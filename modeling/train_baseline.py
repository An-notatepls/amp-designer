"""Fit a training-only E. coli activity baseline; never open held-out CSVs.

Usage: python train_baseline.py PATH/TO/FIXED_SPLIT_DIRECTORY

Inputs: train.csv and the split manifest/group/assignment metadata produced by
prepare_splits.py. Features describe sequence plus C-terminal amidation, not
strain, notes, IDs or assay outcomes. Scaling is fitted inside each CV fold.
Select C by mean five-fold grouped balanced accuracy; ties favor smaller C.
The selected CV score is a development score, not unbiased final performance.
No prediction on test or clinical data is made by this script.
"""

import argparse
import csv
import hashlib
import io
import json
import os
import platform
import shutil
import tempfile
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import scipy
import sklearn
import threadpoolctl
from sklearn.base import clone
from sklearn.dummy import DummyClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    confusion_matrix, matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

VERSION = "physchem-logistic-v1"
SEED, N_FOLDS = 42, 5
C_VALUES = (0.01, 0.1, 1.0, 10.0)
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
FEATURES = {
    "length": "Number of residues",
    "basic_fraction": "Fraction K or R",
    "acidic_fraction": "Fraction D or E",
    "histidine_fraction": "Fraction H",
    "hydrophobic_fraction": "Fraction A, V, I, L, M, F, W or Y (declared grouping)",
    "aromatic_fraction": "Fraction F, W or Y",
    "glycine_fraction": "Fraction G",
    "proline_fraction": "Fraction P",
    "c_amidated": "1 for recorded C-terminal amidation; otherwise 0",
    "charge_proxy": "count(K,R) - count(D,E) + amidation; not a pH-specific charge estimate",
}


def sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)):
            raise ValueError(f"Invalid CSV header: {path}")
        rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f"Malformed CSV: {path}")
    return rows


def identity(row):
    return row["sequence"], row["c_amidated"]


def feature_vector(sequence, c_amidated):
    if not isinstance(sequence, str) or not 10 <= len(sequence) <= 50 or not set(sequence) <= STANDARD:
        raise ValueError("Features require a canonical uppercase sequence of length 10-50")
    if c_amidated not in ("True", "False"):
        raise ValueError("c_amidated must be the CSV string True or False")
    length = len(sequence)
    count = Counter(sequence)
    fraction = lambda residues: sum(count[a] for a in residues) / length
    amidated = int(c_amidated == "True")
    return [length, fraction("KR"), fraction("DE"), fraction("H"),
            fraction("AVILMFWY"), fraction("FWY"), fraction("G"),
            fraction("P"), amidated, count["K"] + count["R"] - count["D"] - count["E"] + amidated]


def feature_matrix(rows):
    # Deliberately ignore labels, group IDs, notes, strain and peptide IDs.
    return np.asarray([feature_vector(r["sequence"], r["c_amidated"]) for r in rows], dtype=float)


def load_training(split_dir):
    """Validate the training partition without opening test/clinical CSVs."""
    split_dir = Path(split_dir)
    manifest = json.loads((split_dir / "split_manifest.json").read_text(encoding="utf-8"))
    required_checks = {
        "identity_partition_unique", "no_family_overlap_across_train_test_and_clinical_reservation",
        "all_qualifying_pairs_protected", "clinical_families_absent_from_training",
        "primary_train_and_test_have_both_classes",
    }
    if any(manifest.get("checks", {}).get(k) is not True for k in required_checks):
        raise ValueError("Split manifest does not confirm all required integrity checks")
    for name in ("train.csv", "sequence_groups.csv", "split_assignments.csv"):
        if sha_file(split_dir / name) != manifest.get("output_sha256", {}).get(name):
            raise ValueError(f"Split file changed or is not recorded in the manifest: {name}")
    rows = read_csv(split_dir / "train.csv")
    required = {"sequence", "c_amidated", "label", "status", "model_eligible", "group_id", "partition", "peptide_ids_json"}
    if not rows or any(not required <= set(r) for r in rows):
        raise ValueError("Training table is empty or missing required columns")
    seen = set()
    for row in rows:
        if row["partition"] != "train" or row["model_eligible"] != "True":
            raise ValueError("Training table contains a non-training or ineligible row")
        if row["label"] not in {"0", "1"} or row["status"] != ("active" if row["label"] == "1" else "inactive"):
            raise ValueError("Training label/status is invalid")
        if not row["group_id"] or identity(row) in seen:
            raise ValueError("Missing family ID or duplicate sequence+amidation identity")
        seen.add(identity(row))
        ids = json.loads(row["peptide_ids_json"])
        if not isinstance(ids, list) or not ids or any(type(v) is not int or v < 1 for v in ids) or len(ids) != len(set(ids)):
            raise ValueError("Invalid peptide ID provenance")
        feature_vector(row["sequence"], row["c_amidated"])
    mapping_rows = read_csv(split_dir / "sequence_groups.csv")
    mapping = {r["sequence"]: r for r in mapping_rows}
    if len(mapping) != len(mapping_rows):
        raise ValueError("Duplicate sequence in group index")
    assignments = read_csv(split_dir / "split_assignments.csv")
    index = {identity(r): r for r in assignments}
    if len(index) != len(assignments) or {identity(r) for r in assignments if r["partition"] == "train"} != seen:
        raise ValueError("Training identity set disagrees with assignments")
    family_roles = {}
    for r in assignments:
        role = "clinical_reservation" if r["partition"].startswith("clinical_") else r["partition"]
        family_roles.setdefault(r["group_id"], set()).add(role)
    if any(len(roles) != 1 for roles in family_roles.values()):
        raise ValueError("A family crosses training/test/clinical boundaries")
    for row in rows:
        entry = mapping.get(row["sequence"], {})
        assigned = index[identity(row)]
        if (entry.get("group_id") != row["group_id"] or entry.get("group_partition") != "train"
                or assigned["group_id"] != row["group_id"] or assigned["partition"] != "train"
                or assigned["peptide_ids_json"] != row["peptide_ids_json"]):
            raise ValueError("Training row disagrees with group/identity metadata")
    actual_counts = {"examples": len(rows), "unique_sequences": len({r["sequence"] for r in rows}),
                     "groups": len({r["group_id"] for r in rows}),
                     "active": sum(r["label"] == "1" for r in rows),
                     "inactive": sum(r["label"] == "0" for r in rows)}
    if actual_counts != manifest.get("counts", {}).get("train"):
        raise ValueError("Training counts disagree with split manifest")
    if {r["label"] for r in rows} != {"0", "1"}:
        raise ValueError("Training requires both classes")
    return rows, manifest


def make_folds(y, groups):
    splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    folds = list(splitter.split(np.zeros((len(y), 1)), y, groups))
    coverage = np.zeros(len(y), dtype=int)
    for train, valid in folds:
        if set(groups[train]) & set(groups[valid]):
            raise ValueError("Family overlap in a cross-validation fold")
        if set(y[train]) != {0, 1} or set(y[valid]) != {0, 1}:
            raise ValueError("A fixed CV fold lacks a class; review the CV design instead of searching seeds")
        coverage[valid] += 1
    if not np.all(coverage == 1):
        raise ValueError("Every training example must be validated exactly once")
    return folds


def logistic_pipeline(c):
    return Pipeline([("scale", StandardScaler()),
                     ("classifier", LogisticRegression(C=c, solver="lbfgs", max_iter=5000,
                                                       class_weight=None, random_state=SEED))])


def metrics(y, predicted, probability):
    tn, fp, fn, tp = (int(v) for v in confusion_matrix(y, predicted, labels=[0, 1]).ravel())
    return {"balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
            "accuracy": float(accuracy_score(y, predicted)),
            "mcc": float(matthews_corrcoef(y, predicted)),
            "roc_auc": float(roc_auc_score(y, probability)),
            "average_precision": float(average_precision_score(y, probability)),
            "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp),
            "tn": tn, "fp": fp, "fn": fn, "tp": tp}


def evaluate_cv(estimator, x, y, folds, model_name, c=None):
    scores, predictions = [], []
    for fold, (train, valid) in enumerate(folds, 1):
        fitted = clone(estimator)
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            fitted.fit(x[train], y[train])
        predicted = fitted.predict(x[valid])
        probability = fitted.predict_proba(x[valid])[:, list(fitted.classes_).index(1)]
        scores.append({"model": model_name, "C": c, "fold": fold, **metrics(y[valid], predicted, probability)})
        for index, label, prob in zip(valid, predicted, probability):
            predictions.append({"model": model_name, "C": c, "fold": fold, "row_index": int(index),
                                "true_label": int(y[index]), "predicted_label": int(label), "probability_active": float(prob)})
    return scores, predictions


def summary(scores):
    keys = ("balanced_accuracy", "accuracy", "mcc", "roc_auc", "average_precision", "sensitivity", "specificity")
    return {key: {"mean": float(np.mean([r[key] for r in scores])),
                  "std_across_folds": float(np.std([r[key] for r in scores], ddof=1))} for key in keys}


def csv_text(rows):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def versions():
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "scikit-learn": sklearn.__version__, "joblib": joblib.__version__, "threadpoolctl": threadpoolctl.__version__}


def train_baseline(split_dir, results_root=None, progress=False):
    split_dir = Path(split_dir).expanduser().resolve()
    rows, split_manifest = load_training(split_dir)
    spec = {"version": VERSION, "split_manifest_sha256": sha_file(split_dir / "split_manifest.json"),
            "train_sha256": sha_file(split_dir / "train.csv"), "script_sha256": sha_file(Path(__file__)),
            "environment": versions(), "features": FEATURES, "C_values": list(C_VALUES),
            "cv": {"method": "StratifiedGroupKFold", "n_splits": N_FOLDS, "shuffle": True, "seed": SEED},
            "selection": "highest mean fold balanced accuracy; exact ties choose smaller C",
            "probability_threshold": 0.5, "class_weight": None,
            "example_weighting": "each sequence+amidation identity equal; assay counts not used as weights"}
    digest = hashlib.sha256(canonical(spec).encode()).hexdigest()[:16]
    root = Path(results_root) if results_root is not None else Path(__file__).resolve().parent / "results" / "baselines"
    output = root / f"{VERSION}-{digest}"
    if output.exists():
        report = json.loads((output / "training_report.json").read_text(encoding="utf-8"))
        if report.get("spec") != spec:
            raise ValueError("Existing result has a different specification")
        for name, expected in report["output_sha256"].items():
            if sha_file(output / name) != expected:
                raise ValueError(f"Existing output changed; preserved without overwrite: {name}")
        return output, report
    x = feature_matrix(rows)
    y = np.asarray([int(r["label"]) for r in rows])
    groups = np.asarray([r["group_id"] for r in rows])
    folds = make_folds(y, groups)
    all_scores, all_predictions = evaluate_cv(DummyClassifier(strategy="most_frequent"), x, y, folds, "majority_dummy")
    dummy_summary = summary(all_scores)
    candidates = []
    for c in C_VALUES:
        scores, predictions = evaluate_cv(logistic_pipeline(c), x, y, folds, "logistic", c)
        all_scores.extend(scores)
        all_predictions.extend(predictions)
        stats = summary(scores)
        candidates.append({"C": c, "cv": stats})
        if progress:
            print(f"C={c:g}: grouped CV balanced accuracy {stats['balanced_accuracy']['mean']:.3f}", flush=True)
    selected = max(candidates, key=lambda r: (r["cv"]["balanced_accuracy"]["mean"], -r["C"]))
    fitted = logistic_pipeline(selected["C"])
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        fitted.fit(x, y)
    assignment_rows = []
    for fold, (_, valid) in enumerate(folds, 1):
        for index in valid:
            assignment_rows.append({"row_index": int(index), "sequence": rows[index]["sequence"],
                                    "c_amidated": rows[index]["c_amidated"], "group_id": rows[index]["group_id"],
                                    "label": int(y[index]), "validation_fold": fold})
    assignment_rows.sort(key=lambda r: r["row_index"])
    for item in all_predictions:
        row = rows[item["row_index"]]
        item.update({"sequence": row["sequence"], "c_amidated": row["c_amidated"], "group_id": row["group_id"]})
    features = [{"row_index": i, "sequence": r["sequence"], "c_amidated_input": r["c_amidated"],
                 "label": int(y[i]), "group_id": r["group_id"], **dict(zip(FEATURES, x[i].tolist()))} for i, r in enumerate(rows)]
    scaler = fitted.named_steps["scale"]
    classifier = fitted.named_steps["classifier"]
    coefficients = [{"feature": name, "coefficient_per_training_sd": float(coef), "training_mean": float(mean),
                     "training_sd_scale": float(scale), "training_variance": float(var)}
                    for name, coef, mean, scale, var in zip(FEATURES, classifier.coef_[0], scaler.mean_, scaler.scale_, scaler.var_)]
    report = {"spec": spec, "trained_at_utc": datetime.now(timezone.utc).isoformat(),
              "training_counts": split_manifest["counts"]["train"],
              "dummy_cv": dummy_summary, "logistic_candidates": candidates,
              "selected_C": selected["C"], "selected_logistic_cv": selected["cv"],
              "intercept": float(classifier.intercept_[0]),
              "heldout_csvs_opened": [], "final_test_evaluated": False,
              "checks": {"split_training_hash_verified": True, "group_disjoint_cv": True,
                         "both_classes_in_each_fold": True, "each_row_validated_once": True,
                         "scaling_inside_cv_pipeline": True},
              "limitations": ["Selected CV scores reuse folds for model selection and are optimistic development estimates.",
                              "Fold standard deviations are descriptive, not confidence intervals.",
                              "Features are correlated; coefficients describe this fitted model, not causal residue effects.",
                              "Probability outputs are uncalibrated model scores, not proven biological success probabilities.",
                              "Similarity groups do not guarantee all biological relationships or shared study effects are isolated.",
                              "The curated threshold is an operational label, not a clinical breakpoint or evidence of bacterial killing."],
              "source_attribution": "Derived from DBAASP data, https://dbaasp.org; see project README for citation and terms."}
    payloads = {"cv_scores.csv": csv_text(all_scores), "cv_predictions.csv": csv_text(all_predictions),
                "cv_folds.csv": csv_text(assignment_rows), "train_features.csv": csv_text(features),
                "coefficients.csv": csv_text(coefficients),
                "feature_definitions.json": json.dumps(FEATURES, indent=2) + "\n"}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".baseline-preparing-", dir=output.parent))
    try:
        for name, content in payloads.items():
            (temporary / name).write_text(content, encoding="utf-8", newline="")
        joblib.dump({"pipeline": fitted, "feature_names": list(FEATURES), "spec": spec}, temporary / "model.joblib")
        report["output_sha256"] = {p.name: sha_file(p) for p in sorted(temporary.iterdir())}
        (temporary / "training_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        for path in temporary.iterdir():
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
        temporary.rename(output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_dir", type=Path)
    parser.add_argument("--results-root", type=Path)
    args = parser.parse_args()
    output, report = train_baseline(args.split_dir, args.results_root, progress=True)
    print("Training examples:", report["training_counts"]["examples"])
    print("Selected C:", report["selected_C"])
    for name, key in (("Majority dummy", "dummy_cv"), ("Selected logistic", "selected_logistic_cv")):
        cv = report[key]
        print(f"{name}: balanced accuracy {cv['balanced_accuracy']['mean']:.3f} "
              f"(fold SD {cv['balanced_accuracy']['std_across_folds']:.3f}); "
              f"ROC AUC {cv['roc_auc']['mean']:.3f}; average precision {cv['average_precision']['mean']:.3f}")
    print("Results:", output)
    print("Development CV scores only. Test and clinical CSVs were not opened.")
    print("The selected model was refitted on all training rows; do not claim final test performance yet.")


if __name__ == "__main__":
    main()
