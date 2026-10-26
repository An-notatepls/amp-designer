"""Compare feature sets on saved training folds; keep held-out CSVs sealed.

Usage: python modeling/compare_features.py SPLIT_DIR --baseline-dir BASELINE_DIR

Reuses cv_folds.csv exactly. C is selected separately for every feature set
using mean fold balanced accuracy. Paired cluster bootstrap intervals describe
fixed, pooled out-of-fold predictions after selection: they do NOT account for
model-selection uncertainty, refitting, or repeated reuse of the same folds.
Grouping sensitivity is training-only and never reallocates the frozen test.
"""

import argparse
import importlib.util
import json
import shutil
import tempfile
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import joblib
import numpy as np
from sklearn.dummy import DummyClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score, roc_auc_score

import train_baseline as base

VERSION = "feature-comparison-v1"
C_VALUES = (0.0001, 0.001, 0.01, 0.1, 1.0, 10.0)
MODEL_ORDER = ("dummy", "charge_uncorrected", "charge_corrected",
               "charge_corrected_amidation", "composition", "physchem")
AA = "ACDEFGHIKLMNPQRSTVWY"
FEATURE_NAMES = {
    "dummy": ["constant"],
    "charge_uncorrected": ["charge_residue_counts"],
    "charge_corrected": ["charge_proxy"],
    "charge_corrected_amidation": ["charge_proxy", "c_amidated"],
    "composition": [f"fraction_{a}" for a in AA] + ["length", "c_amidated"],
    "physchem": list(base.FEATURES),
}
INTERVAL_METRICS = ("balanced_accuracy", "roc_auc", "average_precision")
THRESHOLDS = (Fraction(3, 5), Fraction(7, 10), Fraction(4, 5))


def matrices(rows):
    physchem = base.feature_matrix(rows)
    amidation = physchem[:, list(base.FEATURES).index("c_amidated")]
    corrected = physchem[:, list(base.FEATURES).index("charge_proxy")]
    composition = [[r["sequence"].count(a) / len(r["sequence"]) for a in AA]
                   + [len(r["sequence"]), int(r["c_amidated"] == "True")] for r in rows]
    return {"dummy": np.ones((len(rows), 1)),
            "charge_uncorrected": (corrected - amidation)[:, None],
            "charge_corrected": corrected[:, None],
            "charge_corrected_amidation": np.column_stack((corrected, amidation)),
            "composition": np.asarray(composition, dtype=float), "physchem": physchem}


def load_saved_folds(baseline_dir, split_dir, rows):
    baseline_dir, split_dir = Path(baseline_dir), Path(split_dir)
    report_path = baseline_dir / "training_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report["spec"]["train_sha256"] != base.sha_file(split_dir / "train.csv"):
        raise ValueError("Baseline used a different training table")
    if report["spec"]["split_manifest_sha256"] != base.sha_file(split_dir / "split_manifest.json"):
        raise ValueError("Baseline used a different split manifest")
    path = baseline_dir / "cv_folds.csv"
    if base.sha_file(path) != report["output_sha256"].get("cv_folds.csv"):
        raise ValueError("Saved baseline folds changed")
    assignments = base.read_csv(path)
    if len(assignments) != len(rows):
        raise ValueError("Saved fold count differs from training table")
    fold_id = np.zeros(len(rows), dtype=int)
    seen = set()
    for entry in assignments:
        index, fold = int(entry["row_index"]), int(entry["validation_fold"])
        if index in seen or not 0 <= index < len(rows) or not 1 <= fold <= base.N_FOLDS:
            raise ValueError("Invalid or duplicate saved fold assignment")
        seen.add(index)
        for key in ("sequence", "c_amidated", "group_id", "label"):
            if entry[key] != rows[index][key]:
                raise ValueError(f"Saved fold row disagrees with training data: {key}")
        fold_id[index] = fold
    folds = [(np.flatnonzero(fold_id != f), np.flatnonzero(fold_id == f)) for f in range(1, base.N_FOLDS + 1)]
    y, groups = labels_groups(rows)
    for train, valid in folds:
        if set(groups[train]) & set(groups[valid]):
            raise ValueError("Family crosses saved fold boundaries")
        if set(y[train]) != {0, 1} or set(y[valid]) != {0, 1}:
            raise ValueError("Saved fold lacks a class")
    return folds, report


def labels_groups(rows):
    return np.asarray([int(r["label"]) for r in rows]), np.asarray([r["group_id"] for r in rows])


def compare(x_by_model, y, folds, progress=False):
    selected, all_scores, selected_predictions, fitted_models = {}, [], {}, {}
    for model, x in x_by_model.items():
        candidates = []
        for c in ((None,) if model == "dummy" else C_VALUES):
            estimator = DummyClassifier(strategy="most_frequent") if c is None else base.logistic_pipeline(c)
            scores, predictions = base.evaluate_cv(estimator, x, y, folds, model, c)
            all_scores.extend(scores)
            candidates.append({"C": c, "summary": base.summary(scores), "predictions": predictions})
        best = max(candidates, key=lambda r: (r["summary"]["balanced_accuracy"]["mean"], -(r["C"] or 0)))
        predictions = sorted(best["predictions"], key=lambda r: r["row_index"])
        if [r["row_index"] for r in predictions] != list(range(len(y))):
            raise ValueError("Out-of-fold predictions do not cover training exactly once")
        probability = np.asarray([r["probability_active"] for r in predictions])
        predicted = np.asarray([r["predicted_label"] for r in predictions])
        selected[model] = {"C": best["C"], "n_features": x.shape[1], "cv": best["summary"],
                           "pooled_oof": base.metrics(y, predicted, probability)}
        selected_predictions[model] = predictions
        estimator = DummyClassifier(strategy="most_frequent") if model == "dummy" else base.logistic_pipeline(best["C"])
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            estimator.fit(x, y)
        fitted_models[model] = estimator
        if progress:
            print(f"{model}: C={best['C']}; mean fold balanced accuracy "
                  f"{selected[model]['cv']['balanced_accuracy']['mean']:.3f}; "
                  f"mean fold AUROC {selected[model]['cv']['roc_auc']['mean']:.3f}", flush=True)
    return selected, all_scores, selected_predictions, fitted_models


def bootstrap_weights(groups, draws, seed=42):
    """Each draw samples G whole groups with replacement, keeping all their rows."""
    unique, inverse = np.unique(groups, return_inverse=True)
    if len(unique) < 2:
        raise ValueError("Bootstrap requires at least two groups")
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(unique), np.full(len(unique), 1 / len(unique)), size=draws)
    return counts[:, inverse]


def weighted_metrics(y, predicted, probability, weights):
    positive, negative = y == 1, y == 0
    sensitivity = np.sum(weights * positive * (predicted == 1)) / np.sum(weights * positive)
    specificity = np.sum(weights * negative * (predicted == 0)) / np.sum(weights * negative)
    return np.array([(sensitivity + specificity) / 2,
                     roc_auc_score(y, probability, sample_weight=weights),
                     average_precision_score(y, probability, sample_weight=weights)])


def paired_bootstrap(y, groups, predictions, draws=2000, progress=False):
    if draws < 20:
        raise ValueError("Use at least 20 bootstrap draws (2000 recommended)")
    weights = bootstrap_weights(groups, draws)
    valid = (weights[:, y == 1].sum(axis=1) > 0) & (weights[:, y == 0].sum(axis=1) > 0)
    weights = weights[valid]
    if len(weights) < 0.9 * draws:
        raise ValueError("Too many one-class cluster resamples; intervals are not supported by these data")
    samples = {}
    for model, rows in predictions.items():
        probability = np.asarray([r["probability_active"] for r in rows])
        predicted = np.asarray([r["predicted_label"] for r in rows])
        samples[model] = np.asarray([weighted_metrics(y, predicted, probability, w) for w in weights])
        if progress:
            print(f"Bootstrap complete: {model} ({len(weights)} shared valid draws)", flush=True)
    intervals = []
    deltas = []
    names = list(predictions)
    for index, model in enumerate(names):
        for j, metric in enumerate(INTERVAL_METRICS):
            low, high = np.quantile(samples[model][:, j], [0.025, 0.975])
            intervals.append({"model": model, "metric": metric, "low": float(low), "high": float(high)})
            if index:
                previous = names[index - 1]
                delta = samples[model][:, j] - samples[previous][:, j]
                low, high = np.quantile(delta, [0.025, 0.975])
                current_p = predictions[model]
                previous_p = predictions[previous]
                score = lambda p: weighted_metrics(y, np.array([r["predicted_label"] for r in p]),
                                                  np.array([r["probability_active"] for r in p]), np.ones(len(y)))[j]
                deltas.append({"model": model, "previous_model": previous, "metric": metric,
                               "pooled_oof_delta": float(score(current_p) - score(previous_p)),
                               "low": float(low), "high": float(high)})
    return intervals, deltas, {"requested_draws": draws, "valid_shared_draws": len(weights),
                              "one_class_draws_omitted": int(draws - len(weights)), "seed": 42}


def grouping_module():
    path = Path(__file__).resolve().with_name("prepare_splits.py")
    if not path.is_file():
        path = Path(__file__).resolve().parent.parent / "prepare_splits.py"
    if not path.is_file():
        raise ValueError("Keep prepare_splits.py in the project root, beside the modeling folder")
    spec = importlib.util.spec_from_file_location("amp_grouping", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def cluster_description(sizes):
    values = list(sizes.values())
    return {"n_sequences": sum(values), "n_groups": len(values),
            "largest_group": max(values), "largest_fraction": max(values) / sum(values),
            "singleton_groups": sum(v == 1 for v in values),
            "median_group_size": float(np.median(values)),
            "p90_group_size": float(np.quantile(values, 0.9)),
            "size_histogram": {str(size): count for size, count in sorted(Counter(values).items())}}


def grouping_sensitivity(split_dir, rows, saved_folds, progress=False):
    """Regroup unlabeled audit sequences; evaluate only the fixed training pool."""
    module, _ = grouping_module()
    index = base.read_csv(Path(split_dir) / "sequence_groups.csv")
    by_sequence = {r["sequence"]: r for r in index}
    assignments = base.read_csv(Path(split_dir) / "split_assignments.csv")
    y, _ = labels_groups(rows)
    matrices_by_model = matrices(rows)
    records, mappings, fold_rows = [], [], []
    for threshold in THRESHOLDS:
        if progress:
            print(f"Checking grouping at {float(threshold):.0%} global edit similarity", flush=True)
        mapping, sizes, _ = module.sequence_groups(list(by_sequence), threshold)
        roles = defaultdict(set)
        for row in assignments:
            role = "clinical_reservation" if row["partition"].startswith("clinical_") else row["partition"]
            roles[mapping[row["sequence"]]].add(role)
        crossing = {g for g, parts in roles.items() if "train" in parts and len(parts) > 1}
        description = cluster_description(sizes)
        train_sizes = Counter(mapping[r["sequence"]] for r in rows)
        record = {"similarity": str(threshold), "audit_clusters": description,
                  "training_groups": len(train_sizes), "largest_training_group_examples": max(train_sizes.values()),
                  "groups_linking_training_to_frozen_holdouts": len(crossing),
                  "training_examples_in_those_groups": sum(mapping[r["sequence"]] in crossing for r in rows)}
        mappings.extend({"similarity": str(threshold), "sequence": sequence, "group_id": group,
                         "audit_group_size": sizes[group]} for sequence, group in sorted(mapping.items()))
        groups = np.asarray([mapping[r["sequence"]] for r in rows])
        if threshold == Fraction(4, 5):
            if any(mapping[s] != by_sequence[s]["group_id"] for s in mapping):
                raise ValueError("80% grouping no longer matches the frozen split")
            folds = saved_folds
            record["fold_policy"] = "original saved shared folds"
        else:
            record["fold_policy"] = "fixed-seed StratifiedGroupKFold within original training pool"
            if len(train_sizes) < base.N_FOLDS:
                record["status"] = "not_evaluable: fewer groups than folds"
                records.append(record)
                continue
            try:
                folds = base.make_folds(y, groups)
            except ValueError as error:
                record["status"] = f"not_evaluable: {error}"
                records.append(record)
                continue
        for fold, (_, valid) in enumerate(folds, 1):
            for i in valid:
                fold_rows.append({"similarity": str(threshold), "row_index": int(i),
                                  "sequence": rows[i]["sequence"], "c_amidated": rows[i]["c_amidated"],
                                  "group_id": groups[i], "label": int(y[i]), "validation_fold": fold})
        result, sensitivity_scores, _, _ = compare({"dummy": matrices_by_model["dummy"], "physchem": matrices_by_model["physchem"]}, y, folds)
        record.update({"status": "evaluated_training_only", "models": result, "cv_scores": sensitivity_scores})
        records.append(record)
    return records, mappings, fold_rows


def run_comparison(split_dir, baseline_dir, results_root=None, draws=2000, progress=False):
    split_dir, baseline_dir = Path(split_dir).resolve(), Path(baseline_dir).resolve()
    rows, _ = base.load_training(split_dir)
    folds, _ = load_saved_folds(baseline_dir, split_dir, rows)
    _, group_path = grouping_module()
    spec = {"version": VERSION, "script_sha256": base.sha_file(Path(__file__)),
            "baseline_script_sha256": base.sha_file(Path(base.__file__)),
            "grouping_script_sha256": base.sha_file(group_path),
            "train_sha256": base.sha_file(split_dir / "train.csv"),
            "split_manifest_sha256": base.sha_file(split_dir / "split_manifest.json"),
            "saved_cv_folds_sha256": base.sha_file(baseline_dir / "cv_folds.csv"),
            "baseline_report_sha256": base.sha_file(baseline_dir / "training_report.json"),
            "environment": base.versions(), "feature_names": FEATURE_NAMES,
            "C_values": list(C_VALUES), "bootstrap_draws": draws,
            "bootstrap_seed": 42, "similarity_thresholds": [str(t) for t in THRESHOLDS]}
    digest = __import__("hashlib").sha256(base.canonical(spec).encode()).hexdigest()[:16]
    if results_root is None:
        results_root = baseline_dir.parent / "comparisons"
    output = Path(results_root) / f"{VERSION}-{digest}"
    if output.exists():
        report = json.loads((output / "comparison_report.json").read_text(encoding="utf-8"))
        if report["spec"] != spec:
            raise ValueError("Existing comparison differs; refusing overwrite")
        for name, digest in report["output_sha256"].items():
            if base.sha_file(output / name) != digest:
                raise ValueError(f"Existing comparison output changed: {name}")
        return output, report
    y, groups = labels_groups(rows)
    selected, scores, predictions, models = compare(matrices(rows), y, folds, progress)
    intervals, deltas, boot = paired_bootstrap(y, groups, predictions, draws, progress)
    sensitivity, mapping_rows, sensitivity_folds = grouping_sensitivity(split_dir, rows, folds, progress)
    ci = {(r["model"], r["metric"]): r for r in intervals}
    ladder = []
    for model, result in selected.items():
        row = {"model": model, "n_features": result["n_features"], "selected_C": result["C"],
               "C_at_grid_edge": result["C"] in (C_VALUES[0], C_VALUES[-1])}
        for metric in INTERVAL_METRICS:
            row.update({f"mean_fold_{metric}": result["cv"][metric]["mean"],
                        f"pooled_oof_{metric}": result["pooled_oof"][metric],
                        f"pooled_oof_{metric}_low": ci[model, metric]["low"],
                        f"pooled_oof_{metric}_high": ci[model, metric]["high"]})
        ladder.append(row)
    prediction_rows = []
    for model, items in predictions.items():
        for item in items:
            r = rows[item["row_index"]]
            prediction_rows.append({**item, "sequence": r["sequence"], "c_amidated": r["c_amidated"], "group_id": r["group_id"]})
    report = {"spec": spec, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
              "training_examples": len(rows), "training_groups": len(set(groups)),
              "models": selected, "bootstrap": boot, "grouping_sensitivity": sensitivity,
              "heldout_csvs_opened": [], "final_test_evaluated": False,
              "methods": {
                  "selection": "maximum mean fold balanced accuracy; exact ties favor smaller C",
                  "scaling": "StandardScaler fitted inside each training fold",
                  "intervals": "95% percentile paired whole-family bootstrap of fixed pooled OOF predictions; identical resamples for all models",
                  "sensitivity": "same frozen training examples; regroup all unlabeled audit sequences, then form training CV folds; no holdout reassignment",
                  "primary_similarity": "1 - global unit-cost Levenshtein distance / longer length; 80%; connected components"},
              "limitations": [
                  "Model and C comparisons use development folds; their scores and bootstrap intervals are not final generalization estimates.",
                  "Bootstrap holds selected C and fitted fold predictions fixed; it omits training/selection uncertainty and shared-fold dependence.",
                  "Mean fold scores select C; pooled OOF scores have different weighting and are the targets of bootstrap intervals.",
                  "Charge is a count-based proxy; corrected means +1 for amidation, not a full pH-specific calculation.",
                  "Corrected charge already contains amidation information; dropping its explicit flag does not remove all amidation information.",
                  "Composition includes length and amidation; composition and physchem are alternative feature sets, not nested steps.",
                  "Correlated features and L2 regularization affect comparisons; they do not establish causal chemical effects.",
                  "Lower similarity cutoffs can connect original training and held-out families; diagnostic counts expose this without changing the original split.",
                  "Threshold CV results use different folds and are sensitivity analyses, not a rule for choosing the best-looking split.",
                  "Do not claim charge explains most of the signal unless the measured comparisons support it."],
              "source_attribution": "Derived from DBAASP data, https://dbaasp.org; see project README for citation and terms."}
    payloads = {"feature_ladder.csv": base.csv_text(ladder), "cv_scores.csv": base.csv_text(scores),
                "selected_oof_predictions.csv": base.csv_text(prediction_rows),
                "paired_differences.csv": base.csv_text(deltas),
                "cluster_mapping.csv": base.csv_text(mapping_rows),
                "grouping_cv_folds.csv": base.csv_text(sensitivity_folds),
                "grouping_sensitivity.json": json.dumps(sensitivity, indent=2, allow_nan=False) + "\n"}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".comparison-preparing-", dir=output.parent))
    try:
        for name, content in payloads.items():
            (temporary / name).write_text(content, encoding="utf-8", newline="")
        # Preserve exact source folds without regenerating them.
        shutil.copyfile(baseline_dir / "cv_folds.csv", temporary / "cv_folds.csv")
        joblib.dump({"models": models, "feature_names": FEATURE_NAMES, "spec": spec}, temporary / "models.joblib")
        report["output_sha256"] = {p.name: base.sha_file(p) for p in sorted(temporary.iterdir())}
        (temporary / "comparison_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.rename(output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_dir", type=Path)
    parser.add_argument("--baseline-dir", required=True, type=Path)
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    args = parser.parse_args()
    if args.bootstrap_draws < 20:
        parser.error("At least 20 draws are required; keep the default 2000 for the report")
    output, report = run_comparison(args.split_dir, args.baseline_dir, args.results_root, args.bootstrap_draws, progress=True)
    print("Results:", output)
    print("Reused original folds. Test and clinical CSVs were not opened.")
    print("All scores and intervals describe development predictions, not final test performance.")


if __name__ == "__main__":
    main()
