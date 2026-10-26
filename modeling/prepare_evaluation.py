"""Freeze a pre-test protocol; audit sequences and training OOF predictions only.

Usage: python modeling/prepare_evaluation.py SPLIT_DIR --comparison-dir RUN_DIR

Never opens test.csv or either clinical CSV. Test sequences/identities come
from the label-free split_assignments.csv. No model is fitted or evaluated on
held-out outcomes. Commit the generated protocol before final evaluation.
"""

import argparse
import hashlib
import itertools
import json
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import Bio
import joblib
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

import compare_features as compare
import train_baseline as base

VERSION = "pretest-v1"
PRIMARY = "charge_corrected_amidation"
EXPECTED_C = {"dummy": None, "charge_uncorrected": 0.1, "charge_corrected": 1.0,
              PRIMARY: 0.1, "composition": 0.1, "physchem": 0.01}
KS = (5, 10, 20)
TIE_LIMIT = 100
POLICY = {
    "decision_stage": "chosen after development comparisons, before test outcomes; not retrospective preregistration",
    "primary_model": PRIMARY, "primary_C": 0.1,
    "primary_features": ["charge_proxy", "c_amidated"],
    "charge_proxy": "count(K,R) - count(D,E) + int(c_amidated)",
    "prediction_threshold": 0.5, "prediction_tie": "score exactly 0.5 predicts inactive, matching sklearn logistic decision > 0",
    "primary_metric": "balanced_accuracy",
    "reported_models_C": EXPECTED_C,
    "secondary_models_C": {model: c for model, c in EXPECTED_C.items() if model != PRIMARY},
    "secondary_metrics": ["roc_auc", "average_precision", "mcc", "accuracy", "precision", "sensitivity", "specificity", "confusion_counts"],
    "ranking_budgets": list(KS), "ranking_status": "illustrative candidate budgets fixed before inspecting P@k; not a validated lab budget",
    "ranking_modes": ["all_identities", "one_representative_per_80pct_family"],
    "ranking_ties": "round probabilities to 12 decimals; expected hits from uniform selection within a boundary tie",
    "family_representative": "highest rounded probability; ties choose lexicographically smallest sequence+amidation identity, ignoring labels",
    "local_alignment": {"mode": "local", "matrix": "BLOSUM62", "open_gap_score": -10.0,
                        "extend_gap_score": -0.5, "identity_denominator": "all alignment columns including internal gaps",
                        "coverage_denominator": "aligned residue pairs divided by peptide length; report both peptides",
                        "minimum_aligned_residue_pairs": 10, "minimum_shorter_peptide_coverage": 0.8,
                        "maximum_optimal_alignments_per_pair": TIE_LIMIT,
                        "tie_policy": "scan up to 100 score-optimal alignments; an incomplete tie scan excludes the test identity from local/combined strict subsets"},
    "subsets": {
        "overall": "all original primary test identities",
        "global60_strict": "maximum direct global edit similarity to training is strictly below 0.60",
        "local60_strict": "no detected coverage-qualified score-optimal local match has identity >=0.60; all optimal-alignment tie scans complete",
        "combined60_strict": "intersection of global60_strict and local60_strict",
    },
    "primary_analysis": "primary model, original overall test; subsets and other models are secondary, not selection rules",
    "final_intervals": {"method": "paired percentile whole-80pct-family bootstrap of fixed test predictions",
                        "draws": 5000, "seed": 42, "level": 0.95,
                        "minimum_groups_for_auc_interval": 10, "minimum_examples_per_class_for_auc_interval": 5,
                        "one_class_draws": "omit consistently across compared models; disclose counts; suppress interval if >10% omitted"},
    "small_subsets": "always report n, classes and families; AUROC undefined with one class; suppress AUROC intervals below declared support limits",
    "clinical_analysis": "17 clinical-only identities are exploratory; only three inactive; do not use clinical buffer for training or selection",
    "selection_rule": "evaluate the fixed list once; select no feature set, C, cutoff or grouping threshold from test/clinical performance",
    "power": "71 test identities in 47 families provide a limited final check, not a powered ranking of small feature-set differences; no promised CI width",
}


def alignment_summary(alignment, first, second):
    coordinates = alignment.coordinates
    columns = pairs = matches = 0
    for i in range(coordinates.shape[1] - 1):
        a0, a1 = (int(v) for v in coordinates[0, i:i + 2])
        b0, b1 = (int(v) for v in coordinates[1, i:i + 2])
        da, db = a1 - a0, b1 - b0
        columns += max(da, db)
        if da and db:
            if da != db:
                raise ValueError("Unexpected non-diagonal paired alignment segment")
            pairs += da
            matches += sum(a == b for a, b in zip(first[a0:a1], second[b0:b1]))
    return {"score": float(alignment.score), "matches": matches, "columns": columns,
            "aligned_pairs": pairs, "identity": matches / columns if columns else 0.0,
            "first_coverage": pairs / len(first), "second_coverage": pairs / len(second),
            "shorter_coverage": pairs / min(len(first), len(second))}


def qualifying(local, shorter_length):
    return local["aligned_pairs"] >= 10 and 5 * local["aligned_pairs"] >= 4 * shorter_length


def local_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def local_pair(first, second, aligner, tie_limit=TIE_LIMIT):
    raw = qualified = None
    complete = True
    rank = lambda r: (r["identity"], r["aligned_pairs"], r["score"])
    for i, alignment in enumerate(itertools.islice(aligner.align(first, second), tie_limit + 1)):
        if i == tie_limit:
            complete = False
            break
        item = alignment_summary(alignment, first, second)
        if raw is None or rank(item) > rank(raw):
            raw = item
        if qualifying(item, min(len(first), len(second))) and (qualified is None or rank(item) > rank(qualified)):
            qualified = item
    return raw, qualified, complete


def novelty_flags(distance, longer_length, qualified, complete):
    global_strict = 5 * distance > 2 * longer_length
    local_strict = complete and (qualified is None or 5 * qualified["matches"] < 3 * qualified["columns"])
    return global_strict, local_strict, global_strict and local_strict


def neighbor_audit(training_rows, assignments, progress=False):
    module, _ = compare.grouping_module()
    training = sorted({r["sequence"] for r in training_rows})
    targets = sorted((r for r in assignments if r["partition"] == "test"), key=base.identity)
    if not training or not targets:
        raise ValueError("Need nonempty training and primary-test assignment metadata")
    aligner = local_aligner()
    output = []
    for index, target in enumerate(targets, 1):
        sequence = target["sequence"]
        best_global = None
        best_raw = best_qualified = None
        complete = True
        truncated_pairs = 0
        for neighbor in training:
            distance = module.edit_distance(sequence, neighbor)
            longer = max(len(sequence), len(neighbor))
            if best_global is None or distance * best_global["longer"] < best_global["distance"] * longer:
                best_global = {"distance": distance, "longer": longer, "sequence": neighbor}
            raw, qualified, scanned = local_pair(sequence, neighbor, aligner)
            complete = complete and scanned
            truncated_pairs += int(not scanned)
            for kind, item in (("raw", raw), ("qualified", qualified)):
                if item is None:
                    continue
                item = {**item, "training_sequence": neighbor}
                previous = best_raw if kind == "raw" else best_qualified
                if previous is None or (item["identity"], item["aligned_pairs"], item["score"]) > (previous["identity"], previous["aligned_pairs"], previous["score"]):
                    if kind == "raw":
                        best_raw = item
                    else:
                        best_qualified = item
        if 5 * best_global["distance"] <= best_global["longer"]:
            raise ValueError("Neighbor audit found a train/test pair violating the fixed 80% grouping rule; stop before evaluation")
        g, l, combined = novelty_flags(best_global["distance"], best_global["longer"], best_qualified, complete)
        row = {"sequence": sequence, "c_amidated": target["c_amidated"], "group_id": target["group_id"],
               "nearest_global_training_sequence": best_global["sequence"],
               "maximum_global_edit_similarity": 1 - best_global["distance"] / best_global["longer"],
               "global60_strict": g, "local60_strict": l, "combined60_strict": combined,
               "optimal_alignment_scan_complete": complete, "tie_truncated_training_pairs": truncated_pairs}
        for prefix, best in (("unrestricted_local", best_raw), ("qualified_local", best_qualified)):
            for key in ("training_sequence", "identity", "score", "aligned_pairs", "matches", "columns", "first_coverage", "second_coverage", "shorter_coverage"):
                row[f"{prefix}_{key}"] = best[key] if best is not None else None
        output.append(row)
        if progress and (index % 10 == 0 or index == len(targets)):
            print(f"Sequence-only neighbor audit: {index}/{len(targets)} test identities", flush=True)
    return output


def identity_histogram(values):
    counts = {"[0,.2)": 0, "[.2,.4)": 0, "[.4,.6)": 0, "[.6,.8)": 0, "[.8,1]": 0, "no_qualifying_alignment": 0}
    bins = list(counts)[:-1]
    for value in values:
        if value is None:
            counts["no_qualifying_alignment"] += 1
        else:
            counts[bins[int(np.searchsorted([0.2, 0.4, 0.6, 0.8], value, side="right"))]] += 1
    return counts


def precision_at_budgets(y, probability, budgets=KS, weights=None):
    """Expected precision for uniform boundary-tie selection, with cluster weights."""
    y, probability = np.asarray(y), np.round(np.asarray(probability, dtype=float), 12)
    weights = np.ones(len(y)) if weights is None else np.asarray(weights, dtype=float)
    if (not len(y) or len(y) != len(probability) or len(y) != len(weights)
            or not np.all(np.isin(y, [0, 1])) or not np.all(np.isfinite(probability))
            or np.any((probability < 0) | (probability > 1))
            or not np.all(np.isfinite(weights)) or np.any(weights < 0)):
        raise ValueError("Invalid ranking inputs")
    _, inverse = np.unique(-probability, return_inverse=True)
    count = np.bincount(inverse, weights=weights)
    positives = np.bincount(inverse, weights=weights * y)
    cumulative = np.cumsum(count)
    cumulative_positives = np.cumsum(positives)
    result = {}
    for k in budgets:
        if k < 1 or cumulative[-1] < k:
            result[k] = None
            continue
        j = int(np.searchsorted(cumulative, k, side="left"))
        before = cumulative[j - 1] if j else 0
        before_positives = cumulative_positives[j - 1] if j else 0
        result[k] = float((before_positives + (k - before) * positives[j] / count[j]) / k)
    return result


def family_representatives(rows, probability):
    chosen = {}
    for i, row in enumerate(rows):
        key = row["group_id"]
        old = chosen.get(key)
        rank = (-round(float(probability[i]), 12), base.identity(row))
        if old is None or rank < (-round(float(probability[old]), 12), base.identity(rows[old])):
            chosen[key] = i
    return np.array([chosen[g] for g in sorted(chosen)], dtype=int)


def load_comparison(split_dir, comparison_dir, rows):
    directory = Path(comparison_dir)
    report = json.loads((directory / "comparison_report.json").read_text(encoding="utf-8"))
    spec = report["spec"]
    if spec["train_sha256"] != base.sha_file(Path(split_dir) / "train.csv") or spec["split_manifest_sha256"] != base.sha_file(Path(split_dir) / "split_manifest.json"):
        raise ValueError("Comparison uses a different training table or split")
    for name in ("selected_oof_predictions.csv", "models.joblib", "cv_folds.csv"):
        if base.sha_file(directory / name) != report["output_sha256"].get(name):
            raise ValueError(f"Comparison input changed: {name}")
    bundle = joblib.load(directory / "models.joblib")
    if bundle["spec"] != spec or set(bundle["models"]) != set(EXPECTED_C):
        raise ValueError("Model bundle disagrees with comparison specification")
    folds = base.read_csv(directory / "cv_folds.csv")
    if len(folds) != len(rows):
        raise ValueError("Saved fold table has the wrong row count")
    fold_index = {}
    for entry in folds:
        i = int(entry["row_index"])
        if i in fold_index or not 0 <= i < len(rows):
            raise ValueError("Saved fold has an invalid row index")
        if any(entry[k] != rows[i][k] for k in ("sequence", "c_amidated", "group_id", "label")):
            raise ValueError("Saved fold disagrees with training identity/label")
        fold_index[i] = int(entry["validation_fold"])
    locations = defaultdict(set)
    for i, row in enumerate(rows):
        locations[row["group_id"]].add(fold_index[i])
    if set(fold_index.values()) != set(range(1, 6)) or any(len(v) != 1 for v in locations.values()):
        raise ValueError("Saved validation folds do not preserve groups")
    predictions = defaultdict(list)
    for entry in base.read_csv(directory / "selected_oof_predictions.csv"):
        model = entry["model"]
        i = int(entry["row_index"])
        if model not in EXPECTED_C or not 0 <= i < len(rows):
            raise ValueError("Unexpected prediction model/index")
        if any(entry[k] != rows[i][k] for k in ("sequence", "c_amidated", "group_id")) or entry["true_label"] != rows[i]["label"] or int(entry["fold"]) != fold_index[i]:
            raise ValueError("OOF prediction identity/label/fold mismatch")
        p = float(entry["probability_active"])
        if not np.isfinite(p) or not 0 <= p <= 1 or entry["predicted_label"] not in {"0", "1"}:
            raise ValueError("Invalid OOF probability/prediction")
        c = None if entry["C"] == "" else float(entry["C"])
        if c != EXPECTED_C[model]:
            raise ValueError("Selected C differs from the fixed protocol")
        predictions[model].append(entry)
    y = np.array([int(r["label"]) for r in rows])
    for model, expected_c in EXPECTED_C.items():
        if report["models"][model]["C"] != expected_c or bundle["feature_names"][model] != compare.FEATURE_NAMES[model]:
            raise ValueError("Model configuration differs from the fixed protocol")
        entries = sorted(predictions[model], key=lambda r: int(r["row_index"]))
        if [int(r["row_index"]) for r in entries] != list(range(len(rows))):
            raise ValueError("OOF predictions do not cover training exactly once")
        predictions[model] = entries
        actual = base.metrics(y, np.array([int(r["predicted_label"]) for r in entries]), np.array([float(r["probability_active"]) for r in entries]))
        if any(not np.isclose(value, report["models"][model]["pooled_oof"][metric], atol=1e-12) for metric, value in actual.items()):
            raise ValueError("OOF predictions do not reproduce the recorded pooled metrics")
        if model != "dummy" and bundle["models"][model].named_steps["classifier"].C != expected_c:
            raise ValueError("Fitted estimator C differs from the protocol")
    return dict(predictions), fold_index, bundle, report


def ranking_diagnostics(rows, predictions, folds, draws=2000, progress=False):
    y, groups = compare.labels_groups(rows)
    shared_weights = compare.bootstrap_weights(groups, draws)
    metrics, fold_metrics, operating = [], [], []
    for model in EXPECTED_C:
        entries = predictions[model]
        probability = np.array([float(r["probability_active"]) for r in entries])
        predicted = np.array([int(r["predicted_label"]) for r in entries])
        point = base.metrics(y, predicted, probability)
        point["precision"] = point["tp"] / (point["tp"] + point["fp"]) if point["tp"] + point["fp"] else None
        operating.append({"model": model, **point})
        for mode in POLICY["ranking_modes"]:
            indices = np.arange(len(rows)) if mode == "all_identities" else family_representatives(rows, probability)
            values = precision_at_budgets(y[indices], probability[indices])
            samples = {k: [] for k in KS}
            for weights in shared_weights:
                sample = precision_at_budgets(y[indices], probability[indices], weights=weights[indices])
                for k in KS:
                    if sample[k] is not None:
                        samples[k].append(sample[k])
            for k in KS:
                bounds = np.quantile(samples[k], [0.025, 0.975]) if len(samples[k]) >= draws * 0.9 else (None, None)
                metrics.append({"model": model, "ranking_mode": mode, "k": k, "candidate_pool": len(indices),
                                "expected_precision_at_k": values[k], "conditional_low": bounds[0], "conditional_high": bounds[1],
                                "valid_bootstrap_draws": len(samples[k]), "requested_draws": draws})
        for fold in range(1, 6):
            indices = np.array([i for i in range(len(rows)) if folds[i] == fold])
            for mode in POLICY["ranking_modes"]:
                if mode != "all_identities":
                    reduced = family_representatives([rows[i] for i in indices], probability[indices])
                    indices_mode = indices[reduced]
                else:
                    indices_mode = indices
                values = precision_at_budgets(y[indices_mode], probability[indices_mode])
                fold_metrics.extend({"model": model, "fold": fold, "ranking_mode": mode, "k": k,
                                     "candidate_pool": len(indices_mode), "expected_precision_at_k": values[k]} for k in KS)
        if progress:
            print(f"Training OOF ranking diagnostics complete: {model}", flush=True)
    return metrics, fold_metrics, operating


def protocol_markdown(report):
    spec = report["spec"]
    counts = report["split_counts"]
    subsets = report["neighbor_summary"]["subset_counts"]
    return f"""# Pre-test evaluation protocol

Prepared at UTC: {report['prepared_at_utc']}

This is a development-informed decision recorded before held-out outcomes.
It is not a claim of preregistration before development. Commit this file and
`evaluation_manifest.json` before opening any held-out outcome table.

## Primary model and analysis

- Model: corrected charge plus C-terminal amidation (`{PRIMARY}`).
- Charge: count(K,R) - count(D,E) + amidation; features ordered as charge, amidation.
- L2 logistic regression: C=0.1, prediction threshold 0.5 (exact tie inactive).
- Primary metric: balanced accuracy on the original overall primary test set.
- Secondary: AUROC, average precision, MCC, accuracy, precision, sensitivity,
  specificity, confusion counts, and tie-aware precision@5/10/20.
- Original test: {counts['test']['examples']} identities, {counts['test']['groups']} families;
  class counts are recorded in the existing split manifest. No test predictions
  or outcome-level analysis were generated here.
- Evaluate the five other fixed feature-ladder models once as secondary
  comparisons. Use their already-selected Cs; select nothing from test scores.
- Keep the clinical set exploratory and the clinical-family buffer out of training.

## Fixed artifacts

- Primary model SHA-256: `{report['primary_model_sha256']}`
- Source model bundle SHA-256: `{spec['source_sha256']['models.joblib']}`
- Split manifest SHA-256: `{spec['split_manifest_sha256']}`
- Comparison report SHA-256: `{spec['comparison_report_sha256']}`
- Training table SHA-256: `{spec['train_sha256']}`
- Complete settings, source-code hashes and dependency versions: `evaluation_manifest.json`.

## Sequence-only novelty subsets

All test identities remain in the primary overall analysis. The following are
secondary subsets, defined without reading held-out outcome CSVs:

| Subset | Rule | Identities |
|---|---|---:|
| global60_strict | Maximum direct global edit similarity to training <0.60 | {subsets['global60_strict']} |
| local60_strict | No detected qualifying local identity >=0.60; complete optimal-tie scans | {subsets['local60_strict']} |
| combined60_strict | Intersection of the preceding subsets | {subsets['combined60_strict']} |

Local alignment is Smith-Waterman with BLOSUM62, gap open -10 and extension -0.5.
Identity divides exact residue matches by all alignment columns including gaps.
Qualification requires at least 10 aligned residue pairs and coverage >=80%
of the shorter peptide (paired residues divided by shorter length). Both peptide
coverages and aligned lengths are exported. A short fragment can therefore
qualify despite low coverage of its longer neighbor.

For each pair, scan up to 100 score-optimal alignments. Any incomplete tie scan
excludes that test identity from local/combined strict subsets. Maxima refer to
these score-optimal local alignments, not every possible suboptimal alignment.
Absence of a detected neighbor does not establish biological independence.

Keep the original split fixed. The 60% direct-neighbor screen differs from
transitive 60% connected-component linkage and from alignment percent identity.

## Uncertainty, small subsets and interpretation

Final overall intervals will use 5,000 paired whole-family bootstrap draws,
seed 42, percentile 95% intervals on fixed predictions. Disclose one-class
resample omissions; suppress an interval if more than 10% are omitted.
For AUROC intervals additionally require >=10 families and >=5 examples per
class; otherwise report the support limitation. AUROC is undefined with one
class. Always report subset n, class counts and family counts; never hide an
empty or small subset. The clinical set has only three inactive identities.

The {counts['test']['examples']}-identity primary test provides a limited final check, not a powered
ranking of small feature-set differences. No fixed interval width is promised.
Report all declared outcomes and retain a disappointing result without retuning.

## OOF prioritization diagnostics

Precision@5/10/20 budgets are illustrative and fixed before inspecting P@k;
they are not established laboratory budgets. Round scores to 12 decimal places
and average over random boundary-tie selection. Family-diverse ranking chooses
one highest-scoring identity per original family, breaking within-family ties
by sequence+amidation lexicographic order without labels.

OOF intervals use paired whole-family resampling of fixed predictions and
are conditional development diagnostics. They omit refitting/selection
uncertainty. A resample can contain replicated family copies; it is not an
actual recommendation to synthesize duplicates. Pooled OOF scores come from
different fitted models and are not guaranteed to be mutually calibrated;
fold-specific P@k is also exported. Database P@k does not establish prospective
success for generated peptides or a different candidate population.

## Release conditions

Review the metadata-only diagnostics, then commit this protocol, the membership
table and manifest before test outcomes are opened. Any change to model,
metric, alignment or subset rules requires a dated amendment before outcome
inspection. Model binaries remain generated local artifacts; their hashes
are committed. This script does not perform final test evaluation.

Derived from DBAASP data, https://dbaasp.org/; see project citation and data terms.
"""


def prepare(split_dir, comparison_dir, results_root=None, draws=2000, progress=False):
    split_dir, comparison_dir = Path(split_dir).resolve(), Path(comparison_dir).resolve()
    rows, split = base.load_training(split_dir)
    predictions, folds, bundle, comparison_report = load_comparison(split_dir, comparison_dir, rows)
    module, group_path = compare.grouping_module()
    spec = {"version": VERSION, "policy": POLICY, "train_sha256": base.sha_file(split_dir / "train.csv"),
            "split_manifest_sha256": base.sha_file(split_dir / "split_manifest.json"),
            "comparison_report_sha256": base.sha_file(comparison_dir / "comparison_report.json"),
            "source_sha256": {name: base.sha_file(comparison_dir / name) for name in ("models.joblib", "selected_oof_predictions.csv", "cv_folds.csv")},
            "code_sha256": {"prepare_evaluation.py": base.sha_file(Path(__file__)),
                            "train_baseline.py": base.sha_file(Path(base.__file__)),
                            "compare_features.py": base.sha_file(Path(compare.__file__)),
                            "prepare_splits.py": base.sha_file(group_path)},
            "environment": {**base.versions(), "biopython": Bio.__version__},
            "development_bootstrap_draws": draws, "development_bootstrap_seed": 42}
    digest = hashlib.sha256(base.canonical(spec).encode()).hexdigest()[:16]
    root = Path(results_root) if results_root else comparison_dir.parent / "pretest"
    output = root / f"{VERSION}-{digest}"
    if output.exists():
        report = json.loads((output / "evaluation_manifest.json").read_text(encoding="utf-8"))
        if report["spec"] != spec:
            raise ValueError("Existing pre-test specification differs")
        for name, expected in report["output_sha256"].items():
            if base.sha_file(output / name) != expected:
                raise ValueError(f"Existing pre-test output changed: {name}")
        return output, report
    assignments = base.read_csv(split_dir / "split_assignments.csv")
    neighbors = neighbor_audit(rows, assignments, progress)
    if (len(neighbors) != split["counts"]["test"]["examples"]
            or len({r['sequence'] for r in neighbors}) != split["counts"]["test"]["unique_sequences"]
            or len({r['group_id'] for r in neighbors}) != split["counts"]["test"]["groups"]):
        raise ValueError("Test sequence assignment metadata disagree with the fixed split counts")
    ranking, fold_ranking, operating = ranking_diagnostics(rows, predictions, folds, draws, progress)
    neighbor_summary = {
        "test_identities": len(neighbors), "training_unique_sequences": len({r['sequence'] for r in rows}),
        "subset_counts": {name: sum(r[name] for r in neighbors) for name in ("global60_strict", "local60_strict", "combined60_strict")},
        "global_similarity_histogram": identity_histogram([r["maximum_global_edit_similarity"] for r in neighbors]),
        "unrestricted_local_identity_histogram": identity_histogram([r["unrestricted_local_identity"] for r in neighbors]),
        "qualified_local_identity_histogram": identity_histogram([r["qualified_local_identity"] for r in neighbors]),
        "test_identities_with_incomplete_optimal_tie_scans": sum(not r["optimal_alignment_scan_complete"] for r in neighbors),
        "subset_class_counts": "not computed: held-out outcomes remain sealed"}
    payloads = {"test_neighbor_audit.csv": base.csv_text(neighbors),
                "test_subset_membership.csv": base.csv_text([{k: r[k] for k in ("sequence", "c_amidated", "group_id", "global60_strict", "local60_strict", "combined60_strict")} for r in neighbors]),
                "oof_precision_at_k.csv": base.csv_text(ranking),
                "oof_precision_at_k_by_fold.csv": base.csv_text(fold_ranking),
                "oof_operating_points.csv": base.csv_text(operating),
                "neighbor_summary.json": json.dumps(neighbor_summary, indent=2) + "\n"}
    report = {"spec": spec, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
              "split_counts": split["counts"], "neighbor_summary": neighbor_summary,
              "heldout_csvs_opened": [], "final_test_evaluated": False,
              "protocol_commit_required_before_test": True,
              "notes": ["No fitted model was changed or selected from held-out outcomes.",
                        "This is a development-informed pre-test protocol, not preregistration before training.",
                        "Local similarity is a diagnostic under declared scoring/coverage rules, not proof of homology or independence.",
                        "OOF precision intervals are conditional, selection-affected database diagnostics, not prospective wet-lab success rates."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".pretest-preparing-", dir=output.parent))
    try:
        for name, content in payloads.items():
            (temporary / name).write_text(content, encoding="utf-8", newline="")
        joblib.dump({"pipeline": bundle["models"][PRIMARY], "feature_names": POLICY["primary_features"],
                     "source_bundle_sha256": spec["source_sha256"]["models.joblib"], "policy": POLICY}, temporary / "primary_model.joblib")
        report["primary_model_sha256"] = base.sha_file(temporary / "primary_model.joblib")
        (temporary / "evaluation_protocol.md").write_text(protocol_markdown(report), encoding="utf-8")
        report["output_sha256"] = {p.name: base.sha_file(p) for p in sorted(temporary.iterdir())}
        (temporary / "evaluation_manifest.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.rename(output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_dir", type=Path)
    parser.add_argument("--comparison-dir", required=True, type=Path)
    parser.add_argument("--results-root", type=Path)
    args = parser.parse_args()
    output, report = prepare(args.split_dir, args.comparison_dir, args.results_root, progress=True)
    print("Primary model:", PRIMARY, "C=0.1; prediction threshold=0.5")
    print("Sequence-only novelty subset counts:", json.dumps(report["neighbor_summary"]["subset_counts"], sort_keys=True))
    print("Results:", output)
    print("Commit evaluation_protocol.md, evaluation_manifest.json and the subset membership before opening test outcomes.")
    print("No held-out CSVs were opened; no final test or clinical metrics were computed.")


if __name__ == "__main__":
    main()
