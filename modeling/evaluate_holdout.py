"""Evaluate the committed pre-test protocol without fitting or selecting models.

Usage: python modeling/evaluate_holdout.py SPLIT_DIR --pretest-dir PRETEST_DIR
       --comparison-dir COMPARISON_DIR --results-root results/evaluation

The protocol, manifest, membership table and this evaluator must be committed
and unchanged. Inputs/models/source/dependencies are verified before outcome
CSVs are opened. A protocol-specific run marker prevents silent fresh reruns.
--resume permits an interrupted run only with the identical specification.
"""

import argparse
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import Bio
import joblib
import numpy as np

import train_baseline as base
import compare_features as compare
import prepare_evaluation as prep

VERSION = "holdout-evaluation-v1"
METRICS = ("balanced_accuracy", "roc_auc", "average_precision", "mcc", "accuracy", "precision", "sensitivity", "specificity")


def git_evidence(pretest_dir):
    """Verify exact committed bytes; new evaluator source must be committed too."""
    result = subprocess.run(["git", "-C", str(pretest_dir), "rev-parse", "--show-toplevel"], capture_output=True, check=True)
    root = Path(os.fsdecode(result.stdout).strip()).resolve()
    required = [Path(pretest_dir) / n for n in ("evaluation_protocol.md", "evaluation_manifest.json", "test_subset_membership.csv")]
    required.append(Path(__file__).resolve())
    verified = []
    for path in required:
        relative = path.resolve().relative_to(root).as_posix()
        committed = subprocess.run(["git", "-C", str(root), "show", f"HEAD:{relative}"], capture_output=True)
        if committed.returncode or committed.stdout != path.read_bytes():
            raise ValueError(f"Commit this file unchanged before evaluating outcomes: {relative}")
        verified.append(relative)
    commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, check=True)
    return {"commit": os.fsdecode(commit.stdout).strip(), "verified_committed_files": verified}


def verify_pretest(split_dir, pretest_dir, comparison_dir):
    report = json.loads((pretest_dir / "evaluation_manifest.json").read_text(encoding="utf-8"))
    spec = report["spec"]
    if spec["policy"] != prep.POLICY:
        raise ValueError("Current policy differs from the frozen pre-test policy")
    for name, expected in report["output_sha256"].items():
        if Path(name).name != name or base.sha_file(pretest_dir / name) != expected:
            raise ValueError(f"Frozen pre-test artifact changed: {name}")
    if base.sha_file(split_dir / "split_manifest.json") != spec["split_manifest_sha256"]:
        raise ValueError("Split manifest differs from the pre-test snapshot")
    if base.sha_file(split_dir / "train.csv") != spec["train_sha256"]:
        raise ValueError("Training table differs from the pre-test snapshot")
    if base.sha_file(comparison_dir / "comparison_report.json") != spec["comparison_report_sha256"]:
        raise ValueError("Comparison report changed")
    for name, expected in spec["source_sha256"].items():
        if Path(name).name != name or base.sha_file(comparison_dir / name) != expected:
            raise ValueError(f"Frozen comparison artifact changed: {name}")
    _, grouping_path = compare.grouping_module()
    sources = {"prepare_evaluation.py": Path(prep.__file__), "train_baseline.py": Path(base.__file__),
               "compare_features.py": Path(compare.__file__), "prepare_splits.py": grouping_path}
    if any(base.sha_file(path) != spec["code_sha256"][name] for name, path in sources.items()):
        raise ValueError("Feature/protocol source changed since freezing; preserve the frozen source")
    if {**base.versions(), "biopython": Bio.__version__} != spec["environment"]:
        raise ValueError("Environment differs from the frozen report; use the original environment")
    split = json.loads((split_dir / "split_manifest.json").read_text(encoding="utf-8"))
    if split["counts"] != report["split_counts"]:
        raise ValueError("Split counts disagree with the frozen report")
    return report, split


def read_frozen_outcomes(path, expected_hash, role, assignments, counts):
    """Called only after metadata, commits and run marker are verified."""
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected_hash:
        raise ValueError(f"Frozen outcome table changed: {path.name}")
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8"), newline=""))
    rows = list(reader)
    required = {"sequence", "c_amidated", "label", "status", "model_eligible", "partition", "group_id"}
    if not rows or any(not required <= set(r) or None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("Empty/malformed outcome table")
    index = {base.identity(r): r for r in assignments if r["partition"] == role}
    if len({base.identity(r) for r in rows}) != len(rows) or {base.identity(r) for r in rows} != set(index):
        raise ValueError("Outcome identities differ from frozen assignments")
    for row in rows:
        base.feature_vector(row["sequence"], row["c_amidated"])
        if (row["partition"] != role or row["model_eligible"] != "True" or row["label"] not in {"0", "1"}
                or row["status"] != ("active" if row["label"] == "1" else "inactive")
                or row["group_id"] != index[base.identity(row)]["group_id"]):
            raise ValueError("Invalid outcome identity, class, role or group")
    actual = {"examples": len(rows), "unique_sequences": len({r['sequence'] for r in rows}),
              "groups": len({r['group_id'] for r in rows}), "active": sum(r['label'] == '1' for r in rows),
              "inactive": sum(r['label'] == '0' for r in rows)}
    if actual != counts:
        raise ValueError("Outcome counts differ from the frozen split")
    return rows


def weighted_metrics(y, predicted, probability, weights):
    """Vectorized weighted binary metrics; AUROC/AP handle tied scores exactly."""
    weights = np.atleast_2d(weights).astype(float)
    positive, negative = y == 1, y == 0
    tp, fn = weights @ (positive & (predicted == 1)), weights @ (positive & (predicted == 0))
    tn, fp = weights @ (negative & (predicted == 0)), weights @ (negative & (predicted == 1))
    divide = lambda a, b: np.divide(a, b, out=np.full_like(a, np.nan, dtype=float), where=b > 0)
    sensitivity, specificity = divide(tp, tp + fn), divide(tn, tn + fp)
    precision = np.nan_to_num(divide(tp, tp + fp), nan=0.0)
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = np.nan_to_num(divide(tp * tn - fp * fn, denom), nan=0.0)
    _, inverse = np.unique(probability, return_inverse=True)
    indicator = np.eye(int(inverse.max()) + 1)[inverse]
    pos = weights @ (indicator * positive[:, None])
    neg = weights @ (indicator * negative[:, None])
    before_neg = np.cumsum(neg, axis=1) - neg
    auc = divide(np.sum(pos * (before_neg + 0.5 * neg), axis=1), pos.sum(axis=1) * neg.sum(axis=1))
    descending_pos, descending_neg = pos[:, ::-1], neg[:, ::-1]
    cumulative_pos = np.cumsum(descending_pos, axis=1)
    cumulative_all = cumulative_pos + np.cumsum(descending_neg, axis=1)
    precision_at_block = np.divide(cumulative_pos, cumulative_all, out=np.zeros_like(cumulative_pos), where=cumulative_all > 0)
    ap = divide(np.sum(descending_pos * precision_at_block, axis=1), pos.sum(axis=1))
    return {"balanced_accuracy": (sensitivity + specificity) / 2, "roc_auc": auc, "average_precision": ap,
            "mcc": mcc, "accuracy": divide(tp + tn, tp + tn + fp + fn), "precision": precision,
            "sensitivity": sensitivity, "specificity": specificity, "tn": tn, "fp": fp, "fn": fn, "tp": tp}


def nullable(value):
    return float(value) if np.isfinite(value) else None


def panel_results(name, rows, probabilities, predictions, draws, seed):
    counts = {"n": len(rows), "active": sum(r["label"] == "1" for r in rows),
              "inactive": sum(r["label"] == "0" for r in rows), "groups": len({r["group_id"] for r in rows})}
    metrics, intervals, ranking = [], [], []
    if not rows:
        for model in prep.EXPECTED_C:
            metrics.append({"panel": name, "model": model, **counts, **{k: None for k in (*METRICS, "tn", "fp", "fn", "tp")}})
        return metrics, intervals, ranking, {**counts, "status": "empty subset; metrics undefined"}
    y, groups = compare.labels_groups(rows)
    weights = compare.bootstrap_weights(groups, draws, seed) if len(set(groups)) > 1 else np.ones((draws, len(y)))
    valid = (weights[:, y == 1].sum(axis=1) > 0) & (weights[:, y == 0].sum(axis=1) > 0)
    valid_weights = weights[valid]
    support = {**counts, "status": "evaluated fixed predictions", "bootstrap_draws": draws,
               "two_class_draws": int(valid.sum()), "one_class_draws_omitted": int((~valid).sum())}
    bootstrap_reason = ("fewer than two sampling groups" if counts["groups"] < 2 else
                        "more than 10% one-class draws omitted" if valid.sum() < 0.9 * draws else None)
    enough_draws = bootstrap_reason is None
    for model in prep.EXPECTED_C:
        probability, predicted = probabilities[model], predictions[model]
        point = weighted_metrics(y, predicted, probability, np.ones(len(y)))
        metrics.append({"panel": name, "model": model, **counts, **{key: nullable(value[0]) for key, value in point.items()}})
        sampled = weighted_metrics(y, predicted, probability, valid_weights) if len(valid_weights) else None
        for metric in METRICS:
            reason = None
            if not np.isfinite(point[metric][0]):
                reason = "metric undefined for this panel"
            elif not enough_draws:
                reason = bootstrap_reason
            elif metric == "roc_auc" and (counts["groups"] < 10 or min(counts["active"], counts["inactive"]) < 5):
                reason = "AUROC interval support below frozen minimum (10 groups, 5 examples/class)"
            bounds = np.quantile(sampled[metric], [0.025, 0.975]) if reason is None else (None, None)
            intervals.append({"panel": name, "model": model, "metric": metric, "estimate": nullable(point[metric][0]),
                              "low": bounds[0], "high": bounds[1], "interval_suppressed_reason": reason,
                              "valid_draws": int(valid.sum()), "requested_draws": draws})
        for mode in prep.POLICY["ranking_modes"]:
            idx = np.arange(len(rows)) if mode == "all_identities" else prep.family_representatives(rows, probability)
            points = prep.precision_at_budgets(y[idx], probability[idx])
            samples = {k: [] for k in prep.KS}
            if enough_draws:
                for w in valid_weights:
                    values = prep.precision_at_budgets(y[idx], probability[idx], weights=w[idx])
                    for k in prep.KS:
                        if values[k] is not None:
                            samples[k].append(values[k])
            for k in prep.KS:
                if points[k] is None:
                    reason = "k exceeds candidate pool"
                elif not enough_draws:
                    reason = bootstrap_reason
                elif len(samples[k]) < 0.9 * draws:
                    reason = "too few valid ranking resamples"
                else:
                    reason = None
                bounds = np.quantile(samples[k], [0.025, 0.975]) if reason is None else (None, None)
                ranking.append({"panel": name, "model": model, "ranking_mode": mode, "k": k, "candidate_pool": len(idx),
                                "expected_precision_at_k": points[k], "low": bounds[0], "high": bounds[1],
                                "valid_draws": len(samples[k]), "interval_suppressed_reason": reason})
    return metrics, intervals, ranking, support


def evaluate(split_dir, pretest_dir, comparison_dir, results_root, resume=False, progress=False):
    split_dir, pretest_dir, comparison_dir = (Path(p).resolve() for p in (split_dir, pretest_dir, comparison_dir))
    frozen, split = verify_pretest(split_dir, pretest_dir, comparison_dir)
    evidence = git_evidence(pretest_dir)
    manifest_hash = base.sha_file(pretest_dir / "evaluation_manifest.json")
    output = Path(results_root).resolve() / f"{VERSION}-{manifest_hash[:16]}"
    spec = {"version": VERSION, "pretest_manifest_sha256": manifest_hash,
            "evaluator_sha256": base.sha_file(Path(__file__)), "environment": base.versions(),
            "frozen_outcome_sha256": {n: split["output_sha256"][n] for n in ("test.csv", "clinical_exploratory.csv")},
            "conventions": {"undefined_precision": "0 if no predicted positive", "undefined_mcc": "0 if denominator zero",
                            "one_class_panels": "report undefined metrics as null; retain counts and all panels",
                            "comparisons": "all fixed models reported; no winner selected from outcomes"}}
    if output.exists():
        report = json.loads((output / "evaluation_report.json").read_text(encoding="utf-8"))
        if report["spec"] != spec:
            raise ValueError("Existing evaluation has different code/settings; preserve its results without rerunning")
        for name, expected in report["output_sha256"].items():
            if base.sha_file(output / name) != expected:
                raise ValueError(f"Saved evaluation artifact changed: {name}")
        return output, report
    assignments = base.read_csv(split_dir / "split_assignments.csv")
    membership = base.read_csv(pretest_dir / "test_subset_membership.csv")
    index = {base.identity(r): r for r in membership}
    expected_test = {base.identity(r): r for r in assignments if r["partition"] == "test"}
    if len(index) != len(membership) or set(index) != set(expected_test):
        raise ValueError("Frozen subset membership does not match test assignments")
    for key, row in index.items():
        names = ("global60_strict", "local60_strict", "combined60_strict")
        if (row["group_id"] != expected_test[key]["group_id"] or any(row[k] not in {"True", "False"} for k in names)
                or (row["combined60_strict"] == "True") != (row["global60_strict"] == row["local60_strict"] == "True")):
            raise ValueError("Invalid subset flags or groups")
    primary = joblib.load(pretest_dir / "primary_model.joblib")
    bundle = joblib.load(comparison_dir / "models.joblib")
    if primary["feature_names"] != prep.POLICY["primary_features"] or primary["policy"] != prep.POLICY:
        raise ValueError("Primary model metadata differ from the frozen policy")
    models = {**bundle["models"], prep.PRIMARY: primary["pipeline"]}
    for model, c in prep.EXPECTED_C.items():
        if model not in models or (c is not None and models[model].named_steps["classifier"].C != c):
            raise ValueError("Model or C does not match frozen list")
    marker = pretest_dir / "evaluation_started.json"
    marker_info = {"spec": spec, "output": str(output), "started_at_utc": datetime.now(timezone.utc).isoformat(), "git": evidence}
    try:
        with marker.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(marker_info, indent=2))
    except FileExistsError:
        previous = json.loads(marker.read_text(encoding="utf-8"))
        if previous["spec"] != spec or previous["output"] != str(output):
            raise ValueError("This protocol already has an evaluation attempt with different code/output settings")
        if not resume:
            raise ValueError("An evaluation attempt is recorded. Use --resume only to recover an interrupted identical run")
        marker_info = previous
    all_metrics, all_intervals, all_ranking, prediction_rows, supports = [], [], [], [], {}
    settings = frozen["spec"]["policy"]["final_intervals"]
    for role, filename in (("test", "test.csv"), ("clinical_exploratory", "clinical_exploratory.csv")):
        rows = read_frozen_outcomes(split_dir / filename, split["output_sha256"][filename], role, assignments, split["counts"][role])
        x = compare.matrices(rows)
        probabilities, predictions = {}, {}
        for model in prep.EXPECTED_C:
            fitted = models[model]
            probabilities[model] = fitted.predict_proba(x[model])[:, list(fitted.classes_).index(1)]
            # Explicit binary rule: logistic log-odds >0; dummy retains its fixed class.
            predictions[model] = fitted.predict(x[model])
            if not np.all(np.isfinite(probabilities[model])):
                raise ValueError("Nonfinite model probabilities")
            for row, predicted, probability in zip(rows, predictions[model], probabilities[model]):
                prediction_rows.append({"partition": role, "model": model, "sequence": row["sequence"],
                                        "c_amidated": row["c_amidated"], "group_id": row["group_id"],
                                        "true_label": row["label"], "predicted_label": int(predicted), "probability_active": float(probability)})
        panels = {"overall": np.arange(len(rows))} if role == "test" else {"clinical_exploratory": np.arange(len(rows))}
        if role == "test":
            for name in ("global60_strict", "local60_strict", "combined60_strict"):
                panels[name] = np.array([i for i, row in enumerate(rows) if index[base.identity(row)][name] == "True"], dtype=int)
        for name, idx in panels.items():
            values = panel_results(name, [rows[i] for i in idx], {m: p[idx] for m, p in probabilities.items()},
                                   {m: p[idx] for m, p in predictions.items()}, settings["draws"], settings["seed"])
            metrics, intervals, ranking, supports[name] = values
            all_metrics.extend(metrics); all_intervals.extend(intervals); all_ranking.extend(ranking)
            if progress:
                print(f"Evaluated {name}: n={supports[name]['n']}, groups={supports[name]['groups']}", flush=True)
    report = {"spec": spec, "started_at_utc": marker_info["started_at_utc"], "completed_at_utc": datetime.now(timezone.utc).isoformat(),
              "git": marker_info["git"], "primary_model": prep.PRIMARY, "primary_panel": "overall", "primary_metric": "balanced_accuracy",
              "primary_result": next(r for r in all_metrics if r["panel"] == "overall" and r["model"] == prep.PRIMARY),
              "panels": supports, "outcome_files_opened": ["test.csv", "clinical_exploratory.csv"],
              "model_fitting_performed": False, "model_selection_performed": False,
              "notes": ["Primary result is the fixed model on the original overall test.",
                        "All other models/subsets are secondary; clinical results are exploratory.",
                        "Bootstrap intervals condition on fitted models and this curated dataset; they are not prospective wet-lab guarantees.",
                        "Novelty subsets are diagnostic and do not prove biological independence.",
                        "No code, feature, threshold or grouping choice may be selected using these outcomes."],
              "source_attribution": "Derived from DBAASP data, https://dbaasp.org; see project citation and data terms."}
    payloads = {"metrics.csv": base.csv_text(all_metrics), "intervals.csv": base.csv_text(all_intervals),
                "precision_at_k.csv": base.csv_text(all_ranking), "predictions.csv": base.csv_text(prediction_rows)}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".evaluation-preparing-", dir=output.parent))
    try:
        for name, content in payloads.items():
            (temporary / name).write_text(content, encoding="utf-8", newline="")
        report["output_sha256"] = {p.name: base.sha_file(p) for p in sorted(temporary.iterdir())}
        (temporary / "evaluation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.rename(output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_dir", type=Path)
    parser.add_argument("--pretest-dir", required=True, type=Path)
    parser.add_argument("--comparison-dir", required=True, type=Path)
    parser.add_argument("--results-root", type=Path, default=Path("results/evaluation"))
    parser.add_argument("--resume", action="store_true", help="Recover an interrupted identical run; never changes frozen settings")
    args = parser.parse_args()
    output, report = evaluate(args.split_dir, args.pretest_dir, args.comparison_dir, args.results_root, args.resume, progress=True)
    print("PRIMARY OVERALL TEST RESULT:")
    print(json.dumps(report["primary_result"], indent=2))
    print("Results:", output)
    print("The declared holdouts have now been evaluated. Preserve all results; do not retune from test or clinical scores.")


if __name__ == "__main__":
    main()
