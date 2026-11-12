"""Training-only ESM-2 comparison under docs/embedding_v2_protocol.json.

The encoder is frozen. Every classifier/scaler fits within a development fold.
No v1 test or clinical outcomes are read and the shipped predictor is unchanged.
"""

import argparse
import hashlib
import json
import shutil
import tempfile
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.exceptions import ConvergenceWarning

import compare_features as comparison
import review_training_notes as review
import train_baseline as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("balanced_accuracy", "roc_auc", "average_precision", "accuracy", "sensitivity", "specificity")
ARMS = ("charge_original_training", "charge_expanded_training", "esm2_expanded_training")


def mean_residues(hidden, mask):
    """Mean only residue tokens, excluding start/end/padding tokens."""
    hidden, mask = np.asarray(hidden), np.asarray(mask, dtype=bool)
    if hidden.ndim != 3 or mask.shape != hidden.shape[:2] or np.any(mask.sum(axis=1) == 0):
        raise ValueError("Invalid residue mask or empty sequence")
    if not np.isfinite(hidden).all():
        raise ValueError("Nonfinite hidden representation")
    return (hidden * mask[:, :, None]).sum(axis=1) / mask.sum(axis=1)[:, None]


def load_candidates(protocol):
    dataset = ROOT / protocol["dataset"]
    report_path = ROOT / protocol["review_report"]
    if dataset.name != "training_candidates.csv" or report_path.name != "review_report.json":
        raise ValueError("Protocol must use the training review outputs")
    if base.sha_file(dataset) != protocol["dataset_sha256"] or base.sha_file(report_path) != protocol["review_report_sha256"]:
        raise ValueError("Candidate dataset/report fingerprint differs from protocol")
    report = json.loads(report_path.read_text())
    if report["output_sha256"]["training_candidates.csv"] != protocol["dataset_sha256"]:
        raise ValueError("Review report does not fingerprint the candidate dataset")
    original, _ = base.load_training(review.DEFAULT_SPLIT)
    originals = {base.identity(r): r for r in original}
    groups = {r["sequence"]: r for r in base.read_csv(review.DEFAULT_SPLIT / "sequence_groups.csv")}
    rows = base.read_csv(dataset)
    if len(rows) != protocol["candidate_training_identities"]:
        raise ValueError("Candidate identity count differs from protocol")
    seen = set()
    for row in rows:
        key = base.identity(row)
        if key in seen or row["partition"] != "train" or row["model_eligible"] != "True":
            raise ValueError("Duplicate, non-training or ineligible candidate")
        seen.add(key)
        if row["original_training_identity"] not in {"True", "False"}:
            raise ValueError("Invalid original-training flag")
        if row["label"] not in {"0", "1"} or row["status"] != ("active" if row["label"] == "1" else "inactive"):
            raise ValueError("Invalid candidate label/status")
        group = groups.get(row["sequence"], {})
        if group.get("group_partition") != "train" or group.get("group_id") != row["group_id"]:
            raise ValueError("Candidate is outside frozen training families")
        if (key in originals) != (row["original_training_identity"] == "True"):
            raise ValueError("Candidate original-training membership changed")
        if key in originals and row["label"] != originals[key]["label"]:
            raise ValueError("Original training label changed")
        base.feature_vector(row["sequence"], row["c_amidated"])
    if not set(originals) <= seen or len(originals) != protocol["original_training_identities"]:
        raise ValueError("Original training identities are missing")
    return rows


def extract_embeddings(rows, protocol, cache_root):
    sequences = sorted({r["sequence"] for r in rows})
    spec = {"sequences": sequences, "model_id": protocol["model_id"],
            "revision": protocol["model_revision"], "pooling": protocol["representation"]}
    key = hashlib.sha256(base.canonical(spec).encode()).hexdigest()[:16]
    cache = Path(cache_root) / key
    if cache.exists():
        manifest = json.loads((cache / "embedding_manifest.json").read_text())
        if manifest["spec"] != spec or base.sha_file(cache / "embeddings.npz") != manifest["embeddings_sha256"]:
            raise ValueError("Embedding cache changed")
        with np.load(cache / "embeddings.npz", allow_pickle=False) as saved:
            vectors = saved["vectors"]
            if saved["sequences"].tolist() != sequences:
                raise ValueError("Embedding cache sequence order changed")
        if vectors.shape != (len(sequences), 320) or not np.isfinite(vectors).all():
            raise ValueError("Invalid cached embedding dimensions/values")
        return dict(zip(sequences, vectors)), manifest

    # Optional imports: ordinary offline test suites do not download weights or
    # require PyTorch. Only this explicit command downloads the fixed model.
    import torch
    import transformers
    from huggingface_hub import snapshot_download
    from transformers import AutoTokenizer, EsmModel

    required = ["config.json", "model.safetensors", "tokenizer_config.json", "special_tokens_map.json", "vocab.txt"]
    snapshot = Path(snapshot_download(protocol["model_id"], revision=protocol["model_revision"], allow_patterns=required))
    tokenizer = AutoTokenizer.from_pretrained(snapshot, local_files_only=True, trust_remote_code=False)
    encoder = EsmModel.from_pretrained(snapshot, local_files_only=True, add_pooling_layer=False, use_safetensors=True)
    encoder.eval()
    torch.set_num_threads(2)
    vectors = []
    with torch.inference_mode():
        for start in range(0, len(sequences), 32):
            batch = sequences[start:start + 32]
            tokens = tokenizer([" ".join(s) for s in batch], padding=True,
                               return_special_tokens_mask=True, return_tensors="pt", truncation=False)
            special = tokens.pop("special_tokens_mask")
            mask = (tokens["attention_mask"].bool() & ~special.bool()).cpu().numpy()
            if mask.sum(axis=1).tolist() != list(map(len, batch)):
                raise ValueError("Tokenizer did not preserve exactly one token per residue")
            hidden = encoder(**tokens).last_hidden_state.cpu().numpy()
            vectors.append(mean_residues(hidden, mask).astype(np.float32))
            print(f"Embedded {min(start + 32, len(sequences))}/{len(sequences)} training sequences", flush=True)
    vectors = np.concatenate(vectors)
    if vectors.shape != (len(sequences), 320):
        raise ValueError("Unexpected ESM-2 representation size")
    manifest = {"spec": spec, "weights_sha256": {name: base.sha_file(snapshot / name) for name in required},
                "torch": torch.__version__, "transformers": transformers.__version__, "encoder_fine_tuned": False}
    cache.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".esm-cache-", dir=cache.parent))
    try:
        np.savez_compressed(temporary / "embeddings.npz", vectors=vectors, sequences=np.asarray(sequences))
        manifest["embeddings_sha256"] = base.sha_file(temporary / "embeddings.npz")
        (temporary / "embedding_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        temporary.replace(cache)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return dict(zip(sequences, vectors)), manifest


def feature_matrices(rows, embeddings):
    charge = comparison.matrices(rows)["charge_corrected_amidation"]
    protein = np.asarray([embeddings[r["sequence"]] for r in rows], dtype=float)
    if protein.ndim != 2 or not np.isfinite(protein).all():
        raise ValueError("Invalid embedding matrix")
    return charge, np.column_stack((protein, [int(r["c_amidated"] == "True") for r in rows]))


def fit_arm(x, y, folds, original, model, c, original_only=False):
    scores, predictions = [], []
    for number, (train, valid) in enumerate(folds, 1):
        fitting = train[original[train]] if original_only else train
        if set(y[fitting]) != {0, 1}:
            raise ValueError("Training fold lacks a class; do not search seeds")
        estimator = base.logistic_pipeline(c)
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            estimator.fit(x[fitting], y[fitting])
        probability = estimator.predict_proba(x[valid])[:, list(estimator.classes_).index(1)]
        predicted = (probability > 0.5).astype(int)
        scores.append({"model": model, "C": c, "fold": number, "fit_rows": len(fitting),
                       "validation_rows": len(valid), **base.metrics(y[valid], predicted, probability)})
        predictions.extend({"row_index": int(index), "model": model, "fold": number,
                            "true_label": int(y[index]), "predicted_label": int(pred),
                            "probability_active": float(prob)}
                           for index, pred, prob in zip(valid, predicted, probability))
    predictions.sort(key=lambda row: row["row_index"])
    if [r["row_index"] for r in predictions] != list(range(len(y))):
        raise ValueError("Validation coverage is not exactly once per identity")
    return scores, predictions


def compare(rows, embeddings, protocol):
    if protocol["cv"] != {"folds": 5, "seed": 42, "method": protocol["cv"]["method"]}:
        raise ValueError("This runner implements exactly five folds with seed 42")
    charge, protein = feature_matrices(rows, embeddings)
    y, groups = comparison.labels_groups(rows)
    original = np.asarray([r["original_training_identity"] == "True" for r in rows])
    folds = base.make_folds(y, groups)
    selected, all_scores, predictions = {}, [], {}
    for arm in ARMS:
        x = protein if arm.startswith("esm2") else charge
        candidates = []
        grid = protocol["embedding_C_grid"] if arm.startswith("esm2") else [protocol["baseline_C"]]
        for c in grid:
            scores, predicted = fit_arm(x, y, folds, original, arm, c, arm == ARMS[0])
            all_scores.extend(scores)
            summary = base.summary(scores)
            candidates.append((summary["balanced_accuracy"]["mean"], -c, c, summary, predicted))
        _, _, c, summary, predicted = max(candidates, key=lambda row: row[:2])
        probability = np.asarray([r["probability_active"] for r in predicted])
        labels = np.asarray([r["predicted_label"] for r in predicted])
        selected[arm] = {"C": c, "features": x.shape[1], "cv": summary, "pooled_oof": base.metrics(y, labels, probability)}
        predictions[arm] = predicted
        print(f"{arm}: C={c:g}; mean-fold balanced accuracy {summary['balanced_accuracy']['mean']:.3f}; "
              f"pooled OOF accuracy {selected[arm]['pooled_oof']['accuracy']:.3f}; AUROC {selected[arm]['pooled_oof']['roc_auc']:.3f}", flush=True)
    assignments = []
    for number, (_, valid) in enumerate(folds, 1):
        assignments.extend({"row_index": int(i), "sequence": rows[i]["sequence"], "c_amidated": rows[i]["c_amidated"],
                            "group_id": rows[i]["group_id"], "validation_fold": number,
                            "original_training_identity": rows[i]["original_training_identity"]} for i in valid)
    return selected, predictions, all_scores, sorted(assignments, key=lambda r: r["row_index"])


def intervals(rows, predictions, draws=2000):
    y, groups = comparison.labels_groups(rows)
    weights = comparison.bootstrap_weights(groups, draws, seed=42)
    valid = (weights[:, y == 1].sum(axis=1) > 0) & (weights[:, y == 0].sum(axis=1) > 0)
    weights = weights[valid]
    if len(weights) < draws * 0.9:
        raise ValueError("Insufficient two-class group-bootstrap draws")
    samples, estimates, output, deltas = {}, {}, [], []
    for model, records in predictions.items():
        pred = np.asarray([r["predicted_label"] for r in records])
        prob = np.asarray([r["probability_active"] for r in records])
        first = np.asarray([comparison.weighted_metrics(y, pred, prob, w) for w in weights])
        extra = np.column_stack((weights @ (y == pred) / weights.sum(axis=1),
                                 weights @ ((y == 1) & (pred == 1)) / (weights @ (y == 1)),
                                 weights @ ((y == 0) & (pred == 0)) / (weights @ (y == 0))))
        samples[model] = np.column_stack((first, extra))
        estimates[model] = base.metrics(y, pred, prob)
        for index, metric in enumerate(METRICS):
            low, high = np.quantile(samples[model][:, index], [0.025, 0.975])
            output.append({"model": model, "metric": metric, "estimate": estimates[model][metric],
                           "low": float(low), "high": float(high)})
    for before, after in zip(ARMS, ARMS[1:]):
        for index, metric in enumerate(METRICS):
            low, high = np.quantile(samples[after][:, index] - samples[before][:, index], [0.025, 0.975])
            deltas.append({"model": after, "previous_model": before, "metric": metric,
                           "estimate": estimates[after][metric] - estimates[before][metric],
                           "low": float(low), "high": float(high)})
    return output, deltas, {"requested_draws": draws, "valid_shared_draws": len(weights), "seed": 42}


def run(protocol_path, cache_root, results_root):
    protocol = json.loads(Path(protocol_path).read_text())
    rows = load_candidates(protocol)
    embeddings, manifest = extract_embeddings(rows, protocol, cache_root)
    spec = {"protocol_sha256": base.sha_file(protocol_path), "runner_sha256": base.sha_file(__file__),
            "baseline_helper_sha256": base.sha_file(base.__file__),
            "comparison_helper_sha256": base.sha_file(comparison.__file__),
            "review_helper_sha256": base.sha_file(review.__file__), "environment": base.versions(),
            "embedding_manifest_sha256": hashlib.sha256(base.canonical(manifest).encode()).hexdigest()}
    digest = hashlib.sha256(base.canonical(spec).encode()).hexdigest()[:16]
    output = Path(results_root) / f"{protocol['version']}-{digest}"
    if output.exists():
        report = json.loads((output / "comparison_report.json").read_text())
        if report["spec"] != spec:
            raise ValueError("Existing comparison specification changed")
        for name, expected in report["output_sha256"].items():
            if base.sha_file(output / name) != expected:
                raise ValueError("Existing comparison output changed: " + name)
        return output, report
    selected, predictions, scores, folds = compare(rows, embeddings, protocol)
    ci, deltas, bootstrap = intervals(rows, predictions, protocol["bootstrap"]["draws"])
    report = {"spec": spec, "protocol": protocol, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
              "counts": {"identities": len(rows), "sequences": len(embeddings), "groups": len({r['group_id'] for r in rows}),
                         "active": sum(r['label'] == '1' for r in rows)},
              "selected": selected, "bootstrap": bootstrap, "heldout_classifier_csvs_opened": [],
              "encoder_fine_tuned": False, "v1_model_changed": False,
              "scope": "Selected training-only development predictions; no final test improvement estimated"}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".esm-comparison-", dir=output.parent))
    try:
        tables = {"cv_scores.csv": scores, "cv_folds.csv": folds,
                  "selected_oof_predictions.csv": [r for arm in ARMS for r in predictions[arm]],
                  "intervals.csv": ci, "paired_differences.csv": deltas}
        for name, table in tables.items():
            (temporary / name).write_text(review.csv_text(table, list(table[0])), encoding="utf-8")
        (temporary / "embedding_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        report["output_sha256"] = {p.name: base.sha_file(p) for p in temporary.iterdir()}
        (temporary / "comparison_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        temporary.replace(output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=ROOT / "docs/embedding_v2_protocol.json")
    parser.add_argument("--cache-root", type=Path, default=ROOT / "data/processed/embedding_cache")
    parser.add_argument("--results-root", type=Path, default=ROOT / "results/development")
    args = parser.parse_args()
    output, report = run(args.protocol, args.cache_root, args.results_root)
    print("Results:", output)
    print("Training-only development comparison; v1 held-out accuracy is unchanged.")


if __name__ == "__main__":
    main()
