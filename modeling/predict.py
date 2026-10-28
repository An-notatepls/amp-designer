"""Score a peptide with the frozen primary model; never fit or select a model.

Run from the project root. Load only model bundles produced by this project
that you trust: joblib uses pickle and can execute code when loading.
"""
import argparse
import csv
import io
import json
from pathlib import Path

import joblib
import numpy as np

import train_baseline as base

DEFAULT_PRETEST = Path("results/pretest/pretest-v1-cbd14c16acf65e99")
FEATURES = ["charge_proxy", "c_amidated"]
MODEL = "charge_corrected_amidation"


def load_primary(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "evaluation_manifest.json").read_text(encoding="utf-8"))
    policy = manifest["spec"]["policy"]
    if (policy.get("primary_model") != MODEL
            or policy.get("primary_features") != FEATURES
            or policy.get("primary_C") != 0.1
            or policy.get("prediction_threshold") != 0.5):
        raise ValueError("Manifest differs from the supported frozen primary model")
    expected_env = manifest["spec"]["environment"]
    actual_env = base.versions()
    mismatches = {k: (expected_env.get(k), v) for k, v in actual_env.items()
                  if expected_env.get(k) != v}
    if mismatches:
        raise ValueError(f"Use the original model environment (expected, actual): {mismatches}")
    if base.sha_file(Path(base.__file__)) != manifest["spec"]["code_sha256"]["train_baseline.py"]:
        raise ValueError("Training feature source differs from the frozen source")
    model_path = directory / "primary_model.joblib"
    digest = base.sha_file(model_path)
    if digest != manifest["primary_model_sha256"] or digest != manifest["output_sha256"][model_path.name]:
        raise ValueError("Primary model fingerprint differs from the frozen manifest")
    bundle = joblib.load(model_path)
    if bundle.get("policy") != policy or bundle.get("feature_names") != FEATURES:
        raise ValueError("Saved model metadata differ from the frozen manifest")
    pipeline = bundle["pipeline"]
    if (pipeline.named_steps["classifier"].C != 0.1
            or list(pipeline.classes_) != [0, 1]
            or pipeline.n_features_in_ != 2):
        raise ValueError("Unexpected classifier configuration")
    return pipeline, digest


def prediction_features(rows):
    if not rows:
        raise ValueError("Provide at least one peptide")
    matrix = base.feature_matrix(rows)
    indices = [list(base.FEATURES).index(name) for name in FEATURES]
    return matrix[:, indices]


def score_rows(rows, pipeline, model_sha256):
    x = prediction_features(rows)
    probabilities = np.asarray(pipeline.predict_proba(x))[:, 1]
    # Use the fitted classifier's decision rule, including its inactive tie.
    labels = np.asarray(pipeline.predict(x))
    if not np.isfinite(probabilities).all() or ((probabilities < 0) | (probabilities > 1)).any():
        raise ValueError("Model returned invalid probabilities")
    return [{"sequence": row["sequence"], "c_amidated": row["c_amidated"],
             "charge_proxy": int(features[0]), "model_probability_active": float(probability),
             "predicted_label": "active" if int(label) == 1 else "inactive",
             "activity_definition": "recorded E. coli MIC <= 32 uM",
             "model": MODEL, "model_sha256": model_sha256}
            for row, features, probability, label in zip(rows, x, probabilities, labels)]


def input_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not {"sequence", "c_amidated"} <= set(names):
            raise ValueError("CSV needs unique column names including sequence,c_amidated")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed input CSV")
    return [{"sequence": row["sequence"], "c_amidated": row["c_amidated"]} for row in rows]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="Canonical uppercase sequence, 10-50 residues")
    source.add_argument("--input-csv", type=Path, help="Columns sequence,c_amidated; True/False only")
    chemistry = parser.add_mutually_exclusive_group()
    chemistry.add_argument("--amidated", action="store_true", help="C-terminal amidation")
    chemistry.add_argument("--unmodified", action="store_true", help="Unmodified termini")
    parser.add_argument("--pretest-dir", type=Path, default=DEFAULT_PRETEST)
    parser.add_argument("--output", type=Path, help="Write CSV to a new file; never overwrite")
    args = parser.parse_args(argv)
    if args.sequence is not None and not (args.amidated or args.unmodified):
        parser.error("Specify --amidated or --unmodified for a single sequence")
    if args.input_csv is not None and (args.amidated or args.unmodified):
        parser.error("For CSV input, specify chemistry in each row")
    try:
        rows = (input_rows(args.input_csv) if args.input_csv is not None else
                [{"sequence": args.sequence, "c_amidated": str(args.amidated)}])
        prediction_features(rows)  # Reject invalid input before loading the model.
        pipeline, digest = load_primary(args.pretest_dir)
        results = score_rows(rows, pipeline, digest)
        if args.output:
            buffer = io.StringIO(newline="")
            writer = csv.DictWriter(buffer, fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)
            with args.output.open("x", encoding="utf-8", newline="") as handle:
                handle.write(buffer.getvalue())
            print(f"Wrote {len(results)} predictions to {args.output}")
        else:
            print(json.dumps(results[0] if args.sequence is not None else results, indent=2))
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"Prediction failed: {error}\n")


if __name__ == "__main__":
    main()
