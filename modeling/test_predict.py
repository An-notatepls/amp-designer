"""Synthetic inference tests; no real holdout or database files are opened."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

import compare_features
import predict


def fixture(root):
    policy = {"primary_model": predict.MODEL, "primary_features": predict.FEATURES,
              "primary_C": 0.1, "prediction_threshold": 0.5}
    model = Pipeline([("scaler", StandardScaler()),
                      ("classifier", LogisticRegression(C=0.1))])
    model.fit([[-3, 0], [-2, 1], [2, 0], [3, 1]], [0, 0, 1, 1])
    path = root / "primary_model.joblib"
    joblib.dump({"pipeline": model, "feature_names": predict.FEATURES, "policy": policy}, path)
    digest = predict.base.sha_file(path)
    manifest = {"spec": {"policy": policy, "environment": predict.base.versions(),
                         "code_sha256": {"train_baseline.py": predict.base.sha_file(Path(predict.base.__file__))}},
                "primary_model_sha256": digest, "output_sha256": {path.name: digest}}
    (root / "evaluation_manifest.json").write_text(json.dumps(manifest))
    return model, manifest


class PredictionTests(unittest.TestCase):
    def test_features_match_training_and_amidation(self):
        rows = [{"sequence": "AKRDEAAAAA", "c_amidated": flag} for flag in ("False", "True")]
        actual = predict.prediction_features(rows)
        np.testing.assert_array_equal(actual, [[0, 0], [1, 1]])
        np.testing.assert_array_equal(actual, compare_features.matrices(rows)[predict.MODEL])

    def test_invalid_inputs_are_rejected(self):
        bad = [("AAAA", "False"), ("A" * 51, "False"), ("a" * 10, "False"),
               ("AAAAAXAAAA", "False"), ("A" * 10, "unknown"), ("A" * 10, True)]
        for sequence, state in bad:
            with self.subTest(sequence=sequence, state=state), self.assertRaises(ValueError):
                predict.prediction_features([{"sequence": sequence, "c_amidated": state}])
        with self.assertRaises(ValueError):
            predict.prediction_features([])

    def test_loading_and_scoring_never_fit_and_preserve_input(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            model, _ = fixture(root)
            rows = [{"sequence": "AKRDEAAAAA", "c_amidated": "True"}]
            before = {p.name: p.read_bytes() for p in root.iterdir()}
            with patch.object(Pipeline, "fit", side_effect=AssertionError("Inference must not fit")):
                loaded, digest = predict.load_primary(root)
                result = predict.score_rows(rows, loaded, digest)[0]
            self.assertAlmostEqual(result["model_probability_active"], model.predict_proba([[1, 1]])[0, 1])
            self.assertEqual(result["predicted_label"], "active" if model.predict([[1, 1]])[0] else "inactive")
            self.assertEqual(before, {p.name: p.read_bytes() for p in root.iterdir()})

    def test_tampered_model_rejected_before_deserialization(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            fixture(root)
            with (root / "primary_model.joblib").open("ab") as handle:
                handle.write(b"tampered")
            with patch.object(joblib, "load", side_effect=AssertionError("Must not load")):
                with self.assertRaisesRegex(ValueError, "fingerprint"):
                    predict.load_primary(root)

    def test_wrong_environment_and_policy_rejected_before_loading(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, manifest = fixture(root)
            path = root / "evaluation_manifest.json"
            for kind in ("environment", "policy"):
                modified = json.loads(json.dumps(manifest))
                if kind == "environment":
                    modified["spec"]["environment"]["scikit-learn"] = "wrong"
                else:
                    modified["spec"]["policy"]["primary_C"] = 10
                path.write_text(json.dumps(modified))
                with patch.object(joblib, "load", side_effect=AssertionError("Must not load")):
                    with self.assertRaises(ValueError):
                        predict.load_primary(root)

    def test_csv_cli_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            fixture(root)
            source, output = root / "input.csv", root / "predictions.csv"
            source.write_text("sequence,c_amidated\nAKRDEAAAAA,False\nAKRDEAAAAA,True\n")
            args = ["--input-csv", str(source), "--pretest-dir", str(root), "--output", str(output)]
            with contextlib.redirect_stdout(io.StringIO()):
                predict.main(args)
            before = output.read_bytes()
            self.assertIn(b"model_probability_active", before)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                predict.main(args)
            self.assertEqual(output.read_bytes(), before)
            source.write_text("sequence,sequence,c_amidated\nA,A,False\n")
            with self.assertRaises(ValueError):
                predict.input_rows(source)

    def test_single_cli_requires_chemistry_and_uses_inactive_tie(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            predict.main(["--sequence", "AKRDEAAAAA"])
        class TieModel:
            def predict_proba(self, x):
                return np.array([[0.5, 0.5]])
            def predict(self, x):
                return np.array([0])
        result = predict.score_rows([{"sequence": "AKRDEAAAAA", "c_amidated": "False"}], TieModel(), "test")[0]
        self.assertEqual(result["predicted_label"], "inactive")


if __name__ == "__main__":
    unittest.main(verbosity=2)
