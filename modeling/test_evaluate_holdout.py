"""Synthetic-only tests. Run: python modeling/test_evaluate_holdout.py -v"""

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import joblib
import numpy as np
from sklearn.pipeline import Pipeline

import train_baseline as base
import prepare_evaluation as prep
import evaluate_holdout as evaluate
from test_prepare_evaluation import comparison_fixture, write_csv


def frozen_fixture(parent):
    root, comparison, training = comparison_fixture(parent)
    assignments = base.read_csv(root / "split_assignments.csv")
    index = base.read_csv(root / "sequence_groups.csv")
    test_rows = []
    for i, r in enumerate(a for a in assignments if a["partition"] == "test"):
        label = str(i % 2)
        test_rows.append({**training[0], **r, "label": label, "status": "active" if label == "1" else "inactive"})
    clinical = []
    for i, sequence in enumerate(("QWERTYIPASDFGHKL", "VILMFWYACDEGHKNP")):
        label = str(i)
        row = {**training[0], "sequence": sequence, "c_amidated": "False", "group_id": f"clinical_family_{i}",
               "partition": "clinical_exploratory", "peptide_ids_json": json.dumps([2001 + i]),
               "label": label, "status": "active" if label == "1" else "inactive"}
        clinical.append(row)
        assignments.append({k: row[k] for k in ("sequence", "c_amidated", "group_id", "partition", "peptide_ids_json")})
        index.append({"sequence": sequence, "group_id": row["group_id"], "group_partition": "clinical_family_holdout"})
    write_csv(root / "test.csv", test_rows)
    write_csv(root / "clinical_exploratory.csv", clinical)
    write_csv(root / "sequence_groups.csv", index)
    write_csv(root / "split_assignments.csv", assignments)
    split = json.loads((root / "split_manifest.json").read_text())
    split["counts"]["clinical_exploratory"] = {"examples": 2, "unique_sequences": 2, "groups": 2, "active": 1, "inactive": 1}
    split["counts"]["clinical_family_buffer"] = {"examples": 0, "unique_sequences": 0, "groups": 0, "active": 0, "inactive": 0}
    for name in (*split["output_sha256"], "test.csv", "clinical_exploratory.csv"):
        split["output_sha256"][name] = base.sha_file(root / name)
    (root / "split_manifest.json").write_text(json.dumps(split))
    comparison_report = json.loads((comparison / "comparison_report.json").read_text())
    comparison_report["spec"]["split_manifest_sha256"] = base.sha_file(root / "split_manifest.json")
    models = joblib.load(comparison / "models.joblib")
    models["spec"] = comparison_report["spec"]
    joblib.dump(models, comparison / "models.joblib")
    comparison_report["output_sha256"]["models.joblib"] = base.sha_file(comparison / "models.joblib")
    (comparison / "comparison_report.json").write_text(json.dumps(comparison_report))
    pretest, _ = prep.prepare(root, comparison, parent / "pretest", draws=20)
    return root, pretest, comparison


class HoldoutTests(unittest.TestCase):
    def test_vectorized_metrics_match_independent_sklearn_calls_with_ties(self):
        rng = np.random.default_rng(19)
        y = np.tile([0, 1], 8)
        p = rng.choice([0.1, 0.3, 0.5, 0.8], size=len(y))
        predicted = (p > 0.5).astype(int)
        weights = rng.integers(1, 4, size=(20, len(y)))
        actual = evaluate.weighted_metrics(y, predicted, p, weights)
        for i, w in enumerate(weights):
            repeated = np.repeat(np.arange(len(y)), w)
            expected = base.metrics(y[repeated], predicted[repeated], p[repeated])
            expected["precision"] = expected["tp"] / (expected["tp"] + expected["fp"]) if expected["tp"] + expected["fp"] else 0
            for key, value in expected.items():
                self.assertAlmostEqual(value, actual[key][i], places=12, msg=key)

    def test_empty_and_single_class_panels_remain_visible(self):
        rows = [{"sequence": "ACDEFGHIKLMN", "c_amidated": "False", "label": "1", "group_id": "one"}]
        probs = {m: np.array([0.7]) for m in prep.EXPECTED_C}
        preds = {m: np.array([1]) for m in prep.EXPECTED_C}
        metrics, intervals, ranking, support = evaluate.panel_results("single", rows, probs, preds, 20, 42)
        self.assertTrue(all(r["roc_auc"] is None and r["balanced_accuracy"] is None for r in metrics))
        self.assertTrue(all(r["low"] is None for r in intervals))
        self.assertEqual(support["n"], 1)
        empty, _, _, support = evaluate.panel_results("empty", [], probs, preds, 20, 42)
        self.assertEqual(len(empty), 6)
        self.assertEqual(support["n"], 0)

    def test_auc_support_rule_suppresses_only_the_unsupported_interval(self):
        rows = [{"sequence": "ACDEFGHIKL" + "ACDEFG"[i] + "N", "c_amidated": "False",
                 "label": str(i % 2), "group_id": f"family_{i}"} for i in range(6)]
        p = np.array([0.1, 0.6, 0.2, 0.7, 0.4, 0.8])
        metrics, intervals, _, support = evaluate.panel_results("tiny", rows,
            {m: p for m in prep.EXPECTED_C}, {m: (p > 0.5).astype(int) for m in prep.EXPECTED_C}, 100, 42)
        auc = [r for r in intervals if r["metric"] == "roc_auc"]
        self.assertTrue(all(r["estimate"] == 1 and r["low"] is None for r in auc))
        self.assertTrue(all("support" in r["interval_suppressed_reason"] for r in auc))

    def test_git_gate_checks_committed_bytes_including_evaluator(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / "pretest"
            folder.mkdir()
            for name in ("evaluation_protocol.md", "evaluation_manifest.json", "test_subset_membership.csv"):
                (folder / name).write_text("synthetic committed metadata")
            script = root / "evaluate_holdout.py"
            script.write_text("synthetic evaluator source")
            def git(*args):
                return subprocess.run(["git", "-C", str(root), *args], capture_output=True, check=True)
            git("init", "-q")
            git("add", ".")
            git("-c", "user.name=Synthetic Test", "-c", "user.email=test@example.invalid",
                "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", "commit", "-qm", "Synthetic freeze")
            with patch.object(evaluate, "__file__", str(script)):
                evidence = evaluate.git_evidence(folder)
                self.assertEqual(len(evidence["verified_committed_files"]), 4)
                (folder / "evaluation_protocol.md").write_text("changed")
                with self.assertRaisesRegex(ValueError, "Commit this file unchanged"):
                    evaluate.git_evidence(folder)

    def test_complete_evaluation_never_fits_and_reuses_saved_results(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, pretest, comparison = frozen_fixture(parent)
            gate_passed = False
            original_read = evaluate.read_frozen_outcomes
            opened = []
            def gate(path):
                nonlocal gate_passed
                gate_passed = True
                return {"commit": "synthetic", "verified_committed_files": []}
            def reader(path, *args):
                self.assertTrue(gate_passed)
                self.assertTrue((pretest / "evaluation_started.json").is_file())
                opened.append(path.name)
                return original_read(path, *args)
            with patch.object(evaluate, "git_evidence", gate), patch.object(evaluate, "read_frozen_outcomes", reader), patch.object(Pipeline, "fit", side_effect=AssertionError("Evaluation tried to fit")):
                output, report = evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
            self.assertEqual(opened, ["test.csv", "clinical_exploratory.csv"])
            self.assertFalse(report["model_fitting_performed"])
            self.assertFalse(report["model_selection_performed"])
            self.assertEqual(set(report["panels"]), {"overall", "global60_strict", "local60_strict", "combined60_strict", "clinical_exploratory"})
            self.assertEqual(report["primary_result"]["model"], prep.PRIMARY)
            self.assertEqual(len(base.read_csv(output / "metrics.csv")), 30)
            with patch.object(evaluate, "git_evidence", return_value={"commit": "synthetic"}), patch.object(evaluate, "read_frozen_outcomes", side_effect=AssertionError("Cached evaluation reopened outcomes")):
                repeated, repeated_report = evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
            self.assertEqual((repeated, repeated_report), (output, report))
            (output / "metrics.csv").write_text("changed")
            with patch.object(evaluate, "git_evidence", return_value={"commit": "synthetic"}):
                with self.assertRaisesRegex(ValueError, "evaluation artifact changed"):
                    evaluate.evaluate(split, pretest, comparison, parent / "evaluation")

    def test_changed_model_fails_before_git_gate_or_outcome_access(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, pretest, comparison = frozen_fixture(parent)
            (pretest / "primary_model.joblib").write_bytes(b"changed")
            with patch.object(evaluate, "git_evidence") as gate, patch.object(evaluate, "read_frozen_outcomes") as reader:
                with self.assertRaisesRegex(ValueError, "Frozen pre-test artifact changed"):
                    evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
                gate.assert_not_called()
                reader.assert_not_called()

    def test_changed_outcome_file_is_rejected_without_computing_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, pretest, comparison = frozen_fixture(parent)
            with (split / "test.csv").open("a") as handle:
                handle.write("\n")
            with patch.object(evaluate, "git_evidence", return_value={"commit": "synthetic"}), patch.object(evaluate, "panel_results") as scorer:
                with self.assertRaisesRegex(ValueError, "Frozen outcome table changed"):
                    evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
                scorer.assert_not_called()

    def test_interrupted_run_requires_explicit_identical_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, pretest, comparison = frozen_fixture(parent)
            with patch.object(evaluate, "git_evidence", return_value={"commit": "synthetic"}):
                with patch.object(evaluate, "read_frozen_outcomes", side_effect=RuntimeError("synthetic interruption")):
                    with self.assertRaisesRegex(RuntimeError, "synthetic interruption"):
                        evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
                with self.assertRaisesRegex(ValueError, "evaluation attempt is recorded"):
                    evaluate.evaluate(split, pretest, comparison, parent / "evaluation")
                _, report = evaluate.evaluate(split, pretest, comparison, parent / "evaluation", resume=True)
                self.assertEqual(report["primary_panel"], "overall")


if __name__ == "__main__":
    unittest.main()
