"""Run: python modeling/test_compare_features.py -v

Uses the existing baseline test fixture; keep test_train_baseline.py nearby.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import modeling.compare_features as compare
import train_baseline as base
from test_train_baseline import fixture, write_csv


def grouped_fixture(root):
    rows = fixture(root)
    module, _ = compare.grouping_module()
    mapping, _, _ = module.sequence_groups([r["sequence"] for r in rows], compare.Fraction(4, 5))
    for name in ("train.csv", "sequence_groups.csv", "split_assignments.csv"):
        table = base.read_csv(root / name)
        for row in table:
            row["group_id"] = mapping[row["sequence"]]
        write_csv(root / name, table)
    manifest = json.loads((root / "split_manifest.json").read_text())
    for name in manifest["output_sha256"]:
        manifest["output_sha256"][name] = base.sha_file(root / name)
    (root / "split_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return base.read_csv(root / "train.csv")


class ComparisonTests(unittest.TestCase):
    def test_charge_correction_and_feature_definitions(self):
        rows = [{"sequence": "AKRDEHFWYGP", "c_amidated": amide} for amide in ("False", "True")]
        matrices = compare.matrices(rows)
        np.testing.assert_array_equal(matrices["charge_uncorrected"][:, 0], [0, 0])
        np.testing.assert_array_equal(matrices["charge_corrected"][:, 0], [0, 1])
        np.testing.assert_array_equal(matrices["charge_corrected_amidation"], [[0, 0], [1, 1]])
        np.testing.assert_allclose(matrices["composition"][:, :20].sum(axis=1), 1)
        for name, values in matrices.items():
            self.assertEqual(values.shape[1], len(compare.FEATURE_NAMES[name]))
        changed = [{**r, "label": "1", "note": "clinical", "id": 999} for r in rows]
        for name, values in matrices.items():
            np.testing.assert_array_equal(values, compare.matrices(changed)[name])

    def test_bootstrap_samples_complete_groups_with_same_draws(self):
        groups = np.asarray(["a", "a", "a", "b", "c", "c"])
        weights = compare.bootstrap_weights(groups, 30)
        np.testing.assert_array_equal(weights[:, 0], weights[:, 1])
        np.testing.assert_array_equal(weights[:, 1], weights[:, 2])
        np.testing.assert_array_equal(weights[:, 4], weights[:, 5])
        np.testing.assert_array_equal(weights[:, [0, 3, 4]].sum(axis=1), 3)
        np.testing.assert_array_equal(weights, compare.bootstrap_weights(groups, 30))

    def test_weighted_metrics_equal_explicit_repeated_rows_with_ties(self):
        y = np.asarray([0, 1, 0, 1, 0, 1])
        probability = np.asarray([0.2, 0.8, 0.5, 0.5, 0.8, 0.8])
        predicted = (probability > 0.5).astype(int)
        weights = np.asarray([3, 1, 2, 1, 1, 2])
        indices = np.repeat(np.arange(len(y)), weights)
        explicit = base.metrics(y[indices], predicted[indices], probability[indices])
        np.testing.assert_allclose(compare.weighted_metrics(y, predicted, probability, weights),
                                   [explicit[m] for m in compare.INTERVAL_METRICS])

    def test_identical_models_have_zero_paired_differences(self):
        y = np.tile([0, 1], 20)
        groups = np.repeat(np.arange(20), 2)
        probabilities = np.where(y == 1, 0.7, 0.3)
        predictions = [{"predicted_label": int(label), "probability_active": float(prob)}
                       for label, prob in zip(y, probabilities)]
        intervals, deltas, boot = compare.paired_bootstrap(y, groups, {"one": predictions, "same": predictions}, draws=30)
        self.assertEqual(len(intervals), 6)
        self.assertEqual(boot["valid_shared_draws"], 30)
        for delta in deltas:
            self.assertEqual((delta["pooled_oof_delta"], delta["low"], delta["high"]), (0, 0, 0))

    def test_saved_fold_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            rows = grouped_fixture(root)
            baseline_dir, _ = base.train_baseline(root, Path(directory) / "baseline")
            with (baseline_dir / "cv_folds.csv").open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "folds changed"):
                compare.load_saved_folds(baseline_dir, root, rows)

    def test_cluster_size_distribution_is_explicit(self):
        result = compare.cluster_description({"a": 1, "b": 1, "c": 3, "d": 5})
        self.assertEqual(result["n_sequences"], 10)
        self.assertEqual(result["n_groups"], 4)
        self.assertEqual(result["largest_fraction"], 0.5)
        self.assertEqual(result["singleton_groups"], 2)
        self.assertEqual(result["size_histogram"], {"1": 2, "3": 1, "5": 1})

    def test_complete_training_only_comparison_reuses_folds_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            grouped_fixture(root)
            baseline_dir, _ = base.train_baseline(root, Path(directory) / "baseline")
            before = {p: p.read_bytes() for folder in (root, baseline_dir) for p in folder.iterdir() if p.is_file()}
            original_open = Path.open

            def guard(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Heldout CSV was opened")
                return original_open(path, *args, **kwargs)

            with patch.object(Path, "open", guard):
                output, report = compare.run_comparison(root, baseline_dir, Path(directory) / "comparisons", draws=20)
            self.assertEqual(list(report["models"]), list(compare.MODEL_ORDER))
            self.assertEqual(report["heldout_csvs_opened"], [])
            self.assertFalse(report["final_test_evaluated"])
            self.assertEqual(len(report["grouping_sensitivity"]), 3)
            self.assertEqual((output / "cv_folds.csv").read_bytes(), (baseline_dir / "cv_folds.csv").read_bytes())
            self.assertEqual(len(base.read_csv(output / "feature_ladder.csv")), 6)
            self.assertEqual(len(base.read_csv(output / "paired_differences.csv")), 15)
            for p, content in before.items():
                self.assertEqual(p.read_bytes(), content)
            repeated, repeated_report = compare.run_comparison(root, baseline_dir, Path(directory) / "comparisons", draws=20)
            self.assertEqual(repeated, output)
            self.assertEqual(repeated_report, report)
            (output / "feature_ladder.csv").write_text("changed")
            with self.assertRaisesRegex(ValueError, "output changed"):
                compare.run_comparison(root, baseline_dir, Path(directory) / "comparisons", draws=20)


if __name__ == "__main__":
    unittest.main()
