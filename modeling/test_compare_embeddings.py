"""Offline synthetic tests; never download or import a protein encoder."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import compare_embeddings as experiment
from test_review_training_notes import fixture as review_fixture
from test_train_baseline import fixture


def vectors(rows):
    return {r["sequence"]: np.asarray([r["sequence"].count(a) for a in "KRDEG"], dtype=float)
            for r in rows}


def protocol():
    return {"version": "synthetic-embedding-comparison", "baseline_C": 0.1,
            "embedding_C_grid": [0.01, 0.1], "bootstrap": {"draws": 30},
            "cv": {"folds": 5, "seed": 42, "method": "synthetic grouped folds"}}


class EmbeddingTests(unittest.TestCase):
    def test_pooling_excludes_special_and_padding_tokens(self):
        hidden = np.array([[[100, 100], [2, 4], [4, 8], [900, 900]]], dtype=float)
        actual = experiment.mean_residues(hidden, [[False, True, True, False]])
        np.testing.assert_array_equal(actual, [[3, 6]])
        with self.assertRaises(ValueError):
            experiment.mean_residues(hidden, [[False] * 4])

    def test_amidation_is_separate_from_shared_sequence_embedding(self):
        rows = [{"sequence": "AKRDEAAAAA", "c_amidated": state} for state in ("False", "True")]
        charge, protein = experiment.feature_matrices(rows, vectors(rows))
        np.testing.assert_array_equal(protein[0, :-1], protein[1, :-1])
        np.testing.assert_array_equal(protein[:, -1], [0, 1])
        np.testing.assert_array_equal(charge[:, 0], [0, 1])

    def test_all_arms_share_group_folds_and_original_arm_fits_fewer_rows(self):
        with tempfile.TemporaryDirectory() as temporary:
            rows = fixture(Path(temporary) / "split")
            for i, row in enumerate(rows):
                row["original_training_identity"] = str(i % 2 == 0)
            selected, predictions, scores, folds = experiment.compare(rows, vectors(rows), protocol())
            by_family = {}
            for row in folds:
                by_family.setdefault(row["group_id"], set()).add(row["validation_fold"])
            self.assertTrue(all(len(values) == 1 for values in by_family.values()))
            reference = [r["fold"] for r in predictions[experiment.ARMS[0]]]
            for arm in experiment.ARMS:
                self.assertEqual([r["fold"] for r in predictions[arm]], reference)
                self.assertEqual(len(predictions[arm]), 80)
            for row in scores:
                self.assertEqual(row["fit_rows"], 32 if row["model"] == experiment.ARMS[0] else 64)

    def test_scaling_uses_only_the_arm_training_subset(self):
        with tempfile.TemporaryDirectory() as temporary:
            rows = fixture(Path(temporary) / "split")
            y, groups = experiment.comparison.labels_groups(rows)
            folds = experiment.base.make_folds(y, groups)
            original = np.arange(len(rows)) % 2 == 0
            x = np.arange(len(rows) * 2, dtype=float).reshape(-1, 2)
            observed = []
            real_fit = experiment.base.StandardScaler.fit

            def fit(scaler, values, *args, **kwargs):
                observed.append(values.mean(axis=0))
                return real_fit(scaler, values, *args, **kwargs)

            with patch.object(experiment.base.StandardScaler, "fit", fit):
                experiment.fit_arm(x, y, folds, original, "synthetic", 0.1, original_only=True)
            for mean, (train, _) in zip(observed, folds):
                np.testing.assert_array_equal(mean, x[train[original[train]]].mean(axis=0))

    def test_identical_predictions_have_zero_paired_intervals(self):
        with tempfile.TemporaryDirectory() as temporary:
            rows = fixture(Path(temporary) / "split")
            predicted = [{"predicted_label": int(r["label"]), "probability_active": 0.8 if r["label"] == "1" else 0.2}
                         for r in rows]
            _, deltas, _ = experiment.intervals(rows, {arm: predicted for arm in experiment.ARMS}, draws=30)
            self.assertTrue(all(r["estimate"] == r["low"] == r["high"] == 0 for r in deltas))

    def test_full_offline_run_preserves_inputs_guards_holdouts_and_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split, original, overlay, _ = review_fixture(root)
            output, _ = experiment.review.build(split, original, overlay, root / "reviews")
            p = {**protocol(), "dataset": str(output / "training_candidates.csv"),
                 "dataset_sha256": experiment.base.sha_file(output / "training_candidates.csv"),
                 "review_report": str(output / "review_report.json"),
                 "review_report_sha256": experiment.base.sha_file(output / "review_report.json"),
                 "original_training_identities": 80, "candidate_training_identities": 81}
            path = root / "protocol.json"
            path.write_text(json.dumps(p))
            original_open = Path.open

            def guard(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Held-out classifier CSV opened")
                return original_open(path, *args, **kwargs)

            def fake_extract(rows, *args):
                return vectors(rows), {"fixture": True, "encoder_fine_tuned": False}

            before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
            with patch.object(experiment.review, "DEFAULT_SPLIT", split), patch.object(Path, "open", guard), \
                    patch.object(experiment, "extract_embeddings", fake_extract):
                destination, report = experiment.run(path, root / "cache", root / "results")
                self.assertFalse(report["encoder_fine_tuned"])
                self.assertFalse(report["v1_model_changed"])
                self.assertEqual(report["heldout_classifier_csvs_opened"], [])
                self.assertEqual(experiment.run(path, root / "cache", root / "results"), (destination, report))
                (destination / "cv_folds.csv").write_text("changed")
                with self.assertRaisesRegex(ValueError, "output changed"):
                    experiment.run(path, root / "cache", root / "results")
            for source, content in before.items():
                self.assertEqual(source.read_bytes(), content)
            (output / "training_candidates.csv").write_text("changed")
            with patch.object(experiment.review, "DEFAULT_SPLIT", split):
                with self.assertRaisesRegex(ValueError, "fingerprint"):
                    experiment.load_candidates(p)

    def test_changed_cache_fails_before_loading_a_model(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = [{"sequence": "AAAAAAAAAA"}]
            p = {"model_id": "synthetic", "model_revision": "fixed", "representation": "residue mean"}
            spec = {"sequences": ["AAAAAAAAAA"], "model_id": "synthetic", "revision": "fixed", "pooling": "residue mean"}
            key = hashlib.sha256(experiment.base.canonical(spec).encode()).hexdigest()[:16]
            cache = root / key
            cache.mkdir()
            (cache / "embeddings.npz").write_bytes(b"changed")
            (cache / "embedding_manifest.json").write_text(json.dumps({"spec": spec, "embeddings_sha256": "original"}))
            with self.assertRaisesRegex(ValueError, "cache changed"):
                experiment.extract_embeddings(rows, p, root)


if __name__ == "__main__":
    unittest.main(verbosity=2)
