"""Offline tests: features, grouped validation, training-only access and provenance.

Run: python -m unittest -v test_train_baseline
"""

import csv
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import joblib
import numpy as np

import train_baseline as baseline


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture(root):
    root.mkdir()
    rng = random.Random(123)
    rows, mapping, assignments = [], [], []
    for group in range(40):
        sequence = "".join(rng.choice(sorted(baseline.STANDARD)) for _ in range(16))
        label = group % 2
        group_id = f"fixture_family_{group}"
        mapping.append({"sequence": sequence, "group_id": group_id, "group_partition": "train"})
        for amide in ("False", "True"):
            row = {"sequence": sequence, "c_amidated": amide, "label": str(label),
                   "status": "active" if label else "inactive", "model_eligible": "True",
                   "group_id": group_id, "partition": "train", "peptide_ids_json": json.dumps([len(rows) + 1])}
            rows.append(row)
            assignments.append({k: row[k] for k in ("sequence", "c_amidated", "group_id", "partition", "peptide_ids_json")})
    write_csv(root / "train.csv", rows)
    write_csv(root / "sequence_groups.csv", mapping)
    write_csv(root / "split_assignments.csv", assignments)
    checks = {"identity_partition_unique": True,
              "no_family_overlap_across_train_test_and_clinical_reservation": True,
              "all_qualifying_pairs_protected": True,
              "clinical_families_absent_from_training": True,
              "primary_train_and_test_have_both_classes": True}
    manifest = {"checks": checks, "counts": {"train": {"examples": 80, "unique_sequences": 40,
                 "groups": 40, "active": 40, "inactive": 40}},
                "output_sha256": {name: baseline.sha_file(root / name) for name in
                                   ("train.csv", "sequence_groups.csv", "split_assignments.csv")}}
    (root / "split_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    for name in ("test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"):
        (root / name).write_text("DELIBERATELY INVALID: MUST NEVER BE READ", encoding="utf-8")
    return rows


class BaselineTests(unittest.TestCase):
    def test_features_match_hand_counts_and_amidation(self):
        expected = [11, 2 / 11, 2 / 11, 1 / 11, 4 / 11, 3 / 11, 1 / 11, 1 / 11, 1, 1]
        np.testing.assert_allclose(baseline.feature_vector("AKRDEHFWYGP", "True"), expected)
        plain = baseline.feature_vector("AKRDEHFWYGP", "False")
        self.assertEqual(plain[-2:], [0, 0])
        for seq, amide in (("akrdehfwyGP", "False"), ("AXRDEHFWYGP", "False"), ("AAA", "False"), ("AKRDEHFWYGP", "unknown")):
            with self.assertRaises(ValueError):
                baseline.feature_vector(seq, amide)

    def test_features_cannot_use_label_note_or_id(self):
        row = {"sequence": "AKRDEHFWYGP", "c_amidated": "False", "label": "0", "note": "one", "id": "1"}
        changed = {**row, "label": "1", "note": "opposite", "id": "999"}
        np.testing.assert_array_equal(baseline.feature_matrix([row]), baseline.feature_matrix([changed]))

    def test_cv_families_stay_together_and_validation_covers_once(self):
        y = np.repeat(np.arange(40) % 2, 2)
        groups = np.repeat(np.arange(40), 2)
        folds = baseline.make_folds(y, groups)
        seen = []
        for train, valid in folds:
            self.assertFalse(set(groups[train]) & set(groups[valid]))
            self.assertEqual(set(y[valid]), {0, 1})
            seen.extend(valid.tolist())
        self.assertEqual(sorted(seen), list(range(80)))

    def test_modified_training_bytes_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            fixture(root)
            with (root / "train.csv").open("a", encoding="utf-8") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "changed"):
                baseline.load_training(root)

    def test_inconsistent_family_metadata_rejected_even_with_updated_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            fixture(root)
            rows = baseline.read_csv(root / "split_assignments.csv")
            rows[0]["partition"] = "test"
            write_csv(root / "split_assignments.csv", rows)
            manifest = json.loads((root / "split_manifest.json").read_text())
            manifest["output_sha256"]["split_assignments.csv"] = baseline.sha_file(root / "split_assignments.csv")
            (root / "split_manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "identity set"):
                baseline.load_training(root)

    def test_complete_run_never_opens_holdouts_and_scales_within_each_fold(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            rows = fixture(root)
            x = baseline.feature_matrix(rows)
            y = np.array([int(r["label"]) for r in rows])
            groups = np.array([r["group_id"] for r in rows])
            folds = baseline.make_folds(y, groups)
            observed_means = []
            original_open = Path.open
            original_fit = baseline.StandardScaler.fit

            def guarded_open(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Held-out data was opened")
                return original_open(path, *args, **kwargs)

            def observed_fit(scaler, values, *args, **kwargs):
                fitted = original_fit(scaler, values, *args, **kwargs)
                observed_means.append(fitted.mean_.copy())
                return fitted

            with patch.object(Path, "open", guarded_open), patch.object(baseline.StandardScaler, "fit", observed_fit):
                output, report = baseline.train_baseline(root, Path(directory) / "results")
            expected_means = [x[train].mean(axis=0) for _ in baseline.C_VALUES for train, _ in folds] + [x.mean(axis=0)]
            self.assertEqual(len(observed_means), len(expected_means))
            np.testing.assert_allclose(observed_means, expected_means)
            self.assertFalse(report["final_test_evaluated"])
            self.assertEqual(report["heldout_csvs_opened"], [])
            self.assertAlmostEqual(report["dummy_cv"]["balanced_accuracy"]["mean"], 0.5)
            saved = joblib.load(output / "model.joblib")
            np.testing.assert_allclose(saved["pipeline"].named_steps["scale"].mean_, x.mean(axis=0))
            self.assertEqual(saved["pipeline"].predict_proba(x).shape, (80, 2))
            self.assertEqual(saved["feature_names"], list(baseline.FEATURES))

    def test_result_reuse_and_tampered_export_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "split"
            fixture(root)
            output, first = baseline.train_baseline(root, Path(directory) / "results")
            before = (output / "training_report.json").read_bytes()
            repeated_output, second = baseline.train_baseline(root, Path(directory) / "results")
            self.assertEqual(output, repeated_output)
            self.assertEqual(first, second)
            self.assertEqual(before, (output / "training_report.json").read_bytes())
            (output / "coefficients.csv").write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Existing output changed"):
                baseline.train_baseline(root, Path(directory) / "results")


if __name__ == "__main__":
    unittest.main()
