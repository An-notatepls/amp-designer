"""Run: python modeling/test_prepare_evaluation.py -v

Keep test_train_baseline.py and all modeling scripts in the same folder.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import joblib
import numpy as np
from sklearn.dummy import DummyClassifier

import prepare_evaluation as prep
import train_baseline as base
import compare_features as compare
from test_train_baseline import fixture, write_csv


def comparison_fixture(parent):
    root = parent / "split"
    fixture(root)
    rows = base.read_csv(root / "train.csv")
    index = base.read_csv(root / "sequence_groups.csv")
    assignments = base.read_csv(root / "split_assignments.csv")
    for i, sequence in enumerate(("ACDEFGHIKLMN", "WWWWWWWWWWWW", "MNPQRSTVWYAC", "DEFGHIKLMNPQ")):
        assignments.append({"sequence": sequence, "c_amidated": "False", "group_id": f"test_family_{i}",
                            "partition": "test", "peptide_ids_json": json.dumps([1001 + i])})
        index.append({"sequence": sequence, "group_id": f"test_family_{i}", "group_partition": "test"})
    write_csv(root / "sequence_groups.csv", index)
    write_csv(root / "split_assignments.csv", assignments)
    manifest = json.loads((root / "split_manifest.json").read_text())
    manifest["counts"]["test"] = {"examples": 4, "unique_sequences": 4, "groups": 4, "active": 2, "inactive": 2}
    for name in manifest["output_sha256"]:
        manifest["output_sha256"][name] = base.sha_file(root / name)
    (root / "split_manifest.json").write_text(json.dumps(manifest))
    directory = parent / "comparison"
    directory.mkdir()
    y, groups = compare.labels_groups(rows)
    folds = base.make_folds(y, groups)
    xs = compare.matrices(rows)
    fold_rows = []
    for fold, (_, valid) in enumerate(folds, 1):
        for i in valid:
            fold_rows.append({"row_index": int(i), "sequence": rows[i]["sequence"], "c_amidated": rows[i]["c_amidated"],
                              "group_id": rows[i]["group_id"], "label": rows[i]["label"], "validation_fold": fold})
    write_csv(directory / "cv_folds.csv", fold_rows)
    all_predictions, models, summaries = [], {}, {}
    for model, c in prep.EXPECTED_C.items():
        estimator = DummyClassifier(strategy="most_frequent") if c is None else base.logistic_pipeline(c)
        _, predictions = base.evaluate_cv(estimator, xs[model], y, folds, model, c)
        predictions.sort(key=lambda r: r["row_index"])
        for item in predictions:
            item.update({key: rows[item["row_index"]][key] for key in ("sequence", "c_amidated", "group_id")})
        all_predictions.extend(predictions)
        estimator.fit(xs[model], y)
        models[model] = estimator
        summaries[model] = {"C": c, "pooled_oof": base.metrics(y,
            np.array([r["predicted_label"] for r in predictions]),
            np.array([r["probability_active"] for r in predictions]))}
    write_csv(directory / "selected_oof_predictions.csv", all_predictions)
    spec = {"train_sha256": base.sha_file(root / "train.csv"), "split_manifest_sha256": base.sha_file(root / "split_manifest.json")}
    joblib.dump({"models": models, "feature_names": compare.FEATURE_NAMES, "spec": spec}, directory / "models.joblib")
    report = {"spec": spec, "models": summaries,
              "output_sha256": {p.name: base.sha_file(p) for p in directory.iterdir()}}
    (directory / "comparison_report.json").write_text(json.dumps(report))
    return root, directory, rows


class PretestTests(unittest.TestCase):
    def test_truncation_is_detected_despite_low_global_similarity(self):
        fragment = "ACDEFGHIKLMN"
        longer = "RRRRRR" + fragment + "WWWWWWW"
        module, _ = compare.grouping_module()
        distance = module.edit_distance(fragment, longer)
        self.assertEqual(len(longer), 25)
        self.assertAlmostEqual(1 - distance / len(longer), 0.48)
        raw, qualified, complete = prep.local_pair(fragment, longer, prep.local_aligner())
        self.assertTrue(complete)
        self.assertEqual(qualified["identity"], 1)
        self.assertEqual(qualified["aligned_pairs"], 12)
        self.assertEqual(qualified["shorter_coverage"], 1)
        self.assertAlmostEqual(qualified["second_coverage"], 12 / 25)
        self.assertEqual(prep.novelty_flags(distance, len(longer), qualified, complete), (True, False, False))

    def test_tiny_perfect_local_match_fails_coverage_qualification(self):
        raw, qualified, complete = prep.local_pair("AAAAAAAAAAAA", "CCCCACCCCCCC", prep.local_aligner())
        self.assertEqual(raw["identity"], 1)
        self.assertEqual(raw["aligned_pairs"], 1)
        self.assertIsNone(qualified)
        self.assertTrue(complete)

    def test_internal_gaps_count_in_the_identity_denominator(self):
        first = "ACDEFGHIKLMNPQRSTVWY"
        second = first[:12] + "AAA" + first[12:]
        _, qualified, complete = prep.local_pair(first, second, prep.local_aligner())
        self.assertTrue(complete)
        self.assertEqual(qualified["matches"], 20)
        self.assertEqual(qualified["aligned_pairs"], 20)
        self.assertEqual(qualified["columns"], 23)
        self.assertAlmostEqual(qualified["identity"], 20 / 23)
        self.assertEqual(qualified["first_coverage"], 1)
        self.assertAlmostEqual(qualified["second_coverage"], 20 / 23)

    def test_exact_sixty_percent_and_incomplete_scan_are_not_strict(self):
        exact = {"matches": 6, "columns": 10}
        self.assertEqual(prep.novelty_flags(2, 5, exact, True), (False, False, False))
        self.assertEqual(prep.novelty_flags(3, 5, None, False), (True, False, False))
        self.assertEqual(prep.novelty_flags(3, 5, {"matches": 5, "columns": 10}, True), (True, True, True))
        _, _, complete = prep.local_pair("AAAAAAAAAAAA", "CCCCACCCCCCC", prep.local_aligner(), tie_limit=1)
        self.assertFalse(complete)
        row = {"sequence": "ACDEFGHIKLMN", "c_amidated": "False", "group_id": "train_family"}
        with self.assertRaisesRegex(ValueError, "violating the fixed 80%"):
            prep.neighbor_audit([row], [{**row, "group_id": "wrong_test_family", "partition": "test"}])

    def test_precision_ties_are_order_independent_and_weighted_correctly(self):
        y = np.array([1, 0, 1, 0])
        p = np.array([0.9, 0.5, 0.5, 0.5])
        result = prep.precision_at_budgets(y, p, (1, 2, 3, 5))
        self.assertEqual(result[1], 1)
        self.assertAlmostEqual(result[2], (1 + 1 / 3) / 2)
        self.assertIsNone(result[5])
        order = [3, 2, 1, 0]
        self.assertEqual(result, prep.precision_at_budgets(y[order], p[order], (1, 2, 3, 5)))
        weights = np.array([2, 1, 3, 1])
        duplicated = np.repeat(np.arange(4), weights)
        weighted = prep.precision_at_budgets(y, p, (2, 3, 5), weights)
        explicit = prep.precision_at_budgets(y[duplicated], p[duplicated], (2, 3, 5))
        self.assertEqual(weighted, explicit)
        self.assertEqual(prep.precision_at_budgets(y, np.ones(4), (1, 2)), {1: 0.5, 2: 0.5})

    def test_family_representative_uses_scores_and_identity_not_labels(self):
        rows = [{"sequence": "ACDEFGHIKLMN", "c_amidated": "False", "group_id": "one", "label": "1"},
                {"sequence": "CDEFGHIKLMNP", "c_amidated": "True", "group_id": "one", "label": "0"},
                {"sequence": "MNPQRSTVWYAC", "c_amidated": "False", "group_id": "two", "label": "0"}]
        indices = prep.family_representatives(rows, [0.5, 0.5, 0.6])
        np.testing.assert_array_equal(indices, [0, 2])
        changed = [{**r, "label": "0" if r["label"] == "1" else "1"} for r in rows]
        np.testing.assert_array_equal(indices, prep.family_representatives(changed, [0.5, 0.5, 0.6]))

    def test_histogram_preserves_missing_and_boundary_values(self):
        result = prep.identity_histogram([None, 0.0, 0.2, 0.6, 0.8, 1.0])
        self.assertEqual(sum(result.values()), 6)
        self.assertEqual(result["no_qualifying_alignment"], 1)
        self.assertEqual(result["[.6,.8)"], 1)
        self.assertEqual(result["[.8,1]"], 2)

    def test_complete_preparation_keeps_holdouts_sealed_and_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root, comparison_dir, rows = comparison_fixture(parent)
            before = {p: p.read_bytes() for folder in (root, comparison_dir) for p in folder.iterdir() if p.is_file()}
            original_open = Path.open

            def guard(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Held-out outcome CSV was opened")
                return original_open(path, *args, **kwargs)

            with patch.object(Path, "open", guard):
                output, report = prep.prepare(root, comparison_dir, parent / "pretest", draws=20)
            self.assertEqual(report["heldout_csvs_opened"], [])
            self.assertFalse(report["final_test_evaluated"])
            self.assertEqual(report["neighbor_summary"]["test_identities"], 4)
            membership = base.read_csv(output / "test_subset_membership.csv")
            self.assertEqual(len(membership), 4)
            self.assertTrue(all("label" not in row for row in membership))
            self.assertEqual(len(base.read_csv(output / "oof_precision_at_k.csv")), 36)
            self.assertEqual(len(base.read_csv(output / "oof_precision_at_k_by_fold.csv")), 180)
            model = joblib.load(output / "primary_model.joblib")
            self.assertEqual(model["pipeline"].named_steps["classifier"].C, 0.1)
            self.assertEqual(model["feature_names"], ["charge_proxy", "c_amidated"])
            self.assertEqual(report["primary_model_sha256"], base.sha_file(output / "primary_model.joblib"))
            for path, data in before.items():
                self.assertEqual(path.read_bytes(), data)
            repeated, repeated_report = prep.prepare(root, comparison_dir, parent / "pretest", draws=20)
            self.assertEqual((repeated, repeated_report), (output, report))
            (output / "evaluation_protocol.md").write_text("changed")
            with self.assertRaisesRegex(ValueError, "pre-test output changed"):
                prep.prepare(root, comparison_dir, parent / "pretest", draws=20)
            with (comparison_dir / "selected_oof_predictions.csv").open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "Comparison input changed"):
                prep.load_comparison(root, comparison_dir, rows)


if __name__ == "__main__":
    unittest.main()
