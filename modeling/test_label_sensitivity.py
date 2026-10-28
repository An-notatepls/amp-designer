"""Synthetic-only tests: python modeling/test_label_sensitivity.py -v"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import train_baseline as base
import check_label_sensitivity as sensitivity
from test_train_baseline import fixture, write_csv


def sensitivity_fixture(parent):
    cohort = parent / "processed" / "provenance" / "cohort"
    split = cohort / "splits" / "fixed"
    split.parent.mkdir(parents=True)
    fixture(split)
    training = base.read_csv(split / "train.csv")
    labels, provenance = sensitivity.label_modules()
    measurements = []
    for i, item in enumerate(training):
        pid = json.loads(item["peptide_ids_json"])[0]
        values = [(8.0 if item["label"] == "1" else 96.0, None)]
        if i == 0:
            values = [(48.0, None)]
        elif i == 1:
            values = [(32.0, ">")]
        elif i == 2:
            values = [(24.0, None)]
        elif i == 3:
            values = [(16.0, "range")]
        elif i == 5:
            values = [(48.0, None), (96.0, None)]
        elif i == 7:
            values = [(8.0, None), (24.0, None)]
        for j, (value, op) in enumerate(values):
            note = "laboratory strain" if i == 0 else None
            row = {"peptide_id": pid, "measurement_id": 10000 + 10 * pid + j,
                   "sequence": item["sequence"], "c_amidated": item["c_amidated"],
                   "strain": "Escherichia coli", "medium": "MHB",
                   "medium_raw_json": json.dumps({"name": "MHB"}), "note_raw_json": json.dumps(note),
                   "note_review_status": "pending" if note else "not_needed", "note_review_reason": "",
                   "classification_eligible": "True", "regression_eligible": "False",
                   "exclusion_reasons": "[]", "op": op or "=", "point_uM": value if op is None else None,
                   "lower_uM": value, "upper_uM": None if op == ">" else 32.0 if op == "range" else value,
                   "log10_mic_uM": None}
            row.update({f"{key}_raw_json": "null" for key in labels.ASSAY_KEYS})
            measurements.append(row)
    # An identity mixed at 32 can become newly eligible at another cutoff.
    latent = {**training[0], "sequence": "ACDEFGHIKLMNPQRS", "peptide_ids_json": "[9000]",
              "label": "", "status": "mixed", "model_eligible": "False", "group_id": "family_latent"}
    for j, value in enumerate((24.0, 48.0)):
        measurements.append({**measurements[0], "peptide_id": 9000, "measurement_id": 90000 + j,
                             "sequence": latent["sequence"], "note_raw_json": "null", "note_review_status": "not_needed",
                             "point_uM": value, "lower_uM": value, "upper_uM": value, "op": "="})
    group_rows = base.read_csv(split / "sequence_groups.csv")
    group_rows.append({"sequence": latent["sequence"], "group_id": latent["group_id"], "group_partition": "train"})
    write_csv(split / "sequence_groups.csv", group_rows)
    # This out-of-scope MIC is deliberately unparsable: filtering must happen first.
    measurements.append({**measurements[0], "peptide_id": 999999, "measurement_id": 999999,
                         "sequence": "WWWWWWWWWWWW", "point_uM": "DO NOT PARSE", "lower_uM": "DO NOT PARSE",
                         "upper_uM": "DO NOT PARSE", "classification_eligible": "DO NOT PARSE"})
    processed = cohort.parent.parent
    source = processed / "ecoli_measurements.csv"
    write_csv(source, measurements)
    review = parent / "reviews.csv"
    write_csv(review, [{"peptide_id": measurements[0]["peptide_id"], "measurement_id": measurements[0]["measurement_id"],
                       "note": "laboratory strain", "decision": "include", "review_reason": "Synthetic ordinary origin note"},
                      {"peptide_id": 999999, "measurement_id": 999999, "note": "outside training", "decision": "", "review_reason": "Not processed"}])
    build_report_path = processed / "build_report.json"
    build_report_path.write_text(json.dumps({"policy": {"script_sha256": base.sha_file(Path(labels.__file__))}}))
    report = {"policy": {"source_measurements_sha256": base.sha_file(source),
                          "source_report_sha256": base.sha_file(build_report_path),
                          "note_reviews_sha256": base.sha_file(review),
                          "script_sha256": base.sha_file(Path(provenance.__file__))}}
    (cohort / "provenance_report.json").write_text(json.dumps(report))
    split_report = json.loads((split / "split_manifest.json").read_text())
    source_hashes = {}
    for name in ("broad_candidates.csv", "broad_audit.csv", "clinical_only_candidates.csv"):
        if name == "broad_audit.csv":
            audit = [{k: v for k, v in r.items() if k not in {"group_id", "partition"}} for r in [*training, latent]]
            write_csv(cohort / name, audit)
        else:
            (cohort / name).write_text("synthetic metadata hash anchor\n")
        source_hashes[name] = base.sha_file(cohort / name)
    split_report["spec"] = {"source_sha256": source_hashes}
    split_report["output_sha256"]["sequence_groups.csv"] = base.sha_file(split / "sequence_groups.csv")
    (split / "split_manifest.json").write_text(json.dumps(split_report))
    return split, review, source, training


def selected_rows(source, review, training, directory):
    labels, provenance = sensitivity.label_modules()
    selected, reviews, *_ = sensitivity.select_training_rows(source, review, training, directory)
    rows = provenance.read_measurements(selected)
    labels.apply_note_reviews(rows, reviews)
    return rows, labels


class LabelSensitivityTests(unittest.TestCase):
    def test_complete_training_only_check_replays_32_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, review, source, training = sensitivity_fixture(parent)
            before = {p: p.read_bytes() for p in parent.rglob('*') if p.is_file()}
            original = Path.open
            def guard(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Held-out classifier CSV opened")
                return original(path, *args, **kwargs)
            with patch.object(Path, "open", guard):
                output, report = sensitivity.check(split, review, parent / "results")
            self.assertTrue(report["replayed_32_uM_labels_match"])
            self.assertFalse(report["models_fitted"])
            self.assertFalse(report["frozen_primary_model_changed"])
            self.assertEqual(report["heldout_classifier_csvs_opened"], [])
            self.assertEqual(report["training_source_peptide_ids"], 81)
            at32 = next(r for r in report["summary"] if r["threshold_uM"] == 32)
            self.assertEqual((at32["eligible_identities"], at32["active"], at32["inactive"]), (80, 40, 40))
            self.assertEqual(len(base.read_csv(output / "training_identity_labels.csv")), 81)
            for threshold in (16, 64):
                expanded = next(r for r in report["summary"] if r["threshold_uM"] == threshold and r["scope"] == "all_training_assigned_audit_identities")
                self.assertEqual(expanded["newly_eligible_from_32"], 1)
            for path, payload in before.items():
                self.assertEqual(path.read_bytes(), payload)
            again, repeated = sensitivity.check(split, review, parent / "results")
            self.assertEqual((again, repeated), (output, report))
            (output / "threshold_summary.csv").write_text("changed")
            with self.assertRaisesRegex(ValueError, "sensitivity output changed"):
                sensitivity.check(split, review, parent / "results")

    def test_ranges_and_censored_values_remain_unresolved_when_needed(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            _, review, source, training = sensitivity_fixture(parent)
            rows, labels = selected_rows(source, review, training, parent)
            details, _, _ = sensitivity.rebuild(training, rows, labels)
            by = {(r["sequence"], r["c_amidated"]): r for r in details}
            self.assertEqual(by[base.identity(training[1])]["status_64"], "ambiguous")
            self.assertEqual(by[base.identity(training[3])]["status_16"], "ambiguous")
            self.assertEqual(by[base.identity(training[2])]["label_16"], 0)
            self.assertEqual(by[base.identity(training[0])]["label_64"], 1)

    def test_conflicting_measurements_are_not_resolved_by_majority_vote(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            _, review, source, training = sensitivity_fixture(parent)
            rows, labels = selected_rows(source, review, training, parent)
            details, _, _ = sensitivity.rebuild(training, rows, labels)
            by = {(r["sequence"], r["c_amidated"]): r for r in details}
            self.assertEqual(by[base.identity(training[5])]["status_64"], "mixed")
            self.assertIsNone(by[base.identity(training[5])]["label_64"])
            self.assertEqual(by[base.identity(training[7])]["status_16"], "mixed")

    def test_changed_review_snapshot_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            split, review, _, _ = sensitivity_fixture(parent)
            with review.open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "Input differs"):
                sensitivity.check(split, review, parent / "results")

    def test_selected_identity_mismatch_is_rejected_before_bound_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            _, review, source, training = sensitivity_fixture(parent)
            rows = base.read_csv(source)
            rows[0]["sequence"] = "ACDEFGHIKLMN"
            write_csv(source, rows)
            with self.assertRaisesRegex(ValueError, "identity differs"):
                sensitivity.select_training_rows(source, review, training, parent)

    def test_duplicate_selected_measurement_is_rejected_by_existing_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            _, review, source, training = sensitivity_fixture(parent)
            rows = base.read_csv(source)
            rows.append(rows[0].copy())
            write_csv(source, rows)
            with self.assertRaisesRegex(ValueError, "Duplicate measurement key"):
                selected_rows(source, review, training, parent)

    def test_32_label_drift_stops_the_check(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            _, review, source, training = sensitivity_fixture(parent)
            rows, labels = selected_rows(source, review, training, parent)
            first_id = json.loads(training[0]["peptide_ids_json"])[0]
            for row in rows:
                if row["peptide_id"] == first_id:
                    row.update(point_uM=1.0, lower_uM=1.0, upper_uM=1.0)
            with self.assertRaisesRegex(ValueError, "replay differs from frozen training labels"):
                sensitivity.rebuild(training, rows, labels)


if __name__ == "__main__":
    unittest.main()
