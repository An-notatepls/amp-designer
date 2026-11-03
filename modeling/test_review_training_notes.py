"""Synthetic checks for training-only review overlays and immutable v1 inputs."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import review_training_notes as review
from test_label_sensitivity import sensitivity_fixture
from test_train_baseline import write_csv


def fixture(root):
    split, original_reviews, source, _ = sensitivity_fixture(root)
    measurements = review.base.read_csv(source)
    for row in measurements:
        if row["peptide_id"] == "9000":
            row.update(note_raw_json=json.dumps("Clinical isolate GD004"),
                       note_review_status="pending", point_uM="8.0", lower_uM="8.0", upper_uM="8.0")
    write_csv(source, measurements)
    cohort = split.parent.parent
    audit_path = cohort / "broad_audit.csv"
    audit = review.base.read_csv(audit_path)
    for row in audit:
        if json.loads(row["peptide_ids_json"]) == [9000]:
            row.update(status="needs_note_review", model_eligible="False", label="")
    write_csv(audit_path, audit)
    manifest = json.loads((split / "split_manifest.json").read_text())
    manifest["spec"]["source_sha256"]["broad_audit.csv"] = review.fingerprint(audit_path)
    (split / "split_manifest.json").write_text(json.dumps(manifest))
    report_path = cohort / "provenance_report.json"
    report = json.loads(report_path.read_text())
    report["policy"]["source_measurements_sha256"] = review.fingerprint(source)
    report_path.write_text(json.dumps(report))
    _, roster, rows, labels, _ = review.load_context(split, original_reviews)
    queue = review.pending_queue(rows, roster, labels)
    assert len(queue) == 2
    overlay = root / "new_reviews.csv"
    for row in queue:
        row.update(decision="include", review_reason="Synthetic single-isolate annotation review")
    write_csv(overlay, queue)
    return split, original_reviews, overlay, source


class TrainingReviewTests(unittest.TestCase):
    def test_build_filters_before_parsing_never_opens_holdouts_and_preserves_v1(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split, original, overlay, _ = fixture(root)
            before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
            original_open = Path.open

            def guard(path, *args, **kwargs):
                if path.name in {"test.csv", "clinical_exploratory.csv", "clinical_family_buffer.csv"}:
                    raise AssertionError("Held-out classifier CSV opened")
                return original_open(path, *args, **kwargs)

            with patch.object(Path, "open", guard):
                output, report = review.build(split, original, overlay, root / "results")
            self.assertEqual(report["counts"]["candidate_training_identities"], 81)
            self.assertEqual(report["counts"]["newly_eligible_identities"], 1)
            self.assertEqual(report["counts"]["remaining_pending_measurements"], 0)
            self.assertFalse(report["models_fitted"])
            self.assertFalse(report["accuracy_improvement_estimated"])
            self.assertEqual(report["heldout_classifier_csvs_opened"], [])
            for path, contents in before.items():
                self.assertEqual(path.read_bytes(), contents)
            self.assertEqual(review.build(split, original, overlay, root / "results"), (output, report))
            (output / "training_candidates.csv").write_text("tampered")
            with self.assertRaisesRegex(ValueError, "output changed"):
                review.build(split, original, overlay, root / "results")

    def test_duplicate_nontraining_and_previously_reviewed_keys_are_rejected(self):
        for mode in ("duplicate", "nontraining", "previously_reviewed"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                split, original, overlay, _ = fixture(root)
                rows = review.base.read_csv(overlay)
                if mode == "duplicate":
                    rows.append(dict(rows[0]))
                elif mode == "previously_reviewed":
                    previous = review.base.read_csv(original)[0]
                    rows[0].update(peptide_id=previous["peptide_id"], measurement_id=previous["measurement_id"])
                else:
                    rows[0]["peptide_id"] = "999999"
                write_csv(overlay, rows)
                with self.assertRaisesRegex(ValueError, "Duplicate or non-pending training review"):
                    review.build(split, original, overlay, root / "results")

    def test_stale_note_and_changed_family_are_rejected(self):
        for field in ("note", "group_id"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                split, original, overlay, _ = fixture(root)
                rows = review.base.read_csv(overlay)
                rows[0][field] = "changed"
                write_csv(overlay, rows)
                with self.assertRaisesRegex(ValueError, "changed|differs"):
                    review.build(split, original, overlay, root / "results")

    def test_blank_decision_keeps_identity_pending(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split, original, overlay, _ = fixture(root)
            rows = review.base.read_csv(overlay)
            rows[0]["decision"] = ""
            write_csv(overlay, rows)
            _, report = review.build(split, original, overlay, root / "results")
            self.assertEqual(report["counts"]["candidate_training_identities"], 80)
            self.assertEqual(report["counts"]["identity_statuses"]["needs_note_review"], 1)

    def test_excluding_all_rows_is_not_an_inactive_label(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split, original, overlay, _ = fixture(root)
            rows = review.base.read_csv(overlay)
            for row in rows:
                row["decision"] = "exclude"
            write_csv(overlay, rows)
            _, report = review.build(split, original, overlay, root / "results")
            self.assertEqual(report["counts"]["candidate_training_identities"], 80)
            self.assertEqual(report["counts"]["identity_statuses"]["no_remaining_groups"], 1)

    def test_frozen_source_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split, original, overlay, source = fixture(root)
            source.write_text(source.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "Frozen source changed: measurements"):
                review.build(split, original, overlay, root / "results")


if __name__ == "__main__":
    unittest.main(verbosity=2)
