"""Build a separate v2 training candidate table from frozen training families.

Reviews are explicit, note-matched decisions, not predictions or keyword-based
automatic approvals. Original source files, v1 labels and model stay unchanged.
Only training-selected measurement strings are parsed as outcomes.
"""

import argparse
import csv
import hashlib
import io
import json
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import check_label_sensitivity as sensitivity
import train_baseline as base

VERSION = "training-note-review-v2"
ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "data/processed/ecoli-full-v1-3cec454243584936/45c31cc2e742/provenance/83780b72da84"
DEFAULT_SPLIT = COHORT / "splits/grouped-v1-dbe608485cce9a7e"
QUEUE_HEADER = [
    "peptide_id", "measurement_id", "sequence", "c_amidated", "group_id",
    "note", "strain", "medium", "reference_raw", "decision", "review_reason",
]


def row_key(row):
    return str(row["peptide_id"]), str(row["measurement_id"])


def csv_text(rows, fields):
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue()


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_context(split_dir, original_reviews):
    """Verify the historical inputs and replay v1 before applying new reviews."""
    split_dir, original_reviews = Path(split_dir), Path(original_reviews)
    training, manifest = base.load_training(split_dir)
    cohort, processed = split_dir.parent.parent, split_dir.parent.parent.parent.parent
    policy = json.loads((cohort / "provenance_report.json").read_text())["policy"]
    build_report = processed / "build_report.json"
    report = json.loads(build_report.read_text())
    labels, provenance = sensitivity.label_modules()
    inputs = {
        "split_manifest": split_dir / "split_manifest.json",
        "original_train": split_dir / "train.csv",
        "sequence_groups": split_dir / "sequence_groups.csv",
        "split_assignments": split_dir / "split_assignments.csv",
        "measurements": processed / "ecoli_measurements.csv",
        "original_reviews": original_reviews,
        "build_report": build_report,
        "provenance_report": cohort / "provenance_report.json",
        "label_script": Path(labels.__file__),
        "provenance_script": Path(provenance.__file__),
        "training_helper": Path(base.__file__),
        "sensitivity_helper": Path(sensitivity.__file__),
        "data_utils": Path(labels.__file__).with_name("data_utils.py"),
    }
    expected = {
        "measurements": policy["source_measurements_sha256"],
        "original_reviews": policy["note_reviews_sha256"],
        "build_report": policy["source_report_sha256"],
        "label_script": report["policy"]["script_sha256"],
        "provenance_script": policy["script_sha256"],
    }
    for name, digest in expected.items():
        if fingerprint(inputs[name]) != digest:
            raise ValueError(f"Frozen source changed: {name}")
    for name, digest in manifest["spec"]["source_sha256"].items():
        if Path(name).name != name or fingerprint(cohort / name) != digest:
            raise ValueError(f"Frozen cohort source changed: {name}")
        inputs[f"cohort_{name}"] = cohort / name
    roster = sensitivity.training_roster(cohort, split_dir, training)
    with tempfile.TemporaryDirectory(prefix="amp-v2-training-") as temporary:
        selected, reviews, *_ = sensitivity.select_training_rows(
            inputs["measurements"], original_reviews, roster, Path(temporary))
        rows = provenance.read_measurements(selected)
        labels.apply_note_reviews(rows, reviews)
    # This existing check also compares source IDs and group-level counts.
    sensitivity.rebuild(roster, rows, labels, frozen_training=training)
    return training, roster, rows, labels, inputs


def pending_queue(rows, roster, labels):
    groups = {base.identity(r): r["group_id"] for r in roster}
    queue = []
    for row in rows:
        if not row["classification_eligible"] or row["note_review_status"] != "pending":
            continue
        amidated = str(row["c_amidated"])
        note = labels.name_of(json.loads(row["note_raw_json"]))
        if not note:
            raise ValueError("Pending row has no note")
        queue.append({
            "peptide_id": row["peptide_id"], "measurement_id": row["measurement_id"],
            "sequence": row["sequence"], "c_amidated": amidated,
            "group_id": groups[(row["sequence"], amidated)],
            "note": note, "strain": row["strain"], "medium": row["medium"],
            "reference_raw": row.get("reference_raw", ""),
            "decision": "", "review_reason": "",
        })
    return sorted(queue, key=lambda r: (int(r["peptide_id"]), int(r["measurement_id"])))


def apply_overlay(rows, queue, reviews, labels):
    expected = {row_key(row): row for row in queue}
    seen = set()
    for review in base.read_csv(reviews):
        if not {"peptide_id", "measurement_id", "note", "decision", "review_reason"} <= set(review):
            raise ValueError("New review file lacks required columns")
        key = row_key(review)
        if key in seen or key not in expected:
            raise ValueError(f"Duplicate or non-pending training review: {key}")
        seen.add(key)
        current = expected[key]
        if review["note"] != current["note"]:
            raise ValueError(f"Review note changed: {key}")
        for field in ("sequence", "c_amidated", "group_id"):
            if field in review and review[field] != str(current[field]):
                raise ValueError(f"Review {field} differs from training metadata: {key}")
        decision = review["decision"].strip().casefold()
        if decision and (decision not in {"include", "exclude"} or not review["review_reason"].strip()):
            raise ValueError(f"Review requires include/exclude and a reason: {key}")
    labels.apply_note_reviews(rows, reviews)


def build(split_dir, original_reviews, reviews, results_root):
    training, roster, rows, labels, inputs = load_context(split_dir, original_reviews)
    initial_queue = pending_queue(rows, roster, labels)
    apply_overlay(rows, initial_queue, reviews, labels)
    rebuilt = labels.classifier_table(labels.aggregate(rows, 32.0))
    by_identity = {(r["sequence"], str(r["c_amidated"])): r for r in rebuilt}
    if set(by_identity) - {base.identity(r) for r in roster}:
        raise ValueError("v2 rebuild escaped frozen training families")
    for row in training:
        updated = by_identity.get(base.identity(row))
        if not updated or not updated["model_eligible"] or updated["label"] != int(row["label"]):
            raise ValueError("New reviews changed an original ready training label")
    roster_map = {base.identity(r): r for r in roster}
    original_keys = {base.identity(r) for r in training}
    audit, ready = [], []
    for key, original in sorted(roster_map.items()):
        item = by_identity.get(key)
        if item is None:
            # Excluding every retained measurement is not an inactive label.
            item = {"sequence": key[0], "c_amidated": key[1], "peptide_ids_json": original["peptide_ids_json"],
                    "n_assay_groups": 0, "active_groups": 0, "inactive_groups": 0,
                    "mixed_labels": False, "needs_note_review": False, "has_unresolved_group": False,
                    "label": None, "status": "no_remaining_groups", "model_eligible": False}
        item = {**item, "group_id": original["group_id"], "partition": "train",
                "original_training_identity": key in original_keys,
                "original_status": original["status"]}
        audit.append(item)
        if item["model_eligible"]:
            ready.append(item)
    remaining = pending_queue(rows, roster, labels)
    source_hashes = {name: fingerprint(path) for name, path in inputs.items()}
    source_hashes.update(new_reviews=fingerprint(reviews), review_script=fingerprint(__file__))
    spec = {"version": VERSION, "threshold_uM": 32, "source_sha256": source_hashes,
            "scope": "frozen training-assigned families only; annotation review under unchanged broad labeling rules"}
    digest = hashlib.sha256(base.canonical(spec).encode()).hexdigest()[:16]
    output = Path(results_root) / f"{VERSION}-{digest}"
    if output.exists():
        report = json.loads((output / "review_report.json").read_text())
        if report["spec"] != spec:
            raise ValueError("Existing v2 review specification differs")
        for name, expected in report["output_sha256"].items():
            if fingerprint(output / name) != expected:
                raise ValueError(f"Existing v2 output changed: {name}")
        return output, report
    transitions = Counter((r["original_status"], r["status"]) for r in audit)
    report = {
        "spec": spec, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "counts": {"training_family_audit_identities": len(roster), "original_training_identities": len(training),
                   "initial_pending_measurements": len(initial_queue), "remaining_pending_measurements": len(remaining),
                   "reviewed_measurements": len(initial_queue) - len(remaining),
                   "candidate_training_identities": len(ready),
                   "candidate_training_sequences": len({r["sequence"] for r in ready}),
                   "candidate_training_groups": len({r["group_id"] for r in ready}),
                   "newly_eligible_identities": sum(not r["original_training_identity"] for r in ready),
                   "active": sum(r["label"] == 1 for r in ready), "inactive": sum(r["label"] == 0 for r in ready),
                   "identity_statuses": dict(Counter(r["status"] for r in audit))},
        "status_transitions": [{"before": before, "after": after, "identities": n}
                               for (before, after), n in sorted(transitions.items())],
        "original_32_uM_replay_matches": True, "original_ready_labels_preserved": True,
        "heldout_classifier_csvs_opened": [], "models_fitted": False,
        "accuracy_improvement_estimated": False,
        "review_scope": "DBAASP note annotations; no claim of original-publication verification",
        "evaluation_limit": "v1 holdouts are already evaluated; a future performance claim requires new independent evaluation",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".amp-training-review-", dir=output.parent))
    try:
        tables = {"training_candidates.csv": (ready, list(audit[0])),
                  "training_audit.csv": (audit, list(audit[0])),
                  "initial_review_queue.csv": (initial_queue, QUEUE_HEADER),
                  "remaining_review_queue.csv": (remaining, QUEUE_HEADER)}
        for name, (items, fields) in tables.items():
            (temporary / name).write_text(csv_text(items, fields), encoding="utf-8")
        report["output_sha256"] = {name: fingerprint(temporary / name) for name in tables}
        (temporary / "review_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        temporary.replace(output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-dir", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--original-reviews", type=Path, default=ROOT / "data/curation/full_note_reviews.csv")
    parser.add_argument("--reviews", type=Path, default=ROOT / "data/curation/training_note_reviews_v2.csv")
    parser.add_argument("--results-root", type=Path, default=ROOT / "results/development")
    args = parser.parse_args()
    output, report = build(args.split_dir, args.original_reviews, args.reviews, args.results_root)
    print(json.dumps(report["counts"], indent=2))
    print("Results:", output)
    print("Training candidates only. No model fitted; no accuracy improvement claimed.")


if __name__ == "__main__":
    main()
