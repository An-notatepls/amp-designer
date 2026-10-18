"""Build provenance-defined cohorts without changing the broad E. coli task.

Usage:
  python build_provenance.py PROCESSED_DIR --note-reviews REVIEWS.csv

Reads existing measurement tables; makes no network requests and never edits
its inputs. Provenance is inferred from explicit strain/notes annotations, not
verified against publications. Empty notes and wild-type/resistance descriptions
alone do not establish laboratory/reference origin. No train/test split is made
here: all cohorts must later share sequence-family group assignments.
"""

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from build_labels import (
    aggregate, apply_note_reviews, atomic_text, classifier_table, dump, name_of,
    write_csv,
)

NAMED_LAB_REFERENCE_TARGETS = {
    "Escherichia coli ATCC 25922", "Escherichia coli ATCC 8739",
    "Escherichia coli MG1655", "Escherichia coli K-12",
    "Escherichia coli K12 MG1655",
}
BOOL_KEYS = ("c_amidated", "classification_eligible", "regression_eligible")
NUMBER_KEYS = ("point_uM", "lower_uM", "upper_uM", "log10_mic_uM")
IDENTITY_HEADER = [
    "sequence", "c_amidated", "peptide_ids_json", "n_assay_groups", "active_groups",
    "inactive_groups", "mixed_labels", "needs_note_review", "has_unresolved_group",
    "label", "status", "model_eligible", "mixed_in_broad_dataset",
]


def provenance(row):
    note = name_of(json.loads(row["note_raw_json"]))
    named = row["strain"].strip() in NAMED_LAB_REFERENCE_TARGETS
    clinical = bool(re.search(r"\bclinical\b", note, re.IGNORECASE))
    laboratory = bool(re.search(r"\b(?:lab|laboratory) strain\b", note, re.IGNORECASE))
    if clinical and (named or laboratory):
        return "unspecified", True, "Clinical-origin note conflicts with named laboratory/reference target or laboratory note; source review required"
    if clinical:
        return "clinical_isolate", False, "Explicit clinical-origin note; concentration interpretation and treatment eligibility are separate decisions"
    if named:
        return "laboratory_reference", False, "Target matches the predefined named laboratory/reference list; not a claim of wild-type phenotype"
    if laboratory:
        return "laboratory_reference", False, "Explicit laboratory-strain note; exact strain identifier may be unspecified"
    return "unspecified", False, "No positive laboratory/reference or clinical-origin annotation; empty, WT, MDR, resistant or pathogenic descriptions alone do not establish origin"


def read_measurements(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "peptide_id", "measurement_id", "sequence", "strain", "medium",
            "medium_raw_json", "note_raw_json", "note_review_status",
            "exclusion_reasons", "op", "cfu_raw_json", "ph_raw_json",
            "saltType_raw_json", "ionicStrength_raw_json", "temperature_raw_json",
            *BOOL_KEYS, *NUMBER_KEYS,
        }
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Input measurement table lacks columns from the updated build_labels.py")
        rows = list(reader)
    seen = set()
    for row in rows:
        row["peptide_id"] = int(row["peptide_id"])
        row["measurement_id"] = int(row["measurement_id"])
        key = (row["peptide_id"], row["measurement_id"])
        if key in seen:
            raise ValueError(f"Duplicate measurement key: {key}")
        seen.add(key)
        for name in BOOL_KEYS:
            if row[name] not in {"True", "False"}:
                raise ValueError(f"Unexpected boolean in {name} for {key}")
            row[name] = row[name] == "True"
        for name in NUMBER_KEYS:
            row[name] = float(row[name]) if row[name] else None
            if row[name] is not None and not math.isfinite(row[name]):
                raise ValueError(f"Nonfinite numeric value in {name} for {key}")
        if row["note_review_status"] not in {"pending", "not_needed", "included", "excluded"}:
            raise ValueError(f"Unexpected note-review status for {key}")
        row.setdefault("note_review_reason", "")
    return rows


def identity(row):
    return row["sequence"], row["c_amidated"]


def build_cohorts(rows, threshold):
    for row in rows:
        origin, conflict, evidence = provenance(row)
        row.update(provenance=origin, provenance_conflict=conflict,
                   provenance_evidence=evidence)
    # Preserve the current broad policy: pending notes still block identities.
    broad = classifier_table(aggregate(rows, threshold))
    broadly_mixed = {identity(row) for row in broad if row["mixed_labels"]}
    cohorts = {"broad": broad}
    origin_evidence = defaultdict(set)
    conflicted = set()
    for row in rows:
        if row["classification_eligible"]:
            origin_evidence[identity(row)].add(row["provenance"])
            if row["provenance_conflict"]:
                conflicted.add(identity(row))
    for origin in ("laboratory_reference", "clinical_isolate", "unspecified"):
        selected = [r for r in rows if r["classification_eligible"]
                    and r["provenance"] == origin and not r["provenance_conflict"]
                    and r["note_review_status"] in {"not_needed", "included"}]
        cohorts[origin] = classifier_table(aggregate(selected, threshold))
    # "Clinical-only" means all observed base-eligible origin evidence is
    # clinical; it does not mean the peptide was never tested elsewhere.
    cohorts["clinical_only"] = [r for r in cohorts["clinical_isolate"]
        if origin_evidence[identity(r)] == {"clinical_isolate"}
        and identity(r) not in conflicted]
    for table in cohorts.values():
        for row in table:
            row["mixed_in_broad_dataset"] = identity(row) in broadly_mixed
    return cohorts


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("processed_dir", type=Path)
    parser.add_argument("--note-reviews", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=32.0)
    args = parser.parse_args()
    if not math.isfinite(args.threshold) or args.threshold <= 0:
        parser.error("Threshold must be finite and positive")
    source = args.processed_dir.expanduser().resolve()
    measurements = source / "ecoli_measurements.csv"
    report_path = source / "build_report.json"
    source_report = json.loads(report_path.read_text(encoding="utf-8"))
    rows = read_measurements(measurements)
    if len(rows) != source_report["counts"]["ecoli_mic_rows"]:
        raise ValueError("Measurement count differs from build_report.json")
    apply_note_reviews(rows, args.note_reviews)
    cohorts = build_cohorts(rows, args.threshold)
    policy = {
        "threshold_uM": args.threshold,
        "named_lab_reference_targets": sorted(NAMED_LAB_REFERENCE_TARGETS),
        "origin_rule": "Positive annotated origin; contradictory annotations require review; absence of notes is not reference evidence",
        "narrow_cohort_rule": "Select reviewed, base-eligible measurement rows first, then recompute labels within each provenance category",
        "broad_rule": "Unchanged broad policy including pending-note and mixed-identity exclusions",
        "clinical_only_definition": "Within the captured base-eligible measurements, only clinical-origin evidence is recorded; this is not proof of absence of other experiments",
        "evaluation_rule": "Candidate tables are not train/test partitions; keep sequence families disjoint across the eventual split and reserve clinical test examples before comparing models",
        "source_measurements_sha256": sha(measurements),
        "source_report_sha256": sha(report_path),
        "note_reviews_sha256": sha(args.note_reviews),
        "script_sha256": sha(Path(__file__)),
    }
    fingerprint = hashlib.sha256(dump(policy).encode()).hexdigest()[:12]
    output = source / "provenance" / fingerprint
    write_csv(output / "measurement_provenance.csv", rows,
              ["peptide_id", "measurement_id", "provenance", "provenance_conflict"])
    counts = {}
    for name, table in cohorts.items():
        ready = [r for r in table if r["model_eligible"]]
        write_csv(output / f"{name}_audit.csv", table, IDENTITY_HEADER)
        write_csv(output / f"{name}_candidates.csv", ready, IDENTITY_HEADER)
        counts[name] = {
            "identities": len(table), "statuses": dict(Counter(r["status"] for r in table)),
            "ready_identities": len(ready),
            "ready_sequences": len({r["sequence"] for r in ready}),
            "ready_labels": dict(Counter("active" if r["label"] else "inactive" for r in ready)),
            "ready_but_mixed_in_broad": sum(r["mixed_in_broad_dataset"] for r in ready),
        }
    counts["base_eligible_measurement_origins"] = dict(Counter(
        r["provenance"] for r in rows if r["classification_eligible"]))
    counts["base_eligible_provenance_conflicts"] = sum(
        r["provenance_conflict"] for r in rows if r["classification_eligible"])
    atomic_text(output / "provenance_report.json", dump({"policy": policy, "counts": counts}) + "\n")
    print(json.dumps(counts, indent=2))
    print("Wrote provenance cohorts to:", output)
    print("These are candidate cohorts. No test set is locked and no model is trained yet.")


if __name__ == "__main__":
    main()


