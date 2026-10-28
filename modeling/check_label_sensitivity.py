"""Relabel frozen training identities at 16/32/64 uM without fitting models.

Usage: python modeling/check_label_sensitivity.py SPLIT_DIR
       --note-reviews data/curation/full_note_reviews.csv

Reuses original converted bounds and reviews, within frozen training-assigned
families. Reports both the original ready training identities and the broader
training-family audit pool, including newly eligible identities at another
cutoff. Only selected rows are interpreted/aggregated. No test.csv or clinical
split CSV is opened. No alternative-threshold model performance is estimated.
"""

import argparse
import csv
import hashlib
import importlib
import json
import shutil
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import train_baseline as base

VERSION = "training-label-sensitivity-v1"
THRESHOLDS = (16, 32, 64)


def label_modules():
    root = Path(__file__).resolve().parent
    if not (root / "build_labels.py").is_file():
        root = root.parent
    if not (root / "build_labels.py").is_file() or not (root / "build_provenance.py").is_file():
        raise ValueError("Keep build_labels.py and build_provenance.py in the project root")
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return importlib.import_module("build_labels"), importlib.import_module("build_provenance")


def select_training_rows(measurement_path, review_path, training, temporary):
    """Filter string metadata before the existing parser interprets MIC bounds."""
    allowed = {}
    for row in training:
        for pid in json.loads(row["peptide_ids_json"]):
            key = base.identity(row)
            if pid in allowed and allowed[pid] != key:
                raise ValueError("A source peptide ID belongs to multiple modeled identities")
            allowed[pid] = key
    selected, seen_ids = [], set()
    with measurement_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        if not {"peptide_id", "measurement_id", "sequence", "c_amidated"} <= set(header):
            raise ValueError("Measurement source lacks identity columns")
        for row in reader:
            pid = int(row["peptide_id"])
            if pid not in allowed:
                continue
            if base.identity(row) != allowed[pid]:
                raise ValueError("Source measurement identity differs from frozen training metadata")
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed selected measurement row")
            selected.append(row)
            seen_ids.add(pid)
    if seen_ids != set(allowed):
        raise ValueError("Some frozen training source peptide IDs have no measurement rows")
    selected_path = temporary / "training_measurements.csv"
    with selected_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(selected)
    keys = {(r["peptide_id"], r["measurement_id"]) for r in selected}
    filtered_reviews = []
    with review_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        review_header = reader.fieldnames or []
        if not {"peptide_id", "measurement_id", "note", "decision", "review_reason"} <= set(review_header):
            raise ValueError("Review source lacks required columns")
        for row in reader:
            if int(row["peptide_id"]) not in allowed:
                continue
            if (row["peptide_id"], row["measurement_id"]) not in keys:
                raise ValueError("Unknown review measurement for a selected training peptide")
            filtered_reviews.append(row)
    selected_reviews = temporary / "training_reviews.csv"
    with selected_reviews.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=review_header)
        writer.writeheader()
        writer.writerows(filtered_reviews)
    return selected_path, selected_reviews, len(selected), len(allowed), len(filtered_reviews)


def training_roster(cohort, split_dir, frozen_training):
    groups = {r["sequence"]: r for r in base.read_csv(split_dir / "sequence_groups.csv")}
    selected = []
    with (cohort / "broad_audit.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            group = groups.get(row["sequence"])
            if group is None:
                raise ValueError("Audit sequence absent from the frozen group index")
            if group["group_partition"] != "train":
                continue
            if row["model_eligible"] not in {"True", "False"} or row["c_amidated"] not in {"True", "False"}:
                raise ValueError("Invalid audit identity/eligibility metadata")
            selected.append({**row, "group_id": group["group_id"]})
    ready = {base.identity(r): r for r in selected if r["model_eligible"] == "True"}
    frozen = {base.identity(r): r for r in frozen_training}
    if len({base.identity(r) for r in selected}) != len(selected) or set(ready) != set(frozen):
        raise ValueError("Training-family audit's ready identities differ from frozen training")
    for key, row in ready.items():
        if row["label"] != frozen[key]["label"] or row["status"] != frozen[key]["status"] or row["group_id"] != frozen[key]["group_id"]:
            raise ValueError("Training-family audit disagrees with frozen labels/groups")
        if json.loads(row["peptide_ids_json"]) != json.loads(frozen[key]["peptide_ids_json"]):
            raise ValueError("Training-family audit disagrees with frozen source peptide IDs")
        for field in ("n_assay_groups", "active_groups", "inactive_groups"):
            if field in row and field in frozen[key] and int(row[field]) != int(frozen[key][field]):
                raise ValueError(f"Training-family audit disagrees with frozen {field}")
    return selected


def rebuild(training, measurements, labels, frozen_training=None):
    expected = {(r["sequence"], r["c_amidated"] == "True"): r for r in training}
    frozen_training = training if frozen_training is None else frozen_training
    frozen_keys = {(r["sequence"], r["c_amidated"] == "True") for r in frozen_training}
    tables = {}
    for threshold in THRESHOLDS:
        audit = labels.classifier_table(labels.aggregate(measurements, threshold))
        table = {(r["sequence"], r["c_amidated"]): r for r in audit}
        if set(table) - set(expected):
            raise ValueError("Rebuild produced an identity outside frozen training scope")
        tables[threshold] = table
    if set(tables[32]) != set(expected):
        raise ValueError("32-uM replay did not recover every frozen training identity")
    for key, source in expected.items():
        replay = tables[32][key]
        eligible = source["model_eligible"] == "True"
        source_label = int(source["label"]) if eligible else None
        if replay["model_eligible"] != eligible or replay["label"] != source_label or replay["status"] != source["status"]:
            raise ValueError("32-uM replay differs from frozen training labels; preserve the test seal and inspect the inputs")
        if json.loads(replay["peptide_ids_json"]) != json.loads(source["peptide_ids_json"]):
            raise ValueError("32-uM replay uses different source peptide IDs")
        for field in ("n_assay_groups", "active_groups", "inactive_groups"):
            if field in source and replay[field] != int(source[field]):
                raise ValueError(f"32-uM replay differs in {field}")
    details, summary, transitions = [], [], []
    for key, source in expected.items():
        row = {"sequence": source["sequence"], "c_amidated": source["c_amidated"],
               "group_id": source["group_id"], "source_peptide_ids_json": source["peptide_ids_json"],
               "original_training_identity": key in frozen_keys, "original_status_32": source["status"],
               "frozen_label_32": int(source["label"]) if source["model_eligible"] == "True" else None}
        for threshold in THRESHOLDS:
            result = tables[threshold].get(key)
            row[f"status_{threshold}"] = result["status"] if result else "no_remaining_groups"
            row[f"label_{threshold}"] = result["label"] if result and result["model_eligible"] else None
        # If both labels are definite, raising the MIC cutoff cannot turn active into inactive.
        if row["label_16"] is not None and row["frozen_label_32"] is not None and row["label_16"] > row["frozen_label_32"]:
            raise ValueError("Non-monotone definite label at 16 uM")
        if row["label_64"] is not None and row["frozen_label_32"] is not None and row["label_64"] < row["frozen_label_32"]:
            raise ValueError("Non-monotone definite label at 64 uM")
        details.append(row)
    for scope in ("original_ready_training", "all_training_assigned_audit_identities"):
        scoped = [r for r in details if r["original_training_identity"]] if scope == "original_ready_training" else details
        for threshold in THRESHOLDS:
            counts = Counter(r[f"status_{threshold}"] for r in scoped)
            ready = [r for r in scoped if r[f"label_{threshold}"] is not None]
            summary.append({"scope": scope, "threshold_uM": threshold, "identities_in_scope": len(scoped),
                            "original_training_identities": len(frozen_training),
                            "eligible_identities": len(ready), "eligible_sequences": len({r['sequence'] for r in ready}),
                            "eligible_groups": len({r['group_id'] for r in ready}),
                            "active": counts["active"], "inactive": counts["inactive"],
                            "active_fraction_among_eligible": counts["active"] / len(ready) if ready else None,
                            "mixed": counts["mixed"], "ambiguous": counts["ambiguous"],
                            "needs_note_review": counts["needs_note_review"], "no_remaining_groups": counts["no_remaining_groups"],
                            "definite_label_changes_from_32": sum(r["frozen_label_32"] is not None and r[f"label_{threshold}"] != r["frozen_label_32"] for r in ready),
                            "newly_eligible_from_32": sum(r["frozen_label_32"] is None for r in ready),
                            "excluded_from_original_training": sum(r["original_training_identity"] and r[f"label_{threshold}"] is None for r in scoped)})
            changes = Counter((r["original_status_32"], r[f"status_{threshold}"]) for r in scoped)
            transitions.extend({"scope": scope, "threshold_uM": threshold, "status_at_frozen_32": before,
                                "status_at_threshold": after, "identities": n} for (before, after), n in sorted(changes.items()))
    return details, summary, transitions


def check(split_dir, note_reviews, results_root=Path("results/label_sensitivity")):
    split_dir, note_reviews = Path(split_dir).resolve(), Path(note_reviews).resolve()
    training, split = base.load_training(split_dir)
    cohort = split_dir.parent.parent
    processed = cohort.parent.parent
    for name, expected in split["spec"]["source_sha256"].items():
        if Path(name).name != name or base.sha_file(cohort / name) != expected:
            raise ValueError(f"Original cohort input changed: {name}")
    provenance_report = json.loads((cohort / "provenance_report.json").read_text(encoding="utf-8"))
    policy = provenance_report["policy"]
    source = processed / "ecoli_measurements.csv"
    build_report_path = processed / "build_report.json"
    build_report = json.loads(build_report_path.read_text(encoding="utf-8"))
    checks = {source: policy["source_measurements_sha256"], build_report_path: policy["source_report_sha256"],
              note_reviews: policy["note_reviews_sha256"]}
    for path, expected in checks.items():
        if base.sha_file(path) != expected:
            raise ValueError(f"Input differs from the recorded provenance build: {path.name}")
    labels, provenance = label_modules()
    if base.sha_file(Path(labels.__file__)) != build_report["policy"]["script_sha256"]:
        raise ValueError("build_labels.py differs from the recorded source; use the original labeling code")
    if base.sha_file(Path(provenance.__file__)) != policy["script_sha256"]:
        raise ValueError("build_provenance.py differs from the recorded source")
    roster = training_roster(cohort, split_dir, training)
    spec = {"version": VERSION, "thresholds_uM": list(THRESHOLDS),
            "scope": "original ready training and all audit identities in frozen training-assigned families; no held-out families",
            "split_manifest_sha256": base.sha_file(split_dir / "split_manifest.json"),
            "train_sha256": base.sha_file(split_dir / "train.csv"),
            "measurements_sha256": base.sha_file(source), "review_sha256": base.sha_file(note_reviews),
            "provenance_report_sha256": base.sha_file(cohort / "provenance_report.json"),
            "build_report_sha256": base.sha_file(build_report_path),
            "script_sha256": base.sha_file(Path(__file__)),
            "environment": base.versions(),
            "label_source_sha256": base.sha_file(Path(labels.__file__)),
            "provenance_source_sha256": base.sha_file(Path(provenance.__file__)),
            "data_utils_sha256": base.sha_file(Path(sys.modules['data_utils'].__file__))}
    digest = hashlib.sha256(base.canonical(spec).encode()).hexdigest()[:16]
    output = Path(results_root).resolve() / f"{VERSION}-{digest}"
    if output.exists():
        report = json.loads((output / "sensitivity_report.json").read_text(encoding="utf-8"))
        if report["spec"] != spec:
            raise ValueError("Existing sensitivity specification differs")
        for name, expected in report["output_sha256"].items():
            if base.sha_file(output / name) != expected:
                raise ValueError(f"Existing sensitivity output changed: {name}")
        return output, report
    with tempfile.TemporaryDirectory(prefix="amp-training-labels-") as directory:
        selected, reviews, n_rows, n_ids, n_reviews = select_training_rows(source, note_reviews, roster, Path(directory))
        rows = provenance.read_measurements(selected)
        labels.apply_note_reviews(rows, reviews)
        details, summary, transitions = rebuild(roster, rows, labels, frozen_training=training)
    report = {"spec": spec, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
              "training_source_peptide_ids": n_ids, "training_measurement_rows": n_rows, "training_review_rows": n_reviews,
              "summary": summary, "training_family_audit_identities": len(roster), "replayed_32_uM_labels_match": True,
              "frozen_primary_threshold_uM": 32, "frozen_primary_model_changed": False,
              "heldout_classifier_csvs_opened": [], "models_fitted": False,
              "limitations": ["Restricted to the original audit identities in training-assigned families; no conclusions about held-out family label changes or new sequences outside the index.",
                              "Full-cohort audit/measurement CSVs are streamed for identity/group filtering; only training-family rows are interpreted and aggregated. Full-file hashes bind the source snapshots.",
                              "No test/clinical classifier table is opened and no performance at an alternative MIC threshold is estimated.",
                              "This diagnostic does not select a replacement activity threshold, feature set or model."],
              "source_attribution": "Derived from DBAASP data, https://dbaasp.org; see project citation and data terms."}
    payloads = {"training_identity_labels.csv": base.csv_text(details), "threshold_summary.csv": base.csv_text(summary),
                "label_transitions.csv": base.csv_text(transitions)}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".label-sensitivity-", dir=output.parent))
    try:
        for name, content in payloads.items():
            (temporary / name).write_text(content, encoding="utf-8", newline="")
        report["output_sha256"] = {p.name: base.sha_file(p) for p in sorted(temporary.iterdir())}
        (temporary / "sensitivity_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.rename(output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split_dir", type=Path)
    parser.add_argument("--note-reviews", required=True, type=Path)
    parser.add_argument("--results-root", type=Path, default=Path("results/label_sensitivity"))
    args = parser.parse_args()
    output, report = check(args.split_dir, args.note_reviews, args.results_root)
    print(json.dumps(report["summary"], indent=2))
    print("32-uM training replay: exact label match")
    print("Results:", output)
    print("Training-family label check; primary cutoff remains 32 uM. No models fitted; test.csv and clinical split CSVs remain unopened.")


if __name__ == "__main__":
    main()
