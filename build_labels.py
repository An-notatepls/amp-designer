"""Build auditable E. coli MIC tables from a validated pilot manifest.

Usage: python build_labels.py PATH/TO/manifest.json --threshold 32

This is a conservative starting policy, not a validated activity definition.
Empty chemistry annotations mean 'no modification recorded', not confirmed
unmodified chemistry. Amidation is allowed and represented explicitly.
Regression uses single numeric MICs only. Classification may also use ranges
and ordinary inequality bounds, without replacing them with point estimates.
Assay groups remain separate; these are NOT one-label-per-peptide tables.
The classifier audit has one row per sequence plus amidation state. Only
definite, consistent identities with no pending note review enter the classifier
dataset. This deliberately narrows the first baseline and may introduce bias.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import statistics
import tempfile
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from data_utils import (
    PROJECT_ROOT, is_ecoli_mic, load_json, target_name,
    validate_manifest, validate_record,
)

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 10, 50
ALLOW_AMIDATION = True
AMIDATION_DELTA = -0.9848  # approximate average mass difference, g/mol
CHEMISTRY_KEYS = (
    "unusualAminoAcids", "intrachainBonds", "interchainBonds",
    "coordinationBonds",
)
# A provisional inclusion list, not a claim that these are all wild-type.
ALLOWED_STRAINS = {
    "Escherichia coli", "Escherichia coli ATCC 25922",
    "Escherichia coli ATCC 8739", "Escherichia coli MG1655",
    "Escherichia coli K-12", "Escherichia coli K12 MG1655",
}
ALLOWED_MEDIA = {"MHB", "CAMHB"}  # retained as separate assay groups
UNIT_KIND = {"μm": "molar", "μg/ml": "mass", "mg/l": "mass"}
ASSAY_KEYS = ("cfu", "ph", "saltType", "ionicStrength", "temperature")
NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
SINGLE = re.compile(NUMBER)
BOUND = re.compile(rf"(>=|<=|>|<|≥|≤)\s*({NUMBER})")
RANGE = re.compile(rf"({NUMBER})\s*[-–]\s*({NUMBER})")
PLUS_MINUS = re.compile(rf"({NUMBER})\s*±\s*({NUMBER})")


def dump(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)


def name_of(value):
    if isinstance(value, dict):
        value = value.get("name")
    return value.strip() if isinstance(value, str) else ""


def normalize_unit(value):
    return unicodedata.normalize("NFKC", value).strip().casefold()


def terminal(record, primary, alias):
    """Return annotation and an issue; do not hide conflicting aliases."""
    if primary in record and alias in record and record[primary] != record[alias]:
        return None, f"conflicting {primary} aliases"
    if primary not in record and alias not in record:
        return None, f"missing {primary} metadata"
    return record.get(primary) if primary in record else record[alias], None


def peptide_info(record):
    reasons = []
    sequence = record.get("sequence") or ""
    if not isinstance(sequence, str):
        sequence = ""
    if not sequence:
        reasons.append("no top-level sequence")
    else:
        if not MIN_LEN <= len(sequence) <= MAX_LEN:
            reasons.append("length outside range")
        if not set(sequence) <= STANDARD:
            reasons.append("non-standard residues")
    if name_of(record.get("complexity")) != "Monomer":
        reasons.append("not an annotated monomer")
    nterm, nissue = terminal(record, "NTerminus", "nTerminus")
    cterm, cissue = terminal(record, "CTerminus", "cTerminus")
    reasons.extend(issue for issue in (nissue, cissue) if issue)
    if nterm:
        reasons.append("recorded N-terminal modification")
    amidated = name_of(cterm) == "AMD"
    if cterm and not (ALLOW_AMIDATION and amidated):
        reasons.append("C-terminal modification outside policy")
    for key in CHEMISTRY_KEYS:
        if key not in record:
            reasons.append(f"missing {key} metadata")
        elif record[key]:
            reasons.append(f"recorded {key}")
    return sequence, amidated, reasons


def parse_concentration(raw):
    """Preserve bounds and uncertainty; never impute a range midpoint."""
    text = "" if raw is None else str(raw).strip()
    out = {"format": "unrecognized", "op": None, "point": None,
           "lower": None, "upper": None, "uncertainty": None}
    if not text:
        out["format"] = "missing"
    elif text.startswith((">>", "<<")):
        out["format"] = "ambiguous_bound"
    elif SINGLE.fullmatch(text):
        value = float(text)
        out.update(format="single_number", op="=", point=value,
                   lower=value, upper=value)
    elif match := BOUND.fullmatch(text):
        op = {"≥": ">=", "≤": "<="}.get(match[1], match[1])
        out.update(format="censored", op=op)
        out["lower" if op.startswith(">") else "upper"] = float(match[2])
    elif match := RANGE.fullmatch(text):
        out.update(format="range", lower=float(match[1]), upper=float(match[2]))
    elif match := PLUS_MINUS.fullmatch(text):
        out.update(format="plus_minus", point=float(match[1]),
                   uncertainty=float(match[2]))
    values = [out[key] for key in ("point", "lower", "upper")
              if out[key] is not None]
    uncertainty = out["uncertainty"]
    if (any(not math.isfinite(v) or v <= 0 for v in values)
            or (uncertainty is not None and not math.isfinite(uncertainty))
            or (out["format"] == "range" and out["lower"] > out["upper"])):
        # Avoid nonfinite numbers leaking into JSON reports or calculations.
        out.update(format="invalid_numeric", point=None, lower=None,
                   upper=None, uncertainty=None)
    return out


def measurement_label(row, threshold):
    """1 active, 0 inactive, None unresolved from the reported bounds."""
    upper, lower = row["upper_uM"], row["lower_uM"]
    if upper is not None and upper <= threshold:
        return 1
    if lower is not None:
        if lower > threshold or (lower == threshold and row["op"] == ">"):
            return 0
    return None


def build_measurement(record, measurement, source_path, sequence, amidated,
                      peptide_reasons, weight):
    parsed = parse_concentration(measurement.get("concentration"))
    strain = target_name(measurement).strip()
    medium = name_of(measurement.get("medium"))
    unit = name_of(measurement.get("unit"))
    kind = UNIT_KIND.get(normalize_unit(unit))
    reasons = list(peptide_reasons)
    if strain not in ALLOWED_STRAINS:
        reasons.append("strain outside inclusion list")
    if medium not in ALLOWED_MEDIA:
        reasons.append("medium outside inclusion list")
    if kind is None:
        reasons.append("unrecognized unit")
    # +/- represents reported uncertainty, not an interval of known MIC bounds.
    if parsed["format"] not in {"single_number", "range", "censored"}:
        reasons.append(f"concentration {parsed['format']}")
    row = {
        "peptide_id": record["id"], "measurement_id": measurement.get("id"),
        "source_file": str(source_path), "sequence": sequence,
        "c_amidated": amidated, "mw_g_mol": weight, "strain": strain,
        "medium": medium, "unit_raw": unit,
        "raw_concentration": measurement.get("concentration"),
        "format": parsed["format"], "op": parsed["op"],
        "point_raw": parsed["point"], "lower_raw": parsed["lower"],
        "upper_raw": parsed["upper"], "uncertainty_raw": parsed["uncertainty"],
        "point_uM": None, "lower_uM": None, "upper_uM": None,
        "log10_mic_uM": None, "classification_eligible": not reasons,
        "regression_eligible": not reasons and parsed["format"] == "single_number",
        "exclusion_reasons": dump(reasons),
        "reference_raw": dump(measurement.get("reference")),
        "measurement_raw_json": dump(measurement),
        "medium_raw_json": dump(measurement.get("medium")),
        "note_raw_json": dump(measurement.get("note")),
        "note_review_status": "pending" if name_of(measurement.get("note")) else "not_needed",
        "note_review_reason": "",
    }
    row.update({f"{key}_raw_json": dump(measurement.get(key)) for key in ASSAY_KEYS})
    if not reasons:
        factor = 1.0 if kind == "molar" else 1000.0 / weight
        for key in ("point", "lower", "upper"):
            if parsed[key] is not None:
                row[f"{key}_uM"] = parsed[key] * factor
        converted = [row[f"{key}_uM"] for key in ("point", "lower", "upper")
                     if row[f"{key}_uM"] is not None]
        if any(not math.isfinite(value) or value <= 0 for value in converted):
            reasons.append("non-finite or non-positive converted concentration")
            row.update(classification_eligible=False, regression_eligible=False,
                       point_uM=None, lower_uM=None, upper_uM=None,
                       exclusion_reasons=dump(reasons))
        elif row["regression_eligible"]:
            row["log10_mic_uM"] = math.log10(row["point_uM"])
    return row


def aggregate(rows, threshold):
    """Group within recorded assay context, never silently across all strains."""
    groups = defaultdict(list)
    context_keys = ["medium_raw_json", "note_raw_json"] + [f"{key}_raw_json" for key in ASSAY_KEYS]
    for row in rows:
        if row["classification_eligible"] and row["note_review_status"] != "excluded":
            key = (row["peptide_id"], row["strain"],
                   *(row[k] for k in context_keys))
            groups[key].append(row)
    result = []
    for group in groups.values():
        labels = [measurement_label(row, threshold) for row in group]
        known = set(labels) - {None}
        label = None
        if len(known) > 1:
            status = "conflict"
        elif None in labels:
            status = "ambiguous"
        else:
            label = next(iter(known))
            status = "active" if label == 1 else "inactive"
        points = [r["log10_mic_uM"] for r in group if r["regression_eligible"]]
        first = group[0]
        output = {key: first[key] for key in
                  ["peptide_id", "sequence", "c_amidated", "strain", "medium"]
                  + context_keys}
        output.update(
            n_measurements=len(group), n_single_numeric=len(points),
            median_log10_mic_uM=statistics.median(points) if points else None,
            spread_log10=max(points) - min(points) if len(points) > 1 else None,
            label=label, label_status=status, threshold_uM=threshold,
            measurement_ids_json=dump([r["measurement_id"] for r in group]),
            needs_note_review=any(r["note_review_status"] == "pending" for r in group),
        )
        result.append(output)
    return result


def classifier_table(assay_groups):
    """No majority vote: opposing evidence always produces a mixed flag."""
    identities = defaultdict(list)
    for group in assay_groups:
        identities[(group["sequence"], group["c_amidated"])].append(group)
    rows = []
    for (sequence, amidated), groups in identities.items():
        known = {g["label"] for g in groups if g["label"] is not None}
        mixed = len(known) > 1 or any(g["label_status"] == "conflict" for g in groups)
        pending = any(g["needs_note_review"] for g in groups)
        unresolved = any(g["label"] is None for g in groups)
        if mixed:
            status = "mixed"
        elif pending:
            status = "needs_note_review"
        elif unresolved:
            status = "ambiguous"
        else:
            status = "active" if known == {1} else "inactive"
        eligible = status in {"active", "inactive"}
        rows.append({
            "sequence": sequence, "c_amidated": amidated,
            "peptide_ids_json": dump(sorted({g["peptide_id"] for g in groups})),
            "n_assay_groups": len(groups),
            "active_groups": sum(g["label"] == 1 for g in groups),
            "inactive_groups": sum(g["label"] == 0 for g in groups),
            "mixed_labels": mixed, "needs_note_review": pending,
            "has_unresolved_group": unresolved,
            "label": next(iter(known)) if eligible else None,
            "status": status, "model_eligible": eligible,
        })
    return rows


def apply_note_reviews(rows, review_path):
    """Reviews are keyed to a measurement AND its exact current note."""
    if review_path is None:
        return
    by_key = {(str(r["peptide_id"]), str(r["measurement_id"])): r for r in rows}
    seen = set()
    with review_path.open(encoding="utf-8", newline="") as handle:
        reviews = csv.DictReader(handle)
        required = {"peptide_id", "measurement_id", "note", "decision", "review_reason"}
        if not required <= set(reviews.fieldnames or []):
            raise ValueError("Note review CSV is missing required columns")
        for review in reviews:
            key = (review["peptide_id"], review["measurement_id"])
            if key in seen or key not in by_key:
                raise ValueError(f"Duplicate or unknown note-review key: {key}")
            seen.add(key)
            row = by_key[key]
            if review["note"] != name_of(json.loads(row["note_raw_json"])):
                raise ValueError(f"Note changed since review for {key}")
            decision = review["decision"].strip().casefold()
            if not decision:  # blank rows remain pending
                continue
            if decision not in {"include", "exclude"} or not review["review_reason"].strip():
                raise ValueError(f"Review for {key} needs include/exclude and a reason")
            row["note_review_status"] = "included" if decision == "include" else "excluded"
            row["note_review_reason"] = review["review_reason"].strip()


def note_review_table(rows):
    return [{
        "peptide_id": r["peptide_id"], "measurement_id": r["measurement_id"],
        "note": name_of(json.loads(r["note_raw_json"])), "strain": r["strain"],
        "medium": r["medium"], "reference_raw": r["reference_raw"],
        "raw_concentration": r["raw_concentration"], "unit_raw": r["unit_raw"],
        "current_review_status": r["note_review_status"],
        "decision": "", "review_reason": "",
    } for r in rows if r["classification_eligible"] and r["note_review_status"] == "pending"]


def atomic_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                     dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)  # processed outputs may be regenerated
    finally:
        temporary.unlink(missing_ok=True)


def write_csv(path, rows, empty_header):
    stream = io.StringIO(newline="")
    fields = list(rows[0]) if rows else empty_header
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, stream.getvalue())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--threshold", type=float, default=32.0)
    parser.add_argument("--note-reviews", type=Path,
                        help="Optional completed note review CSV; keep originals unchanged")
    args = parser.parse_args()
    if not math.isfinite(args.threshold) or args.threshold <= 0:
        parser.error("--threshold must be finite and positive")
    manifest_path = args.manifest.expanduser().resolve()
    manifest = load_json(manifest_path)
    ids = validate_manifest(manifest)
    records_dir = manifest_path.parent / "records"
    records, input_hashes, seen_ids = [], {}, set()
    # Validate EVERY selected file before writing anything. Extra files ignored.
    for pid in ids:
        path = records_dir / f"peptide_{pid}.json"
        record = validate_record(load_json(path), pid)
        for measurement in record["targetActivities"]:
            mid = measurement.get("id")
            if mid is None:
                raise ValueError(f"Missing measurement ID in peptide {pid}")
            key = (pid, mid)
            if key in seen_ids:
                raise ValueError(f"Duplicate measurement key {key}")
            seen_ids.add(key)
        records.append((path, record))
        input_hashes[str(pid)] = hashlib.sha256(path.read_bytes()).hexdigest()
    peptide_rows, rows = [], []
    for path, record in records:
        sequence, amidated, reasons = peptide_info(record)
        mics = [m for m in record["targetActivities"] if is_ecoli_mic(m)]
        if not mics:
            reasons.append("no E. coli MIC")
        weight = None
        if not reasons:
            weight = ProteinAnalysis(sequence).molecular_weight()
            if amidated:
                weight += AMIDATION_DELTA
        chemical_keys = ["complexity", "NTerminus", "nTerminus", "CTerminus",
                         "cTerminus", *CHEMISTRY_KEYS]
        peptide_rows.append({
            "peptide_id": record["id"], "sequence": sequence,
            "length": len(sequence), "c_amidated": amidated,
            "mw_g_mol": weight, "n_ecoli_mic": len(mics),
            "eligible": not reasons, "exclusion_reasons": dump(reasons),
            "chemical_metadata_json": dump({k: record[k] for k in chemical_keys
                                            if k in record}),
            "source_file": str(path), "source_sha256": input_hashes[str(record["id"])],
        })
        rows.extend(build_measurement(record, m, path, sequence, amidated, reasons,
                                      weight) for m in mics)
    apply_note_reviews(rows, args.note_reviews)
    groups = aggregate(rows, args.threshold)
    classifier_audit = classifier_table(groups)
    classifier_ready = [r for r in classifier_audit if r["model_eligible"]]
    note_reviews = note_review_table(rows)
    policy = {
        "length": [MIN_LEN, MAX_LEN], "allow_amidation": ALLOW_AMIDATION,
        "amidation_delta_average_g_mol": AMIDATION_DELTA,
        "strains": sorted(ALLOWED_STRAINS), "media": sorted(ALLOWED_MEDIA),
        "threshold_uM": args.threshold, "regression_formats": ["single_number"],
        "classification_formats": ["single_number", "range", "censored"],
        "assay_keys": ["strain", "medium", "note", *ASSAY_KEYS],
        "classifier_policy": "one row per sequence+amidation; exclude mixed, unresolved, or pending-note identities; no majority vote",
        "note_policy": "nonempty notes require manual review; differing notes never pooled; notes do not prove causation or co-treatment",
        "note_review_sha256": hashlib.sha256(args.note_reviews.read_bytes()).hexdigest() if args.note_reviews else None,
        "note_review_source_file": str(args.note_reviews.resolve()) if args.note_reviews else None,
        "empty_chemistry_assumption": "no modification recorded; not verified absent",
        "missing_chemistry_policy": "exclude absent fields; accept present null/empty",
        "classification_policy": "conflict if opposing definite labels; otherwise any unresolved row makes group ambiguous",
        "regression_policy": "median log10 MIC of single numeric rows within assay group; inspect classification conflicts separately",
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "record_sha256": input_hashes,
    }
    fingerprint = hashlib.sha256(dump(policy).encode()).hexdigest()[:12]
    out_dir = PROJECT_ROOT / "data" / "processed" / manifest_path.parent.name / fingerprint
    write_csv(out_dir / "peptide_eligibility.csv", peptide_rows, ["peptide_id", "eligible"])
    write_csv(out_dir / "ecoli_measurements.csv", rows,
              ["peptide_id", "measurement_id", "classification_eligible", "regression_eligible"])
    write_csv(out_dir / "ecoli_assay_labels.csv", groups,
              ["peptide_id", "strain", "medium", "label", "label_status", "median_log10_mic_uM"])
    classifier_header = ["sequence", "c_amidated", "peptide_ids_json", "n_assay_groups",
                         "active_groups", "inactive_groups", "mixed_labels", "needs_note_review",
                         "has_unresolved_group", "label", "status", "model_eligible"]
    write_csv(out_dir / "ecoli_classifier_audit.csv", classifier_audit, classifier_header)
    write_csv(out_dir / "ecoli_classifier_dataset.csv", classifier_ready, classifier_header)
    write_csv(out_dir / "note_review.csv", note_reviews,
              ["peptide_id", "measurement_id", "note", "strain", "medium", "reference_raw",
               "raw_concentration", "unit_raw", "current_review_status", "decision", "review_reason"])
    counts = {
        "selected_records": len(ids), "eligible_peptides": sum(r["eligible"] for r in peptide_rows),
        "ecoli_mic_rows": len(rows),
        "classification_rows": sum(r["classification_eligible"] for r in rows),
        "regression_rows": sum(r["regression_eligible"] for r in rows),
        "assay_group_statuses": dict(Counter(r["label_status"] for r in groups)),
        "classifier_identity_statuses": dict(Counter(r["status"] for r in classifier_audit)),
        "classifier_ready_identities": len(classifier_ready),
        "classifier_ready_sequences": len({r["sequence"] for r in classifier_ready}),
        "pending_note_reviews": len(note_reviews),
        "peptide_exclusions_all_reasons": dict(Counter(
            reason for row in peptide_rows for reason in json.loads(row["exclusion_reasons"]))),
    }
    atomic_text(out_dir / "build_report.json", dump({"policy": policy, "counts": counts}) + "\n")
    print(json.dumps(counts, indent=2))
    print("Wrote tables to:", out_dir)
    print("Rows in the label table are assay groups, not unique peptides.")
    print("Use ecoli_classifier_dataset.csv for the first classifier; keep all rows of a sequence in the same split.")
    print("Regression summaries remain provisional until assay notes are reviewed.")


if __name__ == "__main__":
    main()
