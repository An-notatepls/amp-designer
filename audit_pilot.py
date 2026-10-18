"""Audit only manifest IDs; report completeness, filter yield, and metadata states."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from data_utils import (STANDARD_AMINO_ACIDS, concentration_kind, is_ecoli_mic,
                        load_json, metadata_state, sequence_reasons,
                        target_name, validate_manifest, validate_record)

CHEMICAL_KEYS = ("complexity", "NTerminus", "nTerminus", "CTerminus", "cTerminus",
                 "unusualAminoAcids", "opms", "intrachainBonds",
                 "interchainBonds", "coordinationBonds")
ASSAY_KEYS = ("medium", "cfu", "ph", "ionicStrength", "saltType", "note", "reference")


def summarize(manifest, records_dir):
    ids = validate_manifest(manifest)
    counters = defaultdict(Counter)
    records, missing, invalid = [], [], []
    seq_ids, ecoli_ids, seq_ecoli_ids = set(), set(), set()
    for peptide_id in ids:
        path = records_dir / f"peptide_{peptide_id}.json"
        if not path.exists():
            missing.append(peptide_id)
            continue
        try:
            records.append(validate_record(load_json(path), peptide_id))
        except (OSError, ValueError) as error:
            invalid.append({"id": peptide_id, "error": str(error)})
    all_keys = sorted(set().union(*(record.keys() for record in records)))
    for record in records:
        peptide_id = record["id"]
        for key in all_keys:
            counters["top_level_key_states"][f"{key}: {metadata_state(record, key)}"] += 1
        for key in CHEMICAL_KEYS:
            counters["chemical_field_states"][f"{key}: {metadata_state(record, key)}"] += 1
        complexity = record.get("complexity") or {}
        label = complexity.get("name") if isinstance(complexity, dict) else str(complexity)
        counters["complexity"][label or "Missing"] += 1
        reasons = sequence_reasons(record)
        counters["sequence_rejection_reasons"].update(reasons)
        sequence = record.get("sequence") or ""
        if isinstance(sequence, str):
            counters["nonstandard_character_occurrences"].update(c for c in sequence if c not in STANDARD_AMINO_ACIDS)
        if not reasons:
            seq_ids.add(peptide_id)
        match = False
        for measurement in record["targetActivities"]:
            name = target_name(measurement)
            measure = measurement.get("activityMeasureValue") or "Missing"
            counters["measurement_types"][measure] += 1
            if "coli" in name.casefold():
                counters["all_coli_target_names"][name] += 1
            if not is_ecoli_mic(measurement):
                continue
            match = True
            counters["ecoli_mic_target_names"][name] += 1
            unit = measurement.get("unit") or {}
            counters["ecoli_mic_units"][unit.get("name") or "Missing"] += 1
            value = measurement.get("concentration")
            counters["ecoli_mic_formats"][concentration_kind(value)] += 1
            counters["ecoli_mic_value_types"][type(value).__name__] += 1
            for key in ASSAY_KEYS:
                counters["ecoli_mic_assay_states"][f"{key}: {metadata_state(measurement, key)}"] += 1
        if match:
            ecoli_ids.add(peptide_id)
            if not reasons:
                seq_ecoli_ids.add(peptide_id)
    n, selected, matched = len(records), len(ids), len(ecoli_ids)
    complete = not missing and not invalid
    return {"manifest_spec": manifest.get("spec"),
            "manifest_namespace": manifest.get("namespace"),
            "selected_records": selected, "valid_downloaded_records": n,
            "missing_ids": missing, "invalid_records": invalid, "complete": complete,
            "records_with_ecoli_mic": matched,
            "filter_yield_among_valid_downloads": matched / n if n else None,
            "filter_yield_of_selected_pilot": matched / selected if complete else None,
            "sequence_passing_records": len(seq_ids),
            "sequence_and_ecoli_mic_records": len(seq_ecoli_ids),
            "sequence_and_ecoli_mic_ids": sorted(seq_ecoli_ids),
            "counters": {key: dict(value) for key, value in sorted(counters.items())}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    report = summarize(load_json(args.manifest), args.manifest.parent / "records")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["complete"]:
        print("INCOMPLETE PILOT: resume downloads. The whole-pilot yield is not reported.")


if __name__ == "__main__":
    main()
