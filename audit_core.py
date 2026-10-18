import argparse
import json
from collections import Counter
from pathlib import Path

from data_utils import (
    concentration_kind,
    is_ecoli_mic,
    load_json,
    metadata_state,
    sequence_reasons,
    validate_manifest,
    validate_record,
)


def dump(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def field_label(record, key):
    """Preserve both the metadata state and its actual value."""
    return dump({
        "state": metadata_state(record, key),
        "value": record.get(key),
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()

    ids = validate_manifest(load_json(args.manifest))
    records_dir = args.manifest.parent / "records"

    chemical_keys = (
        "complexity",
        "CTerminus",
        "NTerminus",
        "unusualAminoAcids",
        "opms",
        "intrachainBonds",
        "interchainBonds",
        "coordinationBonds",
    )
    chemical_values = {key: Counter() for key in chemical_keys}
    states = Counter()
    media = Counter()
    alias_conflicts = Counter()
    odd = []

    core_count = 0
    mic_count = 0

    for peptide_id in ids:
        path = records_dir / f"peptide_{peptide_id}.json"
        record = validate_record(load_json(path), peptide_id)

        if sequence_reasons(record):
            continue

        mics = [
            measurement
            for measurement in record["targetActivities"]
            if is_ecoli_mic(measurement)
        ]
        if not mics:
            continue

        core_count += 1

        for key in chemical_keys:
            states[f"{key}: {metadata_state(record, key)}"] += 1
            chemical_values[key][field_label(record, key)] += 1

        for primary, alias in (
            ("CTerminus", "cTerminus"),
            ("NTerminus", "nTerminus"),
        ):
            if (
                primary in record
                and alias in record
                and record[primary] != record[alias]
            ):
                alias_conflicts[f"{primary}/{alias}"] += 1

        for measurement in mics:
            mic_count += 1
            media[field_label(measurement, "medium")] += 1

            value = measurement.get("concentration")
            if concentration_kind(value) == "unrecognized":
                odd.append({
                    "peptide_id": peptide_id,
                    "measurement_id": measurement.get("id"),
                    "concentration": value,
                    "unit": measurement.get("unit"),
                })

    print("Core candidates:", core_count)
    print("E. coli MIC measurements in core:", mic_count)

    print("\nChemical metadata states — record counts:")
    for label, count in states.most_common():
        print(f"  {label}: {count}")

    for key, counter in chemical_values.items():
        print(f"\n{key} values — record counts:")
        # Print every distinct chemical value, not just the top 15.
        for label, count in counter.most_common():
            print(f"  {label}: {count}")

    print("\nMedium — measurement counts, top 15:")
    for label, count in media.most_common(15):
        print(f"  {label}: {count}")
    print("Distinct medium representations:", len(media))

    print("\nTerminal alias conflicts:", dict(alias_conflicts))
    print("\nUnrecognized concentrations:")
    print(json.dumps(odd, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()