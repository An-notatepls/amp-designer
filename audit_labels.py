from collections import Counter

from data_utils import (
    RAW_DIR,
    concentration_kind,
    is_ecoli_mic,
    load_json,
    target_name,
)

coli_names = Counter()
concentration_types = Counter()

for path in RAW_DIR.glob("peptide_*.json"):
    record = load_json(path)

    for measurement in record.get("targetActivities") or []:
        name = target_name(measurement)

        if "coli" in name.casefold():
            coli_names[name] += 1

        value = measurement.get("concentration")
        concentration_types[type(value).__name__] += 1

print("Target names containing 'coli' — measurement counts:")
for name, count in coli_names.most_common():
    print(f"  {name}: {count}")

print("\nConcentration field types:")
for value_type, count in concentration_types.most_common():
    print(f"  {value_type}: {count}")

ecoli_mic_formats = Counter()

for path in RAW_DIR.glob("peptide_*.json"):
    record = load_json(path)

    for measurement in record.get("targetActivities") or []:
        if not is_ecoli_mic(measurement):
            continue

        value = measurement.get("concentration")
        kind = concentration_kind(value)
        ecoli_mic_formats[kind] += 1

        print(
            f"Peptide {record['id']} | "
            f"{target_name(measurement)} | "
            f"{value!r} | {kind}"
        )

print("\nE. coli MIC concentration formats:")
for kind, count in ecoli_mic_formats.most_common():
    print(f"  {kind}: {count}")