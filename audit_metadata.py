import json
from collections import Counter

from data_utils import RAW_DIR, is_ecoli_mic, load_json


records = [
    load_json(path)
    for path in sorted(RAW_DIR.glob("peptide_*.json"))
]

key_counts = Counter(
    key
    for record in records
    for key in record
)

print("Top-level keys — number of records containing each key:")
for key, count in sorted(key_counts.items()):
    print(f"  {key}: {count}/{len(records)}")

# Discover candidate fields without assuming their exact names.
keywords = (
    "terminus", "unusual", "amino", "bond",
    "complexity", "modif", "chiral", "stereo",
)
structure_keys = [
    key for key in sorted(key_counts)
    if any(word in key.casefold() for word in keywords)
]

for record in records:
    print(f"\nPEPTIDE {record.get('id')}")

    for key in structure_keys:
        if key not in record:
            print(f"  {key}: FIELD ABSENT")
        else:
            value = json.dumps(record[key], ensure_ascii=False)
            print(f"  {key}: {value}")

    print("  E. coli MIC assay conditions:")

    for measurement in record.get("targetActivities") or []:
        if not is_ecoli_mic(measurement):
            continue

        fields = (
            "id", "targetSpecies", "concentration", "unit",
            "medium", "cfu", "ph", "ionicStrength",
            "saltType", "note", "reference",
        )
        details = {
            field: measurement[field]
            for field in fields
            if field in measurement
        }
        print(json.dumps(details, indent=2, ensure_ascii=False))


print("\nAdditional identity metadata:")

for record in records:
    sequence = record.get("sequence") or ""

    print(f"\nPeptide {record.get('id')}")
    print("Sequence:", sequence)
    print("Calculated sequence length:", len(sequence))
    print("Reported sequence length:", record.get("sequenceLength"))

    for field in ("opms", "unusualAminoAcids", "noteReference"):
        if field in record:
            print(
                f"{field}:",
                json.dumps(record[field], ensure_ascii=False),
            )
        else:
            print(f"{field}: FIELD ABSENT")