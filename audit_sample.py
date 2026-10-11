import json
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"

measurement_types = Counter()
mic_units = Counter()

record_count = 0
missing_sequences = 0
short_standard_sequences = 0
ecoli_measurements = 0
ecoli_peptide_ids = set()

standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

for path in sorted(RAW_DIR.glob("peptide_*.json")):
    record = json.loads(path.read_text(encoding="utf-8"))
    record_count += 1

    sequence = record.get("sequence") or ""

    if not sequence:
        missing_sequences += 1
    elif (
        10 <= len(sequence) <= 50
        and set(sequence).issubset(standard_amino_acids)
    ):
        short_standard_sequences += 1

    for measurement in record.get("targetActivities") or []:
        measure = measurement.get("activityMeasureValue") or "Missing"
        measurement_types[measure] += 1

        if measure != "MIC":
            continue

        unit = measurement.get("unit") or {}
        mic_units[unit.get("name") or "Missing"] += 1

        target = measurement.get("targetSpecies") or {}
        target_name = target.get("name") or ""

        if (
            target_name == "Escherichia coli"
            or target_name.startswith("Escherichia coli ")
        ):
            ecoli_measurements += 1
            ecoli_peptide_ids.add(record["id"])

print("Peptide records:", record_count)
print("Records without a sequence:", missing_sequences)
print("Records with standard sequences of length 10–50:",
      short_standard_sequences)

print("\nActivity measurements by type:")
for measure, count in measurement_types.most_common():
    print(f"  {measure}: {count}")

print("\nMIC measurements by unit:")
for unit, count in mic_units.most_common():
    print(f"  {unit}: {count}")

print("\nE. coli MIC measurements:", ecoli_measurements)
print("Distinct peptide IDs with E. coli MIC:",
      len(ecoli_peptide_ids))



print("\nStructure of records without a top-level sequence:")
for path in sorted(RAW_DIR.glob("peptide_*.json")):
    record = json.loads(path.read_text(encoding="utf-8"))

    if record.get("sequence"):
        continue

    chains = record.get("monomers") or []

    print(
        f"\nPeptide {record.get('id')} | "
        f"Complexity: {record.get('complexity')} | "
        f"Component chains: {len(chains)}"
    )

    for chain in chains:
        print(
            f"  Chain ID: {chain.get('id')} | "
            f"Sequence: {chain.get('sequence')}"
        )