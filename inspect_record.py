import json

from data_utils import RAW_DIR, is_ecoli_mic, load_json


record = load_json(RAW_DIR / "peptide_107.json")
activities = record.get("targetActivities") or []

print("Peptide ID:", record.get("id"))
print("Sequence:", record.get("sequence"))
print("Number of activity measurements:", len(activities))

print("\nE. coli MIC measurements:")

for measurement in activities:
    if is_ecoli_mic(measurement):
        print(json.dumps(measurement, indent=2, ensure_ascii=False))