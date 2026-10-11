import json
from pathlib import Path
from urllib.request import urlopen

import ssl
import certifi

#Locate the project folder using this script's location
PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_PATH = PROJECT_ROOT / "data" / "raw" / "peptide_107.json"

# Download one record from DBAASP.
url = "https://dbaasp.org/peptides/107"

ssl_context = ssl.create_default_context(cafile=certifi.where())

with urlopen(url, timeout=30, context=ssl_context) as response:
    record = json.load(response)


# Save the original record for inspection.
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH.write_text(
    json.dumps(record, indent=2, ensure_ascii=False),
    encoding="utf-8",
)

print("Saved record to:", OUTPUT_PATH)
print("Peptide ID:", record.get("id"))
print("Sequence:", record.get("sequence"))

activities = record.get("targetActivities", [])
print("Number of activity measurements:", len(activities))

# Inspect the first measurement without assuming its structure.
if activities:
    print("\nFirst activity measurement:")
    print(json.dumps(activities[0], indent=2, ensure_ascii=False))


print("\nAll activity measurements:")

for index, measurement in enumerate(activities, start=1):
    target = measurement.get("targetSpecies") or {}
    unit = measurement.get("unit") or {}

    print(
        f"{index}. "
        f"Target: {target.get('name')} | "
        f"Measure: {measurement.get('activityMeasureValue')} | "
        f"Concentration: {measurement.get('concentration')} | "
        f"Unit: {unit.get('name')}"
    )


# First measurement 

print("\nE. coli MIC measurements:")

for measurement in activities:
    target = measurement.get("targetSpecies") or {}
    target_name = target.get("name", "")
    measure = measurement.get("activityMeasureValue")

    is_ecoli = (
        target_name == "Escherichia coli"
        or target_name.startswith("Escherichia coli ")
    )
    is_mic = measure == "MIC"

    if is_ecoli and is_mic:
        unit = measurement.get("unit") or {}

        print(
            f"{target_name} | "
            f"{measurement.get('concentration')} | "
            f"{unit.get('name')}"
        )