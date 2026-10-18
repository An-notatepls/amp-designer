"""Expand the pilot's frozen listing cohort and reuse its validated records.

No network requests. The result includes every ID in the saved filtered listing,
not every DBAASP entry and not necessarily only E. coli MIC records. Listing
membership is frozen; individual peptide records may have different retrieval
dates. Raw pilot records and the original manifest are never overwritten.
"""

import argparse
import hashlib
import shlex
from pathlib import Path

from data_utils import RAW_DIR, load_json, save_json, validate_manifest, validate_record


def prepare(pilot_manifest_path):
    pilot_manifest_path = Path(pilot_manifest_path).expanduser().resolve()
    pilot = load_json(pilot_manifest_path)
    selected = validate_manifest(pilot)
    population = pilot.get("population_ids")
    if (not isinstance(population, list) or not population
            or any(type(pid) is not int or pid < 1 for pid in population)
            or len(set(population)) != len(population)):
        raise ValueError("Pilot manifest must contain a nonempty, unique integer population_ids list")
    if pilot.get("reported_total") != len(population):
        raise ValueError("population_ids count differs from reported_total")
    if not set(selected) <= set(population):
        raise ValueError("Pilot selected IDs are outside its listing population")
    spec = pilot.get("spec")
    if not isinstance(spec, dict) or not isinstance(spec.get("filters"), dict):
        raise ValueError("Pilot manifest lacks its original query specification")

    # Validate all source records before creating a destination.
    cached = []
    for pid in selected:
        source = pilot_manifest_path.parent / "records" / f"peptide_{pid}.json"
        record = validate_record(load_json(source), pid)
        cached.append((pid, record))

    digest = hashlib.sha256(pilot_manifest_path.read_bytes()).hexdigest()
    namespace = f"ecoli-full-v1-{digest[:16]}"
    manifest_path = RAW_DIR / "full" / namespace / "manifest.json"
    full_spec = dict(spec)
    full_spec.update(sample_size=len(population), seed=None,
                     selection_method="all IDs in the source listing cohort")
    full = {
        "spec": full_spec, "run": "ecoli-full-v1", "namespace": namespace,
        "listing_started_at_utc": pilot.get("listing_started_at_utc"),
        "source_prepared_at_utc": pilot.get("prepared_at_utc"),
        "source_manifest_sha256": digest,
        "source_pilot_namespace": pilot.get("namespace"),
        "reported_total": len(population), "population_ids": sorted(population),
        "selected_ids": sorted(population),
        "cached_pilot_ids": sorted(selected),
        "scope": "All IDs returned by the saved server-filtered listing; client MIC and chemistry filters still required",
        "retrieval_caveat": "Records are reused or fetched at different times; this is not a simultaneous database snapshot",
    }
    validate_manifest(full)
    save_json(manifest_path, full)
    for pid, record in cached:
        save_json(manifest_path.parent / "records" / f"peptide_{pid}.json", record)
    return manifest_path, len(population), len(cached)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot_manifest", type=Path)
    args = parser.parse_args()
    path, total, cached = prepare(args.pilot_manifest)
    print(f"Listing cohort IDs: {total}")
    print(f"Pilot records reused: {cached}")
    print(f"IDs beyond the pilot: {total - cached}")
    print(f"Full manifest: {path}")
    print("Next, run this exact command (rerun it to resume):")
    print(shlex.join(["python", "download_pilot.py", str(path)]))


if __name__ == "__main__":
    main()
