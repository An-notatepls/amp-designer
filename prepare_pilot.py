"""Create a random pilot manifest in a query- and run-specific namespace."""
import argparse
import hashlib
import json
import random
import time
from datetime import datetime, timezone
from urllib.parse import urlencode

from data_utils import RAW_DIR, PILOT_DIR, fetch_json, load_json, save_json, validate_manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", default="ecoli-v1")
    parser.add_argument("--size", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--page-size", type=int, default=100)
    args = parser.parse_args()
    if not args.run or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in args.run):
        parser.error("--run must use letters, numbers, hyphens, or underscores")
    if not 1 <= args.page_size <= 100 or args.size < 1:
        parser.error("Page size must be 1–100 and sample size must be positive")

    spec = {"endpoint": "https://dbaasp.org/peptides",
            "filters": {"targetSpecies.value": "Escherichia coli"},
            "page_size": args.page_size, "sample_size": args.size, "seed": args.seed}
    digest = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]
    namespace = f"{args.run}-{digest}"
    listings = RAW_DIR / "listings" / namespace
    pilot = PILOT_DIR / namespace
    manifest_path = pilot / "manifest.json"
    if manifest_path.exists():
        manifest = load_json(manifest_path)
        validate_manifest(manifest)
        if manifest.get("spec") != spec:
            raise ValueError("Existing manifest configuration does not match")
        print(f"Using frozen manifest: {manifest_path}")
        return

    started = datetime.now(timezone.utc).isoformat()
    ids, offset, total = [], 0, None
    while True:
        query = {**spec["filters"], "limit": args.page_size, "offset": offset}
        url = spec["endpoint"] + "?" + urlencode(query)
        path = listings / f"page_n{args.page_size}_offset{offset}.json"
        if path.exists():
            envelope = load_json(path)
            if envelope.get("query") != query or envelope.get("url") != url:
                raise ValueError("Cached page query does not match; use a new --run")
        else:
            envelope = {"query": query, "url": url,
                        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                        "response": fetch_json(url)}
            save_json(path, envelope)
            time.sleep(0.5)
        page = envelope.get("response")
        if not isinstance(page, dict) or not isinstance(page.get("data"), list):
            raise ValueError(f"Invalid listing page: {path}")
        count = page.get("totalCount")
        if type(count) is not int or count < 0:
            raise ValueError("Invalid totalCount")
        if total is None:
            total = count
        if count != total:
            raise ValueError("The database may have changed. Preserve this run and rerun with a new --run name (or deliberately clear its listing pages).")
        entries = page["data"]
        if len(entries) > args.page_size or (not entries and offset < total):
            raise ValueError("Unexpected pagination; use a new --run and inspect the API")
        for entry in entries:
            peptide_id = entry.get("id") if isinstance(entry, dict) else None
            if type(peptide_id) is not int:
                raise ValueError("Listing entry lacks an integer ID")
            ids.append(peptide_id)
        offset += len(entries)
        print(f"Listing entries: {offset}/{total}", flush=True)
        if offset >= total:
            break
    population = sorted(set(ids))
    if not population or len(ids) != total or len(population) != total:
        raise ValueError("Unique ID count differs from totalCount. The database may have changed; use a new --run (or deliberately clear this run's listing pages) and rerun.")
    manifest = {"spec": spec, "run": args.run, "namespace": namespace,
                "listing_started_at_utc": started,
                "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
                "reported_total": total, "population_ids": population,
                "selected_ids": random.Random(args.seed).sample(population, min(args.size, total))}
    save_json(manifest_path, manifest)
    print(f"Selected IDs: {len(manifest['selected_ids'])}")
    print(f"Manifest: {manifest_path}")
    print("Pass this manifest path to download_pilot.py and audit_pilot.py.")


if __name__ == "__main__":
    main()
