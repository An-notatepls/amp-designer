"""Resume a manifest-defined pilot; validate cached records before skipping."""
import argparse
import time
from pathlib import Path

from data_utils import fetch_json, load_json, save_json, validate_manifest, validate_record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    ids = validate_manifest(load_json(args.manifest))
    records_dir = args.manifest.parent / "records"
    failed, downloaded, cached = [], 0, 0
    try:
        for index, peptide_id in enumerate(ids, 1):
            out = records_dir / f"peptide_{peptide_id}.json"
            if out.exists():
                try:
                    validate_record(load_json(out), peptide_id)
                    cached += 1
                    print(f"{index}/{len(ids)} validated cache: {peptide_id}", flush=True)
                except (OSError, ValueError) as error:
                    failed.append(peptide_id)
                    print(f"Invalid cache preserved: {out}: {error}. Use a new run or inspect/quarantine it manually.", flush=True)
                continue
            for attempt in range(3):
                try:
                    record = validate_record(fetch_json(f"https://dbaasp.org/peptides/{peptide_id}"), peptide_id)
                    save_json(out, record)
                    downloaded += 1
                    print(f"{index}/{len(ids)} downloaded: {peptide_id}", flush=True)
                    break
                except (OSError, ValueError) as error:
                    print(f"ID {peptide_id}, attempt {attempt + 1}: {error}", flush=True)
                    if attempt < 2:
                        time.sleep(2 ** attempt)
            else:
                failed.append(peptide_id)
            time.sleep(0.5)
    except KeyboardInterrupt:
        print("\nInterrupted. Rerun the same command to resume; complete records are retained.")
        return
    print(f"Downloaded now: {downloaded}; validated cache: {cached}; failed: {failed}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
