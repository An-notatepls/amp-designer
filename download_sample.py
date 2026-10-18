import json
import time
from urllib.error import URLError
from urllib.request import urlopen
from data_utils import RAW_DIR, fetch_json, save_json




def main():
    # Retrieve real peptide IDs from the database.
    listing = fetch_json(
        "https://dbaasp.org/peptides?limit=10&offset=0"
    )
    save_json(RAW_DIR / "sample_listing.json", listing)

    peptides = listing["data"]
    available = 0

    for peptide in peptides:
        peptide_id = peptide["id"]
        output_path = RAW_DIR / f"peptide_{peptide_id}.json"

        # Reuse previously downloaded records.
        if output_path.exists():
            print(f"Already saved: peptide {peptide_id}")
            available += 1
            continue

        try:
            record = fetch_json(
                f"https://dbaasp.org/peptides/{peptide_id}"
            )
            save_json(output_path, record)
            print(f"Downloaded: peptide {peptide_id}")
            available += 1

        except (URLError, TimeoutError, json.JSONDecodeError) as error:
            print(f"Failed: peptide {peptide_id}: {error}")

        time.sleep(0.5)

    print(f"\nSample records available: {available}/{len(peptides)}")


if __name__ == "__main__":
    main()