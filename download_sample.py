import json
import ssl
import time
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

import certifi


PROJECT_ROOT = Path(__file__).resolve().parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())


def fetch_json(url):
    """Download JSON using trusted HTTPS certificates."""
    with urlopen(url, timeout=30, context=SSL_CONTEXT) as response:
        return json.load(response)


def save_json(path, data):
    """Save JSON in a readable format."""
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


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