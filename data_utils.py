import json
import ssl
from pathlib import Path
from urllib.request import urlopen

import certifi


PROJECT_ROOT = Path(__file__).resolve().parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"


def load_json(path):
    """Read a saved JSON file."""
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path, data):
    """Save JSON in a readable format."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def fetch_json(url):
    """Download JSON with certificate verification enabled."""
    context = ssl.create_default_context(cafile=certifi.where())

    with urlopen(url, timeout=30, context=context) as response:
        return json.load(response)


def target_name(measurement):
    """Return the original target name without changing it."""
    target = measurement.get("targetSpecies") or {}
    return target.get("name") or ""


def is_ecoli_mic(measurement):
    """Provisional name-based filter; review against target metadata."""
    name = target_name(measurement)
    is_ecoli = (
        name == "Escherichia coli"
        or name.startswith("Escherichia coli ")
    )
    return is_ecoli and measurement.get("activityMeasureValue") == "MIC"