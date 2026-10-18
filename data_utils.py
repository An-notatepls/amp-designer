"""Shared I/O and exploratory filters; raw files are immutable by default."""
import json
import os
import re
import ssl
import tempfile
from pathlib import Path
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PILOT_DIR = RAW_DIR / "pilot"
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save_json(path, data, *, overwrite=False):
    """Publish complete JSON atomically; never replace raw files by default."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        if load_json(path) == data:
            return
        raise FileExistsError(f"Preserving existing file: {path}")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            # Atomic no-clobber publication on the same filesystem.
            os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def fetch_json(url):
    import certifi

    context = ssl.create_default_context(cafile=certifi.where())
    with urlopen(url, timeout=60, context=context) as response:
        return json.load(response)


def validate_record(record, peptide_id):
    if not isinstance(record, dict) or record.get("id") != peptide_id:
        raise ValueError(f"Record ID does not match requested ID {peptide_id}")
    activities = record.get("targetActivities")
    if not isinstance(activities, list) or any(
        not isinstance(item, dict) for item in activities
    ):
        raise ValueError(f"Invalid targetActivities for ID {peptide_id}")
    return record


def target_name(measurement):
    target = measurement.get("targetSpecies") or {}
    return target.get("name") or ""


def is_ecoli_mic(measurement):
    """Provisional client filter; yield is conditional on this name rule."""
    name = target_name(measurement)
    return (name == "Escherichia coli" or name.startswith("Escherichia coli ")) and (
        measurement.get("activityMeasureValue") == "MIC"
    )


def concentration_kind(value):
    if value is None or not str(value).strip():
        return "missing"
    text = str(value).strip()
    number = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    for pattern, kind in (
        (number, "single_number"),
        (rf"(?:>=|<=|>|<|≥|≤)\s*{number}", "censored"),
        (rf"{number}\s*[-–]\s*{number}", "range"),
        (rf"{number}\s*±\s*{number}", "plus_minus"),
    ):
        if re.fullmatch(pattern, text):
            return kind
    return "unrecognized"


def has_content(value):
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, dict)):
        return bool(value)
    return True


def metadata_state(record, key):
    if key not in record:
        return "absent"
    return "populated" if has_content(record[key]) else "empty"


def sequence_reasons(record):
    sequence = record.get("sequence") or ""
    if not isinstance(sequence, str):
        return ["Sequence is not text"]
    if not sequence:
        return ["Component representation" if record.get("monomers") else "Unavailable sequence"]
    reasons = []
    if not 10 <= len(sequence) <= 50:
        reasons.append("Length outside 10–50")
    if set(sequence) - STANDARD_AMINO_ACIDS:
        reasons.append("Nonstandard characters")
    return reasons


def validate_manifest(manifest):
    ids = manifest.get("selected_ids") if isinstance(manifest, dict) else None
    if not isinstance(ids, list) or not ids or any(type(i) is not int for i in ids):
        raise ValueError("Manifest must contain a nonempty list of integer selected_ids")
    if len(ids) != len(set(ids)):
        raise ValueError("Manifest contains duplicate IDs")
    return ids
