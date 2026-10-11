from collections import Counter

from data_utils import (
    RAW_DIR,
    concentration_kind,
    is_ecoli_mic,
    load_json,
    target_name,
)


STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

complexities = Counter()
sequence_status = Counter()
sequence_rejections = Counter()
nonstandard_characters = Counter()
metadata_flags = Counter()
nonempty_keys = Counter()
measurement_types = Counter()
mic_units = Counter()
coli_names = Counter()
ecoli_mic_formats = Counter()

sequence_ids = set()
ecoli_ids = set()
metadata_candidate_ids = set()

record_count = 0
missing_ids = 0
missing_mic_concentrations = 0


def has_content(value):
    """Count populated values without treating numeric zero as missing."""
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, dict)):
        return bool(value)
    return True


for path in sorted(RAW_DIR.glob("peptide_*.json")):
    record = load_json(path)

    if not isinstance(record, dict):
        raise ValueError(f"Expected a JSON object in {path.name}")

    record_count += 1

    peptide_id = record.get("id")
    if peptide_id is None:
        peptide_id = f"file:{path.name}"
        missing_ids += 1

    # Audit every top-level key, without relying on guessed field names.
    for key, value in record.items():
        if has_content(value):
            nonempty_keys[key] += 1

    complexity = record.get("complexity")
    if isinstance(complexity, dict):
        complexity_name = complexity.get("name") or "Missing"
    else:
        complexity_name = complexity or "Missing"

    complexities[complexity_name] += 1

    sequence = record.get("sequence") or ""
    chains = record.get("monomers") or []
    reasons = []

    if not sequence:
        if chains:
            sequence_status["No top-level sequence; components present"] += 1
            reasons.append("Component-chain representation")
        else:
            sequence_status["No top-level or component sequence provided"] += 1
            reasons.append("Sequence unavailable")
    else:
        sequence_status["Top-level sequence present"] += 1

        if not 10 <= len(sequence) <= 50:
            reasons.append("Length outside 10–50")

        invalid = [
            character for character in sequence
            if character not in STANDARD_AMINO_ACIDS
        ]

        if invalid:
            reasons.append("Nonstandard sequence characters")
            nonstandard_characters.update(invalid)

    sequence_rejections.update(reasons)
    sequence_passes = not reasons

    if sequence_passes:
        sequence_ids.add(peptide_id)

    # A conservative screen for recorded chemical metadata.
    # This is NOT proof of complete chemical identity.
    flags = []

    if complexity_name != "Monomer":
        flags.append("Not explicitly classified as Monomer")

    for aliases, label in (
        (("NTerminus", "nTerminus"), "N-terminal annotation"),
        (("CTerminus", "cTerminus"), "C-terminal annotation"),
    ):
        if not any(key in record for key in aliases):
            flags.append(f"Missing {label} fields")
        elif any(has_content(record.get(key)) for key in aliases):
            flags.append(f"Recorded {label}")

    for key in (
        "unusualAminoAcids",
        "opms",
        "intrachainBonds",
        "interchainBonds",
        "coordinationBonds",
    ):
        if key not in record:
            flags.append(f"Missing field: {key}")
        elif has_content(record[key]):
            flags.append(f"Recorded field: {key}")

    metadata_flags.update(flags)

    if sequence_passes and not flags:
        metadata_candidate_ids.add(peptide_id)

    for measurement in record.get("targetActivities") or []:
        measure = measurement.get("activityMeasureValue") or "Missing"
        measurement_types[measure] += 1

        name = target_name(measurement)
        if "coli" in name.casefold():
            coli_names[name] += 1

        if measure != "MIC":
            continue

        unit = measurement.get("unit") or {}
        mic_units[unit.get("name") or "Missing"] += 1

        value = measurement.get("concentration")
        kind = concentration_kind(value)

        if kind == "missing":
            missing_mic_concentrations += 1

        if is_ecoli_mic(measurement):
            ecoli_ids.add(peptide_id)
            ecoli_mic_formats[kind] += 1


def print_counter(title, counts):
    print(f"\n{title}")
    if not counts:
        print("  None")
    for label, count in counts.most_common():
        print(f"  {label}: {count}")


print("Records audited:", record_count)
print("Records with missing IDs:", missing_ids)
print("Sequence-passing IDs:", len(sequence_ids))
print("IDs with E. coli MIC measurements:", len(ecoli_ids))

print(
    "Sequence-passing IDs with E. coli MIC:",
    len(sequence_ids & ecoli_ids),
)
print(
    "Above intersection also passing conservative metadata screen:",
    len(metadata_candidate_ids & ecoli_ids),
)
print("MIC measurements without concentration:",
      missing_mic_concentrations)

print_counter("Complexity — record counts:", complexities)
print_counter("Sequence representation — record counts:", sequence_status)
print_counter("Sequence rejection reasons — record counts:",
              sequence_rejections)
print_counter("Nonstandard characters — occurrence counts:",
              nonstandard_characters)
print_counter("Chemical metadata flags — record counts:", metadata_flags)
print_counter("Populated top-level keys — record counts:", nonempty_keys)
print_counter("Measurement types — measurement counts:", measurement_types)
print_counter("MIC units — measurement counts:", mic_units)
print_counter("Target names containing 'coli' — measurement counts:",
              coli_names)
print_counter("E. coli MIC concentration formats — measurement counts:",
              ecoli_mic_formats)