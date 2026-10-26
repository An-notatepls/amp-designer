"""Prepare immutable sequence-group partitions; never fit a model.

Usage: python prepare_splits.py PATH/TO/PROVENANCE_OUTPUT

Requires broad_candidates.csv, broad_audit.csv, clinical_only_candidates.csv.
Defaults: 80% normalized full-length edit similarity, seed 42, 20% group-hash
probability for the primary test set. Exact test size is not guaranteed.
Clinical-only identities and their entire connected sequence groups are
reserved before primary train/test assignment. Other members of those groups
are buffers, not extra independent clinical test examples.

Similarity = 1 - Levenshtein_distance / max(sequence_lengths). Each substitution,
insertion and deletion costs one. This is not standard biological alignment
identity or proof of homology; it does not protect every shared motif. Families
are connected components and may contain distant members linked by a chain.
Group on all audit sequences, without activity labels or model results.
No random seed or cutoff search is performed. No third-party packages required.
"""

import argparse
import csv
import hashlib
import io
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
INPUT_NAMES = ("broad_candidates.csv", "broad_audit.csv", "clinical_only_candidates.csv")
VERSION = "grouped-v1"


def sha_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def identity(row):
    return row["sequence"], row["c_amidated"]


def read_table(path, candidates):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        required = {"sequence", "c_amidated", "label", "status", "model_eligible", "peptide_ids_json"}
        if not required <= set(header) or len(header) != len(set(header)):
            raise ValueError(f"Unexpected header: {path}")
        rows = list(reader)
    seen = set()
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed CSV row in {path}")
        sequence = row["sequence"]
        if not 10 <= len(sequence) <= 50 or not set(sequence) <= STANDARD:
            raise ValueError(f"Noncanonical/out-of-scope sequence in {path}")
        if row["c_amidated"] not in {"True", "False"} or row["model_eligible"] not in {"True", "False"}:
            raise ValueError(f"Invalid boolean in {path}")
        key = identity(row)
        if key in seen:
            raise ValueError(f"Duplicate sequence+amidation identity in {path}")
        seen.add(key)
        ids = json.loads(row["peptide_ids_json"])
        if (not isinstance(ids, list) or not ids or any(type(pid) is not int or pid < 1 for pid in ids)
                or len(ids) != len(set(ids))):
            raise ValueError(f"Invalid peptide IDs in {path}")
        if candidates or row["model_eligible"] == "True":
            if row["model_eligible"] != "True" or row["label"] not in {"0", "1"}:
                raise ValueError(f"Candidate row is not definitely labeled in {path}")
            if row["status"] != ("active" if row["label"] == "1" else "inactive"):
                raise ValueError(f"Status disagrees with label in {path}")
    return sorted(rows, key=identity)


def validate_inputs(broad, audit, clinical):
    if not broad or not audit or not clinical:
        raise ValueError("This split design requires nonempty broad, audit and clinical-only tables")
    audit_by_id = {identity(row): row for row in audit}
    broad_by_id = {identity(row): row for row in broad}
    ready = {identity(row) for row in audit if row["model_eligible"] == "True"}
    if set(broad_by_id) != ready:
        raise ValueError("Broad candidates do not equal the eligible identities in the audit table")
    for key, row in broad_by_id.items():
        if row["label"] != audit_by_id[key]["label"]:
            raise ValueError("Broad candidate label disagrees with audit label")
    for row in clinical:
        key = identity(row)
        if key not in broad_by_id or row["label"] != broad_by_id[key]["label"]:
            raise ValueError("Clinical-only candidate is absent from broad candidates or its label differs")


def edit_distance(a, b):
    """Exact unit-cost Levenshtein distance using a bit-vector recurrence."""
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return len(b)
    width = len(a)
    mask, top = (1 << width) - 1, 1 << (width - 1)
    equality = {}
    for i, residue in enumerate(a):
        equality[residue] = equality.get(residue, 0) | (1 << i)
    positive, negative, score = mask, 0, width
    for residue in b:
        equal = equality.get(residue, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        plus = negative | ~(horizontal | positive)
        minus = positive & horizontal
        score += bool(plus & top) - bool(minus & top)
        plus, minus = (plus << 1) | 1, minus << 1
        positive = (minus | ~(vertical | plus)) & mask
        negative = (plus & vertical) & mask
    return score


def sequence_groups(sequences, similarity, progress=False):
    sequences = sorted(set(sequences))
    parent = list(range(len(sequences)))

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    denominator = similarity.denominator
    difference = denominator - similarity.numerator
    edges = []
    for i, first in enumerate(sequences):
        for j in range(i + 1, len(sequences)):
            second = sequences[j]
            length = max(len(first), len(second))
            if abs(len(first) - len(second)) * denominator > difference * length:
                continue
            distance = edit_distance(first, second)
            if distance * denominator <= difference * length:
                parent[root(i)] = root(j)
                edges.append((first, second, distance))
        if progress and (i + 1) % 100 == 0:
            print(f"Grouped comparisons for {i + 1}/{len(sequences)} sequences", flush=True)
    components = defaultdict(list)
    for index, sequence in enumerate(sequences):
        components[root(index)].append(sequence)
    mapping, sizes = {}, {}
    for members in components.values():
        group = "family_" + hashlib.sha256("\n".join(members).encode()).hexdigest()[:16]
        if group in sizes:
            raise ValueError("Group-hash collision")
        sizes[group] = len(members)
        mapping.update({sequence: group for sequence in members})
    return mapping, sizes, edges


def allocate_groups(mapping, clinical_sequences, seed, test_probability):
    reserved = {mapping[sequence] for sequence in clinical_sequences}
    assignment = {}
    for group in sorted(set(mapping.values())):
        number = int(hashlib.sha256(f"{seed}:{group}".encode()).hexdigest(), 16)
        selected = number * test_probability.denominator < test_probability.numerator * (1 << 256)
        assignment[group] = ("clinical_family_holdout" if group in reserved else "test" if selected else "train")
    return assignment


def partition_rows(broad, clinical, mapping, assignment):
    clinical_ids = {identity(row) for row in clinical}
    tables = {"train": [], "test": [], "clinical_exploratory": [], "clinical_family_buffer": []}
    for row in broad:
        group = mapping[row["sequence"]]
        part = assignment[group]
        role = part
        if part == "clinical_family_holdout":
            role = "clinical_exploratory" if identity(row) in clinical_ids else "clinical_family_buffer"
        tables[role].append({**row, "group_id": group, "partition": role})
    return tables


def check_partitions(tables, broad, clinical, mapping, assignment, edges):
    records = [row for table in tables.values() for row in table]
    if len(records) != len(broad) or {identity(r) for r in records} != {identity(r) for r in broad}:
        raise ValueError("Partitions do not cover broad candidates exactly once")
    if len({identity(row) for row in records}) != len(records):
        raise ValueError("Identity appears in more than one partition")
    if {identity(r) for r in tables["clinical_exploratory"]} != {identity(r) for r in clinical}:
        raise ValueError("Clinical-only reservation is incomplete")
    locations = defaultdict(set)
    for row in records:
        # Clinical examples and their buffers share the same held-out family.
        locations[row["group_id"]].add("clinical_family_holdout" if row["partition"].startswith("clinical_") else row["partition"])
    if any(len(parts) != 1 for parts in locations.values()):
        raise ValueError("Sequence group crosses primary train/test/clinical boundaries")
    for first, second, _ in edges:
        if mapping[first] != mapping[second] or assignment[mapping[first]] != assignment[mapping[second]]:
            raise ValueError("A close sequence pair crosses partition boundaries")
    for name in ("train", "test"):
        if {r["label"] for r in tables[name]} != {"0", "1"}:
            raise ValueError(f"Fixed split lacks both classes in {name}; inspect the declared design rather than searching seeds for model performance")


def counts(rows):
    return {
        "examples": len(rows), "unique_sequences": len({r["sequence"] for r in rows}),
        "groups": len({r["group_id"] for r in rows}),
        "active": sum(r["label"] == "1" for r in rows),
        "inactive": sum(r["label"] == "0" for r in rows),
    }


def csv_text(rows, header):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=header)
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def existing_result(output, spec):
    manifest_path = output / "split_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"Incomplete existing split directory preserved: {output}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("spec") != spec:
        raise ValueError("Existing split specification differs; refusing to overwrite")
    for name, digest in manifest["output_sha256"].items():
        path = output / name
        if not path.is_file() or sha_file(path) != digest:
            raise ValueError(f"Existing split file changed or missing: {path}")
    return manifest


def prepare(source, similarity=Fraction(4, 5), seed=42, test_probability=Fraction(1, 5), progress=False):
    source = Path(source).expanduser().resolve()
    paths = {name: source / name for name in INPUT_NAMES}
    broad = read_table(paths[INPUT_NAMES[0]], candidates=True)
    audit = read_table(paths[INPUT_NAMES[1]], candidates=False)
    clinical = read_table(paths[INPUT_NAMES[2]], candidates=True)
    validate_inputs(broad, audit, clinical)
    spec = {
        "version": VERSION, "source_sha256": {name: sha_file(path) for name, path in paths.items()},
        "script_sha256": sha_file(Path(__file__)), "seed": seed,
        "similarity_fraction": str(similarity), "test_group_probability": str(test_probability),
        "similarity_method": "1 - unit-cost global Levenshtein distance / longer length",
        "group_method": "connected components over all unique audit sequences; amidation ignored for grouping",
        "allocation_method": "reserve clinical-only families; otherwise SHA256(str(seed)+':'+group_id) threshold; independent of labels",
        "clinical_policy": "all clinical-only identities held out, all related broad examples buffered; exploratory only",
        "limitations": "Not proof of homology; local motifs may remain shared; family chaining is possible; conditional on curated eligibility",
    }
    digest = hashlib.sha256(canonical_json(spec).encode()).hexdigest()[:16]
    output = source / "splits" / f"{VERSION}-{digest}"
    if output.exists():
        return output, existing_result(output, spec)
    mapping, sizes, edges = sequence_groups([row["sequence"] for row in audit], similarity, progress)
    assignment = allocate_groups(mapping, {r["sequence"] for r in clinical}, seed, test_probability)
    tables = partition_rows(broad, clinical, mapping, assignment)
    check_partitions(tables, broad, clinical, mapping, assignment, edges)
    header = list(broad[0]) + ["group_id", "partition"]
    payloads = {f"{name}.csv": csv_text(table, header) for name, table in tables.items()}
    index = [{"sequence": sequence, "group_id": group, "group_sequences_in_audit": sizes[group],
              "group_partition": assignment[group]} for sequence, group in sorted(mapping.items())]
    payloads["sequence_groups.csv"] = csv_text(index, list(index[0]))
    assignment_rows = [{k: r[k] for k in ("sequence", "c_amidated", "peptide_ids_json", "group_id", "partition")}
                       for r in [item for table in tables.values() for item in table]]
    assignment_rows.sort(key=lambda r: (r["sequence"], r["c_amidated"]))
    payloads["split_assignments.csv"] = csv_text(assignment_rows, list(assignment_rows[0]))
    manifest = {
        "spec": spec, "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "model_fitted_by_this_script": False,
        "counts": {name: counts(table) for name, table in tables.items()},
        "audit_unique_sequences": len(mapping), "audit_groups": len(sizes),
        "qualifying_close_pairs": len(edges),
        "largest_audit_group_sequences": max(sizes.values()),
        "checks": {"identity_partition_unique": True,
                   "no_family_overlap_across_train_test_and_clinical_reservation": True,
                   "all_qualifying_pairs_protected": True, "clinical_families_absent_from_training": True,
                   "primary_train_and_test_have_both_classes": True},
        "evaluation_rule": "Use train.csv for all fitting/tuning with group-aware CV. Use test.csv only after selecting the model. Clinical metrics are exploratory. Do not train on clinical_family_buffer.csv.",
        "buffer_note": "Clinical exploratory rows and their buffers may share families, or different amidation states of the same sequence; the entire reservation is disjoint from train and primary test.",
        "output_sha256": {name: hashlib.sha256(content.encode()).hexdigest() for name, content in payloads.items()},
        "source_attribution": "Derived from DBAASP data, https://dbaasp.org; see project README and DBAASP terms/citation.",
    }
    payloads["split_manifest.json"] = json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".split-preparing-", dir=output.parent))
    try:
        for name, content in payloads.items():
            with (temporary / name).open("w", encoding="utf-8", newline="") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
        # Rename a complete directory; never replace a populated split.
        try:
            temporary.rename(output)
        except OSError:
            if output.exists():
                return output, existing_result(output, spec)
            raise
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return output, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("provenance_dir", type=Path)
    parser.add_argument("--similarity", default="0.8")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--test-probability", default="0.2")
    args = parser.parse_args()
    try:
        similarity, probability = Fraction(args.similarity), Fraction(args.test_probability)
    except (ValueError, ZeroDivisionError):
        parser.error("Similarity and test probability must be finite numbers")
    if not 0 < similarity <= 1 or not 0 < probability < 1:
        parser.error("Require 0 < similarity <= 1 and 0 < test probability < 1")
    output, manifest = prepare(args.provenance_dir, similarity, args.seed, probability, progress=True)
    print(json.dumps(manifest["counts"], indent=2))
    print("Checks:", json.dumps(manifest["checks"], sort_keys=True))
    print("Split files:", output)
    print("No model was trained. Use train.csv for fitting and group-aware model selection.")
    print("Preserve split_manifest.json, sequence_groups.csv and split_assignments.csv for reproducibility.")
    print("The primary test and clinical holdout are for final evaluation; do not search seeds after inspecting model performance.")


if __name__ == "__main__":
    main()
