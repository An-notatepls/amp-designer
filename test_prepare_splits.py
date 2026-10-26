"""Offline tests for sequence grouping, partition integrity and split reuse."""

import csv
import io
import json
import random
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

import prepare_splits as split


def basic_distance(a, b):
    previous = list(range(len(b) + 1))
    for i, first in enumerate(a, 1):
        current = [i]
        for j, second in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (first != second)))
        previous = current
    return previous[-1]


def fixture_row(sequence, label, pid, amidated="False"):
    return {"sequence": sequence, "c_amidated": amidated,
            "peptide_ids_json": json.dumps([pid]), "label": str(label),
            "status": "active" if label else "inactive", "model_eligible": "True"}


def write_fixture(folder):
    rng = random.Random(123)
    sequences = set()
    while len(sequences) < 80:
        sequences.add("".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=16)))
    broad = [fixture_row(sequence, i % 2, i + 1) for i, sequence in enumerate(sorted(sequences))]
    # Include a chemical variant: it must always follow its exact sequence.
    variant = dict(broad[0], c_amidated="True", peptide_ids_json="[1001]")
    broad.append(variant)
    clinical = [broad[5], broad[6], broad[7]]
    for name, rows in zip(split.INPUT_NAMES, (broad, broad, clinical)):
        (folder / name).write_text(split.csv_text(rows, list(rows[0])), encoding="utf-8", newline="")
    return broad, clinical


class SplitTests(unittest.TestCase):
    def test_edit_distance_matches_independent_dynamic_programming(self):
        rng = random.Random(7)
        for _ in range(300):
            a = "".join(rng.choices("ACDE", k=rng.randrange(25)))
            b = "".join(rng.choices("ACDE", k=rng.randrange(25)))
            self.assertEqual(split.edit_distance(a, b), basic_distance(a, b))
            self.assertEqual(split.edit_distance(a, b), split.edit_distance(b, a))

    def test_boundary_transitive_groups_and_order_invariance(self):
        sequences = ["AAAAAAAAAA", "CAAAAAAAAA", "CCAAAAAAAA", "WWWWWWWWWW"]
        mapping, sizes, edges = split.sequence_groups(sequences, Fraction(9, 10))
        self.assertEqual(mapping[sequences[0]], mapping[sequences[2]])
        self.assertNotEqual(mapping[sequences[0]], mapping[sequences[3]])
        self.assertEqual(max(sizes.values()), 3)
        self.assertEqual(len(edges), 2)
        self.assertEqual(mapping, split.sequence_groups(reversed(sequences), Fraction(9, 10))[0])

    def test_labels_not_used_for_group_or_partition_allocation(self):
        sequences = ["AAAAAAAAAA", "CAAAAAAAAA", "WWWWWWWWWW", "GGGGGGGGGG"]
        mapping, _, _ = split.sequence_groups(sequences, Fraction(9, 10))
        assignment = split.allocate_groups(mapping, {sequences[0]}, 42, Fraction(1, 5))
        self.assertEqual(assignment[mapping[sequences[0]]], "clinical_family_holdout")
        self.assertEqual(assignment[mapping[sequences[1]]], "clinical_family_holdout")
        self.assertEqual(assignment, split.allocate_groups(mapping, {sequences[0]}, 42, Fraction(1, 5)))

    def test_cross_partition_family_is_rejected(self):
        rows = [dict(fixture_row("AAAAAAAAAA", 1, 1), group_id="one", partition="train"),
                dict(fixture_row("CAAAAAAAAA", 0, 2), group_id="one", partition="test")]
        tables = {"train": [rows[0]], "test": [rows[1]], "clinical_exploratory": [], "clinical_family_buffer": []}
        with self.assertRaisesRegex(ValueError, "crosses"):
            split.check_partitions(tables, rows, [], {r["sequence"]: "one" for r in rows}, {"one": "train"}, [])

    def test_chemical_variants_clinical_reservation_complete_exports_and_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            broad, clinical = write_fixture(root)
            original = {name: (root / name).read_bytes() for name in split.INPUT_NAMES}
            output, manifest = split.prepare(root)
            tables = {}
            for name in ("train", "test", "clinical_exploratory", "clinical_family_buffer"):
                with (output / f"{name}.csv").open(encoding="utf-8", newline="") as handle:
                    tables[name] = list(csv.DictReader(handle))
            self.assertEqual(sum(map(len, tables.values())), len(broad))
            self.assertEqual({split.identity(r) for r in tables["clinical_exploratory"]},
                             {split.identity(r) for r in clinical})
            for source_row in broad:
                variants = [r for table in tables.values() for r in table if r["sequence"] == source_row["sequence"]]
                self.assertEqual(len({r["group_id"] for r in variants}), 1)
                self.assertEqual(len({r["partition"] for r in variants}), 1)
            self.assertFalse(manifest["model_fitted_by_this_script"])
            repeated_output, repeated_manifest = split.prepare(root)
            self.assertEqual(output, repeated_output)
            self.assertEqual(manifest, repeated_manifest)
            self.assertEqual(original, {name: (root / name).read_bytes() for name in split.INPUT_NAMES})
            (output / "train.csv").write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "changed or missing"):
                split.prepare(root)

    def test_inconsistent_or_duplicate_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            broad, _ = write_fixture(root)
            clinical_path = root / split.INPUT_NAMES[2]
            changed = dict(broad[5], label=str(1 - int(broad[5]["label"])))
            changed["status"] = "active" if changed["label"] == "1" else "inactive"
            clinical_path.write_text(split.csv_text([changed], list(changed)), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "clinical|Clinical"):
                split.prepare(root)
            candidate_path = root / split.INPUT_NAMES[0]
            candidate_path.write_text(split.csv_text([broad[0], broad[0]], list(broad[0])), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                split.read_table(candidate_path, candidates=True)


if __name__ == "__main__":
    unittest.main()
