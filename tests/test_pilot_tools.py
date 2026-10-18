"""Offline checks for retrieval integrity and pilot population accounting."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import audit_pilot
import data_utils
import download_pilot
import prepare_pilot


def record(peptide_id, activities=None):
    return {"id": peptide_id, "sequence": "ACDEFGHIKL",
            "targetActivities": activities or []}


class PilotToolsTests(unittest.TestCase):
    def test_atomic_write_failure_and_immutability(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "record.json"
            def interrupted_dump(data, handle, **kwargs):
                handle.write('{"partial":')
                raise KeyboardInterrupt()
            with patch.object(data_utils.json, "dump", side_effect=interrupted_dump):
                with self.assertRaises(KeyboardInterrupt):
                    data_utils.save_json(path, {"id": 1})
            self.assertFalse(path.exists())
            self.assertEqual(list(Path(folder).iterdir()), [])
            data_utils.save_json(path, {"id": 1})
            data_utils.save_json(path, {"id": 1})
            with self.assertRaises(FileExistsError):
                data_utils.save_json(path, {"id": 2})
            self.assertEqual(data_utils.load_json(path), {"id": 1})

    def test_audit_ignores_unselected_files_and_reports_incomplete(self):
        activity = {"targetSpecies": {"name": "Escherichia coli ATCC 25922"},
                    "activityMeasureValue": "MIC", "concentration": ">100"}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            data_utils.save_json(path / "peptide_1.json", record(1, [activity]))
            data_utils.save_json(path / "peptide_99.json", record(99, [activity]))
            report = audit_pilot.summarize({"selected_ids": [1, 2]}, path)
            self.assertEqual(report["valid_downloaded_records"], 1)
            self.assertEqual(report["missing_ids"], [2])
            self.assertIsNone(report["filter_yield_of_selected_pilot"])
            self.assertEqual(report["sequence_and_ecoli_mic_records"], 1)
            data_utils.save_json(path / "peptide_2.json", record(2))
            report = audit_pilot.summarize({"selected_ids": [1, 2]}, path)
            self.assertEqual(report["filter_yield_of_selected_pilot"], 0.5)
            self.assertEqual(report["counters"]["ecoli_mic_formats"], {"censored": 1})

    def test_metadata_states_preserve_unknowns_and_zero(self):
        self.assertEqual(data_utils.metadata_state({}, "x"), "absent")
        for value in [None, "", " ", [], {}]:
            self.assertEqual(data_utils.metadata_state({"x": value}, "x"), "empty")
        self.assertEqual(data_utils.metadata_state({"x": 0}, "x"), "populated")
        self.assertEqual(data_utils.metadata_state({"x": {"name": "AMD"}}, "x"), "populated")

    def test_downloader_validates_cache_without_overwriting_it(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            manifest = root / "manifest.json"
            data_utils.save_json(manifest, {"selected_ids": [1, 2]})
            cached = root / "records" / "peptide_1.json"
            cached.parent.mkdir()
            cached.write_text('{"broken":', encoding="utf-8")
            with patch("sys.argv", ["download_pilot.py", str(manifest)]), \
                 patch.object(download_pilot, "fetch_json", return_value=record(2)) as fetch, \
                 patch.object(download_pilot.time, "sleep"), \
                 contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit):
                    download_pilot.main()
            self.assertEqual(cached.read_text(), '{"broken":')
            self.assertEqual(fetch.call_count, 1)
            self.assertEqual(data_utils.load_json(root / "records" / "peptide_2.json")["id"], 2)

    def test_query_namespace_and_frozen_manifest(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with patch.object(prepare_pilot, "RAW_DIR", root), \
                 patch.object(prepare_pilot, "PILOT_DIR", root / "pilot"), \
                 patch.object(prepare_pilot, "fetch_json", return_value={"totalCount": 3, "data": [{"id": 1}, {"id": 2}, {"id": 3}]}) as fetch, \
                 patch.object(prepare_pilot.time, "sleep"), \
                 contextlib.redirect_stdout(io.StringIO()):
                with patch("sys.argv", ["prepare_pilot.py", "--size", "2"]):
                    prepare_pilot.main()
                    prepare_pilot.main()
                self.assertEqual(fetch.call_count, 1)
                with patch("sys.argv", ["prepare_pilot.py", "--size", "3"]):
                    prepare_pilot.main()
                self.assertEqual(fetch.call_count, 2)
            manifests = list((root / "pilot").glob("*/manifest.json"))
            self.assertEqual(len(manifests), 2)

    def test_database_change_stops_listing(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pages = [{"totalCount": 2, "data": [{"id": 1}]},
                     {"totalCount": 3, "data": [{"id": 2}]}]
            with patch.object(prepare_pilot, "RAW_DIR", root), \
                 patch.object(prepare_pilot, "PILOT_DIR", root / "pilot"), \
                 patch.object(prepare_pilot, "fetch_json", side_effect=pages), \
                 patch.object(prepare_pilot.time, "sleep"), \
                 patch("sys.argv", ["prepare_pilot.py", "--page-size", "1"]), \
                 contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, "database may have changed"):
                    prepare_pilot.main()
            self.assertFalse(list((root / "pilot").glob("*/manifest.json")))


if __name__ == "__main__":
    unittest.main()
