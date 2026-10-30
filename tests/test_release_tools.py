"""Release tests use synthetic fixtures and committed metadata, never new evaluation."""
import csv
import hashlib
import json
import math
import subprocess
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'modeling'))
import export_local_release as exporter
import predict_portable as portable
import reproduce_modeling as reproduce
import verify_release
from release_paths import LABEL_REL, COHORT_REL, SENSITIVITY_REL, RAW_MANIFEST_REL


def toy_model():
    return {'schema': 'amp-primary-linear-v1', 'feature_names': ['charge_proxy', 'c_amidated'],
            'source_model_sha256': portable.EXPECTED_SOURCE, 'classes': [0, 1],
            'mean': [0, 0], 'scale': [1, 1], 'coef': [1, 0], 'intercept': 0,
            'prediction_threshold': 0.5}


class ReleaseTests(unittest.TestCase):
    def test_committed_portable_model_runs_without_site_packages(self):
        result = subprocess.run([sys.executable, '-I', '-S', str(ROOT / 'modeling/predict_portable.py'),
                                 '--sequence', 'AAAAAAAAAA', '--unmodified'],
                                check=True, capture_output=True, text=True)
        prediction = json.loads(result.stdout)
        self.assertEqual(prediction['predicted_label'], 'inactive')
        self.assertEqual(prediction['source_model_sha256'], portable.EXPECTED_SOURCE)
        self.assertAlmostEqual(prediction['model_probability_active'], 0.1696300091588566, places=12)

    def test_portable_features_probability_and_tie(self):
        m = portable.validate_model(toy_model())
        row = portable.predict('AKRDEAAAAA', False, m)
        self.assertEqual(row['charge_proxy'], 0)
        self.assertEqual(row['model_probability_active'], 0.5)
        self.assertEqual(row['predicted_label'], 'inactive')
        amidated = portable.predict('AKRDEAAAAA', True, m)
        self.assertEqual(amidated['charge_proxy'], 1)
        self.assertAlmostEqual(amidated['model_probability_active'], 1/(1+math.exp(-1)))

    def test_portable_rejects_unsupported_chemistry_and_nonfinite_model(self):
        for seq, state in [('AAAA', False), ('a'*10, False), ('AXAAAAAAAA', False), ('A'*10, 'False')]:
            with self.assertRaises(ValueError): portable.predict(seq, state, toy_model())
        for change in [{'scale': [1, 0]}, {'coef': [1, float('nan')]}, {'source_model_sha256': 'wrong'}]:
            with self.assertRaises(ValueError): portable.validate_model({**toy_model(), **change})

    def test_portable_json_checksum_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as f:
            root = Path(f); path = root/'primary_model.json'
            path.write_text(json.dumps(toy_model()))
            (root/'model_manifest.json').write_text(json.dumps({'files': {path.name: exporter.sha(path)}}))
            self.assertEqual(portable.load_model(path)['classes'], [0, 1])
            path.write_text(path.read_text()+' ')
            with self.assertRaisesRegex(ValueError, 'fingerprint'): portable.load_model(path)

    def test_export_matches_original_pipeline_without_refitting(self):
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.linear_model import LogisticRegression
        model = Pipeline([('scale', StandardScaler()), ('classifier', LogisticRegression(C=.1))])
        model.fit([[-3, 0], [-2, 1], [2, 0], [3, 1]], [0, 0, 1, 1])
        with tempfile.TemporaryDirectory() as f, patch.object(exporter.original_predict, 'load_primary', return_value=(model, portable.EXPECTED_SOURCE)):
            destination = Path(f)/'model'
            with patch.object(Pipeline, 'fit', side_effect=AssertionError('Export must never fit')):
                exporter.export_model(Path(f), destination)
            m = portable.load_model(destination/'primary_model.json')
            for sequence, state in [('A'*10, False), ('AKRDEAAAAA', True), ('K'*20, False), ('D'*20, True)]:
                result = portable.predict(sequence, state, m)
                x = [[result['charge_proxy'], int(state)]]
                self.assertAlmostEqual(result['model_probability_active'], model.predict_proba(x)[0, 1], places=12)
                self.assertEqual(result['predicted_label']=='active', bool(model.predict(x)[0]))
            with self.assertRaises(FileExistsError): exporter.export_model(Path(f), destination)

    def test_numerical_comparison_rejects_drift_and_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as f:
            a,b=Path(f)/'a.csv',Path(f)/'b.csv'
            a.write_text('id,value\nx,0.75\n');b.write_text('id,value\nx,0.750000000001\n')
            reproduce.compare_rows(a,b,['id'])
            b.write_text('id,value\nx,0.80\n')
            with self.assertRaises(ValueError): reproduce.compare_rows(a,b,['id'])
            b.write_text('id,value\nx,0.75\nx,0.75\n')
            with self.assertRaises(ValueError): reproduce.compare_rows(a,b,['id'])

    def test_source_export_preserves_timestamps_and_bytes(self):
        with tempfile.TemporaryDirectory() as f:
            root=Path(f)/'source';dest=Path(f)/'export'
            def write(relative, value):
                path=root/relative;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(value));return path
            full=write(RAW_MANIFEST_REL, {'selected_ids':list(range(1,2218)), 'source_pilot_namespace':'pilot',
                'listing_started_at_utc':'2026-10-01T12:00:00+00:00', 'source_prepared_at_utc':'2026-10-01T12:01:00+00:00'})
            build=write(LABEL_REL/'build_report.json', {'policy':{'manifest_sha256':exporter.sha(full)}})
            rows=write(LABEL_REL/'ecoli_measurements.csv', {'fixture':'synthetic'})
            prov=write(COHORT_REL/'provenance_report.json', {'fixture':'synthetic'})
            write(SENSITIVITY_REL/'sensitivity_report.json', {'spec':{'build_report_sha256':exporter.sha(build),
                'measurements_sha256':exporter.sha(rows),'provenance_report_sha256':exporter.sha(prov)}})
            exporter.export_sources(root,dest)
            self.assertEqual(full.read_bytes(), (dest/'data/manifests/ecoli-full-v1/manifest.json').read_bytes())
            dates=json.loads((dest/'data/manifests/ecoli-full-v1/retrieval_metadata.json').read_text())
            self.assertIsNone(dates['record_download_dates'])
            self.assertEqual(dates['observed_metadata'][0]['value'],'2026-10-01T12:00:00+00:00')
            full.write_text('{}')
            with self.assertRaisesRegex(ValueError,'manifest'):exporter.export_sources(root,Path(f)/'bad')

    def test_committed_inputs_and_sources_verify(self):
        report=verify_release.verify()
        self.assertGreaterEqual(report['checked_fingerprints'],35)


if __name__=='__main__':
    unittest.main(verbosity=2)
