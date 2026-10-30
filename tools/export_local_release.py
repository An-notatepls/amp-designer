"""Export verified original local artifacts; never refit or rewrite frozen files.

Run on the original Mac with .venv. Output is a new ZIP suitable for review and
committing. Never load joblib artifacts from an untrusted source.
"""
import argparse
import hashlib
import json
import shutil
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from release_paths import (ROOT, LABEL_REL, COHORT_REL, PRETEST_REL,
                           SENSITIVITY_REL, RAW_MANIFEST_REL)
sys.path.insert(0, str(ROOT / 'modeling'))
import predict as original_predict
import predict_portable as portable


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export_model(pretest_dir, destination):
    import numpy as np
    pipeline, digest = original_predict.load_primary(pretest_dir)
    if digest != portable.EXPECTED_SOURCE:
        raise ValueError('This is not the original evaluated primary model')
    scale, classifier = pipeline.named_steps['scale'], pipeline.named_steps['classifier']
    model = {'schema': 'amp-primary-linear-v1', 'feature_names': ['charge_proxy', 'c_amidated'],
             'classes': [int(v) for v in pipeline.classes_], 'mean': scale.mean_.tolist(),
             'scale': scale.scale_.tolist(), 'coef': classifier.coef_[0].tolist(),
             'intercept': float(classifier.intercept_[0]), 'prediction_threshold': 0.5,
             'source_model_sha256': digest,
             'note': 'Exact parameters exported from evaluated frozen model; uncalibrated database prediction.'}
    portable.validate_model(model)
    # Cover the supported charge/amidation feature range; no outcomes used.
    x = np.array([[charge, state] for charge in range(-50, 52) for state in (0, 1)], dtype=float)
    decisions = model['intercept'] + ((x - scale.mean_) / scale.scale_) @ np.array(model['coef'])
    original_decisions = pipeline.decision_function(x)
    np.testing.assert_allclose(decisions, original_decisions, rtol=1e-12, atol=1e-12)
    if not np.array_equal(decisions > 0, pipeline.predict(x) == 1):
        raise ValueError('Portable classification differs on the feature grid')
    destination.mkdir(parents=True)
    path = destination / 'primary_model.json'
    path.write_text(json.dumps(model, indent=2, allow_nan=False) + '\n')
    (destination / 'model_manifest.json').write_text(json.dumps({
        'source_model_sha256': digest, 'files': {path.name: sha(path)},
        'export_validation': 'A 204-point charge/amidation grid; decision tolerance 1e-12; labels identical.'}, indent=2) + '\n')


def export_sources(root, destination):
    root, destination = Path(root), Path(destination)
    sensitivity = json.loads((root / SENSITIVITY_REL / 'sensitivity_report.json').read_text())['spec']
    build_path = root / LABEL_REL / 'build_report.json'
    if sha(build_path) != sensitivity['build_report_sha256']:
        raise ValueError('Local build report differs from frozen sensitivity report')
    build = json.loads(build_path.read_text())
    full_path = root / RAW_MANIFEST_REL
    if sha(full_path) != build['policy']['manifest_sha256']:
        raise ValueError('Full cohort manifest differs from the original build')
    full = json.loads(full_path.read_text())
    ids = full.get('selected_ids')
    if not isinstance(ids, list) or len(ids) != 2217 or len(set(ids)) != 2217:
        raise ValueError('Unexpected full cohort IDs')
    local_sources = [(full_path, Path('data/manifests/ecoli-full-v1/manifest.json'), build['policy']['manifest_sha256']),
                     (build_path, LABEL_REL / 'build_report.json', sensitivity['build_report_sha256']),
                     (root / LABEL_REL / 'ecoli_measurements.csv', LABEL_REL / 'ecoli_measurements.csv', sensitivity['measurements_sha256']),
                     (root / COHORT_REL / 'provenance_report.json', COHORT_REL / 'provenance_report.json', sensitivity['provenance_report_sha256'])]
    files = {}
    for source, relative, expected in local_sources:
        if sha(source) != expected:
            raise ValueError(f'Original source fingerprint differs: {source}')
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        files[str(relative)] = expected
    observed = []
    keys = ('listing_started_at_utc', 'source_prepared_at_utc', 'prepared_at_utc')
    for key in keys:
        if full.get(key):
            observed.append({'source': 'data/manifests/ecoli-full-v1/manifest.json', 'field': key, 'value': full[key]})
    listings = root / 'data/raw/listings' / str(full.get('source_pilot_namespace', ''))
    if listings.is_dir():
        for path in sorted(listings.glob('*.json')):
            envelope = json.loads(path.read_text())
            if envelope.get('retrieved_at_utc'):
                observed.append({'source': str(path.relative_to(root)), 'source_sha256': sha(path),
                                 'field': 'retrieved_at_utc', 'value': envelope['retrieved_at_utc']})
    dates = {'documented_at_utc': datetime.now(timezone.utc).isoformat(),
             'observed_metadata': observed,
             'record_download_dates': None,
             'caveat': 'Listing/manifest timestamps are preserved as labeled. Per-record retrieval timestamps were not saved by the downloader. File modification times and documentation dates are not substituted.'}
    target = destination / 'data/manifests/ecoli-full-v1/retrieval_metadata.json'
    target.write_text(json.dumps(dates, indent=2) + '\n')
    (destination / 'data/manifests/ecoli-full-v1/export_checksums.json').write_text(json.dumps(files, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'amp-local-release-export.zip')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Preserving existing output; choose a new ZIP filename')
    with tempfile.TemporaryDirectory(prefix='amp-release-export-') as folder:
        temporary = Path(folder)
        export_model(ROOT / PRETEST_REL, temporary / 'artifacts/model-v1')
        export_sources(ROOT, temporary)
        with args.output.open('xb') as handle:
            with zipfile.ZipFile(handle, 'w', zipfile.ZIP_DEFLATED) as archive:
                for path in sorted(temporary.rglob('*')):
                    if path.is_file():
                        archive.write(path, path.relative_to(temporary))
    print('Verified export:', args.output)
    print('Contains portable model, exact original manifest/source tables, and observed date metadata.')
    print('No model was fit; no final evaluation was repeated. Review/upload this ZIP to finish release packaging.')


if __name__ == '__main__':
    main()
