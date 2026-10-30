"""Check committed inputs/summary fingerprints without training or evaluation."""
import hashlib
import json
from datetime import datetime
from pathlib import Path
from release_paths import (ROOT, SPLIT_REL, COHORT_REL, LABEL_REL, BASELINE_REL,
                           COMPARISON_REL, PRETEST_REL, EVALUATION_REL, SENSITIVITY_REL)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(root=ROOT):
    root = Path(root).resolve()
    checked = []
    def check(path, expected):
        if sha(path) != expected:
            raise ValueError(f'Fingerprint mismatch: {path.relative_to(root)}')
        checked.append(str(path.relative_to(root)))
    split = json.loads((root / SPLIT_REL / 'split_manifest.json').read_text())
    for name, digest in split['output_sha256'].items():
        check(root / SPLIT_REL / name, digest)
    for name, digest in split['spec']['source_sha256'].items():
        check(root / COHORT_REL / name, digest)
    pretest = json.loads((root / PRETEST_REL / 'evaluation_manifest.json').read_text())
    check(root / SPLIT_REL / 'split_manifest.json', pretest['spec']['split_manifest_sha256'])
    check(root / COMPARISON_REL / 'comparison_report.json', pretest['spec']['comparison_report_sha256'])
    for name, digest in pretest['spec']['code_sha256'].items():
        path = root / name if name == 'prepare_splits.py' else root / 'modeling' / name
        check(path, digest)
    missing_generated = []
    for directory, report_name in [(BASELINE_REL, 'training_report.json'),
                                   (COMPARISON_REL, 'comparison_report.json'),
                                   (PRETEST_REL, 'evaluation_manifest.json'),
                                   (EVALUATION_REL, 'evaluation_report.json'),
                                   (SENSITIVITY_REL, 'sensitivity_report.json')]:
        report = json.loads((root / directory / report_name).read_text())
        for name, digest in report['output_sha256'].items():
            path = root / directory / name
            if path.exists():
                check(path, digest)
            else:
                missing_generated.append(str(directory / name))
    final = json.loads((root / EVALUATION_REL / 'evaluation_report.json').read_text())
    check(root / PRETEST_REL / 'evaluation_manifest.json', final['spec']['pretest_manifest_sha256'])
    check(root / 'modeling/evaluate_holdout.py', final['spec']['evaluator_sha256'])
    sensitivity = json.loads((root / SENSITIVITY_REL / 'sensitivity_report.json').read_text())['spec']
    for relative, key in [(LABEL_REL / 'build_report.json', 'build_report_sha256'),
                          (LABEL_REL / 'ecoli_measurements.csv', 'measurements_sha256'),
                          (COHORT_REL / 'provenance_report.json', 'provenance_report_sha256')]:
        check(root / relative, sensitivity[key])
    for name, key in [('build_labels.py', 'label_source_sha256'),
                      ('build_provenance.py', 'provenance_source_sha256'),
                      ('data_utils.py', 'data_utils_sha256')]:
        check(root / name, sensitivity[key])
    check(root / 'data/curation/full_note_reviews.csv', sensitivity['review_sha256'])
    build = json.loads((root / LABEL_REL / 'build_report.json').read_text())
    collection = root / 'data/manifests/ecoli-full-v1'
    check(collection / 'manifest.json', build['policy']['manifest_sha256'])
    cohort = json.loads((collection / 'manifest.json').read_text())
    ids = cohort['selected_ids']
    if (len(ids) != build['counts']['selected_records'] or len(ids) != len(set(ids))
            or any(type(value) is not int or value <= 0 for value in ids)):
        raise ValueError('Frozen cohort manifest has invalid ID membership')
    for relative, digest in json.loads((collection / 'export_checksums.json').read_text()).items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError('Export checksum path escapes the repository')
        check(path, digest)
    dates = json.loads((collection / 'retrieval_metadata.json').read_text())
    if dates['record_download_dates'] is not None:
        raise ValueError('Unrecorded peptide retrieval dates must remain unknown')
    if not dates['observed_metadata']:
        raise ValueError('No collection timestamp metadata recorded')
    for item in dates['observed_metadata']:
        if datetime.fromisoformat(item['value']).tzinfo is None:
            raise ValueError('Collection timestamps require a timezone')
    model_dir = root / 'artifacts/model-v1'
    model_manifest = json.loads((model_dir / 'model_manifest.json').read_text())
    check(model_dir / 'primary_model.json', model_manifest['files']['primary_model.json'])
    model = json.loads((model_dir / 'primary_model.json').read_text())
    if (model['source_model_sha256'] != pretest['primary_model_sha256']
            or model_manifest['source_model_sha256'] != pretest['primary_model_sha256']):
        raise ValueError('Portable model is not linked to the evaluated primary artifact')
    return {'checked_fingerprints': len(checked),
            'scope': 'Frozen inputs/sources, present summaries, original cohort metadata and portable model; no training/evaluation.',
            'missing_generated_artifacts': missing_generated}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
