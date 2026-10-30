"""Check committed inputs/summary fingerprints without training or evaluation."""
import hashlib
import json
from pathlib import Path
from release_paths import (ROOT, SPLIT_REL, COHORT_REL, BASELINE_REL,
                           COMPARISON_REL, PRETEST_REL, EVALUATION_REL, SENSITIVITY_REL)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(root=ROOT):
    root = Path(root)
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
    return {'checked_fingerprints': len(checked),
            'scope': 'Frozen modeling inputs, frozen sources, and present summary exports; no training/evaluation.',
            'missing_generated_artifacts': missing_generated}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
