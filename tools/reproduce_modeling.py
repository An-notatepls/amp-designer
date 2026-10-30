"""Rebuild training-only results in a new directory and compare numerical summaries.

Never opens test/clinical outcomes, rewrites reports, or changes model choices.
Reports/model binaries are not expected to have identical timestamped bytes.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path
from release_paths import ROOT, SPLIT_REL, BASELINE_REL, COMPARISON_REL, SENSITIVITY_REL
from verify_release import verify
sys.path.insert(0, str(ROOT / 'modeling'))
import train_baseline as base
import compare_features as compare


def compare_rows(expected, actual, keys, tolerance=1e-10):
    def load(path):
        with Path(path).open(newline='', encoding='utf-8') as handle:
            return list(csv.DictReader(handle))
    left, right = load(expected), load(actual)
    def index(rows):
        result = {tuple(row[key] for key in keys): row for row in rows}
        if len(result) != len(rows):
            raise ValueError('Duplicate comparison keys')
        return result
    left, right = index(left), index(right)
    if left.keys() != right.keys():
        raise ValueError(f'Result row membership differs: {expected}')
    for key in left:
        if left[key].keys() != right[key].keys():
            raise ValueError('Result columns differ')
        for column, value in left[key].items():
            other = right[key][column]
            if value == other:
                continue
            try:
                match = math.isclose(float(value), float(other), rel_tol=tolerance, abs_tol=tolerance)
            except ValueError:
                match = False
            if not match:
                raise ValueError(f'Numerical reproduction mismatch: {expected.name}, {key}, {column}: {value} != {other}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results/reproduced')
    parser.add_argument('--baseline-only', action='store_true')
    args = parser.parse_args()
    verify()
    output = args.output_dir.resolve()
    if any(output == (ROOT / p).resolve() or (ROOT / p).resolve() in output.parents
           for p in (BASELINE_REL, COMPARISON_REL)):
        parser.error('Choose a new output directory outside the original runs')
    old = json.loads((ROOT / BASELINE_REL / 'training_report.json').read_text())
    print('Recorded environment:', json.dumps(old['spec']['environment']), flush=True)
    print('Current environment:', json.dumps(base.versions()), flush=True)
    new_base, report = base.train_baseline(ROOT / SPLIT_REL, output / 'baselines', progress=True)
    if report['selected_C'] != old['selected_C']:
        raise ValueError('Selected baseline C differs')
    compare_rows(ROOT / BASELINE_REL / 'cv_scores.csv', new_base / 'cv_scores.csv', ['model', 'C', 'fold'])
    if (ROOT / BASELINE_REL / 'cv_folds.csv').read_bytes() != (new_base / 'cv_folds.csv').read_bytes():
        raise ValueError('Saved folds differ')
    compare_rows(ROOT / BASELINE_REL / 'coefficients.csv', new_base / 'coefficients.csv', ['feature'])
    print('PASS: baseline folds, selected C, CV metrics and coefficients match.', flush=True)
    if not args.baseline_only:
        new_comp, _ = compare.run_comparison(ROOT / SPLIT_REL, new_base, output / 'comparisons', progress=True)
        compare_rows(ROOT / COMPARISON_REL / 'feature_ladder.csv', new_comp / 'feature_ladder.csv', ['model'])
        compare_rows(ROOT / COMPARISON_REL / 'paired_differences.csv', new_comp / 'paired_differences.csv', ['model', 'previous_model', 'metric'])
        print('PASS: feature ladder and paired intervals match.', flush=True)
        import check_label_sensitivity as sensitivity
        new_labels, _ = sensitivity.check(ROOT / SPLIT_REL,
                                         ROOT / 'data/curation/full_note_reviews.csv',
                                         output / 'label_sensitivity')
        original_labels = json.loads((ROOT / SENSITIVITY_REL / 'sensitivity_report.json').read_text())
        for name, digest in original_labels['output_sha256'].items():
            if base.sha_file(new_labels / name) != digest:
                raise ValueError(f'Training-label sensitivity export differs: {name}')
        print('PASS: training-only 16/32/64-uM sensitivity exports match original fingerprints.', flush=True)
    print('Training-only reproduction complete. Original reports and held-out outcomes were not changed.')


if __name__ == '__main__':
    main()
