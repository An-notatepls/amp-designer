"""Recorded v1 paths; original hash-linked artifacts remain unchanged."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
LABEL_REL = Path('data/processed/ecoli-full-v1-3cec454243584936/45c31cc2e742')
COHORT_REL = LABEL_REL / 'provenance/83780b72da84'
SPLIT_REL = COHORT_REL / 'splits/grouped-v1-dbe608485cce9a7e'
BASELINE_REL = Path('results/baselines/physchem-logistic-v1-b8379a4f32704165')
COMPARISON_REL = Path('results/comparisons/feature-comparison-v1-3c0530a63f1c77d0')
PRETEST_REL = Path('results/pretest/pretest-v1-cbd14c16acf65e99')
EVALUATION_REL = Path('results/evaluation/holdout-evaluation-v1-836296044b94ac62')
SENSITIVITY_REL = Path('results/label_sensitivity/training-label-sensitivity-v1-d536b5d4387dd9f6')
RAW_MANIFEST_REL = Path('data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json')
