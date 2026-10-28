# Reproduction and prediction usage

Run commands from `amp-designer`, not from inside `modeling`. This guide refers
to the recorded local experiment. An ordinary clone does not include ignored
raw/processed data or fitted joblib bundles. Complete clean-checkout data and
model recovery is still a release task; do not claim it is already available.

## Environment and tests

The recorded model environment is Python 3.13.4, numpy 2.3.5, scipy 1.17.0,
scikit-learn 1.8.0, joblib 1.5.3, threadpoolctl 3.6.0 and Biopython 1.86.
Use the existing environment; do not upgrade it for this frozen experiment.
For a separate reconstruction environment with the matching Python:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r modeling/requirements-diagnostics.txt
```

Data collection tools also require their existing dependencies (including
certifi). Complete data-tool dependency pins remain a release task.

```bash
./.venv/bin/python modeling/test_predict.py -v
```

Inference tests use synthetic data only. Existing split, baseline, comparison,
pretest and evaluator tests can be run with their individual test scripts.
The prediction test imports `compare_features.py` to check exact feature parity.

## Predict one sequence

Keep the original, trusted `primary_model.joblib` and `evaluation_manifest.json`
in `results/pretest/pretest-v1-cbd14c16acf65e99/`. Do not replace the frozen
model with a newly fit or all-data model. joblib uses pickle; load only trusted
artifacts you produced. Fingerprint checks detect changes relative to the local
manifest, but do not authenticate an arbitrary model/manifest downloaded together.

```bash
./.venv/bin/python modeling/predict.py --sequence AAAAAAAAAA --unmodified
./.venv/bin/python modeling/predict.py --sequence AAAAAAAAAA --amidated
```

These synthetic inputs demonstrate syntax, not recommended candidates. The
tool requires a known terminal state; it does not infer chemistry from sequence.
Use only canonical uppercase 10–50-residue monomers without additional chemistry.
It cannot verify D-residues, bonds or chemical modifications from sequence alone.

The output includes a charge proxy, model probability, binary prediction,
activity definition and model fingerprint. A model probability is not a
validated confidence estimate or an experimental success rate. The model can
assign the same score to many unrelated peptides with identical features; this
is not a generator or a toxicity/selectivity predictor. It does not calculate
sequence novelty to the training set.

The default pretest directory can be overridden with `--pretest-dir`. The tool
requires the model's original numerical-library environment and training-feature
source fingerprint. Biopython is not used by inference itself. Invalid sequences,
unknown chemistry, changed model bytes and existing output files are rejected.

## Predict a CSV

Create a CSV with explicit chemistry for every row:

```csv
sequence,c_amidated
AAAAAAAAAA,False
AAAAAAAAAA,True
```

```bash
./.venv/bin/python modeling/predict.py --input-csv peptides.csv --output predictions.csv
```

CSV states must be exactly `True` or `False`. Extra input columns are ignored
and not propagated. Output rows preserve input order, including duplicates.
The output path's parent folder must exist; an existing file is not overwritten.

## Recorded data-to-model workflow

The full raw manifest is
`data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json`.
Its selected cached records must exist. Manifest preparation and download
scripts are `prepare_pilot.py`, `prepare_full.py` and `download_pilot.py`;
`audit_pilot.py` accepts full-cohort manifests despite its historical name.
Use their `--help` for collection arguments. A fresh live query does not recreate
the saved historical collection.

```bash
python audit_pilot.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json"
python build_labels.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json" --threshold 32 --note-reviews data/curation/full_note_reviews.csv
```

Use the actual output directory printed by the label builder as `LABEL_DIR`,
then the printed provenance path as `COHORT_DIR`. The following are shell
variables assigned to the recorded local paths; rebuilt fingerprints can differ.

```bash
LABEL_DIR="data/processed/ecoli-full-v1-3cec454243584936/45c31cc2e742"
COHORT_DIR="$LABEL_DIR/provenance/83780b72da84"
SPLIT_DIR="$COHORT_DIR/splits/grouped-v1-dbe608485cce9a7e"
```

The stage order and command shapes are:

```bash
python build_provenance.py "$LABEL_DIR" --note-reviews data/curation/full_note_reviews.csv --threshold 32
python prepare_splits.py "$COHORT_DIR"
python modeling/train_baseline.py "$SPLIT_DIR" --results-root results/baselines
```

Use the baseline output actually printed by that run:

```bash
python modeling/compare_features.py "$SPLIT_DIR" --baseline-dir "results/baselines/physchem-logistic-v1-b8379a4f32704165" --results-root results/comparisons
python modeling/prepare_evaluation.py "$SPLIT_DIR" --comparison-dir "results/comparisons/feature-comparison-v1-3c0530a63f1c77d0" --results-root results/pretest
python modeling/check_label_sensitivity.py "$SPLIT_DIR" --note-reviews data/curation/full_note_reviews.csv
```

These commands explain the pipeline; do not regenerate source bundles merely
to improve the evaluated experiment. The frozen pretest protocol, subset
membership and evaluator were committed before final evaluation. Recreating
the historical chronology requires the preserved commits/artifacts; a new
pretest run made after seeing results cannot retroactively preregister choices.

The final evaluator has already been run for this experiment. Its completed
outputs are the authoritative results. Preserve them; no command is needed to
repeat model selection or reopen the outcomes. Identical completed evaluator
runs verify and reuse cached results. `--resume` is for an interrupted identical
run, not a way to change settings after results are known.

## Tracked versus local artifacts

Commit source, tests, documentation, reviewed decisions, dependency specs,
fixed split definitions, saved CV assignments, pretest protocol/subsets and small
result reports. Raw caches, cleaned data tables and joblib models remain local
under current ignore rules. Preserve original evaluation predictions locally
for audit; decide data redistribution separately under the source terms.

A future release must supply an explicit trusted model/data recovery mechanism,
record actual collection dates from manifests/logs, complete dependency pins,
and add a code license. Dataset access dates and the source data terms are
separate from code licensing and documentation-update dates.
