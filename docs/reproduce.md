# Reproduction and inference

Updated 2026-10-29. Commands run from the repository root. The fixed experiment
is preserved; this guide does not replace original reports, rewrite history,
select another model or reopen final outcomes for tuning.

## Environment roles

`.python-version` records **Python 3.13.4**, the original experiment interpreter.
The model/parameter exporter requires the original frozen environment; the
summary notebook has a separate environment. Dependencies are layered:

| File | Purpose |
|---|---|
| `requirements.txt` | Complete current data/model workflow: certifi plus diagnostic/model requirements |
| `modeling/requirements-baseline.txt` | Pinned numerical/modeling libraries |
| `modeling/requirements-diagnostics.txt` | Baseline requirements plus pinned Biopython for alignment/data chemistry |
| `notebooks/requirements-notebook.txt` | Separate notebook display/execution environment |

The requirements pin direct dependencies, not every transitive dependency.
For a fresh modeling environment with the recorded Python:

```bash
python3 -m venv .venv
./.venv/bin/python -m pip install -r requirements.txt
```

Confirm `python3 --version` is 3.13.4 first if reproducing the original
interpreter exactly. Do not upgrade an existing frozen modeling environment
merely to display the notebook.

## Level 1: offline tests and integrity

```bash
./.venv/bin/python -m unittest discover -s tests -p 'test_*.py' -v
./.venv/bin/python -m unittest test_prepare_splits.py -v
./.venv/bin/python -m unittest discover -s modeling -p 'test_*.py' -v
./.venv/bin/python tools/verify_release.py
```

There are 64 tests, including eight release checks. Fixtures are synthetic
except frozen-artifact integrity checks and an isolated invocation of the
committed portable model with its recorded synthetic syntax-example prediction.
The verifier checks the exact frozen partitions, provenance inputs, source
fingerprints and present report exports. It lists uncommitted generated files
separately; their absence does not mean that every original artifact is present.
Hosted CI on Ubuntu 24.04 / Python 3.13.4 runs all 64 tests, frozen-input
verification and the full modeling/label-sensitivity reproduction; see the
[latest PR checks](https://github.com/An-notatepls/amp-designer/actions).
The earlier 63-test/full-modeling check also passed in
[this recorded run](https://github.com/An-notatepls/amp-designer/actions/runs/40631386045).
Local authoring tests also passed on Python 3.12.14 with the recorded numerical
library versions.

## Level 2: training-only modeling reproduction

The original `train.csv`, `test.csv`, clinical and buffer partitions are now
included as exact bytes matching `split_manifest.json`. The three provenance
tables also match that manifest's source hashes. Do not replace them with newly
queried or relabeled data. The existing `split_assignments.csv` has no labels;
it is assignment metadata rather than a complete training dataset.

```bash
./.venv/bin/python tools/reproduce_modeling.py
```

This uses the frozen training and comparison scripts. It writes a new run under
`results/reproduced/`, checks exact saved CV folds and selected C, and compares
baseline CV scores/coefficients, the feature ladder and paired bootstrap
intervals at numerical tolerance 1e-10. It also rebuilds the training-only
16/32/64-µM label-sensitivity exports and checks their original fingerprints. It also runs the existing grouping
sensitivity diagnostic. No test/clinical outcome files are opened by training.
The initial integrity verifier hashes them as immutable inputs, which is
explicitly different from interpreting their outcomes for model selection.

Full local reproduction matched the committed numeric results on Python
3.12.14 with numpy 2.3.5, scipy 1.17.0, sklearn 1.8.0, joblib 1.5.3 and
threadpoolctl 3.6.0. The same full numerical reproduction also passed hosted CI with Python
3.13.4 on Ubuntu 24.04. Floating-point platform variation can produce differences; the verifier
fails rather than widening tolerance silently. New run directory names,
timestamps and serialized model bytes need not match original report hashes.
Original hash-linked reports are never hand-edited or overwritten.

For a quicker check:

```bash
./.venv/bin/python tools/reproduce_modeling.py --baseline-only
```

The completed final metrics are preserved under
`results/evaluation/holdout-evaluation-v1-836296044b94ac62/`. This reproduction
command does not rerun final evaluation. These public/evaluated outcomes are
not a new blind test for future model development.

## The read-only walkthrough

Open `notebooks/project_walkthrough.ipynb` to read its saved tables and figures.
It runs entirely from `notebooks/data/walkthrough/`, with no database, model,
network access or fitting. For execution, create a separate environment:

```bash
python3 -m venv .venv-notebook
./.venv-notebook/bin/python -m pip install -r notebooks/requirements-notebook.txt
```

Select that interpreter as the VS Code notebook kernel and use Run All. The
notebook has 12 code cells and six figure outputs. Its code was validated in
sequential in-process IPython because authoring-kernel socket binding was
unavailable; local normal-kernel Run All remains to be confirmed.

## Prediction from a fresh clone

The original evaluated model is now committed as portable scaler/logistic
parameters under `artifacts/model-v1/`. Inference needs only Python's standard
library; it reads JSON rather than pickle and never fits a model:

```bash
python3 modeling/predict_portable.py --sequence AAAAAAAAAA --unmodified
```

The synthetic input demonstrates syntax only. It produces an inactive label
with score approximately 0.16963, matching the original local CLI. Use
`--amidated` for C-terminal amidation. Inputs must be canonical uppercase
10–50-residue monomers with the supported chemistry. Model scores are
uncalibrated and do not establish experimental activity, toxicity or efficacy.

`model_manifest.json` verifies the JSON checksum and records the original
primary joblib SHA-256. The original-machine exporter checked the exact source
bundle before extracting parameters and compared decisions across a 204-point
charge/amidation grid. The received export was independently compared with the
previously reproduced training-only model on the same grid within 1e-12.
There was no refitting, model selection or repeated held-out evaluation.

The JSON model/manifest are versioned together; a checksum detects changes
relative to that manifest, rather than authenticating arbitrary downloads.
An isolated-process test runs the real portable CLI with `python -I`, without
third-party site packages.

The original joblib command still works on the original machine:

```bash
./.venv/bin/python modeling/predict.py --sequence AAAAAAAAAA --unmodified
```

Its batch mode takes a CSV with `sequence,c_amidated` and states exactly
`True` or `False`, plus an optional new `--output` CSV path. A fresh clone does
not include the original joblib bundle; use the portable command there.

## Original-machine export and provenance

The supplied export is now integrated: portable model JSON, the unmodified
2,217-ID full cohort manifest, original measurement/build/provenance files,
checksums and actual observed listing timestamps. See
[collection metadata](collection_metadata.md). The listing window is
2026-10-11 19:45:54–19:46:41 UTC. Individual peptide retrieval dates were not
recorded and remain unknown; export/documentation dates are not substitutes.

`tools/export_local_release.py` is retained as the provenance/recovery utility.
It requires the original trusted joblib and raw cache, so it is not a new-clone
setup step and should not be rerun merely to use the committed model. It never
fits a model or rewrites the frozen source reports.

## Level 3: best-effort DBAASP rebuild

**A new live crawl has not been demonstrated from a clean clone.** The frozen
raw manifest and historical record checksums are now committed. Per-record download dates cannot be recovered from the
current cache if they were never recorded. Documentation/export dates and file
modification times must not be substituted for retrieval dates.

The first download commands below are optional and make network requests.
Use a fresh cache for a new crawl; preserve any original local record cache.
Restore the saved manifest to the expected raw namespace without overwriting
an existing file:

```bash
mkdir -p data/raw/full/ecoli-full-v1-3cec454243584936
cp -n data/manifests/ecoli-full-v1/manifest.json data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json
```

Run `python tools/verify_release.py` for committed-source integrity first. The
original download/audit commands are:

```bash
python download_pilot.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json"
python audit_pilot.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json"
python audit_core.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json"
python build_labels.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json" --threshold 32 --note-reviews data/curation/full_note_reviews.csv
```

Use the actual label folder printed, then run `build_provenance.py` with that
folder, `--threshold 32` and the same review file. `prepare_pilot.py` prepares a
new random pilot; `prepare_full.py` expands its saved listing. Those are new
collection workflows, not historical reproduction of the existing cohort.

The current downloader validates cached IDs/structure, resumes, retries and
publishes complete immutable JSON files. Its historical names are retained for
compatibility. Use its `--help`; do not mix exploratory and manifest-defined files.
For a background download, redirect output to a new log and never start two
writers for the same manifest. Ctrl+C preserves completed validated records.

DBAASP records can change or disappear. Compare new record bytes with the
original build report's `policy.record_sha256` mapping (keys are peptide IDs). A mismatch is source drift,
not a reason to replace the frozen experiment or silently alter expected hashes.
A current crawl cannot guarantee a byte-identical historical snapshot.

## History and data terms

The pretest protocol was committed in b884ffe and evaluator in 387a4f5, before
final results in 526f231. The runtime marker was created at evaluation start
and committed with results. Preserve that history; directory names alone do
not establish chronology. No historical commits or report hashes are rewritten.

Original code is MIT licensed. Derived DBAASP tables retain attribution and
source-data terms; see `data/README.md`. The pilot review CSV is retained as
historical curation evidence; the active build uses the full review file.
Obsolete exploratory code was removed from the current tree and remains in Git
history. The active manifest-aware audit commands and test dependencies remain.
