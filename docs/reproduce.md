# Reproduction and inference

Updated 2026-10-28. Commands run from the repository root. The fixed experiment
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

There are 63 tests, including seven new release checks. Test fixtures are
synthetic except a fingerprint consistency check over committed artifacts.
The verifier checks the exact frozen partitions, provenance inputs, source
fingerprints and present report exports. It lists uncommitted generated files
separately; their absence does not mean that every original artifact is present.
Hosted CI on Ubuntu 24.04 / Python 3.13.4 passed all 63 tests, frozen-input
verification and the full modeling reproduction; see [the successful run](https://github.com/An-notatepls/amp-designer/actions/runs/40631386045).
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
intervals at numerical tolerance 1e-10. It also runs the existing grouping
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

## Export the original model and raw-cohort provenance

The exact original primary joblib bundle, raw 2,217-ID manifest, cached listing
timestamps and complete measurement/build reports remain on the original Mac.
They must not be recreated from the now-known test results. The export tool
checks them against frozen fingerprints before producing portable files:

```bash
./.venv/bin/python tools/export_local_release.py
```

This creates a new `amp-local-release-export.zip`. It contains:

- A JSON scaler/logistic parameter model and its checksum manifest, exported
  from the exact evaluated primary model. A 204-point charge/amidation grid covering the supported
  feature combinations are checked against its decision function without fitting.
- The unmodified original full ID manifest under `data/manifests/ecoli-full-v1/`.
- Original measurement, build and provenance files with unchanged fingerprints.
- The listing/manifest timestamp values actually present in the local cache.
  Individual record-download timestamps were not saved and remain unknown.

The exporter only loads the trusted original local joblib artifact. Do not
load pickle/joblib from an untrusted source. Keep the ZIP outside staging;
review and add its extracted intended files explicitly. Uploading this ZIP is
the remaining handoff needed to finish the original-model/raw-cohort release.
There is no invented release model or backdated dataset access date in this branch.

After exporting and placing the portable files under `artifacts/model-v1/`:

```bash
python3 modeling/predict_portable.py --sequence AAAAAAAAAA --unmodified
```

This inference command needs only Python's standard library, reads JSON rather
than pickle, checks the model checksum, and never fits a model. Before that
export is installed, it reports a missing artifact rather than making up a
prediction. Scores are uncalibrated and do not establish experimental activity.
The synthetic sequence demonstrates syntax only. Use `--amidated` when applicable;
canonical uppercase 10–50-residue monomers and the supported chemistry are required.

Your existing local model command remains available unchanged:

```bash
./.venv/bin/python modeling/predict.py --sequence AAAAAAAAAA --unmodified
```

Its batch mode takes a CSV with `sequence,c_amidated` and states exactly
`True` or `False`, plus an optional new `--output` CSV path.

## Level 3: best-effort DBAASP rebuild

**Not yet demonstrated from a clean clone.** The frozen raw manifest must first
be exported as above. Per-record download dates cannot be recovered from the
current cache if they were never recorded. Documentation/export dates and file
modification times must not be substituted for retrieval dates.

Once available, copy the manifest's exact bytes into the original
`data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json` location without
replacing existing files. Then the original commands are:

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
original build report's `record_sha256` mapping. A mismatch is source drift,
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
