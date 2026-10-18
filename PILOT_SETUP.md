# Replacement pilot tools

Copy data_utils.py, prepare_pilot.py, download_pilot.py, and audit_pilot.py
into the existing amp-designer folder. The shared module replaces the old one;
save your current work in Git before copying. Copy tests/test_pilot_tools.py
into tests if you want to run the offline checks. Keep your existing README,
environment, and exploratory records. The existing audit_sample.py,
audit_labels.py, and audit_metadata.py remain exploratory-sample tools;
use audit_pilot.py for pilot statistics. Do not mix their populations.

The scripts require certifi, already installed in your environment. No
Biopython or modeling libraries are required yet.

## Layout and snapshot identity

- data/raw/peptide_*.json: existing exploratory sample, left unchanged.
- data/raw/listings/<run>-<configuration-hash>/: timestamped listing envelopes.
- data/raw/pilot/<run>-<configuration-hash>/manifest.json: frozen selected IDs.
- data/raw/pilot/<run>-<configuration-hash>/records/: full pilot records.
- data/processed/: future cleaned tables and derived reports.

The configuration hash includes endpoint, filters, page size, seed, and sample
size. A new --run name creates a fresh retrieval even with the same parameters.
Cached listing pages validate their stored URL/query. A frozen manifest is
reused unchanged on resume. The manifest records the population IDs, but this
is a snapshot assembled across requests, not a server-side transactional
snapshot. Count and duplicate checks detect some changes, not all changes.

## Commands

Run from amp-designer with .venv activated:

    python -m unittest discover -s tests -p test_pilot_tools.py
    python prepare_pilot.py --run ecoli-v1 --size 300 --seed 42

The preparation command prints an absolute manifest path. Substitute that
exact printed path for MANIFEST_PATH in the commands below (keep quotes):

    python -u download_pilot.py "MANIFEST_PATH"
    python audit_pilot.py "MANIFEST_PATH"

To capture an audit report without modifying raw data:

    python audit_pilot.py "MANIFEST_PATH" > results/pilot_audit.txt

Make sure results exists first. An optional background download on macOS is:

    nohup python -u download_pilot.py "MANIFEST_PATH" > results/pilot_download.log 2>&1 &

Monitor results/pilot_download.log. Do not start concurrent download processes
for the same manifest. Ctrl+C in a foreground run retains completed files.
The 0.5-second delay is additional to network latency: 300 requests require at
least 150 seconds of delays and may take tens of minutes or longer. This bundle
does not claim that the live 300-record crawl has been completed or verified.

## Integrity and failures

save_json publishes complete JSON using a unique temporary file in the same
directory, flushes/fsyncs it, and uses a no-clobber hard link for immutable
writes. This works on the user's normal macOS local filesystem. Unsupported
hard-link filesystems produce an error rather than falling back to unsafe
writes. overwrite=True is an explicit option for future derived outputs;
pilot scripts never use it. A process kill can leave a hidden temporary file,
but the final JSON filename never refers to that partial file.

Cached full records must have the requested ID and a list of measurement
objects. Bad caches are reported and preserved, not silently skipped or
overwritten. Inspect/quarantine a bad cache or create a new run. Transient
download/JSON failures retry three times; exhausted failures produce a nonzero
exit status. The audit reports missing and invalid IDs, and withholds the
whole-pilot yield until every selected ID has a valid record.

If listing counts or unique IDs disagree, the database may have changed.
Prefer a new --run name to preserve provenance. You can deliberately remove
only that failed run's listing folder and rerun, but do not mix old and new
pages. Old offset-only listing files and the old raw-root manifest are ignored.

## Interpretation

The API search is for target species, not specifically for MIC. The audit
reports the fraction of selected records with at least one client-filtered
E. coli MIC. It reports exact names containing 'coli' across all endpoints and
exact names in the MIC subset. The client predicate remains provisional:
prefix matching may miss alternative names. A low yield triggers investigation;
it does not by itself establish a server defect. Validate search coverage,
metadata conventions, and exports before defining the final dataset.

The 300-record random pilot estimates structure within the API-returned
population; it is not yet a modeling dataset. Incomplete downloads can bias
observed yield, so the complete-pilot denominator is withheld. Sequence/MIC
intersection counts are preliminary, not chemical-eligibility counts.

Metadata states are absent, empty, and populated. Empty includes null, blank
strings, and empty containers, but not numeric zero. Neither absent nor empty
is translated into 'unmodified' or 'all L-amino acids'. Known chemistry fields
have state counts, and every discovered top-level key has state counts across
all valid downloaded records. Aliased terminal keys are preserved separately.
Assay states describe structured fields only; descriptions can contain further
conditions. No units are converted, no ranges are collapsed, and no labels
are aggregated by these scripts.

Before a full crawl, inspect DBAASP's CSV/activity export with species and
structure filters. It may cover the necessary measurements with fewer requests;
if it omits required metadata, retrieve only missing full records. Record data
usage terms and retrieval provenance. Settle cleaning/aggregation policy before
the full target-specific retrieval and keep original snapshots immutable.

Add observed pilot results and methods to the existing project README only
after obtaining them. In particular report selected/valid counts, server-query
parameters, client definition, filter yield, and the exclusions used later.
