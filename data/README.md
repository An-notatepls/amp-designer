# Data attribution and release scope

Source: DBAASP, https://www.dbaasp.org/, and the original experimental
publications linked by its records. Please cite Pirtskhalava et al. (2021),
DBAASP v3, Nucleic Acids Research 49(D1):D288–D297,
https://doi.org/10.1093/nar/gkaa991.

The website Data Access Policy, checked 2026-10-28, permits copying, adaptation
and redistribution with public acknowledgement identifying DBAASP as source:
https://www.dbaasp.org/terms-and-conditions . Older downloadable policy text
contains differing language; this release records the current website policy
rather than claiming all historical texts agree. Source data are not relicensed
under this repository's MIT code license.

Committed derived data are filtered assay-derived identity labels, provenance
assignments and frozen evaluation partitions. Cleaning changes include chemistry,
length, target, unit, assay and manual-note filters, label aggregation at 32 µM,
and sequence-family splitting. They are not raw DBAASP records or a complete
database mirror. See docs/data_policy.md for assumptions and limitations.

The four partition CSVs and three provenance tables match the fingerprints in
the original split manifest. They permit modeling reproduction while preserving
original inputs. Test/clinical outcomes are now public/evaluated; do not reuse
them as a fresh blind test for future model development.

The original full 2,217-ID manifest and collection timestamp metadata must be
exported from the original local cache with tools/export_local_release.py.
Dates are not inferred from filesystem modification times. Per-record retrieval
dates were not recorded by the downloader; that missing history is disclosed.
