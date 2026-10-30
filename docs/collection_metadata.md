# Collection metadata and source snapshots

Documentation updated 2026-10-29. These dates come from the original exported
manifest and cached listing metadata. They are not inferred from file modification
times or the date this repository was documented.

| Event | Observed timestamp (UTC) | Meaning |
|---|---|---|
| Listing run started | 2026-10-11 19:45:54.378162 | Start of saved filtered ID-list collection |
| First listing page metadata | 2026-10-11 19:45:54.378398 | First observed listing request timestamp |
| Last listing page metadata | 2026-10-11 19:46:41.624464 | Last observed listing request timestamp |
| Source pilot manifest prepared | 2026-10-11 19:46:42.560386 | Listing/pilot-selection metadata prepared |
| Release export documented | 2026-10-30 03:00:25.164734 | Local export operation; **not a retrieval date** |

There are 23 observed listing-page timestamps for the saved 2,217-ID
server-filtered E. coli collection. The full original manifest is at
[`data/manifests/ecoli-full-v1/manifest.json`](../data/manifests/ecoli-full-v1/manifest.json).
The exported [`retrieval_metadata.json`](../data/manifests/ecoli-full-v1/retrieval_metadata.json)
preserves exact values, source fields, original listing-cache paths and hashes.
The cache page files themselves are not included in this release.

**Individual peptide-record retrieval dates are unknown.** The old downloader
saved raw record JSON without retrieval envelopes. Some records were reused
from the pilot and others downloaded later. Listing timestamps must not be
presented as access dates for every peptide record, and the cohort is not a
simultaneous database snapshot.

## Frozen sources now available

The original measurement CSV, label build report and provenance report are
committed at their recorded processed-data paths. Each matches the fingerprints
in the original sensitivity/provenance reports. The original build report also
records a hash for each selected raw peptide file; those hashes permit later
drift checks but do not restore missing historical raw content.

`export_checksums.json` preserves checksums for the four unmodified source
copies. The portable model is linked to the pretest primary joblib fingerprint
and has its own JSON checksum. No executable pickle or joblib file is required
for inference. See [reproduction](reproduce.md) for the independent scopes of
summary verification, modeling reproduction and a best-effort live raw rebuild.
