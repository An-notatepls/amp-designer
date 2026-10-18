# Data preparation and labeling policy

**Version: v0.1 — draft. Updated: 2026-10-18.**

This document describes the implemented starting rules and the evaluation
decisions still pending. It is not a declaration that note review is complete
or that a v1 dataset/test split is frozen. It accompanies the
[project README](../README.md).

## 1. Task, units of observation and cohort hierarchy

The primary experiment predicts recorded E. coli MIC activity using the
reviewed broad dataset. Clinical-isolate measurements are eligible when they
otherwise meet the rules. This is not a reference-strain-only task.

- A **record** describes a database peptide and may contain many measurements.
- A **measurement row** describes one reported activity entry with its target,
  concentration and recorded experimental context.
- An **assay group** combines retained rows with matching grouping fields below.
  It does not prove that their unrecorded conditions were identical.
- A **peptide identity** is `(sequence, C-terminal amidation state)`. Different
  database IDs with the same modeled identity are combined at this stage.
- A **sequence family group** will keep related sequences together during
  evaluation. Its final definition is not yet frozen.

Laboratory/reference-only and clinical-isolate analyses are secondary. Their
candidate tables are not automatically evaluation partitions.

## 2. Collection, validation and traceability

The pilot manifest selected 300 records from a saved server-side E. coli search
containing 2,217 unique IDs. The full manifest expands to every ID in that saved
listing. Its scope is the returned ID collection, not the whole database or a
guarantee of an E. coli MIC measurement.

Process only manifest-selected records. Validate file readability, matching
peptide ID and activity-list structure. The label build checks selected-record
availability and measurement keys before writing outputs. Extra files are
ignored; missing/invalid selected records must not silently count as complete.

Keep raw data immutable. Use complete atomic JSON writes, validated caches,
separate processed directories, input hashes and recorded review decisions.
Record dates may differ; listing membership is frozen, but this does not create
a simultaneous historical snapshot of all record contents.

## 3. Initial peptide eligibility

Requirements:

1. A top-level sequence containing only `ACDEFGHIKLMNPQRSTVWY`.
2. Length from 10 through 50 residues, inclusive.
3. An explicit `Monomer` annotation.
4. No recorded N-terminal modification.
5. No recorded C-terminal modification other than `AMD` amidation.
6. Empty `unusualAminoAcids`, `intrachainBonds`, `interchainBonds` and
   `coordinationBonds` annotations.

Absent required chemistry fields and conflicting terminal aliases require
exclusion/review under the implemented rules. Present null/empty annotations
are accepted as **no modification recorded**, not verified unmodified chemistry.
`opms` entries link to structural information and are not modification flags.

Do not uppercase sequences or remove noncanonical characters to force a pass.
Lowercase/D-residue representations and unusual-residue annotations are outside
the initial representation. Incomplete annotations can still conceal chemical
differences; the filter is not independent stereochemical verification.

The 10–50 window and simplified chemistry keep the undergraduate baseline
focused. They were not selected using model/test performance and are not claims
about which peptides can have biological activity.

## 4. Target and assay inclusion

The provisional client E. coli MIC check requires `activityMeasureValue == MIC`
and a target name equal to `Escherichia coli` or beginning with that name plus a
space. Broad substring matches such as `coli` are not used for inclusion.

The implemented modeling target-name list is:

- `Escherichia coli`
- `Escherichia coli ATCC 25922`
- `Escherichia coli ATCC 8739`
- `Escherichia coli MG1655`
- `Escherichia coli K-12`
- `Escherichia coli K12 MG1655`

This exact-name list is a scope restriction, not a complete taxonomy or a
wild-type/reference-only definition. Generic species entries may describe
clinical isolates in their notes. No strain aliases are silently inferred.

Allowed medium names are `MHB` and `CAMHB`; they remain separate. Missing/blank
medium names and other media are excluded from this starting policy. A broth
name alone does not establish identical protocols or absence of supplements.
Treatment-related notes and contradictory annotations require review.

## 5. Units and molecular weight

Normalize unit names with Unicode NFKC, whitespace trimming and case folding.
Supported units are µM, µg/mL and mg/L. Unsupported units are not guessed.

For eligible peptide chemistry, compute average molecular weight using
Biopython `ProteinAnalysis(sequence).molecular_weight()`. For C-terminal
amidation, add the approximate average-mass correction **−0.9848 g/mol**.

For mass concentrations:

```text
MIC in µM = 1000 × MIC in µg/mL ÷ molecular weight in g/mol
```

mg/L is numerically equivalent to µg/mL for this conversion. Retain raw units,
concentrations and the molecular weight used. Convert bounds individually.
Nonpositive, reversed-range or nonfinite numeric values are not modeling
targets. The assumed mass applies to the represented peptide chemistry and
does not independently establish purity, counterion content or annotation
completeness.

## 6. Concentration interpretation

| Reported format | Initial treatment |
|---|---|
| Single positive number, including scientific notation | Point value for classification and regression |
| Numeric range, such as `8–16` | Preserve both endpoints; classification only |
| `>`, `>=`, `<`, `<=`, `≥`, `≤` | Preserve censoring direction and bound; classification only |
| Central value plus `±` uncertainty | Preserve central value and uncertainty; exclude from initial modeling |
| `>>`, `<<`, mixed forms such as `1->64`, missing or unrecognized text | Preserve and flag; do not force a point estimate |

The primary classifier threshold is **32 µM**. An exact value at the threshold
is active. A bound/range is classified only if it establishes a definite side:

| Example at 32 µM | Label |
|---|---|
| `32` | Active |
| `8–16` | Active |
| `16–64` | Unresolved |
| `>32` | Inactive |
| `>=32` | Unresolved: the true value could equal 32 |
| `<32` or `<=32` | Active |
| `>8` | Unresolved |

This is an operational research label, not a clinical breakpoint or a guarantee
of activity against every E. coli strain. The concentration string is used;
the separate numeric `activity` field is not silently substituted for it.

## 7. Assay notes and recorded review decisions

Every nonempty activity note starts as pending review. Store decisions in a
tracked review CSV, matched to peptide ID, measurement ID and the exact current
note. Accept only `include` or `exclude` with a reason; a blank decision stays
pending. Stale/duplicate review keys must not be accepted silently.

An ordinary single-isolate or laboratory-strain description may be included
under the broad task. It does not verify the entire assay against its paper.
Multiple-isolate summaries, percentile summaries, conflicting descriptions and
notes mentioning extra conditions may require the original source. Phrases
such as "also with" do not establish which condition produced a specific MIC.

Retain pending cases rather than treating them as invalid measurements.
Reviewing an origin note and establishing the meaning of an aggregated
concentration are separate decisions. The current combined decision file is
`data/curation/full_note_reviews.csv`.

## 8. Grouping and broad classifier labels

The label builder groups by peptide ID, exact target strain and the recorded
JSON representations of medium, note, inoculum (`cfu`), pH, salt type, ionic
strength and temperature. Different notes must not be pooled. Reference IDs
are retained but are not themselves grouping keys. Matching raw metadata does
not prove independent experiments or identical unrecorded conditions.

Explicitly excluded review rows do not enter aggregation. Within a group,
opposing definite labels produce conflict; unresolved evidence prevents a
definite group label under the conservative policy.

Groups are then combined by sequence plus amidation state. The primary broad
classifier includes only identities with consistent definite group evidence
and no pending note review. Conflict/mixed status takes precedence; review
and unresolved-evidence flags are also retained. No majority vote is used.
An unreviewed eligible row can block an identity even when another group has
a definite label. This deliberate conservatism may introduce selection bias.

Assay tables and classifier tables have different row meanings. Never treat
the number of assay groups as the number of independent peptide sequences.

## 9. Provenance-specific experiments

`build_provenance.py` uses explicit annotations to distinguish:

- **Laboratory/reference:** one of the five named targets above, excluding the
  generic species entry, or an explicit laboratory-strain note.
- **Clinical isolate:** an explicit clinical-origin note.
- **Unspecified:** no positive origin evidence. An empty note, `WT`, `MDR`,
  resistance or pathogenicity description alone does not establish origin.

A clinical note conflicting with a named laboratory/reference target or
laboratory note is flagged for source review and excluded from the narrower
provenance branches. These assignments are annotation-based, not publication-
verified classifications or guarantees of wild-type phenotype.

For each narrower branch, select base-eligible, reviewed measurement rows first
and recompute labels from the remaining groups. A narrow-cohort identity may
be consistent while its broader evidence is mixed. Retain that fact explicitly
as `mixed_in_broad_dataset`; do not describe removed measurements as wrong.

Clinical-only candidates require exclusively clinical-origin evidence among
the captured base-eligible rows, with no provenance conflict. This does not
prove the peptide has never been tested against another strain. Candidate
tables are not locked holdouts and must not automatically enter training.

## 10. Regression is a separate, provisional target

Only eligible single numeric MICs contribute to median log10 MIC summaries
within an assay group. Ranges, censored values and `±` central values do not
become point-value regression observations. Preserve the number of numeric
measurements and their spread where available.

Stored regression summaries are not a finalized regression training dataset.
Their assay notes and cohort definition still require review. A classifier
exclusion must not be assumed to automatically remove every corresponding
regression summary already stored in an audit table.

## 11. Evaluation and freeze requirements

Before training or evaluating:

1. Confirm the primary broad candidate table and its input fingerprints.
2. Freeze reviewed decisions and leave unresolved cases explicitly pending.
3. Define sequence-family grouping, including how chemical variants are linked.
4. Record a reproducible, group-disjoint split and training-only model-selection
   procedure. Never select the split based on model performance.
5. Decide whether sufficient independent reference and clinical groups support
   a secondary transfer experiment; identify its held-out examples in advance.
6. Record the version, freeze date, cohort counts, group parameters, random seed
   and exact split membership. The policy remains a draft until then.

An 80% normalized edit-similarity cutoff has been explored, but is not yet a
frozen definition. It measures one minus unit-cost edit distance divided by the
longer sequence length; it is not a standard biological percent-identity claim.
Connected grouping can link sequences through intermediates, and does not
establish homology. Exact sequences and amidation variants must remain together.

Use the same sequence-group assignments for the predefined 16, 32 and 64 µM
sensitivity tasks and across compared training variants. Label eligibility can
change with the threshold; report those changes. Do not select the threshold
or final model using held-out test performance.

For clinical robustness testing, held-out clinical examples and their close
relatives must be absent from both training variants. A model trained with
clinical data can use only the clinical examples assigned to its training
groups. Distinguish transfer to new peptides from testing a familiar peptide
under a new strain/assay context.

## 12. Documentation, dates and data terms

Keep raw measurements, review reasons, source identifiers and build reports.
Commit original code, documentation, review decisions and release metadata;
generate raw caches and processed tables locally under the current repository
convention. Export frozen manifest metadata and split definitions to a tracked
location before claiming fresh-checkout reproducibility.

Transcribe actual listing/record access dates from manifests, cache envelopes
and retrieval logs. Reused records may predate the remaining downloads.
Documentation and terms-check dates are not dataset retrieval dates.

DBAASP's [terms](https://dbaasp.org/terms-and-conditions), checked 2026-10-18,
require public source acknowledgement when its data or tools are used. Its
requested publication is Pirtskhalava et al., *DBAASP v3*, 2021,
[doi:10.1093/nar/gkaa991](https://doi.org/10.1093/nar/gkaa991).
Any future code license is separate from the source data's terms.
