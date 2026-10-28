# Data preparation and labeling policy

**Version: v0.1 — frozen for the reported experiment. Documentation updated: 2026-10-27.**

The fixed sequence-group split, primary model and evaluation protocol have now
been used for final evaluation. Pending reviews remain excluded; the underlying
database is not fully curated. A completed evaluation does not make the MIC
cutoff a validated clinical breakpoint or the repository a complete release.
See [final results](final_results.md), [reproduction](reproduce.md) and the
[README](../README.md).

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
- A **sequence family group** is a connected component under the fixed 80%
  global edit-similarity rule in Section 12. It keeps related sequences and
  amidation variants together across primary training/test boundaries.

Laboratory/reference-only and clinical-isolate analyses are secondary. Their
candidate tables are not automatically evaluation partitions. The current
clinical-only candidates were explicitly reserved by the split procedure.

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
prove the peptide has never been tested against another strain. Cohort tables
and evaluation partitions are distinct: the split reserves the 17 current
clinical-only candidates and their complete 80%-similarity families. These
reserved examples must not subsequently enter primary training.

## 10. Regression is a separate, provisional target

Only eligible single numeric MICs contribute to median log10 MIC summaries
within an assay group. Ranges, censored values and `±` central values do not
become point-value regression observations. Preserve the number of numeric
measurements and their spread where available.

Stored regression summaries are not a finalized regression training dataset.
Their assay notes and cohort definition still require review. A classifier
exclusion must not be assumed to automatically remove every corresponding
regression summary already stored in an audit table.

## 11. Confirmed cohort and recorded run

The full manifest selects 2,217 validated records. Of these, 1,951 contain at
least one locally matching E. coli MIC (88.0%). Sequence plus E. coli MIC checks
retain 1,420 records; chemistry checks retain 1,222. There are 3,955 E. coli MIC
rows before the modeling filters. These are not counts of classifier examples.

The broad audit contains 776 sequence-plus-amidation identities:

| Identity status | Count |
|---|---:|
| Active | 271 |
| Inactive | 218 |
| Needs note review | 232 |
| Mixed evidence | 43 |
| Ambiguous | 12 |

Active and inactive subsets form 489 ready identities across 486 sequences.
The other 287 identities remain outside the classifier. Status categories are
mutually exclusive; secondary flags may overlap. These counts supersede the
earlier 391-example dataset.

Recorded inputs and outputs:

- Full manifest: `data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json`.
- Label output: `data/processed/ecoli-full-v1-3cec454243584936/45c31cc2e742/`.
- Provenance: the label output's `provenance/83780b72da84/` directory.
- Split: the provenance output's `splits/grouped-v1-dbe608485cce9a7e/` directory.
- Baseline: `results/baselines/physchem-logistic-v1-b8379a4f32704165/`.
- Comparison: `results/comparisons/feature-comparison-v1-3c0530a63f1c77d0/`.

The provenance builder applies recorded reviews to measurement evidence and
recomputes cohort labels. An older label output or classifier CSV must not be
assumed to contain the same ready identities. The confirmed provenance
`broad_candidates.csv` and the fixed split's `train.csv` define this experiment's
eligible cohort and training examples.

Reports retain input hashes, specifications, dependency versions and preparation
timestamps. Directory names above identify the confirmed local run, not names
to impose on a different rebuild.

## 12. Fixed sequence-family grouping and partitions

`prepare_splits.py` validates the broad candidate, broad audit and clinical-only
candidate tables. Ready broad identities must equal the audit's eligible set;
clinical-only labels must agree with their broad labels. Modeled identities
must be unique.

### Similarity and connected components

Use unit-cost global Levenshtein distance: each substitution, insertion and
deletion costs one.

```text
similarity = 1 - edit_distance / max(sequence_lengths)
```

Join pairs with similarity at least 0.8 using a rational threshold comparison.
Form connected components over all 772 unique broad-audit sequences, including
unready identities. Ignore amidation for grouping, while retaining it in model
inputs. Component IDs hash their sorted sequence membership. This produces 481
audit groups; the largest contains 44 sequences.

Activity labels and model results do not determine the graph or partition
allocation. Connected components can join distant sequences through
intermediates. This rule is not standard alignment identity or proof of
homology and does not isolate every shared motif or study effect.

### Allocation

1. Reserve every family containing a current clinical-only candidate.
2. Put its clinical-only identities in `clinical_exploratory.csv`; place the
   other broad identities in those families in `clinical_family_buffer.csv`.
   Do not train on the buffer.
3. For remaining groups, hash `str(seed) + ':' + group_id` with SHA-256, using
   seed 42. A group enters the test set when its hash-derived value is below
   the fixed 0.2 threshold; otherwise it enters training.

The 20% setting applies per family and does not guarantee 20% of examples.
Do not search seeds or alter the split based on model performance.

| Partition | Identities | Unique sequences | Groups | Active | Inactive |
|---|---:|---:|---:|---:|---:|
| Training | 394 | 392 | 242 | 217 | 177 |
| Primary test | 71 | 71 | 47 | 34 | 37 |
| Clinical exploratory | 17 | 17 | 14 | 14 | 3 |
| Clinical family buffer | 7 | 7 | 2 | 6 | 1 |

The buffer's two groups are contained in the clinical reservation's 14 groups;
these group counts are not additive. Different amidation states of one
sequence can appear in the clinical and buffer files, so their unique-sequence
counts are not additive either. Both files remain part of the same reservation.

Each broad identity is assigned once. Training, primary test and the entire
clinical reservation are disjoint by 80%-similarity family. Both primary
training and test contain both classes. Checks confirm all qualifying sequence
pairs share their component and partition role. These checks are conditional
on the declared rule and valid input metadata.

The broad training set can contain other eligible clinical-origin examples.
This is not a pure reference-to-clinical transfer experiment. The reserved
clinical set has only 17 identities in 14 families and three inactive examples;
its completed evaluation is exploratory.

## 13. Training-only baseline and feature comparison

The modeling scripts verify training bytes against the split manifest and
check identity/group metadata. They do not open primary test, clinical
exploratory or clinical-buffer CSVs. Split/group metadata can be read, but
held-out outcomes must not guide model choices or preprocessing.

### Cross-validation and selection

- Use five-fold `StratifiedGroupKFold`, shuffling with seed 42, within training.
- Verify no family crosses a fold boundary, every training row is validated
  once, and both classes occur in every fitting and validation subset.
- Save `cv_folds.csv`. Before reuse, verify its recorded hash and agreement
  with the training identities and labels. All primary feature comparisons
  use these same saved folds; do not regenerate them to improve scores.
- Fit `StandardScaler` inside each fitting fold's pipeline, then transform
  that fold's validation rows. Do not scale all rows before cross-validation.
- Use L2 logistic regression with `lbfgs`, `max_iter=5000`, no class weighting
  and equal weights per modeled identity. Assay counts are not sample weights.
  Treat convergence warnings as failures.
- Select C by maximum mean fold balanced accuracy; exact ties favor smaller C.
  The initial grid was `0.01, 0.1, 1, 10`. The comparison extends it to
  `0.0001, 0.001, 0.01, 0.1, 1, 10` for every logistic feature set.
- Keep the prediction threshold at 0.5. This is distinct from the 32 µM MIC
  threshold used to construct the activity labels.

After selecting C, fit each feature model on all 394 training identities and
save it locally. Refitting is not independent evaluation. Corrected charge
plus amidation (C=0.1, prediction threshold 0.5) was subsequently declared
primary before test evaluation; see Section 15.

### Representations

| Model | Features |
|---|---|
| Majority dummy | Most frequent class in each fitting subset |
| Residue charge counts | `count(K,R) - count(D,E)` |
| Corrected charge proxy | Residue charge counts plus 1 for amidation |
| Corrected charge + amidation | Corrected proxy and a separate 0/1 amidation flag |
| Composition | 20 amino-acid fractions, length and amidation: 22 features |
| Physicochemical | Length; K/R, D/E, H, hydrophobic, aromatic, G and P fractions; amidation; corrected charge: 10 features |

The hydrophobic set is `AVILMFWY`; the aromatic set is `FWY`. Charge is a
count-based proxy, not Biopython's `charge_at_pH` or a complete pH-specific
calculation. The amidation correction was already present in the original
ten-feature baseline; do not apply a second +1.

Dropping the separate amidation flag does not remove the amidation information
already encoded by corrected charge. Composition and physicochemical features
are alternatives, not nested additions. Correlation, scaling and regularization
affect comparisons; they do not establish causal chemical effects. Neither
representation captures residue order.

### Recorded development results

| Model | Selected C | Mean fold balanced accuracy | Mean fold AUROC |
|---|---:|---:|---:|
| Majority dummy | — | 0.500 | 0.500 |
| Residue charge counts | 0.1 | 0.684 | 0.722 |
| Corrected charge proxy | 1 | 0.713 | 0.743 |
| Corrected charge + amidation | 0.1 | 0.739 | 0.752 |
| Composition | 0.1 | 0.731 | 0.810 |
| Physicochemical | 0.01 | 0.735 | 0.755 |

Scores reuse validation evidence for selection and are not final test estimates.
Composition's higher AUROC does not establish a reliable balanced-accuracy
improvement. The compact, original and composition models remain comparisons,
not declared final winners.

## 14. Bootstrap intervals and grouping sensitivity

### Paired family bootstrap

For each feature set, retain the selected C's out-of-fold probabilities and
predicted labels, pooled across folds. Use 2,000 draws with seed 42, sampling
the 242 original training families with replacement. Keep every selected
family's identities together with their multiplicity. Share draws across all
models and paired differences.

Compute percentile 95% intervals for balanced accuracy, AUROC, average precision
and differences between adjacent feature-ladder entries. Omit one-class draws
consistently across models, report their count and reject excessive omission.
The confirmed run retained all 2,000 draws. Adjacent differences are directional
as recorded in `paired_differences.csv`; this is not an all-pairs comparison or
a multiple-testing-adjusted significance analysis.

Mean fold scores select C. Bootstrap intervals target pooled predictions, which
have different weighting. Do not attach a pooled interval to a mean fold score
as though they describe the same summary.

Intervals condition on selected C and the already-fitted fold predictions.
They do not refit models, repeat selection or fully capture dependence from
overlapping fold-training sets. Fold standard deviations are descriptive, not
confidence intervals. Neither output establishes final test performance or
equivalence between models.

### Sequence-group sensitivity

Keep the original 394 training identities. Regroup all 772 unlabeled audit
sequences at prespecified 60%, 70% and 80% global edit similarity. At 60% and
70%, construct training-only five-fold grouped validation with seed 42. At 80%,
reuse the saved original folds. Within each threshold the dummy and
physicochemical models share folds; logistic regression uses the extended C
grid. If too few groups or a one-class fold prevents evaluation, report that
failure instead of searching seeds. Save alternative fold assignments and
group mappings.

| Similarity | Audit groups | Largest audit group | Physchem balanced accuracy | AUROC |
|---|---:|---:|---:|---:|
| 60% | 334 | 79/772 sequences | 0.730 | 0.762 |
| 70% | 411 | 47/772 | 0.732 | 0.765 |
| 80% | 481 | 44/772 | 0.735 | 0.755 |

Reports retain the complete cluster-size histograms, singleton counts, median
and 90th-percentile group size. Audit sizes count unique sequences; training
group sizes can differ because amidation variants are separate identities.

Broader grouping connects 148 training identities at 60%, and 68 at 70%, to
groups also represented in the frozen holdouts. This diagnostic reads sequence
and partition metadata, not held-out outcomes. It does not move examples or
open held-out CSVs. Similar training CV scores do not prove the original test
is independent under the broader thresholds. Do not choose the grouping rule
from the best-looking validation score.

## 15. Frozen model, sensitivity and completed evaluation

Corrected charge plus amidation was declared primary after development but
before opening held-out outcomes: L2 logistic regression C=0.1, threshold 0.5,
primary metric balanced accuracy on the original 71-identity test. An exact
0.5 prediction tie is inactive. All other compared models are secondary.
The pretest protocol records artifact hashes, package versions, local-alignment
settings, subset membership, ranking budgets and bootstrap support rules.
It was pushed before evaluation; this is not preregistration before development.

Training-only cutoff sensitivity at 16/32/64 µM exactly replayed the original
394 labels at 32 µM. At 16 µM, 377 remained eligible (170 active, 207 inactive);
at 64 µM, 380 remained eligible (242 active, 138 inactive). New eligibility was
also examined among all 639 training-assigned audit identities, without adding
examples to the frozen training set or fitting alternative-cutoff models.
The primary cutoff stays 32 µM.

Final overall test results: balanced accuracy 0.717 (95% family-bootstrap
interval 0.593–0.874), AUROC 0.751 (0.592–0.888), TP=23, FN=11, TN=28, FP=9.
The 17-identity clinical set remains exploratory; the family buffer was not
used for fitting or evaluation. Novelty subsets and all secondary models were
reported under the same fixed protocol. See [final results](final_results.md)
for the full context, intervals, ranking diagnostics and limitations.

Do not alter model choices, labels, cutoffs or partitions to improve these
reported held-out scores. Future development needs a separate version and a
new independent evaluation resource; the current outcomes are now known.
Release packaging and unresolved curation remain separate unfinished tasks.

## 16. Documentation, dates, committed artifacts and data terms

Retain raw measurements, review reasons, source identifiers, original manifests
and build reports. Commit source, tests, dependency specifications, documentation,
review decisions and selected small result summaries. The fixed split's
`split_manifest.json`, `sequence_groups.csv` and `split_assignments.csv` are
explicitly tracked exceptions under the ignored processed-data tree. Baseline
`cv_folds.csv` and run reports support the shared-fold comparison.

Raw records/listing caches, cleaned datasets, feature/prediction exports and
fitted `.joblib` artifacts are generated locally. Reports may hash generated
files that are not committed; a fresh clone does not provide every run file.
Data-recovery instructions, complete data-tool dependency pins and retrieval
dates remain release tasks. Do not claim full clean-checkout reproduction yet.

Transcribe actual access dates from manifests, cache envelopes and retrieval
logs. Reused records may predate the remaining downloads. Documentation and
terms-check dates are not dataset retrieval dates.

This project acknowledges [DBAASP](https://dbaasp.org/) and the original
experimental publications linked in its records. Consult its
[terms](https://dbaasp.org/terms-and-conditions) when using or sharing source or
derived data. The earlier project terms check was dated 2026-10-18; that date
is not a new terms check or retrieval date for this document.

The requested publication is Pirtskhalava et al., *DBAASP v3*, 2021,
[doi:10.1093/nar/gkaa991](https://doi.org/10.1093/nar/gkaa991).
Any future project code license is separate from the source data's terms.
