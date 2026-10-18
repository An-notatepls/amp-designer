# AMP Designer: E. coli peptide activity predictor

An undergraduate bioinformatics project predicting recorded antimicrobial
activity against *Escherichia coli* from peptide sequence and chemical features.
The first task is classification using minimum inhibitory concentration (MIC).
Candidate generation is planned; it has not been implemented. MIC measures
growth inhibition, not necessarily bacterial killing.

## Status

**Policy v0.1 — draft. Documentation updated: 2026-10-18.**

Data collection, auditing and label-building code are implemented. Additional
note-review decisions have been prepared; their local rebuild and the
provenance-specific cohorts still need confirmation. No model has been trained,
no evaluation results are available, and no test set or v1 policy is locked.

## Confirmed data counts

These are the last confirmed full-cohort results reported on 2026-10-13.
Counts are not interchangeable: records, measurements and peptide identities
represent different things.

| Stage | Count |
|---|---:|
| Downloaded and validated peptide records | 2,217 |
| Records with at least one locally matching E. coli MIC | 1,951 |
| Records passing the sequence and E. coli MIC checks | 1,420 |
| Records also passing the peptide chemistry rules | 1,222 |
| E. coli MIC measurements before modeling filters | 3,955 |
| Classifier-ready peptide identities | 391 |
| Unique sequences in that classifier table | 390 |
| Active / inactive classifier examples | 190 / 201 |

The local MIC-filter yield is 88.0% of the downloaded records. The saved
server-side search is a filtered ID collection, not the whole DBAASP database.
Record contents were retrieved at different times.

Applying the additional reviewed notes to the uploaded assay tables predicts
489 classifier examples across 486 sequences. Those counts remain **expected,
not locally confirmed**, until the corresponding rebuild is checked.

## Primary training cohort

**The first model will use the reviewed broad E. coli classifier dataset.**
It includes eligible laboratory/reference assays, accepted clinical-isolate
assays and assays whose strain origin is unspecified. It is not a
reference-strain-only dataset. The broad model's candidate file is
`ecoli_classifier_dataset.csv` from `build_labels.py`, or the corresponding
`broad_candidates.csv` from `build_provenance.py` with the same input policy.

Each example represents a **sequence plus C-terminal amidation state**. The
model must represent amidation explicitly. Different amidation states of the
same sequence are kept together during evaluation splitting.

Laboratory/reference-only models and transfer to clinical isolates are
secondary experiments, subject to cohort size and positive strain-origin
annotations. Clinical-isolate measurements are not inherently invalid. An
empty note does not establish reference-strain origin.

## How labels are built

The initial peptide scope is canonical sequences of 10–50 amino acids,
annotated as monomers, with no recorded unusual residues or additional bonds.
C-terminal amidation is allowed; other recorded terminal modifications are
excluded. The length window is a practical scope choice, not a biological
definition of an antimicrobial peptide.

MIC values are converted to micromolar (µM) where supported. At the provisional
32 µM threshold, a definite MIC at or below 32 is active and a definite MIC
above 32 is inactive. Ranges and inequalities keep their bounds. Unclear or
contradictory evidence is flagged rather than assigned a label by majority vote.

### Worked example — synthetic, not an actual DBAASP record

Suppose one otherwise eligible peptide has three measurements with the same
strain, medium and other recorded conditions. Its notes are empty or accepted
after review. All concentrations are already in µM.

| Reported MIC | Interpretation at 32 µM |
|---|---|
| `8` | Active: 8 ≤ 32 |
| `8–16` | Active: the entire reported range is ≤32 |
| `≤16` | Active: the upper bound is ≤32 |

All three support the active label, so this peptide identity is included as
active. The range is never replaced with a midpoint. If the last measurement
were `>64`, the evidence would conflict and the identity would be excluded from
the first classifier. A range such as `16–64` would instead leave uncertainty.

See [the complete data policy](docs/data_policy.md) for chemistry assumptions,
units, grouping, note review and provenance-specific rules.

## Baselines and evaluation plan

The first comparison will be a majority-class dummy predictor versus
regularized logistic regression using understandable features: length,
amino-acid composition, charge/hydrophobicity summaries and amidation.
Coefficients describe associations; they do not establish biological causes.
No claim is made about later model architectures.

- Group similar sequences before splitting. Exact sequences and their chemical
  variants must not cross training and evaluation partitions.
- An 80% normalized edit-similarity rule has been explored as a grouping
  diagnostic. It is not yet the frozen split definition or proof of homology.
- Use group-aware cross-validation within training data for model selection.
- Keep a held-out test set out of preprocessing fits and model tuning.
- Report class counts, balanced accuracy, precision/recall and MCC alongside
  the dummy baseline, with uncertainty where feasible.
- Rebuild labels at 16 and 64 µM as predefined sensitivity analyses around
  the primary 32 µM task. Keep sequence-group assignments fixed, report changes
  in eligible examples, and do not choose the cutoff by test performance.
- Consider reference-to-clinical-isolate transfer only if independent groups
  support it. Reserved clinical test examples and their close relatives must
  not enter either training variant.

No split is locked yet. Dataset fingerprints, grouping parameters, seed and
partition counts will be recorded before evaluation.

## Limitations

Published MIC measurements vary across strains, laboratories and conditions.
Unknown conditions remain unknown, even when recorded fields match. Empty
chemistry annotations do not prove an absence of modifications.

The initial representation excludes lowercase/noncanonical sequences and
recorded unusual residues, including D-amino-acid representations caught by
these rules. It cannot establish stereochemistry independently. Results will
not cover all modified, cyclic, multichain or longer peptides.

Excluding mixed and unresolved examples can favor easier cases. Repeated or
closely related measurements are not necessarily independent. The dataset is
not a representative sample of all possible peptides or E. coli isolates.
Predictions will not establish clinical efficacy or experimental activity.

## Reproduce the current data build

Run from the project folder. Initial setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install certifi biopython
```

On later sessions, activate the existing environment rather than recreating it.
Dependency versions must be recorded for the eventual reproducible release.

For the existing downloaded full cohort:

```bash
python audit_pilot.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json"
python build_labels.py "data/raw/full/ecoli-full-v1-3cec454243584936/manifest.json" --threshold 32 --note-reviews data/curation/full_note_reviews.csv
```

Despite its name, `audit_pilot.py` audits the records selected by whichever
manifest is supplied. `download_pilot.py` similarly downloads or resumes any
compatible manifest-defined collection.

The label builder prints its output folder. Pass that folder as the first
argument to `build_provenance.py`, together with the same review file and
threshold. Use `python build_provenance.py --help` for its argument syntax.

These commands assume the local manifest, downloaded records and review file
exist. A newly queried listing is not an exact reproduction of the saved
collection. A clean-checkout reproduction recipe, archived manifest metadata
and pinned dependencies remain release tasks.

### Files committed versus generated

The intended repository convention is:

- **Commit:** source code, tests, documentation, dependency specifications,
  reviewed decision CSVs in `data/curation/`, and the frozen ID-manifest metadata
  and split definitions once prepared for release.
- **Generate locally:** downloaded records and listing caches under `data/raw/`,
  cleaned tables under `data/processed/`, and trained model artifacts.
- **Do not commit:** the virtual environment, credentials or temporary files.

The raw and processed data directories are excluded by the project's ignore
rules. Frozen manifest metadata must therefore be explicitly exported to a
tracked location before claiming reproducibility from a fresh checkout.

## Data source, citation and license

Data source: [DBAASP](https://dbaasp.org/). Dataset retrieval dates must be
transcribed from saved manifest/cache timestamps and retrieval logs before
release. They must not be replaced with the documentation-update date.

DBAASP's [terms and conditions](https://dbaasp.org/terms-and-conditions), checked
2026-10-18, permit access, adaptation and redistribution while requiring public
acknowledgement of DBAASP as the source. This project acknowledges DBAASP and
the original experimental publications linked in its records. Changes made by
this project's cleaning pipeline should be described when sharing derived data.

Please cite:

Pirtskhalava M, Amstrong AA, Grigolava M, et al. **DBAASP v3: database of
antimicrobial/cytotoxic activity and structure of peptides as a resource for
development of new therapeutics.** *Nucleic Acids Research*. 2021;49(D1):D288–D297.
[doi:10.1093/nar/gkaa991](https://doi.org/10.1093/nar/gkaa991).

A project code license has not yet been selected. DBAASP's data terms are
separate from any future license for this repository's original code.
