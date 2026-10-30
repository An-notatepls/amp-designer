# AMP Designer: E. coli peptide activity predictor

An undergraduate bioinformatics project predicting recorded *Escherichia coli*
peptide activity from sequence and C-terminal amidation. The current product
is a binary MIC activity predictor. Candidate generation has not been implemented.
MIC measures growth inhibition, not necessarily bacterial killing.

Read the [project walkthrough](notebooks/project_walkthrough.ipynb) for the
data-to-results story, figures and limitations. Saved outputs are included;
rerunning it reads bundled summaries without training or repeating evaluation.

- The declared primary model achieved **test AUROC 0.751 (95% family-bootstrap
  interval 0.592–0.888)** and **balanced accuracy 0.717 (0.593–0.874)**.
- Simple charge features explain substantial predictive signal. These experiments
  do not establish a richer feature set as superior.
- The fixed test and exploratory clinical sets have now been evaluated. No model
  or cutoff was selected from their results.

## Status

Documentation updated: **2026-10-28**. Data policy v0.1 is frozen for the reported
experiment; unresolved reviews remain excluded. Collection, labeling, grouped
validation, final evaluation and an inference command are implemented. The
inference command checks the frozen model fingerprint and does not retrain it.
The walkthrough is implemented with saved outputs and six figures. Tests and
training-only modeling reproduction are supported by committed, fingerprint-
verified inputs. The original-model export, original raw ID manifest and observed
collection-date metadata still require the local export described below.
Normal-kernel notebook Run All remains to be checked. Hosted CI on Python
3.13.4 passed all 63 tests, fingerprint verification and the full training-only
numerical reproduction. A complete live DBAASP rebuild has not been demonstrated.

## Data and primary cohort

The saved DBAASP search contains 2,217 downloaded records, of which 1,951 have
an E. coli MIC. Chemistry and sequence filtering retain 1,222 records containing
3,955 E. coli MIC rows before assay-level modeling filters. After review and
aggregation, the **reviewed broad cohort** has 489 eligible sequence-plus-amidation
identities across 486 sequences: 271 active and 218 inactive.

The broad cohort includes accepted laboratory/reference, clinical-isolate and
unspecified-origin assays. Empty notes do not establish reference-strain origin.
The primary label is recorded MIC **≤32 µM**, a project cutoff rather than a
validated clinical breakpoint. Eligible peptides are annotated monomers of
10–50 canonical residues, with no recorded additional chemistry except permitted
C-terminal amidation. Ranges and inequalities retain their bounds; conflicting
or unresolved identities are excluded rather than labeled by majority vote.

For example, consistent measurements `8`, `8–16` and `≤16` µM support active at
32 µM. Adding `>64` creates opposing evidence and excludes that identity.
See [data policy](docs/data_policy.md) for the exact assay and chemistry rules.

| Fixed partition | Identities | Families | Active / inactive |
|---|---:|---:|---:|
| Training | 394 | 242 | 217 / 177 |
| Primary test | 71 | 47 | 34 / 37 |
| Clinical exploratory | 17 | 14 | 14 / 3 |
| Clinical family buffer | 7 | 2 | 6 / 1 |

Sequence-family connected components use global edit similarity ≥80%; identical
sequences and amidation variants stay together. Clinical families are reserved
before allocation. The buffer's two families are within the clinical reservation.
This heuristic does not prove homology or remove every shared motif/study effect.

## Final findings

The primary model is standardized L2 logistic regression using the corrected
charge proxy and an amidation flag, with C=0.1 and prediction threshold 0.5.
It was declared after development comparisons and before held-out evaluation.

| Model | Test AUROC (95% interval) | Test balanced accuracy |
|---|---:|---:|
| Majority dummy | 0.500 (0.500–0.500) | 0.500 |
| Uncorrected charge | 0.767 (0.611–0.880) | 0.760 |
| Corrected charge | 0.769 (0.620–0.889) | 0.775 |
| **Corrected charge + amidation — primary** | **0.751 (0.592–0.888)** | **0.717** |
| Composition | 0.780 (0.659–0.898) | 0.731 |
| Physicochemical features | 0.746 (0.589–0.878) | 0.702 |

The primary model correctly identified 23/34 active and 28/37 inactive identities.
Precision was 0.719, sensitivity 0.676 and specificity 0.757. Its test precision
at top 10 was 0.80 (95% interval 0.50–1.00), compared with 47.9% active prevalence
in the test pool. These are retrospective database results, not prospective
wet-lab success rates. Secondary models remain secondary despite their scores.

On the combined strict novelty subset, AUROC was 0.814 (0.659–1.000; n=32).
The clinical analysis had only three inactive identities and cannot establish
reliable clinical transfer. See [complete results and caveats](docs/final_results.md)
for subset scores, development intervals, ranking results and label sensitivity.

## Quick start

Use Python **3.13.4** for the original interpreter; `.python-version` records it.
From a fresh clone:

```bash
python3 -m venv .venv
./.venv/bin/python -m pip install -r requirements.txt
./.venv/bin/python tools/verify_release.py
./.venv/bin/python tools/reproduce_modeling.py
```

The last command reproduces training-only baseline/feature comparisons in a
new output directory and compares numerical results; it does not repeat final
evaluation or change original reports. Locally, all 63 tests and the full
numeric reproduction passed on Python 3.12.14 with the recorded numerical-library
versions. The full check also passed in [hosted CI on Python 3.13.4](https://github.com/An-notatepls/amp-designer/actions/runs/40631386045).

Read the [reproduction guide](docs/reproduce.md) for the three reproducibility
levels, dependency roles, tests and the separate notebook environment.

The original evaluated model remains on the original Mac. Its verified export
produces a portable JSON model; no replacement model is invented. After installing
that export, inference needs only Python's standard library:

```bash
python3 modeling/predict_portable.py --sequence AAAAAAAAAA --unmodified
```

The synthetic sequence demonstrates syntax only; use `--amidated` when appropriate.
Predictions cover only canonical 10–50-residue monomers with supported chemistry;
the probability is an uncalibrated model score. Existing local joblib inference
continues to work with `modeling/predict.py` in the frozen modeling environment.

## Limitations and source

The test is small and cannot reliably rank modest model differences. Published
assays differ across strains, laboratories and conditions; empty annotations
are not verified absence of chemistry. Excluding unresolved/mixed cases can
favor easier examples. Of 232 pending-review identities, 180 have provisional
active-only and 45 inactive-only evidence; the excluded queue leans active.
They are not accepted labels. Predictions do not establish toxicity, selectivity,
clinical efficacy or activity against arbitrary isolates.

Source: [DBAASP](https://dbaasp.org/) and the original publications linked in its
records. Please cite Pirtskhalava et al., *DBAASP v3*, Nucleic Acids Research
49(D1):D288–D297 (2021), [doi:10.1093/nar/gkaa991](https://doi.org/10.1093/nar/gkaa991).
Consult [DBAASP terms](https://dbaasp.org/terms-and-conditions) for source/derived
data. Original code is [MIT licensed](LICENSE). Derived data retain DBAASP attribution
and terms; see [data attribution](data/README.md). Observed collection metadata
will be exported from the original cache. Per-record retrieval dates were not
recorded; this documentation date is not a dataset access date.
