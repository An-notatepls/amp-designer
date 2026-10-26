# Pre-test evaluation protocol

Prepared at UTC: 2026-10-26T02:08:55.843179+00:00

This is a development-informed decision recorded before held-out outcomes.
It is not a claim of preregistration before development. Commit this file and
`evaluation_manifest.json` before opening any held-out outcome table.

## Primary model and analysis

- Model: corrected charge plus C-terminal amidation (`charge_corrected_amidation`).
- Charge: count(K,R) - count(D,E) + amidation; features ordered as charge, amidation.
- L2 logistic regression: C=0.1, prediction threshold 0.5 (exact tie inactive).
- Primary metric: balanced accuracy on the original overall primary test set.
- Secondary: AUROC, average precision, MCC, accuracy, precision, sensitivity,
  specificity, confusion counts, and tie-aware precision@5/10/20.
- Original test: 71 identities, 47 families;
  class counts are recorded in the existing split manifest. No test predictions
  or outcome-level analysis were generated here.
- Evaluate the five other fixed feature-ladder models once as secondary
  comparisons. Use their already-selected Cs; select nothing from test scores.
- Keep the clinical set exploratory and the clinical-family buffer out of training.

## Fixed artifacts

- Primary model SHA-256: `5573f24eb93693b207b8e356c7e842433c1f9f26b08fd1fea6039dc85e31680e`
- Source model bundle SHA-256: `115fdea3c1b8dd6aeb15eb283fa01a6c9f5a4dd2d22205684a69af1ecb2ccbac`
- Split manifest SHA-256: `5bc08d219cb6a86f801eb8d918bc9b6f25ce9875e7d627d31487093567a4e3a3`
- Comparison report SHA-256: `9cd0ae0287a29117db3b7c561944877b32ad47d545449f939a86827da44b7024`
- Training table SHA-256: `9e5461de996dc8bbb34f4bf36ac069b0a8a24ad2f8ecb3f24a5721734a49ffa1`
- Complete settings, source-code hashes and dependency versions: `evaluation_manifest.json`.

## Sequence-only novelty subsets

All test identities remain in the primary overall analysis. The following are
secondary subsets, defined without reading held-out outcome CSVs:

| Subset | Rule | Identities |
|---|---|---:|
| global60_strict | Maximum direct global edit similarity to training <0.60 | 45 |
| local60_strict | No detected qualifying local identity >=0.60; complete optimal-tie scans | 32 |
| combined60_strict | Intersection of the preceding subsets | 32 |

Local alignment is Smith-Waterman with BLOSUM62, gap open -10 and extension -0.5.
Identity divides exact residue matches by all alignment columns including gaps.
Qualification requires at least 10 aligned residue pairs and coverage >=80%
of the shorter peptide (paired residues divided by shorter length). Both peptide
coverages and aligned lengths are exported. A short fragment can therefore
qualify despite low coverage of its longer neighbor.

For each pair, scan up to 100 score-optimal alignments. Any incomplete tie scan
excludes that test identity from local/combined strict subsets. Maxima refer to
these score-optimal local alignments, not every possible suboptimal alignment.
Absence of a detected neighbor does not establish biological independence.

Keep the original split fixed. The 60% direct-neighbor screen differs from
transitive 60% connected-component linkage and from alignment percent identity.

## Uncertainty, small subsets and interpretation

Final overall intervals will use 5,000 paired whole-family bootstrap draws,
seed 42, percentile 95% intervals on fixed predictions. Disclose one-class
resample omissions; suppress an interval if more than 10% are omitted.
For AUROC intervals additionally require >=10 families and >=5 examples per
class; otherwise report the support limitation. AUROC is undefined with one
class. Always report subset n, class counts and family counts; never hide an
empty or small subset. The clinical set has only three inactive identities.

The 71-identity primary test provides a limited final check, not a powered
ranking of small feature-set differences. No fixed interval width is promised.
Report all declared outcomes and retain a disappointing result without retuning.

## OOF prioritization diagnostics

Precision@5/10/20 budgets are illustrative and fixed before inspecting P@k;
they are not established laboratory budgets. Round scores to 12 decimal places
and average over random boundary-tie selection. Family-diverse ranking chooses
one highest-scoring identity per original family, breaking within-family ties
by sequence+amidation lexicographic order without labels.

OOF intervals use paired whole-family resampling of fixed predictions and
are conditional development diagnostics. They omit refitting/selection
uncertainty. A resample can contain replicated family copies; it is not an
actual recommendation to synthesize duplicates. Pooled OOF scores come from
different fitted models and are not guaranteed to be mutually calibrated;
fold-specific P@k is also exported. Database P@k does not establish prospective
success for generated peptides or a different candidate population.

## Release conditions

Review the metadata-only diagnostics, then commit this protocol, the membership
table and manifest before test outcomes are opened. Any change to model,
metric, alignment or subset rules requires a dated amendment before outcome
inspection. Model binaries remain generated local artifacts; their hashes
are committed. This script does not perform final test evaluation.

Derived from DBAASP data, https://dbaasp.org/; see project citation and data terms.
