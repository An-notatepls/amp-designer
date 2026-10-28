# Final evaluation and development results

Documentation updated 2026-10-27. Final results were reported by the evaluator
and checked against the exported metrics, intervals and ranking summaries.
This document does not replace the original run reports or their fingerprints.

## Declared primary analysis

Primary model: corrected charge plus amidation; standardized L2 logistic
regression, C=0.1, prediction threshold 0.5 (an exact tie predicts inactive).
Features are `count(K,R) - count(D,E) + int(c_amidated)` and `int(c_amidated)`.
Charge is a count proxy, not a full pH-specific charge calculation.
The target is recorded E. coli MIC ≤32 µM. Models were fit on 394 training
identities only. Primary metric: balanced accuracy on the original overall test.

Choices were recorded in the development-informed pretest protocol, not
preregistered before all development. The protocol was pushed before evaluation
(commit b884ffe); the evaluator was pushed in 387a4f5. Other features, subsets,
ranking measures and clinical results are secondary. No new feature set,
regularization strength, classification cutoff or split is selected from them.

## Primary model by evaluation panel

Intervals below are 95% percentile whole-family bootstrap intervals from 5,000
shared draws, conditional on fixed models/predictions. They do not include
training-set selection or model-fitting uncertainty.

| Panel | Identities | Active / inactive | Families | AUROC (95% interval) | Balanced accuracy (95% interval) |
|---|---:|---:|---:|---:|---:|
| overall | 71 | 34 / 37 | 47 | 0.751 (0.592–0.888) | 0.717 (0.593–0.874) |
| global60_strict | 45 | 24 / 21 | 29 | 0.729 (0.559–0.893) | 0.646 (0.519–0.844) |
| local60_strict | 32 | 17 / 15 | 19 | 0.814 (0.659–1.000) | 0.686 (0.565–0.929) |
| combined60_strict | 32 | 17 / 15 | 19 | 0.814 (0.659–1.000) | 0.686 (0.565–0.929) |
| clinical_exploratory | 17 | 14 / 3 | 14 | 0.762 (interval suppressed) | 0.595 (0.400–0.929) |

The overall confusion counts are TP=23, FN=11, TN=28 and FP=9. Accuracy=0.718,
precision=0.719, sensitivity=0.676, specificity=0.757, MCC=0.435 and average
precision=0.721. The AUROC interval is wide; 71 examples in 47 families are a
limited final check rather than a powered comparison of small model differences.

The local and combined strict panels contain the same 32 identities. Their
matching scores are not independent confirmations. A higher subset AUROC does
not demonstrate superior generalization; these panels have different membership,
class balance and small samples.

Clinical results are exploratory: TP=12, FN=2, TN=1, FP=2. Specificity is only
1/3. There are too few inactive examples for the declared AUROC interval rule
(minimum five examples/class and ten groups). Of 5,000 bootstrap draws, 4,835
contained both classes. The AUROC interval is suppressed, not zero. Clinical
average precision of 0.930 should be interpreted against 14/17 active prevalence.
This is not a pure reference-to-clinical transfer experiment: the broad training
pool can contain other eligible clinical-origin evidence.

## Sequence-neighbor diagnostic

Grouping used global edit similarity ≥0.8 across all 772 broad-audit sequences,
including unready identities; labels did not define groups or allocation.
Global similarity is `1 - Levenshtein_distance / max(lengths)`. Connected
components can chain sequences and can miss truncation relationships.

The pretest local diagnostic used score-optimal Smith–Waterman alignments,
BLOSUM62, gap-open −10 and gap-extend −0.5. Identity includes internal gap
columns in its denominator. A qualifying match needs at least ten aligned
residue pairs and coverage of at least 80% of the shorter peptide. Up to 100
score-optimal alignment ties were scanned per pair; incomplete scans would
exclude an identity from the strict local subset (none occurred).

Global-strict requires maximum direct global similarity <0.6. Local-strict
requires no detected qualifying score-optimal local identity ≥0.6. Combined is
their intersection. These are diagnostics under declared scoring rules, not
proof of homology, complete local-match searches or biological independence.
There were 45 global-strict and 32 local/combined-strict test identities.

## Development comparison

All models shared saved five-fold group-aware cross-validation assignments.
Scaling was fitted inside each training fold. C was selected by mean fold
balanced accuracy from 0.0001, 0.001, 0.01, 0.1, 1 and 10. These pooled OOF
AUROCs differ from averages of per-fold AUROC. Intervals use 2,000 shared family
resamples and are conditional on selected out-of-fold predictions; they do not
remove hyperparameter-selection effects or fully account for overlapping fits.

| Model | Selected C | Mean fold balanced accuracy | Pooled OOF AUROC (95% interval) |
|---|---:|---:|---:|
| Dummy | — | 0.500 | 0.500 (0.500–0.500) |
| Uncorrected charge | 0.1 | 0.684 | 0.730 (0.657–0.803) |
| Corrected charge | 1 | 0.713 | 0.747 (0.675–0.813) |
| Corrected charge + amidation | 0.1 | 0.739 | 0.750 (0.676–0.817) |
| Composition | 0.1 | 0.731 | 0.794 (0.739–0.853) |
| Physicochemical | 0.01 | 0.735 | 0.758 (0.684–0.821) |

Composition minus the primary representation had a paired development AUROC
difference of +0.044 (95% interval −0.016 to +0.113). This did not establish
superiority or equivalence. The simpler two-feature model was declared primary
before final testing. All six final test results remain reported in the README
and original exports; a secondary model's higher test score does not change that choice.

Training CV under broader grouping thresholds was diagnostic only: physicochemical
mean fold balanced accuracy was 0.730 at 60%, 0.732 at 70%, and 0.735 at 80%.
The largest broad-audit components contained 79/772, 47/772 and 44/772 sequences,
respectively. At 60%, 148/394 training identities lay in components linking to
held-out families; at 70%, 68/394 did. These are transitive component links,
not counts of direct test neighbors. No partitions or primary models were changed.

## Candidate prioritization

Overall primary-model results:

| Ranking pool | Precision@5 (95% interval) | Precision@10 | Precision@20 |
|---|---:|---:|---:|
| All 71 identities | 0.80 (0.40–1.00) | 0.80 (0.50–1.00) | 0.736 (0.520–0.900) |
| One representative per family (47) | 0.80 (0.40–1.00) | 0.80 (0.50–1.00) | 0.683 (0.450–0.875) |

The identity-pool active fraction is 34/71=47.9%. Probabilities are rounded to
12 decimals for ranking ties; boundary ties use expected hits under uniform
selection, so precision can be fractional. A family representative is the
highest-scoring identity, with label-independent lexical tie-breaking. Intervals
resample families with fixed predictions; model scores are not demonstrated
calibrated probabilities or prospective wet-lab hit rates.

## Training-only MIC cutoff sensitivity

The 32 µM reconstruction exactly replayed all 394 original training labels.
No models were fit or held-out partition outcomes used for this diagnostic.

| Threshold | Eligible original training identities | Active | Inactive | Excluded | Definite label changes from 32 |
|---|---:|---:|---:|---:|---:|
| 16 µM | 377 | 170 | 207 | 17 | 30 |
| 32 µM | 394 | 217 | 177 | 0 | 0 |
| 64 µM | 380 | 242 | 138 | 14 | 25 |

Among all 639 audit identities assigned to training families, eligible totals
were 392 at 16 µM, 394 at 32 µM and 394 at 64 µM. Equal totals at 32 and 64 do
not mean equal membership: at 64, 14 originally eligible identities were
excluded and 14 previously unready identities became eligible. At 16, 17 were
excluded and 15 became eligible. These new cases were not added to training.
This measures label/eligibility sensitivity, not model performance at new cutoffs.
The primary MIC cutoff remains 32 µM.

## Recorded run locations

- Split: `grouped-v1-dbe608485cce9a7e`.
- Comparison: `results/comparisons/feature-comparison-v1-3c0530a63f1c77d0`.
- Pretest: `results/pretest/pretest-v1-cbd14c16acf65e99`.
- Sensitivity: `results/label_sensitivity/training-label-sensitivity-v1-d536b5d4387dd9f6`.
- Final evaluation: `results/evaluation/holdout-evaluation-v1-836296044b94ac62`.

Preserve the evaluator's `metrics.csv`, `intervals.csv`, `precision_at_k.csv`,
`predictions.csv`, `evaluation_report.json` and pretest `evaluation_started.json`.
No refitting, cutoff changes or test-driven feature selection are part of this version.
