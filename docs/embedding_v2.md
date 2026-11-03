# Does a protein embedding improve the predictor?

We tested the next idea after reviewing the training notes: replace the simple
charge features with a numerical description from a pretrained protein model.
The first result is a small, uncertain difference, not a demonstrated gain in
final predictive accuracy. The shipped v1 predictor remains unchanged.

## What is an embedding?

An embedding is a list of numbers that summarizes a sequence using patterns
learned from other protein sequences. Here we used the small, frozen
[ESM-2 8M model](https://huggingface.co/facebook/esm2_t6_8M_UR50D), following
the [original ESM project](https://github.com/facebookresearch/esm). Its weights
were downloaded at revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`.
Weights are cached locally, not included in this repository.

We averaged the final-layer representation over actual amino-acid tokens,
excluding start, end and padding tokens. That gives 320 numbers per sequence.
We then added the amidation flag, because changing this terminal modification
does not change the amino-acid letters. The encoder was not fine-tuned; only a
standardized logistic regression classifier was fitted inside each fold.

## A comparison that separates two questions

The [training review](training_v2.md) recovered 48 eligible identities, increasing
the training candidate pool from 394 to 442. To distinguish a dataset-size
effect from a representation effect, we compared three arms on identical
validation rows and identical five-fold sequence-family assignments:

| Arm | Rows available for fitting within each training fold | Pooled validation accuracy | Balanced accuracy | AUROC |
|---|---|---:|---:|---:|
| Charge + amidation, original-data arm | Original 394 identities only | 73.8% | 0.732 | 0.752 |
| Charge + amidation, expanded-data arm | All 442 candidates | 73.8% | 0.732 | 0.753 |
| ESM-2 + amidation, expanded-data arm | All 442 candidates | 74.2% | 0.738 | 0.768 |

Each row is predicted once, by a model that did not fit that identity or its
sequence family. Even the original-data arm predicts the same 442 validation
identities as the expanded arms. Within a training fold, however, it only fits
the rows that belonged to the original ready dataset.

The five-fold grouping and seed are fixed for this experiment. They differ
from the original v1 cross-validation setup because the eligible pool now has
442 identities in 267 families. The table compares arms within this experiment;
it should not be compared directly with the v1 held-out accuracy of 71.8%.

The two charge models have the same thresholded predictions in this run,
although their probability scores differ slightly. The embedding model has a
net two additional correct predictions out of 442. That is an accuracy
difference of about **0.45 percentage points**.

## What the uncertainty tells us

We resampled whole sequence families 2,000 times, using the same draws for
every arm. The paired embedding-minus-expanded-charge differences were:

| Metric | Estimated difference | 95% family-bootstrap interval |
|---|---:|---:|
| Accuracy | +0.45 percentage points | −6.01 to +6.68 percentage points |
| Balanced accuracy | +0.005 | −0.054 to +0.063 |
| AUROC | +0.015 | −0.058 to +0.089 |

Each interval includes zero, so the current evidence does not establish a
reliable improvement. These intervals describe fixed, selected development
predictions. They do not account for every source of uncertainty, such as
hyperparameter selection, refitting, repeated development or uncertain labels.
The encoder's pretraining overlap with published peptides was not audited.
The frozen 80% global-edit sequence families are also a heuristic: truncations,
short shared motifs and shared experimental sources can still cross CV folds.

The charge arms use fixed C=0.1. For the embedding classifier, we searched the
declared C values 0.0001, 0.001, 0.01, 0.1, 1 and 10, selecting by mean-fold
balanced accuracy. C=10 was selected, at the upper edge of that grid. This
does not identify an interior optimum. Model choice and settings were recorded
in [the protocol](embedding_v2_protocol.json), committed in
`029d63078dbaab23373c85308653800cf6196492` before the comparison ran.

## Rerun the experiment

Keep the embedding environment separate from the frozen v1 environment:

```bash
python3 -m venv .venv-embeddings
./.venv-embeddings/bin/python -m pip install -r modeling/requirements-embeddings.txt
./.venv-embeddings/bin/python modeling/compare_embeddings.py
```

On Linux, a CPU-only PyTorch installation can be installed first using the
instructions on [PyTorch's installation page](https://pytorch.org/get-started/locally/).
This experiment used Python 3.12.14 and torch 2.9.1+cpu, with transformers
4.57.1 and the recorded baseline-library versions. The direct dependencies are
pinned; transitive dependencies are not fully locked. The full embedding run
has not been verified on macOS or Python 3.13.4. The numerical baseline remains
in its own original environment.

The first run downloads the specified model files and caches one embedding per
unique training sequence under `data/processed/embedding_cache/`. The script
checks token counts, source fingerprints, cached file hashes and finite values.
An unchanged rerun verifies and reuses the existing result directory. Different
platforms or library versions can produce different floating-point results.

The output folder under `results/development/` contains:

- `comparison_report.json`: selected C values, metrics, environment and hashes.
- `cv_scores.csv`: every arm, C and fold score.
- `cv_folds.csv`: shared validation-fold membership for all 442 identities.
- `selected_oof_predictions.csv`: each arm's selected validation predictions.
- `intervals.csv` and `paired_differences.csv`: uncertainty and paired changes.
- `embedding_manifest.json`: fixed model revision, pooling method, model-file
  hashes, cache hash and library versions.

The downloaded weights and embedding arrays are generated locally. CI runs the
offline synthetic tests without downloading a model. The actual run above
checks that all 439 unique sequences produce 320-dimensional representations.

## What we keep and what comes next

The original simple predictor remains the released model. This experiment adds
evidence about a plausible alternative without establishing better final
accuracy. A new independent evaluation would be needed for that claim because
the v1 test and clinical sets have already been evaluated.

Further work should prioritize source-level review of the remaining assay notes
and an independent collection suitable for evaluation. Trying additional
representations is reasonable, but a more complex model is not by itself
evidence of an improvement.
