# Training-data review for the next experiment

This is the first development step after the completed v1 predictor. The aim
is to recover useful training examples by interpreting unresolved assay notes.
This step does not measure an improvement in prediction accuracy.

## What changed

The original sequence-family assignments include identities whose notes were
still unresolved. We selected only families assigned to training, replayed the
original labeling rules, and verified that their original labels matched.
We then reviewed the notes using a worksheet with no separate MIC-value or
outcome-label columns. Original notes remain verbatim, including any values
mentioned in the note itself. Decisions concern assay meaning, not whether a
row supports an active label. The cutoff remains 32 µM.

| Training-data quantity | Before this review | After this review |
|---|---:|---:|
| Eligible identities | 394 | 442 |
| Unique sequences | 392 | 439 |
| Sequence families | 242 | 267 |
| Active identities | 217 | 250 |
| Inactive identities | 177 | 192 |
| Identities with status `needs_note_review` | 193 | 145 |
| Pending eligible measurement notes | 249 | 180 |

The 48 newly eligible identities comprise 33 active and 15 inactive examples.
All original 394 ready training labels are retained. The training-family audit
pool still contains 41 mixed identities and 11 ambiguous identities; this
review did not resolve those cases by voting or dropping opposing results.

## What the review establishes

The new decision file is `data/curation/training_note_reviews_v2.csv`. It has
249 measurement-keyed rows, including 69 explicit inclusion decisions covering
27 distinct note texts. Each decision is matched to the exact DBAASP note,
sequence, amidation state and frozen training family. The file records that
this is an AI-assisted annotation review; original publications were not
checked. Decisions can be inspected directly rather than inferred from a
keyword rule in the training script.

An ordinary single-isolate description, such as `Clinical isolate GD004`, is
acceptable under the existing broad task. It identifies the tested organism
without stating an added treatment or a summary over several isolates. It
does not establish reference-strain origin or verify unrecorded conditions.
Ordinary species, pathotype or resistance descriptions can also be included
under that same policy. Clinical assays are already part of the v1 broad task;
this review does not introduce a reference-only policy.

Unclear cases stay pending. Examples include:

- `8 MDR Clinical isolates`: the meaning of the combined measurement needs
  the original source.
- `also with 12.5% serum`: the note does not tell us which condition produced
  this MIC.
- `not sufficiently soluble to test at 256 micro g/mL`: the measurement may
  describe a solubility limit rather than a definite activity result.
- `Ciprofloxacin resistant strain and clinical isolate 20411`: the wording
  may describe one strain or two distinct strains.
- `Methicillin resistant wild type`: the unusual organism/resistance
  description needs source verification.

Blank decisions remain pending, with a reason explaining what evidence is
missing. These rows are not silently treated as inactive. The remaining 180
notes include notes on mixed or ambiguous identities as well as identities
whose primary status is `needs_note_review`; those counts have different units.

## Reproduce this training-only build

Use the existing modeling environment, with dependencies from `requirements.txt`.
From the repository root:

```bash
./.venv/bin/python modeling/test_review_training_notes.py -v
./.venv/bin/python modeling/review_training_notes.py
```

The second command prints a new folder under `results/development/`. It writes:

- `training_candidates.csv`: 442 eligible identities, each with a training
  family ID and an indicator of whether it belonged to the original 394.
- `training_audit.csv`: all 639 identities assigned to training families,
  including unresolved and mixed identities.
- `initial_review_queue.csv`: the starting worksheet with no separate MIC-value
  or outcome-label columns; original notes are retained verbatim.
- `remaining_review_queue.csv`: the 180 notes that still need review.
- `review_report.json`: counts, transitions, source/output hashes and scope.

The folder name depends on source and review fingerprints. An unchanged rerun
verifies existing outputs and reuses them; changed output bytes are rejected.
The committed report and tables capture this review pass. For additional
decisions, use a new versioned review file and pass it with `--reviews`.

Keep `training_candidates.csv` separate from the frozen v1 `train.csv`. It is
not a replacement split manifest and should not be passed to the old training
CLI as though it were the original dataset. A v2 modeling runner must validate
this report and construct group-aware development folds explicitly.

## How we will judge the next model

A first charge-and-amidation versus embedding comparison has now been run
using shared training-family validation folds. See [the embedding experiment](embedding_v2.md)
for its results and separate reproduction environment. Amidation remains a
separate feature because amino-acid embeddings do not represent that modification.
The comparison evaluates every arm on the same rows, separating the additional
training-data effect from the representation effect.

The embedding model, pooling method, regularization search and fold allocation
were recorded before viewing validation results. Scaling was fitted inside each
fold. The report includes accuracy, balanced accuracy, AUROC and uncertainty,
and distinguishes selected development scores from final evaluation. More
examples or a larger model do not guarantee better predictions.

The v1 test and clinical sets have already been evaluated. They are historical
results, not fresh blind tests for v2. A convincing final claim of improvement
requires new independent evaluation, with related sequences kept together.
None of this review fits a model, opens held-out classifier CSVs, changes v1
artifacts or estimates a new accuracy.
