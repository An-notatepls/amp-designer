# Antimicrobial Peptide Activity Prediction and Design

## Objective
Build a reproducible machine learning pipeline using DBAASP
data to predict antibacterial activity of short peptides.

After evaluating the predictor, extend the project to generate
and rank unvalidated peptide candidates.

## Initial Scope
- Study one bacterial species, selected after a data audit.
- Preserve strain and experimental metadata.
- Initially consider linear peptides with standard amino acids.
- Use a provisional peptide length range of 10–50 amino acids.
- Select MIC regression or activity classification based on
  the availability and quality of experimental measurements.
- Compare simple baseline models before adding deep learning.

## Success Criteria
- Document data retrieval, cleaning, and exclusions.
- Evaluate on held-out data while keeping related sequences together.
- Compare models against a simple baseline.
- Provide reproducible training and prediction instructions.
- Build a sequence-scoring demo.

## Limitations
Predictions and generated candidates are computational estimates.
They do not demonstrate bacterial killing, clinical effectiveness,
or safety. Experimental validation is outside this project's scope.

## Status
Environment setup complete. Data audit pending.