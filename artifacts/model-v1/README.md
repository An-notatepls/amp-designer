# Frozen primary model (portable JSON)

This is the evaluated charge-plus-amidation logistic model, exported from the
original trusted joblib bundle without refitting. Training used 394 identities,
standardized features, C=0.1, and a probability threshold of 0.5 (an exact tie
is inactive). The target is recorded E. coli MIC <=32 µM.

Source joblib SHA-256:
`5573f24eb93693b207b8e356c7e842433c1f9f26b08fd1fea6039dc85e31680e`.
`model_manifest.json` records the portable file checksum. No pickle is needed
for inference. From the project root:

```bash
python3 modeling/predict_portable.py --sequence AAAAAAAAAA --unmodified
```

Model outputs are uncalibrated scores for a narrowly filtered database task,
not experimental confirmation, toxicity predictions or clinical efficacy.
The synthetic input above demonstrates syntax only. DBAASP-derived training
data retain source attribution; see ../../data/README.md and the main README.
