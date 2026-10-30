"""Predict with the exported frozen two-feature model, using only Python stdlib.

No pickle, sklearn, network, fitting, data download or held-out evaluation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

STANDARD = set('ACDEFGHIKLMNPQRSTVWY')
EXPECTED_SOURCE = '5573f24eb93693b207b8e356c7e842433c1f9f26b08fd1fea6039dc85e31680e'
DEFAULT_MODEL = Path(__file__).resolve().parents[1] / 'artifacts/model-v1/primary_model.json'


def validate_model(model):
    if (model.get('schema') != 'amp-primary-linear-v1'
            or model.get('feature_names') != ['charge_proxy', 'c_amidated']
            or model.get('source_model_sha256') != EXPECTED_SOURCE
            or model.get('classes') != [0, 1]
            or model.get('prediction_threshold') != 0.5):
        raise ValueError('Unsupported model schema or source fingerprint')
    for key in ('mean', 'scale', 'coef'):
        values = model.get(key)
        if not isinstance(values, list) or len(values) != 2:
            raise ValueError(f'Invalid {key}')
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
            raise ValueError(f'Non-finite or nonnumeric {key}')
    if any(v <= 0 for v in model['scale']):
        raise ValueError('Scales must be positive')
    if type(model.get('intercept')) not in (int, float) or not math.isfinite(model['intercept']):
        raise ValueError('Invalid intercept')
    return model


def load_model(path):
    path = Path(path)
    manifest = json.loads((path.parent / 'model_manifest.json').read_text(encoding='utf-8'))
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != manifest['files'][path.name]:
        raise ValueError('Portable model fingerprint mismatch')
    return validate_model(json.loads(content))


def predict(sequence, amidated, model):
    if type(amidated) is not bool:
        raise ValueError('Amidation must be specified as a boolean')
    if not isinstance(sequence, str) or not 10 <= len(sequence) <= 50 or not set(sequence) <= STANDARD:
        raise ValueError('Use canonical uppercase sequences of length 10-50')
    features = [sequence.count('K') + sequence.count('R') - sequence.count('D') - sequence.count('E') + int(amidated), int(amidated)]
    decision = model['intercept'] + sum((v-m)/s*c for v, m, s, c in zip(features, model['mean'], model['scale'], model['coef']))
    if decision >= 0:
        probability = 1 / (1 + math.exp(-decision))
    else:
        e = math.exp(decision)
        probability = e / (1 + e)
    return {'sequence': sequence, 'c_amidated': amidated, 'charge_proxy': features[0],
            'model_probability_active': probability, 'predicted_label': 'active' if decision > 0 else 'inactive',
            'activity_definition': 'recorded E. coli MIC <= 32 uM',
            'model': 'charge_corrected_amidation', 'source_model_sha256': model['source_model_sha256']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequence', required=True)
    chemistry = parser.add_mutually_exclusive_group(required=True)
    chemistry.add_argument('--amidated', action='store_true')
    chemistry.add_argument('--unmodified', action='store_true')
    parser.add_argument('--model', type=Path, default=DEFAULT_MODEL)
    args = parser.parse_args()
    try:
        print(json.dumps(predict(args.sequence, args.amidated, load_model(args.model)), indent=2))
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f'Prediction failed: {error}\n')


if __name__ == '__main__':
    main()
