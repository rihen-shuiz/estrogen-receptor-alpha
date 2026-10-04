"""Small explicit-path CLI/file helpers, shared by the thirteen task adapters."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .run_io import read_table

RESOURCE_ROOT = Path(__file__).resolve().parent / 'resources'

def parser(description):
    result = argparse.ArgumentParser(description=description)
    result.add_argument('--input-root', type=Path, default=RESOURCE_ROOT)
    result.add_argument('--source-folder', type=Path)
    result.add_argument('--settings', type=Path)
    result.add_argument('--output', type=Path, required=True)
    return result

def prepare(args):
    root = args.input_root.resolve()
    source = (args.source_folder or root / 'frozen').resolve()
    output = args.output.resolve()
    if output.is_relative_to(root) or output == source:
        raise ValueError('Outputs must be outside the bundled inputs and separate from the source folder')
    output.mkdir(parents=True, exist_ok=True)
    settings_path = args.settings or root / 'config/analysis_settings.json'
    if not settings_path.is_file():
        raise FileNotFoundError(f'Missing settings: {settings_path}')
    return root, source, output, json.loads(settings_path.read_text())

def records(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))

def json_write(path, value):
    def scalar(item):
        return item.tolist() if hasattr(item, 'tolist') else str(item)
    Path(path).write_text(json.dumps(value, indent=2, default=scalar) + '\n')

def sequences(root, settings):
    result = {}
    for receptor in ['esr1', 'esr2']:
        path = root / 'data' / f"{receptor}_uniprot_{settings[receptor + '_uniprot']}.fasta"
        result[receptor] = ''.join(line.strip() for line in path.read_text().splitlines() if not line.startswith('>'))
    return result

def coordinate_mapping(mapping, chain_key):
    rows = mapping.loc[mapping.chain_key.eq(chain_key) & mapping.has_ca].dropna(subset=['esr1_position'])
    if rows.esr1_position.duplicated().any():
        raise ValueError('Observed chain positions are not unique')
    return {int(row.esr1_position): np.array([row.x_A, row.y_A, row.z_A]) for row in rows.itertuples()}

def aligned_caches(path, ids):
    with np.load(path, allow_pickle=False) as data:
        keys = data['chain_keys'].astype(str).tolist()
        positions, coords = data['positions'], data['coords']
        if coords.shape != (len(keys), len(positions), 3):
            raise ValueError('Aligned coordinate shape disagrees with chain/position index')
        lookup = dict(zip(keys, coords))
        return {key: {int(position): xyz for position, xyz in zip(positions, lookup[key]) if np.isfinite(xyz).all()} for key in ids}
