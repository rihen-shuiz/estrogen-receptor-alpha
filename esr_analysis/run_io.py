"""Keep inputs and saved results compatible. Hash local inputs and settings;
never reuse saved geometry after its cohort, mappings or fitting definition changes.
Outputs are plain CSV/JSON/NPZ files; no scientific calculations happen here.
"""
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def save_table(out, name, rows, columns=None):
    table = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows, columns=columns)
    table.to_csv(Path(out) / name, index=False)
    return table


def read_table(out, name):
    path = Path(out) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing prerequisite {path}. Run the producing stage first.")
    # "NA" is a real sodium CCD code (and can be a label chain ID).
    # Our CSV writer represents missing values with empty fields only.
    return pd.read_csv(path, keep_default_na=False, na_values=[""], low_memory=False)


def positions(ranges):
    result = []
    for start, end in ranges:
        result.extend(range(start, end + 1))
    return result


def geometry_inputs(out, settings):
    """Fingerprint atom correspondence, reference and fitting definition."""
    return {"residue_mapping_sha256":sha256(Path(out)/"residue_mapping.csv"),
            "reference_sha256":sha256(Path(out)/"reference.json"),
            "chain_inventory_sha256":sha256(Path(out)/"receptor_chains.csv"),
            "core_ranges":settings["esr1_core_ranges"],
            "lbd_range":[settings["esr1_lbd_start"],settings["esr1_lbd_end"]],
            "min_core_ca":settings["min_core_ca_for_alignment"],
            "min_common_ca":settings["min_pairwise_common_ca"]}


def validate_saved_geometry(out, settings):
    """Reject a matrix/ensemble after its mapping, reference or fit changed."""
    signature = Path(out)/"geometry_inputs.json"
    if not signature.exists() or json.loads(signature.read_text()) != geometry_inputs(out,settings):
        raise ValueError("Saved geometry no longer matches mapping/reference/core selections. Use a new run and recompute global-rmsd.")
