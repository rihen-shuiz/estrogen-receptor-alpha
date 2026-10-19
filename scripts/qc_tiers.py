"""Tracker K05: graded QC tiers (A / B / C) per entry, instead of a hard cutoff. Labels only; nothing is removed.

Tier A limits = the old flags: analysis_settings.json qc_resolution_A, qc_r_free, qc_mean_B_A2. The script checks that
esr_analysis/resources/config/qc_tiers.json holds the same A values, so the old flags stay the reference.
Tier B limits = relaxed starting values in esr_analysis/resources/config/qc_tiers.json. Edit that file and rerun.
If the file is missing, it is created once with the starting values; an existing file is never overwritten.
Tier C = beyond the tier B limit. A missing metric is tier unknown, never a pass.
A value equal to a limit belongs to the better tier, as in the old flags (fail is strictly above the limit).
qc_label: C if any metric is C; otherwise unknown if any metric is missing; otherwise B if any metric is B; otherwise A.
Fail precedes unknown, as in the project QC rule (flag_quality.quality_status).
Outputs: results/qc_tiers.csv (one row per entry, with a plain-language warning) and qc_tiers.json (created if missing).
The old qc_* columns in results/entries.csv are read only, to check the A tier against them.
Run from the repo root: python scripts/qc_tiers.py
"""
import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
TIERS_PATH = ROOT / "esr_analysis/resources/config/qc_tiers.json"
ENTRIES = ROOT / "results/entries.csv"
OUT = ROOT / "results/qc_tiers.csv"
METRICS = {  # metric key -> (value column, old flag column, unit text, limit key in qc_tiers.json)
    "resolution": ("resolution", "qc_res", "A", "resolution_max_A"),
    "r_free": ("r_free", "qc_rfree", "", "r_free_max"),
    "mean_b": ("mean_b", "qc_b", "A^2", "mean_b_max_A2"),
}
START_VALUES = {
    "name": "qc_tiers",
    "tier_A_limits": {"resolution_max_A": SETTINGS["qc_resolution_A"], "r_free_max": SETTINGS["qc_r_free"],
                      "mean_b_max_A2": SETTINGS["qc_mean_B_A2"]},
    "tier_B_limits": {"resolution_max_A": 3.5, "r_free_max": 0.40, "mean_b_max_A2": 80.0},
    "tier_A_source": "analysis_settings.json: qc_resolution_A, qc_r_free, qc_mean_B_A2 (the old flags)",
    "tier_B_source": "starting values; the team may edit these and rerun scripts/qc_tiers.py",
    "tier_C": "beyond the tier B limit",
    "missing_metric": "unknown, never a pass",
    "qc_label_rule": "C if any metric is C; else unknown if any metric is missing; else B if any metric is B; else A",
}


def load_tiers():
    if not TIERS_PATH.exists():
        TIERS_PATH.write_text(json.dumps(START_VALUES, indent=2) + "\n")
    tiers = json.loads(TIERS_PATH.read_text())
    expected_a = {"resolution_max_A": SETTINGS["qc_resolution_A"], "r_free_max": SETTINGS["qc_r_free"],
                  "mean_b_max_A2": SETTINGS["qc_mean_B_A2"]}
    if tiers["tier_A_limits"] != expected_a:
        raise ValueError("qc_tiers.json tier A limits differ from the old flags in analysis_settings.json")
    return tiers


def tier_of(value, a_limit, b_limit):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "unknown"
    if value <= a_limit:
        return "A"
    if value <= b_limit:
        return "B"
    return "C"


def main():
    tiers = load_tiers()
    entries = pd.read_csv(ENTRIES, dtype=str, keep_default_na=False)
    rows = []
    for row in entries.itertuples(index=False):
        record = {"entry_id": row.entry_id}
        tier_names = {}
        warnings = []
        for metric, (column, flag, unit, limit_key) in METRICS.items():
            text = getattr(row, column)
            value = float(text) if text != "" else math.nan
            tier = tier_of(value, tiers["tier_A_limits"][limit_key], tiers["tier_B_limits"][limit_key])
            tier_names[metric] = tier
            record[metric] = value
            record[f"tier_{metric}"] = tier
            flag_value = getattr(row, flag)
            if (tier == "A") != (flag_value == "pass") or (tier == "unknown") != (flag_value == "unknown"):
                raise ValueError(f"{row.entry_id}: tier {tier} for {metric} disagrees with old flag {flag_value}")
            if tier == "unknown":
                warnings.append(f"{metric.replace('_', '-')} missing (unknown)")
            elif tier in ("B", "C"):
                decimals = {"resolution": 2, "r_free": 3, "mean_b": 1}[metric]
                warnings.append(f"{metric.replace('_', '-')} {value:.{decimals}f}{' ' + unit if unit else ''} (tier {tier})")
        values = list(tier_names.values())
        if "C" in values:
            label = "C"
        elif "unknown" in values:
            label = "unknown"
        elif "B" in values:
            label = "B"
        else:
            label = "A"
        record["qc_label"] = label
        record["warning"] = "; ".join(warnings) if warnings else "none"
        record["old_qc_pass"] = row.qc_pass
        rows.append(record)
    table = pd.DataFrame(rows).rename(columns={"resolution": "resolution_A", "mean_b": "mean_b_A2",
                                               "tier_r_free": "tier_rfree", "tier_mean_b": "tier_b"})
    table = table[["entry_id", "resolution_A", "r_free", "mean_b_A2", "tier_resolution", "tier_rfree", "tier_b",
                   "qc_label", "warning", "old_qc_pass"]]
    table.to_csv(OUT, index=False)
    print(f"K05: {len(table)} entries; qc_label {table.qc_label.value_counts().to_dict()}")
    print(f"K05: tier counts resolution {table.tier_resolution.value_counts().to_dict()}, "
          f"R-free {table.tier_rfree.value_counts().to_dict()}, mean B {table.tier_b.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
