"""Saved-input report figure: quality."""
from ._common import json, pd, frame, save, aligned_annotations

def plot(source, output, resources):
    table = aligned_annotations(source)
    settings = json.loads((resources / 'config/analysis_settings.json').read_text())
    fig, axes = frame(f"Selected-entry quality — n={len(table)}; missing values remain explicit", columns=3, size=(15, 5))
    for axis, column, label, cutoff in zip(axes, ['resolution_A', 'r_free', 'mean_B_factor_A2'], ['Resolution (Å)', 'Depositor R-free', 'Mean B (Å²)'], ['qc_resolution_A', 'qc_r_free', 'qc_mean_B_A2']):
        numeric = pd.to_numeric(table[column], errors='coerce')
        axis.hist(numeric.dropna(), bins=20, color='#4477aa')
        axis.axvline(settings[cutoff], color='#cc3311', ls='--', label=f"Cutoff {settings[cutoff]:g}")
        axis.set(xlabel=label, ylabel='Selected entries', title=f"Known {int(numeric.notna().sum())}/{len(table)}")
        axis.legend(fontsize=8)
    save(fig, output, '03_quality', 'Selected one-chain-per-entry cohort; pass requires all three known metrics at or below the configured cutoffs. Fail takes precedence over missing data. Counts are metadata QC, not evidence of biological stability.')
