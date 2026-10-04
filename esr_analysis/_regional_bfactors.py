"""Reported ESR1 portion of the original regional B producer.

Population SD, valid-B observations, matching positions, coverage and region
aggregation are unchanged. ESR2, matched-pair and sensitivity extensions are
excluded. Annotation inputs retain their source and provisional status.
"""
from pathlib import Path
import numpy as np
import pandas as pd

def run(*, source_folder, output_folder, annotation_file):
    source, out = (Path(source_folder), Path(output_folder))
    out.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(source / 'residue_mapping.csv', low_memory=False)
    entries = pd.read_csv(source / 'entries.csv').set_index('pdb_id')
    primary = (source / 'rmsd_entry_index.txt').read_text().splitlines()
    annotations = pd.read_csv(annotation_file)
    if not annotations.receptor.eq('esr1').all():
        raise ValueError('This reported producer requires explicit ESR1 annotations')
    annotations.to_csv(out / 'helix_annotations.csv', index=False)
    selected = mapping.loc[mapping.chain_key.isin(primary) & mapping.canonical_position.between(305, 550) & mapping.has_ca].copy()
    if selected.duplicated(['chain_key', 'canonical_position']).any():
        raise ValueError('Duplicate observed chain/native positions')
    if not (selected.occupancy.isna() | selected.occupancy.gt(0)).all():
        raise ValueError('Unusable CA occupancy in selected observations')
    selected['valid_B'] = np.isfinite(selected.B_factor_A2) & selected.B_factor_A2.ge(0)
    selected['B_z'] = np.nan
    norms = []
    for key, group in selected.groupby('chain_key', sort=False):
        valid = group.loc[group.valid_B, 'B_factor_A2']
        mean, sd = (valid.mean(), valid.std(ddof=0))
        if len(valid) >= 2 and sd > 0:
            selected.loc[valid.index, 'B_z'] = (valid - mean) / sd
        norms.append({'chain_key': key, 'valid_B_CA': len(valid), 'invalid_B_CA': int((~group.valid_B).sum()), 'zero_B_CA': int(group.B_factor_A2.eq(0).sum()), 'mean_B_A2': mean, 'SD_B_A2_ddof0': sd, 'normalisation_available': len(valid) >= 2 and sd > 0})
    pd.DataFrame(norms).to_csv(out / 'chain_B_normalisation.csv', index=False)
    selected.loc[~selected.valid_B].to_csv(out / 'invalid_CA_B.csv', index=False)
    local = pd.read_csv(source / 'local_displacements.csv')
    selected = selected.merge(local[['chain_key', 'esr1_position', 'displacement_A']], on=['chain_key', 'esr1_position'], how='left', validate='one_to_one')
    selected[['chain_key', 'canonical_position', 'esr1_position', 'B_factor_A2', 'valid_B', 'B_z', 'displacement_A']].to_csv(out / 'matched_B_and_displacement.csv', index=False)
    rows = []
    for key in primary:
        group = selected.loc[selected.chain_key.eq(key)]
        for region in annotations.itertuples():
            observed = group.loc[group.canonical_position.between(region.start, region.end)]
            valid = observed.loc[observed.valid_B]
            matched = valid.loc[np.isfinite(valid.displacement_A)]
            expected = region.end - region.start + 1
            rows.append({'chain_key': key, 'pdb_id': key.split('|')[0], 'receptor': 'esr1', 'QC': entries.loc[key.split('|')[0]].qc_status, 'interval': region.interval, 'start': region.start, 'end': region.end, 'expected_CA': expected, 'observed_CA': len(observed), 'valid_B_CA': len(valid), 'coverage': len(observed) / expected, 'B_coverage': len(valid) / expected, 'median_B_A2': valid.B_factor_A2.median(), 'mean_B_A2': valid.B_factor_A2.mean(), 'median_B_z': valid.B_z.median(), 'matched_B_displacement_CA': len(matched), 'matched_median_B_z': matched.B_z.median(), 'reference_RMSD_A': np.sqrt(np.mean(matched.displacement_A ** 2)) if len(matched) else np.nan, 'shape_coverage': len(matched) / expected})
    summary = pd.DataFrame(rows)
    summary.to_csv(out / 'helix_chain_summaries.csv', index=False)
    distributions = []
    for interval, group in summary.groupby('interval', sort=False):
        for cohort, subset in [('full', group), ('combined_QC', group.loc[group.QC.eq('pass')])]:
            retained = subset.loc[subset.B_coverage.ge(0.8) & subset.valid_B_CA.gt(0)]
            distributions.append({'receptor': 'esr1', 'interval': interval, 'cohort': cohort, 'coverage_threshold': 0.8, 'cohort_chains': len(subset), 'retained_chains': len(retained), 'median_chain_median_B_A2': retained.median_B_A2.median(), 'median_chain_median_z': retained.median_B_z.median(), 'median_coverage': retained.coverage.median()})
    pd.DataFrame(distributions).to_csv(out / 'helix_distribution_and_coverage.csv', index=False)
    profiles = pd.read_csv(source / 'ensemble_profiles.csv')
    ensemble = []
    for region in annotations.itertuples():
        profile = profiles.loc[profiles.esr1_position.between(region.start, region.end)]
        ensemble.append({'interval': region.interval, 'available_positions': int(profile.ensemble_RMSF_A.notna().sum()), 'min_observed_entries': int(profile.n_observed.min()), 'max_observed_entries': int(profile.n_observed.max()), 'RMS_position_ensemble_RMSF_A': float(np.sqrt(np.mean(profile.ensemble_RMSF_A ** 2)))})
    pd.DataFrame(ensemble).to_csv(out / 'helix_ensemble_summary.csv', index=False)
    return summary
