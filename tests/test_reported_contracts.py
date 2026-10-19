"""Focused first-month scientific contracts, using only installed resources."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from esr_analysis._task_io import RESOURCE_ROOT, coordinate_mapping
from esr_analysis.align_structures import apply_transform, fit_to_reference
from esr_analysis.calculate_rmsd import compute_distances
from esr_analysis.cluster_structures import structural_groups
from esr_analysis.flag_quality import quality_status
from esr_analysis.run_io import read_table, positions

SOURCE = RESOURCE_ROOT / 'frozen'
SETTINGS = json.loads((RESOURCE_ROOT / 'config/analysis_settings.json').read_text())


def test_master_subsets_contexts_and_QC_denominators():
    ids = (RESOURCE_ROOT / 'data/pdb_ids_esr1.txt').read_text().splitlines()
    entries = read_table(SOURCE, 'entries.csv')
    chains = read_table(SOURCE, 'receptor_chains.csv')
    target = chains.loc[chains.identity.eq('esr1')]
    selected = read_table(SOURCE, 'entry_annotations.csv')
    assert len(ids) == len(set(ids)) == len(entries) == 487
    assert set(ids) == set(entries.pdb_id)
    assert len(target) == 901 and target.pdb_id.nunique() == 444
    assert target.eligible_historical.sum() == 788
    assert len(selected) == 351
    assert selected.qc_status.value_counts().to_dict() == {'pass': 247, 'fail': 104}
    assert entries.qc_status.value_counts().to_dict() == {'pass': 369, 'fail': 117, 'unknown': 1}
    assert selected.declared_receptor_chains.gt(1).sum() == 319
    assert selected.partner_chains.gt(0).sum() == 211
    assert selected.dna_chains.gt(0).sum() == 0


def test_identity_and_supported_action_units():
    table = read_table(SOURCE, 'ligand_identities.csv')
    assert table.receptor_context.eq('ESR1').all()
    assert table.chemical_identity_key.nunique() == 400
    units = table.drop_duplicates('chemical_identity_key')
    assert units.supported_action_type.fillna('unknown').value_counts().to_dict() == {'unknown': 393, 'degrader': 4, 'agonist': 3}
    known = units.loc[units.evidence_supported]
    assert known.supported_action_evidence.notna().all()
    assert known.supported_action_receptor_context.str.contains('ESR1').all()


def test_reference_and_fitted_coordinate_correspondence():
    reference = json.loads((SOURCE / 'reference.json').read_text())
    assert reference['chain_key'] == '7NFB|B|1'
    chains = read_table(SOURCE, 'receptor_chains.csv').set_index('chain_key')
    assert chains.loc['7NFB|B|1', ['mapped_lbd_ca', 'mapped_core_ca', 'mapped_h12_ca']].tolist() == [244, 84, 10]
    mapping = read_table(SOURCE, 'residue_mapping.csv')
    mobile = coordinate_mapping(mapping, '3Q95|B|1')
    target = coordinate_mapping(mapping, '7NFB|B|1')
    aligned, fitted = fit_to_reference(mobile, target, positions(SETTINGS['esr1_core_ranges']), 60)
    assert fitted['n_common_core_ca'] == 84
    saved = json.loads((SOURCE / 'alignment_transforms.json').read_text())['3Q95|B|1']
    keys = sorted(mobile)
    expected = apply_transform(np.array([mobile[key] for key in keys]), np.asarray(saved['rotation']), np.asarray(saved['mobile_centroid']), np.asarray(saved['reference_centroid']))
    assert np.allclose(np.array([aligned[key] for key in keys]), expected, atol=1e-10, rtol=0)


def test_indexed_matrix_partition_and_pair_medians():
    ids = (SOURCE / 'rmsd_entry_index.txt').read_text().splitlines()
    matrix = np.load(SOURCE / 'rmsd_entry.npy')
    assert matrix.shape == (351, 351) and len(set(ids)) == 351
    assert np.isfinite(matrix).all() and np.allclose(matrix, matrix.T, atol=1e-12)
    assert np.array_equal(np.diag(matrix), np.zeros(351))
    assignments, _ = structural_groups(ids, matrix, [1.5], 1.5)
    saved = read_table(SOURCE, 'entry_clusters.csv')
    assert assignments.chain_key.tolist() == saved.chain_key.tolist()
    assert np.array_equal(assignments['cluster_1.5A'], saved['cluster_1.5A'])
    counts = assignments['cluster_1.5A'].value_counts().sort_index().tolist()
    assert counts == [2, 2, 2, 124, 1, 6, 212, 1, 1]
    labels = assignments['cluster_1.5A'].to_numpy()
    left, right = np.triu_indices(len(ids), 1)
    assert len(left) == 61425
    same = labels[left] == labels[right]
    assert round(float(np.median(matrix[left[same], right[same]])), 3) == .800
    assert round(float(np.median(matrix[left[~same], right[~same]])), 3) == 3.820


def test_fail_precedes_unknown_and_thresholds_are_inclusive():
    assert quality_status(np.nan, .4, 10, SETTINGS) == 'fail'
    assert quality_status(np.nan, .2, 10, SETTINGS) == 'unknown'
    assert quality_status(2.5, .3, 50, SETTINGS) == 'pass'


def test_no_pair_refit_and_missing_overlap_is_not_zero_distance():
    a = {p: np.array([p, 0., 0.]) for p in range(3)}
    b = {p: value + np.array([3., 0., 0.]) for p, value in a.items()}
    ids, matrix, overlap, excluded = compute_distances(['a', 'b'], {'a': a, 'b': b}, min_common=3)
    assert matrix[0, 1] == 3 and overlap[0, 1] == 3
    _, missing, _, _ = compute_distances(['a', 'b'], {'a': a, 'b': b}, min_common=4)
    assert np.isnan(missing[0, 1])


@pytest.mark.parametrize('interval,count,full,qc,variation', [
    ('H6', 328, 1.144, 1.086, .703),
    ('H12_reference', 344, .339, .325, 6.500),
    ('H5', 351, -.846, -.830, .423),
    ('H8', 351, -.860, -.857, .334),
])
def test_reported_regional_values(interval, count, full, qc, variation):
    table = read_table(SOURCE, 'helix_distribution_and_coverage.csv')
    table = table.loc[table.interval.eq(interval) & table.coverage_threshold.eq(.8)].set_index('cohort')
    assert table.loc['full'].retained_chains == count
    assert np.isclose(table.loc['full'].median_chain_median_z, full, atol=.0005, rtol=0)
    assert np.isclose(table.loc['combined_QC'].median_chain_median_z, qc, atol=.0005, rtol=0)
    ensemble = read_table(SOURCE, 'helix_ensemble_summary.csv').set_index('interval')
    assert np.isclose(ensemble.loc[interval].RMS_position_ensemble_RMSF_A, variation, atol=.0005, rtol=0)


def test_provisional_H12_and_population_SD_policy():
    helices = read_table(SOURCE, 'helix_annotations.csv').set_index('interval')
    assert helices.loc['H12_reference', ['start', 'end']].tolist() == [536, 545]
    assert 'provisional' in str(helices.loc['H5'].status)
    norms = read_table(SOURCE, 'chain_B_normalisation.csv').set_index('chain_key')
    mapping = read_table(SOURCE, 'residue_mapping.csv')
    rows = mapping.loc[mapping.chain_key.eq('7NFB|B|1') & mapping.has_ca & mapping.canonical_position.between(305, 550)]
    values = rows.loc[np.isfinite(rows.B_factor_A2) & rows.B_factor_A2.ge(0), 'B_factor_A2']
    assert np.isclose(values.std(ddof=0), norms.loc['7NFB|B|1'].SD_B_A2_ddof0, atol=1e-12)
