"""Saved-input report figure: rmsf_coverage."""
from ._common import read_table, frame, save

def plot(source, output, resources):
    profiles = read_table(source, 'ensemble_profiles.csv')
    annotations = read_table(source, 'helix_annotations.csv')
    fig, axes = frame('Experimental-ensemble CA RMSF and available-entry coverage', rows=2, size=(14, 7))
    axes[0].plot(profiles.esr1_position, profiles.ensemble_RMSF_A, color='#4477aa')
    axes[0].set(ylabel='Population ensemble RMSF (Å)')
    axes[1].plot(profiles.esr1_position, profiles.n_observed, color='#228833')
    axes[1].set(ylabel='Observed entry frames', xlabel='ESR1 canonical/native residue')
    terminal = annotations.loc[annotations.interval.eq('H12_reference')].iloc[0]
    for axis in axes:
        axis.axvspan(terminal.start, terminal.end, color='#ee6677', alpha=0.18, label='Provisional H12*')
    axes[0].legend()
    save(fig, output, '11_rmsf_coverage', 'One experimental chain per entry in the saved core-reference frame; reference first replaces its same-entry selection. Population RMSF uses per-position available observations, with missing positions retained. Experimental structural variation is not a molecular-dynamics trajectory; coverage is not an independent-replicate count.')
