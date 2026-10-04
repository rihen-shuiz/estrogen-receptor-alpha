"""Saved-input report figure: identity_action."""
from ._common import pd, read_table, frame, save

def plot(source, output, resources):
    chemistry = read_table(source, 'ligand_identities.csv')
    units = []
    for identity, group in chemistry.groupby('chemical_identity_key'):
        labels = sorted(set(group.loc[group.evidence_supported, 'supported_action_type'].dropna()))
        units.append(labels[0] if len(labels) == 1 else 'unknown' if not labels else 'conflicting')
    counts = pd.Series(units).value_counts().reindex(['agonist', 'degrader', 'unknown', 'conflicting'], fill_value=0)
    fig, axes = frame(f"ESR1 exact chemical identities and supported action evidence — n={len(units)}", columns=2, size=(13, 5))
    axes[0].bar(['CCD/context rows', 'Exact identities'], [len(chemistry), len(units)], color='#4477aa')
    axes[0].set(ylabel='ESR1 identity population')
    axes[1].barh(counts.index, counts.values, color=['#228833', '#aa3377', '#bbbbbb', '#cc3311'])
    for n, value in enumerate(counts.values):
        axes[1].text(value, n, f" {value}", va='center')
    axes[1].set(xlabel='Exact identities; current evidence-support policy', xlim=(0, len(units) * 1.12))
    save(fig, output, '04_identity_action', 'Full canonical-isomeric-SMILES identity, without salt/charge/tautomer/stereochemical merging. Only source-supported exact-compound receptor-context actions count; unknowns are not inferred from generic activity or scaffold similarity. Later expanded curation is excluded.')
