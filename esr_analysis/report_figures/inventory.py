"""Saved-input report figure: inventory."""
from ._common import pd, read_table, frame, save, aligned_annotations

def plot(source, output, resources):
    entries = read_table(source, 'entries.csv')
    chains = read_table(source, 'receptor_chains.csv')
    target = chains.loc[chains.identity.eq('esr1')]
    selected = aligned_annotations(source)
    values = [len(entries), target.pdb_id.nunique(), target.loc[target.eligible_historical].pdb_id.nunique()]
    fig, axes = frame('ESR1 entry inventory — full master and explicit selected subset', columns=2, size=(14, 6))
    axes[0].barh(['Master entries', 'Verified ESR1 chain present', 'Applicable LBD entry'], values, color='#4477aa')
    for n, value in enumerate(values):
        axes[0].text(value, n, f" {value}", va='center')
    axes[0].set(xlabel='Distinct deposited entries', xlim=(0, max(values) * 1.16))
    contexts = pd.Series({'Multiple receptor copies': selected.declared_receptor_chains.gt(1).sum(), 'Partner-containing': selected.partner_chains.gt(0).sum(), 'DNA-containing': selected.dna_chains.gt(0).sum()})
    axes[1].barh(contexts.index, contexts.values, color='#66ccee')
    for n, value in enumerate(contexts.values):
        axes[1].text(value, n, f" {value}", va='center')
    axes[1].set(xlabel=f"Selected entries (n={len(selected)}); categories overlap", xlim=(0, len(selected) * 1.1))
    save(fig, output, '02_inventory', f"{len(entries)} full master entries; {len(target)} verified ESR1 chain/models and {int(target.eligible_historical.sum())} eligible copies. Fragment entries remain in the master. Context categories overlap and are not summed as attrition.")
