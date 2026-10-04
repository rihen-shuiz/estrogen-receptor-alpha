"""Saved-input report figure: group_counts."""
from ._common import np, frame, save, aligned_annotations

def plot(source, output, resources):
    table = aligned_annotations(source)
    total = table.groupby('struct_cluster').size().sort_index()
    passing = table.loc[table.qc_status.eq('pass')].groupby('struct_cluster').size().reindex(total.index, fill_value=0)
    fig, axes = frame('Original baseline groups and their QC-pass members', size=(12, 6))
    axis = axes[0]
    x = np.arange(len(total))
    axis.bar(x, total.values, color='#bbbbbb', label=f"All selected entries (n={len(table)})")
    axis.bar(x, passing.values, color='#228833', label=f"Known QC pass (n={int(passing.sum())})")
    for i, n in enumerate(total):
        axis.text(i, n, f"{n}", ha='center', va='bottom')
    axis.set_xticks(x, total.index)
    axis.set(xlabel='Saved baseline group ID', ylabel='Selected entries')
    axis.legend()
    save(fig, output, '08_group_counts', 'QC-pass overlay retains membership in existing baseline groups. It is not independent QC reclustering or a new partition. Original labels and all nine groups remain unchanged.')
