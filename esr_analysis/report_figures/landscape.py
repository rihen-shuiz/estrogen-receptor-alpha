"""Saved-input report figure: landscape."""
from ._common import np, GROUP_COLOURS, QC_COLOURS, frame, save, aligned_annotations

def plot(source, output, resources):
    table = aligned_annotations(source)
    matrix = np.load(source / 'rmsd_entry.npy')
    variance = table[['PC1', 'PC2']].var(ddof=1).to_numpy() / np.var(matrix, axis=0, ddof=1).sum()
    fig, axes = frame('Same saved distance-row PCA coordinates — baseline groups and QC', columns=2, size=(14, 6))
    for number, (label, group) in enumerate(table.groupby('struct_cluster', sort=True)):
        axes[0].scatter(group.PC1, group.PC2, s=18, color=GROUP_COLOURS[number], label=f"{label} (n={len(group)})", alpha=0.8)
    for label, group in table.groupby('qc_status', sort=True):
        axes[1].scatter(group.PC1, group.PC2, s=18, color=QC_COLOURS[label], label=f"{label} (n={len(group)})", alpha=0.8)
    for axis in axes:
        axis.set(xlabel=f"PC1 ({variance[0]:.1%} variance)", ylabel=f"PC2 ({variance[1]:.1%} variance)")
        axis.legend(fontsize=8)
    limits = (axes[0].get_xlim(), axes[0].get_ylim())
    axes[1].set_xlim(limits[0])
    axes[1].set_ylim(limits[1])
    save(fig, output, '09_landscape', 'Unscaled RMSD-row PCA, n=351 selected entries. Both panels use identical saved PC1/PC2 values, preserved display signs and shared limits. QC labels do not refit the embedding. PCA is descriptive geometry, not an independent biological validation.')
