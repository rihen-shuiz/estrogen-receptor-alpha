"""Saved-input report figure: distances."""
from ._common import np, frame, save, aligned_annotations

def plot(source, output, resources):
    matrix = np.load(source / 'rmsd_entry.npy')
    labels = aligned_annotations(source).struct_cluster.to_numpy()
    order = np.argsort(labels, kind='stable')
    values = matrix[np.triu_indices(len(matrix), 1)]
    fig, axes = frame(f"Common-reference CA RMSD — {len(matrix)} selected entries, {len(values):,} distinct pairs", columns=2, size=(14, 6))
    image = axes[0].imshow(matrix[np.ix_(order, order)], cmap='viridis', vmin=0, aspect='equal')
    fig.colorbar(image, ax=axes[0], label='CA RMSD (Å)')
    axes[0].set(xlabel='Entry index sorted by baseline group', ylabel='Same ordered entry index')
    axes[1].hist(values, bins=60, color='#4477aa')
    axes[1].set(xlabel='CA RMSD (Å)', ylabel='Distinct entry pairs')
    save(fig, output, '07_distances', 'Symmetric indexed 351-entry matrix; diagonal excluded from the pair histogram. Shared observed LBD positions measured after the saved reference-core fit, without fitting pairs again. Sorting uses original baseline labels and stable within-group entry order.')
