"""Block 12: PCA of unscaled RMSD rows in the common reference frame.
Input: indexed finite distance matrix. Output: two coordinates per row and
variance fractions. Preserve singleton/all-zero conventions and sklearn PCA."""
import numpy as np
from sklearn.decomposition import PCA
from .distance_matrix import validate_distance_matrix


def structural_pca(matrix):
    matrix = validate_distance_matrix(matrix)
    ids = range(len(matrix))
    if len(ids) == 1 or not np.any(matrix):
        coordinates = np.zeros((len(ids), 2), dtype=float)
        variance = np.zeros(2, dtype=float)
    else:
        pca = PCA(n_components=2)
        coordinates = pca.fit_transform(matrix)
        variance = pca.explained_variance_ratio_
    return coordinates, variance



def main():
    from ._task_io import parser, prepare, json_write
    from .run_io import read_table
    args = parser('Block 12: unscaled distance-row PCA; preserve frozen display signs').parse_args()
    root, source, out, settings = prepare(args)
    ids = (source / 'rmsd_entry_index.txt').read_text().splitlines()
    coordinates, variance = structural_pca(np.load(source / 'rmsd_entry.npy'))
    frozen = read_table(source, 'entry_embedding.csv').set_index('chain_key').loc[ids]
    for component, column in enumerate(['PC1', 'PC2']):
        if np.dot(coordinates[:, component], frozen[column]) < 0:
            coordinates[:, component] *= -1
    embedding = frozen[['struct_cluster']].reset_index()
    embedding[['PC1', 'PC2']] = coordinates
    embedding.to_csv(out / 'entry_embedding.csv', index=False)
    json_write(out / 'PCA_variance.json', {'explained_variance_ratio': variance.tolist(), 'features': 'unscaled RMSD rows', 'display_orientation': 'signs anchored to bundled baseline; no rotation/refit of QC panel'})
    print('Block 12: variance ' + str(variance.tolist()))


if __name__ == "__main__":
    main()
