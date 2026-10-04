"""Saved-input report figure: alignment."""
from ._common import json, np, Figure, FigureCanvasAgg, read_table, save

def plot(source, output, resources):
    from .._task_io import coordinate_mapping
    from ..align_structures import apply_transform
    mapping = read_table(source, 'residue_mapping.csv')
    transforms = json.loads((source / 'alignment_transforms.json').read_text())
    mobile = coordinate_mapping(mapping, '3Q95|B|1')
    reference = coordinate_mapping(mapping, '7NFB|B|1')
    settings = json.loads((resources / 'config/analysis_settings.json').read_text())
    core = {p for a, b in settings['esr1_core_ranges'] for p in range(a, b + 1)}
    shared = sorted(core & set(mobile) & set(reference))
    transform = transforms['3Q95|B|1']
    cm = np.asarray(transform['mobile_centroid'])
    cr = np.asarray(transform['reference_centroid'])
    rotation = np.asarray(transform['rotation'])
    xyz = np.array([mobile[p] for p in shared])
    ref = np.array([reference[p] for p in shared])
    stages = [(xyz, ref, 'Deposited coordinates'), (xyz - cm, ref - cr, 'Center shared core'), (apply_transform(xyz, rotation, cm, cr), ref, 'Saved conserved-core fit')]
    fig = Figure(figsize=(15, 5), layout='constrained')
    FigureCanvasAgg(fig)
    fig.suptitle(f"3Q95 B to 7NFB B — {len(shared)} shared core CA positions")
    for number, (mobile_xyz, ref_xyz, title) in enumerate(stages, 1):
        axis = fig.add_subplot(1, 3, number, projection='3d')
        axis.scatter(*ref_xyz.T, s=10, color='#999999', label='7NFB B')
        axis.scatter(*mobile_xyz.T, s=10, color='#4477aa', label='3Q95 B')
        axis.set(title=title, xlabel='x (Å)', ylabel='y (Å)', zlabel='z (Å)')
        axis.view_init(elev=20, azim=35)
        axis.legend(fontsize=8)
    save(fig, output, '06_alignment', 'Three actual coordinate stages from the saved whole-chain conserved-core transform; centering alone is not alignment. Atom correspondence is canonical/native mapping, not residue-index order or a new pairwise refit.')
