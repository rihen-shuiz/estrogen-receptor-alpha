"""Saved-input report figure: core_windows."""
from ._common import json, frame, save

def plot(source, output, resources):
    from matplotlib.patches import Rectangle
    settings = json.loads((resources / 'config/analysis_settings.json').read_text())
    fig, axes = frame('Configured ESR1 fitting windows — blocks 04/09', size=(13, 4))
    axis = axes[0]
    start, end = (settings['esr1_lbd_start'], settings['esr1_lbd_end'])
    axis.add_patch(Rectangle((start, 0.1), end - start + 1, 0.45, color='#bbbbbb', label='LBD envelope'))
    for number, (a, b) in enumerate(settings['esr1_core_ranges']):
        axis.add_patch(Rectangle((a, 0.7), b - a + 1, 0.45, color='#4477aa', label='Conserved fitting core' if number == 0 else None))
        axis.text((a + b) / 2, 1.23, f"{a}–{b}", ha='center')
    a, b = settings['h12_range']
    axis.add_patch(Rectangle((a, 1.65), b - a + 1, 0.45, color='#cc3311', label='Configured H12 analysis window'))
    axis.text((a + b) / 2, 2.16, f"{a}–{b}", ha='center')
    axis.set(xlim=(start - 10, end + 10), ylim=(-0.1, 2.8), xlabel='ESR1 canonical residue', yticks=[])
    axis.legend(loc='upper left', fontsize=9)
    save(fig, output, '01_core_windows', 'Project-defined LBD/core/H12 ranges; 84 configured core positions. The provisional regional H12* 536–545 interval is a distinct annotation, not this configured 538–547 window.')
