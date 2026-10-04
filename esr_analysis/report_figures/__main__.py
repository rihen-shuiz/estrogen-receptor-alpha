"""Thin offline figure orchestration; no extraction, fitting or live retrieval."""
import argparse
from importlib import import_module
from pathlib import Path
from .._task_io import RESOURCE_ROOT
NAMES = ['core_windows', 'inventory', 'quality', 'identity_action', 'reference_context', 'alignment', 'distances', 'group_counts', 'landscape', 'regional_B', 'rmsf_coverage']

def main():
    parser = argparse.ArgumentParser(description='Eleven first-month figures from declared frozen inputs only')
    parser.add_argument('--input-root', type=Path, default=RESOURCE_ROOT)
    parser.add_argument('--source-folder', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-image', type=Path)
    args = parser.parse_args()
    root = args.input_root.resolve()
    source = (args.source_folder or root / 'frozen').resolve()
    output = args.output.resolve()
    if output.is_relative_to(root) or output == source:
        raise ValueError('Figure output must be separate from inputs')
    for name in NAMES:
        module = import_module('.' + name, __package__)
        if name == 'reference_context':
            module.plot(source, output, root, args.reference_image)
        else:
            module.plot(source, output, root)
        print('Saved ' + name)
if __name__ == '__main__':
    main()
