"""Reproduce the actual deposited 7NFB reference scene using explicit PyMOL."""
import argparse
import json
import os
from pathlib import Path
import subprocess
from ._task_io import RESOURCE_ROOT
from .pymol_scene import _PYMOL_SCENE


def render_reference(input_root, output, pymol):
    root, out = Path(input_root).resolve(), Path(output).resolve()
    if out.is_relative_to(root):
        raise ValueError('Render output must be separate from bundled inputs')
    out.mkdir(parents=True, exist_ok=True)
    cif = root / 'data/raw/esr1/7NFB.cif'
    if not cif.is_file():
        raise FileNotFoundError(f'Actual reference mmCIF missing: {cif}')
    scene = json.loads((root / 'reference/reference_scene.json').read_text())
    bootstrap = '''
import os
from pathlib import Path
output = Path(os.environ['ESR1_REFERENCE_OUTPUT']).resolve()
output.mkdir(parents=True, exist_ok=True)
for molecule in scene['molecules']:
    molecule['path'] = str(Path(os.environ['ESR1_REFERENCE_CIF']).resolve())
scene['png'] = str(output / 'reference_context.png')
scene['details'] = str(output / 'reference_context.details.json')
scene['error'] = str(output / 'reference_context.error.txt')
'''
    recipe = out / 'reference_context.pml'
    recipe.write_text('# Relocatable: supply ESR1_REFERENCE_CIF and ESR1_REFERENCE_OUTPUT.\npython\nscene = '
                      + repr(scene) + '\n' + bootstrap + _PYMOL_SCENE + '\npython end\n')
    env = dict(os.environ, ESR1_REFERENCE_CIF=str(cif), ESR1_REFERENCE_OUTPUT=str(out))
    result = subprocess.run([str(Path(pymol).resolve()), '-cq', str(recipe)], env=env,
                            capture_output=True, text=True, timeout=180)
    (out / 'reference_render.log').write_text(result.stdout + result.stderr)
    if result.returncode or not (out / 'reference_context.details.json').is_file():
        raise RuntimeError('Actual PyMOL rendering failed; inspect reference_render.log and error TXT')
    details = json.loads((out / 'reference_context.details.json').read_text())
    if not details['whole_molecule_transform_checked']:
        raise RuntimeError('Saved whole-object transform was not verified')
    return details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=RESOURCE_ROOT)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pymol', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(render_reference(args.input_root, args.output, args.pymol), indent=2))


if __name__ == '__main__':
    main()
