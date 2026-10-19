from pathlib import Path
import importlib
import json
import os
import subprocess
import sys

out = Path('/tmp/esr_firstmonth_verify/trial')
out.mkdir(exist_ok=True)
names = ['query_structures','download_structures','extract_metadata','extract_lbd','flag_quality',
         'build_ligand_table','assign_action_classes','choose_reference','align_structures',
         'calculate_rmsd','cluster_structures','structural_pca','regional_variation']
env = dict(os.environ, PYTHONPATH='', PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR='/tmp/esr_firstmonth_verify/mpl',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
records = []
for number, name in enumerate(names, 1):
    module = importlib.import_module('esr_analysis.' + name)
    assert '/tmp/esr_firstmonth_verify/env/' in module.__file__, module.__file__
    help_result = subprocess.run([sys.executable, '-m', 'esr_analysis.' + name, '--help'], env=env, cwd='/tmp', capture_output=True, text=True)
    assert help_result.returncode == 0, help_result.stderr
    command = [sys.executable, '-m', 'esr_analysis.' + name, '--output', str(out / f'block_{number:02d}')]
    result = subprocess.run(command, env=env, cwd='/tmp', capture_output=True, text=True)
    (out / f'block_{number:02d}.log').write_text(result.stdout + result.stderr)
    records.append({'block': number, 'module': name, 'command': command, 'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr, 'help_returncode': help_result.returncode})
    (out / 'COMMAND_RESULTS.json').write_text(json.dumps(records, indent=2))
    print(number, name, result.returncode, result.stdout.strip(), flush=True)
    if result.returncode:
        print(result.stderr, flush=True)
        raise SystemExit(result.returncode)
command=[sys.executable,'-m','esr_analysis.report_figures','--output',str(out/'figures')]
result=subprocess.run(command,env=env,cwd='/tmp',capture_output=True,text=True)
(out/'figure_command.log').write_text(result.stdout+result.stderr)
print('Figures',result.returncode,result.stdout,result.stderr,flush=True)
raise SystemExit(result.returncode)
