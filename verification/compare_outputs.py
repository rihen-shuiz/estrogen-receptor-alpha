from pathlib import Path
import json
import numpy as np
import pandas as pd
from esr_analysis._task_io import RESOURCE_ROOT
from esr_analysis.run_io import read_table

source=RESOURCE_ROOT/'frozen'
out=Path('/tmp/esr_firstmonth_verify/trial')
checks=[]
def equal(name,actual,expected):
    assert actual==expected,(name,actual,expected)
    checks.append({'check':name,'status':'PASS','actual':actual,'expected':expected})
def arrays(name,actual,expected,atol=1e-12):
    assert np.allclose(actual,expected,atol=atol,rtol=0,equal_nan=True),name
    finite=np.isfinite(actual)&np.isfinite(expected)
    difference=float(np.max(np.abs(actual[finite]-expected[finite]))) if finite.any() else 0.
    checks.append({'check':name,'status':'PASS','max_absolute_difference':difference,'atol':atol,'rtol':0})
equal('01 ordered frozen cohort',(out/'block_01/pdb_ids_esr1.txt').read_text().splitlines(),(RESOURCE_ROOT/'data/pdb_ids_esr1.txt').read_text().splitlines())
raw=read_table(out/'block_03','entries.csv').set_index('pdb_id')
baseline=read_table(source,'entries.csv').set_index('pdb_id').loc[raw.index]
arrays('03 raw deposited metadata',raw[['resolution_A','r_free','mean_B_factor_A2']].to_numpy(float),baseline[['resolution_A','r_free','mean_B_factor_A2']].to_numpy(float))
rawchains=read_table(out/'block_04','receptor_chains.csv').set_index('chain_key')
savedchains=read_table(source,'receptor_chains.csv').set_index('chain_key').loc[rawchains.index]
for column in ['identity','mapped_lbd_ca','mapped_core_ca','mapped_h12_ca','eligible_historical','auth_asym_id']:
    equal('04 raw '+column,rawchains[column].tolist(),savedchains[column].tolist())
quality=read_table(out/'block_05','entries.csv').set_index('pdb_id')
equal('05 recomputed master QC',quality.qc_status.to_dict(),read_table(source,'entries.csv').set_index('pdb_id').qc_status.to_dict())
fresh=read_table(out/'block_06','ligand_identities.csv').set_index(['ligand_code','receptor_context'])
old=read_table(source,'ligand_identities.csv').set_index(['ligand_code','receptor_context']).loc[fresh.index]
equal('06 exact identity population',len(fresh),400)
equal('06 exact identities/stereochemical SMILES',fresh.chemical_identity_key.tolist(),old.chemical_identity_key.tolist())
actions=read_table(out/'block_07','ligand_actions.csv').set_index(['ligand_code','receptor_context']).loc[fresh.index]
equal('07 supported labels',actions.supported_action_type.fillna('unknown').tolist(),old.supported_action_type.fillna('unknown').tolist())
equal('07 evidence support',actions.evidence_supported.tolist(),old.evidence_supported.tolist())
equal('08 reference',json.loads((out/'block_08/reference.json').read_text())['chain_key'],'7NFB|B|1')
new_transforms=json.loads((out/'block_09/alignment_transforms.json').read_text())
old_transforms=json.loads((source/'alignment_transforms.json').read_text())
for key,transform in new_transforms.items():
    equal('09 shared core '+key,transform['n_common_core_ca'],84)
    for column in ['rotation','mobile_centroid','reference_centroid']:
        arrays('09 saved transform '+key+' '+column,np.asarray(transform[column]),np.asarray(old_transforms[key][column]),1e-10)
equal('10 retained matrix order',(out/'block_10/rmsd_entry_index.txt').read_text().splitlines(),(source/'rmsd_entry_index.txt').read_text().splitlines())
arrays('10 reconstructed full RMSD',np.load(out/'block_10/rmsd_entry.npy'),np.load(source/'rmsd_entry.npy'))
arrays('10 overlap',np.load(out/'block_10/overlap_entry.npy'),np.load(source/'overlap_entry.npy'),0)
equal('11 group labels',read_table(out/'block_11','entry_clusters.csv')['cluster_1.5A'].tolist(),read_table(source,'entry_clusters.csv')['cluster_1.5A'].tolist())
arrays('12 PCA coordinates/signs',read_table(out/'block_12','entry_embedding.csv')[['PC1','PC2']].to_numpy(),read_table(source,'entry_embedding.csv')[['PC1','PC2']].to_numpy(),1e-8)
arrays('13 ensemble profiles',read_table(out/'block_13','ensemble_profiles.csv').to_numpy(),read_table(source,'ensemble_profiles.csv').to_numpy(),1e-10)
new=read_table(out/'block_13','helix_chain_summaries.csv').set_index(['chain_key','interval']).sort_index()
old=read_table(source,'helix_chain_summaries.csv').set_index(['chain_key','interval']).sort_index()
equal('13 regional row/index correspondence',list(map(list,new.index.tolist())),list(map(list,old.index.tolist())))
numeric=['valid_B_CA','B_coverage','median_B_A2','mean_B_A2','median_B_z','matched_B_displacement_CA','matched_median_B_z','reference_RMSD_A','shape_coverage']
arrays('13 per-chain regional calculations',new[numeric].to_numpy(),old[numeric].to_numpy(),1e-10)
old_distribution=read_table(source,'helix_distribution_and_coverage.csv')
old_distribution=old_distribution.loc[old_distribution.coverage_threshold.eq(.8)].set_index(['interval','cohort']).sort_index()
new_distribution=read_table(out/'block_13','helix_distribution_and_coverage.csv').set_index(['interval','cohort']).sort_index()
arrays('13 coverage-aware distributions',new_distribution[['retained_chains','median_chain_median_z']].to_numpy(),old_distribution[['retained_chains','median_chain_median_z']].to_numpy(),1e-10)
summary=json.loads((out/'block_13/local_summary.json').read_text())
equal('13 experimental frames',summary['frames'],351)
equal('13 common positions',summary['common_atoms_count'],172)
checks.append({'check':'MDAnalysis new common-atom check','status':'SKIPPED','reason':summary['mdanalysis']})
checks.append({'check':'Exact seven-page report PDF/recipe visual matching','status':'SKIPPED','reason':'Named reference files absent; user explicitly directed proceeding with local folder contents'})
(out/'NUMERICAL_COMPARISONS.json').write_text(json.dumps(checks,indent=2)+'\n')
print('PASS numerical comparisons:',sum(c['status']=='PASS' for c in checks),'SKIPPED:',sum(c['status']=='SKIPPED' for c in checks))
