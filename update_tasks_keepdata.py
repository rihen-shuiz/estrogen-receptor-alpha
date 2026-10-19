#!/usr/bin/env python3
"""Keep-all-data revision of esr1.tasks (professor notes, 2026-10-18).

Run from the repo root:  python3 update_tasks_keepdata.py
- makes tracker/esr1.tasks.bak3
- rewrites S03, S04, S05, S08, S13, S16 (new specs)
- adds K00..K09 (phase keep)
- logs REOPEN (who=human) for every structure task whose output will be rebuilt
Nothing else in the file is touched. Stops without changes if the file is not as expected.
"""
import re, shutil, subprocess, sys
from pathlib import Path

P = Path("tracker/esr1.tasks")
if not P.exists():
    sys.exit("Run from the repo root (tracker/esr1.tasks not found).")
text = P.read_text(encoding="utf-8")
if re.search(r"^T K00 ", text, re.M):
    sys.exit("Already applied (T K00 exists). Nothing changed.")

REPLACE = {}
REPLACE["S03"] = '''T S03 S run "Align EVERY eligible chain (not one per entry); log core RMSD per chain"
  help Apply S02 to all eligible ESR1 chains (the completeness screen stays exactly as it is) on the core ranges in esr_analysis/resources/config/alignment_core_v3.json (K04). Also fit with the old core (340-360, 400-420, 450-470, 500-520) and store both core RMSDs so the choice can be compared. Same reference as the last run; log its chain_key in a NOTE. A large core RMSD is labelled, never a reason to drop the chain.
  meta wk=3 who=B ph=core
  needs S02 Q04 K00 K04
  out results/alignment_core_rmsd.csv
  chk csv results/alignment_core_rmsd.csv rows>=700 cols=entry_id,chain_key,core_rmsd,core_rmsd_old_core,core_used nonempty=chain_key
  orig esr_analysis/align_structures.py
  src professor notes 2026-10-18 (keep as much data as possible; which helices are stable)'''
REPLACE["S04"] = '''T S04 S run "Pairwise Calpha RMSD matrix over ALL eligible chains (H12 included) - save to disk"
  help N x N RMSD over common LBD Calpha positions, H12 positions included, after the core fit. Rows are every eligible chain from S03, not one chain per entry; the file name keeps "entry" so downstream paths still work, and the row-order file (S05) says which chain is which. Keep min_pairwise_common_ca as in the settings (never lower it silently). A chain is left out of the matrix only if it cannot meet it, and then it is listed with the reason in results/alignment/distance_exclusions.json and in chain_status (K06), never silently. Save with np.save, no NaN, never recompute in notebooks. The previous 351-chain matrix is archived in results/archive_v2_strict/ (K00).
  meta wk=3 who=B ph=core
  needs S03
  tool numpy
  out results/rmsd_entry.npy
  chk npy results/rmsd_entry.npy square>=700
  orig esr_analysis/calculate_rmsd.py
  src desc 5.1 step 2; professor notes 2026-10-18'''
REPLACE["S05"] = '''T S05 XS code "Save matrix row order (index file, one row per chain)"
  help Ordered chain_key / entry_id list matching the matrix rows, so every analysis loads labels without guessing. Rows are chains; an entry with several chains appears several times.
  meta wk=3 who=B ph=core
  needs S04
  out results/rmsd_entry_index.csv
  chk csv results/rmsd_entry_index.csv rows>=700 cols=chain_key,entry_id nonempty=chain_key uniq=chain_key
  orig results/current/rmsd_entry_index.txt'''
REPLACE["S08"] = '''T S08 S code "Hierarchical clustering of all chains at cuts 1.0 / 1.5 / 2.0 / 2.5 A"
  help Average linkage on the RMSD matrix (rows = chains); one label column per cut, plus the linkage saved for the dendrogram. Every chain in the matrix gets a label; nothing is filtered before clustering. Copy qc_label (K05) next to the labels for reading only.
  meta wk=3 who=B ph=core
  needs S04 S05 K05
  tool scipy
  out results/entry_clusters.csv
  chk csv results/entry_clusters.csv rows>=700 cols=entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5,qc_label nonempty=chain_key
  orig esr_analysis/cluster_structures.py
  src desc 5.1 step 3'''
REPLACE["S13"] = '''T S13 S code "QC-tier sensitivity: all chains vs tier A only vs tiers A+B"
  help Replaces the old "drop R-free > 0.30" test. The main result stays on ALL chains. For each cut recluster (same method as S08, a slice of the saved S04 matrix, not recomputed) the chains with qc_label A only and with A or B, and tabulate clusters, top-3 share, clusters that have no member in the subset, and ARI against the all-chain labelling. Poor or failing chains stay in every file; this table only shows whether the conclusion depends on them.
  meta wk=6 who=B ph=core
  needs S08 Q06 K05
  out results/qc_clustering_delta.csv
  chk csv results/qc_clustering_delta.csv rows>=8 cols=cut,subset,n_chains,n_clusters,top3_fraction,clusters_without_members,ari_vs_all
  orig esr_analysis/compare_qc_clustering.py
  src desc 5.1 report list; professor notes 2026-10-18'''
REPLACE["S16"] = '''T S16 S code "Measure H12 position per chain after core fit (keep every chain)"
  help RMSD of H12 residues 538-547 versus the reference after the core fit; this is the agonist/antagonist signal. Keep every chain: no shared H12 Calpha -> h12_rmsd blank and h12_label not_modelled; fewer than 10 shared -> partial; otherwise complete. Add h12_mean_b and h12_b_group from K03 so a reader sees how mobile H12 is. Large values are reported next to n_h12_ca and never filtered.
  meta wk=4 who=B ph=core
  needs S03 K03
  out results/h12_metric.csv
  chk csv results/h12_metric.csv rows>=700 cols=entry_id,chain_key,n_h12_ca,h12_rmsd,h12_label,h12_mean_b,h12_b_group nonempty=chain_key,h12_label
  src desc 3.2, 5.1; professor notes 2026-10-18'''

NEW = '''# ===================== KEEP-ALL-DATA REVISION (phase keep) =====================
# Professor notes after the presentation: (1) which helices are most stable by B-factor group, use them for alignment;
# (2) label a structure you do not like, keep as much data as possible. H12 stays in. Bad data is labelled "bad", never removed.
# The completeness screen (LBD / core Calpha counts) is unchanged. Strict outputs are archived, not deleted.

T K00 XS run "Archive the strict (351-chain) outputs before the rewrite"
  help COPY, never move or delete: results/rmsd_entry.npy, rmsd_entry_index.csv, entry_clusters.csv, linkage_entry.npy, cluster_stats.csv, medoids.csv, qc_clustering_delta.csv, h12_metric.csv, alignment_core_rmsd.csv and figures/dendrogram.png, figures/structural_sensitivity.png into results/archive_v2_strict/ (same file names). Write MANIFEST.csv with the sha256 of each copy. These stay as the "strict" comparison.
  meta wk=7 who=B ph=keep
  out results/archive_v2_strict/MANIFEST.csv
  chk csv results/archive_v2_strict/MANIFEST.csv rows>=10 cols=file,sha256 nonempty=file,sha256 uniq=file
  src professor notes 2026-10-18

T K01 S run "Attrition funnel: every step that removes a chain, entry or residue"
  help Start from all downloaded ESR1 chains and list each step that removes or ignores something: completeness screen, apo/no-ligand, one-chain-per-entry selection (frozen index), distance pruning (min_pairwise_common_ca), QC flags, anything else you find in the code. One row per step: n before, n after, n removed, rule, code file, and necessary = yes only when it is identity or mathematics (wrong protein, no overlap), otherwise no. Also write one row per ESR1 chain with the stage where it was lost (empty if it reached the old matrix) and the reason. Numbers come from files, never from memory. Delete nothing.
  meta wk=7 who=A ph=keep
  needs Q04 Q06 K00
  out results/attrition_funnel.csv results/chain_loss_reasons.csv
  chk csv results/attrition_funnel.csv rows>=4 cols=step,n_before,n_after,n_removed,rule,code_file,necessary nonempty=step,rule
  chk csv results/chain_loss_reasons.csv rows>=788 cols=chain_key,entry_id,stage_lost,reason nonempty=chain_key,entry_id
  src professor notes 2026-10-18

T K02 S code "Helix table: residue range of each LBD helix, with its source"
  help One row per helix of the ERalpha LBD (305-550) from the mmCIF helix annotation (struct_conf) of the reference chain in docs/reference.md, or from DSSP if the binary exists; record which. Helix names the source does not give are unknown, never invented. H12 must come out as 538-547 as in the settings; if the source disagrees, report both and do not overwrite the settings.
  meta wk=7 who=B ph=keep
  needs S01
  out results/helix_table.csv
  chk csv results/helix_table.csv rows>=8 cols=helix,start,end,source,name_status nonempty=helix,start,end,source
  src professor note 1

T K03 M code "Normalised B-factor per helix for every eligible chain"
  help For every eligible ESR1 chain take the Calpha B-factors and divide each by that chain's own median Calpha B, so resolution and refinement differences cancel; average per helix from K02. Keep raw means too. Per chain also: core mean B, H12 mean B and the ratio H12 / core, and b_group = low / mid / high tertile of the chain's mean B. A chain missing a helix gets blank, not zero. Do not drop high-B chains.
  meta wk=7 who=B ph=keep
  needs K02 Q04
  out results/helix_bfactor_by_chain.csv
  chk csv results/helix_bfactor_by_chain.csv rows>=700 cols=chain_key,entry_id,b_group,mean_b_chain,h12_mean_b,core_mean_b,h12_to_core_ratio nonempty=chain_key,b_group
  src professor note 1

T K04 S code "Which helices are most stable per B-factor group; choose the alignment core"
  help Per b_group and over all chains: median and IQR of normalised B per helix, rank 1 = most stable, n chains. Report the H12 / core ratio (median and quartiles, raw and normalised) and say whether the team's earlier estimate, H12 about 12 times the core, holds; if not, state what was measured. Recommended core = helices that rank in the better half in every b_group, never H12, at least 60 residues in total; say if it equals the current core (340-360, 400-420, 450-470, 500-520) or how it differs. Save it as core_ranges in the config JSON that S03 reads. Figure: heatmap helix x b_group of median normalised B. Every number in the text comes from the CSV.
  meta wk=7 who=B ph=keep
  needs K03
  out results/helix_stability_by_bgroup.csv docs/helix_stability.md figures/helix_bfactor_stability.png esr_analysis/resources/config/alignment_core_v3.json
  chk csv results/helix_stability_by_bgroup.csv rows>=8 cols=helix,b_group,n_chains,median_norm_b,iqr_norm_b,rank nonempty=helix,b_group
  chk lines docs/helix_stability.md >=15
  chk grep docs/helix_stability.md (?i)recommended core
  chk grep esr_analysis/resources/config/alignment_core_v3.json core_ranges
  chk file figures/helix_bfactor_stability.png
  src professor note 1

T K05 S code "Graded QC tiers instead of a hard cutoff"
  help Tiers A / B / C for resolution, R-free and mean B per entry. A limits = the old flags (2.5 A, 0.30, 50). B limits are relaxed starting values in esr_analysis/resources/config/qc_tiers.json (resolution 3.5 A, R-free 0.40, mean B 80 A^2); the human may edit that file and rerun. qc_label: A = all metrics A; B = some B, no C; C = any C; unknown = a metric is missing (never a pass). Add a short plain-language warning per entry, for example "resolution 3.2 A (tier B)". No entry is removed; the old qc_* columns in results/entries.csv stay untouched.
  meta wk=7 who=A ph=keep
  needs Q06
  out results/qc_tiers.csv esr_analysis/resources/config/qc_tiers.json
  chk csv results/qc_tiers.csv rows>=450 cols=entry_id,tier_resolution,tier_rfree,tier_b,qc_label,warning nonempty=entry_id,qc_label
  chk grep esr_analysis/resources/config/qc_tiers.json resolution
  src professor note 2

T K06 S code "Master chain status table: every chain labelled, none dropped"
  help One row per ESR1 chain: completeness status and reason, qc_label and warning (K05), h12_b_group and h12_to_core_ratio (K03), in_matrix true/false, why it was lost in the old pipeline (K01). in_matrix is false only for a stated mathematical reason (no overlap), and the chain still has its row.
  meta wk=7 who=A ph=keep
  needs S05 K01 K03 K05
  out results/chain_status.csv
  chk csv results/chain_status.csv rows>=788 cols=chain_key,entry_id,completeness_status,qc_label,h12_b_group,in_matrix,warning nonempty=chain_key,entry_id,in_matrix
  src professor note 2

T K07 S code "Entry-level view next to the chain-level result"
  help Chains of one entry are not independent crystals. Pick one chain per entry by a stated rule (most mapped LBD Calpha, then lowest mean B, then chain_key order), slice those rows from the saved S04 matrix (do not recompute), cluster at the same four cuts, and compare with the chain-level clustering: clusters per cut, top-3 share, ARI. Keep both views; neither replaces the other.
  meta wk=7 who=B ph=keep
  needs S08 K06
  out results/entry_level_clusters.csv results/entry_vs_chain_clusters.csv
  chk csv results/entry_level_clusters.csv rows>=350 cols=entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5 nonempty=entry_id,chain_key uniq=entry_id
  chk csv results/entry_vs_chain_clusters.csv rows>=4 cols=cut,n_clusters_chain,n_clusters_entry,top3_chain,top3_entry,ari
  src professor notes 2026-10-18

T K09 M run "Alignment-core sensitivity: old core vs data-chosen core"
  help Rebuild the matrix with the old core fit (stored by S03) as a NEW file results/rmsd_old_core.npy (the primary matrix is not replaced), cluster it at 1.5 A, and compare with the primary result: cluster counts, ARI, and Spearman correlation of all pairwise RMSDs. Say whether the choice of core changes the clusters.
  meta wk=7 who=B ph=keep
  needs S03 S04 S08 K04
  out results/core_sensitivity.csv results/rmsd_old_core.npy
  chk csv results/core_sensitivity.csv rows>=1 cols=core,n_chains,n_clusters_1_5,ari_vs_primary,spearman_pairwise
  chk npy results/rmsd_old_core.npy square>=700
  src professor note 1

T K08 S text "Data-retention note for the professor: what was removed before, what is kept now"
  help Plain markdown, every number from a results file. Old strict set (351 chains) versus now; what each former filter removed (K01); QC tier counts and how many chains carry a warning (K05, K06); the H12 / core B ratio (K04); entry-level versus chain-level clusters (K07); effect of the strict subsets (S13). State plainly: nothing is removed for quality, chains are labelled; the only chains outside the matrix are listed with a mathematical reason. Answer note 1 (most stable helices per B-group, the core used) and note 2 in a few sentences each.
  meta wk=8 who=ALL ph=keep
  needs K01 K04 K06 K07 S13
  out docs/data_retention.md
  chk lines docs/data_retention.md >=25
  chk grep docs/data_retention.md (?i)label
  chk grep docs/data_retention.md (?i)stable
  src professor notes 2026-10-18
'''

lines = text.split("\n")
def block(tid):
    s = next((i for i,l in enumerate(lines) if l.startswith(f"T {tid} ")), None)
    if s is None: sys.exit(f"Task {tid} not found; file is not the expected 115-task version. Nothing changed.")
    e = s
    while e < len(lines) and lines[e].strip() != "": e += 1
    return s, e
for tid in REPLACE:
    pass
if sum(l.startswith("T ") for l in lines) != 115:
    sys.exit("Expected 115 tasks. Nothing changed.")
shutil.copy(P, str(P) + ".bak3")
for tid, new in REPLACE.items():
    s, e = block(tid)
    lines[s:e] = new.split("\n")
out = "\n".join(lines)
marker = "# ===================== STRUCTURE (phase core) ====================="
if marker not in out: sys.exit("Section marker not found. Nothing changed (backup exists).")
out = out.replace(marker, NEW + "\n" + marker, 1)
P.write_text(out, encoding="utf-8")
print("tasks file updated: 6 specs rewritten, 10 tasks added (K00-K09)")

py = sys.executable
note = "keep-all-data rewrite (professor notes 2026-10-18): rerun on all eligible chains; strict outputs archived by K00"
for tid in ["S03","S04","S05","S06","S08","S09","S10","S12","S13","S14","S16"]:
    r = subprocess.run([py, "tracker/tl.py", "log", tid, "REOPEN", "--who", "human", "--note", note], capture_output=True, text=True)
    print(tid, "REOPEN" if r.returncode == 0 else "FAILED: " + (r.stdout + r.stderr).strip()[:200])
print("Now run: python3 tracker/tl.py validate   (expect 0 errors) and python3 tracker/tl.py status")
