#!/usr/bin/env bash
# One command that reproduces the project from raw inputs (tracker Y39).
#   bash run_all.sh            # everything: upstream rebuild + check, original pipeline, extended analyses, documents
#   bash run_all.sh upstream   # data download + upstream rebuild into results/rebuild/ with a frozen-vs-rebuilt report
#   bash run_all.sh core       # original tracker pipeline (QC labels, core, superposition, RMSD, clusters, R/E/X tasks)
#   bash run_all.sh ext        # extended analyses of the final run (scripts/ext_*.py) -> results/ext, figures/ext
#   bash run_all.sh docs       # figure index, gallery, numbers and the generated documents
# Environment: python with requirements-lock.txt (PYTHON=... to choose the interpreter). Optional external tools:
# P2Rank (prank on PATH) for X08; PyMOL (PYMOL_PYTHON or the conda env, see scripts/pymol_safe.py) for F01/F02/Y30.
# Network is needed for the downloads (RCSB, ChEMBL); every downloaded file is checked or cached.
# Order follows the tracker dependencies (tracker/esr1.tasks 'needs'). Seeds: 42 everywhere.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"
SEP="$("$PY" -c "import os; print(os.pathsep)")"   # ";" on Windows, ":" elsewhere
export PYTHONPATH=".${SEP}scripts${PYTHONPATH:+${SEP}$PYTHONPATH}"
export ESR1_ROOT="$(pwd)"
stage="${1:-all}"
run() { echo ">>> $*"; "$PY" "$@"; }

upstream() {
  run scripts/fetch_data.py                 # W04: mmCIFs, CCD files, data/lbd (skips files already present)
  run scripts/rebuild_all.py                # Y31: Q04 Q06 Q08 C01-C07 C10 E03 E04 X01 X13 X14 -> results/rebuild + report
}

core() {
  run scripts/ligand_of_interest.py         # Q18
  run scripts/clean_chain_pdbs.py           # S17
  run scripts/chem_ligands_of_interest.py   # C11
  run scripts/attrition_funnel.py           # K01
  run scripts/helix_table.py                # K02
  run scripts/helix_bfactors.py             # K03
  run scripts/helix_stability.py            # K04 (writes alignment_core_v3.json)
  run scripts/qc_tiers.py                   # K05
  run scripts/align_selected_chains.py      # S03
  run scripts/build_selected_rmsd.py        # S04
  run scripts/save_rmsd_index.py            # S05
  run scripts/h12_metric.py                 # S16
  run scripts/cluster_selected.py           # S08
  run scripts/cluster_stats.py              # S09
  run scripts/figure_dendrogram.py          # S10
  run scripts/medoids.py                    # S12
  run scripts/chain_status.py               # K06
  run scripts/qc_clustering_delta.py        # S13
  run scripts/figure_sensitivity.py         # S14
  run scripts/allatom_rmsd_subset.py        # S07
  run scripts/pair_superposition.py         # E05
  run scripts/cluster_entropy.py            # R01
  run scripts/diversity_curve.py            # R02
  run scripts/saturation.py                 # R04
  run scripts/figure_diversity_curve.py     # R03
  run scripts/receptor_displacements.py     # E06 E07
  run scripts/figure_displacement_heatmap.py  # E08
  run scripts/consistent_positions.py       # E10
  run scripts/qc_outliers.py                # Q14
  run scripts/qc_bfactor_clusters.py        # Q15
  run scripts/core_sensitivity.py           # K09
  run scripts/entropy_qc_compare.py         # R06
  run scripts/entry_level_clusters.py       # K07
  run scripts/representative_set.py         # S15
  run scripts/esr2_rmsd.py                  # X02
  run scripts/esr2_clusters.py              # X03
  run scripts/esr2_pca.py                   # X04
  run scripts/esr2_basin_counts.py          # X05
  run scripts/potency_by_cluster.py         # X15 X16
  if command -v prank >/dev/null 2>&1; then run scripts/p2rank_batch.py; else echo "skip X08 (P2Rank not on PATH)"; fi
  run scripts/pocket_table.py               # X09
  run scripts/pocket_pca.py                 # X12
  if [ -n "${PYMOL:-}" ]; then run scripts/pymol_figures.py --pymol "$PYMOL"; else echo "skip F01/F02 render (PYMOL unset)"; fi
  run scripts/provisional_clusters.py       # P01 (September matrix)
  run scripts/provisional_saturation.py     # P02
  run scripts/provisional_qc_delta.py       # P03
  run scripts/provisional_chem_struct.py    # P04
}

ext() {
  for s in ext_metadata ext_embeddings ext_zoo ext_metric_check ext_compactness ext_agreement ext_dendrogram \
           ext_chem_zoo ext_ligand_descriptors ext_displacement_profile ext_structural_features ext_phipsi \
           ext_supervised ext_split_hypotheses ext_outliers ext_outlier_quality ext_metadata_tests ext_colour_by \
           ext_interactive ext_coord_pca ext_bfactors ext_esr2 ext_joint_esr1_esr2 ext_selectivity \
           ext_chem_struct_link ext_redundancy ext_time ext_robustness ext_esr2_core_fraction ext_pymol_renders \
           ext_ncs_copies ext_h12_contacts ext_cl4_substructure; do
    run "scripts/$s.py"
  done
}

docs() {
  run scripts/ext_build_index.py            # Y33 docs/figures/INDEX.md, docs/gallery.html
  run scripts/collect_numbers.py            # W16 + Y35 (numbers_ext.json, renders docs/templates/*.md)
}

case "$stage" in
  all) upstream; core; ext; docs ;;
  upstream) upstream ;;
  core) core ;;
  ext) ext ;;
  docs) docs ;;
  *) echo "usage: bash run_all.sh [all|upstream|core|ext|docs]"; exit 2 ;;
esac
echo "run_all.sh: stage '$stage' finished"
