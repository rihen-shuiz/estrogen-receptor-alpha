"""VIS1 per-image decisions (classification, cut parameters, visual check verdicts).

Keys are project-relative PNG paths (the first path of a SHA-1 group). Images under figures/sv/ are SV figures
(method a, listed from figures/panels/sv/<ID>/) and are not keyed here except for check verdicts.

Classification follows the VIS1 definition of a graph (one plot with its axes, ticks, axis titles, own legend, own
colorbar and attached marginal axes). Declared reading of that definition, fixed before cutting:
- Axes stacked on one shared x axis whose x tick labels are drawn only under the lowest axes form ONE graph (the
  upper axes cannot be read without the shared axis). Stacked axes that each carry their own x tick labels are
  separate graphs; the shared x-axis title then stays only under the lowest one (noted in the manifest).
- A key that sits next to a graph and is needed to read it (e.g. scaffold drawings for bar labels S1-S12) stays
  with that graph (rule 5).
- A PyMOL render is one graph; a figure holding two renders is multi.

mode:
- "svg"    : cut positions chosen from the PNG's white gutters, guided by the axes layout read from the SVG with
             the same stem (top-level <g id="axes_N"> rectangles + text anchors); panels are crops of the PNG.
- only     : allowed cut directions, declared from the layout ("v" = one row of graphs, "h" = one column).
- top_title: False when the top strip holds the graphs' own titles (no suptitle), so it must not be dropped.
- cut_x    : declared vertical cut (px) where the widest gutter is not the one between graphs (07_distances: the
             colorbar of the heatmap and the histogram's y labels share a column band; cut in the background
             column run 1285-1311 between the colorbar label and the histogram's y-axis title).
- "gutter" : plain recursive XY-cut of the PNG along background gutters (no SVG used); a cut is taken only if
             both sides keep a graph-sized part (>= 12 % of the image height and >= 8 % of its width).
check: visual verdict from figures/panels/_check/ (filled after looking at every sheet).
"""

SINGLE = {
    "esr_analysis/resources/reference/reference_context.png": "PyMOL render with its caption text",
    "figures/01_core_windows.png": "one axes",
    "figures/08_group_counts.png": "one axes",
    "figures/10_regional_B.png": "one axes",
    "figures/dendrogram.png": "one dendrogram",
    "figures/deposits_per_year.png": "one axes",
    "figures/displacement_heatmap.png": "heatmap with its top bars and colorbar (attached)",
    "figures/displacement_structure.png": "PyMOL render",
    "figures/diversity_curve.png": "one axes",
    "figures/esr2_pca_landscape.png": "one axes",
    "figures/ext/Y06_alluvial.png": "one alluvial plot with its legend",
    "figures/ext/Y07_dendrogram_metadata.png": "dendrogram with attached annotation bars and their legends",
    "figures/ext/Y07_rmsd_heatmap_annotated.png": "heatmap with attached side bars, legend and colorbar",
    "figures/ext/Y09_scaffolds.png": "bar chart; the 12 scaffold drawings are the key to its S1-S12 labels (rule 5)",
    "figures/ext/Y10_displacement_profile.png": "3 axes stacked on one shared x axis, tick labels only at the bottom",
    "figures/ext/Y13_decision_tree.png": "one tree diagram",
    "figures/ext/Y20_coord_pca_loadings.png": "3 axes stacked on one shared x axis, tick labels only at the bottom",
    "figures/ext/Y23_receptor_residue_diff.png": "profile + flag track on one shared x axis, tick labels only at the bottom",
    "figures/ext/Y27_timeline.png": "3 axes stacked on one shared x axis, tick labels only at the bottom",
    "figures/ext/pymol_panels/bputty_cl4.png": "PyMOL render",
    "figures/ext/pymol_panels/bputty_cl7.png": "PyMOL render",
    "figures/ext/pymol_panels/medoids_front.png": "PyMOL render",
    "figures/ext/pymol_panels/medoids_side.png": "PyMOL render",
    "figures/helix_bfactor_stability.png": "heatmap with its colorbar",
    "figures/medoids.png": "PyMOL render",
    "figures/pocket_pca.png": "one axes",
    "figures/potency_by_cluster.png": "one axes",
    "figures/provisional_diversity_curve.png": "one axes",
    "reference_render/reference_context.png": "PyMOL render with its caption text",
    "results/archive_v2_strict/dendrogram.png": "one dendrogram",
    "results/figures/01_core_windows.png": "one axes",
    "results/figures/08_group_counts.png": "one axes",
    "results/figures/10_regional_B.png": "one axes",
    "results/reference_render/reference_context.png": "PyMOL render with its caption text",
}

# multi images: expected number of graphs (from looking at the originals), mode and parameters.
# groups: extra axes groupings for svg mode (thin axes such as colorbars are attached to their neighbour
# automatically). titles: graph title or axis label per panel for gutter mode (read from the PNG), in reading order; slug = slugified.
MULTI = {
    # ---- PNG only (older pipeline figures) ----
    "figures/02_inventory.png": dict(n=2, mode="gutter", only="v",
                                     titles=["Distinct deposited entries", "Selected entries (n=351); categories overlap"]),
    "figures/03_quality.png": dict(n=3, mode="gutter", only="v",
                                   titles=["Resolution (Å)", "Depositor R-free", "Mean B (Å²)"]),
    "figures/04_identity_action.png": dict(n=2, mode="gutter", only="v",
                                           titles=["ESR1 identity population", "Exact identities; current evidence-support policy"]),
    "figures/06_alignment.png": dict(n=3, mode="gutter", only="v",
                                     titles=["Deposited coordinates", "Center shared core", "Saved conserved-core fit"]),
    "figures/07_distances.png": dict(n=2, mode="gutter", only="v", cut_x=[1298], titles=["Entry index sorted by baseline group", "CA RMSD (Å)"]),
    "figures/09_landscape.png": dict(n=2, mode="gutter", only="v",
                                     titles=["PC1 (90.6% variance)", "PC1 (90.6% variance)"],
                                     panel_notes={1: "legend: baseline groups 1-9", 2: "legend: QC fail / pass"}),
    "figures/11_rmsf_coverage.png": dict(n=2, mode="gutter", only="h",
                                         titles=["Population ensemble RMSF (Å)", "Observed entry frames"],
                                         note="stacked axes, each with its own x tick labels; the shared x-axis "
                                              "title is printed only under the lower panel"),
    "results/figures/02_inventory.png": dict(n=2, mode="gutter", only="v",
                                             titles=["Distinct deposited entries", "Selected entries (n=351); categories overlap"]),
    "results/figures/03_quality.png": dict(n=3, mode="gutter", only="v",
                                           titles=["Resolution (Å)", "Depositor R-free", "Mean B (Å²)"]),
    "results/figures/04_identity_action.png": dict(n=2, mode="gutter", only="v",
                                                   titles=["ESR1 identity population",
                                                          "Exact identities; current evidence-support policy"]),
    "results/figures/06_alignment.png": dict(n=3, mode="gutter", only="v",
                                             titles=["Deposited coordinates", "Center shared core",
                                                    "Saved conserved-core fit"]),
    "results/figures/07_distances.png": dict(n=2, mode="gutter", only="v", cut_x=[1298],
                                             titles=["Entry index sorted by baseline group", "CA RMSD (Å)"]),
    "results/figures/09_landscape.png": dict(n=2, mode="gutter", only="v",
                                             titles=["PC1 (90.6% variance)", "PC1 (90.6% variance)"],
                                     panel_notes={1: "legend: baseline groups 1-9", 2: "legend: QC fail / pass"}),
    "results/figures/11_rmsf_coverage.png": dict(n=2, mode="gutter", only="h",
                                                 titles=["Population ensemble RMSF (Å)", "Observed entry frames"],
                                                 note="stacked axes, each with its own x tick labels; the shared "
                                                      "x-axis title is printed only under the lower panel"),
    "figures/provisional_threshold_sensitivity.png": dict(n=3, mode="gutter", only="v", top_title=False,
                                                          titles=["Clusters (count)", "Top-3 share of entries",
                                                                 "Entropy (bits)"]),
    "figures/qc_outliers.png": dict(n=2, mode="gutter", only="v", titles=["Resolution (Å)", "R-free"]),
    "figures/structural_sensitivity.png": dict(n=2, mode="gutter", only="v",
                                               titles=["Clusters", "Chains in the 3 largest clusters (%)"]),
    "results/archive_v2_strict/structural_sensitivity.png": dict(n=2, mode="gutter", only="v",
                                                                 titles=["Clusters",
                                                                        "Entries in the 3 largest clusters (%)"]),
    # PyMOL composites: the SVG holds the renders as images, so the plain gutter cut is used
    "figures/ext/Y30_pymol_bputty.png": dict(n=2, mode="gutter", only="v", legend="bottom",
                                             titles=["cl7 medoid 5TM1|A|1", "cl4 medoid 7UJ7|C|1"]),
    "figures/ext/Y30_pymol_medoids.png": dict(n=2, mode="gutter", only="v", legend="bottom",
                                              titles=["front view", "side view (rotated 90 deg about y)"]),
    # ---- matplotlib figures with an SVG twin ----
    "figures/ext/Y02_metadata_overview.png": dict(n=6, mode="svg"),
    "figures/ext/Y03_embeddings.png": dict(n=5, mode="svg"),
    "figures/ext/Y04_zoo_grid.png": dict(n=15, mode="svg"),
    "figures/ext/Y04_zoo_metrics.png": dict(n=3, mode="svg"),
    "figures/ext/Y05_compactness.png": dict(n=3, mode="svg"),
    "figures/ext/Y06_agreement.png": dict(n=3, mode="svg"),
    "figures/ext/Y08_chem_map.png": dict(n=2, mode="svg"),
    "figures/ext/Y08_chem_zoo_metrics.png": dict(n=3, mode="svg"),
    "figures/ext/Y09_physchem_by_cluster.png": dict(n=10, mode="svg"),
    "figures/ext/Y11_contact_map_diff.png": dict(n=3, mode="svg"),
    "figures/ext/Y11_features_by_cluster.png": dict(n=12, mode="svg"),
    "figures/ext/Y12_phipsi_pca.png": dict(n=3, mode="svg",
                                           note="right-hand axes stacked on one x axis but each has its own x tick "
                                                "labels; the x-axis title is printed only under the lower one"),
    "figures/ext/Y13_feature_importance.png": dict(n=3, mode="svg"),
    "figures/ext/Y14_split_hypotheses.png": dict(n=3, mode="svg"),
    "figures/ext/Y15_outliers.png": dict(n=8, mode="svg"),
    "figures/ext/Y16_outlier_quality.png": dict(n=6, mode="svg"),
    "figures/ext/Y17_metadata_effects.png": dict(n=3, mode="svg"),
    "figures/ext/Y18_colour_by_biology.png": dict(n=12, mode="svg"),
    "figures/ext/Y18_colour_by_experiment.png": dict(n=12, mode="svg"),
    "figures/ext/Y18_colour_by_pca.png": dict(n=8, mode="svg"),
    "figures/ext/Y18_colour_by_quality.png": dict(n=12, mode="svg"),
    "figures/ext/Y20_coord_pca.png": dict(n=4, mode="svg"),
    "figures/ext/Y21_core_selection.png": dict(n=2, mode="svg"),
    "figures/ext/Y21_helix_b.png": dict(n=3, mode="svg"),
    "figures/ext/Y22_esr2_parallel.png": dict(n=6, mode="svg"),
    "figures/ext/Y23_joint_embedding.png": dict(n=2, mode="svg"),
    "figures/ext/Y24_potency_esr1_vs_esr2.png": dict(n=3, mode="svg"),
    "figures/ext/Y25_chem_struct_link.png": dict(n=3, mode="svg"),
    "figures/ext/Y25_nn_agreement.png": dict(n=2, mode="svg"),
    "figures/ext/Y26_hill_numbers.png": dict(n=2, mode="svg"),
    "figures/ext/Y26_rarefaction.png": dict(n=2, mode="svg"),
    "figures/ext/Y26_set_cover.png": dict(n=2, mode="svg"),
    "figures/ext/Y27_novelty.png": dict(n=2, mode="svg"),
    "figures/ext/Y28_robustness_grid.png": dict(n=2, mode="svg"),
    "figures/ext/Y42_metric_check.png": dict(n=3, mode="svg"),
    "figures/ext/Y42_region_rmsd.png": dict(n=4, mode="svg"),
    "figures/ext/Y43_ncs_copies.png": dict(n=3, mode="svg"),
    "figures/ext/Y44_h12_contacts.png": dict(n=3, mode="svg"),
    "figures/ext/Y45_cl4_substructure.png": dict(n=3, mode="svg"),
}

# SV figures whose SV-script panels fail the VIS1 check (fragments of a neighbour or a cut-off label at a panel
# edge, seen on the contact sheets and confirmed by the edge-ink check): cut from the original PNG instead (method c).
SV_RECUT = {
    "figures/sv/SV01a/SV01a_forest_completeness.png": dict(
        n=4, mode="svg", reason="the forest_cl4 panel holds the end of the cl7 forest's title, x tick label and "
                                "legend at its left edge, forest_cl7 a fragment of the cl4 tick labels at its right"),
    "figures/sv/SV01b/SV01b_forest_completeness.png": dict(
        n=4, mode="svg", reason="the forest_cl4 panel holds fragments of the cl7 forest at its left edge, "
                                "forest_cl7 fragments of the cl4 tick labels at its right"),
    "figures/sv/SV01b/SV01b_forest_mobility.png": dict(
        n=4, mode="svg", reason="the forest_cl4 panel holds fragments of the cl7 forest at its left edge, "
                                "forest_cl7 fragments of the cl4 tick labels at its right"),
    "figures/sv/SV04/h12_map_partial.png": dict(
        n=2, mode="svg", reason="the legend of panel b is cut off at its right edge ('(n=' without the count)"),
    "figures/sv/SV05/SV05_core_v4.png": dict(
        n=5, mode="svg", reason="the y-axis title of panel c is cut off at its top edge"),
    "figures/sv/SV06/SV06_outliers_map.png": dict(
        n=2, mode="svg", reason="panel A holds the y-axis title of panel B at its right edge"),
    "figures/sv/SV14/SV14_filling.png": dict(
        n=2, mode="svg", reason="by_cluster holds the y-axis title and tick labels of by_class at its right edge, "
                                "by_class the end of the by_cluster title at its left"),
}

# SV figures that are one graph under the declared reading (shared axis whose tick labels are drawn on one axes only;
# the same rule as Y10 / Y20 loadings / Y23 residue diff / Y27 timeline), or single-panel SV figures with a defect.
SV_SINGLE = {
    "figures/sv/SV02/SV02_helix_track.png": "helix-fraction track and observing track on one shared x axis, x tick "
                                            "labels only under the lower track (the SV helix_fraction panel has none)",
    "figures/sv/SV15/SV15_dpx_by_group.png": "by-cluster and by-class strips on one shared y axis, y tick labels only "
                                             "on the left strip (the SV by_class panel has none); one marker legend",
}

# Visual check verdicts per contact sheet (key = image path; SV figures included). "pass" or "fail: <reason>".
# Every sheet in figures/panels/_check/ was opened and looked at; cut figures not listed here were not checked.
_PASS = [
    "figures/02_inventory.png", "figures/03_quality.png", "figures/04_identity_action.png",
    "figures/06_alignment.png", "figures/07_distances.png", "figures/09_landscape.png",
    "figures/11_rmsf_coverage.png", "figures/provisional_threshold_sensitivity.png", "figures/qc_outliers.png",
    "figures/structural_sensitivity.png", "results/archive_v2_strict/structural_sensitivity.png",
    "results/figures/02_inventory.png", "results/figures/03_quality.png", "results/figures/04_identity_action.png",
    "results/figures/06_alignment.png", "results/figures/07_distances.png", "results/figures/09_landscape.png",
    "results/figures/11_rmsf_coverage.png",
    "figures/ext/Y02_metadata_overview.png", "figures/ext/Y03_embeddings.png", "figures/ext/Y04_zoo_grid.png",
    "figures/ext/Y04_zoo_metrics.png", "figures/ext/Y05_compactness.png", "figures/ext/Y06_agreement.png",
    "figures/ext/Y08_chem_map.png", "figures/ext/Y08_chem_zoo_metrics.png", "figures/ext/Y11_contact_map_diff.png",
    "figures/ext/Y11_features_by_cluster.png", "figures/ext/Y12_phipsi_pca.png",
    "figures/ext/Y13_feature_importance.png", "figures/ext/Y14_split_hypotheses.png", "figures/ext/Y15_outliers.png",
    "figures/ext/Y17_metadata_effects.png", "figures/ext/Y18_colour_by_biology.png",
    "figures/ext/Y18_colour_by_experiment.png", "figures/ext/Y18_colour_by_pca.png",
    "figures/ext/Y18_colour_by_quality.png", "figures/ext/Y20_coord_pca.png", "figures/ext/Y21_core_selection.png",
    "figures/ext/Y22_esr2_parallel.png", "figures/ext/Y23_joint_embedding.png",
    "figures/ext/Y24_potency_esr1_vs_esr2.png", "figures/ext/Y25_chem_struct_link.png",
    "figures/ext/Y25_nn_agreement.png", "figures/ext/Y26_hill_numbers.png", "figures/ext/Y26_rarefaction.png",
    "figures/ext/Y27_novelty.png", "figures/ext/Y28_robustness_grid.png", "figures/ext/Y30_pymol_bputty.png",
    "figures/ext/Y30_pymol_medoids.png", "figures/ext/Y42_metric_check.png", "figures/ext/Y42_region_rmsd.png",
    "figures/ext/Y43_ncs_copies.png", "figures/ext/Y44_h12_contacts.png", "figures/ext/Y45_cl4_substructure.png",
    # SV figures whose SV-script panels pass (method a)
    "figures/sv/SV03/SV03_descriptors_by_cluster.png", "figures/sv/SV07/SV07_spike_profile.png",
    "figures/sv/SV13/SV13_forest.png", "figures/sv/SV13/SV13_h12_class_panels.png",
    # SV figures re-cut from the PNG (SV_RECUT) whose PNG panels pass
    "figures/sv/SV04/h12_map_partial.png", "figures/sv/SV05/SV05_core_v4.png",
]
CHECK = {k: "pass" for k in _PASS}
CHECK["figures/ext/Y21_helix_b.png"] = (
    "fail: in the original the colorbar label of the heatmap ('... 1 = most stable') runs into the figure suptitle "
    "(overprinted on 'median'), so any rectangle that holds the whole colorbar label also holds a cut-off piece of "
    "the suptitle")
