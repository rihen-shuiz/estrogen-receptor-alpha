# Helix stability and the alignment core

This note asks which LBD helices of the reference ERalpha structure (7NFB, chain B, model 1) have the most stable C-alpha B-factors, and which helices should form the alignment core. It covers all 788 eligible ESR1 chains. No chain was removed.

## Data and method

The helix ranges come from the mmCIF helix annotation (struct_conf, HELX_P) of 7NFB chain B. The source gives helix IDs only, so every helix is named "unknown" except the one that covers the H12 range in the settings, which is named H12 (results/helix_table.csv). The source places H12 at 537-546, while the settings place it at 538-547. Both ranges are reported. The settings were not changed.

For each chain, every C-alpha B-factor was divided by that chain's own median C-alpha B. The normalised values were averaged per helix. Chains were split into low, mid and high tertiles of their mean C-alpha B (results/helix_bfactor_by_chain.csv). Within each group, helices were ranked by median normalised B, with rank 1 the lowest median. The "better half" is rank 5 or better out of 11 helices. H12 is ranked but cannot enter the core.

## Most stable helices per B-factor group

The table gives the median and interquartile range of normalised B over all chains, and the rank in each tertile and over all chains (results/helix_stability_by_bgroup.csv).

| Helix | Range | Name | Median (all) | IQR (all) | Rank low / mid / high / all |
|---|---|---|---|---|---|
| AB3 | 306-310 | unknown | 1.618 | 0.442 | 11 / 10 / 10 / 10 |
| AB4 | 311-322 | unknown | 1.150 | 0.167 | 8 / 7 / 6 / 7 |
| AB5 | 338-362 | unknown | 0.981 | 0.118 | 3 / 4 / 4 / 4 |
| AB6 | 366-370 | unknown | 1.137 | 0.195 | 6 / 6 / 8 / 6 |
| AB7 | 371-396 | unknown | 0.824 | 0.088 | 2 / 2 / 2 / 2 |
| AB8 | 412-416 | unknown | 1.623 | 0.478 | 10 / 11 / 11 / 11 |
| AB9 | 420-439 | unknown | 1.051 | 0.124 | 5 / 5 / 5 / 5 |
| AC1 | 441-456 | unknown | 0.779 | 0.093 | 1 / 1 / 1 / 1 |
| AC2 | 465-493 | unknown | 1.153 | 0.169 | 7 / 8 / 7 / 8 |
| AC3 | 496-531 | unknown | 0.963 | 0.113 | 4 / 3 / 3 / 3 |
| AC4 | 537-546 | H12 (settings 538-547) | 1.279 | 0.348 | 9 / 9 / 9 / 9 |

AC1, AB7, AB5, AC3 and AB9 are in the better half in every tertile, so they form the recommended core. AB6 is rank 6 in the low and mid groups but rank 8 in the high group, so it is excluded. AB8 and AB3 are the least stable helices in every group. The figure figures/helix_bfactor_stability.png shows the median per helix and group.

## H12 / core B ratio

The ratio is H12 mean B divided by core mean B, with the core and H12 ranges from the settings (results/h12_core_ratio_by_bgroup.csv). Over all chains, the median ratio is 1.255 (quartiles 1.118 to 1.416, n = 765). The tertile medians are 1.299 (low), 1.224 (mid) and 1.252 (high). The normalised ratio is identical to the raw ratio in every row, because the per-chain median cancels in the division.

The earlier estimate, that H12 is about 12 times the core, does not hold for mean C-alpha B on this data. The measured median ratio is 1.255, so H12 runs about 1.25 times the core in mean B.

## Recommended core

The recommended core is five helices, 123 residues in total, never H12 (config/alignment_core_v3.json, core_ranges): AB5 338-362, AB7 371-396, AB9 420-439, AC1 441-456, AC3 496-531.

This does not equal the current core (340-360, 400-420, 450-470, 500-520, 84 residues). The differences are:

- 400-420 is replaced. The source annotates no helix between 397 and 411. AB8 (412-416) is the least stable helix in every group, so it is excluded. AB9 (420-439) is added.
- 450-470 is replaced by AC1 (441-456), the most stable helix in every group. Residues 457-470 are excluded, which covers the start of AC2 (465-493).
- 340-360 widens to AB5 (338-362).
- 500-520 widens to AC3 (496-531).
- AB7 (371-396) is added. It is rank 2 in every group.

The settings file analysis_settings.json still holds the current core. The S03 spec reads core_ranges from alignment_core_v3.json. The completeness screen still uses the current core, as before.

## Limits

B-factors depend on refinement, crystal packing and resolution. The normalisation removes the per-chain scale, but it does not remove differences in packing. The reference is an engineered, PROSS-designed variant with 28 declared substitutions, some inside the LBD. The ranking describes this reference set, not wild-type ERalpha.

## Files

results/helix_table.csv (helix ranges and names), results/helix_bfactor_by_chain.csv (per-chain B), results/helix_stability_by_bgroup.csv (ranks and IQR), results/h12_core_ratio_by_bgroup.csv (ratio), esr_analysis/resources/config/alignment_core_v3.json (recommended core), figures/helix_bfactor_stability.png (heatmap).
