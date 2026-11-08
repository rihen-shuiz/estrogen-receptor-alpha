Not done / substituted: (1) method (b) (render the SVG with the other axes hidden) not run: cairosvg 2.9.1 cannot load libcairo and rsvg-convert is absent (results/sv/VIS1/renderer_check.txt); macOS sips/qlmanage and Chrome can rasterise SVG but are not named by the step and would substitute the missing DejaVu Sans font, so they were not used. Every multi image without SV-script panels was cut by (c), PNG gutters; for the matplotlib figures with an SVG twin the SVG was read (not rendered) for the axes layout to choose the cut lines, so all panels are crops of the original PNG pixels (method column 'png_gutter'). (2) 9 multi images left uncut (not_cut, original used as is) because no rectangle separates two of their graphs without cutting one or taking a fragment of the other (method b would hide the other axes): Y09_physchem_by_cluster.png (axes_9 / axes_10: a title or tick label of one reaches 20 px past the edge of the other); Y16_outlier_quality.png (axes_3 / axes_4: a title or tick label of one reaches 17 px past the edge of the other); Y21_helix_b.png (the heatmap's colorbar label runs into the figure suptitle in the original); Y26_set_cover.png (axes_1 / axes_2: a title or tick label of one reaches 25 px past the edge of the other); SV01a_forest_completeness.png (axes_1 / axes_2: a title or tick label of one reaches 16 px past the edge of the other); SV01b_forest_completeness.png (axes_1 / axes_2: a title or tick label of one reaches 19 px past the edge of the other); SV01b_forest_mobility.png (axes_1 / axes_2: a title or tick label of one reaches 90 px past the edge of the other); SV06_outliers_map.png (axes_1 / axes_3: a title or tick label of one reaches 29 px past the edge of the other); SV14_filling.png (axes_1 / axes_2: a title or tick label of one reaches 37 px past the edge of the other). (3) Declared before cutting (scripts/sv/VIS1_spec.py): axes stacked on a shared axis whose tick labels are drawn on one axes only are one graph (Y10, Y20 loadings, Y23 residue diff, Y27 timeline; applied to SV02 helix track and SV15 dpx_by_group, whose SV-script panels lack those tick labels); a key needed to read a graph stays with it (Y09 scaffold drawings). (4) Check beyond the step: the SV-script panels of 7 SV figures fail the VIS1 check (neighbour fragments or a cut-off label at a panel edge, also found by an edge-ink scan); 2 of them were re-cut from the PNG (h12_map_partial.png, SV05_core_v4.png), the rest are among the not_cut images; the SV files themselves are untouched. (5) Shared legend files carry the same <parent>__ prefix as the panels (<parent>__<stem>__legend.png), because stems repeat across folders. (6) Supporting files not named in the run file: results/sv/VIS1/inventory.csv, cut_log.csv, panel_boxes.csv, renderer_check.txt; scripts/sv/VIS1_spec.py. Nothing installed.

# VIS1 - Cut multi-graph images into single graphs

**Verdict**

- 130 PNG paths, 128 distinct images by SHA-1 (2 duplicate paths); 0 smaller than 200 px.
- Counts over the 128 distinct images: single 58; cut into panels 61 images -> 245 panels (4 SV figures listed with 28 SV-script panels, method sv_script; 57 images cut by VIS1 into 217 panels, method png_gutter, plus 2 shared legend files); not_cut 9.
- Every contact sheet in figures/panels/_check/ was opened and looked at (79 sheets). Pass: 61 (57 VIS1 cuts, 4 SV figures with SV-script panels). Fail: 1 VIS1 cut (Y21_helix_b) and the SV-script panels of 9 SV figures (7 with fragments or cut-off labels, 2 lacking the shared-axis tick labels); failures were fixed where a PNG cut could fix them, else the image is not_cut (or single under the shared-axis rule) with the reason in the manifest.
- Originals unchanged: panels are crops only (no pixel edits, no redrawing); re-hashed after the run, 130 of 130 originals match their inventory SHA-1.

## Step 1 - inventory (results/sv/VIS1/inventory.csv)

- Scope: every PNG under the project except data/, figures/panels/, visuals/, virtual environments (folders with pyvenv.cfg) and hidden folders.
- By folder: figures/ext 54, figures/sv 35, figures/ 25, results/figures 11, results/archive_v2_strict 2, esr_analysis/resources 1, reference_render/ 1, results/reference_render 1.
- Identical by SHA-1 (cut once, all paths in the manifest): a1114ed454: esr_analysis/resources/reference/reference_context.png, figures/05_reference_context.png, results/figures/05_reference_context.png.

## Steps 2-4 - classification and method (results/sv/VIS1/cut_log.csv)

| class / method | images | paths | panels |
|---|---|---|---|
| single (one graph, incl. PyMOL renders; original used as is) | 58 | 60 | - |
| multi, (a) sv_script: SV-script panels listed | 4 | 4 | 28 |
| multi, (c) png_gutter guided by the SVG axes layout | 37 | 37 | 172 |
| multi, (c) png_gutter, PNG only | 20 | 20 | 45 (+2 legend) |
| multi, not_cut (original used as is) | 9 | 9 | - |

- Single SV figures: 24 (one panel saved by the SV script, or one graph under the shared-axis rule: SV02_helix_track, SV15_dpx_by_group); they are listed as single with their SV-script panel paths in the notes.
- (c) cuts: recursive XY-cut at the widest full-height/full-width background gutter. With an SVG twin, a cut may not cross any graph's extent (axes rectangle, legend frame and text boxes estimated from DejaVu Sans metrics, 4 px tolerance) and must have graphs on both sides; thin axes (colorbars, side bars) and overlapping axes (insets, twins) are joined to their graph. PNG-only: a cut must leave a graph-sized part on both sides (>= 12 % of the height, >= 8 % of the width, a straight ink run >= 8 % of the shorter side); declared cut directions and one declared cut column (07_distances) in scripts/sv/VIS1_spec.py. Crops: ink extent + 10 px padding, clipped at the cut lines. Strips without a graph (suptitles, footnotes) are dropped (112 strips in the cut images; panel_boxes.csv).

## Step 5 - shared legends

- 2 figure-level keys sit below all graphs and were saved once as <stem>__legend.png and referenced in needs_legend: figures_ext__Y30_pymol_bputty__legend.png, figures_ext__Y30_pymol_medoids__legend.png. Keys next to a graph stay in that graph's crop (e.g. SV13_forest C2b panel).

## Step 9 - check (figures/panels/_check/)

- One sheet per multi image (original left with numbered crop boxes, panels right); SV figures re-cut from the PNG also keep the sheet of their SV-script panels (<tag>__sv_panels.png); panels of not_cut images are written only to _check/_failed/ (21 files) so the failure can be seen.
- Not cut and why (manifest notes):
  - figures/ext/Y09_physchem_by_cluster.png: no background row or column separates axes_9+axes_10 (their extents from the SVG layout overlap by 20 px in x and 599 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 8 graphs found, 10 expected
  - figures/ext/Y16_outlier_quality.png: no background row or column separates axes_3+axes_4 (their extents from the SVG layout overlap by 17 px in x and 740 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 4 graphs found, 6 expected
  - figures/ext/Y21_helix_b.png: visual check fail: in the original the colorbar label of the heatmap ('... 1 = most stable') runs into the figure suptitle (overprinted on 'median'), so any rectangle that holds the whole colorbar label also holds a cut-off piece of the suptitle
  - figures/ext/Y26_set_cover.png: no background row or column separates axes_1+axes_2 (their extents from the SVG layout overlap by 25 px in x and 1047 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 0 graphs found, 2 expected
  - figures/sv/SV01a/SV01a_forest_completeness.png: no background row or column separates axes_1+axes_2 (their extents from the SVG layout overlap by 16 px in x and 1675 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 2 graphs found, 4 expected; SV-script panels (SV01a_forest_completeness__forest_cl4.png, SV01a_forest_completeness__forest_cl7.png, SV01a_forest_completeness__gains_cl4.png, SV01a_forest_completeness__gains_cl7.png) fail the VIS1 check: the forest_cl4 panel holds the end of the cl7 forest's title, x tick label and legend at its left edge, forest_cl7 a fragment of the cl4 tick labels at its right; the PNG cut fails too
  - figures/sv/SV01b/SV01b_forest_completeness.png: no background row or column separates axes_1+axes_2 (their extents from the SVG layout overlap by 19 px in x and 1743 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 2 graphs found, 4 expected; SV-script panels (SV01b_forest_completeness__forest_cl4.png, SV01b_forest_completeness__forest_cl7.png, SV01b_forest_completeness__gains_cl4.png, SV01b_forest_completeness__gains_cl7.png) fail the VIS1 check: the forest_cl4 panel holds fragments of the cl7 forest at its left edge, forest_cl7 fragments of the cl4 tick labels at its right; the PNG cut fails too
  - figures/sv/SV01b/SV01b_forest_mobility.png: no background row or column separates axes_1+axes_2 (their extents from the SVG layout overlap by 90 px in x and 1743 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 2 graphs found, 4 expected; SV-script panels (SV01b_forest_mobility__forest_cl4.png, SV01b_forest_mobility__forest_cl7.png, SV01b_forest_mobility__gains_cl4.png, SV01b_forest_mobility__gains_cl7.png) fail the VIS1 check: the forest_cl4 panel holds fragments of the cl7 forest at its left edge, forest_cl7 fragments of the cl4 tick labels at its right; the PNG cut fails too
  - figures/sv/SV06/SV06_outliers_map.png: no background row or column separates axes_1+axes_3 (their extents from the SVG layout overlap by 29 px in x and 1159 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 0 graphs found, 2 expected; SV-script panels (SV06_outliers_map__A_mds_map.png, SV06_outliers_map__B_own_vs_other_medoid.png) fail the VIS1 check: panel A holds the y-axis title of panel B at its right edge; the PNG cut fails too
  - figures/sv/SV14/SV14_filling.png: no background row or column separates axes_1+axes_2 (their extents from the SVG layout overlap by 37 px in x and 744 px in y: a tick label or title of one reaches past the edge of the other); a rectangular PNG crop would cut one graph or include a fragment of the other (method b, which hides the other axes, is unavailable); 0 graphs found, 2 expected; SV-script panels (SV14_filling__by_class.png, SV14_filling__by_cluster.png) fail the VIS1 check: by_cluster holds the y-axis title and tick labels of by_class at its right edge, by_class the end of the by_cluster title at its left; the PNG cut fails too
- SV-script panels that fail the check (SV files untouched):
  - figures/sv/SV01a/SV01a_forest_completeness.png: the forest_cl4 panel holds the end of the cl7 forest's title, x tick label and legend at its left edge, forest_cl7 a fragment of the cl4 tick labels at its right -> not_cut
  - figures/sv/SV01b/SV01b_forest_completeness.png: the forest_cl4 panel holds fragments of the cl7 forest at its left edge, forest_cl7 fragments of the cl4 tick labels at its right -> not_cut
  - figures/sv/SV01b/SV01b_forest_mobility.png: the forest_cl4 panel holds fragments of the cl7 forest at its left edge, forest_cl7 fragments of the cl4 tick labels at its right -> not_cut
  - figures/sv/SV04/h12_map_partial.png: the legend of panel b is cut off at its right edge ('(n=' without the count) -> re-cut from the PNG
  - figures/sv/SV05/SV05_core_v4.png: the y-axis title of panel c is cut off at its top edge -> re-cut from the PNG
  - figures/sv/SV06/SV06_outliers_map.png: panel A holds the y-axis title of panel B at its right edge -> not_cut
  - figures/sv/SV14/SV14_filling.png: by_cluster holds the y-axis title and tick labels of by_class at its right edge, by_class the end of the by_cluster title at its left -> not_cut
- Notes from the sheets: crops end exactly at a narrow gutter where graphs sit close (Y06 panel 01, Y21_core_selection panel 01; nothing cut); 06_alignment panel 01 starts at the left edge of the original, whose z-axis title is already clipped there; 11_rmsf_coverage and Y12 right-hand axes keep the shared x-axis title only under the lower panel (noted per panel).

## Files

- figures/panels/manifest.csv (84.0 KB)
- results/sv/VIS1/SUMMARY.md (13.5 KB)
- results/sv/VIS1/env.txt (0.5 KB)
- results/sv/VIS1/renderer_check.txt (1.0 KB)
- results/sv/VIS1/inventory.csv (12.5 KB)
- results/sv/VIS1/cut_log.csv (27.5 KB)
- results/sv/VIS1/panel_boxes.csv (40.8 KB)
- figures/panels/*.png: 219 files (217 panels + 2 legend files; 17.9 MB)
- figures/panels/_check/*.png: 79 contact sheets (30.5 MB); _check/_failed/: 21 files (1.6 MB)
- scripts/sv/VIS1_cut.py (37.6 KB)
- scripts/sv/VIS1_spec.py (17.7 KB)
- scripts/sv/VIS1_summary.py (15.1 KB)
- docs/EXPLORATION.md: appended 'Suggested, not done (VIS1, 2026-11-08)' (rule 1.1)
