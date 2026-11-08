"""VIS1 SUMMARY.md, env.txt and renderer_check.txt, written from the saved outputs only (rule 1.7).

Inputs: results/sv/VIS1/{inventory.csv, cut_log.csv, panel_boxes.csv}, figures/panels/manifest.csv,
figures/panels/*.png, figures/panels/_check/*.png, scripts/sv/VIS1_spec.py (classification, check verdicts).
Run after scripts/sv/VIS1_cut.py: .venv/bin/python -B scripts/sv/VIS1_summary.py
"""
from __future__ import annotations

import platform
import re
import shutil
import sys
import time
from importlib.metadata import version
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/sv/VIS1"
PANELS = ROOT / "figures/panels"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import VIS1_spec as spec  # noqa: E402


def size(p: Path) -> str:
    n = p.stat().st_size
    return f"{n / 1024 / 1024:.1f} MB" if n >= 1024 * 1024 else f"{n / 1024:.1f} KB"


def dir_size(files) -> str:
    n = sum(f.stat().st_size for f in files)
    return f"{n / 1024 / 1024:.1f} MB"


def renderer_check() -> list[str]:
    """Which SVG renderers exist on this machine (nothing installed)."""
    out = [f"# VIS1 SVG renderer check, {time.strftime('%Y-%m-%dT%H:%M%z')}"]
    for tool in ("rsvg-convert", "inkscape", "magick", "convert", "resvg", "cairosvg"):
        out.append(f"{tool}: {shutil.which(tool) or 'not found'}")
    try:
        import cairosvg  # noqa: F401
        out.append("python cairosvg: imports")
    except OSError as e:
        out.append("python cairosvg " + version("cairosvg") + ": import fails, " + str(e).splitlines()[0]
                   + " (system libcairo missing)")
    for tool in ("sips", "qlmanage"):
        out.append(f"{tool}: {shutil.which(tool) or 'not found'} (macOS system tool)")
    chrome = Path("/Applications/Google Chrome.app")
    out.append(f"Google Chrome.app: {'present' if chrome.exists() else 'absent'}")
    svgs = sorted((ROOT / "figures").rglob("*.svg"))
    live = [s for s in svgs if "<text" in s.read_text(errors="ignore")]
    out.append(f"SVG twins under figures/: {len(svgs)}; with live <text> elements (matplotlib svg.fonttype 'none', "
               f"scripts/ext_common.py): {len(live)}")
    font_dirs = [Path("/System/Library/Fonts"), Path("/Library/Fonts"), Path.home() / "Library/Fonts"]
    dv = [f for d in font_dirs if d.exists() for f in d.rglob("*") if "dejavu" in f.name.lower()]
    out.append(f"DejaVu Sans installed system-wide (font used in the SVG text): {'yes' if dv else 'no'}")
    out.append("Decision: method (b) needs cairosvg or rsvg-convert (step 4b); neither works here, so (b) is "
               "unavailable and every multi image without SV-script panels uses (c). sips/qlmanage/Chrome can "
               "rasterise SVG but are not named by the step and would substitute the missing DejaVu Sans font, so the "
               "render would not reproduce the original; not used (rule 1.5). The SVG is only read for the axes "
               "layout (top-level <g id='axes_N'> rectangles and text anchors) to place cuts in the PNG.")
    return out


def main():
    inv = pd.read_csv(RES / "inventory.csv")
    log = pd.read_csv(RES / "cut_log.csv", keep_default_na=False)
    man = pd.read_csv(PANELS / "manifest.csv", keep_default_na=False)
    boxes = pd.read_csv(RES / "panel_boxes.csv", keep_default_na=False)
    panels = sorted(PANELS.glob("*.png"))
    sheets = sorted((PANELS / "_check").glob("*.png"))
    failed_dir = sorted((PANELS / "_check/_failed").glob("*.png"))

    n_paths, n_uniq = len(inv), inv.sha1.nunique()
    dup = inv[inv.duplicated("sha1", keep=False)].sort_values("path")
    n_small = int(inv.small.astype(str).eq("True").sum())
    single = log[log.status == "single"]
    cut = log[log.status == "cut"]
    svl = log[log.status == "sv_listed"]
    nc = log[log.status == "not_cut"]
    sv_mask = log.path.str.startswith("figures/sv/")
    n_cut_panels = int(cut.n_panels.sum())
    n_sv_panels = int(svl.n_panels.sum())
    leg = man[man.panel_title == "shared legend / colour key"]
    gut = cut[cut["mode"] == "gutter"]
    svgm = cut[cut["mode"] == "svg"]
    recut = cut[sv_mask[cut.index]]
    sv_single = single[sv_mask[single.index]]
    chk_pass = sum(1 for v in spec.CHECK.values() if v == "pass")
    chk_fail = sum(1 for v in spec.CHECK.values() if v.startswith("fail"))
    paths_of = lambda df: df.all_paths.str.split(";").str.len().sum()  # noqa: E731

    # renderer check + env first (cited in SUMMARY)
    rc = renderer_check()
    (RES / "renderer_check.txt").write_text("\n".join(rc) + "\n")
    env = [f"# VIS1 environment, {time.strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for p in ("numpy", "pandas", "Pillow", "matplotlib"):
        env.append(f"{p} {version(p)}")
    env += ["matplotlib is used only for its bundled DejaVu Sans font files (text-extent estimates, contact-sheet "
            "labels); no figure is redrawn",
            "SVG renderer: none usable (see results/sv/VIS1/renderer_check.txt); method (b) not run",
            "network: none", "nothing installed, upgraded or downgraded by VIS1"]
    (RES / "env.txt").write_text("\n".join(env) + "\n")

    nc_list = []
    for r in nc.itertuples():
        m = re.search(r"separates (\S+) \(their extents from the SVG layout overlap by (\d+) px in x", r.note)
        if m:
            why = f"{m.group(1).replace('+', ' / ')}: a title or tick label of one reaches {m.group(2)} px past the edge of the other"
        elif "suptitle" in r.note:
            why = "the heatmap's colorbar label runs into the figure suptitle in the original"
        else:
            why = r.note.split(";")[0]
        nc_list.append(f"{Path(r.path).name} ({why})")

    line1 = (
        "Not done / substituted: (1) method (b) (render the SVG with the other axes hidden) not run: cairosvg "
        f"{version('cairosvg')} cannot load libcairo and rsvg-convert is absent (results/sv/VIS1/renderer_check.txt); "
        "macOS sips/qlmanage and Chrome can rasterise SVG but are not named by the step and would substitute the "
        "missing DejaVu Sans font, so they were not used. Every multi image without SV-script panels was cut by (c), "
        "PNG gutters; for the matplotlib figures with an SVG twin the SVG was read (not rendered) for the axes layout "
        "to choose the cut lines, so all panels are crops of the original PNG pixels (method column 'png_gutter'). "
        f"(2) {len(nc)} multi images left uncut (not_cut, original used as is) because no rectangle separates two of "
        "their graphs without cutting one or taking a fragment of the other (method b would hide the other axes): "
        + "; ".join(nc_list) + ". "
        f"(3) Declared before cutting (scripts/sv/VIS1_spec.py): axes stacked on a shared axis whose tick labels are "
        "drawn on one axes only are one graph (Y10, Y20 loadings, Y23 residue diff, Y27 timeline; applied to SV02 "
        "helix track and SV15 dpx_by_group, whose SV-script panels lack those tick labels); a key needed to read a "
        "graph stays with it (Y09 scaffold drawings). (4) Check beyond the step: the SV-script panels of "
        f"{len(spec.SV_RECUT)} SV figures fail the VIS1 check (neighbour fragments or a cut-off label at a panel "
        f"edge, also found by an edge-ink scan); {len(recut)} of them were re-cut from the PNG "
        f"({', '.join(Path(p).name for p in recut.path)}), the rest are among the not_cut images; the SV files "
        "themselves are untouched. (5) Shared legend files carry the same <parent>__ prefix as the panels "
        "(<parent>__<stem>__legend.png), because stems repeat across folders. (6) Supporting files not named in the "
        "run file: results/sv/VIS1/inventory.csv, cut_log.csv, panel_boxes.csv, renderer_check.txt; "
        "scripts/sv/VIS1_spec.py. Nothing installed.")

    L = [line1, "", "# VIS1 - Cut multi-graph images into single graphs", "", "**Verdict**", ""]
    L.append(f"- {n_paths} PNG paths, {n_uniq} distinct images by SHA-1 ({n_paths - n_uniq} duplicate paths); "
             f"{n_small} smaller than 200 px.")
    L.append(f"- Counts over the {n_uniq} distinct images: single {len(single)}; cut into panels "
             f"{len(cut) + len(svl)} images -> {n_cut_panels + n_sv_panels} panels "
             f"({len(svl)} SV figures listed with {n_sv_panels} SV-script panels, method sv_script; "
             f"{len(cut)} images cut by VIS1 into {n_cut_panels} panels, method png_gutter, plus {len(leg)} shared "
             f"legend files); not_cut {len(nc)}.")
    L.append(f"- Every contact sheet in figures/panels/_check/ was opened and looked at ({len(sheets)} sheets). "
             f"Pass: {chk_pass} ({len(cut)} VIS1 cuts, {len(svl)} SV figures with SV-script panels). Fail: "
             f"{chk_fail} VIS1 cut (Y21_helix_b) and the SV-script panels of {len(spec.SV_RECUT) + len(spec.SV_SINGLE)} "
             f"SV figures ({len(spec.SV_RECUT)} with fragments or cut-off labels, {len(spec.SV_SINGLE)} lacking the "
             "shared-axis tick labels); failures were fixed where a PNG cut could fix them, else the image is not_cut "
             "(or single under the shared-axis rule) with the reason in the manifest.")
    import hashlib
    same = sum(hashlib.sha1((ROOT / r.path).read_bytes()).hexdigest() == r.sha1 for r in inv.itertuples())
    L.append(f"- Originals unchanged: panels are crops only (no pixel edits, no redrawing); re-hashed after the run, "
             f"{same} of {n_paths} originals match their inventory SHA-1.")
    L.append("")
    L.append("## Step 1 - inventory (results/sv/VIS1/inventory.csv)")
    L.append("")
    L.append("- Scope: every PNG under the project except data/, figures/panels/, visuals/, virtual environments "
             "(folders with pyvenv.cfg) and hidden folders.")
    by_dir = inv.assign(top=inv.path.map(lambda p: "/".join(p.split("/")[:2]) if p.count("/") >= 2 else
                                         p.split("/")[0] + "/")).top.value_counts()
    L.append("- By folder: " + ", ".join(f"{k} {v}" for k, v in by_dir.items()) + ".")
    L.append("- Identical by SHA-1 (cut once, all paths in the manifest): "
             + "; ".join(f"{g.sha1.iloc[0][:10]}: " + ", ".join(g.path) for _, g in dup.groupby("sha1")) + ".")
    L.append("")
    L.append("## Steps 2-4 - classification and method (results/sv/VIS1/cut_log.csv)")
    L.append("")
    L.append("| class / method | images | paths | panels |")
    L.append("|---|---|---|---|")
    L.append(f"| single (one graph, incl. PyMOL renders; original used as is) | {len(single)} | "
             f"{paths_of(single)} | - |")
    L.append(f"| multi, (a) sv_script: SV-script panels listed | {len(svl)} | {paths_of(svl)} | {n_sv_panels} |")
    L.append(f"| multi, (c) png_gutter guided by the SVG axes layout | {len(svgm)} | {paths_of(svgm)} | "
             f"{int(svgm.n_panels.sum())} |")
    L.append(f"| multi, (c) png_gutter, PNG only | {len(gut)} | {paths_of(gut)} | {int(gut.n_panels.sum())} "
             f"(+{len(leg)} legend) |")
    L.append(f"| multi, not_cut (original used as is) | {len(nc)} | {paths_of(nc)} | - |")
    L.append("")
    L.append(f"- Single SV figures: {len(sv_single)} (one panel saved by the SV script, or one graph under the "
             "shared-axis rule: SV02_helix_track, SV15_dpx_by_group); they are listed as single with their "
             "SV-script panel paths in the notes.")
    L.append("- (c) cuts: recursive XY-cut at the widest full-height/full-width background gutter. With an SVG twin, "
             "a cut may not cross any graph's extent (axes rectangle, legend frame and text boxes estimated from "
             "DejaVu Sans metrics, 4 px tolerance) and must have graphs on both sides; thin axes (colorbars, side "
             "bars) and overlapping axes (insets, twins) are joined to their graph. PNG-only: a cut must leave a "
             "graph-sized part on both sides (>= 12 % of the height, >= 8 % of the width, a straight ink run >= "
             "8 % of the shorter side); declared cut directions and one declared cut column (07_distances) in "
             "scripts/sv/VIS1_spec.py. Crops: ink extent + 10 px padding, clipped at the cut lines. Strips without "
             f"a graph (suptitles, footnotes) are dropped ({int(((boxes.status == 'dropped') & boxes.source.isin(cut.path)).sum())} "
             "strips in the cut images; panel_boxes.csv).")
    L.append("")
    L.append("## Step 5 - shared legends")
    L.append("")
    L.append(f"- {len(leg)} figure-level keys sit below all graphs and were saved once as <stem>__legend.png and "
             "referenced in needs_legend: " + ", ".join(Path(p).name for p in leg.panel_path) + ". Keys next to a "
             "graph stay in that graph's crop (e.g. SV13_forest C2b panel).")
    L.append("")
    L.append("## Step 9 - check (figures/panels/_check/)")
    L.append("")
    L.append("- One sheet per multi image (original left with numbered crop boxes, panels right); SV figures re-cut "
             "from the PNG also keep the sheet of their SV-script panels (<tag>__sv_panels.png); panels of not_cut "
             f"images are written only to _check/_failed/ ({len(failed_dir)} files) so the failure can be seen.")
    L.append("- Not cut and why (manifest notes):")
    for r in nc.itertuples():
        L.append(f"  - {r.path}: {r.note}")
    L.append("- SV-script panels that fail the check (SV files untouched):")
    for k, v in spec.SV_RECUT.items():
        st = log.set_index("path").loc[k, "status"]
        L.append(f"  - {k}: {v['reason']} -> {'re-cut from the PNG' if st == 'cut' else 'not_cut'}")
    L.append("- Notes from the sheets: crops end exactly at a narrow gutter where graphs sit close (Y06 panel 01, "
             "Y21_core_selection panel 01; nothing cut); 06_alignment panel 01 starts at the left edge of the "
             "original, whose z-axis title is already clipped there; 11_rmsf_coverage and Y12 right-hand axes keep "
             "the shared x-axis title only under the lower panel (noted per panel).")
    L.append("")
    L.append("## Files")
    L.append("")
    files = [PANELS / "manifest.csv", RES / "SUMMARY.md", RES / "env.txt", RES / "renderer_check.txt",
             RES / "inventory.csv", RES / "cut_log.csv", RES / "panel_boxes.csv"]
    scripts = [ROOT / "scripts/sv" / x for x in ("VIS1_cut.py", "VIS1_spec.py", "VIS1_summary.py")]
    head = list(L)
    (RES / "SUMMARY.md").write_text("\n".join(head) + "\n")
    for _ in range(3):  # converge the listed size of SUMMARY.md itself
        L = head + [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files]
        L.append(f"- figures/panels/*.png: {len(panels)} files ({n_cut_panels} panels + {len(leg)} legend files; "
                 f"{dir_size(panels)})")
        L.append(f"- figures/panels/_check/*.png: {len(sheets)} contact sheets ({dir_size(sheets)}); "
                 f"_check/_failed/: {len(failed_dir)} files ({dir_size(failed_dir)})")
        L += [f"- {p.relative_to(ROOT)} ({size(p)})" for p in scripts]
        L.append("- docs/EXPLORATION.md: appended 'Suggested, not done (VIS1, 2026-11-08)' (rule 1.1)")
        (RES / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
