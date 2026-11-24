"""Stage 12.5 - Render 10 diverse ERa LBD ligand complexes (focusing on H11, H12, and ligand) using headless PyMOL.

Usage:
  python pipeline/curated_action/05_pymol_10_ligands.py [--download-only | --render-only]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402

# 10 carefully selected complexes spanning all functional classes, drug classes, H12 states, and key outliers
SELECTED_10 = [
    dict(idx=1, pdb="7NEL", chain="A", lig="EST", title="17β-Estradiol (EST)",
         role="Natural full agonist (cl7, Y537S + pep)", lig_rgb="[0.11, 0.62, 0.47]", h12_rgb="[0.16, 0.47, 0.84]"),
    dict(idx=2, pdb="3ERD", chain="A", lig="DES", title="Diethylstilbestrol (DES)",
         role="Non-steroidal agonist (cl7, WT + pep)", lig_rgb="[0.11, 0.62, 0.47]", h12_rgb="[0.16, 0.47, 0.84]"),
    dict(idx=3, pdb="2QZO", chain="A", lig="KN1", title="KN1",
         role="Partial agonist in cl7 (Y537S + pep)", lig_rgb="[0.48, 0.20, 0.58]", h12_rgb="[0.16, 0.47, 0.84]"),
    dict(idx=4, pdb="1XPC", chain="A", lig="AIT", title="AIT (stereoisomer of AIJ)",
         role="Partial agonist in cl4 (WT, no pep)", lig_rgb="[0.48, 0.20, 0.58]", h12_rgb="[0.92, 0.41, 0.20]"),
    dict(idx=5, pdb="7UJC", chain="A", lig="RAL", title="Raloxifene (RAL)",
         role="Clinical SERM antagonist (cl4, Y537S)", lig_rgb="[0.76, 0.09, 0.36]", h12_rgb="[0.92, 0.41, 0.20]"),
    dict(idx=6, pdb="5UFX", chain="A", lig="86Y", title="86Y (oral SERD)",
         role="SERD antagonist (cl4, L536S)", lig_rgb="[0.76, 0.09, 0.36]", h12_rgb="[0.92, 0.41, 0.20]"),
    dict(idx=7, pdb="6OWC", chain="A", lig="ND1", title="H3B-6545 (ND1)",
         role="Covalent SERCA antagonist (Y537S)", lig_rgb="[0.76, 0.09, 0.36]", h12_rgb="[0.92, 0.41, 0.20]"),
    dict(idx=8, pdb="5KD9", chain="B", lig="OBT", title="OBT (cl7 Tukey outlier +1.65 Å)",
         role="Antagonist trapped in cl7 (Y537S + pep)", lig_rgb="[0.76, 0.09, 0.36]", h12_rgb="[0.16, 0.47, 0.84]"),
    dict(idx=9, pdb="6SBO", chain="A", lig="L5B", title="Amcenestrant (L5B)",
         role="Clinical SERD, action=unknown (cl4, L536S)", lig_rgb="[0.45, 0.45, 0.45]", h12_rgb="[0.92, 0.41, 0.20]"),
    dict(idx=10, pdb="9B25", chain="A", lig="A1AIZ", title="A1AIZ (minor intermediate H12)",
         role="Intermediate H12 outlier (minor, L536S)", lig_rgb="[0.45, 0.45, 0.45]", h12_rgb="[0.83, 0.63, 0.09]"),
]

STRUCT_DIR = cc.PYMOL_DIR / "structures"
PML_DIR = cc.PYMOL_DIR / "pml"
SHIM_DIR = Path(__file__).resolve().parent / "pymol_shim"
REF_PDB = "7NEL"  # All 10 structures aligned on core helices onto 7NEL A for identical H11/H12 viewing frame


def download_structures():
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    pdbs = sorted({x["pdb"] for x in SELECTED_10})
    for pdb in pdbs:
        out_cif = STRUCT_DIR / f"{pdb}.cif"
        if out_cif.exists() and out_cif.stat().st_size > 10_000:
            print(f"[cache] {out_cif.name}")
            continue
        url = f"https://files.rcsb.org/download/{pdb}.cif"
        print(f"[download] {url} -> {out_cif.name}")
        urllib.request.urlretrieve(url, out_cif)


def add_caption_banner(raw_png: Path, out_png: Path, item: dict, meta_row: pd.Series | None):
    img = Image.open(raw_png).convert("RGB")
    w, h = img.size
    banner_h = 115
    canvas = Image.new("RGB", (w, h + banner_h), (255, 255, 255))
    canvas.paste(img, (0, banner_h))
    draw = ImageDraw.Draw(canvas)
    try:
        f_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 30)
        f_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
        f_leg = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    except Exception:
        f_title = f_sub = f_leg = ImageFont.load_default()

    header = f"#{item['idx']:02d}  {item['pdb']} (chain {item['chain']}) — {item['title']}"
    draw.text((28, 12), header, fill=(11, 11, 11), font=f_title)
    if meta_row is not None:
        dist_v = meta_row.get("dist_H12c_H11axis_A", float("nan"))
        dist_s = f"{float(dist_v):.2f} Å" if pd.notna(dist_v) else "n/a"
        res_s = f"{float(meta_row['resolution_A']):.2f} Å" if pd.notna(meta_row.get("resolution_A")) else "n/a"
        dc_s = str(meta_row.get("drug_classes", "")) if pd.notna(meta_row.get("drug_classes")) else "—"
        sub = (f"Class: {meta_row['action_class']} | Drug: {dc_s} | H12 group: {meta_row['group']} | "
               f"d(H12–H11)={dist_s} | Res: {res_s} | {item['role']}")
    else:
        sub = item["role"]
    draw.text((28, 52), sub, fill=(75, 75, 72), font=f_sub)
    draw.text((28, 84), "Colors:  H11 (497–526) steel-blue  |  H12 (538–547) blue=cl7 / orange=cl4 / gold=minor  |  Ligand = sticks  |  LBD = grey cartoon",
              fill=(110, 110, 105), font=f_leg)
    draw.line([(20, banner_h - 2), (w - 20, banner_h - 2)], fill=(220, 220, 215), width=2)
    canvas.save(out_png)


def render_one(item: dict, meta_row: pd.Series | None):
    PML_DIR.mkdir(parents=True, exist_ok=True)
    ref_cif = (STRUCT_DIR / f"{REF_PDB}.cif").resolve()
    mob_cif = (STRUCT_DIR / f"{item['pdb']}.cif").resolve()
    raw_png = (PML_DIR / f"raw_{item['idx']:02d}_{item['pdb']}_{item['lig']}.png").resolve()
    final_png = (cc.PYMOL_DIR / f"{item['idx']:02d}_{item['pdb']}_{item['chain']}_{item['lig']}.png").resolve()
    pml_path = PML_DIR / f"{item['idx']:02d}_{item['pdb']}_{item['lig']}.pml"

    ch = item["chain"]
    lig = item["lig"]
    # Align target chain onto reference 7NEL chain A core (342-361 + 372-394 + 422-437 + 442-455 + 497-526)
    pml = f"""
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set cartoon_fancy_helices, 1
set cartoon_highlight_color, grey85
set cartoon_transparency, 0.55
set stick_radius, 0.22
set depth_cue, 0
set orthoscopic, 1

load {ref_cif}, ref_mol
load {mob_cif}, mob_mol

# Keep only target receptor chain + its ligand (and peptide if short chain present)
remove mob_mol and solvent
remove mob_mol and inorganic

select ref_core, ref_mol and chain A and name CA and (resi 342-361 or resi 372-394 or resi 422-437 or resi 442-455 or resi 497-526)
select mob_core, mob_mol and chain {ch} and name CA and (resi 342-361 or resi 372-394 or resi 422-437 or resi 442-455 or resi 497-526)
align mob_core, ref_core

hide everything
delete ref_mol

select rec_chain, mob_mol and polymer.protein and chain {ch}
select h11_sel, rec_chain and resi 497-526
select loop_sel, rec_chain and resi 527-537
select h12_sel, rec_chain and resi 538-547
select lig_sel, mob_mol and resn {lig} and (chain {ch} or (byres (mob_mol and resn {lig}) within 6 of rec_chain))
select pep_sel, mob_mol and polymer.protein and not chain {ch} and (byres (mob_mol and polymer.protein and not chain {ch}) within 6 of rec_chain)

show cartoon, rec_chain
color grey80, rec_chain

set_color h11_col, [0.23, 0.41, 0.60]
set_color h12_col, {item['h12_rgb']}
set_color lig_col, {item['lig_rgb']}
set_color loop_col, [0.55, 0.55, 0.55]

color h11_col, h11_sel
color loop_col, loop_sel
color h12_col, h12_sel
set cartoon_transparency, 0.0, h11_sel or h12_sel or loop_sel

show cartoon, pep_sel
color palegreen, pep_sel
set cartoon_transparency, 0.25, pep_sel

show sticks, lig_sel
color lig_col, lig_sel and elem C
util.cnc lig_sel
show spheres, lig_sel
set sphere_scale, 0.20, lig_sel

# Focus view tightly on H11 + H12 + ligand in consistent reference orientation
orient (h11_sel or h12_sel or lig_sel)
zoom (h11_sel or h12_sel or lig_sel), 6.0

ray 1400, 1000
png {raw_png}
quit
"""
    pml_path.write_text(pml.strip() + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["DISPLAY"] = ""
    env["PYTHONPATH"] = str(SHIM_DIR) + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    cmd = ["/usr/bin/python3", "/usr/lib/python3/dist-packages/pymol/__init__.py", "-ckq", str(pml_path)]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
    if res.returncode != 0 or not raw_png.exists():
        raise RuntimeError(f"PyMOL failed for {item['pdb']}: rc={res.returncode}\n{res.stderr}\n{res.stdout}")
    add_caption_banner(raw_png, final_png, item, meta_row)
    print(f"[pymol] rendered {final_png.name} ({final_png.stat().st_size // 1024} KB)")


def make_contact_sheet():
    """Assemble a 2x5 contact sheet of all 10 PyMOL ligand visualizations."""
    imgs = []
    for item in SELECTED_10:
        p = cc.PYMOL_DIR / f"{item['idx']:02d}_{item['pdb']}_{item['chain']}_{item['lig']}.png"
        imgs.append(Image.open(p).convert("RGB").resize((700, 558), Image.Resampling.LANCZOS))
    sheet = Image.new("RGB", (700 * 5, 558 * 2), (255, 255, 255))
    for i, im in enumerate(imgs):
        sheet.paste(im, ((i % 5) * 700, (i // 5) * 558))
    out = cc.PYMOL_DIR / "00_contact_sheet_10_ligands.png"
    sheet.save(out)
    print(f"[pymol] contact sheet -> {out.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download-only", action="store_true")
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    cc.ensure_dirs()

    if not args.render_only:
        download_structures()
    if args.download_only:
        return

    ch = cc.read_csv_safe(cc.TABLES_DIR / "chains_curated.csv")
    sel_rows = []
    for item in SELECTED_10:
        sub = ch[(ch["pdb_norm"] == item["pdb"]) & (ch["auth_asym_id"] == item["chain"])]
        meta = sub.iloc[0] if len(sub) else None
        render_one(item, meta)
        sel_rows.append({
            "idx": item["idx"],
            "pdb_id": item["pdb"],
            "chain": item["chain"],
            "ligand_id": item["lig"],
            "title": item["title"],
            "selection_role": item["role"],
            "action_class": meta["action_class"] if meta is not None else "",
            "drug_classes": meta["drug_classes"] if meta is not None else "",
            "h12_group": meta["group"] if meta is not None else "",
            "mutation_status": meta["mutation_status"] if meta is not None else "",
            "peptide_bound": meta["peptide_bound"] if meta is not None else "",
            "dist_H12c_H11axis_A": meta["dist_H12c_H11axis_A"] if meta is not None else "",
            "resolution_A": meta["resolution_A"] if meta is not None else "",
            "png_file": f"results/12_curated_action/pymol_ligands/{item['idx']:02d}_{item['pdb']}_{item['chain']}_{item['lig']}.png",
        })
    pd.DataFrame(sel_rows).to_csv(cc.PYMOL_DIR / "ligands_10_selection.csv", index=False)
    make_contact_sheet()
    print("=== Stage 12.5 PyMOL 10-ligand gallery complete ===")


if __name__ == "__main__":
    main()
