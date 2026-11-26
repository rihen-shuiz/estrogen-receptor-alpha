"""Stage 12.5 - PyMOL views of an agonist and an antagonist complex: where H12 goes and what the antagonist blocks.

Structures (Shiau et al. 1998, Cell 95:927, both wild-type LBD): 3ERD chain A, diethylstilbestrol (agonist, H12 closed
over the pocket, GRIP1 coactivator peptide chain C); 3ERT chain A, 4-hydroxytamoxifen (antagonist, H12 in the
coactivator groove). 3ERT is superposed on 3ERD with the C-alpha atoms of the alignment core (core_v3), with outlier
rejection.
mmCIF files are read from $ESR1_MMCIF_DIR or results/12_curated_action/pymol/structures/, and downloaded from RCSB when
missing.
Outputs: results/12_curated_action/pymol/*.png (ray-traced views), results/12_curated_action/tables/pymol_contacts.csv
(closest approach of the antagonist to the agonist H12), figures/pymol_agonist_antagonist.png,
figures/pymol_h12_blocked.png and figures/pymol_esr2_displacement.png (3ERD coloured by the ESR1-ESR2 C-alpha
distance of data/esr2/esr1_esr2_displacement_profile.csv; labelled panels).
"""
from __future__ import annotations

import json
import os
import sys
import textwrap
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from pymol import cmd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402

OUT = cc.OUT_ROOT / "pymol"
CIF_DIR = Path(os.environ.get("ESR1_MMCIF_DIR", OUT / "structures"))
FIGS = cc.ROOT / "figures"
AGO, ANT = ("3ERD", "A", "DES", "C"), ("3ERT", "A", "OHT")
H11, H12 = "497-530", "536-549"
BLUE, ORANGE, GREEN, CRIMSON, PEPTIDE = "0x2a78d6", "0xeb6834", "0x1b9e77", "0xc2185b", "0x8fbf7f"


def cif(pdb: str) -> Path:
    path = CIF_DIR / f"{pdb}.cif"
    if not path.exists():
        CIF_DIR.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.cif", path)
    return path


def load_pair():
    cmd.reinitialize()
    core = json.loads((cc.DATA / "coordinates/core_v3.json").read_text())["core_ranges"]
    core_sel = " or ".join(f"resi {a}-{b}" for a, b in core)
    cmd.load(str(cif(AGO[0])), "ago")
    cmd.load(str(cif(ANT[0])), "ant")
    cmd.remove("solvent or inorganic or resn ACY")
    cmd.remove(f"ago and not (chain {AGO[1]} or chain {AGO[3]})")
    cmd.remove(f"ant and not chain {ANT[1]}")
    cmd.remove("not alt ''+A")  # keep the first alternate conformation only
    cmd.alter("all", "alt=''")
    # core C-alpha fit with outlier rejection (5 cycles): the core's last residues (528-531) follow the lid
    r = cmd.align(f"ant and chain {ANT[1]} and name CA and ({core_sel})", f"ago and chain {AGO[1]} and name CA and ({core_sel})",
                  cycles=5)
    return r[0], r[1]


def style_scene():
    cmd.bg_color("white")
    for k, v in (("ray_opaque_background", 1), ("antialias", 2), ("ray_trace_mode", 1), ("ray_trace_gain", 0.4),
                 ("ray_shadows", 0), ("depth_cue", 0), ("orthoscopic", 1), ("cartoon_fancy_helices", 1),
                 ("cartoon_transparency", 0.0), ("stick_radius", 0.28), ("sphere_transparency", 0.35)):
        cmd.set(k, v)
    cmd.hide("everything")


def set_camera(sel: str, turns=(("x", -90), ("y", 20)), zoom: float | None = None):
    """Absolute camera: orient from a fixed starting view (orient alone depends on the current view)."""
    cmd.set_view((1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -50, 0, 0, 0, 40, 100, -20))
    cmd.orient(sel)
    for axis, deg in turns:
        cmd.turn(axis, deg)
    if zoom is not None:
        cmd.zoom(sel, zoom)
    return cmd.get_view()


def project(sel: str, width: int, height: int) -> tuple[float, float]:
    """Pixel position of a selection's centre in an orthoscopic ray-traced image (for labels drawn afterwards)."""
    v = np.array(cmd.get_view())
    rot = v[:9].reshape(3, 3)
    cam = (np.mean(cmd.get_coords(sel), axis=0) - v[12:15]) @ rot + v[9:12]
    scale = height / (2 * abs(v[11]) * np.tan(np.radians(cmd.get_setting_float("field_of_view")) / 2))
    return width / 2 + cam[0] * scale, height / 2 - cam[1] * scale


def side_camera(axis_sel: str, base_sel: str, tip_sel: str, distance: float = 60.0):
    """Camera with the helix axis left-right and the ligand's base-to-tip direction pointing up the screen."""
    ca = cmd.get_coords(axis_sel)
    c = ca.mean(axis=0)
    x = np.linalg.svd(ca - c)[2][0]
    base = cmd.get_coords(base_sel).mean(axis=0)
    t = cmd.get_coords(tip_sel).mean(axis=0) - base
    y = t - (t @ x) * x
    y /= np.linalg.norm(y)
    rot = np.column_stack([x, y, np.cross(x, y)])
    origin = 0.4 * c + 0.6 * base
    cmd.set_view(list(rot.reshape(-1)) + [0.0, 0.0, -distance] + list(origin) + [distance - 40, distance + 40, -20.0])


def render(png: Path, labels: dict[str, tuple[str, tuple[float, float]]], w=1600, h=1300) -> dict:
    """Ray-trace, then return each label's anchor pixel and its text direction (dx, dy; screen y points down)."""
    cmd.png(str(png), width=w, height=h, dpi=300, ray=1)
    return {text: (*project(sel, w, h), d) for text, (sel, d) in labels.items()}


def complex_view(obj, chain, lig, h12_color, lig_color, peptide_chain, png: Path, camera) -> dict:
    style_scene()
    rec = f"{obj} and chain {chain} and polymer"
    cmd.show("cartoon", rec)
    cmd.color("grey80", rec)
    cmd.set("cartoon_transparency", 0.35, rec)
    cmd.color("0x5b7ea6", f"{rec} and resi {H11}")
    cmd.color(h12_color, f"{rec} and resi {H12}")
    cmd.set("cartoon_transparency", 0.0, f"{rec} and (resi {H11} or resi {H12})")
    ligand = f"{obj} and resn {lig} and chain {chain}"
    cmd.show("sticks", ligand)
    cmd.color(lig_color, f"{ligand} and elem C")
    cmd.util.cnc(ligand)
    labels = {"H12": (f"{rec} and resi 538-547 and name CA", (1.0, 0.25) if peptide_chain else (0.6, 0.9)),
              "H11": (f"{rec} and resi 500-515 and name CA", (0.3, -1.0)), "ligand": (ligand, (-1.0, -0.2))}
    if peptide_chain:
        cmd.show("cartoon", f"{obj} and chain {peptide_chain}")
        cmd.color(PEPTIDE, f"{obj} and chain {peptide_chain}")
        labels["coactivator peptide"] = (f"{obj} and chain {peptide_chain} and name CA", (0.5, 1.0))
    cmd.set_view(camera)
    return render(png, labels)


def blocked_view(png: Path) -> tuple[dict, dict]:
    """Closed (agonist) H12 with the antagonist ligand from the superposed 3ERT: the side chain lands on H12."""
    style_scene()
    rec = f"ago and chain {AGO[1]} and polymer"
    h12 = f"{rec} and resi {H12}"
    oht = f"ant and resn {ANT[2]}"
    des = f"ago and resn {AGO[2]} and chain {AGO[1]}"
    cmd.show("cartoon", rec)
    cmd.color("grey85", rec)
    cmd.set("cartoon_transparency", 0.7, rec)
    cmd.color(BLUE, h12)
    cmd.set("cartoon_transparency", 0.0, h12)
    cmd.show("sticks", f"{oht} or {des}")
    cmd.color(CRIMSON, f"{oht} and elem C")
    cmd.color(GREEN, f"{des} and elem C")
    cmd.util.cnc(f"{oht} or {des}")
    # closest approach of the antagonist to the closed (agonist) H12
    m_lig, m_h12 = cmd.get_model(oht), cmd.get_model(f"{h12} and not hydro")
    a = np.array([x.coord for x in m_lig.atom])
    b = np.array([x.coord for x in m_h12.atom])
    d = np.linalg.norm(a[:, None] - b[None], axis=2)
    i, j = np.unravel_index(d.argmin(), d.shape)
    la, ha = m_lig.atom[i], m_h12.atom[j]
    near = "+".join(x.name for x, dmin in zip(m_lig.atom, d.min(axis=1)) if dmin < 3.0)
    cmd.show("spheres", f"{oht} and name {near}")
    cmd.show("sticks", f"{h12} and resi {ha.resi} and not name N+C+O")
    cmd.show("spheres", f"{h12} and resi {ha.resi} and not name N+C+O")
    cmd.color(BLUE, f"{h12} and resi {ha.resi} and elem C")
    cent = lambda sel: np.mean(cmd.get_coords(sel), axis=0)
    shift = float(np.linalg.norm(cent(f"{rec} and resi 538-547 and name CA") - cent("ant and polymer and resi 538-547 and name CA")))
    side_camera(f"{h12} and name CA", oht, f"{oht} and elem N", distance=62)
    labels = render(png, {"closed H12 (agonist position)": (f"{h12} and name CA", (-0.7, -1.0)),
                          f"{ha.resn.capitalize()}{ha.resi}": (f"{h12} and resi {ha.resi} and not name N+C+O", (1.0, -0.8)),
                          "4-hydroxytamoxifen (antagonist)": (f"{oht} and not name {near}", (-1.0, 0.15)),
                          "diethylstilbestrol (agonist)": (des, (1.0, 0.35))})
    stats = {"antagonist_atom": la.name, "agonist_h12_residue": f"{ha.resn}{ha.resi}", "agonist_h12_atom": ha.name,
             "closest_approach_A": float(d[i, j]), "n_ligand_atoms_within_3A_of_agonist_h12": int((d.min(axis=1) < 3.0).sum()),
             "h12_centre_shift_A": shift}
    return stats, labels


def displacement_view(png: Path) -> dict:
    """3ERD chain A coloured by the median ESR1-ESR2 C-alpha distance at each position (blue 0 -> red 3 A), with the
    two pocket selectivity positions as spheres (data/esr2/esr1_esr2_displacement_profile.csv)."""
    prof = pd.read_csv(cc.DATA / "esr2/esr1_esr2_displacement_profile.csv").set_index("esr1_position")
    style_scene()
    rec = f"ago and chain {AGO[1]} and polymer"
    cmd.alter(rec, "b=-1.0")
    for pos, r in prof.iterrows():
        cmd.alter(f"{rec} and resi {pos}", f"b={min(float(r.median_displacement_A), 3.0)}")
    cmd.show("cartoon", rec)
    cmd.set("cartoon_transparency", 0.0, rec)  # object settings persist from the previous scene
    cmd.spectrum("b", "blue_white_red", f"{rec} and b > -0.5", minimum=0, maximum=3.0)
    cmd.color("grey70", f"{rec} and b < -0.5")
    cmd.set("cartoon_tube_radius", 0.5)
    sel = f"{rec} and resi 384+421 and not name N+C+O"
    cmd.show("spheres", sel)
    cmd.set("sphere_transparency", 0.0, sel)
    cmd.set_view(set_camera(f"ago and chain {AGO[1]} and polymer"))
    sel384, sel421 = f"{rec} and resi 384 and not name N+C+O", f"{rec} and resi 421 and not name N+C+O"
    left = project(sel384, 1600, 1300)[0] < project(sel421, 1600, 1300)[0]  # label each position on its own side
    lab = {f"L384 / M336 ({prof.loc[384].median_displacement_A:.2f} Å)": (sel384, (-1.0, -0.6) if left else (0.3, -1.0)),
           f"M421 / I373 ({prof.loc[421].median_displacement_A:.2f} Å)": (sel421, (0.3, -1.0) if left else (-1.0, -0.6)),
           "H12": (f"{rec} and resi 538-547 and name CA", (1.0, 0.6))}
    return render(png, lab)


def panels(images, titles, labels, out: Path, caption: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    fig, axs = plt.subplots(1, len(images), figsize=(6.4 * len(images), 5.6))
    for ax, img, ttl, lab in zip(np.atleast_1d(axs), images, titles, labels):
        im = Image.open(img)
        ax.imshow(im)
        for text, (x, y, d) in lab.items():
            u = np.array(d) / np.linalg.norm(d)
            tx, ty = np.array([x, y]) + u * 300
            ha = "right" if u[0] < -0.3 else "left" if u[0] > 0.3 else "center"
            ax.annotate(text, (x, y), (min(max(tx, 60), im.width - 60), min(max(ty, 40), im.height - 40)),
                        ha=ha, va="center", fontsize=12, color="#0b0b0b",
                        arrowprops=dict(arrowstyle="-", color="#52514e", lw=1.0, shrinkA=4, shrinkB=6))
        ax.set_title(ttl, fontsize=13, fontweight="bold", loc="left")
        ax.axis("off")
    fig.savefig(OUT / f"{out.stem}_untitled.png", dpi=200, bbox_inches="tight", facecolor="white")  # slide version
    fig.suptitle(textwrap.fill(caption, 75 * len(images)), x=0.01, y=1.0, ha="left", va="bottom", fontsize=12.5,
                 fontweight="bold")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rms, n_fit = load_pair()
    camera = set_camera(f"ago and chain {AGO[1]} and polymer")
    a_png, b_png, c_png = OUT / "agonist_3ERD_DES.png", OUT / "antagonist_3ERT_OHT.png", OUT / "h12_blocked.png"
    lab_a = complex_view("ago", AGO[1], AGO[2], BLUE, GREEN, AGO[3], a_png, camera)
    lab_b = complex_view("ant", ANT[1], ANT[2], ORANGE, CRIMSON, None, b_png, camera)
    lab_a["diethylstilbestrol"] = lab_a.pop("ligand")
    lab_b["4-hydroxytamoxifen"] = lab_b.pop("ligand")
    stats, lab_c = blocked_view(c_png)
    stats.update({"agonist": "3ERD A DES", "antagonist": "3ERT A OHT", "core_fit_rmsd_A": float(rms),
                  "core_fit_atoms": int(n_fit)})
    pd.DataFrame([stats]).to_csv(cc.TABLES_DIR / "pymol_contacts.csv", index=False)
    panels([a_png, b_png], ["Agonist: diethylstilbestrol (3ERD)", "Antagonist: 4-hydroxytamoxifen (3ERT)"],
           [lab_a, lab_b], FIGS / "pymol_agonist_antagonist.png",
           "With the agonist, H12 closes over the pocket and holds the coactivator peptide; with the antagonist, H12 "
           f"sits {stats['h12_centre_shift_A']:.0f} Å away in the coactivator groove.")
    res = stats["agonist_h12_residue"].capitalize()
    panels([c_png], ["Antagonist from 3ERT placed in the agonist structure 3ERD (core superposition)"], [lab_c],
           FIGS / "pymol_h12_blocked.png",
           f"The antagonist side chain overlaps {res} of the closed H12 ({stats['closest_approach_A']:.1f} Å; "
           f"{stats['n_ligand_atoms_within_3A_of_agonist_h12']} atoms within 3 Å): H12 cannot close over it.")
    d_png = OUT / "esr1_esr2_displacement.png"
    lab_d = displacement_view(d_png)
    panels([d_png], ["ESR1 (3ERD) coloured by ESR1-ESR2 C-alpha distance: blue 0 Å, white 1.5 Å, red ≥ 3 Å"], [lab_d],
           FIGS / "pymol_esr2_displacement.png",
           "Over 135 matched pairs the ESR1 and ESR2 backbones differ mainly in loops and at the H12 end; the two pocket "
           "selectivity positions barely move.")
    print(stats)


if __name__ == "__main__":
    main()
