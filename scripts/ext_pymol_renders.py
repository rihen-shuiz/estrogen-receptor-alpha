"""Tracker Y30: PyMOL renders of the two main-cluster medoids (through scripts/pymol_safe.py only).

Chains: the complete-H12 medoids of cl7 and cl4 (results/ext/cluster_displacement_profile.csv, Y10: 5TM1|A|1 and
7UJ7|C|1), loaded from data/raw mmCIFs (absolute paths from ESR1_ROOT), one flat object per chain, moved with their
saved S03 transforms (results/alignment/alignment_transforms.json; the analysis superposition, no refit).
Images: (1) superposition, cl7 blue / cl4 orange (shared palette), body semi-transparent, H12 (538-547, author numbering
via residue_mapping.csv) opaque; two views (front, rotated 90 deg about y). (2) B-factor putty of each medoid: C-alpha B
divided by the chain's median C-alpha B, putty radius and colour (light -> dark blue) by that normalised B.
Every render: pymol_safe.render (faulthandler, 150 s timeout, PNG check). If PyMOL fails, a matplotlib C-alpha trace of
the same superposed coordinates is drawn instead and the figure title says so.
Outputs: figures/ext/Y30_pymol_medoids.png, figures/ext/Y30_pymol_bputty.png (matplotlib composites of the PyMOL
panels, with legends), raw panels in figures/ext/pymol_panels/.
Run from the repo root: PYTHONPATH=. python scripts/ext_pymol_renders.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ext_common import BLUE, FIG, H12, INK, ORANGE, RES, ROOT, coords, mpl_setup, rmsd, save_fig
from pymol_figures import matrix_4x4, author_residues
from pymol_safe import render

PANELS = FIG / "pymol_panels"
PANELS.mkdir(parents=True, exist_ok=True)


def chain_info(keys):
    ch = pd.read_csv(ROOT / "results/lbd_extract/receptor_chains.csv", usecols=["chain_key", "label_asym_id",
                                                                              "auth_asym_id"]).set_index("chain_key")
    mp = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_ca", "auth_seq_id", "insertion_code", "B_factor_A2"])
    return ch.loc[keys], mp[mp.chain_key.isin(keys)]


def load_cmds(obj, key, ch, tr):
    pdb = key.split("|")[0]
    path = (ROOT / f"data/raw/{pdb}.cif").as_posix()
    c = ch.loc[key]
    return [f"load {path}, raw_{obj}", f"create {obj}, raw_{obj} and polymer and chain {c.auth_asym_id} and segi {c.label_asym_id}",
            f"delete raw_{obj}", f"PY:cmd.transform_object('{obj}', {matrix_4x4(tr[key])}, state=1, homogenous=1, transpose=0)"]


BASE = ["set cif_use_auth, 1", "set orthoscopic, 1", "set ray_shadows, 0", "set antialias, 2", "bg_color white",
        "set opaque_background, 1", "set cartoon_fancy_helices, 1"]


def crop(path):
    """Trim the white margin of a rendered panel (keeps 30 px)."""
    im = Image.open(path).convert("RGB")
    a = np.asarray(im)
    m = (a < 245).any(2)
    ys, xs = np.where(m)
    if len(ys):
        im = im.crop((max(xs.min() - 30, 0), max(ys.min() - 30, 0), min(xs.max() + 30, a.shape[1]), min(ys.max() + 30, a.shape[0])))
    return im


def main():
    prof = pd.read_csv(RES / "cluster_displacement_profile.csv")
    keys = {"cl7": prof.medoid_cl7.iloc[0], "cl4": prof.medoid_cl4.iloc[0]}
    tr = json.loads((ROOT / "results/alignment/alignment_transforms.json").read_text())
    ch, mp = chain_info(list(keys.values()))
    h12 = {k: author_residues(mp, v, range(H12[0], H12[1] + 1)) for k, v in keys.items()}
    hexc = {"cl7": "0x2a78d6", "cl4": "0xeb6834"}
    light = {"cl7": "0xb7d3f6", "cl4": "0xf7c4ae"}
    ok_all = True
    views = {"front": [], "side": ["turn y, 90"]}
    for vname, extra in views.items():
        cmds = list(BASE)
        for k, v in keys.items():
            cmds += load_cmds(k, v, ch, tr)
            cmds += [f"color {light[k]}, {k}", f"set cartoon_transparency, 0.55, {k}"]
            if h12[k]:
                sel = f"{k} and resi {'+'.join(h12[k].values())}"
                cmds += [f"create {k}_h12, {sel}", f"color {hexc[k]}, {k}_h12", f"set cartoon_transparency, 0, {k}_h12"]
        cmds += ["hide everything", "show cartoon", "orient cl7", "turn x, -90"] + extra + ["zoom all, 2"]
        ok = render(cmds, PANELS / f"medoids_{vname}.png", 1400, 1100)
        ok_all &= ok
    for k, v in keys.items():
        b = mp[(mp.chain_key == v) & mp.has_ca].dropna(subset=["auth_seq_id"])
        med = b.B_factor_A2.median()
        vals = {f"{int(a)}": float(round(float(x) / float(med), 3)) for a, x in zip(b.auth_seq_id, b.B_factor_A2)}
        other = "cl4" if k == "cl7" else "cl7"
        cmds = list(BASE) + load_cmds(k, v, ch, tr) + load_cmds(other, keys[other], ch, tr) + [
            f"PY:vals = {vals}",
            f"PY:cmd.alter('{k}', 'b = vals.get(resi, 1.0)', space={{'vals': vals}})",
            "set_color s1, [0.80, 0.89, 0.98]", "set_color s2, [0.22, 0.53, 0.90]", "set_color s3, [0.05, 0.21, 0.42]",
            "hide everything", f"show cartoon, {k}", "cartoon putty", "set cartoon_putty_transform, 0",
            "set cartoon_putty_scale_min, 0.5", "set cartoon_putty_scale_max, 3", "set cartoon_putty_radius, 0.35",
            f"spectrum b, s1 s2 s3, {k}, minimum=0.5, maximum=2.0", "orient cl7", "turn x, -90", f"zoom {k}, 2",
            f"delete {other}"]
        ok = render(cmds, PANELS / f"bputty_{k}.png", 1200, 1100)
        ok_all &= ok

    plt = mpl_setup()
    from matplotlib.patches import Patch
    if ok_all:
        fig, axs = plt.subplots(1, 2, figsize=(16, 7.2))
        for ax, v in zip(axs, ("front", "side")):
            ax.imshow(crop(PANELS / f"medoids_{v}.png"))
            ax.axis("off")
            ax.set_title(f"{v} view" + (" (rotated 90 deg about y)" if v == "side" else ""))
        fig.legend(handles=[Patch(color=BLUE, label=f"cl7 medoid {keys['cl7']} (H12 opaque)"),
                            Patch(color=ORANGE, label=f"cl4 medoid {keys['cl4']} (H12 opaque)")], loc="lower center", ncol=2)
        how = "PyMOL 3.1 (open source), headless via pymol_safe.py"
    else:
        x, pos = coords()
        _, idx = rmsd()
        fig, axs = plt.subplots(1, 2, figsize=(16, 7.2))
        for ax, (a, b) in zip(axs, ((0, 1), (2, 1))):
            for k, col in (("cl7", BLUE), ("cl4", ORANGE)):
                i = int(np.where(idx.chain_key == keys[k])[0][0])
                ax.plot(x[i, :, a], x[i, :, b], color=col, lw=1, alpha=0.5)
                m = (pos >= H12[0]) & (pos <= H12[1])
                ax.plot(x[i, m, a], x[i, m, b], color=col, lw=4)
            ax.set_aspect("equal"); ax.axis("off")
        how = "FALLBACK: matplotlib C-alpha trace (PyMOL render failed; see figures/ext/pymol_panels/*.log)"
    fig.suptitle(f"The two states superposed on core v3: H12 of the cl7 medoid caps the pocket, H12 of the cl4 medoid lies "
                 f"~12 A away in the coactivator groove ({how})", x=0.01, ha="left", fontweight="bold", fontsize=11)
    save_fig(fig, "Y30_pymol_medoids", section="N. Structure images", task="Y30",
             shows="PyMOL cartoon of the complete-H12 medoids of cl7 and cl4 in the S03 frame, H12 opaque, two views",
             takeaway="Everything except H12 (and two loops) overlays; H12 occupies two different positions.",
             source="results/ext/cluster_displacement_profile.csv, results/alignment/alignment_transforms.json",
             script="scripts/ext_pymol_renders.py")
    if ok_all:
        fig, axs = plt.subplots(1, 2, figsize=(16, 7.5))
        for ax, k in zip(axs, ("cl7", "cl4")):
            ax.imshow(crop(PANELS / f"bputty_{k}.png"))
            ax.axis("off")
            ax.set_title(f"{k} medoid {keys[k]}")
        import matplotlib
        sm = matplotlib.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
            "p", [(0.80, 0.89, 0.98), (0.22, 0.53, 0.90), (0.05, 0.21, 0.42)]), norm=matplotlib.colors.Normalize(0.5, 2.0))
        fig.colorbar(sm, ax=axs, orientation="horizontal", fraction=0.04, pad=0.02,
                     label="C-alpha B / chain median C-alpha B (putty radius and colour)")
        fig.suptitle("B-factor putty: the helical core is thin and pale in both states; H12, the H11-H12 loop and the H2-H3 "
                     "loop are thick (mobile) in both; same orientation (cl7 frame)", x=0.01, ha="left", fontweight="bold", fontsize=11)
        save_fig(fig, "Y30_pymol_bputty", section="N. Structure images", task="Y30",
                 shows="PyMOL B-factor putty of both complete-H12 medoids (normalised C-alpha B)",
                 takeaway="Mobility concentrates in H12 and the loops that also carry the cl7/cl4 difference.",
                 source="results/lbd_extract/residue_mapping.csv", script="scripts/ext_pymol_renders.py")
    print("PyMOL ok" if ok_all else "fallback used", keys, {k: list(v.values()) for k, v in h12.items()})


if __name__ == "__main__":
    main()
