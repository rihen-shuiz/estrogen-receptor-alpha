"""Short deck: ligand action class and the position of helix 12 (results of the hand-labelled ligands, November 2026).

    python docs/slides_action_class.py     # -> docs/slides_action_class.pptx (with speaker notes) and .pdf

Charts are the figure scripts in pipeline/figures/ redrawn without their titles (ESR1_NO_TITLE=1, written to
results/12_curated_action/slides/); the PyMOL panels come from pipeline/curated_action/05_pymol_agonist_antagonist.py.
Every number on a slide is read from a saved table at build time:
  results/12_curated_action/tables/pairs_core3_with_chem.csv   class x lid position of the 291 labelled pairs
  results/12_curated_action/tables/effects_summary.csv         median shifts inside one lid position
  results/03_lid_outliers/lid_distance_by_group.csv            ordinary variation of the lid distance
  results/12_curated_action/tables/steroid_summary.csv         steroidal vs non-steroidal ligands
  results/12_curated_action/tables/redundancy_hill_summary.csv effective conformations per class
  results/12_curated_action/tables/pymol_contacts.csv          3ERD / 3ERT superposition and contacts
  results/12_curated_action/tables/join_audit_summary.json     labelled pairs, chains and entries
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from PIL import Image  # noqa: E402
from pptx import Presentation  # noqa: E402
from pptx.dml.color import RGBColor  # noqa: E402
from pptx.util import Inches, Pt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results/12_curated_action/tables"
FIGS = ROOT / "results/12_curated_action/slides"
PYMOL = ROOT / "results/12_curated_action/pymol"
OUT = ROOT / "docs/slides_action_class"
REPO = "github.com/rihen-shuiz/estrogen-receptor-alpha"
DATE = "26 November 2026"
AUTHORS = "Yerzhan Takauov · Ayan Duisenbekov · Dana Koshkinbayeva · Ayaulym Bekenova"
W, H = 13.333, 7.5
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
CHARTS = ["curated_class_map", "curated_h12_within", "curated_steroid_class"]


def csv(name):
    return pd.read_csv(T / name, keep_default_na=False, na_values=[""])


def draw_charts():
    """Redraw the charts used on slides without their figure titles."""
    env = {**os.environ, "ESR1_NO_TITLE": "1", "ESR1_FIG_DIR": str(FIGS), "MPLBACKEND": "Agg"}
    for name in CHARTS:
        subprocess.run([sys.executable, "-B", f"fig_{name}.py"], cwd=ROOT / "pipeline/figures", env=env, check=True,
                       stdout=subprocess.DEVNULL)


def numbers() -> dict:
    p = csv("pairs_core3_with_chem.csv")
    n = {c: int((p.action_class == c).sum()) for c in ("agonist", "partial_agonist", "antagonist")}
    k = {(c, g): int(((p.action_class == c) & (p.h12_group_excl == g)).sum()) for c in n for g in ("cl7", "cl4")}
    eff = csv("effects_summary.csv").set_index("contrast_id")
    ov = pd.read_csv(ROOT / "results/03_lid_outliers/lid_distance_by_group.csv").set_index("group").ordinary_threshold
    st = csv("steroid_summary.csv")
    sp = st[(st.level == "pairs_core3") & (st.steroid.astype(str) == "True")].iloc[0]
    so = st[(st.level == "pairs_core3") & (st.steroid.astype(str) == "False")].iloc[0]
    sl = st[(st.level == "unique_ligands") & (st.steroid.astype(str) == "True")].iloc[0]
    hill = csv("redundancy_hill_summary.csv").set_index(["subset", "cut_A"]).hill_q1
    pm = csv("pymol_contacts.csv").iloc[0]
    j = json.loads((T / "join_audit_summary.json").read_text())
    return dict(n=n, k=k, eff=eff, ov=ov, sp=sp, so=so, sl=sl, hill=hill, pm=pm, j=j)


def slides(v) -> list[dict]:
    n, k, eff, ov, sp, so, sl, hill, pm, j = (v[x] for x in ("n", "k", "eff", "ov", "sp", "so", "sl", "hill", "pm", "j"))
    a, b = eff.loc["cl7_Y537S_pep: ant vs ago"], eff.loc["cl7_Y537S_pep: part vs ago"]
    c4 = eff.loc["cl4_L536S_nopep: ant vs part"]
    sgn = lambda x: f"{x:+.2f}".replace("-", "−")
    return [
        dict(kind="title",
             notes="We linked a hand-assigned action class to every ligand in the ESR1 crystal structures and asked what "
                   "the class does to helix 12. The structures fall into two groups by the position of this helix, and "
                   "the question is whether function explains them."),
        dict(title="Agonists close the lid (H12); antagonists mostly leave it in the coactivator groove",
             images=[FIGS / "curated_class_map.png"],
             foot=f"Map of the 788 chains (UMAP of C-alpha coordinates); marker = hand-labelled class. Of "
                  f"{j['core3_strict_pairs']} PDB entry × ligand pairs, {k[('agonist', 'cl7')]}/{n['agonist']} agonist "
                  f"pairs are in cl7 (H12 over the pocket) and {k[('antagonist', 'cl4')]}/{n['antagonist']} antagonist "
                  f"pairs in cl4 (H12 in the groove); partial agonists {k[('partial_agonist', 'cl7')]} cl7, "
                  f"{k[('partial_agonist', 'cl4')]} cl4.",
             notes=f"ERα's ligand-binding domain ends in helix 12, the lid. In cl7 it closes over the ligand and builds "
                   f"the coactivator site; in cl4 it lies in that site. {k[('agonist', 'cl7')]} of {n['agonist']} "
                   f"agonist pairs are in cl7 and {k[('antagonist', 'cl4')]} of {n['antagonist']} antagonist pairs in "
                   f"cl4; partial agonists sit mostly in cl7 ({k[('partial_agonist', 'cl7')]} of "
                   f"{n['partial_agonist']})."),
        dict(title=f"Inside one lid position the class moves H12 by at most {abs(a.diff_a_minus_b):.2f} Å, "
                   "less than copies of one crystal differ",
             images=[FIGS / "curated_h12_within.png"],
             foot=f"Lid-to-H11 distance, median per pair. cl7: Y537S construct with peptide; cl4: L536S constructs "
                  f"without peptide. Dashed lines: ± ordinary variation ({ov['cl7']:.2f} / {ov['cl4']:.2f} Å).",
             notes=f"We compared classes only where constructs and peptide match. In cl7 antagonists sit "
                   f"{sgn(a.diff_a_minus_b)} Å and partial agonists {sgn(b.diff_a_minus_b)} Å from agonists; in cl4 "
                   f"antagonists and partial agonists differ by {sgn(c4.diff_a_minus_b)} Å. All of this is well inside "
                   f"the {ov['cl7']:.2f} Å that copies of one crystal already differ by."),
        dict(title="Steroidal ligands sit together in chemical space, and most of them are agonists, not antagonists",
             images=[FIGS / "curated_steroid_class.png"],
             foot=f"{int(sl.n)} of 240 labelled ligands carry the steroid four-ring nucleus and fall into one "
                  f"K-means cluster of the fingerprints. Antagonists: {int(sp.antagonist)} of {int(sp.n)} steroid pairs "
                  f"vs {int(so.antagonist)} of {int(so.n)} other pairs (Fisher p = {sp.antagonist_fisher_p:.0e}).".replace(
                      "e-0", "e-"),
             notes="Estradiol, estrone and estriol are steroids and agonists, so a steroid scaffold does not make an "
                   "antagonist. The three steroidal antagonists carry a long side chain on the steroid, like "
                   "fulvestrant or tamoxifen: it is the side chain that matters, as the next slides show."),
        dict(title="Agonist vs antagonist: H12 closes over the pocket or moves into the coactivator groove",
             images=[PYMOL / "pymol_agonist_antagonist_untitled.png"],
             foot=f"3ERD (diethylstilbestrol, GRIP1 peptide) and 3ERT (4-hydroxytamoxifen), wild type, Shiau et al. "
                  f"1998. Superposed on the core ({pm.core_fit_rmsd_A:.2f} Å over {int(pm.core_fit_atoms)} C-alpha); "
                  f"the H12 centre moves {pm.h12_centre_shift_A:.0f} Å.",
             notes="Same view for both structures after superposition on the stable core. With the agonist, H12 packs "
                   "over the pocket and the coactivator peptide binds next to it; with the antagonist, H12 lies in "
                   "the groove the peptide would use."),
        dict(title="The antagonist's side chain sits where the closed H12 would be",
             images=[PYMOL / "pymol_h12_blocked_untitled.png"],
             foot=f"4-hydroxytamoxifen from 3ERT placed into 3ERD: its side chain comes within "
                  f"{pm.closest_approach_A:.1f} Å of {pm.agonist_h12_residue.capitalize()} "
                  f"({int(pm.n_ligand_atoms_within_3A_of_agonist_h12)} atoms within 3 Å), so H12 cannot close. "
                  "Labels of less-studied ligands may be wrong; check them case by case.",
             notes="The antagonist protrudes from the pocket towards H12 and physically blocks the closed position, so "
                   "H12 has to move into the groove. This is the textbook mechanism; for 4-hydroxytamoxifen the label "
                   "is certain, but some of our labels come from assay thresholds and could be wrong."),
        dict(kind="text", title="What this gives the project",
             bullets=[f"H12 position reads out the action class: {k[('agonist', 'cl7')]}/{n['agonist']} agonist pairs "
                      f"closed, {k[('antagonist', 'cl4')]}/{n['antagonist']} antagonist pairs in the groove.",
                      f"It is a switch, not a gradient: inside one position the class leaves ≤ "
                      f"{abs(a.diff_a_minus_b):.2f} Å on H12, below the {ov['cl7']:.2f} Å crystal-copy noise.",
                      f"Redundancy: each class holds only {hill[('agonist', 1.5)]:.1f}–{hill[('antagonist', 1.5)]:.1f} "
                      "effective conformations; new structures of a known class rarely add a state.",
                      "A steroid scaffold does not make an antagonist; a side chain that blocks H12 does.",
                      f"Caveats: {k[('antagonist', 'cl7')]} antagonist pairs crystallise closed (mostly Y537S + peptide "
                      "constructs); some labels are derived from assay thresholds.",
                      "Next: check the closed-lid antagonists one by one; relate assay efficacy to H12 geometry; ESR2."],
             notes="The main gain is that function labels now explain the two structural groups. The open questions are "
                   "the exceptions, which mix crystal context with ligand action."),
    ]


def boxes(n_img: int, ratios: list[float]) -> list[tuple[float, float, float, float]]:
    """Image boxes (x, y from top, w, h) in inches inside the content area, aspect kept, centred."""
    x0, y0, aw, ah, gap = 0.5, 1.45, W - 1.0, 4.95, 0.3
    cw = (aw - gap * (n_img - 1)) / n_img
    out = []
    for i, r in enumerate(ratios):
        w = min(cw, ah * r)
        h = w / r
        out.append((x0 + i * (cw + gap) + (cw - w) / 2, y0 + (ah - h) / 2, w, h))
    return out


def build_pptx(sl, path: Path):
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    blank = prs.slide_layouts[6]

    def text(slide, x, y, w, h, s, size, bold=False, color=INK, lines=None):
        tf = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)).text_frame
        tf.word_wrap = True
        for i, line in enumerate(lines or [s]):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.text = line
            p.font.size, p.font.bold = Pt(size), bold
            p.font.color.rgb = RGBColor.from_string(color[1:])
            p.space_after = Pt(10 if lines else 0)

    for i, s in enumerate(sl, 1):
        slide = prs.slides.add_slide(blank)
        if s.get("kind") == "title":
            text(slide, 0.7, 2.0, 11.9, 1.6, "", 36, True,
                 lines=["ESR1 ligand-binding domain:", "what ligand action class does to helix 12"])
            text(slide, 0.7, 3.9, 11.9, 0.5, "Hand-labelled ligands linked to the crystal structures", 20, color=INK2)
            text(slide, 0.7, 4.5, 11.9, 0.5, f"{AUTHORS} · {DATE}", 16, color=INK2)
            text(slide, 0.7, 6.4, 11.9, 0.5, REPO, 14, color=MUTED)
            slide.notes_slide.notes_text_frame.text = s["notes"]
            continue
        text(slide, 0.5, 0.3, W - 1.0, 1.1, s["title"], 26, True)
        if s.get("kind") == "text":
            text(slide, 0.7, 1.5, W - 1.4, 5.3, "", 19, lines=[f"•  {b}" for b in s["bullets"]])
        else:
            ratios = [Image.open(p).width / Image.open(p).height for p in s["images"]]
            for p, (x, y, w, h) in zip(s["images"], boxes(len(ratios), ratios)):
                slide.shapes.add_picture(str(p), Inches(x), Inches(y), Inches(w), Inches(h))
            text(slide, 0.5, 6.5, W - 1.4, 0.8, s["foot"], 13, color=INK2)
        text(slide, W - 0.8, H - 0.45, 0.5, 0.3, str(i), 11, color=MUTED)
        slide.notes_slide.notes_text_frame.text = s["notes"]
    prs.core_properties.title = "ESR1: ligand action class and helix 12"
    prs.core_properties.author = ""
    prs.core_properties.comments = ""
    prs.core_properties.last_modified_by = ""
    prs.save(path)


def build_pdf(sl, path: Path):
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    with PdfPages(path, metadata={"Title": "ESR1: ligand action class and helix 12", "CreationDate": None}) as pdf:
        for i, s in enumerate(sl, 1):
            fig = plt.figure(figsize=(W, H))
            put = lambda x, y, t, **kw: fig.text(x / W, 1 - y / H, t, va="top", **kw)
            if s.get("kind") == "title":
                put(0.7, 2.0, "ESR1 ligand-binding domain:\nwhat ligand action class does to helix 12", fontsize=34,
                    fontweight="bold", color=INK, linespacing=1.2)
                put(0.7, 3.9, "Hand-labelled ligands linked to the crystal structures", fontsize=19, color=INK2)
                put(0.7, 4.5, f"{AUTHORS} · {DATE}", fontsize=15, color=INK2)
                put(0.7, 6.4, REPO, fontsize=13, color=MUTED)
            else:
                put(0.5, 0.3, textwrap.fill(s["title"], 70), fontsize=24, fontweight="bold", color=INK, linespacing=1.15)
                if s.get("kind") == "text":
                    put(0.7, 1.6, "\n\n".join(textwrap.fill(f"•  {b}", 88, subsequent_indent="    ")
                                               for b in s["bullets"]), fontsize=18, color=INK, linespacing=1.25)
                else:
                    ims = [Image.open(p) for p in s["images"]]
                    for im, (x, y, w, h) in zip(ims, boxes(len(ims), [m.width / m.height for m in ims])):
                        ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
                        ax.imshow(im)
                        ax.axis("off")
                    put(0.5, 6.5, textwrap.fill(s["foot"], 140), fontsize=12, color=INK2, linespacing=1.3)
                fig.text((W - 0.35) / W, 0.3 / H, str(i), ha="right", fontsize=11, color=MUTED)
            pdf.savefig(fig)
            plt.close(fig)


def main():
    draw_charts()
    sl = slides(numbers())
    build_pptx(sl, OUT.with_suffix(".pptx"))
    build_pdf(sl, OUT.with_suffix(".pdf"))
    (OUT.parent / f"{OUT.name}_notes.md").write_text(
        "# Speaker notes\n\n" + "\n\n".join(f"**Slide {i}.** {s['notes']}" for i, s in enumerate(sl, 1) if "notes" in s)
        + "\n")
    print("wrote", OUT.with_suffix(".pptx").name, OUT.with_suffix(".pdf").name, f"({len(sl)} slides)")


if __name__ == "__main__":
    main()
