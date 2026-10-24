"""Tracker Y09: physico-chemical descriptors and Murcko scaffolds of the 306 ligands of interest.

Input: results/ligands_of_interest.csv (C11 SMILES), results/ext/chem_labels_best.csv (Y08: chemistry cluster and the
structural majority of each ligand's chains: cl7 / cl4 when >= 80% of its chains are in that S08 cluster).
RDKit descriptors: MW, Crippen cLogP, TPSA, HBD, HBA (Lipinski), rotatable bonds, heavy atoms, rings, aromatic rings,
fraction sp3, formal charge, basic amine present (SMARTS for aliphatic amines not amide), phenol count, Murcko
scaffold and generic scaffold.
Unit of every test = ligand (not chain), so one ligand seen in 47 chains counts once. Mann-Whitney U cl7-only vs
cl4-only ligands with rank-biserial effect size and BH q over the descriptors.
Outputs: results/ext/ligand_descriptors.csv, results/ext/ligand_descriptor_tests.csv,
figures/ext/Y09_physchem_by_cluster.png, figures/ext/Y09_scaffolds.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_ligand_descriptors.py
"""
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import BLUE, GREY, INK2, ORANGE, RES, ROOT, bh, mpl_setup, save_fig, write_numbers

BASIC_AMINE = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N-a);!$(N-C=N);!$(N#*);!$(N=*)]([C;!$(C=*)])")
PHENOL = Chem.MolFromSmarts("[OX2H]c")


def describe(smi):
    m = Chem.MolFromSmiles(smi)
    scaf = MurckoScaffold.GetScaffoldForMol(m)
    try:
        gen = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaf))
    except Exception:
        gen = ""
    return dict(mw=Descriptors.MolWt(m), clogp=Crippen.MolLogP(m), tpsa=rdMolDescriptors.CalcTPSA(m),
                hbd=Lipinski.NumHDonors(m), hba=Lipinski.NumHAcceptors(m), rotb=Lipinski.NumRotatableBonds(m),
                heavy_atoms=m.GetNumHeavyAtoms(), rings=rdMolDescriptors.CalcNumRings(m),
                aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(m), fsp3=rdMolDescriptors.CalcFractionCSP3(m),
                formal_charge=Chem.GetFormalCharge(m), basic_amine=m.HasSubstructMatch(BASIC_AMINE),
                n_phenol=len(m.GetSubstructMatches(PHENOL)),
                murcko=Chem.MolToSmiles(scaf) or "acyclic", generic_scaffold=gen or "acyclic")


def main():
    lig = pd.read_csv(ROOT / "results/ligands_of_interest.csv")
    lab = pd.read_csv(RES / "chem_labels_best.csv")[["ligand_id", "chem_cluster_auto", "structural_majority"]]
    d = pd.DataFrame([describe(s) for s in lig.canonical_isomeric_smiles])
    d.insert(0, "ligand_id", lig.ligand_id)
    d.insert(1, "n_chains", lig.n_chains)
    d = d.merge(lab, on="ligand_id", how="left")
    d.to_csv(RES / "ligand_descriptors.csv", index=False)
    a, b = d[d.structural_majority == "cl7"], d[d.structural_majority == "cl4"]
    vars_ = ["mw", "clogp", "tpsa", "hbd", "hba", "rotb", "heavy_atoms", "rings", "aromatic_rings", "fsp3", "n_phenol"]
    rows = []
    for v in vars_:
        u, p = mannwhitneyu(a[v], b[v])
        rows.append(dict(variable=v, n_cl7=len(a), n_cl4=len(b), median_cl7=a[v].median(), median_cl4=b[v].median(),
                         U=u, rank_biserial=1 - 2 * u / (len(a) * len(b)), p_value=p))
    from scipy.stats import fisher_exact
    tab = [[a.basic_amine.sum(), (~a.basic_amine).sum()], [b.basic_amine.sum(), (~b.basic_amine).sum()]]
    orr, p = fisher_exact(tab)
    rows.append(dict(variable="basic_amine (Fisher)", n_cl7=len(a), n_cl4=len(b), median_cl7=a.basic_amine.mean(),
                     median_cl4=b.basic_amine.mean(), U=np.nan, rank_biserial=np.nan, p_value=p, odds_ratio=orr))
    t = pd.DataFrame(rows)
    t["q_value"] = bh(t.p_value)
    t.to_csv(RES / "ligand_descriptor_tests.csv", index=False)

    plt = mpl_setup()
    show = ["mw", "clogp", "tpsa", "heavy_atoms", "rotb", "hbd", "hba", "aromatic_rings"]
    units = {"mw": "g/mol", "clogp": "", "tpsa": "A^2", "heavy_atoms": "atoms", "rotb": "bonds", "hbd": "donors",
             "hba": "acceptors", "aromatic_rings": "rings"}
    fig, axs = plt.subplots(2, 5, figsize=(17, 7))
    axs = axs.ravel()
    rng = np.random.default_rng(42)
    for ax, v in zip(axs, show):
        for k, (g, col, name) in enumerate(((a, BLUE, "cl7"), (b, ORANGE, "cl4"))):
            parts = ax.violinplot([g[v].values], positions=[k], widths=0.75, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(col); pc.set_alpha(0.35); pc.set_edgecolor("none")
            ax.scatter(k + rng.uniform(-0.15, 0.15, len(g)), g[v], s=5, c=col, lw=0, alpha=0.7)
            ax.hlines(g[v].median(), k - 0.3, k + 0.3, color="black", lw=1.5)
        r = t.set_index("variable").loc[v]
        ax.set_xticks([0, 1], [f"cl7 ligands\nn={len(a)}", f"cl4 ligands\nn={len(b)}"])
        ax.set_ylabel(f"{v} ({units[v]})" if units[v] else v)
        ax.set_title(f"{v}: q={r.q_value:.1g}, r={r.rank_biserial:+.2f}", fontsize=9)
    ax = axs[8]
    fr = [a.basic_amine.mean() * 100, b.basic_amine.mean() * 100]
    ax.bar([0, 1], fr, color=[BLUE, ORANGE], width=0.6)
    for k, f in enumerate(fr):
        ax.text(k, f + 1, f"{f:.0f}%", ha="center", fontsize=9, color=INK2)
    rr = t[t.variable.str.startswith("basic")].iloc[0]
    ax.set_xticks([0, 1], ["cl7 ligands", "cl4 ligands"])
    ax.set(ylabel="ligands with a basic amine (%)", title=f"basic amine: Fisher q={rr.q_value:.1g}")
    ax = axs[9]
    sig = t.sort_values("q_value").head(8)
    ax.barh(sig.variable[::-1], -np.log10(sig.q_value[::-1]), color=GREY)
    ax.axvline(-np.log10(0.05), color="black", ls="--", lw=1)
    ax.set(xlabel="-log10 q (BH)", title="Strongest differences")
    fig.suptitle(f"Ligands of cl4 chains are larger, more flexible and carry a basic amine more often than ligands of "
                 f"cl7 chains (unit = ligand; cl7 n={len(a)}, cl4 n={len(b)})", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    mw = t.set_index("variable").loc["mw"]
    save_fig(fig, "Y09_physchem_by_cluster", section="B. Chemistry clustering", task="Y09",
             shows="physchem descriptor distributions of cl7-only vs cl4-only ligands with Mann-Whitney / Fisher BH q",
             takeaway=(f"Median MW {mw.median_cl4:.0f} vs {mw.median_cl7:.0f} g/mol (cl4 vs cl7, q={mw.q_value:.1g}); "
                       f"basic amine in {fr[1]:.0f}% of cl4 vs {fr[0]:.0f}% of cl7 ligands (q={rr.q_value:.1g})."),
             source="results/ext/ligand_descriptors.csv, results/ext/ligand_descriptor_tests.csv",
             script="scripts/ext_ligand_descriptors.py")

    # scaffolds
    sc = d[d.murcko != "acyclic"].groupby("murcko").agg(n=("ligand_id", "size"),
                                                       n7=("structural_majority", lambda s: (s == "cl7").sum()),
                                                       n4=("structural_majority", lambda s: (s == "cl4").sum()))
    sc = sc.sort_values("n", ascending=False).head(12)
    fig = plt.figure(figsize=(17, 8.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.5])
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(sc))[::-1]
    ax.barh(y, sc.n7, color=BLUE, label="ligands in cl7 chains")
    ax.barh(y, sc.n4, left=sc.n7, color=ORANGE, label="ligands in cl4 chains")
    ax.barh(y, sc.n - sc.n7 - sc.n4, left=sc.n7 + sc.n4, color=GREY, label="mixed / minor")
    ax.set_yticks(y, [f"S{i + 1}" for i in range(len(sc))])
    ax.set(xlabel="ligands with this Murcko scaffold", title=f"12 most common of {d.murcko.nunique()} Murcko scaffolds")
    ax.legend(fontsize=8, loc="lower right")
    ax2 = fig.add_subplot(gs[1])
    img = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in sc.index], molsPerRow=4, subImgSize=(260, 200),
                               legends=[f"S{i + 1}: n={n} (cl7 {a_}, cl4 {b_})" for i, (n, a_, b_) in
                                        enumerate(zip(sc.n, sc.n7, sc.n4))])
    from PIL import Image
    buf = io.BytesIO()
    img.save(buf, format="PNG") if hasattr(img, "save") else buf.write(img.data)
    ax2.imshow(Image.open(io.BytesIO(buf.getvalue())))
    ax2.axis("off")
    pure = int(((sc.n7 == sc.n) | (sc.n4 == sc.n)).sum())
    fig.suptitle(f"Scaffolds are cluster-specific: {pure} of the 12 most common scaffolds occur in only one structural "
                 "cluster", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    sc.reset_index().to_csv(RES / "ligand_top_scaffolds.csv", index=False)
    save_fig(fig, "Y09_scaffolds", section="B. Chemistry clustering", task="Y09",
             shows="12 most common Murcko scaffolds, counts split by the structural cluster of the ligand's chains",
             takeaway=f"{pure}/12 top scaffolds are confined to one structural cluster; {d.murcko.nunique()} distinct scaffolds in 306 ligands.",
             source="results/ext/ligand_descriptors.csv, results/ext/ligand_top_scaffolds.csv",
             script="scripts/ext_ligand_descriptors.py")
    write_numbers("Y09", dict(n_ligands=len(d), n_scaffolds=int(d.murcko.nunique()), n_cl7=len(a), n_cl4=len(b),
                              tests=t.set_index("variable")[["median_cl7", "median_cl4", "q_value"]].round(5).to_dict("index"),
                              top12_pure=pure))
    print(t.round(4).to_string())


if __name__ == "__main__":
    main()
