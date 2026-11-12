"""Rebuild the evidence figures from evidence/data and re-check every number in evidence/INDEX.md.

    python evidence/rebuild.py            # figures from evidence/data -> evidence/figures, then the number checks
    python evidence/rebuild.py --refresh  # first copy the source tables from the analysis outputs into evidence/data

Each claim below carries the number as stated in the index and a function that recomputes it from the tables in
evidence/data. A claim passes when every recomputed value equals the stated one within its tolerance (the rounding
of the stated value). The script writes evidence/checks.csv and evidence/INDEX.md and exits with status 1 if any
check fails.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE / "data"
FIGS = ROOT / "figures"
FIG_SCRIPTS = ROOT / "pipeline/figures"

# evidence/data file <- analysis output it is copied from
SOURCES = {
    "lid_distance_per_chain.csv": "results/03_lid_outliers/lid_distance_per_chain.csv",
    "lid_distance_by_group.csv": "results/03_lid_outliers/lid_distance_by_group.csv",
    "lid_outliers.csv": "results/03_lid_outliers/lid_outliers.csv",
    "lid_spread_drivers.csv": "results/03_lid_outliers/lid_spread_drivers.csv",
    "enrichment_tests.csv": "results/04_outliers_consensus/enrichment_tests.csv",
    "maps_coords.csv": "results/05_maps/maps_coords.csv",
    "maps_quality.csv": "results/05_maps/maps_quality.csv",
    "maps_islands.csv": "results/05_maps/maps_islands.csv",
    "numbers_05_maps.json": "results/05_maps/numbers.json",
    "no_lid_k_scan.csv": "results/06_lid_hidden/no_lid_k_scan.csv",
    "no_lid_separation.csv": "results/06_lid_hidden/no_lid_separation.csv",
    "no_lid_residue_profile.csv": "results/06_lid_hidden/no_lid_residue_profile.csv",
    "numbers_06_lid_hidden.json": "results/06_lid_hidden/numbers.json",
    "geometry_free_features.csv": "results/07_geometry_free/geometry_free_features.csv",
    "geometry_free_per_chain.csv": "results/07_geometry_free/geometry_free_per_chain.csv",
    "geometry_free_kscan.csv": "results/07_geometry_free/geometry_free_kscan.csv",
    "additives_lattice_tests.csv": "results/08_crystal_environment/additives_lattice_tests.csv",
    "lattice_within_crystal.csv": "results/08_crystal_environment/lattice_within_crystal.csv",
    "midpoint_error_by_k.csv": "data/lid/midpoint_error_by_k.csv",
    "midpoint_clustering_ari.csv": "data/lid/midpoint_clustering_ari.csv",
    "h12_window_status.csv": "data/lid/h12_window_status.csv",
    "dihedral_summary.csv": "data/lid/dihedral_summary.csv",
    "lid_order_odds_ratios.csv": "data/lid/lid_order_odds_ratios.csv",
    "chain_class.csv": "data/ligands/chain_class.csv",
    "assay_spread.csv": "data/ligands/assay_spread.csv",
    "druglikeness_by_group.csv": "data/ligands/druglikeness_by_group.csv",
    "pocket_filling_tests.csv": "data/ligands/pocket_filling_tests.csv",
}



def D(name: str) -> pd.DataFrame:
    if name.endswith(".json"):
        return json.loads((DATA / name).read_text())
    return pd.read_csv(DATA / name, keep_default_na=False, na_values=[""])


def _lid_copies():
    g, fl, s = D("lid_distance_per_chain.csv"), D("lid_outliers.csv"), D("lid_distance_by_group.csv").set_index("group")
    flagged = set(fl.chain_key)
    ents = [e for e in fl.entry_id.unique() if (g.entry_id == e).sum() > 1]
    alone = 0
    for e in ents:
        o = g[(g.entry_id == e) & ~g.chain_key.isin(flagged)].lid_to_h11_A.dropna()
        ok = o.apply(lambda x: any(abs(x - s.loc[k, "median"]) <= s.loc[k, "ordinary_threshold"] for k in ("cl7", "cl4")))
        alone += bool(len(o) and ok.all())
    sp = g[g.entry_id == "9B25"].lid_to_h11_A
    return alone, len(ents), sp.min(), sp.max()


def _within(o, grp):
    w = D("lattice_within_crystal.csv")
    return w[(w.contact == "sym_atoms") & (w.outcome == o) & (w.group == grp)].iloc[0]


def _ff(feature):
    return D("geometry_free_features.csv").set_index("feature").loc[feature]


def _cls(cls, grp):
    cc = D("chain_class.csv")
    return cc[(cc.action_class == cls) & (cc.cluster == grp)].entry_id.nunique(), cc[cc.action_class == cls].entry_id.nunique()


def _incomplete(grp, col):
    st = D("h12_window_status.csv")
    st = st[(st.level == "chain") & (st.cluster_group == grp)]
    return int(st[st[col] != "complete"].n.sum())


def _h12_rho(grp):
    from scipy.stats import spearmanr
    p = D("geometry_free_per_chain.csv")
    e = p[p.group == grp].groupby("entry_id")[["len_H12", "mobility_h12_z"]].mean().dropna()
    return spearmanr(e.len_H12, e.mobility_h12_z)[0], len(e)


# id, theme, claim, number (as stated), values (stated), tolerances, check, figure, source, obtained, calculated, caveat
CLAIMS = [
    ("L1", "Lid position", "The lid-to-H11 distance separates the two groups with no overlap.",
     "cl7 12.48 Å, cl4 22.26 Å (medians)", (12.48, 22.26), 0.005,
     lambda: tuple(D("lid_distance_by_group.csv").set_index("group").loc[["cl7", "cl4"], "median"]),
     "lid_outliers", "lid_distance_by_group.csv",
     "H12 centroid to the H11 axis, core-v3 superposition (28 Oct helix geometry)",
     "median over chains with an H12 axis (431 cl7, 307 cl4)", "needs ≥ 6 helical H12 C-alpha; 40 chains have none"),
    ("L2", "Lid position", "Copies of one crystal differ more than chains of a group spread; this sets 'ordinary variation'.",
     "copies 0.91 / 1.08 Å vs SD 0.39 / 0.61 Å (cl7 / cl4)", (0.91, 1.08, 0.39, 0.61), 0.005,
     lambda: tuple(D("lid_distance_by_group.csv").set_index("group").loc[["cl7", "cl4"], c].iloc[i]
                   for c in ("copy_p95", "sd") for i in (0, 1)),
     "lid_outliers", "lid_distance_by_group.csv", "same distance, chains of one entry and one group paired",
     "95th percentile of |difference| between copies; SD within group", "copies sit in different lattice contacts"),
    ("L3", "Lid position", "15 chains exceed ordinary variation; the two cl7 chains named in the meeting sit furthest out.",
     "15 chains; 5KCU B +1.53 Å, 5KD9 B +1.65 Å; cl4 Tukey outliers 0", (15, 1.53, 1.65, 0), 0.005,
     lambda: (int(D("lid_outliers.csv").beyond_ordinary.astype(str).eq("True").sum()),
              *D("lid_outliers.csv").set_index("chain_key").loc[["5KCU|B|1", "5KD9|B|1"], "deviation_A"],
              int(D("lid_distance_by_group.csv").set_index("group").loc["cl4", "n_tukey"])),
     "lid_outliers", "lid_outliers.csv", "Tukey fences per group; every grey chain; cl4 |robust z| ≥ 2",
     "|value − group median| > copy 95th pct (grey chains: vs nearer main group)", "cl4 extremes are not Tukey outliers"),
    ("L4", "Lid position", "Five grey chains hold the lid between the two positions.",
     "5 chains, 14.96-19.04 Å", (5, 14.96, 19.04), 0.005,
     lambda: (lambda f: (len(f), f.lid_to_h11_A.min(), f.lid_to_h11_A.max()))(
         D("lid_outliers.csv").query("group == 'minor' and lid_state == 'between cl7 and cl4'")),
     "lid_outliers", "lid_outliers.csv", "grey = the five minor clusters (14 chains)",
     "above cl7 median + 0.91 Å and below cl4 median − 1.08 Å", "1A52, 6PFM at ≥ 2.8 Å; most lack H12 545-547"),
    ("L5", "Lid position", "The unusual lid usually belongs to one copy, not to the crystal; 9B25 holds both positions.",
     "7 of 15 crystals; 9B25 copies 14.96-23.33 Å", (7, 15, 14.96, 23.33), 0.005, _lid_copies,
     "lid_outlier_copies", "lid_distance_per_chain.csv", "all copies of every crystal holding an examined chain",
     "other copies all within a main group's ordinary range", "copies differ in lattice environment"),
    ("L6", "Lid position", "Inside a group, no covariate shifts the lid beyond ordinary variation.",
     "0 real; largest ligand size in cl4 +0.77 Å (q 6e-7)", (0, 0.77), 0.005,
     lambda: (lambda d: (int(d.real.astype(str).eq("True").sum()), d.effect.abs().max()))(D("lid_spread_drivers.csv")),
     "lid_spread_drivers", "lid_spread_drivers.csv", "entry means within cl7 / cl4; peptide, construct, ligand, class, quality, contacts",
     "Mann-Whitney / Spearman, BH; effect = entry-median difference vs copy 95th pct", "detectable (q < 0.05) is not large"),
    ("L7", "Lid position", "Action classes nearly coincide with the lid groups, so class cannot be tested apart from lid.",
     "agonist 25 of 27 entries in cl7; SERM 19 of 20 in cl4", (25, 27, 19, 20), 0, lambda: (*_cls("agonist", "cl7"), *_cls("SERM", "cl4")),
     "class_vs_group", "chain_class.csv", "uncurated ChEMBL mechanism labels (4 Nov)",
     "entries with ≥ 1 chain of the class in the group", "51 of 362 entry × ligand rows labelled"),
    ("M1", "Maps", "The Sammon map keeps RMSD distances best; UMAP on PCA keeps neighbours.",
     "Sammon Spearman 0.92 (stress 0.041); UMAP-on-PCA trustworthiness 0.86", (0.92, 0.041, 0.86), 0.005,
     lambda: (lambda q, n: (q.loc["sammon", "spearman_dist"], n["sammon_stress"], q.loc["umap_on_pca", "trustworthiness"]))(
         D("maps_quality.csv").set_index("map"), D("numbers_05_maps.json")),
     "maps_quality", "maps_quality.csv", "Sammon on the 788 × 788 RMSD; SOM and UMAP on 6 PCA scores (95 % variance)",
     "trustworthiness (10 neighbours); Spearman of map vs RMSD distance", "SOM positions are a 12 × 12 grid"),
    ("M2", "Maps", "All three new maps show the same two lid groups.", "silhouette 0.76 / 0.45 / 0.76 (Sammon / SOM / UMAP)",
     (0.76, 0.45, 0.76), 0.005,
     lambda: tuple(D("maps_quality.csv").set_index("map").loc[["sammon", "som", "umap_on_pca"], "silhouette_cl7_cl4"]),
     "maps_three", "maps_quality.csv", "same three maps", "silhouette of cl7 vs cl4 in map coordinates", "—"),
    ("M3", "Maps", "No island the maps break into is set apart by more than copies of one crystal.",
     "largest gap 0.60 Å vs copies 1.26 Å", (0.60, 1.26), 0.005,
     lambda: (D("maps_islands.csv").gap_A.max(), D("maps_islands.csv").copy_rmsd_p95.iloc[0]),
     "map_islands", "maps_islands.csv", "HDBSCAN (≥ 10 chains) on each map",
     "median RMSD to other islands minus inside, vs copy RMSD 95th pct", "islands follow crystal form (UMAP ARI 0.68)"),
    ("N1", "Lid hidden", "Hiding residues 528-550 leaves one ensemble: no split exceeds copy variation.",
     "best split 784 / 4; gaps 0.26-0.71 Å vs copies 1.12 Å", (784, 4, 0.26, 0.71, 1.12), 0.005,
     lambda: (lambda k, s_: (*map(int, k.loc[k.silhouette.idxmax(), "sizes"].split(";")), s_.gap.min(), s_.gap.max(),
                             s_.copy_p95.iloc[0]))(D("no_lid_k_scan.csv"), D("no_lid_separation.csv")),
     "recover_lid_groups", "no_lid_separation.csv", "saved RMSD without 528-550 (Oct metric check)",
     "average linkage, K by silhouette; spectral K = 2; median between − within RMSD", "balanced splits overlap the lid groups (ARI ≤ 0.42)"),
    ("N2", "Lid hidden", "What remains is loop 331-340, where copies of one crystal differ as much.",
     "91 % of the squared difference; 1 residue (338) passes", (0.91, 1), 0.005,
     lambda: (D("numbers_06_lid_hidden.json")["splits"]["spectral_k2"]["loop_331_340_share"],
              len(D("numbers_06_lid_hidden.json")["splits"]["spectral_k2"]["residues_beyond_ordinary"])),
     "no_lid_profile", "no_lid_residue_profile.csv", "superposed C-alpha, group means of the spectral split",
     "per-residue distance vs within-group RMS and copy 95th pct", "loop partly unmodelled in many chains"),
    ("G1", "No superposition", "Only distances to H12 change between the groups; all other helix pairs stay put.",
     "8 of 36 real, all with H12; H7-H12 +5.45 Å; others ≤ 0.20 Å", (8, 5.45, 0.20), 0.005,
     lambda: (lambda f: (int(f.real.astype(str).eq("True").sum()), f.loc[f.feature == "d_H7_H12", "diff_cl4_minus_cl7"].iloc[0],
                         round(f[~f.feature.str.contains("H12")].diff_cl4_minus_cl7.abs().max(), 2)))(
         D("geometry_free_features.csv").query("kind == 'distance'")),
     "helix_pair_shift", "geometry_free_features.csv", "deposited coordinates, no fitting; mass-weighted backbone COM per helix and chain",
     "cl4 − cl7 median vs max(pooled SD, copy 95th pct)", "COM needs ≥ 70 % of a helix observed"),
    ("G2", "No superposition", "These features alone recover the lid groups, but only with H12.",
     "ARI 0.94 (K = 2); without H12 ≤ 0.20", (0.94, 0.20), 0.005,
     lambda: (lambda k: (k[(k.features == "all_features") & (k.k == 2)].ari_vs_rmsd_groups.iloc[0],
                         k[k.features == "without_H12"].ari_vs_rmsd_groups.max()))(D("geometry_free_kscan.csv")),
     "recover_lid_groups", "geometry_free_kscan.csv", "36 distances + 8 phi/psi helix lengths, z-scored",
     "Ward linkage, K = 2-10; ARI with the 1.5 Å RMSD groups", "763 chains with every feature"),
    ("G3", "No superposition", "Helix ends shift by a few residues, no more than copies of one crystal differ.",
     "H9 +5, H10/11 −4, H3 −3 residues", (5, -4, -3), 0,
     lambda: tuple(_ff(f"len_{h}").diff_cl4_minus_cl7 for h in ("H9", "H10/11", "H3")),
     "helix_lengths", "geometry_free_features.csv", "helical residue: phi −100..−30°, psi −80..−10°; runs ≥ 4",
     "cl4 − cl7 median length vs copy 95th pct (6, 4, 3)", "the box splits H4/5 at its kink"),
    ("G4", "No superposition", "In cl4 a longer H12 is a less mobile H12; in cl7 there is no link.",
     "ρ −0.37 (123 crystals) vs −0.05 (217)", (-0.37, 123, -0.05, 217), 0.005,
     lambda: (*_h12_rho("cl4"), *_h12_rho("cl7")),
     "h12_length_mobility", "geometry_free_per_chain.csv", "phi/psi H12 length; H12 C-alpha B z-score vs core (28 Oct)",
     "Spearman on entry means", "B also reflects refinement (TLS) and resolution"),
    ("C1", "Crystal environment", "Additives, recipe and solvent content leave lid order and mobility inside ordinary variation.",
     "0 of 22 tests real; largest 0.66 of ordinary", (0, 0.66), 0.005,
     lambda: (lambda t: (int(t.real.astype(str).eq("True").sum()), (t.effect / t.threshold).abs().max()))(
         D("additives_lattice_tests.csv").dropna(subset=["effect"])),
     "crystal_components", "additives_lattice_tests.csv", "mmCIF non-polymer copies, recipe text, solvent %, lattice contacts",
     "entry-level Spearman within groups, BH; effect vs copy 95th pct", "bound additives are rare (median 0 per chain)"),
    ("C2", "Crystal environment", "Within one crystal, the cl4 copy with more lattice contacts at the lid has more of it modelled.",
     "38 vs 14 crystals (q 0.0096)", (38, 14, 0.0096), 0.00005,
     lambda: (lambda r: (int(r.more_contacts_better), int(r.more_contacts_worse), r.q))(_within("lid_order", "cl4")),
     "lattice_within_crystal", "lattice_within_crystal.csv", "copies share drop, additives and solvent",
     "sign test over crystals whose copies differ in contacts", "contacts are counted on modelled atoms"),
    ("C3", "Crystal environment", "Between crystals, symmetry-mate contacts predict a complete H12 in cl4 (28 Oct).",
     "odds ratio 20.2 per SD [8.0-50.9]", (20.2, 8.0, 50.9), 0.05,
     lambda: tuple(D("lid_order_odds_ratios.csv").query(
         "outcome == 'complete' and role == 'primary' and cluster == 4 and level == 'chain' and predictor == 'lattice_sym'"
     )[["estimate", "ci_low", "ci_high"]].iloc[0]),
     "—", "lid_order_odds_ratios.csv", "logistic model within cl4, blocks quality / experiment / biology",
     "odds ratio per SD, entry-clustered errors", "partly mechanical: contacts need modelled atoms"),
    ("W1", "H12 window", "Ending the H12 window at 544 instead of 547 makes an incomplete lid rare (28 Oct).",
     "cl4 209 → 19 chains", (209, 19), 0,
     lambda: (_incomplete("cl4", "h12_status_538_547"), _incomplete("cl4", "h12_status_538_544")),
     "h12_window", "h12_window_status.csv", "consensus secondary structure gives 538-544",
     "chains with any C-alpha of the window unobserved", "window not yet signed off"),
    ("W2", "H12 window", "Residue 537 alone switches between the groups (28 Oct).",
     "psi +156.6° vs −46.8°, overlap 0", (156.6, -46.8, 0), 0.05,
     lambda: tuple(D("dihedral_summary.csv").set_index("position").loc[537, ["psi_mean_cl7", "psi_mean_cl4", "overlap_chain"]]),
     "—", "dihedral_summary.csv", "backbone dihedrals, no superposition", "circular mean per group; share of chains in the other's range",
     "537 carries the Y537S construct mutation"),
    ("W3", "Helix midpoints", "Midpoints tolerate missing ends and reproduce the groups (28 Oct).",
     "≤ 0.91 Å with 4 residues missing; ARI 0.954", (0.91, 0.954), 0.005,
     lambda: (D("midpoint_error_by_k.csv").query("k_removed == 4").median_midpoint_error_A.max(),
              D("midpoint_clustering_ari.csv").query("distance == 'midpoint_rms' and level == 'chain' and K == 2").ari.iloc[0]),
     "midpoint_error", "midpoint_error_by_k.csv", "complete helices truncated 1-6 residues", "median midpoint shift; RMS midpoint distance, average linkage",
     "agreement restates the H12 split"),
    ("O1", "Outliers", "Consensus RMSD outliers are rarer in cl7 than expected (28 Oct).",
     "2 vs 8.0 expected of 14 (q 0.0023)", (2, 8.0, 0.0023), 0.05,
     lambda: tuple(D("enrichment_tests.csv").query("level == 'chain' and test == 'cl7 vs cl4 (main only)'")[
         ["observed", "expected", "q_bh"]].iloc[0]),
     "—", "enrichment_tests.csv", "22 chains flagged by ≥ 3 of 5 detectors (Oct)", "hypergeometric, BH over 8 tests",
     "minor-group enrichment holds by construction"),
    ("P1", "Ligands", "Ligand size explains most drug-likeness differences between the groups.",
     "7 distinct quantities differ beyond size", (7,), 0,
     lambda: (int(D("druglikeness_by_group.csv").query("label == 'adds beyond size'").identical_to.isna().sum()),),
     "druglikeness_size", "druglikeness_by_group.csv", "RDKit descriptors of the crystallised ligands",
     "cl7 vs cl4 with heavy atoms as covariate, q < 0.05 at both levels", "oral drug-likeness, not binding"),
    ("P2", "Ligands", "Pocket filling does not differ between the groups (4 Nov).", "entry-level q 0.14", (0.14,), 0.005,
     lambda: (D("pocket_filling_tests.csv").query("comparison == 'cluster' and level == 'entry'").q_bh.iloc[0],),
     "—", "pocket_filling_tests.csv", "buriedness grid (0.4 Å), 701 chains with a complete ligand",
     "ligand / pocket volume, Mann-Whitney on entry medians", "pocket 100 % = wild-type E2 pocket"),
    ("P3", "Ligands", "Expressing affinity relative to E2 in the same assay makes assays agree better (4 Nov).",
     "4 of 5 ligands; median SD 0.86 → 0.56", (4, 0.86, 0.56), 0.005,
     lambda: (lambda a: (int((a.sd_dpx < a.sd_px_abs).sum()), a.sd_px_abs.median(), a.sd_dpx.median()))(D("assay_spread.csv")),
     "assay_spread", "assay_spread.csv", "ChEMBL binding activities on CHEMBL206", "SD across assays of pX vs dpX = pX − pX(E2)",
     "only 5 ligands have ≥ 2 paired assays"),
]


def refresh():
    DATA.mkdir(exist_ok=True)
    for dst, src in SOURCES.items():
        shutil.copyfile(ROOT / src, DATA / dst)
    print(f"copied {len(SOURCES)} tables into {DATA.relative_to(ROOT)}")


def figures():
    env = dict(os.environ, ESR1_DATA_DIR=str(DATA), ESR1_FIG_DIR=str(FIGS), MPLBACKEND="Agg")
    names = sorted({cl[7] for cl in CLAIMS if cl[7] != "—"})
    scripts = [FIG_SCRIPTS / f"fig_{n}.py" for n in names]
    for sc in scripts:
        r = subprocess.run([sys.executable, "-B", str(sc)], cwd=FIG_SCRIPTS, env=env, capture_output=True, text=True)
        if r.returncode:
            print(r.stderr)
            raise SystemExit(f"figure script failed: {sc.name}")
    print(f"rebuilt {len(scripts)} figures in {FIGS.relative_to(ROOT)}")


def check():
    rows, ok_all = [], True
    for cid, theme, claim, number, stated, tol, fn, fig, src, obt, calc, cav in CLAIMS:
        got = tuple(float(x) for x in fn())
        ok = len(got) == len(stated) and all(abs(g - s) <= tol + 1e-9 for g, s in zip(got, stated))
        ok_all &= ok
        rows.append(dict(id=cid, stated=";".join(map(str, stated)), recomputed=";".join(f"{g:.4g}" for g in got),
                         tolerance=tol, ok=ok))
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "checks.csv", index=False)
    print(out.to_string(index=False))
    return ok_all


def index():
    lines = ["# Evidence index", "",
             "Results added since the October report (28 Oct, 4 Nov and 11 Nov) that carry the conclusions. Every number is recomputed "
             "from the copies in `evidence/data/` by `python evidence/rebuild.py` (result in `checks.csv`), which also redraws the "
             "figures (`figures/`). "
             "'Ordinary variation' = the larger of the within-group SD and the 95th percentile of differences between "
             "copies of one crystal; an effect counts only above it.", "",
             "| # | Claim | Number | Figure | Source file | How obtained | How calculated | Caveat |",
             "|---|---|---|---|---|---|---|---|"]
    for cid, theme, claim, number, stated, tol, fn, fig, src, obt, calc, cav in CLAIMS:
        f = f"[{fig}](../figures/{fig}.png)" if fig != "—" else "—"
        lines.append(f"| {cid} | **{theme}.** {claim} | {number} | {f} | `{SOURCES[src]}` | {obt} | {calc} | {cav} |")
    lines += ["", "`evidence/data/` holds a copy of each source file under the same name; "
              "`python evidence/rebuild.py --refresh` copies them again from `results/` and `data/`."]
    (HERE / "INDEX.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    if "--refresh" in sys.argv:
        refresh()
    figures()
    passed = check()
    index()
    print("all numbers match" if passed else "SOME NUMBERS DO NOT MATCH")
    sys.exit(0 if passed else 1)
