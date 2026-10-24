"""Tracker Y44: crystal contacts at H12 - does packing hold H12 in place?

For every matrix chain: the deposited mmCIF (data/raw) is read with gemmi, waters and ligands removed (polymer atoms
only), a NeighborSearch with crystal symmetry is built (5 A cells), and for each atom of the chain's modelled H12
residues (ESR1 538-547, author numbering from results/lbd_extract/residue_mapping.csv) the polymer atoms within 4.0 A
that do NOT belong to the same chain in the same asymmetric unit are collected. Each contact atom is classified by its
partner: 'peptide' (non-receptor polymer of <= peptide_max_length residues), 'receptor ASU' (another ESR1/ESR2 chain of
the same asymmetric unit) or 'symmetry' (any chain of a symmetry-related copy, including the chain's own image).
Per chain: h12_contact_atoms (all partners except the peptide), h12_contact_residues (distinct H12 residues touching),
contacts_with_symmetry, contacts_receptor_asu, contacts_peptide, plus the number of modelled H12 residues.
Tests (chain level, descriptive; entry-level medians reported too): cl7 vs cl4; within cl4 complete vs partial H12;
Spearman of contacts vs H12 shift within cl4.
Outputs: results/ext/h12_crystal_contacts.csv, figures/ext/Y44_h12_contacts.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_h12_contacts.py
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK, INK2, ORANGE, RES, ROOT, SETTINGS, VIOLET, YELLOW, chain_table, mpl_setup,
                        save_fig, write_numbers)

H12 = range(SETTINGS["h12_range"][0], SETTINGS["h12_range"][1] + 1)


def main():
    t = chain_table()
    rc = pd.read_csv(ROOT / "results/lbd_extract/receptor_chains.csv", low_memory=False)
    rc = rc[rc.model == 1] if "model" in rc else rc
    lc = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False)[["chain_key", "auth_asym_id"]]
    t = t.merge(lc, on="chain_key", how="left")
    mp = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_residue", "auth_seq_id"])
    mp = mp[mp.chain_key.isin(t.chain_key) & mp.has_residue & mp.esr1_position.isin(list(H12))]
    h12res = mp.groupby("chain_key").auth_seq_id.apply(lambda s: set(int(x) for x in s.dropna()))
    rows = []
    for ent, g in t.groupby("entry_id"):
        st = gemmi.read_structure(str(ROOT / f"data/raw/{ent}.cif"))
        st.setup_entities()
        st.remove_ligands_and_waters()
        st.remove_hydrogens()
        model = st[0]
        frame_note = "standard"
        try:
            ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        except RuntimeError:
            st.standardize_crystal_frame()  # non-standard SCALEn: move coordinates to the standard frame (rigid change)
            model = st[0]
            ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
            frame_note = "standardized"
        r = rc[rc.pdb_id == ent]
        kind = {}
        for x in r.itertuples():
            if x.identity in ("esr1", "esr2"):
                kind[x.auth_asym_id] = "receptor"
            elif str(x.polymer_type).startswith("polypeptide") and x.sequence_length <= SETTINGS["peptide_max_length"]:
                kind[x.auth_asym_id] = "peptide"
            else:
                kind[x.auth_asym_id] = "other"
        for c in g.itertuples():
            res = h12res.get(c.chain_key, set())
            counts = {"symmetry": set(), "receptor ASU": set(), "peptide": set(), "other": set()}
            touching = set()
            chain = model.find_chain(c.auth_asym_id) if c.auth_asym_id in [ch.name for ch in model] else None
            if chain is not None:
                for rr in chain:
                    if rr.seqid.num not in res:
                        continue
                    for at in rr:
                        for m in ns.find_atoms(at.pos, "\0", radius=4.0):
                            cra = m.to_cra(model)
                            same = cra.chain.name == c.auth_asym_id
                            if same and m.image_idx == 0:
                                continue
                            if m.image_idx != 0:
                                key = "symmetry" if kind.get(cra.chain.name) != "peptide" else "peptide"
                            else:
                                key = {"receptor": "receptor ASU", "peptide": "peptide"}.get(kind.get(cra.chain.name), "other")
                            counts[key].add((cra.chain.name, m.image_idx, cra.residue.seqid.num, cra.atom.name))
                            touching.add(rr.seqid.num)
            rows.append(dict(chain_key=c.chain_key, entry_id=ent, cluster=c.cluster, h12_label=c.h12_label,
                             n_h12_residues_modelled=len(res), contacts_with_symmetry=len(counts["symmetry"]),
                             contacts_receptor_asu=len(counts["receptor ASU"]), contacts_peptide=len(counts["peptide"]),
                             contacts_other=len(counts["other"]),
                             h12_contact_atoms=len(counts["symmetry"]) + len(counts["receptor ASU"]) + len(counts["other"]),
                             h12_contact_residues=len(touching), space_group=c.space_group, crystal_frame=frame_note))
    df = pd.DataFrame(rows)
    sf = pd.read_csv(RES / "structural_features.csv")[["chain_key", "h12_shift_A"]]
    df = df.merge(sf, on="chain_key", how="left")
    df.to_csv(RES / "h12_crystal_contacts.csv", index=False)
    a, b = df[df.cluster == 7], df[df.cluster == 4]
    p74 = mannwhitneyu(a.h12_contact_atoms, b.h12_contact_atoms).pvalue
    c4c, c4p = b[b.h12_label == "complete"], b[b.h12_label == "partial"]
    pcp = mannwhitneyu(c4c.h12_contact_atoms, c4p.h12_contact_atoms).pvalue
    rho = spearmanr(b.h12_contact_atoms, b.h12_shift_A, nan_policy="omit")[0]
    frac0 = {k: float((g.h12_contact_atoms == 0).mean()) for k, g in (("cl7", a), ("cl4", b))}

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.6))
    ax = axs[0]
    rng = np.random.default_rng(42)
    groups = [("cl7", a, BLUE), ("cl4 complete H12", c4c, ORANGE), ("cl4 partial H12", c4p, YELLOW)]
    for k, (nm, g, col) in enumerate(groups):
        ax.scatter(k + rng.uniform(-.2, .2, len(g)), g.h12_contact_atoms, s=6, c=col, lw=0, alpha=0.7)
        ax.hlines(g.h12_contact_atoms.median(), k - .3, k + .3, color=INK, lw=1.5)
    ax.set_xticks(range(3), [f"{nm}\nn={len(g)}" for nm, g, _ in groups])
    ax.set(ylabel="H12 atoms' crystal-contact partner atoms (< 4 A, no peptide)", yscale="symlog",
           title=f"Packing at H12: cl7 vs cl4 p={p74:.1g}; cl4 complete vs partial p={pcp:.1g}")
    ax = axs[1]
    stack = df.groupby("cluster")[["contacts_with_symmetry", "contacts_receptor_asu", "contacts_peptide"]].mean()
    stack = stack.loc[[7, 4]]
    left = np.zeros(2)
    for col, c in zip(stack.columns, (VIOLET, GREY, YELLOW)):
        ax.barh(["cl7", "cl4"], stack[col], left=left, color=c, label=col.replace("contacts_", "").replace("_", " "))
        left += stack[col].values
    ax.set(xlabel="mean partner atoms within 4 A of H12 per chain", title="Who touches H12?")
    ax.legend(fontsize=8)
    ax = axs[2]
    ax.scatter(b.h12_contact_atoms, b.h12_shift_A, s=10, c=ORANGE, lw=0)
    ax.set(xscale="symlog", xlabel="H12 crystal-contact atoms (cl4 chains)", ylabel="H12 shift from cl7 mean (A)",
           title=f"Within cl4: contacts vs H12 position, Spearman {rho:.2f}")
    fig.suptitle(f"Lattice contacts at H12: median {a.h12_contact_atoms.median():.0f} atoms in cl7 (the peptide sits on H12) vs "
                 f"{b.h12_contact_atoms.median():.0f} in cl4; in cl4 more packing goes with an ordered H12 and a larger H12 shift", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y44_h12_contacts", section="D. Outliers and quality", task="Y44",
             shows="crystal-contact atoms around H12 (symmetry mates, other ASU chains, peptide) by cluster and H12 label; contacts vs H12 shift in cl4",
             takeaway=(f"Median contact atoms cl7 {a.h12_contact_atoms.median():.0f} vs cl4 {b.h12_contact_atoms.median():.0f} (p={p74:.1g}); "
                       f"partial-H12 cl4 chains have {c4p.h12_contact_atoms.median():.0f} vs {c4c.h12_contact_atoms.median():.0f} for complete (p={pcp:.1g})."),
             source="results/ext/h12_crystal_contacts.csv", script="scripts/ext_h12_contacts.py")
    write_numbers("Y44", dict(median_cl7=float(a.h12_contact_atoms.median()), median_cl4=float(b.h12_contact_atoms.median()),
                              p_cl7_cl4=float(p74), median_cl4_complete=float(c4c.h12_contact_atoms.median()),
                              median_cl4_partial=float(c4p.h12_contact_atoms.median()), p_complete_partial=float(pcp),
                              rho_cl4_shift=float(rho), frac_no_contact=frac0,
                              mean_by_partner=stack.round(2).to_dict("index")))
    print(stack, p74, pcp, rho, frac0)


if __name__ == "__main__":
    main()
