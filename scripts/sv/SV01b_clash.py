"""SV01b step 8 (runs only if the re-issued SV01b verdict is "mixed"): medoid H12 clash test.

For every chain with a partial or not-modelled H12 (538-547; SV01a outcome, all clusters), the H12 of the cl7 medoid
5TM1|A|1 and of the cl4 medoid 7UJ7|C|1 (the complete-H12 medoids named in the run file) is superposed onto the
chain's core and its clashes with lattice neighbours are counted: does the crystal block one of the two positions?

Rules fixed before the test was run:
- Superposition: Kabsch (proper rotation) of the medoid C-alpha onto the target C-alpha over the core v3 positions
  (esr_analysis/resources/config/alignment_core_v3.json core_ranges) that both chains have, in the target's crystal
  frame as read for the neighbour search (lattice_contacts.load_structure). One conformer per atom: highest occupancy,
  then blank/A (read_structure policy). The direct fit is compared with the SV00 gate core-v3 transforms
  (data/sv/gate/alignment_transforms.json, medoid -> reference -> target) on the medoid H12 C-alpha (frame check).
- Moved atoms: every heavy atom of the medoid residues 538-547 (selected conformer).
- Lattice neighbours: heavy atoms of other chains' ASU copies and of every symmetry copy (the chain's own image
  included), the chain itself (ASU copy) and its own bound peptide copy left out (lattice_contacts.lattice_hits;
  same neighbour set as the step-6 lattice contacts). Waters and ligands are not in the structure.
- Clash: a moved heavy atom within 3.0 A of a neighbour heavy atom (CLASH_A). Counts: moved atoms in clash (all heavy
  atoms; backbone N, CA, C, O), split ASU / symmetry. [declared] A medoid position is 'blocked' in a chain when >= 1
  backbone atom of the moved H12 clashes.
Outputs: results/sv/SV01b/clash_test.csv (one row per chain x medoid), clash_summary.csv, clash_checks.json.
Run: .venv/bin/python -B scripts/sv/SV01b_clash.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common import lattice_contacts as lc  # noqa: E402

OUT = ROOT / "results/sv/SV01b"
H12 = (538, 547)
CLASH_A = 3.0
MEDOIDS = {"cl7": "5TM1|A|1", "cl4": "7UJ7|C|1"}


def kabsch(P, Q):
    """R, t with P @ R + t ~ Q (proper rotation)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, qc - pc @ R


def chain_atoms(model, auth, lookup):
    """{esr1_position: {atom name: xyz}} with one conformer per atom."""
    ch = lc.get_chain(model, auth)
    out = {}
    for res in ch:
        p = lookup.get((res.seqid.num, res.seqid.icode))
        if p is None:
            continue
        out[p] = {k: np.array(a.pos.tolist()) for k, a in lc.selected_atoms(res).items()}
    return out


def main():
    nums = json.loads((OUT / "numbers.json").read_text())
    if not nums["run_step8"]:
        (OUT / "clash_test_not_run.md").write_text(
            f"Step 8 not run: re-issued SV01b verdict is '{nums['step_verdict']}', not 'mixed'.\n")
        print("step 8 not run")
        return
    t = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    oc = pd.read_csv(ROOT / "results/sv/SV01a/chain_outcomes.csv")
    t = t.merge(oc[["chain_key", "h12_status", "cluster_group", "n_obs_window"]], on="chain_key")
    lat = pd.read_csv(OUT / "lattice_contacts_per_chain.csv", keep_default_na=False).set_index("chain_key")
    targets = t[t.h12_status.isin(["partial", "not_modelled"])].reset_index(drop=True)
    core_pos = set()
    for lo, hi in json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())[
            "core_ranges"]:
        core_pos |= set(range(lo, hi + 1))
    lookup = lc.residue_lookup(list(targets.chain_key) + list(MEDOIDS.values()))
    auth_of = t.set_index("chain_key").auth_asym_id
    med = {}
    for g, ck in MEDOIDS.items():
        model, _, _ = lc.load_structure(ck.split("|")[0])
        at = chain_atoms(model, auth_of[ck], lookup[ck])
        assert all(p in at for p in range(H12[0], H12[1] + 1)), ck
        med[g] = at
    tr = json.loads((ROOT / "data/sv/gate/alignment_transforms.json").read_text())

    def to_ref(ck, xyz):
        d = tr[ck]
        return (xyz - np.array(d["mobile_centroid"])) @ np.array(d["rotation"]) + np.array(d["reference_centroid"])

    def from_ref(ck, xyz):
        d = tr[ck]
        return (xyz - np.array(d["reference_centroid"])) @ np.array(d["rotation"]).T + np.array(d["mobile_centroid"])

    rows, frame_dev = [], []
    for ent, g in targets.groupby("entry_id", sort=False):
        model, ns, frame = lc.load_structure(ent)
        kinds = lc.chain_kinds(ent)
        for r in g.itertuples():
            tat = chain_atoms(model, r.auth_asym_id, lookup[r.chain_key])
            lr = lat.loc[r.chain_key]
            own = ((lr.own_peptide_auth_asym_id, int(float(lr.own_peptide_image_idx)), lr.own_peptide_source)
                   if lr.own_peptide_source != "none" else None)
            for mg, mat in med.items():
                common = sorted(p for p in core_pos if p in tat and "CA" in tat[p] and p in mat and "CA" in mat[p])
                P = np.array([mat[p]["CA"] for p in common])
                Q = np.array([tat[p]["CA"] for p in common])
                R, tvec = kabsch(P, Q)
                rmsd = float(np.sqrt(((P @ R + tvec - Q) ** 2).sum(1).mean()))
                moved = [(p, name, xyz @ R + tvec) for p in range(H12[0], H12[1] + 1) for name, xyz in mat[p].items()]
                if frame == "standard":
                    ca = np.array([mat[p]["CA"] for p in range(H12[0], H12[1] + 1)])
                    via_gate = from_ref(r.chain_key, to_ref(MEDOIDS[mg], ca))
                    direct = ca @ R + tvec
                    frame_dev.append(float(np.sqrt(((via_gate - direct) ** 2).sum(1)).max()))
                clash = {"all": set(), "bb": set(), "asu": set(), "sym": set(), "bb_asu": set(), "bb_sym": set()}
                partners, kinds_hit, dmin = set(), set(), np.inf
                for p, name, xyz in moved:
                    hits = lc.lattice_hits(model, ns, gemmi.Position(*xyz), r.auth_asym_id, kinds, exclude=own,
                                           cutoff=CLASH_A)
                    if not hits:
                        continue
                    a = (p, name)
                    clash["all"].add(a)
                    bb = name in lc.BACKBONE
                    if bb:
                        clash["bb"].add(a)
                    for h in hits:
                        src = "asu" if h[2] == "asu" else "sym"
                        clash[src].add(a)
                        if bb:
                            clash[f"bb_{src}"].add(a)
                        partners.add(h[0] if src == "asu" else f"{h[0]}@{h[1]}")
                        kinds_hit.add(h[3])
                        dmin = min(dmin, h[4])
                rows.append(dict(chain_key=r.chain_key, entry_id=ent, s08_cluster=int(r.s08_cluster),
                                 cluster_group=r.cluster_group, h12_status=r.h12_status,
                                 n_obs_538_547=int(r.n_obs_window), medoid=MEDOIDS[mg], medoid_cluster=mg,
                                 n_core_ca_fit=len(common), core_rmsd_A=rmsd, n_moved_atoms=len(moved),
                                 n_moved_backbone=sum(1 for _, n, _ in moved if n in lc.BACKBONE),
                                 clash_atoms=len(clash["all"]), clash_backbone_atoms=len(clash["bb"]),
                                 clash_atoms_asu=len(clash["asu"]), clash_atoms_sym=len(clash["sym"]),
                                 clash_backbone_asu=len(clash["bb_asu"]), clash_backbone_sym=len(clash["bb_sym"]),
                                 blocked=len(clash["bb"]) > 0, min_distance_A=dmin if np.isfinite(dmin) else np.nan,
                                 clash_partners=";".join(sorted(partners)), clash_partner_kinds=";".join(sorted(kinds_hit)),
                                 own_peptide_excluded=own is not None, crystal_frame=frame))
    ct = pd.DataFrame(rows)
    ct.to_csv(OUT / "clash_test.csv", index=False)

    # summary: per cluster group x H12 status x medoid, and the joint pattern per chain
    summ = []
    w = ct.pivot_table(index=["chain_key", "cluster_group", "h12_status"], columns="medoid_cluster",
                       values="blocked").reset_index()
    for (grp, st), d in list(ct.groupby(["cluster_group", "h12_status"])) + [(("all", "all"), ct)]:
        for mg in MEDOIDS:
            dd = d[d.medoid_cluster == mg]
            summ.append(dict(cluster_group=grp, h12_status=st, measure=f"{mg} medoid H12 blocked", n_chains=len(dd),
                             n_entries=dd.entry_id.nunique(), k=int(dd.blocked.sum()), fraction=dd.blocked.mean(),
                             entries_with_blocked_chain=dd[dd.blocked].entry_id.nunique(),
                             median_clash_atoms=dd.clash_atoms.median(),
                             median_clash_backbone_atoms=dd.clash_backbone_atoms.median()))
        ww = w if grp == "all" else w[(w.cluster_group == grp) & (w.h12_status == st)]
        for name, sel in (("only cl7 position blocked", (ww.cl7 == 1) & (ww.cl4 == 0)),
                          ("only cl4 position blocked", (ww.cl7 == 0) & (ww.cl4 == 1)),
                          ("both positions blocked", (ww.cl7 == 1) & (ww.cl4 == 1)),
                          ("neither position blocked", (ww.cl7 == 0) & (ww.cl4 == 0))):
            summ.append(dict(cluster_group=grp, h12_status=st, measure=name, n_chains=len(ww), k=int(sel.sum()),
                             fraction=sel.mean() if len(ww) else np.nan))
    sm_ = pd.DataFrame(summ)
    sm_.to_csv(OUT / "clash_summary.csv", index=False)
    checks = dict(n_targets=int(targets.shape[0]), n_rows=len(ct),
                  targets_by_group_status=targets.groupby(["cluster_group", "h12_status"]).size().rename(
                      "n").reset_index().to_dict("records"),
                  core_ca_fit_min=int(ct.n_core_ca_fit.min()), core_rmsd_median=float(ct.core_rmsd_A.median()),
                  core_rmsd_max=float(ct.core_rmsd_A.max()),
                  frame_check_direct_vs_gate_h12_ca_max_dev_A=dict(
                      n=len(frame_dev), median=float(np.median(frame_dev)), max=float(np.max(frame_dev))),
                  medoid_h12_atoms={g: sum(len(m[p]) for p in range(H12[0], H12[1] + 1)) for g, m in med.items()},
                  clash_cutoff_A=CLASH_A)
    (OUT / "clash_checks.json").write_text(json.dumps(checks, indent=1))
    print(json.dumps(checks, indent=1))
    print(sm_.to_string())


if __name__ == "__main__":
    main()
