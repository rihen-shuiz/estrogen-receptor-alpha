"""SV02 addition steps 7, 8, 10, 11 (RUN_SV1028): phi/psi per chain and the per-residue Ramachandran summary.

Step 7: phi and psi of every kept chain at every position 305-550, from the SV00 backbone
(data/sv/backbone/backbone_atoms.csv.gz, N/CA/C, UniProt numbering = esr1_position, which equals auth_seq_id + the
SV00 offset for every atom - checked in SV02_dssp.py), NumPy dihedral (IUPAC sign, degrees in (-180, 180]).
phi(i) = C(i-1)-N(i)-CA(i)-C(i) only if C(i-1)-N(i) <= 2.0 A; psi(i) = N(i)-CA(i)-C(i)-N(i+1) only if
C(i)-N(i+1) <= 2.0 A. Neighbours 304 and 551 are taken from the same table when present.
mkdssp did not run (not installed), so the step's mkdssp phi/psi comparison on 5 chains is not possible; as a
substitute sanity check (not the step's comparison) gemmi.calculate_phi_psi on the restored mmCIFs (data/raw, first
model; gemmi reads the first altloc, SV00 kept the highest-occupancy one, so residues with an altloc in residue
i-1, i or i+1 are reported apart) is compared on 5 chains fixed here: the S12/Y10 medoids 5TM1|A|1 and
7UJ7|C|1, the reference 7NFB|B|1, 1A52|A|1 (first chain of the matrix) and 6IAR|A|1 (the one struct_ref_seq offset).
Step 8 [DEFAULT] positions, fixed: 537 (main); 533-536, 538-540 (hinge); results/ext/numbers_Y12.json top_positions
other than 537 (329, 359, 457); 333 and 467 (SV07 spike positions named in the step).
Step 10: per position and cluster (cl7, cl4, minor = S08 clusters 1, 2, 3, 5, 6 pooled, all): chains and entries with
both angles; circular mean and circular SD (sqrt(-2 ln R)) of phi and psi; separation = sqrt(dphi^2 + dpsi^2)
between the cl7 and cl4 circular means, each difference wrapped to +-180; overlap = share of cl7 and cl4 units closer
(same wrapped distance) to the other cluster's mean than to their own - chain level (chain-level means) and entry level
(unit = entry x cluster, angles = circular mean over the entry's chains of that cluster; entry-level means).
Step 11: at 537 the circular-mean psi per cluster vs docs/EXPLORATION.md (+157 cl7, -47 cl4), tolerance 10 deg.
Outputs: results/sv/SV02/phipsi_per_chain.csv, rama_summary.csv, phipsi_check_gemmi.csv, numbers_rama.json.
Run: .venv/bin/python -B -I scripts/sv/SV02_phipsi.py
"""
import json
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV02"
POS = np.arange(305, 551)
BOND_A = 2.0
CHECK_CHAINS = ["5TM1|A|1", "7UJ7|C|1", "7NFB|B|1", "1A52|A|1", "6IAR|A|1"]
EXPECTED_537_PSI = {"cl7": 157.0, "cl4": -47.0}
TOL = 10.0


def positions_fixed():
    top = json.loads((ROOT / "results/ext/numbers_Y12.json").read_text())["top_positions"]
    hinge = [533, 534, 535, 536, 538, 539, 540]
    spikes = [333, 467]
    pos = [537] + hinge + [p for p in top if p != 537] + spikes
    assert len(pos) == len(set(pos)) == 13, pos
    return pos, top


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def wrap(d):
    return (np.asarray(d) + 180.0) % 360.0 - 180.0


def cmean(a):
    a = np.radians(np.asarray(a, float))
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean()))) if len(a) else np.nan


def csd(a):
    a = np.radians(np.asarray(a, float))
    if not len(a):
        return np.nan
    R = np.hypot(np.sin(a).mean(), np.cos(a).mean())
    return float(np.degrees(np.sqrt(-2.0 * np.log(min(R, 1.0)))))


def tdist(phi, psi, m):
    return np.hypot(wrap(np.asarray(phi) - m[0]), wrap(np.asarray(psi) - m[1]))


def cluster_label(c):
    return {7: "cl7", 4: "cl4"}.get(int(c), "minor")


def angles(b, kept):
    rng = np.arange(304, 552)
    row = {k: i for i, k in enumerate(kept)}
    X = np.full((len(kept), len(rng), 3, 3), np.nan)  # N, CA, C
    rn = np.full((len(kept), len(rng)), "", dtype=object)
    s = b[b.esr1_position.between(304, 551) & b.atom.isin(["N", "CA", "C"]) & b.chain_key.isin(row)]
    if s.duplicated(["chain_key", "esr1_position", "atom"]).any():
        raise ValueError("duplicated atom per chain and position")
    r = s.chain_key.map(row).to_numpy()
    p = (s.esr1_position.astype(int) - 304).to_numpy()
    a = s.atom.map({"N": 0, "CA": 1, "C": 2}).to_numpy()
    X[r, p, a] = s[["x", "y", "z"]].to_numpy(float)
    rn[r, p] = s.resname.to_numpy()
    N, CA, C = X[..., 0, :], X[..., 1, :], X[..., 2, :]
    cn_prev = np.linalg.norm(N[:, 1:] - C[:, :-1], axis=-1)  # C(i-1)-N(i) for i = 305..551
    phi = np.full((len(kept), len(rng)), np.nan)
    psi = np.full((len(kept), len(rng)), np.nan)
    with np.errstate(invalid="ignore"):
        ph = dihedral(C[:, :-1], N[:, 1:], CA[:, 1:], C[:, 1:])
        ok = cn_prev <= BOND_A
        phi[:, 1:] = np.where(ok, ph, np.nan)
        ps = dihedral(N[:, :-1], CA[:, :-1], C[:, :-1], N[:, 1:])
        psi[:, :-1] = np.where(ok, ps, np.nan)
    has_res = np.isfinite(X).all(axis=(2, 3))
    sl = slice(1, 1 + len(POS))  # 305-550
    return phi[:, sl], psi[:, sl], has_res[:, sl], rn[:, sl], cn_prev


def gemmi_check(df):
    rows = []
    for key in CHECK_CHAINS:
        pdb, label, model = key.split("|")
        st = gemmi.read_structure(str(ROOT / f"data/raw/{pdb}.cif"))
        st.setup_entities()
        mine = df[df.chain_key == key].set_index("auth_seq_id")
        for ch in st[int(model) - 1]:
            sub = ch.get_subchain(label)
            if len(sub) == 0:
                continue
            res = list(sub)
            for i, r in enumerate(res):
                prev = res[i - 1] if i > 0 and res[i - 1].seqid.num == r.seqid.num - 1 else None
                nxt = res[i + 1] if i + 1 < len(res) and res[i + 1].seqid.num == r.seqid.num + 1 else None
                ang = gemmi.calculate_phi_psi(prev, r, nxt)
                if r.seqid.num not in mine.index or r.seqid.icode.strip():
                    continue
                m = mine.loc[r.seqid.num]
                g_phi, g_psi = np.degrees(ang[0]), np.degrees(ang[1])
                rows.append({"chain_key": key, "auth_seq_id": r.seqid.num, "position": int(m.position),
                             "phi_numpy": m.phi, "phi_gemmi": g_phi, "psi_numpy": m.psi, "psi_gemmi": g_psi,
                             "has_altloc": any(at.has_altloc() for at in r),
                             "altloc_in_triplet": any(at.has_altloc() for x in (prev, r, nxt) if x is not None
                                                      for at in x)})
            break
    c = pd.DataFrame(rows)
    c["dphi"] = np.abs(wrap(c.phi_numpy - c.phi_gemmi))
    c["dpsi"] = np.abs(wrap(c.psi_numpy - c.psi_gemmi))
    c["only_numpy_phi"] = c.phi_numpy.notna() & c.phi_gemmi.isna()
    c["only_gemmi_phi"] = c.phi_numpy.isna() & c.phi_gemmi.notna()
    c["only_numpy_psi"] = c.psi_numpy.notna() & c.psi_gemmi.isna()
    c["only_gemmi_psi"] = c.psi_numpy.isna() & c.psi_gemmi.notna()
    return c


def summarise(df, positions, ht_label):
    rows = []
    for p in positions:
        d = df[(df.position == p) & df.phi.notna() & df.psi.notna()]
        rec = {"position": p, "context": ht_label(p)}
        for g in ("cl7", "cl4", "minor", "all"):
            x = d if g == "all" else d[d.cluster == g]
            rec[f"n_chains_{g}"] = len(x)
            rec[f"n_entries_{g}"] = x.entry_id.nunique()
            rec[f"phi_mean_{g}"] = cmean(x.phi)
            rec[f"phi_csd_{g}"] = csd(x.phi)
            rec[f"psi_mean_{g}"] = cmean(x.psi)
            rec[f"psi_csd_{g}"] = csd(x.psi)
        m7, m4 = (rec["phi_mean_cl7"], rec["psi_mean_cl7"]), (rec["phi_mean_cl4"], rec["psi_mean_cl4"])
        rec["separation_cl7_cl4"] = float(np.hypot(wrap(m7[0] - m4[0]), wrap(m7[1] - m4[1])))
        # chain level
        c7, c4 = d[d.cluster == "cl7"], d[d.cluster == "cl4"]
        o7 = tdist(c7.phi, c7.psi, m4) < tdist(c7.phi, c7.psi, m7)
        o4 = tdist(c4.phi, c4.psi, m7) < tdist(c4.phi, c4.psi, m4)
        rec["overlap_chain_cl7"] = float(o7.mean()) if len(o7) else np.nan
        rec["overlap_chain_cl4"] = float(o4.mean()) if len(o4) else np.nan
        rec["overlap_chain"] = float(np.concatenate([o7, o4]).mean()) if len(o7) + len(o4) else np.nan
        # entry level: unit = entry x cluster
        e = d[d.cluster.isin(["cl7", "cl4"])].groupby(["entry_id", "cluster"]).agg(
            phi=("phi", cmean), psi=("psi", cmean)).reset_index()
        e7, e4 = e[e.cluster == "cl7"], e[e.cluster == "cl4"]
        em7, em4 = (cmean(e7.phi), cmean(e7.psi)), (cmean(e4.phi), cmean(e4.psi))
        rec.update(phi_mean_entry_cl7=em7[0], psi_mean_entry_cl7=em7[1], phi_mean_entry_cl4=em4[0],
                   psi_mean_entry_cl4=em4[1], n_entry_units_cl7=len(e7), n_entry_units_cl4=len(e4))
        q7 = tdist(e7.phi, e7.psi, em4) < tdist(e7.phi, e7.psi, em7)
        q4 = tdist(e4.phi, e4.psi, em7) < tdist(e4.phi, e4.psi, em4)
        rec["overlap_entry_cl7"] = float(q7.mean()) if len(q7) else np.nan
        rec["overlap_entry_cl4"] = float(q4.mean()) if len(q4) else np.nan
        rec["overlap_entry"] = float(np.concatenate([q7, q4]).mean()) if len(q7) + len(q4) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def helix_context():
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from ext_common import HELIX_NAMES
    ht = pd.read_csv(ROOT / "results/helix_table.csv").sort_values("start")

    def nm(r):
        return HELIX_NAMES.get(r.helix, f"{r.helix} (unnamed)")

    def label(p):
        for r in ht.itertuples():
            if r.start <= p <= r.end:
                return f"{nm(r)} (helix table {r.helix} {r.start}-{r.end})"
        before = ht[ht.end < p].iloc[-1] if (ht.end < p).any() else None
        after = ht[ht.start > p].iloc[0] if (ht.start > p).any() else None
        b = nm(before) if before is not None else "N-terminus"
        a = nm(after) if after is not None else "C-terminus"
        return f"loop {b}-{a} (helix table)"
    return label


def main():
    positions, top = positions_fixed()
    z = np.load(ROOT / "data/sv/backbone/backbone.npz")
    kept = [str(k) for k in z["chain_keys"]]
    b = pd.read_csv(ROOT / "data/sv/backbone/backbone_atoms.csv.gz", dtype={"ins_code": str, "altloc": str},
                    keep_default_na=False, na_values=[""])
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv").set_index("chain_key")
    phi, psi, has, rn, cn_prev = angles(b, kept)

    auth = b[b.esr1_position.between(305, 550)].drop_duplicates(["chain_key", "esr1_position"]).set_index(
        ["chain_key", "esr1_position"]).auth_seq_id
    rows = []
    for i, k in enumerate(kept):
        for j, p in enumerate(POS):
            if not has[i, j]:
                continue
            rows.append({"chain_key": k, "entry_id": k.split("|")[0], "cluster": cluster_label(ct.s08_cluster[k]),
                         "s08_cluster": int(ct.s08_cluster[k]), "position": int(p), "resname": rn[i, j],
                         "auth_seq_id": int(auth[(k, float(p))]), "phi": phi[i, j], "psi": psi[i, j]})
    df = pd.DataFrame(rows)
    df.round({"phi": 3, "psi": 3}).to_csv(OUT / "phipsi_per_chain.csv", index=False)

    chk = gemmi_check(df)
    chk.round(4).to_csv(OUT / "phipsi_check_gemmi.csv", index=False)

    label = helix_context()
    rs = summarise(df, positions, label)
    rs.round(3).to_csv(OUT / "rama_summary.csv", index=False)

    r537 = rs[rs.position == 537].iloc[0]
    check = {}
    for g, exp in EXPECTED_537_PSI.items():
        got = float(r537[f"psi_mean_{g}"])
        diff = float(abs(wrap(got - exp)))
        check[g] = {"expected_psi": exp, "psi_circular_mean": round(got, 2), "abs_diff": round(diff, 2),
                    "within_10": diff <= TOL}
    num = {
        "positions": positions, "y12_top_positions": top,
        "rows_phipsi": len(df), "chains": df.chain_key.nunique(),
        "residues_with_phi": int(df.phi.notna().sum()), "residues_with_psi": int(df.psi.notna().sum()),
        "residues_with_both": int((df.phi.notna() & df.psi.notna()).sum()),
        "peptide_links_observed": int(np.isfinite(cn_prev).sum()),
        "peptide_links_gt_2A": int((cn_prev > BOND_A).sum()),
        "check_537_psi": check, "check_537_pass": all(v["within_10"] for v in check.values()),
        "gemmi_check": {"chains": CHECK_CHAINS, "residues_compared": len(chk),
                        "max_dphi": float(chk.dphi.max()), "max_dpsi": float(chk.dpsi.max()),
                        "residues_altloc_in_triplet": int(chk.altloc_in_triplet.sum()),
                        "max_dphi_no_altloc_in_triplet": float(chk[~chk.altloc_in_triplet].dphi.max()),
                        "max_dpsi_no_altloc_in_triplet": float(chk[~chk.altloc_in_triplet].dpsi.max()),
                        "only_numpy_phi": int(chk.only_numpy_phi.sum()), "only_gemmi_phi": int(chk.only_gemmi_phi.sum()),
                        "only_numpy_psi": int(chk.only_numpy_psi.sum()), "only_gemmi_psi": int(chk.only_gemmi_psi.sum())},
        "mkdssp_comparison": "not possible: mkdssp not installed (docs/sv/env_pre.txt)",
    }
    (OUT / "numbers_rama.json").write_text(json.dumps(num, indent=1))
    print(json.dumps(num, indent=1))
    pd.set_option("display.width", 250)
    print(rs[["position", "context", "n_chains_cl7", "n_chains_cl4", "n_chains_minor", "phi_mean_cl7", "psi_mean_cl7",
              "phi_mean_cl4", "psi_mean_cl4", "separation_cl7_cl4", "overlap_chain", "overlap_entry"]].round(1).to_string())


if __name__ == "__main__":
    main()
