"""SV02 steps 1-2 (RUN_SV1028): secondary structure of every kept chain from the SV00 backbone.

Tool: mkdssp is not installed on this machine (docs/sv/env_pre.txt), so the step's fallback pydssp 0.9.1 (numpy
backend, pydssp.assign(..., out_type='c3')) is used. pydssp returns H / E / '-' only. Its 'H' is the union of 3-, 4-
and 5-turn helices (pydssp_numpy.assign: helix = helix3 + helix4 + helix5), so 3_10 (G) and pi (I) cannot be
separated from alpha and are reported as not available.

Input: data/sv/backbone/backbone_atoms.csv.gz (SV00: N/CA/C/O of every residue of the 788 kept chains, altloc-selected,
usable occupancy, deposited frame, label_seq_id + esr1_position = UniProt P03372 number).
Full backbone: all residues of the chain in the table (also those outside 305-550 and residues without a UniProt
number) are passed to pydssp, so helix ends near 305 or 550 see their flanking residues.
Gaps: pydssp has no chain-break handling (the amide H of residue i is placed from C(i-1)), so each chain is split into
segments of consecutive label_seq_id residues with all four backbone atoms and C(i)-N(i+1) <= 2.5 A (mkdssp's peptide
bond limit, kMaxPeptideBondLength); pydssp runs once per segment. Proline residues get donor_mask 0 (no amide H), as
in DSSP. pydssp fails on segments shorter than 6 residues, so those are padded with dummy residues placed 1e4 A away
(no H-bond energy reaches them); the padding is dropped from the output. Strand bridges between different segments are
therefore not detected (E only; H and 'other' = not H are unaffected).
UniProt numbering: esr1_position of the SV00 table; checked against auth_seq_id + uniprot_author_offset of
results/sv/SV00/chain_offsets.csv for every residue.

Outputs: data/sv/SV02/dssp_residues.csv.gz (one row per residue passed to pydssp), data/sv/SV02/dssp_segments.csv,
results/sv/SV02/sse_codes_per_chain.csv (788 chains x positions 305-550: H / E / - ; empty = not observed with all
four backbone atoms), results/sv/SV02/dssp_run.json.
Run: .venv/bin/python -B -I scripts/sv/SV02_dssp.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pydssp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV02"
DATA = ROOT / "data/sv/SV02"
ATOMS = ["N", "CA", "C", "O"]
POS = np.arange(305, 551)
MAX_PEPTIDE_A = 2.5   # mkdssp kMaxPeptideBondLength
MIN_PYDSSP_LEN = 6    # pydssp_numpy.assign raises on shorter inputs


def dummy_residues(n, start_k):
    """n far-away residues with a sane internal geometry (no H-bond energy reaches real atoms or each other)."""
    tmpl = np.array([[0.0, 0.0, 0.0], [1.46, 0.0, 0.0], [2.0, 1.4, 0.0], [1.5, 2.5, 0.0]])
    return np.stack([tmpl + np.array([1e4 * (start_k + k + 1), 0.0, 0.0]) for k in range(n)])


def run_segment(xyz, donor):
    L = len(xyz)
    if L >= MIN_PYDSSP_LEN:
        return pydssp.assign(xyz, donor_mask=donor, out_type="c3"), 0
    pad = MIN_PYDSSP_LEN - L
    x = np.concatenate([xyz, dummy_residues(pad, 0)])
    d = np.concatenate([donor, np.zeros(pad)])
    return pydssp.assign(x, donor_mask=d, out_type="c3")[:L], pad


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    z = np.load(ROOT / "data/sv/backbone/backbone.npz")
    kept = [str(k) for k in z["chain_keys"]]
    b = pd.read_csv(ROOT / "data/sv/backbone/backbone_atoms.csv.gz", dtype={"ins_code": str, "altloc": str},
                    keep_default_na=False, na_values=[""])
    off = pd.read_csv(ROOT / "results/sv/SV00/chain_offsets.csv")[["chain_key", "uniprot_author_offset"]]

    # UniProt numbering check: esr1_position == auth_seq_id + SV00 offset
    m = b[b.esr1_position.notna()].merge(off, on="chain_key", how="left")
    numbering_equal = int(((m.auth_seq_id + m.uniprot_author_offset) == m.esr1_position).sum())
    numbering_differ = int(len(m) - numbering_equal)

    # one row per residue, coordinates of the four backbone atoms
    piv = b.pivot_table(index=["chain_key", "label_seq_id"], columns="atom", values=["x", "y", "z"], aggfunc="first")
    res = b.drop_duplicates(["chain_key", "label_seq_id"])[["chain_key", "entry_id", "label_seq_id", "esr1_position",
                                                             "auth_seq_id", "resname"]].set_index(["chain_key",
                                                                                                   "label_seq_id"])
    xyz = np.stack([np.stack([piv[(c, a)].to_numpy(float) for c in ("x", "y", "z")], -1) for a in ATOMS], 1)
    res = res.loc[piv.index].reset_index()
    full = np.isfinite(xyz).all(axis=(1, 2))
    res["all_backbone"] = full

    rows, segs = [], []
    n_padded_res = 0
    for key in kept:
        idx = np.where((res.chain_key == key).to_numpy() & full)[0]
        sub = res.iloc[idx].sort_values("label_seq_id")
        idx = sub.index.to_numpy()
        seq = sub.label_seq_id.to_numpy()
        X = xyz[idx]
        brk = np.ones(len(idx), bool)  # True = a segment starts here
        if len(idx) > 1:
            cn = np.linalg.norm(X[1:, 0] - X[:-1, 2], axis=1)
            brk[1:] = (np.diff(seq) != 1) | (cn > MAX_PEPTIDE_A)
        seg_id = np.cumsum(brk) - 1
        for s in np.unique(seg_id):
            sel = seg_id == s
            donor = (sub.resname.to_numpy()[sel] != "PRO").astype(float)
            sse, pad = run_segment(X[sel], donor)
            n_padded_res += int(sel.sum()) if pad else 0
            ss = sub[sel]
            segs.append({"chain_key": key, "segment": int(s), "label_seq_start": int(ss.label_seq_id.iloc[0]),
                         "label_seq_end": int(ss.label_seq_id.iloc[-1]), "esr1_start": ss.esr1_position.min(),
                         "esr1_end": ss.esr1_position.max(), "n_residues": int(sel.sum()), "padded": bool(pad)})
            for r, code in zip(ss.itertuples(), sse):
                rows.append({"chain_key": key, "entry_id": r.entry_id, "label_seq_id": int(r.label_seq_id),
                             "esr1_position": r.esr1_position, "auth_seq_id": r.auth_seq_id, "resname": r.resname,
                             "segment": int(s), "sse": str(code)})
    d = pd.DataFrame(rows)
    d.to_csv(DATA / "dssp_residues.csv.gz", index=False)
    sg = pd.DataFrame(segs)
    sg.to_csv(DATA / "dssp_segments.csv", index=False)

    # chain x position matrix 305-550 (empty = not observed with all four backbone atoms)
    lbd = d[d.esr1_position.between(305, 550)]
    if lbd.duplicated(["chain_key", "esr1_position"]).any():
        raise ValueError("duplicated UniProt position within a chain")
    mat = lbd.pivot(index="chain_key", columns="esr1_position", values="sse").reindex(index=kept, columns=POS)
    mat.columns = [str(int(c)) for c in mat.columns]
    mat.index.name = "chain_key"
    mat.to_csv(OUT / "sse_codes_per_chain.csv")

    # cross-check against the SV00 C-alpha observation mask
    mask = np.load(ROOT / "results/sv/SV00/observed_mask.npy")
    obs_dssp = mat.notna().to_numpy()
    info = {
        "tool": "pydssp 0.9.1 (numpy backend), fallback for mkdssp (not installed)",
        "chains": len(kept), "residues_passed": int(len(d)),
        "residues_in_table": int(len(res)), "residues_missing_backbone_atoms": int((~full).sum()),
        "segments": int(len(sg)), "segments_padded_lt6": int(sg.padded.sum()), "residues_in_padded_segments": n_padded_res,
        "segments_per_chain_median": float(sg.groupby("chain_key").size().median()),
        "segment_rule": f"consecutive label_seq_id, all of N/CA/C/O, C(i)-N(i+1) <= {MAX_PEPTIDE_A} A",
        "numbering_atoms_equal_auth_plus_offset": numbering_equal, "numbering_atoms_differ": numbering_differ,
        "lbd_positions_observed_dssp": int(obs_dssp.sum()), "lbd_positions_observed_ca_mask": int(mask.sum()),
        "ca_observed_but_incomplete_backbone": int((mask & ~obs_dssp).sum()),
        "dssp_observed_not_in_ca_mask": int((~mask & obs_dssp).sum()),
        "sse_counts_lbd": {k: int(v) for k, v in lbd.sse.value_counts().items()},
    }
    (OUT / "dssp_run.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
