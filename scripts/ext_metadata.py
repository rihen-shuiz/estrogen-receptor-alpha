"""Tracker Y02: experimental metadata from mmCIF for every ESR1 entry with a chain in the 788-chain matrix (351 entries).

Source order: the local mmCIF (data/raw/<id>.cif, sha256 in data/raw/sha256sums.txt); when a value is missing
there, the RCSB Data API (https://data.rcsb.org/rest/v1/core/entry/<id>) is asked once per entry and cached in
data/cache/rcsb_entry/<id>.json. A value that neither source reports is written as "not_reported" (text) or blank
(numbers) and the row is kept. Nothing is dropped.

Per entry: exp_method, cryst_method, ph, temp_K (crystal growth), cryst_details (raw _exptl_crystal_grow text),
precipitant flags (PEG, ammonium sulfate, MPD, other salt, DMSO, glycerol) and one additive_class, peg_mw, bound
additive HET codes (settings solvent_additive_codes found in the entry), space_group, crystal_system, cell a/b/c (A),
alpha/beta/gamma (deg), cell_volume_A3, chains_asu (ESR1 chains in the asymmetric unit), polymer_chains_asu,
release_year, deposit_year, expression_system, construct entity mutation text (_entity.pdbx_mutation), mutations
(from the residue-level mapping: results/lbd_chains.csv 'mismatches' with cause engineered mutation / conflict), named
hot-spot mutations (Y537S, D538G, L536S, Y537N, E380Q, S463P), peptide_present (any non-receptor polymer of
<= peptide_max_length residues), peptide_names, oligomer state (author/software assembly), r_work, completeness,
wilson_b, diffraction temperature, source (synchrotron / home), data_resolution.
Also a chain-level table results/ext/chain_metadata.csv (788 matrix chains, mutations per chain).
Figure: figures/ext/Y02_metadata_overview.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_metadata.py
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK2, ORANGE, RES, ROOT, SETTINGS, mpl_setup, save_fig, write_numbers,
                        chain_table)

CACHE = ROOT / "data/cache/rcsb_entry"
HOTSPOTS = ["Y537S", "D538G", "L536S", "Y537N", "Y537C", "E380Q", "S463P", "L536Q", "C381S", "C417S", "C530S"]
ADDITIVES = set(SETTINGS["solvent_additive_codes"])


def first(block, tag):
    v = block.find_values(tag)
    if len(v) == 0:
        return None
    s = v.str(0)
    return None if s in ("?", ".", "") else s


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def rcsb(entry):
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{entry}.json"
    if p.exists():
        return json.loads(p.read_text())
    try:
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{entry}", timeout=30)
        js = r.json() if r.ok else {}
    except Exception:
        js = {}
    p.write_text(json.dumps(js))
    time.sleep(0.1)
    return js


def precipitant_flags(text: str) -> dict:
    t = (text or "").upper()
    f = dict(
        has_peg=bool(re.search(r"PEG|POLYETHYLENE ?GLYCOL|PEG ?MME|P\.E\.G", t)),
        has_ammonium_sulfate=bool(re.search(r"AMMONIUM ?SULF|AMMONIUM ?SULPH|\(NH4\)2 ?SO4", t)),
        has_mpd=bool(re.search(r"\bMPD\b|METHYL-?2,4-?PENTANEDIOL", t)),
        has_dmso=bool(re.search(r"DMSO|DIMETHYL ?SULF|DIMETHYL ?SULPH", t)),
        has_glycerol=bool(re.search(r"GLYCEROL", t)),
        has_other_salt=bool(re.search(r"SODIUM|MAGNESIUM|LITHIUM|POTASSIUM|CALCIUM|AMMONIUM|CHLORIDE|FORMATE|CITRATE|"
                                      r"MALONATE|ACETATE|NACL|MGCL|LICL", t)),
    )
    m = re.search(r"PEG\s*(?:MME\s*)?-?\s*(\d{3,5})|POLYETHYLENE GLYCOL\s*(\d{3,5})", t)
    f["peg_mw"] = float(next(g for g in m.groups() if g)) if m else np.nan
    if not t.strip():
        cls = "not_reported"
    elif f["has_peg"]:
        cls = "PEG"
    elif f["has_ammonium_sulfate"]:
        cls = "ammonium_sulfate"
    elif f["has_mpd"]:
        cls = "MPD"
    elif f["has_other_salt"]:
        cls = "other_salt"
    else:
        cls = "other"
    f["additive_class"] = cls
    return f


def peptide_class(names: str) -> str:
    t = names.upper()
    if t == "NONE":
        return "none"
    if re.search(r"COACTIVATOR 2|NCOA2|GRIP|TIF2|TRANSCRIPTION INTERMEDIARY FACTOR-2|SRC-?2", t):
        return "NCOA2/SRC-2/GRIP1"
    if re.search(r"COACTIVATOR 1|NCOA1|SRC-?1", t):
        return "NCOA1/SRC-1"
    if re.search(r"COACTIVATOR 3|NCOA3|SRC-?3|AIB1", t):
        return "NCOA3/SRC-3"
    return "other peptide"


def crystal_system(sg: str | None) -> str:
    if not sg:
        return "not_reported"
    try:
        n = gemmi.SpaceGroup(sg).number
    except Exception:
        return "not_reported"
    for hi, name in ((2, "triclinic"), (15, "monoclinic"), (74, "orthorhombic"), (142, "tetragonal"),
                     (167, "trigonal"), (194, "hexagonal"), (230, "cubic")):
        if n <= hi:
            return name
    return "not_reported"


def parse_entry(path: Path, receptor_chains: pd.DataFrame, entry_ligands: str, receptor="esr1") -> dict:
    entry = path.stem.upper()
    b = gemmi.cif.read(str(path)).sole_block()
    details = first(b, "_exptl_crystal_grow.pdbx_details") or ""
    row = dict(entry_id=entry)
    row["exp_method"] = first(b, "_exptl.method") or "not_reported"
    row["cryst_method"] = (first(b, "_exptl_crystal_grow.method") or "not_reported").upper().strip()
    row["ph"] = num(first(b, "_exptl_crystal_grow.pH"))
    if np.isnan(row["ph"]):
        m = re.search(r"PH\s*=?\s*(\d{1,2}(?:\.\d+)?)", details.upper())
        row["ph"] = float(m.group(1)) if m and 2 <= float(m.group(1)) <= 12 else np.nan
        row["ph_source"] = "parsed from crystal-growth text" if m else "not_reported"
    else:
        row["ph_source"] = "_exptl_crystal_grow.pH"
    row["temp_K"] = num(first(b, "_exptl_crystal_grow.temp"))
    row["cryst_details"] = details.replace("\n", " ").strip() or "not_reported"
    row.update(precipitant_flags(details))
    sg = first(b, "_symmetry.space_group_name_H-M")
    row["space_group"] = sg or "not_reported"
    row["crystal_system"] = crystal_system(sg)
    for k in ("a", "b", "c"):
        row[f"cell_{k}"] = num(first(b, f"_cell.length_{k}"))
    for k in ("alpha", "beta", "gamma"):
        row[f"cell_{k}"] = num(first(b, f"_cell.angle_{k}"))
    try:
        row["cell_volume_A3"] = gemmi.UnitCell(row["cell_a"], row["cell_b"], row["cell_c"], row["cell_alpha"],
                                               row["cell_beta"], row["cell_gamma"]).volume
    except Exception:
        row["cell_volume_A3"] = np.nan
    rev = list(b.find_values("_pdbx_audit_revision_history.revision_date"))
    row["release_year"] = int(min(rev)[:4]) if rev else np.nan
    dep = first(b, "_pdbx_database_status.recvd_initial_deposition_date")
    row["deposit_year"] = int(dep[:4]) if dep else np.nan
    hosts = [h for h in b.find_values("_entity_src_gen.pdbx_host_org_scientific_name") if h not in ("?", ".")]
    row["expression_system"] = (gemmi.cif.as_string(hosts[0]) if hosts else "not_reported")
    # mutation text of the receptor entity
    rc = receptor_chains[receptor_chains.pdb_id == entry]
    rec = rc[rc.identity == receptor]
    rec_entities = set(rec.entity_id.astype(str))
    ent_ids = list(b.find_values("_entity.id"))
    muts = list(b.find_values("_entity.pdbx_mutation"))
    mtext = sorted({gemmi.cif.as_string(m) for e, m in zip(ent_ids, muts) if e in rec_entities and m not in ("?", ".")})
    row["entity_mutation_text"] = "; ".join(mtext) if mtext else "not_reported"
    row["chains_asu"] = int(len(rec))
    row["polymer_chains_asu"] = int(len(rc))
    pep = rc[(rc.identity != receptor) & (rc.polymer_type.astype(str).str.contains("polypeptide")) &
             (rc.sequence_length <= SETTINGS["peptide_max_length"])]
    row["peptide_present"] = bool(len(pep))
    row["n_peptide_chains"] = int(len(pep))
    row["peptide_names"] = "; ".join(sorted(set(pep.description.astype(str)))) if len(pep) else "none"
    row["peptide_class"] = peptide_class(row["peptide_names"])
    olig = list(b.find_values("_pdbx_struct_assembly.oligomeric_details"))
    row["oligomer_state"] = gemmi.cif.as_string(olig[0]).lower() if olig else "not_reported"
    row["r_work"] = num(first(b, "_refine.ls_R_factor_R_work"))
    row["r_free_cif"] = num(first(b, "_refine.ls_R_factor_R_free"))
    row["completeness_pct"] = num(first(b, "_reflns.percent_possible_obs") or first(b, "_refine.ls_percent_reflns_obs"))
    row["wilson_b"] = num(first(b, "_reflns.B_iso_Wilson_estimate"))
    row["diffraction_temp_K"] = num(first(b, "_diffrn.ambient_temp"))
    row["source"] = (first(b, "_diffrn_source.source") or "not_reported").upper()
    row["data_resolution_A"] = num(first(b, "_refine.ls_d_res_high"))
    codes = set(str(entry_ligands).split(";")) if isinstance(entry_ligands, str) else set()
    add = sorted(codes & ADDITIVES)
    row["bound_additives"] = ";".join(add) if add else "none"
    row["n_bound_additive_types"] = len(add)
    return row


def fill_from_rcsb(df):
    gaps = []
    for i, r in df.iterrows():
        need = [c for c in ("ph", "temp_K", "release_year") if pd.isna(r[c])]
        if r.expression_system == "not_reported":
            need.append("expression_system")
        if not need:
            continue
        js = rcsb(r.entry_id)
        g = (js.get("exptl_crystal_grow") or [{}])[0]
        if "ph" in need and g.get("p_h") is not None:
            df.at[i, "ph"] = g["p_h"]; df.at[i, "ph_source"] = "RCSB API exptl_crystal_grow.pH"
        if "temp_K" in need and g.get("temp") is not None:
            df.at[i, "temp_K"] = g["temp"]
        if "release_year" in need:
            d = (js.get("rcsb_accession_info") or {}).get("initial_release_date")
            if d:
                df.at[i, "release_year"] = int(d[:4])
        gaps.append(r.entry_id)
    return gaps


def chain_level(meta):
    t = chain_table()
    lc = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False)[["chain_key", "mismatches", "mismatch_count"]]
    t = t.merge(lc, on="chain_key", how="left", validate="one_to_one")

    def mut_list(s):
        out = []
        for item in str(s).split(";") if isinstance(s, str) else []:
            m = re.match(r"(\d+)([A-Z])>([A-Z]):(.*)", item)
            if m and m.group(3) != "X" and ("engineered" in m.group(4) or "conflict" in m.group(4) or "variant" in m.group(4)):
                out.append(f"{m.group(2)}{m.group(1)}{m.group(3)}")
        return out
    t["mutations"] = t.mismatches.map(lambda s: ";".join(mut_list(s)) or "none")
    for h in HOTSPOTS:
        t[f"mut_{h}"] = t.mutations.str.split(";").map(lambda L, h=h: h in L)
    lbd_hot = ["Y537S", "D538G", "L536S", "Y537N", "Y537C", "E380Q", "S463P", "L536Q"]
    t["activating_mutation"] = t[[f"mut_{h}" for h in lbd_hot]].any(axis=1)
    t["cys_triple_mutant"] = t[["mut_C381S", "mut_C417S", "mut_C530S"]].all(axis=1)
    def status(r):
        L = set(r.mutations.split(";")) if r.mutations != "none" else set()
        if L & {"Y537S", "Y537N", "Y537C"}:
            return "Y537S"
        if "D538G" in L:
            return "D538G"
        if L & {"L536S", "L536Q"}:
            return "L536S"
        if not L:
            return "wild type"
        if L <= {"C381S", "C417S", "C530S", "C530A"}:
            return "Cys->Ser/Ala only"
        return "other"
    t["mutation_status"] = t.apply(status, axis=1)
    cols = ["chain_key", "entry_id", "mutations", "mutation_status", "activating_mutation", "cys_triple_mutant"] + \
           [f"mut_{h}" for h in HOTSPOTS]
    m = t[cols].merge(meta[["entry_id", "peptide_present", "peptide_names", "peptide_class", "space_group", "additive_class", "ph",
                            "release_year", "exp_method", "cryst_method", "chains_asu"]], on="entry_id", how="left")
    return m


def figure(meta, cm):
    plt = mpl_setup()
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
    n = len(meta)
    ax = axs[0, 0]
    ph = meta.ph.dropna()
    ax.hist(ph, bins=np.arange(4, 10.01, 0.25), color=BLUE, edgecolor="white")
    ax.set(xlabel="crystallisation pH", ylabel="entries", title=f"pH reported for {len(ph)}/{n} entries")
    ax = axs[0, 1]
    vc = meta.additive_class.value_counts()
    ax.barh(vc.index[::-1], vc.values[::-1], color=[GREY if i == "not_reported" else BLUE for i in vc.index[::-1]])
    for y, v in enumerate(vc.values[::-1]):
        ax.text(v, y, f" {v}", va="center", fontsize=8, color=INK2)
    ax.set(xlabel="entries", title=f"Main precipitant class (n={n})")
    ax = axs[0, 2]
    vc = meta.space_group.value_counts()
    top = vc.head(8)
    other = vc.iloc[8:].sum()
    lab = list(top.index) + ([f"other ({len(vc) - 8} groups)"] if other else [])
    val = list(top.values) + ([other] if other else [])
    ax.barh(lab[::-1], val[::-1], color=[GREY if l.startswith("other") else BLUE for l in lab[::-1]])
    for y, v in enumerate(val[::-1]):
        ax.text(v, y, f" {v}", va="center", fontsize=8, color=INK2)
    ax.set(xlabel="entries", title=f"{len(vc)} space groups; top: {vc.index[0]}")
    ax = axs[1, 0]
    yr = meta.release_year.dropna().astype(int)
    ax.hist(yr, bins=np.arange(yr.min(), yr.max() + 2) - 0.5, color=BLUE, edgecolor="white")
    ax.set(xlabel="release year", ylabel="entries", title=f"Releases {yr.min()}-{yr.max()}")
    ax = axs[1, 1]
    pp = meta.peptide_present.mean() * 100
    vc = cm.mutation_status.value_counts()
    ax.barh(vc.index[::-1], vc.values[::-1], color=ORANGE)
    for y, v in enumerate(vc.values[::-1]):
        ax.text(v, y, f" {v}", va="center", fontsize=8, color=INK2)
    ax.set(xlabel="chains (of 788)", title="LBD mutation status per chain")
    ax = axs[1, 2]
    vc = meta.chains_asu.value_counts().sort_index()
    ax.bar(vc.index.astype(str), vc.values, color=BLUE)
    ax.set(xlabel="ESR1 chains in the asymmetric unit", ylabel="entries",
           title=f"Peptide in {pp:.0f}% of entries; ASU 1-{meta.chains_asu.max()} chains")
    fig.suptitle(f"Experimental metadata of the {n} ESR1 entries (local mmCIF, RCSB API for gaps)", x=0.01,
                 ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


def main():
    rc = pd.read_csv(ROOT / "results/lbd_extract/receptor_chains.csv", low_memory=False)
    ent = pd.read_csv(ROOT / "results/entries.csv")[["entry_id", "ligand_ids"]]
    lig = dict(zip(ent.entry_id, ent.ligand_ids))
    matrix_entries = set(pd.read_csv(ROOT / "results/rmsd_entry_index.csv").entry_id)  # the 351 entries with matrix chains
    rows = [parse_entry(p, rc, lig.get(p.stem.upper())) for p in sorted((ROOT / "data/raw").glob("*.cif"))
            if p.stem.upper() in matrix_entries]
    meta = pd.DataFrame(rows)
    gaps = fill_from_rcsb(meta)
    meta["ph_bin"] = pd.cut(meta.ph, [0, 6.5, 7.5, 8.5, 14], labels=["<6.5", "6.5-7.5", "7.5-8.5", ">8.5"]).astype(str)
    meta.loc[meta.ph.isna(), "ph_bin"] = "not_reported"
    meta.to_csv(RES / "entry_metadata.csv", index=False)
    cm = chain_level(meta)
    em = cm.groupby("entry_id").mutations.agg(lambda x: ";".join(sorted({m for v in x for m in v.split(";") if m != "none"})) or "none")
    ms = cm.groupby("entry_id").mutation_status.agg(lambda x: x.mode().iloc[0])
    meta["mutations"] = meta.entry_id.map(em).fillna("not_in_matrix")
    meta["mutation_status"] = meta.entry_id.map(ms).fillna("not_in_matrix")
    meta.to_csv(RES / "entry_metadata.csv", index=False)
    cm.to_csv(RES / "chain_metadata.csv", index=False)
    fig = figure(meta, cm)
    n = len(meta)
    save_fig(fig, "Y02_metadata_overview", section="E. Experimental metadata", task="Y02",
             shows="pH, precipitant class, space groups, release years, LBD mutation status, chains per ASU",
             takeaway=(f"{meta.ph.notna().sum()}/{n} entries report a pH (median {meta.ph.median():.1f}); "
                       f"{(meta.additive_class == 'PEG').sum()} PEG-grown; most common space group "
                       f"{meta.space_group.value_counts().index[0]} ({meta.space_group.value_counts().iloc[0]} entries); "
                       f"coactivator/other peptide in {meta.peptide_present.sum()} entries."),
             source="results/ext/entry_metadata.csv, results/ext/chain_metadata.csv", script="scripts/ext_metadata.py")
    write_numbers("Y02", dict(n_entries=n, n_ph=int(meta.ph.notna().sum()), ph_median=float(meta.ph.median()),
                              n_peg=int((meta.additive_class == "PEG").sum()), n_peptide=int(meta.peptide_present.sum()),
                              n_space_groups=int(meta.space_group.nunique()),
                              top_space_group=meta.space_group.value_counts().index[0],
                              n_rcsb_gap_queries=len(gaps),
                              chains_mutation_status=cm.mutation_status.value_counts().to_dict(),
                              n_xray=int((meta.exp_method == "X-RAY DIFFRACTION").sum())))
    print(meta.describe(include="all").T[["count", "unique", "top"]].to_string())
    print(cm.mutation_status.value_counts())


if __name__ == "__main__":
    main()
