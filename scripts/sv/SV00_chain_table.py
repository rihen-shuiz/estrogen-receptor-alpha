"""SV00 step 6-7: shared chain table (one row per kept chain) and the C-alpha observation mask.

Rows: the 788 kept chains in matrix order (results/rmsd_entry_index.csv, tracker S05 out).
Existing columns are reused (sources in results/sv/SV00/column_sources.md, written by this script):
S08 cluster (results/entry_clusters.csv cut_1.5), QC tier (results/qc_tiers.csv, K05), H12 label (results/h12_metric.csv,
S16), entry-level peptide and construct mutations (results/ext/chain_metadata.csv, Y02), ligand codes and status
(results/chain_ligand_of_interest.csv, Q18), primary-ligand heavy atoms and basic-amine flag (results/ext/
structural_features.csv from results/ext/ligand_descriptors.csv, Y09/Y11; the existing SMARTS definition), resolution,
R-free and entry mean B (results/qc_tiers.csv), chain mean C-alpha B (results/lbd_chains.csv, Q04), space group
(results/entries.csv, Q05), ESR1 chains in the ASU and release year (results/ext/entry_metadata.csv, Y02).
SV00 columns: pH, growth / ambient temperature, R-work, TLS flag, multi_valued (SV00_metadata.py, mmCIF + section 3
cleaning), UniProt-author offset (SIFTS, SV00_metadata.py), remediated (SV00_gate.py), chain-level bound peptide and
dimer partner (SV00_contacts.py), latest revision date (data/sv/manifest.csv).
New chemistry columns: carboxylic-acid flag = RDKit match of SMARTS [CX3](=O)[OX1H0-,OX2H1] on the primary ligand's
canonical isomeric SMILES (results/ligands_of_interest.csv, the SMILES the basic-amine flag used); ligand_chemotype =
amine (basic amine only) / acid (acid only) / both / other (ligand with neither) / none (no ligand of interest).
Primary ligand = the chain's first ligand-of-interest row (the existing convention of scripts/ext_common.chain_table).
Mask: results/sv/SV00/observed_mask.npy, bool [788 chains x 246 positions 305-550], True = C-alpha observed in the
regenerated residue mapping (data/sv/lbd_extract/residue_mapping.csv has_ca: occupancy > 0 or unknown), with
observed_mask_rows.csv (row, chain_key, entry_id) and observed_mask_cols.csv (col, esr1_position).
Outputs: results/sv/SV00/chain_table.csv, missing_by_column.csv, column_sources.md, observed_mask*.
Run: .venv/bin/python -I scripts/sv/SV00_chain_table.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
OUT = RES / "sv/SV00"
ACID = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
LBD = np.arange(305, 551)


def rd(path, **kw):
    return pd.read_csv(path, keep_default_na=False, na_values=[""], low_memory=False, **kw)


def has_acid(smiles):
    if not isinstance(smiles, str) or not smiles:
        return np.nan
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    return bool(mol.HasSubstructMatch(ACID))


def truth(value):
    """True / False from a bool or its CSV text; anything else (blank) is False."""
    return value is True or value == "True" or (isinstance(value, (bool, np.bool_)) and bool(value))


def chemotype(amine, acid, has_ligand):
    if not has_ligand:
        return "none"
    amine, acid = truth(amine), truth(acid)
    return "both" if amine and acid else "amine" if amine else "acid" if acid else "other"


def main():
    idx = rd(RES / "rmsd_entry_index.csv")
    t = idx[["chain_key", "entry_id"]].copy()
    lbd_old = rd(RES / "lbd_chains.csv").set_index("chain_key")
    t["auth_asym_id"] = t.chain_key.map(lbd_old.auth_asym_id)

    cl = rd(RES / "entry_clusters.csv").set_index("chain_key")
    t["s08_cluster"] = t.chain_key.map(cl["cut_1.5"]).astype(int)
    qc = rd(RES / "qc_tiers.csv").set_index("entry_id")
    for c in ["qc_label", "tier_resolution", "tier_rfree", "tier_b"]:
        t[c] = t.entry_id.map(qc[c])
    h12 = rd(RES / "h12_metric.csv").set_index("chain_key")
    t["h12_label"] = t.chain_key.map(h12.h12_label)
    t["n_h12_ca"] = t.chain_key.map(h12.n_h12_ca)

    cm = rd(RES / "ext/chain_metadata.csv").set_index("chain_key")
    t["peptide_present_entry"] = t.chain_key.map(cm.peptide_present)
    t["peptide_class_entry"] = t.chain_key.map(cm.peptide_class)
    pep = rd(OUT / "chain_peptide.csv").set_index("chain_key")
    t["peptide_bound"] = t.chain_key.map(pep.peptide_bound)
    t["peptide_bound_source"] = t.chain_key.map(pep.peptide_source)
    t["peptide_bound_auth_asym_id"] = t.chain_key.map(pep.peptide_auth_asym_id)
    t["peptide_bound_length"] = t.chain_key.map(pep.peptide_length)
    t["construct_mutations"] = t.chain_key.map(cm.mutations)
    t["mutation_status"] = t.chain_key.map(cm.mutation_status)

    loi = rd(RES / "chain_ligand_of_interest.csv")
    codes = loi.groupby("chain_key").ligand_id.apply(lambda s: list(dict.fromkeys(x for x in s if isinstance(x, str))))
    t["ligand_codes"] = t.chain_key.map(lambda k: ";".join(codes.get(k, [])))
    t["n_ligand_codes"] = t.chain_key.map(lambda k: len(codes.get(k, [])))
    first = loi.drop_duplicates("chain_key").set_index("chain_key")
    t["ligand_status"] = t.chain_key.map(first.status)
    sf = rd(RES / "ext/structural_features.csv").set_index("chain_key")
    t["ligand_id_primary"] = t.chain_key.map(sf.ligand_id)
    if not (t.ligand_id_primary.fillna("") == t.chain_key.map(first.ligand_id).fillna("")).all():
        raise ValueError("structural_features ligand_id is not the first ligand-of-interest row")
    t["ligand_heavy_atoms"] = t.chain_key.map(sf.ligand_heavy_atoms)
    t["ligand_basic_amine"] = t.chain_key.map(sf.ligand_basic_amine)
    lo = rd(RES / "ligands_of_interest.csv").set_index("ligand_id")
    acid = {lig: has_acid(lo.loc[lig, "canonical_isomeric_smiles"]) for lig in lo.index}
    t["ligand_carboxylic_acid"] = t.ligand_id_primary.map(lambda x: acid.get(x, np.nan) if isinstance(x, str) else np.nan)
    has_lig = t.ligand_id_primary.map(lambda x: isinstance(x, str) and x != "")
    t["ligand_chemotype"] = [chemotype(a, c, h) for a, c, h in zip(t.ligand_basic_amine, t.ligand_carboxylic_acid, has_lig)]

    t["resolution_A"] = t.entry_id.map(qc.resolution_A)
    t["r_free"] = t.entry_id.map(qc.r_free)
    meta = rd(OUT / "entry_metadata_mmcif.csv").set_index("entry_id")
    t["r_work"] = t.entry_id.map(meta.r_work)
    t["mean_b_entry_A2"] = t.entry_id.map(qc.mean_b_A2)
    t["mean_ca_b_chain_A2"] = t.chain_key.map(lbd_old.mean_CA_B_A2)
    t["ph"] = t.entry_id.map(meta.ph)
    t["growth_temp_K"] = t.entry_id.map(meta.growth_temp_K)
    t["ambient_temp_K"] = t.entry_id.map(meta.ambient_temp_K)
    ent = rd(RES / "entries.csv").set_index("entry_id")
    t["space_group"] = t.entry_id.map(ent.space_group)
    em = rd(RES / "ext/entry_metadata.csv").set_index("entry_id")
    t["esr1_chains_asu"] = t.entry_id.map(em.chains_asu)
    t["release_year"] = t.entry_id.map(em.release_year)
    t["tls_refined"] = t.entry_id.map(meta.tls_refined)
    off = rd(OUT / "chain_offsets.csv").set_index("chain_key")
    t["uniprot_author_offset"] = t.chain_key.map(off.uniprot_author_offset)
    t["offset_constant"] = t.chain_key.map(off.offset_sifts_constant)
    t["offset_source"] = t.chain_key.map(off.offset_source)
    t["multi_valued"] = t.entry_id.map(meta.multi_valued)
    t["multi_valued_fields"] = t.entry_id.map(meta.multi_valued_fields)
    man = rd(ROOT / "data/sv/manifest.csv").set_index("entry")
    t["latest_revision_date"] = t.entry_id.map(man.latest_revision_date)
    gate = json.loads((OUT / "gate.json").read_text())
    t["remediated"] = t.chain_key.isin(gate["remediated_chains"])

    dim = rd(OUT / "dimer_partner.csv").set_index("chain_key")
    for c in ["dimer_partner_chain_key", "dimer_partner_auth_asym_id", "dimer_partner_source",
              "dimer_partner_h11_atoms", "dimer_second_best_h11_atoms", "dimer_partner_self_image"]:
        t[c] = t.chain_key.map(dim[c])
    t["dimer_partner_peptide_bound"] = t.dimer_partner_chain_key.map(pep.peptide_bound)
    t["dimer_partner_cluster"] = t.dimer_partner_chain_key.map(cl["cut_1.5"])
    has_partner = t.dimer_partner_source != "none"
    for c in ["dimer_partner_peptide_bound", "dimer_partner_cluster"]:
        t[c] = t[c].astype(object).where(has_partner, np.nan)
    if len(t) != 788 or t.chain_key.duplicated().any():
        raise ValueError(f"chain table has {len(t)} rows")
    t.to_csv(OUT / "chain_table.csv", index=False)

    miss = pd.DataFrame({"column": t.columns, "n_missing": [int(t[c].isna().sum() + (t[c].astype(str) == "").sum()
                                                                if t[c].dtype == object else t[c].isna().sum())
                                                            for c in t.columns]})
    miss["n_rows"] = len(t)
    miss["fraction_missing"] = (miss.n_missing / len(t)).round(4)
    miss.to_csv(OUT / "missing_by_column.csv", index=False)

    # observation mask (regenerated mapping) and its check against the frozen mapping
    cols = ["chain_key", "esr1_position", "has_ca"]
    mask = np.zeros((len(t), len(LBD)), dtype=bool)
    frozen = np.zeros_like(mask)
    row = {k: i for i, k in enumerate(t.chain_key)}
    for path, arr in ((ROOT / "data/sv/lbd_extract/residue_mapping.csv", mask),
                      (RES / "lbd_extract/residue_mapping.csv", frozen)):
        m = rd(path, usecols=cols)
        m = m[m.chain_key.isin(row) & m.has_ca.astype(bool) & m.esr1_position.between(LBD[0], LBD[-1])]
        arr[m.chain_key.map(row).to_numpy(), (m.esr1_position.astype(int) - LBD[0]).to_numpy()] = True
    np.save(OUT / "observed_mask.npy", mask)
    pd.DataFrame({"row": range(len(t)), "chain_key": t.chain_key, "entry_id": t.entry_id}).to_csv(
        OUT / "observed_mask_rows.csv", index=False)
    pd.DataFrame({"col": range(len(LBD)), "esr1_position": LBD}).to_csv(OUT / "observed_mask_cols.csv", index=False)
    n_rows_diff = int((mask != frozen).any(axis=1).sum())
    print(f"chain table {t.shape}; mask {mask.shape}, {int(mask.sum())} observed CA; rows differing from frozen "
          f"mapping: {n_rows_diff}")
    print("chemotype:", t.ligand_chemotype.value_counts().to_dict())
    print("partner source:", t.dimer_partner_source.value_counts().to_dict())
    print(miss[miss.n_missing > 0].to_string(index=False))
    (OUT / "mask_check.json").write_text(json.dumps({"observed_ca": int(mask.sum()), "rows_differ_from_frozen_mapping":
                                                     n_rows_diff, "shape": list(mask.shape)}, indent=1) + "\n")
    write_sources(t)
    crosscheck(t)


def crosscheck(t):
    """Entry-level agreement of the chain-table values with the other existing / mmCIF values of the same quantity."""
    e = t.drop_duplicates("entry_id").set_index("entry_id")
    em = rd(RES / "ext/entry_metadata.csv").set_index("entry_id")
    meta = rd(OUT / "entry_metadata_mmcif.csv").set_index("entry_id")
    lbd = rd(ROOT / "data/sv/lbd_chains.csv")
    n_esr1 = lbd[(lbd.identity == "esr1") & (lbd.model.astype(str) == "1")].groupby("entry_id").size()
    pairs = [("ph", e.ph, em.ph, "Y02 ph (also parsed from crystal-growth text)"),
             ("growth_temp_K", e.growth_temp_K, em.temp_K, "Y02 temp_K"),
             ("ambient_temp_K", e.ambient_temp_K, em.diffraction_temp_K, "Y02 diffraction_temp_K"),
             ("r_work", e.r_work, em.r_work, "Y02 r_work"),
             ("resolution_A", e.resolution_A, meta.resolution_mmcif, "mmCIF _refine.ls_d_res_high"),
             ("r_free", e.r_free, meta.r_free_mmcif, "mmCIF _refine.ls_R_factor_R_free"),
             ("esr1_chains_asu", e.esr1_chains_asu, n_esr1, "regenerated Q04 table: ESR1 chains in model 1"),
             ("release_year", e.release_year, meta.initial_release_date.astype(str).str[:4].astype(float),
              "mmCIF initial release (revision ordinal 1)"),
             ("space_group", e.space_group, meta.space_group_mmcif, "mmCIF _symmetry.space_group_name_H-M")]
    rows = []
    for name, a, b, other in pairs:
        b = b.reindex(e.index)
        if name == "space_group":
            a2, b2 = a.astype(str).str.replace(" ", ""), b.astype(str).str.replace(" ", "")
            both = a.notna() & b.notna()
            agree = (a2[both] == b2[both]).sum()
        else:
            a2, b2 = pd.to_numeric(a, errors="coerce"), pd.to_numeric(b, errors="coerce")
            both = a2.notna() & b2.notna()
            agree = (np.abs(a2[both] - b2[both]) < 1e-3).sum()
            a, b = a2, b2
        rows.append({"column": name, "compared_with": other, "n_entries": len(e), "n_both": int(both.sum()),
                     "n_agree": int(agree), "n_only_chain_table": int((a.notna() & b.isna()).sum()),
                     "n_only_other": int((a.isna() & b.notna()).sum())})
    pd.DataFrame(rows).to_csv(OUT / "value_crosscheck.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def write_sources(t):
    src = {
        "chain_key, entry_id": "results/rmsd_entry_index.csv (S05), matrix row order",
        "auth_asym_id": "results/lbd_chains.csv (Q04)",
        "s08_cluster": "results/entry_clusters.csv cut_1.5 (S08, 1.5 A average linkage)",
        "qc_label, tier_resolution, tier_rfree, tier_b": "results/qc_tiers.csv (K05)",
        "h12_label, n_h12_ca": "results/h12_metric.csv (S16)",
        "peptide_present_entry, peptide_class_entry": "results/ext/chain_metadata.csv (Y02; entry level: any non-receptor "
                                                      "polypeptide of <= peptide_max_length = 50 residues in the entry)",
        "peptide_bound, peptide_bound_source, peptide_bound_auth_asym_id, peptide_bound_length":
            "SV00 results/sv/SV00/chain_peptide.csv (SV00_contacts.py): peptide copy (same definition, <= 50 residues) "
            "with the most receptor heavy atoms within 4.0 A, ASU or symmetry copy",
        "construct_mutations, mutation_status": "results/ext/chain_metadata.csv mutations, mutation_status (Y02, from "
                                                "the residue-level mapping of results/lbd_chains.csv)",
        "ligand_codes, n_ligand_codes, ligand_status": "results/chain_ligand_of_interest.csv (Q18); codes = distinct "
                                                       "ligand_id per chain; status of the first row",
        "ligand_id_primary": "first ligand-of-interest row per chain (= results/ext/structural_features.csv ligand_id, "
                             "the convention of scripts/ext_common.chain_table)",
        "ligand_heavy_atoms, ligand_basic_amine": "results/ext/structural_features.csv (Y11) from "
            "results/ext/ligand_descriptors.csv (Y09): RDKit heavy atoms; basic amine = SMARTS "
            "[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N-a);!$(N-C=N);!$(N#*);!$(N=*)]([C;!$(C=*)]) (existing definition)",
        "ligand_carboxylic_acid": "SV00: RDKit SMARTS [CX3](=O)[OX1H0-,OX2H1] on results/ligands_of_interest.csv "
                                  "canonical_isomeric_smiles of ligand_id_primary (acid or carboxylate only)",
        "ligand_chemotype": "SV00: amine / acid / both / other / none from ligand_basic_amine and ligand_carboxylic_acid; "
                            "none = no ligand of interest",
        "resolution_A, r_free, mean_b_entry_A2": "results/qc_tiers.csv (K05; values of results/entries.csv, Q05)",
        "r_work": "SV00 results/sv/SV00/entry_metadata_mmcif.csv: mmCIF _refine.ls_R_factor_R_work (first row); "
                  "Data API only where the mmCIF has none (r_work_source)",
        "mean_ca_b_chain_A2": "results/lbd_chains.csv mean_CA_B_A2 (Q04, chain mean C-alpha B)",
        "ph, growth_temp_K": "SV00 entry_metadata_mmcif.csv: _exptl_crystal_grow.pH / .temp, crystal_id 1, '?'/'.' -> "
                             "NaN, [DEFAULT] plausible pH 3-11, 270-315 K, else NaN; Data API fill only where missing",
        "ambient_temp_K": "SV00 entry_metadata_mmcif.csv: _diffrn.ambient_temp, diffrn id 1, [DEFAULT] 50-320 K",
        "space_group": "results/entries.csv space_group (Q05)",
        "esr1_chains_asu, release_year": "results/ext/entry_metadata.csv chains_asu, release_year (Y02)",
        "tls_refined": "SV00 entry_metadata_mmcif.csv: category _pdbx_refine_tls exists in the mmCIF",
        "uniprot_author_offset, offset_constant, offset_source": "SV00 results/sv/SV00/chain_offsets.csv: SIFTS (PDBe "
            "API mappings/uniprot, P03372 segments, UniProt 305-550) minus author number (_pdbx_poly_seq_scheme); "
            "offset = most common per-residue value, offset_constant = one value over the LBD",
        "multi_valued, multi_valued_fields": "SV00 entry_metadata_mmcif.csv: several rows in _exptl_crystal_grow, _diffrn "
                                             "or _refine (id 1 / first row used)",
        "latest_revision_date": "data/sv/manifest.csv (_pdbx_audit_revision_history of the restored file)",
        "remediated": "SV00 results/sv/SV00/gate.json remediated_chains (only filled when the gate is 'tolerated')",
        "dimer_partner_*": "SV00 results/sv/SV00/dimer_partner.csv (SV00_contacts.py): receptor copy with the most H11 "
                           "(497-530) heavy atoms within 4.0 A; source asu / symmetry / none; self_image = the chain's "
                           "own symmetry image",
        "dimer_partner_peptide_bound": "peptide_bound of the partner chain (chain_peptide.csv)",
        "dimer_partner_cluster": "results/entry_clusters.csv cut_1.5 of the partner chain (blank if the partner is not "
                                 "a kept chain)",
    }
    lines = ["# SV00 chain_table.csv: column sources", "",
             f"{len(t)} rows (kept chains, matrix order), {len(t.columns)} columns. Missing values per column: "
             "results/sv/SV00/missing_by_column.csv.", "", "| column(s) | source |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in src.items()]
    (OUT / "column_sources.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
