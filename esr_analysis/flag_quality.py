"""How do coverage and crystallographic quality affect inclusion? Join one
entry's depositor/fallback quality to its independently mapped chains. Compare
the historical selector with all receptor instances; preserve missing QC as
unknown. Outputs explain exclusions instead of assigning biology from a count.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .map_residues import residue_letter, missing_stretches
from .read_structure import read_cif
from .run_io import read_table, save_table


def quality_status(resolution, rfree, mean_b, settings):
    """Existing combined cutoff status; fail takes precedence over missing values."""
    values = [resolution, rfree, mean_b]
    limits = [settings["qc_resolution_A"], settings["qc_r_free"], settings["qc_mean_B_A2"]]
    status = "pass"
    if any(np.isfinite(v) and v > limit for v, limit in zip(values, limits)):
        status = "fail"
    elif any(not np.isfinite(v) for v in values):
        status = "unknown"
    return status


def normalize_cached_coordinates(root, out, settings):
    """Repair cached CA presence and missing residue IDs before fitting.

    Existing occupancy identifies unusable CA atoms. Only entries with affected
    rows are read again for residue metadata and chain observation totals; atom
    correspondence, assembly/contact distances and thresholds are unchanged.
    Repeated quality runs leave the repaired coordinate tables unchanged.
    """
    from .extract_lbd import residue_observations, coordinate_metadata

    chains = read_table(out,"receptor_chains.csv")
    mapping = read_table(out,"residue_mapping.csv")
    entries = read_table(out,"entries.csv")
    invalid_ca = mapping.has_ca & mapping.occupancy.notna() & (mapping.occupancy <= 0)
    missing_ids = mapping.has_residue & ~mapping.has_ca & mapping.auth_seq_id.isna()
    affected = mapping[invalid_ca | missing_ids]
    if affected.empty:
        return
    for field in ["auth_seq_id","insertion_code","altloc","observed_aa","residue_status"]:
        mapping[field] = mapping[field].astype(object)
    for field in ["missing_lbd_stretches","missing_core_stretches","missing_h12_stretches",
                  "eligibility_reason","mismatches"]:
        chains[field] = chains[field].astype(object)
    chain_index = chains.set_index("chain_key")
    source_paths = entries.set_index("pdb_id").source_path
    first_models = {}
    for pdb_id, entry_rows in affected.groupby("pdb_id",sort=False):
        path = Path(source_paths.loc[pdb_id])
        if not path.is_absolute():
            path = Path(root)/path
        data, atoms = read_cif(path)
        first_models[pdb_id] = str(atoms.pdbx_PDB_model_num.iloc[0])
        for key, rows in entry_rows.groupby("chain_key",sort=False):
            chain = chain_index.loc[key]
            model = key.split("|")[2]
            observed = atoms[(atoms.label_asym_id==chain.label_asym_id) & (atoms.pdbx_PDB_model_num==model)]
            ca, residues, raw_ca = residue_observations(observed)
            for index, row in rows.iterrows():
                query = int(row.label_seq_id)
                residue, metadata = coordinate_metadata(query,ca,residues,raw_ca)
                if residue is None:
                    raise ValueError(f"Cached observed residue absent from raw mmCIF: {key}, label position {query}")
                for field, value in metadata.items():
                    mapping.at[index,field] = value
                actual_aa = residue_letter(residue.label_comp_id,data)
                mapping.at[index,"observed_aa"] = actual_aa if metadata["has_residue"] else ""
                difference = row.difference_annotation if pd.notna(row.difference_annotation) else ""
                mapping.at[index,"residue_status"] = ("match" if actual_aa==row.expected_aa else
                    "deposited_difference" if difference else "unexplained_difference")

            chain_row = chains.index[chains.chain_key==key][0]
            chains.at[chain_row,"observed_ca"] = len(ca)
            chains.at[chain_row,"observed_residues"] = len(residues)
            chains.at[chain_row,"mean_CA_B_A2"] = pd.Series([atom.B_iso_or_equiv for atom in ca.values()],dtype=float).mean()
            block = mapping[mapping.chain_key==key]
            mapped_ca = set(block.loc[block.has_ca,"canonical_position"].astype(int))
            mapped_residues = set(block.loc[block.has_residue,"canonical_position"].astype(int))
            declared = set(block.canonical_position.astype(int))
            equivalent = json.loads(chain.canonical_to_esr1)
            regions = {
                "lbd": set(range(settings["esr1_lbd_start"],settings["esr1_lbd_end"]+1)),
                "core": {p for start,end in settings["esr1_core_ranges"] for p in range(start,end+1)},
                "h12": set(range(settings["h12_range"][0],settings["h12_range"][1]+1)),
            }
            for region, esr1_positions in regions.items():
                expected = esr1_positions if chain.identity=="esr1" else {
                    int(native) for native,e1 in equivalent.items() if e1 in esr1_positions}
                count = len(mapped_ca & expected)
                chains.at[chain_row,f"mapped_{region}_ca"] = count
                chains.at[chain_row,f"{region}_fraction"] = count/len(expected) if expected else float("nan")
                chains.at[chain_row,f"missing_{region}_stretches"] = missing_stretches(expected,mapped_ca)
                if region=="lbd":
                    chains.at[chain_row,"observed_lbd_residues"] = len(mapped_residues & expected)
                    declared_lbd = declared & expected
            eligible = (chains.at[chain_row,"mapped_lbd_ca"] >= settings["min_chain_lbd_ca"] and
                        chains.at[chain_row,"mapped_core_ca"] >= settings["min_core_ca_for_alignment"])
            chains.at[chain_row,"eligible_historical"] = eligible
            chains.at[chain_row,"eligibility_reason"] = ("eligible_historical_coverage" if eligible else
                "non_LBD_fragment" if not declared_lbd else "short_construct_LBD" if
                len(declared_lbd)<settings["min_chain_lbd_ca"] else "incomplete_observed_LBD")
            differences = block[block.residue_status!="match"]
            chains.at[chain_row,"mismatch_count"] = len(differences)
            chains.at[chain_row,"mismatches"] = ";".join(
                f"{int(row.canonical_position)}{row.expected_aa}>{row.observed_aa if pd.notna(row.observed_aa) and row.observed_aa!='' else row.declared_aa}:{row.difference_annotation if pd.notna(row.difference_annotation) else 'unresolved cause'}"
                for row in differences.itertuples())
    for index, entry in entries.iterrows():
        if entry.pdb_id not in first_models:
            continue
        first = chains[(chains.pdb_id==entry.pdb_id) & (chains.model.astype(str).str.replace(r"\.0$","",regex=True)==first_models[entry.pdb_id])]
        entries.at[index,"eligible_receptor_chains"] = int(first.eligible_historical.sum())
    save_table(out,"residue_mapping.csv",mapping)
    save_table(out,"receptor_chains.csv",chains)
    save_table(out,"entries.csv",entries)
    print(f"Coordinate normalization: excluded {int(invalid_ca.sum())} zero/negative-occupancy CA atoms; restored {int(missing_ids.sum())} observed-residue IDs from {affected.pdb_id.nunique()} raw entries.")


def validate_inventory_labels(out,settings):
    """Check cached declared correspondence and synchronize chemical categories.

    Geometry is unchanged: this joins already observed contact instances to
    their canonical category. Low sequence agreement halts cached-input reuse;
    it is never repaired by silently relabeling a protein as a receptor.
    """
    chains = read_table(out,"receptor_chains.csv")
    mapping = read_table(out,"residue_mapping.csv")
    mapping["sequence_match"] = mapping.declared_aa==mapping.expected_aa
    agreement = mapping.groupby("chain_key").sequence_match.mean()
    chains["declared_sequence_agreement"] = chains.chain_key.map(agreement)
    receptor = chains.identity.isin(["esr1","esr2"])
    failed = receptor & (chains.declared_sequence_agreement.isna() | (chains.declared_sequence_agreement<settings["annotated_min_identity"]))
    if failed.any():
        raise ValueError("Cached receptor accession contradicts its deposited sequence: " + ", ".join(chains.loc[failed,"chain_key"]) + ". Inspect raw UniProt segments and rerun inventory into a new directory.")
    chains["identity_sequence_checked"] = receptor
    save_table(out,"receptor_chains.csv",chains)
    contacts = read_table(out,"contacts.csv")
    ligands = read_table(out,"ligand_instances.csv")
    lookup = ligands.set_index("ligand_instance_key").category
    ligand_contacts = contacts.contact_type=="receptor_ligand"
    corrected = contacts.loc[ligand_contacts,"ligand_instance_key"].map(lookup)
    if corrected.isna().any():
        raise ValueError("A ligand contact has no unique inventory instance")
    contacts.loc[ligand_contacts,"ligand_category"] = corrected
    excluded = ligand_contacts & (contacts.ligand_category!="candidate_small_molecule")
    contacts.loc[excluded,"site"] = "excluded_chemical_category"
    save_table(out,"contacts.csv",contacts)
    print(f"Identity check: {int(receptor.sum())} receptor chains; minimum declared sequence agreement {chains.loc[receptor,'declared_sequence_agreement'].min():.3f}. Categories joined to ligand contacts.")


def run_quality(root,out,settings):
    normalize_cached_coordinates(root,out,settings)
    validate_inventory_labels(out,settings)
    chains = read_table(out,"receptor_chains.csv")
    entries = read_table(out,"entries.csv")
    counts = []
    comparisons = []
    for receptor in ("esr1","esr2"):
        historical = read_table(root/"data", "receptor_chain_choices.csv")
        historical = historical[historical.receptor == receptor].rename(columns={
            "author_chain":"chain_id", "historical_eligible":"usable_for_structural_analysis",
            "historical_lbd_ca":"mapped_lbd_ca"})
        active_entries = set(entries.loc[entries.query_memberships.str.contains(receptor),"pdb_id"])
        historical = historical[historical.pdb_id.isin(active_entries)]
        for row in historical.itertuples():
            observed = chains[(chains.pdb_id==row.pdb_id)&(chains.auth_asym_id==row.chain_id)]
            receptor_rows = chains[(chains.pdb_id==row.pdb_id)&(chains.identity==receptor)]
            reason = "selected_chain_verified" if (observed.identity==receptor).any() else "wrong_or_unresolved_selected_chain"
            if reason=="selected_chain_verified" and not receptor_rows.eligible_historical.any():
                reason = ";".join(sorted(set(receptor_rows.eligibility_reason)))
            comparisons.append({"pdb_id":row.pdb_id,"receptor":receptor,"historical_author_chain":row.chain_id,
                                "historical_eligible":row.usable_for_structural_analysis,"historical_lbd_ca":row.mapped_lbd_ca,
                                "corrected_eligible_chains":int(receptor_rows.eligible_historical.sum()),"reason":reason})
        current = chains[chains.identity==receptor]
        eligible = current[current.eligible_historical]
        counts.append({"receptor":receptor,"query_entries":len(active_entries),"historical_eligible_entries":int(historical.usable_for_structural_analysis.sum()),
                       "verified_chains":len(current),"eligible_chains":len(eligible),"eligible_entries":eligible.pdb_id.nunique(),
                       "eligible_qc_pass_entries":eligible.loc[eligible.quality=="pass","pdb_id"].nunique(),"qc_unknown_chains":int((eligible.quality=="unknown").sum())})
    save_table(out,"filtering_counts.csv",counts)
    save_table(out,"historical_selection_comparison.csv",comparisons)
    print(f"Quality: coverage and QC remain separate. Counts and exclusions: {out}/filtering_counts.csv")


def main():
    from ._task_io import parser, prepare
    args = parser('Block 05: recompute metadata QC only; do not refit or recluster').parse_args()
    root, source, out, settings = prepare(args)
    entries = read_table(source, 'entries.csv')
    entries['qc_status'] = [quality_status(row.resolution_A, row.r_free, row.mean_B_factor_A2, settings) for row in entries.itertuples()]
    save_table(out, 'entries.csv', entries)
    chains = read_table(source, 'receptor_chains.csv')
    chains['quality'] = chains.pdb_id.map(entries.set_index('pdb_id').qc_status)
    save_table(out, 'receptor_chains.csv', chains)
    print('Block 05: ' + str(entries.qc_status.value_counts().to_dict()))


if __name__ == "__main__":
    main()
