"""Block 04: canonical receptor mapping, domain applicability and observations.

Input: parsed mmCIF/atoms, frozen sequences, subtype equivalence and settings.
Output: entry composition, all declared chain/model rows and residue mapping.
Retain fragments, missing coordinates, occupancy and mutation evidence.
No master-entry filtering; eligibility remains the existing 180/60 screen."""
import numpy as np
import pandas as pd
from .read_structure import category_rows
from .map_residues import (identify_entity, map_entity, residue_letter, missing_stretches, receptor_accessions)
from .extract_metadata import entry_quality
from .run_io import positions


def residue_observations(atoms):
    """Separate usable CA atoms, observed residues, and deposited CA occupancy.

    Zero/negative occupancy coordinates are not observations. Unknown occupancy
    remains usable, matching the contact atom policy. Residue metadata comes
    from another observed atom when CA is absent, without inventing a CA position.
    """
    usable = atoms[atoms.occupancy.isna() | (atoms.occupancy > 0)]
    residues = {}
    for atom in usable.itertuples():
        if str(atom.label_seq_id).isdigit():
            residues.setdefault(int(atom.label_seq_id), atom)
    raw_ca = {int(atom.label_seq_id): atom for atom in atoms[atoms.label_atom_id=="CA"].itertuples()
              if str(atom.label_seq_id).isdigit()}
    ca = {position: atom for position, atom in raw_ca.items()
          if pd.isna(atom.occupancy) or atom.occupancy > 0}
    return ca, residues, raw_ca


def coordinate_metadata(query, ca_lookup, residue_lookup, raw_ca_lookup):
    """Return residue identifiers independently of an observed CA coordinate."""
    ca = ca_lookup.get(query)
    raw_ca = raw_ca_lookup.get(query)
    residue = ca if ca is not None else residue_lookup.get(query, raw_ca)
    return residue, {
        "has_residue": query in residue_lookup, "has_ca": ca is not None,
        "auth_seq_id": residue.auth_seq_id if residue is not None else "",
        "insertion_code": residue.pdbx_PDB_ins_code if residue is not None else "",
        "altloc": ca.label_alt_id if ca is not None else "",
        # Preserve deposited zero occupancy as evidence for excluding its CA.
        "occupancy": float(raw_ca.occupancy) if raw_ca is not None else np.nan,
        "B_factor_A2": float(ca.B_iso_or_equiv) if ca is not None else np.nan,
        "x_A": float(ca.Cartn_x) if ca is not None else np.nan,
        "y_A": float(ca.Cartn_y) if ca is not None else np.nan,
        "z_A": float(ca.Cartn_z) if ca is not None else np.nan,
    }


def inventory_entry(pdb_id, data, atoms, sequences, equivalence, settings):
    entities = {r["id"]: r for r in category_rows(data, "_entity")}
    polymers = {r["entity_id"]: r for r in category_rows(data, "_entity_poly")}
    asym = [r for r in category_rows(data, "_struct_asym") if r["entity_id"] in polymers]
    protein_labels = [r["id"] for r in asym if "polypeptide" in polymers[r["entity_id"]]["type"]]
    entry = entry_quality(pdb_id, data, atoms, protein_labels, settings)
    esr1_lbd = set(range(settings["esr1_lbd_start"], settings["esr1_lbd_end"]+1))
    core = set(positions(settings["esr1_core_ranges"]))
    h12 = set(range(settings["h12_range"][0], settings["h12_range"][1]+1))
    dbd = set(range(settings["dna_binding_range"][0], settings["dna_binding_range"][1]+1))
    chain_rows, mapping_rows = [], []
    for instance in asym:
        label, entity = instance["id"], instance["entity_id"]
        polymer = polymers[entity]
        sequence = "".join(polymer["pdbx_seq_one_letter_code_can"].split())
        identity, identity_source, accessions = identify_entity(entity, sequence, data, sequences, settings, polymer["type"])
        chain_atoms = atoms[atoms["label_asym_id"] == label]
        models = chain_atoms["pdbx_PDB_model_num"].unique().tolist() or [atoms["pdbx_PDB_model_num"].iloc[0]]
        pairs, mapping_source = [], "not a verified receptor"
        if identity in ("esr1", "esr2"):
            pairs, mapping_source = map_entity(entity, sequence, data, sequences[identity], receptor_accessions(identity,settings))
        native_to_esr1 = equivalence if identity == "esr2" else {c:c for q,c in pairs}
        expected_lbd = {native for native, e1 in native_to_esr1.items() if e1 in esr1_lbd}
        # ESR1 denominator is the configured canonical domain, even for a short fragment.
        if identity == "esr1":
            expected_lbd = esr1_lbd
        expected_core = {native for native,e1 in native_to_esr1.items() if e1 in core}
        expected_h12 = {native for native,e1 in native_to_esr1.items() if e1 in h12}
        if identity == "esr1":
            expected_core, expected_h12 = core, h12
        for model in models:
            observed = chain_atoms[chain_atoms["pdbx_PDB_model_num"] == model]
            auth = observed["auth_asym_id"].iloc[0] if not observed.empty else "?"
            if observed.empty:
                scheme = [r for r in category_rows(data, "_pdbx_poly_seq_scheme") if r.get("asym_id") == label]
                if scheme:
                    auth = scheme[0].get("pdb_strand_id", "?")
            chain_key = f"{pdb_id}|{label}|{model}"
            ca_lookup, residue_lookup, raw_ca_lookup = residue_observations(observed)
            residue_observed = set(residue_lookup)
            mapped_ca, mapped_residues, mismatches = set(), set(), []
            local_rows = []
            differences = category_rows(data, "_struct_ref_seq_dif")
            for query, canonical in pairs:
                residue, atom_metadata = coordinate_metadata(query, ca_lookup, residue_lookup, raw_ca_lookup)
                e1 = native_to_esr1.get(canonical)
                expected_aa = sequences[identity][canonical-1]
                actual_aa = residue_letter(residue.label_comp_id, data) if residue is not None else sequence[query-1]
                mismatch = actual_aa != expected_aa
                deposited_diff = [r for r in differences if r.get("seq_num") == str(query) and r.get("pdbx_pdb_strand_id", auth) == auth]
                detail = ";".join(r.get("details", "") for r in deposited_diff)
                status = "match"
                if mismatch:
                    status = "deposited_difference" if detail else "unexplained_difference"
                    mismatches.append(f"{canonical}{expected_aa}>{actual_aa}:{detail or 'unresolved cause'}")
                if atom_metadata["has_ca"]:
                    mapped_ca.add(canonical)
                if query in residue_observed:
                    mapped_residues.add(canonical)
                row = {"chain_key": chain_key,"pdb_id":pdb_id,"identity":identity,"label_seq_id":query,
                       "canonical_position":canonical,"esr1_position":e1,"declared_aa":sequence[query-1],
                       "expected_aa":expected_aa,"observed_aa":actual_aa if query in residue_observed else "",
                       "residue_status":status,"difference_annotation":detail,
                       **atom_metadata}
                local_rows.append(row)
            mapping_rows.extend(local_rows)
            declared = {c for q,c in pairs}
            esr1_declared = {native_to_esr1[c] for c in declared if c in native_to_esr1}
            lbd_count = len(mapped_ca & expected_lbd)
            core_count = len(mapped_ca & expected_core)
            domain = "not_receptor"
            if identity in ("esr1", "esr2"):
                in_lbd, in_dbd = bool(esr1_declared & esr1_lbd), bool(esr1_declared & dbd)
                domain = "multiple_domains" if in_lbd and in_dbd else "LBD" if in_lbd else "DNA_binding" if in_dbd else "other_receptor_fragment"
            eligible = identity in ("esr1", "esr2") and lbd_count >= settings["min_chain_lbd_ca"] and core_count >= settings["min_core_ca_for_alignment"]
            reason = "eligible_historical_coverage"
            if not eligible:
                reason = "not_receptor" if identity not in ("esr1", "esr2") else "non_LBD_fragment" if not declared & expected_lbd else "short_construct_LBD" if len(declared & expected_lbd) < settings["min_chain_lbd_ca"] else "incomplete_observed_LBD"
            chain_rows.append({"chain_key":chain_key,"pdb_id":pdb_id,"entity_id":entity,"label_asym_id":label,
                "auth_asym_id":auth,"model":str(model),"description":entities[entity].get("pdbx_description",""),
                "identity":identity,"identity_source":identity_source,"accessions":";".join(accessions),
                "polymer_type":polymer["type"],"sequence_length":len(sequence),"observed_residues":len(residue_observed),
                "observed_ca":len(ca_lookup),"domain":domain,"domain_source":settings["domain_annotation_source"],
                "mapping_source":mapping_source,"construct_lbd_positions":len(declared & expected_lbd),
                "observed_lbd_residues":len(mapped_residues & expected_lbd),"mapped_lbd_ca":lbd_count,
                "mapped_core_ca":core_count,"mapped_h12_ca":len(mapped_ca & expected_h12),
                "lbd_expected":len(expected_lbd),"core_expected":len(expected_core),"h12_expected":len(expected_h12),
                "lbd_fraction":lbd_count/len(expected_lbd) if expected_lbd else np.nan,
                "core_fraction":core_count/len(expected_core) if expected_core else np.nan,
                "h12_fraction":len(mapped_ca & expected_h12)/len(expected_h12) if expected_h12 else np.nan,
                "missing_lbd_stretches":missing_stretches(expected_lbd,mapped_ca),
                "missing_core_stretches":missing_stretches(expected_core,mapped_ca),
                "missing_h12_stretches":missing_stretches(expected_h12,mapped_ca),
                "mismatch_count":len(mismatches),"mismatches":";".join(mismatches),
                "eligible_historical":eligible,"eligibility_reason":reason,"quality":entry["qc_status"],
                "mean_CA_B_A2":pd.Series([ca.B_iso_or_equiv for ca in ca_lookup.values()],dtype=float).mean(),
                "label_to_canonical":dict(pairs),"canonical_to_esr1":native_to_esr1})
    first_model = str(atoms["pdbx_PDB_model_num"].iloc[0])
    first_chains = [r for r in chain_rows if r["model"] == first_model]
    receptors = [r for r in first_chains if r["identity"] in ("esr1","esr2")]
    proteins = [r for r in first_chains if "polypeptide" in r["polymer_type"]]
    entry.update({"declared_protein_chains":len(proteins),"observed_protein_chains":sum(r["observed_residues"]>0 for r in proteins),
                  "declared_receptor_chains":len(receptors),"observed_receptor_chains":sum(r["observed_residues"]>0 for r in receptors),
                  "eligible_receptor_chains":sum(r["eligible_historical"] for r in receptors),
                  "dna_chains":sum(r["identity"]=="dna" for r in first_chains),"rna_chains":sum(r["identity"]=="rna" for r in first_chains),
                  "partner_chains":sum(r["identity"] in ("partner","peptide") for r in first_chains),
                  "model_count":atoms["pdbx_PDB_model_num"].nunique(),"status":"ok"})
    return entry, chain_rows, mapping_rows



def main():
    from ._task_io import parser, prepare, sequences
    from .read_structure import read_cif
    from .map_residues import receptor_equivalence
    from .run_io import save_table, sha256
    p = parser('Block 04: identify and map all declared chains/models in explicit raw entries')
    p.add_argument('--ids', nargs='+', default=['3Q95', '7NFB', '7BAA', '6CHZ'])
    p.add_argument('--raw-folder', type=__import__('pathlib').Path)
    args = p.parse_args()
    root, source, out, settings = prepare(args)
    seq = sequences(root, settings)
    equivalence = receptor_equivalence(seq)
    entries, chains, mapping = [], [], []
    for pdb in args.ids:
        path = (args.raw_folder or root / 'data/raw/esr1') / (pdb + '.cif')
        data, atoms = read_cif(path)
        entry, local_chains, local_mapping = inventory_entry(pdb, data, atoms, seq, equivalence, settings)
        entry.update(source_path=str(path), source_sha256=sha256(path), query_memberships='esr1')
        entries.append(entry)
        chains.extend(local_chains)
        mapping.extend(local_mapping)
    save_table(out, 'entries.csv', entries)
    save_table(out, 'receptor_chains.csv', chains)
    save_table(out, 'residue_mapping.csv', mapping)
    print(f'Block 04: {len(chains)} declared chain/models; fragments retained')


if __name__ == "__main__":
    main()
