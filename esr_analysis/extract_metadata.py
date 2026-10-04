"""Block 03: deposited entry metadata with calculated mean-B fallback.
Input: parsed mmCIF, selected atoms, protein labels, settings. Output: metadata
with source fields and existing combined QC status. Missing values stay unknown."""
import numpy as np
from .read_structure import category_rows, first_value
from .flag_quality import quality_status


def numeric_field(data, names):
    value, source = first_value(data, names)
    return (float(value) if value is not None else np.nan), source


def entry_quality(pdb_id, data, atoms, polymer_labels, settings):
    resolution, resolution_source = numeric_field(data, ["_refine.ls_d_res_high", "_reflns.d_resolution_high"])
    rfree, rfree_source = numeric_field(data, ["_refine.ls_R_factor_R_free"])
    mean_b, b_source = numeric_field(data, ["_refine.B_iso_mean"])
    model = atoms["pdbx_PDB_model_num"].iloc[0]
    if np.isnan(mean_b):
        protein_atoms = atoms[atoms["label_asym_id"].isin(polymer_labels) & (atoms["pdbx_PDB_model_num"] == model)]
        mean_b = protein_atoms["B_iso_or_equiv"].mean()
        b_source = "computed first-model protein atoms; occupancy-selected alternates"
    status = quality_status(resolution, rfree, mean_b, settings)
    citation = category_rows(data, "_citation")
    primary = next((r for r in citation if r.get("id") == "primary"), citation[0] if citation else {})
    return {"pdb_id": pdb_id, "resolution_A": resolution, "resolution_source": resolution_source,
            "r_free": rfree, "r_free_source": rfree_source + "; depositor", "mean_B_factor_A2": mean_b,
            "B_source": b_source, "qc_status": status,
            "citation_title": primary.get("title", "unknown"), "citation_doi": primary.get("pdbx_database_id_DOI", ""),
            "citation_pmid": primary.get("pdbx_database_id_PubMed", ""), "citation_journal": primary.get("journal_abbrev", "unknown"),
            "space_group": first_value(data, ["_space_group.name_H-M_alt", "_symmetry.space_group_name_H-M"])[0],
            "deposit_date": first_value(data, ["_pdbx_database_status.recvd_initial_deposition_date"])[0]}



def main():
    import pandas as pd
    from ._task_io import parser, prepare
    from .read_structure import read_cif
    from .run_io import sha256
    p = parser('Block 03: deposited metadata and declared molecular contents of local mmCIFs')
    p.add_argument('--ids', nargs='+', default=['3Q95', '7NFB', '7BAA', '6CHZ'])
    p.add_argument('--raw-folder', type=__import__('pathlib').Path)
    args = p.parse_args()
    root, source, out, settings = prepare(args)
    rows = []
    for pdb in args.ids:
        path = (args.raw_folder or root / 'data/raw/esr1') / (pdb + '.cif')
        data, atoms = read_cif(path)
        polymers = {row['entity_id']: row for row in category_rows(data, '_entity_poly')}
        labels = [row['id'] for row in category_rows(data, '_struct_asym') if row['entity_id'] in polymers and 'polypeptide' in polymers[row['entity_id']]['type']]
        record = entry_quality(pdb, data, atoms, labels, settings)
        record.update(declared_protein_chains=len(labels), model_count=atoms.pdbx_PDB_model_num.nunique(), source_path=str(path), source_sha256=sha256(path))
        rows.append(record)
    pd.DataFrame(rows).to_csv(out / 'entries.csv', index=False)
    print(f'Block 03: metadata for {len(rows)} raw entries')


if __name__ == "__main__":
    main()
