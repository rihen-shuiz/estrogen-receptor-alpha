"""Which verified chain is a defensible coordinate reference? Require 95% LBD,
all configured core/H12 positions and known passing historical QC, then prefer
best resolution. Save both completeness-first and resolution-first shortlists,
the chosen context and its actual deposited primary-citation identifiers.
"""
from .run_io import read_table, save_table, write_json


def run_reference(root,out,settings):
    chains = read_table(out,"receptor_chains.csv")
    entries = read_table(out,"entries.csv")
    candidates = chains[(chains.identity=="esr1") & chains.eligible_historical].merge(
        entries[["pdb_id","resolution_A","r_free","mean_B_factor_A2","citation_title","citation_doi","citation_pmid","citation_journal"]],on="pdb_id",validate="many_to_one")
    contacts = read_table(out,"contacts.csv")
    ligands = read_table(out,"ligand_instances.csv")
    ligand_lookup = ligands.set_index("ligand_instance_key")
    ligand_context, partner_context = {}, {}
    for key in candidates.chain_key:
        rows = contacts[contacts.chain_key==key]
        instances = rows.loc[rows.contact_type=="receptor_ligand","ligand_instance_key"].dropna().unique()
        codes = []
        for instance in instances:
            ligand = ligand_lookup.loc[instance]
            if ligand.category=="candidate_small_molecule":
                codes.append(ligand.ligand_code)
        ligand_context[key] = ";".join(sorted(set(codes))) or "no candidate ASU contacts observed"
        partners = set()
        if "partner_chain_key" in contacts:
            partner_rows = contacts[(contacts.chain_key==key) | (contacts.partner_chain_key==key)]
            partners = set(partner_rows.loc[partner_rows.contact_type.isin(["receptor_receptor","receptor_partner"]),"partner_chain_key"].dropna())
            partners |= set(partner_rows.loc[partner_rows.partner_chain_key==key,"chain_key"].dropna())
        partners.discard(key)
        instances = set()
        if "partner_instance_id" in contacts:
            polymer_rows = contacts[(contacts.chain_key==key) & contacts.contact_type.isin(["receptor_receptor","receptor_partner"])]
            instances = set(polymer_rows.partner_instance_id.dropna())
            reverse_rows = contacts[(contacts.partner_chain_key==key) & contacts.contact_type.isin(["receptor_receptor","receptor_partner"])]
            instances |= set(reverse_rows.instance_id.dropna())
        partner_context[key] = ";".join(sorted(instances or partners)) or "no polymer contacts observed"
    candidates["contacting_candidate_ligand_codes"] = candidates.chain_key.map(ligand_context)
    candidates["contacting_polymer_chain_keys"] = candidates.chain_key.map(partner_context)
    candidates["high_completeness"] = ((candidates.lbd_fraction>=settings["reference_min_lbd_fraction"]) &
        (candidates.core_fraction>=settings["reference_min_core_fraction"]) &
        (candidates.h12_fraction>=settings["reference_min_h12_fraction"]) & (candidates.quality=="pass"))
    candidates = candidates.sort_values(["mapped_lbd_ca","mapped_core_ca","mapped_h12_ca","resolution_A","chain_key"],ascending=[False,False,False,True,True])
    candidates["completeness_rank"] = range(1,len(candidates)+1)
    by_resolution = candidates.sort_values(["resolution_A","mapped_lbd_ca","mismatch_count","chain_key"],ascending=[True,False,True,True])
    ranks = {key:rank for rank,key in enumerate(by_resolution.chain_key,1)}
    candidates["resolution_rank"] = candidates.chain_key.map(ranks)
    save_table(out,"reference_candidates.csv",candidates)
    eligible = by_resolution[by_resolution.high_completeness]
    if eligible.empty:
        raise RuntimeError("No reference satisfies explicit high-completeness + known-QC policy. Inspect reference_candidates.csv; do not silently relax.")
    chosen = select_reference_chain(candidates, settings)
    record = {key:chosen[key] for key in ["chain_key","pdb_id","label_asym_id","auth_asym_id","model"]}
    record["model"] = str(record["model"])
    record["policy"] = (f"known historical QC pass; LBD >={settings['reference_min_lbd_fraction']:.0%}, "
        f"core >={settings['reference_min_core_fraction']:.0%}, H12 >={settings['reference_min_h12_fraction']:.0%}; "
        "best resolution, then coverage, fewer differences, chain key")
    record["historical_reference"] = settings["historical_reference"]
    write_json(out/"reference.json",record)
    selected_contacts = contacts[contacts.chain_key==chosen.chain_key]
    save_table(out,"reference_contacts.csv",selected_contacts)
    text = ["REFERENCE DOSSIER — computed from frozen mmCIF data",f"Chosen: {chosen.chain_key}; author chain {chosen.auth_asym_id}",record["policy"],
        f"LBD {chosen.mapped_lbd_ca}/{chosen.lbd_expected}, core {chosen.mapped_core_ca}/{chosen.core_expected}, H12 {chosen.mapped_h12_ca}/{chosen.h12_expected} CA atoms.",
        f"Resolution {chosen.resolution_A} Å; depositor R-free {chosen.r_free}; entry mean B {chosen.mean_B_factor_A2} Å²; chain CA mean B {chosen.mean_CA_B_A2} Å².",
        f"Declared differences: {chosen.mismatches}",f"Primary citation: {chosen.citation_title}; {chosen.citation_journal}; DOI {chosen.citation_doi}; PMID {chosen.citation_pmid}",
        f"Contacting candidate ligand codes (ASU): {chosen.contacting_candidate_ligand_codes}; polymer source-chain contacts: {chosen.contacting_polymer_chain_keys}.",
        "Contacts and site confidence: reference_contacts.csv. Assembly evidence: assemblies.csv.",
        "Publication reading, molecular-viewer and validation inspection are recorded separately in reference_inspection.txt.",
        "3Q95|B|1 remains the historical reference; new matrices will be recomputed with this chosen chain.",
        "Resolution alone does not establish stability, pharmacological action or biological superiority."]
    (out/"REFERENCE_DOSSIER.txt").write_text("\n".join(text)+"\n")
    print(f"Reference: {chosen.chain_key}, {chosen.resolution_A} Å; {len(eligible)} high-completeness candidates. Dossier: {out}/REFERENCE_DOSSIER.txt")
    return record

def select_reference_chain(chains, settings, *, receptor="esr1", policy=None):
    """Shared completeness/QC/ranking policy; explicit subtype-specific limits.

    Chains include entry resolution metadata. ESR1 defaults retain the existing
    stage policy. Other receptors require an explicit policy rather than silently
    relaxing the ESR1 reference thresholds.
    """
    if receptor != "esr1" and policy is None:
        raise ValueError("A subtype-specific reference policy must be explicit")
    policy = policy or {"lbd_fraction": settings['reference_min_lbd_fraction'],
                        "core_fraction": settings['reference_min_core_fraction'],
                        "h12_fraction": settings['reference_min_h12_fraction']}
    candidates = chains.loc[
        chains.identity.eq(receptor) & chains.eligible_historical & chains.quality.eq('pass')
    ].copy()
    for field, limit in policy.items():
        candidates = candidates.loc[candidates[field].ge(limit)]
    candidates = candidates.sort_values(
        ['resolution_A', 'mapped_lbd_ca', 'mismatch_count', 'chain_key'],
        ascending=[True, False, True, True],
    )
    candidates = candidates.loc[candidates.resolution_A.notna()]
    if candidates.empty:
        raise ValueError('No reference satisfies the declared subtype completeness and known-QC policy')
    return candidates.iloc[0]


def main():
    from ._task_io import parser, prepare, json_write
    args = parser('Block 08: original known-QC completeness and resolution ranking').parse_args()
    root, source, out, settings = prepare(args)
    chains = read_table(source, 'receptor_chains.csv')
    entries = read_table(source, 'entries.csv')
    candidates = chains.merge(entries[['pdb_id', 'resolution_A', 'r_free']], on='pdb_id', validate='many_to_one')
    chosen = select_reference_chain(candidates, settings)
    candidates.to_csv(out / 'reference_candidates.csv', index=False)
    record = chosen.to_dict()
    record['model'] = str(record['model'])
    record['policy'] = 'known QC pass; configured LBD/core/H12 completeness; resolution, coverage, differences, key'
    json_write(out / 'reference.json', record)
    print(f"Block 08: {chosen.chain_key}")


if __name__ == "__main__":
    main()
