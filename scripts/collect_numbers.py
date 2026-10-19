"""Collect the numbers for the team presentation: results/numbers.json and docs/numbers.md.

Reads current results/ and data/ files and the config. A number read from results/task_outputs/
is marked (Sep). docs/numbers.md has one line per number: "name: value (source file)".
Run from anywhere: python scripts/collect_numbers.py
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import ligand_of_interest as loi  # noqa: E402  (pocket residues and contact cutoff of task Q18)

SETTINGS = "esr_analysis/resources/config/analysis_settings.json"
NUMBERS = []


def read(path):
    return pd.read_csv(ROOT / path, dtype=str, keep_default_na=False)


def ids(path):
    return [line.strip() for line in (ROOT / path).read_text().splitlines() if line.strip()]


def add(name, value, *sources, sep=False):
    if any(item["name"] == name for item in NUMBERS):
        raise ValueError(f"duplicate number name {name}")
    NUMBERS.append({"name": name, "value": value, "source": list(sources), "sep": sep})


def text(value):
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, list):
        return ",".join(str(item) for item in value)
    return str(value)


def dataset(settings):
    esr1, esr2 = ids("data/pdb_ids_esr1.txt"), ids("data/pdb_ids_esr2.txt")
    entries, lbd = read("results/entries.csv"), read("results/lbd_chains.csv")
    extract, manifest = read("results/lbd_extract/entries.csv"), read("data/raw/download_manifest.csv")
    add("query_esr1_entries", len(esr1), "data/pdb_ids_esr1.txt")
    add("query_esr2_entries", len(esr2), "data/pdb_ids_esr2.txt")
    add("entries_master_table", len(entries), "results/entries.csv")
    add("raw_cif_download_ok", int((manifest.download_ok == "True").sum()), "data/raw/download_manifest.csv")
    add("lbd_extract_status_ok", int((extract.status == "ok").sum()), "results/lbd_extract/entries.csv")
    add("esr1_entries_with_esr1_chain", lbd.loc[lbd.identity == "esr1", "entry_id"].nunique(),
        "results/lbd_chains.csv")
    no_receptor = extract[extract.pdb_id.isin(esr1) & (extract.declared_receptor_chains == "0")]
    add("esr1_entries_without_receptor_chain", len(no_receptor), "results/lbd_extract/entries.csv")
    add("esr1_chain_models", int((lbd.identity == "esr1").sum()), "results/lbd_chains.csv")
    eligible = lbd[(lbd.identity == "esr1") & (lbd.eligible_historical == "True")]
    add("esr1_eligible_chains", len(eligible), "results/lbd_chains.csv")
    add("esr1_selected_entries", eligible.entry_id.nunique(), "results/lbd_chains.csv")
    add("lbd_chain_rows", len(lbd), "results/lbd_chains.csv")
    add("eligibility_min_lbd_ca", settings["min_chain_lbd_ca"], SETTINGS)
    add("eligibility_min_core_ca", settings["min_core_ca_for_alignment"], SETTINGS)


def qc(settings):
    entries, lbd = read("results/entries.csv"), read("results/lbd_chains.csv")
    selected = set(lbd.loc[(lbd.identity == "esr1") & (lbd.eligible_historical == "True"), "entry_id"])
    sets = [("all", entries, ["results/entries.csv"]),
            ("selected", entries[entries.entry_id.isin(selected)], ["results/entries.csv", "results/lbd_chains.csv"])]
    flags = [("qc_res", "resolution"), ("qc_rfree", "r_free"), ("qc_b", "mean_B"), ("qc_pass", "combined")]
    for set_name, frame, sources in sets:
        for column, label in flags:
            counts = frame[column].value_counts()
            for value in ["pass", "fail", "unknown"]:
                add(f"qc_{set_name}_{label}_{value}", int(counts.get(value, 0)), *sources)
    add("qc_cutoff_resolution_A", settings["qc_resolution_A"], SETTINGS)
    add("qc_cutoff_r_free", settings["qc_r_free"], SETTINGS)
    add("qc_cutoff_mean_B_A2", settings["qc_mean_B_A2"], SETTINGS)


def reference(settings):
    entries, lbd = read("results/entries.csv"), read("results/lbd_chains.csv")
    key = re.search(r"Chain key `([^`]+)`", (ROOT / "docs/reference.md").read_text()).group(1)
    row = lbd[lbd.chain_key == key].iloc[0]
    qc_row = entries[entries.entry_id == key.split("|")[0]].iloc[0]
    add("reference_chain_key", key, "docs/reference.md")
    add("reference_resolution_A", float(qc_row.resolution), "results/entries.csv")
    add("reference_r_free", float(qc_row.r_free), "results/entries.csv")
    add("reference_mean_B_A2", float(qc_row.mean_b), "results/entries.csv")
    lbd_source = "results/lbd_chains.csv"
    add("reference_lbd_ca_mapped", int(row.mapped_lbd_ca), lbd_source)
    add("reference_lbd_ca_expected", int(row.lbd_expected), lbd_source)
    add("reference_core_ca_mapped", int(row.mapped_core_ca), lbd_source)
    add("reference_core_ca_expected", int(row.core_expected), lbd_source)
    add("reference_h12_ca_mapped", int(row.mapped_h12_ca), lbd_source)
    add("reference_h12_ca_expected", int(row.h12_expected), lbd_source)
    add("reference_declared_differences", int(row.mismatch_count), lbd_source)

    # Selection step 2 of docs/reference.md, applied to the eligible ESR1 chains that pass entry QC.
    eligible = lbd[(lbd.identity == "esr1") & (lbd.eligible_historical == "True")]
    eligible = eligible.merge(entries[["entry_id", "qc_pass", "resolution"]], on="entry_id")
    passing = eligible[eligible.qc_pass == "pass"].copy()
    lbd_ca, lbd_expected = pd.to_numeric(passing.mapped_lbd_ca), pd.to_numeric(passing.lbd_expected)
    core_ca, core_expected = pd.to_numeric(passing.mapped_core_ca), pd.to_numeric(passing.core_expected)
    h12_ca, h12_expected = pd.to_numeric(passing.mapped_h12_ca), pd.to_numeric(passing.h12_expected)
    keep = (lbd_ca / lbd_expected >= 0.95) & (core_ca == core_expected) & (h12_ca == h12_expected)
    step2 = passing[keep].copy()
    step2["res"] = step2.resolution.astype(float)
    step2["lbd_ca"] = lbd_ca[keep]
    step2["diffs"] = step2.mismatch_count.astype(int)
    step2 = step2.sort_values(["res", "lbd_ca", "diffs", "chain_key"], ascending=[True, False, True, True])
    add("reference_step2_chains", len(step2), "results/lbd_chains.csv", "results/entries.csv")
    add("reference_step2_entries", step2.entry_id.nunique(), "results/lbd_chains.csv", "results/entries.csv")
    add("reference_step2_first_chain", step2.chain_key.iloc[0], "results/lbd_chains.csv", "results/entries.csv")
    add("reference_step2_runner_up_chain", step2.chain_key.iloc[1], "results/lbd_chains.csv", "results/entries.csv")
    add("reference_step2_runner_up_resolution_A", float(step2.res.iloc[1]), "results/entries.csv")

    sep_reference = json.loads((ROOT / "results/task_outputs/block_08/reference.json").read_text())
    add("reference_chain_key_sep", sep_reference["chain_key"], "results/task_outputs/block_08/reference.json",
        sep=True)


def ligands_of_interest():
    chains, table = read("results/chain_ligand_of_interest.csv"), read("results/ligands_of_interest.csv")
    with_ligand = chains[chains.ligand_id != ""]
    source = "results/chain_ligand_of_interest.csv"
    add("loi_chains_evaluated", chains.chain_key.nunique(), source)
    add("loi_chain_ligand_rows", len(chains), source)
    add("loi_chains_with_ligand", with_ligand.chain_key.nunique(), source)
    add("loi_chains_apo_or_unclear", int((chains.status == "apo_or_unclear").sum()), source)
    add("loi_ligand_pairs", len(with_ligand), source)
    add("loi_multi_ligand_rows", int((chains.status == "multi").sum()), source)
    add("loi_source_subject_of_investigation", int((chains.source == "subject_of_investigation").sum()), source)
    add("loi_source_pocket_contact", int((chains.source == "pocket_contact").sum()), source)
    add("loi_entries_with_ligand", with_ligand.entry_id.nunique(), source)
    add("loi_distinct_ligand_ids", with_ligand.ligand_id.nunique(), source)
    add("loi_ligands_in_clustering", int((table.in_clustering == "True").sum()), "results/ligands_of_interest.csv")
    add("loi_pocket_residues", len(loi.POCKET_POSITIONS), "scripts/ligand_of_interest.py")
    add("loi_contact_cutoff_A", loi.CONTACT_CUTOFF_A, "scripts/ligand_of_interest.py")
    add("loi_min_pocket_contacts", loi.MIN_POCKET_CONTACTS, "scripts/ligand_of_interest.py")


def chemistry():
    sets = [("480", "results/ligand_clusters.csv", "results/kmeans_scores.csv"),
            ("462", "results/ligand_clusters_noadd.csv", "results/kmeans_scores_noadd.csv"),
            ("306", "results/ligand_clusters_loi.csv", "results/kmeans_scores_loi.csv")]
    for tag, clusters_path, scores_path in sets:
        scores, clusters = read(scores_path), read(clusters_path)
        add(f"chem_{tag}_ligands", len(clusters), clusters_path)
        silhouette = {int(k): float(value) for k, value in zip(scores.k, scores.silhouette)}
        for k in sorted(silhouette):
            add(f"chem_{tag}_silhouette_k{k}", silhouette[k], scores_path)
        best = max(silhouette, key=silhouette.get)
        add(f"chem_{tag}_best_k", best, scores_path)
        add(f"chem_{tag}_best_silhouette", silhouette[best], scores_path)
        sizes = clusters[f"cluster_k{best}"].astype(int).value_counts().sort_index()
        add(f"chem_{tag}_cluster_sizes_best_k", [int(size) for size in sizes], clusters_path)
    settings = (ROOT / "docs/chemistry_k_choice.md").read_text()
    add("chem_kmeans_n_init", int(re.search(r"n_init=(\d+)", settings).group(1)), "docs/chemistry_k_choice.md")
    add("chem_kmeans_seed", int(re.search(r"seed (\d+)", settings).group(1)), "docs/chemistry_k_choice.md")


def chembl():
    activities = read("data/chembl_esr1_activities.csv")
    source = "data/chembl_esr1_activities.csv"
    add("chembl_activity_rows", len(activities), source)
    for standard_type in ["IC50", "Ki", "Kd"]:
        add(f"chembl_standard_type_{standard_type}", int((activities.standard_type == standard_type).sum()), source)
    add("chembl_pchembl_present", int((activities.pchembl_value != "").sum()), source)
    add("chembl_molecules_distinct", activities.molecule_chembl_id.nunique(), source)
    blank = activities.loc[activities.canonical_smiles == "", "molecule_chembl_id"]
    add("chembl_molecules_blank_canonical_smiles", blank.nunique(), source)
    matches = read("results/chembl_pdb_matches.csv")
    source = "results/chembl_pdb_matches.csv"
    add("chembl_pdb_match_rows", len(matches), source)
    add("chembl_pdb_match_ligand_codes", matches.ligand_id.nunique(), source)
    add("chembl_pdb_match_molecules", matches.molecule_chembl_id.nunique(), source)
    add("chembl_pdb_match_exact", int((matches.match_type == "exact").sum()), source)
    add("chembl_pdb_match_tanimoto_ge_0_9", int((matches.match_type == "tanimoto>=0.9").sum()), source)
    add("chembl_pdb_match_min_tanimoto", float(matches.tanimoto.astype(float).min()), source)


def esr2():
    entries = read("results/esr2_entries.csv")
    source = "results/esr2_entries.csv"
    add("esr2_entries", len(entries), source)
    add("esr2_chain_rows", int(entries.esr2_chains.astype(int).sum()), source)
    add("esr2_eligible_chain_rows", int(entries.esr2_eligible_chains.astype(int).sum()), source)
    for column, label in [("qc_res", "resolution"), ("qc_rfree", "r_free"), ("qc_b", "mean_B"), ("qc_pass", "combined")]:
        counts = entries[column].value_counts()
        for value in ["pass", "fail"]:
            add(f"esr2_qc_{label}_{value}", int(counts.get(value, 0)), source)
    numbering = read("results/numbering_map_esr1_esr2.csv")
    source = "results/numbering_map_esr1_esr2.csv"
    add("esr2_numbering_positions", len(numbering), source)
    add("esr2_numbering_same_residue", int((numbering.same_residue == "True").sum()), source)
    add("esr2_numbering_in_lbd", int((numbering.in_lbd == "True").sum()), source)
    pairs = read("results/matched_pairs.csv")
    source = "results/matched_pairs.csv"
    add("esr2_matched_pairs", len(pairs), source)
    add("esr2_matched_esr1_entries", pairs.esr1_entry.nunique(), source)
    add("esr2_matched_esr2_entries", pairs.esr2_entry.nunique(), source)
    add("esr2_matched_identical_codes", int((pairs.match_type == "identical").sum()), source)
    add("esr2_matched_tanimoto_1", int((pairs.tanimoto.astype(float) == 1.0).sum()), source)


def labels():
    current = "results/ligand_action_labels.csv"
    add("labels_current_file", "present" if (ROOT / current).exists() else "absent", current)
    identities = read("results/task_outputs/block_06/ligand_identities.csv")
    actions = read("results/task_outputs/block_07/ligand_actions.csv")
    ident_source = "results/task_outputs/block_06/ligand_identities.csv"
    action_source = "results/task_outputs/block_07/ligand_actions.csv"
    current_codes = set(read("results/ligands.csv").ligand_id)
    sep_codes = set(identities.ligand_code)
    add("labels_sep_identities", len(sep_codes), ident_source, sep=True)
    add("labels_sep_codes_in_current_table", len(sep_codes & current_codes), ident_source,
        "results/ligands.csv", sep=True)
    add("labels_current_codes_not_in_sep", len(current_codes - sep_codes), "results/ligands.csv",
        ident_source, sep=True)
    add("labels_sep_supported", int((actions.evidence_supported == "True").sum()), action_source, sep=True)
    add("labels_sep_unsupported", int((actions.evidence_supported == "False").sum()), action_source, sep=True)
    supported = actions.loc[actions.supported_action_type != "", "supported_action_type"]
    for action_type, count in supported.value_counts().sort_index().items():
        add(f"labels_sep_supported_{action_type}", int(count), action_source, sep=True)
    provisional = actions[actions.action_status == "provisional"]
    add("labels_sep_provisional_rows", len(provisional), action_source, sep=True)
    for action_type, count in provisional.action_type.value_counts().sort_index().items():
        add(f"labels_sep_provisional_{action_type}", int(count), action_source, sep=True)


def sorted_sizes(values):
    return [int(size) for size in pd.Series(values).value_counts().tolist()]


def provisional_tests():
    """P01-P04 from results/provisional/. The September 1.5 Å groups are read from results/task_outputs/."""
    prov = "results/provisional"
    sensitivity = f"{prov}/entropy_sensitivity.csv"
    for row in read(sensitivity).itertuples(index=False):
        add(f"p01_n_clusters_cut_{row.cut}A", int(row.n_clusters), sensitivity)
        add(f"p01_top3_share_cut_{row.cut}A", float(row.top3_fraction), sensitivity)
        add(f"p01_entropy_bits_cut_{row.cut}A", float(row.entropy_bits), sensitivity)
    clusters_path = f"{prov}/clusters.csv"
    clusters = read(clusters_path)
    sep_path = "results/task_outputs/block_11/entry_clusters.csv"
    sep = read(sep_path)
    ours = clusters["cut_1.5"].astype(int)
    add("p01_sizes_cut_1.5A", sorted_sizes(ours), clusters_path)
    add("p01_sep_sizes_cut_1.5A", sorted_sizes(sep["cluster_1.5A"].astype(int)), sep_path, sep=True)
    sep_label = dict(zip(sep.chain_key, sep["cluster_1.5A"].astype(int)))
    agreement = adjusted_rand_score([sep_label[key] for key in clusters.chain_key], ours)
    add("p01_matches_sep_cut_1.5A", "yes" if agreement == 1.0 else "no", clusters_path, sep_path, sep=True)

    saturation_path = f"{prov}/saturation.csv"
    for row in read(saturation_path).itertuples(index=False):
        add(f"p02_n90_cut_{row.cut}A", int(row.n_90), saturation_path)

    delta_path = f"{prov}/qc_delta.csv"
    for row in read(delta_path).itertuples(index=False):
        tag = row.set.replace("_", "")
        if row.cut == "1.0":
            add(f"p03_{tag}_n_entries", int(row.n_entries), delta_path)
        add(f"p03_{tag}_n_clusters_cut_{row.cut}A", int(row.n_clusters), delta_path)
        add(f"p03_{tag}_top3_share_cut_{row.cut}A", float(row.top3_fraction), delta_path)
        add(f"p03_{tag}_entropy_bits_cut_{row.cut}A", float(row.entropy_bits), delta_path)
        add(f"p03_{tag}_n90_cut_{row.cut}A", int(row.n_90), delta_path)
        if row.set != "all_351":
            add(f"p03_clusters_absent_count_cut_{row.cut}A", int(row.n_clusters_absent), delta_path)
            add(f"p03_clusters_absent_cut_{row.cut}A", row.clusters_absent or "none", delta_path)

    test_path = f"{prov}/chem_struct_test.csv"
    test = read(test_path)
    chi = test[test.test == "chi_squared_independence"].iloc[0]
    perm = test[test.test == "permutation_ligand_level"].iloc[0]
    add("p04_n", int(chi.n), test_path)
    add("p04_chi2_statistic", float(chi.statistic), test_path)
    add("p04_chi2_df", int(chi.df), test_path)
    add("p04_chi2_p_value", float(chi.p_value), test_path)
    add("p04_cells_expected_below_5", int(chi.cells_expected_below_5), test_path)
    add("p04_cells_total", int(chi.n_cells), test_path)
    add("p04_permutation_p_value", float(perm.p_value), test_path)
    add("p04_permutations", int(perm.n_permutations), test_path)
    add("p04_ligands", int(chi.n_ligands), test_path)
    add("p04_excluded_no_ligand", int(chi.n_excluded_no_ligand), test_path)
    add("p04_multi_ligand_resolved", int(chi.n_multi_resolved), test_path)


def write():
    settings = json.loads((ROOT / SETTINGS).read_text())
    dataset(settings)
    qc(settings)
    reference(settings)
    ligands_of_interest()
    chemistry()
    chembl()
    esr2()
    labels()
    provisional_tests()
    (ROOT / "results/numbers.json").write_text(json.dumps({"numbers": NUMBERS}, indent=2) + "\n")
    lines = [f"{item['name']}: {text(item['value'])} ({', '.join(item['source'])})" + (" (Sep)" if item["sep"] else "")
             for item in NUMBERS]
    (ROOT / "docs/numbers.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(NUMBERS)} numbers ({sum(item['sep'] for item in NUMBERS)} marked Sep) to "
          "results/numbers.json and docs/numbers.md")


if __name__ == "__main__":
    write()
