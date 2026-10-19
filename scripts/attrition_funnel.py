"""Tracker K01: attrition funnel. Every step that removes or ignores a chain, entry or residue.

Reads files only (no recomputation of coordinates). Inputs:
  data/pdb_ids_esr1.txt, data/pdb_ids_esr2.txt        entry origin (ESR1 query vs ESR2 validation)
  results/lbd_extract/receptor_chains.csv             all declared chain-models of the 526 entries
  results/lbd_chains.csv                              ESR1 chains with eligibility and mapped CA counts (Q04)
  results/lbd_extract/residue_mapping.csv             CA coordinates and B-factors per canonical position
  esr_analysis/resources/data/receptor_chain_choices.csv   historical author chain per entry (one per entry)
  esr_analysis/resources/frozen/rmsd_entry_index.txt  old strict matrix index (351 chains)
  esr_analysis/resources/config/analysis_settings.json    thresholds (min_chain_lbd_ca, min_pairwise_common_ca, ...)
  results/chain_ligand_of_interest.csv                apo labels (Q18)
  results/task_outputs/block_13/helix_distribution_and_coverage.csv   old regional-B coverage cohorts
Outputs: results/attrition_funnel.csv (one row per step) and results/chain_loss_reasons.csv (one row per ESR1 chain).
Nothing is deleted. Chains kept out of a matrix keep their row and a reason.
Run from the repo root: python scripts/attrition_funnel.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
MIN_LBD = SETTINGS["min_chain_lbd_ca"]
MIN_CORE = SETTINGS["min_core_ca_for_alignment"]
MIN_COMMON = SETTINGS["min_pairwise_common_ca"]
LBD_START, LBD_END = SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"]
R_FREE_LIMIT = SETTINGS["qc_r_free"]


def read(path, **kwargs):
    return pd.read_csv(ROOT / path, keep_default_na=False, na_values=[""], low_memory=False, **kwargs)


def truthy(series):
    return series.astype(str).str.lower().eq("true")


def main():
    esr1_ids = {line.strip().upper() for line in (ROOT / "data/pdb_ids_esr1.txt").read_text().splitlines() if line.strip()}
    esr2_ids = {line.strip().upper() for line in (ROOT / "data/pdb_ids_esr2.txt").read_text().splitlines() if line.strip()}
    chains_all = read("results/lbd_extract/receptor_chains.csv", dtype={"pdb_id": str})
    chains_all["pdb_id"] = chains_all.pdb_id.str.upper()
    esr1_chains = read("results/lbd_chains.csv", dtype={"entry_id": str})
    esr1_chains = esr1_chains[esr1_chains.identity == "esr1"].copy()
    esr1_chains["entry_id"] = esr1_chains.entry_id.str.upper()
    if set(esr1_chains.chain_key) != set(chains_all.loc[chains_all.identity == "esr1", "chain_key"]):
        raise ValueError("ESR1 chain sets differ between lbd_chains.csv and receptor_chains.csv")

    entries = read("results/entries.csv", dtype=str)
    entries["entry_id"] = entries.entry_id.str.upper()
    selected_mask = esr1_chains.eligible_historical.pipe(truthy)
    eligible = esr1_chains[selected_mask].copy()

    # One chain per entry: the historical author chain (receptor_chain_choices.csv, historical_eligible True).
    choices = read("esr_analysis/resources/data/receptor_chain_choices.csv", dtype=str)
    picks = choices[(choices.receptor == "esr1") & truthy(choices.historical_eligible)]
    pick_keys = {}
    for row in picks.itertuples():
        match = eligible[(eligible.entry_id == row.pdb_id.upper()) & (eligible.auth_asym_id.astype(str) == row.author_chain)]
        if len(match) != 1:
            raise ValueError(f"Historical pick {row.pdb_id} {row.author_chain} does not match exactly one eligible chain")
        pick_keys[row.pdb_id.upper()] = match.chain_key.iloc[0]
    frozen = [line.strip() for line in (ROOT / "esr_analysis/resources/frozen/rmsd_entry_index.txt").read_text().splitlines() if line.strip()]
    if set(frozen) != set(pick_keys.values()):
        raise ValueError("Frozen index is not the historical one-chain-per-entry choice")
    esr1_chains["historical_pick"] = esr1_chains.chain_key.isin(set(pick_keys.values()))

    # Residue level: CA coverage and B-factor validity for the eligible chains.
    residues = read("results/lbd_extract/residue_mapping.csv", dtype={"chain_key": str},
                    usecols=["chain_key", "esr1_position", "has_ca", "B_factor_A2"])
    residues = residues[residues.chain_key.isin(set(eligible.chain_key)) & truthy(residues.has_ca)].dropna(subset=["esr1_position"])
    residues["position"] = residues.esr1_position.astype(float).astype(int)
    residues = residues[residues.position.between(LBD_START, LBD_END)]
    positions = LBD_END - LBD_START + 1
    keys = sorted(eligible.chain_key)
    index = {key: i for i, key in enumerate(keys)}
    observed = np.zeros((len(keys), positions), dtype=np.int32)
    for key, pos in zip(residues.chain_key, residues.position):
        observed[index[key], pos - LBD_START] = 1
    overlap = observed @ observed.T
    off_diagonal = overlap[np.triu_indices(len(keys), k=1)]
    pairs_below = int((off_diagonal < MIN_COMMON).sum())
    min_overlap = int(off_diagonal.min())
    invalid_b = int((~np.isfinite(residues.B_factor_A2.astype(float)) | (residues.B_factor_A2.astype(float) < 0)).sum())

    # Old block-13 coverage cohorts (full cohort, coverage threshold 0.8).
    coverage = read("results/task_outputs/block_13/helix_distribution_and_coverage.csv", dtype=str)
    full = coverage[(coverage.cohort == "full") & coverage.coverage_threshold.astype(float).eq(0.8)]
    b_dropped = int(full.cohort_chains.astype(int).sum() - full.retained_chains.astype(int).sum())

    qc_fail_selected = int((entries.set_index("entry_id").loc[list(pick_keys), "r_free"].astype(float) > R_FREE_LIMIT).sum())
    apo = read("results/chain_ligand_of_interest.csv", dtype=str)
    apo_chains = int(apo.loc[apo.status == "apo_or_unclear", "chain_key"].nunique())

    # Completeness-screen reasons, as written by extract_lbd.py.
    ineligible = esr1_chains[~esr1_chains.eligible_historical.pipe(truthy)].copy()

    def completeness_reason(row):
        if row.eligibility_reason == "non_LBD_fragment":
            return "no declared LBD positions in this construct: no LBD overlap, nothing to measure"
        if row.eligibility_reason == "short_construct_LBD":
            return f"declared LBD positions {row.construct_lbd_positions} below min_chain_lbd_ca {MIN_LBD}"
        return (f"mapped LBD CA {row.mapped_lbd_ca} (min {MIN_LBD}) or core CA {row.mapped_core_ca} (min {MIN_CORE}); "
                f"declared LBD positions {row.construct_lbd_positions}")

    ineligible["detail"] = ineligible.apply(completeness_reason, axis=1)
    reasons = {row.chain_key: row.detail for row in ineligible.itertuples()}
    n_non_lbd = int((ineligible.eligibility_reason == "non_LBD_fragment").sum())
    n_threshold = len(ineligible) - n_non_lbd

    in_esr1_entries = chains_all[chains_all.pdb_id.isin(esr1_ids)]
    counts = {identity: int(n) for identity, n in in_esr1_entries.identity.value_counts().items()}
    n_all = len(chains_all)
    n_esr2_entry_chains = int(chains_all.pdb_id.isin(esr2_ids).sum())
    n_after_entry_split = n_all - n_esr2_entry_chains
    n_peptide, n_partner, n_dna = counts.get("peptide", 0), counts.get("partner", 0), counts.get("dna", 0)
    if n_after_entry_split != sum(counts.values()) or counts.get("esr2", 0) != 0:
        raise ValueError("Unexpected identity mix in ESR1 entries")
    lbd_expected_total = int(eligible.lbd_expected.astype(int).sum())
    lbd_observed_total = int(eligible.mapped_lbd_ca.astype(int).sum())
    entries_no_esr1 = len(esr1_ids - set(esr1_chains.entry_id))
    if pairs_below != 0:
        raise ValueError(f"{pairs_below} pairs below min_pairwise_common_ca; run S04 pruning before this funnel")
    old_pruned = len(json.loads((ROOT / "results/archive_v2_strict/alignment/distance_exclusions.json").read_text()))

    steps = [
        ("chain-models of ESR2 validation entries", n_all, n_after_entry_split, n_esr2_entry_chains,
         "chains of the 39 ESR2 entries (data/pdb_ids_esr2.txt) are not ESR1 receptor chains", "esr_analysis/extract_lbd.py; data/pdb_ids_esr2.txt", "yes"),
        ("non-receptor peptide chains in ESR1 entries", n_after_entry_split, n_after_entry_split - n_peptide, n_peptide,
         "identity peptide: not the receptor protein", "esr_analysis/extract_lbd.py identify_entity", "yes"),
        ("non-receptor partner chains in ESR1 entries", n_after_entry_split - n_peptide, n_after_entry_split - n_peptide - n_partner, n_partner,
         "identity partner: not the receptor protein", "esr_analysis/extract_lbd.py identify_entity", "yes"),
        ("DNA chains in ESR1 entries", n_after_entry_split - n_peptide - n_partner, counts.get("esr1", 0), n_dna,
         "identity dna: not the receptor protein", "esr_analysis/extract_lbd.py identify_entity", "yes"),
        ("ESR1 entries without an ESR1 receptor chain", len(esr1_ids), len(esr1_ids) - entries_no_esr1, entries_no_esr1,
         "entry has no ESR1 protein chain (entry-level step)", "esr_analysis/extract_lbd.py", "yes"),
        ("completeness: non-LBD fragment", len(esr1_chains), len(esr1_chains) - n_non_lbd, n_non_lbd,
         "eligibility_reason non_LBD_fragment: no LBD overlap (mathematics)", "esr_analysis/extract_lbd.py; esr_analysis/flag_quality.py", "yes"),
        ("completeness: LBD/core CA below screen", len(esr1_chains) - n_non_lbd, len(eligible), n_threshold,
         f"screen unchanged: mapped LBD CA >= {MIN_LBD} and core CA >= {MIN_CORE}", "esr_analysis/extract_lbd.py; esr_analysis/flag_quality.py", "no"),
        ("one chain per entry (historical author chain)", len(eligible), len(pick_keys), len(eligible) - len(pick_keys),
         "project choice: historical author chain per entry (receptor_chain_choices.csv); equals frozen index",
         "esr_analysis/resources/data/receptor_chain_choices.csv; esr_analysis/resources/frozen/rmsd_entry_index.txt", "no"),
        ("distance pruning, old 351-chain matrix", len(pick_keys), len(pick_keys) - old_pruned, old_pruned,
         f"min_pairwise_common_ca {MIN_COMMON}; chains dropped in the strict run (results/archive_v2_strict/alignment/distance_exclusions.json)",
         "esr_analysis/calculate_rmsd.py", "no"),
        ("distance pruning, keep-all matrix (all eligible chains)", len(eligible), len(eligible), 0,
         f"min_pairwise_common_ca {MIN_COMMON}; pairs below {MIN_COMMON}: {pairs_below}; minimum pairwise overlap {min_overlap}",
         "esr_analysis/calculate_rmsd.py", "no"),
        ("residue: LBD positions with no CA coordinate", lbd_expected_total, lbd_observed_total,
         lbd_expected_total - lbd_observed_total, "positions without modelled CA cannot be fitted or measured; chains kept",
         "esr_analysis/extract_lbd.py; esr_analysis/align_structures.py", "yes"),
        ("residue: CA B-factor invalid (NaN or negative), eligible chains", int(len(residues)), int(len(residues) - invalid_b), invalid_b,
         "B-factor only; coordinates kept; K03 leaves these blank", "esr_analysis/_regional_bfactors.py (old)", "no"),
        ("apo / no pocket ligand (labelled, Q18)", len(eligible), len(eligible), 0,
         f"status apo_or_unclear labelled, not removed ({apo_chains} chains)", "scripts/ligand_of_interest.py", "no"),
        ("QC flags (resolution, R-free, mean B)", len(entries), len(entries), 0,
         "flags only; combined pass needs all three; nothing deleted", "esr_analysis/flag_quality.py", "no"),
        ("old sensitivity drop: R-free above limit (strict set)", len(pick_keys), len(pick_keys) - qc_fail_selected, qc_fail_selected,
         f"R-free above {R_FREE_LIMIT}; used only by the strict S13 comparison", "scripts/qc_clustering_delta.py (archived)", "no"),
        ("regional-B coverage at 0.8 (old block 13, full cohort)", int(full.cohort_chains.astype(int).sum()),
         int(full.retained_chains.astype(int).sum()), b_dropped,
         "unit: chain x interval pairs (10 intervals); within-chain regional B distributions only; the chain stays in every matrix",
         "esr_analysis/_regional_bfactors.py", "no"),
    ]
    funnel = pd.DataFrame(steps, columns=["step", "n_before", "n_after", "n_removed", "rule", "code_file", "necessary"])
    funnel.to_csv(ROOT / "results/attrition_funnel.csv", index=False)

    loss = esr1_chains[["chain_key", "entry_id", "chain_id", "auth_asym_id", "eligible_historical", "historical_pick",
                        "mapped_lbd_ca", "mapped_core_ca"]].copy()
    loss["stage_lost"] = ""
    loss["reason"] = ""
    is_eligible = loss.eligible_historical.pipe(truthy)
    loss.loc[~is_eligible, "stage_lost"] = "completeness_screen"
    loss.loc[~is_eligible, "reason"] = loss.loc[~is_eligible, "chain_key"].map(reasons)
    dropped_pick = is_eligible & ~loss.historical_pick
    loss.loc[dropped_pick, "stage_lost"] = "one_chain_per_entry"
    loss.loc[dropped_pick, "reason"] = loss.loc[dropped_pick, "entry_id"].map(
        {entry: f"eligible; historical chain of entry is {key}" for entry, key in pick_keys.items()}) + "; this chain is not it"
    loss["in_old_matrix"] = loss.historical_pick
    loss["stage_lost_keep_all"] = np.where(is_eligible, "", "completeness_screen")
    loss["reason_keep_all"] = np.where(is_eligible, "", loss.chain_key.map(reasons).fillna(""))
    loss = loss.sort_values("chain_key").reset_index(drop=True)
    loss.to_csv(ROOT / "results/chain_loss_reasons.csv", index=False)

    print(funnel[["step", "n_before", "n_after", "n_removed", "necessary"]].to_string(index=False))
    print(f"chain_loss_reasons.csv rows: {len(loss)}; in old matrix: {int(loss.in_old_matrix.sum())}; "
          f"completeness out: {int((loss.stage_lost == 'completeness_screen').sum())}; one-per-entry out: {int((loss.stage_lost == 'one_chain_per_entry').sum())}")


if __name__ == "__main__":
    main()
