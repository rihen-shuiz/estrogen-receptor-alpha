"""Tracker Y35: every number used in the final-run documents, read from result files (extends collect_numbers.py).

1. Collect: every results/ext/numbers_<task>.json written by the ext_* scripts is flattened into keys such as
   "Y04.best_run" or "Y14.peptide_entry.effect"; a set of headline numbers is recounted here directly from the CSVs
   (independent of the producing scripts) under the prefix "R." (e.g. "R.n_chains", "R.cl7_size").
   Output: results/ext/numbers_ext.json ({key: {"value": v, "source": file}}) and docs/numbers_ext.md.
2. Render: every docs/templates/<NAME>.md is rendered to docs/<NAME>.md. A placeholder {{key}} or {{key|fmt}}
   (fmt = a Python format spec such as .2f, .0%, d) is replaced by the number; an unknown key stops the script, so no
   document can contain a number that is not in numbers_ext.json.
Called at the end of scripts/collect_numbers.py as well.
Run from the repo root: PYTHONPATH=. python scripts/collect_numbers_ext.py
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/ext"


def flatten(prefix, obj, out, source):
    if isinstance(obj, dict):
        for k, v in obj.items():
            flatten(f"{prefix}.{k}", v, out, source)
    elif isinstance(obj, list) and obj and all(isinstance(x, (int, float, str)) for x in obj) and len(obj) <= 12:
        out[prefix] = {"value": ", ".join(str(x) for x in obj), "source": source}
    elif isinstance(obj, (int, float, str, bool)) or obj is None:
        out[prefix] = {"value": obj, "source": source}


def recount(out):
    def put(k, v, src):
        out[f"R.{k}"] = {"value": v.item() if hasattr(v, "item") else v, "source": src}
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    d = np.load(ROOT / "results/rmsd_entry.npy")
    cl = pd.read_csv(ROOT / "results/entry_clusters.csv")
    put("n_chains", len(idx), "results/rmsd_entry_index.csv")
    put("n_entries", idx.entry_id.nunique(), "results/rmsd_entry_index.csv")
    vc = cl["cut_1.5"].value_counts()
    put("cl7_size", int(vc[7]), "results/entry_clusters.csv")
    put("cl4_size", int(vc[4]), "results/entry_clusters.csv")
    put("minor_size", int(len(cl) - vc[7] - vc[4]), "results/entry_clusters.csv")
    put("n_clusters_15", int(cl["cut_1.5"].nunique()), "results/entry_clusters.csv")
    m7, m4 = (cl["cut_1.5"] == 7).values, (cl["cut_1.5"] == 4).values
    put("median_cl7_cl4_A", round(float(np.median(d[np.ix_(m7, m4)])), 2), "results/rmsd_entry.npy")
    cs = pd.read_csv(ROOT / "results/chain_status.csv")
    put("n_chains_all_esr1", len(cs), "results/chain_status.csv")
    put("n_chains_out_of_matrix", int((~cs.in_matrix).sum()), "results/chain_status.csv")
    em = pd.read_csv(RES / "entry_metadata.csv")
    put("n_peptide_entries", int(em.peptide_present.sum()), "results/ext/entry_metadata.csv")
    put("n_peg_entries", int((em.additive_class == "PEG").sum()), "results/ext/entry_metadata.csv")
    cm = pd.read_csv(RES / "chain_metadata.csv").merge(cl[["chain_key", "cut_1.5"]], on="chain_key")
    put("cl7_with_peptide", int(cm[(cm["cut_1.5"] == 7)].peptide_present.sum()), "results/ext/chain_metadata.csv")
    put("cl4_with_peptide", int(cm[(cm["cut_1.5"] == 4)].peptide_present.sum()), "results/ext/chain_metadata.csv")
    put("cl7_Y537S", int(((cm["cut_1.5"] == 7) & (cm.mutation_status == "Y537S")).sum()), "results/ext/chain_metadata.csv")
    put("cl4_L536S", int(((cm["cut_1.5"] == 4) & (cm.mutation_status == "L536S")).sum()), "results/ext/chain_metadata.csv")
    put("cl7_L536S", int(((cm["cut_1.5"] == 7) & (cm.mutation_status == "L536S")).sum()), "results/ext/chain_metadata.csv")
    z = pd.read_csv(RES / "zoo_runs.csv")
    put("zoo_n_runs", len(z), "results/ext/zoo_runs.csv")
    put("zoo_n_families", int(z.method.nunique()), "results/ext/zoo_runs.csv")
    cz = pd.read_csv(RES / "chem_zoo_runs.csv")
    put("chem_n_runs", len(cz), "results/ext/chem_zoo_runs.csv")
    cards = pd.read_csv(RES / "figure_cards.csv")
    put("n_ext_figures", int(cards.file.str.endswith(".png").sum()), "results/ext/figure_cards.csv")
    put("n_all_figures", int(cards.file.str.endswith(".png").sum() + len(list((ROOT / "figures").glob("*.png")))),
        "results/ext/figure_cards.csv + figures/*.png")
    st = pd.read_csv(RES / "structural_features.csv")
    for c in (7, 4):
        g = st[st.cluster == c]
        put(f"h12_shift_cl{c}", round(float(g.h12_shift_A.median()), 1), "results/ext/structural_features.csv")
        put(f"k362_e542_cl{c}", round(float(g.k362_e542_A.median()), 1), "results/ext/structural_features.csv")
        put(f"h12_h3_angle_cl{c}", round(float(g.h12_h3_angle_deg.median()), 0), "results/ext/structural_features.csv")
    sh = pd.read_csv(RES / "split_hypotheses.csv").set_index(["hypothesis", "level"])
    for key, hyp, lev in (("or_peptide_entry", "coactivator/other peptide present", "entry"),
                          ("or_Y537S_entry", "Y537S construct", "entry"), ("or_L536S_entry", "L536S construct", "entry"),
                          ("or_P21_entry", "space group P21", "entry"), ("or_amine_ligand", "ligand has a basic amine", "ligand")):
        put(key, round(float(sh.loc[(hyp, lev), "effect"]), 1), "results/ext/split_hypotheses.csv")
    q = sh.loc[("quality adds to peptide + mutation", "chain")]
    put("quality_gain_r2", round(float(q.effect), 3), "results/ext/split_hypotheses.csv")
    put("quality_gain_p", float(q.p_value), "results/ext/split_hypotheses.csv")
    for nm, hyp in (("res", "resolution (A)"), ("meanb", "mean B (A^2)"), ("ph", "crystallisation pH")):
        det = sh.loc[(hyp, "entry"), "detail"]
        a, b = re.findall(r"cl7 ([\d.]+), cl4 ([\d.]+)", det)[0]
        put(f"{nm}_cl7_entry", float(a), "results/ext/split_hypotheses.csv")
        put(f"{nm}_cl4_entry", float(b), "results/ext/split_hypotheses.csv")
    dt = pd.read_csv(RES / "ligand_descriptor_tests.csv").set_index("variable")
    put("amine_cl7", round(float(dt.loc["basic_amine (Fisher)", "median_cl7"]) * 100), "results/ext/ligand_descriptor_tests.csv")
    put("amine_cl4", round(float(dt.loc["basic_amine (Fisher)", "median_cl4"]) * 100), "results/ext/ligand_descriptor_tests.csv")
    put("heavy_cl7", float(dt.loc["heavy_atoms", "median_cl7"]), "results/ext/ligand_descriptor_tests.csv")
    put("heavy_cl4", float(dt.loc["heavy_atoms", "median_cl4"]), "results/ext/ligand_descriptor_tests.csv")
    put("mw_cl7", round(float(dt.loc["mw", "median_cl7"])), "results/ext/ligand_descriptor_tests.csv")
    put("mw_cl4", round(float(dt.loc["mw", "median_cl4"])), "results/ext/ligand_descriptor_tests.csv")
    sat = pd.read_csv(ROOT / "results/saturation.csv").set_index("cut")
    for c in sat.index:
        put(f"n90_cut{str(c).replace('.', 'p')}", int(sat.loc[c, "n_90"]), "results/saturation.csv")
    ent = pd.read_csv(ROOT / "results/entropy.csv").set_index("cut")
    for c in ent.index:
        put(f"entropy_cut{str(c).replace('.', 'p')}", round(float(ent.loc[c, "entropy_bits"]), 2), "results/entropy.csv")
    cst = pd.read_csv(ROOT / "results/cluster_stats.csv").set_index("cut")
    for c in cst.index:
        put(f"nclusters_cut{str(c).replace('.', 'p')}", int(cst.loc[c, "n_clusters"]), "results/cluster_stats.csv")
    sr = pd.read_csv(RES / "selectivity_ratios.csv")
    put("sel_cl7_median", round(float(sr[sr.structural == "cl7"].selectivity.median()), 2), "results/ext/selectivity_ratios.csv")
    put("sel_cl4_median", round(float(sr[sr.structural == "cl4"].selectivity.median()), 2), "results/ext/selectivity_ratios.csv")
    tiers = pd.read_csv(ROOT / "results/chain_status.csv")
    tiers = tiers[tiers.in_matrix]
    for lev in ("A", "B", "C"):
        put(f"tier_{lev}_chains", int((tiers.qc_label == lev).sum()), "results/chain_status.csv")
    put("rebuild_identical", int((pd.read_csv(ROOT / "results/rebuild/rebuild_report.csv").status == "identical").sum()),
        "results/rebuild/rebuild_report.csv")
    put("rebuild_equivalent", int((pd.read_csv(ROOT / "results/rebuild/rebuild_report.csv").status == "equivalent").sum()),
        "results/rebuild/rebuild_report.csv")
    put("rebuild_differs", int((pd.read_csv(ROOT / "results/rebuild/rebuild_report.csv").status == "differs").sum()),
        "results/rebuild/rebuild_report.csv")
    settings = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
    for k in ("min_chain_lbd_ca", "min_core_ca_for_alignment", "min_pairwise_common_ca", "random_seed", "qc_resolution_A",
              "qc_r_free", "qc_mean_B_A2", "contact_cutoff_A", "peptide_max_length", "esr1_lbd_start", "esr1_lbd_end",
              "primary_structural_threshold_A"):
        put(f"set_{k}", settings[k], "esr_analysis/resources/config/analysis_settings.json")
    put("set_h12_range", f"{settings['h12_range'][0]}-{settings['h12_range'][1]}", "esr_analysis/resources/config/analysis_settings.json")
    v3 = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    put("core_v3_ranges", "; ".join(f"{a}-{b}" for a, b in v3["core_ranges"]), "esr_analysis/resources/config/alignment_core_v3.json")
    put("core_v3_n", int(v3["residues_total"]), "esr_analysis/resources/config/alignment_core_v3.json")
    put("old_core_ranges", "; ".join(f"{a}-{b}" for a, b in settings["esr1_core_ranges"]),
        "esr_analysis/resources/config/analysis_settings.json")
    log = (ROOT / "tracker/esr1.log").read_text(encoding="utf8").splitlines()
    state = {}
    for line in log:
        p = line.split(" | ")
        if len(p) > 2 and p[1] != "-":
            state[p[1]] = p[2]
    tasks = re.findall(r"^T (\w\d\d) ", (ROOT / "tracker/esr1.tasks").read_text(encoding="utf8"), re.M)
    for s in ("DONE", "PART", "BLOCK", "MANUAL"):
        put(f"tasks_{s}", sum(state.get(t) == s for t in tasks), "tracker/esr1.log")
    put("tasks_TODO", sum(t not in state for t in tasks), "tracker/esr1.log")
    put("tasks_total", len(tasks), "tracker/esr1.tasks")


def fmt(v, spec):
    if spec is None or v is None or isinstance(v, str):
        return str(v)
    if spec == "d":
        return str(int(round(float(v))))
    return format(float(v), spec)


def render(nums):
    tdir = ROOT / "docs/templates"
    done = []
    for tmpl in sorted(tdir.glob("*.md")):
        txt = tmpl.read_text(encoding="utf8")
        missing = []

        def sub(m):
            key, spec = m.group(1), m.group(2)
            if key not in nums:
                missing.append(key)
                return m.group(0)
            return fmt(nums[key]["value"], spec)
        out = re.sub(r"\{\{([A-Za-z0-9_.\-|/ ]+?)(?:\|([^}]+))?\}\}", sub, txt)
        if missing:
            raise KeyError(f"{tmpl.name}: unknown number keys {sorted(set(missing))}")
        header = f"<!-- generated by scripts/collect_numbers_ext.py from docs/templates/{tmpl.name}; numbers from results/ext/numbers_ext.json -->\n"
        (ROOT / "docs" / tmpl.name).write_text(header + out, encoding="utf8")
        done.append(tmpl.name)
    return done


def main():
    nums = {}
    for f in sorted(RES.glob("numbers_Y*.json")):
        flatten(f.stem.replace("numbers_", ""), json.loads(f.read_text()), nums, f"results/ext/{f.name}")
    recount(nums)
    (RES / "numbers_ext.json").write_text(json.dumps(nums, indent=1, default=str))
    lines = ["# Numbers for the final-run documents", "",
             "Generated by `scripts/collect_numbers_ext.py`. Keys `R.*` are recounted directly from the CSVs; the others",
             "are written by the producing `scripts/ext_*.py`.", "", "| key | value | source |", "|---|---|---|"]
    for k, v in nums.items():
        val = str(v["value"]).replace("|", "/")
        lines.append(f"| {k} | {val[:120]} | {v['source']} |")
    (ROOT / "docs/numbers_ext.md").write_text("\n".join(lines) + "\n", encoding="utf8")
    done = render(nums) if (ROOT / "docs/templates").exists() else []
    print(f"{len(nums)} numbers -> results/ext/numbers_ext.json; rendered: {done}")


if __name__ == "__main__":
    main()
