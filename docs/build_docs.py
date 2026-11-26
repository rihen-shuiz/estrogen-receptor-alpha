"""Build the PDF documents in docs/ from their HTML sources.

    python docs/build_docs.py      # -> docs/report_month3.pdf, docs/project_summary_full.pdf

The figures are redrawn without their titles (ESR1_NO_TITLE=1) into results/doc_figures/, where the HTML looks for
them; the PyMOL panels come from results/12_curated_action/pymol/. The PDF is printed by headless Google Chrome
($CHROME, or the default macOS location).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG_OUT = ROOT / "results/doc_figures"
DOCS = ["report_month3", "project_summary_full"]
FIGURES = ["h12_window", "lattice_within_crystal", "midpoint_error", "lid_outliers", "maps_three", "recover_lid_groups",
           "helix_pair_shift", "crystal_components", "assay_spread", "curated_labels", "curated_class_h12",
           "curated_class_map", "curated_h12_within", "curated_effects", "curated_same_ligand", "curated_pocket",
           "curated_steroid_map", "curated_steroid_class", "curated_redundancy", "curated_qc", "esr2_displacement",
           "esr2_map", "dendrogram", "rmsd_blocks", "where_groups_differ", "set_cover", "ligand_size", "esr2_blocks",
           "druglikeness_size", "lid_outlier_copies", "h12_length_mobility"]
CHROME = os.environ.get("CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")


def main():
    env = {**os.environ, "ESR1_NO_TITLE": "1", "ESR1_FIG_DIR": str(FIG_OUT), "MPLBACKEND": "Agg"}
    for name in FIGURES:
        subprocess.run([sys.executable, "-B", f"fig_{name}.py"], cwd=ROOT / "pipeline/figures", env=env, check=True,
                       stdout=subprocess.DEVNULL)
    for doc in DOCS:
        html, pdf = ROOT / "docs" / f"{doc}.html", ROOT / "docs" / f"{doc}.pdf"
        subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer", "--allow-file-access-from-files",
                        f"--print-to-pdf={pdf}", html.as_uri()], check=True, capture_output=True)
        print("wrote", pdf.relative_to(ROOT))


if __name__ == "__main__":
    main()
