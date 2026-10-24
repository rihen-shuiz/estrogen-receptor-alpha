"""Rebuild tracker E03 (ESR1<->ESR2 numbering map) and E04 (ligands matched across receptors).

Inputs: UniProt sequences P03372 / Q92731 in esr_analysis/resources/data, settings; for E04
results/rebuild/entries.csv (ligand_ids), results/rebuild/fp_index.csv + tanimoto_ecfp4.npy (C02/C04
rebuild), data/pdb_ids_esr1.txt and data/pdb_ids_esr2.txt.
Method E03: esr_analysis.map_residues.receptor_equivalence (Biopython global alignment, BLOSUM62, gap
-10/-0.5) of the full canonical sequences; one row per aligned pair, sorted by ESR1 position; in_lbd =
ESR1 position in esr1_lbd_start..esr1_lbd_end (305-550).
Method E04: for every ESR1-list entry x ESR2-list entry, compare their ligands except the
solvent_additive_codes of analysis_settings.json; keep the pair if the best ligand-ligand ECFP4 Tanimoto
is > 0.9 and record that best ligand pair (first in sorted ligand order on ties); match_type 'identical'
when the CCD codes are equal, else 'similar'.
Outputs: results/rebuild/numbering_map_esr1_esr2.csv, results/rebuild/matched_pairs.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_esr1_esr2_maps.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from esr_analysis._task_io import RESOURCE_ROOT, sequences
from esr_analysis.map_residues import receptor_equivalence
from rebuild_common import OUT, ROOT, SETTINGS, esr2_ids, read_csv


def numbering_map():
    seq = sequences(RESOURCE_ROOT, SETTINGS)
    pairs = sorted((e1, e2) for e2, e1 in receptor_equivalence(seq).items())
    lbd = range(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)
    rows = [{"esr1_pos": e1, "esr2_pos": e2, "esr1_aa": seq["esr1"][e1 - 1], "esr2_aa": seq["esr2"][e2 - 1],
             "same_residue": seq["esr1"][e1 - 1] == seq["esr2"][e2 - 1], "in_lbd": e1 in lbd} for e1, e2 in pairs]
    pd.DataFrame(rows).to_csv(OUT / "numbering_map_esr1_esr2.csv", index=False)
    print(f"E03 rebuild: {len(rows)} aligned positions, {sum(r['in_lbd'] for r in rows)} in LBD")


def matched_pairs():
    entries = read_csv(OUT / "entries.csv")
    index = read_csv(OUT / "fp_index.csv")
    tanimoto = np.load(OUT / "tanimoto_ecfp4.npy")
    row = {code: i for i, code in enumerate(index.ligand_id.astype(str))}
    skip = set(SETTINGS["solvent_additive_codes"])
    ligands = {e.entry_id: [c for c in str(e.ligand_ids if isinstance(e.ligand_ids, str) else "").split(";")
                            if c and c not in skip and c in row] for e in entries.itertuples()}
    esr1 = [line.strip() for line in (ROOT / "data/pdb_ids_esr1.txt").read_text().splitlines() if line.strip()]
    esr2 = esr2_ids()
    rows = []
    for a in esr1:
        for b in esr2:
            best = None
            for x in ligands[a]:
                for y in ligands[b]:
                    t = tanimoto[row[x], row[y]]
                    if t > 0.9 and (best is None or t > best[2]):
                        best = (x, y, t)
            if best:
                rows.append({"esr1_entry": a, "esr2_entry": b, "esr1_ligand": best[0], "esr2_ligand": best[1],
                             "tanimoto": float(best[2]), "match_type": "identical" if best[0] == best[1] else "similar"})
    table = pd.DataFrame(rows, columns=["esr1_entry", "esr2_entry", "esr1_ligand", "esr2_ligand", "tanimoto", "match_type"])
    table.to_csv(OUT / "matched_pairs.csv", index=False)
    print(f"E04 rebuild: {len(table)} entry pairs; {table.match_type.value_counts().to_dict()}")


def main():
    numbering_map()
    matched_pairs()


if __name__ == "__main__":
    main()
