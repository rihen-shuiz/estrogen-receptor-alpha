"""Tracker Y19: interactive, offline HTML scatter of the RMSD embeddings (plotly, library embedded in the file).

Input: results/ext/colour_by_table.csv (Y18; every matrix chain with its labels and metrics).
One trace per structural cluster (shared palette), hover shows PDB id, chain, cluster, tier, resolution, R-free,
mean B, H12 label, ligand, peptide, construct mutation, space group, outlier count. A dropdown switches the map
(metric MDS, PCA of RMSD rows, classical MDS, t-SNE, UMAP).
Output: docs/interactive_embedding.html (self-contained, opens without network).
Run from the repo root: PYTHONPATH=. python scripts/ext_interactive.py
"""
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import CLUSTER_COLOR, RES, ROOT, cluster_name, register_card

MAPS = [("metric MDS", "mds1", "mds2"), ("PCA of RMSD rows", "pca1", "pca2"), ("classical MDS", "cmds1", "cmds2"),
        ("t-SNE", "tsne1", "tsne2"), ("UMAP", "umap1", "umap2")]
SYMBOL = {7: "circle", 4: "circle", 1: "triangle-up", 2: "square", 3: "diamond", 5: "triangle-down", 6: "cross"}


def main():
    t = pd.read_csv(RES / "colour_by_table.csv")
    t["chain"] = t.chain_key.str.split("|").str[1]
    hover = ("<b>%{customdata[0]}</b> chain %{customdata[1]}<br>cluster %{customdata[2]} | tier %{customdata[3]}"
             "<br>resolution %{customdata[4]} A | R-free %{customdata[5]} | mean B %{customdata[6]} A^2"
             "<br>H12 %{customdata[7]} | ligand %{customdata[8]}<br>peptide %{customdata[9]}"
             "<br>construct %{customdata[10]} | %{customdata[11]}<br>outlier detectors: %{customdata[12]}<extra></extra>")
    cols = ["entry_id", "chain", "cluster", "tier", "resolution_A", "r_free", "mean_b_A2", "h12_label", "ligand_id",
            "peptide_class", "mutation_status", "space_group", "n_methods"]
    fig = go.Figure()
    clusters = sorted(t.cluster.unique(), key=lambda c: (c not in (7, 4), c))
    for c in clusters:
        d = t[t.cluster == c]
        fig.add_trace(go.Scattergl(x=d.mds1, y=d.mds2, mode="markers", name=f"{cluster_name(c)} (n={len(d)})",
                                   marker=dict(color=CLUSTER_COLOR[c], symbol=SYMBOL[c], size=6 if c in (7, 4) else 10,
                                               line=dict(width=0.5, color="white")),
                                   customdata=d[cols].round(3).astype(str).values, hovertemplate=hover))
    buttons = []
    for name, x, y in MAPS:
        buttons.append(dict(label=name, method="update",
                            args=[{"x": [t[t.cluster == c][x] for c in clusters],
                                   "y": [t[t.cluster == c][y] for c in clusters]},
                                  {"xaxis.title.text": f"{name} 1", "yaxis.title.text": f"{name} 2"}]))
    fig.update_layout(
        title=dict(text="ESR1 LBD conformations: 788 chains, 5 maps of the same C-alpha RMSD matrix "
                        "(hover for details; click legend to hide a cluster)", x=0.01),
        updatemenus=[dict(buttons=buttons, direction="down", x=1.0, xanchor="right", y=1.12, yanchor="top")],
        xaxis_title="metric MDS 1", yaxis_title="metric MDS 2", template="simple_white", height=760,
        legend=dict(title="S08 cluster (1.5 A)"), font=dict(family="Arial", size=12))
    out = ROOT / "docs/interactive_embedding.html"
    fig.write_html(out, include_plotlyjs=True, full_html=True)
    register_card("Y19_interactive_embedding", section="F. Embeddings and colour-by", task="Y19",
                  shows="interactive scatter (plotly, offline) of five embeddings with hover: PDB id, chain, cluster, quality, ligand, peptide, mutation",
                  takeaway="Open in a browser to look up any chain; the dropdown switches between MDS, PCA, t-SNE and UMAP.",
                  source="results/ext/colour_by_table.csv", script="scripts/ext_interactive.py",
                  file="docs/interactive_embedding.html")
    print(out, out.stat().st_size)


if __name__ == "__main__":
    main()
