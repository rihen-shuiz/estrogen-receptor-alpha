"""Action-class styles for SV13-SV15 figures (RUN_SV1028 rule 1.8).

Classes come only from results/sv/SV12/action_class.csv (rule 1.9). Unlabelled ligands are grey. Class colours never
reuse the cluster colours (cl7 blue, cl4 orange, minor clusters aqua; scripts/sv/common/sv_fig.CLUSTER_COLOR), which
is checked at import. Created by SV14 (it did not exist before).

Usage:
    from common.class_style import CLASS_ORDER, CLASS_COLOR, CLASS_MARKER, CLASS_LABEL, style
    st = style("SERM")   # {'color': ..., 'marker': ..., 'label': ...}
"""
from __future__ import annotations

from .sv_fig import CLUSTER_COLOR, GREY

CLASS_ORDER = ("agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other", "unlabelled")
CLASS_COLOR = {
    "agonist": "#4a3aa7",          # violet (project palette slot VIOLET)
    "partial_agonist": "#9e8fd6",  # light violet
    "SERM": "#c2185b",             # dark magenta
    "SERD": "#5d4037",             # brown
    "antagonist": "#0b0b0b",       # ink
    "other": "#8a8a00",            # olive
    "unlabelled": GREY,            # grey, rule 1.8
}
CLASS_MARKER = {"agonist": "^", "partial_agonist": "<", "SERM": "s", "SERD": "D", "antagonist": "v", "other": "P",
                "unlabelled": "o"}
CLASS_LABEL = {"agonist": "agonist", "partial_agonist": "partial agonist", "SERM": "SERM", "SERD": "SERD",
               "antagonist": "antagonist", "other": "other", "unlabelled": "unlabelled (no class)"}

_cluster = {c.lower() for c in CLUSTER_COLOR.values()}
_clash = [k for k, c in CLASS_COLOR.items() if c.lower() in _cluster]
if _clash:
    raise ValueError(f"class colours reuse cluster colours: {_clash}")


def style(cls: str) -> dict:
    """Colour, marker and legend label of an action class; anything unknown is styled as unlabelled."""
    key = cls if cls in CLASS_COLOR else "unlabelled"
    return {"color": CLASS_COLOR[key], "marker": CLASS_MARKER[key], "label": CLASS_LABEL[key]}
