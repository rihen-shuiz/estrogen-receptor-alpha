"""VIS1 - cut multi-graph PNGs into one PNG per graph (crops of the original pixels only).

Steps (RUN_SV1028 VIS1):
1. Inventory: every PNG except data/, figures/panels/, visuals/, virtual environments and hidden folders;
   SHA-1 groups are cut once; images < 200 px on a side are skipped.  -> results/sv/VIS1/inventory.csv
2. Classification and parameters: scripts/sv/VIS1_spec.py.
3. Method (a) SV figures: panels already saved by the SV scripts in figures/panels/sv/<ID>/ are listed.
   Method (b) (render the SVG with other axes hidden) is unavailable: cairosvg cannot load libcairo,
   rsvg-convert is absent (see results/sv/VIS1/renderer_check.txt). All other multi images use (c): recursive
   XY-cut along background gutters of the PNG. For matplotlib figures with an SVG twin the cut positions are chosen
   with the axes layout read from the SVG (top-level <g id="axes_N"> rectangles and text anchors): a cut must not
   cross any graph, and the widest gutter that separates graphs is used. Strips without a graph (suptitle) are dropped.
4. Panels: figures/panels/<parent>__<stem>__<nn>_<slug>.png (<parent> = parent folder path, '/' -> '_'),
   manifest figures/panels/manifest.csv, contact sheets figures/panels/_check/.
Run: .venv/bin/python -B scripts/sv/VIS1_cut.py
"""
from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import VIS1_spec as spec  # noqa: E402

PANELS = ROOT / "figures/panels"
CHECK = PANELS / "_check"
RES = ROOT / "results/sv/VIS1"
MIN_SIDE = 200
PAD = 10
BG_TOL = 8
NS = "{http://www.w3.org/2000/svg}"
FONT = Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans-Bold.ttf"
PALETTE = [(214, 39, 40), (31, 119, 180), (44, 160, 44), (148, 103, 189), (255, 127, 14), (23, 190, 207),
           (227, 119, 194), (140, 86, 75), (188, 189, 34)]


# ----------------------------------------------------------------------------------------------- inventory
def inventory():
    rows = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        rel = Path(dirpath).relative_to(ROOT)
        keep = []
        for d in dirnames:
            r = (rel / d).as_posix()
            if (d.startswith(".") or r in ("data", "visuals", "figures/panels")
                    or (Path(dirpath) / d / "pyvenv.cfg").exists()):
                continue
            keep.append(d)
        dirnames[:] = sorted(keep)
        for f in sorted(filenames):
            if not f.lower().endswith(".png") or f.startswith("."):
                continue
            p = Path(dirpath) / f
            w, h = Image.open(p).size
            rows.append(dict(path=p.relative_to(ROOT).as_posix(), sha1=hashlib.sha1(p.read_bytes()).hexdigest(),
                             width=w, height=h, svg_twin=p.with_suffix(".svg").exists(),
                             small=min(w, h) < MIN_SIDE))
    rows.sort(key=lambda r: r["path"])
    return rows


# ----------------------------------------------------------------------------------------------- pixels
def load_rgb(p: Path) -> np.ndarray:
    im = Image.open(p)
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        bg.alpha_composite(im)
        im = bg
    return np.asarray(im.convert("RGB")).astype(np.int16)


def blank_mask(a: np.ndarray):
    """Background = white plus light colours (all channels >= 240) covering >= 3 % of the image."""
    sub = a[::4, ::4].reshape(-1, 3)
    cols, cnt = np.unique(sub, axis=0, return_counts=True)
    bgs = [np.array([255, 255, 255], dtype=np.int16)]
    for c, n in zip(cols, cnt):
        if n >= 0.03 * len(sub) and (c >= 240).all() and not (c == 255).all():
            bgs.append(c.astype(np.int16))
    m = np.zeros(a.shape[:2], bool)
    for c in bgs:
        m |= np.abs(a - c).max(axis=2) <= BG_TOL
    return m, [tuple(int(v) for v in c) for c in bgs]


def runs(b: np.ndarray):
    """(start, end) of True runs in a 1-D boolean array (end exclusive)."""
    if not b.any():
        return []
    d = np.diff(np.concatenate([[0], b.astype(np.int8), [0]]))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def trim(blank, box):
    x0, y0, x1, y1 = box
    ink = ~blank[y0:y1, x0:x1]
    if not ink.any():
        return None
    r = np.where(ink.any(1))[0]
    c = np.where(ink.any(0))[0]
    return (x0 + int(c[0]), y0 + int(r[0]), x0 + int(c[-1]) + 1, y0 + int(r[-1]) + 1)


# ----------------------------------------------------------------------------------------------- SVG layout
def _nums(s):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*(?:e-?\d+)?", s)]


def _rect_patch(g):
    """Rectangle of the first patch child of a group (axes background / legend frame), or None."""
    for ch in g:
        cid = ch.get("id", "")
        if cid.startswith("patch"):
            path = ch.find(f"{NS}path")
            if path is None:
                return None
            v = _nums(path.get("d"))
            if len(v) < 8:
                return None
            xs, ys = v[0::2], v[1::2]
            return (min(xs), min(ys), max(xs), max(ys))
        return None
    return None


def _anchor(t):
    x, y = t.get("x"), t.get("y")
    if x is None:
        ts = t.find(f"{NS}tspan")
        if ts is not None:
            x, y = ts.get("x"), ts.get("y")
    try:
        return float(x.split()[0]), float(y.split()[0])
    except (AttributeError, ValueError):
        return None


def _text(t):
    return " ".join("".join(t.itertext()).split())


def _translate(el):
    """Translation part of an SVG transform attribute (translate or matrix); rotation/scale ignored."""
    tr = el.get("transform") or ""
    m = re.search(r"translate\(([^)]*)\)", tr)
    if m:
        v = _nums(m.group(1))
        return (v[0], v[1] if len(v) > 1 else 0.0)
    m = re.search(r"matrix\(([^)]*)\)", tr)
    if m:
        v = _nums(m.group(1))
        return (v[4], v[5])
    return (0.0, 0.0)


_FONTS = {}


def _font(bold: bool, size_px: float):
    key = (bold, max(1, round(size_px)))
    if key not in _FONTS:
        name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
        _FONTS[key] = ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts/ttf" / name), key[1])
    return _FONTS[key], size_px / key[1]


def _style(el, name, default=None):
    m = re.search(name + r"\s*:\s*([^;]+)", el.get("style") or "")
    return m.group(1).strip() if m else default


def _text_box(text, x, y, fs, bold, anchor, angle, scale):
    """Approximate pixel box of one text run (DejaVu Sans metrics, as matplotlib used for the PNG)."""
    sx, sy = scale
    font, k = _font(bold, fs * sx)
    if not text.strip():
        return None
    w = font.getlength(text) * k
    l, t, r, b = (v * k for v in font.getbbox(text, anchor="ls"))  # ink box relative to the baseline start
    f = {"start": 0.0, "middle": 0.5, "end": 1.0}.get(anchor, 0.0)
    a = math.radians(angle)
    pts = []
    for cx, cy in ((l - w * f, t), (r - w * f, t), (l - w * f, b), (r - w * f, b)):
        pts.append((x * sx + cx * math.cos(a) - cy * math.sin(a), y * sy + cx * math.sin(a) + cy * math.cos(a)))
    xs, ys = [q[0] for q in pts], [q[1] for q in pts]
    return (min(xs) - 2, min(ys) - 2, max(xs) + 2, max(ys) + 2)


def _texts(el, off=(0.0, 0.0), scale=None):
    """Yield (text element, anchor, pixel box or None) for every <text> below el, with accumulated group
    translations (anchor in SVG units; box needs scale = (sx, sy))."""
    if el.tag == f"{NS}text":
        a = _anchor(el)
        tr = _translate(el)  # multi-line texts are placed by transform="translate(x y)" without x/y
        if a is None and (tr != (0.0, 0.0) or el.findall(f"{NS}tspan")):
            a = (0.0, 0.0)
        if a:
            off = (off[0] + tr[0], off[1] + tr[1])
            box = None
            if scale:
                bold = _style(el, "font-weight", "") in ("700", "bold")
                m = re.search(r"rotate\((-?[\d.]+)", el.get("transform") or "")
                ang = float(m.group(1)) if m else 0.0
                spans = el.findall(f"{NS}tspan")
                boxes = []
                if spans:  # mathtext: one box per tspan (start anchored)
                    for ts in spans:
                        try:
                            tx, ty = float(ts.get("x")), float(ts.get("y"))
                        except (TypeError, ValueError):
                            continue
                        fs = _nums(_style(ts, "font-size", _style(el, "font-size", "10")))[0]
                        bx = _text_box(ts.text or "", tx + off[0], ty + off[1], fs, bold, "start", ang, scale)
                        if bx:
                            boxes.append(bx)
                else:
                    fs = _nums(_style(el, "font-size", "10"))[0]
                    bx = _text_box(_text(el), a[0] + off[0], a[1] + off[1], fs, bold,
                                   _style(el, "text-anchor", "start"), ang, scale)
                    if bx:
                        boxes.append(bx)
                if boxes:
                    box = (min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes),
                           max(b[3] for b in boxes))
            yield el, (a[0] + off[0], a[1] + off[1]), box
        return
    t = _translate(el) if el.tag == f"{NS}g" else (0.0, 0.0)
    off = (off[0] + t[0], off[1] + t[1])
    for ch in el:
        yield from _texts(ch, off, scale)


def _first_text(el):
    for t, a, _ in _texts(el):
        if _text(t):
            return _text(t), a
    return None, None


def svg_layout(svg: Path, png_w: int, png_h: int):
    root = ET.parse(svg).getroot()
    sw, sh = _nums(root.get("width"))[0], _nums(root.get("height"))[0]
    sx, sy = png_w / sw, png_h / sh
    fig = root.find(f"{NS}g")
    axes, fig_texts, fig_legends = [], [], []
    for g in fig:
        gid = g.get("id", "")
        if gid.startswith("axes"):
            rect = _rect_patch(g)
            tboxes = [b for _, _, b in _texts(g, scale=(sx, sy)) if b]
            legs = [_rect_patch(lg) for lg in g.iter(f"{NS}g") if lg.get("id", "").startswith("legend")]
            # title candidates: direct text children above the axes rectangle; labels: axis label texts
            title, ylabel, xlabel = None, None, None
            if rect:
                cands = []
                for ch in g:
                    if ch.get("id", "").startswith("text"):
                        txt, a = _first_text(ch)
                        if a and a[1] < rect[1]:
                            cands.append(txt)
                if cands:
                    title = max(cands, key=len)
            for ch in g:
                if ch.get("id", "").startswith("matplotlib.axis"):
                    for sub in ch:
                        if sub.get("id", "").startswith("text"):
                            txt, _ = _first_text(sub)
                            if txt:
                                if ch.get("id").endswith("_2") and ylabel is None:
                                    ylabel = txt
                                elif ch.get("id").endswith("_1") and xlabel is None:
                                    xlabel = txt
            sc = lambda r: (r[0] * sx, r[1] * sy, r[2] * sx, r[3] * sy)  # noqa: E731
            parts = list(tboxes) + ([sc(rect)] if rect else []) + [sc(lr) for lr in legs if lr]
            if not parts:
                continue
            box = (min(b[0] for b in parts), min(b[1] for b in parts), max(b[2] for b in parts),
                   max(b[3] for b in parts))
            axes.append(dict(id=gid, rect=sc(rect) if rect else None, box=box,
                             title=title, ylabel=ylabel, xlabel=xlabel))
        elif gid.startswith("text"):
            txt, a = _first_text(g)
            bx = [b for _, _, b in _texts(g, scale=(sx, sy)) if b]
            if txt:
                fig_texts.append(dict(id=gid, text=txt, xy=(a[0] * sx, a[1] * sy) if a else None,
                                      box=(min(b[0] for b in bx), min(b[1] for b in bx), max(b[2] for b in bx),
                                           max(b[3] for b in bx)) if bx else None))
        elif gid.startswith("legend"):
            lr = _rect_patch(g)
            if lr:
                fig_legends.append(dict(id=gid, box=(lr[0] * sx, lr[1] * sy, lr[2] * sx, lr[3] * sy)))
    return axes, fig_texts, fig_legends, (sx, sy)


def group_axes(axes, extra_groups=()):
    """Attach thin axes (colorbars, side bars) and overlapping axes (twins, insets) to their neighbour."""
    parent = {a["id"]: a["id"] for a in axes}

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def union(i, j):
        parent[find(i)] = find(j)

    R = {a["id"]: (a["rect"] or a["box"]) for a in axes}
    for a in axes:
        x0, y0, x1, y1 = R[a["id"]]
        w, h = x1 - x0, y1 - y0
        best = None
        for b in axes:
            if b is a:
                continue
            u0, v0, u1, v1 = R[b["id"]]
            bw, bh = u1 - u0, v1 - v0
            ov = max(0, min(x1, u1) - max(x0, u0)) * max(0, min(y1, v1) - max(y0, v0))
            if ov > 0.5 * min(w * h, bw * bh):  # twin / inset
                union(a["id"], b["id"])
                continue
            if w < 70 and h > 2 * w and bw >= 70:  # vertical thin axes (colorbar)
                vov = max(0, min(y1, v1) - max(y0, v0))
                gap = min(abs(x0 - u1), abs(u0 - x1))
                if vov >= 0.5 * h and gap <= 150 and (best is None or gap < best[0]):
                    best = (gap, b["id"])
            if h < 70 and w > 2 * h and bh >= 70:  # horizontal thin axes
                hov = max(0, min(x1, u1) - max(x0, u0))
                gap = min(abs(y0 - v1), abs(v0 - y1))
                if hov >= 0.5 * w and gap <= 150 and (best is None or gap < best[0]):
                    best = (gap, b["id"])
        if best:
            union(a["id"], best[1])
    for grp in extra_groups:
        for i in grp[1:]:
            union(i, grp[0])
    groups = {}
    for a in axes:
        groups.setdefault(find(a["id"]), []).append(a)
    items = []
    for members in groups.values():
        bx = [m["box"] for m in members]
        box = (min(b[0] for b in bx), min(b[1] for b in bx), max(b[2] for b in bx), max(b[3] for b in bx))
        main = max(members, key=lambda m: (lambda r: (r[2] - r[0]) * (r[3] - r[1]))(m["rect"] or m["box"]))
        label = main["title"] or main["ylabel"] or main["xlabel"] or main["id"]
        items.append(dict(kind="graph", ids=[m["id"] for m in members], box=box, title=label))
    return items


# ----------------------------------------------------------------------------------------------- cutting
TOL = 4  # px: estimated text boxes include a 2 px margin and overshoot the ink by 2-9 px (checked on Y16, Y09, Y26)


def _side(box, axis, mid, tol=0):
    lo, hi = (box[1], box[3]) if axis == "h" else (box[0], box[2])
    if hi - tol <= mid:
        return -1
    if lo + tol >= mid:
        return 1
    return 0


def guided_cut(blank, region, items, leaves, drops, fails):
    """Recursive XY-cut; a cut may not cross an item box. Prefers cuts with items on both sides (widest gutter)."""
    t = trim(blank, region)
    if t is None:
        return
    x0, y0, x1, y1 = t
    inside = [it for it in items if it["box"][0] < x1 and it["box"][2] > x0 and it["box"][1] < y1
              and it["box"][3] > y0]
    if not inside:
        drops.append(t)
        return
    sub = blank[y0:y1, x0:x1]
    best2, best1 = None, None
    for axis, line, off in (("h", sub.all(1), y0), ("v", sub.all(0), x0)):
        for a, b in runs(line):
            mid = off + (a + b) // 2
            width = b - a
            cand = (width, axis, mid)
            # cut between graphs: the line may graze an estimated text box by <= TOL px (it is background anyway)
            sides = [_side(it["box"], axis, mid, TOL) for it in inside]
            if 0 not in sides and -1 in sides and 1 in sides:
                if best2 is None or width > best2[0]:
                    best2 = cand
                continue
            # cut that only removes ink outside every graph (suptitle): strict, no tolerance
            sides = [_side(it["box"], axis, mid) for it in inside]
            if 0 not in sides and len(set(sides)) == 1 and (best1 is None or width > best1[0]):
                best1 = cand
    pick = best2 or best1
    if pick is None:
        if len(inside) == 1:
            leaves.append(dict(ink=t, allow=region, item=inside[0]))
        else:
            fails.append(dict(box=t, ids=[it["ids"] for it in inside]))
        return
    _, axis, mid = pick
    if axis == "h":
        parts = [(region[0], region[1], region[2], mid), (region[0], mid, region[2], region[3])]
    else:
        parts = [(region[0], region[1], mid, region[3]), (mid, region[1], region[2], region[3])]
    for p in parts:
        guided_cut(blank, p, items, leaves, drops, fails)


def gutter_cut(blank, region, g, leaves, big, only="hv"):
    """Plain XY-cut: split at the widest background gutter (>= g px) that leaves a graph-sized part on both
    sides (big(box) True); otherwise the region is one leaf. only: allowed cut directions ('h', 'v', 'hv')."""
    t = trim(blank, region)
    if t is None:
        return
    x0, y0, x1, y1 = t
    sub = blank[y0:y1, x0:x1]
    best = None
    for axis, line, off in (("h", sub.all(1), y0), ("v", sub.all(0), x0)):
        if axis not in only:
            continue
        for a, b in runs(line):
            if b - a < g or (best is not None and b - a <= best[0]):
                continue
            mid = off + (a + b) // 2
            if axis == "h":
                parts = [(x0, y0, x1, mid), (x0, mid, x1, y1)]
            else:
                parts = [(x0, y0, mid, y1), (mid, y0, x1, y1)]
            tp = [trim(blank, q) for q in parts]
            if all(q is not None and big(q) for q in tp):
                best = (b - a, axis, mid)
    if best is None:
        leaves.append(dict(ink=t, allow=region, item=None))
        return
    _, axis, mid = best
    if axis == "h":
        parts = [(region[0], region[1], region[2], mid), (region[0], mid, region[2], region[3])]
    else:
        parts = [(region[0], region[1], mid, region[3]), (mid, region[1], region[2], region[3])]
    for q in parts:
        gutter_cut(blank, q, g, leaves, big, only)


def longest_run(ink: np.ndarray) -> int:
    """Longest straight run of ink pixels along any row or column."""
    best = 0
    for arr in (ink, ink.T):
        d = np.diff(np.pad(arr.astype(np.int8), ((0, 0), (1, 1))), axis=1)
        st = np.argwhere(d == 1)
        en = np.argwhere(d == -1)
        if len(st):
            best = max(best, int((en[:, 1] - st[:, 1]).max()))
    return best


def edge_ink(path: Path, e: int = 4):
    """Ink pixels within e px of each edge of a panel PNG (own ink of a padded crop never touches the edge)."""
    a = load_rgb(path)
    ink = ~blank_mask(a)[0]
    d = dict(top=int(ink[:e].sum()), bottom=int(ink[-e:].sum()), left=int(ink[:, :e].sum()),
             right=int(ink[:, -e:].sum()))
    return {k: v for k, v in d.items() if v}


def foreign_ink(blank, crop, own, others, grow=6):
    """Ink pixels inside the crop that lie in another graph's or a figure-level text's estimated box and outside
    the panel's own estimated box (grown by `grow` px): fragments of neighbours."""
    x0, y0, x1, y1 = crop
    ink = ~blank[y0:y1, x0:x1]
    mine = np.zeros_like(ink)
    oth = np.zeros_like(ink)

    def paint(m, b, g):
        a0, b0 = max(0, int(b[0] - g) - x0), max(0, int(b[1] - g) - y0)
        a1, b1 = min(x1 - x0, int(math.ceil(b[2] + g)) - x0), min(y1 - y0, int(math.ceil(b[3] + g)) - y0)
        if a1 > a0 and b1 > b0:
            m[b0:b1, a0:a1] = True

    paint(mine, own, grow)
    for b in others:
        paint(oth, b, 0)
    return int((ink & oth & ~mine).sum())


def crop_box(leaf):
    (x0, y0, x1, y1), (a0, b0, a1, b1) = leaf["ink"], leaf["allow"]
    return (max(a0, x0 - PAD), max(b0, y0 - PAD), min(a1, x1 + PAD), min(b1, y1 + PAD))


def reading_order(leaves):
    """Rows by top edge (a graph joins the row if its top lies within the upper half of the row's first graph),
    then left to right."""
    leaves = sorted(leaves, key=lambda L: (L["ink"][1], L["ink"][0]))
    rows = []
    for L in leaves:
        if rows:
            f = rows[-1][0]["ink"]
            if L["ink"][1] < f[1] + 0.5 * (f[3] - f[1]):
                rows[-1].append(L)
                continue
        rows.append([L])
    return [L for r in rows for L in sorted(r, key=lambda L: (L["ink"][0], L["ink"][1]))]


def slugify(s: str, n: int = 40) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")
    if len(s) > n:
        cut = s[:n]
        s = cut[:cut.rfind("_")] if "_" in cut[n // 2:] else cut
    return s.strip("_") or "panel"


# ----------------------------------------------------------------------------------------------- contact sheet
def contact_sheet(orig: Path, panels, boxes, drops, out: Path, legend=None):
    font = ImageFont.truetype(str(FONT), 26)
    small = ImageFont.truetype(str(FONT), 20)
    im = Image.open(orig).convert("RGB")
    W, H = 2600, 1500
    LW = 1150
    s = min((LW - 20) / im.width, (H - 80) / im.height)
    left = im.resize((max(1, int(im.width * s)), max(1, int(im.height * s))), Image.LANCZOS)
    sheet = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(sheet)
    d.text((10, 8), f"original: {orig.relative_to(ROOT)}  ({im.width}x{im.height})", fill="black", font=small)
    sheet.paste(left, (10, 40))
    for k, b in enumerate(boxes):
        c = PALETTE[k % len(PALETTE)]
        bb = [10 + b[0] * s, 40 + b[1] * s, 10 + b[2] * s, 40 + b[3] * s]
        d.rectangle(bb, outline=c, width=3)
        d.text((bb[0] + 4, bb[1] + 2), f"{k + 1:02d}", fill=c, font=font)
    for b in drops:
        bb = [10 + b[0] * s, 40 + b[1] * s, 10 + b[2] * s, 40 + b[3] * s]
        d.rectangle(bb, outline=(150, 150, 150), width=1)
    items = list(panels) + ([legend] if legend else [])
    n = len(items)
    cols = max(1, math.ceil(math.sqrt(n * 1.3)))
    rows = max(1, math.ceil(n / cols))
    cw, ch = (W - LW - 20) // cols, (H - 40) // rows
    for k, p in enumerate(items):
        pim = Image.open(p).convert("RGB")
        ps = min((cw - 16) / pim.width, (ch - 46) / pim.height, 1.0)
        th = pim.resize((max(1, int(pim.width * ps)), max(1, int(pim.height * ps))), Image.LANCZOS)
        cx, cy = LW + 10 + (k % cols) * cw, 40 + (k // cols) * ch
        lab = "legend" if (legend and k == n - 1) else f"{k + 1:02d}"
        c = (0, 0, 0) if lab == "legend" else PALETTE[k % len(PALETTE)]
        d.text((cx + 4, cy), f"{lab}  {pim.width}x{pim.height}", fill=c, font=small)
        sheet.paste(th, (cx + 8, cy + 30))
        d.rectangle([cx + 7, cy + 29, cx + 8 + th.width, cy + 30 + th.height], outline=c, width=2)
    out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out)


# ----------------------------------------------------------------------------------------------- main
def parent_tag(rel: str) -> str:
    return Path(rel).parent.as_posix().replace("/", "_")


def main():
    RES.mkdir(parents=True, exist_ok=True)
    CHECK.mkdir(parents=True, exist_ok=True)
    for p in list(PANELS.glob("*.png")) + list(CHECK.glob("*.png")) + list(CHECK.glob("_failed/*.png")):  # earlier VIS1 runs
        p.unlink()
    inv = inventory()
    with open(RES / "inventory.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(inv[0].keys()))
        w.writeheader()
        w.writerows(inv)
    groups = {}
    for r in inv:
        groups.setdefault(r["sha1"], []).append(r)

    manifest, log, box_rows = [], [], []
    for sha, members in sorted(groups.items(), key=lambda kv: kv[1][0]["path"]):
        rep = members[0]
        rel = rep["path"]
        paths = ";".join(m["path"] for m in members)
        base = dict(sha1=sha)
        L = dict(path=rel, all_paths=paths, sha1=sha, width=rep["width"], height=rep["height"],
                 svg_twin=rep["svg_twin"], cls="", mode="", expected_n="", n_panels=0, n_dropped=0,
                 status="", check=spec.CHECK.get(rel, "not checked"), foreign_ink="", note="")
        if rep["small"]:
            L.update(cls="skipped", status="skipped_small", note="< 200 px on a side")
            log.append(L)
            continue
        src = ROOT / rel
        stem = src.stem
        # ---- (a) SV figures
        cfg = None
        if rel.startswith("figures/sv/"):
            sv_id = rel.split("/")[2]
            pan = sorted((ROOT / "figures/panels/sv" / sv_id).glob(f"{stem}__*.png"))
            pan_rel = [q.relative_to(ROOT).as_posix() for q in pan]
            edges = {q.relative_to(ROOT).as_posix(): edge_ink(q) for q in pan}
            edge_note = "; ".join(f"{k.rsplit('/', 1)[1]}: ink at the {'/'.join(v)} edge" for k, v in edges.items()
                                  if v)
            L["note"] = ("edge check: " + edge_note) if edge_note else ""
            if len(pan) > 1:
                sheet_name = f"{parent_tag(rel)}__{stem}" + ("__sv_panels" if rel in spec.SV_RECUT else "")
                contact_sheet(src, pan, [], [], CHECK / f"{sheet_name}.png")
            if rel in spec.SV_RECUT:
                cfg = dict(spec.SV_RECUT[rel])
                cfg["sv_note"] = (f"SV-script panels ({', '.join(q.rsplit('/', 1)[1] for q in pan_rel)}) fail the "
                                  f"VIS1 check: {cfg['reason']}; panels cut from the original PNG instead")
            elif not pan or rel in spec.SV_SINGLE or len(pan) == 1:
                reason = spec.SV_SINGLE.get(rel, "one panel saved by the SV script" if pan else
                                            "no panel saved by the SV script")
                L.update(cls="single", mode="single", status="single", note="; ".join(
                    x for x in (reason, L["note"]) if x))
                manifest.append({"source_path(s)": paths, **base, "panel_path": "", "method": "single",
                                 "panel_title": "", "needs_legend": "",
                                 "notes": f"one graph ({reason}); original used as is"
                                          + (f"; SV-script panel(s) {', '.join(pan_rel)}" if pan else "")
                                          + (f"; {L['note']}" if L["note"] else "")})
                log.append(L)
                continue
            else:
                L.update(cls="multi", mode="sv_script", n_panels=len(pan), status="sv_listed", expected_n=len(pan))
                for q in pan_rel:
                    manifest.append({"source_path(s)": paths, **base, "panel_path": q, "method": "sv_script",
                                     "panel_title": q.rsplit("/", 1)[1][:-4].split("__", 1)[1], "needs_legend": "",
                                     "notes": f"saved per graph by the {sv_id} figure script (rule 1.8); visual "
                                              f"check {spec.CHECK.get(rel, 'not checked')}"})
                log.append(L)
                continue
        # ---- single
        if cfg is None and rel in spec.SINGLE:
            L.update(cls="single", mode="single", status="single", note=spec.SINGLE[rel])
            manifest.append({"source_path(s)": paths, **base, "panel_path": "", "method": "single",
                             "panel_title": "", "needs_legend": "",
                             "notes": f"one graph ({spec.SINGLE[rel]}); original used as is"})
            log.append(L)
            continue
        if cfg is None and rel not in spec.MULTI:
            L.update(cls="unclassified", status="not_cut", note="not classified in VIS1_spec.py")
            manifest.append({"source_path(s)": paths, **base, "panel_path": "", "method": "not_cut",
                             "panel_title": "", "needs_legend": "", "notes": L["note"]})
            log.append(L)
            continue
        cfg = cfg or spec.MULTI[rel]
        L.update(cls="multi", mode=cfg["mode"], expected_n=cfg["n"])
        a = load_rgb(src)
        blank, bgs = blank_mask(a)
        H, W = blank.shape
        leaves, drops, fails = [], [], []
        legend_leaf = None
        if cfg["mode"] == "svg":
            axes, ftexts, flegs, _ = svg_layout(src.with_suffix(".svg"), W, H)
            items = group_axes(axes, cfg.get("groups", ()))
            for fl in flegs:
                items.append(dict(kind="legend", ids=[fl["id"]], box=fl["box"], title="shared legend"))
            guided_cut(blank, (0, 0, W, H), items, leaves, drops, fails)
            graph_leaves = [lf for lf in leaves if lf["item"]["kind"] == "graph"]
            leg = [lf for lf in leaves if lf["item"]["kind"] == "legend"]
            legend_leaf = leg[0] if leg else None
        else:
            region, title_drop = (0, 0, W, H), []
            t = trim(blank, region)
            hr = [r for r in runs(blank[t[1]:t[3], t[0]:t[2]].all(1)) if r[1] - r[0] >= 3]
            if cfg.get("top_title", True) and hr and hr[0][0] <= 0.1 * H:
                # a suptitle strip above the first horizontal gutter is text only (not a graph)
                title_drop = [(t[0], t[1], t[2], t[1] + hr[0][0])]
                region = (0, t[1] + (hr[0][0] + hr[0][1]) // 2, W, H)
            bottom = [r for r in hr if (t[3] - t[1]) - r[1] <= 0.2 * H]
            if cfg.get("legend") == "bottom" and bottom:
                # a figure-level key below every graph (rule 5: saved once as <stem>__legend.png); the widest
                # horizontal gutter in the bottom 20 % separates it from the graphs
                r = max(bottom, key=lambda r: r[1] - r[0])
                ymid = t[1] + (r[0] + r[1]) // 2
                legend_leaf = dict(ink=trim(blank, (0, ymid, W, H)), allow=(0, ymid, W, H), item=None)
                region = (region[0], region[1], W, ymid)
            regions = [region]
            for xc in cfg.get("cut_x", []):
                # declared cut position (from the column gutter list); it must be background over the region
                new_regions = []
                for q in regions:
                    if q[0] < xc < q[2]:
                        if not blank[q[1]:q[3], xc].all():
                            raise SystemExit(f"{rel}: cut_x {xc} is not a background column")
                        new_regions += [(q[0], q[1], xc, q[3]), (xc, q[1], q[2], q[3])]
                    else:
                        new_regions.append(q)
                regions = new_regions
            L_run = max(60, int(0.08 * min(H, W)))

            def big_fn(q):
                """Graph-sized part: >= 12 % of the image height, >= 8 % of its width, and holding plot structure
                (a straight run of >= L_run ink pixels in some row or column: spine, bar, heatmap, render),
                which text-only strips (tick labels, titles) lack."""
                if (q[3] - q[1]) < 0.12 * H or (q[2] - q[0]) < 0.08 * W:
                    return False
                ink = ~blank[q[1]:q[3], q[0]:q[2]]
                return longest_run(ink) >= L_run
            for q in regions:
                gutter_cut(blank, q, cfg.get("g", 3), leaves, big_fn, cfg.get("only", "hv"))
            big, small = [], []
            for lf in leaves:
                (big if big_fn(lf["ink"]) else small).append(lf)
            graph_leaves = big
            drops = title_drop + [lf["ink"] for lf in small]
        graph_leaves = reading_order(graph_leaves)
        L["n_dropped"] = len(drops)
        problems = []
        if fails:
            msgs = []
            for f_ in fails:
                its = [it for it in items if it["ids"] in f_["ids"]]
                ov = ""
                if len(its) == 2:
                    a_, b_ = its[0]["box"], its[1]["box"]
                    ox = min(a_[2], b_[2]) - max(a_[0], b_[0])
                    oy = min(a_[3], b_[3]) - max(a_[1], b_[1])
                    ov = (f" (their extents from the SVG layout overlap by {ox:.0f} px in x and {oy:.0f} px in y:"
                          " a tick label or title of one reaches past the edge of the other)")
                msgs.append("+".join("/".join(i) for i in f_["ids"]) + ov)
            problems.append("no background row or column separates " + "; ".join(msgs)
                            + "; a rectangular PNG crop would cut one graph or include a fragment of the other "
                              "(method b, which hides the other axes, is unavailable)")
        if len(graph_leaves) != cfg["n"]:
            problems.append(f"{len(graph_leaves)} graphs found, {cfg['n']} expected")
        chk = spec.CHECK.get(rel, "not checked")
        if chk.startswith("fail"):
            problems.append("visual check " + chk)
        boxes = [crop_box(lf) for lf in graph_leaves]
        if cfg["mode"] == "svg":
            fr = []
            for k, (lf, b) in enumerate(zip(graph_leaves, boxes)):
                others = [it["box"] for it in items if it is not lf["item"]]
                others += [t["box"] for t in ftexts if t.get("box")]
                n_f = foreign_ink(blank, b, lf["item"]["box"], others)
                if n_f:
                    fr.append(f"{k + 1:02d}:{n_f}")
            L["foreign_ink"] = " ".join(fr)
        # write panels (also for failed figures, to the check folder only, so the failure can be seen)
        tag = f"{parent_tag(rel)}__{stem}"
        im = Image.open(src)
        out_paths, titles = [], []
        given = cfg.get("titles")
        for k, (lf, b) in enumerate(zip(graph_leaves, boxes)):
            title = lf["item"]["title"] if lf["item"] else (given[k] if given and k < len(given) else "panel")
            slug = slugify(title)
            titles.append(title)
            out_paths.append(PANELS / f"{tag}__{k + 1:02d}_{slug}.png")
        leg_path = PANELS / f"{tag}__legend.png" if legend_leaf else None
        if problems:
            tmp = CHECK / "_failed"
            tmp.mkdir(exist_ok=True)
            out_paths = [tmp / p.name for p in out_paths]
            leg_path = tmp / leg_path.name if leg_path else None
        for p, b in zip(out_paths, boxes):
            im.crop(b).save(p)
            box_rows.append(dict(source=rel, panel=p.relative_to(ROOT).as_posix(), x0=b[0], y0=b[1], x1=b[2],
                                 y1=b[3], status="failed" if problems else "panel"))
        if leg_path:
            lb = crop_box(legend_leaf)
            im.crop(lb).save(leg_path)
            box_rows.append(dict(source=rel, panel=leg_path.relative_to(ROOT).as_posix(), x0=lb[0], y0=lb[1],
                                 x1=lb[2], y1=lb[3], status="failed" if problems else "legend"))
        for d in drops:
            box_rows.append(dict(source=rel, panel="", x0=d[0], y0=d[1], x1=d[2], y1=d[3], status="dropped"))
        contact_sheet(src, out_paths, boxes, drops, CHECK / f"{tag}.png", legend=leg_path)
        if problems:
            if cfg.get("sv_note"):
                problems.append(cfg["sv_note"].replace("panels cut from the original PNG instead",
                                                       "the PNG cut fails too"))
            L.update(status="not_cut", note="; ".join(problems))
            manifest.append({"source_path(s)": paths, **base, "panel_path": "", "method": "not_cut",
                             "panel_title": "", "needs_legend": "",
                             "notes": "not cut: " + "; ".join(problems) + "; original used as is"})
            log.append(L)
            continue
        L.update(status="cut", n_panels=len(out_paths), note=cfg.get("note", ""))
        method_note = ("cut positions from the SVG axes layout; crop of the original PNG" if cfg["mode"] == "svg"
                       else "XY-cut along background gutters (PNG only)")
        legend_for = cfg.get("legend_for", "all" if leg_path else "")
        for k, (p, t) in enumerate(zip(out_paths, titles)):
            nl = ""
            if leg_path and (legend_for == "all" or (k + 1) in legend_for):
                nl = leg_path.relative_to(ROOT).as_posix()
            extra = cfg.get("panel_notes", {}).get(k + 1, "")
            manifest.append({"source_path(s)": paths, **base, "panel_path": p.relative_to(ROOT).as_posix(),
                             "method": "png_gutter", "panel_title": t, "needs_legend": nl,
                             "notes": "; ".join(x for x in (method_note, cfg.get("note", ""), extra,
                                                            cfg.get("sv_note", "")) if x)})
        if leg_path:
            manifest.append({"source_path(s)": paths, **base, "panel_path": leg_path.relative_to(ROOT).as_posix(),
                             "method": "png_gutter", "panel_title": "shared legend / colour key",
                             "needs_legend": "",
                             "notes": "figure-level key below the graphs, not next to one graph; saved once (rule 5)"})
        log.append(L)

    with open(PANELS / "manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source_path(s)", "sha1", "panel_path", "method", "panel_title",
                                          "needs_legend", "notes"])
        w.writeheader()
        w.writerows(manifest)
    with open(RES / "panel_boxes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "panel", "x0", "y0", "x1", "y1", "status"])
        w.writeheader()
        w.writerows(box_rows)
    with open(RES / "cut_log.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(log[0].keys()))
        w.writeheader()
        w.writerows(log)
    for L in log:
        if L["status"] in ("cut", "not_cut") and L["mode"] != "sv_script":
            print(f"{L['status']:8s} {L['n_panels']:3d}/{L['expected_n']!s:3s} drop={L['n_dropped']} {L['path']} "
                  f"foreign=[{L['foreign_ink']}] {L['note'][:160]}")
    print(f"images {len(inv)}, unique {len(groups)}, manifest rows {len(manifest)}")


if __name__ == "__main__":
    main()
