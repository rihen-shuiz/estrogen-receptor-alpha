#!/usr/bin/env python3
"""tl.py - task/log tracker (task language v1 + log language v1). Stdlib only.

Layout: this file lives in <repo>/tracker/ next to *.tasks and *.log files.
Each X.tasks pairs with X.log (use one pair per subsystem on big projects).
Checks run relative to the repo root = parent of tracker/ (override: TL_ROOT).

  tl.py validate                      parse everything, report errors/warnings
  tl.py status [-p PREFIX] [-n 10]    progress board, ready tasks, blocked tasks
  tl.py next [-p PREFIX]              print the next ready task card
  tl.py show ID                       task card + history
  tl.py check ID|all [--cmd]          run chk lines (no logging)
  tl.py audit [--write] [--cmd]       check every TODO task; --write logs DONE/PART as who=audit
  tl.py log ID STATE [--note ..] [--ev a,b] [--who ..] [--ref ..]
  tl.py log - NOTE --note "..." [--file main]
  tl.py bundle [--since ISO|last|all] [--ids A,B] [--all] [--sample N] [--budget CHARS]
               [--no-code] [--out FILE]        digest zip for a reviewer (no data dumps)
  tl.py pull PATH [--lines a:b] [--grep RX -C n] [--sym NAME] [--rows a:b] [--where col=v] [--cols a,b]
                                      exact bounded piece of one file, for the reviewer
  tl.py map [PATH ...]                code map: files, classes, functions
Env: TL_ROOT, TL_DIR, TL_TASKS (comma list), TL_WHO (default who), TL_SKIP (extra dirs to skip),
     TL_CMD_TIMEOUT (seconds, default 300)
"""
import argparse, ast, collections, contextlib, csv, datetime as dt, glob, hashlib, io, json, os
import platform, random, re, subprocess, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
TDIR = os.path.abspath(os.environ.get("TL_DIR", HERE))
ROOT = os.path.abspath(os.environ.get("TL_ROOT") or os.environ.get("ESR_ROOT") or os.path.dirname(HERE))
PROV = os.path.join(TDIR, "prov")
LAST_BUNDLE = os.path.join(TDIR, ".last_bundle")

SIZES = {"XS", "S", "M"}
KINDS = {"code", "run", "fig", "text", "decide", "manual"}
STATES = {"WIP", "DONE", "PART", "BLOCK", "MANUAL", "SKIP", "REOPEN"}
NOTE_REQ = {"PART", "BLOCK", "MANUAL", "SKIP", "WIP", "REOPEN"}
LOG_KEYS = {"who", "chk", "ev", "note", "ref"}
CHK_TYPES = {"file", "glob", "lines", "csv", "npy", "grep", "cmd", "man"}
SINGLE = {"file", "lines", "csv", "npy"}          # checks that read ONE file
SHOW_STATES = ["DONE", "PART", "WIP", "BLOCK", "MANUAL", "SKIP", "TODO"]
HASH_MAX = 256 * 1024 * 1024
SENT = {"NA", "N/A", "nan", "NaN", "None", "null", "NULL"}


def task_files():
    env = os.environ.get("TL_TASKS") or os.environ.get("ESR_TASKS")
    if env: return [p.strip() for p in env.split(",") if p.strip()]
    return sorted(glob.glob(os.path.join(TDIR, "*.tasks")))


def log_for(taskfile):
    return os.path.splitext(taskfile)[0] + ".log"


def bn(p): return os.path.basename(p)


# ---------- task language ----------
def parse_tasks(files):
    tasks, errs, conf = {}, [], {"phases": [], "who": []}
    for path in files:
        cur = None
        for n, raw in enumerate(open(path, encoding="utf-8"), 1):
            line = raw.rstrip("\n"); s = line.strip(); at = f"{bn(path)}:{n}"
            if not s or s.startswith("#") or s.startswith("LANG "): continue
            if line.startswith("CONF "):
                for kv in line[5:].split():
                    k, _, v = kv.partition("=")
                    if k in conf: conf[k] += [x for x in v.split(",") if x and x not in conf[k]]
                    else: errs.append(f"{at}: unknown CONF key '{k}'")
                continue
            if line.startswith("T "):
                m = re.match(r'^T (\S+) (\S+) (\S+) "(.*)"$', line)
                if not m: errs.append(f"{at}: bad T line"); cur = None; continue
                tid, size, kind, title = m.groups()
                if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", tid): errs.append(f"{at}: bad id '{tid}'")
                if tid in tasks: errs.append(f"{at}: duplicate id {tid} (also in {bn(tasks[tid]['file'])})")
                cur = dict(id=tid, size=size, kind=kind, title=title, help="", meta={}, needs=[], tool=[],
                           out=[], chk=[], orig=[], src="", line=n, file=path, raw=[line])
                tasks[tid] = cur
                continue
            if not line.startswith("  ") or cur is None:
                errs.append(f"{at}: field outside a task"); continue
            cur["raw"].append(line)
            key, _, rest = s.partition(" "); rest = rest.strip()
            if key == "help": cur["help"] = rest
            elif key == "meta":
                for kv in rest.split():
                    k, _, v = kv.partition("="); cur["meta"][k] = v
            elif key == "needs": cur["needs"] += [x for x in re.split(r"[ ,]+", rest) if x]
            elif key == "tool": cur["tool"] += [x.strip() for x in rest.split(",") if x.strip()]
            elif key == "out": cur["out"].append(rest)
            elif key == "chk":
                t, _, a = rest.partition(" "); cur["chk"].append((t, a.strip()))
            elif key == "orig": cur["orig"].append(rest)
            elif key == "src": cur["src"] = rest
            else: errs.append(f"{at}: unknown field '{key}'")
    for t in tasks.values():
        p = f"{bn(t['file'])}:{t['line']}: {t['id']}"
        m = t["meta"]
        m.setdefault("wk", "0"); m.setdefault("ph", "-"); m.setdefault("who", "-")
        if t["size"] not in SIZES: errs.append(f"{p}: bad size {t['size']}")
        if t["kind"] not in KINDS: errs.append(f"{p}: bad kind {t['kind']}")
        if not m["wk"].isdigit(): errs.append(f"{p}: meta wk not an integer")
        if conf["who"] and m["who"] not in conf["who"] + ["ALL", "-"]: errs.append(f"{p}: meta who '{m['who']}' not in CONF who")
        if conf["phases"] and m["ph"] not in conf["phases"] + ["-"]: errs.append(f"{p}: meta ph '{m['ph']}' not in CONF phases")
        if not t["help"]: errs.append(f"{p}: missing help")
        if not t["chk"]: errs.append(f"{p}: no chk")
        for x in t["needs"]:
            if x not in tasks: errs.append(f"{p}: needs unknown {x}")
        for ty, ar in t["chk"]:
            if ty not in CHK_TYPES: errs.append(f"{p}: bad chk type {ty}")
            elif ty in SINGLE and re.search(r"[*?\[]", ar.split(" ")[0]):
                errs.append(f"{p}: warn: wildcard in single-file chk '{ty} {ar.split(' ')[0]}' (first match wins; use the exact path)")
            elif ty == "cmd" and "tl.py" in ar:
                errs.append(f"{p}: warn: cmd chk calls tl.py (circular)")
        if t["kind"] in ("manual", "decide") and not any(ty == "man" for ty, _ in t["chk"]):
            errs.append(f"{p}: manual/decide task needs a 'chk man'")
        if t["kind"] in ("code", "run", "fig", "text") and not t["out"]:
            errs.append(f"{p}: warn: no out path")
        if len(t["chk"]) == 1 and t["chk"][0][0] == "file" and t["kind"] in ("code", "run"):
            errs.append(f"{p}: warn: only a 'file' check on a {t['kind']} task (existence says little; add csv/npy/grep/cmd)")
    color = {}
    def dfs(u, stack):
        color[u] = 1
        for v in tasks[u]["needs"]:
            if v not in tasks: continue
            if color.get(v) == 1: errs.append(f"cycle: {' > '.join(stack + [u, v])}")
            elif v not in color: dfs(v, stack + [u])
        color[u] = 2
    for u in tasks:
        if u not in color: dfs(u, [])
    return tasks, errs, conf


# ---------- checks ----------
def matches(pat):
    return sorted(glob.glob(os.path.join(ROOT, pat), recursive=True))


def first_match(pat):
    hits = matches(pat)
    return (hits[0] if hits else None), len(hits)


def amb(k): return f" [{k} matches, used first]" if k > 1 else ""


def run_chk(ctype, a, skip_cmd=False):
    """-> (PASS|FAIL|MAN|SKIP, detail)"""
    try:
        if ctype == "man": return "MAN", a
        if ctype == "file":
            p, k = first_match(a)
            return ("PASS" if p and os.path.getsize(p) > 0 else "FAIL"), (os.path.relpath(p, ROOT) + amb(k) if p else "missing")
        if ctype == "glob":
            m = re.match(r"^(.*\S)\s+>=(\d+)$", a)
            pat, n = (m.group(1), int(m.group(2))) if m else (a, 1)
            k = len(matches(pat))
            return ("PASS" if k >= n else "FAIL"), f"{k} match (need {n})"
        if ctype == "lines":
            path, _, rest = a.partition(" ")
            n = int(re.search(r">=(\d+)", rest).group(1))
            p, k = first_match(path)
            if not p: return "FAIL", "missing"
            c = sum(1 for l in open(p, encoding="utf-8", errors="replace") if l.strip())
            return ("PASS" if c >= n else "FAIL"), f"{c} lines (need {n})" + amb(k)
        if ctype == "csv":
            toks = a.split()
            p, k = first_match(toks[0])
            if not p: return "FAIL", "missing"
            need_rows, need_cols, ne, uq = 0, [], [], []
            for t in toks[1:]:
                if t.startswith("rows>="): need_rows = int(t[6:])
                elif t.startswith("cols="): need_cols = t[5:].split(",")
                elif t.startswith("nonempty="): ne = t[9:].split(",")
                elif t.startswith("uniq="): uq = t[5:].split(",")
                else: return "FAIL", f"unknown csv option '{t}'"
            rows, empty, dup, seen = 0, collections.Counter(), collections.Counter(), {c: set() for c in uq}
            with open(p, encoding="utf-8-sig", newline="") as f:
                rd = csv.DictReader(f); cols = rd.fieldnames or []
                for r in rd:
                    rows += 1
                    for c in ne:
                        if (r.get(c) or "").strip() == "": empty[c] += 1
                    for c in uq:
                        v = r.get(c)
                        if v in seen[c]: dup[c] += 1
                        seen[c].add(v)
            miss = [c for c in need_cols + ne + uq if c not in cols]
            ok = rows >= need_rows and not miss and not empty and not dup
            d = f"{rows} rows (need {need_rows})"
            if miss: d += f", missing cols {sorted(set(miss))}"
            if empty: d += ", empty cells " + ",".join(f"{c}={n}" for c, n in empty.items())
            if dup: d += ", duplicates " + ",".join(f"{c}={n}" for c, n in dup.items())
            return ("PASS" if ok else "FAIL"), d + amb(k)
        if ctype == "npy":
            toks = a.split(); path = toks[0]
            n = int(re.search(r">=(\d+)", a).group(1)) if ">=" in a else 1
            p, k = first_match(path)
            if not p: return "FAIL", "missing"
            try: import numpy as np
            except ImportError: return "FAIL", "numpy not installed"
            x = np.load(p, mmap_mode="r"); sh = x.shape
            ok = len(sh) == 2 and sh[0] == sh[1] and sh[0] >= n
            d = f"shape {sh} (need square >= {n})"
            if "nonan" in toks and x.dtype.kind == "f":
                nn = int(np.isnan(x).sum()); ok = ok and nn == 0; d += f", nan={nn}"
            return ("PASS" if ok else "FAIL"), d + amb(k)
        if ctype == "grep":
            pat, _, rx = a.partition(" ")
            hits = [h for h in matches(pat) if os.path.isfile(h)]
            for p in hits:
                if re.search(rx, open(p, encoding="utf-8", errors="replace").read(), re.M):
                    return "PASS", os.path.relpath(p, ROOT)
            return "FAIL", f"no match in {len(hits)} file(s)"
        if ctype == "cmd":
            if skip_cmd: return "SKIP", "cmd not run (use --cmd)"
            r = subprocess.run(a, shell=True, cwd=ROOT, capture_output=True, text=True,
                               timeout=int(os.environ.get("TL_CMD_TIMEOUT", "300")))
            tail = (r.stderr.strip() or r.stdout.strip()).splitlines()[-1:] or [""]
            return ("PASS" if r.returncode == 0 else "FAIL"), f"exit {r.returncode} {tail[0][:100]}"
    except Exception as e:
        return "FAIL", f"{type(e).__name__}: {e}"
    return "FAIL", f"unknown chk type {ctype}"


def evaluate(task, skip_cmd=False):
    return [(t, a) + run_chk(t, a, skip_cmd) for t, a in task["chk"]]


# ---------- log language ----------
def parse_logs(tasks, files):
    ents, errs = [], []
    for path in files:
        lf = log_for(path)
        if not os.path.exists(lf): continue
        prev = None
        for n, raw in enumerate(open(lf, encoding="utf-8"), 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#") or line.startswith("LANG "): continue
            parts = line.split(" | "); p = f"{bn(lf)}:{n}"
            if len(parts) < 3: errs.append(f"{p}: need 'time | id | STATE | ...'"); continue
            ts, tid, state = parts[0].strip(), parts[1].strip(), parts[2].strip()
            try:
                t = dt.datetime.fromisoformat(ts)
                if t.tzinfo is None: t = t.astimezone()
                if prev and t < prev: errs.append(f"{p}: warn: time goes backwards")
                prev = t
            except ValueError:
                errs.append(f"{p}: bad time '{ts}'"); continue
            kv = {}
            for x in parts[3:]:
                k, _, v = x.partition("=")
                if k.strip() not in LOG_KEYS: errs.append(f"{p}: unknown key '{k.strip()}'"); continue
                kv[k.strip()] = v.strip()
            if state == "NOTE":
                if tid != "-": errs.append(f"{p}: NOTE needs id '-'")
                if not kv.get("note"): errs.append(f"{p}: NOTE needs note=")
            else:
                if state not in STATES: errs.append(f"{p}: bad state {state}")
                if tid not in tasks: errs.append(f"{p}: unknown task {tid}")
                if not kv.get("who"): errs.append(f"{p}: missing who=")
                if state == "DONE":
                    m = re.match(r"^(\d+)/(\d+)$", kv.get("chk", ""))
                    if not m or m.group(1) != m.group(2): errs.append(f"{p}: DONE needs chk=n/n all passing")
                    if not kv.get("ev"): errs.append(f"{p}: DONE needs ev=")
                if state in NOTE_REQ and not kv.get("note"):
                    errs.append(f"{p}: {'warn: ' if state == 'REOPEN' else ''}{state} needs note=")
            ents.append(dict(ts=ts, t=t, id=tid, state=state, kv=kv, file=lf, line=n))
    ents.sort(key=lambda e: e["t"])        # stable: same-minute lines keep file order
    return ents, errs


def current_states(tasks, ents):
    st = {t: "TODO" for t in tasks}
    for e in ents:
        if e["id"] in st:
            st[e["id"]] = "TODO" if e["state"] == "REOPEN" else e["state"]
    return st


def last_lines(tasks, ents):
    last = {}
    for e in ents:
        if e["id"] in tasks: last[e["id"]] = e
    return last


def append_log(logfile, tid, state, who, chk=None, ev=None, note=None, ref=None):
    clean = lambda s: re.sub(r"\s+", " ", str(s).replace(" | ", " / ")).strip()
    ts = dt.datetime.now().astimezone().isoformat(timespec="minutes")
    parts = [ts, tid, state, f"who={who}"]
    for k, v in (("chk", chk), ("ev", ev), ("ref", ref), ("note", note)):
        if v: parts.append(f"{k}={clean(v)}")
    new = not os.path.exists(logfile)
    with open(logfile, "a", encoding="utf-8") as f:
        if new: f.write("LANG log v1\n")
        f.write(" | ".join(parts) + "\n")
    return " | ".join(parts)


# ---------- provenance: hash of the evidence at the moment of DONE ----------
def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()[:16]


def sh(cmd, timeout=60):
    try:
        r = subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        return (r.stdout + r.stderr).strip()
    except Exception as e:
        return f"({type(e).__name__})"


def in_git(): return sh("git rev-parse --is-inside-work-tree") == "true"


def write_prov(tid, who, ev):
    items = []
    for p in ev:
        fp = os.path.join(ROOT, p)
        if os.path.isfile(fp):
            sz = os.path.getsize(fp)
            items.append(dict(path=p, size=sz, sha=sha(fp) if sz <= HASH_MAX else None))
    d = dict(task=tid, ts=dt.datetime.now().astimezone().isoformat(timespec="minutes"), who=who,
             git=(sh("git rev-parse --short HEAD 2>/dev/null") or None) if in_git() else None, tl_sha=sha(os.path.abspath(__file__)), ev=items)
    os.makedirs(PROV, exist_ok=True)
    with open(os.path.join(PROV, f"{tid}.json"), "w", encoding="utf-8") as f: json.dump(d, f, indent=1)


def prov_drift(tid):
    """-> (status, [lines]); status in none|ok|drift"""
    pf = os.path.join(PROV, f"{tid}.json")
    if not os.path.exists(pf): return "none", []
    d = json.load(open(pf, encoding="utf-8")); bad = []
    for it in d["ev"]:
        fp = os.path.join(ROOT, it["path"])
        if not os.path.isfile(fp): bad.append(f"{it['path']}: missing"); continue
        if it["sha"] and sha(fp) != it["sha"]: bad.append(f"{it['path']}: changed after DONE ({d['ts']})")
    return ("drift" if bad else "ok"), bad


# ---------- commands ----------
def load():
    files = task_files()
    if not files: return {}, [], ["no *.tasks file found in " + TDIR], [], files, {"phases": [], "who": []}
    tasks, terr, conf = parse_tasks(files)
    ents, lerr = parse_logs(tasks, files)
    return tasks, ents, terr, lerr, files, conf


def hard(errs): return [e for e in errs if "warn:" not in e]


def cmd_validate(_):
    tasks, ents, terr, lerr, files, _c = load()
    errs = terr + lerr
    for e in errs: print(e)
    print(f"{len(tasks)} tasks in {len(files)} file(s), {len(ents)} log lines, {len(hard(errs))} errors, {len(errs) - len(hard(errs))} warnings")
    return 1 if hard(errs) else 0


def ready_tasks(tasks, st, prefix=None):
    r = [t for t in tasks.values() if st[t["id"]] in ("TODO", "PART", "WIP")
         and all(st[x] in ("DONE", "SKIP") for x in t["needs"]) and (not prefix or t["id"].startswith(prefix))]
    r.sort(key=lambda t: (int(t["meta"]["wk"]), t["id"]))
    return r


def cmd_status(a):
    tasks, ents, terr, lerr, files, conf = load()
    if hard(terr): print("TASK FILE ERRORS - run validate"); return 1
    st = current_states(tasks, ents)
    phases = list(conf["phases"]) + [p for p in dict.fromkeys(t["meta"]["ph"] for t in tasks.values()) if p not in conf["phases"]]
    w = max([len(p) for p in phases] + [5]) + 1
    print(f"{'phase':{w}}{'total':>6}" + "".join(f"{c:>8}" for c in SHOW_STATES))
    for ph in phases:
        ids = [t for t in tasks if tasks[t]["meta"]["ph"] == ph]
        if not ids: continue
        print(f"{ph:{w}}{len(ids):>6}" + "".join(f"{sum(st[i] == c for i in ids):>8}" for c in SHOW_STATES))
    done = sum(v == "DONE" for v in st.values())
    print(f"TOTAL {done}/{len(tasks)} DONE ({100 * done // max(1, len(tasks))}%)")
    ready = ready_tasks(tasks, st, a.prefix)
    auto = [t for t in ready if t["kind"] not in ("manual", "decide")]
    human = [t for t in ready if t["kind"] in ("manual", "decide")]
    print(f"\nREADY for Claude (top {a.n}):")
    for t in auto[:a.n]: print(f"  {t['id']:6} {t['size']:2} {t['kind']:6} wk{t['meta']['wk']:>2} {st[t['id']]:5} {t['title']}")
    print("\nNEEDS HUMAN (manual/decide, ready):")
    for t in human[:a.n]: print(f"  {t['id']:6} {t['size']:2} {t['kind']:6} wk{t['meta']['wk']:>2} {t['title']}")
    last = last_lines(tasks, ents)
    blocked = [t for t in tasks if st[t] in ("BLOCK", "MANUAL")]
    if blocked:
        print("\nBLOCK/MANUAL:")
        for t in blocked: print(f"  {t:6} {st[t]:6} {last[t]['kv'].get('note', '')[:90]}")
    return 0


def cmd_next(a):
    tasks, ents, terr, lerr, files, conf = load()
    if hard(terr): print("TASK FILE ERRORS - run validate"); return 1
    st = current_states(tasks, ents)
    auto = [t for t in ready_tasks(tasks, st, a.prefix) if t["kind"] not in ("manual", "decide")]
    if not auto: print("nothing ready for Claude"); return 0
    t = auto[0]; print("\n".join(t["raw"])); print(f"  # state: {st[t['id']]}")
    return 0


def cmd_show(a):
    tasks, ents, *_ = load()
    if a.id not in tasks: print(f"unknown task {a.id}"); return 1
    print("\n".join(tasks[a.id]["raw"]))
    hist = [e for e in ents if e["id"] == a.id]
    if hist:
        print("  # history:")
        for e in hist: print(f"  #   {e['ts']} {e['state']} who={e['kv'].get('who', '')} {e['kv'].get('note', '')[:80]}")
    return 0


def print_results(res):
    for ty, ar, r, d in res: print(f"  {r:5} chk {ty} {ar}  [{d}]")


def cmd_check(a):
    tasks, *_ = load()
    ids = list(tasks) if a.id == "all" else [a.id]
    for i in ids:
        if i not in tasks: print(f"unknown task {i}"); return 1
        print(f"{i} {tasks[i]['title']}"); print_results(evaluate(tasks[i], not a.cmd))
    return 0


def derive(res):
    states = [r[2] for r in res]
    if "MAN" in states or "SKIP" in states: return None
    p = states.count("PASS")
    return "DONE" if p == len(states) else ("PART" if p else None)


def cmd_audit(a):
    tasks, ents, terr, lerr, files, _c = load()
    if hard(terr): print("TASK FILE ERRORS - run validate"); return 1
    st = current_states(tasks, ents)
    n_done = n_part = 0
    for t in tasks.values():
        if st[t["id"]] != "TODO": continue
        res = evaluate(t, not a.cmd); d = derive(res)
        passed = sum(r[2] == "PASS" for r in res)
        print(f"{t['id']:6} {d or '-':5} {passed}/{len(res)}  {t['title']}")
        if a.write and d:
            have = [os.path.relpath(m, ROOT) for o in t["out"] for m in matches(o)[:3] if os.path.isfile(m)]
            ev = ",".join(dict.fromkeys(have)) or "chk"
            lf = log_for(t["file"])
            if d == "DONE":
                append_log(lf, t["id"], "DONE", "audit", chk=f"{passed}/{len(res)}", ev=ev, note="audit: all checks pass (shape only)")
                write_prov(t["id"], "audit", [] if ev == "chk" else ev.split(",")); n_done += 1
            else:
                failed = "; ".join(f"{r[0]} {r[1]}" for r in res if r[2] == "FAIL")
                append_log(lf, t["id"], "PART", "audit", chk=f"{passed}/{len(res)}", ev=ev, note=f"audit: failing: {failed}"); n_part += 1
    if a.write: print(f"\nwrote {n_done} DONE, {n_part} PART")
    return 0


def cmd_log(a):
    tasks, ents, terr, lerr, files, _c = load()
    if hard(terr): print("TASK FILE ERRORS - run validate"); return 1
    who = a.who or os.environ.get("TL_WHO", "claude")
    if a.id == "-":
        if a.state != "NOTE" or not a.note: print("usage: log - NOTE --note TEXT [--file NAME]"); return 1
        tf = next((f for f in files if bn(f)[:-6] == a.file), files[0]) if a.file else files[0]
        print(append_log(log_for(tf), "-", "NOTE", who, note=a.note)); return 0
    if a.id not in tasks: print(f"unknown task {a.id}"); return 1
    if a.state not in STATES: print(f"state must be one of {sorted(STATES)}"); return 1
    t, st, chk = tasks[a.id], current_states(tasks, ents), None
    if a.state in NOTE_REQ and not a.note: print(f"{a.state} needs --note"); return 1
    ev = [x.strip() for x in (a.ev or "").split(",") if x.strip()]
    if a.state == "DONE":
        if not ev: print("DONE needs --ev (files produced or inspected, comma separated, repo-relative)"); return 1
        if any(ty == "man" for ty, _ in t["chk"]) and who != "human":
            print("task has a manual check: log MANUAL and let the human log DONE (--who human)"); return 1
        if not a.audit:
            open_needs = [x for x in t["needs"] if st[x] not in ("DONE", "SKIP")]
            if open_needs: print(f"needs not DONE: {open_needs}"); return 1
        badev = [p for p in ev if os.path.isabs(p) or p.startswith("..") or not os.path.isfile(os.path.join(ROOT, p))]
        if badev: print(f"REFUSED: ev must be existing repo-relative files: {badev}"); return 1
        res = evaluate(t); auto = [r for r in res if r[2] != "MAN"]
        print_results(res)
        if [r for r in auto if r[2] != "PASS"]:
            print("REFUSED: checks failing. Fix, or log PART/BLOCK with --note."); return 1
        if t["out"] and not any(os.path.relpath(m, ROOT) in ev for o in t["out"] for m in matches(o)):
            print("warn: ev contains none of the task's out files")
        chk = f"{len(auto)}/{len(auto)}"
    line = append_log(log_for(t["file"]), a.id, a.state, who, chk=chk, ev=a.ev, note=a.note, ref=a.ref)
    if a.state == "DONE": write_prov(a.id, who, ev)
    print(line)
    return 0


# ---------- digests: short text that lets a reviewer judge a file without having it ----------
SKIP_DIRS = {".git", ".venv", "venv", "env", "node_modules", "__pycache__", "build", "dist", ".mypy_cache",
             ".pytest_cache", ".ipynb_checkpoints", "site-packages"}
CODE_EXT = {".py", ".sh", ".r", ".sql", ".pml", ".js", ".ts", ".go", ".rs", ".java", ".c", ".cpp", ".h", ".jl", ".m"}
TEXT_EXT = CODE_EXT | {".md", ".txt", ".json", ".yaml", ".yml", ".toml", ".cfg", ".ini", ".ipynb", ".log", ".list", ".tasks"}
SECRET = re.compile(r"(?i)((api[_-]?key|token|secret|passw(or)?d)\s*[=:]\s*)\S+")


def redact(s): return SECRET.sub(r"\1[REDACTED]", s)


def skip_dirs():
    return SKIP_DIRS | {x for x in os.environ.get("TL_SKIP", "").split(",") if x}


def digest_csv(path, ext):
    delim = "\t" if ext == ".tsv" else ","
    with open(path, encoding="utf-8-sig", newline="") as f:
        rd = csv.reader(f, delimiter=delim); hdr = next(rd, None)
        if not hdr: return "empty file"
        k = len(hdr); n = ragged = 0
        empty, sent = [0] * k, [0] * k
        num = [[0, None, None, 0.0] for _ in range(k)]       # count, min, max, sum
        head, tail = [], collections.deque(maxlen=2)
        for r in rd:
            n += 1
            if len(r) != k: ragged += 1
            for i, v in enumerate(r[:k]):
                if v == "": empty[i] += 1; continue
                if v in SENT: sent[i] += 1
                try:
                    x = float(v); c = num[i]
                    c[0] += 1; c[1] = x if c[1] is None or x < c[1] else c[1]; c[2] = x if c[2] is None or x > c[2] else c[2]; c[3] += x
                except ValueError: pass
            if n <= 5: head.append(r)
            else: tail.append(r)
    out = [f"rows={n} cols={k}" + (f" RAGGED_ROWS={ragged}" if ragged else ""),
           "columns (empty cells): " + ", ".join(f"{c}({e})" for c, e in zip(hdr, empty))]
    s = [f"{hdr[i]}({sent[i]})" for i in range(k) if sent[i]]
    if s: out.append("literal NA-like strings (NA/nan/None/null) per column: " + ", ".join(s))
    nm = [f"{hdr[i]}[min={c[1]:.4g} max={c[2]:.4g} mean={c[3] / c[0]:.4g}]" for i, c in enumerate(num)
          if c[0] and c[0] == n - empty[i]][:15]
    if nm: out.append("numeric: " + " ".join(nm))
    out += [",".join(r)[:300] for r in head]
    if tail: out += ["..."] + [",".join(r)[:300] for r in tail]
    return "\n".join(out)


def digest_file(path, cap=8000):
    size = os.path.getsize(path)
    head = f"## {os.path.relpath(path, ROOT)}  {size} bytes  sha256:{sha(path)}\n"
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext in (".csv", ".tsv"): return head + digest_csv(path, ext) + "\n"
        if ext in (".npy", ".npz"):
            import numpy as np
            if ext == ".npz":
                z = np.load(path); return head + "arrays: " + ", ".join(f"{k}{z[k].shape}" for k in z.files) + "\n"
            x = np.load(path, mmap_mode="r"); out = [f"shape={x.shape} dtype={x.dtype}"]
            if x.size and x.size <= 50_000_000 and x.dtype.kind in "fiu":
                out.append(f"nan={int(np.isnan(x).sum()) if x.dtype.kind == 'f' else 0} min={np.nanmin(x):.6g} max={np.nanmax(x):.6g} mean={np.nanmean(x):.6g}")
                if x.ndim == 2 and x.shape[0] == x.shape[1]:
                    out.append(f"symmetric={bool(np.allclose(x, x.T, equal_nan=True))} max|diag|={float(np.abs(np.diag(x)).max()):.3g}")
            return head + "\n".join(out) + "\n"
        if ext == ".png":
            import struct
            with open(path, "rb") as f: b = f.read(24)
            if b[:8] == b"\x89PNG\r\n\x1a\n": return head + "png %dx%d px\n" % struct.unpack(">II", b[16:24])
            return head
        if ext in TEXT_EXT or ext == "":
            txt = open(path, encoding="utf-8", errors="replace").read()
            lines = txt.splitlines()
            if ext in (".txt", ".list", ".log") and len(lines) > 60:
                return head + f"{len(lines)} lines\n" + redact("\n".join(lines[:25] + ["..."] + lines[-10:]))[:cap] + "\n"
            if len(txt) > cap:
                return head + redact(txt[: int(cap * .75)]) + f"\n[... {len(txt) - cap} chars omitted ...]\n" + redact(txt[-int(cap * .25):]) + "\n"
            return head + redact(txt) + "\n"
    except Exception as e:
        return head + f"digest failed: {type(e).__name__}: {e}\n"
    return head


# ---------- code map ----------
def code_files():
    skip = skip_dirs()
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = sorted(d for d in dn if d not in skip and not d.startswith(".") and not d.endswith(".egg-info")
                       and os.path.abspath(os.path.join(dp, d)) != TDIR)
        for f in sorted(fn):
            if os.path.splitext(f)[1].lower() in CODE_EXT: yield os.path.join(dp, f)


def first_doc(node):
    d = ast.get_docstring(node)
    return (" - " + d.strip().splitlines()[0][:90]) if d else ""


def sig(fn):
    a = fn.args; names = [x.arg for x in a.posonlyargs + a.args]
    if a.vararg: names.append("*" + a.vararg.arg)
    names += [x.arg for x in a.kwonlyargs]
    if a.kwarg: names.append("**" + a.kwarg.arg)
    return ", ".join(names)


def file_map(path):
    rel = os.path.relpath(path, ROOT)
    try: src = open(path, encoding="utf-8", errors="replace").read()
    except Exception as e: return f"{rel}  (unreadable: {e})\n"
    n = src.count("\n") + 1
    if not path.endswith(".py"):
        pat = re.compile(r"^\s*(?:export\s+)?(?:async\s+)?(?:pub\s+)?(?:def|class|function|func|fn|struct|interface|impl)\s+\w+[^\n{;]*", re.M)
        items = [f"  :{src.count(chr(10), 0, m.start()) + 1} {m.group(0).strip()[:100]}" for m in pat.finditer(src)][:60]
        return f"{rel} ({n} lines)\n" + "\n".join(items) + ("\n" if items else "")
    try: tree = ast.parse(src)
    except SyntaxError as e: return f"{rel} ({n} lines) SYNTAX ERROR line {e.lineno}\n"
    out = [f"{rel} ({n} lines){first_doc(tree)}"]
    imps = sorted({(x.names[0].name if isinstance(x, ast.Import) else (x.module or "")).split(".")[0]
                   for x in ast.walk(tree) if isinstance(x, (ast.Import, ast.ImportFrom)) and (isinstance(x, ast.Import) or x.level == 0)} - {""})
    if imps: out.append("  imports: " + ", ".join(imps))
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append(f"  def {node.name}({sig(node)}) :{node.lineno}{first_doc(node)}")
        elif isinstance(node, ast.ClassDef):
            out.append(f"  class {node.name} :{node.lineno}{first_doc(node)}")
            for m in node.body:
                if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    out.append(f"    def {m.name}({sig(m)}) :{m.lineno}")
    return "\n".join(out) + "\n"


def cmd_map(a):
    paths = []
    for p in (a.paths or [ROOT]):
        fp = os.path.join(ROOT, p)
        if os.path.isdir(fp):
            paths += [x for x in code_files() if x.startswith(os.path.abspath(fp))]
        elif os.path.isfile(fp): paths.append(fp)
    for p in paths: print(file_map(p), end="")
    return 0


# ---------- pull: exact bounded piece of a file ----------
def parse_range(s, default_end=10 ** 9):
    a, _, b = s.partition(":")
    return int(a or 1), int(b) if b else default_end


def cmd_pull(a):
    p = os.path.realpath(os.path.join(ROOT, a.path))
    if not p.startswith(os.path.realpath(ROOT) + os.sep) or not os.path.isfile(p):
        print("file not found inside the repo"); return 1
    ext = os.path.splitext(p)[1].lower(); sel = ""
    if a.sym:
        src = open(p, encoding="utf-8", errors="replace").read(); tree = ast.parse(src); lines = src.splitlines()
        want = a.sym.split("."); node = tree
        for w in want:
            node = next((n for n in getattr(node, "body", []) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.name == w), None)
            if node is None: print(f"symbol {a.sym} not found"); return 1
        body = "\n".join(lines[node.lineno - 1: node.end_lineno]); sel = f"symbol {a.sym} lines {node.lineno}-{node.end_lineno}"
    elif ext in (".csv", ".tsv") and (a.rows or a.where or a.cols):
        r0, r1 = parse_range(a.rows or "1:") ; out = []; kept = n = 0
        wk, _, wv = (a.where or "").partition("=")
        with open(p, encoding="utf-8-sig", newline="") as f:
            rd = csv.DictReader(f, delimiter="\t" if ext == ".tsv" else ",")
            cols = a.cols.split(",") if a.cols else rd.fieldnames
            out.append(",".join(cols))
            for r in rd:
                n += 1
                if not (r0 <= n < r1): continue
                if a.where and (r.get(wk) or "") != wv: continue
                kept += 1; out.append(",".join((r.get(c) or "") for c in cols))
                if sum(len(x) for x in out) > a.cap: out.append("[cap reached]"); break
        body = "\n".join(out); sel = f"rows {a.rows or 'all'} where {a.where or '-'} -> {kept} rows"
    else:
        lines = open(p, encoding="utf-8", errors="replace").read().splitlines()
        if a.grep:
            rx = re.compile(a.grep); keep = set()
            for i, l in enumerate(lines):
                if rx.search(l): keep.update(range(max(0, i - a.C), min(len(lines), i + a.C + 1)))
            body = "\n".join(f"{i + 1}: {lines[i]}" for i in sorted(keep)); sel = f"grep {a.grep} ({len(keep)} lines)"
        else:
            l0, l1 = parse_range(a.lines or "1:40", 40)
            body = "\n".join(f"{i + 1}: {l}" for i, l in enumerate(lines[l0 - 1:l1], l0 - 1)); sel = f"lines {l0}-{min(l1, len(lines))} of {len(lines)}"
    body = redact(body)
    if len(body) > a.cap: body = body[:a.cap] + f"\n[truncated at {a.cap} chars]"
    print(f"## {os.path.relpath(p, ROOT)}  {sel}\n{body}")
    return 0


# ---------- review bundle ----------
def cmd_bundle(a):
    tasks, ents, terr, lerr, files_, conf = load()
    if hard(terr): print("TASK FILE ERRORS - run validate"); return 1
    st = current_states(tasks, ents); last = last_lines(tasks, ents)
    now = dt.datetime.now().astimezone()
    if a.since in (None, "all") or (a.since == "last" and not os.path.exists(LAST_BUNDLE)):
        since = ents[0]["t"] if ents else now          # default: whole history (queue stays small anyway)
    elif a.since == "last": since = dt.datetime.fromisoformat(open(LAST_BUNDLE).read().strip())
    else: since = dt.datetime.fromisoformat(a.since)
    if since.tzinfo is None: since = since.astimezone()
    budget = a.budget; used = [0]; files = {}

    def add(name, text):
        if used[0] + len(text) > budget:
            files.setdefault("TRUNCATED.txt", ""); files["TRUNCATED.txt"] += f"not included (size budget): {name}\n"; return
        used[0] += len(text); files[name] = text

    # ---- who needs review
    queue = collections.OrderedDict()
    for tid, e in sorted(last.items(), key=lambda kv: kv[1]["t"], reverse=True):
        if e["state"] in ("PART", "BLOCK", "MANUAL", "WIP", "REOPEN") and len(queue) < 40:
            queue[tid] = f"{e['state']} (always reviewed)"
    done_win = [t for t in tasks if st[t] == "DONE" and last[t]["t"] >= since]
    for t in done_win:
        if tasks[t]["kind"] in ("text", "decide") or a.all: queue.setdefault(t, "DONE text/decide (100% review)" if not a.all else "DONE (--all)")
    rest = sorted(t for t in done_win if t not in queue)
    rng = random.Random(f"{len(ents)}{now:%Y%m%d}")
    for t in rng.sample(rest, min(a.sample, len(rest))): queue[t] = "DONE sampled" + (" (audit: shape-only)" if last[t]["kv"].get("who") == "audit" else "")
    for t in (a.ids.split(",") if a.ids else []):
        if t in tasks: queue.setdefault(t, "requested")
    for tid, why in queue.items():
        e = last[tid] if tid in last else None
        if not e: continue
        paths = [p for p in e["kv"].get("ev", "").split(",") if p and p != "chk"]
        for o in tasks[tid]["out"]: paths += [os.path.relpath(m, ROOT) for m in matches(o)[:5]]
        seen, txt = set(), f"# {tid} {tasks[tid]['title']}\nreview reason: {why}\nstate={st[tid]} who={e['kv'].get('who')} note={e['kv'].get('note', '')}\n\n"
        for p in paths:
            fp = os.path.join(ROOT, p)
            if p in seen or not os.path.isfile(fp): continue
            seen.add(p); txt += digest_file(fp) + "\n"
        add(f"evidence/{tid}.txt", txt[:60000])

    # ---- changed code
    n_code = 0
    if not a.no_code:
        cf = [(os.path.getmtime(p), p) for p in code_files()]
        changed = sorted([x for x in cf if dt.datetime.fromtimestamp(x[0]).astimezone() >= since], reverse=True)
        for _, p in changed[:a.max_code]: add("code/" + os.path.relpath(p, ROOT), digest_file(p, cap=12000)); n_code += 1
        for _, p in changed[a.max_code:]: files.setdefault("TRUNCATED.txt", ""); files["TRUNCATED.txt"] += f"changed code not included (--max-code): {os.path.relpath(p, ROOT)}\n"
        cm = "".join(file_map(p) for _, p in sorted(cf, key=lambda x: x[1]))
        add("codemap.txt", cm[:250000] + ("\n[codemap truncated]\n" if len(cm) > 250000 else ""))

    # ---- state, re-checks, drift
    counts = {s: sum(v == s for v in st.values()) for s in SHOW_STATES}
    stale, noev, drift, noprov = [], [], [], []
    for tid in [t for t in tasks if st[t] == "DONE"]:
        bad = [f"{r[0]} {r[1]} [{r[3]}]" for r in evaluate(tasks[tid], True) if r[2] == "FAIL"]
        if bad: stale.append(f"- {tid}: " + "; ".join(bad))
        for p in last[tid]["kv"].get("ev", "").split(","):
            if p and p != "chk" and not os.path.isfile(os.path.join(ROOT, p)): noev.append(f"- {tid}: ev path missing: {p}")
        s, lines = prov_drift(tid)
        if s == "drift": drift += [f"- {tid}: {x}" for x in lines]
        elif s == "none": noprov.append(tid)
    rows = [f"| {t} | {st[t]} | {e['kv'].get('who', '')} | {e['kv'].get('note', '')[:140]} |" for t, e in sorted(last.items())
            if st[t] != "DONE" or e["t"] >= since][:400]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): cmd_validate(None)
    vl = buf.getvalue().strip().splitlines(); groups = collections.OrderedDict(); keep = []
    for l in vl:
        if "warn:" in l:
            k = re.sub(r"'.*", "", l.split("warn: ", 1)[1]).strip(); groups.setdefault(k, []).append(l)
        else: keep.append(l)
    validate = "\n".join(keep[:-1] + [f"WARNINGS BY KIND ({sum(len(v) for v in groups.values())} total; `tl.py validate` lists all):"]
                         + [f"- {k}: {len(v)}  e.g. {v[0].split(': warn: ')[0]}" for k, v in groups.items()] + keep[-1:])
    sbuf = io.StringIO()
    with contextlib.redirect_stdout(sbuf): cmd_status(argparse.Namespace(n=10, prefix=None))
    try:
        from importlib import metadata
        pk = sorted({f"{d.metadata['Name']}=={d.version}" for d in metadata.distributions() if d.metadata["Name"]}, key=str.lower)[:300]
    except Exception: pk = []
    shas = " | ".join(f"{bn(f)} {sha(f)}" for f in files_ + [os.path.abspath(__file__)])
    summary = f"""# Review bundle {now.isoformat(timespec='minutes')}
repo root: {ROOT}
sha256: {shas}
window: log lines since {since.isoformat(timespec='minutes')}; code files changed in window: {n_code}
How to read: this is a digest, not the project. Ask for exact pieces with `tl.py pull` (see REVIEW_PROTOCOL.md).

## State counts (from the log)
{', '.join(f'{k} {v}' for k, v in counts.items())}

## Review queue (digests in evidence/)
{chr(10).join(f'- {t}: {w}' for t, w in queue.items()) or 'empty'}
Not digested: {len([t for t in done_win if t not in queue])} other DONE tasks in the window (their checks are re-run below).

## DONE tasks whose checks fail NOW (cmd checks skipped)
{chr(10).join(stale) or 'none'}

## DONE tasks whose evidence changed after DONE (hash drift)
{chr(10).join(drift) or 'none'}
DONE tasks with no provenance file (logged before prov/ existed): {len(noprov)}

## DONE tasks whose evidence path does not exist
{chr(10).join(noev) or 'none'}

## Validate (errors + warnings)
{validate.strip()}

## Latest line per task (non-DONE, plus anything changed in the window)
| task | state | who | note |
|---|---|---|---|
""" + "\n".join(rows) + "\n"
    files["SUMMARY.md"] = summary
    files["status.txt"] = sbuf.getvalue()
    files["log.txt"] = "".join(f"## {bn(log_for(f))}\n" + open(log_for(f), encoding="utf-8").read() for f in files_ if os.path.exists(log_for(f)))[-300000:]
    files["git.txt"] = ("\n\n".join(f"$ {c}\n{sh(c)}" for c in ["git status --short | head -80", "git log --oneline -20", "git diff --stat | tail -30"])
                        if in_git() else "not a git repository")
    files["env.txt"] = f"python {sys.version.split()[0]} ({sys.executable}) on {platform.platform()}\n" + "\n".join(pk) + "\n"
    out = a.out or os.path.join(TDIR, f"bundle_{now.strftime('%Y%m%d_%H%M%S')}.zip")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for n, t in files.items(): z.writestr(n, t)
    open(LAST_BUNDLE, "w").write(now.isoformat(timespec="minutes"))
    print(f"wrote {out} ({os.path.getsize(out) // 1024} KB, {len(files)} files, {used[0] // 1000}k chars)")
    print(f"DONE {counts['DONE']}/{len(tasks)}; queue {len(queue)}; stale DONE {len(stale)}; drift {len(drift)}; missing evidence {len(noev)}; truncated: {'yes' if 'TRUNCATED.txt' in files else 'no'}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="c", required=True)
    sub.add_parser("validate").set_defaults(f=cmd_validate)
    s = sub.add_parser("status"); s.add_argument("-n", type=int, default=10); s.add_argument("-p", "--prefix"); s.set_defaults(f=cmd_status)
    s = sub.add_parser("next"); s.add_argument("-p", "--prefix"); s.set_defaults(f=cmd_next)
    s = sub.add_parser("show"); s.add_argument("id"); s.set_defaults(f=cmd_show)
    s = sub.add_parser("check"); s.add_argument("id"); s.add_argument("--cmd", action="store_true"); s.set_defaults(f=cmd_check)
    s = sub.add_parser("audit"); s.add_argument("--write", action="store_true"); s.add_argument("--cmd", action="store_true"); s.set_defaults(f=cmd_audit)
    s = sub.add_parser("log"); s.add_argument("id"); s.add_argument("state")
    s.add_argument("--note"); s.add_argument("--ev"); s.add_argument("--who"); s.add_argument("--ref")
    s.add_argument("--file"); s.add_argument("--audit", action="store_true"); s.set_defaults(f=cmd_log)
    s = sub.add_parser("bundle"); s.add_argument("--since"); s.add_argument("--ids"); s.add_argument("--all", action="store_true")
    s.add_argument("--sample", type=int, default=6); s.add_argument("--budget", type=int, default=250_000)
    s.add_argument("--no-code", action="store_true"); s.add_argument("--max-code", type=int, default=30); s.add_argument("--out"); s.set_defaults(f=cmd_bundle)
    s = sub.add_parser("pull"); s.add_argument("path"); s.add_argument("--lines"); s.add_argument("--grep"); s.add_argument("-C", type=int, default=2)
    s.add_argument("--sym"); s.add_argument("--rows"); s.add_argument("--where"); s.add_argument("--cols"); s.add_argument("--cap", type=int, default=6000); s.set_defaults(f=cmd_pull)
    s = sub.add_parser("map"); s.add_argument("paths", nargs="*"); s.set_defaults(f=cmd_map)
    a = ap.parse_args()
    sys.exit(a.f(a))


if __name__ == "__main__":
    main()
