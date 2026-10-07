#!/usr/bin/env python3
"""固定幅の独立 A/B (plan architecture-solver-host-memory §6.3、2026-10-07 事前登録)。AWS で回す。

    # 段階 1: 既存の base 3 本だけから診断幅 T = 2·S0 を凍結する (new は使わない)
    python3 fixedwidth_eval.py freeze --base RUN RUN RUN --out fixedwidth_c44dual_ckpt100
    # 段階 2: 新しい 12 本 (plan.json の順・名前) を、凍結した T と既存 base 3 本で評価する
    python3 fixedwidth_eval.py eval --plan fixedwidth_c44dual_ckpt100/plan.json [--out fixedwidth_c44dual_ckpt100/result]

段階 1 (freeze):
  量は compare_runs.py の追加診断 B と同じ集合 (最終出力 h5 の全データセット・CSV の全値列) に、表の行キー
  (step・inner・phase 等) と初期出力 (最初の res_*.h5) の保存量・原始量・幾何量を足したもの。
  S0 = max_{i<j} ||B_i − B_j||_∞ (float64 の絶対 L∞)、T = 2·S0。S0 = 0 の量は「差 0 を要求」。
  整数・文字列・行キー・初期出力は厳密一致。base 3 本の間で構造 (ファイル集合・データセット集合・shape・dtype・属性) が
  一致しなければ凍結しない (止める)。T_frozen.tsv と、その sha256 を T_frozen.sha256 に書く。

段階 2 (eval):
  plan.json の T_frozen の sha256 を照合してから評価する (T を後から変えると止まる)。新しい各 run X について
  - 出力ファイル集合・各 h5 のデータセット集合・shape・dtype・属性 (CHECKPOINT の totalTime 等 = 物理時刻を含む) が既存 B1 と厳密一致
  - 表の行キー (step・inner・phase) が B1 と厳密一致
  - 初期出力の保存量・原始量・幾何量と整数・文字列の量が B1 とビット一致
  - 幅の量: v = max_i ||X − B_i||_∞ ≤ T (差 0 を要求する量は v = 0)。非有限・shape 違い・欠落は失敗
  新規 base 6 本 (b) と new 6 本 (n) を分けて集計し、plan §6.3 の 3 行の表で判定する。
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_runs as cr  # noqa: E402

MODES = ("幅", "差0を要求", "厳密一致(整数・文字列)", "厳密一致(行キー)", "厳密一致(初期出力)")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def first_res(d):
    st = cr.res_steps(d)
    return f"res_{st[0]}.h5" if len(st) > 1 else None


def structure(d):
    """run の構造: 最終ファイル・表・初期出力の名前と、各 h5 の (データセット → shape/dtype/属性, グループ属性)、各表の (列, 行キー)。"""
    finals, tables = cr.final_files(d)
    f0 = first_res(d)
    h5 = {fn: cr.h5_items(os.path.join(d, fn)) for fn in finals + ([f0] if f0 else [])}
    tb = {}
    for fn in tables:
        t, keys = cr.read_table(os.path.join(d, fn))
        tb[fn] = (list(t.keys()), keys)
    return dict(finals=finals, tables=tables, init=f0, h5=h5, tb=tb)


def struct_diff(sa, sb):
    """2 つの run の構造の違いを文字列の list で返す (空なら一致)。"""
    out = []
    for k in ("finals", "tables", "init"):
        if sa[k] != sb[k]:
            out.append(f"{k}: {sa[k]} vs {sb[k]}")
    for fn in set(sa["h5"]) & set(sb["h5"]):
        da, ga = sa["h5"][fn]
        db, gb = sb["h5"][fn]
        if set(da) != set(db):
            out.append(f"{fn}: データセット集合 (片側だけ {sorted(set(da) ^ set(db))[:6]})")
        for n in set(da) & set(db):
            if da[n] != db[n]:
                out.append(f"{fn}:{n}: shape/dtype/属性 {da[n]} vs {db[n]}")
        if ga != gb:
            out.append(f"{fn}: グループ属性 {ga} vs {gb}")
    for fn in set(sa["tb"]) & set(sb["tb"]):
        if sa["tb"][fn][0] != sb["tb"][fn][0]:
            out.append(f"{fn}: 列が違う")
        if sa["tb"][fn][1] != sb["tb"][fn][1]:
            out.append(f"{fn}: 行キーが対応しない")
    return out


def load_quantity(d, fn, name):
    """量の値 (h5 のデータセットか表の列)。無ければ None。"""
    if fn.endswith(".h5"):
        its, _ = cr.h5_items(os.path.join(d, fn))
        return cr.h5_data(os.path.join(d, fn), [name])[name] if name in its else None
    t, keys = cr.read_table(os.path.join(d, fn))
    if name == "(行キー)":
        return np.array(["|".join(k) for k in keys], dtype=object)
    return t.get(name)


def quantities(st):
    """(ファイル, 量, 種類) の list。種類は float / exact / key / init。"""
    q = []
    for fn in st["finals"]:
        for n in sorted(st["h5"][fn][0]):
            dt = st["h5"][fn][0][n][1]
            q.append((fn, n, "float" if dt.startswith("float") else "exact"))
    for fn in st["tables"]:
        cols = st["tb"][fn][0]
        for n in sorted(set(cols) - set(cr.TABLE_KEYS)):
            q.append((fn, n, "table"))
        q.append((fn, "(行キー)", "key"))
    if st["init"]:
        for n in sorted(st["h5"][st["init"]][0]):
            if (not n.startswith("VALUE/")) or cr.T0_VALUE.match(n.split("/", 1)[1]):
                q.append((st["init"], n, "init"))
    return q


def exact_equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and a.dtype == b.dtype and (a.tobytes() == b.tobytes() if a.dtype != object
                                                          else bool(np.all(a == b)))


# ------------------------------------------------------------------ 段階 1
def cmd_freeze(a):
    base = [os.path.abspath(os.path.expanduser(x)) for x in a.base]
    if len(base) != 3:
        sys.exit("[freeze] base はちょうど 3 本")
    sts = [structure(d) for d in base]
    for k in (1, 2):
        dif = struct_diff(sts[0], sts[k])
        if dif:
            sys.exit(f"[freeze] base 3 本の構造が一致しない ({os.path.basename(base[k])}): {dif[:5]} — 凍結しない")
    tag = {i: f"B{i + 1}" for i in range(3)}
    rows = []
    for fn, n, kind in quantities(sts[0]):
        vals = [load_quantity(d, fn, n) for d in base]
        arrs = [None if v is None else np.asarray(v) for v in vals]
        if any(x is None for x in arrs):
            sys.exit(f"[freeze] {fn}:{n} が base の一部に無い")
        isfloat = all(x.dtype.kind == "f" for x in arrs)
        if kind in ("float", "table") and isfloat:
            best, pair = -1.0, None
            for i in range(3):
                for j in range(i + 1, 3):
                    v = cr.dabs(arrs[i], arrs[j])
                    if v is None:
                        sys.exit(f"[freeze] {fn}:{n} の base {tag[i]}-{tag[j]} が比較不能 (非有限か shape 違い) — 凍結しない")
                    if v > best:
                        best, pair = v, (i, j)
            mode = "幅" if best > 0 else "差0を要求"
            rows.append((fn, n, mode, best, 2.0 * best, f"{tag[pair[0]]}-{tag[pair[1]]}" if best > 0 else "全組で差 0"))
        else:
            same = all(exact_equal(arrs[0], x) for x in arrs[1:])
            mode = {"key": "厳密一致(行キー)", "init": "厳密一致(初期出力)"}.get(kind, "厳密一致(整数・文字列)")
            if kind == "init" and isfloat:
                s0 = max(cr.dabs(arrs[i], arrs[j]) for i in range(3) for j in range(i + 1, 3))
            else:
                s0 = 0.0 if same else float("nan")
            if not same:
                sys.exit(f"[freeze] {fn}:{n} ({mode}) が base 3 本で一致しない (S0 {s0}) — 凍結しない")
            rows.append((fn, n, mode, s0, 0.0, "-"))
    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, "T_frozen.tsv")
    nw = sum(1 for r in rows if r[2] == "幅")
    nz = sum(1 for r in rows if r[2] == "差0を要求")
    ne = sum(1 for r in rows if r[2] == "厳密一致(整数・文字列)")
    nk = sum(1 for r in rows if r[2] == "厳密一致(行キー)")
    ni = sum(1 for r in rows if r[2] == "厳密一致(初期出力)")
    with open(p, "w") as f:
        f.write(f"# 固定幅の独立 A/B の診断幅 (plan architecture-solver-host-memory §6.3)。凍結 {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")
        f.write(f"# 既存 base 3 本だけから算定 (new は使わない): " + " ".join(f"{tag[i]}={os.path.basename(d)}" for i, d in enumerate(base)) + "\n")
        f.write("# S0 = max_{i<j} ||B_i − B_j||_∞ (float64 の絶対 L∞)、T = 2·S0。mode: 幅 = |X − B_i| ≤ T、差0を要求 = S0 = 0 なので X − B_i = 0 を要求、"
                "厳密一致 = ビット・値の一致 (S0 は情報、T は使わない)\n")
        f.write(f"# 量の数: 最終出力の量 {nw + nz + ne} (幅 {nw}、差0を要求 {nz}、厳密一致(整数・文字列) {ne}) + 行キー {nk} + 初期出力 {ni}\n")
        f.write("# 構造 (ファイル集合・データセット集合・shape・dtype・属性・CHECKPOINT の totalTime 等) は base 3 本で一致を確認済み。段階 2 では B1 との厳密一致を要求する\n")
        f.write("file\tquantity\tmode\tS0\tT\tS0_pair\n")
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.9e}\t{r[4]:.9e}\t{r[5]}\n")
    h = sha256(p)
    open(os.path.join(a.out, "T_frozen.sha256"), "w").write(f"{h}  T_frozen.tsv\n")
    print(f"[freeze] {p}: 最終出力の量 {nw + nz + ne} (幅 {nw}、差0を要求 {nz}、厳密一致 {ne}) + 行キー {nk} + 初期出力 {ni}")
    print(f"[freeze] sha256 {h}")


# ------------------------------------------------------------------ 段階 2
def read_T(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or line.startswith("file\t"):
            continue
        fn, n, mode, s0, t, pair = line.rstrip("\n").split("\t")
        rows.append(dict(file=fn, name=n, mode=mode, S0=float(s0), T=float(t), pair=pair))
    return rows


def eval_run(x, base, base_st, trows):
    """新しい run x を凍結した T と既存 base 3 本で評価。返り値 (超過の list, 構造の違いの list)。"""
    st = structure(x)
    sdiff = struct_diff(base_st, st)
    over = []
    for r in trows:
        xv = load_quantity(x, r["file"], r["name"])
        if xv is None:
            over.append((r, float("nan"), "比較不能: 欠落"))
            continue
        if r["mode"].startswith("厳密一致"):
            b1 = load_quantity(base[0], r["file"], r["name"])
            if not exact_equal(b1, xv):
                v = cr.dabs(b1, xv) if np.asarray(xv).dtype.kind == "f" else float("nan")
                over.append((r, v if v is not None else float("nan"), "厳密一致でない (B1 と)"))
            continue
        worst, wi = -1.0, None
        bad = ""
        for i, b in enumerate(base):
            v = cr.dabs(load_quantity(b, r["file"], r["name"]), xv)
            if v is None:
                bad = f"比較不能 (B{i + 1} と: 非有限か shape 違い)"
                break
            if v > worst:
                worst, wi = v, i
        if bad:
            over.append((r, float("nan"), bad))
        elif worst > r["T"]:
            over.append((r, worst, f"幅超過 (最悪 B{wi + 1})"))
    return over, sdiff


def cmd_eval(a):
    plan = json.load(open(a.plan))
    pdir = os.path.dirname(os.path.abspath(a.plan))
    tpath = os.path.join(pdir, "T_frozen.tsv")
    if sha256(tpath) != plan["T_frozen_sha256"]:
        sys.exit("[eval] T_frozen.tsv の sha256 が plan.json と違う — 凍結後に T が変わった。評価しない")
    root = os.path.join(HERE)
    base = [os.path.join(root, r) for r in plan["base_existing"]]
    base_st = structure(base[0])
    trows = read_T(tpath)
    out = a.out or os.path.join(pdir, "result")
    os.makedirs(out, exist_ok=True)
    lines = [f"固定幅の独立 A/B (plan §6.3) 段階 2 の評価 {time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
             f"T_frozen.tsv sha256 {plan['T_frozen_sha256']} (照合済み)、既存 base: {', '.join(plan['base_existing'])}", ""]
    groups = {"b": [], "n": []}
    per_q = {}
    for job in plan["runs"]:
        x = os.path.join(root, job["run"])
        if os.path.exists(os.path.join(x, ".exclude")) and os.path.isdir(x + "_retry"):
            lines.append(f"{job['run']}: インフラ中断で除外 → {job['run']}_retry を使う")
            x = x + "_retry"
        if not os.path.isdir(x):
            lines.append(f"{job['run']}: run が無い — 判定不能")
            groups[job["group"]].append(None)
            continue
        over, sdiff = eval_run(x, base, base_st, trows)
        nan = open(os.path.join(x, "NANCHECK.txt")).read().strip().splitlines()[-1] if os.path.exists(os.path.join(x, "NANCHECK.txt")) else "NANCHECK: (無し)"
        bad = len(over) + len(sdiff) + (0 if "PASS" in nan else 1)
        groups[job["group"]].append(bad)
        lines.append(f"{job['run']} ({job['group']}): 超過・不一致 {len(over)}、構造の違い {len(sdiff)}、{nan}")
        for s in sdiff[:10]:
            lines.append(f"    構造: {s}")
        for idx, (r, v, why) in enumerate(over):
            if idx < 40:
                lines.append(f"    {r['file']}:{r['name']} [{r['mode']}] v {v:.3e} T {r['T']:.3e} — {why}")
            k = (r["file"], r["name"])
            per_q.setdefault(k, {"b": 0, "n": 0})[job["group"]] += 1
    b_ok = all(g == 0 for g in groups["b"]) and len(groups["b"]) == 6
    n_ok = all(g == 0 for g in groups["n"]) and len(groups["n"]) == 6
    b_any_none = any(g is None for g in groups["b"] + groups["n"])
    n_plan = {g: sum(1 for j in plan["runs"] if j["group"] == g) for g in ("b", "n")}
    if b_any_none or n_plan != {"b": 6, "n": 6}:
        verdict = f"判定不能 (run が欠けているか本数が計画 [b 6・n 6] と違う: {n_plan})"
    elif not b_ok:
        verdict = "新規 base も超過 → 判定不能 (基準 3 本では再現性を捉えられていない。T や反復数をその場で増やして合格にしない)"
    elif not n_ok:
        verdict = "新規 base は幅内、new だけ超過 → 変更起因の差を優先して追う (諮る)"
    else:
        verdict = "新規 base・new とも全量で幅内 → 「追加標本でもこの固定幅を超える差は検出されなかった」(限定した支持材料として result レビューへ)"
    lines += ["", "超過した量 (b 本数 / n 本数):"] + [f"  {k[0]}:{k[1]}  b {v['b']} / n {v['n']}" for k, v in sorted(per_q.items())]
    lines += ["", f"判定 (plan §6.3): {verdict}"]
    open(os.path.join(out, "fixedwidth_result.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("freeze")
    s.add_argument("--base", nargs=3, required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("eval")
    s.add_argument("--plan", required=True)
    s.add_argument("--out")
    a = ap.parse_args()
    {"freeze": cmd_freeze, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
