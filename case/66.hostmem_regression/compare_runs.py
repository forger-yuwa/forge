#!/usr/bin/env python3
"""case/66 回帰の比較 (plan architecture-solver-host-memory §6 の (a)(b)(c) と出力互換性)。AWS 上で回す (h5py・numpy)。

    # 構成ごと (registry.tsv から run を選ぶ。base だけなら同ビルド内のばらつきだけを出す)
    python3 compare_runs.py --cfg c36node [--base-build base] [--new-build new] [--out report.txt]
    python3 compare_runs.py --all [--new-build new] --out-dir compare/       # 全構成 (構成ごとの報告 + summary.txt)
    # run を直接指定
    python3 compare_runs.py --base RUN RUN RUN [--new RUN RUN RUN]
    # 2 ファイルの単純比較 (dual-time の分割 res_100 と連続 res_200 など)
    python3 compare_runs.py --diff2 A.h5 B.h5

判定の定義 (§6。結果を見てから変えない):
  (a) step 0 の残差行 (residual_history.csv の最初の step の全行) を全列ビット一致で比べる。base 反復の間で値が割れる列
      (既知の 1 ulp の非決定性、atomicAdd) は「new の値が base で観測した値のどれか」で可。
      初期出力 (最初の res_*.h5 = writeInitialOutputs、残差の組立より前) では**保存量・原始量・幾何量だけ**をビット一致で比べる
      (残差・補正量・勾配・リミッタ・診断は初期出力では比べない。level 2 の初期出力には未初期化の残差が入る: §6 M4)。
  (b) 最後の出力 (最終 res_*.h5・境界出力・CSV 出力・残差履歴) の各データセット (列) について
      m(A,B) = max|A−B| / max|A| (A は組の前側、max|A| = 0 なら差が 0 のとき 0、そうでなければ inf)。
      S = 同ビルド内ペア差 (base 3 回の 3 対 + new 3 回の 3 対) の最大、D = base×new の 9 対の最大。合格は D ≤ 2·S
      (S = 0 のときは D = 0)。new が無いときは S_base だけを出す (= 基準のばらつき)。
  (c) 各 run の NaN/Inf (run_matrix.py が残した NANCHECK.txt と、ここでの再検査)。
  出力互換: 全 h5 のデータセット集合・shape・dtype・属性 (値も) が一致すること。属性は全 run で比べる。
  ログ: 起動時の名前登録・checkpoint 復元・遷移初期化・ψ 退避件数・警告の行が全 run で一致すること。
"""
import argparse
import csv
import glob
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# 初期出力で比べる量 (保存量・原始量・幾何量)。/VALUE 以外 (/MESH・/CHECKPOINT・/BCONDS 等) は幾何・状態として全部比べる。
T0_VALUE = re.compile(
    r"^(ro|roUx|roUy|roUz|roe|roK|roOmega|roY\d+|roXi|roGamma|roReth|ro[gQ]\w*_\d+|"
    r"P|T|Ux|Uy|Uz|k|omega|Y\d+|Xi|h0|[gQ][0-2]?_\d+|wall_dist|volume|ccx|ccy|ccz)$")
# 比較対象のファイル (run が書いたもの)。入力 (INPUT_FILES) と xmf・ログは除く。
OUT_GLOBS = ["*.h5", "*.csv", "*.out"]
SKIP_FILES = {"residual_history.png"}
# ログから拾う行 (時間・速度を含む行は拾わない)
LOG_PAT = re.compile(r"(registered|restored|history|\[variables\]|snapshot:|\[FORGE_OUT_RESIDUALS\]|\[FORGE_RESID_SNAP\]|"
                     r"psi-dualeval\] ON|WARNING|warning|警告|ignored|無視|lineImplicit|line-implicit|"
                     r"\[transition\]|transition .* (read|initiali)|output: |'output')", re.I)
LOG_DROP = re.compile(r"(ms/step|elapsed|eta |wall|Time = |sec|秒)", re.I)


# ------------------------------------------------------------------ run の選択
def read_registry():
    p = os.path.join(HERE, "registry.tsv")
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def runs_for(cfg, build):
    out = []
    for r in read_registry():
        if r["cfg"] != cfg or r["build"] != build:
            continue
        d = os.path.join(HERE, r["run"])
        st = open(os.path.join(d, ".state")).read() if os.path.exists(os.path.join(d, ".state")) else ""
        if os.path.exists(os.path.join(d, ".exclude")):   # 破棄予定 (起動失敗など) は比較に使わない
            continue
        if st.startswith("done"):
            out.append(d)
    return out


def input_files(d):
    p = os.path.join(d, "INPUT_FILES")
    return set(open(p).read().split()) if os.path.exists(p) else set()


def out_files(d):
    inp = input_files(d)
    fs = set()
    for g in OUT_GLOBS:
        for p in glob.glob(os.path.join(d, g)):
            b = os.path.basename(p)
            if b not in inp and b not in SKIP_FILES:
                fs.add(b)
    return fs


def res_steps(d):
    st = []
    for p in glob.glob(os.path.join(d, "res_*.h5")):
        m = re.match(r"res_(\d+)\.h5$", os.path.basename(p))
        if m:
            st.append(int(m.group(1)))
    return sorted(st)


# ------------------------------------------------------------------ h5
def h5_items(path):
    """{データセット名: (shape, dtype, attrs)} と {グループ名: attrs}。"""
    ds, grp = {}, {}
    with h5py.File(path, "r") as f:
        grp["/"] = {k: _attr(v) for k, v in f.attrs.items()}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                ds[name] = (obj.shape, str(obj.dtype), {k: _attr(v) for k, v in obj.attrs.items()})
            else:
                grp[name] = {k: _attr(v) for k, v in obj.attrs.items()}
        f.visititems(visit)
    return ds, grp


def _attr(v):
    if isinstance(v, bytes):
        return v.decode(errors="replace")
    if isinstance(v, np.ndarray):
        return tuple(v.ravel().tolist())
    if isinstance(v, np.generic):
        return v.item()
    return v


def h5_data(path, names):
    with h5py.File(path, "r") as f:
        return {n: f[n][()] for n in names}


def metric(a, b):
    """max|a−b| / max|a| (整数は不一致の個数を返さず、一致なら 0、不一致なら inf)。"""
    if a.shape != b.shape:
        return float("inf")
    if a.dtype.kind in "iub":
        return 0.0 if np.array_equal(a, b) else float("inf")
    a64 = np.asarray(a, dtype=np.float64)
    b64 = np.asarray(b, dtype=np.float64)
    if a64.size == 0:
        return 0.0
    d = np.nanmax(np.abs(a64 - b64)) if np.all(np.isfinite(a64)) and np.all(np.isfinite(b64)) else float("inf")
    s = np.max(np.abs(a64))
    if s == 0.0:
        return 0.0 if d == 0.0 else float("inf")
    return float(d / s)


def bit_equal(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


# ------------------------------------------------------------------ CSV
def read_table(path):
    """数値表を {列名: ndarray} にする。区切りは ',' (probe は ' , ')。行キー (step,inner,phase) を別に返す。"""
    with open(path) as f:
        lines = [x.rstrip("\n") for x in f if x.strip()]
    if not lines:
        return {}, []
    sep = ","
    hdr = [h.strip() for h in lines[0].split(sep)]
    rows = [[c.strip() for c in ln.split(sep)] for ln in lines[1:]]
    keys = []
    cols = {h: [] for h in hdr}
    keycols = [i for i, h in enumerate(hdr) if h in ("step", "inner", "phase", "Step", "var", "physID")]
    for r in rows:
        keys.append(tuple(r[i] if i < len(r) else "" for i in keycols))
        for i, h in enumerate(hdr):
            cols[h].append(r[i] if i < len(r) else "")
    out = {}
    for h, v in cols.items():
        try:
            out[h] = np.array([float(x) for x in v])
        except ValueError:
            out[h] = np.array(v, dtype=object)
    return out, keys


def step0_rows(path):
    with open(path) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        rows = list(rd)
    if not rows:
        return hdr, []
    s0 = rows[0][0]
    return hdr, [r for r in rows if r[0] == s0]


# ------------------------------------------------------------------ 比較本体
class Report:
    def __init__(self):
        self.lines = []
        self.verdicts = {}

    def p(self, s=""):
        self.lines.append(s)

    def v(self, key, val):
        self.verdicts[key] = val
        self.p(f"  >> {key}: {val}")


def compare(base, new, rep, cfgname=""):
    allruns = base + new
    tag = {d: ("B" if d in base else "N") + str((base if d in base else new).index(d) + 1) for d in allruns}
    rep.p(f"=== 構成 {cfgname}")
    for d in allruns:
        rep.p(f"  {tag[d]}: {os.path.relpath(d, HERE)}")

    # ---- (c) NaN/Inf
    rep.p("\n[c] NaN/Inf (各 run)")
    nan_bad = 0
    for d in allruns:
        p = os.path.join(d, "NANCHECK.txt")
        v = open(p).read().strip().splitlines()[-1] if os.path.exists(p) else "NANCHECK: (無し)"
        if "PASS" not in v:
            nan_bad += 1
        rep.p(f"  {tag[d]}: {v}")
    rep.v("NAN", "PASS" if nan_bad == 0 else f"FAIL ({nan_bad} run)")

    # ---- 出力ファイル集合
    rep.p("\n[構造] 出力ファイルの集合")
    fsets = {d: out_files(d) for d in allruns}
    ref = fsets[allruns[0]]
    fs_bad = [d for d in allruns if fsets[d] != ref]
    for d in fs_bad:
        rep.p(f"  {tag[d]}: 余分 {sorted(fsets[d] - ref)} / 欠け {sorted(ref - fsets[d])}")
    rep.p(f"  {len(ref)} ファイル: {', '.join(sorted(ref))}")
    common = set.intersection(*fsets.values())
    h5s = sorted(f for f in common if f.endswith(".h5"))

    # ---- h5 の構造 (集合・shape・dtype・属性)
    rep.p("\n[構造] h5 のデータセット集合・shape・dtype・属性")
    st_bad = len(fs_bad)
    items = {}
    for fn in h5s:
        items[fn] = {d: h5_items(os.path.join(d, fn)) for d in allruns}
        r_ds, r_grp = items[fn][allruns[0]]
        for d in allruns[1:]:
            ds, grp = items[fn][d]
            if set(ds) != set(r_ds):
                st_bad += 1
                rep.p(f"  {fn} {tag[d]}: データセット集合が違う 余分 {sorted(set(ds) - set(r_ds))[:8]} 欠け {sorted(set(r_ds) - set(ds))[:8]}")
            for n in set(ds) & set(r_ds):
                if ds[n][:2] != r_ds[n][:2]:
                    st_bad += 1
                    rep.p(f"  {fn}:{n} {tag[d]}: shape/dtype {ds[n][:2]} vs {r_ds[n][:2]}")
                if ds[n][2] != r_ds[n][2]:
                    st_bad += 1
                    rep.p(f"  {fn}:{n} {tag[d]}: 属性 {ds[n][2]} vs {r_ds[n][2]}")
            for g in set(grp) | set(r_grp):
                if grp.get(g) != r_grp.get(g):
                    st_bad += 1
                    rep.p(f"  {fn} グループ {g} {tag[d]}: 属性 {grp.get(g)} vs {r_grp.get(g)}")
        rep.p(f"  {fn}: {len(r_ds)} データセット")
    rep.v("STRUCT", "PASS" if st_bad == 0 else f"FAIL ({st_bad} 件)")

    # ---- ログの行
    rep.p("\n[ログ] 名前登録・復元・初期化・警告の行")
    logs = {}
    for d in allruns:
        p = os.path.join(d, "forge_run.log") if os.path.exists(os.path.join(d, "forge_run.log")) else os.path.join(d, "convert.log")
        seen = []
        if os.path.exists(p):
            with open(p, errors="replace") as f:
                for ln in f:
                    if LOG_PAT.search(ln) and not LOG_DROP.search(ln):
                        x = ln.rstrip()
                        if x not in seen:
                            seen.append(x)
        logs[d] = seen
    lg_bad = 0
    for d in allruns[1:]:
        if logs[d] != logs[allruns[0]]:
            lg_bad += 1
            a, b = set(logs[allruns[0]]), set(logs[d])
            rep.p(f"  {tag[d]}: B1 に無い行 {sorted(b - a)[:5]} / B1 にだけある行 {sorted(a - b)[:5]}")
    for x in logs[allruns[0]][:40]:
        rep.p(f"    | {x[:200]}")
    rep.v("LOG", "PASS" if lg_bad == 0 else f"FAIL ({lg_bad} run)")

    # ---- (a) step 0 の残差行
    rep.p("\n[a] step 0 の残差行 (全列ビット一致; base で割れた列は new が base の観測値のどれか)")
    rc = [d for d in allruns if os.path.exists(os.path.join(d, "residual_history.csv"))]
    if len(rc) == len(allruns) and rc:
        hdr0, _ = step0_rows(os.path.join(rc[0], "residual_history.csv"))
        rows = {d: step0_rows(os.path.join(d, "residual_history.csv"))[1] for d in allruns}
        nrows = {len(v) for v in rows.values()}
        a_bad = 0
        split_cols = []
        if len(nrows) != 1:
            a_bad += 1
            rep.p(f"  step 0 の行数が違う: {[len(rows[d]) for d in allruns]}")
        else:
            for i in range(nrows.pop()):
                for j, h in enumerate(hdr0):
                    bvals = {rows[d][i][j] for d in base}
                    if len(bvals) > 1:
                        split_cols.append((i, h, len(bvals)))
                    for d in new:
                        if rows[d][i][j] not in bvals:
                            a_bad += 1
                            rep.p(f"  行 {i} 列 {h}: {tag[d]} {rows[d][i][j]} ∉ base {sorted(bvals)}")
        rep.p(f"  base で値が割れた (行, 列, 値の数): {split_cols if split_cols else 'なし (base 3 回でビット一致)'}")
        rep.v("STEP0", ("PASS" if a_bad == 0 else f"FAIL ({a_bad} 値)") if new else
              f"BASE-ONLY (割れた列 {len(split_cols)})")
    else:
        rep.p("  residual_history.csv が無い (変換器など)")

    # ---- 初期出力 (保存量・原始量・幾何量のビット一致)
    steps = res_steps(allruns[0])
    if steps:
        fn0 = f"res_{steps[0]}.h5"
        rep.p(f"\n[a] 初期出力 {fn0}: 保存量・原始量・幾何量のビット一致 (他の量は情報として不一致の数だけ)")
        ds0 = items.get(fn0, {}).get(allruns[0], ({}, {}))[0]
        cmp_names = [n for n in ds0 if (not n.startswith("VALUE/")) or T0_VALUE.match(n.split("/", 1)[1])]
        info_names = [n for n in ds0 if n not in cmp_names]
        ref = h5_data(os.path.join(allruns[0], fn0), cmp_names + info_names)
        t0_bad, info_diff = [], {}
        for d in allruns[1:]:
            x = h5_data(os.path.join(d, fn0), cmp_names + info_names)
            for n in cmp_names:
                if not bit_equal(ref[n], x[n]):
                    t0_bad.append((tag[d], n, metric(ref[n], x[n])))
            for n in info_names:
                if not bit_equal(ref[n], x[n]):
                    info_diff.setdefault(n, []).append(tag[d])
        rep.p(f"  比べた量 {len(cmp_names)}、比べない量 {len(info_names)}")
        for t, n, m in t0_bad[:20]:
            rep.p(f"  不一致 {t} {n}: m={m:.3e}")
        if info_diff:
            rep.p(f"  (情報) 比べない量で B1 と違うもの: {sorted(info_diff)[:20]}")
        rep.v("INIT_OUT", "PASS" if not t0_bad else f"FAIL ({len(t0_bad)} 件)")
        if len(steps) == 1:
            rep.p("  (出力が 1 回だけ)")

    # ---- (b) 最後の出力・境界出力・CSV・残差履歴
    rep.p("\n[b] N step 後: m = max|A−B|/max|A|。S = 同ビルド内ペアの最大、D = base×new の最大、合格 D ≤ 2·S")
    finals = []
    if steps:
        last = steps[-1]
        finals = [f for f in h5s if re.search(rf"_{last}\.h5$", f) and (len(steps) > 1 or True)]
        if len(steps) > 1:
            finals = [f for f in finals if f != f"res_{steps[0]}.h5"]
    if not steps:
        finals = h5s   # 変換器
    tables = sorted(f for f in common if f.endswith((".csv", ".out")))
    pairs_b = [(base[i], base[j]) for i in range(len(base)) for j in range(i + 1, len(base))]
    pairs_n = [(new[i], new[j]) for i in range(len(new)) for j in range(i + 1, len(new))]
    pairs_x = [(b, n) for b in base for n in new]
    rows_out = []   # (ファイル, 量, S_base, S_new, D, verdict)

    def eval_set(getter, names, fname):
        cache = {d: getter(d) for d in allruns}
        for n in names:
            def m(p):
                a, b = cache[p[0]].get(n), cache[p[1]].get(n)
                if a is None or b is None:
                    return float("inf")
                return metric(a, b)
            sb = max((m(p) for p in pairs_b), default=float("nan"))
            sn = max((m(p) for p in pairs_n), default=float("nan"))
            dd = max((m(p) for p in pairs_x), default=float("nan"))
            if new:
                s = max(sb, sn)
                ok = (dd == 0.0) if s == 0.0 else (dd <= 2.0 * s)
                vd = "PASS" if ok else "FAIL"
            else:
                vd = "-"
            rows_out.append((fname, n, sb, sn, dd, vd))

    for fn in finals:
        names = sorted(items[fn][allruns[0]][0])
        eval_set(lambda d, fn=fn, names=names: h5_data(os.path.join(d, fn), names), names, fn)
    for fn in tables:
        def get(d, fn=fn):
            t, keys = read_table(os.path.join(d, fn))
            return {k: v for k, v in t.items() if v.dtype != object and k not in ("step", "inner", "Step")}
        t0, _ = read_table(os.path.join(allruns[0], fn))
        names = [k for k, v in t0.items() if v.dtype != object and k not in ("step", "inner", "Step")]
        eval_set(get, names, fn)

    b_fail = [r for r in rows_out if r[5] == "FAIL"]
    rep.p(f"  {'ファイル':28s} {'量':28s} {'S_base':>10s} {'S_new':>10s} {'D':>10s}  判定")
    # 主要なものから: 最終場の保存量、残差履歴、その他。S_base の大きい順に上位を出し、FAIL は全部出す
    def fmt(x):
        return "      -   " if x != x else f"{x:10.3e}"
    shown = 0
    for r in sorted(rows_out, key=lambda r: (r[5] != "FAIL", -(r[2] if r[2] == r[2] else 0))):
        if r[5] == "FAIL" or shown < 25:
            rep.p(f"  {r[0][:28]:28s} {r[1][:28]:28s} {fmt(r[2])} {fmt(r[3])} {fmt(r[4])}  {r[5]}")
            shown += 1
    nz = sum(1 for r in rows_out if r[2] == r[2] and r[2] > 0)
    rep.p(f"  ({len(rows_out)} 量。base 内で差が 0 でない量 {nz})")
    key = {}
    for r in rows_out:
        if r[0].startswith("res_") and r[1] in ("VALUE/ro", "VALUE/roUx", "VALUE/roe", "VALUE/P", "VALUE/T"):
            key[r[1]] = r[2]
        if r[0] == "residual_history.csv" and r[1] in ("rms_ro", "rms_roe"):
            key["csv:" + r[1]] = r[2]
    rep.p("  主要量の S_base: " + ", ".join(f"{k}={v:.2e}" for k, v in key.items()))
    if new:
        rep.v("NSTEP", "PASS" if not b_fail else f"FAIL ({len(b_fail)} 量)")
    else:
        smax = max((r[2] for r in rows_out if r[2] == r[2]), default=float("nan"))
        rep.v("NSTEP", f"BASE-ONLY (S_base 最大 {smax:.3e})")
    rep.rows = rows_out
    return rep


def diff2(a, b):
    da, _ = h5_items(a)
    db, _ = h5_items(b)
    names = sorted(set(da) & set(db))
    xa, xb = h5_data(a, names), h5_data(b, names)
    print(f"A={a}\nB={b}\n集合の差: A のみ {sorted(set(da) - set(db))[:10]} / B のみ {sorted(set(db) - set(da))[:10]}")
    for n in names:
        m = metric(xa[n], xb[n])
        if m != 0.0:
            print(f"  {n:40s} m={m:.3e}")
    print(f"  ({len(names)} 共通データセット、差 0 でないものだけ表示)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--base", nargs="*", default=[])
    ap.add_argument("--new", nargs="*", default=[])
    ap.add_argument("--base-build", default="base")
    ap.add_argument("--new-build")
    ap.add_argument("--out")
    ap.add_argument("--out-dir")
    ap.add_argument("--diff2", nargs=2)
    a = ap.parse_args()
    if a.diff2:
        diff2(*a.diff2)
        return
    jobs = []
    if a.all:
        import matrix_spec as ms
        for c in ms.CONFIGS:
            jobs.append((c, runs_for(c, a.base_build), runs_for(c, a.new_build) if a.new_build else []))
    elif a.cfg:
        jobs.append((a.cfg, runs_for(a.cfg, a.base_build), runs_for(a.cfg, a.new_build) if a.new_build else []))
    else:
        jobs.append(("(指定)", [os.path.abspath(x) for x in a.base], [os.path.abspath(x) for x in a.new]))
    summary = []
    for c, base, new in jobs:
        if len(base) < 2:
            summary.append(f"{c:22s} base {len(base)} 本 (2 本未満のため比較しない)")
            continue
        rep = compare(base, new, Report(), c)
        txt = "\n".join(rep.lines) + "\n"
        if a.out_dir:
            os.makedirs(a.out_dir, exist_ok=True)
            open(os.path.join(a.out_dir, f"{c}.txt"), "w").write(txt)
        elif a.out:
            open(a.out, "w").write(txt)
        else:
            print(txt)
        summary.append(f"{c:22s} base {len(base)} new {len(new)} | " +
                       " | ".join(f"{k} {v}" for k, v in rep.verdicts.items()))
    s = "\n".join(summary) + "\n"
    if a.out_dir:
        open(os.path.join(a.out_dir, "summary.txt"), "w").write(s)
    print(s)


if __name__ == "__main__":
    main()
