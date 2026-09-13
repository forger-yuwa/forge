#!/usr/bin/env python3
"""凝縮 run の軸 (対称線) 上の全エンタルピー保存を時系列で判定する (plans/active/condensation-float-speedup.md §4.3 (d), codex plan M6)。

  metric(step) = max_x |h0(x) − h0_in|  [J/kg]   (h0 = VALUE/h0, h0_in = 軸の最上流 3 % の平均)
  VERDICT: 末尾 3 枚の metric の変化 ≤ --steady (既定 10 J/kg) で STEADY、そうでなければ DRIFTING (新旧とも ≥3 枚・有限・共通 step が無ければ FAIL)。
  cell 場 (VALUE/h0 の長さ = セル数) は CELLS/centCoords、node 場は MESH/COORD で位置を取る。
  --ref <run> を与えると同じ step のスナップショットで metric の差 (new − ref) を出し、--allow (既定 100 J/kg) 以下と
  new ≤ max(--frac (既定 1e-3)·--h0scale (既定 3e5 J/kg ≈ c_p·300 K), ref + --steady) を要求する (twophase face-T fix の 0.1 % 規約; 参照が既に超える cell 場では悪化しないこと。datum 付きの h0 自体は ~0 なので尺度にしない)。

軸ノード: |y| ≤ --ytol (既定 1e-9 m) の節点。無ければ (中心線に節点の無い偶数分割メッシュ) x 列ごとに |y| 最小の節点。--xmin/--xmax で評価範囲 (凝縮帯) を絞る。
使い方: python3 solver_density_cuda/tools/cond_axis_h0.py <run_dir> [--ref <run_dir>] [--xmin X --xmax X]
"""
import argparse, glob, os, re, sys
import h5py, numpy as np

def mesh_file(run):
    cfg = open(os.path.join(run, "solverConfig.yaml")).read()
    m = re.search(r'meshFileName:\s*"?([^"\s]+)"?', cfg)
    return os.path.join(run, m.group(1))

def coords_for_values(run):
    """VALUE/* の位置座標: node 場 (長さ = 節点数) は MESH/COORD、cell 場 (長さ = セル数) は CELLS/centCoords (codex result M3)。"""
    with h5py.File(mesh_file(run), "r") as f:
        C = f["MESH/COORD"][...]; C = C.reshape(-1, 3) if C.ndim == 1 else C
        cc = None
        if "CELLS" in f:
            for k in ("centCoords", "centroid", "CENT"):
                if k in f["CELLS"]: cc = f["CELLS"][k][...]; cc = cc.reshape(-1, 3) if cc.ndim == 1 else cc; break
    n_val = None
    for step, fn in snapshots(run)[:1]:
        with h5py.File(fn, "r") as h:
            if "h0" in h["VALUE"]: n_val = h["VALUE/h0"].shape[0]
    if n_val is None: sys.exit("no res_*.h5 with VALUE/h0 in %s" % run)
    if n_val == C.shape[0]: return C, "node"
    if cc is not None and n_val == cc.shape[0]: return cc, "cell"
    sys.exit("VALUE/h0 length %d matches neither nodes (%d) nor cells (%s) in %s" % (n_val, C.shape[0], cc.shape[0] if cc is not None else "n/a", run))

def axis_nodes(run, ytol):
    """軸 (対称線) の点列: |y| ≤ ytol の点。無ければ (中心線に点が無い偶数分割メッシュ) x 列ごとに |y| 最小の点を取る。node/cell は場の長さで自動判別。"""
    C, kind = coords_for_values(run)
    y = np.abs(C[:, 1]); idx = np.where(y <= ytol)[0]
    if idx.size < 10:
        xr = np.round(C[:, 0], 6)
        order = np.lexsort((y, xr))               # x 列ごとに |y| 昇順
        xs_sorted = xr[order]
        first = np.concatenate(([True], xs_sorted[1:] != xs_sorted[:-1]))
        idx = order[first]
    return C[idx, 0], idx

def snapshots(run):
    fs = glob.glob(os.path.join(run, "res_*.h5"))
    out = []
    for f in fs:
        m = re.match(r"res_(\d+)\.h5$", os.path.basename(f))
        if m: out.append((int(m.group(1)), f))
    return sorted(out)

def metric(run, x, idx, xmin, xmax):
    rows = []
    for step, f in snapshots(run):
        with h5py.File(f, "r") as h:
            if "h0" not in h["VALUE"]: continue
            h0 = h["VALUE/h0"][...][idx].astype(np.float64)
        order = np.argsort(x); xs = x[order]; hs = h0[order]
        nin = max(3, int(0.03 * xs.size)); h0_in = hs[:nin].mean()
        sel = np.ones(xs.size, bool)
        if xmin is not None: sel &= xs >= xmin
        if xmax is not None: sel &= xs <= xmax
        d = np.abs(hs[sel] - h0_in); j = int(np.argmax(d))
        rows.append((step, float(d[j]), float(xs[sel][j]), float(h0_in)))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--ref"); ap.add_argument("--ytol", type=float, default=1.0e-9)
    ap.add_argument("--xmin", type=float); ap.add_argument("--xmax", type=float)
    ap.add_argument("--steady", type=float, default=10.0, help="末尾 3 枚の metric 変化の許容 [J/kg]")
    ap.add_argument("--allow", type=float, default=100.0, help="new − ref の許容 [J/kg]")
    ap.add_argument("--frac", type=float, default=1.0e-3, help="metric ≤ frac·h0scale (両 run)")
    ap.add_argument("--h0scale", type=float, default=3.0e5, help="h0 の尺度 [J/kg] (既定 3e5 ≈ c_p·300 K; datum 付き h0 は ~0 なので h0_in は尺度に使えない)")
    a = ap.parse_args()
    x, idx = axis_nodes(a.run, a.ytol)
    rows = metric(a.run, x, idx, a.xmin, a.xmax)
    if not rows: sys.exit("no res_*.h5 with VALUE/h0 in %s" % a.run)
    kind = coords_for_values(a.run)[1]
    print("run: %s  axis points: %d (%s field)  x %.4f..%.4f" % (a.run, idx.size, kind, x.min(), x.max()))
    if len(rows) < 3: print("FAIL: need >= 3 snapshots with VALUE/h0 (have %d)" % len(rows)); sys.exit(1)
    if not all(np.isfinite(r[1]) and np.isfinite(r[3]) for r in rows): print("FAIL: non-finite h0 metric"); sys.exit(1)
    print("%8s %14s %10s %14s" % ("step", "max|h0-h0in|", "at x", "h0_in"))
    for s, d, xx, hin in rows: print("%8d %14.3f %10.5f %14.1f" % (s, d, xx, hin))
    tail = [r[1] for r in rows[-3:]]
    steady = (max(tail) - min(tail)) <= a.steady if len(tail) >= 2 else False
    last = rows[-1]; ok = last[1] <= a.frac * a.h0scale   # 絶対水準 (参照が既に超えている場合は下の new ≤ ref+allow で判定; 情報として出す)
    verdict = "STEADY" if steady else "DRIFTING"
    print("h0 conservation: metric %.3f J/kg = %.2e of h0scale %.1e (limit %.0e) -> %s; tail change %.3f J/kg (limit %.1f) -> %s"
          % (last[1], last[1] / a.h0scale, a.h0scale, a.frac, "ok" if ok else "EXCEED", max(tail) - min(tail), a.steady, verdict))
    fail = (not steady) or ((not ok) and not a.ref)   # --ref があれば絶対水準の超過は「参照より悪化しない」で判定する
    if a.ref:
        xr, idr = axis_nodes(a.ref, a.ytol); rr = metric(a.ref, xr, idr, a.xmin, a.xmax)
        if len(rr) < 3 or not all(np.isfinite(r[1]) for r in rr): print("FAIL: ref needs >= 3 finite snapshots (have %d)" % len(rr)); sys.exit(1)
        tailr = [r[1] for r in rr[-3:]]; steady_r = (max(tailr) - min(tailr)) <= a.steady
        print("ref %s: last metric %.3f J/kg, tail change %.3f -> %s" % (a.ref, rr[-1][1], max(tailr) - min(tailr), "STEADY" if steady_r else "DRIFTING"))
        dr = {s: d for s, d, _, _ in rr}
        common = [(s, d, dr[s]) for s, d, _, _ in rows if s in dr]
        if not common: print("FAIL: no common snapshot step between run and ref"); sys.exit(1)
        s, dn, dref = common[-1]
        delta = dn - dref
        # 合否: 参照が定常、new − ref ≤ allow、かつ new ≤ max(frac·h0scale, ref) (参照が既に絶対水準を超える cell 場では「悪化しない」を要求)
        okd = steady_r and delta <= a.allow and dn <= max(a.frac * a.h0scale, dref + a.steady)   # 参照側の定常揺らぎ (--steady) を許す
        print("vs ref @step %d: new %.3f  ref %.3f  new-ref %+.3f J/kg (allow %.0f; new <= max(%.0f, ref+%.0f)) -> %s" % (s, dn, dref, delta, a.allow, a.frac * a.h0scale, a.steady, "ok" if okd else "EXCEED"))
        fail = fail or (not okd)
    print("VERDICT: %s" % ("FAIL" if fail else "PASS"))
    sys.exit(1 if fail else 0)

if __name__ == "__main__": main()
