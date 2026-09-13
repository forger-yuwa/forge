#!/usr/bin/env python3
"""凝縮 run の軸 (対称線) 上の全エンタルピー保存を時系列で判定する (plans/active/condensation-float-speedup.md §4.3 (d), codex plan M6)。

  metric(step) = max_x |h0(x) − h0_in|  [J/kg]   (h0 = VALUE/h0, h0_in = 軸の最上流 3 % の平均)
  VERDICT: 末尾 3 枚の metric の変化 ≤ --steady (既定 10 J/kg) で STEADY、そうでなければ DRIFTING。
  --ref <run> を与えると同じ step のスナップショットで metric の差 (new − ref) を出し、--allow (既定 100 J/kg) 以下と
  両方 ≤ --frac (既定 1e-3)·--h0scale (既定 3e5 J/kg ≈ c_p·300 K) を要求する (twophase face-T fix の 0.1 % 規約。datum 付きの h0 自体は ~0 なので尺度にしない)。

軸ノード: |y| ≤ --ytol (既定 1e-9 m) の節点。無ければ (中心線に節点の無い偶数分割メッシュ) x 列ごとに |y| 最小の節点。--xmin/--xmax で評価範囲 (凝縮帯) を絞る。
使い方: python3 solver_density_cuda/tools/cond_axis_h0.py <run_dir> [--ref <run_dir>] [--xmin X --xmax X]
"""
import argparse, glob, os, re, sys
import h5py, numpy as np

def mesh_file(run):
    cfg = open(os.path.join(run, "solverConfig.yaml")).read()
    m = re.search(r'meshFileName:\s*"?([^"\s]+)"?', cfg)
    return os.path.join(run, m.group(1))

def axis_nodes(run, ytol):
    """軸 (対称線) の節点列: |y| ≤ ytol の節点。無ければ (中心線に節点が無い偶数分割メッシュ) x 列ごとに |y| 最小の節点を取る。"""
    with h5py.File(mesh_file(run), "r") as f:
        C = f["MESH/COORD"][...]
    C = C.reshape(-1, 3) if C.ndim == 1 else C
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
    print("run: %s  axis nodes: %d  x %.4f..%.4f" % (a.run, idx.size, x.min(), x.max()))
    print("%8s %14s %10s %14s" % ("step", "max|h0-h0in|", "at x", "h0_in"))
    for s, d, xx, hin in rows: print("%8d %14.3f %10.5f %14.1f" % (s, d, xx, hin))
    tail = [r[1] for r in rows[-3:]]
    steady = (max(tail) - min(tail)) <= a.steady if len(tail) >= 2 else False
    last = rows[-1]; ok = last[1] <= a.frac * a.h0scale
    verdict = "STEADY" if steady else "DRIFTING"
    print("h0 conservation: metric %.3f J/kg = %.2e of h0scale %.1e (limit %.0e) -> %s; tail change %.3f J/kg (limit %.1f) -> %s"
          % (last[1], last[1] / a.h0scale, a.h0scale, a.frac, "ok" if ok else "EXCEED", max(tail) - min(tail), a.steady, verdict))
    fail = (not ok) or (not steady)
    if a.ref:
        xr, idr = axis_nodes(a.ref, a.ytol); rr = metric(a.ref, xr, idr, a.xmin, a.xmax)
        dr = {s: d for s, d, _, _ in rr}
        common = [(s, d, dr[s]) for s, d, _, _ in rows if s in dr]
        if common:
            s, dn, dref = common[-1]
            delta = dn - dref
            okd = delta <= a.allow and dref <= a.frac * a.h0scale
            print("vs ref %s @step %d: new %.3f  ref %.3f  new-ref %+.3f J/kg (allow %.0f) -> %s" % (a.ref, s, dn, dref, delta, a.allow, "ok" if okd else "EXCEED"))
            fail = fail or (not okd)
    print("VERDICT: %s" % ("FAIL" if fail else "PASS"))
    sys.exit(1 if fail else 0)

if __name__ == "__main__": main()
