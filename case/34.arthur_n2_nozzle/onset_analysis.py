#!/usr/bin/env python3
"""Arthur ノズル凝縮 run の onset 解析 (plan condensation-air §4.4/§6):
  - dry 基準 (restart_dry.h5 の保存量から p_dry = (γ−1)(ρe − ½ρ|u|²)) との中心線 p/p_dry − 1 が 1 % を超える最初の x を onset とし、
    その点の dry 状態 (p, T) と、g>1e-4 の最初の x を出す。
  - 膨張率 Ṗ = −(u/p) dp/dx [1/s] を dry 中心線で評価 (Daum & Gyarmathy の相関パラメータ)。
  - Daum & Gyarmathy 最小 onset 曲線 (daum_gyarmathy_min_onset_n2.csv) の T_onset(p) と比べた ΔT を出す。
usage: onset_analysis.py RUN_DIR [--res res_N.h5] [--dry restart_dry.h5] [--gamma 1.4] [--R 296.8]"""
import sys, os, glob, argparse, h5py, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--res", default=None); ap.add_argument("--dry", default="restart_dry.h5")
ap.add_argument("--gamma", type=float, default=1.4); ap.add_argument("--R", type=float, default=296.8); ap.add_argument("--thr", type=float, default=0.01)
ap.add_argument("--series", action="store_true", help="全 res_*.h5 で報告量を出し末尾窓 (最大 4 枚, 最低 3 枚) の振幅で STEADY/DRIFTING/OSCILLATING/UNSETTLED 判定 (不合格は非ゼロ終了)")
ap.add_argument("--stations", default="3,4,5", help="cond/dry 比を報告する軸位置 [in]")
a = ap.parse_args()
def latest(run):
    fs = [f for f in glob.glob(os.path.join(run, "res_*.h5")) if "nan" not in f]; return max(fs, key=lambda f: int(os.path.basename(f)[4:-3]))
import yaml as _yaml
_meshname = _yaml.safe_load(open(os.path.join(a.run, "solverConfig.yaml")))["mesh"]["meshFileName"]
mesh = h5py.File(os.path.join(a.run, _meshname), "r")
d = h5py.File(os.path.join(a.run, a.dry), "r")["VALUE"]
f_res_coord = None
def analyze(res, quiet=False):
    global f_res_coord
    f = h5py.File(res, "r"); V = f["VALUE"]
    coord = np.array(f["MESH/COORD"]).reshape(-1, 3); f_res_coord = coord
    # node (median-dual): 値の位置はノード座標 (MESH/COORD)。cell: 入力メッシュの CELLS/centCoords。
    c = coord[:len(V["P"])] if len(V["P"]) == len(coord) else np.array(mesh["CELLS/centCoords"]).reshape(-1, 3)[:len(V["P"])]
    return _analyze(res, V, c, quiet)
def _analyze(res, V, c, quiet):
    ro = np.array(d["ro"], float); u = np.array(d["roUx"], float)/ro; v = np.array(d["roUy"], float)/ro; w = np.array(d["roUz"], float)/ro; roe = np.array(d["roe"], float)
    p_dry = (a.gamma - 1)*(roe - 0.5*ro*(u*u + v*v + w*w)); T_dry = p_dry/(ro*a.R)
    P = np.array(V["P"], float); g = np.array(V["g_0"], float) if "g_0" in V else np.zeros_like(P)
    # 中心線: source-flow メッシュはセル中心 x が列ごとに揺れるので x をビン分け (~400 列) し、各ビンで |y| 最小のセルを取る
    nb = 400; edges = np.linspace(c[:, 0].min(), c[:, 0].max() + 1e-9, nb + 1); ib = np.clip(np.digitize(c[:, 0], edges) - 1, 0, nb - 1)
    bins = [np.where(ib == k)[0] for k in range(nb)]; bins = [b for b in bins if len(b)]
    is_node = (len(V["P"]) == len(np.array(f_res_coord)))
    if is_node:
        # node (値=ノード): 中心線は y=0 のノード列、壁は wall_dist=0 かつ y>0 のノード列 (中心線外の内部点を拾わない; codex 2026-09-13 M4)
        C = np.where(np.abs(c[:, 1]) < 1e-9)[0]                      # 中心線ノード (y=0)
        W = np.array([b[np.argmax(c[b, 1])] for b in bins])            # 壁ノード列 (各 x ビンで y 最大 = 輪郭壁上のノード; 全 slip なので wall_dist は使えない)
    else:
        C = np.array([b[np.argmin(np.abs(c[b, 1]))] for b in bins])   # cell: 各 x ビンで |y| 最小 (第一セル列)
        W = np.array([b[np.argmax(c[b, 1])] for b in bins])           # 壁セル列 (y 最大 = 輪郭壁の第一セル; Arthur は上半分メッシュ)
    C = C[np.argsort(c[C, 0])]; W = W[np.argsort(c[W, 0])]
    x = c[C, 0]; pd = p_dry[C]; Td = T_dry[C]; ud = u[C]; ratio = P[C]/pd - 1.0; gc = g[C]
    xw = c[W, 0]; ratio_w = P[W]/p_dry[W]        # 壁 p/p_dry (Arthur が測ったのは壁面静圧)
    Pdot = -ud/pd*np.gradient(pd, x)   # [1/s] (膨張で正)
    def first(mask):
        i = np.where(mask)[0]; return int(i[0]) if len(i) else None
    io = first(ratio > a.thr); ig = first(gc > 1e-4)
    on = np.genfromtxt(os.path.join(a.run, "..", "daum_gyarmathy_min_onset_n2.csv"), delimiter=",", comments="#", skip_header=1)
    def T_dg(p): return np.interp(np.log10(p), np.log10(on[:, 0]), on[:, 1])
    if not quiet: print(f"res={os.path.basename(res)}  dry={a.dry}  cells on centerline={len(C)}")
    for lab, i in (("Δp/p_dry > %.0f %%" % (a.thr*100), io), ("g > 1e-4", ig)):
        if i is None:
            if not quiet: print(f"  onset ({lab}): not reached")
            continue
        if not quiet: print(f"  onset ({lab}): x={x[i]*1e3:.2f} mm = {x[i]/0.0254:.2f} in,  p_dry={pd[i]:.1f} Pa, T_dry={Td[i]:.2f} K, Pdot={Pdot[i]:.3g} 1/s;  D&G min onset T({pd[i]:.0f} Pa)={T_dg(pd[i]):.1f} K -> ΔT = {Td[i]-T_dg(pd[i]):+.1f} K (正=D&G より高温=過冷却が浅い)")
    if not quiet: print(f"  exit: x={x[-1]/0.0254:.2f} in  p_dry={pd[-1]:.1f} Pa T_dry={Td[-1]:.1f} K  g_exit={gc[-1]:.4f}  max(p/p_dry-1)={ratio.max():.3f}  Pdot range along nozzle [{Pdot[np.isfinite(Pdot)].min():.3g}, {Pdot[np.isfinite(Pdot)].max():.3g}]")
    # 報告量 (series 用): onset x [in]/T/p, 各 station の cond/dry 比, g_exit, 出口面 u_n/c の最小 (超音速判定; sonic があれば)
    st = [float(v) for v in a.stations.split(",")]
    ratio_st = {v: float(np.interp(v*0.0254, xw, ratio_w)) for v in st}   # 壁セル列の p/p_dry (codex 2026-09-12 M5: 中心線ではなく壁)
    unc = np.nan
    if "sonic" in V and "Ux" in V:
        # 出口境界に隣接するセル (最後の x ビン) の法線 Mach。Arthur の出口面は x 法線なので u_n=Ux (近似: 出口面の法線ベクトルは使わない)
        ex = bins[np.argmax([c[b, 0].mean() for b in bins])]
        unc = float(np.min(np.array(V["Ux"], float)[ex]/np.array(V["sonic"], float)[ex]))
    return dict(onset_in=(x[io]/0.0254 if io is not None else np.nan), T_on=(Td[io] if io is not None else np.nan), p_on=(pd[io] if io is not None else np.nan),
                Pdot_on=(Pdot[io] if io is not None else np.nan), g_exit=float(gc[-1]), un_c_min=unc, **{f"r{v:g}": ratio_st[v] for v in st})

if a.series:
    fs = sorted(((int(os.path.basename(f)[4:-3]), f) for f in glob.glob(os.path.join(a.run, "res_*.h5")) if "nan" not in f)); fs = [(k, f) for k, f in fs if k > 0]
    st = [float(v) for v in a.stations.split(",")]
    rows = [(k, analyze(f, quiet=True)) for k, f in fs]
    import yaml
    cfgc = yaml.safe_load(open(os.path.join(a.run, "solverConfig.yaml"))).get("condensation", {})
    dry = int(cfgc.get("condensation", 0)) == 0     # dry/cond は config で決める (g_exit=0 から推定しない; codex M5)
    # 閾値 (plan condensation-air §5 8): onset 0.02 in / T_on 0.3 K / p_on 1 % / Pdot 5 % / 壁圧比 0.2 % / g_exit 0.5 %
    # onset は x ビン (400 列, 1 ビン ≈ 0.33 mm ≈ 0.013 in) の離散位置なので、1 ビン分の交差ゆらぎ (Δp/p ≈ Ṗ Δx/u ≈ 0.9 %, ΔT ≈ 0.25 K) を許容:
    #   onset 0.02 in / T_on 0.4 K / p_on 2 % / Pdot 10 % / g_exit 0.5 % / 壁圧比 0.2 %
    TOL = ([] if dry else [("onset_in", "abs", 0.02), ("T_on", "abs", 0.4), ("p_on", "rel", 2.0), ("Pdot_on", "rel", 10.0), ("g_exit", "rel", 0.5)]) + [(f"r{v:g}", "rel", 0.2) for v in st]
    # 出口の超音速判定は全保存時刻で独立の合否 (u_n/c > 1)
    sup_ok = all(np.isfinite(m["un_c_min"]) and m["un_c_min"] > 1.0 for _, m in rows)
    print(f"== {a.run.rstrip('/')} ==  ({len(rows)} snapshots)")
    print("| step | onset [in] | T_on [K] | p_on [Pa] | Pdot [1/s] | " + " | ".join(f"cond/dry @{v:g} in" for v in st) + " | g_exit | min u_n/c (exit) |")
    print("|---|---|---|---|---|" + "---|"*len(st) + "---|---|")
    for k, m in rows:
        print(f"| {k} | {m['onset_in']:.3f} | {m['T_on']:.2f} | {m['p_on']:.1f} | {m['Pdot_on']:.3g} | " + " | ".join(f"{m[f'r{v:g}']:.4f}" for v in st) + f" | {m['g_exit']:.4f} | {m['un_c_min']:.3f} |")
    print(f"exit supersonic (u_n/c>1) at all saved times: {'YES' if sup_ok else 'NO  <-- FAIL'}  (min over snapshots = {min(m['un_c_min'] for _, m in rows):.3f})")
    if len(rows) < 3: print("VERDICT(series): TRANSIENT-UNSETTLED (snapshots < 3)"); sys.exit(1)
    win = rows[-4:]; bad = []; osc = []
    for q, kind, tol in TOL:
        v = np.array([m[q] for _, m in win])
        if not np.all(np.isfinite(v)): print(f"  {q:9s}: NON-FINITE {v} <-- UNEVALUABLE"); osc.append(q + "(non-finite)"); continue
        amp = v.max() - v.min(); lim = tol if kind == "abs" else tol/100*abs(v[-1]); dv = np.diff(v); mono = np.all(dv >= 0) or np.all(dv <= 0)
        flag = "" if amp <= lim else (" DRIFT" if mono else " OSC"); print(f"  {q:9s}: window {win[0][0]}..{win[-1][0]} amp={amp:.4g} (tol {lim:.4g}){flag}")
        if amp > lim: (bad if mono else osc).append(q)
    if not bad and not osc: v = "STEADY"
    elif any(q.endswith("(non-finite)") for q in osc): v = "UNEVALUABLE (" + ",".join(osc + bad) + ")"
    elif osc: v = "OSCILLATING (" + ",".join(osc + bad) + ")"
    else: v = "DRIFTING (" + ",".join(bad) + ")"
    if not sup_ok: v += " + EXIT-NOT-SUPERSONIC"
    print("VERDICT(series):", v); sys.exit(0 if (v == "STEADY") else 1)
else:
    analyze(a.res or latest(a.run))
