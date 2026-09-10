#!/usr/bin/env python3
"""Kantrowitz γ_v / 二相 frozen 音速の修正 (plan condensation-kantrowitz-gamma-twophase-sonic) の A/B 比較。
usage: compare_condfix.py RUN_DIR [RUN_DIR ...] [--step N] [--out PREFIX]
  最初の run を参照にして場の差 (max |Δ|/max|ref|) を出し、壁 p/p0 vs Fig.3 実験・onset・出口諸量を表にする。"""
import sys, glob, os, argparse, h5py, yaml, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for fp in glob.glob("/home/sano/.fonts/*Noto*CJK*"): font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = "Noto Sans CJK JP"
ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--step", type=int, default=None)
ap.add_argument("--out", default="compare_condfix"); ap.add_argument("--p0", type=float, default=59070.0)
ap.add_argument("--series", action="store_true", help="各 run の全 res_*.h5 で報告量の時系列を出し、末尾 2 点の変化で STEADY 判定")
a = ap.parse_args()
exp = np.genfromtxt("wyslouzil_fig3_pp0.csv", delimiter=",", skip_header=1)[:, :3]
xe, iso, cond = exp[:, 0]*10, exp[:, 1], exp[:, 2]

def latest(run):
    fs = glob.glob(os.path.join(run, "res_*.h5")); st = {int(os.path.basename(f)[4:-3]): f for f in fs}
    return st[a.step] if a.step in st else st[max(st)]

def load(run, fn=None):
    fn = fn or latest(run); f = h5py.File(fn, "r"); c = np.array(f["MESH/COORD"]).reshape(-1, 3); V = f["VALUE"]
    d = {k: np.array(V[k]) for k in ("P", "T", "ro", "Ux", "Uy", "wall_dist") if k in V}
    d["cell"] = False
    if len(d["P"]) != len(c):   # cell 離散化: VALUE はセル中心量。入力メッシュ h5 の CELLS/centCoords を座標に使う (codex M4)
        mesh = yaml.safe_load(open(os.path.join(run, "solverConfig.yaml")))["mesh"]["meshFileName"]
        c = np.array(h5py.File(os.path.join(run, mesh), "r")["CELLS/centCoords"]).reshape(-1, 3)[:len(d["P"])]
        d["cell"] = True
    d["g"] = np.array(V["g_0"]) if "g_0" in V else np.zeros(len(c))
    d["sonic"] = np.array(V["sonic"]) if "sonic" in V else None
    d["Y1"] = np.array(V["Y1"]) if "Y1" in V else None
    d["h0"] = np.array(V["h0"]) if "h0" in V else None
    d["xyz"] = c; d["file"] = fn; return d

def lines(d):
    c = d["xyz"]; xs = np.round(c[:, 0], 6); ux = np.unique(xs)
    W, C = [], []
    if d["cell"]:
        # cell: 各 x 列 (y>0) で wall_dist 最小のセル = 輪郭壁の第一セル (extract_wall_pp0_cell.py と同じ定義)
        for xv in ux:
            cand = np.where(xs == xv)[0]
            j = cand[np.argmin(np.abs(c[cand, 1]))]; C.append(j)
            cy = cand[c[cand, 1] > 0]
            if len(cy): W.append(cy[np.argmin(d["wall_dist"][cy])])
        return np.array(W), np.array(C)
    wall = (d["wall_dist"] <= 0) & (c[:, 1] > 0)
    for xv in ux:
        cand = np.where(xs == xv)[0]
        j = cand[np.argmin(np.abs(c[cand, 1]))]; C.append(j)
        cw = cand[wall[cand]]
        if len(cw): W.append(cw[np.argmax(c[cw, 1])])
    return np.array(W), np.array(C)

def onset(x_mm, g, thr):
    # 閾値交差を線形補間 (セル内の移動も拾う; codex M3)
    i = np.where(g > thr)[0]
    if not len(i): return np.nan
    k = i[0]
    if k == 0 or g[k] == g[k-1]: return x_mm[k]
    return x_mm[k-1] + (thr - g[k-1])/(g[k] - g[k-1])*(x_mm[k] - x_mm[k-1])

def metrics(d):
    W, C = lines(d); xw = d["xyz"][W, 0]*1e3; xc = d["xyz"][C, 0]*1e3
    pw = d["P"][W]/a.p0; gc = d["g"][C]
    pi = np.interp(xe, xw, pw); dev = (pi - cond)/cond*100; sel = xe >= 10
    je = C[np.argmax(xc)]
    M = (np.hypot(d["Ux"][je], d["Uy"][je])/d["sonic"][je]) if d["sonic"] is not None else np.nan
    h0err = np.nan
    if d["h0"] is not None:   # 中心線 h0 の入口値からの最大偏差 [kJ/kg] (h0 は 298.15 K 基準の sensible なので相対値は使わない)
        h0c = d["h0"][C]; h0in = h0c[np.argmin(xc)]; h0err = np.max(np.abs(h0c - h0in))/1e3
    return dict(onset=onset(xc, gc, 1e-3), dev=dev[sel].mean(), dev21=np.interp(21.0, xe, dev), dev42=np.interp(42.0, xe, dev),
                dev52=np.interp(52.0, xe, dev), g_exit=gc[np.argmax(xc)], M_exit=M,
                c_exit=(d["sonic"][je] if d["sonic"] is not None else np.nan), h0err=h0err)

if a.series:
    # 末尾窓 (最大 4 枚, 最低 3 枚) で各報告量の振幅 (max−min) と単調性を見る (codex M3)。
    #   全量が許容内 → STEADY / 単調で超過 → DRIFTING / 非単調で超過 → OSCILLATING / 枚数不足 → TRANSIENT-UNSETTLED。
    # 許容値: (量, 絶対 or 相対 [%], 値)
    TOL = [("onset", "abs", 0.1), ("dev", "abs", 0.1), ("dev21", "abs", 0.1), ("dev42", "abs", 0.1), ("dev52", "abs", 0.1),
           ("g_exit", "rel", 0.5), ("M_exit", "rel", 0.1), ("c_exit", "rel", 0.1), ("h0err", "abs", 0.01)]
    UNIT = dict(onset="mm", dev="%pt", dev21="%pt", dev42="%pt", dev52="%pt", g_exit="", M_exit="", c_exit="m/s", h0err="kJ/kg")
    allok = True
    for run in a.runs:
        fs = sorted(((int(os.path.basename(f)[4:-3]), f) for f in glob.glob(os.path.join(run, "res_*.h5"))))
        fs = [(st, f) for st, f in fs if st > 0]
        print(f"\n== {run.rstrip('/')} ==")
        print("| step | onset g=1e-3 [mm] | 壁偏差 x≥10 平均 [%] | @21 | @42 | @52 | g_exit | M_exit | c_exit [m/s] | h0 中心線 max 偏差 [kJ/kg] |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        ser = []
        for st, f in fs:
            m = metrics(load(run, f)); ser.append((st, m))
            print(f"| {st} | {m['onset']:.3f} | {m['dev']:+.3f} | {m['dev21']:+.2f} | {m['dev42']:+.2f} | {m['dev52']:+.2f} | {m['g_exit']:.5f} | {m['M_exit']:.4f} | {m['c_exit']:.2f} | {m['h0err']:.4f} |")
        if len(ser) < 3:
            print("VERDICT(series): TRANSIENT-UNSETTLED (snapshots < 3)"); allok = False; continue
        win = ser[-4:]; steps = [s for s, _ in win]
        bad = []; mono_bad = []
        for q, kind, tol in TOL:
            v = np.array([m[q] for _, m in win])
            if not np.all(np.isfinite(v)):   # 欠落/NaN/Inf は評価不能 = 不合格 (codex result M4)
                print(f"  {q:7s}: window {steps[0]}..{steps[-1]} NON-FINITE {v}  <-- UNEVALUABLE"); mono_bad.append(q + "(non-finite)"); continue
            amp = v.max() - v.min(); lim = tol if kind == "abs" else tol/100*abs(v[-1])
            dv = np.diff(v); mono = np.all(dv >= 0) or np.all(dv <= 0)
            flag = "" if amp <= lim else (" DRIFT" if mono else " OSC")
            print(f"  {q:7s}: window {steps[0]}..{steps[-1]} amp={amp:.4g} {UNIT[q]} (tol {lim:.4g}){flag}")
            if amp > lim: (bad if mono else mono_bad).append(q)
        if not bad and not mono_bad: v = "STEADY"
        elif any(q.endswith("(non-finite)") for q in mono_bad): v = "UNEVALUABLE (" + ",".join(mono_bad + bad) + ")"
        elif mono_bad: v = "OSCILLATING (" + ",".join(mono_bad + bad) + ")"
        else: v = "DRIFTING (" + ",".join(bad) + ")"
        print(f"VERDICT(series): {v}")
        if v != "STEADY": allok = False
    sys.exit(0 if allok else 1)

runs = [(r.rstrip("/"), load(r)) for r in a.runs]
ref = runs[0][1]
rows = []; fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
ax[0].plot(xe, cond, "ks", ms=6, label="実験 凝縮 1.00 kPa (Fig.3)"); ax[0].plot(xe, iso, "k.", ms=4, alpha=0.4, label="実験 dry isentrope")
for name, d in runs:
    W, C = lines(d); xw = d["xyz"][W, 0]*1e3; xc = d["xyz"][C, 0]*1e3
    pw = d["P"][W]/a.p0; gw = d["g"][W]; gc = d["g"][C]
    m = (xw >= -1) & (xw <= 95); lab = os.path.basename(name)
    ax[0].plot(xw[m], pw[m], "-", lw=1.3, label=lab); ax[1].plot(xc, gc, "-", lw=1.3, label=lab)
    ax[2].plot(xc, d["T"][C], "-", lw=1.3, label=lab)
    pi = np.interp(xe, xw, pw); dev = (pi - cond)/cond*100; sel = xe >= 10
    dv = {xx: np.interp(xx, xe, dev) for xx in (21.0, 42.0, 52.0)}
    je = C[np.argmax(xc)]  # 出口中心
    M = (np.hypot(d["Ux"][je], d["Uy"][je])/d["sonic"][je]) if d["sonic"] is not None else np.nan
    diff = {}
    for k in ("ro", "P", "T", "Ux", "g"):
        if k in d and k in ref and len(d[k]) == len(ref[k]):
            den = np.max(np.abs(ref[k])); diff[k] = np.max(np.abs(d[k] - ref[k]))/den if den > 0 else 0.0
    rows.append(dict(run=lab, file=os.path.basename(d["file"]),
                     onset_c_1e4=onset(xc, gc, 1e-4), onset_c_1e3=onset(xc, gc, 1e-3), onset_w_1e3=onset(xw, gw, 1e-3),
                     dev_mean_x10=dev[sel].mean(), dev21=dv[21.0], dev42=dv[42.0], dev52=dv[52.0],
                     g_exit=gc[np.argmax(xc)], T_exit=d["T"][je], c_exit=(d["sonic"][je] if d["sonic"] is not None else np.nan), M_exit=M,
                     pp0_exit=d["P"][je]/a.p0, gmax=d["g"].max(), Tmin=d["T"].min(), nan=int(np.isnan(d["P"]).sum()), diff=diff))
ax[0].set(xlabel="x [mm]", ylabel="p_wall/p0", title="輪郭壁 p/p0 vs 実験", xlim=(-1, 95)); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7)
ax[1].set(xlabel="x [mm]", ylabel="g (中心線)", title="液相質量分率 (onset)"); ax[1].grid(alpha=.3); ax[1].legend(fontsize=7)
ax[2].set(xlabel="x [mm]", ylabel="T [K] (中心線)", title="中心線 静温"); ax[2].grid(alpha=.3); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(a.out + ".png", dpi=130)
hdr = "| run | res | onset c g>1e-4 [mm] | onset c g>1e-3 | onset wall g>1e-3 | 壁 p/p0 偏差 x≥10 平均 [%] | @21 | @42 | @52 | g_exit | T_exit [K] | c_exit [m/s] | M_exit | p/p0 exit | g_max | T_min | NaN | max|Δ|/max|ref| vs 1st (ro/P/T/Ux/g) |"
out = [hdr, "|" + "---|"*19]
for r in rows:
    dd = "/".join(f"{r['diff'].get(k, np.nan):.1e}" for k in ("ro", "P", "T", "Ux", "g"))
    out.append(f"| {r['run']} | {r['file']} | {r['onset_c_1e4']:.2f} | {r['onset_c_1e3']:.2f} | {r['onset_w_1e3']:.2f} | {r['dev_mean_x10']:+.2f} | {r['dev21']:+.1f} | {r['dev42']:+.1f} | {r['dev52']:+.1f} | {r['g_exit']:.5f} | {r['T_exit']:.2f} | {r['c_exit']:.2f} | {r['M_exit']:.4f} | {r['pp0_exit']:.5f} | {r['gmax']:.5f} | {r['Tmin']:.2f} | {r['nan']} | {dd} |")
txt = "\n".join(out); print(txt); open(a.out + ".txt", "w").write(txt + "\n")
