#!/usr/bin/env python3
"""カウル衝撃足の再同定と flag0/flag1 比較 (plan convection-slau-wall-normal-chi-usage-rule §4.4 / §4.2 / §6)。

**結果を見る前に書いて commit したもの** (plan §5.1 #3)。手順・閾値はすべて plan §4.4 の登録どおりで、ここで変えない。
旧抽出 (v3sern_foot_perrun.py / v3sern_series.py) は max|dp/dx| をランプ全域から取り圧縮・膨張を区別せず、
x≈0 のランプ膨張角を拾っていた。また CELLS/centCoords を使っていた (node の値位置はノード座標)。

手順 (plan §4.4):
  1. 探索区間: ランプ壁ノード (physID 4) のうち x ∈ [x_lip, L_ramp]。x_lip = カウル (physID 5 = cowl_in) の最大 x、
     y_lip はそのノードの y。L_ramp = ランプ壁ノードの最大 x。
  2. 候補: dp_w/dx > 0 の局所最大すべて。上昇比 r = p(x_c + 2t) / p(x_c − 2t) (壁ノード線形補間)、登録条件 r ≥ 1.05。
  3. 唇由来の識別: 候補 (x_c, y_w) → 唇 (x_lip, y_lip) の線分上 50 点で、場の ∇p を線分法線へ射影した符号が 80 % 以上同符号。
     ∇p·n は法線方向の中心差分 (刻み DELTA = t/4)、場の P は節点の線形補間 (scipy LinearNDInterpolator)。
  4. 分岐表 (上から優先) は plan §4.4 の表。
比較 (plan §4.2): 共通固定窓 W = [x_f − 2t, x_f + 5t] (x_f = flag0 末尾 dump 群の識別候補位置の中央値)、共通格子 = W 内の
  ランプ壁ノード (同一メッシュ)、台形則の重み。末尾 40 % の全 dump 対で l_ij = ||p1_j − p0_i||_w / ||p̄0||_w。
  帯内: max l_ij ≤ 0.01 / 差が残る: min l_ij > 0.01 / 他: 判定不能。位置差はスカラー差区間 (許容 = W 内ノード間隔の中央値)。
準定常: 各 run の系列 (x_f, 窓平均 p, p(x_f − t), p(x_f + t), p(x_f + 3t); 座標は flag0 の x_f 中央値に固定) を
  check_quasisteady.py --series-csv --drift 0.002 --osc 0.005 --tail 0.4 で判定。STEADY / OSCILLATING のときだけ比較する。

usage: v3sern_shockfoot.py RUN0 RUN1 --out-dir DIR [--tail 0.4] [--t-cowl 0.002] [--label m6_on]
"""
import argparse, glob, json, re, subprocess, sys
from pathlib import Path
import h5py, numpy as np
from scipy.interpolate import LinearNDInterpolator

RAMP_ID, COWL_IN_ID = 4, 5
R_MIN, LIP_FRAC, N_LINE = 1.05, 0.80, 50
L2_TOL = 0.01
QS_DRIFT, QS_OSC = 0.002, 0.005
TOOLS = Path(__file__).resolve().parents[2] / "solver_density_cuda" / "tools"


def dumps(d):
    out = {}
    for f in glob.glob(str(Path(d) / "res_[0-9]*.h5")):
        m = re.fullmatch(r"res_(\d+)\.h5", Path(f).name)
        if m:
            out[int(m.group(1))] = f
    return {k: out[k] for k in sorted(out) if k > 0}


def bc_nodes(f, pid):
    ids = np.asarray(f["BCONDS"][str(pid)]["iCells"]).ravel()
    return np.unique(ids[ids >= 0])


def geometry(mesh):
    with h5py.File(mesh) as f:
        X = np.asarray(f["MESH/COORD"], np.float64).reshape(-1, 3)
        ramp = bc_nodes(f, RAMP_ID)
        cowl = bc_nodes(f, COWL_IN_ID)
    ramp = ramp[np.argsort(X[ramp, 0])]
    j = cowl[np.argmax(X[cowl, 0])]
    return X, ramp, float(X[j, 0]), float(X[j, 1]), float(X[ramp, 0].max())


def lip_test(interp, xc, yc, xl, yl, delta):
    d = np.array([xl - xc, yl - yc]); L = np.hypot(*d)
    if L == 0:
        return 0.0
    d /= L
    n = np.array([d[1], -d[0]])
    if n[0] < 0:                                   # 法線は下流 (+x) 向きに揃える (符号の一致だけを見るので向きは任意)
        n = -n
    s = np.linspace(0.02, 0.98, N_LINE)[:, None]
    pts = np.array([xc, yc]) + s * L * d
    g = (interp(pts + delta * n) - interp(pts - delta * n)) / (2 * delta)
    g = g[np.isfinite(g)]
    if g.size < 0.9 * N_LINE:                      # 有効点が 45/50 未満なら判定不能 = 不成立 (codex result m5)
        return 0.0
    return max((g > 0).mean(), (g < 0).mean())


def detect(P, X, ramp, x_lip, y_lip, L_ramp, t, xy_interp_nodes):
    """1 dump の識別済み候補 [(x_c, r)] を r 降順で返す。"""
    x = X[ramp, 0]; p = P[ramp]
    sel = (x >= x_lip) & (x <= L_ramp)
    xs, ps = x[sel], p[sel]
    if xs.size < 5:
        return []
    dp = np.gradient(ps, xs)
    cand = [k for k in range(1, len(xs) - 1) if dp[k] > 0 and dp[k] >= dp[k - 1] and dp[k] >= dp[k + 1]]
    interp = None
    out = []
    for k in cand:
        xc = xs[k]
        r = np.interp(xc + 2 * t, x, p) / np.interp(xc - 2 * t, x, p)
        if r < R_MIN:
            continue
        if interp is None:
            interp = LinearNDInterpolator(X[xy_interp_nodes, :2], P[xy_interp_nodes])
        yc = np.interp(xc, x, X[ramp, 1])
        frac = lip_test(interp, xc, yc, x_lip, y_lip, t / 4)
        if frac >= LIP_FRAC:
            out.append((float(xc), float(r), float(frac)))
    return sorted(out, key=lambda c: -c[1])


def quasisteady(csv, cols):
    r = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(csv),
                        "--series-cols", ",".join(cols), "--drift", str(QS_DRIFT), "--osc", str(QS_OSC),
                        "--tail", "1.0"], capture_output=True, text=True)   # 系列は既に末尾窓 (二重に切らない、codex result m5)
    txt = r.stdout + r.stderr
    labels = dict(re.findall(r"^\s*(\S+)\s*:.*?\s(STEADY|OSCILLATING|DRIFTING|TRANSIENT-UNSETTLED|NONFINITE)\s*$", txt, re.M))
    return txt, labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run0"); ap.add_argument("run1")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tail", type=float, default=0.4)
    ap.add_argument("--t-cowl", type=float, default=0.002)
    ap.add_argument("--label", default="m6_on")
    ap.add_argument("--mesh", default=None)
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    t = a.t_cowl
    X, ramp, x_lip, y_lip, L_ramp = geometry(a.mesh or str(Path(a.run0) / "sern.h5"))
    nodes_all = np.arange(len(X))
    res = {"label": a.label, "x_lip": x_lip, "y_lip": y_lip, "L_ramp": L_ramp, "t": t,
           "registered": {"R_MIN": R_MIN, "LIP_FRAC": LIP_FRAC, "N_LINE": N_LINE, "L2_TOL": L2_TOL,
                          "QS": [QS_DRIFT, QS_OSC], "tail": a.tail}}
    runs = {"flag0": a.run0, "flag1": a.run1}
    P, det = {}, {}
    for tag, rd in runs.items():
        d = dumps(rd); steps = list(d)
        ntail = max(3, int(round(len(steps) * a.tail)))
        steps = steps[-ntail:]
        P[tag] = {}
        det[tag] = {}
        for s in steps:
            with h5py.File(d[s]) as g:
                P[tag][s] = np.asarray(g["/VALUE/P"], np.float64)
            det[tag][s] = detect(P[tag][s], X, ramp, x_lip, y_lip, L_ramp, t, nodes_all)
    def summarize(tag):
        found = [c[0][0] for c in det[tag].values() if c]
        multi = sum(1 for c in det[tag].values() if len(c) > 1)
        return {"n_dumps": len(det[tag]), "n_found": len(found), "n_multi": multi,
                "x_c_median": float(np.median(found)) if found else None,
                "candidates_last": det[tag][max(det[tag])]}
    s0, s1 = summarize("flag0"), summarize("flag1")
    res["detect"] = {"flag0": s0, "flag1": s1}
    xr = X[ramp, 0]
    # 分岐表 (plan §4.4、上から優先)
    # run の「検出」= 末尾窓のどれかの dump で識別済み候補がある (登録表に無い「半数以上」条件は削除、codex result m5)。
    # 未検出の主張が最も強くなる側 (検出を最も広く取る) に揃えた。
    half0 = s0["n_found"] > 0
    half1 = s1["n_found"] > 0
    if not half0 and not half1:
        res["verdict"] = f"登録条件では同定できない ({a.label})"
        res["next"] = "m6_on なら m4_off の flag 0/1 対で再実行。m4_off でも未検出なら 判定保留 (衝撃が当たらないとは書かない)"
    elif half0 != half1:
        res["verdict"] = "判定不能 (片側のみ検出)"
    else:
        xf = s0["x_c_median"]
        win = (xr >= xf - 2 * t) & (xr <= xf + 5 * t)
        wx = xr[win]
        spacing = float(np.median(np.diff(wx))) if wx.size > 1 else float("nan")
        res["x_f"] = xf; res["window_nodes"] = int(win.sum()); res["local_spacing"] = spacing
        full_len = 7 * t
        cov = (min(xf + 5 * t, xr.max()) - max(xf - 2 * t, xr.min())) / full_len
        if abs(s1["x_c_median"] - xf) > spacing:
            res["verdict"] = "判定不能 (候補不一致)"
        elif cov < 0.5:
            res["verdict"] = "判定保留 (窓がランプ端を越え残り 50 % 未満)"
        else:
            # 準定常 (固定座標)
            qs = {}
            for tag in runs:
                rows = []
                for s, Pd in P[tag].items():
                    p = Pd[ramp]
                    xc = det[tag][s][0][0] if det[tag][s] else float("nan")
                    rows.append((s, xc, p[win].mean(), *[np.interp(xf + c * t, xr, p) for c in (-1, 1, 3)]))
                csv = out / f"{a.label}_{tag}_series.csv"
                with open(csv, "w") as fcsv:
                    fcsv.write("step,xfoot,pwin_mean,p_m1t,p_p1t,p_p3t\n")
                    for r_ in rows:
                        fcsv.write(",".join(f"{v:.9g}" for v in r_) + "\n")
                txt, lab = quasisteady(csv, ["xfoot", "pwin_mean", "p_m1t", "p_p1t", "p_p3t"])
                (out / f"{a.label}_{tag}_QS.txt").write_text(txt)
                qs[tag] = lab
            res["quasisteady"] = qs
            ok = all(v in ("STEADY", "OSCILLATING") for lab in qs.values() for v in lab.values()) and all(len(l) == 5 for l in qs.values())
            if not ok:
                res["verdict"] = "判定保留 (準定常でない量がある)"
            else:
                w = np.zeros(wx.size)
                dx = np.diff(wx); w[:-1] += dx / 2; w[1:] += dx / 2
                pw0 = [Pd[ramp][win] for Pd in P["flag0"].values()]
                pw1 = [Pd[ramp][win] for Pd in P["flag1"].values()]
                pbar = np.mean(pw0, axis=0)
                nrm = np.sqrt(np.sum(w * pbar ** 2))
                l = np.array([[np.sqrt(np.sum(w * (b - c) ** 2)) / nrm for b in pw1] for c in pw0])
                x0 = [c[0][0] for c in det["flag0"].values() if c]; x1 = [c[0][0] for c in det["flag1"].values() if c]
                dlo, dhi = min(x1) - max(x0), max(x1) - min(x0)
                l_v = "帯内" if l.max() <= L2_TOL else ("差が残る" if l.min() > L2_TOL else "判定不能")
                x_v = "帯内" if (-spacing <= dlo and dhi <= spacing) else ("差が残る" if (dlo > spacing or dhi < -spacing) else "判定不能")
                res["compare"] = {"l_min": float(l.min()), "l_max": float(l.max()), "l_verdict": l_v,
                                  "dx_interval": [dlo, dhi], "dx_tol": spacing, "dx_verdict": x_v}
                order = {"帯内": 0, "判定不能": 1, "差が残る": 2}
                res["verdict"] = max((l_v, x_v), key=lambda v: order[v])
    (out / f"{a.label}_VERDICT.json").write_text(json.dumps(res, ensure_ascii=False, indent=2))
    print(json.dumps({k: res[k] for k in res if k not in ("registered",)}, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
