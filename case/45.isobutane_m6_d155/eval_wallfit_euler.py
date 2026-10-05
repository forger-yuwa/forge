"""補間壁 (A) vs 同時当てはめ壁 (B) の Euler A/B を固定評価器で判定する。plan verification-m6-axis-wave-mesh-su2 §5.1 #15 (事前登録)。
A: run_0053〜0055 / B: run_0056〜0058 (各腕 3 回の独立再実行)。
延長 (準定常未達時の 6000 step × 1 回, run_0059〜0064) があれば、本段 (step 1000〜12000) に延長 (12000+1000〜) を連結して判定する。
v2 (2026-10-05, codex 6 回目): 分布差の U に 2T・2E を入れ、E ≤ Δ/10 を前提ゲートに、残差 VERDICT は全文の行を保存。登録時の出力 (wallfit_euler_ab.json) は上書きしない。
usage: python3 eval_wallfit_euler.py [case_dir] [--fixed-coef]   (→ _band_ab/wallfit_euler_ab_v2.json / _band_ab/wallfit_euler_ab_diag_fixedcoef.json) → _band_ab/wallfit_euler_ab.json, 各 run の wallfit_series.csv と QUASISTEADY_wallfit.txt
"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
from scipy.interpolate import BSpline
FIXED = "--fixed-coef" in sys.argv   # codex 6 回目の診断 B: 刻み感度 E を、0.05 で当てた平滑化曲線を固定して 0.025 で評価する
_args = [a for a in sys.argv[1:] if not a.startswith("--")]
C = Path(_args[0]).resolve() if _args else Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402

X_E, X_F = 39.82004263, 95.22667765
WIN_T = (41.82004263, 94.22667765); WIN_O = (24.82004263, 95.22667765)
P_REF, MD = 2237.0, 6.0
DELTA = dict(M_wave=0.001, M_resid_diff=0.001, overshoot=0.003, exit_core_M=0.00018, P_wave=0.010, P_resid_diff=0.010)
ARMS = {"A": [f"run_00{n}_euler_wallfit_interp_r{k}" for n, k in ((53, 1), (54, 2), (55, 3))],
        "B": [f"run_00{n}_euler_wallfit_fit_r{k}" for n, k in ((56, 1), (57, 2), (58, 3))]}


def grid(h):
    g = np.arange(0.0, X_F + 1e-9, h)
    return np.unique(np.r_[g, WIN_T, WIN_O, X_E, X_F])


def pspline_fixed(x, v, lam=1.0, k=3, coef=False):
    """3 次 B-spline、内部ノット x = 10k (固定)、2 階差分罰則。フィット区間 [X_E, X_F] 固定。coef=True なら BSpline も返す。"""
    m = (x >= X_E) & (x <= X_F)
    ki = np.arange(10.0, X_F, 10.0); ki = ki[(ki > X_E) & (ki < X_F)]
    t = np.r_[[X_E] * (k + 1), ki, [X_F] * (k + 1)]; nc = len(t) - k - 1
    B = BSpline.design_matrix(x[m], t, k).toarray()
    D = np.diff(np.eye(nc), 2, axis=0)
    c = np.linalg.solve(B.T @ B + lam * D.T @ D, B.T @ v[m])
    out = np.full_like(v, np.nan); out[m] = B @ c
    return (out, BSpline(t, c, k)) if coef else out


def quantities(F, h, fixed_from=None):
    """fixed_from=h0 なら、平滑化曲線を h0 格子で当てはめて固定し、h 格子ではその曲線を評価するだけ (再フィットしない; codex 6 回目の診断 B)。"""
    xq = grid(h); wt = (xq >= WIN_T[0]) & (xq <= WIN_T[1]); wo = (xq >= WIN_O[0]) & (xq <= WIN_O[1])
    q, dist = {}, {}
    for eta in (0.0, 0.1):
        dm = 100 * (eta_line(F, "M", eta, xq) / MD - 1); dp = 100 * (eta_line(F, "P", eta, xq) / P_REF - 1)
        if fixed_from is None:
            rm = dm - pspline_fixed(xq, dm); rp = dp - pspline_fixed(xq, dp)
        else:
            x0 = grid(fixed_from)
            _, sm = pspline_fixed(x0, 100 * (eta_line(F, "M", eta, x0) / MD - 1), coef=True)
            _, sp = pspline_fixed(x0, 100 * (eta_line(F, "P", eta, x0) / P_REF - 1), coef=True)
            rm = dm - sm(xq); rp = dp - sp(xq)
        q[f"M_wave_eta{eta}"] = float(np.abs(rm[wt]).max()); q[f"P_wave_eta{eta}"] = float(np.abs(rp[wt]).max())
        q[f"overshoot_eta{eta}"] = float(dm[wo].max())
        cf = np.polyfit(xq[wt], dp[wt], 1); q[f"P_slope_eta{eta}"] = float(cf[0] * (WIN_T[1] - WIN_T[0]))
        dist[f"M_resid_eta{eta}"] = (xq[wt], rm[wt]); dist[f"P_resid_eta{eta}"] = (xq[wt], rp[wt]); dist[f"P_raw_eta{eta}"] = (xq[wt], dp[wt])
    eta_last = F["R"][-1] / F["R"][-1, -1]; eg = np.linspace(0.05, 0.7, 131)
    Me = np.interp(eg, eta_last, F["V"]["M"][-1])
    q["exit_core_M"] = float(np.trapezoid(Me * eg, eg) / np.trapezoid(eg, eg))
    return q, dist


def snaps(run):
    fs = [f.name for f in (C / run).glob("res_[0-9]*.h5")]
    return sorted(fs, key=lambda f: int(re.findall(r"\d+", f)[0]))


def ext_of(run):
    m = sorted(C.glob("run_00[0-9][0-9]_" + run[9:] + "_ext6k"))
    return m[0] if m else None


def series_files(run):
    """[(run_dir, res 名, 通算 step)]。延長があれば本段 12000 の後に連結 (延長の step 0 は本段 12000 と同じ場なので除く)。"""
    rd = C / run; out = [(rd, f, int(re.findall(r"\d+", f)[0])) for f in snaps(run) if int(re.findall(r"\d+", f)[0]) > 0]
    e = ext_of(run)
    if e is not None:
        out += [(e, f, 12000 + int(re.findall(r"\d+", f)[0])) for f in snaps(e.name) if int(re.findall(r"\d+", f)[0]) > 0]
    return out


R = {}
for arm, runs in ARMS.items():
    for run in runs:
        rd = C / run
        sf = series_files(run)
        series = []
        for d_, f, st_ in sf:
            F = load_field(d_, f); q, _ = quantities(F, 0.05); q["step"] = st_; series.append(q)
        cols = [k for k in series[0] if k != "step"]
        with open(rd / "wallfit_series.csv", "w") as fh:
            fh.write("step," + ",".join(cols) + "\n")
            for q in series:
                fh.write(f"{q['step']}," + ",".join(f"{q[c]:.10g}" for c in cols) + "\n")
        qs = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(rd / "wallfit_series.csv"),
                             "--series-cols", ",".join(cols)], capture_output=True, text=True).stdout
        (rd / "QUASISTEADY_wallfit.txt").write_text(qs)
        verd = {c: re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", qs, re.M).group(1) for c in cols}
        tail = series[-5:]
        rep = {c: float(np.mean([q[c] for q in tail])) for c in cols}
        T = {}
        for c in cols:
            v = np.array([q[c] for q in tail]); st = np.array([q["step"] for q in tail], float)
            T[c] = max(float(np.ptp(v)), float(abs(np.polyfit(st, v, 1)[0]) * (st[-1] - st[0])))
        last = sf[-5:]
        dists = [quantities(load_field(d_, f), 0.05)[1] for d_, f, _ in last]
        dmean = {k: (dists[0][k][0], np.mean([d[k][1] for d in dists], axis=0)) for k in dists[0]}
        qE, _ = quantities(load_field(sf[-1][0], sf[-1][1]), 0.025, fixed_from=0.05 if FIXED else None); q5 = series[-1]
        E = {c: abs(qE[c] - q5[c]) for c in cols if c in qE}
        cv = (sf[-1][0] / "CONVERGENCE_VERDICT.txt")
        conv = cv.read_text() if cv.exists() else "missing"
        conv_overall = [l.strip() for l in conv.splitlines() if "-> " in l and l.startswith("===")]
        R[run] = dict(arm=arm, rep=rep, T=T, E=E, quasisteady=verd, dist=dmean, n_snaps=len(series), last_step=series[-1]["step"],
                      convergence=conv_overall[-1] if conv_overall else conv.strip().splitlines()[-1] if conv != "missing" else "missing",
                      convergence_file=str(cv.relative_to(C)) if cv.exists() else None)
        # 分布の時間変動 T_dist (同じ x ごとの末尾 5 枚の ptp と線形ドリフトの最大) と刻み感度 E_dist (0.025 で同じ曲線を評価し 0.05 の点へ戻した差)
        st5 = np.array([st_ for _, _, st_ in last], float)
        R[run]["T_dist"] = {k: float(max(np.ptp(np.array([d[k][1] for d in dists]), axis=0).max(),
                                         (np.abs(np.polyfit(st5, np.array([d[k][1] for d in dists]), 1)[0]) * (st5[-1] - st5[0])).max()))
                            for k in dists[0]}
        _, dE = quantities(load_field(sf[-1][0], sf[-1][1]), 0.025, fixed_from=0.05 if FIXED else None)
        R[run]["E_dist"] = {k: float(np.abs(np.interp(dists[-1][k][0], dE[k][0], dE[k][1]) - dists[-1][k][1]).max()) for k in dE}

# メッシュ照合: B が使われた証拠 (準備ディレクトリの壁節点)
mesh_check = {}
try:
    W = {}
    for arm in ("interp", "fit"):
        pd = C / f"_prep_wallfit_{arm}"; info = json.loads((pd / "prepare_info.json").read_text())
        with h5py.File(pd / "nozzle.h5") as f:
            nc = f["/MESH/COORD"][:].reshape(-1, 3)
        ni = int(info["mesh"]["ni"]); S = float(info["scale_m"]); nj = nc.shape[0] // ni
        W[arm] = ((nc[:, 0] / S).reshape(ni, nj)[:, -1], (nc[:, 1] / S).reshape(ni, nj)[:, -1])
    xa, ra = W["interp"]; xb, rb = W["fit"]
    mesh_check = dict(same_x=bool(np.allclose(xa, xb, atol=1e-9)), wall_dr_max=float(np.abs(rb - ra).max()),
                      x_at_max=float(xa[np.argmax(np.abs(rb - ra))]))
except Exception as e:  # noqa: BLE001
    mesh_check = {"error": str(e)}

scal = list(next(iter(R.values()))["rep"].keys())
out = {"runs": {k: {kk: vv for kk, vv in v.items() if kk != "dist"} for k, v in R.items()}, "mesh_check": mesh_check, "delta": DELTA}
dkey = lambda c: ("M_wave" if c.startswith("M_wave") else "P_wave" if c.startswith("P_wave") else "overshoot" if c.startswith("overshoot")
                  else "exit_core_M" if c == "exit_core_M" else None)
judge = {}
for c in scal:
    dk = dkey(c)
    a = np.array([R[r]["rep"][c] for r in ARMS["A"]]); b = np.array([R[r]["rep"][c] for r in ARMS["B"]])
    T = max(R[r]["T"][c] for r in ARMS["A"] + ARMS["B"]); E = max(R[r]["E"].get(c, 0.0) for r in ARMS["A"] + ARMS["B"])
    row = dict(A_mean=float(a.mean()), B_mean=float(b.mean()), B_minus_A=float(b.mean() - a.mean()), R_A=float(np.ptp(a)), R_B=float(np.ptp(b)), T=T, E=E)
    if dk:
        D = DELTA[dk]; U = max(3 * row["R_A"], 3 * row["R_B"], 2 * T, 2 * E, D / 10)
        row.update(Delta=D, U=U, tail_ok=bool(T <= D / 4), step_ok=bool(E <= D / 10))
        row["judge"] = ("small" if (U <= D / 2 and abs(row["B_minus_A"]) + U <= D) else "different" if abs(row["B_minus_A"]) - U > D else "hold")
    judge[c] = row
for k in ("M_resid_eta0.0", "M_resid_eta0.1", "P_resid_eta0.0", "P_resid_eta0.1"):
    dk = "M_resid_diff" if k.startswith("M") else "P_resid_diff"; D = DELTA[dk]
    A_ = np.array([R[r]["dist"][k][1] for r in ARMS["A"]]); B_ = np.array([R[r]["dist"][k][1] for r in ARMS["B"]])
    diff = float(np.abs(B_.mean(0) - A_.mean(0)).max()); RA_ = float(np.ptp(A_, 0).max()); RB_ = float(np.ptp(B_, 0).max())
    Td = max(R[r]["T_dist"][k] for r in ARMS["A"] + ARMS["B"]); Ed = max(R[r]["E_dist"][k] for r in ARMS["A"] + ARMS["B"])
    U = max(3 * RA_, 3 * RB_, 2 * Td, 2 * Ed, D / 10)
    judge[k + "_dist"] = dict(B_minus_A_maxabs=diff, R_A=RA_, R_B=RB_, T=Td, E=Ed, Delta=D, U=U, tail_ok=bool(Td <= D / 4), step_ok=bool(Ed <= D / 10),
                              judge=("small" if (U <= D / 2 and diff + U <= D) else "different" if diff - U > D else "hold"))
for k in ("P_raw_eta0.0", "P_raw_eta0.1"):
    A_ = np.array([R[r]["dist"][k][1] for r in ARMS["A"]]); B_ = np.array([R[r]["dist"][k][1] for r in ARMS["B"]])
    judge[k + "_dist"] = dict(B_minus_A_maxabs=float(np.abs(B_.mean(0) - A_.mean(0)).max()), note="併記のみ")
qs_all = {r: all(v == "STEADY" for v in R[r]["quasisteady"].values()) for r in R}
out["judge"] = judge; out["all_steady"] = qs_all
gated = [v["judge"] for v in judge.values() if "judge" in v]
tail_ok = all(v.get("tail_ok", True) and v.get("step_ok", True) for v in judge.values())
out["overall"] = ("保留 (準定常・刻み精度の前提未達)" if not (all(qs_all.values()) and tail_ok) else
                  "支持: 両壁の差は小さい" if all(g == "small" for g in gated) else
                  "棄却: 差あり (向きを確認)" if any(g == "different" for g in gated) else "保留")
(C / "_band_ab").mkdir(exist_ok=True)
(C / ("_band_ab/wallfit_euler_ab_diag_fixedcoef.json" if FIXED else "_band_ab/wallfit_euler_ab_v2.json")).write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps({"mesh_check": mesh_check, "all_steady": qs_all, "overall": out["overall"],
                  "judge": {k: {kk: (round(vv, 7) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in judge.items()}}, indent=1, ensure_ascii=False))
