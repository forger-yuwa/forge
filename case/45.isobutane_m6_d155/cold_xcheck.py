"""冷却ノズルの forge と SU2 の照合 (plan tooling-nozzle-isothermal-wall-chain §5.1 #26) の後処理。両コードの場から同じ式で同じ量を作る。
SU2 は forge と同じ節点 (格子の節点の順番 = msh の順番 = forge の node 番号) を解くので、場を同じ (ni, nj) の構造に並べて扱う。

量 (各スナップショット):
  - δ_loc・θ_r: 一定 x の断面 (各行 j を x 方向に補間) で、帯の外縁 y_b(x) の値 (ρ_e, u_e) を縁として面積で等価な厚さ
    (cold_pair.theta_diag と同じ式)。y_b = 各 x で TP の断熱 (run_0181) と 300 K (run_0183) の抽出の band_y_b の大きい方、感度に ×1.25。
  - 壁の τ_w・q_w: 壁の法線に沿った壁節点・第 1・第 2 内部節点から 2 次の片側差分で ∂u_t/∂n・∂T/∂n、μ は Sutherland
    (1.716e-5/273/111)、λ = μ c_p/Pr (Pr 0.72)。q_w は流体から壁へ向かう向きを正。Q_w = Σ q_w 2π r ds (台形)。
  - 入口・出口の質量流量と全エンタルピー流量 (h0 = c_p T + ½|u|²、k は含めない): 列 i = 0 と i = ni−1 で 2π∫ρu_x (·) r dr。
  - x = 40・70・94 の断面の ρu・T の分布。
CPG の定数は γ 1.27354、c_p 1360、R = c_p(γ−1)/γ。TP の run を読むと h0 は CPG の式になるので、収支は CPG の run だけで見る。

usage:
  python3 cold_xcheck.py reduce-forge <run> [--last 5]      (AWS: 最後の 5 枚の res → _band_ab/cold_pair/xcheck_<run>.npz)
  python3 cold_xcheck.py reduce-su2 <run> [--last 5]        (手元: 最後の 5 つの restart_flow_<iter>.csv → 同上)
  python3 cold_xcheck.py gates-forge <run>                   (AWS: NaN の全件検査 + check_convergence --segment → xcheck_gates_<run>.json)
  python3 cold_xcheck.py gates-su2 <run>                     (手元: history.csv を forge の残差 CSV の列名に変換して check_convergence、
                                                              履歴と restart の CSV の非有限値 → xcheck_gates_<run>.json)
  python3 cold_xcheck.py compare                             (手元: 判定 → _band_ab/cold_pair/xcheck.json)
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUTD = HERE / "_band_ab" / "cold_pair"
GAM, CP, PR = 1.27354, 1360.0, 0.72
RGAS = CP * (GAM - 1.0) / GAM
XS = (40.0, 70.0, 94.0)
XE = np.arange(40.0, 94.0 + 1e-9, 0.25)


def mu_suth(T):
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


def mesh_info(n_nodes):
    z = np.load(OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")
    ni = len(z["x"]); S = float(z["scale"])
    if n_nodes % ni:
        raise SystemExit(f"節点数 {n_nodes} が ni {ni} で割り切れない (格子が違う)")
    return ni, n_nodes // ni, S


def common_yb():
    a = np.load(OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz"); b = np.load(OUTD / "theta_run_0183_ns_coldmesh_tw300_ext_100000.npz")
    x = np.asarray(a["x"]); yb = np.maximum(np.asarray(a["band_y_b"]), np.interp(x, b["x"], b["band_y_b"]))
    return x, yb


def reduce_fields(xy, ro, ux, uy, T, k, h0_with_k, ni, nj, S, yb_x, yb):
    x = xy[:, 0].reshape(ni, nj) / S; r = xy[:, 1].reshape(ni, nj) / S
    RO = ro.reshape(ni, nj); UX = ux.reshape(ni, nj); UY = uy.reshape(ni, nj); TT = T.reshape(ni, nj); KK = k.reshape(ni, nj)
    xt = x[:, 0].copy()
    R = np.empty_like(r); Q = {"ro": np.empty_like(RO), "ux": np.empty_like(UX), "T": np.empty_like(TT)}
    for j in range(nj):
        R[:, j] = np.interp(xt, x[:, j], r[:, j])
        Q["ro"][:, j] = np.interp(xt, x[:, j], RO[:, j]); Q["ux"][:, j] = np.interp(xt, x[:, j], UX[:, j]); Q["T"][:, j] = np.interp(xt, x[:, j], TT[:, j])
    out = {"x": xt}
    yb_i = np.interp(xt, yb_x, yb)
    for fac, tag in ((1.0, ""), (1.25, "_s125")):
        th = np.full(ni, np.nan); dl = np.full(ni, np.nan); qe = np.full(ni, np.nan)
        for i in range(ni):
            rr = R[i]; rw = rr[-1]; rb = rw - fac * yb_i[i]
            if not np.isfinite(rb) or rb <= rr[0]:
                continue
            rf = np.linspace(rb, rw, 4001)
            rho = np.interp(rf, rr, Q["ro"][i]); u = np.interp(rf, rr, Q["ux"][i])
            re_, ue_ = rho[0], u[0]
            qm = np.trapezoid((re_ * ue_ - rho * u) * rf, rf); qq = np.trapezoid(rho * u * (ue_ - u) * rf, rf)
            dl[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qm / (re_ * ue_), 0.0))
            th[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qq / (re_ * ue_ ** 2), 0.0))
            qe[i] = 0.5 * re_ * ue_ ** 2
        out["delta_loc" + tag] = dl; out["theta_r" + tag] = th; out["qdyn_e" + tag] = qe
    # 壁の τ_w・q_w (実際の節点、壁法線の線)
    P0 = np.stack([x[:, -1], r[:, -1]], 1) * S; P1 = np.stack([x[:, -2], r[:, -2]], 1) * S; P2 = np.stack([x[:, -3], r[:, -3]], 1) * S
    d1 = np.linalg.norm(P1 - P0, axis=1); d2 = np.linalg.norm(P2 - P0, axis=1)
    tv = np.gradient(P0, axis=0); tv /= np.linalg.norm(tv, axis=1)[:, None]
    ut = lambda k: UX[:, k] * tv[:, 0] + UY[:, k] * tv[:, 1]  # noqa: E731
    c0 = -(d1 + d2) / (d1 * d2); c1 = d2 / (d1 * (d2 - d1)); c2 = -d1 / (d2 * (d2 - d1))
    dudn = c0 * ut(-1) + c1 * ut(-2) + c2 * ut(-3)
    dTdn = c0 * TT[:, -1] + c1 * TT[:, -2] + c2 * TT[:, -3]
    Tw = TT[:, -1]; muw = mu_suth(Tw)
    out["tau_w"] = muw * np.abs(dudn); out["q_w"] = muw * CP / PR * dTdn; out["T_w"] = Tw; out["x_w"] = x[:, -1]
    ok = np.isfinite(out["qdyn_e"])                       # C_f は帯の外縁 (y_b) の ρ_e u_e²/2 で割る (壁の x へ補間)
    out["c_f"] = out["tau_w"] / np.interp(x[:, -1], xt[ok], out["qdyn_e"][ok], left=np.nan, right=np.nan)
    ds = np.linalg.norm(np.diff(P0, axis=0), axis=1); rwm = P0[:, 1]
    fq = 2.0 * np.pi * out["q_w"] * rwm
    out["Q_w"] = np.array(float(np.sum(0.5 * (fq[1:] + fq[:-1]) * ds)))
    # 入口・出口の流量: 列の折れ線を通る流束 2π∫(F_x dr − F_r dx) r (出口の列は壁法線の層で一定 x でない)
    for i, tag in ((0, "in"), (ni - 1, "out")):
        xx = x[i] * S; rr = r[i] * S
        h0 = CP * TT[i] + 0.5 * (UX[i] ** 2 + UY[i] ** 2) + (KK[i] if h0_with_k else 0.0)
        fl = lambda q: 2 * np.pi * float(np.sum(0.5 * ((RO[i] * q * UX[i] * rr)[1:] + (RO[i] * q * UX[i] * rr)[:-1]) * np.diff(rr)  # noqa: E731
                                                - 0.5 * ((RO[i] * q * UY[i] * rr)[1:] + (RO[i] * q * UY[i] * rr)[:-1]) * np.diff(xx)))
        out["mdot_" + tag] = np.array(fl(1.0)); out["Hdot_" + tag] = np.array(fl(h0)); out["Kdot_" + tag] = np.array(fl(KK[i]))
    for xs in XS:
        i = int(np.argmin(np.abs(xt - xs)))
        out[f"prof_r_{int(xs)}"] = R[i]; out[f"prof_rhou_{int(xs)}"] = Q["ro"][i] * Q["ux"][i]; out[f"prof_T_{int(xs)}"] = Q["T"][i]
    return out


def load_forge(run, step):
    import h5py
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.array(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
    with h5py.File(run / f"res_{step}.h5", "r") as h:
        return xy, *(np.array(h["VALUE/" + k], dtype=float) for k in ("ro", "Ux", "Uy", "T", "k"))


def load_su2(run, it):
    A = np.loadtxt(run / f"restart_flow_{it:06d}.csv", delimiter=",", skiprows=1)
    xy = A[:, 1:3]; ro = A[:, 3]; ux = A[:, 4] / ro; uy = A[:, 5] / ro
    e = A[:, 6] / ro - 0.5 * (ux ** 2 + uy ** 2) - A[:, 7]            # SU2 の Energy は k を含む
    return xy, ro, ux, uy, e * (GAM - 1.0) / RGAS, A[:, 7]


SU2_COLS = {"rms[Rho]": "rms_ro", "rms[RhoU]": "rms_roUx", "rms[RhoV]": "rms_roUy", "rms[RhoE]": "rms_roe", "rms[k]": "rms_roK", "rms[w]": "rms_roOmega"}


def _conv(target: str, extra=()) -> dict:
    import subprocess
    tools = HERE.parents[1] / "solver_density_cuda/tools"
    cc = subprocess.run([sys.executable, str(tools / "check_convergence.py"), target, *extra], capture_output=True, text=True)
    lines = cc.stdout.splitlines()
    verdict = [l for l in lines if "->" in l]
    return {"cmd": f"check_convergence.py {target} {' '.join(extra)}".strip(), "rc": cc.returncode, "verdict": verdict[-1] if verdict else "(判定行なし)",
            "rising": [l.strip() for l in lines if "RISING" in l], "diverged": any("DIVERGED" in l for l in lines), "tail": lines[-12:]}


def gates_forge(run: Path) -> Path:
    """AWS: 全段の残差と全スナップショットの NaN 検査 (ns_n012.nan_scan) と、本段区間の check_convergence。"""
    import ns_n012 as NS
    nan = NS.nan_scan(run)
    NS.jdump(run / "NAN_SCAN.json", nan)
    out = {"run": run.name, "kind": "forge", "nan_verdict": nan.get("VERDICT"), "nan_first": nan.get("first_nonfinite"),
           "convergence": _conv(str(run), ("--segment",))}
    p = OUTD / f"xcheck_gates_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return p


def gates_su2(run: Path) -> Path:
    """手元: SU2 の history.csv (log10 の rms) を forge の残差 CSV の列名・線形値に直して check_convergence にかける (2D なので
    rms_roUz は 0)。履歴の非有限値と、全 restart の CSV の非有限値も数える。"""
    import csv
    with open(run / "history.csv") as fh:
        rd = csv.reader(fh); head = [h.strip().strip('"') for h in next(rd)]; rows = [r for r in rd if r]
    ic = {h: i for i, h in enumerate(head)}
    it = [int(float(r[ic["Inner_Iter"]])) for r in rows]
    cols = {dst: np.array([float(r[ic[src]]) for r in rows]) for src, dst in SU2_COLS.items()}
    nonfin_hist = int(sum(np.count_nonzero(~np.isfinite(v)) for v in cols.values()))
    conv_csv = run / "residual_history_su2.csv"
    with open(conv_csv, "w") as fh:
        fh.write("step," + ",".join(["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]) + "\n")
        for n in range(len(it)):
            v = [10.0 ** cols[c][n] for c in ("rms_ro", "rms_roUx", "rms_roUy")] + [0.0] + [10.0 ** cols[c][n] for c in ("rms_roe", "rms_roK", "rms_roOmega")]
            fh.write(f"{it[n]}," + ",".join(f"{x:.9e}" for x in v) + "\n")
    nf = {}
    for f in sorted(run.glob("restart_flow_[0-9]*.csv")):
        A = np.loadtxt(f, delimiter=",", skiprows=1)
        nf[f.name] = int(np.count_nonzero(~np.isfinite(A)))
    itmax = None
    for line in (run / "sst.cfg").read_text().splitlines():
        if line.strip().startswith("ITER="):
            itmax = int(line.split("=")[1])
    out = {"run": run.name, "kind": "su2", "last_iter": it[-1], "ITER": itmax, "reached_ITER": itmax is not None and it[-1] >= itmax - 1,
           "nonfinite_history": nonfin_hist, "nonfinite_restart": nf,
           "nan_verdict": "CLEAN" if nonfin_hist == 0 and not any(nf.values()) else "NONFINITE",
           "convergence": _conv(str(conv_csv))}
    p = OUTD / f"xcheck_gates_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return p


def reduce_run(run: Path, kind: str, last: int = 5) -> Path:
    yb_x, yb = common_yb()
    if kind == "forge":
        steps = sorted(int(re.match(r"res_(\d+)\.h5$", p.name).group(1)) for p in run.glob("res_[0-9]*.h5"))[-last:]
        loader = lambda s: load_forge(run, s)  # noqa: E731
    else:
        steps = sorted(int(re.match(r"restart_flow_(\d+)\.csv$", p.name).group(1)) for p in run.glob("restart_flow_[0-9]*.csv"))[-last:]
        loader = lambda s: load_su2(run, s)  # noqa: E731
    rec = {}
    for k, st in enumerate(steps):
        xy, ro, ux, uy, T, tke = loader(st)
        rec[f"s{k}_nonfinite"] = np.array(int(sum(np.count_nonzero(~np.isfinite(a)) for a in (ro, ux, uy, T, tke))))
        ni, nj, S = mesh_info(len(ro))
        o = reduce_fields(xy, ro, ux, uy, T, tke, kind == "su2", ni, nj, S, yb_x, yb)   # SU2 のエネルギーは k を含む
        for key, v in o.items():
            rec[f"s{k}_{key}"] = np.asarray(v)
        print(f"[reduce] {run.name} {st}: Q_w {float(o['Q_w']) / 1e6:.4f} MW、ṁ 入口 {float(o['mdot_in']):.4f} 出口 {float(o['mdot_out']):.4f}", flush=True)
    p = OUTD / f"xcheck_{run.name}.npz"
    np.savez(p, steps=np.array(steps), kind=np.array(kind), **rec)
    return p


def summarize(name):
    """xcheck_<run>.npz → 判定窓の平均と時系列 (x = 40/70/94 の δ_loc・θ_r・τ_w、Q_w、断熱の T_w、収支)。"""
    Z = np.load(OUTD / f"xcheck_{name}.npz"); n = len(Z["steps"])
    g = lambda k, key: np.interp(XE, Z[f"s{k}_x"], Z[f"s{k}_{key}"])  # noqa: E731
    gw = lambda k, key: np.interp(XE, Z[f"s{k}_x_w"], Z[f"s{k}_{key}"])  # noqa: E731
    d = {key: np.array([g(k, key) for k in range(n)]) for key in ("delta_loc", "theta_r", "delta_loc_s125", "theta_r_s125")}
    d.update({key: np.array([gw(k, key) for k in range(n)]) for key in ("tau_w", "c_f", "q_w", "T_w")})
    d["Q_w"] = np.array([float(Z[f"s{k}_Q_w"]) for k in range(n)])
    for key in ("mdot_in", "mdot_out", "Hdot_in", "Hdot_out"):
        d[key] = np.array([float(Z[f"s{k}_{key}"]) for k in range(n)])
    d["steps"] = np.array(Z["steps"], dtype=float)
    d["nonfinite"] = np.array([int(Z[f"s{k}_nonfinite"]) for k in range(n)])
    return d


def compare() -> dict:
    sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda/tools"))
    from check_quasisteady import classify
    runs = {"forge_plain_ad": "run_0184_ns_coldmesh_cpg_ad_plain", "forge_plain_tw": "run_0185_ns_coldmesh_cpg_tw300_plain",
            "forge_prod_ad": "run_0186_ns_coldmesh_cpg_ad_dilat2", "forge_prod_tw": "run_0187_ns_coldmesh_cpg_tw300_dilat2",
            "su2_ad": "run_0188_su2_coldmesh_cpg_ad", "su2_tw": "run_0189_su2_coldmesh_cpg_tw300"}
    D = {k: summarize(v) for k, v in runs.items() if (OUTD / f"xcheck_{v}.npz").is_file()}
    out = {"runs": runs, "available": sorted(D), "gates": {}, "metrics": {}}
    ix = [int(np.argmin(np.abs(XE - x))) for x in XS]
    # 前提: 準定常 (5 枚、drift・osc 0.1 %) と収支
    for k, d in D.items():
        bad = []
        for key in ("delta_loc", "theta_r", "c_f"):
            for i in ix:
                v = classify(d["steps"], d[key][:, i], 1.0, 0.001, 0.001, 5)[0]
                if v != "STEADY":
                    bad.append(f"{key} x={XE[i]:.0f} {v}")
        if k.endswith("_ad"):
            for i in ix:
                v = classify(d["steps"], d["T_w"][:, i], 1.0, 0.001, 0.001, 5)[0]
                if v != "STEADY":
                    bad.append(f"T_w x={XE[i]:.0f} {v}")
        else:
            v = classify(d["steps"], d["Q_w"], 1.0, 0.001, 0.001, 5)[0]
            if v != "STEADY":
                bad.append(f"Q_w {v}")
        mb = float(np.mean(d["mdot_in"] - d["mdot_out"]) / np.mean(d["mdot_in"]))
        hb = float(np.mean(d["Hdot_in"] - d["Hdot_out"] - d["Q_w"]) / np.mean(d["Hdot_in"]))
        gp = OUTD / f"xcheck_gates_{runs[k]}.json"
        g = json.loads(gp.read_text()) if gp.is_file() else None
        nonfin = int(sum(int(Z) for Z in d["nonfinite"]))
        g_ok = (g is not None and g["nan_verdict"] == "CLEAN" and not g["convergence"]["diverged"] and not g["convergence"]["rising"]
                and nonfin == 0 and (g["kind"] == "forge" or g["reached_ITER"]))
        out["gates"][k] = {"not_steady": bad, "mass_balance": mb, "enthalpy_balance": hb, "nonfinite_window": nonfin,
                           "gates_file": gp.name if g else "(無い → 判定不能)", "nan": g and g["nan_verdict"],
                           "convergence": g and g["convergence"]["verdict"], "rising": g and g["convergence"]["rising"],
                           "ok": g_ok and (not bad) and abs(mb) <= 1e-3 and abs(hb) <= 1e-3}
    m = lambda k, key: D[k][key].mean(0)  # noqa: E731
    t = (XE >= 40) & (XE <= 94)
    if all(k in D for k in ("forge_plain_ad", "forge_plain_tw", "su2_ad", "su2_tw")):
        Rf = m("forge_plain_tw", "delta_loc") / m("forge_plain_ad", "delta_loc"); Rs = m("su2_tw", "delta_loc") / m("su2_ad", "delta_loc")
        met = {"R_loc_max_rel": float(np.max(np.abs(Rf / Rs - 1.0)[t])), "R_loc_forge": Rf[ix].tolist(), "R_loc_su2": Rs[ix].tolist()}
        for arm in ("ad", "tw"):
            for key in ("delta_loc", "theta_r", "c_f", "q_w"):
                if key == "q_w" and arm == "ad":
                    continue
                a = m(f"forge_plain_{arm}", key); b = m(f"su2_{arm}", key)
                good = t & np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
                met[f"{key}_{arm}_max_rel"] = float(np.max(np.abs(a[good] / b[good] - 1.0))) if good.any() else None
        met["Q_w_tw_rel"] = float(np.mean(D["forge_plain_tw"]["Q_w"]) / np.mean(D["su2_tw"]["Q_w"]) - 1.0)
        lim = {"R_loc_max_rel": 0.02, "delta_loc_ad_max_rel": 0.03, "delta_loc_tw_max_rel": 0.03, "theta_r_ad_max_rel": 0.03,
               "theta_r_tw_max_rel": 0.03, "c_f_ad_max_rel": 0.03, "c_f_tw_max_rel": 0.03, "q_w_tw_max_rel": 0.05}
        met["pass"] = {k2: (met[k2] is not None and met[k2] <= v) for k2, v in lim.items()}
        gates_ok = all(out["gates"][k]["ok"] for k in ("forge_plain_ad", "forge_plain_tw", "su2_ad", "su2_tw"))
        out["VERDICT"] = ("判定不能 (前提不成立)" if not gates_ok else
                          ("forge の解き方は冷却ノズルでも SU2 と一致 (この CPG 条件・格子・比較量で)" if all(met["pass"].values())
                           else "不一致 (forge 側の原因を調べる — 諮問)"))
        out["metrics"]["plain_vs_su2"] = met
    for tag, (a, b) in {"prod_minus_plain_R": ("forge_prod", "forge_plain")}.items():
        if all(f"{p_}_{w}" in D for p_ in (a, b) for w in ("ad", "tw")):
            Ra = m(f"{a}_tw", "delta_loc") / m(f"{a}_ad", "delta_loc"); Rb = m(f"{b}_tw", "delta_loc") / m(f"{b}_ad", "delta_loc")
            out["metrics"][tag] = {"R_prod": Ra[ix].tolist(), "R_plain": Rb[ix].tolist(), "max_rel": float(np.max(np.abs(Ra / Rb - 1.0)[t]))}
    (OUTD / "xcheck.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    for c in ("reduce-forge", "reduce-su2"):
        p = sp.add_parser(c); p.add_argument("run"); p.add_argument("--last", type=int, default=5)
    for c in ("gates-forge", "gates-su2"):
        p = sp.add_parser(c); p.add_argument("run")
    sp.add_parser("compare")
    a = ap.parse_args()
    if a.cmd == "reduce-forge":
        print(reduce_run(HERE / a.run, "forge", a.last))
    elif a.cmd == "reduce-su2":
        print(reduce_run(HERE / a.run, "su2", a.last))
    elif a.cmd == "gates-forge":
        sys.path.insert(0, str(HERE))
        print(gates_forge(HERE / a.run))
    elif a.cmd == "gates-su2":
        print(gates_su2(HERE / a.run))
    else:
        compare()
