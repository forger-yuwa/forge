#!/usr/bin/env python3
r"""③ベル (M) L100: δ* 補正壁の判定 — plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-4) の事前登録 (run を回す前に commit で凍結)。

組 (NS, 狙いの Euler):
  U (補正なし)  = (run_0019_tp_m_ns_l100, run_0018_tp_m_euler_l100_pin)            — 同じ非粘性壁
  Vext・Vzero  = (補正壁の NS [problem_bell_tp_m_ns_l100_corr_{vext,vzero}.yaml], 狙い [problem_bell_tp_m_euler_l100_inv.yaml])
  格子の項      = 狙いを ×1.5 格子 (problem_bell_tp_m_euler_l100_inv_fine.yaml) に替えた値
場: 各 run の最終の res (NS 本段 60000、Euler 本段 30000)。run のゲート (GATES.json の全 PASS) を前提に読む (落ちていたら判定不能)。

指標 (組ごと):
  pw    : 壁圧の誤差 = (p_NS − p_E)/p_E の Δx 重みの RMS、x ∈ [1.17, L]。壁の節点 (j = nj−1) の p を共通の x 格子 (4001 点) へ補間
  pw_t  : 同じ量の x ∈ [0, 1.17] (推定区間、記録のみ)
  axM   : 軸上 Mach の差 M_NS − M_E の Δx 重みの RMS、x ∈ [0, L]
  exM   : 出口面 (最終列) の Mach の差の面積重み (r dr) の RMS、r ∈ [0, R_C] (R_C = 6.5 r_t — 補正なしの NS の出口の境界層の縁
          r ≈ 7.07 − 0.40 = 6.67 の内側、両方の組で同じ)
  edM   : 境界層の外縁の Mach — 点 (x, r_s(x))、r_s = r_wall,NS(x) − y_b(x) (y_b は run_0019 の edge 抽出の帯の位置、凍結:
          dstar_l100_ns_run0019_constx_edge.npz の band_y_b) で NS と Euler の Mach を比べた差の Δx 重みの RMS、x ∈ [1.17, L]。
          一定 x の断面に直してから (各行を x 方向に補間) 半径方向に補間
  eta   : η_exit (metrics.json)
不確かさ (許容幅) u_k = max_V |m_k(V; 狙い fine) − m_k(V; 狙い)| (狙いの格子の差だけ)。η は u_η = |η(狙い fine) − η(狙い)|。
推定区間への感度 s_k = |m_k(Vext) − m_k(Vzero)|、s_η = |η_Vext − η_Vzero| は別に報告し、許容幅には使わない
  (codex 2026-10-12 `notes/reviews/2026-10-12-bell-m-corrected-wall-prereg-diagnose.md`: 幅に入れると片方の悪化を隠す)。
判定 (各補正壁 V について):
  主     : pw(V) ≤ 0.7·pw(U)  (壁圧の誤差が補正なしから 30 % 以上減る)
  併記   : axM・exM・edM について m(V) ≤ m(U) + u (不確かさ以上に悪化しない)
  η      : 改善は必須にしない — Δη = η_V − η_U と u_η を報告
結論: 両方の V で主・併記とも成立 → 「補正は効く (推定区間 x < 1.17 の採り方に依らない)」、片方だけ → 「推定区間の採り方で結果が
  変わる」、どちらも不成立 → 「補正の効きを支持しない」。出力 corrected_wall_judgement.json。

使い方: python judge_corrected_wall.py [--runs vext=<run_dir> vzero=<run_dir> target=<run_dir> target_fine=<run_dir>]
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = HERE.parents[2]
X_CAL = 1.17
R_C = 6.5
NX = 4001
PW_REDUCTION = 0.30
RUN_U = ("run_0019_tp_m_ns_l100", "run_0018_tp_m_euler_l100_pin")
EDGE = HERE / "dstar_l100_ns_run0019_constx_edge.npz"


def load(run: Path) -> dict:
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj, S = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]), float(info["scale_m"])
    res = sorted(run.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
    with h5py.File(run / "nozzle.h5") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res) as f:
        V = f["VALUE"]
        P = V["P"][:]
        M = np.hypot(V["Ux"][:], V["Uy"][:]) / V["sonic"][:]
    if len(P) != ni * nj:
        raise SystemExit(f"{run}: node 構造格子でない")
    g = json.loads((run / "GATES.json").read_text()) if (run / "GATES.json").exists() else None
    met = json.loads((run / "metrics.json").read_text())
    return dict(run=run, res=res.name, x=c[:, 0].reshape(ni, nj) / S, r=c[:, 1].reshape(ni, nj) / S,
                P=P.reshape(ni, nj), M=M.reshape(ni, nj), gates=g, eta=float(met["eta_exit"]))


def gates_ok(d: dict) -> bool:
    g = d["gates"]
    if not g:
        return False
    return g.get("status") == "PASS"


def const_x(d: dict):
    """各行 j を x 方向に補間して列を軸の x の一定 x の断面にそろえた (x_t, r, M)。"""
    xt = d["x"][:, 0]
    nj = d["x"].shape[1]
    r = np.empty_like(d["r"]); M = np.empty_like(d["M"])
    for j in range(nj):
        r[:, j] = np.interp(xt, d["x"][:, j], d["r"][:, j]); M[:, j] = np.interp(xt, d["x"][:, j], d["M"][:, j])
    return xt, r, M


def wrms(v, w):
    return float(np.sqrt(np.sum(w * v * v) / np.sum(w)))


def trap_w(x):
    w = np.zeros_like(x); w[1:-1] = 0.5 * (x[2:] - x[:-2]); w[0] = 0.5 * (x[1] - x[0]); w[-1] = 0.5 * (x[-1] - x[-2])
    return w


def metrics(ns: dict, eu: dict, yb_x, yb) -> dict:
    L = min(ns["x"][-1, -1], eu["x"][-1, -1])
    xg = np.linspace(0.0, L, NX); w = trap_w(xg)
    pN = np.interp(xg, ns["x"][:, -1], ns["P"][:, -1]); pE = np.interp(xg, eu["x"][:, -1], eu["P"][:, -1])
    e = pN / pE - 1.0
    cal, thr = xg >= X_CAL, xg < X_CAL
    out = {"pw": wrms(e[cal], w[cal]), "pw_t": wrms(e[thr], w[thr]), "pw_max_cal": float(np.max(np.abs(e[cal])))}
    mN = np.interp(xg, ns["x"][:, 0], ns["M"][:, 0]); mE = np.interp(xg, eu["x"][:, 0], eu["M"][:, 0])
    out["axM"] = wrms(mN - mE, w)
    rg = np.linspace(0.0, R_C, 2001)
    MN = np.interp(rg, ns["r"][-1, :], ns["M"][-1, :]); ME = np.interp(rg, eu["r"][-1, :], eu["M"][-1, :])
    out["exM"] = wrms(MN - ME, trap_w(rg) * np.maximum(rg, 1e-12))
    xtN, rN, MNc = const_x(ns); xtE, rE, MEc = const_x(eu)
    xs = xg[cal]
    rwall = np.interp(xs, xtN, rN[:, -1]); rs = rwall - np.interp(xs, yb_x, yb)
    dM, ok = np.full(len(xs), np.nan), np.zeros(len(xs), bool)
    for k, (xx, rr) in enumerate(zip(xs, rs)):
        iN = int(np.clip(np.searchsorted(xtN, xx), 1, len(xtN) - 1)); iE = int(np.clip(np.searchsorted(xtE, xx), 1, len(xtE) - 1))
        aN = (xx - xtN[iN - 1]) / (xtN[iN] - xtN[iN - 1]); aE = (xx - xtE[iE - 1]) / (xtE[iE] - xtE[iE - 1])
        rcN = (1 - aN) * rN[iN - 1] + aN * rN[iN]; McN = (1 - aN) * MNc[iN - 1] + aN * MNc[iN]
        rcE = (1 - aE) * rE[iE - 1] + aE * rE[iE]; McE = (1 - aE) * MEc[iE - 1] + aE * MEc[iE]
        if rr <= rcE[-1]:
            dM[k] = np.interp(rr, rcN, McN) - np.interp(rr, rcE, McE); ok[k] = True
    out["edM"] = wrms(dM[ok], w[cal][ok]); out["edM_n_outside_euler"] = int((~ok).sum())
    out["eta_ns"], out["eta_euler"] = ns["eta"], eu["eta"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=[])
    a = ap.parse_args()
    runs = dict(kv.split("=", 1) for kv in a.runs)
    need = ("vext", "vzero", "target", "target_fine")
    if any(k not in runs for k in need):
        raise SystemExit(f"--runs に {need} を与える (例 vext=run_0024_...)")
    ed = np.load(EDGE, allow_pickle=True)
    yb_x, yb = ed["x"], ed["band_y_b"]
    D = {k: load(CASE / v) for k, v in runs.items()}
    D["U_ns"], D["U_eu"] = load(CASE / RUN_U[0]), load(CASE / RUN_U[1])
    gate = {k: {"run": str(d["run"].relative_to(ROOT)), "res": d["res"], "gates_ok": gates_ok(d)} for k, d in D.items()}
    res = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #7 (7b-4)", "frozen_script": Path(__file__).name,
           "settings": {"X_CAL": X_CAL, "R_C": R_C, "NX": NX, "PW_REDUCTION": PW_REDUCTION, "edge": EDGE.name}, "gate": gate}
    if not all(g["gates_ok"] for g in gate.values()):
        res["verdict"] = "判定不能 (run のゲートが PASS でない)"
        (HERE / "corrected_wall_judgement.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
        print(json.dumps(res, indent=1, ensure_ascii=False)); return
    m = {"U": metrics(D["U_ns"], D["U_eu"], yb_x, yb),
         "Vext": metrics(D["vext"], D["target"], yb_x, yb), "Vzero": metrics(D["vzero"], D["target"], yb_x, yb),
         "Vext_fine": metrics(D["vext"], D["target_fine"], yb_x, yb), "Vzero_fine": metrics(D["vzero"], D["target_fine"], yb_x, yb)}
    keys = ("pw", "axM", "exM", "edM")
    # 許容幅 u は狙いの格子の差だけ (codex 2026-10-12: V-ext/V-zero の差を幅に入れると、片方の悪化をその差で隠して
    # 「採り方に依らない」と判定できてしまう)。V-ext/V-zero の差は推定区間への感度 s として別に報告する
    u = {k: max(abs(m["Vext_fine"][k] - m["Vext"][k]), abs(m["Vzero_fine"][k] - m["Vzero"][k])) for k in keys + ("pw_t",)}
    s_var = {k: abs(m["Vext"][k] - m["Vzero"][k]) for k in keys + ("pw_t",)}
    u_eta = abs(D["target_fine"]["eta"] - D["target"]["eta"])
    s_eta = abs(D["vext"]["eta"] - D["vzero"]["eta"])
    judge = {}
    for V in ("Vext", "Vzero"):
        prim = m[V]["pw"] <= (1.0 - PW_REDUCTION) * m["U"]["pw"]
        sec = {k: bool(m[V][k] <= m["U"][k] + u[k]) for k in ("axM", "exM", "edM")}
        judge[V] = {"pw_reduction": 1.0 - m[V]["pw"] / m["U"]["pw"], "primary": bool(prim), "secondary": sec,
                    "pass": bool(prim and all(sec.values())), "d_eta": D[V.lower()]["eta"] - D["U_ns"]["eta"]}
    npass = sum(judge[V]["pass"] for V in judge)
    res.update({"metrics": m, "uncertainty": u, "u_eta": u_eta, "variant_sensitivity": s_var, "s_eta": s_eta, "judge": judge,
                "verdict": {2: "補正は効く (推定区間 x < 1.17 の採り方に依らない)", 1: "推定区間の採り方で結果が変わる",
                            0: "補正の効きを支持しない"}[npass]})
    (HERE / "corrected_wall_judgement.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps({k: res[k] for k in ("uncertainty", "u_eta", "variant_sensitivity", "s_eta", "judge", "verdict")}, indent=1, ensure_ascii=False))
    for k, v in m.items():
        print(k, {kk: round(vv, 6) if isinstance(vv, float) else vv for kk, vv in v.items()})


if __name__ == "__main__":
    main()
