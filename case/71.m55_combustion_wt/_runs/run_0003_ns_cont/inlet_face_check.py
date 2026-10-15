#!/usr/bin/env python3
r"""入口の面の実効の流入とエネルギー流束の確認 (forge 0 step)。plan tooling-nozzle-card-local-mesh §5.1 #3d「入口の面の実効の流入と
エネルギー流束は延長の序盤に確かめる」、codex 諮問 `notes/reviews/2026-10-15-m55-ns-unsettled-3-diagnose.md`「該当 9 点に接続する
面の流入・逆流の別、境界状態、実流入のエンタルピー流束を確認する。実効入口条件の違反なら修正が先」。

ソルバの経路 (node、`inlet_Pressure`、TP 多成分) をそのまま写す:
  - 入口の境界の半割面 (nozzle.h5 の `BCONDS/<入口>/iPlanes`・`iCells`) ごとに、節点の状態から境界の状態 (bvar) を作る
    (`boundaryCond_d.cu` の `inlet_Pressure_d`): Un_c = U·n (n は外向きの単位法線)、Un_c < 0 なら流入 — mach_c = Un_c/sonic、
    rf = 0 なので mach_new = |mach_c|、(Pt, Tt, 入口の組成 Yin, mach_new) から `thermo_isentropic_from_total_mix` (同じ Newton・
    同じ打ち切り) で Ts・Ps・ρ・|u|、速度は −|u| n。Un_c ≥ 0 は逆流の分岐 (U = 0、P = Pt、T = Tt、ρ = Pt/(R Tt))。
  - 境界の流束 (`convectiveFlux_boundary_d.inc.cuh`): mdot = W·(ρ_b U_b·n) (W は面のベクトルの大きさ; 軸対称は /PLANES/rSurfVect)、
    エネルギーは流入 (mdot < 0) なら mdot·Ht_R、Ht_R = (roe_b + P_b)/ρ_b = e(Ts) + |u|²/2 + Ps/ρ_b = h(Ts) + |u|²/2。
    入口の半割面には粘性・熱伝導の流束は無い (`viscousFlux_d.cu` は壁の半割面だけ)。
  - 節点の状態は 2 通り: A = 出力の原始量 (Ux・Uy・P・T・sonic・Y)、B = 保存量から thermoFloat 0 の経路 (double の組成の正規化・
    e_mix(T) = e の Newton・P = ρ R T・sonic = √(γ P/ρ)) で作り直したもの (延長の最初の step が見る状態)。

報告: 流入・逆流の面の数と位置、境界の状態と節点の状態の差、面ごとの Ht_R − h0_in (h0_in = h_mix(Tt, Yin))、実効の流入の
質量流束・エネルギー流束の合計 (2π 倍) と ṁ·h0_in の差、節点の値で取った ṁ_in (停止規則の r の分子と同じ式) との比、
`total_state` で T0 > Tt + 1 K の入口の点 (run_0002_ns の res_183000 で 9 点) の上の全部。
違反の判定 (呼び出し側が事前に置いた目安; 0 step の整合検査で、収束の判定ではない): (v1) 逆流の面がある、(v2) 流入の面で
|Ht_R − h0_in| > 1e−6 h0_in、(v3) 合計のエネルギー流束と ṁ·h0_in の差 > 1e−6、(v4) 流入の面で mach_new ≥ 1 (亜音速の入口の前提)。

使い方 (AWS A): python3 inlet_face_check.py --tools ~/forge-cardmesh/solver_density_cuda/tools --run <run_0002_ns> --step 183000 \
    --out-json <出力先>/inlet_face_check.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

TMIN, TMAX = 50.0, 6000.0
TOL = 1.0e-6


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _sc(f, Y, T) -> float:
    """_TPGas の h・cp・s0 (配列の T を取る) をスカラの T で呼ぶ。"""
    return float(f(Y, np.atleast_1d(np.asarray(T, dtype=float)))[0])


def isentropic_from_total(gas, Yin, Pt, Tt, mach):
    """`thermo_isentropic_from_total_mix` (thermo_d.cuh:893) の写し (スカラ)。"""
    R = float(gas.Rmix(Yin))
    h0 = _sc(gas.h, Yin, Tt)
    Ts = max(Tt / (1.0 + 0.2 * mach * mach), 10.0)
    for _ in range(25):
        cp = _sc(gas.cp, Yin, Ts)
        g = cp / max(cp - R, 1.0e-6)
        a2 = g * R * Ts
        f = _sc(gas.h, Yin, Ts) + 0.5 * mach * mach * a2 - h0
        df = cp + 0.5 * mach * mach * g * R
        dT = min(max(f / df, -0.4 * Ts), 0.4 * Ts)
        Ts = min(max(Ts - dT, 10.0), Tt)
        if abs(dT) < 1.0e-3 + 1.0e-6 * Ts:
            break
    cp = _sc(gas.cp, Yin, Ts)
    g = cp / max(cp - R, 1.0e-6)
    a = np.sqrt(g * R * Ts)
    Ps = Pt * np.exp(-(_sc(gas.s0, Yin, Tt) - _sc(gas.s0, Yin, Ts)) / R)
    return Ts, Ps, Ps / (R * Ts), mach * a


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.tools).resolve()))
    import forge_species as fsp
    from total_quantities import total_state

    run = Path(a.run).resolve()
    res = run / f"res_{a.step}.h5"
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    pp = cfg["physProp"]
    roMin, pMin = float(pp.get("roMin", 1.0e-4)), float(pp.get("pMin", 1.0))
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    name, inl = next((k, v) for k, v in bc.items() if v["kind"] == "inlet_Pressure")
    pid = int(inl["physID"])
    fl = inl["floats"]
    Pt, Tt = float(fl["Pt"]), float(fl["Tt"])
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])

    th = fsp.run_thermo(str(run), res_path=str(res), source="auto")
    gas = fsp.thermo_gas(th)
    nsp = len(th["names"])
    yin = np.array([max(float(fl[f"Y{s}"]), 0.0) for s in range(nsp)])
    yin = yin / yin.sum()
    Yin = list(yin)
    h0_in = _sc(gas.h, Yin, Tt)
    R_in = float(gas.Rmix(Yin))

    with h5py.File(run / "nozzle.h5", "r") as h:
        coord = np.array(h["/MESH/COORD"][:]).reshape(-1, 3)
        g = h[f"BCONDS/{pid}"]
        kind = g.attrs["bcondKind"]
        kind = kind.decode() if isinstance(kind, bytes) else str(kind)
        ipl = np.array(g["iPlanes"][:]); icl = np.array(g["iCells"][:])
        sv = np.array(h["PLANES/surfVect"][:]).reshape(-1, 3)[ipl]
        rsv = np.array(h["PLANES/rSurfVect"][:]).reshape(-1, 3)[ipl]
        pc = np.array(h["PLANES/centCoords"][:]).reshape(-1, 3)[ipl]
    if kind != "inlet_Pressure":
        raise SystemExit(f"BCONDS/{pid} の種類 {kind} が inlet_Pressure でない")
    with h5py.File(res, "r") as h:
        V = {k: np.array(h["VALUE"][k][:], dtype=np.float64) for k in
             ("ro", "roUx", "roUy", "roUz", "roe", "Ux", "Uy", "Uz", "P", "T", "sonic", "h0")}
        roY = [np.array(h["VALUE"][f"roY{s}"][:], dtype=np.float64) for s in range(nsp)]
        Yout = [np.array(h["VALUE"][f"Y{s}"][:], dtype=np.float64) for s in range(nsp)]

    # 面の外向きの単位法線 (surfVect の向き = ic0 → ghost の外向き; 入口は x の最小の列なので −x のはず)
    n = sv / np.linalg.norm(sv, axis=1, keepdims=True)
    W = np.linalg.norm(rsv, axis=1)                       # 軸対称の面の重み (1 rad あたり; ソルバの W_f = Σ r_k S_k)
    col_x = coord[icl, 0]
    i_idx, j_idx = icl // nj, icl % nj

    def node_state(src):
        if src == "A":
            ro = V["ro"][icl]; U = np.stack([V["Ux"][icl], V["Uy"][icl], V["Uz"][icl]], 1)
            return ro, U, V["P"][icl], V["T"][icl], V["sonic"][icl]
        ro = np.maximum(V["ro"][icl], roMin)
        U = np.stack([V["roUx"][icl], V["roUy"][icl], V["roUz"][icl]], 1) / ro[:, None]
        e = V["roe"][icl] / ro - 0.5 * np.sum(U * U, 1)
        ys = [np.maximum(r[icl] / ro, 0.0) for r in roY]
        s = np.maximum(sum(ys), 1.0e-30)
        Y = [y / s for y in ys]
        R = gas.Rmix(Y)
        T = np.clip(V["T"][icl].copy(), TMIN, TMAX)
        for _ in range(60):
            cv = np.maximum(gas.cp(Y, T) - R, 1.0e-2 * R)
            dT = np.clip((gas.h(Y, T) - R * T - e) / cv, -0.5 * T, 0.5 * T)
            T = np.clip(T - dT, TMIN, TMAX)
            if np.max(np.abs(dT)) < 1.0e-10 * TMAX:
                break
        P = np.maximum(ro * R * T, pMin)
        cp = gas.cp(Y, T)
        gam = cp / (cp - R)
        return ro, U, P, T, np.sqrt(gam * P / ro)

    sA = total_state(str(run), str(res))
    T0A = np.asarray(sA["T0"], float)[icl]
    hot9 = T0A > Tt + 1.0
    out = {"run": str(run), "res": res.name, "res_sha256": sha256(res), "inlet": {"name": name, "physID": pid, "Pt": Pt, "Tt": Tt,
           "Yin": yin.tolist(), "h0_in": h0_in, "R_in": R_in},
           "faces": int(len(ipl)), "nodes_column_i0": bool(np.all(i_idx == 0)), "j_range": [int(j_idx.min()), int(j_idx.max())],
           "x_column_spread_m": float(col_x.max() - col_x.min()),
           "normal": {"nx_min": float(n[:, 0].min()), "nx_max": float(n[:, 0].max()), "max_abs_ny": float(np.abs(n[:, 1]).max())},
           "states": {}, "rule": ("v1 逆流の面がある / v2 流入の面で |Ht_R − h0_in| > 1e−6 h0_in / v3 合計のエネルギー流束と ṁ·h0_in の差 > 1e−6 / "
                                  "v4 流入の面で mach_new ≥ 1 (0 step の整合検査、収束の判定ではない)")}
    for src in ("A", "B"):
        ro, U, P, T, a_s = node_state(src)
        Un = np.sum(U * n, 1)
        inflow = Un < 0.0
        rows = []
        mdot = np.zeros(len(ipl)); eflux = np.zeros(len(ipl)); dH = np.full(len(ipl), np.nan); mach_new = np.zeros(len(ipl))
        Ts_b = np.zeros(len(ipl)); Ps_b = np.zeros(len(ipl)); ro_b = np.zeros(len(ipl)); um_b = np.zeros(len(ipl))
        for k in range(len(ipl)):
            if inflow[k]:
                mach_new[k] = abs(Un[k] / a_s[k])
                Ts_b[k], Ps_b[k], ro_b[k], um_b[k] = isentropic_from_total(gas, Yin, Pt, Tt, mach_new[k])
                Ht = _sc(gas.h, Yin, Ts_b[k]) + 0.5 * um_b[k] ** 2
                dH[k] = Ht - h0_in
                mdot[k] = -W[k] * ro_b[k] * um_b[k]              # 外向きを正 (ソルバの mdot と同じ符号)
                eflux[k] = mdot[k] * Ht
            else:
                Ts_b[k], Ps_b[k], ro_b[k], um_b[k] = Tt, Pt, Pt / (R_in * Tt), 0.0
        m_in = -2.0 * np.pi * mdot.sum()                           # 流入を正
        e_in = -2.0 * np.pi * eflux.sum()
        # 節点の値の ṁ_in (停止規則の mass_flows と同じ: 列の折れ線を通る 2π∫ρ u_x r dr; 入口は x 一定なので dx 項 0)
        r_col = coord[icl, 1]
        o = np.argsort(r_col)
        q = (ro * U[:, 0])[o]
        m_node = 2.0 * np.pi * float(np.trapezoid(q * r_col[o], r_col[o]))
        st = {
            "n_inflow": int(inflow.sum()), "n_reverse": int((Un > 0.0).sum()), "n_zero": int((Un == 0.0).sum()),
            "reverse_j_from_wall": sorted(int(nj - 1 - j) for j in j_idx[Un > 0.0]),
            "zero_j_from_wall": sorted(int(nj - 1 - j) for j in j_idx[Un == 0.0]),
            "mach_new_max": float(mach_new.max()),
            "dH_over_h0in": {"max_abs": float(np.nanmax(np.abs(dH))) / h0_in if inflow.any() else None},
            "mdot_in_face_kg_s": m_in, "energy_in_face_W": e_in, "energy_minus_mdot_h0in_rel": (e_in - m_in * h0_in) / (m_in * h0_in),
            "mdot_in_node_kg_s": m_node, "face_over_node_mdot": m_in / m_node,
            "node_minus_boundary": {"P_rel_max_abs": float(np.max(np.abs(P[inflow] - Ps_b[inflow]) / Ps_b[inflow])) if inflow.any() else None,
                                    "T_max_abs_K": float(np.max(np.abs(T[inflow] - Ts_b[inflow]))) if inflow.any() else None,
                                    "T_node_max_minus_Tt": float((T - Tt).max())},
        }
        for k in np.where(hot9)[0]:
            rows.append({"node": int(icl[k]), "j": int(j_idx[k]), "j_from_wall": int(nj - 1 - j_idx[k]), "r_m": float(coord[icl[k], 1]),
                         "T0_A": float(T0A[k]), "T_node": float(T[k]), "Un_c": float(Un[k]), "inflow": bool(inflow[k]),
                         "mach_new": float(mach_new[k]), "Ts_b": float(Ts_b[k]), "Ps_b": float(Ps_b[k]), "P_node": float(P[k]),
                         "Ht_R_minus_h0in": (float(dH[k]) if inflow[k] else None),
                         "mdot_face_per_rad": float(mdot[k]), "share_of_inflow": float(mdot[k] / mdot.sum()) if mdot.sum() != 0 else None})
        st["hot_T0_points"] = rows
        # 面ごとの表 (壁からの番号の順): 境界の状態と節点の状態、面の質量流束 (1 rad あたり、外向きを正)
        o2 = np.argsort(-j_idx)
        st["faces_table"] = {"j_from_wall": [int(nj - 1 - j) for j in j_idx[o2]], "Un_c": Un[o2].tolist(), "T_node": T[o2].tolist(),
                             "Ts_b": Ts_b[o2].tolist(), "P_node": P[o2].tolist(), "Ps_b": Ps_b[o2].tolist(),
                             "mach_new": mach_new[o2].tolist(), "mdot_per_rad": mdot[o2].tolist(), "W_per_rad": W[o2].tolist()}
        # v1 は真の逆流 (Un_c > 0) だけ。壁との角の節点 (壁の Dirichlet で U = 0 → Un_c = 0) はソルバでも逆流の分岐 (mdot = 0) に入るが、
        # 壁の節点として想定どおりなので zero に分けて記録する
        viol = {"v1_reverse": int((Un > 0.0).sum()) > 0,
                "v2_dH": bool(inflow.any() and np.nanmax(np.abs(dH)) > TOL * h0_in),
                "v3_total_energy": bool(abs((e_in - m_in * h0_in) / (m_in * h0_in)) > TOL),
                "v4_supersonic": bool(mach_new.max() >= 1.0)}
        st["violations"] = viol
        out["states"][src] = st
        print(f"[inlet_face_check] {src}: 流入 {st['n_inflow']} / 逆流 {st['n_reverse']} (壁から {st['reverse_j_from_wall']}) / "
              f"零 {st['n_zero']} (壁から {st['zero_j_from_wall']}), "
              f"|Ht_R−h0_in|/h0_in 最大 {st['dH_over_h0in']['max_abs']}, ṁ_in 面 {m_in:.6f} 節点 {m_node:.6f} kg/s, "
              f"E−ṁh0 相対 {st['energy_minus_mdot_h0in_rel']:.3e}, mach 最大 {st['mach_new_max']:.4f}, 違反 {viol}", flush=True)
    out["hot_T0_n"] = int(hot9.sum())
    Path(a.out_json).write_text(json.dumps(out, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
