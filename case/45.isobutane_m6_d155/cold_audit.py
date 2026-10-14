"""冷却壁の腕の質量の収支の監査 (plan tooling-nozzle-isothermal-wall-chain §5.1 #27 (iv)、codex diagnose 2026-10-08)。
監査 run (cold_cfl.py prep で親の最終場からビット一致で restart、2 step、FORGE_DUMP_MASSFLUX と extraFields [res_ro, volume]) の
最初の評価の面の質量流束 (massflux.bin、flow_float = FP64) と res_ro (res_1.h5) から、節点集合ごとの離散の収支を組む。

  (1) 節点ごとの恒等式: res_ro(n) と、n に接する面の流束の符号付きの和 (内部面 + 境界面) が一致するか (符号・状態の対応の確認)。
  (2) 領域 Ω_I = {列 i ≤ I の節点} について: 境界 (入口・壁・軸) と切断面 (Ω_I と外を分ける内部面) を通る流出の和と −Σ_{Ω_I} res_ro。
  (3) 切断面の数値流束と、同じ列の節点値の台形 (cold_xcheck と同じ 2π∫(F_x dr − F_r dx) r) の比。
節点の番号 = 双対 CV の番号 (n = i·nj + j、i: 入口 → 出口、j: 軸 → 壁) を CELLS/centCoords と MESH/COORD で確かめる。
usage (AWS): python3 cold_audit.py <audit_run> [--cuts 0,200,1000,2000,3000,4000,4717]  → _band_ab/cold_pair/audit_<run>.json
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
OUTD = HERE / "_band_ab" / "cold_pair"


def planes(h):
    st = np.asarray(h["PLANES/STRUCT"][:], dtype=np.int64)
    n = int(h["PLANES/surfArea"].shape[0])
    c0 = np.empty(n, np.int64); c1 = np.full(n, -1, np.int64)
    k = 0
    for ip in range(n):                       # [節点数, 節点…, セル数, セル…]
        k += 1 + st[k]
        nc = st[k]; c0[ip] = st[k + 1]
        if nc >= 2:
            c1[ip] = st[k + 2]
        k += 1 + nc
    if k != len(st):
        raise SystemExit(f"PLANES/STRUCT の読み取りが末尾と合わない ({k} != {len(st)})")
    return c0, c1


def main(run: Path, cuts):
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
        cc = np.asarray(h["CELLS/centCoords"][:], float).reshape(-1, 3)[:, :2]
        c0, c1 = planes(h)
        bpl = {pid: np.asarray(h[f"BCONDS/{pid}/iPlanes"][:], np.int64) for pid in h["BCONDS"].keys()}
    nn = len(xy); ni = 4719; nj = nn // ni
    if ni * nj != nn:
        raise SystemExit("節点数が ni の倍数でない")
    # 節点 = 双対 CV の番号の確認: 双対の重心は節点の近く (最近接の格子間隔より近い)
    dd = np.linalg.norm(cc[:nn] - xy, axis=1)
    mf = np.fromfile(run / "massflux.bin", dtype=np.float64)
    if len(mf) != len(c0):
        raise SystemExit(f"massflux.bin の面数 {len(mf)} が PLANES {len(c0)} と違う (型が float64 でない?)")
    with h5py.File(run / "res_1.h5", "r") as h:
        res = np.asarray(h["VALUE/res_ro"][:], float)
    with h5py.File(run / "res_0.h5", "r") as h:
        ro = np.asarray(h["VALUE/ro"][:], float); ux = np.asarray(h["VALUE/Ux"][:], float); uy = np.asarray(h["VALUE/Uy"][:], float)
    st = np.fromfile(str(run / "massflux.bin") + ".state", dtype=np.float64).reshape(6, -1)[:, :nn]
    out = {"run": run.name, "nodes": nn, "ni": ni, "nj": nj, "planes": int(len(c0)),
           "dual_centroid_offset_max": float(dd.max()), "state_ro_vs_res0_maxrel": float(np.max(np.abs(st[0] / ro - 1.0)))}
    internal = c1 >= 0
    internal &= c1 < nn
    bnd = ~internal
    # (1) 節点ごとの流束の和 (c0 から c1 へ正と仮定; 境界面は c0 の外向き)
    div = np.zeros(nn)
    np.add.at(div, c0, mf)
    np.add.at(div, c1[internal], -mf[internal])
    # res_ro の符号を最小二乗で当てる (res = −div を期待)
    a = float(np.dot(res, -div) / np.dot(div, div))
    r1 = res + div
    out["identity"] = {"res_over_minus_div_ls": a, "max_abs_res_plus_div": float(np.max(np.abs(r1))),
                       "max_abs_res": float(np.max(np.abs(res))), "rel": float(np.max(np.abs(r1)) / max(np.max(np.abs(res)), 1e-300)),
                       "worst_nodes": [int(v) for v in np.argsort(-np.abs(r1))[:5]]}
    # 境界ごとの流出の和
    out["boundary_outflow"] = {pid: float(mf[ip].sum()) for pid, ip in bpl.items()}
    out["sum_res_all"] = float(res.sum())
    # (2)(3) 領域 Ω_I
    col = np.arange(nn) // nj
    x = xy[:, 0].reshape(ni, nj); r = xy[:, 1].reshape(ni, nj)
    Fx = (ro * ux).reshape(ni, nj) * r; Fr = (ro * uy).reshape(ni, nj) * r
    rows = []
    for I in cuts:
        inside = col <= I
        cut = internal & (inside[c0] != inside[np.where(internal, c1, 0)])
        sgn = np.where(inside[c0], 1.0, -1.0)
        cut_out = float(np.sum(sgn[cut] * mf[cut]))
        bout = {pid: float(mf[ip][inside[c0[ip]]].sum()) for pid, ip in bpl.items()}
        tot_out = cut_out + sum(bout.values())
        sres = float(res[inside].sum())
        trap = 2 * np.pi * float(np.sum(0.5 * (Fx[I, 1:] + Fx[I, :-1]) * np.diff(r[I])) - np.sum(0.5 * (Fr[I, 1:] + Fr[I, :-1]) * np.diff(x[I])))
        rows.append({"I": int(I), "x_w": float(x[I, -1]), "cut_outflow": cut_out, "boundary_outflow": bout, "total_outflow": tot_out,
                     "sum_res": sres, "closure": tot_out + sres, "trapezoid_mdot": trap})
    out["regions"] = rows
    OUTD.mkdir(parents=True, exist_ok=True)
    p = OUTD / f"audit_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return p


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--cuts", default="0,200,1000,2000,3000,4000,4717")
    a = ap.parse_args()
    main(HERE / a.run, [int(v) for v in a.cuts.split(",")])
