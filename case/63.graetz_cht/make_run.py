#!/usr/bin/env python3
r"""case/63 Graetz の run ディレクトリを組む (実行はしない)。

plan [`boundary-cht-axisymmetric-graetz.md`](../../plans/accepted/boundary-cht-axisymmetric-graetz.md) §6。

    kind = dry : 乾式 1 step (加熱壁は非連成の等温 T_in。共役壁の節点確定と起動確認用)
    kind = iso : 非連成の等温壁 T_w = T_in + ΔT (§6 の G1)
    kind = cht : 共役 (固体殻、外面 Robin T_c = T_in + ΔT) (§6 の G2)

- IC: 十分発達した Poiseuille 流 (u = 2U_m(1−r²/R²)、線形の圧力勾配 dp/dx = −8μU_m/R²、T = T_in) を
  MESH/COORD (node 値の位置) から組んで VALUE に書く。上流・下流の断熱区間も同じ。
- 入口: `inlet_profile_1.csv` (列 y, ro, Ux, Uy, Uz, Ps。y で 1D 線形補間) に放物速度を書き、`ints: {inletProfile: 1}`。
  ρ・Ps は一様 (亜音速ではエントロピー Ps/ρ^γ として使われる。`boundaryCond_d.cu` inlet_uniformVelocity_d)。

    python3 case/63.graetz_cht/make_run.py dry case/63.graetz_cht/run_0001_dry_r32 --mesh case/63.graetz_cht/mesh/graetz_r32.h5
    python3 case/63.graetz_cht/make_run.py cht case/63.graetz_cht/run_00NN_... --mesh ... --solid ... --dT 10
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import h5py
import numpy as np

import graetz_common as gc
from graetz_common import axcht


def poiseuille_ic(h5path: Path, x_out: float):
    """Poiseuille + 線形圧力 + T_in を VALUE に書く (dtype は変換器の出力に揃える)。"""
    with h5py.File(h5path, "a") as f:
        xyz = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
        g = f["VALUE"]
        n = g["ro"].shape[0]
        if len(xyz) != n:
            raise SystemExit(f"MESH/COORD の点数 {len(xyz)} と VALUE の長さ {n} が違う (node 変換でない?)")
        x, r = xyz[:, 0], xyz[:, 1]
        u = 2.0 * gc.U_M * np.clip(1.0 - (r / gc.R) ** 2, 0.0, None)
        p = gc.P_OUT + 8.0 * gc.MU * gc.U_M / gc.R ** 2 * (x_out - x)
        ro = p / (gc.RGAS * gc.T_IN)
        roe = p / (gc.GAMMA - 1.0) + 0.5 * ro * u ** 2
        for name, val in (("ro", ro), ("roUx", ro * u), ("roUy", 0 * u), ("roUz", 0 * u), ("roe", roe)):
            dt = g[name].dtype
            del g[name]
            g.create_dataset(name, data=np.asarray(val, dtype=dt))
        f["VALUE"].attrs["graetz_ic"] = (f"Poiseuille U_m={gc.U_M} R={gc.R} dp/dx=-8 mu U_m/R^2 "
                                         f"p(x_out={x_out})={gc.P_OUT} T={gc.T_IN} (make_run.poiseuille_ic)")
    return ro.min(), ro.max(), p.min(), p.max()


def inlet_csv(run: Path, n: int = 201):
    y = np.linspace(0.0, gc.R, n)
    u = 2.0 * gc.U_M * (1.0 - (y / gc.R) ** 2)
    # 入口の Ps と ρ は T_in・p(入口) の一様値 (エントロピーの指定)。p(入口) は Poiseuille の線形圧力
    p_in = gc.P_OUT + 8.0 * gc.MU * gc.U_M / gc.R ** 2 * (gc.L_UP + gc.L_HEAT + gc.L_DOWN)
    ro = p_in / (gc.RGAS * gc.T_IN)
    rows = np.c_[y, np.full(n, ro), u, np.zeros(n), np.zeros(n), np.full(n, p_in)]
    np.savetxt(run / "inlet_profile_1.csv", rows, delimiter=" ", header="y ro Ux Uy Uz Ps", comments="",
               fmt="%.12e")
    return p_in, ro


BC = """# case/63 Graetz ({kind}, dT={dT} K)。make_run.py が生成。physID は gen_mesh.py の Physical Curve。
inlet:     {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {inl_ints},
            floats: {{ro: {ro:.10g}, Ux: {um:.10g}, Uy: 0.0, Uz: 0.0, Ps: {pin:.10g}}}}}
outlet:    {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {pout}, Pt: {ptb:.10g}, Tt: {tin}}}}}
wall_up:   {{physID: 3, kind: wall, outputHDFflg: 1, ints: , floats: }}
wall_heat: {{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: {heat_ints}, floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw}}}}}
wall_down: {{physID: 5, kind: wall, outputHDFflg: 1, ints: , floats: }}
axis:      {{physID: 6, kind: axis, outputHDFflg: 0, ints: , floats: }}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["dry", "iso", "cht"])
    ap.add_argument("run")
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--solid", default=None)
    ap.add_argument("--dT", type=float, default=0.0, help="T_c − T_in [K] (iso は壁温、cht は Robin の外温)")
    ap.add_argument("--nstep", type=int, default=None)
    ap.add_argument("--cfl", type=float, default=None, help="cfl / cfl_pseudo を差し替え (起動 A/B 用。既定はテンプレートの 2)")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    run = Path(a.run)
    if run.exists():
        raise SystemExit(f"{run} は既にある (run を使い回さない)")
    if a.kind == "cht" and not a.solid:
        raise SystemExit("cht には --solid が要る")
    if a.kind == "dry" and a.dT != 0.0:
        raise SystemExit("dry は dT=0 のみ")
    tdir = gc.HERE / "template"
    run.mkdir(parents=True)
    cfg = (tdir / ("solverConfig_dry.yaml" if a.kind == "dry" else f"solverConfig_{a.kind}.yaml")).read_text()
    if a.nstep is not None:
        cfg, n = re.subn(r"last: \{nStepOuter: \d+\}", f"last: {{nStepOuter: {a.nstep}}}", cfg)
        if n != 1:
            raise SystemExit("nStepOuter を差し替えられない")
    if a.cfl is not None:
        cfg, n = re.subn(r"cfl: [0-9.]+, cfl_pseudo: [0-9.]+,", f"cfl: {a.cfl}, cfl_pseudo: {a.cfl},", cfg)
        if n != 1:
            raise SystemExit("cfl を差し替えられない")
    (run / "solverConfig.yaml").write_text(cfg)
    shutil.copy(tdir / "probe.yaml", run / "probe.yaml")
    shutil.copy(a.mesh, run / "mesh.h5")
    x_out = gc.L_HEAT + gc.L_DOWN
    romin, romax, pmin, pmax = poiseuille_ic(run / "mesh.h5", x_out)
    p_in, ro_in = inlet_csv(run)
    tw = gc.T_IN + a.dT
    (run / "bcondConfig.yaml").write_text(BC.format(
        kind=a.kind, dT=a.dT, inl_ints="{inletProfile: 1}", ro=ro_in, um=gc.U_M, pin=p_in,
        pout=gc.P_OUT, ptb=gc.P_OUT + 0.5 * gc.RHO * gc.U_M ** 2, tin=gc.T_IN,
        heat_ints="{conjugate: 1}" if a.kind == "cht" else "", tw=tw))
    prov = [f"mesh.h5  <- {a.mesh} (sha256 {axcht.sha256(a.mesh)}; VALUE を Poiseuille IC にパッチ)",
            f"cfl_pseudo {a.cfl if a.cfl is not None else 'テンプレート既定'}",
            f"kind {a.kind}, dT {a.dT} K (T_w/T_c {tw} K), IC p {pmin:.3f}..{pmax:.3f} Pa, ro {romin:.6f}..{romax:.6f}",
            f"inlet_profile_1.csv: 放物 Ux (U_m {gc.U_M:.6f}), ro {ro_in:.8f}, Ps {p_in:.4f}",
            gc.summary()]
    if a.kind == "cht":
        # 固体の外面温度はソルバが solid.h5 の ROBIN/TC から読む (solidFem2d.cpp)。ΔT と食い違う固体を拒否する
        # (2026-09-27 codex plan レビュー m5)
        with h5py.File(a.solid, "r") as sf:
            tc = np.asarray(sf["ROBIN/TC"][:], float)
            ks = np.asarray(sf["SOLID/K_V"][:], float)
        if not np.allclose(tc, tw, rtol=0, atol=1e-9):
            shutil.rmtree(run)
            raise SystemExit(f"solid.h5 の ROBIN/TC ({tc.min()}..{tc.max()} K) が T_in+ΔT = {tw} K と違う (ΔT を取り違えた固体?)")
        shutil.copy(a.solid, run / "solid.h5")
        gj = Path(str(a.solid)[:-3] + ".grid.json")
        if gj.exists():
            shutil.copy(gj, run / "solid.grid.json")
        prov.append(f"solid.h5 <- {a.solid} (sha256 {axcht.sha256(a.solid)}; 読み出した ROBIN/TC {tc.min():.9g}..{tc.max():.9g} K, "
                    f"k_s {ks.tolist()})")
    if a.note:
        prov.append(f"note: {a.note}")
    (run / "RUN_INPUTS.txt").write_text("\n".join(prov) + "\n")
    print(f"[make_run] {run}\n  " + "\n  ".join(prov))


if __name__ == "__main__":
    main()
