#!/usr/bin/env python3
r"""case/64 の run ディレクトリを組む (実行はしない)。kind = dry (乾式 1 step、壁は等温 300 K) / cht (共役 A1/A2)。
IC は Poiseuille (放物速度・線形圧力・300 K)、入口は放物速度の inletProfile (case/63 と同じ作法)。
固体の ROBIN/TC・H と k_s を読み、登録値 (T_c = T_in + 10 K、h_o 570、k_s = 比 × k_f) と違えば停止する。

    python3 make_run.py dry run_0001_dry_r32 --mesh mesh/pipe_r32.h5
    python3 make_run.py cht run_00NN_a1_r32 --mesh mesh/pipe_r32.h5 --solid mesh/solid_run_0001_dry_r32_A1.h5 --case A1
"""
from __future__ import annotations
import argparse, re, shutil
from pathlib import Path
import h5py
import numpy as np
import pipe_common as pc
from pipe_common import axcht, gc

X_IN, X_OUT = -pc.L_UP, pc.L_HEAT + pc.L_DOWN
BC = """# case/64 ({kind}{case})。make_run.py が生成。
inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {{inletProfile: 1}},
          floats: {{ro: {ro:.10g}, Ux: {um:.10g}, Uy: 0.0, Uz: 0.0, Ps: {pin:.10g}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {pout}, Pt: {ptb:.10g}, Tt: {tin}}}}}
wall:   {{physID: 3, kind: wall_isothermal, outputHDFflg: 1, ints: {wints}, floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tin}}}}}
axis:   {{physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["dry", "cht"]); ap.add_argument("run")
    ap.add_argument("--mesh", required=True); ap.add_argument("--solid"); ap.add_argument("--case", choices=sorted(pc.KS_RATIO))
    ap.add_argument("--nstep", type=int); ap.add_argument("--out", type=int, help="outStepInterval (50 の倍数にしない)")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    run = Path(a.run)
    if run.exists():
        raise SystemExit(f"{run} は既にある")
    if a.kind == "cht" and not (a.solid and a.case):
        raise SystemExit("cht には --solid と --case が要る")
    if a.out and a.out % 50 == 0:
        raise SystemExit("outStepInterval は連成間隔 50 の倍数にしない (plan §4.7)")
    tdir = pc.HERE / "template"
    run.mkdir(parents=True)
    cfg = (tdir / ("solverConfig_dry.yaml" if a.kind == "dry" else "solverConfig_cht.yaml")).read_text()
    if a.nstep is not None:
        cfg, n = re.subn(r"last: \{nStepOuter: \d+\}", f"last: {{nStepOuter: {a.nstep}}}", cfg); assert n == 1
    if a.out is not None:
        cfg, n = re.subn(r"outStepInterval: \d+", f"outStepInterval: {a.out}", cfg); assert n == 1
    (run / "solverConfig.yaml").write_text(cfg)
    shutil.copy(tdir / "probe.yaml", run / "probe.yaml"); shutil.copy(a.mesh, run / "mesh.h5")
    # IC: Poiseuille + 線形圧力 + T_in
    with h5py.File(run / "mesh.h5", "a") as f:
        xyz = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3); g = f["VALUE"]
        x, r = xyz[:, 0], xyz[:, 1]
        u = 2 * gc.U_M * np.clip(1 - (r / pc.R) ** 2, 0, None)
        p = gc.P_OUT + 8 * gc.MU * gc.U_M / pc.R ** 2 * (X_OUT - x)
        ro = p / (gc.RGAS * gc.T_IN)
        for name, val in (("ro", ro), ("roUx", ro * u), ("roUy", 0 * u), ("roUz", 0 * u), ("roe", p / (gc.GAMMA - 1) + 0.5 * ro * u * u)):
            dt = g[name].dtype; del g[name]; g.create_dataset(name, data=np.asarray(val, dtype=dt))
    p_in = gc.P_OUT + 8 * gc.MU * gc.U_M / pc.R ** 2 * (X_OUT - X_IN); ro_in = p_in / (gc.RGAS * gc.T_IN)
    y = np.linspace(0, pc.R, 201)
    np.savetxt(run / "inlet_profile_1.csv", np.c_[y, np.full_like(y, ro_in), 2 * gc.U_M * (1 - (y / pc.R) ** 2), 0 * y, 0 * y, np.full_like(y, p_in)],
               header="y ro Ux Uy Uz Ps", comments="", fmt="%.12e")
    (run / "bcondConfig.yaml").write_text(BC.format(kind=a.kind, case=f" {a.case}" if a.case else "", ro=ro_in, um=gc.U_M, pin=p_in,
                                                     pout=gc.P_OUT, ptb=gc.P_OUT + 0.5 * gc.RHO * gc.U_M ** 2, tin=gc.T_IN,
                                                     wints="{conjugate: 1}" if a.kind == "cht" else ""))
    prov = [f"mesh.h5 <- {a.mesh} (sha256 {axcht.sha256(a.mesh)}; Poiseuille IC)", f"kind {a.kind} {a.case or ''}", gc.summary()]
    if a.kind == "cht":
        with h5py.File(a.solid, "r") as sf:
            tc = np.asarray(sf["ROBIN/TC"][:]); h = np.asarray(sf["ROBIN/H"][:]); ks = float(np.asarray(sf["SOLID/K_V"][:])[0])
        want_ks = pc.KS_RATIO[a.case] * gc.K_F
        if not (np.allclose(tc, gc.T_IN + pc.DT_C, atol=1e-9) and np.allclose(h, pc.H_O, atol=1e-9) and abs(ks - want_ks) < 1e-9 * want_ks):
            shutil.rmtree(run); raise SystemExit(f"solid.h5 の登録値が違う (TC {tc.min()}..{tc.max()}, H {h.min()}..{h.max()}, k_s {ks} / 期待 {want_ks})")
        shutil.copy(a.solid, run / "solid.h5")
        gj = Path(str(a.solid)[:-3] + ".grid.json")
        if gj.exists(): shutil.copy(gj, run / "solid.grid.json")
        prov.append(f"solid.h5 <- {a.solid} (sha256 {axcht.sha256(a.solid)}; ROBIN/TC {tc[0]} K、H {h[0]}、Robin 辺 {len(tc)}、k_s {ks:.6g})")
    if a.note: prov.append(f"note: {a.note}")
    (run / "RUN_INPUTS.txt").write_text("\n".join(prov) + "\n"); print(f"[make_run] {run}\n  " + "\n  ".join(prov))


if __name__ == "__main__":
    main()
