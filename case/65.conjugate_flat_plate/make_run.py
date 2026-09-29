#!/usr/bin/env python3
r"""case/65 の run ディレクトリ (kind = dry / cht)。IC は一様流 (U_∞・300 K・101325 Pa)、入口は一様速度 (エントロピー固定)。
固体の ROBIN/TC・H・k_s を登録値 (T_h = T_∞ + 10 K、h 1e8、k_s = 比 × k_f) と照合する。

    python3 make_run.py dry run_0001_dry_n32 --mesh mesh/plate_n32.h5
    python3 make_run.py cht run_00NN_c1_n32 --mesh mesh/plate_n32.h5 --solid mesh/solid_run_0002_dry_n32_C1.h5 --case C1
"""
from __future__ import annotations
import argparse, re, shutil
from pathlib import Path
import h5py
import numpy as np
import plate_common as pc
from plate_common import axcht, gc

BC = """# case/65 ({kind}{case})。make_run.py が生成。
inlet:     {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.10g}, Ux: {U:.10g}, Uy: 0.0, Uz: 0.0, Ps: {p:.10g}}}}}
outlet:    {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.10g}, Pt: {pt:.10g}, Tt: {tt:.10g}}}}}
top:       {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
slip_up:   {{physID: 4, kind: slip, outputHDFflg: 0, ints: , floats: }}
plate:     {{physID: 5, kind: wall_isothermal, outputHDFflg: 1, ints: {wints}, floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tin}}}}}
slip_down: {{physID: 6, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["dry", "cht"]); ap.add_argument("run")
    ap.add_argument("--mesh", required=True); ap.add_argument("--solid"); ap.add_argument("--case", choices=sorted(pc.KS_RATIO))
    ap.add_argument("--nstep", type=int); ap.add_argument("--out", type=int); ap.add_argument("--note", default="")
    a = ap.parse_args()
    run = Path(a.run)
    if run.exists(): raise SystemExit(f"{run} は既にある")
    if a.kind == "cht" and not (a.solid and a.case): raise SystemExit("cht には --solid と --case が要る")
    if a.out and a.out % 50 == 0: raise SystemExit("outStepInterval は 50 の倍数にしない")
    tdir = pc.HERE / "template"; run.mkdir(parents=True)
    cfg = (tdir / ("solverConfig_dry.yaml" if a.kind == "dry" else "solverConfig_cht.yaml")).read_text()
    if a.nstep is not None:
        cfg, n = re.subn(r"last: \{nStepOuter: \d+\}", f"last: {{nStepOuter: {a.nstep}}}", cfg); assert n == 1
    if a.out is not None:
        cfg, n = re.subn(r"outStepInterval: \d+", f"outStepInterval: {a.out}", cfg); assert n == 1
    (run / "solverConfig.yaml").write_text(cfg)
    shutil.copy(pc.HERE.parent / "63.graetz_cht" / "template" / "probe.yaml", run / "probe.yaml"); shutil.copy(a.mesh, run / "mesh.h5")
    ro = gc.P_OUT / (gc.RGAS * gc.T_IN)
    with h5py.File(run / "mesh.h5", "a") as f:
        g = f["VALUE"]; n = g["ro"].shape[0]
        for name, val in (("ro", ro), ("roUx", ro * pc.U_INF), ("roUy", 0.0), ("roUz", 0.0), ("roe", gc.P_OUT / (gc.GAMMA - 1) + 0.5 * ro * pc.U_INF ** 2)):
            dt = g[name].dtype; del g[name]; g.create_dataset(name, data=np.full(n, val, dtype=dt))
    pt = gc.P_OUT + 0.5 * ro * pc.U_INF ** 2; tt = gc.T_IN + pc.U_INF ** 2 / (2 * gc.CP)
    (run / "bcondConfig.yaml").write_text(BC.format(kind=a.kind, case=f" {a.case}" if a.case else "", ro=ro, U=pc.U_INF, p=gc.P_OUT, pt=pt, tt=tt,
                                                     tin=gc.T_IN, wints="{conjugate: 1}" if a.kind == "cht" else ""))
    prov = [f"mesh.h5 <- {a.mesh} (sha256 {axcht.sha256(a.mesh)}; 一様流 IC)", f"kind {a.kind} {a.case or ''}", f"U {pc.U_INF:.6f} m/s, L {pc.L*1e3:.4f} mm, b {pc.B*1e3:.4f} mm"]
    if a.kind == "cht":
        with h5py.File(a.solid, "r") as sf:
            tc = np.asarray(sf["ROBIN/TC"][:]); h = np.asarray(sf["ROBIN/H"][:]); ks = float(np.asarray(sf["SOLID/K_V"][:])[0])
        want = pc.KS_RATIO[a.case] * gc.K_F
        if not (np.allclose(tc, gc.T_IN + pc.DT_H, atol=1e-9) and np.allclose(h, pc.H_BACK) and abs(ks - want) < 1e-9 * want):
            shutil.rmtree(run); raise SystemExit("solid.h5 の登録値が違う")
        shutil.copy(a.solid, run / "solid.h5")
        prov.append(f"solid.h5 <- {a.solid} (sha256 {axcht.sha256(a.solid)}; ROBIN/TC {tc[0]}、H {h[0]:g}、辺 {len(tc)}、k_s {ks:.6g})")
    if a.note: prov.append(f"note: {a.note}")
    (run / "RUN_INPUTS.txt").write_text("\n".join(prov) + "\n"); print(f"[make_run] {run}\n  " + "\n  ".join(prov))


if __name__ == "__main__":
    main()
