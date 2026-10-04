#!/usr/bin/env python3
"""1 step の場がビット再現するかを測る (plan convection-slau-wall-normal-chi #10d の前提確認)。

同一メッシュ・同一 IC・同一設定・同一ブロックサイズで 1 step を 2 回回し、res_1.h5 をビット比較する。
node の残差 gather は atomicAdd なので、これが同一でなければ「非対象面のビット不変」は測定できない。
"""
import shutil, subprocess, sys, os
from pathlib import Path
import h5py, numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_runs as G

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh" / "fp_y1_12um.h5"
TOOLS = G.TOOLS

def make(rd, flag):
    if rd.exists(): shutil.rmtree(rd)
    rd.mkdir(parents=True)
    shutil.copy(MESH, rd / "mesh.h5")
    G.patch_ic(rd / "mesh.h5", None)
    txt = G.cfg(nsteps=1, cfl=0.2, relax=1.0, conv=0, lim=0, ninner=5, outint=1, plain=True, lam=True)
    if flag is not None:
        txt = txt.replace("space: {convMethod: 0, limiter: 0}",
                          f"space: {{convMethod: 0, limiter: 0, slauWallNormalChi: {flag}}}")
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    (rd / "solverConfig.yaml").write_text(txt)
    (rd / "bcondConfig.yaml").write_text(G.BC.format(ro=G.RO_INF, u=G.U_INF, p=G.P_INF, t=G.T_INF,
                                                     k=G.K_INF, om=G.OM_INF, wallline=G.wall_line("adiabatic")))

def run(rd):
    env = dict(G.ENV, FORGE_CUDA_BLOCKSIZE="128")
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=env, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode

def cmp(a, b, label):
    with h5py.File(a) as fa, h5py.File(b) as fb:
        keys = sorted(k for k in fa["/VALUE"] if fa[f"/VALUE/{k}"].shape == fb[f"/VALUE/{k}"].shape)
        print(f"--- {label} ---")
        allsame = True
        for k in keys:
            x = fa[f"/VALUE/{k}"][:]; y = fb[f"/VALUE/{k}"][:]
            same = np.array_equal(x.view(np.uint8), y.view(np.uint8))
            ndiff = int((x != y).sum())
            if not same: allsame = False
            nan = int((~np.isfinite(x)).sum()) + int((~np.isfinite(y)).sum())
            print(f"  {k:12s} bit-identical={str(same):5s} ndiff={ndiff:7d}/{x.size} nonfinite={nan}")
        print(f"  => ALL BIT-IDENTICAL: {allsame}")
        return allsame

if __name__ == "__main__":
    dirs = {}
    for name, flag in (("run_0030_bitrep_a", None), ("run_0031_bitrep_b", None)):
        rd = HERE / name; make(rd, flag); rc = run(rd)
        res = rd / "res_1.h5"
        print(f"{name}: rc={rc} res_1={res.exists()}")
        if rc != 0 or not res.exists():
            print((rd / "run_case_stdout.log").read_text()[-1500:]); raise SystemExit(1)
        dirs[name] = res
    cmp(dirs["run_0030_bitrep_a"], dirs["run_0031_bitrep_b"], "同一設定の反復 (flag 未指定 x2), 1 step")
