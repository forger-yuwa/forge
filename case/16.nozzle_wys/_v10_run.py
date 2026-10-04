#!/usr/bin/env python3
"""#10: case/16 を node + no-slip (層流 NS) で段階起動し、flag0/flag1 の A/B を作る。

**入口は Pt 101325 / Tt 293.15** — 組込み初期場 `nozzle_wys` (setInitial.hpp:370, P0/T0 がハードコード) と
**総圧・総温を厳密に揃える**ため。V3 は flag0/flag1 の同一設定どうしの比較なので作動点の絶対値は無関係で、
IC/BC 不整合 (divergence-and-startup の発散主因の筆頭) を完全に消せる。
初回は `uniform_p101325_u10` (入口全圧の 1.7 倍の圧力) で組んで層流段から残差が rising、2 次段で発散した。

plan convection-slau-wall-normal-chi §5.1 #10 / §6 V2・V3。
既存の node run は全部 Euler (slip 壁) で chi_n の対象外なので、**新規に NS 設定を組む**。
段階: Euler slip (1次) -> 層流 no-slip (1次) -> 2次ランプ -> 本段。
出口は **outflow** (このノズルの出口は超音速。outlet_statPress は node の亜音速壁列に背圧を課す)。

usage: _v10_run.py --mesh M.h5 --out RUNDIR --flag 0|1 [--main-steps 12000]
"""
import argparse, os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH",""),
           FORGE_CUDA_BLOCKSIZE="128")

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: {visc}, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1039.0, gamma: 1.4, pMin: 20.0}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {n}}}
  deltaT: {{control: 1, dt: 1e-6, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: 1.0, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-10, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {oi}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {conv}, limiter: {lim}, slauWallNormalChi: {flag}}}
turbulence: {{model: "none"}}
output: {{level: 1}}
initial: "nozzle_wys"
"""
BC = """inlet:  {{physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {{Pt: 101325.0, Tt: 293.15, k: 1.0, omega: 1000.0}}}}
outlet: {{physID: 2, kind: outflow, outputHDFflg: 0, ints: , floats: }}
wall:   {{physID: 3, kind: {wk}, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
"""

def stage(rd, tag, n, cfl, conv, lim, visc, wk, flag, oi):
    (rd/"solverConfig.yaml").write_text(CFG.format(n=n, cfl=cfl, conv=conv, lim=lim, visc=visc, flag=flag, oi=oi))
    (rd/"solverConfig_main.yaml").write_text((rd/"solverConfig.yaml").read_text())
    (rd/"bcondConfig.yaml").write_text(BC.format(wk=wk))
    r = subprocess.run([str(TOOLS/"run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd/f"run_{tag}.log").write_text(r.stdout + r.stderr)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    last = int(res[-1].stem.split("_")[1]) if res else -1
    print(f"  [{tag}] rc={r.returncode} last_dump={last}/{n}")
    if r.returncode != 0 or last < n:
        print((rd/f"run_{tag}.log").read_text()[-1200:]); raise SystemExit(f"stage {tag} 失敗")
    if tag != "main":
        subprocess.run([sys.executable, str(TOOLS/"restart_field.py"), str(res[-1]), str(rd/"m.h5")],
                       env=ENV, check=True, capture_output=True, text=True)
        for f in rd.glob("res_*"): f.unlink()
        for f in ("residual_history.csv","residual_history.png","CONVERGENCE_VERDICT.txt","forge_run.log"):
            if (rd/f).exists(): shutil.move(str(rd/f), str(rd/f"{Path(f).stem}_{tag}{Path(f).suffix}"))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--flag", type=int, required=True); ap.add_argument("--main-steps", type=int, default=12000)
    a = ap.parse_args()
    rd = Path(a.out)
    if rd.exists(): shutil.rmtree(rd)
    rd.mkdir(parents=True)
    shutil.copy(a.mesh, rd/"m.h5")
    (rd/"probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    print(f"== {rd} (flag {a.flag}) ==")
    # **1 段で 1 つだけ変える** (初回は 2 次化と CFL 倍増を同時にやって ramp 段で発散した)
    stage(rd, "euler", 6000, 0.3, 0, 0, 0, "slip", a.flag, 6000)   # 非粘性 slip 暖機 (1 次)
    stage(rd, "lam",   6000, 0.3, 0, 0, 1, "wall", a.flag, 6000)   # 壁を no-slip に (CFL 据え置き)
    stage(rd, "ramp1", 4000, 0.3, 1, 2, 1, "wall", a.flag, 4000)   # 2 次化 (CFL 据え置き)
    stage(rd, "ramp2", 4000, 0.6, 1, 2, 1, "wall", a.flag, 4000)   # CFL を上げる
    stage(rd, "main",  a.main_steps, 1.0, 1, 2, 1, "wall", a.flag, max(a.main_steps//24, 1))
