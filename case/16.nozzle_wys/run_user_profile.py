#!/usr/bin/env python3
"""ユーザ指定形状 (2026-08-19, mesh/nozzle_user_profile.py) の Wyslouzil ノズルを TP split
(MIXDRY=N2 擬似種 + H2O 独立種, thermoHrefTemp 298.15) で回す run の生成・投入。

  Euler (slip, visc=0) / NS (SST 低 Re no-slip 断熱壁, Sutherland + const-Pr) × 凝縮 off/on (Kw+HK)
  離散化: cell (押し出し 1 層 mesh/nozzle_user_2d.msh, frontback slip) / node (平面 mesh/nozzle_user_2d_planar[_inv].msh)
  IC: 1D 等エントロピー (壁 y(x) から A/A*, TP semi-perfect, datum 298.15)
  起動: soft (1 次 cfl0.5, 3000) → mid (1 次 cfl1, 3000) → 本段 (2 次 cfl_main, nsteps)
  node NS (2026-09-07, 最新の超音速ノズル node レシピ): 変換時 `mesh.nodeInletCornerWall: 1` (入口∩壁角の質量蓄積を根治),
    `nodeWallDirichlet: 1`, 低 Re SST + `katoLaunder: 1`, soft/mid 段は nStepInner 10。
  --mesh3d: 12.7 mm 押し出し 3D (mesh/nozzle_user_3d.msh, 4 壁 no-slip)。--ic-from RES.h5 で 2D 場を最近傍移植 (warm start)。

usage: python3 run_user_profile.py RUN_DIR --disc cell|node --phys euler|sst [--cond] [--cfl 2] [--nsteps 12000]
                                   [--prepare-only] [--mesh3d] [--ic-from RES.h5] [--stage-steps 3000]
"""
import subprocess, sys, shutil, os, re, argparse
from pathlib import Path
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))   # リポジトリ相対 (AWS でも動く)
import yaml, h5py
from forge_design.gas.semiperfect import mixture_pseudo_species_split, GasSemiPerfect
from forge_design.evaluate.ic import paste_isentropic_ic, zero_wall_velocity_ic
from forge_design.evaluate.runner import run_forge, FORGE_TOOLS, FORGE_BUILD, _ENV

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE / "mesh"))
from nozzle_user_profile import y_user_mm

ap = argparse.ArgumentParser()
ap.add_argument("run_dir")
ap.add_argument("--disc", choices=["cell", "node"], required=True)
ap.add_argument("--phys", choices=["euler", "sst"], required=True)
ap.add_argument("--cond", action="store_true")
ap.add_argument("--cfl", type=float, default=2.0)
ap.add_argument("--nsteps", type=int, default=12000)
ap.add_argument("--prepare-only", action="store_true")
ap.add_argument("--no-stage", action="store_true")
ap.add_argument("--mesh3d", action="store_true", help="12.7 mm 押し出し 3D メッシュ (4 壁 no-slip)")
ap.add_argument("--ic-from", default=None, help="warm start 元 res_*.h5 (interp_field.py 最近傍移植; 2D→3D 可)")
ap.add_argument("--stage-steps", type=int, default=3000, help="soft/mid 各段の step 数")
ap.add_argument("--out-interval", type=int, default=3000)
ap.add_argument("--cpg", action="store_true", help="診断用: CPG 単成分 (thermalMethod 0, species 無し)")
ap.add_argument("--laminar", action="store_true", help="診断用: NS 層流 (turbulence none)")
ap.add_argument("--reuse-h5", default=None, help="変換済み h5 を流用 (convert をスキップ; 同じ msh/config の run から)")
ap.add_argument("--msh", default=None, help="mesh/ 内の .msh 名を上書き (壁厚感度試験用)")
ap.add_argument("--cfg-sub", action="append", default=[], help="診断用: solverConfig 文字列置換 OLD=NEW (複数可)")
ap.add_argument("--bc-sub", action="append", default=[], help="診断用: bcondConfig 文字列置換 OLD=NEW (複数可)")
ap.add_argument("--soft-diag", type=int, default=0, help="診断用: soft 段の設定 (1 次 cfl0.5) だけを N step, 100 step 毎出力で回して終了")
a = ap.parse_args()
run_dir = CASE / a.run_dir
run_dir.mkdir(exist_ok=False)
MESH_H5 = "nozzle_user_3d.h5" if a.mesh3d else "nozzle_user_2d.h5"

# 1. species DB (N2 擬似種 MIXDRY + H2O)
db, Ys, order = mixture_pseudo_species_split({"N2": 0.98905, "H2O": 0.01095}, keep=("H2O",))
(run_dir / "species_db.yaml").write_text(yaml.safe_dump(db, sort_keys=False))
Y = [Ys[k] for k in order]

# 2. config (変換前に書く: node は discretization を convertGmshToForge が読む、SST は壁 no-slip で wall_dist)
euler = a.phys == "euler"
node = a.disc == "node"
cond = "" if not a.cond else """condensation:
  condensation: 1
  nCondSpecies: 1
  condModel: 1
  condGasSpecies: 1
  condKantrowitz: 1
  condGrowthModel: 0
"""
phys = ("  viscMethod: 0\n  visc: 0.0\n  thermCond: 0.0\n" if euler else
        "  viscMethod: 1\n  visc: 0.0000175\n  thermCond: 0.0241\n  thermCondMethod: 1\n  prandtlLam: 0.72\n")
turb = ('turbulence: {model: "none"}' if (euler or a.laminar) else
        'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInf: 1.0, omegaInf: 1000.0}')
# node NS: 入口∩壁の角ノード CV の入口側半割面を壁へ帰属 (変換時に焼き込み。run_0195 の壁ノード P>Pt 暴走の根治)
corner = "  nodeInletCornerWall: 1\n" if (node and not euler) else ""
species_cfg = "" if a.cpg else '  species: ["MIXDRY", "H2O"]\n  speciesDBFile: "species_db.yaml"\n  thermoHrefTemp: 298.15\n'
cfg = f"""mesh:
  meshFormat: "hdf5"
  discretization: "{a.disc}"
  isAxisymmetric: 0
  nodeWallDirichlet: 1
{corner}  meshFileName: "{MESH_H5}"
  valueFileName: "{MESH_H5}"
gpu: 1
solver: "SLAU"
physProp:
  isCompressible: 1
  thermalMethod: {0 if a.cpg else 2}
{phys}  ro: 1.2
  cp: 1039.0
  gamma: 1.4
{species_cfg}time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {a.nsteps}}}
  deltaT: {{control: 1, dt: 0.00001, cfl: {a.cfl}, cfl_pseudo: {a.cfl}, implicitRelax: 1.0, blockDPLUR: 1,
           dt_min: 0.00000001, dt_max: 1.0, detectNaN: 1}}
  outStepInterval: {a.out_interval}
  outStepStart: 0
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2}}
{turb}
initial: "nozzle_wys"
{cond}"""
for kv in a.cfg_sub:
    o, n = kv.split("=", 1); assert o in cfg, f"--cfg-sub: {o!r} not in config"; cfg = cfg.replace(o, n)
(run_dir / "solverConfig.yaml").write_text(cfg)
wall_kind = "slip" if euler else "wall"
species_bc = "" if a.cpg else f"    Y0: {Y[0]:.5f}\n    Y1: {Y[1]:.5f}\n"
bc = f"""inlet:
  physID: 1
  kind: inlet_Pressure
  outputHDFflg: 0
  ints:
  floats:
    Pt: 59070.0
    Tt: 286.65
{species_bc}    k: 1.0
    omega: 1000.0

outlet:
  physID: 2
  kind: outlet_statPress
  outputHDFflg: 0
  ints:
  floats:
    Ps: 2000.0
    Pt: 2000.0
    Tt: 286.65

wall:
  physID: 3
  kind: {wall_kind}          # Euler: slip / NS: no-slip 断熱壁
  outputHDFflg: 0
  ints:
  floats:
"""
if a.mesh3d:
    bc += """
sym:
  physID: 4
  kind: slip          # 半幅モデルの対称面 (z=6.35 mm)
  outputHDFflg: 0
  ints:
  floats:
"""
if not node and not a.mesh3d:
    bc += """
frontback:
  physID: 4
  kind: slip          # pseudo-2D 対称面
  outputHDFflg: 0
  ints:
  floats:
"""
for kv in a.bc_sub:
    o, n = kv.split("=", 1); assert o in bc, f"--bc-sub: {o!r} not in bcond"; bc = bc.replace(o, n)
(run_dir / "bcondConfig.yaml").write_text(bc)
shutil.copy(CASE / "run_0050_fig3_2d_sst_kwhk" / "probe.yaml", run_dir / "probe.yaml")

# 3. mesh 変換
msh = {("cell", True): "nozzle_user_2d.msh", ("cell", False): "nozzle_user_2d.msh",
       ("node", True): "nozzle_user_2d_planar_inv.msh", ("node", False): "nozzle_user_2d_planar.msh"}[(a.disc, euler)]
if a.mesh3d:
    msh = "nozzle_user_3d.msh"
if a.msh:
    msh = a.msh
if a.reuse_h5:
    shutil.copy(Path(a.reuse_h5).resolve(), run_dir / MESH_H5)
    (run_dir / "convert.log").write_text(f"reused {Path(a.reuse_h5).resolve()}\n")
    class _R: returncode = 0
    r = _R()
else:
    shutil.copy(CASE / "mesh" / msh, run_dir / msh)
    r = subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), msh, MESH_H5], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "convert.log").write_text(r.stdout + r.stderr)
# AWS (CUDA 13 ホスト) では convertGmshToForge が h5 を書き切った後の cudaFree で非零終了する (既知・無害,
# procedures/cloud-aws-gpu.md「移設で必ず踏む 4 点」③) → 成否は exit code でなく出力ファイルで判定する
if not (run_dir / MESH_H5).exists() or (run_dir / MESH_H5).stat().st_size < 1024:
    raise RuntimeError(f"convertGmshToForge failed (rc={r.returncode}); see convert.log")
# 品質検査: check_mesh_quality.py は node 変換 h5 の /VIZMESH (primal セル) を直接測れる (2026-09-08,
# ベクトル化済み) ので cell 再変換は不要
q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), str(run_dir / MESH_H5)],
                   capture_output=True, text=True)
(run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
print(q.stdout.strip().splitlines()[-1] if q.stdout.strip() else q.stderr[-300:])

# 4. IC: 1D 等エントロピー (平面: A/A* = y/y_t。paste_isentropic_ic は (r/r_t)^2 なので r=sqrt(y/y_t)*r_t を渡す)
class _Wall:
    x_throat, r_throat = 0.0, 0.0025
    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        y = np.array([y_user_mm(v * 1e3) * 1e-3 for v in np.ravel(x)]).reshape(x.shape)
        return np.sqrt(y / 0.0025) * 0.0025
gas = GasSemiPerfect({"N2": 0.98905, "H2O": 0.01095}, Tt=286.65)
info = (paste_isentropic_ic(str(run_dir / MESH_H5), _Wall(), 1.0, 59070.0, 286.65, 1.4, 1039.0, k_init=1.0, omega_init=1000.0)
        if a.cpg else
        paste_isentropic_ic(str(run_dir / MESH_H5), _Wall(), 1.0, 59070.0, 286.65, 1.4, 1039.0,
                            k_init=1.0, omega_init=1000.0, gas=gas, h_ref_T=298.15, species_Y=Y))
print("IC:", info)
if a.ic_from:
    # warm start: 過去 res を最近傍移植 (interp_field は (x,y) で照会するので 2D 場→3D メッシュにも使える)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"), str(Path(a.ic_from).resolve()), str(run_dir / MESH_H5)],
                   env=_ENV, check=True, capture_output=True, text=True)
    (run_dir / "IC_FROM.txt").write_text(str(Path(a.ic_from).resolve()) + "\n")
    print("IC interpolated from", a.ic_from)
if node and not euler:
    print("wall nodes zeroed:", zero_wall_velocity_ic(str(run_dir / MESH_H5)))
print("prepared", run_dir)
if a.prepare_only:
    sys.exit(0)

# 5. 段階起動 (soft 1次 cfl0.5 → mid 1次 cfl1 → 本段 2次 cfl_main)
cfg_main = (run_dir / "solverConfig.yaml").read_text()

def _stage(c, nsteps, label):
    c = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", c)
    c = re.sub(r"outStepInterval: \d+", f"outStepInterval: {nsteps}", c)
    (run_dir / "solverConfig.yaml").write_text(c)
    rc = run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if rc != 0 or not res or int(res[-1].stem.split("_")[1]) < nsteps:
        raise RuntimeError(f"{label} 段が失敗 (rc={rc}); res_nan_*.h5 を見る")
    # 同一メッシュの段間引き継ぎは index コピー (interp_field は (x,y) 最近傍なので 3D では z 列を混同し
    # k/ω が step 2 で爆発する — AWS run_0219 で実害, 2026-09-08。[[interp-field-coincident-nodes-trap]] と同趣旨)
    with h5py.File(res[-1], "r") as src, h5py.File(run_dir / MESH_H5, "r+") as dst:
        n = len(dst["VALUE/ro"])
        for k in src["VALUE"]:
            if k == "wall_dist" or k not in dst["VALUE"] or len(src["VALUE"][k]) != n:
                continue
            dst["VALUE"][k][:] = src["VALUE"][k][:]
    shutil.copy(run_dir / "residual_history.csv", run_dir / f"residual_history_{label}.csv")
    for f in run_dir.glob("res_*"):
        f.unlink()
    print(f"[{label}] done")

if a.soft_diag:
    inner = "nStepInner: 10" if (node and not euler) else "nStepInner: 5"
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", cfg_main).replace("convMethod: 1, limiter: 2", "convMethod: 0, limiter: 0").replace("nStepInner: 5", inner)
    soft = re.sub(r"nStepOuter: \d+", f"nStepOuter: {a.soft_diag}", soft); soft = re.sub(r"outStepInterval: \d+", "outStepInterval: 100", soft)
    (run_dir / "solverConfig.yaml").write_text(soft)
    print("soft-diag rc", run_forge(run_dir)); sys.exit(0)
if not a.no_stage:
    inner = "nStepInner: 10" if (node and not euler) else "nStepInner: 5"   # node NS の起動レシピ (case/40 y+1)
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", cfg_main).replace("convMethod: 1, limiter: 2", "convMethod: 0, limiter: 0").replace("nStepInner: 5", inner)
    _stage(soft, a.stage_steps, "soft")
    mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 1.0, cfl_pseudo: 1.0", cfg_main).replace("convMethod: 1, limiter: 2", "convMethod: 0, limiter: 0").replace("nStepInner: 5", inner)
    _stage(mid, a.stage_steps, "mid")
(run_dir / "solverConfig.yaml").write_text(cfg_main)
rc = run_forge(run_dir)
print("main rc", rc)
print((run_dir / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-1] if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else "")
