"""バッチ評価: 問題定義 YAML + dv → run ディレクトリ準備 → forge 実行 → 判定・抽出。

パイプライン (Phase 0, thruster_bell):
  1. NozzleWall 構築 + validate (幾何フィルタ)
  2. 構造化 quad メッシュ → nozzle.msh (msh4.1 直書き)
  3. solverConfig/bcondConfig/probe を生成
  4. convertGmshToForge → nozzle.h5 (wall_dist・幾何量は変換器が計算)
  5. check_mesh_quality.py ゲート (FAIL → 評価失敗)
  6. 準 1D 等エントロピー IC 貼り付け
  7. run_case.sh (forge 実行 + 収束判定 + 残差 png)  ※hook 準拠の起動経路
  8. CONVERGENCE_VERDICT / metrics.json
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.wall import NozzleWall, WallParams
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic

# リポジトリ位置から導く (AWS など別マシンでも動くように。FORGE_ROOT で上書き可)
FORGE_ROOT = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
FORGE_TOOLS = FORGE_ROOT / "solver_density_cuda" / "tools"
FORGE_BUILD = FORGE_ROOT / "solver_density_cuda" / "build"
_ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")


def build_wall(p: Problem) -> tuple:
    g = p.geometry
    theta_a = np.deg2rad(dv_value(p, "theta_a_deg"))
    L_div = dv_value(p, "L_over_rt")
    # theta_e は dv 宣言があれば dv (TOP 幾何 dv — 親計画 §4.6③)、なければ geometry 固定値
    theta_e_deg = (dv_value(p, "theta_e_deg") if "theta_e_deg" in p.dv
                   else float(g.get("theta_e_deg", 8.0)))
    wp = WallParams(
        r_inlet=float(g.get("r_inlet", 2.5)),
        L_pipe=float(g.get("L_pipe", 0.5)),
        L_contract=float(g.get("L_contract", 3.0)),
        Ru=float(g.get("Ru", 1.5)),
        Rd=float(g.get("Rd", 0.4)),
        phi_u=np.deg2rad(float(g.get("phi_u_deg", 25.0))),
        theta_a=theta_a,
        L_div=L_div,
        r_exit=float(np.sqrt(p.spec["eps"])),
        theta_e=np.deg2rad(theta_e_deg),
        bell_type=str(g.get("bell_type", "hermite")),
    )
    wall = NozzleWall(wp)
    msgs = wall.validate()
    if msgs:
        raise ValueError("幾何フィルタ不合格: " + "; ".join(msgs))
    return wall, wp


def _solver_config(p: Problem, nsteps: int, out_interval: int) -> str:
    """離散化別の solverConfig を生成する。

    node は実証レシピに従う (plans/active/boundary-node-nozzle-wall-outlet-stability.md):
    **全域 2 次 (convMethod 1, bndFirstOrder なし) + 陰解法 (timeIntegration 11,
    blockDPLUR, cfl_pseudo 4)** + nodeWallDirichlet/axisCentroidShift (nodeAxisDirichlet は 2026-08-16 撤去) +
    [2026-08-16 更新: `nodeAxisDirichlet` は撤去。node 既定が「値=ノード座標 + 軸 u_r 三点セット + エッジ中点再構成 + LSQ」
     になり (case/43)、軸ノードは通常 DOF として解く。]
    katoLaunder。**収束場からの warm start (--warm-from) が前提** (冷間 IC は発散)。
    壁第一セルは細かすぎ禁物 (wall_first_frac ≳ 5e-3)。かつての 2 次・陰解法発散は
    いずれも軸行真空化が種で、旧 nodeAxisDirichlet で対症 (現在は値=ノード整合セットで根治) (bndFirstOrder による境界 1 次化は
    軸病変のマスキングだったため撤去。同 plan §2.6–2.7, case/40 run_0009–0022)。
    残差は入口/収縮部の微小リミットサイクルで ~1e-7 プラトー (cell 全域 2 次と同格)、
    計量は quasisteady STEADY を確認して使う。
    """
    sst = p.evaluate.get("turbulence", "sst") == "sst"
    disc = p.mesh.get("discretization", "node")  # 既定 node — ユーザ方針 2026-08-03
    if disc == "node":
        node_keys = ', isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1'
        tint, cfl, conv, dplur, kato = 11, 4.0, 1, 1, ", katoLaunder: 1"
    else:
        node_keys = ", isAxisymmetric: 1"
        tint, cfl, conv, dplur, kato = 11, 4.0, 1, 1, ""
    turb = (f'turbulence: {{model: "sst", scalarDiffusion: 1, dilatationCorrection: 2{kato}}}'
            if sst else 'turbulence: {model: "none"}')
    return f"""mesh: {{meshFormat: "hdf5", discretization: "{disc}"{node_keys}, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 1,
           ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: {dplur}, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_interval}
  timeIntegration: {tint}
  nStepInner: 5
space: {{convMethod: {conv}, limiter: 2}}
{turb}
initial: "uniform_p101325_u10"
"""


def _bcond_config(p: Problem) -> str:
    Pt = float(p.spec["Pt"])
    Tt = float(p.spec["Tt"])
    pa = float(p.spec["p_ambient"])
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{Pt: {Pt}, Tt: {Tt}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}
wall:   {{physID: 3, kind: wall,             outputHDFflg: 1, ints: , floats: }}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


PROBE_STUB = "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n"


def prepare(problem_path, run_dir, nsteps=None, warm_from=None) -> dict:
    """メッシュ・config・IC まで準備する (forge は起動しない)。

    warm_from: 過去の収束 res_*.h5 (cell/node 不問)。指定時は準 1D IC の上から
    `interp_field.py` (最近傍 cross-mesh restart) で場を移植する。**node は必須**
    (冷間 IC は発散 — 実証レシピ参照)。
    """
    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    wall, wp = build_wall(p)
    scale = float(p.spec["r_throat"])

    mp = Mesh2DParams(
        ni=int(p.mesh.get("ni", 241)),
        nj=int(p.mesh.get("nj", 81)),
        wall_first_frac=float(p.mesh.get("wall_first_frac", 2.0e-3)),
        throat_refine=float(p.mesh.get("throat_refine", 3.0)),
        scale=scale,
    )
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 壁座標も記録 (IC・後処理・帰還エンジンの入力)
    xs = np.linspace(wall.x_in, wall.x_e, 800)
    np.savetxt(run_dir / "contour.csv", np.c_[xs * scale, wall.r(xs) * scale],
               delimiter=",", header="x_m,r_m", comments="")

    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    (run_dir / "solverConfig.yaml").write_text(_solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)

    disc = p.mesh.get("discretization", "node")
    if disc == "node":
        # 品質ゲートは primal (cell 変換) で測る: 双対 h5 は品質ツールが読めないし、
        # 品質の実体は primal quad の AR/skew。cell で一時変換 → 検査 → node で本変換。
        (run_dir / "solverConfig.yaml").write_text(
            (run_dir / "solverConfig.yaml").read_text().replace(
                f'discretization: "{disc}"', 'discretization: "cell"').replace(
                ", axisCentroidShift: 1, nodeWallDirichlet: 1", ""))
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                       cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
        q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "nozzle_qc.h5"],
                           cwd=run_dir, env=_ENV, capture_output=True, text=True)
        (run_dir / "MESH_QUALITY.txt").write_text("(primal/cell 変換で検査)\n" + q.stdout + q.stderr)
        (run_dir / "nozzle_qc.h5").unlink()
        for f in run_dir.glob("nozzle_qc.xmf"):
            f.unlink()
        if q.returncode != 0:
            raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
        (run_dir / "solverConfig.yaml").write_text(_solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                       cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    else:
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                       cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
        q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "nozzle.h5"],
                           cwd=run_dir, env=_ENV, capture_output=True, text=True)
        (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
        if q.returncode != 0:
            raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")

    exit1d = paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                                 float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    if warm_from is not None:
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(warm_from), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        if disc == "node":
            _smooth_wall_entropy(run_dir / "nozzle.h5", p.gamma, p.cp)
    elif disc == "node":
        raise ValueError("node は warm start 必須 (--warm-from に収束 res_*.h5 を指定。"
                         "冷間 IC は発散 — plans/active/boundary-node-nozzle-wall-outlet-stability.md)")
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n,
            "scale_m": scale, "exit_1d": exit1d,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "cells": int(quads.shape[0])}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def _smooth_wall_entropy(h5path, gamma: float, cp: float, npass: int = 10) -> None:
    """cell→node 最近傍 interp が壁ノード列に刻む等圧エントロピー市松 (T/ρ 逆相ジグザグ) を除去する。

    壁ノードは u=0 で移流平滑化が働かず、P は音響で均されるが T/ρ の市松は接線熱伝導でしか
    消えない (実測: 静的凍結、ベル部 |ΔΔT|~500K@interp直後 → 収束後も ~200K 残存)。
    ここで壁列の T を x ソートで 1-2-1 フィルタし、P 不変で ro/roe を再構成する
    (u は壁 Dirichlet が 0 にするため運動量も 0 化)。plan
    boundary-node-nozzle-wall-outlet-stability §2.8 参照。
    """
    import h5py

    R = cp * (gamma - 1.0) / gamma
    with h5py.File(h5path, "r+") as f:
        wall = None
        for bid in f["BCONDS"]:
            kind = f["BCONDS"][bid].attrs.get("bcondKind", b"")
            kind = kind.decode() if isinstance(kind, bytes) else str(kind)
            if kind in ("wall", "wall_isothermal"):
                wall = np.array(sorted(set(int(i) for i in f["BCONDS"][bid]["iCells"][:])))
        if wall is None or wall.size < 5:
            return
        cc = f["CELLS/centCoords"][:].reshape(-1, 3)
        order = np.argsort(cc[wall, 0])
        w = wall[order]
        v = f["VALUE"]
        ro = v["ro"][:]; roe = v["roe"][:]
        rux = v["roUx"][:]; ruy = v["roUy"][:]; ruz = v["roUz"][:]
        ke = 0.5 * (rux[w]**2 + ruy[w]**2 + ruz[w]**2) / np.maximum(ro[w], 1e-30)
        P = (gamma - 1.0) * (roe[w] - ke)
        T = P / np.maximum(ro[w] * R, 1e-30)
        for _ in range(npass):
            T[1:-1] = 0.25 * T[:-2] + 0.5 * T[1:-1] + 0.25 * T[2:]
        ro_new = P / np.maximum(R * T, 1e-30)
        roe_new = P / (gamma - 1.0)  # u=0 (壁 Dirichlet と整合、運動量も 0 化)
        ro[w] = ro_new; roe[w] = roe_new
        rux[w] = 0.0; ruy[w] = 0.0; ruz[w] = 0.0
        v["ro"][:] = ro; v["roe"][:] = roe
        v["roUx"][:] = rux; v["roUy"][:] = ruy; v["roUz"][:] = ruz


def run_forge(run_dir) -> int:
    """run_case.sh 経由で forge を実行 (収束判定・残差 png も生成される)。"""
    r = subprocess.run([str(FORGE_TOOLS / "run_case.sh"), str(Path(run_dir).resolve())],
                       env=_ENV, capture_output=True, text=True)
    (Path(run_dir) / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def collect(problem_path, run_dir) -> dict:
    """収束 VERDICT + 最終 res からメトリクスを抽出し metrics.json に保存。"""
    from ..metrics.extract import thrust_metrics

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text() if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else "(no verdict)"
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    out = {"convergence_verdict": verdict.strip().splitlines()[-2:] if verdict else []}
    if res:
        out["res_file"] = res[-1].name
        out.update(thrust_metrics(run_dir / "nozzle.h5", res[-1], p.gamma,
                                  float(p.spec["Pt"]), float(p.spec["p_ambient"]),
                                  float(p.spec["r_throat"])))
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse

    ap = argparse.ArgumentParser(description="forge_design バッチ評価 (Phase 0)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--warm-from", default=None,
                    help="収束 res_*.h5 から cross-mesh warm start (node は必須)")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, a.steps, warm_from=a.warm_from)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    rc = run_forge(a.run_dir)
    out = collect(a.problem, a.run_dir)
    print(json.dumps(out, indent=1))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
