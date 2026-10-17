"""wall-driven チェーン W3: 多項式スロート域壁 + θ_D 円錐延長ドメインの node Euler 評価。

plan: plans/active/tooling-nozzle-walldriven-chain.md W3。既存経路 (runner_wt の
モード F チェーン) には触れない。W4 はこの run の場から D 発の C⁻ をデータ線として
抽出する — 円錐延長は「D の C⁻ より上流の場は延長壁に依存しない」(超音速の依存域)
ことを利用した打ち切りで、θ_D 一定なので D で新たな膨張も作らない。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_walldriven \
      case/41.wind_tunnel_design/problem_m4_walldriven.yaml run_dir [--prepare-only]
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.wall_walldriven import WallDrivenCFDWall, WallDrivenThroatRegion
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import run_converter, FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge
from .runner_wt import _bcond, _config_euler, _config_euler_node


def build_wall(p: Problem) -> tuple[WallDrivenCFDWall, dict]:
    """dv (theta_D, R_t, L_D_frac) + geometry → CFD 壁。決定的。

    L_D は独立入力にせず `L_D_frac × 3 R_t tanθ_D` (変曲なし条件の上限比) で与える。
    円錐長は `L_ext_factor × r_D` (既定 4.0 — D 発 C⁻ の軸着地を覆う。着地位置は
    run 後に実測で検証し、足りなければ factor を上げて作り直す)。終端は既定で
    円錐のまま (`L_turn`/`L_cyl` は既定 0) — 段階起動さえすれば戻しテーパは
    不要と run_0043 で確定済み (WallDrivenCFDWall docstring 参照)。
    """
    theta_D = np.deg2rad(float(dv_value(p, "theta_D")))
    R_t = float(dv_value(p, "R_t"))
    L_D_frac = float(dv_value(p, "L_D_frac", 0.9))
    L_D = L_D_frac * 3.0 * R_t * np.tan(theta_D)
    region = WallDrivenThroatRegion(
        r_U=float(p.geometry.get("r_inlet", 2.5)), R_t=R_t,
        L_U=float(p.geometry.get("L_U", 3.5)), L_D=L_D, theta_D=theta_D)
    v = region.validate()
    if v:
        raise ValueError("wall-driven 壁フィルタ不合格: " + "; ".join(v))
    wall = WallDrivenCFDWall(
        region, L_pipe=float(p.geometry.get("L_pipe", 0.5)),
        L_cone=float(p.geometry.get("L_ext_factor", 4.0)) * region.dn.r_D,
        L_turn=float(p.geometry.get("L_turn", 0.0)),
        L_cyl=float(p.geometry.get("L_cyl", 0.0)))
    diag = region.diagnostics()
    diag.update({"theta_D_deg": float(np.rad2deg(theta_D)), "R_t": R_t,
                 "L_D": L_D, "L_D_frac": L_D_frac, "L_U": region.up.L_U,
                 "x_D": wall.x_D, "x_turn": wall.x_turn, "x_cyl": wall.x_cyl,
                 "x_in": wall.x_in, "x_e": wall.x_e})
    return wall, diag


def prepare(problem_path, run_dir, nsteps=None) -> dict:
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_walldriven":
        raise ValueError("runner_walldriven は wind_tunnel_axisym_walldriven 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    wall, diag = build_wall(p)
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 幾何記録 (W4 のデータ線抽出・診断の入力)
    xs = np.linspace(wall.x_in, wall.x_e, 1200)
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[xs * scale, wall.r(xs) * scale, wall.theta(xs), wall.kappa(xs)],
               delimiter=",", header="x_m,r_m,theta_rad,kappa_per_rt", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 16000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 4, 1)))
    (run_dir / "bcondConfig.yaml").write_text(_bcond(p, euler=True))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質検査は cell 変換の一時コピーで (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(_config_euler(p, n, out_int, 4.0, 1))
    run_converter([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"], cwd=run_dir, env=_ENV)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node 変換 (node config を書いてから変換 — 必須) + 等エントロピー IC
    (run_dir / "solverConfig.yaml").write_text(_config_euler_node(p, n, out_int, 4.0, 1))
    run_converter([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"], cwd=run_dir, env=_ENV)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    info = {"chain": "walldriven_w3", "discretization": "node", **diag,
            "r_D": diag["r_D"], "nStepOuter": n, "scale_m": scale,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(run_dir, problem_path=None) -> int:
    """soft 段 (1 次 + cfl0.5, 3000 step) → 本段 (runner_wt.run_staged と同方式)。

    cold 単段 (conv1+cfl4 直投入) は run_0040/0041 で step ~9 発散 (D 発膨張波の
    軸反射が壁へ戻る位置で準 1D IC との不整合過渡が最大になる)。**段階起動だけで
    十分** (run_0043: 円錐のみ・戻しテーパなしでも段階起動なら NaN なし完走)。
    divergence-and-startup の標準手順 (易しい条件で立ち上げ→引き継ぎ) に従う。

    段間は段終了ゲート + `restart_field.py` (`runner_wt._staged_handoff`; 2026-10-14 まで `interp_field.py`)。段ごとの実効設定を
    `stage_manifest.json` に、段の残差履歴を `residual_history_<tag>.csv` (tag = soft / main) に残す (plan
    tooling-generator-defaults-and-launch-gate §5.1 #9)。段の config の文字列は従来と同じ。
    2026-10-15 (plan §5.1 #11): soft 段の値は YAML 上の位置で書き換え (`config_builder.set_scalars`; 旧: 文字列置換・正規表現、
    本段からの文字列は同じ)、各段は明示すべきキーを足して宣言を書く (機種 walldriven、`problem_path` は逸脱の出所の記録用)。
    """
    from . import config_builder as CB
    from .runner_axismach import _stage_manifest_cls
    from .runner_wt import _staged_handoff, _staged_launch
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    bc_path = run_dir / "bcondConfig.yaml"
    bc_text = bc_path.read_text() if bc_path.exists() else ""
    sm = _stage_manifest_cls()(run_dir)
    recipe = {"family": "walldriven", "problem": problem_path, "generator": "runner_walldriven.run_staged"}
    cfg_soft = CB.set_scalars(cfg_main, {"time.deltaT.cfl": "0.5", "time.deltaT.cfl_pseudo": "0.5", "space.convMethod": "0",
                                         "time.last.nStepOuter": "3000", "time.outStepInterval": "3000"})
    # 計算の後の判定 (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12): 仕様は最初の起動の前に、総合は本段の後に 1 度だけ
    from . import judge as J
    J.write_spec(run_dir, J.spec_snapshot("runner_walldriven.run_staged"))
    with J.Staged(run_dir, "runner_walldriven.run_staged"):
        rc = _staged_launch(run_dir, sm, bc_text, cfg_soft, "soft", launch=run_forge, recipe=recipe)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < 3000:
            raise RuntimeError("soft 段が失敗 (発散切り分けは res_nan_*.h5 を見る)")
        _staged_handoff(run_dir, res[-1], cfg_soft, "soft")
        for f in run_dir.glob("res_*"):
            f.unlink()
        rc = _staged_launch(run_dir, sm, bc_text, cfg_main, "main", launch=run_forge, recipe=recipe)
    return rc


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="wall-driven W3 評価 (node Euler)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, nsteps=a.steps)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    rc = run_staged(a.run_dir, problem_path=a.problem)
    print(f"forge exit={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
