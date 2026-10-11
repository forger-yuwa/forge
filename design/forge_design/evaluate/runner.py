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

TP 経路 (2026-10-11, plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2a): 問題 YAML が `gas.model: semiperfect` +
`gas.exhaust` (CEA の排気ファイル; `gas.rocket_exhaust`) のときだけ、ガスを燃焼室の平衡組成を凍結した TP (単一の lump
擬似種 EXH、`thermalMethod: 2`・`thermoHrefTemp`・`viscMethod: 2` + `physProp.transport`) にし、IC・壁の平滑化・
1D 理想値を同じ契約 (`RocketExhaust`) の気体で作り、真空の出口を `outflow` にする。config は推奨レシピ
(procedures/recommended-settings.md §1・§1.2・§2・§3) の値で書き、段 (`stage`: soft | mid | main) を選べる。
既定の CPG 経路は従来の config・BC・IC・平滑化のまま (この節の追加で文字単位で変わらない)。
TP の段階起動・同一メッシュの継続・推力の CV 会計 (2026-10-11, 同 plan §5.1 #2b): `run_staged_tp` (soft → mid → 本段を同じ run・
同じメッシュで、段の間は restart_field.py、stage_manifest.json)、`continue_tp` (本段の継続 = polish)、`collect` が本段の全スナップショットで
`metrics.extract.thrust_cv` を出して thrust_history.csv を書く。CLI は `--staged` / `--continue-from <run>`。
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


def converter_path() -> Path:
    """convertGmshToForge の場所: 環境変数 FORGE_CONVERTER > FORGE_BIN (run_case.sh が使う forge) と同じビルドの変換器 >
    solver_density_cuda/build。変換器も solverConfig.yaml を読むので、forge と同じビルドのものを使う
    (lump 記法の physProp.species は plan thermophysics-solver-owned-species-db #6a 以降の変換器でないと読めない)。"""
    cand = os.environ.get("FORGE_CONVERTER")
    if cand:
        return Path(cand)
    fb = os.environ.get("FORGE_BIN")
    if fb and (Path(fb).parent / "convertGmshToForge").exists():
        return Path(fb).parent / "convertGmshToForge"
    return FORGE_BUILD / "convertGmshToForge"


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
    return f"""mesh: {{discretization: "{disc}"{node_keys}, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1,
           visc: 1.8e-5, thermCond: 0.0257, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
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
    if _is_tp(p):
        return _bcond_config_tp(p)
    Pt = float(p.spec["Pt"])
    Tt = float(p.spec["Tt"])
    pa = float(p.spec["p_ambient"])
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{Pt: {Pt}, Tt: {Tt}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}
wall:   {{physID: 3, kind: wall,             outputHDFflg: 1, ints: , floats: }}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


# ---------------------------------------------------------------- TP 経路 (plan tooling-nozzle-bell-two-mode §4.1, §5.1 #2a)

# 段階起動の段 (procedures/recommended-settings.md §1.2: soft → mid → 本段)。implicitRelax 0.7 は §1・§3 の基本設定で全段に付ける。
# main の cfl は evaluate.cfl_main (既定 4.0; §1 の node NS ノズル 4〜6、§3 の TP + relax 0.7 は 6〜8 まで可)
TP_STAGES = {
    "soft": {"convMethod": 0, "limiter": 0, "cfl": 0.5, "nStepInner": 10},
    "mid": {"convMethod": 0, "limiter": 0, "cfl": 1.0, "nStepInner": 10},
    "main": {"convMethod": 1, "limiter": 2, "cfl": None, "nStepInner": 4},
}
TP_NODE_KEYS = ", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, nodeInletCornerWall: 1"
TP_K_INIT, TP_OMEGA_INIT = 1.0, 18000.0     # IC・入口の k/ω (CPG 経路と同じ値)
# 推力の CV 会計 (metrics.extract.thrust_cv) が入口・出口面の粘性応力に使うソルバの勾配 (output.extraFields; 場は不変)
TP_EXTRA_FIELDS = ("dUxdx", "dUxdy", "dUydx", "dUydy", "axisym_divU")
# 段階起動の既定の step 数 (recommended-settings §1.2: soft 2000 → mid 2000 → 本段 8000〜12000)。evaluate.stage_steps で上書き。
# 本段は evaluate.nStepOuter、res の間隔は evaluate.outStepInterval (既定は本段の 1/12 — 推力の時系列の準定常判定に 4 点以上)
TP_STAGE_STEPS = {"soft": 2000, "mid": 2000}


def _is_tp(p: Problem) -> bool:
    """TP 経路か (`gas.model: semiperfect` で `evaluate.cfd_gas` が cpg でない)。③ベルの TP は組成・全温・全圧を
    `gas.exhaust` (CEA の排気ファイル) で与える (plan bell-two-mode §4.1)。CPG (既定) は False。"""
    if not p.uses_tp_cfd:
        return False
    if p.exhaust is None:
        raise ValueError("③ベルの TP は gas.exhaust (CEA の排気ファイル; python -m forge_design.gas.rocket_exhaust で作る) で "
                         "組成・燃焼室温度・燃焼室圧を与える (plan tooling-nozzle-bell-two-mode §4.1)")
    return True


def _outlet_kind(p: Problem) -> str:
    """`evaluate.outlet_kind`: outflow | outlet_statPress。既定は TP = outflow (真空の設計点、出口は超音速の流出;
    plan §4.1)、CPG = outlet_statPress (従来)。"""
    k = str(p.evaluate.get("outlet_kind", "outflow" if _is_tp(p) else "outlet_statPress"))
    if k not in ("outflow", "outlet_statPress"):
        raise ValueError(f"evaluate.outlet_kind '{k}' は未知 (outflow | outlet_statPress)")
    return k


def tp_gas_record(p: Problem) -> dict:
    """TP の気体の来歴と CFD の参照値 (prepare_info.json・config 用)。cp・gamma は燃焼室温度の凍結値
    (TP では NASA-9 が優先され、R = cp(γ−1)/γ が混合気の気体定数に一致するように置く参照値)。
    visc・thermCond は CEA の燃焼室の凍結輸送物性 (visc は局所 dt と陰解法対角の剛性見積り: 領域の最大 ≈ 燃焼室)。"""
    ex = p.exhaust
    gas = ex.frozen_gas()
    cp_c = float(gas.cp_mass(ex.Tc)[0])
    ch = ex.data["chamber"]
    return {"exhaust": ex.summary(), "R": float(gas.R), "cp_ref": cp_c, "gamma_ref": cp_c / (cp_c - float(gas.R)),
            "visc_ref": float(ex.visc_estimate() or 1.0e-4),
            "thermCond_ref": float(ch.get("cond_W_mK_cea_frozen") or 0.5), "thermoHrefTemp": float(ex.href_T)}


def _solver_config_tp(p: Problem, nsteps: int, out_interval: int, stage: str = "main",
                      discretization: str | None = None) -> str:
    """TP (凍結組成の lump 擬似種 EXH) の solverConfig。推奨レシピ: §1 (node・SLAU・陰解法 + implicitRelax 0.7・
    output level 1)、§1.2 (段)、§2 (SST、`wallTreatmentSST: 0`)、§3 (`thermalMethod: 2`・`thermoHrefTemp`)、
    `viscMethod: 2` + `physProp.transport` (solver-settings)。`mesh.bndFirstOrder` は書かない。
    `discretization="cell"` は品質検査用の一時変換だけに使う (node 専用キーを落とす)。"""
    from ..gas.composition import physprop_species_flow, physprop_transport_flow, solver_species_config
    if stage not in TP_STAGES:
        raise ValueError(f"stage '{stage}' は未知 ({list(TP_STAGES)})")
    disc = discretization or p.mesh.get("discretization", "node")
    if disc != "node" and discretization is None:
        raise ValueError("③ベルの TP は node (median-dual) のみ (recommended-settings §1)")
    if int(p.evaluate.get("wall_treatment_sst", 0)) != 0:
        raise ValueError("③ベルの TP は低 Re 壁 (wallTreatmentSST: 0) のみ。壁関数は使わない (recommended-settings §2)")
    layout = p.species_layout()
    items, external = solver_species_config(layout)
    _db = ', speciesDBFile: "species_db_external.yaml"' if external else ""
    transport = p.transport_for_ns(layout)
    rec = tp_gas_record(p)
    sst = p.evaluate.get("turbulence", "sst") == "sst"
    st = TP_STAGES[stage]
    cfl = float(p.evaluate.get("cfl_main", 4.0)) if st["cfl"] is None else st["cfl"]
    node_keys = TP_NODE_KEYS if disc == "node" else ", isAxisymmetric: 1"
    turb = ('turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, '
            f'turbulentPrandtl: 0.9, kInit: {TP_K_INIT}, omegaInit: {TP_OMEGA_INIT}}}' if sst else 'turbulence: {model: "none"}')
    return f"""mesh: {{discretization: "{disc}"{node_keys}, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 2, viscMethod: 2, visc: {rec['visc_ref']!r}, thermCond: {rec['thermCond_ref']!r}, prandtlLam: 0.72,
           cp: {rec['cp_ref']!r}, gamma: {rec['gamma_ref']!r},
           species: {physprop_species_flow(items)}{_db}, thermoHrefTemp: {rec['thermoHrefTemp']!r},
           transport: {physprop_transport_flow(transport)}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1.0e-5, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: 0.7,
           dt_min: 1.0e-8, dt_max: 1.0, blockDPLUR: 1, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_interval}
  timeIntegration: 11
  nStepInner: {st['nStepInner']}
space: {{convMethod: {st['convMethod']}, limiter: {st['limiter']}}}
{turb}
output: {{level: 1, extraFields: [{", ".join(TP_EXTRA_FIELDS)}]}}
initial: "uniform_p101325_u10"
"""


def _bcond_config_tp(p: Problem) -> str:
    """TP の bcondConfig。入口 = 燃焼室のよどみ (Pt = Pc, Tt = Tc; 排気ファイル)、壁 = spec.wall_thermal (既定 断熱 wall)、
    出口 = `_outlet_kind` (既定 outflow)。単一輸送種 (EXH) なので入口の Y は書かない (多種なら Y0.. を足す)。"""
    Pt, Tt = float(p.spec["Pt"]), float(p.spec["Tt"])
    layout = p.species_layout()
    ys = "" if layout.n == 1 else "".join(f"Y{i}: {y!r}, " for i, y in enumerate(layout.Y_transport("inflow")))
    if _outlet_kind(p) == "outflow":
        outlet = "{physID: 2, kind: outflow,          outputHDFflg: 1, ints: , floats: }"
    else:
        pa = float(p.spec["p_ambient"])
        outlet = f"{{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}"
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{{ys}Pt: {Pt!r}, Tt: {Tt!r}, k: {TP_K_INIT}, omega: {TP_OMEGA_INIT}}}}}
outlet: {outlet}
wall:   {p.wall_bcond_line(False, phys_id=3, output=1)}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


def ideal_1d_tp(p: Problem) -> dict:
    """TP の 1D 理想値 (同じ凍結 TP で ε まで等エントロピー膨張): ṁ_1D・F_1D,vac(ε)・F_1D(ε; p_a) [物理単位、A_t = π r_t²]。"""
    rt = float(p.spec["r_throat"])
    return p.exhaust.ideal_1d(float(p.spec["eps"]), p_a=float(p.spec["p_ambient"]), A_t=float(np.pi * rt * rt))


PROBE_STUB = "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n"

# 第一セル厚の x 依存 (スロートだけ詰める) と細分幅のキー (meshing.mesh2d.Mesh2DParams)。**書いてあるときだけ**渡す
# (書いていない問題は従来の 4 キーだけの格子とビット同一; plan bell-two-mode §5.1 #2b 引き継ぎ (c): スロートの y1+ ≤ 1)。
MESH_OPTIONAL_KEYS = ("throat_width", "wall_first_frac_throat", "wall_first_blend_x0", "wall_first_blend_x1",
                      "wall_first_up_x0", "wall_first_up_x1")


def mesh_params(p: Problem, scale: float) -> Mesh2DParams:
    """問題の mesh ブロック → Mesh2DParams。ni・nj・wall_first_frac・throat_refine は従来の既定、`MESH_OPTIONAL_KEYS` は
    書いてあるときだけ渡す (以前は黙って無視していた)。"""
    kw = {k: float(p.mesh[k]) for k in MESH_OPTIONAL_KEYS if p.mesh.get(k) is not None}
    return Mesh2DParams(
        ni=int(p.mesh.get("ni", 241)),
        nj=int(p.mesh.get("nj", 81)),
        wall_first_frac=float(p.mesh.get("wall_first_frac", 2.0e-3)),
        throat_refine=float(p.mesh.get("throat_refine", 3.0)),
        scale=scale,
        **kw,
    )


def prepare(problem_path, run_dir, nsteps=None, warm_from=None, stage=None) -> dict:
    """メッシュ・config・IC まで準備する (forge は起動しない)。

    warm_from: 過去の収束 res_*.h5 (cell/node 不問)。指定時は準 1D IC の上から
    `interp_field.py` (最近傍 cross-mesh restart) で場を移植する。**node は必須**
    (冷間 IC は発散 — 実証レシピ参照)。
    stage: TP 経路だけ (soft | mid | main; 既定は warm_from なしで soft、ありで main)。TP の冷間 IC は soft 段に限って
    warm_from なしで起動できる (recommended-settings §1.2)。CPG 経路では指定不可 (従来の 1 段 config)。
    注意: warm_from は形状が変わる cross-mesh の引き継ぎ (`interp_field.py`) で、**同一メッシュの段の継続には使わない**
    (AGENTS「メッシュ変更後の restart」: 同一メッシュは `restart_field.py`)。TP の段の継続 (`restart_field.py` + `StageManifest`)
    は `run_staged_tp`・`continue_tp` (plan bell-two-mode §5.1 #2b)。
    """
    p = load_problem(problem_path)
    tp = _is_tp(p)
    if stage is not None and not tp:
        raise ValueError("stage は TP 経路 (gas.model: semiperfect + gas.exhaust) だけで指定できる")
    if tp:
        stage = _tp_stage(stage, warm_from)     # run dir を作る前に段を検査する
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    wall, wp = build_wall(p)
    scale = float(p.spec["r_throat"])

    mp = mesh_params(p, scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 壁座標も記録 (IC・後処理・帰還エンジンの入力)
    xs = np.linspace(wall.x_in, wall.x_e, 800)
    np.savetxt(run_dir / "contour.csv", np.c_[xs * scale, wall.r(xs) * scale],
               delimiter=",", header="x_m,r_m", comments="")
    if tp:
        return _prepare_tp(p, problem_path, run_dir, wall, scale, mp, quads, nsteps, warm_from, stage)

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


def _tp_stage(stage, warm_from) -> str:
    """TP の段を決めて検査する: 既定は warm_from なしで soft、ありで main。冷間 IC (準 1D 等エントロピー) からは soft だけ
    (recommended-settings §1.2: 一様・準 1D の IC から 2 次・本段 CFL を直接始めない)。"""
    stage = stage or ("main" if warm_from is not None else "soft")
    if stage not in TP_STAGES:
        raise ValueError(f"stage '{stage}' は未知 ({list(TP_STAGES)})")
    if warm_from is None and stage != "soft":
        raise ValueError("TP の冷間 IC (準 1D 等エントロピー) からは soft 段で起動する (recommended-settings §1.2)。"
                         "mid / main は前段の res_*.h5 を warm_from に与える")
    return stage


def _environment_record() -> dict:
    import platform

    import scipy
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "forge_bin": os.environ.get("FORGE_BIN"), "converter": str(converter_path())}


def _prepare_tp(p: Problem, problem_path, run_dir: Path, wall, scale: float, mp, quads, nsteps, warm_from, stage) -> dict:
    """TP 経路の prepare (メッシュ生成の後; plan bell-two-mode §4.1 #2a)。config は `_solver_config_tp`、BC は
    `_bcond_config_tp`、IC は排気の契約の気体 (`p.gas_model` = `RocketExhaust.moc_gas`) で `thermoHrefTemp` の基準の
    内部エネルギー、新規 IC の化学種属性は宛先の `forge --resolve-species` と照合して付ける。変換器は `converter_path()`
    (FORGE_BIN と同じビルド; lump 記法の physProp.species を読める変換器が要る)。"""
    from ..gas.composition import species_db_raw_yaml, solver_species_config, write_species_meta
    from .ic import stamp_isentropic_ic_species

    stage = _tp_stage(stage, warm_from)
    ex = p.exhaust
    layout = p.species_layout()
    transport = p.transport_for_ns(layout)
    _, external = solver_species_config(layout)
    write_species_meta(layout, run_dir, transport)
    if external:
        (run_dir / "species_db_external.yaml").write_text(species_db_raw_yaml(external))
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
    if n % out_int:
        # forge は outStepInterval の倍数でしか res を書かない → 最終 step の res が残るよう n を割り切る間隔に直す
        out_int = next(n // k for k in (3, 2, 4, 6, 1) if n % k == 0)
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config_tp(p))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    conv = str(converter_path())
    # 品質ゲートは primal (cell 変換) で測る (CPG 経路と同じ理由)。AR の上限は mesh.ar_max (既定 1000; 壁法線の構造層は 5000 まで可)
    ar_max = int(p.mesh.get("ar_max", 1000))
    (run_dir / "solverConfig.yaml").write_text(_solver_config_tp(p, n, out_int, stage, discretization="cell"))
    subprocess.run([conv, "nozzle.msh", "nozzle_qc.h5"], cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "nozzle_qc.h5", "--ar-max", str(ar_max)],
                       cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text("(primal/cell 変換で検査)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    for f in run_dir.glob("nozzle_qc.xmf"):
        f.unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    (run_dir / "solverConfig.yaml").write_text(_solver_config_tp(p, n, out_int, stage))
    subprocess.run([conv, "nozzle.msh", "nozzle.h5"], cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)

    rec = tp_gas_record(p)
    gas = p.gas_model                       # = ex.moc_gas() (同じ Y・種 DB・Tc)
    species = list(layout.species)
    Ys = None if layout.n == 1 else layout.Y_transport("inflow")
    exit1d = paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale, float(p.spec["Pt"]), float(p.spec["Tt"]),
                                 rec["gamma_ref"], rec["cp_ref"], k_init=TP_K_INIT, omega_init=TP_OMEGA_INIT,
                                 gas=gas, h_ref_T=ex.href_T, species_Y=Ys)
    # 新規初期場の化学種属性: IC と同じ気体・基準・輸送種順序・MW を宛先の --resolve-species の記録と照合 (違えば例外)
    stamp = stamp_isentropic_ic_species(run_dir / "nozzle.h5", run_dir, gas, ex.href_T, species,
                                        [float(layout.entries[s].MW) for s in species], Ys)
    if warm_from is not None:
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(warm_from), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        _smooth_wall_entropy(run_dir / "nozzle.h5", rec["gamma_ref"], rec["cp_ref"], gas=ex.frozen_gas())
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n, "outStepInterval": out_int,
            "stage": stage, "warm_from": (str(warm_from) if warm_from is not None else None),
            "scale_m": scale, "exit_1d": exit1d, "ic_species": stamp,
            "gas": {**rec, "species": species, "transport": transport, "mode": layout.mode},
            "ideal_1d": ideal_1d_tp(p), "outlet_kind": _outlet_kind(p),
            "mesh": {"ni": mp.ni, "nj": mp.nj, "cells": int(quads.shape[0]), "wall_first_frac": mp.wall_first_frac,
                     "ar_max": ar_max, **{k: getattr(mp, k) for k in MESH_OPTIONAL_KEYS if p.mesh.get(k) is not None}},
            "environment": _environment_record()}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    return info


def _tp_T_of_e(gas, e) -> np.ndarray:
    """TP の比内部エネルギー e (forge の thermoHrefTemp 基準) → T。`gas` は FrozenGas (e_sens は単調増加)。"""
    Tg = np.linspace(50.0, 5999.0, 60000)
    eg = gas.e_sens(Tg)
    return np.interp(np.asarray(e, dtype=float), eg, Tg)


def _smooth_wall_entropy(h5path, gamma: float, cp: float, npass: int = 10, gas=None) -> None:
    """cell→node 最近傍 interp が壁ノード列に刻む等圧エントロピー市松 (T/ρ 逆相ジグザグ) を除去する。

    壁ノードは u=0 で移流平滑化が働かず、P は音響で均されるが T/ρ の市松は接線熱伝導でしか
    消えない (実測: 静的凍結、ベル部 |ΔΔT|~500K@interp直後 → 収束後も ~200K 残存)。
    ここで壁列の T を x ソートで 1-2-1 フィルタし、P 不変で ro/roe を再構成する
    (u は壁 Dirichlet が 0 にするため運動量も 0 化)。plan
    boundary-node-nozzle-wall-outlet-stability §2.8 参照。
    `gas` (TP 経路; FrozenGas): T を TP の内部エネルギー e(T) から引き、roe を ρ e(T) で組み直す (CPG の
    P/(γ−1) は使わない)。None (CPG) は従来の式のまま。
    """
    import h5py

    R = cp * (gamma - 1.0) / gamma if gas is None else float(gas.R)
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
        if gas is None:
            P = (gamma - 1.0) * (roe[w] - ke)
            T = P / np.maximum(ro[w] * R, 1e-30)
        else:
            T = _tp_T_of_e(gas, (roe[w] - ke) / np.maximum(ro[w], 1e-30))
            P = ro[w] * R * T
        for _ in range(npass):
            T[1:-1] = 0.25 * T[:-2] + 0.5 * T[1:-1] + 0.25 * T[2:]
        ro_new = P / np.maximum(R * T, 1e-30)
        if gas is None:
            roe_new = P / (gamma - 1.0)  # u=0 (壁 Dirichlet と整合、運動量も 0 化)
        else:
            roe_new = ro_new * gas.e_sens(T)  # u=0 (同上)、TP の内部エネルギー (thermoHrefTemp 基準)
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


# ---------------------------------------------------------------- TP の段階起動と同一メッシュの継続 (plan bell-two-mode §4.1, §5.1 #2b)

def _res_sorted(run_dir) -> list:
    """run の主スナップショット res_<step>.h5 を step 順に (境界ダンプ res_<群>_<id>_<step>.h5 は含めない)。"""
    import re
    rs = [f for f in Path(run_dir).glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)]
    return sorted(rs, key=lambda f: int(f.stem.split("_")[1]))


def tp_stage_plan(p: Problem, stages=None) -> list:
    """TP の段の計画 [(tag, 段, step 数, outStepInterval)]。段は `stages` > `evaluate.stages` > (soft, mid, main)。
    soft・mid の step 数は `evaluate.stage_steps` (既定 TP_STAGE_STEPS)、本段は `evaluate.nStepOuter` (既定 12000) で
    res を `evaluate.outStepInterval` (既定は本段の 1/12) ごとに書く。前段は最終 step だけを書く。最後の段は main、
    本段の step 数が outStepInterval の倍数でなければ ValueError (最終 step の res が書かれず判定がずれるため; 黙って直さない)。"""
    stages = tuple(stages or p.evaluate.get("stages") or ("soft", "mid", "main"))
    bad = [s for s in stages if s not in TP_STAGES]
    if bad or not stages or stages[-1] != "main" or "main" in stages[:-1]:
        raise ValueError(f"TP の段 {stages}: 段は {list(TP_STAGES)} から選び、main で終える (main は最後の 1 回だけ)")
    steps = {**TP_STAGE_STEPS, **(p.evaluate.get("stage_steps") or {})}
    n_main = int(p.evaluate.get("nStepOuter", 12000))
    out_main = int(p.evaluate.get("outStepInterval", max(n_main // 12, 1)))
    if n_main <= 0 or out_main <= 0 or n_main % out_main:
        raise ValueError(f"本段の nStepOuter {n_main} が outStepInterval {out_main} の正の倍数でない")
    plan = []
    for i, s in enumerate(stages):
        n = n_main if s == "main" else int(steps[s])
        if n <= 0:
            raise ValueError(f"段 {s} の step 数 {n} が正でない")
        plan.append((f"S{i + 1}_{s}", s, n, out_main if s == "main" else n))
    return plan


def _tp_stage_run(p: Problem, run_dir: Path, sm, bc_text: str, tag: str, stage: str, n: int, out_int: int) -> Path:
    """1 段を回す: config を書き、forge を起動し、段の残差履歴を `residual_history_<tag>.csv` に写して manifest に足す。
    最終 res が step n まで届かない・forge が非 0・段終了ゲート (必要保存量が有限で ρ>0; runner_axismach.stage_gate) に
    落ちたら RuntimeError (次段へ進まない)。戻り値は段の最終 res。"""
    import shutil

    from .runner_axismach import stage_gate
    cfg = _solver_config_tp(p, n, out_int, stage)
    (run_dir / "solverConfig.yaml").write_text(cfg)
    (run_dir / "residual_history.csv").unlink(missing_ok=True)   # 前段の履歴を別段の名前で写さない
    rc = run_forge(run_dir)
    hist = run_dir / "residual_history.csv"
    if hist.exists():
        shutil.copy(hist, run_dir / f"residual_history_{tag}.csv")
    sm.add(tag, cfg, bc_text, history=f"residual_history_{tag}.csv")
    sm.write()
    res = _res_sorted(run_dir)
    last = int(res[-1].stem.split("_")[1]) if res else -1
    if rc != 0 or last < n:
        raise RuntimeError(f"TP の段 {tag} が失敗 (rc={rc}, 最終 res step {last} < {n}; 発散の切り分けは res_nan_*.h5)")
    probs = stage_gate(res[-1], cfg)
    if probs:
        raise RuntimeError(f"TP の段 {tag} の最終場 {res[-1].name} が段終了ゲートで不合格 — 次段へ進まない:\n    "
                           + "\n    ".join(probs))
    return res[-1]


def run_staged_tp(problem_path, run_dir, warm_from=None, stages=None) -> dict:
    """TP の段階起動 (recommended-settings §1.2): `prepare` (冷間は準 1D 等エントロピー IC、`warm_from` は形状が変わる
    cross-mesh の `interp_field.py`) の後、**同じ run・同じメッシュ**で段 (`tp_stage_plan`) を順に回す。段の間は
    `restart_field.py` (保存量の index コピー、ビット一致の検査; 引き継ぎ (a))。段ごとの実効設定を `stage_manifest.json`
    (StageManifest) に記録するので、`check_convergence.py --segment` が本段の区間 (方程式・BC・離散化が同じ段の連結) を
    判定する。前段の res は消し、本段の res と境界ダンプを残す (推力の時系列用)。戻り値は prepare_info (段の計画を追記)。"""
    from .runner_axismach import _restart_same_mesh, _stage_manifest_cls
    p = load_problem(problem_path)
    if not _is_tp(p):
        raise ValueError("run_staged_tp は TP 経路 (gas.model: semiperfect + gas.exhaust) 専用")
    plan = tp_stage_plan(p, stages)
    if warm_from is None and plan[0][1] != "soft":
        raise ValueError("TP の冷間 IC は soft 段で起動する (recommended-settings §1.2)")
    info = prepare(problem_path, run_dir, nsteps=plan[0][2], warm_from=warm_from, stage=plan[0][1])
    run_dir = Path(run_dir)
    bc_text = (run_dir / "bcondConfig.yaml").read_text()
    sm = _stage_manifest_cls()(run_dir)
    for k, (tag, st, n, oi) in enumerate(plan):
        last = _tp_stage_run(p, run_dir, sm, bc_text, tag, st, n, oi)
        if k == len(plan) - 1:
            break
        _restart_same_mesh(last, run_dir / "nozzle.h5")          # 同一メッシュ: index コピー
        for f in run_dir.glob("res_*"):
            f.unlink()
    info["stages"] = [{"tag": t, "stage": s, "nStepOuter": n, "outStepInterval": oi} for t, s, n, oi in plan]
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    return info


# 同一メッシュの継続で親 run から写すもの (メッシュと入力・化学種の記録・段の記録)。res・ログ・VERDICT は写さない
TP_CONTINUE_COPY = ("nozzle.h5", "nozzle.xmf", "nozzle.msh", "contour.csv", "bcondConfig.yaml", "probe.yaml",
                    "species_meta.yaml", "species_db_external.yaml", "MESH_QUALITY.txt", "forge_launches.jsonl")


def continue_tp(problem_path, src_run, new_run, nsteps=None) -> dict:
    """TP の本段を**同一メッシュで継続**する (polish 用; plan §4.1「polish は restart_field.py による同一メッシュの継続」)。
    親 run のメッシュ・入力・化学種の記録・段の履歴と manifest を新しい run に写し、親の最終 res を `restart_field.py` で
    初期場にして本段をもう 1 段 (`S<k>_main_p`) 回す。親の本段と同じ実効設定なので `check_convergence.py --segment` は
    親の本段とこの段を 1 区間として判定する。`nsteps` 既定は親の本段と同じ。戻り値は新しい run の prepare_info。"""
    import shutil

    from .runner_axismach import _restart_same_mesh, _stage_manifest_cls
    p = load_problem(problem_path)
    if not _is_tp(p):
        raise ValueError("continue_tp は TP 経路専用 (CPG の polish は従来どおり)")
    src_run, new_run = Path(src_run), Path(new_run)
    res = _res_sorted(src_run)
    if not res:
        raise FileNotFoundError(f"{src_run}: 継続元の res_*.h5 が無い")
    man_p = src_run / "stage_manifest.json"
    if not man_p.exists():
        raise FileNotFoundError(f"{src_run}: stage_manifest.json が無い (段の記録の無い run は継続しない)")
    man = json.loads(man_p.read_text())
    if not man["stages"] or not man["stages"][-1]["tag"].split("_", 1)[-1].startswith("main"):
        raise ValueError(f"{src_run}: 最後の段が本段でない ({[s['tag'] for s in man['stages']]})")
    info_src = json.loads((src_run / "prepare_info.json").read_text())
    st_main = [s for s in info_src.get("stages", []) if s["stage"] == "main"]
    if not st_main:
        raise ValueError(f"{src_run}: prepare_info に本段の記録が無い")
    n = int(nsteps or st_main[-1]["nStepOuter"])
    oi = int(st_main[-1]["outStepInterval"])
    if n % oi:
        raise ValueError(f"継続の step 数 {n} が outStepInterval {oi} の倍数でない")
    new_run.mkdir(parents=True, exist_ok=False)
    for name in TP_CONTINUE_COPY:
        if (src_run / name).exists():
            shutil.copy2(src_run / name, new_run / name)
    for f in src_run.glob("resolved_species_*.yaml"):
        shutil.copy2(f, new_run / f.name)
    for st in man["stages"]:
        h = src_run / st["history"]
        if h.exists():
            shutil.copy2(h, new_run / st["history"])
    sm = _stage_manifest_cls()(new_run)
    sm.stages = [dict(s) for s in man["stages"]]
    k = len(sm.stages) + 1
    tag = f"S{k}_main_p"
    bc_text = (new_run / "bcondConfig.yaml").read_text()
    (new_run / "solverConfig.yaml").write_text(_solver_config_tp(p, n, oi, "main"))   # restart_field が宛先を解決する
    _restart_same_mesh(res[-1], new_run / "nozzle.h5")
    info = {**info_src, "run_dir": str(new_run), "continued_from": str(res[-1].resolve()),
            "stages": info_src.get("stages", []) + [{"tag": tag, "stage": "main", "nStepOuter": n, "outStepInterval": oi}]}
    (new_run / "prepare_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    _tp_stage_run(p, new_run, sm, bc_text, tag, "main", n, oi)
    return info


def tp_iso_k(run_dir) -> bool:
    """run の solverConfig の `turbulence.sstIsotropicStress` (−(2/3)ρk δij を応力に足すか; 既定 0)。重複キーは拒否。"""
    from .runner_axismach import _cfg_load
    doc = _cfg_load((Path(run_dir) / "solverConfig.yaml").read_text())
    t = doc.get("turbulence") or {}
    return int((t.get("sstIsotropicStress", 0) if isinstance(t, dict) else 0) or 0) == 1


THRUST_HISTORY_COLS = ("thrust_exit_N", "thrust_wall_N", "eta", "eta_wall", "mdot_exit_kg_s", "mdot_in_kg_s",
                       "mdot_ratio_1d", "cv_residual", "friction_N")
# 準定常ゲート (ii) で判定する列 (両推力会計・η・ṁ・CV 残差; plan §4.1)
THRUST_QS_COLS = ("thrust_exit_N", "thrust_wall_N", "eta", "mdot_exit_kg_s", "cv_residual")


def thrust_history_tp(p: Problem, run_dir) -> list:
    """本段の全スナップショット (res_<step>.h5 と同じ step の壁ダンプがあるもの) で推力の CV 会計を出す。
    壁ダンプの無い step (forge は step 0 の境界ダンプを書かない) は摩擦が評価できないので時系列に入れない。"""
    from ..metrics.extract import CV_IDS, thrust_cv
    run_dir = Path(run_dir)
    ideal = ideal_1d_tp(p)
    gas = p.exhaust.frozen_gas()
    iso_k = tp_iso_k(run_dir)
    rows = []
    for r in _res_sorted(run_dir):
        step = int(r.stem.split("_")[1])
        wd = run_dir / f"res_wall_{CV_IDS['wall']}_{step}.h5"     # bcondConfig の壁 (名前 wall・physID 3)
        if not wd.exists():
            continue
        o = thrust_cv(run_dir / "nozzle.h5", r, float(p.spec["p_ambient"]), float(ideal["F"]), wall_dump=wd, gas=gas,
                      iso_k=iso_k, mdot_ideal=float(ideal["mdot"]))
        o["step"] = step
        o["eta"] = o["eta_exit"]
        rows.append(o)
    return rows


def write_thrust_history(path, rows) -> Path:
    """thrust_history_tp の結果を `step` + THRUST_HISTORY_COLS の CSV に (check_quasisteady --series-csv 用)。"""
    import csv
    path = Path(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("step",) + THRUST_HISTORY_COLS)
        for o in rows:
            w.writerow([o["step"]] + [("" if o.get(c) is None else repr(float(o[c]))) for c in THRUST_HISTORY_COLS])
    return path


def collect(problem_path, run_dir) -> dict:
    """収束 VERDICT + 最終 res からメトリクスを抽出し metrics.json に保存。"""
    from ..metrics.extract import thrust_metrics

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text() if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else "(no verdict)"
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    out = {"convergence_verdict": verdict.strip().splitlines()[-2:] if verdict else []}
    if res and _is_tp(p):
        # TP: 推力は閉じた CV の 2 会計 (metrics.extract.thrust_cv; CPG の圧力復元は使わない)。本段の全スナップショットで
        # 時系列 thrust_history.csv を書き (準定常ゲート用)、最終 step の値を metrics に入れる。η は出口面会計 (eta_cf)
        res = _res_sorted(run_dir)
        rows = thrust_history_tp(p, run_dir)
        write_thrust_history(run_dir / "thrust_history.csv", rows)
        out["res_file"] = res[-1].name
        out["ideal_1d"] = ideal_1d_tp(p)
        if not rows or rows[-1]["step"] != int(res[-1].stem.split("_")[1]):
            out["thrust"] = f"最終 step {res[-1].name} の壁ダンプが無く推力の CV 会計を出せない"
        else:
            last = rows[-1]
            out.update({k: last[k] for k in ("thrust_exit_N", "thrust_wall_N", "friction_N", "eta_exit", "eta_wall",
                                             "cv_residual", "cv_closable", "mdot_exit_kg_s", "mdot_in_kg_s",
                                             "mass_imbalance", "mdot_ratio_1d", "pressure_source", "open_viscous",
                                             "wall_traction_convention", "exit_area_m2")})
            out["eta_cf"] = last["eta_exit"]
            out["thrust_N"] = last["thrust_exit_N"]
            out["mdot_kg_s"] = last["mdot_exit_kg_s"]
            out["cv_groups"] = last["groups"]
            out["thrust_history"] = {"file": "thrust_history.csv", "n": len(rows),
                                     "steps": [rows[0]["step"], rows[-1]["step"]]}
    elif res:
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
    ap.add_argument("--stage", default=None, choices=sorted(TP_STAGES),
                    help="TP 経路の段 (既定: warm-from なし soft / あり main)。CPG では指定不可")
    ap.add_argument("--staged", action="store_true",
                    help="TP 経路: 段階起動 (evaluate.stages、既定 soft → mid → main) を同じ run・同じメッシュで回す "
                         "(段の間は restart_field.py、stage_manifest.json に記録)")
    ap.add_argument("--continue-from", default=None,
                    help="TP 経路: 指定 run の本段を同一メッシュで継続する (restart_field.py; polish と同じ)")
    a = ap.parse_args(argv)
    if a.staged or a.continue_from:
        if (a.staged and a.continue_from) or a.prepare_only or a.stage or (a.steps and not a.continue_from) \
                or (a.warm_from and a.continue_from):
            ap.error("--staged と --continue-from は排他。--prepare-only・--stage とは併用しない "
                     "(--steps は --continue-from の継続 step 数だけ、--warm-from は --staged の cross-mesh 起点だけ)")
        if a.continue_from:
            info = continue_tp(a.problem, a.continue_from, a.run_dir, nsteps=a.steps)
        else:
            info = run_staged_tp(a.problem, a.run_dir, warm_from=a.warm_from)
        out = collect(a.problem, a.run_dir)
        print(json.dumps({k: v for k, v in out.items() if k != "cv_groups"}, indent=1, ensure_ascii=False))
        return 0
    info = prepare(a.problem, a.run_dir, a.steps, warm_from=a.warm_from, stage=a.stage)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    rc = run_forge(a.run_dir)
    out = collect(a.problem, a.run_dir)
    print(json.dumps(out, indent=1))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
