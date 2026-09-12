"""① 風洞 (モード F) の評価: 逆設計 → メッシュ → forge → 軸 M 達成度 (Phase 3)。

v1 チェーン: dv (Bézier 自由 CP + 軸長 L_ax + 円弧 R) → Sauer starting line →
逆 MOC マーチ → 壁 (Euler 用は非粘性壁 / NS 用は +δ* 経験式) → ModeFWall →
TFI メッシュ → forge (Euler=cell/slip/visc0 | RANS-SST=node) → 軸 M vs 目標。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_wt \
      case/41.wind_tunnel_design/problem_m4.yaml run_dir [--euler] [--steps N]
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.bezier import MachBezier
from ..geometry.cminus_map import build_cminus_map
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import SauerThroat
from ..geometry.wall_modef import ModeFWall
from ..feedback.deltastar import deltastar_offset
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge


def design_chain(p: Problem, viscous: bool) -> dict:
    r"""dv → 目標 Bézier → 逆設計壁 (+δ*) → ModeFWall。決定的。

    **B7 (plan §9.2)**: 目標曲線は $x_{\mathrm{reach}}$ (設計壁始点 $x_0$ から出て
    最初に成功する C⁻ の軸着地点 — 壁の帰還がそもそも届かない境界) で分割する。
    $x_0\le x<x_{\mathrm{reach}}$ は実測 (CFD アンカー時は `axis_segment`、Sauer
    時は解析式) をそのまま目標にし、自由 CP は $x\ge x_{\mathrm{reach}}$ だけを
    担当する。Bézier のグローバル台 (Bernstein 基底は局所支持を持たない) が遠方
    CP の形を到達不能域まで押し戻す問題を解消する (実測: 分割前は
    $x/r_t\approx1.38$ で CP2 が既に 11% の重みを持ち $\Delta M\approx+0.05$ の
    主残差を作っていた)。$x_{\mathrm{reach}}$ 自体は下流形状にほぼ依存しない
    ため、**まず旧来の単一 Bézier ("ブートストラップ") で仮設計して MOC マップの
    被覆から $x_{\mathrm{reach}}$ を得てから、本目標を組んで MOC をやり直す**。
    """
    g = p.gamma
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))   # 物理定数でなく数値上の選択
    # **B6 CFD アンカー** (plan §9.2): `geometry.throat_anchor_run` があれば
    # Sauer をその run の CFD 抽出値で置き換える (アンカーは抽出元 run に凍結)
    anchor_run = p.geometry.get("throat_anchor_run")
    if anchor_run:
        from ..feedback.cfd_anchor import from_run
        st = from_run(anchor_run, float(p.spec["r_throat"]), R, g)
    else:
        st = SauerThroat(R=R, gamma=g)
    n_start = int(p.geometry.get("n_start", 41))
    x0, rr0, MM0, tt0 = st.starting_line(M_start=M_start, n=n_start)
    h = 1e-4
    M0 = float(st.mach(x0, 0.0))
    dM0 = float((st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h))
    d2M0 = float((st.mach(x0 + h, 0.0) - 2.0 * M0 + st.mach(x0 - h, 0.0)) / h ** 2)
    # **B1: ブートストラップ段の中心線接続次数のみに効く** (plan §9.2)。B7 導入後は
    # 本目標側の x_reach 接続が常に C2 (M,M',M'') を課すため、ここは x_reach の
    # 推定に使う仮設計の質にしか影響しない。
    order = int(p.geometry.get("centerline_start_order", 1))
    if order not in (1, 2):
        raise ValueError("geometry.centerline_start_order は 1 か 2")
    start_bc = (M0, dM0) if order == 1 else (M0, dM0, d2M0)
    xd = dv_value(p, "L_ax")
    cps = [dv_value(p, f"mc_cp{i+1}") for i in range(3)]
    eps = (1.0 / Md) * ((2.0 + (g - 1.0) * Md * Md) / (g + 1.0)) ** (
        (g + 1.0) / (2.0 * (g - 1.0)))
    x_end = xd + 2.3 * np.sqrt(eps) / np.tan(np.arcsin(1.0 / Md))
    n_axis_inv = int(p.geometry.get("n_axis_inv", 500))
    th_wall0 = float(np.arcsin(min(x0 / R, 1.0)))

    # === B10 (2026-08-17): **target と自由 dv の起点を $x_{reach,CFD}$ にする** ===
    # `geometry.x_reach_cfd` が与えられていればその方式。無ければ従来 B8。
    #
    # 思想: 自由壁が軸へ影響できる最初の位置は $J$ (固定円弧終端) 発の C⁻ が軸へ
    # 着地する $x_{reach,CFD}$ であり、$x_0$ ではない。そこより上流に target を
    # 置くと「直しようのないズレ」を設計に組み込むことになる。B7 は上流に実測こぶ
    # 曲線をコピーする帳簿で回避しようとしたが、それは帳簿操作であって設計ではない。
    # B10 は**上流に target を置かない** — 上流の軸 M・こぶ・圧力勾配は target 誤差
    # ではなく**独立の物理診断**として記録する (plan §9.2 B10)。
    #
    # $x_{reach,CFD}$ と境界値 $(M,M',M'')$ は正式ゲートを通した **node** 基準 run
    # から一度取得して固定する (cell は離散化雑音で $M''$ が測定不能 — 実測で
    # stencil 間 6.09、node は 0.049。plan §9.2 B10-a/b)。
    x_reach_cfd = p.geometry.get("x_reach_cfd")
    if x_reach_cfd is not None:
        x_reach = float(x_reach_cfd)
        anc = p.geometry.get("x_reach_anchor")     # [M, M', M''] (node 基準 run 由来)
        if anc is None or len(anc) != 3:
            raise ValueError("geometry.x_reach_cfd 指定時は geometry.x_reach_anchor: "
                             "[M, M', M''] が必須 (node 基準 run から取得した値)")
        start_r = (float(anc[0]), float(anc[1]), float(anc[2]))
        # 自由 CP は**拘束付き最小二乗で決定的に射影**する (手調整・旧値の
        # インデックスずらし流用は禁止 — plan §9.2 B10)。参照曲線は現行下流 target。
        bz_ref = MachBezier.from_constraints(x0, xd, start=start_bc, free_cp=cps,
                                             end=(Md, 0.0, 0.0))
        bz = MachBezier.fit_free_cp(x_reach, xd, start_r, len(cps), (Md, 0.0, 0.0),
                                    lambda x: float(np.atleast_1d(bz_ref(float(x)))[0]))

        def target(x: float) -> float:      # $x \ge x_{reach}$ にのみ定義
            x = float(x)
            if x < x_reach:
                raise ValueError(f"target は x >= x_reach ({x_reach:.4f}) にのみ定義 "
                                 f"(要求 x={x:.4f})。上流は物理診断であって target 誤差ではない")
            return Md if x >= xd else float(np.atleast_1d(bz(x))[0])

        # 逆設計へは $[x_0, x_{reach})$ の軸データも要る (MOC の Cauchy 入力)。
        # ここは**設計 target ではなく**、基準 run の実測軸 M をそのまま渡す
        # (壁は $J$ 上流で固定なので実測が「達成される値」そのもの)。
        seg = st.axis_segment_curve(x0, x_reach) if hasattr(st, "axis_segment_curve") else None

        def target_moc(x: float) -> float:
            x = float(x)
            if x < x_reach:
                return (float(seg(np.float64(x))) if seg is not None
                        else float(st.mach(x, 0.0)))
            return target(x)

        target_eval = target        # 評価も $x \ge x_{reach}$ のみ (帳簿コピー廃止)
        res = inverse_design(st, target_moc, x_axis_end=float(x_end), n_axis=n_axis_inv,
                             n_start=n_start, gamma=g, th_wall0=th_wall0, M_start=M_start)
    else:
        # === [旧] B8 (2026-08-16): 設計 = B6 グローバル Bézier / 評価 = B7 分割帳簿 ===
        bz = MachBezier.from_constraints(x0, xd, start=start_bc, free_cp=cps,
                                         end=(Md, 0.0, 0.0))

        def target(x: float) -> float:      # 設計用 = B6 グローバル (滑らか)
            return Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])

        res = inverse_design(st, target, x_axis_end=float(x_end), n_axis=n_axis_inv,
                             n_start=n_start, gamma=g, th_wall0=th_wall0, M_start=M_start)
        cmap = build_cminus_map(res["pts"], res["wall"], g)
        if len(cmap) < 4:
            raise ValueError("C⁻ マップ点が不足 (x_reach 未確定)")
        x_reach = float(cmap[:, 1].min())
        if hasattr(st, "axis_segment_curve"):
            seg = st.axis_segment_curve(x0, x_reach)

            def target_eval(x: float) -> float:
                x = float(x)
                if x < x_reach:
                    return float(seg(np.float64(x)))
                return target(x)
        else:
            target_eval = target

    wall_inv = res["wall"]
    pts = (deltastar_offset(wall_inv, float(p.spec["r_throat"]),
                            float(p.spec["Pt"]), float(p.spec["Tt"]), g, p.cp)
           if viscous else wall_inv[:, :2])
    wall = ModeFWall(pts, R=R,
                     r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                     L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                     L_contract=float(p.geometry.get("L_contract", 3.0)),
                     phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))),
                     blend_len=float(p.geometry.get("blend_len", 0.0)))
    msgs = wall.validate()
    if msgs:
        raise ValueError("モード F 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": wall_inv, "target": target,
            "target_eval": target_eval, "x0": float(x0),
            "xd": float(xd), "Md": Md, "pts": res["pts"], "R": R,
            "centerline_start_order": order, "sauer_d2M0": d2M0,
            "x_reach": x_reach, "x_reach_mode": ("cfd" if x_reach_cfd is not None else "moc"),
            "free_cp": [float(v) for v in bz.cp[len(bz.cp) - 3 - len(cps):len(bz.cp) - 3]],
            "target_d2M0": float(np.atleast_1d(bz.deriv(bz.x0, 2))[0]),
            "mdot_ratio_moc": res["mdot_exit"] / res["mdot_start"]}


def _config_euler(p: Problem, nsteps: int, out_int: int, cfl: float, conv: int) -> str:
    return f"""mesh: {{meshFormat: "hdf5", discretization: "cell", isAxisymmetric: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0,
           ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {conv}, limiter: 2}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""


def _config_euler_node(p: Problem, nsteps: int, out_int: int, cfl: float,
                       conv: int = 1) -> str:
    """**node (median-dual) Euler** — B10 の基準 run 用 (2026-08-17)。

    動機: cell モードは atomicAdd 由来のセルレベル雑音床 ~6e-4 を持ち
    ([[cell-atomicadd-nondeterminism]])、軸 $M$ の 2 階微分がこの雑音に埋もれて
    $M''$ を測定できない (実測: 雑音 $\\epsilon/dx^2$ が細分で 0.097→1.53 と増え、
    測定値 −1.33/+0.51/−5.59 と同オーダー)。node は雑音床 ~1e-6 で約 600 倍良い。

    設定の注意:
    - `nodeWallDirichlet` は**書かない**。既定は 1 だが、これが作用する `wall_flag`
      は `bcondKind` が `wall`/`wall_isothermal` の CV にしか立たない
      (`mesh/mesh.cpp`) ため、**slip 壁では no-op** で害はない (2026-08-17 に
      ソース確認。「slip で接線速度まで殺す」と書いた旧記述は誤りで訂正)。
      no-slip へ切り替えるときに初めて意味を持つ
      ([[node-wall-velocity-needs-dirichlet-flag]])。
    - 軸対称 node の旧 `nodeAxisDirichlet` (軸行真空化の対症, 2026-08-16 撤去) と
    [2026-08-16 更新: `nodeAxisDirichlet` は撤去。node 既定が「値=ノード座標 + 軸 u_r 三点セット + エッジ中点再構成 + LSQ」
     になり (case/43)、軸ノードは通常 DOF として解く。]
      `axisCentroidShift` が要る ([[node-axisym-axis-dirichlet]])。
    - `bndFirstOrder` は**使用禁止** (AGENTS)。
    - メッシュは**この config で変換**すること — cell 用 h5 を流用すると黙って
      primal cell で解かれる ([[node-mesh-must-convert-with-node-config]])。
    - 既知の未解決事項: node slip 壁 + 接線密度勾配の市松スプリアス流
      ([[node-slip-spurious-flow]])。本ケースは slip 壁なので、同じ市松診断で
      cell と比較して確認すること。
    """
    return f"""mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0,
           ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {conv}, limiter: 2}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""


def _config_sst_node(p: Problem, nsteps: int, out_int: int, cfl: float) -> str:
    """NS v1 (δ* 段) の node 壁解像 y+~1 SST 設定。

    case/40 `run_0050_node_yp1_outletic` (node y+1 の生産構成、「正 (最終形)」判定)
    の設定を踏襲: node + `nodeWallDirichlet` (壁速度 Dirichlet — 既定 0 だと壁速度が
    非零 [[node-wall-velocity-needs-dirichlet-flag]]) +
    `axisCentroidShift` + 低 Re SST (`wallTreatmentSST: 0` — **壁関数は使わない**。
    node 壁関数経路には既知の Cf 欠損 [[node-wallfunction-pk-convention-deficit]]、
    ユーザ指示 2026-08-16)。粘性は Sutherland (viscMethod 1)。"""
    return f"""mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 1,
           ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2}}
turbulence: {{model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, turbulentPrandtl: 0.9}}
initial: "uniform_p101325_u10"
"""


def _bcond(p: Problem, euler: bool) -> str:
    Pt, Tt = float(p.spec["Pt"]), float(p.spec["Tt"])
    pa = float(p.spec.get("p_ambient", 1000.0))
    # 壁行は Problem.wall_bcond_line (spec.wall_thermal が単一ソース: 断熱 wall / 等温 wall_isothermal+Ts)
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{Pt: {Pt}, Tt: {Tt}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}
wall:   {p.wall_bcond_line(euler, phys_id=3, output=1)}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


def prepare(problem_path, run_dir, euler: bool = True, nsteps=None,
            ic_from=None) -> dict:
    """run ディレクトリ準備。euler=False で **NS v1** (δ* オフセット壁 +
    node 壁解像 y+~1 SST — `_config_sst_node`)。

    ic_from: warm start 元の run dir (最終 res を `interp_field.py` で cross-mesh
    移植。AGENTS のメッシュ変更後 restart 必須ルール)。Euler 場を SST へ移植する
    場合、Euler res の k/omega は 0 なので **interp 後に roK/roOmega を入口相当値
    で貼り直す** (SST init 0 は step 0 で omega 発散 [[sst-mesh-walldist-gotcha]])。
    """
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym":
        raise ValueError("runner_wt は wind_tunnel_axisym 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p, viscous=not euler)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 目標軸分布と設計壁の記録 (帰還 v2 とΔM 評価の入力)
    xs = np.linspace(d["x0"], d["xd"], 400)
    np.savetxt(run_dir / "target_axis_M.csv",
               np.c_[xs * scale, [d["target"](x) for x in xs]], delimiter=",",
               header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
    # solverConfig は convertGmshToForge が cwd から読む — node/SST の場合は
    # **変換前に**書いておくことで node (median-dual) メッシュ + no-slip 壁の
    # wall_dist が作られる ([[node-mesh-must-convert-with-node-config]])
    cfg = (_config_euler(p, n, out_int, 4.0, 1) if euler
           else _config_sst_node(p, n, out_int, 4.0))
    (run_dir / "bcondConfig.yaml").write_text(_bcond(p, euler))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    if not euler:
        # 品質ツールは node 形式の CONNE を読めない (cell 前提)。品質は primal
        # メッシュの性質なので **cell 変換の一時コピーで検査**してから node 変換。
        (run_dir / "solverConfig.yaml").write_text(_config_euler(p, n, out_int, 4.0, 1))
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                       cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
        q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                            "nozzle_qc.h5"], cwd=run_dir, env=_ENV,
                           capture_output=True, text=True)
        (run_dir / "MESH_QUALITY.txt").write_text("# cell 変換コピーで検査 (品質は primal の性質)\n"
                                                  + q.stdout + q.stderr)
        (run_dir / "nozzle_qc.h5").unlink()
        if q.returncode != 0:
            raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    (run_dir / "solverConfig.yaml").write_text(cfg)
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    if euler:
        q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                            "nozzle.h5"], cwd=run_dir, env=_ENV, capture_output=True, text=True)
        (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
        if q.returncode != 0:
            raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        # ソースが Euler (k/omega≈0) のときだけ乱流を入口相当で貼り直す。
        # ソースが SST 中継場なら**実 k/ω 構造を保持** (これを潰すと y+~1 の
        # 剛直壁で 2 次化直後に出口コーナー ω NaN — run_0028 で実測)
        import h5py as _h5
        with _h5.File(src) as fs:
            om_src_max = float(np.max(fs["/VALUE/omega"][:])) if "/VALUE/omega" in fs else 0.0
        if om_src_max < 1.0:
            with _h5.File(run_dir / "nozzle.h5", "r+") as f:
                ro = f["/VALUE/ro"][:]
                f["/VALUE/roK"][:] = ro * 1.0
                f["/VALUE/roOmega"][:] = ro * 18000.0
        if not euler:
            # **近壁 ω の粘性底層フロア** (2026-08-16, run_0030 初回で実測):
            # 粗→細 (y+50→y+1.5) の最近傍 interp は細メッシュ近壁の ω を
            # 粗第一セル値 (~1e7) のまま貼るが、平衡値は ~1e12 — 5 桁不足の
            # 調整過渡で本段 (2次+cfl4) が壁帯 ω NaN。標準の底層解
            # ω = 6ν/(β1 y²) (β1=0.075) を床として適用する (壁遠方は床が
            # 微小になり無影響)。
            from ..feedback.deltastar import _sutherland
            g = p.gamma
            R_gas = p.cp * (g - 1.0) / g
            with _h5.File(run_dir / "nozzle.h5", "r+") as f:
                ro = f["/VALUE/ro"][:]
                roUx, roUy = f["/VALUE/roUx"][:], f["/VALUE/roUy"][:]
                roe = f["/VALUE/roe"][:]
                wd = f["/VALUE/wall_dist"][:]
                T = np.maximum((roe - 0.5 * (roUx ** 2 + roUy ** 2) / np.maximum(ro, 1e-12))
                               * (g - 1.0) / (np.maximum(ro, 1e-12) * R_gas), 50.0)
                nu = _sutherland(T) / np.maximum(ro, 1e-12)
                om_floor = 6.0 * nu / (0.075 * np.maximum(wd, 1e-9) ** 2)
                f["/VALUE/roOmega"][:] = np.maximum(f["/VALUE/roOmega"][:], ro * om_floor)
    info = {"x0": d["x0"], "xd": d["xd"], "Md": d["Md"], "x_reach": d["x_reach"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "scale_m": scale,
            "euler": euler, "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(problem_path, run_dir, nsteps=None) -> None:
    """soft 段 (1次+cfl0.5, 3000 step) → 本段 — bell と同じ 2 段起動を config
    書換えで同一 run dir 内実行 (soft 完走後に本 config へ戻し restart なしで
    冷間から回すのでなく、soft の最終場を IC に採用する)。"""
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    cfg_soft = cfg_main.replace("cfl: 4.0, cfl_pseudo: 4.0", "cfl: 0.5, cfl_pseudo: 0.5")
    cfg_soft = cfg_soft.replace("convMethod: 1", "convMethod: 0")
    # node y+~1 SST の確定レシピ ([[node-yp1-dual-geometry-float32-fix]]):
    # soft 段は nStepInner 10 (壁 ω ~1e10 の剛性対策。cfl4 直投入は step ~69 で
    # roOmega NaN — run_0028 初回で再確認)
    cfg_soft = cfg_soft.replace("nStepInner: 5", "nStepInner: 10")
    import re
    cfg_soft = re.sub(r"nStepOuter: \d+", "nStepOuter: 3000", cfg_soft)
    cfg_soft = re.sub(r"outStepInterval: \d+", "outStepInterval: 3000", cfg_soft)
    (run_dir / "solverConfig.yaml").write_text(cfg_soft)
    run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not res or int("".join(c for c in res[-1].stem if c.isdigit())) < 3000:
        raise RuntimeError("soft 段が失敗")
    # soft 最終場を nozzle.h5 の IC に移植 (同一メッシュ → interp は index 一致)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                    str(res[-1]), str(run_dir / "nozzle.h5")],
                   env=_ENV, check=True, capture_output=True, text=True)
    for f in run_dir.glob("res_*"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    """達成軸 M vs 目標の ΔM (固定格子)・質量流量比を評価。"""
    from ..metrics.extract import axis_mach

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    tgt = np.loadtxt(run_dir / "target_axis_M.csv", delimiter=",", skiprows=1)
    x_ach, M_ach = axis_mach(run_dir / "nozzle.h5", res[-1],
                             axis_band=0.9 * float(p.spec["r_throat"]) * 0.1)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    scale = info["scale_m"]
    # 評価窓: M>1.05 の少し先 〜 目標終端手前 (§4.7(c) のマスク方針)
    xs = np.linspace(info["x0"] * scale + 0.5 * scale, info["xd"] * scale - 0.3 * scale, 200)
    Mt = np.interp(xs, tgt[:, 0], tgt[:, 1])
    Ma = np.interp(xs, x_ach, M_ach)
    dM = Ma - Mt
    out = {"res_file": res[-1].name,
           "dM_max": float(np.max(np.abs(dM))),
           "dM_max_rel_Md": float(np.max(np.abs(dM)) / info["Md"]),
           "dM_rms": float(np.sqrt(np.mean(dM ** 2))),
           "M_axis_exit": float(np.interp(info["xd"] * scale, x_ach, M_ach)),
           "mdot_ratio_moc": info["mdot_ratio_moc"]}
    np.savetxt(run_dir / "achieved_vs_target.csv", np.c_[xs, Mt, Ma],
               delimiter=",", header="x_m,M_target,M_achieved", comments="")
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="モード F 風洞評価 (Phase 3 v1)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, euler=True, nsteps=a.steps)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    run_staged(a.problem, a.run_dir, a.steps)
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
