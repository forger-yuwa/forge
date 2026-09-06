"""⑤ SERN 評価 runner (S3): 逆設計 → 2 バンド構造メッシュ → forge 平面 2D (Euler/SST) → 力係数。

plan: plans/active/tooling-nozzle-sern-chain.md §4.2, §4.7。問題型 `sern_2d`。
1 run = 1 作動点。作動点は spec.external で与える (多作動点束ねは S6 で driver 側)。
段階起動: soft (1 次 + cfl 0.5, 3000 step) → 本段 (2 次 + cfl_main)。IC は領域別一様
(中間線より上 = 燃焼器出口状態、下 = 外部流)。
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

from ..geometry.moc_sern import PlanarMOC, SernKernelSpec, wall_forces
from ..geometry.rao_planar import ideal_gross_thrust
from ..meshing.mesh_sern import PHYS_SERN, SernMeshParams, generate_sern_mesh, write_msh41_named
from ..metrics.sern_forces import force_history, steadiness
from ..probdef import Problem, dv_value, load_problem

# リポジトリ位置から導く (AWS など別マシンでも動くように。FORGE_ROOT で上書き可)
FORGE_ROOT = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
FORGE_TOOLS = FORGE_ROOT / "solver_density_cuda" / "tools"
FORGE_BUILD = FORGE_ROOT / "solver_density_cuda" / "build"
_ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
MESH = "sern.h5"


def _dv(p: Problem, name, default=None) -> float:
    v = dv_value(p, name, default)
    return float(v["value"] if isinstance(v, dict) else v)


def design_snapshot(p: Problem) -> dict:
    """作動点で上書きされる**前**の設計点 (入口・外部流・ガス) を控える。逆設計はこれで固定する
    (plan §4.10: 形状は設計点で 1 つに決まる。作動点は CFD の BC/IC だけを変える)。"""
    return {"inflow": dict(p.spec["inflow"]), "external": dict(p.spec["external"]),
            "gamma": float(p.gamma), "cp": float(p.cp)}


def select_operating_point(p: Problem, op: str | None) -> dict:
    """spec.operating_points[] から名前で 1 点を選び spec.external / inflow / ガスに反映する。
    op=None で operating_points が無ければ spec.external をそのまま使う。戻り値 = 選んだ点 (重み込み)。

    作動点は飛行条件 (external) だけでなく**燃焼器出口 (inflow) とガス (gamma/cp) も動く** —
    飛行 M が変わればインレット圧縮と燃焼加熱が変わり、powered と power-off では γ が 1.18 と 1.39 で違う
    (plan §4.10、NASA TM X-71972 TABLE 1 + CEA2)。"""
    ops = p.spec.get("operating_points")
    if not ops:
        if op not in (None, "", "default"):
            raise ValueError("spec.operating_points が無いのに --op が指定された")
        return {"name": "default", "weight": 1.0, "external": dict(p.spec["external"]),
                "inflow": dict(p.spec["inflow"]), "gas": {"gamma": p.gamma, "cp": p.cp}}
    names = [o["name"] for o in ops]
    if op in (None, "", "default"):
        op = names[0]
    if op not in names:
        raise ValueError(f"作動点 '{op}' が無い (候補: {names})")
    o = ops[names.index(op)]
    p.spec["external"] = dict(o["external"])
    if "inflow" in o:
        p.spec["inflow"] = {**p.spec["inflow"], **o["inflow"]}
    if "gas" in o:                      # 作動点ごとの γ / cp (powered 1.18 vs power-off 1.39)
        gas = {k: float(v) for k, v in o["gas"].items()}
        unknown = set(gas) - {"gamma", "cp"}
        if unknown:
            raise ValueError(f"operating_points[{op}].gas の未知キー: {sorted(unknown)} (gamma | cp のみ)")
        p.gamma = gas.get("gamma", p.gamma)
        p.cp = gas.get("cp", p.cp)
        p.raw.setdefault("gas", {}).update(gas)
    return {"name": op, "weight": float(o.get("weight", 1.0)), "external": dict(p.spec["external"]),
            "inflow": dict(p.spec["inflow"]), "gas": {"gamma": p.gamma, "cp": p.cp}}


def gas_states(p: Problem) -> dict:
    g, cp = p.gamma, p.cp
    R = cp * (g - 1.0) / g
    fi, ex = p.spec["inflow"], p.spec["external"]
    if fi.get("mode", "supersonic") != "supersonic":
        raise ValueError("inflow.mode は現状 supersonic のみ (sonic_throat は S1 接続が未実装)")

    def st(M, P, T):
        ro = P / (R * T)
        u = M * np.sqrt(g * R * T)
        k = 1.5 * (0.01 * u) ** 2
        omega = ro * k / (1.8e-5 * 10.0)
        return {"M": float(M), "P": float(P), "T": float(T), "ro": float(ro), "u": float(u),
                "k": float(k), "omega": float(omega)}
    return {"exhaust": st(fi["M_in"], fi["p_in"], fi["T_in"]),
            "ext": st(ex["M_inf"], ex["p_inf"], ex["T_inf"]), "R": R}


def design_from_problem(p: Problem, design: dict | None = None):
    """逆設計は**設計点**で行う (`design` = 作動点適用前の `design_snapshot(p)`)。作動点 (operating_points) は
    CFD の境界条件・IC だけを変え、形状は変えない (2026-09-05 修正: それまで作動点ごとに p_ext が変わり kernel/形状が
    作動点依存になっていた — run_0010/0017 は作動点間で形状が一致していない)。
    2026-09-05 追補 (§4.10): 作動点が inflow / gas も動かすようになったので、入口状態と γ も設計点で固定する。"""
    geo = p.geometry
    d0 = design or design_snapshot(p)
    fi_d, ext_d, g_d = d0["inflow"], d0["external"], float(d0["gamma"])
    p_ext_ratio = float(ext_d["p_inf"]) / float(fi_d["p_in"])
    spec = SernKernelSpec(M_in=float(fi_d["M_in"]),
                          theta_r0=np.deg2rad(_dv(p, "theta_r0_deg")), theta_c0=np.deg2rad(_dv(p, "theta_c0_deg")),
                          L_cowl=_dv(p, "L_cowl"), gamma=g_d, p_ext_over_p_in=p_ext_ratio,
                          x_max=float(geo.get("x_max_kernel", 10.0)), nj=int(geo.get("nj_moc", 301)),
                          dx=float(geo.get("dx_moc", 2e-3)))
    if str(geo.get("mode", "keypoint")) == "straight":
        k = PlanarMOC(spec).march()
        d = k.straight_design(_dv(p, "L_ramp"))
    else:
        k = PlanarMOC(spec).march(stop_at=(_dv(p, "f"), _dv(p, "M_c")))   # c の少し先で打ち切る
        d = k.design_ramp(M_c=_dv(p, "M_c"), f=_dv(p, "f"))
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    fr = wall_forces(d, spec.M_in, g_d, pa_over_pin=p_ext_ratio, x_ref=float(xr), y_ref=float(yr))
    theta_b = float(k.TH[-1, 0])   # 自由境界の終端角 (せん断層に格子線を沿わせる)
    return k, d, fr, theta_b


def _solver_config(p: Problem, nsteps: int, out_int: int, cfl: float, p_ref: float) -> str:
    """`evaluate.implicit_relax` / `evaluate.p_min` を deltaT / space に挿入できる (2026-09-06)。

    SERN は元から `blockDPLUR: 1` (陰解法) だが `implicitRelax` を設定しておらず既定 1.0 だった。
    風洞チェーンは同じ陰解法に対し「cfl 6 + implicitRelax 0.7 が生産推奨」(case/45 run_0018) と配管済みで、
    メモリ [[implicit-cfl-ceiling-eos-floor]] も「NS 陰解法の上限は P 床洗浄律速、効くのは implicitRelax のみ」
    と記録している。run_0075 の発散 (M∞10 の boat-tail 膨張で 18 % のノードが pMin=1 Pa に着地 → 負密度) は
    まさにこの指紋なので、relax を効かせる。"""
    _lim = int(p.evaluate.get("limiter", 2))   # 2=Venkatakrishnan (既定), 1=Barth
    ir = p.evaluate.get("implicit_relax")
    _relax = f", implicitRelax: {float(ir)}" if ir is not None else ""
    pm = p.evaluate.get("p_min")
    _pmin = f", pMin: {float(pm)}" if pm is not None else ""   # physProp 配下 (space ではない)
    disc = p.mesh.get("discretization", "cell")
    model = p.evaluate.get("model", "euler")
    node_keys = ", nodeWallDirichlet: 1" if (disc == "node" and model != "euler") else ""
    if model == "euler":
        phys = f"physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}{_pmin}}}"
        turb = 'turbulence: {model: "none"}'
    else:
        phys = (f"physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, "
                f"thermCondMethod: 1, prandtlLam: 0.72, cp: {p.cp}, gamma: {p.gamma}{_pmin}}}")
        turb = 'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 1}'
    return f"""mesh: {{meshFormat: "hdf5", discretization: "{disc}", isAxisymmetric: 0{node_keys}, meshFileName: "{MESH}", valueFileName: "{MESH}"}}
gpu: 1
solver: "SLAU"
{phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1{_relax}, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: {_lim}, pRef: {p_ref}}}
{turb}
initial: "uniform_p101325_u10"
"""


def _bcond_config(p: Problem, st: dict) -> str:
    model = p.evaluate.get("model", "euler")
    wall_kind = "slip" if model == "euler" else "wall"
    ex, en = st["exhaust"], st["ext"]

    def inlet(name, pid, s):
        return (f"{name}: {{physID: {pid}, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , "
                f"floats: {{ro: {s['ro']:.6g}, Ux: {s['u']:.6g}, Uy: 0.0, Uz: 0.0, Ps: {s['P']:.6g}, k: {s['k']:.6g}, omega: {s['omega']:.6g}}}}}\n")

    def outlet(name, pid):
        return (f"{name}: {{physID: {pid}, kind: outlet_statPress, outputHDFflg: 0, ints: , "
                f"floats: {{Ps: {en['P']:.6g}, Pt: {en['P']:.6g}, Tt: {en['T']:.6g}}}}}\n")

    def wall(name, pid):
        return f"{name}: {{physID: {pid}, kind: {wall_kind}, outputHDFflg: 1, ints: , floats: }}\n"
    P = PHYS_SERN
    return (inlet("inlet_nozzle", P["inlet_nozzle"], ex) + inlet("inlet_ext", P["inlet_ext"], en)
            + outlet("outlet", P["outlet"]) + wall("ramp", P["ramp"]) + wall("cowl_in", P["cowl_in"])
            + wall("cowl_out", P["cowl_out"]) + outlet("bottom", P["bottom"])
            + (outlet("top_out", P["top_out"]) if p.evaluate.get("top_out_kind", "outlet") == "outlet"
               else f"top_out: {{physID: {P['top_out']}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n")
            # 機体上面 + base (§4.11): 機体の力なので帳簿外だが base 圧の診断のため壁出力する。Euler/SST とも slip
            + (f"vehicle: {{physID: {P['vehicle']}, kind: slip, outputHDFflg: 1, ints: , floats: }}\n"
               if int(p.mesh.get("ext_top", 0)) else ""))


def apply_wall_offset(design, wall_offset: dict, H: float):
    """壁を法線方向 (流体と反対側) に dn(x) [m] だけ動かした SernDesign を返す (無次元化して適用)。
    ramp: 上壁なので +n = 上、cowl: 下壁 (内面が上向き) なので +n = 下。表は (x_m, dn_m) の 2 列。"""
    import copy
    d = copy.deepcopy(design)
    for name, arr in (("ramp", d.ramp_xy), ("cowl", d.cowl_xy)):
        tbl = wall_offset.get(name)
        if tbl is None:
            continue
        tbl = np.asarray(tbl, dtype=float)
        dn = np.interp(arr[:, 0], tbl[:, 0] / H, tbl[:, 1] / H, left=tbl[0, 1] / H, right=tbl[-1, 1] / H)
        t = np.gradient(arr, axis=0); t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-30)[:, None]
        nrm = np.column_stack([-t[:, 1], t[:, 0]]) if name == "ramp" else np.column_stack([t[:, 1], -t[:, 0]])
        arr += dn[:, None] * nrm
    return d


def paste_region_ic(h5path, y_mid, y_top, scale: float, st: dict, gamma: float) -> None:
    """領域別一様 IC: 中間線とランプ/プルーム上線の間 = 燃焼器出口状態、それ以外 (カウル下・ランプ側外部流) = 外部流。"""
    ex, en = st["exhaust"], st["ext"]
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn, yn = cc[:, 0] / scale, cc[:, 1] / scale
        upper = (yn > y_mid(xn)) & (yn < y_top(xn))
        ro = np.where(upper, ex["ro"], en["ro"]); u = np.where(upper, ex["u"], en["u"]); P = np.where(upper, ex["P"], en["P"])
        v = f["/VALUE"]
        v["ro"][:] = ro; v["roUx"][:] = ro * u; v["roUy"][:] = 0.0; v["roUz"][:] = 0.0
        v["roe"][:] = P / (gamma - 1.0) + 0.5 * ro * u * u
        if "roK" in v:
            v["roK"][:] = ro * np.where(upper, ex["k"], en["k"]); v["roOmega"][:] = ro * np.where(upper, ex["omega"], en["omega"])


def convert_mesh(run_dir, msh: str, out: str) -> None:
    """gmsh msh → forge h5。**exit code で判定しない**: 一部の環境 (AWS g5 / CUDA 13) で converter は
    h5 を書き切ってから終了時に `GPUassert: invalid argument` を出して非零で抜ける (既知・無害)。
    成否は出力ファイルの存在とサイズで見る。"""
    r = subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), msh, out], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    f = Path(run_dir) / out
    if not f.exists() or f.stat().st_size < 1024:
        raise RuntimeError(f"convertGmshToForge が {out} を作れなかった (rc={r.returncode})\n"
                           + (r.stdout or "")[-1500:] + (r.stderr or "")[-1500:])


def prepare(problem_path, run_dir, nsteps=None, op: str | None = None, wall_offset=None) -> dict:
    """op: 作動点名 (spec.operating_points)。wall_offset: {"ramp": (x_m, dn_m), "cowl": (x_m, dn_m)} の
    法線オフセット表 [m] (S5 δ* 一発補正。壁を流体と反対側へ dn だけ動かす)。"""
    p = load_problem(problem_path)
    if p.type != "sern_2d":
        raise ValueError("runner_sern は sern_2d 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d0 = design_snapshot(p)                        # 設計点 (作動点で上書きされる前に保存)
    opinfo = select_operating_point(p, op)
    st = gas_states(p)
    kern, design, fr_moc, theta_b = design_from_problem(p, design=d0)
    H = float(p.spec["H_m"])
    m = p.mesh
    mp = SernMeshParams(ni_up=int(m.get("ni_up", 16)), ni_noz=int(m.get("ni_noz", 120)), ni_plume=int(m.get("ni_plume", 220)),
                        nj_top=int(m.get("nj_top", 101)), nj_bot=int(m.get("nj_bot", 61)), L_up=float(m.get("L_up", 0.5)),
                        x_out_extra=float(m.get("x_out_extra", 2.0)), bot_depth=float(m.get("bot_depth", 3.0)),
                        first_wall_frac=float(m.get("first_wall_frac", 2e-3)),
                        cowl_thickness=float(m.get("cowl_thickness", 2e-3 if m.get("discretization", "cell") == "node" else 0.0)),
                        interface_angle=float(m.get("interface_angle_rad", theta_b)),
                        top_ext_angle=float(np.deg2rad(m.get("top_ext_angle_deg", np.rad2deg(design.info["theta_e"])))),
                        ext_top=bool(int(m.get("ext_top", 0))), top_depth=float(m.get("top_depth", 2.0)),
                        nj_ext_top=int(m.get("nj_ext_top", 41)), nj_wake=int(m.get("nj_wake", 9)),
                        vehicle_clearance=float(m.get("vehicle_clearance", 0.02)), first_top_frac=float(m.get("first_top_frac", 0.02)),
                        vehicle_taper=float(m.get("vehicle_taper", 0.0)), ramp_fillet=float(m.get("ramp_fillet", 0.0)),
                        scale=H)
    if wall_offset:
        design = apply_wall_offset(design, wall_offset, H)
    coords, quads, bedges, minfo, y_mid, y_top = generate_sern_mesh(design, mp)
    write_msh41_named(run_dir / "sern.msh", coords, quads, bedges, PHYS_SERN)
    np.savetxt(run_dir / "ramp_contour.csv", design.ramp_xy * H, delimiter=",", header="x_m,y_m", comments="")
    np.savetxt(run_dir / "cowl_contour.csv", design.cowl_xy * H, delimiter=",", header="x_m,y_m", comments="")
    n = int(nsteps or p.evaluate.get("nStepOuter", 6000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl = float(p.evaluate.get("cfl_main", 4.0))
    cfg = _solver_config(p, n, out_int, cfl, st["ext"]["P"])
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p, st))
    (run_dir / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    disc = p.mesh.get("discretization", "cell")
    # 品質ゲートは primal (cell) 変換で
    (run_dir / "solverConfig.yaml").write_text(cfg.replace(f'discretization: "{disc}"', 'discretization: "cell"').replace(", nodeWallDirichlet: 1", ""))
    convert_mesh(run_dir, "sern.msh", "sern_qc.h5")
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "sern_qc.h5", "--mode", "2d"], cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    if disc == "cell":
        (run_dir / "sern_qc.h5").rename(run_dir / MESH)
    else:
        (run_dir / "sern_qc.h5").unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg)
        convert_mesh(run_dir, "sern.msh", MESH)
    for f in run_dir.glob("sern_qc.xmf"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg)
    paste_region_ic(run_dir / MESH, y_mid, y_top, H, st, p.gamma)
    ex = st["exhaust"]
    F_ideal_nd, M_e_id = ideal_gross_thrust(ex["M"], st["ext"]["P"] / ex["P"], p.gamma)
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n, "H_m": H, "states": st,
            "operating_point": opinfo, "wall_offset": bool(wall_offset), "design_point": d0,
            "design": {"key_point": list(design.key_point), "foot_a": list(design.foot_a), "lip_e": list(design.lip_e),
                       "L_ramp": design.L_ramp, "mass_fraction_check": design.mass_fraction_check,
                       "theta_e_deg": float(np.rad2deg(design.info["theta_e"])), "theta_b_deg": float(np.rad2deg(theta_b)),
                       "p_te_over_p_in": kern.p_te_over_p_in, "warnings": design.info["warnings"]},
            "moc_forces": fr_moc, "F_ideal_N_per_m": F_ideal_nd * ex["P"] * H, "M_e_ideal": M_e_id,
            "mesh": minfo, "discretization": disc, "model": p.evaluate.get("model", "euler")}
    (run_dir / "solverConfig_main.yaml").write_text(cfg)
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def restart_by_index(res_h5, mesh_h5) -> None:
    """同一メッシュの stage 間移植: VALUE を index でコピーする (座標最近傍の `interp_field.py` は使わない)。
    理由 (2026-09-04, case/46 run_0009): スリットカウルの上下壁ノードは座標が一致し、最近傍補間が双子を同じ元
    ノードに写す → 排気側の壁ノードが外部流の圧力を持ち 2 次で発散した (interp_field の全 134 station で誤写像を確認)。"""
    with h5py.File(res_h5, "r") as src, h5py.File(mesh_h5, "r+") as dst:
        n = len(dst["VALUE/ro"])
        for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"):   # 状態量のみ (wall_dist は触らない)
            if k in src["VALUE"] and k in dst["VALUE"] and len(src["VALUE"][k]) == n:
                dst["VALUE"][k][:] = src["VALUE"][k][:]


def warm_from_run(dst_run_dir, src_run_dir) -> dict:
    """**別作動点**の収束場を同一メッシュの run へ移植する (plan §4.7 / §5.1-1c、codex レビューの A′)。

    保存量の素コピー (`restart_by_index`) は**作動点間では使えない**: 同じ `roe` を別の γ で読むと
    P = (γ−1)(roe − ½ρ|u|²) が (γ_dst−1)/(γ_src−1) 倍ずれる (m6_on→m10_on で 1.24 倍、→m4_off で 2.18 倍)。
    そこで入口状態の比で相似スケールし、**目標作動点の γ で roe を組み直す**:

        ρ' = ρ·sρ,  (ρu)' = (ρu)·sρ·su,  P' = P·sP,  roe' = P'/(γ_dst−1) + ½|ρu'|²/ρ'
        (ρk)' = (ρk)·sρ·su²   (k ~ u²),   (ρω)' = (ρω)·sρ·su   (ω ~ u/L, L は同一メッシュで不変)

    sρ, su, sP は入口 (燃焼器出口) 状態の比。NPR が違えば内部の波の当たり方は変わるので、
    移植後に**短い適応段**を必ず入れる (run_staged の warm_adapt_steps)。戻り値 = 使ったスケール。
    """
    dst_run_dir, src_run_dir = Path(dst_run_dir), Path(src_run_dir)
    si = json.loads((src_run_dir / "prepare_info.json").read_text())
    di = json.loads((dst_run_dir / "prepare_info.json").read_text())
    if si["mesh"]["cells"] != di["mesh"]["cells"] or si["design"]["L_ramp"] != di["design"]["L_ramp"]:
        raise ValueError("warm_from_run: メッシュが同一でない (同一設計・同一 dv の run 間でのみ使う)")
    se, de = si["states"]["exhaust"], di["states"]["exhaust"]
    g_s = float(si["operating_point"]["gas"]["gamma"]); g_d = float(di["operating_point"]["gas"]["gamma"])
    s_ro, s_u, s_P = de["ro"] / se["ro"], de["u"] / se["u"], de["P"] / se["P"]
    res = sorted(src_run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not res:
        raise RuntimeError(f"warm_from_run: {src_run_dir} に res_*.h5 が無い")
    with h5py.File(res[-1], "r") as src, h5py.File(dst_run_dir / MESH, "r+") as dst:
        n = len(dst["VALUE/ro"])
        if len(src["VALUE/ro"]) != n:
            raise ValueError("warm_from_run: VALUE 長が一致しない")
        ro = src["VALUE/ro"][:].astype(np.float64)
        mom = [src[f"VALUE/{k}"][:].astype(np.float64) for k in ("roUx", "roUy", "roUz")]
        roe = src["VALUE/roe"][:].astype(np.float64)
        P = (g_s - 1.0) * (roe - 0.5 * sum(m * m for m in mom) / np.maximum(ro, 1e-30))
        ro_n = ro * s_ro
        mom_n = [m * (s_ro * s_u) for m in mom]
        roe_n = P * s_P / (g_d - 1.0) + 0.5 * sum(m * m for m in mom_n) / np.maximum(ro_n, 1e-30)
        dst["VALUE/ro"][:] = ro_n
        for k, m in zip(("roUx", "roUy", "roUz"), mom_n):
            dst["VALUE"][k][:] = m
        dst["VALUE/roe"][:] = roe_n
        if "roK" in src["VALUE"] and "roK" in dst["VALUE"]:
            dst["VALUE/roK"][:] = src["VALUE/roK"][:] * (s_ro * s_u * s_u)
            dst["VALUE/roOmega"][:] = src["VALUE/roOmega"][:] * (s_ro * s_u)
    return {"src": str(src_run_dir), "s_ro": s_ro, "s_u": s_u, "s_P": s_P, "gamma": [g_s, g_d]}


def run_forge(run_dir) -> int:
    r = subprocess.run([str(FORGE_TOOLS / "run_case.sh"), str(Path(run_dir).resolve())], env=_ENV, capture_output=True, text=True)
    (Path(run_dir) / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def run_staged(run_dir, stages: str = "full", soft_steps: int = 3000, soft_cfl: float = 0.5, soft_conv: int = 0,
               warm_lam_steps: int = 0, warm_lam_cfl: float = 0.2, mid_steps: int = 0,
               warm_src=None, warm_adapt_steps: int = 500) -> int:
    """soft_cfl / soft_conv: soft 段の CFL と convMethod (既定 0.5 / 1 次)。3D SST の後縁 3 重点など、1 次でも
    立ち上がりが厳しいケースで下げる。

    warm_lam_steps > 0 (SST のみ): soft 段の**前に層流暖機段** (`turbulence: none` + 粘性あり、1 次、`warm_lam_cfl`) を
    入れる (plan §4.7)。カウル後縁のせん断層は排気と外部流の密度・温度比が大きく、平均場が立つ前に SST を回すと
    `roOmega` が発散する — 板厚では作動点ごとに要求が逆転して解けなかった (m4_off は 2e-3、m10_on は 5e-3 が必要)。
    暖機で平均場を作ってから乱流方程式を入れると 3 作動点とも通り、力係数は暖機なしの成功例と一致する (case/46 run_0048)。

    mid_steps > 0: soft と本段の間に **2 次 + soft CFL** の段を入れる。次数と CFL を同時に上げるとランプ膨張角部で
    `roOmega` が発散する (case/46 run_0054 doe_000, θ_r0 18.3°)。`procedures/solver-settings.md` の段階戦略に沿う。

    warm_src: **別作動点の収束 run ディレクトリ**。与えると層流暖機と soft を飛ばし、`warm_from_run` で
    熱力学整合リマップ → **適応段** (1 次, soft CFL, warm_adapt_steps) → mid → 本段 とする。
    NPR が大きく違う作動点間 (m6_on 35.4 ↔ m4_off 2.8) では使わないこと (§5.1-1c)。"""
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig_main.yaml").read_text() if (run_dir / "solverConfig_main.yaml").exists() \
        else (run_dir / "solverConfig.yaml").read_text()
    if stages == "main":  # soft 段の res_*.h5 が残っている状態から本段だけ回す
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg_main)
        return run_forge(run_dir)
    if stages == "none":
        return run_forge(run_dir)
    if warm_src is not None:      # 別作動点からの warm start: 暖機と soft を飛ばし、適応段だけ入れる
        info_warm = warm_from_run(run_dir, warm_src)
        (run_dir / "WARM_START.json").write_text(json.dumps(info_warm, indent=1))
        ad = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {soft_cfl}, cfl_pseudo: {soft_cfl}", cfg_main)
        ad = ad.replace("convMethod: 1", f"convMethod: {soft_conv}")
        ad = re.sub(r"nStepOuter: \d+", f"nStepOuter: {warm_adapt_steps}", ad)
        ad = re.sub(r"outStepInterval: \d+", f"outStepInterval: {warm_adapt_steps}", ad)
        (run_dir / "solverConfig.yaml").write_text(ad)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError("warm start の適応段が失敗 (res_nan_*.h5 / forge_run.log を見る)")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
        warm_lam_steps = 0      # 以降は mid → 本段
    if warm_lam_steps > 0 and 'model: "sst"' in cfg_main:      # 層流暖機段 (SST を後から入れる)
        lam = re.sub(r'turbulence: \{model: "sst"[^}]*\}', 'turbulence: {model: "none"}', cfg_main)
        if 'model: "none"' not in lam:
            raise RuntimeError("層流暖機: turbulence 行の置換に失敗 (solverConfig の書式が変わった)")
        lam = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {warm_lam_cfl}, cfl_pseudo: {warm_lam_cfl}", lam)
        lam = lam.replace("convMethod: 1", "convMethod: 0")
        lam = re.sub(r"nStepOuter: \d+", f"nStepOuter: {warm_lam_steps}", lam)
        lam = re.sub(r"outStepInterval: \d+", f"outStepInterval: {warm_lam_steps}", lam)
        (run_dir / "solverConfig.yaml").write_text(lam)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError("層流暖機段が失敗 (res_nan_*.h5 / forge_run.log を見る)")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
    if warm_src is not None:      # soft は適応段で代替済み → mid へ
        return _run_mid_and_main(run_dir, cfg_main, soft_cfl, mid_steps)
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {soft_cfl}, cfl_pseudo: {soft_cfl}", cfg_main)
    soft = soft.replace("convMethod: 1", f"convMethod: {soft_conv}")
    soft = re.sub(r"nStepOuter: \d+", f"nStepOuter: {soft_steps}", soft)
    soft = re.sub(r"outStepInterval: \d+", f"outStepInterval: {soft_steps}", soft)
    (run_dir / "solverConfig.yaml").write_text(soft)
    rc = run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if rc != 0 or not res:
        raise RuntimeError("soft 段が失敗 (res_nan_*.h5 / forge_run.log を見る)")
    restart_by_index(res[-1], run_dir / MESH)
    for f in run_dir.glob("res_*"):
        f.unlink()
    return _run_mid_and_main(run_dir, cfg_main, soft_cfl, mid_steps)


def _run_mid_and_main(run_dir, cfg_main: str, soft_cfl: float, mid_steps: int) -> int:
    """mid 段 (2 次 + soft CFL) → 本段。soft/暖機/warm start の後段として共用する。"""
    if mid_steps > 0:      # mid 段: 2 次に上げるが CFL は soft のまま (次数と CFL を同時に上げない)
        mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {soft_cfl}, cfl_pseudo: {soft_cfl}", cfg_main)
        mid = re.sub(r"nStepOuter: \d+", f"nStepOuter: {mid_steps}", mid)
        mid = re.sub(r"outStepInterval: \d+", f"outStepInterval: {mid_steps}", mid)
        (run_dir / "solverConfig.yaml").write_text(mid)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError("mid 段 (2 次 + soft CFL) が失敗 (res_nan_*.h5 / forge_run.log を見る)")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    st = info["states"]; ex, en = st["exhaust"], st["ext"]; H = info["H_m"]
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    hist = force_history(run_dir, p_a=en["P"], F_ideal=info["F_ideal_N_per_m"], H=H, x_ref=float(xr) * H, y_ref=float(yr) * H,
                         mdot_u_in=ex["ro"] * ex["u"] ** 2 * H, p_in=ex["P"],
                         twall_on_fluid=(info.get("discretization", "cell") == "cell"))
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-2:] if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else []
    out = {"convergence_verdict": verdict, "n_snapshots": len(hist), "history": hist,
           "operating_point": info.get("operating_point"), "L_ramp": info["design"]["L_ramp"]}
    if hist:
        last = hist[-1]
        out.update({k: last[k] for k in ("step", "C_T", "C_T_wall", "C_L", "C_M", "T_wall", "L", "M_noseup")})
        for k in ("C_T_with_shear", "C_T_friction", "sep_frac_ramp", "sep_x_min_ramp"):
            if k in last:
                out[k] = last[k]
        out["steadiness"] = {k: steadiness([h[k] for h in hist]) for k in ("C_T", "C_L", "C_M")}
        out["moc_forces"] = info["moc_forces"]
        # MOC は**設計点**の値なので、作動点が設計点と一致するときだけ差を出す (2026-09-05, plan §4.10:
        # 作動点は inflow/gas も動かすので、オフデザイン run で差を取ると意味の無い数になる)。
        d0, o = info.get("design_point"), info.get("operating_point") or {}
        on_design = bool(d0) and all(o.get(k) == d0.get(k) for k in ("inflow", "external")) \
            and (o.get("gas", {}).get("gamma") == d0.get("gamma"))
        out["on_design_point"] = on_design
        if on_design:
            out["cfd_vs_moc"] = {k: (last[k] - info["moc_forces"][k]) for k in ("C_T", "C_L", "C_M")}
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="forge_design ⑤ SERN 評価 (1 作動点)")
    ap.add_argument("problem"); ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--stages", default="full", choices=["full", "none", "main"])
    ap.add_argument("--op", default=None, help="作動点名 (spec.operating_points)")
    ap.add_argument("--wall-offset", default=None, help="δ* オフセット JSON ({ramp: [[x_m, dn_m],...], cowl: [...]})")
    a = ap.parse_args(argv)
    wo = json.loads(Path(a.wall_offset).read_text()) if a.wall_offset else None
    info = prepare(a.problem, a.run_dir, a.steps, op=a.op, wall_offset=wo)
    print(json.dumps({k: info[k] for k in ("design", "moc_forces", "mesh")}, indent=1))
    if a.prepare_only:
        return 0
    rc = run_staged(a.run_dir, a.stages)
    out = collect(a.problem, a.run_dir)
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
