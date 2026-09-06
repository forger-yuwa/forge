"""⑤ SERN の 3D 確認 runner (S7): 2D 設計 (key point 逆 MOC) を有限スパン + 側壁の 3D hex メッシュで評価する。

問題型 `sern_2d` の YAML に `mesh3d:` 節 (W, Z_ext, L_sw, nz_in, nz_out, ni_*, nj_*, first_wall_frac) を足して使う。
IC は y > y_mid(x) かつ z ≤ W/2 を排気状態、他を外部流。力は壁面出力の quad ごとに (p − p_a) n A を積分し、
半スパン (W/2) あたりの 2D 換算係数 C_T, C_L, C_M を出す (ランプは z ≤ W/2 [ノズル幅内] と z > W/2 [外側] に分けて記録)。
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

from ..geometry.rao_planar import ideal_gross_thrust
from ..meshing.mesh_sern3d import PHYS_SERN3D, SernMesh3DParams, generate_sern_mesh3d, write_msh41_3d
from ..metrics.sern_forces import steadiness
from ..probdef import load_problem
from . import runner_sern as R2

MESH = "sern.h5"


def _solver_config(p, nsteps, out_int, cfl, p_ref):
    cfg = R2._solver_config(p, nsteps, out_int, cfl, p_ref)
    return cfg  # 2D と同じ (isAxisymmetric: 0, discretization は mesh.discretization)


def _bcond_config(p, st):
    model = p.evaluate.get("model", "euler"); wall_kind = "slip" if model == "euler" else "wall"
    ex, en = st["exhaust"], st["ext"]; P = PHYS_SERN3D

    def inlet(name, s):
        return (f"{name}: {{physID: {P[name]}, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , "
                f"floats: {{ro: {s['ro']:.6g}, Ux: {s['u']:.6g}, Uy: 0.0, Uz: 0.0, Ps: {s['P']:.6g}, k: {s['k']:.6g}, omega: {s['omega']:.6g}}}}}\n")

    def outlet(name):
        return f"{name}: {{physID: {P[name]}, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {en['P']:.6g}, Pt: {en['P']:.6g}, Tt: {en['T']:.6g}}}}}\n"

    def wall(name, kind=None):
        return f"{name}: {{physID: {P[name]}, kind: {kind or wall_kind}, outputHDFflg: 1, ints: , floats: }}\n"
    return (inlet("inlet_nozzle", ex) + inlet("inlet_ext", en) + outlet("outlet") + wall("ramp") + wall("cowl_in") + wall("cowl_out")
            + outlet("bottom") + (outlet("top_out") if p.evaluate.get("top_out_kind", "outlet") == "outlet"
                                  else f"top_out: {{physID: {P['top_out']}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n")
            + f"sym: {{physID: {P['sym']}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n"
            + f"side_far: {{physID: {P['side_far']}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n"
            + wall("sidewall_in") + wall("sidewall_out"))


def paste_region_ic3d(h5path, y_mid, scale, half_W_m, st, gamma, minfo=None):
    """領域別一様 IC。node では CELLS/centCoords が双対重心で、側壁面 (z=W/2) 上のノードは内外どちらとも
    決めにくい (ランプ∩側壁の共有ノードなど) ので、`minfo` があれば**ノード index** で判定する:
    排気 = 上バンド (j ≥ jm) かつ k ≤ k_sw の base ノード + カウル上コピー (dup1)、側壁外コピー (dup2) と他は外部流。"""
    ex, en = st["exhaust"], st["ext"]
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn, yn = cc[:, 0] / scale, cc[:, 1] / scale
        if minfo is not None and len(cc) == minfo["nodes"]:
            NJ, nz, jm, ksw = minfo["NJ"], minfo["nz"], minfo["jm"], minfo["k_sw"]
            N_base = minfo["ni"] * NJ * nz; ids = np.arange(len(cc))
            k = ids % nz; j = (ids // nz) % NJ; i = ids // (nz * NJ); ite = minfo["i_te"]
            # 中間線 (j = jm) の元ノードはカウル区間 (i < i_te) では cowl_out = 外部流。排気側は dup1 コピーが担う。
            # カウル下流 (i ≥ i_te) の共有中間線ノードだけ排気に入れる (プリューム境界)
            upper = (ids < N_base) & (k <= ksw) & ((j > jm) | ((j == jm) & (i >= ite)))
            upper |= (ids >= N_base) & (ids < N_base + minfo["n_dup_cowl"])
        else:
            upper = (yn > y_mid(xn)) & (cc[:, 2] <= half_W_m * (1 + 1e-9))
        ro = np.where(upper, ex["ro"], en["ro"]); u = np.where(upper, ex["u"], en["u"]); P = np.where(upper, ex["P"], en["P"])
        v = f["/VALUE"]
        v["ro"][:] = ro; v["roUx"][:] = ro * u; v["roUy"][:] = 0.0; v["roUz"][:] = 0.0
        v["roe"][:] = P / (gamma - 1.0) + 0.5 * ro * u * u
        if "roK" in v:
            v["roK"][:] = ro * np.where(upper, ex["k"], en["k"]); v["roOmega"][:] = ro * np.where(upper, ex["omega"], en["omega"])


def prepare(problem_path, run_dir, nsteps=None, op=None) -> dict:
    p = load_problem(problem_path)
    run_dir = Path(run_dir); run_dir.mkdir(parents=True, exist_ok=False)
    d0 = R2.design_snapshot(p); opinfo = R2.select_operating_point(p, op); st = R2.gas_states(p)
    kern, design, fr_moc, theta_b = R2.design_from_problem(p, design=d0)
    H = float(p.spec["H_m"]); m = p.raw.get("mesh3d", {}); m2 = p.mesh
    prm = SernMesh3DParams(ni_up=int(m.get("ni_up", 10)), ni_noz=int(m.get("ni_noz", 60)), ni_plume=int(m.get("ni_plume", 110)),
                           nj_top=int(m.get("nj_top", 49)), nj_bot=int(m.get("nj_bot", 31)), nz_in=int(m.get("nz_in", 25)), nz_out=int(m.get("nz_out", 17)),
                           W=float(m.get("W", 2.0)), Z_ext=float(m.get("Z_ext", 1.5)), L_sw=m.get("L_sw"), L_up=float(m2.get("L_up", 0.5)),
                           x_out_extra=float(m2.get("x_out_extra", 2.0)), bot_depth=float(m2.get("bot_depth", 3.0)),
                           first_wall_frac=float(m.get("first_wall_frac", m2.get("first_wall_frac", 4e-3))), first_z_frac=float(m.get("first_z_frac", 4e-3)),
                           interface_angle=float(m2.get("interface_angle_rad", theta_b)),
                           top_ext_angle=float(np.deg2rad(m2.get("top_ext_angle_deg", np.rad2deg(design.info["theta_e"])))), scale=H)
    coords, hexes, B, minfo, y_mid = generate_sern_mesh3d(design, prm)
    write_msh41_3d(run_dir / "sern.msh", coords, hexes, B, PHYS_SERN3D)
    np.savetxt(run_dir / "ramp_contour.csv", design.ramp_xy * H, delimiter=",", header="x_m,y_m", comments="")
    np.savetxt(run_dir / "cowl_contour.csv", design.cowl_xy * H, delimiter=",", header="x_m,y_m", comments="")
    n = int(nsteps or p.evaluate.get("nStepOuter", 4000)); out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl = float(p.evaluate.get("cfl_main", 1.0)); cfg = _solver_config(p, n, out_int, cfl, st["ext"]["P"])
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p, st)); (run_dir / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    disc = p.mesh.get("discretization", "node")
    (run_dir / "solverConfig.yaml").write_text(cfg.replace(f'discretization: "{disc}"', 'discretization: "cell"').replace(", nodeWallDirichlet: 1", ""))
    R2.convert_mesh(run_dir, "sern.msh", "sern_qc.h5")
    q = subprocess.run([sys.executable, str(R2.FORGE_TOOLS / "check_mesh_quality.py"), "sern_qc.h5", "--mode", "3d"], cwd=run_dir, env=R2._ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    if disc == "cell":
        (run_dir / "sern_qc.h5").rename(run_dir / MESH)
    else:
        (run_dir / "sern_qc.h5").unlink(); (run_dir / "solverConfig.yaml").write_text(cfg)
        R2.convert_mesh(run_dir, "sern.msh", MESH)
    for f in run_dir.glob("sern_qc.xmf"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg); (run_dir / "solverConfig_main.yaml").write_text(cfg)
    paste_region_ic3d(run_dir / MESH, y_mid, H, 0.5 * prm.W * H, st, p.gamma, minfo=minfo if disc == "node" else None)
    ex = st["exhaust"]; F_ideal_nd, M_e_id = ideal_gross_thrust(ex["M"], st["ext"]["P"] / ex["P"], p.gamma)
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n, "H_m": H, "states": st, "operating_point": opinfo,
            "design_point": d0, "design": {"key_point": list(design.key_point), "L_ramp": design.L_ramp, "theta_e_deg": float(np.rad2deg(design.info["theta_e"])),
            "theta_b_deg": float(np.rad2deg(theta_b)), "warnings": design.info["warnings"]}, "moc_forces": fr_moc,
            "F_ideal_N_per_m": F_ideal_nd * ex["P"] * H, "mesh": minfo, "discretization": disc, "model": p.evaluate.get("model", "euler"), "dim": 3,
            "half_W_m": 0.5 * prm.W * H}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(run_dir, stages="full", soft_steps=2000, soft_cfl=0.5, soft_conv=0,
               warm_lam_steps=0, warm_lam_cfl=0.2, mid_steps=0):
    return R2.run_staged(run_dir, stages, soft_steps, soft_cfl=soft_cfl, soft_conv=soft_conv,
                         warm_lam_steps=warm_lam_steps, warm_lam_cfl=warm_lam_cfl, mid_steps=mid_steps)


def warm_from_same_mesh(run_dir, res_h5, keep_turb: bool = True) -> None:
    """同一メッシュの別 run (例: Euler) の場を index コピーして暖機起動にする。keep_turb=True なら roK/roOmega は
    現在の IC を保つ (Euler 場には無い)。座標一致ノードがあるので最近傍補間は使わない。"""
    keys = ["ro", "roUx", "roUy", "roUz", "roe"] + ([] if keep_turb else ["roK", "roOmega"])
    # 状態量だけをコピーする。res には wall_dist も入っており、Euler run (slip 壁) の wall_dist は番兵値なので
    # 丸ごとコピーすると SST の壁距離が全滅する (run_0026 で実害、2026-09-05)
    with h5py.File(res_h5, "r") as src, h5py.File(Path(run_dir) / MESH, "r+") as dst:
        n = len(dst["VALUE/ro"])
        for k in keys:
            if k in src["VALUE"] and k in dst["VALUE"] and len(src["VALUE"][k]) == n:
                dst["VALUE"][k][:] = src["VALUE"][k][:]


# ---------------------------------------------------------------------- 3D 壁面力
def _parse_conne(flat):
    """壁面出力 MESH/CONNE を面ごとの節点 id リストにする。
    3D (quad): 1 面 5 整数 [5, n0, n1, n2, n3] (case/46 run_0023 で確認)。2D (line): 1 面 4 整数 [2, 2, a, b]。"""
    flat = flat.astype(int)
    if len(flat) % 5 == 0 and np.all(flat[0::5] == 5):
        return list(flat.reshape(-1, 5)[:, 1:5])
    if len(flat) % 4 == 0 and np.all(flat[0::4] == 2):
        return list(flat.reshape(-1, 4)[:, 2:4])
    out = []; i = 0
    while i < len(flat):
        nn = int(flat[i + 1]); out.append(flat[i + 2:i + 2 + nn]); i += 2 + nn
    return out


def surface_forces(h5file, expect_dir, p_a, x_ref, y_ref, z_split=None):
    """quad ごとに (p − p_a)·n·A。n は流体から壁へ向く向きに expect_dir (概ねの向きベクトル) で揃える。
    戻り: 全体と (z_split があれば) z ≤ z_split / z > z_split の分割。"""
    with h5py.File(h5file, "r") as f:
        xyz = f["MESH/COORD"][:].reshape(-1, 3).astype(float); faces = _parse_conne(f["MESH/CONNE"][:]); p = f["VALUE/Ps"][:].astype(float)
        tw = np.column_stack([f["VALUE/twall_x"][:], f["VALUE/twall_y"][:], f["VALUE/twall_z"][:]]).astype(float) if "VALUE/twall_x" in f else None
    nf = len(faces); F = np.zeros((nf, 3)); cen = np.zeros((nf, 3)); area = np.zeros(nf)
    per_node = len(p) == len(xyz)
    for k, fc in enumerate(faces):
        q = xyz[fc]; c = q.mean(0)
        nrm = 0.5 * np.cross(q[2] - q[0], q[3] - q[1]) if len(fc) == 4 else 0.5 * np.cross(q[1] - q[0], q[2] - q[0])
        A = np.linalg.norm(nrm); nvec = nrm / max(A, 1e-30)
        if np.dot(nvec, expect_dir) < 0:
            nvec = -nvec
        pf = p[fc].mean() if per_node else p[k]
        F[k] = (pf - p_a) * nvec * A; cen[k] = c; area[k] = A
    Mz = (cen[:, 0] - x_ref) * F[:, 1] - (cen[:, 1] - y_ref) * F[:, 0]
    out = {"Fx_p": float(F[:, 0].sum()), "Fy_p": float(F[:, 1].sum()), "Fz_p": float(F[:, 2].sum()), "M_noseup_p": float(-Mz.sum()), "area": float(area.sum()), "n_faces": nf}
    if z_split is not None:
        m = cen[:, 2] <= z_split * (1 + 1e-9)
        out["inside"] = {"Fx_p": float(F[m, 0].sum()), "Fy_p": float(F[m, 1].sum()), "M_noseup_p": float(-Mz[m].sum()), "area": float(area[m].sum())}
        out["outside"] = {"Fx_p": float(F[~m, 0].sum()), "Fy_p": float(F[~m, 1].sum()), "M_noseup_p": float(-Mz[~m].sum()), "area": float(area[~m].sum())}
    if tw is not None:
        tf = np.array([tw[fc].mean(0) if per_node else tw[k] for k, fc in enumerate(faces)]) * area[:, None]
        out["Fx_tau"] = float(tf[:, 0].sum()); out["Fy_tau"] = float(tf[:, 1].sum())
    return out


def forces3d(run_dir, step, p_a, F_ideal_per_m, half_W, H, x_ref=0.0, y_ref=0.0, mdot_u_in=0.0, p_in=0.0, twall_on_fluid=False):
    P = PHYS_SERN3D; run_dir = Path(run_dir)
    spec = {"ramp": (P["ramp"], (0, 1, 0)), "cowl_in": (P["cowl_in"], (0, -1, 0)), "cowl_out": (P["cowl_out"], (0, 1, 0)),
            "sidewall_in": (P["sidewall_in"], (0, 0, 1)), "sidewall_out": (P["sidewall_out"], (0, 0, -1))}
    parts = {}
    for name, (pid, d) in spec.items():
        c = list(run_dir.glob(f"res_{name}_{pid}_{step}.h5"))
        if c:
            parts[name] = surface_forces(c[0], np.array(d, float), p_a, x_ref, y_ref, z_split=half_W if name == "ramp" else None)
    Fx = sum(v["Fx_p"] for v in parts.values()); Fy = sum(v["Fy_p"] for v in parts.values()); Mn = sum(v["M_noseup_p"] for v in parts.values())
    F_ideal = F_ideal_per_m * half_W; inlet = (mdot_u_in + (p_in - p_a) * H) * half_W
    out = {"step": step, "T_wall": -Fx, "L": Fy, "M_noseup": Mn, "F_gross": inlet - Fx, "F_ideal": F_ideal, "C_T": (inlet - Fx) / F_ideal,
           "C_T_wall": -Fx / F_ideal, "C_L": Fy / F_ideal, "C_M": Mn / (F_ideal * H), "parts": parts}
    if "ramp" in parts and "inside" in parts["ramp"]:
        ri, ro = parts["ramp"]["inside"], parts["ramp"]["outside"]
        out["C_L_ramp_inside"] = ri["Fy_p"] / F_ideal; out["C_L_ramp_outside"] = ro["Fy_p"] / F_ideal
        out["C_T_ramp_inside"] = -ri["Fx_p"] / F_ideal; out["C_T_ramp_outside"] = -ro["Fx_p"] / F_ideal
    if any("Fx_tau" in v for v in parts.values()):
        sgn = 1.0 if twall_on_fluid else -1.0
        Fx_t = sgn * sum(v.get("Fx_tau", 0.0) for v in parts.values())
        out["C_T_with_shear"] = (inlet - Fx + Fx_t) / F_ideal; out["C_T_friction"] = Fx_t / F_ideal
    return out


def collect(problem_path, run_dir) -> dict:
    p = load_problem(problem_path); run_dir = Path(run_dir); info = json.loads((run_dir / "prepare_info.json").read_text())
    st = info["states"]; ex, en = st["exhaust"], st["ext"]; H = info["H_m"]; xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    steps = sorted({int(f.stem.split("_")[-1]) for f in run_dir.glob("res_ramp_*_[0-9]*.h5")})
    hist = [forces3d(run_dir, s, en["P"], info["F_ideal_N_per_m"], info["half_W_m"], H, float(xr) * H, float(yr) * H,
                     ex["ro"] * ex["u"] ** 2 * H, ex["P"], twall_on_fluid=(info["discretization"] == "cell")) for s in steps]
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-2:] if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else []
    out = {"convergence_verdict": verdict, "n_snapshots": len(hist), "history": hist, "dim": 3, "moc_forces": info["moc_forces"]}
    if hist:
        last = hist[-1]
        out.update({k: last[k] for k in last if k != "parts"})
        out["steadiness"] = {k: steadiness([h[k] for h in hist]) for k in ("C_T", "C_L", "C_M")}
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="⑤ SERN 3D 確認 run")
    ap.add_argument("problem"); ap.add_argument("run_dir"); ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true"); ap.add_argument("--stages", default="full", choices=["full", "none", "main"])
    ap.add_argument("--op", default=None); ap.add_argument("--soft-steps", type=int, default=2000)
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, a.steps, op=a.op); print(json.dumps({k: info[k] for k in ("design", "mesh", "moc_forces")}, indent=1))
    if a.prepare_only:
        return 0
    # 2D と同じく段階起動は problem の `opt:` を正本にする (§4.12 のレシピを 3D にそのまま移植するため)
    o = (load_problem(a.problem).raw.get("opt") or {})
    rc = run_staged(a.run_dir, a.stages, int(o.get("soft_steps", a.soft_steps)),
                    soft_cfl=float(o.get("soft_cfl", 0.5)), soft_conv=int(o.get("soft_conv", 0)),
                    warm_lam_steps=int(o.get("warm_lam_steps", 0)), warm_lam_cfl=float(o.get("warm_lam_cfl", 0.2)),
                    mid_steps=int(o.get("mid_steps", 0)))
    out = collect(a.problem, a.run_dir)
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1)); return rc


if __name__ == "__main__":
    raise SystemExit(main())
