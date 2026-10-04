#!/usr/bin/env python3
r"""C3X / Mark II **翼列 1 ピッチ分の流体メッシュ**を作る (V5 の流体側)。

- 翼面は `ref/vane_*_profile.csv` を**固体メッシュと同じ弧長等間隔 N 点**に再標本化し、
  各区間を Line + `Transfinite Curve = 2` にする → **界面節点が固体と 1 対 1**で揃う。
- 周期境界は **キャンバ線を ±ピッチ/2 だけ平行移動した 2 本**。上流は軸方向、下流は出口空気角方向に延長する。
  forge の周期は**面重心を (dx,dy,dz) 平行移動して最近傍照合**する (`mesh::setPeriodicPartner`) ので、
  2 本の離散が平行移動で一致していることが要件 (gmsh の `Periodic Curve` で保証する)。
- 壁には境界層 (`Field[BoundaryLayer]`) を張る。既定 hwall は $y_1^+\approx1$ 相当の 2 µm。

usage: python3 case/53.c3x_vane_cht/tools/gen_fluid_mesh.py [--vane c3x] [--n-wall 240]
       [--hwall 2e-6] [--ratio 1.15] [--bl-thick 1.2e-3] [--lc 2.5e-3]
"""
import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_solid_mesh import load_profile, resample_closed, parse_msh41   # noqa: E402

GEOM = {   # 表 IV
    "c3x":    dict(pitch=11.773, exit_angle=72.38, case="case/53.c3x_vane_cht"),
    "markii": dict(pitch=12.974, exit_angle=70.96, case="case/54.markii_vane_cht"),
}


def camber_line(P, m=60):
    """前縁→後縁のキャンバ線 (両面を正規化弧長で対応させた中点)。単位は入力と同じ [cm]。"""
    i_le = int(np.argmin(P[:, 0])); i_te = int(np.argmax(P[:, 0]))
    n = len(P)
    fwd = [(i_le + k) % n for k in range(((i_te - i_le) % n) + 1)]
    bwd = [(i_le - k) % n for k in range(((i_le - i_te) % n) + 1)]

    def arc(idx):
        q = P[idx]
        d = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(q, axis=0).T))])
        return q, d / d[-1]

    qa, sa = arc(fwd)
    qb, sb = arc(bwd)
    t = np.linspace(0.0, 1.0, m)
    A = np.column_stack([np.interp(t, sa, qa[:, 0]), np.interp(t, sa, qa[:, 1])])
    B = np.column_stack([np.interp(t, sb, qb[:, 0]), np.interp(t, sb, qb[:, 1])])
    return 0.5 * (A + B)


def periodic_path(cam, x_in, x_out, exit_angle_deg):
    """キャンバ線を上流 (軸方向) と下流 (出口空気角) に延長した 1 本の折れ線 [cm]。"""
    le, te = cam[0], cam[-1]
    up = np.array([[x_in, le[1]]])
    # 出口空気角は軸方向からの角度。翼列は y 方向にピッチを持つので dy/dx = -tan(angle)
    th = math.radians(exit_angle_deg)
    dy = -math.tan(th) * (x_out - te[0])
    dn = np.array([[x_out, te[1] + dy]])
    return np.vstack([up, cam, dn])


def write_geo(path, wall, low, pitch, x_in, x_out, hwall, ratio, thick, lc):
    """wall: 翼面点 (cm), low: 下側周期線 (cm), pitch [cm]"""
    L = [f"lc = {lc};", "// C3X cascade passage (cm -> m)"]
    pid = 1
    wid = []
    for (x, y) in wall:
        L.append(f"Point({pid}) = {{{x/100:.9f}, {y/100:.9f}, 0, lc}};"); wid.append(pid); pid += 1
    lid = []
    for (x, y) in low:
        L.append(f"Point({pid}) = {{{x/100:.9f}, {y/100:.9f}, 0, lc}};"); lid.append(pid); pid += 1
    uid = []
    for (x, y) in low:
        L.append(f"Point({pid}) = {{{x/100:.9f}, {(y+pitch)/100:.9f}, 0, lc}};"); uid.append(pid); pid += 1

    cid = 1
    wall_c = []
    for k in range(len(wid)):
        L.append(f"Line({cid}) = {{{wid[k]}, {wid[(k+1) % len(wid)]}}};"); wall_c.append(cid); cid += 1
    low_c, up_c = [], []
    for k in range(len(lid) - 1):
        L.append(f"Line({cid}) = {{{lid[k]}, {lid[k+1]}}};"); low_c.append(cid); cid += 1
    for k in range(len(uid) - 1):
        L.append(f"Line({cid}) = {{{uid[k]}, {uid[k+1]}}};"); up_c.append(cid); cid += 1
    L.append(f"Line({cid}) = {{{lid[0]}, {uid[0]}}};"); inlet_c = cid; cid += 1
    L.append(f"Line({cid}) = {{{lid[-1]}, {uid[-1]}}};"); outlet_c = cid; cid += 1

    L.append("Curve Loop(1) = {" + ",".join(str(c) for c in low_c) + f",{outlet_c},"
             + ",".join(f"-{c}" for c in reversed(up_c)) + f",-{inlet_c}}};")
    L.append("Curve Loop(2) = {" + ",".join(str(c) for c in wall_c) + "};")
    L.append("Plane Surface(1) = {1, 2};")
    L.append("Transfinite Curve {" + ",".join(str(c) for c in wall_c) + "} = 2;")
    # 周期: 下側 -> 上側 を (0, pitch) 平行移動
    L.append("Periodic Curve {" + ",".join(str(c) for c in up_c) + "} = {"
             + ",".join(str(c) for c in low_c) + f"}} Translate{{0, {pitch/100:.9f}, 0}};")
    L.append(f"Field[1] = BoundaryLayer;")
    L.append("Field[1].EdgesList = {" + ",".join(str(c) for c in wall_c) + "};")
    L.append(f"Field[1].hwall_n = {hwall};")
    L.append(f"Field[1].ratio = {ratio};")
    L.append(f"Field[1].thickness = {thick};")
    L.append("Field[1].Quads = 1;")
    L.append("BoundaryLayer Field = 1;")
    L.append('Physical Curve("inlet", 1) = {' + str(inlet_c) + "};")
    L.append('Physical Curve("outlet", 2) = {' + str(outlet_c) + "};")
    L.append('Physical Curve("per_low", 3) = {' + ",".join(str(c) for c in low_c) + "};")
    L.append('Physical Curve("per_up", 4) = {' + ",".join(str(c) for c in up_c) + "};")
    L.append('Physical Curve("wall", 5) = {' + ",".join(str(c) for c in wall_c) + "};")
    L.append('Physical Surface("fluid", 6) = {1};')
    L.append("Mesh.MshFileVersion = 4.1;")
    path.write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="c3x", choices=list(GEOM))
    ap.add_argument("--n-wall", type=int, default=240)
    ap.add_argument("--hwall", type=float, default=2.0e-6)
    ap.add_argument("--ratio", type=float, default=1.15)
    ap.add_argument("--bl-thick", type=float, default=1.2e-3)
    ap.add_argument("--lc", type=float, default=2.5e-3)
    ap.add_argument("--x-in", type=float, default=None, help="入口 x [cm] (既定: 前縁 - 1 軸弦長)")
    ap.add_argument("--no-convert", action="store_true")
    ap.add_argument("--x-out", type=float, default=None, help="出口 x [cm] (既定: 後縁 + 1.5 軸弦長)")
    # GEOM の exit_angle は報告 表 IV の設計値。**報告の翼型座標そのものから測った
    # 喉/ピッチは別の角度を与える** (C3X: o/pitch 0.2881 -> 73.26° に対し表 IV は 72.38°)。
    # **この 0.88° の食い違いは解に効かない** (2026-09-20 A/B: run_0025 72.38° と
    # run_0026 73.26° で壁圧の bias/rms が 4 桁まで同一、出口気流角も 73.07 / 73.06°)。
    # 周期の出口ブロックは流れの向きを決めず、流れが自分で出口角を決める。
    # したがって後縁側に残る壁圧欠損の原因ではない (case/53 README「負圧面後縁の壁圧」)。
    ap.add_argument("--curv-smooth", type=int, default=0,
                    help="壁節点の再標本化に使う曲率を平滑化する回数 (0=従来)。"
                         "曲率推定のノイズが接線間隔の 2 節点交番を作り、解の市松を強制する "
                         "(case/53 README「メッシュ由来の強制」)。80 で 2.78%% -> 0.86%%")
    ap.add_argument("--exit-angle", type=float, default=None,
                    help="出口気流角 [deg] を上書き (既定は GEOM の表 IV 設計値)")
    a = ap.parse_args()

    g = GEOM[a.vane]
    P = load_profile(a.vane)
    wall, arc = resample_closed(P, a.n_wall, curv_smooth=a.curv_smooth)
    bx = P[:, 0].max() - P[:, 0].min()
    x_in = a.x_in if a.x_in is not None else P[:, 0].min() - 1.0 * bx
    x_out = a.x_out if a.x_out is not None else P[:, 0].max() + 1.5 * bx
    cam = camber_line(P)
    exit_angle = a.exit_angle if a.exit_angle is not None else g["exit_angle"]
    print(f"[mesh] exit_angle = {exit_angle:.2f} deg"
          f"{' (表 IV 既定)' if a.exit_angle is None else ' (上書き)'}")
    path = periodic_path(cam, x_in, x_out, exit_angle)
    low = path - np.array([0.0, g["pitch"] / 2])

    mesh_dir = ROOT / g["case"] / "mesh"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    geo = mesh_dir / f"fluid_{a.vane}.geo"
    msh = mesh_dir / f"fluid_{a.vane}.msh"
    write_geo(geo, wall, low, g["pitch"], x_in, x_out, a.hwall, a.ratio, a.bl_thick, a.lc)
    print(f"[{a.vane}] domain x {x_in:.2f} .. {x_out:.2f} cm, pitch {g['pitch']} cm, "
          f"wall {a.n_wall} nodes (spacing {arc/a.n_wall*10:.2f} mm), hwall {a.hwall*1e6:.1f} um")
    r = subprocess.run(["gmsh", "-2", str(geo), "-o", str(msh), "-format", "msh41", "-v", "2"],
                       capture_output=True, text=True)
    sys.stdout.write(r.stdout[-1500:]); sys.stderr.write(r.stderr[-1500:])
    if r.returncode != 0:
        sys.exit(f"gmsh failed ({r.returncode})")
    xy, tris, phys = parse_msh41(msh)
    print(f"[{a.vane}] mesh: {len(xy)} nodes, {len(tris)} triangles; "
          + ", ".join(f"{k}:{len(v)} edges" for k, v in sorted(phys.items())))

    if a.no_convert:
        return
    # ---- forge 形式へ変換 (node 離散化。wall_dist は no-slip 壁から取る) ----
    import os
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:"
               + os.environ.get("LD_LIBRARY_PATH", ""))
    pitch_m = g["pitch"] / 100.0
    (mesh_dir / "solverConfig.yaml").write_text(
        'mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, '
        f'meshFileName: "fluid_{a.vane}.h5", valueFileName: "fluid_{a.vane}.h5"}}\n'
        "gpu: 1\nsolver: \"SLAU\"\n"
        "physProp: {thermalMethod: 0, viscMethod: 1, visc: 3.3e-5, thermCond: 0.05, "
        "thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}\n"
        "time:\n  unsteady: 0\n  dualTime: 0\n  last: {nStepOuter: 1}\n"
        "  deltaT: {control: 1, dt: 1e-9, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, "
        "dt_min: 1e-12, dt_max: 1.0}\n  outStepStart: 0\n  outStepInterval: 1\n"
        "  timeIntegration: 11\n  nStepInner: 5\n"
        "space: {convMethod: 0, limiter: 0}\nturbulence: {model: \"sst\", wallTreatmentSST: 0}\n"
        'initial: "uniform_p101325_u10"\n')
    (mesh_dir / "bcondConfig.yaml").write_text(
        "inlet:   {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , "
        "floats: {Pt: 319500.0, Tt: 786.0, k: 57.4, omega: 1384.0}}\n"
        "outlet:  {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , "
        "floats: {Ps: 188760.0, Pt: 319500.0, Tt: 786.0}}\n"
        f"per_low: {{physID: 3, kind: periodic, outputHDFflg: 0, ints: {{type: 0, partnerBCID: 4}}, "
        f"floats: {{dx: 0.0, dy: {pitch_m:.9f}, dz: 0.0}}}}\n"
        f"per_up:  {{physID: 4, kind: periodic, outputHDFflg: 0, ints: {{type: 0, partnerBCID: 3}}, "
        f"floats: {{dx: 0.0, dy: {-pitch_m:.9f}, dz: 0.0}}}}\n"
        "wall:    {physID: 5, kind: wall_isothermal, outputHDFflg: 1, ints: , "
        "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 560.0}}\n")
    conv = ROOT / "solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge"
    r2 = subprocess.run([str(conv), msh.name, f"fluid_{a.vane}.h5"], cwd=mesh_dir, env=env,
                        capture_output=True, text=True)
    sys.stdout.write(r2.stdout[-800:]); sys.stderr.write(r2.stderr[-800:])
    if r2.returncode != 0:
        sys.exit(f"converter failed ({r2.returncode})")
    subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_mesh_quality.py"),
                    str(mesh_dir / f"fluid_{a.vane}.h5")], env=env)


if __name__ == "__main__":
    main()
