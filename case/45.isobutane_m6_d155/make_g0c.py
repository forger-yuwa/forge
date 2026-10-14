"""G0c (G0 の縮流部 x −8〜−2 r_t だけ流れ方向の密度を max(1, d/2) に並べ直した格子、ni 4467) の .msh を作る。
plan architecture-float-state-double-geometry §5.1 #22 (forge-25 の Grefc の作り方を G0 に当てる。ユーザ決定: 診断に限り AR ≤ 10,000)。
forge-25 の関数 (feature/nozzle-core-grid の case/45.isobutane_m6_d155/core_grid_mesh.py、commit を引数で固定) を git show で読み込み、
PROBLEM・WALL_REF をこの作業ツリーの同じファイルに向ける。G0 (ni 4719) も同じ手順で作り、区間の外の列が G0 とビット一致することを確かめる。
使い方 (design の venv で): /home/sano/work/forge/design/.venv-opt/bin/python make_g0c.py <commit> <出力ディレクトリ>
出力: <出力>/nozzle.msh (G0c)、<出力>/g0c_info.json (列の数・区間・ビット一致の検査・第一層の厚さ)"""
import importlib.util, json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
commit, outdir = sys.argv[1], Path(sys.argv[2]); outdir.mkdir(parents=True, exist_ok=True)
src = subprocess.run(["git", "-C", str(HERE), "show", f"{commit}:case/45.isobutane_m6_d155/core_grid_mesh.py"], capture_output=True, text=True, check=True).stdout
mod_path = outdir / "_core_grid_mesh_from_forge25.py"; mod_path.write_text(src)
sys.argv = ["core_grid_mesh.py", "select"]
spec = importlib.util.spec_from_file_location("core_grid_mesh_f25", mod_path); CG = importlib.util.module_from_spec(spec)
sys.path.insert(0, str(HERE.parents[1] / "design"))
spec.loader.exec_module(CG)
CG.PROBLEM = HERE / "problem_d155_ns_prod_coldmesh_tw300.yaml"; CG.WALL_REF = HERE / "_band_ab" / "prod_confirm" / "prep"
from forge_design.meshing import mesh2d  # noqa: E402
wall = CG.Wall(CG.WALL_REF)
m = CG.mesh_block()
NI0, NJ = int(m["ni"]), int(m["nj"])
assert NI0 == 4719 and NJ == 121, (NI0, NJ)
# G0 (元の station)
CG.use_grefc_stations(False)
prm0 = CG._mesh_params_from(m, wall.scale, NI0, NJ, float(m["wall_first_frac"]))
c0, q0, b0 = mesh2d.generate_axisym_mesh(wall, prm0)
# G0c
m1 = dict(m); m1["ni"] = 4467
CG.use_grefc_stations(True)
prm1 = CG._mesh_params_from(m1, wall.scale, 4467, NJ, float(m1["wall_first_frac"]))
c1, q1, b1 = mesh2d.generate_axisym_mesh(wall, prm1)
CG.use_grefc_stations(False)
P0 = c0[:, :2].reshape(NI0, NJ, 2) / wall.scale; P1 = c1[:, :2].reshape(4467, NJ, 2) / wall.scale
x0, x1 = P0[:, 0, 0], P1[:, 0, 0]
# 区間の外の列のビット一致: 前 ia 列と後ろの列
ia = int(np.argmin(np.abs(x0 - (-8.0)))); ib = int(np.argmin(np.abs(x0 - (-2.0)))); nseg = 4467 - (NI0 - (ib - ia))
pre_ok = bool(np.array_equal(c0.reshape(NI0, NJ, -1)[:ia], c1.reshape(4467, NJ, -1)[:ia]))
post_ok = bool(np.array_equal(c0.reshape(NI0, NJ, -1)[ib + 1:], c1.reshape(4467, NJ, -1)[ia + nseg + 1:]))
dx1 = np.diff(x1); dx0 = np.diff(x0)
info = dict(commit=commit, msh_digits=int(m["msh_digits"]), ni_G0=NI0, ni_G0c=4467, nj=NJ, seg_cols_G0=int(ib - ia), seg_cols_G0c=int(nseg), seg_x=[float(x0[ia]), float(x0[ib])],
            pre_cols_bit_equal=pre_ok, post_cols_bit_equal=post_ok, n_pre=int(ia), n_post=int(NI0 - ib - 1),
            joint_dx_ratio=[float(dx1[ia] / dx1[ia - 1]), float(dx1[ia + nseg] / dx1[ia + nseg - 1])],
            seg_dx_max_ratio_vs_G0=float(np.max(dx1[ia:ia + nseg]) / np.max(dx0[ia:ib])))
# 座標の桁数は生産 (runner_axismach.py の write_msh41_2d) と同じく問題 YAML の mesh.msh_digits (冷却壁は 17)。
# 2026-10-11: これを渡さず既定 10 桁で書き、実際の G0 の入力と区間の外の座標が相対 ±5e-10 ずれた (壁際で局所の間隔の 0.13 %、g0c_check.py で発見)
digits = int(m["msh_digits"]); assert digits == 17, digits
mesh2d.write_msh41_2d(str(outdir / "nozzle.msh"), c1, q1, b1, digits=digits)
(outdir / "g0c_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
print(json.dumps(info, ensure_ascii=False))
if not (pre_ok and post_ok): sys.exit("区間の外の列が G0 とビット一致しない")
