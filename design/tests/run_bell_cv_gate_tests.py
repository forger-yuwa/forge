#!/usr/bin/env python3
"""③ベル 2 モードの共通評価器 — 推力の CV 会計とゲート (plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2b)。

forge を起動しない単体試験 (全流れ解での確認は case/68 の run):
- `metrics.extract.thrust_cv`: 合成した node メッシュ (TOP ベルの壁・直管) と場で、
  静止気体の閉曲面で Σ(p − p_a) n_x dA = 0 (CV 残差 0、p_a に依らない)、面積 πr²、外向き法線の向き、一様流の直管で
  壁会計 = 出口面会計・ṁ_in = ṁ_out、twall の規約 (node: 流体に働く traction = −twall、cell: +twall) と摩擦の推力寄与の符号、
  P の出所 (VALUE/P 優先、無ければ TP は e_sens(T) の逆算・CPG は γ)、粘性応力・壁ダンプが無いときは閉じない扱い。
  Euler (slip 壁、§5.1 #6 (W1)): `wall_model="slip"` で内壁の traction 0・壁ダンプを読まない・μ_eff = 0 の確認、開境界の粘性項が
  viscMethod 0 (vis_lam = 0) で 0、NS で壁ダンプに twall が無ければ止める、粘性の場に slip を使えば止める。
- `opt.driver` のゲート (i)〜(v) (2026-10-11 組み替え、plan §4.1・§5.1 #2c'): 合成した残差履歴 + StageManifest で (i) は PASS・
  STALLED・still converging を通し (STALLED は併記)、DIVERGED・有意な末尾上昇 (rising)・manifest 無し (判定不能) を落とす、
  (ii) 推力系 4 列 (1.5e-4) の DRIFTING・振動・単調で漸近値が離れた列・NONFINITE・列欠損・末尾窓 8 点未満を落とし、ṁ は 1e-3 で判定、
  cv_residual は列にしない、(iii) 物理範囲、(iv) 最終 |CV 残差| 1e-4 超過・閉じない会計を落とす、(v) メッシュ品質の VERDICT の読みと
  AR > 1000 の格子の例外の条件 (MESH_AR_RELAX.json)、総合の分類 (RETRYABLE / SUSPECT / INFEASIBLE)。
- 壁法線の近壁層 (#2d): mesh2d の `wall_normal_end_stations` (入口・出口面を平面に)、runner の格子キーの受け渡し、
  `mesh_ar_relax_check` (AR > 1000 のスキュー・面の平面性) の負例。
- (D) 直接モード (#3f): `geometry.bell_mode: direct_bspline` の壁 (`DirectBellWall`)、成立しない点の `WallInfeasible` (理由つき)、
  TOP のキー・dv の誤りを止める、driver の設計変数の並び。
- (M) 特性曲線モード (#6 (W1)): `geometry.bell_mode: moc_rao` の壁 (`RaoBellWall` = 生成器の壁ファイルの `RaoWall` + (D) と同じ
  収縮部) が npz の壁点を通る (kernel は円弧の上、B・E は厳密、BE の内点は当てはめの残差以内)、T4 の評価と同じ作り方、スロートで C2・
  B で C1、L・ε・R_t・排気ファイルの不一致で止める、prepare の変換器の手前まで (偽の変換器) で Euler 格子の品質 PASS と
  WALL_DESIGN.json の V2 参照値 (生成器の推力を N・η に換算)。合成の Euler run で推力の時系列・metrics (η − η_cs)・ゲート
  (壁解像は対象外と書いて ok、他は NS と同じ)、geometry.rao_wall_file の相対パスの書き換え。
- 引き継ぎ (b): 評価ごとに書き出す問題 YAML の相対パス (gas.exhaust・gas.species_db) を絶対パスに直す。
- runner: TP の段の計画 (`tp_stage_plan`) の検査、TP の config の出力 (勾配の extraFields・pRef を書かない)、
  格子のスロート細分キー (`mesh_params`) は書いたときだけ渡る。

実行: design/.venv-opt の Python で `python design/tests/run_bell_cv_gate_tests.py`。
"""
import csv
import json
import math
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate import runner  # noqa: E402
from forge_design.metrics import extract as EX  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh  # noqa: E402
from forge_design.opt import driver as DR  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CASE68 = ROOT / "case" / "68.bell_two_mode"
PROBLEM_TP = CASE68 / "problem_bell_tp_top.yaml"
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))
from stage_manifest import StageManifest  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:
        return str(e) or type(e).__name__
    return None


# ---------------------------------------------------------------- 合成 node メッシュ / 場 / 壁ダンプ
IDS = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4}


def write_mesh_h5(path, coords, quads, bedges, viz=True, planes=True, flip_planes=False):
    """forge の node 変換 h5 と同じ並びの最小構成 (thrust_cv が読む分だけ)。入口の壁角ノードは壁 bcond に入れる
    (nodeInletCornerWall と同じ帰属)。`planes`: 境界双対面 (bcond の節点ごとに外向きの半辺ベクトルの和) を
    PLANES/surfVect・BCONDS/<id>/iPlanes に書く (`flip_planes` で向きを逆にした壊れたメッシュを作る)。"""
    # node の CELLS/centCoords は双対 CV 重心 (境界ノードでは領域の内側へずれる)。隣接四角形の重心の平均で近似する
    qc = coords[quads].mean(axis=1)
    cc = np.zeros_like(coords); cnt = np.zeros(len(coords))
    for k in range(4):
        np.add.at(cc, quads[:, k], qc); np.add.at(cnt, quads[:, k], 1.0)
    cc /= np.maximum(cnt, 1.0)[:, None]
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float64).ravel()
        f["CELLS/centCoords"] = cc.astype(np.float64).ravel()
        if viz:
            f["VIZMESH/CONNE"] = np.column_stack([np.full(len(quads), 5), quads]).astype(np.int32).ravel()
        wall_nodes = {n for e in bedges["wall"] for n in e}
        # 境界辺の外向き (辺を含む四角形の重心から離れる向き) を試験側で独立に作る
        ekey = {}
        for qi, qd in enumerate(quads):
            for k in range(4):
                ekey[frozenset((int(qd[k]), int(qd[(k + 1) % 4])))] = qi
        sv_all = []
        for name, pid in IDS.items():
            e = np.array(bedges[name], dtype=np.int32)
            nodes = sorted({int(n) for n in e.ravel()})
            if name in ("inlet", "outlet"):
                nodes = [n for n in nodes if not (name == "inlet" and n in wall_nodes)]
            f[f"BCONDS/{pid}/vizBfaceNodes"] = e.ravel()
            f[f"BCONDS/{pid}/vizBfaceSizes"] = np.full(len(e), 2, dtype=np.int32)
            f[f"BCONDS/{pid}/iCells"] = np.array(nodes, dtype=np.int32)
            sv = {nd: np.zeros(2) for nd in nodes}
            for ea, eb in e:
                dd = coords[eb, :2] - coords[ea, :2]
                nv = np.array([dd[1], -dd[0]])
                if np.dot(nv, 0.5 * (coords[ea, :2] + coords[eb, :2]) - qc[ekey[frozenset((int(ea), int(eb)))], :2]) < 0:
                    nv = -nv
                for nd in (int(ea), int(eb)):
                    if nd in sv:
                        sv[nd] += 0.5 * nv * (-1.0 if flip_planes else 1.0)
            f[f"BCONDS/{pid}/iPlanes"] = np.arange(len(sv_all), len(sv_all) + len(nodes), dtype=np.int32)
            sv_all += [np.r_[sv[nd], 0.0] for nd in nodes]
        if planes:
            f["PLANES/surfVect"] = np.array(sv_all).ravel()
        else:
            for pid in IDS.values():
                del f[f"BCONDS/{pid}/iPlanes"]


def write_res_h5(path, n, ro, ux, ur, P=None, grads=True, roe=None):
    with h5py.File(path, "w") as f:
        V = f.create_group("VALUE")
        V["ro"] = np.broadcast_to(ro, n).astype(np.float64)
        V["Ux"] = np.broadcast_to(ux, n).astype(np.float64)
        V["Uy"] = np.broadcast_to(ur, n).astype(np.float64)
        V["Uz"] = np.zeros(n)
        V["roUx"] = V["ro"][:] * V["Ux"][:]
        V["roUy"] = V["ro"][:] * V["Uy"][:]
        V["roUz"] = np.zeros(n)
        if P is not None:
            V["P"] = np.broadcast_to(P, n).astype(np.float64)
        if roe is not None:
            V["roe"] = np.broadcast_to(roe, n).astype(np.float64)
        if grads:
            for g in ("dUxdx", "dUxdy", "dUydx", "dUydy", "axisym_divU"):
                V[g] = np.zeros(n)
            V["vis_lam"] = np.full(n, 1.0e-4)
            V["vis_turb"] = np.zeros(n)


def write_wall_dump(path, coords, wall_nodes, tx, ty):
    wn = np.array(wall_nodes)
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords[wn].astype(np.float64).ravel()
        f["VALUE/twall_x"] = np.broadcast_to(tx, len(wn)).astype(np.float64)
        f["VALUE/twall_y"] = np.broadcast_to(ty, len(wn)).astype(np.float64)


def pipe_mesh(L=2.0, R=0.5, ni=41, nj=17):
    xs, rs = np.linspace(0, L, ni), np.linspace(0, R, nj) ** 1.0
    X, Rr = np.meshgrid(xs, rs, indexing="ij")
    coords = np.column_stack([X.ravel(), Rr.ravel(), np.zeros(X.size)])
    nid = lambda i, j: i * nj + j  # noqa: E731
    quads = np.array([(nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1))
                      for i in range(ni - 1) for j in range(nj - 1)])
    bedges = {"axis": [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)],
              "wall": [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)],
              "inlet": [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)],
              "outlet": [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)]}
    return coords, quads, bedges


tmp = Path(tempfile.mkdtemp(prefix="bell_cv_"))

# ---- (1) TOP ベルの壁 (runner と同じ壁・同じメッシュ生成) で静止気体: 閉曲面なので Σ p n_x dA = 0
p = load_problem(PROBLEM_TP)
wall, _ = runner.build_wall(p)
coords, quads, bedges = generate_axisym_mesh(wall, Mesh2DParams(ni=121, nj=31, wall_first_frac=2e-3, scale=0.14))
N = len(coords)
write_mesh_h5(tmp / "bell.h5", coords, quads, bedges)
p0 = 2.0e5
write_res_h5(tmp / "res_100.h5", N, 1.0, 0.0, 0.0, P=p0)
wall_nodes = sorted({n for e in bedges["wall"] for n in e})
write_wall_dump(tmp / "res_wall_3_100.h5", coords, wall_nodes, 0.0, 0.0)
r_in, r_e = float(wall.r(np.array([wall.x_in]))[0]) * 0.14, float(np.sqrt(p.spec["eps"])) * 0.14
F_id = 1.0e6
for pa in (0.0, p0 / 3.0, 1.7e6):
    o = EX.thrust_cv(tmp / "bell.h5", tmp / "res_100.h5", pa, F_id)
    g = o["groups"]
    check(f"(1) 静止気体 p_a={pa:g}: CV 残差 = 0 (閉曲面 Σ(p−p_a) n_x dA = 0)", abs(o["cv_residual"]) < 1e-12,
          f"{o['cv_residual']:.3e}")
    check(f"(1) 静止気体 p_a={pa:g}: 出口面会計 = (p−p_a)·A_e", math.isclose(o["thrust_exit_N"], (p0 - pa) * math.pi * r_e ** 2,
                                                                       rel_tol=1e-12, abs_tol=1e-6))
check("(1) 入口面積 = π r_in² (r 重みの半辺積分が厳密)", math.isclose(g["inlet"]["area"], math.pi * r_in ** 2, rel_tol=1e-12),
      f"{g['inlet']['area']:.12g} vs {math.pi * r_in ** 2:.12g}")
check("(1) 出口面積 = π r_e²", math.isclose(g["outlet"]["area"], math.pi * r_e ** 2, rel_tol=1e-12))
check("(1) 入口の外向き法線は −x (圧力項 = −(p−p_a) A_in)", math.isclose(g["inlet"]["pressure"], -(p0 - 1.7e6) * math.pi * r_in ** 2,
                                                               rel_tol=1e-12))
check("(1) 壁の圧力項 = (p−p_a)(A_in − A_e) (収縮部 +、拡大部 −、正味は面積差)",
      math.isclose(g["wall"]["pressure"], (p0 - 1.7e6) * math.pi * (r_in ** 2 - r_e ** 2), rel_tol=1e-10))
check("(1) VALUE/P を使う", o["pressure_source"] == "VALUE/P")
check("(1) 勾配と壁ダンプがあれば閉じる会計", o["cv_closable"] and o["open_viscous"] == "evaluated")
# VIZMESH が無い場合 (双対 CV 重心で向きを決める) も同じ
write_mesh_h5(tmp / "bell_noviz.h5", coords, quads, bedges, viz=False)
o2 = EX.thrust_cv(tmp / "bell_noviz.h5", tmp / "res_100.h5", 0.0, F_id)
check("(1) VIZMESH 無し (向きは境界双対面のベクトルから) でも同じ値",
      math.isclose(o2["thrust_exit_N"], p0 * math.pi * r_e ** 2, rel_tol=1e-12) and abs(o2["cv_residual"]) < 1e-12)
write_mesh_h5(tmp / "bell_flip.h5", coords, quads, bedges, flip_planes=True)
check("(1) 四角形と境界双対面で外向きが食い違うメッシュは例外 (黙って片方を信じない)",
      "食い違う" in (raises(lambda: EX.thrust_cv(tmp / "bell_flip.h5", tmp / "res_100.h5", 0.0, F_id), ValueError) or ""))
write_mesh_h5(tmp / "bell_none.h5", coords, quads, bedges, viz=False, planes=False)
check("(1) 向きを決める情報が無いメッシュは例外",
      "決められない" in (raises(lambda: EX.thrust_cv(tmp / "bell_none.h5", tmp / "res_100.h5", 0.0, F_id), ValueError) or ""))

# ---- (2) 直管の一様流: 壁会計 = 出口面会計、ṁ_in = ṁ_out = ρ u A
pc, pq, pb = pipe_mesh()
NP = len(pc)
write_mesh_h5(tmp / "pipe.h5", pc, pq, pb)
ro, u, P = 1.3, 250.0, 1.0e5
write_res_h5(tmp / "res_200.h5", NP, ro, u, 0.0, P=P)
pwall = sorted({n for e in pb["wall"] for n in e})
write_wall_dump(tmp / "res_wall_3_200.h5", pc, pwall, 0.0, 0.0)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_200.h5", 0.0, F_id)
A = math.pi * 0.5 ** 2
check("(2) 直管の一様流: 出口面会計 = (ρu² + p)A", math.isclose(o["thrust_exit_N"], (ro * u * u + P) * A, rel_tol=1e-12))
check("(2) 直管の一様流: 壁会計 = 出口面会計 (CV 残差 0)", abs(o["cv_residual"]) < 1e-12, f"{o['cv_residual']:.2e}")
check("(2) ṁ_in = ṁ_out = ρuA", math.isclose(o["mdot_in_kg_s"], ro * u * A, rel_tol=1e-12)
      and math.isclose(o["mdot_exit_kg_s"], ro * u * A, rel_tol=1e-12) and abs(o["mass_imbalance"]) < 1e-14)
# 壁の摩擦: node の twall_x = τ0 > 0 (流れ +x に引かれる壁の力) → 流体に働く traction は −τ0 → 推力を τ0·A_wall 下げる
tau0 = 37.0
write_wall_dump(tmp / "res_wall_3_200.h5", pc, pwall, tau0, 0.0)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_200.h5", 0.0, F_id)
A_wall = 2 * math.pi * 0.5 * 2.0
check("(2) node の twall: 摩擦の推力寄与 = −τ0·A_wall (流体に働く traction = −twall)",
      math.isclose(o["friction_N"], -tau0 * A_wall, rel_tol=1e-12), f"{o['friction_N']:.6g} vs {-tau0 * A_wall:.6g}")
check("(2) 摩擦で壁会計だけが下がる (出口面会計は不変)", math.isclose(o["thrust_wall_N"], (ro * u * u + P) * A - tau0 * A_wall, rel_tol=1e-12)
      and math.isclose(o["thrust_exit_N"], (ro * u * u + P) * A, rel_tol=1e-12))
check("(2) 内壁の法線速度 0 (no-slip の場で u·n を積まない根拠の検査)", o["groups"]["wall"]["max_abs_u_n"] == 0.0)
# 規約の関数: node は −twall、cell は +twall
dof_xyz = pc[:, :2]
tr_n, conv_n = EX._wall_traction_on_fluid(tmp / "res_wall_3_200.h5", dof_xyz, np.array(pwall), True)
check("(2) _wall_traction_on_fluid: node は −twall", np.allclose(tr_n[:, 0], -tau0) and "node" in conv_n)
# cell の壁ダンプ (面値・CONNE で面中心) を DOF = 面中心に対応づけ +twall
with h5py.File(tmp / "cellwall.h5", "w") as f:
    f["MESH/COORD"] = pc.ravel()
    ed = np.array(pb["wall"])
    f["MESH/CONNE"] = np.column_stack([np.full(len(ed), 2), np.full(len(ed), 2), ed]).ravel()
    f["VALUE/twall_x"] = np.full(len(ed), -tau0)
    f["VALUE/twall_y"] = np.zeros(len(ed))
cdof = np.vstack([0.5 * (pc[ed[:, 0], :2] + pc[ed[:, 1], :2])])
tr_c, conv_c = EX._wall_traction_on_fluid(tmp / "cellwall.h5", cdof, np.arange(len(ed)), False)
check("(2) _wall_traction_on_fluid: cell は +twall (残差に足す粘性流束 = 流体に働く)", np.allclose(tr_c[:, 0], -tau0) and "cell" in conv_c)
check("(2) 壁ダンプの点を DOF に対応づけられないなら例外",
      raises(lambda: EX._wall_traction_on_fluid(tmp / "res_wall_3_200.h5", dof_xyz + 0.01, np.array(pwall), True)) is not None)
# 壁ダンプ無し・勾配無しは閉じない扱い
(tmp / "res_wall_3_200.h5").unlink()
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_200.h5", 0.0, F_id)
check("(2) 壁ダンプ無し: 摩擦は None・cv_closable False", o["friction_N"] is None and not o["cv_closable"])
write_res_h5(tmp / "res_300.h5", NP, ro, u, 0.0, P=P, grads=False)
write_wall_dump(tmp / "res_wall_3_300.h5", pc, pwall, 0.0, 0.0)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_300.h5", 0.0, F_id)
check("(2) 勾配が res に無い: open_viscous missing・cv_closable False", o["open_viscous"] == "missing" and not o["cv_closable"])
# 開境界の粘性応力: 出口で ∂u/∂x = s 一様 → (τn)_x = μ(2s − 2/3 s) = 4/3 μ s → 出口面会計を 4/3 μ s A_e 下げる
with h5py.File(tmp / "res_200.h5", "r+") as f:
    f["VALUE/dUxdx"][:] = 1.0e4
    f["VALUE/axisym_divU"][:] = 1.0e4
    f["VALUE/vis_lam"][:] = 2.0e-4
write_wall_dump(tmp / "res_wall_3_200.h5", pc, pwall, 0.0, 0.0)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_200.h5", 0.0, F_id)
tv = 4.0 / 3.0 * 2.0e-4 * 1.0e4
check("(2) 出口面の粘性応力 τ_xx = μ(2∂u/∂x − 2/3 ∇·u) を出口面会計から引く",
      math.isclose(o["groups"]["outlet"]["viscous_tau_n_x"], tv * A, rel_tol=1e-12)
      and math.isclose(o["thrust_exit_N"], (ro * u * u + P) * A - tv * A, rel_tol=1e-12))
check("(2) 入口面の粘性応力 (外向き −x) も同じ式で壁会計に入り、一様勾配なら CV 残差 0", abs(o["cv_residual"]) < 1e-12)

# ---- (3) P の出所
gas = p.exhaust.frozen_gas()
T, rho, uu = 2000.0, 0.8, 300.0
roe = rho * (float(gas.e_sens(np.array([T]))[0]) + 0.5 * uu * uu)
P_tp = EX.tp_pressure_from_conserved(np.array([rho]), np.array([roe]), np.array([0.5 * rho * uu * uu]), gas)[0]
check("(3) TP の P: e_sens(T) の逆算で ρRT (相対 1e-6)", math.isclose(P_tp, rho * gas.R * T, rel_tol=1e-6), f"{P_tp:.6g}")
write_res_h5(tmp / "res_400.h5", NP, rho, uu, 0.0, P=None, roe=roe)
with h5py.File(tmp / "res_400.h5", "r") as f:
    Pv, src = EX._pressure(f["VALUE"], np.arange(3), gas=gas)
    check("(3) VALUE/P が無ければ TP は e_sens の逆算", "TP" in src and np.allclose(Pv, rho * gas.R * T, rtol=1e-6))
    Pc_, srcc = EX._pressure(f["VALUE"], np.arange(3), gamma=1.4)
    check("(3) CPG は (γ−1)(ρE − ρu²/2)", "CPG" in srcc and np.allclose(Pc_, 0.4 * (roe - 0.5 * rho * uu * uu)))
    check("(3) 気体も γ も無ければ例外", raises(lambda: EX._pressure(f["VALUE"], np.arange(3))) is not None)


# ---- (4) Euler (slip 壁) の会計 (plan §5.1 #6 (W1)): wall_model="slip" は内壁の traction 0・壁ダンプを読まない・μ_eff = 0 を確認
def write_slip_dump(path, coords, wall_nodes):
    """forge の slip 境界のダンプと同じく twall_* が無い (case/45 run_0001_euler_shortest_dry/res_wall_3_12000.h5 の VALUE と同じ並び)。"""
    wn = np.array(wall_nodes)
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords[wn].astype(np.float64).ravel()
        for k in ("Ps", "Pt", "Ts", "Tt", "Ux", "Uy", "Uz", "ro", "roUx", "roUy", "roUz", "roe"):
            f[f"VALUE/{k}"] = np.zeros(len(wn))


def write_euler_res(path, n, ro, ux, ur, P, grad=0.0, grads=True):
    """Euler の res: vis_lam = vis_turb = 0 (forge の viscMethod 0 は vis_lam = visc = 0、乱流なし)。勾配は `grad` (非零でも粘性項は 0)。"""
    write_res_h5(path, n, ro, ux, ur, P=P, grads=grads)
    with h5py.File(path, "r+") as f:
        if grads:
            f["VALUE/vis_lam"][:] = 0.0
            for g in ("dUxdx", "dUxdy", "dUydx", "dUydy", "axisym_divU"):
                f[f"VALUE/{g}"][:] = grad
        else:
            f["VALUE/vis_lam"] = np.zeros(n)
            f["VALUE/vis_turb"] = np.zeros(n)


write_euler_res(tmp / "res_500.h5", NP, ro, u, 0.0, P, grad=1.0e4)
write_slip_dump(tmp / "res_wall_3_500.h5", pc, pwall)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_500.h5", 0.0, F_id, wall_model="slip")
check("(4) Euler (slip): 直管の一様流で壁会計 = 出口面会計 (CV 残差 0)・出口面会計 = (ρu² + p)A",
      abs(o["cv_residual"]) < 1e-12 and math.isclose(o["thrust_exit_N"], (ro * u * u + P) * A, rel_tol=1e-12), f"{o['cv_residual']:.1e}")
check("(4) Euler (slip): 摩擦 0・cv_closable True・使った経路 (wall_model slip、規約の文字列) を戻り値に残す",
      o["friction_N"] == 0.0 and o["cv_closable"] and o["wall_model"] == "slip" and "slip" in o["wall_traction_convention"]
      and "wall_dump" not in o, o["wall_traction_convention"])
check("(4) Euler: 勾配が非零 (1e4) でも開境界の粘性項は 0 (viscMethod 0 の vis_lam = visc = 0)",
      o["groups"]["outlet"]["viscous_tau_n_x"] == 0.0 and o["groups"]["inlet"]["viscous_tau_n_x"] == 0.0 and o["open_viscous"] == "evaluated")
m = raises(lambda: EX.thrust_cv(tmp / "pipe.h5", tmp / "res_500.h5", 0.0, F_id), ValueError) or ""
check("(4) NS (既定 wall_model viscous) で壁ダンプに twall_x が無ければ止める (摩擦を黙って 0 にしない)", "twall" in m, m[:90])
m = raises(lambda: EX.thrust_cv(tmp / "pipe.h5", tmp / "res_500.h5", 0.0, F_id, wall_dump=tmp / "res_wall_3_500.h5"), ValueError) or ""
check("(4) NS で壁ダンプを明示しても twall_x が無ければ止める", "twall" in m, m[:90])
m = raises(lambda: EX.thrust_cv(tmp / "pipe.h5", tmp / "res_200.h5", 0.0, F_id, wall_model="slip"), ValueError) or ""
check("(4) 粘性の場 (vis_lam 2e-4) に slip を使うと止める (μ_eff ≠ 0 — 壁の摩擦を黙って 0 にしない)", "μ_eff" in m, m[:90])
write_res_h5(tmp / "res_600.h5", NP, ro, u, 0.0, P=P, grads=False)
m = raises(lambda: EX.thrust_cv(tmp / "pipe.h5", tmp / "res_600.h5", 0.0, F_id, wall_model="slip"), ValueError) or ""
check("(4) vis_lam・vis_turb の無い res に slip を使うと止める (μ_eff = 0 を確かめられない)", "vis_lam" in m, m[:90])
write_euler_res(tmp / "res_700.h5", NP, ro, u, 0.0, P, grads=False)
o = EX.thrust_cv(tmp / "pipe.h5", tmp / "res_700.h5", 0.0, F_id, wall_model="slip")
check("(4) Euler で勾配が res に無くても μ_eff = 0 なので閉じる会計 (open_viscous に理由)",
      o["cv_closable"] and o["open_viscous"].startswith("zero") and abs(o["cv_residual"]) < 1e-12, o["open_viscous"])
check("(4) 未知の wall_model は止める",
      raises(lambda: EX.thrust_cv(tmp / "pipe.h5", tmp / "res_500.h5", 0.0, F_id, wall_model="inviscid"), ValueError) is not None)
write_euler_res(tmp / "res_800.h5", N, 1.0, 0.0, 0.0, p0)
for pa in (0.0, p0 / 3.0):
    o = EX.thrust_cv(tmp / "bell.h5", tmp / "res_800.h5", pa, F_id, wall_model="slip")
    check(f"(4) Euler (slip): TOP ベルの静止気体 p_a={pa:g} で CV 残差 0・出口面会計 = (p−p_a)A_e",
          abs(o["cv_residual"]) < 1e-12 and math.isclose(o["thrust_exit_N"], (p0 - pa) * math.pi * r_e ** 2, rel_tol=1e-12))

# ---------------------------------------------------------------- ゲート (i): 合成の残差履歴 + StageManifest
def make_run(d: Path, main_hist, pre_hist=None, manifest=True):
    d.mkdir(parents=True)
    cfg_soft = runner._solver_config_tp(p, 2000, 2000, "soft")
    cfg_main = runner._solver_config_tp(p, len(main_hist), max(len(main_hist) // 12, 1), "main")
    bc = runner._bcond_config_tp(p)
    cols = ["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]

    def wr(name, hist):
        with open(d / name, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(cols)
            for s, v in enumerate(hist):
                w.writerow([s, 4, "outer_end", v, 1e3 * v, 5e2 * v, 0.0, 1e7 * v, 10 * v, 1e8 * v])
    sm = StageManifest(d)
    if pre_hist is not None:
        wr("residual_history_S1_soft.csv", pre_hist)
        sm.add("S1_soft", cfg_soft, bc, history="residual_history_S1_soft.csv")
    wr("residual_history_S3_main.csv", main_hist)
    sm.add("S3_main", cfg_main, bc, history="residual_history_S3_main.csv")
    if manifest:
        sm.write()
    return d


n = 3000
s = np.arange(n)
falling = list(3e-4 * 10 ** (-4.0 * s / n))                                           # 4 桁落ちる
stalled = list(3e-4 * np.maximum(10 ** (-1.5 * s / 1000), 10 ** -1.5) * (1 + 0.05 * np.sin(s / 7.0)))   # 1.5 桁で横ばい
pre = list(3e-2 * 10 ** (-2.0 * s / n))
g = DR.gate_convergence(make_run(tmp / "conv_pass", falling, pre))
check("(i) 本段で 4 桁落ちる履歴は PASS (前段の別区間を基準にしない)", g["ok"] and g["class"] == "PASS", g["verdict"])
check("(i) 判定区間は本段 (S3_main) だけ", g["segment"] is not None and "S3_main" in g["segment"] and "S1_soft" not in g["segment"],
      g["segment"])
g = DR.gate_convergence(make_run(tmp / "conv_stall", stalled, pre))
check("(i) STALLED (1.5 桁で横ばい) は記録して通す (FAIL は DIVERGED と有意な末尾上昇だけ; §4.1 (i))",
      g["ok"] and g["class"] == "STALLED" and "併記" in g["note"], g["verdict"])
slow = list(3e-4 * 10 ** (-1.5 * s / n))                                              # 1.5 桁、末尾まで低下中
g = DR.gate_convergence(make_run(tmp / "conv_slow", slow, pre))
check("(i) まだ低下中 (still converging) も通す", g["ok"] and g["class"] == "CONVERGING", g["verdict"])
rise = list(np.where(s < 0.6 * n, 3e-4 * 10 ** (-3.0 * s / (0.6 * n)), 3e-7 * 10 ** (1.0 * (s - 0.6 * n) / (0.4 * n))))
g = DR.gate_convergence(make_run(tmp / "conv_rise", rise, pre))
check("(i) 有意な末尾上昇 (3 桁落ちてから 1 桁戻る) は RISING で落とす", (not g["ok"]) and g["class"] == "RISING"
      and "rms_ro" in g["rising"], f"{g['verdict']} rising={g['rising']}")
g = DR.gate_convergence(make_run(tmp / "conv_nomani", falling, pre, manifest=False))
check("(i) stage_manifest.json が無ければ判定不能 = 不合格", (not g["ok"]) and g["class"] == "UNDECIDABLE", g["verdict"])
nanh = falling[:-5] + [float("nan")] * 5
g = DR.gate_convergence(make_run(tmp / "conv_nan", nanh, pre))
check("(i) NaN は DIVERGED で落とす", (not g["ok"]) and g["class"] == "DIVERGED", g["verdict"])


# ---------------------------------------------------------------- ゲート (ii): 推力の時系列
def write_series(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(("step",) + runner.THRUST_HISTORY_COLS)
        for r in rows:
            w.writerow([r["step"]] + [repr(float(r.get(c, 1.0))) for c in runner.THRUST_HISTORY_COLS])


def qs_case(name, rows):
    (tmp / name).mkdir()
    write_series(tmp / name / "thrust_history.csv", rows)
    return DR.gate_quasisteady(tmp / name / "thrust_history.csv")


steps = np.arange(1, 21) * 3000                                                     # 本段の 1/20 ごと (末尾窓 8 点)
base = dict(thrust_exit_N=1.135e6, thrust_wall_N=1.135e6, eta=0.9792, eta_wall=0.9792, mdot_exit_kg_s=268.8,
            mdot_in_kg_s=268.8, mdot_ratio_1d=0.9925, cv_residual=4.0e-5, friction_N=-1.18e4)
steady = [{"step": int(st), **{k: v * (1 + 1e-6 * math.sin(i)) for k, v in base.items()}} for i, st in enumerate(steps)]
gq = qs_case("qs_ok", steady)
check("(ii) 平坦な系列 (20 点) は ALL STEADY で通る (推力系 4 列 + ṁ)", gq["ok"] and set(gq["per_column"].values()) == {"STEADY"}
      and set(gq["per_column"]) == set(runner.THRUST_QS_COLS + runner.THRUST_QS_MDOT_COLS) and gq["n_tail"] == 8, gq["overall"])
check("(ii) cv_residual は列にしない (0 を横切る; iv で見る)", "cv_residual" not in gq["per_column"])
gq = qs_case("qs_short", steady[:12])
check("(ii) 末尾窓 8 点未満 (12 点 → 5 点) は落とす (本段の出力間隔は本段長の 1/20 以下)", (not gq["ok"]) and gq["n_tail"] == 5
      and any("末尾窓" in x for x in gq["problems"]), gq["overall"])
drift = [{"step": int(st), **base, "thrust_exit_N": base["thrust_exit_N"] * (1 + 1e-3 * i), "eta": base["eta"] * (1 + 1e-3 * i)}
         for i, st in enumerate(steps)]
gq = qs_case("qs_drift", drift)
check("(ii) DRIFTING (推力が単調に 1e-3/点 動く) を落とす", (not gq["ok"]) and gq["per_column"].get("thrust_exit_N") == "DRIFTING",
      str(gq["per_column"]))
ew = [{"step": int(st), **base, "eta_wall": base["eta_wall"] * (1 + 2e-4 * (-1) ** i)} for i, st in enumerate(steps)]
gq = qs_case("qs_etawall", ew)
check("(ii) η (壁会計) の振動 ±2e-4 (> 1.5e-4) も落とす (壁会計の η も列)", (not gq["ok"]) and gq["per_column"].get("eta_wall") == "OSCILLATING",
      str(gq["per_column"]))
cvd = [{"step": int(st), **base, "cv_residual": 4.0e-5 * (1 + 0.5 * i)} for i, st in enumerate(steps)]
gq = qs_case("qs_cv", cvd)
check("(ii) CV 残差の時系列が動いても (ii) は落とさない (列でない)", gq["ok"], gq["overall"])
mo = [{"step": int(st), **base, "mdot_exit_kg_s": base["mdot_exit_kg_s"] * (1 + 4e-4 * math.sin(1.3 * i))} for i, st in enumerate(steps)]
gq = qs_case("qs_mdot_ok", mo)
check("(ii) ṁ の振動 ±4e-4 は ṁ の閾値 1e-3 で通る", gq["ok"] and gq["per_column"].get("mdot_exit_kg_s") == "STEADY", gq["overall"])
mo2 = [{"step": int(st), **base, "thrust_exit_N": base["thrust_exit_N"] * (1 + 4e-4 * math.sin(1.3 * i))} for i, st in enumerate(steps)]
gq = qs_case("qs_thrust_osc", mo2)
check("(ii) 同じ ±4e-4 の振動が推力にあれば 1.5e-4 で落とす", (not gq["ok"]) and gq["per_column"].get("thrust_exit_N") == "OSCILLATING",
      str(gq["per_column"]))
# 単調で増分がゆっくり減衰 (r 0.98): 末尾窓のドリフトは 6e-5 (< 実効の単調許容 7.5e-5) で STEADY だが、漸近値は最終値から 3.3e-4
asy = [{"step": int(st), **base, "eta": base["eta"] * (1 - 5e-4 * 0.98 ** i)} for i, st in enumerate(steps)]
gq = qs_case("qs_asym", asy)
check("(ii) 単調列は漸近値の最終比も 1.5e-4 以内 (STEADY でも漸近値が離れていれば落とす)", (not gq["ok"])
      and gq["per_column"].get("eta") == "STEADY" and not gq["asymptote"]["eta"]["ok"]
      and abs(gq["asymptote"]["eta"]["final_ratio"]) > 1.5e-4, str(gq["asymptote"].get("eta")))
nanr = [dict(r) for r in steady]
nanr[-1]["thrust_wall_N"] = float("nan")
gq = qs_case("qs_nan", nanr)
check("(ii) 非有限値は NONFINITE で落とす", (not gq["ok"]) and gq["per_column"].get("thrust_wall_N") == "NONFINITE",
      str(gq["per_column"]))
(tmp / "qs_cols").mkdir()
with open(tmp / "qs_cols" / "thrust_history.csv", "w") as fh:
    fh.write("step,thrust_exit_N\n" + "".join(f"{int(st)},1.0\n" for st in steps))
gq = DR.gate_quasisteady(tmp / "qs_cols" / "thrust_history.csv")
check("(ii) 列が欠けていれば判定不能 = 不合格", not gq["ok"], gq["overall"])
gq = DR.gate_quasisteady(tmp / "qs_ok" / "thrust_history.csv", qs={"drift": 1e-9, "osc": 1e-9})
check("(ii) evaluate.quasisteady の閾値を渡せる (極端に厳しくすれば落ちる)", not gq["ok"])

# ---------------------------------------------------------------- ゲート (iii)(iv)(v) と総合の分類
m_ok = {"eta_cf": 0.979, "mdot_ratio_1d": 0.9925, "cv_residual": 4e-5, "cv_closable": True}
check("(iii) 物理範囲内は通る", DR.gate_physics(m_ok)["ok"])
check("(iii) ṁ/ṁ_1D 0.90 (不始動の偽アトラクタ) を落とす", not DR.gate_physics({**m_ok, "mdot_ratio_1d": 0.90})["ok"])
check("(iii) η 1.05 を落とす", not DR.gate_physics({**m_ok, "eta_cf": 1.05})["ok"])
check("(iii) 欠損・非有限は落とす", not DR.gate_physics({**m_ok, "eta_cf": None})["ok"]
      and not DR.gate_physics({**m_ok, "mdot_ratio_1d": float("nan")})["ok"])
check("(iv) 最終 |CV 残差| 4e-5 は通る", DR.gate_cv(m_ok)["ok"] and DR.gate_cv(m_ok)["tol"] == 1e-4)
check("(iv) |CV 残差| 2e-4 (> 1e-4) を落とす (旧 1e-3 では通っていた)", not DR.gate_cv({**m_ok, "cv_residual": -2e-4})["ok"])
check("(iv) run_0002 の最終値 −6.6e-4 は落とす", not DR.gate_cv({**m_ok, "cv_residual": -6.6e-4})["ok"])
check("(iv) ちょうど 1e-4 は通る (≤)", DR.gate_cv({**m_ok, "cv_residual": 1e-4})["ok"])
check("(iv) 閉じない会計 (cv_closable False) は落とす", not DR.gate_cv({**m_ok, "cv_closable": False})["ok"])
check("(iv) NaN は落とす", not DR.gate_cv({**m_ok, "cv_residual": float("nan")})["ok"])
(tmp / "mq").mkdir()
(tmp / "mq" / "MESH_QUALITY.txt").write_text("(primal/cell 変換で検査)\nmax AR 836\nVERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
check("(v) MESH_QUALITY の PASS を読む", DR.gate_mesh_quality(tmp / "mq")["ok"])
(tmp / "mq" / "MESH_QUALITY.txt").write_text("VERDICT: FAIL (AR 1200 > 1000)\n")
check("(v) MESH_QUALITY の FAIL を落とす", not DR.gate_mesh_quality(tmp / "mq")["ok"])
check("(v) MESH_QUALITY.txt が無ければ落とす", not DR.gate_mesh_quality(tmp / "nonexistent")["ok"])
(tmp / "mq" / "MESH_QUALITY.txt").write_text("max AR 2700\nVERDICT: PASS (AR<= 5000, skew<= 0.90)\nAR_RELAX: PASS (...)\n")
gm = DR.gate_mesh_quality(tmp / "mq")
check("(v) AR ≤ 5000 の格子で MESH_AR_RELAX.json が無ければ判定不能 = 不合格", (not gm["ok"]) and gm["ar_limit"] == 5000.0
      and "判定不能" in str(gm.get("ar_relax")), str(gm))
(tmp / "mq" / "MESH_AR_RELAX.json").write_text('{"ok": false, "n_high_ar": 10, "max_skew_high_ar": 0.3, "n_violations": 4, '
                                              '"inlet_planar": true, "outlet_planar": true}')
check("(v) 例外の条件 (AR > 1000 のスキュー ≤ 0.1) を満たさない格子は落とす", not DR.gate_mesh_quality(tmp / "mq")["ok"])
(tmp / "mq" / "MESH_AR_RELAX.json").write_text('{"ok": true, "n_high_ar": 10, "max_skew_high_ar": 0.004, "n_violations": 0, '
                                              '"inlet_planar": true, "outlet_planar": true}')
check("(v) 例外の条件を満たせば AR ≤ 5000 の格子を通す", DR.gate_mesh_quality(tmp / "mq")["ok"])
(tmp / "wr").mkdir()
gw = DR.gate_wall_resolution(tmp / "wr")
check("(v) 壁ダンプの無い run の y1+ は判定不能 = 不合格", not gw["ok"], gw["verdict"])
# (v-b) は受け渡しの関門と同じ規則 (plan tooling-nozzle-isothermal-wall-chain §5.1 #36s): bcond の全 no-slip 壁・停止の step
# (既定は最終の主スナップショット)・面積の重み。合成の run は test_gate_bad_input.make_wall_run (2 枚の no-slip 壁 + 出口の出力)
import test_gate_bad_input as TGB  # noqa: E402

_WC = TGB.wall_cases()
_both = {"wall_a": [1000, 2000], "wall_b": [1000, 2000]}


def _wall_run(name, walls, steps, snaps=(1000, 2000)):
    d = Path(TGB.make_wall_run(str(tmp / name), walls, steps))
    for s_ in snaps:                                   # 主スナップショット (vis_lam の無い空の HDF5)
        h5py.File(d / f"res_{s_}.h5", "w").close()
    return d


gw = DR.gate_wall_resolution(_wall_run("wr_ctrl", _WC["control"], _both))
check("(v-b) 対照: 両壁 (wall_a・wall_b) とも最終の step 2000 に出力 → PASS (群 wall に固定しない・面積)",
      gw["ok"] and gw["step"] == 2000 and sorted(gw["walls"]) == ["wall_a", "wall_b"] and gw["weight"] == "area", str(gw))
gw = DR.gate_wall_resolution(_wall_run("wr_b_old", _WC["control"], {"wall_a": [1000, 2000], "wall_b": [1000]}))
check("(v-b) wall_b の最終の step の出力だけ欠く → 判定不能 = 不合格 (最新の step に戻らない)",
      not gw["ok"] and gw["kind"] == "UNDETERMINED" and "wall_b" in (gw["why"] or ""), gw["verdict"])
gw = DR.gate_wall_resolution(_wall_run("wr_area", _WC["inv_area_fail"], _both))
check("(v-b) 点数の割合 1/102 では通るが面積の割合 約 48 % の入力 → FAIL (面積で判定)", not gw["ok"] and gw["kind"] == "FAIL", gw["verdict"])
gw = DR.gate_wall_resolution(_wall_run("wr_point", _WC["inv_point_fail"], _both))
check("(v-b) 点数の割合 10/31 では落ちるが面積の割合 約 0.5 % の入力 → PASS", gw["ok"], gw["verdict"])
gw = DR.gate_wall_resolution(_wall_run("wr_step", _WC["control"], _both), step=1000)
check("(v-b) step を渡せばその step (推力の会計の res_file) で判定", gw["ok"] and gw["step"] == 1000, str(gw.get("step")))

okg = {"ok": True}
allok = {k: okg for k in ("convergence", "quasisteady", "physics", "cv", "mesh_quality", "wall_resolution")}
check("総合: 全ゲート通過で PASS", DR.combine_tp_gates(allok)[:2] == ("PASS", None))
check("総合: STALLED は RETRYABLE", DR.combine_tp_gates({**allok, "convergence": {"ok": False}})[1] == "RETRYABLE")
check("総合: DRIFTING は RETRYABLE", DR.combine_tp_gates({**allok, "quasisteady": {"ok": False}})[1] == "RETRYABLE")
check("総合: CV 残差超過は SUSPECT", DR.combine_tp_gates({**allok, "cv": {"ok": False}})[1] == "SUSPECT")
check("総合: y1+ 不合格は SUSPECT", DR.combine_tp_gates({**allok, "wall_resolution": {"ok": False}})[1] == "SUSPECT")
phf = DR.gate_physics({**m_ok, "mdot_ratio_1d": 0.90})
check("総合: 収束・定常で物理ゲート外は INFEASIBLE", DR.combine_tp_gates({**allok, "physics": phf})[1] == "INFEASIBLE")
phc = DR.gate_physics({**m_ok, "mdot_ratio_1d": 1.01})
check("総合: Cd > 1 は SUSPECT", DR.combine_tp_gates({**allok, "physics": phc})[1] == "SUSPECT")
check("総合: ゲートの欠落は不合格", DR.combine_tp_gates({k: v for k, v in allok.items() if k != "cv"})[0] == "FAIL")

# ---------------------------------------------------------------- 引き継ぎ (b): 評価ごとの YAML の相対パス
raw = yaml.safe_load(PROBLEM_TP.read_text())
ab = DR.absolutize_problem_paths(raw, PROBLEM_TP.parent)
check("(b) gas.exhaust の相対パスを問題 YAML の所在で絶対パスに", Path(ab["gas"]["exhaust"]).is_absolute()
      and Path(ab["gas"]["exhaust"]).exists() and raw["gas"]["exhaust"] == "exhaust_lox_lh2_pc100_of6.yaml")
check("(b) 絶対パスはそのまま", DR.absolutize_problem_paths(ab, "/nonexistent")["gas"]["exhaust"] == ab["gas"]["exhaust"])
check("(b) gas.species_db も直す", DR.absolutize_problem_paths({"gas": {"species_db": "db.yaml"}}, "/a/b")["gas"]["species_db"] == "/a/b/db.yaml")
check("(b) 空・無しは触らない", DR.absolutize_problem_paths({"gas": {"exhaust": ""}}, "/a")["gas"]["exhaust"] == ""
      and DR.absolutize_problem_paths({"spec": {}}, "/a") == {"spec": {}})
check("(b) 文字列でない値は例外", raises(lambda: DR.absolutize_problem_paths({"gas": {"exhaust": 3}}, "/a"), ValueError) is not None)
cdir = tmp / "campaign"
camp = DR.Campaign(PROBLEM_TP, cdir, None)
check("(b) TP の Campaign: --seed-res なしで作れる・1D 理想値が分母", camp.tp and camp.seed_res is None
      and math.isclose(camp.mdot_1d, runner.ideal_1d_tp(p)["mdot"]) and math.isclose(camp.F_1d, runner.ideal_1d_tp(p)["F"]))
out_yaml = camp._write_problem([30.0, 9.5, 17.0], cdir / "e_000.yaml")
p2 = load_problem(out_yaml)
check("(b) 書き出した評価 YAML を別ディレクトリから読める (排気ファイルが解決する)", p2.exhaust is not None
      and math.isclose(p2.exhaust.Tc, p.exhaust.Tc) and p2.dv["theta_a_deg"]["value"] == 30.0)
check("(b) CPG の Campaign は --seed-res 必須", raises(lambda: DR.Campaign(ROOT / "case" / "40.nozzle_design_tool" / "problem_bell_top.yaml",
                                                                    tmp / "camp_cpg", None), ValueError) is not None)

# ---------------------------------------------------------------- runner: 段の計画・TP の config・格子のキー
plan = runner.tp_stage_plan(p)
check("runner: 既定の段 soft → mid → main (main は evaluate.nStepOuter・前段は最終 step だけ書く)",
      [t[1] for t in plan] == ["soft", "mid", "main"] and plan[0][2] == 2000 and plan[0][3] == 2000
      and plan[-1][2] == int(p.evaluate["nStepOuter"]) and plan[-1][3] == int(p.evaluate["outStepInterval"]), str(plan))
p_bad = load_problem(PROBLEM_TP)
p_bad.evaluate = {**p.evaluate, "nStepOuter": 1000, "outStepInterval": 300}
check("runner: 本段の step 数が outStepInterval の倍数でなければ例外 (黙って直さない)", raises(lambda: runner.tp_stage_plan(p_bad), ValueError) is not None)
check("runner: main で終わらない段は例外", raises(lambda: runner.tp_stage_plan(p, ("soft", "mid")), ValueError) is not None)
check("runner: 未知の段は例外", raises(lambda: runner.tp_stage_plan(p, ("soft", "fast", "main")), ValueError) is not None)
check("runner: main の重複は例外", raises(lambda: runner.tp_stage_plan(p, ("main", "main")), ValueError) is not None)
p_bad.evaluate = {**p.evaluate, "stage_steps": {"soft": 0}}
check("runner: 0 step の段は例外", raises(lambda: runner.tp_stage_plan(p_bad), ValueError) is not None)
p_bad.evaluate = {**p.evaluate, "stage_steps": {"mid": 500}}
check("runner: stage_steps で前段の step 数を上書き", runner.tp_stage_plan(p_bad)[1][2] == 500)
p_bad.evaluate = {k: v for k, v in p.evaluate.items() if k != "outStepInterval"}
p_bad.evaluate["nStepOuter"] = 60000
check("runner: 本段の出力間隔の既定は本段長の 1/20 (準定常ゲートの末尾窓 8 点)", runner.tp_stage_plan(p_bad)[-1][3] == 3000)
cfg = runner._solver_config_tp(p, 12000, 1000, "main")
doc = yaml.safe_load(cfg)
check("runner: TP の config は推力の粘性応力用の勾配を extraFields で出す (level 1 のまま)",
      doc["output"]["level"] == 1 and list(doc["output"]["extraFields"]) == list(runner.TP_EXTRA_FIELDS))
check("runner: TP の config に pRef と bndFirstOrder を書かない (引き継ぎ (d) は上位に諮るまで現状)",
      "pRef" not in cfg and "bndFirstOrder" not in cfg)
mp0 = runner.mesh_params(p, 0.14)
check("runner: スロート細分キーを書いていない問題は従来の 4 キーだけの格子 (既定値のまま)",
      mp0.wall_first_frac_throat is None and mp0.throat_width == Mesh2DParams().throat_width
      and mp0.wall_first_blend_x1 == Mesh2DParams().wall_first_blend_x1)
p_m = load_problem(PROBLEM_TP)
p_m.mesh = {**p.mesh, "wall_first_frac_throat": 3.0e-6, "wall_first_blend_x1": 8.0, "wall_first_up_x0": -3.0, "wall_first_up_x1": -1.0}
mp1 = runner.mesh_params(p_m, 0.14)
check("runner: wall_first_frac_throat・ブレンドのキーを読む (引き継ぎ (c))", mp1.wall_first_frac_throat == 3.0e-6
      and mp1.wall_first_blend_x1 == 8.0 and mp1.wall_first_up_x0 == -3.0 and mp1.wall_first_up_x1 == -1.0)

p_m.mesh = {**p.mesh, "wall_normal_layer": [0.02, 0.8], "wall_normal_end_stations": 20}
mp2 = runner.mesh_params(p_m, 0.14)
check("runner: wall_normal_layer (組)・wall_normal_end_stations (整数) を渡す (#2d)", mp2.wall_normal_layer == (0.02, 0.8)
      and mp2.wall_normal_end_stations == 20 and isinstance(mp2.wall_normal_end_stations, int))
p_m.mesh = {**p.mesh, "wall_normal_layer": 0.02}
check("runner: wall_normal_layer が組でなければ例外", raises(lambda: runner.mesh_params(p_m, 0.14), ValueError) is not None)
p_m.mesh = {**p.mesh, "wall_normal_layer": [0.02, 0.8], "wall_normal_end_stations": 2.5}
check("runner: wall_normal_end_stations が整数でなければ例外", raises(lambda: runner.mesh_params(p_m, 0.14), ValueError) is not None)


# ---------------------------------------------------------------- 壁法線の近壁層 (#2d): mesh2d の端の断面・AR ≤ 5000 の例外の条件
def write_qc_h5(path, coords, quads):
    """primal の四角形 (VIZMESH/CONNE、XDMF の四角形 = 5) と座標だけの h5 (mesh_ar_relax_check が読む分)。"""
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float64).ravel()
        f["VIZMESH/CONNE"] = np.column_stack([np.full(len(quads), 5), quads]).astype(np.int32).ravel()
        f["VIZMESH"].attrs["nVizCells"] = len(quads)


wn_base = dict(ni=241, nj=61, wall_first_frac=2.5e-5, wall_first_frac_throat=3.5e-6, wall_first_blend_x0=0.3,
               wall_first_blend_x1=14.0, throat_refine=8.0, throat_width=3.0, scale=0.14)
c_r, q_r, b_r = generate_axisym_mesh(wall, Mesh2DParams(**wn_base))
c_n, q_n, b_n = generate_axisym_mesh(wall, Mesh2DParams(**wn_base, wall_normal_layer=(0.02, 0.8)))
c_e, q_e, b_e = generate_axisym_mesh(wall, Mesh2DParams(**wn_base, wall_normal_layer=(0.02, 0.8), wall_normal_end_stations=8))
oute = lambda c, b: float(np.ptp(c[sorted({n for e in b["outlet"] for n in e}), 0]))  # noqa: E731
check("mesh2d: 法線の層だけだと θ_e ≠ 0 の出口面は平面でない (節点が下流へずれる)", oute(c_n, b_n) > 1e-3, f"{oute(c_n, b_n):.3e} m")
check("mesh2d: wall_normal_end_stations で入口・出口面が平面 (x 一定) に戻る", oute(c_e, b_e) == 0.0
      and float(np.ptp(c_e[sorted({n for e in b_e['inlet'] for n in e}), 0])) == 0.0)
ni_, nj_ = 241, 61
mid = slice(20, ni_ - 20)
check("mesh2d: 端の断面の外 (内部) は法線の層と同じ節点", np.array_equal(c_e.reshape(ni_, nj_, 3)[mid], c_n.reshape(ni_, nj_, 3)[mid]))
check("mesh2d: wall_normal_end_stations だけ (層なし) は例外", raises(lambda: generate_axisym_mesh(
    wall, Mesh2DParams(**wn_base, wall_normal_end_stations=8)), ValueError) is not None)
check("mesh2d: wall_normal_end_stations の範囲外は例外", raises(lambda: generate_axisym_mesh(
    wall, Mesh2DParams(**wn_base, wall_normal_layer=(0.02, 0.8), wall_normal_end_stations=0)), ValueError) is not None)
c_0, _, _ = generate_axisym_mesh(wall, Mesh2DParams(**wn_base, wall_normal_layer=None, wall_normal_end_stations=None))
check("mesh2d: 新しいキーを書かなければ従来の格子とビット同一", np.array_equal(c_0, c_r))
# 例外の条件の合格例は ni 481 (出口近くの AR を 1000 未満にして、半径方向へ戻す端の断面に AR > 1000 のセルを置かない)
c_g, q_g, b_g = generate_axisym_mesh(wall, Mesh2DParams(**{**wn_base, "ni": 481}, wall_normal_layer=(0.02, 0.8),
                                                        wall_normal_end_stations=10))
write_qc_h5(tmp / "qc_r.h5", c_r, q_r)
write_qc_h5(tmp / "qc_e.h5", c_e, q_e)
write_qc_h5(tmp / "qc_n.h5", c_n, q_n)
write_qc_h5(tmp / "qc_g.h5", c_g, q_g)
ck_r = runner.mesh_ar_relax_check(tmp / "qc_r.h5", c_r, b_r, 0.14, 5000.0)
ck_e = runner.mesh_ar_relax_check(tmp / "qc_e.h5", c_e, b_e, 0.14, 5000.0)
ck_n = runner.mesh_ar_relax_check(tmp / "qc_n.h5", c_n, b_n, 0.14, 5000.0)
ck_g = runner.mesh_ar_relax_check(tmp / "qc_g.h5", c_g, b_g, 0.14, 5000.0)
check("AR 緩和の条件: 半径方向の線の格子は AR > 1000 のセルにスキューがあり不合格 (違反を数える)",
      ck_r["n_high_ar"] > 0 and ck_r["n_violations"] > 0 and not ck_r["ok"] and ck_r["required"],
      f"AR>1000 {ck_r['n_high_ar']}・違反 {ck_r['n_violations']}・最大スキュー {ck_r['max_skew_high_ar']:.3f}")
check("AR 緩和の条件: 壁法線の層 + 端の断面で合格 (AR > 1000 のスキュー ≤ 0.1・面は平面)",
      ck_g["ok"] and ck_g["n_high_ar"] > 0 and ck_g["max_skew_high_ar"] <= 0.1 and ck_g["inlet_planar"] and ck_g["outlet_planar"],
      f"AR>1000 {ck_g['n_high_ar']}・最大スキュー {ck_g['max_skew_high_ar']:.4f}")
check("AR 緩和の条件: 半径方向へ戻す端の断面に AR > 1000 のセルがあれば不合格 (ni 241: 出口近くで AR 1400・スキュー 0.101)",
      (not ck_e["ok"]) and ck_e["n_violations"] > 0 and ck_e["outlet_planar"] and ck_e["worst_high_ar"]["x_over_rt"] > 17.0,
      f"違反 {ck_e['n_violations']}・最悪 x {ck_e['worst_high_ar']['x_over_rt']:.2f} スキュー {ck_e['max_skew_high_ar']:.4f}")
check("AR 緩和の条件: 出口面が平面でなければ不合格", (not ck_n["ok"]) and not ck_n["outlet_planar"] and ck_n["n_violations"] == 0)
check("AR 緩和の条件: ar_max ≤ 1000 なら適用なし (required False)",
      not runner.mesh_ar_relax_check(tmp / "qc_g.h5", c_g, b_g, 0.14, 1000.0)["required"])

# ---------------------------------------------------------------- (D) 直接モード (#3f)
from forge_design.geometry.wall import DirectBellWall  # noqa: E402
from forge_design.geometry.wall_bell_c2 import BellC2Spec, cone15_length, design_bell_c2  # noqa: E402

PROBLEM_D = CASE68 / "problem_bell_tp_d_wn.yaml"
pd = load_problem(PROBLEM_D)
wd, des = runner.build_wall(pd)
Lc = cone15_length(50.0)
ref = design_bell_c2(0.8 * Lc, 10.0, 25.0, BellC2Spec())
check("(D): bell_mode direct_bspline は design_bell_c2 の壁 (同じ入力で同じ x_m・係数)", isinstance(wd, DirectBellWall)
      and des.ok and des.x_m == ref.x_m and np.array_equal(des.wall.coef, ref.wall.coef), f"x_m {des.x_m}")
check("(D): 壁の口 (x_in −L_U−L_pipe・x_e = L・スロート (0, 1)・r(x, 1))", wd.x_in == -4.0 and wd.x_e == pd.dv["L_over_rt"]["value"]
      and abs(float(wd.r(np.array([0.0]))[0]) - 1.0) < 1e-12 and abs(float(wd.r(np.array([0.0]), 1)[0])) < 1e-12
      and (wd.x_throat, wd.r_throat) == (0.0, 1.0)
      and np.array_equal(wd.drdx(np.array([1.0, 5.0])), des.wall.r(np.array([1.0, 5.0]), 1)))
check("(D): 端の丸め (1e-12) は端の値、端から大きく外れた x は例外",
      float(wd.r(np.array([wd.x_e + 1e-12]))[0]) == float(wd.r(np.array([wd.x_e]))[0])
      and raises(lambda: wd.r(np.array([wd.x_e + 1e-3])), ValueError) is not None)
cdd, _, _ = generate_axisym_mesh(wd, Mesh2DParams(**wn_base, wall_normal_layer=(0.02, 0.8), wall_normal_end_stations=8))
check("(D): mesh2d が (D) の壁で格子を作れる (壁の節点は壁の上)",
      np.allclose(cdd.reshape(ni_, nj_, 3)[:, -1, 1] / 0.14, wd.r(cdd.reshape(ni_, nj_, 3)[:, -1, 0] / 0.14), atol=1e-12))
pd_bad = load_problem(PROBLEM_D)
pd_bad.dv = {**pd.dv, "theta_max_deg": {**pd.dv["theta_max_deg"], "value": 20.0},
             "theta_e_deg": {**pd.dv["theta_e_deg"], "value": 5.0},
             "L_over_rt": {**pd.dv["L_over_rt"], "value": 0.6 * Lc}}
dbad = design_bell_c2(0.6 * Lc, 5.0, 20.0, BellC2Spec())   # V0-b の D_L60_tm20_te5 (成立区間が空)
e = None
try:
    runner.build_wall(pd_bad)
except runner.WallInfeasible as ex:
    e = ex
if dbad.ok:
    check("(D): 成立しない点の負例 (θ_max 20°・θ_e 5°・L60 が成立してしまった — 負例を選び直す)", False, f"x_m {dbad.x_m}")
else:
    check("(D): 成立区間が空の点は WallInfeasible (ValueError の派生) で理由と設計の記録を返す",
          e is not None and isinstance(e, ValueError) and "成立しない" in e.reason and e.record.get("design_status") == "INFEASIBLE",
          (e.reason[:120] if e else "例外なし"))
pd_bad = load_problem(PROBLEM_D)
pd_bad.geometry = {**pd.geometry, "Rd": 0.4}
check("(D): TOP のキー (Rd) を書いたら止める (黙って無視しない)", "効かない" in (raises(lambda: runner.build_wall(pd_bad), ValueError) or ""))
pd_bad = load_problem(PROBLEM_D)
pd_bad.dv = {**pd.dv, "theta_a_deg": {"value": 30.0, "min": 20.0, "max": 40.0}}
check("(D): dv の余分 (theta_a_deg) は止める", raises(lambda: runner.build_wall(pd_bad), ValueError) is not None)
pd_bad = load_problem(PROBLEM_D)
pd_bad.geometry = {**pd.geometry, "bell_mode": "tic"}          # moc_rao は 2026-10-12 に実装 (下の (M) の節)
check("(D): 未実装の bell_mode は止める", "未実装" in (raises(lambda: runner.build_wall(pd_bad), ValueError) or ""))
check("(D): TOP の問題 (bell_mode 無し) は従来の NozzleWall のまま", type(runner.build_wall(p)[0]).__name__ == "NozzleWall")
check("driver: 設計変数の並び (TOP は θ_a・θ_e・L、(D) は θ_max・θ_e・L)", DR.dv_order_of(p.raw) == DR.DV_ORDER
      and DR.dv_order_of(pd.raw) == ("theta_max_deg", "theta_e_deg", "L_over_rt"))
campd = DR.Campaign(PROBLEM_D, tmp / "camp_d", None)
check("driver: (D) の Campaign は dv の範囲を (D) の並びで読む", campd.dv_order == DR.DV_ORDER_DIRECT
      and campd.bounds[0].tolist() == [20.0, 45.0] and campd.feasible([25.0, 10.0, 0.8 * Lc]))
check("driver: (D) の成立しない点は事前フィルタで落ちる", (not dbad.ok) and not campd.feasible([20.0, 5.0, 0.6 * Lc]))
row = campd.evaluate([20.0, 5.0, 0.6 * Lc], "d_infeasible")          # forge は起動しない (壁の構築で止まる)
check("driver: (D) の成立しない点の評価は INFEASIBLE・理由 (note)・設計の記録 (wall_design) を台帳に返す",
      row["status"] == "FAIL" and row["fail_class"] == "INFEASIBLE" and "成立しない" in row["note"]
      and row.get("wall_design", {}).get("design_status") == "INFEASIBLE"
      and json.loads((tmp / "camp_d" / "ledger.jsonl").read_text().splitlines()[-1])["fail_class"] == "INFEASIBLE",
      row["note"][:120])

# ---------------------------------------------------------------- (M) 特性曲線モード (plan §4.4・§5.1 #6 (W1))
# geometry.bell_mode moc_rao: 生成器の壁ファイル (rao_m/rao_t6_parts/*_n321_mq2_wall.npz) を読み、(D) と同じ収縮部につなぐ。
# Euler の問題 YAML (problem_bell_tp_m_euler_l100/l70.yaml) で、壁・照合・接合・prepare の変換器の手前まで (変換器は偽物で止める)
from forge_design.geometry import rao_axisym as RA  # noqa: E402
from forge_design.geometry.wall import RaoBellWall  # noqa: E402
from forge_design.meshing.mesh2d import write_msh41_2d  # noqa: E402

RAO_DIR = CASE68 / "rao_m" / "rao_t6_parts"


class _StopAtConverter(Exception):
    pass


def _yaml_abs(src: Path, dst: Path, mut=None) -> Path:
    """問題 YAML を別の場所へ写す (gas.exhaust・geometry.rao_wall_file を絶対パスに)。`mut(d)` で中身を変える。"""
    d = yaml.safe_load(src.read_text())
    d["gas"]["exhaust"] = str((src.parent / d["gas"]["exhaust"]).resolve())
    d["geometry"]["rao_wall_file"] = str((src.parent / d["geometry"]["rao_wall_file"]).resolve())
    if mut:
        mut(d)
    dst.write_text(yaml.safe_dump(d, allow_unicode=True))
    return dst


for tag, lab, Lfrac in (("l100", "L100_n321_mq2", 1.0), ("l70", "L70_n321_mq2", 0.7)):
    PROBLEM_M = CASE68 / f"problem_bell_tp_m_euler_{tag}.yaml"
    pm = load_problem(PROBLEM_M)
    wm, recm = runner.build_wall(pm)
    jrec = json.loads((RAO_DIR / f"{lab}.json").read_text())
    with np.load(RAO_DIR / f"{lab}_wall.npz") as z:
        BE, KW = np.asarray(z["BE"]), np.asarray(z["kernel_wall"])
    check(f"(M) {lab}: bell_mode moc_rao は RaoBellWall (x_in −L_U−L_pipe・x_e = 壁ファイルの L・スロート (0, 1))",
          isinstance(wm, RaoBellWall) and wm.x_in == -4.0 and abs(wm.x_e - Lfrac * Lc) <= 1e-9 * Lc
          and (wm.x_throat, wm.r_throat) == (0.0, 1.0) and wm.validate() == [], f"x_e {wm.x_e!r}")
    # 同じ作り方 (run_rao_t4h.py の _rao_wall: BE の末尾を記録の E に、RaoWall の pspline) と係数・節点がビット一致
    P = BE[:, :2].copy()
    P[-1] = (jrec["x_E"], jrec["r_E"])
    ref_w = RA.RaoWall(1.0, np.deg2rad(jrec["theta_B_deg"]), P, np.deg2rad(jrec["theta_E_deg"]))
    check(f"(M) {lab}: 壁は T4 の評価 (run_rao_t4h) と同じ作り方で B-spline の係数・節点がビット一致",
          np.array_equal(ref_w.spl.c, wm.rao.spl.c) and np.array_equal(ref_w.spl.t, wm.rao.spl.t))
    # npz の壁点を通る: kernel の点は円弧 (閉じた式) の上、B と E は厳密、BE の内点は当てはめの残差 (記録の resid_max) 以内
    fi = recm["fit_info"]
    dk = float(np.max(np.abs(wm.r(KW[:, 0]) - KW[:, 1])))
    dBE = np.abs(wm.r(BE[1:-1, 0]) - BE[1:-1, 1])
    eB = abs(float(wm.r(np.array([wm.x_B]))[0]) - jrec["r_B"])
    eE = abs(float(wm.r(np.array([wm.x_e]))[0]) - jrec["r_E"])
    check(f"(M) {lab}: npz の kernel の壁点 (801 点) が円弧の上 (≤ 1e-9 r_t)", dk <= 1e-9, f"{dk:.2e}")
    check(f"(M) {lab}: 記録の B と E を通る (≤ 1e-12 r_t)", eB <= 1e-12 and eE <= 1e-12, f"B {eB:.1e}・E {eE:.1e}")
    check(f"(M) {lab}: npz の BE の内点を当てはめの残差以内で通る (最大 ≤ resid_max ≤ 1e-4 r_t)",
          float(dBE.max()) <= fi["resid_max"] + 1e-15 and fi["resid_max"] <= 1e-4,
          f"最大 {float(dBE.max()):.2e}・rms {float(np.sqrt(np.mean(dBE ** 2))):.2e} r_t")
    # 接合: スロートで C2、B で C1 (曲率は跳ぶ)。式の検算と、片側の値による独立の検算
    tj, bj = wm.throat_jump(), wm.B_jump()
    check(f"(M) {lab}: スロートで C2 (収縮部と円弧の r・r′・r″ の差 ≤ 1e-12)", max(abs(tj[k]) for k in (0, 1, 2)) <= 1e-12,
          f"{[f'{tj[k]:.1e}' for k in (0, 1, 2)]}")
    check(f"(M) {lab}: B で C1 (円弧と BE の r・r′ の差 ≤ 1e-12)", max(abs(bj[0]), abs(bj[1])) <= 1e-12,
          f"{bj[0]:.1e}・{bj[1]:.1e}")
    h = 1e-7
    dl, dr_ = (float(wm.r(np.array([v]), 1)[0]) for v in (wm.x_B - h, wm.x_B + h))
    sl, sr = (float(wm.r(np.array([v]), 2)[0]) for v in (wm.x_B - h, wm.x_B + h))
    t0l, t0r = (float(wm.r(np.array([v]), 2)[0]) for v in (-h, h))
    check(f"(M) {lab}: 片側の値で検算 (B の r′ の差 ≤ 1e-5・r″ の跳びが B_jump と 1e-5 以内、スロートの r″ の差 ≤ 1e-5)",
          abs(dr_ - dl) <= 1e-5 and abs((sr - sl) - bj[2]) <= 1e-5 and abs(t0r - t0l) <= 1e-5,
          f"Δr′ {dr_ - dl:.1e}・Δr″ {sr - sl:+.4f} (B_jump {bj[2]:+.4f})・スロート Δr″ {t0r - t0l:.1e}")
    check(f"(M) {lab}: 曲率の跳び Δκ_B < 0 を記録 (生成器の 6 点当てはめの dkappa_B と 5e-3 以内)",
          bj["dkappa"] < 0 and abs(bj["dkappa"] - jrec["dkappa_B"]) <= 5e-3, f"{bj['dkappa']:+.5f} / 記録 {jrec['dkappa_B']:+.5f}")
    ck = recm["checks"]
    check(f"(M) {lab}: 照合 (L・ε・R_t・排気ファイル) の記録", abs(ck["L_rel_diff"]) <= 1e-9 and abs(ck["eps_rel_diff"]) <= 1e-9
          and ck["R_t_wall"] == 1.0 and ck["gas"]["verified"], f"L {ck['L_rel_diff']:+.1e}・ε {ck['eps_rel_diff']:+.1e}")
    check(f"(M) {lab}: 来歴 (壁ファイル・記録の sha256、生成条件 n 321・m = q = 320・円弧・det_D > 0、キャンペーンの commit)",
          len(recm["wall_file_sha256"]) == 64 and len(recm["record_sha256"]) == 64 and recm["generator"]["n"] == 321
          and recm["generator"]["m"] == 320 and recm["generator"]["q"] == 320 and recm["generator"]["kernel"] == "arc"
          and recm["generator"]["accept_rule"] == "det_D > 0" and bool((recm["campaign"] or {}).get("commit")))

    # 照合の負例 (黙って合わせない)
    def _bad_m(_pm=PROBLEM_M, **kw):
        q = load_problem(_pm)
        for k, v in kw.items():
            setattr(q, k, v)
        return raises(lambda: runner.build_wall(q), (ValueError, FileNotFoundError)) or ""
    m = _bad_m(dv={"L_over_rt": {"value": Lfrac * Lc * (1 + 1e-6), "fixed": True}})
    check(f"(M) {lab}: L の不一致 (相対 1e-6) で止める", "L:" in m, m[:100])
    m = _bad_m(spec={**pm.spec, "eps": 49.0})
    check(f"(M) {lab}: ε の不一致で止める", "ε:" in m, m[:100])
    m = _bad_m(geometry={**pm.geometry, "R_t": 0.9})
    check(f"(M) {lab}: R_t の不一致で止める", "R_t" in m, m[:100])
    other = "L70_n321_mq2_wall.npz" if tag == "l100" else "L100_n321_mq2_wall.npz"
    m = _bad_m(geometry={**pm.geometry, "rao_wall_file": f"rao_m/rao_t6_parts/{other}"})
    check(f"(M) {lab}: 別の長さの壁ファイルは L で止める", "L:" in m, m[:100])
    m = _bad_m(geometry={k: v for k, v in pm.geometry.items() if k != "rao_wall_file"})
    check(f"(M) {lab}: rao_wall_file が無ければ止める", "rao_wall_file" in m, m[:100])
    m = _bad_m(geometry={**pm.geometry, "rao_wall_file": "rao_m/rao_t6_parts/missing_wall.npz"})
    check(f"(M) {lab}: 壁ファイルが無ければ止める", "無い" in m, m[:100])
    m = _bad_m(geometry={**pm.geometry, "Rd": 0.4})
    check(f"(M) {lab}: TOP のキー (Rd) を書いたら止める", "効かない" in m, m[:100])
    m = _bad_m(dv={**pm.dv, "theta_e_deg": {"value": 10.0, "min": 2.0, "max": 20.0}})
    check(f"(M) {lab}: dv の余分 (theta_e_deg) は止める", "dv" in m, m[:100])

# 排気ファイルが違う (中身の違うファイル) なら止める / 照合に落ちたら run dir を作らない / CPG では止める
dd = tmp / "m_gas"
dd.mkdir()
ex_other = dd / "exhaust_other.yaml"
ex_other.write_text((CASE68 / "exhaust_lox_lh2_pc100_of6.yaml").read_text() + "\n# 別のファイル (中身が違う)\n")
y = _yaml_abs(CASE68 / "problem_bell_tp_m_euler_l100.yaml", dd / "p.yaml", lambda d: d["gas"].update(exhaust=str(ex_other)))
m = raises(lambda: runner.build_wall(load_problem(y)), ValueError) or ""
check("(M): 壁ファイルを作った排気ファイルと問題の gas.exhaust が違えば止める", "排気ファイル" in m, m[:100])
y = _yaml_abs(CASE68 / "problem_bell_tp_m_euler_l100.yaml", dd / "p_L.yaml",
              lambda d: d["dv"]["L_over_rt"].update(value=0.9 * Lc))
m = raises(lambda: runner.prepare(y, dd / "run_bad"), ValueError) or ""
check("(M): 照合に落ちたら run dir を作らない (prepare)", "L:" in m and not (dd / "run_bad").exists(), m[:100])
y = _yaml_abs(CASE68 / "problem_bell_tp_m_euler_l100.yaml", dd / "p_cpg.yaml",
              lambda d: (d["gas"].pop("exhaust"), d["gas"].update(model="cpg", gamma=1.2, cp=3700.0),
                         d["spec"].update(Pt=1.0e7, Tt=3500.0), d["evaluate"].pop("model")))
m = raises(lambda: runner.prepare(y, dd / "run_cpg"), ValueError) or ""
check("(M): CPG (TP でない) の moc_rao は止める", "moc_rao" in m and not (dd / "run_cpg").exists(), m[:100])

# prepare の変換器の手前まで (変換器は偽物: 品質検査用の cell 変換は生成した格子の四角形の h5 を書き、本変換で止める)。
# 品質の判定は本物の check_mesh_quality.py (既定 AR ≤ 1000・skew ≤ 0.9、AR 緩和なし)
_orig_conv = runner.run_converter
for tag, lab in (("l100", "L100_n321_mq2"), ("l70", "L70_n321_mq2")):
    PROBLEM_M = CASE68 / f"problem_bell_tp_m_euler_{tag}.yaml"
    pm = load_problem(PROBLEM_M)
    wm, recm = runner.build_wall(pm)
    mpm = runner.mesh_params(pm, 0.14)
    cm, qm, bm = generate_axisym_mesh(wm, mpm)
    seen = []

    def _fake_conv(args, cwd, env=None, _c=cm, _q=qm, _seen=seen):
        out, run_dir = str(args[-1]), cwd
        _seen.append((out, (Path(run_dir) / "solverConfig.yaml").read_text()))
        if out == "nozzle_qc.h5":
            write_qc_h5(Path(run_dir) / out, _c, _q)
            return
        raise _StopAtConverter(out)
    runner.run_converter = _fake_conv
    rd = tmp / f"prep_m_{tag}"
    try:
        e = raises(lambda: runner.prepare(PROBLEM_M, rd), _StopAtConverter)
    finally:
        runner.run_converter = _orig_conv
    check(f"(M) prepare {tag}: 変換器の手前まで通る (品質検査用の cell 変換 → 本変換で停止)",
          e == "nozzle.h5" and [s[0] for s in seen] == ["nozzle_qc.h5", "nozzle.h5"], str(e))
    write_msh41_2d(tmp / f"ref_{tag}.msh", cm, qm, bm)
    check(f"(M) prepare {tag}: nozzle.msh は壁と mesh_params の格子 (バイト一致)",
          (rd / "nozzle.msh").read_bytes() == (tmp / f"ref_{tag}.msh").read_bytes())
    mq = (rd / "MESH_QUALITY.txt").read_text()
    vq = [ln for ln in mq.splitlines() if ln.startswith("VERDICT")]
    check(f"(M) prepare {tag}: Euler 格子の品質 PASS (check_mesh_quality.py 既定 AR ≤ 1000・skew ≤ 0.9)",
          bool(vq) and vq[-1].startswith("VERDICT: PASS") and "AR<=1000" in vq[-1].replace(" ", ""), vq[-1] if vq else mq[:200])
    rx = json.loads((rd / "MESH_AR_RELAX.json").read_text())
    check(f"(M) prepare {tag}: AR 緩和は使わない (required False・AR > 1000 のセル 0)",
          (not rx["required"]) and rx["n_high_ar"] == 0 and rx["inlet_planar"] and rx["outlet_planar"],
          f"AR 最大 {rx['max_ar']:.2f}・skew 最大 {rx['max_skew']:.3f}")
    Rg = cm[:, 1].reshape(mpm.ni, mpm.nj) / 0.14
    Xg = cm[:, 0].reshape(mpm.ni, mpm.nj)
    fr = (Rg[:, -1] - Rg[:, -2]) / Rg[:, -1]
    check(f"(M) prepare {tag}: 壁際を細分しない (全断面で第一セル / 局所半径 = 0.005、j 方向の線は x 一定)",
          float(np.max(np.abs(fr - 0.005))) <= 1e-12 and float(np.max(np.ptp(Xg, axis=1))) == 0.0,
          f"第一セル比 {fr.min():.6f}〜{fr.max():.6f}")
    sc = yaml.safe_load(seen[-1][1])
    check(f"(M) prepare {tag}: 本変換の solverConfig は Euler の node の soft 段 (viscMethod 0・turbulence none)",
          sc["mesh"]["discretization"] == "node" and sc["physProp"]["viscMethod"] == 0 and sc["turbulence"]["model"] == "none"
          and sc["space"]["convMethod"] == 0 and "transport" not in sc["physProp"])
    check(f"(M) prepare {tag}: 品質検査用の変換は cell", yaml.safe_load(seen[0][1])["mesh"]["discretization"] == "cell")
    check(f"(M) prepare {tag}: bcondConfig の物理壁は slip",
          yaml.safe_load((rd / "bcondConfig.yaml").read_text())["wall"]["kind"] == "slip")
    wdj = json.loads((rd / "WALL_DESIGN.json").read_text())
    ref = wdj["moc_reference"]
    i1 = pm.exhaust.ideal_1d(50.0, A_t=np.pi * 0.14 ** 2)
    jr = json.loads((RAO_DIR / f"{lab}.json").read_text())
    F_cs_N = jr["F_cs"] * i1["rho_star"] * i1["a_star"] ** 2 * 0.14 ** 2
    check(f"(M) prepare {tag}: WALL_DESIGN.json の V2 参照値 F_cs [N] = F_cs·ρ*a*²r_t² (相対 1e-12)、η = F/F_1D",
          abs(ref["F_cs_N"] / F_cs_N - 1) <= 1e-12 and abs(ref["eta_F_cs"] - F_cs_N / i1["F"]) <= 1e-12
          and wdj["bell_mode"] == "moc_rao" and wdj["wall_file_sha256"] == recm["wall_file_sha256"],
          f"F_cs {ref['F_cs_N']:.1f} N・η {ref['eta_F_cs']:.6f}・F_wall {ref['F_wall_N']:.1f} N・η {ref['eta_F_wall']:.6f}")
    check(f"(M) prepare {tag}: 換算の検算 (MOC の気体の表だけの η と相対 1e-6、ρ*a*² の 2 通りが相対 1e-6)",
          abs(ref["check"]["eta_rel_diff"]) <= 1e-6 and abs(ref["check"]["scale_rel_diff"]) <= 1e-6,
          f"η {ref['check']['eta_rel_diff']:+.1e}・ρ*a*² {ref['check']['scale_rel_diff']:+.1e}")
    check(f"(M) prepare {tag}: 変換器の前に止めたので nozzle.h5・prepare_info.json は無い",
          not (rd / "nozzle.h5").exists() and not (rd / "prepare_info.json").exists())

# ---------------------------------------------------------------- Euler の run の推力の時系列・metrics・ゲート (plan §5.1 #6 (W1))
# 合成の Euler run (直管の格子・一様流、ṁ/ṁ_1D ≈ 0.99・η ≈ 0.98 になる状態) を problem_bell_tp_m_euler_l100.yaml で集計し、
# runner.collect → driver.tp_gates (flow_model euler) を通す。forge・変換器は使わない
PROBLEM_EU = CASE68 / "problem_bell_tp_m_euler_l100.yaml"
pe = load_problem(PROBLEM_EU)
id_e = runner.ideal_1d_tp(pe)
A_p = math.pi * 0.5 ** 2
u_e = 4000.0
ro_e = 0.99 * id_e["mdot"] / (u_e * A_p)                          # ṁ/ṁ_1D = 0.99
P_e = 0.98 * id_e["F"] / A_p - ro_e * u_e * u_e                   # η = 0.98
eu = tmp / "eu_run"
eu.mkdir()
cols_e = ["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe"]   # Euler: k・ω の列なし
sm_e = StageManifest(eu)
for tag_e, st_e, hist_e in (("S1_soft", "soft", pre), ("S2_mid", "mid", pre), ("S3_main", "main", falling)):
    with open(eu / f"residual_history_{tag_e}.csv", "w", newline="") as fh:
        w_ = csv.writer(fh); w_.writerow(cols_e)
        for s_, v_ in enumerate(hist_e):
            w_.writerow([s_, 4, "outer_end", v_, 1e3 * v_, 5e2 * v_, 0.0, 1e7 * v_])
    sm_e.add(tag_e, runner._solver_config_tp(pe, len(hist_e), len(hist_e) // 20 or 1, st_e), runner._bcond_config_tp(pe),
             history=f"residual_history_{tag_e}.csv")
sm_e.write()
(eu / "solverConfig.yaml").write_text(runner._solver_config_tp(pe, 20000, 1000, "main"))
write_mesh_h5(eu / "nozzle.h5", pc, pq, pb)
for st_ in range(0, 20001, 1000):
    write_euler_res(eu / f"res_{st_}.h5", NP, ro_e, u_e, 0.0, P_e, grad=1.0e3)
    if st_ > 0:
        write_slip_dump(eu / f"res_wall_3_{st_}.h5", pc, pwall)
(eu / "MESH_QUALITY.txt").write_text("(primal/cell 変換で検査)\nVERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
(eu / "WALL_DESIGN.json").write_text((tmp / "prep_m_l100" / "WALL_DESIGN.json").read_text())    # 上の prepare が書いたもの
rows_e = runner.thrust_history_tp(pe, eu)
check("Euler run: 推力の時系列は step 0 を除く全スナップショット (NS と同じ step の集合)・全行 slip",
      [r_["step"] for r_ in rows_e] == list(range(1000, 20001, 1000)) and all(r_["wall_model"] == "slip" for r_ in rows_e))
check("Euler run: 合成の一様流で推力の 2 会計が閉じる (|cv| < 1e-12)・η 0.98・ṁ/ṁ_1D 0.99",
      max(abs(r_["cv_residual"]) for r_ in rows_e) < 1e-12 and math.isclose(rows_e[-1]["eta_exit"], 0.98, rel_tol=1e-12)
      and math.isclose(rows_e[-1]["mdot_ratio_1d"], 0.99, rel_tol=1e-12))
me = runner.collect(PROBLEM_EU, eu)
mr = me.get("moc_reference", {})
check("Euler run: metrics に flow_model euler・wall_model slip・cv_closable True・thrust_history (20 点)",
      me.get("flow_model") == "euler" and me.get("wall_model") == "slip" and me.get("cv_closable") is True
      and me.get("thrust_history", {}).get("n") == 20 and (eu / "thrust_history.csv").exists(), str({k: me.get(k) for k in ("flow_model", "wall_model", "cv_closable")}))
check("Euler run: metrics に MOC の参照値との差 η − η_cs (出口面・壁) と推力の相対差・ṁ/ṁ_1D の差",
      math.isclose(mr.get("eta_exit_minus_eta_cs", 1.0), me["eta_exit"] - mr["eta_F_cs"], rel_tol=0, abs_tol=0)
      and math.isclose(mr["eta_wall_minus_eta_cs"], me["eta_wall"] - mr["eta_F_cs"], rel_tol=0, abs_tol=0)
      and math.isclose(mr["thrust_exit_rel_diff_cs"], me["thrust_exit_N"] / mr["F_cs_N"] - 1.0, rel_tol=0, abs_tol=0)
      and math.isclose(mr["mdot_ratio_1d_minus_moc"], me["mdot_ratio_1d"] - mr["mdot_ratio_1d"], rel_tol=0, abs_tol=0),
      f"η − η_cs {mr.get('eta_exit_minus_eta_cs', float('nan')):+.6f}")
ge = DR.tp_gates(eu, me, flow_model="euler")
gj = json.loads((eu / "GATES.json").read_text())
wr_e = gj["gates"]["wall_resolution"]
check("Euler の GATES: 壁解像 (y1+) は回さず対象外と理由を書いて ok (WALL_RESOLUTION_VERDICT.txt を作らない)",
      wr_e["ok"] and wr_e["skipped"] and "対象外" in wr_e["verdict"] and "Euler" in wr_e["verdict"]
      and not (eu / "WALL_RESOLUTION_VERDICT.txt").exists() and gj.get("flow_model") == "euler", wr_e["verdict"][:60])
check("Euler の GATES: (i) 残差 (k・ω の列なし)・(ii) 準定常・(iii) 物理・(iv) CV・(v) メッシュ品質は NS と同じ判定で通り PASS",
      ge["status"] == "PASS" and all(ge["gates"][k]["ok"] for k in ("convergence", "quasisteady", "physics", "cv", "mesh_quality")),
      f"{ge['status']} {ge['note']} {ge['gates']['convergence'].get('verdict')}")
m = raises(lambda: DR.tp_gates(eu, me, flow_model="ns"), ValueError) or ""
check("Euler の run を NS として判定しようとすると止める (wall_model と flow_model の食い違い)", "食い違う" in m, m[:80])
m = raises(lambda: DR.tp_gates(eu, {**me, "wall_model": "viscous"}, flow_model="euler"), ValueError) or ""
check("NS の会計 (viscous) を Euler として判定しようとすると止める", "食い違う" in m, m[:80])
check("tp_gates: 未知の flow_model は止める", raises(lambda: DR.tp_gates(eu, me, flow_model="rans"), ValueError) is not None)
me_nothrust = {k: v for k, v in me.items() if k not in ("eta_exit", "wall_model", "eta_cf", "cv_residual", "cv_closable")}
ge2 = DR.tp_gates(eu, me_nothrust, flow_model="euler")
check("tp_gates: 推力の会計が無い Euler の metrics は照合せず、(iii)(iv) が欠損で落とす (例外にしない)",
      ge2["status"] == "FAIL" and not ge2["gates"]["physics"]["ok"] and not ge2["gates"]["cv"]["ok"]
      and ge2["gates"]["wall_resolution"]["skipped"], ge2["note"])

# 評価ごとの YAML の相対パス: geometry.rao_wall_file も直す (PROBLEM_PATH_KEYS)
check("(b) geometry.rao_wall_file も問題 YAML の所在で絶対パスに",
      DR.absolutize_problem_paths({"geometry": {"rao_wall_file": "rao_m/x_wall.npz"}}, "/a/b")["geometry"]["rao_wall_file"]
      == "/a/b/rao_m/x_wall.npz")
raw_m = yaml.safe_load(PROBLEM_EU.read_text())
ab_m = DR.absolutize_problem_paths(raw_m, PROBLEM_EU.parent)
(tmp / "elsewhere").mkdir()
(tmp / "elsewhere" / "m.yaml").write_text(yaml.safe_dump(ab_m, allow_unicode=True))
wm2, rec2 = runner.build_wall(load_problem(tmp / "elsewhere" / "m.yaml"))
check("(b) 書き出した (M) の問題 YAML を別ディレクトリから読める (壁ファイル・排気ファイルが解決し、同じ壁)",
      Path(ab_m["geometry"]["rao_wall_file"]).is_absolute() and Path(ab_m["geometry"]["rao_wall_file"]).exists()
      and rec2["wall_file_sha256"] == json.loads((tmp / "prep_m_l100" / "WALL_DESIGN.json").read_text())["wall_file_sha256"]
      and np.array_equal(wm2.rao.spl.c, runner.build_wall(load_problem(PROBLEM_EU))[0].rao.spl.c))

# ---------------------------------------------------------------- driver: TP は全点冷間・停止規則 (2026-10-11 diagnostician)
_seen = {}
_orig_staged = runner.run_staged_tp
def _fake_staged(prob, run_dir, warm_from=None, **kw):      # forge は起動しない: 渡された warm_from を記録して段の失敗にする
    _seen["warm_from"] = warm_from
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    (Path(run_dir) / "res_nan_41558.h5").write_bytes(b"")
    raise RuntimeError("TP の段 S3_main が失敗 (rc=1)")
runner.run_staged_tp = _fake_staged
try:
    campw = DR.Campaign(PROBLEM_D, tmp / "camp_cold", tmp / "nonexistent_seed.h5")
    campw.rows.append({"tag": "doe_000", "x": [25.0, 10.0, 0.8 * Lc], "status": "PASS", "run_dir": str(tmp / "x"),
                       "res_file": "res_60000.h5", "fail_class": None})
    r1 = campw._evaluate_tp([26.0, 10.0, 0.8 * Lc], "doe_001")
    check("driver: TP は PASS 点や --seed-res があっても warm start を渡さない (全点冷間)",
          "warm_from" in _seen and _seen["warm_from"] is None and r1["warm_from"] == "", str(_seen))
    check("driver: 段の失敗は RETRYABLE、note に最初の NaN の step", r1["fail_class"] == "RETRYABLE" and "NaN step 41558" in r1["note"],
          r1["note"])
    campw.rows.append({"tag": "doe_002", "fail_class": "RETRYABLE"})
    check("driver: RETRYABLE 2 連続では止めない", raises(lambda: campw._stop_on_repeated_failures(), SystemExit) is None)
    campw.rows.append({"tag": "doe_003", "fail_class": "RETRYABLE"})
    check("driver: RETRYABLE 3 連続で止める", raises(lambda: campw._stop_on_repeated_failures(), SystemExit) is not None)
    campw.rows[-1] = {"tag": "doe_003", "fail_class": "SUSPECT"}
    check("driver: 種類の違う失敗の 3 連続は止めない", raises(lambda: campw._stop_on_repeated_failures(), SystemExit) is None)
    rows0 = list(campw.rows)
    campw.rows = [{"tag": f"doe_{i:03d}", "fail_class": ("RETRYABLE" if i < 8 else None)} for i in range(20)]
    check("driver: DOE 40 点で失敗 8 点 (= 20 %) は止めない", raises(lambda: campw._stop_on_doe_failure_rate(40), SystemExit) is None)
    campw.rows.append({"tag": "doe_020", "fail_class": "SUSPECT"})
    check("driver: DOE 40 点で失敗 9 点 (> 20 %) で止める", raises(lambda: campw._stop_on_doe_failure_rate(40), SystemExit) is not None)
    campw.rows = [{"tag": f"inf_00_{i}", "fail_class": "RETRYABLE"} for i in range(9)]
    check("driver: infill の失敗は DOE の割合に数えない", raises(lambda: campw._stop_on_doe_failure_rate(40), SystemExit) is None)
    campw.rows = rows0
finally:
    runner.run_staged_tp = _orig_staged

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
