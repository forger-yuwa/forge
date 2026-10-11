#!/usr/bin/env python3
"""③ベル 2 モードの共通評価器 — 推力の CV 会計とゲート (plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2b)。

forge を起動しない単体試験 (全流れ解での確認は case/68 の run):
- `metrics.extract.thrust_cv`: 合成した node メッシュ (TOP ベルの壁・直管) と場で、
  静止気体の閉曲面で Σ(p − p_a) n_x dA = 0 (CV 残差 0、p_a に依らない)、面積 πr²、外向き法線の向き、一様流の直管で
  壁会計 = 出口面会計・ṁ_in = ṁ_out、twall の規約 (node: 流体に働く traction = −twall、cell: +twall) と摩擦の推力寄与の符号、
  P の出所 (VALUE/P 優先、無ければ TP は e_sens(T) の逆算・CPG は γ)、粘性応力・壁ダンプが無いときは閉じない扱い。
- `opt.driver` のゲート (i)〜(v) (2026-10-11 組み替え、plan §4.1・§5.1 #2c'): 合成した残差履歴 + StageManifest で (i) は PASS・
  STALLED・still converging を通し (STALLED は併記)、DIVERGED・有意な末尾上昇 (rising)・manifest 無し (判定不能) を落とす、
  (ii) 推力系 4 列 (1.5e-4) の DRIFTING・振動・単調で漸近値が離れた列・NONFINITE・列欠損・末尾窓 8 点未満を落とし、ṁ は 1e-3 で判定、
  cv_residual は列にしない、(iii) 物理範囲、(iv) 最終 |CV 残差| 1e-4 超過・閉じない会計を落とす、(v) メッシュ品質の VERDICT の読みと
  AR > 1000 の格子の例外の条件 (MESH_AR_RELAX.json)、総合の分類 (RETRYABLE / SUSPECT / INFEASIBLE)。
- 壁法線の近壁層 (#2d): mesh2d の `wall_normal_end_stations` (入口・出口面を平面に)、runner の格子キーの受け渡し、
  `mesh_ar_relax_check` (AR > 1000 のスキュー・面の平面性) の負例。
- (D) 直接モード (#3f): `geometry.bell_mode: direct_bspline` の壁 (`DirectBellWall`)、成立しない点の `WallInfeasible` (理由つき)、
  TOP のキー・dv の誤りを止める、driver の設計変数の並び。
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
pd_bad.geometry = {**pd.geometry, "bell_mode": "moc_rao"}
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
