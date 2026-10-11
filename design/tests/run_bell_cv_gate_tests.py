#!/usr/bin/env python3
"""③ベル 2 モードの共通評価器 — 推力の CV 会計とゲート (plan tooling-nozzle-bell-two-mode §4.1・§5.1 #2b)。

forge を起動しない単体試験 (全流れ解での確認は case/68 の run):
- `metrics.extract.thrust_cv`: 合成した node メッシュ (TOP ベルの壁・直管) と場で、
  静止気体の閉曲面で Σ(p − p_a) n_x dA = 0 (CV 残差 0、p_a に依らない)、面積 πr²、外向き法線の向き、一様流の直管で
  壁会計 = 出口面会計・ṁ_in = ṁ_out、twall の規約 (node: 流体に働く traction = −twall、cell: +twall) と摩擦の推力寄与の符号、
  P の出所 (VALUE/P 優先、無ければ TP は e_sens(T) の逆算・CPG は γ)、粘性応力・壁ダンプが無いときは閉じない扱い。
- `opt.driver` のゲート (i)〜(v): 合成した残差履歴 + StageManifest で (i) が STALLED・未収束・manifest 無しを落とし
  PASS を通す、(ii) 推力の時系列の DRIFTING・NONFINITE・列欠損を落とす、(iii) 物理範囲、(iv) CV 残差超過・閉じない会計を落とす、
  (v) メッシュ品質の VERDICT の読み、総合の分類 (RETRYABLE / SUSPECT / INFEASIBLE)。
- 引き継ぎ (b): 評価ごとに書き出す問題 YAML の相対パス (gas.exhaust・gas.species_db) を絶対パスに直す。
- runner: TP の段の計画 (`tp_stage_plan`) の検査、TP の config の出力 (勾配の extraFields・pRef を書かない)、
  格子のスロート細分キー (`mesh_params`) は書いたときだけ渡る。

実行: design/.venv-opt の Python で `python design/tests/run_bell_cv_gate_tests.py`。
"""
import csv
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
check("(i) STALLED (1.5 桁で横ばい) を落とす", (not g["ok"]) and g["class"] == "STALLED", g["verdict"])
g = DR.gate_convergence(make_run(tmp / "conv_nomani", falling, pre, manifest=False))
check("(i) stage_manifest.json が無ければ判定不能 = 不合格", not g["ok"], g["verdict"])
nanh = falling[:-5] + [float("nan")] * 5
g = DR.gate_convergence(make_run(tmp / "conv_nan", nanh, pre))
check("(i) NaN は DIVERGED で落とす", (not g["ok"]) and g["class"] == "DIVERGED", g["verdict"])


# ---------------------------------------------------------------- ゲート (ii): 推力の時系列
def write_series(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(("step",) + runner.THRUST_HISTORY_COLS)
        for r in rows:
            w.writerow([r["step"]] + [repr(float(r.get(c, 1.0))) for c in runner.THRUST_HISTORY_COLS])


steps = np.arange(1, 13) * 1000
base = dict(thrust_exit_N=1.135e6, thrust_wall_N=1.135e6, eta=0.9792, eta_wall=0.9792, mdot_exit_kg_s=268.8,
            mdot_in_kg_s=268.8, mdot_ratio_1d=0.9925, cv_residual=4.0e-4, friction_N=-1.18e4)
steady = [{"step": int(st), **{k: v * (1 + 1e-6 * math.sin(i)) for k, v in base.items()}} for i, st in enumerate(steps)]
(tmp / "qs_ok").mkdir()
write_series(tmp / "qs_ok" / "thrust_history.csv", steady)
gq = DR.gate_quasisteady(tmp / "qs_ok" / "thrust_history.csv")
check("(ii) 平坦な系列は ALL STEADY で通る", gq["ok"] and set(gq["per_column"].values()) == {"STEADY"}, gq["overall"])
drift = [{"step": int(st), **base, "thrust_exit_N": base["thrust_exit_N"] * (1 + 0.03 * i), "eta": base["eta"] * (1 + 0.03 * i)}
         for i, st in enumerate(steps)]
(tmp / "qs_drift").mkdir()
write_series(tmp / "qs_drift" / "thrust_history.csv", drift)
gq = DR.gate_quasisteady(tmp / "qs_drift" / "thrust_history.csv")
check("(ii) DRIFTING (推力が単調に 33 % 動く; 末尾窓で 9 %) を落とす", (not gq["ok"]) and gq["per_column"].get("thrust_exit_N") == "DRIFTING",
      str(gq["per_column"]))
cvd = [{"step": int(st), **base, "cv_residual": 4.0e-4 * (1 + 0.2 * i)} for i, st in enumerate(steps)]
(tmp / "qs_cv").mkdir()
write_series(tmp / "qs_cv" / "thrust_history.csv", cvd)
gq = DR.gate_quasisteady(tmp / "qs_cv" / "thrust_history.csv")
check("(ii) CV 残差の時系列の DRIFTING も落とす", (not gq["ok"]) and gq["per_column"].get("cv_residual") == "DRIFTING",
      str(gq["per_column"]))
nanr = [dict(r) for r in steady]
nanr[-1]["thrust_wall_N"] = float("nan")
(tmp / "qs_nan").mkdir()
write_series(tmp / "qs_nan" / "thrust_history.csv", nanr)
gq = DR.gate_quasisteady(tmp / "qs_nan" / "thrust_history.csv")
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
m_ok = {"eta_cf": 0.979, "mdot_ratio_1d": 0.9925, "cv_residual": 4e-4, "cv_closable": True}
check("(iii) 物理範囲内は通る", DR.gate_physics(m_ok)["ok"])
check("(iii) ṁ/ṁ_1D 0.90 (不始動の偽アトラクタ) を落とす", not DR.gate_physics({**m_ok, "mdot_ratio_1d": 0.90})["ok"])
check("(iii) η 1.05 を落とす", not DR.gate_physics({**m_ok, "eta_cf": 1.05})["ok"])
check("(iii) 欠損・非有限は落とす", not DR.gate_physics({**m_ok, "eta_cf": None})["ok"]
      and not DR.gate_physics({**m_ok, "mdot_ratio_1d": float("nan")})["ok"])
check("(iv) |CV 残差| 4e-4 は通る", DR.gate_cv(m_ok)["ok"])
check("(iv) |CV 残差| 2e-3 (> 0.1 %) を落とす", not DR.gate_cv({**m_ok, "cv_residual": -2e-3})["ok"])
check("(iv) ちょうど 1e-3 は通る (≤)", DR.gate_cv({**m_ok, "cv_residual": 1e-3})["ok"])
check("(iv) 閉じない会計 (cv_closable False) は落とす", not DR.gate_cv({**m_ok, "cv_closable": False})["ok"])
check("(iv) NaN は落とす", not DR.gate_cv({**m_ok, "cv_residual": float("nan")})["ok"])
(tmp / "mq").mkdir()
(tmp / "mq" / "MESH_QUALITY.txt").write_text("(primal/cell 変換で検査)\nmax AR 836\nVERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
check("(v) MESH_QUALITY の PASS を読む", DR.gate_mesh_quality(tmp / "mq")["ok"])
(tmp / "mq" / "MESH_QUALITY.txt").write_text("VERDICT: FAIL (AR 1200 > 1000)\n")
check("(v) MESH_QUALITY の FAIL を落とす", not DR.gate_mesh_quality(tmp / "mq")["ok"])
check("(v) MESH_QUALITY.txt が無ければ落とす", not DR.gate_mesh_quality(tmp / "nonexistent")["ok"])
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

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
