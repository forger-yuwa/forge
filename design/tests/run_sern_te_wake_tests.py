#!/usr/bin/env python3
"""格子 B (後縁下流の中間線の局所変形、`mesh3d.te_wake_blend_H`) の試験
(plan convection-zero-thickness-edge-reconstruction §4.2・§5.1 #2c、codex diagnose 2026-10-08 te-grid-kink)。forge は回さない。

合格 (測る前に固定):
  (a) 既定 (無効) は**変更前のメッシャ** (commit ae856141 の mesh_sern3d.py) と座標・ヘキサ・境界面がビット一致。0 を明示しても同じ
  (b) 有効: 節点数・ヘキサ・境界面の接続・x と z の配列がビット一致、y が動くのは変形区間 (L_cowl, L_cowl + L_b) の station の
      主ブロックだけで下端 (j=0)・上線 (j=NJ−1) を含まない、外部境界と壁の節点は動かない (z 一定の面内 sym・side_far だけ面内で動く)、
      閉性・ヤコビアンの符号・float32 での新たな節点衝突なし
  (c) 中間線の端点: L_cowl で位置が後縁・勾配が**下側壁面の最後の辺の傾き** (座標から測る。上側壁面とは別)、
      L_cowl + L_b で旧中間線の位置と勾配の両方に一致。実際の最初の下流の辺の傾き = 弦の式 m0 + (m1 − m0)(2t₁ − t₁²)
  (d) 不正な指定は生成を失敗させる (負・出口に届く・区間に station が無い)
  (e) runner_sern3d.prepare が問題 YAML の `mesh3d.te_wake_blend_H` をメッシャに渡す (無ければ 0)
  (f) 格子の確認の道具 (te_wake_grid_check.py): msh と info から (v) 同一性・動いた範囲、(i) 後縁の最初の下流の辺、(ii) 折れ角を出す
  (g) 初期場の道具 (restart_field_deformed.py): 不動節点は SRC の保存量とビット一致、動いた節点は interp_field.py の結果と
      ビット一致、wall_dist・/AUX は DST のまま、双子への最寄り・接続違い・節点数違いは拒否、全節点不動なら restart_field.py と同じ
  (h) 曲線版の定数 (来歴): 有効時だけ info に te_wake_curve_version (= mesh_sern3d.TE_WAKE_CURVE_VERSION) が入り、既定 0 の info には
      入らない (既定の格子・info のビット一致は (a))
  (i) --admission (plan §6.0 の投入条件の表) をメモリ上の g1 (生産 YAML problem_3d_prod_m6on_g1.yaml の設計 MOC、L_b 1.0 H、座標は変換後と
      同じ float32 に丸める) に掛ける: 幾何の 6 行 (最初の下流の辺・後縁の折れ・復帰区間の折れ・固定すべき量・層の Δy・primal の品質) が PASS、
      codex の g3 の再生成値 (後縁の折れ 2.20°・復帰 5.08°・スケール済みヤコビアン 0.67245・skew 0.44706 → 0.53048) と一致、
      入力が最終 node 格子でない (前提)・品質と双対のファイルが無い → 判定不能、総合は UNDECIDABLE (PASS にならない)
  (j) 変換した粗い格子 (g1 の x・y 分布で nz・ni_noz を減らした格子を runner_sern3d.prepare で作り、最終 node の h5 にする) の
      run ディレクトリに CLI の --admission を掛ける: check_dual_closure の出力を渡せば ADMISSION VERDICT: PASS (終了コード 0)、
      渡さなければ UNDECIDABLE (2)、L_b 0.3 H の B は FAIL (1)、info・MESH_QUALITY が無い・msh 入力は UNDECIDABLE
  (k) 意図的に壊した格子で該当の行が FAIL・入力の欠落で判定不能: 復帰区間の折れ・x・上流・壁・外部境界・接続・第 1 層・反転ヘキサ・
      節点数・境界面の接続なし・msh の座標、品質の VERDICT の順位・閾値違い・セル数違い、変形領域の新しい AR 超過、双対の FAIL・tol 違い・CV 数違い
  (l) --provenance-out: 格子署名 = mark_zero_thickness_edges.mesh_signature の再計算値、曲線版・L_b・設定の sha256・commit の記録の有無、
      判定の行。--same-mesh: /VALUE を書き換えた同じ格子は YES (0)、A と B は NO (1)、読めない入力は UNDECIDABLE (2)

  design/.venv-opt/bin/python design/tests/run_sern_te_wake_tests.py
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
CHECK = ROOT / "case" / "46.sern_design" / "diag" / "te_wake_grid_check.py"
RFD = TOOLS / "restart_field_deformed.py"
sys.path.insert(0, str(HERE.parent))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
from forge_design.meshing.mesh_sern3d import (PHYS_SERN3D, TE_WAKE_CURVE_VERSION, SernMesh3DParams,  # noqa: E402
                                              generate_sern_mesh3d, write_msh41_3d)

REF_COMMIT = "ae8561411db51598b9c3b7fe0f6a15c9a2518ceb"   # 本変更の直前に mesh_sern3d.py を変えた commit (L_sw_exact)
FAIL = 0
tmp = Path(tempfile.mkdtemp(prefix="te_wake_tests_"))


def check(name, cond, detail="", only_on_fail=False):
    global FAIL
    show = detail and (not only_on_fail or not cond)
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if show else ""))
    if not cond:
        FAIL += 1


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    except Exception as e:  # noqa: BLE001
        print(f"    (別の例外 {type(e).__name__}: {e})")
    return False


def load_ref_mesher():
    """変更前の mesh_sern3d.py を git から読んで別名のモジュールにする (相対 import は forge_design.meshing で解決)。"""
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{REF_COMMIT}:design/forge_design/meshing/mesh_sern3d.py"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    f = tmp / "_mesh_sern3d_ref.py"; f.write_text(r.stdout)
    spec = importlib.util.spec_from_file_location("forge_design.meshing._mesh_sern3d_ref", f)
    m = importlib.util.module_from_spec(spec); m.__package__ = "forge_design.meshing"
    sys.modules[spec.name] = m          # dataclass (from __future__ import annotations) が自分のモジュールを引く
    spec.loader.exec_module(m)
    return m


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def same_B(B1, B2):
    return set(B1) == set(B2) and all(np.array_equal(np.asarray(B1[k], dtype=np.int64).reshape(-1, 4) if len(B1[k]) else np.zeros((0, 4)),
                                                     np.asarray(B2[k], dtype=np.int64).reshape(-1, 4) if len(B2[k]) else np.zeros((0, 4)))
                                      for k in B1)


def closure(hexes, B):
    from collections import Counter
    FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    cnt = Counter(frozenset(h[i] for i in f) for h in hexes for f in FACES)
    bnd = Counter(frozenset(q) for faces in B.values() for q in faces)
    unshared = [f for f, c in cnt.items() if c == 1]
    return sum(1 for f in unshared if bnd.get(f, 0) != 1), sum(1 for f in bnd if cnt.get(f, 0) != 1)


# ------------------------------------------------------------------------------------------------ 共通の格子
k = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                             x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
d = k.design_ramp(M_c=3.3, f=0.45)
THB = float(k.TH[-1, 0])
# 生産に近い構成: 有限厚 (下側と上側の壁面の傾きが違う)・カウル板の自由な側端・外部流ブロック・機体側面
P0 = SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=60, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5, L_sw=0.67,
                      interface_angle=THB, top_ext_angle=d.info["theta_e"], cowl_thickness=0.005, ext_top=True, nj_ext_top=9,
                      vehicle_taper=0.3, scale=0.1)
LB = 0.3
cA, hA, BA, iA, ymA = generate_sern_mesh3d(d, P0)
cB, hB, BB, iB, ymB = generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=LB))
NJ, nz, jm, it, ni = iA["NJ"], iA["nz"], iA["jm"], iA["i_te"], iA["ni"]
NB = ni * NJ * nz
base = lambda i, j, kk: (i * NJ + j) * nz + kk
H = P0.scale
print(f"--- 格子: nodes {iA['nodes']} hex {len(hA)} ni {ni} NJ {NJ} nz {nz} i_te {it}; B: L_b {LB}, station {iB['te_wake_n_stations']} "
      f"[i {iB['te_wake_i_first']}..{iB['te_wake_i_last']}], m0 {iB['te_wake_slope_te_deg']:.4f}°, m1 {iB['te_wake_slope_end_deg']:.4f}°")

# ------------------------------------------------------------------------------------------------ (a)
ref = load_ref_mesher()
if ref is None:
    check("(a) 変更前のメッシャを git から読めた", False, REF_COMMIT)
else:
    for label, prm in (("生産に近い構成", P0),
                       ("厚さ 0・外部流ブロックなし", replace(P0, cowl_thickness=0.0, ext_top=False, L_sw=None)),
                       ("外側空間なし (nz_out=0)", replace(P0, nz_out=0, ext_top=False))):
        prm_ref = ref.SernMesh3DParams(**{f: getattr(prm, f) for f in ref.SernMesh3DParams.__dataclass_fields__})
        c0, h0, B0, i0, y0 = ref.generate_sern_mesh3d(d, prm_ref)
        c1, h1, B1, i1, y1 = generate_sern_mesh3d(d, prm)
        check(f"(a) {label}: 既定は変更前と座標・ヘキサ・境界面がビット一致",
              np.array_equal(c0, c1) and np.array_equal(h0, h1) and same_B(B0, B1))
        xq = np.linspace(-0.5, float(i0["x_out"]), 997)
        check(f"(a) {label}: 返す中間線の関数も同じ値", np.array_equal(y0(xq), y1(xq)))
        check(f"(a) {label}: info は te_wake_blend_H 0 が増えただけ",
              {kk: v for kk, v in i1.items() if kk != "te_wake_blend_H"} == i0 and i1["te_wake_blend_H"] == 0.0)
c0b, h0b, B0b, _, _ = generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=0.0))
check("(a) 0 を明示しても既定とビット一致", np.array_equal(c0b, cA) and np.array_equal(h0b, hA) and same_B(B0b, BA))

# ------------------------------------------------------------------------------------------------ (b)
check("(b) 節点数・ヘキサの接続・境界面の接続が同一", cA.shape == cB.shape and np.array_equal(hA, hB) and same_B(BA, BB))
check("(b) x・z 配列がビット一致", np.array_equal(cA[:, 0], cB[:, 0]) and np.array_equal(cA[:, 2], cB[:, 2]))
mv = np.where(cA[:, 1] != cB[:, 1])[0]
ii, jj = mv // (NJ * nz), (mv // nz) % NJ
check("(b) y が動くのは主ブロックだけ", mv.size > 0 and np.all(mv < NB), f"{mv.size} 節点、主ブロック外 {int(np.sum(mv >= NB))}")
check("(b) 動く station = 変形区間 (L_cowl, L_cowl + L_b) の station",
      set(np.unique(ii).tolist()) == set(range(iB["te_wake_i_first"], iB["te_wake_i_last"] + 1)), f"i {sorted(set(ii.tolist()))}")
check("(b) 下端 (j=0) と上線 (j=NJ−1) は動かない", jj.min() >= 1 and jj.max() <= NJ - 2, f"j {jj.min()}..{jj.max()}")
xs = cA[:NB:NJ * nz, 0] / H
check("(b) 動く station の x は (L_cowl, L_cowl + L_b) の中", np.all((xs[np.unique(ii)] > iA["L_cowl"]) & (xs[np.unique(ii)] < iA["L_cowl"] + LB)))
check("(b) y_bot は旧中間線の x_out の値のまま", iB["y_bot"] == iA["y_bot"], f"{iB['y_bot']} vs {iA['y_bot']}")
mset = np.zeros(len(cA), bool); mset[mv] = True
moved_tags = {nm: int(mset[np.unique(np.asarray(q))].sum()) for nm, q in BA.items() if len(q)}
in_plane = {"sym", "side_far"}
check("(b) 外部境界・壁の節点は動かない (動くのは z 一定の面 sym・side_far の面内だけ)",
      all(n == 0 for nm, n in moved_tags.items() if nm not in in_plane), {nm: n for nm, n in moved_tags.items() if n})
mi, xi = closure(hB, BB)
check("(b) 境界の閉性", mi == 0 and xi == 0, f"missing {mi} extra {xi}")
_P = cB[hB]
_J = np.einsum("ij,ij->i", _P[:, 1] - _P[:, 0], np.cross(_P[:, 3] - _P[:, 0], _P[:, 4] - _P[:, 0]))
check("(b) ヘキサの符号付きヤコビアンが同符号・非退化", np.all(_J > 1e-18) or np.all(_J < -1e-18))
check("(b) float32 で新たな節点衝突が起きない",
      np.unique(cB, axis=0).shape[0] - np.unique(cB.astype(np.float32), axis=0).shape[0]
      == np.unique(cA, axis=0).shape[0] - np.unique(cA.astype(np.float32), axis=0).shape[0])
check("(b) info に変形の記録 (区間の終端・端点の傾き・station 数・最大移動量)",
      all(kk in iB for kk in ("te_wake_x_end", "te_wake_slope_te_deg", "te_wake_slope_end_deg", "te_wake_n_stations", "te_wake_max_dy"))
      and iB["te_wake_blend_H"] == LB)

# ------------------------------------------------------------------------------------------------ (c)
Lc, yte = float(iA["L_cowl"]), float(d.cowl_xy[-1, 1]); x1 = Lc + LB
m0 = np.tan(np.radians(iB["te_wake_slope_te_deg"])); m1 = np.tan(THB)
G = cB[:NB].reshape(ni, NJ, nz, 3) / H
lower_last = (G[it, jm, 0, 1] - G[it - 1, jm, 0, 1]) / (G[it, jm, 0, 0] - G[it - 1, jm, 0, 0])
# 上側壁面 (cowl_in) の最後の辺: 後縁の 1 つ上流の station のカウル上コピー (番号はメッシャの列挙と同じ規則で数える:
# i < i_te・k ≤ k_sw を順に、側壁より下流の側端 (k = k_sw, i ≥ i_sw) は上下で共有するので飛ばす)
_dup1, _nid = {}, NB
for _i in range(it):
    for _k in range(iA["k_sw"] + 1):
        if _k == iA["k_sw"] and _i >= iA["i_sw"] and P0.nz_out > 0 and P0.share_cowl_free_edge:
            continue
        _dup1[(_i, _k)] = _nid; _nid += 1
assert _nid - NB == iA["n_dup_cowl"]
_up = cB[_dup1[(it - 1, 0)]] / H
upper_last = (G[it, jm, 0, 1] - _up[1]) / (G[it, jm, 0, 0] - _up[0])
check("(c) 後縁の勾配 m0 = 下側壁面の最後の辺の傾き (座標から測った値と一致)", abs(m0 - lower_last) < 1e-9,
      f"m0 {np.degrees(np.arctan(m0)):.5f}° vs 辺 {np.degrees(np.arctan(lower_last)):.5f}°")
check("(c) 上側壁面の最後の辺とは別 (有限厚で傾きが違う)", np.isfinite(upper_last) and abs(m0 - upper_last) > 1e-3,
      f"上側 {np.degrees(np.arctan(upper_last)):.4f}°")
eps = 1e-7
yA = lambda x: float(ymA(x)); yB = lambda x: float(ymB(x))
check("(c) L_cowl で位置 = 後縁 (旧中間線と同じ)", yB(Lc) == yA(Lc) and abs(yB(Lc + 1e-12) - yte) < 1e-10)
check("(c) L_cowl の右側の勾配 = m0", abs((yB(Lc + eps) - yB(Lc)) / eps - m0) < 1e-5, f"{(yB(Lc + eps) - yB(Lc)) / eps:.6f} vs {m0:.6f}")
check("(c) L_cowl + L_b で位置 = 旧中間線", abs(yB(x1 - 1e-12) - yA(x1)) < 1e-10 and yB(x1) == yA(x1))
check("(c) L_cowl + L_b の左側の勾配 = 旧中間線の勾配 (tan θ_b)", abs((yB(x1) - yB(x1 - eps)) / eps - m1) < 1e-5,
      f"{(yB(x1) - yB(x1 - eps)) / eps:.6f} vs {m1:.6f}")
check("(c) 区間の外 (上流・下流) は旧中間線と同じ値", all(yB(x) == yA(x) for x in (-0.3, 0.0, 0.5, Lc, x1, x1 + 0.2, 3.0)))
t1 = (G[it + 1, jm, 0, 0] - G[it, jm, 0, 0]) / LB
chord = (G[it + 1, jm, 0, 1] - G[it, jm, 0, 1]) / (G[it + 1, jm, 0, 0] - G[it, jm, 0, 0])
check("(c) 実際の最初の下流の辺の傾き = 弦の式 m0 + (m1 − m0)(2t₁ − t₁²)", abs(chord - (m0 + (m1 - m0) * (2 * t1 - t1 * t1))) < 1e-9,
      f"辺 {np.degrees(np.arctan(chord)):.4f}°、接線 {np.degrees(np.arctan(m0)):.4f}°、t₁ {t1:.4f}")

# ------------------------------------------------------------------------------------------------ (d)
check("(d) 負の L_b は失敗", raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=-0.1))))
check("(d) 出口に届く L_b は失敗", raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=float(iA["x_out"]) - Lc))))
check("(d) 区間に station が無い L_b は失敗 (変形が効かない)", raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=0.01))))

# ------------------------------------------------------------------------------------------------ (e) runner の受け渡し
import yaml  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402


class _Stop(Exception):
    pass


_seen = {}
_saved = (R2.design_from_problem, R3.generate_sern_mesh3d)
R2.design_from_problem = lambda p, design=None: (None, d, {}, THB)        # MOC を回さない (受け渡しだけ見る)


def _capture(design, prm):
    _seen["prm"] = prm
    raise _Stop


R3.generate_sern_mesh3d = _capture
try:
    src = yaml.safe_load(open(ROOT / "case" / "46.sern_design" / "problem_3d_prod_3op_wallres_lswx08.yaml"))
    for label, extra, want in (("キー無し", None, 0.0), ("0.3", 0.3, 0.3)):
        y = json.loads(json.dumps(src))
        if extra is not None:
            y["mesh3d"]["te_wake_blend_H"] = extra
        pth = tmp / f"prob_{label}.yaml"; pth.write_text(yaml.safe_dump(y, allow_unicode=True))
        _seen.clear()
        try:
            R3.prepare(pth, tmp / f"run_{label}", op="m10_on")
        except _Stop:
            pass
        check(f"(e) runner: {label} → te_wake_blend_H {want}", "prm" in _seen and _seen["prm"].te_wake_blend_H == want,
              str(getattr(_seen.get("prm"), "te_wake_blend_H", None)))
finally:
    R2.design_from_problem, R3.generate_sern_mesh3d = _saved

# ------------------------------------------------------------------------------------------------ (f) 格子の確認の道具
TC = load_module(CHECK, "te_wake_grid_check")
gA, gB = tmp / "gA", tmp / "gB"
for g, c, h, B, info in ((gA, cA, hA, BA, iA), (gB, cB, hB, BB, iB)):
    g.mkdir()
    write_msh41_3d(g / "sern.msh", c, h, B, PHYS_SERN3D)
    (g / "prepare_info.json").write_text(json.dumps({"H_m": H, "mesh": info}))
A_, B_ = TC.load(str(gA)), TC.load(str(gB))
check("(f) msh の読み込み: 節点数・ヘキサが生成と一致", A_["coords"].shape == cA.shape and np.array_equal(A_["hexes"], hA))
S = TC.Struct(iA)
idn = TC.identity(A_, B_, S, iB)
check("(f) (v) 同一性: x・z・接続・境界面が同一", idn["x_equal"] and idn["z_equal"] and idn["hex_equal"] and idn["bface_equal"])
# msh は 10 桁で書くので、上線の直下で丸め程度 (~1e-17 m) しか動かない節点は「動いていない」に入る (float32 の forge 入力も同じ)
_lost = np.setdiff1d(mv, idn["moved"])
check("(f) (v) 動いた節点の範囲が変形区間の station・主ブロック (取りこぼしは丸め程度の移動だけ)",
      idn["moved_i"] == (iB["te_wake_i_first"], iB["te_wake_i_last"]) and idn["moved_nonbase"] == 0
      and np.all(np.isin(idn["moved"], mv)) and (_lost.size == 0 or np.abs(cA[_lost, 1] - cB[_lost, 1]).max() < 1e-12),
      f"i {idn['moved_i']}、道具 {idn['moved'].size} / 生成 {mv.size}、丸め程度 {_lost.size}")
check("(f) (v) 外部境界・壁の節点は動いていない (sym・side_far の面内は可)",
      all(n == 0 for p, n in idn["moved_by_tag"].items() if p not in TC.IN_PLANE_TAGS))
check("(f) カウル上コピーの番号の再現が座標と合う", all(np.allclose(A_["coords"][S.dup1[(i_, k_)], [0, 2]], A_["coords"][S.base(i_, jm, k_), [0, 2]])
                                                  for (i_, k_) in list(S.dup1)[:50]))
rA, rB = TC.te_column(A_, S, 5), TC.te_column(B_, S, 5)
key = "後縁節点 j=jm  z≤W/2  上流=下側壁面 (cowl_out)"
check("(f) (i) A の最初の下流の辺 = θ_b", abs(rA[key][0][0] - np.degrees(THB)) < 1e-6 and abs(rA[key][0][2] - np.degrees(THB)) < 1e-6,
      f"{rA[key][0]}")
check("(f) (i) B の最初の下流の辺 = 弦の角度 (全 z で同じ)", abs(rB[key][0][0] - np.degrees(np.arctan(chord))) < 1e-6
      and abs(rB[key][0][2] - rB[key][0][0]) < 1e-9, f"{rB[key][0]} vs {np.degrees(np.arctan(chord)):.4f}")
check("(f) (i) 上流の辺 (下側壁面) の中央値 = m0", abs(rA[key][1][1] - iB["te_wake_slope_te_deg"]) < 1e-6, f"{rA[key][1]}")
tA, tB = TC.turn_maxima(A_, S, iB["te_wake_i_last"] + 1), TC.turn_maxima(B_, S, iB["te_wake_i_last"] + 1)
kte = ("中間線/下側壁面 (j=jm)", "後縁 (i=i_te)")
check("(f) (ii) 後縁の折れ角: B < A", tB[kte][0] < tA[kte][0], f"A {tA[kte][0]:.2f}° B {tB[kte][0]:.2f}°")
koth = ("下側 (j<jm)", "それ以外")
check("(f) (ii) 変形区間の外の折れ角の最大は A と B で同じ", abs(tA[koth][0] - tB[koth][0]) < 1e-9)
JA, _ = TC.corner_jacobians(A_["coords"], A_["hexes"]); JB, _ = TC.corner_jacobians(B_["coords"], B_["hexes"])
sg = 1.0 if np.median(JA) > 0 else -1.0
check("(f) (iv) A・B とも全ヘキサの 8 頂点のヤコビアンが正", np.all(sg * JA > 0) and np.all(sg * JB > 0))
r = subprocess.run([sys.executable, str(CHECK), str(gA), str(gB), "--n-layers", "5"], capture_output=True, text=True)
check("(f) 道具の CLI が完走し、品質・双対幾何の回し方を出す", r.returncode == 0 and "check_dual_closure.py" in r.stdout and "MESH_QUALITY.txt" in r.stdout,
      r.stderr[-300:], only_on_fail=True)
r = subprocess.run([sys.executable, str(CHECK), str(gA), str(gA)], capture_output=True, text=True)
check("(f) A と A: 動いた節点 0", r.returncode == 0 and "y が動いた節点        : 0" in r.stdout, r.stdout[-300:] + r.stderr[-300:], only_on_fail=True)

# ------------------------------------------------------------------------------------------------ (g) 初期場の道具
CONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]
xdmf = lambda h: np.column_stack([np.full(len(h), 9, np.int64), h]).ravel().astype(np.int32)


def write_mesh_h5(path, coords, hexes, B, ic=1.0, wall=0.5, aux=True):
    """合成の forge 入力 h5: MESH/COORD (float32)・VIZMESH/CONNE・BCONDS/*/vizBface*・VALUE (一様 IC + wall_dist)・/AUX。"""
    n = len(coords)
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float32).ravel()
        f["VIZMESH/CONNE"] = xdmf(hexes)
        for nm, q in B.items():
            if len(q):
                q = np.asarray(q, np.int32)
                f[f"BCONDS/{PHYS_SERN3D[nm]}/vizBfaceSizes"] = np.full(len(q), 4, np.int32)
                f[f"BCONDS/{PHYS_SERN3D[nm]}/vizBfaceNodes"] = q.ravel()
        for kk in CONS:
            f["VALUE/" + kk] = np.full(n, ic, np.float32)
        f["VALUE/wall_dist"] = np.full(n, wall, np.float32)
        if aux:
            f["AUX/w_recon_vel"] = np.linspace(0.0, 1.0, n)


def write_res_h5(path, coords, hexes):
    """合成の res: 座標の滑らかな関数の原始量と、それと**わずかにずらした**保存量 (forge の res は原始量が更新前・
    密度が更新後なので ρ·原始量 は保存量と一致しない。不動節点が保存量を、動いた節点が原始量を使うことを見分ける)。"""
    n = len(coords); x, y, z = (coords.astype(np.float32).astype(np.float64).T / H)
    ro = (1.0 + 0.3 * x + 0.2 * y + 0.1 * z).astype(np.float32)
    prim = {"Ux": 100.0 + 50 * x - 30 * y, "Uy": -20.0 + 10 * y, "Uz": 5.0 * z, "P": 1e5 * (1 + 0.1 * x), "T": 300 + 10 * y,
            "k": 1.0 + 0.1 * x, "omega": 1e3 * (1 + y * y), "Y0": 0.3 + 0.1 * np.tanh(y), "Y1": 0.7 - 0.1 * np.tanh(y)}
    prim = {kk: np.asarray(v, np.float32) for kk, v in prim.items()}
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float32).ravel()
        f["MESH/CONNE"] = xdmf(hexes)
        f["VALUE/ro"] = ro
        for c, p_ in (("roUx", "Ux"), ("roUy", "Uy"), ("roUz", "Uz"), ("roK", "k"), ("roOmega", "omega"), ("roY0", "Y0"), ("roY1", "Y1")):
            f["VALUE/" + c] = (ro * prim[p_] * np.float32(1.0 + 1e-4)).astype(np.float32)
        f["VALUE/roe"] = (prim["P"] / 0.4 + 0.5 * ro * prim["Ux"] ** 2).astype(np.float32)
        for kk, v in prim.items():
            f["VALUE/" + kk] = v
        f["VALUE/wall_dist"] = np.full(n, 9.0, np.float32)


def run_rfd(*args):
    return subprocess.run([sys.executable, str(RFD), *map(str, args), "--force-species"], capture_output=True, text=True)


rd = tmp / "rfd"; rd.mkdir()
src_res, src_mesh, dst = rd / "res_A.h5", rd / "mesh_A.h5", rd / "mesh_B.h5"
write_res_h5(src_res, cA, hA); write_mesh_h5(src_mesh, cA, hA, BA); write_mesh_h5(dst, cB, hB, BB)
dst_interp = rd / "mesh_B_interp.h5"; shutil.copy(dst, dst_interp)
r = run_rfd(src_res, dst, "--src-mesh", src_mesh, "--list-moved", rd / "moved.csv")
check("(g) 完走して VERDICT: OK", r.returncode == 0 and "VERDICT: OK" in r.stdout, (r.stdout + r.stderr)[-600:], only_on_fail=True)
moved_f32 = np.where(np.any(cA.astype(np.float32) != cB.astype(np.float32), axis=1))[0]
check("(g) 動いた節点の数を出す (float32 の座標で数えた値と一致)", f"**動いた {moved_f32.size}**" in r.stdout, f"{moved_f32.size}")
ri = subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(src_res), str(dst_interp), "--force-species"],
                    capture_output=True, text=True)
check("(g) 比較用の interp_field.py が完走", ri.returncode == 0, ri.stderr[-300:], only_on_fail=True)
still = np.setdiff1d(np.arange(len(cA)), moved_f32)
with h5py.File(src_res, "r") as s, h5py.File(dst, "r") as g_, h5py.File(dst_interp, "r") as gi:
    ok_still = all(np.array_equal(g_["VALUE/" + kk][...][still], s["VALUE/" + kk][...][still]) for kk in CONS)
    ok_moved = all(np.array_equal(g_["VALUE/" + kk][...][moved_f32], gi["VALUE/" + kk][...][moved_f32]) for kk in CONS)
    differs = not np.array_equal(g_["VALUE/roUx"][...][moved_f32], s["VALUE/roUx"][...][moved_f32])
    wd = np.array_equal(g_["VALUE/wall_dist"][...], np.full(len(cA), 0.5, np.float32))
    aux = np.array_equal(g_["AUX/w_recon_vel"][...], np.linspace(0.0, 1.0, len(cA)))
check("(g) 不動節点: 全保存量が SRC の保存量とビット一致 (index コピー)", ok_still)
check("(g) 動いた節点: 全保存量が interp_field.py の結果とビット一致 (原始量の最近傍補間)", ok_moved)
check("(g) 動いた節点は元の番号の保存量の index コピーではない", differs)
check("(g) wall_dist は DST のまま", wd)
check("(g) /AUX (処置の重み w) は DST のまま", aux)
lst = np.loadtxt(rd / "moved.csv", delimiter=",", skiprows=1)
check("(g) 動いた節点の一覧 CSV (id・位置・移動量・最寄り)", lst.shape[0] == moved_f32.size and np.array_equal(lst[:, 0].astype(int), moved_f32))
# 全節点不動なら restart_field.py と同じ
dst_same, dst_same_rf = rd / "mesh_A_rfd.h5", rd / "mesh_A_rf.h5"
write_mesh_h5(dst_same, cA, hA, BA); shutil.copy(dst_same, dst_same_rf)
r1 = run_rfd(src_res, dst_same)
r2 = subprocess.run([sys.executable, str(TOOLS / "restart_field.py"), str(src_res), str(dst_same_rf), "--force-species"],
                    capture_output=True, text=True)
with h5py.File(dst_same, "r") as g1, h5py.File(dst_same_rf, "r") as g2:
    same = all(np.array_equal(g1["VALUE/" + kk][...], g2["VALUE/" + kk][...]) for kk in CONS + ["wall_dist"])
check("(g) 座標が全節点で一致: restart_field.py と同じ結果 (と案内)", r1.returncode == 0 and r2.returncode == 0 and same
      and "restart_field.py を使うのが本来" in r1.stdout, (r1.stdout + r2.stdout + r2.stderr)[-400:], only_on_fail=True)
# 拒否: 接続違い・節点数違い・境界面違い・双子への最寄り
bad_conn = rd / "bad_conn.h5"; write_mesh_h5(bad_conn, cB, hB[:, [1, 2, 3, 0, 5, 6, 7, 4]], BB)
r = run_rfd(src_res, bad_conn)
check("(g) 接続が違えば拒否 (interp_field を案内)", r.returncode != 0 and "interp_field" in (r.stdout + r.stderr))
bad_n = rd / "bad_n.h5"; write_mesh_h5(bad_n, cB[:-1], hB[:10], {})
r = run_rfd(src_res, bad_n)
check("(g) 節点数が違えば拒否", r.returncode != 0 and "節点数" in (r.stdout + r.stderr))
BB2 = dict(BB); BB2["cowl_in"] = [tuple(q[::-1]) for q in BB["cowl_in"]]
bad_bf = rd / "bad_bf.h5"; write_mesh_h5(bad_bf, cB, hB, BB2)
r = run_rfd(src_res, bad_bf, "--src-mesh", src_mesh)
check("(g) 境界面の接続が --src-mesh と違えば拒否", r.returncode != 0 and "境界面" in (r.stdout + r.stderr))
# 双子: 厚さ 0 の格子 (カウル上流が座標一致の双子) で、動いた節点 1 つを双子の座標の極近くへ置く
p0t = replace(P0, cowl_thickness=0.0)
cA0, hA0, BA0, iA0, _ = generate_sern_mesh3d(d, p0t)
cB0, _, _, _, _ = generate_sern_mesh3d(d, replace(p0t, te_wake_blend_H=LB))
twin_id = next(iter(range(ni * NJ * nz, len(cA0))))            # 最初のカウル上コピー (= 下側の base と座標一致)
assert np.unique(cA0[np.all(cA0 == cA0[twin_id], axis=1)], axis=0).shape[0] == 1
cB0t = cB0.copy(); mv0 = np.where(cA0[:, 1] != cB0[:, 1])[0]
cB0t[mv0[0]] = cA0[twin_id] + np.array([1e-9, 0.0, 0.0])
sr0, dt0 = rd / "res_A0.h5", rd / "twin.h5"
write_res_h5(sr0, cA0, hA0); write_mesh_h5(dt0, cB0t, hA0, BA0)
r = run_rfd(sr0, dt0)
check("(g) 動いた節点の最寄りが座標一致の双子なら拒否", r.returncode != 0 and "双子" in (r.stdout + r.stderr), (r.stdout + r.stderr)[-300:], only_on_fail=True)
nos = rd / "mesh_B_nospecies.h5"; write_mesh_h5(nos, cB, hB, BB)
r = subprocess.run([sys.executable, str(RFD), str(src_res), str(nos)], capture_output=True, text=True,
                   env={kk: v for kk, v in os.environ.items() if kk != "FORGE_ALLOW_UNVERIFIED_SPECIES"})
with h5py.File(nos, "r") as g_:
    untouched_nos = np.array_equal(g_["VALUE/ro"][...], np.ones(len(cB), np.float32))
check("(g) 化学種を照合できない (未検証・設定なし) ときは既定で拒否し何も書かない (restart_field と同じ)",
      r.returncode != 0 and "REFUSED" in (r.stdout + r.stderr) and untouched_nos, (r.stdout + r.stderr)[-300:], only_on_fail=True)
dry = rd / "mesh_B_dry.h5"; write_mesh_h5(dry, cB, hB, BB)
r = run_rfd(src_res, dry, "--dry-run")
with h5py.File(dry, "r") as g_:
    untouched = np.array_equal(g_["VALUE/ro"][...], np.ones(len(cB), np.float32))
check("(g) --dry-run は書かない", r.returncode == 0 and untouched)

# ------------------------------------------------------------------------------------------------ (h) 曲線版の定数
check("(h) 有効時の info に曲線版 te_wake_curve_version (= mesh_sern3d.TE_WAKE_CURVE_VERSION)",
      iB.get("te_wake_curve_version") == TE_WAKE_CURVE_VERSION, str(iB.get("te_wake_curve_version")))
check("(h) 既定 0 の info に曲線版が入らない (既定の格子と info のビット一致・同一は (a))", "te_wake_curve_version" not in iA)


# ------------------------------------------------------------------------------------------------ (i)–(l) --admission・来歴
GEOM_ROWS = ("最初の下流の辺", "後縁の折れ", "復帰区間の折れ", "固定すべき量", "層の Δy", "primal の品質")


def _row_status(R, row):
    return R.status_of(row)


def _item(R, row, prefix):
    hit = [r for r in R.rows if r["row"] == row and r["item"].startswith(prefix)]
    return hit[0] if len(hit) == 1 else None


def _admission_tests():
    import hashlib
    import yaml
    from forge_design.evaluate import runner_sern as R2
    from forge_design.evaluate import runner_sern3d as R3
    from forge_design.probdef import load_problem
    TC = load_module(CHECK, "te_wake_grid_check_adm")
    MZ = load_module(TOOLS / "mark_zero_thickness_edges.py", "mark_zte_for_te_wake_tests")
    prob_g1 = ROOT / "case" / "46.sern_design" / "problem_3d_prod_m6on_g1.yaml"
    p = load_problem(prob_g1)
    cached = R2.design_from_problem(p)                   # 生産の設計 MOC (~10 s)。以下の prepare はこれを使い回す
    Hg = float(p.spec["H_m"])

    class _Stop(Exception):
        pass
    seen = {}

    def _cap(design, prm):
        seen["prm"] = prm
        raise _Stop
    saved = (R2.design_from_problem, R3.generate_sern_mesh3d)
    R2.design_from_problem = lambda p_, design=None: cached
    try:
        # ---------------------------------------------------------------- (i) メモリ上の g1
        R3.generate_sern_mesh3d = _cap
        try:
            R3.prepare(prob_g1, tmp / "cap_g1")
        except _Stop:
            pass
        R3.generate_sern_mesh3d = saved[1]
        prm, des = seen["prm"], cached[1]
        c0, h0, B0, i0, _ = generate_sern_mesh3d(des, prm)
        c1, h1, B1, i1, _ = generate_sern_mesh3d(des, replace(prm, te_wake_blend_H=1.0))
        f32 = lambda c: c.astype(np.float32).astype(np.float64)     # 変換後の h5 と同じ座標の表現 (float32)
        gA = TC.from_mesher(f32(c0), h0, B0, i0, Hg, "g1 A"); gB = TC.from_mesher(f32(c1), h1, B1, i1, Hg, "g1 B (L_b 1.0 H)")
        gA["coord_repr"] = gB["coord_repr"] = "float32"
        v, R = TC.admission(gA, gB)
        print(f"--- (i) g1 (メモリ上): 節点 {i0['nodes']}、ヘキサ {len(h0)}")
        for row in GEOM_ROWS:
            check(f"(i) g1: {row} PASS", _row_status(R, row) == TC.PASS,
                  "; ".join(f"{r['item']} {r['status']} {r['value']}" for r in R.rows if r["row"] == row), only_on_fail=True)
        check("(i) g1: 前提の「入力が最終 node 格子」はメモリ上なので判定不能",
              [r["status"] for r in R.rows if r["row"] == "前提" and "最終 node 格子" in r["item"]] == [TC.UND, TC.UND])
        check("(i) g1: 品質・双対 (MESH_QUALITY・check_dual_closure の出力が無い) は判定不能",
              all(r["status"] == TC.UND for r in R.rows if r["row"] == "品質・双対の判定"))
        check("(i) g1: 総合は UNDECIDABLE (判定不能を合格にしない)", v == TC.UND, TC.VERDICT_WORD[v])
        te = max(r["num"] for r in R.rows if r["row"] == "後縁の折れ")
        ret = _item(R, "復帰区間の折れ", "主ブロック")["num"]
        fe = _item(R, "最初の下流の辺", "下側")["num"]
        jmin = _item(R, "primal の品質", "スケール済み")["num"]
        sk = _item(R, "primal の品質", "skew 最大 (B")
        check("(i) g1: 後縁の折れ・最初の下流の辺の m0 からの差 = codex の g3 再生成値 2.1975° (x station が g3 と同じ)",
              abs(te - 2.1975) < 1e-3 and abs(fe - 2.1975) < 1e-3, f"後縁 {te:.5f}° 最初の辺 {fe:.5f}°")
        check("(i) g1: 復帰区間の最大の折れ = codex の再生成値 5.08°", abs(ret - 5.0795) < 2e-3, f"{ret:.5f}°")
        # codex の再生成は倍精度の座標。float32 に丸めると A の skew 最大は 0.44927 になる (B は 0.53048 のまま) ので、照合は倍精度で
        _, jA64, _, skA64 = TC.hex_quality(c0, h0, 1.0)
        _, jB64, _, skB64 = TC.hex_quality(c1, h1, 1.0)
        check("(i) g1 (倍精度の座標): スケール済みヤコビアンの最小 B 0.67245・skew 最大 0.44706 → 0.53048 (codex の g3 再生成値と同じ)",
              abs(jB64 - 0.67245) < 5e-5 and abs(float(skA64.max()) - 0.44706) < 5e-5 and abs(float(skB64.max()) - 0.53048) < 5e-5,
              f"J A {jA64:.5f} B {jB64:.5f} skew A {skA64.max():.5f} B {skB64.max():.5f}")
        check("(i) g1 (float32): 判定に使う B のスケール済みヤコビアン・skew 最大も同じ値、skew の増加 ≤ 0.10",
              abs(jmin - 0.67245) < 5e-5 and abs(sk["num"] - 0.53048) < 5e-5 and _item(R, "primal の品質", "skew 最大の A")["num"] <= 0.10,
              f"J {jmin:.5f} skew A {sk['ref_A']} B {sk['num']:.5f}")
        del c0, c1, h0, h1, gA, gB, skA64, skB64

        # ---------------------------------------------------------------- (j) 変換した粗い格子 (runner の prepare)
        def prep(lab, lb):
            y = yaml.safe_load(open(prob_g1))
            # g1 の x・y 分布 (ni_plume 153・nj_top 49・nj_bot 31・first_wall_frac) のまま、z と上流のノズル区間を粗くする。
            # 基部の後流の第一間隔・z の第一間隔・外部流バンドは品質ゲート (AR ≤ 5000) を通す粗さに
            y["mesh3d"].update({"nz_in": 4, "nz_out": 3, "ni_noz": 30, "nj_vside": 5, "first_wake_frac": 0.004, "first_z_frac": 0.02})
            y["mesh"].update({"nj_ext_top": 9, "first_top_frac": 0.004})
            if lb is not None:
                y["mesh3d"]["te_wake_blend_H"] = lb
            pth = tmp / f"adm_{lab}.yaml"
            pth.write_text(yaml.safe_dump(y, allow_unicode=True))
            R3.prepare(pth, tmp / f"adm_run_{lab}")
            return tmp / f"adm_run_{lab}"
        ra, rb, rb3 = prep("A", None), prep("B", 1.0), prep("B03", 0.3)

        def closure(run, name):
            r_ = subprocess.run([sys.executable, str(TOOLS / "check_dual_closure.py"), str(run / "sern.h5"), "--tol", "1e-5"],
                                capture_output=True, text=True)
            (tmp / name).write_text(r_.stdout + r_.stderr)
            return tmp / name
        clo_b, clo_b3 = closure(rb, "closure_B.txt"), closure(rb3, "closure_B03.txt")

        def cli(*args):
            return subprocess.run([sys.executable, str(CHECK), *map(str, args)], capture_output=True, text=True)
        r = cli(ra, rb, "--admission", "--dual-closure-b", clo_b)
        print(f"--- (j) 変換した粗い格子: {json.loads((rb / 'prepare_info.json').read_text())['mesh']['nodes']} 節点")
        check("(j) run ディレクトリ (最終 node の h5・MESH_QUALITY・info) + 双対の出力: ADMISSION VERDICT: PASS、終了コード 0",
              r.returncode == 0 and r.stdout.strip().endswith("ADMISSION VERDICT: PASS"), (r.stdout + r.stderr)[-1500:], only_on_fail=True)
        check("(j) 行ごとに測定値と閾値を出す", "≤ 3°" in r.stdout and "≤ 6°" in r.stdout and "≥ 0.65" in r.stdout and "B ≥ A" in r.stdout)
        r = cli(ra, rb, "--admission")
        check("(j) check_dual_closure の出力を渡さない: UNDECIDABLE、終了コード 2",
              r.returncode == 2 and "ADMISSION VERDICT: UNDECIDABLE" in r.stdout, r.stdout[-400:], only_on_fail=True)
        r = cli(ra, rb3, "--admission", "--dual-closure-b", clo_b3, "--expect-blend-H", "0.3")
        check("(j) L_b 0.3 H の B (事前登録の値を 0.3 にしても): 後縁の折れ・最初の下流の辺・復帰区間の折れが FAIL、終了コード 1",
              r.returncode == 1 and "ADMISSION VERDICT: FAIL" in r.stdout
              and all(f"{row} FAIL" in r.stdout for row in ("最初の下流の辺", "後縁の折れ", "復帰区間の折れ")), r.stdout[-600:], only_on_fail=True)
        r = cli(ra, rb3, "--admission", "--dual-closure-b", clo_b3)
        check("(j) L_b が事前登録の 1.0 H でない B は前提が FAIL", r.returncode == 1 and "前提 FAIL" in r.stdout, r.stdout[-300:], only_on_fail=True)
        r = cli(ra, rb, "--admission", "--dual-closure-b", clo_b, "--info-b", tmp / "no_such_info.json")
        check("(j) info が無い: UNDECIDABLE、終了コード 2", r.returncode == 2 and "ADMISSION VERDICT: UNDECIDABLE" in r.stdout,
              r.stdout[-300:] + r.stderr[-300:], only_on_fail=True)
        r = cli(ra, rb, "--admission", "--dual-closure-b", clo_b, "--quality-b", tmp / "no_such_quality.txt")
        check("(j) B の MESH_QUALITY が無い: UNDECIDABLE", r.returncode == 2 and "ADMISSION VERDICT: UNDECIDABLE" in r.stdout,
              r.stdout[-300:], only_on_fail=True)
        r = cli(ra / "sern.msh", rb / "sern.msh", "--info-a", ra / "prepare_info.json", "--info-b", rb / "prepare_info.json", "--admission",
                "--quality-a", ra / "MESH_QUALITY.txt", "--quality-b", rb / "MESH_QUALITY.txt", "--dual-closure-b", clo_b)
        check("(j) msh (変換前) を渡す: 前提・層の Δy が判定不能で UNDECIDABLE (合格にしない)",
              r.returncode == 2 and "ADMISSION VERDICT: UNDECIDABLE" in r.stdout and "前提 判定不能" in r.stdout and "層の Δy 判定不能" in r.stdout,
              r.stdout[-500:], only_on_fail=True)

        # ---------------------------------------------------------------- (k) 壊した格子・入力の欠落
        A_, B_ = TC.load(str(ra)), TC.load(str(rb))
        qa, qb = str(ra / "MESH_QUALITY.txt"), str(rb / "MESH_QUALITY.txt")

        def adm(Bx, Ax=None, **kw):
            kw.setdefault("quality_a", qa); kw.setdefault("quality_b", qb); kw.setdefault("dual_b", str(clo_b))
            return TC.admission(Ax or A_, Bx, **kw)
        v0, R0 = adm(B_)
        check("(k) 基準 (変換した B そのまま、Python から): PASS", v0 == TC.PASS, TC.VERDICT_WORD[v0])
        S_ = TC.Struct(B_["info"]); jm_, NJ_ = S_.jm, S_.NJ
        i1_ = int(B_["info"]["te_wake_i_first"])
        cB_ = B_["coords"]

        def mut(fn, hexes=False):
            Bx = dict(B_)
            Bx["coords"] = cB_.copy()
            if hexes:
                Bx["hexes"] = B_["hexes"].copy()
            fn(Bx)
            Bx["coords"] = Bx["coords"].astype(np.float32).astype(np.float64)      # 変換後の表現のまま
            return Bx

        def statuses(R):
            return {row: R.status_of(row) for row in dict.fromkeys(r["row"] for r in R.rows)}

        # 復帰区間の折れ: 上側の第 20 層 (判定する 10 層の外)・k=1 (sym・side_far でない) の 1 節点を y にずらす
        nk = S_.base(i1_ + 3, jm_ + 20, 1)
        dx_ = cB_[S_.base(i1_ + 4, jm_ + 20, 1), 0] - cB_[nk, 0]
        dy_ = min(abs(cB_[S_.base(i1_ + 3, jm_ + 21, 1), 1] - cB_[nk, 1]), abs(cB_[nk, 1] - cB_[S_.base(i1_ + 3, jm_ + 19, 1), 1]))
        shift = min(0.15 * dx_, 0.3 * dy_)

        def f_kink(Bx):
            Bx["coords"][nk, 1] += shift
        v, R = adm(mut(f_kink)); st = statuses(R)
        check("(k) 復帰区間に折れを作る: 復帰区間の折れ FAIL、最初の下流の辺・後縁の折れ・固定すべき量・層の Δy は PASS",
              st["復帰区間の折れ"] == TC.FAIL and all(st[x] == TC.PASS for x in ("最初の下流の辺", "後縁の折れ", "固定すべき量", "層の Δy"))
              and v == TC.FAIL, f"{st} 折れ {_item(R, '復帰区間の折れ', '主ブロック')['num']:.3f}° (ずらし {shift:.3e} m)")

        def f_x(Bx):
            Bx["coords"][nk, 0] += 1e-6
        v, R = adm(mut(f_x))
        check("(k) x を 1 節点変える: 固定すべき量 (x) FAIL", _item(R, "固定すべき量", "x 座標")["status"] == TC.FAIL and v == TC.FAIL)
        nu = S_.base(S_.i_te - 5, jm_ - 20, 1)

        def f_up(Bx):
            Bx["coords"][nu, 1] += 1e-5
        v, R = adm(mut(f_up))
        check("(k) 上流の節点を動かす: 固定すべき量 (上流の座標) FAIL", _item(R, "固定すべき量", "上流の座標")["status"] == TC.FAIL and v == TC.FAIL)
        nr = S_.base(i1_ + 3, NJ_ - 1, 1)            # 上線 (j = NJ−1)、k ≤ k_sw = ramp

        def f_ramp(Bx):
            Bx["coords"][nr, 1] -= 1e-5
        v, R = adm(mut(f_ramp))
        check("(k) 復帰区間の壁 (ramp) の節点を動かす: 固体形状 FAIL", _item(R, "固定すべき量", "固体形状")["status"] == TC.FAIL and v == TC.FAIL,
              _item(R, "固定すべき量", "固体形状")["value"])
        nb_ = S_.base(i1_ + 3, 0, 1)                 # 下端 (j = 0) = bottom

        def f_bot(Bx):
            Bx["coords"][nb_, 1] += 1e-5
        v, R = adm(mut(f_bot))
        check("(k) 外部境界 (bottom) の節点を動かす: 外部境界の形状 FAIL", _item(R, "固定すべき量", "外部境界")["status"] == TC.FAIL and v == TC.FAIL)
        check("(k) 基準で sym・side_far の面内の移動を数として記録 (FAIL にしない)",
              "sym" in _item(R0, "固定すべき量", "外部境界")["note"] and _item(R0, "固定すべき量", "外部境界")["status"] == TC.PASS)

        def f_hex(Bx):
            Bx["hexes"][0] = Bx["hexes"][0][[1, 2, 3, 0, 5, 6, 7, 4]]
        v, R = adm(mut(f_hex, hexes=True))
        check("(k) ヘキサの接続を 1 つ変える: 固定すべき量 (接続) FAIL", _item(R, "固定すべき量", "ヘキサ")["status"] == TC.FAIL and v == TC.FAIL)
        n1 = S_.base(i1_ + 1, jm_ - 1, 2)
        dy1 = cB_[S_.base(i1_ + 1, jm_, 2), 1] - cB_[n1, 1]

        def f_l1(Bx):
            Bx["coords"][n1, 1] += 1e-3 * dy1
        v, R = adm(mut(f_l1))
        check("(k) 復帰区間の下側の第 1 層を 1e-3 だけ変える: 第 1 層 FAIL・第 1〜第 10 層 PASS",
              _item(R, "層の Δy", "下側の第 1 層")["status"] == TC.FAIL and _item(R, "層の Δy", "下側の第 1〜")["status"] == TC.PASS and v == TC.FAIL,
              f"{_item(R, '層の Δy', '下側の第 1 層')['value']}")
        n9 = S_.base(i1_ + 2, jm_ - 12, 1)

        def f_inv(Bx):
            Bx["coords"][n9, 1] = cB_[S_.base(i1_ + 2, jm_ - 14, 1), 1] - 1e-6
        v, R = adm(mut(f_inv))
        check("(k) ヘキサを反転させる: 非正の頂点ヤコビアン・スケール済みの最小が FAIL",
              _item(R, "primal の品質", "非正")["status"] == TC.FAIL and _item(R, "primal の品質", "スケール済み")["status"] == TC.FAIL and v == TC.FAIL)

        def f_n(Bx):
            Bx["coords"] = np.vstack([Bx["coords"], Bx["coords"][-1:]])
        v, R = adm(mut(f_n))
        check("(k) 節点数が違う: 固定すべき量 FAIL、他の幾何の行は判定不能",
              R.status_of("固定すべき量") == TC.FAIL and R.status_of("復帰区間の折れ") == TC.UND and v == TC.FAIL)
        Bn = dict(B_); Bn["bfaces"] = None
        v, R = adm(Bn)
        check("(k) 境界面の接続が無い入力: 固定すべき量が判定不能、総合 UNDECIDABLE",
              R.status_of("固定すべき量") == TC.UND and v == TC.UND, str(statuses(R)))
        Bm = dict(B_); Bm["source"], Bm["coord_repr"] = "msh", "msh10g"
        v, R = adm(Bm)
        check("(k) msh の座標 (丸めの上限が未登録): 層の Δy が判定不能", R.status_of("層の Δy") == TC.UND and v == TC.UND)

        # 品質・双対の判定 (合成した出力ファイル)
        MA, MB = A_["hexes"].shape[0], B_["hexes"].shape[0]

        def qfile(name, verdict, cells, ar=5000.0, sk=0.9):
            pth = tmp / name
            pth.write_text(f"mode: 3D (per-face skew / edge AR)\ncells: {cells}\nVERDICT: {verdict} (AR<= {ar:.0f}, skew<= {sk:.2f})\n")
            return str(pth)
        v, R = adm(B_, quality_a=qfile("qa_pass", "PASS", MA), quality_b=qfile("qb_soft", "SOFT-PASS (<0.1% outliers)", MB))
        check("(k) 品質の VERDICT が A PASS → B SOFT-PASS: (a) FAIL", _item(R, "品質・双対の判定", "(a)")["status"] == TC.FAIL and v == TC.FAIL)
        v, R = adm(B_, quality_a=qfile("qa_soft", "SOFT-PASS (<0.1% outliers)", MA), quality_b=qfile("qb_pass", "PASS", MB))
        check("(k) 品質の VERDICT が A SOFT-PASS → B PASS: (a) PASS", _item(R, "品質・双対の判定", "(a)")["status"] == TC.PASS and v == TC.PASS)
        v, R = adm(B_, quality_a=qfile("qa_1000", "PASS", MA, ar=1000), quality_b=qfile("qb_5000", "PASS", MB))
        check("(k) 品質の閾値が A と B で違う: (a) 判定不能", _item(R, "品質・双対の判定", "(a)")["status"] == TC.UND and v == TC.UND)
        v, R = adm(B_, quality_a=qfile("qa_cells", "PASS", MA + 1), quality_b=qfile("qb_cells", "PASS", MB))
        check("(k) MESH_QUALITY のセル数が格子と違う (別の格子の出力): (a) 判定不能", _item(R, "品質・双対の判定", "(a)")["status"] == TC.UND)
        # 変形領域の AR: A より B で AR が上がるヘキサの間に閾値を置くと、新しい超過になる
        _, _, arsA, _ = TC.hex_quality(A_["coords"], A_["hexes"], 1.0)
        _, _, arsB, _ = TC.hex_quality(B_["coords"], B_["hexes"], 1.0)
        mvd = np.zeros(len(cB_), bool); mvd[np.where(A_["coords"][:, 1] != cB_[:, 1])[0]] = True
        touch = mvd[B_["hexes"]].any(axis=1)
        gain = np.where(touch & (arsB > arsA), arsB - arsA, -1.0)
        hx = int(np.argmax(gain))
        thr = 0.5 * (arsA[hx] + arsB[hx])
        v, R = adm(B_, quality_a=qfile("qa_thr", "PASS", MA, ar=thr), quality_b=qfile("qb_thr", "PASS", MB, ar=thr))
        it_b = _item(R, "品質・双対の判定", "(b)")
        check("(k) 変形領域で A は閾値以下・B は超過するヘキサがある閾値: (b) FAIL",
              gain[hx] > 0 and it_b["status"] == TC.FAIL and it_b["num"] >= 1, f"閾値 {thr:.3f}: {it_b['value']}")
        v, R = adm(B_, ar_max=4000.0)
        check("(k) --ar-max が MESH_QUALITY の AR<= と違う: (b) 判定不能", _item(R, "品質・双対の判定", "(b)")["status"] == TC.UND)
        v, R = adm(B_, quality_a=None, quality_b=None, ar_max=5000.0)
        check("(k) MESH_QUALITY が無く --ar-max だけ: (a) 判定不能・(b) は判定する (PASS)",
              _item(R, "品質・双対の判定", "(a)")["status"] == TC.UND and _item(R, "品質・双対の判定", "(b)")["status"] == TC.PASS and v == TC.UND)
        nB_ = len(cB_)

        def dfile(name, verdict, n_cv, tol="1e-05", clo=True):
            pth = tmp / name
            pth.write_text(f"CV {n_cv} / 双対面 10 (内部 8, 境界 2)\n体積: min 1e-12  max 1e-3  非正または非有限 0\n"
                           + (f"閉性 |ΣS|/Σ|S|: median 1e-09  p99 1e-08  max 1e-07   (> {tol}: 0 CV, 面の無い CV: 0)\n" if clo else "")
                           + f"VERDICT: {verdict}\n")
            return str(pth)
        v, R = adm(B_, dual_b=dfile("d_fail", "FAIL", nB_))
        check("(k) 双対の VERDICT FAIL: (c) FAIL", _item(R, "品質・双対の判定", "(c)")["status"] == TC.FAIL and v == TC.FAIL)
        v, R = adm(B_, dual_b=dfile("d_tol", "PASS", nB_, tol="1e-04"))
        check("(k) 双対の tol が 1e-4: (c) 判定不能", _item(R, "品質・双対の判定", "(c)")["status"] == TC.UND and v == TC.UND)
        v, R = adm(B_, dual_b=dfile("d_ncv", "PASS", nB_ + 7))
        check("(k) 双対の CV 数が B の節点数と違う: (c) 判定不能", _item(R, "品質・双対の判定", "(c)")["status"] == TC.UND and v == TC.UND)
        v, R = adm(B_, dual_b=dfile("d_und", "FAIL (判定不能: PLANES/STRUCT が無い)", nB_, clo=False))
        check("(k) check_dual_closure 自身の判定不能: (c) 判定不能", _item(R, "品質・双対の判定", "(c)")["status"] == TC.UND and v == TC.UND)
        v, R = adm(B_, dual_b=None)
        check("(k) 双対の出力が無い: (c) 判定不能", _item(R, "品質・双対の判定", "(c)")["status"] == TC.UND and v == TC.UND)

        # ---------------------------------------------------------------- (l) 来歴・--same-mesh
        prov = tmp / "prov.json"
        att = tmp / "ic_transfer.log"; att.write_text("restart_field_deformed: VERDICT: OK\n")
        r = cli(ra, rb, "--admission", "--dual-closure-b", clo_b, "--provenance-out", prov, "--code-commit-b", "0123abcd", "--attach-b", att)
        doc = json.loads(prov.read_text()) if prov.exists() else {}
        pa_, pb_ = doc.get("runs", {}).get("A", {}), doc.get("runs", {}).get("B", {})
        sigA, sigB = MZ.mesh_signature(str(ra / "sern.h5")), MZ.mesh_signature(str(rb / "sern.h5"))
        check("(l) 来歴: 格子署名 = mark_zero_thickness_edges.mesh_signature を実入力格子 (sern.h5) で再計算した値",
              pa_.get("mesh_signature") == sigA and pb_.get("mesh_signature") == sigB and pa_.get("mesh_signature_version") == MZ.SIG_VERSION)
        check("(l) 来歴: A と B の署名は違う (A_B_same_mesh_signature false)", sigA != sigB and doc.get("A_B_same_mesh_signature") is False)
        check("(l) 来歴: L_b (A 0・B 1.0) と曲線版 (B = 現行の定数、A は変形なし)",
              pa_.get("te_wake_blend_H") == 0.0 and pb_.get("te_wake_blend_H") == 1.0 and pb_.get("te_wake_curve_version") == TE_WAKE_CURVE_VERSION
              and pa_.get("te_wake_curve_version") is None and doc.get("mesher_curve_version_now") == TE_WAKE_CURVE_VERSION)
        check("(l) 来歴: commit は prepare_info に記録が無い → A は記録なし、B は申告を「未検証」として記録",
              pa_.get("code_commit") is None and pb_.get("code_commit") == "0123abcd" and "未検証" in pb_.get("code_commit_source", ""))
        sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
        check("(l) 来歴: 設定ファイル・問題 YAML・添付の sha256",
              pb_.get("config_sha256", {}).get("solverConfig.yaml") == sha(rb / "solverConfig.yaml")
              and pb_.get("config_sha256", {}).get("bcondConfig.yaml") == sha(rb / "bcondConfig.yaml")
              and pb_.get("config_sha256", {}).get("prepare_info.json") == sha(rb / "prepare_info.json")
              and pb_.get("problem_yaml", {}).get("sha256") == sha(tmp / "adm_B.yaml")
              and pb_.get("attached", [{}])[0].get("sha256") == sha(att))
        check("(l) 来歴: 判定の行と総合 (PASS) を含む、設計 DB に取り込まない旨", doc.get("admission", {}).get("verdict") == "PASS"
              and len(doc.get("admission", {}).get("rows", [])) > 20 and "取り込まない" in doc.get("design_db", ""))
        old_info = json.loads((ra / "prepare_info.json").read_text())
        old_info["mesh"] = {kk: vv for kk, vv in old_info["mesh"].items() if not kk.startswith("te_wake_")}
        (tmp / "old_info.json").write_text(json.dumps(old_info))
        e_old = TC.provenance_entry(str(ra), str(tmp / "old_info.json"))
        check("(l) 来歴: オプション導入前の info (te_wake_blend_H のキー無し) は 0 として出所を書く",
              e_old["te_wake_blend_H"] == 0.0 and "キー無し" in e_old["te_wake_blend_H_source"])
        cp = tmp / "A_copy_new_values.h5"
        shutil.copy(ra / "sern.h5", cp)
        with h5py.File(cp, "r+") as f:
            f["VALUE/ro"][...] = 2.0 * f["VALUE/ro"][...]          # 初期場の書き込み (restart_field 相当) では格子は変わらない
        r = cli("--same-mesh", ra, cp)
        check("(l) --same-mesh: run ディレクトリ A と /VALUE を書き換えた A の h5 → SAME MESH: YES、終了コード 0",
              r.returncode == 0 and "SAME MESH: YES" in r.stdout, r.stdout[-300:] + r.stderr[-300:], only_on_fail=True)
        r = cli("--same-mesh", ra / "sern.h5", rb / "sern.h5")
        check("(l) --same-mesh: A と B → SAME MESH: NO、終了コード 1", r.returncode == 1 and "SAME MESH: NO" in r.stdout)
        r = cli("--same-mesh", ra / "sern.h5", tmp / "no_such.h5")
        check("(l) --same-mesh: 読めない入力 → SAME MESH: UNDECIDABLE、終了コード 2", r.returncode == 2 and "SAME MESH: UNDECIDABLE" in r.stdout)
    finally:
        R2.design_from_problem, R3.generate_sern_mesh3d = saved


_admission_tests()

shutil.rmtree(tmp, ignore_errors=True)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
