#!/usr/bin/env python3
"""後縁下流の局所変形の変位の局所化 (`mesh3d.te_wake_mode: local`) の試験
(plan tooling-sern-te-wake-grid §5.1 #4c、事前登録 2026-10-11、codex diagnose g4-grid-redesign)。forge は回さない。

合格 (測る前に固定):
  (a) 既定 (キー無し)・te_wake_blend_H 0 の明示・0 + local・band の明示は、変更前のメッシャ (commit aa712de0 の mesh_sern3d.py) と
      座標・ヘキサ・境界面・返す中間線の関数がビット一致。info は L_b > 0 の band に te_wake_mode "band" が増えるだけ (L_b 0 は不変)
  (b) local の変位の式: 0 格子 (L_b 0) と x・z・接続・境界面がビット一致、動くのは変形区間の station の主ブロックだけ、
      各節点の y の差 = δ(x)·φ(d) (δ = 現方式の中間線 − 旧中間線、d = 0 格子での旧中間線からの鉛直距離 /H)、d ≤ 0.02 H で φ = 1 (= δ)、
      d ≥ 0.30 H でビット不変、間は φ = 1 − 10t³ + 15t⁴ − 6t⁵ (t = (d/H − 0.02)/0.28)、上下の帯の両方に同じ式が成り立つ、
      中間線は現方式の Hermite と丸め内で一致、下端・上線・壁・外部境界の節点が動かない、閉性・ヤコビアンの符号・float32 の衝突
  (c) φ の端点・遷移の値と単調性、定数が問題 YAML から変えられない (格子パラメータに無い)、不正な te_wake_mode・上線/下端が
      支持領域に入る指定は生成を失敗させる
  (d) info と識別: local の info に te_wake_mode・版 (= TE_WAKE_MODE_VERSION)・定数、評価方式の識別 (eval_method) に版が入り
      band と別の id、prepare の記録 (メッシャの info) とキャンペーンの要求 (問題 YAML) の id が一致、band の識別子の文字列は従来のまま
  (e) runner_sern3d の受け渡し: 問題 YAML の mesh3d.te_wake_mode (無ければ "band")
  (f) te_wake_hex_compare.compare / rows_4c: A = 0 格子・B = local で支持領域の外の不変が PASS、B = band では判定不能、
      支持領域の外を 1 節点動かした B は FAIL
  (g) 変換した粗い格子 (g1 の x・y 分布、run_sern_te_wake_tests (j) と同じ) を runner_sern3d.prepare で作り、prepare_info に mode・版、
      評価方式の id = 問題 YAML からの要求、CLI の te_wake_hex_compare が TE_WAKE_LOCAL VERDICT: PASS (変換器が無ければ FAIL として数える)
  --production [g3] [g4] (重い。g4 + g3 で約 4 分・最大 RSS 2.2 GB [2026-10-11 実測]): 問題 YAML (problem_3d_prod_3op_wallres_lswx08_tewake_<g>_B10.yaml) から
      0 格子・現方式・局所化をメモリ上で作り、座標をメッシャの .10g 出力 → float32 に丸めて te_wake_hex_compare と
      te_wake_grid_check.admission を掛ける。g4 は codex diagnose 2026-10-11 の予備値 (現方式 4,180 / 6,105 / 44,880、局所化 0 / 0 /
      40,700、変形領域の最小 0.651909、復帰の折れ 5.999851° → 5.170986°) を照合する。g3 は値を出すだけ (事前の値が無い)

  design/.venv-opt/bin/python design/tests/run_sern_te_wake_local_tests.py [--production g3 g4]
"""
import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIAG = ROOT / "case" / "46.sern_design" / "diag"
CASE = ROOT / "case" / "46.sern_design"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(DIAG))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
from forge_design.meshing import mesh_sern3d as MS  # noqa: E402
from forge_design.meshing.mesh_sern3d import (TE_WAKE_CURVE_VERSION, TE_WAKE_MODE_VERSION, SernMesh3DParams,  # noqa: E402
                                              generate_sern_mesh3d, te_wake_local_weight)

REF_COMMIT = "aa712de060df009f714b4a9c7608dec8e6aad4b7"   # 本変更の直前に mesh_sern3d.py を変えた commit
FAIL = 0
tmp = Path(tempfile.mkdtemp(prefix="te_wake_local_tests_"))


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
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{REF_COMMIT}:design/forge_design/meshing/mesh_sern3d.py"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    f = tmp / "_mesh_sern3d_ref_local.py"; f.write_text(r.stdout)
    spec = importlib.util.spec_from_file_location("forge_design.meshing._mesh_sern3d_ref_local", f)
    m = importlib.util.module_from_spec(spec); m.__package__ = "forge_design.meshing"
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
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


# ------------------------------------------------------------------------------------------------ 共通の格子 (run_sern_te_wake_tests と同じ)
k_ = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                              x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
d = k_.design_ramp(M_c=3.3, f=0.45)
THB = float(k_.TH[-1, 0])
P0 = SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=60, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5, L_sw=0.67,
                      interface_angle=THB, top_ext_angle=d.info["theta_e"], cowl_thickness=0.005, ext_top=True, nj_ext_top=9,
                      vehicle_taper=0.3, scale=0.1)
LB = 0.3


def unit_tests():
    cA, hA, BA, iA, ymA = generate_sern_mesh3d(d, P0)
    cB, hB, BB, iB, ymB = generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=LB))
    cL, hL, BL, iL, ymL = generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=LB, te_wake_mode="local"))
    NJ, nz, jm, ni = iA["NJ"], iA["nz"], iA["jm"], iA["ni"]
    NB = ni * NJ * nz
    H = P0.scale
    i0, i1 = iL["te_wake_i_first"], iL["te_wake_i_last"]
    print(f"--- 格子: nodes {iA['nodes']} ni {ni} NJ {NJ} nz {nz}; 変形区間 i {i0}..{i1}")

    # ---------------------------------------------------------------------------------------- (a)
    ref = load_ref_mesher()
    if ref is None:
        check("(a) 変更前のメッシャを git から読めた", False, REF_COMMIT)
    else:
        for label, prm in (("生産に近い構成", P0),
                           ("厚さ 0・外部流ブロックなし", replace(P0, cowl_thickness=0.0, ext_top=False, L_sw=None)),
                           ("外側空間なし (nz_out=0)", replace(P0, nz_out=0, ext_top=False))):
            for sub, kw, want_info in (("キー無し", {}, {}), ("L_b 0 明示", {"te_wake_blend_H": 0.0}, {}),
                                       ("L_b 0 + local", {"te_wake_blend_H": 0.0, "te_wake_mode": "local"}, {}),
                                       ("band L_b 0.3 明示", {"te_wake_blend_H": LB, "te_wake_mode": "band"}, {"te_wake_mode": "band"})):
                p1 = replace(prm, **kw)
                prm_ref = ref.SernMesh3DParams(**{f: getattr(p1, f) for f in ref.SernMesh3DParams.__dataclass_fields__})
                c0, h0, B0, i0r, y0 = ref.generate_sern_mesh3d(d, prm_ref)
                c1, h1, B1, i1r, y1 = generate_sern_mesh3d(d, p1)
                xq = np.linspace(-0.5, float(i0r["x_out"]), 997)
                check(f"(a) {label} / {sub}: 変更前と座標・ヘキサ・境界面・中間線の関数がビット一致",
                      np.array_equal(c0, c1) and np.array_equal(h0, h1) and same_B(B0, B1) and np.array_equal(y0(xq), y1(xq)))
                extra = {kk: v for kk, v in i1r.items() if kk not in i0r}
                check(f"(a) {label} / {sub}: info は変更前と同じ (増えるのは {want_info or 'なし'})",
                      extra == want_info and all(i1r[kk] == v for kk, v in i0r.items()), str(extra))

    # ---------------------------------------------------------------------------------------- (b)
    check("(b) local: 節点数・ヘキサ・境界面の接続が 0 格子と同一", cA.shape == cL.shape and np.array_equal(hA, hL) and same_B(BA, BL))
    check("(b) local: x・z 配列がビット一致", np.array_equal(cA[:, 0], cL[:, 0]) and np.array_equal(cA[:, 2], cL[:, 2]))
    mv = np.flatnonzero(np.any(cA != cL, axis=1))
    ii, jj = mv // (NJ * nz), (mv // nz) % NJ
    check("(b) local: 動くのは主ブロック (base 節点) だけ", mv.size > 0 and np.all(mv < NB), f"{mv.size} 節点")
    check("(b) local: 動く station = 変形区間 (te_wake_i_first..te_wake_i_last)、区間外はビット不変",
          set(np.unique(ii).tolist()) <= set(range(i0, i1 + 1)) and np.array_equal(cA[:NB].reshape(ni, -1)[:i0], cL[:NB].reshape(ni, -1)[:i0])
          and np.array_equal(cA[:NB].reshape(ni, -1)[i1 + 1:], cL[:NB].reshape(ni, -1)[i1 + 1:]) and np.array_equal(cA[NB:], cL[NB:]),
          f"i {sorted(set(ii.tolist()))}")
    # 式の照合 (前処理なしの倍精度の座標、/H)
    GA = cA[:NB].reshape(ni, NJ, nz, 3) / H
    GL = cL[:NB].reshape(ni, NJ, nz, 3) / H
    xs = GA[:, 0, 0, 0]
    dlt = np.array([float(ymB(x)) - float(ymA(x)) for x in xs])               # 現方式の中間線 − 旧中間線 (各 station の x)
    check("(b) δ(x) は変形区間の外で厳密に 0", np.all(dlt[:i0] == 0.0) and np.all(dlt[i1 + 1:] == 0.0) and np.all(dlt[i0:i1 + 1] != 0.0))
    dd = np.abs(GA[..., 1] - GA[:, jm:jm + 1, :, 1])                          # (ni, NJ, nz) 0 格子の旧中間線からの鉛直距離
    want = dlt[:, None, None] * te_wake_local_weight(dd)
    got = GL[..., 1] - GA[..., 1]
    tol = 4.0 * np.finfo(float).eps * (np.abs(GA[..., 1]) + np.abs(want)) + 1e-15
    err = np.abs(got - want)
    check("(b) 全 base 節点で y の差 = δ(x)·φ(d) (丸め内)", bool(np.all(err <= tol)), f"max 誤差 {float(err.max()):.3e}")
    sec = (slice(i0, i1 + 1),)
    full = dd[sec] <= MS.TE_WAKE_LOCAL_D_FULL_H
    zero = dd[sec] >= MS.TE_WAKE_LOCAL_D_ZERO_H
    trans = ~full & ~zero
    dl3 = np.broadcast_to(dlt[sec][:, None, None], full.shape)
    check("(b) d ≤ 0.02 H の節点は δ だけ動く (φ = 1、中間線を含む)",
          full.sum() > 0 and np.all(np.abs(got[sec][full] - dl3[full]) <= tol[sec][full]), f"{int(full.sum())} 節点")
    check("(b) d ≥ 0.30 H の節点はビット不変 (φ = 0)", zero.sum() > 0 and np.array_equal(GA[sec][..., 1][zero], GL[sec][..., 1][zero]),
          f"{int(zero.sum())} 節点")
    phi_obs = got[sec][trans] / dl3[trans]
    t_ = (dd[sec][trans] - 0.02) / 0.28
    phi_lit = 1 - 10 * t_ ** 3 + 15 * t_ ** 4 - 6 * t_ ** 5
    check("(b) 遷移帯 (0.02 < d < 0.30 H) の φ = 1 − 10t³ + 15t⁴ − 6t⁵、t = (d/H − 0.02)/0.28 (0 < φ < 1)",
          trans.sum() > 0 and np.all(np.abs(phi_obs - phi_lit) < 1e-9) and np.all((phi_obs > 0) & (phi_obs < 1)),
          f"{int(trans.sum())} 節点、φ {float(phi_obs.min()):.4f}〜{float(phi_obs.max()):.4f}")
    jidx = np.broadcast_to(np.arange(NJ)[None, :, None], full.shape)
    for lab, side in (("下側の帯 (j < jm)", jidx < jm), ("上側の帯 (j > jm)", jidx > jm)):
        m_ = trans & side
        check(f"(b) {lab}にも遷移帯の節点があり同じ式が成り立つ", m_.sum() > 0 and np.all(np.abs(got[sec][m_] - (dl3 * te_wake_local_weight(dd[sec]))[m_])
                                                                     <= tol[sec][m_]), f"{int(m_.sum())} 節点")
    yb = cB[:NB].reshape(ni, NJ, nz, 3) / H
    check("(b) 中間線 (j = jm) は現方式の Hermite と丸め内で一致",
          np.all(np.abs(GL[:, jm, :, 1] - yb[:, jm, :, 1]) <= 4 * np.finfo(float).eps * np.abs(yb[:, jm, :, 1]) + 1e-15))
    check("(b) 下端 (j = 0)・上線 (j = NJ−1) は動かない", np.array_equal(GA[:, 0], GL[:, 0]) and np.array_equal(GA[:, NJ - 1], GL[:, NJ - 1]))
    check("(b) y_bot は旧中間線の x_out の値のまま", iL["y_bot"] == iA["y_bot"])
    mset = np.zeros(len(cA), bool); mset[mv] = True
    moved_tags = {nm: int(mset[np.unique(np.asarray(q))].sum()) for nm, q in BA.items() if len(q)}
    check("(b) 外部境界・壁の節点は動かない (動くのは z 一定の面 sym・side_far の面内だけ)",
          all(n == 0 for nm, n in moved_tags.items() if nm not in ("sym", "side_far")), {nm: n for nm, n in moved_tags.items() if n})
    mi, xi = closure(hL, BL)
    check("(b) 境界の閉性", mi == 0 and xi == 0, f"missing {mi} extra {xi}")
    _P = cL[hL]
    _J = np.einsum("ij,ij->i", _P[:, 1] - _P[:, 0], np.cross(_P[:, 3] - _P[:, 0], _P[:, 4] - _P[:, 0]))
    check("(b) ヘキサの符号付きヤコビアンが同符号・非退化", np.all(_J > 1e-18) or np.all(_J < -1e-18))
    check("(b) float32 で新たな節点衝突が起きない",
          np.unique(cL, axis=0).shape[0] - np.unique(cL.astype(np.float32), axis=0).shape[0]
          == np.unique(cA, axis=0).shape[0] - np.unique(cA.astype(np.float32), axis=0).shape[0])
    check("(b) 動いた節点の数は現方式より少ない (帯全体を再配分しない)",
          mv.size < int(np.count_nonzero(np.any(cA != cB, axis=1))), f"local {mv.size} / band {int(np.count_nonzero(np.any(cA != cB, axis=1)))}")

    # ---------------------------------------------------------------------------------------- (c)
    w = te_wake_local_weight(np.array([0.0, 0.01, 0.02, 0.16, 0.30, 0.31, 1.0]))
    check("(c) φ の端点・中点: φ(0) = φ(0.01) = φ(0.02) = 1、φ(0.16) = 0.5、φ(0.30) = φ(0.31) = φ(1) = 0 (厳密)",
          w[0] == 1.0 and w[1] == 1.0 and w[2] == 1.0 and abs(w[3] - 0.5) < 1e-15 and w[4] == 0.0 and w[5] == 0.0 and w[6] == 0.0, str(w))
    dg = np.linspace(0.0, 0.4, 4001); wg = te_wake_local_weight(dg)
    h_ = 1e-6
    der = lambda x: (te_wake_local_weight(x + h_) - te_wake_local_weight(x - h_)) / (2 * h_)
    check("(c) φ は単調非増加で、両端 (0.02・0.30 H) で傾きが 0", np.all(np.diff(wg) <= 0.0) and abs(float(der(0.02))) < 1e-5
          and abs(float(der(0.30))) < 1e-5, f"φ'(0.02) {float(der(0.02)):.2e} φ'(0.30) {float(der(0.30)):.2e}")
    check("(c) 定数 0.02 H / 0.30 H と版の文字列", MS.TE_WAKE_LOCAL_D_FULL_H == 0.02 and MS.TE_WAKE_LOCAL_D_ZERO_H == 0.30
          and TE_WAKE_MODE_VERSION == "local-0.02-0.30-quintic-v1")
    check("(c) 幅は格子パラメータに無い (問題 YAML から変えられない)",
          [f for f in SernMesh3DParams.__dataclass_fields__ if f.startswith("te_wake")] == ["te_wake_blend_H", "te_wake_mode"])
    check("(c) 不正な te_wake_mode は失敗", raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=LB, te_wake_mode="Local")))
          and raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_mode="foo"))))
    saved = MS.TE_WAKE_LOCAL_D_ZERO_H
    try:
        MS.TE_WAKE_LOCAL_D_ZERO_H = iL["te_wake_local_clear_top_H"] + 0.01      # 上線が支持領域に入る状況を作る (試験だけ)
        top_in = raises(lambda: generate_sern_mesh3d(d, replace(P0, te_wake_blend_H=LB, te_wake_mode="local")))
    finally:
        MS.TE_WAKE_LOCAL_D_ZERO_H = saved
    check("(c) 上線 (壁・プルーム線) か下端が支持領域に入る station があれば生成を失敗させる (壁を動かさない)", top_in)
    check("(c) info に上線・下端の余裕 (いまの構成では 0.30 H より遠い)",
          iL["te_wake_local_clear_top_H"] >= 0.30 and iL["te_wake_local_clear_bot_H"] >= 0.30,
          f"上線 {iL['te_wake_local_clear_top_H']:.4f} H、下端 {iL['te_wake_local_clear_bot_H']:.4f} H")

    # ---------------------------------------------------------------------------------------- (d)
    check("(d) local の info: te_wake_mode・版・定数・曲線版",
          iL.get("te_wake_mode") == "local" and iL.get("te_wake_mode_version") == TE_WAKE_MODE_VERSION
          and iL.get("te_wake_local_d_full_H") == 0.02 and iL.get("te_wake_local_d_zero_H") == 0.30
          and iL.get("te_wake_curve_version") == TE_WAKE_CURVE_VERSION)
    check("(d) band の info: te_wake_mode band・版なし", iB.get("te_wake_mode") == "band" and "te_wake_mode_version" not in iB)
    check("(d) local と band の中間線の記録 (端点の傾き・区間・最大移動) は同じ",
          all(iL[kk] == iB[kk] for kk in ("te_wake_x_end", "te_wake_slope_te_deg", "te_wake_slope_end_deg", "te_wake_n_stations",
                                         "te_wake_i_first", "te_wake_i_last", "te_wake_max_dy")))
    from forge_design.evaluate import runner_sern as R2
    from forge_design.evaluate import runner_sern3d as R3
    from forge_design.probdef import load_problem
    check("(d) 識別子の文字列: band は従来のまま、local は版を足す",
          R2.te_wake_effective_of(1.0, TE_WAKE_CURVE_VERSION) == f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION}"
          and R2.te_wake_effective_of(1.0, TE_WAKE_CURVE_VERSION, TE_WAKE_MODE_VERSION)
          == f"L_b=1.0;curve={TE_WAKE_CURVE_VERSION};mode={TE_WAKE_MODE_VERSION}"
          and R2.te_wake_effective_of(0.0, None, TE_WAKE_MODE_VERSION) == "off")
    check("(d) te_wake_mode_spec: local + L_b > 0 だけ版、band・キー無し・L_b 0・2D は None、不正は失敗",
          R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 1.0, "te_wake_mode": "local"}}, 3) == TE_WAKE_MODE_VERSION
          and R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 1.0, "te_wake_mode": "band"}}, 3) is None
          and R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 1.0}}, 3) is None
          and R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 0.0, "te_wake_mode": "local"}}, 3) is None
          and R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 1.0, "te_wake_mode": "local"}}, 2) is None
          and raises(lambda: R2.te_wake_mode_spec({"mesh3d": {"te_wake_blend_H": 1.0, "te_wake_mode": "x"}}, 3)))
    conv = {"sha256": "c" * 64}
    pB = load_problem(CASE / "problem_3d_prod_3op_wallres_lswx08_tewake_g4_B10.yaml")
    pL = load_problem(CASE / "problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml")
    reqB, reqL = R2.eval_method_required(pB, 3, conv), R2.eval_method_required(pL, 3, conv)
    check("(d) 評価方式: local は版の項目が入り band と別の id、band の項目は従来どおり",
          reqL["eval_method"].get("te_wake_mode_version") == TE_WAKE_MODE_VERSION and "te_wake_mode_version" not in reqB["eval_method"]
          and reqL[R2.EVAL_METHOD_ID] != reqB[R2.EVAL_METHOD_ID] and reqL[R2.EVAL_METHOD_ID] is not None
          and set(reqB["eval_method"]) == {"version", "dim", "discretization", "mesh_recipe_sha256", "mesher_source_sha256",
                                           "te_wake_blend_H", "te_wake_curve_version", "converter_sha256"})
    check("(d) 評価方式: 格子生成のレシピに te_wake_mode が入る",
          R2.mesh_recipe(R3.sern_mesh3d_params(pL, 0.0, 0.0), pL.mesh).get("te_wake_mode") == "local"
          and R2.mesh_recipe(R3.sern_mesh3d_params(pB, 0.0, 0.0), pB.mesh).get("te_wake_mode") == "band")
    # prepare の記録 (メッシャの info の実効値) とキャンペーンの要求が一致すること (MOC を解かない: 設計由来の角度はレシピで "design")
    for lab, p_, mi in (("local", pL, {"te_wake_blend_H": 1.0, "te_wake_curve_version": TE_WAKE_CURVE_VERSION, "te_wake_mode": "local",
                                      "te_wake_mode_version": TE_WAKE_MODE_VERSION}),
                        ("band", pB, {"te_wake_blend_H": 1.0, "te_wake_curve_version": TE_WAKE_CURVE_VERSION, "te_wake_mode": "band"})):
        rec = R2.eval_method_record(p_, 3, R3.sern_mesh3d_params(p_, -0.7, 0.05), mi, "node", conv)
        req = R2.eval_method_required(p_, 3, conv)
        check(f"(d) {lab}: prepare の記録の id・識別子 = キャンペーンの要求", rec[R2.EVAL_METHOD_ID] == req[R2.EVAL_METHOD_ID]
              and rec[R2.TE_WAKE_EFFECTIVE] == req[R2.TE_WAKE_EFFECTIVE], f"{rec[R2.TE_WAKE_EFFECTIVE]} / {req[R2.TE_WAKE_EFFECTIVE]}")

    # ---------------------------------------------------------------------------------------- (e)
    check("(e) runner: mesh3d.te_wake_mode local → prm.te_wake_mode local、キー無し → band",
          R3.sern_mesh3d_params(pL, 0.0, 0.0).te_wake_mode == "local" and R3.sern_mesh3d_params(pB, 0.0, 0.0).te_wake_mode == "band")

    # ---------------------------------------------------------------------------------------- (f)
    import te_wake_grid_check as TC
    import te_wake_hex_compare as HX
    gA = TC.from_mesher(cA, hA, BA, iA, H, "A 0"); gL = TC.from_mesher(cL, hL, BL, iL, H, "B local"); gB = TC.from_mesher(cB, hB, BB, iB, H, "B band")
    oL = HX.compare(gA, gL)
    vL, rowsL = HX.rows_4c(oL)
    check("(f) hex_compare A = 0・B = local: 支持領域の外の座標の不変が PASS (動いた節点 = 支持領域の中)",
          oL["support"]["status"] == TC.PASS and oL["support"]["n_moved_outside_support"] == 0 and oL["support"]["n_moved"] == mv.size,
          json.dumps(oL["support"], ensure_ascii=False))
    check("(f) hex_compare: 変形領域の外のヘキサは変わらない・#4c の行が 6 つ",
          oL["outside_region_changed_hex"] == 0 and len(rowsL) == 6, str([(r, s) for r, s, _ in rowsL]))
    oB = HX.compare(gA, gB)
    check("(f) hex_compare A = 0・B = band: 支持領域の行は判定不能 (B が local でない)", oB["support"]["status"] == TC.UND
          and HX.rows_4c(oB)[0] != TC.PASS)
    cT = cL.copy()
    jf = int(np.flatnonzero(dd[i0, :, 0] > 0.5)[0])                          # 変形区間の最初の station で旧中間線から 0.5 H より遠い節点
    cT[(i0 * NJ + jf) * nz + 0, 1] += 1e-6
    oT = HX.compare(gA, TC.from_mesher(cT, hL, BL, iL, H, "B local (支持領域の外を 1 節点動かした)"))
    check("(f) 支持領域の外を 1 節点動かした B → FAIL・総合 FAIL", oT["support"]["status"] == TC.FAIL
          and oT["support"]["n_moved_outside_support"] == 1 and HX.rows_4c(oT)[0] == TC.FAIL)
    iX = dict(iL); iX["te_wake_mode_version"] = "local-0.02-0.30-quintic-v0"
    oX = HX.compare(gA, TC.from_mesher(cL, hL, BL, iX, H, "B local (版違い)"))
    check("(f) B の info の版が現行の定数と違えば判定不能", oX["support"]["status"] == TC.UND)

    # ---------------------------------------------------------------------------------------- (g) 変換した run ディレクトリ (CLI)
    if not R2.CONVERTER.exists():
        check("(g) 変換器が無いので CLI の試験を飛ばす (判定不能。合格扱いにしない)", False, str(R2.CONVERTER))
        return
    import os
    import yaml
    src = CASE / "problem_3d_prod_m6on_g1.yaml"
    cached = R2.design_from_problem(load_problem(src))
    saved_d = R2.design_from_problem
    R2.design_from_problem = lambda p_, design=None: cached
    os.environ.setdefault("FORGE_ALLOW_UNVERIFIED_SPECIES", "1")
    runs = {}
    try:
        for lab, lb, mode in (("A", 0.0, None), ("L", 1.0, "local")):
            y = yaml.safe_load(open(src))
            # run_sern_te_wake_tests (j) と同じ粗い格子 (g1 の x・y 分布のまま z と上流のノズル区間を粗くする)
            y["mesh3d"].update({"nz_in": 4, "nz_out": 3, "ni_noz": 30, "nj_vside": 5, "first_wake_frac": 0.004, "first_z_frac": 0.02,
                                "te_wake_blend_H": lb})
            y["mesh"].update({"nj_ext_top": 9, "first_top_frac": 0.004})
            if mode:
                y["mesh3d"]["te_wake_mode"] = mode
            pth = tmp / f"cli_{lab}.yaml"; pth.write_text(yaml.safe_dump(y, allow_unicode=True))
            runs[lab] = (tmp / f"cli_run_{lab}", pth, R3.prepare(pth, tmp / f"cli_run_{lab}"))
    finally:
        R2.design_from_problem = saved_d
    (rA, _, _), (rL, pL_, infoL) = runs["A"], runs["L"]
    mi = infoL["mesh"]
    check("(g) prepare_info の mesh に te_wake_mode・版、評価方式に版・識別子に mode",
          mi.get("te_wake_mode") == "local" and mi.get("te_wake_mode_version") == TE_WAKE_MODE_VERSION
          and infoL["eval_method"]["eval_method"].get("te_wake_mode_version") == TE_WAKE_MODE_VERSION
          and infoL["eval_method"][R2.TE_WAKE_EFFECTIVE].endswith(f";mode={TE_WAKE_MODE_VERSION}"))
    check("(g) prepare の評価方式の id = 問題 YAML からのキャンペーンの要求 (変換器も同じ)",
          infoL["eval_method"][R2.EVAL_METHOD_ID] == R2.eval_method_required(load_problem(pL_), 3)[R2.EVAL_METHOD_ID] is not None)
    r = subprocess.run([sys.executable, str(DIAG / "te_wake_hex_compare.py"), str(rA), str(rL), "--json", str(tmp / "cli_hx.json")],
                       capture_output=True, text=True)
    hx = json.loads((tmp / "cli_hx.json").read_text()) if (tmp / "cli_hx.json").exists() else {}
    check("(g) CLI: 最終 node 格子の run ディレクトリ A (0)・B (local) で TE_WAKE_LOCAL VERDICT: PASS、終了コード 0、JSON に行",
          r.returncode == 0 and "TE_WAKE_LOCAL VERDICT: PASS" in r.stdout and hx.get("te_wake_local_verdict") == "PASS"
          and len(hx.get("te_wake_local_rows", [])) == 6 and hx.get("support", {}).get("n_moved_outside_support") == 0,
          (r.stdout[-600:] + r.stderr[-600:]), only_on_fail=True)


# ================================================================================================ --production
def to_f32_msh(c):
    """メッシャの msh 出力 (.10g) → 変換器の float32 の座標 (を倍精度で持つ)。"""
    s = np.array([float(f"{v:.10g}") for v in np.asarray(c, float).ravel()]).reshape(c.shape)
    return s.astype(np.float32).astype(np.float64)


def production(label):
    """問題 YAML (B10) からメモリ上で 0 格子・現方式・局所化を作り、対応比較と投入条件の行を出す。戻り = 結果の dict。"""
    import time
    from forge_design.evaluate import runner_sern as R2
    from forge_design.evaluate import runner_sern3d as R3
    from forge_design.probdef import load_problem
    import te_wake_grid_check as TC
    import te_wake_hex_compare as HX
    prob = CASE / f"problem_3d_prod_3op_wallres_lswx08_tewake_{label}_B10.yaml"
    probL = CASE / f"problem_3d_prod_3op_wallres_lswx08_tewake_{label}_L10.yaml"
    p = load_problem(prob)
    t0 = time.time()
    seen = {}

    class _Stop(Exception):
        pass

    def _cap(design, prm):
        seen["design"], seen["prm"] = design, prm
        raise _Stop
    saved = R3.generate_sern_mesh3d
    R3.generate_sern_mesh3d = _cap
    try:
        for pp, key in ((prob, "prm"), (probL, "prmL")):
            try:
                R3.prepare(pp, tmp / f"cap_{label}_{key}", op="m10_on")      # runner と同じ手順で格子パラメータを組む (格子は作らない)
            except _Stop:
                pass
            seen[key + "_got"] = seen["prm"]
    finally:
        R3.generate_sern_mesh3d = saved
    design, prm, prmL = seen["design"], seen["prm_got"], seen["prmL_got"]
    H = float(p.spec["H_m"])
    check(f"(P) {label}: B10 は band・L_b 1.0、L10 は local で他の格子パラメータは同じ",
          prm.te_wake_mode == "band" and prm.te_wake_blend_H == 1.0 and prmL == replace(prm, te_wake_mode="local"))
    print(f"--- {label}: 設計 {time.time() - t0:.1f} s", flush=True)
    res = {"label": label, "H_m": H}
    c0, h0, B0, i0, _ = generate_sern_mesh3d(design, replace(prm, te_wake_blend_H=0.0))
    print(f"--- {label}: 0 格子 {i0['nodes']} 節点・{len(h0)} ヘキサ ({time.time() - t0:.1f} s)", flush=True)
    c0 = to_f32_msh(c0)
    cb, hb, Bb, ib, _ = generate_sern_mesh3d(design, prm)
    same_b = np.array_equal(hb, h0) and same_B(Bb, B0); del hb, Bb
    cb = to_f32_msh(cb)
    cl, hl, Bl, il, _ = generate_sern_mesh3d(design, prmL)
    same_l = np.array_equal(hl, h0) and same_B(Bl, B0); del hl, Bl
    cl = to_f32_msh(cl)
    check(f"(P) {label}: 3 格子の接続・境界面が同一", same_b and same_l)
    print(f"--- {label}: 3 格子 ({time.time() - t0:.1f} s); local の上線・下端の余裕 {il['te_wake_local_clear_top_H']:.4f} / "
          f"{il['te_wake_local_clear_bot_H']:.4f} H", flush=True)
    res["clear_top_H"], res["clear_bot_H"] = il["te_wake_local_clear_top_H"], il["te_wake_local_clear_bot_H"]
    gA = TC.from_mesher(c0, h0, B0, i0, H, f"{label} 0 格子 (float32)")
    gb = TC.from_mesher(cb, h0, B0, ib, H, f"{label} 現方式 (float32)")
    gl = TC.from_mesher(cl, h0, B0, il, H, f"{label} 局所化 (float32)")
    for g in (gA, gb, gl):
        g["coord_repr"] = "float32"
    keys = ("n_lt_thr_A", "n_lt_thr_B", "new_lt_thr", "worsened_existing_lt_thr", "min_A", "min_B", "region_min_A", "region_min_B",
            "n_region_hex", "n_moved_nodes", "new_ar_gt_max", "ar_gt5000_A", "ar_gt5000_B", "region_ar_max_B", "nonpos_jac_hex_B",
            "outside_region_changed_hex")
    for tag, g in (("band", gb), ("local", gl)):
        o = HX.compare(gA, g)
        v4, rows = HX.rows_4c(o)
        res[tag] = {kk: o[kk] for kk in keys}
        res[tag]["support"] = o["support"]; res[tag]["rows_4c"] = [(r, s, val) for r, s, val in rows]; res[tag]["verdict_4c"] = v4
        print(f"--- {label} {tag} ({time.time() - t0:.1f} s): " + ", ".join(f"{kk} {o[kk]}" for kk in keys), flush=True)
        print(f"    支持領域: {json.dumps(o['support'], ensure_ascii=False)}")
        for r, s, val in rows:
            print(f"    {s:6s} | {r} | {val}")
        del o
    # 投入条件の行 (後縁・復帰の折れ、層厚、固定座標と境界、primal の品質)。品質・双対のファイルは無いので判定不能
    for tag, g in (("band", gb), ("local", gl)):
        v, R = TC.admission(gA, g)
        res[tag]["admission_rows"] = [{kk: r[kk] for kk in ("row", "item", "status", "value", "ref_A", "num")} for r in R.rows]
        res[tag]["return_fold_deg"] = next(r["num"] for r in R.rows if r["row"] == "復帰区間の折れ" and r["num"] is not None)
        print(f"--- {label} {tag} admission ({time.time() - t0:.1f} s): " + ", ".join(f"{row} {R.status_of(row)}"
                                                                              for row in dict.fromkeys(r["row"] for r in R.rows)))
        for r in R.rows:
            if r["row"] not in ("前提", "品質・双対の判定"):
                print(f"    {r['status']:6s} | {r['row']} | {r['item']} | {r['value']} | A {r['ref_A'] or '-'}")
    # 遷移帯の層厚比 (局所化 / 0 格子、変形区間の station の主ブロック) と、第 1〜第 10 層の旧中間線からの最大距離
    S = TC.Struct(il)
    G0, GL = S.grid(c0), S.grid(cl)
    sl = slice(il["te_wake_i_first"], il["te_wake_i_last"] + 1)
    dy0, dyl = np.diff(G0[sl, :, :, 1], axis=1), np.diff(GL[sl, :, :, 1], axis=1)
    ratio = dyl / dy0
    res["layer_ratio_min"], res["layer_ratio_max"] = float(ratio.min()), float(ratio.max())
    jm = S.jm
    # 第 1〜第 10 層の距離は投入条件の層厚の列 (i = i_te .. i_last + 1、te_wake_grid_check の「層の Δy」と同じ) で測る (codex の値の定義)
    y0 = G0[S.i_te:il["te_wake_i_last"] + 2][..., 1]                          # (列, NJ, nz)
    res["layers_1_10_max_dist_H"] = float(np.max(np.abs(y0[:, [jm - 10, jm + 10], :] - y0[:, jm:jm + 1, :])) / H)
    print(f"--- {label}: 遷移帯の層厚比 (局所化/0 格子) {res['layer_ratio_min']:.4f}〜{res['layer_ratio_max']:.4f}、"
          f"第 10 層の旧中間線からの最大距離 (列 i_te..i_last+1) {res['layers_1_10_max_dist_H']:.6f} H ({time.time() - t0:.1f} s)", flush=True)
    return res


CODEX_G4 = {"band": {"new_lt_thr": 4180, "worsened_existing_lt_thr": 6105, "n_lt_thr_B": 44880, "return_fold_deg": 5.999851},
            "local": {"new_lt_thr": 0, "worsened_existing_lt_thr": 0, "n_lt_thr_B": 40700, "region_min_B": 0.651909,
                      "return_fold_deg": 5.170986, "new_ar_gt_max": 0, "nonpos_jac_hex_B": 0},
            "min_B": 0.452350, "layer_ratio_min": 0.2495, "layers_1_10_max_dist_H": 0.002963}


def check_g4(res):
    for tag in ("band", "local"):
        for kk, v in CODEX_G4[tag].items():
            got = res[tag][kk]
            ok = (got == v) if isinstance(v, int) else abs(got - v) < 5e-7 * max(1.0, abs(v)) + 5e-7
            check(f"(P) g4 {tag}: {kk} = codex の予備値 {v}", ok, f"{got}")
        check(f"(P) g4 {tag}: 全域最小 = codex の予備値 0.452350", abs(res[tag]["min_B"] - CODEX_G4["min_B"]) < 5e-7, f"{res[tag]['min_B']}")
    check("(P) g4 local: 遷移帯の層厚比の最小 = codex の 0.2495", abs(res["layer_ratio_min"] - 0.2495) < 5e-5, f"{res['layer_ratio_min']:.6f}")
    check("(P) g4: 第 1〜第 10 層の旧中間線からの最大距離 = codex の 0.002963 H", abs(res["layers_1_10_max_dist_H"] - 0.002963) < 5e-7,
          f"{res['layers_1_10_max_dist_H']:.7f}")
    check("(P) g4 local: 支持領域の外の座標はビット不変", res["local"]["support"]["status"] == "PASS"
          and res["local"]["support"]["n_moved_outside_support"] == 0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--production", nargs="*", choices=("g3", "g4"), default=None)
    ap.add_argument("--json", help="--production の結果の書き出し先")
    a = ap.parse_args()
    unit_tests()
    out = {}
    for lab in (a.production or []):
        out[lab] = production(lab)
        if lab == "g4":
            check_g4(out[lab])
    if a.json and out:
        json.dump(out, open(a.json, "w"), indent=1, ensure_ascii=False, default=str)
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
    sys.exit(1 if FAIL else 0)
