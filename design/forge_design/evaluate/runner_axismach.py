r"""axis-Mach チェーン A4: Hall + 5次 Hermite 軸 Mach → 逆 MOC → node Euler 評価。

plan: plans/accepted/tooling-nozzle-axismach-chain.md §6 A4。既存経路 (runner_wt の
モード F / runner_walldriven) には触れない。

チェーン: HallThroat (初期値線 + 軸アンカー) → QuinticHermiteAxisLaw
(自由度 L_c のみ) → inverse_design (E→F 閉包込み) → wall_qa →
AxisMachCFDWall (直管 + U→T Hermite + 設計壁 spline。旧・縦線構成のみ
[T, x0] に骨接放物線が入る) →
TFI メッシュ → node Euler (段階起動)。

CFD-in-the-loop アンカー更新 (A5) は problem YAML の geometry キーで受ける:
  x_reach_cfd:     壁始点発 C⁻ の軸着地点 (node 基準 run から)
  x_reach_anchor:  [M, M', M''] (同 run の軸平滑化フィットから)
  axis_segment_run: [x0, x_reach) の実測軸 M を target_moc に渡す元 run
これらが無ければ初回 (Hall アンカー、x_A = x0)。
初期値線は geometry.start_line で選ぶ ('throat_char' = スロート特性線 [A8] /
'vertical' = M_start の縦線 [旧構成])。壁の決め方は geometry.wall_mode
('flux' = 断面の質量流束閉包 [A9] / 'streamline' = 流線積分 [旧構成])。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
      case/41.wind_tunnel_design/problem_m4_axismach.yaml run_dir [--prepare-only]
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.axis_law import QuinticHermiteAxisLaw, KnotQuinticAxisLaw
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import HallThroat
from ..geometry.wall_axismach import (AxisMachCFDWall, area_ratio_isentropic,
                                      wall_qa)
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge
from .runner_wt import (_bcond, _config_euler, _config_euler_node,
                        _config_sst_node)


def axis_curve_node(run_dir, scale: float, lam: float = 1e-5, mode: str = "evenfit"):
    r"""node run の軸列から平滑化スプライン M(x) を作る (§15 の実装)。

    mode="evenfit" (既定): 軸ノードを除外し r>0 の 4 点で M=a₀+a₂r² を最小二乗した a₀ を軸値に
    使う (`axis_extract.axis_curve_evenfit`)。旧生産設定 (`nodeAxisDirichlet: 1`, 撤去済) の軸ノードは第一内点の
    コピーで ½M_rr r₁² のバイアス (case/41 で ~0.08% M_d) を持っていたため。現行 (軸 DOF) でも外挿の方が精度が良い。
    mode="node": 従来の軸ノード直読。末尾 3 スナップ平均 → `make_smoothing_spline`。
    M, M', M'' は**同一スプライン**から評価する (生の 2 階差分禁止)。
    戻り値: (spline, x_min, x_max) — x は r_t 単位。"""
    if mode == "evenfit":
        from .axis_extract import axis_curve_evenfit
        return axis_curve_evenfit(run_dir, scale, lam=lam)
    import h5py
    from scipy.interpolate import make_smoothing_spline
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) > 0]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    with h5py.File(rd / "nozzle.h5") as nz:
        nc = nz["/MESH/COORD"][:].reshape(-1, 3)
    Ms = []
    for rf in res[-3:]:
        with h5py.File(rf) as f:
            Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
        if len(Ux) != len(nc):
            raise ValueError("node run でない (VALUE 長 != 節点数)")
        Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    M = np.mean(Ms, axis=0)
    ax = nc[:, 1] < 1e-12 * max(float(nc[:, 1].max()), 1.0)
    x = nc[ax, 0] / scale
    o = np.argsort(x)
    x, Ma = x[o], M[ax][o]
    spl = make_smoothing_spline(x, Ma, lam=lam)
    return spl, float(x.min()), float(x.max())


def _gam_or_gas(p: Problem):
    """MOC/幾何関数に渡す γ: cpg なら float、semiperfect ならガスモデル
    (cfd_gas: cpg のときは設計も CPG に落とす — 熱力学の一致)。"""
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return p.gamma
    return p.gas_model if p.is_semiperfect else p.gamma


def _apply_gas_to_config(cfg: str, p: Problem, run_dir) -> str:
    """semi-perfect のとき forge config を **単一擬似種 TP** (thermalMethod 2) に書き換え、
    NASA-9 混合擬似種を run_dir/species_db.yaml へ出す。cpg なら無変更。
    設計 (MOC) と同一係数の熱力学で CFD が回る (forge 内蔵 DB と同じ CEA 値)。"""
    # 切り分け用: evaluate.axisym_method で CPG でも SU2 流軸対称に切替可
    if int(p.evaluate.get("axisym_method", 0)) == 1:
        cfg = cfg.replace("isAxisymmetric: 1", "isAxisymmetric: 1, axisymMethod: 1", 1)
    if not p.is_semiperfect or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        # cfd_gas: cpg = 設計は semi-perfect のまま CFD だけ CPG(γ*, cp 参照値) で回す
        # (TP × node 軸対称の forge 側発散 [case/42 run_0001] の回避。相対比較には十分)
        return cfg
    import yaml as _yaml
    from ..gas.semiperfect import mixture_pseudo_species, mixture_pseudo_species_split
    Y = dict(p.raw["gas"]["species"])
    species_list = _tp_species_list(p)
    if species_list == ["MIX"]:
        db = mixture_pseudo_species(Y, "MIX")
    else:
        # split_h2o (plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md): H₂O 以外を擬似種
        # MIXDRY、H₂O を独立種にした 2 種 TP。凝縮 (condGasSpecies=1) が指せる。
        db, _, _ = mixture_pseudo_species_split(Y, keep=(str(p.evaluate.get("tp_keep_species", "H2O")),))
    (Path(run_dir) / "species_db.yaml").write_text(_yaml.safe_dump(db, sort_keys=False))
    # thermalMethod 0 → 2、species/speciesDBFile を physProp に追加 (cp/gamma は参照値のまま
    # 残すが TP では NASA-9 が優先される)
    cfg = cfg.replace("thermalMethod: 0", "thermalMethod: 2", 1)
    # [2026-08-16] nodeAxisDirichlet は撤去済み (node は軸ノードを DOF として解く整合セットが常時 ON、
    # TP スカラーの軸ピン問題は消滅)。axisRFloor は evaluate.axis_r_floor 指定時のみ従来どおり付ける。
    floor = float(p.evaluate.get("axis_r_floor", 0.0))
    if floor > 0:
        cfg = cfg.replace("isAxisymmetric: 1", f"isAxisymmetric: 1, axisRFloor: {floor}", 1)
    # thermoHrefTemp (sensible-enthalpy datum): 絶対基準 (生成エンタルピー込み) のままだと
    # 陰解法 Jacobian の χ_eos = c² − κh が桁違いになり block-DPLUR が軸近傍で発散する
    # (case/42 run_0020–0025 で切り分け: 一定 cp 種/陽解法は完走、実 NASA-9 + 陰解法だけ発散、
    #  thermoHrefTemp 298.15 で完走)。IC の roe も同じ datum で作る (paste_isentropic_ic の h_ref)。
    href = float(p.evaluate.get("thermo_href_temp", 298.15))
    sp_txt = "[" + ", ".join(species_list) + "]"
    cfg = cfg.replace("cp: %s, gamma: %s}" % (p.cp, p.gamma),
                      "cp: %s, gamma: %s,\n           species: %s, speciesDBFile: \"species_db.yaml\", thermoHrefTemp: %s}"
                      % (p.cp, p.gamma, sp_txt, href), 1)
    if f"species: {sp_txt}" not in cfg:
        raise RuntimeError("_apply_gas_to_config: physProp の書き換えに失敗 (テンプレート変更?)")
    # 凝縮 (evaluate.condensation: dict) — forge の condensation ブロックをそのまま通す
    cond = p.evaluate.get("condensation")
    if cond:
        cond = dict(cond)
        if len(species_list) >= 2 and "condGasSpecies" not in cond:
            cond["condGasSpecies"] = species_list.index(str(p.evaluate.get("tp_keep_species", "H2O")).upper())
        cfg = cfg.rstrip("\n") + "\ncondensation: {" + ", ".join(f"{k}: {v}" for k, v in cond.items()) + "}\n"
    return cfg


def _tp_species_list(p: Problem) -> list:
    """TP の species リスト。`evaluate.tp_species`: 'pseudo' (既定, [MIX]) / 'split_h2o'
    ([MIXDRY, H2O]; `tp_keep_species` で残す種を変えられる)。"""
    mode = str(p.evaluate.get("tp_species", "pseudo"))
    if mode == "pseudo":
        return ["MIX"]
    if mode == "split_h2o":
        return ["MIXDRY", str(p.evaluate.get("tp_keep_species", "H2O")).upper()]
    raise ValueError("evaluate.tp_species は 'pseudo' か 'split_h2o'")


def _tp_species_Y(p: Problem):
    """IC/BC 用の組成 [Y_s] (species リスト順)。pseudo なら None (roY 不要)。"""
    if _tp_species_list(p) == ["MIX"] or not p.is_semiperfect \
            or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return None
    from ..gas.semiperfect import mixture_pseudo_species_split
    _, Ys, order = mixture_pseudo_species_split(dict(p.raw["gas"]["species"]),
                                                keep=(str(p.evaluate.get("tp_keep_species", "H2O")),))
    return [Ys[k] for k in order]


def _bcond_with_species(txt: str, Ys) -> str:
    """bcond の inlet floats に Y0.. を追記 (多成分 TP の入口組成)。"""
    if Ys is None:
        return txt
    ins = ", ".join(f"Y{i}: {y:.8f}" for i, y in enumerate(Ys))
    return txt.replace("floats: {Pt:", "floats: {" + ins + ", Pt:", 1)



def design_chain(p: Problem) -> dict:
    """Hall (+CFD アンカー) → Hermite law → 逆 MOC → 壁 QA → CFD 壁。決定的。

    ガス: `p.gas_model` (cpg なら float γ と等価 / semiperfect なら NASA-9 テーブル)。
    Hall 遷音速級数は定数 γ 前提なので**スロートの局所 γ* を渡す。MOC の ν↔M・
    質量流束・面積比はガスモデル経由 (`moc_kernel._is_gas` 規約)。"""
    gas = p.gas_model
    # cfd_gas: cpg のときは**設計も** CPG(γ 参照値) で作る。設計だけ semi-perfect にすると
    # 壁 (A/A*=13.1) と CFD (CPG γ=1.309 なら A/A*=15.3) の熱力学が食い違い、出口 M が
    # 3.86 で止まる (case/42 run_0003–0011 で実測、破棄)。設計と CFD の熱力学は必ず一致。
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        from ..gas import GasCPG
        gas = GasCPG(p.gamma, p.cp)
    use_gas = getattr(gas, "kind", "cpg") == "semiperfect"
    g = gas if use_gas else p.gamma                   # MOC/面積比に渡す「γ or ガス」
    g_hall = float(gas.gamma_throat(float(p.spec["Tt"])))
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))
    n_start = int(p.geometry.get("n_start", 41))
    ht = HallThroat(R=R, gamma=g_hall)
    # 初期値線: 'throat_char' = スロート壁点発 C⁻ (CONTUR 流、壁が T から MOC 出力に
    # なる) / 'vertical' = M_start の縦線 (旧構成、回帰対照)。x0 = その軸着地点。
    start_line = str(p.geometry.get("start_line", "vertical"))
    if start_line == "throat_char":
        x0 = float(ht.throat_characteristic(n=n_start)[0][0])
    elif start_line == "vertical":
        x0 = ht.x_axis_of_mach(M_start)
    else:
        raise ValueError("geometry.start_line は 'throat_char' か 'vertical'")

    # --- アンカー: 初回 = Hall 解析値 (x_A = x0)。反復 = CFD 実測 (x_A = x_reach) ---
    x_reach_cfd = p.geometry.get("x_reach_cfd")
    seg = None
    if x_reach_cfd is not None:
        x_A = float(x_reach_cfd)
        anc = p.geometry.get("x_reach_anchor")
        if anc is None or len(anc) != 3:
            raise ValueError("x_reach_cfd 指定時は x_reach_anchor: [M,M',M''] が必須")
        M_A, Mp_A, Mpp_A = float(anc[0]), float(anc[1]), float(anc[2])
        seg_run = p.geometry.get("axis_segment_run")
        if seg_run:
            seg, _, _ = axis_curve_node(seg_run, float(p.spec["r_throat"]))
    else:
        x_A = x0
        M_A, Mp_A, Mpp_A = ht.axis_anchor(x0)

    # --- 軸 M 則: 'quintic' = 単一 5 次 Hermite (DOF L_c) / 'knot' = 内部 knot 1 個の
    # 区分 C² (A6, DOF L_c + M_knot)。高マッハ (M6) では単一 quintic の L_c 上限
    # (単調性) が短すぎて壁角 θ_w > μ_w の fold を起こすので knot を使う。
    # 'bspline_M' (B) / 'bspline_dnu' (C): axislaw-smoothness 比較 (plans/accepted/
    # tooling-nozzle-axislaw-smoothness.md)。$L_c$ に上限はなく常に Lc_mode: explicit。
    axis_law = str(p.geometry.get("axis_law", "quintic"))
    M_K = None
    if axis_law == "quintic":
        lo, hi = QuinticHermiteAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md)
    elif axis_law == "knot":
        # 既定 knot Mach = 急膨張の終わり。M_A + 0.25ΔM (M4.2 で 1.9、M6 で 2.4)
        M_K = float(p.geometry.get("M_knot", M_A + 0.25 * (Md - M_A)))
        lo, hi = KnotQuinticAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md, M_K)
    elif axis_law in ("bspline_M", "bspline_dnu", "onepoint"):
        lo, hi = 1e-2, float("inf")
    else:
        raise ValueError("geometry.axis_law は 'quintic'/'knot'/'bspline_M'/'bspline_dnu'/'onepoint'")
    # --- L_c の決め方 (Lc_mode): 許容窓は常に計算する ---
    #   'explicit'    = dv.L_c を直接与える (窓内検査)
    #   'max'         = 窓上限×0.98
    #   'from_length' = dv.L_total (スロート x=0 → 物理出口 x_F の長さ, r_t 単位) を
    #                   設計変数とし、軸 M 極大点 x_E (= x_A + L_c) は終端特性線込みの
    #                   設計パスから逆算する (E→F 一様化区間は物理で決まるため)。
    #                   plans/accepted/tooling-nozzle-axismach-length-dv.md

    def _make_law(L_c: float):
        if axis_law == "knot":
            return KnotQuinticAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md, M_K)
        if axis_law == "bspline_M":
            from ..geometry.axis_law_bspline import MonotoneBSplineAxisLaw
            return MonotoneBSplineAxisLaw(
                x_A, x_A + L_c, M_A, Mp_A, Mpp_A, Md,
                n_interior=int(p.geometry.get("bspline_n_interior", 15)),
                exit_curvature=str(p.geometry.get("bspline_exit_curvature", "hard")),
                spread_x_frac=float(p.geometry.get("bspline_spread_x_frac", 0.5)),
                spread_M_frac=p.geometry.get("bspline_spread_M_frac", 0.75))
        if axis_law == "bspline_dnu":
            from ..geometry.axis_law_bspline import NonnegDnuBSplineAxisLaw
            return NonnegDnuBSplineAxisLaw(
                x_A, x_A + L_c, M_A, Mp_A, Mpp_A, Md, gas,
                n_interior=int(p.geometry.get("bspline_n_interior", 20)),
                exit_slope=str(p.geometry.get("bspline_exit_curvature", "hard")),
                spread_x_frac=float(p.geometry.get("bspline_spread_x_frac", 0.5)),
                spread_M_frac=p.geometry.get("bspline_spread_M_frac", 0.75))
        if axis_law == "onepoint":
            # D 案 (plans/accepted/tooling-nozzle-axislaw-onepoint.md): 端点アンカー固定 +
            # 内部補間点 1 点 (ξ_P, η_P) の C⁴ 区分 5 次。実験用選択肢。
            from ..geometry.axis_law_onepoint import OnePointC4AxisLaw
            return OnePointC4AxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md,
                                     float(p.geometry["onepoint_xi_P"]),
                                     float(p.geometry["onepoint_eta_P"]))
        return QuinticHermiteAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md)

    def _checked_law(L_c: float):
        law = _make_law(L_c)
        v = law.gates()["violations"]
        if v:
            raise ValueError("軸 Mach law ゲート不合格: " + "; ".join(v))
        return law

    rF_pred = float(np.sqrt(area_ratio_isentropic(Md, g)))

    def _run_inverse(law):
        # target: [x0, x_A) は実測 (反復時) / [x_A, x_E] law / 以降 M_d
        def target_moc(x: float) -> float:
            x = float(x)
            if x < x_A:
                if seg is not None:
                    return float(seg(np.float64(x)))
                return float(ht.mach(x, 0.0)) if x_reach_cfd is None else float(law(x_A))
            return float(law(x))
        x_end = law.x_E + float(p.geometry.get("x_end_margin", 2.3)) \
            * rF_pred * float(np.sqrt(Md * Md - 1.0))
        return inverse_design(ht, target_moc, x_axis_end=float(x_end),
                              n_axis=int(p.geometry.get("n_axis_inv", 500)),
                              n_start=n_start, gamma=g,
                              dx_wall=float(p.geometry.get("dx_wall", 0.02)),
                              th_wall0=float(np.arctan(x0 / R)), M_start=M_start,
                              exit_mode=str(p.geometry.get("exit_mode", "characteristic")),
                              x_E=law.x_E, M_d=Md, start_line=start_line,
                              wall_mode=str(p.geometry.get("wall_mode", "streamline")),
                              blend_width=float(p.geometry.get("wall_blend_width", 1.0)),
                              axis_dx0=p.geometry.get("axis_dx0"))

    mode = str(p.geometry.get("Lc_mode", "explicit"))
    solve_diag = None
    if mode == "max":
        if axis_law in ("knot", "bspline_M", "bspline_dnu", "onepoint") and hi >= 400.0:
            raise ValueError(f"{axis_law} 則の L_c に上限はない (窓が開放) — "
                             "Lc_mode: explicit で L_c を与えること")
        L_c = hi * 0.98
        law = _checked_law(L_c)
        res = _run_inverse(law)
    elif mode == "explicit":
        L_c = float(dv_value(p, "L_c"))
        if not (lo <= L_c <= hi):
            raise ValueError(f"L_c = {L_c:.4g} が許容窓 ({lo:.4g}, {hi:.4g}) 外 "
                             "(M'≥0 単調ゲート)")
        law = _checked_law(L_c)
        res = _run_inverse(law)
    elif mode == "from_length":
        # 逆問題: x_F(L_c) = L_total を L_c について解く。x_F − x_E ≈ r_F√(Md²−1)
        # (終端 Mach line 近似, 実測 0.04% 一致) なので写像はほぼ**傾き 1** の単調。
        # ただし離散 MOC の x_F には解像度依存のノイズ床がある (n_axis 500 で
        # ~0.05 r_t — 隣接 L_c 間で階段状。secant の局所勾配推定は壊れる) ため、
        # ステップ = −残差 の固定点反復 (縮小率 |1 − dx_F/dL_c| ≈ 0.05) を使い、
        # 最良反復点を採用する。既定 tol 0.05 はこのノイズ床相当 — より詰めるなら
        # n_axis_inv を上げて L_total_tol を下げる。
        L_target = float(dv_value(p, "L_total"))
        tol = float(p.geometry.get("L_total_tol", 0.05))
        n_iter_max = int(p.geometry.get("L_total_maxiter", 8))
        # 窓端は数値マージンを取ってクランプ (単調ゲート境界での law 破綻を避ける)
        lo_c = lo * 1.001
        hi_c = hi * 0.999 if np.isfinite(hi) else hi

        def _clamp(v: float) -> float:
            return float(min(max(v, lo_c), hi_c))

        n_eval = 0

        def _eval(L_c: float):
            nonlocal n_eval
            law = _checked_law(L_c)
            res = _run_inverse(law)
            xF = float(res["exit"].get("x_F", float("nan")))
            if not np.isfinite(xF):
                raise ValueError("Lc_mode: from_length — 終端特性線の x_F 追跡に失敗 "
                                 "(場が短い: geometry.x_end_margin を増やす)")
            n_eval += 1
            return xF - L_target, law, res

        L_c = _clamp(L_target - x_A - rF_pred * float(np.sqrt(Md * Md - 1.0)))
        f0, law, res = _eval(L_c)
        best = (abs(f0), L_c, f0, law, res)
        for _ in range(n_iter_max):
            if abs(f0) <= tol:
                break
            L_new = _clamp(L_c - f0)
            if L_new == L_c:               # 窓端に張り付いた — 窓内で到達不能
                raise ValueError(
                    f"dv.L_total = {L_target:.4g} は L_c 許容窓 ({lo:.4g}, {hi:.4g}) "
                    f"内で実現不能 (窓端 L_c = {L_c:.4g} で x_F = {L_target + f0:.4g})")
            f0, law, res = _eval(L_new)
            L_c = L_new
            if abs(f0) < best[0]:
                best = (abs(f0), L_c, f0, law, res)
            else:                          # 改善停止 = ノイズ床に到達
                break
        _, L_c, f0, law, res = best
        if abs(f0) > tol:
            raise ValueError(
                f"Lc_mode: from_length 未収束 (|x_F − L_total| = {abs(f0):.3g} > "
                f"tol {tol:.3g}) — 逆 MOC の離散化ノイズ床の可能性: "
                "geometry.n_axis_inv を上げるか L_total_tol を緩める")
        solve_diag = {"L_total_target": float(L_target),
                      "xF_residual": float(f0), "n_design_evals": int(n_eval),
                      "tol": float(tol)}
    else:
        raise ValueError("geometry.Lc_mode は 'explicit' / 'max' / 'from_length'")
    gates = law.gates()
    qa = wall_qa(res["wall"], Md, law.x_E, g, R=R)
    if qa["violations"]:
        raise ValueError("壁 QA 不合格: " + "; ".join(qa["violations"]))
    # 壁表現 (A14): 'interp' = 補間 5 次 B-spline (現行・比較基準) /
    # 'lsq' = 制約付き最小二乗 B-spline (n_cp は誤差ゲートで自動、または wall_ncp)
    wall_repr = str(p.geometry.get("wall_repr", "interp"))
    wkw = dict(R=R, r_U=float(p.geometry.get("r_inlet", 2.5)),
               L_U=float(p.geometry.get("L_U", 3.5)),
               L_pipe=float(p.geometry.get("L_pipe", 0.5)))
    if wall_repr == "lsq":
        from ..geometry.wall_axismach import LSQBsplineCFDWall
        wall = LSQBsplineCFDWall(res["wall"], n_cp=p.geometry.get("wall_ncp"),
                                 tol_dr=float(p.geometry.get("wall_lsq_tol_dr", 5e-4)),
                                 tol_dtheta_deg=float(p.geometry.get("wall_lsq_tol_dtheta", 0.05)),
                                 **wkw)
    elif wall_repr == "interp":
        wall = AxisMachCFDWall(res["wall"], **wkw)
    else:
        raise ValueError("geometry.wall_repr は 'interp' か 'lsq'")
    msgs = wall.validate()
    if msgs:
        raise ValueError("axis-Mach 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": res["wall"], "law": law, "qa": qa,
            "exit": {k: v for k, v in res["exit"].items() if k != "term_path"},
            "x0": float(x0), "x_A": float(x_A), "x_E": float(law.x_E),
            "L_c": float(L_c), "Lc_window": (float(lo), float(hi)),
            "Lc_mode": mode, "Lc_solve": solve_diag,
            "axis_law": axis_law, "M_knot": M_K,
            "x_K": (float(law.x_K) if axis_law == "knot" else None),
            "anchor": (float(M_A), float(Mp_A), float(Mpp_A)),
            "anchor_source": ("cfd" if x_reach_cfd is not None else "hall"),
            "start_line": start_line, "wall_mode": res["wall_mode"],
            "wall_repr": wall_repr,
            "gas": (gas.summary() if hasattr(gas, "summary")
                    else {"kind": "cpg", "gamma": p.gamma, "cp": p.cp}),
            "gamma_hall": g_hall,
            "wall_fit": getattr(wall, "fit_diag", None),
            "Md": Md, "R": R, "gates": gates,
            # cplus 閉包では構成的に 1 になる循環指標なので出さない (A9 の教訓)
            "mdot_ratio_moc": (None if not np.isfinite(res["mdot_exit"])
                               else float(res["mdot_exit"] / res["mdot_start"])),
            "cd_series": float(ht.cd_series())}


def prepare(problem_path, run_dir, nsteps=None, ic_from=None, cfl_main=None, implicit_relax=None) -> dict:
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 記録: 目標軸分布 (x0 → x_E) と設計壁
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
    (run_dir / "bcondConfig.yaml").write_text(_bcond_with_species(_bcond(p, euler=True), _tp_species_Y(p)))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質検査は cell 変換の一時コピー (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node 変換 (node config を書いてから変換 — 必須) + 等エントロピー IC
    cfg_e = _apply_gas_to_config(_config_euler_node(p, n, out_int, 4.0, 1), p, run_dir)
    if cfl_main is not None:
        cfg_e = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {float(cfl_main)}, cfl_pseudo: {float(cfl_main)}", cfg_e)
    if implicit_relax is not None:
        cfg_e = cfg_e.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    (run_dir / "solverConfig.yaml").write_text(cfg_e)
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model),
                        h_ref_T=float(p.evaluate.get('thermo_href_temp', 298.15)),
                        species_Y=_tp_species_Y(p))
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
    info = {"chain": "axismach", "discretization": "node",
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "Lc_window": list(d["Lc_window"]), "anchor": list(d["anchor"]),
            "Lc_mode": d["Lc_mode"], "Lc_solve": d["Lc_solve"],
            "anchor_source": d["anchor_source"], "start_line": d["start_line"],
            "wall_mode": d["wall_mode"], "wall_repr": d["wall_repr"],
            "axis_law": d["axis_law"], "M_knot": d["M_knot"], "x_K": d["x_K"],
            "L_U": float(p.geometry.get("L_U", 3.5)),
            "gas": d["gas"], "gamma_hall": d["gamma_hall"],
            "wall_fit": d["wall_fit"],
            "Md": d["Md"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "exit": d["exit"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "cd_series": d["cd_series"],
            "nStepOuter": n, "scale_m": scale, "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(run_dir, cfl_main: float | None = None, mid_stage: bool = False, stages: str = "full",
               mesh_h5: str = "nozzle.h5") -> int:
    """soft 段 (1次+cfl0.5, 3000 step) → [mid 段 (2次+cfl1, 3000 step)] → 本段。
    walldriven W3 と同方式・段階起動必須。`cfl_main` で本段 CFL を上書き
    (semi-perfect TP は cfl4 で本段 step ~60 に爆発 — case/42 実測。TP は cfl≤2 が実績
    [[cutler-cpg-vs-tp-dplur-sst]])。`mid_stage=True` で 2 次化と CFL 上げを分離する。"""
    import re
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    if cfl_main is not None:
        cfg_main = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+",
                          f"cfl: {cfl_main}, cfl_pseudo: {cfl_main}", cfg_main)

    def _stage(cfg, nsteps, label):
        cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", cfg)
        cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {nsteps}", cfg)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"{label} 段が失敗 (発散切り分けは res_nan_*.h5 を見る)")
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(res[-1]), str(run_dir / mesh_h5)],
                       env=_ENV, check=True, capture_output=True, text=True)
        for f in run_dir.glob("res_*"):
            f.unlink()

    if stages == "none":
        (run_dir / "solverConfig.yaml").write_text(cfg_main)
        return run_forge(run_dir)
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", cfg_main)
    soft = soft.replace("convMethod: 1", "convMethod: 0")
    _stage(soft, 3000, "soft")
    if mid_stage and stages == "full":
        mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 1.0, cfl_pseudo: 1.0", cfg_main)
        _stage(mid, 3000, "mid")
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    """軸 M vs 目標 (x_A ≤ x ≤ x_E−0.3 マスク) + 出口一様性 + overshoot。"""
    from ..metrics.extract import axis_mach, exit_uniformity

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    info = json.loads((run_dir / "prepare_info.json").read_text())
    scale = info["scale_m"]
    Md = float(info["Md"])
    x_ach, M_ach = axis_mach(run_dir / "nozzle.h5", res[-1],
                             axis_band=0.9 * scale * 0.1)
    tgt = np.loadtxt(run_dir / "target_axis_M.csv", delimiter=",", skiprows=1)
    ok = np.isfinite(tgt[:, 1])
    xs = np.linspace(max(info["x_A"], tgt[ok, 0].min() / scale) * scale + 1e-9,
                     (info["x_E"] - 0.3) * scale, 200)
    Mt = np.interp(xs, tgt[ok, 0], tgt[ok, 1])
    Ma = np.interp(xs, x_ach, M_ach)
    dM = Ma - Mt
    # overshoot: x_E 以降〜出口手前の軸 M の最大値
    m_dn = x_ach > info["x_E"] * scale
    M_max_dn = float(np.max(M_ach[m_dn])) if m_dn.any() else float("nan")
    try:
        uni = exit_uniformity(run_dir / "nozzle.h5", res[-1], Md,
                              x_d=info["x_E"] * scale, gamma=p.gamma)
    except Exception as e:  # noqa: BLE001 — 一様性は診断 (核心は軸 M)
        uni = {"error": str(e)}
    mf = None
    if info.get("euler_ref"):
        try:
            from ..metrics.deltastar import massflow_ratio
            mf = massflow_ratio(run_dir, info["euler_ref"])
        except Exception as e:  # noqa: BLE001 — 帳簿は診断
            mf = {"error": str(e)}
    out = {"res_file": res[-1].name,
           "massflow": mf,
           "dM_max": float(np.max(np.abs(dM))),
           "dM_max_rel_Md": float(np.max(np.abs(dM)) / Md),
           "dM_rms": float(np.sqrt(np.mean(dM ** 2))),
           "M_axis_exit": float(np.interp(info["x_E"] * scale, x_ach, M_ach)),
           "M_axis_max_downstream": M_max_dn,
           "overshoot_rel": float((M_max_dn - Md) / Md) if np.isfinite(M_max_dn) else None,
           "exit_uniformity": uni,
           "mdot_ratio_moc": info["mdot_ratio_moc"]}
    np.savetxt(run_dir / "achieved_vs_target.csv", np.c_[xs, Mt, Ma],
               delimiter=",", header="x_m,M_target,M_achieved", comments="")
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="axis-Mach チェーン評価 (node Euler)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--cfl", type=float, default=None, help="本段 cfl (YAML evaluate.cfl_main を上書き)")
    ap.add_argument("--implicit-relax", type=float, default=None, help="implicitRelax を deltaT に挿入")
    ap.add_argument("--stages", default="full", choices=("full", "soft", "none"), help="full=soft/mid/本段, soft=soft+本段, none=本段のみ")
    ap.add_argument("--ic-from", default=None)
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, nsteps=a.steps, ic_from=a.ic_from, cfl_main=a.cfl, implicit_relax=a.implicit_relax)
    info["stages"] = a.stages
    (Path(a.run_dir) / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps(info, indent=1, default=str))
    if a.prepare_only:
        return 0
    pp = load_problem(a.problem)
    import time as _t
    t0 = _t.time()
    rc = run_staged(a.run_dir, cfl_main=(a.cfl if a.cfl is not None else pp.evaluate.get("cfl_main")),
                    mid_stage=bool(pp.evaluate.get("mid_stage", pp.is_semiperfect)), stages=a.stages)
    print(f"forge exit={rc} (wall {_t.time() - t0:.0f} s, stages={a.stages})")
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return rc


if __name__ == "__main__":
    sys.exit(main())


# --- A12: 粘性 δ* 補正 (RANS 経路) ------------------------------------------------
def prepare_ns(problem_path, run_dir, nsteps=None, ic_from=None,
               dstar_csv=None, dstar_blend=(6.0, 9.0),
               delta_r_csv=None, offset: str = "normal", euler_ref=None,
               omega: float | None = None, prev_run=None, initializer=None,
               cfl_main: float | None = None, implicit_relax: float | None = None) -> dict:
    r"""**物理壁 (inviscid + δ*) の RANS run** を準備する (A12)。

    plan: plans/active/tooling-nozzle-axismach-viscous-deltastar.md。

    - `dstar_csv` なし → v1: `feedback.deltastar.deltastar_offset` (相関) で法線オフセット
    - `dstar_csv` あり → v2: CSV (x_rt, dstar_rt) を補間して法線オフセット
      (CFD 抽出 δ* の固定点反復)。`dstar_blend=(x_lo, x_hi)` [r_t] で相関→CSV の
      smoothstep 区間を指定 (既定 (6, 9) = case/41 の規約。CSV が全域を覆うなら
      (-1, -0.5) 等で CSV を全域採用 [case/44 v3])
    - メッシュ/段階起動レシピは B8 系 NS v1 (run_0028-0030) で確立したものを流用。
      coarse 中継 (y+~50) は YAML の mesh/evaluate 設定だけの違いで同じ関数で作る
    - **`delta_r_csv` (生産経路, plans/active/tooling-nozzle-deltastar-core-matched-euler.md)**:
      半径方向補正 δ_r(x) [r_t] の CSV (列 x_rt, delta_r; `feedback.deltastar_loop` が作る
      `delta_r_next.csv`) を全域そのまま使い、`offset="radial"` で壁を作る。`dstar_csv`/`dstar_blend`
      (旧 v3 継ぎはぎ) とは排他。`euler_ref` (固定 Euler 参照 run) は帳簿用に prepare_info へ記録。
    """
    from ..feedback.deltastar import _sutherland
    from ..geometry.wall_axismach import PhysicalNozzleWall
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p)
    scale = float(p.spec["r_throat"])
    wall_inv = d["wall_inv"]                      # (n,4) [x,r,th,M] r_t 単位
    # 物理壁 (A13): 上流履歴込み δ* + 真のスロート探索 + 上流 Hermite 再生成。
    # v2 (dstar_csv) は「下流は CFD 抽出、x<6 は履歴相関、[6,9] smoothstep」の合成。
    dstar_x = None
    if dstar_csv is not None:
        tbl = np.loadtxt(dstar_csv, delimiter=",", skiprows=1)
        base = PhysicalNozzleWall(d["wall"], wall_inv, scale, float(p.spec["Pt"]),
                                  float(p.spec["Tt"]), _gam_or_gas(p), p.cp)

        b_lo, b_hi = float(dstar_blend[0]), float(dstar_blend[1])

        def dstar_x(x, _tbl=tbl, _corr=base._dstar_hist):
            x = np.asarray(x, dtype=float)
            w = np.clip((x - b_lo) / max(b_hi - b_lo, 1e-9), 0.0, 1.0)
            w = w * w * (3.0 - 2.0 * w)
            csv = np.interp(np.clip(x, _tbl[0, 0], _tbl[-1, 0]), _tbl[:, 0], _tbl[:, 1])
            return (1.0 - w) * _corr(x) + w * csv
    delta_r_x = None
    init_info = None
    init_cfg = initializer if initializer is not None else p.raw.get("deltastar_initializer")
    if init_cfg and delta_r_csv is None and dstar_csv is None:
        # 積分法初期壁 (plan §4.1): 初回 NS 専用。断熱 / 指定壁温は thermal_bc で。
        from ..feedback.deltastar_integral import integral_bl, delta_r_function
        model = str(init_cfg.get("model", "contur"))
        if model not in ("contur", "contur_momentum_integral"):
            raise ValueError(f"deltastar_initializer.model={model!r} は未対応 (contur のみ)")
        res_init = integral_bl(d["wall"], wall_inv, _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]),
                               scale, thermal_bc=init_cfg.get("thermal_bc"),
                               theta0_m=init_cfg.get("theta0_m"), x_virtual_m=init_cfg.get("x_virtual_m"),
                               a_crocco=float(init_cfg.get("a_crocco", 1.0)), closure=str(init_cfg.get("closure", "contur")))
        # 積分法の出力も同じ 5 次 P-spline で平滑化 (N(Re) テーブルの折れ目などを壁曲率に持ち込まない)
        from ..metrics.deltastar import smooth_delta_quintic
        f_s, sm_diag = smooth_delta_quintic(res_init["x"], res_init["delta_r"], knot_spacing=2.0, lam=1.0)
        res_init["delta_r_raw_integral"] = res_init["delta_r"].copy()
        res_init["delta_r"] = f_s(res_init["x"])
        delta_r_x = delta_r_function(res_init)
        offset = "radial"
        init_info = dict(res_init["settings"])
        init_info["smooth"] = {"kind": "quintic_pspline", **sm_diag}
        init_info["delta_r_throat"] = float(np.interp(0.0, res_init["x"], res_init["delta_r"]))
        init_info["delta_r_exit"] = float(res_init["delta_r"][-1])
        (run_dir / "delta_r_initial.csv").write_text("")   # 後で上書き (run_dir は下で作る)
    if delta_r_csv is not None:
        if dstar_csv is not None:
            raise ValueError("delta_r_csv と dstar_csv は排他")
        tbl_r = np.loadtxt(delta_r_csv, delimiter=",", skiprows=1)
        if not np.all(np.isfinite(tbl_r[:, 1])):
            raise ValueError("delta_r_csv に非有限値がある (deltastar_loop.extract_and_merge で前回値保持済みの CSV を渡す)")
        delta_r_x = lambda x, _t=tbl_r: np.interp(x, _t[:, 0], _t[:, 1])
        offset = "radial"
    wall = PhysicalNozzleWall(d["wall"], wall_inv, scale, float(p.spec["Pt"]),
                              float(p.spec["Tt"]), _gam_or_gas(p), p.cp, dstar_x=dstar_x,
                              offset=offset, delta_r_x=delta_r_x)
    if init_info is not None:
        np.savetxt(run_dir / "delta_r_initial.csv",
                   np.c_[res_init["x"], res_init["delta_r"], res_init["dstar_n"], res_init["theta"] * scale,
                         res_init["H"], res_init["N"], res_init["Cf"], res_init["M"], res_init["Tw"], res_init["Taw"]],
                   delimiter=",", comments="",
                   header="x_rt,delta_r_rt,dstar_n_rt,theta_m,H,N,Cf,M_e,Tw,Taw")
        (run_dir / "delta_r_initial.json").write_text(json.dumps(init_info, indent=1, default=str))
    if init_info is not None:
        dstar_src = f"integral_bl:{init_info['model']} thermal_bc={init_info['thermal_bc']} (radial)"
    elif delta_r_csv is not None:
        dstar_src = f"delta_r_csv:{delta_r_csv} (radial, core-matched Euler ref {euler_ref}, omega {omega})"
    else:
        dstar_src = ("correlation_hist_v1" if dstar_csv is None else f"{dstar_csv} (blend {dstar_blend})") + f" ({offset})"
    msgs = wall.validate()
    if msgs:
        raise ValueError("物理壁フィルタ不合格: " + "; ".join(msgs))
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 561)), nj=int(p.mesh.get("nj", 97)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 4.5e-5)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    xs_p = np.linspace(wall.x_throat, wall.x_e, 1500)
    np.savetxt(run_dir / "wall_physical.csv", np.c_[xs_p * scale, wall.r(xs_p) * scale],
               delimiter=",", header="x_m,r_m", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 48000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    if n % out_int:
        # forge は outStepInterval の倍数でしか res を書かない → 最終 step の res が残るよう n を割り切る間隔に直す
        out_int = next(n // k for k in (3, 2, 4, 6, 1) if n % k == 0)
    cfl_main = float(cfl_main if cfl_main is not None else p.evaluate.get("cfl_main", 1.0))
    implicit_relax = implicit_relax if implicit_relax is not None else p.evaluate.get("implicit_relax")
    (run_dir / "bcondConfig.yaml").write_text(_bcond_with_species(_bcond(p, euler=False), _tp_species_Y(p)))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質は cell 変換コピーで検査 (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node/SST 変換 (config を先に書く — wall_dist は no-slip 壁で作られる)
    cfg_ns = _apply_gas_to_config(_config_sst_node(p, n, out_int, cfl_main), p, run_dir)
    if implicit_relax is not None:
        # 陰解法の緩和 (cfl 6 + implicitRelax 0.7 が生産推奨: case/45 run_0018)。deltaT ブロックに挿入
        cfg_ns = cfg_ns.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    (run_dir / "solverConfig.yaml").write_text(cfg_ns)
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model),
                        h_ref_T=float(p.evaluate.get('thermo_href_temp', 298.15)),
                        species_Y=_tp_species_Y(p))
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        import h5py as _h5
        # ソースが Euler (k/omega~0) のときだけ入口相当を貼る (runner_wt と同じ規約)
        with _h5.File(src) as fs:
            om_src_max = float(np.max(fs["/VALUE/omega"][:])) if "/VALUE/omega" in fs else 0.0
        if om_src_max < 1.0:
            with _h5.File(run_dir / "nozzle.h5", "r+") as f:
                ro = f["/VALUE/ro"][:]
                f["/VALUE/roK"][:] = ro * 1.0
                f["/VALUE/roOmega"][:] = ro * 18000.0
        # 近壁 omega の粘性底層フロア omega = 6 nu / (beta1 y^2) (run_0030 の教訓)
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
    info = {"chain": "axismach_ns", "discretization": "node", "viscous": True,
            "dstar_source": dstar_src, "mdot_ratio_moc": None,
            "throat_physical": {"x": wall.x_throat, "r": wall.r_throat,
                                "kappa": wall.kappa_throat,
                                "dstar_throat_correlation": float(wall._dstar_hist(0.0)),
                                "delta_r_throat_applied": float(wall.r_throat - 1.0)},
            "offset": wall.offset_mode, "euler_ref": (str(euler_ref) if euler_ref else None),
            "initializer": init_info,
            "omega": omega, "prev_run": (str(prev_run) if prev_run else None),
            "delta_r_csv": (str(delta_r_csv) if delta_r_csv else None),
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "Lc_mode": d["Lc_mode"], "Lc_solve": d["Lc_solve"],
            "anchor": list(d["anchor"]), "anchor_source": d["anchor_source"],
            "start_line": d["start_line"], "wall_mode": d["wall_mode"],
            "Md": d["Md"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "nStepOuter": n, "cfl_main": cfl_main, "implicit_relax": implicit_relax, "scale_m": scale,
            "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged_ns(run_dir, stages: str = "full", ramp=None, ramp_steps: int = 1000) -> int:
    """NS の起動。stages:
    - "full" (既定・run_0030 レシピ): soft (1次 cfl0.5 ni10, 3000 step) → mid (1次 cfl1, 3000) → 本段 (2次 cfl_main)。
    - "none": 本段だけ (収束済み NS 場からの warm start 用)。
    - "ramp": 2 次のまま cfl を `ramp` (例 (1, 2, 3.5)) の順に各 ramp_steps だけ回して本段 cfl_main へ
      (forge に CFL ランプ機能は無いので restart で段階化する。本段の step 数はランプ分を差し引く)。
    各段の最終場を IC に引き継ぐ。"""
    import re
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    n_main = int(re.search(r"nStepOuter: (\d+)", cfg_main).group(1))

    def _stage(cfg, nsteps):
        cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", cfg)
        cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {nsteps}", cfg)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"段階起動が失敗 (rc={rc}, res={res[-1].name if res else None})")
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(res[-1]), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        for f in run_dir.glob("res_*"):
            f.unlink()

    if stages == "full":
        soft = cfg_main
        soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", soft)
        soft = soft.replace("convMethod: 1", "convMethod: 0")
        soft = soft.replace("nStepInner: 5", "nStepInner: 10")
        _stage(soft, 3000)
        mid = cfg_main
        mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 1.0, cfl_pseudo: 1.0", mid)
        mid = mid.replace("convMethod: 1", "convMethod: 0")
        mid = mid.replace("nStepInner: 5", "nStepInner: 10")
        _stage(mid, 3000)
    elif stages == "ramp":
        ramp = tuple(ramp or (1.0, 2.0, 3.5))
        for c in ramp:
            st = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {c}, cfl_pseudo: {c}", cfg_main)
            _stage(st, int(ramp_steps))
        n_main = max(n_main - int(ramp_steps) * len(ramp), int(ramp_steps))
        # 最終 res が書かれるよう outStepInterval の倍数に丸める (forge は outStepInterval の倍数でしか res を書かない)
        out_int = int(re.search(r"outStepInterval: (\d+)", cfg_main).group(1))
        n_main = max((n_main // out_int) * out_int, out_int)
        cfg_main = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n_main}", cfg_main)
    elif stages != "none":
        raise ValueError("stages は 'full' / 'none' / 'ramp'")
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)
