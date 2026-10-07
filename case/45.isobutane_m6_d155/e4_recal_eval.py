"""plan verification-case45-euler-total-enthalpy §6 E4 (2026-10-07 登録 99431498。結果を見る前に実装) の評価器: 出口較正のやり直し。
段 1 = run_0163_euler_e4_recal_d0 (δ₀ = Md_moc_offset +3.770e-4)、段 2 = run_0164_euler_e4_recal_d1 (段 1 が「更新」のときだけ、δ₁)。
起動は run_e4_recal.sh (準備・実行は e4_recal.py)。この評価器は保存場を読むだけで、forge は起動しない。

固定 (登録): 単調壁 (wall_fit_mono_r2 [0, 1.5])・MOC は legacy + fixed2 (MOC のキー無し)・凍結の初期線 run_0062 (res_6000)・r_t・BC・
  熱物性・バイナリ (forge の sha256 = run_0143 の RUN_PROVENANCE; E2 と同じ)。格子は E3 の mesh_euler (2000 × 97・全域 0.005・スロートの
  別指定なし・cap なし)。IC は新しい格子の上の等エントロピー IC (起動前に全節点 |T0 − 1600| ≤ 1 K; e4_recal.py prep が E4_PREP.json に記録)。
長さ: soft 3000 (1 次・cfl 0.5) + 本段 54000 (2 次・cfl 2・implicitRelax 0.7)、出力 1000 ごと。判定窓は本段 42000〜54000 の 13 枚と
  50000〜54000 の 5 枚 (予算であり収束の予測ではない)。本段の step 数・出力間隔・窓は下の定数から作り、実行側 (e4_recal.py) も同じ定数を読む。

出口の評価量 (登録): 各設計の実際の最終断面 (i = ni − 1) で、η = r/r_w ∈ [0.05, 0.7] について、NS の基準格子 G1 の最終断面の η の列
  (固定; NS の問題 problem_d155_ns_finemesh_recal_final_mono.yaml の mesh から mesh2d で求めた 12 点 0.08876〜0.67224) に、自格子の
  M(η) を同じ線形補間 (np.interp) で移して単純平均した M_common。自格子の帯内の節点の単純平均 M_own も別の列に残す
  (nozzle_report.metrics・exitM_sampling_ab の腕 A と同じ)。M = |U|/sonic (res の Ux・Uy・sonic; nozzle_report.load_field と同じ)。
  軸上の M や x_E の値に置き換えない。η の列は登録の数 (ETA_COMMON_N) と sha256 (ETA_COMMON_SHA256) と照合し、AWS に NS の基準 run
  (run_0147_ns_mono_final) の nozzle.h5 があれば、その最終断面の η (float32 の座標) とも照合する (|Δη| ≤ 1e-6 かつ帯内の点の数が同じ)。

前提 (登録):
  P1 残差: 本段区間 (stage_manifest の最後の区間 = main) の check_convergence が PASS、または全不合格列が停滞 (STALLED/plateau) だけ
     (E2 の評価器 parse_convergence の列ごとの内訳)。停滞のみは今回の実務の較正に限って許し、NOT CONVERGED を記録に残す。RISING・still
     converging・DIVERGED・入力不備は保留 (前提の未達)。
  P2 全温の健全性 (E2 の評価器の主指標と幅の条件を流用): 全温の復元は E1・E2 と同じ (euler_t0_stage_ab.recon)。評価点はこの run の節点
     (E2 の評価点の作り方で A = B = 自格子; 補間値 = 節点値)。登録の 9 領域で、(a) 判定窓 13 枚すべてで全領域の T0 − 1600 が
     最大 ≤ +1 K かつ最小 ≥ −1 K、(b) 27 列 (9 領域 × 最大・最小・|偏差| の 99 % 点) の時間方向の最大 − 最小 ≤ 0.1 K を 13 枚・5 枚の
     両窓で (E2 の window_check。check_quasisteady の VERDICT は併記)。
  P3 出口 M: M_common の系列について、13 枚・5 枚の両窓で check_quasisteady (--series-csv --tail 1) が STEADY、かつ 13 枚の幅
     (最大 − 最小) ≤ 5e-5、かつ |末尾 5 枚 (50000〜54000) の平均 − その直前 5 枚 (45000〜49000) の平均| ≤ 5e-5。
  固定の条件: 上の「固定」と、本段の step 数・出力間隔の実効値 (solverConfig.yaml・stage_manifest) と残差の最終 step、段の cfl
     (stage_manifest の soft 0.5・main 2)、段 2 では段 1 との bcondConfig・solverConfig・probe・species_meta の同一と Md_moc_offset = δ₁。
判定 (登録):
  段 1: 前提を満たし、13 枚すべてで |M_common − 6| ≤ 1e-4 → 据え置き (δ₀)。外れたら δ₁ = δ₀ − (平均 M_common − 6) (係数 1、#11f と同じ式)
        を出し、段 2 を 1 回だけ。
  段 2: 段 1 が「更新」で、run_0164 の Md_moc_offset が段 1 の δ₁ と一致し、前提を満たし、13 枚すべてで |M_common − 6| ≤ 1e-4 → 合格
        (δ₁ を採用)。外れたら補正を重ねず保留 (係数 1 の 1 回補正の則を棄却)。
  前提の未達は判別不能 (保留)。窓や補正の回数を変えて救済しない。
  実装時の読み (結果を見る前): 「平均 M_common」は判定窓 13 枚の平均 (DELTA1_BASIS; 末尾 5 枚の平均と、それによる δ₁ も記録する)。
  比較は float の ≤ そのまま (許容を足さない; E2 と同じ)。|·| の 1e-4 は M_common − 6 に当てる。
出力: `_band_ab/e4_recal_eval.json` (評価器の sha256・登録の commit 99431498・前提の定義と実効の判定・各時点の出口 M・全温の主指標)、
  各 run の `e4_recal_series.csv` (全時点の出口 M と全温の主指標)・`e4_recal_exitM_{win13,tail5}.csv`/`.txt`・`e4_recal_t0_{win13,tail5}.csv`/`.txt`。
usage: python3 e4_recal_eval.py [case_dir] [--out PATH]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda/tools"
sys.path.insert(0, str(HERE))
import euler_t0_e2_eval as EV2  # noqa: E402  (全温の主指標・幅の条件・残差の列ごとの内訳・check_quasisteady の読みを流用する)

PLAN = "plans/active/verification-case45-euler-total-enthalpy.md §6 E4"
PLAN_REG_COMMIT = "99431498"
STAGES = ("d0", "d1")
RUNS = {"d0": "run_0163_euler_e4_recal_d0", "d1": "run_0164_euler_e4_recal_d1"}
PROBLEMS = {"d0": "problem_d155_euler_e4_recal_d0.yaml", "d1": "problem_d155_euler_e4_recal_d1.yaml"}
DELTA0 = 3.770e-4                                # 開始値 δ₀ (登録: Md_moc_offset = +3.770e-4)
M_TARGET = 6.0
M_TOL = 1e-4                                     # 13 枚すべてで |M_common − 6| ≤ 1e-4
EXIT_RANGE_TOL = 5e-5                            # 13 枚の幅
EXIT_MEANDIFF_TOL = 5e-5                         # |末尾 5 枚の平均 − 直前 5 枚の平均|
ETA_BAND = (0.05, 0.7)
NS_REF_PROBLEM = "problem_d155_ns_finemesh_recal_final_mono.yaml"   # NS の基準格子 G1 の問題 (mesh ブロック)
NS_REF_RUN = "run_0147_ns_mono_final"            # その NS の run (AWS にあれば最終断面の η と照合)
ETA_COMMON_N = 12                                # 登録時に求めた η の列の数と sha256 (float64 のバイト列)。食い違えば固定の条件の不成立
ETA_COMMON_SHA256 = "96a495c2e80e57da7d5e43b6e5225aa9906eb12d978ad27056e0a0d1b66b8223"
ETA_XCHECK_TOL = 1e-6                            # run_0147 の float32 の座標から求めた η との照合の許容 (登録時のローカルの NS prep で 3.8e-8)
DELTA1_BASIS = "win13"                           # 実装時の読み: δ₁ の「平均 M_common」は判定窓 13 枚の平均
SOFT_STEPS = 3000
MAIN_NSTEPS = 54000
OUT_INTERVAL = 1000
MAIN_STEPS = tuple(range(0, MAIN_NSTEPS + 1, OUT_INTERVAL))                        # 本段の出力 (55 枚)
WIN13 = tuple(range(MAIN_NSTEPS - 12 * OUT_INTERVAL, MAIN_NSTEPS + 1, OUT_INTERVAL))  # 判定窓 42000〜54000 (13 枚)
TAIL5 = tuple(range(MAIN_NSTEPS - 4 * OUT_INTERVAL, MAIN_NSTEPS + 1, OUT_INTERVAL))   # 末尾 50000〜54000 (5 枚)
PREV5 = WIN13[-10:-5]                                                              # 末尾 5 枚の直前 5 枚 45000〜49000
WINDOWS = {"win13": WIN13, "tail5": TAIL5}
CFL_SOFT, CFL_MAIN, RELAX = 0.5, 2.0, 0.7
REF_RUN_BIN = "run_0143_euler_wallfit_monoG1_r1"  # forge の sha256 の基準 (E2 と同じ)
TT_REG = 1600.0
SCALE_M = 0.0768075
SOFT_DIR = "_soft_stage"
PREP_JSON = "E4_PREP.json"
SEGMENT_VERDICT_FILE = "CONVERGENCE_VERDICT_segment.txt"
SEGMENT_MAIN = "main"
STAGES_EXPECTED = ["soft", "main"]
MOC_EXPECT = {"axis_limit": "legacy", "corrector": "fixed2"}
MONO_R2 = [0.0, 1.5]
INITIAL_LINE = ("run_0062_euler_wallfit_fit_r1_ext6k", "res_6000.h5")
MESH_EULER_EXPECT = {"ni": 2000, "nj": 97, "wall_first_frac": 0.005, "wall_first_frac_throat": None, "axis_cap_frac": None,
                     "axis_gap_frac": None, "throat_refine": 4.0, "throat_width": 3.0}
SAME_AS_D0 = ("bcondConfig.yaml", "solverConfig.yaml", "probe.yaml", "species_meta.yaml")   # 段 2 の prep が段 1 の prep と同一であるべきファイル
SERIES_CSV = "e4_recal_series.csv"
OUT_JSON = "_band_ab/e4_recal_eval.json"
EXIT_COLS = ("exitM_common", "exitM_own")

LBL_KEEP = "据え置き: 13 枚すべてで |M_common − 6| ≤ 1e-4 → δ₀ = +3.770e-4 を据え置く"
LBL_UPDATE = "更新: δ₁ = δ₀ − (平均 M_common − 6) で 1 回だけ更新し、作り直した壁の独立の run (段 2) で同じ条件を確かめる"
LBL_PASS2 = "合格: δ₁ の壁で 13 枚すべて |M_common − 6| ≤ 1e-4 → δ₁ を採用"
LBL_HOLD2 = "保留: δ₁ の壁でも |M_common − 6| ≤ 1e-4 を満たさない — 補正を重ねない (係数 1 の 1 回補正の則を棄却)"
LBL_UNDET = "判別不能 (前提の未達・保留)"
LBL_NOTRUN = "未実施"
SCOPE_NOTE = ("Euler の出口較正 (登録の計算手順・固定の予算 soft 3000 + 本段 54000・窓 42000〜54000) についての判定。停滞のみの NOT CONVERGED は"
              "収束の証明ではない。Euler の較正の合格を NS に移せるとは限らない (配点も粘性も違う) — NS の出口 M は NS で判定する")


def _sha(p) -> str | None:
    return EV2._sha(p)


# --- 出口の共通の η の列 -------------------------------------------------------------------------------------------------------
class _FlatWall:
    """η = r/r_w を取り出すための仮の壁 (r ≡ 1)。半径方向の配点は x_e ≥ wall_first_blend_x1 の断面では x に依らない。"""
    x_in, x_e = -20.0, 95.0

    def r(self, x):
        return np.ones_like(np.asarray(x, dtype=float))


def eta_common_from_problem(problem: Path) -> np.ndarray:
    """NS の問題の mesh (NS の経路 runner_axismach.mesh_params) から mesh2d で最終断面の η を作り、η ∈ [0.05, 0.7] を返す (float64)。"""
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.evaluate import runner_axismach as RA
    from forge_design.meshing.mesh2d import generate_axisym_mesh
    p = RA.load_problem(problem)
    mp = RA.mesh_params(p, 1.0, 561, 97, 4.5e-5)
    if mp.axis_gap_frac is not None or mp.axis_cap_frac is not None:
        raise ValueError("NS の基準格子に軸側の cap がある (登録の G1 でない)")
    if not (_FlatWall.x_e >= mp.wall_first_blend_x1 and (mp.wall_first_up_x1 is None or _FlatWall.x_e > mp.wall_first_up_x1)):
        raise ValueError("仮の壁の出口がブレンドの区間の外にない")
    c, _, _ = generate_axisym_mesh(_FlatWall(), mp)
    R = c[:, 1].reshape(mp.ni, mp.nj)
    eta = R[-1] / R[-1, -1]
    return np.ascontiguousarray(eta[(eta >= ETA_BAND[0]) & (eta <= ETA_BAND[1])], dtype=np.float64)


def eta_digest(eta) -> str:
    return hashlib.sha256(np.ascontiguousarray(eta, dtype=np.float64).tobytes()).hexdigest()


def final_section_eta(run_dir: Path) -> np.ndarray:
    """run の nozzle.h5 の最終断面の η (float32 の座標から float64 で)。並びは VIZMESH/CONNE で確かめる。"""
    import h5py
    from ic_index_map import structured_shape
    with h5py.File(run_dir / "nozzle.h5", "r") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3).astype(np.float64)
        ni, nj = structured_shape(f["VIZMESH/CONNE"], nc.shape[0])
    R = nc[:, 1].reshape(ni, nj)
    return R[-1] / R[-1, -1]


def eta_common(case: Path) -> tuple:
    """(η の列, 記録)。登録の数・sha256 と照合し、NS の基準 run があればその最終断面とも照合する。不成立は記録の problems に入れる。"""
    rec = {"problem": NS_REF_PROBLEM, "problem_sha256": _sha(case / NS_REF_PROBLEM), "band": list(ETA_BAND), "problems": []}
    try:
        eta = eta_common_from_problem(case / NS_REF_PROBLEM)
    except Exception as e:  # noqa: BLE001
        rec["problems"].append(f"η の列を作れない: {type(e).__name__}: {e}")
        return None, rec
    rec.update(n=int(eta.size), sha256=eta_digest(eta), eta=[float(v) for v in eta])
    if eta.size != ETA_COMMON_N or rec["sha256"] != ETA_COMMON_SHA256:
        rec["problems"].append(f"η の列 ({eta.size} 点、{rec['sha256'][:16]}) が登録 ({ETA_COMMON_N} 点、{ETA_COMMON_SHA256[:16]}) と違う")
    ref = case / NS_REF_RUN
    if (ref / "nozzle.h5").is_file():
        try:
            er = final_section_eta(ref)
            m = (er >= ETA_BAND[0]) & (er <= ETA_BAND[1])
            x = {"run": NS_REF_RUN, "n_in_band": int(m.sum())}
            if int(m.sum()) == eta.size:
                x["max_abs_deta"] = float(np.abs(er[m] - eta).max())
            x["ok"] = bool(int(m.sum()) == eta.size and x.get("max_abs_deta", np.inf) <= ETA_XCHECK_TOL)
        except Exception as e:  # noqa: BLE001
            x = {"run": NS_REF_RUN, "ok": False, "error": f"{type(e).__name__}: {e}"}
        rec["xcheck_ns_run"] = x
        if not x["ok"]:
            rec["problems"].append(f"{NS_REF_RUN} の最終断面の η と照合が不成立: {x}")
    else:
        rec["xcheck_ns_run"] = {"run": NS_REF_RUN, "status": "未確認 (nozzle.h5 が無い)"}
    return eta, rec


# --- 出口 M -----------------------------------------------------------------------------------------------------------------------
def exit_values(M_col, eta_col, eta_c) -> dict:
    """最終断面の M(η) (自格子の節点、η は 0〜1 の狭義単調増加) → M_common (共通の η の列への線形補間の単純平均) と M_own
    (帯内の自格子の節点の単純平均)。"""
    M_col = np.asarray(M_col, dtype=np.float64)
    e = np.asarray(eta_col, dtype=np.float64)
    if not (e.ndim == 1 and e.size == M_col.size and e[0] == 0.0 and np.all(np.diff(e) > 0)):
        raise ValueError("最終断面の η が 0 から始まる狭義単調増加でない")
    m = (e >= ETA_BAND[0]) & (e <= ETA_BAND[1])
    return {"exitM_common": float(np.mean(np.interp(eta_c, e, M_col))), "exitM_own": float(np.mean(M_col[m])), "n_own": int(m.sum())}


def exit_snapshot(h5: Path, ni: int, nj: int, eta_last, eta_c) -> dict:
    import h5py
    rec = {"file": h5.name, "problems": []}
    if not h5.is_file():
        rec["problems"].append("ファイルが無い")
        return rec
    try:
        with h5py.File(h5, "r") as f:
            U = {k: f[f"/VALUE/{k}"][:].astype(np.float64) for k in ("Ux", "Uy", "sonic")}
        if U["Ux"].size != ni * nj:
            raise ValueError(f"節点数 {U['Ux'].size} が格子 {ni}×{nj} と違う")
        M = (np.hypot(U["Ux"], U["Uy"]) / U["sonic"]).reshape(ni, nj)[-1]
        if not np.all(np.isfinite(M)):
            raise ValueError("最終断面の M に非有限値")
        rec.update(exit_values(M, eta_last, eta_c))
    except Exception as e:  # noqa: BLE001
        rec["problems"].append(f"{type(e).__name__}: {e}")
    return rec


def exit_window_check(series: dict, qs_texts: dict) -> dict:
    """P3: series = {step: M_common}。両窓の check_quasisteady (qs_texts = {窓: 出力}) が exitM_common で STEADY、13 枚の幅 ≤ 5e-5、
    |mean(末尾 5) − mean(直前 5)| ≤ 5e-5。時点の欠け・非有限は不成立。"""
    rec = {"problems": [], "verdicts": {}}
    miss = [s for s in WIN13 if s not in series or not math.isfinite(series[s])]
    if miss:
        rec["problems"].append(f"判定窓の時点が欠ける・非有限: {miss}")
    for wn, txt in qs_texts.items():
        v = EV2.parse_quasisteady(txt or "", ["exitM_common"])["exitM_common"]
        rec["verdicts"][wn] = v
        m = re.search(r"\[(\d+) rows, steps", txt or "")
        if not m or int(m.group(1)) != len(WINDOWS[wn]):
            rec["problems"].append(f"{wn}: check_quasisteady の行数 {m.group(1) if m else '読めない'} が {len(WINDOWS[wn])} でない")
        if v != "STEADY":
            rec["problems"].append(f"{wn}: check_quasisteady の VERDICT が {v} (STEADY でない)")
    if not miss:
        w = np.array([series[s] for s in WIN13])
        rec["range13"] = float(w.max() - w.min())
        rec["mean_tail5"] = float(np.mean([series[s] for s in TAIL5]))
        rec["mean_prev5"] = float(np.mean([series[s] for s in PREV5]))
        rec["meandiff_tail5_prev5"] = rec["mean_tail5"] - rec["mean_prev5"]
        if not rec["range13"] <= EXIT_RANGE_TOL:
            rec["problems"].append(f"13 枚の幅 {rec['range13']:.3e} > {EXIT_RANGE_TOL:g}")
        if not abs(rec["meandiff_tail5_prev5"]) <= EXIT_MEANDIFF_TOL:
            rec["problems"].append(f"|末尾 5 枚の平均 − 直前 5 枚の平均| = {abs(rec['meandiff_tail5_prev5']):.3e} > {EXIT_MEANDIFF_TOL:g}")
    rec["ok"] = not rec["problems"]
    return rec


# --- 全温の健全性 (E2 の主指標と幅の条件) ----------------------------------------------------------------------------------------
def t0_within_check(ind_by_step: dict) -> dict:
    """P2 (a): 判定窓 13 枚すべてで全 9 領域の主指標が最大 ≤ +1 K かつ最小 ≥ −1 K (E2 の within_K)。"""
    rec = {"problems": [], "per_step": {}}
    for s in WIN13:
        ind = ind_by_step.get(s)
        if ind is None or not EV2._complete(ind):
            rec["problems"].append(f"判定窓の {s} の主指標が欠ける・非有限")
            continue
        rec["per_step"][s] = EV2.within_K(ind)
    if rec["per_step"]:
        rec["max_over_steps_K"] = max(rec["per_step"].values())
        bad = [s for s, v in rec["per_step"].items() if not v <= EV2.WITHIN_K]
        if bad:
            rec["problems"].append(f"全領域の |T0 − 1600| > {EV2.WITHIN_K:g} K の時点 {bad} (最大 {rec['max_over_steps_K']:.4g} K)")
    rec["ok"] = not rec["problems"]
    return rec


# --- 判定 ---------------------------------------------------------------------------------------------------------------------------
def _window_values(series: dict) -> list | None:
    if any(s not in series or not math.isfinite(series[s]) for s in WIN13):
        return None
    return [series[s] for s in WIN13]


def delta1_of(series: dict, delta0: float = DELTA0) -> dict:
    """δ₁ = δ₀ − (平均 M_common − 6)。平均は DELTA1_BASIS (判定窓 13 枚)、末尾 5 枚の平均による値も記録する。"""
    w = _window_values(series)
    m13 = float(np.mean(w))
    m5 = float(np.mean([series[s] for s in TAIL5]))
    d13, d5 = delta0 - (m13 - M_TARGET), delta0 - (m5 - M_TARGET)
    d = d13 if DELTA1_BASIS == "win13" else d5
    return {"basis": DELTA1_BASIS, "mean_win13": m13, "mean_tail5": m5, "delta1": d, "delta1_repr": repr(d),
            "delta1_if_tail5": d5, "delta0": delta0}


def within_target(series: dict) -> tuple:
    """(13 枚すべてで |M_common − 6| ≤ 1e-4 か, 各時点の |M_common − 6|)。"""
    w = _window_values(series)
    dev = {s: abs(series[s] - M_TARGET) for s in WIN13}
    return all(v <= M_TOL for v in dev.values()) and w is not None, dev


def judge_stage1(series: dict, pre_ok: bool, pre_reasons=()) -> dict:
    """段 1 (δ₀)。series = {本段の step: M_common}。"""
    if _window_values(series) is None:
        return {"verdict": LBL_UNDET, "reasons": list(pre_reasons) + ["判定窓の M_common が欠ける・非有限"], "facts": {}}
    ok, dev = within_target(series)
    facts = {"abs_dev_by_step": dev, "max_abs_dev": max(dev.values()), "all_within_1e-4": ok, **delta1_of(series)}
    if not pre_ok:
        return {"verdict": LBL_UNDET, "reasons": ["前提の未達"] + list(pre_reasons), "facts": facts}
    if ok:
        return {"verdict": LBL_KEEP, "reasons": [], "facts": facts, "adopted_delta": DELTA0}
    return {"verdict": LBL_UPDATE, "reasons": [f"|M_common − 6| の最大 {facts['max_abs_dev']:.3e} > {M_TOL:g}"], "facts": facts,
            "delta1": facts["delta1"], "delta1_repr": facts["delta1_repr"]}


def judge_stage2(series: dict, pre_ok: bool, pre_reasons, stage1: dict, delta_run) -> dict:
    """段 2 (δ₁)。段 1 が「更新」で、run の Md_moc_offset (delta_run) が段 1 の δ₁ と一致することを要求する。"""
    why = list(pre_reasons)
    if (stage1 or {}).get("verdict") != LBL_UPDATE:
        return {"verdict": LBL_UNDET, "reasons": why + [f"段 1 が「更新」でない ({(stage1 or {}).get('verdict')!r}) — 段 2 は登録外"], "facts": {}}
    if delta_run is None or float(delta_run) != float(stage1["delta1"]):
        why.append(f"run の Md_moc_offset {delta_run!r} が段 1 の δ₁ {stage1['delta1_repr']} と一致しない")
        pre_ok = False
    if _window_values(series) is None:
        return {"verdict": LBL_UNDET, "reasons": why + ["判定窓の M_common が欠ける・非有限"], "facts": {}}
    ok, dev = within_target(series)
    facts = {"abs_dev_by_step": dev, "max_abs_dev": max(dev.values()), "all_within_1e-4": ok,
             "mean_win13": float(np.mean(_window_values(series))), "delta1": stage1["delta1"]}
    if not pre_ok:
        return {"verdict": LBL_UNDET, "reasons": ["前提の未達"] + why, "facts": facts}
    if ok:
        return {"verdict": LBL_PASS2, "reasons": [], "facts": facts, "adopted_delta": stage1["delta1"]}
    return {"verdict": LBL_HOLD2, "reasons": [f"|M_common − 6| の最大 {facts['max_abs_dev']:.3e} > {M_TOL:g}"], "facts": facts}


# --- 本段の step 数の同期 (実行側 e4_recal.py と共有) ------------------------------------------------------------------------------
def window_problems() -> list:
    bad = []
    if MAIN_STEPS != tuple(range(0, MAIN_NSTEPS + 1, OUT_INTERVAL)):
        bad.append(f"MAIN_STEPS が 0〜{MAIN_NSTEPS} の {OUT_INTERVAL} ごとでない")
    if len(WIN13) != 13 or WIN13[-1] != MAIN_NSTEPS or any(b - a != OUT_INTERVAL for a, b in zip(WIN13, WIN13[1:])):
        bad.append(f"判定窓 {WIN13[0]}〜{WIN13[-1]} が本段の最後の 13 枚でない (本段 {MAIN_NSTEPS})")
    if len(TAIL5) != 5 or TAIL5 != WIN13[-5:] or len(PREV5) != 5 or PREV5 != WIN13[-10:-5]:
        bad.append("末尾 5 枚・直前 5 枚の窓が判定窓の最後の 10 枚でない")
    return bad


def config_steps(cfg_text: str) -> dict:
    return EV2.config_steps(cfg_text)


def step_problems(eff: dict) -> list:
    bad = list(window_problems())
    if eff.get("nStepOuter") != MAIN_NSTEPS:
        bad.append(f"本段の nStepOuter {eff.get('nStepOuter')!r} が評価器の MAIN_NSTEPS {MAIN_NSTEPS} と違う")
    if eff.get("outStepInterval") != OUT_INTERVAL:
        bad.append(f"本段の outStepInterval {eff.get('outStepInterval')!r} が評価器の OUT_INTERVAL {OUT_INTERVAL} と違う")
    if "manifest_main_nStepOuter" in eff and eff["manifest_main_nStepOuter"] != MAIN_NSTEPS:
        bad.append(f"stage_manifest の main の nStepOuter {eff['manifest_main_nStepOuter']!r} が {MAIN_NSTEPS} と違う")
    if "residual_last_step" in eff and eff["residual_last_step"] != MAIN_NSTEPS - 1:
        bad.append(f"本段の残差の最終 step {eff['residual_last_step']!r} が {MAIN_NSTEPS - 1} でない")
    return bad


def _num(v):
    m = re.match(r"[-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?", str(v if v is not None else ""))
    return float(m.group(0)) if m else None


def manifest_problems(rd: Path) -> tuple:
    """stage_manifest の段 (soft → main) と段の cfl・cfl_pseudo (soft 0.5・main 2)・本段の implicitRelax 0.7。戻り (記録, 不成立の理由)。"""
    bad = []
    try:
        man = json.loads((rd / "stage_manifest.json").read_text())
        st = man.get("stages") if isinstance(man, dict) else man
        tags = [s.get("tag") for s in st]
        soft = {s.get("tag"): (s.get("soft") or {}) for s in st}
    except (OSError, ValueError, AttributeError, TypeError) as e:
        return {"error": f"{type(e).__name__}: {e}"}, [f"stage_manifest を読めない ({type(e).__name__})"]
    if tags != STAGES_EXPECTED:
        bad.append(f"stage_manifest の段が {tags} ({STAGES_EXPECTED} でない)")
    for tag, cfl in (("soft", CFL_SOFT), ("main", CFL_MAIN)):
        for k in ("cfl", "cfl_pseudo"):
            if _num(soft.get(tag, {}).get(k)) != cfl:
                bad.append(f"段 {tag} の {k} {soft.get(tag, {}).get(k)!r} が {cfl} でない")
    if _num(soft.get("main", {}).get("implicitRelax")) != RELAX:
        bad.append(f"本段の implicitRelax {soft.get('main', {}).get('implicitRelax')!r} が {RELAX} でない")
    return {"tags": tags, "soft": soft}, bad


def run_steps(rd: Path, conv: dict | None = None) -> dict:
    return EV2.run_steps(rd, conv)


# --- 固定の条件 ---------------------------------------------------------------------------------------------------------------------
def _psha(p: Path):
    return EV2._psha(p)


def prepare_info_problems(info: dict, expected_delta) -> list:
    """prepare_info.json の固定の条件: MOC legacy + fixed2・mono_r2・mesh_euler の実効値・Md_moc_offset・凍結の初期線・r_t。"""
    bad = []
    moc = info.get("moc") or {}
    if {k: moc.get(k) for k in MOC_EXPECT} != MOC_EXPECT:
        bad.append(f"MOC が {MOC_EXPECT} でない ({ {k: moc.get(k) for k in MOC_EXPECT} })")
    from throat_mono_judge import mono_r2_matches
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf or not mono_r2_matches(wf["mono_r2"], MONO_R2):
        bad.append(f"壁の mono_r2 {wf.get('mono_r2')!r} が {MONO_R2} でない")
    mesh = info.get("mesh") or {}
    prm = mesh.get("params") or {}
    if mesh.get("source") != "mesh_euler":
        bad.append(f"格子の採用元が {mesh.get('source')!r} (mesh_euler でない)")
    got = {k: prm.get(k) for k in MESH_EULER_EXPECT}
    if got != MESH_EULER_EXPECT:
        bad.append(f"格子の実効値 {got} が登録 {MESH_EULER_EXPECT} と違う")
    if expected_delta is None or info.get("Md_moc_offset") is None or float(info["Md_moc_offset"]) != float(expected_delta):
        bad.append(f"Md_moc_offset {info.get('Md_moc_offset')!r} が {expected_delta!r} でない")
    il = info.get("initial_line") or {}
    res = il.get("res")
    res = res if isinstance(res, list) else [res]                  # prepare_info は snapshot の列 (['res_6000.h5']) で持つ
    if Path(str(il.get("run") or "")).name != INITIAL_LINE[0] or res != [INITIAL_LINE[1]]:
        bad.append(f"凍結の初期線 {il.get('run')!r}/{il.get('res')!r} が {INITIAL_LINE} でない")
    if info.get("scale_m") is None or float(info["scale_m"]) != SCALE_M:
        bad.append(f"r_t (scale_m) {info.get('scale_m')!r} が {SCALE_M} でない")
    return bad


def fixed_conditions(case: Path, stage: str, expected_delta) -> dict:
    rd = case / RUNS[stage]
    probs, rec = [], {}
    ref = _psha(case / REF_RUN_BIN / "RUN_PROVENANCE.txt")
    shas = {"main": _psha(rd / "RUN_PROVENANCE.txt"), "soft": _psha(rd / SOFT_DIR / "RUN_PROVENANCE.txt")}
    rec.update(forge_sha256=shas, forge_sha256_ref=ref)
    if not ref:
        probs.append(f"{REF_RUN_BIN}/RUN_PROVENANCE.txt の forge_sha256 を読めない")
    bad = {k: v for k, v in shas.items() if not v or v != ref}
    if bad:
        probs.append(f"forge の sha256 が {REF_RUN_BIN} ({ref}) と違う・読めない: {bad}")
    man, mp = manifest_problems(rd)
    rec["manifest"] = man
    probs += mp
    try:
        pr = json.loads((rd / PREP_JSON).read_text())
    except (OSError, ValueError):
        pr = {}
    rec["ic_check"] = (pr.get("ic_check") or {}).get("VERDICT")
    if rec["ic_check"] != "OK":
        probs.append(f"起動前の IC の検査 ({PREP_JSON}) が OK でない ({rec['ic_check']!r})")
    if pr.get("stage") != stage:
        probs.append(f"{PREP_JSON} の段 {pr.get('stage')!r} が {stage} でない")
    try:
        import moc_v5c_thermo_ab as V
        rec["Tt"] = V.read_Tt(rd)
    except Exception as e:  # noqa: BLE001
        rec["Tt"] = f"読めない: {type(e).__name__}"
    if rec["Tt"] != TT_REG:
        probs.append(f"入口の Tt が {rec['Tt']!r} ({TT_REG} でない)")
    try:
        info = json.loads((rd / "prepare_info.json").read_text())
    except (OSError, ValueError) as e:
        info = {}
        probs.append(f"prepare_info.json を読めない ({type(e).__name__})")
    rec["Md_moc_offset"] = info.get("Md_moc_offset")
    probs += prepare_info_problems(info, expected_delta)
    rec["files_sha256"] = {n: _sha(rd / n) for n in SAME_AS_D0}
    if stage == "d1":
        try:
            p0 = json.loads((case / RUNS["d0"] / PREP_JSON).read_text()).get("files_sha256") or {}
            p1 = pr.get("files_sha256") or {}
            diff = [n for n in SAME_AS_D0 if (p0.get(n) or p1.get(n)) and p0.get(n) != p1.get(n)]
            rec["same_as_d0"] = {"compared": list(SAME_AS_D0), "differ": diff}
            if diff:
                probs.append(f"段 2 の prep の {diff} が段 1 の prep と違う (変えてよいのは Md_moc_offset と壁だけ)")
        except (OSError, ValueError) as e:
            probs.append(f"段 1 の {PREP_JSON} を読めない ({type(e).__name__})")
    rec["problems"] = probs
    rec["ok"] = not probs
    return rec


# --- 1 段の評価 -------------------------------------------------------------------------------------------------------------------
def series_csv(path: Path, rows: list, cols: list) -> None:
    EV2.write_series_csv(path, rows, cols)


def evaluate_stage(case: Path, stage: str, expected_delta, eta_c, eta_rec) -> dict:
    """1 段 (1 本の run) の前提・出口 M の系列。判定は呼び出し側 (段 1・段 2 で別)。"""
    rd = case / RUNS[stage]
    out = {"run": RUNS[stage], "problem": PROBLEMS[stage], "expected_delta": expected_delta, "snapshots": {}}
    pre = []
    out["fixed"] = fixed_conditions(case, stage, expected_delta)
    pre += [f"[固定] {p}" for p in out["fixed"]["problems"]]
    pre += [f"[固定] η の列: {p}" for p in eta_rec.get("problems", [])]
    pre += [f"[固定] 評価器の窓: {p}" for p in window_problems()]
    segf = rd / SEGMENT_VERDICT_FILE
    conv = EV2.parse_convergence(segf.read_text() if segf.is_file() else None)
    conv["label"] = EV2.convergence_label(conv)
    out["convergence"] = conv
    if not EV2.convergence_ok(conv):
        pre.append(f"[P1] {conv['label']} — {conv.get('reasons')}")
    eff = run_steps(rd, conv)
    sp = [p for p in step_problems(eff) if p not in window_problems()]
    out["steps"] = {"effective": eff, "problems": sp}
    pre += [f"[固定] {p}" for p in sp]
    # 格子
    try:
        X, R, S = EV2.load_geometry(rd)
        eta = EV2.column_eta(R)
        ep = EV2.eval_points(X, R, X, R)       # 評価点 = この run の節点 (E2 の作り方で A = B = 自格子)
    except Exception as e:  # noqa: BLE001
        out["geometry_error"] = f"{type(e).__name__}: {e}"
        pre.append(f"[固定] 格子を読めない: {out['geometry_error']}")
        out["precondition"] = {"ok": False, "reasons": pre}
        out["exit_series"] = {}
        return out
    ni, nj = X.shape
    out["final_section"] = {"x_rt": float(X[-1, 0]), "r_w_rt": float(R[-1, -1]), "n_own_in_band": int(np.sum((eta[-1] >= ETA_BAND[0]) & (eta[-1] <= ETA_BAND[1])))}
    exitM, ind, rows = {}, {}, []
    snaps = [("main", s, rd / f"res_{s}.h5") for s in MAIN_STEPS]
    for stg, step, h5 in snaps:
        rec = {"exit": exit_snapshot(h5, ni, nj, eta[-1], eta_c) if eta_c is not None else {"problems": ["η の列が無い"]}}
        t0 = EV2.eval_snapshot(rd, h5, X, R, eta, ep)
        rec["t0"] = {k: t0.get(k) for k in ("sha256", "problems", "thermo_source", "newton_ok", "raw", "h0_minus_hmix_Tt", "indicators")}
        out["snapshots"][f"{stg}:{step}"] = rec
        row = {"stage": stg, "step": step}
        if not rec["exit"]["problems"]:
            exitM[step] = rec["exit"]["exitM_common"]
            row.update({c: rec["exit"][c] for c in EXIT_COLS})
        elif step in WIN13:
            pre.append(f"[P3] 判定窓の res_{step}.h5 の出口 M を評価できない ({rec['exit']['problems']})")
        if t0.get("indicators") is not None and not t0["problems"]:
            ind[step] = t0["indicators"]
            row.update(EV2.series_row(stg, step, t0["indicators"]))
        elif step in WIN13:
            pre.append(f"[P2] 判定窓の res_{step}.h5 の全温を評価できない ({t0.get('problems')})")
        rows.append(row)
    cols_all = list(EXIT_COLS) + EV2.all_columns()
    series_csv(rd / SERIES_CSV, rows, cols_all)
    out["series_csv"] = {"file": f"{RUNS[stage]}/{SERIES_CSV}", "sha256": _sha(rd / SERIES_CSV)}
    out["exit_series"] = {s: exitM[s] for s in sorted(exitM)}
    out["exit_own_series"] = {s: out["snapshots"][f"main:{s}"]["exit"]["exitM_own"] for s in sorted(exitM)}
    # P3 出口 M (窓ごとに CSV を書いて check_quasisteady に渡す)
    qs = {}
    for wn, wsteps in WINDOWS.items():
        wcsv = rd / f"e4_recal_exitM_{wn}.csv"
        series_csv(wcsv, [r for r in rows if r["step"] in wsteps and "exitM_common" in r], list(EXIT_COLS))
        qs[wn] = EV2.run_quasisteady(wcsv, list(EXIT_COLS), rd / f"e4_recal_exitM_{wn}.txt")
    ex = exit_window_check(exitM, qs)
    ex["own_verdicts_recorded"] = {wn: EV2.parse_quasisteady(t, ["exitM_own"])["exitM_own"] for wn, t in qs.items()}
    out["exit_precondition"] = ex
    if not ex["ok"]:
        pre += [f"[P3] {p}" for p in ex["problems"]]
    # P2 全温の健全性
    tw = t0_within_check(ind)
    out["t0_within"] = tw
    if not tw["ok"]:
        pre += [f"[P2a] {p}" for p in tw["problems"]]
    cols = EV2.qs_columns()
    out["t0_range"] = {}
    for wn, wsteps in WINDOWS.items():
        wcsv = rd / f"e4_recal_t0_{wn}.csv"
        wrows = [r for r in rows if r["step"] in wsteps and all(c in r for c in cols)]
        series_csv(wcsv, wrows, cols)
        txt = EV2.run_quasisteady(wcsv, cols, rd / f"e4_recal_t0_{wn}.txt")
        steps, vals = EV2.read_csv_cols(wcsv, cols)
        wc = EV2.window_check(steps, vals, txt, cols, wsteps)
        out["t0_range"][wn] = wc
        if not wc["ok"]:
            pre.append(f"[P2b] {wn}: {wc['problems']} 不成立の列 {wc['n_bad']}/{len(cols)} "
                       f"(例 {[(c, wc['columns'][c]['range_K']) for c in wc['bad'][:3]]})")
    out["precondition"] = {"ok": not pre, "reasons": pre}
    return out


def evaluate(case: Path) -> dict:
    case = Path(case).resolve()
    out = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "evaluator": Path(__file__).name, "evaluator_sha256": _sha(Path(__file__)),
           "deps_sha256": {**{n: _sha(HERE / n) for n in ("euler_t0_e2_eval.py", "euler_t0_stage_ab.py", "moc_v5c_thermo_ab.py",
                                                         "throat_mono_judge.py", "ic_index_map.py")},
                           **{n: _sha(TOOLS / n) for n in ("check_quasisteady.py", "check_convergence.py", "forge_species.py")},
                           **{f"design/{n}": _sha(ROOT / "design" / n) for n in ("forge_design/evaluate/runner_axismach.py",
                                                                                "forge_design/meshing/mesh2d.py")}},
           "runs": RUNS, "problems": PROBLEMS, "delta0": DELTA0, "delta1_basis": DELTA1_BASIS,
           "definitions": {"exit": f"最終断面の η ∈ {list(ETA_BAND)}、NS の基準格子 G1 の固定の η の列への線形補間の単純平均 (M_common)。自格子の平均 M_own は別列",
                           "target": f"判定窓 {WIN13[0]}〜{WIN13[-1]} の 13 枚すべてで |M_common − {M_TARGET:g}| ≤ {M_TOL:g}",
                           "exit_precondition": f"両窓 STEADY、13 枚の幅 ≤ {EXIT_RANGE_TOL:g}、|末尾 5 枚 − 直前 5 枚 ({PREV5[0]}〜{PREV5[-1]}) の平均| ≤ {EXIT_MEANDIFF_TOL:g}",
                           "t0": f"判定窓 13 枚で全 9 領域の T0 − 1600 が ±{EV2.WITHIN_K:g} K 以内、27 列の時間の幅 ≤ {EV2.RANGE_TOL_K:g} K (両窓)",
                           "convergence": EV2.CONV_DEFINITION, "delta1": "δ₁ = δ₀ − (判定窓 13 枚の M_common の平均 − 6) (末尾 5 枚の平均による値も記録)",
                           "windows": {k: list(v) for k, v in WINDOWS.items()}, "prev5": list(PREV5)},
           "scope_note": SCOPE_NOTE}
    eta_c, eta_rec = eta_common(case)
    out["eta_common"] = eta_rec
    st = {}
    s1 = evaluate_stage(case, "d0", DELTA0, eta_c, eta_rec) if (case / RUNS["d0"]).is_dir() else None
    if s1 is None:
        st["d0"] = {"run": RUNS["d0"], "judgment": {"verdict": LBL_NOTRUN, "reasons": [f"{RUNS['d0']} が無い"], "facts": {}}}
    else:
        s1["judgment"] = judge_stage1(s1["exit_series"], s1["precondition"]["ok"], s1["precondition"]["reasons"])
        st["d0"] = s1
    j1 = st["d0"]["judgment"]
    if (case / RUNS["d1"]).is_dir():
        d1 = j1.get("delta1")
        s2 = evaluate_stage(case, "d1", d1, eta_c, eta_rec)
        info2 = {}
        try:
            info2 = json.loads((case / RUNS["d1"] / "prepare_info.json").read_text())
        except (OSError, ValueError):
            pass
        s2["judgment"] = judge_stage2(s2.get("exit_series", {}), s2["precondition"]["ok"], s2["precondition"]["reasons"], j1,
                                      info2.get("Md_moc_offset"))
        st["d1"] = s2
    else:
        st["d1"] = {"run": RUNS["d1"], "judgment": {"verdict": LBL_NOTRUN, "reasons": [f"{RUNS['d1']} が無い"
                                                                                       + (" (段 1 が「更新」なので段 2 が必要)" if j1["verdict"] == LBL_UPDATE else "")],
                                                    "facts": {}}}
    out["stages"] = st
    v1, v2 = st["d0"]["judgment"]["verdict"], st["d1"]["judgment"]["verdict"]
    if v1 == LBL_KEEP:
        final = {"status": "据え置き", "Md_moc_offset": DELTA0, "reference_run": RUNS["d0"]}
    elif v1 == LBL_UPDATE and v2 == LBL_PASS2:
        final = {"status": "更新して合格", "Md_moc_offset": st["d0"]["judgment"]["delta1"], "reference_run": RUNS["d1"]}
    elif v1 == LBL_UPDATE and v2 == LBL_NOTRUN:
        final = {"status": "段 2 待ち", "delta1": st["d0"]["judgment"]["delta1"], "delta1_repr": st["d0"]["judgment"]["delta1_repr"]}
    else:
        final = {"status": "保留", "stage1": v1, "stage2": v2}
    out["final"] = final
    out["VERDICT"] = {"d0": v1, "d1": v2}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("case", nargs="?", default=str(HERE))
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case = Path(a.case).resolve()
    out = evaluate(case)
    op = Path(a.out) if a.out else case / OUT_JSON
    op.parent.mkdir(parents=True, exist_ok=True)
    op.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    er = out["eta_common"]
    print(f"η の列: {er.get('n')} 点 ({(er.get('sha256') or '')[:16]}) 照合 {er.get('xcheck_ns_run')}" + (f" 問題 {er['problems']}" if er.get("problems") else ""))
    for stg in STAGES:
        s = out["stages"][stg]
        j = s["judgment"]
        print(f"[{stg}] {s.get('run')}: {j['verdict']}")
        if "convergence" in s:
            ex = s.get("exit_precondition", {})
            print(f"  収束 {s['convergence']['label']} / 出口 M の前提 {'成立' if ex.get('ok') else '不成立'} "
                  f"(VERDICT {ex.get('verdicts')}、13 枚の幅 {ex.get('range13')}、末尾 − 直前 {ex.get('meandiff_tail5_prev5')}) / "
                  f"全温 ±1 K {'成立' if s.get('t0_within', {}).get('ok') else '不成立'} (最大 {s.get('t0_within', {}).get('max_over_steps_K')})")
        f = j.get("facts", {})
        if "max_abs_dev" in f:
            print(f"  窓の M_common: 平均 {f.get('mean_win13')}、max|M_common − 6| {f['max_abs_dev']:.3e}" +
                  (f"、δ₁ = {f['delta1_repr']} (末尾 5 枚なら {f['delta1_if_tail5']!r})" if "delta1_repr" in f else ""))
        for r in j.get("reasons", [])[:12]:
            print("   -", r)
    print("FINAL:", json.dumps(out["final"], ensure_ascii=False))
    print("注:", SCOPE_NOTE)
    print("→", op)
    return 0


if __name__ == "__main__":
    sys.exit(main())
