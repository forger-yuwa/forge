"""plan discretization-moc-axis-limit-and-corrector §6 V5d (2026-10-07 登録 99431498、結果を見る前。諮問
notes/reviews/2026-10-07-euler-grid-switch-plan-diagnose.md) の評価器: 新しい Euler の格子 (mesh_euler) での MOC の比較。
起動は run_moc_v5d.sh (準備・実行は moc_v5d.py)。この評価器は保存場と時系列を読むだけで、forge は起動しない。

腕と run (登録・依頼の固定):
  腕 B = 単調壁・legacy + fixed2 (problem_d155_euler_v5d_B.yaml)、腕 M = 単調壁・analytic + converge (problem_d155_euler_v5d_M.yaml)。
  両腕とも同じ出口較正の値 (plan verification-case45-euler-total-enthalpy の E4/E4V の結果) を Md_moc_offset に書く (--md-offset)。
  格子は E3 の mesh_euler (2000 × 97・全域 0.005)。各腕 3 本、等エントロピー IC (腕ごとに 1 つの prep を 3 本に複製)。
  既定 (共用しない): 腕 B run_0171〜0173_euler_v5d_B_r{1,2,3}、腕 M run_0174〜0176_euler_v5d_M_r{1,2,3} (6 本を新しく回す)。
  --share-e4v: E4V の run (run_0164_euler_e4_recal_d1) を腕 B の 1 本として共用する。腕 B = run_0171・run_0172・run_0164
    (run_0173 は作らない)。共用の条件: E4 の評価の最終が採用で参照 run が run_0164、その Md_moc_offset が --md-offset と同じ、
    腕 B の新しい run と問題 (name を除く)・壁 (保存 spline)・格子 (座標と接続のハッシュ)・設定ファイル・forge のバイナリが同じ。
  旧 G1 の run は対照に混ぜない (全 run で格子の採用元が mesh_euler・登録の実効値であることを前提にする)。
長さ・窓 (E4 と同じ): soft 3000 (1 次・cfl 0.5) + 本段 54000 (2 次・cfl 2・implicitRelax 0.7)、出力 1000 ごと。
  判定窓は本段 42000〜54000 の 13 枚 (win13) と 50000〜54000 の 5 枚 (tail5)。窓は動かさない・延長しない。
量 (V5 と同じ): M 波・P 波・オーバーシュート・出口規格化オーバーシュート・|P 傾き| (いずれも r/r_w = 0.1) と |出口コア M − 6|
  (exit_M_dev)。Δq = 0.001 / 0.010 / 0.003 / 0.003 / 0.03 %pt、1.8e-4。出口コア M (exit_core_M、η 重みの面積平均) は記録用。
  評価量の時系列は eval_wallfit_euler.py --series (tag v5d) が書く。X_E・X_F は腕 B の r1 (run_0171) から取り全 run に同じ値。
統計・判定 (V5 と同じ定義、moc_v5_euler_eval.arm_stats・judge_one_sided): run i の窓平均 m_i、s_i = sd_i/√13、
  腕の SE = max(sd(m_i)/√k, √(Σ s_i²)/k)、D = M_M − M_B、SE_D = √(SE_B² + SE_M²)。
  全量で D + 2·SE_D ≤ Δq → 許容幅内 / D − 2·SE_D ≥ Δq の量があれば悪化 / それ以外は保留。SE は自己相関を補正しない実務の指標。
前提 (全 run に同じ条件。1 つでも不成立なら総合は保留。結果を見て窓を動かしたり延長したりしない):
  残差 (E4 と同じ): 本段区間の check_convergence が PASS、または全不合格列が停滞 (STALLED/plateau) だけ (euler_t0_e2_eval.parse_convergence)。
  全温 (E4 と同じ = E2 の健全性): 判定窓 13 枚すべてで全 9 領域の T0 − 1600 が ±1 K 以内、27 列の時間の幅 ≤ 0.1 K (両窓)。
  準定常: 6 量と出口コア M を、両窓 (13 枚・5 枚) の check_quasisteady --tail 1 (drift 0.05・osc 0.10) で判定し、両方 STEADY。
    例外は exit_M_dev だけ (V5 と同じ): 符号付き exit_core_M − 6 が窓の末尾 10 枚で E4 の「絶対許容内」(許容 1.8e-5) を満たし、
    かつ 13 枚の最大絶対値と幅がともに 1.8e-5 以下なら比較に使う (表示は「絶対許容内」。STEADY と書かない)。
  幅 (V5d で追加): 各量の 13 枚の幅 (最大 − 最小) ≤ Δq/10、符号付きの出口 M (exit_core_M − 6) の 13 枚の幅 ≤ 1e-5。
  固定の条件: forge の sha256 (本段・soft 段) が run_0143 の RUN_PROVENANCE と同じ、段 soft → main・cfl 0.5 / 2・relax 0.7、
    本段 54000・出力 1000・残差の最終 step 53999、IC は等エントロピーで起動前の検査が OK、Tt 1600、MOC が腕どおり (腕 M はゲート合格)、
    単調壁 [0, 1.5] (保存 spline の単調性も)、格子の採用元 mesh_euler・実効値、Md_moc_offset = --md-offset (全 run 同じ値)、
    凍結の初期線 run_0062、r_t、壁節点が自腕の spline に一致 (腕 M は腕 B と見分けられる節点がすべて自腕に一致)、
    腕の中の run の壁・格子が同じ、全 run の設定ファイル (prep 時の bcondConfig・solverConfig・probe・species_meta) と本段の実効設定が同じ。
    時系列の記録 (_band_ab/wallfit_series_v5d.json): tag v5d、基準 run_0171、共通の評価座標、各 run 自身の X_E・X_F との差 ≤ 1e-6 r_t
    (座標の整合の許容差であって評価誤差の保証ではない)、CSV・準定常ファイルの sha256、54 枚・最終 54000。
    出口の共通の η の列 (E4 の 12 点、e4_recal_eval.eta_common の sha256 照合)。
    E4 の結果 (_band_ab/e4_recal_eval.json): 評価器の sha256 が今の e4_recal_eval.py と同じ、最終が採用、Md_moc_offset = --md-offset。
    起動の記録 (_band_ab/moc_v5d_launch.json): 較正値・共用の有無・腕の構成が評価の引数と同じ。
出口較正と M = 6 の達成は分ける (登録):
  出口較正の据え置き (V5 の 3 区分、出口コア M の D): |D| + 2SE ≤ 1e-4 → 据え置き / |D| − 2SE > 1e-4 → やり直す / それ以外は保留。
    共通の η の列の M_common による同じ 3 区分は参考値として記録する。
  M = 6 の達成 (腕 M 自身、E4 の絶対目標): 腕 M の 3 本とも判定窓 13 枚すべてで |M_common − 6| ≤ 1e-4。前提は共通の前提と腕 M の
    3 本の前提に加え、E4 の P3 (M_common が両窓 STEADY、13 枚の幅 ≤ 5e-5、|末尾 5 枚 − 直前 5 枚の平均| ≤ 5e-5) を腕 M の各 run に課す。
    未達なら NS を始めず、新しい共通の較正値で両腕を確かめ直す (§6 V5d)。腕 B の同じ値は記録だけ。
結論は「指定の等エントロピーの初期化での比較」に限る (IC に依存しないとは言わない)。
記録だけ (判定に使わない): 評価座標の感度 (moc_v5_euler_eval.coordinate_record、窓の最後の枚 54000)、eval_wallfit_euler の記録の準定常
  ファイル (系列全体の末尾 5 枚) の VERDICT、M_own (自格子の帯内平均)。
出力: _band_ab/moc_v5d_eval.json (評価器の sha256・登録の commit・定義・各 run の前提と系列・判定)。各 run に moc_v5d_q_{win13,tail5}.csv/.txt
  (6 量 + 出口コア M の窓)、moc_v5d_t0_{win13,tail5}.csv/.txt (全温の主指標)、moc_v5d_exitM_{win13,tail5}.csv/.txt (M_common)。
usage: python3 moc_v5d_eval.py [case_dir] --md-offset <repr> [--share-e4v] [--out PATH]
"""
from __future__ import annotations

import argparse
import csv
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
import euler_t0_e2_eval as EV2  # noqa: E402  (残差の列ごとの内訳・全温の主指標・check_quasisteady の呼び方)
import e4_recal_eval as EV4  # noqa: E402  (共通の η の列・M_common・E4 の P3・段の manifest・step 数の同期)
import moc_v5_euler_eval as V5  # noqa: E402  (V5 の 6 量・Δq・腕の統計・壁の照合・評価座標の記録・実効設定)

PLAN = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5d"
PLAN_REG_COMMIT = "99431498"                      # V5d の登録 (2026-10-07、結果を見る前)
TAG = "v5d"
RUNS_NEW = {"B": ("run_0171_euler_v5d_B_r1", "run_0172_euler_v5d_B_r2", "run_0173_euler_v5d_B_r3"),
            "M": ("run_0174_euler_v5d_M_r1", "run_0175_euler_v5d_M_r2", "run_0176_euler_v5d_M_r3")}
E4V_RUN = EV4.RUNS["d1"]                          # run_0164_euler_e4_recal_d1 (共用するときの腕 B の 1 本)
GEOM_REF = RUNS_NEW["B"][0]                       # X_E・X_F と実効設定の基準 (腕 B の r1)
PROBLEMS = {"B": "problem_d155_euler_v5d_B.yaml", "M": "problem_d155_euler_v5d_M.yaml"}
PROBLEM_NAMES = {"B": "isobutane_m6_d155_euler_v5d_B", "M": "isobutane_m6_d155_euler_v5d_M"}
PREPS = {"B": "_prep_v5d_B", "M": "_prep_v5d_M"}
PREP_JSON = "V5D_PREP.json"
LAUNCH_JSON = "_band_ab/moc_v5d_launch.json"
SERIES_JSON = f"_band_ab/wallfit_series_{TAG}.json"
E4_EVAL_JSON = EV4.OUT_JSON                       # _band_ab/e4_recal_eval.json
OUT_JSON = "_band_ab/moc_v5d_eval.json"
SERIES_CSV = f"wallfit_series_{TAG}.csv"
SERIES_QS = f"QUASISTEADY_wallfit_{TAG}.txt"

# 長さ・窓 (E4 と同じ。window_problems で E4 の定数とも照合する)
SOFT_STEPS = 3000
MAIN_NSTEPS = 54000
OUT_INTERVAL = 1000
MAIN_STEPS = tuple(range(0, MAIN_NSTEPS + 1, OUT_INTERVAL))
WIN13 = tuple(range(MAIN_NSTEPS - 12 * OUT_INTERVAL, MAIN_NSTEPS + 1, OUT_INTERVAL))   # 42000〜54000
TAIL5 = WIN13[-5:]                                                                    # 50000〜54000
WINDOWS = {"win13": WIN13, "tail5": TAIL5}
CFL_SOFT, CFL_MAIN, RELAX = 0.5, 2.0, 0.7
SERIES_N_SNAPS = MAIN_NSTEPS // OUT_INTERVAL      # eval_wallfit_euler の系列は本段の step > 0 (54 枚)

# 量と許容幅 (V5 と同じ)
DQ = dict(V5.DQ)
RECORD_ONLY = ("exit_core_M",)
QS_COLS = list(DQ) + list(RECORD_ONLY)
WIDTH_FRAC = 0.1                                  # 各量の 13 枚の幅 ≤ Δq/10 (V5d で追加)
EXIT_SIGNED_WIDTH_TOL = 1e-5                      # 符号付きの出口 M (exit_core_M − 6) の 13 枚の幅
EXIT_EXC_QTY = "exit_M_dev"                       # 絶対許容内の例外を認める唯一の量 (V5 と同じ)
EXIT_EXC_FULL_TOL = 1.8e-5                        # 例外の追加条件: 13 枚の最大絶対値と幅
EXIT_CAL_TOL = 1e-4
M_TARGET = 6.0
M6_TOL = EV4.M_TOL                                # E4 の絶対目標 |M_common − 6| ≤ 1e-4 (13 枚すべて)
GEOM_TOL = V5.GEOM_TOL                            # 1e-6 r_t (座標の整合の許容差)
MOC = {"B": {"axis_limit": "legacy", "corrector": "fixed2"}, "M": {"axis_limit": "analytic", "corrector": "converge"}}
SAME_FILES = ("bcondConfig.yaml", "solverConfig.yaml", "probe.yaml", "species_meta.yaml")   # prep 時に全 run で同じであるべきファイル
# E4 の評価 (e4_recal_eval.evaluate の final.status) のうち、較正値が採用されたもの
E4_ADOPTED = ("据え置き", "更新して合格", "E4V 合格 (独立の検証。段 1 は判別不能のまま)")
MD_ABS_MAX = 0.01                                 # 出口較正の補正量としての範囲 (桁の取り違えを拒む; NS の生成器と同じ)

LBL_WIN, LBL_WORSE, LBL_HOLD = V5.WIN_LABEL, V5.WORSE_LABEL, V5.HOLD_LABEL   # 許容幅内 / 悪化 / 保留
LBL_ABS = "絶対許容内"
CAL_NOTE = "「据え置き」は両腕に共通の較正値を変える必要がないという意味で、M = 6 の達成の証明ではない (達成は armM_M6 で別に判定)"
LBL_CAL_KEEP = f"据え置き (Md_moc_offset を変えない。{CAL_NOTE})"
LBL_CAL_REDO = "較正をやり直す (|D| − 2SE > 1e-4)"
LBL_CAL_HOLD = "保留 (判別不能; 較正は変えない)"
LBL_M6_OK = "達成: 腕 M の 3 本とも判定窓 13 枚すべてで |M_common − 6| ≤ 1e-4 (E4 の絶対目標)"
LBL_M6_NG = ("未達: 腕 M は E4 の絶対目標を外れる → 再較正が要る (NS を始めない。新しい共通の較正値で両腕を確かめ直す; "
             "共通の値で両者を通せなければ MOC 単独の比較と各設計の最適な較正を同じ試験として扱わない)")
LBL_M6_UNDET = "判別不能 (前提の未達)"
SCOPE_NOTE = ("指定の等エントロピーの初期化 (各腕 1 つの prep を 3 本に複製) での比較であり、IC に依存しないとは言わない。"
              "停滞のみの NOT CONVERGED は収束の証明ではない。SE は自己相関を補正しない実務の指標で、保証された信頼限界・厳密な非劣化を主張しない")


def _sha(p) -> str | None:
    return EV2._sha(p)


def _json(p: Path) -> dict:
    return json.loads(Path(p).read_text())


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))


def arms(share: bool) -> dict:
    """腕の構成 (登録・依頼の固定)。share なら腕 B = run_0171・run_0172・run_0164 (run_0173 は作らない)。"""
    return {"B": (list(RUNS_NEW["B"][:2]) + [E4V_RUN]) if share else list(RUNS_NEW["B"]), "M": list(RUNS_NEW["M"])}


def new_runs(share: bool) -> dict:
    """このスクリプトが新しく作る run (腕ごと)。"""
    a = arms(share)
    return {arm: [r for r in runs if r != E4V_RUN] for arm, runs in a.items()}


def parse_md(text) -> tuple:
    """--md-offset の文字列 → (float, repr)。有限・|値| ≤ 0.01 でなければ ValueError。repr を正本の書き方にする。"""
    try:
        v = float(str(text).strip())
    except (TypeError, ValueError):
        raise ValueError(f"--md-offset {text!r} は数でない") from None
    if not math.isfinite(v):
        raise ValueError(f"--md-offset {text!r} は有限でない")
    if abs(v) > MD_ABS_MAX:
        raise ValueError(f"--md-offset {v!r} は |値| ≤ {MD_ABS_MAX:g} の外 (出口較正の補正量として大きすぎる。単位・桁を確かめる)")
    return v, repr(v)


# --- 長さ・窓の照合 -------------------------------------------------------------------------------------------------------------
def window_problems() -> list:
    """V5d の窓の定数が本段の定数から作られた形で、E4 の定数 (登録「E4 と同じ」) と同じか。"""
    bad = []
    if len(WIN13) != 13 or WIN13[-1] != MAIN_NSTEPS or any(b - a != OUT_INTERVAL for a, b in zip(WIN13, WIN13[1:])):
        bad.append(f"判定窓 {WIN13[0]}〜{WIN13[-1]} が本段の最後の 13 枚でない")
    if len(TAIL5) != 5 or TAIL5 != WIN13[-5:]:
        bad.append("末尾の窓が判定窓の最後の 5 枚でない")
    e4 = {"SOFT_STEPS": EV4.SOFT_STEPS, "MAIN_NSTEPS": EV4.MAIN_NSTEPS, "OUT_INTERVAL": EV4.OUT_INTERVAL, "WIN13": EV4.WIN13,
          "TAIL5": EV4.TAIL5, "CFL": (EV4.CFL_SOFT, EV4.CFL_MAIN, EV4.RELAX)}
    mine = {"SOFT_STEPS": SOFT_STEPS, "MAIN_NSTEPS": MAIN_NSTEPS, "OUT_INTERVAL": OUT_INTERVAL, "WIN13": WIN13, "TAIL5": TAIL5,
            "CFL": (CFL_SOFT, CFL_MAIN, RELAX)}
    bad += [f"{k} が E4 と違う ({mine[k]!r} / E4 {e4[k]!r})" for k in mine if mine[k] != e4[k]]
    bad += [f"E4 の評価器の窓: {p}" for p in EV4.window_problems()]
    return bad


def step_problems(eff: dict) -> list:
    """実効値 {nStepOuter, outStepInterval, manifest_main_nStepOuter, residual_last_step} と V5d の定数の食い違い。"""
    bad = list(window_problems())
    if eff.get("nStepOuter") != MAIN_NSTEPS:
        bad.append(f"本段の nStepOuter {eff.get('nStepOuter')!r} が {MAIN_NSTEPS} と違う")
    if eff.get("outStepInterval") != OUT_INTERVAL:
        bad.append(f"本段の outStepInterval {eff.get('outStepInterval')!r} が {OUT_INTERVAL} と違う")
    if "manifest_main_nStepOuter" in eff and eff["manifest_main_nStepOuter"] != MAIN_NSTEPS:
        bad.append(f"stage_manifest の main の nStepOuter {eff['manifest_main_nStepOuter']!r} が {MAIN_NSTEPS} と違う")
    if "residual_last_step" in eff and eff["residual_last_step"] != MAIN_NSTEPS - 1:
        bad.append(f"本段の残差の最終 step {eff['residual_last_step']!r} が {MAIN_NSTEPS - 1} でない")
    return bad


# --- E4 の結果・起動の記録 ------------------------------------------------------------------------------------------------------
def e4_adoption(case: Path, md: float, share: bool) -> dict:
    """E4 の評価 (_band_ab/e4_recal_eval.json) が較正値を採用していて、その値が md と同じか。share なら参照 run が run_0164。"""
    p = Path(case) / E4_EVAL_JSON
    rec = {"file": E4_EVAL_JSON, "problems": []}
    try:
        ev = _json(p)
    except (OSError, ValueError) as e:
        rec["problems"].append(f"{E4_EVAL_JSON} を読めない ({type(e).__name__})")
        return rec
    fin = ev.get("final") or {}
    rec.update(evaluator_sha256=ev.get("evaluator_sha256"), final=fin)
    if ev.get("evaluator_sha256") != _sha(HERE / "e4_recal_eval.py"):
        rec["problems"].append(f"E4 の評価器の sha256 ({str(ev.get('evaluator_sha256'))[:16]}) が今の e4_recal_eval.py と違う")
    if fin.get("status") not in E4_ADOPTED:
        rec["problems"].append(f"E4 の最終 {fin.get('status')!r} が採用でない ({' / '.join(E4_ADOPTED)} のどれでもない)")
    v = fin.get("Md_moc_offset")
    if not (_is_num(v) and float(v) == md):
        rec["problems"].append(f"E4 の採用値 {v!r} が --md-offset {md!r} と違う")
    if share and fin.get("reference_run") != E4V_RUN:
        rec["problems"].append(f"共用: E4 の参照 run {fin.get('reference_run')!r} が {E4V_RUN} でない")
    rec["ok"] = not rec["problems"]
    return rec


def launch_problems(case: Path, md_repr: str, share: bool) -> tuple:
    """起動の記録 (投入スクリプトが forge の前に書く) が評価の引数と同じ構成か。戻り (記録, 不成立の理由)。"""
    try:
        lr = _json(Path(case) / LAUNCH_JSON)
    except (OSError, ValueError) as e:
        return {}, [f"{LAUNCH_JSON} を読めない ({type(e).__name__})"]
    bad = []
    if lr.get("dry"):
        bad.append("起動の記録が乾式確認 (DRY) のもの")
    if lr.get("md_offset_repr") != md_repr:
        bad.append(f"起動の記録の較正値 {lr.get('md_offset_repr')!r} が {md_repr!r} と違う")
    if lr.get("share_e4v") is not share:
        bad.append(f"起動の記録の共用 {lr.get('share_e4v')!r} が {share} と違う")
    if lr.get("arms") != arms(share):
        bad.append(f"起動の記録の腕の構成 {lr.get('arms')} が {arms(share)} と違う")
    return lr, bad


# --- 評価量の時系列 -------------------------------------------------------------------------------------------------------------
def load_window(case: Path, run: str) -> dict:
    """wallfit_series_v5d.csv の判定窓 13 枚 → {量: ndarray}。欠けた枚・非有限・列の欠損は ValueError (0 で埋めない)。"""
    f = Path(case) / run / SERIES_CSV
    if not f.is_file():
        raise FileNotFoundError(f"{SERIES_CSV} が無い (eval_wallfit_euler.py --series を先に回すこと)")
    with open(f) as fh:
        rows = list(csv.DictReader(fh))
    win = [r for r in rows if WIN13[0] <= int(float(r["step"])) <= WIN13[-1]]
    steps = [int(float(r["step"])) for r in win]
    if steps != list(WIN13):
        raise ValueError(f"判定窓 {WIN13[0]}〜{WIN13[-1]} の 13 枚がそろっていない (実際 {steps})")
    out = {}
    for k in QS_COLS:
        if k not in win[0]:
            raise ValueError(f"列 {k} が無い")
        v = np.array([float(r[k]) if r[k] not in ("", None) else np.nan for r in win])
        if not np.all(np.isfinite(v)):
            raise ValueError(f"{k} に非有限値・空欄")
        out[k] = v
    return out


def series_record_problems(srec: dict, runs: list) -> tuple:
    """時系列の記録 (_band_ab/wallfit_series_v5d.json) の共通の検査と run ごとの検査。戻り ({"共通": [...], run: [...]}, 記録)。"""
    common, per = [], {}
    if srec.get("tag") != TAG:
        common.append(f"時系列の記録の tag {srec.get('tag')!r} が {TAG!r} でない")
    if srec.get("geom_ref") != GEOM_REF:
        common.append(f"時系列の記録の X_E・X_F の基準 {srec.get('geom_ref')!r} が {GEOM_REF} (腕 B の r1) でない")
    common += V5.geometry_failures(srec.get("geometry"))
    gr = ((srec.get("runs") or {}).get(GEOM_REF)) or {}
    if not (gr.get("dX_E") == 0.0 and gr.get("dX_F") == 0.0):
        common.append(f"評価座標が基準 {GEOM_REF} 自身の座標でない・記録が無い (差 {gr.get('dX_E')!r}, {gr.get('dX_F')!r})")
    rec = {}
    for run in runs:
        why = []
        r = (srec.get("runs") or {}).get(run)
        if r is None:
            per[run] = [f"時系列の記録 ({SERIES_JSON}) に無い"]
            continue
        if r.get("status") != "ok":
            why.append(f"時系列の評価が {r.get('status')!r} ({r.get('reason')})")
        rec[run] = {k: r.get(k) for k in ("dX_E", "dX_F", "n_snaps", "last_step", "status", "csv_sha256", "quasisteady_sha256")}
        for k in ("dX_E", "dX_F"):
            v = r.get(k)
            if not (_is_num(v) and abs(v) <= GEOM_TOL):
                why.append(f"自身の {k[1:]} と評価の基準の差 {v!r} が座標の整合の許容差 {GEOM_TOL:g} r_t を超える・欠損")
        if r.get("last_step") != MAIN_NSTEPS or r.get("n_snaps") != SERIES_N_SNAPS:
            why.append(f"時系列が本段 {SERIES_N_SNAPS} 枚・最終 {MAIN_NSTEPS} でない ({r.get('n_snaps')!r} 枚・最終 {r.get('last_step')!r})")
        per[run] = why
    if common:
        per["共通"] = common
    return per, rec


def series_file_problems(case: Path, run: str, r: dict) -> list:
    """CSV・準定常ファイルの sha256 が時系列の記録と同じか (古い成果物・別の呼び出しを取り違えない)。"""
    why = []
    for key, name in (("csv_sha256", SERIES_CSV), ("quasisteady_sha256", SERIES_QS)):
        got = _sha(Path(case) / run / name)
        if not r.get(key) or got != r.get(key):
            why.append(f"{name} が時系列の記録と違う (sha256 記録 {str(r.get(key))[:12]} / 実際 {str(got)[:12]})")
    return why


# --- 準定常・幅・exit_M_dev の例外 -----------------------------------------------------------------------------------------------
def window_csvs_quasisteady(rd: Path, vals: dict) -> dict:
    """両窓の CSV (6 量と出口コア M) を書いて check_quasisteady --tail 1 を回す。戻り {窓: 出力}。"""
    out = {}
    for wn, ws in WINDOWS.items():
        idx = [WIN13.index(s) for s in ws]
        rows = [{"stage": "main", "step": WIN13[i], **{c: float(vals[c][i]) for c in QS_COLS}} for i in idx]
        p = rd / f"moc_v5d_q_{wn}.csv"
        EV2.write_series_csv(p, rows, QS_COLS)
        out[wn] = EV2.run_quasisteady(p, QS_COLS, rd / f"moc_v5d_q_{wn}.txt")
    return out


def exit_exception(signed) -> dict:
    """exit_M_dev だけの例外 (V5 と同じ条件、窓は V5d の 13 枚): 符号付き exit_core_M − 6 が、末尾 10 枚で E4 の「絶対許容内」
    (throat_mono_judge.abs_tolerance_verdict、許容 Δq/10 = 1.8e-5) を満たし、かつ 13 枚の最大絶対値と幅がともに 1.8e-5 以下。"""
    from throat_mono_judge import ABS_TOL_N, ABS_TOL_WITHIN, DELTA_Q, abs_tolerance_verdict
    v = np.asarray(signed, dtype=float)
    rec = {"tail_n": ABS_TOL_N, "tail_tol": DELTA_Q["exit_M_dev"] / 10.0, "full_tol": EXIT_EXC_FULL_TOL, "applies": False}
    if v.shape != (len(WIN13),) or not np.all(np.isfinite(v)):
        rec["reason"] = f"窓の 13 枚がそろわない・非有限 (形 {v.shape})"
        return rec
    e4 = abs_tolerance_verdict(list(WIN13), [float(x) for x in v], DELTA_Q["exit_M_dev"] / 10.0, n=ABS_TOL_N)
    maxabs, width = float(np.abs(v).max()), float(np.ptp(v))
    why = []
    if e4.get("status") != ABS_TOL_WITHIN:
        why.append(f"末尾 {ABS_TOL_N} 枚の E4 が {e4.get('status')} ({e4.get('reason')})")
    if not maxabs <= EXIT_EXC_FULL_TOL:
        why.append(f"13 枚の最大絶対値 {maxabs:.3g} > {EXIT_EXC_FULL_TOL:g}")
    if not width <= EXIT_EXC_FULL_TOL:
        why.append(f"13 枚の幅 {width:.3g} > {EXIT_EXC_FULL_TOL:g}")
    rec.update(e4_tail=e4, full_max_abs=maxabs, full_width=width, applies=not why, reason=("; ".join(why) or None),
               label=(LBL_ABS if not why else None))
    return rec


def judge_quasisteady(verdicts: dict, exc: dict) -> tuple:
    """verdicts = {窓: {量: VERDICT}}。両窓 STEADY なら STEADY、exit_M_dev は例外が成り立てば「絶対許容内」、それ以外は保留。
    戻り (不成立の理由 [文], {量: STEADY / 絶対許容内 / 保留})。"""
    status, bad = {}, []
    for c in QS_COLS:
        src = {wn: (verdicts.get(wn) or {}).get(c, "UNKNOWN") for wn in WINDOWS}
        if all(v == "STEADY" for v in src.values()):
            status[c] = "STEADY"
        elif c == EXIT_EXC_QTY and exc.get("applies") is True:
            status[c] = LBL_ABS
        else:
            status[c] = LBL_HOLD
            bad.append(", ".join(f"{c} {v} ({wn})" for wn, v in src.items() if v != "STEADY")
                       + (f" [絶対許容内の例外に当たらない: {exc.get('reason')}]" if c == EXIT_EXC_QTY else ""))
    return (["準定常でない量 (両窓 STEADY でなく、exit_M_dev の絶対許容内の例外にも当たらない): " + "; ".join(bad)] if bad else []), status


def width_checks(vals: dict) -> tuple:
    """各量の 13 枚の幅 ≤ Δq/10、符号付きの出口 M (exit_core_M − 6) の幅 ≤ 1e-5。戻り (不成立の理由, 記録)。"""
    rec, bad = {}, []
    for c, dq in DQ.items():
        w = float(np.ptp(vals[c]))
        ok = w <= dq * WIDTH_FRAC
        rec[c] = {"width": w, "tol": dq * WIDTH_FRAC, "ok": bool(ok)}
        if not ok:
            bad.append(f"{c} の 13 枚の幅 {w:.3g} > Δq/10 = {dq * WIDTH_FRAC:g}")
    w = float(np.ptp(np.asarray(vals["exit_core_M"], dtype=float) - M_TARGET))
    ok = w <= EXIT_SIGNED_WIDTH_TOL
    rec["exit_core_M_signed"] = {"width": w, "tol": EXIT_SIGNED_WIDTH_TOL, "ok": bool(ok)}
    if not ok:
        bad.append(f"符号付きの出口 M の 13 枚の幅 {w:.3g} > {EXIT_SIGNED_WIDTH_TOL:g}")
    return bad, rec


def window_conditions(case: Path, run: str, vals: dict) -> tuple:
    """1 run の窓の前提 (準定常・幅・exit_M_dev の例外)。戻り (不成立の理由, 記録)。"""
    rd = Path(case) / run
    qs = window_csvs_quasisteady(rd, vals)
    verd, rows_bad = {}, []
    for wn, txt in qs.items():
        verd[wn] = EV2.parse_quasisteady(txt, QS_COLS)
        m = re.search(r"\[(\d+) rows, steps", txt or "")
        if not m or int(m.group(1)) != len(WINDOWS[wn]):
            rows_bad.append(f"{wn}: check_quasisteady の行数 {m.group(1) if m else '読めない'} が {len(WINDOWS[wn])} でない")
    exc = exit_exception(np.asarray(vals["exit_core_M"], dtype=float) - M_TARGET)
    why_q, status = judge_quasisteady(verd, exc)
    why_w, wrec = width_checks(vals)
    qp = rd / SERIES_QS
    rec = {"status": status, "verdicts": verd, "exit_M_dev_exception": exc, "widths": wrec,
           "record_file_tail5_record_only": EV2.parse_quasisteady(qp.read_text() if qp.is_file() else "", QS_COLS),
           "exit_M_dev_vs_signed_maxabs": float(np.abs(np.asarray(vals["exit_M_dev"]) - np.abs(np.asarray(vals["exit_core_M"]) - M_TARGET)).max())}
    return rows_bad + why_q + why_w, rec


# --- 全温 (P2) と出口の M_common ------------------------------------------------------------------------------------------------
def t0_and_exit(case: Path, run: str, eta_c) -> dict:
    """判定窓 13 枚の全温の主指標 (E2 の復元・評価点 = この run の節点) と M_common (E4 の定義)。保存場を読む
    (試験では EV2.load_geometry・EV2.eval_snapshot・EV4.exit_snapshot を差し替える)。"""
    rd = Path(case) / run
    out = {"t0": {}, "exit": {}, "problems": []}
    try:
        X, R, _ = EV2.load_geometry(rd)
        eta = EV2.column_eta(R)
        ep = EV2.eval_points(X, R, X, R)
    except Exception as e:  # noqa: BLE001 — 格子が読めないことは前提の不成立
        out["problems"].append(f"格子を読めない: {type(e).__name__}: {e}")
        return out
    ni, nj = X.shape
    for s in WIN13:
        h5 = rd / f"res_{s}.h5"
        t0 = EV2.eval_snapshot(rd, h5, X, R, eta, ep)
        out["t0"][s] = {k: t0.get(k) for k in ("sha256", "problems", "thermo_source", "newton_ok", "indicators", "raw")}
        ex = EV4.exit_snapshot(h5, ni, nj, eta[-1], eta_c) if eta_c is not None else {"problems": ["共通の η の列が無い"]}
        out["exit"][s] = ex
    return out


def t0_conditions(rd: Path, te: dict) -> tuple:
    """P2: 13 枚すべて ±1 K (EV4.t0_within_check) と、27 列の時間の幅 ≤ 0.1 K を両窓で (EV2.window_check)。"""
    why, rec = [], {}
    ind = {s: v["indicators"] for s, v in te["t0"].items() if v.get("indicators") is not None and not v.get("problems")}
    bad_steps = [s for s in WIN13 if s not in ind]
    if bad_steps:
        why.append(f"[P2] 判定窓の全温を評価できない時点 {bad_steps} ({[te['t0'].get(s, {}).get('problems') for s in bad_steps][:2]})")
    tw = EV4.t0_within_check(ind)
    rec["within"] = tw
    if not tw["ok"]:
        why += [f"[P2a] {p}" for p in tw["problems"]]
    cols = EV2.qs_columns()
    rec["range"] = {}
    for wn, ws in WINDOWS.items():
        p = rd / f"moc_v5d_t0_{wn}.csv"
        rows = [EV2.series_row("main", s, ind[s]) for s in ws if s in ind]
        EV2.write_series_csv(p, rows, cols)
        txt = EV2.run_quasisteady(p, cols, rd / f"moc_v5d_t0_{wn}.txt")
        steps, vals = EV2.read_csv_cols(p, cols)
        wc = EV2.window_check(steps, vals, txt, cols, ws)
        rec["range"][wn] = {k: wc[k] for k in ("ok", "problems", "n_bad", "bad", "verdict_counts_recorded")} | {
            "max_range_K": max((v["range_K"] for v in wc["columns"].values() if v["range_K"] is not None), default=None)}
        if not wc["ok"]:
            why.append(f"[P2b] {wn}: {wc['problems']} 不成立の列 {wc['n_bad']}/{len(cols)} "
                       f"(例 {[(c, wc['columns'][c]['range_K']) for c in wc['bad'][:3]]})")
    return why, rec


def exit_common_conditions(rd: Path, te: dict) -> tuple:
    """M_common の系列 (判定窓 13 枚) と E4 の P3 (両窓 STEADY・13 枚の幅 ≤ 5e-5・|末尾 5 − 直前 5| ≤ 5e-5)。
    戻り (系列 {step: M_common}, P3 の記録, M_own の系列)。"""
    series, own, rows = {}, {}, []
    for s in WIN13:
        ex = te["exit"].get(s) or {}
        if not ex.get("problems") and _is_num(ex.get("exitM_common")):
            series[s] = float(ex["exitM_common"])
            own[s] = ex.get("exitM_own")
            rows.append({"stage": "main", "step": s, "exitM_common": series[s], "exitM_own": ex.get("exitM_own")})
    qs = {}
    for wn, ws in WINDOWS.items():
        p = rd / f"moc_v5d_exitM_{wn}.csv"
        EV2.write_series_csv(p, [r for r in rows if r["step"] in ws], list(EV4.EXIT_COLS))
        qs[wn] = EV2.run_quasisteady(p, list(EV4.EXIT_COLS), rd / f"moc_v5d_exitM_{wn}.txt")
    p3 = EV4.exit_window_check(series, qs)
    return series, p3, own


def m6_target(series: dict) -> dict:
    """E4 の絶対目標: 判定窓 13 枚すべてで |M_common − 6| ≤ 1e-4 (float の ≤ そのまま)。"""
    if any(s not in series or not math.isfinite(series[s]) for s in WIN13):
        return {"complete": False, "within": False}
    dev = {s: abs(series[s] - M_TARGET) for s in WIN13}
    return {"complete": True, "within": all(v <= M6_TOL for v in dev.values()), "max_abs_dev": max(dev.values()),
            "mean_win13": float(np.mean([series[s] for s in WIN13])), "abs_dev_by_step": dev}


# --- 固定の条件 (run ごと) ------------------------------------------------------------------------------------------------------
def prepare_info_problems(info: dict, arm: str, md: float) -> list:
    """prepare_info.json の固定の条件: MOC が腕どおり (腕 M はゲート合格)・単調壁 [0, 1.5]・mesh_euler の実効値・Md_moc_offset・
    凍結の初期線・r_t・IC の作り方 (等エントロピー)。乾式確認の印が無いこと。moc_v5d.py の prep・run も同じ関数で検査する。"""
    from throat_mono_judge import mono_r2_matches
    bad = []
    if info.get("DRY") or info.get("DRY_NO_IC"):
        bad.append("乾式確認 (DRY) の準備から作った run")
    moc = info.get("moc")
    if not isinstance(moc, dict):
        bad.append("prepare_info に moc が無い (MOC の実装より前の準備・古い design パッケージ)")
    else:
        got = {k: moc.get(k) for k in MOC[arm]}
        if got != MOC[arm]:
            bad.append(f"腕 {arm} の moc が {MOC[arm]} でない ({got})")
        gate = moc.get("gate") or {}
        if arm == "M" and not (gate.get("applicable") is True and gate.get("pass") is True):
            bad.append(f"腕 M の MOC のゲートが合格でない (applicable {gate.get('applicable')!r}, pass {gate.get('pass')!r})")
        if arm == "B" and gate.get("applicable") not in (False, None):
            bad.append(f"腕 B の MOC のゲートが legacy + fixed2 の「適用外」でない ({gate})")
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf or not mono_r2_matches(wf.get("mono_r2"), V5.MONO_R2):
        bad.append(f"wall_fit.mono_r2 {wf.get('mono_r2', '(記録なし)')!r} が {V5.MONO_R2} と一致しない")
    mesh = info.get("mesh") or {}
    if mesh.get("source") != "mesh_euler":
        bad.append(f"格子の採用元が {mesh.get('source')!r} (mesh_euler でない。旧 G1 の run を混ぜない)")
    got = {k: (mesh.get("params") or {}).get(k) for k in EV4.MESH_EULER_EXPECT}
    if got != EV4.MESH_EULER_EXPECT:
        bad.append(f"格子の実効値 {got} が登録 {EV4.MESH_EULER_EXPECT} と違う")
    v = info.get("Md_moc_offset")
    if not (_is_num(v) and float(v) == md):
        bad.append(f"Md_moc_offset {v!r} が --md-offset {md!r} でない")
    il = info.get("initial_line") or {}
    res = il.get("res")
    res = res if isinstance(res, list) else [res]
    if Path(str(il.get("run") or "")).name != EV4.INITIAL_LINE[0] or res != [EV4.INITIAL_LINE[1]]:
        bad.append(f"凍結の初期線 {il.get('run')!r}/{il.get('res')!r} が {EV4.INITIAL_LINE} でない")
    if not (_is_num(info.get("scale_m")) and float(info["scale_m"]) == EV4.SCALE_M):
        bad.append(f"r_t (scale_m) {info.get('scale_m')!r} が {EV4.SCALE_M} でない")
    ic = info.get("ic") or {}
    if ic.get("mode") != "isentropic" or ic.get("VERDICT") != "OK":
        bad.append(f"IC が等エントロピー・OK でない (mode {ic.get('mode')!r}, VERDICT {ic.get('VERDICT')!r})")
    return bad


def prep_record_problems(rd: Path, arm: str, md_repr: str, share: bool) -> tuple:
    """起動前の準備の記録: 新しい run は V5D_PREP.json (腕・run 名・較正値・IC の検査 OK)、共用の run_0164 は E4_PREP.json (段 d1・IC OK)。
    戻り (不成立の理由, 記録, prep 時のファイルの sha256)。"""
    why = []
    if rd.name == E4V_RUN:
        try:
            pr = _json(rd / EV4.PREP_JSON)
        except (OSError, ValueError) as e:
            return [f"{EV4.PREP_JSON} を読めない ({type(e).__name__})"], {}, {}
        if not share:
            why.append(f"{E4V_RUN} は共用しない構成の腕に入らない")
        if pr.get("stage") != "d1" or pr.get("run") != E4V_RUN or pr.get("dry"):
            why.append(f"{EV4.PREP_JSON} の段・run・乾式 ({pr.get('stage')!r}, {pr.get('run')!r}, {pr.get('dry')!r}) が d1・{E4V_RUN}・本番でない")
        if (pr.get("ic_check") or {}).get("VERDICT") != "OK":
            why.append(f"起動前の IC の検査 ({EV4.PREP_JSON}) が OK でない")
        return why, {"file": EV4.PREP_JSON, "stage": pr.get("stage"), "Md_moc_offset": pr.get("Md_moc_offset")}, pr.get("files_sha256") or {}
    try:
        pr = _json(rd / PREP_JSON)
    except (OSError, ValueError) as e:
        return [f"{PREP_JSON} を読めない ({type(e).__name__})"], {}, {}
    if pr.get("arm") != arm or rd.name not in (pr.get("runs") or []) or rd.name not in new_runs(share)[arm]:
        why.append(f"{PREP_JSON} の腕・run ({pr.get('arm')!r}, {pr.get('runs')}) が {arm}・{rd.name} と合わない")
    if pr.get("md_offset_repr") != md_repr:
        why.append(f"{PREP_JSON} の較正値 {pr.get('md_offset_repr')!r} が {md_repr!r} と違う")
    if pr.get("share_e4v") is not share:
        why.append(f"{PREP_JSON} の共用 {pr.get('share_e4v')!r} が {share} と違う")
    if pr.get("dry"):
        why.append(f"{PREP_JSON} が乾式確認のもの")
    if (pr.get("ic_check") or {}).get("VERDICT") != "OK":
        why.append(f"起動前の IC の検査 ({PREP_JSON}) が OK でない ({(pr.get('ic_check') or {}).get('VERDICT')!r})")
    if share and arm == "B" and (pr.get("same_as_e4v") or {}).get("status") != "OK":
        why.append(f"共用: 腕 B の prep と {E4V_RUN} の照合が OK でない ({(pr.get('same_as_e4v') or {}).get('status')!r})")
    return why, {"file": PREP_JSON, "arm": pr.get("arm"), "md_offset_repr": pr.get("md_offset_repr")}, pr.get("files_sha256") or {}


def fixed_conditions(case: Path, run: str, arm: str, md: float, md_repr: str, share: bool) -> tuple:
    """1 run の固定の条件。戻り (記録 {..., "problems", "ok"}, prepare_info)。"""
    rd = Path(case) / run
    probs, rec = [], {}
    rc = (rd / "RUN_RC").read_text().strip() if (rd / "RUN_RC").is_file() else None
    rec["RUN_RC"] = rc
    if rc != "0":
        probs.append(f"RUN_RC が 0 でない・無い ({rc!r})")
    if (rd / "EARLY_STOP.txt").exists():
        probs.append("EARLY_STOP.txt がある (早期停止した run)")
    ref = EV2._psha(Path(case) / EV4.REF_RUN_BIN / "RUN_PROVENANCE.txt")
    shas = {"main": EV2._psha(rd / "RUN_PROVENANCE.txt"), "soft": EV2._psha(rd / EV4.SOFT_DIR / "RUN_PROVENANCE.txt")}
    rec.update(forge_sha256=shas, forge_sha256_ref=ref)
    if not ref:
        probs.append(f"{EV4.REF_RUN_BIN}/RUN_PROVENANCE.txt の forge_sha256 を読めない")
    bad = {k: v for k, v in shas.items() if not v or v != ref}
    if bad:
        probs.append(f"forge の sha256 が {EV4.REF_RUN_BIN} ({ref}) と違う・読めない: {bad}")
    man, mp = EV4.manifest_problems(rd)
    rec["manifest"] = man
    probs += mp
    why, prec, files = prep_record_problems(rd, arm, md_repr, share)
    rec["prep"] = prec
    rec["prep_files_sha256"] = {n: files.get(n) for n in SAME_FILES}
    probs += why
    try:
        import moc_v5c_thermo_ab as VC
        rec["Tt"] = VC.read_Tt(rd)
    except Exception as e:  # noqa: BLE001
        rec["Tt"] = f"読めない: {type(e).__name__}"
    if rec["Tt"] != EV4.TT_REG:
        probs.append(f"入口の Tt が {rec['Tt']!r} ({EV4.TT_REG} でない)")
    try:
        info = _json(rd / "prepare_info.json")
    except (OSError, ValueError) as e:
        info = {}
        probs.append(f"prepare_info.json を読めない ({type(e).__name__})")
    rec["Md_moc_offset"] = info.get("Md_moc_offset")
    rec["moc"] = {k: (info.get("moc") or {}).get(k) for k in ("axis_limit", "corrector")} if isinstance(info.get("moc"), dict) else None
    rec["mesh_hashes"] = (info.get("mesh") or {}).get("hashes")
    rec["x_E"] = info.get("x_E")
    probs += prepare_info_problems(info, arm, md)
    if info:
        probs += V5.spline_shape_failures(info)
    rec["problems"] = probs
    rec["ok"] = not probs
    return rec, info


# --- 腕をまたぐ照合 -------------------------------------------------------------------------------------------------------------
def cross_run_problems(a: dict, infos: dict, fixed: dict, case: Path) -> list:
    """腕の中の run の壁・格子が同じ (同じ prep の複製; 共用の run_0164 も腕 B と同じ)、全 run の prep 時の設定ファイルが同じ、
    本段の実効設定 (V5 の effective_settings) が基準 run と同じ。"""
    bad = []
    for arm, runs in a.items():
        sp = [((infos.get(r) or {}).get("wall_fit") or {}).get("spline") for r in runs]
        mh = [((infos.get(r) or {}).get("mesh") or {}).get("hashes") for r in runs]
        if any(s is None for s in sp) or any(s != sp[0] for s in sp):
            bad.append(f"腕 {arm} の run の保存 spline がそろわない (同じ壁でない・欠損)")
        if any(h is None for h in mh) or any(h != mh[0] for h in mh):
            bad.append(f"腕 {arm} の run の格子のハッシュ (prepare_info.mesh.hashes) がそろわない (同じ格子でない・欠損)")
    sb = ((infos.get(a["B"][0]) or {}).get("wall_fit") or {}).get("spline")
    sm = ((infos.get(a["M"][0]) or {}).get("wall_fit") or {}).get("spline")
    if sb is not None and sb == sm:
        bad.append("腕 B と腕 M の保存 spline が同じ (MOC の変更が壁に入っていない)")
    for n in SAME_FILES:
        vals = {r: (fixed.get(r) or {}).get("prep_files_sha256", {}).get(n) for r in a["B"] + a["M"]}
        absent = [r for r, v in vals.items() if not v]
        # species_meta.yaml だけは全 run で無いことを許す (化学種の属性の記録が無い版の prep)。それ以外の欠けと、一部だけの欠けは不成立
        if absent and (n != "species_meta.yaml" or len(absent) != len(vals)):
            bad.append(f"prep 時の {n} の記録が欠ける run: {absent}")
        if len({v for v in vals.values() if v}) > 1:
            bad.append(f"prep 時の {n} が run の間で違う ({ {r: str(v)[:12] for r, v in vals.items()} })")
    eff = {}
    for run in a["B"] + a["M"]:
        try:
            eff[run] = V5.effective_settings(Path(case) / run)
        except Exception as e:  # noqa: BLE001 — 読めないことは不成立
            eff[run] = {"error": f"{type(e).__name__}: {e}"}
    ref = eff.get(GEOM_REF)
    if ref is None or "error" in ref:
        bad.append(f"実効設定の基準 {GEOM_REF} を読めない ({(ref or {}).get('error')})")
    for run, e in eff.items():
        if "error" in e:
            bad.append(f"{run}: 実効設定を読めない ({e['error']})")
        elif ref is not None and "error" not in ref:
            diff = [k for k in e if e[k] != ref.get(k)]
            if diff:
                bad.append(f"{run}: 本段の実効設定が {GEOM_REF} と違う ({', '.join(diff)})")
    return bad


# --- 判定 ---------------------------------------------------------------------------------------------------------------------------
def judge_exit_calibration(D: float, SE: float, tol: float = EXIT_CAL_TOL) -> str:
    """出口較正の 3 区分 (V5 と同じ境界)。「据え置き」は M = 6 の達成の証明ではない。"""
    if abs(D) + 2.0 * SE <= tol:
        return LBL_CAL_KEEP
    if abs(D) - 2.0 * SE > tol:
        return LBL_CAL_REDO
    return LBL_CAL_HOLD


def overall_verdict(verdicts: dict, rows: list, pre_bad: dict, missing: dict) -> str:
    if pre_bad:
        return "保留 (前提不成立): " + "; ".join(f"{r}: {', '.join(w)}" for r, w in pre_bad.items())
    if missing:
        return "保留 (欠損): " + "; ".join(f"{r}: {w}" for r, w in missing.items())
    if verdicts and all(v == LBL_WIN for v in verdicts.values()):
        det = [r["qty"] for r in rows if r.get("detected") and r["qty"] in DQ]
        return "許容幅内 (新しい Euler の格子で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし")
    if any(v == LBL_WORSE for v in verdicts.values()):
        return "悪化 (D − 2·SE ≥ Δq の量: " + ", ".join(k for k, v in verdicts.items() if v == LBL_WORSE) + ")"
    return "保留 (ユーザ判断: " + ", ".join(k for k, v in verdicts.items() if v != LBL_WIN) + ")"


def judge_m6(armM_runs: list, m6: dict, p3: dict, pre_ok_M: bool, pre_reasons) -> dict:
    """腕 M の M = 6 の達成 (E4 の絶対目標)。m6 = {run: m6_target}、p3 = {run: E4 の P3}。"""
    why = list(pre_reasons)
    for r in armM_runs:
        if not (p3.get(r) or {}).get("ok"):
            why.append(f"{r}: E4 の P3 (M_common の準定常・幅・平均差) が不成立 ({(p3.get(r) or {}).get('problems')})")
        if not (m6.get(r) or {}).get("complete"):
            why.append(f"{r}: 判定窓の M_common が欠ける・非有限")
    facts = {r: {k: (m6.get(r) or {}).get(k) for k in ("within", "max_abs_dev", "mean_win13")} for r in armM_runs}
    if not pre_ok_M or why:
        return {"verdict": LBL_M6_UNDET, "reasons": why, "facts": facts}
    if all(m6[r]["within"] for r in armM_runs):
        return {"verdict": LBL_M6_OK, "reasons": [], "facts": facts}
    return {"verdict": LBL_M6_NG, "reasons": [f"{r}: max|M_common − 6| {m6[r]['max_abs_dev']:.3e} > {M6_TOL:g}"
                                             for r in armM_runs if not m6[r]["within"]], "facts": facts}


def default_wall_check(case: Path, run: str, other: str | None) -> dict:
    return V5.default_wall_check(Path(case), run, other)


def default_coord_sensitivity(case: Path, run: str, geom: dict, shifts: dict) -> dict:
    """記録用 (判定に使わない): 窓の最後の枚 (step 54000) で評価座標を動かしたときの評価量の差 (V5 の同名の関数と同じ作り、step だけ違う)。"""
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.report.nozzle_report import eta_line, load_field
    ns = V5._wallfit_defs()
    ns["eta_line"] = eta_line
    F = load_field(Path(case) / run, f"res_{WIN13[-1]}.h5")

    def q_at(dE, dF):
        V5._set_geometry(ns, geom, dE, dF)
        q, _ = ns["quantities"](F, V5.SERIES_H)
        return {c: float(q[c]) for c in QS_COLS}
    base = q_at(0.0, 0.0)
    return {"step": WIN13[-1], "base": base, "diff": {n: {c: v - base[c] for c, v in q_at(*sh).items()} for n, sh in shifts.items()},
            "eval_wallfit_euler_sha256": ns["_source_sha256"]}


def evaluate(case: Path, md_text, share: bool = False, wall_check=None, coord_sensitivity=None) -> dict:
    """§6 V5d の判定。wall_check・coord_sensitivity は nozzle.h5・res_*.h5 を読む部分 (試験では差し替える)。"""
    case = Path(case).resolve()
    wall_check = wall_check or default_wall_check
    coord_sensitivity = coord_sensitivity or default_coord_sensitivity
    md, md_repr = parse_md(md_text)
    a = arms(share)
    runs_all = a["B"] + a["M"]
    role = {r: arm for arm, rs in a.items() for r in rs}
    pre_bad, metas, infos, fixed = {}, {}, {}, {}
    common = [f"評価器の窓: {p}" for p in window_problems()]
    # E4 の結果・起動の記録
    e4 = e4_adoption(case, md, share)
    common += [f"[E4] {p}" for p in e4["problems"]]
    launch, lp = launch_problems(case, md_repr, share)
    common += [f"[起動の記録] {p}" for p in lp]
    # 共通の η の列 (E4 の 12 点)
    eta_c, eta_rec = EV4.eta_common(case)
    common += [f"[η の列] {p}" for p in eta_rec.get("problems", [])]
    # 時系列の記録
    try:
        srec = _json(case / SERIES_JSON)
    except (OSError, ValueError) as e:
        srec = {}
        common.append(f"時系列の記録 {SERIES_JSON} を読めない ({type(e).__name__})")
    sp, srun = series_record_problems(srec, runs_all)
    common += sp.pop("共通", [])
    S, missing, m6, p3, own_exit = {}, {}, {}, {}, {}
    for run in runs_all:
        arm = role[run]
        why = list(sp.get(run, []))
        meta = {"arm": arm, "series_record": srun.get(run)}
        if not (case / run).is_dir():
            pre_bad[run] = ["run dir が無い (欠損)"]
            missing[run] = "run dir が無い"
            metas[run] = meta
            continue
        try:
            if run in srun:
                why += series_file_problems(case, run, srun[run])
            # 固定の条件
            fx, info = fixed_conditions(case, run, arm, md, md_repr, share)
            fixed[run], infos[run] = fx, info
            meta["fixed"] = fx
            why += [f"[固定] {p}" for p in fx["problems"]]
            # 残差 (P1)
            segf = case / run / EV4.SEGMENT_VERDICT_FILE
            conv = EV2.parse_convergence(segf.read_text() if segf.is_file() else None)
            conv["label"] = EV2.convergence_label(conv)
            meta["convergence"] = conv
            if not EV2.convergence_ok(conv):
                why.append(f"[P1] {conv['label']} — {conv.get('reasons')}")
            eff = EV2.run_steps(case / run, conv)
            meta["steps"] = eff
            why += [f"[固定] {p}" for p in step_problems(eff) if p not in window_problems()]
            # 壁
            other = a["B"][0] if arm == "M" else None
            try:
                wc = wall_check(case, run, other)
            except Exception as e:  # noqa: BLE001 — 証拠が取れないことは不成立
                wc = {"own": {"status": "error", "reason": f"{type(e).__name__}: {e}"}}
            meta["wall"] = wc
            if (wc.get("own") or {}).get("status") != "consistent":
                why.append(f"壁節点が自腕の当てはめ後 spline に一致しない ({(wc.get('own') or {}).get('status')})")
            if other is not None:
                vo = wc.get("vs_other") or {}
                nd, nown = vo.get("n_discriminable"), vo.get("n_discriminable_matching_own")
                if vo.get("status") != "ok":
                    why.append(f"腕 B との壁の照合が {vo.get('status')!r}")
                elif not (isinstance(nd, int) and nd > 0):
                    why.append("腕 B と見分けられる壁節点が無い (同じ壁の可能性)")
                elif nown != nd:
                    why.append(f"見分けられる壁節点 {nd} のうち自腕の当てはめに一致するのは {nown!r} (壁の取り違えの可能性)")
            # 全温 (P2) と M_common
            te = t0_and_exit(case, run, eta_c)
            why += [f"[P2] {p}" for p in te["problems"]]
            if not te["problems"]:
                w2, t0rec = t0_conditions(case / run, te)
                meta["t0"] = t0rec
                why += w2
                ser, p3r, own = exit_common_conditions(case / run, te)
                m6[run], p3[run], own_exit[run] = m6_target(ser), p3r, own
                meta["exit_common"] = {"series": ser, "p3_E4_record": p3r, "m6": m6[run], "exitM_own": own}
        except Exception as e:  # noqa: BLE001 — 証拠が読めないことは不成立
            why.append(f"証拠を読めない: {type(e).__name__}: {e}")
        # 評価量の窓 (準定常・幅・例外)
        try:
            S[run] = load_window(case, run)
        except Exception as e:  # noqa: BLE001 — 欠損は欠損として記録 (0 で埋めない)
            missing[run] = f"{type(e).__name__}: {e}"
        if run in S:
            try:
                ww, wrec = window_conditions(case, run, S[run])
            except Exception as e:  # noqa: BLE001
                ww, wrec = [f"準定常・幅を判定できない: {type(e).__name__}: {e}"], {"error": f"{type(e).__name__}: {e}"}
            meta["window"] = wrec
            why += ww
        metas[run] = meta
        if why:
            pre_bad[run] = why
    common += cross_run_problems(a, infos, fixed, case)
    if common:
        pre_bad["共通"] = common
    complete = not missing
    rows, verdicts, exit_cal, exit_cal_common = [], {}, None, None
    if complete:
        for k in QS_COLS:
            st = {arm: V5.arm_stats([S[r][k] for r in runs]) for arm, runs in a.items()}
            D = st["M"]["mean"] - st["B"]["mean"]
            SE = float(np.hypot(st["B"]["se"], st["M"]["se"]))
            row = {"qty": k, "B": st["B"], "M": st["M"], "D_M_minus_B": float(D), "SE_D": SE, "detected": bool(abs(D) > 2 * SE)}
            if k in DQ:
                row.update(dq=DQ[k], verdict=V5.judge_one_sided(D, SE, DQ[k]), D_plus_2SE_over_dq=float((D + 2 * SE) / DQ[k]))
                verdicts[k] = row["verdict"]
            else:
                row["verdict"] = "記録のみ (出口較正に使う)"
            rows.append(row)
        ex = next(r for r in rows if r["qty"] == "exit_core_M")
        exit_cal = {"qty": "exit_core_M (η 重みの面積平均、V5 と同じ)", "D_M_minus_B": ex["D_M_minus_B"], "SE_D": ex["SE_D"],
                    "tol": EXIT_CAL_TOL, "absD_plus_2SE": abs(ex["D_M_minus_B"]) + 2 * ex["SE_D"],
                    "absD_minus_2SE": abs(ex["D_M_minus_B"]) - 2 * ex["SE_D"],
                    "verdict": judge_exit_calibration(ex["D_M_minus_B"], ex["SE_D"]), "note": CAL_NOTE}
    # 参考: M_common (共通の η の列) による同じ 3 区分 (判定に使わない)
    if all(m6.get(r, {}).get("complete") for r in runs_all):
        stc = {arm: V5.arm_stats([np.array([metas[r]["exit_common"]["series"][s] for s in WIN13]) for r in runs]) for arm, runs in a.items()}
        Dc = stc["M"]["mean"] - stc["B"]["mean"]
        SEc = float(np.hypot(stc["B"]["se"], stc["M"]["se"]))
        exit_cal_common = {"qty": "M_common (E4 の共通の η の列、参考値)", "B": stc["B"], "M": stc["M"], "D_M_minus_B": float(Dc), "SE_D": SEc,
                           "absD_plus_2SE": abs(Dc) + 2 * SEc, "verdict_record_only": judge_exit_calibration(Dc, SEc)}
    pre_ok = not pre_bad
    overall = overall_verdict(verdicts, rows, pre_bad, missing)
    if exit_cal is not None:
        if not pre_ok:
            exit_cal["verdict_if_preconditions_held"] = exit_cal["verdict"]
            exit_cal["verdict"] = f"{LBL_CAL_HOLD} — 前提不成立"
        exit_cal["decision_complete"] = not exit_cal["verdict"].startswith("保留")
    # 腕 M の M = 6 の達成 (共通の前提と腕 M の 3 本の前提 + E4 の P3)
    pre_M = [f"{k}: {', '.join(v)}" for k, v in pre_bad.items() if k == "共通" or k in a["M"]]
    m6_judge = judge_m6(a["M"], m6, p3, not pre_M, pre_M)
    m6_B = {r: {k: (m6.get(r) or {}).get(k) for k in ("within", "max_abs_dev", "mean_win13")} | {"p3_ok": (p3.get(r) or {}).get("ok")}
            for r in a["B"]}
    coord = V5.coordinate_record(case, srec, runs_all, S, coord_sensitivity) if srec else {"error": "時系列の記録が無い"}
    out = {"plan": PLAN, "plan_registration_commit": PLAN_REG_COMMIT, "evaluator": Path(__file__).name,
           "evaluator_sha256": _sha(Path(__file__).resolve()),
           "deps_sha256": {**{n: _sha(HERE / n) for n in ("moc_v5_euler_eval.py", "e4_recal_eval.py", "euler_t0_e2_eval.py",
                                                         "euler_t0_stage_ab.py", "moc_v5c_thermo_ab.py", "throat_mono_judge.py",
                                                         "throat_mono_ab.py", "eval_wallfit_euler.py", "ic_index_map.py")},
                           **{n: _sha(TOOLS / n) for n in ("check_quasisteady.py", "check_convergence.py")}},
           "md_offset": md, "md_offset_repr": md_repr, "share_e4v": share, "arms": a, "geom_ref": GEOM_REF,
           "windows": {k: list(v) for k, v in WINDOWS.items()}, "window_fixed": True,
           "definitions": {"quantities": DQ, "width": f"各量の 13 枚の幅 ≤ Δq·{WIDTH_FRAC:g}、符号付きの出口 M の幅 ≤ {EXIT_SIGNED_WIDTH_TOL:g}",
                           "quasisteady": "6 量と出口コア M が両窓 (13 枚・5 枚) の check_quasisteady --tail 1 で STEADY (exit_M_dev だけ絶対許容内の例外)",
                           "convergence": EV2.CONV_DEFINITION, "t0": "判定窓 13 枚で全 9 領域の T0 − 1600 が ±1 K、27 列の時間の幅 ≤ 0.1 K (両窓)",
                           "exit_calibration": "出口コア M の D の 3 区分 (V5 と同じ)。M_common の同じ 3 区分は参考値",
                           "armM_M6": "腕 M の 3 本とも 13 枚すべてで |M_common − 6| ≤ 1e-4 (E4 の絶対目標、前提に E4 の P3)"},
           "e4_adoption": e4, "launch_record": launch, "eta_common": eta_rec, "series_record": SERIES_JSON,
           "series_geometry": srec.get("geometry"), "series_evaluator_sha256": srec.get("evaluator_sha256"),
           "rows": rows, "verdicts": verdicts, "overall": overall, "exit_calibration": exit_cal,
           "exit_calibration_M_common_record_only": exit_cal_common,
           "armM_M6": m6_judge, "armB_M6_record_only": m6_B,
           "preconditions_failed": pre_bad, "missing": missing, "run_meta": metas, "evaluation_coordinates": coord,
           "scope_note": SCOPE_NOTE,
           "limits": ["自己相関は補正していない (SE は実務の指標)", "窓 42000〜54000 は登録で固定 (動かさない・延長しない)",
                      "各腕 3 本は同じ prep (同じ等エントロピー IC) の複製で、run の間の差は forge の再実行の揺れだけを見る",
                      f"「{LBL_ABS}」は {EXIT_EXC_QTY} だけの例外で、STEADY ではない (元の VERDICT は run_meta に残す)",
                      f"{GEOM_TOL:g} r_t は座標の整合の許容差であって評価誤差の保証ではない (感度は evaluation_coordinates に記録)",
                      "軸 (η0) の量は判定対象外", "出口較正の据え置きと M = 6 の達成は別の判定 (armM_M6)"]}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("case", nargs="?", default=str(HERE))
    ap.add_argument("--md-offset", required=True, help="両腕に共通の Md_moc_offset (E4 の採用値の repr)")
    ap.add_argument("--share-e4v", action="store_true", help=f"腕 B に {E4V_RUN} を共用した構成")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case = Path(a.case).resolve()
    out = evaluate(case, a.md_offset, a.share_e4v)
    op = Path(a.out) if a.out else case / OUT_JSON
    op.parent.mkdir(parents=True, exist_ok=True)
    op.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(f"評価器 {out['evaluator']} sha256 {out['evaluator_sha256']}  (登録 {PLAN_REG_COMMIT})  較正値 {out['md_offset_repr']}  "
          f"共用 {out['share_e4v']}  腕 B {out['arms']['B']}  腕 M {out['arms']['M']}")
    print(f"{'量':28s} {'Δq':>8s} {'B 平均':>11s} {'M 平均':>11s} {'D=M−B':>10s} {'2·SE_D':>9s} {'(D+2SE)/Δq':>10s} {'検出':>4s}  判定")
    for r in out["rows"]:
        print(f"{r['qty']:28s} {r.get('dq', float('nan')):8.2e} {r['B']['mean']:11.7g} {r['M']['mean']:11.7g} {r['D_M_minus_B']:+10.2e} "
              f"{2 * r['SE_D']:9.2e} {r.get('D_plus_2SE_over_dq', float('nan')):10.2f} {'あり' if r['detected'] else 'なし':>4s}  {r['verdict']}")
    for run, m in out["run_meta"].items():
        st = (m.get("window") or {}).get("status") or {}
        ns = {c: v for c, v in st.items() if v != "STEADY"}
        print(f"準定常 {run}: " + ("全量 STEADY (両窓)" if st and not ns else ", ".join(f"{c} {v}" for c, v in ns.items()) if ns else "判定できない")
              + f" | 残差 {(m.get('convergence') or {}).get('label')}")
    for r, w in out["preconditions_failed"].items():
        print(f"前提不成立 {r}: " + " / ".join(w[:8]) + (f" / … (他 {len(w) - 8} 件)" if len(w) > 8 else ""))
    if out["exit_calibration"]:
        ec = out["exit_calibration"]
        print(f"出口較正 (出口コア M): D {ec['D_M_minus_B']:+.3e}  |D| + 2SE {ec['absD_plus_2SE']:.3e} (≤ {ec['tol']:g})  → {ec['verdict']}")
    if out["exit_calibration_M_common_record_only"]:
        ec = out["exit_calibration_M_common_record_only"]
        print(f"  参考 (M_common): D {ec['D_M_minus_B']:+.3e}  |D| + 2SE {ec['absD_plus_2SE']:.3e}  → {ec['verdict_record_only']}")
    print(f"腕 M の M = 6: {out['armM_M6']['verdict']}  {out['armM_M6']['facts']}")
    print("総合:", out["overall"])
    print("注:", SCOPE_NOTE)
    print("->", op)
    return 0


if __name__ == "__main__":
    sys.exit(main())
