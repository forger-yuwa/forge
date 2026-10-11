"""⑤ SERN の多作動点 MOO ドライバ (S6)。

dv = (M_c, f, theta_r0_deg, theta_c0_deg, L_cowl) [YAML の dv.min/max]。1 評価 = 作動点数分の forge run
(`spec.operating_points[]`、重み w_k)。目的 (最小化):
    f1 = −Σ_k w_k C_T^(k)      (Kriging。RANS では摩擦込み C_T)
    f2 = L_ramp / H            (逆設計の決定的結果。候補点では Kriging で近似)
ゲート (plan §5.1 R1, `metrics/sern_gates.py`): 逆設計成立 (kernel 内に key point、C⁻ が直線ランプに着地、
`L_ramp_max` は最終輪郭で再検査)、メッシュ品質 PASS、各作動点で **rc == 0 かつ 保存場が有限・正値、全残差に NaN/rising 無し、
実目的量 (RANS は C_T_with_shear) と C_T/C_L/C_M が STEADY**。発散した評価の採用 (旧 §4.13-3) は撤回した:
ゲート不合格はサロゲート学習と Pareto から外す (accepted `tooling-nozzle-moo-loop.md` と同じ)。学習・Pareto に使うのは現行の受理検査の版
(sern_gates.ACCEPTANCE_VERSION) で判定し、残差の入力の検査済み要約が合格の行だけ (plan tooling-sern-te-wake-grid §5.1 R2k。旧い版の行は
台帳に残り、`--rejudge` で付け直す)。
status: PASS / INFEASIBLE (物理: 設計不成立 `DESIGN`, `L_RAMP_MAX`, `CM_WINDOW`) / FAIL (数値: `DIVERGED`,
`RESIDUAL_RISING`, `NOT_CONVERGED`, `UNSTEADY`, `NO_FORCES`, `ERROR`)。`degraded` = 緩レシピ (再試行) で通った評価
(台帳と pareto.json の両方に残す)。C_M (機体 CG 基準) は台帳に記録し制約は任意 (`opt.cm_min/cm_max` = 加重平均、
`opt.cm_window: {op: [lo, hi]}` = 作動点別の窓 [R6(d)])。
既存 `opt/` (lhs / KrigingSet / propose_infill / hypervolume2d) を再利用。使い方:

  design/.venv-opt/bin/python -m forge_design.opt.driver_sern <problem.yaml> <campaign_dir> [--n-doe 12 --n-iter 3 --batch 2]
  design/.venv-opt/bin/python -m forge_design.opt.driver_sern <problem.yaml> <campaign_dir> --rejudge <out_dir>
      (CFD を回さず既存 run を現行ゲートで再判定: out_dir に ledger_rejudged.jsonl / pareto_rejudged.json / rejudge_summary.md)
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

from ..evaluate import runner_sern as R
from ..geometry.moc_sern import PlanarMOC, SernKernelSpec
from ..probdef import load_problem
from .doe import lhs
from .ehvi import hypervolume2d, nondominated_mask
from .moo import propose_infill
from .surrogate import KrigingSet

DV_ORDER = ("M_c", "f", "theta_r0_deg", "theta_c0_deg", "L_cowl")
GATE_KEYS = ("C_T", "C_L", "C_M")


# slauWallNormalChi の既定変更 (runner_sern.FLAG_POLICY と同じ値)
from ..evaluate.runner_sern import FLAG_POLICY  # noqa: E402
# 処置の識別 (import 時に束縛する: 試験は R を偽の runner に差し替えてからキャンペーンを作る)
from ..evaluate.runner_sern import ZTE_EFFECTIVE, zte_signature_of, zte_spec_from_evaluate  # noqa: E402
from ..metrics.sern_gates import floor_events_summary  # noqa: E402  毎更新の床事象の要約 (te-wake-grid §5.1 #2)
from ..metrics.sern_gates import floor_events_options  # noqa: E402  キャンペーンの床判定の要件 (te-wake-grid §5.1 R1)
from ..metrics.sern_gates import JUDGMENT_MIN_SAMPLES, TE_REGIONS, judgment_summary  # noqa: E402  判定区間の §6 #2・#3 (te-wake-grid §5.1 R2b)
from ..evaluate.runner_sern import judgment_plan  # noqa: E402  キャンペーンの判定区間の要件 (import 時に束縛する)
# 受理検査の版と残差の入力の検査済み要約 (te-wake-grid §5.1 R2k)
from ..metrics.sern_gates import ACCEPTANCE_VERSION, ACCEPTANCE_VERSION_KEY, RESIDUAL_FLUID_COLS, residual_evidence_summary  # noqa: E402
# 評価方式の識別 (plan tooling-sern-te-wake-grid §4、codex plan レビュー 2026-10-08 M3。import 時に束縛する)
from ..evaluate.runner_sern import EVAL_METHOD_ID, TE_WAKE_EFFECTIVE, eval_method_required, te_wake_effective_of, te_wake_spec  # noqa: E402

class DesignInfeasible(ValueError):
    """物理的に成立しない候補 (逆設計不成立 / L_ramp_max 超過)。数値失敗と区別する (R1)。"""


class EvalFailure(RuntimeError):
    """CFD 評価の数値失敗 (ゲート不合格)。fail_class を持つ。"""

    def __init__(self, fail_class: str, msg: str) -> None:
        super().__init__(msg); self.fail_class = fail_class


class _KrgBoth:
    def __init__(self, krg: KrigingSet) -> None:
        self.krg = krg

    def predict(self, X):
        return self.krg.predict(np.atleast_2d(np.asarray(X, dtype=float)))


REQUIRED_SCALAR_GRADIENT = "lsq"
REQUIRED_WALL_NORMAL_CHI = 1   # 2026-09-26 から node+SLAU の既定 (auto)。明示 0 の評価は別の応答関数なので混ぜない
# 厚さ 0 の板の自由端の速度再構成の処置 (plan convection-zero-thickness-edge-reconstruction §4「設計チェーンへの配線」):
# 作動点ごとの実効値 ("off" / "tags=...;rings=N" / None = 不明) がキャンペーンの要求 (基準 YAML の
# evaluate.zero_thickness_edge_velocity) と一致する行だけを学習・Pareto に使う。処置の導入前の行 (キーが無い) は "off"。
# **処置なしの旧評価を処置ありのキャンペーンへ持ち越すのは、作動点を名指しした明示の指定
# (`opt.zero_thickness_edge_velocity_carry_over: [m6_on, m4_off]`) があるときだけ** (可否の判定は plan §6.2 #6・#7)
ZTE_CARRY_OVER_KEY = "zero_thickness_edge_velocity_carry_over"


def _zte_effective_of(o: dict):
    """作動点要約の処置の実効値。キーが無い行 (処置の導入前の評価) は "off"、キーがあって None なら不明。"""
    return o[ZTE_EFFECTIVE] if ZTE_EFFECTIVE in o else "off"


def _zte_compatible(op: str, o: dict, required: str = "off", carry_over=()) -> bool:
    eff = _zte_effective_of(o)
    if eff is None:
        return False
    return eff == required or (eff == "off" and op in carry_over)


def _learnable_base(r: dict) -> bool:
    """処置以外の採否: PASS・現行の flag_policy・全作動点の実効 scalarGradient が lsq かつ実効 chi が 1。
    不明 (None) は除外 (codex result 2026-09-27 chi-default M4: 日付と scalarGradient だけでは chi 0/1/不明が混ざる)。"""
    if r.get("status") != "PASS" or r.get("flag_policy") != FLAG_POLICY:
        return False
    ops = r.get("ops") or {}
    return bool(ops) and all(o.get("scalar_gradient_effective") == REQUIRED_SCALAR_GRADIENT
                             and o.get("slau_wall_normal_chi_effective") in (REQUIRED_WALL_NORMAL_CHI, str(REQUIRED_WALL_NORMAL_CHI))
                             for o in ops.values())


def _zte_row_ok(r: dict, zte_required: str = "off", zte_carry_over=()) -> bool:
    """全作動点の処置の実効値がキャンペーンの要求と一致 (処置なしは明示の持ち越し指定の作動点だけ可)。不明は不可。"""
    ops = r.get("ops") or {}
    return bool(ops) and all(_zte_compatible(name, o, zte_required, zte_carry_over) for name, o in ops.items())


# --- 評価方式の識別 (plan tooling-sern-te-wake-grid §4「格子の来歴と評価方式の識別を分ける」、codex plan レビュー 2026-10-08 M3) ---
# 学習可否は**評価方式の識別** (runner_sern.eval_method_*: 離散化・格子生成のレシピ・メッシャの版・後縁下流の格子変形 L_b と
# 曲線版・変換器の版) で決める。実格子の署名・双対幾何のハッシュは設計点ごとに変わるので来歴として残すだけで照合しない。
#   - 後縁下流の格子変形を要求しないキャンペーン (基準 YAML の te_wake_blend_H が無い・0 = 修正前の系列): 従来どおりで、
#     変形を入れた評価・実効値が不明な評価だけを外す (識別の導入前の行は変形なし = off)
#   - 要求するキャンペーン (> 0 = 修正後の学習系列): 全作動点の評価方式の id がキャンペーンの要求 (基準 YAML・現在のコード・
#     現在の変換器) と一致する行だけ。旧方式・変換器の版違い・来歴不明 (識別の無い旧行・None) は外す (台帳には残す)
# 旧評価の持ち越しは**識別の一致とは独立の承認**: `opt.eval_method_carry_over` に作動点・持ち越す旧評価の評価方式
# (id か識別の導入前を表す "legacy")・plan §6 #6 (旧評価の持ち越し) と §6 #7 (設計差) の判定・根拠を書いた項目のうち、
# 両判定が PASS のものだけが効く。承認は処置以外の採否 (_learnable_base) と処置の一致 (_zte_row_ok) を緩めない
EVAL_METHOD_CARRY_OVER_KEY = "eval_method_carry_over"
_CARRY_KEYS = ("op", "from_eval_method", "check_6", "check_7", "evidence")
_CARRY_VERDICTS = ("PASS", "FAIL", "UNDECIDABLE")
LEGACY_METHOD = "legacy"


def _te_wake_of(o: dict):
    """作動点要約の後縁下流の格子変形の実効値。キーが無い行 (識別の導入前の評価) は "off"、キーがあって None なら不明。"""
    return o[TE_WAKE_EFFECTIVE] if TE_WAKE_EFFECTIVE in o else "off"


def _method_carried(op: str, o: dict, approved=()) -> bool:
    """承認済みの持ち越し (両判定 PASS の項目だけが approved に入っている) に当たるか。"""
    for a in approved:
        if a["op"] != op:
            continue
        if a["from_eval_method"] == LEGACY_METHOD:
            if EVAL_METHOD_ID not in o:          # 識別の導入前の評価 (None = 記録はあるが不明、は持ち越さない)
                return True
        elif o.get(EVAL_METHOD_ID) == a["from_eval_method"]:
            return True
    return False


def _method_row_ok(r: dict, te_wake_required: str = "off", method_required=None, approved=()) -> bool:
    """評価方式の一致。要求なし (off) は格子変形を入れた評価・不明だけを外す。要求あり は全作動点の id が一致
    (または承認済みの持ち越し)。要求ありで要求の id が不明 (None) なら、持ち越し以外は採らない。"""
    ops = r.get("ops") or {}
    if not ops:
        return False
    if te_wake_required == "off":
        return all(_te_wake_of(o) == "off" for o in ops.values())
    return all((method_required is not None and o.get(EVAL_METHOD_ID) == method_required) or _method_carried(name, o, approved)
               for name, o in ops.items())


# --- 床事象の判定の要件 (plan tooling-sern-te-wake-grid §5.1 R1、codex result レビュー 2026-10-11 M1) ---
# 基準 YAML の `evaluate.floor_events: 1` のキャンペーンは、各作動点に**保存された**床事象の判定 (作動点要約の floor_events =
# sern_gates.floor_events_summary) が PASS・判定区間の長さ ≥ 要件・検査した更新数 ≥ 要件の行だけを学習・Pareto に使う。
# 要件の長さは `evaluate.floor_events_window_steps` (省略時は区間が空でないことだけ。区間は末尾の割合で決まる)。
# 評価時にゲートが任意 (required: false) だった行でも、保存された判定が要件を満たせば採る (評価時の設定ではなく判定を照合する)。
# 欠測 (要約が無い・None・区間や更新数が読めない)・短い区間・FAIL・判定不能は外す。台帳には残り、`--rejudge` で
# 現行ゲートの判定を付け直せば戻せる。旧評価の持ち越しの承認 (eval_method_carry_over) はこの要件を緩めない。
# floor_events: 0 (またはキー無し) のキャンペーンは要件なし = 従来どおり
def floor_requirement_of(evaluate) -> dict | None:
    """基準 YAML の evaluate 節からキャンペーンの床判定の要件を作る。必須でなければ None (要件なし)。
    戻り = {"required": True, "window_steps": 判定区間の長さの下限 [step] (None = 指定なし), "min_updates": 検査した更新数の下限}。
    定常陰解法は 1 step 1 更新なので、更新数の下限は区間長と同じ (指定なしは 1)。"""
    o = floor_events_options(evaluate)          # 値の検査もここで (0/1・正の整数)
    if not o["require_floor_events"]:
        return None
    w = o["floor_window_steps"]
    return {"required": True, "window_steps": w, "min_updates": int(w) if w is not None else 1}


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _floor_op_ok(o: dict, req: dict) -> bool:
    """作動点要約の floor_events が要件を満たすか: verdict PASS・区間 [a, b] (整数) の長さ b − a ≥ 要件 (かつ > 0)・
    検査した更新数 ≥ 要件。読めない・欠けた要約は満たさない (記録が無いことを 0 件として扱わない)。"""
    fe = o.get("floor_events")
    if not isinstance(fe, dict) or fe.get("verdict") != "PASS":
        return False
    win, n = fe.get("window"), fe.get("n_updates_checked")
    if not (isinstance(win, (list, tuple)) and len(win) == 2 and all(_is_int(v) for v in win) and _is_int(n)):
        return False
    need_len = req.get("window_steps") or 1
    return (win[1] - win[0]) >= max(int(need_len), 1) and n >= max(int(req.get("min_updates") or 1), 1)


def _floor_row_ok(r: dict, req: dict | None = None) -> bool:
    """全作動点の床事象の判定がキャンペーンの要件を満たす。要件なし (None) は従来どおり (見ない)。"""
    if req is None:
        return True
    ops = r.get("ops") or {}
    return bool(ops) and all(_floor_op_ok(o, req) for o in ops.values())


# --- 判定区間の §6 #2・#3 の要件 (plan tooling-sern-te-wake-grid §5.1 R2b、codex diagnose 2026-10-12 M2) ---
# 基準 YAML に `evaluate.judgment_steps` があるキャンペーンは、各作動点に**保存された**判定区間の判定 (作動点要約の judgment と
# judgment_summary = sern_gates.judgment_summary) が PASS で、判定窓 ≥ W (力の窓・温度の窓とも区間長 ≥ W、標本 20 以上)、
# キャンペーンの必須量 (evaluate.model から: SST 4 量・Euler 3 量) が全部判定されて合格、温度の 4 領域がそろう行だけを学習・Pareto に使う。
# #1 (床事象) も要件に含める (judgment_steps は evaluate.floor_events: 1 と窓を要求するので、上の床の要件と同じものを持つ)。
# 判定の無い行 (判定区間の段の導入前・旧 GATES PASS) は欠測として外す。**旧評価の持ち越しの承認 (eval_method_carry_over) は
# この要件を緩めない** (旧 ID の承認だけでは学習に戻らない。保存系列を現行の判定で再判定 [--rejudge] して要約を付け直す)。
# judgment_steps が無いキャンペーンは要件なし = 従来どおり
def judgment_requirement_of(evaluate) -> dict | None:
    """基準 YAML の evaluate 節からキャンペーンの判定区間の要件を作る。judgment_steps が無ければ None (要件なし)。
    戻り = {"window_steps": W, "quantities": 必須量, "regions": 温度の監視領域, "min_samples": 窓あたりの標本の下限,
    "floor": #1 の床事象の要件 (floor_requirement_of)}。"""
    jp = judgment_plan(evaluate)              # 値の検査もここで (runner の collect と同じ予定)
    if jp is None:
        return None
    return {"window_steps": int(jp["window_steps"]), "quantities": list(jp["required"]), "regions": list(TE_REGIONS),
            "min_samples": int(JUDGMENT_MIN_SAMPLES), "floor": floor_requirement_of(evaluate)}


def _window_ok(win, need: int) -> bool:
    return (isinstance(win, (list, tuple)) and len(win) == 2 and all(_is_int(v) for v in win) and (win[1] - win[0]) >= max(int(need), 1))


def _judgment_op_ok(o: dict, req: dict) -> bool:
    """作動点要約の判定区間の判定が要件を満たすか: #1 床事象 (_floor_op_ok)・judgment PASS・要約の verdict PASS・窓の長さ ≥ W・
    力の窓 PASS (区間長 ≥ W、両窓の標本 ≥ 下限、必須量がすべて判定され合格)・温度の窓 PASS (区間長 ≥ W、標本 ≥ 下限、領域がそろう)。
    欠けた要約は満たさない。"""
    if not (isinstance(req.get("floor"), dict) and _floor_op_ok(o, req["floor"])):
        return False
    js = o.get("judgment_summary")
    if o.get("judgment") != "PASS" or not isinstance(js, dict) or js.get("verdict") != "PASS":
        return False
    W, nmin = int(req["window_steps"]), int(req.get("min_samples") or 1)
    if not (_is_int(js.get("window_steps")) and js["window_steps"] >= W):
        return False
    fw, tp = js.get("force_window"), js.get("temperature")
    if not (isinstance(fw, dict) and isinstance(tp, dict)):
        return False
    oks = fw.get("ok") if isinstance(fw.get("ok"), dict) else {}
    if not (fw.get("verdict") == "PASS" and _window_ok(fw.get("window"), W) and _is_int(fw.get("n_tail")) and _is_int(fw.get("n_prev"))
            and fw["n_tail"] >= nmin and fw["n_prev"] >= nmin and all(oks.get(q) is True for q in req["quantities"])):
        return False
    regions = tp.get("regions") if isinstance(tp.get("regions"), list) else []
    return (tp.get("verdict") == "PASS" and _window_ok(tp.get("window"), W) and _is_int(tp.get("n_window")) and tp["n_window"] >= nmin
            and set(req["regions"]) <= set(regions))


def _judgment_row_ok(r: dict, req: dict | None = None) -> bool:
    """全作動点の判定区間の判定がキャンペーンの要件を満たす。要件なし (None) は従来どおり (見ない)。"""
    if req is None:
        return True
    ops = r.get("ops") or {}
    return bool(ops) and all(_judgment_op_ok(o, req) for o in ops.values())


# --- 受理検査の版と残差の入力の証拠 (plan tooling-sern-te-wake-grid §5.1 R2k、codex result 再レビュー 3 回目 2026-10-12 M1) ---
# 学習・Pareto は保存された判定 (status PASS) を使うので、受理検査を厳しくしても、旧い検査で PASS した評価が再判定なしで入っていた
# (9e5f0ffd の判定で PASS した、化学種の残差列の欠落・1 行だけ・途中の欠測/重複の残差。codex の再現)。そこで各作動点に
#   (1) 受理検査の版 (作動点要約の acceptance_version = 評価時の gates の値) が現行の sern_gates.ACCEPTANCE_VERSION と一致し、
#   (2) 残差の入力の検査済み要約 (residual_evidence = sern_gates.residual_evidence_summary) が合格: 入力が妥当 (input_ok True・判定不能でない)、
#       要求した残差列が流体 5 列を含み欠けた列が無い、時系列が session の nStepOuter と最終場に対して完全 (outer_end が step 0..N−1 の N 行・最終場 res_N)
# であることを要求する。**キャンペーンの設定によらず一律**: 床事象・判定区間 (R1・R2b) は問題 YAML が要求するときだけ回す検査で、窓の長さも
# キャンペーンが決めるのでキャンペーンの要求に置いたが、ここで照合する検査は evaluate_gates が全評価に課す必須の検査なので、要求しない
# キャンペーン (床事象・判定区間を持たない 2D など) を作れると旧い検査の PASS が黙って戻る。R1 と同じく、照合するのは評価時の設定でなく
# 保存された判定と要約、版・要約の無い行 (R2k の前の台帳の全行)・旧い版の行は台帳に残して外し、旧評価の持ち越しの承認
# (eval_method_carry_over・処置の持ち越し) は緩めない。戻すのは `--rejudge` (現行ゲートで判定し直した行は現行の版と要約を持つ)
def _residual_evidence_ok(ev) -> bool:
    """残差の入力の検査済み要約が合格か (保存された要約の欄だけを見る。run は読み直さない)。欄の欠け・None は不合格。"""
    if not isinstance(ev, dict) or ev.get("input_ok") is not True or ev.get("unverifiable") is not False:
        return False
    req, miss = ev.get("required_columns"), ev.get("missing_columns")
    if not (isinstance(req, list) and all(isinstance(c, str) for c in req) and set(RESIDUAL_FLUID_COLS) <= set(req)
            and isinstance(miss, list) and not miss):
        return False
    s = ev.get("series")
    if not isinstance(s, dict):
        return False
    n = s.get("n_steps")
    return (_is_int(n) and n > 0 and _is_int(s.get("n_outer_end")) and s["n_outer_end"] == n
            and _is_int(s.get("first")) and s["first"] == 0 and _is_int(s.get("last")) and s["last"] == n - 1
            and s.get("final_field") == f"res_{n}.h5")


def _acceptance_op_ok(o: dict) -> bool:
    """作動点要約が現行の受理検査の版で判定され、残差の入力の検査済み要約が合格か。"""
    return o.get(ACCEPTANCE_VERSION_KEY) == ACCEPTANCE_VERSION and _residual_evidence_ok(o.get("residual_evidence"))


def _acceptance_row_ok(r: dict) -> bool:
    """全作動点が現行の受理検査の版と残差の入力の要約の合格を持つ (キャンペーンの設定によらない)。"""
    ops = r.get("ops") or {}
    return bool(ops) and all(_acceptance_op_ok(o) for o in ops.values())


def _learnable(r: dict, zte_required: str = "off", zte_carry_over=(), method=None, floor=None, judgment=None) -> bool:
    """学習・Pareto に使ってよい行 = 処置以外の採否 (_learnable_base) かつ受理検査の版と残差の入力の要約 (_acceptance_row_ok、一律)
    かつ処置の一致 (_zte_row_ok) かつ評価方式の一致 (_method_row_ok) かつ床事象の判定の要件 (_floor_row_ok) かつ判定区間の要件
    (_judgment_row_ok)。method = (格子変形の要求, 評価方式の id の要求, 承認済みの持ち越し)、None は要求なし。
    floor = floor_requirement_of・judgment = judgment_requirement_of の戻り (None は要件なし)。"""
    return (_learnable_base(r) and _acceptance_row_ok(r) and _zte_row_ok(r, zte_required, zte_carry_over)
            and _method_row_ok(r, *(method or ("off", None, ()))) and _floor_row_ok(r, floor) and _judgment_row_ok(r, judgment))


def _zte_policy(obj):
    """キャンペーン (または試験の偽物) の処置の要求と持ち越し指定。属性が無ければ処置なし・持ち越しなし。"""
    return getattr(obj, "zte_required", "off"), tuple(getattr(obj, "zte_carry_over", ()))


def _method_policy(obj):
    """キャンペーン (または試験の偽物) の評価方式の要求: (格子変形の要求, 評価方式の id の要求, 承認済みの持ち越し)。
    属性が無ければ要求なし (off)。"""
    return (getattr(obj, "te_wake_required", "off"), getattr(obj, "eval_method_required", None),
            tuple(getattr(obj, "eval_method_carry_over", ())))


def _floor_policy(obj):
    """キャンペーン (または試験の偽物) の床判定の要件。属性が無ければ要件なし (None)。"""
    return getattr(obj, "floor_required", None)


def _judgment_policy(obj):
    """キャンペーン (または試験の偽物) の判定区間の要件。属性が無ければ要件なし (None)。"""
    return getattr(obj, "judgment_required", None)


def parse_eval_method_carry_over(items, op_names, te_wake_required: str) -> tuple:
    """`opt.eval_method_carry_over` を検査する。戻り = (全項目, 承認済み [check_6・check_7 とも PASS] の項目)。
    項目 = {op, from_eval_method (16 進 64 桁か "legacy"), check_6 (plan §6 #6 旧評価の持ち越し), check_7 (§6 #7 設計差),
    evidence (判定の記録の所在)}。判定は PASS / FAIL / UNDECIDABLE。格子変形を要求しないキャンペーンでは書けない。"""
    items = items or []
    if not isinstance(items, list) or not all(isinstance(a, dict) for a in items):
        raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY} は項目 (dict) のリスト: {items!r}")
    if items and te_wake_required == "off":
        raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY} は後縁下流の格子変形を要求するキャンペーン (te_wake_blend_H > 0) でだけ意味がある")
    for a in items:
        unknown, missing = sorted(set(a) - set(_CARRY_KEYS)), [k for k in _CARRY_KEYS if k not in a]
        if unknown or missing:
            raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY} の項目のキー (未知 {unknown}・欠け {missing}): {a!r}")
        if a["op"] not in op_names:
            raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY}: 無い作動点 {a['op']!r}")
        fm = a["from_eval_method"]
        if not (fm == LEGACY_METHOD or (isinstance(fm, str) and re.fullmatch(r"[0-9a-f]{64}", fm))):
            raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY}: from_eval_method は評価方式の id (16 進 64 桁) か {LEGACY_METHOD!r}: {fm!r}")
        for k in ("check_6", "check_7"):
            if a[k] not in _CARRY_VERDICTS:
                raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY}: {k} は {_CARRY_VERDICTS} のどれか: {a[k]!r}")
        if not (isinstance(a["evidence"], str) and a["evidence"].strip()):
            raise ValueError(f"opt.{EVAL_METHOD_CARRY_OVER_KEY}: evidence (判定の記録の所在) が空: {a!r}")
    items = [dict(a) for a in items]
    return items, tuple(a for a in items if a["check_6"] == "PASS" and a["check_7"] == "PASS")


class SernCampaign:
    def __init__(self, base_yaml, campaign_dir, ref=(-0.90, 20.0), seed: int = 0) -> None:
        self.base_yaml = Path(base_yaml); self.dir = Path(campaign_dir); self.dir.mkdir(parents=True, exist_ok=True)
        self.ref = tuple(ref); self.seed = int(seed)
        self.base_raw = yaml.safe_load(self.base_yaml.read_text())
        self.bounds = np.array([[float(self.base_raw["dv"][k]["min"]), float(self.base_raw["dv"][k]["max"])] for k in DV_ORDER])
        self.ops = self.base_raw["spec"].get("operating_points") or [{"name": "default", "weight": 1.0, "external": self.base_raw["spec"]["external"]}]
        self.optcfg = self.base_raw.get("opt", {})
        self.ledger = self.dir / "ledger.jsonl"
        self.rows = [json.loads(l) for l in self.ledger.read_text().splitlines() if l.strip()] if self.ledger.exists() else []
        # 処置の要求 (基準 YAML) と、処置なしの旧評価の持ち越し指定 (作動点名のリスト、既定なし)
        self.zte_required = zte_signature_of(zte_spec_from_evaluate(self.base_raw.get("evaluate") or {}))
        _co = self.optcfg.get(ZTE_CARRY_OVER_KEY) or []
        if not isinstance(_co, list) or not all(isinstance(o, str) for o in _co):
            raise ValueError(f"opt.{ZTE_CARRY_OVER_KEY} は作動点名のリスト: {_co!r}")
        _unknown = sorted(set(_co) - {o["name"] for o in self.ops})
        if _unknown:
            raise ValueError(f"opt.{ZTE_CARRY_OVER_KEY} に無い作動点: {_unknown}")
        if _co and self.zte_required == "off":
            raise ValueError(f"opt.{ZTE_CARRY_OVER_KEY} は処置ありのキャンペーン (evaluate.zero_thickness_edge_velocity) でだけ意味がある")
        self.zte_carry_over = tuple(_co)
        # 評価方式の要求 (driver は 2D runner)。後縁下流の格子変形を要求するときだけ id を照合する (修正後の学習系列)
        self.te_wake_required = te_wake_effective_of(*te_wake_spec(self.base_raw, 2))
        self.eval_method_required = None
        if self.te_wake_required != "off":
            self.eval_method_required = eval_method_required(load_problem(self.base_yaml), 2)[EVAL_METHOD_ID]
            if self.eval_method_required is None:
                print("[campaign] 要求する評価方式の id を作れない (変換器が無い等)。承認済みの持ち越し以外は学習しない", flush=True)
        self.eval_method_carry_over_records, self.eval_method_carry_over = parse_eval_method_carry_over(
            self.optcfg.get(EVAL_METHOD_CARRY_OVER_KEY), {o["name"] for o in self.ops}, self.te_wake_required)
        # 床事象の判定の要件 (基準 YAML の evaluate.floor_events / floor_events_window_steps、plan tooling-sern-te-wake-grid §5.1 R1)
        self.floor_required = floor_requirement_of(self.base_raw.get("evaluate") or {})
        # 判定区間の §6 #2・#3 の要件 (基準 YAML の evaluate.judgment_steps、plan tooling-sern-te-wake-grid §5.1 R2b)
        self.judgment_required = judgment_requirement_of(self.base_raw.get("evaluate") or {})
        _old = sum(1 for r in self.rows if r.get("status") == "PASS" and r.get("flag_policy") != FLAG_POLICY)
        if _old:
            print(f"[campaign] flag_policy が {FLAG_POLICY} でない PASS 行 {_old} 件を学習から除外 "
                  f"(既定変更前の評価: slauWallNormalChi 2026-09-26 / scalarGradient 2026-09-27)", flush=True)
        _nolsq = sum(1 for r in self.rows if r.get("status") == "PASS" and r.get("flag_policy") == FLAG_POLICY and not _learnable_base(r))
        if _nolsq:
            print(f"[campaign] 実効 scalarGradient (lsq) / slauWallNormalChi (1) が全作動点で確認できない PASS 行 {_nolsq} 件を学習から除外", flush=True)
        _noacc = sum(1 for r in self.rows if _learnable_base(r) and not _acceptance_row_ok(r))
        if _noacc:
            print(f"[campaign] 受理検査の版が現行 ({ACCEPTANCE_VERSION}) でない・残差の入力の検査済み要約が無いか不合格の PASS 行 {_noacc} 件を"
                  f"学習から除外 (旧い受理検査の判定。台帳には残る。--rejudge で現行ゲートの判定を付け直せる)", flush=True)
        _nozte = sum(1 for r in self.rows if _learnable_base(r) and not _zte_row_ok(r, self.zte_required, self.zte_carry_over))
        if _nozte:
            print(f"[campaign] 処置 (zero_thickness_edge_velocity) の実効値が要求 {self.zte_required} と違う・不明な PASS 行 {_nozte} 件を学習から除外"
                  f" (持ち越し指定: {list(self.zte_carry_over) or 'なし'})", flush=True)
        _nometh = sum(1 for r in self.rows if _learnable_base(r) and _zte_row_ok(r, self.zte_required, self.zte_carry_over)
                      and not _method_row_ok(r, *_method_policy(self)))
        if _nometh:
            print(f"[campaign] 評価方式 (格子変形 {self.te_wake_required}・id {self.eval_method_required}) が要求と違う・不明な PASS 行 "
                  f"{_nometh} 件を学習から除外 (承認済みの持ち越し: {[a['op'] for a in self.eval_method_carry_over] or 'なし'})", flush=True)
        _nofloor = sum(1 for r in self.rows if _learnable(r, self.zte_required, self.zte_carry_over, method=_method_policy(self))
                       and not _floor_row_ok(r, self.floor_required))
        if _nofloor:
            print(f"[campaign] 床事象の判定が要件 {self.floor_required} を満たさない (欠測・短い区間・FAIL・判定不能) PASS 行 "
                  f"{_nofloor} 件を学習から除外 (台帳には残る。--rejudge で付け直せる)", flush=True)
        _nojudg = sum(1 for r in self.rows if _learnable(r, self.zte_required, self.zte_carry_over, method=_method_policy(self),
                                                         floor=self.floor_required)
                      and not _judgment_row_ok(r, self.judgment_required))
        if _nojudg:
            print(f"[campaign] 判定区間の §6 #2・#3 が要件 {self.judgment_required} を満たさない (判定なし・未達・判定不能・短い窓・"
                  f"必須量の欠け) PASS 行 {_nojudg} 件を学習から除外 (台帳には残る。持ち越しの承認では戻らない。--rejudge で付け直せる)",
                  flush=True)

    def _write_problem(self, x, path: Path) -> Path:
        raw = json.loads(json.dumps(self.base_raw))
        for k, v in zip(DV_ORDER, np.asarray(x, dtype=float)):
            raw["dv"][k]["value"] = float(v)
        path.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False)); return path

    # -- 逆設計の成立判定 (粗い kernel で高速化) ----------------------------------
    def feasible(self, x) -> bool:
        try:
            self._design_probe(x); return True
        except Exception:
            return False

    def _design_probe(self, x):
        raw = self.base_raw; geo = raw["geometry"]; st_in = raw["spec"]["inflow"]
        ex0 = raw["spec"]["external"]      # 設計点 (作動点ではない)
        g = float(raw["gas"].get("gamma", 1.4)); cp = float(raw["gas"].get("cp", 1004.5)); Rg = cp * (g - 1) / g
        p_ratio = float(ex0["p_inf"]) / float(st_in["p_in"])
        M_c, f, tr, tc, Lc = [float(v) for v in x]
        spec = SernKernelSpec(M_in=float(st_in["M_in"]), theta_r0=np.deg2rad(tr), theta_c0=np.deg2rad(tc), L_cowl=Lc, gamma=g,
                              p_ext_over_p_in=p_ratio, x_max=float(geo.get("x_max_kernel", 10.0)),
                              nj=int(geo.get("nj_moc_probe", 151)), dx=float(geo.get("dx_moc_probe", 4e-3)))
        k = PlanarMOC(spec).march(stop_at=(f, M_c)); d = k.design_ramp(M_c=M_c, f=f, ds=2e-3)
        if d.info["warnings"]:
            raise ValueError(str(d.info["warnings"]))
        Lmax = float(geo.get("L_ramp_max", 1e9))
        if d.L_ramp > Lmax:
            raise ValueError(f"L_ramp {d.L_ramp:.2f} > {Lmax}")
        return d

    def _eval_op(self, prob, tag: str, op: str, row: dict):
        """1 作動点を標準レシピで回し、ゲート不合格なら緩レシピ (step 倍・cfl 半分) で 1 回だけ再試行する (plan §4.13-1,2)。
        戻り値 (run_dir, rc, metrics, degraded)。**採用条件は rc == 0 かつ gates.verdict == PASS のみ** (旧 §4.13-3 の
        「rc != 0 でも力係数 STEADY なら採用」は撤回、R1)。degraded = 再試行 (緩レシピ) で通った場合 True。"""
        oc = self.optcfg
        # YAML の operating_points[].warm_from で「どの作動点の収束場から立ち上げるか」を指定する
        # (同一メッシュ・熱力学整合リマップ。NPR が遠い作動点には付けない — plan §5.1-1c)
        wf = next((o.get("warm_from") for o in self.ops if o["name"] == op), None)
        warm_src = self.dir / f"{tag}_{wf}" if wf else None
        if warm_src is not None and not (warm_src / "prepare_info.json").exists():
            print(f"   [{tag}/{op}] warm_from={wf} の run が無いので cold start", flush=True)
            warm_src = None
        ladder = [(dict(soft_steps=int(oc.get("soft_steps", 1500)), soft_cfl=float(oc.get("soft_cfl", 0.5)),
                        warm_lam_steps=int(oc.get("warm_lam_steps", 0)), warm_lam_cfl=float(oc.get("warm_lam_cfl", 0.2)),
                        mid_steps=int(oc.get("mid_steps", 0)), warm_src=warm_src,
                        warm_adapt_steps=int(oc.get("warm_adapt_steps", 500))), ""),
                   (dict(soft_steps=2 * int(oc.get("soft_steps", 1500)), soft_cfl=0.5 * float(oc.get("soft_cfl", 0.5)),
                         warm_lam_steps=2 * int(oc.get("warm_lam_steps", 0)), warm_lam_cfl=0.5 * float(oc.get("warm_lam_cfl", 0.2)),
                         mid_steps=2 * int(oc.get("mid_steps", 0)), warm_src=None), "_retry")]   # 再試行は cold start に落とす
        last = None
        for kw, suffix in ladder:
            rd = self.dir / f"{tag}_{op}{suffix}"
            if rd.exists():
                continue
            info = R.prepare(prob, rd, op=op)
            self._check_design(info)
            try:
                rc = R.run_staged(rd, "full", **kw)
            except RuntimeError as e:      # 起動段で落ちた: 壁出力が無いので次の梯子へ
                last = (rd, 1, {"gates": {"verdict": "FAIL", "fail_class": "DIVERGED", "reasons": [str(e)[:200]]}}, True)
                print(f"   [{tag}/{op}{suffix}] {e}", flush=True)
                continue
            out = R.collect(prob, rd, rc=rc, require_residual_pass=bool(oc.get("require_residual_pass", False)))
            g = out["gates"]
            if rc == 0 and g["verdict"] == "PASS":
                return rd, rc, out, suffix == "_retry"
            print(f"   [{tag}/{op}{suffix}] gate {g['fail_class']}: {'; '.join(g['reasons'])[:200]}", flush=True)
            last = (rd, rc, out, True)
        if last is None:                    # 梯子の全段が既存ディレクトリでスキップされた (resume 時)
            raise RuntimeError(f"{tag}/{op}: 梯子の全段が既存ディレクトリ。再開するなら該当 run を消すこと")
        return last

    def _check_design(self, info: dict) -> None:
        """設計段の物理的成立性 (INFEASIBLE)。L_ramp_max は probe (粗い kernel) でなく**最終輪郭**で再検査する (R6(d)/§5.1-8b)。"""
        if info["design"]["warnings"]:
            raise DesignInfeasible("design warnings: " + "; ".join(info["design"]["warnings"]))
        Lmax = float(self.base_raw["geometry"].get("L_ramp_max", 1e9))
        if float(info["design"]["L_ramp"]) > Lmax * (1 + float(self.optcfg.get("l_ramp_tol", 0.0))):
            raise DesignInfeasible(f"L_ramp {info['design']['L_ramp']:.3f} > L_ramp_max {Lmax} (final contour)")

    @staticmethod
    def _op_summary(out: dict, rd) -> dict:
        g = out.get("gates", {})
        # C_T 欄 = 目的量の値。目的量 (gates.objective) かその値が無ければ None (欠測)。圧力の C_T は C_T_p 欄だけに載せ、
        # 目的量の欄に代わりに入れない (SST で C_T_with_shear が欠けた行が圧力の C_T を目的量として表示していた。§5.1 R2g)
        obj_key = g.get("objective")
        return {"C_T": out.get(obj_key) if obj_key else None, "C_T_p": out.get("C_T"), "C_L": out.get("C_L"), "C_M": out.get("C_M"),
                "step": out.get("step"), "run_dir": str(rd), "sep_frac_ramp": out.get("sep_frac_ramp"), "sep_x_min_ramp": out.get("sep_x_min_ramp"),
                "forge_rc": out.get("forge_rc"), "gate": g.get("verdict"), "gate_fail_class": g.get("fail_class"),
                "slau_wall_normal_chi_effective": out.get("slau_wall_normal_chi_effective"), "flag_policy": out.get("flag_policy"),
                "scalar_gradient_effective": out.get("scalar_gradient_effective"),
                # 処置の実効値 (設計 DB の識別) と、その評価が読んだ w の格子署名・ハッシュ (来歴)
                ZTE_EFFECTIVE: out.get(ZTE_EFFECTIVE),
                "zero_thickness_edge_velocity_mesh_signature": (out.get("zero_thickness_edge_velocity") or {}).get("mesh_signature"),
                "zero_thickness_edge_velocity_field_hash": (out.get("zero_thickness_edge_velocity") or {}).get("field_hash"),
                "residual": g.get("residual", {}).get("verdict"), "objective": g.get("objective"),
                # 受理検査の版と残差の入力の検査済み要約 (学習の採否が現行の版と要約の合格を要求する、§5.1 R2k)。
                # **gates に版が無ければ None のまま** (旧い判定に現行の版を付けない)。要約も gates の残差の判定から作り、欠けた欄を補わない
                ACCEPTANCE_VERSION_KEY: g.get(ACCEPTANCE_VERSION_KEY),
                "residual_evidence": residual_evidence_summary(g.get("residual")),
                # 毎更新の EOS 床事象 (必須か・判定・区間・件数)。必須で不合格なら gate_fail_class が FLOOR_EVENT / FLOOR_UNVERIFIABLE
                "floor_events": floor_events_summary(g.get("floor_events")),
                # 判定区間の段の §6 #2・#3 の判定 (問題 YAML が evaluate.judgment_steps を持つ評価だけ。plan tooling-sern-te-wake-grid §5.1 R2)。
                # 学習の採否 (_judgment_op_ok) は要約 (窓の長さ・必須量の充足・温度の領域) まで照合する (§5.1 R2b)
                **({"judgment": (g.get("judgment") or {}).get("verdict"), "judgment_summary": judgment_summary(g.get("judgment"))}
                   if "judgment" in g else {}),
                "steadiness": {k: v.get("verdict") for k, v in g.get("steadiness", {}).get("series", {}).items()},
                # 評価方式の識別 (学習可否) と格子の来歴 (照合しない)。collect が出さない評価 (識別の導入前) はキーを作らない
                **({TE_WAKE_EFFECTIVE: out[TE_WAKE_EFFECTIVE]} if TE_WAKE_EFFECTIVE in out else {}),
                **({EVAL_METHOD_ID: out[EVAL_METHOD_ID], "eval_method": (out.get("eval_method") or {}).get("eval_method"),
                    "eval_method_reason": (out.get("eval_method") or {}).get("reason"),
                    # 評価方式から外した設計変数の値 (L_sw など。設計点として台帳に残す。plan tooling-sern-te-wake-grid §5.1 R3)
                    "design_values": (out.get("eval_method") or {}).get("design_values")} if EVAL_METHOD_ID in out else {}),
                **({"mesh_signature": ((out.get("mesh_provenance") or {}).get("now") or {}).get("mesh_signature"),
                    "dual_hash": ((out.get("mesh_provenance") or {}).get("now") or {}).get("dual_hash")} if "mesh_provenance" in out else {})}

    def _cm_check(self, row: dict) -> None:
        """C_M 制約: 加重平均の窓 (`cm_min/cm_max`) と作動点別の窓 (`cm_window: {op: [lo, hi]}`, R6(d))。物理的 INFEASIBLE。"""
        lo, hi = self.optcfg.get("cm_min"), self.optcfg.get("cm_max")
        if (lo is not None and row["C_M_w"] < lo) or (hi is not None and row["C_M_w"] > hi):
            raise DesignInfeasible(f"C_M_w {row['C_M_w']:.3f} outside [{lo}, {hi}]")
        for op, win in (self.optcfg.get("cm_window") or {}).items():
            cm = row["ops"].get(op, {}).get("C_M")
            if cm is None:
                continue
            wlo, whi = (win + [None, None])[:2] if isinstance(win, list) else (win.get("min"), win.get("max"))
            if (wlo is not None and cm < wlo) or (whi is not None and cm > whi):
                raise DesignInfeasible(f"C_M[{op}] {cm:.3f} outside [{wlo}, {whi}]")

    # -- 1 点評価 -----------------------------------------------------------------
    def evaluate(self, x, tag: str) -> dict:
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time(); prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "status": "FAIL", "fail_class": None, "ops": {}, "note": "", "degraded": False, "degraded_ops": [],
               "flag_policy": FLAG_POLICY, "zero_thickness_edge_velocity": self.zte_required,
               "te_wake": self.te_wake_required, "eval_method_required": self.eval_method_required}
        try:
            ct_w, L_ramp, cm_w, wsum = 0.0, None, 0.0, 0.0
            for o in self.ops:
                rd, rc, out, degraded = self._eval_op(prob, tag, o["name"], row)
                info = json.loads((rd / "prepare_info.json").read_text())
                L_ramp = info["design"]["L_ramp"]
                g = out.get("gates", {"verdict": "FAIL", "fail_class": "ERROR", "reasons": ["no gates"]})
                row["ops"][o["name"]] = self._op_summary(out, rd)
                if rc != 0 or g["verdict"] != "PASS":
                    raise EvalFailure(g.get("fail_class") or "DIVERGED", f"{o['name']}: rc={rc} {'; '.join(g.get('reasons', []))[:200]}")
                ct = out[g["objective"]]
                if ct is None or not np.isfinite(ct):
                    raise EvalFailure("DIVERGED", f"{o['name']}: objective {g['objective']}={ct}")
                if degraded:
                    row["degraded_ops"].append(o["name"])
                w = float(o.get("weight", 1.0)); wsum += w
                ct_w += w * ct; cm_w += w * out["C_M"]
                # 容量節約: 場は最終ステップのみ残す。**.xmf も一緒に消す** (h5 だけ消すと XDMF が消えた
                # h5 を指し続けて「中身の無い xmf」が残る — 2026-09-05 修正)
                vol = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
                for f_ in vol[:-1]:
                    f_.unlink()
                    f_.with_suffix(".xmf").unlink(missing_ok=True)
            row.update({"status": "PASS", "C_T_w": ct_w / wsum, "C_M_w": cm_w / wsum, "L_ramp": float(L_ramp),
                        "degraded": bool(row["degraded_ops"])})
            self._cm_check(row)
        except DesignInfeasible as e:       # 物理的に不成立 (学習対象外だが「失敗」ではない)
            row["note"] = str(e)[:300]
            row["status"] = "INFEASIBLE"
            row["fail_class"] = "CM_WINDOW" if "C_M" in str(e) else ("L_RAMP_MAX" if "L_ramp_max" in str(e) else "DESIGN")
        except EvalFailure as e:            # 数値失敗 (ゲート不合格)
            row["note"] = str(e)[:300]; row["status"] = "FAIL"; row["fail_class"] = e.fail_class
        except Exception as e:              # 想定外 (I/O 等)。INFEASIBLE と混同しない
            row["note"] = f"{type(e).__name__}: {str(e)[:280]}"; row["status"] = "FAIL"; row["fail_class"] = "ERROR"
        row["elapsed_s"] = time.time() - t0
        with open(self.ledger, "a") as fh:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
        self.rows.append(row)
        print(f"[{tag}] {row['status']}{'/' + row['fail_class'] if row['fail_class'] else ''}{' (degraded)' if row['degraded'] else ''} "
              f"x={np.round(x, 3).tolist()} C_T_w={row.get('C_T_w')} L={row.get('L_ramp')} C_M_w={row.get('C_M_w')} ({row['elapsed_s']:.0f}s) {row['note']}", flush=True)
        return row

    def _XF(self):
        # slauWallNormalChi の既定変更 (2026-09-26) 前後の評価を同じ応答関数として学習しない (plan
        # convection-slau-wall-normal-chi-default §4.4、codex plan M3)。flag_policy の無い旧行と不一致の行は除外する。
        # さらに mesh.scalarGradient の node 既定 lsq 化 (2026-09-27) 以降は、**全作動点の実効値が lsq と確認できた行だけ**を使う
        # (日付の一致だけでは gg 評価・不明が混ざる。codex diagnose 2026-09-27、plan gradient-scalar-lsq-unification #6)。
        ok = [r for r in self.rows if _learnable(r, *_zte_policy(self), method=_method_policy(self), floor=_floor_policy(self),
                                                 judgment=_judgment_policy(self))]
        X = np.array([r["x"] for r in ok]); F = np.array([[-r["C_T_w"], r["L_ramp"]] for r in ok])
        return X, F

    def run(self, n_doe: int, n_iter: int, batch: int) -> None:
        done = {r["tag"] for r in self.rows}
        cand = lhs(self.bounds, 3 * n_doe, seed=self.seed)
        doe = []
        for x in cand:
            if len(doe) >= n_doe:
                break
            if self.feasible(x):
                doe.append(x)
        print(f"DOE {len(doe)} 点", flush=True)
        for i, x in enumerate(doe):
            tag = f"doe_{i:03d}"
            if tag not in done:
                self.evaluate(x, tag)
        for it in range(n_iter):
            X, F = self._XF()
            if len(X) < 4:
                raise RuntimeError("PASS 評価が少なすぎる (<4)")
            print(f"--- iter {it}: PASS {len(X)} 点, HV={hypervolume2d(F, self.ref):.4f} ---", flush=True)
            krg = KrigingSet(self.bounds).fit(X, F)
            xs = propose_infill(_KrgBoth(krg), X, F, self.ref, self.bounds, n_infill=batch, seed=self.seed + 1000 + it, feasible=self.feasible)
            for j, x in enumerate(xs):
                tag = f"inf_{it:02d}_{j}"
                if tag not in done:
                    self.evaluate(x, tag)
        summary = self.summary()
        (self.dir / "pareto.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
        print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)

    def summary(self, rows=None) -> dict:
        """Pareto 要約。**degraded / tag / 作動点ごとのゲート要約を落とさない** (R1: pareto.json でも追える)。"""
        rows = self.rows if rows is None else rows
        # Pareto・HV の母集団も学習と同じ選別 (codex result 2026-09-27 M2: 学習から外した gg・旧方針の行が Pareto に混ざっていた)
        ok = [r for r in rows if _learnable(r, *_zte_policy(self), method=_method_policy(self), floor=_floor_policy(self),
                                            judgment=_judgment_policy(self))]
        X = np.array([r["x"] for r in ok]); F = np.array([[-r["C_T_w"], r["L_ramp"]] for r in ok])
        pareto = []
        if len(ok):
            for i in np.where(nondominated_mask(F))[0]:
                r = ok[i]
                pareto.append({"tag": r["tag"], "x": dict(zip(DV_ORDER, X[i].tolist())), "C_T_w": float(-F[i, 0]), "L_ramp": float(F[i, 1]),
                               "C_M_w": r["C_M_w"], "degraded": bool(r.get("degraded")), "degraded_ops": r.get("degraded_ops", []),
                               "flag_policy": r.get("flag_policy"),
                               "ops": {op: {**{k: v.get(k) for k in ("C_T", "C_M", "gate", "residual", "steadiness", "scalar_gradient_effective")},
                                            ZTE_EFFECTIVE: _zte_effective_of(v), TE_WAKE_EFFECTIVE: _te_wake_of(v),
                                            EVAL_METHOD_ID: v.get(EVAL_METHOD_ID), ACCEPTANCE_VERSION_KEY: v.get(ACCEPTANCE_VERSION_KEY)}
                                       for op, v in r.get("ops", {}).items()}})
        pareto.sort(key=lambda r: r["L_ramp"])
        classes = {}
        for r in rows:
            k = r["status"] + ("/" + r["fail_class"] if r.get("fail_class") else "")
            classes[k] = classes.get(k, 0) + 1
        n_status_pass = sum(1 for r in rows if r["status"] == "PASS")
        return {"n_eval": len(rows), "n_pass": int(len(ok)), "n_pass_excluded_by_policy": int(n_status_pass - len(ok)),
                "flag_policy": FLAG_POLICY, "required_scalar_gradient": REQUIRED_SCALAR_GRADIENT,
                "required_zero_thickness_edge_velocity": _zte_policy(self)[0], ZTE_CARRY_OVER_KEY: list(_zte_policy(self)[1]),
                "required_te_wake": _method_policy(self)[0], "required_eval_method_id": _method_policy(self)[1],
                "required_floor_events": _floor_policy(self), "required_judgment": _judgment_policy(self),
                "required_acceptance_version": ACCEPTANCE_VERSION,     # 一律 (キャンペーンの設定によらない、§5.1 R2k)
                EVAL_METHOD_CARRY_OVER_KEY: list(getattr(self, "eval_method_carry_over_records", ())),
                "n_degraded": int(sum(1 for r in ok if r.get("degraded"))),
                "hv": (hypervolume2d(F, self.ref) if len(ok) else 0.0), "ref": self.ref, "status_counts": classes,
                "gate_policy": "R1: rc==0 + finite field + residual no NaN/rising + objective & C_T/C_L/C_M STEADY (no divergent adoption)",
                "operating_points": self.ops, "pareto": pareto}

    # -- 既存キャンペーンの再判定 (CFD なし) ----------------------------------------
    def rejudge(self, out_dir) -> dict:
        """台帳の各評価を現行ゲートで判定し直す。元 run は読むだけ (metrics/force_history は out_dir/<tag>_<op>/ に書く)。
        rc は台帳に無いので run_case_stdout.log から復元する。判定し直した行は現行の受理検査の版と残差の入力の要約を持つ
        (pareto_rejudged の学習の母集団に入りうる。元の ledger.jsonl は書き換えない)。"""
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        new_rows = []; lines = ["| tag | 旧 status | 新 status | fail_class | degraded | C_T_w (旧→新) | 理由 |", "| --- | --- | --- | --- | --- | --- | --- |"]
        for r0 in self.rows:
            tag = r0["tag"]; prob = self.dir / f"{tag}.yaml"
            row = {"tag": tag, "x": r0["x"], "status": "FAIL", "fail_class": None, "ops": {}, "note": "", "degraded": False, "degraded_ops": [],
                   "old_status": r0["status"], "old_fail_class": r0.get("fail_class"), "old_C_T_w": r0.get("C_T_w"),
                   "flag_policy": r0.get("flag_policy"),   # 評価時の方針を引き継ぐ (再判定で現行方針に書き換えない)
                   "zero_thickness_edge_velocity": r0.get("zero_thickness_edge_velocity", "off"),   # 同上 (処置の導入前の行は off)
                   "te_wake": r0.get("te_wake", "off"), "eval_method_required": r0.get("eval_method_required")}   # 同上 (評価方式)
            try:
                if not prob.exists():
                    raise EvalFailure("ERROR", "problem yaml missing")
                ct_w, cm_w, wsum, L_ramp = 0.0, 0.0, 0.0, None
                for o in self.ops:
                    rd_std, rd_retry = self.dir / f"{tag}_{o['name']}", self.dir / f"{tag}_{o['name']}_retry"
                    cands = [d for d in (rd_std, rd_retry) if (d / "prepare_info.json").exists()]
                    if not cands:
                        raise EvalFailure("NO_FORCES", f"{o['name']}: no run dir")
                    chosen = None
                    for rd in cands:       # 標準 → 再試行の順に、現行ゲートを通る最初の run を採る
                        out = R.collect(prob, rd, out_dir=out_dir / rd.name, require_residual_pass=bool(self.optcfg.get("require_residual_pass", False)))
                        row["ops"][o["name"]] = self._op_summary(out, rd)
                        if out["forge_rc"] == 0 and out["gates"]["verdict"] == "PASS":
                            chosen = (rd, out); break
                    if chosen is None:
                        g = out["gates"]
                        raise EvalFailure(g.get("fail_class") or "DIVERGED", f"{o['name']}: rc={out['forge_rc']} {'; '.join(g['reasons'])[:200]}")
                    rd, out = chosen
                    info = json.loads((rd / "prepare_info.json").read_text()); self._check_design(info); L_ramp = info["design"]["L_ramp"]
                    if rd.name.endswith("_retry"):
                        row["degraded_ops"].append(o["name"])
                    w = float(o.get("weight", 1.0)); wsum += w; ct_w += w * out[out["gates"]["objective"]]; cm_w += w * out["C_M"]
                row.update({"status": "PASS", "C_T_w": ct_w / wsum, "C_M_w": cm_w / wsum, "L_ramp": float(L_ramp), "degraded": bool(row["degraded_ops"])})
                self._cm_check(row)
            except DesignInfeasible as e:
                row["note"] = str(e)[:300]; row["status"] = "INFEASIBLE"
                row["fail_class"] = "CM_WINDOW" if "C_M" in str(e) else ("L_RAMP_MAX" if "L_ramp_max" in str(e) else "DESIGN")
            except EvalFailure as e:
                row["note"] = str(e)[:300]; row["status"] = "FAIL"; row["fail_class"] = e.fail_class
            except Exception as e:
                row["note"] = f"{type(e).__name__}: {str(e)[:280]}"; row["status"] = "FAIL"; row["fail_class"] = "ERROR"
            new_rows.append(row)
            ctn = row.get("C_T_w"); cto = r0.get("C_T_w")
            lines.append(f"| {tag} | {r0['status']}{'/' + r0['fail_class'] if r0.get('fail_class') else ''} | {row['status']} | {row['fail_class'] or ''} | "
                         f"{'yes' if row['degraded'] else ''} | {'' if cto is None else f'{cto:.4f}'} → {'' if ctn is None else f'{ctn:.4f}'} | {row['note'][:120]} |")
        with open(out_dir / "ledger_rejudged.jsonl", "w") as fh:
            for r in new_rows:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        summ = self.summary(new_rows); summ["source_campaign"] = str(self.dir)
        old_ok = [r for r in self.rows if r["status"] == "PASS"]
        summ["old"] = {"n_pass": len(old_ok), "hv": (hypervolume2d(np.array([[-r["C_T_w"], r["L_ramp"]] for r in old_ok]), self.ref) if old_ok else 0.0)}
        (out_dir / "pareto_rejudged.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False))
        md = [f"# 再判定 (R1 ゲート): {self.dir}", "", f"旧: PASS {summ['old']['n_pass']} / HV {summ['old']['hv']:.4f} → "
              f"新: PASS {summ['n_pass']} (degraded {summ['n_degraded']}) / HV {summ['hv']:.4f}", "", f"status 内訳: {summ['status_counts']}", ""] + lines
        (out_dir / "rejudge_summary.md").write_text("\n".join(md) + "\n")
        print("\n".join(md), flush=True)
        return summ


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="⑤ SERN 多作動点 MOO キャンペーン")
    ap.add_argument("problem"); ap.add_argument("campaign_dir")
    ap.add_argument("--n-doe", type=int, default=12); ap.add_argument("--n-iter", type=int, default=3)
    ap.add_argument("--batch", type=int, default=2); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref", type=float, nargs=2, default=(-0.90, 20.0))
    ap.add_argument("--rejudge", default=None, metavar="OUT_DIR", help="CFD を回さず既存キャンペーンを現行ゲートで再判定して OUT_DIR に書く")
    a = ap.parse_args(argv)
    c = SernCampaign(a.problem, a.campaign_dir, ref=tuple(a.ref), seed=a.seed)
    if a.rejudge:
        c.rejudge(a.rejudge); return 0
    c.run(a.n_doe, a.n_iter, a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
