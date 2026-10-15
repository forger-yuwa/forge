#!/usr/bin/env python3
r"""仮想の停止の後の同じ軌道の検証 — 停止規則 v2 が止めたはずの格子点 t_s から span の間、目的の量が t_s の値から許容の内に
留まり、各格子点の判定が SETTLED のままかを、停止規則なしで回した run で確かめる
(plan `plans/active/tooling-settle-convergence-and-monitor.md` §6「V3(b)・V4 の実施」の (1) `run_0012_ns_free`・(1′)・(1″)
(`wt_nozzle_v2` の再生は `--criteria`); codex diagnose `notes/reviews/2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md`
の Major 1)。保持 0 (`wt_nozzle_v2`) では停止候補 t_s は最初に SETTLED になった格子点 (hold_step が保持 0 で即到達)。

判定の部品は `settle.py` (`RunData`・`judge`・量を作る関数) をそのまま使う。design には依存しない (停止規則 v2 の保持の数え方は
`design/forge_design/evaluate/stop_rule.py` から写し、出典の行を書いた; 同じ判定になることは `test_settle_virtual_stop.py` の (vi) が
本物の v2 の見張りと偽の forge で確かめる)。

## 手順

1. 再生: 最初に判定できる格子点 first = 開始 + 2W + 間隔 (stop_rule.py:1181) から入力のデータの最後まで、間隔ごとに
   `settle.judge(run, s)` (RunData は 1 つを使い回す)。
2. 停止候補 t_s: 停止規則 v2 の見張り (`watch_main_stage_v2`) と同じ規則で、最初に STOP_REQUEST を置く格子点を求める。
   - SETTLED の格子点で保持を進め、条件を満たし始めた格子点 s₀ から s − s₀ ≥ hold になった格子点で要求 (kind settled)。
     SETTLED でない格子点で保持をやり直す (`breaks`)。hold は settle_criteria.yaml の hold、無ければ 40,000。
   - SETTLED でない格子点が DIVERGED、または INDETERMINATE で理由の分類に「未設定」(unset) を含むとき、v2 はその格子点で要求
     (kind diverged / criteria_unset) を置いて止める (FAILED) → 停止候補は無い。
   - 出典: 保持は `_hold_step` (stop_rule.py:624–637) を `hold_step` に写した。要求の分岐は watch_main_stage_v2 の判定の繰り返し
     (stop_rule.py:1558–1588)、理由の分類と why は `judge_summary` (stop_rule.py:1297–1310)、hold の既定は stop_rule.py:110・1124。
   - 上限の格子点 (§5.1 #5c (A)、stop_rule.watch_main_stage_v2 と同じ): 保持を満たした格子点 s が本段の上限 cap (本段の nStepOuter)
     以上なら settled の要求を出さない (上限の後に step が無く forge は要求を見ないので v2 は CENSORED) → 停止候補なし。DIVERGED・
     未設定の要求は cap に依らず今どおり。cap は `main_stage_cap` (MAIN_STAGE_STARTED の main_config.nStepOuter、無ければ
     prepare_info.json の stop_rule.cap) から。どちらも無ければ cap なしで再生し、その旨を note に書く。
   - 再生では forge が cap まで走っているとみなす (v2 は forge が上限より前に終わった後の格子点でも要求を置かない (§5.1 #5c (B)) が、
     これは時刻の印なので再生には持ち込まない。cap が分からず最後の格子点で保持を満たしたときは note に書く)。
3. t_s から t_s + span (既定 80,000) の格子点 (両端を含む) で:
   (a) 目的の量の t_s の値からの差: 出口コア M は |M − M(t_s)| ≤ tol、線は max_k |d_k − d_k(t_s)| ≤ tol_pos (位置ごと)、
       出口の流れ角は max_i |θ_i − θ_i(t_s)| ≤ tol (節点ごと)、壁は |Q_i − Q_i(t_s)| ≤ q_tol·A_i + eps_rel·|Q_i(t_s)| (節点ごと;
       A_i は t_s の壁の h5 から `settle.wall_node_area`、0 の許容は settle と同じく差が 0 のときだけ内)。許容は
       settle_criteria.yaml の各量の閾値そのもの。値は settle.py の関数 (`exit_core_mach_of`・`eta_line_nodes` + `line_quantities`・
       `exit_flow_angle_of`・`RunData.wall`) で作る (演算の順も同じ; 線の軸の節点の重みの検査も同じ)。線の maxabs・halfamp の差は
       記録だけ (判定に使わない)。
   (b) 各格子点の全体の verdict が SETTLED。
4. 判定 (この順に決める):
   - 未成立: 基準が無い・誤り、停止候補が無い、t_s > ts_max (既定 120,000)。
   - B: span の格子点で (a) の超過 (節点の値の非有限を含む) か (b) の NOT_SETTLED・DIVERGED が 1 つでもある。どの量・位置・step
     かを書く。**データの欠け・判定不能が同時にあっても B にする** (超過はそれ自体が反例; 欠けは `incomplete` に並べる。
     2026-10-16 に結果を見る前に決めた順; plan に無いので呼び出し側の確認待ち)。
   - 未成立: 超過は無いが、span の終わりまでデータが無い、(b) の INDETERMINATE、(a) の値を作れない (データ不足・格子と線の不整合・
     series_csv の入力・壁の値の非有限) のどれかがある。
   - A: 上のどれでもない (span の全期間で (a)・(b) とも内)。span は有限の検証範囲で、永久の整定の保証ではない。
5. 追加の記録: 目的の量 (exit_core_mach・line_mach・exit_flow_angle・wall_heat_nodes) がすべて SETTLED なのに残差か流量の収支が
   SETTLED でない格子点の数と step (NOT_SETTLED・DIVERGED と INDETERMINATE を分け、t_s の前・span・その後に分ける)。
   残差の条件が独立の情報を足すかの確認の材料 (AND に入っているだけでは足したと言わない)。
6. 参考の記録 (A/B の判定には使わない; 2026-10-16 主セッションの指示): 残差と流量の収支の条件を外し、目的の量の verdict の AND
   (優先は settle と同じ) に手順 2 と同じ保持の規則・要求の分岐を当てた停止候補 t_s′。t_s′ から span の (a) は手順 3 と同じ差、
   (b) は全体でなく目的の量の verdict。t_s との関係 (t_s′ なし / t_s′ = t_s / 何 step 早い・遅い) と、t_s′ からの検証の結果
   (超過なし / 超過あり / 判定できない; 順は手順 4 と同じ) を書く。t_s′ が早く、そこからの検証が「超過あり」で本検証が A なら、
   残差・流量の条件が早すぎる停止を防いだことの材料になる。

## 出力 (`<run>/settle_virtual_stop/`、--out で変える)

- `virtual_stop.json`: 結果 (基準の sha256・道具と settle.py の sha256・t_s・保持の経過・span の量ごとの最大の比・B の根拠・欠け・
  追加の記録・所要・使ったファイルの大きさと sha256)。手順 6 の記録は `target_only` (`used_for_verdict: false`)。
- `VERDICT.txt`: 1 行目 `VERDICT: A` / `VERDICT: B` / `VERDICT: 未成立`、続けて理由と、「参考 (判定に使わない)」の t_s′ の行。
  B のときは「B の内訳」の行で (a) 目的の量の累積の超過と (b) 全条件の持続の違反 (崩れた条件 = 量ごとの格子点の数・verdict) を分けて
  書く (JSON の `b_breakdown`; codex diagnose 2026-10-16 settle-v2-span10k — (b) だけの B は早すぎる停止の証明ではない)。
- `grid_verdicts.csv`: 格子点ごとの全体と量ごとの verdict・v2 の保持の状態と事象・区間・追加の記録の印、目的の量だけの verdict と
  その保持の状態・事象・t_s′ からの区間。
- `grid_checks.csv`: 格子点ごと・量ごとの settle の検査の値と許容 (統計量の時系列) と、検査の無い理由。
- `span_deviation.csv`: 格子点ごと・量ごとの t_s からの差・許容・比・位置 (t_s の前と span の後も参考に出す; 判定は region が span の行だけ)。
- `span_deviation_ts_prime.csv`: 同じ形で t_s′ からの差 (参考; t_s′ が無ければ見出しだけ)。
- `deviation_from_ts.png` (量ごとの t_s からの差と許容の線)、`grid_verdicts.png` (格子点ごとの量ごとの verdict)。t_s′ が t_s と
  違えば紫の破線で示す。

使い方:
  python3 settle_virtual_stop.py <run_dir> [--span 80000] [--ts-max 120000] [--out DIR] [--criteria FILE] [--no-plot] [--json]
  --criteria FILE: run の settle_criteria.yaml の代わりに FILE を基準にする (§6 (1″) の `wt_nozzle_v2` の再生など; run の基準は
    読まず書き換えない。基準の中の相対のパスは run から解く)。run の基準の結果を上書きしないよう --out が要り、既定の
    <run>/settle_virtual_stop は使えない。結果の criteria.path・files は FILE を指す。
終了コード: 0 A / 1 B / 3 未成立 (基準・引数の誤りも 3)。
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import settle  # noqa: E402

SCHEMA = "forge.settle_virtual_stop/v1"
PLAN = ("plans/active/tooling-settle-convergence-and-monitor.md §6「V3(b)・V4 の実施」(1) run_0012_ns_free "
        "(codex diagnose 2026-10-16-settle-v3b-v4-preregistration の Major 1)")
OUT_DIR = "settle_virtual_stop"
SPAN_DEFAULT = 80_000
TS_MAX_DEFAULT = 120_000
# = stop_rule.HOLD_STEPS_DEFAULT (stop_rule.py:110)。settle_criteria.yaml に hold が無いときの v2 の保持 (stop_rule.py:1124)
HOLD_STEPS_DEFAULT = 40_000
# 日本語のラベルのフォント (手元・AWS とも ~/.fonts に置いてある; 無ければ既定のフォントで描く)
FONT = str(Path.home() / ".fonts" / "NotoSansCJKjp-Regular.otf")
RESULT_A, RESULT_B, RESULT_NA = "A", "B", "未成立"
EXIT_CODE = {RESULT_A: 0, RESULT_B: 1, RESULT_NA: 3}
REQ_SETTLED, REQ_DIVERGED, REQ_UNSET = "settled", "diverged", "criteria_unset"
REGION_BEFORE, REGION_SPAN, REGION_AFTER, REGION_NO_TS = "before_ts", "span", "after_span", "no_ts"
FAIL_VERDICTS = ("NOT_SETTLED", "DIVERGED")
UNIT = {"exit_core_mach": "[-]", "line_mach": "[%pt]", "exit_flow_angle": "[deg]", "wall_heat_nodes": "[比]"}
LABEL = {"exit_core_mach": "出口コア M", "line_mach": "線の各位置の d", "exit_flow_angle": "出口の流れ角 (節点ごと)",
         "wall_heat_nodes": "壁の局所の熱流 (節点ごと)"}


class _GridMismatch(Exception):
    """線の補間に軸の節点が入る (settle の CAT_GRID と同じ条件)。"""


# ====================================================================== 停止規則 v2 の保持 (写し)
def hold_step(hs: dict, ok: bool, step: int, hold: int, why=None) -> bool:
    """停止規則 v2 の保持の 1 格子点ぶんの更新。`design/forge_design/evaluate/stop_rule.py:624–637` の `_hold_step` の写し
    (design に依存しないため; 同じ判定になることは test_settle_virtual_stop.py の (vi) が元の関数と比べる)。
    hs = {start, reached, breaks}。条件を満たした最初の格子点 s₀ から数え、step − s₀ ≥ hold になった格子点で True。途中で条件が
    外れたら到達を取り消して breaks に記録し、次に満たした格子点から数え直す。hold 0 は条件を満たした格子点で True。"""
    if not ok:
        if hs.get("start") is not None:
            hs["breaks"].append({"start": hs["start"], "broken_at": int(step), "why": why})
            hs["start"] = None
        return False
    if hs.get("start") is None:
        hs["start"] = int(step)
        hs["reached"].append(int(step))
    return int(step) - hs["start"] >= int(hold)


def summarize(result: dict, kinds=None) -> dict:
    """`settle.judge` の結果の要約。categories と why は `judge_summary` (stop_rule.py:1297–1310) と同じ作り方 (digest は作らない)。
    kinds を渡すとその種類の量だけで要約する (verdict は settle と同じ優先 DIVERGED > INDETERMINATE > NOT_SETTLED > SETTLED の AND;
    t_s′ の記録用)。"""
    qs = result.get("quantities") or []
    whole = kinds is None
    if not whole:
        qs = [q for q in qs if q["kind"] in kinds]
    cats = {c for q in qs for c in q.get("categories", [])}
    if whole and result.get("reason_category"):
        cats.add(result["reason_category"])
    out = {"step": result.get("terminal_step"), "end_step": result.get("end_step"),
           "verdict": result["verdict"] if whole else settle._overall(qs),
           "by_quantity": {q["name"]: q["verdict"] for q in qs}, "kinds": {q["name"]: q["kind"] for q in qs},
           "categories": sorted(cats)}
    why = ([str(result["reason"])] if whole and result.get("reason") else []) + [
        f"{q['name']}: [{r['label']}] {r['text']}" for q in qs if q["verdict"] != "SETTLED" for r in q.get("reasons", [])[:1]]
    if why:
        out["why"] = why[:8]
    return out


def emulate_v2(rows: list, hold: int, cap: int | None = None) -> tuple[dict | None, dict]:
    """格子点の順の要約 rows (first から) に、停止規則 v2 の見張りの判定の分岐 (stop_rule.py の watch_main_stage_v2 の判定の繰り返し)
    を当てて、最初の STOP_REQUEST を返す。戻り (request, hs)。request = {kind, step, settled_since, reason} か None。v2 と同じく要求の
    後も保持の状態は更新し続ける。各 row に v2 の事象 `v2_event` と、その格子点の後の保持の始まり `hold_start` を書き足す。
    cap (本段の上限) を渡すと、保持を満たした格子点が cap 以上なら settled の要求を出さず (v2 は CENSORED; §5.1 #5c (A))、その格子点を
    hs["reached_at_cap"] に書く (無ければ None)。DIVERGED・未設定の要求は cap に依らない。"""
    hs = {"hold_steps": int(hold), "start": None, "reached": [], "breaks": [], "reached_at_cap": None}
    req = None
    for row in rows:
        s, v = int(row["step"]), row["verdict"]
        ev = []
        if v == "SETTLED":
            new = hs["start"] is None
            if new:
                ev.append("hold_start")
            if hold_step(hs, True, s, hold) and req is None:
                if cap is not None and s >= cap:
                    hs["reached_at_cap"] = s          # 上限の後に step が無く要求は果たされない (v2 は CENSORED)
                    ev.append("hold_at_cap:no_request")
                else:
                    req = {"kind": REQ_SETTLED, "step": s, "settled_since": hs["start"],
                           "reason": f"3 条件が step {hs['start']} から {s - hs['start']} step 続いた (保持 {hold})"}
                    ev.append("request:settled")
        else:
            had = hs["start"]
            hold_step(hs, False, s, hold, "; ".join(row.get("why") or [v])[:400])
            if had is not None:
                ev.append("hold_break")
            if req is None:
                if v == "DIVERGED":
                    req = {"kind": REQ_DIVERGED, "step": s, "settled_since": None, "reason": f"判定が DIVERGED (step {s})"}
                    ev.append("request:diverged")
                elif v == "INDETERMINATE" and settle.CAT_UNSET in row["categories"]:
                    req = {"kind": REQ_UNSET, "step": s, "settled_since": None,
                           "reason": f"判定が INDETERMINATE (未設定; 基準と run の出力が合わない) (step {s})"}
                    ev.append("request:criteria_unset")
        row["hold_start"] = hs["start"]
        row["v2_event"] = ",".join(ev)
    return req, hs


def main_stage_cap(run) -> tuple[int | None, str]:
    """本段の上限 cap (= 本段の nStepOuter; 停止規則 v2 の見張りの上限) と出所。`<run>/MAIN_STAGE_STARTED` (run_staged_ns が本段の forge の
    直前に書く) の main_config.nStepOuter、無ければ `<run>/prepare_info.json` の stop_rule.cap。どちらも無い・読めない・正の整数でない
    → (None, 理由)。"""
    run = Path(run)
    why = []
    for name, get in (("MAIN_STAGE_STARTED", lambda d: (d.get("main_config") or {}).get("nStepOuter")),
                      ("prepare_info.json", lambda d: (d.get("stop_rule") or {}).get("cap"))):
        p = run / name
        if not p.is_file():
            why.append(f"{name} が無い")
            continue
        try:
            v = get(json.loads(p.read_text()))
        except Exception as e:  # noqa: BLE001
            why.append(f"{name} を読めない ({type(e).__name__})")
            continue
        if isinstance(v, int) and not isinstance(v, bool) and v > 0:
            return v, (f"{name} の main_config.nStepOuter" if name == "MAIN_STAGE_STARTED" else f"{name} の stop_rule.cap")
        why.append(f"{name} に正の整数の値が無い ({v!r})")
    return None, "; ".join(why)


def first_judgeable(crit: dict) -> int:
    """v2 が最初に判定する格子点 = 開始 + 2W + 間隔 (stop_rule.py:1181)。"""
    return crit["grid"]["start_step"] + 2 * crit["window"]["W"] + crit["grid"]["out_step_interval"]


def v2_hold(crit: dict) -> tuple[int, str]:
    """v2 の保持 (stop_rule.py:1124)。"""
    if "hold" in crit:
        return int(crit["hold"]), "settle_criteria.yaml の hold"
    return HOLD_STEPS_DEFAULT, f"既定 {HOLD_STEPS_DEFAULT} (stop_rule.HOLD_STEPS_DEFAULT)"


# ====================================================================== 目的の量の値 (settle.py の関数で作る)
def target_values(data, crit: dict, q: dict, s: int) -> dict:
    """格子点 s の目的の量の値。settle.py の関数と同じ演算の順で作る (settle._eval_* と同じ呼び方)。
    読めない・揃わない → settle._DataMissing、節点の値の非有限 → settle._NonFinite、線に軸の節点が入る → _GridMismatch。"""
    kind = q["kind"]
    if "series_csv" in q:
        raise settle._DataMissing(f"{q['name']}: series_csv (再生用の系列) の入力は位置ごと・節点ごとの値を持たない — t_s からの差を作れない")
    ctx = {"data": data}
    if kind in settle.NODE_KINDS:
        st = crit["structured"]
        obs = settle._node_obs(ctx, q, s, None)
        if kind == "exit_core_mach":
            return {"value": settle.exit_core_mach_of(obs, st["ni"], st["nj"], q["eta_range"])}
        if kind == "exit_flow_angle":
            return {"theta": settle.exit_flow_angle_of(obs, st["ni"], st["nj"], q["eta_range"])}
        S = st["scale_m"]
        xq = np.linspace(q["x_range_m"][0] / S, q["x_range_m"][1] / S, q["n"])     # settle._eval_line と同じ
        dg = {}
        lq = settle.line_quantities(settle.eta_line_nodes(obs, "M", q["eta"], xq, st["ni"], st["nj"], diag=dg), q["Md"])
        wax = dg.get("axis_weight_max", 0.0)
        if wax > settle.AXIS_WEIGHT_MAX:
            raise _GridMismatch(f"線 η {q['eta']} の補間に軸の節点が重み {wax:.3g} で入る (> {settle.AXIS_WEIGHT_MAX:g})")
        return {"d": lq["d"], "maxabs": lq["maxabs"], "halfamp": lq["halfamp"], "x_m": xq * S}
    if kind == "wall_heat_nodes":
        walls = []
        for wl in q["walls"]:
            ob = data.wall(q["file_pattern"], q["dataset"], wl, s)
            if not np.all(np.isfinite(ob["Q"])):      # settle._eval_wall と同じく判定不能 (適用外の NaN と発散を区別できない)
                nbad = int(np.count_nonzero(~np.isfinite(ob["Q"])))
                raise settle._DataMissing(f"{Path(ob['path']).name}: {q['dataset']} に非有限 {nbad}/{len(ob['Q'])} 点")
            walls.append({"name": wl["name"], "physID": wl["physID"], "Q": ob["Q"], "coord": ob["coord"], "conne": ob["conne"]})
        return {"walls": walls}
    raise ValueError(f"目的の量でない種類 {kind}")


def _ratio(diff: float, tol: float) -> float:
    if tol > 0:
        return float(diff / tol)
    return math.inf if diff > 0 else 0.0


def deviation(q: dict, ref: dict, cur: dict, area_cache: dict) -> dict:
    """t_s の値 ref からの差 (§手順 3 (a))。戻り {diff, tol, ratio, exceeded, where, ...}。exceeded は settle と同じ `差 <= 許容` の否定。"""
    th = q["thresholds"]
    kind = q["kind"]
    if kind == "exit_core_mach":
        diff = abs(float(cur["value"]) - float(ref["value"]))
        tol = th["tol"]
        return {"diff": diff, "tol": tol, "ratio": _ratio(diff, tol), "exceeded": not diff <= tol, "where": "出口コアの平均",
                "signed": float(cur["value"]) - float(ref["value"]), "value": float(cur["value"]), "ref": float(ref["value"])}
    if kind == "line_mach":
        dd = np.abs(cur["d"] - ref["d"])
        k = int(np.argmax(dd))
        tol = th["tol_pos"]
        return {"diff": float(dd[k]), "tol": tol, "ratio": _ratio(float(dd[k]), tol), "exceeded": not float(dd[k]) <= tol,
                "where": f"xq の位置 {k} (x {float(ref['x_m'][k]):.6g} m)", "index": k, "x_m": float(ref["x_m"][k]),
                "signed": float(cur["d"][k] - ref["d"][k]), "value": float(cur["d"][k]), "ref": float(ref["d"][k]),
                "maxabs_diff": float(cur["maxabs"] - ref["maxabs"]), "halfamp_diff": float(cur["halfamp"] - ref["halfamp"])}
    if kind == "exit_flow_angle":
        if cur["theta"].shape != ref["theta"].shape:
            raise settle._DataMissing(f"{q['name']}: 出口コアの節点の数が t_s と違う ({cur['theta'].shape} と {ref['theta'].shape})")
        dd = np.abs(cur["theta"] - ref["theta"])
        k = int(np.argmax(dd))
        tol = th["tol"]
        return {"diff": float(dd[k]), "tol": tol, "ratio": _ratio(float(dd[k]), tol), "exceeded": not float(dd[k]) <= tol,
                "where": f"コアの {k} 番目の節点", "index": k, "signed": float(cur["theta"][k] - ref["theta"][k]),
                "value": float(cur["theta"][k]), "ref": float(ref["theta"][k])}
    if kind == "wall_heat_nodes":
        worst = None
        for wr, wc in zip(ref["walls"], cur["walls"]):
            if not np.array_equal(wc["coord"], wr["coord"]):
                raise settle._DataMissing(f"{q['name']}: 壁 {wr['name']} の節点の COORD が t_s と違う (節点の順を対応させられない)")
            key = (q["name"], wr["name"], wr["physID"])
            if key not in area_cache:
                area_cache[key] = settle.wall_node_area(wr["coord"], wr["conne"], q["axisymmetric"])
            A, _how = area_cache[key]
            tol = th["q_tol"] * A + th["eps_rel"] * np.abs(wr["Q"])
            dQ = np.abs(wc["Q"] - wr["Q"])
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(tol > 0, dQ / tol, np.where(dQ > 0, np.inf, 0.0))   # settle._eval_wall と同じ
            k = int(np.argmax(ratio))
            cand = {"diff": float(dQ[k]), "tol": float(tol[k]), "ratio": float(ratio[k]), "exceeded": not float(ratio[k]) <= 1.0,
                    "where": f"壁 {wr['name']} 節点 {k} (x {float(wr['coord'][k, 0]):.6g} m, r {float(wr['coord'][k, 1]):.6g} m)",
                    "wall": wr["name"], "index": k, "x_m": float(wr["coord"][k, 0]), "r_m": float(wr["coord"][k, 1]),
                    "signed": float(wc["Q"][k] - wr["Q"][k]), "value": float(wc["Q"][k]), "ref": float(wr["Q"][k])}
            if worst is None or cand["ratio"] > worst["ratio"]:
                worst = cand
        return worst
    raise ValueError(kind)


def _tol_text(q: dict) -> str:
    th = q["thresholds"]
    k = q["kind"]
    if k == "exit_core_mach":
        return f"|M − M(t_s)| ≤ tol {th['tol']:.6g}"
    if k == "line_mach":
        return f"max_k |d_k − d_k(t_s)| ≤ tol_pos {th['tol_pos']:.6g} %pt"
    if k == "exit_flow_angle":
        return f"max_i |θ_i − θ_i(t_s)| ≤ tol {th['tol']:.6g} deg"
    return f"|Q_i − Q_i(t_s)| ≤ q_tol {th['q_tol']:.6g}·A_i + eps_rel {th['eps_rel']:.6g}·|Q_i(t_s)|"


# ====================================================================== 本体
def _sha256_file(p: Path) -> str | None:
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _region(s: int, ts, span: int) -> str:
    if ts is None:
        return REGION_NO_TS
    if s < ts:
        return REGION_BEFORE
    return REGION_SPAN if s <= ts + span else REGION_AFTER


def _span_check(ts, span: int, targets: list, vals: dict, grid: list, by_step: dict, last, verdict_key: str, why_key: str,
                b_label: str, ts_label: str, b_kinds=None) -> dict:
    """停止候補 ts から ts + span の格子点の (a)・(b) (§手順 3)。値 vals は全格子点で作ったものを使い回す。
    verdict_key・why_key は (b) に使う row の verdict と理由のキー (本検証は全体の "verdict"・"why"、t_s′ の記録は目的の量だけの
    "target_verdict"・"target_why")。b_kinds は (b) の違反の内訳に数える量の種類 (None = すべて)。
    戻り {a, b_not_settled, evidence, incomplete, dev_rows, series, grid_points}。b_not_settled の各要素に、SETTLED でなかった量と
    その verdict (`quantities`) を付ける (B の内訳)。
    evidence は超過 (B の根拠)、incomplete はデータの欠け・判定不能。ts が None なら量ごとに「{ts_label} なし」だけを書く。"""
    out = {"a": {}, "b_not_settled": [], "evidence": [], "incomplete": [], "dev_rows": [],
           "series": {q["name"]: {"kind": q["kind"], "steps": [], "signed": [], "diff": [], "tol": [], "ratio": []} for q in targets},
           "grid_points": []}
    span_end = ts + span if ts is not None else None
    span_steps = [s for s in grid if ts is not None and ts <= s <= span_end]
    out["grid_points"] = span_steps
    area_cache = {}
    for q in targets:
        nm = q["name"]
        summ = {"kind": q["kind"], "tolerance": _tol_text(q), "n_points": 0, "max_ratio": None, "at_step": None, "where": None,
                "diff": None, "tol": None, "n_exceeded": 0, "first_exceeded": None, "exceeded": [], "status": None}
        out["a"][nm] = summ
        if ts is None:
            summ["status"] = f"{ts_label} なし"
            continue
        ref = vals[nm].get(ts)
        if ref is None or not ref["ok"]:
            txt = ref["text"] if ref else f"{ts_label} の格子点のデータが無い"
            summ["status"] = f"{ts_label} の値を作れない"
            out["incomplete"].append(f"(a) {nm}: {ts_label} {ts} の値を作れない — {txt}")
            continue
        for s in grid:
            cur = vals[nm][s]
            reg = _region(s, ts, span)
            dr = {"step": s, "region": reg, "quantity": nm, "kind": q["kind"]}
            if not cur["ok"]:
                dr.update(status=cur["category"], detail=cur["text"])
                out["dev_rows"].append(dr)
                if reg == REGION_SPAN:
                    if cur["category"] == settle.CAT_NONFINITE:
                        out["evidence"].append(f"(a) {nm}: step {s} の節点の値が非有限 — {cur['text']}")
                        summ["n_exceeded"] += 1
                        summ["first_exceeded"] = summ["first_exceeded"] or s
                    else:
                        out["incomplete"].append(f"(a) {nm}: step {s} の値を作れない [{settle.CAT_LABEL[cur['category']]}] — {cur['text']}")
                continue
            try:
                d = deviation(q, ref["v"], cur["v"], area_cache)
            except settle._DataMissing as e:
                dr.update(status=settle.CAT_DATA, detail=str(e))
                out["dev_rows"].append(dr)
                if reg == REGION_SPAN:
                    out["incomplete"].append(f"(a) {nm}: step {s} の差を作れない — {e}")
                continue
            dr.update(status="ok", **{k: d.get(k) for k in ("diff", "tol", "ratio", "exceeded", "where", "signed", "value", "ref",
                                                             "x_m", "r_m", "maxabs_diff", "halfamp_diff")})
            out["dev_rows"].append(dr)
            sr = out["series"][nm]
            sr["steps"].append(s)
            sr["signed"].append(d["signed"])
            sr["diff"].append(d["diff"])
            sr["tol"].append(d["tol"])
            sr["ratio"].append(d["ratio"])
            if reg != REGION_SPAN:
                continue
            summ["n_points"] += 1
            if summ["max_ratio"] is None or d["ratio"] > summ["max_ratio"]:
                summ.update(max_ratio=d["ratio"], at_step=s, where=d["where"], diff=d["diff"], tol=d["tol"])
            if d["exceeded"]:
                summ["n_exceeded"] += 1
                if summ["first_exceeded"] is None:
                    summ["first_exceeded"] = s
                if len(summ["exceeded"]) < 20:
                    summ["exceeded"].append({"step": s, "ratio": d["ratio"], "diff": d["diff"], "tol": d["tol"], "where": d["where"]})
        if summ["status"] is None:
            summ["status"] = "超過" if summ["n_exceeded"] else "内"
        if summ["n_exceeded"] and summ["max_ratio"] is not None:
            out["evidence"].append(
                f"(a) {nm} ({LABEL.get(q['kind'], q['kind'])}): span の {summ['n_exceeded']} 格子点で超過 (最初 step {summ['first_exceeded']})、"
                f"最大は step {summ['at_step']} の {summ['where']} で差 {summ['diff']:.6g} > 許容 {summ['tol']:.6g} (比 {summ['max_ratio']:.4g})")
    # (b) span の格子点の verdict
    if ts is not None:
        if last is None or last < span_end:
            out["incomplete"].append(f"データが span の終わり (step {span_end}) まで無い (データの最後 {last})")
        for s in span_steps:
            r = by_step[s]
            v = r[verdict_key]
            if v in FAIL_VERDICTS:
                qs = {nm: qv for nm, qv in r["by_quantity"].items()
                      if qv != "SETTLED" and (b_kinds is None or r["kinds"][nm] in b_kinds)}
                out["b_not_settled"].append({"step": s, "verdict": v, "why": r.get(why_key), "quantities": qs})
            elif v == "INDETERMINATE":
                out["incomplete"].append(f"(b) 格子点 {s} の{b_label}が INDETERMINATE — " + "; ".join((r.get(why_key) or [])[:2]))
        bb = out["b_not_settled"]
        if bb:
            out["evidence"].append(f"(b) span の {len(bb)} 格子点で{b_label}が SETTLED でない: "
                                   + ", ".join(f"{b['step']} {b['verdict']}" for b in bb[:10]) + (" …" if len(bb) > 10 else "")
                                   + f" (最初の理由: {(bb[0]['why'] or ['-'])[0]})")
    return out


def analyse(run_dir, *, span: int = SPAN_DEFAULT, ts_max: int = TS_MAX_DEFAULT, out_dir=None, write: bool = True,
            plot: bool = True, criteria=None) -> dict:
    """run を再生して仮想の停止の後を検証する。戻りは結果の dict (`virtual_stop.json` の中身)。write で出力を書く。
    criteria を渡すと run の settle_criteria.yaml の代わりにそのファイルを基準にする (run の基準は読まず書き換えない; 基準の中の
    相対のパスは run から解く)。"""
    t0 = time.perf_counter()
    run = Path(run_dir)
    if span <= 0 or ts_max < 0:
        raise ValueError(f"span {span} は正、ts_max {ts_max} は 0 以上")
    crit_path = Path(criteria) if criteria is not None else run / settle.CRITERIA_FILE
    res = {"schema": SCHEMA, "plan": PLAN, "run_dir": str(run.resolve()),
           "created": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
           "tool": "solver_density_cuda/tools/settle_virtual_stop.py", "tool_sha256": _sha256_file(Path(__file__)),
           "settle_sha256": _sha256_file(Path(settle.__file__)), "params": {"span": int(span), "ts_max": int(ts_max)},
           "reasons": [], "incomplete": [], "b_evidence": [], "notes": []}
    if criteria is not None:
        res["params"]["criteria"] = str(crit_path.resolve())
        res["notes"].append(f"基準は --criteria {crit_path} (run の {settle.CRITERIA_FILE} の代わり; run の基準は使わない)")
    timing = {}
    try:
        lc = settle.load_criteria(run, path=crit_path if criteria is not None else None)
    except settle.CriteriaError as e:
        res.update(verdict=RESULT_NA, reasons=[f"基準のファイルの誤り ({crit_path}): {e}"])
        return _finish(res, run, out_dir, write, plot, None, [], [], {}, timing, t0)
    if lc["status"] != "ok":
        res.update(verdict=RESULT_NA, reasons=[lc["message"]], criteria={"path": lc["path"], "sha256": lc.get("sha256")})
        return _finish(res, run, out_dir, write, plot, None, [], [], {}, timing, t0)
    crit = lc["criteria"]
    dt = crit["grid"]["out_step_interval"]
    start = crit["grid"]["start_step"]
    W = crit["window"]["W"]
    hold, hold_src = v2_hold(crit)
    first = first_judgeable(crit)
    targets = [q for q in crit["quantities"] if q["kind"] in settle.TARGET_KINDS]
    res["criteria"] = {"path": lc["path"], "sha256": lc["sha256"], "provisional": settle.any_provisional(crit), "W": W,
                       "min_points": crit["window"]["min_points"], "out_step_interval": dt, "start_step": start, "hold": hold,
                       "hold_source": hold_src, "first_judgeable": first,
                       "targets": [{"name": q["name"], "kind": q["kind"], "tolerance": _tol_text(q)} for q in targets],
                       "content": crit}
    if span % dt != 0:
        res["notes"].append(f"span {span} が格子の間隔 {dt} の倍数でない (span の終わりは t_s + span 以下の最後の格子点)")

    # ---------------------------------------------------------------- 1. 再生
    data = settle.RunData(run, crit)
    last, cands = data.last_step()
    res["data_last"] = {"step": last, "sources": [{"step": c[0], "what": c[1]} for c in cands]}
    grid = list(range(first, last + 1, dt)) if last is not None and last >= first else []
    rows, checks = [], []
    t1 = time.perf_counter()
    for s in grid:
        r = settle.judge(run, s, data=data, loaded=lc)
        row = summarize(r)
        row["step"] = s
        rt = summarize(r, kinds=settle.TARGET_KINDS)
        row["target_verdict"] = rt["verdict"]
        row["target_categories"] = rt["categories"]
        row["target_why"] = rt.get("why")
        rows.append(row)
        for q in r.get("quantities", []):
            base = {"step": s, "quantity": q["name"], "kind": q["kind"], "verdict": q["verdict"], "categories": "|".join(q.get("categories", []))}
            for c in q.get("checks", []):
                checks.append({**base, "check": c["check"], "value": c["value"], "limit": c["limit"], "ok": c["ok"],
                               "detail": c.get("detail", "")})
            for rr in q.get("reasons", []):
                if rr["category"] != settle.CAT_THRESHOLD:
                    checks.append({**base, "check": "", "value": "", "limit": "", "ok": "", "detail": f"[{rr['label']}] {rr['text']}"})
    timing["replay_s"] = time.perf_counter() - t1
    res["grid"] = {"first": first, "last": grid[-1] if grid else None, "n_points": len(grid)}
    if not grid:
        res["reasons"].append(f"判定できる格子点が無い (最初に判定できる格子点 {first} = 開始 + 2W + 間隔; データの最後 {last})")

    # ---------------------------------------------------------------- 2. 停止候補 t_s
    cap, cap_src = main_stage_cap(run)
    if cap is None:
        res["notes"].append(f"本段の上限が分からない ({cap_src}) — 上限の格子点の規則 (保持を上限の格子点で満たしたら要求なし) を当てずに再生")
    req, hs = emulate_v2(rows, hold, cap)
    res["v2_emulation"] = {"request": req, "hold": hs, "cap": {"value": cap, "source": cap_src},
                           "source": "design/forge_design/evaluate/stop_rule.py:624–637 (_hold_step)・1558–1588 (要求の分岐)・"
                                     "1297–1310 (judge_summary)・1181 (first)・1124 (hold)"}
    ts = req["step"] if req is not None and req["kind"] == REQ_SETTLED else None
    res["t_s"] = ts
    res["settled_since"] = req["settled_since"] if ts is not None else None
    span_end = ts + span if ts is not None else None
    res["span_end"] = span_end
    if req is None:
        if hs["reached_at_cap"] is not None:
            res["reasons"].append(f"停止候補が無い: 保持は上限の格子点 {hs['reached_at_cap']} で満たした (v2 は要求を置けない — CENSORED; "
                                  f"上限 {cap} は {cap_src})")
        elif grid:
            res["reasons"].append(f"停止候補が無い: v2 は格子点 {grid[0]}〜{grid[-1]} で要求を置かない"
                                  + (f" (保持の途中: 条件を満たし始めた {hs['start']}、保持 {hold})" if hs["start"] is not None else ""))
    elif ts is None:
        res["reasons"].append(f"停止候補が無い: v2 は step {req['step']} で {req['kind']} の要求を置いて止める (FAILED) — {req['reason']}")
    else:
        if grid and ts == grid[-1]:
            res["notes"].append(f"t_s {ts} はデータの最後の格子点"
                                + (" (本段の上限が分からないので上限の格子点の規則を当てていない; 上限の格子点なら v2 は要求を置かない)"
                                   if cap is None else ""))
        if ts > ts_max:
            res["reasons"].append(f"t_s {ts} > {ts_max} (検証の範囲 t_s + {span} がデータの範囲に収まる前提を満たさない)")

    # ---------------------------------------------------------------- 3. t_s からの差 (値は全格子点で作り、判定は span だけ)
    t2 = time.perf_counter()
    vals = {q["name"]: {} for q in targets}
    for s in grid:
        for q in targets:
            try:
                vals[q["name"]][s] = {"ok": True, "v": target_values(data, crit, q, s)}
            except settle._NonFinite as e:
                vals[q["name"]][s] = {"ok": False, "category": settle.CAT_NONFINITE, "text": str(e)}
            except settle._DataMissing as e:
                vals[q["name"]][s] = {"ok": False, "category": settle.CAT_DATA, "text": str(e)}
            except _GridMismatch as e:
                vals[q["name"]][s] = {"ok": False, "category": settle.CAT_GRID, "text": str(e)}
    by_step = {r["step"]: r for r in rows}
    sc = _span_check(ts, span, targets, vals, grid, by_step, last, "verdict", "why", "全体の判定", "t_s")
    res["b_evidence"] += sc["evidence"]
    res["incomplete"] += sc["incomplete"]
    res["span"] = {"t_s": ts, "end": span_end, "grid_points": sc["grid_points"], "n_grid_points": len(sc["grid_points"]),
                   "a": sc["a"], "b_not_settled": sc["b_not_settled"]}

    # ---------------------------------------------------------------- 3′. 目的の量だけで保持した停止候補 t_s′ (記録だけ; 判定に使わない)
    rows_t = [{"step": r["step"], "verdict": r["target_verdict"], "categories": r["target_categories"], "why": r.get("target_why")}
              for r in rows]
    req_t, hs_t = emulate_v2(rows_t, hold, cap)
    for r, rt in zip(rows, rows_t):
        r["hold_start_target_only"] = rt["hold_start"]
        r["v2_event_target_only"] = rt["v2_event"]
    tsp = req_t["step"] if req_t is not None and req_t["kind"] == REQ_SETTLED else None
    sct = _span_check(tsp, span, targets, vals, grid, by_step, last, "target_verdict", "target_why", "目的の量の判定", "t_s′",
                      b_kinds=settle.TARGET_KINDS)
    to_notes = []
    if tsp is None:
        rel = "t_s′ なし" + ((f" (目的の量だけの判定でも step {req_t['step']} で {req_t['kind']} の要求になる — {req_t['reason']})")
                            if req_t is not None else
                            (f" (保持は上限の格子点 {hs_t['reached_at_cap']} で満たした — v2 は要求を置けない)" if hs_t["reached_at_cap"] is not None
                             else (f" (格子点 {grid[0]}〜{grid[-1]} で保持を満たさない)" if grid else "")))
        outcome = "t_s′ なし"
    else:
        if ts is None:
            rel = f"t_s なし・t_s′ {tsp}"
        elif tsp == ts:
            rel = "t_s′ = t_s"
        elif tsp < ts:
            rel = f"t_s′ < t_s ({ts - tsp} step 早い)"
        else:
            rel = f"t_s′ > t_s ({tsp - ts} step 遅い)"
        outcome = "超過あり" if sct["evidence"] else ("判定できない" if sct["incomplete"] else "超過なし")
        if tsp > ts_max:
            to_notes.append(f"t_s′ {tsp} > {ts_max} (本検証なら未成立の条件)")
    timing["values_s"] = time.perf_counter() - t2
    res["target_only"] = {
        "used_for_verdict": False,
        "definition": ("残差 (residual) と流量の収支 (mass_balance) の条件を外し、目的の量 (" + ", ".join(settle.TARGET_KINDS)
                       + ") の verdict の AND (優先は settle と同じ) に停止規則 v2 と同じ保持の規則・要求の分岐を当てた停止候補 t_s′。"
                       "t_s′ から span の (a) は本検証と同じ差、(b) は全体でなく目的の量の verdict。残差の条件が独立の情報を足すかの確認の材料"),
        "request": req_t, "hold": hs_t, "t_s_prime": tsp, "settled_since": req_t["settled_since"] if tsp is not None else None,
        "span_end": tsp + span if tsp is not None else None, "relation": rel, "outcome": outcome,
        "evidence": sct["evidence"], "incomplete": sct["incomplete"],
        "span": {"grid_points": sct["grid_points"], "n_grid_points": len(sct["grid_points"]), "a": sct["a"],
                 "b_not_settled": sct["b_not_settled"]},
        "notes": to_notes}

    # ---------------------------------------------------------------- 4. 判定
    if res["reasons"]:
        verdict = RESULT_NA
    elif res["b_evidence"]:
        verdict = RESULT_B
    elif res["incomplete"]:
        verdict = RESULT_NA
    else:
        verdict = RESULT_A
    res["verdict"] = verdict
    # B の内訳 (codex diagnose 2026-10-16 settle-v2-span10k): (a) 目的の量の累積の超過と (b) 全条件の持続の違反を分けて書く
    kinds_of = {q["name"]: q["kind"] for q in crit["quantities"]}
    res["b_breakdown"] = b_breakdown(sc["a"], sc["b_not_settled"], kinds_of)
    res["target_only"]["b_breakdown"] = b_breakdown(sct["a"], sct["b_not_settled"], kinds_of)

    # ---------------------------------------------------------------- 5. 追加の記録
    rm_kinds = (settle.KIND_RESIDUAL, settle.KIND_MASS)
    extra = {"definition": "目的の量 (" + ", ".join(settle.TARGET_KINDS) + ") がすべて SETTLED なのに residual か mass_balance が SETTLED でない格子点",
             "failed": {"count": 0, "steps": [], "by_region": {}}, "indeterminate": {"count": 0, "steps": [], "by_region": {}},
             "target_only_settled_points": 0}
    for r in rows:
        s = r["step"]
        tv = [v for nm, v in r["by_quantity"].items() if r["kinds"][nm] in settle.TARGET_KINDS]
        rmv = {nm: v for nm, v in r["by_quantity"].items() if r["kinds"][nm] in rm_kinds}
        r["targets_all_settled"] = bool(tv) and all(v == "SETTLED" for v in tv)
        r["res_mass_failed"] = any(v in FAIL_VERDICTS for v in rmv.values())
        r["res_mass_indeterminate"] = any(v == "INDETERMINATE" for v in rmv.values())
        r["region"] = _region(s, ts, span)
        r["region_ts_prime"] = _region(s, res["target_only"]["t_s_prime"], span)
        if not r["targets_all_settled"]:
            continue
        extra["target_only_settled_points"] += 1
        for key, flag in (("failed", r["res_mass_failed"]), ("indeterminate", r["res_mass_indeterminate"] and not r["res_mass_failed"])):
            if flag:
                e = extra[key]
                e["count"] += 1
                e["steps"].append({"step": s, "region": r["region"], "residual_mass": rmv})
                e["by_region"][r["region"]] = e["by_region"].get(r["region"], 0) + 1
    res["extra_record"] = extra
    return _finish(res, run, out_dir, write, plot, data, rows, checks,
                   {"dev_rows": sc["dev_rows"], "series": sc["series"], "dev_rows_tp": sct["dev_rows"]}, timing, t0)


def b_breakdown(span_a: dict, b_not: list, kinds_of: dict) -> dict:
    """B の内訳 (codex diagnose 2026-10-16 settle-v2-span10k の「B の解釈」): (a) 目的の量の t_s からの累積の超過 (量ごとの超過の
    格子点の数・最初・最大の比) と、(b) span の格子点の全条件の持続の違反 (SETTLED でなかった条件 = 量ごとの格子点の数・verdict・step)
    を分ける。(b) だけの B は目的の量を早く止めすぎた証明ではない (残差・流量の条件の再逸脱でも B になる)。"""
    a_q = {}
    for nm, sm in span_a.items():
        if sm.get("n_exceeded"):
            a_q[nm] = {"kind": sm.get("kind"), "n_exceeded": sm["n_exceeded"], "first_exceeded": sm.get("first_exceeded"),
                       "max_ratio": sm.get("max_ratio"), "at_step": sm.get("at_step"), "where": sm.get("where")}
    by_q = {}
    for b in b_not:
        for nm, v in (b.get("quantities") or {}).items():
            e = by_q.setdefault(nm, {"kind": kinds_of.get(nm), "n_grid_points": 0, "verdicts": {}, "first_step": b["step"], "steps": []})
            e["n_grid_points"] += 1
            e["verdicts"][v] = e["verdicts"].get(v, 0) + 1
            if len(e["steps"]) < 20:
                e["steps"].append(b["step"])
    return {"a_cumulative_exceedance": {"present": bool(a_q), "quantities": a_q,
                                        "definition": "目的の量の t_s の値からの差が許容を超えた (累積のドリフト; 節点の値の非有限を含む)"},
            "b_persistence_violation": {"present": bool(b_not), "n_grid_points": len(b_not),
                                        "first_step": b_not[0]["step"] if b_not else None, "by_quantity": by_q,
                                        "definition": "span の格子点で全条件の判定が SETTLED でない (どの条件が崩れたか)。これだけの B は"
                                                      "目的の量を早く止めすぎた証明ではない"}}


def _breakdown_lines(bd: dict, head: str) -> list:
    a, b = bd["a_cumulative_exceedance"], bd["b_persistence_violation"]
    out = []
    if a["present"]:
        out.append(f"{head} (a) 目的の量の累積の超過: あり — " + "; ".join(
            f"{nm} ({x['kind']}) 超過 {x['n_exceeded']} 格子点、最初 step {x['first_exceeded']}、最大の比 "
            + (f"{x['max_ratio']:.4g}" if isinstance(x.get("max_ratio"), (int, float)) else "-") + f" (step {x['at_step']})"
            for nm, x in a["quantities"].items()))
    else:
        out.append(f"{head} (a) 目的の量の累積の超過: なし")
    if b["present"]:
        out.append(f"{head} (b) 全条件の持続の違反: あり — span の {b['n_grid_points']} 格子点 (最初 step {b['first_step']}); 崩れた条件: " + "; ".join(
            f"{nm} ({x['kind']}) {x['n_grid_points']} 格子点 {x['verdicts']} 最初 step {x['first_step']}" for nm, x in b["by_quantity"].items()))
        if not a["present"]:
            out.append(f"{head}の注: (a) の超過は無く (b) だけ — 目的の量を早く止めすぎた証明ではない (条件の再逸脱)")
    else:
        out.append(f"{head} (b) 全条件の持続の違反: なし")
    return out


DEV_COLS = ["step", "region", "quantity", "kind", "status", "diff", "tol", "ratio", "exceeded", "where", "signed", "value", "ref",
            "x_m", "r_m", "maxabs_diff", "halfamp_diff", "detail"]


def _finish(res, run, out_dir, write, plot, data, rows, checks, dev, timing, t0) -> dict:
    res["provisional"] = bool((res.get("criteria") or {}).get("provisional"))
    od = Path(out_dir) if out_dir else run / OUT_DIR
    if write:
        t3 = time.perf_counter()
        od.mkdir(parents=True, exist_ok=True)
        _write_grid_csv(od / "grid_verdicts.csv", rows)
        _write_rows(od / "grid_checks.csv", checks, ["step", "quantity", "kind", "verdict", "categories", "check", "value", "limit", "ok",
                                                       "detail"])
        _write_rows(od / "span_deviation.csv", dev.get("dev_rows", []), DEV_COLS)
        _write_rows(od / "span_deviation_ts_prime.csv", dev.get("dev_rows_tp", []), DEV_COLS)
        figs = []
        if plot and rows:
            try:
                figs = _plots(od, res, rows, dev.get("series", {}))
            except Exception as e:  # noqa: BLE001 — 図の失敗は判定を変えない (記録する)
                res["notes"].append(f"図を作れない ({type(e).__name__}: {e})")
        res["figures"] = [str(p) for p in figs]
        timing["outputs_s"] = time.perf_counter() - t3
        timing["total_s"] = time.perf_counter() - t0
        res["timing_s"] = timing
        cp = Path((res.get("criteria") or {}).get("path") or (res.get("params") or {}).get("criteria") or run / settle.CRITERIA_FILE)
        extra = [cp] if cp.exists() else []         # 使った基準のファイル (--criteria ならそれ; run の基準ではない)
        res["files"] = settle.files_record(data, extra)
        (od / "VERDICT.txt").write_text("\n".join(verdict_lines(res)) + "\n")
        tmp = od / "virtual_stop.json.tmp"
        tmp.write_text(json.dumps(settle._jsonable(res), ensure_ascii=False, indent=1) + "\n")
        tmp.replace(od / "virtual_stop.json")
        res["out_dir"] = str(od)
    else:
        timing["total_s"] = time.perf_counter() - t0
        res["timing_s"] = timing
    return res


def verdict_lines(res: dict) -> list:
    """VERDICT.txt の中身 (1 行目が VERDICT の行)。"""
    out = [f"VERDICT: {res['verdict']}"]
    out.append(f"run: {res['run_dir']}")
    c = res.get("criteria") or {}
    if res.get("t_s") is not None:
        out.append(f"t_s: {res['t_s']} (条件を満たし始めた step {res['settled_since']}、保持 {c.get('hold')}) / span {res['params']['span']} "
                   f"({res['t_s']}〜{res['span_end']}) / t_s の上限 {res['params']['ts_max']}")
    else:
        req = (res.get("v2_emulation") or {}).get("request")
        out.append("t_s: なし" + (f" (v2 は step {req['step']} で {req['kind']} の要求)" if req else ""))
    for k, lab in (("reasons", "未成立の理由"), ("b_evidence", "超過"), ("incomplete", "欠け・判定不能")):
        for x in res.get(k) or []:
            out.append(f"{lab}: {x}")
    if res.get("verdict") == RESULT_B and res.get("b_breakdown"):
        out += _breakdown_lines(res["b_breakdown"], "B の内訳")
    a = ((res.get("span") or {}).get("a")) or {}
    for nm, s in a.items():
        if s.get("max_ratio") is not None:
            out.append(f"(a) {nm}: {s['status']}、span の最大の比 {s['max_ratio']:.4g} (step {s['at_step']}、{s['where']}; {s['tolerance']})")
        else:
            out.append(f"(a) {nm}: {s.get('status')}")
    ex = res.get("extra_record")
    if ex:
        out.append(f"追加の記録: 目的の量だけなら SETTLED なのに残差か流量で落ちた格子点 {ex['failed']['count']} "
                   f"(区間別 {ex['failed']['by_region']})、判定不能だった格子点 {ex['indeterminate']['count']} "
                   f"(目的の量がすべて SETTLED の格子点 {ex['target_only_settled_points']})")
    to = res.get("target_only")
    if to:
        hd = "参考 (判定に使わない) 目的の量だけで保持した停止候補"
        if to["t_s_prime"] is not None:
            out.append(f"{hd} t_s′: {to['t_s_prime']} (条件を満たし始めた step {to['settled_since']}、span {to['t_s_prime']}〜{to['span_end']}) — "
                       f"{to['relation']}、t_s′ からの検証: {to['outcome']}")
        else:
            out.append(f"{hd} t_s′: {to['relation']}")
        for k, lab in (("evidence", "超過"), ("incomplete", "欠け・判定不能")):
            for x in to.get(k) or []:
                out.append(f"参考 t_s′ の{lab}: {x}")
        for nm, sm in ((to.get("span") or {}).get("a") or {}).items():
            if sm.get("max_ratio") is not None:
                out.append(f"参考 t_s′ (a) {nm}: {sm['status']}、span の最大の比 {sm['max_ratio']:.4g} (step {sm['at_step']}、{sm['where']})")
        for x in to.get("notes") or []:
            out.append(f"参考 t_s′ の note: {x}")
    for x in res.get("notes") or []:
        out.append(f"note: {x}")
    if res.get("provisional"):
        src = ((c.get("content") or {}).get("source") or {}) if isinstance(c.get("content"), dict) else {}
        ver = f"{src['preset_version']} の値" if src.get("preset_version") else "基準のファイルの値"
        out.append(f"閾値: {settle.PROVISIONAL_TEXT} — この検証では {ver}に固定 (基準 {c.get('path')}; plan §6 V3(b)・(1″))")
    return out


# ====================================================================== 出力
def _cell(v):
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, (list, tuple, dict)):
        return json.dumps(settle._jsonable(v), ensure_ascii=False)
    return "" if v is None else str(v)


def _write_rows(p: Path, rows: list, cols: list):
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([_cell(r.get(c)) for c in cols])


def _write_grid_csv(p: Path, rows: list):
    names = []
    for r in rows:
        for nm in r["by_quantity"]:
            if nm not in names:
                names.append(nm)
    cols = ["step", "verdict", "categories", "region", "hold_start", "v2_event", "targets_all_settled", "res_mass_failed",
            "res_mass_indeterminate", "target_verdict", "region_ts_prime", "hold_start_target_only", "v2_event_target_only"] + [
        f"v:{nm}" for nm in names]
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            base = [r["step"], r["verdict"], "|".join(r["categories"]), r.get("region", ""), _cell(r.get("hold_start")),
                    r.get("v2_event", ""), r.get("targets_all_settled", ""), r.get("res_mass_failed", ""),
                    r.get("res_mass_indeterminate", ""), r.get("target_verdict", ""), r.get("region_ts_prime", ""),
                    _cell(r.get("hold_start_target_only")), r.get("v2_event_target_only", "")]
            w.writerow(base + [r["by_quantity"].get(nm, "") for nm in names])


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    if Path(FONT).is_file():
        font_manager.fontManager.addfont(FONT)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def _span_marks(ax, res):
    ts, se = res.get("t_s"), res.get("span_end")
    if ts is not None:
        ax.axvspan(ts, se, color="tab:blue", alpha=0.08, lw=0)
        ax.axvline(ts, color="tab:blue", lw=0.8, ls=":")
    ss = res.get("settled_since")
    if ss is not None:
        ax.axvline(ss, color="tab:green", lw=0.8, ls=":")
    tsp = (res.get("target_only") or {}).get("t_s_prime")
    if tsp is not None and tsp != ts:
        ax.axvline(tsp, color="tab:purple", lw=0.9, ls="--")


def _plots(od: Path, res: dict, rows: list, series: dict) -> list:
    plt = _pyplot()
    figs = []
    title = f"{Path(res['run_dir']).name}: VERDICT {res['verdict']}" + (f"  t_s {res['t_s']}・span {res['params']['span']}"
                                                                       if res.get("t_s") is not None else "  (t_s なし)")
    names = [nm for nm, s in series.items() if s["steps"]]
    if names:
        fig, axes = plt.subplots(len(names), 1, figsize=(10, 2.6 * len(names)), sharex=True, squeeze=False)
        for ax, nm in zip(axes[:, 0], names):
            s = series[nm]
            x = np.asarray(s["steps"])
            kind = s["kind"]
            if kind == "exit_core_mach":
                ax.plot(x, s["signed"], "o-", ms=3, lw=1, label="M − M(t_s)")
                tol = np.asarray(s["tol"])
                ax.plot(x, tol, "r--", lw=1, label="± 許容")
                ax.plot(x, -tol, "r--", lw=1)
                ax.set_ylabel(f"ΔM {UNIT[kind]}")
            elif kind == "wall_heat_nodes":
                ax.plot(x, s["ratio"], "o-", ms=3, lw=1, label="max_i |ΔQ_i| / (q_tol·A_i + ε·|Q_i(t_s)|)")
                ax.axhline(1.0, color="r", ls="--", lw=1, label="許容")
                ax.set_ylabel(f"比 {UNIT[kind]}")
            else:
                ax.plot(x, s["diff"], "o-", ms=3, lw=1, label=("max_k |d_k − d_k(t_s)|" if kind == "line_mach" else "max_i |θ_i − θ_i(t_s)|"))
                ax.plot(x, s["tol"], "r--", lw=1, label="許容")
                ax.set_ylabel(f"差 {UNIT[kind]}")
            _span_marks(ax, res)
            ax.set_title(f"{nm} ({LABEL.get(kind, kind)})", fontsize=9, loc="left")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7, loc="upper left")
        axes[-1, 0].set_xlabel("step (青の帯 = t_s〜t_s + span、緑の点線 = 条件を満たし始めた step、紫の破線 = 目的の量だけの t_s′ (参考))")
        fig.suptitle("t_s の値からの差と許容の線 — " + title, fontsize=10)
        fig.tight_layout()
        p = od / "deviation_from_ts.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        figs.append(p)
    # 格子点ごとの verdict
    qn = []
    for r in rows:
        for nm in r["by_quantity"]:
            if nm not in qn:
                qn.append(nm)
    labels = ["全体"] + qn
    colors = {"SETTLED": "tab:green", "NOT_SETTLED": "tab:orange", "INDETERMINATE": "0.6", "DIVERGED": "tab:red"}
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(labels) + 1.6))
    for v, col in colors.items():
        xs, ys = [], []
        for r in rows:
            vals = [r["verdict"]] + [r["by_quantity"].get(nm) for nm in qn]
            for k, vv in enumerate(vals):
                if vv == v:
                    xs.append(r["step"])
                    ys.append(k)
        if xs:
            ax.scatter(xs, ys, s=14, marker="s", color=col, label=v)
    _span_marks(ax, res)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("格子点の step (終端)")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("格子点ごとの量ごとの判定 — " + title, fontsize=9)
    fig.tight_layout()
    p = od / "grid_verdicts.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    figs.append(p)
    return figs


# ====================================================================== CLI
class _Parser(argparse.ArgumentParser):
    def error(self, message):  # 引数の誤りは未成立 (3)
        self.print_usage(sys.stderr)
        print(f"{self.prog}: 引数の誤り: {message}", file=sys.stderr)
        print(f"VERDICT: {RESULT_NA}")
        sys.exit(EXIT_CODE[RESULT_NA])


def main(argv=None) -> int:
    ap = _Parser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--span", type=int, default=SPAN_DEFAULT, help=f"t_s から検証する step 数 (既定 {SPAN_DEFAULT})")
    ap.add_argument("--ts-max", type=int, default=TS_MAX_DEFAULT, help=f"t_s の上限 (超えたら未成立; 既定 {TS_MAX_DEFAULT})")
    ap.add_argument("--out", default=None, help=f"出力の置き場所 (既定 <run>/{OUT_DIR}; --criteria のときは必須)")
    ap.add_argument("--criteria", default=None,
                    help=f"run の {settle.CRITERIA_FILE} の代わりに使う基準のファイル (run の基準は読まず書き換えない; --out が要る)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--json", action="store_true", help="結果の JSON を標準出力に出す")
    a = ap.parse_args(argv)
    run = Path(a.run_dir)
    if not run.is_dir():
        print(f"ERROR: {run} はディレクトリでない")
        print(f"VERDICT: {RESULT_NA}")
        return EXIT_CODE[RESULT_NA]
    if a.criteria is not None:
        # 別の基準の結果で run の基準の結果 (<run>/settle_virtual_stop/) を上書きしない
        if a.out is None:
            print(f"ERROR: --criteria には --out が要る (既定の {run / OUT_DIR} は run の基準の結果の置き場所)")
            print(f"VERDICT: {RESULT_NA}")
            return EXIT_CODE[RESULT_NA]
        if Path(a.out).resolve() == (run / OUT_DIR).resolve():
            print(f"ERROR: --out {a.out} は run の基準の結果の置き場所 ({run / OUT_DIR}) — 別の場所にする")
            print(f"VERDICT: {RESULT_NA}")
            return EXIT_CODE[RESULT_NA]
    try:
        res = analyse(run, span=a.span, ts_max=a.ts_max, out_dir=a.out, write=True, plot=not a.no_plot, criteria=a.criteria)
    except ValueError as e:
        print(f"ERROR: {e}")
        print(f"VERDICT: {RESULT_NA}")
        return EXIT_CODE[RESULT_NA]
    lines = verdict_lines(res)
    for ln in lines[1:]:
        print(ln)
    print(f"所要: {', '.join(f'{k} {v:.2f} s' for k, v in (res.get('timing_s') or {}).items())}")
    print(f"結果: {res.get('out_dir')}")
    if a.json:
        print(json.dumps(settle._jsonable(res), ensure_ascii=False, indent=1))
    print(lines[0])
    return EXIT_CODE[res["verdict"]]


if __name__ == "__main__":
    sys.exit(main())
