r"""風洞ノズルの設計の判定 (qualification) (plan tooling-nozzle-card-local-mesh §5.1 #6b・§4.3・§6 V5; orchestrator plan
tooling-design-orchestrator §4.2 の先行)。

計算の成立 (停止規則・関門・`judge_run`) とは別に、**設計の要求を満たしたか**を 1 本の NS run で判定する。
判定 (`result`):
  accepted            要求の公差 (V5b・V5c)・壁解像・寸法の関門・必要な成果物をすべて満たし、計算の成立 (停止規則の到達) も確かめた
  back_to_correction  評価済みで公差 (V5b・V5c) が未達、または壁解像・寸法の関門が FAIL (補正へ戻す)
  rejected            back_to_correction に当たり、反復の上限 (`--iteration` ≥ `--max-iterations`) に達した
  incomplete          必要な判定 (V5b・V5c・壁解像・計算の成立) か成果物が欠けた
次の手 (`next`): accepted → done、back_to_correction → recalibrate (下の規則をすべて満たす) か escalate (人に返す)、
rejected → escalate、incomplete → make_deliverables (流れの判定はすべて合格で成果物だけが無い) か escalate。

判定の中身 (定義は既存の部品と同じ):
  - 場: `--step` → 停止規則の状態 (`STOP_RULE_STATE.json`) の stop_step の全場 → 最後の全場 `res_<N>.h5` (res_0 を除く)。
  - V5b 出口コア M: 最後の列の η ∈ [0.05, 0.7] の節点の M の平均 (`stop_rule.recorded_quantities` = `nozzle_report.metrics` の
    exit_core.M_mean)。許容 |M/M_d − 1| ≤ tol (既定 1e−3 = ±0.1 %、plan §6 V5b; 2026-10-16 ユーザ決定、前は 5e−4)。
  - V5c 近軸の線の一様性: 線 η (既定 0.03) の M を `nozzle_report.eta_line` で x 窓に n 点 (既定 2001) 取り、
    最大偏差 max|100(M/M_d − 1)| と片振幅 100(max M − min M)/(2 M_d) (`settle.line_quantities` = `v5c_series.py` の式)。
    許容 0.1 %・0.1 % (plan §6 V5c; 片振幅は 2026-10-16 ユーザ決定、前は 0.05 %)。η 0.1 と軸 (η 0) は比較用に併記 (判定に使わない)。補間に軸の節点が重み > 1e−6 で入る
    (第 1 層の η が線より上) か、窓が格子の外・非有限なら判定不能 (plan §6 V5c、codex M3)。
    **この量は「近軸の線 (η) の一様性」であって試験部全体の一様性ではない** (plan §6 V5c)。
    線の定常の保証 (plan §5.1 #7a、2026-10-17 codex result M1): 評価の線 (η・x 窓・点数) が停止規則の監視の line_mach と違う・
    line_mach が基準に無い・評価の場の step が停止の step と違う、のどれかなら判定不能 (→ incomplete)。**監視の線は停止の時の判定の
    記録 (下の「計算の成立」の結び付き) の基準の中身から取り、その量が SETTLED であることも要る**。停止の時の判定・基準・入力が今の run と
    結び付かない (照合できない・記録が無い) ときは監視による保証にしない (2026-10-17 codex result 2 回目 M1: 今の settle_criteria.yaml の
    読み取り値だけを η 0.06 に変えると、古い停止の記録のまま accepted になった)。保証が無い線は、本段の全場 (MAIN_STAGE_STARTED の
    start_step < N ≤ 評価の場の step) の系列の最大偏差・片振幅を `check_quasisteady.py --series-csv` が両方 STEADY と判定したときだけ
    合否に進める (VERDICT・系列の CSV を記録)。全場が 4 枚未満なら判定不能。
  - 壁解像: `check_wall_resolution.py <run> --step <場の step> --groups <bcond の wall* 全部> --weight area` の VERDICT
    (y₁⁺ > 1 の面積 ≤ 2 %、plan §6 V4)。ツールが評価した壁と step (EVALUATED 行) が bcond の全部の壁・場の step・停止の step と
    一致しなければ判定不能 (plan §5.1 #7b)。
  - 計算の成立: `settle_criteria.yaml` のある run は停止規則 v2 の REACHED と停止の step の再評価 SETTLED。**その前に停止の時の
    判定・基準・入力が今の run と結び付くこと** (受け渡しの関門と同じ照合 `deltastar_loop._v2_binding`: 状態・停止の時の判定の専用
    ファイル (古い run は停止の時のままの SETTLE_VERDICT.json)・基準の sha256・本段の起動の記録・再評価が使った入力の sha256・
    STOPPED_AT・停止の場)。結び付かない・状態や停止の時の判定の記録が無いときは INDETERMINATE (→ incomplete; plan §5.1 #7a の続き、
    2026-10-17 codex result 2 回目 M1)。旧方式は `STOP_RULE_STATE.json` の REACHED。状態が無い旧い run (記録の写し等) は UNVERIFIED
    (accepted にも自動の較正し直しにもしない)。
  - 寸法の関門と成果物: `<run>/deliverables/` の `DIMENSIONS.json` (verdict PASS、カードを渡したらその sha256 と一致) と `MANIFEST.json`
    (verdict PASS・失敗/飛ばした部品なし・必要な部品 `deliverables.PARTS` がそろい・列挙したファイルの sha256 が一致)。無ければ ABSENT。
    **成果物の生成元の照合** (plan §5.1 #6b の続き、2026-10-17 codex result 3 回目 Major: run_0014 の評価に run_0010 の成果物を渡すと
    accepted になった): MANIFEST.json の `provenance` (生成元の run の壁・格子・設定、場の step と sha256・壁ダンプ、停止の step、カード・
    問題 YAML の sha256) を今の評価の対象 (評価の場の step と `res_<step>.h5`、run の壁・格子・設定、渡したカード・問題 YAML) と内容で
    照合する (`provenance_check`; 置き場所・run の名前は見ない — 別の置き場所への写しは中身が同じなら通す)。一致しない・生成元の記録が
    無い (古い成果物)・照合できない → 成果物 INDETERMINATE・寸法の関門も INDETERMINATE (その成果物の寸法はこの run の判定にしない)。
    結果は incomplete (次の手は既存の規則どおり escalate; 理由に作り直しを書く)。run の名前の不一致は警告のまま (合否は内容の照合で決める)。
    **配る msh と計算の格子の元の msh の照合** (§5.1 #6c2 の再オープン、2026-10-17 codex result 4 回目 Major: 記録していた格子は h5 だけで、
    配る msh・Fluent・寸法の関門の入力の run の `nozzle.msh` と結ぶ検査が無く、MANIFEST の `msh_sha256_matches` が false でも、
    内部の配点だけ違う msh でも accepted だった): 準備の記録 (`prepare_info.json` の `mesh.hashes.msh_sha256`; h5 はこの msh から変換
    された) を基準に、今の run の `nozzle.msh`・配る `mesh/nozzle.msh` (実ファイルと MANIFEST の記録)・生成元の記録 (v2 の `msh`・
    `prepare_msh`; v1 は MANIFEST の `parts.mesh.msh41_sha256` などの既存の記録) の sha256 が一致することを要る (`_msh_binding`)。
    書き出しのときの一致の印 (`msh_sha256_matches`・`matches_prepare`) が false なら不一致。不一致 → MISMATCH、準備の記録に msh の
    sha256 が無い・run の `nozzle.msh` が無い・書き出しのときの msh の記録が無い → UNCHECKABLE (どちらも成果物 INDETERMINATE → incomplete)。

公差・窓の出所 (どれを使ったかを記録の `tolerances.*.source` に書く; 上が優先):
  1. CLI (`--exit-tol-rel`・`--line-eta`・`--line-x-range-m X0 X1`・`--line-tol-maxabs`・`--line-tol-halfamp`)
  2. カードの `qualification:` (キー exit_core_M_tol_rel・line_eta・line_x_range_m・line_tol_maxabs_pct・line_tol_halfamp_pct;
     **要求カードの書式 forge.wt_card/v1 はまだこのキーを受け付けない** — 生成器に通すカードには書けない)
  3. 問題 YAML (`--problem`) の `qualification:` (同じキー)、線 (η・x 窓) は `evaluate.settle.line_eta`・`x_range_m`
     (監視の line_mach と同じ線; plan §5.1 #5c)
  4. run の `settle_criteria.yaml` の line_mach の η・x 窓 (準備のときに解決した値)
  5. 既定: 公差は plan §6 V5b・V5c、線は雛形 axismach_final_v1 の evaluate.settle (preset wt_nozzle_v2 の線: η 0.03・
     x = (x_E + 2, x_F − 1) r_t・2001 点)。M_d はカードの M_exit (無ければ prepare_info の Md; 両方あれば一致を要る)。

較正し直しの規則 (plan §4.3 の末尾、codex diagnose 2026-10-15 m55-three-coef M2; 1 つでも欠けたら escalate):
  (a) V5b・V5c の両方を判定した (判定不能でない)
  (b) 前の更新で δ の誤差が縮んだ: この run の max|δ_in − δ_E| (x ≥ 0) < 前の run (`--prev-run`) の同じ量
  (c) 同じ 3 係数で形の誤差をさらに減らせる: この run で較正し直した射影 (`contur_calibrate`) の最大誤差 ≤ (1 − r_min) × 今の最大誤差
      (r_min = `--recal-min-reduction`、既定 0.5。**plan に数値が無い暫定の値** — 2026-10-15 の判断 (0.624 → 0.612 mm は回さない) を
      越えない側に置いた)
  (d) 較正し直しはまだ (`--recalibrations-done` < `--max-recalibrations`、既定 1 回まで)
  (e) 計算の成立を確かめた (上の PASS; UNVERIFIED の run は人に返す)
壁解像・寸法の FAIL は δ の較正で直す問題でない (壁解像は plan §4.2 の表の作り直し、寸法は寸法の解き直し) ので escalate。

使い方 (design/ で、または PYTHONPATH=design):
  design/.venv-opt/bin/python -m forge_design.feedback.wt_qualify RUN_DIR --card CARD.yaml [--problem P.yaml]
      [--prev-run PREV_RUN] [--iteration K --max-iterations N] [--recalibrations-done R]
      [--delta-next CSV] [--prev-delta-next CSV] [--step S] [--deliverables DIR] [--out-dir DIR | --no-write] [--json]
出力 (`--out-dir`、既定 `<run>/wt_qualify/`): `WT_QUALIFICATION.json`・`check_wall_resolution.txt`、recalibrate なら
`contur_overrides.yaml` (新しい係数のカードの overrides の断片)、監視されていない V5c の線の系列を作ったら `v5c_line_series.csv`。
終了コード: 0 accepted、2 back_to_correction で recalibrate、3 escalate (back_to_correction・rejected)、4 incomplete、1 入力の不備。
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from . import contur_calibrate as CC

FORMAT = "forge.wt_qualify/v1"
PLAN = ("plans/accepted/tooling-nozzle-card-local-mesh.md §5.1 #6b・§4.3・§6 V5; "
        "plans/active/tooling-design-orchestrator.md §4.2")
OUT_DIRNAME = "wt_qualify"
OUT_JSON = "WT_QUALIFICATION.json"
REPO = Path(__file__).resolve().parents[3]
TOOLS = REPO / "solver_density_cuda" / "tools"

EXIT_ETA_RANGE = (0.05, 0.7)            # 出口コアの η (stop_rule.recorded_quantities・settle の exit_eta_range と同じ)
COMPARE_ETAS = (0.1, 0.0)               # 比較用に併記する線 (判定に使わない; plan §6 V5c)
DEFAULT_TOLS = {                        # plan §6 V5b・V5c (ユーザ了承 2026-10-15「よい」)
    "exit_core_M_tol_rel": (1.0e-3, "plan §6 V5b (5.5 ± 0.1 %; 2026-10-16 ユーザ決定、前は ± 0.05 %)"),
    "line_tol_maxabs_pct": (0.1, "plan §6 V5c (最大絶対偏差 ≤ 0.1 %)"),
    "line_tol_halfamp_pct": (0.1, "plan §6 V5c (片振幅 ≤ 0.1 %; 2026-10-16 ユーザ決定、前は 0.05 %)"),
}
QUAL_KEYS = ("exit_core_M_tol_rel", "line_eta", "line_x_range_m", "line_tol_maxabs_pct", "line_tol_halfamp_pct")
DEFAULT_PRESET = "wt_nozzle_v3"         # 雛形 axismach_final_v1 の evaluate.settle の既定 (2026-10-16 から v3; 線の定義は v2 と同じ)
RECAL_MIN_REDUCTION = 0.5               # 規則 (c) の暫定の値 (plan に数値が無い; 上の docstring)
MAX_RECALIBRATIONS = 1                  # plan §4.3「1 回まで」
MAX_ITERATIONS = 3                      # 手順書 design-chain-quickstart §8 の --max-iterations 3 (承認の既定; plan tooling-design-orchestrator §4.7.1)
QUAL_KEYS_V2 = ("line_n", "optional_parts")   # 要求カードの版 2 の qualification だけが書けるキー (公差の解決とは別に読む)
RESULTS = ("accepted", "back_to_correction", "rejected", "incomplete")
EXIT_CODES = {"accepted": 0, "recalibrate": 2, "escalate": 3, "incomplete": 4}


def _settle():
    from ..evaluate import settle_presets as SP
    return SP._settle()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ------------------------------------------------------------------ 入力
def _load_yaml(path):
    return yaml.safe_load(Path(path).read_text()) if path is not None else None


def _qual_block(doc, where: str) -> dict:
    """`qualification:` のブロックを読む (無ければ {})。知らないキー・形の違いは ValueError。"""
    q = (doc or {}).get("qualification") if isinstance(doc, dict) else None
    if q is None:
        return {}
    if not isinstance(q, dict):
        raise ValueError(f"{where} の qualification が mapping でない ({q!r})")
    allowed = QUAL_KEYS + (QUAL_KEYS_V2 if where == "カード" else ())
    bad = [k for k in q if k not in allowed]
    if bad:
        raise ValueError(f"{where} の qualification の未知のキー {bad} (書けるのは {list(allowed)})")
    out = {}
    for k, v in q.items():
        if k == "optional_parts":
            out[k] = list(v)
        elif k == "line_n":
            out[k] = int(v)
        elif k == "line_x_range_m":
            if not (isinstance(v, (list, tuple)) and len(v) == 2):
                raise ValueError(f"{where} の qualification.line_x_range_m が [x0, x1] でない ({v!r})")
            out[k] = [float(v[0]), float(v[1])]
        else:
            out[k] = float(v)
    return out


def _settle_line_from_criteria(run_dir) -> dict | None:
    p = Path(run_dir) / "settle_criteria.yaml"
    if not p.is_file():
        return None
    c = yaml.safe_load(p.read_text()) or {}
    for q in c.get("quantities") or []:
        if isinstance(q, dict) and q.get("kind") == "line_mach":
            return {"eta": float(q["eta"]), "x_range_m": [float(q["x_range_m"][0]), float(q["x_range_m"][1])],
                    "n": int(q.get("n", 2001)), "Md": q.get("Md")}
    return None


def resolve_tolerances(run_dir, info: dict, card: dict | None, problem: dict | None, cli: dict) -> dict:
    """公差・窓を出所の順 (CLI → カード → 問題 YAML → run の settle_criteria.yaml → 既定) で決める。戻り: {キー: {value, source}}。"""
    from ..evaluate import settle_presets as SP
    cq = _qual_block(card, "カード")
    pq = _qual_block(problem, "問題 YAML")
    es = (((problem or {}).get("evaluate") or {}).get("settle")) if isinstance(problem, dict) else None
    es = es if isinstance(es, dict) else {}
    sc = _settle_line_from_criteria(run_dir)
    line_def = SP.PRESETS[DEFAULT_PRESET]["line"]
    out = {}

    def pick(key, default_value, default_src, from_settle=None, from_criteria=None):
        for val, src in ((cli.get(key), "CLI"), (cq.get(key), "カードの qualification"), (pq.get(key), "問題 YAML の qualification"),
                         (from_settle, "問題 YAML の evaluate.settle (監視の line_mach と同じ線)"),
                         (from_criteria, "run の settle_criteria.yaml の line_mach")):
            if val is not None:
                out[key] = {"value": val, "source": src}
                return
        out[key] = {"value": default_value, "source": default_src}

    for k, (v, src) in DEFAULT_TOLS.items():
        pick(k, v, src)
    pick("line_eta", float(line_def["eta"]), f"既定 (雛形 axismach_final_v1 の evaluate.settle preset {DEFAULT_PRESET} の線 η)",
         from_settle=None if es.get("line_eta") is None else float(es["line_eta"]),
         from_criteria=None if sc is None else sc["eta"])
    xr_settle = None if es.get("x_range_m") is None else [float(es["x_range_m"][0]), float(es["x_range_m"][1])]
    dom = (info.get("physical_wall") or {}).get("domain")
    if info.get("x_E") is not None and isinstance(dom, (list, tuple)) and len(dom) == 2:
        S = float(info["scale_m"])
        o0, o1 = line_def["x_offsets_rt"]
        xr_def = [(float(info["x_E"]) + o0) * S, (float(dom[1]) + o1) * S]
        src_def = (f"既定 (preset {DEFAULT_PRESET} の線: (x_E + {o0:g}, x_F {o1:+g}) r_t × scale_m; x_E {float(info['x_E'])!r}、"
                   f"x_F = prepare_info.physical_wall.domain[1] {float(dom[1])!r})")
    else:
        xr_def, src_def = None, "既定を作れない (prepare_info に x_E か physical_wall.domain が無い)"
    pick("line_x_range_m", xr_def, src_def, from_settle=xr_settle, from_criteria=None if sc is None else sc["x_range_m"])
    if "line_n" in cq:
        out["line_n"] = {"value": int(cq["line_n"]), "source": "カードの qualification"}
    else:
        out["line_n"] = {"value": int(line_def["n"]) if sc is None else int(sc["n"]),
                         "source": "run の settle_criteria.yaml" if sc is not None else f"既定 (preset {DEFAULT_PRESET})"}
    out["settle_line"] = sc
    return out


# ------------------------------------------------------------------ 場・計算の成立
def _res_steps(run_dir) -> list[int]:
    return sorted(int(m.group(1)) for f in Path(run_dir).iterdir() if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)) and m.group(1) != "0")


def select_field(run_dir, step: int | None = None) -> dict:
    """判定に使う全場。戻り: {step, res, source}。"""
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    if step is not None:
        if not (run_dir / f"res_{int(step)}.h5").is_file():
            raise ValueError(f"--step {step}: res_{int(step)}.h5 が無い")
        return {"step": int(step), "res": f"res_{int(step)}.h5", "source": "CLI (--step)"}
    st = run_dir / SR.STATE
    if st.is_file():
        s = json.loads(st.read_text())
        k = s.get("stop_step")
        if isinstance(k, int) and (run_dir / f"res_{k}.h5").is_file():
            return {"step": k, "res": f"res_{k}.h5", "source": f"{SR.STATE} の stop_step (status {s.get('status')})"}
    steps = _res_steps(run_dir)
    if not steps:
        raise ValueError(f"{run_dir} に全場 res_<N>.h5 (N > 0) が無い")
    why = "停止規則の stop_step の全場が無い" if st.is_file() else "停止規則の止めた場の記録が無い"
    return {"step": steps[-1], "res": f"res_{steps[-1]}.h5", "source": f"最後の全場 ({why})"}


def stop_binding(run_dir) -> dict:
    """停止の時の判定・基準・入力と今の run の結び付き (plan §5.1 #7a の続き、2026-10-17 codex result 2 回目 M1)。受け渡しの関門の照合
    `deltastar_loop._v2_binding` をそのまま使う (状態・停止の時の判定の専用ファイル・基準の sha256・本段の起動・再評価の入力の sha256・
    STOPPED_AT・停止の場・見張りの記録)。照合が例外で止まったときも結び付かないとする (合格にしない)。
    戻り: _v2_binding の戻り {pass, problems, verdict (停止の時の判定の中身), verdict_file, stop_step, criteria_sha256, ...}。"""
    from ..evaluate import stop_rule as SR
    from . import deltastar_loop as DL
    run_dir = Path(run_dir)
    try:
        return DL._v2_binding(run_dir, SR.read_state(run_dir), SR._settle())
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
        return {"pass": False, "problems": [f"照合が例外で止まった ({type(e).__name__}: {e})"], "verdict": None, "verdict_file": None,
                "stop_step": None, "criteria_sha256": None}


def bound_lines(verdict: dict | None) -> list[dict]:
    """停止の時の判定 (結び付いた記録) の基準の中身 (criteria.content) の line_mach と、その量の停止の時の verdict。
    今の settle_criteria.yaml を読み直した値は使わない (読み取りだけを変えても監視の線にならない; codex result 2 回目 M1)。"""
    if not isinstance(verdict, dict):
        return []
    content = (verdict.get("criteria") or {}).get("content") or {}
    by_name = {str(q.get("name")): q for q in verdict.get("quantities") or [] if isinstance(q, dict)}
    out = []
    for q in content.get("quantities") or []:
        if not (isinstance(q, dict) and q.get("kind") == "line_mach"):
            continue
        try:
            line = {"name": str(q.get("name")), "eta": float(q["eta"]), "x_range_m": [float(q["x_range_m"][0]), float(q["x_range_m"][1])],
                    "n": int(q.get("n", 2001)), "Md": q.get("Md")}
        except (KeyError, TypeError, ValueError, IndexError):
            continue
        line["verdict"] = (by_name.get(line["name"]) or {}).get("verdict")
        out.append(line)
    return out


def flow_check(run_dir) -> dict:
    """計算の成立: PASS / FAIL / INDETERMINATE / UNVERIFIED (判定は停止規則の記録だけを読む; 収束・関門の判定は judge_run・受け渡しの関門の役目)。
    v2 は停止の時の判定・基準・入力が今の run と結び付いた (`stop_binding`) ときだけ PASS・FAIL を決め、結び付かなければ INDETERMINATE。
    結び付いたときは停止の時の判定の監視の線 (`bound_lines`) を monitored_lines に入れる (V5c の線の保証が使う)。"""
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    st_path = run_dir / SR.STATE
    v2 = (run_dir / "settle_criteria.yaml").is_file()
    if not st_path.is_file():
        if v2:
            # 2026-10-17 (codex result 2 回目 M1) から、記録が無い v2 の run は FAIL でなく INDETERMINATE (確かめられない; 結果は同じ incomplete)
            return {"status": "INDETERMINATE", "method": "v2", "monitored_lines": [],
                    "why": f"settle_criteria.yaml のある run に {SR.STATE} が無い (本段が終わっていない・記録が失われた; 停止の時の判定を照合できない)"}
        return {"status": "UNVERIFIED", "method": None, "why": f"{SR.STATE} が無い (記録の写し・停止規則を使わない run) — 計算の成立を確かめられない"}
    s = json.loads(st_path.read_text())
    rec = {"method": s.get("method", "v1" if s.get("version") == 1 else None), "stop_status": s.get("status"), "stop_step": s.get("stop_step")}
    if v2 or s.get("method") == "v2":
        rv = (s.get("reeval") or {}).get("verdict")
        rec["reeval_verdict"] = rv
        b = stop_binding(run_dir)
        rec["binding"] = {"pass": bool(b.get("pass")), "problems": list(b.get("problems") or [])[:10], "n_problems": len(b.get("problems") or []),
                          "verdict_file": b.get("verdict_file"), "stop_step": b.get("stop_step"), "criteria_sha256": b.get("criteria_sha256"),
                          "tool": "design/forge_design/feedback/deltastar_loop.py _v2_binding (受け渡しの関門と同じ照合)"}
        if not rec["binding"]["pass"]:
            rec.update(status="INDETERMINATE", monitored_lines=[],
                       why=(f"停止の時の判定・基準・入力が今の run と結び付かない (状態 {s.get('status')}・再評価 {rv} を確かめられない): "
                            + "; ".join(rec["binding"]["problems"])[:1200]))
            return rec
        bv = (b.get("verdict") or {}).get("verdict")
        rec["bound_verdict"] = bv
        rec["monitored_lines"] = bound_lines(b.get("verdict"))
        ok = s.get("status") == "REACHED" and bv == "SETTLED"
        rec.update(status="PASS" if ok else "FAIL",
                   why=None if ok else f"停止規則 v2 の状態 {s.get('status')}・停止の step の再評価 {bv} (REACHED・SETTLED が要る)")
        return rec
    ok = s.get("status") == "REACHED"
    rec.update(status="PASS" if ok else "FAIL", why=None if ok else f"停止規則の状態 {s.get('status')} (REACHED が要る)")
    return rec


# ------------------------------------------------------------------ V5b・V5c
def v5b_exit_core(F: dict, Md: float, tol_rel: float, x_E: float) -> dict:
    """出口コア M (stop_rule.recorded_quantities と同じ定義) と許容の判定。"""
    from ..evaluate import stop_rule as SR
    M = float(SR.recorded_quantities(F["X"], F["R"], F["V"], float(x_E))["exit_core_M"])
    lo, hi = Md * (1.0 - tol_rel), Md * (1.0 + tol_rel)
    if not np.isfinite(M):
        return {"verdict": "INDETERMINATE", "M": M, "why": "非有限"}
    return {"verdict": "PASS" if lo <= M <= hi else "FAIL", "M": M, "dev_pct": 100.0 * (M / Md - 1.0), "Md": Md, "tol_rel": tol_rel,
            "range": [lo, hi], "eta_range": list(EXIT_ETA_RANGE),
            "definition": "最後の列の η ∈ [0.05, 0.7] の節点の M の平均 (stop_rule.recorded_quantities = nozzle_report.metrics の exit_core.M_mean)"}


def _axis_weight(F: dict, eta: float, xq) -> float:
    """eta_line の補間に入る軸の節点 (j = 0) の重みの最大 (settle.eta_line_nodes の axis_weight_max と同じ数え方: 使う列ごと)。"""
    X, R = F["X"], F["R"]
    xa = X[:, 0]
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2)
    wmax = 0.0
    for col in np.unique(np.concatenate([i, i + 1])):
        xp = R[col] / R[col, -1]
        k = int(np.searchsorted(xp, eta, side="right")) - 1
        if k == 0 and len(xp) > 1:
            wmax = max(wmax, float(1.0 - (eta - xp[0]) / (xp[1] - xp[0])))
    return wmax


def v5c_line(F: dict, Md: float, eta: float, x_range_m, n: int, tol_maxabs: float, tol_halfamp: float) -> dict:
    """近軸の線の一様性 (v5c_series.py・settle の line_mach と同じ式) と比較用の線。"""
    from ..report.nozzle_report import eta_line
    st = _settle()
    S = float(F["S"])
    x0, x1 = float(x_range_m[0]), float(x_range_m[1])
    xa = F["X"][:, 0]
    rec = {"eta": float(eta), "x_range_m": [x0, x1], "n": int(n), "Md": Md, "tol_maxabs_pct": tol_maxabs, "tol_halfamp_pct": tol_halfamp,
           "definition": "max|100(M/M_d − 1)|・100(max M − min M)/(2 M_d)、線は nozzle_report.eta_line (v5c_series.py と同じ)",
           "scope": "近軸の線 (η) の一様性であって試験部全体の一様性ではない (plan §6 V5c)"}
    if not (x0 < x1 and x0 / S >= xa[0] and x1 / S <= xa[-1]):
        rec.update(verdict="INDETERMINATE", why=f"窓 [{x0}, {x1}] m が格子の x [{xa[0] * S:.6g}, {xa[-1] * S:.6g}] m の外か空")
        return rec
    xq = np.linspace(x0 / S, x1 / S, int(n))
    M = eta_line(F, "M", float(eta), xq)
    wax = _axis_weight(F, float(eta), xq)
    rec["axis_weight_max"] = wax
    rec["axis_weight_limit"] = float(st.AXIS_WEIGHT_MAX)
    if not np.all(np.isfinite(M)):
        rec.update(verdict="INDETERMINATE", why="線の M に非有限がある")
        return rec
    q = st.line_quantities(M, Md)
    d = q["d"]
    rec.update(maxabs_pct=q["maxabs"], halfamp_pct=q["halfamp"], mean_pct=float(np.mean(d)),
               x_m_at_maxabs=float(xq[int(np.argmax(np.abs(d)))] * S), x_m_at_min=float(xq[int(np.argmin(d))] * S),
               x_m_at_max=float(xq[int(np.argmax(d))] * S))
    if wax > float(st.AXIS_WEIGHT_MAX):
        rec.update(verdict="INDETERMINATE", why=f"線 η {eta} の補間に軸の節点が重み {wax:.3g} で入る (> {st.AXIS_WEIGHT_MAX:g}; 第 1 層の η が線より上)")
    else:
        rec["checks"] = {"maxabs": q["maxabs"] <= tol_maxabs, "halfamp": q["halfamp"] <= tol_halfamp}
        rec["verdict"] = "PASS" if all(rec["checks"].values()) else "FAIL"
    cmp_ = {}
    for e in COMPARE_ETAS:
        Mc = eta_line(F, "M", float(e), xq)
        qc = st.line_quantities(Mc, Md)
        cmp_[f"eta{e:g}"] = {"maxabs_pct": qc["maxabs"], "halfamp_pct": qc["halfamp"], "mean_pct": float(np.mean(qc["d"]))}
    rec["compare"] = cmp_
    return rec


# ------------------------------------------------------------------ V5c の線の定常の保証 (§5.1 #7a)
def monitored_line(flow: dict, v5c: dict) -> dict | None:
    """停止の時の判定に結び付いた監視の線 (flow の monitored_lines) のうち、評価の線 (η・x 窓・点数) と同じもの。無ければ最初の線
    (違いを書くため)、線が無ければ None。"""
    lines = flow.get("monitored_lines") or []
    for ln in lines:
        if "eta" in v5c and _same_line(ln, v5c):
            return ln
    return lines[0] if lines else None


def _same_line(a: dict, b: dict) -> bool:
    return (abs(float(a["eta"]) - float(b["eta"])) <= 1e-12
            and all(abs(float(x) - float(y)) <= 1e-9 for x, y in zip(a["x_range_m"], b["x_range_m"]))
            and int(a["n"]) == int(b["n"]))


def line_monitor_gaps(v5c: dict, settle_line: dict | None, fld: dict, flow: dict) -> list[str]:
    """V5c の線と場が停止規則の監視で保証されていない理由 (空 = 保証される: 停止の時の判定・基準・入力が今の run と結び付き
    (flow の binding)、その判定の監視の line_mach (settle_line = `monitored_line`; 今のファイルを読み直した値ではない) と同じ線
    (η・x 窓・点数) で、その量が停止の時に SETTLED、場が停止の step)。2026-10-17 codex result M1: 違っても警告だけで accepted に
    できた (run_0014 に line_eta 0.06)。同 2 回目 M1: 停止の時の基準と結び付けず、今の settle_criteria.yaml の読み取り値と比べていた。"""
    why = []
    bnd = flow.get("binding")
    if not (isinstance(bnd, dict) and bnd.get("pass") is True):
        why.append("停止の時の判定・基準・入力が今の run と結び付かない (" + ("; ".join(bnd.get("problems") or [])[:400]
                                                                   if isinstance(bnd, dict) else "照合の記録が無い") + ")")
    if settle_line is None:
        why.append("監視の量 line_mach が停止の時の判定の基準に無い (停止規則が線を監視していない・結び付いた判定が無い)")
    elif "eta" in v5c:
        if not _same_line(settle_line, v5c):
            why.append(f"評価の線 (η {v5c['eta']}, x {v5c['x_range_m']} m, {v5c['n']} 点) が監視の line_mach (η {settle_line['eta']}, "
                       f"x {settle_line['x_range_m']} m, {settle_line['n']} 点) と違う")
        elif settle_line.get("verdict") != "SETTLED":
            why.append(f"停止の時の判定で監視の線 {settle_line.get('name')} が SETTLED でない ({settle_line.get('verdict')!r})")
    stop = flow.get("stop_step")
    if not isinstance(stop, int):
        why.append("停止の step の記録が無い (STOP_RULE_STATE.json の stop_step)")
    elif stop != int(fld["step"]):
        why.append(f"評価の場の step {fld['step']} が停止の step {stop} と違う")
    return why


def main_stage_fields(run_dir, end_step: int) -> dict:
    """本段の全場 res_<N>.h5 (本段の開始の記録 MAIN_STAGE_STARTED の start_step < N ≤ end_step、開始の時に残っていたファイルを除く)。
    戻り: {steps, start_step, why (作れないとき)}。"""
    run_dir = Path(run_dir)
    p = run_dir / "MAIN_STAGE_STARTED"
    if not p.is_file():
        return {"steps": [], "why": "本段の開始の記録 MAIN_STAGE_STARTED が無い (系列の区間を特定できない)"}
    try:
        ms = json.loads(p.read_text())
    except ValueError as e:
        return {"steps": [], "why": f"MAIN_STAGE_STARTED を読めない ({e})"}
    s0 = ms.get("start_step")
    if not isinstance(s0, int):
        return {"steps": [], "why": "MAIN_STAGE_STARTED に start_step が無い (系列の区間を特定できない)"}
    stale = {Path(str(f)).name for f in (ms.get("res_files_present") or [])}
    steps = [s for s in _res_steps(run_dir) if s0 < s <= int(end_step) and f"res_{s}.h5" not in stale]
    return {"steps": steps, "start_step": s0}


QS_COLS = ("maxabs_pct", "halfamp_pct")   # 系列で判定する量 = V5c が合否に使う 2 量
QS_MIN_SNAPS = 4                          # check_quasisteady.py の --min-snaps の既定と同じ (明示して渡す)


def line_series_steadiness(run_dir, fld: dict, v5c: dict, Md: float, out_dir=None) -> dict:
    """監視されていない V5c の線について、本段の全場の時系列 (評価の場で終わる) の最大偏差・片振幅を
    `check_quasisteady.py --series-csv` で判定する。status: STEADY / NOT_STEADY / INDETERMINATE (系列を作れない)。
    系列の CSV は out_dir (無ければ一時ディレクトリで、値は記録に埋め込む) に書く。"""
    import csv
    import hashlib
    import tempfile
    from ..report.nozzle_report import load_field
    run_dir = Path(run_dir)
    rec = {"tool": "solver_density_cuda/tools/check_quasisteady.py --series-csv", "cols": list(QS_COLS), "min_snaps": QS_MIN_SNAPS}
    ms = main_stage_fields(run_dir, int(fld["step"]))
    rec["steps"] = ms["steps"]
    if ms.get("why"):
        rec.update(status="INDETERMINATE", why=ms["why"])
        return rec
    if not ms["steps"] or ms["steps"][-1] != int(fld["step"]):
        rec.update(status="INDETERMINATE", why=f"本段の全場の系列が評価の場 (step {fld['step']}) で終わらない")
        return rec
    if len(ms["steps"]) < QS_MIN_SNAPS:
        rec.update(status="INDETERMINATE", why=f"本段の全場が {len(ms['steps'])} 枚 (< {QS_MIN_SNAPS}) で系列を作れない "
                                               f"(step {ms['start_step']} < N ≤ {fld['step']}: {ms['steps']})")
        return rec
    rows = []
    for s in ms["steps"]:
        r = v5c_line(load_field(run_dir, f"res_{s}.h5"), Md, float(v5c["eta"]), v5c["x_range_m"], int(v5c["n"]),
                     float(v5c["tol_maxabs_pct"]), float(v5c["tol_halfamp_pct"]))
        rows.append({"step": s, **{c: r.get(c, float("nan")) for c in QS_COLS}})
    rec["values"] = rows
    with tempfile.TemporaryDirectory() as td:
        od = Path(out_dir) if out_dir is not None else Path(td)
        csv_path = od / "v5c_line_series.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["step", *QS_COLS])
            w.writeheader()
            w.writerows(rows)
        cmd = [sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(csv_path), "--series-cols", ",".join(QS_COLS),
               "--min-snaps", str(QS_MIN_SNAPS)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        rec.update(csv=str(csv_path) if out_dir is not None else None, csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                   cmd=" ".join(cmd[1:]), returncode=r.returncode, output=(r.stdout + r.stderr)[-3000:])
    txt = r.stdout + r.stderr
    verd = {c: (m.group(1) if (m := re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", txt, re.M)) else None) for c in QS_COLS}
    rec["verdicts"] = verd
    ok = r.returncode == 0 and all(v == "STEADY" for v in verd.values())
    rec["status"] = "STEADY" if ok else "NOT_STEADY"
    if not ok:
        rec["why"] = f"check_quasisteady の判定 {verd} (終了コード {r.returncode})"
    return rec


def line_guarantee(v5c: dict, settle_line: dict | None, fld: dict, flow: dict, run_dir, Md: float, series_dir=None) -> dict:
    """V5c の線の定常の保証 (§5.1 #7a; v5c をその場で書き換えて返す): 停止の時の判定と結び付かない・監視の線と違う・監視の量が無い・
    その量が停止の時に SETTLED でない・停止の step でない、のどれかなら判定不能 (線だけの合否は line_verdict に残す)。ただし本段の
    全場の系列で check_quasisteady が STEADY を返したときだけ合否に進める (その VERDICT と系列のファイルを v5c["steadiness"] に記録)。
    settle_line は停止の時の判定に結び付いた監視の線 (`monitored_line`)。v5c が既に判定不能 (窓・軸の重み・非有限) なら何もしない。"""
    if v5c.get("verdict") not in ("PASS", "FAIL"):
        return v5c
    gaps = line_monitor_gaps(v5c, settle_line, fld, flow)
    if not gaps:
        v5c["steadiness"] = {"basis": "停止規則の監視 (停止の時の判定に結び付いた基準の line_mach と同じ線・その量が SETTLED、停止の step の場)",
                             "monitored_line": settle_line, "verdict_file": (flow.get("binding") or {}).get("verdict_file")}
        return v5c
    qs = line_series_steadiness(run_dir, fld, v5c, Md, series_dir)
    v5c["steadiness"] = {"basis": "本段の全場の系列の check_quasisteady (監視で保証されない)", "monitor_gaps": gaps, "series": qs}
    if qs["status"] != "STEADY":
        v5c["line_verdict"] = v5c["verdict"]
        v5c.update(verdict="INDETERMINATE",
                   why="線の定常を確かめていない: " + "; ".join(gaps) + f" — 全場の系列の準定常 {qs['status']} ({qs.get('why')})")
    return v5c


# ------------------------------------------------------------------ 壁解像・成果物
def wall_resolution(run_dir, step: int, out_txt=None, stop_step=None) -> dict:
    """`check_wall_resolution.py --step --groups <wall*> --weight area` を回して VERDICT を読む。PASS・FAIL でも、ツールが評価した
    壁と step (出力の EVALUATED 行) が要求した全部の壁 (bcond の名前・physID)・要求した step・停止の step (あれば) と一致しなければ
    判定不能にする (2026-10-17 codex result M2: ツールが最新の出力に代えた・存在しない壁を飛ばした PASS を受理していた)。"""
    run_dir = Path(run_dir)
    bc = yaml.safe_load((run_dir / "bcondConfig.yaml").read_text()) or {}
    walls = [k for k, v in bc.items() if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")]
    if not walls:
        return {"verdict": "INDETERMINATE", "why": "bcondConfig に no-slip 壁 (kind wall*) が無い"}
    want_pid = {k: (int(bc[k]["physID"]) if bc[k].get("physID") is not None else None) for k in walls}
    cmd = [sys.executable, str(TOOLS / "check_wall_resolution.py"), str(run_dir), "--step", str(int(step)), "--groups", ",".join(walls),
           "--weight", "area"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    txt = r.stdout + r.stderr
    if out_txt is not None:
        Path(out_txt).write_text(txt + f"\nCHECK_WALL_RESOLUTION_EXIT: {r.returncode}\n")
    m = re.search(r"^VERDICT: (\S+)(.*)$", txt, re.M)
    v = m.group(1) if m else "INDETERMINATE"
    out = {"verdict": v if v in ("PASS", "FAIL") else "INDETERMINATE", "verdict_line": m.group(0).strip() if m else None,
           "returncode": r.returncode, "groups": walls, "cmd": " ".join(cmd[1:]),
           "tool": "solver_density_cuda/tools/check_wall_resolution.py"}
    mm = re.search(r"目標 ([0-9.eE+-]+) を超える面積割合 最大 ([0-9.]+) % \(許容 ([0-9.eE+-]+) %\)", txt)
    if mm:
        out.update(target=float(mm.group(1)), over_area_pct=float(mm.group(2)), over_area_allow_pct=float(mm.group(3)))
    if out["verdict"] in ("PASS", "FAIL"):
        # 評価した壁と step の照合 (ツールの出力の EVALUATED 行; 書式 名前/physID@step)
        ev = re.search(r"^EVALUATED: (.*)$", txt, re.M)
        got = [(g.group(1), int(g.group(2)), int(g.group(3))) for g in re.finditer(r"([^\s,/]+)/(\d+)@(\d+)", ev.group(1))] if ev else []
        out["evaluated"] = [{"name": n, "physID": p, "step": s} for n, p, s in got]
        probs = []
        if not got:
            probs.append("ツールの出力に評価した壁と step (EVALUATED 行) が無い")
        names = [n for n, _, _ in got]
        miss = [w for w in walls if w not in names]
        extra = sorted({n for n in names if n not in walls})
        if miss:
            probs.append(f"要求した壁 {miss} を評価していない")
        if extra:
            probs.append(f"要求していない壁 {extra} を評価した")
        bad_pid = [f"{n} (physID {p}、bcond は {want_pid[n]})" for n, p, _ in got if n in want_pid and want_pid[n] is not None and p != want_pid[n]]
        if bad_pid:
            probs.append(f"bcond と違う physID の出力を評価した {bad_pid}")
        bad_step = sorted({s for _, _, s in got if s != int(step)})
        if bad_step:
            probs.append(f"要求した step {int(step)} でない出力 (step {bad_step}) を評価した")
        if isinstance(stop_step, int) and int(step) != stop_step:
            probs.append(f"評価した step {int(step)} が停止の step {stop_step} でない (V4 の判定ではない)")
        if probs:
            out.update(tool_verdict=out["verdict"], verdict="INDETERMINATE", why="; ".join(probs))
    if out["verdict"] == "INDETERMINATE" and not m:
        out["tail"] = txt[-1500:]
    return out


PROV_FIX = "この run から成果物を作り直す (python -m forge_design.export.deliverables <run> --card <カード> --problem <問題 YAML>)"


def _cmp_sha(name: str, recorded, current, mism: list, unchk: list, compared: dict) -> None:
    """記録と今の sha256 を比べて、違えば mism、どちらかが無ければ unchk に理由を足す (compared に両方を残す)。"""
    compared[name] = {"recorded": recorded, "current": current}
    if recorded is None or current is None:
        unchk.append(f"{name} の sha256 が" + ("生成元の記録" if recorded is None else "今の評価の対象") + "に無い (照合できない)")
    elif recorded != current:
        mism.append(f"{name} の sha256 が違う (成果物の生成元 {str(recorded)[:16]}…、今の評価の対象 {str(current)[:16]}…)")


def _find_key(o, key: str, path: str = "") -> list:
    """入れ子の dict・list の中の key の値を (場所, 値) で全部集める。"""
    out = []
    if isinstance(o, dict):
        for k, v in o.items():
            p = f"{path}.{k}" if path else str(k)
            if k == key:
                out.append((p, v))
            out += _find_key(v, key, p)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out += _find_key(v, key, f"{path}[{i}]")
    return out


MSH_DELIVERED = "mesh/nozzle.msh"          # part_mesh が配る msh (書き出しのときの run の nozzle.msh の写し)


def _msh_binding(man: dict, pv: dict, cur: dict, deliv_dir, mism: list, unchk: list, compared: dict) -> None:
    """配る msh と計算の格子 (h5) の元の msh の照合 (plan §5.1 #6c2 の再オープン、2026-10-17 codex result 4 回目 Major: 記録していた格子は
    h5 だけで、配る msh・Fluent・寸法の関門の入力 (run の nozzle.msh) と結ぶ検査が無かった)。
    基準は今の準備の記録 (`prepare_info.json` の `mesh.hashes.msh_sha256`; h5 は準備のときにこの msh から変換された)。基準と一致を要るもの:
    今の run の `nozzle.msh`、配る `mesh/nozzle.msh` の実ファイルと MANIFEST の記録 (`parts.mesh.msh41_sha256`・files_sha256; 記録か
    ファイルがあれば)、書き出しのときの msh (v2: 生成元の記録の `msh`・`prepare_msh`。v1 (msh の項目が無い古い形) は MANIFEST の
    `parts.mesh.msh41_sha256` か files_sha256 の `mesh/nozzle.msh` — どちらも無ければ書き出しのときの msh が分からず照合できない)。
    書き出しのときの一致の印 (MANIFEST の `msh_sha256_matches`・v2 の `msh.matches_prepare`) が false なら不一致、true でなければ照合できない。
    基準が無い・要るものが無い → unchk、違う → mism に理由を足す。"""
    from ..export import deliverables as DV
    pr = cur["prepare_msh"]
    ref = pr["sha256"]
    ref_name = f"計算の格子 ({cur['mesh']['file']}) の元の msh (準備の記録 {pr['file']} の {pr['key']})"
    v2 = pv.get("format") == DV.PROVENANCE_FORMAT
    parts = man.get("parts") if isinstance(man.get("parts"), dict) else {}
    pm = parts.get("mesh") if isinstance(parts.get("mesh"), dict) else {}
    fs = man.get("files_sha256") if isinstance(man.get("files_sha256"), dict) else {}
    items = [(f"今の run の {DV.MSH_FILE}", cur["msh"]["sha256"])]
    if v2:
        items += [(f"生成元の記録の msh (書き出しのときの run の {DV.MSH_FILE}; 配る msh・Fluent・寸法の関門の入力)", (pv.get("msh") or {}).get("sha256")),
                  (f"生成元の記録の準備の記録の msh (書き出しのときの {DV.PREPARE_FILE})", (pv.get("prepare_msh") or {}).get("sha256"))]
    elif "msh41_sha256" not in pm and MSH_DELIVERED not in fs:
        unchk.append("書き出しのときの msh の記録が無い (v1 の生成元の記録に msh の項目が無く、MANIFEST の parts.mesh.msh41_sha256・files_sha256 の "
                     f"{MSH_DELIVERED} も無い; 照合できない)")
    if "msh41_sha256" in pm:
        items.append(("配る msh の記録 (MANIFEST の parts.mesh.msh41_sha256)", pm.get("msh41_sha256")))
    if MSH_DELIVERED in fs:
        items.append((f"配る msh の記録 (MANIFEST の files_sha256 の {MSH_DELIVERED})", fs.get(MSH_DELIVERED)))
    dfile = None if deliv_dir is None else Path(deliv_dir) / MSH_DELIVERED
    if dfile is not None and (dfile.is_file() or "msh41_sha256" in pm or MSH_DELIVERED in fs):
        items.append((f"配る msh (成果物の {MSH_DELIVERED})", DV._sha_file(dfile)))
    if ref is None:
        unchk.append(f"準備の記録に msh の sha256 が無い ({pr.get('why')}; 計算の格子の元の msh が分からない — 照合できない)")
    for nm, s in items:
        if not isinstance(s, str) or not s:
            unchk.append(f"{nm} が無い (照合できない)")
        elif ref is not None and s != ref:
            mism.append(f"{nm} の sha256 {s[:16]}… が{ref_name} {ref[:16]}… と違う")
    flags = [(f"MANIFEST の parts.mesh.{p}", v) for p, v in _find_key(pm, "msh_sha256_matches")]
    if v2:
        flags.append(("生成元の記録の msh.matches_prepare", (pv.get("msh") or {}).get("matches_prepare")))
    for nm, v in flags:
        if v is False:
            mism.append(f"書き出しのときに msh が準備の記録と違った ({nm} = false)")
        elif v is not True:
            unchk.append(f"書き出しのときの msh と準備の記録の一致の印が true でない ({nm} = {v!r}; 照合できない)")
    compared["msh"] = {"reference": ref, "reference_source": ref_name, "items": {nm: s for nm, s in items},
                       "flags": {nm: v for nm, v in flags}, "provenance_format": pv.get("format")}


def _h5_binding(run_dir: Path, cur: dict, mism: list, unchk: list, compared: dict) -> None:
    """計算の格子 (h5) が準備のときに msh から変換したものと同じかの照合 (2026-10-17、result 段レビュー 4 回目の後に implementer が見つけた穴:
    準備の後に h5 だけを別の msh から変換し直すと、msh・準備の記録・生成元の記録 (h5 の sha256) がそろったまま通る)。
    準備の記録 `prepare_info.json` の `mesh.hashes.coord_sha256`・`topology_sha256` (変換の直後の h5 の座標と接続) と、今の h5 から同じ式
    (`runner_axismach.h5_mesh_hashes`) で計算し直した値が一致することを要る。記録が無い・h5 を読めない → 照合できない。"""
    from ..evaluate.runner_axismach import h5_mesh_hashes
    h5 = Path(run_dir) / cur["mesh"]["file"]
    try:
        hashes = (json.loads((Path(run_dir) / "prepare_info.json").read_text()).get("mesh") or {}).get("hashes") or {}
    except Exception as e:  # noqa: BLE001
        unchk.append(f"準備の記録 prepare_info.json を読めない ({e}; 計算の格子の座標・接続を照合できない)")
        return
    rc, rt = hashes.get("coord_sha256"), hashes.get("topology_sha256")
    if not (isinstance(rc, str) and isinstance(rt, str)):
        unchk.append("準備の記録に計算の格子の座標・接続の sha256 (mesh.hashes.coord_sha256・topology_sha256) が無い (照合できない)")
        return
    try:
        now = h5_mesh_hashes(h5)
    except Exception as e:  # noqa: BLE001
        unchk.append(f"計算の格子 {h5.name} を読めない ({e}; 照合できない)")
        return
    compared["h5"] = {"coord": {"recorded": rc, "current": now["coord_sha256"]}, "topology": {"recorded": rt, "current": now["topology_sha256"]}}
    if now["coord_sha256"] != rc:
        mism.append(f"計算の格子 {h5.name} の座標が準備の記録と違う (変換の後に h5 が作り直された; {now['coord_sha256'][:16]}… と {rc[:16]}…)")
    if now["topology_sha256"] != rt:
        mism.append(f"計算の格子 {h5.name} の接続が準備の記録と違う ({now['topology_sha256'][:16]}… と {rt[:16]}…)")


def provenance_check(run_dir, man: dict, field: dict | None, card_path=None, problem_path=None, deliv_dir=None) -> dict:
    """成果物の生成元の記録 (MANIFEST.json の `provenance`) と今の評価の対象の照合 (plan §5.1 #6b の続き、2026-10-17 codex result
    3 回目 Major)。比べるのは内容 (sha256) と step だけ — 壁・格子・設定 (`deliverables.run_identity`)、評価の場 (`field` の step と
    `res_<step>.h5` の sha256、生成元が読んだ壁ダンプ)、停止の step (`deliverables.stop_identity`)、カード・問題 YAML (両方に無ければ
    比べるものが無いので一致、片方だけなら照合できない)、配る msh と計算の格子の元の msh (`_msh_binding`; §5.1 #6c2 の再オープン、
    4 回目 Major。deliv_dir は配る `mesh/nozzle.msh` の実ファイルを読むのに使う)。書式は v2 と v1 (msh の項目が無い古い形) を読む。
    status: MATCH / MISMATCH (どれかが違う) / ABSENT (生成元の記録が無い古い成果物) / UNCHECKABLE (照合できない)。"""
    from ..export import deliverables as DV
    run_dir = Path(run_dir)
    pv = man.get("provenance") if isinstance(man, dict) else None
    rec = {"tool": "design/forge_design/feedback/wt_qualify.py provenance_check (deliverables.run_identity・field_identity・stop_identity と同じ定義)"}
    if not isinstance(pv, dict):
        rec.update(status="ABSENT", problems=["成果物に生成元の記録 (MANIFEST.json の provenance) が無い (古い成果物; 生成元の run・場を確かめられない)"])
        return rec
    if pv.get("format") not in (DV.PROVENANCE_FORMAT, DV.PROVENANCE_FORMAT_V1):
        rec.update(status="UNCHECKABLE", problems=[f"生成元の記録の書式 {pv.get('format')!r} が {DV.PROVENANCE_FORMAT} か "
                                                   f"{DV.PROVENANCE_FORMAT_V1} でない (照合できない)"])
        return rec
    mism, unchk, cmp_ = [], [], {}
    cur = DV.run_identity(run_dir)
    _cmp_sha(f"壁 ({cur['wall']['file']})", (pv.get("wall") or {}).get("sha256"), cur["wall"]["sha256"], mism, unchk, cmp_)
    _cmp_sha(f"格子 ({cur['mesh']['file']})", (pv.get("mesh") or {}).get("sha256"), cur["mesh"]["sha256"], mism, unchk, cmp_)
    for n in DV.PROVENANCE_CONFIGS:
        _cmp_sha(f"設定 ({n})", (pv.get("config") or {}).get(n), cur["config"][n], mism, unchk, cmp_)
    # 配る msh と計算の格子の元の msh (§5.1 #6c2 の再オープン、2026-10-17 codex result 4 回目 Major)
    _msh_binding(man, pv, cur, deliv_dir, mism, unchk, cmp_)
    _h5_binding(run_dir, cur, mism, unchk, cmp_)
    pf = pv.get("field")
    if not isinstance(field, dict) or not isinstance(field.get("step"), int):
        unchk.append("評価の場の step が渡されていない (照合できない)")
    elif not isinstance(pf, dict) or not isinstance(pf.get("step"), int):
        unchk.append(f"生成元の場の記録が無い ({pv.get('field_why') or '理由の記録なし'}; 照合できない)")
    else:
        cmp_["場の step"] = {"recorded": pf["step"], "current": int(field["step"])}
        if pf["step"] != int(field["step"]):
            mism.append(f"場の step が違う (成果物の生成元 {pf['step']}、評価の場 {int(field['step'])})")
        else:
            fi = DV.field_identity(run_dir, int(field["step"]))
            _cmp_sha(f"場 ({fi['file']})", pf.get("sha256"), fi["sha256"], mism, unchk, cmp_)
        wd = pf.get("wall_dumps")
        if not isinstance(wd, dict) or not wd:
            unchk.append("生成元が読んだ壁ダンプの記録が無い (照合できない)")
        else:
            for n, s in sorted(wd.items()):
                _cmp_sha(f"壁ダンプ ({n})", s, DV._sha_file(run_dir / n), mism, unchk, cmp_)
    ps = (pv.get("stop_rule") or {}).get("stop_step") if isinstance(pv.get("stop_rule"), dict) else None
    cs = (DV.stop_identity(run_dir) or {}).get("stop_step")
    cmp_["停止の step"] = {"recorded": ps, "current": cs}
    if ps != cs:
        mism.append(f"停止の step が違う (成果物の生成元 {ps}、今の run {cs})")
    for nm, rp, p in (("カード", pv.get("card"), card_path), ("問題 YAML", pv.get("problem"), problem_path)):
        rs = rp.get("sha256") if isinstance(rp, dict) else None
        cs_ = None if p is None else DV._sha_file(Path(p))
        if rp is None and p is None:
            cmp_[nm] = {"recorded": None, "current": None, "note": "書き出し・評価のどちらにも渡していない"}
            continue
        _cmp_sha(nm, rs, cs_, mism, unchk, cmp_)
    rec["compared"] = cmp_
    rec["recorded_run_name"] = pv.get("run_name")
    rec["status"] = "MISMATCH" if mism else "UNCHECKABLE" if unchk else "MATCH"
    rec["problems"] = mism + unchk
    return rec


def deliverables_check(run_dir, deliv_dir=None, card_path=None, optional_parts=(), field=None, problem_path=None, gate=None) -> dict:
    """成果物 (MANIFEST.json) と寸法の関門 (DIMENSIONS.json)。status: PASS / FAIL / INDETERMINATE / ABSENT。
    field (評価の場 {step, ...})・problem_path (評価に渡した問題 YAML) は生成元の照合 (`provenance_check`) に使う。照合が MATCH で
    なければ成果物と寸法の関門を INDETERMINATE にする (2026-10-17 codex result 3 回目 Major)。"""
    from ..export.deliverables import PARTS, PASS_REQUIRED, REQUIRED_FILES
    run_dir = Path(run_dir)
    od = Path(deliv_dir) if deliv_dir is not None else run_dir / "deliverables"
    rec = {"dir": str(od), "required_parts": [p for p in PARTS if p not in optional_parts]}
    if not od.is_dir():
        rec.update(status="ABSENT", dimensions={"status": "ABSENT"}, why=f"{od} が無い (成果物を作っていない)")
        return rec
    probs = []
    try:
        man = json.loads((od / "MANIFEST.json").read_text())
        dim = json.loads((od / "DIMENSIONS.json").read_text())
    except (OSError, ValueError) as e:
        rec.update(status="INDETERMINATE", dimensions={"status": "INDETERMINATE"}, why=f"MANIFEST.json・DIMENSIONS.json を読めない ({e})")
        return rec
    # 寸法の関門
    dv = dim.get("verdict")
    dprobs = []
    if dv != "PASS":
        dprobs.append(f"DIMENSIONS.json の verdict {dv!r}")
    if card_path is not None:
        sha = CC._sha256(card_path)
        if dim.get("card_sha256") != sha:
            dprobs.append(f"寸法の関門のカード ({dim.get('card')}, sha256 {str(dim.get('card_sha256'))[:16]}…) が渡したカード "
                          f"({Path(card_path).name}, {sha[:16]}…) と違う")
    rec["dimensions"] = {"status": "PASS" if not dprobs else "FAIL", "verdict": dv, "items": {k: (v or {}).get("pass") for k, v in
                                                                                          (dim.get("items") or {}).items()},
                         "card": dim.get("card"), "card_sha256": dim.get("card_sha256"), "problems": dprobs}
    # 成果物
    if man.get("verdict") != "PASS":
        probs.append(f"MANIFEST.json の verdict {man.get('verdict')!r}")
    if man.get("failed_parts"):
        probs.append(f"失敗した部品 {man.get('failed_parts')}")
    parts = man.get("parts") or {}
    miss = [p for p in rec["required_parts"] if not isinstance(parts.get(p), dict) or parts[p].get("ok") is not True]
    if miss:
        probs.append(f"必要な部品が欠けた・失敗した {miss} (skipped {man.get('skipped')}, errors {list((man.get('errors') or {}).keys())})")
    # 部品の検査の総合と MANIFEST の総合の矛盾 (plan tooling-design-orchestrator §4.7.1 (6); codex 2026-10-17: cad.pass=false でも PASS だった)
    nopass = [p for p in rec["required_parts"] if p in PASS_REQUIRED and isinstance(parts.get(p), dict) and parts[p].get("pass") is not True]
    failed = [p for p, v in parts.items() if isinstance(v, dict) and (v.get("ok") is False or v.get("pass") is False)]
    if nopass:
        probs.append(f"部品の検査の総合 (pass) が true でない {[(p, parts[p].get('pass')) for p in nopass]} (MANIFEST の verdict {man.get('verdict')!r})")
    if failed and man.get("verdict") == "PASS":
        probs.append(f"MANIFEST の verdict PASS と部品の不合格 {failed} が矛盾する")
    if dv != "PASS" and man.get("verdict") == "PASS":
        probs.append(f"MANIFEST の verdict PASS と寸法の関門 {dv!r} が矛盾する")
    # 必須のファイルは MANIFEST の列挙 (自己申告) でなく、この版のコードの一覧 (REQUIRED_FILES) で見る
    fl = man.get("files_sha256") or {}
    req_missing = [f for p in rec["required_parts"] for f in REQUIRED_FILES.get(p, ()) if f not in fl or not (od / f).is_file()]
    if req_missing:
        probs.append(f"必須のファイルが無い・MANIFEST の列挙に無い {req_missing[:10]}" + (f" ほか {len(req_missing) - 10}" if len(req_missing) > 10 else ""))
    rec["required_files"] = {p: list(REQUIRED_FILES.get(p, ())) for p in rec["required_parts"]}
    bad_files = []
    for rel, sha in (man.get("files_sha256") or {}).items():
        p = od / rel
        if not p.is_file():
            bad_files.append(f"{rel} (無い)")
        elif CC._sha256(p) != sha:
            bad_files.append(f"{rel} (sha256 が違う)")
    if bad_files:
        probs.append(f"MANIFEST の列挙と違うファイル {bad_files[:10]}" + (f" ほか {len(bad_files) - 10}" if len(bad_files) > 10 else ""))
    if not man.get("files_sha256"):
        probs.append("MANIFEST.json に files_sha256 が無い")
    warns = []
    if man.get("run") and Path(str(man["run"])).name != run_dir.name:
        # 名前は警告だけ (合否は下の生成元の照合 = 内容で決める; 別の置き場所への写しは名前が違っても中身が同じなら通す)
        warns.append(f"MANIFEST の run {man.get('run')} がこの run ({run_dir.name}) と名前が違う")
    rec.update(status="PASS" if not probs and not dprobs else "FAIL", problems=probs, warnings=warns,
               manifest={"verdict": man.get("verdict"), "created": man.get("created"), "git_commit": man.get("git_commit"),
                         "n_files": len(man.get("files_sha256") or {})})
    # 生成元の照合 (2026-10-17 codex result 3 回目 Major): 一致しなければ、その成果物の中身 (ファイルの sha256・寸法の関門) が正しくても
    # この run の成果物・寸法の判定にしない (寸法の FAIL を補正 (寸法の解き直し) にも回さない)
    pv = provenance_check(run_dir, man, field, card_path, problem_path, deliv_dir=od)
    rec["provenance"] = pv
    if pv["status"] != "MATCH":
        what = {"ABSENT": "の記録が無い (古い成果物)", "MISMATCH": "が今の評価の対象と一致しない",
                "UNCHECKABLE": "を今の評価の対象と照合できない"}.get(pv["status"], f" {pv['status']}")
        rec["problems"] = probs + [f"生成元: {p}" for p in pv["problems"]]
        rec["dimensions"].update(status="INDETERMINATE", status_before_provenance=rec["dimensions"]["status"],
                                 why=f"成果物の生成元{what} — この寸法の関門はこの run の判定にしない")
        rec.update(status="INDETERMINATE", status_before_provenance=rec["status"],
                   why=f"成果物の生成元{what} ({pv['status']}): " + "; ".join(pv["problems"])[:800] + f" — {PROV_FIX}")
    if gate is not None:
        # 正式な入口 (plan tooling-design-orchestrator §4.7.1 (6)): 成果物の部品の検査は、判定専用の環境で作ったもの (MANIFEST の gate が
        # 承認の記録と合う) だけを信じる。合わなければ部品の自己申告を判定に使わない
        g = man.get("gate") or {}
        gp = [f"{k} {str(g.get(k))[:16]}… ≠ 承認 {str(v)[:16]}…" for k, v in gate.items() if g.get(k) != v]
        rec["gate"] = {"manifest": g, "expected": gate, "problems": gp}
        if gp and rec["status"] != "ABSENT":
            rec["dimensions"].update(status="INDETERMINATE", why="成果物を判定専用の環境で作っていない")
            rec.update(status="INDETERMINATE", status_before_gate=rec["status"],
                       why="成果物を判定専用の環境で作っていない (MANIFEST の gate: " + "; ".join(gp) +
                           ") — python -m forge_design.feedback.gate_env deliverables --approval … で作り直す")
    return rec


def approval_rechecks(run_dir, flow: dict, apr: dict) -> dict:
    """正式な入口で、run の側で作られた判定 (停止の判定・計算の格子の品質) を判定専用の環境で回し直す (plan tooling-design-orchestrator
    §4.7.1 (6))。run は作業側のコードで回るので、その判定の記録 (SETTLE_VERDICT・MESH_QUALITY) の自己申告を使わない。"""
    out, probs = {}, []
    stop = flow.get("stop_step")
    if isinstance(stop, int) and not isinstance(stop, bool):
        cmd = [sys.executable, str(TOOLS / "check_settled.py"), str(run_dir), "--end-step", str(stop), "--no-write"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        vl = next((ln for ln in reversed(r.stdout.splitlines()) if ln.startswith("VERDICT:")), None)
        out["settle"] = {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdict_line": vl}
        if not (vl or "").startswith("VERDICT: SETTLED"):
            probs.append(f"停止の判定を固定した版で回し直すと SETTLED でない ({vl or r.stderr.strip()[-200:]}; 停止の step {stop})")
    else:
        probs.append("停止の step が無い (停止の判定を固定した版で回し直せない)")
    mq = (apr["protected"].get("evaluation") or {}).get("mesh_quality")
    h5 = Path(run_dir) / "nozzle.h5"
    if not mq:
        probs.append("承認の記録にメッシュ品質の閾値が無い")
    elif not h5.is_file():
        probs.append(f"計算の格子 {h5.name} が無い (品質を回し直せない)")
    else:
        cmd = [sys.executable, str(TOOLS / "check_mesh_quality.py"), str(h5), "--ar-max", f"{mq['ar_max']:g}", "--skew-max", f"{mq['skew_max']:g}"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        vl = next((ln for ln in reversed(r.stdout.splitlines()) if ln.startswith("VERDICT:")), None)
        out["mesh_quality"] = {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdict_line": vl}
        if not (vl or "").startswith("VERDICT: PASS"):
            probs.append(f"計算の格子の品質を固定した版で回し直すと PASS でない ({vl or r.stderr.strip()[-200:]})")
    out["problems"] = probs
    return out


# ------------------------------------------------------------------ 較正し直しの規則
def prev_error(prev_run, prev_delta_next=None, x_far_m: float = CC.X_FAR_M) -> dict:
    """前の run の max|δ_in − δ_E| (x ≥ 0; δ_in は表の delta_in_prev = 格子で測った値)。"""
    prev_run = Path(prev_run)
    loc = CC.locate_delta_e(prev_run, prev_delta_next)
    t = CC.read_delta_next(loc["path"])
    S = float(json.loads((prev_run / "prepare_info.json").read_text())["scale_m"])
    x_m = np.asarray(t["x_rt"], dtype=float) * S
    use = x_m >= 0.0
    if np.any(np.asarray(t["held"])[use] != 0) or not np.isfinite(np.asarray(t["delta_E"])[use]).all():
        raise ValueError(f"前の run の δ_E の表 {loc['path']} の x ≥ 0 に held・非有限がある")
    e = CC.error_summary(x_m, np.asarray(t["delta_in_prev"]) * S, np.asarray(t["delta_E"]) * S, x_far_m)
    return {"run": str(prev_run.resolve()), "delta_E": loc, "scale_m": S, "error": e}


def recalibration_rule(run_dir, v5b: dict, v5c: dict, flow: dict, *, delta_next=None, prev_run=None, prev_delta_next=None,
                       recalibrations_done: int = 0, max_recalibrations: int = MAX_RECALIBRATIONS,
                       min_reduction: float = RECAL_MIN_REDUCTION, x_far_m: float = CC.X_FAR_M) -> dict:
    """§4.3 の較正し直しの規則 (a)〜(e)。戻り: {apply, conditions, calibration (contur_calibrate の記録 or None), why}。"""
    cond = {}
    cond["a_both_judged"] = {"ok": v5b.get("verdict") in ("PASS", "FAIL") and v5c.get("verdict") in ("PASS", "FAIL"),
                             "v5b": v5b.get("verdict"), "v5c": v5c.get("verdict")}
    cal, cal_err = None, None
    try:
        cal = CC.calibrate_run(run_dir, delta_next=delta_next, x_far_m=x_far_m)
    except (ValueError, OSError, KeyError) as e:
        cal_err = str(e)
    cur = None if cal is None else cal["errors"]["current_wall"]["max_abs_mm"]
    proj = None if cal is None else cal["errors"]["projected_after_calibration"]["max_abs_mm"]
    if prev_run is None:
        cond["b_prev_update_reduced_error"] = {"ok": False, "why": "--prev-run が無い (前の更新で誤差が縮んだかを確かめられない)", "current_max_abs_mm": cur}
    else:
        try:
            pe = prev_error(prev_run, prev_delta_next, x_far_m)
            pv = pe["error"]["max_abs_mm"]
            cond["b_prev_update_reduced_error"] = {"ok": cur is not None and cur < pv, "prev_max_abs_mm": pv, "current_max_abs_mm": cur,
                                                   "prev_run": pe["run"], "prev_delta_E": pe["delta_E"]["path"],
                                                   "prev_x_m_at_max": pe["error"]["x_m_at_max"]}
        except (ValueError, OSError, KeyError) as e:
            cond["b_prev_update_reduced_error"] = {"ok": False, "why": f"前の run の誤差を出せない ({e})", "current_max_abs_mm": cur}
    if cal is None:
        cond["c_further_reducible"] = {"ok": False, "why": f"この run で較正し直しの射影を作れない ({cal_err})"}
    else:
        red = 1.0 - proj / cur if cur > 0 else 0.0
        k = f"max_abs_x_ge_{x_far_m:g}m_mm"
        cond["c_further_reducible"] = {"ok": proj <= (1.0 - min_reduction) * cur, "current_max_abs_mm": cur, "projected_max_abs_mm": proj,
                                       "reduction": red, "min_reduction": min_reduction,
                                       f"current_{k}": cal["errors"]["current_wall"][k],
                                       f"projected_{k}": cal["errors"]["projected_after_calibration"][k],
                                       "note": "min_reduction は plan に数値が無い暫定の値 (--recal-min-reduction)"}
    cond["d_count"] = {"ok": int(recalibrations_done) < int(max_recalibrations), "recalibrations_done": int(recalibrations_done),
                       "max_recalibrations": int(max_recalibrations)}
    cond["e_flow_verified"] = {"ok": flow.get("status") == "PASS", "flow": flow.get("status")}
    apply = all(c["ok"] for c in cond.values())
    why = [f"{k}: " + (c.get("why") or json.dumps({kk: vv for kk, vv in c.items() if kk != "ok"}, ensure_ascii=False, default=str))
           for k, c in cond.items() if not c["ok"]]
    return {"apply": apply, "conditions": cond, "calibration": cal, "why_not": why,
            "rule": "plan §4.3 の末尾 (codex diagnose 2026-10-15 m55-three-coef M2): V5b・V5c の両方を判定・前の更新で誤差が縮んだ・"
                    "同じ 3 係数でさらに減らせる・1 回まで (+ 計算の成立を確かめた)"}


# ------------------------------------------------------------------ 判定
def decide(flow: dict, v5b: dict, v5c: dict, wr: dict, dl: dict, *, iteration=None, max_iterations=None, rule_fn=None) -> dict:
    """各判定から結果と次の手を決める (ファイルを読まない純粋な部分; 規則の評価だけ rule_fn() を呼ぶ)。
    順: 判定不能・計算の成立の FAIL・INDETERMINATE (停止の時の判定と結び付かない) → incomplete (escalate)。公差・壁解像・寸法の FAIL → 反復の上限なら rejected、でなければ
    back_to_correction (壁の形だけなら較正し直しの規則 → recalibrate / escalate、壁解像・寸法を含めば escalate)。すべて合格で成果物が
    無い・不完全、または計算の成立が UNVERIFIED → incomplete (成果物だけが無く成立を確かめていれば make_deliverables)。残りは accepted。
    戻り: {result, next, missing, corrections, rule, iteration}。"""
    corrections, missing = [], []
    for nm, c in (("V5b", v5b), ("V5c", v5c), ("wall_resolution", wr)):
        if c.get("verdict") not in ("PASS", "FAIL"):
            missing.append(f"{nm} が判定不能 ({c.get('why') or c.get('verdict_line') or c.get('verdict')})")
    if flow.get("status") in ("FAIL", "INDETERMINATE"):     # INDETERMINATE: 停止の時の判定と結び付かない (§5.1 #7a の続き)
        missing.append(f"計算の成立 ({flow.get('status')}): {flow.get('why')}")
    if v5b.get("verdict") == "FAIL" or v5c.get("verdict") == "FAIL":
        corrections.append({"kind": "wall_shape", "by": [n for n, c in (("V5b", v5b), ("V5c", v5c)) if c.get("verdict") == "FAIL"],
                            "procedure": "plan §4.3 の較正し直し (contur_calibrate)"})
    if wr.get("verdict") == "FAIL":
        corrections.append({"kind": "mesh_first_layer", "by": ["wall_resolution"], "procedure": "plan §4.2 (NS の壁量から第一層の表を作り直す)"})
    dim_status = (dl.get("dimensions") or {}).get("status")
    if dim_status == "FAIL":
        corrections.append({"kind": "dimensions", "by": ["dimension_gate"], "procedure": "寸法の解き直し (カード・r_t)"})
    rule = None
    it = {"iteration": iteration, "max_iterations": max_iterations, "evaluated": iteration is not None and max_iterations is not None}
    if not it["evaluated"]:
        it["note"] = "反復の上限は判定していない (--iteration・--max-iterations が無い)"
    if missing:
        result, nxt = "incomplete", "escalate"
    elif corrections:
        if it["evaluated"] and int(iteration) >= int(max_iterations):
            result, nxt = "rejected", "escalate"
        else:
            result = "back_to_correction"
            if {c["kind"] for c in corrections} == {"wall_shape"} and rule_fn is not None:
                rule = rule_fn()
                nxt = "recalibrate" if rule["apply"] else "escalate"
            else:
                nxt = "escalate"
    elif dl.get("status") != "PASS" or dim_status != "PASS":
        result = "incomplete"
        if dl.get("status") == "ABSENT":
            missing.append("成果物 (deliverables) が無い — 流れの判定はすべて合格")
            nxt = "make_deliverables" if flow.get("status") == "PASS" else "escalate"
        elif dl.get("status") == "INDETERMINATE" and dl.get("why"):
            # 読めない・生成元が一致しない・記録が無い・照合できない (2026-10-17 codex result 3 回目 Major): 次の手は今までの規則
            # (ABSENT 以外は escalate) のまま、理由 (why) に作り直しを書く
            missing.append(f"成果物をこの run の判定に使えない ({dl['why']}; 寸法 {dim_status})")
            nxt = "escalate"
        else:
            missing.append(f"成果物が不完全 ({dl.get('problems') or dl.get('why')}; 寸法 {dim_status})")
            nxt = "escalate"
        if flow.get("status") != "PASS":
            missing.append(f"計算の成立を確かめていない ({flow.get('why')})")
    elif flow.get("status") != "PASS":
        result, nxt = "incomplete", "escalate"
        missing.append(f"計算の成立を確かめていない ({flow.get('why')})")
    else:
        result, nxt = "accepted", "done"
    return {"result": result, "next": nxt, "missing": missing, "corrections": corrections, "rule": rule, "iteration": it}


def qualify(run_dir, card=None, problem=None, *, prev_run=None, delta_next=None, prev_delta_next=None, step=None, deliverables=None,
            iteration=None, max_iterations=None, recalibrations_done: int = 0, max_recalibrations: int | None = None,
            recal_min_reduction: float = RECAL_MIN_REDUCTION, optional_parts=(), cli_tols: dict | None = None,
            x_far_m: float = CC.X_FAR_M, wall_txt=None, series_dir=None, approval=None) -> dict:
    """1 本の NS run を判定して記録の dict を返す (ファイルは書かない; 壁解像の出力だけ wall_txt に、監視されていない V5c の線の
    全場の系列だけ series_dir (無ければ一時ディレクトリ) に書く)。"""
    from ..report.nozzle_report import load_field
    run_dir = Path(run_dir).resolve()
    info = json.loads((run_dir / "prepare_info.json").read_text())
    card_doc = _load_yaml(card)
    prob_doc = _load_yaml(problem)
    warns = []
    # M_d
    Md_info = info.get("Md")
    Md_card = None if card_doc is None else card_doc.get("M_exit")
    if Md_card is not None and Md_info is not None and abs(float(Md_card) - float(Md_info)) > 1e-12:
        raise ValueError(f"カードの M_exit {Md_card!r} と prepare_info の Md {Md_info!r} が違う (別の設計の run とカード)")
    if Md_card is None and Md_info is None:
        raise ValueError("M_d が無い (カードの M_exit・prepare_info の Md)")
    Md = float(Md_card if Md_card is not None else Md_info)
    if prob_doc is not None:
        rt = ((prob_doc.get("spec") or {}).get("r_throat")) if isinstance(prob_doc, dict) else None
        if rt is not None and abs(float(rt) / float(info["scale_m"]) - 1.0) > 1e-12:
            warns.append(f"問題 YAML の spec.r_throat {rt!r} が run の scale_m {info['scale_m']!r} と違う (別の設計の問題 YAML?)")
    apr, gate_run = None, None
    if approval is not None:
        # 正式な入口 (plan tooling-design-orchestrator §4.7.1): 判定の条件は承認の記録の値だけ。CLI・問題 YAML の相違は拒否する
        apr, max_iterations, max_recalibrations, optional_parts = approval_mode(
            approval, card, problem, cli_tols or {}, optional_parts, iteration, max_iterations, max_recalibrations)
        from . import gate_env as GE
        gate_run = GE.assert_in_gate_env(apr)     # 判定専用の環境 (plan §4.7.1 (5)) で回っているか
        tol = approved_tolerances(run_dir, info, apr, approval)
    else:
        if max_recalibrations is None:
            max_recalibrations = MAX_RECALIBRATIONS
        cq = _qual_block(card_doc, "カード")
        if not optional_parts and cq.get("optional_parts"):
            optional_parts = tuple(cq["optional_parts"])
            warns.append(f"無くてもよい成果物はカードの qualification.optional_parts {list(optional_parts)} (CLI の指定なし)")
        tol = resolve_tolerances(run_dir, info, card_doc, prob_doc, cli_tols or {})
    fld = select_field(run_dir, step)
    flow = flow_check(run_dir)
    rechecks = approval_rechecks(run_dir, flow, apr) if apr is not None else None
    F = load_field(run_dir, fld["res"])
    v5b = v5b_exit_core(F, Md, float(tol["exit_core_M_tol_rel"]["value"]), float(info["x_E"]))
    if tol["line_x_range_m"]["value"] is None:
        v5c = {"verdict": "INDETERMINATE", "why": tol["line_x_range_m"]["source"]}
    else:
        v5c = v5c_line(F, Md, float(tol["line_eta"]["value"]), tol["line_x_range_m"]["value"], int(tol["line_n"]["value"]),
                       float(tol["line_tol_maxabs_pct"]["value"]), float(tol["line_tol_halfamp_pct"]["value"]))
    sl = tol.get("settle_line")
    if sl is not None and "eta" in v5c:
        same = abs(sl["eta"] - v5c["eta"]) <= 1e-12 and all(abs(a - b) <= 1e-9 for a, b in zip(sl["x_range_m"], v5c["x_range_m"]))
        v5c["matches_settle_line"] = same
        if not same:
            warns.append(f"V5c の線 (η {v5c['eta']}, {v5c['x_range_m']} m) が監視の line_mach (η {sl['eta']}, {sl['x_range_m']} m) と違う — "
                         "線の定常は停止規則の判定で保証されていない")
    if isinstance(flow.get("stop_step"), int) and flow["stop_step"] != fld["step"]:
        warns.append(f"判定した場 (step {fld['step']}) が停止規則の止めた場 (step {flow['stop_step']}) でない")
    # 線の保証は停止の時の判定に結び付いた監視の線で見る (今の settle_criteria.yaml の読み取り値 sl は出所と警告にだけ使う;
    # 2026-10-17 codex result 2 回目 M1)
    line_guarantee(v5c, monitored_line(flow, v5c), fld, flow, run_dir, Md, series_dir)
    wr = wall_resolution(run_dir, fld["step"], wall_txt, stop_step=flow.get("stop_step"))
    gate_bind = None
    if apr is not None:
        from ..intake import approval as _AP
        gate_bind = {"gate_sha256": apr.get("gate_sha256"), "approval_sha256": _AP.sha256_file(approval)}
    dl = deliverables_check(run_dir, deliverables, card, optional_parts, field=fld, problem_path=problem, gate=gate_bind)
    checks = {"flow": flow, "V5b": v5b, "V5c": v5c, "wall_resolution": wr, "deliverables": dl}
    dec = decide(flow, v5b, v5c, wr, dl, iteration=iteration, max_iterations=max_iterations,
                 rule_fn=lambda: recalibration_rule(run_dir, v5b, v5c, flow, delta_next=delta_next, prev_run=prev_run,
                                                    prev_delta_next=prev_delta_next, recalibrations_done=recalibrations_done,
                                                    max_recalibrations=max_recalibrations, min_reduction=recal_min_reduction,
                                                    x_far_m=x_far_m))
    result, nxt, missing, corrections, rule, iter_rec = (dec["result"], dec["next"], dec["missing"], dec["corrections"], dec["rule"],
                                                        dec["iteration"])
    apr_rec = None
    if apr is not None:
        from ..intake import approval as AP
        rc = AP.check_run(apr, run_dir, problem)
        rwg = (apr["protected"].get("requirements_without_gate") or [])
        apr_rec = {"path": str(Path(approval).resolve()), "sha256": AP.sha256_file(approval), "protected_sha256": apr["protected_sha256"],
                   "gate": {k: v for k, v in (apr.get("gate") or {}).items() if k != "env"}, "gate_run": gate_run,
                   "user_ref": apr.get("user_ref"), "run_check": rc,
                   "requirements_without_gate": [r.get("name") for r in rwg]}
        if rc["status"] != "OK":
            # その run は承認した問題の評価になっていない (plan §4.2 要求の一致) → 結果によらず incomplete
            missing = list(missing) + [f"run の設定 (停止規則の基準・物理) が承認と違う ({'; '.join(rc['problems'] + rc['diffs'])[:400]})"]
            result, nxt = "incomplete", "escalate"
        apr_rec["rechecks"] = rechecks
        if rechecks and rechecks["problems"]:
            missing = list(missing) + [f"固定した版での回し直し: {p}" for p in rechecks["problems"]]
            result, nxt = "incomplete", "escalate"
        if rwg:
            missing = list(missing) + [f"判定器の無い要求が未評価 {[r.get('name') for r in rwg]} (plan §4.2: 評価していない要求があれば "
                                       "accepted にしない)"]
            if result == "accepted":
                result, nxt = "incomplete", "escalate"
    if rule is not None and rule.get("calibration") is not None:
        dst = (rule["calibration"]["inputs"]["delta_E"].get("state") or {})
        if dst.get("delta_E_step") is not None and dst["delta_E_step"] != fld["step"]:
            warns.append(f"δ_E の表の抽出の step {dst['delta_E_step']} が判定した場の step {fld['step']} と違う")
    show = {"V5b_exit_core_M": v5b.get("M"), "V5b_range": v5b.get("range"),
            "V5c_maxabs_pct": v5c.get("maxabs_pct"), "V5c_halfamp_pct": v5c.get("halfamp_pct"),
            "V5c_tol": [v5c.get("tol_maxabs_pct"), v5c.get("tol_halfamp_pct")], "V5c_eta": v5c.get("eta"),
            "V5c_compare": v5c.get("compare"), "wall_resolution": wr.get("verdict_line"), "flow": flow.get("status")}
    if rule is not None:
        c = rule["conditions"]
        show.update(delta_error_current_mm=c["c_further_reducible"].get("current_max_abs_mm"),
                    delta_error_projected_mm=c["c_further_reducible"].get("projected_max_abs_mm"),
                    delta_error_reduction=c["c_further_reducible"].get("reduction"),
                    delta_error_prev_mm=c["b_prev_update_reduced_error"].get("prev_max_abs_mm"))
        if rule["calibration"] is not None:
            show["recalibrated_coefficients"] = rule["calibration"]["coefficients"]
    return {"format": FORMAT, "plan": PLAN, "created": _now(), "git_head": CC._git_head(), "run": str(run_dir),
            "inputs": {"card": None if card is None else str(Path(card).resolve()), "card_sha256": None if card is None else CC._sha256(card),
                       "problem": None if problem is None else str(Path(problem).resolve()),
                       "problem_sha256": None if problem is None else CC._sha256(problem),
                       "prev_run": None if prev_run is None else str(Path(prev_run).resolve()), "Md": Md,
                       "Md_source": "カードの M_exit" if Md_card is not None else "prepare_info の Md"},
            "field": fld, "tolerances": {k: v for k, v in tol.items() if k != "settle_line"},
            "result": result, "next": nxt, "missing": missing, "corrections": corrections, "checks": checks,
            "recalibration_rule": None if rule is None else {k: v for k, v in rule.items() if k != "calibration"},
            "recalibration": None if rule is None or rule["calibration"] is None else
            {"coefficients": rule["calibration"]["coefficients"], "overrides": rule["calibration"]["overrides"],
             "overrides_yaml": rule["calibration"]["overrides_yaml"], "errors": rule["calibration"]["errors"],
             "delta_E": rule["calibration"]["inputs"]["delta_E"]},
            "iteration": iter_rec, "show": show, "warnings": warns, "approval": apr_rec}


def approval_mode(approval, card, problem, cli: dict, optional_parts, iteration, max_iterations, max_recalibrations):
    """正式な入口の検査 (plan tooling-design-orchestrator §4.7.1 (1)・(2))。戻り: (承認の記録, 反復の上限, 較正し直しの上限, 無くてもよい部品)。
    CLI の公差・窓・--optional-parts・反復の上限・較正し直しの上限は値が同じでも拒否 (入口を 1 つにする)。カード・問題 YAML が承認の値と
    違えば拒否。--iteration (台帳の反復の番号) は必須。"""
    from ..intake import approval as AP
    if card is None or problem is None:
        raise ValueError("--approval には --card と --problem が要る (承認の値と照合する)")
    given = sorted(k for k, v in (cli or {}).items() if v is not None)
    if given:
        raise ValueError(f"--approval のときは CLI の公差・窓の指定を受け付けない ({given}; 値が同じでも。承認の記録の値を使う)")
    if optional_parts:
        raise ValueError("--approval のときは --optional-parts を受け付けない (承認の記録の qualification.optional_parts を使う)")
    if max_iterations is not None or max_recalibrations is not None:
        raise ValueError("--approval のときは --max-iterations・--max-recalibrations を受け付けない (承認の記録の run_policy を使う)")
    if iteration is None:
        raise ValueError("--approval のときは --iteration (台帳の反復の番号) が要る")
    ck = AP.check(approval, card, problem)
    if ck["status"] != "OK":
        raise ValueError(f"カード・問題 YAML が承認の記録と違う ({ck['kind']}): " + "; ".join(ck["problems"] + ck["diffs"])[:800])
    apr = AP.load_approval(approval)
    rp = apr["protected"]["run_policy"]
    q = apr["protected"]["qualification"]
    return apr, int(rp["max_iterations"]), int(rp["max_recalibrations"]), tuple(q.get("optional_parts") or ())


def approved_tolerances(run_dir, info: dict, apr: dict, approval_path) -> dict:
    """公差・窓を承認の記録の値で決める (出所は承認の記録)。x の窓が規則 (プリセット) なら prepare_info から作る。"""
    base = resolve_tolerances(run_dir, info, None, None, {})
    q = apr["protected"]["qualification"]
    src = f"承認の記録 ({Path(approval_path).name})"
    out = dict(base)
    for k in ("exit_core_M_tol_rel", "line_eta", "line_tol_maxabs_pct", "line_tol_halfamp_pct"):
        out[k] = {"value": float(q[k]), "source": src}
    out["line_n"] = {"value": int(q["line_n"]), "source": src}
    xr = q.get("line_x_range_m")
    if isinstance(xr, list):
        out["line_x_range_m"] = {"value": [float(xr[0]), float(xr[1])], "source": src}
    elif isinstance(xr, dict) and xr.get("rule") == "preset":
        dom = (info.get("physical_wall") or {}).get("domain")
        if info.get("x_E") is None or not (isinstance(dom, (list, tuple)) and len(dom) == 2):
            out["line_x_range_m"] = {"value": None, "source": f"{src} の規則を作れない (prepare_info に x_E か physical_wall.domain が無い)"}
        else:
            S = float(info["scale_m"])
            o0, o1 = xr["x_offsets_rt"]
            out["line_x_range_m"] = {"value": [(float(info["x_E"]) + o0) * S, (float(dom[1]) + o1) * S],
                                     "source": f"{src} の規則 (x_E + {o0:g}, x_F {o1:+g}) r_t × scale_m"}
    else:
        raise ValueError(f"承認の記録の qualification.line_x_range_m が読めない ({xr!r})")
    return out


def summary_lines(rec: dict) -> list[str]:
    c = rec["checks"]
    v5b, v5c, wr, fl, dl = c["V5b"], c["V5c"], c["wall_resolution"], c["flow"], c["deliverables"]
    L = [f"run: {rec['run']}  (場 {rec['field']['res']}: {rec['field']['source']})",
         f"判定: {rec['result']}  次の手: {rec['next']}"]
    if v5b.get("M") is not None:
        L.append(f"V5b 出口コア M {v5b['M']:.6f} (許容 {v5b['range'][0]:.5f}〜{v5b['range'][1]:.5f}; {rec['tolerances']['exit_core_M_tol_rel']['source']}) "
                 f"→ {v5b['verdict']}")
    if "maxabs_pct" in v5c:
        L.append(f"V5c 線 η {v5c['eta']} x {v5c['x_range_m'][0]}〜{v5c['x_range_m'][1]} m: 最大偏差 {v5c['maxabs_pct']:.4f} % (≤ {v5c['tol_maxabs_pct']})・"
                 f"片振幅 {v5c['halfamp_pct']:.4f} % (≤ {v5c['tol_halfamp_pct']}) → {v5c['verdict']}  [窓 {rec['tolerances']['line_x_range_m']['source']}]")
        for k, v in (v5c.get("compare") or {}).items():
            L.append(f"    比較 {k}: 最大偏差 {v['maxabs_pct']:.4f} %・片振幅 {v['halfamp_pct']:.4f} % (判定に使わない)")
        stq = v5c.get("steadiness") or {}
        if stq:
            ser = stq.get("series") or {}
            L.append(f"    線の定常: {stq['basis']}" + (f" → {ser.get('status')} (系列 {ser.get('csv') or ('一時ファイル' if ser.get('csv_sha256') else '作っていない')}; "
                                                       f"{ser.get('verdicts') or ser.get('why')})" if ser else ""))
        if v5c["verdict"] == "INDETERMINATE":
            L.append(f"    判定不能の理由: {v5c.get('why')}" + (f" (線だけの合否 {v5c['line_verdict']})" if v5c.get("line_verdict") else ""))
    else:
        L.append(f"V5c → {v5c['verdict']} ({v5c.get('why')})")
    L.append(f"壁解像 → {wr['verdict']} ({wr.get('verdict_line')})" + (f" — 判定不能: {wr.get('why')}" if wr.get("tool_verdict") else ""))
    L.append(f"計算の成立 → {fl['status']}" + (f" ({fl.get('why')})" if fl.get("why") else ""))
    L.append(f"寸法の関門 → {dl.get('dimensions', {}).get('status')}・成果物 → {dl['status']}"
             + (f"・生成元 → {dl['provenance']['status']}" if dl.get("provenance") else "") + (f" ({dl.get('why')})" if dl.get("why") else ""))
    for m in rec["missing"]:
        L.append(f"欠け: {m}")
    rr = rec.get("recalibration_rule")
    if rr is not None:
        L.append(f"較正し直しの規則: {'当たる → recalibrate' if rr['apply'] else '当たらない → escalate'}")
        for k, cnd in rr["conditions"].items():
            L.append(f"    {k}: {'ok' if cnd['ok'] else 'NO'}  " + json.dumps({kk: vv for kk, vv in cnd.items() if kk not in ("ok", "prev_delta_E", "prev_run")},
                                                                         ensure_ascii=False, default=str))
    if rec.get("recalibration") and rec["next"] == "recalibrate":
        L += ["新しい係数 (カードの overrides の断片):", rec["recalibration"]["overrides_yaml"].rstrip()]
    L += [f"警告: {w}" for w in rec["warnings"]]
    return L


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m forge_design.feedback.wt_qualify",
        description="風洞ノズルの設計の判定 (V5b 出口コア M・V5c 近軸の線の一様性・壁解像・寸法の関門・成果物) と、補正へ戻すときの "
                    "較正し直しの規則 (plan tooling-nozzle-card-local-mesh §4.3・§5.1 #6b)。",
        epilog="判定: accepted / back_to_correction / rejected / incomplete、次の手: done / recalibrate / escalate / make_deliverables。"
               "終了コード 0 accepted、2 recalibrate、3 escalate、4 incomplete、1 入力の不備。出力 <out-dir>/WT_QUALIFICATION.json。")
    ap.add_argument("run", help="NS の run のディレクトリ")
    ap.add_argument("--card", default=None, help="要求カード (M_exit・寸法の関門のカードの照合・qualification の公差)")
    ap.add_argument("--problem", default=None, help="問題 YAML (qualification・evaluate.settle の線)")
    ap.add_argument("--prev-run", default=None, help="壁を更新する前の NS run (規則 (b) 前の更新で δ の誤差が縮んだか)")
    ap.add_argument("--delta-next", default=None, help="この run の δ_E の表を明示 (delta_r_next.csv の形式、.gz 可)")
    ap.add_argument("--prev-delta-next", default=None, help="前の run の δ_E の表を明示")
    ap.add_argument("--step", type=int, default=None, help="判定する全場の step (既定: 停止規則の stop_step → 最後の全場)")
    ap.add_argument("--deliverables", default=None, help="成果物の置き場 (既定 <run>/deliverables)")
    ap.add_argument("--optional-parts", default="", help="無くてもよい成果物の部品 (カンマ区切り; 既定は全部必要)")
    ap.add_argument("--iteration", type=int, default=None, help="この run の反復の番号 (1 = 最初の NS)")
    ap.add_argument("--max-iterations", type=int, default=None, help="反復の上限 (iteration ≥ これで rejected)")
    ap.add_argument("--recalibrations-done", type=int, default=0, help="済んだ較正し直しの回数 (既定 0)")
    ap.add_argument("--max-recalibrations", type=int, default=None, help=f"較正し直しの上限 (既定 {MAX_RECALIBRATIONS}、plan §4.3)")
    ap.add_argument("--approval", default=None,
                    help="承認の記録 (card_approval.json; 正式な入口)。公差・窓・必須成果物・反復と較正し直しの上限は承認の値だけを使い、"
                         "CLI の指定は値が同じでも拒否する。--card・--problem・--iteration が要る (plan tooling-design-orchestrator §4.7.1)")
    ap.add_argument("--recal-min-reduction", type=float, default=RECAL_MIN_REDUCTION,
                    help=f"規則 (c): 較正し直しの射影で最大誤差が減る割合の下限 (既定 {RECAL_MIN_REDUCTION}、plan に数値が無い暫定の値)")
    ap.add_argument("--exit-tol-rel", type=float, default=None, help="V5b の許容 |M/M_d − 1| (既定 1e-3)")
    ap.add_argument("--line-eta", type=float, default=None, help="V5c の線の r/r_w (既定 0.03)")
    ap.add_argument("--line-x-range-m", type=float, nargs=2, default=None, metavar=("X0", "X1"), help="V5c の x 窓 [m] (設計スロートから)")
    ap.add_argument("--line-tol-maxabs", type=float, default=None, help="V5c の最大偏差の許容 [%%] (既定 0.1)")
    ap.add_argument("--line-tol-halfamp", type=float, default=None, help="V5c の片振幅の許容 [%%] (既定 0.1)")
    ap.add_argument("--out-dir", default=None, help=f"出力の置き場 (既定 <run>/{OUT_DIRNAME})")
    ap.add_argument("--no-write", action="store_true", help="ファイルを書かない (壁解像の出力も残さない)")
    ap.add_argument("--json", action="store_true", help="記録の JSON を標準出力に出す")
    a = ap.parse_args(argv)
    cli = {"exit_core_M_tol_rel": a.exit_tol_rel, "line_eta": a.line_eta,
           "line_x_range_m": None if a.line_x_range_m is None else [float(a.line_x_range_m[0]), float(a.line_x_range_m[1])],
           "line_tol_maxabs_pct": a.line_tol_maxabs, "line_tol_halfamp_pct": a.line_tol_halfamp}
    od = None if a.no_write else Path(a.out_dir if a.out_dir else Path(a.run) / OUT_DIRNAME)
    if od is not None:
        od.mkdir(parents=True, exist_ok=True)
    try:
        rec = qualify(a.run, a.card, a.problem, prev_run=a.prev_run, delta_next=a.delta_next, prev_delta_next=a.prev_delta_next,
                      step=a.step, deliverables=a.deliverables, iteration=a.iteration, max_iterations=a.max_iterations,
                      recalibrations_done=a.recalibrations_done, max_recalibrations=a.max_recalibrations,
                      recal_min_reduction=a.recal_min_reduction,
                      optional_parts=tuple(s for s in a.optional_parts.split(",") if s), cli_tols=cli,
                      wall_txt=None if od is None else od / "check_wall_resolution.txt", series_dir=od, approval=a.approval)
    except (ValueError, OSError, KeyError) as e:
        if od is not None and od.is_dir() and not any(od.iterdir()):
            od.rmdir()
        print(f"wt_qualify: 止めた — {e}", file=sys.stderr)
        return 1
    if a.json:
        print(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    else:
        print("\n".join(summary_lines(rec)))
    if od is not None:
        (od / OUT_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
        if rec["next"] == "recalibrate" and rec.get("recalibration"):
            (od / CC.OUT_YAML).write_text(rec["recalibration"]["overrides_yaml"])
        print(f"書いた: {od / OUT_JSON}", file=sys.stderr if a.json else sys.stdout)
    return EXIT_CODES["accepted"] if rec["result"] == "accepted" else EXIT_CODES["incomplete"] if rec["result"] == "incomplete" \
        else EXIT_CODES[rec["next"]]


if __name__ == "__main__":
    sys.exit(main())
