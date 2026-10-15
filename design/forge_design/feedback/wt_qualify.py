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
    exit_core.M_mean)。許容 |M/M_d − 1| ≤ tol (既定 5e−4 = ±0.05 %、plan §6 V5b)。
  - V5c 近軸の線の一様性: 線 η (既定 0.03) の M を `nozzle_report.eta_line` で x 窓に n 点 (既定 2001) 取り、
    最大偏差 max|100(M/M_d − 1)| と片振幅 100(max M − min M)/(2 M_d) (`settle.line_quantities` = `v5c_series.py` の式)。
    許容 0.1 %・0.05 % (plan §6 V5c)。η 0.1 と軸 (η 0) は比較用に併記 (判定に使わない)。補間に軸の節点が重み > 1e−6 で入る
    (第 1 層の η が線より上) か、窓が格子の外・非有限なら判定不能 (plan §6 V5c、codex M3)。
    **この量は「近軸の線 (η) の一様性」であって試験部全体の一様性ではない** (plan §6 V5c)。
  - 壁解像: `check_wall_resolution.py <run> --step <場の step> --groups <bcond の wall* 全部> --weight area` の VERDICT
    (y₁⁺ > 1 の面積 ≤ 2 %、plan §6 V4)。
  - 計算の成立: `settle_criteria.yaml` のある run は停止規則 v2 の REACHED と停止の step の再評価 SETTLED、旧方式は
    `STOP_RULE_STATE.json` の REACHED。状態が無い run (記録の写し等) は UNVERIFIED (accepted にも自動の較正し直しにもしない)。
  - 寸法の関門と成果物: `<run>/deliverables/` の `DIMENSIONS.json` (verdict PASS、カードを渡したらその sha256 と一致) と `MANIFEST.json`
    (verdict PASS・失敗/飛ばした部品なし・必要な部品 `deliverables.PARTS` がそろい・列挙したファイルの sha256 が一致)。無ければ ABSENT。

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
`contur_overrides.yaml` (新しい係数のカードの overrides の断片)。
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
PLAN = ("plans/active/tooling-nozzle-card-local-mesh.md §5.1 #6b・§4.3・§6 V5; "
        "plans/active/tooling-design-orchestrator.md §4.2")
OUT_DIRNAME = "wt_qualify"
OUT_JSON = "WT_QUALIFICATION.json"
REPO = Path(__file__).resolve().parents[3]
TOOLS = REPO / "solver_density_cuda" / "tools"

EXIT_ETA_RANGE = (0.05, 0.7)            # 出口コアの η (stop_rule.recorded_quantities・settle の exit_eta_range と同じ)
COMPARE_ETAS = (0.1, 0.0)               # 比較用に併記する線 (判定に使わない; plan §6 V5c)
DEFAULT_TOLS = {                        # plan §6 V5b・V5c (ユーザ了承 2026-10-15「よい」)
    "exit_core_M_tol_rel": (5.0e-4, "plan §6 V5b (5.5 ± 0.05 %、M6 の生産のゲートと同じ幅)"),
    "line_tol_maxabs_pct": (0.1, "plan §6 V5c (最大絶対偏差 ≤ 0.1 %)"),
    "line_tol_halfamp_pct": (0.05, "plan §6 V5c (片振幅 ≤ 0.05 %)"),
}
QUAL_KEYS = ("exit_core_M_tol_rel", "line_eta", "line_x_range_m", "line_tol_maxabs_pct", "line_tol_halfamp_pct")
DEFAULT_PRESET = "wt_nozzle_v2"         # 雛形 axismach_final_v1 の evaluate.settle の既定
RECAL_MIN_REDUCTION = 0.5               # 規則 (c) の暫定の値 (plan に数値が無い; 上の docstring)
MAX_RECALIBRATIONS = 1                  # plan §4.3「1 回まで」
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
    bad = [k for k in q if k not in QUAL_KEYS]
    if bad:
        raise ValueError(f"{where} の qualification の未知のキー {bad} (書けるのは {list(QUAL_KEYS)})")
    out = {}
    for k, v in q.items():
        if k == "line_x_range_m":
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


def flow_check(run_dir) -> dict:
    """計算の成立: PASS / FAIL / UNVERIFIED (判定は停止規則の記録だけを読む; 収束・関門の判定は judge_run・受け渡しの関門の役目)。"""
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    st_path = run_dir / SR.STATE
    v2 = (run_dir / "settle_criteria.yaml").is_file()
    if not st_path.is_file():
        if v2:
            return {"status": "FAIL", "method": "v2", "why": f"settle_criteria.yaml のある run に {SR.STATE} が無い (本段が終わっていない・記録が失われた)"}
        return {"status": "UNVERIFIED", "method": None, "why": f"{SR.STATE} が無い (記録の写し・停止規則を使わない run) — 計算の成立を確かめられない"}
    s = json.loads(st_path.read_text())
    rec = {"method": s.get("method", "v1" if s.get("version") == 1 else None), "stop_status": s.get("status"), "stop_step": s.get("stop_step")}
    if v2 or s.get("method") == "v2":
        rv = (s.get("reeval") or {}).get("verdict")
        rec["reeval_verdict"] = rv
        ok = s.get("status") == "REACHED" and rv == "SETTLED"
        rec.update(status="PASS" if ok else "FAIL",
                   why=None if ok else f"停止規則 v2 の状態 {s.get('status')}・停止の step の再評価 {rv} (REACHED・SETTLED が要る)")
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


# ------------------------------------------------------------------ 壁解像・成果物
def wall_resolution(run_dir, step: int, out_txt=None) -> dict:
    """`check_wall_resolution.py --step --groups <wall*> --weight area` を回して VERDICT を読む。"""
    run_dir = Path(run_dir)
    bc = yaml.safe_load((run_dir / "bcondConfig.yaml").read_text()) or {}
    walls = [k for k, v in bc.items() if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")]
    if not walls:
        return {"verdict": "INDETERMINATE", "why": "bcondConfig に no-slip 壁 (kind wall*) が無い"}
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
    if out["verdict"] == "INDETERMINATE" and not m:
        out["tail"] = txt[-1500:]
    return out


def deliverables_check(run_dir, deliv_dir=None, card_path=None, optional_parts=()) -> dict:
    """成果物 (MANIFEST.json) と寸法の関門 (DIMENSIONS.json)。status: PASS / FAIL / INDETERMINATE / ABSENT。"""
    from ..export.deliverables import PARTS
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
    miss = [p for p in rec["required_parts"] if p not in parts or (isinstance(parts.get(p), dict) and parts[p].get("ok") is False)]
    if miss:
        probs.append(f"必要な部品が欠けた・失敗した {miss} (skipped {man.get('skipped')}, errors {list((man.get('errors') or {}).keys())})")
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
        warns.append(f"MANIFEST の run {man.get('run')} がこの run ({run_dir.name}) と名前が違う")
    rec.update(status="PASS" if not probs and not dprobs else "FAIL", problems=probs, warnings=warns,
               manifest={"verdict": man.get("verdict"), "created": man.get("created"), "git_commit": man.get("git_commit"),
                         "n_files": len(man.get("files_sha256") or {})})
    return rec


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
    順: 判定不能・計算の成立の FAIL → incomplete (escalate)。公差・壁解像・寸法の FAIL → 反復の上限なら rejected、でなければ
    back_to_correction (壁の形だけなら較正し直しの規則 → recalibrate / escalate、壁解像・寸法を含めば escalate)。すべて合格で成果物が
    無い・不完全、または計算の成立が UNVERIFIED → incomplete (成果物だけが無く成立を確かめていれば make_deliverables)。残りは accepted。
    戻り: {result, next, missing, corrections, rule, iteration}。"""
    corrections, missing = [], []
    for nm, c in (("V5b", v5b), ("V5c", v5c), ("wall_resolution", wr)):
        if c.get("verdict") not in ("PASS", "FAIL"):
            missing.append(f"{nm} が判定不能 ({c.get('why') or c.get('verdict_line') or c.get('verdict')})")
    if flow.get("status") == "FAIL":
        missing.append(f"計算の成立: {flow.get('why')}")
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
            iteration=None, max_iterations=None, recalibrations_done: int = 0, max_recalibrations: int = MAX_RECALIBRATIONS,
            recal_min_reduction: float = RECAL_MIN_REDUCTION, optional_parts=(), cli_tols: dict | None = None,
            x_far_m: float = CC.X_FAR_M, wall_txt=None) -> dict:
    """1 本の NS run を判定して記録の dict を返す (ファイルは書かない; 壁解像の出力だけ wall_txt に書く)。"""
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
    tol = resolve_tolerances(run_dir, info, card_doc, prob_doc, cli_tols or {})
    fld = select_field(run_dir, step)
    flow = flow_check(run_dir)
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
    wr = wall_resolution(run_dir, fld["step"], wall_txt)
    dl = deliverables_check(run_dir, deliverables, card, optional_parts)
    checks = {"flow": flow, "V5b": v5b, "V5c": v5c, "wall_resolution": wr, "deliverables": dl}
    dec = decide(flow, v5b, v5c, wr, dl, iteration=iteration, max_iterations=max_iterations,
                 rule_fn=lambda: recalibration_rule(run_dir, v5b, v5c, flow, delta_next=delta_next, prev_run=prev_run,
                                                    prev_delta_next=prev_delta_next, recalibrations_done=recalibrations_done,
                                                    max_recalibrations=max_recalibrations, min_reduction=recal_min_reduction,
                                                    x_far_m=x_far_m))
    result, nxt, missing, corrections, rule, iter_rec = (dec["result"], dec["next"], dec["missing"], dec["corrections"], dec["rule"],
                                                        dec["iteration"])
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
            "iteration": iter_rec, "show": show, "warnings": warns}


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
    else:
        L.append(f"V5c → {v5c['verdict']} ({v5c.get('why')})")
    L.append(f"壁解像 → {wr['verdict']} ({wr.get('verdict_line')})")
    L.append(f"計算の成立 → {fl['status']}" + (f" ({fl.get('why')})" if fl.get("why") else ""))
    L.append(f"寸法の関門 → {dl.get('dimensions', {}).get('status')}・成果物 → {dl['status']}" + (f" ({dl.get('why')})" if dl.get("why") else ""))
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
    ap.add_argument("--max-recalibrations", type=int, default=MAX_RECALIBRATIONS, help=f"較正し直しの上限 (既定 {MAX_RECALIBRATIONS}、plan §4.3)")
    ap.add_argument("--recal-min-reduction", type=float, default=RECAL_MIN_REDUCTION,
                    help=f"規則 (c): 較正し直しの射影で最大誤差が減る割合の下限 (既定 {RECAL_MIN_REDUCTION}、plan に数値が無い暫定の値)")
    ap.add_argument("--exit-tol-rel", type=float, default=None, help="V5b の許容 |M/M_d − 1| (既定 5e-4)")
    ap.add_argument("--line-eta", type=float, default=None, help="V5c の線の r/r_w (既定 0.03)")
    ap.add_argument("--line-x-range-m", type=float, nargs=2, default=None, metavar=("X0", "X1"), help="V5c の x 窓 [m] (設計スロートから)")
    ap.add_argument("--line-tol-maxabs", type=float, default=None, help="V5c の最大偏差の許容 [%%] (既定 0.1)")
    ap.add_argument("--line-tol-halfamp", type=float, default=None, help="V5c の片振幅の許容 [%%] (既定 0.05)")
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
                      wall_txt=None if od is None else od / "check_wall_resolution.txt")
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
