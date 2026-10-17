"""計算の後の判定の確定と判定の仕様 `JUDGE_SPEC.json` の生成 (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12 後半;
ユーザ 2026-10-15「計算後の収束判定は完全に機械化してほしい」「あほなAIでもね」)。

**確定の時点** (§4.6-2): 設計チェーンの runner (`runner_axismach` の Euler `run_staged`・NS `run_staged_ns`、`runner_wt.run_staged`・
`runner_walldriven.run_staged`・`runner_sern.run_staged` (2D・3D・driver_sern)・ベルの `run_staged_tp`・`continue_tp`・`continue_cpg` の後の
起動 (opt/driver・opt/polish・runner.main)・`feedback/euler_loop`) は、**全部の段** (S1・stages=none の本段を含む) を `run_forge`
(`FORGE_JUDGE_DEFER=1` を渡す) で起動する。`run_case.sh` の `judge_run.py --from-run-case` は段の結果と健全性だけを記録し
(`state: deferred`)、親が段の履歴・manifest・停止の記録の確定・出力の整理 (`stop_rule.postprocess()` を含む; 判定が読む系列の
CSV もここで書く) を終えた後に `finalize` で総合の判定を**一度だけ**確定する (`judge_run.py <run>`、段の記録は参考として残る)。
段階起動の親は `Staged` (with 文) で包む: 入るときに `JUDGE_PENDING.json` (pid・host) を置き (Stop フックは親が生きている間の
`state: deferred` を確定待ちとして見逃す)、出るときに forge を起動していれば (例外で抜けるときも) `finalize` する。

**判定の仕様** (§4.6-1): 生成器は実行の前に機種の既定の仕様を run の `JUDGE_SPEC.json` に書く (`write_spec`; `judge_run.normalize_spec`
で検査してから書く)。量の定義・単位・抽出器・領域・窓・閾値・必須かを書き、設計の受理の条件 (停止規則の到達・受け渡しの関門・
SERN のゲート・ベルの tp_gates) は各ゲートに残す (`judge_run.py` は共通の検査だけ)。機種の既定:

| 機種 (spec の source) | 量 (系列の出所) | 閾値 |
| --- | --- | --- |
| axis-Mach NS (`spec_axismach_ns`) | 出口コア M・質量の収支 r・θ_r (停止規則を使う run だけ) (`judge_series.csv`: 停止規則の記録の行、無ければ本段の res から同じ式) | 停止規則を使う run: M・r・θ_r は criterion stop_rule (停止規則の `judge_rows` そのもの、窓・閾値は run の `stop_rule.json` の settings; §5.1 #14a (2))。使わない run: M は check_quasisteady の既定、r は停止規則の既定の絶対の閾値 (|平均|・振幅) |
| axis-Mach NS、停止規則 v2 (`settle_criteria.yaml` のある run; `spec_axismach_settle`) | 基準の量の全部 (残差・流量の収支・出口コア M・線の M・流れ角・壁の熱流など、基準のとおり) | criterion settle: convergence (`check_settled`) と同じ 1 回の判定の各量の verdict (閾値・窓・区間は基準のもの; plan tooling-design-orchestrator §4.2「解析の成立の読み方」、judge plan §4.6 の 8)。仕様は名前・種類と基準の sha256 を持ち、判定のときに照合する |
| axis-Mach Euler (`spec_axismach_euler`) | 出口コア M・r (本段の res から) | 同上。本段の出力が 4 枚に届かなければ準定常は INDETERMINATE (窓の不足) |
| SERN 2D・3D (`spec_sern`) | 目的量 (粘性: C_T_with_shear、非粘性: C_T)・C_T・C_L・C_M (`force_history.csv`; sern_gates の定義) | sern_forces の steadiness (tail 0.4・drift 0.02・osc 0.05) |
| ベル TP (`spec_bell`) | 推力系 4 列・ṁ (`thrust_history.csv`; tp_gates (ii) の列) | tp_gates (ii) の drift・osc (推力系 1.5e-4、ṁ 1e-3) |
| ベル CPG・風洞 (WT・walldriven・euler_loop) (`spec_snapshot`) | shock・machmax・pmax (全 res から; health.gate・opt/driver の check_quasisteady の既定の量、軸対称なので asym は除く) | check_quasisteady の既定 |

収束は auto (manifest があれば `--segment`)。同じ数値設定の継続 (opt/polish の CPG のチャンク・euler_loop の継続チャンク) は
`--from-floor <継続元>` を仕様に書く (`spec_*(..., from_floor=)`)。壁解像は実効の粘性と境界条件で決まる (仕様に書かない)。
"""
from __future__ import annotations

import csv
import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

DEFER = "FORGE_JUDGE_DEFER"
PENDING = "JUDGE_PENDING.json"
SPEC = "JUDGE_SPEC.json"
SERIES = "judge_series.csv"            # axis-Mach の系列 (出口コア M・質量の収支・θ_r)
SERIES_META = "judge_series.json"      # 系列の出所の記録
PLAN = "plans/active/tooling-convergence-and-wall-resolution-gates.md §4.6・§5.1 #12"
LAUNCH_MARGIN_NS = 50_000_000          # forge_run.log の mtime と入った時刻の比較の余裕 (ファイルの時刻は粗いクロック; 50 ms)


def _tools() -> Path:
    from .runner import FORGE_TOOLS
    return FORGE_TOOLS


def _judge_mod():
    t = str(_tools())
    if t not in sys.path:
        sys.path.insert(0, t)
    import judge_run
    return judge_run


def defer_env(env: dict | None = None) -> dict:
    """`run_case.sh` に渡す環境: `env` (既定 runner._ENV) に `FORGE_JUDGE_DEFER=1` を足したもの (段の結果だけを記録させる)。"""
    from .runner import _ENV
    return {**(_ENV if env is None else env), DEFER: "1"}


# ------------------------------------------------------------------------------------------------ 確定

def finalize(run_dir, *, solver_rc: int | None = None, series=None, generator: str | None = None) -> dict:
    """段の記録・manifest・停止の記録の確定・出力の整理が終わった後に**一度だけ**呼び、`judge_run.py <run>` で総合の判定を確定する。
    series: 判定の仕様が読む系列の CSV を書く呼び出し (出力の整理の最後; 失敗は記録して判定は続ける = その量は INDETERMINATE)。
    solver_rc: forge の終了コードが分かっていれば渡す (無ければ judge_run が前の段の記録・forge_run.log から決める)。
    判定の FAIL・判定不能・判定器の例外でも例外にしない。`JUDGE_PENDING.json` を消す。戻り {rc, overall, state, series_error}。"""
    run_dir = Path(run_dir)
    out = {"run": str(run_dir), "generator": generator, "series_error": None}
    if series is not None:
        try:
            series()
        except Exception as e:  # noqa: BLE001 — 系列が書けなければ判定がその量を INDETERMINATE にする
            out["series_error"] = f"{type(e).__name__}: {e}"
            print(f"[judge] {run_dir.name}: 判定の系列を書けない ({out['series_error']}) — その量は判定不能になる", flush=True)
    from .runner import _ENV
    env = {k: v for k, v in _ENV.items() if k != DEFER}
    cmd = ["python3", str(_tools() / "judge_run.py"), str(run_dir.resolve())]
    if solver_rc is not None:
        cmd += ["--solver-rc", str(int(solver_rc))]
    try:
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        out["rc"] = r.returncode
        out["output"] = (r.stdout + r.stderr)[-4000:]
    except OSError as e:
        out["rc"] = 3
        out["output"] = f"judge_run.py を起動できない: {e}"
    try:
        v = json.loads((run_dir / "RUN_VERDICT.json").read_text())
        out["overall"], out["state"] = v.get("overall"), v.get("state")
    except (OSError, ValueError):
        out["overall"], out["state"] = None, None
    (run_dir / PENDING).unlink(missing_ok=True)
    print(f"[judge] {run_dir.name}: OVERALL {out['overall']} (state {out['state']}, judge_run exit {out['rc']}; "
          f"{run_dir / 'RUN_VERDICT.txt'})", flush=True)
    return out


class Staged:
    """段階起動の親の判定の確定 (with 文)。入るときに `JUDGE_PENDING.json` を置き、出るときに、この呼び出しの中で forge を起動していれば
    (`forge_run.log` が入った後に書かれた; 関門で止まった起動は数えない) 例外で抜けるときも `finalize` する。例外は握りつぶさない。
    `solver_rc`・`series` は抜ける前に設定してよい (例: 本段の停止規則の forge の終了コード)。"""

    def __init__(self, run_dir, generator: str, series=None):
        self.run_dir = Path(run_dir)
        self.generator = generator
        self.series = series
        self.solver_rc = None
        self.result = None

    def __enter__(self):
        self.t0 = time.time_ns()
        tmp = self.run_dir / (PENDING + ".tmp")
        tmp.write_text(json.dumps({"pid": os.getpid(), "host": socket.gethostname(), "generator": self.generator,
                                   "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "plan": PLAN}, ensure_ascii=False))
        os.replace(tmp, self.run_dir / PENDING)
        return self

    def launched(self) -> bool:
        try:
            return (self.run_dir / "forge_run.log").stat().st_mtime_ns >= self.t0 - LAUNCH_MARGIN_NS
        except OSError:
            return False

    def __exit__(self, et, ev, tb):
        try:
            if self.launched():
                self.result = finalize(self.run_dir, solver_rc=self.solver_rc, series=self.series, generator=self.generator)
            else:
                (self.run_dir / PENDING).unlink(missing_ok=True)
        except Exception as e:  # noqa: BLE001 — 確定の失敗で元の例外・戻り値を壊さない
            print(f"[judge] {self.run_dir.name}: 判定の確定に失敗 ({type(e).__name__}: {e})", flush=True)
        return False


# ------------------------------------------------------------------------------------------------ 仕様

def write_spec(run_dir, spec: dict) -> Path:
    """仕様を `judge_run.normalize_spec` で検査してから run の `JUDGE_SPEC.json` に書く (一時ファイル → rename)。誤りは ValueError。"""
    JR = _judge_mod()
    try:
        JR.normalize_spec(spec, True)
    except JR.SpecError as e:
        raise ValueError(f"生成器の判定の仕様が誤り: {e}") from e
    p = Path(run_dir) / SPEC
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n")
    os.replace(tmp, p)
    return p


def _qty(name, definition, unit, region, extractor, *, criterion="relative", thresholds=None, tail=0.4, min_points=4,
         accept=("STEADY",), required=True, role=None) -> dict:
    """仕様の量 1 つ。criterion stop_rule (停止規則の judge_rows で判定、role = stop_rule_role) は窓・閾値を書かない (run の
    stop_rule.json の settings から読む; judge_run.normalize_spec が書いたものを誤りにする)。"""
    if criterion == "stop_rule":
        return {"name": name, "definition": definition, "unit": unit, "region": region, "extractor": extractor,
                "criterion": criterion, "stop_rule_role": role, "accept": ["STEADY"], "required": required}
    q = {"name": name, "definition": definition, "unit": unit, "region": region, "extractor": extractor,
         "window": {"tail_frac": tail, "min_points": min_points}, "criterion": criterion,
         "thresholds": dict(thresholds or {"drift": 0.05, "osc": 0.10}), "accept": list(accept), "required": required}
    return q


def _conv(from_floor: str | None) -> dict:
    if from_floor:
        return {"mode": "from_floor", "reference": str(from_floor)}
    return {"mode": "auto"}


def spec_axismach(generator: str, *, stop_rule: dict | None = None, euler: bool = False, from_floor=None) -> dict:
    """axis-Mach の NS (euler=False) と Euler (True) の仕様。stop_rule: 本段で停止規則を使う run の prepare_info の stop_rule (閾値・断面)、
    使わない run は None (θ_r は帯の位置 y_b が決まらないので量に入れない)。
    停止規則を使う NS の run は、出口コア M・質量の収支 r・θ_r を criterion stop_rule (停止規則の judge_rows そのもの、窓・閾値は run の
    stop_rule.json の settings) で判定する (機種の受理のゲートと PASS の意味を揃える; plan tooling-convergence-and-wall-resolution-gates
    §5.1 #12b (2)・#14a (2))。使わない NS と Euler (ゲートが無い) は今までどおり M が check_quasisteady の既定、r が停止規則の既定の絶対の閾値。"""
    from . import stop_rule as SR
    th = SR.thresholds(**{k: stop_rule[k] for k in ("window", "min_points", "r_mean_max", "r_amp_max", "theta_drift_max",
                                                     "theta_amp_max", "x_fracs") if k in (stop_rule or {})}) if stop_rule else SR.thresholds()
    gate = bool(stop_rule) and not euler        # 機種のゲート (本段の停止規則) がある run
    csvx = lambda col: {"kind": "csv_column", "file": SERIES, "column": col, "step_column": "step"}   # noqa: E731
    m_kw = {"criterion": "stop_rule", "role": "relative"} if gate else {}
    r_kw = ({"criterion": "stop_rule", "role": "r"} if gate else
            {"criterion": "absolute", "thresholds": {"abs_mean": float(th["r_mean_max"]), "abs_amp": float(th["r_amp_max"])}})
    qs = [_qty("exit_core_M", "出口コア M: 最後の列 (i = ni − 1) の η = r/r_w ∈ [0.05, 0.7] の節点の M = |U|/c の平均 "
               "(stop_rule.recorded_quantities; nozzle_report.metrics・case/70 qs_series.py と同じ定義)", "-", "出口の列 η ∈ [0.05, 0.7]",
               csvx("exit_core_M"), **m_kw),
          _qty("mass_imbalance_r", "質量の収支 r = (ṁ_in − ṁ_out)/ṁ_in: 入口の列 (i = 0)・出口の列 (i = ni − 1) の折れ線を通る流束 "
               "2π∫(ρu_x dr − ρu_r dx) r (stop_rule.mass_flows; M6 の cold_xcheck.reduce_fields と同じ式)", "-", "入口の列と出口の列",
               csvx("r"), **r_kw)]
    notes = [f"系列 {SERIES} は確定の時点に親が書く (judge.write_axismach_series): 停止規則の記録 stop_rule.json の行 (STOP_RULE_STATE.json の "
             "record_sha256 と一致するときだけ = 定義・出典・区間が停止規則と一致するとき)、無ければ本段の res_<N>.h5 (N > 0) から同じ式で"]
    if gate:
        notes.append("閾値: exit_core_M・r・θ_r は criterion stop_rule = 停止規則の judge_rows そのもの (窓・閾値・断面は run の stop_rule.json "
                     f"の settings; 準備の時点の値は窓 {th['window']} step・{th['min_points']} 点、r は |窓の平均| ≤ {th['r_mean_max']:g}・"
                     f"振幅 ≤ {th['r_amp_max']:g}、θ_r と出口コア M は D ≤ {th['theta_drift_max']:g}・振幅/|平均| ≤ {th['theta_amp_max']:g})。"
                     "系列の最後の出力の時点で判定し、停止規則の到達 (stop_rule の項目) とは別の判定として並べる (照合は stop_rule_cross)")
    else:
        notes.append("閾値: exit_core_M は check_quasisteady の既定 (drift 0.05・osc 0.10; 共通の検査)、r は停止規則の既定の絶対の閾値 "
                     f"(|平均| ≤ {th['r_mean_max']:g}・振幅 ≤ {th['r_amp_max']:g}、末尾窓 4 割)。設計の受理 (停止規則の到達・受け渡しの関門) は各ゲートに残る")
    if euler:
        notes.append("Euler: 本段の出力が 4 枚 (min_points) に届かなければ準定常は INDETERMINATE (窓の不足; 理由つきの適用外にはしない — "
                     "軸 M などの報告する量があるので)")
    elif stop_rule:
        for f in th["x_fracs"]:
            k = SR.theta_key(f)
            qs.append(_qty(k, f"x = {f:g}·x_w (設計壁の出口 physical_wall.domain[1] に対する比) の一定 x の断面の運動量厚さ θ_r (帯の外縁 y_b は "
                           "本段の最初の出力 × 固定の Euler 参照で固定; stop_rule.theta_r_sections)", "r_t", f"x/x_w = {f:g} の断面", csvx(k),
                           criterion="stop_rule", role="theta"))
    else:
        notes.append("θ_r は量に入れない: 本段で停止規則を使わない (帯の位置 y_b を決める Euler 参照が無い) run")
    return {"judge_spec_version": 1, "source": f"generator:{generator}", "convergence": _conv(from_floor),
            "quasisteady": {"quantities": qs}, "notes": notes}


def spec_axismach_settle(generator: str, run_dir) -> dict:
    """停止規則 v2 の axis-Mach NS (run に `settle_criteria.yaml` がある) の仕様 (judge plan §4.6 の 8; 2026-10-18 codex 諮問 judge-v2)。
    準定常の量は基準の量の全部を criterion settle で書く: 数値は judge_run の convergence と同じ 1 回の `check_settled` の判定から写し、
    名前・種類と基準の sha256 を仕様に持って判定のときに照合する (欠け・余り・基準の書き換えは判定不能)。出口コア M・r の
    `judge_series.csv` の量は書かない (定義・窓・閾値が基準と別の判定を二重に持たない)。基準が読めなければ ValueError (起動しない)。"""
    from . import stop_rule as SR
    S = SR._settle()
    lc = S.load_criteria(run_dir)                       # 誤りは CriteriaError
    if lc["status"] != "ok":
        raise ValueError(f"{run_dir}: {S.CRITERIA_FILE} の評価対象が未設定 ({lc.get('reason')}) — 判定の仕様を書けない")
    crit = lc["criteria"]
    qs = []
    for q in crit["quantities"]:
        qs.append({"name": q["name"], "definition": f"settle_criteria.yaml の量 {q['name']} (種類 {q['kind']}; 判定は check_settled、"
                                                    "閾値・窓は基準のもの)",
                   "unit": "-", "region": "settle_criteria.yaml のとおり",
                   "extractor": {"kind": "settle", "quantity": q["name"], "settle_kind": q["kind"]},
                   "criterion": "settle", "accept": ["SETTLED"], "required": True})
    notes = [f"停止規則 v2 の run: convergence は check_settled (judge_run が run_case.sh と同じ道具の選択で決める; PASS = 設定した量・区間で "
             f"SETTLED、残差の低下は意味しない)。準定常は同じ判定の量 {[q['name'] for q in crit['quantities']]} を criterion settle で写す",
             "停止の時の再評価 (stop_rule_v2/SETTLE_VERDICT_stop_<k>.json) と今の入力・起動・終端・量の数値を照合し、通らなければ判定不能"]
    if crit.get("excluded"):
        notes.append("基準が判定から外した量: " + "; ".join(f"{e.get('kind')} ({e.get('reason')})" for e in crit["excluded"]))
    return {"judge_spec_version": 1, "source": f"generator:{generator}", "convergence": {"mode": "auto"},
            "quasisteady": {"quantities": qs, "settle_criteria": {"file": S.CRITERIA_FILE, "sha256": lc["sha256"]}}, "notes": notes}


def spec_sern(generator: str, viscous: bool, from_floor=None) -> dict:
    """SERN 2D・3D: sern_gates.steadiness_gate の量 (目的量 + C_T・C_L・C_M) を force_history.csv から、sern_forces.steadiness の閾値で。"""
    from ..metrics import sern_forces as SF
    from ..metrics.sern_gates import FORCE_KEYS
    obj = "C_T_with_shear" if viscous else "C_T"
    defs = {"C_T_with_shear": ("推力係数 (壁の圧力 + 摩擦): (F_gross + 摩擦の推力寄与)/F_ideal", "ランプ・カウル内外面"),
            "C_T": ("推力係数 (壁の圧力): F_gross/F_ideal、F_gross = 入口の運動量 + (p_in − p_a)H + 壁の推力", "ランプ・カウル内外面"),
            "C_L": ("揚力係数 (壁の圧力): L/F_ideal", "ランプ・カウル内外面"),
            "C_M": ("頭上げモーメント係数 (壁の圧力): M/(F_ideal H)、基準点 spec.moment_ref", "ランプ・カウル内外面")}
    keys = [obj] + [k for k in FORCE_KEYS if k != obj]
    th = {"drift": SF.STEADY_DRIFT, "osc": SF.STEADY_OSC}
    qs = [_qty(k, defs[k][0] + " (metrics.sern_forces.sern_forces; sern_gates.steadiness_gate の量)", "-", defs[k][1],
               {"kind": "csv_column", "file": "force_history.csv", "column": k, "step_column": "step"},
               thresholds=th, tail=SF.STEADY_TAIL, min_points=SF.STEADY_MIN_SNAPS) for k in keys]
    return {"judge_spec_version": 1, "source": f"generator:{generator}", "convergence": _conv(from_floor),
            "quasisteady": {"quantities": qs},
            "notes": ["系列 force_history.csv は確定の時点に親が書く (runner_sern.write_force_series; collect と同じ関数)",
                      f"閾値は sern_forces.steadiness (tail {SF.STEADY_TAIL}・drift {SF.STEADY_DRIFT}・osc {SF.STEADY_OSC}); "
                      "受理 (rc・床・目的量の窓) は sern_gates に残る"]}


def spec_bell(generator: str, from_floor=None) -> dict:
    """ベル TP: tp_gates (ii) の列 (推力系 4 列・ṁ) を thrust_history.csv から、(ii) の閾値で。"""
    from ..opt.driver import TP_QS_DEFAULT, TP_QS_MDOT
    from .runner import THRUST_QS_COLS, THRUST_QS_MDOT_COLS
    defs = {"thrust_exit_N": ("推力 (出口面の会計)", "N"), "thrust_wall_N": ("推力 (壁の会計)", "N"),
            "eta": ("推力効率 η (出口面の会計) = 推力/F_1D", "-"), "eta_wall": ("推力効率 (壁の会計)", "-"),
            "mdot_exit_kg_s": ("出口の質量流量", "kg/s")}
    qs = []
    for cols, q in ((THRUST_QS_COLS, TP_QS_DEFAULT), (THRUST_QS_MDOT_COLS, {**TP_QS_DEFAULT, **TP_QS_MDOT})):
        for c in cols:
            d, u = defs.get(c, (c, "-"))
            qs.append(_qty(c, d + " (metrics.extract.thrust_cv の閉じた CV の会計; runner.thrust_history_tp)", u, "ノズルの閉じた CV",
                           {"kind": "csv_column", "file": "thrust_history.csv", "column": c, "step_column": "step"},
                           thresholds={"drift": float(q["drift"]), "osc": float(q["osc"])}, tail=float(q["tail"]),
                           min_points=int(q["min_snaps"])))
    return {"judge_spec_version": 1, "source": f"generator:{generator}", "convergence": _conv(from_floor),
            "quasisteady": {"quantities": qs},
            "notes": ["系列 thrust_history.csv は確定の時点に親が書く (runner.write_thrust_history; collect と同じ関数)",
                      "閾値は opt/driver の tp_gates (ii) と同じ。(ii) の末尾窓 8 点・漸近値の最終比と (i)(iii)〜(v) はゲートに残る"]}


def spec_snapshot(generator: str, from_floor=None) -> dict:
    """ベル CPG・風洞 (WT・walldriven・euler_loop): health.gate・opt/driver の CPG が回す check_quasisteady の既定の量 (shock・machmax・pmax;
    軸対称の格子では asym は適用外) を全 res_<N>.h5 から、check_quasisteady の既定の閾値で。"""
    defs = {"shock": ("軸近く (|y| < 1 mm) の P が上流の値 + 0.15 (最大 − 上流) を最初に超える x [mm] (check_quasisteady.q_shock)", "mm",
                      "軸の近く"),
            "machmax": ("全域の M = |U|/c の最大 (check_quasisteady.q_machmax)", "-", "全域"),
            "pmax": ("全域の静圧の最大 (check_quasisteady.q_pmax)", "Pa", "全域")}
    qs = [_qty(k, d, u, r, {"kind": "snapshot", "quantity": k}) for k, (d, u, r) in defs.items()]
    return {"judge_spec_version": 1, "source": f"generator:{generator}", "convergence": _conv(from_floor),
            "quasisteady": {"quantities": qs},
            "notes": ["量は health.gate・opt/driver の CPG の check_quasisteady の既定 (BUILTIN から軸対称で適用外の asym を除いた); "
                      "軸 M の系列の判定 (health.axis_M_series_steady) と設計の受理はゲートに残る"]}


def same_numerics(cfg_a: str, cfg_b: str, bcond_a: str, bcond_b: str) -> bool:
    """2 つの段 (config・bcond) が同じ数値設定 (`stage_manifest.stage_key` の hard キー: 方程式・BC・空間離散化) か。継続の run の
    `--from-floor` を仕様に書くかの判定 (違えば継続元の床は基準にならない; 段の連結と同じ基準)。読めなければ False (auto に倒す)。"""
    t = str(_tools())
    if t not in sys.path:
        sys.path.insert(0, t)
    import stage_manifest as SMm
    try:
        return SMm.stage_key(cfg_a, bcond_a) == SMm.stage_key(cfg_b, bcond_b)
    except Exception:  # noqa: BLE001
        return False


# ------------------------------------------------------------------------------------------------ 系列 (確定の時点に書く)

def _write_rows(path: Path, rows: list, cols: list) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step"] + cols)
        for r in rows:
            w.writerow([int(r["step"])] + [("nan" if r.get("nonfinite_all") else
                                            ("" if r.get(c) is None else repr(float(r[c])))) for c in cols])
    os.replace(tmp, path)


def _snapshot_rows(run_dir: Path) -> list:
    """本段の res_<N>.h5 (N > 0) から出口コア M・流量・質量の収支 (停止規則と同じ式) の行を作る。格子は prepare_info の mesh (ni・nj)・
    scale_m・x_E と nozzle.h5 の座標 (stop_rule.RunGeometry と同じ並び; 設計壁の出口は要らない)。"""
    import h5py
    import numpy as np

    from . import stop_rule as SR
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S, x_E = float(info["scale_m"]), float(info["x_E"])
    with h5py.File(run_dir / "nozzle.h5", "r") as h:
        nc = np.array(h["/MESH/COORD"][:]).reshape(-1, 3)
    if nc.shape[0] != ni * nj:
        raise ValueError(f"節点数 {nc.shape[0]} が ni × nj = {ni * nj} でない (構造格子の node run でない)")
    X = (nc[:, 0].astype(float) / S).reshape(ni, nj)
    R = (nc[:, 1].astype(float) / S).reshape(ni, nj)
    Xf = (nc[:, 0] / S).reshape(ni, nj)
    Rf = (nc[:, 1] / S).reshape(ni, nj)
    rows = []
    for n in SR.res_steps(run_dir):
        if n <= 0:
            continue
        V = SR.read_values(run_dir / f"res_{n}.h5")
        row = {"step": n, "nonfinite_all": SR.nonfinite_count(V)}
        if not row["nonfinite_all"]:
            sh = (ni, nj)
            RO = V["ro"].astype(float).reshape(sh)
            mi, mo = SR.mass_flows(X, R, RO, V["Ux"].astype(float).reshape(sh), V["Uy"].astype(float).reshape(sh), S)
            row.update(mdot_in=mi, mdot_out=mo, r=(mi - mo) / mi if mi != 0 else float("nan"))
            row.update(SR.recorded_quantities(Xf, Rf, {k: V[k].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}, x_E))
        rows.append(row)
    return rows


def write_axismach_series(run_dir) -> Path:
    """axis-Mach の系列 `judge_series.csv` (step・exit_core_M・mdot_rt・r・θ_r_<比>) を書く。停止規則の記録 stop_rule.json の行を、
    本段の終了の状態 STOP_RULE_STATE.json の record_sha256 が stop_rule.json と一致するときだけ使う (定義・出典・区間が停止規則と
    一致するとき)。それ以外は本段の res_<N>.h5 (N > 0) から同じ式で作る (θ_r は帯の位置が無いので出さない)。出所は judge_series.json。
    列は 1 行でも値のある量だけ (無い量の列は書かない = 判定は「列が無い」で INDETERMINATE)。非有限の行は全列 nan (判定は FAIL)。"""
    from . import stop_rule as SR
    run_dir = Path(run_dir)
    meta = {"plan": PLAN}
    rows = None
    rec_p, st = run_dir / SR.RECORD, SR.read_state(run_dir)
    if rec_p.exists() and st is not None and st.get("record_sha256"):
        sha = SR._sha256(rec_p)
        if sha == st["record_sha256"]:
            rec = json.loads(rec_p.read_text())
            rows = [r for r in rec.get("rows") or [] if int(r.get("step", 0)) > 0]
            meta["source"] = {"kind": "stop_rule.json", "sha256": sha, "status": st.get("status"), "n_rows": len(rows)}
        else:
            meta["stop_rule_rejected"] = f"stop_rule.json の sha256 {sha[:16]} が STOP_RULE_STATE.json の record_sha256 と違う"
    elif rec_p.exists():
        meta["stop_rule_rejected"] = "STOP_RULE_STATE.json (本段の終了の状態) が無いか record_sha256 が無い"
    if rows is None:
        rows = _snapshot_rows(run_dir)
        meta["source"] = {"kind": "res_<N>.h5 (N > 0)", "n_rows": len(rows)}
    keys = ["exit_core_M", "mdot_rt", "r"] + sorted({k for r in rows for k in r if k.startswith("theta_r_") and not k.startswith("x_")})
    cols = [k for k in keys if any(k in r for r in rows)]
    meta["columns"] = cols
    path = run_dir / SERIES
    _write_rows(path, rows, cols)
    (run_dir / SERIES_META).write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    return path
