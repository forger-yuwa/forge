#!/usr/bin/env python3
r"""収束の判定の新方式への呼び出し側の移行 (plan tooling-settle-convergence-and-monitor §4.9 の 1〜3、§5.1 #7; codex diagnose
2026-10-16 settle-migration-7 の M1・M2・M5・Minor と判別の A/B) の試験。GPU・forge・変換器は使わない。

合成の run (h370 の試験と同じ構造格子の場 56 × 48・固定の Euler の参照) に、停止規則 v2 が書く記録 (STOP_RULE_STATE.json・stop_rule.json・
SETTLE_VERDICT.json・STOPPED_AT・`stop_rule_v2/delta_E_<step>/`) を**本物の関数で**作る: 判定は `settle.judge`・`settle.write_verdict`、
状態は `stop_rule._v2_state`・`write_state`、δ_E は `stop_rule.extract_delta_e` (本物の抽出器)。本段の起動の記録 (MAIN_STAGE_STARTED・
forge_launches.jsonl・solverConfig.yaml) は runner_axismach と forge と同じ形。基準は wt_nozzle_v2 の形 (窓 1 万・保持 0・残差は毎 step)。

  (a) 関門の新方式 (`deltastar_loop.handoff_gate`、settle_criteria.yaml のある run):
      - A (正例): REACHED (停止 26,001、非格子) + 再評価 SETTLED → 進む (method v2・settled true)、run の直下の delta_r_next.csv は停止規則の
        抽出の表とバイト同一、**関門の中で抽出し直さない** (extract_and_merge・extract_delta_e・check_convergence を呼ぶと失敗する偽物)。
        `_after_ns` も同じ表から fixed_point.json を作る。
      - 判別の A/B (codex): B = A と同じ入力で MAIN_STAGE_STARTED の main_config.fnv だけを不一致 → 止める (上書きなし、理由に照合、
        delta_r_next.csv を公開しない)。同じ形の負例: 基準の sha256 の不一致・STOPPED_AT の step の不一致・δ_E の status が ok でない・
        起動の記録が無い・再評価の終端が停止の step と違う・抽出の設定の不一致。どれも上書き 3 つを全部渡しても通さない。
      - (a) CENSORED・FAILED は止める; CENSORED は allow_not_reached で通る (settled false)、FAILED は上書きでも通らない。
      - (b)(c) 再評価の全体: 残差の閾値の超過 (CENSORED の run) は allow_unconverged、流量の収支の超過は allow_unsettled で通り、
        目的の量の判定不能 (節点の出力の欠け) は上書きでも通らない。
      - (d)(f)(g)(h) が新方式でも効く: y₁⁺ FAIL・Euler の参照の違い・ω の違い・使用域の held (表と記録の sha を揃えた上で)・次の壁の不合格。
  (b) 移行: 旧方式の run (基準のファイルなし) は今の関門 (method v1) のまま。同じ場で旧方式の関門が作る delta_r_next.csv と新方式の表が
      同じ形式 (列の見出し・行・値がビット同一)。既存の関門の試験 (run_h370_prelaunch_tests・run_handoff_diag_mark_tests) は別に回す。
  (c) run_case.sh の道具の選択 (偽の forge = forge_launches.jsonl に起動を 1 行足すだけ): 本段 + 基準 → check_settled.py (1 行目 TOOL)、
      基準なし・前の段と確かめられた run (MAIN_STAGE_STARTED なし・stage_manifest.json に本段の段なし) → check_convergence.py、
      手で組んだ run (基準が本段の記録を指さない) → check_settled.py、照合の不一致・基準の誤り・**本段でない起動 (最後の起動の cfg_fnv が
      本段の記録と違う)・起動の記録の欠け・記録が無く前の段と確かめられない run** → check_settled.py の INDETERMINATE (旧の道具へ落とさない;
      §5.1 #7c、result 段レビュー M1)、選択に失敗 → check_settled.py、check_settled.py が VERDICT を出さない → INDETERMINATE・終了コード 3。
      停止規則の後処理 `stop_rule.postprocess` も同じ選択。
  (d) Stop フック (`hook_convergence_gate.py` を一時の根に写して): 両方の道具の CONVERGENCE_VERDICT.txt で block しない、無い・古いと
      今どおり block、案内文に選択の規則と「通過は合格ではない」。
  (e) #7b 残りの読み手 (同じ選択の規則): 成果物の `wall_tables.thermal_use` (新方式の run は check_settled.py --bind・停止規則の
      再評価と照合・道具の名前を verdict_line に; 照合の不一致・再評価との食い違いは未評価、NOT_SETTLED は FAIL、前の段と確かめられた
      run は今どおり、確かめられない run は未評価)、
      旧方式の run は変更前 (HEAD の wall_tables を git から読む) と出力がビット同一。`deliverables.part_excel` の status.json に載る。
      設計の報告の収束の欄 (`nozzle_report.convergence_gate_row`) は 1 行目の TOOL で読み分け、旧い形・旧方式は変更前と同じ。
      `hook_forge_guard.py` の案内文。
  (e′) result 段レビュー 2 回目 M3 (§5.1 #11): 成果物の収束と停止規則の再評価の照合 (`stop_rule.match_reeval`) は、同じ run を相対・絶対・
      別の場所のパスで渡しても、記録の浮動小数を全部 1 ulp ずらしても一致 (修正前 = PRE_REVIEW2_COMMIT の wall_tables は「未評価」; 検出力)、
      各検査の合否・量の verdict・許容を超える値・全体の verdict・入力の中身が違えば不一致 (未評価)、記録のファイルの一覧が無い・記録が無い
      旧い形は例外にせず理由つきの未評価。照合の部品の単体 (所在のパス・表示の文は比べない、絶対の床 1e−12 の扱い)。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_settle_migration_tests.py
所要: 約 1 分 (本物の δ_E の抽出 3 回)。出力は一時ディレクトリ (全部合格なら消す)。
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, str(TOOLS))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")

from forge_design.evaluate import settle_presets as SP  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.metrics import deltastar as DS  # noqa: E402
import stage_manifest as SM  # noqa: E402

ST = SR._settle()
TMP = Path(tempfile.mkdtemp(prefix="settle_migration_"))
FAIL = 0
NCHK = 0


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


# ====================================================================================== 合成の構造格子 (run_h370_prelaunch_tests と同じ)
S_SC = 0.0766
NI, NJ = 56, 48
X0, X1 = -2.0, 20.0
DT, W = 2000, 10000
ALL_OV = dict(allow_not_reached="試験: 未到達", allow_unconverged="試験: 未収束", allow_unsettled="試験: 落ち着いていない")


def rwE_f(x):
    return 1.0 + 0.002 * x * x


def dl_f(x):
    return 0.01 + 0.002 * (x + 2.0)


def core_q(x, r, rw):
    return 100.0 * (1.0 + 0.02 * np.tanh(x / 10.0)) * (1.0 - 0.01 * (r / rw) ** 2)


def radial(n):
    return 1.0 - (1.0 - np.linspace(0.0, 1.0, n)) ** 2


def write_field(rd: Path, X, R, Q, steps, scale=S_SC, perturb=None):
    """node 構造格子 (nozzle.h5 の COORD [m]・res_<N>.h5 の ro・Ux・Uy・sonic)。perturb(step) → Q の倍率。"""
    ni, nj = X.shape
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * scale, R.ravel() * scale, np.zeros(ni * nj)].ravel()
    for st in steps:
        k = 1.0 if perturb is None else perturb(st)
        with h5py.File(rd / f"res_{st}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(ni * nj)
            f["/VALUE/Ux"] = (Q * k).ravel()
            f["/VALUE/Uy"] = np.zeros(ni * nj)
            f["/VALUE/sonic"] = np.full(ni * nj, 300.0)


def synth_ns():
    x = np.linspace(X0, X1, NI)
    rw = rwE_f(x) + dl_f(x)
    R = rw[:, None] * radial(NJ)[None, :]
    d = rw[:, None] - R
    slope = (0.004 * x + 0.002) / np.sqrt(1.0 + (0.004 * x + 0.002) ** 2)
    beta = np.clip(1.0 - d / (0.1 * rw[:, None]), 0.0, 1.0)
    X = x[:, None] + beta * d * slope[:, None]
    y = rw[:, None] - R
    Q = core_q(X, R, rw[:, None]) * np.minimum(1.0, np.maximum(y / (5.0 * dl_f(X)), 0.0) ** (1.0 / 7.0))
    return X, R, Q


def make_euler():
    x = np.linspace(X0 - 0.5, X1 + 0.5, 70)
    rw = rwE_f(x)
    R = rw[:, None] * radial(NJ)[None, :]
    X = np.repeat(x[:, None], NJ, axis=1)
    rd = TMP / "euler"
    rd.mkdir()
    write_field(rd, X, R, core_q(X, R, rw[:, None]), steps=(1000, 2000), perturb=lambda st: 1.0 if st == 1000 else 1.003)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": 70, "nj": NJ}, "scale_m": S_SC, "x_E": 15.0}))
    return rd, {"run": str(rd), "res": "res_1000.h5", "sha256": DS._sha256_file(rd / "res_1000.h5")}


def ns_info(pin):
    sr = SR.resolve_settings({"stop_rule": True}, nsteps_arg=None, n_main=200000, n_main_source="試験", out_int=DT, euler_ref=pin)
    return {"mesh": {"ni": NI, "nj": NJ}, "scale_m": S_SC, "x_E": 15.0, "wall_thermal": {"mode": "isothermal", "Tw": 370.0},
            "physical_wall": {"domain": [X0, X1]}, "throat_physical": {"r": 1.012}, "euler_ref": SR.euler_ref_spec(pin), "stop_rule": sr}


# ====================================================================================== 新方式の run (停止規則 v2 の記録を本物の関数で)
RES_COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe"]


def criteria(segment=True, residual="judge"):
    """residual: judge = 残差を判定する (wt_nozzle_v2 の形)、watch = 残差を外し非有限の見張りだけ (wt_nozzle_v3 の形; plan §4.4・§5.1 #8)。"""
    q_res = ({"kind": "residual", "file": "residual_history.csv", "series": "raw", "required": ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe"],
              "inactive_ok": ["rms_roUz"],
              "thresholds": {"slope_dec": math.log10(1.05), "amp_dec": 0.1, "rise_dec": math.log10(1.10), "excursion_dec": 0.5}}
             if residual == "judge" else {"kind": "residual_finite", "file": "residual_history.csv", "thresholds": {}})
    q = [q_res,
         {"kind": "mass_balance", "file": "boundary_massflux.csv", "inlet": [1], "outlet": [2],
          "thresholds": {"mean_abs": 5e-4, "swing": 5e-3, "mdot_change_rel": 1e-4}},
         {"kind": "exit_core_mach", "eta_range": [0.05, 0.7], "thresholds": {"tol": 1e-4}}]
    return {"schema": ST.SCHEMA, "source": {"made_by": "run_settle_migration_tests.py"},
            "structured": {"ni": NI, "nj": NJ, "scale_m": S_SC},
            "grid": {"out_step_interval": DT, "start_step": 0, "segment": SP.main_stage_segment() if segment else None},
            "window": {"W": W, "min_points": W // DT + 1}, "hold": 0, "quantities": q,
            "excluded": ([] if residual == "judge" else [{"kind": "residual", "reason": "試験: 残差を止める判定から外した (wt_nozzle_v3 の形)"}]),
            "preconditions": []}


CFG_TEXT = "time: {last: {nStepOuter: %d}, outStepInterval: 2000}\noutput: {boundaryMassflux: 1, nodeProbeFile: settle_nodes.txt}\n"


def write_launch(rd: Path, cfg_text: str | None = None):
    """forge と同じ起動の記録 (forge_launches.jsonl に cfg_fnv = solverConfig.yaml の FNV-1a 64 の 1 行)。"""
    cfg = (rd / "solverConfig.yaml").read_bytes() if cfg_text is None else cfg_text.encode()
    with open(rd / "forge_launches.jsonl", "a") as f:
        f.write(json.dumps({"time": int(time.time()), "cfg_fnv": SM.fnv1a64(cfg), "exe_size": 0}) + "\n")


def write_main_notice(rd: Path, cap: int):
    """runner_axismach._write_main_stage_started と同じ鍵 (event・main_config の sha256・nStepOuter・fnv・settle_criteria の sha256)。"""
    cfg = (rd / "solverConfig.yaml").read_bytes()
    rec = {"event": "main_stage_started", "time_unix": time.time(), "stages": "none",
           "main_config": {"file": "solverConfig.yaml", "sha256": DS._sha256_file(rd / "solverConfig.yaml"), "nStepOuter": cap,
                           "fnv": SM.fnv1a64(cfg)},
           "start_step": 0, "settle_criteria": {"file": ST.CRITERIA_FILE, "sha256": DS._sha256_file(rd / ST.CRITERIA_FILE)}}
    (rd / ST.MAIN_STAGE_RECORD).write_text(json.dumps(rec, indent=1, ensure_ascii=False))


def write_data(rd: Path, k: int, *, resid_fn=None, mdot_fn=None, skip_nodes=()):
    """forge の出力の形: 残差の outer_end の行 0〜k−1 (毎 step)、流量の行と res_nodes_<s>.h5 を格子点 2,000〜k と停止の step に。
    節点は出口の列 i = ni−1 の全節点 (exit_core_mach の要件)。"""
    X, R, _ = synth_ns()
    resid_fn = resid_fn or (lambda s: 1e-4)
    with open(rd / "residual_history.csv", "w") as f:
        f.write("step,inner,phase," + ",".join(RES_COLS) + "\n")
        for s in range(0, k):
            v = resid_fn(s)
            f.write(f"{s},-1,outer_end,{v!r},{v!r},{v!r},0.0,{v!r}\n")
    obs = sorted(set(range(DT, (k // DT) * DT + 1, DT)) | {k})
    mdot_fn = mdot_fn or (lambda s: (10.0, 10.0))
    with open(rd / "boundary_massflux.csv", "w") as f:
        f.write("# 試験\nstep,mdot_inlet_1,raw_inlet_1,mdot_outlet_2,raw_outlet_2,mdot_wall_3,raw_wall_3\n")
        for s in obs:
            mi, mo = mdot_fn(s)
            f.write(f"{s},{-mi!r},{-mi!r},{mo!r},{mo!r},0.0,0.0\n")
    ids = np.array([(NI - 1) * NJ + j for j in range(NJ)], dtype=np.int64)
    C = np.c_[X.ravel() * S_SC, R.ravel() * S_SC, np.zeros(NI * NJ)][ids]
    for s in obs:
        if s in skip_nodes:
            continue
        with h5py.File(rd / f"res_nodes_{s}.h5", "w") as h:
            h.attrs["step"] = s
            h["MESH/COORD"] = C.ravel()
            h["MESH/NODE_ID"] = ids
            for key, val in (("ro", 0.1), ("Ux", 5.5 * 340.0), ("Uy", 0.0), ("Uz", 0.0), ("P", 1e4), ("T", 300.0), ("sonic", 340.0)):
                h["VALUE/" + key] = np.full(len(ids), val)


def build(name: str, pin, *, kind="REACHED", resid_fn=None, mdot_fn=None, skip_nodes=(), delta_from: Path | None = None,
          residual="judge") -> Path:
    """新方式の run。REACHED: 上限 100,000・要求の格子点 26,000・協調停止 26,001 (STOPPED_AT)。CENSORED: 上限 = 停止の step 26,000
    (STOPPED_AT なし)。delta_from: 同じ場の別の run の δ_E の抽出を写す (場は同じなので抽出も同じ; 試験の時間を縮める)。"""
    rd = TMP / name
    rd.mkdir()
    k, cap = steps_of(kind)
    X, R, Q = synth_ns()
    write_field(rd, X, R, Q, steps=(0, k))
    info = ns_info(pin)
    (rd / "prepare_info.json").write_text(json.dumps(info))
    (rd / "solverConfig.yaml").write_text(CFG_TEXT % cap)
    return finish_v2(rd, kind=kind, resid_fn=resid_fn, mdot_fn=mdot_fn, skip_nodes=skip_nodes,
                     delta=("from", delta_from) if delta_from is not None else ("extract", info), residual=residual)


def steps_of(kind):
    """(停止の step, 上限): REACHED は協調停止 26,001・上限 100,000、CENSORED は上限 = 停止の step 26,000。"""
    return (26001, 100000) if kind == "REACHED" else (26000, 26000)


def finish_v2(rd: Path, *, kind="REACHED", resid_fn=None, mdot_fn=None, skip_nodes=(), delta=None, residual="judge") -> Path:
    """run (場・solverConfig.yaml は置いてある) に停止規則 v2 の記録一式を書く: 基準・起動の記録・本段の開始の記録・forge の出力の形の
    データ・STOPPED_AT・停止の step の再評価 (SETTLE_VERDICT.json)・stop_rule.json・STOP_RULE_STATE.json。
    delta: ("extract", prepare_info) = 本物の抽出、("from", run) = 同じ場の別の run の抽出を写す、None = δ_E の記録なし (関門には使えない)。"""
    k, cap = steps_of(kind)
    (rd / ST.CRITERIA_FILE).write_text(yaml.safe_dump(criteria(residual=residual), allow_unicode=True, sort_keys=False))
    with open(rd / "forge_launches.jsonl", "w") as f:                  # 前の段 (S1) の起動 → 本段の起動
        f.write(json.dumps({"time": 1, "cfg_fnv": "0123456789abcdef"}) + "\n")
    write_launch(rd)
    write_main_notice(rd, cap)
    write_data(rd, k, resid_fn=resid_fn, mdot_fn=mdot_fn, skip_nodes=skip_nodes)
    if kind == "REACHED":
        files = [{"path": f, "bytes": (rd / f).stat().st_size} for f in ("boundary_massflux.csv", "residual_history.csv")
                 + tuple(f_ for f_ in (f"res_{k}.h5", f"res_nodes_{k}.h5") if (rd / f_).exists())]
        (rd / SR.STOPPED_AT).write_text(json.dumps({"status": "ok", "step": k, "residual_step": k - 1, "request": "{}",
                                                    "nStepOuter": cap, "files": files}, indent=2) + "\n")
    # 停止の step の再評価 (watch_main_stage_v2 と同じ: 新しい RunData・end_step = 停止の step → SETTLE_VERDICT.json)
    lc = ST.load_criteria(rd)
    data = ST.RunData(rd, lc["criteria"])
    r = ST.judge(rd, k, data=data, loaded=lc)
    pv = ST.write_verdict(rd, r, data)
    reeval = {**SR.judge_summary(r), "end_source": r.get("end_source"), "file": ST.VERDICT_FILE, "sha256": SR._sha256(pv)}
    od = rd / SR.V2_DIR / f"delta_E_{k}"
    sr_info = None
    if delta is None:
        de = None
    elif delta[0] == "extract":
        sr_info = delta[1]["stop_rule"]
        de = SR.extract_delta_e(rd, k, SR.euler_ref_spec(sr_info["euler_ref"]), od)
    else:
        shutil.copytree(delta[1] / SR.V2_DIR / f"delta_E_{k}", od)
        de = json.loads((delta[1] / SR.RECORD).read_text())["delta_E"]
        de["dir"] = str(od)
    if sr_info is None and (rd / "prepare_info.json").exists():
        sr_info = json.loads((rd / "prepare_info.json").read_text()).get("stop_rule")
    res_k = rd / f"res_{k}.h5"
    rec = {"version": 2, "method": "v2", "plan": SR.V2_PLAN, "status": kind, "started": SR._now(), "prepare_stop_rule": sr_info,
           "cap": cap, "criteria": {"path": lc["path"], "sha256": lc["sha256"], "out_step_interval": DT, "start_step": 0, "W": W,
                                    "min_points": W // DT + 1, "hold": 0, "hold_source": "settle_criteria.yaml の hold", "first_judgeable": 2 * W + DT},
           "note": None, "judgements": [], "hold": {"hold_steps": 0, "start": 26000, "reached": [26000], "breaks": []},
           "request": ({"kind": "settled", "grid_step": 26000, "settled_since": 26000} if kind == "REACHED" else None),
           "forge_rc": 0, "run_case_rc": 0,
           "stopped_at": ({"present": True, "step": k, "sha256": SR._sha256(rd / SR.STOPPED_AT), "problems": []} if kind == "REACHED"
                          else {"present": False, "step": None, "sha256": None, "problems": []}),
           "reeval": reeval, "stop_sha256": SR._sha256(res_k) if res_k.exists() else None, "status_note": "試験", "failures": [],
           "reason": None, "stop_step": k, "last_output": k, "delta_E": de, "finished": SR._now()}
    SR.save_record(rd, rec)
    SR.write_state(rd, SR._v2_state(rd, rec, 0))
    return rd


def clone(src: Path, name: str) -> Path:
    dst = TMP / name
    shutil.copytree(src, dst, symlinks=True)
    return dst


def rewrite_records(rd: Path, mutate):
    """stop_rule.json を mutate(rec) で書き換えて保存し、STOP_RULE_STATE.json を _v2_state で作り直す (sha256 の結び付きを保つ)。"""
    rec = json.loads((rd / SR.RECORD).read_text())
    mutate(rec)
    SR.save_record(rd, rec)
    SR.write_state(rd, SR._v2_state(rd, rec, 0))


# ====================================================================================== 関門の呼び出し (偽物で抽出し直しを検出)
FORBIDDEN = []


def forbidden(name):
    def f(*a, **k):
        FORBIDDEN.append(name)
        raise RuntimeError(f"試験: 新方式の関門で {name} を呼んだ (抽出し直し・旧の判定)")
    return f


def stub_wr(kind="PASS"):
    return lambda run_dir, step=None: {"cmd": ["check_wall_resolution.py"], "rc": {"PASS": 0, "FAIL": 1}.get(kind, 2), "kind": kind,
                                       "verdict": [f"VERDICT: {kind}"], "output": f"VERDICT: {kind} (試験の偽物; step {step})\n"}


def stub_nw(ok=True):
    return lambda problem, delta_next: ({"pass": True, "wall": "試験の偽物"} if ok else {"pass": False, "why": "試験の偽物: 次の壁が非単調"})


REAL = {k: getattr(DL, k) for k in ("convergence_verdict", "wall_resolution_verdict", "next_wall_check", "extract_and_merge")}
REAL_SR = SR.extract_delta_e


def gate(rd, pin, wr="PASS", nw=True, eref=None, **kw):
    DL.wall_resolution_verdict, DL.next_wall_check = stub_wr(wr), stub_nw(nw)
    DL.convergence_verdict = forbidden("check_convergence (convergence_verdict)")
    DL.extract_and_merge = forbidden("extract_and_merge")
    SR.extract_delta_e = forbidden("extract_delta_e")
    try:
        return DL.handoff_gate(rd, problem="試験.yaml", euler_ref=(pin if eref is None else eref), **kw)
    finally:
        for k_, v_ in REAL.items():
            setattr(DL, k_, v_)
        SR.extract_delta_e = REAL_SR


def blocked_by(g, key):
    return any(key in b for b in g["blocked_by"])


def not_published(rd):
    return not (rd / "delta_r_next.csv").exists() and not (rd / "fixed_point.json").exists()


# ====================================================================================== (a) 関門の新方式
def test_gate(pin, eu):
    print("== (a) 関門の新方式", flush=True)
    A = build("A_reached", pin)
    FORBIDDEN.clear()
    g = gate(A, pin)
    saved = json.loads((A / DL.HANDOFF_GATE).read_text())
    k = 26001
    srcf = A / SR.V2_DIR / f"delta_E_{k}" / "delta_r_next.csv"
    check("A: REACHED + 再評価 SETTLED + 照合が全部一致 → 進む (method v2・settled true・上書きなし)",
          g["proceed"] and g["method"] == "v2" and saved["method"] == "v2" and g["settled"] is True and g["override"]["applied"] == []
          and not g["blocked_by"], "; ".join(g["blocked_by"])[:300])
    check("A: 関門の中で抽出し直さない・旧の判定を回さない (extract_and_merge・extract_delta_e・check_convergence を呼んでいない)",
          FORBIDDEN == [], str(FORBIDDEN))
    check("A: run の直下の delta_r_next.csv は停止規則の抽出 (stop_rule_v2/delta_E_26001) とバイト同一、HANDOFF_GATE の sha も同じ",
          (A / "delta_r_next.csv").read_bytes() == srcf.read_bytes() and g["extract"]["delta_r_next_sha256"] == SR._sha256(srcf)
          and g["extract"]["made"] and "delta_r_next.csv" in g["extract"]["made"])
    nx = DL.read_delta_r_next(A / "delta_r_next.csv")
    check("A: delta_r_next.csv は旧方式と同じ名前の列 (read_delta_r_next が named で読める)・使用域に held なし",
          nx["header_kind"] == "named" and not np.any(np.asarray(nx["held"])[(nx["x_rt"] >= 0) & (nx["x_rt"] <= X1)]))
    c = g["checks"]
    check("A: 各検査の記録 (照合・再評価の全体・残差・流量・目的の量・y₁⁺・δ_E・参照・ω・被覆・次の壁) が合格",
          all(c[n]["pass"] for n in ("binding", "reeval", "residual", "mass_balance", "targets", "wall_resolution", "delta_E", "euler_ref",
                                      "omega", "coverage", "next_wall")) and c["stop_step"] == k
          and c["reeval"]["by_quantity"] == {"residual": "SETTLED", "mass_balance": "SETTLED", "exit_core_mach": "SETTLED"},
          json.dumps({n: c[n].get("pass") for n in c if isinstance(c[n], dict)}, ensure_ascii=False))
    check("A: δ_E の作り方 = 停止規則 v2 の停止の step の 1 枚 (dir・場と表の sha256)、帯の感度の警告なし",
          g["delta_e_source"]["kind"] == "stop_rule_v2_stop_snapshot" and g["delta_e_source"]["stop_step"] == k
          and g["delta_e_source"]["dir"] == f"{SR.V2_DIR}/delta_E_{k}" and g["warnings"] == [], f"{g['warnings']}")
    summ = json.loads((A / "delta_r_extract_summary.json").read_text())
    check("A: 抽出の要約に落ち着きの印 (settled true・δ_E の作り方) を写す、収束の判定の記録は SETTLE_VERDICT.json の残差",
          summ["settled"] is True and summ["delta_e_source"]["kind"] == "stop_rule_v2_stop_snapshot" and g["convergence"]["kind"] == "PASS"
          and "SETTLE_VERDICT.json" in g["convergence"]["tool"] and "OVERALL: SETTLED" in (A / DL.HANDOFF_CONVERGENCE).read_text())

    # _after_ns (run_pass0_integral・run_from_table・run_pass の後段) も同じ表から
    A2 = clone(A, "A_after_ns")
    for f in ("delta_r_next.csv", DL.HANDOFF_GATE, "fixed_point.json", "delta_r_extract_summary.json"):
        (A2 / f).unlink(missing_ok=True)
    DL.wall_resolution_verdict, DL.next_wall_check = stub_wr(), stub_nw()
    DL.convergence_verdict, DL.extract_and_merge = forbidden("convergence_verdict"), forbidden("extract_and_merge")
    SR.extract_delta_e = forbidden("extract_delta_e")
    FORBIDDEN.clear()
    try:
        m = DL._after_ns(A2, pin, None, {}, "extract(new)", None, None, problem="試験.yaml")
    finally:
        for k_, v_ in REAL.items():
            setattr(DL, k_, v_)
        SR.extract_delta_e = REAL_SR
    fp = json.loads((A2 / "fixed_point.json").read_text())
    check("A: _after_ns は新方式の関門の表から fixed_point.json (settled true・δ_E の作り方) と delta_r_next.csv を作る (抽出し直さない)",
          m["handoff_gate"]["proceed"] and fp["settled"] is True and fp["delta_e_source"]["kind"] == "stop_rule_v2_stop_snapshot"
          and (A2 / "delta_r_next.csv").read_bytes() == srcf.read_bytes() and FORBIDDEN == [] and "massflow" in fp)

    # --- 判別の A/B と照合の負例 (上書きなし・全部の上書き) ---
    def mut_notice_fnv(rd):
        n = json.loads((rd / ST.MAIN_STAGE_RECORD).read_text())
        n["main_config"]["fnv"] = "deadbeefdeadbeef"
        (rd / ST.MAIN_STAGE_RECORD).write_text(json.dumps(n))

    def mut_criteria(rd):
        with open(rd / ST.CRITERIA_FILE, "a") as f:
            f.write("# 準備の後に書き換えた\n")

    def mut_stopped_at(rd):
        sa = json.loads((rd / SR.STOPPED_AT).read_text())
        sa["step"], sa["residual_step"] = 26003, 26002
        (rd / SR.STOPPED_AT).write_text(json.dumps(sa))

    def mut_delta_status(rd):
        rewrite_records(rd, lambda rec: rec["delta_E"].update(status="error", error="試験: 抽出の失敗"))

    def mut_no_launch(rd):
        (rd / "forge_launches.jsonl").unlink()

    def mut_rejudged(rd):
        # 後で check_settled.py --end-step 26000 を回して SETTLE_VERDICT.json を書き換えた (停止の step の再評価でない)
        lc = ST.load_criteria(rd)
        d = ST.RunData(rd, lc["criteria"])
        ST.write_verdict(rd, ST.judge(rd, 26000, data=d, loaded=lc), d)

    def mut_extract_settings(rd):
        rewrite_records(rd, lambda rec: rec["delta_E"]["extract"].update(band_select=repr("adaptive")))

    def mut_launch_other(rd):
        write_launch(rd, cfg_text="time: {last: {nStepOuter: 3000}}\n")    # 本段の後に別の設定で forge を起動した

    cases = [("B: MAIN_STAGE_STARTED.main_config.fnv だけ不一致 (判別の A/B)", mut_notice_fnv, "照合"),
             ("基準 settle_criteria.yaml の sha256 の不一致", mut_criteria, "照合"),
             ("STOPPED_AT の step の不一致", mut_stopped_at, "照合"),
             ("δ_E の status が ok でない (状態・見張りの記録とも)", mut_delta_status, "δ_E"),
             ("起動の記録 forge_launches.jsonl が無い (照合できない)", mut_no_launch, "照合"),
             ("SETTLE_VERDICT.json が停止の step の再評価でない (終端 26,000 で書き直された)", mut_rejudged, "照合"),
             ("δ_E の抽出の設定が受け渡しの固定の条件と違う (帯 adaptive)", mut_extract_settings, "δ_E"),
             ("最後の forge の起動が本段でない (本段の後に別の設定で起動)", mut_launch_other, "照合")]
    for i, (lab, mut, key) in enumerate(cases):
        rd = clone(A, f"neg_{i}")
        for f in ("delta_r_next.csv", DL.HANDOFF_GATE, "delta_r_extract_summary.json", "delta_r_equiv.csv", "delta_r_equiv_diag.json",
                  "delta_r_edge.csv"):
            (rd / f).unlink(missing_ok=True)
        mut(rd)
        g = gate(rd, pin)
        g2 = gate(rd, pin, **ALL_OV)
        check(f"負例 {lab} → 止める (上書きなし・上書き 3 つでも)、理由に「{key}」、delta_r_next.csv を公開しない",
              not g["proceed"] and blocked_by(g, key) and not g2["proceed"] and blocked_by(g2, key) and not_published(rd)
              and g["method"] == "v2", "; ".join(g["blocked_by"])[:260])
    check("判別の A/B: A (整合の入力) は進み B (fnv だけ不一致) は止める — 陽性の対照が成立", True)

    # --- (a) CENSORED・FAILED ---
    Cn = build("C_censored", pin, kind="CENSORED")
    g = gate(Cn, pin)
    check("CENSORED (上限 26,000 で停止) → 止める (停止の状態)、照合は合格", not g["proceed"] and blocked_by(g, "停止の状態 CENSORED")
          and g["checks"]["binding"]["pass"] and not_published(Cn), "; ".join(g["blocked_by"])[:200])
    g = gate(Cn, pin, allow_not_reached="試験: 上限で止めた場から作る")
    check("CENSORED + allow_not_reached → 進む、settled false・unsettled_overrides に停止の状態 (理由つき)",
          g["proceed"] and g["settled"] is False and g["override"]["applied"] == ["allow_not_reached"]
          and g["unsettled_overrides"][0]["condition"] == "stop_state" and (Cn / "delta_r_next.csv").exists(), "; ".join(g["blocked_by"])[:200])
    Fd = clone(A, "F_failed")
    for f in ("delta_r_next.csv", DL.HANDOFF_GATE):
        (Fd / f).unlink(missing_ok=True)
    rewrite_records(Fd, lambda rec: rec.update(status="FAILED", failures=["reeval_not_settled"], delta_E=None))
    g = gate(Fd, pin, **ALL_OV)
    check("FAILED (δ_E なし) → 上書き 3 つでも止める (停止の状態・δ_E の契約)", not g["proceed"] and blocked_by(g, "停止の状態 FAILED")
          and blocked_by(g, "δ_E") and not_published(Fd), "; ".join(g["blocked_by"])[:200])

    # --- (b)(c) 再評価の全体 (CENSORED の run で; REACHED は再評価が SETTLED でなければ停止規則が FAILED にする) ---
    Cr = build("C_resid", pin, kind="CENSORED", resid_fn=lambda s: 1e-4 * 10.0 ** (s / 20000.0), delta_from=Cn)
    g = gate(Cr, pin, allow_not_reached="試験")
    check("(b) 残差が上がり続ける (閾値の超過) → allow_not_reached だけでは止める (残差の再評価 NOT_SETTLED)",
          not g["proceed"] and blocked_by(g, "(b) 残差の再評価 NOT_SETTLED") and g["convergence"]["kind"] == "NOT_CONVERGED",
          "; ".join(g["blocked_by"])[:200])
    g = gate(Cr, pin, allow_not_reached="試験", allow_unconverged="試験: 残差は横ばいでない")
    check("(b) + allow_unconverged → 進む (settled false)", g["proceed"] and g["settled"] is False and "allow_unconverged" in g["override"]["applied"],
          "; ".join(g["blocked_by"])[:200])
    Cm = build("C_mass", pin, kind="CENSORED", mdot_fn=lambda s: (10.0, 9.99), delta_from=Cn)
    g = gate(Cm, pin, allow_not_reached="試験", allow_unconverged="試験")
    check("(c) 流量の収支 r = 1e−3 (> 5e−4) → allow_unsettled が無ければ止める", not g["proceed"] and blocked_by(g, "(c) 流量の収支"),
          "; ".join(g["blocked_by"])[:200])
    g = gate(Cm, pin, **ALL_OV)
    check("(c) + allow_unsettled → 進む、unsettled_overrides に流量の収支 (量と不合格の検査)",
          g["proceed"] and any(u["condition"] == "mass_balance" and u["quantities"] == {"mass_balance": "NOT_SETTLED"}
                               for u in g["unsettled_overrides"]), "; ".join(g["blocked_by"])[:200])
    Ct = build("C_target_missing", pin, kind="CENSORED", skip_nodes=(20000,), delta_from=Cn)
    g = gate(Ct, pin, **ALL_OV)
    check("目的の量の判定不能 (res_nodes_20000.h5 が無い) → 上書き 3 つでも止める", not g["proceed"] and blocked_by(g, "目的の量")
          and g["checks"]["targets"]["kind"] == "UNDETERMINED" and not_published(Ct), "; ".join(g["blocked_by"])[:200])

    # --- (d)(f)(g)(h) は新方式でも効く ---
    def fresh(name):
        rd = clone(A, name)
        for f in ("delta_r_next.csv", DL.HANDOFF_GATE, "delta_r_extract_summary.json"):
            (rd / f).unlink(missing_ok=True)
        return rd
    rd = fresh("d_y1")
    g = gate(rd, pin, wr="FAIL", **ALL_OV)
    check("(d) 局所 y₁⁺ FAIL → 上書きでも止める", not g["proceed"] and blocked_by(g, "局所 y₁⁺") and not_published(rd))
    rd = fresh("d_y1_und")
    g = gate(rd, pin, wr="UNDETERMINED", **ALL_OV)
    check("(d) 局所 y₁⁺ 判定不能 → 上書きでも止める", not g["proceed"] and blocked_by(g, "局所 y₁⁺"))
    rd = fresh("f_eref")
    g = gate(rd, pin, eref=str(eu), **ALL_OV)
    check("(f) 呼び出し側の Euler の参照が本段の記録 (固定の res_1000) と違う → 上書きでも止める",
          not g["proceed"] and not g["checks"]["euler_ref"]["pass"] and not_published(rd))
    rd = fresh("f_omega")
    g = gate(rd, pin, omega=0.4, **ALL_OV)
    check("(f) ω 0.4 (固定・抽出の 0.5 と違う) → 上書きでも止める", not g["proceed"] and not g["checks"]["omega"]["pass"] and not_published(rd))
    rd = fresh("g_held")
    dn = rd / SR.V2_DIR / f"delta_E_{k}" / "delta_r_next.csv"
    nx = DL.read_delta_r_next(dn)
    h = np.asarray(nx["held"]).copy()
    h[np.argmin(np.abs(nx["x_rt"] - 10.0))] = 1
    np.savetxt(dn, np.c_[nx["x_rt"], nx["delta_target"], nx["delta_in_prev"], nx["delta_E"], h], delimiter=",",
               header=",".join(DL.DELTA_R_NEXT_COLUMNS), comments="")
    new_sha = SR._sha256(dn)
    rewrite_records(rd, lambda rec: rec["delta_E"].update(delta_r_next_sha256=new_sha))
    g = gate(rd, pin, **ALL_OV)
    check("(g) 使用域の held 1 点 (表と記録の sha256 を揃えて δ_E の契約は通す) → 上書きでも止める (被覆)",
          not g["proceed"] and g["checks"]["delta_E"]["pass"] and blocked_by(g, "被覆") and g["checks"]["coverage"]["n_held_use"] == 1
          and not (rd / "delta_r_next.csv").exists(), "; ".join(g["blocked_by"])[:200])
    rd = fresh("h_wall")
    g = gate(rd, pin, nw=False, **ALL_OV)
    check("(h) 次の壁の幾何が不合格 → 上書きでも止める", not g["proceed"] and blocked_by(g, "次の壁") and not_published(rd))
    rd = fresh("e_sens_missing")
    (rd / SR.V2_DIR / f"delta_E_{k}" / "band_sensitivity.json").unlink()
    g = gate(rd, pin)
    check("帯の位置の感度の記録が欠けても合否は変えず、警告として HANDOFF_GATE.json に残す (「問題なし」と扱わない)",
          g["proceed"] and any("帯の位置の感度" in w for w in g["warnings"])
          and any("帯の位置の感度" in w for w in json.loads((rd / DL.HANDOFF_GATE).read_text())["warnings"]))
    rd = clone(A, "x_extractor")
    try:
        DL.handoff_gate(rd, problem="試験.yaml", euler_ref=pin, extractor=lambda *a, **k: {}, delta_e_source={"kind": "x", "n_snapshots": 5})
        ok = False
    except ValueError as e:
        ok = "停止規則 v2" in str(e)
    check("新方式で呼び出し側の extractor を渡すと ValueError (抽出し直さない)", ok)
    return A


# ====================================================================================== (b) 移行: 旧方式の run は今の関門
def test_migration(pin, A):
    print("== (b) 移行: 旧方式 (基準のファイルなし) は今の関門・表の形式は同じ", flush=True)
    X, R, Q = synth_ns()
    rd = TMP / "old_twin"
    rd.mkdir()
    s_stop = 62500
    write_field(rd, X, R, Q, steps=(s_stop,))
    info = ns_info(pin)
    sr = SR.resolve_settings({"stop_rule": True}, nsteps_arg=None, n_main=200000, n_main_source="試験", out_int=2500, euler_ref=pin)
    info["stop_rule"] = sr
    (rd / "prepare_info.json").write_text(json.dumps(info))
    uc = SR.use_columns(rd)
    rows = [{"step": st, "nonfinite_all": 0, "r": 0.0,
             "delta": {"delta_E": [1e-3] * len(uc["x"]), "delta_target": [2e-3] * len(uc["x"]), "held": [0] * len(uc["x"]),
                       "missing": [], "n_missing": 0, "n_held": 0, "seconds": 0.0}} for st in range(2500, s_stop + 1, 2500)]
    SR.save_record(rd, {"version": 1, "settings": sr, "status": "REACHED", "rows": rows, "stop_step": s_stop,
                        "delta_series": {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"],
                                         "extract": sr["extract"], "euler_ref": sr["euler_ref"], "span_steps": 40000}})
    SR.write_state(rd, {"method": "v1", "status": "REACHED", "stop_step": s_stop, "forge_rc": 143})
    DL.wall_resolution_verdict, DL.next_wall_check = stub_wr(), stub_nw()
    DL.convergence_verdict = lambda run_dir: {"cmd": ["check_convergence.py"], "rc": 0, "kind": "PASS", "verdict": ["OVERALL: PASS (偽物)"],
                                              "output": "OVERALL: PASS (試験の偽物)\n"}
    try:
        g = DL.handoff_gate(rd, problem="試験.yaml", euler_ref=pin)
    finally:
        for k_, v_ in REAL.items():
            setattr(DL, k_, v_)
    check("旧方式の run (settle_criteria.yaml なし) は今の関門: method v1・version 3・旧の検査 (保存則・2 窓) で進む",
          g["proceed"] and g["method"] == "v1" and g["version"] == 3 and "two_windows" in g["checks"] and "conservation" in g["checks"]
          and g["delta_e_source"]["kind"] == "single_snapshot", "; ".join(g["blocked_by"])[:200])
    a_old = (rd / "delta_r_next.csv").read_text().splitlines()
    a_new = (A / "delta_r_next.csv").read_text().splitlines()
    o, n = DL.read_delta_r_next(rd / "delta_r_next.csv"), DL.read_delta_r_next(A / "delta_r_next.csv")
    check("同じ場で旧方式の関門の抽出 (extract_and_merge) と新方式の表 (停止規則の extract_delta_e) が同じ形式: 見出しの行・行の数・"
          "全列の値がビット同一", a_old[0] == a_new[0] and len(a_old) == len(a_new)
          and all(np.array_equal(o[c], n[c], equal_nan=True) for c in DL.DELTA_R_NEXT_COLUMNS), a_old[0][:120])
    check("旧方式の delta_r_next.csv と新方式は同じファイルの名前の集合を run の直下に写す",
          set(g["extract"]["made"]) == set(json.loads((A / DL.HANDOFF_GATE).read_text())["extract"]["made"]),
          f"{g['extract']['made']}")


# ====================================================================================== (a′) 残差を外した基準 (wt_nozzle_v3) の関門
def test_gate_v3(pin, A):
    """plan §4.4「残差を止める判定から外す」(2026-10-16 ユーザ決定)・§5.1 #8: 基準に残差が無く非有限の見張り (residual_finite) がある run は
    関門の (b) が対象外 (止めない)、見張りの DIVERGED は上書きでも止める。対照として同じ残差を v2 の形の基準で判定すると (b) で止まる。"""
    print("== (a′) 残差を外した基準 (wt_nozzle_v3 の形) の関門", flush=True)
    osc = lambda s: 1e-4 * (10.0 ** 0.6 if (s // 7) % 2 else 1.0)       # 振幅 0.6 桁 (float の run_0015 の ρ の 0.585 桁ほど)
    V3 = build("A_v3_osc", pin, resid_fn=osc, delta_from=A, residual="watch")
    lc = ST.load_criteria(V3)
    kinds = [q["kind"] for q in lc["criteria"]["quantities"]]
    g = gate(V3, pin)
    c = g["checks"]
    check("v3: 残差の振幅 0.6 桁でも、再評価は SETTLED (残差は判定しない・非有限の見張りは SETTLED)・関門は上書きなしで進む、"
          "(b) 残差は対象外 (NOT_APPLICABLE)・(b′) 見張りは合格",
          lc["status"] == "ok" and "residual" not in kinds and "residual_finite" in kinds and g["proceed"] and g["settled"] is True
          and g["override"]["applied"] == [] and not g["blocked_by"] and c["residual"]["kind"] == "NOT_APPLICABLE"
          and c["residual_watch"]["pass"] and c["reeval"]["pass"] and g["convergence"]["kind"] == "NOT_APPLICABLE"
          and c["reeval"]["by_quantity"] == {"residual_finite": "SETTLED", "mass_balance": "SETTLED", "exit_core_mach": "SETTLED"},
          "; ".join(g["blocked_by"])[:300] + " " + json.dumps(c.get("residual"), ensure_ascii=False)[:200])
    check("v3: HANDOFF_CONVERGENCE.txt に (b) 対象外と (b′) 見張りを書く",
          "(b) 残差: NOT_APPLICABLE" in (V3 / DL.HANDOFF_CONVERGENCE).read_text()
          and "(b′) 残差の非有限の見張り: PASS" in (V3 / DL.HANDOFF_CONVERGENCE).read_text(), (V3 / DL.HANDOFF_CONVERGENCE).read_text()[-300:])
    V2c = build("C_v2_osc", pin, kind="CENSORED", resid_fn=osc)   # δ_E は本物の抽出 (CENSORED の 26,000; 下の v3 の非有限の run が写す)
    g = gate(V2c, pin, allow_not_reached="試験")
    check("v3 の対照: 同じ残差を v2 の形の基準 (残差を判定) で → (b) 残差の再評価 NOT_SETTLED で止まる",
          not g["proceed"] and blocked_by(g, "(b) 残差の再評価 NOT_SETTLED"), "; ".join(g["blocked_by"])[:200])
    nan_at = lambda s: float("nan") if s >= 25000 else 1e-4
    V3n = build("C_v3_nan", pin, kind="CENSORED", resid_fn=nan_at, delta_from=V2c, residual="watch")
    g = gate(V3n, pin, **ALL_OV)
    check("v3: 残差に非有限 (25,000 から) → 見張りが DIVERGED、上書き 3 つでも止める ((b′) 残差の非有限の見張り)",
          not g["proceed"] and blocked_by(g, "(b′) 残差の非有限の見張り DIVERGED") and g["checks"]["residual"]["kind"] == "NOT_APPLICABLE"
          and not_published(V3n), "; ".join(g["blocked_by"])[:300])


# ====================================================================================== (c) run_case.sh の道具の選択
FAKE_FORGE = TMP / "fake_forge.py"
FAKE_SRC = f'''#!{sys.executable}
import json, os, sys
sys.path.insert(0, {str(TOOLS)!r})
import stage_manifest as SM
if not os.environ.get("FAKE_NO_LAUNCH"):
    cfg = open("solverConfig.yaml", "rb").read()
    open("forge_launches.jsonl", "a").write(json.dumps({{"time": 0, "cfg_fnv": SM.fnv1a64(cfg)}}) + "\\n")
print("fake forge (run_settle_migration_tests)")
sys.exit(0)
'''


def run_case(rd, extra_env=None):
    env = {**os.environ, "FORGE_BIN": str(FAKE_FORGE), **(extra_env or {})}
    p = subprocess.run(["bash", str(TOOLS / "run_case.sh"), str(rd)], capture_output=True, text=True, env=env, timeout=600)
    t = (rd / "CONVERGENCE_VERDICT.txt").read_text() if (rd / "CONVERGENCE_VERDICT.txt").exists() else ""
    return p.returncode, t


def tool_line(t):
    return t.splitlines()[0] if t else ""


def exit_code(t):
    ln = [x for x in t.splitlines() if x.startswith("CONVERGENCE_EXIT:")]
    return int(ln[-1].split(":")[1]) if ln else None


def py_wrapper(name, script):
    d = TMP / name
    d.mkdir()
    real = shutil.which("python3")
    (d / "python3").write_text("#!/bin/sh\n" + script.replace("@REAL@", real))
    (d / "python3").chmod(0o755)
    return {"PATH": f"{d}:{os.environ['PATH']}"}


def test_run_case(A):
    print("== (c) run_case.sh の道具の選択 (偽の forge)", flush=True)
    FAKE_FORGE.write_text(FAKE_SRC)
    FAKE_FORGE.chmod(0o755)

    def fresh(name):
        rd = clone(A, name)
        (rd / "CONVERGENCE_VERDICT.txt").unlink(missing_ok=True)
        return rd
    rd = fresh("rc_main")
    rc, t = run_case(rd)
    check("本段 (MAIN_STAGE_STARTED と最後の起動・設定・基準が一致) + 基準 → check_settled.py --bind: 1 行目 TOOL、SETTLED、終了コード 0",
          tool_line(t).startswith("TOOL: check_settled.py (本段の起動") and "VERDICT: SETTLED" in t and exit_code(t) == 0
          and "照合:" in t and "OVERALL" not in t, t[:300])
    rd = fresh("rc_pre")
    (rd / ST.MAIN_STAGE_RECORD).unlink()
    (rd / "stage_manifest.json").write_text(json.dumps({"manifest_version": 3, "stages": [{"order": 0, "tag": "S1_soft"}]}))
    rc, t = run_case(rd)
    check("前の段と確かめられた run (MAIN_STAGE_STARTED なし・stage_manifest.json の段 S1_soft に本段なし) → check_convergence.py (1 行目 TOOL)",
          tool_line(t).startswith("TOOL: check_convergence.py (本段の開始の記録") and "前の段と確かめた" in tool_line(t)
          and "CONVERGENCE_EXIT:" in t and "VERDICT: SETTLED" not in t, t[:200])
    rd = fresh("rc_pre_unconfirmed")
    (rd / ST.MAIN_STAGE_RECORD).unlink()
    (rd / "stage_manifest.json").unlink(missing_ok=True)
    rc, t = run_case(rd)
    check("MAIN_STAGE_STARTED が無く stage_manifest.json も無い (前の段と確かめられない) → check_settled.py の INDETERMINATE (記録の欠け)・"
          "終了コード 3 (以前は check_convergence.py; result 段レビュー M1)",
          tool_line(t).startswith("TOOL: check_settled.py") and "VERDICT: INDETERMINATE (本段の開始の記録の欠け)" in t and exit_code(t) == 3,
          t[:300])
    rd = fresh("rc_legacy")
    (rd / ST.CRITERIA_FILE).unlink()
    rc, t = run_case(rd)
    check("基準のファイルの無い run → check_convergence.py (今どおり; 1 行目 TOOL)",
          tool_line(t) == "TOOL: check_convergence.py (settle_criteria.yaml が無い)" and "CONVERGENCE_EXIT:" in t, t[:200])
    rd = fresh("rc_not_main")
    with open(rd / "solverConfig.yaml", "a") as f:
        f.write("# 本段の後に別の設定で起動\n")
    rc, t = run_case(rd)
    check("最後の起動の cfg_fnv が本段の記録と違う (以前は「本段でない起動」で check_convergence.py) → check_settled.py の INDETERMINATE "
          "(照合の不一致)・終了コード 3 (result 段レビュー M1)",
          tool_line(t).startswith("TOOL: check_settled.py (本段の起動") and "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" in t
          and "最後の forge の起動" in t and exit_code(t) == 3, t[:300])
    rd = fresh("rc_no_launch")
    (rd / "forge_launches.jsonl").unlink(missing_ok=True)
    rc, t = run_case(rd, {"FAKE_NO_LAUNCH": "1"})
    check("本段の記録があり起動の記録 forge_launches.jsonl が無い (設定・基準は記録と一致) → check_settled.py の INDETERMINATE "
          "(照合の不一致; 以前は問題なしで判定)", tool_line(t).startswith("TOOL: check_settled.py")
          and "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" in t and "起動の記録" in t and exit_code(t) == 3, t[:300])
    rd = fresh("rc_bind_bad")
    with open(rd / "solverConfig.yaml", "a") as f:
        f.write("# 起動の記録の無い forge で、本段の記録と違う設定\n")
    (rd / "forge_launches.jsonl").unlink()
    rc, t = run_case(rd, {"FAKE_NO_LAUNCH": "1"})
    check("起動の記録が無く設定が本段の記録と違う → check_settled.py の INDETERMINATE (照合の不一致; 旧の道具へ落とさない)・終了コード 3、"
          "SETTLE_VERDICT.json を書き換えない",
          tool_line(t).startswith("TOOL: check_settled.py") and "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" in t and exit_code(t) == 3
          and (rd / ST.VERDICT_FILE).read_bytes() == (A / ST.VERDICT_FILE).read_bytes(), t[:300])
    rd = fresh("rc_crit_bad")
    with open(rd / ST.CRITERIA_FILE, "a") as f:
        f.write("# 準備の後に書き換えた\n")
    rc, t = run_case(rd)
    check("基準の sha256 が本段の記録と違う → check_settled.py の INDETERMINATE・終了コード 3",
          tool_line(t).startswith("TOOL: check_settled.py") and "INDETERMINATE (本段の起動との照合の不一致)" in t and exit_code(t) == 3, t[:300])
    rd = fresh("rc_invalid")
    c = yaml.safe_load((rd / ST.CRITERIA_FILE).read_text())
    c["quantities"][0]["thresholds"]["bogus"] = 1
    (rd / ST.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    rc, t = run_case(rd)
    check("基準のファイルの誤り (知らないキー) → check_settled.py の INDETERMINATE (check_convergence.py へ落とさない)",
          tool_line(t).startswith("TOOL: check_settled.py") and "VERDICT: INDETERMINATE" in t and exit_code(t) == 3
          and "OVERALL" not in t, t[:300])
    rd = fresh("rc_manual")
    (rd / ST.MAIN_STAGE_RECORD).unlink()
    c = yaml.safe_load((rd / ST.CRITERIA_FILE).read_text())
    c["grid"]["segment"] = None
    (rd / ST.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    rc, t = run_case(rd)
    check("手で組んだ run (基準が本段の記録を指さない・MAIN_STAGE_STARTED なし) → check_settled.py で判定 (SETTLED)",
          tool_line(t).startswith("TOOL: check_settled.py (手で組んだ run") and "VERDICT: SETTLED" in t and exit_code(t) == 0, t[:300])
    rd = fresh("rc_select_fail")
    env = py_wrapper("py_select_fail", 'case "$*" in *--select*) exit 1;; esac\nexec @REAL@ "$@"\n')
    rc, t = run_case(rd, env)
    check("道具の選択に失敗 (check_settled.py --select が落ちる) → 基準のある run は check_settled.py (旧の道具へ落とさない)",
          tool_line(t).startswith("TOOL: check_settled.py (道具の選択に失敗") and "VERDICT: SETTLED" in t, t[:300])
    rd = fresh("rc_no_verdict")
    env = py_wrapper("py_settled_fail", 'case "$1" in *check_settled.py) echo "Traceback (試験)"; exit 1;; esac\nexec @REAL@ "$@"\n')
    rc, t = run_case(rd, env)
    check("check_settled.py が VERDICT を出さずに終わる → VERDICT: INDETERMINATE を足し終了コード 3 (1 = NOT_SETTLED と取り違えない)",
          "VERDICT: INDETERMINATE (check_settled.py が VERDICT を出さずに終わった; 終了コード 1)" in t and exit_code(t) == 3, t[-300:])
    # 停止規則の後処理 (run_case.sh が後処理をしなかったときの代わり) も同じ選択
    rd = fresh("pp_main")
    write_launch(rd)
    pp = SR.postprocess(rd)
    t = (rd / "CONVERGENCE_VERDICT.txt").read_text()
    check("stop_rule.postprocess: 本段 + 基準 → check_settled.py --bind (1 行目 TOOL・SETTLED・CONVERGENCE_EXIT 0)",
          pp["tool"] == "check_settled.py" and pp["convergence_exit"] == 0 and tool_line(t).startswith("TOOL: check_settled.py")
          and "VERDICT: SETTLED" in t, t[:200])
    rd = fresh("pp_legacy")
    (rd / ST.CRITERIA_FILE).unlink()
    pp = SR.postprocess(rd)
    t = (rd / "CONVERGENCE_VERDICT.txt").read_text()
    check("stop_rule.postprocess: 基準なし → check_convergence.py (今どおり)", pp["tool"] == "check_convergence.py"
          and tool_line(t) == "TOOL: check_convergence.py (settle_criteria.yaml が無い)" and "CONVERGENCE_EXIT" in t)
    return TMP / "rc_main", TMP / "rc_legacy"


# ====================================================================================== (d) Stop フック
def test_hook(rc_main: Path, rc_legacy: Path):
    print("== (d) Stop フック (hook_convergence_gate.py)", flush=True)
    root = TMP / "hookroot"
    (root / "solver_density_cuda" / "tools").mkdir(parents=True)
    for f in ("hook_convergence_gate.py", "forge_gate_claims.py"):
        shutil.copy2(TOOLS / f, root / "solver_density_cuda" / "tools" / f)
    case = root / "case" / "99.hooktest"
    now = time.time()

    def mk(name, verdict_src: Path | None, verdict_age: float | None):
        d = case / name
        d.mkdir(parents=True)
        (d / "residual_history.csv").write_text("step,inner,phase,rms_ro\n0,-1,outer_end,1e-4\n")
        os.utime(d / "residual_history.csv", (now - 600, now - 600))
        if verdict_src is not None:
            shutil.copy2(verdict_src, d / "CONVERGENCE_VERDICT.txt")
            os.utime(d / "CONVERGENCE_VERDICT.txt", (now - verdict_age, now - verdict_age))
        return d

    def hook():
        p = subprocess.run([sys.executable, str(root / "solver_density_cuda" / "tools" / "hook_convergence_gate.py")], input="{}",
                           capture_output=True, text=True, timeout=60)
        return json.loads(p.stdout.strip() or "{}")

    mk("run_0001_settled", rc_main / "CONVERGENCE_VERDICT.txt", 300)
    mk("run_0002_convergence", rc_legacy / "CONVERGENCE_VERDICT.txt", 300)
    out = hook()
    check("check_settled.py・check_convergence.py の両方の形式の CONVERGENCE_VERDICT.txt (residual より新しい) で block しない",
          out == {} and "TOOL: check_settled.py" in (rc_main / "CONVERGENCE_VERDICT.txt").read_text()
          and "TOOL: check_convergence.py" in (rc_legacy / "CONVERGENCE_VERDICT.txt").read_text(), str(out)[:200])
    mk("run_0003_none", None, None)
    out = hook()
    check("CONVERGENCE_VERDICT.txt が無い run があれば今どおり block (その run だけを挙げる)",
          out.get("decision") == "block" and "run_0003_none" in out.get("reason", "") and "run_0001" not in out.get("reason", ""),
          str(out)[:200])
    check("block の案内文: 道具の選択 (check_settled.py --select・--bind / check_convergence.py) と「通過は収束の合格ではない」",
          "check_settled.py --select" in out.get("reason", "") and "--bind" in out.get("reason", "")
          and "収束の合格ではない" in out.get("reason", ""))
    shutil.rmtree(case / "run_0003_none")
    mk("run_0004_old", rc_main / "CONVERGENCE_VERDICT.txt", 900)
    out = hook()
    check("CONVERGENCE_VERDICT.txt が residual より古い → 今どおり block", out.get("decision") == "block" and "run_0004_old" in out.get("reason", ""),
          str(out)[:200])


# ====================================================================================== (e) #7b 残りの読み手
_TMPDIR_RE = __import__("re").compile(r"/tmp/thermal_use_(conv|qs)_[A-Za-z0-9_]+")


def _norm(o) -> str:
    return _TMPDIR_RE.sub(lambda m: f"/tmp/thermal_use_{m.group(1)}_<TMP>", json.dumps(o, ensure_ascii=False, default=str))


def _load_ref(relpath: str, modname: str, commit: str = "HEAD"):
    """変更前 (既定 HEAD; 検出力の確かめは固定の commit) のモジュールを git から読み、同じパッケージの中の名前で読み込む
    (相対 import を同じに解く)。"""
    import importlib.util
    src = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{relpath}"], capture_output=True, text=True, check=True).stdout
    f = TMP / "ref_tree" / Path(relpath).parent / (modname.rsplit(".", 1)[-1] + ".py")   # 同じ深さに置く (REPO = parents[3])
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(src)
    spec = importlib.util.spec_from_file_location(modname, f)
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    for k_ in ("REPO", "TOOLS"):                          # 道具は今の作業ツリーのもの (比べるのは wall_tables の中身だけ)
        if hasattr(m, k_):
            setattr(m, k_, getattr(sys.modules[modname.rsplit(".", 1)[0] + "." + Path(relpath).stem], k_))
    return m


PRE_REVIEW2_COMMIT = "b64dbfdd"     # result 段レビュー 2 回目 M3 (§5.1 #11) の修正の前 (検出力の確かめに使う)


def _rebind_verdict(rd: Path, mutate):
    """停止規則の再評価の記録 SETTLE_VERDICT.json を mutate(v) で書き換え、STOP_RULE_STATE の reeval (digest・sha256・verdict) を合わせる
    (別の機械・別の所在で書いた記録の形; 記録どうしは結び付いたまま)。"""
    p = rd / ST.VERDICT_FILE
    v = json.loads(p.read_text())
    mutate(v)
    p.write_text(json.dumps(v, ensure_ascii=False, indent=1) + "\n")
    dig, sha = SR.quantities_digest(v), SR._sha256(p)
    rewrite_records(rd, lambda rec: rec["reeval"].update(digest=dig, sha256=sha, verdict=v["verdict"]))


def test_reeval_match(V: Path):
    """result 段レビュー 2 回目 M3 (plan §5.1 #11): 成果物の収束 (`wall_tables._settled_convergence`) が停止規則 v2 の再評価の記録と
    今の判定を照合するとき、所在のパス (相対・絶対・別の場所) と演算の丸め (1 ulp) では一致、各検査の合否・量の verdict・許容を超える値・
    入力の中身が違えば不一致 (未評価)。修正前 (PRE_REVIEW2_COMMIT の wall_tables; 量の辞書全体の sha256 の完全一致) が所在と丸めで
    未評価だったこと (検出力) も確かめる。"""
    print("== (e′) 成果物の収束と停止規則の再評価の照合 (result 段レビュー 2 回目 M3、§5.1 #11)", flush=True)
    from forge_design.export import wall_tables as WT
    try:
        PRE = _load_ref("design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_review2", commit=PRE_REVIEW2_COMMIT)
    except Exception as e:  # noqa: BLE001
        PRE = None
        check(f"修正前 ({PRE_REVIEW2_COMMIT}) の wall_tables を読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))

    def conv(run, mod=WT):
        return mod.thermal_use(run)["convergence"]

    # 照合の部品の単体 (quantities_mismatch)
    a = [{"name": "q", "verdict": "SETTLED", "checks": [{"check": "c", "value": 0.1, "limit": 1.0, "ok": True, "detail": "x 0.1"}],
          "values": {"file": "/a/run/residual_history.csv", "n": 3, "v": 1e-9, "w": 0.5, "z": 3.3e-16, "nan": float("nan")},
          "reasons": [{"category": "data", "label": "データ不足", "text": "/a/run/x が無い"}]}]
    b = json.loads(json.dumps(SR._settle()._jsonable(a)))
    b[0]["values"]["file"] = "rel/run/residual_history.csv"
    b[0]["checks"][0]["detail"] = "x 0.1000000001"
    b[0]["reasons"][0]["text"] = "rel/run/x が無い"
    b[0]["checks"][0]["value"] = math.nextafter(0.1, 1.0)
    b[0]["values"]["v"] = 1e-9 * (1 + 5e-10)
    b[0]["values"]["w"] = 0.5 * (1 + 5e-10)
    b[0]["values"]["z"] = 0.0
    check("照合の単体: 所在のパス (file)・表示の文 (detail・text) の違い、1 ulp・相対 5e−10 の差、丸めの大きさの量 (3e−16 と 0)、"
          "非有限 (nan と nan) は一致", SR.quantities_mismatch(a, b) == [], str(SR.quantities_mismatch(a, b)))
    c_ = json.loads(json.dumps(b))
    c_[0]["values"]["v"] = 1e-9 * (1 + 1e-6)
    check(f"照合の単体: 絶対の床 (REEVAL_ATOL {SR.REEVAL_ATOL:g}) より小さい差は一致とみなす (1e−9 の量の相対 1e−6 = 絶対 1e−15; "
          "許容の決め方を試験に残す — 入力の中身は sha256 で別に揃える)", SR.quantities_mismatch(a, c_) == [], str(SR.quantities_mismatch(a, c_)))
    neg = []
    for label, mut in (("合否 ok", lambda x: x[0]["checks"][0].update(ok=False)), ("量の verdict", lambda x: x[0].update(verdict="NOT_SETTLED")),
                       ("整数", lambda x: x[0]["values"].update(n=4)), ("相対 1e−6 の差 (0.5 の量)", lambda x: x[0]["values"].update(w=0.5 * (1 + 1e-6))),
                       ("キーの欠け", lambda x: x[0]["values"].pop("n")), ("検査の名前", lambda x: x[0]["checks"][0].update(check="d")),
                       ("理由の分類", lambda x: x[0]["reasons"][0].update(category="threshold")),
                       ("真偽と整数", lambda x: x[0]["checks"][0].update(ok=1))):
        c_ = json.loads(json.dumps(b))
        mut(c_)
        if not SR.quantities_mismatch(a, c_):
            neg.append(label)
    check("照合の単体: 合否 ok・量の verdict・整数・許容を超える差 (相対 1e−6)・キーの欠け・検査の名前・理由の分類・真偽と整数の取り違えは不一致",
          not neg, str(neg))

    # 相対のパス (同じ run を作業ディレクトリからの相対で渡す)
    cwd = os.getcwd()
    try:
        os.chdir(V.parent)
        c = conv(Path(V.name))
        c_old = conv(Path(V.name), PRE) if PRE is not None else None
    finally:
        os.chdir(cwd)
    check("同じ run を相対のパスで渡す → 停止規則の再評価と一致・収束 PASS",
          c["result"] == "PASS" and c["stop_rule_reeval"]["matches"] is True, json.dumps(c.get("stop_rule_reeval"), ensure_ascii=False)[:300])
    check("検出力: 修正前は相対のパスで「未評価」(量の sha256 に所在のパスが入る)",
          c_old is not None and c_old["result"] == "未評価" and c_old["stop_rule_reeval"]["matches"] is False, str((c_old or {}).get("reason")))
    c = conv(V.resolve())
    check("同じ run を絶対のパスで渡す → 一致・収束 PASS", c["result"] == "PASS" and c["stop_rule_reeval"]["matches"] is True)
    # 別の場所 (記録を書いた run と違う所在; AWS からの写しの形)
    Vm = clone(V, "wt_v2_moved")
    c = conv(Vm)
    check("run を別の場所へ写す (記録の所在と違う) → 一致・収束 PASS (入力の中身は記録の相対のパスで引き直して sha256 が一致)",
          c["result"] == "PASS" and c["stop_rule_reeval"]["matches"] is True, json.dumps(c["stop_rule_reeval"].get("compare"), ensure_ascii=False)[:300])
    if PRE is not None:
        c_old = conv(Vm, PRE)
        check("検出力: 修正前は別の場所へ写した run で「未評価」", c_old["result"] == "未評価", str(c_old.get("reason")))
    # 演算の丸め (記録の浮動小数をすべて 1 ulp ずらし、記録どうしの結び付きは保つ)
    V1 = clone(V, "wt_v2_ulp")
    nb = [0]

    def bump(o):
        if isinstance(o, dict):
            return {k: bump(x) for k, x in o.items()}
        if isinstance(o, list):
            return [bump(x) for x in o]
        if isinstance(o, float) and math.isfinite(o) and o != 0.0:
            nb[0] += 1
            return math.nextafter(o, math.inf)
        return o
    _rebind_verdict(V1, lambda v: v.update(quantities=bump(v["quantities"])))
    c = conv(V1)
    check(f"記録の浮動小数を全部 1 ulp ずらす ({nb[0]} 個; 別の数値ライブラリで書いた記録の形) → 一致・収束 PASS",
          nb[0] > 10 and c["result"] == "PASS" and c["stop_rule_reeval"]["matches"] is True, json.dumps(c["stop_rule_reeval"].get("compare"), ensure_ascii=False)[:300])
    if PRE is not None:
        c_old = conv(V1, PRE)
        check("検出力: 修正前は 1 ulp の差で「未評価」", c_old["result"] == "未評価", str(c_old.get("reason")))

    # 不一致の負例 (未評価にする)
    def first_check(v, kind="exit_core_mach"):
        return next(q for q in v["quantities"] if q["kind"] == kind)["checks"][0]
    for name, mut, word in (
            ("ok", lambda v: first_check(v).update(ok=False), ".ok"),
            ("qverdict", lambda v: next(q for q in v["quantities"] if q["kind"] == "exit_core_mach").update(verdict="NOT_SETTLED"), ".verdict"),
            ("value", lambda v: first_check(v).update(value=first_check(v)["value"] * (1 + 1e-6) + 1e-9), ".value"),
            ("verdict", lambda v: v.update(verdict="NOT_SETTLED"), "verdict")):
        Vx = clone(V, f"wt_v2_neg_{name}")
        _rebind_verdict(Vx, mut)
        c = conv(Vx)
        check(f"記録の {name} だけが今の判定と違う (記録どうしは結び付いたまま) → 不一致・収束 未評価 (理由に食い違いの箇所)",
              c["result"] == "未評価" and c["stop_rule_reeval"]["matches"] is False and "食い違う" in (c.get("reason") or "")
              and word in (c.get("reason") or ""), str(c.get("reason"))[:300])
    Vi = clone(V, "wt_v2_input")
    p = Vi / "boundary_massflux.csv"
    p.write_text("# 再評価の後に足した行 (判定には効かない)\n" + p.read_text())
    c = conv(Vi)
    check("入力の中身が記録と違う (判定に効かない注記の行を足した流量の CSV) → 判定は同じでも不一致・収束 未評価 (入力の sha256 の照合を残す)",
          c["result"] == "未評価" and c["stop_rule_reeval"]["matches"] is False and "boundary_massflux.csv" in (c.get("reason") or ""),
          str(c.get("reason"))[:300])
    Vf = clone(V, "wt_v2_nofiles")
    _rebind_verdict(Vf, lambda v: v.pop("files"))
    c = conv(Vf)
    check("記録に使ったファイルの一覧が無い (旧い記録の形) → 例外にせず理由つきの未評価", c["result"] == "未評価"
          and "files" in (c.get("reason") or ""), str(c.get("reason"))[:300])
    Vv = clone(V, "wt_v2_noverdict")
    (Vv / ST.VERDICT_FILE).unlink()
    c = conv(Vv)
    check("再評価の記録 SETTLE_VERDICT.json が無い → 例外にせず理由つきの未評価", c["result"] == "未評価"
          and ST.VERDICT_FILE in (c.get("reason") or ""), str(c.get("reason"))[:300])


def _old_gate_row(run: Path):
    """変更前の nozzle_report.conditions の収束の行 (同じ式を写したもの)。"""
    cv = run / "CONVERGENCE_VERDICT.txt"
    if not cv.exists():
        return None
    t = cv.read_text().splitlines()
    overall = next((l.strip() for l in t if "OVERALL" in l), "?")
    worst = [l.strip() for l in t if l.strip().startswith("rms_") and "<--" in l][:3]
    return ["収束 (check_convergence)", overall + (" — " + "; ".join(worst) if worst else "")]


def test_readers():
    print("== (e) #7b 新方式の run を読む残りの読み手 (成果物の thermal_use・設計の報告の収束の欄・forge の直接実行の案内)", flush=True)
    sys.path.insert(0, str(DESIGN / "tests"))
    import run_wall_tables_tests as RWT
    from forge_design.export import deliverables as DLV
    from forge_design.export import wall_tables as WT
    from forge_design.report import nozzle_report as NR
    ref = _load_ref("design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_ref")
    wsteps = (18000, 20000, 22000, 24000, 26000, 26001)

    # --- 旧方式 (基準なし): thermal_use は変更前とビット同一 (一時ディレクトリの道だけ置き換え) ---
    L = TMP / "wt_legacy"
    RWT.make_synth(L, steps=(100, 200, 300, 400, 500), tw=1e-5, residual=True)
    new, old = WT.thermal_use(L), ref.thermal_use(L)
    check("thermal_use (旧方式の run、基準なし): 変更前の wall_tables と出力がビット同一 (一時ディレクトリの道を除く)・check_convergence.py・PASS",
          _norm(new) == _norm(old) and new["convergence"]["tool"] == "check_convergence.py" and new["verdict"] == "PASS",
          f"{new['verdict']} / {old['verdict']}")
    Lf = TMP / "wt_legacy_fail"
    RWT.make_synth(Lf, steps=(100, 200, 300), tw=1e-5, residual=False)
    with open(Lf / "residual_history.csv", "w") as f:
        f.write("step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n")
        for i in range(300):
            f.write(f"{i},-1,outer_end,1e-4,1e-4,1e-4,0,1e-4\n")
    new, old = WT.thermal_use(Lf), ref.thermal_use(Lf)
    check("thermal_use (旧方式の run、残差が下がらない): 変更前とビット同一 (収束 FAIL)",
          _norm(new) == _norm(old) and new["convergence"]["result"] == old["convergence"]["result"] == "FAIL")

    # --- 新方式 (本段 + 基準): check_settled.py --bind・停止規則の再評価と照合 ---
    V = TMP / "wt_v2"
    RWT.make_synth(V, steps=wsteps, tw=1e-5, residual=False)
    finish_v2(V, kind="REACHED", delta=None)
    th = WT.thermal_use(V)
    c = th["convergence"]
    check("thermal_use (新方式の run): 収束は check_settled.py --bind (終端 26,000・停止の観測 26,001)、SETTLED → PASS、停止規則 v2 の"
          "再評価と量の判定が一致、check_convergence を回さない",
          c["tool"] == "check_settled.py --bind" and c["result"] == "PASS" and c["judged"]["end_step"] == 26001
          and c["judged"]["terminal_step"] == 26000 and c["stop_rule_reeval"]["matches"] is True
          and "VERDICT: SETTLED" in c["verdict_line"] and "[check_settled.py --bind、終端 26000・停止の観測 26001]" in c["verdict_line"]
          and "check_convergence" not in json.dumps(c, ensure_ascii=False) and c["selection"]["mode"] == "main",
          json.dumps({k: c.get(k) for k in ("result", "verdict_line", "reason")}, ensure_ascii=False))
    check(f"thermal_use (新方式の run): 準定常・壁解像は今どおり ({th['quasisteady']['result']}・{th['wall_resolution']['result']}) で総合 {th['verdict']}",
          th["verdict"] == "PASS" and th["quasisteady"]["tool"] == "check_quasisteady.py --series-csv")
    check("thermal_use (新方式の run): run ディレクトリに書かない (SETTLE_VERDICT.json は停止規則の再評価のまま)",
          json.loads((V / ST.VERDICT_FILE).read_text())["end_step"] == 26001
          and SR._sha256(V / ST.VERDICT_FILE) == SR.read_state(V)["reeval"]["sha256"])
    # 成果物 (deliverables.part_excel): status.json・xlsx の注記に新方式の判定と道具が載る
    out = TMP / "deliv_v2"
    rx = DLV.part_excel(V, out)
    sj = json.loads((out / "excel" / "wall_cfd_nodes.status.json").read_text())
    tc = sj["thermal_use"]["convergence"]
    check("成果物 (part_excel): wall_cfd_nodes.status.json の thermal_use の収束が check_settled.py の判定 (道具の名前つき)・PASS",
          tc["result"] == "PASS" and "check_settled.py --bind" in (tc["verdict_line"] or "") and sj["use_for_thermal_design"] is True
          and rx["thermal_use"]["convergence"]["tool"] == "check_settled.py --bind", json.dumps(tc, ensure_ascii=False))
    # 照合の不一致・再評価との食い違い・未収束・前の段
    Vb = clone(V, "wt_v2_bind")
    with open(Vb / ST.CRITERIA_FILE, "a") as f:
        f.write("# 準備の後に書き換えた\n")
    c = WT.thermal_use(Vb)["convergence"]
    check("thermal_use: 基準の sha256 が本段の記録と違う → check_settled.py の INDETERMINATE → 収束 未評価 (check_convergence へ落とさない)",
          c["tool"] == "check_settled.py --bind" and c["result"] == "未評価" and "照合の不一致" in (c.get("reason") or ""), str(c.get("reason")))
    Vr = clone(V, "wt_v2_reeval")
    rewrite_records(Vr, lambda rec: rec["reeval"].update(digest="0" * 64))
    c = WT.thermal_use(Vr)["convergence"]
    check("thermal_use: 停止規則 v2 の再評価 (STOP_RULE_STATE の digest) と今の判定が食い違う → 未評価",
          c["result"] == "未評価" and c["stop_rule_reeval"]["matches"] is False and "食い違う" in (c.get("reason") or ""), str(c.get("reason")))
    test_reeval_match(V)
    Vn = TMP / "wt_v2_rising"
    RWT.make_synth(Vn, steps=wsteps, tw=1e-5, residual=False)
    finish_v2(Vn, kind="REACHED", resid_fn=lambda s: 1e-4 * 10.0 ** (s / 20000.0), delta=None)
    th = WT.thermal_use(Vn)
    check("thermal_use: 残差が上がり続ける (NOT_SETTLED) → 収束 FAIL・総合 FAIL",
          th["convergence"]["result"] == "FAIL" and th["verdict"] == "FAIL" and "NOT_SETTLED" in th["convergence"]["verdict_line"])
    Vp = clone(V, "wt_v2_premain")
    (Vp / ST.MAIN_STAGE_RECORD).unlink()
    (Vp / "stage_manifest.json").write_text(json.dumps({"manifest_version": 3, "stages": [{"order": 0, "tag": "S1_soft"}]}))
    c = WT.thermal_use(Vp)["convergence"]
    check("thermal_use: 前の段と確かめられた run (MAIN_STAGE_STARTED なし・stage_manifest.json に本段なし) → 今どおり check_convergence.py "
          "(選択の記録つき)", c["tool"] == "check_convergence.py" and c["selection"]["mode"] == "pre_main")
    Vu = clone(V, "wt_v2_unconfirmed")
    (Vu / ST.MAIN_STAGE_RECORD).unlink()
    (Vu / "stage_manifest.json").unlink(missing_ok=True)
    c = WT.thermal_use(Vu)["convergence"]
    check("thermal_use: MAIN_STAGE_STARTED が無く前の段と確かめられない → check_settled.py の INDETERMINATE → 収束 未評価 "
          "(check_convergence へ落とさない; result 段レビュー M1)", c["tool"] == "check_settled.py --bind" and c["result"] == "未評価"
          and c["selection"]["mode"] == "unconfirmed", str(c.get("reason")))

    # --- 設計の報告 (nozzle_report.convergence_gate_row): 道具に応じて VERDICT を読む ---
    Vc = clone(V, "rep_v2")
    write_launch(Vc)
    SR.postprocess(Vc)
    row = NR.convergence_gate_row(Vc)
    check("報告の収束の欄 (新方式の run、CONVERGENCE_VERDICT.txt の 1 行目 TOOL: check_settled.py) → 見出し「収束 (check_settled)」・VERDICT: SETTLED",
          row == ["収束 (check_settled)", "VERDICT: SETTLED"], str(row))
    Vcn = clone(Vn, "rep_v2_rising")
    write_launch(Vcn)
    SR.postprocess(Vcn)
    row = NR.convergence_gate_row(Vcn)
    check("報告の収束の欄 (新方式・NOT_SETTLED) → VERDICT: NOT_SETTLED と SETTLED でない量", row[0] == "収束 (check_settled)"
          and row[1].startswith("VERDICT: NOT_SETTLED") and "[NOT_SETTLED  ] residual" in row[1], str(row))
    Lc = clone(L, "rep_legacy")
    SR.postprocess(Lc)                                   # 1 行目 TOOL: check_convergence.py の形
    row = NR.convergence_gate_row(Lc)
    check("報告の収束の欄 (旧方式、TOOL: check_convergence.py の形) → 今までと同じ見出しと値 (OVERALL の行)",
          row == _old_gate_row(Lc) and row[0] == "収束 (check_convergence)" and row[1].startswith("OVERALL"), str(row))
    Lo = clone(L, "rep_legacy_old_format")
    txt = (Lc / "CONVERGENCE_VERDICT.txt").read_text().split("\n", 1)[1]      # TOOL の行の無い 2026-10-16 より前の形
    (Lo / "CONVERGENCE_VERDICT.txt").write_text(txt)
    check("報告の収束の欄 (TOOL の行の無い旧い CONVERGENCE_VERDICT.txt) → 変更前と同じ", NR.convergence_gate_row(Lo) == _old_gate_row(Lo))
    check("報告の収束の欄: CONVERGENCE_VERDICT.txt が無い → 行を出さない (今どおり)", NR.convergence_gate_row(TMP / "euler") is None)

    # --- forge の直接実行を弾くフックの案内 ---
    p = subprocess.run([sys.executable, str(TOOLS / "hook_forge_guard.py")], input=json.dumps({"tool_input": {"command": "build/forge"}}),
                       capture_output=True, text=True, timeout=60)
    r = json.loads(p.stdout)["hookSpecificOutput"]["permissionDecisionReason"]
    check("hook_forge_guard: 直接の forge の実行を弾き、案内に道具の選択 (check_settled.py・check_convergence.py・--select) を書く",
          json.loads(p.stdout)["hookSpecificOutput"]["permissionDecision"] == "deny" and "check_settled.py" in r and "--select" in r, r[:200])


def main() -> int:
    print(f"出力: {TMP}", flush=True)
    eu, pin = make_euler()
    A = test_gate(pin, eu)
    test_gate_v3(pin, A)
    test_migration(pin, A)
    rc_main, rc_legacy = test_run_case(A)
    test_hook(rc_main, rc_legacy)
    test_readers()
    print(f"\n{NCHK - FAIL}/{NCHK} ok — {'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (一時ディレクトリ {TMP})", flush=True)
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
