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
  (f) result 段レビュー 3 回目・4 回目 M1 (§5.1 #16): 停止の時の再評価は専用ファイル `stop_rule_v2/SETTLE_VERDICT_stop_<step>.json`
      (`stop_rule.write_stop_verdict`; 合成の run も同じ関数で書く。`legacy=True` は専用ファイルの無い古い形)。通常の CLI
      (`check_settled.py`、相対のパス・既定で SETTLE_VERDICT.json を上書き) で判定し直した後も、成果物の収束 (相対・絶対のパス) が
      PASS、受け渡しの照合 (`_v2_binding`) が通り、関門が進む — 新しい形 (専用ファイル) と古い形 (CLI が上書きの直前に停止の時の
      原本を同じ名前の専用ファイルへ移す; STOP_RULE_STATE は変えない) の両方。古い形を別の終端で判定し直しても原本が残るので通る
      (4 回目 M1 の前は負例)。--no-write・別の --out は移行しない。負例: 入力 (流量の CSV) を変えて通常の CLI (両方の形)、古い形で
      移行の前に SETTLE_VERDICT.json が上書きされていた (原本なし)、移した専用ファイルの書き換え・名前の衝突 (CLI は上書きしない)、
      新しい形で専用ファイルを書き換えた → 未評価・照合の失敗・関門は上書きでも止める。CLI の --out で専用ファイルを指すと引数の誤り。
      検出力: 修正前 (PRE_M1_COMMIT) は古い形で CLI の後に「未評価」・照合の失敗、修正前 (PRE_M1B_COMMIT = 61fc3d59) は入力を変えて
      通常の CLI を回すと古い形の成果物・受け渡し、新しい形の受け渡しが通った。
  (g) check_settled.py が写して持つ停止規則 v2 の名前 (置き場・状態のファイル・専用ファイルの名前) が stop_rule.py と同じ。

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
        # 投入の関門 (run_case.sh → launch_gate.py、plan tooling-generator-defaults-and-launch-gate §5.1 #10) がメッシュの品質
        # (check_mesh_quality.py --bcond) を見るので、node 変換と同じ置き方の primal の接続と壁 (physID 3 = 上辺) の面も書く
        # (座標だけの h5 は検査できない入力として止まる; run_stop_rule_tests.py の write_mesh と同じ形)
        nid = lambda i, j: i * nj + j                                                   # noqa: E731
        quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
        wall = [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)]
        f["MESH/CONNE"] = np.array([x for n in range(ni * nj) for x in (1, n)], dtype=np.int64)
        f["VIZMESH/CONNE"] = np.array([x for q in quads for x in [5] + q], dtype=np.int64)
        f["VIZMESH"].attrs["nVizCells"] = np.int64(len(quads))
        b = f.create_group("BCONDS/3")
        b["vizBfaceNodes"] = np.array([n for e in wall for n in e], dtype=np.int64)
        b["vizBfaceSizes"] = np.full(len(wall), 2, dtype=np.int64)
        b.attrs["bcondKind"] = "wall_isothermal"
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


# mesh・bcondConfig.yaml は投入の関門 (run_case.sh → launch_gate.py) が読む (meshFileName が無い・bcond が無い run は検査できない入力として止まる)
CFG_TEXT = ('mesh: {discretization: "node", isAxisymmetric: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}\n'
            "time: {last: {nStepOuter: %d}, outStepInterval: 2000}\noutput: {boundaryMassflux: 1, nodeProbeFile: settle_nodes.txt}\n")
BCOND_TEXT = ("inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 3500000.0, Tt: 1250.0}}\n"
              "outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {Ps: 2770.0}}\n"
              "wall:   {physID: 3, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Tw: 370.0}}\n")


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
          residual="judge", legacy=False) -> Path:
    """新方式の run。REACHED: 上限 100,000・要求の格子点 26,000・協調停止 26,001 (STOPPED_AT)。CENSORED: 上限 = 停止の step 26,000
    (STOPPED_AT なし)。delta_from: 同じ場の別の run の δ_E の抽出を写す (場は同じなので抽出も同じ; 試験の時間を縮める)。
    legacy: 停止の時の再評価の専用ファイルの無い古い形 (2026-10-16 の V4・V6 の形; `finish_v2` の legacy)。"""
    rd = TMP / name
    rd.mkdir()
    k, cap = steps_of(kind)
    X, R, Q = synth_ns()
    write_field(rd, X, R, Q, steps=(0, k))
    info = ns_info(pin)
    (rd / "prepare_info.json").write_text(json.dumps(info))
    (rd / "solverConfig.yaml").write_text(CFG_TEXT % cap)
    (rd / "bcondConfig.yaml").write_text(BCOND_TEXT)
    return finish_v2(rd, kind=kind, resid_fn=resid_fn, mdot_fn=mdot_fn, skip_nodes=skip_nodes,
                     delta=("from", delta_from) if delta_from is not None else ("extract", info), residual=residual, legacy=legacy)


def steps_of(kind):
    """(停止の step, 上限): REACHED は協調停止 26,001・上限 100,000、CENSORED は上限 = 停止の step 26,000。"""
    return (26001, 100000) if kind == "REACHED" else (26000, 26000)


def finish_v2(rd: Path, *, kind="REACHED", resid_fn=None, mdot_fn=None, skip_nodes=(), delta=None, residual="judge", legacy=False) -> Path:
    """run (場・solverConfig.yaml は置いてある) に停止規則 v2 の記録一式を書く: 基準・起動の記録・本段の開始の記録・forge の出力の形の
    データ・STOPPED_AT・停止の step の再評価 (停止の時の専用ファイル stop_rule_v2/SETTLE_VERDICT_stop_<k>.json と SETTLE_VERDICT.json;
    `stop_rule.write_stop_verdict`、§5.1 #16)・stop_rule.json・STOP_RULE_STATE.json。
    delta: ("extract", prepare_info) = 本物の抽出、("from", run) = 同じ場の別の run の抽出を写す、None = δ_E の記録なし (関門には使えない)。
    legacy: 専用ファイルの無い古い形 (§5.1 #16 より前の見張りと同じ: SETTLE_VERDICT.json だけに書き、reeval に stop_file が無い)。"""
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
    # 停止の step の再評価 (watch_main_stage_v2 と同じ: 新しい RunData・end_step = 停止の step → 専用ファイルと SETTLE_VERDICT.json)
    lc = ST.load_criteria(rd)
    data = ST.RunData(rd, lc["criteria"])
    r = ST.judge(rd, k, data=data, loaded=lc)
    if legacy:
        pv = ST.write_verdict(rd, r, data)
        reeval = {**SR.judge_summary(r), "end_source": r.get("end_source"), "file": ST.VERDICT_FILE, "sha256": SR._sha256(pv)}
    else:
        reeval = SR.write_stop_verdict(rd, k, r, data)
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
    # 収束の判定の記録は停止の時の専用ファイル (§5.1 #16、result 段レビュー 3 回目 M1; 以前は SETTLE_VERDICT.json — 通常の CLI が
    # 上書きする最新の判定なので、停止の時の照合に使わなくなった)
    stop_rel = SR.stop_verdict_relpath(k)
    check(f"A: 抽出の要約に落ち着きの印 (settled true・δ_E の作り方) を写す、収束の判定の記録は停止の時の専用ファイル {stop_rel} の残差",
          summ["settled"] is True and summ["delta_e_source"]["kind"] == "stop_rule_v2_stop_snapshot" and g["convergence"]["kind"] == "PASS"
          and g["convergence"]["tool"].startswith(stop_rel + " ") and c["binding"]["items"]["reeval"]["form"] == "stop_file"
          and stop_rel in (A / DL.HANDOFF_CONVERGENCE).read_text() and "OVERALL: SETTLED" in (A / DL.HANDOFF_CONVERGENCE).read_text(),
          g["convergence"]["tool"])

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
        # 停止の時の専用ファイルを終端 26,000 の判定で書き換えた (停止の step の再評価でない)。§5.1 #16 より前はこの負例を
        # SETTLE_VERDICT.json の書き換えで作っていたが、SETTLE_VERDICT.json は通常の CLI が上書きする最新の判定になり、停止の時の照合に
        # 使わない (書き換えても関門は進む — 正例は (f) の試験) ので、専用ファイルを書き換える
        lc = ST.load_criteria(rd)
        d = ST.RunData(rd, lc["criteria"])
        ST.write_verdict(rd, ST.judge(rd, 26000, data=d, loaded=lc), d, rd / SR.stop_verdict_relpath(26001))

    def mut_rejudged_bound(rd):
        # 同じ書き換えで、状態の記録 (要約・専用ファイルの sha256) も書き換えに合わせた: 停止の step との終端の照合だけで止まる
        lc = ST.load_criteria(rd)
        d = ST.RunData(rd, lc["criteria"])
        r = ST.judge(rd, 26000, data=d, loaded=lc)
        p_ = ST.write_verdict(rd, r, d, rd / SR.stop_verdict_relpath(26001))
        rewrite_records(rd, lambda rec: rec["reeval"].update(**SR.judge_summary(r), stop_file_sha256=SR._sha256(p_)))

    def mut_extract_settings(rd):
        rewrite_records(rd, lambda rec: rec["delta_E"]["extract"].update(band_select=repr("adaptive")))

    def mut_launch_other(rd):
        write_launch(rd, cfg_text="time: {last: {nStepOuter: 3000}}\n")    # 本段の後に別の設定で forge を起動した

    cases = [("B: MAIN_STAGE_STARTED.main_config.fnv だけ不一致 (判別の A/B)", mut_notice_fnv, "照合"),
             ("基準 settle_criteria.yaml の sha256 の不一致", mut_criteria, "照合"),
             ("STOPPED_AT の step の不一致", mut_stopped_at, "照合"),
             ("δ_E の status が ok でない (状態・見張りの記録とも)", mut_delta_status, "δ_E"),
             ("起動の記録 forge_launches.jsonl が無い (照合できない)", mut_no_launch, "照合"),
             ("停止の時の専用ファイルが停止の step の再評価でない (終端 26,000 で書き直された; sha256 が状態の記録と違う)", mut_rejudged, "照合"),
             ("停止の時の専用ファイルが停止の step の再評価でない (終端 26,000; 状態の記録も合わせた → 終端の照合)", mut_rejudged_bound,
              "再評価の終端"),
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
    """Stop フックは 2026-10-18 の forge-gate の取り込みから RUN_VERDICT.json (judge_run.py) の有無・鮮度・確定を見る
    (plan tooling-convergence-and-wall-resolution-gates §4.6-5; 以前の CONVERGENCE_VERDICT.txt の有無と時刻の検査は置き換わった)。
    ここでは run_case.sh が書いた基準のある run・基準の無い run の判定でフックが止めないこと、判定の無い run・判定の後に入力が
    変わった run で止めることを見る (FAIL・判定不能のまま終えるのは許す設計; 詳しい場合分けは tests/unit/test_judge_run.py)。"""
    print("== (d) Stop フック (hook_convergence_gate.py、RUN_VERDICT.json)", flush=True)
    root = TMP / "hookroot"
    case = root / "case" / "99.hooktest"
    case.mkdir(parents=True)

    def mk(name, src: Path | None):
        d = case / name
        if src is not None:
            shutil.copytree(src, d, symlinks=True)       # copy2: 判定の入力の識別 (sha256・h5 の mtime) を保つ
        else:
            d.mkdir()
            (d / "residual_history.csv").write_text("step,inner,phase,rms_ro\n0,-1,outer_end,1e-4\n")
        old = time.time() - 600                          # 実行中 (更新から 120 s 未満) として飛ばされないように
        os.utime(d / "residual_history.csv", (old, old))
        return d

    def hook():
        env = {**os.environ, "FORGE_GATE_ROOT": str(root), "HOME": str(TMP / "hookhome")}
        p = subprocess.run([sys.executable, str(TOOLS / "hook_convergence_gate.py")], input="{}", capture_output=True, text=True,
                           timeout=120, env=env)
        return json.loads(p.stdout.strip() or "{}")

    mk("run_0001_settled", rc_main)
    mk("run_0002_convergence", rc_legacy)
    out = hook()
    check("run_case.sh が判定を書いた run (基準あり・基準なし) では block しない (判定の総合が FAIL・判定不能でも終えてよい)",
          out == {} and (rc_main / "RUN_VERDICT.json").is_file() and (rc_legacy / "RUN_VERDICT.json").is_file(), str(out)[:300])
    mk("run_0003_none", None)
    out = hook()
    check("RUN_VERDICT.json が無い run があれば block (その run だけを挙げる)、案内は judge_run.py",
          out.get("decision") == "block" and "run_0003_none" in out.get("reason", "") and "run_0001" not in out.get("reason", "")
          and "judge_run.py" in out.get("reason", ""), str(out)[:300])
    shutil.rmtree(case / "run_0003_none")
    d = mk("run_0004_stale", rc_main)
    with open(d / "residual_history.csv", "a") as f:
        f.write("# 判定の後に書き足した\n")
    old = time.time() - 600
    os.utime(d / "residual_history.csv", (old, old))
    out = hook()
    check("判定の後に残差の履歴が変わった run → block (STALE)", out.get("decision") == "block" and "run_0004_stale" in out.get("reason", ""),
          str(out)[:300])


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
    """停止規則の停止の時の再評価の記録 (専用ファイル stop_rule_v2/SETTLE_VERDICT_stop_<k>.json; 古い形の run は SETTLE_VERDICT.json) を
    mutate(v) で書き換え、STOP_RULE_STATE の reeval (digest・sha256・stop_file_sha256・verdict) を合わせる (別の機械・別の所在で書いた
    記録の形; 記録どうしは結び付いたまま)。専用ファイルのある run は、停止規則と同じく同じバイトを SETTLE_VERDICT.json にも書く。"""
    rel = (SR.read_state(rd)["reeval"] or {}).get("stop_file")
    p = rd / (rel or ST.VERDICT_FILE)
    v = json.loads(p.read_text())
    mutate(v)
    txt = json.dumps(v, ensure_ascii=False, indent=1) + "\n"
    p.write_text(txt)
    if rel:
        (rd / ST.VERDICT_FILE).write_text(txt)
    upd = {"digest": SR.quantities_digest(v), "sha256": SR._sha256(rd / ST.VERDICT_FILE), "verdict": v["verdict"]}
    if rel:
        upd["stop_file_sha256"] = SR._sha256(p)
    rewrite_records(rd, lambda rec: rec["reeval"].update(**upd))


def test_reeval_match(V: Path, VL: Path):
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
    # 停止の時の記録が無い: 新しい形は専用ファイル、古い形は SETTLE_VERDICT.json (§5.1 #16 より前は SETTLE_VERDICT.json を消す負例だけ
    # だった — 新しい形では SETTLE_VERDICT.json は最新の判定で、消しても停止の時の照合は専用ファイルで成り立つ)
    Vv = clone(V, "wt_v2_noverdict")
    stop_rel = SR.read_state(Vv)["reeval"]["stop_file"]
    (Vv / stop_rel).unlink()
    c = conv(Vv)
    check("停止の時の再評価の専用ファイルが無い → 例外にせず理由つきの未評価 (SETTLE_VERDICT.json へ落とさない)", c["result"] == "未評価"
          and stop_rel in (c.get("reason") or "") and (Vv / ST.VERDICT_FILE).is_file(), str(c.get("reason"))[:300])
    Vv2 = clone(V, "wt_v2_noverdict_latest")
    (Vv2 / ST.VERDICT_FILE).unlink()
    c = conv(Vv2)
    check("最新の判定 SETTLE_VERDICT.json が無くても停止の時の専用ファイルと照合できる → 収束 PASS", c["result"] == "PASS"
          and c["stop_rule_reeval"]["compare"]["form"] == "stop_file", str(c.get("reason"))[:300])
    Vl = clone(VL, "wt_v2_legacy_noverdict")
    (Vl / ST.VERDICT_FILE).unlink()
    c = conv(Vl)
    check("古い形の run (専用ファイルなし) で再評価の記録 SETTLE_VERDICT.json が無い → 例外にせず理由つきの未評価", c["result"] == "未評価"
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
    VL = TMP / "wt_v2_legacy"                       # 停止の時の専用ファイルの無い古い形 (§5.1 #16 より前の見張りの形)
    RWT.make_synth(VL, steps=wsteps, tw=1e-5, residual=False)
    finish_v2(VL, kind="REACHED", delta=None, legacy=True)
    test_reeval_match(V, VL)
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



# ====================================================================================== (f) 通常の CLI で判定し直した後の照合 (§5.1 #16)
PRE_M1_COMMIT = "8cf2a284"          # result 段レビュー 3 回目 M1 (§5.1 #16) の修正の前 (検出力の確かめに使う)
PRE_M1B_COMMIT = "61fc3d59"         # result 段レビュー 4 回目 M1 (§5.1 #16 の続き) の修正の前: 古い形で原本が失われたとき要約の一致で
                                    # 代わりの記録 (form overwritten) にし、今の入力から作り直した判定と比べた (検出力の確かめに使う)


def _cli_rejudge(rd: Path, *args, mod=None):
    """通常の CLI (check_settled.py; 既定で SETTLE_VERDICT.json に書く) を run の親から相対のパスで回す。mod: 別の版の check_settled
    (検出力の確かめ; 既定は今の作業ツリー)。戻り (終了コード, 最後の行, 出力の全文)。"""
    import io
    from contextlib import redirect_stdout
    import check_settled
    mod = mod or check_settled
    cwd = os.getcwd()
    buf = io.StringIO()
    try:
        os.chdir(rd.parent)
        with redirect_stdout(buf):
            try:
                rc = mod.main([rd.name, *args])
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    out = buf.getvalue()
    return rc, (out.strip().splitlines() or [""])[-1], out


def _with_stop_rule(mod, fn):
    """`from ..evaluate import stop_rule` を mod (修正前の stop_rule) に差し替えて fn() を呼ぶ (検出力の確かめ; 呼んだ後に戻す)。"""
    import forge_design.evaluate as EV
    key = "forge_design.evaluate.stop_rule"
    saved_m, saved_a = sys.modules[key], EV.stop_rule
    sys.modules[key], EV.stop_rule = mod, mod
    try:
        return fn()
    finally:
        sys.modules[key], EV.stop_rule = saved_m, saved_a


def _bump_massflux(rd: Path, dm: float = 1.0):
    """流量の CSV の流入・流出の絶対値を各 dm [kg/s] 増やす (result 段レビュー 4 回目 M1 の再現と同じ入力の書き換え; 収支 r は変わらず
    SETTLED のまま)。写しの実ファイルを消してから書く (symlink を通して元を書き換えない)。"""
    p = rd / "boundary_massflux.csv"
    lines = p.read_text().splitlines()
    hi = next(i for i, ln in enumerate(lines) if ln.startswith("step,"))
    cols = lines[hi].split(",")
    sign = {c: (-1.0 if "inlet" in c else 1.0) for c in cols if c.startswith(("mdot_inlet", "raw_inlet", "mdot_outlet", "raw_outlet"))}
    out = lines[:hi + 1]
    for ln in lines[hi + 1:]:
        if not ln or ln.startswith("#"):
            out.append(ln)
            continue
        v = ln.split(",")
        for c, sg in sign.items():
            j = cols.index(c)
            v[j] = repr(float(v[j]) + sg * dm)
        out.append(",".join(v))
    p.unlink()
    p.write_text("\n".join(out) + "\n")


def test_cli_rejudge(pin, A: Path):
    """result 段レビュー 3 回目・4 回目 M1 (plan §5.1 #16): 通常の CLI (`check_settled.py`、既定で SETTLE_VERDICT.json を上書き) で判定し
    直した後も、成果物の収束 (`wall_tables.thermal_use`、相対・絶対のパス) と受け渡しの照合 (`_v2_binding`・関門) が停止の時の記録と
    照合できること。新しい形 (停止の時の専用ファイルあり) と古い形 (専用ファイルなし、§5.1 #16 より前の見張りの形; 通常の CLI が上書きの
    直前に停止の時の原本を専用ファイルへ移す) の両方。負例: 入力を変えて通常の CLI (両方の形)、古い形で移行の前に SETTLE_VERDICT.json が
    上書きされていた (原本なし)、移した専用ファイルの書き換え・名前の衝突、新しい形で専用ファイルを書き換えた、CLI の --out で専用ファイル
    を指す。検出力: 修正前 (PRE_M1_COMMIT) は古い形で CLI の後に照合が失敗、修正前 (PRE_M1B_COMMIT) は入力を変えて CLI を回すと通った。"""
    print("== (f) 通常の CLI で判定し直す → 成果物の照合 → 受け渡しの照合 (result 段レビュー 3 回目・4 回目 M1、§5.1 #16)", flush=True)
    sys.path.insert(0, str(DESIGN / "tests"))
    import run_wall_tables_tests as RWT
    from forge_design.export import wall_tables as WT
    import check_settled as CS
    wsteps = (18000, 20000, 22000, 24000, 26000, 26001)
    k = 26001
    stop_rel = SR.stop_verdict_relpath(k)

    # (g) check_settled.py が写して持つ停止規則 v2 の名前 (design のパッケージを読まないため) が stop_rule.py と同じ
    check("(g) check_settled.py の停止規則 v2 の名前 (置き場・状態のファイル・専用ファイルの名前) が stop_rule.py と同じ",
          (CS.STOP_RULE_V2_DIR, CS.STOP_RULE_STATE, CS.STOP_VERDICT_NAME) == (SR.V2_DIR, SR.STATE, SR.V2_STOP_VERDICT)
          and SR.legacy_stop_verdict_relpath({"end_step": k}) == stop_rel,
          f"{(CS.STOP_RULE_V2_DIR, CS.STOP_RULE_STATE, CS.STOP_VERDICT_NAME)}")

    def wall_run(name, legacy):
        """壁のダンプのある run (成果物の収束) に、停止の場 res_26001.h5 (照合は sha256 だけを見る) と停止規則 v2 の記録一式。"""
        rd = TMP / name
        RWT.make_synth(rd, steps=wsteps, tw=1e-5, residual=False)
        with h5py.File(rd / f"res_{k}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(4)
        return finish_v2(rd, kind="REACHED", delta=None, legacy=legacy)

    def unpublish(rd):
        for f in ("delta_r_next.csv", DL.HANDOFF_GATE, "fixed_point.json", "delta_r_extract_summary.json", "delta_r_equiv.csv",
                  "delta_r_equiv_diag.json", "delta_r_edge.csv"):
            (rd / f).unlink(missing_ok=True)
        return rd

    def thermal(rd, mod=WT):
        cwd = os.getcwd()
        try:
            os.chdir(rd.parent)
            c_rel = mod.thermal_use(Path(rd.name))["convergence"]
        finally:
            os.chdir(cwd)
        return c_rel, mod.thermal_use(rd.resolve())["convergence"]

    def conv_ok(c, form):
        return (c["result"] == "PASS" and c["stop_rule_reeval"].get("matches") is True
                and (c["stop_rule_reeval"].get("compare") or {}).get("form") == form)

    def chain_ok(rd, form):
        """成果物の収束 (相対・絶対) が PASS・受け渡しの照合が通る (記録の形 form)。戻り (合否, 詳細)。"""
        c_rel, c_abs = thermal(rd)
        b = DL._v2_binding(rd, SR.read_state(rd), ST)
        ok = (conv_ok(c_rel, form) and conv_ok(c_abs, form) and b["pass"] and b["items"]["reeval"]["form"] == form
              and b["items"]["reeval"]["exact"] is True and "rejudge_match" not in b["items"]["reeval"])
        return ok, f"{c_rel.get('reason')} / {c_abs.get('reason')} / {b['problems'][:2]}"

    W_new, W_old = wall_run("cli_wt_new", False), wall_run("cli_wt_legacy", True)
    A_old = build("cli_A_legacy", pin, legacy=True, delta_from=A)
    for label, Wr, Ar, form in (("新しい形 (停止の時の専用ファイルあり)", W_new, A, "stop_file"),
                                ("古い形 (専用ファイルなし; V4・V6 の形)", W_old, A_old, "migrated")):
        # 成果物と受け渡しの照合 (同じ run)
        rd = clone(Wr, f"cli_chain_{form}")
        st = SR.read_state(rd)
        ve = st["reeval"]
        state_bytes = (rd / SR.STATE).read_bytes()
        stop_bytes = (rd / stop_rel).read_bytes() if form == "stop_file" else None
        c_rel, c_abs = thermal(rd)
        b = DL._v2_binding(rd, st, ST)
        check(f"(f) {label}: CLI の前 — 成果物の収束 PASS・受け渡しの照合が通る (記録の形 "
              + ("stop_file" if form == "stop_file" else "settle_verdict = 移行していない原本") + ")",
              c_rel["result"] == c_abs["result"] == "PASS" and b["pass"]
              and b["items"]["reeval"]["form"] == ("stop_file" if form == "stop_file" else "settle_verdict")
              and (form == "stop_file" or not (rd / stop_rel).exists()), f"{c_rel.get('reason')} / {b['problems'][:2]}")
        orig = (rd / ST.VERDICT_FILE).read_bytes()
        rc, last, out = _cli_rejudge(rd)
        if form == "stop_file":
            check(f"(f) {label}: 通常の CLI (相対のパス・既定で書く) → SETTLED・終了コード 0、SETTLE_VERDICT.json は停止の時と違うバイトに"
                  "上書き、停止の時の専用ファイルは変わらない、移行の文を出さない",
                  rc == 0 and last.startswith("VERDICT: SETTLED") and SR._sha256(rd / ST.VERDICT_FILE) != ve["sha256"]
                  and (rd / stop_rel).read_bytes() == stop_bytes and "移行:" not in out, f"rc {rc} {last}")
        else:
            check(f"(f) {label}: 通常の CLI → SETTLED・終了コード 0、上書きの前に停止の時の原本を {stop_rel} へ移した (バイト同一・sha256 = "
                  "状態の記録)、SETTLE_VERDICT.json は上書き、STOP_RULE_STATE.json は変えない",
                  rc == 0 and last.startswith("VERDICT: SETTLED") and "移行:" in out and (rd / stop_rel).read_bytes() == orig
                  and SR._sha256(rd / stop_rel) == ve["sha256"] and SR._sha256(rd / ST.VERDICT_FILE) != ve["sha256"]
                  and (rd / SR.STATE).read_bytes() == state_bytes, f"rc {rc} {last} / {out[-300:]}")
        ok, det = chain_ok(rd, form)
        check(f"(f) {label}: CLI の後の成果物の収束 (相対・絶対のパス) PASS・受け渡しの照合 (_v2_binding) が通る (記録の形 {form}、"
              "停止の時の原本そのもの = exact、今の判定との代わりの照合はしない)", ok and DL._v2_binding(rd, SR.read_state(rd), ST)["verdict_file"]
              == stop_rel, det)
        rc2, last2, out2 = _cli_rejudge(rd)
        ok, det = chain_ok(rd, form)
        check(f"(f) {label}: もう一度通常の CLI → 移した専用ファイルは変わらない (移行は 1 回だけ)・照合は通るまま",
              rc2 == 0 and "移行:" not in out2 and (form == "stop_file" or (rd / stop_rel).read_bytes() == orig) and ok, det)
        # 受け渡しの関門 (δ_E のある run)
        ag = unpublish(clone(Ar, f"cli_gate_{form}"))
        rc, last, _ = _cli_rejudge(ag)
        g = gate(ag, pin)
        check(f"(f) {label}: CLI の後の受け渡しの関門 → 進む (上書きなし・settled true)、run の直下に delta_r_next.csv を公開",
              rc == 0 and g["proceed"] and g["settled"] is True and g["override"]["applied"] == []
              and g["checks"]["binding"]["items"]["reeval"]["form"] == form and (ag / "delta_r_next.csv").is_file(),
              "; ".join(g["blocked_by"])[:300])

    # 古い形で別の終端 (--end-step 26000) の通常の CLI: 4 回目 M1 の修正の前は「要約が停止の時と違う」で未評価になる負例だった。
    # 修正後は上書きの前に原本を移すので、別の終端で判定し直しても停止の時の記録は残り、照合は通る (期待値を負例 → 正例に変えた)
    rd = clone(W_old, "cli_old_other_end")
    rc, last, out = _cli_rejudge(rd, "--end-step", "26000")
    ok, det = chain_ok(rd, "migrated")
    check("(f) 古い形: 通常の CLI を別の終端 (--end-step 26000) で回す → 上書きの前に原本を移し、成果物の収束 PASS・受け渡しの照合が通る "
          "(4 回目 M1 の前は要約の不一致で未評価の負例; 停止の時の記録が残るので正例に変えた)",
          "移行:" in out and json.loads((rd / ST.VERDICT_FILE).read_text())["end_step"] == 26000 and ok, det)
    # --no-write・別の --out は SETTLE_VERDICT.json を上書きしないので移行もしない (原本は SETTLE_VERDICT.json のまま)
    rd = clone(W_old, "cli_old_nowrite")
    _cli_rejudge(rd, "--no-write")
    _, _, out_o = _cli_rejudge(rd, "--out", str(Path(rd.name) / "other_verdict.json"))   # _cli_rejudge は run の親で回す
    c_rel, _ = thermal(rd)
    check("(f) 古い形: --no-write・別の --out は移行しない (専用ファイルなし、SETTLE_VERDICT.json は原本のまま → 記録の形 settle_verdict で PASS)",
          not (rd / stop_rel).exists() and "移行:" not in out_o and (rd / "other_verdict.json").is_file() and conv_ok(c_rel, "settle_verdict"),
          str(c_rel.get("reason"))[:300])

    # --- 負例 ---
    try:
        cs_old = _load_ref("solver_density_cuda/tools/check_settled.py", "_check_settled_pre_m1b", commit=PRE_M1B_COMMIT)
    except Exception as e:  # noqa: BLE001
        cs_old = None
        check(f"(f) 修正前 ({PRE_M1B_COMMIT}) の check_settled.py を git から読めた (読めなければ原本なしの負例・検出力は判定不能 — 合格にしない)",
              False, str(e))

    def overwrite_before_migration(rd):
        # 移行の前の版の通常の CLI (PRE_M1B_COMMIT; 移行しない) で SETTLE_VERDICT.json を上書きした古い形の run (原本は残っていない)。
        # その後の今の CLI は原本でないので移行しない
        _cli_rejudge(rd, mod=cs_old)
        _cli_rejudge(rd)

    def conflict(rd):
        # 古い形で、移行の前に同じ名前の専用ファイルに違う中身がある → CLI は原本を失わないよう SETTLE_VERDICT.json を上書きしない
        (rd / stop_rel).parent.mkdir(parents=True, exist_ok=True)
        (rd / stop_rel).write_text("{}\n")
        orig = (rd / ST.VERDICT_FILE).read_bytes()
        rc_, _, out_ = _cli_rejudge(rd)
        if not (rc_ == 0 and "書かない" in out_ and (rd / ST.VERDICT_FILE).read_bytes() == orig):
            raise AssertionError(f"衝突で SETTLE_VERDICT.json を上書きした・終了コード {rc_}: {out_[-200:]}")

    negs = (("古い形: 入力 (流量の CSV の流入・流出を各 +1 kg/s) を変えて通常の CLI (レビューの再現; 原本は移した)", "legacy",
             lambda rd: (_bump_massflux(rd), _cli_rejudge(rd)), "中身 (sha256) が記録と違う"),
            ("新しい形: 入力 (流量の CSV の流入・流出を各 +1 kg/s) を変えて通常の CLI", "new",
             lambda rd: (_bump_massflux(rd), _cli_rejudge(rd)), "中身 (sha256) が記録と違う"),
            ("古い形: 移行の前に SETTLE_VERDICT.json が上書きされていた (原本なし) → 今の CLI", "legacy",
             overwrite_before_migration, "原本が失われた"),
            ("古い形: 移した専用ファイルを停止の後に書き換えた", "legacy",
             lambda rd: (_cli_rejudge(rd), (rd / stop_rel).write_text((rd / stop_rel).read_text().replace('"verdict": "SETTLED"',
                                                                                                         '"verdict": "SETTLED" ', 1))),
             "停止の時の原本でない"),
            ("古い形: 移行の前に同じ名前の専用ファイルに違う中身がある (CLI は SETTLE_VERDICT.json を上書きしない)", "legacy",
             conflict, "停止の時の原本でない"),
            ("新しい形: 停止の時の専用ファイルを停止の後に書き換えた (sha256 が状態の記録と違う)", "new",
             lambda rd: (_cli_rejudge(rd), (rd / stop_rel).write_text((rd / stop_rel).read_text().replace('"verdict": "SETTLED"',
                                                                                                         '"verdict": "SETTLED" ', 1))),
             "停止の後に書き換えられた"))
    for i, (lab, form, mut, key) in enumerate(negs):
        if cs_old is None and mut is overwrite_before_migration:
            continue
        rd = clone(W_old if form == "legacy" else W_new, f"cli_neg_wt_{i}_{form}")
        ag = unpublish(clone(A_old if form == "legacy" else A, f"cli_neg_gate_{i}_{form}"))
        errs = []
        for r_ in (rd, ag):
            try:
                mut(r_)
            except AssertionError as e:
                errs.append(str(e))
        c_rel, _ = thermal(rd)
        b = DL._v2_binding(rd, SR.read_state(rd), ST)
        g = gate(ag, pin, **ALL_OV)
        check(f"(f) 負例 {lab} → 成果物の収束は未評価、受け渡しの照合は失敗、関門は上書き 3 つでも止める (理由に「{key}」)",
              not errs and c_rel["result"] == "未評価" and key in (c_rel.get("reason") or "") and not b["pass"]
              and any(key in p_ for p_ in b["problems"]) and not g["proceed"] and blocked_by(g, "照合") and not_published(ag),
              f"{errs[:1]} / {c_rel.get('reason')} / {b['problems'][:2]} / {g['blocked_by'][:1]}")
    rd = clone(W_new, "cli_out_guard")
    before = (rd / stop_rel).read_bytes()
    rc, last, _ = _cli_rejudge(rd, "--out", str(Path(rd.name) / stop_rel))
    check("(f) CLI の --out で停止の時の専用ファイルを指す → 引数の誤り (INDETERMINATE・終了コード 3)、専用ファイルは変わらない",
          rc == 3 and last == "VERDICT: INDETERMINATE (引数の誤り)" and (rd / stop_rel).read_bytes() == before, f"rc {rc} {last}")

    # --- 検出力 1: 修正前 (PRE_M1_COMMIT) は古い形で CLI の後に照合が失敗した ---
    try:
        sr_old = _load_ref("design/forge_design/evaluate/stop_rule.py", "forge_design.evaluate._stop_rule_pre_m1", commit=PRE_M1_COMMIT)
        dl_old = _load_ref("design/forge_design/feedback/deltastar_loop.py", "forge_design.feedback._deltastar_loop_pre_m1", commit=PRE_M1_COMMIT)
        wt_old = _load_ref("design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_m1", commit=PRE_M1_COMMIT)
    except Exception as e:  # noqa: BLE001
        check(f"(f) 修正前 ({PRE_M1_COMMIT}) のモジュールを git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
    else:
        rd = clone(W_old, "cli_detect")
        pre = _with_stop_rule(sr_old, lambda: (thermal(rd, wt_old)[0], dl_old._v2_binding(rd, sr_old.read_state(rd), ST)))
        _cli_rejudge(rd)
        post = _with_stop_rule(sr_old, lambda: (thermal(rd, wt_old)[0], dl_old._v2_binding(rd, sr_old.read_state(rd), ST)))
        check("(f) 検出力: 修正前 (3 回目 M1) は古い形の run で、CLI の前は成果物の収束 PASS・受け渡しの照合が通り、通常の CLI の後は成果物の"
              "収束「未評価」・受け渡しの照合が失敗 (レビューの再現)",
              pre[0]["result"] == "PASS" and pre[1]["pass"] and post[0]["result"] == "未評価" and not post[1]["pass"],
              f"前 {pre[0]['result']}・{pre[1]['pass']} / 後 {post[0]['result']}・{post[1]['problems'][:1]}")

    # --- 検出力 2: 修正前 (PRE_M1B_COMMIT) は「入力を変える → 通常の CLI (その版) → 成果物・受け渡し」が通った ---
    try:
        sr_b = _load_ref("design/forge_design/evaluate/stop_rule.py", "forge_design.evaluate._stop_rule_pre_m1b", commit=PRE_M1B_COMMIT)
        dl_b = _load_ref("design/forge_design/feedback/deltastar_loop.py", "forge_design.feedback._deltastar_loop_pre_m1b", commit=PRE_M1B_COMMIT)
        wt_b = _load_ref("design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_m1b", commit=PRE_M1B_COMMIT)
        if cs_old is None:
            raise RuntimeError("check_settled.py の修正前を読めない")
    except Exception as e:  # noqa: BLE001
        check(f"(f) 修正前 ({PRE_M1B_COMMIT}) のモジュールを git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
        return
    res = {}
    for form, Wr in (("legacy", W_old), ("new", W_new)):
        rd = clone(Wr, f"cli_detect_b_{form}")
        _bump_massflux(rd)
        _cli_rejudge(rd, mod=cs_old)
        res[form] = _with_stop_rule(sr_b, lambda: (thermal(rd, wt_b)[0], dl_b._v2_binding(rd, sr_b.read_state(rd), ST)))
    lg, nw = res["legacy"], res["new"]
    check(f"(f) 検出力: 修正前 ({PRE_M1B_COMMIT}) は入力を変えて通常の CLI を回すと、古い形で成果物の収束 PASS (記録の形 overwritten)・"
          "受け渡しの照合が通り、新しい形でも受け渡しの照合が通った (今は両方とも拒否)",
          lg[0]["result"] == "PASS" and (lg[0]["stop_rule_reeval"].get("compare") or {}).get("form") == "overwritten" and lg[1]["pass"]
          and nw[1]["pass"],
          f"古い形 {lg[0]['result']}・{lg[1]['pass']} / 新しい形 {nw[0]['result']}・{nw[1]['pass']}")


def main() -> int:
    print(f"出力: {TMP}", flush=True)
    eu, pin = make_euler()
    A = test_gate(pin, eu)
    test_gate_v3(pin, A)
    test_migration(pin, A)
    rc_main, rc_legacy = test_run_case(A)
    test_hook(rc_main, rc_legacy)
    test_readers()
    test_cli_rejudge(pin, A)
    print(f"\n{NCHK - FAIL}/{NCHK} ok — {'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (一時ディレクトリ {TMP})", flush=True)
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
