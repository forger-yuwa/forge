#!/usr/bin/env python3
"""停止規則 v2 (`evaluate/stop_rule.py` の `watch_main_stage_v2`・`runner_axismach.run_staged_ns` の分岐) の試験
(plan tooling-settle-convergence-and-monitor §4.4・§4.8 の 1・§6 V4、§5.1 #5)。

偽の forge (`FORGE_BIN` を下の台本に差し替え、本物の `run_case.sh` 経由) は forge の出力の契約
(procedures/solver-settings.md「本段の監視用の出力と協調停止」・`solver_density_cuda/output/monitorOutputs.cpp`) を写す:
step i = 1..nStepOuter を回し、毎 step STOP_REQUEST を確かめる。出力の step (outStepInterval の倍数・最後・停止の step) に
res_nodes_i.h5 → 流量の行 i の順に書き (流量の行が完成の印)、その後に残差の行 i−1 (outer_end) を書く (forge と同じ順)。
全場 res_i.h5 は最後と停止の step。STOP_REQUEST を見た step で出力を書き、残差を flush して STOPPED_AT (forge と同じ JSON) を置いて
終了コード 0。起動時に STOP_REQUEST・STOPPED_AT・SNAPSHOT_REQUEST があれば終了コード 1。
試験を決定的にするため、偽の forge は判定できる格子点ごとに見張りの記録 (stop_rule.json の judgements) にその格子点が載るまで待つ
(lockstep; 見張りは判定 → 要求を置く → 記録の順に書くので、要求は格子点 s の直後の step s+1 で見える)。

  (a) 部品: 最後の完成した行・残差の行の完成・流量の最後の step (書きかけの最後の行を読まない)、オンラインの RunData が書きかけの行を
      読まず path を元のファイルに戻す、STOP_REQUEST の書き方、STOPPED_AT の検査、起動前の検査 (制御のファイルの残り・基準の誤用・
      未設定・出力の間隔・節点の表・出力の設定)。
  (b) 偽の forge で `run_staged_ns`:
      - REACHED (M1 の正例): 格子 2,000・窓 20,000・保持 38,000 で 42,000 から SETTLED → 80,000 で保持を満たして STOP_REQUEST →
        **非出力の step 80,001 で止まる** → STOPPED_AT の検査 → 80,001 で再評価 SETTLED → REACHED。δ_E の抽出 (偽物) は 80,001 で 1 回。
        オンラインの格子点ごとの判定が終了後の再生 (`settle.judge(run, s)`) と同じ (verdict と量の sha256)、止めた step の再評価
        (SETTLE_VERDICT.json) が `check_settled.py --end-step 80001` と同じ (終了コード 0)。
      - 保持の既定 (基準のファイルに hold が無い → 40,000): 82,000 で要求、82,001 で REACHED。
      - 要求の後に崩れる: (i) 停止の step の観測だけが跳ねる → 再評価 NOT_SETTLED → FAILED (reeval_not_settled)、
        (ii) 要求を見てから 4,000 step 遅れて止まり、その間に出口の M が変わる → 要求の後の格子点が NOT_SETTLED → FAILED
        (broken_after_request)。
      - forge が協調停止で終了コード 3・STOPPED_AT なし → FAILED (forge_exit・no_stopped_at)。STOPPED_AT の step がずれている → FAILED
        (stopped_at_invalid)。
      - 上限で止まる → CENSORED (保持の途中と記録、δ_E は上限の step で 1 回)。
      - 要求を無視する forge → stop_wait_s の後に kill → FAILED (timeout; REACHED にしない)。
      - forge は要求で止まったが run_case.sh の後処理 (check_convergence.py) が stop_wait_s を超える → 時間切れにせず REACHED
        (§5.1 #5b; forge の終了は run_case_stdout.log の `[run_case] forge exit=`)。後処理は PATH の先頭に置いた python3 の包み
        (check_convergence.py のときだけ眠ってから本物を exec) で遅らせる。
      - 残差が非有限 → DIVERGED の判定で STOP_REQUEST → FAILED (diverged)。基準に無い残差の列が値を持つ → 未設定の判定で止めて FAILED
        (criteria_unset)。
      - 起動前に STOP_REQUEST が残っている → forge を起動せず RuntimeError (ファイルは消さない)。基準の誤り・出力の設定の不足も起動しない。
      - settle_criteria.yaml の無い run は旧方式 (STOP_RULE_STATE.json の method v1; 旧方式の試験は run_stop_rule_tests.py)。
        stop_rule=False は v2 の run でも NOT_USED (method v2)。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_stop_rule_v2_tests.py [--quick]
所要: (a) 数秒、(b) 約 2〜3 分 (偽の forge 約 14 本)。--quick は (b) を省く。出力は一時ディレクトリ (全部合格なら消す)。
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
QUICK = "--quick" in sys.argv[1:]
TMP = Path(tempfile.mkdtemp(prefix="stop_rule_v2_tests_"))
FAKE_FORGE = TMP / "fake_forge_v2.py"
os.environ["FORGE_BIN"] = str(FAKE_FORGE)            # runner._ENV は import 時の環境を写すので import の前に
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import settle_presets as SP  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.evaluate import runner as RN  # noqa: E402

ST = SR._settle()
FAIL = 0
NI, NJ, SCALE = 6, 5, 0.1
DT, W = 2000, 20000
TEMPLATE_CFG = '''mesh: {discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 10000}
  deltaT: {control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0, dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 2000
  fieldOutStepInterval: 20000
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 1, limiter: 2}
turbulence: {model: "none"}
initial: "uniform_p101325_u10"
output: {level: 1, nodeProbeFile: nodes.txt, boundaryMassflux: 1}
'''
TEMPLATE_BC = ("inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 3500000.0, Tt: 1250.0}}\n"
               "outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {Ps: 2770.0}}\n"
               "wall:   {physID: 3, kind: slip, outputHDFflg: 1, ints: , floats: }\n")


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception, must=""):
    try:
        fn()
    except exc as e:  # noqa: BLE001
        return (must in str(e)), str(e)
    return False, "例外にならない"


# ====================================================================================== 合成の構造格子・基準
def coords():
    x = np.linspace(0.0, 1.0, NI) * SCALE * 10
    X = np.repeat(x[:, None], NJ, axis=1)
    R = (0.5 + 0.1 * np.linspace(0, 1, NI))[:, None] * np.linspace(0.0, 1.0, NJ)[None, :] * SCALE * 10
    return np.c_[X.ravel(), R.ravel(), np.zeros(NI * NJ)]


def criteria(hold: int | None = 38000, dt: int = DT, quantities: list | None = None) -> dict:
    q = quantities if quantities is not None else [
        {"kind": "residual", "file": "residual_history.csv", "required": ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe"],
         "inactive_ok": ["rms_roUz"],
         "thresholds": {"slope_dec": math.log10(1.05), "amp_dec": 0.1, "rise_dec": math.log10(1.10), "excursion_dec": 0.5}},
        {"kind": "mass_balance", "file": "boundary_massflux.csv", "inlet": [1], "outlet": [2],
         "thresholds": {"mean_abs": 5e-4, "swing": 5e-3, "mdot_change_rel": 1e-4}},
        {"kind": "exit_core_mach", "eta_range": [0.05, 0.7], "thresholds": {"tol": 1e-4}},
        {"kind": "exit_flow_angle", "eta_range": [0.05, 0.7], "thresholds": {"tol": 0.01}},
    ]
    c = {"schema": ST.SCHEMA, "source": {"made_by": "run_stop_rule_v2_tests.py"},
         "structured": {"ni": NI, "nj": NJ, "scale_m": SCALE, "mesh_file": "nozzle.h5"},
         "grid": {"out_step_interval": dt, "start_step": 0, "segment": None}, "window": {"W": W, "min_points": 11},
         "quantities": q, "excluded": [], "preconditions": []}
    if hold is not None:
        c["hold"] = hold
    return c


FAKE_SRC = r'''#!{py}
"""偽の forge (停止規則 v2 の試験用; docstring は run_stop_rule_v2_tests.py)。"""
import json, os, re, sys, time
import h5py, numpy as np
for f in ("STOP_REQUEST", "STOPPED_AT", "SNAPSHOT_REQUEST"):
    if os.path.exists(f):
        print(f"[output] 起動時の検査で停止: 作業ディレクトリに {{f}} が既にある", flush=True); sys.exit(1)
cfg = open("solverConfig.yaml").read()
n_outer = int(re.search(r"nStepOuter:\s*(\d+)", cfg).group(1))
dt = int(re.search(r"outStepInterval:\s*(\d+)", cfg).group(1))
sc = json.load(open("fake_scenario.json"))
coord = np.load("fake_coord.npy")
ids = [int(l.split("#")[0]) for l in open("nodes.txt") if l.split("#")[0].strip()]
cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe"] + (["rms_roK"] if sc.get("extra_col") else [])
hist = open("residual_history.csv", "w")
hist.write("step,inner,phase," + ",".join(cols) + "\n")
mf = open("boundary_massflux.csv", "w")
mf.write("# fake boundary_massflux.csv\nstep,mdot_inlet_1,raw_inlet_1,mdot_outlet_2,raw_outlet_2,mdot_wall_3,raw_wall_3\n"); mf.flush()
A = 340.0
def mach(step, stop):
    m = 5.5
    if sc.get("jump_at") is not None and step >= sc["jump_at"]:
        m += sc.get("jump", 0.01)
    if stop:
        m += sc.get("stop_jump", 0.0)
    return m
def write_field(step):
    with h5py.File(f"res_{{step}}.h5", "w") as h:
        h["VALUE/ro"] = np.ones(len(coord))
def write_outputs(step, stop):
    if step == n_outer or stop:
        write_field(step)
    n = len(ids); m = mach(step, stop)
    with h5py.File(f"res_nodes_{{step}}.h5", "w") as h:
        h.attrs["step"] = step
        h["MESH/COORD"] = coord[ids].ravel()
        h["MESH/NODE_ID"] = np.asarray(ids, dtype=np.int64)
        for k, v in (("ro", 0.1), ("Ux", m * A), ("Uy", 0.0), ("Uz", 0.0), ("P", 1e4), ("T", 300.0), ("sonic", A)):
            h["VALUE/" + k] = np.full(n, v)
    mf.write(f"{{step}},-1.0,-1.0,1.0,1.0,0.0,0.0\n"); mf.flush()
def resid_row(r):
    v = 1e-4 if sc.get("nan_at") is None or r < sc["nan_at"] else float("nan")
    vals = [v, v, v, 0.0, v] + ([1e-3] if sc.get("extra_col") else [])
    return f"{{r}},-1,outer_end," + ",".join(repr(float(x)) for x in vals) + "\n"
def wait_watcher(s):
    first = sc["first_judgeable"]
    if s < first or s % dt:
        return
    if sc.get("pace"):
        time.sleep(float(sc["pace"]))
    t0 = time.time()
    while time.time() - t0 < 30.0:
        try:
            rec = json.load(open("stop_rule.json"))
        except Exception:
            rec = None
        if rec and any(j.get("step") is not None and j["step"] >= s for j in rec.get("judgements", [])):
            return
        time.sleep(0.005)
write_field(0)
req_seen = None
i = 0
stop_now = False
print("fake forge v2 start", n_outer, dt, flush=True)
while i < n_outer:
    i += 1
    if not sc.get("ignore_stop") and os.path.exists("STOP_REQUEST"):
        req_seen = i if req_seen is None else req_seen
        stop_now = i - req_seen >= int(sc.get("stop_delay", 0))
    out = (i % dt == 0) or i == n_outer or stop_now
    if out:
        write_outputs(i, stop_now)
    hist.write(resid_row(i - 1))
    if out or (i - 1) % 7 == 0:
        hist.flush()
    if out and not stop_now:
        wait_watcher(i)
    if stop_now:
        break
hist.flush()
if stop_now:
    print(f"[coopStop] step {{i}}: STOP_REQUEST を検出", flush=True)
    if sc.get("stop_rc"):
        sys.exit(int(sc["stop_rc"]))
    k = i + int(sc.get("bad_stopped_at", 0))
    files = [f"res_{{i}}.h5", f"res_nodes_{{i}}.h5", "boundary_massflux.csv", "residual_history.csv"]
    js = {{"status": "ok", "step": k, "residual_step": k - 1, "request": open("STOP_REQUEST", encoding="utf-8").read()[:4096].rstrip(),
          "nStepOuter": n_outer, "files": [{{"path": f, "bytes": os.path.getsize(f)}} for f in files]}}
    open("STOPPED_AT.tmp", "w").write(json.dumps(js, indent=2, ensure_ascii=False) + "\n")
    os.replace("STOPPED_AT.tmp", "STOPPED_AT")
print("fake forge v2 end", flush=True)
sys.exit(0)
'''

EULER_DUMMY = TMP / "euler_dummy"
DELTA_CALLS = []


def fake_extract_delta_e(run_dir, step, euler_ref, out_dir, sens_factor=SR.YB_SENS_FACTOR):
    """δ_E の抽出 (`stop_rule.extract_delta_e`) の偽物: 呼ばれた step を記録して形だけ返す。"""
    DELTA_CALLS.append({"run": str(run_dir), "step": int(step), "euler_ref": str(euler_ref), "out_dir": str(out_dir)})
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    return {"status": "ok", "step": int(step), "dir": str(out_dir), "res_sha256": SR._sha256(Path(run_dir) / f"res_{step}.h5"),
            "delta_r_next_sha256": "試験の偽物", "seconds": 0.0,
            "band_sensitivity": {"factor": float(sens_factor), "max_abs_B_minus_A_mm": 0.0, "x_m_at": 0.0}}


def make_run(name: str, scenario: dict, *, hold: int | None = 38000, cap: int = 100000, crit: dict | None = None,
             cfg_updates: dict | None = None, cfg_text: str | None = None, nodes: list | None = None) -> Path:
    """偽の forge 用の v2 の run (prepare_ns の出力の形: prepare_info.json の stop_rule・solverConfig・settle_criteria.yaml・節点の表)。"""
    rd = TMP / name
    rd.mkdir()
    C = coords()
    np.save(rd / "fake_coord.npy", C)
    with h5py.File(rd / "nozzle.h5", "w") as h:
        h["MESH/COORD"] = C.ravel()
    c = crit if crit is not None else criteria(hold)
    first = c["grid"]["start_step"] + 2 * c["window"]["W"] + c["grid"]["out_step_interval"]
    (rd / "fake_scenario.json").write_text(json.dumps({"first_judgeable": first, **scenario}))
    cfg = cfg_text if cfg_text is not None else TEMPLATE_CFG
    cfg = RA._cfg_set(cfg, {"nStepOuter": 10000, **(cfg_updates or {})})
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(TEMPLATE_BC)
    (rd / "MESH_QUALITY.txt").write_text("VERDICT: PASS (試験の偽物)\n")
    (rd / ST.CRITERIA_FILE).write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    if nodes is None:
        X, R = SP.mesh_xr(rd / "nozzle.h5", NI, NJ, SCALE)
        nodes = SP.probe_node_ids(NI, NJ, X, R, ST.validate_criteria(c)).tolist()
    SP.write_node_probe_file(rd / "nodes.txt", nodes, header="試験の節点の表")
    EULER_DUMMY.mkdir(exist_ok=True)
    sr = SR.resolve_settings({}, nsteps_arg=None, n_main=10000, n_main_source="試験", out_int=DT, euler_ref=str(EULER_DUMMY), requested=None)
    sr["cap"] = cap
    info = {"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": 1.0, "physical_wall": {"domain": [0.0, 1.0]},
            "wall_thermal": {"mode": "adiabatic"}, "euler_ref": str(EULER_DUMMY), "nStepOuter": 10000, "stop_rule": sr}
    (rd / "prepare_info.json").write_text(json.dumps(info))
    return rd


# 見張りの時間 (stop_rule.watch_main_stage_v2)。stop_wait_s は「要求を置いてから run_case.sh (子プロセス) が終わるまで」の上限で、
# forge の終了の後の run_case.sh の後処理 (check_convergence.py・残差の図; python の起動 2 回) も含む。2.0 s では負荷のあるとき
# (CPU を 16 本の空回しで埋めた手元の再現、2026-10-16) に後処理が間に合わず、正しく止まった run まで timeout に分類された
# (REACHED・保持の既定・崩れ・終了コード 3 の 6 件)。正常に止まる場合は run_case.sh が終わった時点で見張りが抜けるので、上限を
# 大きくしても試験は長くならない。要求を無視する forge の試験 (timeout) だけは短い上限 TIMEOUT_STOP_WAIT_S を使う
# (偽の forge は要求の後も 1 格子点ごとに 0.4 s 以上かけて 10 格子点を回るので、負荷があっても先に上限が来る)。
WATCH = {"poll_s": 0.02, "stop_wait_s": 120.0, "term_wait_s": 3.0, "kill_wait_s": 5.0, "post_wait_s": 120.0}
TIMEOUT_STOP_WAIT_S = 2.0
# 後処理が stop_wait_s を超える試験 (§5.1 #5b): 要求から forge の終了までは偽の forge で 1 s 未満 (負荷があっても数 s)、後処理は
# check_convergence.py の前に SLOW_POST_S 眠るので、stop_wait_s を超えた時点で forge は終わっていて run_case.sh は後処理の途中
SLOW_STOP_WAIT_S = 6.0
SLOW_POST_S = 15


class slow_postprocess:
    """run_case.sh の後処理を遅らせる: PATH の先頭に python3 の包みを置き、check_convergence.py のときだけ seconds 眠ってから本物の
    python3 を exec する (見張りが run_case.sh に渡す環境 runner._ENV を一時的に書き換える; 偽の forge は sys.executable で動くので影響しない)。"""

    def __init__(self, seconds: float):
        self.seconds = seconds

    def __enter__(self):
        real = shutil.which("python3", path=RN._ENV.get("PATH"))
        d = TMP / "slow_python3"
        d.mkdir(exist_ok=True)
        w = d / "python3"
        w.write_text("#!/bin/sh\ncase \"$1\" in *check_convergence.py) sleep %s;; esac\nexec %s \"$@\"\n" % (self.seconds, real))
        w.chmod(0o755)
        self.saved = RN._ENV.get("PATH")
        RN._ENV["PATH"] = f"{d}:{self.saved}"
        return self

    def __exit__(self, *exc):
        RN._ENV["PATH"] = self.saved
        return False


def run_case(rd, watch: dict | None = None, **kw):
    rc = RA.run_staged_ns(rd, stages="none", stop_rule_watch={**WATCH, **(watch or {})}, **kw)
    rec = json.loads((rd / SR.RECORD).read_text()) if (rd / SR.RECORD).exists() else None
    return rc, rec, SR.read_state(rd)


def state_ok(rc, st, status):
    """戻り値の付加情報と STOP_RULE_STATE.json が同じ状態を持つか (method v2)。"""
    return (st is not None and st["status"] == status and st.get("method") == "v2" and st.get("version") == 2
            and getattr(rc, "stop_rule", None) is not None and rc.stop_rule["status"] == status
            and rc.stop_rule == {k: v for k, v in st.items() if k not in ("version", "plan")})


def replay_matches(rd) -> tuple:
    """オンラインの格子点ごとの判定 (stop_rule.json の judgements) が、終了後に新しい RunData で同じ終端を判定した結果と同じか
    (verdict と量の sha256)。戻り (全部同じか, 比べた数, 違いの例)。"""
    rec = json.loads((rd / SR.RECORD).read_text())
    lc = ST.load_criteria(rd)
    data = ST.RunData(rd, lc["criteria"])
    bad = []
    for j in rec["judgements"]:
        r = ST.judge(rd, j["step"], data=data, loaded=lc)
        if r["verdict"] != j["verdict"] or SR.quantities_digest(r) != j["digest"]:
            bad.append((j["step"], j["verdict"], r["verdict"]))
    return (not bad and len(rec["judgements"]) > 0), len(rec["judgements"]), bad[:3]


def cli_replay(rd, end_step: int) -> tuple:
    """`check_settled.py <run> --end-step N --json --no-write` の終了コードと JSON。"""
    p = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda" / "tools" / "check_settled.py"), str(rd), "--end-step", str(end_step),
                        "--json", "--no-write"], capture_output=True, text=True)
    lines = p.stdout.splitlines()
    try:
        a = lines.index("{")
        b = max(i for i, ln in enumerate(lines) if ln == "}")
        return p.returncode, json.loads("\n".join(lines[a:b + 1]))
    except ValueError:
        return p.returncode, None


# ====================================================================================== (a) 部品
def test_parts():
    print("== (a) 部品: 完成した行・オンラインの読み込み・STOP_REQUEST・STOPPED_AT・起動前の検査", flush=True)
    d = TMP / "parts"
    d.mkdir()
    p = d / "residual_history.csv"
    p.write_bytes(b"step,inner,phase,rms_ro\n0,-1,outer_begin,1\n0,-1,outer_end,1\n1,-1,outer_begin,1\n1,-1,oute")
    check("最後の完成した行は改行で終わる行 (書きかけの行を除く)", SR._last_complete_line(p) == b"1,-1,outer_begin,1")
    check("残差の行の完成: 行 0 の outer_end は完成、行 1 は outer_begin だけなので未完成",
          SR.residual_row_ready(p, 0) and not SR.residual_row_ready(p, 1))
    q = d / "boundary_massflux.csv"
    q.write_bytes(b"# c\nstep,mdot_inlet_1\n")
    check("流量の CSV: 見出しだけなら最後の step は None", SR.massflux_last_step(q) is None)
    q.write_bytes(b"# c\nstep,mdot_inlet_1\n2000,-1\n4000,-1\n6")
    check("流量の CSV: 書きかけの最後の行 (6…) を読まず最後の step は 4000", SR.massflux_last_step(q) == 4000)

    # オンラインの RunData: 書きかけの最後の行を読まない・path を元のファイルに戻す
    rd = d / "online_run"
    rd.mkdir()
    c = ST.validate_criteria(criteria())
    rows = "".join(f"{s},-1,outer_end,1e-4,1e-4,1e-4,0,1e-4\n" for s in range(0, 50))
    (rd / "residual_history.csv").write_text("step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n" + rows + "50,-1,outer_e")
    (rd / "boundary_massflux.csv").write_text("# c\nstep,mdot_inlet_1,raw_inlet_1,mdot_outlet_2,raw_outlet_2\n2000,-1,-1,1,1\n4")
    od = SR._online_data(ST, rd, c, rd / "tmp")
    R = od.residual("residual_history.csv")
    M = od.massflux("boundary_massflux.csv")
    check("オンラインの RunData: 書きかけの行を読まない (残差 50 行・読めない行 0、流量 1 行・読めない行 0)",
          len(R["steps"]) == 50 and R["bad"] == [] and list(M["steps"]) == [2000] and M["bad"] == [], f"{len(R['steps'])} {R['bad']} {M['bad']}")
    check("オンラインの RunData: 結果の path と使ったファイルの一覧は元のファイル (写しでない)",
          R["path"] == str(rd / "residual_history.csv") and M["path"] == str(rd / "boundary_massflux.csv")
          and set(od.used) == {str(rd / "residual_history.csv"), str(rd / "boundary_massflux.csv")}, str(sorted(od.used)))
    plain = ST.RunData(rd, c).massflux("boundary_massflux.csv")
    check("(対照) 親の RunData は書きかけの行を読めない行として数える (オンラインで判定不能にならないよう写しを使う理由)", plain["bad"] == [4])
    (rd / "residual_history.csv").write_text("step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n" + rows
                                             + "50,-1,outer_end,1e-4,1e-4,1e-4,0,1e-4\n")
    check("refresh の前は残差の読み込みを使い回す (50 行)", len(od.residual("residual_history.csv")["steps"]) == 50)
    od.refresh()
    check("refresh の後は読み直す (51 行)", len(od.residual("residual_history.csv")["steps"]) == 51)
    od2 = SR._online_data(ST, rd, c, rd / "tmp2")
    ok, e = raises(lambda: od2.residual("no_such.csv"), Exception, str(rd / "no_such.csv"))
    check("オンラインの RunData: 無いファイルは親と同じ判定不能 (文言は元の path)", ok, e[:100])

    # STOP_REQUEST の書き方
    t = SR.write_stop_request(d, {"kind": "settled", "reason": "試験 — 日本語", "request_step": 80000})
    check("STOP_REQUEST: 一時ファイル → rename、中身は JSON 1 行 (照合に使う文字列)",
          (d / SR.STOP_REQUEST).read_text() == t and json.loads(t)["request_step"] == 80000 and not (d / ".STOP_REQUEST.tmp").exists())
    ok, e = raises(lambda: SR.write_stop_request(d / "x", {"reason": "a" * 5000}), ValueError, "4 KiB")
    check("STOP_REQUEST: 4 KiB を越える中身は書かない", ok, e[:80])
    check("制御のファイルの残りを見つける", SR.leftover_control_files(d) == [SR.STOP_REQUEST])

    # STOPPED_AT の検査
    sd = d / "stopped"
    sd.mkdir()
    for f, n in (("res_80001.h5", 10), ("res_nodes_80001.h5", 5), ("boundary_massflux.csv", 7), ("residual_history.csv", 9)):
        (sd / f).write_bytes(b"x" * n)
    files = [{"path": f, "bytes": (sd / f).stat().st_size} for f in ("res_80001.h5", "res_nodes_80001.h5", "boundary_massflux.csv",
                                                                       "residual_history.csv")]
    good = {"status": "ok", "step": 80001, "residual_step": 80000, "request": t, "nStepOuter": 100000, "files": files}
    kw = dict(request_text=t, request_step=80000, cap=100000, crit=c)
    check("STOPPED_AT の検査: 正しい中身は問題なし", SR.check_stopped_at(sd, good, **kw) == [], str(SR.check_stopped_at(sd, good, **kw)))
    cases = {"無い": None, "residual_step のずれ": {**good, "residual_step": 80001}, "別の要求": {**good, "request": "{}"},
             "上限の違い": {**good, "nStepOuter": 90000}, "要求より前の step": {**good, "step": 79999, "residual_step": 79998},
             "大きさの違い": {**good, "files": [{**files[0], "bytes": 11}] + files[1:]},
             "節点の h5 が一覧に無い": {**good, "files": [f for f in files if f["path"] != "res_nodes_80001.h5"]}}
    for nm, st in cases.items():
        check(f"STOPPED_AT の検査: {nm} → 問題あり", len(SR.check_stopped_at(sd, st, **kw)) > 0, str(SR.check_stopped_at(sd, st, **kw))[:120])

    # 起動前の検査
    rd = make_run("pre_ok", {})
    cfg = yaml.safe_load((rd / "solverConfig.yaml").read_text())
    pf = SR.preflight_v2(rd, cfg, 100000)
    check("起動前の検査: そろった run は問題なし・最初に判定できる格子点 42,000・保持 38,000 (基準のファイル)",
          pf["problems"] == [] and pf["first_judgeable"] == 42000 and pf["hold"] == 38000, str(pf["problems"]))
    pf = SR.preflight_v2(rd, cfg, 60000)
    check("起動前の検査: 上限 60,000 では保持を終えられない (note、問題ではない)", pf["problems"] == [] and "必ず CENSORED" in (pf["note"] or ""))
    pf = SR.preflight_v2(rd, {**cfg, "time": {**cfg["time"], "outStepInterval": 1000}}, 100000)
    check("起動前の検査: 出力の間隔が格子と違う → 問題", any("outStepInterval" in x for x in pf["problems"]), str(pf["problems"])[:120])
    pf = SR.preflight_v2(rd, {**cfg, "output": {"level": 1, "nodeProbeFile": "nodes.txt"}}, 100000)
    check("起動前の検査: boundaryMassflux が無い → 問題", any("boundaryMassflux" in x for x in pf["problems"]))
    pf = SR.preflight_v2(rd, {**cfg, "output": {"level": 1, "boundaryMassflux": 1}}, 100000)
    check("起動前の検査: nodeProbeFile が無い → 問題", any("nodeProbeFile" in x for x in pf["problems"]))
    rd2 = make_run("pre_nodes_short", {}, nodes=[(NI - 1) * NJ + j for j in range(NJ - 1)])
    pf = SR.preflight_v2(rd2, cfg, 100000)
    check("起動前の検査: 節点の表に出口の列の節点が足りない → 問題", any("足りない" in x for x in pf["problems"]), str(pf["problems"])[:120])
    c_unset = criteria()
    c_unset["quantities"] = c_unset["quantities"][:2]
    rd3 = make_run("pre_unset", {}, crit=c_unset)
    pf = SR.preflight_v2(rd3, cfg, 100000)
    check("起動前の検査: 目的の量が無い → 評価対象が未設定 (問題)", any("評価対象が未設定" in x for x in pf["problems"]), str(pf["problems"])[:120])
    c_bad = criteria()
    c_bad["quantities"][0]["thresholds"]["bogus"] = 1
    rd4 = make_run("pre_bad", {}, crit=c_bad, nodes=[(NI - 1) * NJ + j for j in range(NJ)])
    pf = SR.preflight_v2(rd4, cfg, 100000)
    check("起動前の検査: 知らないキー → 基準の誤り (問題)", any("settle_criteria.yaml の誤り" in x for x in pf["problems"]), str(pf["problems"])[:120])
    c_ser = criteria()
    c_ser["quantities"][2] = {"kind": "exit_core_mach", "eta_range": [0.05, 0.7], "thresholds": {"tol": 1e-4},
                              "series_csv": {"path": "x.csv", "step_col": "step", "value_col": "M"}}
    rd5 = make_run("pre_series", {}, crit=c_ser)
    pf = SR.preflight_v2(rd5, cfg, 100000)
    check("起動前の検査: series_csv (再生用) はオンラインに使えない (問題)", any("series_csv" in x for x in pf["problems"]))
    (rd / SR.SNAPSHOT_REQUEST).write_text("")
    pf = SR.preflight_v2(rd, cfg, 100000)
    check("起動前の検査: SNAPSHOT_REQUEST の残り → 問題 (消さない)",
          any("SNAPSHOT_REQUEST" in x for x in pf["problems"]) and (rd / SR.SNAPSHOT_REQUEST).exists())


# ====================================================================================== (b) 偽の forge
def test_watcher():
    print("== (b) 偽の forge で停止規則 v2 の見張り", flush=True)
    real = SR.extract_delta_e
    SR.extract_delta_e = fake_extract_delta_e
    try:
        # --- REACHED (M1 の正例): 非出力の step 80,001 で止める ---
        DELTA_CALLS.clear()
        rd = make_run("reach_80001", {})
        rc, rec, st = run_case(rd)
        check("REACHED: 42,000 から SETTLED → 80,000 で保持 38,000 を満たして要求 → 非出力の step 80,001 で止まる (rc 0)",
              rc == 0 and state_ok(rc, st, "REACHED") and st["settled_since"] == 42000 and st["request_step"] == 80000
              and st["stop_step"] == 80001 and st["forge_rc"] == 0 and st["failures"] == [],
              json.dumps({k: (st or {}).get(k) for k in ("status", "settled_since", "request_step", "stop_step", "failures", "reason")},
                         ensure_ascii=False)[:300])
        req = json.loads((rd / SR.STOP_REQUEST).read_text())
        check("REACHED: STOP_REQUEST は JSON (理由・条件を満たし始めた step・保持・要求の格子点)",
              req["kind"] == "settled" and req["settled_since"] == 42000 and req["hold_steps"] == 38000 and req["request_step"] == 80000
              and "3 条件" in req["reason"])
        sa = json.loads((rd / SR.STOPPED_AT).read_text())
        check("REACHED: STOPPED_AT の request は置いた STOP_REQUEST と同じ・検査に通った",
              sa["request"] == (rd / SR.STOP_REQUEST).read_text() and st["stopped_at"]["problems"] == [] and st["stopped_at"]["step"] == 80001)
        ve = st["reeval"]
        check("REACHED: 再評価は 80,001 を end_step に固定 (格子の終端 80,000 + 停止の観測 80,001) で SETTLED、SETTLE_VERDICT.json の sha256",
              ve["verdict"] == "SETTLED" and ve["step"] == 80000 and ve["stop_step"] == 80001 and ve["end_step"] == 80001
              and ve["sha256"] == SR._sha256(rd / ST.VERDICT_FILE), json.dumps(ve, ensure_ascii=False)[:200])
        sv = json.loads((rd / ST.VERDICT_FILE).read_text())
        stop_checks = [c_["check"] for q in sv["quantities"] for c_ in q["checks"] if "停止の観測" in c_["check"]]
        check("REACHED: 再評価に停止の観測の別検査が入っている (残差・流量・出口コア M・流れ角)", len(stop_checks) >= 4, str(len(stop_checks)))
        ok, n, bad = replay_matches(rd)
        check("REACHED: オンラインの格子点ごとの判定 (42,000〜80,000) が終了後の再生と同じ (verdict・量の sha256)",
              ok and n == 20, f"{n} 点、違い {bad}")
        crc, cj = cli_replay(rd, 80001)
        check("REACHED: check_settled.py --end-step 80001 の判定が見張りの再評価と同じ (終了コード 0・verdict・量の sha256)",
              crc == 0 and cj is not None and cj["verdict"] == "SETTLED" and SR.quantities_digest(cj) == SR.quantities_digest(sv),
              f"rc {crc}")
        check("REACHED: δ_E の抽出は止めた step 80,001 の全場から 1 回だけ (Euler の参照つき)",
              [c_["step"] for c_ in DELTA_CALLS] == [80001] and st["delta_E"]["status"] == "ok" and st["delta_E"]["step"] == 80001,
              str(DELTA_CALLS))
        check("REACHED: θ_r・δ の時系列は作らない (stop_rule.json に rows・delta_series が無い)、全場は res_0 と止めた step だけ",
              "rows" not in rec and "delta_series" not in rec and SR.res_steps(rd) == [0, 80001], str(SR.res_steps(rd)))
        check("REACHED: forge を kill していない (stop_forge なし)・run_case.sh の後処理が走った",
              rec.get("stop_forge") is None and rec["postprocess"]["by"] == "run_case.sh")
        check("REACHED: solverConfig と stage_manifest の本段の nStepOuter = 上限 100,000",
              RA._cfg_value((rd / "solverConfig.yaml").read_text(), "nStepOuter") == 100000)

        # --- 保持の既定 (基準のファイルに hold が無い → 40,000) ---
        DELTA_CALLS.clear()
        rd = make_run("reach_default_hold", {}, hold=None)
        rc, rec, st = run_case(rd)
        check("保持の既定 40,000: 82,000 で要求・82,001 で REACHED", state_ok(rc, st, "REACHED") and st["hold_steps"] == 40000
              and st["request_step"] == 82000 and st["stop_step"] == 82001, f"{st and (st['status'], st['request_step'], st['stop_step'])}")

        # --- 要求の後に崩れる (i): 停止の観測だけが跳ねる ---
        DELTA_CALLS.clear()
        rd = make_run("break_at_stop", {"stop_jump": 0.01})
        rc, rec, st = run_case(rd)
        check("要求の後に崩れる (停止の step の出口 M が +0.01): 再評価 NOT_SETTLED → FAILED (reeval_not_settled)、δ_E は抽出しない (rc ≠ 0)",
              rc != 0 and state_ok(rc, st, "FAILED") and st["failures"] == ["reeval_not_settled"] and st["reeval"]["verdict"] == "NOT_SETTLED"
              and st["stop_step"] == 80001 and DELTA_CALLS == [] and st["delta_E"] is None, f"{st and (st['failures'], st['reeval'])}")

        # --- 要求の後に崩れる (ii): 要求を見てから 4,000 step 遅れて止まり、その間に出口 M が変わる ---
        rd = make_run("break_after_request", {"stop_delay": 4000, "jump_at": 80001})
        rc, rec, st = run_case(rd)
        after = [j for j in rec["judgements"] if j.get("after_request")]
        check("要求の後に崩れる (要求の後の格子点 82,000・84,000 が NOT_SETTLED): FAILED (broken_after_request・reeval_not_settled)、停止 84,001",
              state_ok(rc, st, "FAILED") and "broken_after_request" in st["failures"] and "reeval_not_settled" in st["failures"]
              and st["stop_step"] == 84001 and [j["step"] for j in after] == [82000, 84000] and all(j["verdict"] != "SETTLED" for j in after),
              f"{st and st['failures']} {[(j['step'], j['verdict']) for j in after]}")

        # --- forge が終了コード 3・STOPPED_AT なし ---
        rd = make_run("rc3", {"stop_rc": 3})
        rc, rec, st = run_case(rd)
        check("協調停止で forge が終了コード 3・STOPPED_AT なし → FAILED (forge_exit・no_stopped_at; rc 3)",
              rc == 3 and state_ok(rc, st, "FAILED") and set(st["failures"]) == {"forge_exit", "no_stopped_at"} and st["stop_step"] is None
              and st["forge_rc"] == 3, f"{st and (st['failures'], st['forge_rc'])}")

        # --- STOPPED_AT の中身が合わない ---
        rd = make_run("bad_stopped_at", {"bad_stopped_at": 1})
        rc, rec, st = run_case(rd)
        check("STOPPED_AT の step が書いた出力とずれている → FAILED (stopped_at_invalid)",
              state_ok(rc, st, "FAILED") and st["failures"] == ["stopped_at_invalid"] and st["stopped_at"]["problems"],
              f"{st and (st['failures'], st['stopped_at'])}")

        # --- 上限で止まる → CENSORED ---
        DELTA_CALLS.clear()
        rd = make_run("censored", {}, cap=60000)
        rc, rec, st = run_case(rd)
        check("上限 60,000 で保持 38,000 を終えられない → CENSORED (rc 0、保持の途中と記録、STOP_REQUEST を置かない)",
              rc == 0 and state_ok(rc, st, "CENSORED") and st["stop_step"] == 60000 and "保持の途中" in (st.get("reason") or "")
              and not (rd / SR.STOP_REQUEST).exists() and st["request_step"] is None, f"{st and (st['status'], st.get('reason'))}")
        check("CENSORED: δ_E は上限の step 60,000 の全場から 1 回、再評価は 60,000 で SETTLED (記録)",
              [c_["step"] for c_ in DELTA_CALLS] == [60000] and st["reeval"]["end_step"] == 60000 and st["reeval"]["verdict"] == "SETTLED")

        # --- 要求を無視する forge → 時間切れで kill ---
        rd = make_run("timeout", {"ignore_stop": True, "pace": 0.4})
        rc, rec, st = run_case(rd, watch={"stop_wait_s": TIMEOUT_STOP_WAIT_S})
        check("要求を無視する forge → stop_wait_s の後に kill → FAILED (timeout; REACHED にしない)",
              state_ok(rc, st, "FAILED") and "timeout" in st["failures"] and rec["stop_forge"]["stopped"]
              and rec["stop_forge"]["signals"][:1] == ["SIGTERM"], f"{st and st['failures']} {rec.get('stop_forge')}")
        check("要求を無視する forge: forge の終了を見ていない (stop_wait の記録なし) — 止まらない forge は今どおり時間切れ",
              "stop_wait" not in rec, str(rec.get("stop_wait")))

        # --- forge は止まったが run_case.sh の後処理が stop_wait_s を超える → 時間切れにしない (§5.1 #5b) ---
        DELTA_CALLS.clear()
        rd = make_run("slow_post", {})
        with slow_postprocess(SLOW_POST_S):
            rc, rec, st = run_case(rd, watch={"stop_wait_s": SLOW_STOP_WAIT_S})
        sw = rec.get("stop_wait") or {}
        check(f"後処理が遅い (check_convergence.py の前に {SLOW_POST_S} s) が forge は止まった → stop_wait_s {SLOW_STOP_WAIT_S} s を超えても "
              "REACHED (rc 0・80,001・failures なし)",
              rc == 0 and state_ok(rc, st, "REACHED") and st["stop_step"] == 80001 and st["failures"] == [] and st["forge_rc"] == 0,
              f"{st and (st['status'], st['failures'], st.get('reason'))}")
        check("後処理が遅い: 新しい経路を通った (stop_wait の記録: forge の終了 rc 0 を stop_wait_s の後に見た)・forge を kill していない",
              sw.get("forge_rc") == 0 and sw.get("waited_s", 0) > SLOW_STOP_WAIT_S and sw.get("stop_wait_s") == SLOW_STOP_WAIT_S
              and rec.get("stop_forge") is None, json.dumps(sw, ensure_ascii=False)[:300])
        check("後処理が遅い: run_case.sh の後処理が終わるまで待った (CONVERGENCE_VERDICT.txt を run_case.sh が書いた)・δ_E は 80,001 で 1 回",
              rec["postprocess"] == {"by": "run_case.sh"} and "CONVERGENCE_EXIT" in (rd / "CONVERGENCE_VERDICT.txt").read_text()
              and [c_["step"] for c_ in DELTA_CALLS] == [80001] and not rec.get("run_case_killed"), str(rec.get("postprocess")))

        # --- 非有限の残差 → DIVERGED の判定で止める ---
        rd = make_run("diverged", {"nan_at": 50000})
        rc, rec, st = run_case(rd)
        check("残差が 50,000 から非有限 → 52,000 の判定が DIVERGED → 要求 (diverged) → FAILED (diverged)",
              state_ok(rc, st, "FAILED") and st["failures"] == ["diverged"] and st["request_kind"] == "diverged" and st["request_step"] == 52000
              and st["stop_step"] == 52001 and st["reeval"]["verdict"] == "DIVERGED", f"{st and (st['failures'], st['request_step'])}")

        # --- 基準に無い残差の列が値を持つ → 未設定の判定で止める ---
        rd = make_run("unset_col", {"extra_col": True})
        rc, rec, st = run_case(rd)
        check("基準に無い残差の列 rms_roK が値を持つ → 42,000 の判定が INDETERMINATE (未設定) → 要求で止めて FAILED (criteria_unset)",
              state_ok(rc, st, "FAILED") and st["failures"] == ["criteria_unset"] and st["request_step"] == 42000 and st["stop_step"] == 42001,
              f"{st and (st['failures'], st['request_step'])}")

        # --- 起動前の制御のファイルの残り・設定の不足 ---
        rd = make_run("leftover", {})
        (rd / SR.STOP_REQUEST).write_text("前の run の残り")
        ok, msg = raises(lambda: run_case(rd), RuntimeError, "STOP_REQUEST")
        check("起動前に STOP_REQUEST が残っている → forge を起動せず RuntimeError、ファイルは消さない・状態も書かない",
              ok and (rd / SR.STOP_REQUEST).read_text() == "前の run の残り" and not (rd / "run_case_stdout.log").exists()
              and SR.read_state(rd) is None, msg[:160])
        rd = make_run("no_massflux_cfg", {}, cfg_text=TEMPLATE_CFG.replace(", boundaryMassflux: 1", ""))
        ok, msg = raises(lambda: run_case(rd), RuntimeError, "boundaryMassflux")
        check("本段の設定に output.boundaryMassflux が無い → forge を起動せず RuntimeError", ok and not (rd / "run_case_stdout.log").exists(),
              msg[:160])

        # --- 分岐: settle_criteria.yaml の無い run は旧方式・stop_rule=False ---
        rd = make_run("no_criteria", {})
        (rd / ST.CRITERIA_FILE).unlink()
        check("settle_criteria.yaml の無い run → v2 を使わない (`_stop_rule_v2_for_main` が None、旧方式が決める)",
              RA._stop_rule_v2_for_main(rd, None, (rd / "solverConfig.yaml").read_text()) is None)
        rd = make_run("arg_false", {})
        rc, rec, st = run_case(rd, stop_rule=False)
        check("v2 の run でも stop_rule=False → 見張らない・状態 NOT_USED (method v2)",
              st is not None and st["status"] == "NOT_USED" and st["method"] == "v2" and rec is None and not (rd / SR.STOP_REQUEST).exists())
    finally:
        SR.extract_delta_e = real


def main() -> int:
    FAKE_FORGE.write_text(FAKE_SRC.format(py=sys.executable))
    FAKE_FORGE.chmod(0o755)
    test_parts()
    if not QUICK:
        test_watcher()
    print(f"\n{'全部合格' if FAIL == 0 else f'不合格 {FAIL} 件'} (一時ディレクトリ {TMP})", flush=True)
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
