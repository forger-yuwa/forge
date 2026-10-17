#!/usr/bin/env python3
"""計算の後の判定の確定 (`evaluate/judge.py`) の統合試験 (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12 後半、
plan 段 judge の m8)。本物の forge は起動しない: 本物の `run_case.sh` に偽の forge (`FORGE_BIN`) を渡す。

(a) 段階起動 (偽の forge): 親が `judge.Staged` で包み、全部の段を `runner.run_forge` (`FORGE_JUDGE_DEFER=1`) で起動する。S1 は manifest の
    無い新しい run でも段の記録だけ (`state: deferred`)、2 段目も段の記録、抜けたときに総合を 1 度だけ確定 (`state: final`・段の記録 2 つを
    参考に残す・`JUDGE_PENDING.json` を消す)。起動の後のエコー・残差のプロット・収束の副判定は残る。
(b) 段の途中の例外: 1 段目の後に例外 → 例外は外へ出るが、起動していたので確定する。
(c) 関門で未起動: 投入の関門が止めた (bcondConfig.yaml が無い) → 起動していないので確定しない (RUN_VERDICT を書かない)、
    `JUDGE_PENDING.json` は消す。粘性の壁の outputHDFflg 0 も関門で止まる (確定しない)。
(d) 系列 `judge_series.csv` (`write_axismach_series`): 停止規則の記録 stop_rule.json の行を STOP_RULE_STATE.json の record_sha256 と
    一致するときだけ使い (θ_r の列も)、一致しなければ本段の res_<N>.h5 (N > 0) から同じ式で作る (θ_r の列は無い)。量の値は停止規則の
    `snapshot` と一致。非有限の出力の行は nan。
(e) 機種の既定の仕様 (`spec_*`) が `judge_run.normalize_spec` を通り、量の名前・抽出器・閾値が期待どおり。誤った仕様は書かない。
(f) `same_numerics` (継続の run の `--from-floor` を書くか): 同じ config は True、convMethod (hard キー) が違えば False。
(g) 判定器の例外: judge_run.py が例外で終わっても `finalize` は例外にせず、RUN_VERDICT は state error・INDETERMINATE。
(h) 古い判定: 確定の後に系列を書き換えると `--check-fresh` が STALE (Stop フックが止める側)。
(i) axis-Mach Euler の準定常 (§5.1 #14b): 本段の出力 6 枚 (12,000 step を 2,000 の間隔) で `spec_axismach(euler=True)` の窓
    (tail 0.4・min_points 4) を満たし STEADY・DRIFTING を判定できる。3 枚は窓の不足 (判定不能)。
(j) criterion stop_rule (§5.1 #14b): 停止規則の記録 stop_rule.json → judge_series.csv → 停止規則の `judge_rows` そのもので判定した値・合否が
    記録の行に judge_rows を直に当てた結果と同じ (到達 = first_stop が最後の行で REACHED、未到達 = 届かない)。停止規則の到達との照合
    (stop_rule_cross) が一致。stop_rule.json・settings が無ければ判定不能。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_judge_finalize_tests.py   (約 1 分)
"""
from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
TOOLS = ROOT / "solver_density_cuda" / "tools"
TMP = Path(tempfile.mkdtemp(prefix="judge_finalize_tests_"))
FAKE = TMP / "fake_forge.py"
os.environ["FORGE_BIN"] = str(FAKE)                    # runner._ENV は import 時の環境を写すので import の前に
os.environ.pop("FORGE_JUDGE_DEFER", None)
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, str(TOOLS))

from forge_design.evaluate import judge as J  # noqa: E402
from forge_design.evaluate import runner as RR  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
import judge_run as JR  # noqa: E402
import stage_manifest as SM  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""), flush=True)
    if not cond:
        FAIL += 1


FAKE.write_text('''#!%s
# 偽の forge (試験用): 残差の履歴・体積の出力 (h0 つき) を書き、正常終了の行を出す
import os, sys
import h5py, numpy as np
open("FORGE_STARTED", "a").write("started\\n")
with open("residual_history.csv", "w") as f:
    f.write("step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\\n")
    for i in range(200):
        v = 10.0 ** (-5.0 * i / 199.0)
        for ph in ("outer_begin", "outer_end"):
            f.write("%%d,-1,%%s,%%s\\n" %% (i, ph, ",".join(["%%.6e" %% v] * 5)))
n = len(h5py.File("mesh.h5", "r")["MESH/COORD"]) // 3
for s in (100, 200):
    with h5py.File("res_%%d.h5" %% s, "w") as h:
        for k, val in (("ro", 1.2), ("P", 1e5), ("T", 300.0), ("roUx", 1.0), ("h0", 1004.5 * 300.0)):
            h.create_dataset("VALUE/" + k, data=np.full(n, val))
print("Time = 0.010 s (wall, 200 steps, 0.05 ms/step)")
sys.exit(int(os.environ.get("FAKE_FORGE_RC", "0")))
''' % sys.executable)
FAKE.chmod(FAKE.stat().st_mode | stat.S_IXUSR)
QS_NA = {"judge_spec_version": 1, "quasisteady": {"not_applicable": "試験: 報告する派生量が無い"}}
CFG = ('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node"}\n'
       'physProp: {thermalMethod: 0, cp: 1004.5, gamma: 1.4, viscMethod: 0, visc: %s, thermCond: 0.0}\n')


def small_run(tag, viscous=False, wall_out=0, bcond=True) -> Path:
    """関門を通る node の小さな run (test_judge_run の rc_run と同じ置き方、上辺が壁)。"""
    d = TMP / tag
    d.mkdir()
    ni, nj = 11, 4
    X, Y = np.meshgrid(np.linspace(0.0, 1.0, ni), np.linspace(0.0, 0.2, nj), indexing="ij")
    nid = lambda i, j: i * nj + j  # noqa: E731
    quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    wall = [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)]
    with h5py.File(d / "mesh.h5", "w") as f:
        f.create_dataset("MESH/COORD", data=np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)].ravel())
        f.create_dataset("MESH/CONNE", data=np.array([x for n in range(ni * nj) for x in (1, n)], np.int64))
        f.create_dataset("VIZMESH/CONNE", data=np.array([x for q in quads for x in [5] + q], np.int64))
        f["VIZMESH"].attrs["nVizCells"] = len(quads)
        f.create_dataset("VALUE/wall_dist", data=0.2 - Y.ravel())
        f.create_dataset("BCONDS/3/vizBfaceNodes", data=np.array(wall, np.int64).ravel())
        f.create_dataset("BCONDS/3/vizBfaceSizes", data=np.full(len(wall), 2, np.int64))
        f["BCONDS/3"].attrs["bcondKind"] = "wall"
    (d / "solverConfig.yaml").write_text(CFG % ("1.8e-5" if viscous else "0.0"))
    if bcond:
        (d / "bcondConfig.yaml").write_text("wall: {physID: 3, kind: wall, outputHDFflg: %d, ints: , floats: }\n" % wall_out)
    (d / "JUDGE_SPEC.json").write_text(json.dumps(QS_NA))
    return d


def verdict(d: Path):
    p = d / "RUN_VERDICT.json"
    return json.loads(p.read_text()) if p.exists() else None


# ------------------------------------------------------------------------------------------------ (a) 段階起動
print("== (a) 段階起動 (偽の forge、FORGE_JUDGE_DEFER=1、親が 1 度だけ確定)", flush=True)
d = small_run("staged")
cfg = (d / "solverConfig.yaml").read_text()
bc = (d / "bcondConfig.yaml").read_text()
sm = SM.StageManifest(str(d))
seen = {}
with J.Staged(d, "試験 (a)") as S:
    check("入ったら JUDGE_PENDING.json (pid・host) を置く", json.loads((d / J.PENDING).read_text())["pid"] == os.getpid())
    rc1 = RR.run_forge(d)
    v1 = verdict(d)
    seen["S1"] = v1
    shutil.copy(d / "residual_history.csv", d / "residual_history_S1.csv")
    sm.add("S1", cfg, bc, history="residual_history_S1.csv")
    sm.write()
    for f in d.glob("res_*"):
        f.unlink()
    rc2 = RR.run_forge(d)
    seen["main"] = verdict(d)
    shutil.copy(d / "residual_history.csv", d / "residual_history_main.csv")
    sm.add("main", cfg, bc, history="residual_history_main.csv")
    sm.write()
v = verdict(d)
check("S1 (manifest の無い新しい run) でも run_case.sh は段の記録だけ (state deferred; FORGE_JUDGE_DEFER=1)",
      rc1 == 0 and seen["S1"] and seen["S1"]["state"] == "deferred" and "FORGE_JUDGE_DEFER" in seen["S1"]["state_note"],
      json.dumps(seen["S1"] and seen["S1"].get("state_note"), ensure_ascii=False))
check("2 段目も段の記録 (段 2 つ)", rc2 == 0 and seen["main"]["state"] == "deferred" and len(seen["main"]["stages"]) == 2)
check("抜けたときに総合を 1 度だけ確定 (state final・段の記録 2 つを参考に残す・収束は --segment で PASS・総合 PASS)",
      v["state"] == "final" and len(v.get("stages") or []) == 2 and v["checks"]["convergence"]["status"] == "PASS"
      and "segment" in v["checks"]["convergence"]["reason"] and v["overall"] == "PASS" and S.result and S.result["overall"] == "PASS",
      v.get("overall_reason"))
check("確定の後は JUDGE_PENDING.json が無い", not (d / J.PENDING).exists())
check("エコー LAUNCH_ECHO.txt・残差のプロット・収束の副判定 CONVERGENCE_VERDICT.txt は今のまま",
      (d / "LAUNCH_ECHO.txt").exists() and (d / "residual_history.png").exists() and (d / "CONVERGENCE_VERDICT.txt").exists())
check("run_case.sh の終了コードは solver のまま (run_forge の戻り 0; 判定の終了コードを混ぜない)", rc1 == 0 and rc2 == 0)
r = subprocess.run(["python3", str(TOOLS / "judge_run.py"), str(d), "--check-fresh"], capture_output=True, text=True)
check("確定の直後は --check-fresh が FRESH", r.returncode == 0, r.stdout)

# ------------------------------------------------------------------------------------------------ (b) 例外
print("== (b) 段の途中の例外", flush=True)
d = small_run("exc")
raised = False
try:
    with J.Staged(d, "試験 (b)"):
        RR.run_forge(d)
        raise RuntimeError("段終了ゲートで不合格 (試験)")
except RuntimeError:
    raised = True
v = verdict(d)
check("例外は外へ出る・起動していたので確定する (state final)・JUDGE_PENDING.json は消す",
      raised and v and v["state"] == "final" and not (d / J.PENDING).exists(), v and v.get("state"))

# ------------------------------------------------------------------------------------------------ (c) 関門で未起動
print("== (c) 関門で未起動", flush=True)
d = small_run("gate", bcond=False)
with J.Staged(d, "試験 (c)") as S:
    rc = RR.run_forge(d)
check("bcondConfig.yaml が無い → 関門が止める (run_forge ≠ 0・forge 未起動)・確定しない (RUN_VERDICT を書かない)・JUDGE_PENDING.json を消す",
      rc != 0 and not (d / "FORGE_STARTED").exists() and verdict(d) is None and S.result is None and not (d / J.PENDING).exists(), f"rc {rc}")
d = small_run("gate_wall", viscous=True, wall_out=0)
with J.Staged(d, "試験 (c)") as S:
    rc = RR.run_forge(d)
gate = (d / "LAUNCH_GATE.txt").read_text() if (d / "LAUNCH_GATE.txt").exists() else ""
check("定数粘性の壁の outputHDFflg 0 → 関門が止める ([壁の出力] の FAIL)・確定しない",
      rc != 0 and not (d / "FORGE_STARTED").exists() and "[壁の出力]" in gate and verdict(d) is None, f"rc {rc}")

# ------------------------------------------------------------------------------------------------ (d) axis-Mach の系列
print("== (d) axis-Mach の系列 judge_series.csv", flush=True)
NI, NJ, SCALE, X_E = 41, 9, 0.01, 30.0


def axis_run(tag) -> Path:
    d = TMP / tag
    d.mkdir()
    x = np.linspace(-5.0, 35.0, NI)
    r_w = 1.0 + 0.02 * (x + 5.0)
    X = np.repeat(x[:, None], NJ, axis=1)
    R = r_w[:, None] * np.linspace(0.0, 1.0, NJ)[None, :]
    with h5py.File(d / "nozzle.h5", "w") as f:
        f.create_dataset("MESH/COORD", data=np.c_[X.ravel() * SCALE, R.ravel() * SCALE, np.zeros(NI * NJ)].ravel())
    (d / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": X_E,
                                                     "physical_wall": {"domain": [-5.0, 34.0]}}))
    for k, s in enumerate((0, 2500, 5000, 7500, 10000)):
        a = 1.0 + 0.01 * k
        with h5py.File(d / f"res_{s}.h5", "w") as f:
            n = NI * NJ
            f["VALUE/ro"] = np.full(n, 1.2 * a)
            f["VALUE/Ux"] = np.full(n, 500.0)
            f["VALUE/Uy"] = np.zeros(n)
            f["VALUE/sonic"] = np.full(n, 300.0)
            f["VALUE/P"] = np.full(n, 1e5)
    return d


d = axis_run("series_snap")
J.write_axismach_series(d)
meta = json.loads((d / J.SERIES_META).read_text())
rows = list(__import__("csv").DictReader(open(d / J.SERIES)))
geo = SR.RunGeometry(d)
want = SR.snapshot(geo, 5000, None)
got = [r for r in rows if int(r["step"]) == 5000][0]
check("停止規則の記録が無い → 本段の res_<N>.h5 (N > 0) から (step 0 を入れない・θ_r の列は無い)",
      meta["source"]["kind"].startswith("res_") and [int(r["step"]) for r in rows] == [2500, 5000, 7500, 10000]
      and meta["columns"] == ["exit_core_M", "mdot_rt", "r"], json.dumps(meta, ensure_ascii=False))
check("量の値は停止規則の snapshot と同じ式 (出口コア M・流量・質量の収支)",
      all(abs(float(got[k]) - want[k]) <= 1e-12 * max(1.0, abs(want[k])) for k in ("exit_core_M", "mdot_rt", "r")),
      json.dumps({k: (got[k], want[k]) for k in ("exit_core_M", "mdot_rt", "r")}))
rows_sr = [{"step": 2500, "nonfinite_all": 0, "r": 1e-4, "exit_core_M": 2.0, "mdot_rt": 3.0, "theta_r_0.4": 0.01, "x_theta_r_0.4": 13.6},
           {"step": 5000, "nonfinite_all": 3}]
SR.save_record(d, {"status": "REACHED", "rows": rows_sr})
SR.write_state(d, {"status": "REACHED", "record_sha256": SR._sha256(d / SR.RECORD)})
J.write_axismach_series(d)
meta = json.loads((d / J.SERIES_META).read_text())
rows = list(__import__("csv").DictReader(open(d / J.SERIES)))
check("停止規則の記録 (record_sha256 が一致) → 記録の行から (θ_r の列も、x_θ の列は入れない)・非有限の行は nan",
      meta["source"]["kind"] == "stop_rule.json" and meta["columns"] == ["exit_core_M", "mdot_rt", "r", "theta_r_0.4"]
      and rows[0]["theta_r_0.4"] == repr(0.01) and rows[1]["r"] == "nan", json.dumps({"meta": meta, "rows": rows}, ensure_ascii=False))
SR.save_record(d, {"status": "REACHED", "rows": rows_sr + [{"step": 7500, "nonfinite_all": 0, "r": 0.0}]})
J.write_axismach_series(d)
meta = json.loads((d / J.SERIES_META).read_text())
check("stop_rule.json が状態の記録と違う (sha256 の不一致) → 使わずに res から (出典が停止規則と一致しない)",
      meta["source"]["kind"].startswith("res_") and "sha256" in meta.get("stop_rule_rejected", ""), json.dumps(meta, ensure_ascii=False))

# ------------------------------------------------------------------------------------------------ (e) 機種の既定の仕様
print("== (e) 機種の既定の仕様 (JUDGE_SPEC.json)", flush=True)
specs = {"axismach_ns": J.spec_axismach("g", stop_rule=SR.resolve_settings({}, nsteps_arg=None, n_main=60000, n_main_source="試験",
                                                                           out_int=2500, euler_ref="x")),
         "axismach_ns_nosr": J.spec_axismach("g"), "axismach_euler": J.spec_axismach("g", euler=True),
         "sern_rans": J.spec_sern("g", True), "sern_euler": J.spec_sern("g", False), "bell_tp": J.spec_bell("g"),
         "snapshot": J.spec_snapshot("g"), "snapshot_ff": J.spec_snapshot("g", from_floor="../cont_1")}
names = {}
for k, sp in specs.items():
    try:
        JR.normalize_spec(sp, True)
        names[k] = [q["name"] for q in sp["quasisteady"]["quantities"]]
    except Exception as e:  # noqa: BLE001
        names[k] = f"誤り: {e}"
check("全機種の仕様が normalize_spec を通る", all(isinstance(v, list) for v in names.values()), json.dumps(names, ensure_ascii=False))
qns = specs["axismach_ns"]["quasisteady"]["quantities"]
check("axis-Mach NS (停止規則あり): 出口コア M・質量の収支・θ_r × 3 を criterion stop_rule (役 relative・r・theta、窓・閾値は書かない; #14a (2))",
      names["axismach_ns"] == ["exit_core_M", "mass_imbalance_r", "theta_r_0.4", "theta_r_0.7", "theta_r_0.95"]
      and all(q["criterion"] == "stop_rule" and "thresholds" not in q and "window" not in q and q["accept"] == ["STEADY"] for q in qns)
      and [q["stop_rule_role"] for q in qns] == ["relative", "r", "theta", "theta", "theta"],
      json.dumps([{k: q.get(k) for k in ("name", "criterion", "stop_rule_role")} for q in qns], ensure_ascii=False))
q0 = {q["name"]: q for q in specs["axismach_ns_nosr"]["quasisteady"]["quantities"]}
check("axis-Mach NS (停止規則なし)・Euler: θ_r を入れない (理由を notes に)、M は relative・r は absolute (停止規則の既定の閾値) のまま",
      names["axismach_ns_nosr"] == names["axismach_euler"] == ["exit_core_M", "mass_imbalance_r"]
      and any("θ_r は量に入れない" in n for n in specs["axismach_ns_nosr"]["notes"])
      and any("INDETERMINATE" in n for n in specs["axismach_euler"]["notes"])
      and q0["exit_core_M"]["criterion"] == "relative" and q0["mass_imbalance_r"]["criterion"] == "absolute"
      and q0["mass_imbalance_r"]["thresholds"] == {"abs_mean": SR.R_MEAN_MAX, "abs_amp": SR.R_AMP_MAX}
      and all(q["criterion"] != "stop_rule" for q in specs["axismach_euler"]["quasisteady"]["quantities"]))
check("SERN: 粘性は C_T_with_shear + C_T・C_L・C_M、非粘性は C_T・C_L・C_M (sern_forces の閾値 0.02/0.05)",
      names["sern_rans"] == ["C_T_with_shear", "C_T", "C_L", "C_M"] and names["sern_euler"] == ["C_T", "C_L", "C_M"]
      and specs["sern_rans"]["quasisteady"]["quantities"][0]["thresholds"] == {"drift": 0.02, "osc": 0.05})
check("ベル TP: tp_gates (ii) の列と閾値 (推力系 1.5e-4・ṁ 1e-3)",
      names["bell_tp"] == ["thrust_exit_N", "thrust_wall_N", "eta", "eta_wall", "mdot_exit_kg_s"]
      and specs["bell_tp"]["quasisteady"]["quantities"][0]["thresholds"] == {"drift": 1.5e-4, "osc": 1.5e-4}
      and specs["bell_tp"]["quasisteady"]["quantities"][-1]["thresholds"] == {"drift": 1e-3, "osc": 1e-3})
check("CPG・風洞: snapshot の shock・machmax・pmax、継続は from_floor",
      names["snapshot"] == ["shock", "machmax", "pmax"] and specs["snapshot_ff"]["convergence"] == {"mode": "from_floor", "reference": "../cont_1"})
d = TMP / "spec_bad"
d.mkdir()
try:
    J.write_spec(d, {"judge_spec_version": 1, "quasisteady": {"quantities": []}})
    bad = False
except ValueError:
    bad = True
check("誤った仕様は書かない (ValueError)", bad and not (d / J.SPEC).exists())

from forge_design.report import nozzle_report as NR  # noqa: E402
vp = TMP / "verdicts.json"
vp.write_text(json.dumps({"exit_core_M": {"verdict": "INDETERMINATE", "note": "窓の不足"}}))
ok_ind = NR.load_verdicts(vp)["items"]["exit_core_M"]["verdict"] == "INDETERMINATE"
vp.write_text(json.dumps({"exit_core_M": {"verdict": "STEADYISH"}}))
try:
    NR.load_verdicts(vp)
    ok_bad = False
except ValueError:
    ok_bad = True
check("nozzle_report --verdicts は INDETERMINATE (判定不能) を受ける・未知の綴りは拒否のまま", ok_ind and ok_bad)

# ------------------------------------------------------------------------------------------------ (f) same_numerics
print("== (f) 継続の run の --from-floor (same_numerics)", flush=True)
c1 = 'space: {convMethod: 1, limiter: 2}\nturbulence: {model: "sst", wallTreatmentSST: 0}\ntime: {last: {nStepOuter: 1000}}\n'
c2 = c1.replace("nStepOuter: 1000", "nStepOuter: 5000")
c3 = c1.replace("convMethod: 1", "convMethod: 0")
check("同じ数値設定 (反復数だけ違う) → True、convMethod が違う (soft → 本段) → False",
      J.same_numerics(c1, c2, "w: {}\n", "w: {}\n") and not J.same_numerics(c1, c3, "w: {}\n", "w: {}\n"))

# ------------------------------------------------------------------------------------------------ (g) 判定器の例外
print("== (g) 判定器の例外", flush=True)
d = small_run("judge_exc")
RR.run_forge(d)
(d / "JUDGE_SPEC.json").write_text(json.dumps(QS_NA))
orig_tools = J._tools
boom = TMP / "boom_tools"
boom.mkdir()
for f in TOOLS.glob("*.py"):
    (boom / f.name).symlink_to(f)
(boom / "judge_run.py").unlink()
(boom / "judge_run.py").write_text((TOOLS / "judge_run.py").read_text().replace(
    "def check_physical(run):", "def check_physical(run):\n    raise RuntimeError('試験の例外')\n\n\ndef _unused_check_physical(run):", 1))
J._tools = lambda: boom
res = J.finalize(d, generator="試験 (g)")
J._tools = orig_tools
v = verdict(d)
check("judge_run.py が例外 (終了 3) でも finalize は例外にしない・RUN_VERDICT は state error・INDETERMINATE",
      res["rc"] == 3 and v["state"] == "error" and v["overall"] == "INDETERMINATE", json.dumps({k: res.get(k) for k in ("rc", "overall")}))

# ------------------------------------------------------------------------------------------------ (h) 古い判定
print("== (h) 古い判定", flush=True)
d = small_run("stale")
(d / "JUDGE_SPEC.json").write_text(json.dumps({"judge_spec_version": 1, "quasisteady": {"quantities": [
    {"name": "r", "definition": "試験", "unit": "-", "region": "全体",
     "extractor": {"kind": "csv_column", "file": "judge_series.csv", "column": "r"}, "criterion": "absolute",
     "thresholds": {"abs_mean": 5e-4, "abs_amp": 5e-3}}]}}))


def series():
    with open(d / "judge_series.csv", "w") as f:
        f.write("step,r\n" + "".join("%d,%r\n" % (s, 1e-5) for s in range(100, 2001, 100)))


with J.Staged(d, "試験 (h)", series=series) as S:
    RR.run_forge(d)
v = verdict(d)
check("系列を確定の時点に書く (series) → 準定常の量を判定 (PASS)",
      v["state"] == "final" and v["checks"]["quasisteady"]["status"] == "PASS", v["checks"]["quasisteady"]["reason"])
with open(d / "judge_series.csv", "a") as f:
    f.write("2100,0.5\n")
r = subprocess.run(["python3", str(TOOLS / "judge_run.py"), str(d), "--check-fresh"], capture_output=True, text=True)
check("確定の後に系列が変わった → --check-fresh が STALE (Stop フックが止める)", r.returncode == 1 and "judge_series.csv" in r.stdout, r.stdout)

# ------------------------------------------------------------------------------------------------ (i) Euler の 6 枚
print("== (i) axis-Mach Euler の準定常: 本段の出力 6 枚で窓を満たす (spec_axismach(euler=True) の tail 0.4・min_points 4; #14b)", flush=True)


def write_series_csv(d: Path, rows: list, cols: list) -> None:
    with open(d / J.SERIES, "w") as f:
        f.write("step," + ",".join(cols) + "\n")
        for r_ in rows:
            f.write("%d," % r_["step"] + ",".join(repr(float(r_[c])) for c in cols) + "\n")


def qs(d: Path) -> dict:
    """準定常の項目 (judge_run.check_quasisteady を作業ツリーの仕様・停止の状態で)。"""
    return JR.check_quasisteady(str(d), JR.load_spec(str(d)), JR.read_stop_state(str(d)))


d = TMP / "euler6"
d.mkdir()
J.write_spec(d, J.spec_axismach("試験 (i)", euler=True))
steps6 = [2000 * k for k in range(1, 7)]                 # 12,000 step の本段を n/6 = 2,000 の間隔で
write_series_csv(d, [{"step": s, "exit_core_M": 4.2 * (1 + 1e-4 * (-1) ** k), "r": 2e-4 * (-1) ** k} for k, s in enumerate(steps6)],
                 ["exit_core_M", "r"])
c = qs(d)
qq = c["quantities"]
check("6 枚 (2,000〜12,000 step) で出口コア M・r が判定できる (窓の不足でない; 末尾窓 3 点・STEADY = PASS)",
      c["status"] == "PASS" and all(qq[k]["n"] == 6 and qq[k]["tail_points"] == 3 and qq[k]["verdict"] == "STEADY" for k in qq),
      json.dumps({k: (q["status"], q["reason"][:120]) for k, q in qq.items()}, ensure_ascii=False))
write_series_csv(d, [{"step": s, "exit_core_M": 4.2 * (1 + 0.05 * k), "r": 2e-4 * (-1) ** k} for k, s in enumerate(steps6)],
                 ["exit_core_M", "r"])
c = qs(d)
check("6 枚でドリフトも判定できる (出口コア M が 5 %/出力で増える → DRIFTING = FAIL; 判定不能でない)",
      c["status"] == "FAIL" and c["quantities"]["exit_core_M"]["verdict"] == "DRIFTING" and c["quantities"]["mass_imbalance_r"]["status"] == "PASS",
      c["reason"])
write_series_csv(d, [{"step": s, "exit_core_M": 4.2, "r": 2e-4} for s in (4000, 8000, 12000)], ["exit_core_M", "r"])
c = qs(d)
check("3 枚 (n/3 の間隔) は窓の不足 → 判定不能 (#12b (1) の状態)",
      c["status"] == "INDETERMINATE" and all("窓の不足" in q["reason"] for q in c["quantities"].values()), c["reason"])

# ------------------------------------------------------------------------------------------------ (j) criterion stop_rule
print("== (j) criterion stop_rule: 停止規則の記録 → judge_series.csv → judge_rows そのもの (到達と未到達; #14b)", flush=True)
SRS = SR.resolve_settings({}, nsteps_arg=None, n_main=60000, n_main_source="試験", out_int=2500, euler_ref="x")
TH = {k: SRS[k] for k in SR.thresholds()}
TKEYS = [SR.theta_key(f) for f in TH["x_fracs"]]


def sr_rows(n: int, drift: float = 0.0, m_osc: float = 0.0) -> list:
    """停止規則の記録の行 (本段の出力 2,500 step ごと): r は ±1e-4、θ_r は断面ごとの値 × (1 ± 1e-4 + drift·i)、出口コア M は
    5.5 × (1 ± m_osc) (i = 22 だけ 0: 25 行の最後の窓 (i = 16〜24) の前半・後半の平均を揃え、D でなく振幅だけが外れるように)。"""
    rows = []
    for i in range(n):
        p = 0 if i == 22 else (-1) ** i
        row = {"step": 2500 * (i + 1), "nonfinite_all": 0, "r": 1e-4 * (-1) ** i, "exit_core_M": 5.5 * (1 + m_osc * p), "mdot_rt": 3.0}
        for f in TH["x_fracs"]:
            row[SR.theta_key(f)] = 0.02 * f * (1 + 1e-4 * (-1) ** i + drift * i)
            row["x_" + SR.theta_key(f)] = 30.0 * f
        rows.append(row)
    return rows


def sr_run(tag: str, rows: list, status: str) -> Path:
    d = TMP / tag
    d.mkdir()
    J.write_spec(d, J.spec_axismach("試験 (j)", stop_rule=SRS))
    SR.save_record(d, {"version": 1, "settings": SRS, "status": status, "rows": rows})
    SR.write_state(d, {"status": status, "record_sha256": SR._sha256(d / SR.RECORD)})
    J.write_axismach_series(d)
    return d


def same_as_judge_rows(qq: dict, rows: list) -> tuple:
    """criterion stop_rule の各量 = 停止規則の judge_rows (記録の行そのもの; 出口コア M は relative_keys) の値・合否。戻り (一致, 内訳)。"""
    j = SR.judge_rows(rows, **TH, relative_keys=("exit_core_M",))
    j0 = SR.judge_rows(rows, **TH)                       # 停止規則そのものの条件 (r・θ_r)
    rq = qq["mass_imbalance_r"]
    ok = (rq["stop_rule_judge"]["r_mean"] == j["r_mean"] and rq["stop_rule_judge"]["r_amp"] == j["r_amp"]
          and (rq["status"] == "PASS") == j["r_ok"])
    for k in TKEYS:
        ok &= (qq[k]["stop_rule_judge"]["D"] == j["theta"][k]["D"] and qq[k]["stop_rule_judge"]["amp"] == j["theta"][k]["amp"]
               and (qq[k]["status"] == "PASS") == j["theta"][k]["ok"])
    m = j["relative"]["exit_core_M"]
    ok &= (qq["exit_core_M"]["stop_rule_judge"]["D"] == m["D"] and qq["exit_core_M"]["stop_rule_judge"]["amp"] == m["amp"]
           and (qq["exit_core_M"]["status"] == "PASS") == m["ok"])
    ok &= all(qq[k]["status"] == "PASS" for k in ["mass_imbalance_r"] + TKEYS) == j0["ok"]
    return ok, {"judge_rows_ok": j0["ok"], "why": j0.get("why"), "M": m}


rows_ok = sr_rows(25)                                    # 22,500 で最初に条件を満たし、保持 40,000 の後の 62,500 (最後の行) で REACHED
fs = SR.first_stop(rows_ok, hold_steps=int(SRS["hold_steps"]), **TH)
d = sr_run("sr_reached", rows_ok, "REACHED")
meta = json.loads((d / J.SERIES_META).read_text())
c = qs(d)
ok, det = same_as_judge_rows(c["quantities"], rows_ok)
check("到達: 停止規則の first_stop (保持 40,000) が最後の行で REACHED、系列は stop_rule.json の行から、criterion stop_rule の 5 量は全 PASS",
      fs and fs["status"] == "REACHED" and fs["k"] == len(rows_ok) - 1 and meta["source"]["kind"] == "stop_rule.json" and c["status"] == "PASS",
      json.dumps({"first_stop": fs and fs["status"], "reason": c["reason"]}, ensure_ascii=False))
check("到達: 各量の値・合否が judge_rows (記録の行そのもの) と同じ (r・θ_r × 3・出口コア M = relative_keys)", ok, json.dumps(det, ensure_ascii=False))
check("到達: 停止規則の到達と照合して一致 (stop_rule_cross agree)、閾値の出所は stop_rule.json の settings",
      c["stop_rule_cross"]["agree"] is True and c["stop_rule_settings"]["settings"] == SR.thresholds(**TH)
      and c["stop_rule_settings"]["sha256"] == SR._sha256(d / SR.RECORD), json.dumps(c.get("stop_rule_cross"), ensure_ascii=False))
rows_bad = sr_rows(25, drift=2e-4, m_osc=1e-2)           # θ_r が 1 出力あたり 2e-4 ずつ動く・出口コア M が ±1 % で振れる
fs = SR.first_stop(rows_bad, hold_steps=int(SRS["hold_steps"]), **TH)
d = sr_run("sr_censored", rows_bad, "CENSORED")
c = qs(d)
ok, det = same_as_judge_rows(c["quantities"], rows_bad)
qq = c["quantities"]
check("未到達: first_stop は届かない、θ_r は DRIFTING・出口コア M は OSCILLATING・r は STEADY (準定常 FAIL)",
      fs is None and c["status"] == "FAIL" and all(qq[k]["verdict"] == "DRIFTING" for k in TKEYS)
      and qq["exit_core_M"]["verdict"] == "OSCILLATING" and qq["mass_imbalance_r"]["verdict"] == "STEADY",
      json.dumps({k: q.get("verdict") for k, q in qq.items()}))
check("未到達: 各量の値・合否が judge_rows と同じ", ok, json.dumps(det, ensure_ascii=False, default=str))
check("未到達: 停止規則 CENSORED と criterion stop_rule の FAIL は一致 (stop_rule_cross agree)", c["stop_rule_cross"]["agree"] is True,
      json.dumps(c.get("stop_rule_cross"), ensure_ascii=False))
(d / SR.RECORD).write_text(json.dumps({"version": 1, "status": "CENSORED", "rows": rows_bad}))
c = qs(d)
check("stop_rule.json に settings が無い → criterion stop_rule の量はすべて判定不能 (閾値を決めない)",
      c["status"] == "INDETERMINATE" and all(q["status"] == "INDETERMINATE" and "settings" in q["reason"] for q in c["quantities"].values()),
      c["reason"])
(d / SR.RECORD).unlink()
c = qs(d)
check("stop_rule.json が無い → 判定不能", c["status"] == "INDETERMINATE"
      and all("stop_rule.json が無い" in q["reason"] for q in c["quantities"].values()), c["reason"])

print()
if FAIL:
    print(f"{FAIL} FAIL (出力を残した: {TMP})")
    sys.exit(1)
shutil.rmtree(TMP, ignore_errors=True)
print("ALL PASS")
