#!/usr/bin/env python3
"""判定区間の段の接続試験 (GPU・forge が要る。plan tooling-sern-te-wake-grid §5.1 R2、codex result レビュー 2026-10-11 M2)。

生産 YAML (case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml) から、格子を粗くし (約 2.4 万節点) 段の step 数だけを
小さくした試験用の YAML を作り、**標準の実行経路** `python -m forge_design.evaluate.runner_sern3d` (prepare → 段階起動 →
判定区間の段 → 未定常なら同一設定で 1 回延長) を 1 作動点で通す。物理の検証ではなく配線の確認 (段が短いので窓・定常性の
合否は問わない)。生産 YAML の値 (judgment_steps 20000・窓 10000) は変えない (試験の中で確かめる)。

合格 (測る前に固定):
  (a) runner が終わり JUDGMENT.json を書く。段の数 = 1 + 延長の有無 (延長は 1 回まで)
  (b) 各判定区間で床事象の判定が「段階起動の最後の session に窓を当てた判定不能」にならない: verdict は PASS か FAIL、
      区間 = (N − W, N] = (40, 80]、検査した更新数 = W = 40 (レビューの再現の解消)。floor_events.csv の最後の session の長さ = N
  (c) 判定区間の config は本段と nStepOuter だけが違う。全場のスナップショットは step 0 を別枠に、他は TE_MONITOR.csv に 1 行ずつ
      (step 2..80 の 40 行) 記録され、最後の 1 枚 (res_80.h5) だけが残る
  (d) metrics.json の gates に判定区間 (#2・#3) が入る: 力の窓の標本 20/20、温度の窓の標本 20、監視の領域 4 つ。
      延長したら 1 回目の記録が *_judgment1.* に残る
  (e) 延長の枝 (新しく prepare・段階起動した run で、完走した 1 回目に限り延長の判断だけを強制): 2 回目は 1 回目の最終場から継続し、1 回目の記録を退避、
      floor_events.csv に判定区間の session が 2 つ足され、床事象の判定は区間 (N − W, N] (2 回目は完走したとき)

  FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1 \\
    /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_judgment_connection_test.py [--out DIR] [--keep]
"""
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

import check_floor_events as CFE  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

PROD = ROOT / "case" / "46.sern_design" / "problem_3d_prod_3op_wallres_lswx08.yaml"
COARSE3D = {"ni_up": 6, "ni_noz": 30, "ni_plume": 40, "nj_top": 21, "nj_bot": 15, "nz_in": 4, "nz_out": 3, "nj_vside": 5,
            "first_wall_frac": 2.0e-3, "first_z_frac": 0.02, "first_wake_frac": 0.004}
COARSE2 = {"nj_ext_top": 9, "first_top_frac": 0.004}
N, W, OI = 80, 40, 2
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


ap = argparse.ArgumentParser()
ap.add_argument("--out", default=None, help="試験の作業ディレクトリ (既定: 一時ディレクトリ)")
ap.add_argument("--keep", action="store_true", help="作業ディレクトリを残す")
ap.add_argument("--op", default="m10_on")
a = ap.parse_args()
work = Path(a.out) if a.out else Path(tempfile.mkdtemp(prefix="sern_judgment_conn_"))
work.mkdir(parents=True, exist_ok=True)

jo_prod = R2.judgment_options(load_problem(PROD).evaluate)
check("生産 YAML の値は変えない (judgment_steps 20000・延長 1 回・窓 10000)", jo_prod == {"steps": 20000, "extend_once": True, "window_steps": 10000})

y = yaml.safe_load(PROD.read_text())
y["name"] = y["name"] + "_judgment_conn_test"
y["mesh3d"].update(COARSE3D); y["mesh"].update(COARSE2)
y["evaluate"].update({"nStepOuter": 20, "outStepInterval": OI, "floor_events_window_steps": W, "judgment_steps": N})
y["opt"].update({"warm_lam_steps": 20, "soft_steps": 20, "mid_steps": 20})
prob = work / "problem_conn.yaml"
prob.write_text(yaml.safe_dump(y, allow_unicode=True, sort_keys=False))
rd = work / "run_conn"
if rd.exists():
    shutil.rmtree(rd)
env = dict(os.environ)
env.setdefault("FORGE_CUDA_BLOCKSIZE", "128")          # node SLAU のレジスタ上限 (検証の手順と同じ)
env.setdefault("FORGE_ALLOW_UNVERIFIED_SPECIES", "1")  # 検証の手順と同じ (領域別の一様な初期場は化学種の属性を持たない)
t0 = time.time()
r = subprocess.run([sys.executable, "-m", "forge_design.evaluate.runner_sern3d", str(prob), str(rd), "--op", a.op, "--poll-s", "1"],
                   cwd=HERE.parent, env=env, capture_output=True, text=True, timeout=3600)
(work / "runner_stdout.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
print(f"runner rc={r.returncode} ({time.time() - t0:.0f} s)、ログ {work / 'runner_stdout.log'}", flush=True)
print("\n".join(l for l in r.stdout.splitlines() if l.startswith(("JUDGMENT", "GATES"))), flush=True)

jf = rd / R2.JUDGMENT_FILE
check("(a) runner が JUDGMENT.json を書く", jf.exists(), r.stderr[-800:])
if jf.exists():
    j = json.loads(jf.read_text())
    att = [x for x in j["attempts"] if "gates" in x]
    print(f"JUDGMENT {j['verdict']}: {j['reason']}")
    for x in att:
        print(f"  attempt {x['attempt']}: rc {x['rc']} gates {x['gates']}/{x['fail_class']} floor {x['floor_events'] and x['floor_events']['verdict']} "
              f"window {x['force_window']['verdict']} temp {x['temperature']['verdict']} snapshots {x['snapshots']}")
    check("(a) 判定区間の段が 1 回か 2 回 (延長は 1 回まで)、延長は未定常のときだけ",
          len(att) in (1, 2) and (len(att) == 1 or att[0]["fail_class"] in ("JUDGMENT_UNSETTLED", "UNSTEADY")), [(x["attempt"], x["fail_class"]) for x in att])
    for x in att:
        fe = x["floor_events"] or {}
        check(f"(b) {x['attempt']} 回目: 床事象の判定は PASS か FAIL (段階起動の session に窓を当てた判定不能にならない)・区間 ({N - W}, {N}]・更新 {W}",
              fe.get("verdict") in ("PASS", "FAIL") and fe.get("window") == [N - W, N] and fe.get("n_updates_checked") == W, fe)
        check(f"(d) {x['attempt']} 回目: 力の窓の標本 20/20・温度の窓の標本 20・領域 4 つ",
              x["force_window"]["n_tail"] == 20 and x["force_window"]["n_prev"] == 20 and x["temperature"]["n_window"] == 20
              and sorted(x["temperature"]["regions"]) == sorted(("R_TE", "R_SE", "RET", "ALL")), (x["force_window"], x["temperature"]["n_window"]))
    rows, errs = CFE.load(str(rd / CFE.FILE_NAME))
    ends = [rr for rr in rows if rr["kind"] == "session_end"]
    check("(b) floor_events.csv: 最後の session の長さ = N、段階起動の session (暖機・soft・mid・本段) + 判定区間の session が並ぶ",
          not errs and ends and ends[-1]["step"] == N and len(ends) == 4 + len(att), [(rr["session"], rr["step"]) for rr in ends])
    jr = CFE.judge(str(rd), window_steps=W)
    check("(b) check_floor_events (手作業と同じ道具) でも判定区間の段の最後の session に窓 40 が当たる",
          jr["verdict"] in ("PASS", "FAIL") and jr["window"] == [N - W, N], jr["reasons"])
    cm, cj = (rd / "solverConfig_main.yaml").read_text(), (rd / R2.JUDGMENT_CONFIG).read_text()
    dl = [(p, q) for p, q in zip(cm.splitlines(), cj.splitlines()) if p != q]
    check("(c) 判定区間の config は本段と nStepOuter だけが違う", len(dl) == 1 and "nStepOuter: 20" in dl[0][0] and f"nStepOuter: {N}" in dl[0][1], dl)
    mon = list(csv.DictReader(open(rd / R2.MONITOR_CSV)))
    check(f"(c) TE_MONITOR.csv に step {OI}..{N} の {N // OI} 行、step 0 は別枠",
          [int(float(m["step"])) for m in mon] == list(range(OI, N + 1, OI)) and (rd / R2.MONITOR_STEP0_CSV).exists(), len(mon))
    vols = sorted(f.name for f in rd.glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", f.name))
    check("(c) 全場のスナップショットは最後の 1 枚 (res_80.h5) だけ残る", vols == [f"res_{N}.h5"], vols)
    m = json.loads((rd / "metrics.json").read_text())
    check("(d) metrics.json の gates に判定区間 (§6 #2・#3) が入り、JUDGMENT.json の最後の段と同じ判定",
          "judgment" in m["gates"] and m["gates"]["fail_class"] == att[-1]["fail_class"] and m["gates"]["verdict"] == att[-1]["gates"])
    if len(att) == 2:
        check("(d) 延長したら 1 回目の記録を *_judgment1.* に残す",
              all((rd / f).exists() for f in ("metrics_judgment1.json", "force_history_judgment1.csv", "TE_MONITOR_judgment1.csv",
                                               "residual_judgment1.csv", "residual_judgment2.csv")))
    print(f"run: {rd}")

    # (e) 延長の配線 (GPU): 短い段は窓より先に残差・定常性で落ちて延長の枝を通らないことがあるので、**新しく prepare・段階起動した
    #     run** で、1 回目が完走 (rc 0) したときに限り**延長の判断だけ**を強制する (判定そのものは変えない)。実際の forge で
    #     1 回目の最終場からの継続・記録の退避・2 つ目の session・最後の判定区間の metrics を確かめる
    from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
    rd2 = work / "run_conn_ext"
    if rd2.exists():
        shutil.rmtree(rd2)
    os.environ.update({k: env[k] for k in ("FORGE_CUDA_BLOCKSIZE", "FORGE_ALLOW_UNVERIFIED_SPECIES")})
    o = y["opt"]
    R3.prepare(prob, rd2, op=a.op)
    rc0 = R3.run_staged(rd2, "full", int(o["soft_steps"]), soft_cfl=float(o["soft_cfl"]), soft_conv=int(o.get("soft_conv", 0)),
                        warm_lam_steps=int(o["warm_lam_steps"]), warm_lam_cfl=float(o["warm_lam_cfl"]), mid_steps=int(o["mid_steps"]))
    n_beg0 = len([rr for rr in CFE.load(str(rd2 / CFE.FILE_NAME))[0] if rr["kind"] == "session_begin"])
    real_ext = R2.judgment_extendable
    forced = []

    def _force_once(g):
        if not forced and g.get("rc") == 0:
            forced.append(g.get("fail_class"))
            return True, "試験: 延長の配線を通すため 1 回目 (完走) だけ強制"
        return real_ext(g)
    R2.judgment_extendable = _force_once
    try:
        R2.run_judgment(prob, rd2, R3.collect, rc0, monitor=R3.te_monitor_snapshot, precheck=R3.te_monitor_unavailable, poll_s=1.0)
    finally:
        R2.judgment_extendable = real_ext
    j2 = json.loads((rd2 / R2.JUDGMENT_FILE).read_text())
    a2 = [x for x in j2["attempts"] if "gates" in x]
    print(f"(e) 段階起動 rc {rc0}、JUDGMENT {j2['verdict']}: {j2['reason']}、強制した延長: {forced}")
    for x in a2:
        print(f"  attempt {x['attempt']}: restart {x['restart_from']} rc {x['rc']} gates {x['gates']}/{x['fail_class']} "
              f"floor {(x['floor_events'] or {}).get('verdict')} {(x['floor_events'] or {}).get('window')} snapshots {x['snapshots']}")
    check("(e) 1 回目が完走し、延長の判断を 1 回だけ強制した", rc0 == 0 and len(forced) == 1 and a2 and a2[0]["rc"] == 0, (rc0, forced))
    check("(e) 延長の枝: 段が 2 回、1 回目は本段の最終場 (res_20.h5)・2 回目は 1 回目の最終場 (res_80.h5) から継続",
          [x.get("restart_from") for x in j2["attempts"]] == ["res_20.h5", f"res_{N}.h5"], [(x.get("attempt"), x.get("restart_from")) for x in j2["attempts"]])
    check("(e) 1 回目の記録を *_judgment1.* に退避 (監視 40 行)、2 回目は新しい監視の CSV と metrics.json (最後の判定区間が正本)",
          all((rd2 / f).exists() for f in ("metrics_judgment1.json", "force_history_judgment1.csv", "TE_MONITOR_judgment1.csv",
                                            "TE_MONITOR_step0_judgment1.csv", "residual_judgment1.csv", "residual_judgment2.csv"))
          and len(list(csv.DictReader(open(rd2 / "TE_MONITOR_judgment1.csv")))) == N // OI
          and (rd2 / R2.MONITOR_CSV).exists()
          and json.loads((rd2 / "metrics.json").read_text())["gates"]["fail_class"] == a2[-1]["fail_class"])
    rows2 = CFE.load(str(rd2 / CFE.FILE_NAME))[0]
    begins2 = [rr["session"] for rr in rows2 if rr["kind"] == "session_begin"]
    ends2 = {rr["session"]: rr["step"] for rr in rows2 if rr["kind"] == "session_end"}
    check("(e) floor_events.csv に判定区間の 2 つの session が足され (session_begin)、1 回目の session は N step で閉じ、"
          "1 回目の床事象の判定は区間 (40, 80]",
          len(begins2) == n_beg0 + 2 and ends2.get(begins2[-2]) == N and (a2[0]["floor_events"] or {}).get("window") == [N - W, N],
          (begins2, ends2))
    if a2 and a2[-1]["rc"] == 0:
        check("(e) 2 回目も完走したら、2 回目の床事象の判定も区間 (40, 80]・監視 40 行",
              ends2.get(begins2[-1]) == N and (a2[-1]["floor_events"] or {}).get("window") == [N - W, N]
              and len(list(csv.DictReader(open(rd2 / R2.MONITOR_CSV)))) == N // OI)
    else:
        # 粗い格子・短い起動の流れが続きで発散した (配線とは別)。発散した段は合格にしない: 床事象の判定は PASS にならず、
        # 総合も FAIL (延長の後の段も延長しない)
        print("  (2 回目は完走しなかった: 粗い格子・短い起動の流れの発散。延長の配線の確認には影響しない)")
        check("(e) 2 回目が発散したら、その段の床事象の判定は PASS にならず (session が閉じていない)、JUDGMENT は FAIL",
              (a2[-1]["floor_events"] or {}).get("verdict") != "PASS" and j2["verdict"] == "FAIL" and a2[-1]["fail_class"] == "DIVERGED")
    print(f"run (延長の配線): {rd2}")
if not a.keep and not a.out:
    shutil.rmtree(work, ignore_errors=True)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
