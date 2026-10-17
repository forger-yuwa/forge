#!/usr/bin/env python3
"""`wt_qualify` の解析の成立に `judge_run` の回し直しを入れた試験 (plan tooling-design-orchestrator §5.1 #5c、§4.2「解析の成立の読み方」;
2026-10-18 codex 諮問 judge-v2 の (6))。GPU・forge なし。

  (D) decide() の分岐 (合成の判定):
      D1 回し直し PASS・すべて合格・成果物 PASS → accepted
      D2 回し直しで壁解像だけ FAIL (この判定の壁解像も FAIL) → back_to_correction・補正 mesh_first_layer (解析の総合は FAIL のまま)
      D3 回し直しの quasisteady が判定不能 (convergence は PASS) → incomplete
      D4 回し直しの壁解像 FAIL・この判定の壁解像 PASS (食い違い) → incomplete
      D5 古い run (NOT_APPLICABLE) → 今までどおり accepted
  (Q) 本物の run の写し (case/71 run_0014・試走 A case/90 run_0003_ns_p0、h5 は symlink):
      Q1 run_0014 の写し → 回し直し PASS (7 項目)・accepted (正例は accepted のまま)
      Q2 試走 A run_0003_ns_p0 の写し (公差未達でも整定) → 回し直し PASS・back_to_correction (V5b・V5c FAIL)
      Q3 run_0014 の写しで run の中の仕様 JUDGE_SPEC.json の準定常の量が判定の側の仕様と違う → incomplete
  必須の CAD の欠け (cad.pass=false) → incomplete は run_wt_qualify_tests.py の (4c) が見る。

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wt_qualify_rv_tests.py
所要: 数分。run が無ければその項目を飛ばす。
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import judge as J  # noqa: E402
from forge_design.feedback import wt_qualify as WQ  # noqa: E402

C71 = ROOT / "case" / "71.m55_combustion_wt"
R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", C71 / "run_0014_chain_v2"))
CARD14, PROB14 = C71 / "card_m55_3coef.yaml", C71 / "problem_m55_p5_settle_v2.yaml"
QA = ROOT / "case" / "90.qs_trial_a"
RA3 = Path(os.environ.get("FORGE_QS_TRIAL_A_RUN_0003", QA / "run_0003_ns_p0"))
FAIL = 0
SKIPPED = []


def check(name, ok, detail=""):
    global FAIL
    print(("ok  " if ok else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok:
        FAIL += 1


def dec(rv, wr="PASS", dl=None):
    d = WQ.decide({"status": "PASS"}, {"verdict": "PASS"}, {"verdict": "PASS"}, {"verdict": wr},
                  dl or {"status": "PASS", "dimensions": {"status": "PASS"}}, rv=rv)
    return d


print("== (D) decide() の分岐")
rvP = {"status": "PASS", "wall_resolution": "PASS", "why": "7 項目すべて PASS"}
rvW = {"status": "WALL_ONLY_FAIL", "wall_resolution": "FAIL", "why": "壁解像だけ FAIL"}
rvQ = {"status": "INCOMPLETE", "wall_resolution": "PASS", "why": "解析の成立の項目が PASS でない: quasisteady INDETERMINATE (仕様の量の余り)"}
d = dec(rvP)
check("D1 回し直し PASS・すべて合格・成果物 PASS → accepted", (d["result"], d["next"]) == ("accepted", "done"), str(d["missing"]))
d = dec(rvW, wr="FAIL")
check("D2 回し直しで壁解像だけ FAIL (この判定も FAIL) → back_to_correction・mesh_first_layer",
      d["result"] == "back_to_correction" and [c["kind"] for c in d["corrections"]] == ["mesh_first_layer"] and not d["missing"],
      f"{d['result']} {d['corrections']} {d['missing']}")
d = dec(rvQ)
check("D3 回し直しの quasisteady 判定不能 (convergence は PASS) → incomplete",
      d["result"] == "incomplete" and any("quasisteady" in m for m in d["missing"]), f"{d['result']} {d['missing']}")
d = dec(rvW, wr="PASS")
check("D4 回し直しの壁解像 FAIL・この判定の壁解像 PASS (食い違い) → incomplete",
      d["result"] == "incomplete" and any("食い違う" in m for m in d["missing"]), f"{d['result']} {d['missing']}")
d = dec({"status": "NOT_APPLICABLE", "why": "古い run"})
check("D5 古い run (NOT_APPLICABLE) → 今までどおり accepted", d["result"] == "accepted", str(d["missing"]))


def copy_run(src: Path, dst: Path) -> Path:
    """run の写し (h5・xmf は symlink、ディレクトリは複製 (中の h5 も symlink))。"""
    run = dst / src.name
    run.mkdir(parents=True)
    for f in src.iterdir():
        if f.is_dir():
            if f.name in ("deliverables", "wt_qualify", "wt_qualify_tol01", "wt_qualify_final", "handoff_extract", "contur_calibrate"):
                continue
            shutil.copytree(f, run / f.name, symlinks=True, copy_function=lambda s, d: os.symlink(Path(s).resolve(), d)
                            if s.endswith(".h5") else shutil.copy2(s, d))
        elif f.suffix in (".h5", ".xmf"):
            (run / f.name).symlink_to(f.resolve())
        else:
            shutil.copy2(f, run / f.name)
    return run


print("== (Q) 本物の run の写し")
if (R14 / "nozzle.h5").exists() and (R14 / "deliverables").is_dir():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        run = copy_run(R14, td)
        rec = WQ.qualify(run, CARD14, PROB14, deliverables=R14 / "deliverables", wall_txt=td / "wr.txt", series_dir=td, judge_dir=td / "jr")
        rv = rec["checks"]["run_verdict"]
        check("Q1 run_0014 の写し → 回し直し PASS (7 項目)・accepted のまま",
              rv["status"] == "PASS" and all(v == "PASS" for v in rv["items"].values()) and rec["result"] == "accepted"
              and (td / "jr" / "RUN_VERDICT.json").is_file() and not (run / "RUN_VERDICT.json").exists(),
              f"{rv['status']} {rv.get('items')} {rec['result']} {rec['missing']}")
        sp = J.spec_axismach_settle("試験", run)
        sp["quasisteady"]["quantities"].pop()
        (run / J.SPEC).write_text(json.dumps(sp, ensure_ascii=False))
        rec = WQ.qualify(run, CARD14, PROB14, deliverables=R14 / "deliverables", wall_txt=td / "wr2.txt", series_dir=td, judge_dir=td / "jr2")
        rv = rec["checks"]["run_verdict"]
        check("Q3 run の中の仕様の準定常の量が判定の側の仕様と違う → 回し直し INCOMPLETE・incomplete",
              rv["status"] == "INCOMPLETE" and rec["result"] == "incomplete" and "判定の側の仕様と違う" in (rv.get("why") or ""),
              f"{rv['status']} {rec['result']} {rv.get('why')}")
else:
    SKIPPED.append(f"Q1・Q3 ({R14} か成果物が無い)")

if (RA3 / "STOP_RULE_STATE.json").is_file() and (QA / "card_p0.yaml").is_file():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        run = copy_run(RA3, td)
        rec = WQ.qualify(run, QA / "card_p0.yaml", QA / "problem_p0.yaml", wall_txt=td / "wr.txt", series_dir=td, judge_dir=td / "jr")
        rv = rec["checks"]["run_verdict"]
        c = rec["checks"]
        check("Q2 試走 A run_0003_ns_p0 の写し (公差未達でも整定) → 回し直し PASS・back_to_correction (V5b・V5c FAIL)",
              rv["status"] == "PASS" and rec["result"] == "back_to_correction" and c["V5b"]["verdict"] == "FAIL"
              and c["V5c"]["verdict"] == "FAIL", f"{rv['status']} {rv.get('items')} {rec['result']} V5b {c['V5b'].get('verdict')} "
                                                   f"V5c {c['V5c'].get('verdict')} {rec['missing']}")
else:
    SKIPPED.append(f"Q2 ({RA3} が無い)")

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
