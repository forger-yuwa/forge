#!/usr/bin/env python3
"""正式な入口の回し直し (`wt_qualify.approval_rechecks`) の試験 (plan tooling-design-orchestrator §4.7.1 (6)、§5.1 #2c)。
run の側で作られた判定 (停止の判定・計算の格子の品質) を、判定の側のコードで回し直す。手元の case/71 run_0014_chain_v2 を読む (書かない)。

  (P) 停止の step 46171・格子の AR ≤ 5000 → 回し直しに問題なし (停止の判定 SETTLED・格子の品質 PASS)
  (N1) 停止の step を 20000 にした流れの記録 → 停止の判定を回し直すと SETTLED でない
  (N2) 承認の格子の AR の上限 10 → 格子の品質が PASS でない
  (N3) 停止の step が無い → 回し直せない
  (R14) run の物理の照合 (plan §5.1 #3) の 2 例目: case/71 の承認と run_0014 の本段の設定が一致、停止規則の基準 (wt_nozzle_v2) も一致

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_approval_rechecks_tests.py
所要: 1〜2 分 (check_settled を 2 回)。run が無ければ飛ばす (rc 0 と「飛ばした」)。
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.feedback import wt_qualify as Q  # noqa: E402

R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", Path(__file__).resolve().parents[2] / "case" / "71.m55_combustion_wt" / "run_0014_chain_v2"))
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


if not (R14 / "SETTLE_VERDICT.json").is_file() or not (R14 / "nozzle.h5").is_file():
    print(f"飛ばした: {R14} に run が無い")
    sys.exit(0)
APR = {"protected": {"evaluation": {"mesh_quality": {"ar_max": 5000.0, "skew_max": 0.9}}}}
r = Q.approval_rechecks(R14, {"stop_step": 46171}, APR)
check("(P) 停止の step 46171・AR ≤ 5000 → 問題なし", not r["problems"], f"{r.get('settle', {}).get('verdict_line')} / {r.get('mesh_quality', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {"stop_step": 20000}, APR)
check("(N1) 停止の step 20000 → 停止の判定が SETTLED でない", any("SETTLED でない" in p for p in r["problems"]), f"{r.get('settle', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {"stop_step": 46171}, {"protected": {"evaluation": {"mesh_quality": {"ar_max": 10.0, "skew_max": 0.9}}}})
check("(N2) AR の上限 10 → 格子の品質が PASS でない", any("品質" in p for p in r["problems"]), f"{r.get('mesh_quality', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {}, APR)
check("(N3) 停止の step が無い → 回し直せない", any("停止の step が無い" in p for p in r["problems"]))

# run の物理の照合 (plan §5.1 #3) の正例: case/71 の承認 (card_m55_3coef を版 2 にしたもの + run_0014 を作った問題 YAML) と run_0014 の本段の設定
import json  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import yaml  # noqa: E402

from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

C71 = R14.parent
CARD, PROB = C71 / "card_m55_3coef.yaml", C71 / "problem_m55_p5_settle_v2.yaml"
if CARD.is_file() and PROB.is_file():
    with tempfile.TemporaryDirectory() as tdn:
        td = Path(tdn)
        repo = td / "repo"
        repo.mkdir()
        for args in (["init", "-q"], ["config", "user.email", "t@t"], ["config", "user.name", "t"]):
            subprocess.run(["git", "-C", str(repo), *args], check=True)
        (repo / "x").write_text("x")
        subprocess.run(["git", "-C", str(repo), "add", "x"], check=True)
        subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "x"], check=True)
        d = yaml.safe_load(CARD.read_text())
        src = {k: {"kind": "user", "ref": "case/71"} for k in ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall",
                                                             "tolerances", "p_ambient", "budget")}
        d.update(format=W.CARD_FORMAT_V2, budget={"gpu_hours": 4.0, "disk_gb": 50.0}, mesh_limits={"max_nodes_ns": 443000, "max_nodes_euler": 467000}, sources={**src, "mesh_limits": {"kind": "user", "ref": "case/71 のユーザ決定「a」2026-10-15 (手組みの格子の節点数)"}})  # §5.1 #8
        card = td / "card.yaml"
        card.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
        pv = AP.preview(card, PROB, repo_dir=repo, judge_env=False)
        (td / "pv.json").write_text(json.dumps(pv, ensure_ascii=False))
        check("(R14-前提) case/71 の承認の preview に問題が無い", not pv["problems"], str(pv["problems"])[:300])
        if not pv["problems"]:
            AP.write(td / "pv.json", "試験", td / "ap.json")
            ph = AP.check_run_physics(AP.load_approval(td / "ap.json"), R14, PROB)
            check("(R14) case/71 run_0014 の本段の設定が承認の物理と一致 (組成・寸法の違う 2 例目)", ph["status"] == "OK", str(ph["problems"])[:300])
            rr = AP.check_run(td / "ap.json", R14, PROB)
            # run_0014 はプリセット wt_nozzle_v2 (残差を判定する版) で回っていて、v2 の中身は今も同じなので、基準も一致する
            # (試走 A の run_0003 は v3 で、v3 に後から足した residual_finite.required が無いので MISMATCH になる; run_approval_tests (N15))
            check("(R14) check_run (停止規則の基準 wt_nozzle_v2 + 物理) が OK", rr["status"] == "OK", str(rr["problems"] + rr["diffs"])[:300])
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
