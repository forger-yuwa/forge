#!/usr/bin/env python3
"""`judge_run.py` の停止規則 v2 対応の試験 (plan tooling-design-orchestrator §5.1 #5b・§6 J、judge plan §4.6 の 8;
2026-10-18 codex 諮問 judge-v2)。GPU・forge なし。

case/71 の run_0014_chain_v2・run_0011_ns_v2 (手元) の写し (h5 は symlink、他は複製) に、判定の仕様 `JUDGE_SPEC.json` を
#5b の生成器の関数 (`judge.spec_axismach_settle`) で**後から**書いて `judge_run.py` を回す (移行試験。生成器が起動の前に仕様を書く
経路は (G1) で別に見る)。

  §6 J (事前登録した判別の A/B):
  (J-A) 変更なし → convergence PASS (SETTLED)・stop_rule PASS (REACHED)・quasisteady PASS
  (J-B) boundary_massflux.csv の先頭のコメント行の 1 文字を変える → check_settled の数値と SETTLED は A と同じ、停止の時の再評価との照合が
        通らず convergence INDETERMINATE (総合は PASS にならない)。A が PASS でなければ試験の不成立
  そのほか:
  (R11) run_0011 の写しも convergence・stop_rule・quasisteady が PASS
  (N1)  停止規則 v2 の状態 FAILED → stop_rule FAIL (途中の記録なら判定不能)
  (N2)  仕様から量を 1 つ除く (余り) → quasisteady 判定不能
  (N3)  仕様に基準に無い量を足す (欠け) → quasisteady 判定不能
  (N4)  仕様の基準の sha256 が違う → quasisteady 判定不能
  (N5)  仕様の量の種類が基準と違う → quasisteady 判定不能
  (N6)  STOP_RULE_STATE.json に再評価の記録が無い → convergence 判定不能 (照合から外さない)
  (N7)  流量の数値を変えて判定し直しが NOT_SETTLED → convergence FAIL のまま、照合の不能を併記 (既知の FAIL を隠さない)
  (F1)  鮮度: 判定の直後は FRESH、流量の CSV のコメントを変えると STALE
  (F2)  鮮度: 判定のときに無かった固定の入力 (settle_nodes.txt) が後から現れると STALE
  (G1)  生成器の経路: run_staged_ns (ダミーの forge) が最初の起動の前に criterion settle の仕様と基準の sha256 を書く

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_judge_v2_tests.py
所要: 数分 (judge_run を 12 回)。run が無ければ飛ばす。
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(REPO / "solver_density_cuda" / "tools"))

from forge_design.evaluate import judge as J  # noqa: E402
from forge_design.evaluate import runner as RN  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
import judge_run as JR  # noqa: E402

C71 = REPO / "case" / "71.m55_combustion_wt"
R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", C71 / "run_0014_chain_v2"))
R11 = Path(os.environ.get("FORGE_CASE71_RUN_0011", C71 / "run_0011_ns_v2"))
JUDGE = REPO / "solver_density_cuda" / "tools" / "judge_run.py"
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""), flush=True)
    if not cond:
        FAIL += 1


if not (R14 / "nozzle.h5").exists() or not (R14 / "settle_criteria.yaml").is_file():
    print(f"飛ばした: {R14} に run が無い")
    sys.exit(0)

td = Path(tempfile.mkdtemp(prefix="judge_v2_tests_"))


def copy_run(src: Path, tag: str) -> Path:
    """run の写し: h5 は symlink、ディレクトリ (stop_rule_v2 など) とテキストは複製 (mtime も)。"""
    d = td / tag
    d.mkdir()
    for f in src.iterdir():
        if f.is_dir():
            shutil.copytree(f, d / f.name, symlinks=True)
        elif f.suffix == ".h5":
            (d / f.name).symlink_to(f.resolve())
        else:
            shutil.copy2(f, d / f.name)
    return d


def write_spec(d: Path, mutate=None) -> dict:
    spec = J.spec_axismach_settle("試験 (後付けの仕様; 移行試験 run_judge_v2_tests)", d)
    if mutate:
        mutate(spec)
    (d / J.SPEC).write_text(json.dumps(spec, ensure_ascii=False, indent=1) + "\n")
    return spec


def judge(d: Path) -> dict:
    r = subprocess.run([sys.executable, str(JUDGE), str(d)], capture_output=True, text=True, timeout=1800)
    v = json.loads((d / JR.VERDICT_JSON).read_text())
    v["_rc"] = r.returncode
    return v


def st(v, k):
    return (v.get("checks") or {}).get(k, {}).get("status")


def statuses(v):
    return {k: st(v, k) for k in JR.CHECK_ORDER}


def edit_comment(d: Path):
    """boundary_massflux.csv の先頭のコメント行の 1 文字を変える (数値の行は変えない)。"""
    p = d / "boundary_massflux.csv"
    lines = p.read_text().split("\n")
    assert lines[0].startswith("# boundary_massflux.csv")
    lines[0] = lines[0].replace("# boundary_massflux.csv", "# Boundary_massflux.csv", 1)
    p.write_text("\n".join(lines))


try:
    print("== §6 J: 判別の A/B (run_0014 の写し)")
    a = copy_run(R14, "J_A")
    write_spec(a)
    va = judge(a)
    ca = va["checks"]["convergence"]
    check("(J-A) 変更なし → convergence PASS (SETTLED)・stop_rule PASS (REACHED)・quasisteady PASS",
          st(va, "convergence") == "PASS" and ca.get("verdict") == "SETTLED" and ca.get("tool") == "check_settled.py"
          and st(va, "stop_rule") == "PASS" and st(va, "quasisteady") == "PASS" and (ca.get("stop_match") or {}).get("pass") is True,
          f"{statuses(va)} {ca.get('reason', '')[:120]}")
    b = copy_run(R14, "J_B")
    write_spec(b)
    edit_comment(b)
    vb = judge(b)
    cb = vb["checks"]["convergence"]
    probs = (cb.get("stop_match") or {}).get("problems") or []
    check("(J-B) コメントの 1 文字 → 数値と SETTLED は A と同じ、照合が通らず convergence INDETERMINATE (総合は PASS でない)",
          cb.get("verdict") == "SETTLED" and cb.get("by_quantity") == ca.get("by_quantity") and st(vb, "convergence") == "INDETERMINATE"
          and vb["overall"] != "PASS" and any("boundary_massflux.csv" in p for p in probs),
          f"{statuses(vb)} {cb.get('reason', '')[:200]}")
    check("(J) 判別に使わない項目 (solver・nonfinite・physical・wall_resolution) は A と B で同じ",
          all(st(va, k) == st(vb, k) for k in ("solver", "nonfinite", "physical", "wall_resolution")),
          f"A {statuses(va)} / B {statuses(vb)}")

    if (R11 / "nozzle.h5").exists() and (R11 / "settle_criteria.yaml").is_file():
        print("== run_0011 の写し")
        r11 = copy_run(R11, "R11")
        write_spec(r11)
        v11 = judge(r11)
        check("(R11) run_0011 の写し: convergence・stop_rule・quasisteady が PASS",
              all(st(v11, k) == "PASS" for k in ("convergence", "stop_rule", "quasisteady")), f"{statuses(v11)}")
    else:
        print("飛ばした: (R11) run_0011 が無い")

    print("== 停止規則の状態・仕様の照合・再評価の記録")
    c1 = JR.check_stop_rule({"raw": "FAILED", "final": True, "source": "STOP_RULE_STATE.json", "record": {"method": "v2"}})
    c2 = JR.check_stop_rule({"raw": "FAILED", "final": False, "source": "stop_rule.json", "record": {"status": "FAILED"}})
    check("(N1) 停止規則 v2 の FAILED → FAIL (途中の記録なら判定不能)", c1["status"] == "FAIL" and c2["status"] == "INDETERMINATE",
          f"{c1['status']} / {c2['status']}")
    cases = {
        "N2": ("仕様から量を 1 つ除く (余り)", lambda s: s["quasisteady"]["quantities"].pop(), "仕様に無い"),
        "N3": ("仕様に基準に無い量を足す (欠け)", lambda s: s["quasisteady"]["quantities"].append(
            {**s["quasisteady"]["quantities"][0], "name": "bogus", "extractor": {"kind": "settle", "quantity": "bogus", "settle_kind": "residual"}}),
            "無い"),
        "N4": ("仕様の基準の sha256 が違う", lambda s: s["quasisteady"]["settle_criteria"].update(sha256="0" * 64), "sha256"),
        "N5": ("仕様の量の種類が基準と違う", lambda s: s["quasisteady"]["quantities"][1]["extractor"].update(settle_kind="residual"), "種類"),
    }
    for key, (what, mut, word) in cases.items():
        d = copy_run(R14, key)
        write_spec(d, mut)
        v = judge(d)
        q = v["checks"]["quasisteady"]
        txt = q.get("reason", "") + json.dumps(q.get("quantities", {}), ensure_ascii=False)
        check(f"({key}) {what} → quasisteady 判定不能、convergence は PASS のまま (数値の判定は同じ)",
              st(v, "quasisteady") == "INDETERMINATE" and st(v, "convergence") == "PASS" and word in txt and v["overall"] != "PASS",
              f"{statuses(v)} {q.get('reason', '')[:160]}")
    d = copy_run(R14, "N6")
    write_spec(d)
    s_ = json.loads((d / "STOP_RULE_STATE.json").read_text())
    s_.pop("reeval", None)
    (d / "STOP_RULE_STATE.json").write_text(json.dumps(s_, ensure_ascii=False, indent=1))
    v = judge(d)
    check("(N6) STOP_RULE_STATE.json に再評価の記録が無い → convergence 判定不能 (照合から外さない)",
          st(v, "convergence") == "INDETERMINATE" and "再評価" in v["checks"]["convergence"].get("reason", ""),
          f"{statuses(v)} {v['checks']['convergence'].get('reason', '')[:160]}")
    d = copy_run(R14, "N7")
    write_spec(d)
    p = d / "boundary_massflux.csv"
    lines = p.read_text().split("\n")
    hdr_i = next(i for i, ln in enumerate(lines) if ln.startswith("step,"))
    cols = lines[hdr_i].split(",")
    jo = cols.index("mdot_outlet_2")
    for i in range(hdr_i + 1, len(lines)):
        if lines[i].strip():
            f = lines[i].split(",")
            if int(f[0]) >= 30000:
                f[jo] = repr(float(f[jo]) * 1.01)       # 質量の収支を 1 % 崩す (許容 |平均| 5e-4)
            lines[i] = ",".join(f)
    p.write_text("\n".join(lines))
    v = judge(d)
    cv = v["checks"]["convergence"]
    check("(N7) 判定し直しが NOT_SETTLED (流量を変えた) → convergence FAIL のまま、照合の不能を併記",
          st(v, "convergence") == "FAIL" and cv.get("verdict") == "NOT_SETTLED" and "照合も通らない" in cv.get("reason", "")
          and v["overall"] == "FAIL", f"{statuses(v)} {cv.get('reason', '')[:200]}")

    print("== 鮮度")
    fr = subprocess.run([sys.executable, str(JUDGE), str(a), "--check-fresh"], capture_output=True, text=True)
    edit_comment(a)
    fr2 = subprocess.run([sys.executable, str(JUDGE), str(a), "--check-fresh"], capture_output=True, text=True)
    check("(F1) 判定の直後は FRESH、流量の CSV のコメントを変えると STALE",
          fr.returncode == 0 and fr2.returncode == 1 and "boundary_massflux.csv" in fr2.stdout, f"{fr.stdout.strip()[:80]} / {fr2.stdout.strip()[:200]}")
    d = copy_run(R14, "F2")
    write_spec(d)
    nodes = (d / "settle_nodes.txt").read_bytes()
    (d / "settle_nodes.txt").unlink()
    v = judge(d)
    (d / "settle_nodes.txt").write_bytes(nodes)
    fr3 = subprocess.run([sys.executable, str(JUDGE), str(d), "--check-fresh"], capture_output=True, text=True)
    check("(F2) 判定のときに無かった固定の入力 (settle_nodes.txt) が後から現れると STALE",
          fr3.returncode == 1 and "settle_nodes.txt" in fr3.stdout, f"判定 {v['overall']}; {fr3.stdout.strip()[:200]}")

    print("== 生成器の経路")
    dummy = td / "dummy_forge.sh"
    dummy.write_text("#!/bin/bash\nexit 0\n")
    dummy.chmod(0o755)
    RN._ENV["FORGE_BIN"] = str(dummy)
    os.environ["FORGE_BIN"] = str(dummy)
    g = td / "G1"
    g.mkdir()
    skip = {"STOP_REQUEST", "STOPPED_AT", "SNAPSHOT_REQUEST", "MAIN_STAGE_STARTED", "STOP_RULE_STATE.json", "execution_plan.json", J.SPEC}
    for f in R14.iterdir():
        if f.name in skip or (f.name.startswith("res_") and f.suffix in (".h5", ".xmf")) or f.is_dir():
            continue
        if f.suffix == ".h5":
            (g / f.name).symlink_to(f.resolve())
        else:
            shutil.copy2(f, g / f.name)
    try:
        RA.run_staged_ns(g)
    except (RuntimeError, ValueError):
        pass                                     # ダミーの forge は場を書かないので最初の段で止まる (試験の対象は仕様を書くところまで)
    sp = json.loads((g / J.SPEC).read_text()) if (g / J.SPEC).is_file() else {}
    qs = (sp.get("quasisteady") or {})
    crit_sha = hashlib.sha256((g / "settle_criteria.yaml").read_bytes()).hexdigest()
    names = [q["name"] for q in qs.get("quantities") or []]
    check("(G1) run_staged_ns は最初の起動の前に criterion settle の仕様 (基準の量の全部・基準の sha256) を書く",
          names == ["residual", "mass_balance", "exit_core_mach", "line_mach", "exit_flow_angle", "wall_heat_nodes"]
          and all(q["criterion"] == "settle" for q in qs.get("quantities") or [])
          and (qs.get("settle_criteria") or {}).get("sha256") == crit_sha, f"{names} {str((qs.get('settle_criteria') or {}))[:120]}")
finally:
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
