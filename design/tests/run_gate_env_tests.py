#!/usr/bin/env python3
"""判定専用の環境 (`forge_design.feedback.gate_env`) の試験 (plan tooling-design-orchestrator §4.7.1 (5)、§5.1 #2b)。
正例と負例を並べる (全部止める実装を合格にしない)。

**このリポジトリの HEAD を承認の commit に使う**ので、追跡ファイルに変更の無い状態 (commit した後) で回す。checkout は一時の
ディレクトリに作り (FORGE_GATE_ROOT)、最後に `git worktree remove` で消す。

  正例
  (P1) 承認の記録 (判定の Python・外部の道具つき) の check が OK (sparse checkout を作る)
  (P2) probe: 判定専用の環境で forge_design の import の元が checkout の中、印の検査が通る
  (P3) judge: wt_qualify --approval が判定専用の環境で門を通る (固定の入力には場が無いので、その先の場の読み込みで止まる —
       門の誤りでも checkout の欠け (物性の DB など) でもないこと)
  負例
  (N1) 同じ HEAD・追跡ファイルでも、判定の Python が違う (/usr/bin/python3: prefix・配布物) → check が MISMATCH
  (N2) 印の無い環境 (この試験のプロセス) の assert_in_gate_env → 止める
  (N3) 印はあるが checkout が import の元と違う → 止める
  (N4) checkout のコードの探索先に .py を足す → check が止める (消すと通る)
  (N5) 外部の道具 (変換器の代わりの一時のファイル) の中身が変わる → compare_env が違いを出す
  (N6) wt_qualify --approval を作業ツリーから直接回す → 「判定専用の環境で回していない」で止まる
  (N7) 承認の記録の gate を書き換える → load_approval が止める
  (N8) 作業ツリーの追跡ファイルに変更がある状態の preview → 問題にする (承認の commit を固定できない)

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_gate_env_tests.py
所要: 約 30 秒。
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from forge_design.feedback import gate_env as GE  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

FIX = HERE / "fixtures" / "approval_trial_a"
REPO = HERE.parents[1]
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def raises(fn, needle):
    try:
        fn()
    except ValueError as e:
        return needle in str(e)
    return False


def git(*args, cwd=REPO):
    return subprocess.run(["git", "-C", str(cwd), *args], capture_output=True, text=True)


if git("status", "--porcelain", "--untracked-files=no").stdout.strip():
    print("止めた: 作業ツリーの追跡ファイルに変更がある (この試験は HEAD を承認の commit に使う; commit してから回す)")
    sys.exit(2)

SOURCES = {s: {"kind": "user", "ref": "依頼の原文"} for s in ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall")}
SOURCES.update({"p_ambient": {"kind": "recipe_default"}, "tolerances": {"kind": "recipe_default"},
                "budget": {"kind": "user", "ref": "見積もりの承認"}})

tdn = tempfile.mkdtemp(prefix="gate_env_tests_")
td = Path(tdn)
os.environ[GE.GATE_ROOT_ENV] = str(td / "gates")
head = git("rev-parse", "HEAD").stdout.strip()
co = td / "gates" / head
try:
    d = yaml.safe_load((FIX / "card_p0.yaml").read_text())
    d.update(format=W.CARD_FORMAT_V2, budget={"gpu_hours": 4.0, "disk_gb": 50.0}, sources=SOURCES)
    card = td / "card.yaml"
    card.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    prob = td / "problem.yaml"
    shutil.copy(FIX / "problem_p0.yaml", prob)
    pv = AP.preview(card, prob)               # 既定: このリポジトリの HEAD・今の Python と道具を記録
    check("前提: preview に問題が無い (HEAD がきれい・判定の環境を記録)", not pv["problems"] and pv["gate"]["commit"] == head
          and pv["gate"].get("env") is not None, str(pv["problems"])[:300])
    (td / "preview.json").write_text(json.dumps(pv, ensure_ascii=False))
    apath = td / "card_approval.json"
    AP.write(td / "preview.json", "ユーザ「よい」(試験の仮の承認)", apath)

    print("== 正例")
    r = GE.check(apath)
    check("(P1) check が OK (sparse checkout を作る)", r["status"] == "OK" and Path(r["checkout"]) == co, str(r["problems"])[:300])
    sz = sum(f.stat().st_size for f in co.rglob("*") if f.is_file() and ".git" not in f.parts)
    check("(P1) checkout は design・solver_density_cuda/tools・data だけ (50 MB 未満)",
          sz < 50e6 and (co / "design" / "forge_design").is_dir() and (co / "solver_density_cuda" / "data" / "species").is_dir()
          and not (co / "case").exists(), f"{sz / 1e6:.1f} MB")
    pr = GE.probe(apath)
    check("(P2) probe: import の元が checkout の中、印の検査が通る",
          pr["status"] == "OK" and pr["forge_design"].startswith(str(co)) and pr["assert"]["commit"] == head,
          str({k: pr.get(k) for k in ("status", "forge_design", "problems")})[:300])
    run = td / "run"
    shutil.copytree(FIX / "run_0003_ns_p0", run)
    q = subprocess.run([sys.executable, "-m", "forge_design.feedback.gate_env", "judge", "--approval", str(apath), "--",
                        str(run), "--card", str(card), "--problem", str(prob), "--iteration", "1", "--no-write"],
                       capture_output=True, text=True, cwd=str(HERE.parent))
    gate_words = ("判定専用の環境", "checkout", "import の元", "承認の記録と違う", "CLI の公差")
    check("(P3) judge: 門を通り、その先の場の読み込み (res_*.h5 が無い) で止まる",
          q.returncode != 0 and "止めた" in q.stderr and not any(w in q.stderr for w in gate_words)
          and "No such file" not in q.stderr and ("res_" in q.stderr or "全場" in q.stderr),
          q.stderr.strip()[-300:])

    print("== 負例")
    r = GE.check(apath, python="/usr/bin/python3")
    check("(N1) 判定の Python が違う (/usr/bin/python3) → MISMATCH",
          r["status"] == "MISMATCH" and any("prefix" in p or "配布物" in p for p in r["problems"]), str(r["problems"])[:300])
    apr = AP.load_approval(apath)
    os.environ.pop(GE.MARK_ENV, None)
    check("(N2) 印の無い環境の assert_in_gate_env → 止める", raises(lambda: GE.assert_in_gate_env(apr), "判定専用の環境で回していない"))
    os.environ[GE.MARK_ENV] = "1"
    os.environ[GE.CHECKOUT_ENV] = str(co)
    check("(N3) 印はあるが import の元 (作業ツリー) が checkout の外 → 止める", raises(lambda: GE.assert_in_gate_env(apr), "import の元"))
    os.environ.pop(GE.MARK_ENV, None)
    os.environ.pop(GE.CHECKOUT_ENV, None)
    stray = co / "design" / "forge_design" / "zz_stray.py"
    stray.write_text("X = 1\n")
    check("(N4) checkout のコードの探索先に .py を足す → check が止める", raises(lambda: GE.check(apath), "未追跡"))
    stray.unlink()
    check("(N4) 消すと通る", GE.check(apath)["status"] == "OK")
    tool = td / "fake_converter"
    tool.write_text("a")
    snap = GE.snapshot(tools={"FORGE_CONVERTER": str(tool)})
    tool.write_text("b")
    now = GE.snapshot(tools={"FORGE_CONVERTER": str(tool)})
    probs = GE.compare_env(snap, now)
    check("(N5) 外部の道具の中身が変わる → 違いを出す (同じなら出さない)",
          any("FORGE_CONVERTER" in p for p in probs) and not GE.compare_env(snap, snap), str(probs)[:200])
    q = subprocess.run([sys.executable, "-m", "forge_design.feedback.wt_qualify", str(run), "--card", str(card), "--problem", str(prob),
                        "--iteration", "1", "--no-write", "--approval", str(apath)], capture_output=True, text=True, cwd=str(HERE.parent))
    check("(N6) wt_qualify --approval を作業ツリーから直接回す → 止まる",
          q.returncode == 1 and "判定専用の環境で回していない" in q.stderr, q.stderr.strip()[-200:])
    bad = json.loads(apath.read_text())
    bad["gate"]["commit"] = "0" * 40
    (td / "approval_bad.json").write_text(json.dumps(bad, ensure_ascii=False))
    check("(N7) 承認の記録の gate を書き換える → load_approval が止める",
          raises(lambda: AP.load_approval(td / "approval_bad.json"), "合わない"))
    rr = td / "dirty_repo"
    rr.mkdir()
    for args in (["init", "-q"], ["config", "user.email", "t@t"], ["config", "user.name", "t"]):
        git(*args, cwd=rr)
    (rr / "x.txt").write_text("x")
    git("add", "x.txt", cwd=rr)
    git("commit", "-q", "-m", "x", cwd=rr)
    (rr / "x.txt").write_text("y")
    pv2 = AP.preview(card, prob, repo_dir=rr, judge_env=False)
    check("(N8) 追跡ファイルに変更がある状態の preview → 問題にする", any("追跡ファイルに変更" in p for p in pv2["problems"]), str(pv2["problems"])[:200])
finally:
    if co.exists():
        git("worktree", "remove", "--force", str(co))
    git("worktree", "prune")
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
