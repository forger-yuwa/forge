#!/usr/bin/env python3
"""判定専用の環境で成果物を書き出す通しの試験 (plan tooling-design-orchestrator §4.7.1 (5)・(6)、§5.1 #2c)。

case/71 の run_0014_chain_v2 (手元) の写しから、承認の記録 (版 2 のカード・問題 YAML・判定の Python・外部の道具) を作り、
`gate_env deliverables` で成果物を書き出す。作業側の環境変数を渡さず、承認のときに記録した道具 (gmsh 用の Python・Excel 用の Python・
FreeCAD・変換器) だけで部品がそろうか、MANIFEST に承認と合う gate が書かれるか、判定の側 (deliverables_check に gate を渡す) が
その成果物を PASS にするかを見る。

**このリポジトリの HEAD を承認の commit に使う**ので、追跡ファイルに変更の無い状態で回す。checkout は一時のディレクトリに作り最後に消す。
変換器は FORGE_CONVERTER (無ければ solver_density_cuda/.build-native/card4a/convertGmshToForge)。

  (P1) gate_env deliverables が rc 0 で、MANIFEST の verdict PASS・部品の pass がそろう
  (P2) MANIFEST の gate が承認の記録 (gate の sha256・承認の記録の sha256・承認の commit) と合う
  (P3) deliverables_check (gate つき) が PASS
  (N1) 同じ成果物を gate の違う承認 (別の記録) で見ると INDETERMINATE

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_gate_deliverables_tests.py
所要: 約 2〜3 分 (FreeCAD・gmsh)。run・道具が無ければ飛ばす (rc 0 と「飛ばした」)。
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
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

from forge_design.feedback import gate_env as GE  # noqa: E402
from forge_design.feedback import wt_qualify as Q  # noqa: E402
from forge_design.intake import approval as AP  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402

C71 = REPO / "case" / "71.m55_combustion_wt"
R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", C71 / "run_0014_chain_v2"))
CARD, PROB = C71 / "card_m55_3coef.yaml", C71 / "problem_m55_p5_settle_v2.yaml"
CONV = os.environ.get("FORGE_CONVERTER") or str(REPO / "solver_density_cuda" / ".build-native" / "card4a" / "convertGmshToForge")
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def git(*args, cwd=REPO):
    return subprocess.run(["git", "-C", str(cwd), *args], capture_output=True, text=True)


need = [R14 / "nozzle.h5", R14 / "wall_repr.json", CARD, PROB, Path(CONV)]
miss = [str(p) for p in need if not p.exists()]
if miss:
    print(f"飛ばした: 無い {miss}")
    sys.exit(0)
if git("status", "--porcelain", "--untracked-files=no").stdout.strip():
    print("止めた: 作業ツリーの追跡ファイルに変更がある (HEAD を承認の commit に使う; commit してから回す)")
    sys.exit(2)

SOURCES = {s: {"kind": "user", "ref": "case/71 の要求"} for s in ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall",
                                                                "tolerances", "p_ambient")}
SOURCES["budget"] = {"kind": "user", "ref": "試験の仮の承認"}
td = Path(tempfile.mkdtemp(prefix="gate_deliv_tests_"))
os.environ[GE.GATE_ROOT_ENV] = str(td / "gates")
head = git("rev-parse", "HEAD").stdout.strip()
try:
    d = yaml.safe_load(CARD.read_text())
    d.update(format=W.CARD_FORMAT_V2, budget={"gpu_hours": 4.0, "disk_gb": 50.0}, sources=SOURCES)
    card = td / "card.yaml"
    card.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    prob = td / "problem.yaml"
    shutil.copy(PROB, prob)
    os.environ["FORGE_CONVERTER"] = CONV          # 承認のときに記録する変換器 (current_tools が拾う)
    pv = AP.preview(card, prob)
    check("前提: preview に問題が無い", not pv["problems"], str(pv["problems"])[:400])
    if pv["problems"]:
        raise SystemExit(1)
    (td / "preview.json").write_text(json.dumps(pv, ensure_ascii=False))
    apath = td / "card_approval.json"
    AP.write(td / "preview.json", "試験の仮の承認", apath)
    # run の写し: h5 は symlink、ほかは複製 (元の run に書かない)
    run = td / R14.name
    run.mkdir()
    for f in R14.iterdir():
        if f.is_file():
            if f.suffix == ".h5":
                (run / f.name).symlink_to(f)
            else:
                shutil.copy2(f, run / f.name)
    out = td / "deliverables"
    q = subprocess.run([sys.executable, "-m", "forge_design.feedback.gate_env", "deliverables", "--approval", str(apath), "--",
                        str(run), "--card", str(card), "--problem", str(prob), "--out", str(out)],
                       capture_output=True, text=True, cwd=str(HERE.parent))
    man = json.loads((out / "MANIFEST.json").read_text()) if (out / "MANIFEST.json").is_file() else {}
    parts = man.get("parts") or {}
    check("(P1) gate_env deliverables が rc 0・MANIFEST の verdict PASS・cad/mesh/fluent の pass がそろう",
          q.returncode == 0 and man.get("verdict") == "PASS" and all((parts.get(p) or {}).get("pass") is True for p in ("cad", "mesh", "fluent")),
          f"rc {q.returncode}, verdict {man.get('verdict')}, failed {man.get('failed_parts')}, errors {list((man.get('errors') or {}))} "
          f"{q.stderr.strip()[-300:]}")
    g = man.get("gate") or {}
    apd = json.loads(apath.read_text())
    check("(P2) MANIFEST の gate が承認の記録と合う (gate の sha256・承認の記録の sha256・承認の commit)",
          g.get("gate_sha256") == apd["gate_sha256"] and g.get("approval_sha256") == AP.sha256_file(apath) and g.get("commit") == head,
          str(g)[:300])
    G = {"gate_sha256": apd["gate_sha256"], "approval_sha256": AP.sha256_file(apath)}
    fld = {"step": int(sorted(int(p.stem.split("_")[1]) for p in run.glob("res_*.h5") if p.stem.split("_")[1].isdigit())[-1])}
    fld["res"] = f"res_{fld['step']}.h5"
    r = Q.deliverables_check(run, out, card, field=fld, problem_path=prob, gate=G)
    check("(P3) deliverables_check (gate つき) が PASS", r["status"] == "PASS", f"{r.get('why') or r.get('problems')}")
    r = Q.deliverables_check(run, out, card, field=fld, problem_path=prob, gate={**G, "gate_sha256": "0" * 64})
    check("(N1) gate の違う承認で見ると INDETERMINATE", r["status"] == "INDETERMINATE", r.get("why"))
finally:
    co = td / "gates" / head
    if co.exists():
        git("worktree", "remove", "--force", str(co))
    git("worktree", "prune")
    shutil.rmtree(td, ignore_errors=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
