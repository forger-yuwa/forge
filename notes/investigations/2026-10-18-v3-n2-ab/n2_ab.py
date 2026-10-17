#!/usr/bin/env python3
"""V3 13a の N2 の FAIL の判別の A/B (plan tooling-design-orchestrator §6「N2 の FAIL の扱いと判別の A/B」、codex 諮問 v3-n2-fail、2026-10-18)。

変える因子は数値の生成のコードの版だけ: A = a36a87ff (試走 A の p0 の版)、B = b1a02511 (今の版)。固定: カード `case/90.qs_trial_a/card_p0.yaml`
(v1、sha256 23caf0b7…)、Python と依存 (design/.venv-opt)、`--axis-cap-frac 0.03`・`--radial-ratio-max 1.2`。カードの生成 (wt_card) → 壁 → 配点の表 →
節点数の決定 (hand_mesh build) まで各 1 回。forge・全格子の変換・品質はしない。

B の今の hand_mesh build は節点数の上限で止まるので、上限の照合 (`runner_axismach.check_node_budget`) だけを、同じ関数を非常に大きな上限で呼ぶ形に
差し替えて表と節点数を取り出す (数え方は変えない)。A の旧版は上限を自動で広げるのでそのまま回す。両方の NS・Euler の節点数は、出力の問題 YAML から
今の数え方 (check_node_budget) で数え直す。

期待: A が旧の値 411,240 / 433,590、B が今の値 412,160 / 434,560 を再現すれば、コードの版の変更で差を説明できる (#6e 単独への帰属は保留)。
A と B が一致すれば、コードの版だけの仮説を棄却する。どちらかが記録を再現しなければ対照の不成立。

実行: /home/sano/work/forge/design/.venv-opt/bin/python notes/investigations/2026-10-18-v3-n2-ab/n2_ab.py   (結果は同じディレクトリの result.json)
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PY = sys.executable
CARD = REPO / "case" / "90.qs_trial_a" / "card_p0.yaml"
A_COMMIT, B_COMMIT = "a36a87ff", "b1a02511"
td = Path(tempfile.mkdtemp(prefix="n2_ab_"))


def worktree(commit: str) -> Path:
    d = td / f"wt_{commit}"
    subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", str(d), commit], check=True, capture_output=True)
    return d


def gen(wt: Path, out: Path):
    env = dict(os.environ, PYTHONPATH=str(wt / "design"), PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([PY, "-m", "forge_design.intake.wt_card", str(CARD), "--out", str(out)], cwd=wt / "design", env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"wt_card ({wt.name}): {r.stderr[-600:]}")


def summarize(wt: Path, pc: Path, prob: Path, rec: Path) -> dict:
    """今の版の数え方で NS・Euler の節点数を数え直し、比べる量をそろえる。"""
    sys.path.insert(0, str(REPO / "design"))
    from forge_design.evaluate import runner_axismach as RA
    from forge_design.probdef import load_problem
    p = load_problem(prob)
    big = {"max_nodes": 10 ** 9, "max_nodes_reason": "A/B の数え直し"}
    ns = RA.check_node_budget(RA.mesh_params(p, 1.0, ni=561, nj=97, wall_first_frac=4.5e-5), {**dict(p.mesh), **big}, "mesh")
    eu = RA.check_node_budget(RA.mesh_params_euler(p, 1.0), {**dict(RA.mesh_euler_block(p)), **big}, "mesh_euler")
    d = yaml.safe_load(prob.read_text())
    dc = yaml.safe_load(pc.read_text())
    r = json.loads(rec.read_text())
    g = dc.get("geometry") or {}
    return {"r_throat": (dc.get("spec") or {}).get("r_throat"), "M_knot": g.get("M_knot"), "n_axis_inv": g.get("n_axis_inv"),
            "physical_wall": r.get("physical_wall"), "design_wall": r.get("design_wall"), "ni_min": r.get("ni_min"), "ni": r.get("ni"),
            "nj": r.get("nj"), "nodes_record": r.get("nodes"), "nodes_ns": ns["n_nodes"], "nodes_euler": eu["n_nodes"],
            "wall_first_frac_table": d["mesh"].get("wall_first_frac_table"), "x_density_table": d["mesh"].get("x_density_table")}


try:
    out = {}
    wa = worktree(A_COMMIT)
    (td / "A").mkdir()
    gen(wa, td / "A" / "problem_p0_card.yaml")
    hm_a = wa / "case" / "71.m55_combustion_wt" / "_mesh" / "hand_mesh.py"
    r = subprocess.run([PY, str(hm_a), "build", "--design", str(wa / "design"), "--card", str(CARD), "--problem", str(td / "A" / "problem_p0_card.yaml"),
                        "--out-yaml", str(td / "A" / "problem_p0.yaml"), "--record", str(td / "A" / "build_record.json"), "--axis-cap-frac", "0.03",
                        "--radial-ratio-max", "1.2", "--ratio-summary", str(wa / "case/71.m55_combustion_wt/_mesh/c45/h_over_est_ratio_summary.json"),
                        "--plan-item", "N2 A/B"], capture_output=True, text=True, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    if r.returncode != 0:
        raise RuntimeError(f"A の hand_mesh build: {r.stderr[-600:]}")
    out["A"] = summarize(wa, td / "A" / "problem_p0_card.yaml", td / "A" / "problem_p0.yaml", td / "A" / "build_record.json")
    # B: 今の版 (このリポジトリの作業ツリー = b1a02511 の後の版でもよいが、版を固定するため worktree)
    wb = worktree(B_COMMIT)
    (td / "B").mkdir()
    gen(wb, td / "B" / "problem_p0_card.yaml")
    hm_b = wb / "case" / "71.m55_combustion_wt" / "_mesh" / "hand_mesh.py"
    drv = td / "B" / "run_b.py"
    drv.write_text(f'''import importlib.util, sys, argparse
spec = importlib.util.spec_from_file_location("hm_b", {str(hm_b)!r}); hm = importlib.util.module_from_spec(spec); spec.loader.exec_module(hm)
hm._imports({str(wb / "design")!r})
orig = hm.RA.check_node_budget
hm.RA.check_node_budget = lambda mp, m, block, allow_hint=None: orig(mp, {{**dict(m), "max_nodes": 10 ** 9, "max_nodes_reason": "A/B"}}, block)
a = argparse.Namespace(mode="build", design={str(wb / "design")!r}, card={str(CARD)!r}, problem={str(td / "B" / "problem_p0_card.yaml")!r},
    ratio_summary={str(wb / "case/71.m55_combustion_wt/_mesh/c45/h_over_est_ratio_summary.json")!r}, out_yaml={str(td / "B" / "problem_p0.yaml")!r},
    record={str(td / "B" / "build_record.json")!r}, run_dir=None, axis_cap_frac=0.03, radial_ratio_max=1.2, plan_item="N2 A/B", approval=None, sanity=None)
hm.build(a)
''')
    r = subprocess.run([PY, str(drv)], capture_output=True, text=True, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    if r.returncode != 0:
        raise RuntimeError(f"B の hand_mesh build: {r.stderr[-600:]}")
    out["B"] = summarize(wb, td / "B" / "problem_p0_card.yaml", td / "B" / "problem_p0.yaml", td / "B" / "build_record.json")
    A, B = out["A"], out["B"]
    rep_A = (A["nodes_ns"], A["nodes_euler"]) == (411240, 433590)
    rep_B = (B["nodes_ns"], B["nodes_euler"]) == (412160, 434560)
    same = (A["nodes_ns"], A["nodes_euler"]) == (B["nodes_ns"], B["nodes_euler"])
    verdict = ("コードの版の変更で差を説明できる (A は旧の値、B は今の値を再現; #6e 単独への帰属は保留)" if rep_A and rep_B else
               "A と B が一致 → コードの版だけの仮説を棄却 (旧の実行の環境・入力・オプションを候補へ)" if same else
               "対照の不成立 (どちらかが記録を再現しない)")
    diffs = sorted(k for k in A if A[k] != B[k])
    res = {"A_commit": A_COMMIT, "B_commit": B_COMMIT, "verdict": verdict, "rep_A": rep_A, "rep_B": rep_B, "differs": diffs,
           "A": {k: v for k, v in A.items() if not k.endswith("_table")}, "B": {k: v for k, v in B.items() if not k.endswith("_table")},
           "tables_equal": {k: A[k] == B[k] for k in ("wall_first_frac_table", "x_density_table")}}
    (HERE / "result.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str) + "\n")
    print(json.dumps({k: res[k] for k in ("verdict", "differs", "tables_equal")}, ensure_ascii=False, indent=1))
    print(json.dumps({t: {k: res[t][k] for k in ("r_throat", "M_knot", "n_axis_inv", "ni_min", "ni", "nj", "nodes_ns", "nodes_euler")} for t in ("A", "B")},
                     ensure_ascii=False, indent=1))
finally:
    for c in (A_COMMIT, B_COMMIT):
        subprocess.run(["git", "-C", str(REPO), "worktree", "remove", "--force", str(td / f"wt_{c}")], capture_output=True)
    subprocess.run(["git", "-C", str(REPO), "worktree", "prune"], capture_output=True)
