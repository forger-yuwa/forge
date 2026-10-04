#!/usr/bin/env python3
"""SERN 2D/3D runner の準備経路の先頭 (設計点の控え・作動点の選択・気体状態・形状生成) を生産 YAML で通す回帰試験。

2026-09-27 に runner_sern.py の編集で `_dv` を消してしまい、`design_from_problem` が NameError で止まった
(codex result chi-default-2 M1)。選別の単体試験では検出できなかったので、実際の問題 YAML で形状生成まで呼ぶ。
メッシュ生成・forge 実行はしない。
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from forge_design.probdef import load_problem  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402,F401

CASE = os.path.join(HERE, "..", "..", "case", "46.sern_design")
YAMLS = ["problem_moo_frozen_tp_cycle3op.yaml", "problem_3d_prod_m6on_g1.yaml"]
fails = 0
for y in YAMLS:
    p0 = load_problem(os.path.join(CASE, y))
    ops = [o["name"] for o in (p0.spec.get("operating_points") or [])] or [None]
    for op in ops:
        p = load_problem(os.path.join(CASE, y))   # runner と同じく作動点ごとに読み直す (select_operating_point は p を上書きする)
        try:
            d0 = R2.design_snapshot(p); R2.select_operating_point(p, op); R2.gas_states(p)
            kern, design, fr, th = R2.design_from_problem(p, design=d0)
            ok = design is not None
        except Exception as e:
            ok = False; print(f"  {type(e).__name__}: {e}")
        fails += not ok
        print(("PASS " if ok else "FAIL ") + f"{y} op={op}: 設計点の控え → 作動点 → 気体 → 形状生成")
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
