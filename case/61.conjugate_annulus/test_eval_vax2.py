#!/usr/bin/env python3
r"""V-ax2 評価器の自己試験 — 正例・負例・判別 A/B (plan §6 V-ax2 / V-ax2b の M2)。

流体計算をしない。登録格子 (流体 32×4・固体 16 層) の `mesh.h5` と `solid.h5` から、評価器が読むファイル一式を
**合成**した run を一時ディレクトリに作る (axcht.synth_run):

  - 共役壁の荷重 $Q_{f,i}=q_*A^r_{{\rm fluid},i}$、`iface_q_eff`=$q_*$、固体は $Ku=b+E^{\rm T}Q_f$ の直接解、`Ts`=固体界面温度
  - 内壁の荷重 $Q_{f,i}=-(Q'/r_a)A^r_{{\rm fluid},i}$

期待 (結果を見る前に書いた):

  1. **正例**: `eval_vax2.py --no-qs` が (a)〜(e) すべて PASS (rc 0)
  2. **負例** (`--simulate-old-denominator` = $q_{\rm bad}=q_{\rm eff}A^r_{\rm fluid}/A_{\rm planar}$): (d)(e) が FAIL、(a)(b)(c) は PASS のまま (rc 1)
  3. **判別 A/B**: 合成荷重 $Q_i=q_*A_i^r$ を A = $A_i^r$ で割ると PASS、B = $A_{{\rm planar},i}$ で割ると FAIL (固体側・流体側とも)

`--run <run_dir>` を付けると、実 run でも「無改変の (d) が PASS し、模擬すると FAIL する」ことを確かめる
(無改変の (d) が FAIL の run では判定不能と出す)。

    python3 case/61.conjugate_annulus/test_eval_vax2.py [--run case/61.conjugate_annulus/run_0002_annulus]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import axcht
from eval_vax2 import PID_CJ, PID_IN, NAME_CJ, NAME_IN, RA, analytic

EVAL = axcht.HERE / "eval_vax2.py"


def evaluate(run: Path, simold: bool):
    cmd = [sys.executable, str(EVAL), str(run), "--no-qs"] + (["--simulate-old-denominator"] if simold else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    js = run / ("vax2_eval_simold.json" if simold else "vax2_eval.json")
    if not js.exists():
        print(p.stdout + p.stderr)
        raise SystemExit(f"評価器が結果を書かなかった: {' '.join(cmd)}")
    res = json.loads(js.read_text())
    return {it["name"][:3]: it["pass"] for it in res["items"]}, res["verdict"], p.returncode, p.stdout


def check(label, cond, detail=""):
    print(f"  [{'OK' if cond else 'NG'}] {label}" + (f"  ({detail})" if detail else ""))
    return cond


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mesh", default=str(axcht.HERE / "mesh" / "annulus_r32_x4.h5"))
    ap.add_argument("--solid", default=str(axcht.HERE / "mesh" / "solid_s16_x4.h5"))
    ap.add_argument("--run", default=None, help="実 run でも負例を確かめる")
    a = ap.parse_args()
    an = analytic()
    good = True

    with tempfile.TemporaryDirectory() as td:
        run = Path(td) / "synth"
        Tw = axcht.synth_run(run, Path(a.mesh), Path(a.solid), axcht.HERE / "template" / "solverConfig.yaml",
                             axcht.HERE / "template" / "bcondConfig.yaml",
                             walls=[(PID_IN, NAME_IN, -an["Qp"] / RA)], pid_cj=PID_CJ, name_cj=NAME_CJ,
                             q_cj=an["q"])
        print(f"合成 run: 固体だけの直接解で T_w {Tw.min():.5f} .. {Tw.max():.5f} K (解析 {an['Tw']:.5f}、"
              f"誤差 max {abs(Tw - an['Tw']).max() / an['drop'] * 100:.4f} % of 降下) — 荷重移送+固体離散化の誤差")
        print("1. 正例")
        it, v, rc, out = evaluate(run, False)
        good &= check("VERDICT PASS / rc 0", v == "PASS" and rc == 0, f"{v}, rc {rc}")
        if v != "PASS":
            print(out)
        print("2. 負例 (分母修正を外した q_eff の模擬)")
        it, v, rc, out = evaluate(run, True)
        good &= check("(d)(e) FAIL", (not it["(d)"]) and (not it["(e)"]))
        good &= check("(a)(b)(c) は PASS のまま", it["(a)"] and it["(b)"] and it["(c)"])
        good &= check("VERDICT FAIL / rc 1", v == "FAIL" and rc == 1, f"{v}, rc {rc}")

    print("3. 判別 A/B (合成荷重 Q_i = q* A_i^r)")
    ab = axcht.synthetic_ab(Path(a.mesh), PID_CJ, Path(a.solid), an["q"])
    for side, r in ab.items():
        good &= check(f"{side}: A (÷A^r) PASS", r["A_pass"], f"max|q-q*| {r['A_err']:.3e} W/m², 許容 {r['tol']:.4f}")
        good &= check(f"{side}: B (÷A_planar) FAIL", not r["B_pass"], f"max|q-q*| {r['B_err']:.4f} W/m²")

    if a.run:
        run = Path(a.run)
        print(f"4. 実 run: {run}")
        it0, v0, _, _ = evaluate(run, False)
        it1, v1, _, _ = evaluate(run, True)
        if not it0["(d)"]:
            print("  判定不能: 無改変の (d) が FAIL (負例の検出力を示せない)")
            good = False
        else:
            good &= check("無改変 (d) PASS → 模擬 (d)(e) FAIL", (not it1["(d)"]) and (not it1["(e)"]))

    print(f"\nVERDICT: {'PASS' if good else 'FAIL'}")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
