#!/usr/bin/env python3
r"""V-ax2b 評価器の自己試験 (流体計算なし) — 正例・荷重移送の予測との照合・判別 A/B。

plan §6 V-ax2b。登録格子の `mesh.h5` と `solid.h5` から合成 run を作る (case/61 `axcht.synth_run`):
共役壁の荷重 $Q_{f,i}=q\,A^r_{{\rm fluid},i}$ (流体半辺の厳密な面積分)、熱い壁 $Q_{f,i}=-q\,A^r_{{\rm fluid},i}$、
固体は $Ku=b+E^{\rm T}Q_f$ の直接解、`Ts` = 固体界面温度。

期待 (結果を見る前に書いた):

  1. 各格子で (b)(c) PASS、(d) は合成出力が倍精度なので PASS (恒等式は作り方から自明 — 評価器の経路の確認だけ)
  2. 界面温度の最大誤差 (% of 降下) が plan §6 V-ax2b の**参考予測** (codex の独立計算、流体半辺積分の荷重を固体円板に
     与えた荷重移送だけの試験) 0.9421 / 0.2983 / 0.08796 % (N_r = 8 / 16 / 32 一様) と一致すること (相対 1e-3 以内)。
     **これは連成計算の予測ではない** (登録文どおり)
  3. 判別 A/B: 合成荷重 $Q_i=q_*A_i^r$ を A = $A_i^r$ で割ると PASS、B = $A_{{\rm planar},i}$ で割ると FAIL (固体側・流体側)

    python3 case/62.conjugate_disk/test_eval_vax2b.py
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "61.conjugate_annulus"))
import axcht  # noqa: E402
from eval_vax2b import NAME_CJ, NAME_HOT, PID_CJ, PID_HOT, analytic  # noqa: E402

PRED = {"r8_u": 0.9421, "r16_u": 0.2983, "r32_u": 0.08796}


def main():
    an = analytic()
    good = True
    with tempfile.TemporaryDirectory() as td:
        for tag in ("r8_u", "r16_u", "r32_u", "r32_g1p1"):
            mesh, solid = HERE / "mesh" / f"disk_{tag}.h5", HERE / "mesh" / f"solid_disk_{tag}.h5"
            run = Path(td) / tag
            axcht.synth_run(run, mesh, solid, HERE / "template" / "solverConfig.yaml",
                            HERE / "template" / "bcondConfig.yaml",
                            walls=[(PID_HOT, NAME_HOT, -an["q"])], pid_cj=PID_CJ, name_cj=NAME_CJ, q_cj=an["q"])
            p = subprocess.run([sys.executable, str(HERE / "eval_vax2b.py"), str(run), "--no-qs"],
                               capture_output=True, text=True)
            res = json.loads((run / "vax2b_eval.json").read_text())
            items = {it["name"][:3]: it["pass"] for it in res["items"]}
            e = res["Tw_err_pct_max"]
            ok1 = items["(b)"] and items["(c)"] and res["identity_pass"]
            line = f"{tag:10s}: T_w 誤差 max {e:.5f} % of 降下  (b)(c)(d) {'OK' if ok1 else 'NG'}"
            good &= ok1
            if tag in PRED:
                okp = abs(e / PRED[tag] - 1) <= 1e-3
                good &= okp
                line += f"  予測 {PRED[tag]} % と {'一致' if okp else '**不一致**'} (比 {e/PRED[tag]:.5f})"
            print(line)
            if not ok1:
                print(p.stdout[-3000:])
    print("判別 A/B (合成荷重 Q_i = q* A_i^r、N_r=32 非一様)")
    ab = axcht.synthetic_ab(HERE / "mesh" / "disk_r32_g1p1.h5", PID_CJ, HERE / "mesh" / "solid_disk_r32_g1p1.h5", an["q"])
    for side, r in ab.items():
        okA, okB = r["A_pass"], not r["B_pass"]
        good &= okA and okB
        print(f"  {side}: A (÷A^r) {'PASS' if r['A_pass'] else 'FAIL'} ({r['A_err']:.2e}) / "
              f"B (÷A_planar) {'PASS' if r['B_pass'] else 'FAIL'} ({r['B_err']:.4f} W/m²) -> {'OK' if okA and okB else 'NG'}")
    print(f"\nVERDICT: {'PASS' if good else 'FAIL'}")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
