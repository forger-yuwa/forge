#!/usr/bin/env python3
"""1 つの run について**報告してよいか**を全ゲートで判定する (AGENTS.md の各ルールの束ね)。

これが通らない run の数値は報告しない。個別ツールを呼び出して VERDICT を集約する:

| # | 何を | ツール | 合格条件 |
| --- | --- | --- | --- |
| 1 | NaN/発散 | 残差 csv + 最終 res | 非有限ゼロ |
| 2 | 残差の収束 | `check_convergence.py --segment` | 判定区間で PASS |
| 3 | 結論量の準定常 | `check_cavity_steady.py` | 全量 STEADY |
| 4 | 壁解像 | `check_wall_resolution.py` | y1+ 超過面積 <= --yplus-frac |
| 5 | 保存性 | `cavity_eval.py` | 正味/片道 <= --mass-tol, CV 収支 <= --budget-tol |

usage:
  python3 tools/check_case_gates.py <run_dir> [--yplus-frac 2] [--mass-tol 0.01] [--budget-tol 0.05]

終了コード 0 = 全ゲート PASS。1 = どれか不合格 (理由を出力)。2 = 判定不能。
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
STOOLS = ROOT / "solver_density_cuda" / "tools"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--yplus-frac", type=float, default=2.0, help="y1+>1 を許す面積 [%]")
    ap.add_argument("--mass-tol", type=float, default=0.01, help="開口の |正味/片道| 許容")
    ap.add_argument("--budget-tol", type=float, default=0.05, help="CV エネルギー収支 残差 許容")
    a = ap.parse_args()
    rd = Path(a.run)
    fails, notes = [], []
    py = sys.executable

    print("======== ゲート判定: %s ========" % a.run)

    # --- 1. NaN / 発散 ---
    bad = 0
    for f in list(rd.glob("residual_history*.csv")):
        for ln in f.read_text().splitlines()[1:]:
            if re.search(r"(?i)\b(nan|inf)\b", ln):
                bad += 1
    if list(rd.glob("res_nan_*.h5")):
        bad += 1
    print("[1] NaN/発散       : %s" % ("OK" if bad == 0 else "**FAIL** (%d 箇所 / NaN ダンプ有り)" % bad))
    if bad:
        fails.append("NaN/発散")

    # --- 2. 残差の収束 (判定区間) ---
    seg = ["--segment"] if (rd / "stage_manifest.json").exists() else []
    rc, out = run([py, str(STOOLS / "check_convergence.py"), str(rd)] + seg)
    head = next((l for l in out.splitlines() if l.startswith("===")), "(出力なし)")
    print("[2] 残差の収束     : %s" % ("OK" if rc == 0 else "**FAIL**"))
    print("      %s" % head.strip())
    if rc != 0:
        fails.append("残差")
        notes.append(out)

    # --- 3. 結論量の準定常 ---
    rc3, out3 = run([py, str(HERE / "check_cavity_steady.py"), str(rd)])
    v3 = next((l for l in out3.splitlines() if l.startswith("VERDICT")), "")
    ok3 = "STEADY (全量)" in v3
    print("[3] 結論量の準定常 : %s   %s" % ("OK" if ok3 else "**FAIL**", v3.strip()))
    if not ok3:
        fails.append("準定常")
        notes.append(out3)

    # --- 4. 壁解像 ---
    rc4, out4 = run([py, str(STOOLS / "check_wall_resolution.py"), str(rd),
                     "--over-frac", str(a.yplus_frac)])
    v4 = next((l for l in out4.splitlines() if l.startswith("VERDICT")), "")
    print("[4] 壁解像 y1+     : %s   %s" % ("OK" if rc4 == 0 else "**FAIL**", v4.strip()))
    if rc4 != 0:
        fails.append("壁解像")

    # --- 5. 保存性 (開口の質量収支と CV エネルギー収支) ---
    j = rd / "cavity_eval.json"
    if not j.exists():
        run([py, str(HERE / "cavity_eval.py"), str(rd)])
    if j.exists():
        d = json.loads(j.read_text())
        imb = abs(float(d["field"].get("mdot_imbalance", float("nan"))))
        bud = abs(float(d.get("budget_residual", float("nan")))) if "budget_residual" in d else None
        ok5 = imb <= a.mass_tol
        print("[5] 保存性         : %s   開口の正味/片道 %.3e (許容 %.3g)"
              % ("OK" if ok5 else "**FAIL**", imb, a.mass_tol))
        if not ok5:
            fails.append("質量収支")
        if bud is not None:
            ok5b = bud <= a.budget_tol
            print("                      CV エネルギー収支 残差 %.3f (許容 %.3g) %s"
                  % (bud, a.budget_tol, "OK" if ok5b else "**FAIL**"))
            if not ok5b:
                fails.append("エネルギー収支")
    else:
        print("[5] 保存性         : **判定不能** (cavity_eval.json が作れない)")
        fails.append("保存性(判定不能)")

    print("\nGATES: %s" % ("PASS (報告してよい)" if not fails
                           else "FAIL — " + " / ".join(fails) + " (数値を報告しないこと)"))
    if fails and notes:
        print("--- 不合格の詳細 ---")
        for n in notes:
            print("\n".join(n.strip().splitlines()[-12:]))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
