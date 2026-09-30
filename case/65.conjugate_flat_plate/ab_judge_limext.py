#!/usr/bin/env python3
r"""case/65 B3 延長の判定 (plan `boundary-cht-conjugate-flat-plate.md` §4.5、事前登録どおり)。

    python3 ab_judge_limext.py <run_A (limiter 2, run_0017)> <B 前半 (run_0018)> <B 延長 (run_0019)>

- B の判定区間: 累計 step [36000,40000)・[40000,44000)・[44000,48000) (= 延長 run の [24000,36000))。A の比較区間は [8000,12000)。
- 能動残差の全列の区間平均と下流 slip の規格化 P 振幅が、最後の 2 区間で A の 1/10 以下、かつ 3 区間で厳密減少
  → 「同程度の停滞が無制限再構成でも残る」を棄却 (再構成依存を支持)。
- 未達の列が残り、3 区間の max/min ≤ 1.1 → 「累計 48000 step までの延長で過渡が抜ける」を棄却。
- それ以外 (なお下降中など) は判定不能。自動延長しない。
- 別判定 (併記): `check_convergence.py` を limiter 0 の同一設定区間 (run_0018 + run_0019 を連結) で、
  `check_quasisteady.py --series-csv --abs-scale 1` を延長 run の下流 slip P (q∞ 規格化、閾値 0.002) と界面 q
  (評価窓平均 |q| 規格化、閾値 0.006) で。領域別の残差の絶対二乗和と割合を 4000 step 区間ごと。
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ab_judge_cfl import COLS, Q_INF, fmt  # noqa: E402
from ab_judge_lim import down_amp, resid  # noqa: E402
from ab_series import q_scale  # noqa: E402

TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
LI = 4000


def main():
    rA, rB1, rB2 = (Path(p) for p in sys.argv[1:4])
    print(f"=== B3 延長の判定  A={rA.name} [8000,12000)  B={rB1.name}+{rB2.name} 累計 [36000,48000)")
    RA, _ = resid(rA, 8000, 12000, LI)
    DA, _ = down_amp(rA, 8000, 12000, LI)
    RB, _ = resid(rB2, 24000, 36000, LI)
    DB, dB = down_amp(rB2, 24000, 36000, LI)
    allA = {**RA, "down_P": DA}; allB = {**RB, "down_P": DB}
    print("  区間平均 (A は 1 区間、B は 3 区間):")
    for k in allB:
        print(f"    {k:9s} A {allA[k][0]:.3e} | B {fmt(allB[k])}  B/A {'/'.join(f'{v:.4f}' for v in allB[k] / allA[k][0])}  (1/10 閾値 {allA[k][0]/10:.3e})")
    low = all((allB[k][1:] <= allA[k][0] / 10).all() and allB[k][0] > allB[k][1] > allB[k][2] for k in allB)
    flat_fail = [k for k in allB if not (allB[k][1:] <= allA[k][0] / 10).all()]
    flat = bool(flat_fail) and all(allB[k].max() / allB[k].min() <= 1.1 for k in flat_fail)
    # 領域別の残差 (4000 step 区間ごと)
    rr = dB["region_res2"]; st = dB["step"]; names = [str(x) for x in dB["region_names"]]
    print("  領域別の残差二乗和 (区間ごと、絶対値と割合):")
    for j, nm in enumerate(dB["res_names"]):
        parts = []
        for k in range(3):
            m = (st > 24000 + k * LI) & (st <= 24000 + (k + 1) * LI)
            tot = rr[m, j, :].sum(axis=0)
            parts.append("[" + ", ".join(f"{names[i]} {tot[i]:.2e} ({tot[i]/tot.sum()*100:.0f}%)" for i in range(len(names))) + "]")
        print(f"    {str(nm):9s} " + " ".join(parts))
    # 収束 (limiter 0 の同一設定区間を連結)
    with tempfile.TemporaryDirectory() as td:
        t = Path(td) / "run_limiter0_joined"; t.mkdir()
        with open(t / "residual_history.csv", "w") as out:
            with open(rB1 / "residual_history.csv") as f:
                out.write(f.read())
            n1 = int(np.load(rB1 / "ab_series.npz", allow_pickle=True)["step"][-1])
            with open(rB2 / "residual_history.csv") as f:
                next(f)
                for line in f:
                    s, rest = line.split(",", 1)
                    out.write(f"{int(s) + n1},{rest}")
        p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(t)], capture_output=True, text=True)
        print("  check_convergence (run_0018 + run_0019 連結、limiter 0 の同一設定区間):")
        print("    " + "\n    ".join((p.stdout + p.stderr).strip().splitlines()[-7:]))
    # 準定常 (延長 run の後半)
    qs = q_scale()
    g = dB["group"]
    Pn = dB["P"][:, g == "down"] / Q_INF; qn = dB["q"] / qs
    for tag, arr, thr in (("下流 slip の P", Pn, 0.002), ("界面 q", qn, 0.006)):
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            cols = [f"c{i}" for i in range(arr.shape[1])]
            f.write("step," + ",".join(cols) + "\n")
            for s_, row in zip(st, arr):
                f.write(f"{s_}," + ",".join(f"{v:.15e}" for v in row) + "\n")
            path = f.name
        p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", path, "--series-cols", ",".join(cols),
                            "--tail", "0.5", "--drift", str(thr), "--osc", str(thr), "--abs-scale", "1"], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l]
        print(f"  check_quasisteady ({tag}、延長 run の後半、閾値 {thr}): {ov[-1].strip() if ov else '出力なし'} (rc {p.returncode})")
        Path(path).unlink()
    if low:
        print("VERDICT: 「同程度の停滞が無制限再構成でも残る」を棄却 (再構成依存を支持)")
    elif flat:
        print(f"VERDICT: 「累計 48000 step までの延長で過渡が抜ける」を棄却 (未達: {', '.join(flat_fail)})")
    else:
        print(f"VERDICT: 判定不能 (未達: {', '.join(flat_fail) or 'なし (減少が厳密でない)'}、なお変化中)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
