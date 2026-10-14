#!/usr/bin/env python3
r"""③ベル (M) L100: δ* 補正壁の判定 (`judge_corrected_wall.py`、凍結 `58199b1c`) の内訳 — 記録用 (判定には使わない)。
plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-4)。組は判定と同じ (U = run_0019/run_0018、V-ext = run_0026/run_0024、
V-zero = run_0027/run_0024、格子 = run_0024 vs run_0025)。出力: corrected_wall_breakdown.json (区間ごとの壁圧の相対差の RMS・最大・
平均、軸上 Mach の差の RMS、出口面の Mach の半径分布) と図 corrected_wall_breakdown.png (git 管理外)。"""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import judge_corrected_wall as j  # noqa: E402

RUNS = {"U_ns": "run_0019_tp_m_ns_l100", "U_eu": "run_0018_tp_m_euler_l100_pin", "vext": "run_0026_tp_m_ns_l100_corr_vext",
        "vzero": "run_0027_tp_m_ns_l100_corr_vzero", "target": "run_0024_tp_m_euler_l100_inv",
        "target_fine": "run_0025_tp_m_euler_l100_inv_fine"}
PAIRS = {"U": ("U_ns", "U_eu"), "Vext": ("vext", "target"), "Vzero": ("vzero", "target"), "grid": ("target_fine", "target")}
REGS = [(0, 0.41), (0.41, 1.17), (1.17, 2), (2, 4), (4, 8), (8, 13), (13, 18), (18, 99)]
R_EXIT = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5]


def main():
    D = {k: j.load(j.CASE / v) for k, v in RUNS.items()}
    L = float(D["U_ns"]["x"][-1, -1]); xg = np.linspace(0.0, L, 4001)
    out = {"runs": RUNS, "regions": REGS, "pairs": {}}
    curves = {}
    for n, (a, b) in PAIRS.items():
        A, B = D[a], D[b]
        ep = np.interp(xg, A["x"][:, -1], A["P"][:, -1]) / np.interp(xg, B["x"][:, -1], B["P"][:, -1]) - 1.0
        em = np.interp(xg, A["x"][:, 0], A["M"][:, 0]) - np.interp(xg, B["x"][:, 0], B["M"][:, 0])
        rows = []
        for lo, hi in REGS:
            m = (xg >= lo) & (xg < hi)
            rows.append({"x": [lo, min(hi, L)], "pw_rms": float(np.sqrt(np.mean(ep[m] ** 2))), "pw_max": float(np.max(np.abs(ep[m]))),
                         "pw_mean": float(np.mean(ep[m])), "axM_rms": float(np.sqrt(np.mean(em[m] ** 2)))})
        out["pairs"][n] = {"regions": rows,
                           "exit_M": {"r": R_EXIT, "a": np.interp(R_EXIT, A["r"][-1], A["M"][-1]).tolist(),
                                      "b": np.interp(R_EXIT, B["r"][-1], B["M"][-1]).tolist()}}
        curves[n] = (ep, em, (A["r"][-1], A["M"][-1], B["r"][-1], B["M"][-1]))
    (HERE / "corrected_wall_breakdown.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in font_manager.findSystemFonts([os.path.expanduser("~/.fonts")]):
        font_manager.fontManager.addfont(f)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
    col = {"U": "#888888", "Vext": "#2f6db5", "Vzero": "#c9721f", "grid": "#3a8f5c"}
    lab = {"U": "補正なし (NS run_0019 − Euler run_0018)", "Vext": "補正 V-ext (NS run_0026 − 狙い run_0024)",
           "Vzero": "補正 V-zero (NS run_0027 − 狙い run_0024)", "grid": "狙いの格子の差 (run_0025 − run_0024)"}
    fig, ax = plt.subplots(3, 1, figsize=(10.5, 11))
    for n in ("U", "Vext", "Vzero"):
        ax[0].plot(xg, curves[n][0] * 100, color=col[n], lw=1.3, label=lab[n])
        ax[1].plot(xg, curves[n][1], color=col[n], lw=1.2, label=lab[n])
    ax[1].plot(xg, curves["grid"][1], color=col["grid"], lw=1.0, ls="--", label=lab["grid"])
    for a_ in ax[:2]:
        a_.axhline(0, color="0.4", lw=0.8); a_.axvspan(0, 1.17, color="0.9", zorder=0); a_.grid(alpha=0.3); a_.set_xlabel("x [r_t]")
    ax[0].set_ylabel("壁圧の相対差 [%] (NS / 狙い − 1)"); ax[0].legend(fontsize=8.5)
    ax[0].set_title("壁圧: 灰色の帯 = 推定区間 x < 1.17 (較正の外)", fontsize=11)
    ax[1].set_ylabel("軸上の Mach の差"); ax[1].legend(fontsize=8.5)
    for n in ("U", "Vext", "Vzero"):
        rA, MA, rB, MB = curves[n][2]
        ax[2].plot(rA, MA, color=col[n], lw=1.3, label=lab[n].split(" (")[0] + " の NS")
    rA, MA, rB, MB = curves["U"][2]; ax[2].plot(rB, MB, color="0.3", lw=1.0, ls=":", label="狙い (補正なし) = Euler run_0018")
    rA, MA, rB, MB = curves["Vext"][2]; ax[2].plot(rB, MB, color="k", lw=1.0, ls="--", label="狙い (補正後) = Euler run_0024")
    ax[2].axvline(j.R_C, color="0.5", lw=0.8); ax[2].set_xlabel("r [r_t] (出口面)"); ax[2].set_ylabel("Mach"); ax[2].legend(fontsize=8.5)
    ax[2].grid(alpha=0.3); ax[2].set_title("出口面の Mach (縦線 r = 6.5 が判定の範囲の端)", fontsize=11)
    fig.tight_layout(); fig.savefig(HERE / "corrected_wall_breakdown.png", dpi=130); plt.close(fig)
    for n, v in out["pairs"].items():
        print(n, " ".join(f"{r['pw_rms']*100:.2f}/{r['pw_mean']*100:+.2f}" for r in v["regions"]))


if __name__ == "__main__":
    main()
