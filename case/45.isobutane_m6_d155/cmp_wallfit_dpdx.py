"""補間壁 (A) と同時当てはめ壁 (B) の Euler 解で、無次元圧力勾配 (r_t/p)∂p/∂x を比べる (plan verification-m6 §5.1 #15 の補足)。
各腕 3 回の延長 run (run_0059〜0064) の最終場で、軸 (η=0)・η=0.1・壁の 1 つ内側の列に沿って比べ、B−A と腕内レンジを出す。図 _band_ab/wallfit_dpdx.png。
usage: python3 cmp_wallfit_dpdx.py
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager; font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402
A = ["run_0059_euler_wallfit_interp_r1_ext6k", "run_0060_euler_wallfit_interp_r2_ext6k", "run_0061_euler_wallfit_interp_r3_ext6k"]
B = ["run_0062_euler_wallfit_fit_r1_ext6k", "run_0063_euler_wallfit_fit_r2_ext6k", "run_0064_euler_wallfit_fit_r3_ext6k"]
xq = np.arange(0.0, 95.2, 0.02)


def lines(run):
    F = load_field(C / run, "res_6000.h5")
    out = {f"eta{e}": eta_line(F, "dpdx_nd", e, xq) for e in (0.0, 0.1, 0.5)}
    X, D = F["X"][:, -2], F["V"]["dpdx_nd"][:, -2]
    out["wall_row"] = np.interp(xq, X, D)
    return out


LA = [lines(r) for r in A]; LB = [lines(r) for r in B]
res = {}
for k in LA[0]:
    a = np.array([l[k] for l in LA]); b = np.array([l[k] for l in LB]); d = b.mean(0) - a.mean(0)
    rng = np.maximum(np.ptp(a, 0), np.ptp(b, 0))
    res[k] = {}
    for lo, hi in ((0, 2), (2, 10), (10, 42), (42, 94)):
        m = (xq >= lo) & (xq < hi)
        res[k][f"[{lo},{hi})"] = dict(dpdx_A_range=[float(a.mean(0)[m].min()), float(a.mean(0)[m].max())],
                                      B_minus_A_maxabs=float(np.abs(d[m]).max()), x_at=float(xq[m][np.argmax(np.abs(d[m]))]),
                                      repeat_range_max=float(rng[m].max()))
(C / "_band_ab/wallfit_dpdx.json").write_text(json.dumps(res, indent=1))
fig, ax = plt.subplots(4, 1, figsize=(11, 12))
for i, (k, ttl) in enumerate((("eta0.0", "軸 (η=0)"), ("eta0.1", "η=0.1"), ("wall_row", "壁の 1 つ内側の列"))):
    ax[i].plot(xq, np.mean([l[k] for l in LA], 0), "C0", lw=1, label="補間壁 (3 回平均)")
    ax[i].plot(xq, np.mean([l[k] for l in LB], 0), "C3--", lw=1, label="当てはめ壁 (3 回平均)")
    ax[i].set_title(f"(r_t/p)∂p/∂x — {ttl}"); ax[i].grid(alpha=.3); ax[i].set_xlim(0, 95)
    ax[i].set_ylim(np.percentile(np.mean([l[k] for l in LA], 0), [0.5, 99.5]) * np.array([1.1, 1.1]))
for k, c in (("eta0.0", "k"), ("eta0.1", "C2"), ("wall_row", "C1")):
    ax[3].plot(xq, np.mean([l[k] for l in LB], 0) - np.mean([l[k] for l in LA], 0), c, lw=.8, label=k)
ax[3].set_title("差 (当てはめ − 補間)"); ax[3].grid(alpha=.3); ax[3].set_xlim(0, 95); ax[3].set_xlabel("x / r_t")
ax[0].legend(); ax[3].legend()
for a_ in ax:
    a_.set_ylabel("(r_t/p)∂p/∂x")

fig.tight_layout(); fig.savefig(C / "_band_ab/wallfit_dpdx.png", dpi=110)
for k, v in res.items():
    print(k, json.dumps(v))
