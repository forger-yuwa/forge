"""スロート直後の r″ の山は避けられるか (形状のみ, CFD 0 step)。plan verification-m6-axis-wave-mesh-su2 §5.1 #16 の検討材料。
同じ MOC 点群 (n_axis 2400) で同時当てはめの始点条件だけを変える:
  V0: 現行 (r″(0)=1/R=0.5) / V1: r″(0)=κ_H (Hall 壁流線の曲率 0.4735) / V2: r″(0) 自由 (r′(0)=0 は保つ)
  / V3: V1 + スロート直後 (x<x_w) の点の重みを 0 (MOC の始点の約 0.1° のずれを無理に追わない)
出力: x∈[0,0.3] の r″ の最大・最小と r‴ の最大、点上の位置・角度のずれ (x<0.1 と x≥0.1 に分けて)。
usage: design/.venv-opt/bin/python throat_r2_variants.py → _band_ab/throat_r2_variants.json, _band_ab/throat_r2_variants.png
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); plt.rcParams["font.family"] = "Noto Sans CJK JP"
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain, load_problem  # noqa: E402

H0, H1, XG, LAM, SR, ST = 0.0125, 0.5, 6.0, 1e-9, 1e-6, 1e-4


def fit(tb, kappa0, w_cut=0.0, k=5):
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]; x0, xe = x[0], x[-1]
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (XG - x0), 1.0); xs.append(xs[-1] + H0 + (H1 - H0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1]); xi = xi[xi < xe - 0.5 * H1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]; nc = len(t) - k - 1; E = np.eye(nc)
    D = lambda xq, d: np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1); w = np.gradient(x); w = w / w.mean()
    w = np.where((x > x0) & (x < w_cut), 0.0, w)
    xq = np.linspace(x0, xe, 8000); B3 = D(xq, 3)
    A = (B0.T * (w / SR ** 2)) @ B0 + (B1.T * (w / ST ** 2)) @ B1 + LAM * (B3.T * np.gradient(xq)) @ B3 / SR ** 2
    b = B0.T @ (w * r / SR ** 2) + B1.T @ (w * np.tan(th) / ST ** 2)
    rows = [D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]]; dv = [r[0], 0.0, r[-1], np.tan(th[-1])]
    if kappa0 is not None:
        rows.append(D(np.r_[x0], 2)[0]); dv.append(kappa0)
    Cm = np.array(rows); m = len(dv)
    c = np.linalg.solve(np.block([[A, Cm.T], [Cm, np.zeros((m, m))]]), np.r_[b, dv])[:nc]
    return BSpline(t, c, k)


p = load_problem(C / "problem_d155_euler_c2final_n2400.yaml"); d = design_chain(p); tb = d["wall_inv"]; R = d["R"]
KH = 0.4735
V = {"V0 r″(0)=1/R": fit(tb, 1 / R), "V1 r″(0)=κ_H 0.4735": fit(tb, KH), "V2 r″(0) 自由": fit(tb, None), "V3 κ_H + x<0.1 の点を無視": fit(tb, KH, 0.1)}
xg = np.linspace(0, 0.3, 30001); out = {}
x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]; near = x < 0.1
for nm, s in V.items():
    r2 = s(xg, 2); dr = s(x) - r; dth = np.degrees(np.arctan(s(x, 1)) - th)
    out[nm] = dict(r2_at0=float(s(np.r_[0.0], 2)[0]), r2_max_0_0p3=float(r2.max()), r2_min_0_0p1=float(r2[xg <= 0.1].min()),
                   r3_absmax_0_0p3=float(np.abs(s(xg, 3)).max()),
                   r2_overshoot_vs_monotone=float(r2.max() - max(r2[0], r2[-1])),
                   dr_max_xlt0p1=float(np.abs(dr[near]).max()), dth_max_xlt0p1_deg=float(np.abs(dth[near]).max()),
                   dr_max_xge0p1=float(np.abs(dr[~near]).max()), dth_max_xge0p1_deg=float(np.abs(dth[~near]).max()),
                   r2_extrema_count_0_95=int(np.count_nonzero(np.diff(np.sign(np.diff(s(np.linspace(0, x[-1], 200001), 2)))) != 0)))
(C / "_band_ab/throat_r2_variants.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
xx = np.linspace(-0.3, 1.0, 20000)
ax.plot(xx, d["wall"].r(xx, 2), color="C0", lw=1, alpha=.6, label="補間壁 (現行の生産)")
for (nm, s), c in zip(V.items(), ("k", "C3", "C2", "C1")):
    xp = xx[xx >= 0]; ax.plot(xp, s(xp, 2), color=c, lw=1.2, label="当てはめ " + nm)
xu = xx[xx < 0]; ax.plot(xu, d["wall"].r(xu, 2), color="gray", lw=1.2, label="縮流部 Hermite (r″(0)=1/R)")
ax.plot(x[x < 1], np.full((x < 1).sum(), 0.32), "|", color="C3", ms=8, label="MOC 壁点")
ax.set_xlim(-0.3, 1.0); ax.set_ylim(0.3, 0.7); ax.grid(alpha=.3); ax.set_xlabel("x / r_t"); ax.set_ylabel("r″ [1/r_t]")
ax.set_title("スロート付近の r″: 始点条件の違い (MOC 2400 点)"); ax.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(C / "_band_ab/throat_r2_variants.png", dpi=120)
for nm, v in out.items():
    print(nm, {k: round(x_, 6) if isinstance(x_, float) else x_ for k, x_ in v.items()})
