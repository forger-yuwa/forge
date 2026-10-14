import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import glob
for fp in glob.glob(str(__import__("pathlib").Path.home() / ".fonts" / "*CJK*")):
    fm.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Noto Sans CJK JP"]
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
R = "remote/"
r = json.load(open(R + "m55_result.json"))
rt = r["r_t_m"]; x_in, x_e = r["x_in"], r["x_e"]
d = np.genfromtxt(R + "m55_y1req.csv", delimiter=",", names=True)
x, rw, y1req, ftab = d["x_rt"], d["r_w_rt"], d["y1_req_m"], d["first_frac_table"]
def blendA(xs, frac=2.4e-6, ft=7.4e-7, b0=1.0, b1=17.0):
    t = np.clip((xs - b0) / (b1 - b0), 0, 1); t = t * t * (3 - 2 * t); return ft + (frac - ft) * t
def stations(ni, refine, width, dtab=None):
    xs = np.linspace(x_in, x_e, 8001)
    dens = 1 + (refine - 1) * np.exp(-(xs / width) ** 2)
    if dtab is not None:
        dens = dens * np.interp(xs, dtab[:, 0], dtab[:, 1])
    cum = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xs))]); cum /= cum[-1]
    return np.interp(np.linspace(0, 1, ni), cum, xs)
xa = stations(r["A"]["ni"], 4, 3); xb = stations(r["B"]["ni"], 4, 3)
xc = np.load(R + "m55_C_local_stations.npz")["xs"]
fig, ax = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True)
m = x > x_in + 0.02
ax[0].semilogy(x[m], y1req[m] * 1e6, "k-", lw=1.6, label="要求 y₁,req (y₁⁺ 0.8、estimate_mesh の式)")
ax[0].semilogy(x, blendA(x) * rw * rt * 1e6, "C0--", lw=1.4, label="(A) 2 値 + ブレンド (2.4e−6 / 7.4e−7、上流は全域スロート値)")
ax[0].semilogy(x, ftab * rw * rt * 1e6, "C3-", lw=1.0, label="(B)(C) 第 1 層の表 (要求の下包絡)")
ax[0].set_ylabel("第 1 層 y₁ [µm]"); ax[0].legend(fontsize=9, loc="lower right"); ax[0].grid(alpha=0.3, which="both")
ax[0].set_title(f"case/71 M5.5 (r_t {rt * 1e3:.2f} mm): 壁の第 1 層と流れ方向の間隔 (AR ≤ 0.97×5000)")
for xs, lab, c in ((xa, f"(A) ni {r['A']['ni']}", "C0"), (xb, f"(B) ni {r['B']['ni']}", "C2"), (xc, f"(C) ni {r['C_local']['ni']}", "C3")):
    ax[1].semilogy(0.5 * (xs[1:] + xs[:-1]), np.diff(xs), color=c, lw=1.3, label=lab)
href = np.array(r["C_local"]["x_density_table"])
xx = np.linspace(x_in, x_e, 4000)
dens = 1 + 3 * np.exp(-(xx / 3) ** 2); I = float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(xx))); ax[1].semilogy(xx, I / (1999 * dens), "k:", lw=1.2, label="上限: 雛形 ni 2000 の間隔")
ax[1].semilogy(x[m], 0.97 * 5000 * ftab[m] * rw[m] / np.sqrt(1 + np.gradient(rw, x)[m] ** 2), color="0.5", lw=0.8, label="上限: 4850·y₁·cosθ (表)")
ax[1].set_ylabel("Δx [r_t]"); ax[1].set_xlabel("x [r_t] (設計スロート 0)"); ax[1].legend(fontsize=9, ncol=2); ax[1].grid(alpha=0.3, which="both")
for a_ in ax:
    a_.axvline(-r["L_U"], color="0.7", lw=0.6); a_.axvline(0, color="0.7", lw=0.6)
fig.tight_layout(); fig.savefig("m55_mesh_est.png", dpi=130)
print("ok")
