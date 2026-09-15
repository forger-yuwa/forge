"""va3 M4.19 (Euler, node) の入口全温分布 run を uniform run と比較する後処理。
軸線・壁線 (slip 壁ノード)・出口断面 (x_max−2 r_t) の T0/T/M/S/g、質量流量、凝縮 onset、2D の g/S/T0 コンター。
usage: python3 analyze_inletTt_va3.py  → figs/va3_inletTt_*.png, study_va3_inletTt.json
"""
import json, sys
from pathlib import Path
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
sys.path.insert(0, "/home/sano/work/forge/solver_density_cuda/tools")
from total_quantities import total_state

CASE = Path(__file__).resolve().parent
RUNS = [("uniform Tt 1161 K, noneq (run_0126)", "run_0126_va3_M4.19_Lc8_noneq_rerun_merged", "C0", "-"),
        ("inlet Tt profile, dry (run_0128, 24000 step)",   "run_0128_va3_M4.19_Lc8_dry_inletTt_merged", "C2", "-."),
        ("inlet Tt profile, noneq (run_0127, 24000 step)", "run_0127_va3_M4.19_Lc8_noneq_inletTt_merged", "C1", "--")]
RUNS = [r for r in RUNS if (CASE / r[1] / "res_12000.h5").exists() or list((CASE / r[1]).glob("res_[0-9]*.h5"))]


def psat_l(T):
    Tc = np.maximum(T, 120.0)
    return np.exp(54.842763 - 6763.22 / Tc - 4.210 * np.log(Tc) + 0.000367 * Tc
                  + np.tanh(0.0415 * (Tc - 218.8)) * (53.878 - 1331.22 / Tc - 9.44523 * np.log(Tc) + 0.014025 * Tc))


def load(run):
    rd = CASE / run
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    info = json.loads((rd / "prepare_info.json").read_text()); S = float(info["scale_m"])
    with h5py.File(rd / "nozzle.h5") as n:
        nc = n["MESH/COORD"][:].reshape(-1, 3)
        wall = np.unique(n["BCONDS/3/vizBfaceNodes"][:])
        conn = n["VIZMESH/CONNE"][:] if "VIZMESH/CONNE" in n else None
    f = h5py.File(res)
    q = {k: f["VALUE/" + k][:] for k in ("T", "P", "ro", "Ux", "Uy", "sonic", "Y1")}
    q["g"] = f["VALUE/g_0"][:] if "VALUE/g_0" in f else np.zeros_like(q["T"])
    q["M"] = np.hypot(q["Ux"], q["Uy"]) / q["sonic"]
    pv = q["ro"] * (q["Y1"] - q["g"]) * 461.5 * q["T"]
    q["S"] = f["VALUE/condS_0"][:] if "VALUE/condS_0" in f else pv / psat_l(q["T"])
    for k in ("condR30_0", "condDrdt_0", "condTsat_0"):
        if "VALUE/" + k in f: q[k] = f["VALUE/" + k][:]
    st = total_state(rd, res); q["T0"] = st["T0"]; q["P0"] = st["P0"]
    x, r = nc[:, 0] / S, nc[:, 1] / S
    return dict(x=x, r=r, q=q, wall=wall, S=S, res=res.name, x_E=info["x_E"], x_A=info["x_A"])


def line(d, idx):
    o = np.argsort(d["x"][idx]); idx = idx[o]
    return {"x": d["x"][idx], "r": d["r"][idx], **{k: v[idx] for k, v in d["q"].items()}}


def plane(d, x_target):
    cols = np.unique(d["x"]); xc = cols[np.argmin(abs(cols - x_target))]
    idx = np.where(np.abs(d["x"] - xc) < 1e-9)[0]; idx = idx[np.argsort(d["r"][idx])]
    p = {"x": xc, "r": d["r"][idx], **{k: v[idx] for k, v in d["q"].items()}}
    w = p["ro"] * p["Ux"] * p["r"]
    p["mdot"] = float(2 * np.pi * np.trapz(w, p["r"]) * d["S"] ** 2)
    for k in ("M", "T", "T0", "g", "S", "P"):
        p[k + "_avg"] = float(np.trapz(w * p[k], p["r"]) / np.trapz(w, p["r"]))
    return p


def onset(l, thr):
    i = np.where(l["g"] > thr)[0]
    return (float(l["x"][i[0]]), float(l["T"][i[0]]), float(l["S"][i[0]])) if len(i) else (None, None, None)


out = {}; D = {}
for lab, run, c, ls in RUNS:
    d = load(run); D[run] = d
    ax = line(d, np.where(d["r"] < 1e-9)[0]); wl = line(d, d["wall"])
    inlet = plane(d, d["x"].min()); ex2 = plane(d, d["x"].max() - 2.0); exF = plane(d, d["x"].max())
    Yw = float(np.median(d["q"]["Y1"][d["x"] < d["x"].min() + 0.1]))
    thr = 1e-3 * Yw
    rec = dict(run=run, res=d["res"], Y_H2O_inlet=Yw,
               mdot_inlet=inlet["mdot"], mdot_exit=exF["mdot"],
               axis=dict(T_min=float(ax["T"].min()), S_max=float(ax["S"].max()), g_max=float(ax["g"].max()),
                         M_exit=float(ax["M"][-1]), T_exit=float(ax["T"][-1]), T0_exit=float(ax["T0"][-1]),
                         onset=onset(ax, thr)),
               wall=dict(T_min=float(wl["T"].min()), S_max=float(wl["S"].max()), g_max=float(wl["g"].max()),
                         x_g_max=float(wl["x"][np.argmax(wl["g"])]), M_exit=float(wl["M"][-1]),
                         T_exit=float(wl["T"][-1]), T0_exit=float(wl["T0"][-1]), onset=onset(wl, thr)),
               exit_plane_xmax_minus_2rt=dict(x_rt=float(ex2["x"]), M_avg=ex2["M_avg"], T_avg=ex2["T_avg"], T0_avg=ex2["T0_avg"],
                                              g_avg=ex2["g_avg"], g_frac_avg=ex2["g_avg"] / Yw, g_max=float(ex2["g"].max()),
                                              S_max=float(ex2["S"].max()), T0_min=float(ex2["T0"].min()), T0_max=float(ex2["T0"].max()),
                                              M_min=float(ex2["M"].min()), M_max=float(ex2["M"].max())),
               field=dict(g_max=float(d["q"]["g"].max()), g_frac_max=float(d["q"]["g"].max() / Yw),
                          x_g_max=float(d["x"][np.argmax(d["q"]["g"])]), r_g_max=float(d["r"][np.argmax(d["q"]["g"])]),
                          S_max=float(d["q"]["S"].max()), T_min=float(d["q"]["T"].min())))
    out[run] = rec; d["ax"], d["wl"], d["ex2"], d["inlet"] = ax, wl, ex2, inlet
(CASE / "study_va3_inletTt.json").write_text(json.dumps(out, indent=1))
for r, v in out.items():
    print(r, json.dumps({k: v[k] for k in ("mdot_inlet", "mdot_exit")}), "\n  axis", v["axis"], "\n  wall", v["wall"], "\n  exit", v["exit_plane_xmax_minus_2rt"], "\n  field", v["field"])

# ---- 図 1: 線分布 (軸 / 壁 / 出口断面) ----
fig, axs = plt.subplots(3, 4, figsize=(20, 12))
for lab, run, c, ls in RUNS:
    d = D[run]
    for j, (k, yl) in enumerate([("M", "M"), ("T", "T [K]"), ("S", "S (liquid)"), ("g", "g (liquid mass frac.)")]):
        axs[0, j].plot(d["ax"]["x"], d["ax"][k], c, ls=ls, label=lab); axs[0, j].set_ylabel("axis: " + yl)
        axs[1, j].plot(d["wl"]["x"], d["wl"][k], c, ls=ls, label=lab); axs[1, j].set_ylabel("wall: " + yl)
    for j, (k, yl) in enumerate([("T0", "T0 [K]"), ("T", "T [K]"), ("M", "M"), ("g", "g")]):
        axs[2, j].plot(d["ex2"][k], d["ex2"]["r"], c, ls=ls, label=lab); axs[2, j].set_xlabel(f"exit plane (x={d['ex2']['x']:.2f} r_t): " + yl); axs[2, j].set_ylabel("r / r_t")
for a in axs[:2].ravel():
    a.set_xlabel("x / r_t"); a.grid(alpha=.3)
    for xv in (D[RUNS[0][1]]["x_A"], D[RUNS[0][1]]["x_E"]): a.axvline(xv, color="k", lw=.5, ls=":")
for a in axs[2]: a.grid(alpha=.3)
axs[0, 2].axhline(1, color="k", lw=.5); axs[1, 2].axhline(1, color="k", lw=.5)
axs[0, 2].set_yscale("log"); axs[1, 2].set_yscale("log")
axs[0, 0].legend(fontsize=8)
fig.suptitle("va3 M4.19 L_c8 node Euler: uniform Tt vs inlet Tt profile (925-1180 K), noneq condensation")
fig.tight_layout(); fig.savefig(CASE / "figs/va3_inletTt_lines.png", dpi=130); plt.close(fig)

# ---- 図 2: 2D コンター (T0, S, g) : profile noneq と profile dry ----
def tri(d):
    return mtri.Triangulation(d["x"], d["r"])
runs2 = [r for r in RUNS if "inletTt" in r[1]]
fig, axs = plt.subplots(len(runs2) * 3, 1, figsize=(18, 3.2 * 3 * len(runs2)))
i = 0
for lab, run, c, ls in runs2:
    d = D[run]; t = tri(d)
    for k, lv in (("T0", np.linspace(920, 1185, 54)), ("S", None), ("g", None)):
        a = axs[i]; i += 1
        v = d["q"][k]
        if k == "S": v = np.log10(np.maximum(v, 1e-3)); lv = np.linspace(-3, max(0.5, float(v.max())), 36); ttl = "log10 S"
        elif k == "g": lv = np.linspace(0, max(1e-6, float(v.max())), 36); ttl = "g"
        else: ttl = "T0 [K]"
        cf = a.tricontourf(t, v, levels=lv, cmap="viridis" if k != "S" else "coolwarm"); plt.colorbar(cf, ax=a, pad=.01)
        if k == "S": a.tricontour(t, v, levels=[0.0], colors="k", linewidths=.8)
        a.set_aspect("equal"); a.set_title(f"{lab}: {ttl}", fontsize=10); a.set_xlabel("x / r_t"); a.set_ylabel("r / r_t")
fig.tight_layout(); fig.savefig(CASE / "figs/va3_inletTt_contours.png", dpi=110); plt.close(fig)
print("figs written")
