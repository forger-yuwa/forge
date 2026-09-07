"""[非推奨 2026-09-04] 旧 v3 継ぎはぎ (ρU_x 最大縁 + x<x_lo は相関×比)。生産経路は
feedback/deltastar_loop.py (固定 Euler 基準・帯局所抽出) に置換済み — plans/active/tooling-nozzle-deltastar-core-matched-euler.md。
記録用に残す。

v3 用 δ*(x) の構成 (case/44): NS run から x≥x_meas_lo で δ* を抽出 (x<8 も測れる — コアの半径方向
一様性が十分な本 case では ok)、x<x_meas_lo は相関を測定端の比でスケールして接続、全域を平滑化して
x∈[-0.8, x_F] を覆う CSV を出す (prepare_ns の --dstar-csv + --dstar-blend -1,-0.5 で全域採用)。
使い方: design/.venv-opt/bin/python build_dstar_v3_va.py <run_dir> [--x-lo 3.0] [--dx 0.25]
出力: <run_dir>/dstar_v3.csv, <run_dir>/dstar_v3.png"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
try: _fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.metrics.deltastar import deltastar_from_run
from scipy.interpolate import make_smoothing_spline
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--x-lo", type=float, default=3.0)
ap.add_argument("--dx", type=float, default=0.25); ap.add_argument("--lam", type=float, default=1e-2)
ap.add_argument("--x-hi-trust", type=float, default=None, help="この x より下流は使わず、[x_hi_trust-6, x_hi_trust] の線形フィット勾配で外挿 (出口 BC の影響で末端の δ* が急増するのを避ける)")
ap.add_argument("--out", default="dstar_v3.csv")
a = ap.parse_args(); rd = Path(a.run_dir)
info = json.loads((rd / "prepare_info.json").read_text()); S = float(info["scale_m"])
res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
wp = np.loadtxt(rd / "wall_physical.csv", delimiter=",", skiprows=1) / S
wi = np.loadtxt(rd / "wall_design.csv", delimiter=",", skiprows=1) / [S, S, 1, 1]
xs = np.arange(a.x_lo, wp[-1, 0] - 0.5, a.dx)
d = deltastar_from_run(rd / "nozzle.h5", res, wp, S, xs)
ok = np.asarray(d["ok"]) & np.isfinite(d["dstar"])
# この run が使った δ*(x) (= phys − inv の法線距離) を「前回入力」として比を取る
x_full = np.concatenate([np.linspace(-0.8, wi[0, 0], 40, endpoint=False), wi[:, 0]])
ds_prev_tbl = (np.interp(wi[:, 0], wp[:, 0], wp[:, 1]) - wi[:, 1]) * np.cos(wi[:, 2])
# 上流 (x<0) は前回入力の関数形が無いので相関 (PhysicalNozzleWall と同じ) を作る
from forge_design.evaluate.runner_axismach import design_chain, _gam_or_gas
from forge_design.probdef import load_problem
from forge_design.geometry.wall_axismach import PhysicalNozzleWall
prob = load_problem(Path("/home/sano/work/forge/case/44.vitiated_air_wt/problem_va_R2_LU6_Lc8_ns.yaml"))
dc = design_chain(prob)
base = PhysicalNozzleWall(dc["wall"], dc["wall_inv"], S, float(prob.spec["Pt"]), float(prob.spec["Tt"]), _gam_or_gas(prob), prob.cp)
corr = base._dstar_hist(x_full)
spl = make_smoothing_spline(xs[ok], d["dstar"][ok], lam=a.lam)
x_lo, x_hi = xs[ok][0], xs[ok][-1]
slope_hi, val_hi = float(spl.derivative()(x_hi)), float(spl(x_hi))
if a.x_hi_trust is not None and a.x_hi_trust < x_hi:
    x_hi = float(a.x_hi_trust); mfit = ok & (xs >= x_hi - 6.0) & (xs <= x_hi)
    slope_hi, icpt = np.polyfit(xs[mfit], d["dstar"][mfit], 1); val_hi = slope_hi * x_hi + icpt
    print(f"x_hi_trust={x_hi}: linear-fit slope over [{x_hi-6},{x_hi}] = {slope_hi:.5f} (spline end slope {float(spl.derivative()(xs[ok][-1])):.5f})")
ratio_lo = float(spl(x_lo)) / float(base._dstar_hist(x_lo))
ds = np.empty_like(x_full)
m_in = (x_full >= x_lo) & (x_full <= x_hi); ds[m_in] = spl(x_full[m_in])
m_hi = x_full > x_hi; ds[m_hi] = val_hi + slope_hi * (x_full[m_hi] - x_hi)
m_lo = x_full < x_lo; ds[m_lo] = corr[m_lo] * ratio_lo        # 測定端の比で相関をスケール (連続接続)
ds = np.clip(make_smoothing_spline(x_full, ds, lam=1e-3)(x_full), 0.0, None)
np.savetxt(rd / a.out, np.c_[x_full, ds], delimiter=",",
           header=f"x_rt,dstar_rt (CFD x in [{x_lo},{x_hi}] / corr*{ratio_lo:.3f} x<{x_lo} / lin-extrap x>{x_hi})", comments="")
print(f"measured {ok.sum()}/{len(xs)}; ratio_lo(x={x_lo})={ratio_lo:.3f}; throat dstar new={ds[np.argmin(abs(x_full))]:.5f} (corr {float(base._dstar_hist(0.0)):.5f}, prev-run {ds_prev_tbl[0]:.5f})")
prev_x = wi[:, 0]
r_prev = np.interp(prev_x, x_full, ds) / np.maximum(ds_prev_tbl, 1e-9)
print(f"new/prev-input ratio over design wall: median {np.median(r_prev):.3f} range [{r_prev.min():.3f}, {r_prev.max():.3f}]")
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(xs[ok], d["dstar"][ok], "o", ms=3, label="δ* CFD 抽出 (x≥%g)" % a.x_lo)
ax.plot(x_full, corr, "-", label="δ* 相関 (v1)"); ax.plot(prev_x, ds_prev_tbl, ":", label="この run の入力 δ*")
ax.plot(x_full, ds, "--", label="次 run へ渡す δ* (v3)"); ax.set_xlabel("x / r_t"); ax.set_ylabel("δ* / r_t"); ax.legend(); ax.grid(alpha=.3)
ax.set_title(f"{rd.name}: δ* (new/prev median {np.median(r_prev):.3f})"); fig.tight_layout(); fig.savefig(rd / (Path(a.out).stem + ".png"), dpi=110)
print("saved", rd / a.out)
