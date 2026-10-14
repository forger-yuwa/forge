"""plan tooling-nozzle-core-grid §4.38: B0 の遅いモードの形と種類 (計算なし、D1 の場から)。
run の res_2000 と res_2032 (level 0 の保存量) の差 ΔQ = Q_2032 − Q_2000 (32 step のあいだの動き) を、もとの状態
(base_res: 同じ格子の level 1 の場、P・T・sonic・Ux・Uy) で線形化して原始量の変化に直す:
  Δρ、Δu = (Δ(ρu) − uΔρ)/ρ、Δv 同様、Δe = Δ(ρE/ρ) − u·Δu (組成は固定とみなす)、
  R = P/(ρT)、γ = c²ρ/P、c_v = R/(γ−1)、ΔT = Δe/c_v、Δp = R(TΔρ + ρΔT)
成分 (無次元): 圧力 a = Δp/(ρc²)、エントロピー s = Δp/(γP) − Δρ/ρ、速度 Δu/c・Δv/c。
縮流部 (x∈[−12.5, −2] r_t) と島 (x∈[−6, −3.25]、壁から 60〜100 節点) で RMS と、x の向き・半径の向きの符号のそろい方 (隣り合う節点の相関) を出す。
図: 縮流部の (x, r) に a・s・Δu/c・Δv/c を描く。
usage (AWS): python3 cg_mode_shape.py <run> <base_res> <tag> → _band_ab/core_grid/mode_<tag>.json と mode_<tag>.png"""
import json, sys
from pathlib import Path
import h5py
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
run, base, tag = sys.argv[1], sys.argv[2], sys.argv[3]
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def L(path, names):
    with h5py.File(path, "r") as h:
        return {n: np.asarray(h["VALUE/" + n][:], float) for n in names}, np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
Q0, xy = L(HERE / run / "res_2000.h5", ["ro", "roUx", "roUy", "roe"])
Q1, _ = L(HERE / run / "res_2032.h5", ["ro", "roUx", "roUy", "roe"])
B, xyb = L(HERE / base, ["ro", "Ux", "Uy", "P", "T", "sonic"])
if not np.array_equal(xy, xyb):
    raise SystemExit("base の格子が run と違う")
n = len(Q0["ro"]); nj = n // ni
rho, u, v, P, T, c = B["ro"], B["Ux"], B["Uy"], B["P"], B["T"], B["sonic"]
dr = Q1["ro"] - Q0["ro"]; du = (Q1["roUx"] - Q0["roUx"] - u * dr) / rho; dv = (Q1["roUy"] - Q0["roUy"] - v * dr) / rho
de = (Q1["roe"] / Q1["ro"] - Q0["roe"] / Q0["ro"]) - (u * du + v * dv)
R = P / (rho * T); gam = c * c * rho / P; cv = R / (gam - 1)
dT = de / cv; dp = R * (T * dr + rho * dT)
comp = {"a_pressure": dp / (rho * c * c), "s_entropy": dp / (gam * P) - dr / rho, "du_over_c": du / c, "dv_over_c": dv / c}
x = xy[:, 0] / S; r = xy[:, 1] / S; jw = nj - 1 - np.tile(np.arange(nj), ni)
cont = (x >= -12.5) & (x <= -2); isl = (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100)
out = {"run": run, "base": base, "ni": ni, "nj": nj, "steps": [2000, 2032], "comp": {}}
for k, f in comp.items():
    F = f.reshape(ni, nj)
    def adjcorr(A, axis):
        a = np.take(A, range(A.shape[axis] - 1), axis=axis); b = np.take(A, range(1, A.shape[axis]), axis=axis)
        return float(np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b)))
    m2 = cont.reshape(ni, nj); cols = np.where(m2.any(1))[0]
    sub = F[cols]                                                     # 縮流部の列
    out["comp"][k] = {"rms_contraction": float(np.sqrt(np.mean(f[cont] ** 2))), "rms_island": float(np.sqrt(np.mean(f[isl] ** 2))),
                      "rms_all": float(np.sqrt(np.mean(f ** 2))), "mean_island": float(np.mean(f[isl])),
                      "frac_positive_island": float(np.mean(f[isl] > 0)),
                      "adjcorr_x_contraction": adjcorr(sub, 0), "adjcorr_r_contraction": adjcorr(sub, 1)}
# 速度の変化の向き (島): 流れの向きとそろうか
um = np.hypot(u, v); ut = (du * u + dv * v) / np.maximum(um, 1e-30); un = (-du * v + dv * u) / np.maximum(um, 1e-30)
out["island_velocity"] = {"rms_dut_over_c": float(np.sqrt(np.mean((ut[isl] / c[isl]) ** 2))), "rms_dun_over_c": float(np.sqrt(np.mean((un[isl] / c[isl]) ** 2))),
                          "mach_island_median": float(np.median(um[isl] / c[isl]))}
(OUTD / f"mode_{tag}.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
for k, v_ in out["comp"].items():
    print(f"{k:12s} RMS 縮流部 {v_['rms_contraction']:.3e} 島 {v_['rms_island']:.3e} 全域 {v_['rms_all']:.3e} | 島の平均 {v_['mean_island']:+.3e} 正の割合 {v_['frac_positive_island']:.2f} | 隣の相関 x {v_['adjcorr_x_contraction']:+.3f} r {v_['adjcorr_r_contraction']:+.3f}")
print("島の速度の変化 (流れの向き/直交)", out["island_velocity"])
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    X = x.reshape(ni, nj); Rr = r.reshape(ni, nj); sel = (X[:, 0] >= -12.5) & (X[:, 0] <= 2)
    for a_, (k, f) in zip(ax, comp.items()):
        F = f.reshape(ni, nj)[sel]; vmax = np.percentile(np.abs(F), 99) or 1.0
        pc = a_.pcolormesh(X[sel], Rr[sel], F, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"); fig.colorbar(pc, ax=a_)
        a_.set_title(f"{tag}: {k} (change over 32 steps)"); a_.set_ylabel("r / r_t")
        a_.plot([-6, -3.25, -3.25, -6, -6], [2.65, 2.65, 4.77, 4.77, 2.65], "k--", lw=0.8)
    ax[-1].set_xlabel("x / r_t"); fig.tight_layout(); fig.savefig(OUTD / f"mode_{tag}.png", dpi=110)
    print(f"[mode] {OUTD / f'mode_{tag}.png'}")
except Exception as e:
    print("図は描けなかった:", e)
