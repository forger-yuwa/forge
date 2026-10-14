#!/usr/bin/env python3
"""③ベル (M) L100 (CFD ピン版): NS (run_0019) と Euler (run_0018) の圧力・Mach の分布を比べ、波の位置がずれているかを確かめる。
2026-10-12、plan §5.1 #7 (コアで合わせる方法の不採用の根拠の確認 — ユーザの条件 1)。軸上 M(x)・p(x) の比較と、局所の窓での
x 方向のずれ (相互相関で、勾配 dM/dx どうしを ±1 r_t の範囲でずらして最も合う Δx)、断面 x 4・10.89・21.51 での M・p の半径分布。
出力: wave_shift_l100.json・wave_shift_axis.png・wave_shift_radial.png (case/68.bell_two_mode/rao_m/)。"""
import json, os, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
NS = HERE.parent / "run_0019_tp_m_ns_l100"; EU = HERE.parent / "run_0018_tp_m_euler_l100_pin"

def load(run, res):
    info = json.loads((run / "prepare_info.json").read_text()); ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]); S = float(info["scale_m"])
    with h5py.File(run / "nozzle.h5") as f: c = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(run / res) as f:
        V = f["VALUE"]; P = V["P"][:]; M = np.hypot(V["Ux"][:], V["Uy"][:]) / V["sonic"][:]
    return dict(x=c[:, 0].reshape(ni, nj) / S, r=c[:, 1].reshape(ni, nj) / S, P=P.reshape(ni, nj), M=M.reshape(ni, nj))

def main():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in font_manager.findSystemFonts([os.path.expanduser("~/.fonts")]): font_manager.fontManager.addfont(f)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
    n = load(NS, "res_60000.h5"); e = load(EU, "res_30000.h5")
    # 軸 (j = 0) の分布: 両方とも軸の節点は x 一定の列なので x[:,0] で比べる
    xn, xe = n["x"][:, 0], e["x"][:, 0]
    Mn, Me = n["M"][:, 0], np.interp(xn, xe, e["M"][:, 0]); Pn, Pe = n["P"][:, 0], np.interp(xn, xe, e["P"][:, 0])
    m = xn >= 0
    # 局所のずれ: 窓 [xc − 2, xc + 2] で dM/dx (NS) と dM/dx (Euler を Δx ずらしたもの) の相関が最大の Δx
    xg = np.linspace(0, xn.max(), 4001); dMn = np.gradient(np.interp(xg, xn, Mn), xg); dMe_f = lambda s: np.gradient(np.interp(xg - s, xe, e["M"][:, 0]), xg)
    shifts = np.linspace(-1.0, 1.0, 401); cen = np.arange(2.0, xn.max() - 2.0 + 1e-9, 1.0); rows = []
    for xc in cen:
        w = (xg >= xc - 2) & (xg <= xc + 2)
        best = max(((float(np.corrcoef(dMn[w], dMe_f(s)[w])[0, 1]), s) for s in shifts), key=lambda t: t[0])
        rows.append({"x_center": float(xc), "shift_x": float(best[1]), "corr": float(best[0]), "corr_at_0": float(np.corrcoef(dMn[w], dMe_f(0.0)[w])[0, 1])})
    out = {"ns": str(NS.relative_to(ROOT)), "euler": str(EU.relative_to(ROOT)),
           "axis_M_maxdiff": float(np.max(np.abs(Mn[m] - Me[m]))), "axis_P_maxreldiff": float(np.max(np.abs(Pn[m] / Pe[m] - 1))),
           "axis_local_shift": rows}
    # 半径分布
    st = {}
    for xs in (3.99, 10.89, 21.51):
        i = int(np.argmin(abs(xn - xs))); xi = xn[i]
        # NS の列は壁近くで x がずれるので、各行で x = xi に補間
        rN = np.array([np.interp(xi, n["x"][:, j], n["r"][:, j]) for j in range(n["x"].shape[1])])
        MN = np.array([np.interp(xi, n["x"][:, j], n["M"][:, j]) for j in range(n["x"].shape[1])])
        PN = np.array([np.interp(xi, n["x"][:, j], n["P"][:, j]) for j in range(n["x"].shape[1])])
        k = int(np.clip(np.searchsorted(xe, xi), 1, len(xe) - 1)); w = (xi - xe[k - 1]) / (xe[k] - xe[k - 1])
        rE = (1 - w) * e["r"][k - 1] + w * e["r"][k]; ME = (1 - w) * e["M"][k - 1] + w * e["M"][k]; PE = (1 - w) * e["P"][k - 1] + w * e["P"][k]
        R = rN[-1]
        # 半径方向の局所のずれ: dM/dr の最大の位置 (波の位置) を NS と Euler で
        rg = np.linspace(0, 0.9 * R, 2001); gN = np.gradient(np.interp(rg, rN, MN), rg); gE = np.gradient(np.interp(rg, rE, ME), rg)
        jN, jE = int(np.argmax(np.abs(gN))), int(np.argmax(np.abs(gE)))
        st[f"{xi:.2f}"] = {"R": float(R), "r_maxgrad_ns": float(rg[jN] / R), "r_maxgrad_euler": float(rg[jE] / R),
                           "maxgrad_ns": float(gN[jN]), "maxgrad_euler": float(gE[jE]),
                           "M_maxdiff_core": float(np.max(np.abs(np.interp(rg, rN, MN) - np.interp(rg, rE, ME))))}
        st[f"{xi:.2f}"]["_prof"] = (rN / R, MN, PN, rE / R, ME, PE)
    out["radial"] = {k: {kk: vv for kk, vv in v.items() if kk != "_prof"} for k, v in st.items()}
    (HERE / "wave_shift_l100.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    # 図: 軸
    fig, ax = plt.subplots(3, 1, figsize=(10, 9.5), sharex=True)
    ax[0].plot(xn[m], Mn[m], color="#222", lw=1.4, label="NS"); ax[0].plot(xn[m], Me[m], color="#3a8f5c", lw=1.2, ls="--", label="Euler (同じ壁)")
    ax[0].set_ylabel("軸上の Mach"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(xn[m], Mn[m] - Me[m], color="#c9721f", lw=1.2); ax[1].axhline(0, color="0.5", lw=0.8); ax[1].set_ylabel("Mach の差 (NS − Euler)"); ax[1].grid(alpha=0.3)
    xc = [r_["x_center"] for r_ in rows]; ax[2].plot(xc, [r_["shift_x"] for r_ in rows], "o-", color="#2f6db5", ms=3, label="窓 ±2 r_t で dM/dx が最も合う Euler のずらし量 Δx")
    ax2 = ax[2].twinx(); ax2.plot(xc, [r_["corr"] for r_ in rows], color="0.5", lw=0.8, label="そのときの相関"); ax2.plot(xc, [r_["corr_at_0"] for r_ in rows], color="0.5", lw=0.8, ls=":", label="ずらさないときの相関"); ax2.set_ylim(0, 1.05); ax2.set_ylabel("相関")
    ax[2].axhline(0, color="0.5", lw=0.8); ax[2].set_ylabel("Δx [r_t] (正 = Euler を下流へ)"); ax[2].set_xlabel("x [r_t]"); ax[2].grid(alpha=0.3)
    h1, l1 = ax[2].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax[2].legend(h1 + h2, l1 + l2, fontsize=8.5, loc="lower left")
    fig.suptitle("軸上の Mach: NS と Euler (L100、CFD ピン版の同じ壁)", fontsize=12); fig.tight_layout(); fig.savefig(HERE / "wave_shift_axis.png", dpi=130); plt.close(fig)
    # 図: 半径
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for c_, (k, v) in enumerate(st.items()):
        rn, MN, PN, re, ME, PE = v["_prof"]
        ax[0, c_].plot(rn, MN, color="#222", lw=1.4, label="NS"); ax[0, c_].plot(re, ME, color="#3a8f5c", lw=1.2, ls="--", label="Euler")
        ax[0, c_].set_title(f"x = {k} r_t", fontsize=11); ax[0, c_].grid(alpha=0.3); ax[0, c_].legend(fontsize=8)
        ax[1, c_].plot(rn, PN / 1e3, color="#222", lw=1.4); ax[1, c_].plot(re, PE / 1e3, color="#3a8f5c", lw=1.2, ls="--"); ax[1, c_].grid(alpha=0.3); ax[1, c_].set_xlabel("r / R (軸 0、壁 1)")
        for a_ in ax[:, c_]:
            a_.axvline(v["r_maxgrad_ns"], color="#222", lw=0.8, ls=":"); a_.axvline(v["r_maxgrad_euler"], color="#3a8f5c", lw=0.8, ls=":")
    ax[0, 0].set_ylabel("Mach"); ax[1, 0].set_ylabel("静圧 [kPa]")
    fig.suptitle("断面の Mach と静圧の半径分布 (点線 = |dM/dr| が最大の位置、黒 NS・緑 Euler)", fontsize=12); fig.tight_layout(); fig.savefig(HERE / "wave_shift_radial.png", dpi=125); plt.close(fig)
    print(json.dumps({k: v for k, v in out.items() if k != "axis_local_shift"}, ensure_ascii=False, indent=1))
    print("local shifts:", [(r_["x_center"], round(r_["shift_x"], 3), round(r_["corr"], 3), round(r_["corr_at_0"], 3)) for r_ in rows])

if __name__ == "__main__":
    main()
