#!/usr/bin/env python3
r"""壁節点のエネルギー残差の内訳 (対流 / ソース / 粘性) を弧長に沿って分解する。

`output.interfaceDiag: 1` の壁ダンプにある
  ifaceRconv (対流流束の直後の res_roe) / ifaceRpre (粘性流束の直前) / ifaceRraw (壁ピンの直前) / ifaceFw (壁半割面の寄与)
から、面積あたり [W/m²]・固体向き正で
  conv = Rconv/A,  src = (Rpre-Rconv)/A,  visc_int = (Rraw-Rpre-Fw)/A,  wall = -Fw/A,  q_eff = (Rraw-Fw)/A
を作り、負圧面の区間ごとに「滑らかな成分を引いた残り」の rms を出す
(plan boundary-conjugate-heat-transfer §5.1 #57/#58)。

usage: python3 case/53.c3x_vane_cht/tools/resid_split.py <run_dir> [--step N] [--out fig.png]
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).parent))
from compare_h import arc_map


def smooth(y, w):
    k = np.ones(w) / w
    yp = np.r_[y[w:0:-1], y, y[-2:-w - 2:-1]]
    return np.convolve(yp, k, "same")[w:-w]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--step", type=int); ap.add_argument("--phys", type=int, default=5)
    ap.add_argument("--out"); ap.add_argument("--win", type=int, default=41)
    a = ap.parse_args()
    rd = Path(a.run)
    fs = sorted(rd.glob(f"res_wall_{a.phys}_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    fn = fs[-1] if a.step is None else rd / f"res_wall_{a.phys}_{a.step}.h5"
    with h5py.File(fn, "r") as f:
        C = np.array(f["MESH/COORD"]).reshape(-1, 3)
        v = {k: np.array(f["VALUE"][k], dtype=float) for k in
             ("ifaceRconv", "ifaceRpre", "ifaceRraw", "ifaceFw", "iface_q_eff", "iface_q_compact")}
    s, ss = arc_map(C[:, :2])
    A = (v["ifaceRraw"] - v["ifaceFw"]) / v["iface_q_eff"]          # 面積 (q_eff の定義から逆算)
    comp = {"conv": v["ifaceRconv"] / A, "src": (v["ifaceRpre"] - v["ifaceRconv"]) / A,
            "visc_int": (v["ifaceRraw"] - v["ifaceRpre"] - v["ifaceFw"]) / A,
            "wall(-Fw)": -v["ifaceFw"] / A, "q_eff": v["iface_q_eff"], "q_compact": v["iface_q_compact"]}
    o = np.argsort(s[ss]); sS = s[ss][o]
    print(f"[resid_split] {fn}  suction side n={ss.sum()}")
    print(f"{'component':<12}" + "".join(f"{r:>22}" for r in ("0.30-0.40 (shock)", "0.45-0.95 (waves)")) + "   [kW/m2: mean / rms of (y - smooth)]")
    for k, y in comp.items():
        y = y[ss][o]; d = y - smooth(y, a.win)
        row = ""
        for lo, hi in ((0.30, 0.40), (0.45, 0.95)):
            m = (sS >= lo) & (sS <= hi)
            row += f"{y[m].mean()/1e3:12.2f} /{np.sqrt((d[m]**2).mean())/1e3:7.3f}"
        print(f"{k:<12}{row}")
    # 成分間の相関 (うねり区間)
    m = (sS >= 0.45) & (sS <= 0.95)
    D = {k: (comp[k][ss][o] - smooth(comp[k][ss][o], a.win))[m] for k in comp}
    print("corr of detrended q_eff with: " + "  ".join(f"{k} {np.corrcoef(D['q_eff'], D[k])[0,1]:+.2f}" for k in D if k != "q_eff"))
    if a.out:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        for k in ("q_eff", "wall(-Fw)", "q_compact"):
            ax[0].plot(sS, comp[k][ss][o] / 1e3, lw=1, label=k)
        for k in ("conv", "src", "visc_int"):
            ax[1].plot(sS, comp[k][ss][o] / 1e3, lw=1, label=k)
            ax[2].plot(sS, (comp[k][ss][o] - smooth(comp[k][ss][o], a.win)) / 1e3, lw=1, label=k + " - smooth")
        ax[2].plot(sS, (comp["q_eff"][ss][o] - smooth(comp["q_eff"][ss][o], a.win)) / 1e3, "k", lw=1, label="q_eff - smooth")
        for x in ax: x.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8); x.grid(alpha=.3); x.set_ylabel("kW/m$^2$")
        ax[2].set_xlabel("s/S (suction side)"); ax[1].set_ylim(-15, 15); ax[2].set_ylim(-6, 6); ax[0].set_xlim(0.2, 1.0)
        fig.tight_layout(); fig.savefig(a.out, dpi=110)


if __name__ == "__main__":
    main()
