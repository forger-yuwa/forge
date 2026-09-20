#!/usr/bin/env python3
r"""forge / SU2 / 実測 の**壁面熱伝達係数**を突き合わせる (plan §5.1 #28)。

目的は「負圧面前縁の層流域で $h$ が +40〜75 % 過大」が **forge 固有か、遷移モデルを持たない
低 Re SST に共通か**の切り分け。したがって **forge と SU2 は同一 `.geo`・同一 BC・
一様壁温 566 K** で回す (実測分布を課すと壁温の内挿差が交絡する)。

usage: python3 case/53.c3x_vane_cht/tools/compare_su2.py <forge_run_dir> <su2_dir> [--run run108]
"""
import argparse, csv, sys
from pathlib import Path
import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_data import TABLES          # noqa: E402
from compare_h import arc_map, H0    # noqa: E402

TG = {"run108": 786.0, "run42": 788.0}


def su2_wall(d: Path):
    """SU2 の `surface_flow.csv` から壁面座標と熱流束を読む (固体向き正に揃える)。"""
    f = d / "surface_flow.csv"
    if not f.exists():
        cand = sorted(d.glob("surface_flow*.csv"))
        if not cand:
            sys.exit(f"{d}: surface_flow.csv が無い")
        f = cand[-1]
    rows = list(csv.DictReader(open(f)))
    key = {k.strip().strip('"'): k for k in rows[0]}
    def col(*names):
        for n in names:
            if n in key:
                return np.array([float(r[key[n]]) for r in rows])
        sys.exit(f"{f}: 列 {names} が無い ({list(key)[:12]})")
    x = col("x"); y = col("y")
    q = col("Heat_Flux", "Heat_Flux_(W/m^2)")
    T = col("Temperature") if "Temperature" in key else None
    return np.column_stack([x, y]), q, T, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("forge_run"); ap.add_argument("su2_dir")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--Tw", type=float, default=566.0, help="一様壁温 [K] (両者で同じ)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    Tg = TG[a.run]

    # ---- forge ----
    rd = Path(a.forge_run if Path(a.forge_run).is_absolute() else ROOT / a.forge_run)
    fs = sorted(rd.glob("res_wall_5_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    with h5py.File(fs[-1], "r") as f:
        Cf = np.array(f["MESH/COORD"]).reshape(-1, 3)
        qf = np.array(f["VALUE"]["iface_q_compact"])
    sf, ssf = arc_map(Cf[:, :2])
    hf = qf / (Tg - a.Tw)

    # ---- SU2 ----
    Cs, qs, Ts, csvf = su2_wall(Path(a.su2_dir if Path(a.su2_dir).is_absolute() else ROOT / a.su2_dir))
    ss_, sss = arc_map(Cs)
    # SU2 の Heat_Flux は壁から流体へ出る向きが正 (冷却壁では負)。固体向き正に揃える。
    if np.mean(qs) < 0:
        qs = -qs
    hs = qs / (Tg - a.Tw)

    # ---- 実測 ----
    rows = [r for r in TABLES[a.run]["rows"] if r[3] is not None]
    sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows]) * H0
    i = int(np.argmin(sd))
    exp = {"PS": (sd[:i + 1][::-1], hd[:i + 1][::-1]), "SS": (sd[i:], hd[i:])}

    print(f"[compare_su2] forge {fs[-1].relative_to(ROOT)} ({len(Cf)} nodes)")
    print(f"              SU2   {csvf.relative_to(ROOT)} ({len(Cs)} nodes)")
    print(f"              Tw = {a.Tw} K (両者一様), Tg = {Tg} K\n")
    print(f"{'region':<26}{'forge vs exp':>18}{'SU2 vs exp':>16}{'forge vs SU2':>16}")
    stats = {}
    for name, sel in (("PS", lambda s, t: ~t & (s <= 0.87)),
                      ("SS laminar (s/S<0.25)", lambda s, t: t & (s < 0.25)),
                      ("SS post-transition", lambda s, t: t & (s >= 0.25) & (s <= 0.87)),
                      ("all (s/S<=0.87)", lambda s, t: s <= 0.87)):
        mf = sel(sf, ssf); ms = sel(ss_, sss)
        def vs_exp(s, h, m):
            he = np.array([np.interp(s[k], *exp["SS" if (k in np.where(sss if h is hs else ssf)[0]) else "PS"])
                           for k in range(len(s))])
            return None
        # 実測点上で比べる (点数を揃える)
        side_of = lambda t: "SS" if t else "PS"
        def on_exp(s, tt, h, want_ss):
            se, he = exp["SS" if want_ss else "PS"]
            m = (tt == want_ss) & (s <= 0.87)
            o = np.argsort(s[m])
            hi = np.interp(se[se <= 0.87], s[m][o], h[m][o])
            return hi, he[se <= 0.87], se[se <= 0.87]
        res = []
        for want_ss, lo, hi_ in (((name != "PS"), 0.0, 1.0),):
            pass
        # 領域ごとに forge/SU2 を実測 station に載せて比較
        want_ss = (name != "PS")
        hff, hee, sse = on_exp(sf, ssf, hf, want_ss)
        hsu, _, _ = on_exp(ss_, sss, hs, want_ss)
        if name == "SS laminar (s/S<0.25)":
            k = sse < 0.25
        elif name == "SS post-transition":
            k = sse >= 0.25
        elif name == "PS":
            k = np.ones(len(sse), bool)
        else:
            hff2, hee2, sse2 = on_exp(sf, ssf, hf, False)
            hsu2, _, _ = on_exp(ss_, sss, hs, False)
            hff = np.r_[hff, hff2]; hsu = np.r_[hsu, hsu2]; hee = np.r_[hee, hee2]
            k = np.ones(len(hee), bool)
        def bias(u, v):
            r = (u[k] - v[k]) / v[k]
            return f"{100*r.mean():+6.1f}/{100*np.sqrt((r**2).mean()):5.1f}%"
        stats[name] = (bias(hff, hee), bias(hsu, hee), bias(hff, hsu))
        print(f"{name:<26}{stats[name][0]:>18}{stats[name][1]:>16}{stats[name][2]:>16}")

    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 5.8))
        for side, sgn, col in (("SS", +1, "tab:red"), ("PS", -1, "tab:blue")):
            se, he = exp[side]
            ax.plot(sgn * se, he / H0, "o", color=col, ms=4, label="exp " + side, zorder=5)
        for s_, t_, h_, col, lab in ((sf, ssf, hf, "k", "forge"), (ss_, sss, hs, "tab:green", "SU2")):
            for side, sgn in (("SS", +1), ("PS", -1)):
                m = (t_ if side == "SS" else ~t_)
                o = np.argsort(s_[m])
                ax.plot(sgn * s_[m][o], (h_[m] / H0)[o], "-", color=col, lw=1.8,
                        label=(lab if side == "SS" else None), alpha=0.85)
        ax.axvspan(0, 0.25, color="0.9", zorder=0)
        ax.set_xlabel("$-s/S$ (PS)   |   $+s/S$ (SS)"); ax.set_ylabel("$h/h_0$")
        ax.set_title(f"{TABLES[a.run]['vane']} {a.run} — forge vs SU2 vs experiment "
                     f"(uniform $T_w$={a.Tw:.0f} K)", fontsize=10)
        ax.grid(alpha=0.3); ax.legend(fontsize=9)
        out = a.out or str(rd / "h_forge_su2.png")
        fig.tight_layout(); fig.savefig(out, dpi=110); print(f"\n-> {out}")
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
