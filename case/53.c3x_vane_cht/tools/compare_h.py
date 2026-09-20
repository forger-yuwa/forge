#!/usr/bin/env python3
r"""V5 段 (a): **実測壁温を課した CFD の熱伝達係数を実測 $h$ と比べる**。

報告 (NASA CR-168015 p.125) の定義に合わせる:
  - $h \equiv q_w / (T_g - T_w)$、$T_g$ は**ガス全温** (run 108 は $T_{T1}$=786 K)。
  - $h$ は $h_0$=1135 W/m²K、$T_w$ は 811 K で正規化して報告されている。

壁熱流束は**壁ダンプの界面診断** (`iface_q_compact` = $k_{eff}(T_1-T_w)/d_1$) を一次に取り、
`iface_q_2nd` (2 次片側差分) を感度として併記する。`qwall` は低 Re 経路のソルバ出力。

usage: python3 case/53.c3x_vane_cht/tools/compare_h.py <run_dir> [--run run108] [--step N]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_data import TABLES                     # noqa: E402

H0, TREF = 1135.0, 811.0
TG = {"run108": 786.0}


def arc_map(P):
    """壁節点列 P (順不同) を弧長順に並べ、LE 起点の正規化弧長と SS/PS 判定を返す。

    `setup_run.wall_profile_csv` と**同じ規約** (LE=x 最小, TE=x 最大, 長い側が SS)。
    壁ダンプは節点順が不定なので、最近傍でループを作り直す。
    """
    n = len(P)
    # 最近傍で閉ループを作る (壁は単一閉曲線)
    used = np.zeros(n, bool)
    order = [int(np.argmin(P[:, 0]))]
    used[order[0]] = True
    for _ in range(n - 1):
        d = np.hypot(*(P - P[order[-1]]).T)
        d[used] = np.inf
        k = int(np.argmin(d))
        order.append(k)
        used[k] = True
    order = np.array(order)
    Q = P[order]
    seg = np.hypot(*np.diff(np.vstack([Q, Q[:1]]), axis=0).T)
    i_le = int(np.argmin(Q[:, 0]))
    i_te = int(np.argmax(Q[:, 0]))
    s = np.zeros(n)
    acc = 0.0
    for k in range(n):
        idx = (i_le + k) % n
        s[idx] = acc
        acc += seg[idx]
    total, s_te = acc, s[i_te]
    fwd = s <= s_te
    arc_f, arc_b = s_te, total - s_te
    ss_is_fwd = arc_f > arc_b
    s_norm = np.where(fwd, s / arc_f, (total - s) / arc_b)
    is_ss = fwd if ss_is_fwd else ~fwd
    inv = np.empty(n, int)
    inv[order] = np.arange(n)
    return s_norm[inv], is_ss[inv]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--phys", type=int, default=5)
    a = ap.parse_args()

    rd = Path(a.run_dir if Path(a.run_dir).is_absolute() else ROOT / "case/53.c3x_vane_cht" / a.run_dir)
    fs = sorted(rd.glob(f"res_wall_{a.phys}_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    if not fs:
        sys.exit(f"no wall dump in {rd}")
    fn = fs[-1] if a.step is None else rd / f"res_wall_{a.phys}_{a.step}.h5"
    with h5py.File(fn, "r") as f:
        C = np.array(f["MESH/COORD"]).reshape(-1, 3)
        v = f["VALUE"]
        Tw = np.array(v["Ts"])
        q = {k: np.array(v[k]) for k in ("iface_q_compact", "iface_q_2nd", "iface_q_recon", "qwall")}
        ok = np.array(v["iface_ok"]) > 0.5
    Tg = TG[a.run]
    s_norm, is_ss = arc_map(C[:, :2])

    rows = TABLES[a.run]["rows"]
    sd = np.array([r[0] for r in rows])
    hd = np.array([r[3] for r in rows]) * H0
    i_stag = int(np.argmin(sd))
    exp = {"PS": (sd[:i_stag + 1][::-1], hd[:i_stag + 1][::-1]), "SS": (sd[i_stag:], hd[i_stag:])}

    print(f"[compare_h] {fn.relative_to(ROOT)}  Tg={Tg} K  nodes={len(C)} (iface_ok {ok.sum()})")
    print(f"{'source':<18}{'h mean':>9}{'h max':>9}  |  vs exp: bias% rms% (points in data range)")
    out = {}
    for key, qq in q.items():
        h = qq / (Tg - Tw)
        out[key] = h
        m = ok & (s_norm <= 0.87)          # 後縁側はデータ範囲外 (README 参照)
        he = np.array([np.interp(s_norm[k], *exp["SS" if is_ss[k] else "PS"]) for k in range(len(C))])
        r = (h[m] - he[m]) / he[m]
        print(f"{key:<18}{h[ok].mean():9.1f}{h[ok].max():9.1f}  |  "
              f"{100*r.mean():+7.1f}{100*np.sqrt((r**2).mean()):8.1f}  (n={m.sum()})")

    np.savez(rd / "h_compare.npz", s_norm=s_norm, is_ss=is_ss, Tw=Tw, ok=ok, Tg=Tg,
             **{k: v for k, v in out.items()})
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5.5))
        for side, sgn, col in (("SS", +1, "tab:red"), ("PS", -1, "tab:blue")):
            se, he_ = exp[side]
            ax.plot(sgn * se, he_ / H0, "o-", color=col, ms=3, lw=1, label=f"exp {side}")
            sel = ok & (is_ss if side == "SS" else ~is_ss)
            o = np.argsort(s_norm[sel])
            ax.plot(sgn * s_norm[sel][o], out["iface_q_compact"][sel][o] / H0, "-", color=col,
                    alpha=0.55, lw=2, label=f"forge {side}")
        ax.set_xlabel("$-s/S$ (PS)   |   $+s/S$ (SS)")
        ax.set_ylabel("$h/h_0$   ($h_0$=1135 W/m²K)")
        ax.grid(alpha=0.3); ax.legend(); ax.set_title(f"C3X {a.run} — h (measured $T_w$ imposed)")
        fig.tight_layout(); fig.savefig(rd / "h_compare.png", dpi=110)
        print(f"[compare_h] -> {(rd/'h_compare.png').relative_to(ROOT)}")
    except Exception as e:                      # 図は補助
        print("[compare_h] plot skipped:", e)


if __name__ == "__main__":
    main()
