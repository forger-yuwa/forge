#!/usr/bin/env python3
r"""壁ダンプ時系列から、報告する熱伝達量 (領域別の $h$ 偏差・うねり rms・2 節点振幅) の CSV を作る。

`check_quasisteady.py --series-csv` に渡して、量そのものが頭打ちしているかを判定するための入力。
翼の run は残差が落ちきらない (`NOT CONVERGED (stalled/plateau)`) ので、引用してよいのは
ここで `STEADY` と判定された量だけ (AGENTS.md「準定常確認」)。

列: step, PS, SS_lam, SS_post, all  … 実測に対する $h$ の平均偏差 [%] ($s/S\le0.87$)
    wave_rms  … 負圧面 $s/S$ 0.45–0.95 の (q − 41 点移動平均) の rms [W/m²]
    zig_rms   … 同区間の 2 節点振幅 (q_i − (q_{i-1}+q_{i+1})/2)/2 の rms [W/m²]

usage: python3 case/53.c3x_vane_cht/tools/h_series.py <run_dir> --run run108|run42 [--flux iface_q_eff]
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).parent))
import run_data as rd_mod
from compare_h import arc_map, TABLES, H0, TG
from resid_split import smooth


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--run", required=True)
    ap.add_argument("--flux", default="iface_q_eff"); ap.add_argument("--phys", type=int, default=5)
    a = ap.parse_args(); rd = Path(a.run_dir)
    coarse = {i for i, d in getattr(rd_mod, "RUN42_PARTIAL", {}).items() if a.run == "run42" and d["tol"] > 0.01}
    rows = [r for i, r in enumerate(TABLES[a.run]["rows"]) if r[3] is not None and i not in coarse]
    sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows]) * H0; i0 = int(np.argmin(sd))
    exp = {"PS": (sd[:i0 + 1][::-1], hd[:i0 + 1][::-1]), "SS": (sd[i0:], hd[i0:])}
    fs = sorted(rd.glob(f"res_wall_{a.phys}_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    out = ["step,PS,SS_lam,SS_post,all,wave_rms,zig_rms"]
    geo = None
    for fn in fs:
        step = int(fn.stem.split("_")[-1])
        if step == 0: continue
        with h5py.File(fn) as f:
            v = f["VALUE"]
            if a.flux not in v: sys.exit(f"[h_series] {fn} に {a.flux} が無い")
            q = np.array(v[a.flux], float); Tw = np.array(v["Ts"], float); ok = np.array(v["iface_ok"]) > 0.5
            if geo is None:
                C = np.array(f["MESH/COORD"], float).reshape(-1, 3)[:, :2]; s, ss = arc_map(C)
                he = np.array([np.interp(s[k], *exp["SS" if ss[k] else "PS"]) for k in range(len(C))])
                o = np.argsort(s[ss]); sS = s[ss][o]; mw = (sS >= 0.45) & (sS <= 0.95); geo = True
        h = q / (TG[a.run] - Tw); r = 100 * (h - he) / he
        reg = [ok & ~ss & (s <= 0.87), ok & ss & (s < 0.25), ok & ss & (s >= 0.25) & (s <= 0.87), ok & (s <= 0.87)]
        y = q[ss][o]; w = y - smooth(y, 41); z = np.r_[0, (y[1:-1] - 0.5 * (y[:-2] + y[2:])) / 2, 0]
        out.append(f"{step}," + ",".join(f"{np.nanmean(r[m]):.4f}" for m in reg)
                   + f",{np.sqrt(np.nanmean(w[mw]**2)):.2f},{np.sqrt(np.nanmean(z[mw]**2)):.2f}")
    dst = rd / f"h_series_{a.flux}.csv"; dst.write_text("\n".join(out) + "\n")
    print(f"[h_series] {len(out)-1} dumps -> {dst}"); print("\n".join(out[:1] + out[-3:]))


if __name__ == "__main__":
    main()
