#!/usr/bin/env python3
r"""**角帯を除いた熱回路係数**と、帯幅への感度を出す (plan §4.8.6.4 / 残作業 #51)。

壁温が不連続な run (mixA/mixB) では、底面と側壁が出会う角に $q''\propto s^{-1}$ の特異点が立つ
(§4.8.6.4)。$1/s$ の面積分は**対数発散**するので、**角近傍の壁入熱は格子収束しない**。
$G_{ij}$ は mixA/mixB の**総**壁入熱から同定しているので、この格子依存な帯を含んでいる。

そこで §4.4.2 のリップ帯と同じ扱いにする: **角から幅 $\varepsilon$ の帯を除いた入熱**で
同定し直し、$\varepsilon$ を振って係数がどれだけ動くかを見る。動かなければ係数は特異点に
依存していない、動くならその分が不確かさである。

除く帯:
  * 側壁 (`cav_outer`, `cyl_side`): $z < -\mathrm{depth}+\varepsilon$ を除く
  * 底面 (`cav_floor`): 内外の縁から $\varepsilon$ 以内 ($\min(r-R_i,\,R_o-r)<\varepsilon$) を除く

いずれも `cavity_eval.face_cut_integral` の**面切断**で除く (ノード採否だと除去範囲が格子依存)。

usage:
  python3 tools/corner_band_fit.py --runs R20 R500 R1000 --mixed RA RB [--eps 0,0.05,0.1,0.2,0.5,1.0]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import cavity_eval as ce                                   # noqa: E402
import geom_common as gc                                    # noqa: E402
from fit_wall_network import WALLS, PAIRS, JP, read_run, build   # noqa: E402


def wall_Q_nocorner(run, eps_m):
    """角帯 (幅 eps) を除いた壁ごとの入熱 [W]。eps=0 なら総入熱そのもの。"""
    run = Path(run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    step = int(ce.snapshots(run)[-1].stem.split("_")[1])
    wh = ce.wall_heat(run, step, man, D)
    G = man["geometry"]
    dep, Ri, Ro = float(G["depth"]), float(G["Ri"]), float(G["Ro"])
    out = {}
    for g in WALLS:
        d = wh.get(g)
        if d is None:
            continue
        if eps_m <= 0.0 or "_dump" not in d:
            out[g] = float(d["Q_W"])
            continue
        w = d["_dump"]
        qin = -np.asarray(w["qwall"], float)
        if g == "cav_floor":
            # 縁までの距離で切る (内外どちらの縁からも eps より遠い部分を残す)
            r = np.hypot(w["xyz"][:, 0], w["xyz"][:, 1])
            sc = np.minimum(r - Ri, Ro - r)
            out[g] = ce.face_cut_integral(w, qin, eps_m, below=False, scalar=sc)
        else:
            # 底面から eps より上を残す
            out[g] = ce.face_cut_integral(w, qin, -dep + eps_m, below=False)
    return out, wh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="一様壁温の run (3 本以上)")
    ap.add_argument("--mixed", nargs="+", required=True, help="非一様壁温の run (2 本以上)")
    ap.add_argument("--eps", default="0,0.05,0.1,0.2,0.5,1.0", help="角帯の幅 [mm] のリスト")
    a = ap.parse_args()
    eps_list = [float(x) for x in a.eps.split(",") if x != ""]
    runs = list(a.runs) + list(a.mixed)

    print("=== 角帯を除いた熱回路係数の感度 (plan §4.8.6.4 / 残作業 #51) ===")
    print("  除く帯: 側壁 = 底面から eps 以内 / 底面 = 内外の縁から eps 以内 (いずれも面切断)")
    print()
    base = {r: read_run(r) for r in runs}
    res = []
    for eps in eps_list:
        rows = []
        for r in runs:
            q, _ = wall_Q_nocorner(r, eps * 1e-3)
            d = dict(base[r])
            d["Q"] = np.array([float(q[w]) for w in WALLS])
            rows.append(d)
        A, b = build(rows)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        res.append((eps, x, rows))
        mixa = next(d for d in rows if d["name"] == Path(a.mixed[0]).name)
        print("  eps = %.2f mm  (mixA の残る壁入熱 %s W, 合計 %.3f)"
              % (eps, " / ".join("%.3f" % v for v in mixa["Q"]), float(mixa["Q"].sum())))
        print("     G_i0 : " + "  ".join("%s %.5f" % (JP[w], x[i]) for i, w in enumerate(WALLS)))
        print("     G_ij : " + "  ".join("%s-%s %.5f" % (JP[WALLS[p]], JP[WALLS[q]], x[3 + k])
                                         for k, (p, q) in enumerate(PAIRS)))
    print()
    print("  === eps=0 (総入熱) からの移動 [%] ===")
    x0 = res[0][1]
    hdr = [JP[w] for w in WALLS] + ["%s-%s" % (JP[WALLS[p]], JP[WALLS[q]]) for p, q in PAIRS]
    print("  %-12s %s" % ("eps [mm]", " ".join("%10s" % h for h in hdr)))
    for eps, x, _ in res[1:]:
        f = lambda i: (100.0 * (x[i] - x0[i]) / x0[i]) if abs(x0[i]) > 1e-12 else float("nan")
        print("  %-12.2f %s" % (eps, " ".join("%+10.1f" % f(i) for i in range(6))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
