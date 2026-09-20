#!/usr/bin/env python3
"""リップ帯を均す幅をどこまで広げてよいか (plan §5.1 #25, codex 2026-09-19)。

開口リップは幾何的特異点で $q''\\propto s^{-1/2}$ と発散するため、ピーク $q''$ は格子収束しない
(§4.4.2)。FEM には「リップ帯の積分入熱を保存したまま幅 $\\varepsilon$ に均した値」を渡すが、
**熱量を保存しても局所温度勾配は保存されない**。均し幅が構造側の熱拡散長より大きいと
勾配を潰してしまう (= 熱応力に対して非安全側になり得る)。

判定は単純で、**均し幅 $\\varepsilon$ を評価したい時刻の熱拡散長 $\\sqrt{\\alpha t}$ と比べる**。
$\\varepsilon \\lesssim \\sqrt{\\alpha t}$ なら、構造自身がその幅を均してしまうので害がない。

`cavity_eval.py --lip-scan` が出す Q(ε) を入力に、各 ε での均し熱流束を出す。

usage: python3 tools/lip_smear.py --qtot 211.8 --scan 0.1:193.7,0.2:181.3,... [--alpha 3.09e-6]
"""
import argparse
import math


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qtot", type=float, required=True, help="全周の総入熱 [W]")
    ap.add_argument("--scan", required=True,
                    help="`eps_mm:Q_W` をカンマ区切り (cavity_eval --lip-scan の出力)")
    ap.add_argument("--perim-mm", type=float, default=295.3,
                    help="リップ帯が乗る周長の合計 [mm] (外筒 2πRo + 円柱 2πRi)")
    ap.add_argument("--alpha", type=float, default=3.09e-6,
                    help="構造の熱拡散率 [m2/s] (既定 Inconel 718 相当)")
    ap.add_argument("--times", default="0.01,0.1,1,10",
                    help="評価したい時刻 [s]")
    a = ap.parse_args()

    pts = []
    for kv in a.scan.split(","):
        e, _, q = kv.partition(":")
        pts.append((float(e) * 1e-3, float(q)))
    pts.sort()

    print("リップ帯に集中する入熱と、幅 ε に均したときの熱流束")
    print("  (帯の入熱 = 総入熱 %.1f W − ε を除いた入熱。周長 %.1f mm)"
          % (a.qtot, a.perim_mm))
    print("  %-10s %12s %10s %14s" % ("ε [mm]", "帯の入熱 [W]", "帯の割合", "均した q'' [kW/m2]"))
    for e, q in pts:
        Qb = a.qtot - q
        A = a.perim_mm * 1e-3 * e
        print("  %-10.2f %12.2f %9.1f %% %14.1f"
              % (e * 1e3, Qb, 100 * Qb / a.qtot, Qb / A / 1e3))

    print("\n構造の熱拡散長 √(αt)   (α = %.3g m2/s)" % a.alpha)
    print("  %-12s %14s %s" % ("t [s]", "√(αt) [mm]", "この時刻で許される均し幅"))
    for t in (float(x) for x in a.times.split(",")):
        L = math.sqrt(a.alpha * t) * 1e3
        ok = [e * 1e3 for e, _ in pts if e * 1e3 <= L]
        print("  %-12g %14.3f %s" % (t, L, ("ε ≤ %.2f mm" % max(ok)) if ok else "**どの ε も広すぎる**"))
    print("\n  ε ≲ √(αt) なら構造自身がその幅を均すので、熱量さえ合っていれば害はない。")
    print("  逆に ε > √(αt) の時刻を評価するなら、均さずリップ形状 (フィレット半径) を入れること。")


if __name__ == "__main__":
    main()
