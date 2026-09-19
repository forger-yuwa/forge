#!/usr/bin/env python3
"""偏心スイープの比較。**h_ref の分母 (基準温度) が偏心で動く**ので、分母を固定した比較も出す。

usage:
  python3 tools/compare_offsets.py run_0200_off000 run_0201_off050 ... [--ref run_0200_off000]

h_ref = q'' / (T0_ref - T_w) は「すきま中央面の最近傍点の全温」を駆動温度に使う (ユーザ指定)。
偏心するとすきま内部の温度そのものが変わる (同心 11.9 K -> 偏心 1.0 mm で 5.3 K) ため、
**q'' がほぼ同じでも h_ref が上がる**。これは熱伝達の改善ではなく分母の縮小である。
そこで 3 通り並べる:

  q''        : 分子そのもの (モデル非依存。**まずこれを見る**)
  h_ref      : 各 run 自身の基準温度で割ったもの (ユーザ指定の定義)
  h_ref@ref  : **基準 run の基準温度**で全 run を割ったもの (分母を固定した比較)

`--ref` の run の T_ref (壁面積平均) を固定分母に使う。
"""
import argparse
import json
from pathlib import Path

CAV = ("cav_outer", "cyl_side", "cav_floor")
LAB = {"cav_outer": "外筒壁", "cyl_side": "内円柱側面", "cav_floor": "底面"}


def load(run):
    j = json.loads((Path(run) / "cavity_eval.json").read_text())
    m = json.loads((Path(run) / "manifest.json").read_text())
    return j, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--ref", default=None, help="分母を固定するときの基準 run (既定は先頭)")
    ap.add_argument("--wall-T", type=float, default=500.0)
    a = ap.parse_args()
    D = {r: load(r) for r in a.runs}
    ref = a.ref or a.runs[0]
    Tref0 = {g: D[ref][0]["wall"][g]["Tref_mean"] for g in CAV if g in D[ref][0]["wall"]}

    print("基準 run = %s  (固定分母に使う T_ref: %s)"
          % (ref, ", ".join("%s %.1f K" % (LAB[g], v) for g, v in Tref0.items())))
    for g in CAV:
        print("\n=== %s ===" % LAB[g])
        print("  %-18s %9s %11s %11s %11s %11s %11s"
              % ("run", "gap min", "q'' [W/m2]", "T_ref [K]", "h_ref", "h_ref@ref", "Q [W] 全周"))
        for r in a.runs:
            j, m = D[r]
            w = j["wall"].get(g)
            if w is None:
                continue
            Tr = w["Tref_mean"]
            h_fixed = w["qpp_mean"] / max(Tref0[g] - a.wall_T, 1e-30)
            print("  %-18s %9.2f %11.0f %11.1f %11.2f %11.2f %11.3f"
                  % (r, m["geometry"]["gap_min"] * 1e3, w["qpp_mean"], Tr,
                     w["h_ref"], h_fixed, 2 * w["Q_W"]))
        b = D[a.runs[0]][0]["wall"].get(g)
        if b:
            print("  --- 先頭 run からの変化 [%] ---")
            for r in a.runs[1:]:
                w = D[r][0]["wall"].get(g)
                if w is None:
                    continue
                hf0 = b["qpp_mean"] / max(Tref0[g] - a.wall_T, 1e-30)
                hf = w["qpp_mean"] / max(Tref0[g] - a.wall_T, 1e-30)
                print("  %-18s %9s %+10.1f %+11.1f %+11.1f %+11.1f %+11.1f"
                      % (r, "", 100 * (w["qpp_mean"] / b["qpp_mean"] - 1),
                         100 * (w["Tref_mean"] / b["Tref_mean"] - 1),
                         100 * (w["h_ref"] / b["h_ref"] - 1),
                         100 * (hf / hf0 - 1), 100 * (w["Q_W"] / b["Q_W"] - 1)))
    print("\n** q'' は分子そのもの。h_ref は各 run 自身の基準温度で割ったもの (分母が偏心で動く)。")
    print("** h_ref@ref は分母を基準 run の T_ref に固定したもの = 分子の変化だけを見る指標。")


if __name__ == "__main__":
    main()
