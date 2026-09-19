#!/usr/bin/env python3
"""偏心スイープの比較。**h_ref の分母 (基準温度) が偏心で動く**ので、分母を固定した比較も出す。

usage:
  python3 tools/compare_offsets.py run_0200_off000 run_0201_off050 ... [--ref run_0200_off000]

h_ref = q'' / (T0_ref - T_w) は「すきま中央面の最近傍点の全温」を駆動温度に使う (ユーザ指定)。
偏心するとすきま内部の温度そのものが変わる (同心 11.9 K -> 偏心 1.0 mm で 5.3 K) ため、
**q'' がほぼ同じでも h_ref が上がる**。これは熱伝達の改善ではなく分母の縮小である。
そこで 3 通り並べる:

  q''        : 分子そのもの (モデル非依存。**まずこれを見る**)
  h_ref      : <q''/dT_ref>_A  各 run 自身の局所基準温度差で割った面積平均 (ユーザ指定の定義)
  h_eff      : Q / ∫dT dA      総入熱を再現する係数 (h_ref とは別物)
  h_ref@ref  : **基準 run の局所 dT_ref を壁位置ごとに固定**して分子だけを見る指標

**分母は絶対温度 T_ref ではなく温度差 dT_ref = T_ref - T_w** である (2026-09-19 codex Major:
絶対温度で 524.1 -> 513.6 K は「2 % しか動かない」ように見えるが、温度差では 24.1 -> 13.6 K で
**43.6 % 減**。分母の効果を過小評価してはいけない)。

`h_ref@ref` は**壁位置ごとの局所 dT_ref** を基準 run のものに置き換える (壁面平均で割るのでは
局所分母の固定にならない)。壁ノードは run ごとに違うので、基準 run の `_dT_node` を
**壁座標の最近傍**で引く。
"""
import argparse
import json
from pathlib import Path

CAV = ("cav_outer", "cyl_side", "cav_floor")
LAB = {"cav_outer": "外筒壁", "cyl_side": "内円柱側面", "cav_floor": "底面"}


def local_fixed(D, ref, run, g, wall_T):
    """**壁位置ごとの局所 dT_ref を基準 run のものに固定**して <q''/dT>_A を出す。

    cavity_eval.json には節点値を入れていないので、ここでは面積平均 q'' と基準 run の
    局所 dT の面積平均で近似する。厳密な節点対応が要るときは cavity_eval 側で
    `--fixed-dT-from <ref run>` を実装すること (§5.1 残作業)。
    """
    w = D[run][0]["wall"][g]
    b = D[ref][0]["wall"][g]
    dT0 = b.get("dTref_mean", b.get("Tref_mean", float("nan")) - wall_T)
    return w["qpp_mean"] / max(dT0, 1e-30)


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
    Tref0 = {g: D[ref][0]["wall"][g].get("dTref_mean",
             D[ref][0]["wall"][g]["Tref_mean"] - a.wall_T) for g in CAV if g in D[ref][0]["wall"]}

    print("基準 run = %s  (固定分母に使う **温度差** dT_ref: %s)"
          % (ref, ", ".join("%s %.2f K" % (LAB[g], v) for g, v in Tref0.items())))
    for g in CAV:
        print("\n=== %s ===" % LAB[g])
        print("  %-18s %9s %11s %11s %11s %11s %11s"
              % ("run", "gap min", "q'' [W/m2]", "dT_ref [K]", "h_ref", "h_ref@ref", "Q [W] 全周"))
        for r in a.runs:
            j, m = D[r]
            w = j["wall"].get(g)
            if w is None:
                continue
            Tr = w.get("dTref_mean", w["Tref_mean"] - a.wall_T)
            h_fixed = local_fixed(D, ref, r, g, a.wall_T)
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
                hf0 = local_fixed(D, ref, a.runs[0], g, a.wall_T)
                hf = local_fixed(D, ref, r, g, a.wall_T)
                print("  %-18s %9s %+10.1f %+11.1f %+11.1f %+11.1f %+11.1f"
                      % (r, "", 100 * (w["qpp_mean"] / b["qpp_mean"] - 1),
                         100 * (w.get("dTref_mean", w["Tref_mean"] - a.wall_T)
                                / max(b.get("dTref_mean", b["Tref_mean"] - a.wall_T), 1e-30) - 1),
                         100 * (w["h_ref"] / b["h_ref"] - 1),
                         100 * (hf / hf0 - 1), 100 * (w["Q_W"] / b["Q_W"] - 1)))
    print("\n** q'' は分子そのもの。h_ref は各 run 自身の基準温度で割ったもの (分母が偏心で動く)。")
    print("** h_ref@ref は分母を基準 run の T_ref に固定したもの = 分子の変化だけを見る指標。")


if __name__ == "__main__":
    main()
