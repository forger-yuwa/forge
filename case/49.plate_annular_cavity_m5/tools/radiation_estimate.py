#!/usr/bin/env python3
"""壁間放射の寄与を、同定した対流コンダクタンスと同じ土俵で見積もる (plan §5.1 #23)。

§4.7.7 で「ガスを介した壁間コンダクタンス $G_{ij}$」を同定した。壁に数百 K の差が付く
運用では**放射も同じ経路で熱を運ぶ**ので、両者を同じ単位 [W/K] に揃えて比べる。

外筒内壁と円柱側面は**同軸の 2 円筒**で、すきま 3 mm に対し深さ 25 mm・周長 157 mm なので
軸方向・周方向には十分長い。長い同軸円筒の放射交換は

    q_1 = sigma (T1^4 - T2^4) / ( 1/eps1 + (A1/A2)(1/eps2 - 1) )      [W/m^2, 内側基準]

これを $G_{rad} = q_1 A_1 / (T_1-T_2)$ として [W/K] に直す (非線形なので温度差に依存する)。

**上限側の見積もり**である: 有限深さのキャビティでは開口から外へ逃げる分があり、
実際の壁間交換はこれより小さい。逆に言えば「これより小さいなら無視してよい」の判定に使える。

usage: python3 tools/radiation_estimate.py [--run RUN] [--eps 0.2,0.5,0.8]
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402

SIGMA = 5.670374419e-8


def g_rad_cyl(T1, T2, A1, A2, e1, e2):
    """同軸 2 円筒 (1 = 内側 = 円柱側面) の放射コンダクタンス [W/K]。"""
    if abs(T1 - T2) < 1e-9:
        return 0.0, 0.0
    den = 1.0 / e1 + (A1 / A2) * (1.0 / e2 - 1.0)
    q1 = SIGMA * (T1 ** 4 - T2 ** 4) / den          # [W/m2] 内側基準・正 = 内側が放出
    return q1, q1 * A1 / (T1 - T2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run_0408_s2_off000_mixB_outer1000_cyl20_ext1",
                    help="壁温を読む run (非一様壁温のもの)")
    ap.add_argument("--eps", default="0.2,0.4,0.6,0.8,0.9")
    ap.add_argument("--g-conv", type=float, default=0.04384,
                    help="比較する対流の壁間コンダクタンス [W/K] (§4.7.7 の外筒-円柱)")
    a = ap.parse_args()

    run = Path(a.run)
    man = gc.load_manifest(run=run)
    G = man["geometry"]
    D = json.loads((run / "conditions.json").read_text())
    by = D.get("wall_T_by_group", {}) or {}
    T_out = float(by.get("cav_outer", D["wall_T"]))
    T_cyl = float(by.get("cyl_side", D["wall_T"]))
    ev = json.loads((run / "cavity_eval.json").read_text())
    A_out = float(ev["wall"]["cav_outer"]["area_m2"])
    A_cyl = float(ev["wall"]["cyl_side"]["area_m2"])

    print("run %s" % run.name)
    print("  外筒内壁 %.2f K (A %.6f m2, 半割) / 円柱側面 %.2f K (A %.6f m2)"
          % (T_out, A_out, T_cyl, A_cyl))
    print("  すきま %.1f mm / 深さ %.1f mm / 周長 %.0f mm -> 同軸 2 円筒として扱える"
          % ((G["Ro"] - G["Ri"]) * 1e3, G["depth"] * 1e3, 2 * 3.14159 * G["Ri"] * 1e3))
    print("\n放射コンダクタンス G_rad [W/K] (半割) と、同定した対流 G_ij = %.5f W/K との比"
          % a.g_conv)
    print("  %-8s %14s %14s %12s" % ("eps", "q1 [kW/m2]", "G_rad [W/K]", "G_rad/G_conv"))
    for e in (float(x) for x in a.eps.split(",")):
        q1, g = g_rad_cyl(T_cyl, T_out, A_cyl, A_out, e, e)
        print("  %-8.2f %14.2f %14.5f %12.2f" % (e, q1 / 1e3, abs(g), abs(g) / a.g_conv))
    print("\n  (負号は内側=円柱側面が受け取る向き。ここでは外筒が高温なので円柱が受ける)")
    print("  **上限側の見積もり**: 有限深さなので開口から逃げる分だけ実際は小さい。")
    print("  温度差に対して非線形なので、FEM では線形化せず放射要素として与えること。")


if __name__ == "__main__":
    main()
