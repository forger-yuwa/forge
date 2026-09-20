#!/usr/bin/env python3
"""流入境界層の厚さ $\\delta$ を変えた `inletProfile` CSV を作る (plan §5.1 #11 の残り)。

前駆平板計算が与える入口分布は $\\delta_{99}\\simeq5.06$ mm の 1 本しか無いので、
**キャビティの入熱が流入 BL 厚にどれだけ依存するか**が測れていない。ここでは前駆の
プロファイル形状を保ったまま壁法線座標を $k$ 倍して $\\delta\\to k\\delta$ の分布を作る
(自己相似な引き伸ばし)。

  z      -> k z            (同じ z/δ で同じ値になる)
  ro,Ux,Uy,Uz,Ps,k  -> そのまま   (速度・温度スケールは変えない)
  omega  -> omega / k      (長さスケールが k 倍 -> ω ~ 1/ℓ)

**厳密な相似解ではない** (Re_θ が k 倍になるので形状係数もわずかに変わる)。
入熱の $\\delta$ 依存の**桁**を見るための感度用であることを明記して使うこと。

usage: python3 tools/scale_inlet_delta.py --in mesh/inlet_profile_m9.csv --k 0.7 \
           --out mesh/inlet_profile_m9_d07.csv
"""
import argparse
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=float, required=True, help="δ の倍率")
    a = ap.parse_args()

    lines = Path(a.src).read_text().strip().splitlines()
    hdr = lines[0].split()
    data = np.array([[float(x) for x in ln.split()] for ln in lines[1:]])
    iz = hdr.index("z")
    iw = hdr.index("omega")
    z0 = data[:, iz].copy()
    data[:, iz] = z0 * a.k
    data[:, iw] = data[:, iw] / a.k
    out = [" ".join(hdr)]
    for row in data:
        out.append(" ".join("%.9g" % v for v in row))
    Path(a.out).write_text("\n".join(out) + "\n")
    print("δ x %.3g:  z 範囲 %.4g -> %.4g mm  (%d 点)  -> %s"
          % (a.k, z0.max() * 1e3, data[:, iz].max() * 1e3, len(data), a.out))


if __name__ == "__main__":
    main()
