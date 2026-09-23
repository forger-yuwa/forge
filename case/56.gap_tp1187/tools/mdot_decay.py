#!/usr/bin/env python3
"""すきま断面を通る正味の質量流束 |∫ρU_y dx| の step 依存を出す。

**なぜこの量か**: 閉じた袋小路の定常解ではどの深さでもゼロ。ゼロからの隔たりが
「まだ収束していない量」そのものであり、壁熱流束のような派生量と違って
離散化の細部に依らない。減り方が

  - **べき乗則** (log-log 直線) なら step を積んでも届かない
  - **幾何級数** (semi-log 直線) なら有限 step で収束する

なので、収束の「質」をこの傾きで判定できる (§4.7-6c / 4.7-6d-2)。

usage: mdot_decay.py RUN_DIR [RUN_DIR ...] [--zw 15] [--offset N]
"""
import argparse, glob, os
import numpy as np
import h5py

W = 0.18e-2
XW = 0.9e-3


def mdot(res, zW, signed=False):
    """帯 z/W の質量流束 |∫ ρU_y dx|。

    **保存量 `roUy` を直接積分する** (2026-09-24 に修正、codex plan-3 m3)。
    以前は `ro * Uy` = 「commit 後の保存量 × commit **前**の原始量」を掛けていた。
    `res_*.h5` の原始量は commit 前の評価値なので、**更新時点の違う量を混ぜていた**。
    床を論じる指標に不要な誤差を持ち込むので直した (実測の差は run_0026/res_400000 で
    5.90398e-10 → 5.88065e-10、0.4 %。減衰率の結論は変わらない)。

    `signed=True` で符号つきを返す。**床付近は符号が零をまたぐ**ので、
    絶対値だけ見ていると「床に張り付いている」ように見える (実際は零の周りの振れ)。
    """
    with h5py.File(res) as h:
        c = h["/MESH/COORD"][:].reshape(-1, 3)
        if "/VALUE/roUy" in h:
            f = h["/VALUE/roUy"][:].astype(float)          # 保存量をそのまま
        else:                                              # 古い出力に roUy が無い場合の後方互換
            f = h["/VALUE/ro"][:].astype(float) * h["/VALUE/Uy"][:].astype(float)
    m = (np.abs(c[:, 1] + zW * W) < 0.12e-3) & (np.abs(c[:, 0]) <= XW + 1e-9)
    if m.sum() < 5:
        return np.nan
    o = np.argsort(c[m, 0])
    v = float(np.trapz(f[m][o], c[m, 0][o]))
    return v if signed else abs(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--zw", type=float, default=15.0)
    ap.add_argument("--offset", type=int, default=0, help="通算 step の下駄")
    a = ap.parse_args()
    for rd in a.runs:
        fs = sorted(glob.glob(os.path.join(rd, "res_[0-9]*.h5")),
                    key=lambda p: int(os.path.basename(p)[4:-3]))
        pts = [(int(os.path.basename(f)[4:-3]) + a.offset, mdot(f, a.zw)) for f in fs]
        pts = [(s, v) for s, v in pts if s > 0 and np.isfinite(v)]
        if len(pts) < 3:
            print(f"{rd}: 点が {len(pts)} 個しかない"); continue
        s = np.array([p[0] for p in pts], float); v = np.array([p[1] for p in pts])
        print(f"\n=== {rd}  (z/W={a.zw:g}, {len(pts)} 点, step {s[0]:.0f}..{s[-1]:.0f}) ===")
        for i in range(0, len(pts), max(1, len(pts) // 6)):
            print(f"   step {s[i]:9.0f}   |mdot| = {v[i]:.4e}")
        print(f"   step {s[-1]:9.0f}   |mdot| = {v[-1]:.4e}")
        # べき乗則 (log-log) と 幾何級数 (semi-log) のどちらが説明するか。
        # **相関で比べない**: 相関は説明変数のスケールに依るので、step 範囲が狭いと
        # log(step) がほぼ定数になり log-log が不当に勝つ (2026-09-23 に実際に誤判定した:
        # step 1.60M→1.65M の 3 点で「べき乗則 step^-67.6、相関 -1.0000」と出た)。
        # 同じ log(|mdot|) 空間での残差 RMS で比べ、step 幅が足りなければ判定を保留する。
        # **減衰区間と末尾区間を分ける** (codex plan M3)。床付近まで一括 fit すると
        # 「べき乗則」「DRIFTING」と出るが、これは判定が不適切なのであって指数減衰が
        # 無いのではない。末尾の頭打ちを検出し、そこから前を減衰区間として扱う。
        n = len(v)
        knee = n
        if n >= 6:
            lv = np.log(v)
            # 末尾から遡り、単調減少が崩れた (= 床に入った) 最初の点を探す
            for k in range(n - 1, 2, -1):
                if lv[k] >= lv[k - 1]:
                    knee = k
                else:
                    break
        if knee < n:
            tail = v[knee - 1:]
            print(f"   **末尾 {n-knee+1} 点は頭打ち**: {v[knee-1]:.4e} .. {v[-1]:.4e}"
                  f"  (平均 {tail.mean():.4e} ± {tail.std():.1e}, 振れ幅 {tail.max()/tail.min():.2f} 倍)")
            print(f"   -> **床としてはこの平均±振幅を報告する。最小値を床と呼ばないこと。**")
            if knee >= 4:
                bd, ad = np.polyfit(s[:knee], np.log(v[:knee]), 1)
                print(f"   減衰区間 (先頭 {knee} 点): |mdot| ∝ exp({bd:+.3e}·step)"
                      f"  e 折り {abs(1/bd):.3e} step、{v[0]/v[knee-1]:.3g} 倍落ちた")
            s, v = s[:knee], v[:knee]   # 以降の判定は減衰区間だけで行う
        span = np.log10(s[-1] / s[0]) if len(s) > 1 else 0.0
        bp, ap_ = np.polyfit(np.log(s), np.log(v), 1)
        be, ae = np.polyfit(s, np.log(v), 1)
        rp = float(np.sqrt(np.mean((np.log(v) - (bp * np.log(s) + ap_)) ** 2)))
        re = float(np.sqrt(np.mean((np.log(v) - (be * s + ae)) ** 2)))
        print(f"   べき乗則  |mdot| ∝ step^{bp:+.4f}       残差 RMS {rp:.4f}"
              f"   (step 倍加あたり {(1-2**bp)*100:+.1f} %)")
        print(f"   幾何級数  |mdot| ∝ exp({be:+.3e}·step)  残差 RMS {re:.4f}"
              f"   (1e5 step あたり {(1-np.exp(be*1e5))*100:+.1f} %、e 折り {abs(1/be):.3e} step)")
        if span < 0.3:
            print(f"   -> **判定保留**: step 範囲が {span:.2f} decade しかない "
                  f"(log-log の判別には 0.3 decade 以上要る)。幾何級数側の数値だけ使うこと。")
        else:
            print(f"   -> {'**べき乗則** (step を積んでも届かない)' if rp < re else '**幾何級数** (有限 step で収束する)'}"
                  f"   [step 範囲 {span:.2f} decade]")


if __name__ == "__main__":
    main()
