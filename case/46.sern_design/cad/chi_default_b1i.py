#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B1 (i) 接続模型の救済: chi 省略 (auto → 1) の run を前 plan の V1 基準で判定する。

基準 (前 plan convection-slau-wall-normal-chi §6、測る前に固定済み):
  V1-e/V1-a: 監視 3 CV (153797・153880・189814) が step 200 までに 10 ρMin (= 1e-3、roMin 既定 1e-4) を超え、以後終了まで維持
  V1-d: 全壁ノードで床到達 0 (保存スナップショット)。NaN 0。初回面流束が明示 1 とビット同一 (別途比較済みの結果を引数で渡す)
プローブは T・P しか出さないので、ρ は **下限 P/(R_max T)** (R_max = 排気の 340.4、空気は 287.1) で評価する (安全側)。
スナップショットでは VALUE/ro を直接読む。

  python3 chi_default_b1i.py RUN_DIR [--flux-mismatch N] [--out FILE]
"""
import argparse
import glob
import os
import re

import h5py
import numpy as np

MON = {"p0_153797": 153797, "p1_153880": 153880, "p2_189814": 189814}
R_MAX, R_MIN = 340.4, 287.1
THR = 1.0e-3


def probe(run, i):
    rows = []
    for line in open(os.path.join(run, f"point_probe_{i}.out")):
        p = [x.strip() for x in line.split(",")]
        if not p[0].isdigit():
            continue
        rows.append((int(p[0]), float(p[2]), float(p[3])))
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--flux-mismatch", type=int, default=None); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    names = [l.split(":")[0].strip() for l in open(os.path.join(a.run, "probe.yaml")) if re.match(r"\s+p\d_", l)]
    L = [f"# chi 既定化 B1 (i) 接続模型の救済: {a.run}", ""]
    ok = True
    L.append("| 監視点 | 10ρMin 到達 step (下限 ρ) | 到達後の下限 ρ 最小 | 最終 step の ρ 範囲 [下限, 上限] | V1-e (≤ 200 で到達し維持) |")
    L.append("| --- | --- | --- | --- | --- |")
    for i, nm in enumerate(names):
        if nm not in MON:
            continue
        d = probe(a.run, i)
        st, T, P = d[:, 0], d[:, 1], d[:, 2]
        lo, hi = P / (R_MAX * T), P / (R_MIN * T)
        k = np.where(lo > THR)[0]
        first = int(st[k[0]]) if k.size else None
        after = lo[k[0]:].min() if k.size else float("nan")
        good = first is not None and first <= 200 and after > THR
        ok &= good
        L.append(f"| {nm} | {first} | {after:.3e} | [{lo[-1]:.3e}, {hi[-1]:.3e}] | {'PASS' if good else 'FAIL'} |")
    L.append("")
    L.append("| スナップショット | 監視 3 節点の ρ | 壁ノード ρ 最小 | 床到達 (ρ ≤ roMin 1e-4 の壁ノード) | NaN |")
    L.append("| --- | --- | --- | --- | --- |")
    for f in sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda s: int(re.search(r"res_(\d+)", s).group(1))):
        with h5py.File(f, "r") as h:
            ro = np.array(h["VALUE/ro"]); wd = np.array(h["VALUE/wall_dist"])
        w = wd <= 0
        nfl = int(np.count_nonzero(ro[w] <= 1.0e-4)); nn = int(np.isnan(ro).sum())
        if not f.endswith("res_0.h5"):
            ok &= nfl == 0 and nn == 0
        L.append(f"| {os.path.basename(f)} | {', '.join(f'{ro[n]:.3e}' for n in MON.values())} | {ro[w].min():.3e} | {nfl} | {nn} |")
    if a.flux_mismatch is not None:
        ok &= a.flux_mismatch == 0
        L.append(f"\n初回面流束 (省略 vs 明示 1) の不一致面数: {a.flux_mismatch}")
    L.append(f"\nVERDICT B1 (i): {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(L) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
