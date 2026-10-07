#!/usr/bin/env python3
"""dual-time の「分割 (ckpt100 → 再開 100)」と「連続 200」の比較 (plan §6 の構成表: 分割と連続の差)。AWS で回す。

    python3 split_vs_cont.py [--split-cfg c09restart100] [--cont-cfg c09cont200] [--builds base new]

registry.tsv から各ビルドの再開 run (res_100.h5 = 物理時刻 200 step 目) と連続 run (res_200.h5) を選び、
主要データセットの m = max|A−B|/max|A| を「分割 vs 連続」(3×3 組の最大) と、同じ側の反復どうし (3 組の最大) で並べる。
判定はしない (分割と連続は再開時に保存量以外の状態が移らないので一致しない。差の大きさをビルド間で比べるための表)。
"""
import argparse
import itertools
import os

import compare_runs as c

NAMES = ["VALUE/ro", "VALUE/roUx", "VALUE/roUy", "VALUE/roUz", "VALUE/roe", "VALUE/roXi", "VALUE/roY0", "VALUE/roY1",
         "VALUE/P", "VALUE/T", "CHECKPOINT/roN", "CHECKPOINT/roUxN", "CHECKPOINT/roeN", "CHECKPOINT/roXiP"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-cfg", default="c09restart100")
    ap.add_argument("--cont-cfg", default="c09cont200")
    ap.add_argument("--split-file", default="res_100.h5")
    ap.add_argument("--cont-file", default="res_200.h5")
    ap.add_argument("--builds", nargs="+", default=["base", "new"])
    a = ap.parse_args()
    data, names = {}, None
    for b in a.builds:
        for kind, cfg, fn in (("分割", a.split_cfg, a.split_file), ("連続", a.cont_cfg, a.cont_file)):
            runs = c.runs_for(cfg, b)
            if names is None:
                ds, _ = c.h5_items(os.path.join(runs[0], fn))
                names = [n for n in NAMES if n in ds]
            data[(b, kind)] = [c.h5_data(os.path.join(r, fn), names) for r in runs]
            print(f"{b} {kind}: {', '.join(os.path.basename(r) for r in runs)} ({fn})")

    def cross(x, y, n):
        return max(c.metric(p[n], q[n]) for p in x for q in y)

    def within(x, n):
        return max((c.metric(p[n], q[n]) for p, q in itertools.combinations(x, 2)), default=float("nan"))

    hdr = f"{'量':20s}"
    for b in a.builds:
        hdr += f" | {b}: 分割vs連続  連続内     分割内"
    print("\nm = max|A−B|/max|A| (分割vs連続は 3×3 組の最大、内は同じ側の反復 3 組の最大)")
    print(hdr)
    for n in names:
        line = f"{n:20s}"
        for b in a.builds:
            line += (f" | {cross(data[(b, '分割')], data[(b, '連続')], n):11.3e}"
                     f" {within(data[(b, '連続')], n):9.3e} {within(data[(b, '分割')], n):9.3e}")
        print(line)


if __name__ == "__main__":
    main()
