#!/usr/bin/env python3
"""FORGE_MEMLOG=1 の工程別メモリ (forge_run.log の [memlog] 行) と run_matrix.py --memwatch の採取 (mem_samples.csv) を表にする。

    python3 memlog_summary.py RUN [RUN ...] [--nodes N]

- 工程ごとに VmRSS / VmHWM [MB] と GPU (cudaMemGetInfo の total−free = **他プロセス込み**、と最初の memlog からの増分) を run ごとに並べる。
- --nodes (省略時はログの "Number of Cells") で割った B/節点 も出す (1 格子なので傾きではなく「切片込みの平均」)。
- mem_samples.csv からは自分の forge の VmHWM 最大・GPU (プロセス別、nvidia-smi) 最大を出す。
"""
import argparse
import csv
import os
import re


def parse(run):
    txt = open(os.path.join(run, "forge_run.log"), errors="replace").read()
    m = re.search(r"Number of Cells: (\d+)", txt)
    n = int(m.group(1)) if m else None
    rows, seen = [], {}
    for m in re.finditer(r"\[memlog\] (.+?) @(\S+) VmRSS=(\d+)MB VmHWM=(\d+)MB(?: \| (.*))?", txt):
        lab = m.group(1)
        seen[lab] = seen.get(lab, 0) + 1
        key = lab if seen[lab] == 1 else f"{lab}#{seen[lab]}"
        ex = m.group(5) or ""
        g = re.search(r"GPU used=(\d+)MB \(since first memlog \+(-?\d+)MB\)", ex)
        items = dict((a, float(b)) for a, b in re.findall(r"(\S+?)\[\d+\]=([\d.e+-]+)MB", ex))
        rows.append(dict(key=key, rss=int(m.group(3)), hwm=int(m.group(4)),
                         gpu=int(g.group(1)) if g else None, gpud=int(g.group(2)) if g else None, items=items))
    samp = []
    p = os.path.join(run, "mem_samples.csv")
    if os.path.exists(p):
        with open(p) as f:
            samp = list(csv.DictReader(f))
    return n, rows, samp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--nodes", type=int)
    ap.add_argument("--items", action="store_true", help="最後の工程のコンテナ別内訳も出す")
    a = ap.parse_args()
    data = {r: parse(r) for r in a.runs}
    names = [os.path.basename(r.rstrip("/")) for r in a.runs]
    n = a.nodes or next((v[0] for v in data.values() if v[0]), None)
    print(f"節点数 {n}")
    keys = [x["key"] for x in data[a.runs[0]][1]]
    hdr = " | ".join(f"{nm[-12:]:>28s}" for nm in names)
    print(f"{'工程':52s} | {hdr}")
    print(f"{'':52s} | " + " | ".join(f"{'RSS/HWM MB  GPU(Δ) MB':>28s}" for _ in names))
    for k in keys:
        cells = []
        for r in a.runs:
            x = next((y for y in data[r][1] if y["key"] == k), None)
            if x is None:
                cells.append(f"{'-':>28s}")
            else:
                g = f"{x['gpu']}({x['gpud']:+d})" if x["gpu"] is not None else "-"
                cells.append(f"{x['rss']:6d}/{x['hwm']:6d} {g:>14s}")
        print(f"{k[:52]:52s} | " + " | ".join(cells))
    print()
    for r, nm in zip(a.runs, names):
        nn, rows, samp = data[r]
        hwm = max((x["hwm"] for x in rows), default=0)
        last = rows[-1] if rows else None
        s_hwm = max((int(s["VmHWM_kB"]) for s in samp), default=-1)
        s_rss = max((int(s["VmRSS_kB"]) for s in samp), default=-1)
        s_gpu = max((int(s["gpu_proc_MiB"]) for s in samp), default=-1)
        s_gt = max((int(s["gpu_total_MiB"]) for s in samp), default=-1)
        bpn = f"{hwm * 1048576 / nn:.0f} B/節点" if nn else "-"
        gd = max((x["gpud"] for x in rows if x["gpud"] is not None), default=None)
        print(f"{nm}: memlog VmHWM 最大 {hwm} MB ({bpn})、GPU 増分最大 {gd} MB; "
              f"採取 (1 s): VmHWM {s_hwm / 1024:.0f} MB・VmRSS 最大 {s_rss / 1024:.0f} MB・GPU (プロセス) {s_gpu} MiB・GPU 合計 {s_gt} MiB "
              f"({len(samp)} 点)")
        if a.items and last:
            for k2, v in last["items"].items():
                print(f"    {k2:36s} {v:9.1f} MB")


if __name__ == "__main__":
    main()
