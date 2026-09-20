#!/usr/bin/env python3
"""#39 の単一因子 A/B — 入口分布だけを変えたときの深部の活発さ。

当初 (2026-09-20) は `run_0005_gap_t8_deep` (厚い BL) と `run_0006_xcheck_forge` (一様) を
比べて「入口を一様にすると深部が 60 倍活発」と書いたが、この 2 run は **5 因子が同時に
変わっていた** (入口分布 / 熱力学 NASA-9↔理想気体 / 輸送 Chapman-Enskog↔Sutherland /
5 成分↔単一成分 / SST 補正)。→ codex result-2 M2・result-3 で交絡ありに訂正。

ここでは `run_0005` から index コピーで継続し **`inletProfile` だけを外した**
`run_0009_gap_t8_uniform` と比べる (solverConfig は完全同一)。

指標:
  - 深さ別の壁熱流束 q_w [W/m²] (絶対値。比では深部で意味を持たない)
  - 深部 (z/W > 4) の T のばらつきを **float32 ULP 単位**で (丸めと物理を分けるため)
  - 深部の |U| の代表値
"""
import sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
DEPTHS_CM = [0.25, 0.51, 0.76, 1.52, 2.54, 3.81]
W = 0.18e-2


def latest(run):
    fs = sorted((CASE / run).glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
    if not fs:
        raise SystemExit(f"{run}: res_*.h5 が無い")
    return fs[-1]


def metrics(run):
    f = latest(run)
    with h5py.File(f) as h:
        c = h["/MESH/COORD"][:].reshape(-1, 3)
        T = h["/VALUE/T"][:].astype(np.float32)
        U = np.sqrt(h["/VALUE/Ux"][:].astype(float) ** 2 + h["/VALUE/Uy"][:].astype(float) ** 2)
    # すきま内 = y < 0 (上面が y=0)。深部は z/W > 4 すなわち深さ > 4W
    y = c[:, 1]
    deep = (y < 0) & (np.abs(y) > 4.0 * W)
    if deep.sum() == 0:
        raise SystemExit(f"{run}: 深部ノードが無い")
    Td = T[deep]
    # float32 ULP は値そのもので決まる。代表温度での ULP で割る。
    ulp = np.spacing(np.float32(float(np.mean(Td))))
    q = None
    csv = CASE / run / "gap_wall_q.csv"
    if csv.exists():
        q = np.loadtxt(csv, delimiter=",", skiprows=1)
    return dict(run=run, res=f.name, step=int(f.stem.split("_")[1]),
                n_deep=int(deep.sum()), T_mean=float(np.mean(Td)),
                T_span_K=float(np.ptp(Td)), T_span_ulp=float(np.ptp(Td) / ulp),
                U_deep_max=float(np.max(U[deep])), U_deep_rms=float(np.sqrt(np.mean(U[deep] ** 2))),
                q=q)


def main():
    a, b = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else \
           ("run_0005_gap_t8_deep", "run_0009_gap_t8_uniform")
    ma, mb = metrics(a), metrics(b)
    print(f"A = {ma['run']}  ({ma['res']}, step {ma['step']})   入口 = inletProfile 1 (厚い BL)")
    print(f"B = {mb['run']}  ({mb['res']}, step {mb['step']})   入口 = 一様")
    print("  ** 変えたのは bcondConfig の inletProfile のみ。solverConfig は完全同一。\n")
    print(f"{'指標 (深部 z/W>4, ノード数 A/B = ' + str(ma['n_deep']) + '/' + str(mb['n_deep']) + ')':52} "
          f"{'A':>14} {'B':>14} {'B/A':>8}")
    for k, lbl, fmt in (("T_mean", "T 平均 [K]", "{:14.4f}"),
                        ("T_span_K", "T のばらつき (最大-最小) [K]", "{:14.3e}"),
                        ("T_span_ulp", "同上 [float32 ULP]", "{:14.1f}"),
                        ("U_deep_max", "|U| 最大 [m/s]", "{:14.4e}"),
                        ("U_deep_rms", "|U| RMS [m/s]", "{:14.4e}")):
        va, vb = ma[k], mb[k]
        r = vb / va if va else float("nan")
        print(f"{lbl:52} " + fmt.format(va) + " " + fmt.format(vb) + f" {r:8.2f}")

    if ma["q"] is not None and mb["q"] is not None:
        print(f"\n{'深さ [cm]':>10} {'z/W':>7} {'q_w A [W/m²]':>14} {'q_w B [W/m²]':>14} {'B/A':>8}")
        for dcm in DEPTHS_CM:
            d = dcm * 1e-2
            qa = float(np.interp(d, ma["q"][:, 0], ma["q"][:, 1]))
            qb = float(np.interp(d, mb["q"][:, 0], mb["q"][:, 1]))
            r = qb / qa if abs(qa) > 0 else float("nan")
            print(f"{dcm:10.2f} {d/W:7.2f} {qa:14.4e} {qb:14.4e} {r:8.2f}")
    else:
        print("\n(gap_wall_q.csv が両方に無い。先に tools/gap_eval.py を各 run で回すこと)")
    print("\n注: B が未完了なら step を確認すること。深部は最も遅く落ち着くので、"
          "両者が同じ step で、かつ量が STEADY になってから比べる。")


if __name__ == "__main__":
    main()
