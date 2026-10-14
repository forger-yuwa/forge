"""残差の場所 (codex diagnose 2026-10-09「延長の判断より先に res_ro の局所ノルム・寄与領域・移動/固定・振幅成長を調べる」)。
各スナップショットの VALUE/res_ro (extraFields) と volume から:
  - 欠損 Σ res_ro × 2π (符号付き) と Σ|res_ro| × 2π (絶対値)、RMS(res_ro/V)
  - 壁の x (列の壁節点の x / r_t) の区間ごとの Σ|res_ro| と、壁からの層の帯 (0、1〜5、6〜30、31〜) ごとの Σ|res_ro|
  - |res_ro|/V の上位 5 節点の (列 i、層 j = 壁から、x_w)
節点の番号 n = i·nj + j (i: 入口 → 出口、j: 軸 → 壁)。
usage (AWS の case dir): python3 cold_resloc.py <run> [<run> ...] → _band_ab/cold_pair/resloc_<run>.json
"""
import json
import math
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402

XB = [-1e9, -5.0, 0.0, 5.0, 20.0, 60.0, 1e9]
LB = [(0, 0), (1, 5), (6, 30), (31, 10_000)]


def main(run: Path):
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    ni, nj, S = XC.mesh_info(len(xy))
    xw = xy[:, 0].reshape(ni, nj)[:, -1] / S
    steps = sorted(int(m.group(1)) for p in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", p.name)))
    rows = []
    for st in steps:
        with h5py.File(run / f"res_{st}.h5", "r") as h:
            if "VALUE/res_ro" not in h:
                continue
            r = np.asarray(h["VALUE/res_ro"][:], float).reshape(ni, nj)
            v = np.asarray(h["VALUE/volume"][:], float).reshape(ni, nj) if "VALUE/volume" in h else None
        a = np.abs(r)
        rec = {"step": st, "deficit": float(r.sum() * 2 * math.pi), "abs_sum": float(a.sum() * 2 * math.pi)}
        if v is not None:
            q = a / np.maximum(v, 1e-300)
            rec["rms_res_over_V"] = float(np.sqrt(np.mean((r / np.maximum(v, 1e-300)) ** 2)))
            top = np.argsort(q, axis=None)[::-1][:5]
            rec["top5"] = [{"i": int(k // nj), "layer_from_wall": int(nj - 1 - k % nj), "x_w": round(float(xw[k // nj]), 3),
                            "res_over_V": float(q.flat[k])} for k in top]
        rec["abs_by_x"] = {f"{XB[k]:g}..{XB[k + 1]:g}": float(a[(xw >= XB[k]) & (xw < XB[k + 1]), :].sum() * 2 * math.pi) for k in range(len(XB) - 1)}
        lay = nj - 1 - np.arange(nj)
        rec["abs_by_layer"] = {f"{lo}-{hi}": float(a[:, (lay >= lo) & (lay <= hi)].sum() * 2 * math.pi) for lo, hi in LB}
        rows.append(rec)
        t = rec.get("top5", [{}])[0]
        print(f"[resloc] {run.name} {st}: 欠損 {rec['deficit']:.4f}・Σ|res| {rec['abs_sum']:.4f} kg/s、RMS(res/V) {rec.get('rms_res_over_V', float('nan')):.3e}、"
              f"最大 (列 {t.get('i')}・層 {t.get('layer_from_wall')}・x_w {t.get('x_w')})、x 別 " +
              " ".join(f"{k}:{w:.3f}" for k, w in rec["abs_by_x"].items()), flush=True)
    (XC.OUTD / f"resloc_{run.name}.json").write_text(json.dumps({"run": run.name, "rows": rows}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    for r in sys.argv[1:]:
        main(HERE / r)
