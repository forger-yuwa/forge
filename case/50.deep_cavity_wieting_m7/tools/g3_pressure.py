#!/usr/bin/env python3
"""G3 — キャビティ内圧 / 模型面圧 の照合 ([W70] Fig 4: p_c/p_m = 1.0 ± 0.1)。

分母問題 (q_fp の再構成) の影響を受けない検証点なので、すきま熱伝達とは独立に
「開口キャビティの基本挙動が出ているか」を確かめられる。
[W70] の圧力孔は前縁から 3.5 in = 88.9 mm の**模型面上**、キャビティ側は床近傍。
"""
import json
from pathlib import Path
import h5py
import numpy as np

CASE = Path(__file__).resolve().parents[1]
RUNS = {0.063: "run_0016_T1_wd0063_fine_settle", 0.211: "run_0008_T1_wd0211_long",
        0.383: "run_0009_T1_wd0383_long", 0.524: "run_0011_T1_wd0524_long"}
X_PM = 88.9e-3   # 模型圧力孔 (前縁から 3.5 in)


def main():
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    d = geom["cavity"]["depth"] * 1e-3
    xr = geom["cavity"]["x_rear_wall_from_le"] * 1e-3
    print(f"{'w/d':>6} {'p_c [Pa]':>10} {'p_m [Pa]':>10} {'p_c/p_m':>9} {'判定':>6}   run")
    out = {}
    for wd, run in RUNS.items():
        w = geom["cavity"]["widths"][str(wd)] * 1e-3
        f = sorted((CASE / run).glob("res_[0-9]*.h5"),
                   key=lambda p: int(p.stem.split("_")[1]))[-1]
        with h5py.File(f) as h:
            c = h["/MESH/COORD"][:].reshape(-1, 3)
            P = h["/VALUE/P"][:].astype(float)
        # キャビティ内: 床から 1.27 mm、幅の中央に最も近いノード
        yt, xt = -d + 1.27e-3, xr - 0.5 * w
        m = (c[:, 1] < -0.5 * d) & (np.abs(c[:, 0] - xt) < w)
        k = np.argmin((c[m, 0] - xt) ** 2 + (c[m, 1] - yt) ** 2)
        p_c = float(P[m][k])
        # 模型面 (y=0 の板上、x=88.9 mm) 近傍の平均
        m2 = (np.abs(c[:, 1]) < 1e-5) & (np.abs(c[:, 0] - X_PM) < 2e-3)
        p_m = float(np.mean(P[m2])) if m2.sum() else float("nan")
        r = p_c / p_m
        ok = "PASS" if abs(r - 1.0) <= 0.1 else "FAIL"
        print(f"{wd:6.3f} {p_c:10.1f} {p_m:10.1f} {r:9.3f} {ok:>6}   {run}")
        out[str(wd)] = dict(run=run, p_c=p_c, p_m=p_m, ratio=r, verdict=ok, n_model_nodes=int(m2.sum()))
    v = "PASS" if all(o["verdict"] == "PASS" for o in out.values()) else "FAIL"
    print(f"\nG3 VERDICT: {v}   (実験 [W70] Fig 4: p_c/p_m = 1.0 ± 0.1, Re・スパン長に依らず)")
    (CASE / "g3_pressure.json").write_text(
        json.dumps(dict(rows=out, verdict=v, ref="W70 Fig 4: 1.0 +- 0.1"), indent=2, ensure_ascii=False),
        encoding="utf-8")
    print(f"  → {CASE/'g3_pressure.json'}")


if __name__ == "__main__":
    main()
