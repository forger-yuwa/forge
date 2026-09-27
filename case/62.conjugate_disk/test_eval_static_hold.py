#!/usr/bin/env python3
"""eval_static_hold.py の熱流束判定の負例 (plan axisymmetric-graded-grid-static-gas、2026-09-27 plan レビュー M1)。

合成 run ディレクトリ (mesh.h5・res_<step>.h5・両壁ダンプ) を一時ディレクトリに作り、**評価器本体を実行**する。
正しい向き (加熱壁 −120.5 / 共役壁 +120.5 W/m²) は PASS、両壁の符号を入れ替えると FAIL になること
(旧実装は絶対値比較なので入れ替えても PASS した)。判定は VERDICT 行と終了コードの組で照合する。
"""
import os, subprocess, sys, tempfile
import h5py, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); EVAL = os.path.join(HERE, "eval_static_hold.py")
Q = 0.0241 * 25 / 0.005

def make(d, qhot, qcj):
    n = 20
    with h5py.File(f"{d}/mesh.h5", "w") as m:
        m["CELLS/volume"] = np.ones(n); m["CELLS/centCoords"] = np.column_stack([np.zeros(n), np.linspace(.005, .02, n), np.zeros(n)]).ravel()
    for st in (10000, 15000, 20000):
        with h5py.File(f"{d}/res_{st}.h5", "w") as h:
            for k in ("Ux", "Uy"): h[f"VALUE/{k}"] = np.zeros(n)
            h["VALUE/P"] = np.full(n, 1000.0)
        for stem, q in (("res_wall_hot_3", qhot), ("res_wall_cj_4", qcj)):
            with h5py.File(f"{d}/{stem}_{st}.h5", "w") as w:
                r = np.linspace(.005, .02, 9); w["MESH/COORD"] = np.column_stack([np.zeros(9), r, np.zeros(9)]).ravel()
                w["VALUE/iface_q_eff"] = np.full(9, q)

def run(qhot, qcj):
    with tempfile.TemporaryDirectory() as d:
        make(d, qhot, qcj)
        p = subprocess.run([sys.executable, EVAL, d, "--win0", "10000", "--win1", "20000"], capture_output=True, text=True)
        v = [l for l in p.stdout.splitlines() if l.startswith("VERDICT")]
        return (v[-1].split()[2] if v else "(none)"), p.returncode

cases = [("正しい向き (hot −120.5 / cj +120.5)", (-Q, Q), ("PASS", 0)),
         ("両壁の符号を入れ替え", (Q, -Q), ("FAIL", 1)),
         ("共役壁だけ符号反転", (-Q, -Q), ("FAIL", 1))]
bad = 0
for nm, qq, exp in cases:
    got = run(*qq); good = got == exp; bad += not good
    print(f"[{'OK' if good else 'NG'}] {nm}: 期待 {exp} -> {got}")
print(f"\nVERDICT: {'PASS' if not bad else 'FAIL'} ({len(cases)-bad}/{len(cases)})")
raise SystemExit(bad)
