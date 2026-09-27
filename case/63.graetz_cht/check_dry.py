#!/usr/bin/env python3
r"""乾式 1 step の壁ダンプで境界の帰属を検査する (壁 3/4/5 は r = R、x 範囲が区間どおり)。

2026-09-27: gen_mesh.py の初版は Extrude の戻り値の対応を取り違え、壁と軸が入れ替わっていた (壁ダンプが y = 0)。

    python3 case/63.graetz_cht/check_dry.py case/63.graetz_cht/run_0002_dry_r16
"""
import sys
from pathlib import Path

import h5py
import numpy as np

import graetz_common as gc

run = Path(sys.argv[1])
want = {"wall_up_3": (-gc.L_UP, 0.0), "wall_heat_4": (0.0, gc.L_HEAT), "wall_down_5": (gc.L_HEAT, gc.L_HEAT + gc.L_DOWN)}
ok = True
for name, (x0, x1) in want.items():
    f = sorted(run.glob(f"res_{name}_*.h5"))
    if not f:
        print(f"FAIL {name}: 壁ダンプが無い"); ok = False; continue
    with h5py.File(f[-1], "r") as h:
        c = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
        okf = h["VALUE/iface_ok"][:]
    good = (np.allclose(c[:, 1], gc.R, rtol=0, atol=1e-9) and abs(c[:, 0].min() - x0) < 1e-8
            and abs(c[:, 0].max() - x1) < 1e-8 and bool(np.all(okf == 1)))
    print(f"{'PASS' if good else 'FAIL'} {name}: {len(c)} 節点、y {c[:,1].min():.9g}..{c[:,1].max():.9g}、"
          f"x {c[:,0].min():.6g}..{c[:,0].max():.6g}、iface_ok 全 1: {bool(np.all(okf == 1))}")
    ok &= good
print("VERDICT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
