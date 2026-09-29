#!/usr/bin/env python3
r"""case/65 の板 (0 ≤ x ≤ L、−b ≤ y ≤ 0)。下面 y = −b に Robin (h 1e8、T_h = T_∞ + 10 K)、両端面は断熱。
界面節点は乾式 1 step の壁ダンプ `res_plate_5_1.h5` の座標そのもの。層数は --ns (既定 = 流体の N)。

    python3 gen_solid.py --wall run_0002_dry_n32/res_plate_5_1.h5 --case C1
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import plate_common as pc
from plate_common import axcht, gc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True); ap.add_argument("--case", choices=sorted(pc.KS_RATIO), required=True); ap.add_argument("--ns", type=int)
    a = ap.parse_args()
    c = axcht.read_wall_coords(Path(a.wall))
    if np.abs(c[:, 1]).max() > 1e-12 or abs(c[:, 0].min()) > 1e-9 or abs(c[:, 0].max() - pc.L) > 1e-7:
        raise SystemExit("壁ダンプが y=0、0 ≤ x ≤ L の板でない")
    iface = c[np.argsort(c[:, 0]), :2]
    tag = Path(a.wall).resolve().parent.name
    ns = a.ns or int(tag.split("_n")[-1])
    nodes, tris, outer, back = axcht.solid_strip(iface, normal_axis=1, back_coord=-pc.B, ns=ns)
    k_s = pc.KS_RATIO[a.case] * gc.K_F; Th = gc.T_IN + pc.DT_H
    out = pc.HERE / "mesh" / f"solid_{tag}_{a.case}"
    meta = {"case": f"case/65 {a.case} 板", "wall_dump": str(Path(a.wall)), "wall_dump_sha256": axcht.sha256(a.wall),
            "n_iface": int(len(iface)), "n_layers": int(ns), "n_axial": int(len(iface) - 1), "k_solid": k_s,
            "robin": {"edge": "下面 y = −b", "h": pc.H_BACK, "T_c": Th}, "adiabatic": "両端面"}
    h5, meta, log = axcht.write_solid(out, nodes, tris, outer, back, k_s, pc.H_BACK, Th, f"case/65 {a.case} 板 k_s {k_s:.4g}", meta)
    print(f"[gen_solid] {h5}: {len(nodes)} 節点 / 界面 {len(iface)} / Robin 辺 {len(back)} / 対角 {meta['diagonal']}")


if __name__ == "__main__":
    main()
