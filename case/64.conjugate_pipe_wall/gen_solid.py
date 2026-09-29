#!/usr/bin/env python3
r"""case/64 の固体殻 (R ≤ r ≤ r_o、全長)。外面の Robin (h_o、T_c = T_in + 10 K) は**加熱区間 0 ≤ x ≤ L_h の辺だけ**、それ以外の外面と両端面は断熱。
界面節点は乾式 1 step の壁ダンプ `res_wall_3_1.h5` の座標そのもの。帯トポロジは case/61・63 と同じ (行 j=0 が外面)。
半径方向の層数は流体と同じ間隔になるよう N_r (--ns 既定 = 流体の N_r)。

    python3 gen_solid.py --wall run_0001_dry_r32/res_wall_3_1.h5 --case A1
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pipe_common as pc
from pipe_common import axcht, gc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wall", required=True); ap.add_argument("--case", choices=sorted(pc.KS_RATIO), required=True)
    ap.add_argument("--ns", type=int, default=None)
    a = ap.parse_args()
    c = axcht.read_wall_coords(Path(a.wall))
    if not np.allclose(c[:, 1], pc.R, rtol=0, atol=1e-9):
        raise SystemExit("壁ダンプが r = R の直線でない")
    order = np.argsort(c[:, 0]); iface = c[order, :2]
    if abs(iface[0, 0] + pc.L_UP) > 1e-7 or abs(iface[-1, 0] - pc.L_HEAT - pc.L_DOWN) > 1e-7:
        raise SystemExit("壁ダンプの x 範囲が全長でない")
    nax = len(iface) - 1
    tag = Path(a.wall).resolve().parent.name
    nr = int(tag.split("_r")[-1]) if "_r" in tag else 32
    ns = a.ns or nr
    rc = float(iface[0, 1]) + (pc.R_O - pc.R)
    nodes, tris, outer, back = axcht.solid_strip(iface, normal_axis=1, back_coord=rc, ns=ns)
    if axcht.diag_direction(nodes, len(iface), ns) != "左上→右下":
        raise SystemExit("対角が登録と違う")
    # 外面の辺のうち加熱区間 (両端点が 0 ≤ x ≤ L_h) だけを Robin の孔に
    xb = nodes[:, 0]
    heat = np.array([(-1e-9 <= xb[e[0]] <= pc.L_HEAT + 1e-7) and (-1e-9 <= xb[e[1]] <= pc.L_HEAT + 1e-7) for e in back])
    back_h = back[heat]
    k_s = pc.KS_RATIO[a.case] * gc.K_F
    Tc = gc.T_IN + pc.DT_C
    out = pc.HERE / "mesh" / f"solid_{tag}_{a.case}"
    meta = {"case": f"case/64 {a.case} 固体殻 (全長)", "wall_dump": str(Path(a.wall)), "wall_dump_sha256": axcht.sha256(a.wall),
            "n_iface": int(len(iface)), "n_layers": int(ns), "n_axial": int(nax),
            "radial": {"from": float(iface[0, 1]), "to": rc, "layers": ns},
            "k_solid": k_s, "robin": {"edges": "外面のうち 0 ≤ x ≤ L_h", "n_edges": int(len(back_h)), "h": pc.H_O, "T_c": Tc},
            "adiabatic": "加熱区間外の外面と両端面"}
    note = f"case/64 {a.case}: k_s {k_s:.6g}、外面 Robin h {pc.H_O} / T_c {Tc} K を加熱区間の {len(back_h)} 辺だけ"
    h5, meta, log = axcht.write_solid(out, nodes, tris, outer, back_h, k_s, pc.H_O, Tc, note, meta)
    print(f"[gen_solid] {h5}: {len(nodes)} 節点 / 界面 {len(iface)} / Robin 辺 {len(back_h)} (外面 {len(back)} 辺のうち) / k_s {k_s:.4g}")


if __name__ == "__main__":
    main()
