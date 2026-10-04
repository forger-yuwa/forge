#!/usr/bin/env python3
r"""CHT 段 (c) 用の `solid.json` (`cht_loop.py --solid-mode fem2d` の入力) を組む。

中身はすべて**既に検証済みの成果物からの写し**で、ここで新しい値は作らない:
  - 固体メッシュ  : `mesh/solid_<vane>.npz` (**流体の壁節点を外周に使ったもの**を渡すこと。
                    `gen_solid_mesh.py --outer-from <run>/wall_nodes_ordered.csv`)
  - 熱伝導率      : `ref/material_astm310.json` の $k_s(T)$ 表 (外部出典。報告に数値は無い)
  - 孔の Robin    : `ref/run108_internal_bc.json` の**逆算値** ($h_c$, $T_c$)。**公開値ではない**

usage: python3 case/53.c3x_vane_cht/tools/make_solid_json.py [--vane c3x] [--run run108]
                                                             [--npz mesh/solid_c3x.npz] [--Tc 300]
"""
import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASE = {"c3x": ROOT / "case/53.c3x_vane_cht", "markii": ROOT / "case/54.markii_vane_cht"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="c3x")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--npz", default=None)
    ap.add_argument("--Tc", type=float, default=None, help="冷却剤温度 [K] (既定 = 逆算に使った値)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    case = CASE[a.vane]
    bc = json.loads((case / f"ref/{a.run}_internal_bc.json").read_text())
    mat = json.loads((case / "ref/material_astm310.json").read_text())
    npz = Path(a.npz) if a.npz else case / f"mesh/solid_{a.vane}.npz"
    if not npz.is_absolute():
        npz = ROOT / npz
    if not npz.exists():
        raise SystemExit(f"solid mesh not found: {npz}")

    # **孔ごとの T_c** (公開値を使った同定なら 10 個違う値が入る)
    Tc_list = bc["T_c_K"] if isinstance(bc["T_c_K"], list) else [bc["T_c_K"]] * len(bc["h_c_W_m2K"])
    if a.Tc is not None:
        Tc_list = [a.Tc] * len(Tc_list)
    Tc = float(np.mean(Tc_list)) if False else Tc_list
    spec = {
        "_note": "冷却孔の h_c / T_c は報告に無いので公開量から逆算した**推定値**",
        "_source": {"internal_bc": f"ref/{a.run}_internal_bc.json",
                    "material": "ref/material_astm310.json (外部出典)",
                    "mesh": str(npz.relative_to(ROOT))},
        "mesh_npz": str(npz),
        # 出典の表は °C なので K に直す (Fem2DOperator は K で内挿する)
        "k_table": [[t + 273.15 for t in mat["k_table"]["T"]], mat["k_table"]["k"]],
        "holes": [{"h": float(h), "T_c": float(t)}
                  for h, t in zip(bc["h_c_W_m2K"], Tc_list)],
        "T_init": float(np.mean(Tc_list)),
    }
    out = Path(a.out) if a.out else case / f"ref/solid_{a.vane}_{a.run}.json"
    out.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n")
    print(f"[make_solid_json] {len(spec['holes'])} holes, T_c={min(Tc_list):.1f}-{max(Tc_list):.1f} K, "
          f"k_s {min(spec['k_table'][1]):.1f}..{max(spec['k_table'][1]):.1f} W/mK -> "
          f"{out}")


if __name__ == "__main__":
    main()
