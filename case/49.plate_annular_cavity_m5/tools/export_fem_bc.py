#!/usr/bin/env python3
"""非定常熱構造 FEM へ渡す境界条件一式を書き出す (plan §4.7.7 の正式案)。

§4.7.6 で「等温壁 3 本の $h_{aw}$ 包絡は上限にならない (符号すら外す)」ことが確定したので、
**熱回路**で渡す:

    Q_i = G_i0 (T_aw - T_i) + sum_j G_ij (T_j - T_i)   [W]   + 放射

FEM 側の入れ方:
  1. 各壁に T_gas = T_aw 固定・h_i = G_i0 / A_i の対流 BC
  2. 壁面どうしに線形コンダクタンス G_ij (面間熱伝達要素)
  3. 壁面どうしに放射 (T^4。線形化しない)
  4. リップ帯は帯積分入熱を保存したまま幅 eps に均す。eps <= sqrt(alpha t) で選ぶ

出力 (既定 `fem_bc/`):
  network.json      G_i0 / G_ij / T_aw / 面積 / 放射の形態係数パラメータ / リップ帯
  qpp_<壁>.csv      壁温ごとの q''(theta, z) 分布 (FEM に直接マップする用)
  README.md         使い方と適用範囲

usage:
  python3 tools/export_fem_bc.py --runs RUN_20 RUN_500 RUN_1000 \
      --mixed RUN_MIXA RUN_MIXB [--out fem_bc]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import cavity_eval as ce  # noqa: E402
import geom_common as gc  # noqa: E402
from fit_wall_network import WALLS, PAIRS, JP, read_run, build  # noqa: E402

SIGMA = 5.670374419e-8


def qpp_grid(run, nth=37, nz=26):
    """壁ごとに q''(theta, z) を (theta, z) 格子へ束ねる。半割なので theta は 0..180 度。"""
    run = Path(run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    step = int(ce.snapshots(run)[-1].stem.split("_")[1])
    wh = ce.wall_heat(run, step, man, D)
    dep = float(man["geometry"]["depth"])
    te = np.linspace(0.0, 180.0, nth + 1)
    out = {}
    for g in WALLS + ["cyl_top"]:
        if g not in wh:
            continue
        d = wh[g]
        th = np.degrees(np.arctan2(np.abs(d["_y_node"]), -d["_x_node"]))
        it = np.clip(np.digitize(th, te) - 1, 0, nth - 1)
        if g == "cav_floor" or g == "cyl_top":
            # 水平面は半径方向に束ねる
            rr = np.hypot(d["_x_node"], d["_y_node"])
            ze = np.linspace(rr.min(), rr.max(), nz + 1)
            iz = np.clip(np.digitize(rr, ze) - 1, 0, nz - 1)
            axis, unit = "r", "m"
        else:
            ze = np.linspace(-dep, 0.0, nz + 1)
            iz = np.clip(np.digitize(d["_z_node"], ze) - 1, 0, nz - 1)
            axis, unit = "z", "m"
        ib = it * nz + iz
        w = d["_w_node"]
        den = np.bincount(ib, weights=w, minlength=nth * nz)
        num = np.bincount(ib, weights=d["_qin_node"] * w, minlength=nth * nz)
        q = (num / np.maximum(den, 1e-30)).reshape(nth, nz)
        q[den.reshape(nth, nz) <= 0] = np.nan
        out[g] = dict(theta_deg=(0.5 * (te[1:] + te[:-1])).tolist(),
                      axis=axis, axis_unit=unit,
                      axis_vals=(0.5 * (ze[1:] + ze[:-1])).tolist(),
                      qpp_W_m2=q.tolist(), area_m2=float(d["area_m2"]),
                      Q_W_half=float(d["Q_W"]))
    return D, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="一様壁温の run (3 点以上)")
    ap.add_argument("--mixed", nargs="+", required=True, help="非一様壁温の run (2 本以上)")
    ap.add_argument("--out", default="fem_bc")
    ap.add_argument("--eps-mm", type=float, default=1.0, help="リップ帯の幅 [mm]")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(exist_ok=True)

    rows = [read_run(r) for r in list(a.runs) + list(a.mixed)]
    A, b = build(rows)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    G0, Gij = x[:3], x[3:]
    Taw = rows[0]["Taw"]

    # 面積は一様 run の 1 本から (半割)
    D0, grids0 = qpp_grid(a.runs[0])
    areas = {g: grids0[g]["area_m2"] for g in grids0}

    net = dict(
        _comment="case/49 環状深キャビティの FEM 受け渡し。Q_i = G_i0 (T_aw - T_i) + sum_j G_ij (T_j - T_i) + 放射。"
                 "**包絡 (一様壁の h_aw の最大) は上限にならない** (plan §4.7.6)。",
        T_aw_K=Taw, T_t_K=float(D0["Tt_tp"] if D0.get("gas_used") == "TP" else D0["Tt_cpg"]),
        mach=float(D0["mach"]), half_model=True,
        note_half="G は半割の値。全周は 2 倍にすること",
        G_i0_W_per_K={WALLS[i]: float(G0[i]) for i in range(3)},
        h_i0_W_per_m2K={WALLS[i]: float(G0[i] / areas[WALLS[i]]) for i in range(3)},
        G_ij_W_per_K={"%s-%s" % (WALLS[p], WALLS[q]): float(Gij[k])
                      for k, (p, q) in enumerate(PAIRS)},
        area_m2_half=areas,
        radiation=dict(
            form="長い同軸 2 円筒: q1 = sigma (T1^4 - T2^4) / (1/e1 + (A1/A2)(1/e2 - 1))",
            inner="cyl_side", outer="cav_outer",
            A_inner_m2=areas.get("cyl_side"), A_outer_m2=areas.get("cav_outer"),
            note="eps>=0.4 で対流の G_ij を上回る (eps=0.8 で 4 倍)。線形化せず T^4 で入れること",
        ),
        lip=dict(band_mm=a.eps_mm,
                 rule="eps <= sqrt(alpha t) で選ぶ。秒オーダーなら 1 mm、0.1 s 以下は 0.5 mm 以下",
                 note="帯積分入熱を保存したまま均す。ピーク q'' は格子収束しない (plan §4.4.2)"),
        validity=dict(
            geometry="Ro %.1f / Ri %.1f / depth %.1f / x_off %.1f mm"
                     % tuple(float(gc.load_manifest(run=Path(a.runs[0]))["geometry"][k]) * 1e3
                             for k in ("Ro", "Ri", "depth", "x_off")),
            wall_T_range_K=[min(r["T"].min() for r in rows), max(r["T"].max() for r in rows)],
            note="この形状・マッハ数・同心/偏心の別でのみ有効。外挿しない",
        ),
        source_runs=[r["name"] for r in rows],
    )
    (out / "network.json").write_text(json.dumps(net, indent=2, ensure_ascii=False))

    for r in a.runs:
        D, grids = qpp_grid(r)
        tw = float(D["wall_T"])
        for g, gd in grids.items():
            f = out / ("qpp_%s_Tw%03.0fK.csv" % (g, tw))
            lines = ["# %s  Tw=%.2f K  q'' [W/m2] (壁に入る側が正)  半割" % (g, tw),
                     "# 行 = theta [deg] (0=上流, 180=下流), 列 = %s [%s]"
                     % (gd["axis"], gd["axis_unit"]),
                     "theta_deg," + ",".join("%.6g" % v for v in gd["axis_vals"])]
            for i, th in enumerate(gd["theta_deg"]):
                lines.append("%.4g," % th + ",".join(
                    ("" if not np.isfinite(v) else "%.6g" % v) for v in gd["qpp_W_m2"][i]))
            f.write_text("\n".join(lines) + "\n")
    print("wrote %s/network.json と q'' 分布 CSV" % out)
    print("  G_i0 [W/K]:", {k: round(v, 5) for k, v in net["G_i0_W_per_K"].items()})
    print("  G_ij [W/K]:", {k: round(v, 5) for k, v in net["G_ij_W_per_K"].items()})
    print("  h_i0 [W/m2K]:", {k: round(v, 3) for k, v in net["h_i0_W_per_m2K"].items()})


if __name__ == "__main__":
    main()
