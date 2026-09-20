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


def qpp_grid(run, nth=37):
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
        # **等間隔 bin にしない**。節点は層状に並んでいて層厚が 2 桁違うので、等間隔だと
        # セルの 19〜42 % が空欄になり、受け渡しファイルとして使えない (2026-09-20 実測)。
        # **実在する層 (丸めた座標の一意値) で束ねる**ので空欄が出ない。
        if g in ("cav_floor", "cyl_top"):
            # **水平面は半径が層になっていない** (バタフライ / O グリッドなので r が θ で連続に
            # 変わる)。一意値で束ねると 1 点しか入らない列が並ぶ (cyl_top で 999 列・空欄 89 %)。
            # さらに**偏心するとギャップ幅が θ で 5->1 mm と変わる**ので、絶対半径の分位でも
            # θ ごとに範囲が違って埋まらない (偏心の底面で空欄 20 %)。
            # 底面は**すきま内の正規化半径** (0=内円柱壁, 1=外筒壁)、円柱上面は**円柱軸からの
            # 半径を Ri で正規化**して、どの θ でも 0..1 を張る。
            G_ = man["geometry"]
            if g == "cav_floor":
                rr = np.hypot(d["_x_node"], d["_y_node"])
                thc = np.arctan2(np.abs(d["_y_node"]), -d["_x_node"])
                ri = gc.inner_radius_at(np.clip(thc, 0.0, np.pi), man)
                coord = np.clip((rr - ri) / np.maximum(G_["Ro"] - ri, 1e-12), 0.0, 1.0)
                axis, unit = "gap_frac", "0=内円柱壁 / 1=外筒壁"
            else:
                rc = np.hypot(d["_x_node"] - G_["x_off"], d["_y_node"])
                coord = np.clip(rc / max(G_["Ri"], 1e-12), 0.0, 1.0)
                axis, unit = "r_over_Ri", "0=円柱中心 / 1=リップ"
            # 正規化しただけでは足りない: 壁際に節点が密集するので等間隔だと中央が空く。
            # **分位 bin** にし、さらに**空欄が消えるまで bin 数を自動で下げる**
            # (受け渡しファイルに穴を残さないことを優先する)。
            for nb in (16, 12, 10, 8, 6, 4):
                eg = np.unique(np.quantile(coord, np.linspace(0.0, 1.0, nb + 1)))
                if len(eg) < 3:
                    continue
                iz = np.clip(np.digitize(coord, eg) - 1, 0, len(eg) - 2)
                cnt = np.bincount(it * (len(eg) - 1) + iz, minlength=nth * (len(eg) - 1))
                if cnt.min() > 0:
                    break
            zc = 0.5 * (eg[1:] + eg[:-1])
        else:
            # 側壁は押し出しの z 層そのもの (空欄ゼロ)
            coord = d["_z_node"]
            axis, unit = "z", "m"
            zc, iz = np.unique(np.round(coord.astype(np.float64), 9), return_inverse=True)
        nz = len(zc)
        ib = it * nz + iz
        w = d["_w_node"]
        den = np.bincount(ib, weights=w, minlength=nth * nz)
        num = np.bincount(ib, weights=d["_qin_node"] * w, minlength=nth * nz)
        q = (num / np.maximum(den, 1e-30)).reshape(nth, nz)
        q[den.reshape(nth, nz) <= 0] = np.nan
        out[g] = dict(theta_deg=(0.5 * (te[1:] + te[:-1])).tolist(),
                      axis=axis, axis_unit=unit,
                      axis_vals=zc.tolist(),
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
    # **README も書く** (docstring で書くと言っているのに出していなかった)
    G = gc.load_manifest(run=Path(a.runs[0]))["geometry"]
    rd = """# case/49 環状深キャビティ — FEM 受け渡し

形状: 外径 %.1f / 内径 %.1f / 深さ %.1f mm、偏心 %.1f mm。マッハ %.1f、半割モデル。

## 使い方

`Q_i = G_i0 (T_aw - T_i) + sum_j G_ij (T_j - T_i)` + 放射 で渡す。

1. **対流 BC**: 各壁に `T_gas = T_aw = %.1f K` 固定、`h_i = G_i0 / A_i`
   (外筒内壁 %.3f / 円柱側面 %.3f / 底面 %.3f W/m2K)。
2. **壁間結合**: 壁面どうしに線形コンダクタンス `G_ij` を張る (面間熱伝達要素)。
   **これを省くと符号を間違える** — 外筒 1000 degC / 円柱 20 degC のとき外筒の
   正味入熱は **-12.7 W (放熱側)** だが、壁間項を落とすと +37.5 W になる。
3. **放射**: `network.json` の `radiation` の式で面間放射を張る。**線形化しない**。
   放射率 0.4 以上で対流の壁間結合を上回る (0.8 で 4 倍)。
4. **リップ**: `qpp_*.csv` の開口端 %.1f mm は帯積分入熱を保存したまま均してある。
   評価時刻 t に対し `eps <= sqrt(alpha t)` を満たす幅か確認すること
   (秒オーダーなら 1 mm、0.1 s 以下なら 0.5 mm 以下)。

## やってはいけないこと

- **一様壁温の h の包絡を上限として使う**。側壁で 2 倍過大・底面で 5.8 倍過小になり、
  条件によっては符号も逆になる (plan §4.7.6)。
- 壁温 %.0f-%.0f K の外を外挿する。
- 同心の係数を偏心に (またはその逆に) 流用する。底面の `G_i0` が 15 倍違う。

## ファイル

- `network.json` … G_i0 / G_ij / T_aw / 面積 / 放射 / リップ / 適用範囲
- `qpp_<壁>_Tw<温度>K.csv` … q''(theta, z) [W/m2]。行 = theta [deg] (0=上流, 180=下流)、
  列 = z [m] (底面・円柱上面は半径 r [m])。**壁に入る側が正、半割**。

出典 run: %s
""" % (2 * G["Ro"] * 1e3, 2 * G["Ri"] * 1e3, G["depth"] * 1e3, G["x_off"] * 1e3,
       net["mach"], Taw, net["h_i0_W_per_m2K"]["cav_outer"],
       net["h_i0_W_per_m2K"]["cyl_side"], net["h_i0_W_per_m2K"]["cav_floor"],
       a.eps_mm, net["validity"]["wall_T_range_K"][0], net["validity"]["wall_T_range_K"][1],
       ", ".join(net["source_runs"]))
    (out / "README.md").write_text(rd)
    print("wrote %s/network.json, README.md と q'' 分布 CSV" % out)
    print("  G_i0 [W/K]:", {k: round(v, 5) for k, v in net["G_i0_W_per_K"].items()})
    print("  G_ij [W/K]:", {k: round(v, 5) for k, v in net["G_ij_W_per_K"].items()})
    print("  h_i0 [W/m2K]:", {k: round(v, 3) for k, v in net["h_i0_W_per_m2K"].items()})


if __name__ == "__main__":
    main()
