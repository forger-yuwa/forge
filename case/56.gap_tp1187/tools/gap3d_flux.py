#!/usr/bin/env python3
"""3D すきまの給気の測定 — 口と衝突域の直前の断面で、質量流束と過剰全エンタルピー流束を面積積分する (case/56)。

    python3 tools/gap3d_flux.py <run_dir> [--step N] [--r-edge 2.5e-3] [--n 120]

acceptance.json の T4-3D-RADIUS (2026-09-29) の測定量。**A/B の両腕で同じ断面・同じ範囲**を使う:

- P1 (縦すきまの口): 平面 y = 0、x ∈ [x2d, xu] (上流横すきまの下流壁 〜 主横すきまの上流壁)、
  z ∈ [0, W/2 + R] (R = --r-edge。鋭い腕では z > W/2 はタイル上面なので流束 0)。
  流入 = すきまへ入る向き (−y)、流出 = すきまから出る向き (+y)。
- P2 (衝突域の直前): 平面 x = 0 (主横すきまの中心)、y ∈ [−D, 0]、z ∈ [0, W/2 + R]。
  流入 = 前向き壁へ向かう向き (+x)。
- 過剰全エンタルピー流束 = ∫ ρ u_n (h0 − h_w) dA、h_w = 壁温・入口組成のエンタルピー (NASA-9、run の species_db.yaml)。
  h0 は出力の VALUE/h0 (全エンタルピー、thermoHrefTemp 基準)。h_w も同じ基準 (h(Tw) − h(Tref)) にそろえる。

場は節点値なので、断面上の規則格子へ k 近傍の距離逆数重みで写す (近似。両腕で同じ手順なので差の比較には使える)。
対称面 z = 0 で切った半ピッチの領域なので、積分値は半分 (z ≥ 0) の量。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py
import yaml
from scipy.spatial import cKDTree

RU = 8.314462618


def h_mix(db, species, Y, T):
    h = 0.0
    for s, y in zip(species, Y):
        c = db[s]; a = c["nasa9_low"] if T < c["Tmid"] else c["nasa9_high"]
        hs = RU / c["MW"] * T * (-a[0] * T ** -2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3
                                 + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + a[7] / T)
        h += y * hs
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--r-edge", type=float, default=2.5e-3, help="断面の z 範囲 W/2 + R の R (両腕で同じ値)")
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2, help="横すきまの縁の半径")
    ap.add_argument("--l-long", type=float, default=15.24e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--n", type=int, default=120, help="断面の格子の各方向の点数")
    ap.add_argument("--wall-r", type=float, default=None,
                    help="その run の縦すきま口の実際の半径 [m] (鋭い腕 0、丸み腕 2.5e-3)。P1 面のうち"
                         "タイル上面 (固体) に当たる点の流束を 0 にするのに使う。省略時は mesh 名から推定しない → 必須")
    ap.add_argument("--k", type=int, default=8)
    a = ap.parse_args()
    if a.wall_r is None:
        raise SystemExit("--wall-r (その run の縦すきま口の半径) を指定する")
    run = Path(a.run)
    steps = sorted(int(re.search(r"res_(\d+)\.h5$", f).group(1)) for f in glob.glob(str(run / "res_[0-9]*.h5")))
    step = a.step if a.step is not None else steps[-1]
    with h5py.File(run / f"res_{step}.h5", "r") as h:
        X = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        V = h["VALUE"]
        ro = np.asarray(V["ro"], float); U = np.column_stack([np.asarray(V[k], float) for k in ("Ux", "Uy", "Uz")])
        h0 = np.asarray(V["h0"], float)
    su = json.loads((run / "case_setup.json").read_text())
    db = yaml.safe_load(open(run / "species_db.yaml"))
    species = list(su["Y"].keys()) if isinstance(su["Y"], dict) else None
    Y = [su["Y"][s] for s in species]
    # forge の h0 は thermoHrefTemp (既定 298.15 K) 基準 (NASA 絶対値 − h(Tref))。壁の h も同じ基準にそろえる
    tref = 298.15
    m_ = re.search(r"thermoHrefTemp:\s*([0-9.eE+-]+)", (run / "solverConfig.yaml").read_text())
    if m_:
        tref = float(m_.group(1))
    hw = h_mix(db, species, Y, float(su.get("Tw", 300.0))) - h_mix(db, species, Y, tref)
    W, R, D = a.w, a.r_edge, a.depth
    xu, x2d = -0.5 * W, -a.l_long
    zmax = 0.5 * W + R
    tree = cKDTree(X)

    def sample(P):
        d, i = tree.query(P, k=a.k)
        w = 1.0 / np.maximum(d, 1e-12) ** 2; w /= w.sum(1, keepdims=True)
        return (w * ro[i]).sum(1), np.einsum("pk,pkc->pc", w, U[i]), (w * h0[i]).sum(1)

    out = {}
    # P1: y = 0、x ∈ [x2d, xu]、z ∈ [0, zmax]  (x は端を少し外して横すきまの円弧を避けない — 同じ範囲なら差は比べられる)
    xs = np.linspace(x2d, xu, 4 * a.n); zs = np.linspace(0.0, zmax, a.n)
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    P = np.column_stack([XX.ravel(), np.zeros(XX.size), ZZ.ravel()])
    r_, u_, h_ = sample(P)
    dA = (xs[1] - xs[0]) * (zs[1] - zs[0])
    mdot = r_ * (-u_[:, 1])                       # すきまへ入る向きを正
    # 固体面の除外 (codex 2026-09-29 Minor: 近傍補間は壁の不透過を保証せず、鋭い腕で P1 流入質量の 1.04 % が
    # 固体上面を通る偽の流入だった)。y = 0 の面がタイル上面 (固体) に当たるのは、上面が平らな区間
    # x ∈ [x2vd, xvu] で z ≥ W/2 + R_arm のところ。横すきまの円弧区間では上面が y < 0 なので面は流体。
    xvu_, x2vd_ = xu - a.r, x2d + a.r            # 生成器と同じ (x2vd = x2d + r、xvu = xu − r)
    solid = (P[:, 0] >= x2vd_) & (P[:, 0] <= xvu_) & (P[:, 2] >= 0.5 * W + a.wall_r - 1e-12)
    mdot = np.where(solid, 0.0, mdot)
    n_solid = int(solid.sum())
    out["P1_mouth"] = dict(n_solid_points_zeroed=n_solid, m_in=float(mdot[mdot > 0].sum() * dA), m_out=float(-mdot[mdot < 0].sum() * dA),
                           H_in=float((mdot * (h_ - hw))[mdot > 0].sum() * dA),
                           H_out=float(-(mdot * (h_ - hw))[mdot < 0].sum() * dA))
    # P2: x = 0、y ∈ [−D, 0]、z ∈ [0, zmax]
    ys = np.linspace(-D, 0.0, 4 * a.n)
    YY, ZZ = np.meshgrid(ys, zs, indexing="ij")
    P = np.column_stack([np.zeros(YY.size), YY.ravel(), ZZ.ravel()])
    r_, u_, h_ = sample(P)
    dA = (ys[1] - ys[0]) * (zs[1] - zs[0])
    mdot = r_ * u_[:, 0]                          # 前向き壁へ向かう向きを正
    prof = mdot.reshape(len(ys), len(zs)).sum(1) * (zs[1] - zs[0])     # 深さごとの単位高さあたり質量流量
    out["P2_preimpinge"] = dict(m_in=float(mdot[mdot > 0].sum() * dA), m_out=float(-mdot[mdot < 0].sum() * dA),
                                H_in=float((mdot * (h_ - hw))[mdot > 0].sum() * dA),
                                H_out=float(-(mdot * (h_ - hw))[mdot < 0].sum() * dA),
                                depth_profile_m_per_m={f"{-yy*1e3:.2f}mm": float(v) for yy, v in
                                                       list(zip(ys, prof))[::-max(1, len(ys) // 16)]})
    out["meta"] = dict(run=str(run), step=step, h_w=hw, W=W, R_section=R, zmax=zmax, x2d=x2d, xu=xu, depth=D,
                       note="半ピッチ (z ≥ 0) の量、単位 kg/s・W")
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
