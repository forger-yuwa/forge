"""case/45 冷却壁の NS 格子の作り直し (plan tooling-nozzle-core-grid) の格子の設計と 0 step の検査。

格子の系列 (§4.4): 300 K の生産の問題 problem_d155_ns_prod_coldmesh_tw300.yaml の mesh ブロックのうち nj と axis_cap_frac だけを変える。
  G0 = 今 (nj 121、上限なし) / G1 = c 0.03・q ≤ 1.2 / Gc = c 0.015・q ≤ 1.2 / G2 = c 0.015・q ≤ 1.1
  (axis_cap_frac: 上限 c を固定し、列ごとに比 q_i を解く。nj は全列で q_i ≤ q_max になる最小の値)。
壁は元の物理壁 (wall_repr.json) から mesh2d.generate_axisym_mesh で作る (prepare_ns と同じ経路)。

usage (design/.venv-opt の python):
  python core_grid_mesh.py check-g0 [--wall-ref DIR] [--g0-res RES_H5]   G0 が既存の格子の座標を再現するか
  python core_grid_mesh.py select   [--wall-ref DIR]                     G1・Gc・G2 の nj を選び、格子の指標を出す
  python core_grid_mesh.py geom-ab  [--wall-ref DIR]                     §4.9 の 0 step の幾何 A/B (比を固定する形 G1q 対 G1)
  python core_grid_mesh.py metric-ab [--wall-ref DIR] [--g0-res RES_H5]  §4.6 の出口の主流 M の指標の 0 step の A/B
  python core_grid_mesh.py bl-count  [--wall-ref DIR] [--g0-res RES_H5]  G0 の解の境界層の外縁 (μt/μ が最大の 1 %) より内側の節点数
出力: _band_ab/core_grid/<cmd>.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import _mesh_params_from  # noqa: E402
from forge_design.geometry.wall_axismach import load_wall_file  # noqa: E402
from forge_design.meshing import mesh2d  # noqa: E402

PROBLEM = HERE / "problem_d155_ns_prod_coldmesh_tw300.yaml"
WALL_REF = HERE / "_band_ab" / "prod_confirm" / "prep"
OUT = HERE / "_band_ab" / "core_grid"
ARMS = {"G1": (0.03, 1.2, 92), "Gc": (0.015, 1.2, 122), "G2": (0.015, 1.1, 170)}   # (c, q_max, 見積もりの nj)
KINK_DEG = 2.0          # §4.9 の事前基準 (格子の幾何の比較の基準で、CFD 精度の保証ではない)


class Wall:
    """物理壁 (wall_repr.json) を mesh2d の壁の形にする (cold_pair.mesh_checks と同じ)。"""

    def __init__(self, ref: Path):
        W = load_wall_file(ref)
        self.ph = W["physical"]; self.x_in, self.x_e = (float(v) for v in W["domain"]); self.scale = float(W["scale_m"])

    def r(self, x, d=0):
        return self.ph.r(np.asarray(x), d) if d else self.ph.r(np.asarray(x))


def mesh_block(nj=None, cap=None):
    m = dict(yaml.safe_load(open(PROBLEM))["mesh"])
    if nj is not None:
        m["nj"] = int(nj)
    if cap is not None:
        m["axis_cap_frac"] = float(cap)
    return m


def generate(wall: Wall, m: dict):
    prm = _mesh_params_from(m, wall.scale, int(m["ni"]), int(m["nj"]), float(m["wall_first_frac"]))
    coords, quads, _ = mesh2d.generate_axisym_mesh(wall, prm)
    ni, nj = prm.ni, prm.nj
    P = coords[:, :2].reshape(ni, nj, 2) / wall.scale            # r_t 単位、j = 0 が軸
    return P, prm


def column_ratio_max(prm, wall: Wall):
    """各列の分布の隣接比の最大 (変形前、mesh2d が捨てる値を同じ関数で取り直す)。"""
    xs = mesh2d._x_stations(wall.x_in, wall.x_e, prm.ni, prm.throat_refine, prm.throat_width,
                            prm.local_center, prm.local_refine, prm.local_width, prm.x_density_table)
    fr = mesh2d._first_frac_profile(xs, prm)
    q = np.empty(prm.ni)
    for i in range(prm.ni):
        if prm.axis_cap_frac is None:
            g = np.diff(mesh2d._radial_fracs(prm.nj, float(fr[i])))
            q[i] = float(np.max(g[:-1] / g[1:]))
        else:
            _, q[i] = mesh2d._radial_fracs_capfixed(prm.nj, float(fr[i]), float(prm.axis_cap_frac))
    return xs, fr, q


def metrics(P, x_band=(-12.5, -4.0)):
    X, R = P[..., 0], P[..., 1]
    rw = R[:, -1]
    ang = np.arctan2(np.diff(R, axis=0), np.diff(X, axis=0))
    kink = np.degrees(np.abs(np.diff(ang, axis=0)))                # (ni-2, nj) j 一定の格子線の折れ角
    xk = X[1:-1, -1]
    band = (xk >= x_band[0]) & (xk <= x_band[1])
    i, j = np.unravel_index(np.argmax(kink), kink.shape)
    c = [P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]]
    e = np.stack([np.linalg.norm(c[(k + 1) % 4] - c[k], axis=-1) for k in range(4)], -1)
    ar = e.max(-1) / e.min(-1)
    ang4 = []
    for k in range(4):
        a = c[(k - 1) % 4] - c[k]; b = c[(k + 1) % 4] - c[k]
        cs = (a * b).sum(-1) / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1))
        ang4.append(np.degrees(np.arccos(np.clip(cs, -1, 1))))
    A = np.stack(ang4, -1); skew = np.maximum((A.max(-1) - 90) / 90, (90 - A.min(-1)) / 90)
    area = 0.5 * ((c[2][..., 0] - c[0][..., 0]) * (c[3][..., 1] - c[1][..., 1]) - (c[3][..., 0] - c[1][..., 0]) * (c[2][..., 1] - c[0][..., 1]))
    gr = np.hypot(np.diff(X, axis=1), np.diff(R, axis=1))           # 半径方向の辺 (変形後) [r_t]
    gw = gr[:, ::-1]                                               # 壁側から
    qphys = (gw[:, 1:20] / gw[:, :19]).max(1)
    xw = X[:, -1]
    def at(xq, a):
        return float(a[int(np.argmin(np.abs(xw - xq)))])
    return {
        "nj": int(P.shape[1]), "nodes": int(P.shape[0] * P.shape[1]),
        "kink_max_deg": float(kink.max()), "kink_at": [float(xk[i]), int(j)], "kink_p999_deg": float(np.percentile(kink, 99.9)),
        "kink_over_2deg": int((kink > KINK_DEG).sum()),
        "kink_band_max_deg": float(kink[band].max()), "kink_band_over_2deg": int((kink[band] > KINK_DEG).sum()), "band": list(x_band),
        "ar_max": float(ar.max()), "ar_p99": float(np.percentile(ar, 99)), "skew_max": float(skew.max()),
        "area_min_rt2": float(area.min()), "lines_cross": int((np.diff(X, axis=0) <= 0).sum() + (np.diff(R, axis=1) <= 0).sum()),
        "radial_edge_over_rw_max": float((gr / rw[:, None]).max()), "radial_edge_over_rw_axis5_max": float((gr[:, :5] / rw[:, None]).max()),
        "wall_ratio_phys": {"x0": at(0.0, qphys), "x40": at(40.0, qphys), "x94": at(94.0, qphys), "max": float(qphys.max())},
        "axis_gap_rt": {"x0": at(0.0, R[:, 1] - R[:, 0]), "x40": at(40.0, R[:, 1] - R[:, 0]), "x94": at(94.0, R[:, 1] - R[:, 0])},
        "axis_dr_over_dx_max": float(((R[:-1, 1] - R[:-1, 0]) / np.diff(X[:, 0])).max()),
    }


def fixed_q_fracs(q, cmax):
    """§4.9 の A 腕 (比を固定し上限 c_i を列ごとに解く形、廃案の G1q) の分布関数。mesh2d._radial_fracs_capfixed の代わりに差し込む。"""
    def f(nj, first_frac, cap):
        n = nj - 1; k = np.arange(n); fr = float(first_frac)
        tot = lambda c: float(np.minimum(fr * q ** k, c).sum())  # noqa: E731
        if tot(cmax) < 1.0:
            raise ValueError(f"G1q: nj {nj} では c_max {cmax} で和が 1 に届かない (f {fr})")
        lo, hi = fr, cmax
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if tot(mid) < 1.0 else (lo, mid)
        g = np.minimum(fr * q ** k, hi); g = g / g.sum()
        s = np.concatenate([[0.0], np.cumsum(g[::-1])]); s /= s[-1]
        return s, q
    return f


def cmd_check_g0(a):
    wall = Wall(Path(a.wall_ref)); P, prm = generate(wall, mesh_block())
    import h5py
    with h5py.File(a.g0_res) as h:
        C = np.asarray(h["MESH/COORD"][:], dtype=np.float64).reshape(prm.ni, prm.nj, 3)[..., :2] / wall.scale
    d = np.abs(C - P)
    out = {"max_abs_diff_rt": float(d.max()), "first_layer_rel_err_max": float(np.max(np.abs(
        np.linalg.norm(C[:, -1] - C[:, -2], axis=1) / np.linalg.norm(P[:, -1] - P[:, -2], axis=1) - 1))), "g0_res": str(a.g0_res)}
    out["metrics"] = metrics(P)
    return out


def cmd_select(a):
    wall = Wall(Path(a.wall_ref)); out = {}
    P0, prm0 = generate(wall, mesh_block()); out["G0"] = metrics(P0)
    _, fr, q0 = column_ratio_max(prm0, wall); out["G0"]["q_col_max"] = float(q0.max())
    def qcol(nj, c):
        prm = _mesh_params_from(mesh_block(nj, c), wall.scale, 4719, nj, 1.0)
        try:
            return column_ratio_max(prm, wall)[2]
        except ValueError:
            return None
    for name, (c, qmax, nj0) in ARMS.items():
        ok = lambda q: q is not None and q.max() <= qmax + 1e-12  # noqa: E731
        nj = nj0
        while not ok(qcol(nj, c)):
            nj += 1
        while nj > 3 and ok(qcol(nj - 1, c)):
            nj -= 1
        q = qcol(nj, c)
        P, _ = generate(wall, mesh_block(nj, c))
        out[name] = metrics(P); out[name].update(c=c, q_max_target=qmax, q_col_max=float(q.max()), q_col_min=float(q.min()))
        print(f"[select] {name}: nj {nj}、列の比 {q.min():.4f}〜{q.max():.4f}、折れ角 最大 {out[name]['kink_max_deg']:.2f}°", flush=True)
    return out


def cmd_geom_ab(a):
    wall = Wall(Path(a.wall_ref)); out = {}
    P, _ = generate(wall, mesh_block(92, 0.03)); out["B_G1"] = metrics(P)
    orig = mesh2d._radial_fracs_capfixed
    try:
        mesh2d._radial_fracs_capfixed = fixed_q_fracs(1.2, 0.03)
        P, _ = generate(wall, mesh_block(92, 0.03)); out["A_G1q"] = metrics(P)
    finally:
        mesh2d._radial_fracs_capfixed = orig
    A, B = out["A_G1q"], out["B_G1"]
    support = A["kink_band_over_2deg"] > 0 and B["kink_over_2deg"] == 0
    out["verdict"] = ("支持 (A に 2° 超があり、B は全域で 2° 以下)" if support else
                      "不支持 (採用の根拠を再検証)")
    return out


def cmd_metric_ab(a):
    """G0 の出口の M(η) を連続の分布 (PCHIP) にし、各格子の最終断面の節点の η に載せて、節点平均 A と線平均 B を連続の平均と比べる。"""
    import h5py
    from scipy.interpolate import PchipInterpolator
    wall = Wall(Path(a.wall_ref)); NI, NJ0 = 4719, 121
    with h5py.File(a.g0_res) as h:
        Rx = np.asarray(h["MESH/COORD"][:], dtype=np.float64).reshape(NI, NJ0, 3)[-1, :, 1]
        V = {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64).reshape(NI, NJ0)[-1] for k in ("Ux", "Uy", "sonic")}
    f = PchipInterpolator(Rx / Rx[-1], np.hypot(V["Ux"], V["Uy"]) / V["sonic"])
    e = np.linspace(0.05, 0.7, 200001); true = float(np.trapezoid(f(e), e) / 0.65)
    out = {"true_mean": true}
    for name, nj, c in (("G0", None, None), ("G1", 92, 0.03), ("Gc", 122, 0.015), ("G2", 170, 0.015)):
        P, _ = generate(wall, mesh_block(nj, c))
        eta = P[-1, :, 1] / P[-1, -1, 1]; M = f(eta); core = (eta >= 0.05) & (eta <= 0.7)
        ee = np.unique(np.concatenate([[0.05], eta[core], [0.7]]))
        A = float(M[core].mean()); B = float(np.trapezoid(np.interp(ee, eta, M), ee) / 0.65)
        out[name] = {"n_core": int(core.sum()), "A_node_mean": A, "B_line_mean": B, "A_minus_true_pct": 100 * (A / true - 1), "B_minus_true_pct": 100 * (B / true - 1)}
    return out


def cmd_bl_count(a):
    """G0 の解 (res) で列ごとに μt/μ が最大値の 1 % に落ちる外縁を取り、各格子の同じ列で壁からその距離までの節点を数える。
    x ≤ −4 は入口から乱れが入っていて外縁が決まらない (壁から軸まで 1 % を下回らない列は数えない)。"""
    import h5py
    wall = Wall(Path(a.wall_ref)); NI, NJ0 = 4719, 121
    with h5py.File(a.g0_res) as h:
        C = np.asarray(h["MESH/COORD"][:], dtype=np.float64).reshape(NI, NJ0, 3)[..., :2] / wall.scale
        q = (np.asarray(h["VALUE/vis_turb"][:], dtype=np.float64) / np.asarray(h["VALUE/vis_lam"][:], dtype=np.float64)).reshape(NI, NJ0)
    dist0 = np.linalg.norm(C - C[:, -1:, :], axis=2)
    delta = np.full(NI, np.nan)
    for i in range(NI):
        jp = int(np.argmax(q[i])); out = np.where(q[i, :jp] < 0.01 * q[i, jp])[0]
        if len(out):
            delta[i] = dist0[i, out.max()]
    xw = C[:, -1, 0]; res = {"defined_columns": int(np.isfinite(delta).sum())}
    for name, nj, c in (("G0", None, None), ("G1", 92, 0.03), ("Gc", 122, 0.015), ("G2", 170, 0.015)):
        P, _ = generate(wall, mesh_block(nj, c))
        d = np.linalg.norm(P - P[:, -1:, :], axis=2)
        n = np.where(np.isfinite(delta), (d <= delta[:, None]).sum(1), -1)
        row = {}
        for xq in (-2.0, 0.0, 5.0, 40.0, 94.0):
            i = int(np.argmin(np.abs(xw - xq))); row[f"x{xq:g}"] = int(n[i])
        ok = n >= 0
        row["min"], row["max"] = int(n[ok].min()), int(n[ok].max())
        res[name] = row
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["check-g0", "select", "geom-ab", "metric-ab", "bl-count"])
    ap.add_argument("--wall-ref", default=str(WALL_REF)); ap.add_argument("--g0-res", default=str(HERE / "run_0353_m9_L5" / "res_115000.h5"))
    a = ap.parse_args()
    out = {"check-g0": cmd_check_g0, "select": cmd_select, "geom-ab": cmd_geom_ab, "metric-ab": cmd_metric_ab, "bl-count": cmd_bl_count}[a.cmd](a)
    out["problem"] = PROBLEM.name; out["wall_ref"] = a.wall_ref
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / f"{a.cmd}.json"; p.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "metrics"}, ensure_ascii=False, indent=1)[:6000])
    print(f"[core_grid_mesh] {p}")


if __name__ == "__main__":
    main()
