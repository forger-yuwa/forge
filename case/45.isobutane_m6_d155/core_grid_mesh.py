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
  python core_grid_mesh.py yaml                                          G1・Gc・G2 の問題 YAML (problem_d155_ns_prod_coldmesh_tw300_cg{1,c,2}.yaml) を書く
  python core_grid_mesh.py post-ab   [--wall-ref DIR] [--g0-res RES_H5] [--band-dir DIR]
                                     §5.1 #2(c) の 0 step の後処理の A/B: G0 の場を列ごとに壁からの距離 d の PCHIP にし、
                                     各格子と参照の細かい格子 Gref に載せて、判定の量の後処理だけの差を測る
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


GREF = (560, 0.004)     # 後処理の A/B の参照格子 (nj, c)。比 ≤ 約 1.03
WINDOWS = ((40.0, 50.0), (65.0, 75.0), (84.0, 94.0))   # θ_r・δ_loc の x の窓 (plan tooling-nozzle-core-grid §4.10)


def radial_fracs_all(prm, wall: Wall):
    """各列の変形前の分布 s (0 = 軸, 1 = 壁) を mesh2d と同じ関数で作り直す。"""
    xs = mesh2d._x_stations(wall.x_in, wall.x_e, prm.ni, prm.throat_refine, prm.throat_width,
                            prm.local_center, prm.local_refine, prm.local_width, prm.x_density_table)
    fr = mesh2d._first_frac_profile(xs, prm)
    S = np.empty((prm.ni, prm.nj)); q = np.empty(prm.ni)
    for i in range(prm.ni):
        if prm.axis_cap_frac is None:
            S[i] = mesh2d._radial_fracs(prm.nj, float(fr[i])); g = np.diff(S[i]); q[i] = float(np.max(g[:-1] / g[1:]))
        else:
            S[i], q[i] = mesh2d._radial_fracs_capfixed(prm.nj, float(fr[i]), float(prm.axis_cap_frac))
    return S, q


def cmd_post_ab(a):
    """G0 の場 (res) を列ごとに d = r_w·(1 − s) の PCHIP で連続にし、各格子の節点に載せて後処理だけを比べる。
    どの格子でも列 i の節点は、壁法線の層の変形 (d だけの関数) で同じ曲線の上に乗るので、d の 1 次元の補間で同じ場を標本化できる。"""
    import h5py
    from scipy.interpolate import Akima1DInterpolator, PchipInterpolator
    sys.path.insert(0, str(HERE))
    import cold_xcheck as XC
    Truth = {"pchip": PchipInterpolator, "akima": Akima1DInterpolator}[a.truth]
    sys.path.insert(0, str(HERE.parents[1] / "design" / "forge_design" / "report"))
    from forge_design.report.nozzle_report import eta_line
    XC.OUTD = Path(a.band_dir)
    yb_x, yb = XC.common_yb()
    wall = Wall(Path(a.wall_ref)); NI, NJ0 = 4719, 121
    keys = ("ro", "Ux", "Uy", "T", "k", "sonic")
    with h5py.File(a.g0_res) as h:
        V0 = {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64).reshape(NI, NJ0) for k in keys}
    prm0 = _mesh_params_from(mesh_block(), wall.scale, NI, NJ0, 1.0)
    S0, _ = radial_fracs_all(prm0, wall)
    xs0 = mesh2d._x_stations(wall.x_in, wall.x_e, NI, prm0.throat_refine, prm0.throat_width, density_table=prm0.x_density_table)
    rw = wall.r(xs0)
    d0 = (rw[:, None] * (1.0 - S0))[:, ::-1]                  # 壁 (d = 0) から軸へ増える
    interp = [{k: Truth(d0[i], V0[k][i, ::-1]) for k in keys} for i in range(NI)]
    x_E = float(json.loads((Path(a.wall_ref) / "prepare_info.json").read_text())["x_E"])
    out = {"x_E": x_E, "ref": {"nj": GREF[0], "c": GREF[1]}, "truth": a.truth}
    arrays = {}
    for name, nj, c in (("G0", None, None), ("G1", 92, 0.03), ("Gc", 122, 0.015), ("G2", 170, 0.015), ("Gref",) + GREF):
        m = mesh_block(nj, c); P, prm = generate(wall, m)
        S, q = radial_fracs_all(prm, wall)
        d = rw[:, None] * (1.0 - S)
        V = {k: np.empty((NI, prm.nj)) for k in keys}
        for i in range(NI):
            for k in keys:
                V[k][i] = interp[i][k](d[i])
        xy = (P.reshape(-1, 2) * wall.scale)
        rec = {"nj": prm.nj, "q_col_max": float(q.max())}
        for prof in ("linear", "pchip"):
            o = XC.reduce_fields(xy, V["ro"].ravel(), V["Ux"].ravel(), V["Uy"].ravel(), V["T"].ravel(), V["k"].ravel(), False,
                                 NI, prm.nj, wall.scale, yb_x, yb, profile=prof)
            sfx = "" if prof == "linear" else "_pchip"
            rec["Q_w" + sfx] = float(o["Q_w"])
            for xq in (40.0, 70.0, 94.0):
                i = int(np.argmin(np.abs(o["x"] - xq)))
                rec[f"theta_r_{int(xq)}{sfx}"] = float(o["theta_r"][i]); rec[f"delta_loc_{int(xq)}{sfx}"] = float(o["delta_loc"][i])
            for (lo, hi) in WINDOWS:              # x の窓の平均 (列の x は不等間隔なので台形の線平均)
                w = (o["x"] >= lo) & (o["x"] <= hi)
                for k in ("theta_r", "delta_loc"):
                    rec[f"{k}_w{int(lo)}_{int(hi)}{sfx}"] = float(np.trapezoid(o[k][w], o["x"][w]) / (o["x"][w][-1] - o["x"][w][0]))
            arrays[(name, prof)] = (o["x"], o["theta_r"], o["delta_loc"])
        M = np.hypot(V["Ux"], V["Uy"]) / V["sonic"]
        F = {"X": P[..., 0], "R": P[..., 1], "V": {"M": M}}
        eta = P[-1, :, 1] / P[-1, -1, 1]; core = (eta >= 0.05) & (eta <= 0.7)
        ee = np.unique(np.concatenate([[0.05], eta[core], [0.7]]))
        rec["exit_M_node_mean"] = float(M[-1][core].mean())
        rec["exit_M_line_mean"] = float(np.trapezoid(np.interp(ee, eta, M[-1]), ee) / 0.65)
        xq = np.linspace(float(P[0, 0, 0]), float(P[-1, 0, 0]), 2401); w = (xq >= x_E + 2) & (xq <= float(P[-1, 0, 0]) - 1)
        for et in (0.0, 0.1):
            dd = 100 * (eta_line(F, "M", et, xq) / 6.0 - 1)
            rec[f"eta{et}_mean_pct"] = float(dd[w].mean())
            rec[f"eta{et}_overshoot_pct"] = float(dd[(xq >= x_E - 15)].max())
        out[name] = rec
        print(f"[post-ab] {name} nj {prm.nj}: Q_w {rec['Q_w']:.6e}、θ_r(70) {rec['theta_r_70']:.6e}、出口 M (線) {rec['exit_M_line_mean']:.7f}", flush=True)
    ref = out["Gref"]; tab = {}
    for name in ("G0", "G1", "Gc", "G2"):
        r = out[name]; tab[name] = {}
        wk = [f"{k}_w{int(lo)}_{int(hi)}{sfx}" for sfx in ("", "_pchip") for k in ("theta_r", "delta_loc") for (lo, hi) in WINDOWS]
        for k in ["Q_w", "theta_r_40", "theta_r_70", "theta_r_94", "delta_loc_40", "delta_loc_70", "delta_loc_94",
                  "theta_r_40_pchip", "theta_r_70_pchip", "theta_r_94_pchip", "delta_loc_40_pchip", "delta_loc_70_pchip", "delta_loc_94_pchip",
                  "exit_M_node_mean", "exit_M_line_mean"] + wk:
            tab[name][k + "_rel_pct"] = 100 * (r[k] / ref[k] - 1)
        for k in ("eta0.0_mean_pct", "eta0.0_overshoot_pct", "eta0.1_mean_pct", "eta0.1_overshoot_pct"):
            tab[name][k + "_diff_pctpt"] = r[k] - ref[k]
    out["vs_ref"] = tab
    # x 方向の誤差の分布 (試験部 [40, 94]、Gref 比): 中央値と最大
    xr = arrays[("Gref", "linear")][0]; w = (xr >= 40) & (xr <= 94); px = {}
    for name in ("G0", "G1", "Gc", "G2"):
        for prof in ("linear", "pchip"):
            for k, idx in (("theta_r", 1), ("delta_loc", 2)):
                e = 100 * (arrays[(name, prof)][idx][w] / arrays[("Gref", prof)][idx][w] - 1)
                px[f"{name}_{prof}_{k}"] = {"median_pct": float(np.median(e)), "max_abs_pct": float(np.abs(e).max()),
                                            "x_at_max": float(xr[w][np.argmax(np.abs(e))])}
    out["x_profile_error"] = px
    return out


YAML_NAMES = {"G1": "cg1", "Gc": "cgc", "G2": "cg2"}
YAML_NJ = {"G1": 92, "Gc": 122, "G2": 170}      # select.json の結果 (2026-10-10、物理壁で全列の q_i ≤ q_max)


def cmd_yaml(a):
    """生産の 300 K の YAML の本文を残し、name・mesh の nj・axis_cap_frac だけを変えた YAML を書く。読んだ辞書が
    その 3 キー以外で生産と一致することを確かめる (違えば書かずに止める)。"""
    src = PROBLEM.read_text(); base = yaml.safe_load(src)
    if "axis_cap_frac" in base["mesh"]:
        raise SystemExit("生産の YAML に axis_cap_frac がある — 想定と違うので止める")
    out = {}
    for arm, tag in YAML_NAMES.items():
        c, qmax, _ = ARMS[arm]; nj = YAML_NJ[arm]
        name = f"{base['name']}_{tag}"
        head = (f"# plan tooling-nozzle-core-grid §4.3・§4.4 の格子 {arm} (c {c}・q_max {qmax}・nj {nj})。core_grid_mesh.py yaml が\n"
                f"# {PROBLEM.name} から作る (手で編集しない)。違いは name・mesh.nj・mesh.axis_cap_frac だけ。\n"
                f"# axis_cap_frac は変形前の分布のパラメータ (壁法線の層のつなぎ目では物理の間隔が最大 25 % 超える)。nj は全列で比 ≤ {qmax} になる最小値\n")
        lines = src.splitlines(); n_name = n_nj = 0; body = []
        for ln in lines:
            if ln.startswith("name: "):
                ln = f"name: {name}"; n_name += 1
            elif ln.startswith("  nj: "):
                body.append(f"  nj: {nj}"); ln = f"  axis_cap_frac: {c}   # 主流と軸付近の間隔 (/局所半径) の上限 (plan tooling-nozzle-core-grid §4.2)"; n_nj += 1
            body.append(ln)
        if n_name != 1 or n_nj != 1:
            raise SystemExit(f"{arm}: name ({n_name}) か nj ({n_nj}) の行が 1 つでない — 止める")
        text = head + "\n".join(body) + "\n"
        new = yaml.safe_load(text)
        b = json.loads(json.dumps(base)); nn = json.loads(json.dumps(new))
        for d in (b, nn):
            d.pop("name"); d["mesh"].pop("nj"); d["mesh"].pop("axis_cap_frac", None)
        if b != nn or new["mesh"]["nj"] != nj or new["mesh"]["axis_cap_frac"] != c or new["name"] != name:
            raise SystemExit(f"{arm}: 生産の YAML との差が name・nj・axis_cap_frac 以外にもある — 止める")
        dst = HERE / f"{name}.yaml"; dst.write_text(text)
        out[arm] = {"file": dst.name, "nj": nj, "axis_cap_frac": c, "identical_except_name_nj_cap": True}
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["check-g0", "select", "geom-ab", "metric-ab", "bl-count", "post-ab", "yaml"])
    ap.add_argument("--wall-ref", default=str(WALL_REF)); ap.add_argument("--g0-res", default=str(HERE / "run_0353_m9_L5" / "res_115000.h5"))
    ap.add_argument("--truth", choices=["pchip", "akima"], default="pchip", help="post-ab: 連続の場の作り方 (G0 の節点の補間)")
    ap.add_argument("--band-dir", default=str(HERE / "_band_ab" / "cold_pair"), help="帯の外縁の npz (theta_run_0181/0183_*.npz) の場所")
    a = ap.parse_args()
    out = {"check-g0": cmd_check_g0, "select": cmd_select, "geom-ab": cmd_geom_ab, "metric-ab": cmd_metric_ab, "bl-count": cmd_bl_count, "post-ab": cmd_post_ab, "yaml": cmd_yaml}[a.cmd](a)
    out["problem"] = PROBLEM.name; out["wall_ref"] = a.wall_ref
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / (f"{a.cmd}_{a.truth}.json" if a.cmd == "post-ab" else f"{a.cmd}.json"); p.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "metrics"}, ensure_ascii=False, indent=1)[:6000])
    print(f"[core_grid_mesh] {p}")


if __name__ == "__main__":
    main()
