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
  python core_grid_mesh.py prep <G1|Gc|G2>                              (AWS) §5.1 #3: _prepare_ns で格子・変換・初期値の移送を
                                     _band_ab/core_grid/prep_<格子>/ に作り、品質・第一層・壁・法線・壁距離・初期値を検査する
                                     (FORGE_BIN・REAL_CONVERTER・FORGE_CONVERTER は run_cold_pair.sh と同じに設定して呼ぶ)
  python core_grid_mesh.py prep-compare                                  (AWS) prep_Gc と prep_G2 で座標が一致する範囲 (列ごと)
  python core_grid_mesh.py view [--g0-res RES_H5]                       ParaView で見る表示用ファイル (G0・G1・G1x) を _band_ab/core_grid/view/ に書く
  python core_grid_mesh.py prep-ic                                       (AWS) 移送後の初期値 (nozzle.h5 の保存量) の非有限・正値と、
                                     新しい nj で reduce_fields が通ること (θ_r・δ_loc を G0 の res_100000 と並べる)
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
ARMS = {"G1": (0.03, 1.3, 75), "Gc": (0.015, 1.2, 122), "G2": (0.015, 1.1, 170)}   # (c, q_max, 見積もりの nj)。G1 は 2026-10-10 に比 1.2 → 1.3 (ユーザ「壁際の比も 1.2〜1.3」)
KINK_DEG = 2.0          # §4.9 の事前基準 (格子の幾何の比較の基準で、CFD 精度の保証ではない)
# G1x (2026-10-10 ユーザ決定「1 万上限にしましょう」): 壁際の AR の上限 1 万。今の格子の x の間隔は AR の目標 4500 (上限 5000) で決まっているので、
# 同じ余裕で目標 9000 にする = AR で決まっている範囲の密度を 1/2 にする (今の表の相対密度 d に対し max(1, d/XFACTOR); 1 は試験部の上限 0.06 r_t)
XFACTOR = 2.0
AR_MAX_G1X = 10000


class Wall:
    """物理壁 (wall_repr.json) を mesh2d の壁の形にする (cold_pair.mesh_checks と同じ)。"""

    def __init__(self, ref: Path):
        W = load_wall_file(ref)
        self.ph = W["physical"]; self.x_in, self.x_e = (float(v) for v in W["domain"]); self.scale = float(W["scale_m"])

    def r(self, x, d=0):
        return self.ph.r(np.asarray(x), d) if d else self.ph.r(np.asarray(x))


def x_table_coarse(m, factor=XFACTOR):
    """今の x_density_table を、AR で決まっている範囲 (相対密度 > 1) だけ 1/factor にした表と、密度の積分から決め直した ni。"""
    t = np.asarray(m["x_density_table"], dtype=float)
    new = np.c_[t[:, 0], np.maximum(1.0, t[:, 1] / factor)]
    xs = np.linspace(t[0, 0], t[-1, 0], 200001)
    ratio = np.trapezoid(np.interp(xs, new[:, 0], new[:, 1]), xs) / np.trapezoid(np.interp(xs, t[:, 0], t[:, 1]), xs)
    ni = int(np.ceil((int(m["ni"]) - 1) * ratio)) + 1
    return [[float(a), float(f"{b:.6g}")] for a, b in new], ni


def mesh_block(nj=None, cap=None, xcoarse=False):
    m = dict(yaml.safe_load(open(PROBLEM))["mesh"])
    if nj is not None:
        m["nj"] = int(nj)
    if cap is not None:
        m["axis_cap_frac"] = float(cap)
    if xcoarse:
        m["x_density_table"], m["ni"] = x_table_coarse(m)
        m["ar_max"] = AR_MAX_G1X
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


def cmd_select_g1x(a):
    """G1x (G1 の半径方向 + x 方向を粗く) の格子の指標。列の比 ≤ q_max と AR ≤ 1 万を確かめる。"""
    wall = Wall(Path(a.wall_ref)); c, qmax, _ = ARMS["G1"]; nj = YAML_NJ["G1"]
    m = mesh_block(nj, c, xcoarse=True)
    prm = _mesh_params_from(m, wall.scale, int(m["ni"]), nj, 1.0)
    xs, fr, q = column_ratio_max(prm, wall)
    P, _ = generate(wall, m); out = {"G1x": metrics(P)}
    out["G1x"].update(ni=int(m["ni"]), nj=nj, c=c, q_col_min=float(q.min()), q_col_max=float(q.max()), xfactor=XFACTOR,
                      ar_ok=bool(out["G1x"]["ar_max"] <= AR_MAX_G1X), q_ok=bool(q.max() <= qmax + 1e-12))
    dx = np.diff(P[:, -1, 0]); xw = P[:-1, -1, 0]
    out["G1x"]["dx_wall_rt"] = {f"x{v:g}": float(dx[int(np.argmin(np.abs(xw - v)))]) for v in (-8.0, -5.0, -2.0, 0.0, 2.0, 5.0, 12.0, 40.0)}
    return out


REF_CAP = 0.015          # Gref の主流の上限 (plan §4.16)


def cmd_select_ref(a):
    """Gref (plan §4.16): 主流の上限 0.015 r_w、壁際の比は各列で G0 のその列の比以下 (境界層は G0 と同じかそれより細かい)。条件を満たす最小の nj。"""
    wall = Wall(Path(a.wall_ref))
    prm0 = _mesh_params_from(mesh_block(), wall.scale, 4719, 121, 1.0)
    _, _, q0 = column_ratio_max(prm0, wall)
    def ok(nj):
        prm = _mesh_params_from(mesh_block(nj, REF_CAP), wall.scale, 4719, nj, 1.0)
        try:
            q = column_ratio_max(prm, wall)[2]
        except ValueError:
            return False, None
        return bool(np.all(q <= q0 + 1e-12)), q
    nj = 160
    while not ok(nj)[0]:
        nj += 1
    while ok(nj - 1)[0]:
        nj -= 1
    _, q = ok(nj)
    P, _ = generate(wall, mesh_block(nj, REF_CAP)); m = metrics(P)
    m.update(nj=nj, c=REF_CAP, q_col_min=float(q.min()), q_col_max=float(q.max()), q_over_G0_max=float((q / q0).max()))
    print(f"[select-ref] Gref: nj {nj}、節点 {m['nodes']}、列の比 {q.min():.4f}〜{q.max():.4f} (G0 比の最大 {m['q_over_G0_max']:.4f})、折れ角 {m['kink_max_deg']:.2f}°、AR {m['ar_max']:.0f}、skew {m['skew_max']:.3f}", flush=True)
    return {"Gref": m}


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
    # 2026-10-10 の記録 (旧 G1 = c 0.03・比 ≤ 1.2・nj 92 で、分布の形だけを替えた比較)。今の G1 (比 ≤ 1.3・nj 75) ではない
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
    for name, nj, c in grid_list():
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
    for name, nj, c in grid_list():
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
    for name, nj, c in grid_list(ref=True):
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


YAML_NAMES = {"G1": "cg1", "Gc": "cgc", "G2": "cg2", "Gref": "cgref"}
YAML_NJ = {"G1": 75, "Gc": 122, "G2": 170, "Gref": 160}      # select.json の結果 (2026-10-10、物理壁で全列の q_i ≤ q_max)


def grid_list(ref=False):
    g = [("G0", None, None)] + [(a, YAML_NJ[a], ARMS[a][0]) for a in ("G1", "Gc", "G2")]
    return g + [("Gref",) + GREF] if ref else g


def cmd_yaml(a):
    """生産の 300 K の YAML の本文を残し、name・mesh の nj・axis_cap_frac だけを変えた YAML を書く。読んだ辞書が
    その 3 キー以外で生産と一致することを確かめる (違えば書かずに止める)。"""
    src = PROBLEM.read_text(); base = yaml.safe_load(src)
    if "axis_cap_frac" in base["mesh"]:
        raise SystemExit("生産の YAML に axis_cap_frac がある — 想定と違うので止める")
    out = {}
    for arm, tag in YAML_NAMES.items():
        c, qmax, _ = ARMS[arm] if arm in ARMS else (REF_CAP, "各列で G0 以下", None); nj = YAML_NJ[arm]
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
    # G1x: G1 の半径方向 + x 方向を粗く (ni・x_density_table・ar_max も変える)
    c, qmax, _ = ARMS["G1"]; nj = YAML_NJ["G1"]; mx = mesh_block(nj, c, xcoarse=True); name = f"{base['name']}_cg1x"
    head = (f"# plan tooling-nozzle-core-grid §5.1 #10 の格子 G1x (G1 = c {c}・q_max {qmax}・nj {nj} に加え、x 方向を粗く: 壁際の AR の上限 1 万・目標 9000、"
            f"AR で決まる範囲の x の密度を 1/{XFACTOR:g}、ni {mx['ni']})。\n# core_grid_mesh.py yaml が {PROBLEM.name} から作る (手で編集しない)。"
            f"違いは name・mesh.ni・nj・axis_cap_frac・x_density_table・ar_max だけ。AR ≤ 1 万は 2026-10-10 ユーザ決定 (本 plan の試験に限る)\n")
    body = []; cnt = {"name": 0, "ni": 0, "nj": 0, "xd": 0, "ar": 0}
    for ln in src.splitlines():
        if ln.startswith("name: "):
            ln = f"name: {name}"; cnt["name"] += 1
        elif ln.startswith("  ni: "):
            ln = f"  ni: {mx['ni']}"; cnt["ni"] += 1
        elif ln.startswith("  nj: "):
            body.append(f"  nj: {nj}"); ln = f"  axis_cap_frac: {c}   # 主流と軸付近の間隔 (/局所半径) の上限 (plan tooling-nozzle-core-grid §4.2)"; cnt["nj"] += 1
        elif ln.startswith("  x_density_table: "):
            ln = "  x_density_table: " + json.dumps(mx["x_density_table"]); cnt["xd"] += 1
        elif ln.startswith("  ar_max: "):
            ln = f"  ar_max: {AR_MAX_G1X}"; cnt["ar"] += 1
        body.append(ln)
    if any(v != 1 for v in cnt.values()):
        raise SystemExit(f"G1x: 置き換える行の数が 1 でない {cnt} — 止める")
    text = head + "\n".join(body) + "\n"; new = yaml.safe_load(text)
    b = json.loads(json.dumps(base)); nn = json.loads(json.dumps(new))
    for d in (b, nn):
        d.pop("name")
        for k in ("ni", "nj", "axis_cap_frac", "x_density_table", "ar_max"):
            d["mesh"].pop(k, None)
    if b != nn or new["mesh"]["ni"] != mx["ni"] or new["mesh"]["x_density_table"] != mx["x_density_table"]:
        raise SystemExit("G1x: 生産の YAML との差が想定のキー以外にもある — 止める")
    (HERE / f"{name}.yaml").write_text(text)
    out["G1x"] = {"file": f"{name}.yaml", "ni": mx["ni"], "nj": nj, "axis_cap_frac": c, "ar_max": AR_MAX_G1X, "xfactor": XFACTOR}
    return out


IC_RUN, IC_RES = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5"    # §4.4 の初期値 (B0 と同じ出発点)
G0_MESH_RUN = "run_0183_ns_coldmesh_tw300_ext"                        # 今の格子 G0 の nozzle.h5 (照合の基準)
SAME_RTOL = 1e-12


def _h5(path, *keys):
    import h5py
    with h5py.File(path, "r") as h:
        return [np.array(h[k]) for k in keys]


def cmd_prep(a):
    """§5.1 #3。cold_pair.prep と同じ _prepare_ns (300 K の問題 YAML・同じ壁の δ_r・同じ Euler 参照) で、run ではない準備のディレクトリを作る。"""
    import os
    sys.path.insert(0, str(HERE))
    import cold_pair as CP
    import ns_n012 as NS
    NS.check_dry_env(False)
    conv = Path(os.environ.get("REAL_CONVERTER", ""))
    if not conv.is_file() or NS.sha256_file(conv) != CP.CONV_SHA:
        raise SystemExit(f"REAL_CONVERTER ({conv}) が FP64 の変換器 (sha256 {CP.CONV_SHA[:16]}…) でない — 止める")
    arm = a.arm
    problem = HERE / f"{yaml.safe_load(PROBLEM.read_text())['name']}_{YAML_NAMES.get(arm, 'cg1x')}.yaml"
    out = OUT / f"prep_{arm}"
    if out.exists():
        raise SystemExit(f"{out} が既にある — 止める")
    src = HERE / IC_RUN; rs = NS.res_files(src)
    if not rs or rs[-1].name != IC_RES:
        raise SystemExit(f"IC のドナー {src.name} の最後の res が {rs[-1].name if rs else None} ({IC_RES} であること) — 止める")
    dr_csv = HERE / CP.WALL_REF / "delta_r_initial.csv"
    info = NS._prepare_ns(problem, out, nsteps=12000, ic_from=src, delta_r_csv=str(dr_csv), offset="radial",
                          euler_ref=str(HERE / CP.EULER_REF), cfl_main=5.0, implicit_relax=CP.RELAX)
    NS.jdump(out / "prepare_info.json", info)
    rec = {"plan": "plans/active/tooling-nozzle-core-grid.md §5.1 #3", "tool": "core_grid_mesh.py prep", "created": NS.now(),
           "arm": arm, "problem": problem.name, "problem_sha256": NS.sha256_file(problem),
           "code_commit": (HERE.parents[1] / "COMMIT").read_text().strip() if (HERE.parents[1] / "COMMIT").is_file() else None,
           "forge_bin": os.environ.get("FORGE_BIN"), "forge_sha256": NS.sha256_file(Path(os.environ["FORGE_BIN"])),
           "converter": str(conv), "converter_sha256": CP.CONV_SHA,
           "ic": {"src_run": src.name, "src_res": IC_RES, "src_res_sha256": NS.sha256_file(rs[-1])},
           "delta_r_csv_sha256": NS.sha256_file(dr_csv)}
    rec["mesh_quality"] = CP.mesh_quality_strict(out)
    # 生成時の倍精度座標との照合 (cold_pair.mesh_checks は Mesh2DParams を決め打ちのキーで組むので、ここは _mesh_params_from で組む)
    wall = Wall(out); m = yaml.safe_load(problem.read_text())["mesh"]
    P, prm = generate(wall, m)
    C, dt = _h5(out / "nozzle.h5", "MESH/COORD")[0], None
    import h5py
    with h5py.File(out / "nozzle.h5", "r") as h:
        dt = str(h["MESH/COORD"].dtype); wd = np.array(h["VALUE/wall_dist"])
        V = {k: np.array(h["VALUE/" + k]) for k in ("ro", "Ux", "Uy", "T", "P", "k", "omega") if "VALUE/" + k in h}
    C = C.reshape(prm.ni, prm.nj, 3)[..., :2] / wall.scale
    rel1 = np.abs(np.linalg.norm(C[:, -1] - C[:, -2], axis=1) / np.linalg.norm(P[:, -1] - P[:, -2], axis=1) - 1.0)
    rec["vs_generated"] = {"coord_dtype": dt, "max_abs_diff_rt": float(np.abs(C - P).max()), "first_layer_rel_err_max": float(rel1.max())}
    rec["vs_generated"]["metrics"] = metrics(C)
    # 物理壁 (wall_repr.json) が生産の壁と一致
    rec["geometry_vs_production"] = CP.geom_check(out, wall.scale)
    # 壁距離を変換し直して一致
    wd_new, coord_new, _ = NS.reconvert_wall_dist(out)
    with h5py.File(out / "nozzle.h5", "r") as h:
        coord0 = np.array(h["MESH/COORD"])
    rec["wall_dist_reconvert_rel_max"] = float(np.max(np.abs(wd - wd_new) / np.maximum(np.abs(wd_new), 1e-30))) if np.array_equal(coord0, coord_new) else None
    # G0 (今の格子) との照合: x の station・壁節点・第一内部節点・壁の法線・壁距離 (壁と第一内部節点)
    C0, wd0 = _h5(HERE / G0_MESH_RUN / "nozzle.h5", "MESH/COORD", "VALUE/wall_dist")
    if C0.size // 3 % prm.ni or (C0.size // 3) // prm.ni not in (121,) or prm.ni != 4719:
        # G1x: x の station が G0 と違う → 節点の一致ではなく、両方の格子の第一層の厚さ (壁節点と第一内部節点の距離) が
        # 同じ設計の関数 f(x)·r_w(x) (第一セルの表 × 局所半径) に一致することを見る (G0 の列から補間すると表の折れ点で 0.4 % ずれる)
        C0 = C0.reshape(4719, -1, 3)[..., :2] / wall.scale
        def d1_vs_design(Cg, m_):
            prm_ = _mesh_params_from(m_, wall.scale, int(m_["ni"]), int(m_["nj"]), float(m_["wall_first_frac"]))
            xs_ = mesh2d._x_stations(wall.x_in, wall.x_e, prm_.ni, prm_.throat_refine, prm_.throat_width, density_table=prm_.x_density_table)
            fr_ = mesh2d._first_frac_profile(xs_, prm_)
            d1_ = np.linalg.norm(Cg[:, -1] - Cg[:, -2], axis=1)
            return float(np.max(np.abs(d1_ / (fr_ * wall.r(xs_)) - 1.0)))
        rec["vs_G0"] = {"note": "x の station が違うので節点の一致は見ない",
                        "first_layer_vs_design_rel_max": d1_vs_design(C, m), "G0_first_layer_vs_design_rel_max": d1_vs_design(C0, mesh_block())}
    else:
        rec["vs_G0"] = None
    if rec["vs_G0"] is None:
      C0 = C0.reshape(prm.ni, -1, 3)[..., :2] / wall.scale; nj0 = C0.shape[1]; wd0 = wd0.reshape(prm.ni, nj0); wdn = wd.reshape(prm.ni, prm.nj)
      rel = lambda A, B: float(np.max(np.abs(A - B) / np.maximum(np.abs(B), 1e-30)))  # noqa: E731
      t = np.gradient(C[:, -1], axis=0); t /= np.linalg.norm(t, axis=1)[:, None]; t0 = np.gradient(C0[:, -1], axis=0); t0 /= np.linalg.norm(t0, axis=1)[:, None]
      rec["vs_G0"] = {"wall_nodes_rel_max": rel(C[:, -1], C0[:, -1]), "first_interior_rel_max": rel(C[:, -2], C0[:, -2]),
                      "x_station_rel_max": rel(C[:, -1, 0], C0[:, -1, 0]), "wall_tangent_max_abs_diff": float(np.abs(t - t0).max()),
                      "wall_dist_first_interior_rel_max": rel(wdn[:, -2], wd0[:, -2])}
    # 初期値 (移送後の nozzle.h5 の VALUE): 非有限・物理性
    rec["ic_values"] = {k: {"nonfinite": int(np.count_nonzero(~np.isfinite(v))), "min": float(np.nanmin(v)), "max": float(np.nanmax(v))} for k, v in V.items()}
    ok = {
        "mesh_quality_pass": rec["mesh_quality"].startswith("VERDICT: PASS"),
        "coord_float64": dt == "float64", "first_layer_exact": rec["vs_generated"]["first_layer_rel_err_max"] <= CP.FIRST_LAYER_TOL,
        "matches_generated": rec["vs_generated"]["max_abs_diff_rt"] <= 1e-9,
        "wall_same_as_production": bool(rec["geometry_vs_production"]["ok"]),
        "wall_dist_reconvert": rec["wall_dist_reconvert_rel_max"] is not None and rec["wall_dist_reconvert_rel_max"] <= NS.WALLDIST_RTOL,
        "same_wall_and_first_layer_as_G0": (max(rec["vs_G0"]["first_layer_vs_design_rel_max"], rec["vs_G0"]["G0_first_layer_vs_design_rel_max"]) <= 1e-9) if "note" in rec["vs_G0"] else
                                           (max(rec["vs_G0"]["wall_nodes_rel_max"], rec["vs_G0"]["first_interior_rel_max"], rec["vs_G0"]["x_station_rel_max"]) <= SAME_RTOL
                                            and rec["vs_G0"]["wall_tangent_max_abs_diff"] <= SAME_RTOL),
        "wall_dist_first_interior_as_G0": True if "note" in rec["vs_G0"] else rec["vs_G0"]["wall_dist_first_interior_rel_max"] <= NS.WALLDIST_RTOL,
        "kinks_le_2deg": rec["vs_generated"]["metrics"]["kink_over_2deg"] == 0,
        # CFD ピンの初期線が生産の準備と同じ (run_0062 の nozzle.h5 は AWS で消されたので手元の複製を使う。初期線のハッシュで同一性を確かめる)
        "initial_line_same_as_production": (info.get("initial_line") or {}).get("sha256_16") == json.loads(
            (HERE / CP.WALL_REF / "prepare_info.json").read_text())["initial_line"]["sha256_16"],
        "ic_finite": all(v["nonfinite"] == 0 for v in rec["ic_values"].values()),
        "ic_positive": all(rec["ic_values"][k]["min"] > 0 for k in ("ro", "T", "P") if k in rec["ic_values"]),
    }
    rec["checks"] = ok; rec["VERDICT"] = "PASS" if all(ok.values()) else "FAIL"
    rec["nozzle_sha256"] = NS.sha256_file(out / "nozzle.h5")
    NS.jdump(out / "CORE_GRID_PREP.json", rec)
    return rec


def cmd_prep_compare(a):
    """prep_Gc と prep_G2 の節点座標が軸側から一致する範囲 (列ごとの本数)。上限 c に達した区間が共通なら一致する。"""
    out = {}
    A = _h5(OUT / "prep_Gc" / "nozzle.h5", "MESH/COORD")[0].reshape(4719, -1, 3)[..., :2]
    B = _h5(OUT / "prep_G2" / "nozzle.h5", "MESH/COORD")[0].reshape(4719, -1, 3)[..., :2]
    n = min(A.shape[1], B.shape[1]); same = np.all(np.isclose(A[:, :n], B[:, :n], rtol=0, atol=1e-12), axis=2)
    cnt = np.argmin(np.concatenate([same, np.zeros((same.shape[0], 1), bool)], 1), axis=1)   # 軸から連続して一致する本数
    xw = A[:, -1, 0] / Wall(Path(a.wall_ref)).scale
    for xq in (-8.0, -2.0, 0.0, 5.0, 40.0, 94.0):
        i = int(np.argmin(np.abs(xw - xq))); out[f"x{xq:g}"] = {"common_from_axis": int(cnt[i]), "Gc_nj": A.shape[1], "G2_nj": B.shape[1]}
    out["min_common"] = int(cnt.min()); out["max_common"] = int(cnt.max())
    return out


def cmd_prep_ic(a):
    import h5py
    sys.path.insert(0, str(HERE))
    import cold_xcheck as XC
    yb_x, yb = XC.common_yb(); scale = Wall(OUT / "prep_G1").scale; NI = 4719
    def reduce(xy, ro, ux, uy, k, nj, ni=NI):
        o = XC.reduce_fields(xy, ro, ux, uy, np.ones_like(ro), k, False, ni, nj, scale, yb_x, yb, profile="pchip")
        r = {}
        for (lo, hi) in WINDOWS:
            w = (o["x"] >= lo) & (o["x"] <= hi)
            for q in ("theta_r", "delta_loc"):
                r[f"{q}_w{int(lo)}_{int(hi)}"] = float(np.trapezoid(o[q][w], o["x"][w]) / (o["x"][w][-1] - o["x"][w][0]))
        return r
    out = {}
    with h5py.File(HERE / IC_RUN / IC_RES, "r") as h:
        xy = np.array(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
        g0 = {k: np.array(h["VALUE/" + k], dtype=float) for k in ("ro", "Ux", "Uy", "k")}
        S = {k: np.array(h["VALUE/" + k], dtype=float) for k in ("ro", "roUx", "roUy", "roe", "roK", "roOmega")}
    def prim(D):
        r = D["ro"]
        return {"ro": r, "u": D["roUx"] / r, "v": D["roUy"] / r, "e": D["roe"] / r, "k": D["roK"] / r, "omega": D["roOmega"] / r}
    P0 = prim(S); nj0 = len(S["ro"]) // NI
    def wall_omega(P, nj, ni=NI):
        w = P["omega"].reshape(ni, nj)
        return {"wall_median": float(np.median(w[:, -1])), "wall_max": float(w[:, -1].max()),
                "first_interior_median": float(np.median(w[:, -2])), "max_excl_wall": float(w[:, :-1].max())}
    out["G0_res_100000_omega"] = wall_omega(P0, nj0)
    out["G0_res_100000"] = reduce(xy, g0["ro"], g0["Ux"], g0["Uy"], g0["k"], len(g0["ro"]) // NI)
    for arm in [a_ for a_ in ("G1", "Gc", "G2", "G1x") if (OUT / f"prep_{a_}" / "nozzle.h5").is_file()]:
        with h5py.File(OUT / f"prep_{arm}" / "nozzle.h5", "r") as h:
            keys = [k for k in h["VALUE"].keys() if h["VALUE/" + k].shape and h["VALUE/" + k].ndim == 1]
            V = {k: np.array(h["VALUE/" + k], dtype=float) for k in keys}
            xy = np.array(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
        n = len(V["ro"]); nonf = {k: int(np.count_nonzero(~np.isfinite(v))) for k, v in V.items() if len(v) == n}
        pos = {k: float(V[k].min()) for k in ("ro", "roe", "roK", "roOmega") if k in V}
        # roe は燃焼ガスのエネルギーの基準 (生成エンタルピー) で負になりうる (G0 の res_100000 でも最小 −5.4e6) ので正値は見ない。
        # interp_field は原始量を最近傍で移して ρ を掛け直すので、保存量 (ρ × 別の節点の量) は元の値域を超えうる → 原始量の値域で見る。
        # ω は prepare_ns が移送の後に近壁の下限 6ν/(β₁ d²) を掛ける (runner_axismach.py:1443–1454) ので、壁節点を除いて見る。
        ni_a = int(json.loads((OUT / f"prep_{arm}" / "prepare_info.json").read_text())["mesh"]["ni"])
        P = prim({k: V[k] for k in S}); excess = {}
        for k in P:
            lo, hi = P0[k].min(), P0[k].max(); x = P[k] if k != "omega" else P[k].reshape(ni_a, -1)[:, :-1].ravel()
            excess[k] = float(max(lo - x.min(), x.max() - hi, 0.0) / (hi - lo))
        rec = {"value_keys": sorted(V), "nonfinite": nonf, "min": pos, "primitive_excess_over_range": excess,
               "omega": wall_omega(P, n // ni_a, ni_a),
               "ok": all(v == 0 for v in nonf.values()) and pos.get("ro", 1) > 0 and pos.get("roK", 0) >= 0 and pos.get("roOmega", 1) > 0
                     and all(v <= 1e-6 for v in excess.values())}
        ro = V["ro"]; rec["reduce"] = reduce(xy, ro, V["roUx"] / ro, V["roUy"] / ro, V["roK"] / ro, n // ni_a, ni_a)
        rec["reduce_rel_vs_G0_pct"] = {k: 100 * (rec["reduce"][k] / out["G0_res_100000"][k] - 1) for k in rec["reduce"]}
        out[arm] = rec
    return out


def cmd_view(a):
    """ParaView 用の表示ファイル (mesh_view.py と同じ形式: XDMF + HDF5、四角形の一次要素)。座標は物理壁から生成 (準備の nozzle.h5 と差 0)。
    節点: y_over_rt (列の壁節点までの直線距離 / r_t)、y_plus_est (G0 の解の壁の u_τ/ν_w を同じ x に補間した見積もり)、eta (r/r_w)。
    セル: aspect_ratio、skew、dr_over_dx、dx_over_rt、dr_over_rt。"""
    import csv
    import h5py
    wall = Wall(Path(a.wall_ref)); RT = wall.scale
    rows = list(csv.DictReader(open(Path(a.g0_res).parent / "y1p_115000_wall.csv")))
    xw0 = np.array([float(r["x"]) for r in rows]) / RT
    ut = np.array([np.sqrt(abs(float(r["tau_t"])) / float(r["rho_w"])) for r in rows]); nu = np.array([float(r["mu_w"]) / float(r["rho_w"]) for r in rows])
    o = np.argsort(xw0); xw0, ut_nu = xw0[o], (ut / nu)[o]
    vdir = OUT / "view"; vdir.mkdir(parents=True, exist_ok=True); out = {}
    c1, _, _ = ARMS["G1"]; nj1 = YAML_NJ["G1"]
    for name, m in (("G0", mesh_block()), ("G1", mesh_block(nj1, c1)), ("G1x", mesh_block(nj1, c1, xcoarse=True))):
        P, prm = generate(wall, m); ni, nj = P.shape[:2]
        X = np.concatenate([P.reshape(-1, 2), np.zeros((ni * nj, 1))], 1) * RT
        y = np.linalg.norm(P - P[:, -1:, :], axis=2)                                  # [r_t]
        yplus = y * RT * np.interp(P[:, -1, 0], xw0, ut_nu)[:, None]
        I, J = np.meshgrid(np.arange(ni - 1), np.arange(nj - 1), indexing="ij")
        n00 = (I * nj + J).ravel(); quad = np.stack([n00, n00 + nj, n00 + nj + 1, n00 + 1], 1).astype(np.int32)
        Pq = P.reshape(-1, 2)[quad]
        e = [np.linalg.norm(Pq[:, (k + 1) % 4] - Pq[:, k], axis=1) for k in range(4)]
        E = np.stack(e, 1); ar = E.max(1) / E.min(1)
        ang = []
        for k in range(4):
            u = Pq[:, (k - 1) % 4] - Pq[:, k]; v = Pq[:, (k + 1) % 4] - Pq[:, k]
            ang.append(np.degrees(np.arccos(np.clip((u * v).sum(1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1))))
        A = np.stack(ang, 1); skew = np.maximum((A.max(1) - 90) / 90, (90 - A.min(1)) / 90)
        dx = 0.5 * (e[0] + e[2]); dr = 0.5 * (e[1] + e[3])
        h5 = vdir / f"meshview_{name}.h5"
        with h5py.File(h5, "w") as h:
            h["geom/xyz"] = X; h["topo/quad"] = quad
            h["node/y_over_rt"] = y.ravel(); h["node/y_plus_est"] = yplus.ravel(); h["node/eta"] = (P[..., 1] / P[:, -1:, 1]).ravel()
            h["cell/aspect_ratio"] = ar; h["cell/skew"] = skew; h["cell/dr_over_dx"] = dr / dx; h["cell/dx_over_rt"] = dx; h["cell/dr_over_rt"] = dr
        nn, nc = X.shape[0], quad.shape[0]
        def attr(nm, center, n):
            grp = "node" if center == "Node" else "cell"
            return (f'      <Attribute Name="{nm}" AttributeType="Scalar" Center="{center}">\n'
                    f'        <DataItem Dimensions="{n}" NumberType="Float" Precision="8" Format="HDF">{h5.name}:/{grp}/{nm}</DataItem>\n      </Attribute>\n')
        xmf = ('<?xml version="1.0" ?>\n<Xdmf Version="3.0">\n  <Domain>\n    <Grid Name="case45_' + name + '" GridType="Uniform">\n'
               f'      <Topology TopologyType="Quadrilateral" NumberOfElements="{nc}">\n        <DataItem Dimensions="{nc} 4" NumberType="Int" Precision="4" Format="HDF">{h5.name}:/topo/quad</DataItem>\n      </Topology>\n'
               f'      <Geometry GeometryType="XYZ">\n        <DataItem Dimensions="{nn} 3" NumberType="Float" Precision="8" Format="HDF">{h5.name}:/geom/xyz</DataItem>\n      </Geometry>\n'
               + "".join(attr(k, "Node", nn) for k in ("y_over_rt", "y_plus_est", "eta"))
               + "".join(attr(k, "Cell", nc) for k in ("aspect_ratio", "skew", "dr_over_dx", "dx_over_rt", "dr_over_rt"))
               + "    </Grid>\n  </Domain>\n</Xdmf>\n")
        (vdir / f"meshview_{name}.xmf").write_text(xmf)
        out[name] = {"file": str(vdir / f"meshview_{name}.xmf"), "ni": ni, "nj": nj, "nodes": nn, "cells": nc, "ar_max": float(ar.max()), "skew_max": float(skew.max())}
        print(f"[view] {name}: {vdir / f'meshview_{name}.xmf'} (節点 {nn}、セル {nc})", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["check-g0", "select", "geom-ab", "metric-ab", "bl-count", "post-ab", "yaml", "prep", "prep-compare", "prep-ic", "select-g1x", "view", "select-ref"])
    ap.add_argument("--wall-ref", default=str(WALL_REF)); ap.add_argument("--g0-res", default=str(HERE / "run_0353_m9_L5" / "res_115000.h5"))
    ap.add_argument("arm", nargs="?", choices=["G1", "Gc", "G2", "G1x", "Gref"], help="prep: 格子")
    ap.add_argument("--truth", choices=["pchip", "akima"], default="pchip", help="post-ab: 連続の場の作り方 (G0 の節点の補間)")
    ap.add_argument("--band-dir", default=str(HERE / "_band_ab" / "cold_pair"), help="帯の外縁の npz (theta_run_0181/0183_*.npz) の場所")
    a = ap.parse_args()
    out = {"check-g0": cmd_check_g0, "select": cmd_select, "geom-ab": cmd_geom_ab, "metric-ab": cmd_metric_ab, "bl-count": cmd_bl_count, "post-ab": cmd_post_ab, "yaml": cmd_yaml, "prep": cmd_prep, "prep-compare": cmd_prep_compare, "prep-ic": cmd_prep_ic, "select-g1x": cmd_select_g1x, "view": cmd_view, "select-ref": cmd_select_ref}[a.cmd](a)
    out["problem"] = PROBLEM.name; out["wall_ref"] = a.wall_ref
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / (f"{a.cmd}_{a.truth}.json" if a.cmd == "post-ab" else f"prep_{a.arm}.json" if a.cmd == "prep" else f"{a.cmd}.json"); p.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "metrics"}, ensure_ascii=False, indent=1)[:6000])
    print(f"[core_grid_mesh] {p}")


if __name__ == "__main__":
    main()
