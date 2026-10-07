r"""物理壁の STEP の書き出しと読み直しの検査 (plans/accepted/tooling-nozzle-wall-single-bspline.md §4.3・§6 W5)。

run の壁ファイル (`wall_repr.json`、物理壁の表現 `single_bspline`) の 1 本の 5 次 B-spline $r(x)$ [r_t] を、平面の
B-spline 曲線 $C(u) = (x(u), r(u), 0)$ として STEP (`B_SPLINE_CURVE_WITH_KNOTS`) に 1 本だけ書く。

- 助変数 $u = x$ [mm]。ノットも mm にする。換算 $s = 1000\cdot$`scale_m` (壁ファイルの `units.scale_m` = 問題の
  `spec.r_throat`; ケースの定数を書かない)。$u$ での微分がそのまま $x$ [mm] での微分になる。
- 制御点 $P_i = (\bar x_i, s\,c_i)$。$\bar x_i$ はグレビル点 (ノット $k$ 個の平均)。B-spline は 1 次関数を正確に表すので
  実数演算では $x(u) = u$。次数 5、重み 1 (非有理)。
- 単位は mm。$x$ 軸が流れ方向 (軸)、原点は設計スロート (設計壁の $x = 0$)。物理スロート ($r' = 0$) の位置は添え書きに記す。
  曲線は上半分 ($r \ge 0$)。回転体 (内面) にするのは CAD 側の作業。
- 書き出し・読み直しは FreeCAD (OpenCascade) の `freecadcmd` で行う (`freecad_wall_step_job.py` をファイルで渡す)。
- 添え書き (JSON): 単位・原点・入口端と出口端の座標・物理スロート・CAD の形と CFD の形の関係 (run の `nozzle.h5` の壁節点を
  直線でつないだ多角形と曲線の、全壁辺の符号付きの差「弦 − 曲線」の分布)。

`de_boor` は scipy の BSpline を使わない自前の評価器 (W5 の転送誤差の検査で、読み直した STEP と保存した B-spline を
独立に評価するため)。

usage:
  design/.venv-opt/bin/python -m forge_design.export.wall_step RUN_DIR [--out DIR] [--freecadcmd PATH]
    → DIR/wall_physical.step と DIR/wall_physical_step.json (既定 DIR = RUN_DIR)
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

FREECADCMD = os.environ.get("FREECADCMD", "/home/sano/opt/squashfs-root/usr/bin/freecadcmd")
JOB_SCRIPT = Path(__file__).with_name("freecad_wall_step_job.py")
# 転送誤差の許容 (plan §6 W5): 位置 [mm]・接線方向の角度 [rad]・曲率 (絶対 [1/mm] または相対)
TRANSFER_TOL = {"pos_mm": 1e-6, "angle_rad": 1e-9, "kappa_abs_per_mm": 1e-9, "kappa_rel": 1e-9}


def greville(t, k: int) -> np.ndarray:
    """グレビル点 x̄ᵢ = (t_{i+1} + … + t_{i+k}) / k。"""
    t = np.asarray(t, dtype=float)
    n = len(t) - k - 1
    return np.array([t[i + 1:i + k + 1].sum() / k for i in range(n)])


def step_curve_data(rec: dict) -> dict:
    """壁ファイルの記録 (`load_wall_file(...)["record"]`、表現 single_bspline) → STEP に書く曲線のデータ (mm)。"""
    if rec.get("physical_wall_repr") != "single_bspline":
        raise ValueError(f"STEP は物理壁の表現 single_bspline の壁ファイルだけ (受け取った: {rec.get('physical_wall_repr')!r})")
    ph = rec["physical_wall"]
    k = int(ph["k"])
    s = 1000.0 * float(rec["units"]["scale_m"])           # r_t → mm
    t_mm = np.asarray(ph["t"], dtype=float) * s
    c_mm = np.asarray(ph["c"], dtype=float) * s
    xbar = greville(t_mm, k)
    dist, mult = np.unique(t_mm, return_counts=True)
    return {"degree": k, "scale_mm_per_rt": s, "t_mm": t_mm, "knots_mm": dist, "mults": mult.astype(int),
            "poles_mm": np.c_[xbar, c_mm], "n_poles": int(len(c_mm)),
            "domain_mm": [float(t_mm[0]), float(t_mm[-1])], "min_knot_gap_mm": float(np.diff(dist).min())}


class DeBoor:
    """自前の de Boor 評価器 (scipy を使わない)。t: ノット列、P: (n, d) の制御点、k: 次数。
    `__call__(u, side)` → (nder+1, d): 0..nder 階の u 微分。side = 'right' / 'left' はノット上で右 / 左の区間を使う
    (重複ノットでの右 / 左の極限)。導関数は制御点の差分 (k (Q_{j+1} − Q_j) / (t_{j+k+1} − t_{j+1})) で次数を下げて評価する。"""

    def __init__(self, t, P, k: int, nder: int = 2) -> None:
        tm, Q, p = np.asarray(t, dtype=float), np.asarray(P, dtype=float), int(k)
        if Q.ndim == 1:
            Q = Q[:, None]
        if len(tm) != len(Q) + p + 1:
            raise ValueError("DeBoor: ノット数 ≠ 制御点数 + 次数 + 1")
        self.levels = []
        for m in range(nder + 1):
            self.levels.append((tm.tolist(), Q.tolist(), p))
            if m == nder:
                break
            n = len(Q)
            den = tm[p + 1:p + n] - tm[1:n]
            Q = np.where(den[:, None] > 0.0, p * (Q[1:] - Q[:-1]) / np.where(den > 0.0, den, 1.0)[:, None], 0.0)
            tm, p = tm[1:-1], p - 1
        self.lo, self.hi = float(t[k]), float(t[len(t) - k - 1])

    @staticmethod
    def _span(t, n, p, u, side):
        import bisect
        i = (bisect.bisect_right(t, u) if side == "right" else bisect.bisect_left(t, u)) - 1
        return min(max(i, p), n - 1)

    def __call__(self, u: float, side: str = "right") -> np.ndarray:
        u = float(u)
        out = []
        for t, Q, p in self.levels:
            n = len(Q)
            i = self._span(t, n, p, u, side)
            d = [list(Q[j + i - p]) for j in range(p + 1)]
            for r in range(1, p + 1):
                for j in range(p, r - 1, -1):
                    a0 = t[j + i - p]
                    den = t[j + 1 + i - r] - a0
                    al = (u - a0) / den if den != 0.0 else 0.0
                    d[j] = [(1.0 - al) * x0 + al * x1 for x0, x1 in zip(d[j - 1], d[j])]
            out.append(d[p])
        return np.asarray(out)


def tangent_angle_curvature(d1, d2):
    """平面曲線の接線方向の角度 atan2(y′, x′) と曲率 (x′y″ − y′x″)/(x′² + y′²)^{3/2} (u での微分から)。"""
    d1, d2 = np.asarray(d1, dtype=float), np.asarray(d2, dtype=float)
    ang = np.arctan2(d1[..., 1], d1[..., 0])
    kap = (d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0]) / (d1[..., 0] ** 2 + d1[..., 1] ** 2) ** 1.5
    return ang, kap


def angle_diff(a, b):
    """角度の差 (−π, π] を atan2 で (丸めに弱い acos の内積は使わない)。"""
    return np.arctan2(np.sin(np.asarray(a) - np.asarray(b)), np.cos(np.asarray(a) - np.asarray(b)))


def _freecad(job: dict, freecadcmd: str | None = None, timeout: int = 600) -> dict:
    exe = freecadcmd or FREECADCMD
    if not Path(exe).is_file():
        raise FileNotFoundError(f"freecadcmd が無い: {exe} (FREECADCMD か --freecadcmd で指定)")
    with tempfile.TemporaryDirectory(prefix="wall_step_") as td:
        jp, rp = Path(td) / "job.json", Path(td) / "result.json"
        job = dict(job, result=str(rp))
        jp.write_text(json.dumps(job))
        env = dict(os.environ, WALL_STEP_JOB=str(jp), QT_QPA_PLATFORM="offscreen")
        q = subprocess.run([exe, str(JOB_SCRIPT)], env=env, capture_output=True, text=True, timeout=timeout, cwd=td)
        if not rp.exists():
            raise RuntimeError(f"freecadcmd が結果を返さない (rc={q.returncode}):\n{q.stdout[-2000:]}\n{q.stderr[-2000:]}")
        res = json.loads(rp.read_text())
    if not res.get("ok"):
        raise RuntimeError(f"freecadcmd の作業が失敗:\n{res.get('error')}")
    return res


def write_step(data: dict, step_path, freecadcmd: str | None = None) -> dict:
    """曲線のデータ (`step_curve_data`) を STEP に書く (辺 1 本)。戻り: FreeCAD の報告 (次数・制御点数・非有理)。"""
    job = {"mode": "write", "step": str(Path(step_path).resolve()), "degree": int(data["degree"]),
           "poles": np.asarray(data["poles_mm"], dtype=float).tolist(), "mults": [int(m) for m in data["mults"]],
           "knots": [float(u) for u in data["knots_mm"]]}
    return _freecad(job, freecadcmd)


def read_step(step_path, u=(), revolve: bool = False, freecadcmd: str | None = None) -> dict:
    """STEP を読み直す。u で位置・1 階・2 階微分を評価し、revolve なら x 軸まわりの回転面を作って妥当性と面積を返す。"""
    job = {"mode": "read", "step": str(Path(step_path).resolve()), "u": [float(v) for v in u], "revolve": bool(revolve)}
    return _freecad(job, freecadcmd)


def sample_params(data: dict, n_per_interval: int = 3, eps_rel: float = 1e-12) -> list:
    """転送誤差の検査点: 各ノット区間の内部 n_per_interval 点 + 全ノットの左右の極限 (u = ノット ∓ δ、δ は u の丸めより十分大きい)。
    戻り: [(u, side, label)]。"""
    kn = np.asarray(data["knots_mm"], dtype=float)
    scale = max(abs(kn[0]), abs(kn[-1]))
    dl = max(1e-7, eps_rel * scale * 1e4)
    pts = []
    fr = (np.arange(n_per_interval) + 0.5) / n_per_interval
    for a, b in zip(kn[:-1], kn[1:]):
        for f in fr:
            pts.append((float(a + f * (b - a)), "right", "interior"))
    mult = np.asarray(data["mults"], dtype=int)
    for j, uk in enumerate(kn):
        if j > 0:
            pts.append((float(uk - dl), "left", f"knot_left_mult{mult[j]}"))
        if j < len(kn) - 1:
            pts.append((float(uk + dl), "right", f"knot_right_mult{mult[j]}"))
    return pts


def transfer_check(data: dict, rd: dict, pts: list) -> dict:
    """読み直した STEP (`read_step` の結果) と保存した B-spline (data) の構造と転送誤差 (plan §6 W5)。
    構造: 次数・異なるノット (数と値; 値は STEP の実数の桁数で丸まるので位置の許容差 [mm] 以内)・重複度・制御点数・非有理・辺 1 本・
    平面性 (全制御点の |z| ≤ 位置の許容差)・辺が定義域の全体 (辺の助変数の範囲 = 最初と最後のノット)・実際の端点 (辺の頂点 2 つ) が
    保存した曲線の両端と 3 次元で一致。
    転送誤差: 各点の位置 [mm、z を含む 3 次元]・接線方向の角度 [rad]・曲率 (絶対 [1/mm] または相対)。保存した側は自前の de Boor
    (`DeBoor`) で評価する (2026-10-07 result 段レビュー M1: z の移動・途中で切れた辺を検出していなかった)。"""
    tol = TRANSFER_TOL
    st = {"degree": rd.get("degree") == data["degree"],
          "n_poles": rd.get("n_poles") == data["n_poles"],
          "mults": list(rd.get("mults", [])) == [int(m) for m in data["mults"]],
          "n_knots": len(rd.get("knots", [])) == len(data["knots_mm"]),
          "non_rational": rd.get("is_rational") is False,
          "single_edge": rd.get("n_edges") == 1, "curve_type": rd.get("curve_type") == "BSplineCurve"}
    kn_rd = np.asarray(rd.get("knots", []), dtype=float)
    po_rd = np.asarray(rd.get("poles", []), dtype=float)
    st["knot_values"] = bool(len(kn_rd) == len(data["knots_mm"]) and np.abs(kn_rd - data["knots_mm"]).max() <= tol["pos_mm"])
    st["planar"] = bool(po_rd.ndim == 2 and po_rd.shape[0] > 0 and po_rd.shape[1] == 3 and np.abs(po_rd[:, 2]).max() <= tol["pos_mm"])
    k0, k1 = float(data["knots_mm"][0]), float(data["knots_mm"][-1])
    ef, el = rd.get("edge_first"), rd.get("edge_last")
    st["edge_full_range"] = bool(ef is not None and el is not None and abs(float(ef) - k0) <= tol["pos_mm"] and abs(float(el) - k1) <= tol["pos_mm"])
    db_end = DeBoor(data["t_mm"], data["poles_mm"], data["degree"], nder=0)
    ends = [np.r_[db_end(k0, "right")[0], 0.0], np.r_[db_end(k1, "left")[0], 0.0]]
    vt = np.asarray(rd.get("vertices", []), dtype=float)
    e_end = None
    if vt.shape == (2, 3):
        e_end = max(min(float(np.linalg.norm(v - e)) for v in vt) for e in ends)   # 頂点の順は問わない、両端それぞれに最も近い頂点
    st["edge_endpoints"] = bool(e_end is not None and e_end <= tol["pos_mm"])
    out = {"structure": st,
           "edge_range": [ef, el], "edge_endpoints_max_err_mm": e_end,
           "knots_max_abs_diff_mm": (float(np.abs(kn_rd - data["knots_mm"]).max()) if len(kn_rd) == len(data["knots_mm"]) else None),
           "poles_max_abs_diff_mm": (float(np.abs(po_rd[:, :2] - data["poles_mm"]).max()) if po_rd.shape[0] == data["n_poles"] else None),
           "poles_max_abs_z": (float(np.abs(po_rd[:, 2]).max()) if po_rd.size else None)}
    db = DeBoor(data["t_mm"], data["poles_mm"], data["degree"], nder=2)
    ev = rd.get("eval", [])
    if len(ev) != len(pts):
        out["error"] = f"評価点の数が合わない ({len(ev)} ≠ {len(pts)})"
        out["pass"] = False
        return out
    pos, ang, kap, rows = [], [], [], []
    for (u, side, lab), (p0, d1, d2) in zip(pts, ev):
        ref = db(u, side)
        e_pos = float(np.linalg.norm(np.asarray(p0, dtype=float) - np.r_[ref[0, 0], ref[0, 1], 0.0]))   # z を含む 3 次元
        a_rd, k_rd = tangent_angle_curvature(np.array(d1[:2]), np.array(d2[:2]))
        a_rf, k_rf = tangent_angle_curvature(ref[1], ref[2])
        e_ang = float(abs(angle_diff(a_rd, a_rf)))
        e_kap = float(abs(k_rd - k_rf))
        e_kap_rel = e_kap / abs(float(k_rf)) if k_rf != 0.0 else float("inf")
        pos.append(e_pos); ang.append(e_ang); kap.append((e_kap, e_kap_rel))
        rows.append((u, lab, e_pos, e_ang, e_kap, e_kap_rel, p0[0] - u))
    pos, ang = np.asarray(pos), np.asarray(ang)
    kap_ok = np.array([(a <= tol["kappa_abs_per_mm"]) or (r <= tol["kappa_rel"]) for a, r in kap])
    ip, ia = int(np.argmax(pos)), int(np.argmax(ang))
    ik = int(np.argmax([a for a, _ in kap]))
    out.update({"n_points": len(pts), "tol": dict(tol),
                "pos_max_mm": float(pos.max()), "u_pos_max": rows[ip][0],
                "angle_max_rad": float(ang.max()), "u_angle_max": rows[ia][0],
                "kappa_abs_max_per_mm": float(kap[ik][0]), "kappa_rel_at_abs_max": float(kap[ik][1]), "u_kappa_max": rows[ik][0],
                "kappa_points_failing_both": int((~kap_ok).sum()),
                "x_of_u_max_abs_diff_mm": float(max(abs(r[6]) for r in rows))})
    st_ok = all(st.values())
    out["pass"] = bool(st_ok and pos.max() <= tol["pos_mm"] and ang.max() <= tol["angle_rad"] and kap_ok.all())
    return out


def revolve_ok(rv) -> bool:
    """x 軸まわりの回転面 (内面) が作れたか: 妥当 (`isValid`)・面が 1 枚以上・面積が正 (2026-10-07 result 段レビュー M1:
    CLI の終了判定に含めていなかった)。"""
    return bool(isinstance(rv, dict) and rv.get("is_valid") is True and int(rv.get("n_faces") or 0) >= 1 and float(rv.get("area_mm2") or 0.0) > 0.0)


def chord_deviation(rec: dict, run_dir) -> dict:
    """CFD が解いた多角形 (run の nozzle.h5 の no-slip 壁の節点を x 順に直線でつないだもの) と曲線の半径方向の差
    (弦 − 曲線、符号付き、µm) を全壁辺で測る。正 = 弦が曲線より外 (流路が広い側)。節点座標は float32 [m]。
    端の節点が丸めで定義域の外に出た分は端で評価する。"""
    import h5py
    import yaml
    from ..geometry.wall_axismach import load_wall_file
    run_dir = Path(run_dir)
    h5, bc = run_dir / "nozzle.h5", run_dir / "bcondConfig.yaml"
    if not (h5.exists() and bc.exists()):
        return {"status": "未計算", "reason": f"{h5} または {bc} が無い"}
    W = load_wall_file(run_dir)
    S = float(rec["units"]["scale_m"])
    x_in, x_e = (float(v) for v in rec["domain"])
    walls = [str(v["physID"]) for v in (yaml.safe_load(bc.read_text()) or {}).values()
             if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")]
    with h5py.File(h5) as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3).astype(float)
        ids = np.unique(np.concatenate([f["BCONDS"][w]["vizBfaceNodes"][:].ravel() for w in walls if w in f["BCONDS"]]))
    xy = nc[ids, :2] / S
    xy = xy[np.argsort(xy[:, 0])]
    um = S * 1e6
    rows = []
    for (xa, ra), (xb, rb) in zip(xy[:-1], xy[1:]):
        q = np.linspace(xa, xb, 201)
        chord = ra + (rb - ra) * (q - xa) / (xb - xa)
        e = (chord - W["physical"].r(np.clip(q, x_in, x_e))) * um
        rows.append((xa, xb, float(e.min()), float(q[np.argmin(e)]), float(e.max()), float(q[np.argmax(e)])))
    rows = np.asarray(rows)
    i0, i1 = int(np.argmin(rows[:, 2])), int(np.argmax(rows[:, 4]))
    allv = np.r_[rows[:, 2], rows[:, 4]]
    out = {"status": "ok", "source": f"{h5} の壁節点 (bcond physID {walls}, float32 [m])", "n_wall_nodes": int(len(xy)),
           "n_edges": int(len(rows)), "sign": "弦 − 曲線 (正 = 弦が外 = 流路が広い側)",
           "min_um": float(rows[i0, 2]), "x_of_min_rt": float(rows[i0, 3]), "edge_dx_at_min_rt": float(rows[i0, 1] - rows[i0, 0]),
           "max_um": float(rows[i1, 4]), "x_of_max_rt": float(rows[i1, 5]), "edge_dx_at_max_rt": float(rows[i1, 1] - rows[i1, 0]),
           "percentiles_um": {str(p): float(np.percentile(allv, p)) for p in (1, 5, 50, 95, 99)},
           "node_offset_um": {"max_abs": float(np.abs(xy[:, 1] - W["physical"].r(np.clip(xy[:, 0], x_in, x_e))).max() * um),
                              "note": "壁節点そのもの (float32) と曲線の差"},
           "by_region_um": {}}
    for name, lo, hi in (("pipe+contraction [x_in, 0)", -np.inf, 0.0), ("throat [0, 1)", 0.0, 1.0), ("expansion [1, x_e]", 1.0, np.inf)):
        m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
        if m.any():
            out["by_region_um"][name] = {"min": float(rows[m, 2].min()), "max": float(rows[m, 4].max()),
                                         "edge_dx_max_rt": float((rows[m, 1] - rows[m, 0]).max())}
    return out


def sidecar(rec: dict, data: dict, run_dir=None) -> dict:
    """STEP の添え書き: 単位・原点・入口端と出口端・物理スロート・表現の誤差の保証・CAD と CFD の形の関係。"""
    ph = rec["physical_wall"]
    s = data["scale_mm_per_rt"]
    db = DeBoor(data["t_mm"], data["poles_mm"], data["degree"], nder=0)
    u0, u1 = data["domain_mm"]
    p_in, p_out = db(u0, "right")[0], db(u1, "left")[0]
    out = {"units": "mm", "axis": "x = 流れ方向 (回転軸)、y = 半径 (上半分 r ≥ 0)、z = 0",
           "origin": "設計スロート (設計壁の x = 0)", "scale_mm_per_rt": s, "scale_m": float(rec["units"]["scale_m"]),
           "curve": {"type": "B_SPLINE_CURVE_WITH_KNOTS (非有理、重み 1)", "degree": data["degree"], "n_poles": data["n_poles"],
                     "n_distinct_knots": int(len(data["knots_mm"])), "parameter": "u = x [mm] (ノットも mm)",
                     "min_knot_gap_mm": data["min_knot_gap_mm"], "poles": "Pᵢ = (グレビル点 x̄ᵢ, 係数 cᵢ) [mm]"},
           "inlet_end_mm": [float(p_in[0]), float(p_in[1])], "outlet_end_mm": [float(p_out[0]), float(p_out[1])],
           "physical_throat": {"x_mm": float(ph["throat"]["x"]) * s, "r_mm": float(ph["throat"]["r"]) * s,
                               "note": "物理スロート (r′ = 0)。原点 (設計スロート) からの位置"},
           "joints_mm": [float(v) * s for v in ph["joints"]],
           "representation": {"note": "この曲線は今の物理壁 (区間ごとの式の和) を 1 本の B-spline に作り直したもの。保証は許容誤差まで",
                              "tol_rt": ph.get("tol"), "max_err_vs_legacy_rt": ph.get("max_err_vs_legacy"),
                              "radius_tol_mm": float(ph["tol"]["r"]) * s if ph.get("tol") else None},
           "cad_vs_cfd": {"note": "CAD に渡すのは滑らかな曲線。CFD はその壁の上に置いた節点を直線でつないだ多角形を解いている。"
                                  "加工公差との比較は指定された公差の値で行う (未確認)。受け取り側の CAD での読み込みは未確認"}}
    if run_dir is not None:
        out["cad_vs_cfd"]["chord_minus_curve"] = chord_deviation(rec, run_dir)
    return out


def export_run(run_dir, out_dir=None, freecadcmd: str | None = None, n_per_interval: int = 3) -> dict:
    """run の壁ファイルから STEP と添え書きを書き、読み直して構造と転送誤差を検査する。戻り: 添え書き (検査結果を含む)。"""
    from ..geometry.wall_axismach import load_wall_file
    run_dir = Path(run_dir)
    out_dir = Path(out_dir) if out_dir else run_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    W = load_wall_file(run_dir)
    data = step_curve_data(W["record"])
    step = out_dir / "wall_physical.step"
    wr = write_step(data, step, freecadcmd)
    pts = sample_params(data, n_per_interval)
    rd = read_step(step, [u for u, _, _ in pts], revolve=True, freecadcmd=freecadcmd)
    sc = sidecar(W["record"], data, run_dir)
    sc["step_file"] = step.name
    sc["write"] = wr
    sc["readback"] = {"transfer": transfer_check(data, rd, pts), "revolve": rd.get("revolve"),
                      "revolve_ok": revolve_ok(rd.get("revolve")), "freecad_version": rd.get("freecad_version")}
    (out_dir / "wall_physical_step.json").write_text(json.dumps(sc, indent=1, ensure_ascii=False, default=float))
    return sc


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--freecadcmd", default=None)
    a = ap.parse_args(argv)
    sc = export_run(a.run_dir, a.out, a.freecadcmd)
    tr = sc["readback"]["transfer"]
    rv_ok = revolve_ok(sc["readback"]["revolve"])
    print(json.dumps({"step": sc["step_file"], "transfer_pass": tr["pass"], "pos_max_mm": tr.get("pos_max_mm"),
                      "angle_max_rad": tr.get("angle_max_rad"), "kappa_abs_max_per_mm": tr.get("kappa_abs_max_per_mm"),
                      "revolve": sc["readback"]["revolve"], "revolve_ok": rv_ok}, indent=1, ensure_ascii=False))
    return 0 if (tr["pass"] and rv_ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
