r"""寸法の出荷ゲート (plans/active/tooling-nozzle-wt-deliverables.md §4.2a、codex 諮問 2026-10-13)。

配布する STEP (2D の壁曲線) とメッシュから寸法を**測り直し**、要求カード (`forge.wt_card/v1`) の値と公差で比べる。
YAML の値 × $r_t$ の往復は単位換算の検査であって形状の検査ではないので、ここでは形状そのものを測る。

測り方 (STEP は FreeCAD で読み直した制御点・ノットから scipy の B-spline で評価する。$x$ の原点は設計スロート):

- 入口半径: 曲線の始点の $r$ (端点補間の B-spline なので最初の制御点)。
- 直管長: 始点から、影響する制御点の $r$ がすべて入口半径 (差 ≤ `POLE_EQ_MM`) であるノット区間の最後のノット $u^*$ までの $x$ の長さ。
  半径の閾値では探さない (接合点は $r'=r''=0$ なので閾値は下流にずれる; codex plan レビュー M2)。
- 縮流部長: $x(u^*)$ から設計スロート $x=0$ まで。
- スロート径: 物理 = 曲線の $\min r \times 2$。設計スロート径で指定したカードは設計壁 (`wall_repr.json`) の最小半径 × 2 と比べ、物理は参考値。
- スロート→出口長: **物理出口 = 最終 STEP の下流端** (設計スロートが原点)。メッシュの `outlet` 境界の $x$ と端点の転送・丸めの範囲
  (`END_TOL_MM`) で一致することを別に検査する (codex 諮問 2026-10-13: MOC の $x_F$ で代用しない、途中の弦誤差で許容しない)。

メッシュ側は壁節点から同じ量を測る (直管の終端は接合点を挟む 2 節点の区間)。メッシュと STEP の差は弦の形状誤差
(`wall_step.chord_deviation`) の範囲で見る。

usage:
  design/.venv-opt/bin/python -m forge_design.export.dimension_gate STEP --card CARD [--msh MSH] [--wall-file WALL_REPR] [--out JSON]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

POLE_EQ_MM = 1e-9          # 制御点の半径が入口半径と「等しい」とみなす差 [mm]
NODE_EQ_MM = 1e-6          # メッシュの壁節点が入口半径と等しいとみなす差 [mm] (msh は 10 有効桁)
END_TOL_MM = 1e-3          # 出口端 (STEP の下流端とメッシュの outlet の x) の一致の許容 [mm]: 端点の転送・丸めの検査 (弦誤差では許容しない)
DENSE = 20001              # 最小半径の粗探索の点数


def _bspline(rd: dict):
    from scipy.interpolate import BSpline
    k = int(rd["degree"])
    t = np.repeat(np.asarray(rd["knots"], dtype=float), np.asarray(rd["mults"], dtype=int))
    P = np.asarray(rd["poles"], dtype=float)
    if rd.get("is_rational"):
        raise ValueError("有理 B-spline の STEP は未対応 (wall_step が書くのは非有理)")
    return BSpline(t, P, k, extrapolate=False), t, P, k


def step_dimensions(step_path, freecadcmd: str | None = None) -> dict:
    """STEP (2D の壁曲線 1 本) から寸法を測る [mm]。"""
    from .wall_step import read_step
    rd = read_step(step_path, freecadcmd=freecadcmd)
    if int(rd.get("n_edges", 0)) != 1 or rd.get("curve_type") != "BSplineCurve":
        raise ValueError(f"STEP は B-spline の辺 1 本を想定 (辺 {rd.get('n_edges')}・型 {rd.get('curve_type')})")
    C, t, P, k = _bspline(rd)
    u0, u1 = float(rd["first"]), float(rd["last"])
    if np.abs(P[:, 2]).max() > 1e-9:
        raise ValueError(f"STEP の曲線が z=0 の平面に無い (|z| 最大 {np.abs(P[:, 2]).max():.3g} mm)")
    p0, p1 = C(u0), C(u1)
    r_in = float(P[0, 1])
    # 直管の終端: 区間 m (= [t_m, t_{m+1}]) に効く制御点は P[m-k..m]。始点から入口半径と等しい制御点の個数を j とすると、
    # 区間 m の制御点がすべて等しいのは m ≤ j-1 のとき → 直管の終端のノットは t[j] (j ≥ k+1 のとき; それ未満は直管なし)。
    eq = np.abs(P[:, 1] - r_in) <= POLE_EQ_MM
    j = int(np.argmin(eq)) if not eq.all() else len(P)
    if j >= k + 1:
        u_star = float(t[j])
        x_star = float(C(u_star)[0])
        straight = {"u_end": u_star, "x_end_mm": x_star, "n_equal_poles": j}
    else:
        x_star = float(p0[0])
        straight = {"u_end": None, "x_end_mm": x_star, "n_equal_poles": j, "note": "直管なし (等しい制御点が次数 + 1 未満)"}
    # 最小半径: 密に評価 → 近傍を黄金分割で詰める
    us = np.linspace(u0, u1, DENSE)
    xy = C(us)
    if np.any(np.diff(xy[:, 0]) <= 0):
        raise ValueError("STEP の曲線の x が単調増加でない (折り返し・多価の曲線は寸法を測れない)")
    i = int(np.argmin(xy[:, 1]))
    a, b = us[max(i - 1, 0)], us[min(i + 1, DENSE - 1)]
    g = (np.sqrt(5.0) - 1.0) / 2.0
    for _ in range(200):
        c, d = b - g * (b - a), a + g * (b - a)
        if C(c)[1] < C(d)[1]:
            b = d
        else:
            a = c
        if b - a < 1e-12 * max(1.0, abs(u1 - u0)):
            break
    um = 0.5 * (a + b)
    pm = C(um)
    return {"source": str(step_path), "units": "mm", "x_origin": "設計スロート (x=0)",
            "x_start_mm": float(p0[0]), "x_end_mm": float(p1[0]),
            "inlet_radius_mm": r_in, "straight": straight,
            "L_pipe_mm": float(x_star - p0[0]), "L_U_mm": float(0.0 - x_star),
            "throat_physical": {"r_mm": float(pm[1]), "x_mm": float(pm[0]), "diameter_mm": 2.0 * float(pm[1])},
            "L_div_mm": float(p1[0]), "exit_radius_mm": float(p1[1]),
            "curve": {"degree": k, "n_poles": int(len(P)), "u_range": [u0, u1]}}


def design_throat_diameter_mm(wall_file) -> dict:
    """設計壁 (wall_repr.json の設計壁) の最小半径 × 2 [mm] (設計スロート径で指定したカードの比較用)。"""
    from ..geometry.wall_axismach import load_wall_file
    W = load_wall_file(wall_file)
    s = float(W["scale_m"]) * 1000.0
    x0, x1 = (float(v) for v in W["domain"])
    xs = np.linspace(max(x0, -3.0), min(x1, 3.0), 60001)
    r = np.asarray(W["design"].r(xs), dtype=float)
    i = int(np.argmin(r))
    return {"r_mm": float(r[i] * s), "x_mm": float(xs[i] * s), "diameter_mm": float(2.0 * r[i] * s)}


def mesh_dimensions(msh_path) -> dict:
    """解析メッシュ (mesh2d の msh4.1) の壁節点から寸法を測る [mm]。直管の終端は接合点を挟む 2 節点の区間で返す。"""
    from .mesh_structure import read_msh41_mesh2d
    coords, _quads, bedges, _tags = read_msh41_mesh2d(msh_path)
    e = np.asarray(bedges["wall"], dtype=np.int64)
    # 壁の辺を x の順につないで節点列を作る
    nodes = np.unique(e.ravel())
    xw = coords[nodes, 0] * 1000.0
    rw = coords[nodes, 1] * 1000.0
    o = np.argsort(xw)
    xw, rw = xw[o], rw[o]
    r_in = float(rw[0])
    eq = np.abs(rw - r_in) <= NODE_EQ_MM
    j = int(np.argmin(eq)) if not eq.all() else len(rw)
    if j >= 2:
        interval = [float(xw[j - 1]), float(xw[j])] if j < len(xw) else [float(xw[-1]), float(xw[-1])]
    else:
        interval = [float(xw[0]), float(xw[0])]
    i = int(np.argmin(rw))
    return {"source": str(msh_path), "units": "mm", "n_wall_nodes": int(len(xw)),
            "x_start_mm": float(xw[0]), "x_end_mm": float(xw[-1]), "inlet_radius_mm": r_in,
            "straight_end_interval_mm": interval,
            "L_pipe_interval_mm": [interval[0] - float(xw[0]), interval[1] - float(xw[0])],
            "throat_physical": {"r_mm": float(rw[i]), "x_mm": float(xw[i]), "diameter_mm": 2.0 * float(rw[i])},
            "L_div_mm": float(xw[-1]), "exit_radius_mm": float(rw[-1])}


def card_targets(card: dict) -> dict:
    """要求カードから比べる寸法 [mm] と公差 [mm] を取り出す。"""
    tol = card.get("tolerances") or {}
    if "length_m" not in tol or "diameter_m" not in tol:
        raise ValueError("カードに tolerances.length_m と tolerances.diameter_m が要る")
    th = card["throat"]
    if ("diameter_physical_m" in th) == ("diameter_design_m" in th):
        raise ValueError("カードの throat は diameter_physical_m か diameter_design_m のどちらか 1 つ")
    kind = "physical" if "diameter_physical_m" in th else "design"
    return {"inlet_radius_mm": 1000.0 * float(card["inlet"]["radius_m"]),
            "L_pipe_mm": 1000.0 * float(card["inlet"]["pipe_length_m"]),
            "L_U_mm": 1000.0 * float(card["inlet"]["contraction_length_m"]),
            "L_div_mm": 1000.0 * float(card["divergent_length_m"]),
            "throat_kind": kind,
            "throat_diameter_mm": 1000.0 * float(th[f"diameter_{kind}_m"]),
            "tol_length_mm": 1000.0 * float(tol["length_m"]),
            "tol_diameter_mm": 1000.0 * float(tol["diameter_m"])}


def evaluate_gate(card: dict, step_dims: dict, mesh_dims: dict | None = None, design_throat: dict | None = None,
                  chord: dict | None = None) -> dict:
    """カードと測った寸法を比べる。戻り: 項目ごとの {card, step, mesh, diff, tol, pass} と総合 verdict (PASS / FAIL)。"""
    T = card_targets(card)
    tl, td = T["tol_length_mm"], T["tol_diameter_mm"]
    items = {}

    def add(name, card_v, step_v, tol, mesh_v=None, note=None):
        d = float(step_v - card_v)
        items[name] = {"card_mm": card_v, "step_mm": step_v, "diff_mm": d, "tol_mm": tol, "pass": bool(abs(d) <= tol)}
        if mesh_v is not None:
            items[name]["mesh_mm"] = mesh_v
        if note:
            items[name]["note"] = note

    m = mesh_dims or {}
    add("inlet_radius", T["inlet_radius_mm"], step_dims["inlet_radius_mm"], td / 2.0, m.get("inlet_radius_mm"))
    add("L_pipe", T["L_pipe_mm"], step_dims["L_pipe_mm"], tl, m.get("L_pipe_interval_mm"))
    add("L_U", T["L_U_mm"], step_dims["L_U_mm"], tl)
    add("L_div", T["L_div_mm"], step_dims["L_div_mm"], tl, m.get("L_div_mm"))
    phys = step_dims["throat_physical"]["diameter_mm"]
    if T["throat_kind"] == "physical":
        add("throat_diameter_physical", T["throat_diameter_mm"], phys, td,
            (m.get("throat_physical") or {}).get("diameter_mm"))
    else:
        if design_throat is None:
            raise ValueError("設計スロート径のカードには設計壁の最小径 (design_throat_diameter_mm) が要る")
        add("throat_diameter_design", T["throat_diameter_mm"], design_throat["diameter_mm"], td,
            note=f"物理スロート径 (参考) {phys:.6f} mm")
    # メッシュと STEP の整合 (弦の形状誤差の範囲か): 出口位置・入口半径は節点が曲線上なので丸めの範囲、最小半径は弦で取り逃しうる
    consistency = {}
    if mesh_dims:
        consistency["x_end_diff_mm"] = float(m["x_end_mm"] - step_dims["x_end_mm"])
        consistency["x_end_pass"] = bool(abs(consistency["x_end_diff_mm"]) <= END_TOL_MM)
        consistency["inlet_radius_diff_mm"] = float(m["inlet_radius_mm"] - step_dims["inlet_radius_mm"])
        consistency["throat_r_mesh_minus_curve_mm"] = float(m["throat_physical"]["r_mm"] - step_dims["throat_physical"]["r_mm"])
        lo, hi = m["straight_end_interval_mm"]
        consistency["straight_end_in_mesh_interval"] = bool(lo - 1e-6 <= step_dims["straight"]["x_end_mm"] <= hi + 1e-6)
    if chord:
        consistency["chord"] = chord
    ok = all(v["pass"] for v in items.values()) and consistency.get("x_end_pass", True)
    return {"format": "forge.dimension_gate/v1", "verdict": "PASS" if ok else "FAIL", "items": items,
            "throat_kind": T["throat_kind"], "mesh_vs_step": consistency,
            "note": "STEP とメッシュから測り直した寸法をカードと比べた (YAML の往復ではない)。出口マッハ数は CFD の判定で、ここでは測らない"}


def main(argv=None) -> int:
    import yaml
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step")
    ap.add_argument("--card", required=True)
    ap.add_argument("--msh", default=None)
    ap.add_argument("--wall-file", default=None, help="設計スロート径のカードのときの wall_repr.json (または run ディレクトリ)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--freecadcmd", default=None)
    a = ap.parse_args(argv)
    card = yaml.safe_load(Path(a.card).read_text())
    sd = step_dimensions(a.step, a.freecadcmd)
    md = mesh_dimensions(a.msh) if a.msh else None
    dt = design_throat_diameter_mm(a.wall_file) if a.wall_file else None
    res = evaluate_gate(card, sd, md, dt)
    res["step_dimensions"], res["mesh_dimensions"], res["design_throat"] = sd, md, dt
    txt = json.dumps(res, indent=1, ensure_ascii=False)
    if a.out:
        Path(a.out).write_text(txt)
    print(txt)
    print(f"VERDICT: {res['verdict']}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
