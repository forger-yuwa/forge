r"""設計壁 (Euler で回す壁) と物理壁 (実際に作る壁) の座標の比較の表と図 (plan tooling-nozzle-card-local-mesh §5.1 #3u、
ユーザ 2026-10-15「オイラーで作った壁と、物理壁の座標をくらべたデータとプロットも、最終的な成果物としてほしい」)。

- 設計壁: `wall_repr.json` の `design_wall` (直管 + 上流 Hermite + MOC の壁 S。境界層を無視した形で、Euler はこの壁で回す)。
- 物理壁: `wall_repr.json` の `physical_wall` (設計壁 + 境界層の厚み δ_r。STEP と同じ曲線)。
どちらも `geometry.wall_axismach.load_wall_file(run)` の評価関数 (r_t 単位) を `scale_m` で mm に換算して評価する (外挿しない)。

- `compare_table`: 表 `wall_design_vs_physical_mm.csv` (列 `x_mm, r_design_mm, r_physical_mm, dr_mm, angle_design_deg,
  angle_physical_deg`)。x の点は形状表 `wall_shape_mm.csv` と同じ (`wall_tables.table_x_mm`)。片方の壁の定義域の外の x は空欄。
- `delta_record_check`: 表の Δr のスロート (x = 0) と出口が run の `prepare_info.json` の δ_r の記録と一致するか
  (許容 `DR_RECORD_TOL_MM`; 記録をその作り方で再現して比べる。再現できない点は照合しない = 合格にしない)。
- `plot_compare`: 図 3 枚 (全長の壁の形・Δr の分布・スロート付近の拡大)。
- `export_compare`: 上の 3 つを OUT_DIR に書き、説明の JSON `wall_design_vs_physical.json` (列の意味・区間・照合) も書く。

usage (design/.venv-opt):
  python -m forge_design.export.wall_compare RUN_DIR OUT_DIR
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

from .wall_tables import _write_csv, table_x_mm

CSV_NAME = "wall_design_vs_physical_mm.csv"
JSON_NAME = "wall_design_vs_physical.json"
PNG_FULL = "wall_design_vs_physical_full.png"
PNG_DR = "wall_design_vs_physical_dr.png"
PNG_THROAT = "wall_design_vs_physical_throat.png"
COLUMNS = ("x_mm", "r_design_mm", "r_physical_mm", "dr_mm", "angle_design_deg", "angle_physical_deg")
DR_RECORD_TOL_MM = 1e-6        # 表の Δr と prepare の δ_r の記録 (同じ作り方で再現した値) の一致の許容 [mm]
ZOOM_MM = (-300.0, 500.0)      # スロート付近の拡大の x の範囲 [mm] (定義域で切る)
JP_FONT = os.environ.get("FORGE_JP_FONT", os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf"))


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _load(wall) -> dict:
    """`load_wall_file` の戻り (dict) はそのまま、run ディレクトリ・壁ファイルの道は読む。"""
    if isinstance(wall, dict) and "design" in wall and "physical" in wall:
        return wall
    from ..geometry.wall_axismach import load_wall_file
    return load_wall_file(wall)


def _segments(W: dict, s_mm: float) -> dict:
    """設計壁・物理壁の区間 [mm] と、区間ごとの Δr の意味 (注記の行)。"""
    d, p = W["design"], W["physical"]
    rec = W["record"]["physical_wall"]
    xu_d = -float(d.up.L_U)
    mm = lambda v: float(v) * s_mm                                          # noqa: E731
    des = [{"part": "直管 (半径 r_U 一定)", "x_mm": [mm(d.x_in), mm(xu_d)]},
           {"part": "縮流部 (上流 Hermite)", "x_mm": [mm(xu_d), 0.0]},
           {"part": "拡大部 (MOC の壁 S)", "x_mm": [0.0, mm(d.x_e)]}]
    if W["pw_upstream"] == "poly":
        xu_p = float(rec["joints"][0]) if W["kind"] == "single_bspline" else -float(p.L_U)
        phy = [{"part": "直管 (半径 r_U 一定)", "x_mm": [mm(p.x_in), mm(xu_p)]},
               {"part": "縮流部 (多項式 Q: 上流端で直管、x = 0 で S + δ_r に値・1 階・2 階の微分が連続)", "x_mm": [mm(xu_p), 0.0]},
               {"part": "拡大部 (設計壁 S + δ_r)", "x_mm": [0.0, mm(p.x_e)]}]
        x_in = float(max(d.x_in, p.x_in))
        dr_pipe = float((p.r(np.array([x_in])) - d.r(np.array([x_in])))[0]) * s_mm
        meaning = ["x ≥ 0 (設計スロートから出口): Δr = δ_r = 境界層の補正 (物理壁 = 設計壁 S + δ_r)。冷却壁ではスロート付近で負になりうる",
                   "縮流部 (x < 0): 物理壁は多項式 Q、設計壁は Hermite で別の曲線。Δr は境界層の厚みではなく 2 つの上流の形の違い"
                   " (x = 0 で δ_r(0) に連続)",
                   f"直管: Δr = 物理壁の r_U − 設計壁の r_U (この壁では {dr_pipe!r} mm)"]
        upstream_is_bl = False
    else:
        lo, hi = (float(v) for v in rec["ramp"])
        phy = [{"part": f"設計壁 + s(x)·δ_r (s は 5 次 smoothstep、x ∈ [{mm(lo)!r}, {mm(hi)!r}] mm で 0 → 1)",
                "x_mm": [mm(p.x_in), mm(p.x_e)]}]
        meaning = [f"Δr = s(x)·δ_r: x ≥ {mm(hi)!r} mm は境界層の補正 δ_r、ランプの区間は δ_r を 0 から入れる途中、"
                   f"x ≤ {mm(lo)!r} mm は 0"]
        upstream_is_bl = True
    repr_note = None
    if W["kind"] == "single_bspline":
        e0 = (rec.get("max_err_vs_legacy") or {}).get("d0")
        repr_note = ("物理壁は 5 次 B-spline 1 本 (STEP と同じ曲線)。区分の式 (直管・Q・S + δ_r) との差は値で最大 "
                     + (f"{float(e0) * s_mm:.3g} mm" if e0 is not None else "不明 (記録なし)"))
    return {"design": des, "physical": phy, "meaning": meaning, "repr_note": repr_note, "upstream_is_bl": upstream_is_bl}


def compare_table(wall, out_csv, coarse_mm: float = 1.0, fine_per_mm: int = 10, throat_half_width_mm: float | None = None,
                  run_name: str | None = None) -> dict:
    """設計壁と物理壁の比較の表 (plan §5.1 #3u)。x の点 = 2 つの壁の定義域の和の範囲で形状表と同じ刻み (`table_x_mm`)。
    各壁は自分の定義域 [x_in, x_e] の中だけ評価し (外は空欄)、r [mm]・壁角 atan(r′) [deg] を出す。dr_mm = r_physical_mm − r_design_mm
    (片方が空欄なら空欄)。戻り: 説明 (列の意味・単位・原点・刻み・区間・注記) の dict と、作図・照合用の配列 (キー `_x` ほか)。"""
    from ..geometry.wall_axismach import DESIGN_UPSTREAM_KIND
    W = _load(wall)
    if W["kind"] == DESIGN_UPSTREAM_KIND:
        raise ValueError("壁ファイルの physical が H の上流を付けた設計壁 (design_upstream = Euler の参照の壁) で物理壁でない — 比較しない")
    s_mm = 1000.0 * float(W["scale_m"])
    walls = {"design": W["design"], "physical": W["physical"]}
    dom = {k: (float(w.x_in), float(w.x_e)) for k, w in walls.items()}
    lo, hi = min(v[0] for v in dom.values()), max(v[1] for v in dom.values())
    xm, hw = table_x_mm(lo * s_mm, hi * s_mm, s_mm, coarse_mm, fine_per_mm, throat_half_width_mm)
    ev = {}
    for k, w in walls.items():
        a, b = dom[k]
        m = (xm >= a * s_mm) & (xm <= b * s_mm)
        xr = np.clip(xm[m] / s_mm, a, b)                     # mm → r_t の往復の丸めで定義域から出た分は端で評価 (形状表と同じ)
        r, ang = np.full(len(xm), np.nan), np.full(len(xm), np.nan)
        r[m] = w.r(xr) * s_mm
        ang[m] = np.degrees(np.arctan(w.r(xr, 1)))
        ev[k] = (r, ang, m)
    dr = ev["physical"][0] - ev["design"][0]
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_csv, COLUMNS, [xm, ev["design"][0], ev["physical"][0], dr, ev["design"][1], ev["physical"][1]])
    seg = _segments(W, s_mm)
    dom_mm = {k: [v[0] * s_mm, v[1] * s_mm] for k, v in dom.items()}
    n_blank = {k: int((~ev[k][2]).sum()) for k in walls}
    th = W["throat"]
    i0 = np.flatnonzero(xm == 0.0)
    r_d0 = float(ev["design"][0][i0[0]]) if len(i0) else float("nan")
    fin = np.isfinite(dr)
    imin, imax = (int(np.nanargmin(dr)), int(np.nanargmax(dr))) if fin.any() else (None, None)
    fmt = lambda a: "[" + ", ".join(f"{v:.6g}" for v in a) + "]"                     # noqa: E731
    notes = [f"{('run ' + run_name + ': ') if run_name else ''}原点 = 設計スロート (設計壁の x = 0)、x = 流れ方向 (回転軸)、r = 半径。"
             f"単位 mm (角度は deg)。x の点は形状表 wall_shape_mm.csv と同じ ({coarse_mm:g} mm ごと + |x| ≤ {hw:.6g} mm は "
             f"{1.0 / fine_per_mm:g} mm ごと)",
             "r_design_mm: 設計壁 = Euler で回す壁 (境界層を無視した MOC の形; wall_repr.json の design_wall)",
             "r_physical_mm: 物理壁 = 実際に作る壁 (設計壁 + 境界層の厚み; wall_repr.json の physical_wall = STEP と同じ曲線)",
             "dr_mm = r_physical_mm − r_design_mm。angle_design_deg・angle_physical_deg = atan(dr/dx) (正 = 下流で広がる)",
             "設計壁の区間: " + "・".join(f"{s_['part']} {fmt(s_['x_mm'])} mm" for s_ in seg["design"]),
             "物理壁の区間: " + "・".join(f"{s_['part']} {fmt(s_['x_mm'])} mm" for s_ in seg["physical"])]
    notes += ["Δr の意味 — " + t for t in seg["meaning"]]
    if seg["repr_note"]:
        notes.append(seg["repr_note"])
    notes.append(f"定義域: 設計壁 {fmt(dom_mm['design'])} mm・物理壁 {fmt(dom_mm['physical'])} mm。片方の定義域の外の x はその壁の列と "
                 f"dr_mm を空欄にした (空欄: 設計壁 {n_blank['design']} 行・物理壁 {n_blank['physical']} 行)")
    notes.append(f"物理スロート (物理壁の最小径): x = {float(th['x']) * s_mm:.6g} mm・r = {float(th['r']) * s_mm:.9g} mm "
                 f"(設計スロートの r = {r_d0:.9g} mm)")
    return {"csv": str(out_csv), "n_rows": int(len(xm)), "columns": list(COLUMNS), "source": str(W["path"]),
            "physical_wall_repr": W["kind"], "pw_upstream": W["pw_upstream"], "scale_mm_per_rt": s_mm,
            "units": "mm (角度は deg)", "origin": "設計スロート (設計壁の x = 0)。x = 流れ方向 (回転軸)、r = 半径",
            "spacing": {"coarse_mm": coarse_mm, "fine_mm": 1.0 / fine_per_mm, "fine_zone_mm": [max(-hw, lo * s_mm), min(hw, hi * s_mm)],
                        "fine_zone_rule": "|x| ≤ max(50 mm, 2 r_t) (形状表と同じ)"},
            "domain_mm": dom_mm, "n_blank": n_blank, "segments": {"design": seg["design"], "physical": seg["physical"]},
            "dr_meaning": seg["meaning"], "upstream_dr_is_boundary_layer": seg["upstream_is_bl"],
            "physical_throat_mm": {"x": float(th["x"]) * s_mm, "r": float(th["r"]) * s_mm},
            "dr_range_mm": ({"min": float(dr[imin]), "x_at_min": float(xm[imin]), "max": float(dr[imax]), "x_at_max": float(xm[imax])}
                            if imin is not None else None),
            "column_notes": {
                "x_mm": "軸方向の位置 [mm]。原点 = 設計スロート",
                "r_design_mm": "設計壁 (Euler で回す壁、境界層なし) の半径 [mm]",
                "r_physical_mm": "物理壁 (設計壁 + 境界層の厚み、実際に作る壁) の半径 [mm]",
                "dr_mm": "r_physical_mm − r_design_mm [mm] (区間で意味が違う: notes)",
                "angle_design_deg": "設計壁の傾き atan(dr/dx) [deg]",
                "angle_physical_deg": "物理壁の傾き atan(dr/dx) [deg] (形状表の wall_angle_deg と同じ)"},
            "notes": notes,
            "_x": xm, "_r_design": ev["design"][0], "_r_physical": ev["physical"][0], "_dr": dr,
            "_x_design_end_mm": dom_mm["design"][1], "_x_design_end_rt": dom["design"][1]}


# ---------------------------------------------------------------------------------------------------------------------
# prepare の δ_r の記録との照合
# ---------------------------------------------------------------------------------------------------------------------

def _row(tab: dict, x_mm: float):
    i = np.flatnonzero(tab["_x"] == x_mm)
    return float(tab["_dr"][i[0]]) if len(i) and np.isfinite(tab["_dr"][i[0]]) else None


def _find_table(run: Path, info: dict):
    """delta_r_csv の表を探す: 記録の道 → case ディレクトリ (run の親) → run の中、の順で同じ名前のもの。sha256 の記録と照合する。
    戻り (道, 説明)。見つからない・sha256 が合わないときは (None, 理由)。"""
    rec_path = str(info.get("delta_r_csv"))
    tbl = info.get("delta_r_table") or {}
    want = tbl.get("sha256")
    name = Path(rec_path).name
    tried = []
    for c in (Path(rec_path), run.parent / name, run / name):
        if c.is_file():
            h = _sha(c)
            tried.append(f"{c} (sha256 {h[:12]}…)")
            if want is None:
                return c, f"{c} (sha256 の記録が無いので照合していない)"
            if h == want:
                return c, f"{c} (sha256 が prepare_info.json の delta_r_table.sha256 と一致)"
    return None, (f"delta_r_csv の表 {name} が見つからないか sha256 が記録 ({(want or '記録なし')[:12]}…) と合わない "
                  f"(探した: 記録の道・case ディレクトリ・run の中; 見つかったもの {tried or 'なし'})")


def delta_record_check(run, tab: dict, W: dict, tol_mm: float = DR_RECORD_TOL_MM) -> dict:
    """表の Δr のスロート (x = 0) と出口 (設計壁の終点) を run の `prepare_info.json` の δ_r の記録と比べる (plan §5.1 #3u)。
    記録は作り方ごとに再現して比べる (再現できない点は照合しない = `pass` None):

    - `initializer` (積分法の初期壁、`delta_r_throat`・`delta_r_exit` は r_t 単位): 出口の記録は積分の表の最後の点の値で、壁はその値を通る
      (表より下流は端値) ので直接比べる。スロートの記録は `np.interp(0, 表)` (直線補間; runner_axismach の integral_delta_r) で、壁は
      同じ表の点を通る 5 次補間 (`runner_axismach.delta_r_from_table`) なので、点の間の x = 0 では直線補間の誤差だけ違う。run に表
      (`delta_r_initial.csv`) があれば、記録が表の直線補間と一致することを確かめ、表を壁と同じ 5 次補間で x = 0 に評価した値と比べる。
      無ければ照合しない (差だけ記録)。
    - `delta_r_csv` (δ_r の表を渡した run): 表を探して sha256 を記録と照合し、prepare と同じ作り方 (表の端を `delta_r_table_ends` の記録
      どおり揃え、5 次補間) で x = 0 と出口に評価した値と比べる。
    戻り: 点ごとの値・差・合否と、全体の `status` (PASS / FAIL / 一部未照合 / 未照合)・`pass` (True / False / None)。"""
    run = Path(run)
    s_mm = 1000.0 * float(W["scale_m"])
    x_e = float(tab["_x_design_end_rt"])
    pts = {"throat": {"x_mm": 0.0, "dr_table_mm": _row(tab, 0.0)},
           "exit": {"x_mm": float(tab["_x_design_end_mm"]), "dr_table_mm": _row(tab, float(tab["_x_design_end_mm"]))}}
    out = {"tol_mm": tol_mm, "points": pts}

    def finish(status_if_none="未照合"):
        ps = [v.get("pass") for v in pts.values()]
        if any(v is False for v in ps):
            out["status"], out["pass"] = "FAIL", False
        elif all(v is True for v in ps):
            out["status"], out["pass"] = "PASS", True
        else:
            out["status"], out["pass"] = ("一部未照合" if any(v is True for v in ps) else status_if_none), None
        return out

    def judge(key, ref_rt, how, record_rt=None):
        pt = pts[key]
        pt["how"] = how
        if record_rt is not None:
            pt["record_mm"] = float(record_rt) * s_mm
            if pt["dr_table_mm"] is not None:
                pt["diff_vs_record_mm"] = pt["dr_table_mm"] - pt["record_mm"]
        if ref_rt is None:
            pt["pass"] = None
            return
        pt["reference_mm"] = float(ref_rt) * s_mm
        if pt["dr_table_mm"] is None:
            pt["pass"], pt["reason"] = None, "表にこの x の行が無いか空欄"
            return
        pt["diff_mm"] = pt["dr_table_mm"] - pt["reference_mm"]
        pt["pass"] = bool(abs(pt["diff_mm"]) <= tol_mm)

    p = run / "prepare_info.json"
    if not p.exists():
        out["reason"] = "prepare_info.json が無い"
        return finish()
    info = json.loads(p.read_text())
    out["prepare_info_sha256"] = _sha(p)
    sc = info.get("scale_m")
    if sc is not None and float(sc) != float(W["scale_m"]):
        out["reason"] = f"prepare_info.json の scale_m {sc!r} が壁ファイルの scale_m {W['scale_m']!r} と違う (別の run の記録)"
        for v in pts.values():
            v["pass"] = False
        return finish()
    ini = info.get("initializer")
    from ..evaluate.runner_axismach import delta_r_from_table
    if isinstance(ini, dict) and ini.get("delta_r_throat") is not None and ini.get("delta_r_exit") is not None:
        out["source"] = "initializer"
        rt, re_ = float(ini["delta_r_throat"]), float(ini["delta_r_exit"])
        judge("exit", re_, "記録 initializer.delta_r_exit (積分の表の最後の点の値。壁は表の点を通り、表より下流は端値) と直接比べる", re_)
        csv = run / "delta_r_initial.csv"
        if csv.is_file():
            t = np.loadtxt(csv, delimiter=",", skiprows=1, ndmin=2)
            xt, dt = t[:, 0], t[:, 1]
            lin = float(np.interp(0.0, xt, dt))
            pts["throat"]["record_equals_table_linear"] = bool(abs(lin - rt) * s_mm <= 1e-12)
            pts["throat"]["table"] = {"file": csv.name, "sha256": _sha(csv), "n": int(len(xt))}
            judge("throat", float(delta_r_from_table(xt, dt)(np.array([0.0]))[0]),
                  "記録 initializer.delta_r_throat は積分の表の直線補間 np.interp(0, 表)。run の delta_r_initial.csv を壁と同じ 5 次補間 "
                  "(runner_axismach.delta_r_from_table) で x = 0 に評価した値と比べる", rt)
            if not pts["throat"]["record_equals_table_linear"]:
                pts["throat"]["pass"] = False
                pts["throat"]["reason"] = "delta_r_initial.csv の直線補間が記録 delta_r_throat と合わない (記録を作った表でない)"
        else:
            judge("throat", None, "記録 initializer.delta_r_throat は積分の表 (1500 点) の直線補間 np.interp(0, 表) で、壁は同じ表の点を通る "
                  "5 次補間。表 (delta_r_initial.csv) が run に無いので同じ作り方で再現できない → 照合しない (差は diff_vs_record_mm に記録)", rt)
            pts["throat"]["reason"] = "delta_r_initial.csv が run に無い (直線補間の記録を再現できない)"
        return finish()
    if info.get("delta_r_csv"):
        out["source"] = "delta_r_csv"
        path, why = _find_table(run, info)
        out["table"] = why
        if path is None:
            for v in pts.values():
                v["pass"] = None
                v["reason"] = why
            return finish()
        t = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        xt, dt = t[:, 0].copy(), t[:, 1]
        ends = info.get("delta_r_table_ends") or {}
        for key, idx in (("start", 0), ("end", -1)):
            e = ends.get(key) or {}
            if e.get("snapped"):
                xt[idx] = float(e["x_design_rt"])
        f = delta_r_from_table(xt, dt)
        how = ("δ_r の表 (delta_r_csv) を prepare と同じ作り方 (表の端を delta_r_table_ends の記録どおり揃え、5 次補間 "
               "runner_axismach.delta_r_from_table) で評価した値と比べる")
        judge("throat", float(f(np.array([0.0]))[0]), how)
        judge("exit", float(f(np.array([x_e]))[0]), how)
        return finish()
    out["reason"] = "prepare_info.json に δ_r の記録 (initializer の delta_r_throat・delta_r_exit / delta_r_csv) が無い"
    return finish()


# ---------------------------------------------------------------------------------------------------------------------
# 図
# ---------------------------------------------------------------------------------------------------------------------

def _font_rc() -> tuple:
    """日本語フォント (Noto Sans CJK JP) を登録した rc の辞書と、使ったフォントの説明。無ければ既定のフォント (日本語は豆腐になる)。"""
    from matplotlib import font_manager
    rc = {"axes.unicode_minus": False}
    if Path(JP_FONT).is_file():
        if not any(getattr(f, "fname", None) == JP_FONT for f in font_manager.fontManager.ttflist):
            font_manager.fontManager.addfont(JP_FONT)
        name = font_manager.FontProperties(fname=JP_FONT).get_name()
        rc["font.family"] = [name, "DejaVu Sans"]
        return rc, f"{name} ({JP_FONT})"
    return rc, f"既定 (日本語フォント {JP_FONT} が無い)"


def _legend_below(ax, handles=None, ncol=2):
    """凡例は軸の下の外に置く (データに重ねない)。"""
    kw = {"loc": "upper center", "bbox_to_anchor": (0.5, -0.16), "ncol": ncol, "frameon": False, "fontsize": 9}
    return ax.legend(handles=handles, **kw) if handles is not None else ax.legend(**kw)


def plot_compare(tab: dict, out_dir, title: str, zoom_mm=ZOOM_MM, dpi: int = 150) -> dict:
    """図 3 枚 (PNG): (a) 全長の設計壁と物理壁、(b) Δr(x)、(c) スロート付近 (`zoom_mm`、定義域で切る) の 2 つの壁 (左軸) と Δr (右軸)。
    題に `title` (run の名前)。凡例は軸の下の外。"""
    import matplotlib
    from matplotlib.figure import Figure
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    x, rd, rp, dr = tab["_x"], tab["_r_design"], tab["_r_physical"], tab["_dr"]
    th = tab["physical_throat_mm"]
    lab_d = "設計壁 (Euler で回す壁、境界層なし)"
    lab_p = "物理壁 (設計壁 + 境界層の厚み = 作る壁)"
    lab_0 = "設計スロート x = 0"
    xl = "x [mm] (原点 = 設計スロート)"
    rc, font = _font_rc()
    files = {}
    with matplotlib.rc_context(rc):
        # (a) 全長
        fig = Figure(figsize=(11.0, 4.8))
        ax = fig.add_subplot()
        ax.plot(x, rd, "--", color="C0", lw=1.4, label=lab_d)
        ax.plot(x, rp, "-", color="C3", lw=1.0, label=lab_p)
        ax.axvline(0.0, color="0.5", ls=":", lw=0.9, label=lab_0)
        ax.set_xlabel(xl)
        ax.set_ylabel("半径 r [mm]")
        ax.set_title(f"{title}: 設計壁と物理壁 (全長)")
        ax.grid(True, lw=0.3, alpha=0.6)
        _legend_below(ax, ncol=3)
        fig.savefig(out_dir / PNG_FULL, dpi=dpi, bbox_inches="tight")
        files["full"] = PNG_FULL
        # (b) Δr
        fig = Figure(figsize=(11.0, 4.8))
        ax = fig.add_subplot()
        if not tab["upstream_dr_is_boundary_layer"] and np.nanmin(x) < 0.0:
            ax.axvspan(float(np.nanmin(x)), 0.0, color="0.92", lw=0,
                       label="x < 0 (直管・縮流部): 物理壁の上流は設計壁と別の曲線 (Δr は境界層の厚みではない)")
        ax.axhline(0.0, color="0.3", lw=0.7)
        ax.plot(x, dr, "-", color="C2", lw=1.2, label="Δr = r_物理 − r_設計")
        ax.axvline(0.0, color="0.5", ls=":", lw=0.9, label=lab_0)
        ax.set_xlabel(xl)
        ax.set_ylabel("Δr = r_物理 − r_設計 [mm]")
        ax.set_title(f"{title}: 物理壁と設計壁の半径の差 Δr")
        ax.grid(True, lw=0.3, alpha=0.6)
        _legend_below(ax, ncol=3)
        fig.savefig(out_dir / PNG_DR, dpi=dpi, bbox_inches="tight")
        files["dr"] = PNG_DR
        # (c) スロート付近の拡大
        z0, z1 = max(float(zoom_mm[0]), float(np.nanmin(x))), min(float(zoom_mm[1]), float(np.nanmax(x)))
        m = (x >= z0) & (x <= z1)
        fig = Figure(figsize=(11.0, 5.2))
        ax = fig.add_subplot()
        h1 = ax.plot(x[m], rd[m], "--", color="C0", lw=1.4, label=lab_d)[0]
        h2 = ax.plot(x[m], rp[m], "-", color="C3", lw=1.0, label=lab_p)[0]
        hs = [h1, h2]
        if z0 <= th["x"] <= z1:
            hs.append(ax.plot([th["x"]], [th["r"]], "o", ms=4, color="C3",
                              label=f"物理スロート (最小径) x = {th['x']:.3g} mm")[0])
        hs.append(ax.axvline(0.0, color="0.5", ls=":", lw=0.9, label=lab_0))
        ax.set_xlabel(xl)
        ax.set_ylabel("半径 r [mm]")
        ax.set_xlim(z0, z1)
        ax.grid(True, lw=0.3, alpha=0.6)
        ax2 = ax.twinx()
        hs.append(ax2.plot(x[m], dr[m], "-", color="C2", lw=1.2, label="Δr = r_物理 − r_設計 (右の軸)")[0])
        ax2.axhline(0.0, color="C2", lw=0.5, ls="--", alpha=0.6)
        ax2.set_ylabel("Δr = r_物理 − r_設計 [mm]", color="C2")
        ax2.tick_params(axis="y", colors="C2")
        ax.set_title(f"{title}: スロート付近の拡大 (x {z0:g}〜{z1:g} mm)")
        _legend_below(ax, handles=hs, ncol=3)
        fig.savefig(out_dir / PNG_THROAT, dpi=dpi, bbox_inches="tight")
        files["throat"] = PNG_THROAT
    return {"files": files, "font": font, "zoom_mm": [z0, z1], "dpi": dpi, "matplotlib": matplotlib.__version__}


# ---------------------------------------------------------------------------------------------------------------------
# まとめて書く
# ---------------------------------------------------------------------------------------------------------------------

def export_compare(run, out_dir, rel_to=None) -> dict:
    """表・照合・図を `out_dir` に書き、説明の JSON (`JSON_NAME`) も書く。戻り: 成果物の記録 (道は `rel_to` からの相対; 配列は含めない)。
    照合 (`delta_record_check`) が FAIL なら `pass` False、PASS なら True、照合できない点があれば None (合格にしない・不合格にもしない)。"""
    run, out_dir = Path(run), Path(out_dir)
    rel = (lambda p: str(Path(p).relative_to(rel_to))) if rel_to is not None else (lambda p: str(Path(p).name))  # noqa: E731
    W = _load(run)
    tab = compare_table(W, out_dir / CSV_NAME, run_name=run.name)
    chk = delta_record_check(run, tab, W)
    figs = plot_compare(tab, out_dir, run.name)
    pub = {k: v for k, v in tab.items() if not k.startswith("_")}
    pub["csv"] = rel(out_dir / CSV_NAME)
    doc = {"format": "forge.wall_design_vs_physical/v1", "run": str(run), "wall_file": W["path"],
           "wall_file_sha256": _sha(Path(W["path"])), **pub, "delta_record_check": chk,
           "figures": {k: rel(out_dir / v) for k, v in figs["files"].items()}, "figure_font": figs["font"],
           "zoom_mm": figs["zoom_mm"]}
    (out_dir / JSON_NAME).write_text(json.dumps(doc, indent=1, ensure_ascii=False, default=float))
    tp, ex = chk["points"]["throat"], chk["points"]["exit"]
    return {"csv": pub["csv"], "json": rel(out_dir / JSON_NAME), "figures": doc["figures"], "n_rows": pub["n_rows"],
            "domain_mm": pub["domain_mm"], "n_blank": pub["n_blank"], "dr_throat_mm": tp["dr_table_mm"], "dr_exit_mm": ex["dr_table_mm"],
            "dr_range_mm": pub["dr_range_mm"], "physical_throat_mm": pub["physical_throat_mm"],
            "delta_record_check": chk, "notes": pub["notes"], "figure_font": figs["font"], "pass": chk["pass"]}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("out_dir")
    a = ap.parse_args(argv)
    r = export_compare(a.run_dir, a.out_dir)
    print(json.dumps({k: r[k] for k in ("csv", "figures", "dr_throat_mm", "dr_exit_mm", "pass")}, indent=1, ensure_ascii=False))
    print(f"照合: {r['delta_record_check']['status']}")
    return 0 if r["pass"] is not False else 1


if __name__ == "__main__":
    sys.exit(main())
