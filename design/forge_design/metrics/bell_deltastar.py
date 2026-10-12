r"""③ベル (M) の排除厚さ: NS からの $\delta_r(x)$ 抽出の受理と、CONTUR 積分法 (`contur_v2`) の未較正の分布。

計画: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #7 (7b-1)・(7b-3) 前半、codex diagnose 2026-10-12
(notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md の M1〜M5・m1)。

- **M1 一定 x の断面**: ③の NS 格子は近壁が壁法線の層 (`mesh.wall_normal_layer`) で、列の中の x が最大 0.28 r_t ずれる。
  ③の呼び出しはここを通し、`remap_constant_x=True` を既定にする (①風洞の `deltastar_from_core_matched_euler` の既定は変えない)。
  流量の帳簿 (`massflow_ratio`) も一定 x で取り直す (③の prepare_info に x_E が無い件は `deltastar._exit_x` が wall_design.L で扱う)。
- **M2 帯位置の感度**: 同じ壁 (delta_in = 0) では k_band の感度系列が主系列とビット一致して感度にならない。各断面で実際に
  選ばれた帯位置 $y_b$ を倍率 `YB_SENS_FACTORS` 倍に固定して測り直し、主系列との相対差の最大を**抽出不確かさ**とする
  (受理線 `ACCEPT_REL` = 5 %、事前登録)。
- **M3 時系列**: 途中のスナップショットを同じ抽出器で読み、代表 x の $\delta_r$ と全体の L2 を CSV にする
  (`check_quasisteady.py --series-csv` で判定)。
- **M5 半径方向の定義**: 抽出値は一定 x の半径方向等価厚さ (環状の質量等価)。CONTUR の $\delta^*_n$ (法線方向、横曲率重み;
  欠損流量 $2\pi\rho_e u_e r_w\delta^*_n$) を同じ定義に直す: 壁での軸方向の流量密度 $q_w=\rho_e u_e\cos\phi_w$ で
  $2\pi q_w(r_w\delta_r-\delta_r^2/2)=2\pi\rho_e u_e r_w\delta^*_n$ →
  $\delta_r=r_w-\sqrt{r_w^2-2r_w\delta^*_n/\cos\phi_w}$ (1 次は従来の `integral_bl` の $\delta^*_n/\cos\phi_w$)。
- **m1 `contur_v2` の接続**: ③の壁 (`RaoBellWall`)・排気 (`RocketExhaust` の MOC 気体)・NS と同じ輸送物性の μ(T)・断熱壁、
  縁の状態は **Euler 基準の run の壁の節点の M を全域で** (`integral_bl(edge_M_full=...)`)。較正 (kδ) はしない —
  比 $\delta_{r,NS}/\delta_{r,CONTUR}$ の分布を数値で出すだけ (`ratio_stats`)。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import h5py
import numpy as np

from .deltastar import deltastar_from_core_matched_euler, massflow_ratio

YB_SENS_FACTORS = (0.75, 1.5, 2.0)          # 帯位置の倍率 (codex M2 の受理の判定に使う系列; 2026-10-12 主セッション指定)
ACCEPT_REL = 0.05                           # 抽出不確かさの受理線 (plan §5.1 #7 の事前登録)
STATIONS = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)  # 代表 x [r_t] (+ 出口)。codex の再計算の表と同じ
UNC_BANDS = ((-np.inf, 0.0), (0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, 16.0), (16.0, np.inf))
RATIO_REGIONS = ((0.5, None), (1.0, None), (2.0, None), (4.0, None))   # 比の統計の区間 [x_lo, 出口]


# ------------------------------------------------------------------------------------------------ 抽出 (M1・M2)
def extract_bell(ns_run, euler_run, *, ns_res: str | None = None, yb_sens_factors=YB_SENS_FACTORS,
                 remap_constant_x: bool = True, accept_rel: float = ACCEPT_REL, **kw) -> dict:
    """③の NS run から δ_r(x) を抽出する (band 法・adaptive、既定は一定 x の断面)。

    yb_sens_factors: 帯位置の倍率 (None で計算しない — 時系列用)。抽出器の soft の理由 `yb_sens` を相対だけで判定するため
    `gate_sens_floor` の既定を 0 にする (①の 1e-3 r_t の床は③のスロート近くの δ_r ≈ 0.002 r_t には大きすぎる)。
    戻り値は `deltastar_from_core_matched_euler` の dict に `uncertainty` (`extraction_uncertainty`) を足したもの。"""
    kw.setdefault("gate_sens_floor", 0.0)
    out = deltastar_from_core_matched_euler(ns_run, euler_run, remap_constant_x=remap_constant_x, ns_res=ns_res,
                                            yb_sens_factors=(None if yb_sens_factors is None else tuple(yb_sens_factors)), **kw)
    if remap_constant_x:
        # 流量の帳簿も一定 x の断面で (列のままは③ L100 で 3 % 小さい — massflow_ratio の docstring)。列のままの値は記録に残す
        out["massflow_column"] = out["massflow"]
        out["massflow"] = massflow_ratio(ns_run, euler_run, ns_res=ns_res, remap_constant_x=True)
    out["uncertainty"] = (extraction_uncertainty(out, accept_rel=accept_rel) if yb_sens_factors is not None
                          else {"computed": False, "reason": "yb_sens_factors=None (時系列などで帯位置の感度を省いた)"})
    return out


def extraction_uncertainty(out: dict, accept_rel: float = ACCEPT_REL, x_lo: float = 0.0, bands=UNC_BANDS) -> dict:
    r"""帯位置の感度から各断面の抽出不確かさ $u(x)=\max_f|\delta_{r,f}-\delta_r|/\delta_r$ を作り、受理を判定する。

    判定の区間は $x\ge x_{lo}$ (既定 0 = スロートから出口; 収縮部は補正しないので記録だけ)。受理 = 区間の全断面で
    $u\le$ accept_rel かつ非有限なし。参考として倍率ごとの符号つきの差と、外向き (f > 1) だけの最大も返す
    (判定には使わない — 登録の系列は `settings.yb_sens_factors` の全部)。"""
    x = np.asarray(out["x"], dtype=float)
    d = np.asarray(out["delta_r_raw"], dtype=float)
    S = np.asarray(out["delta_r_sens_yb"], dtype=float)
    facs = [float(f) for f in out["settings"]["yb_sens_factors"]]
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = (S - d[:, None]) / d[:, None]
    finite_all = np.all(np.isfinite(rel), axis=1)
    u = np.where(finite_all, np.max(np.abs(np.where(np.isfinite(rel), rel, 0.0)), axis=1), np.nan)
    outward = [k for k, f in enumerate(facs) if f > 1.0]
    u_out = (np.where(np.all(np.isfinite(rel[:, outward]), axis=1), np.max(np.abs(np.nan_to_num(rel[:, outward])), axis=1), np.nan)
             if outward else np.full(len(x), np.nan))
    reg = x >= x_lo
    n_nan = int(np.sum(~np.isfinite(u[reg])))
    u_reg = u[reg]
    i_max = int(np.nanargmax(np.where(reg, u, -np.inf))) if np.any(np.isfinite(u_reg)) else -1

    def band_row(a, b):
        m = (x >= a) & (x < b)
        if not m.any():
            return None
        row = {"x_lo": (None if not np.isfinite(a) else float(a)), "x_hi": (None if not np.isfinite(b) else float(b)),
               "n": int(m.sum()), "u_max": float(np.nanmax(u[m])) if np.any(np.isfinite(u[m])) else float("nan"),
               "frac_over": float(np.mean(~(u[m] <= accept_rel))), "u_outward_max": (float(np.nanmax(u_out[m]))
                                                                                   if np.any(np.isfinite(u_out[m])) else float("nan"))}
        for k, f in enumerate(facs):
            rk = rel[m, k]
            row[f"rel_f{f:g}_min"] = float(np.nanmin(rk)) if np.any(np.isfinite(rk)) else float("nan")
            row[f"rel_f{f:g}_max"] = float(np.nanmax(rk)) if np.any(np.isfinite(rk)) else float("nan")
        return row

    rows = [r for r in (band_row(a, b) for a, b in bands) if r is not None]
    accepted = bool(n_nan == 0 and len(u_reg) > 0 and np.nanmax(u_reg) <= accept_rel)
    return {"computed": True, "accept_rel": accept_rel, "region": f"x >= {x_lo:g} (スロート〜出口)", "factors": facs,
            "accepted": accepted, "n_region": int(reg.sum()), "n_nonfinite": n_nan,
            "n_over": int(np.sum(~(u_reg <= accept_rel))),
            "u_max": (float(u[i_max]) if i_max >= 0 else float("nan")), "x_at_u_max": (float(x[i_max]) if i_max >= 0 else float("nan")),
            "u_outward_max_region": (float(np.nanmax(u_out[reg])) if np.any(np.isfinite(u_out[reg])) else float("nan")),
            "bands": rows, "u": u, "u_outward": u_out, "rel": rel}


def _stop_ok(slope, rms, share, a, alpha, d, reason, settings, a_frac: float = 0.95) -> dict:
    """帯の停止条件 (`band_local_deficit` の適応の再試行の条件と同じ量) を満たすか — 要素ごとの bool の dict。
    slope: 帯の中の比の変化 |Δ| ≤ slope_max、a_core: 帯の比の切片 a ≥ 0.95 α (帯が境界層の深部でない)、clean: 帯の中の残留欠損の比
    ≤ 0.1、fit: 帯の当てはめの RMS ≤ gate_band_rms、root: δ_r が有限で no_root / nan / 帯が取れない (none) でない。"""
    slope, rms, share, a, alpha, d = (np.asarray(v, dtype=float) for v in (slope, rms, share, a, alpha, d))
    hard = np.array([any(h in str(r).split(",") for h in ("no_root", "nan", "none")) for r in np.ravel(reason)],
                    dtype=bool).reshape(np.shape(reason))
    with np.errstate(invalid="ignore"):
        return {"slope": np.abs(slope) <= float(settings["slope_max"]),
                "a_core": ~(np.isfinite(alpha) & (a < a_frac * alpha)),
                "clean": ~(np.isfinite(share) & (np.abs(share) > 0.1)),
                "fit": rms <= float(settings["gate_band_rms"]),
                "root": np.isfinite(d) & ~hard}


def band_acceptance(out: dict, accept_rel: float = ACCEPT_REL, x_lo: float = 0.0) -> dict:
    r"""帯の選び方の A/B (plan §5.1 #7「次の A/B (事前登録)」、codex diagnose 2026-10-12 の受理の諮問) の判定に使う断面ごとの条件。

    区間 $x\ge x_{lo}$ の各断面で:
      - `sens`: 帯位置の感度 $u\le$ accept_rel (登録の倍率すべてで有限)
      - `main_stop`: 主系列の帯が停止条件 (`_stop_ok`) をすべて満たす
      - `sens_stop`: 感度帯 (全倍率) が停止条件をすべて満たす
      - (band_select="edge" のときだけ) `edge_found`: 局所縁 $y^\ast$ (`edge_y_star`、平滑化の前) が有限、
        `band_outside_edge`: 最終の帯の下端 $y_b\ge y^\ast$ (帯が局所縁の外)
    `all_pass` = 全条件を区間の全断面で満たす (A/B の当てはめは親が事前登録に対して行う)。edge のときは平滑化の前後の縁の位置
    (局所 $y^\ast$、移動中央値、P-spline + 端の保持、$y_b=c\cdot$平滑縁) と端の保持の区間 (`held`) も断面ごとに返す
    (`edge_detail` が要る — `return_ladders=True` で抽出する)。"""
    st = out["settings"]
    x = np.asarray(out["x"], dtype=float)
    reg = x >= x_lo
    unc = out["uncertainty"] if (isinstance(out.get("uncertainty"), dict) and out["uncertainty"].get("computed")) else \
        extraction_uncertainty(out, accept_rel=accept_rel, x_lo=x_lo)
    u = np.asarray(unc["u"], dtype=float)
    crit = {"sens": np.isfinite(u) & (u <= accept_rel)}
    main = _stop_ok(out["band_slope"], out["band_rms"], out["band_deficit_share"], out["band_ratio_a"], out["alpha"],
                    out["delta_r_raw"], out["reason"], st)
    crit["main_stop"] = np.logical_and.reduce(list(main.values()))
    facs = [float(f) for f in st["yb_sens_factors"]]
    sens_parts = {}
    for k, f in enumerate(facs):
        sens_parts[f] = _stop_ok(out["yb_sens_band_slope"][:, k], out["yb_sens_band_rms"][:, k], out["yb_sens_band_deficit_share"][:, k],
                                 out["yb_sens_band_ratio_a"][:, k], out["yb_sens_alpha"][:, k], out["delta_r_sens_yb"][:, k],
                                 out["yb_sens_reason"][:, k], st)
    crit["sens_stop"] = np.logical_and.reduce([np.logical_and.reduce(list(sp.values())) for sp in sens_parts.values()])
    edge = None
    if st.get("band_select") == "edge":
        ys = np.asarray(out["edge_y_star"], dtype=float)
        yb = np.asarray(out["band_y_b"], dtype=float)
        crit["edge_found"] = np.isfinite(ys)
        crit["band_outside_edge"] = np.isfinite(ys) & (yb >= ys)
        ed = out.get("edge_detail")
        if ed is not None:
            ex = np.asarray(ed["x"], dtype=float)
            pos = np.clip(np.searchsorted(ex, x), 0, len(ex) - 1)
            if not np.array_equal(ex[pos], x):
                raise ValueError("edge_detail の x が抽出の断面と対応しない")
            y_ref = np.asarray(ed["y_ref"], dtype=float)[pos]
            edge = {"x": x, "y_star": ys, "y_ref": y_ref, "y_edge_median": y_ref * np.exp(np.asarray(ed["kappa_med"], dtype=float)[pos]),
                    "y_edge_smooth": y_ref * np.exp(np.asarray(ed["kappa_s"], dtype=float)[pos]), "y_b": yb,
                    "c": float(st["edge"]["c"]),
                    "held": (x < 3.0) | (x > float(ex.max()) - 3.0)}       # edge_band_positions の端の保持 (hold_lo 3、出口から 3)
    keys = list(crit)
    all_c = np.logical_and.reduce([crit[k] for k in keys])

    def where_fail(m):
        bad = reg & ~np.asarray(m, dtype=bool)
        if not bad.any():
            return {"n_fail": 0}
        xb = x[bad]
        return {"n_fail": int(bad.sum()), "x_first": float(xb.min()), "x_last": float(xb.max())}
    summary = {k: where_fail(crit[k]) for k in keys}
    summary["main_stop_parts"] = {k: where_fail(v) for k, v in main.items()}
    summary["sens_stop_parts"] = {f"f{f:g}": {k: where_fail(v) for k, v in sp.items()} for f, sp in sens_parts.items()}
    summary["all"] = where_fail(all_c)
    return {"region": f"x >= {x_lo:g}", "n_region": int(reg.sum()), "accept_rel": accept_rel, "factors": facs,
            "criteria": keys, "all_pass": bool(np.all(all_c[reg])), "summary": summary,
            "crit": crit, "all": all_c, "main_parts": main, "sens_parts": sens_parts, "edge": edge}


def station_values(x, v, stations=STATIONS) -> dict:
    """代表 x (線形補間) と出口 (最後の断面) の値。"""
    x = np.asarray(x, dtype=float); v = np.asarray(v, dtype=float)
    out = {f"x{s:g}": float(np.interp(s, x, v)) for s in stations}
    out["exit"] = float(v[-1])
    return out


def remap_comparison(ns_run, euler_run, stations=STATIONS, **kw) -> dict:
    """旧 (列のまま、`remap_constant_x=False`) と新 (一定 x) の抽出の差 (codex M1 の再計算の表)。帯位置の感度は省く。"""
    old = deltastar_from_core_matched_euler(ns_run, euler_run, remap_constant_x=False, **kw)
    new = deltastar_from_core_matched_euler(ns_run, euler_run, remap_constant_x=True, **kw)
    if not np.array_equal(old["x"], new["x"]):
        raise ValueError("旧と新で断面の x が違う (列の x = 軸の x のはず)")
    rows = []
    for key in list(stations) + ["exit"]:
        xs = float(old["x"][-1]) if key == "exit" else float(key)
        ro = float(np.interp(xs, old["x"], old["delta_r_raw"])); rn = float(np.interp(xs, new["x"], new["delta_r_raw"]))
        so = float(np.interp(xs, old["x"], old["delta_r_smooth"])); sn = float(np.interp(xs, new["x"], new["delta_r_smooth"]))
        rows.append({"station": ("exit" if key == "exit" else f"{key:g}"), "x": xs, "raw_column": ro, "raw_constx": rn,
                     "raw_rel_change": rn / ro - 1.0, "smooth_column": so, "smooth_constx": sn, "smooth_rel_change": sn / so - 1.0})
    m = new["x"] >= 0.0
    rel = new["delta_r_raw"][m] / old["delta_r_raw"][m] - 1.0
    return {"rows": rows, "max_rel_change_expansion": float(np.max(np.abs(rel))),
            "x_at_max": float(new["x"][m][int(np.argmax(np.abs(rel)))]),
            "remap_max_column_dx": float(new["settings"]["remap_max_column_dx"]),
            "x": new["x"], "raw_column": old["delta_r_raw"], "raw_constx": new["delta_r_raw"],
            "smooth_column": old["delta_r_smooth"], "smooth_constx": new["delta_r_smooth"]}


# ------------------------------------------------------------------------------------------------ 時系列 (M3)
def snapshot_files(run_dir) -> list:
    """`res_<n>.h5` (段の res だけ; res_wall_* 等は除く) を step 順に [(step, name)]。"""
    out = []
    for f in Path(run_dir).glob("res_*.h5"):
        m = re.fullmatch(r"res_(\d+)\.h5", f.name)
        if m:
            out.append((int(m.group(1)), f.name))
    return sorted(out)


def deltastar_series(ns_run, euler_run, out_csv, stations=STATIONS, min_step: int | None = None,
                     out_npz=None, **kw) -> dict:
    r"""NS run の途中のスナップショットを順に抽出 (`extract_bell`、帯位置の感度は省く) して $\delta_r$ の時系列を CSV にする。

    列: step, res, dr_x<s> (代表 x の生値), dr_exit, dr_l2 (スロート〜出口 [x ≥ 0] の生値の RMS)、dr_smooth_l2 (平滑値の RMS)、
    n_hard_fail (hard 不合格の断面数)。Euler 基準は最終場で固定 (`_load_structured` の最新)。out_npz を与えると全断面の
    δ_r(x) を step ごとに保存する。戻り値: dict(rows, csv, cols, command) — command は `check_quasisteady.py` の雛形
    (drift・osc の閾値は呼び出し側が決める)。"""
    snaps = [(s, n) for s, n in snapshot_files(ns_run) if (min_step is None or s >= min_step)]
    if not snaps:
        raise FileNotFoundError(f"{ns_run}: res_<n>.h5 が無い")
    cols = [f"dr_x{s:g}".replace(".", "p") for s in stations] + ["dr_exit", "dr_l2", "dr_smooth_l2"]
    rows, X, D, DS = [], None, [], []
    for step, name in snaps:
        o = extract_bell(ns_run, euler_run, ns_res=name, yb_sens_factors=None, **kw)
        if X is None:
            X = o["x"]
        elif not np.array_equal(X, o["x"]):
            raise ValueError(f"{name}: 断面の x が前のスナップショットと違う")
        m = o["x"] >= 0.0
        sv = station_values(o["x"], o["delta_r_raw"], stations)
        rows.append({"step": step, "res": name, **{c: sv[k] for c, k in zip(cols, [f"x{s:g}" for s in stations] + ["exit"])},
                     "dr_l2": float(np.sqrt(np.mean(o["delta_r_raw"][m] ** 2))),
                     "dr_smooth_l2": float(np.sqrt(np.mean(o["delta_r_smooth"][m] ** 2))),
                     "n_hard_fail": int(np.sum(~o["hard_ok"]))})
        D.append(o["delta_r_raw"]); DS.append(o["delta_r_smooth"])
    out_csv = Path(out_csv)
    hdr = ["step", "res"] + cols + ["n_hard_fail"]
    with open(out_csv, "w") as fh:
        fh.write(",".join(hdr) + "\n")
        for r in rows:
            fh.write(",".join([str(r["step"]), r["res"]] + [repr(float(r[c])) for c in cols] + [str(r["n_hard_fail"])]) + "\n")
    if out_npz is not None:
        np.savez(out_npz, step=np.array([r["step"] for r in rows]), x=X, delta_r_raw=np.array(D), delta_r_smooth=np.array(DS))
    cmd = (f"python3 solver_density_cuda/tools/check_quasisteady.py --series-csv {out_csv} --series-cols {','.join(cols)} "
           "--drift <D> --osc <O>")
    return {"rows": rows, "csv": str(out_csv), "cols": cols, "command": cmd}


# ------------------------------------------------------------------------------------------------ CONTUR (m1・M5)
def euler_wall_edge(euler_run, gas) -> dict:
    """Euler 基準の run (node 構造格子、slip 壁) の壁の節点 (j = nj−1) の状態: x, r [r_t], M = |U|/√(γ(T) R T) (MOC 気体)。
    検算: 壁の T と等エントロピーの T(M) (超音速側) の差、壁の h0 の相対の振れ。"""
    run_dir = Path(euler_run)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj, S = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]), float(info["scale_m"])
    snaps = snapshot_files(run_dir)
    if not snaps:
        raise FileNotFoundError(f"{run_dir}: res_<n>.h5 が無い")
    res = snaps[-1][1]
    with h5py.File(run_dir / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(run_dir / res) as f:
        V = {k: f["/VALUE/" + k][:] for k in ("Ux", "Uy", "T", "h0")}
    if len(nc) != ni * nj or any(len(v) != ni * nj for v in V.values()):
        raise ValueError(f"{run_dir}: node 構造格子 run でない")
    sl = slice(nj - 1, ni * nj, nj)                      # 壁の行 (index = i*nj + j)
    x = nc[sl, 0] / S; r = nc[sl, 1] / S
    T = V["T"][sl]; U = np.hypot(V["Ux"][sl], V["Uy"][sl])
    M = U / np.sqrt(np.asarray(gas.gamma(T), dtype=float) * float(gas.R) * T)
    sup = M >= 1.0
    dT = (np.asarray(gas.T_of_M(np.maximum(M[sup], 1.0)), dtype=float) - T[sup]) if sup.any() else np.zeros(0)
    h0 = V["h0"][sl]
    return {"run": str(run_dir), "res": res, "x": x, "r": r, "M": M, "T": T,
            "check": {"max_abs_T_minus_T_of_M_supersonic_K": float(np.max(np.abs(dT))) if dT.size else float("nan"),
                      "h0_rel_spread": float(np.ptp(h0) / np.mean(h0)), "x_first": float(x[0]), "x_last": float(x[-1]),
                      "M_at_x0": float(np.interp(0.0, x, M))}}


def delta_r_annular(dstar_n, cos_phi, rw):
    r"""CONTUR の $\delta^*_n$ を一定 x の環状の質量等価厚さへ (M5; docstring の式)。根が無い (2δ*/cos φ > r_w) は NaN。"""
    dstar_n = np.asarray(dstar_n, dtype=float); cos_phi = np.asarray(cos_phi, dtype=float); rw = np.asarray(rw, dtype=float)
    disc = rw * rw - 2.0 * rw * dstar_n / cos_phi
    return np.where(disc >= 0.0, rw - np.sqrt(np.maximum(disc, 0.0)), np.nan)


def contur_bell(problem_yaml, euler_run, *, x_out=None, x_virtual_m: float | None = None, theta0_m: float | None = None,
                edge: str = "euler_full", edge_n_grid: int = 4000, rtol: float = 1e-6) -> dict:
    r"""③の壁・排気・断熱壁・Euler 基準の縁で `contur_v2` を回し、未較正の $\delta_{r,CONTUR}(x)$ (M5 の環状の定義と従来の 1 次) を返す。

    problem_yaml: NS の問題 (壁ファイル・gas.exhaust・gas.transport)。壁は `runner.build_wall` (= NS run の壁; 照合は呼び出し側)。
    edge: "euler_full" (既定; 縁の M を全域で Euler の壁の節点から) / "euler_1d_upstream" (感度用: 従来の EdgeConditions の
    上流 1D 面積比 + x > 0 は Euler の壁の M — スロートで M が跳ぶ)。x_virtual_m: None は入口直管長 L_pipe·r_t
    (`integral_bl` の文書の既定。`RaoBellWall` は L_pipe 属性を持たないので明示する)。x_out: 出す x [r_t] (積分区間に切る)。"""
    from ..evaluate.runner import build_wall, rao_bell_spec
    from ..evaluate.runner_axismach import contur_mu_fn
    from ..feedback.deltastar_integral import integral_bl
    from ..probdef import load_problem
    if edge not in ("euler_full", "euler_1d_upstream"):
        raise ValueError(f"edge {edge!r} は euler_full | euler_1d_upstream")
    p = load_problem(problem_yaml)
    wall, rec = build_wall(p)
    gas = p.gas_model
    rt_m = float(p.spec["r_throat"])
    tbc = p.wall_thermal_bc_integral
    if str(tbc.get("mode", "adiabatic")) != "adiabatic":
        raise ValueError(f"③の較正の基準ケースは断熱壁 (問題の壁は {tbc})")
    ed = euler_wall_edge(euler_run, gas)
    x_F = min(float(ed["x"][-1]), float(wall.x_e))
    m = ed["x"] <= x_F
    tbl_full = np.c_[ed["x"][m], ed["M"][m]]
    wall_tbl = np.c_[ed["x"][m], ed["r"][m], np.zeros(int(m.sum())), ed["M"][m]]
    if x_virtual_m is None and theta0_m is None:
        x_virtual_m = float(rao_bell_spec(p).L_pipe) * rt_m
    x0 = float(wall.x_in) + 1e-6
    if x_out is None:
        xo = np.linspace(x0, x_F, 1500)
    else:
        # NS の断面の x をそのまま渡せるよう、端の丸め (入口 x_in と積分の始点 x_in + 1e-6、出口の数 ulp) だけは端へ寄せる
        xo = np.asarray(x_out, dtype=float)
        xo = np.clip(xo[(xo >= x0 - 2e-6) & (xo <= x_F + 1e-9 * max(1.0, abs(x_F)))], x0, x_F)
    res = integral_bl(wall, wall_tbl, gas, p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt_m, thermal_bc=tbc,
                      theta0_m=theta0_m, x_virtual_m=x_virtual_m, x_out=xo, rtol=rtol, closure_version="contur_v2",
                      mu_fn=contur_mu_fn(p), edge_n_grid=edge_n_grid,
                      edge_M_full=(tbl_full if edge == "euler_full" else None))
    rw = np.asarray(wall.r(res["x"]), dtype=float)
    res["r_wall"] = rw
    res["delta_r_lin"] = res["delta_r"]                              # 従来の 1 次 δ*_n / cos φ
    res["delta_r_ann"] = delta_r_annular(res["dstar_n"], res["cos_phi"], rw)
    res["edge"] = {"mode": edge, "euler": {k: v for k, v in ed.items() if k in ("run", "res", "check")}}
    res["wall_record"] = {k: rec.get(k) for k in ("wall_file", "wall_file_sha256") if k in rec}
    res["problem"] = str(problem_yaml)
    return res


def ratio_stats(x, d_ns, d_c, regions=RATIO_REGIONS, x_exit: float | None = None) -> list:
    r"""比 $R=\delta_{r,NS}/\delta_{r,CONTUR}$ の区間ごとの分布 (較正ではない — 定数 1 個で表せそうかの数値だけ)。

    各区間 [x_lo, 出口]: n、中央値・平均・最小・最大・最大/最小・相対標準偏差、参考の最良定数
    $k^\ast=\arg\min_k\sum(k\,\delta_C/\delta_{NS}-1)^2=\sum a/\sum a^2$ ($a=\delta_C/\delta_{NS}$) とそのときの相対 RMS・局所の最大。"""
    x = np.asarray(x, dtype=float); R = np.asarray(d_ns, dtype=float) / np.asarray(d_c, dtype=float)
    xe = float(x[-1]) if x_exit is None else float(x_exit)
    out = []
    for lo, hi in regions:
        hi_ = xe if hi is None else float(hi)
        m = (x >= lo) & (x <= hi_) & np.isfinite(R) & (R > 0)
        if not m.any():
            continue
        r = R[m]; a = 1.0 / r
        k = float(np.sum(a) / np.sum(a * a))
        e = k * a - 1.0
        out.append({"x_lo": float(lo), "x_hi": hi_, "n": int(m.sum()), "median": float(np.median(r)), "mean": float(np.mean(r)),
                    "min": float(np.min(r)), "max": float(np.max(r)), "max_over_min": float(np.max(r) / np.min(r)),
                    "rel_std": float(np.std(r) / np.mean(r)), "x_at_min": float(x[m][int(np.argmin(r))]),
                    "x_at_max": float(x[m][int(np.argmax(r))]),
                    "ref_best_const": k, "ref_const_rel_rms": float(np.sqrt(np.mean(e * e))), "ref_const_rel_maxabs": float(np.max(np.abs(e)))})
    return out
