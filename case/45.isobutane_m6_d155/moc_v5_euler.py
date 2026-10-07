"""plan discretization-moc-axis-limit-and-corrector §5.1 #5 の 5a・§6 V5 (2026-10-07): Euler の確認の run 準備・実行と IC の検査。
腕 M (mocG1) = problem_d155_euler_pin_G1_recal_mono_moc.yaml (腕 B の単調壁の問題 + geometry.moc_axis_limit: analytic・
geometry.moc_corrector: converge)。生産 Euler 格子 G1。比べる相手は腕 B (run_0143〜0145、throat_mono_ab.py の monoG1)。
  r1〜r3 (run_0150〜0152): IC = run_0114 の最終場 (res_6000) の保存済み保存量を、検証付き番号写像 (ic_index_map --mode index) で写す。
    設計壁が最大約 6.2 µm 動き、前回の番号写像の上限 (1 µm) を超える。§6 V5 のとおり上限の数値だけを上げて通さない:
    先に変換後メッシュを検査 (inspect_ic) して IC_INSPECTION.json に記録し、下の条件 C1〜C5 がすべて成立したときだけ、
    その構造から決まる上限で写す。
  IC 依存の確認 (run_0153): 腕 M の prep を IC を入れる前に複製し (同じ問題・同じ格子)、RA.prepare の等エントロピー IC のまま回す。
起動は throat_mono_ab.py と同じ RA.run_staged(stages="soft"): soft (1 次 cfl 0.5、3000 step) → 本段 2 次 cfl 2・implicitRelax 0.7・
18000 step・1000 ごと出力。

番号写像を使う条件 (IC の検査、inspect_ic。結果を見て変えない):
  C1 ic_index_map と同じ検査 (節点数・接続、論理位置 (i, j)、境界種別と節点集合・境界値、座標系と単位、化学種とエネルギー基準、
     要素の反転・ねじれ・退化なし) がすべて成立。移動量は測るだけ (検査のための上限 1 m は写像には使わない)。
  C2 列の x が不変: 全節点で |Δx| ≤ ulp32(x) + |ΔX_F|·r_t (格子の x 配置は入口〜設計の出口 x_F の関数なので、x_F の差だけ動きうる)。
  C3 半径方向の移動が、列ごとの壁の移動を η = r/r_w を保って縮尺したものだけ:
     全節点で |r_DST − η_SRC·r_w,DST| ≤ ½(ulp(r_SRC) + ulp(r_DST)) + ½·η_SRC·(ulp(r_w,SRC) + ulp(r_w,DST)) (float32 の丸め)。
     これから |Δr| ≤ |Δr_w| + 2·ulp(r_w) (内側の節点は同じ列の壁より動かない)。
  C4 壁の移動が設計壁の変化で説明される: 腕 M の壁節点が腕 M の当てはめ後 spline に一致し (M4)、腕 B の参照 (run_0143) の壁節点が
     腕 B の spline に一致する (= 腕 M − 腕 B の壁の差は 2 本の設計壁の差そのもの)。両 spline の範囲外では腕 M と腕 B の壁が
     float32 の丸め以内で同じ。腕 B と IC 格子の壁の差は、腕 B の写像で受けた上限 (腕 B の IC_MAP.json の limit_m) 以下で、
     腕 B の写像は同じ IC (run_0114 の res_6000) から取ったもの。
  C5 腕 M の変換後メッシュの品質が PASS (MESH_QUALITY.txt)。
上限: limit_m = max_i √(tol_x,i² + (|Δr_w,i| + 2·ulp(r_w,i))²) (C2・C3 が成り立つときの列ごとの移動の上界)。数値を先に決めて
  広げたものではない。番号写像は各節点を同じ論理位置 (i, j) かつ同じ正規化位置 η に写す。最近傍は壁の移動 (µm) が壁際の
  セル厚 (0.35 µm) を超える所で層を取り違えるので使わない。番号写像で残る差 (体積比、Δr_w/r_w 程度) は IC 依存の run で確かめる。

usage: [CASE_RUNS=<run_0062・run_0114 のある case dir>] python3 moc_v5_euler.py prep <prep_dir> [--with-isen <prep_dir_isen>]
           [--armB-ref <腕 B の run (既定 run_0143_euler_wallfit_monoG1_r1)>] [--dry]
       python3 moc_v5_euler.py run <run_dir>
       python3 moc_v5_euler.py verify-prep <prep_dir> [<run_dir> ...]
--dry: ローカルの乾式確認 (forge を起動しない)。FORGE_BIN を存在しない道に向け、FORGE_ALLOW_UNVERIFIED_SPECIES=1 にする
  (化学種の属性は付かない)。ic_index_map は --no-species-resolve。その prep は prepare_info に DRY の印が付き、run・verify-prep が拒否する。
  IC の検査の記録は _band_ab/moc_v5_ic_inspection_dry.json (本番は _band_ab/moc_v5_ic_inspection.json)。
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import throat_mono_ab as TM  # noqa: E402  (run_0114 の場所・壁の証拠・prep の照合を共有する)

PLAN = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5"
PROBLEM_M = "problem_d155_euler_pin_G1_recal_mono_moc.yaml"
ARM_M, ARM_ISEN = "mocG1", "mocG1_isen"
IC_RES = "res_6000.h5"
ARMB_REF = "run_0143_euler_wallfit_monoG1_r1"
MONO_R2 = [0.0, 1.5]
MOC_EXPECT = {"axis_limit": "analytic", "corrector": "converge"}
INSPECT_DISP_M = 1.0                      # 検査で移動量を測るための ic_index_map の上限 (写像には使わない)
C1_CHECKS = ("connectivity", "logical_ij", "boundaries", "coordinate_system_units", "species_energy", "no_inversion", "displacement")
DISP_BINS_M = [0.0, 1e-9, 1e-8, 1e-7, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 1e-5, math.inf]
X_BANDS_RT = [-math.inf, -1.0, 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, math.inf]
V4_DESIGN_WALL_CHANGE = {"max_um": 6.18, "x_rt": 0.40, "source": "plan §9 2026-10-07 V4 (腕 M − 腕 B の設計壁、記録用)"}


def _ulp(a):
    """float32 の 1 ulp (座標は float32 で保存されている)。float64 で返す。"""
    return np.spacing(np.abs(np.asarray(a, dtype=np.float64)).astype(np.float32)).astype(np.float64)


def _coords(h5path):
    with h5py.File(h5path, "r") as f:
        raw = f["/MESH/COORD"][:]
    return raw.reshape(-1, 3).astype(np.float64), str(raw.dtype)


def _quality(path: Path) -> dict:
    """MESH_QUALITY.txt → {verdict, ar_max, skew_max} (読めなければ verdict None)。"""
    out = {"file": str(path), "verdict": None, "ar_max": None, "skew_max": None}
    if not path.is_file():
        return out
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("aspect ratio") and "max=" in s:
            out["ar_max"] = float(s.split("max=", 1)[1].split()[0])
        elif s.startswith("skewness") and "max=" in s:
            out["skew_max"] = float(s.split("max=", 1)[1].split()[0])
        elif s.startswith("VERDICT:"):
            out["verdict"] = s.split(":", 1)[1].strip()
    return out


def _x_last_design(run_dir: Path) -> float:
    """設計の出口 x [m] (wall_design.csv の最後の x; float64)。"""
    return float(np.loadtxt(run_dir / "wall_design.csv", delimiter=",", skiprows=1)[-1, 0])


def _bands(x_rt, d, edges):
    out = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x_rt >= lo) & (x_rt < hi)
        out[f"[{lo:g}, {hi:g})"] = {"n": int(m.sum()), "n_moved": int(np.count_nonzero(d[m] > 0)),
                                    "max_um": (float(d[m].max() * 1e6) if m.any() else None)}
    return out


def inspect_ic(src_res: Path, src_run: Path, dst: Path, armB_ref: Path) -> dict:
    """変換後メッシュ (dst/nozzle.h5) と IC 格子 (src_run/nozzle.h5) の検査 (C1〜C5)。何も書き換えない。戻り値 = 記録
    (VERDICT "OK" / "REFUSED"、条件ごとの合否と測定値、上限 limit_m とその理由)。"""
    from ic_index_map import MapRefused, run_checks, structured_shape
    rec = {"plan": PLAN, "tool": "moc_v5_euler.inspect_ic", "src_res": str(src_res), "src_mesh": str(src_run / "nozzle.h5"),
           "dst": str(dst / "nozzle.h5"), "armB_ref": str(armB_ref), "conditions": {}}
    cond = rec["conditions"]

    def put(name, fails, detail):
        cond[name] = {"ok": not fails, "failures": list(fails), "detail": detail}

    # --- C1: ic_index_map と同じ検査 (移動量は測るだけ) ---
    try:
        chk, _ = run_checks(str(src_res), str(dst / "nozzle.h5"), str(src_run / "nozzle.h5"), str(src_run), str(dst),
                            max_disp_m=INSPECT_DISP_M, resolve_species=False)
        c1_fail = []
        checks = chk["checks"]
    except MapRefused as e:
        checks = e.record.get("checks", {})
        c1_fail = list(e.failures)
    disp_rec = (checks.get("displacement") or {}).get("detail")
    # 検査の項目がそろって合格していること (項目の欠落を合格にしない)
    c1_fail += [f"ic_index_map の検査 {n} が{'無い' if n not in checks else '不合格'}" for n in C1_CHECKS
                if (checks.get(n) or {}).get("ok") is not True and not any(x.startswith(f"[{n}]") for x in c1_fail)]
    put("C1_ic_index_map_checks", c1_fail,
        {"checks": {k: {"ok": v.get("ok"), "failures": v.get("failures")} for k, v in checks.items()},
         "no_inversion_detail": (checks.get("no_inversion") or {}).get("detail"),
         "boundaries_detail": (checks.get("boundaries") or {}).get("detail")})
    # --- 座標 (float32 → float64) ---
    cs, dts = _coords(src_run / "nozzle.h5")
    cd, dtd = _coords(dst / "nozzle.h5")
    with h5py.File(dst / "nozzle.h5", "r") as f:
        ni, nj = structured_shape(f["VIZMESH/CONNE"], cd.shape[0])
    S = float(json.loads((dst / "prepare_info.json").read_text())["scale_m"])
    X0, R0 = cs[:, 0].reshape(ni, nj), cs[:, 1].reshape(ni, nj)
    X1, R1 = cd[:, 0].reshape(ni, nj), cd[:, 1].reshape(ni, nj)
    dX, dR = X1 - X0, R1 - R0
    d = np.sqrt(dX ** 2 + dR ** 2)
    # 移動の分布 (ic_index_map の記録より細かい刻み・x の帯・壁からの層)
    k = int(np.argmax(d))
    i_, j_ = divmod(k, nj)
    moved = d[d > 0]
    hist = {f"({lo:g}, {hi:g}]": int(np.count_nonzero((d > lo) & (d <= hi))) for lo, hi in zip(DISP_BINS_M[:-1], DISP_BINS_M[1:])}
    rec["displacement"] = {
        "max_m": float(d.max()), "max_um": float(d.max() * 1e6), "max_rt": float(d.max() / S),
        "argmax": {"i": i_, "j": j_, "layer_from_wall": nj - 1 - j_, "x_rt": float(X1[i_, j_] / S), "r_rt": float(R1[i_, j_] / S)},
        "n_nodes": int(d.size), "n_moved": int(moved.size),
        "quantiles_um": ({q: float(np.quantile(moved, float(q)) * 1e6) for q in ("0.5", "0.9", "0.99")} if moved.size else None),
        "histogram_m": hist, "by_x_band_rt": _bands((X1 / S).ravel(), d.ravel(), X_BANDS_RT),
        "by_layer_from_wall_max_um": [float(d[:, nj - 1 - L].max() * 1e6) for L in range(nj)],
        "dx_max_m": float(np.abs(dX).max()), "dx_n_nonzero": int(np.count_nonzero(dX)),
        "dr_max_m": float(np.abs(dR).max()), "dr_n_nonzero": int(np.count_nonzero(dR)),
        "ic_index_map_record": ({kk: disp_rec.get(kk) for kk in ("max_um", "argmax", "n_moved", "moved_quantiles_m", "histogram_m")}
                                if disp_rec else None)}
    # --- C2: 列の x が不変 ---
    dXF = abs(_x_last_design(dst) - _x_last_design(src_run))
    tol_x = _ulp(np.maximum(np.abs(X0), np.abs(X1))) + dXF
    rx = np.abs(dX) / tol_x
    col_x_const = bool(np.all(X1 == X1[:, :1]) and np.all(X0 == X0[:, :1]))
    f2 = []
    if not rx.max() <= 1.0:
        kk = np.unravel_index(int(np.argmax(rx)), rx.shape)
        f2.append(f"|Δx| が ulp32(x) + |ΔX_F| を超える節点 {int(np.count_nonzero(rx > 1))} 個 (最大 {rx.max():.3g} 倍、i {kk[0]}, j {kk[1]})")
    if not col_x_const:
        f2.append("列の x が j に一定でない (構造格子の列でない)")
    put("C2_columns_x_unchanged", f2, {"dX_F_design_m": dXF, "dx_over_tol_max": float(rx.max()),
                                       "n_dx_nonzero": int(np.count_nonzero(dX)), "dx_max_m": float(np.abs(dX).max())})
    # --- C3: 半径方向は列ごとの壁の移動の縮尺 (η を保つ) ---
    Rw0, Rw1 = R0[:, -1:], R1[:, -1:]
    eta0 = R0 / Rw0
    res3 = np.abs(R1 - eta0 * Rw1)
    tol3 = 0.5 * (_ulp(R0) + _ulp(R1)) + 0.5 * eta0 * (_ulp(Rw0) + _ulp(Rw1))
    r3 = res3 / tol3
    dRw = (Rw1 - Rw0)[:, 0]
    excess = np.abs(dR) - np.abs(dRw)[:, None]
    tol_w = 2.0 * _ulp(np.maximum(Rw0, Rw1))[:, 0]
    f3 = []
    if not r3.max() <= 1.0:
        kk = np.unravel_index(int(np.argmax(r3)), r3.shape)
        f3.append(f"η = r/r_w が float32 の丸めを超えて変わる節点 {int(np.count_nonzero(r3 > 1))} 個 (最大 {r3.max():.3g} 倍、i {kk[0]}, j {kk[1]})")
    if not np.all(excess <= tol_w[:, None]):
        f3.append(f"同じ列の壁より大きく動く節点 {int(np.count_nonzero(excess > tol_w[:, None]))} 個")
    if not (np.all(R0[:, 0] == 0.0) and np.all(R1[:, 0] == 0.0)):
        f3.append("j = 0 が軸 (r = 0) でない")
    # 壁際のセル厚との比 (最近傍なら何層ずれるかの目安)。セル高さの変化の丸め由来の度合い (面積比の極値の説明)
    h0 = np.diff(R0, axis=1)
    hloc = np.minimum(np.c_[h0[:, :1], h0], np.c_[h0, h0[:, -1:]])
    over_h = np.abs(dR) / hloc
    dh = np.diff(R1, axis=1) - h0
    dh_exp = h0 * (dRw / Rw0[:, 0])[:, None]
    dh_tol = (0.5 * (_ulp(R0[:, 1:]) + _ulp(R0[:, :-1]) + _ulp(R1[:, 1:]) + _ulp(R1[:, :-1]))
              + 0.5 * (h0 / Rw0) * (_ulp(Rw0) + _ulp(Rw1)))       # 4 座標と壁の移動の丸め (記録用)
    put("C3_radial_scaling_of_wall_shift", f3, {
        "eta_resid_over_f32_max": float(r3.max()), "abs_dr_minus_abs_drw_max_m": float(excess.max()),
        "wall_shift_max_m": float(np.abs(dRw).max()), "wall_shift_max_x_rt": float(X1[int(np.argmax(np.abs(dRw))), -1] / S),
        "wall_rel_shift_max": float(np.abs(dRw / Rw0[:, 0]).max()), "n_columns_wall_moved": int(np.count_nonzero(dRw)),
        "wall_shift_by_x_band_rt": _bands(X1[:, -1] / S, np.abs(dRw), X_BANDS_RT),
        "disp_over_local_cell_height_max": float(over_h.max()),
        "n_nodes_disp_over_half_cell": int(np.count_nonzero(over_h > 0.5)),
        "cell_height_change_minus_scaling_over_f32_max": float((np.abs(dh - dh_exp) / dh_tol).max()),
        "cell_height_rel_change_max": float(np.abs(dh / h0).max())})
    # --- C4: 壁の移動が設計壁の変化で説明される ---
    f4, d4 = [], {"v4_preestimate": V4_DESIGN_WALL_CHANGE}
    info_d = json.loads((dst / "prepare_info.json").read_text())
    try:
        ev_d = TM.wall_evidence(dst)
    except (OSError, ValueError, KeyError) as e:          # 構造格子でない等 — 不成立として記録する (例外で記録を失わない)
        ev_d = {"status": "error", "reason": f"{type(e).__name__}: {e}"}
    d4["dst_wall_evidence"] = ev_d
    if ev_d.get("status") != "consistent":
        f4.append(f"腕 M の壁節点が腕 M の当てはめ後 spline と一致しない ({ev_d.get('status')})")
    try:
        info_b = json.loads((armB_ref / "prepare_info.json").read_text())
        ev_b = TM.wall_evidence(armB_ref)
        d4["armB_wall_evidence"] = ev_b
        if ev_b.get("status") != "consistent":
            f4.append(f"腕 B の参照の壁節点が腕 B の spline と一致しない ({ev_b.get('status')})")
        from throat_mono_judge import mono_r2_matches
        if not mono_r2_matches((info_b.get("wall_fit") or {}).get("mono_r2"), MONO_R2):
            f4.append(f"腕 B の参照の mono_r2 {(info_b.get('wall_fit') or {}).get('mono_r2')!r} が {MONO_R2} でない")
        if info_b.get("moc") is not None and {k: info_b["moc"].get(k) for k in ("axis_limit", "corrector")} != \
                {"axis_limit": "legacy", "corrector": "fixed2"}:
            f4.append(f"腕 B の参照の moc が legacy・fixed2 でない ({info_b['moc'].get('axis_limit')}, {info_b['moc'].get('corrector')})")
        spl_d, spl_b = TM.wall_spline(info_d), TM.wall_spline(info_b)
        cb, _ = _coords(armB_ref / "nozzle.h5")
        Xb, Rb = cb[:, 0].reshape(ni, nj), cb[:, 1].reshape(ni, nj)
        xw, rw, rwb, xwb = X1[:, -1], R1[:, -1], Rb[:, -1], Xb[:, -1]
        tol_xb = _ulp(np.maximum(np.abs(xw), np.abs(xwb))) + abs(_x_last_design(dst) - _x_last_design(armB_ref))
        if not np.all(np.abs(xw - xwb) <= tol_xb):
            f4.append("腕 M と腕 B の参照の壁節点の x が float32 の丸め以内で一致しない")
        lo = max(spl_d.t[0], spl_b.t[0]) * S
        hi = min(spl_d.t[-1], spl_b.t[-1]) * S
        inr = (xw >= lo) & (xw <= hi)
        d_an = S * (spl_d(xw[inr] / S) - spl_b(xw[inr] / S))
        d_meas = rw[inr] - rwb[inr]
        out_tol = 2.0 * _ulp(np.maximum(rw, rwb))[~inr]
        out_diff = np.abs(rw - rwb)[~inr]
        if out_diff.size and not np.all(out_diff <= out_tol):
            f4.append(f"両 spline の範囲外で腕 M と腕 B の壁が違う (最大 {out_diff.max():.3e} m)")
        kk = int(np.argmax(np.abs(d_an)))
        d4.update(analytic_M_minus_B_max_m=float(np.abs(d_an).max()), analytic_M_minus_B_max_um=float(np.abs(d_an).max() * 1e6),
                  analytic_M_minus_B_x_rt=float(xw[inr][kk] / S), measured_M_minus_B_max_m=float(np.abs(d_meas).max()),
                  measured_minus_analytic_max_m=float(np.abs(d_meas - d_an).max()),
                  outside_spline_range_max_diff_m=(float(out_diff.max()) if out_diff.size else 0.0),
                  B_minus_IC_wall_max_m=float(np.abs(rwb - R0[:, -1]).max()))
        mb = json.loads((armB_ref / "IC_MAP.json").read_text())
        lim_b = (((mb.get("checks") or {}).get("displacement") or {}).get("detail") or {}).get("limit_m")
        d4.update(armB_ic_map_VERDICT=mb.get("VERDICT"), armB_ic_map_limit_m=lim_b, armB_ic_map_src_res=mb.get("src_res"),
                  armB_ic_map_mode=mb.get("mode"))
        if mb.get("VERDICT") != "OK" or mb.get("mode") != "index":
            f4.append(f"腕 B の参照の IC 写像が OK・index でない ({mb.get('VERDICT')!r}, {mb.get('mode')!r})")
        if Path(str(mb.get("src_res"))).name != src_res.name or Path(str(mb.get("src_res"))).parent.name != src_run.name:
            f4.append(f"腕 B の参照の IC が同じ保存場でない ({mb.get('src_res')} / {src_res})")
        if not (isinstance(lim_b, (int, float)) and np.abs(rwb - R0[:, -1]).max() <= float(lim_b)):
            f4.append(f"腕 B の格子と IC 格子の壁の差が腕 B の写像の上限 {lim_b!r} を超える (または上限の記録が無い)")
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
        f4.append(f"腕 B の参照 {armB_ref} を読めない: {type(e).__name__}: {e}")
    put("C4_wall_shift_explained_by_design", f4, d4)
    # --- C5: メッシュ品質 ---
    q_d, q_s = _quality(dst / "MESH_QUALITY.txt"), _quality(src_run / "MESH_QUALITY.txt")
    f5 = [] if (q_d["verdict"] or "").startswith("PASS") else [f"腕 M のメッシュ品質が PASS でない ({q_d['verdict']!r})"]
    put("C5_mesh_quality", f5, {"dst": q_d, "src": q_s})
    # --- 上限と結論 ---
    tol_x_col = tol_x.max(axis=1)
    limit = float(np.max(np.sqrt(tol_x_col ** 2 + (np.abs(dRw) + tol_w) ** 2)))
    ok = all(v["ok"] for v in cond.values())
    rec.update(coord_dtype={"src": dts, "dst": dtd}, ni=ni, nj=nj, scale_m=S, limit_m=limit, limit_um=limit * 1e6,
               limit_reason=("C2 (列の x は float32 の丸め + 設計の出口 x の差以内で不変) と C3 (半径方向は列ごとの壁の移動を η = r/r_w を"
                             "保って縮尺したもの) が成り立つときの列ごとの移動の上界 max_i √(tol_x,i² + (|Δr_w,i| + 2·ulp(r_w,i))²)。"
                             "壁の移動は C4 で 2 本の設計壁の差 (と腕 B の写像で受けた差) に一致することを確かめた。数値を先に決めて広げたものではない"),
               VERDICT=("OK" if ok else "REFUSED"),
               failures=[f"[{k}] {x}" for k, v in cond.items() for x in v["failures"]])
    return rec


def _map_ic(prep_dir: Path, info: dict, src: Path, ic_run: Path, armB_ref: Path, dry: bool) -> dict:
    """IC の検査 (inspect_ic) → 成立なら、その上限で ic_index_map --mode index → prepare_info.json・IC_MAP.json・IC_INSPECTION.json。"""
    from ic_index_map import _sha_file
    try:
        ins = inspect_ic(src, ic_run, prep_dir, armB_ref)
    except Exception as e:  # noqa: BLE001 — 検査できないことも不成立として記録する
        ins = {"plan": PLAN, "tool": "moc_v5_euler.inspect_ic", "VERDICT": "REFUSED", "conditions": {},
               "failures": [f"検査を完了できない: {type(e).__name__}: {e}"], "limit_m": None}
    (prep_dir / "IC_INSPECTION.json").write_text(json.dumps(ins, indent=1, ensure_ascii=False, default=float))
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / f"_band_ab/moc_v5_ic_inspection{'_dry' if dry else ''}.json").write_text(json.dumps(ins, indent=1, ensure_ascii=False, default=float))
    if ins["VERDICT"] != "OK":
        raise RuntimeError("IC の検査が不成立 — 番号写像を使わない:\n  " + "\n  ".join(ins["failures"]))
    cmd = [sys.executable, C / "ic_index_map.py", src, prep_dir / "nozzle.h5", "--mode", "index", "--src-mesh", ic_run / "nozzle.h5",
           "--src-run", ic_run, "--dst-run", prep_dir, "--max-disp-m", repr(ins["limit_m"]), "--record", prep_dir / "IC_MAP.json"]
    if dry:
        cmd.append("--no-species-resolve")
    r = TM._run_tool(cmd, prep_dir / "ic_index_map.log")
    rec = json.loads((prep_dir / "IC_MAP.json").read_text()) if (prep_dir / "IC_MAP.json").exists() else {}
    if r.returncode != 0 or rec.get("VERDICT") != "OK":
        raise RuntimeError(f"ic_index_map (index) が失敗 (rc {r.returncode}, VERDICT {rec.get('VERDICT')}):\n{(r.stdout + r.stderr)[-3000:]}")
    dr = rec["checks"]["displacement"]["detail"]
    if dr.get("limit_m") != ins["limit_m"]:
        raise RuntimeError(f"IC 写像の上限 {dr.get('limit_m')!r} が検査の上限 {ins['limit_m']!r} と違う")
    nv = rec["nearest_vs_index"]
    ic = {"run": str(ic_run), "res": src.name, "tool": "ic_index_map", "mode": "index", "log": "ic_index_map.log", "record": "IC_MAP.json",
          "VERDICT": rec["VERDICT"], "species_resolved": rec["species_resolved"], "transferred": rec["transferred"],
          "displacement_max_um": dr["max_um"], "n_moved": dr["n_moved"], "limit_m": dr["limit_m"],
          "inspection": {"record": "IC_INSPECTION.json", "VERDICT": ins["VERDICT"], "limit_m": ins["limit_m"],
                         "limit_reason": ins["limit_reason"]},
          "nearest_mismatch": nv["n_mismatch"], "nearest_mismatch_dj": nv["dj_counts"],
          "dst_sha256_after": rec["dst_sha256_after"], "dst_mesh_digest": rec["dst_mesh_digest"],
          "tool_last_line": (r.stdout.strip().splitlines() or [""])[-1]}
    ev = TM.wall_evidence(prep_dir)
    if ev.get("status") != "consistent":
        raise RuntimeError(f"壁の証拠 (M4) が不一致: {json.dumps(ev, ensure_ascii=False)}")
    if _sha_file(prep_dir / "nozzle.h5") != rec["dst_sha256_after"]:
        raise RuntimeError("写像の後に nozzle.h5 が変わった")
    info.update(stages="soft", wall_arm=ARM_M, ic=ic, wall_evidence_prep=ev, plan=PLAN)
    info.pop("DRY_NO_IC", None)
    if dry:
        info["DRY"] = True
    else:
        info.pop("DRY", None)
    (prep_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(ARM_M, "| IC index", src.name, "VERDICT", ic["VERDICT"], ("(dry)" if dry else ""), f"| 上限 {ins['limit_um']:.4f} µm",
          f"| 移動 最大 {dr['max_um']:.4f} µm・動いた節点 {dr['n_moved']} | 最近傍の食い違い {nv['n_mismatch']} (j の差 {nv['dj_counts']})",
          "| M4", ev["status"], "|", (prep_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
    return info


def _isen_record(isen_dir: Path, mesh_digest_m: str, dry: bool) -> dict:
    """等エントロピー IC の prep (腕 M の prep を IC 前に複製したもの) の IC 記録。格子が腕 M と同じことを確かめる。"""
    from ic_index_map import _sha_file, mesh_digest   # (solver_density_cuda/tools を sys.path に入れる)
    import forge_species as fsp
    md = mesh_digest(isen_dir / "nozzle.h5")
    if md != mesh_digest_m:
        raise RuntimeError(f"等エントロピー IC の prep の格子が腕 M と同じでない ({md[:16]} / {mesh_digest_m[:16]})")
    st = fsp.source_species_state(str(isen_dir / "nozzle.h5"))
    ok = (st["state"] == "verified") and not dry
    rec = {"plan": PLAN, "tool": "paste_isentropic_ic (runner_axismach.prepare)", "mode": "isentropic",
           "VERDICT": ("OK" if ok else ("DRY (化学種の属性なし)" if dry else f"化学種の属性が検証されていない ({st['state']}: {st['why']})")),
           "species_state": st["state"], "species_why": st["why"], "grid_copied_from": "腕 M の prep (IC を入れる前)",
           "dst_mesh_digest": md, "dst_sha256_after": _sha_file(isen_dir / "nozzle.h5")}
    return rec


def prep(prep_dir: Path, isen_dir: Path | None, armB_ref: Path, dry: bool = False) -> dict:
    if dry:
        # 乾式確認: forge を起動しない (--resolve-species も含めて)。runner の _ENV は import 時に環境を写すので import より前に設定する
        os.environ["FORGE_BIN"] = str(Path("/nonexistent/forge-dry-run-guard"))
        os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
    if isen_dir is not None and Path(isen_dir).exists():
        raise SystemExit(f"{isen_dir} が既にある (消してから作る)")
    if not (armB_ref / "prepare_info.json").is_file() or not (armB_ref / "nozzle.h5").is_file():
        raise SystemExit(f"腕 B の参照 {armB_ref} に prepare_info.json・nozzle.h5 が無い (--armB-ref)")
    RA = TM._load_problem_with_runs()
    ic_run = TM.RUNS / TM.IC_RUN
    src = TM._last_res(ic_run)
    if src.name != IC_RES:
        raise RuntimeError(f"IC run の最後の res が {src.name} (期待 {IC_RES})")
    info = RA.prepare(C / PROBLEM_M, prep_dir, nsteps=TM.NSTEPS, ic_from=None, cfl_main=TM.CFL_MAIN, implicit_relax=TM.RELAX)
    from throat_mono_judge import mono_r2_matches
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf or not mono_r2_matches(wf["mono_r2"], MONO_R2):
        raise RuntimeError(f"腕 M の壁の mono_r2 {wf.get('mono_r2')!r} が {MONO_R2} と完全一致しない")
    moc = info.get("moc") or {}
    got = {k: moc.get(k) for k in MOC_EXPECT}
    gate = moc.get("gate") or {}
    if got != MOC_EXPECT or gate.get("applicable") is not True or gate.get("pass") is not True:
        # 古い design パッケージはキーを読まない (moc が無い) — 黙って腕 B と同じ壁で回さない
        raise RuntimeError(f"prepare_info の moc が {MOC_EXPECT}・ゲート合格でない: {got}, gate {gate}")
    if isen_dir is not None:
        isen_dir = Path(isen_dir).resolve()
        shutil.copytree(prep_dir, isen_dir)          # IC を入れる前 (= RA.prepare の等エントロピー IC) の腕 M の prep を複製
    info = _map_ic(prep_dir, info, src, ic_run, armB_ref, dry)
    if isen_dir is not None:
        rec = _isen_record(isen_dir, info["ic"]["dst_mesh_digest"], dry)
        (isen_dir / "IC_MAP.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))
        ii = json.loads((isen_dir / "prepare_info.json").read_text())
        ev = TM.wall_evidence(isen_dir)
        if ev.get("status") != "consistent":
            raise RuntimeError(f"等エントロピー IC の prep の壁の証拠 (M4) が不一致: {ev}")
        ii.update(stages="soft", wall_arm=ARM_ISEN, ic=rec, wall_evidence_prep=ev, plan=PLAN)
        if dry:
            ii["DRY"] = True
        (isen_dir / "prepare_info.json").write_text(json.dumps(ii, indent=1, default=str))
        print(ARM_ISEN, "| IC isentropic VERDICT", rec["VERDICT"], "| species", rec["species_state"], "| M4", ev["status"],
              "| 格子 = 腕 M", rec["dst_mesh_digest"][:16])
    return info


def run(rd: Path) -> int:
    from forge_design.evaluate import runner_axismach as RA
    info = json.loads((rd / "prepare_info.json").read_text())
    if info.get("DRY") or info.get("DRY_NO_IC"):
        raise SystemExit(f"{rd} は乾式確認 (--dry) の prep から作られている — 回さない")
    if info.get("wall_arm") not in (ARM_M, ARM_ISEN):
        raise SystemExit(f"{rd}: wall_arm {info.get('wall_arm')!r} はこのスクリプトの腕 ({ARM_M} / {ARM_ISEN}) でない")
    if (info.get("ic") or {}).get("VERDICT") != "OK":
        raise SystemExit(f"{rd}: IC の VERDICT が OK でない ({(info.get('ic') or {}).get('VERDICT')!r}) — 回さない")
    gate = (info.get("moc") or {}).get("gate") or {}
    if gate.get("pass") is not True:
        raise SystemExit(f"{rd}: MOC のゲートが合格でない ({gate}) — 回さない")
    rc = RA.run_staged(rd, cfl_main=TM.CFL_MAIN, mid_stage=False, stages="soft")
    print(f"forge exit={rc}")
    return rc


def main(argv) -> int:
    if argv and argv[0] == "prep":
        ap = argparse.ArgumentParser(prog="moc_v5_euler.py prep")
        ap.add_argument("prep_dir")
        ap.add_argument("--with-isen", help="IC 依存の確認用の等エントロピー IC の prep (腕 M の prep を IC 前に複製)")
        ap.add_argument("--armB-ref", default=str(C / ARMB_REF), help=f"腕 B の run (既定 {ARMB_REF}; C4 の照合に使う)")
        ap.add_argument("--dry", action="store_true", help="乾式確認 (forge を起動しない)")
        a = ap.parse_args(argv[1:])
        prep(Path(a.prep_dir).resolve(), (Path(a.with_isen) if a.with_isen else None), Path(a.armB_ref).resolve(), dry=a.dry)
        return 0
    if len(argv) == 2 and argv[0] == "run":
        return run(Path(argv[1]).resolve())
    if len(argv) >= 2 and argv[0] == "verify-prep":
        bad = TM.verify_prep(Path(argv[1]).resolve(), [Path(x).resolve() for x in argv[2:]])
        print("VERIFY-PREP: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
