"""⑤ SERN の δ* 抽出 (S5 一発補正の入力)。

node 構造格子 run から、上バンドの各 x station の列 (j 方向) に沿って ρu 分布を取り、
ランプ (j = nj_top−1 側) とカウル内面 (j = 0 側、x ≤ L_cowl) の排除厚さ
    δ* = ∫ (1 − ρu/(ρu)_e) dn,  縁 = 壁から見て最初に edge_frac·max(ρu) に達する点 (探索窓は列の一部)
を計算する。列は鉛直なので dn = dy·cosθ_w で法線長に直す。cell run は非対応 (VALUE 長 ≠ 節点数)。
戻り値は物理単位 [m] の (x, δ*) 表 (runner_sern.prepare の wall_offset にそのまま渡せる)。
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


def _structured_upper(run_dir: Path):
    info = json.loads((run_dir / "prepare_info.json").read_text())
    m = info["mesh"]; ni, njt, njb, ite = m["ni"], m["nj_top"], m["nj_bot"], m["i_te"]
    with h5py.File(run_dir / "sern.h5") as f:
        cc = f["CELLS/centCoords"][:].reshape(-1, 3)
    N_low = ni * njb; N_up = ni * (njt - 1)
    if info.get("discretization", "cell") == "node":
        dup = m.get("dup_stations", list(range(ite)))
        dup_idx = {int(i): k for k, i in enumerate(dup)}

        def up(i, j):
            if j == 0:
                return N_low + N_up + dup_idx[i] if i in dup_idx else i * njb + (njb - 1)
            return N_low + i * (njt - 1) + (j - 1)
        idx = np.array([[up(i, j) for j in range(njt)] for i in range(ni)])
    else:
        # cell: mesh_sern の quad 生成順 = i ごとに下バンド (njb−1) → 上バンド (njt−1)
        per_i = (njb - 1) + (njt - 1)
        idx = np.array([[i * per_i + (njb - 1) + j for j in range(njt - 1)] for i in range(ni - 1)])
    return info, idx, cc


def deltastar_sern(run_dir, res_h5=None, edge_frac: float = 0.99, search_frac: float = 0.5, min_pts: int = 4) -> dict:
    run_dir = Path(run_dir)
    if res_h5 is None:
        res_h5 = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    info, idx, cc = _structured_upper(run_dir)
    with h5py.File(res_h5) as f:
        ro = f["VALUE/ro"][:]; rux = f["VALUE/roUx"][:]
    if len(ro) != len(cc):
        raise ValueError("VALUE 長 ≠ CELLS 数")
    H = info["H_m"]; ni, njt = idx.shape; ite = info["mesh"]["i_te"]
    # cell モードは壁の座標が列に含まれない (最外セル中心まで)。壁面までの残りは第一セル半分幅の一様欠損として近似:
    # n を壁からの距離に直すため、列端の値と壁の差を n_wall_gap として加える
    X = cc[idx, 0]; Y = cc[idx, 1]; Q = rux[idx]
    ramp = np.loadtxt(run_dir / "ramp_contour.csv", delimiter=",", skiprows=1)
    cowl = np.loadtxt(run_dir / "cowl_contour.csv", delimiter=",", skiprows=1)
    is_cell = info.get("discretization", "cell") != "node"
    out = {"ramp": [], "cowl": []}
    for i in range(1, ni):
        x = float(X[i, 0])
        # --- ramp: 壁 (j=njt-1 側) から下へ。cell では壁 y を輪郭から取り、壁面値 q=0 を先頭に足す
        if is_cell:
            yw = float(np.interp(x, ramp[:, 0], ramp[:, 1])); th = np.arctan(np.gradient(ramp[:, 1], ramp[:, 0])[np.searchsorted(ramp[:, 0], x) - 1])
            n = np.concatenate([[0.0], (yw - Y[i, ::-1]) * np.cos(th)]); q = np.concatenate([[0.0], Q[i, ::-1]])
        else:
            yw = Y[i, -1]; th = np.arctan2(Y[i, -1] - Y[i - 1, -1], X[i, -1] - X[i - 1, -1])
            n = (yw - Y[i, ::-1]) * np.cos(th); q = Q[i, ::-1]
        out["ramp"].append((x, _dstar(n, q, edge_frac, search_frac, min_pts)))
        # --- cowl 内面: 壁 (j=0 側) から上へ (x ≤ L_cowl のみ)
        if (i <= ite - (1 if is_cell else 0)) and x >= 0.0:
            if is_cell:
                yw = float(np.interp(x, cowl[:, 0], cowl[:, 1])); th = np.arctan(-(cowl[-1, 1]) / max(cowl[-1, 0], 1e-30))
                n = np.concatenate([[0.0], (Y[i, :] - yw) * np.cos(th)]); q = np.concatenate([[0.0], Q[i, :]])
            else:
                yw = Y[i, 0]; th = np.arctan2(Y[i, 0] - Y[i - 1, 0], X[i, 0] - X[i - 1, 0])
                n = (Y[i, :] - yw) * np.cos(th); q = Q[i, :]
            out["cowl"].append((x, _dstar(n, q, edge_frac, search_frac, min_pts)))
    res = {k: np.asarray(v) for k, v in out.items()}
    for k in res:
        good = np.isfinite(res[k][:, 1])
        res[k] = res[k][good]
    res["H_m"] = H; res["res_file"] = str(res_h5)
    return res


def _dstar(n, q, edge_frac, search_frac, min_pts) -> float:
    n = np.asarray(n, dtype=float); q = np.asarray(q, dtype=float)
    k = max(int(len(n) * search_frac), min_pts + 1)
    qs = q[:k]
    qe = float(np.max(qs))
    if qe <= 0:
        return np.nan
    ie = int(np.argmax(qs >= edge_frac * qe))
    if ie < min_pts:
        return np.nan
    m = slice(0, ie + 1)
    return float(np.trapezoid(1.0 - q[m] / qe, n[m]))


def smooth_table(tbl, win: int = 9):
    """(x, δ*) 表の移動平均 (端は縮小窓)。"""
    tbl = np.asarray(tbl, dtype=float)
    if len(tbl) < 3:
        return tbl
    d = tbl[:, 1].copy(); out = d.copy()
    h = win // 2
    for i in range(len(d)):
        a, b = max(0, i - h), min(len(d), i + h + 1)
        out[i] = np.mean(d[a:b])
    return np.column_stack([tbl[:, 0], out])


def write_offset_json(res: dict, path, win: int = 9) -> dict:
    wo = {"ramp": smooth_table(res["ramp"], win).tolist(), "cowl": smooth_table(res["cowl"], win).tolist()}
    Path(path).write_text(json.dumps(wo, indent=1))
    return wo


# --- Euler 基準の排除厚さ (非一様な非粘性コアを持つ SERN 用。単独抽出は膨張扇の ρu 変化を欠損と誤認する) ----
def deltastar_sern_vs_euler(run_ns, run_euler, res_ns=None, res_euler=None, thresh: float = 0.02, min_pts: int = 3) -> dict:
    """δ* = ∫ (1 − (ρu)_NS/(ρu)_Euler) dn を壁から、欠損がピーク後に初めて thresh を下回る点まで積分。
    Euler 場は別メッシュでもよい (NS 列点へ線形補間)。cell/node どちらも可。戻りは [m] の (x, δ*) 表。"""
    from scipy.interpolate import LinearNDInterpolator
    run_ns, run_euler = Path(run_ns), Path(run_euler)
    if res_ns is None:
        res_ns = sorted(run_ns.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    if res_euler is None:
        res_euler = sorted(run_euler.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    info, idx, cc = _structured_upper(run_ns)
    with h5py.File(res_ns) as f:
        Q = f["VALUE/roUx"][:][idx]
    with h5py.File(run_euler / "sern.h5") as f:
        ce = f["CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res_euler) as f:
        qe_all = f["VALUE/roUx"][:]
    itp = LinearNDInterpolator(ce[:, :2], qe_all)
    X = cc[idx, 0]; Y = cc[idx, 1]
    QE = itp(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
    H = info["H_m"]; ni, njt = idx.shape; ite = info["mesh"]["i_te"]
    ramp = np.loadtxt(run_ns / "ramp_contour.csv", delimiter=",", skiprows=1)
    cowl = np.loadtxt(run_ns / "cowl_contour.csv", delimiter=",", skiprows=1)
    is_cell = info.get("discretization", "cell") != "node"
    out = {"ramp": [], "cowl": []}

    def integ(n, d):
        good = np.isfinite(d)
        n, d = n[good], d[good]
        if len(d) < min_pts + 1:
            return np.nan
        ip = int(np.argmax(d[: max(min_pts, len(d) // 3)]))
        after = np.where(d[ip:] < thresh)[0]
        ie = ip + (int(after[0]) if len(after) else len(d) - ip - 1)
        if ie < min_pts:
            return np.nan
        return float(np.trapezoid(np.clip(d[: ie + 1], 0.0, None), n[: ie + 1]))

    for i in range(1, ni):
        x = float(X[i, 0])
        d = 1.0 - Q[i, ::-1] / np.maximum(QE[i, ::-1], 1e-30)
        if is_cell:
            yw = float(np.interp(x, ramp[:, 0], ramp[:, 1])); th = float(np.arctan(np.gradient(ramp[:, 1], ramp[:, 0])[min(np.searchsorted(ramp[:, 0], x), len(ramp) - 1)]))
            n = np.concatenate([[0.0], (yw - Y[i, ::-1]) * np.cos(th)]); d = np.concatenate([[1.0], d])
        else:
            yw = Y[i, -1]; th = np.arctan2(Y[i, -1] - Y[i - 1, -1], X[i, -1] - X[i - 1, -1])
            n = (yw - Y[i, ::-1]) * np.cos(th)
        out["ramp"].append((x, integ(n, d)))
        if (i <= ite - (1 if is_cell else 0)) and x >= 0.0:
            d = 1.0 - Q[i, :] / np.maximum(QE[i, :], 1e-30)
            if is_cell:
                yw = float(np.interp(x, cowl[:, 0], cowl[:, 1])); th = float(np.arctan(-(cowl[-1, 1]) / max(cowl[-1, 0], 1e-30)))
                n = np.concatenate([[0.0], (Y[i, :] - yw) * np.cos(th)]); d = np.concatenate([[1.0], d])
            else:
                yw = Y[i, 0]; th = np.arctan2(Y[i, 0] - Y[i - 1, 0], X[i, 0] - X[i - 1, 0])
                n = (Y[i, :] - yw) * np.cos(th)
            out["cowl"].append((x, integ(n, d)))
    res = {k: np.asarray(v) for k, v in out.items()}
    for k in res:
        res[k] = res[k][np.isfinite(res[k][:, 1])]
    res["H_m"] = H; res["res_file"] = str(res_ns); res["euler_ref"] = str(res_euler)
    return res
