r"""RANS 場からの排除厚 $\delta^*(x)$ 抽出 (A12 = axis-Mach 粘性補正)。

計画: plans/active/tooling-nozzle-axismach-viscous-deltastar.md §4。

定義 (圧縮性・**軸対称の質量収支** — 2026-09-01 に平面近似から置換):

$$\delta^*\left(1-\frac{\delta^*}{2r_w}\right) = \int_0^{y_e}
  \left(1 - \frac{\rho u}{(\rho u)_e}\right)\left(1-\frac{y}{r_w}\right) dy
  \;\Rightarrow\; \delta^* = r_w\left(1-\sqrt{1-2I/r_w}\right)$$

壁から $y$ 入った欠損は半径 $r_w-y$ の円環を流れるため $(1-y/r_w)$ の重みが付く。
旧・平面近似 ($\delta^*=\int(1-\rho u/(\rho u)_e)dy$) は欠損重心が壁近くに集中する
ため過大評価は小さい (M6 出口 $\delta/r_w\approx0.14$ で **+1.4 %**、出口 M 換算
+0.05 % — case/45 run_0007 実測)。帳簿の正確さのため軸対称形を正とする。

**教訓 (run_0030 で実測、繰り返し禁止)**:
- 探索窓はフリーストリームに届くまで取る (窓 0.35 r* は下流で BL 縁の手前で切れ、
  「δ* が下流ほど減る」という非物理な結論を出した — 撤回済み)。既定 1.5 r*。
- スロート近傍 (x ≲ 8 r_t) はコア流が半径方向に未一様で「ρu 極大 = BL 縁」の
  判定が破綻する。呼び出し側で相関へフォールバックすること。
"""
from __future__ import annotations

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def deltastar_from_run(mesh_h5, res_h5, wall_xy, scale: float,
                       x_stations, window: float = 1.5, n_prof: int = 400,
                       edge_frac: float = 0.995) -> dict:
    """node RANS run から各 x ステーションの δ* を測る。

    wall_xy: (n,2) 物理壁 [m 単位 or r_t 単位 — scale と整合していれば可]。
    x_stations: r_t 単位。戻り: dict(x, dstar, ue_rel, y_edge, ok) 全て r_t 単位。
    """
    with h5py.File(mesh_h5) as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res_h5) as f:
        ro = f["/VALUE/ro"][:]
        Ux = f["/VALUE/Ux"][:]
    if len(ro) != len(nc):
        raise ValueError("node run でない (VALUE 長 != 節点数)")
    x_n, r_n = nc[:, 0] / scale, nc[:, 1] / scale
    itp = LinearNDInterpolator(np.c_[x_n, r_n], ro * Ux)
    wx, wr = np.asarray(wall_xy)[:, 0], np.asarray(wall_xy)[:, 1]
    out = {"x": [], "dstar": [], "y_edge": [], "ok": []}
    for xs in np.atleast_1d(x_stations):
        rw = float(np.interp(xs, wx, wr))
        y = np.linspace(0.0, min(window, 0.95 * rw), n_prof)   # 壁 → 内側
        q = itp(np.full(n_prof, xs), rw - y)
        good = np.isfinite(q)
        if good.sum() < n_prof // 2:
            out["x"].append(float(xs)); out["dstar"].append(np.nan)
            out["y_edge"].append(np.nan); out["ok"].append(False)
            continue
        q = np.where(good, q, 0.0)
        qe = float(np.max(q))
        # 縁 = 壁から見て最初に edge_frac·(ρu)_max へ達する y
        idx = np.argmax(q >= edge_frac * qe)
        y_e = float(y[idx])
        m = y <= y_e
        # 軸対称の質量収支: 欠損に円環重み (1 - y/r_w) を掛け、
        # δ*(1 - δ*/(2 r_w)) = I を δ* について解く (docstring 参照)
        d = 1.0 - q[m] / max(qe, 1e-30)
        I = float(np.trapezoid(d * (1.0 - y[m] / rw), y[m]))
        integ = float(rw * (1.0 - np.sqrt(max(1.0 - 2.0 * I / rw, 0.0))))
        out["x"].append(float(xs)); out["dstar"].append(integ)
        out["y_edge"].append(y_e); out["ok"].append(bool(idx > 3))
    return {k: np.asarray(v) for k, v in out.items()}


# --- 固定 Euler 基準のコア整合質量欠損抽出 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) ----
def _load_structured(run_dir):
    """node 構造格子 run (index = i*nj + j) を (info, x[ni,nj], r[ni,nj], q=ρUx[ni,nj]) で返す (r_t 単位)。"""
    import json
    from pathlib import Path
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
    with h5py.File(res) as f:
        q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
    if len(q) != ni * nj or len(nc) != ni * nj:
        raise ValueError(f"{run_dir}: node 構造格子 run でない (VALUE {len(q)}, COORD {len(nc)}, ni*nj {ni*nj})")
    return dict(info=info, S=S, res=res.name,
                x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S, q=q.reshape(ni, nj))


def massflow_ratio(ns_run, euler_run) -> dict:
    """NS/Euler 質量流量比 → 有効スロート半径・質量流量由来の等価スロート補正量 (r_t 単位)。
    列ごとに 2π∫ρUx r dr (台形則)、x∈(−2.5, x_E) の中央値を採用。"""
    out = {}
    for tag, rd in (("ns", ns_run), ("euler", euler_run)):
        d = _load_structured(rd)
        S = d["S"]
        m = 2 * np.pi * np.trapezoid(d["q"] * d["r"] * S * S, d["r"] * S, axis=1)
        xs = d["x"][:, 0]
        sel = (xs > -2.5) & (xs < d["info"]["x_E"])
        out[f"mdot_{tag}"] = float(np.median(m[sel]))
        out[f"mdot_{tag}_spread"] = float(m[sel].std() / np.median(m[sel]))
        out[f"res_{tag}"] = d["res"]
        if tag == "ns":
            r_tW = float(d["info"].get("throat_physical", {}).get("r", np.nan)) if d["info"].get("throat_physical") else float(np.min(d["r"][:, -1]))
    ratio = out["mdot_ns"] / out["mdot_euler"]
    out["mdot_ratio"] = float(ratio)
    out["r_throat_eff"] = float(np.sqrt(ratio))            # 設計 r_t 単位
    out["r_throat_phys"] = r_tW
    out["cd_rel"] = float(ratio / r_tW ** 2)
    out["delta_r_throat_eff"] = float(r_tW - np.sqrt(ratio))
    out["delta_r_throat_given"] = float(r_tW - 1.0)
    return out



def core_matched_deficit(r, q_ns, q_e, rw_e: float, core_frac: float = 0.30,
                         n_axis_skip: int = 1, outer_frac: float = 0.25, refine: int = 4):
    r"""1 断面の**コア整合質量欠損 → 半径方向等価排除厚** (純関数、単体試験対象)。

    r: 軸→壁の半径 (NS 格子, 最後が NS 壁)、q_ns: NS の ρUx、q_e: 同じ r に写した Euler の ρUx
    (Euler 壁 rw_e より外は壁値) — 配列 (r と同長) か callable(r)。
    refine>1 で NS 区間を細分 (q_ns は区分線形、q_e は callable なら細分点で評価、配列なら線形補間)。
    壁集中メッシュはコアに節点が ~5 個しか無いので、積分と α を安定させるために使う。
    戻り: dict または None (コア点数不足)。"""
    r0 = np.asarray(r, dtype=float); q0 = np.asarray(q_ns, dtype=float)
    if refine > 1:
        t = np.linspace(0.0, 1.0, refine, endpoint=False)
        r = np.concatenate([(r0[:-1, None] + (r0[1:] - r0[:-1])[:, None] * t[None, :]).ravel(), [r0[-1]]])
        qN = np.interp(r, r0, q0)
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.interp(r, r0, np.asarray(q_e, dtype=float))
    else:
        r, qN = r0, q0
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.asarray(q_e, dtype=float)
    rwN = float(r[-1])
    mc = r <= core_frac * rw_e
    if mc.sum() < n_axis_skip + 3:
        return None
    num = np.trapezoid(qN[mc] * r[mc], r[mc]); den = np.trapezoid(qE[mc] * r[mc], r[mc])
    alpha = float(num / den) if den > 0 else np.nan
    qref = alpha * qE
    rel = qN[mc] / np.maximum(qref[mc], 1e-30) - 1.0
    core_rms = float(np.sqrt(np.mean(rel ** 2)))
    core_rms_noaxis = float(np.sqrt(np.mean(rel[n_axis_skip:] ** 2)))
    core_maxdev = float(np.max(np.abs(rel[n_axis_skip:])))
    integrand = (qref - qN) * r
    D = float(2 * np.pi * np.trapezoid(integrand, r))
    # 壁からの累積 F(y) = 2π∫_{rw-y}^{rw} q_ref r' dr' (y = 壁からの距離, 単調増加)
    seg = 0.5 * (qref[1:] * r[1:] + qref[:-1] * r[:-1]) * np.diff(r)
    F_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    r_desc = r[::-1]
    reason = []
    if D < 0:
        q_scale = max(float(np.mean(qref[mc])), 1e-30)
        delta_r = D / (2 * np.pi * rwN * q_scale)                          # 線形化 (負値, 診断用)
        reason.append("negative_deficit")
    elif D > F_from_wall[-1]:
        delta_r = float("nan"); reason.append("no_root")
    else:
        r_eff = float(np.interp(D, F_from_wall, r_desc))
        delta_r = rwN - r_eff
    segD = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)
    G_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(segD[::-1])])
    y_desc = rwN - r_desc
    G_outer = float(np.interp(outer_frac * rwN, y_desc, G_from_wall))
    outer_share = G_outer / D if D > 0 else float("nan")
    return dict(r_wall_euler=float(rw_e), r_wall_ns=rwN, delta_in=rwN - float(rw_e), alpha=alpha,
                core_rms=core_rms, core_rms_noaxis=core_rms_noaxis, core_maxdev=core_maxdev,
                mass_deficit=D, delta_r=float(delta_r), outer_share=outer_share, reason=reason)



def band_local_deficit(r, q_ns, q_e, rw_e: float, delta_in: float | None = None,
                       k_band: float = 1.5, band_min_frac: float = 0.02, band_width: float = 0.5,
                       n_axis_skip: int = 1, core_frac: float = 0.30, refine: int = 4, _retry: int = 0,
                       slope_max: float = 0.01, band_growth: float = 1.25, max_retry: int = 14,
                       fit_from: float = 1.0, y_b_fixed: float | None = None):
    r"""1 断面の**帯局所参照・質量欠損 → 半径方向等価排除厚** (生産抽出、純関数)。

    plan §4.3 の修正 (2026-09-04 V1 pass1 の知見): 内側 30 % コアの単一倍率 α では、上流の壁 δ 誤差が
    特性線で下流の軸へ運ぶ波 (コア形状差) を欠損に取り込む。そこで参照 $q_{ref}$ は
    **境界層のすぐ外側の帯** $y\in[y_b, (1+w)y_b]$ で比 $q_{NS}/q_E$ を $y$ の 1 次で LSQ フィットし、
    境界層域 $y<y_b$ へ延長して作る。$y_b$ は $\max(k\,\delta_{in},\ f_{min} r_w)$ (既定 k=1.5 = 境界層の中) から始め、
    帯内の比の変化が slope_max (1 %) を下回るまで band_growth (1.25) 倍ずつ広げる (適応帯: 4 つの異なる壁の NS 場で
    δ_r のばらつき 3 %、固定 k=4/w=1 では 6.5 % — case/45 run_0004/0019/0020/0021)。
    欠損 $D = 2\pi\int_{y<y_b}(q_{ref}-q_{NS})r\,dr$、$\delta_r = r_w - r_{eff}$。
    内側 30 % の α・コア RMS は診断として併記する。q_e は配列または callable(r)。"""
    r0 = np.asarray(r, dtype=float); q0 = np.asarray(q_ns, dtype=float)
    if refine > 1:
        t = np.linspace(0.0, 1.0, refine, endpoint=False)
        rr = np.concatenate([(r0[:-1, None] + (r0[1:] - r0[:-1])[:, None] * t[None, :]).ravel(), [r0[-1]]])
        qN = np.interp(rr, r0, q0)
        qE = np.asarray(q_e(rr), dtype=float) if callable(q_e) else np.interp(rr, r0, np.asarray(q_e, dtype=float))
    else:
        rr, qN = r0, q0
        qE = np.asarray(q_e(rr), dtype=float) if callable(q_e) else np.asarray(q_e, dtype=float)
    rwN = float(rr[-1]); y = rwN - rr
    d_in = float(rwN - rw_e) if delta_in is None else float(delta_in)
    if y_b_fixed is not None:                      # 2 パス目: 隣接列で平滑化した帯位置を固定して使う (帯のトグルを防ぐ)
        y_b0 = float(y_b_fixed); max_retry = 0
    else:
        y_b0 = max(k_band * max(d_in, 0.0), band_min_frac * rwN) * (band_growth ** _retry)   # 適応探索する「縁」候補
    y_b = fit_from * y_b0                                                                 # 参照フィット帯の下端 = 積分域の上端
    y_t = min((1.0 + band_width) * y_b, 0.85 * rwN)
    if y_b >= 0.8 * rwN:
        return None
    mb = (y >= y_b) & (y <= y_t)
    if mb.sum() < 4:
        return None
    ratio = qN[mb] / np.maximum(qE[mb], 1e-30)
    A = np.c_[np.ones(mb.sum()), y[mb] - y_b]
    coef, *_ = np.linalg.lstsq(A, ratio, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    band_rms = float(np.sqrt(np.mean((ratio - (a + b * (y[mb] - y_b))) ** 2)))
    # 縁候補帯 [y_b0, (1+w) y_b0] での比の変化 (適応判定用; fit_from=1 なら上と同じ帯)
    mb0 = (y >= y_b0) & (y <= min((1.0 + band_width) * y_b0, 0.85 * rwN))
    if mb0.sum() >= 4:
        c0, *_ = np.linalg.lstsq(np.c_[np.ones(mb0.sum()), y[mb0] - y_b0], qN[mb0] / np.maximum(qE[mb0], 1e-30), rcond=None)
        slope0 = float(c0[1]) * (min((1.0 + band_width) * y_b0, 0.85 * rwN) - y_b0)
        a0 = float(c0[0])
    else:
        slope0, a0 = 0.0, a
    # 内側コアの単一倍率 (診断 + 「帯が境界層の深部にある」検出: 縁候補の比がコア比より 5 % 以上低い)
    mc = rr <= core_frac * rw_e
    if mc.sum() >= n_axis_skip + 3:
        num = np.trapezoid(qN[mc] * rr[mc], rr[mc]); den = np.trapezoid(qE[mc] * rr[mc], rr[mc])
        alpha = float(num / den) if den > 0 else float("nan")
        rel = qN[mc] / np.maximum(alpha * qE[mc], 1e-30) - 1.0
        core_rms = float(np.sqrt(np.mean(rel[n_axis_skip:] ** 2)))
    else:
        alpha, core_rms = float("nan"), float("nan")
    if np.isfinite(alpha) and a0 < 0.95 * alpha and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    # 帯が境界層の中にあると比 q_NS/q_E が帯の幅で数 % 以上変わる (コアの波は < 1 %)。
    # 線形フィットは BL の形も吸収してしまい残留欠損では検出できないので、帯内の比の変化量で判定し y_b を倍にする。
    if abs(slope0) > slope_max and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    ratio_fn = a + b * (y - y_b)
    qref = ratio_fn * qE
    mz = y <= y_b                                     # 境界層域 (壁側)
    integrand = (qref - qN) * rr
    D = float(2 * np.pi * np.trapezoid(integrand[mz], rr[mz]))
    # 帯の中の残留欠損 (参照の整合診断): 本来 ~0。大きければ帯が境界層の中にある → y_b を 2 倍して再試行
    # (δ99/δ* は 2〜10 と場所で変わるので k_band 固定では足りない。閾値でなく自己整合で縁を外す)
    D_band = float(2 * np.pi * np.trapezoid(integrand[mb], rr[mb]))
    if D > 0 and abs(D_band) > 0.1 * D and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    seg = 0.5 * (qref[1:] * rr[1:] + qref[:-1] * rr[:-1]) * np.diff(rr)
    F_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    r_desc = rr[::-1]
    reason = []
    if D < 0:
        q_scale = max(float(np.mean(qref[mz])), 1e-30)
        delta_r = D / (2 * np.pi * rwN * q_scale)
        reason.append("negative_deficit")
    else:
        F_zone = float(np.interp(rwN - y_b, r_desc[::-1], F_from_wall[::-1]))  # 帯下端までの容量
        if D > F_zone:
            if _retry < max_retry:
                return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                          n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                          band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
            delta_r = float("nan"); reason.append("no_root")
        else:
            r_eff = float(np.interp(D, F_from_wall, r_desc))
            delta_r = rwN - r_eff
    return dict(r_wall_euler=float(rw_e), r_wall_ns=rwN, delta_in=rwN - float(rw_e), alpha=alpha,
                core_rms=core_rms, core_rms_noaxis=core_rms, core_maxdev=float("nan"),
                mass_deficit=D, delta_r=float(delta_r), outer_share=float("nan"),
                band_y_b=y_b0, band_fit_from=fit_from, band_ratio_a=a, band_slope=slope0, band_rms=band_rms,
                band_deficit_share=(D_band / D if D > 0 else float("nan")), band_retry=_retry,
                reason=reason)


def deltastar_from_core_matched_euler(ns_run, euler_run, core_frac: float = 0.30,
                                      core_frac_sens=(0.25, 0.35), n_axis_skip: int = 1,
                                      outer_frac: float = 0.25, smooth_lam: float = 1e-3,
                                      gate_core_rms: float = 0.01, gate_sens: float = 0.05,
                                      gate_sens_floor: float = 1e-3, refine: int = 4, out_dir=None,
                                      method: str = "band", k_band: float = 1.5, k_band_sens=(1.125, 2.25),
                                      band_min_frac: float = 0.02, band_width: float = 0.5,
                                      gate_band_rms: float = 0.005, slope_max: float = 0.01,
                                      band_growth: float = 1.25, max_retry: int = 14, fit_from: float = 1.0,
                                      band_smooth_x: float = 0.0) -> dict:
    r"""**固定 Euler 基準・コア整合**の半径方向等価排除厚 $\delta_r(x)$ を NS 全列で抽出する。

    定義 (plan §4.2–4.4):

    - 対応は**同じ物理 $r$** (設計 $r_t$ 単位)。NS の有効壁 $r_{ref}=r_{w,NS}-\delta_{in}$ は、
      $\delta_{in}$ を「NS 壁 − Euler 壁」で定義するので $r_{ref}=r_{w,E}$ となり、η 写像は恒等になる。
      Euler 壁の外側 ($r>r_{w,E}$) は Euler 壁値で定数外挿。
    - コア $r \le \eta_c r_{w,E}$ の面積重み流量を一致させる倍率 $\alpha$、$q_{ref}=\alpha q_E$。
    - 符号付き欠損 $D = 2\pi\int_0^{r_{w,NS}}(q_{ref}-q_{NS})\,r\,dr$ (クリップなし)。
      $2\pi\int_{r_{eff}}^{r_{w,NS}} q_{ref}\,r\,dr = D$ の根から $\delta_r = r_{w,NS}-r_{eff}$。
    - ゲート (診断; 値は捨てない): コア相対 RMS、コア範囲感度、$D<0$、根なし、欠損のコア分布。

    戻り値: dict of arrays (x, r_wall_euler, r_wall_ns, delta_in, alpha, core_rms, core_rms_noaxis,
    core_maxdev, mass_deficit, delta_r_raw, delta_r_smooth, delta_r_sens (n,2), ok, reason) と mdot 帳簿。
    `out_dir` を与えると `delta_r_equiv.csv` / `delta_r_equiv_diag.json` を書く。"""
    import json
    from pathlib import Path
    from scipy.interpolate import make_smoothing_spline
    N = _load_structured(ns_run)
    E = _load_structured(euler_run)
    xE = E["x"][:, 0]
    rwE_col = E["r"][:, -1]

    def euler_q_at(x, r):
        """Euler 場を (x, r) へ: 列間線形 × 列内線形 (r>壁は壁値)。"""
        k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
        w = (x - xE[k]) / max(xE[k + 1] - xE[k], 1e-30)
        w = float(np.clip(w, 0.0, 1.0))
        q0 = np.interp(r, E["r"][k], E["q"][k])          # np.interp は範囲外を端値で埋める
        q1 = np.interp(r, E["r"][k + 1], E["q"][k + 1])
        rw = (1 - w) * rwE_col[k] + w * rwE_col[k + 1]
        return (1 - w) * q0 + w * q1, rw

    def extract_one(i, fc, y_b_fixed=None):
        x = float(N["x"][i, 0]); r0 = N["r"][i]; q0 = N["q"][i]
        if x < xE[0] - 1e-9 or x > xE[-1] + 1e-9:
            return None
        _, rwE = euler_q_at(x, r0[:1])
        if method == "band":
            res = band_local_deficit(r0, q0, lambda rr: euler_q_at(x, rr)[0], float(rwE), k_band=fc,
                                     band_min_frac=band_min_frac, band_width=band_width,
                                     n_axis_skip=n_axis_skip, core_frac=core_frac, refine=refine,
                                     slope_max=slope_max, band_growth=band_growth, max_retry=max_retry, fit_from=fit_from,
                                     y_b_fixed=y_b_fixed)
        else:
            res = core_matched_deficit(r0, q0, lambda rr: euler_q_at(x, rr)[0], float(rwE), core_frac=fc,
                                       n_axis_skip=n_axis_skip, outer_frac=outer_frac, refine=refine)
        if res is None:
            return None
        res["x"] = x
        return res

    main_par = k_band if method == "band" else core_frac
    sens_par = list(k_band_sens) if method == "band" else list(core_frac_sens)
    yb_fixed = None
    if method == "band" and band_smooth_x > 0:
        # 1 パス目: 各列の適応 y_b。隣接列で帯が 1.25 倍刻みにトグルすると δ_r が ±1 % ジグザグになる
        # (run_0022 x=90–94 で実測) ので、log y_b を x 方向に移動中央値 (窓 band_smooth_x [r_t]) で平滑化して固定する。
        xs_all = N["x"][:, 0]; yb_all = np.full(len(xs_all), np.nan)
        for i in range(len(xs_all)):
            b0 = extract_one(i, main_par)
            if b0 is not None:
                yb_all[i] = b0.get("band_y_b", np.nan)
        fin = np.isfinite(yb_all)
        lyb = np.log(np.where(fin, yb_all, np.nan))
        sm = np.full_like(lyb, np.nan)
        for i in range(len(xs_all)):
            w = fin & (np.abs(xs_all - xs_all[i]) <= 0.5 * band_smooth_x)
            if w.sum() >= 3:
                sm[i] = np.median(lyb[w])
        yb_fixed = np.exp(sm)
    rows = []
    for i in range(N["x"].shape[0]):
        if yb_fixed is not None and np.isfinite(yb_fixed[i]):
            base = extract_one(i, main_par, y_b_fixed=float(yb_fixed[i]))
            sens = [extract_one(i, main_par, y_b_fixed=float(yb_fixed[i]) * f) for f in (0.8, 1.25)]
        else:
            base = extract_one(i, main_par)
            sens = [extract_one(i, fc) for fc in sens_par]
        if base is None:
            continue
        base["delta_r_sens"] = [s["delta_r"] if s else np.nan for s in sens]
        rows.append(base)
    x = np.array([r_["x"] for r_ in rows])
    draw = np.array([r_["delta_r"] for r_ in rows])
    sens = np.array([r_["delta_r_sens"] for r_ in rows])
    HARD = ("negative_deficit", "no_root", "nan")
    ok = np.ones(len(rows), bool); hard_ok = np.ones(len(rows), bool); reasons = []
    for k, r_ in enumerate(rows):
        rs = list(r_["reason"])
        if method == "band":
            if r_.get("band_rms", 0.0) > gate_band_rms: rs.append("band_fit")
            if np.isfinite(r_.get("band_deficit_share", np.nan)) and abs(r_["band_deficit_share"]) > 0.1: rs.append("band_not_clean")
        elif r_["core_rms_noaxis"] > gate_core_rms: rs.append("core_shape")
        if np.isfinite(draw[k]) and draw[k] > 0 and np.all(np.isfinite(sens[k])):
            if np.max(np.abs(sens[k] - draw[k])) > max(gate_sens * draw[k], gate_sens_floor): rs.append("core_frac_sens")
        if np.isfinite(r_["outer_share"]) and r_["outer_share"] < 0.8: rs.append("deficit_not_near_wall")
        if not np.isfinite(draw[k]): rs.append("nan")
        ok[k] = len(rs) == 0
        hard_ok[k] = not any(h in rs for h in HARD)
        reasons.append(",".join(rs))
    # 平滑化: hard-ok の生値だけを使い、範囲外は端値 (外挿しない)。soft 不合格の値は使う (plan §4.4)。
    dsm = np.full_like(draw, np.nan)
    if hard_ok.sum() > 10:
        xs_, ds_ = x[hard_ok], draw[hard_ok]
        spl = make_smoothing_spline(xs_, ds_, lam=smooth_lam)
        dsm = spl(np.clip(x, xs_[0], xs_[-1]))
    # 壁更新に渡す値: hard-ok なら平滑化値、hard 不合格なら NaN (呼び出し側が前回値を保持する)
    duse = np.where(hard_ok, dsm, np.nan)
    out = dict(x=x, r_wall_euler=np.array([r_["r_wall_euler"] for r_ in rows]),
               r_wall_ns=np.array([r_["r_wall_ns"] for r_ in rows]),
               delta_in=np.array([r_["delta_in"] for r_ in rows]),
               alpha=np.array([r_["alpha"] for r_ in rows]),
               core_rms=np.array([r_["core_rms"] for r_ in rows]),
               core_rms_noaxis=np.array([r_["core_rms_noaxis"] for r_ in rows]),
               core_maxdev=np.array([r_["core_maxdev"] for r_ in rows]),
               mass_deficit=np.array([r_["mass_deficit"] for r_ in rows]),
               outer_share=np.array([r_["outer_share"] for r_ in rows]),
               band_y_b=np.array([r_.get("band_y_b", np.nan) for r_ in rows]),
               band_slope=np.array([r_.get("band_slope", np.nan) for r_ in rows]),
               band_rms=np.array([r_.get("band_rms", np.nan) for r_ in rows]),
               band_deficit_share=np.array([r_.get("band_deficit_share", np.nan) for r_ in rows]),
               band_retry=np.array([r_.get("band_retry", 0) for r_ in rows]),
               delta_r_raw=draw, delta_r_smooth=dsm, delta_r_use=duse, delta_r_sens=sens,
               ok=ok, hard_ok=hard_ok, reason=np.array(reasons),
               settings=dict(method=method, k_band=k_band, k_band_sens=list(k_band_sens), band_min_frac=band_min_frac,
                             band_width=band_width, slope_max=slope_max, band_growth=band_growth, max_retry=max_retry, fit_from=fit_from,
                             band_smooth_x=band_smooth_x,
                             gate_band_rms=gate_band_rms,
                             core_frac=core_frac, core_frac_sens=list(core_frac_sens), n_axis_skip=n_axis_skip,
                             outer_frac=outer_frac, smooth_lam=smooth_lam, gate_core_rms=gate_core_rms,
                             gate_sens=gate_sens, gate_sens_floor=gate_sens_floor, refine=refine,
                             ns_res=N["res"], euler_res=E["res"],
                             ns_run=str(ns_run), euler_run=str(euler_run)))
    out["massflow"] = massflow_ratio(ns_run, euler_run)
    if out_dir is not None:
        od = Path(out_dir); od.mkdir(parents=True, exist_ok=True)
        np.savetxt(od / "delta_r_equiv.csv",
                   np.c_[x, draw, dsm, duse, out["delta_in"], out["alpha"], out["core_rms_noaxis"], sens,
                         out["band_y_b"], out["band_slope"], out["band_rms"], out["band_deficit_share"],
                         ok.astype(int), hard_ok.astype(int)],
                   delimiter=",", comments="",
                   header="x_rt,delta_r_raw,delta_r_smooth,delta_r_use,delta_in,alpha,core_rms_noaxis,delta_r_sens1,delta_r_sens2,band_y_b,band_slope,band_rms,band_deficit_share,ok,hard_ok")
        diag = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in out.items()}
        (od / "delta_r_equiv_diag.json").write_text(json.dumps(diag))
    return out


# --- δ_r(x) の 5 次 P-spline 平滑化 (壁曲率を滑らかにする; 2026-09-04 ユーザ指摘) ----------------------------
def smooth_delta_quintic(x, d, weights=None, knot_spacing: float = 2.0, lam: float = 1.0,
                         x_stretch: float = 4.0, penalty_order: int = 3, positive: bool = True):
    r"""抽出した $\delta_r(x)$ を **5 次 B-spline + 3 階差分ペナルティ (P-spline)** で平滑化する。

    壁は $r_{inv}(x)+\delta_r(x)$ を補間 5 次 B-spline で通すので、δ_r の 2 階微分のノイズがそのまま壁曲率に乗る。
    3 次平滑化スプラインは曲率まで滑らかにしないため、$\int(\delta_r''')^2$ に相当する係数 3 階差分ペナルティで
    **曲率が滑らかな** δ_r を作る (ユーザ提案)。ノットは $\xi=\sqrt{x-x_0+x_s}$ で等間隔 (スロート側で密:
    δ_r は x=0→8 で 20 倍に育つ)。weights=0 の点 (hard 不合格) は無視。戻り: (callable δ_s(x), 診断 dict)。"""
    from scipy.interpolate import BSpline
    x = np.asarray(x, dtype=float); d = np.asarray(d, dtype=float)
    w = np.ones_like(d) if weights is None else np.asarray(weights, dtype=float)
    ok = np.isfinite(d) & (w > 0)
    x0 = float(x[ok].min())
    xi = lambda xx: np.sqrt(np.maximum(xx - x0, 0.0) + x_stretch)
    k = 5
    xi_all = xi(x[ok]); xi_lo, xi_hi = float(xi_all.min()), float(xi_all.max())
    # ノット間隔: 物理 x で knot_spacing 相当 (ξ 空間では dξ = dx/(2ξ) → 中央値で換算)
    dxi = knot_spacing / (2.0 * np.median(xi_all))
    n_int = max(int(np.ceil((xi_hi - xi_lo) / dxi)) - 1, 4)
    inner = np.linspace(xi_lo, xi_hi, n_int + 2)[1:-1]
    t = np.concatenate([[xi_lo] * (k + 1), inner, [xi_hi] * (k + 1)])
    n_cp = len(t) - k - 1
    B = BSpline.design_matrix(xi_all, t, k).toarray()
    D = np.diff(np.eye(n_cp), n=penalty_order, axis=0)
    W = w[ok][:, None]
    A = B.T @ (W * B) + lam * (D.T @ D)
    c = np.linalg.solve(A, B.T @ (w[ok] * d[ok]))
    spl = BSpline(t, c, k, extrapolate=False)

    def f(xx):
        xx = np.asarray(xx, dtype=float)
        q = np.clip(xi(xx), xi_lo, xi_hi)
        v = spl(q)
        return np.maximum(v, 0.0) if positive else v

    fit = f(x[ok])
    rel = (fit - d[ok]) / np.maximum(np.abs(d[ok]), 1e-6)
    diag = dict(n_cp=int(n_cp), knot_spacing=knot_spacing, lam=lam, penalty_order=penalty_order,
                resid_rel_rms=float(np.sqrt(np.mean(rel ** 2))), resid_abs_max=float(np.max(np.abs(fit - d[ok]))),
                resid_rel_p95=float(np.percentile(np.abs(rel), 95)))
    return f, diag
