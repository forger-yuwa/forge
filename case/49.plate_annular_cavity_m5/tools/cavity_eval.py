#!/usr/bin/env python3
"""case/49 キャビティ内部の評価 (plan §4.7 / §4.8)。

usage:
  python3 tools/cavity_eval.py RUN                    # 最終スナップショットの要約
  python3 tools/cavity_eval.py RUN --series           # 全スナップショットの時系列 CSV
  python3 tools/cavity_eval.py RUN --plot             # 深さプロファイル図

出す量 (plan §4.7 の定義):
  dT_floor / dT_mid / dT_mouth   すきま中央・指定深さ・周方向平均の T - Tw
  dT_up / dT_dn                  流れ方向 上流(θ=0°)/下流(θ=180°) の dT_mid  (偏心の指標)
  dT_left / dT_right             θ=±90°  (**全周計算でのみ意味**。半割では恒等)
  zpen(eps)                      T-Tw > eps を満たす最深点 [m] (eps は 10/25/50 K)
  mdot_in / mdot_out             開口面の ro*uz の負側/正側 絶対積分 [kg/s]
  mdot_imbalance                 net / max(in,out)   ← 検査列 (生の net は 0 へ行くので使わない)
  q_outer / q_cylside / q_floor  等温壁の入熱 [W]  (**物理量**: q_w = lambda_w dT/dn)
  q_cyltop                       内円柱上面 (CV 外なので別枠)
  budget                         (Σq_3壁 + 開口エンタルピー流束) / Σq_3壁

注意 (plan §4.8, codex plan-2 M6):
  ここで出す q_w は**物理量の推定**であり、node の壁 Dirichlet (壁ノードのエネルギー残差を 0 化)
  による離散的なエネルギー授受と同一ではない。離散保存の検算は別実装 (残作業 #4)。
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402
from setup import load as load_conditions  # noqa: E402

MU0, T0S, SS = 1.716e-5, 273.0, 111.0


def mu_suth(T):
    return MU0 * (T / T0S) ** 1.5 * (T0S + SS) / (T + SS)


class Field:
    """node 値の IDW 補間 (壁ピン値を拾わないよう、必要に応じマスクした部分集合で作る)。"""

    def __init__(self, coord, vals, mask=None, k=6):
        self.k = k
        idx = np.where(mask)[0] if mask is not None else np.arange(len(coord))
        self.idx = idx
        self.tree = cKDTree(coord[idx])
        self.vals = {n: v[idx] for n, v in vals.items()}

    def at(self, pts, name):
        d, j = self.tree.query(pts, k=self.k)
        d = np.maximum(d, 1e-12)
        w = 1.0 / d ** 2
        return np.sum(self.vals[name][j] * w, axis=-1) / np.sum(w, axis=-1)


def run_conditions(run):
    """**run 自身に固定された作動条件**を読む。無ければ共有 case.json に落ちる (旧 run 用)。

    共有 `case.json` を読み直すと、M や壁温を変えた時点で**過去 run の評価が変わる**
    (2026-09-19 codex Major 10)。`gen_runs.py` は run 作成時に `conditions.json` を書く。
    """
    p = Path(run) / "conditions.json"
    if p.exists():
        return json.loads(p.read_text())
    print("  WARNING: %s に conditions.json が無い -> 共有 case.json を使う "
          "(共有入力を変えるとこの run の評価が変わる)" % run)
    return load_conditions()


def taw_of(D):
    """**run が実際に使った EOS** の回復温度。CPG run で TP の Taw を使うと h_aw が過大になる
    (2026-09-19 codex Major 10: 本ケースの温度差で約 9.8 %)。"""
    gas = str(D.get("gas_used", D.get("gas", "CPG"))).upper()
    if gas == "TP" and D.get("Taw_tp") is not None:
        return float(D["Taw_tp"])
    return float(D["Taw_cpg"])


def read(res):
    with h5py.File(res, "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        keys = [k for k in ("ro", "Ux", "Uy", "Uz", "P", "T", "h0", "k", "omega",
                            "vis_lam", "vis_turb", "wall_dist") if "VALUE/" + k in f]
        v = {k: np.array(f["VALUE/" + k]) for k in keys}
    return c, v


def wall_dump(run, physid, step, name=None):
    """res_wall_<physID>_<step>.h5 を読む。ソルバが残差に入れたのと同じ q_w / tau_w / u_tau。
    要素面積 (CONNE) からノード重みを作り、面積分と面積平均を返す。"""
    # 壁面ダンプは res_<境界名>_<physID>_<step>.h5
    p = Path(run) / ("res_%s_%d_%d.h5" % (name, physid, step))
    if not p.exists():
        cands = sorted(Path(run).glob("res_*_%d_%d.h5" % (physid, step)))
        if not cands:
            return None
        p = cands[0]
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "solver_density_cuda" / "tools"))
    from res_h5_to_vtu import parse_conne
    with h5py.File(p, "r") as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        conne = np.array(f["MESH/CONNE"])
        v = {k: np.array(f["VALUE/" + k]) for k in f["VALUE"]}
    n = len(xyz)
    # 要素数は CONNE をなめて数える (Mixed)
    ncell, i = 0, 0
    while i < len(conne):
        code = int(conne[i])
        nn = {2: 2, 4: 3, 5: 4, 6: 4, 16: 6, 9: 8}.get(code)
        if nn is None:
            break
        i += 1 + nn
        ncell += 1
    conn, offs, _ = parse_conne(conne, ncell)
    warea = np.zeros(n)
    s0 = 0
    for e in range(ncell):
        nd = conn[s0:offs[e]]
        s0 = offs[e]
        pts = xyz[nd]
        if len(nd) == 3:
            a = 0.5 * np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0]))
        elif len(nd) == 4:
            a = (0.5 * np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0]))
                 + 0.5 * np.linalg.norm(np.cross(pts[2] - pts[0], pts[3] - pts[0])))
        else:
            continue
        warea[nd] += a / len(nd)
    return dict(xyz=xyz, w=warea, area=float(warea.sum()), **v)


def snapshots(run):
    return sorted(Path(run).glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))


# ---------------------------------------------------------------- 量
def eval_snapshot(c, v, man, D):
    G, E = man["geometry"], man["eval"]
    Tw, shrink = D["wall_T"], E["shrink_m"]
    out = {}
    # 壁ピン値を避けたガス場だけで補間器を作る
    m_gas = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man, shrink=shrink)
    if m_gas.sum() < 50:
        raise SystemExit("キャビティ内のガスノードが少なすぎる (%d)" % m_gas.sum())
    fg = Field(c, v, mask=m_gas)
    fa = Field(c, v)                                   # 全域 (開口面の流束用)

    th = np.linspace(0.0, np.pi, E["n_theta"])
    w = gc.azimuth_weights(th, man)
    for nm, frac in E["probe_depth_frac"].items():
        pts, _ = gc.probe_points(man, frac, theta=th)
        T = fg.at(pts, "T")
        out["dT_" + nm] = float(np.sum(w * (T - Tw)))
    # 方位別 (深さは mid)
    pts_mid, _ = gc.probe_points(man, E["probe_depth_frac"]["mid"], theta=th)
    T_mid = fg.at(pts_mid, "T") - Tw
    # θ=0 が上流、θ=180° が下流 (2026-09-19 ユーザ指定の向き)
    for nm, ang in (("up", 0.0), ("dn", np.pi), ("left", 0.5 * np.pi), ("right", 0.5 * np.pi)):
        i = int(np.argmin(np.abs(th - ang)))
        out["dT_" + nm] = float(T_mid[i])
    # 侵入深さ: 深さ方向に細かく取り、周方向平均 dT がしきい値を超える最深点
    nz = 201
    fr = np.linspace(0.0, 1.0, nz)
    prof = np.empty(nz)
    for i, f_ in enumerate(fr):
        pts, _ = gc.probe_points(man, f_, theta=th)
        prof[i] = np.sum(w * (fg.at(pts, "T") - Tw))
    out["_depth_frac"] = fr
    out["_dT_profile"] = prof
    for eps in E["zpen_eps_K"]:
        idx = np.where(prof > eps)[0]
        out["zpen_%g" % eps] = float(fr[idx[-1]] * G["depth"]) if idx.size else 0.0

    # 開口面 (z = -shrink) の質量流束: 極座標グリッドで面積分。
    # **線形補間を使う** (IDW は符号が変わる u_z を平滑化して流入/流出を潰す)。
    # 開口近傍のノードだけで Delaunay を張る (全域だと重い)。
    # 開口面の評価深さ: 開口直下すぎると壁ピン (u=0) のノードを拾って流入/流出を潰す。
    # 実測 (2026-09-19): z=-0.2 mm で質量収支 +33 %、z=-0.5〜-1 mm で ±5 % に収まる。
    # **補間は 2D (x,y) で張る**。評価面は 1 枚の平面なので 3 次元は要らないうえ、
    # 厚さ 3 mm の薄いスラブの 3D Delaunay は退化点群になって病的に遅い
    # (2026-09-19 実測: 2.74M 節点の run で 40 分以上終わらなかった)。
    # 構造化 (ヘキサ) では評価面付近に**節点層**があるので、その 1 層だけを取れば
    # z 方向の混ぜ込みも無くなる。非構造では層が無いので薄いスラブで代用する。
    zf = E.get("flux_depth_frac", 0.02) * G["depth"]
    inr = np.hypot(c[:, 0], c[:, 1]) < G["Ro"] + 3.0e-3
    lin = {}
    zc = c[inr & (np.abs(c[:, 2] + zf) < 3.0e-3), 2]
    slab = None
    if zc.size > 200:
        zu = np.unique(np.round(zc, 9))
        zpick = zu[np.argmin(np.abs(zu + zf))]
        layer = inr & (np.abs(c[:, 2] - zpick) < 1e-9)
        if layer.sum() > 200:                     # 節点層がある (構造化)
            slab = layer
            print("  開口評価: z=%.4f mm の節点層 %d 点で 2D 補間"
                  % (zpick * 1e3, int(layer.sum())), flush=True)
        else:                                     # 層が無い (非構造) -> 薄いスラブ
            slab = inr & (np.abs(c[:, 2] + zf) < 1.0e-3)
            print("  開口評価: 厚さ ±1 mm のスラブ %d 点で 2D 補間 (節点層なし)"
                  % int(slab.sum()), flush=True)
    if slab is not None and slab.sum() > 200:
        tri_pts = c[slab][:, :2]                  # (x, y) のみ
        for nm in ("ro", "Uz", "h0"):
            if nm in v:
                lin[nm] = LinearNDInterpolator(tri_pts, v[nm][slab])
    nth, nr = 361, 60
    tg = np.linspace(0.0, np.pi, nth)
    ri = gc.inner_radius_at(tg, man)
    frac_r = (np.arange(nr) + 0.5) / nr
    TH = np.repeat(tg[:, None], nr, axis=1)
    RR = ri[:, None] + frac_r[None, :] * (G["Ro"] - ri[:, None])
    dA = (G["Ro"] - ri[:, None]) / nr * RR * (np.pi / (nth - 1))
    DXo, DYo = gc.ray_dir(TH)          # **θ=0 が上流**。中央面と同じ符号規約にする
    pts = np.stack([RR * DXo, RR * DYo, np.full_like(RR, -zf)], axis=-1).reshape(-1, 3)
    def samp(nm):
        if nm in lin:
            val = lin[nm](pts[:, :2])             # 2D 補間 (x, y)
            bad = ~np.isfinite(val)
            if bad.any():
                val[bad] = fa.at(pts[bad], nm)
            return val.reshape(nth, nr)
        return fa.at(pts, nm).reshape(nth, nr)
    ro = samp("ro")
    uz = samp("Uz")
    h0 = samp("h0") if "h0" in v else None
    flux = ro * uz * dA
    out["mdot_in"] = float(-flux[flux < 0].sum())
    out["mdot_out"] = float(flux[flux > 0].sum())
    out["mdot_net"] = float(flux.sum())
    den = max(out["mdot_in"], out["mdot_out"], 1e-30)
    out["mdot_imbalance"] = out["mdot_net"] / den
    out["_H_open"] = float(np.sum(flux * h0)) if h0 is not None else float("nan")
    out["_flux_z_m"] = float(-zf)

    # NOTE: 場の勾配から q_w を組む案は**この種のメッシュでは使えない** (2026-09-19 実測)。
    # 壁から数十 µm の点を IDW/線形で補間すると、接線 1-2 mm 間隔のノードを拾って
    # dT/dn をソルバ値の 7 倍に出す。**壁熱流束はソルバ出力 (wall_heat) を一次情報とし**、
    # 離散保存の独立検算は残作業 #4 (非拘束 CV 群の流束収支) で行う。
    return out


def mid_surface(man, n_th=181, n_z=200, with_axes=False):
    """すきま中央面の点群 [m]。**底面は「すきま幅の半分」だけ切り取る** (ユーザ指定 2026-09-19):
    z は 0 から -(depth - gap_local(θ)/2) まで。こうすると床の壁点に対する最近傍点が
    床から gap/2 上になり、側壁の壁点が中央面まで gap/2 なのと整合する。"""
    G = man["geometry"]
    th = np.linspace(0.0, np.pi, n_th)
    rc = gc.gap_center_radius(th, man)
    gap = gc.gap_at(th, man)
    fr = np.linspace(0.0, 1.0, n_z)
    z = -(G["depth"] - 0.5 * gap)[:, None] * fr[None, :]    # (n_th, n_z)
    R = np.repeat(rc[:, None], n_z, axis=1)
    TH = np.repeat(th[:, None], n_z, axis=1)
    # **θ=0 は上流 (-x)**。`gc.ray_dir` と符号を揃える。+cos で組むと偏心時に面が x 鏡映され、
    # 固体内に入る (2026-09-19 codex Critical: x_off=1.5 mm で 181 方位中 56 方位が固体内、
    # 底の切り取りも下流 0.5 mm 必要なところ 2.0 mm と逆になっていた)。
    DX, DY = gc.ray_dir(TH)
    pts = np.stack([R * DX, R * DY, z], axis=-1).reshape(-1, 3)
    return (pts, th, fr) if with_axes else pts


def wall_y1plus(wh, man, c, v):
    r"""**第一内部ノードの** $y_1^+ = \rho u_\tau y_1/\mu$ を後処理で出す。

    ソルバの壁面ダンプ `ypls` は `ρ u_τ dcc/μ` で、`dcc` は**ゴーストセル重心と内点セル重心の
    距離**である ([`viscousFlux_d.cu`](../../../solver_density_cuda/cuda_forge/viscousFlux_d.cu))。
    **node 方式では壁ノードが壁面上に乗るので dcc がほぼ 0 に退化し**、y⁺ が実際より桁違いに
    小さく出る (ソルバ自身、流束計算ではこの退化を避けている)。よって壁解像の根拠に使えない
    (2026-09-19 codex Major 9)。

    ここでは場の `wall_dist` から**キャビティ気体ノードの最小正値** $y_1$ を取り、
    壁面ダンプの $u_\tau$ と壁の $\rho,\mu$ で $y_1^+$ を組む。構造化ヘキサなので $y_1$ は
    第一層厚そのものになる。
    """
    if "wall_dist" not in v:
        return
    m = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man, shrink=0.0)
    wd = np.asarray(v["wall_dist"], float)[m]
    pos = wd[wd > 1e-9]
    if pos.size == 0:
        return
    y1 = float(np.percentile(pos, 0.1))          # 最小値はガベージを拾い得るので下位 0.1 %
    for g, d in wh.items():
        ut = d.get("_utau_node")
        ro = d.get("_ro_node")
        mu = d.get("_mu_node")
        if ut is None or ro is None or mu is None:
            continue
        w = d["_w_node"]
        yp = ro * ut * y1 / np.maximum(mu, 1e-30)
        d["y1"] = y1
        d["y1plus_mean"] = float(np.sum(yp * w) / max(np.sum(w), 1e-30))
        d["y1plus_max"] = float(np.max(yp))


def wall_href(wh, man, D, c, v, T0):
    """基準温度 = すきま中央面の最近傍点の**総温 T0** で熱伝達率を出す (ユーザ指定 2026-09-19)。
        h_ref = q'' / (T0_ref - T_w)
    T0 は `tools/total_quantities.py` の `total_state` (h0 の逆算) で作った node 値を使う
    (スクリプトで T + u^2/2cp を組まない: AGENTS.md「出力と後処理の原則」)。"""
    Tw = D["wall_T"]
    n_th, n_z = 181, 200
    mid, th_ax, fr_ax = mid_surface(man, n_th, n_z, with_axes=True)
    m_gas = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man, shrink=man["eval"]["shrink_m"])
    f = Field(c, {"T0": T0}, mask=m_gas)
    T0_mid = f.at(mid, "T0").reshape(n_th, n_z)
    # **最近傍でなく双線形補間で引く**。最近傍だと基準温度が中央面の離散点に量子化され、
    # 壁に沿って階段状になる (2026-09-19 実測: 1 列 564 点でユニーク値 364)。
    # 中央面は (θ, 深さ比 fr) の規則格子なので、壁点をその 2 座標に写して補間する。
    from scipy.interpolate import RegularGridInterpolator
    rgi = RegularGridInterpolator((th_ax, fr_ax), T0_mid, bounds_error=False, fill_value=None)
    G = man["geometry"]
    for g, d in wh.items():
        if "_z_node" not in d:
            continue
        xw, yw, zw = d["_x_node"], d["_y_node"], d["_z_node"]
        thw = np.arctan2(np.abs(yw), -xw)                 # θ=0 が上流、0..pi
        gapw = gc.gap_at(np.clip(thw, 0.0, np.pi), man)
        zmax = np.maximum(G["depth"] - 0.5 * gapw, 1e-12)
        frw = np.clip(-zw / zmax, 0.0, 1.0)
        Tref = rgi(np.stack([np.clip(thw, th_ax[0], th_ax[-1]), frw], axis=1))
        qn = d["_qin_node"]
        w = d["_w_node"]
        # **分母を切り上げない** (旧実装は max(dT, 1 K) で、負の温度差まで +1 K に置換していた。
        # 2026-09-19 codex Major: 壁面積の 21 % がこの補正に掛かっており、ゼロ割ガードの域を
        # 越えて値を作っていた)。代わりに **係数を定義できない面積割合**を報告する。
        dT = Tref - Tw
        # **閾値は手で置かず、T0 の数値ノイズから決める** (2026-09-19)。
        #
        # h = q''/ΔT は ΔT→0 でも破綻しない: 深部では q'' も一緒に小さくなり、比は 30 W/m2K
        # 前後の妥当な値に落ち着く (実測 θ=140°: z=-8.1 mm で ΔT 0.35 K / q'' 11 W/m2 -> h 32)。
        # 実測では **ΔT が負になる面積は 0.00 %** なので、符号の問題も無い。
        # 閾値が要る唯一の理由は「T0 の数値ノイズに埋もれる点を平均に入れない」ことだけ。
        #
        # T0 は h0/cp を float32 の h0 から作るので、分解能は ulp(h0)/cp
        # (実測 h0 ~ 1.25e6 J/kg -> ulp 0.125 -> 1.24e-4 K)。その `href_noise_factor` 倍を閾値にする。
        # 旧既定 1.0 K は側壁面積の 23 %、0.05 K でも 3.2 % を切り落としていた
        # (ノイズ基準なら 0.16 %)。閾値を大きく取ると図に**階段状の境界**が出る。
        cp = float(D.get("cp", 1004.5))
        # 分解能は**実際に使った h0 の大きさ**から測る (T0 = h0/cp、h0 は float32 出力)
        h0f = np.asarray(v.get("h0", np.array([1.25e6])), float)
        h0rep = float(np.median(np.abs(h0f[np.isfinite(h0f)]))) if h0f.size else 1.25e6
        t0_noise = float(np.spacing(np.float32(h0rep))) / max(cp, 1e-30)
        nf = float(man["eval"].get("href_noise_factor", 10.0))
        dT_min = man["eval"].get("href_dT_min_K", None)
        if dT_min is None:
            dT_min = nf * t0_noise                        # 既定は**測って決める**
        d["href_dT_min_basis"] = ("T0 ノイズ %.3e K × %.3g" % (t0_noise, nf)
                                  if man["eval"].get("href_dT_min_K") is None else "手動指定")
        ok = dT > dT_min
        d["href_undef_area_frac"] = float(np.sum(w[~ok]) / max(np.sum(w), 1e-30))
        d["href_neg_area_frac"] = float(np.sum(w[dT <= 0.0]) / max(np.sum(w), 1e-30))
        d["href_dT_min_K"] = dT_min
        d["Tref_mean"] = float(np.sum(Tref * w) / max(np.sum(w), 1e-30))
        d["dTref_mean"] = d["Tref_mean"] - Tw              # **分母は温度差**。絶対温度で語らない
        d["Tref_min"] = float(Tref.min())
        d["Tref_max"] = float(Tref.max())
        # (1) 局所係数の面積平均 <q''/dT>_A  (ユーザ指定の定義。定義可能な面積のみで平均)
        d["h_ref"] = (float(np.sum((qn[ok] / dT[ok]) * w[ok]) / max(np.sum(w[ok]), 1e-30))
                      if ok.any() else float("nan"))
        # (2) 総入熱を再現する係数 h_eff = Q / ∫dT dA  (別物なので別名で出す)
        den = float(np.sum(dT * w))
        d["h_eff"] = float(np.sum(qn * w) / den) if abs(den) > 1e-30 else float("nan")
        hn = np.full_like(qn, np.nan, dtype=float)
        hn[ok] = qn[ok] / dT[ok]
        d["_href_node"] = hn
        d["_Tref_node"] = Tref
        d["_dT_node"] = dT
        if "_z" in d:                        # 深さ分布 (面積重み付き)
            zb = d["_z"]
            edges = np.linspace(d["_z_node"].min(), d["_z_node"].max(), len(zb) + 1)
            ib = np.clip(np.digitize(d["_z_node"], edges) - 1, 0, len(zb) - 1)
            num = np.bincount(ib, weights=(qn / dT) * w, minlength=len(zb))
            den = np.bincount(ib, weights=w, minlength=len(zb))
            d["_href"] = num / np.maximum(den, 1e-30)
            numT = np.bincount(ib, weights=Tref * w, minlength=len(zb))
            d["_Tref"] = numT / np.maximum(den, 1e-30)
    return wh


def wall_Q_below(wh, groups, z_cut):
    """z < z_cut の壁だけの入熱 [W] (開口面 z_cut での流束と収支を組むため)。"""
    tot = 0.0
    for g in groups:
        d = wh.get(g)
        if d is None or "_z_node" not in d:
            continue
        m = d["_z_node"] < z_cut
        tot += float(np.sum(d["_qin_node"][m] * d["_w_node"][m]))
    return tot


def wall_heat(run, step, man, D, prof_n=40):
    """壁グループごとに **ソルバ出力の q_w** (res_wall) から入熱・熱流束・熱伝達率を出す。

    符号: forge の qwall は「壁→流体が正」。ガスが壁を加熱する (壁に入る) 側を正にしたいので
    q_in = -qwall とする。熱伝達率は 2 通り出す:
      h_aw  = q'' / (T_aw - T_w)   … 外部流基準 (設計で使いやすい。T_aw は採用 EOS の値)
      h_loc = q'' / (T_gas - T_w)  … その深さのすきま中央ガス温度基準 (キャビティ内部で物理的)
    """
    PID = man["phys_id"]
    Taw = taw_of(D)
    Tw = D["wall_T"]
    res = {}
    for g in ("cav_outer", "cyl_side", "cav_floor", "cyl_top", "plate", "plate_in"):
        if g not in PID:
            continue
        w = wall_dump(run, PID[g], step, name=g)
        if w is None or "qwall" not in w:
            continue
        qin = -np.asarray(w["qwall"], float)          # 壁に入る側を正
        wt = w["w"]
        Q = float(np.sum(qin * wt))                   # [W] (半割)
        A = float(np.sum(wt))
        qpp = Q / max(A, 1e-30)                       # 面積平均 [W/m2]
        tw = np.sqrt(np.asarray(w["twall_x"], float) ** 2 + np.asarray(w["twall_y"], float) ** 2
                     + np.asarray(w["twall_z"], float) ** 2)
        # **開口リップ (90° の鋭角) は幾何的特異点**で、細分すると q'' が h^-1/2 で発散する
        # (2026-09-19 実測: q''max 290 -> 355 -> 419 kW/m2、増分比が r^-0.5 に一致)。
        # そのため総入熱 Q は観測次数 0.24 の遅い収束になり Richardson 外挿の不確かさが 32 % になる。
        # **リップ帯を除いた入熱**を併記して、残りが収束していることを示す
        # (帯幅は `eval.lip_band_m`、既定 1 mm = 深さの 2 %)。
        lip = man["eval"].get("lip_band_m", 1.0e-3)
        zn = w["xyz"][:, 2]
        deep = zn < -lip
        Q_nolip = float(np.sum(qin[deep] * wt[deep])) if deep.any() else float("nan")
        d = dict(Q_W=Q, Q_nolip_W=Q_nolip, lip_band_m=lip,
                 area_m2=A, qpp_mean=qpp, qpp_max=float(np.max(qin)),
                 qpp_min=float(np.min(qin)),
                 h_aw=qpp / max(Taw - Tw, 1e-30),
                 tau_mean=float(np.sum(tw * wt) / max(A, 1e-30)),
                 ypls_max=float(np.max(w["ypls"])) if "ypls" in w else float("nan"),
                 ypls_mean=float(np.sum(np.asarray(w["ypls"], float) * wt) / max(A, 1e-30))
                 if "ypls" in w else float("nan"),
                 all_zero=bool(np.all(qin == 0.0)))
        # y₁⁺ を後処理で組むための壁面値 (ソルバの `ypls` は dcc 基準で node では退化する)
        if "utau" in w:
            d["_utau_node"] = np.asarray(w["utau"], float)
        if "ro" in w:
            d["_ro_node"] = np.asarray(w["ro"], float)
        if "Ts" in w:
            d["_mu_node"] = mu_suth(np.asarray(w["Ts"], float))
        d["_x_node"] = w["xyz"][:, 0]
        d["_y_node"] = w["xyz"][:, 1]
        d["_z_node"] = w["xyz"][:, 2]
        d["_qin_node"] = qin
        d["_w_node"] = wt
        # 深さ (z) 方向の分布 (側壁) / 半径方向 (床)
        z = w["xyz"][:, 2]
        if g in ("cav_outer", "cyl_side"):
            edges = np.linspace(z.min(), z.max(), prof_n + 1)
            ib = np.clip(np.digitize(z, edges) - 1, 0, prof_n - 1)
            num = np.bincount(ib, weights=qin * wt, minlength=prof_n)
            den = np.bincount(ib, weights=wt, minlength=prof_n)
            d["_z"] = 0.5 * (edges[1:] + edges[:-1])
            d["_qpp"] = num / np.maximum(den, 1e-30)
        res[g] = d
    return res


SERIES_COLS = ["dT_floor", "dT_mid", "dT_mouth", "dT_up", "dT_dn",
               "zpen_25", "mdot_in", "mdot_out", "mdot_imbalance",
               "q_outer", "q_cylside", "q_floor",   # q_* は wall_heat で埋める
               "q_wall_sum", "h_ref"]               # 準定常判定 (check_cavity_steady.py) の対象量


def lip_scan(wh, eps):
    r"""開口リップ帯の幅 ε を振って Q(ε) を出し、√ε 外挿で総入熱を押さえる。

    鋭角凸エッジの $q''\sim s^{-1/2}$ が正しければ、リップから距離 ε までの寄与は
    $\int_0^\varepsilon q''ds \sim 2C\sqrt{\varepsilon}$ なので
        Q(ε) = Q_tot - A√ε        (Q(ε) = 上端 ε を除いた入熱)
    になる。**この直線性が成り立つこと自体が特異性の確認**であり、切片が総入熱の推定値。
    格子を変えても切片が一致すれば、総入熱は格子に依らず決まっている。
    """
    cav = [g for g in ("cav_outer", "cyl_side", "cav_floor") if g in wh]
    if not cav:
        return
    print("  --- リップ帯スキャン: Q(ε) = 開口から ε を除いた入熱 (全周) ---")
    xs, ys = [], []
    for e in sorted(eps):
        q = 0.0
        for g in cav:
            zn, qin, w = wh[g]["_z_node"], wh[g]["_qin_node"], wh[g]["_w_node"]
            m = zn < -e
            q += float(np.sum(qin[m] * w[m])) if m.any() else 0.0
        print("    ε = %6.3f mm   Q = %8.4g W" % (e * 1e3, 2 * q))
        xs.append(np.sqrt(e)); ys.append(2 * q)
    def fit(xs, ys, tag):
        A = np.polyfit(xs, ys, 1)
        f = np.polyval(A, xs)
        r2 = 1.0 - np.sum((np.asarray(ys) - f) ** 2) / max(np.var(ys) * len(ys), 1e-30)
        print("    %s: Q(ε) = %.4g - %.4g*√ε   (R² = %.5f)  -> ε→0 外挿 **%.4g W**"
              % (tag, A[1], -A[0], r2, A[1]))
        return A[1]
    if len(xs) >= 3:
        small = [(x, y) for x, y in zip(xs, ys) if x * x <= 0.5e-3]
        if len(small) >= 3:
            fit([x for x, _ in small], [y for _, y in small], "小 ε 域 (≤0.5 mm) のみ")
        A = np.polyfit(xs, ys, 1)
        fit = np.polyval(A, xs)
        r2 = 1.0 - np.sum((np.asarray(ys) - fit) ** 2) / max(np.var(ys) * len(ys), 1e-30)
        Qtot = sum(wh[g]["Q_W"] for g in cav) * 2
        print("    √ε 直線あてはめ: Q(ε) = %.4g - %.4g*√ε   (R² = %.5f)" % (A[1], -A[0], r2))
        print("    -> ε→0 外挿 **Q_tot = %.4g W**   (格子上の実測 総 Q = %.4g W, 差 %+.2f %%)"
              % (A[1], Qtot, 100 * (Qtot / A[1] - 1)))
        print("    R² が 1 に近いほど q'' ~ s^-1/2 の仮定が効いている。")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--series", action="store_true")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--lip-scan", default=None,
                    help="リップ帯幅 [mm] をカンマ区切りで与え、Q(ε) と √ε 外挿を出す。"
                         "例 --lip-scan 0.1,0.2,0.5,1,2,5")
    ap.add_argument("--flux-depth", type=float, default=None,
                    help="開口流束の評価深さ [mm] (既定は manifest の flux_depth_frac)")
    a = ap.parse_args()
    man = gc.load_manifest(run=a.run)      # run が自分の manifest を持っていればそれを使う
    if a.flux_depth is not None:
        man["eval"]["flux_depth_frac"] = a.flux_depth * 1e-3 / man["geometry"]["depth"]
    D = run_conditions(a.run)
    snaps = snapshots(a.run)
    if not snaps:
        raise SystemExit("res_*.h5 が無い: %s" % a.run)

    if a.series:
        rows = []
        for res in snaps:
            c, v = read(res)
            q = eval_snapshot(c, v, man, D)
            st = int(res.stem.split("_")[1])
            wq = wall_heat(a.run, st, man, D)
            if not wq:                     # 壁ダンプの無いスナップショット (step 0 等) は飛ばす
                print("  step %6d  壁ダンプ無し -> skip" % st, flush=True)
                continue
            for nm, g in (("q_outer", "cav_outer"), ("q_cylside", "cyl_side"), ("q_floor", "cav_floor")):
                q[nm] = wq[g]["Q_W"] if g in wq else float("nan")
            q["q_wall_sum"] = q["q_outer"] + q["q_cylside"] + q["q_floor"]
            # **時系列でも h_ref を出す** (これが無いと準定常判定でスナップショット 0 点になり、
            # ユーザに報告する主量の定常性が確認できないまま通ってしまう。2026-09-19 修正)
            try:
                sys.path.insert(0, str(Path(__file__).resolve().parents[3]
                                       / "solver_density_cuda" / "tools"))
                from total_quantities import total_state
                T0s = np.asarray(total_state(a.run, str(res))["T0"], float)
                wall_href(wq, man, D, c, v, T0s)
            except Exception as e:                       # noqa: BLE001
                print("  step %6d  総温 T0 を作れない (%s) -> h_ref なし" % (st, e), flush=True)
            # 3 壁をまとめた h_ref (面積重み)。**判定量は run 間で同じ定義**にする
            aw = [(wq[g].get("area_m2", float("nan")), wq[g].get("h_ref", float("nan")))
                  for g in ("cav_outer", "cyl_side", "cav_floor") if g in wq]
            aw = [(A, h) for A, h in aw if A == A and h == h]
            q["h_ref"] = (sum(A * h for A, h in aw) / sum(A for A, h in aw)) if aw else float("nan")
            rows.append((st, q))
            print("  step %6d  dT_mid %8.2f  dT_mouth %8.2f  zpen25 %6.2f mm  q_sum %9.3f W"
                  % (rows[-1][0], q["dT_mid"], q["dT_mouth"], q["zpen_25"] * 1e3, q["q_wall_sum"]),
                  flush=True)
        out = Path(a.out or (Path(a.run) / "cavity_series.csv"))
        with open(out, "w") as f:
            f.write("step," + ",".join(SERIES_COLS) + "\n")
            for st, q in rows:
                f.write("%d," % st + ",".join("%.9g" % q[k] for k in SERIES_COLS) + "\n")
        print("wrote", out)
        print("  use: check_quasisteady.py --series-csv %s --series-cols %s"
              % (out, ",".join(SERIES_COLS)))
        return

    c, v = read(snaps[-1])
    step = int(snaps[-1].stem.split("_")[1])
    q = eval_snapshot(c, v, man, D)
    Tw = D["wall_T"]
    Taw = taw_of(D)
    wh = wall_heat(a.run, step, man, D)
    # 総温 T0 (h0 の逆算; CPG/TP 両対応) — AGENTS.md「出力と後処理の原則」
    T0 = None
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "solver_density_cuda" / "tools"))
        from total_quantities import total_state
        T0 = np.asarray(total_state(a.run, str(snaps[-1]))["T0"], float)
    except Exception as e:                       # noqa: BLE001
        print("  WARNING: 総温 T0 を作れない (%s) -> 基準温度基準の h は出せない" % e)
    if wh:
        wall_y1plus(wh, man, c, v)
    if T0 is not None and wh:
        wall_href(wh, man, D, c, v, T0)
    print("=== %s  (%s) ===" % (a.run, snaps[-1].name))
    if wh:
        print("  --- 壁面 (ソルバ出力 q_w, 壁に入る側が正。半割) ---")
        print("  %-10s %9s %10s %10s %8s %9s %9s %8s" %
              ("group", "Q[W]", "q''[W/m2]", "q''max", "h_aw", "dT_ref[K]", "h_ref", "y+ mean"))
        for g, d in wh.items():
            flag = ("  (断熱壁なので 0 が正しい)" if d["all_zero"] and g in ("plate", "plate_in")
                    else ("  (全点 0: qwall 診断の無いバイナリ)" if d["all_zero"] else ""))
            print("  %-10s %9.4g %10.4g %10.4g %8.4g %9.4g %9.4g %8.3g%s" %
                  (g, d["Q_W"], d["qpp_mean"], d["qpp_max"], d["h_aw"],
                   d.get("dTref_mean", float("nan")), d.get("h_ref", float("nan")),
                   d["ypls_mean"], flag))
        for g, d in wh.items():
            if "h_eff" in d:
                print("    %-10s h_eff = Q/∫dT dA = %8.4g W/m2K   (dT_ref 平均 %.2f K; "
                      "dT<=%.3g K = %.2f %% [%s]; dT<0 = %.2f %%)"
                      % (g, d["h_eff"], d.get("dTref_mean", float("nan")),
                         d.get("href_dT_min_K", 0.0),
                         100 * d.get("href_undef_area_frac", 0.0),
                         d.get("href_dT_min_basis", "-"),
                         100 * d.get("href_neg_area_frac", 0.0)))
        y1 = next((d["y1"] for d in wh.values() if "y1" in d), None)
        if y1 is not None:
            print("    --- 壁解像 (**参考値**。正式判定は "
                  "`solver_density_cuda/tools/check_wall_resolution.py <run>` を使うこと: "
                  "ここの u_τ はソルバ出力で node 経路では高せん断域が stale、y1 も全域共通) ---")
            print("    第一層厚 y1 = %.4g m" % y1)
            for g, d in wh.items():
                if "y1plus_mean" in d:
                    print("    %-10s y1+ 平均 %8.3g / 最大 %8.3g   (ソルバ ypls 平均 %8.3g)"
                          % (g, d["y1plus_mean"], d["y1plus_max"], d["ypls_mean"]))
        print("    h_aw  = q''/(T_aw - T_w)       … 外部流の回復温度基準")
        print("    h_ref = <q''/(T0_ref - T_w)>_A … 局所係数の面積平均。**分母は温度差 dT_ref**"
              " (絶対温度ではない)。h_eff = Q/∫dT dA は総入熱を再現する別の係数。")
        print("    h_ref の元定義 = q''/(T0_ref - T_w)  … **基準温度** = すきま中央面の最近傍点の総温"
              " (底面はすきま幅の半分で切取り)")
        cav = [g for g in ("cav_outer", "cyl_side", "cav_floor") if g in wh]
        if cav:
            Qc = sum(wh[g]["Q_W"] for g in cav)
            Ac = sum(wh[g]["area_m2"] for g in cav)
            print("  キャビティ 3 壁 合計 Q = %.4g W (半割) = %.4g W (全周),  平均 q'' = %.4g W/m2,"
                  "  h_aw = %.4g W/m2K   [T_aw-T_w = %.1f K]"
                  % (Qc, 2 * Qc, Qc / Ac, Qc / Ac / (Taw - Tw), Taw - Tw))
            Qn = sum(wh[g].get("Q_nolip_W", float("nan")) for g in cav)
            lipb = wh[cav[0]].get("lip_band_m", 1e-3)
            print("    うち開口リップ帯 (上端 %.1f mm) を除く Q = %.4g W (半割) = %.4g W (全周)"
                  "   [リップ帯 %.4g W = %.1f %%]"
                  % (lipb * 1e3, Qn, 2 * Qn, Qc - Qn, 100 * (Qc - Qn) / max(Qc, 1e-30)))
            print("    ** リップは 90 度の鋭角 = 幾何的特異点で局所 q'' が h^-1/2 で発散する。"
                  "ただし ∫q'' ds ~ 2√ε なので **総 Q 自体は有限で収束する** (収束が遅いだけ)。"
                  "実測の観測次数 0.24 は『総 Q の格子不確かさが大きく精度を確定できない』であって"
                  "『発散する』ではない。リップ帯を除いた Q は限定領域の別指標 **")
    if a.lip_scan and wh:
        lip_scan(wh, [float(x) * 1e-3 for x in a.lip_scan.split(",") if x])
    print("  Tw = %.1f K,  Taw(CPG/TP) = %.1f / %.1f K  -> **使用 %s の %.1f K** (gas=%s)"
          % (Tw, D["Taw_cpg"], D.get("Taw_tp", float("nan")),
             str(D.get("gas_used", D.get("gas", "CPG"))).upper(), taw_of(D),
             D.get("gas_used", "(run に記録なし)")))
    for k in ("dT_mouth", "dT_mid", "dT_floor", "dT_up", "dT_dn"):
        print("  %-10s %9.2f K   (T = %8.2f K)" % (k, q[k], Tw + q[k]))
    for eps in man["eval"]["zpen_eps_K"]:
        print("  zpen(%2g K)  %9.2f mm" % (eps, q["zpen_%g" % eps] * 1e3))
    print("  mdot in/out %.4e / %.4e kg/s   imbalance %.2e" %
          (q["mdot_in"], q["mdot_out"], q["mdot_imbalance"]))
    if wh and np.isfinite(q["_H_open"]):
        zc = q["_flux_z_m"]
        qs = wall_Q_below(wh, ("cav_outer", "cyl_side", "cav_floor"), zc)
        qall = sum(wh[g]["Q_W"] for g in ("cav_outer", "cyl_side", "cav_floor") if g in wh)
        print("  CV 収支 (評価面 z=%.2f mm 以深): Σq_壁 %.3f W  vs  開口からの正味エンタルピー流入 %.3f W"
              "  -> 残差比 %.3f   [全深さの Σq は %.3f W]"
              % (zc * 1e3, qs, -q["_H_open"], (qs + q["_H_open"]) / max(abs(qs), 1e-30), qall))
    Path(Path(a.run) / "cavity_eval.json").write_text(json.dumps(
        {"field": {k: val for k, val in q.items() if not k.startswith("_")},
         "wall": {g: {k: val for k, val in d.items() if not k.startswith("_")}
                  for g, d in wh.items()}}, indent=2))
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
        for fp in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
            font_manager.fontManager.addfont(str(fp))
            matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()
        dep = man["geometry"]["depth"] * 1e3
        fig, ax = plt.subplots(1, 3, figsize=(13.5, 6.0), sharey=True)
        # (a) 温度プロファイル
        z = -q["_depth_frac"] * dep
        ax[0].plot(q["_dT_profile"] + Tw, z, lw=2.2, color="#b91c1c", label="静温 (すきま中央)")
        d0 = wh.get("cav_outer")
        if d0 and "_Tref" in d0:
            ax[0].plot(d0["_Tref"], d0["_z"] * 1e3, lw=1.8, color="#7c3aed",
                       label="基準温度 T0_ref (総温)")
        ax[0].axvline(Tw, ls="--", c="#1d4ed8", lw=1.4, label="壁温 %g K" % Tw)
        ax[0].axvline(Taw, ls=":", c="#111", lw=1.4, label="回復温度 %.0f K" % Taw)
        ax[0].set_xlabel("ガス温度 [K] (周方向平均・すきま中央)")
        ax[0].set_ylabel("深さ z [mm]")
        ax[0].legend(fontsize=9, loc="lower right")
        ax[0].set_title("(a) 温度の深さ分布", fontsize=12, loc="left")
        # (b) 壁熱流束 q'' の深さ分布
        for g, c_, lab in (("cav_outer", "#dc2626", "外筒壁 (r=%.0f mm)" % (man["geometry"]["Ro"] * 1e3)),
                           ("cyl_side", "#ea580c", "内円柱側面")):
            d = wh.get(g)
            if d and "_z" in d:
                ax[1].plot(d["_qpp"] * 1e-3, d["_z"] * 1e3, lw=2.0, color=c_, label=lab)
        ax[1].axvline(0, c="#999", lw=.8)
        ax[1].set_xlabel("壁熱流束 q'' [kW/m²]  (壁に入る側が正)")
        ax[1].legend(fontsize=9, loc="lower right")
        ax[1].set_title("(b) 熱流束の深さ分布", fontsize=12, loc="left")
        # (c) 熱伝達率: 基準温度 (すきま中央面の総温) 基準を主、回復温度基準を従で
        has_ref = any("_href" in (wh.get(g) or {}) for g in ("cav_outer", "cyl_side"))
        for g, c_, lab in (("cav_outer", "#dc2626", "外筒壁"), ("cyl_side", "#ea580c", "内円柱側面")):
            d = wh.get(g)
            if not d or "_z" not in d:
                continue
            if "_href" in d:
                ax[2].plot(d["_href"], d["_z"] * 1e3, lw=2.2, color=c_, label=lab + " (基準温度基準)")
            ax[2].plot(d["_qpp"] / (Taw - Tw), d["_z"] * 1e3, lw=1.4, ls="--", color=c_,
                       alpha=.75, label=lab + " (回復温度基準)")
        ax[2].axvline(0, c="#999", lw=.8)
        ax[2].set_xlabel("熱伝達率 h [W/m²K]")
        ax[2].legend(fontsize=8, loc="lower right")
        ax[2].set_title("(c) 熱伝達率  実線: h=q''/(T0_ref−T_w) / 破線: q''/(T_aw−T_w)",
                        fontsize=11, loc="left")
        for a_ in ax:
            a_.grid(alpha=.3)
            a_.set_ylim(-dep * 1.02, 2)
        qs = sum(wh[g]["Q_W"] for g in ("cav_outer", "cyl_side", "cav_floor") if g in wh)
        fig.suptitle("case/49  M%.2g 環状深キャビティ (すきま %.2f mm × 深さ %.0f mm, 壁 %g K)  — "
                     "キャビティ 3 壁 総入熱 %.1f W (全周)"
                     % (D["mach"], man["geometry"]["gap_nom"] * 1e3, dep, Tw, 2 * qs), fontsize=13)
        fig.tight_layout()
        fig.savefig(Path(a.run) / "cavity_profile.png", dpi=130, bbox_inches="tight")
        print("wrote", Path(a.run) / "cavity_profile.png")


if __name__ == "__main__":
    main()
