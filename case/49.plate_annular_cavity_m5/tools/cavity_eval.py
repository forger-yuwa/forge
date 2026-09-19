#!/usr/bin/env python3
"""case/49 キャビティ内部の評価 (plan §4.7 / §4.8)。

usage:
  python3 tools/cavity_eval.py RUN                    # 最終スナップショットの要約
  python3 tools/cavity_eval.py RUN --series           # 全スナップショットの時系列 CSV
  python3 tools/cavity_eval.py RUN --plot             # 深さプロファイル図

出す量 (plan §4.7 の定義):
  dT_floor / dT_mid / dT_mouth   すきま中央・指定深さ・周方向平均の T - Tw
  dT_up / dT_dn                  流れ方向 上流(θ=180°)/下流(θ=0°) の dT_mid  (偏心の指標)
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
    for nm, ang in (("up", np.pi), ("dn", 0.0), ("left", 0.5 * np.pi), ("right", 0.5 * np.pi)):
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

    # 開口面 (z = -shrink) の質量流束: 極座標グリッドで面積分
    nth, nr = 361, 60
    tg = np.linspace(0.0, np.pi, nth)
    ri = gc.inner_radius_at(tg, man)
    frac_r = (np.arange(nr) + 0.5) / nr
    TH = np.repeat(tg[:, None], nr, axis=1)
    RR = ri[:, None] + frac_r[None, :] * (G["Ro"] - ri[:, None])
    dA = (G["Ro"] - ri[:, None]) / nr * RR * (np.pi / (nth - 1))
    pts = np.stack([RR * np.cos(TH), RR * np.sin(TH), np.full_like(RR, -shrink)], axis=-1).reshape(-1, 3)
    ro = fa.at(pts, "ro").reshape(nth, nr)
    uz = fa.at(pts, "Uz").reshape(nth, nr)
    h0 = fa.at(pts, "h0").reshape(nth, nr) if "h0" in v else None
    flux = ro * uz * dA
    out["mdot_in"] = float(-flux[flux < 0].sum())
    out["mdot_out"] = float(flux[flux > 0].sum())
    out["mdot_net"] = float(flux.sum())
    den = max(out["mdot_in"], out["mdot_out"], 1e-30)
    out["mdot_imbalance"] = out["mdot_net"] / den
    out["_H_open"] = float(np.sum(flux * h0)) if h0 is not None else float("nan")

    # 等温壁の入熱 (照合用の**物理量推定**): q_w = lambda_w dT/dn, 片側 2 次差分。
    # **壁近傍を除外した fg ではなく全ノードの fa を使う** (shrink したマスクで壁際を
    # 評価すると勾配を数倍過大に出す。2026-09-19 実測で solver q_w の 7 倍になった)。
    # 差分間隔は VL 第一層に見合わせる (メッシュが解けていない間隔で差分しない)。
    lam_w = mu_suth(Tw) * D["cp"] / D["prandtl_lam"]
    gmin = G["gap_min"]
    h = max(2.0 * man["mesh"]["vl_first_m"], 0.01 * gmin)
    for nm, area_key in (("q_outer", "cav_outer"), ("q_cylside", "cyl_side"), ("q_floor", "cav_floor")):
        zs = np.linspace(-G["depth"] * 0.999, -G["depth"] * 0.001, 60)
        if nm == "q_floor":
            tg2 = np.linspace(0, np.pi, 91)
            r2 = np.linspace(0.02, 0.98, 12)
            ri2 = gc.inner_radius_at(tg2, man)
            RR2 = ri2[:, None] + r2[None, :] * (G["Ro"] - ri2[:, None])
            TH2 = np.repeat(tg2[:, None], len(r2), axis=1)
            base = np.stack([RR2 * np.cos(TH2), RR2 * np.sin(TH2),
                             np.full_like(RR2, -G["depth"])], axis=-1).reshape(-1, 3)
            nrm = np.tile(np.array([0.0, 0.0, 1.0]), (len(base), 1))
        else:
            tg2 = np.linspace(0, np.pi, 91)
            TH2 = np.repeat(tg2[:, None], len(zs), axis=1)
            ZZ = np.tile(zs, (len(tg2), 1))
            if nm == "q_outer":
                R = G["Ro"]
                base = np.stack([R * np.cos(TH2), R * np.sin(TH2), ZZ], axis=-1).reshape(-1, 3)
                nrm = np.stack([-np.cos(TH2), -np.sin(TH2), np.zeros_like(TH2)], axis=-1).reshape(-1, 3)
            else:
                R, off = G["Ri"], G["x_off"]
                base = np.stack([off + R * np.cos(TH2), R * np.sin(TH2), ZZ], axis=-1).reshape(-1, 3)
                nrm = np.stack([np.cos(TH2), np.sin(TH2), np.zeros_like(TH2)], axis=-1).reshape(-1, 3)
        T1 = fa.at(base + nrm * h, "T")
        T2 = fa.at(base + nrm * (2.0 * h), "T")
        dTdn = (-3.0 * Tw + 4.0 * T1 - T2) / (2.0 * h)
        out[nm] = float(-lam_w * np.mean(dTdn) * gc.wall_area(man, area_key))  # 壁に入る側を正
    out["q_wall_sum"] = out["q_outer"] + out["q_cylside"] + out["q_floor"]
    if np.isfinite(out["_H_open"]):
        out["budget"] = float((out["q_wall_sum"] + out["_H_open"]) / max(abs(out["q_wall_sum"]), 1e-30))
    return out


def wall_heat(run, step, man, D, prof_n=40):
    """壁グループごとに **ソルバ出力の q_w** (res_wall) から入熱・熱流束・熱伝達率を出す。

    符号: forge の qwall は「壁→流体が正」。ガスが壁を加熱する (壁に入る) 側を正にしたいので
    q_in = -qwall とする。熱伝達率は 2 通り出す:
      h_aw  = q'' / (T_aw - T_w)   … 外部流基準 (設計で使いやすい。T_aw は採用 EOS の値)
      h_loc = q'' / (T_gas - T_w)  … その深さのすきま中央ガス温度基準 (キャビティ内部で物理的)
    """
    PID = man["phys_id"]
    Taw = D.get("Taw_tp", D["Taw_cpg"])
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
        d = dict(Q_W=Q, area_m2=A, qpp_mean=qpp, qpp_max=float(np.max(qin)),
                 qpp_min=float(np.min(qin)),
                 h_aw=qpp / max(Taw - Tw, 1e-30),
                 tau_mean=float(np.sum(tw * wt) / max(A, 1e-30)),
                 ypls_max=float(np.max(w["ypls"])) if "ypls" in w else float("nan"),
                 ypls_mean=float(np.sum(np.asarray(w["ypls"], float) * wt) / max(A, 1e-30))
                 if "ypls" in w else float("nan"),
                 all_zero=bool(np.all(qin == 0.0)))
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
               "q_outer", "q_cylside", "q_floor"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--series", action="store_true")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    man = gc.load_manifest()
    D = load_conditions()
    snaps = snapshots(a.run)
    if not snaps:
        raise SystemExit("res_*.h5 が無い: %s" % a.run)

    if a.series:
        rows = []
        for res in snaps:
            c, v = read(res)
            q = eval_snapshot(c, v, man, D)
            rows.append((int(res.stem.split("_")[1]), q))
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
    Taw = D.get("Taw_tp", D["Taw_cpg"])
    wh = wall_heat(a.run, step, man, D)
    print("=== %s  (%s) ===" % (a.run, snaps[-1].name))
    if wh:
        print("  --- 壁面 (ソルバ出力 q_w, 壁に入る側が正。半割) ---")
        print("  %-10s %10s %10s %10s %10s %10s %9s" %
              ("group", "Q[W]", "q''[W/m2]", "q''max", "h_aw", "tau[Pa]", "y+ mean"))
        for g, d in wh.items():
            flag = ("  (断熱壁なので 0 が正しい)" if d["all_zero"] and g in ("plate", "plate_in")
                    else ("  (全点 0: qwall 診断の無いバイナリ)" if d["all_zero"] else ""))
            print("  %-10s %10.4g %10.4g %10.4g %10.4g %10.4g %9.3g%s" %
                  (g, d["Q_W"], d["qpp_mean"], d["qpp_max"], d["h_aw"], d["tau_mean"],
                   d["ypls_mean"], flag))
        cav = [g for g in ("cav_outer", "cyl_side", "cav_floor") if g in wh]
        if cav:
            Qc = sum(wh[g]["Q_W"] for g in cav)
            Ac = sum(wh[g]["area_m2"] for g in cav)
            print("  キャビティ 3 壁 合計 Q = %.4g W (半割) = %.4g W (全周),  平均 q'' = %.4g W/m2,"
                  "  h_aw = %.4g W/m2K   [T_aw-T_w = %.1f K]"
                  % (Qc, 2 * Qc, Qc / Ac, Qc / Ac / (Taw - Tw), Taw - Tw))
    print("  Tw = %.1f K,  Taw(CPG/TP) = %.1f / %.1f K" % (Tw, D["Taw_cpg"], D.get("Taw_tp", float("nan"))))
    for k in ("dT_mouth", "dT_mid", "dT_floor", "dT_up", "dT_dn"):
        print("  %-10s %9.2f K   (T = %8.2f K)" % (k, q[k], Tw + q[k]))
    for eps in man["eval"]["zpen_eps_K"]:
        print("  zpen(%2g K)  %9.2f mm" % (eps, q["zpen_%g" % eps] * 1e3))
    print("  mdot in/out %.4e / %.4e kg/s   imbalance %.2e" %
          (q["mdot_in"], q["mdot_out"], q["mdot_imbalance"]))
    print("  q_outer %9.3f W  q_cylside %9.3f W  q_floor %9.3f W  -> sum %9.3f W (半割)"
          % (q["q_outer"], q["q_cylside"], q["q_floor"], q["q_wall_sum"]))
    print("  budget (Σq + H_open)/Σq = %.3f" % q.get("budget", float("nan")))
    Path(Path(a.run) / "cavity_eval.json").write_text(json.dumps(
        {"field": {k: val for k, val in q.items() if not k.startswith("_")},
         "wall": {g: {k: val for k, val in d.items() if not k.startswith("_")}
                  for g, d in wh.items()}}, indent=2))
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.2, 6.4))
        ax.plot(q["_dT_profile"] + Tw, -q["_depth_frac"] * man["geometry"]["depth"] * 1e3, lw=2)
        ax.axvline(Tw, ls="--", c="r", label="wall %g K" % Tw)
        ax.axvline(D.get("Taw_tp", D["Taw_cpg"]), ls=":", c="k", label="T_aw")
        ax.set_xlabel("T [K] (周方向平均・すきま中央)"); ax.set_ylabel("z [mm]")
        ax.grid(alpha=.3); ax.legend()
        fig.tight_layout(); fig.savefig(Path(a.run) / "cavity_profile.png", dpi=130)
        print("wrote", Path(a.run) / "cavity_profile.png")


if __name__ == "__main__":
    main()
