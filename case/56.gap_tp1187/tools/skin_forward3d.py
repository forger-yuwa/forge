#!/usr/bin/env python3
"""3D の壁熱流束を**TP-1187 の測定手順と同じ演算**に通してから実測と比べる (観測モデル)。

TP-1187 は薄板 0.08 cm・304SS の温度上昇率から q = ρcτ dT/dt を出しており、薄板内の
横方向伝導を補正していない。**CFD の生の q_conv を実測と直接比べてはいけない**
(`acceptance.json` の T4-3D observation_model)。

## 薄板 (面内伝導つき 2 次元)

濡れ面を展開した (ℓ, z) 平面で

    ρcτ ∂T/∂t = λτ (∂²T/∂ℓ² + ∂²T/∂z²) + a(t) q₀(ℓ, z)

を解く。q₀ は CFD の定常壁熱流束 (壁へ入る向きを正、`qwall` の符号を反転)。
√(α t) ≈ 2.8 mm (t=2 s) なので **z の帯を先に平均して 1 次元化しない** (codex diagnose ①)。
円弧 (半径 r の円柱面) は可展面なので (ℓ, z) への展開は長さを保つ。4 辺は断熱で打ち切る
(TC から十分離す。打ち切り位置の感度は `--l-wall/--l-top/--zmax` で見る)。

## 時間演算 (原報 p.7 の逐語どおり)

- 解析時刻 t_a = 模型がトンネル中心線に達した時刻 (挿入開始から「a little over 2 s」)。
- 温度を Δ=0.25 s の区間で平均し、**区間間の差**から上昇率を出す:
  q̂ = ρcτ (T̄₊ − T̄₋)/Δ、T̄₋ = [t_a−Δ, t_a] 平均、T̄₊ = [t_a, t_a+Δ] 平均。
- 露出履歴 a(t) は分からないので 2 腕で挟む (codex diagnose ①):
  `step` (主仮定) = t=0 で a=1 / `ramp` (感度) = 0→t_a で a を 0→1 に線形、以後 1。
  分離形 q = a(t)·q₀ は近似であって挿入中の流れを再現したものではない。
- **分母も同じ演算に通す**: 実測 q_FP は平滑校正板を同じ手順で整理した値なので、
  一様熱流束に同じ演算を掛けた係数 f_FP で割る (step では 1、ramp では 0.979)。

## 熱電対の座標台帳 (原報 Fig 5 側面図、2026-09-26 に 8 倍描画で読み直し)

経路 ℓ は円弧の**壁側接点** (深さ r = 上面から 0.25 cm) を 0 とし、上面側を正に取る:

    前向き壁  ℓ = −(d − r)          d = 上面からの鉛直深さ
    円弧      ℓ = r θ,  θ∈[0, π/2]   (θ=0 が壁側接点、π/2 が上面側接点)
    上面      ℓ = r π/2 + (s − r)    s = 前向き壁面からの流れ方向距離

- 前向き壁 TC 92/91/90/89/88/87 = 深さ 0.25/0.51/0.76/1.52/2.54/3.81 cm (図に寸法明記)。
  TC92 は壁側接点そのもの。
- 上面 TC 94/95/96 = 壁面から 0.25/0.51/1.02 cm (図の寸法線は前向き壁面を基準に
  引かれている)。TC94 は上面側接点、以下 97/98/99 は 96 から 1.52/+1.52/+1.91 cm。
- TC93 は円弧上 (引出線は 45° 付近を指すが寸法が無い)。**93/94 は基準点が図から
  一意に読めない**ので名目位置 ±0.25 cm の区間の最小・最大も出す (handoff §1-a)。
- z は全点 0 (縦すきま中心線。原報 p.5–6 の ±0.03 cm は中心線に対するばらつき)。

TC 番号は L=15.24 cm のもの。L=30.48 cm 側は +25 (`--arm L30`)。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

CASE = Path(__file__).resolve().parents[1]
LAM_S, RHO_C = 15.0, 7900 * 500.0        # 304SS (全 TC・両 L で共通、当てはめで変えない)
ALPHA_S = LAM_S / RHO_C
TAU_TP1187 = 0.08e-2                     # 薄板厚 [m]
R_EDGE = 0.25e-2
# (TC, 種別, 値 [cm], 位置不確かさ [cm])
LEDGER = [
    (87, "wall", 3.81, 0.0), (88, "wall", 2.54, 0.0), (89, "wall", 1.52, 0.0),
    (90, "wall", 0.76, 0.0), (91, "wall", 0.51, 0.0), (92, "wall", 0.25, 0.0),
    (93, "arc", 45.0, 0.25),
    (94, "top", 0.25, 0.25), (95, "top", 0.51, 0.0), (96, "top", 1.02, 0.0),
    (97, "top", 2.54, 0.0), (98, "top", 4.06, 0.0), (99, "top", 5.97, 0.0),
]


def tc_ell(kind, v, r):
    if kind == "wall":
        return -(v * 1e-2 - r)
    if kind == "arc":
        return r * np.deg2rad(v)
    return r * np.pi / 2 + (v * 1e-2 - r)


def surface_points(run, step, W, r, D, tol=1e-6):
    """濡れ面の節点を (ℓ, z, q_in) で返す。q_in = −qwall (壁へ入る向きが正)。"""
    xd, xvd = 0.5 * W, 0.5 * W + r
    with h5py.File(CASE / run / f"res_gap_6_{step}.h5") as h:
        cg = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
        qg = -np.asarray(h["VALUE/qwall"], dtype=float)
    with h5py.File(CASE / run / f"res_plate_4_{step}.h5") as h:
        cp = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
        qp = -np.asarray(h["VALUE/qwall"], dtype=float)
    out = []
    m = (np.abs(cg[:, 0] - xd) < tol) & (cg[:, 1] <= -r + tol) & (cg[:, 1] >= -D - tol)
    out.append(np.c_[cg[m, 1] + r, cg[m, 2], qg[m]])                  # ℓ = y + r (≤0)
    rad = np.hypot(cg[:, 0] - xvd, cg[:, 1] + r)
    m = (np.abs(rad - r) < 5e-5) & (cg[:, 0] <= xvd + tol) & (cg[:, 1] >= -r - tol)
    th = np.arctan2(cg[m, 1] + r, -(cg[m, 0] - xvd))                   # 壁側 0 → 上面側 π/2
    out.append(np.c_[r * np.clip(th, 0, np.pi / 2), cg[m, 2], qg[m]])
    m = (np.abs(cp[:, 1]) < tol) & (cp[:, 0] >= xvd - tol)
    out.append(np.c_[r * np.pi / 2 + (cp[m, 0] - xvd), cp[m, 2], qp[m]])
    pts = np.vstack(out)
    # 接点の重複 (壁/円弧/上面の継ぎ目) を落とす
    key = np.round(pts[:, :2] / 1e-7).astype(np.int64)
    _, iu = np.unique(key, axis=0, return_index=True)
    return pts[np.sort(iu)], {k: int(len(o)) for k, o in zip(("wall", "arc", "top"), out)}


def exposure(kind, t_a):
    if kind == "step":
        return lambda t: 1.0
    return lambda t: min(t / t_a, 1.0)


def observe(q0, h, t_a, dwin, tau, kind, cfl=0.2):
    """q₀ [W/m²] の 2D 配列 (等間隔 h) に測定演算を掛けた q̂ を返す。q0 がスカラでもよい。"""
    q0 = np.asarray(q0, dtype=float)
    two_d = q0.ndim == 2
    dt0 = cfl * h * h / ALPHA_S if two_d else dwin / 200
    nwin = max(8, int(np.ceil(dwin / dt0)))
    dt = dwin / nwin                               # 窓の境界をステップ境界に揃える
    n_tot = int(round((t_a + dwin) / dt))
    i_m0, i_a = int(round((t_a - dwin) / dt)), int(round(t_a / dt))
    a = exposure(kind, t_a)
    src = q0 / (RHO_C * tau)
    T = np.zeros_like(q0)
    sm = np.zeros_like(q0); sp = np.zeros_like(q0)
    ih2 = ALPHA_S / (h * h)

    def lap(T):
        if not two_d:
            return 0.0
        P = np.pad(T, 1, mode="reflect")            # 4 辺断熱 (節点中心の鏡像: z=0 対称面が節点上)
        return ih2 * (P[2:, 1:-1] + P[:-2, 1:-1] + P[1:-1, 2:] + P[1:-1, :-2] - 4.0 * T)

    for n in range(n_tot + 1):
        # 台形則で区間平均を積む
        if i_m0 <= n <= i_a:
            sm += (0.5 if n in (i_m0, i_a) else 1.0) * T
        if i_a <= n <= n_tot:
            sp += (0.5 if n in (i_a, n_tot) else 1.0) * T
        if n == n_tot:
            break
        t = n * dt
        T = T + dt * (a(t + 0.5 * dt) * src + lap(T))
    Tm, Tp = sm / nwin, sp / nwin
    return RHO_C * tau * (Tp - Tm) / dwin


def load_ref(arm):
    acc = json.loads((CASE / "acceptance.json").read_text())
    g = next(x for x in acc["gates"] if x["id"] == "T4-3D")["reference"]
    if arm == "L15":
        d = g["L15_q_over_qFP"]
        return {int(k): (d["run8"][k], d["run14"][k]) for k in d["run8"]}
    d = g["L30_q_over_qFP"]
    return {int(k) - 25: (d["run8"][k], d["run14"][k]) for k in d["run8"]}


def snaps(run):
    return sorted(int(re.search(r"res_gap_6_(\d+)\.h5$", p).group(1))
                  for p in glob.glob(str(CASE / run / "res_gap_6_*.h5")))


def evaluate(a, step):
    pts, segn = surface_points(a.run, step, a.w, a.r, a.depth)
    l_lo, l_hi = -(a.l_wall - a.r), a.r * np.pi / 2 + (a.l_top - a.r)
    sel = (pts[:, 0] >= l_lo - 5 * a.h) & (pts[:, 0] <= l_hi + 5 * a.h) & (pts[:, 1] <= a.zmax + 5 * a.h)
    P = pts[sel]
    ll = np.arange(l_lo, l_hi + 0.5 * a.h, a.h)
    zz = np.arange(0.0, a.zmax + 0.5 * a.h, a.h)
    L, Z = np.meshgrid(ll, zz, indexing="ij")
    q0 = LinearNDInterpolator(P[:, :2], P[:, 2])(L, Z)
    bad = ~np.isfinite(q0)
    if bad.any():
        q0[bad] = NearestNDInterpolator(P[:, :2], P[:, 2])(L[bad], Z[bad])
    res = {k: observe(q0, a.h, a.ta, a.dwin, a.tau, k) for k in ("step", "ramp")}
    f_fp = {k: float(observe(1.0, a.h, a.ta, a.dwin, a.tau, k)) for k in ("step", "ramp")}

    def at(F, lv):                      # z=0 線で ℓ に線形補間
        return float(np.interp(lv, ll, F[:, 0]))

    rows = []
    for tc, kind, v, unc in LEDGER:
        lv = tc_ell(kind, v, a.r)
        if not (ll[0] <= lv <= ll[-1]):
            continue
        row = {"tc": tc, "kind": kind, "ell": lv, "raw": at(q0, lv)}
        for k in ("step", "ramp"):
            row[k] = at(res[k], lv) / f_fp[k]
            if unc > 0:
                m = (ll >= lv - unc * 1e-2) & (ll <= lv + unc * 1e-2)
                row[k + "_rng"] = (float(res[k][m, 0].min() / f_fp[k]), float(res[k][m, 0].max() / f_fp[k]))
        rows.append(row)
    info = {"segn": segn, "n_grid": q0.shape, "f_fp": f_fp,
            "z0_nodes": int((np.abs(P[:, 1]) < 1e-9).sum())}
    return rows, info


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--arm", choices=("L15", "L30"), default="L15", help="照合する実測の側")
    ap.add_argument("--step", type=int, default=None, help="省略時は最新")
    ap.add_argument("--series", action="store_true", help="--min-step 以上の全スナップショットを評価")
    ap.add_argument("--min-step", type=int, default=0)
    ap.add_argument("--csv", default=None, help="--series の時系列 CSV (check_quasisteady 用)")
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=R_EDGE)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--tau", type=float, default=TAU_TP1187)
    ap.add_argument("--ta", type=float, default=2.0, help="解析時刻 t_a [s] (暫定 2 s)")
    ap.add_argument("--dwin", type=float, default=0.25, help="平均化窓 Δ [s]")
    ap.add_argument("--h", type=float, default=1.0e-4, help="薄板格子 [m]")
    ap.add_argument("--l-wall", type=float, default=5.0e-2, help="前向き壁側の打ち切り深さ [m]")
    ap.add_argument("--l-top", type=float, default=5.5e-2, help="上面側の打ち切り (壁面から) [m]")
    ap.add_argument("--zmax", type=float, default=1.5e-2, help="z 方向の打ち切り [m]")
    ap.add_argument("--qfp-forge", type=float, default=77.34e3)
    ap.add_argument("--qfp-meas", type=float, default=63.44e3)
    a = ap.parse_args()

    steps = snaps(a.run)
    if not steps:
        raise SystemExit(f"{a.run}: res_gap_6_*.h5 が無い")
    if a.series:
        steps = [s for s in steps if s >= a.min_step]
    else:
        steps = [a.step if a.step is not None else steps[-1]]

    ref = load_ref(a.arm)
    off = 0 if a.arm == "L15" else 25
    series = []
    for step in steps:
        rows, info = evaluate(a, step)
        series.append((step, rows))
    step, rows = series[-1]

    print(f"[{a.run}] step {step}  腕 {a.arm}  薄板格子 {a.h*1e3:.2f} mm × {info['n_grid']}  "
          f"節点 {info['segn']}  z=0 節点 {info['z0_nodes']}")
    print(f"  τ={a.tau*1e3:.2f} mm  t_a={a.ta:.2f} s  Δ={a.dwin:.2f} s  "
          f"√(α t_a)={np.sqrt(ALPHA_S*a.ta)*1e3:.2f} mm  f_FP step {info['f_fp']['step']:.4f} / ramp {info['f_fp']['ramp']:.4f}")
    print("  q̂ は f_FP で割った値 (分母と同じ演算)。帯 = run8/run14 の反復区間 (未丸めの台帳値)。\n")
    print("   TC  位置            ℓ[mm]    生 q₀[kW/m²]  q̂ step   q̂ ramp  | q̂/q_FP(実測分母) step  ramp  | forge 分母 step | 実測帯        判定(step/ramp)")
    for r_ in rows:
        tc = r_["tc"] + off
        band = ref.get(r_["tc"])
        lab = {"wall": "壁", "arc": "円弧", "top": "上面"}[r_["kind"]]
        ns, nr = r_["step"] / a.qfp_meas, r_["ramp"] / a.qfp_meas
        nf = r_["step"] / a.qfp_forge
        if band:
            lo, hi = min(band), max(band)
            judge = "/".join("IN" if lo <= x <= hi else ("LOW" if x < lo else "HIGH") for x in (ns, nr))
            bs = f"{lo:.2f}–{hi:.2f}"
        else:
            judge, bs = "—", "(実測なし)"
        print(f"  {tc:3d}  {lab:4s}          {r_['ell']*1e3:7.2f}   {r_['raw']/1e3:10.2f}  "
              f"{r_['step']/1e3:8.2f} {r_['ramp']/1e3:8.2f}  |   {ns:7.3f}  {nr:7.3f}   |   {nf:7.3f}      | {bs:12s}  {judge}")
        if "step_rng" in r_:
            s0, s1 = (x / a.qfp_meas for x in r_["step_rng"])
            print(f"        (位置 ±0.25 cm の区間: step {s0:.3f}–{s1:.3f}、ramp "
                  f"{r_['ramp_rng'][0]/a.qfp_meas:.3f}–{r_['ramp_rng'][1]/a.qfp_meas:.3f})")

    # D2: 前向き壁の単調性 (深い→浅いで増加。隣接比 1.2 倍以内の反転は許容)
    wall = sorted([r_ for r_ in rows if r_["kind"] == "wall"], key=lambda x: x["ell"])
    for k in ("step", "ramp"):
        v = [x[k] for x in wall]
        ok = all(v[i + 1] >= v[i] / 1.2 for i in range(len(v) - 1))
        print(f"  D2 (前向き壁の単調性, {k}): {'PASS' if ok else 'FAIL'}  "
              + " ".join(f"{x['tc']+off}:{x[k]/a.qfp_meas:.3f}" for x in wall))

    if a.csv:
        tcs = [r_["tc"] for r_ in series[-1][1]]
        cols = ["step"] + [f"qs_tc{t+off}" for t in tcs] + [f"qr_tc{t+off}" for t in tcs] + [f"q0_tc{t+off}" for t in tcs]
        data = []
        for s, rr in series:
            d = {x["tc"]: x for x in rr}
            data.append([s] + [d[t]["step"] for t in tcs] + [d[t]["ramp"] for t in tcs] + [d[t]["raw"] for t in tcs])
        np.savetxt(a.csv, np.array(data), delimiter=",", header=",".join(cols), comments="", fmt="%.8g")
        print(f"\n  -> {a.csv} ({len(data)} スナップショット)")


if __name__ == "__main__":
    main()
