#!/usr/bin/env python3
"""全周 (360°) 計算の左右対称性を測る — 半割 (y=0 対称面) で計算してよいかのゲート (plan §4.9)。

キャビティ内の気体ノードを y -> -y で鏡像に写し、**最近傍の相手と場を比べる**。
偏心は x 方向だけなので幾何は y=0 対称であり、**対称な解が存在する**。問われるのは
その解が**安定か** (非対称モードが育たないか) なので、
  ① 定常 run が対称に留まっているか (擾乱を入れなくても離散誤差で非対称モードが育つなら失敗)
  ② 非対称擾乱を入れた URANS で非対称が**減衰するか、成長するか**
の 2 つを見る。①→②の順に効かせる。

usage:
  python3 tools/check_symmetry.py <run_dir> [...] [--tol 0.02] [--series]

指標は **方位ビン平均の左右差**を主に使う:
  A(theta) = 帯 (r, z) ごとの theta ビン平均       ->  asym_bin = rms(A(+t) - A(-t)) / 振幅
節点ごとの鏡像最近傍差 (asym_rms) も出すが、**非構造メッシュでは鏡像の相手が遠く
(テトラ全周で最大 1.65 mm = すきま幅の 2/3)、メッシュの非対称が指標に混ざるので
判定には使わない** (診断のみ。鏡像最近傍距離を必ず併記する)。ビン平均は節点の散らばりを
平均してしまうのでメッシュ非対称に鈍い。判定は
  asym_bin <= --tol かつ (時系列なら) 末尾で増加していない -> PASS
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geom_common as gc                      # noqa: E402

WALLS = ("cav_outer", "cyl_side", "cav_floor")


def snapshots(run):
    return sorted(Path(run).glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))


def mirror_pairs(xyz):
    """y -> -y の鏡像相手を最近傍で引く。戻り値 (j, dist)。"""
    m = xyz.copy()
    m[:, 1] = -m[:, 1]
    d, j = cKDTree(xyz).query(m)
    return j, d


def asym_bins(xyz, vals, nb=36, zb=5):
    """方位ビン平均の左右差。theta は y>=0 側を 0..pi に取り、-y 側を折り返して比べる。
    深さ方向 zb 段に分けてから方位 nb ビンで平均する (深さ方向の勾配で薄まらないように)。"""
    th = np.arctan2(np.abs(xyz[:, 1]), -xyz[:, 0])      # theta=0 が上流 (-x)
    side = np.sign(xyz[:, 1])
    z = xyz[:, 2]
    # z 方向にほとんど広がりの無い面 (床) を 5 段に切ると 1 ビンあたりの節点が足りなくなる
    rext = float(np.max(np.hypot(xyz[:, 0], xyz[:, 1])) - np.min(np.hypot(xyz[:, 0], xyz[:, 1])))
    if float(z.max() - z.min()) < 0.2 * max(rext, 1e-12):
        zb = 1
    zi = np.clip(((z - z.min()) / max(z.max() - z.min(), 1e-12) * zb).astype(int), 0, zb - 1)
    ti = np.clip((th / np.pi * nb).astype(int), 0, nb - 1)
    amp = float(np.max(vals) - np.min(vals))
    if amp <= 0:
        return 0.0, 0.0, amp, 0
    d, n = [], 0
    for k in range(zb):
        for t in range(nb):
            m = (zi == k) & (ti == t)
            a, b = m & (side > 0), m & (side < 0)
            if a.sum() < 3 or b.sum() < 3:
                continue
            d.append(vals[a].mean() - vals[b].mean()); n += 1
    if not d:
        return float("nan"), float("nan"), amp, 0
    d = np.asarray(d)
    return float(np.sqrt(np.mean(d ** 2)) / amp), float(np.max(np.abs(d)) / amp), amp, n


def asym(vals, j):
    a = vals - vals[j]
    amp = float(np.max(vals) - np.min(vals))
    if amp <= 0:
        return 0.0, 0.0, amp
    return float(np.sqrt(np.mean(a ** 2)) / amp), float(np.max(np.abs(a)) / amp), amp


def field_asym(res, man):
    with h5py.File(res, "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        T = np.array(f["VALUE/T"])
        U = [np.array(f["VALUE/U%s" % k]) for k in ("x", "y", "z")]
    m = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man, shrink=man["eval"]["shrink_m"])
    if m.sum() < 100:
        return None
    cc = c[m]
    j, d = mirror_pairs(cc)
    out = {"n_cav": int(m.sum()), "pair_dist_max": float(d.max())}
    out["T"] = asym(T[m], j)
    out["bin"] = {"T": asym_bins(cc, T[m]), "Ux": asym_bins(cc, U[0][m]),
                  "Uz": asym_bins(cc, U[2][m])}
    # **この指標のメッシュ由来の床**: 幾何だけで決まる (厳密に左右対称な) 場を同じ指標にかける。
    # 非構造メッシュでは節点の並びが左右で違うので、ビン平均にも有限の差が残る。
    # 物理の非対称がこの床と同程度なら**判定不能**であって「非対称」ではない。
    floor = {"radius": asym_bins(cc, np.hypot(cc[:, 0], cc[:, 1]))}
    with h5py.File(res, "r") as f:
        if "VALUE/wall_dist" in f:
            floor["wall_dist"] = asym_bins(cc, np.array(f["VALUE/wall_dist"])[m])
    out["floor"] = floor
    # **旋回 (周方向循環)**: 鏡像で符号が反転するので、対称解では厳密に 0。
    # メッシュの左右非対称では符号が決まらないため**この指標は汚染されない**。
    # 半割の可否は「旋回擾乱が減衰して 0 に戻るか」で決まる (plan §4.9)。
    r = np.hypot(cc[:, 0], cc[:, 1])
    ut = (-cc[:, 1] * U[0][m] + cc[:, 0] * U[1][m]) / np.maximum(r, 1e-12)
    uref = float(np.sqrt(np.mean(U[0][m] ** 2 + U[1][m] ** 2 + U[2][m] ** 2)))
    out["swirl"] = (float(np.mean(ut) / max(uref, 1e-30)),
                    float(np.sqrt(np.mean(ut ** 2)) / max(uref, 1e-30)), uref)
    # Uy は鏡像で符号が反転する量なので、+ で比べる (反対称成分の残りを見る)
    uy = U[1][m]
    amp = float(np.max(np.abs(uy))) or 1.0
    out["Uy_anti"] = (float(np.sqrt(np.mean((uy + uy[j]) ** 2)) / amp),
                      float(np.max(np.abs(uy + uy[j])) / amp), amp)
    for nm, k in (("Ux", 0), ("Uz", 2)):
        out[nm] = asym(U[k][m], j)
    return out


def wall_asym(run, step, man):
    out = {}
    for g in WALLS:
        pid = man["phys_id"].get(g)
        p = Path(run) / ("res_%s_%d_%d.h5" % (g, pid, step))
        if not p.exists():
            continue
        with h5py.File(p, "r") as f:
            c = np.array(f["MESH/COORD"]).reshape(-1, 3)
            if "qwall" not in f["VALUE"]:
                continue
            q = np.array(f["VALUE/qwall"])
        j, _ = mirror_pairs(c)
        out[g] = (asym(q, j), asym_bins(c, q))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tol", type=float, default=0.02, help="asym_rms の許容 (既定 2 %)")
    ap.add_argument("--series", action="store_true", help="全スナップショットで時系列を見る")
    ap.add_argument("--ref-run", default=None,
                    help="擾乱なしの定常 run。旋回の**減衰率**で判定する (これが本来の判定)")
    ap.add_argument("--decay", type=float, default=0.2,
                    help="残った旋回 / 初期擾乱 の許容 (既定 0.2 = 8 割減衰していれば PASS)")
    a = ap.parse_args()
    man = gc.load_manifest(run=a.runs[0])
    if man["geometry"].get("half_model", True):
        print("注意: manifest が半割 (half_model: true)。全周 run には CASE49_MANIFEST=manifest_full.json")
    worst, grow, floor_lo, floor_hi = 0.0, False, float("inf"), 0.0
    swirl_hist = []
    for run in a.runs:
        snaps = snapshots(run)
        print("\n=== %s  (%d スナップショット) ===" % (run, len(snaps)))
        if not snaps:
            print("  res_*.h5 が無い"); return 1
        use = snaps if a.series else snaps[-1:]
        hist, swirl_hist = [], []
        for res in use:
            step = int(res.stem.split("_")[1])
            fa = field_asym(res, man)
            if fa is None:
                print("  step %6d  キャビティ節点が取れない" % step); continue
            wa = wall_asym(run, step, man)
            print("  step %6d  n_cav %d  鏡像最近傍距離 max %.2e m" % (step, fa["n_cav"], fa["pair_dist_max"]))
            for k in ("T", "Ux", "Uz", "Uy_anti"):
                r, mx, amp = fa[k]
                extra = ""
                if k in fa["bin"]:
                    br, bmx, _, nb = fa["bin"][k]
                    extra = "   |  ビン asym %8.2e (max %8.2e, %d ビン)" % (br, bmx, nb)
                    worst = max(worst, br)
                print("     %-8s 最近傍 asym %8.2e (max %8.2e, 振幅 %.4g)%s"
                      % (k, r, mx, amp, extra))
            for g, ((r, mx, amp), (br, bmx, _, nb)) in wa.items():
                print("     q''[%-10s] 最近傍 asym %8.2e  |  ビン asym %8.2e (max %8.2e, %d ビン, 振幅 %.4g W/m2)"
                      % (g, r, br, bmx, nb, amp))
                worst = max(worst, br)
            sw, swr, uref = fa["swirl"]
            print("     旋回      平均 u_theta/u_rms %+9.3e   rms u_theta/u_rms %8.3e  (u_rms %.4g m/s)"
                  % (sw, swr, uref))
            for k, (r, mx, amp, nb) in fa["floor"].items():
                print("     [床] %-10s ビン asym %8.2e (max %8.2e)  <- 幾何のみ = メッシュ由来"
                      % (k, r, mx))
            fl = max(v[0] for v in fa["floor"].values())
            floor_hi = max(floor_hi, fl)
            floor_lo = min(floor_lo, min(v[0] for v in fa["floor"].values()))
            print("     -> 物理の最大ビン asym / メッシュ床 = %.2f 倍"
                  % (max(fa["bin"][k][0] for k in fa["bin"]) / max(fl, 1e-30)))
            hist.append((step, abs(sw)))
            swirl_hist.append((step, sw))
        if len(hist) >= 3:
            h = [v for _, v in hist]
            n = max(2, len(h) // 3)
            grow = grow or (np.mean(h[-n:]) > 1.3 * np.mean(h[:n]))
            print("  |平均旋回|: 前期 %.2e -> 末期 %.2e  (%s)"
                  % (np.mean(h[:n]), np.mean(h[-n:]), "成長" if grow else "減衰/横ばい"))
    # --- 旋回の減衰率による判定 (主) ---
    if a.ref_run and swirl_hist:
        rs = field_asym(snapshots(a.ref_run)[-1], man)
        sref = rs["swirl"][0] if rs else 0.0
        steps = np.array([st for st, _ in swirl_hist], float)
        dev = np.abs(np.array([v for _, v in swirl_hist], float) - sref)
        d0 = dev[0]
        # **終点 1 点で判定しない** (codex result-1 M8, 2026-09-20)。
        # 減衰したあと再成長しても、振動して偶然終点が基準近傍でも合格してしまっていた。
        # 末尾窓の**最大**で測り、さらに末尾が上昇していないかを別に見る。
        nt = max(2, len(dev) // 3)
        tail_max = float(dev[-nt:].max())
        tail_mean = float(dev[-nt:].mean())
        mid_mean = float(dev[-2 * nt:-nt].mean()) if len(dev) >= 2 * nt else tail_mean
        frac = tail_max / max(d0, 1e-30)
        rising = tail_mean > 1.10 * mid_mean
        # 指数減衰率 (e-folding step)。遅い減衰と速い減衰を区別する
        efold = float("nan")
        ok = dev > 0
        if ok.sum() >= 3:
            sl = np.polyfit(steps[ok], np.log(dev[ok]), 1)[0]
            if sl < 0:
                efold = -1.0 / sl
        print("\n旋回の減衰: 擾乱 %+.3e -> 末尾窓 最大 %.3e / 平均 %.3e   (基準 %+.3e)"
              % (swirl_hist[0][1], tail_max, tail_mean, sref))
        print("  |末尾窓の最大-基準| / |擾乱-基準| = %.4f   (許容 %.3g)   末尾 %d 点"
              % (frac, a.decay, nt))
        print("  末尾窓平均 / 中間窓平均 = %.3f  (%s)   e-folding %s step"
              % (tail_mean / max(mid_mean, 1e-30), "**上昇**" if rising else "非上昇",
                 ("%.3g" % efold) if np.isfinite(efold) else "減衰していない"))
        if len(dev) < 3:
            # **1〜2 点で減衰は測れない**。旧実装はこの場合でも frac を出して判定していた
            # (`--series` を付け忘れると最終スナップショット 1 点だけになる)。
            print("VERDICT: INCONCLUSIVE (時系列が %d 点しかない — `--series` を付けて"
                  "スナップショット列で見ること)" % len(dev))
            return 2
        if d0 < 1e-6:
            print("VERDICT: NO-PERTURBATION (擾乱が基準と同じ — 擾乱が効いていない)")
            return 2
        if rising:
            print("VERDICT: ASYMMETRIC (末尾で再成長 — 終点だけ見ると見逃す)")
            return 1
        if frac >= 1.0:
            print("VERDICT: ASYMMETRIC (反対称モードが減衰しない/成長 — 半割は不可)")
            return 1
        if frac <= a.decay:
            print("VERDICT: SYMMETRIC (反対称モードが減衰 = 半割で可。減衰 %.1f %%、"
                  "末尾窓で再成長なし)" % (100 * (1 - frac)))
            return 0
        print("VERDICT: PARTIAL-DECAY (%.1f %% しか減衰していない — 観測窓を延ばす)" % (100 * (1 - frac)))
        return 2

    print("\n最大 ビン asym = %.3e  (許容 %.3g)   メッシュ床 %.3e 〜 %.3e%s"
          % (worst, a.tol, floor_lo, floor_hi, "   非対称が成長" if grow else ""))
    if grow:
        v, rc = "ASYMMETRIC (非対称が時間的に成長 — 半割は不可)", 1
    elif worst <= a.tol and worst <= floor_lo:
        v, rc = "SYMMETRIC (半割で可)", 0
    elif worst <= max(a.tol, floor_hi):
        # **非構造メッシュの節点の並びが左右で違うだけでこの値が出る**。「非対称」とは言えない。
        v, rc = ("INCONCLUSIVE (メッシュ床 %.1e〜%.1e の中 — 鏡像対称メッシュか URANS の減衰で判定せよ)"
                 % (floor_lo, floor_hi)), 2
    else:
        v, rc = "ASYMMETRIC (床より明確に大 — 半割は不可、全周で計算すること)", 1
    print("VERDICT: %s" % v)
    return rc


if __name__ == "__main__":
    sys.exit(main())
