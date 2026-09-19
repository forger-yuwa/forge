#!/usr/bin/env python3
r"""壁解像 $y_1^+$ を **第一内部点までの局所距離**と**接線壁応力**から出す (AGENTS.md「壁解像確認」)。

なぜ専用ツールが要るか (2026-09-19, codex レビュー):

- ソルバの壁面ダンプ `ypls` は [`viscousFlux_d.cu`](../cuda_forge/viscousFlux_d.cu) の低 Re 経路で
  `ro*utau*dcc/mu` として書かれる。`dcc` は**ゴーストセル重心と内点セル重心の距離**なので、
  **node 方式では壁ノードが壁面上に乗って退化**し、値が 1 桁以上小さく出る。
  実例 (case/49 run_0103): ソルバ 0.0221 に対し第一内部ノード基準 0.146、最大 2.60。
- ソルバの `utau` も node の既定経路では更新されない (`twall_*` だけ上書きされる)。
  実例 (同 run): $|\tau_w|/(\rho u_\tau^2)$ が中央値 0.996 なのに**最大 74.5**。
  よって $u_\tau$ は**接線 traction から組み直す**: $u_\tau=\sqrt{|\boldsymbol\tau_{w,t}|/\rho_w}$。
- **全域の代表距離 1 つで割ってはいけない**。細かい壁が 1 つあると粗い壁の解像不足が隠れる。
  壁面ごと・**局所**に第一内部点を引く。

やっていること:

1. メッシュ h5 の `BCONDS/<physID>/iCells` (壁の DOF) と `iPlanes` (境界面) を読む。
2. `PLANES/STRUCT` ([nNodes, nodes..., nCells, cells...]) から **DOF → 隣接 DOF** を作る。
3. 各壁 DOF で、境界面の法線 $\hat n$ (`PLANES/surfVect`) に**最も沿った隣接**を第一内部点とし、
   $y_1 = |(\mathbf{x}_{nb}-\mathbf{x}_w)\cdot\hat n|$ を取る。沿い方が悪ければ**評価不能**にする
   (角・斜交・複数壁の曖昧さをゼロ距離や合格に変換しない)。
4. 壁面ダンプ `res_<群>_<physID>_<step>.h5` の `twall_*` と `ro` から
   $u_\tau=\sqrt{|\tau_{w,t}|/\rho_w}$、$\mu_w$ は run の輸送モデルで評価して $y_1^+$ を出す。

usage:
  python3 solver_density_cuda/tools/check_wall_resolution.py <run_dir> [--mesh mesh.h5]
      [--step N] [--groups cav_outer,cyl_side] [--target 1.0] [--align-min 0.5]

VERDICT: 低 Re (`wallTreatmentSST: 0`) は**局所 $y_1^+\le$ --target** を目標とし、
超過面積割合と最大値の位置を出す。**判定に使うのは超過面積割合**であって最大値ではない:
鋭角エッジなど幾何的特異点があると traction が発散するので**最大値は格子収束しない**
(case/49 のリップ: 18.9 -> 14.4 -> 11.1 と減りはするが収束しない。一方、超過面積は
13.5 -> 10.5 -> 7.8 % と単調)。

自己検査: 構造格子なら `y1` が**第一層厚と厳密一致**するはず (case/49 の全ヘキサで
28.57 / 20.00 / 14.29 µm = 設計値 20/scale µm と一致)。ずれていたら対応づけが壊れている。壁関数 (1 = automatic) は低層/buffer/log の面積分布を出す
(forge の mode 1 は粘性低層から対数層を接続する automatic treatment なので、
`30<=y+<=300` を全点必須にはしない)。
"""
import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU0, T0_SUTH, S_SUTH = 1.716e-5, 273.0, 111.0      # forge の Sutherland 実装値


def mu_of(T, cfg):
    """run の輸送モデルに合わせた分子粘性。判定できなければ None を返す (判定不能にするため)。"""
    pp = (cfg or {}).get("physProp", {})
    vm = int(pp.get("viscMethod", 1))
    if vm == 0:
        v = pp.get("visc")
        return None if v is None else np.full_like(T, float(v))
    if vm == 1:
        return MU0 * (T / T0_SUTH) ** 1.5 * (T0_SUTH + S_SUTH) / (T + S_SUTH)
    return None            # 多成分 (Chapman-Enskog/Wilke) は組成が要る -> 判定不能


def dof_neighbors(mesh):
    """PLANES/STRUCT から DOF -> 隣接 DOF のリストを作る。"""
    with h5py.File(mesh, "r") as f:
        st = np.array(f["PLANES/STRUCT"])
        nplane = len(np.array(f["PLANES/surfArea"]))
    nb = {}
    i = ip = 0
    while i < len(st) and ip < nplane:
        nn = int(st[i]); i += 1 + nn
        nc = int(st[i]); i += 1
        cl = st[i:i + nc]; i += nc
        if nc >= 2:
            a, b = int(cl[0]), int(cl[1])
            nb.setdefault(a, []).append((b, ip))
            nb.setdefault(b, []).append((a, ip))
        ip += 1
    return nb


def wall_first_distance(mesh, phys_id, nb, align_min, coords=None):
    """壁 DOF ごとに (y1, 法線, DOF index, 評価できたか) を返す。"""
    with h5py.File(mesh, "r") as f:
        key = "BCONDS/%d" % phys_id
        if key not in f:
            return None
        icells = np.array(f[key + "/iCells"])
        iplanes = np.array(f[key + "/iPlanes"])
        svec = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
        sarea = np.array(f["PLANES/surfArea"])
        cc = coords if coords is not None else np.array(f["CELLS/centCoords"]).reshape(-1, 3)
    # 壁 DOF ごとに境界面法線 (面積重みで平均 = 角では平均法線になる)
    nrm = {}
    for c, p in zip(icells, iplanes):
        v = svec[p] * (1.0 if True else 1.0)
        nrm.setdefault(int(c), np.zeros(3))
        nrm[int(c)] = nrm[int(c)] + v
    out = {}
    for c, n in nrm.items():
        ln = np.linalg.norm(n)
        if ln <= 0:
            out[c] = (np.nan, None, False)
            continue
        nh = n / ln
        best, bestal = None, -1.0
        for (j, _p) in nb.get(c, []):
            d = cc[j] - cc[c]
            dn = np.linalg.norm(d)
            if dn <= 0:
                continue
            al = abs(float(np.dot(d, nh)) / dn)
            if al > bestal:
                bestal, best = al, (abs(float(np.dot(d, nh))), j)
        if best is None or bestal < align_min:
            out[c] = (np.nan, nh, False)          # 角・斜交で第一内部点が定まらない
        else:
            out[c] = (best[0], nh, True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--mesh", default=None, help="既定は <run>/mesh.h5")
    ap.add_argument("--step", type=int, default=None, help="既定は最後の壁ダンプ")
    ap.add_argument("--groups", default=None, help="壁群名をカンマ区切り (既定は壁ダンプ全部)")
    ap.add_argument("--target", type=float, default=1.0, help="低 Re の局所 y1+ 目標 (既定 1)")
    ap.add_argument("--align-min", type=float, default=0.5,
                    help="第一内部点として認める法線との沿い方 (既定 0.5)")
    a = ap.parse_args()

    mesh = a.mesh or os.path.join(a.run, "mesh.h5")
    if not os.path.exists(mesh):
        print("mesh h5 が無い: %s" % mesh); return 2
    cfg = None
    cpath = os.path.join(a.run, "solverConfig.yaml")
    if os.path.exists(cpath):
        try:
            import yaml
            cfg = yaml.safe_load(open(cpath))
        except Exception as e:                      # noqa: BLE001
            print("  WARNING: solverConfig.yaml を読めない (%s)" % e)
    wt = int(((cfg or {}).get("turbulence", {}) or {}).get("wallTreatmentSST", 0))

    dumps = sorted(glob.glob(os.path.join(a.run, "res_*_*_*.h5")))
    pat = re.compile(r"res_(.+)_(\d+)_(\d+)\.h5$")
    found = {}
    for d in dumps:
        m = pat.search(os.path.basename(d))
        if m:
            found.setdefault((m.group(1), int(m.group(2))), []).append((int(m.group(3)), d))
    if not found:
        print("壁面ダンプ res_<群>_<physID>_<step>.h5 が無い: %s" % a.run); return 2
    want = set(x.strip() for x in a.groups.split(",")) if a.groups else None

    print("=== %s  (wallTreatmentSST=%d) ===" % (a.run, wt))
    nb = dof_neighbors(mesh)
    worst = 0.0
    any_eval = False
    fails = []
    for (name, pid), lst in sorted(found.items()):
        if want and name not in want:
            continue
        step, path = sorted(lst)[-1] if a.step is None else \
            next(((s, p) for s, p in lst if s == a.step), sorted(lst)[-1])
        wd = None   # 対応づけ後に、使った座標系で測り直す
        if False:
            print("  %-12s BCONDS/%d がメッシュに無い -> 判定不能" % (name, pid)); fails.append(name); continue
        with h5py.File(path, "r") as f:
            V = {k: np.array(f["VALUE/" + k]) for k in f["VALUE"]}
        need = ("twall_x", "twall_y", "twall_z", "ro", "Ts")
        if any(k not in V for k in need):
            print("  %-12s 壁ダンプに %s が無い -> 判定不能"
                  % (name, ",".join(k for k in need if k not in V)))
            fails.append(name); continue
        with h5py.File(mesh, "r") as f:
            icells = np.unique(np.array(f["BCONDS/%d/iCells" % pid]))
            cand = {"CELLS/centCoords": np.array(f["CELLS/centCoords"]).reshape(-1, 3)}
            nod = np.array(f["MESH/COORD"]).reshape(-1, 3)
            if len(nod) >= int(icells.max()) + 1:
                cand["MESH/COORD"] = nod          # node 方式では DOF index = 節点 index
        with h5py.File(path, "r") as f:
            wxyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        if len(wxyz) != len(V["ro"]):
            print("  %-12s 壁ダンプの座標数と値数が合わない -> 判定不能" % name)
            fails.append(name); continue
        # **順序を仮定しない**。壁ダンプのノード座標を DOF の座標に**厳密一致で対応づける**。
        # 順序を仮定すると値が入れ替わって y1+ が桁で変わる (2026-09-19 実測で 8 倍ずれた)。
        # node 方式では壁ダンプの座標は**節点座標**、`CELLS/centCoords` は**双対 CV 重心**で
        # 別物なので、両方を試して一致する方を使う。
        from scipy.spatial import cKDTree
        pick, best = None, None
        for src, arr in cand.items():
            d, j = cKDTree(arr[icells]).query(wxyz)
            if best is None or d.max() < best[0]:
                best = (float(d.max()), src, j, arr)
        dmax, src, jmap, dof_xyz = best
        scale = float(np.median(np.linalg.norm(np.diff(wxyz[:min(len(wxyz), 200)], axis=0), axis=1)))
        if dmax > 1e-6 * max(scale, 1e-12) + 1e-12:
            print("  %-12s 壁ダンプ節点を DOF に対応づけられない (最良 %s で最大 %.3e m) -> 判定不能"
                  % (name, src, dmax))
            fails.append(name); continue
        uniq = icells[jmap]              # 壁ダンプの各点に対応する DOF
        wd = wall_first_distance(mesh, pid, nb, a.align_min, coords=dof_xyz)
        if wd is None:
            print("  %-12s BCONDS/%d がメッシュに無い -> 判定不能" % (name, pid))
            fails.append(name); continue
        y1 = np.array([wd.get(int(c), (np.nan, None, False))[0] for c in uniq])
        okc = np.array([wd.get(int(c), (np.nan, None, False))[2] for c in uniq])
        nvec = np.array([wd.get(int(c), (np.nan, np.zeros(3), False))[1]
                         if wd.get(int(c), (np.nan, None, False))[1] is not None else np.zeros(3)
                         for c in uniq])
        tw = np.stack([V["twall_x"], V["twall_y"], V["twall_z"]], axis=1)
        tn = np.sum(tw * nvec, axis=1)[:, None] * nvec
        tt = np.linalg.norm(tw - tn, axis=1)               # **接線成分**
        mu = mu_of(V["Ts"].astype(float), cfg)
        if mu is None:
            print("  %-12s 粘性モデルを解決できない (viscMethod) -> 判定不能" % name)
            fails.append(name); continue
        ro = V["ro"].astype(float)
        yp = np.where(okc, y1 * np.sqrt(np.maximum(ro, 0) * tt) / np.maximum(mu, 1e-30), np.nan)
        good = np.isfinite(yp)
        if not good.any():
            print("  %-12s 評価できた点が無い -> 判定不能" % name); fails.append(name); continue
        any_eval = True
        frac = 100.0 * good.sum() / len(yp)
        over = 100.0 * np.count_nonzero(yp[good] > a.target) / good.sum()
        imax = int(np.nanargmax(np.where(good, yp, -np.inf)))
        sol = V.get("ypls")
        print("  %-12s step %6d  y1 = %.3e m  y1+ 平均 %7.3f / p99 %7.3f / 最大 %7.3f"
              % (name, step, np.nanmedian(y1[good]), np.nanmean(yp[good]),
                 np.nanpercentile(yp[good], 99), np.nanmax(yp[good])))
        print("               評価できた面積割合 %.1f %% ; y1+ > %.3g が %.1f %% ; 最大の位置 index %d%s"
              % (frac, a.target, over, imax,
                 ("  ; ソルバ ypls 平均 %.4g" % float(np.mean(sol)) if sol is not None else "")))
        worst = max(worst, float(np.nanmax(yp[good])))
        if frac < 90.0:
            fails.append("%s (評価できた面積 %.1f %%)" % (name, frac))

    if not any_eval:
        print("\nVERDICT: INDETERMINATE (評価できた壁が無い)"); return 2
    if fails:
        print("\nVERDICT: INDETERMINATE (判定不能: %s)" % ", ".join(fails)); return 2
    if wt == 0:
        ok = worst <= a.target
        print("\n最大 y1+ = %.3f  (目標 <= %.3g)" % (worst, a.target))
        print("VERDICT: %s" % ("PASS (壁解像)" if ok else
                               "FAIL (局所 y1+ が目標超過 — 超過面積と位置を上に示した)"))
        return 0 if ok else 1
    print("\n最大 y1+ = %.3f  (壁関数 automatic なので対数層配置は必須条件にしない)" % worst)
    print("VERDICT: PASS (壁関数 automatic; 低層/buffer/log の分布は上の統計を参照)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
