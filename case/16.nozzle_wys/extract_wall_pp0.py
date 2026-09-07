#!/usr/bin/env python3
"""node res_*.h5 から壁圧 p_wall/p0 (x 方向分布) を抽出する (2D 平面 / 3D 半幅対称面 共通)。
  contour: 輪郭壁 (wall_dist=0, y>0, 3D は対称面 z=z_max 上) の各 x 列の壁ノード
  side   : 3D のみ、側壁 (z=0) の中央線 (各 x 列で |y| 最小の壁ノード)
  center : 対称面/中心の |y| 最小ノード (参考)
usage: extract_wall_pp0.py RES.h5 OUT.csv [p0=59070]"""
import sys, h5py, numpy as np
fn, out = sys.argv[1], sys.argv[2]; p0 = float(sys.argv[3]) if len(sys.argv) > 3 else 59070.0
f = h5py.File(fn, "r"); c = np.array(f["MESH/COORD"]).reshape(-1, 3); P = np.array(f["VALUE/P"]); wd = np.array(f["VALUE/wall_dist"])
G = np.array(f["VALUE/g_0"]) if "VALUE/g_0" in f else None      # 凝縮液滴質量分率 (あれば)
TT = np.array(f["VALUE/T"])
zmax = c[:, 2].max(); is3d = zmax > 1e-6
xs = np.round(c[:, 0], 7); ux = np.unique(xs)
def line(mask, pick, arr=None):
    xx, pp = [], []
    for xv in ux:
        cand = np.where(mask & (xs == xv))[0]
        if len(cand) == 0: continue
        j = pick(cand); xx.append(c[j, 0] * 1e3); pp.append((P[j] / p0) if arr is None else arr[j])
    return np.array(xx), np.array(pp)
plane = (np.abs(c[:, 2] - zmax) < 1e-7) if is3d else np.ones(len(c), bool)
cx, cp = line(plane & (wd <= 0) & (c[:, 1] > 0), lambda cand: cand[np.argmax(c[cand, 1])])
mx, mp = line(plane, lambda cand: cand[np.argmin(np.abs(c[cand, 1]))])
cols = {"x_mm": cx, "contour_wall": cp, "center": np.interp(cx, mx, mp)}
pickc = lambda cand: cand[np.argmin(np.abs(c[cand, 1]))]
cols["center_T"] = np.interp(cx, *line(plane, pickc, TT))
if G is not None:
    cols["center_g"] = np.interp(cx, *line(plane, pickc, G))
    cols["contour_g"] = line(plane & (wd <= 0) & (c[:, 1] > 0), lambda cand: cand[np.argmax(c[cand, 1])], G)[1]
if is3d:
    sx, sp = line((np.abs(c[:, 2]) < 1e-9), lambda cand: cand[np.argmin(np.abs(c[cand, 1]))])
    cols["side_wall_mid"] = np.interp(cx, sx, sp)
hdr = ",".join(cols); np.savetxt(out, np.column_stack(list(cols.values())), delimiter=",", header=hdr, comments="")
print("wrote", out, "3D" if is3d else "2D", len(cx), "stations")
