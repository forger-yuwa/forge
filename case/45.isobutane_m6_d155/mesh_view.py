"""case/45 の格子を ParaView で見るためのファイルを作る (2026-10-10、ユーザ「y+ とアスペクト比を paraview で見たい」)。
入力: res_<N>.h5 (MESH と VALUE) と check_wall_resolution.py の壁の分布 CSV (y1p_<N>_wall.csv)。
構造格子 (列 i = 0..NI−1、j = 0 軸 .. NJ−1 壁) の一次の四角形を作り、次を書く。
  節点: y_plus (壁節点までの直線距離 × u_τ/ν_w、u_τ = √(τ_t/ρ_w)、ν_w = μ_w/ρ_w は列の壁節点の値)、y1plus_col (列の y₁⁺、ツールの値)、
        y_over_rt (壁からの距離 / r_t)、Mach、T
  セル: aspect_ratio (辺の長さの最大/最小)、dr_over_dx (半径方向の辺の平均 / 軸方向の辺の平均)、dx_over_rt・dr_over_rt
usage: python3 mesh_view.py <run_dir> <N> [--NI 4719 --NJ 121 --rt 0.07667]"""
import argparse, csv, h5py, numpy as np
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("N", type=int)
ap.add_argument("--NI", type=int, default=4719); ap.add_argument("--NJ", type=int, default=121); ap.add_argument("--rt", type=float, default=0.07667)
a = ap.parse_args(); run = Path(a.run); NI, NJ, RT = a.NI, a.NJ, a.rt
with h5py.File(run / f"res_{a.N}.h5") as h:
    X = np.asarray(h["MESH/COORD"][:], dtype=np.float64).reshape(-1, 3)
    V = {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64) for k in ("Ux", "Uy", "sonic", "T")}
assert X.shape[0] == NI * NJ
P = X.reshape(NI, NJ, 3)
# 壁の分布 (列ごと) を、壁節点の座標の最近傍で列に対応させる
rows = list(csv.DictReader(open(run / f"y1p_{a.N}_wall.csv")))
cx = np.array([[float(r["x"]), float(r["y"])] for r in rows])
wall = P[:, NJ - 1, :2]
order = np.lexsort((wall[:, 1], wall[:, 0])); corder = np.lexsort((cx[:, 1], cx[:, 0]))
col = np.empty(len(rows), dtype=int); col[corder] = order
err = np.max(np.linalg.norm(cx - wall[col], axis=1))
if err > 1e-9: raise SystemExit(f"壁節点と CSV の対応が取れない (最大のずれ {err:.3g} m)")
utau = np.empty(NI); nuw = np.empty(NI); y1p = np.empty(NI)
for r, i in zip(rows, col):
    rho, tau, mu = float(r["rho_w"]), float(r["tau_t"]), float(r["mu_w"])
    utau[i] = np.sqrt(abs(tau) / rho); nuw[i] = mu / rho; y1p[i] = float(r["y1plus"])
y = np.linalg.norm(P - P[:, NJ - 1:NJ, :], axis=2)                       # 列の壁節点までの直線距離
yplus = y * (utau / nuw)[:, None]
mach = np.hypot(V["Ux"], V["Uy"]) / np.maximum(V["sonic"], 1e-30)
# 一次の四角形 (i, j) - (i+1, j) - (i+1, j+1) - (i, j+1)
I, J = np.meshgrid(np.arange(NI - 1), np.arange(NJ - 1), indexing="ij")
n00 = (I * NJ + J).ravel(); quad = np.stack([n00, n00 + NJ, n00 + NJ + 1, n00 + 1], 1).astype(np.int32)
Pq = X[quad][:, :, :2]
ax1 = np.linalg.norm(Pq[:, 1] - Pq[:, 0], axis=1); ax2 = np.linalg.norm(Pq[:, 2] - Pq[:, 3], axis=1)
rd1 = np.linalg.norm(Pq[:, 3] - Pq[:, 0], axis=1); rd2 = np.linalg.norm(Pq[:, 2] - Pq[:, 1], axis=1)
edges = np.stack([ax1, ax2, rd1, rd2], 1)
ar = edges.max(1) / np.maximum(edges.min(1), 1e-300)
dx = 0.5 * (ax1 + ax2); dr = 0.5 * (rd1 + rd2)
out = run / f"meshview_{a.N}.h5"
with h5py.File(out, "w") as h:
    h["geom/xyz"] = X; h["topo/quad"] = quad
    h["node/y_plus"] = yplus.ravel(); h["node/y1plus_col"] = np.repeat(y1p, NJ); h["node/y_over_rt"] = (y / RT).ravel()
    h["node/Mach"] = mach; h["node/T"] = V["T"]
    h["cell/aspect_ratio"] = ar; h["cell/dr_over_dx"] = dr / dx; h["cell/dx_over_rt"] = dx / RT; h["cell/dr_over_rt"] = dr / RT
nn, nc = X.shape[0], quad.shape[0]
def attr(name, center, n):
    return (f'      <Attribute Name="{name}" AttributeType="Scalar" Center="{center}">\n'
            f'        <DataItem Dimensions="{n}" NumberType="Float" Precision="8" Format="HDF">{out.name}:/{"node" if center == "Node" else "cell"}/{name}</DataItem>\n      </Attribute>\n')
xmf = ('<?xml version="1.0" ?>\n<Xdmf Version="3.0">\n  <Domain>\n    <Grid Name="case45_mesh" GridType="Uniform">\n'
       f'      <Topology TopologyType="Quadrilateral" NumberOfElements="{nc}">\n        <DataItem Dimensions="{nc} 4" NumberType="Int" Precision="4" Format="HDF">{out.name}:/topo/quad</DataItem>\n      </Topology>\n'
       f'      <Geometry GeometryType="XYZ">\n        <DataItem Dimensions="{nn} 3" NumberType="Float" Precision="8" Format="HDF">{out.name}:/geom/xyz</DataItem>\n      </Geometry>\n'
       + "".join(attr(k, "Node", nn) for k in ("y_plus", "y1plus_col", "y_over_rt", "Mach", "T"))
       + "".join(attr(k, "Cell", nc) for k in ("aspect_ratio", "dr_over_dx", "dx_over_rt", "dr_over_rt"))
       + "    </Grid>\n  </Domain>\n</Xdmf>\n")
(run / f"meshview_{a.N}.xmf").write_text(xmf)
print(f"[mesh_view] {out} と {out.with_suffix('.xmf')} を書いた (節点 {nn}、セル {nc})")
print(f"  y1+ (列ごと、ツールの値): 平均 {y1p.mean():.3f}・最大 {y1p.max():.3f}")
print(f"  アスペクト比: 最大 {ar.max():.0f}・p99 {np.percentile(ar, 99):.0f}・平均 {ar.mean():.0f}")
xw = P[:-1, NJ - 1, 0] / RT
for xs in (-8, -5, -2, 0, 2, 10, 40, 94):
    i = int(np.argmin(np.abs(xw - xs))); c = i * (NJ - 1) + (NJ - 2); cax = i * (NJ - 1)
    print(f"  x/r_t {xs:4d}: 軸の列のセル dx {dx[cax]/RT:.4f}・dr {dr[cax]/RT:.3f} r_t (dr/dx {dr[cax]/dx[cax]:.1f})、壁の列のセル dx {dx[c]/RT:.4f}・dr {dr[c]/RT:.2e} r_t (AR {ar[c]:.0f})")
