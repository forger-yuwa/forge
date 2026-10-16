# 出口の列 (i = ni−1) と途中の列の M(r) を構造格子の index (i·nj + j、j = 0 軸) で読む
import sys, h5py, numpy as np
run, step, ni, nj = sys.argv[1], sys.argv[2], 4810, 92
with h5py.File(f"{run}/res_{step}.h5") as f:
    V = f["VALUE"]; U = np.hypot(V["Ux"][:], V["Uy"][:]); a = V["sonic"][:]; M = U / a
    C = f["MESH/COORD"][:] if "MESH/COORD" in f else None
if C is None:
    with h5py.File(f"{run}/nozzle.h5") as g: C = g["MESH/COORD"][:]
C = C.reshape(-1, 3) if C.ndim == 1 else C
for i in (ni - 1, ni - 400, ni - 1200):
    idx = i * nj + np.arange(nj); r = C[idx, 1]; x = C[idx, 0]; m = M[idx]
    Max = m[0]; rw = r[-1]
    j99 = np.argmax(m < 0.99 * Max) if np.any(m < 0.99 * Max) else nj
    jdev = np.argmax(np.abs(m - Max) > 0.005 * Max)
    core = slice(0, j99)
    print(f"i {i} x {x[0]:.4f} m  r_w {rw:.4f} m  M_axis {Max:.5f}  δ99 (M<0.99 M_axis) from r {r[j99]:.4f} m → BL {rw - r[j99]:.4f} m"
          f"  |M−M_axis|>0.5 % from r {r[jdev]:.4f} m ({r[jdev]/rw:.3f} r_w)")
    for fr in (0.5, 0.7, 0.8, 0.85, 0.9):
        jj = np.searchsorted(r, fr * rw); print(f"    r/r_w {fr:.2f}: M {m[jj]:.5f}")
    c = r <= 0.8 * rw
    print(f"    r ≤ 0.8 r_w: M mean {m[c].mean():.5f}  min {m[c].min():.5f}  max {m[c].max():.5f}  (max−min)/mean {(m[c].max()-m[c].min())/m[c].mean()*100:.3f} %")
