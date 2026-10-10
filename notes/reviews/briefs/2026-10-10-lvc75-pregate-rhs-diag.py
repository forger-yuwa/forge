# rhs_s0 の外れの中身 (記録用): どの節点・行か、値の大きさと float の ulp、rhs_s0 = float(残差) が各書き出しで成り立つか
import numpy as np, h5py
runs = {"m7": "run_0574_lvc75_m7_dump", "m7b": "run_0576_lvc75_m7_dump2", "m5": "run_0575_lvc75_m5_dump"}
def A(r, n):
    d = r + "/linedump"
    meta = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in open(d + "/meta.txt") if l.strip() and not l.startswith("#")}
    return np.fromfile(d + "/" + n + ".f64").reshape(*meta[n])
rhs = {k: A(r, "rhs_s0") for k, r in runs.items()}
nl = A(runs["m7"], "node_line"); idx = nl[:, 0].astype(int); fl = A(runs["m7"], "flags_wall_iso_axis")
for a, b in (("m5", "m7"), ("m7b", "m7")):
    d = np.abs(rhs[a] - rhs[b]); k, r = np.unravel_index(np.argmax(d), d.shape)
    v = rhs[b][k, r]
    print(f"{a}−{b}: 最大 {d[k, r]:.3e} 節点 {idx[k]} 行 {r} 値 {v:.6e} ulp32 {np.spacing(np.float32(abs(v))):.3e} → {d[k, r] / np.spacing(np.float32(abs(v))):.1f} ulp; "
          f"非零の差の要素 {int(np.sum(d > 0))}; 差/ulp32 の最大 {float(np.max(d / np.maximum(np.spacing(np.abs(rhs[b]).astype(np.float32)).astype(float), 1e-45))):.1f}")
# rhs_s0 = float(残差) の照合 (行 0..4 = res_ro, res_roUx, res_roUy, res_roUz?, res_roe)。拘束の行は 0
for k_, r_ in runs.items():
    with h5py.File(r_ + "/res_1.h5", "r") as h:
        R = {f: h["VALUE"][f][...][idx] for f in ("res_ro", "res_roUx", "res_roUy", "res_roe") if f in h["VALUE"]}
    out = []
    for row, f in ((0, "res_ro"), (1, "res_roUx"), (2, "res_roUy"), (4, "res_roe")):
        exp = R[f].astype(np.float32).astype(np.float64)
        if row in (1, 2, 3):
            exp = np.where(fl[:, 0] == 1, 0.0, exp)
        if row == 2:
            exp = np.where(fl[:, 2] == 1, 0.0, exp)
        if row == 4:
            exp = np.where(fl[:, 1] == 1, 0.0, exp)
        out.append((row, int(np.sum(rhs[k_][:, row] != exp))))
    print(k_, "rhs_s0 と float(残差) の不一致の数 (行, 数):", out, " 行 3 の非零", int(np.sum(rhs[k_][:, 3] != 0)))
