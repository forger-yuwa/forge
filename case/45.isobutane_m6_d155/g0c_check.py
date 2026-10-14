"""§6.52 の前提 (G0c の入力の照合、codex diagnose 2026-10-11 Major 3): 実際の G0 の入力 _tr/q32_init.h5 と G0c の入力 _g0c/q32_init.h5 を照合する。
- 区間の外 (G0 の列 0〜99 と 689〜4718 ↔ G0c の列 0〜99 と 437〜4466) の節点座標のビット一致
- 区間の外の 9 保存量のビット一致 (interp_field の最近傍の移植が正しいか)
- 区間の中の最近傍の対応距離 (r_t 単位) と、壁法線の列 j の対応 (同じ j に写っているか)
→ _band_ab/cold_pair/g0c_check.json"""
import json, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
def load(p):
    with h5py.File(p, "r") as h:
        xy = np.asarray(h["MESH/COORD"][:]).reshape(-1, 3)
        V = {k: np.asarray(h["VALUE/" + k][:]) for k in Q9}
    return xy, V
xy0, V0 = load(HERE / "_tr/q32_init.h5"); xy1, V1 = load(HERE / "_g0c/q32_init.h5")
NJ = 121; ni0, ni1 = len(xy0) // NJ, len(xy1) // NJ
_, _, S = XC.mesh_info(len(xy0))
X0 = xy0[:, :2].reshape(ni0, NJ, 2); X1 = xy1[:, :2].reshape(ni1, NJ, 2)
ia = 100; nseg0, nseg1 = 589, 337
pre0, pre1 = slice(0, ia), slice(0, ia)
post0, post1 = slice(ia + nseg0, ni0), slice(ia + nseg1, ni1)
out = {"ni": [ni0, ni1], "coord_dtype": [str(xy0.dtype), str(xy1.dtype)]}
out["pre_coords_bit_equal"] = bool(np.array_equal(X0[pre0], X1[pre1]))
out["post_coords_bit_equal"] = bool(np.array_equal(X0[post0], X1[post1]))
out["post_cols"] = [ni0 - (ia + nseg0), ni1 - (ia + nseg1)]
out["max_abs_coord_diff_outside_rt"] = float(max(np.max(np.abs(X0[pre0] - X1[pre1])), np.max(np.abs(X0[post0] - X1[post1]))) / S)
for nm, s0, s1 in (("pre", pre0, pre1), ("post", post0, post1)):
    bad = {}
    for k in Q9:
        a = V0[k].reshape(ni0, NJ)[s0]; b = V1[k].reshape(ni1, NJ)[s1]
        n = int(np.count_nonzero(a != b))
        if n: bad[k] = n
    out[f"{nm}_values_not_bit_equal"] = bad
# 区間の中: 最近傍の対応
from scipy.spatial import cKDTree  # noqa: E402
tree = cKDTree(xy0[:, :2])
seg = X1[ia:ia + nseg1].reshape(-1, 2)
d, idx = tree.query(seg)
j_src = idx % NJ; j_dst = np.tile(np.arange(NJ), nseg1)
out["seg_nearest_dist_rt"] = dict(max=float(d.max() / S), p99=float(np.percentile(d, 99) / S), median=float(np.median(d) / S))
out["seg_same_j_fraction"] = float(np.mean(j_src == j_dst))
out["seg_j_mismatch_max"] = int(np.max(np.abs(j_src - j_dst)))
ok = out["pre_coords_bit_equal"] and out["post_coords_bit_equal"] and not out["pre_values_not_bit_equal"] and not out["post_values_not_bit_equal"]
out["verdict"] = "照合 OK (区間の外は座標と 9 保存量がビット一致)" if ok else "照合 NG"
print(json.dumps(out, ensure_ascii=False, indent=1))
(HERE / "_band_ab/cold_pair/g0c_check.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("VERDICT:", out["verdict"])
sys.exit(0 if ok else 1)
