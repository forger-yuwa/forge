"""局所収縮率 g = |dq_final| / |dq_sub<m>| の解析 (line-implicit v2 診断, 2026-09-03)。

データ: run_diag_cflmap_cp8_g*/res_100.h5 (`FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=<m>`)。
  - dq_block_old_{0..4}: 最終 subiter (19) の補正 (swap 後)
  - dq_{ro,roUx,roUy,roUz,roe}_new: subiter <m> の補正 (退避スロット)

**マスクの注意 (Codex 指摘 2026-09-03)**: 分母 (初期 dq) に活動点マスク (例: ini>p50) を
入れると「初期が小さく後から成長する点」を定義から除外し、増幅点の不存在証明には使えない。
本スクリプトは ① マスクなし ② 最終側マスク (|dq_final| > 1% of max) の両方を明示的に出す。
usage: python3 analyze_g.py [run_dir=run_diag_cflmap_cp8_g]
"""
import sys
import h5py
import numpy as np

CASE = "/home/sano/work/forge/case/39.periodic_hills"
run = sys.argv[1] if len(sys.argv) > 1 else "run_diag_cflmap_cp8_g"
f = h5py.File(f"{CASE}/{run}/res_100.h5")
m = h5py.File(f"{CASE}/run_0023_prod_ddes_ny160_cont/hill_xc_160x160x60.h5")
n = len(f["VALUE/ro"])
co = np.array(m["MESH/COORD"][:]).reshape(-1, 3)[:n]
x, y = co[:, 0], co[:, 1]
h = 0.028
wd = f["VALUE/wall_dist"][:]

eqs = ["ro", "roUx", "roUy", "roUz", "roe"]
snap_names = ["dq_ro_new", "dq_roUx_new", "dq_roUy_new", "dq_roUz_new", "dq_roe_new"]

print(f"=== {run}: g = |dq_final(subiter19)| / |dq_snap| ===")
for i, eq in enumerate(eqs):
    fin = np.abs(f[f"VALUE/dq_block_old_{i}"][:])
    ini = np.abs(f[f"VALUE/{snap_names[i]}"][:])
    g = fin / np.maximum(ini, 1e-30)
    grow = g > 1
    print(f"\n[{eq}] マスクなし: med={np.median(g):.3g} p99={np.percentile(g,99):.3g} "
          f"g>1: {grow.sum()} ノード")
    mask_fin = fin > 0.01 * fin.max()
    gm = grow & mask_fin
    print(f"  最終側マスク (|dq_fin|>1%max, n={mask_fin.sum()}): g>1 = {gm.sum()} ノード")
    if gm.sum():
        print(f"    wall_dist/h: med={np.median(wd[gm])/h:.3f} p90={np.percentile(wd[gm],90)/h:.3f}")
        print(f"    x/h ヒスト [0,1,3,5,7,9]:", np.histogram(x[gm] / h, bins=[0, 1, 3, 5, 7, 9])[0])
