"""case/45 H の格子の配点 (YAML の x_density_table・wall_first_frac_table、ni 4719) を Hall 初期線の壁で置いたときの壁際の AR と、
第 1 層の表 / 見積もりの要求 の区間ごとの比 (使い捨て)。"""
import sys, json, pickle
from pathlib import Path
import numpy as np, yaml
sys.path.insert(0, str(Path.home() / "forge-mknot/design"))
from forge_design.meshing import mesh2d as M2
D = Path.home() / "m55_mesh_est" / "c45"
wall = pickle.load(open(D / "wall.pkl", "rb"))
m = yaml.safe_load(open(D / "problem_h_hall.yaml"))["mesh"]
ft = np.asarray(m["wall_first_frac_table"], float); dt = np.asarray(m["x_density_table"], float)
x0, x1 = float(wall.x_in), float(wall.x_e)
ft[-1, 0] = max(ft[-1, 0], x1); dt[-1, 0] = max(dt[-1, 0], x1); ft[0, 0] = min(ft[0, 0], x0); dt[0, 0] = min(dt[0, 0], x0)
xs = M2._x_stations(x0, x1, int(m["ni"]), 1.0, 1.0, density_table=dt)
fr = np.exp(np.interp(xs, ft[:, 0], np.log(ft[:, 1])))
r = wall.r(xs); y1 = fr * r; ch = np.hypot(np.diff(xs), np.diff(r))
ar = np.maximum(ch, np.maximum(y1[:-1], y1[1:])) / np.minimum(ch, np.minimum(y1[:-1], y1[1:]))
k = int(np.argmax(ar)); dx = np.diff(xs); q = dx[1:] / dx[:-1]
print(f"H 配点 ni {m['ni']}: 壁際 AR 最大 {ar[k]:.1f} @ x {0.5 * (xs[k] + xs[k + 1]):.3f}; AR>4000 の列 {(ar > 4000).mean() * 100:.1f} %; "
      f"dx {dx.min():.5f}–{dx.max():.4f} r_t; dx 比最大 {np.max(np.maximum(q, 1 / q)):.4f}")
sec = np.load(D / "sections.npz"); x = sec["x_rt"]; f = sec["y1_over_rw"]
rat = np.exp(np.interp(x, ft[:, 0], np.log(ft[:, 1]))) / np.where(f > 0, f, np.nan)
for nm, lo, hi in (("入口の配管", -12.5, -9.0), ("縮流部 x<-2", -9.0, -2.0), ("スロート付近 |x|≤2", -2.0, 2.0), ("拡大部 2<x≤x_e/2", 2.0, x1 / 2),
                   ("拡大部 x>x_e/2", x1 / 2, x1)):
    mm = (x >= lo) & (x <= hi) & np.isfinite(rat)
    print(f"  {nm:<18s} H の第 1 層 / 見積もりの要求: {np.nanmin(rat[mm]):.3f}–{np.nanmax(rat[mm]):.3f} (中央 {np.nanmedian(rat[mm]):.3f})")
