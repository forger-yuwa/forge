"""plan tooling-nozzle-core-grid §4.98 P-hi の格子の資格 (事前登録): H の nozzle.h5 を Gref (run_0608 の nozzle.h5) と比べる。
(1) x ≥ 0 の節点座標が Gref と相対 1e−12 以内 (r_t で割った座標の差 ≤ 1e−12)。(2) x の station (壁の列の x) が全列で Gref と一致。
(3) 第一層の実寸 (壁の節点と 1 つ内側の節点の距離) が全列で Gref と相対 1 % 以内。(4) check_mesh_quality.py --ar-max 5000 が PASS。
(5) 壁の折れ角 (隣り合う壁の線分の向きの差) の最大 ≤ 2°。(6) 初期値の保存量に非有限なし・ρ > 0 (記録として入口・縮流部の主流の M も出す)。
usage (AWS の case/45 で): python3 cg_hi_mesh_check.py <H の準備のディレクトリ> → <dir>/hi_mesh_check.json"""
import json, subprocess, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent
D = Path(sys.argv[1]).resolve(); GREF = HERE / "run_0608_cg_gref_cont_k1"
S = json.loads((GREF / "prepare_info.json").read_text())["scale_m"]; NI = int(json.loads((GREF / "prepare_info.json").read_text())["mesh"]["ni"])
def coords(p):
    with h5py.File(p, "r") as h: return np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2] / S
c0, c1 = coords(GREF / "nozzle.h5"), coords(D / "nozzle.h5")
if c0.shape != c1.shape: raise SystemExit(f"節点数が違う ({c0.shape} 対 {c1.shape})")
NJ = c0.shape[0] // NI
P0, P1 = c0.reshape(NI, NJ, 2)[:, ::-1], c1.reshape(NI, NJ, 2)[:, ::-1]     # [i, jw] (jw = 0 が壁)
out = {"plan": "plans/active/tooling-nozzle-core-grid.md §4.98", "dir": str(D)}
down = P0[:, 0, 0] >= 0.0
out["downstream_max_diff_rt"] = float(np.abs(P0[down] - P1[down]).max())
out["x_station_max_diff_rt"] = float(np.abs(P0[:, 0, 0] - P1[:, 0, 0]).max())
t0 = np.linalg.norm(P0[:, 0] - P0[:, 1], axis=1); t1 = np.linalg.norm(P1[:, 0] - P1[:, 1], axis=1)
out["first_layer_ratio_min_max"] = [float((t1 / t0).min()), float((t1 / t0).max())]
w = P1[:, 0]; v = np.diff(w, axis=0); ang = np.degrees(np.abs(np.diff(np.arctan2(v[:, 1], v[:, 0]))))
out["wall_kink_max_deg"] = float(ang.max())
q = subprocess.run([sys.executable, str(HERE.parents[1] / "solver_density_cuda" / "tools" / "check_mesh_quality.py"), str(D / "nozzle.h5"), "--ar-max", "5000"],
                   capture_output=True, text=True)
out["mesh_quality_tail"] = q.stdout.strip().splitlines()[-6:]; out["mesh_quality_pass"] = ("VERDICT: PASS" in q.stdout)
with h5py.File(D / "nozzle.h5", "r") as h:
    V = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("ro", "roUx", "roUy", "roe") if "VALUE/" + k in h}
out["ic_nonfinite"] = int(sum(np.count_nonzero(~np.isfinite(a)) for a in V.values())); out["ic_min_ro"] = float(V["ro"].min())
qual = {"downstream_same": out["downstream_max_diff_rt"] <= 1e-12, "x_station_same": out["x_station_max_diff_rt"] <= 1e-12,
        "first_layer_within_1pct": abs(out["first_layer_ratio_min_max"][0] - 1) <= 0.01 and abs(out["first_layer_ratio_min_max"][1] - 1) <= 0.01,
        "mesh_quality_pass": out["mesh_quality_pass"], "wall_kink_le_2deg": out["wall_kink_max_deg"] <= 2.0,
        "ic_finite_positive": out["ic_nonfinite"] == 0 and out["ic_min_ro"] > 0}
out["qualification"] = qual; out["pass"] = all(qual.values())
(D / "hi_mesh_check.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print(json.dumps(out, ensure_ascii=False, indent=1))
