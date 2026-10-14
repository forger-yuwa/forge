"""plan tooling-nozzle-core-grid §4.80 の記述量: 1 出力の停滞帯の格子間隔の上下 (§4.78 と同じ定義、cg_lowmach_saw.py の field_stats と同じ計算)。
線: 島 = 壁の列の x/r_t ∈ [−6, −3.25]・jw 60〜100。壁法線の線 = 島の各列の jw 58〜102、流れ方向の線 = jw 60〜100 の各行の島の列 ± 2 列。
指標 (線ごとの中央値): η = (RMS(δ⁴φ)/16)/RMS(δ¹φ)、f = δ²φ の符号が隣どうしで反転する割合、P の a/q (q = ½ρ|u|²)。
usage: python3 cg_lmp1_saw.py <run> <step> → _band_ab/core_grid/lmp1_saw_<run>.jsonl に 1 行足す"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, step = sys.argv[1], int(sys.argv[2])
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
with h5py.File(HERE / run / f"res_{step}.h5", "r") as h:
    v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("P", "ro", "T", "Ux", "Uy")}
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
nj = v["P"].size // ni
g = {k: a.reshape(ni, nj)[:, ::-1] for k, a in v.items()}   # [i, jw] (jw = 0 が壁)
g["umag"] = np.hypot(g["Ux"], g["Uy"]); g["q"] = 0.5 * g["ro"] * g["umag"] ** 2
X = (xy[:, 0] / S).reshape(ni, nj)[:, ::-1]; Y = (xy[:, 1] / S).reshape(ni, nj)[:, ::-1]

def line_stats(phi):
    d1 = np.diff(phi); d2 = phi[2:] - 2 * phi[1:-1] + phi[:-2]
    d4 = phi[4:] - 4 * phi[3:-1] + 6 * phi[2:-2] - 4 * phi[1:-3] + phi[:-4]
    a = np.sqrt(np.mean(d4 ** 2)) / 16.0; b = np.sqrt(np.mean(d1 ** 2)); c2 = d2[1:-1]
    f = float(np.mean(np.sign(c2[1:]) * np.sign(c2[:-1]) < 0)) if c2.size > 1 else float("nan")
    return a, b, (a / b if b > 0 else float("inf")), f

cols = np.where((X[:, 0] >= -6.0) & (X[:, 0] <= -3.25))[0]
L = {"wallnormal": [(i, slice(58, 103)) for i in cols], "streamwise": [(slice(int(cols[0]) - 2, int(cols[-1]) + 3), jw) for jw in range(60, 101)]}
rec = {"step": step, "nonfinite": int(sum(np.count_nonzero(~np.isfinite(g[k])) for k in ("P", "ro", "T", "Ux", "Uy")))}
for d, ls in L.items():
    acc = {}
    for sel in ls:
        xs, ys = X[sel], Y[sel]; s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))])
        acc.setdefault("eta_s", []).append(line_stats(s)[2])
        for k in ("P", "ro", "T", "umag"):
            a, b, eta, f = line_stats(g[k][sel]); acc.setdefault(f"eta_{k}", []).append(eta); acc.setdefault(f"f_{k}", []).append(f)
            if k == "P": acc.setdefault("aP_over_q", []).append(a / float(np.median(g["q"][sel])))
    rec[d] = {k: float(np.median(x)) for k, x in acc.items()}
with open(HERE / "_band_ab" / "core_grid" / f"lmp1_saw_{run}.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
