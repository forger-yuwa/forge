"""plan tooling-nozzle-core-grid §4.62: 双子の状態の差 (1 出力)。D_X(n) = √[Σ V (X⁺ − X)² / Σ V] (X = T・ρ・ρE、体積重み付き、拘束の節点を除く) を
島 (初期状態で固定、x/r_t ∈ [−6, −3.25]・壁から 60〜100 節点目)・縮流部 (x ∈ [−8, −2])・全域で出し、max|ΔT|・max|ΔT/T| とその場所 (x/r_t・壁からの節点番号・
局所の Δx/r_t)、両 run の dt_local の比の中央値 (島) も記録する。体積は base の出力の `volume` (無ければ止める、§4.65)。
上流 (x < −8)・粗化の遷移区間 (−8〜−7.5・−2.5〜−2) の D_T と、列ごとの max|ΔT| (壁から 50〜110 節点目) も記録する (§4.65)。
usage: python3 cg_twin_diff.py <base> <pert> <step> <tag> → _band_ab/core_grid/twin_<tag>.jsonl に 1 行足す"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; b, p, step, tag = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
info = json.loads((HERE / b / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def load(run, n):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in ("T", "ro", "roe", "dt_local", "volume") if "VALUE/" + k in h}
        xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    return v, xy
B, xy = load(b, step); P, _ = load(p, step)
if "volume" not in B: raise SystemExit(f"{b} の res_{step} に volume が無い")
_, xyp = load(p, step)
if not np.array_equal(xy, xyp): raise SystemExit("双子の座標が違う")
V = B["volume"]; N = V.size; nj = N // ni
idx = np.arange(N); i = idx // nj; jj = idx % nj; jw = nj - 1 - jj; x = xy[:, 0] / S
xc = x.reshape(ni, nj)[:, nj - 1]; dxc = np.gradient(xc)   # 壁の列の x 間隔 (r_t)
inner = (jw >= 1) & (jw <= nj - 2)
reg = {"island": (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100), "contraction": (x >= -8.0) & (x <= -2.0) & inner, "all": inner,
       "upstream": (x < -8.0) & inner, "transition": (((x >= -8.0) & (x <= -7.5)) | ((x >= -2.5) & (x <= -2.0))) & inner}
rec = {"step": step, "ni": ni, "nonfinite": int(sum(np.count_nonzero(~np.isfinite(X[k])) for X in (B, P) for k in ("T", "ro", "roe")))}
for rn, m in reg.items():
    for k in ("T", "ro", "roe"):
        d = P[k][m] - B[k][m]
        rec[f"{rn}_D_{k}"] = float(np.sqrt(np.sum(V[m] * d * d) / np.sum(V[m])))
dT = np.abs(P["T"] - B["T"]); dT[~inner] = 0.0; k = int(np.argmax(dT))
rec.update({"max_dT": float(dT[k]), "max_dT_rel": float(np.max(dT / B["T"])), "max_x": float(x[k]), "max_jw": int(jw[k]), "max_dx": float(dxc[i[k]])})
band = (jw >= 50) & (jw <= 110)
colmax = np.zeros(ni); np.maximum.at(colmax, i[band], np.abs(P["T"][band] - B["T"][band]))
rec["col_x"] = [round(float(v), 4) for v in xc]
rec["col_max_dT"] = [float(f"{v:.3e}") for v in colmax]
if "dt_local" in B:   # 列ごとの擬似 dt (帯 壁から 70〜80 節点目の中央値、基準側) と列の間隔 (§4.70、波束の上の Δτ ∝ Δx の確認)
    bd = (jw >= 70) & (jw <= 80); rows = np.where(bd)[0]
    colv = [[] for _ in range(ni)]
    for q in rows: colv[i[q]].append(B["dt_local"][q])
    rec["col_dt_band"] = [float(f"{np.median(v):.4e}") if v else None for v in colv]
    rec["col_dx"] = [round(float(v), 5) for v in dxc]
if "dt_local" in B and "dt_local" in P:
    m = reg["island"]; rec["island_dt_ratio_median"] = float(np.median(P["dt_local"][m] / B["dt_local"][m]))
with open(HERE / "_band_ab" / "core_grid" / f"twin_{tag}.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
