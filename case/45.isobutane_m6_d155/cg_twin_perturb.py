"""plan tooling-nozzle-core-grid §4.62・§4.65: 双子の摂動。run の入力 h5 (restart_field 済み) の /VALUE/roe に ρE·(1 + a·p) を掛ける。
p の定義 (Gref の格子で): 島 (x/r_t ∈ [−6, −3.25]、壁から 60〜100 節点目) で
  p(i, j) = [sin(2π i/8) + sin(2π i/31 + 1.1) + sin(2π i/97 + 2.3)] · [sin(π (jw − 60)/40) + 0.5 sin(3π (jw − 60)/40)] を max|p| = 1 に正規化 (i は Gref の列)。
--ref-coords を与えたときは、その (Gref の) 座標で p を作り、各 jw で x/r_t に沿って区分線形補間して対象の格子の x で評価する (§4.65、正規化は共通)。
usage: python3 cg_twin_perturb.py <run> <a> <対象の座標の h5> [--ref-coords <Gref の座標の h5>] → <run>/perturb.json"""
import argparse, hashlib, json
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("a", type=float); ap.add_argument("coords"); ap.add_argument("--ref-coords")
A = ap.parse_args()
info = json.loads((HERE / A.run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def geom(path, ni_):
    with h5py.File(HERE / path, "r") as g:
        xy = np.asarray(g["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    N = len(xy); nj = N // ni_
    if ni_ * nj != N: raise SystemExit(f"{path}: 節点数 {N} が ni {ni_} の倍数でない")
    idx = np.arange(N); return xy, idx // nj, nj - 1 - (idx % nj), nj
def pattern(xy, i, jw):
    x = xy[:, 0] / S; isl = (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100)
    p = (np.sin(2 * np.pi * i / 8) + np.sin(2 * np.pi * i / 31 + 1.1) + np.sin(2 * np.pi * i / 97 + 2.3)) * \
        (np.sin(np.pi * (jw - 60) / 40) + 0.5 * np.sin(3 * np.pi * (jw - 60) / 40))
    p = np.where(isl, p, 0.0); return p / np.max(np.abs(p))
xy, i, jw, nj = geom(A.coords, ni)
if A.ref_coords:
    with h5py.File(HERE / A.ref_coords, "r") as g: Nr = len(np.asarray(g["MESH/COORD"][:]).reshape(-1, 3))
    nir = Nr // nj
    xr, ir, jwr, njr = geom(A.ref_coords, nir)
    if njr != nj: raise SystemExit("参照と対象で nj が違う")
    pr = pattern(xr, ir, jwr); p = np.zeros(len(xy))
    for j in range(nj):
        mr = jwr == j; mt = jw == j
        xs = xr[mr, 0] / S; o = np.argsort(xs)
        p[mt] = np.interp(xy[mt, 0] / S, xs[o], pr[mr][o], left=0.0, right=0.0)
    mode = f"Gref の p を x/r_t に沿って補間 (参照 {A.ref_coords}、ni {nir})"
else:
    p = pattern(xy, i, jw); mode = "列の番号で直接"
with h5py.File(HERE / A.run / "nozzle.h5", "r+") as h:
    roe = np.asarray(h["VALUE/roe"][:], dtype=np.float64)
    if roe.size != len(xy): raise SystemExit(f"入力の節点数 {roe.size} と座標 {len(xy)} が違う")
    h["VALUE/roe"][...] = (roe * (1.0 + A.a * p)).astype(h["VALUE/roe"].dtype)
    stored = np.asarray(h["VALUE/roe"][:], dtype=np.float64); dt = str(h["VALUE/roe"].dtype)
rec = {"a": A.a, "mode": mode, "nodes": int(np.count_nonzero(p)), "max_rel": float(np.max(np.abs(stored - roe) / np.abs(roe))), "max_abs_p": float(np.max(np.abs(p))),
       "dtype": dt, "sha256_delta": hashlib.sha256((stored - roe).tobytes()).hexdigest()}
(HERE / A.run / "perturb.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1)); print(json.dumps(rec, ensure_ascii=False))
