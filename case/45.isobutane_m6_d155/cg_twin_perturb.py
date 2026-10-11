"""plan tooling-nozzle-core-grid §4.62 (事前登録、codex 2026-10-12 twin): 双子の摂動。run の入力 h5 (restart_field 済み) の /VALUE/roe に、
島 (x/r_t ∈ [−6, −3.25]、壁から 60〜100 節点目、拘束の節点なし) で決まった形の微小な摂動を掛ける: roe ← roe·(1 + a·p)、
p(i, j) = [sin(2π i/8) + sin(2π i/31 + 1.1) + sin(2π i/97 + 2.3)] · [sin(π (jw − 60)/40) + 0.5 sin(3π (jw − 60)/40)] を max|p| = 1 に正規化 (i は列、jw は壁からの節点番号)。
同じ a・同じ格子なら両上限で同じ配列になる (sha256 を記録)。usage: python3 cg_twin_perturb.py <run> <a> → <run>/perturb.json"""
import hashlib, json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, a = sys.argv[1], float(sys.argv[2])
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
with h5py.File(HERE / run / "nozzle.h5", "r+") as h:
    roe = np.asarray(h["VALUE/roe"][:], dtype=np.float64)
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2] if "MESH/COORD" in h else None
    N = roe.size; nj = N // ni
    if ni * nj != N: raise SystemExit("節点数が ni の倍数でない")
    idx = np.arange(N); i = idx // nj; jw = nj - 1 - (idx % nj)
    if xy is None:
        with h5py.File(HERE / "run_0619_cg_gref_b0_ctl" / "res_60000.h5", "r") as g:
            xy = np.asarray(g["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
    x = xy[:, 0] / S
    isl = (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100)
    p = (np.sin(2 * np.pi * i / 8) + np.sin(2 * np.pi * i / 31 + 1.1) + np.sin(2 * np.pi * i / 97 + 2.3)) * \
        (np.sin(np.pi * (jw - 60) / 40) + 0.5 * np.sin(3 * np.pi * (jw - 60) / 40))
    p = np.where(isl, p, 0.0); p /= np.max(np.abs(p))
    new = roe * (1.0 + a * p)
    h["VALUE/roe"][...] = new.astype(h["VALUE/roe"].dtype)
    stored = np.asarray(h["VALUE/roe"][:], dtype=np.float64)
rec = {"a": a, "nodes": int(isl.sum()), "max_rel": float(np.max(np.abs(stored - roe) / np.abs(roe))), "dtype": str(h5py.File(HERE / run / "nozzle.h5", "r")["VALUE/roe"].dtype),
       "sha256_delta": hashlib.sha256((stored - roe).tobytes()).hexdigest()}
(HERE / run / "perturb.json").write_text(json.dumps(rec, indent=1)); print(json.dumps(rec))
