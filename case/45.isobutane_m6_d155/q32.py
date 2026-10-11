"""plan architecture-float-state-double-geometry §6.8 (V3): run の場のファイル (nozzle.h5) の VALUE の浮動小数点の量を float に丸め、
double のまま書き戻す (float のビルドも FP64 のビルドも同じ値 Q32 を読む)。丸めた後の値の sha256 を出す (腕どうしの一致の確認用)。
usage: python3 q32.py <run_dir>"""
import hashlib, sys, h5py, numpy as np
from pathlib import Path
run = Path(sys.argv[1]); f = run / "nozzle.h5"
hs = hashlib.sha256(); n = 0
with h5py.File(f, "r+") as h:
    V = h["VALUE"]
    for k in sorted(V.keys()):
        a = np.asarray(V[k][:])
        if a.dtype.kind != "f": continue
        if not np.all(np.isfinite(a)): raise SystemExit(f"{f}: VALUE/{k} に非有限")
        b = a.astype(np.float32).astype(a.dtype)
        V[k][...] = b; n += 1
        hs.update(k.encode()); hs.update(np.ascontiguousarray(b).tobytes())
print(f"[q32] {run.name}: VALUE の {n} 量を float に丸めた (型は {a.dtype} のまま)、sha256 {hs.hexdigest()[:16]}")
