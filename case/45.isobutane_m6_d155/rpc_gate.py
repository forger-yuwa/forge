"""plan architecture-float-state-double-geometry §6.55 の関門 (事前登録、codex diagnose 2026-10-12): G0c の実際の面ベクトルと hoop の面積で、pRef が作用素を実質変えないこと。
§6.44 の cp_gate.py の 1. と同じ式: 診断 run (FP64 + pRef、1 step、FORGE_DIAG_HOOP_CLOSURE) の hoop.h5 から、各実 CV で
  g = sqrt(sum64_x² + (sum64_y − A_planar)²) / (abs64_x + abs64_y + A_planar)。全 CV で g ≤ 1e-12。旧 G0 の結果は流用しない。
通れば終了コード 0、通らなければ 1。→ _band_ab/cold_pair/rpc_gate.json"""
import json, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
OUT = {"fail": [], "closure": {}}
hp = HERE / "run_0751_pc_hoop" / "hoop.h5"
if not hp.exists(): OUT["fail"].append("診断 run の hoop.h5 が無い")
else:
    with h5py.File(hp, "r") as h:
        n = len(h["cells/sum64_x"])
        sx, sy = np.asarray(h["cells/sum64_x"][:], float), np.asarray(h["cells/sum64_y"][:], float)
        ax, ay = np.asarray(h["cells/abs64_x"][:], float), np.asarray(h["cells/abs64_y"][:], float)
        A = np.asarray(h["cells/A_planar"][:], float)[:n]
        attrs = {k: str(h.attrs[k])[:80] for k in h.attrs}
    g = np.hypot(sx, sy - A) / np.maximum(ax + ay + np.abs(A), 1e-300)
    OUT["closure"] = dict(n=int(n), g_max=float(g.max()), g_p99=float(np.percentile(g, 99)), n_over=int((g > 1e-12).sum()), attrs=attrs)
    print(f"閉性 (G0c): g の最大 {g.max():.3e}、99 % {np.percentile(g, 99):.3e}、1e-12 を超える CV {int((g > 1e-12).sum())} / {n}")
    if not np.isfinite(g).all() or g.max() > 1e-12: OUT["fail"].append(f"閉性 g の最大 {g.max():.3e} > 1e-12")
(HERE / "_band_ab/cold_pair/rpc_gate.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False))
print("関門:", "通過" if not OUT["fail"] else "不成立 " + "; ".join(OUT["fail"]))
sys.exit(0 if not OUT["fail"] else 1)
