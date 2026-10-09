"""ライン解の後退誤差 η (plan time_integration-line-implicit-speed §6.2 (2)、2026-10-09)。

forge のライン行列の書き出し (FORGE_LINE_DUMP_DIR) から、ラインごとに block 三重対角の密行列
  A: 対角 D_k、下 −Kprev_k、上 −Knext_k
を組み、各 sweep の b = rhs_s、緩和前の解 x = dqnew_s / implicitRelax について、尺度 S = diag(ρ_ref, ρ_ref a_ref ×3, ρ_ref a_ref²) で無次元化した
  η = ‖b̂ − Âx̂‖∞ / (‖Â‖∞ ‖x̂‖∞ + ‖b̂‖∞)、Â = S⁻¹AS、x̂ = S⁻¹x、b̂ = S⁻¹b
を出す。合格の基準 (§6.2) は逆行列の経路で全ライン・全 sweep η ≤ 1e-11。
usage: python3 line_eta.py <dump_dir> [<dump_dir> ...] [--ro-ref 0.8739869154 --a-ref 359.7712185] [--limit 1e-11]
       → 標準出力と <dump_dir>/eta.json。--limit を超えるものがあれば終了コード 1。
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jprobe as J  # noqa: E402


def eta_dir(d: Path, ro_ref: float, a_ref: float, nsweep: int = 5):
    A = J.load_dump(d)
    relax = A["_relax"]
    nodes, lines = A["_nodes"], A["_lines"]
    n_all = len(nodes)
    D = A["D"].reshape(n_all, 5, 5); Kp = A["Kprev"].reshape(n_all, 5, 5); Kn = A["Knext"].reshape(n_all, 5, 5)
    sc1 = np.array([ro_ref, ro_ref * a_ref, ro_ref * a_ref, ro_ref * a_ref, ro_ref * a_ref ** 2])
    out = {"dir": str(d), "relax": relax, "ro_ref": ro_ref, "a_ref": a_ref, "lines": []}
    worst = 0.0
    for L in np.unique(lines):
        idx = np.flatnonzero(lines == L); n = len(idx)
        M = np.zeros((5 * n, 5 * n))
        for k in range(n):
            M[5 * k:5 * k + 5, 5 * k:5 * k + 5] = D[idx[k]]
            if k > 0:
                M[5 * k:5 * k + 5, 5 * (k - 1):5 * k] = -Kp[idx[k]]
            if k < n - 1:
                M[5 * k:5 * k + 5, 5 * (k + 1):5 * k + 10] = -Kn[idx[k]]
        sc = np.tile(sc1, n)
        Mh = (M * sc[None, :]) / sc[:, None]
        normA = float(np.max(np.sum(np.abs(Mh), axis=1)))
        rec = {"line": int(L), "n": n, "eta": []}
        for s in range(nsweep):
            if f"rhs_s{s}" not in A:
                break
            b = A[f"rhs_s{s}"][idx].reshape(-1)
            x = A[f"dqnew_s{s}"][idx].reshape(-1) / relax
            bh, xh = b / sc, x / sc
            r = bh - Mh @ xh
            eta = float(np.max(np.abs(r)) / (normA * np.max(np.abs(xh)) + np.max(np.abs(bh))))
            rec["eta"].append(eta)
            worst = max(worst, eta)
        out["lines"].append(rec)
        print(f"{d.name} ライン {L}: η = " + " ".join(f"{v:.2e}" for v in rec["eta"]))
    out["worst"] = worst
    (d / "eta.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    return worst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--ro-ref", type=float, default=0.8739869154)
    ap.add_argument("--a-ref", type=float, default=359.7712185)
    ap.add_argument("--limit", type=float, default=None)
    a = ap.parse_args()
    w = max(eta_dir(Path(x), a.ro_ref, a.a_ref) for x in a.dirs)
    print(f"[line_eta] 最大の η = {w:.3e}" + (f" (基準 {a.limit:.0e}: {'合格' if w <= a.limit else '不合格'})" if a.limit else ""))
    sys.exit(1 if (a.limit is not None and w > a.limit) else 0)
