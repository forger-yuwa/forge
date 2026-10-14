"""plan tooling-nozzle-core-grid §4.55: M-c の資格 (c) の事後の再評価 (codex 2026-10-12 linefaces、結果を見た後の基準)。
旧登録の判定「判定不能 (資格 c)」は維持し、別出力に残す。ε_D = [(D−D_s1) + (D−D_s2) − (D−D_s3)]δ を成分別に、本判定と同じ重み・領域・拘束の除外で
‖ε_D‖ ≤ 0.01‖E‖ かつ |⟨ε_D,R₀⟩| ≤ 0.01|⟨E,R₀⟩| (射影の分母がほぼ 0 の成分は射影比を判定に使わない)。旧 (c) の式 max|dsum − d3| / max|d3| も記録する。
usage: python3 cg_conv_reeval.py <base> → _band_ab/core_grid/conv_reeval.json"""
import json, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"; B = sys.argv[1]
meta = {}
for ln in (HERE / B / "linedump" / "meta.txt").read_text().splitlines():
    if ln.startswith("#") or not ln.strip():
        continue
    n, r, c = ln.split(); meta[n] = (int(r), int(c))
g = lambda n: np.fromfile(HERE / B / "linedump" / f"{n}.f64").reshape(*meta[n])
nodes = g("node_line")[:, 0].astype(int); N = len(nodes); V = g("dt_vol")[:, 1]; con = np.any(g("flags_wall_iso_axis") != 0, axis=1)
D5 = lambda n: g(n).reshape(N, 5, 5)
D, D1, D2, D3 = D5("D"), D5("D_s1"), D5("D_s2"), D5("D_s3")
dsum = (D - D1) + (D - D2); d3 = D - D3
eps = np.einsum("nij,nj->ni", dsum - d3, g("dqold_final"))[:, [0, 1, 2, 4]]
S = json.loads((HERE / B / "prepare_info.json").read_text())["scale_m"]
with h5py.File(HERE / B / "res_0.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
with h5py.File(HERE / B / "res_1.h5", "r") as h:
    R0 = np.stack([np.asarray(h["VALUE/res_" + q][:], float)[nodes] for q in ("ro", "roUx", "roUy", "roe")], axis=1)
x = xy[:, 0] / S; jw = 159 - (nodes % 160)
isl = (x >= -6) & (x <= -3.25) & (jw >= 60) & (jw <= 100) & ~con
cj = json.loads((OUT / "conv_judge.json").read_text())
def vn(a, j): w = V[isl]; b = a[isl, j] / V[isl]; return float(np.sqrt(np.sum(w * b * b) / np.sum(w)))
def vd(a, b, j): w = V[isl]; return float(np.sum(w * (a[isl, j] / V[isl]) * (b[isl, j] / V[isl])) / np.sum(w))
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.55 (事後の再評価、結果を見た後の基準)", "base": B,
       "old_c_formula": "max|dsum − d3| / max|d3| (全要素の最大どうしの比、要素ごとの相対差の最大ではない)",
       "old_c_value": float(np.max(np.abs(dsum - d3)) / np.max(np.abs(d3))), "components": {}}
for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 3)):
    r0 = vn(R0, j); En = cj["region"]["island"][c]["E"]["norm"] * r0; Ep = cj["region"]["island"][c]["E"]["proj"] * r0 ** 2
    pr = abs(vd(eps, R0, j)) / abs(Ep) if abs(Ep) > 1e-3 * En * r0 else None
    out["components"][c] = {"eps_over_E_norm": vn(eps, j) / En, "eps_over_E_proj": pr, "ok": vn(eps, j) / En <= 0.01 and (pr is None or pr <= 0.01)}
out["new_c"] = all(out["components"][c]["ok"] for c in ("roUx", "roe"))
others = cj["qualification"]; out["other_quals"] = {k: others[k] for k in "abdef"}
I = cj["region"]["island"]; rs = [I[c][k] for c in ("roe", "roUx") for k in ("ratio_norm", "ratio_proj")]
out["reevaluated_result"] = ("B (事後の再評価): ライン外の面の対流は ρE・ρuₓ の残るずれの主要部分ではない" if out["new_c"] and all(others[k] for k in "abdef") and any(r >= 0.75 for r in rs)
                             else "再評価でも判定できない")
(OUT / "conv_reeval.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print(json.dumps(out, ensure_ascii=False, indent=1))
