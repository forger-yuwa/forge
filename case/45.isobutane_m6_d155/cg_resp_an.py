"""plan tooling-nozzle-core-grid §4.46 M-r (事前登録、codex 2026-10-11 response): 十分解いた補正に対する実際の残差の応答。
停滞場から 2 外反復 (nStepInner K = 42、SST の保存量の更新は凍結)。最初の外反復の書き出し (linedump) と、res_1 (更新前の残差 R₀ と更新後の状態 Q₁)・
res_2 (更新後の状態の残差 R₁)・res_0 (状態 Q₀) を使う。
記号: L はライン行列 (D、−Kprev、−Knext)、T = diag(V/Δτ)、ω = implicitRelax、c = 1/ω − 1、b_k = rhs_s{k}、δ_k = dqnew_s{k−1}。
  確定した補正 δ = δ_K (= dqnew_s{K−1})。Nδ ≈ b₀ − b_{K−1} (= Nδ_{K−1})。その誤差 N(δ_K − δ_{K−1}) は測れないので、
  測れる ‖b_{K−1} − b_{K−2}‖ = ‖N(δ_{K−1} − δ_{K−2})‖ を資格で抑える (codex 案の「最後の swap の後に右辺だけ再評価」はコード変更が要るため、sweep を 1 回足す形に替えた)。
  H = L/ω + N、A = L + N − T、g = Hδ − R₀、ΔR_lhs = −Aδ、ΔR_act = R₁ − R₀、e = ΔR_act − ΔR_lhs、恒等式 R₁ = Tδ + cLδ − g + e。
指標 (成分 ρ・ρu・ρv・ρE、体積重み付きの R/V のノルム、拘束の節点を除く、領域は島と高さの帯):
  m = ‖e‖/‖ΔR_lhs‖、ρ_act = ‖R₁‖/‖R₀‖、ρ_lhs = ‖R₀ + ΔR_lhs‖/‖R₀‖、減衰 d_act = 1 − ρ_act²、d_lhs = 1 − ρ_lhs²、不足 = d_lhs − d_act、
  Tδ・cLδ・−g・e の ‖·‖/‖R₀‖ と R₀ への符号付き射影 ⟨X, R₀⟩/‖R₀‖²。
資格: (q1) Q₁ − Q₀ (res_1 − res_0 の保存量) と δ の一致 (成分別の相対差 ≤ ε)、(q2) R₀ (res_1) と b₀ の一致 (成分別の相対差 ≤ ε)、
  (q3) ρk・ρω が res_0 と res_1 でビット一致 (凍結の実効)、(q4) 組成の変化 max|ΔY_s| (島) ≤ 0.01·max|δρ/ρ| (島)、
  (q5) ‖b_{K−1} − b_{K−2}‖/‖b₀‖ ≤ 0.1ε (Nδ の近似誤差)、(q6) ‖g‖ ≤ 0.01‖ΔR_lhs‖、(q7) 信号 ‖ΔR_lhs‖/‖R₀‖ ≥ 100ε と |d_lhs| ≥ 100ε。
  ε は成分ごとに max(1e−7, q1・q2 の実測の相対差)。
usage: python3 cg_resp_an.py <run> <tag> → _band_ab/core_grid/resp_<tag>.json"""
import json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; run, tag = sys.argv[1], sys.argv[2]
DD = HERE / run / "linedump"
meta = {}
for ln in (DD / "meta.txt").read_text().splitlines():
    if ln.startswith("#") or not ln.strip():
        continue
    n, r, c = ln.split(); meta[n] = (int(r), int(c))
get = lambda name: np.fromfile(DD / f"{name}.f64", dtype=np.float64).reshape(*meta[name])
K = max(int(n[len("dqnew_s"):]) for n in meta if n.startswith("dqnew_s")) + 1
if K != 42:
    raise SystemExit(f"sweep が {K} 回 (42 回のはず)")
cfg = (HERE / run / "solverConfig.yaml").read_text()
key = lambda k: float(re.search(k + r":\s*([0-9.eE+-]+)", cfg).group(1))
om = key("implicitRelax"); c_ = 1.0 / om - 1.0
info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
nl = get("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64); N = len(nodes)
prv = np.arange(N) - 1; nxt = np.arange(N) + 1
prv[np.r_[True, line[1:] != line[:-1]]] = -1; nxt[np.r_[line[1:] != line[:-1], True]] = -1
Dm = get("D").reshape(N, 5, 5); Kp = get("Kprev").reshape(N, 5, 5); Kn = get("Knext").reshape(N, 5, 5)
dv = get("dt_vol"); dt, V = dv[:, 0], dv[:, 1]
con = np.any(get("flags_wall_iso_axis") != 0, axis=1)
def Lmul(d):
    out = np.einsum("nij,nj->ni", Dm, d)
    m = prv >= 0; out[m] -= np.einsum("nij,nj->ni", Kp[m], d[prv[m]])
    m = nxt >= 0; out[m] -= np.einsum("nij,nj->ni", Kn[m], d[nxt[m]])
    return out
b0, bK1, bK2 = get("rhs_s0"), get(f"rhs_s{K - 1}"), get(f"rhs_s{K - 2}")
delta = get(f"dqnew_s{K - 1}")                                    # 確定した補正 δ_K
Td = (V / dt)[:, None] * delta
Ld = Lmul(delta)
Nd = b0 - bK1                                                     # ≈ Nδ (誤差 N(δ_K − δ_{K−1}))
H_d = Ld / om + Nd
A_d = Ld + Nd - Td
CMP = {"ro": 0, "roUx": 1, "roUy": 2, "roe": 4}
def load(n, names):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], float)[nodes] for k in names}, np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
Q = ["ro", "roUx", "roUy", "roe"]; R = ["res_" + q for q in Q]
s0, xy = load(0, Q + ["roK", "roOmega", "roY0", "roY1"])
s1, _ = load(1, Q + R + ["roK", "roOmega", "roY0", "roY1"])
s2, _ = load(2, R)
x = xy[:, 0] / S; nj = 160; jw = nj - 1 - (nodes % nj)
as5 = lambda d: np.stack([d[k] if k in d else np.zeros(N) for k in ("ro", "roUx", "roUy", "roUz", "roe")], axis=1)
R0 = np.stack([s1["res_" + q] for q in Q], axis=1); R1 = np.stack([s2["res_" + q] for q in Q], axis=1)
dQ = np.stack([s1[q] - s0[q] for q in Q], axis=1)
ci = [0, 1, 2, 4]
sel = lambda a: a[:, ci]
dl, b0s, Tds, Lds, Hs, As = map(sel, (delta, b0, Td, Ld, H_d, A_d))
g = Hs - R0; dRl = -As; dRa = R1 - R0; e = dRa - dRl
reg = {"island": (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100) & ~con, "lines": ~con}
for lo, hi in ((60, 69), (70, 79), (80, 89), (90, 100)):
    reg[f"jw{lo}_{hi}"] = (x >= -6.0) & (x <= -3.25) & (jw >= lo) & (jw <= hi) & ~con
def vn(a, m, j):
    w = V[m]; b = a[m, j] / V[m]
    return float(np.sqrt(np.sum(w * b * b) / np.sum(w)))
def vdot(a, b, m, j):
    w = V[m]
    return float(np.sum(w * (a[m, j] / V[m]) * (b[m, j] / V[m])) / np.sum(w))
def rel(a, b, m, j):                                             # ‖a − b‖/‖b‖ (a・b は同じ単位、V で割らない量にも使う)
    return float(np.linalg.norm(a[m, j] - b[m, j]) / max(np.linalg.norm(b[m, j]), 1e-300))
out = {"run": run, "tag": tag, "K": K, "omega": om, "nodes": N, "constrained": int(con.sum()), "region": {}}
isl = reg["island"]
q = {}
q["q3_sst_frozen"] = bool(np.array_equal(s0["roK"], s1["roK"]) and np.array_equal(s0["roOmega"], s1["roOmega"]))
Y0 = {k: s0[k] / s0["ro"] for k in ("roY0", "roY1")}; Y1 = {k: s1[k] / s1["ro"] for k in ("roY0", "roY1")}
dY = max(float(np.max(np.abs(Y1[k][isl] - Y0[k][isl]))) for k in Y0); drr = float(np.max(np.abs(dQ[isl, 0] / s0["ro"][isl])))
q["q4_dY_island"] = dY; q["q4_drho_rel_island"] = drr; q["q4_ok"] = dY <= 0.01 * drr
for cn, j in CMP.items():
    jj = ci.index(j)
    eps_q1 = rel(dQ, dl, ~con, jj); eps_q2 = rel(R0, b0s, ~con, jj)
    eps = max(1e-7, eps_q1, eps_q2)
    q[f"eps_{cn}"] = eps; q[f"q1_{cn}"] = eps_q1; q[f"q2_{cn}"] = eps_q2
    q[f"q5_{cn}"] = float(np.linalg.norm(bK1[~con, j] - bK2[~con, j]) / np.linalg.norm(b0[~con, j]))
out["qual_raw"] = q
for rn, m in reg.items():
    o = {"n": int(m.sum())}
    for cn, j in CMP.items():
        jj = ci.index(j); r0 = vn(R0, m, jj); nl_ = vn(dRl, m, jj)
        rho_act = vn(R1, m, jj) / r0; rho_lhs = vn(R0 + dRl, m, jj) / r0
        o[cn] = {"R0": r0, "dRlhs_over_R0": nl_ / r0, "m": vn(e, m, jj) / nl_, "g_over_dRlhs": vn(g, m, jj) / nl_,
                 "rho_act": rho_act, "rho_lhs": rho_lhs, "d_act": 1 - rho_act ** 2, "d_lhs": 1 - rho_lhs ** 2,
                 "deficit": (1 - rho_lhs ** 2) - (1 - rho_act ** 2),
                 "terms_norm": {k: vn(X, m, jj) / r0 for k, X in (("Td", Tds), ("cLd", c_ * Lds), ("-g", -g), ("e", e))},
                 "terms_proj": {k: vdot(X, R0, m, jj) / r0 ** 2 for k, X in (("Td", Tds), ("cLd", c_ * Lds), ("-g", -g), ("e", e), ("R1", R1))},
                 "identity_check": vn(R1 - (Tds + c_ * Lds - g + e), m, jj) / r0}
    out["region"][rn] = o
(HERE / "_band_ab" / "core_grid" / f"resp_{tag}.json").write_text(json.dumps(out, indent=1))
print(f"{tag}: 節点 {N}、拘束 {out['constrained']}、SST 凍結 {q['q3_sst_frozen']}、組成 max|ΔY| {dY:.1e} (max|δρ/ρ| {drr:.1e})")
print("  ε・q1・q2・q5:", {cn: (f"{q['eps_' + cn]:.1e}", f"{q['q1_' + cn]:.1e}", f"{q['q2_' + cn]:.1e}", f"{q['q5_' + cn]:.1e}") for cn in CMP})
for rn in ("island", "jw60_69", "jw70_79", "jw80_89", "jw90_100"):
    o = out["region"][rn]
    print(f"  {rn:9s} " + "  ".join(f"{cn}: ΔRlhs/R0 {o[cn]['dRlhs_over_R0']:.2e} m {o[cn]['m']:.3f} ρact {o[cn]['rho_act']:.4f} ρlhs {o[cn]['rho_lhs']:.4f}" for cn in ("ro", "roe")))
