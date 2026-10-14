"""plan tooling-nozzle-core-grid §4.48 M-o の解析と判定 (事前登録、codex 2026-10-11 mismatch): 応答のずれ e のうち、右辺の再構成次数の寄与。
入力: base (診断入りバイナリで M-r 上限 50 を回し直した run、linedump に rhs_final) と、eval の 4 run (Q₀・Q₁ × convMethod 1・0)、比較用の旧 run_0632。
記号: L はライン行列 (D、−Kprev、−Knext)、T = diag(V/Δτ)、ω = implicitRelax、δ = 確定した補正 (dqold_final = dqnew_s41)、
  b(δ) = rhs_final = R₀ − Nδ (確定した δ に対する直接の評価) → Nδ = b₀ − b(δ)、Aδ = Lδ + Nδ − Tδ、真の g = Lδ/ω − b(δ)。
  ΔR₂ = R₂(Q₁) − R₂(Q₀)、ΔR₁ = R₁(Q₁) − R₁(Q₀) (下付きは再構成の次数)、e₂ = ΔR₂ + Aδ、e₁ = ΔR₁ + Aδ、C_rec = ΔR₂ − ΔR₁ = e₂ − e₁。
ノルムは体積重み付きの R/V、拘束の節点を除く、Lδ は全ラインで計算してから領域を切り出す。主判定は島の ρ・ρE、ρu・ρv は併記。
資格: (a) 診断入りバイナリが旧 run_0632 を再現 (rhs_s0 の成分別の差 ≤ 1e−8、δ の差 ≤ 1e−6)、(b) dqold_final = dqnew_s41 (ビット一致)、
  (c) 島の真の ‖g‖/‖Aδ‖ ≤ 0.01、(d) 再現: ‖R₂(Q₀) − R₀‖ と ‖R₂(Q₁) − R₁‖ が島で ‖ΔR₂‖ の 1 % 以下 (R₀・R₁ は base の res_1・res_2)、
  (e) base で ρk・ρω がビット一致 (SST 凍結) と組成の変化 ≤ 0.01·max|δρ/ρ|。
結果 A: 島の ρ・ρE の両方で ‖e₁‖ ≤ 0.5‖e₂‖ かつ |⟨e₁,R₀⟩| ≤ 0.5|⟨e₂,R₀⟩| → 再構成の変更が不整合の主要部分を除く (1 次化を対策としては採らない)。
結果 B: 島の ρ・ρE の両方で ‖C_rec‖ ≤ 0.1‖e₂‖ かつ |⟨C_rec,R₀⟩| ≤ 0.1|⟨e₂,R₀⟩| → この補正の向きでは再構成の寄与は小さく、第 1 仮説を退ける。
それ以外・資格不足は判定不能。
usage: python3 cg_mis_an.py <base> <o2q0> <o2q1> <o1q0> <o1q1> <旧 run_0632> → _band_ab/core_grid/mis_judge.json"""
import json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
base, o2q0, o2q1, o1q0, o1q1, old = sys.argv[1:7]
def dumpget(run):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    return lambda name: np.fromfile(HERE / run / "linedump" / f"{name}.f64", dtype=np.float64).reshape(*meta[name]), meta
get, meta = dumpget(base); gold, _ = dumpget(old)
cfg = (HERE / base / "solverConfig.yaml").read_text()
om = float(re.search(r"implicitRelax:\s*([0-9.eE+-]+)", cfg).group(1))
info = json.loads((HERE / base / "prepare_info.json").read_text()); S = float(info["scale_m"])
nl = get("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64); N = len(nodes)
prv = np.arange(N) - 1; nxt = np.arange(N) + 1
prv[np.r_[True, line[1:] != line[:-1]]] = -1; nxt[np.r_[line[1:] != line[:-1], True]] = -1
Dm = get("D").reshape(N, 5, 5); Kp = get("Kprev").reshape(N, 5, 5); Kn = get("Knext").reshape(N, 5, 5)
dv = get("dt_vol"); dt, V = dv[:, 0], dv[:, 1]
con = np.any(get("flags_wall_iso_axis") != 0, axis=1)
def Lmul(d):
    o = np.einsum("nij,nj->ni", Dm, d)
    m = prv >= 0; o[m] -= np.einsum("nij,nj->ni", Kp[m], d[prv[m]])
    m = nxt >= 0; o[m] -= np.einsum("nij,nj->ni", Kn[m], d[nxt[m]])
    return o
delta = get("dqold_final"); b0 = get("rhs_s0"); bf = get("rhs_final")
Nd = b0 - bf; Ld = Lmul(delta); Td = (V / dt)[:, None] * delta
Ad = Ld + Nd - Td; g = Ld / om - bf
Nd_old = b0 - get("rhs_s41")                                       # 旧近似 (M-r の Nδ ≈ b₀ − b₄₁)
Q = ["ro", "roUx", "roUy", "roe"]; ci = [0, 1, 2, 4]
def res(run, n):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return np.stack([np.asarray(h["VALUE/res_" + q][:], float)[nodes] for q in Q], axis=1)
def state(run, n, names):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], float)[nodes] for k in names}, np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
R0, R1 = res(base, 1), res(base, 2)
R2q0, R2q1, R1q0, R1q1 = res(o2q0, 1), res(o2q1, 1), res(o1q0, 1), res(o1q1, 1)
s0, xy = state(base, 0, Q + ["roK", "roOmega", "roY0", "roY1"]); s1, _ = state(base, 1, Q + ["roK", "roOmega", "roY0", "roY1"])
x = xy[:, 0] / S; nj = 160; jw = nj - 1 - (nodes % nj)
sel = lambda a: a[:, ci]
A_, g_, dlt, Nd_, Ndo = map(sel, (Ad, g, delta, Nd, Nd_old))
dR2 = R2q1 - R2q0; dR1 = R1q1 - R1q0
e2 = dR2 + A_; e1 = dR1 + A_; Crec = dR2 - dR1
reg = {"island": (x >= -6.0) & (x <= -3.25) & (jw >= 60) & (jw <= 100) & ~con, "lines": ~con}
for lo, hi in ((60, 69), (70, 79), (80, 89), (90, 100)):
    reg[f"jw{lo}_{hi}"] = (x >= -6.0) & (x <= -3.25) & (jw >= lo) & (jw <= hi) & ~con
def vn(a, m, j):
    w = V[m]; b = a[m, j] / V[m]
    return float(np.sqrt(np.sum(w * b * b) / np.sum(w)))
def vdot(a, b, m, j):
    w = V[m]
    return float(np.sum(w * (a[m, j] / V[m]) * (b[m, j] / V[m])) / np.sum(w))
CMP = {"ro": 0, "roUx": 1, "roUy": 2, "roe": 3}
isl = reg["island"]
qual = {}
# (a) 旧 run_0632 の再現
rb0, ro0 = get("rhs_s0"), gold("rhs_s0")
qual["a_rhs_s0_rel"] = {c: float(np.linalg.norm(rb0[~con, j] - ro0[~con, j]) / np.linalg.norm(ro0[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
dd_old = gold("dqnew_s41")
qual["a_delta_rel"] = {c: float(np.linalg.norm(delta[~con, j] - dd_old[~con, j]) / np.linalg.norm(dd_old[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
qual["a"] = all(v <= 1e-8 for v in qual["a_rhs_s0_rel"].values()) and all(v <= 1e-6 for v in qual["a_delta_rel"].values())
qual["b"] = bool(np.array_equal(get("dqold_final"), get("dqnew_s41")))
qual["c_g_over_Ad"] = {c: vn(g_, isl, j) / vn(A_, isl, j) for c, j in CMP.items()}
qual["c"] = all(qual["c_g_over_Ad"][c] <= 0.01 for c in ("ro", "roe"))
qual["d_repro"] = {c: max(vn(R2q0 - R0, isl, j), vn(R2q1 - R1, isl, j)) / vn(dR2, isl, j) for c, j in CMP.items()}
qual["d"] = all(qual["d_repro"][c] <= 0.01 for c in ("ro", "roe"))
Y0 = {k: s0[k] / s0["ro"] for k in ("roY0", "roY1")}; Y1 = {k: s1[k] / s1["ro"] for k in ("roY0", "roY1")}
dY = max(float(np.max(np.abs(Y1[k][isl] - Y0[k][isl]))) for k in Y0); drr = float(np.max(np.abs((s1["ro"][isl] - s0["ro"][isl]) / s0["ro"][isl])))
qual["e_sst_frozen"] = bool(np.array_equal(s0["roK"], s1["roK"]) and np.array_equal(s0["roOmega"], s1["roOmega"]))
qual["e_dY"] = dY; qual["e_drho"] = drr
qual["e"] = qual["e_sst_frozen"] and dY <= 0.01 * drr
qual["Nd_old_vs_direct"] = {c: vn(Ndo - Nd_, isl, j) / vn(A_, isl, j) for c, j in CMP.items()}   # 旧近似の誤差 (‖Aδ‖ 比、記録)
ok = qual["a"] and qual["b"] and qual["c"] and qual["d"] and qual["e"]
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.48", "runs": dict(zip(("base", "o2q0", "o2q1", "o1q0", "o1q1", "old"), sys.argv[1:7])),
       "qualification": qual, "region": {}}
for rn, m in reg.items():
    o = {}
    r0 = {c: vn(R0, m, j) for c, j in CMP.items()}
    for c, j in CMP.items():
        o[c] = {"Ad_over_R0": vn(A_, m, j) / r0[c],
                "e2": vn(e2, m, j) / r0[c], "e1": vn(e1, m, j) / r0[c], "Crec": vn(Crec, m, j) / r0[c],
                "e2_proj": vdot(e2, R0, m, j) / r0[c] ** 2, "e1_proj": vdot(e1, R0, m, j) / r0[c] ** 2, "Crec_proj": vdot(Crec, R0, m, j) / r0[c] ** 2,
                "dR2_proj": vdot(dR2, R0, m, j) / r0[c] ** 2, "dR1_proj": vdot(dR1, R0, m, j) / r0[c] ** 2, "Ad_proj": vdot(-A_, R0, m, j) / r0[c] ** 2}
    out["region"][rn] = o
I = out["region"]["island"]
if not ok:
    res_ = "判定不能 (資格不足: " + ", ".join(k for k in ("a", "b", "c", "d", "e") if not qual[k]) + ")"
elif all(I[c]["e1"] <= 0.5 * I[c]["e2"] and abs(I[c]["e1_proj"]) <= 0.5 * abs(I[c]["e2_proj"]) for c in ("ro", "roe")):
    res_ = "A: 再構成の変更が応答のずれの主要部分を除く"
elif all(I[c]["Crec"] <= 0.1 * I[c]["e2"] and abs(I[c]["Crec_proj"]) <= 0.1 * abs(I[c]["e2_proj"]) for c in ("ro", "roe")):
    res_ = "B: この補正の向きでは再構成の寄与は小さい (第 1 仮説を退ける)"
else:
    res_ = "判定不能 (中間・成分の不一致)"
out["result"] = res_
(OUT / "mis_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("資格:", {k: qual[k] for k in ("a", "b", "c", "d", "e")}, " g/Aδ", {k: f"{v:.1e}" for k, v in qual["c_g_over_Ad"].items()}, " 再現", {k: f"{v:.1e}" for k, v in qual["d_repro"].items()})
print("  旧近似の Nδ の誤差 (‖Aδ‖ 比)", {k: f"{v:.1e}" for k, v in qual["Nd_old_vs_direct"].items()}, " a:", qual["a_rhs_s0_rel"], qual["a_delta_rel"])
for rn in ("island", "jw60_69", "jw70_79", "jw80_89", "jw90_100"):
    o = out["region"][rn]
    print(f"  {rn:9s} " + "  ".join(f"{c}: Aδ {o[c]['Ad_over_R0']:.3f} e2 {o[c]['e2']:.3f}({o[c]['e2_proj']:+.3f}) e1 {o[c]['e1']:.3f}({o[c]['e1_proj']:+.3f}) C {o[c]['Crec']:.3f}({o[c]['Crec_proj']:+.3f})" for c in CMP))
print("結果:", res_)
