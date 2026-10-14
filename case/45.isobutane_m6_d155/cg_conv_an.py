"""plan tooling-nozzle-core-grid §4.52 M-c の解析と判定 (事前登録、codex 2026-10-12 dxlink): ライン外の面の対流項の、左辺の予測と右辺の応答の差。
入力: base (diag3 で M-o の base を再現、linedump に D_s1・D_s2・D_s3・D_check・rhs_s{1,2,3}final・rhs_checkfinal・scalar_visc_line/off)、
  eval 8 本 (Q₀・Q₁ × 面の種類 m = 0 全面・1 ライン外だけ・2 ライン面だけ・3 なし、いずれも convMethod 0)、比較用の run_0634 と run_0637。
記号 (M-o・M-v と同じ): L はライン行列、T = diag(V/Δτ)、δ = dqold_final、b₀ = rhs_s0、b(δ) = rhs_final、Aδ = Lδ + (b₀ − b(δ)) − Tδ、
  Dν = scalar_visc_line + scalar_visc_off (全 5 行に同じ量)。
  左辺の対流のライン外の面の作用 A_cx δ = (D − D_s1)δ − (b(δ) − b_s1(δ))、ライン面の作用 A_cl δ = (D − D_s2)δ − (Kprev δ_prev + Knext δ_next)
  (B0 では Kprev・Knext とライン外の近傍の項は対流だけ)。
  右辺: ΔR_m = R_m(Q₁) − R_m(Q₀)、ライン外の面の対流の応答 ΔR_cx = ΔR₁ − ΔR₃、ライン面 ΔR_cl = ΔR₂ − ΔR₃。
  e₁ = ΔR₀ + Aδ、E = e₁ − Dνδ、Eₓ = ΔR_cx + A_cx δ、Eₗ = ΔR_cl + A_cl δ。
資格: (a) base が run_0634 を再現 (rhs_s0 ≤ 1e−8、δ ≤ 1e−6)、(b) D_check = D・rhs_checkfinal = rhs_final (最大相対差 ≤ 1e−14)、
  (c) 左辺の加法性 (D − D_s3) = (D − D_s1) + (D − D_s2) (最大相対差 ≤ 1e−12)、
  (d) 右辺の加法性: 島の ρE・ρuₓ で ‖(ΔR₀ − ΔR₃) − (ΔR_cx + ΔR_cl)‖ ≤ 0.01‖E‖ かつ |⟨同,R₀⟩| ≤ 0.01|⟨E,R₀⟩|、
  (e) 全面の 1 次の R(Q₀) が M-o の run_0637 と一致 (島の成分別の相対差 ≤ 1e−6)、(f) base の SST 凍結 (ρk・ρω ビット一致)。
結果 A: 島の ρE・ρuₓ の両方で ‖E − Eₓ‖/‖E‖ ≤ 0.5 かつ |⟨E − Eₓ,R₀⟩|/|⟨E,R₀⟩| ≤ 0.5 → ライン外の面の対流が、この補正の向きで残る応答の差の主要部分を担う。
結果 B: どれかの比が ≥ 0.75 → この仮説を退ける。その間は判定不能。ρ・ρu_y と Eₗ は符号つきで併記。
usage: python3 cg_conv_an.py <base> <q0m0> <q1m0> <q0m1> <q1m1> <q0m2> <q1m2> <q0m3> <q1m3> <run_0634> <run_0637> → _band_ab/core_grid/conv_judge.json"""
import json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
base = sys.argv[1]; ev = sys.argv[2:10]; old, old1q0 = sys.argv[10], sys.argv[11]
def dumpget(run):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    return lambda name: np.fromfile(HERE / run / "linedump" / f"{name}.f64", dtype=np.float64).reshape(*meta[name])
get, gold = dumpget(base), dumpget(old)
om = float(re.search(r"implicitRelax:\s*([0-9.eE+-]+)", (HERE / base / "solverConfig.yaml").read_text()).group(1))
S = float(json.loads((HERE / base / "prepare_info.json").read_text())["scale_m"])
nl = get("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64); N = len(nodes)
prv = np.arange(N) - 1; nxt = np.arange(N) + 1
prv[np.r_[True, line[1:] != line[:-1]]] = -1; nxt[np.r_[line[1:] != line[:-1], True]] = -1
D5 = lambda name: get(name).reshape(N, 5, 5)
Dm, Kp, Kn = D5("D"), D5("Kprev"), D5("Knext")
dv = get("dt_vol"); dt, V = dv[:, 0], dv[:, 1]
con = np.any(get("flags_wall_iso_axis") != 0, axis=1)
mv = lambda M, d: np.einsum("nij,nj->ni", M, d)
def Kline(d):
    o = np.zeros_like(d)
    m = prv >= 0; o[m] += np.einsum("nij,nj->ni", Kp[m], d[prv[m]])
    m = nxt >= 0; o[m] += np.einsum("nij,nj->ni", Kn[m], d[nxt[m]])
    return o
delta = get("dqold_final"); b0 = get("rhs_s0"); bf = get("rhs_final")
Ld = mv(Dm, delta) - Kline(delta)
Ad = Ld + (b0 - bf) - (V / dt)[:, None] * delta
Dnu = get("scalar_visc_line")[:, 0] + get("scalar_visc_off")[:, 0]
Acx = mv(Dm - D5("D_s1"), delta) - (bf - get("rhs_s1final"))
Acl = mv(Dm - D5("D_s2"), delta) - Kline(delta)
ci = [0, 1, 2, 4]; Q = ["ro", "roUx", "roUy", "roe"]
def res(run, n=1):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return np.stack([np.asarray(h["VALUE/res_" + q][:], float)[nodes] for q in Q], axis=1)
def state(run, n, names):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], float)[nodes] for k in names}, np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
R = [res(r) for r in ev]                        # q0m0, q1m0, q0m1, q1m1, q0m2, q1m2, q0m3, q1m3
dR = [R[2 * m + 1] - R[2 * m] for m in range(4)]
R0 = res(base, 1)                               # 2 次の R(Q₀) (射影の基準、M-o・M-v と同じ)
s0, xy = state(base, 0, ["roK", "roOmega"]); s1, _ = state(base, 1, ["roK", "roOmega"])
x = xy[:, 0] / S; jw = 159 - (nodes % 160)
sel = lambda a: a[:, ci]
A_, dl, Acx_, Acl_ = sel(Ad), sel(delta), sel(Acx), sel(Acl)
e1 = dR[0] + A_; E = e1 - Dnu[:, None] * dl
dRcx = dR[1] - dR[3]; dRcl = dR[2] - dR[3]
Ex = dRcx + Acx_; El = dRcl + Acl_
addR = (dR[0] - dR[3]) - (dRcx + dRcl)
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
q = {}
q["a_rhs_s0"] = {c: float(np.linalg.norm(b0[~con, j] - gold("rhs_s0")[~con, j]) / np.linalg.norm(gold("rhs_s0")[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["a_delta"] = {c: float(np.linalg.norm(delta[~con, j] - gold("dqold_final")[~con, j]) / np.linalg.norm(gold("dqold_final")[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["a"] = all(v <= 1e-8 for v in q["a_rhs_s0"].values()) and all(v <= 1e-6 for v in q["a_delta"].values())
mr = lambda a, b: float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))
q["b_D"] = mr(D5("D_check"), Dm); q["b_rhs"] = mr(get("rhs_checkfinal"), bf)
q["b"] = q["b_D"] <= 1e-14 and q["b_rhs"] <= 1e-14
dsum = (Dm - D5("D_s1")) + (Dm - D5("D_s2")); d3 = Dm - D5("D_s3")
q["c_lhs_add"] = float(np.max(np.abs(dsum - d3)) / max(np.max(np.abs(d3)), 1e-300))
q["c"] = q["c_lhs_add"] <= 1e-12
q["d_norm"] = {c: vn(addR, isl, CMP[c]) / vn(E, isl, CMP[c]) for c in CMP}
q["d_proj"] = {c: abs(vdot(addR, R0, isl, CMP[c])) / abs(vdot(E, R0, isl, CMP[c])) for c in CMP}
q["d"] = all(q["d_norm"][c] <= 0.01 and q["d_proj"][c] <= 0.01 for c in ("roe", "roUx"))
r37 = res(old1q0)
q["e_R0_vs_0637"] = {c: vn(R[0] - r37, isl, j) / vn(r37, isl, j) for c, j in CMP.items()}
q["e"] = all(v <= 1e-6 for v in q["e_R0_vs_0637"].values())
q["f"] = bool(np.array_equal(s0["roK"], s1["roK"]) and np.array_equal(s0["roOmega"], s1["roOmega"]))
ok = all(q[k] for k in "abcdef")
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.52", "runs": {"base": base, "eval": ev, "old": old, "old_o1q0": old1q0},
       "qualification": q, "region": {}}
for rn, m in reg.items():
    o = {}
    for c, j in CMP.items():
        r0 = vn(R0, m, j); r02 = r0 ** 2; eN = vn(E, m, j); eP = vdot(E, R0, m, j)
        pk = lambda X: {"norm": vn(X, m, j) / r0, "proj": vdot(X, R0, m, j) / r02}
        o[c] = {"E": pk(E), "Ex": pk(Ex), "El": pk(El), "E_minus_Ex": pk(E - Ex), "E_minus_El": pk(E - El),
                "dRcx": pk(dRcx), "Acx": pk(-Acx_), "dRcl": pk(dRcl), "Acl": pk(-Acl_),
                "ratio_norm": vn(E - Ex, m, j) / eN, "ratio_proj": abs(vdot(E - Ex, R0, m, j)) / abs(eP),
                "ratio_norm_line": vn(E - El, m, j) / eN, "ratio_proj_line": abs(vdot(E - El, R0, m, j)) / abs(eP)}
    out["region"][rn] = o
I = out["region"]["island"]
rs = [I[c][k] for c in ("roe", "roUx") for k in ("ratio_norm", "ratio_proj")]
if not ok:
    res_ = "判定不能 (資格不足: " + ", ".join(k for k in "abcdef" if not q[k]) + ")"
elif all(r <= 0.5 for r in rs):
    res_ = "A: ライン外の面の対流が、この補正の向きで残る応答の差の主要部分を担う"
elif any(r >= 0.75 for r in rs):
    res_ = "B: ライン外の面の対流は ρE・ρuₓ の残る応答の差の主要部分ではない (仮説を退ける)"
else:
    res_ = "判定不能 (比が 0.5〜0.75)"
out["result"] = res_
(OUT / "conv_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("資格:", {k: q[k] for k in "abcdef"}, " D_check", f"{q['b_D']:.1e}", " 左辺の加法", f"{q['c_lhs_add']:.1e}",
      " 右辺の加法 (ノルム)", {k: f"{v:.1e}" for k, v in q["d_norm"].items()}, " R(Q₀) と 0637", {k: f"{v:.1e}" for k, v in q["e_R0_vs_0637"].items()})
for rn in ("island", "jw60_69", "jw70_79", "jw80_89", "jw90_100"):
    o = out["region"][rn]; print(" ", rn)
    for c in CMP:
        d = o[c]
        print(f"     {c:5s} E {d['E']['norm']:.3f}({d['E']['proj']:+.3f})  Eₓ {d['Ex']['norm']:.3f}({d['Ex']['proj']:+.3f}) [ΔR_cx {d['dRcx']['norm']:.3f}({d['dRcx']['proj']:+.3f}) −A_cxδ {d['Acx']['norm']:.3f}({d['Acx']['proj']:+.3f})]"
              f"  E−Eₓ 比 {d['ratio_norm']:.2f}/{d['ratio_proj']:.2f}  Eₗ {d['El']['norm']:.3f}({d['El']['proj']:+.3f}) E−Eₗ 比 {d['ratio_norm_line']:.2f}/{d['ratio_proj_line']:.2f}")
print("結果:", res_)
