"""plan tooling-nozzle-core-grid §4.56 M-t の解析と判定 (事前登録、codex 2026-10-12 linefaces): ライン面の熱伝導の近傍作用を固定 δ の予測に足す。
入力: base (diag5 で M-o の base を再現、linedump に rhs_s4final (右辺にライン面の熱伝導の近傍作用 K_θδ を足した組み直し) と D_s4)、
  eval 2 本 (Q₀・Q₁、全面・convMethod 0)、比較用の M-c の base `run_0649` と conv_judge.json。
記号: Aδ・Dν は M-o・M-v と同じ。e₁ = ΔR₀ + Aδ、E = e₁ − Dνδ (腕 A、M-c の E と同じ定義)。
  K_θδ = rhs_s4final − rhs_final (ライン面の熱伝導の近傍の ∂R_i/∂Q_j·δ_j、エネルギー行だけ)。予測に足す作用は Cθδ = −K_θδ。
  腕 B: E_B = E + Cθδ。
資格: (a) base の D_s4 = D (近傍作用は対角を変えない、二重加算なし) と、rhs_s4final と rhs_final の差がエネルギー行だけにあること、
  (b) base が run_0634 を再現 (rhs_s0 ≤ 1e−8、δ ≤ 1e−6)、(c) この E が M-c の E を再現 (島の各成分で ‖E − E_Mc‖ ≤ 0.01‖E‖ と射影で ≤ 0.01|⟨E,R₀⟩|、
  E_Mc は run_0649 と run_0650・0651 から同じ式で作る)、(d) 非有限なし、(e) SST 凍結。
結果 A: 島の ρE で ‖E_B‖/‖E‖ ≤ 0.5 かつ |⟨E_B,R₀⟩|/|⟨E,R₀⟩| ≤ 0.5 → ライン面の熱伝導の近傍項の欠落が ρE のずれの主要部分を担う。
結果 B: どちらか ≥ 0.75 → この仮説を退ける (熱伝導係数・μt の状態依存や有限振幅の非線形の確定ではない)。その間・資格不足は判定不能。
usage: python3 cg_therm_an.py <base> <q0m0> <q1m0> <run_0634> <run_0649> <run_0650> <run_0651> → _band_ab/core_grid/therm_judge.json"""
import json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
base, q0, q1, old, mcb, mc0, mc1 = sys.argv[1:8]
def dumpget(run):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    return lambda name: np.fromfile(HERE / run / "linedump" / f"{name}.f64", dtype=np.float64).reshape(*meta[name])
S = float(json.loads((HERE / base / "prepare_info.json").read_text())["scale_m"])
ci = [0, 1, 2, 4]; Q = ["ro", "roUx", "roUy", "roe"]
def build(b, e0, e1):
    """b の書き出しと e0・e1 (Q₀・Q₁ の 1 次の全面の評価) から E を作る"""
    g = dumpget(b)
    nl = g("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64); N = len(nodes)
    prv = np.arange(N) - 1; nxt = np.arange(N) + 1
    prv[np.r_[True, line[1:] != line[:-1]]] = -1; nxt[np.r_[line[1:] != line[:-1], True]] = -1
    D5 = lambda n: g(n).reshape(N, 5, 5)
    Dm, Kp, Kn = D5("D"), D5("Kprev"), D5("Knext")
    dv = g("dt_vol"); dt, V = dv[:, 0], dv[:, 1]
    def Lmul(d):
        o = np.einsum("nij,nj->ni", Dm, d)
        m = prv >= 0; o[m] -= np.einsum("nij,nj->ni", Kp[m], d[prv[m]])
        m = nxt >= 0; o[m] -= np.einsum("nij,nj->ni", Kn[m], d[nxt[m]])
        return o
    delta = g("dqold_final"); b0 = g("rhs_s0"); bf = g("rhs_final")
    Ad = Lmul(delta) + (b0 - bf) - (V / dt)[:, None] * delta
    Dnu = g("scalar_visc_line")[:, 0] + g("scalar_visc_off")[:, 0]
    def res(run, n=1):
        with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
            return np.stack([np.asarray(h["VALUE/res_" + q][:], float)[nodes] for q in Q], axis=1)
    dR = res(e1) - res(e0)
    E = dR + Ad[:, ci] - Dnu[:, None] * delta[:, ci]
    return g, nodes, V, E, res(b, 1), bf, Dm
g, nodes, V, E, R0, bf, Dm = build(base, q0, q1)
_, nodes_mc, _, E_mc, _, _, _ = build(mcb, mc0, mc1)
if not np.array_equal(nodes, nodes_mc):
    raise SystemExit("書き出しの節点が M-c と違う")
N = len(nodes)
con = np.any(g("flags_wall_iso_axis") != 0, axis=1)
Kth = g("rhs_s4final") - bf
Ct = -Kth[:, ci]
EB = E + Ct
with h5py.File(HERE / base / "res_0.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
    k0, w0 = np.asarray(h["VALUE/roK"][:])[nodes], np.asarray(h["VALUE/roOmega"][:])[nodes]
with h5py.File(HERE / base / "res_1.h5", "r") as h:
    k1, w1 = np.asarray(h["VALUE/roK"][:])[nodes], np.asarray(h["VALUE/roOmega"][:])[nodes]
x = xy[:, 0] / S; jw = 159 - (nodes % 160)
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
isl = reg["island"]; q = {}
q["a_D_s4_equal"] = bool(np.array_equal(g("D_s4").reshape(N, 5, 5), Dm))
q["a_rows_0to3_unchanged"] = bool(np.array_equal(g("rhs_s4final")[:, :4], bf[:, :4]))
q["a"] = q["a_D_s4_equal"] and q["a_rows_0to3_unchanged"]
go = dumpget(old)
q["b_rhs_s0"] = {c: float(np.linalg.norm(g("rhs_s0")[~con, j] - go("rhs_s0")[~con, j]) / np.linalg.norm(go("rhs_s0")[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["b_delta"] = {c: float(np.linalg.norm(g("dqold_final")[~con, j] - go("dqold_final")[~con, j]) / np.linalg.norm(go("dqold_final")[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["b"] = all(v <= 1e-8 for v in q["b_rhs_s0"].values()) and all(v <= 1e-6 for v in q["b_delta"].values())
q["c_norm"] = {c: vn(E - E_mc, isl, j) / vn(E, isl, j) for c, j in CMP.items()}
q["c_proj"] = {c: abs(vdot(E - E_mc, R0, isl, j)) / max(abs(vdot(E, R0, isl, j)), 1e-300) for c, j in CMP.items()}
q["c"] = all(q["c_norm"][c] <= 0.01 for c in CMP) and q["c_proj"]["roe"] <= 0.01
q["d"] = bool(np.all(np.isfinite(EB)) and np.all(np.isfinite(Kth)))
q["e"] = bool(np.array_equal(k0, k1) and np.array_equal(w0, w1))
ok = all(q[k] for k in "abcde")
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.56", "runs": {"base": base, "q0": q0, "q1": q1, "old": old, "mc": [mcb, mc0, mc1]},
       "qualification": q, "region": {}}
for rn, m in reg.items():
    o = {}
    for c, j in CMP.items():
        r0 = vn(R0, m, j); r02 = r0 ** 2
        pk = lambda X: {"norm": vn(X, m, j) / r0, "proj": vdot(X, R0, m, j) / r02}
        o[c] = {"E": pk(E), "Ct": pk(Ct), "EB": pk(EB),
                "ratio_norm": vn(EB, m, j) / vn(E, m, j), "ratio_proj": abs(vdot(EB, R0, m, j)) / max(abs(vdot(E, R0, m, j)), 1e-300)}
    out["region"][rn] = o
I = out["region"]["island"]["roe"]
if not ok:
    res_ = "判定不能 (資格不足: " + ", ".join(k for k in "abcde" if not q[k]) + ")"
elif I["ratio_norm"] <= 0.5 and I["ratio_proj"] <= 0.5:
    res_ = "A: ライン面の熱伝導の近傍項の欠落が ρE のずれの主要部分を担う"
elif I["ratio_norm"] >= 0.75 or I["ratio_proj"] >= 0.75:
    res_ = "B: ライン面の熱伝導の近傍項の欠落は ρE のずれの主要部分ではない (仮説を退ける)"
else:
    res_ = "判定不能 (比が 0.5〜0.75)"
out["result"] = res_
(OUT / "therm_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("資格:", {k: q[k] for k in "abcde"}, " E の再現 (ノルム)", {k: f"{v:.1e}" for k, v in q["c_norm"].items()})
for rn in ("island", "jw60_69", "jw70_79", "jw80_89", "jw90_100"):
    d = out["region"][rn]["roe"]
    print(f"  {rn:9s} ρE: E {d['E']['norm']:.3f}({d['E']['proj']:+.3f})  Cθδ {d['Ct']['norm']:.3f}({d['Ct']['proj']:+.3f})  E_B {d['EB']['norm']:.3f}({d['EB']['proj']:+.3f})  比 {d['ratio_norm']:.2f}/{d['ratio_proj']:.2f}")
print("結果:", res_)
