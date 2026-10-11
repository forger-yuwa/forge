"""plan tooling-nozzle-core-grid §4.50 M-v (事前登録、codex 2026-10-12 rest): 左辺の粘性スカラー対角 Dν の作用 Dνδ を、固定した δ への予測から差し引く。
入力: 新しい base (Dν の書き出し入りのバイナリで M-o の base を 1 回再現した run、linedump に scalar_visc_line・scalar_visc_off)、
  M-o の base `run_0634` (δ・左辺・Q₀・Q₁)、M-o の eval 4 本 (R₂・R₁ の Q₀・Q₁)。
Dν は節点ごとの全対象面の `viscous_diag` (2ν_eff δ/dcc) の和 (ライン面 + ライン外の面)。B0 (ITJ 5、ビット 4) では全 5 行の対角に同じ量が入る。
Dν は状態 Q₀ だけで決まるので、新しい base の Dν を run_0634 の δ に当てる (資格 a で両者の D の一致を確かめる)。
  Aδ = Lδ + Nδ − Tδ (M-o と同じ、Nδ = b₀ − b(δ))。e₁ = ΔR₁ + Aδ、e₁⁻ = e₁ − Dνδ (A から Dν を外した予測とのずれ)。
  内訳: e₁ − Dν,line δ、e₁ − Dν,off δ。e₂ も同様に併記。
資格: (a) 新しい base の D (5×5) が run_0634 と一致 (成分の最大相対差 ≤ 1e−12) かつ scalar_visc_line がビット一致、
  (b) 新しい base が run_0634 を再現 (rhs_s0 の成分別の差 ≤ 1e−8、δ の差 ≤ 1e−6)、(c) M-o の資格 (mis_judge.json) がすべて真。
結果 A: 島の ρ・ρE の両方で ‖e₁⁻‖/‖e₁‖ ≤ 0.5 かつ |⟨e₁⁻,R₀⟩|/|⟨e₁,R₀⟩| ≤ 0.5 → スカラー対角がこの補正への予測誤差の主要部分を説明する。
結果 B: 資格を満たし、上の条件を明確に満たさない (どれかの比が 0.51 を超える) → 第 1 仮説を退ける。比が 0.49〜0.51 にかかって分かれるときは判定不能。
usage: python3 cg_dnu_an.py <新しい base> <run_0634> <o2q0> <o2q1> <o1q0> <o1q1> → _band_ab/core_grid/dnu_judge.json"""
import json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
nb, base, o2q0, o2q1, o1q0, o1q1 = sys.argv[1:7]
def dumpget(run):
    meta = {}
    for ln in (HERE / run / "linedump" / "meta.txt").read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        n, r, c = ln.split(); meta[n] = (int(r), int(c))
    return lambda name: np.fromfile(HERE / run / "linedump" / f"{name}.f64", dtype=np.float64).reshape(*meta[name])
get, gnew = dumpget(base), dumpget(nb)
om = float(re.search(r"implicitRelax:\s*([0-9.eE+-]+)", (HERE / base / "solverConfig.yaml").read_text()).group(1))
S = float(json.loads((HERE / base / "prepare_info.json").read_text())["scale_m"])
nl = get("node_line"); nodes = nl[:, 0].astype(np.int64); line = nl[:, 1].astype(np.int64); N = len(nodes)
if not np.array_equal(nl, gnew("node_line")):
    raise SystemExit("書き出しの節点が違う")
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
Ad = Lmul(delta) + (b0 - bf) - (V / dt)[:, None] * delta
Dl, Do = gnew("scalar_visc_line")[:, 0], gnew("scalar_visc_off")[:, 0]
Dnu = Dl + Do
ci = [0, 1, 2, 4]; Q = ["ro", "roUx", "roUy", "roe"]
def res(run, n):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return np.stack([np.asarray(h["VALUE/res_" + q][:], float)[nodes] for q in Q], axis=1)
with h5py.File(HERE / base / "res_0.h5", "r") as h:
    xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[nodes, :2]
x = xy[:, 0] / S; jw = 159 - (nodes % 160)
R0 = res(base, 1)
dR2 = res(o2q1, 1) - res(o2q0, 1); dR1 = res(o1q1, 1) - res(o1q0, 1)
A_ = Ad[:, ci]; dl = delta[:, ci]
DnD, DlD, DoD = Dnu[:, None] * dl, Dl[:, None] * dl, Do[:, None] * dl
e1 = dR1 + A_; e2 = dR2 + A_
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
q = {}
Dn = gnew("D")
q["a_D_maxrel"] = float(np.max(np.abs(Dn - get("D")) / np.maximum(np.abs(get("D")), 1e-300)))
q["a_visc_line_equal"] = bool(np.array_equal(Dl, get("scalar_visc_line")[:, 0]))
q["a"] = q["a_D_maxrel"] <= 1e-12 and q["a_visc_line_equal"]
q["b_rhs_s0"] = {c: float(np.linalg.norm(gnew("rhs_s0")[~con, j] - b0[~con, j]) / np.linalg.norm(b0[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["b_delta"] = {c: float(np.linalg.norm(gnew("dqold_final")[~con, j] - delta[~con, j]) / np.linalg.norm(delta[~con, j])) for c, j in (("ro", 0), ("roUx", 1), ("roUy", 2), ("roe", 4))}
q["b"] = all(v <= 1e-8 for v in q["b_rhs_s0"].values()) and all(v <= 1e-6 for v in q["b_delta"].values())
mo = json.loads((OUT / "mis_judge.json").read_text())["qualification"]
q["c"] = all(mo[k] for k in ("a", "b", "c", "d", "e"))
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.50", "runs": dict(zip(("new_base", "base", "o2q0", "o2q1", "o1q0", "o1q1"), sys.argv[1:7])),
       "qualification": q, "region": {}}
for rn, m in reg.items():
    o = {"Dnu_line_share_median": float(np.median(Dl[m] / np.maximum(Dnu[m], 1e-300))),
         "Dnu_over_T_median": float(np.median(Dnu[m] / (V[m] / dt[m])))}
    for c, j in CMP.items():
        r0 = vn(R0, m, j); r02 = r0 ** 2
        def pack(e):
            return {"norm": vn(e, m, j) / r0, "proj": vdot(e, R0, m, j) / r02}
        o[c] = {"Ad": pack(-A_), "Dnu_d": pack(-DnD), "e1": pack(e1), "e1_minus": pack(e1 - DnD), "e1_minus_line": pack(e1 - DlD), "e1_minus_off": pack(e1 - DoD),
                "e2": pack(e2), "e2_minus": pack(e2 - DnD)}
        o[c]["ratio_norm"] = o[c]["e1_minus"]["norm"] / o[c]["e1"]["norm"]
        o[c]["ratio_proj"] = abs(o[c]["e1_minus"]["proj"]) / abs(o[c]["e1"]["proj"])
    out["region"][rn] = o
I = out["region"]["island"]
rs = [I[c][k] for c in ("ro", "roe") for k in ("ratio_norm", "ratio_proj")]
if not (q["a"] and q["b"] and q["c"]):
    res_ = "判定不能 (資格不足: " + ", ".join(k for k in ("a", "b", "c") if not q[k]) + ")"
elif all(r <= 0.5 for r in rs) and not any(0.49 <= r <= 0.51 for r in rs):
    res_ = "A: スカラー対角がこの補正への予測誤差の主要部分を説明する"
elif any(r > 0.51 for r in rs):
    res_ = "B: スカラー対角だけでは予測誤差の主要部分を説明できない (第 1 仮説を退ける)"
else:
    res_ = "判定不能 (比が閾値 0.5 の近くにかかる)"
out["result"] = res_
(OUT / "dnu_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("資格:", {k: q[k] for k in ("a", "b", "c")}, " D の差", f"{q['a_D_maxrel']:.1e}", " rhs_s0", {k: f"{v:.1e}" for k, v in q["b_rhs_s0"].items()}, " δ", {k: f"{v:.1e}" for k, v in q["b_delta"].items()})
for rn in ("island", "jw60_69", "jw70_79", "jw80_89", "jw90_100"):
    o = out["region"][rn]
    print(f"  {rn:9s} Dν のライン面の割合 {o['Dnu_line_share_median']:.3f}  Dν/(V/Δτ) {o['Dnu_over_T_median']:.2f}")
    for c in CMP:
        d = o[c]
        print(f"     {c:5s} −Aδ {d['Ad']['norm']:.3f}({d['Ad']['proj']:+.3f})  −Dνδ {d['Dnu_d']['norm']:.3f}({d['Dnu_d']['proj']:+.3f})  e1 {d['e1']['norm']:.3f}({d['e1']['proj']:+.3f})"
              f"  e1⁻ {d['e1_minus']['norm']:.3f}({d['e1_minus']['proj']:+.3f})  ライン面だけ {d['e1_minus_line']['norm']:.3f}  ライン外だけ {d['e1_minus_off']['norm']:.3f}  比 {d['ratio_norm']:.2f}/{d['ratio_proj']:.2f}")
print("結果:", res_)
