"""plan condensation-two-phase-default §5.1 #5p (G2) の判定 (事前登録 c87d66cb のとおり)。
  python3 g2_judge.py <case/16 dir>
腕: P = run_0561–0563 (ON 生産)、C = run_0564–0566 (ON cfl2)、O = run_0567–0569 (OFF)。
各 run: 7 量 STEADY (QUASISTEADY_7q.txt)、非有限・負値 0 (series の検査列)、全残差列 RISING 0 (CONVERGENCE_VERDICT.txt)、ON は corr-gate (commit/floor/clamp) ≤ κ。
比較値 = 判定窓 (系列の末尾 50 %) の平均、σ = 3 本の窓平均の標本標準偏差、s_d = √(s_a²/3 + s_b²/3)。
レシピ感度: abs(ON−OFF) > 3·s_{ON−OFF} の量 → abs(P−C) ≤ 0.1·abs(ON−OFF) (ON は C を使う)、区間 [abs(P−C) ± 2 s_PC] が線をまたげば判定不能。
それ以外: abs(P−C) ≤ max(3·s_PC, 絶対許容) (onset 0.01 mm のみ絶対許容)。
"""
import sys, glob, os, re, numpy as np, pandas as pd

Q = ["onset_c_g1e4_mm", "g_exit_mw", "pw_mean_x10", "pw42", "pw52", "Tw_mean_x10_K", "T_cond0_vmean_K"]
ABS_TOL = {"onset_c_g1e4_mm": 0.01}
KAPPA = 4.7683716e-07
C = sys.argv[1]
arms = {"P": range(561, 564), "C": range(564, 567), "O": range(567, 570)}

def rundir(n):
    m = glob.glob(os.path.join(C, f"run_0{n}_g2_*"))
    return m[0] if m else None

ok_all = True; means = {a: [] for a in arms}
for a, ns in arms.items():
    for n in ns:
        d = rundir(n); name = os.path.basename(d)
        s = pd.read_csv(os.path.join(d, "twophase_series.csv")).sort_values("step")
        w = s[s.step >= s.step.max() * 0.5]
        means[a].append(w[Q].mean().values)
        qs = open(os.path.join(d, "QUASISTEADY_7q.txt")).read()
        steady = "ALL STEADY" in qs
        nf = int(s[["nonfinite_cons"]].max().iloc[0]); neg = int(s[["neg_rv", "neg_rg", "neg_rQ"]].max().max())
        cv = open(os.path.join(d, "CONVERGENCE_VERDICT.txt")).read()
        rising = len(re.findall(r"RISING", cv))
        gate = "n/a"
        if a in "PC":
            log = open(os.path.join(d, "forge_run.log"), errors="replace").read().split("\n")
            # codex 2026-10-04 (G3 設計の諮問) Major: FINAL 行の直後の段ブロックだけを読み、段の欠落はエラー、
            # 射影・正式 VERDICT・非有限も記録する。合否は事前登録 (#4j (D)) どおり commit/floor/clamp ≤ κ、射影は記録のみ。
            fi = [k for k, l in enumerate(log) if l.startswith("[twophase-corr-gate] FINAL")]
            gate = "missing"
            if fi:
                blk = log[fi[-1]:fi[-1] + 8]
                formal = re.search(r"VERDICT: (\w+)", blk[0]); nonf = re.search(r"nonfinite (\d+)", blk[0])
                st = {}
                for l in blk[1:]:
                    m = re.search(r"stage (\S+)", l)
                    if m: st[m.group(1)] = [float(x) for x in re.findall(r"rho\w+ ([0-9.eE+-]+)", l)]
                need = ("commit", "passive_floor", "clamp")
                if all(s in st and len(st[s]) == 6 for s in need) and nonf and int(nonf.group(1)) == 0:
                    gated = max(max(st[s]) for s in need)
                    proj = max(st.get("projection", [float("nan")]))
                    gate = ("ok" if gated <= KAPPA else "FAIL") + f" (commit/floor/clamp max {gated:.2e}; projection {proj:.2e} recorded; formal VERDICT {formal.group(1) if formal else '?'})"
                else:
                    gate = "FAIL (stage missing or non-finite)"
        good = steady and nf == 0 and neg == 0 and rising == 0 and (gate == "n/a" or gate.startswith("ok"))
        ok_all &= good
        print(f"{a} {name}: STEADY {steady} | nonfinite {nf} neg {neg} | RISING {rising} | corr-gate {gate} -> {'ok' if good else 'NOT OK'}")

M = {a: np.array(v) for a, v in means.items()}
mean = {a: M[a].mean(0) for a in M}; sd = {a: M[a].std(0, ddof=1) for a in M}
def sdiff(a, b): return np.sqrt(sd[a] ** 2 / 3 + sd[b] ** 2 / 3)
print("\nquantity | P mean (sd) | C mean (sd) | O mean (sd) | ON-OFF (C-O) ± s | P-C ± s | judgement")
undecided = False; fail = False
for i, q in enumerate(Q):
    on_off = mean["C"][i] - mean["O"][i]; s_oo = sdiff("C", "O")[i]
    pc = mean["P"][i] - mean["C"][i]; s_pc = sdiff("P", "C")[i]
    if abs(on_off) > 3 * s_oo:
        lim = 0.1 * abs(on_off)
        lo, hi = abs(pc) - 2 * s_pc, abs(pc) + 2 * s_pc
        if hi <= lim: j = "PASS (significant ON-OFF; |P-C| <= 0.1|ON-OFF|)"
        elif lo > lim: j = "FAIL (recipe sensitivity > 0.1|ON-OFF|)"; fail = True
        else: j = "UNDECIDED (interval straddles the 0.1|ON-OFF| line)"; undecided = True
    else:
        lim = max(3 * s_pc, ABS_TOL.get(q, 0.0))
        if abs(pc) <= lim: j = f"PASS (ON-OFF not significant; |P-C| <= {lim:.3g})"
        else: j = f"FAIL (|P-C| > {lim:.3g})"; fail = True
    print(f"{q} | {mean['P'][i]:.6g} ({sd['P'][i]:.2g}) | {mean['C'][i]:.6g} ({sd['C'][i]:.2g}) | {mean['O'][i]:.6g} ({sd['O'][i]:.2g}) | {on_off:+.4g} ± {s_oo:.2g} | {pc:+.4g} ± {s_pc:.2g} | {j}")
v = "PASS" if (ok_all and not fail and not undecided) else ("FAIL" if (fail or not ok_all) else "UNDECIDED")
print(f"\nVERDICT G2: {v} (per-run checks all ok: {ok_all})")
sys.exit(0 if v == "PASS" else (1 if v == "FAIL" else 2))
