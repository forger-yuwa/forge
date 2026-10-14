"""plan tooling-nozzle-core-grid §4.73 の判定 (事前登録、codex 2026-10-12 gref-track-interp): Grefc の後半の増幅の鈍化が摂動の振幅に依存するか。
入力: A = 既存の Grefc (a = 1e−7、twin_main_grefc.jsonl、run_0678 の perturb.json)、B = 新しい Grefc (a = 1e−10、twin_main_grefc_small.jsonl、run_0682 の perturb.json)。
正規化: 各 run の実現振幅 (perturb.json の max_rel) で割る。
見る量: (1) 正規化した追従窓の振幅 F(n) = √(窓の列の (a_k/a)² の平均)、窓幅 W = 島の列数、中心は各系列で §4.67 と同じ方法 (形の相互相関) で追う
  (step 1500 の x ∈ [−9.5, −2] の a² の重心から始め、6000 まで)。(2) 正規化した列ごとの分布の相対差 ‖p_B − p_A‖/‖p_A‖ (x ∈ [−9.5, −2] の列、p = a_k/a)。
資格: B の全出力で非有限なし、B の島の D_T ≥ 100 L (L は §4.66 の Grefc の雑音の床 2.87e−11 K を振幅比で使う必要はない — B の信号そのものと比べる)、
  B の追跡の相互相関の最大 ≥ 0.5 (1500〜4750)、1500〜2500 の 5 点で F_B/F_A が 0.8〜1.2。
結果 A: 3500〜4750 の 6 点すべてで F_B/F_A が 0.8〜1.2 かつ分布の相対差 ≤ 0.2 → 鈍化を有限振幅の飽和で説明する仮説を退ける (増幅域を出たかは未確定)。
結果 B: 1500〜2500 は再現し、3500 以降に F_B/F_A ≥ 2 の点が 3 出力以上続く → 振幅依存の抑制 (有限振幅の飽和) を支持。中間・資格不足は判定不能。
usage: python3 cg_twin_amp.py → _band_ab/core_grid/twin_amp.json"""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
def rd(tag):
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(OUT / f"twin_{tag}.jsonl")) if r.get("col_max_dT")}.values(), key=lambda r: r["step"])
aA = json.loads((HERE / "run_0678_twc50_pert" / "perturb.json").read_text())["max_rel"]
aB = json.loads((HERE / "run_0682_twc50s_pert" / "perturb.json").read_text())["max_rel"]
MA, MB = rd("main_grefc"), rd("main_grefc_small")
xc = np.array(MA[0]["col_x"]); ni = len(xc)
if not np.allclose(xc, np.array(MB[0]["col_x"])): raise SystemExit("列の x が A と B で違う")
reg = np.where((xc >= -9.5) & (xc <= -2.0))[0]; W = int(np.count_nonzero((xc >= -6.0) & (xc <= -3.25)))
def zs(v):
    v = np.log(np.maximum(v, 1e-300)); return (v - v.mean()) / (v.std() + 1e-300)
def track(M, a):
    S = [r for r in M if 1500 <= r["step"] <= 6000]; st = [r["step"] for r in S]; A = np.array([r["col_max_dT"] for r in S]) / a
    a2 = A[0, reg] ** 2; c = [float(np.sum(reg * a2) / np.sum(a2))]; corr = []
    for n in range(1, len(S)):
        p0, p1 = zs(A[n - 1, reg]), zs(A[n, reg]); best, bs = -2, 0
        for sft in range(-10, 31):
            if sft >= 0: u, v = p0[:len(p0) - sft], p1[sft:]
            else: u, v = p0[-sft:], p1[:len(p1) + sft]
            r = float(np.mean(u * v))
            if r > best: best, bs = r, sft
        corr.append(best); c.append(c[-1] + bs)
    F = []
    for n in range(len(S)):
        lo = int(round(c[n] - W / 2)); cols = np.arange(max(lo, 0), min(lo + W, ni)); F.append(float(np.sqrt(np.mean(A[n, cols] ** 2))))
    return st, np.array(F), np.array(c), corr, A
sA, FA, cA, corrA, PA = track(MA, aA); sB, FB, cB, corrB, PB = track(MB, aB)
if sA != sB: raise SystemExit("出力の step が A と B で違う")
st = np.array(sA); ratio = FB / FA
prof = [float(np.linalg.norm(PB[n, reg] - PA[n, reg]) / np.linalg.norm(PA[n, reg])) for n in range(len(st))]
nfB = sum(r["nonfinite"] for r in MB)
early = (st >= 1500) & (st <= 2500); late = (st >= 3500) & (st <= 4750)
qual = {"nonfinite_B": nfB == 0, "xcorr_B": min(c for s, c in zip(st[1:], corrB) if s <= 4750) >= 0.5,
        "early_match": bool(np.all((ratio[early] >= 0.8) & (ratio[early] <= 1.2))), "signal_B": min(r["island_D_T"] for r in MB if r["step"] >= 1500) >= 100 * 2.87e-11}
if not all(qual.values()):
    res = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif np.all((ratio[late] >= 0.8) & (ratio[late] <= 1.2)) and np.all(np.array(prof)[late] <= 0.2):
    res = "A: 小さい摂動でも後半の鈍化を再現 (有限振幅の飽和の仮説を退ける)"
else:
    big = ratio >= 2; run, best = 0, 0
    for k in range(len(st)):
        if st[k] >= 3500 and big[k]: run += 1; best = max(best, run)
        else: run = 0
    res = "B: 後半に小さい摂動の正規化振幅が 2 倍以上 (振幅依存の抑制を支持)" if best >= 3 else "判定不能 (中間)"
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.73", "a_A": aA, "a_B": aB, "W": W, "steps": st.tolist(), "F_A": FA.tolist(), "F_B": FB.tolist(),
       "ratio_B_over_A": ratio.tolist(), "profile_rel_diff": prof, "center_x_A": [float(np.interp(v, np.arange(ni), xc)) for v in cA],
       "center_x_B": [float(np.interp(v, np.arange(ni), xc)) for v in cB], "xcorr_A": corrA, "xcorr_B": corrB, "qualification": qual, "result": res}
(OUT / "twin_amp.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print(f"実現振幅 A {aA:.3e}  B {aB:.3e}  W {W}")
for k in range(len(st)):
    print(f"  {int(st[k]):5d}  F_A {FA[k]:.3e}  F_B {FB[k]:.3e}  B/A {ratio[k]:.3f}  分布の差 {prof[k]:.3f}  中心 A {out['center_x_A'][k]:.2f} B {out['center_x_B'][k]:.2f}")
print("資格:", qual, " 結果:", res)
