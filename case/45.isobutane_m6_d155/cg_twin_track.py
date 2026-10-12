"""plan tooling-nozzle-core-grid §4.67 (事前登録、codex 2026-10-12 mechanism-interp): Grefc の双子の列ごとの振幅 (twin_main_grefc.jsonl の col_max_dT、
壁から 50〜110 節点目の列ごとの max|ΔT|) で、固定した窓の成長と、同じ波束に追従する窓の成長を比べる。追加の計算はしない。
列ごとの振幅 a_k(n) から、窓の振幅 A(n) = √(窓の列の a_k² の平均) (固定窓・追従窓で同じ定義・同じ窓幅 W 列)。
固定窓: 島 x/r_t ∈ [−6, −3.25] の列 (W = その列数)。
追従窓: 波束の中心 c(n) を、隣り合う出力の振幅の形 (x ∈ [−9.5, −2] の列の ln a_k を平均 0・分散 1 に直したもの) の相互相関が最大になる列のずれで
  つないで追う (最大点の座標だけでは追わない)。起点は区間の最初の出力での x ∈ [−9.5, −2] の a_k² の重心。窓は [c − W/2, c + W/2]。
資格: 区間は step 1500〜2750 (6 点、1000 step 以上)、各点で max|ΔT/T| < 1e−3 (島の記録)、追従窓が x ∈ [−9.5, −2] の a² の 50 % 以上を含む、
  W を 0.67 倍・1.33 倍にしても追従窓の λ が ±20 % 以内、隣り合う出力の相互相関の最大が 0.5 以上 (形の対応がつく)。
結果 A: 追従窓の λ ≥ 1e−3 /step・R² ≥ 0.95・固定窓の λ の ±20 % 以内 → 「流入が主成分」を退け、波束自身の増幅を支持 (列数による格子差の説明までは確定しない)。
結果 B: 固定窓は増幅するが追従窓の λ が固定窓の 20 % 以下 → 流入が主成分 (固定窓の λ を使った e¹⁶/e⁸ の説明を退ける)。中間・資格不足は判定不能。
usage: python3 cg_twin_track.py → _band_ab/core_grid/twin_track.json"""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUT = HERE / "_band_ab" / "core_grid"
M = sorted({r["step"]: r for r in (json.loads(l) for l in open(OUT / "twin_main_grefc.jsonl")) if r.get("col_max_dT")}.values(), key=lambda r: r["step"])
xc = np.array(M[0]["col_x"]); ni = len(xc)
LO, HI = 1500, 2750
P = [r for r in M if LO <= r["step"] <= HI]
st = np.array([r["step"] for r in P], float)
A = np.array([r["col_max_dT"] for r in P])
reg = np.where((xc >= -9.5) & (xc <= -2.0))[0]
fixed = np.where((xc >= -6.0) & (xc <= -3.25))[0]; W0 = len(fixed)
def zs(v):
    v = np.log(np.maximum(v, 1e-300)); return (v - v.mean()) / (v.std() + 1e-300)
def fit(s, y):
    X = np.vstack([s, np.ones_like(s)]).T; c, *_ = np.linalg.lstsq(X, y, rcond=None); yh = X @ c; ss = np.sum((y - y.mean()) ** 2)
    return float(c[0]), float(1 - np.sum((y - yh) ** 2) / ss)
# 波束の中心の追跡 (列の番号)
a2 = A[0, reg] ** 2; c = [float(np.sum(reg * a2) / np.sum(a2))]; corrs = []
for n in range(1, len(P)):
    p0, p1 = zs(A[n - 1, reg]), zs(A[n, reg]); best, bs = -2, 0
    for sft in range(-5, 31):          # 下流向きのずれを主に、上流向きも少し許す
        if sft >= 0: a, b = p0[:len(p0) - sft], p1[sft:]
        else: a, b = p0[-sft:], p1[:len(p1) + sft]
        r = float(np.mean(a * b))
        if r > best: best, bs = r, sft
    corrs.append(best); c.append(c[-1] + bs)
c = np.array(c)
def amp(n, cols): return float(np.sqrt(np.mean(A[n, cols] ** 2)))
res = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.67", "steps": st.tolist(), "center_col": c.tolist(), "center_x": [float(np.interp(v, np.arange(ni), xc)) for v in c],
       "xcorr_max": corrs, "W_fixed": W0, "fixed": {}, "follow": {}}
yf = np.log([amp(n, fixed) for n in range(len(P))]); lf, r2f = fit(st, yf); res["fixed"] = {"lambda": lf, "R2": r2f, "amp": np.exp(yf).tolist()}
for scale in (0.67, 1.0, 1.33):
    W = int(round(W0 * scale)); lam_amp = []; frac = []
    for n in range(len(P)):
        lo = int(round(c[n] - W / 2)); cols = np.arange(max(lo, 0), min(lo + W, ni))
        lam_amp.append(amp(n, cols)); frac.append(float(np.sum(A[n, cols] ** 2 * np.isin(cols, reg)) / np.sum(A[n, reg] ** 2)))
    l, r2 = fit(st, np.log(lam_amp))
    res["follow"][f"W{W}"] = {"lambda": l, "R2": r2, "amp": lam_amp, "capture_frac_min": min(frac)}
rel = [r["max_dT_rel"] for r in M if LO <= r["step"] <= HI]
F = res["follow"]; base = F[f"W{W0}"]; lams = [v["lambda"] for v in F.values()]
qual = {"points": len(P) >= 5 and st[-1] - st[0] >= 1000, "linear": max(rel) < 1e-3, "capture": base["capture_frac_min"] >= 0.5,
        "W_robust": all(abs(l - base["lambda"]) <= 0.2 * abs(base["lambda"]) for l in lams), "xcorr": min(corrs) >= 0.5}
res["qualification"] = qual
if not all(qual.values()):
    out = "判定不能 (資格不足: " + ", ".join(k for k, v in qual.items() if not v) + ")"
elif base["lambda"] >= 1e-3 and base["R2"] >= 0.95 and abs(base["lambda"] - lf) <= 0.2 * abs(lf):
    out = "A: 追従窓でも固定窓と同程度に増幅 (波束自身の増幅を支持、流入が主成分ではない)"
elif lf > 0 and base["lambda"] <= 0.2 * lf:
    out = "B: 追従窓では増幅が弱い (固定窓の増大の主成分は流入)"
else:
    out = "判定不能 (中間)"
res["result"] = out
(OUT / "twin_track.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
print("中心 x:", [round(v, 2) for v in res["center_x"]], " 相互相関の最大:", [round(v, 2) for v in corrs])
print(f"固定窓 (W {W0}): λ {lf:.2e} R² {r2f:.3f}  振幅 {[f'{v:.2e}' for v in res['fixed']['amp']]}")
for k, v in F.items(): print(f"追従窓 {k}: λ {v['lambda']:.2e} R² {v['R2']:.3f} 捕捉の最小 {v['capture_frac_min']:.2f}  振幅 {[f'{a:.2e}' for a in v['amp']]}")
print("資格:", qual, " 結果:", out)
