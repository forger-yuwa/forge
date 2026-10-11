"""plan architecture-float-state-double-geometry §6.44 の前段の関門 (事前登録、codex diagnose 2026-10-11): FP64 で pRef が作用素を実質変えないこと。
通れば終了コード 0 (rpref.sh が Fp に進む)、通らなければ 1 (前提不成立。幅を広げて救済しない)。結果は _band_ab/cold_pair/cp_gate.json。
1. 作用素 (幾何): Cp1 の hoop.h5 (FORGE_DIAG_HOOP_CLOSURE、FP64) から、各実 CV で
   g = sqrt(sum64_x² + (sum64_y − T_y)²) / (abs64_x + abs64_y + A)   (T_y = A_planar、A = A_planar)
   pRef による残差の変化は ΔR = pRef·(Σ±S_x, Σ±S_y − A_planar) なので、その相対の大きさ。全 CV で g ≤ 1e-12。
   (分母は Σ‖W‖ + ‖A‖ の代わりに Σ|S_x| + Σ|S_y| + A を使う。最大で √2 倍ゆるい)
2. 有限期間: mean(a_Cp) − mean(a_C) の 95 % Welch 区間 ⊂ [−0.002, +0.002]。step 500・1,000・1,500・2,000 の Cp_i − C_j の全 9 組で、
   温度差 RMS ≤ 0.05 K かつ最大 |ΔT| ≤ 1 K を、上流域 (x/r_t < −2 の内点) と全域 (全節点) で別々に満たす。"""
import json, math, sys
from pathlib import Path
import numpy as np, h5py
from scipy import stats
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
C = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3")
CP = ("run_0700_rp_cp1", "run_0701_rp_cp2", "run_0702_rp_cp3")
REF = {"A": "run_0487_tr_fp64", "B": "run_0488_tr_f32", "Bp": "run_0489_tr_f32b"}
S = (1000, 1500, 2000)
OUT = {"fail": [], "closure": {}, "a": {}, "pairs": {}}
def fail(m): OUT["fail"].append(m); print("  [関門で不成立]", m)
# 1. 閉性
hp = HERE / CP[0] / "hoop.h5"
if not hp.exists(): fail("Cp1 の hoop.h5 が無い")
else:
    with h5py.File(hp, "r") as h:
        n = len(h["cells/sum64_x"])
        sx, sy = np.asarray(h["cells/sum64_x"][:], float), np.asarray(h["cells/sum64_y"][:], float)
        ax, ay = np.asarray(h["cells/abs64_x"][:], float), np.asarray(h["cells/abs64_y"][:], float)
        A = np.asarray(h["cells/A_planar"][:], float)[:n]
        attrs = {k: (h.attrs[k].item() if hasattr(h.attrs[k], "item") else str(h.attrs[k])) for k in h.attrs}
    g = np.hypot(sx, sy - A) / np.maximum(ax + ay + np.abs(A), 1e-300)
    OUT["closure"] = dict(n=int(n), g_max=float(g.max()), g_p99=float(np.percentile(g, 99)), n_over=int((g > 1e-12).sum()), attrs={k: str(v)[:80] for k, v in attrs.items()})
    print(f"閉性: g の最大 {g.max():.3e}、99 % {np.percentile(g, 99):.3e}、1e-12 を超える CV {int((g > 1e-12).sum())} / {n}")
    if g.max() > 1e-12: fail(f"閉性 g の最大 {g.max():.3e} > 1e-12")
# 2. 有限期間
with h5py.File(HERE / REF["A"] / "nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
ni, nj, Sc = XC.mesh_info(xy0.size // 3); X = (xy0.reshape(-1, 3)[:, 0] / Sc).reshape(ni, nj)
I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
DOM = ((X < -2) & (I >= 1) & (J >= 1) & (J <= nj - 2)).ravel(); ALL = np.ones(DOM.size, bool)
def T(run, s):
    with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: return np.asarray(h["VALUE/T"][:], float)
TA0 = T(REF["A"], 0); z = lambda run, s: (T(run, s) - T(run, 0)) / TA0  # noqa: E731
V = {s: 0.5 * (z(REF["B"], s) + z(REF["Bp"], s)) - z(REF["A"], s) for s in S}
a = lambda run: float(np.mean([np.dot(V[s][DOM], z(run, s)[DOM]) / np.dot(V[s][DOM], V[s][DOM]) for s in S]))  # noqa: E731
aC, aCp = [a(r) for r in C], [a(r) for r in CP]
x, c = np.array(aCp), np.array(aC); e = x.mean() - c.mean(); vx, vc = x.var(ddof=1) / 3, c.var(ddof=1) / 3; se = math.sqrt(vx + vc)
if se == 0: L = U = e
else:
    df = (vx + vc) ** 2 / ((vx ** 2 / 2 if vx else 0) + (vc ** 2 / 2 if vc else 0)); t = float(stats.t.ppf(0.975, df)); L, U = e - t * se, e + t * se
OUT["a"] = dict(a_C=aC, a_Cp=aCp, e=e, L=L, U=U)
print(f"a: C {[round(v, 6) for v in aC]}、Cp {[round(v, 6) for v in aCp]}、差 {e:+.6f} 区間 [{L:+.6f}, {U:+.6f}]")
if not (-0.002 <= L and U <= 0.002): fail(f"a の差の区間 [{L:+.6f}, {U:+.6f}] が [−0.002, 0.002] に入らない")
for s in (500, 1000, 1500, 2000):
    for cp in CP:
        for cr in C:
            d = T(cp, s) - T(cr, s)
            for nm, m in (("上流域", DOM), ("全域", ALL)):
                rms, mx = float(np.sqrt(np.mean(d[m] ** 2))), float(np.max(np.abs(d[m])))
                OUT["pairs"][f"{s}|{cp}|{cr}|{nm}"] = [rms, mx]
                if rms > 0.05 or mx > 1.0: fail(f"step {s} {cp}−{cr} {nm}: RMS {rms:.4f} K・最大 {mx:.3f} K")
worst = {nm: max(v[1] for k, v in OUT["pairs"].items() if k.endswith(nm)) for nm in ("上流域", "全域")}
worst_rms = {nm: max(v[0] for k, v in OUT["pairs"].items() if k.endswith(nm)) for nm in ("上流域", "全域")}
print(f"Cp − C の温度差: RMS の最大 {worst_rms}、|ΔT| の最大 {worst}")
OUT["verdict"] = "関門を通った" if not OUT["fail"] else "前提不成立: " + "; ".join(OUT["fail"][:3])
print("== 関門 §6.44:", OUT["verdict"])
(HERE / "_band_ab/cold_pair/cp_gate.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
sys.exit(0 if not OUT["fail"] else 1)
