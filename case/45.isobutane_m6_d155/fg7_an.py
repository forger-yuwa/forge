"""plan architecture-float-state-double-geometry §6.12 (事前登録) の集計。fg7.sh の後に AWS の case/45 で回す。
ΔR = R32 − R64 (V3 の run_0442 − run_0444、res_1 の commit 前の残差) のうち、拡散の面の流束の入力の丸めで説明できる分
C = Σ±(F_B − F_A) (A = FP64 が自分で作った入力、B = float の入力を移したもの) を、本番と同じ符号 (ic0 に +F、ic1 に −F) で集める。
r = ‖ΔR − C‖₂/‖ΔR‖₂ (壁際 j ≥ 117)。ω が主判定、k は記録。必要な内部面の欠損・非有限は判定不能。結果は fg7_verdict.json。"""
import json, h5py, numpy as np
from pathlib import Path
C45 = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); NJ = 121
class Undecidable(Exception): pass
def res1(r, k):
    with h5py.File(C45 / r / "res_1.h5", "r") as h:
        if "VALUE/" + k not in h: raise Undecidable(f"{r} に {k} が無い")
        a = np.asarray(h["VALUE/" + k][:], np.float64)
    if not np.all(np.isfinite(a)): raise Undecidable(f"{r} の {k} に非有限")
    return a
def flux(path, grp, arm):
    with h5py.File(path, "r") as h:
        key = f"/{grp}/{arm}/face_flux"
        if key not in h: raise Undecidable(f"{path}: {key} が無い")
        F = np.asarray(h[key][:], np.float64)
        pc = np.asarray(h["/mesh/plane_cells"][:]).reshape(-1, 2) if "/mesh/plane_cells" in h else None
        nN = int(h.attrs["nNormalPlanes"]) if "nNormalPlanes" in h.attrs else None
    return F, pc, nN
V = {}
try:
    FA, pc, nN = flux(C45 / "run_0448_v3_self64a/geomab.h5", "scalar", "new")
    FA2, pc2, _ = flux(C45 / "run_0449_v3_self64b/geomab.h5", "scalar", "new")
    FB, _, _ = flux(C45 / "run_0446_v3_dump32/geomab_ref.h5", "scalar", "ref")
    FB2, _, _ = flux(C45 / "run_0450_v3_ref64dump2/geomab_src_ref.h5", "scalar", "ref")
    F32, pc32, nN32 = flux(C45 / "run_0446_v3_dump32/geomab.h5", "scalar", "new")
    if pc is None or nN is None: raise Undecidable("plane_cells か nNormalPlanes が無い")
    if not (np.array_equal(pc, pc2) and np.array_equal(pc, pc32) and nN == nN32): raise Undecidable("接続が腕で違う")
    nP = pc.shape[0]
    for name, F in (("A", FA), ("A2", FA2), ("B", FB), ("B2", FB2), ("F32", F32)):
        if F.size != 2 * nP: raise Undecidable(f"{name} の大きさ {F.size} ≠ 2·nPlanes {2 * nP}")
        bad = ~np.isfinite(F.reshape(2, nP)[:, :nN])
        if bad.any(): raise Undecidable(f"{name}: 内部面に非有限 {int(bad.sum())}")
    n = res1("run_0444_v3_ref64a", "res_roOmega").size; jj = np.arange(n) % NJ; W = jj >= 117
    a, b = pc[:nN, 0], pc[:nN, 1]; ma, mb = a < n, b < n
    def scatter(dF):
        acc = np.zeros(n); np.add.at(acc, a[ma], dF[:nN][ma]); np.add.at(acc, b[mb], -dF[:nN][mb]); return acc
    for s, (k, label) in enumerate((("res_roK", "k"), ("res_roOmega", "ω"))):
        dR = res1("run_0442_v3_new32a", k) - res1("run_0444_v3_ref64a", k)
        n32 = res1("run_0443_v3_new32b", k) - res1("run_0442_v3_new32a", k); n64 = res1("run_0445_v3_ref64b", k) - res1("run_0444_v3_ref64a", k)
        fa, fa2, fb, fb2, f32 = (F.reshape(2, nP)[s] for F in (FA, FA2, FB, FB2, F32))
        C = scatter(fb - fa); C2 = scatter(fb2 - fa2); Call = scatter(f32 - fa)
        nd = np.linalg.norm(dR[W])
        if nd == 0: raise Undecidable(f"{label}: ‖ΔR‖ = 0")
        r = np.linalg.norm(dR[W] - C[W]) / nd; rall = np.linalg.norm(dR[W] - Call[W]) / nd
        rep = max(np.linalg.norm(n32[W]), np.linalg.norm(n64[W]), np.linalg.norm(C2[W] - C[W])) / nd
        verdict = ("判定不能 (再実行の差が ‖ΔR‖ の 10 % を超える)" if rep > 0.1 else
                   "支持 (入力の丸めが過半を説明)" if r <= 0.5 else "棄却" if r >= 0.9 else "判別不能")
        tag = "主判定" if label == "ω" else "記録"
        print(f"  [{tag}] {label} 壁際 j≥117: r = ‖ΔR − C‖/‖ΔR‖ = {r:.3f}、‖C‖/‖ΔR‖ = {np.linalg.norm(C[W]) / nd:.3f}、"
              f"演算の差まで含めた r_all = {rall:.3f}、再実行の差 / ‖ΔR‖ = {rep:.3e} → {verdict}")
        rg = {"全域": np.ones(n, bool)}
        for nm, m in rg.items():
            ndg = np.linalg.norm(dR[m])
            print(f"     (記録) {label} {nm}: r = {np.linalg.norm(dR[m] - C[m]) / ndg:.3f}、r_all = {np.linalg.norm(dR[m] - Call[m]) / ndg:.3f}")
        V[label] = dict(r=float(r), r_all=float(rall), C_over_dR=float(np.linalg.norm(C[W]) / nd), rerun=float(rep), verdict=verdict, primary=(label == "ω"))
except Undecidable as e:
    V["error"] = str(e); print(f"  [判定不能] {e}")
json.dump(V, open(C45 / "fg7_verdict.json", "w"), ensure_ascii=False, indent=1)
