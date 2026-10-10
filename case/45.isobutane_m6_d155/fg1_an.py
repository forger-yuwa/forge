import csv, h5py, numpy as np, sys
C45 = "/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"; C56 = "/home/ubuntu/forge-wallfit/case/56.gap_tp1187"
def vals(r, n):
    with h5py.File(f"{r}/res_{n}.h5") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"].keys()}
def hist(r):
    return [l for l in open(f"{r}/residual_history.csv")]
print("== T1 既定の経路の同一性 (ビット単位)")
for tag, base, n, (o, nw, o2) in (("FP64 case/45", C45, 20, ("run_0390_fg_id64_old", "run_0391_fg_id64_new", "run_0392_fg_id64_old2")),
                                  ("float case/56", C56, 20, ("run_0079_fg_id32_old", "run_0080_fg_id32_new", "run_0081_fg_id32_old2"))):
    try:
        A, B, A2 = vals(f"{base}/{o}", n), vals(f"{base}/{nw}", n), vals(f"{base}/{o2}", n)
        dnew = {k: int(np.count_nonzero(A[k].view(np.uint8) != B[k].view(np.uint8))) for k in A if A[k].dtype.kind == "f"}
        drep = {k: int(np.count_nonzero(A[k].view(np.uint8) != A2[k].view(np.uint8))) for k in A if A[k].dtype.kind == "f"}
        hn = sum(1 for x, y in zip(hist(f"{base}/{o}"), hist(f"{base}/{nw}")) if x != y); hr = sum(1 for x, y in zip(hist(f"{base}/{o}"), hist(f"{base}/{o2}")) if x != y)
        print(f"  {tag}: 新 vs 旧 の場の不一致バイト {sum(dnew.values())} ({[k for k, v in dnew.items() if v][:5]})、残差履歴の不一致行 {hn} | 旧 vs 旧 (再実行) {sum(drep.values())}・{hr}")
    except Exception as e:
        print(f"  {tag}: 読めない ({e})")
print("== V1・V0")
d = h5py.File(f"{C45}/run_0393_fg_v0_dump/geomab.h5"); d2 = h5py.File(f"{C45}/run_0394_fg_v0_dump2/geomab.h5"); rf = h5py.File(f"{C45}/run_0393_fg_v0_dump/geomab_ref.h5")
print("  属性: e の不一致 (旧の差との比較)", d.attrs.get("e_mismatch_vs_current_difference"), "非有限", d.attrs.get("e_nonfinite"))
pc = np.asarray(d["/mesh/plane_cells"][:]).reshape(-1, 2); nP = pc.shape[0]
NJ = 121; nC = 4719 * NJ
cc = [np.asarray(d["/mesh/cc" + c][:], dtype=np.float64) for c in "xyz"]
e64 = np.stack([np.asarray(d["/mesh/ge64_" + c][:]) for c in "xyz"], 1)
eref = np.stack([np.asarray(rf["/mesh/ge_" + c][:]) for c in "xyz"], 1)
e32 = np.stack([np.asarray(d["/mesh/ge_" + c][:], dtype=np.float64) for c in "xyz"], 1)
eold = np.stack([(np.asarray(d["/mesh/cc" + c][:])[pc[:, 1]] - np.asarray(d["/mesh/cc" + c][:])[pc[:, 0]]).astype(np.float64) for c in "xyz"], 1)
ok = np.all(np.isfinite(e64), 1) & (np.linalg.norm(e64, axis=1) > 0) & (pc[:, 0] < nC) & (pc[:, 1] < nC)
print(f"  e64 (dump の double) と REF の e の差の最大 {np.nanmax(np.abs(e64[ok] - eref[ok])):.3g}")
nrm = np.linalg.norm(e64[ok], axis=1)
rnew = np.linalg.norm(e32[ok] - e64[ok], axis=1) / nrm; rold = np.linalg.norm(eold[ok] - e64[ok], axis=1) / nrm
print(f"  V1: ‖e − e64‖/‖e64‖ 新 (e32) 最大 {rnew.max():.3g}・99% {np.percentile(rnew, 99):.3g} | 旧 (float の座標の差) 最大 {rold.max():.3g}・99% {np.percentile(rold, 99):.3g}・>1e-3 の面 {np.sum(rold > 1e-3)}")
jj = np.where(pc < nC, pc % NJ, -1)
first = ok & (np.max(jj, 1) >= 119); three = ok & (np.max(jj, 1) >= 117)
def cmp(grp, comps, label):
    Fo = np.asarray(d[f"/{grp}/old/face_flux"][:], dtype=np.float64); Fn = np.asarray(d[f"/{grp}/new/face_flux"][:], dtype=np.float64)
    Fr = np.asarray(rf[f"/{grp}/ref/face_flux"][:], dtype=np.float64); Fo2 = np.asarray(d2[f"/{grp}/old/face_flux"][:], dtype=np.float64)
    for name, cs in comps:
        for mname, m in (("第一層", first), ("3 層", three), ("全域", ok)):
            mm = m & np.all(np.isfinite(Fr[cs]), 0) & np.all(np.isfinite(Fo[cs]), 0)
            ref = np.sqrt(np.sum(Fr[cs][:, mm] ** 2))
            if ref == 0: print(f"  {label} {name} {mname}: 参照 0 で比較外"); continue
            eo = np.sqrt(np.sum((Fo[cs][:, mm] - Fr[cs][:, mm]) ** 2)) / ref; en = np.sqrt(np.sum((Fn[cs][:, mm] - Fr[cs][:, mm]) ** 2)) / ref
            ez = np.sqrt(np.sum((Fo[cs][:, mm] - Fo2[cs][:, mm]) ** 2)) / ref
            print(f"  {label} {name:6s} {mname:4s}: 旧 {eo:.3e}  新 {en:.3e}  比 {en / eo if eo > 0 else float('nan'):.3f}  再実行の差 {ez:.1e}  (面 {mm.sum()})")
cmp("viscous", [("運動量", [0, 1, 2]), ("エネルギー", [3]), ("熱", [4]), ("仕事", [5])], "粘性")
cmp("scalar", [("k", [0]), ("ω", [1])], "拡散")
