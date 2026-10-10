"""plan architecture-float-state-double-geometry §6.3 の 3' (closure を使う経路の一様な静止場) の集計。fg2d.sh の後に AWS で回す。"""
import json, h5py, numpy as np
from pathlib import Path
C45 = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); P0 = 1.0e5; NJ = 121
Ap = {32: np.asarray(h5py.File(C45 / "run_0396_fg2_dump32/stage2.h5", "r")["/axisym/A_planar"][:], np.float64),
      64: np.asarray(h5py.File(C45 / "run_0397_fg2_dump64/stage2.h5", "r")["/axisym/A_planar"][:], np.float64)}
def stats(run, prec):
    with h5py.File(C45 / run / "res_1.h5", "r") as h: q = np.asarray(h["VALUE/res_roUy"][:], np.float64)
    with h5py.File(C45 / run / "res_10.h5", "r") as h:
        V = h["VALUE"]; u = np.sqrt(sum(np.asarray(V[k][:], np.float64) ** 2 for k in ("Ux", "Uy", "Uz")))
        bad = sum(int(np.count_nonzero(~np.isfinite(np.asarray(V[k][:])))) for k in V.keys())
    q = q / (P0 * Ap[prec][:q.size]); j = np.arange(q.size) % NJ
    im, iu = int(np.argmax(np.abs(q))), int(np.argmax(u))
    prof = {jj: float(np.max(np.abs(q[j == jj]))) for jj in (1, 2, 3, 4, 5, 60, 118, 119)}
    print(f"  {run}: 最大 {np.max(np.abs(q)):.3e} (i {im // NJ}・j {im % NJ})・RMS {np.sqrt(np.mean(q ** 2)):.3e}、10 step 後の最大 |u| {u.max():.3e} (i {iu // NJ}・j {iu % NJ})、非有限 {bad}")
    print(f"     j ごとの最大: " + "  ".join(f"j{k} {v:.2e}" for k, v in prof.items()))
    return dict(qmax=float(np.max(np.abs(q))), qrms=float(np.sqrt(np.mean(q ** 2))), umax=float(u.max()), bad=bad)
print("== 3' closure を使う経路 (hoopAreaFromClosure 1) の一様な静止場")
S = {r: stats(r, 64 if r.endswith("64") else 32) for r in ("run_0405_fg2_stillhc_old32", "run_0406_fg2_stillhc_new32", "run_0407_fg2_stillhc_old32b", "run_0408_fg2_stillhc_new64")}
o, n, ob = S["run_0405_fg2_stillhc_old32"], S["run_0406_fg2_stillhc_new32"], S["run_0407_fg2_stillhc_old32b"]
ok_r = n["qmax"] <= 1.1 * o["qmax"] and n["qrms"] <= 1.1 * o["qrms"]; ok_u = n["umax"] <= 1.1 * o["umax"]
print(f"  [{'PASS' if ok_r else 'FAIL'}] 3'.res: 新 / 旧 最大 {n['qmax'] / o['qmax']:.3f}・RMS {n['qrms'] / o['qrms']:.3f} (旧の再実行 / 旧 {ob['qmax'] / o['qmax']:.3f}・{ob['qrms'] / o['qrms']:.3f})")
print(f"  [{'PASS' if ok_u else 'FAIL'}] 3'.u10: 新 / 旧 {n['umax'] / o['umax']:.3f} (旧の再実行 / 旧 {ob['umax'] / o['umax']:.3f})")
json.dump({"3p.res": ok_r, "3p.u10": ok_u, "stats": S}, open(C45 / "fg2d_verdict.json", "w"), ensure_ascii=False, indent=1)
