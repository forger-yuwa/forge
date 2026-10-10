"""plan architecture-float-state-double-geometry §6.8 (V3、事前登録) の集計。fg6.sh の後に AWS の case/45 で回す。
- E = ‖R32 − R64‖₂ / ‖R64‖₂ (量・領域ごと)。R64 は参照の 1 本目 (ref64a)。R32 は各腕の 1 本目。
- 再実行の差 N_arm = max(‖R32a − R32b‖, ‖R64a − R64b‖) / ‖R64‖ (腕・量・領域ごと)。
- 悪化しない: E_new ≤ 1.1·E_old + 10·N_new。改善: 熱・運動量・乱流で E_old > 10·N_old なら E_new ≤ 0.5·E_old。
- 記録: 最大の誤差の位置、粘性の面の流束の誤差 (GEOMAB の新腕 − REF) を節点に集めたものの大きさ (§6.5 の H3 の切り分け)。
欠損・非有限は判定不能。結果を標準出力と v3_verdict.json に書く。"""
import json, h5py, numpy as np
from pathlib import Path
C = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); NJ = 121; RT = 0.076807
RK = ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega")
TARGET = {"res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega"}   # 熱・運動量・乱流 (改善を求める量)
def load(r):
    with h5py.File(C / r / "res_1.h5", "r") as h:
        out = {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in RK}
    for k, v in out.items():
        if not np.all(np.isfinite(v)): raise SystemExit(f"判定不能: {r} の {k} に非有限")
    return out
R = {t: load(f"run_{n:04d}_v3_{t}") for n, t in ((440, "old32a"), (441, "old32b"), (442, "new32a"), (443, "new32b"), (444, "ref64a"), (445, "ref64b"))}
with h5py.File(C / "run_0444_v3_ref64a/res_0.h5", "r") as h:
    X = np.asarray(h["MESH/COORD"][:], np.float64).reshape(-1, 3) if "MESH/COORD" in h else None
n = R["ref64a"]["res_ro"].size; jj = np.arange(n) % NJ
REG = {"全域": np.ones(n, bool), "壁際 j≥117": jj >= 117}
V = {"rows": [], "fail": []}
print("量          領域       E_old      E_new      N_old      N_new      新/旧    判定")
for k in RK:
    for rg, m in REG.items():
        nr = np.linalg.norm(R["ref64a"][k][m])
        if nr == 0: print(f"{k:12s} {rg}: 参照が 0 で比較外"); continue
        Eo = np.linalg.norm(R["old32a"][k][m] - R["ref64a"][k][m]) / nr; En = np.linalg.norm(R["new32a"][k][m] - R["ref64a"][k][m]) / nr
        n64 = np.linalg.norm(R["ref64b"][k][m] - R["ref64a"][k][m]) / nr
        No = max(np.linalg.norm(R["old32b"][k][m] - R["old32a"][k][m]) / nr, n64); Nn = max(np.linalg.norm(R["new32b"][k][m] - R["new32a"][k][m]) / nr, n64)
        ok1 = En <= 1.1 * Eo + 10 * Nn
        need = (k in TARGET) and Eo > 10 * No
        ok2 = (En <= 0.5 * Eo) if need else True
        tag = ("PASS" if ok1 and ok2 else "FAIL") + (" (改善を要求)" if need else "")
        print(f"{k:12s} {rg:10s} {Eo:.3e}  {En:.3e}  {No:.3e}  {Nn:.3e}  {En / Eo if Eo > 0 else float('nan'):.3f}   {tag}")
        im = np.flatnonzero(m)[np.argmax(np.abs(R["new32a"][k][m] - R["ref64a"][k][m]))]
        io = np.flatnonzero(m)[np.argmax(np.abs(R["old32a"][k][m] - R["ref64a"][k][m]))]
        loc = lambda i: f"j {i % NJ}" + (f"・x/r_t {X[i, 0] / RT:.2f}" if X is not None else "")
        V["rows"].append(dict(k=k, region=rg, E_old=Eo, E_new=En, N_old=No, N_new=Nn, improve_required=need, verdict=tag, max_new_at=loc(im), max_old_at=loc(io)))
        if not (ok1 and ok2): V["fail"].append(f"{k} {rg}")
print("== 最大の誤差の位置 (新 / 旧)")
for r in V["rows"]: print(f"  {r['k']:12s} {r['region']:10s} 新 {r['max_new_at']}  旧 {r['max_old_at']}")
print("== 記録: 粘性の面の流束の誤差を節点に集めたもの (GEOMAB の新腕 − REF) と、残差全体の誤差 (新) の比較")
try:
    d = h5py.File(C / "run_0446_v3_dump32/geomab.h5", "r"); rf = h5py.File(C / "run_0446_v3_dump32/geomab_ref.h5", "r")
    pc = np.asarray(d["/mesh/plane_cells"][:]).reshape(-1, 2); nN = int(d.attrs["nNormalPlanes"])
    Fn = np.asarray(d["/viscous/new/face_flux"][:], np.float64); Fr = np.asarray(rf["/viscous/ref/face_flux"][:], np.float64)
    dF = np.where(np.isfinite(Fn) & np.isfinite(Fr), Fn - Fr, 0.0)
    for comp, k in ((0, "res_roUx"), (1, "res_roUy"), (3, "res_roe")):
        acc = np.zeros(n)
        a, b = pc[:nN, 0], pc[:nN, 1]
        ma, mb = a < n, b < n
        np.add.at(acc, a[ma], dF[comp, :nN][ma]); np.add.at(acc, b[mb], -dF[comp, :nN][mb])
        for rg, m in REG.items():
            nr = np.linalg.norm(R["ref64a"][k][m])
            print(f"  {k:9s} {rg:10s}: 粘性の面の流束の誤差を集めた大きさ / ‖R64‖ = {np.linalg.norm(acc[m]) / nr:.3e}、残差全体の誤差 E_new = {np.linalg.norm(R['new32a'][k][m] - R['ref64a'][k][m]) / nr:.3e}")
except Exception as e:
    print(f"  (集計できない: {e})")
json.dump(V, open(C / "v3_verdict.json", "w"), ensure_ascii=False, indent=1, default=float)
print("== 総合:", "PASS" if not V["fail"] else f"FAIL {V['fail']}")
