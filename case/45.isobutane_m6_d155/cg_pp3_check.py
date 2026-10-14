"""plan tooling-nozzle-core-grid §4.82 P-m3 の資格の検査 (事前登録、codex 諮問 2026-10-13 point-precond)。
(1) 前処理が実際に効いている: 島 (x/r_t −6〜−3.25・壁から 60〜100 節点目) の step 2500 の dt_local の比 P/C を、P の状態から setDT_d.cu の
    setDTlocal_precond_scale_d の式 max((|u|+c)/(½(1+β)|u| + c′), 1) (β = (U_r/c)²、U_r = min(c, max(|u|, εc))、c′ = ½√((1−β)²|u|² + 4U_r²)、ε = 0.15)
    で求めた期待の倍率で割った値の中央値が 0.85〜1.15。設定の写しが P で lowMachPrecond 3、C で 0。
(2) 両腕の restart_field の VERDICT がビット一致、step 0 (inner −1) の全 rms_ 列が両腕で相対 1e−6 以内に一致。
(3) 境界節点が凍結していない: 流入 (i = 0)・流出 (i = ni−1)・軸 (jw = nj−1)・壁 (jw = 0) の各群 (角を除く) で、
    ρ の step 0 → 2500 の変化の中央値 / 1 節点内側の変化の中央値 = r_g が、P で C の 0.1 倍以上。
(4) の一部: step 20000 の Y0 のばらつき ≤ 1e−9 (両腕)。speciesImplicitCoupling・speciesPrecondDt の実効値を記録する。
usage: python3 cg_pp3_check.py <C の run> <P の run> → _band_ab/core_grid/pp3_check.json"""
import csv, json, re, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
C, P = sys.argv[1:3]; EPS = 0.15

def load(run, n, keys):
    info = json.loads((HERE / run / "prepare_info.json").read_text()); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        v = {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
        x = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, 0] / S
    nj = x.size // ni
    return {k: a.reshape(ni, nj)[:, ::-1] for k, a in v.items()}, x.reshape(ni, nj)[:, ::-1], ni, nj   # [i, jw]

def cfg_val(run, key):
    m = re.search(r"'%s' in '[^']*': ([^\s]+)" % key, (HERE / run / "forge_run.log").read_text(errors="replace"))
    return m.group(1) if m else None

out = {"C": C, "P": P, "config": {r: {k: cfg_val(r, k) for k in ("lowMachPrecond", "precondEps", "speciesImplicitCoupling", "speciesPrecondDt", "lineImplicit", "cfl_pseudo")} for r in (C, P)}}
# (1)
gc, X, ni, nj = load(C, 2500, ("dt_local",)); gp, _, _, _ = load(P, 2500, ("dt_local", "Ux", "Uy", "sonic"))
isl = np.zeros((ni, nj), bool); cols = (X[:, 0] >= -6.0) & (X[:, 0] <= -3.25); isl[np.ix_(cols, np.arange(60, 101))] = True
c = np.maximum(gp["sonic"], 1e-8); vm = np.hypot(gp["Ux"], gp["Uy"]); Ur = np.minimum(c, np.maximum(vm, EPS * c)); beta = (Ur / c) ** 2
cp = 0.5 * np.sqrt((1 - beta) ** 2 * vm ** 2 + 4 * Ur ** 2); expect = np.maximum((vm + c) / np.maximum(0.5 * (1 + beta) * vm + cp, 1e-30), 1.0)
meas = gp["dt_local"] / gc["dt_local"]; q = meas[isl] / expect[isl]
out["dt"] = {"island_median_measured_ratio": float(np.median(meas[isl])), "island_median_expected_ratio": float(np.median(expect[isl])),
             "island_median_measured_over_expected": float(np.median(q)), "island_p10_p90_measured_over_expected": [float(np.percentile(q, 10)), float(np.percentile(q, 90))],
             "island_median_M": float(np.median(vm[isl] / c[isl]))}
q1 = out["config"][P]["lowMachPrecond"] == "3" and out["config"][C]["lowMachPrecond"] == "0" and 0.85 <= out["dt"]["island_median_measured_over_expected"] <= 1.15
# (2)
def step0(run):
    rows = list(csv.reader(open(HERE / run / "residual_history.csv"))); h = rows[0]
    r0 = next(r for r in rows[1:] if r[0] == "0" and r[1] == "-1")
    return {h[i]: float(r0[i]) for i in range(len(h)) if h[i].startswith("rms_")}
s0c, s0p = step0(C), step0(P)
rel = {k: (abs(s0p[k] - s0c[k]) / max(abs(s0c[k]), abs(s0p[k])) if max(abs(s0c[k]), abs(s0p[k])) > 0 else 0.0) for k in s0c}
bit = all("VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / r / "restart.log").read_text() for r in (C, P))
out["step0_rel_diff"] = rel; out["restart_bitwise"] = bit
q2 = bit and all(v <= 1e-6 for v in rel.values())
# (3)
def rgroups(run):
    a0, _, _, _ = load(run, 0, ("ro",)); a1, _, _, _ = load(run, 2500, ("ro",)); d = np.abs(a1["ro"] - a0["ro"])
    g = {"inlet": (d[0, 1:-1], d[1, 1:-1]), "outlet": (d[-1, 1:-1], d[-2, 1:-1]), "axis": (d[1:-1, -1], d[1:-1, -2]), "wall": (d[1:-1, 0], d[1:-1, 1])}
    return {k: (float(np.median(b) / np.median(n)) if np.median(n) > 0 else None) for k, (b, n) in g.items()}
rc, rp = rgroups(C), rgroups(P)
out["boundary_r"] = {"C": rc, "P": rp}
q3 = all(rc[k] is not None and rp[k] is not None and rp[k] >= 0.1 * rc[k] for k in rc)
# (4) の一部
ycomp = {}
for r in (C, P):
    try:
        gy, _, _, _ = load(r, 20000, ("Y0",)); ycomp[r] = float(np.ptp(gy["Y0"]))
    except (FileNotFoundError, OSError, KeyError):
        ycomp[r] = None
out["Y0_ptp_20000"] = ycomp
q4y = all(v is not None and v <= 1e-9 for v in ycomp.values())
out["qualification"] = {"precond_effective": bool(q1), "same_start_and_rhs": bool(q2), "boundary_not_frozen": bool(q3), "composition_uniform": bool(q4y)}
(OUTD / "pp3_check.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print(json.dumps({k: out[k] for k in ("config", "dt", "boundary_r", "Y0_ptp_20000", "qualification")}, ensure_ascii=False, indent=1))
print("step 0 の相対差の最大:", max(rel.values()))
