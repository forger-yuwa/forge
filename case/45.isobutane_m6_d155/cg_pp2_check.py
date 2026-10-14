"""plan tooling-nozzle-core-grid §4.84 P-m2 の資格の検査 (事前登録、codex 諮問 2026-10-13 point-precond2)。
(1) P2・C′ の restart_field の VERDICT がビット一致。設定の差: P2 と C は lowMachPrecond の値 (0 → 2) の 1 行だけ、C′ と C は nStepOuter の 1 行だけ。
(2) 設定の写しの lowMachPrecond が P2 で 2。precondEps・speciesImplicitCoupling・speciesPrecondDt が P2 の設定に無いこと (ソースの既定値 0.15・0・1 を記録)。
(3) 刻みの効き: 島 (x/r_t −6〜−3.25・壁から 60〜100 節点目) の step 2500 で、dt_local(P2)/dt_local(C′) を P2 の状態から
    setDTlocal_precond_scale_d の式 max((|u|+c)/(½(1+β)|u| + c′), 1) で求めた期待の倍率で割った値の中央値が 0.85〜1.15、
    かつ島の |c_P2/c_C′ − 1| の中央値 ≤ 0.01 (状態が近い)。
(4) 境界節点の凍結なし: 流入・流出・軸・壁 (角を除く) の ρ の step 0 → 2500 の変化の中央値 / 1 節点内側の変化の中央値 = r_g が、P2 で C′ の 0.1 倍以上。
    step 20000 の P2 の Y0 のばらつき ≤ 1e−9。
(5) C (流用) の同一性: C・P2・C′ の RUN_PROVENANCE の forge_sha256 が同じ (d58ed787…)、C の restart_field の VERDICT がビット一致、C の領域の残差の系列に縮流部の外の列がある。
usage: python3 cg_pp2_check.py <C′ の run> <P2 の run> [C の run (既定 run_0688_pp_c)] → _band_ab/core_grid/pp2_check.json"""
import json, re, subprocess, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
CP, P2 = sys.argv[1:3]; C = sys.argv[3] if len(sys.argv) > 3 else "run_0688_pp_c"; EPS = 0.15

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

def diff_added(a, b):
    out = subprocess.run(["diff", str(HERE / a / "solverConfig.yaml"), str(HERE / b / "solverConfig.yaml")], capture_output=True, text=True).stdout
    return [l for l in out.splitlines() if l.startswith(">")]

def prov_sha(run):
    m = re.search(r"forge_sha256\s*:\s*([0-9a-f]+)", (HERE / run / "RUN_PROVENANCE.txt").read_text())
    return m.group(1) if m else None

out = {"C": C, "Cprime": CP, "P2": P2}
# (1)
bit = {r: "VERDICT: OK (9 量を移した、SRC とビット一致)" in (HERE / r / "restart.log").read_text() for r in (C, CP, P2)}
dP, dC = diff_added(C, P2), diff_added(C, CP)
out["restart_bitwise"] = bit; out["config_diff_P2_vs_C"] = dP; out["config_diff_Cprime_vs_C"] = dC
q1 = all(bit.values()) and len(dP) == 1 and "lowMachPrecond: 2, " in dP[0] and len(dC) == 1 and "nStepOuter: 2500}" in dC[0]
# (2)
y = (HERE / P2 / "solverConfig.yaml").read_text()
absent = all(k not in y for k in ("precondEps", "speciesImplicitCoupling", "speciesPrecondDt"))
out["effective"] = {"lowMachPrecond_P2": cfg_val(P2, "lowMachPrecond"), "lowMachPrecond_Cprime": cfg_val(CP, "lowMachPrecond"),
                    "keys_absent_in_P2_yaml": absent, "source_defaults": {"precondEps": 0.15, "speciesImplicitCoupling": 0, "speciesPrecondDt": 1}}
q2 = out["effective"]["lowMachPrecond_P2"] == "2" and out["effective"]["lowMachPrecond_Cprime"] == "0" and absent
# (3)
gc, X, ni, nj = load(CP, 2500, ("dt_local", "sonic")); gp, _, _, _ = load(P2, 2500, ("dt_local", "Ux", "Uy", "sonic"))
isl = np.zeros((ni, nj), bool); cols = (X[:, 0] >= -6.0) & (X[:, 0] <= -3.25); isl[np.ix_(cols, np.arange(60, 101))] = True
c = np.maximum(gp["sonic"], 1e-8); vm = np.hypot(gp["Ux"], gp["Uy"]); Ur = np.minimum(c, np.maximum(vm, EPS * c)); beta = (Ur / c) ** 2
cp = 0.5 * np.sqrt((1 - beta) ** 2 * vm ** 2 + 4 * Ur ** 2); expect = np.maximum((vm + c) / np.maximum(0.5 * (1 + beta) * vm + cp, 1e-30), 1.0)
meas = gp["dt_local"] / gc["dt_local"]; q = meas[isl] / expect[isl]; close = float(np.median(np.abs(gp["sonic"][isl] / gc["sonic"][isl] - 1)))
out["dt"] = {"island_median_measured_ratio": float(np.median(meas[isl])), "island_median_expected_ratio": float(np.median(expect[isl])),
             "island_median_measured_over_expected": float(np.median(q)), "island_p10_p90_measured_over_expected": [float(np.percentile(q, 10)), float(np.percentile(q, 90))],
             "island_median_abs_sonic_rel_diff": close, "island_median_M": float(np.median(vm[isl] / c[isl]))}
q3 = 0.85 <= out["dt"]["island_median_measured_over_expected"] <= 1.15 and close <= 0.01
# (4)
def rgroups(run):
    a0, _, _, _ = load(run, 0, ("ro",)); a1, _, _, _ = load(run, 2500, ("ro",)); d = np.abs(a1["ro"] - a0["ro"])
    g = {"inlet": (d[0, 1:-1], d[1, 1:-1]), "outlet": (d[-1, 1:-1], d[-2, 1:-1]), "axis": (d[1:-1, -1], d[1:-1, -2]), "wall": (d[1:-1, 0], d[1:-1, 1])}
    return {k: (float(np.median(b) / np.median(n)) if np.median(n) > 0 else None) for k, (b, n) in g.items()}
rc, rp = rgroups(CP), rgroups(P2)
out["boundary_r"] = {"Cprime": rc, "P2": rp}
q4 = all(rc[k] is not None and rp[k] is not None and rp[k] >= 0.1 * rc[k] for k in rc)
try:
    gy, _, _, _ = load(P2, 20000, ("Y0",)); out["Y0_ptp_20000_P2"] = float(np.ptp(gy["Y0"]))
except (FileNotFoundError, OSError, KeyError):
    out["Y0_ptp_20000_P2"] = None
q4 = q4 and out["Y0_ptp_20000_P2"] is not None and out["Y0_ptp_20000_P2"] <= 1e-9
# (5)
shas = {r: prov_sha(r) for r in (C, CP, P2)}
rowsC = [json.loads(l) for l in open(OUTD / f"region_res_{C}.jsonl")]
out["C_identity"] = {"forge_sha256": shas, "C_has_outside_columns": all("outside_res_ro" in r for r in rowsC), "C_outputs": sorted(r["step"] for r in rowsC)}
q5 = len(set(shas.values())) == 1 and all(s and s.startswith("d58ed787") for s in shas.values()) and bit[C] and out["C_identity"]["C_has_outside_columns"]
out["qualification"] = {"same_start_and_config": bool(q1), "effective_values": bool(q2), "precond_dt_effective": bool(q3),
                        "boundary_not_frozen_and_composition": bool(q4), "control_identity": bool(q5)}
(OUTD / "pp2_check.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print(json.dumps({k: out[k] for k in ("config_diff_P2_vs_C", "config_diff_Cprime_vs_C", "effective", "dt", "boundary_r", "Y0_ptp_20000_P2", "C_identity", "qualification")}, ensure_ascii=False, indent=1))
