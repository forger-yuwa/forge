#!/usr/bin/env python3
"""Python 側 NASA-9 (tools/total_quantities._TPGas) の範囲外処理がソルバ (cuda_forge/thermo_d.cuh) と一致することの単体試験
(codex 2026-09-16 result-3 M2)。

参照実装はこのファイル内で thermo_d.cuh の式 (thermo_pick_coeffs: T<Tmid で low、thermo_cp_molar: Tlo/Thi でクランプ、
thermo_h_molar: 範囲外は h(Tb)+cp(Tb)(T−Tb)、thermo_s0_mass: s°(Tb)+cp(Tb) ln(T/Tb)、datum は a7 のシフト = h(Tref) の引き算)
から独立に書き、_TPGas は呼ばない。
  (1) 内蔵 N2/H2O 係数で T=100, 150, 199, 250, 2500, 6500 K の cp/h/s° を _TPGas と参照で比較 (rtol 1e-12)。
  (2) codex の反例: [N2,H2O]=[0.95,0.05] の 100 K / 150 K セルを reinit で [0.8,0.2] に変換 (roe += ρ[e_dst(T)−e_src(T)]) した
      後の roe を参照 EOS (Newton) で反転し、T が元の値と 1e-6 K 以内 (修正前は 101.38 K / 150.10 K だった)。
  (3) 変換器本体 (tools/convert_species_field.py --mode reinit) を 1 セル入力で実際に回し、書かれた roe を参照 EOS で反転して同じ判定。
usage: python3 tests/unit/test_tpgas_lowT.py     ([PASS]/[FAIL], 失敗があれば非ゼロ終了)
"""
import os, shutil, subprocess, sys, tempfile
import numpy as np, h5py, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, TOOLS)
from total_quantities import _TPGas  # noqa: E402

RU = 8.314462618
DB = {
    "N2": {"MW": 0.0280134, "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
           "nasa9_low": [2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03, 1.384646189e-05,
                         -9.625793620e-09, 2.519705809e-12, 7.108460860e+02, -1.076003744e+01],
           "nasa9_high": [5.877124060e+05, -2.239249073e+03, 6.066949220e+00, -6.139685500e-04, 1.491806679e-07,
                          -1.923105485e-11, 1.061954386e-15, 1.283210415e+04, -1.586640027e+01]},
    "H2O": {"MW": 0.0180153, "LJ_sigma": 2.605, "LJ_eps_kB": 572.4, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
            "nasa9_low": [-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03, -7.342557370e-06,
                          4.955043490e-09, -1.336933246e-12, -3.303974310e+04, 1.724205775e+01],
            "nasa9_high": [1.034972096e+06, -2.412698562e+03, 4.646110780e+00, 2.291998307e-03, -6.836830480e-07,
                           9.426468930e-11, -4.822380530e-15, -1.384286509e+04, -7.978148510e+00]},
}
NAMES = ["N2", "H2O"]
TREF = 298.15
g_fail = 0


def check(ok, what):
    global g_fail
    print(("[PASS] " if ok else "[FAIL] ") + what)
    if not ok:
        g_fail += 1


# ----------------------------------------------------------------------------- 参照実装 (thermo_d.cuh の式をそのまま; スカラ double)
def _pick(sp, Tc):
    return sp["nasa9_low"] if Tc < sp["Tmid"] else sp["nasa9_high"]


def _cp_molar_clamped(sp, Tc):
    a = _pick(sp, Tc); Ti = 1.0 / Tc
    return RU * (a[0] * Ti * Ti + a[1] * Ti + a[2] + a[3] * Tc + a[4] * Tc ** 2 + a[5] * Tc ** 3 + a[6] * Tc ** 4)


def _h_molar_clamped(sp, Tc):
    a = _pick(sp, Tc); Ti = 1.0 / Tc
    hRT = (-a[0] * Ti * Ti + a[1] * np.log(Tc) * Ti + a[2] + a[3] * Tc / 2.0 + a[4] * Tc ** 2 / 3.0
           + a[5] * Tc ** 3 / 4.0 + a[6] * Tc ** 4 / 5.0 + a[7] * Ti)
    return RU * Tc * hRT


def _s_molar_clamped(sp, Tc):
    a = _pick(sp, Tc); Ti = 1.0 / Tc
    return RU * (-a[0] * Ti * Ti / 2.0 - a[1] * Ti + a[2] * np.log(Tc) + a[3] * Tc + a[4] * Tc ** 2 / 2.0
                 + a[5] * Tc ** 3 / 3.0 + a[6] * Tc ** 4 / 4.0 + a[8])


def ref_cp_mass(sp, T):
    Tc = min(max(T, sp["Tlo"]), sp["Thi"])
    return _cp_molar_clamped(sp, Tc) / sp["MW"]


def ref_h_mass(sp, T, Tref=TREF):
    if T < sp["Tlo"]:
        h = _h_molar_clamped(sp, sp["Tlo"]) + _cp_molar_clamped(sp, sp["Tlo"]) * (T - sp["Tlo"])
    elif T > sp["Thi"]:
        h = _h_molar_clamped(sp, sp["Thi"]) + _cp_molar_clamped(sp, sp["Thi"]) * (T - sp["Thi"])
    else:
        h = _h_molar_clamped(sp, T)
    # sensible datum: thermo_init_db は a7 を −h(Tref)/Ru だけシフト = 全 T で h(Tref) を引く
    href = _h_molar_clamped(sp, Tref) if Tref > 0 else 0.0
    return (h - href) / sp["MW"]


def ref_s0_mass(sp, T):
    if T < sp["Tlo"]:
        s = _s_molar_clamped(sp, sp["Tlo"]) + _cp_molar_clamped(sp, sp["Tlo"]) * np.log(T / sp["Tlo"])
    elif T > sp["Thi"]:
        s = _s_molar_clamped(sp, sp["Thi"]) + _cp_molar_clamped(sp, sp["Thi"]) * np.log(T / sp["Thi"])
    else:
        s = _s_molar_clamped(sp, T)
    return s / sp["MW"]


def ref_mix(Y, T):
    sps = [DB[n] for n in NAMES]
    R = sum(y * RU / sp["MW"] for y, sp in zip(Y, sps))
    cp = sum(y * ref_cp_mass(sp, T) for y, sp in zip(Y, sps))
    h = sum(y * ref_h_mass(sp, T) for y, sp in zip(Y, sps))
    s0 = sum(y * ref_s0_mass(sp, T) for y, sp in zip(Y, sps))
    return R, cp, h, s0


def ref_T_from_e(Y, e, T0=300.0):
    """thermo_T_from_e 相当 (Newton, 収束 1e-12)。"""
    T = T0
    for _ in range(200):
        R, cp, h, _ = ref_mix(Y, T)
        f = h - R * T - e; df = cp - R
        dT = f / df
        T -= dT
        if abs(dT) < 1e-12 * T:
            break
    return T


def main():
    gas = _TPGas(DB, NAMES, TREF)
    # (1) 単体物性
    Ts = [100.0, 150.0, 199.0, 250.0, 2500.0, 6500.0]
    Y = [0.95, 0.05]
    worst = 0.0; worst_what = ""
    ok = True
    for T in Ts:
        Ta = np.array([T]); Yl = [np.array([Y[0]]), np.array([Y[1]])]
        R, cp, h, s0 = ref_mix(Y, T)
        # 判定は rel 1e-12 + abs (h: 1e-8 J/kg, cp/s0/R: 1e-10 J/kg/K); datum 近傍で h≈0 になる点は相対差が意味を持たない
        for name, mine, ref, atol in (("cp", gas.cp(Yl, Ta)[0], cp, 1e-10), ("h", gas.h(Yl, Ta)[0], h, 1e-8),
                                      ("s0", gas.s0(Yl, Ta)[0], s0, 1e-10), ("R", gas.Rmix(Yl), R, 1e-10)):
            d = abs(mine - ref)
            ok = ok and (d <= 1e-12 * abs(ref) + atol)
            rel = d / max(abs(ref), 1e-300)
            if rel > worst:
                worst, worst_what = rel, f"{name}@{T:.0f}K (mine {mine!r}, ref {ref!r}, |Δ| {d:.2e})"
    check(ok, f"_TPGas cp/h/s0/R vs reference (thermo_d.cuh formulas) at T={Ts}: worst rel diff {worst:.2e} at {worst_what}")

    # (2) codex counterexample in memory: reinit [0.95,0.05] -> [0.8,0.2] at 100 K / 150 K
    for T in (100.0, 150.0):
        Ys = [np.array([0.95]), np.array([0.05])]; Yd = [np.array([0.8]), np.array([0.2])]; Ta = np.array([T])
        e_src = ref_mix([0.95, 0.05], T)[2] - ref_mix([0.95, 0.05], T)[0] * T         # ソルバ側の保存量 (参照)
        e_new = e_src + (gas.h(Yd, Ta) - gas.Rmix(Yd) * Ta)[0] - (gas.h(Ys, Ta) - gas.Rmix(Ys) * Ta)[0]   # 変換器の差分形
        Tchk = ref_T_from_e([0.8, 0.2], e_new, T0=200.0)
        check(abs(Tchk - T) < 1e-6, f"reinit [0.95,0.05]->[0.8,0.2] at {T:.0f} K: solver-formula inversion of converted e -> {Tchk:.9f} K (|ΔT| {abs(Tchk - T):.2e})")

    # (3) the real tool on a 1-cell input
    root = tempfile.mkdtemp(prefix="test_tpgas_lowT_")
    try:
        for T in (100.0, 150.0):
            for tag, Yt in (("src", [0.95, 0.05]), ("dst", [0.8, 0.2])):
                d = os.path.join(root, f"{tag}_{int(T)}"); os.makedirs(d)
                open(os.path.join(d, "solverConfig.yaml"), "w").write(
                    "mesh: {meshFormat: \"hdf5\", discretization: \"node\", meshFileName: \"in.h5\", valueFileName: \"in.h5\"}\ngpu: 1\nsolver: \"SLAU\"\n"
                    "physProp: {isCompressible: 1, thermalMethod: 2, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1000.0, gamma: 1.4,\n"
                    "           species: [N2, H2O], speciesDBFile: \"species_db.yaml\", thermoHrefTemp: 298.15" + (", tracer: exhaust" if tag == "src" else "") + "}\n"
                    "time: {unsteady: 0, dualTime: 0, last: {control: 0, nStepOuter: 1}, deltaT: {control: 1, dt: 1e-8, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-9, dt_max: 1e-3},"
                    " outStepStart: 0, outStepInterval: 1, timeIntegration: 11, nStepInner: 5}\n"
                    "space: {convMethod: 1, limiter: 2}\nturbulence: {model: \"none\"}\ninitial: \"uniform_p101325_u10\"\noutput: {level: 1}\n")
                yaml.safe_dump(DB, open(os.path.join(d, "species_db.yaml"), "w"), sort_keys=False)
                meta = {"mode": "full", "species": NAMES, "keep": [], "condensing_species": None, "condensing_index": None,
                        "tracer": {"enabled": tag == "src", "name": "roXi" if tag == "src" else None, "definition": None}, "lumps": {},
                        "expansion": {s: {s: 1.0} for s in NAMES},
                        "streams": {"inflow": {"X": {}, "Y": dict(zip(NAMES, Yt)), "Y_transport": list(Yt), "sum_input": 1.0},
                                    "external": {"X": {}, "Y": {"N2": 1.0}, "Y_transport": [1.0, 0.0], "sum_input": 1.0}},
                        "MW": {s: DB[s]["MW"] for s in NAMES}}
                yaml.safe_dump(meta, open(os.path.join(d, "species_meta.yaml"), "w"), sort_keys=False)
                ro = 0.8; R, cp, h, _ = ref_mix(Yt, T); e = h - R * T
                with h5py.File(os.path.join(d, "in.h5"), "w") as f:
                    f.create_dataset("MESH/COORD", data=np.zeros(3))
                    for k, v in {"ro": ro, "roUx": 0.0, "roUy": 0.0, "roUz": 0.0, "roe": ro * e, "roY0": ro * Yt[0], "roY1": ro * Yt[1],
                                 "roXi": ro * 1.0}.items():
                        if k == "roXi" and tag != "src":
                            continue
                        f.create_dataset("VALUE/" + k, data=np.array([v], np.float64))
            src = os.path.join(root, f"src_{int(T)}"); dst = os.path.join(root, f"dst_{int(T)}")
            p = subprocess.run([sys.executable, os.path.join(TOOLS, "convert_species_field.py"), os.path.join(src, "in.h5"), os.path.join(dst, "in.h5"),
                                "--meta", os.path.join(dst, "species_meta.yaml"), "--src-meta", os.path.join(src, "species_meta.yaml"), "--mode", "reinit"],
                               capture_output=True, text=True)
            ok_rc = p.returncode == 0
            with h5py.File(os.path.join(dst, "in.h5")) as f:
                ro = float(f["VALUE/ro"][0]); roe = float(f["VALUE/roe"][0]); Yw = [float(f["VALUE/roY0"][0]) / ro, float(f["VALUE/roY1"][0]) / ro]
            Tchk = ref_T_from_e(Yw, roe / ro, T0=200.0)
            check(ok_rc and abs(Tchk - T) < 1e-6 and abs(Yw[1] - 0.2) < 1e-12,
                  f"tool --mode reinit 1-cell {T:.0f} K -> written Y {np.round(Yw, 6).tolist()}, solver-formula T {Tchk:.9f} K (|ΔT| {abs(Tchk - T):.2e}, rc {p.returncode})")
            if not ok_rc:
                print("      " + "\n      ".join((p.stdout + p.stderr).strip().splitlines()[-4:]))
    finally:
        shutil.rmtree(root, ignore_errors=True)
    print(("ALL PASS" if g_fail == 0 else "FAILED") + f" ({g_fail} failures)")
    return 0 if g_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
