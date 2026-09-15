#!/usr/bin/env python3
"""tools/convert_species_field.py の失敗系試験 (codex 2026-09-16 result-2 M3)。

合成の 2 種 (N2, H2O; 内蔵 DB と同じ NASA-9 係数) の run を一時ディレクトリに作り、変換器を子プロセスで呼んで
  - 正常入力は通る (exit 0, "all checks passed")
  - roe=NaN セル / 負の ρ / ΣY≠1 / 負の Y / T が括弧端 (異常に低い e) / roY データセット欠落 /
    species_meta と solverConfig の矛盾 (species 順序, tracer.enabled; source 側・destination 側) / reinit で destination の
    Y_transport の和≠1 / roK=NaN / roOmega=Inf / 負の roK / roOmega=1e39 (float32 への変換で overflow) / tracer 必須なのに ξ が導けない
  がそれぞれ**書き込み前に**明確なメッセージで拒否され非ゼロ終了することを確認する。
usage: python3 tests/unit/test_convert_species_field_fail.py     ([PASS]/[FAIL], 失敗があれば非ゼロ終了)
"""
import os, shutil, subprocess, sys, tempfile
import numpy as np, h5py, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, TOOLS)
from total_quantities import _TPGas  # noqa: E402

TOOL = os.path.join(TOOLS, "convert_species_field.py")
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
N = 8
g_fail = 0


def check(ok, what):
    global g_fail
    print(("[PASS] " if ok else "[FAIL] ") + what)
    if not ok:
        g_fail += 1


def make_run(d, tracer=False, meta_tracer=None):
    os.makedirs(d, exist_ok=True)
    cfg = ("mesh: {meshFormat: \"hdf5\", discretization: \"node\", meshFileName: \"in.h5\", valueFileName: \"in.h5\"}\ngpu: 1\nsolver: \"SLAU\"\n"
           "physProp: {isCompressible: 1, thermalMethod: 2, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1000.0, gamma: 1.4,\n"
           "           species: [N2, H2O], speciesDBFile: \"species_db.yaml\", thermoHrefTemp: 298.15" + (", tracer: exhaust" if tracer else "") + "}\n"
           "time: {unsteady: 0, dualTime: 0, last: {control: 0, nStepOuter: 1}, deltaT: {control: 1, dt: 1e-8, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-9, dt_max: 1e-3},"
           " outStepStart: 0, outStepInterval: 1, timeIntegration: 11, nStepInner: 5}\n"
           "space: {convMethod: 1, limiter: 2}\nturbulence: {model: \"none\"}\ninitial: \"uniform_p101325_u10\"\noutput: {level: 1}\n")
    open(os.path.join(d, "solverConfig.yaml"), "w").write(cfg)
    yaml.safe_dump(DB, open(os.path.join(d, "species_db.yaml"), "w"), sort_keys=False)
    meta = {"mode": "full", "species": NAMES, "keep": [], "condensing_species": None, "condensing_index": None,
            "tracer": {"enabled": (tracer if meta_tracer is None else meta_tracer), "name": "roXi" if tracer else None, "definition": None},
            "lumps": {}, "expansion": {s: {s: 1.0} for s in NAMES},
            "streams": {"inflow": {"X": {}, "Y": {"N2": 0.9, "H2O": 0.1}, "Y_transport": [0.9, 0.1], "sum_input": 1.0},
                        "external": {"X": {}, "Y": {"N2": 1.0}, "Y_transport": [1.0, 0.0], "sum_input": 1.0}},
            "MW": {s: DB[s]["MW"] for s in NAMES}}
    yaml.safe_dump(meta, open(os.path.join(d, "species_meta.yaml"), "w"), sort_keys=False)
    return meta


def write_h5(path, ro, Y, T, u=100.0, roXi=None, drop=None, roK=None, roOmega=None):
    gas = _TPGas(DB, NAMES, 298.15)
    Yl = [Y[0], Y[1]]
    e = gas.h(Yl, T) - gas.Rmix(Yl) * T
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.zeros(3 * N))
        V = {"ro": ro, "roUx": ro * u, "roUy": np.zeros(N), "roUz": np.zeros(N), "roe": ro * (e + 0.5 * u * u),
             "roY0": ro * Y[0], "roY1": ro * Y[1]}
        if roXi is not None:
            V["roXi"] = roXi
        if roK is not None:
            V["roK"] = roK
        if roOmega is not None:
            V["roOmega"] = roOmega
        for k, v in V.items():
            if k != drop:
                # 乱流量だけ float64 で書く (float64 では有限だが float32 で overflow する値の試験用; 宛先 in.h5 は float32)
                f.create_dataset("VALUE/" + k, data=np.asarray(v, np.float64 if k in ("roK", "roOmega") else np.float32))


def run_tool(src, dst, extra=()):
    cmd = [sys.executable, TOOL, os.path.join(src, "in.h5"), os.path.join(dst, "in.h5"),
           "--meta", os.path.join(dst, "species_meta.yaml"), "--src-meta", os.path.join(src, "species_meta.yaml"), "--dry-run", *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def main():
    root = tempfile.mkdtemp(prefix="test_convert_species_field_")
    ro0 = np.full(N, 0.8); Y0 = np.array([np.full(N, 0.95), np.full(N, 0.05)]); T0 = np.linspace(300.0, 1500.0, N)
    dst = os.path.join(root, "dst"); make_run(dst); write_h5(os.path.join(dst, "in.h5"), ro0, Y0, T0)

    def case(name, expect_ok, keyword, ro=ro0, Y=Y0, T=T0, drop=None, roXi=None, dst_dir=dst, extra=(), src_tracer=False, patch=None,
             roK=None, roOmega=None, meta_tracer=None):
        src = os.path.join(root, "src_" + name); make_run(src, tracer=src_tracer, meta_tracer=meta_tracer)
        write_h5(os.path.join(src, "in.h5"), ro, Y, T, roXi=roXi, drop=drop, roK=roK, roOmega=roOmega)
        if patch:
            patch(src)
        rc, out = run_tool(src, dst_dir, extra)
        ok = (rc == 0) == expect_ok and (keyword in out)
        check(ok, f"{name}: expected {'OK' if expect_ok else 'REFUSE'} with '{keyword}' -> rc={rc}" + ("" if keyword in out else " (keyword NOT found)"))
        if not ok:
            print("      " + "\n      ".join(out.strip().splitlines()[-6:]))

    case("good", True, "all checks passed")
    case("good-reinit", True, "composition re-initialized", extra=("--mode", "reinit"), src_tracer=True, roXi=ro0 * 0.3)
    # roe = NaN in one cell
    def nan_roe(src):
        with h5py.File(os.path.join(src, "in.h5"), "r+") as f: a = f["VALUE/roe"][:]; a[3] = np.nan; f["VALUE/roe"][...] = a
    case("roe-nan", False, "非有限", patch=nan_roe)
    case("negative-rho", False, "ρ>0", ro=np.where(np.arange(N) == 2, -0.5, ro0))
    case("sumY-1.5", False, "|ΣY−1|", Y=np.array([np.full(N, 0.95), np.full(N, 0.55)]))
    case("negative-Y", False, "負値", Y=np.array([np.full(N, 1.01), np.full(N, -0.01)]))
    # T at the bracket end: absurdly low internal energy (below e(T_min)) -> Newton pins T at 50 K
    def low_e(src):
        with h5py.File(os.path.join(src, "in.h5"), "r+") as f: a = f["VALUE/roe"][:]; a[5] = -1.0e7 * 0.8; f["VALUE/roe"][...] = a
    case("T-bracket-end", False, "括弧端", patch=low_e)
    case("missing-roY1", False, "必須データセット", drop="roY1")
    # species_meta contradicts solverConfig (order swapped) in the source
    def swap_meta(src):
        m = yaml.safe_load(open(os.path.join(src, "species_meta.yaml"))); m["species"] = ["H2O", "N2"]
        yaml.safe_dump(m, open(os.path.join(src, "species_meta.yaml"), "w"))
    case("meta-config-contradiction", False, "矛盾" , patch=swap_meta)
    # destination meta with Y_transport not summing to 1 (reinit)
    dst_bad = os.path.join(root, "dst_bad"); make_run(dst_bad); write_h5(os.path.join(dst_bad, "in.h5"), ro0, Y0, T0)
    m = yaml.safe_load(open(os.path.join(dst_bad, "species_meta.yaml"))); m["streams"]["inflow"]["Y_transport"] = [0.9, 0.2]
    yaml.safe_dump(m, open(os.path.join(dst_bad, "species_meta.yaml"), "w"))
    case("reinit-Yt-sum", False, "和≠1", dst_dir=dst_bad, extra=("--mode", "reinit"))
    # destination requires a tracer but the source has neither roXi nor a stream label
    dst_tr = os.path.join(root, "dst_tracer"); make_run(dst_tr, tracer=True); write_h5(os.path.join(dst_tr, "in.h5"), ro0, Y0, T0)
    case("tracer-underivable", False, "roXi", dst_dir=dst_tr)
    # source tracer declared but roXi dataset missing
    case("missing-roXi", False, "必須データセット", src_tracer=True)
    # turbulence conserved variables: NaN / Inf must be refused (result-3 M3); finite roK/roOmega pass through
    case("roK-nan", False, "非有限", roK=np.where(np.arange(N) == 1, np.nan, 0.5), roOmega=np.full(N, 100.0))
    case("roOmega-inf", False, "非有限", roK=np.full(N, 0.5), roOmega=np.where(np.arange(N) == 6, np.inf, 100.0))
    case("roK-negative", False, "負値", roK=np.where(np.arange(N) == 4, -1.0, 0.5), roOmega=np.full(N, 100.0))
    case("good-turb", True, "all checks passed", roK=np.full(N, 0.5), roOmega=np.full(N, 100.0))
    # cast overflow: finite float64 1e39 becomes Inf in the float32 destination -> refused (result-4 M1); 1e30 stays finite -> pass
    case("roOmega-cast-overflow-1e39", False, "overflow", roK=np.full(N, 0.5), roOmega=np.full(N, 1.0e39))
    case("roOmega-near-limit-1e30", True, "all checks passed", roK=np.full(N, 0.5), roOmega=np.full(N, 1.0e30))
    # species_meta tracer.enabled contradicts solverConfig physProp.tracer (source side, and destination side) -> refused (result-3 M1)
    case("meta-tracer-contradiction-src", False, "矛盾", src_tracer=True, roXi=ro0 * 0.3, meta_tracer=False)
    dst_ct = os.path.join(root, "dst_tracer_contra"); make_run(dst_ct, tracer=True, meta_tracer=False); write_h5(os.path.join(dst_ct, "in.h5"), ro0, Y0, T0)
    case("meta-tracer-contradiction-dst", False, "矛盾", dst_dir=dst_ct, src_tracer=True, roXi=ro0 * 0.3)
    # destination config has tracer: source roXi is carried (same layout), never dropped
    dst_tr2 = os.path.join(root, "dst_tracer_ok"); make_run(dst_tr2, tracer=True); write_h5(os.path.join(dst_tr2, "in.h5"), ro0, Y0, T0)
    case("tracer-carried", True, "carried from source roXi", dst_dir=dst_tr2, src_tracer=True, roXi=ro0 * 0.3)

    print(("ALL PASS" if g_fail == 0 else "FAILED") + f" ({g_fail} failures); tmp {root}")
    shutil.rmtree(root, ignore_errors=True)
    return 0 if g_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
