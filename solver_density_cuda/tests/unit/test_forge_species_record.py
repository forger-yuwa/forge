#!/usr/bin/env python3
"""forge_species.py の内容照合と解決済み記録の検証 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #3a, §6 V1)。

  python3 solver_density_cuda/tests/unit/test_forge_species_record.py [--forge BIN]

(1) compare_signatures: 両側 builtin でも係数差を検出する (旧実装は source=builtin 同士を省略して [] を返した)。
    source だけの違い (同一係数を外部 DB で与える) は一致。係数不明 (記録の無い内蔵種) は「照合不能」。
(2) ソルバの記録 (forge --resolve-species) を Python で再計算: 互換性ハッシュが C++ と一致 (正規化の一字一句一致)、
    完全性ハッシュ = 全文 SHA-256。
(3) 記録の取り違え (別 run の記録を置く) と記録の改竄を load_record が検出する。
(2)(3) は forge バイナリが要る (--forge か環境変数 FORGE_BIN; 無ければ SKIP と表示し失敗扱いにしない)。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, copy, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import forge_species as fs  # noqa: E402

FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what)
    if not ok:
        FAIL += 1


# 内蔵 N2 (speciesDB.cpp) と同じ値の署名要素
N2_LOW = [2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03, 1.384646189e-05, -9.625793620e-09,
          2.519705809e-12, 7.108460860e+02, -1.076003744e+01]
N2_HIGH = [5.877124060e+05, -2.239249073e+03, 6.066949220e+00, -6.139685500e-04, 1.491806679e-07, -1.923105485e-11,
           1.061954386e-15, 1.283210415e+04, -1.586640027e+01]


def sig(source, low2=None):
    lo = list(N2_LOW)
    if low2 is not None:
        lo[2] = low2
    e = {"MW": 0.0280134, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0, "LJ_sigma": 3.621, "LJ_eps_kB": 97.53,
         "nasa9_low": lo, "nasa9_high": list(N2_HIGH), "source": source}
    return {"thermalMethod": 2, "names": ["N2"], "MW": [0.0280134], "thermoHrefTemp": 298.15, "tracer": None, "species": {"N2": e}}


def test_compare():
    a, b = sig("builtin"), sig("builtin", N2_LOW[2] + 1.0)   # codex の再現: 内蔵 N2.nasa9_low[2] を変える
    bad = fs.compare_signatures(a, b)
    check(any("N2.nasa9_low[2]" in x for x in bad), f"builtin vs builtin coefficient difference detected: {bad}")
    b = sig("builtin", N2_LOW[2] + 0.001)
    check(any("N2.nasa9_low[2]" in x for x in fs.compare_signatures(a, b)), "builtin vs builtin +0.001 detected")
    check(fs.compare_signatures(sig("builtin"), sig("file")) == [], "source-only difference (same coefficients) -> match")
    c = sig("builtin")
    c["species"]["N2"] = dict(c["species"]["N2"], nasa9_low=None, nasa9_high=None)
    bad = fs.compare_signatures(a, c)
    check(len(bad) == 1 and "unverifiable" in bad[0], f"unknown coefficients -> unverifiable: {bad}")
    d = copy.deepcopy(a); d["species"]["N2"]["LJ_sigma"] = 3.7
    check(any("LJ_sigma" in x for x in fs.compare_signatures(a, d)), "LJ difference detected")
    e = copy.deepcopy(a); e["thermoHrefTemp"] = 0.0
    check(any("thermoHrefTemp" in x for x in fs.compare_signatures(a, e)), "thermoHrefTemp difference detected")


def resolve(forge, run_dir):
    p = subprocess.run([forge, "--resolve-species"], cwd=run_dir, capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    return p.returncode, (lines[-1] if lines else ""), p.stderr


def write_case(d, db_low2=None):
    os.makedirs(d, exist_ok=True)
    db = ""
    if db_low2 is not None:
        lo = list(N2_LOW); lo[2] = db_low2
        with open(os.path.join(d, "species_db.yaml"), "w") as f:
            f.write("N2:\n  MW: 0.0280134\n  LJ_sigma: 3.621\n  LJ_eps_kB: 97.53\n  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n")
            f.write("  nasa9_low: [" + ", ".join("%.17g" % x for x in lo) + "]\n")
            f.write("  nasa9_high: [" + ", ".join("%.17g" % x for x in N2_HIGH) + "]\n")
        db = ', speciesDBFile: "species_db.yaml"'
    # 構成は case/44 run_0509 の solverConfig.yaml (実在の有効な config) から種と DB だけを差し替えたもの
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write('mesh: {discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: "m.h5", valueFileName: "m.h5"}\n'
                'gpu: 1\nsolver: "SLAU"\n'
                f'physProp: {{thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1220.7, gamma: 1.31526, '
                f'species: ["N2", "H2O"]{db}, thermoHrefTemp: 298.15}}\n'
                'time:\n  unsteady: 0\n  dualTime: 0\n  last: {nStepOuter: 1}\n'
                '  deltaT: {control: 1, dt: 1e-8, cfl: 6.0, cfl_pseudo: 6.0, dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, '
                'implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}\n'
                '  outStepStart: 0\n  outStepInterval: 1\n  timeIntegration: 11\n  nStepInner: 4\n'
                'space: {convMethod: 1, limiter: 2}\nturbulence: {model: "none"}\ninitial: "uniform_p101325_u10"\n')


def test_records(forge):
    tmp = tempfile.mkdtemp(prefix="forge_species_record_")
    try:
        A, F, M = (os.path.join(tmp, x) for x in ("A", "F", "M"))
        write_case(A); write_case(F, db_low2=N2_LOW[2]); write_case(M, db_low2=N2_LOW[2] + 0.001)
        res = {k: resolve(forge, d) for k, d in (("A", A), ("F", F), ("M", M))}
        for k, (rc, h, err) in res.items():
            check(rc == 0 and len(h) == 64, f"resolve-only {k}: rc={rc} hash={h[:16]}" + ("" if rc == 0 else f"\n{err[-800:]}"))
        if any(rc != 0 for rc, _, _ in res.values()):
            return
        check(res["A"][1] == res["F"][1], "source-only difference -> same compat hash (C++)")
        check(res["A"][1] != res["M"][1], "N2 nasa9_low[2] +0.001 -> different compat hash (C++)")
        recA = os.path.join(A, f"resolved_species_{res['A'][1][:16]}.yaml")
        recM = os.path.join(M, f"resolved_species_{res['M'][1][:16]}.yaml")
        ra = fs.load_record(recA)
        check(ra["consistent"] and ra["compat_recomputed"] == res["A"][1],
              f"Python recomputes the C++ compat hash ({ra['compat_recomputed'][:16]}) problems={ra['problems']}")
        with open(recA, "rb") as f:
            import hashlib
            check(ra["integrity"] == hashlib.sha256(f.read()).hexdigest(), "integrity = sha256 of the whole record file")
        check([s["source"] for s in ra["species"]] == ["builtin", "builtin"], "record keeps source as provenance")
        rm = fs.load_record(recM)
        bad = fs.compare_signatures(fs.signature_from_record(ra), fs.signature_from_record(rm))
        check(bad == ["N2.nasa9_low[2] 6.0827383600000003 vs 6.0837383600000004"] or
              (len(bad) == 1 and bad[0].startswith("N2.nasa9_low[2]")), f"record-based comparison shows N2.nasa9_low[2]: {bad}")
        # (b) 取り違え: M の記録を A の記録名で置く
        X = os.path.join(tmp, "X"); os.makedirs(X)
        swapped = os.path.join(X, os.path.basename(recA))
        shutil.copy(recM, swapped)
        rx = fs.load_record(swapped)
        check(not rx["consistent"] and any("mixed-up" in p for p in rx["problems"]), f"(b) swapped record detected: {rx['problems']}")
        check(rx["integrity"] != ra["integrity"], "(b) integrity hash differs from the original record")
        # 改竄: 係数を書き換えて compat_hash はそのまま
        with open(recA) as f:
            t = f.read()
        edited = os.path.join(X, "edited.yaml")
        with open(edited, "w") as f:
            f.write(t.replace("nasa9_low: [", "nasa9_low: [1", 1))
        re_ = fs.load_record(edited)
        check(not re_["consistent"] and any("recomputed" in p for p in re_["problems"]), f"edited record detected: {re_['problems']}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    a = ap.parse_args()
    test_compare()
    if a.forge and os.path.exists(a.forge):
        test_records(a.forge)
    else:
        print("[SKIP] record tests need a forge binary (--forge BIN or FORGE_BIN)")
    print("ALL PASSED" if FAIL == 0 else f"FAILED ({FAIL})")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
