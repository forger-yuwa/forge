#!/usr/bin/env python3
"""区間可変 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #13-1) のソルバ + Python 鏡像の試験。GPU 不要。

  python3 solver_density_cuda/tests/unit/test_thermo_intervals.py --forge BIN [--seed-run RUN_0509_DIR] [--keep]

`forge --resolve-species` (と起動時の拒否を見るための素の `forge`) だけを使う (計算はしない)。seed run は config を読むだけ。

  (F)  G1-f: 3 区間種 (CEA thermo.inp の N2 3 区間 200/1000/6000/20000 を外部 DB の Tbounds/nasa9_intervals で与える) と
       区切りの違う lump (O2 + Tmid 1500 の試験種 → 和集合 200/1000/1500/6000) を含む config で、ソルバが書いた記録を Python が読み
       (load_record)、互換性ハッシュの再計算がソルバの値と一致、schema = forge_resolved_species_v1_nint・外挿規約 = 区間可変。
       記録を署名にして比べると同一なら差なし、第 3 区間の係数を 1 つ変えると N2C3.nasa9_intervals[2][0] を示す。
  (T)  Python の _TPGas (total_quantities.py) が 3 区間種・和集合 lump をソルバと同じ規約で評価する: 独立に書いた参照
       (区間 k = Tb[k] <= T < Tb[k+1]、端で cp 一定・h 線形・s° 対数の外挿) と 100–25000 K・全区切り ±1e-9 K で相対 1e-12。
  (S)  speciesDBFile に区間可変の書式を書いた run 同士の species_signature/compare_signatures (interp_field・restart の署名) が
       同一 → 差なし、第 3 区間の係数違い → 該当キーを示す。
  (U)  2 区間だけの config (seed run_0509) は schema forge_resolved_species_v1・互換性ハッシュ 4378b7d78339ba27 のまま。
  (R)  上限超過の起動時拒否: 外部 DB に 4 区間の種 → `forge --resolve-species` と素の `forge` の起動がともに非ゼロ終了し
       THERMO_MAX_INTERVALS を示す; 和集合が 4 区間になる lump (N2 + Tlo 298.15 の 3 区間種) も拒否。Python の読込も拒否。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, math, os, re, shutil, subprocess, sys, tempfile

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TOOLS)
import forge_species as fs  # noqa: E402
from total_quantities import _TPGas  # noqa: E402

RU = 8.314462618
FAIL = 0
# CEA thermo.inp (SHA-256 b1bc0707...) の N2・O2
N2C = [[22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12, 710.846086, -10.76003744],
       [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86640027],
       [831013916.0, -642073.354, 202.0264635, -0.03065092046, 2.486903333e-06, -9.70595411e-11, 1.437538881e-15, 4938707.04, -1672.09974]]
O2C = [[-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12, -3391.45487, 18.4969947],
       [-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16, -16890.10929, 17.38716506]]


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def test_db():
    return {
        "N2C3": {"MW": 0.0280134, "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tbounds": [200.0, 1000.0, 6000.0, 20000.0], "nasa9_intervals": N2C},
        "TMID1500": {"MW": 0.0280134, "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tlo": 200.0, "Tmid": 1500.0, "Thi": 6000.0,
                     "nasa9_low": N2C[0], "nasa9_high": N2C[1]},
        "N2T298": {"MW": 0.0280134, "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tbounds": [298.15, 1000.0, 6000.0, 20000.0], "nasa9_intervals": N2C},
    }


# ---- 独立な参照 (C++ thermo_d.cuh の規約を式から書き直したもの; forge_species / _TPGas は使わない) ----
def ref_props(Tb, co, MW, T):
    """区間 k = Tb[k] <= T < Tb[k+1] (区切りちょうどは上)、[Tlo, Thi] の外は端の cp 一定・h 線形・s° 対数。質量基準。"""
    def raw(Tc):
        k = 0
        for j in range(1, len(co)):
            if not Tc < Tb[j]:
                k = j
        a = co[k]
        cp = a[0] / Tc**2 + a[1] / Tc + a[2] + a[3] * Tc + a[4] * Tc**2 + a[5] * Tc**3 + a[6] * Tc**4
        h = (-a[0] / Tc**2 + a[1] * math.log(Tc) / Tc + a[2] + a[3] * Tc / 2 + a[4] * Tc**2 / 3 + a[5] * Tc**3 / 4
             + a[6] * Tc**4 / 5 + a[7] / Tc) * Tc
        s = (-a[0] / (2 * Tc**2) - a[1] / Tc + a[2] * math.log(Tc) + a[3] * Tc + a[4] * Tc**2 / 2 + a[5] * Tc**3 / 3
             + a[6] * Tc**4 / 4 + a[8])
        return cp * RU / MW, h * RU / MW, s * RU / MW
    Tc = min(max(T, Tb[0]), Tb[-1])
    cp, h, s = raw(Tc)
    if T != Tc:
        h += cp * (T - Tc)
        s += cp * math.log(T / Tc)
    return cp, h, s


class Ctx:
    def __init__(self, a, root):
        self.a, self.root, self.n = a, root, 0

    def make(self, species_yaml, db=None):
        self.n += 1
        d = os.path.join(self.root, f"c{self.n:02d}")
        os.makedirs(d)
        t = open(os.path.join(self.a.seed_run, "solverConfig.yaml")).read()
        t, k = re.subn(r"species: *\[[^\]]*\]", "species: " + species_yaml, t, count=1)
        assert k == 1
        t = re.sub(r",? *speciesDBFile: *\"?[^,}\"]*\"?", "", t, count=1)
        if db is not None:
            t = t.replace("thermoHrefTemp:", "speciesDBFile: \"test_db.yaml\", thermoHrefTemp:", 1)
            with open(os.path.join(d, "test_db.yaml"), "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
        with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
            f.write(t)
        return d

    def resolve(self, d, args=("--resolve-species",)):
        p = subprocess.run([self.a.forge, *args], cwd=d, capture_output=True, text=True, timeout=300)
        lines = p.stdout.strip().splitlines()
        h = lines[-1].strip() if lines else ""
        m = re.search(r"\[species\] record (\S+) \(sha256 ([0-9a-f]{64})\)", p.stderr)
        rec = fs.load_record(os.path.join(d, m.group(1))) if (p.returncode == 0 and m) else None
        return p.returncode, h, rec, p.stderr + p.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"), required=os.environ.get("FORGE_BIN") is None)
    ap.add_argument("--seed-run", default=os.path.join(REPO, "case", "44.vitiated_air_wt", "run_0509_va3_M4.19_Lc8_dry_lumpX"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    a.forge = os.path.abspath(a.forge)
    root = tempfile.mkdtemp(prefix="forge_thermo_intervals_")
    C = Ctx(a, root)
    try:
        db = test_db()
        lump = "{name: LMP, lump: {O2: 0.3, TMID1500: 0.7}, basis: mole}"
        # ---- F ----
        d = C.make(f"[{lump}, N2C3, H2O]", db=db)
        rc, h, rec, err = C.resolve(d)
        check(rc == 0 and rec is not None, f"F resolve with a 3-interval species and a union lump (rc={rc})" + ("" if rc == 0 else f": {err[-600:]}"))
        if rec is not None:
            check(rec["consistent"] and rec["compat_recomputed"] == h, f"F Python load_record: consistent, recomputed compat hash == solver ({h[:16]}): {rec['problems']}")
            check(rec["schema"] == fs.SPECIES_RECORD_SCHEMA_NINT and rec["extrapolation"] == fs.SPECIES_RECORD_EXTRAPOLATION_NINT,
                  f"F schema {rec['schema']!r} and the n-interval extrapolation convention")
            sp = {e["name"]: e for e in rec["species"]}
            n3 = sp["N2C3"]
            check(n3.get("Tbounds") == [200.0, 1000.0, 6000.0, 20000.0] and n3["nasa9_intervals"] == N2C and "nasa9_low" not in n3,
                  "F record keeps N2C3's bounds and all 3 intervals bit-exact (Tbounds/nasa9_intervals)")
            lm = sp["LMP"]
            check(lm.get("Tbounds") == [200.0, 1000.0, 1500.0, 6000.0] and "union" in lm["lump"]["synthesis"],
                  f"F union lump: bounds {lm.get('Tbounds')} and synthesis rule '{lm['lump']['synthesis'][:50]}...'")
            check(sp["H2O"].get("Tmid") == 1000.0 and "Tbounds" not in sp["H2O"], "F 2-interval H2O keeps the 2-interval keys in the same record")
            sig = fs.signature_from_record(rec)
            check(fs.compare_signatures(sig, fs.signature_from_record(rec)) == [], "F compare_signatures: same record -> no difference")
            rec2 = fs.load_record(rec["path"])
            for e in rec2["species"]:
                if e["name"] == "N2C3":
                    e["nasa9_intervals"][2][0] *= (1 + 1e-9)
            diff = fs.compare_signatures(sig, fs.signature_from_record(rec2))
            check(len(diff) == 1 and diff[0].startswith("N2C3.nasa9_intervals[2][0]"), f"F compare_signatures shows the changed 3rd-interval coefficient: {diff}")

            # ---- T: _TPGas vs 独立参照 ----
            th = fs._thermo_from_record(rec, "test")
            Ts = sorted(set([100.0, 150.0, 199.9] + list(np.arange(200.0, 25000.0, 37.0))
                            + [b + dd for b in (200.0, 1000.0, 1500.0, 6000.0, 20000.0) for dd in (-1e-9, 0.0, 1e-9)]))
            T = np.array(Ts)
            for name in ("N2C3", "LMP", "H2O"):
                e = th["species"][name]
                Tb, co = fs.nasa9_intervals(e)
                g = _TPGas({name: e}, [name], 0.0)
                cp, hh, s0 = g._props(e, T)
                ref = np.array([ref_props(Tb, co, e["MW"], float(t)) for t in T])
                hmax = np.max(np.abs(ref[:, 1]))
                w = max(np.max(np.abs(cp - ref[:, 0]) / np.abs(ref[:, 0])), np.max(np.abs(hh - ref[:, 1])) / hmax,
                        np.max(np.abs(s0 - ref[:, 2]) / np.abs(ref[:, 2])))
                check(w <= 1e-12, f"T _TPGas {name} ({len(co)} intervals) vs independent reference, 100-25000 K and breakpoints +-1e-9 K: max rel {w:.2e}")

        # ---- S: speciesDBFile の区間可変の書式での署名 (interp_field / restart の照合) ----
        dA = C.make("[N2C3, TMID1500]", db=db)
        dbB = test_db()
        dbB["N2C3"]["nasa9_intervals"] = [list(x) for x in N2C]
        dbB["N2C3"]["nasa9_intervals"][2][0] = N2C[2][0] * (1 + 1e-9)
        dB = C.make("[N2C3, TMID1500]", db=dbB)
        dA2 = C.make("[N2C3, TMID1500]", db=test_db())
        sA, sSame, sMod = fs.species_signature(dA), fs.species_signature(dA2), fs.species_signature(dB)
        check(fs.compare_signatures(sA, sSame) == [], "S species_signature from a Tbounds/nasa9_intervals speciesDBFile: same DB -> no difference")
        dd = fs.compare_signatures(sA, sMod)
        check(len(dd) == 1 and dd[0].startswith("N2C3.nasa9_intervals[2][0]"), f"S 3rd-interval coefficient difference shown: {dd}")

        # ---- U: 2 区間だけ ----
        dU = os.path.join(root, "u0509")
        os.makedirs(dU)
        for f in ("solverConfig.yaml", "species_db.yaml"):
            shutil.copy(os.path.join(a.seed_run, f), dU)
        rcu, hu, recu, _ = C.resolve(dU)
        check(rcu == 0 and hu.startswith("4378b7d78339ba27") and recu is not None and recu["schema"] == fs.SPECIES_RECORD_SCHEMA
              and recu["extrapolation"] == fs.SPECIES_RECORD_EXTRAPOLATION, f"U run_0509 config: hash {hu[:16]}, schema {recu and recu['schema']!r}")

        # ---- R: 上限超過の起動時拒否 ----
        bad = test_db()
        bad["FOURINT"] = {"MW": 0.0280134, "Tbounds": [200.0, 1000.0, 3000.0, 6000.0, 20000.0], "nasa9_intervals": N2C + [N2C[2]]}
        dR = C.make("[N2, H2O]", db=bad)
        rcr, _, _, errr = C.resolve(dR)
        check(rcr != 0 and "THERMO_MAX_INTERVALS" in errr and "FOURINT" in errr, f"R --resolve-species with a 4-interval DB species -> refused (rc={rcr})")
        p = subprocess.run([a.forge], cwd=dR, capture_output=True, text=True, timeout=300)
        out = p.stdout + p.stderr
        check(p.returncode != 0 and "THERMO_MAX_INTERVALS" in out, f"R plain forge start with a 4-interval DB species -> refused at startup (rc={p.returncode})")
        dL = C.make("[{name: LMP, lump: {N2: 0.9, N2T298: 0.1}, basis: mole}, H2O]", db=db)
        rcl, _, _, errl = C.resolve(dL)
        check(rcl != 0 and "makes 4 intervals" in errl, f"R lump whose union has 4 intervals (N2 + Tlo 298.15 3-interval) -> refused (rc={rcl})")
        try:
            fs._thermo_entry(bad["FOURINT"])
            check(False, "R Python reader refuses a 4-interval entry")
        except ValueError as ex:
            check("THERMO_MAX_INTERVALS" in str(ex), f"R Python reader refuses a 4-interval entry: {ex}")
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("ALL PASSED" if FAIL == 0 else f"FAILED ({FAIL})")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
