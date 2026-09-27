#!/usr/bin/env python3
"""新規初期場 (設計 runner の IC 生成) の化学種属性 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #3b, §6 V1 (f))。

  FORGE_CUDA_BLOCKSIZE=256 design/.venv-opt/bin/python design/tests/run_species_attrs_ic_tests.py --forge BIN [--keep]
  (BIN は --resolve-species 対応の forge。1 step の CFD を 2 回回す)

case/44 va3 dry lumpX の problem で `runner_axismach.prepare` (prepare-only 相当) を一時ディレクトリに作り、
  (f-A) 宛先 thermoHrefTemp 298.15、IC の h_ref_T = 298.15 → IC に属性 (species_input_unverified=0) が付き、
        既定 (環境変数なし) のソルバが「species_hash matches」で 1 step 走り、res_1 も species_input_unverified=0
  (f-B) 同じ宛先で IC の h_ref_T だけ 0 → IC 生成 (stamp) が datum 不一致でエラー、属性は付かない。
        既定のソルバは UNVERIFIABLE で起動前に停止 (#3c で属性なしは既定で停止; 以前は FORGE_REQUIRE_VERIFIED_SPECIES=1 で先取りしていた)
  (MW)  組成を作った MW が記録と違う (MIXDRY ×(1+1e-6)) → エラー・属性なし
  (ord) 輸送種の順序が記録と違う → エラー・属性なし
  (e)   照合関数のエネルギー式: IC と同じ datum シフト → 問題なし、シフトを落とした e(T) → 内部エネルギー不一致を検出
規約: ok/FAIL、失敗があれば非ゼロ終了。
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
REPO = Path(__file__).resolve().parents[2]
PROBLEM = REPO / "case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_dry_lumpX.yaml"

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def attrs(h5):
    with h5py.File(h5, "r") as f:
        return {k: (f.attrs[k].decode() if isinstance(f.attrs[k], bytes) else f.attrs[k]) for k in f.attrs if k.startswith("species")}


def run_forge(forge, d):
    env = dict(os.environ)
    env.setdefault("FORGE_CUDA_BLOCKSIZE", "256")
    env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
    env.pop("FORGE_REQUIRE_VERIFIED_SPECIES", None)   # 撤去済み (#3c)。残っていても効かないが試験の前提を明確にする
    cfg = (Path(d) / "solverConfig.yaml").read_text()
    cfg = re.sub(r"nStepOuter: *\d+", "nStepOuter: 1", cfg)
    cfg = re.sub(r"outStepInterval: *\d+", "outStepInterval: 1", cfg)
    (Path(d) / "solverConfig.yaml").write_text(cfg)
    p = subprocess.run([forge], cwd=d, capture_output=True, text=True, env=env)
    (Path(d) / "forge_run.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    return p.returncode, p.stdout + p.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    if not a.forge:
        raise SystemExit("--forge BIN (or FORGE_BIN) is required")
    os.environ["FORGE_BIN"] = os.path.abspath(a.forge)
    from forge_design.evaluate import runner_axismach as RA
    from forge_design.evaluate.ic import _forge_species, paste_isentropic_ic, stamp_isentropic_ic_species
    from forge_design.probdef import load_problem
    fs = _forge_species()
    root = Path(tempfile.mkdtemp(prefix="forge_species_ic_"))
    try:
        p = load_problem(PROBLEM)
        dA = root / "A"
        RA.prepare(PROBLEM, dA, nsteps=1)
        atA = attrs(dA / "nozzle.h5")
        check("(f-A) runner prepare: IC stamped (species_input_unverified=0, record in the run)",
              atA.get("species_input_unverified") == 0 and len(str(atA.get("species_hash", ""))) == 64
              and (dA / str(atA.get("species_record_file"))).exists(), str({k: str(v)[:16] for k, v in atA.items()}))
        rc, out = run_forge(os.environ["FORGE_BIN"], dA)
        at1 = attrs(dA / "res_1.h5") if (dA / "res_1.h5").exists() else {}
        check("(f-A) solver (default) starts verified, res_1 species_input_unverified=0",
              rc == 0 and "species_hash matches" in out and at1.get("species_input_unverified") == 0, f"rc={rc}")

        # (f-B) 同じ宛先 (thermoHrefTemp 298.15) で IC の h_ref_T だけ 0
        design = RA.design_chain(p)
        wall, scale = design["wall"], float(p.spec["r_throat"])
        dB = root / "B"; dB.mkdir()
        for fn in ("solverConfig.yaml", "bcondConfig.yaml", "species_meta.yaml", "probe.yaml", "nozzle.h5"):
            shutil.copy(dA / fn, dB / fn)
        for fn in dA.glob("species_db*.yaml"):      # 外部 DB の生エントリがあるときだけ (lump は config に書く; plan #9)
            shutil.copy(fn, dB / fn.name)
        for fn in dA.glob("resolved_species_*.yaml"):
            shutil.copy(fn, dB / fn.name)
        L = p.species_layout(); species = list(L.species); MW = [float(L.entries[s].MW) for s in species]
        Y = RA._tp_species_Y(p)
        paste_isentropic_ic(dB / "nozzle.h5", wall, scale, float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                            gas=p.gas_model, h_ref_T=0.0, species_Y=Y)
        try:
            stamp_isentropic_ic_species(dB / "nozzle.h5", dB, p.gas_model, 0.0, species, MW, Y)
            err = ""
        except fs.SpeciesCheckError as e:
            err = str(e)
        check("(f-B) IC h_ref_T=0 vs destination thermoHrefTemp 298.15 -> IC generation error, no attributes",
              "datum" in err and attrs(dB / "nozzle.h5") == {}, err.strip().splitlines()[-1] if err else "no error raised")
        rc, out = run_forge(os.environ["FORGE_BIN"], dB)
        check("(f-B) solver (default) stops before computing (UNVERIFIABLE)",
              rc != 0 and "UNVERIFIABLE" in out and not (dB / "res_1.h5").exists(), f"rc={rc}")

        # (MW) / (ord): 組成を作った MW・順序が記録と違う
        for label, sp, mw in (("MW", species, [MW[0] * (1 + 1e-6)] + MW[1:]), ("ord", species[::-1], MW[::-1])):
            h5 = dB / "nozzle.h5"
            try:
                stamp_isentropic_ic_species(h5, dB, p.gas_model, 298.15, sp, mw, Y)
                err = ""
            except fs.SpeciesCheckError as e:
                err = str(e)
            key = "MW[" if label == "MW" else "species order"
            check(f"({label}) IC {label} differs from the destination record -> error, no attributes",
                  key in err and attrs(h5) == {}, err.strip().splitlines()[-1] if err else "no error raised")

        # (e) 照合関数のエネルギー式 (datum の比較とは独立に e(T) の食い違いを検出する)
        rec = fs.load_record(dA / str(atA["species_record_file"]))
        gas = p.gas_model
        ok_mix = [("inflow", Y, gas.R, lambda T: gas.h_mass(T) - gas.R * T - float(gas.h_mass(298.15)[0]))]
        bad_mix = [("inflow", Y, gas.R, lambda T: gas.h_mass(T) - gas.R * T)]
        pr_ok = fs.check_ic_against_record(rec, 298.15, species, MW, ok_mix)
        pr_bad = fs.check_ic_against_record(rec, 298.15, species, MW, bad_mix)
        check("(e) energy check: same datum -> no problem; energy without the datum shift -> internal energy mismatch",
              pr_ok == [] and any("internal energy" in x for x in pr_bad), f"{pr_ok} / {pr_bad[:1]}")
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("\nALL PASS" if FAIL == 0 else f"\n{FAIL} FAILED")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
