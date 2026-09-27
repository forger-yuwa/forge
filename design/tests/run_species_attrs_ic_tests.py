#!/usr/bin/env python3
"""新規初期場 (設計 runner の IC 生成) の化学種属性 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #3b, §6 V1 (f))。

  FORGE_CUDA_BLOCKSIZE=256 design/.venv-opt/bin/python design/tests/run_species_attrs_ic_tests.py --forge BIN [--keep]
  (BIN は --resolve-species 対応の forge。1 step の CFD を 2 回回す)

case/44 va3 dry lumpX の problem で `runner_axismach.prepare` (prepare-only 相当) を一時ディレクトリに作り、
  (f-A) 宛先 thermoHrefTemp 298.15、IC の h_ref_T = 298.15 → IC に属性 (species_input_unverified=0) が付き、
        既定 (環境変数なし) のソルバが「species_hash matches」で 1 step 走り、res_1 も species_input_unverified=0
  (f-B) 同じ宛先で IC の h_ref_T だけ 0 → IC 生成 (stamp) が datum 不一致でエラー、属性は付かない。
        既定のソルバは UNVERIFIABLE で起動前に停止 (#3c で属性なしは既定で停止)
  (MW)  組成を作った MW が記録と違う (MIXDRY ×(1+1e-6)) → エラー・属性なし
  (ord) 輸送種の順序が記録と違う → エラー・属性なし
  (e)   照合関数のエネルギー式: IC と同じ datum シフト → 問題なし、シフトを落とした e(T) → 内部エネルギー不一致を検出
  (f-C) 宛先を解決できない (FORGE_BIN が --resolve-species を持たない旧バイナリ) → IC の属性付与は既定で停止 (#3c)、
        FORGE_ALLOW_UNVERIFIED_SPECIES=1 で 'unverified' (属性なし)。旧バイナリは起動しない
  (stage) runner_axismach の段間継承 (_restart_same_mesh = restart_field.py): 検証済みの SRC → 属性を継承、
        未検証の SRC (属性なし) → 既定で停止し宛先不変、FORGE_ALLOW_UNVERIFIED_SPECIES=1 で属性なしで写す
  --no-cfd: ソルバを起動する 2 項目 ((f-A) の 1 step と (f-B) の起動拒否) を省く (Python 側の変更だけを確かめるとき)
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
    ap.add_argument("--no-cfd", action="store_true", help="ソルバを起動する項目 ((f-A) 1 step・(f-B) 起動拒否) を省く")
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
        dA0 = root / "A0"                        # (stage) 用に 1 step 前の状態を残す (resolve 記録・config・IC)
        shutil.copytree(dA, dA0)
        if a.no_cfd:
            print("SKIP (f-A) solver 1 step (--no-cfd)")
        else:
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
        if a.no_cfd:
            print("SKIP (f-B) solver launch (--no-cfd)")
        else:
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

        # (f-C) 宛先を解決できない (旧バイナリ) → IC の属性付与は既定で停止、その実行だけの許可で属性なし
        marker = root / "old_binary_launched"
        fake = root / "fake_old_forge"
        fake.write_text(f"#!/bin/sh\ntouch '{marker}'\n"); fake.chmod(0o755)
        saved = {k: os.environ.get(k) for k in ("FORGE_BIN", "FORGE_ALLOW_UNVERIFIED_SPECIES")}
        h5 = dB / "nozzle.h5"
        try:
            os.environ["FORGE_BIN"] = str(fake)
            os.environ.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
            try:
                stamp_isentropic_ic_species(h5, dB, p.gas_model, 298.15, species, MW, Y)
                err = ""
            except fs.SpeciesCheckError as e:
                err = str(e)
            check("(f-C) destination unresolvable (old binary) -> IC stamping stops by default, no attributes, 2-way guidance",
                  "旧バイナリ" in err and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in err and "regenerate the initial field" in err
                  and attrs(h5) == {} and not marker.exists(), err.strip().splitlines()[0][:160] if err else "no error raised")
            os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
            try:
                st_ = stamp_isentropic_ic_species(h5, dB, p.gas_model, 298.15, species, MW, Y)
            except fs.SpeciesCheckError as e:
                st_ = f"error: {e}"
            check("(f-C) same with FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> 'unverified', no attributes, old binary not launched",
                  st_ == "unverified" and attrs(h5) == {} and not marker.exists(), str(st_)[:160])
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

        # (stage) runner_axismach の段間継承 (restart_field.py を子プロセスで呼ぶ; 子の環境は RA._ENV)
        S1 = root / "S1"; S1.mkdir()                       # 検証済みの SRC: 属性付きの IC を res として置き、記録を隣に
        shutil.copy(dA0 / "nozzle.h5", S1 / "res_1.h5")
        for fn in dA0.glob("resolved_species_*.yaml"):
            shutil.copy(fn, S1 / fn.name)
        S2 = root / "S2"; S2.mkdir()                       # 未検証の SRC: 属性なし
        shutil.copy(dA0 / "nozzle.h5", S2 / "res_1.h5")
        fs.write_species_attrs(S2 / "res_1.h5", None)
        env_saved = RA._ENV.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
        try:
            dC = root / "C"; shutil.copytree(dA0, dC)
            fs.write_species_attrs(dC / "nozzle.h5", None)
            RA._restart_same_mesh(S1 / "res_1.h5", dC / "nozzle.h5")
            atC = attrs(dC / "nozzle.h5")
            check("(stage) verified SRC -> restart_field inherits (same species_hash, species_input_unverified=0)",
                  atC.get("species_hash") == atA.get("species_hash") and atC.get("species_input_unverified") == 0, str(atC)[:120])
            dD = root / "D"; shutil.copytree(dA0, dD)
            atD0 = attrs(dD / "nozzle.h5")
            try:
                RA._restart_same_mesh(S2 / "res_1.h5", dD / "nozzle.h5")
                err = ""
            except RuntimeError as e:
                err = str(e)
            check("(stage) unverified SRC -> stops by default, destination attributes unchanged",
                  "REFUSED" in err and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in err and attrs(dD / "nozzle.h5") == atD0,
                  err.strip().splitlines()[0][:160] if err else "no error raised")
            RA._ENV["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
            try:
                RA._restart_same_mesh(S2 / "res_1.h5", dD / "nozzle.h5")
                err = ""
            except RuntimeError as e:
                err = str(e)
            check("(stage) unverified SRC + FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> copied without species attributes",
                  err == "" and attrs(dD / "nozzle.h5") == {}, err[:160])
        finally:
            RA._ENV.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
            if env_saved is not None:
                RA._ENV["FORGE_ALLOW_UNVERIFIED_SPECIES"] = env_saved
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("\nALL PASS" if FAIL == 0 else f"\n{FAIL} FAILED")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
