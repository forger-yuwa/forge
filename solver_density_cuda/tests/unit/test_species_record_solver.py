#!/usr/bin/env python3
"""ソルバ入口 (`valueFileName`) の化学種照合と解決済み記録の実機試験 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #3a, §6 V1)。

  FORGE_CUDA_BLOCKSIZE=256 python3 solver_density_cuda/tests/unit/test_species_record_solver.py --forge BIN \
      [--seed-run case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX] [--seed-res res_24000.h5] [--keep]

実在の TP run (既定 case/44 va3 dry lumpX: 種 MIXDRY (外部 DB) + H2O (内蔵)) の入力を一時ディレクトリへ複製し、
seed の res を `restart_field.py` で nozzle.h5 に写して **1 step** だけ回す (拒否される場合は GPU 計算に入る前に終了する)。
seed run の中身は読むだけで書かない (h5 もリンクせず複製する)。

  (0)  旧場 (属性なし) → 既定で「照合不能」(UNVERIFIABLE) で停止 (非ゼロ終了)。エラー文が IC の作り直しと
       その実行だけの FORGE_ALLOW_UNVERIFIED_SPECIES=1 の 2 通りを案内する (#3c: 既定を厳密へ切り替え)
  (0e) 旧場 + FORGE_ALLOW_UNVERIFIED_SPECIES=1 → 通る。res (境界出力を含む) に species_hash / species_record_sha256 /
       species_record_file / species_input_unverified=1、記録ファイルが run に書かれ、Python の再計算 (load_record) と一致
  (0r) (0e) の res (印付き) を restart_field.py で**許可なしに**写す → 通り、DST に同じ species_hash と印 (species_input_unverified=1)
       (#3d: ツールもソルバの fieldHash == own 分岐と同じ規則)。その場を既定 (環境変数なし) のソルバで 1 step → 通り、出力も印を継承
  (R)  `forge --resolve-species` の標準出力最終行 = ソルバ起動時の species_hash
  (A)  同一 config で作った res を restart_field.py で写し属性をコピー (unverified=0) → 通る、出力 unverified=0
  (Ai) 同上で入力 unverified=1 → 通る、出力に未検証の印を継承
  (B)  外部 DB の MIXDRY nasa9_low[2] +0.001 → 停止し MIXDRY.nasa9_low[2] を表示。env を付けても停止
  (c)  種 [N2, H2O] の場 (resolve-only の記録で属性を付けた試験用の場) に外部 DB で N2 nasa9_low[2] +0.001 → 停止し N2.nasa9_low[2] を表示
  (b)  記録の取り違え (入力場の記録名に別 run の記録を置く) → 完全性ハッシュ不一致を表示して停止; Python find_record も検出
  (S)  source だけ違う (外部 DB に内蔵と同一係数の H2O を置いた場 → H2O を外して内蔵にする) → 互換性ハッシュ一致で通る (記録の source は builtin)。
       段 3 (#13-3) 前は seed の外部 DB の H2O がそのまま内蔵と同一係数だったが、段 3 で内蔵 H2O の MW が CEA の 0.01801528 に
       なったので、同一係数の H2O は内蔵の記録から作る
  (S3) 段 3 の旧記録からの継続: seed の場 (外部 DB の H2O は旧 MW 0.0180153) を、H2O を外した config (内蔵 = CEA の MW) で読む
       → 係数不一致で停止し H2O.MW を示す (env でも通さない; 移行は convert_species_field.py / restart_field.py --force-species)
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, os, re, shutil, subprocess, sys, tempfile

import h5py
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TOOLS)
import forge_species as fs  # noqa: E402

ATTRS = ("species_hash", "species_record_sha256", "species_record_file", "species_input_unverified")
FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def attrs(h5):
    with h5py.File(h5, "r") as f:
        return {k: (f.attrs[k].decode() if isinstance(f.attrs[k], bytes) else f.attrs[k]) for k in ATTRS if k in f.attrs}


def stamp(h5, a, unverified=None):
    """試験用: 属性を h5 に書く (実運用では生成処理が付ける; §4.3)。"""
    with h5py.File(h5, "r+") as f:
        for k in ATTRS:
            if k in f.attrs:
                del f.attrs[k]
        for k, v in a.items():
            if k == "species_input_unverified":
                f.attrs.create(k, int(v if unverified is None else unverified), dtype="int32")
            else:
                f.attrs[k] = str(v)


class Case:
    def __init__(self, args, root):
        self.a, self.root = args, root

    def make(self, name, seed_res=None, species=None, db_edit=None):
        """seed run の入力を複製し seed_res の保存量を nozzle.h5 に写す。species で physProp.species を、db_edit(db) で外部 DB を差し替える。"""
        d = os.path.join(self.root, name)
        os.makedirs(d)
        for fn in ("solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "species_meta.yaml", "probe.yaml", "nozzle.h5"):
            src = os.path.join(self.a.seed_run, fn)
            if os.path.exists(src):
                shutil.copy(src, d)
        # 保存量の複製だけに使う (属性は各試験が stamp() で明示的に書く)。seed の旧場・(0e) の res は未検証なので、#3c 以降の
        # restart_field は既定で停止する → その実行だけ FORGE_ALLOW_UNVERIFIED_SPECIES=1 で写す (検証済みの SRC なら従来どおり継承)
        env = dict(os.environ, FORGE_ALLOW_UNVERIFIED_SPECIES="1")
        p = subprocess.run([sys.executable, os.path.join(TOOLS, "restart_field.py"), seed_res or os.path.join(self.a.seed_run, self.a.seed_res),
                            os.path.join(d, "nozzle.h5")], capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise SystemExit(f"restart_field failed: {p.stdout}{p.stderr}")
        cfgp = os.path.join(d, "solverConfig.yaml")
        t = open(cfgp).read()
        t = re.sub(r"nStepOuter: *\d+", "nStepOuter: 1", t)
        t = re.sub(r"outStepInterval: *\d+", "outStepInterval: 1", t)
        if species is not None:
            t = re.sub(r"species: *\[[^\]]*\]", "species: [" + ", ".join(f'"{s}"' for s in species) + "]", t)
        open(cfgp, "w").write(t)
        if db_edit is not None:
            dbp = os.path.join(d, "species_db.yaml")
            db = yaml.load(open(dbp).read(), Loader=fs._StrSafeLoader) if os.path.exists(dbp) else {}
            db = db_edit(db or {})
            with open(dbp, "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
        return d

    def run(self, d, allow=False, extra_env=None):
        env = dict(os.environ)
        env.setdefault("FORGE_CUDA_BLOCKSIZE", "256")
        env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
        if allow:
            env["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
        env.update(extra_env or {})
        p = subprocess.run([self.a.forge], cwd=d, capture_output=True, text=True, env=env)
        with open(os.path.join(d, "forge_run.log"), "w") as f:
            f.write(p.stdout + "\n--- stderr ---\n" + p.stderr)
        return p.returncode, p.stdout + p.stderr

    def resolve(self, d):
        p = subprocess.run([self.a.forge, "--resolve-species"], cwd=d, capture_output=True, text=True)
        lines = p.stdout.strip().splitlines()
        return p.returncode, (lines[-1] if lines else "")


def bump_low2(name, delta):
    def f(db):
        db[name]["nasa9_low"][2] = float(db[name]["nasa9_low"][2]) + delta
        return db
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    ap.add_argument("--seed-run", default=os.path.join(REPO, "case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX"))
    ap.add_argument("--seed-res", default="res_24000.h5")
    ap.add_argument("--keep", action="store_true", help="一時ディレクトリを残す")
    a = ap.parse_args()
    a.forge = os.path.abspath(a.forge) if a.forge else a.forge
    a.seed_run = os.path.abspath(a.seed_run)
    if not a.forge or not os.path.exists(a.forge):
        raise SystemExit("--forge BIN (or FORGE_BIN) is required")
    if not os.path.exists(os.path.join(a.seed_run, a.seed_res)):
        raise SystemExit(f"seed not found: {a.seed_run}/{a.seed_res}")
    root = tempfile.mkdtemp(prefix="forge_species_solver_")
    C = Case(a, root)
    try:
        # (0) 旧場
        d0 = C.make("old")
        rc, out = C.run(d0)
        check(rc != 0 and "UNVERIFIABLE" in out and not os.path.exists(os.path.join(d0, "res_1.h5")),
              f"(0) old field without attributes -> UNVERIFIABLE, stops by default (rc={rc})")
        check("regenerate the initial field" in out and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out
              and "THIS invocation only" in out and "Allowed for this invocation" not in out,
              "(0) the error names both remedies (regenerate the IC / allow this invocation only) and does not claim it was allowed")
        # (0e) 旧場 + env
        d1 = C.make("old_env")
        rc, out = C.run(d1, allow=True)
        r1 = os.path.join(d1, "res_1.h5")
        ok = rc == 0 and os.path.exists(r1)
        check(ok, f"(0e) old field + FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> runs (rc={rc})")
        if not ok:
            print(out[-3000:])
            raise _Abort()
        check("Allowed for this invocation by FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out,
              "(0e) the warning says it was allowed by the env var that was actually set")
        at = attrs(r1)
        check(at.get("species_input_unverified") == 1 and len(at.get("species_hash", "")) == 64,
              f"(0e) res_1.h5 carries species_hash and species_input_unverified=1: {at}")
        bfiles = [f for f in os.listdir(d1) if re.fullmatch(r"res_.+_\d+_1\.h5", f)]
        check(bool(bfiles) and all(attrs(os.path.join(d1, f)) == at for f in bfiles),
              f"(0e) boundary outputs carry the same attributes: {bfiles}")
        H = at["species_hash"]
        rec, why = fs.find_record(r1)
        check(rec is not None and rec["consistent"] and rec["compat_recomputed"] == H
              and rec["provenance"].get("input_status") == "unverified_env",
              f"(0e) record found and verified by Python (integrity + recomputed compat hash): {why}")
        # (0r) 印付きの res → restart_field (許可なし) → 印を継承 → 既定のソルバで 1 step (#3d)
        dr = os.path.join(root, "old_env_restart")
        os.makedirs(dr)
        for fn in ("solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "species_meta.yaml", "probe.yaml", "nozzle.h5"):
            if os.path.exists(os.path.join(d1, fn)):
                shutil.copy(os.path.join(d1, fn), dr)
        env = dict(os.environ); env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
        p = subprocess.run([sys.executable, os.path.join(TOOLS, "restart_field.py"), r1, os.path.join(dr, "nozzle.h5"), "--forge", a.forge],
                           capture_output=True, text=True, env=env)
        atr = attrs(os.path.join(dr, "nozzle.h5"))
        check(p.returncode == 0 and atr.get("species_hash") == H and atr.get("species_input_unverified") == 1,
              f"(0r) restart_field from the marked res without permission -> passes, hash and mark inherited (rc={p.returncode}, {atr})")
        if p.returncode != 0:
            print((p.stdout + p.stderr)[-2000:])
        rc, out = C.run(dr)
        atr1 = attrs(os.path.join(dr, "res_1.h5")) if rc == 0 else {}
        check(rc == 0 and "species_hash matches" in out and atr1.get("species_input_unverified") == 1 and atr1.get("species_hash") == H,
              f"(0r) default solver (no env) on the inherited field -> runs, output keeps the mark (rc={rc}, {atr1.get('species_input_unverified')})")
        # (R) resolve-only
        dR = C.make("resolve")
        rc, h = C.resolve(dR)
        check(rc == 0 and h == H, f"(R) --resolve-species prints the solver's species_hash ({h[:16]} vs {H[:16]})")
        # (A) 同一 config の res から restart (属性をコピー、unverified=0)
        dA = C.make("A", seed_res=r1)
        stamp(os.path.join(dA, "nozzle.h5"), at, unverified=0)
        rc, out = C.run(dA)
        atA = attrs(os.path.join(dA, "res_1.h5")) if rc == 0 else {}
        check(rc == 0 and "species_hash matches" in out and atA.get("species_input_unverified") == 0 and atA.get("species_hash") == H,
              f"(A) same species -> runs, output verified (rc={rc}, {atA.get('species_input_unverified')})")
        dAi = C.make("A_inherit", seed_res=r1)
        stamp(os.path.join(dAi, "nozzle.h5"), at, unverified=1)
        rc, out = C.run(dAi)
        atAi = attrs(os.path.join(dAi, "res_1.h5")) if rc == 0 else {}
        check(rc == 0 and atAi.get("species_input_unverified") == 1, f"(Ai) unverified mark on the input is inherited (rc={rc})")
        # (B) 外部 DB の MIXDRY low[2] +0.001 (入力場の記録を隣に置く)
        dB = C.make("B", seed_res=r1, db_edit=bump_low2("MIXDRY", 0.001))
        stamp(os.path.join(dB, "nozzle.h5"), at)
        shutil.copy(os.path.join(d1, at["species_record_file"]), dB)
        rc, out = C.run(dB)
        check(rc != 0 and "MIXDRY.nasa9_low[2]" in out and not os.path.exists(os.path.join(dB, "res_1.h5")),
              f"(B) external DB MIXDRY nasa9_low[2] +0.001 -> stops, coefficient shown (rc={rc})")
        m = re.search(r"^ +MIXDRY\.nasa9_low\[2\].*$", out, re.M)
        print("       " + (m.group(0).strip() if m else "(no line)"))
        rc, out = C.run(dB, allow=True)
        check(rc != 0 and "never allowed" in out, f"(B) mismatch is not allowed by FORGE_ALLOW_UNVERIFIED_SPECIES (rc={rc})")
        # (c) N2: 種 [N2, H2O] の記録を resolve-only で作り、試験用の場に属性を付ける → N2 low[2] +0.001 の外部 DB で読む
        dN = C.make("N2_src", species=["N2", "H2O"], db_edit=lambda db: {})
        rc, hN = C.resolve(dN)
        recN = [f for f in os.listdir(dN) if f.startswith("resolved_species_")][0]
        with open(os.path.join(dN, recN), "rb") as f:
            import hashlib
            shaN = hashlib.sha256(f.read()).hexdigest()
        recNsp = {e["name"]: e for e in fs.load_record(os.path.join(dN, recN))["species"]}
        builtinN2 = recNsp["N2"]
        if "Tbounds" in builtinN2:   # 段 3 (#13-3) から内蔵 N2 は CEA の 3 区間 → 外部 DB も区間可変の書式で全区間を書く
            n2db = {k: builtinN2[k] for k in ("MW", "LJ_sigma", "LJ_eps_kB", "Tbounds")}
            n2db["nasa9_intervals"] = [list(r_) for r_ in builtinN2["nasa9_intervals"]]
            n2db["nasa9_intervals"][0][2] += 0.001
            n2key = "N2.nasa9_intervals[0][2]"
        else:
            n2db = {k: builtinN2[k] for k in ("MW", "LJ_sigma", "LJ_eps_kB", "Tlo", "Tmid", "Thi", "nasa9_low", "nasa9_high")}
            n2db["nasa9_low"] = list(n2db["nasa9_low"]); n2db["nasa9_low"][2] += 0.001
            n2key = "N2.nasa9_low[2]"
        dc = C.make("c_N2", species=["N2", "H2O"], db_edit=lambda db: {"N2": n2db})
        stamp(os.path.join(dc, "nozzle.h5"), {"species_hash": hN, "species_record_sha256": shaN, "species_record_file": recN,
                                              "species_input_unverified": 0})
        shutil.copy(os.path.join(dN, recN), dc)
        rc, out = C.run(dc)
        diffs = re.findall(r"^ +(\S+): field .* vs current .*$", out, re.M)
        check(rc != 0 and diffs == [n2key], f"(c) external DB N2 low[2] +0.001 -> stops, only {n2key} shown: {diffs}")
        # (b) 記録の取り違え: 入力場の記録名に別 run (MIXDRY 変更 config の resolve-only) の記録を置く
        dX = C.make("b_swap_src", db_edit=bump_low2("MIXDRY", 0.001))
        C.resolve(dX)
        other = [f for f in os.listdir(dX) if f.startswith("resolved_species_")][0]
        db_ = C.make("b_swap", seed_res=r1, db_edit=bump_low2("MIXDRY", 0.001))
        stamp(os.path.join(db_, "nozzle.h5"), at)
        # 正しい記録の SHA 付きの別名 (resolved_species_<hash16>_<sha16>.yaml; seed から複製される) が残っていると、ソルバと find_record は
        # そちらを内容の SHA で見つけて正しい記録を使うので、すり替えを検知する経路を通らない (2026-10-06: 熱物性の不一致で止まるだけだった)。
        # すり替えの検知を試すため、場が名指す記録以外の記録ファイルを消す。
        for f in os.listdir(db_):
            if f.startswith("resolved_species_") and f != at["species_record_file"]:
                os.remove(os.path.join(db_, f))
        shutil.copy(os.path.join(dX, other), os.path.join(db_, at["species_record_file"]))
        rc, out = C.run(db_)
        check(rc != 0 and "integrity" in out and "cannot be identified" in out,
              f"(b) swapped record -> integrity mismatch reported, stops (rc={rc})")
        rec, why = fs.find_record(os.path.join(db_, "nozzle.h5"))
        check(rec is None and "integrity mismatch" in why, f"(b) Python find_record detects the swapped record: {why[:120]}")
        # (S) source だけ違う: seed の外部 DB は H2O も持つ (source=file)。H2O を外部 DB から外して内蔵 (同一係数) にする
        check({e["name"]: e["source"] for e in rec_species(d1, at)}.get("H2O") == "file", "(S) seed record: H2O source is file")

        def drop_h2o(db):
            db.pop("H2O", None)
            return db
        # 内蔵と同一係数の H2O を外部 DB に置いた場を作る (記録は resolve-only で書き、場に属性を付ける)
        h2o_b = recNsp["H2O"]
        h2o_db = {k: h2o_b[k] for k in ("MW", "LJ_sigma", "LJ_eps_kB", "Tlo", "Tmid", "Thi", "nasa9_low", "nasa9_high")}

        def h2o_as_builtin(db):
            db["H2O"] = dict(h2o_db)
            return db
        dS0 = C.make("S_src", db_edit=h2o_as_builtin)
        rc, hS0 = C.resolve(dS0)
        recS0 = [f for f in os.listdir(dS0) if f.startswith("resolved_species_")][0]
        with open(os.path.join(dS0, recS0), "rb") as f:
            shaS0 = hashlib.sha256(f.read()).hexdigest()
        check({e["name"]: e["source"] for e in fs.load_record(os.path.join(dS0, recS0))["species"]}.get("H2O") == "file",
              "(S) source record: H2O comes from the external DB (file) with the built-in coefficients")
        dS = C.make("S", seed_res=r1, db_edit=drop_h2o)
        stamp(os.path.join(dS, "nozzle.h5"), {"species_hash": hS0, "species_record_sha256": shaS0, "species_record_file": recS0,
                                              "species_input_unverified": 0})
        shutil.copy(os.path.join(dS0, recS0), dS)
        rc, out = C.run(dS)
        recS = fs.find_record(os.path.join(dS, "res_1.h5"))[0] if rc == 0 else None
        check(rc == 0 and "species_hash matches" in out and recS is not None
              and {e["name"]: e["source"] for e in recS["species"]}.get("H2O") == "builtin",
              f"(S) source-only difference (H2O file -> builtin, same coefficients) -> same compat hash, runs (rc={rc})")
        # (S3) 段 3 の旧記録からの継続: seed の場 (外部 DB の H2O = 旧 MW) を内蔵 H2O (CEA の MW) の config で読む → 停止、H2O.MW を示す
        dS3 = C.make("S3", seed_res=r1, db_edit=drop_h2o)
        stamp(os.path.join(dS3, "nozzle.h5"), at, unverified=0)
        shutil.copy(os.path.join(d1, at["species_record_file"]), dS3)
        rc, out = C.run(dS3)
        diffs = re.findall(r"^ +(\S+): field .* vs current .*$", out, re.M)
        seed_mw = float(yaml.safe_load(open(os.path.join(a.seed_run, "species_db.yaml")))["H2O"]["MW"])
        want = ["H2O.MW"] if seed_mw != float(h2o_db["MW"]) else []
        check(rc != 0 and want and diffs == want and not os.path.exists(os.path.join(dS3, "res_1.h5")),
              f"(S3) stage-3 continuation: seed field (external H2O MW {seed_mw!r}) read with the built-in H2O "
              f"(MW {float(h2o_db['MW'])!r}) -> stops (rc={rc}), differences {diffs}")
        m = re.search(r"^ +H2O\.MW.*$", out, re.M)
        print("       " + (m.group(0).strip() if m else "(no line)"))
        rc, out = C.run(dS3, allow=True)
        check(rc != 0 and "never allowed" in out, f"(S3) the stage-3 mismatch is not allowed by FORGE_ALLOW_UNVERIFIED_SPECIES (rc={rc})")
    except _Abort:
        pass
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("ALL PASSED" if FAIL == 0 else f"FAILED ({FAIL})")
    sys.exit(1 if FAIL else 0)


class _Abort(Exception):
    pass


def rec_species(run_dir, at):
    return fs.load_record(os.path.join(run_dir, at["species_record_file"]))["species"]


if __name__ == "__main__":
    main()
