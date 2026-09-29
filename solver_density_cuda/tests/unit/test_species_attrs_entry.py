#!/usr/bin/env python3
"""引き継ぎ・種変換の入口での化学種属性の扱い (plans/active/thermophysics-solver-owned-species-db.md §5.1 #3b, §6 V1 (a)(d))。

  python3 solver_density_cuda/tests/unit/test_species_attrs_entry.py --forge BIN     (GPU 不要; BIN は --resolve-species 対応)

合成の 2 種 run (N2, H2O; test_convert_species_field_fail.py と同じ config / DB) を一時ディレクトリに作り、ツールを子プロセスで呼ぶ。
SRC の res には「ソルバが書いたのと同じ」属性を付ける (記録は SRC run の `forge --resolve-species` で作る試験用の記録)。

  (a)  restart_field / interp_field: 同一物性の宛先 → 継承 (species_input_unverified=0、記録が DST の隣に複製され find_record で引ける)
  (a)  restart_field / interp_field: 宛先の外部 DB で N2 nasa9_low[2] +0.001 → 宛先の --resolve-species と不一致で停止し
       N2.nasa9_low[2] を表示、DST の保存量と属性は書き換わらない (d: 起動前の宛先解決)
  (a)  保存場だけを別ディレクトリへ (記録なし) → 照合不能で停止; --force-species で写すが属性は付かない
  (a)  不一致 + --force-species → 写すが属性は付かない
  (i)  SRC 属性なし / species_input_unverified=1 (宛先 TP) → 既定で停止し DST の保存量・属性は不変 (#3c: ソルバと同じ規約)、
       案内は「IC を属性を付ける処理で作り直す / その実行だけ FORGE_ALLOW_UNVERIFIED_SPECIES=1」の 2 通り。
       FORGE_ALLOW_UNVERIFIED_SPECIES=1 / --force-species → 写し、DST の既存属性を消す (species_input_unverified も付けない)
  (n)  宛先に solverConfig.yaml が無い (CPG か TP か判定できない) + 未検証 SRC → 既定で停止
  (cpg) 宛先が CPG (thermalMethod 0) + 属性なし SRC → 属性の対象外なので既定で通る (属性なし)
  (x)  interp_field: 宛先が H2O を内蔵 DB で持つ (同一係数) → 記録で照合して継承 (旧: 設定の署名で「照合不能」拒否)
  (o)  --resolve-species を持たない旧バイナリを FORGE_BIN にしても起動しない。宛先を解決できないので既定で停止 (DST 不変)、
       FORGE_ALLOW_UNVERIFIED_SPECIES=1 で属性なしで写す (#3c 以前は既定で警告して通していた)。
       旧バイナリが手元に無ければ --resolve-species の文字列を含まない偽バイナリ (起動されたら印を残すスクリプト) で代える
  (c)  convert_species_field: 検証済み入力 → 変換後に宛先のハッシュ (記録は宛先 run に書かれ完全性一致)、未検証入力 → 既定で停止
       (DST 不変)、FORGE_ALLOW_UNVERIFIED_SPECIES=1 で属性なし、
       場の記録と SRC run の設定が違う (N2 low[2]) → 書き込み前に拒否
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, hashlib, os, shutil, subprocess, sys, tempfile

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TOOLS)
sys.path.insert(0, HERE)
import forge_species as fs  # noqa: E402
from test_convert_species_field_fail import DB, N, make_run, write_h5  # noqa: E402

FAIL = 0
# 現行ソルバが受け付ける最小の TP config (--resolve-species は physProp だけを使う; 種・DB・datum は変換器試験と同じ)
CFG = ('mesh: {discretization: "node", meshFileName: "in.h5", valueFileName: "in.h5"}\ngpu: 1\nsolver: "SLAU"\n'
       'physProp: {thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1000.0, gamma: 1.4,\n'
       '           species: [N2, H2O], speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15}\n'
       'time: {unsteady: 0, dualTime: 0, last: {nStepOuter: 1}, deltaT: {control: 1, dt: 1e-8, cfl: 1.0, cfl_pseudo: 1.0, '
       'dt_min: 1e-9, dt_max: 1e-3}, outStepStart: 0, outStepInterval: 1, timeIntegration: 11, nStepInner: 5}\n'
       'space: {convMethod: 1, limiter: 2}\nturbulence: {model: "none"}\ninitial: "uniform_p101325_u10"\noutput: {level: 1}\n')


def check(ok, what, detail=""):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what + (f"\n       {detail}" if (detail and not ok) else ""), flush=True)
    if not ok:
        FAIL += 1


def attrs(h5):
    with h5py.File(h5, "r") as f:
        return {k: (f.attrs[k].decode() if isinstance(f.attrs[k], bytes) else f.attrs[k]) for k in fs.SPECIES_ATTRS if k in f.attrs}


def values(h5):
    with h5py.File(h5, "r") as f:
        return {k: np.asarray(f["VALUE/" + k]) for k in f["VALUE"]}


def same_values(a, b):
    return a.keys() == b.keys() and all(np.array_equal(a[k], b[k]) for k in a)


ALLOW = {"FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}


def tool(name, *args, env_extra=None, drop_env=("FORGE_ALLOW_UNVERIFIED_SPECIES",)):
    env = dict(os.environ)
    for k in drop_env:
        env.pop(k, None)
    env.update(env_extra or {})
    p = subprocess.run([sys.executable, os.path.join(TOOLS, name), *args], capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def bump_db(d, name="N2", delta=0.001):
    db = yaml.safe_load(open(os.path.join(d, "species_db.yaml")))
    db[name]["nasa9_low"][2] = float(db[name]["nasa9_low"][2]) + delta
    yaml.safe_dump(db, open(os.path.join(d, "species_db.yaml"), "w"), sort_keys=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    if not a.forge or not fs.forge_supports_resolve(a.forge):
        raise SystemExit("--forge BIN (or FORGE_BIN) with --resolve-species is required")
    os.environ["FORGE_BIN"] = os.path.abspath(a.forge)
    root = tempfile.mkdtemp(prefix="forge_species_attrs_")
    ro0 = np.full(N, 0.8); Y0 = np.array([np.full(N, 0.95), np.full(N, 0.05)]); T0 = np.linspace(300.0, 1500.0, N)
    ro1 = np.full(N, 0.6)

    def run_dir(name, bump=False, h2o_builtin=False, h5="in.h5", ro=ro1):
        d = os.path.join(root, name)
        make_run(d)
        open(os.path.join(d, "solverConfig.yaml"), "w").write(CFG)   # 変換器試験の config は廃止キーを含むので現行形式に置き換える
        if bump:
            bump_db(d)
        if h2o_builtin:
            db = yaml.safe_load(open(os.path.join(d, "species_db.yaml"))); db.pop("H2O")
            yaml.safe_dump(db, open(os.path.join(d, "species_db.yaml"), "w"), sort_keys=False)
        write_h5(os.path.join(d, h5), ro, Y0, T0)
        return d

    try:
        # --- SRC: 属性付きの保存場 (記録は SRC run の resolve-only 記録) ---
        src = run_dir("src", h5="res_10.h5", ro=ro0)
        r = fs.resolve_species(src)
        srcres = os.path.join(src, "res_10.h5")
        src_attrs = {"species_hash": r["hash"], "species_record_sha256": r["record"]["integrity"],
                     "species_record_file": r["record_file"], "species_input_unverified": 0}
        fs.write_species_attrs(srcres, src_attrs)
        st = fs.source_species_state(srcres)
        check(st["state"] == "verified", f"fixture: SRC res verified by its record ({st['why'][:80]})")

        # (a) 同一物性 → 継承
        for t in ("restart_field.py", "interp_field.py"):
            d = run_dir("same_" + t[:-3])
            rc, out = tool(t, srcres, os.path.join(d, "in.h5"))
            at = attrs(os.path.join(d, "in.h5"))
            rec, why = fs.find_record(os.path.join(d, "in.h5"))
            check(rc == 0 and at.get("species_hash") == r["hash"] and at.get("species_input_unverified") == 0
                  and at.get("species_record_sha256") == r["record"]["integrity"] and rec is not None
                  and "forge --resolve-species" in out,
                  f"(a) {t}: same species -> inherited (hash {at.get('species_hash', '')[:16]}, record next to DST: {why[:60]})",
                  out[-600:])

        # (a)(d) 宛先の外部 DB で N2 nasa9_low[2] +0.001 → 宛先解決と不一致で停止、DST 不変
        for t in ("restart_field.py", "interp_field.py"):
            d = run_dir("bump_" + t[:-3], bump=True)
            dst = os.path.join(d, "in.h5")
            fs.write_species_attrs(dst, {"species_hash": "0" * 64, "species_record_sha256": "1" * 64,
                                         "species_record_file": "x.yaml", "species_input_unverified": 0})
            before, at0 = values(dst), attrs(dst)
            rc, out = tool(t, srcres, dst)
            check(rc != 0 and "N2.nasa9_low[2]" in out and "REFUSED" in out and "forge --resolve-species" in out
                  and same_values(before, values(dst)) and attrs(dst) == at0,
                  f"(a)(d) {t}: destination DB N2 nasa9_low[2] +0.001 -> stops before writing, shows N2.nasa9_low[2] (rc={rc})",
                  out[-800:])
            rc, out = tool(t, srcres, dst, "--force-species")
            check(rc == 0 and attrs(dst) == {} and not same_values(before, values(dst)),
                  f"(a) {t} --force-species on a mismatch -> copies, no species attributes (rc={rc})", out[-600:])

        # (a) 保存場だけコピー (記録なし) → 照合不能で停止
        lone = os.path.join(root, "lone"); os.makedirs(lone); shutil.copy(srcres, lone)
        d = run_dir("lone_dst")
        rc, out = tool("restart_field.py", os.path.join(lone, "res_10.h5"), os.path.join(d, "in.h5"))
        check(rc != 0 and "UNVERIFIABLE" in out, f"(a) field copied without its record -> UNVERIFIABLE, stops (rc={rc})", out[-600:])
        rc, out = tool("restart_field.py", os.path.join(lone, "res_10.h5"), os.path.join(d, "in.h5"), "--force-species")
        check(rc == 0 and attrs(os.path.join(d, "in.h5")) == {}, f"(a) same with --force-species -> copies without attributes (rc={rc})")

        # (i) 未検証の SRC → 既定で停止 (DST 不変); 許可 (環境変数 / --force-species) で写し、DST の既存属性を消す
        for unv, label in ((None, "no attributes"), (1, "species_input_unverified=1")):
            s2 = run_dir(f"src_unv{unv}", h5="res_10.h5", ro=ro0)
            if unv is not None:
                fs.write_species_attrs(os.path.join(s2, "res_10.h5"), dict(src_attrs, species_input_unverified=1))
                shutil.copy(os.path.join(src, r["record_file"]), s2)
            for t in ("restart_field.py", "interp_field.py"):
                d = run_dir(f"unv{unv}_" + t[:-3])
                dst = os.path.join(d, "in.h5")
                fs.write_species_attrs(dst, src_attrs)      # 古い属性が残っていたとする
                before, at0 = values(dst), attrs(dst)
                rc, out = tool(t, os.path.join(s2, "res_10.h5"), dst)
                check(rc != 0 and "REFUSED" in out and "UNVERIFIED" in out and "regenerate the initial field" in out
                      and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out and same_values(before, values(dst)) and attrs(dst) == at0,
                      f"(i) {t}: SRC {label} -> stops by default, DST values/attributes unchanged, 2-way guidance (rc={rc})",
                      out[-800:])
                rc, out = tool(t, os.path.join(s2, "res_10.h5"), dst, env_extra=ALLOW)
                check(rc == 0 and attrs(dst) == {} and not same_values(before, values(dst))
                      and "Allowed for this invocation by FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out,
                      f"(i) {t}: SRC {label} + FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> copies, DST attributes removed, "
                      f"no species_input_unverified (rc={rc})", out[-600:])
                fs.write_species_attrs(dst, src_attrs)
                rc, out = tool(t, os.path.join(s2, "res_10.h5"), dst, "--force-species")
                check(rc == 0 and attrs(dst) == {}, f"(i) {t}: SRC {label} + --force-species -> copies without attributes (rc={rc})",
                      out[-600:])

        # (n) 宛先に solverConfig.yaml が無い → CPG か TP か判定できないので既定で停止
        s2 = os.path.join(root, "src_unvNone", "res_10.h5")
        nocfg = os.path.join(root, "nocfg"); os.makedirs(nocfg)
        write_h5(os.path.join(nocfg, "in.h5"), ro1, Y0, T0)
        before = values(os.path.join(nocfg, "in.h5"))
        rc, out = tool("restart_field.py", s2, os.path.join(nocfg, "in.h5"))
        check(rc != 0 and "no solverConfig.yaml" in out and same_values(before, values(os.path.join(nocfg, "in.h5"))),
              f"(n) restart_field: unverified SRC, destination without solverConfig.yaml -> stops (rc={rc})", out[-600:])

        # (cpg) 宛先が CPG → 属性の対象外。属性なしの SRC から既定で通る
        dc = run_dir("cpg_dst")
        open(os.path.join(dc, "solverConfig.yaml"), "w").write(CFG.replace("thermalMethod: 2", "thermalMethod: 0"))
        rc, out = tool("restart_field.py", s2, os.path.join(dc, "in.h5"))
        check(rc == 0 and attrs(os.path.join(dc, "in.h5")) == {},
              f"(cpg) restart_field: SRC without attributes -> CPG destination passes by default (rc={rc})", out[-600:])

        # (x) interp_field: 宛先の H2O が内蔵 (同一係数なら同じ互換性ハッシュ) → 記録で照合して継承
        d = run_dir("h2o_builtin", h2o_builtin=True)
        rd = fs.resolve_species(d, inplace=False)
        rc, out = tool("interp_field.py", srcres, os.path.join(d, "in.h5"))
        if rd["hash"] == r["hash"]:
            check(rc == 0 and attrs(os.path.join(d, "in.h5")).get("species_hash") == r["hash"],
                  f"(x) interp_field: H2O file vs built-in with the same coefficients -> inherited by record (rc={rc})", out[-600:])
        else:
            diff = fs.compare_signatures(fs.signature_from_record(r["record"]), fs.signature_from_record(rd["record"]))
            check(rc != 0 and "REFUSED" in out,
                  f"(x) interp_field: test DB H2O differs from the built-in ({diff[:2]}) -> refused by record (rc={rc})", out[-600:])

        # (o) 旧バイナリ (--resolve-species なし) を FORGE_BIN にしても起動しない。宛先を解決できないので既定で停止
        old = os.path.join(REPO, "solver_density_cuda", "build", "forge")
        marker = os.path.join(root, "old_binary_launched")
        if not (os.path.exists(old) and not fs.forge_supports_resolve(old)):
            old = os.path.join(root, "fake_old_forge")     # 起動されたら印を残す偽の旧バイナリ
            with open(old, "w") as f:
                f.write(f"#!/bin/sh\ntouch '{marker}'\n")
            os.chmod(old, 0o755)
            print(f"[INFO] (o) no old binary at solver_density_cuda/build/forge; using a fake old binary {old}")
        d = run_dir("oldbin")
        dst = os.path.join(d, "in.h5")
        fs.write_species_attrs(dst, src_attrs)
        before, at0 = values(dst), attrs(dst)
        rc, out = tool("restart_field.py", srcres, dst, env_extra={"FORGE_BIN": old})
        check(rc != 0 and "REFUSED" in out and "旧バイナリ" in out and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out
              and same_values(before, values(dst)) and attrs(dst) == at0
              and not os.path.exists(os.path.join(d, "res_1.h5")) and not os.path.exists(marker),
              f"(o) old binary as FORGE_BIN -> not launched; destination unresolvable -> stops by default, DST unchanged (rc={rc})",
              out[-600:])
        rc, out = tool("restart_field.py", srcres, dst, env_extra=dict(ALLOW, FORGE_BIN=old))
        check(rc == 0 and attrs(dst) == {} and "旧バイナリ" in out and not os.path.exists(marker),
              f"(o) same with FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> not launched; copied without attributes + warning (rc={rc})",
              out[-600:])

        # (c) convert_species_field: 検証済み入力 → 宛先のハッシュを付ける / 未検証入力 → 属性なし
        conv = os.path.join(TOOLS, "convert_species_field.py")

        def convert(srcdir, srch5, dstdir, *extra, env_extra=None):
            env = dict(os.environ); env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None); env.update(env_extra or {})
            p = subprocess.run([sys.executable, conv, srch5, os.path.join(dstdir, "in.h5"), "--meta", os.path.join(dstdir, "species_meta.yaml"),
                                "--src-meta", os.path.join(srcdir, "species_meta.yaml"), *extra], capture_output=True, text=True, env=env)
            return p.returncode, p.stdout + p.stderr
        d = run_dir("conv_ok")
        rc, out = convert(src, srcres, d)
        at = attrs(os.path.join(d, "in.h5"))
        rec, why = fs.find_record(os.path.join(d, "in.h5"))
        check(rc == 0 and at.get("species_input_unverified") == 0 and rec is not None and rec["compat_recomputed"] == at.get("species_hash")
              and os.path.dirname(rec["path"]) == d,
              f"(c) convert: verified input -> destination hash {at.get('species_hash', '')[:16]} with its record in the destination run (rc={rc})",
              out[-800:])
        s3 = run_dir("conv_src_unv", h5="res_10.h5", ro=ro0)
        d = run_dir("conv_unv")
        fs.write_species_attrs(os.path.join(d, "in.h5"), src_attrs)
        before, at0 = values(os.path.join(d, "in.h5")), attrs(os.path.join(d, "in.h5"))
        rc, out = convert(s3, os.path.join(s3, "res_10.h5"), d)
        check(rc != 0 and "REFUSED" in out and "FORGE_ALLOW_UNVERIFIED_SPECIES=1" in out
              and same_values(before, values(os.path.join(d, "in.h5"))) and attrs(os.path.join(d, "in.h5")) == at0,
              f"(c) convert: unverified input -> stops by default, DST unchanged (rc={rc})", out[-600:])
        rc, out = convert(s3, os.path.join(s3, "res_10.h5"), d, env_extra=ALLOW)
        check(rc == 0 and attrs(os.path.join(d, "in.h5")) == {},
              f"(c) convert: unverified input + FORGE_ALLOW_UNVERIFIED_SPECIES=1 -> output unverified (no attributes) (rc={rc})",
              out[-600:])
        # 入力の記録と SRC run の設定が違う (場を作った物性と変換器が読む物性が違う) → 拒否
        s4 = run_dir("conv_src_mismatch", h5="res_10.h5", ro=ro0)
        shutil.copy(os.path.join(src, r["record_file"]), s4)
        fs.write_species_attrs(os.path.join(s4, "res_10.h5"), src_attrs)
        bump_db(s4)
        d = run_dir("conv_mis_dst")
        before = values(os.path.join(d, "in.h5"))
        rc, out = convert(s4, os.path.join(s4, "res_10.h5"), d)
        check(rc != 0 and "N2.nasa9_low[2]" in out and same_values(before, values(os.path.join(d, "in.h5"))),
              f"(c) convert: field record != source run config (N2 low[2]) -> refused before writing (rc={rc})", out[-800:])
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print("ALL PASSED" if FAIL == 0 else f"FAILED ({FAIL})")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
