#!/usr/bin/env python3
"""plan convection-slau-wall-normal-chi-default §6 B0 (1 step の面流束のビット不変、判定は測る前に固定)。

  (a) node + SLAU で chi 省略 (auto → 1) vs 明示 1  → massflux ビット同一
  (b) node + SLAU で明示 0 (現行) vs 実装前 commit 22976398 のバイナリ (明示 0) → massflux ビット同一
      現行は 2026-09-27 に scalarGradient の node 既定を lsq に変えたので、条件を揃えるため現行側は scalarGradient: gg を明記する
  (c) cell・非 SLAU (ROE)・nodeWallDirichlet 0 で省略 → 起動エコーが 0 (auto) (設定解決の確認のみ)

入力は case/48 の 1 step 用準備物 (node・SLAU・nodeWallDirichlet 1、壁あり)。AWS で実行:
  python3 chi_default_b0.py --prep ~/sglsq/s1h/s1_prep/case48 --new ~/sglsq/forge_2fa3826c --old ~/sglsq/forge_22976398 --scratch ~/sglsq/chib0
"""
import argparse
import os
import shutil
import subprocess

import numpy as np
import yaml


def mk(prep, d, edit):
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(prep, d)
    for f in os.listdir(d):
        if f.startswith("res_") or f.endswith(".log") or f.startswith("massflux"):
            os.remove(os.path.join(d, f))
    p = os.path.join(d, "solverConfig.yaml")
    c = yaml.safe_load(open(p))
    c["time"]["last"]["nStepOuter"] = 1
    edit(c)
    yaml.safe_dump(c, open(p, "w"), sort_keys=False, default_flow_style=None)


def run(d, fbin, dump=True, timeout=600):
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_CUDA_BLOCKSIZE_SMALL="128")
    if dump:
        env["FORGE_DUMP_MASSFLUX"] = os.path.join(d, "massflux.bin")
    r = subprocess.run([fbin], cwd=d, env=env, capture_output=True, text=True, timeout=timeout)
    open(os.path.join(d, "forge_stdout.txt"), "w").write(r.stdout + r.stderr)
    echo = [l for l in (r.stdout + r.stderr).splitlines() if "slauWallNormalChi' effective" in l]
    return r.returncode, echo


def cmp(a, b):
    x = np.fromfile(os.path.join(a, "massflux.bin"), dtype=np.float32)
    y = np.fromfile(os.path.join(b, "massflux.bin"), dtype=np.float32)
    if x.shape != y.shape:
        return None
    return int(np.count_nonzero(x.view(np.uint32) != y.view(np.uint32))), x.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep", required=True); ap.add_argument("--new", required=True); ap.add_argument("--old", required=True)
    ap.add_argument("--scratch", required=True); ap.add_argument("--out", default="CHI_DEFAULT_B0.txt")
    a = ap.parse_args()
    os.makedirs(a.scratch, exist_ok=True)
    S = lambda n: os.path.join(a.scratch, n)
    def chi(v):
        def f(c):
            c.setdefault("space", {})
            if v is None:
                c["space"].pop("slauWallNormalChi", None)
            else:
                c["space"]["slauWallNormalChi"] = v
        return f
    def both(*fs):
        def f(c):
            for g in fs:
                g(c)
        return f
    sg_gg = lambda c: c["mesh"].__setitem__("scalarGradient", "gg")
    L = ["# chi 既定化 B0 (plan convection-slau-wall-normal-chi-default §6)", f"new {a.new}", f"old {a.old}", f"prep {a.prep}", ""]
    ok = True
    # (a)
    mk(a.prep, S("a_omit"), chi(None)); rc1, e1 = run(S("a_omit"), a.new)
    mk(a.prep, S("a_exp1"), chi(1)); rc2, e2 = run(S("a_exp1"), a.new)
    r = cmp(S("a_omit"), S("a_exp1")) if rc1 == 0 and rc2 == 0 else None
    good = r is not None and r[0] == 0; ok &= good
    L.append(f"(a) 省略 vs 明示 1: rc {rc1}/{rc2}、エコー {e1} / {e2}、massflux 不一致 {r} → {'PASS' if good else 'FAIL'}")
    # (b)
    mk(a.prep, S("b_new0"), both(chi(0), sg_gg)); rc1, e1 = run(S("b_new0"), a.new)
    mk(a.prep, S("b_old0"), chi(0)); rc2, e2 = run(S("b_old0"), a.old)
    r = cmp(S("b_new0"), S("b_old0")) if rc1 == 0 and rc2 == 0 else None
    good = r is not None and r[0] == 0; ok &= good
    L.append(f"(b) 明示 0 (現行 + scalarGradient gg) vs 実装前 22976398: rc {rc1}/{rc2}、エコー {e1} / {e2}、massflux 不一致 {r} → {'PASS' if good else 'FAIL'}")
    # (c) 設定解決のみ (実行が途中で止まってもエコーだけ見る)
    for tag, ed in (("cell", lambda c: c["mesh"].__setitem__("discretization", "cell")),
                    ("roe", lambda c: c.__setitem__("solver", "ROE")),
                    ("nwd0", lambda c: c["mesh"].__setitem__("nodeWallDirichlet", 0))):
        mk(a.prep, S("c_" + tag), both(chi(None), ed))
        try:
            rc, e = run(S("c_" + tag), a.new, dump=False, timeout=120)
        except subprocess.TimeoutExpired:
            rc, e = "timeout", []
            e = [l for l in open(os.path.join(S("c_" + tag), "forge_stdout.txt")).read().splitlines() if "effective" in l] if os.path.exists(os.path.join(S("c_" + tag), "forge_stdout.txt")) else []
        good = bool(e) and all(": 0 (auto" in x for x in e); ok &= good
        L.append(f"(c) {tag} で省略: rc {rc}、エコー {e} → {'PASS' if good else 'FAIL'}")
    L.append(f"\nVERDICT B0: {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(L) + "\n"
    open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
