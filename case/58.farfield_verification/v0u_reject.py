#!/usr/bin/env python3
"""V0u (v) 非対応構成の起動拒否 (plan boundary-node-farfield-characteristic §6 V0u (v))。
V1 の run (TP: run_0005_v1b_fix、SST: run_0004_v1d) の入力を一時ディレクトリに複製し、構成を 1 つずつ変えて forge を起動する。
合格 = 各構成で forge が非 0 終了し、stderr/stdout に farfield の拒否メッセージ (「kind farfield」または「周期境界と共有」) が出る。
  python3 v0u_reject.py [--bin FORGE] [--keep]
"""
import os, re, shutil, subprocess, sys, tempfile
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = sys.argv[sys.argv.index("--bin") + 1] if "--bin" in sys.argv else os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff/forge")
SRC_CPG, SRC_TP, SRC_SST = os.path.join(HERE, "run_0001_v1a"), os.path.join(HERE, "run_0005_v1b_fix"), os.path.join(HERE, "run_0004_v1d")
FILES = ("solverConfig.yaml", "bcondConfig.yaml", "box.h5", "probe.yaml", "species_db.yaml", "species_meta.yaml")


def edit_cell(c):
    c["mesh"]["discretization"] = "cell"
def edit_roe(c):
    c["solver"] = "ROE"
def edit_cond(c):
    # CPG 担体 + N2 選択凝縮 (空気凝縮の生産形) (蒸気質量分率を与える形)。凝縮の設定自体は有効で、farfield の検査で落ちることを見る
    c["condensation"] = {"condensation": 1, "nCondSpecies": 1, "condModel": 0, "condVaporMassFraction": 0.7671}
def edit_tracer(c):
    c["physProp"]["tracer"] = "exhaust"
def edit_transition(c):
    c["turbulence"]["transition"] = "lm2009"
def edit_axisym(c):
    c["mesh"]["isAxisymmetric"] = 1

CASES = [("cell", SRC_TP, edit_cell), ("ROE", SRC_TP, edit_roe), ("凝縮", SRC_CPG, edit_cond), ("トレーサ", SRC_TP, edit_tracer),
         ("遷移", SRC_SST, edit_transition), ("軸対称", SRC_TP, edit_axisym), ("周期と共有", SRC_TP, "periodic")]
def periodic_prep(d):
    """zmin/zmax を周期にして変換し直す (xmin..ymax は farfield のまま → 周期と節点を共有する)"""
    b = yaml.safe_load(open(os.path.join(d, "bcondConfig.yaml")))
    b["zmin"] = {"physID": 5, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 6}, "floats": {"dx": 0.0, "dy": 0.0, "dz": 0.5}}
    b["zmax"] = {"physID": 6, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 5}, "floats": {"dx": 0.0, "dy": 0.0, "dz": -0.5}}
    yaml.safe_dump(b, open(os.path.join(d, "bcondConfig.yaml"), "w"), sort_keys=False)
    os.remove(os.path.join(d, "box.h5"))
    shutil.copy(os.path.join(SRC_SST, "box.msh"), d)
    subprocess.run([os.path.join(os.path.dirname(BIN), "convertGmshToForge"), "box.msh", "box.h5"], cwd=d,
                   stdout=open(os.path.join(d, "convert.log"), "w"), stderr=subprocess.STDOUT)

PAT = re.compile(r"kind farfield|周期境界と共有")

fails = 0
tmp = tempfile.mkdtemp(prefix="v0u_reject_", dir=HERE)
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
for name, src, edit in CASES + [("対照 (変更なし)", SRC_TP, None)]:
    if ONLY and name != ONLY:
        continue
    d = os.path.join(tmp, re.sub(r"\W", "_", name))
    os.makedirs(d)
    for f in FILES:
        if os.path.exists(os.path.join(src, f)):
            shutil.copy(os.path.abspath(os.path.join(src, f)), os.path.join(d, f))
    if edit == "periodic":
        periodic_prep(d)
    c = yaml.safe_load(open(os.path.join(d, "solverConfig.yaml")))
    c["time"]["last"]["nStepOuter"] = 1
    c["time"]["outStepInterval"] = 1000000
    if callable(edit):
        edit(c)
    yaml.safe_dump(c, open(os.path.join(d, "solverConfig.yaml"), "w"), allow_unicode=True, sort_keys=False)
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128")
    p = subprocess.run([BIN], cwd=d, env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=600)
    out = p.stdout + p.stderr
    m = PAT.search(out)
    line = next((l for l in out.splitlines() if PAT.search(l)), "")
    if edit is None:
        ok = p.returncode == 0 and not m
        print(f"{name}: 終了コード {p.returncode} → {'PASS' if ok else 'FAIL'} (対照は起動して完走すること)")
    else:
        ok = p.returncode != 0 and bool(m)
        err = line if m else next((l for l in out.splitlines() if "rror" in l), out.strip().splitlines()[-1] if out.strip() else "")
        print(f"{name}: 終了コード {p.returncode}、{err.strip()[:160]} → {'PASS' if ok else 'FAIL'}")
    fails += 0 if ok else 1
if "--keep" not in sys.argv:
    shutil.rmtree(tmp)
print(f"VERDICT: {'PASS' if fails == 0 else 'FAIL'} ({fails} 件不合格)")
sys.exit(1 if fails else 0)
