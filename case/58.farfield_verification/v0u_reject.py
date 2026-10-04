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

def bc_nocomp(b):
    # 多成分 (TP 2 種) の farfield から組成 Y/X を消す → 起動時に拒否されること (codex result 2026-10-03 M1)
    for v in b.values():
        if v.get("kind") == "farfield":
            v["floats"] = {k: x for k, x in (v.get("floats") or {}).items() if not re.fullmatch(r"[XY]\d+", str(k))}
def bc_xcomp(b):
    # 正常系: farfield の組成を Y ではなく X (モル分率) で明示 → 起動して完走すること
    for v in b.values():
        if v.get("kind") == "farfield":
            fl = v.get("floats") or {}
            ys = {k: x for k, x in fl.items() if re.fullmatch(r"Y\d+", str(k))}
            assert all(float(x) in (0.0, 1.0) for x in ys.values()), "純成分の外気だけを X に書き換える (換算なしで同値)"
            v["floats"] = {**{k: x for k, x in fl.items() if k not in ys}, **{"X" + k[1:]: x for k, x in ys.items()}}
BC_EDIT = {"組成省略": bc_nocomp, "組成 X で明示 (正常系)": bc_xcomp}

CASES = [("cell", SRC_TP, edit_cell), ("ROE", SRC_TP, edit_roe), ("凝縮", SRC_CPG, edit_cond), ("トレーサ", SRC_TP, edit_tracer),
         ("遷移", SRC_SST, edit_transition), ("軸対称", SRC_TP, edit_axisym), ("周期と共有", SRC_TP, "periodic"),
         ("組成省略", SRC_TP, "bc"), ("組成 X で明示 (正常系)", SRC_TP, "bc")]
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
POSITIVE = {"組成 X で明示 (正常系)"}   # 起動して完走すべき追加ケース

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
    if edit == "bc":
        b = yaml.safe_load(open(os.path.join(d, "bcondConfig.yaml")))
        BC_EDIT[name](b)
        yaml.safe_dump(b, open(os.path.join(d, "bcondConfig.yaml"), "w"), allow_unicode=True, sort_keys=False)
    c = yaml.safe_load(open(os.path.join(d, "solverConfig.yaml")))
    c["time"]["last"]["nStepOuter"] = 1
    c["time"]["outStepInterval"] = 1000000
    if callable(edit):
        edit(c)
    yaml.safe_dump(c, open(os.path.join(d, "solverConfig.yaml"), "w"), allow_unicode=True, sort_keys=False)
    # 種 DB 取り込み後のバイナリは属性の無い入力場を拒否するので、この起動試験では明示許可を付ける (拒否の対象は farfield の検査)
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_ALLOW_UNVERIFIED_SPECIES="1")
    p = subprocess.run([BIN], cwd=d, env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=600)
    out = p.stdout + p.stderr
    m = PAT.search(out)
    line = next((l for l in out.splitlines() if PAT.search(l)), "")
    if edit is None or name in POSITIVE:
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
