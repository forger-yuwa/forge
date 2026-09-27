#!/usr/bin/env python3
"""輸送物性の表引き (#5t2-3) の性能測定。plans/active/thermophysics-solver-owned-species-db.md §5.1 #5t2-3 の事前固定の測り方:

  (新物性時間 − 旧物性時間)/旧 step 時間 ≤ 0.10   (実種 n = 5・12 で判定、32 は記録のみ)
  同じ入力・GPU・ビルド条件で、暖機後 300 step × 5 反復。反復ごとの比の最大値で判定。

  python3 solver_density_cuda/tests/unit/bench_transport_table.py --forge BIN [--reps 5] [--steps 300] [--warm 20] [--keep]

各構成を 3 通りで回す (同じバイナリ・同じ初期場・FORGE_PROFILE=1・FORGE_CUDA_BLOCKSIZE=256):
  legacy : physProp.transport なし (viscMethod 2 の従来の Wilke 経路) = 「旧」
  table  : physProp.transport あり、表引き (既定) = 「新」
  double : physProp.transport あり、FORGE_TRANSPORT_TABLE=0 (段 2 の double 評価; 記録のみ)
暖機の除外: 同じ入力で nStep = warm と warm + steps の 2 本を回し、FORGE_PROFILE の区間合計の差を steps で割る
(1 本目の step に入るモジュール読み込み・キャッシュの影響を差し引く)。step 時間は FORGE_PROFILE の share から逆算
(最大 share の区間の total_ms / share)。
構成: n5 = case/44 va3 の lump 記法 MIXDRY + H2O (実種 N2・O2・AR・CO2 は cea、H2O は custom:h2o_iapws_cea_v1) を run_0509 の res_24000 から、
      n12・n32 = test_transport_gpu.py の real12・real32 の種構成を、res_24000 の流れ場に固定組成の ρY を載せた場から
      (組成は空気に近い lump を主にし他の輸送種は各 1 %; 全実種が 0 でないので全組が評価される)。
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_transport_gpu as tg  # noqa: E402
from transport_reference import Reference  # noqa: E402

SEED_REL = tg.SEED_REL
MIXDRY = {"N2": 0.7088730183520131, "O2": 0.23037573447245455, "AR": 0.008503874694843204, "CO2": 0.05224737248068922}


def configs():
    c = {}
    c["n5"] = dict(species=[{"name": "MIXDRY", "lump": MIXDRY, "basis": "mole"}, "H2O"],
                   transport={"N2": "cea", "O2": "cea", "AR": "cea", "CO2": "cea", "H2O": "custom:h2o_iapws_cea_v1"},
                   db=None, X=None)
    db12 = {f"XS{i:02d}": tg.ext_species(i, "fit" if i % 2 else "kinetic") for i in range(1, 7)}
    sp12 = [{"name": "MIXA", "lump": {"N2": 0.70, "O2": 0.22, "AR": 0.01, "CO2": 0.07}, "basis": "mole"},
            {"name": "MIXB", "lump": {"N2": 0.4, "H2O": 0.35, "He": 0.25}, "basis": "mass"}] + list(db12)
    tr12 = {"N2": "cea", "O2": "kinetic", "AR": "cea", "CO2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1", "He": "cea"}
    tr12.update({k: ("fit" if int(k[2:]) % 2 else "kinetic") for k in db12})
    # 組成: 空気に近い lump を主にし、他の種は各 1 % (全実種が 0 でない = 全組が評価される; va3 の流れ場のエネルギーと大きく食い違わない)
    c["n12"] = dict(species=sp12, transport=tr12, db=db12, X=[0.93] + [0.01] * 7)
    db32 = {f"XS{i:02d}": tg.ext_species(i, "fit" if i % 2 else "kinetic") for i in range(1, 27)}
    lumpA = {"N2": 0.5, "O2": 0.2, "AR": 0.05, "CO2": 0.05}
    lumpA.update({f"XS{i:02d}": 0.025 for i in range(1, 9)})
    lumpB = {"H2O": 0.3, "He": 0.2}
    lumpB.update({f"XS{i:02d}": 0.05 for i in range(9, 19)})
    sp32 = [{"name": "LA", "lump": lumpA, "basis": "mole"}, {"name": "LB", "lump": lumpB, "basis": "mass"}, "N2"] + \
           [f"XS{i:02d}" for i in range(19, 27)]
    tr32 = {"N2": "cea", "O2": "cea", "AR": "kinetic", "CO2": "cea", "H2O": "custom:h2o_iapws_cea_v1", "He": "kinetic"}
    tr32.update({k: ("fit" if int(k[2:]) % 2 else "kinetic") for k in db32})
    c["n32"] = dict(species=sp32, transport=tr32, db=db32, X=[1.0 - 0.01 * (len(sp32) - 1)] + [0.01] * (len(sp32) - 1))
    return c


def make_value_file(seed, dst, ref, X):
    """res_24000 の流れ場に、固定組成 (輸送種のモル分率 X) の ρY を載せた初期場。X が None なら res_24000 をそのまま。"""
    src = os.path.join(seed, "res_24000.h5")
    shutil.copy(src, dst)
    if X is None:
        return
    w = [x * m for x, m in zip(X, ref.mw)]
    Y = [v / sum(w) for v in w]
    with h5py.File(dst, "r+") as f:
        g = f["VALUE"]
        ro = np.asarray(g["ro"])
        for k in [k for k in g.keys() if k.startswith("roY") or re.fullmatch(r"Y\d+", k)]:
            del g[k]
        for s, y in enumerate(Y):
            g.create_dataset(f"roY{s}", data=(ro * y).astype(ro.dtype))


def parse_profile(log):
    sec = {}
    for m in re.finditer(r"^(\S+)\s+total_ms=([\d.]+)\s+avg_ms=([\d.]+)\s+share=([\d.]+)\s*%\s+calls=(\d+)", log, re.M):
        sec[m.group(1)] = (float(m.group(2)), float(m.group(4)), int(m.group(5)))
    steps = int(re.search(r"Profiled steps: (\d+)", log).group(1))
    big = max(sec.values(), key=lambda v: v[1])
    step_total = big[0] / (big[1] / 100.0)
    return {"gas": sec.get("gas_properties", (0.0, 0.0, 0))[0], "step": step_total, "steps": steps}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", required=True)
    ap.add_argument("--seed-run", default="")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--warm", type=int, default=20)
    ap.add_argument("--configs", default="n5,n12,n32")
    ap.add_argument("--variants", default="legacy,table,double")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    a.forge = os.path.abspath(a.forge)
    if not a.seed_run:
        for c in (os.path.join(tg.REPO, SEED_REL), os.path.join(os.path.dirname(tg.REPO), "forge", SEED_REL)):
            if os.path.exists(os.path.join(c, "res_24000.h5")):
                a.seed_run = c
                break
    root = tempfile.mkdtemp(prefix="transport_bench_")
    print(f"work dir: {root}\nseed: {a.seed_run}\nforge: {a.forge}", flush=True)
    cfg0 = yaml.safe_load(open(os.path.join(a.seed_run, "solverConfig.yaml")))
    bc0 = yaml.safe_load(open(os.path.join(a.seed_run, "bcondConfig.yaml")))
    C = configs()
    results = {}
    for name in a.configs.split(","):
        c = C[name]
        ref = Reference(c["species"], c["transport"], c["db"])
        d0 = os.path.join(root, name)
        os.makedirs(d0)
        make_value_file(a.seed_run, os.path.join(d0, "init.h5"), ref, c["X"])
        Xbc = c["X"] if c["X"] is not None else [bc0["inlet"]["floats"]["X0"], bc0["inlet"]["floats"]["X1"]]
        for var in a.variants.split(","):
            for rep in range(a.reps):
                for nst in (a.warm, a.warm + a.steps):
                    d = os.path.join(d0, f"{var}_r{rep}_n{nst}")
                    os.makedirs(d)
                    for fn in ("nozzle.h5", "probe.yaml"):
                        os.symlink(os.path.join(a.seed_run, fn), os.path.join(d, fn))
                    os.symlink(os.path.join(d0, "init.h5"), os.path.join(d, "init.h5"))
                    cfg = json.loads(json.dumps(cfg0))
                    ph = cfg["physProp"]
                    ph["viscMethod"] = 2
                    ph["visc"], ph["thermCond"] = 1.8e-5, 0.026
                    ph["species"] = c["species"]
                    ph.pop("speciesDBFile", None)
                    ph.pop("transport", None)
                    if c["db"]:
                        with open(os.path.join(d, "db.yaml"), "w") as f:
                            yaml.safe_dump(c["db"], f, sort_keys=False)
                        ph["speciesDBFile"] = "db.yaml"
                    if var != "legacy":
                        ph["transport"] = c["transport"]
                    cfg["mesh"]["valueFileName"] = "init.h5"
                    cfg["time"]["last"]["nStepOuter"] = nst
                    cfg["time"]["outStepInterval"] = 1000000
                    cfg["output"] = {"level": 1}
                    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
                        yaml.safe_dump(cfg, f, sort_keys=False)
                    bc = json.loads(json.dumps(bc0))
                    fl = bc["inlet"]["floats"]
                    for k in [k for k in fl if k.startswith("X")]:
                        del fl[k]
                    for s, x in enumerate(Xbc):
                        fl[f"X{s}"] = float(x)
                    with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
                        yaml.safe_dump(bc, f, sort_keys=False)
                    e = dict(os.environ)
                    e["FORGE_PROFILE"] = "1"
                    e["FORGE_CUDA_BLOCKSIZE"] = "256"
                    if var == "double":
                        e["FORGE_TRANSPORT_TABLE"] = "0"
                    r = subprocess.run([a.forge], cwd=d, env=e, capture_output=True, text=True, timeout=3600)
                    open(os.path.join(d, "forge.log"), "w").write(r.stdout + r.stderr)
                    if r.returncode != 0 or "Profiled steps" not in r.stdout:
                        print(f"[FAIL] {name} {var} rep {rep} n {nst}: rc {r.returncode}\n{(r.stdout + r.stderr)[-2000:]}")
                        return 1
                    if "nan" in r.stdout.lower().split("=== runtime profile summary ===")[0][-5000:]:
                        print(f"[WARN] {name} {var} rep {rep} n {nst}: 'nan' in log tail")
                    pr = parse_profile(r.stdout)
                    if pr["steps"] != nst or "Non-finite value detected" in r.stdout:
                        print(f"[FAIL] {name} {var} rep {rep}: profiled {pr['steps']} of {nst} steps (NaN halt?)")
                        return 1
                    results[(name, var, rep, nst)] = pr
                p1, p0 = results[(name, var, rep, a.warm + a.steps)], results[(name, var, rep, a.warm)]
                gas = (p1["gas"] - p0["gas"]) / a.steps
                step = (p1["step"] - p0["step"]) / a.steps
                results[(name, var, rep)] = (gas, step)
                print(f"  {name:4s} {var:6s} rep {rep}: gas_properties {gas:.4f} ms/step, step {step:.4f} ms", flush=True)
        # 判定
        if "legacy" in a.variants and "table" in a.variants:
            ratios = []
            for rep in range(a.reps):
                g_old, s_old = results[(name, "legacy", rep)]
                g_new, _ = results[(name, "table", rep)]
                ratios.append((g_new - g_old) / s_old)
            rd = ""
            if "double" in a.variants:
                rds = [(results[(name, "double", rep)][0] - results[(name, "legacy", rep)][0]) / results[(name, "legacy", rep)][1]
                       for rep in range(a.reps)]
                rd = f"; stage-2 double path (record) max {max(rds):+.4f}"
            gate = name in ("n5", "n12")
            ok = max(ratios) <= 0.10
            tag = ("[PASS] " if ok else "[FAIL] ") if gate else "[INFO] "
            print(f"{tag}{name}: (new gas - old gas)/old step per rep {['%+.4f' % x for x in ratios]}, max {max(ratios):+.4f} "
                  f"({'<= 0.10' if gate else 'record only'}){rd}", flush=True)
    if not a.keep:
        shutil.rmtree(root, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
