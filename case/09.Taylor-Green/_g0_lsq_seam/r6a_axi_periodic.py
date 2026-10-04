#!/usr/bin/env python3
"""#6a (plan boundary-node-periodic-gradient-fix §5.1 #6a、codex result M1): 軸対称 × 並進周期で化学種・受動種の勾配経路が旧新で不変。

旧 (`1266aba1`) では `species_gradient_d` の周期半割面除外は死に条件 (§4.2a) で一度も除外しなかった。
新 (`bd22376d`) は除外を `periodicSeamMergeActive` (軸対称では false) に限ったので、軸対称 × 周期でも除外しない
→ 1 step 後の roY*・roXi (と全 VALUE) は旧新で同じはず。

メッシュ: 軸対称 2D の上半分 x∈[0,0.1]、r∈[0,0.02] (61×41 節点)。x は並進周期 (left ↔ right)、r=0 は軸、r=0.02 は slip。
場: R3 の axi (`r3_prepare.py` bake_axi) と同じ作り方で、x 方向は周期 (cos 2πx/L) にする。2 成分 (MIXDRY/H2O) +
受動トレーサ ξ、`time.deltaT.speciesFaceReconstruction: 1`、乱流なし (SST の F1 初期値 m3 の変化を入れない)。

判定 (§5.1 #6a 合格「旧新の 1 step 後 roY/roXi がビット同一」、揺れの扱いは §6 R3 の規則):
  res_1 の全 VALUE を旧 2・新 2 で比較。旧同士がビット一致する配列は旧新もビット一致を要求。
  旧同士でも一致しない配列は、旧新の最大差 ≤ 旧同士の最大差 × 2 かつ不一致数の桁が同じ。
  VERDICT は roY*・roXi を含む全 VALUE (res_1) が上の規則を満たすこと。

使い方:
  python3 r6a_axi_periodic.py prepare <scratch> --species-db PATH
  FORGE_BIN=... python3 r6a_axi_periodic.py run <scratch> <old_a|old_b|new_a|new_b>
  python3 r6a_axi_periodic.py compare <scratch> [--out R6a_axi_periodic.txt]
  python3 r6a_axi_periodic.py addruns <scratch> old_c old_d ... ; (run) ; python3 r6a_axi_periodic.py noise <scratch> --out ...
    (参考: 反復 run でノイズ床の分布を出す。VERDICT は変えない)
"""
import argparse
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r3_compare import cmp, judge, fmt   # noqa: E402  (§6 R3 と同じ規則)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CONVERTER = os.path.join(REPO, "solver_density_cuda", "build", "convertGmshToForge")
QUALITY = os.path.join(REPO, "solver_density_cuda", "tools", "check_mesh_quality.py")
RUN_CASE = os.path.join(REPO, "solver_density_cuda", "tools", "run_case.sh")
RUNS = ("old_a", "old_b", "new_a", "new_b")
PAIRS = (("old_a", "old_b"), ("new_a", "new_b"), ("old_a", "new_a"))
LX, RMAX = 0.1, 0.02

GEO = """// #6a 用の軸対称 × 並進周期 node メッシュ (上半分 2D、r>=0)。physical: left 1 / right 2 / top 3 / axis 4 / fluid 5
Mesh.MshFileVersion = 4.1;
Point(1) = {0,   0,    0, 1};
Point(2) = {%(LX)g, 0,    0, 1};
Point(3) = {%(LX)g, %(R)g, 0, 1};
Point(4) = {0,   %(R)g, 0, 1};
Line(1) = {1, 2};   // axis
Line(2) = {2, 3};   // right (periodic)
Line(3) = {3, 4};   // top (slip)
Line(4) = {4, 1};   // left (periodic)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 3} = 61;
Transfinite Curve {2, 4} = 41;
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface{1};
Physical Curve("left", 1)  = {4};
Physical Curve("right", 2) = {2};
Physical Curve("top", 3)   = {3};
Physical Curve("axis", 4)  = {1};
Physical Surface("fluid", 5) = {1};
""" % {"LX": LX, "R": RMAX}

BCOND = """left:  {physID: 1, kind: periodic, outputHDFflg: 0, ints: {type: 0, partnerBCID: 2}, floats: {dx: %(LX)r,  dy: 0.0, dz: 0.0}}
right: {physID: 2, kind: periodic, outputHDFflg: 0, ints: {type: 0, partnerBCID: 1}, floats: {dx: -%(LX)r, dy: 0.0, dz: 0.0}}
top:   {physID: 3, kind: slip,     outputHDFflg: 0, ints: , floats: }
axis:  {physID: 4, kind: axis,     outputHDFflg: 0, ints: , floats: }
""" % {"LX": LX}

# R3 の axi 設定 (case/44 run_0162 から凝縮を外したもの、`r3_prepare.py` prep_axi) + 受動トレーサ
CFG = {
    "mesh": {"discretization": "node", "isAxisymmetric": 1, "axisCentroidShift": 1,
             "meshFileName": "axi.h5", "valueFileName": "axi.h5"},
    "gpu": 1, "solver": "SLAU",
    "physProp": {"thermalMethod": 2, "viscMethod": 0, "visc": 0.0, "thermCond": 0.0, "cp": 1220.7, "gamma": 1.31526,
                 "species": ["MIXDRY", "H2O"], "speciesDBFile": "species_db.yaml", "thermoHrefTemp": 298.15,
                 "tracer": "exhaust"},
    "time": {"unsteady": 0, "dualTime": 0, "last": {"nStepOuter": 1},
             "deltaT": {"control": 1, "dt": 1e-8, "cfl": 0.5, "cfl_pseudo": 0.5, "dt_min": 1e-9, "dt_max": 0.001,
                        "blockDPLUR": 1, "lowMachPrecond": 0, "detectNaN": 1, "speciesFaceReconstruction": 1},
             "outStepStart": 0, "outStepInterval": 1, "timeIntegration": 11, "nStepInner": 5},
    "space": {"convMethod": 1, "limiter": 2},
    "turbulence": {"model": "none"},
    "initial": "uniform_p101325_u10",
    "output": {"level": 2},
}


def bake(h5):
    """一様状態の e0 を保って ρ・u・Y・ξ を滑らかに変える。x は周期 (cos/sin 2πx/L)、u_r は軸で 0。"""
    with h5py.File(h5, "r+") as f:
        V = f["VALUE"]
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        n = V["ro"].shape[0]
        if xyz.shape[0] != n:
            raise SystemExit(f"MESH/COORD と VALUE の長さが合わない ({xyz.shape}, {n})")
        th = 2.0 * np.pi * xyz[:, 0].astype(np.float64) / LX
        r = xyz[:, 1].astype(np.float64) / RMAX
        ro0 = float(np.array(V["ro"])[0]); u0 = float(np.array(V["roUx"])[0]) / ro0
        e0 = float(np.array(V["roe"])[0]) / ro0 - 0.5 * u0 ** 2
        ro = ro0 * (1.0 + 0.10 * np.cos(th) * (1.0 - 0.3 * r ** 2))
        ux = 100.0 * (1.0 + 0.3 * np.sin(th)) * (1.0 - 0.4 * r ** 2)
        ur = 8.0 * np.sin(th) * r
        y1 = 0.02 + 0.01 * np.cos(th) * (1.0 - r ** 2)
        xi = 0.5 + 0.3 * np.sin(th + 0.7) * (1.0 - 0.5 * r ** 2)
        roe = ro * (e0 + 0.5 * (ux ** 2 + ur ** 2))
        out = {"ro": ro, "roUx": ro * ux, "roUy": ro * ur, "roUz": 0.0 * ro, "roe": roe,
               "roY0": ro * (1.0 - y1), "roY1": ro * y1, "roXi": ro * xi}
        dt = V["ro"].dtype
        for k, a in out.items():
            if k in V:
                V[k][...] = a.astype(V[k].dtype)
            else:                                # 変換器は一様 IC に roY/roXi を書かない
                V.create_dataset(k, data=a.astype(dt))


def cmd_prepare(scratch, species_db):
    prep = os.path.join(scratch, "r6a_prep")
    os.makedirs(prep, exist_ok=False)
    with open(os.path.join(prep, "solverConfig.yaml"), "w") as fp:
        yaml.safe_dump(CFG, fp, sort_keys=False, default_flow_style=None)
    with open(os.path.join(prep, "bcondConfig.yaml"), "w") as fp:
        fp.write(BCOND)
    with open(os.path.join(prep, "probe.yaml"), "w") as fp:
        fp.write("outStepInterval: 1\noutStepStart: 0\npoints:\nsurfaces:\n")   # forge は probe.yaml が無いと起動しない
    shutil.copyfile(species_db, os.path.join(prep, "species_db.yaml"))
    with open(os.path.join(prep, "axi.geo"), "w") as fp:
        fp.write(GEO)
    subprocess.run(["gmsh", "-2", "axi.geo", "-o", "axi.msh", "-format", "msh4"], cwd=prep, check=True,
                   capture_output=True)
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
    r = subprocess.run([CONVERTER, "axi.msh", "axi.h5"], cwd=prep, capture_output=True, text=True, env=env)
    with open(os.path.join(prep, "convert.log"), "w") as fp:
        fp.write(r.stdout + r.stderr)
    tail = (r.stdout + r.stderr).strip().splitlines()[-1:] or [""]
    # AWS g5 の既知の罠: 変換器は終了時の cudaFree で GPUassert (exit≠0) になるが出力 h5 は完全 (memory aws-p1-instance-state)
    known_exit = "writeInputH5: wrote" in r.stdout and tail[0].startswith("GPUassert") and os.path.exists(os.path.join(prep, "axi.h5"))
    if r.returncode != 0 and not known_exit:
        raise SystemExit(f"convert failed: see {prep}/convert.log")
    q = subprocess.run(["python3", QUALITY, "axi.h5"], cwd=prep, capture_output=True, text=True)
    with open(os.path.join(prep, "quality.txt"), "w") as fp:
        fp.write(q.stdout + q.stderr)
    if "VERDICT: PASS" not in q.stdout:
        raise SystemExit(f"mesh quality not PASS: see {prep}/quality.txt")
    bake(os.path.join(prep, "axi.h5"))
    for v in RUNS:
        rd = os.path.join(scratch, f"r6a_axiper_{v}")
        os.makedirs(rd, exist_ok=False)
        for fn in ("axi.h5", "solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml"):
            shutil.copyfile(os.path.join(prep, fn), os.path.join(rd, fn))
    print("prepared", prep)


def cmd_run(scratch, variant):
    rd = os.path.join(scratch, f"r6a_axiper_{variant}")
    if os.path.exists(os.path.join(rd, "res_1.h5")):
        raise SystemExit(f"既に回した跡がある: {rd}")
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_CUDA_BLOCKSIZE_SMALL="128")
    r = subprocess.run(["bash", RUN_CASE, rd], env=env, capture_output=True, text=True)
    print(rd, "rc", r.returncode)


def read_values(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def cmd_compare(scratch, out):
    lines = []
    P = lines.append
    P("#6a: 軸対称 × 並進周期で 1 step 後の VALUE が旧新で不変 (plan boundary-node-periodic-gradient-fix §5.1 #6a、codex result M1)")
    P(f"run root: {os.path.abspath(scratch)}/r6a_axiper_<old|new>_<a|b>")
    for v in ("old_a", "new_a"):
        p = os.path.join(scratch, f"r6a_axiper_{v}", "RUN_PROVENANCE.txt")
        if os.path.exists(p):
            P(f"{v}: " + " ".join(l.strip() for l in open(p) if l.startswith(("forge_sha256", "git_head"))))
    q = os.path.join(scratch, "r6a_prep", "quality.txt")
    if os.path.exists(q):
        P("mesh quality: " + ([l.strip() for l in open(q) if "VERDICT" in l] or ["?"])[-1])
    log = os.path.join(scratch, "r6a_axiper_new_a", "forge_run.log")
    if os.path.exists(log):
        for l in open(log):
            if any(s in l for s in ("setPeriodicPartner", "buildPeriodicNodeGroups", "periodic seam", "isAxisymmetric", "[passive] initial")):
                P("new_a log: " + l.rstrip())
    P("列: 不一致数 (ビット比較) / 最大絶対差。判定は oo=old_a/old_b と on=old_a/new_a に §6 R3 の規則を当てる。")
    vals = {s: {r: read_values(os.path.join(scratch, f"r6a_axiper_{r}", f"res_{s}.h5")) for r in RUNS} for s in (0, 1)}
    bad = []
    for s in (0, 1):
        for r in RUNS:
            for k, a in vals[s][r].items():
                if a.dtype.kind == "f" and not np.isfinite(a).all():
                    bad.append(f"res_{s} {r} {k}: 非有限 {int((~np.isfinite(a)).sum())}")
    P("NaN/Inf: " + ("なし (全 run・res_0/res_1 の全 VALUE)" if not bad else "; ".join(bad)))
    # 入力 (res_0 の保存量) が 4 run で同一であることの確認
    ok_all = not bad
    for s in (1, 0):
        names = sorted(vals[s]["old_a"])
        P("")
        P(f"-- res_{s}.h5 全 VALUE ({len(names)} 配列)" + ("  [VERDICT 対象]" if s == 1 else "  [参考: step 前]"))
        P(f"{'array':<14} {'old_a/old_b':>18} {'new_a/new_b':>18} {'old_a/new_a':>18}  判定")
        for k in names:
            if any(k not in vals[s][r] for r in RUNS):
                P(f"{k:<14} 一部の run に無い (" + ",".join(r for r in RUNS if k not in vals[s][r]) + ")")
                if s == 1:
                    ok_all = False
                continue
            res = {p: cmp(vals[s][p[0]][k], vals[s][p[1]][k]) for p in PAIRS}
            j = judge(res[PAIRS[0]], res[PAIRS[2]])
            if s == 1 and not j.startswith("PASS"):
                ok_all = False
            mark = "  <- roY/roXi" if (k.startswith("roY") or k == "roXi") else ""
            P(f"{k:<14} {fmt(res[PAIRS[0]])} {fmt(res[PAIRS[1]])} {fmt(res[PAIRS[2]])}  {j}{mark}")
    need = [k for k in ("roY0", "roY1", "roXi") if k not in vals[1]["old_a"]]
    if need:
        P(f"注: res_1 に {need} が無い → 判定不能")
        ok_all = False
    P("")
    P(f"VERDICT #6a (軸対称 × 並進周期、res_1 全 VALUE の旧新): {'PASS' if ok_all else 'FAIL'}")
    txt = "\n".join(lines) + "\n"
    with open(out, "w") as fp:
        fp.write(txt)
    sys.stdout.write(txt)


def cmd_addruns(scratch, names):
    """ノイズ床の分布を測る追加の反復 run (参考。VERDICT は old_a/old_b/new_a/new_b の §6 R3 規則のまま)。"""
    prep = os.path.join(scratch, "r6a_prep")
    for v in names:
        rd = os.path.join(scratch, f"r6a_axiper_{v}")
        os.makedirs(rd, exist_ok=False)
        for fn in ("axi.h5", "solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml"):
            shutil.copyfile(os.path.join(prep, fn), os.path.join(rd, fn))


def cmd_noise(scratch, out):
    """参考: 全反復 run の res_1 を総当たりで比べ、旧同士・新同士・旧新の不一致数と最大差の分布を出す (判定なし)。"""
    import glob
    import itertools
    runs = sorted(os.path.basename(p)[len("r6a_axiper_"):] for p in glob.glob(os.path.join(scratch, "r6a_axiper_*"))
                  if os.path.exists(os.path.join(p, "res_1.h5")))
    olds = [r for r in runs if r.startswith("old_")]; news = [r for r in runs if r.startswith("new_")]
    vals = {r: read_values(os.path.join(scratch, f"r6a_axiper_{r}", "res_1.h5")) for r in runs}
    groups = {"旧同士": list(itertools.combinations(olds, 2)), "新同士": list(itertools.combinations(news, 2)),
              "旧新": [(a, b) for a in olds for b in news]}
    lines = [f"\n## 参考: 反復 run 総当たりの res_1 不一致 (判定なし)。旧 {len(olds)} 本 ({','.join(olds)})・新 {len(news)} 本 ({','.join(news)})",
             "列: 不一致数 min/中央値/max、最大絶対差の max。旧同士でも一致しない配列だけ出す (他は全組ビット一致)"]
    lines.append(f"{'array':<14} " + " ".join(f"{g + ' 不一致数 (min/med/max)':>30} {'最大差':>10}" for g in groups))
    for k in sorted(vals[runs[0]]):
        st = {}
        for g, prs in groups.items():
            ns, ms = [], []
            for a, b in prs:
                n, m = cmp(vals[a][k], vals[b][k]); ns.append(n); ms.append(m)
            st[g] = (min(ns), int(np.median(ns)), max(ns), max(ms))
        if all(v[2] == 0 for v in st.values()):
            continue
        lines.append(f"{k:<14} " + " ".join(f"{f'{v[0]}/{v[1]}/{v[2]}':>30} {v[3]:10.3e}" for v in st.values()))
    txt = "\n".join(lines) + "\n"
    with open(out, "a") as fp:
        fp.write(txt)
    sys.stdout.write(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prepare", "run", "compare", "addruns", "noise"))
    ap.add_argument("scratch")
    ap.add_argument("variant", nargs="*")
    ap.add_argument("--species-db")
    ap.add_argument("--out", default=os.path.join(HERE, "R6a_axi_periodic.txt"))
    a = ap.parse_args()
    if a.cmd == "prepare":
        cmd_prepare(os.path.abspath(a.scratch), a.species_db)
    elif a.cmd == "run":
        cmd_run(os.path.abspath(a.scratch), a.variant[0])
    elif a.cmd == "addruns":
        cmd_addruns(os.path.abspath(a.scratch), a.variant)
    elif a.cmd == "noise":
        cmd_noise(os.path.abspath(a.scratch), a.out)
    else:
        cmd_compare(os.path.abspath(a.scratch), a.out)


if __name__ == "__main__":
    main()
