#!/usr/bin/env python3
"""R3 (plan boundary-node-periodic-gradient-fix §6) の入力を作る: 3 ケース × (旧 2 回・新 2 回) の run ディレクトリ。

各ケースで入力 h5 と config を 1 組だけ作り、4 run へ cp する (ビット同一)。run は FORGE_BIN だけを変えて回す。
1 step (nStepOuter 1)、outStepStart 0 / outStepInterval 1 (res_0・res_1)、output.level 2。

- case48: `case/48.flat_plate_cooled_m4/mesh/fp_y1_12um.h5` (node・非周期) + `run_0030_bitrep_a` の設定を SST に。
  メッシュの VALUE は roK=roOmega=0 の一様場なので roK/roOmega を落とし、kInit/omegaInit = 入口 bcond の k/ω で入れる。
- case16: `case/16.nozzle_wys/run_0507_node_lam_flag1_br` の h5 (発達場) と設定 + `time.deltaT.speciesFaceReconstruction: 1`
  (化学種 GG `species_gradient_d` を通す。キーの位置は `procedures/solver-settings.md`「speciesFaceReconstruction」)。
- axi: case/44 の入力メッシュは削除済みのため、同じ軸対称 node 経路 (`mesh.isAxisymmetric: 1`、`axisCentroidShift: 1`、
  bcond `axis`) を通す小さな平面 2D 構造格子 (x∈[0,0.1]、r∈[0,0.02]、101×41 節点) を作る。設定は case/44 run_0162 から
  凝縮を外したもの + speciesFaceReconstruction 1。一様場では 1 step の勾配が境界近傍にしか立たないので、変換後の一様状態から
  比内部エネルギー e0 を保ったまま ρ・u・Y を滑らかに変えた場を焼き込む (T は組成の差だけ動く)。
  e0 は変換器の一様 IC (CPG 式の roe) をそのまま TP で読むので T は ~700 K 級になる (1 step の勾配比較には無関係)。

使い方: python3 r3_prepare.py <scratch_dir>    (出力: <scratch_dir>/r3_prep/<case>/、<scratch_dir>/r3_<case>_<old|new>_<a|b>/)
"""
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CONVERTER = os.path.join(REPO, "solver_density_cuda", "build", "convertGmshToForge")
QUALITY = os.path.join(REPO, "solver_density_cuda", "tools", "check_mesh_quality.py")
CASES = ("case48", "case16", "axi")
VARIANTS = ("old_a", "old_b", "new_a", "new_b")


def _load(path):
    with open(path) as fp:
        return yaml.safe_load(fp)


def _dump(obj, path):
    with open(path, "w") as fp:
        yaml.safe_dump(obj, fp, sort_keys=False, default_flow_style=None)


def _one_step(cfg):
    cfg["time"]["last"] = {"nStepOuter": 1}
    cfg["time"]["outStepStart"] = 0
    cfg["time"]["outStepInterval"] = 1
    cfg["output"] = {"level": 2}
    return cfg


def prep_case48(d):
    src = os.path.join(REPO, "case", "48.flat_plate_cooled_m4")
    base = os.path.join(src, "run_0030_bitrep_a")
    shutil.copyfile(os.path.join(src, "mesh", "fp_y1_12um.h5"), os.path.join(d, "mesh.h5"))
    with h5py.File(os.path.join(d, "mesh.h5"), "r+") as f:
        for k in ("roK", "roOmega"):          # 一様 0 → kInit/omegaInit で入れる (ω=0 は SST の初手発散)
            if k in f["VALUE"]:
                del f["VALUE"][k]
    for fn in ("bcondConfig.yaml", "probe.yaml"):
        shutil.copyfile(os.path.join(base, fn), os.path.join(d, fn))
    inlet = _load(os.path.join(base, "bcondConfig.yaml"))["inlet"]["floats"]
    cfg = _one_step(_load(os.path.join(base, "solverConfig.yaml")))
    cfg["turbulence"]["model"] = "sst"
    cfg["turbulence"]["wallTreatmentSST"] = 0
    cfg["turbulence"]["kInit"] = float(inlet["k"])
    cfg["turbulence"]["omegaInit"] = float(inlet["omega"])
    _dump(cfg, os.path.join(d, "solverConfig.yaml"))
    return "mesh.h5", ["probe.yaml"]


def prep_case16(d):
    base = os.path.join(REPO, "case", "16.nozzle_wys", "run_0507_node_lam_flag1_br")
    shutil.copyfile(os.path.join(base, "nozzle_user_2d.h5"), os.path.join(d, "nozzle_user_2d.h5"))
    for fn in ("bcondConfig.yaml", "species_db.yaml", "probe.yaml"):
        shutil.copyfile(os.path.join(base, fn), os.path.join(d, fn))
    cfg = _one_step(_load(os.path.join(base, "solverConfig.yaml")))
    cfg["time"]["deltaT"]["speciesFaceReconstruction"] = 1
    _dump(cfg, os.path.join(d, "solverConfig.yaml"))
    return "nozzle_user_2d.h5", ["species_db.yaml", "probe.yaml"]


AXI_GEO = """// R3 用の小さな軸対称 node メッシュ (上半分 2D、r>=0)。physical: inlet 1 / outlet 2 / wall 3 / axis 4 / fluid 5
Mesh.MshFileVersion = 4.1;
Point(1) = {0,   0,    0, 1};
Point(2) = {0.1, 0,    0, 1};
Point(3) = {0.1, 0.02, 0, 1};
Point(4) = {0,   0.02, 0, 1};
Line(1) = {1, 2};   // axis
Line(2) = {2, 3};   // outlet
Line(3) = {3, 4};   // wall
Line(4) = {4, 1};   // inlet
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 3} = 101;
Transfinite Curve {2, 4} = 41;
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface{1};
Physical Curve("inlet", 1)  = {4};
Physical Curve("outlet", 2) = {2};
Physical Curve("wall", 3)   = {3};
Physical Curve("axis", 4)   = {1};
Physical Surface("fluid", 5) = {1};
"""

AXI_BCOND = """inlet:  {physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {Y0: 0.98, Y1: 0.02, Pt: 106000.0, Tt: 305.0, k: 1.0, omega: 1000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 95000.0, Pt: 95000.0, Tt: 300.0}}
wall:   {physID: 3, kind: slip,             outputHDFflg: 0, ints: , floats: }
axis:   {physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }
"""


def prep_axi(d):
    base44 = os.path.join(REPO, "case", "44.vitiated_air_wt", "run_0162_va3_M4.19_Lc8_noneq_inletTt_lim1d_cfl2")
    cfg = _one_step(_load(os.path.join(base44, "solverConfig.yaml")))
    cfg.pop("condensation", None)
    cfg["mesh"]["meshFileName"] = "axi.h5"
    cfg["mesh"]["valueFileName"] = "axi.h5"
    cfg["physProp"]["speciesDBFile"] = "species_db.yaml"
    cfg["time"]["deltaT"]["cfl"] = 0.5
    cfg["time"]["deltaT"]["cfl_pseudo"] = 0.5
    cfg["time"]["deltaT"]["speciesFaceReconstruction"] = 1
    _dump(cfg, os.path.join(d, "solverConfig.yaml"))
    with open(os.path.join(d, "bcondConfig.yaml"), "w") as fp:
        fp.write(AXI_BCOND)
    shutil.copyfile(os.path.join(REPO, "case", "16.nozzle_wys", "run_0507_node_lam_flag1_br", "species_db.yaml"),
                    os.path.join(d, "species_db.yaml"))
    shutil.copyfile(os.path.join(REPO, "case", "48.flat_plate_cooled_m4", "run_0030_bitrep_a", "probe.yaml"),
                    os.path.join(d, "probe.yaml"))     # forge は probe.yaml が無いと起動しない (点・面は空)
    with open(os.path.join(d, "axi.geo"), "w") as fp:
        fp.write(AXI_GEO)
    subprocess.run(["gmsh", "-2", "axi.geo", "-o", "axi.msh", "-format", "msh4"], cwd=d, check=True,
                   capture_output=True)
    r = subprocess.run([CONVERTER, "axi.msh", "axi.h5"], cwd=d, capture_output=True, text=True)
    with open(os.path.join(d, "convert.log"), "w") as fp:
        fp.write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise SystemExit(f"convert failed: see {d}/convert.log")
    q = subprocess.run(["python3", QUALITY, "axi.h5"], cwd=d, capture_output=True, text=True)
    with open(os.path.join(d, "quality.txt"), "w") as fp:
        fp.write(q.stdout + q.stderr)
    if "VERDICT: PASS" not in q.stdout:
        raise SystemExit(f"mesh quality not PASS: see {d}/quality.txt")
    bake_axi(os.path.join(d, "axi.h5"))
    return "axi.h5", ["species_db.yaml", "probe.yaml"]


def bake_axi(h5):
    """一様状態の e0 (比内部エネルギー) を保って ρ・u・Y を滑らかに変える。u_r は軸で 0。"""
    with h5py.File(h5, "r+") as f:
        V = f["VALUE"]
        xyz = np.array(f["MESH"]["COORD"]).reshape(-1, 3) if "COORD" in f["MESH"] else None
        n = V["ro"].shape[0]
        if xyz is None or xyz.shape[0] != n:
            raise SystemExit(f"MESH/COORD と VALUE の長さが合わない ({None if xyz is None else xyz.shape}, {n})")
        x, r = xyz[:, 0].astype(np.float64) / 0.1, xyz[:, 1].astype(np.float64) / 0.02
        ro0 = float(np.array(V["ro"])[0]); u0 = np.array(V["roUx"])[0] / ro0
        e0 = float(np.array(V["roe"])[0]) / ro0 - 0.5 * u0 ** 2
        ro = ro0 * (1.0 + 0.10 * np.cos(np.pi * x) * (1.0 - 0.3 * r ** 2))
        ux = 100.0 * (1.0 + 0.5 * x) * (1.0 - 0.4 * r ** 2)
        ur = 8.0 * np.sin(np.pi * x) * r
        y1 = 0.02 + 0.01 * np.cos(np.pi * x) * (1.0 - r ** 2)
        roe = ro * (e0 + 0.5 * (ux ** 2 + ur ** 2))
        out = {"ro": ro, "roUx": ro * ux, "roUy": ro * ur, "roUz": 0.0 * ro, "roe": roe,
               "roY0": ro * (1.0 - y1), "roY1": ro * y1}
        dt = V["ro"].dtype
        for k, a in out.items():
            if k in V:
                V[k][...] = a.astype(V[k].dtype)
            else:                                # 変換器は一様 IC に roY を書かない (無いと Y0=1 で起動する)
                V.create_dataset(k, data=a.astype(dt))


def main():
    scratch = os.path.abspath(sys.argv[1])
    prep = os.path.join(scratch, "r3_prep")
    for case, fn in (("case48", prep_case48), ("case16", prep_case16), ("axi", prep_axi)):
        d = os.path.join(prep, case)
        os.makedirs(d, exist_ok=False)
        h5name, extra = fn(d)
        for v in VARIANTS:
            rd = os.path.join(scratch, f"r3_{case}_{v}")
            os.makedirs(rd, exist_ok=False)
            for f in [h5name, "solverConfig.yaml", "bcondConfig.yaml"] + extra:
                shutil.copyfile(os.path.join(d, f), os.path.join(rd, f))
        print(case, "->", d)


if __name__ == "__main__":
    main()
