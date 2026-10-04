#!/usr/bin/env python3
"""G0/G1/G2' 用の小さな構造格子 (gmsh → 節点ジッタ → node 変換 → 品質チェック)。

- `box`: TGV と同じ三重周期立方体 [0,2π]^3、N^3 セル (`mesh/Taylor-Green.geo` と同じ面番号)。
- `channel`: x・z 周期、y 両壁の箱 [0,1]x[0,1]x[0,0.5]。x は等比 (継ぎ目の両側で間隔が違う = 非対称 stencil)、
  y は両壁へ Bump。壁∩継ぎ目の節点を作る。
- ジッタ (G2'): 節点を ±frac·h の決定論的擬似乱数で動かす。擬似乱数は格子 index を周期で折り返した値
  (i mod N, j mod N, k mod N) のハッシュなので、**周期像は同じ量だけ動く**。境界面上の節点は面内方向だけ動かす
  (箱の面を平面のまま保つ。角の節点は動かない)。

使い方 (関数として g_suite.py / g2p_jitter.py から呼ぶ):
    make_box_h5(workdir, N, jitter=0.0)      -> (h5, bcond)
    make_channel_h5(workdir)                 -> (h5, bcond)
"""
import os
import shutil
import subprocess

import numpy as np
import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CONVERTER = os.path.join(REPO, "solver_density_cuda", "build", "convertGmshToForge")
QUALITY = os.path.join(REPO, "solver_density_cuda", "tools", "check_mesh_quality.py")
HERE = os.path.dirname(os.path.abspath(__file__))
TGV_BCOND = os.path.join(HERE, "bcondConfig.yaml")

GEO = """Mesh.MshFileVersion = 4.1;
nx={nx}; ny={ny};
Point(1) = {{0   ,0   ,0,1}};
Point(2) = {{{Lx},0   ,0,1}};
Point(3) = {{{Lx},{Ly},0,1}};
Point(4) = {{0   ,{Ly},0,1}};
Line(1) = {{1,2}}; Transfinite Line {{1}} = nx {px1};
Line(2) = {{2,3}}; Transfinite Line {{2}} = ny {py};
Line(3) = {{3,4}}; Transfinite Line {{3}} = nx {px3};
Line(4) = {{4,1}}; Transfinite Line {{4}} = ny {py};
Line Loop(1) = {{1,2,3,4}};
Plane Surface(1) = {{1}};
Transfinite Surface {{1}};
Recombine Surface(1);
Extrude {{0, 0, {Lz}}} {{ Surface{{1}}; Layers{{{nz}}}; Recombine; }}
Physical Surface("bottom", 1) = {{13}};
Physical Surface("upper", 2) = {{21}};
Physical Surface("left", 3) = {{26}};
Physical Surface("right", 4) = {{1}};
Physical Surface("inlet", 5) = {{25}};
Physical Surface("outlet", 6) = {{17}};
Physical Volume("fluid", 7) = {{1}};
"""

CONV_CFG = {"mesh": {"discretization": "node", "meshFileName": "mesh.h5", "valueFileName": "mesh.h5"}, "gpu": 1, "solver": "SLAU",
            "physProp": {"thermalMethod": 0, "viscMethod": 0, "visc": 0.0, "thermCond": 0.0, "cp": 1038.8, "gamma": 1.4},
            "time": {"unsteady": 0, "dualTime": 0, "last": {"nStepOuter": 1},
                     "deltaT": {"control": 1, "dt": 1e-9, "cfl": 0.01, "cfl_pseudo": 0.01, "dt_min": 1e-12, "dt_max": 1e-7},
                     "outStepStart": 0, "outStepInterval": 1, "timeIntegration": 11, "nStepInner": 5},
            "space": {"convMethod": 0, "limiter": 0}, "turbulence": {"model": "sst"},
            "initial": "uniform_p101325_u10"}


def _hash01(i, j, k, axis, seed=12345):
    """整数格子 index の決定論的ハッシュ → [-1, 1)。"""
    h = (np.asarray(i, np.uint64) * np.uint64(73856093)) ^ (np.asarray(j, np.uint64) * np.uint64(19349663)) \
        ^ (np.asarray(k, np.uint64) * np.uint64(83492791)) ^ np.uint64(seed + 1000003 * (axis + 1))
    h = (h ^ (h >> np.uint64(13))) * np.uint64(0x5bd1e995) & np.uint64(0xFFFFFFFF)
    h = h ^ (h >> np.uint64(15))
    return (h.astype(np.float64) / 2.0 ** 32) * 2.0 - 1.0


def _jitter_msh(path, N, L, frac):
    """msh 4.1 ASCII の $Nodes を書き換える (等間隔 box 専用: index = round(x/h))。"""
    lines = open(path).read().split("\n")
    s = lines.index("$Nodes")
    nblk = int(lines[s + 1].split()[0])
    p = s + 2
    h = L / N
    for _ in range(nblk):
        _, _, param, nn = map(int, lines[p].split()); p += 1
        assert param == 0
        p += nn                                          # tags
        for q in range(p, p + nn):
            x = np.array(list(map(float, lines[q].split())))
            idx = np.rint(x / h).astype(np.int64)
            onb = (idx == 0) | (idx == N)
            d = np.array([_hash01(idx[0] % N, idx[1] % N, idx[2] % N, a) for a in range(3)]) * frac * h
            d[onb] = 0.0
            x = x + d
            lines[q] = "%.17g %.17g %.17g" % tuple(x)
        p += nn
    open(path, "w").write("\n".join(lines))


def _convert(workdir, geo_text, bcond, jitter_fn=None):
    os.makedirs(workdir, exist_ok=True)
    open(os.path.join(workdir, "mesh.geo"), "w").write(geo_text)
    if shutil.which("gmsh") is None and os.path.exists(os.path.join(workdir, "mesh.msh.jittered")):
        # gmsh の無い機械 (AWS g5) 用: 別の機械で gmsh + ジッタまで済ませた mesh.msh.jittered を置いておけば変換だけ行う
        shutil.copy(os.path.join(workdir, "mesh.msh.jittered"), os.path.join(workdir, "mesh.msh"))
    else:
        subprocess.run(["gmsh", "-3", "mesh.geo", "-o", "mesh.msh", "-format", "msh4"], cwd=workdir, check=True,
                       capture_output=True)
        if jitter_fn:
            jitter_fn(os.path.join(workdir, "mesh.msh"))
        shutil.copy(os.path.join(workdir, "mesh.msh"), os.path.join(workdir, "mesh.msh.jittered"))
    with open(os.path.join(workdir, "solverConfig.yaml"), "w") as fp:
        yaml.safe_dump(CONV_CFG, fp, sort_keys=False)
    if isinstance(bcond, str):
        shutil.copy(bcond, os.path.join(workdir, "bcondConfig.yaml"))
    else:
        with open(os.path.join(workdir, "bcondConfig.yaml"), "w") as fp:
            yaml.safe_dump(bcond, fp, sort_keys=False)
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
    r = subprocess.run([CONVERTER, "mesh.msh", "mesh.h5"], cwd=workdir, capture_output=True, text=True, env=env)
    open(os.path.join(workdir, "convert.log"), "w").write(r.stdout + r.stderr)
    # AWS g5 の既知の罠: 変換器は終了時の cudaFree で GPUassert (exit≠0) になるが出力 h5 は完全 (memory aws-p1-instance-state)
    tail = (r.stdout + r.stderr).strip().splitlines()[-1:] or [""]
    known_exit = "writeInputH5: wrote" in r.stdout and tail[0].startswith("GPUassert")
    if (r.returncode != 0 and not known_exit) or not os.path.exists(os.path.join(workdir, "mesh.h5")):
        raise SystemExit("convert failed: " + workdir)
    q = subprocess.run(["python3", QUALITY, "mesh.h5"], cwd=workdir, capture_output=True, text=True)
    open(os.path.join(workdir, "quality.txt"), "w").write(q.stdout + q.stderr)
    verdict = [l for l in (q.stdout + q.stderr).splitlines() if "VERDICT" in l]
    return os.path.join(workdir, "mesh.h5"), os.path.join(workdir, "bcondConfig.yaml"), (verdict[-1] if verdict else "?")


def make_box_h5(workdir, N, jitter=0.0, bcond=TGV_BCOND):
    L = 2.0 * np.pi
    geo = GEO.format(nx=N + 1, ny=N + 1, nz=N, Lx="2*Pi", Ly="2*Pi", Lz="2*Pi", px1="", px3="", py="")
    fn = (lambda p: _jitter_msh(p, N, L, jitter)) if jitter > 0 else None
    return _convert(workdir, geo, bcond, fn)


def slip_box_bcond():
    """周期の無い箱 (全 6 面 slip)。plan gradient-scalar-lsq-unification S0-c: 周期 gather が no-op なので
    res_1 の最終配列同士を直接比べられる。"""
    names = ("bottom", "upper", "left", "right", "inlet", "outlet")
    return {nm: {"physID": i + 1, "kind": "slip", "outputHDFflg": 0, "ints": None, "floats": None} for i, nm in enumerate(names)}


def channel_bcond(Lx=1.0, Lz=0.5):
    wall = {"kind": "wall", "outputHDFflg": 0, "ints": None, "floats": {"Ux": 0.0, "Uy": 0.0, "Uz": 0.0}}
    return {
        "bottom": dict(wall, physID=1), "upper": dict(wall, physID=2),
        "left": {"physID": 3, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 4},
                 "floats": {"dx": 0.0, "dy": 0.0, "dz": -Lz}},
        "right": {"physID": 4, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 3},
                  "floats": {"dx": 0.0, "dy": 0.0, "dz": Lz}},
        "inlet": {"physID": 5, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 6},
                  "floats": {"dx": Lx, "dy": 0.0, "dz": 0.0}},
        "outlet": {"physID": 6, "kind": "periodic", "outputHDFflg": 0, "ints": {"type": 0, "partnerBCID": 5},
                   "floats": {"dx": -Lx, "dy": 0.0, "dz": 0.0}},
    }


def make_channel_h5(workdir, nx=17, ny=13, nz=8, r=1.1):
    geo = GEO.format(nx=nx, ny=ny, nz=nz, Lx=1.0, Ly=1.0, Lz=0.5,
                     px1="Using Progression %.6f" % r, px3="Using Progression %.6f" % (1.0 / r), py="Using Bump 0.4")
    return _convert(workdir, geo, channel_bcond())
