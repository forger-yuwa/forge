#!/usr/bin/env python3
"""V2c 停滞の次数 A/B (plan boundary-node-farfield-characteristic §5.1 #3c、codex diagnose 2026-10-03 ②)。
旧形状 B (`run_0064_v2c_B`) の最終場を同一格子の新規 2 run へ保存量のままコピーし、`space.convMethod` だけ変える
(A2 = 1 [2 次]、B1 = 0 [1 次])。形状・BC・CFL・リミッタ基準値・バイナリは固定。各 6000 step、全場を 500 step ごと。
判定量 (事前固定): W = max_x(max_t p − min_t p)/Δp (末尾 2000 step、下面の評価線) と、初期・最終の局所残差 Rmax (`v2c_residual_ab.py`、射影なし)。
  python3 v2c_order_ab.py setup BIN | wall
"""
import glob, os, re, shutil, subprocess, sys
import h5py, numpy as np

SRC = "run_0064_v2c_B"; SNAP = f"{SRC}/res_6000.h5"
RUNS = {"run_0150_v2c_ord_A2": 1, "run_0151_v2c_ord_B1": 0}


def setup(binp):
    for dst, cm in RUNS.items():
        os.makedirs(dst)
        for f in ("bcondConfig.yaml", "ramp.h5", "probe.yaml", "GEOM.txt"):
            shutil.copy(os.path.join(SRC, f), dst)
        s = open(os.path.join(SRC, "solverConfig.yaml")).read()
        s = re.sub(r"convMethod: \d+", f"convMethod: {cm}", s)
        s = re.sub(r"outStepInterval: \d+", "outStepInterval: 500", s)
        s = re.sub(r"nStepOuter: \d+", "nStepOuter: 6000", s)
        open(os.path.join(dst, "solverConfig.yaml"), "w").write(s)
        with h5py.File(SNAP) as s5, h5py.File(os.path.join(dst, "ramp.h5"), "r+") as f:
            assert s5["VALUE/ro"].shape[0] == f["VALUE/ro"].shape[0]
            for q in ("ro", "roUx", "roUy", "roUz", "roe"):
                f["VALUE"][q][...] = s5["VALUE"][q][:]
        open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"{SNAP} の保存量を index コピー。convMethod {cm} 以外は {SRC} と同一 (#3c 次数 A/B)\n")
        env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128")
        r = subprocess.run([binp], cwd=dst, env=env, stdin=subprocess.DEVNULL, stdout=open(os.path.join(dst, "forge_run.log"), "w"), stderr=subprocess.STDOUT)
        print(dst, "exit", r.returncode)


def wall():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import eval_v2c as E
    b, p2 = E.oblique(2.5, E.TH); dp = (p2 - 1.0) * 101325.0
    for dst in RUNS:
        fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(dst + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
        last = fs[-1][0]; sel = [p for n, p in fs if n >= last - 2000]
        with h5py.File(os.path.join(dst, "ramp.h5")) as m:
            xyz = np.array(m["MESH/COORD"]).reshape(-1, 3)
        g = open(os.path.join(dst, "GEOM.txt")).read().strip()
        ybx = (np.clip(xyz[:, 0], E.XC, E.XE_OF[g]) - E.XC) * np.tan(E.TH)
        sel_ = (np.abs(xyz[:, 1] - ybx) < 1e-6) & (np.abs(xyz[:, 2] - E.DZ) < 1e-6) & (xyz[:, 0] >= E.XC - 1e-9)
        idx = np.where(sel_)[0][np.argsort(xyz[sel_, 0])]; xs = xyz[idx, 0]   # eval_v2c と同じ評価線 (下面、z = 5 mm)
        P = np.stack([h5py.File(p)["VALUE/P"][:][idx] for p in sel])
        W = (P.max(0) - P.min(0)) / dp
        fin = all(np.all(np.isfinite(h5py.File(p)["VALUE/P"][:])) for p in sel)
        print(f"{dst}: 末尾 {sel[0].split('res_')[-1]}–{last} の {len(sel)} 枚、W = {W.max():.4g} (x {xs[int(W.argmax())]:.3f})、有限 {fin}")


if __name__ == "__main__":
    {"setup": lambda: setup(sys.argv[2]), "wall": wall}[sys.argv[1]]()
