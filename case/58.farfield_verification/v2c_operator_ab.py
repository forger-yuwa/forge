#!/usr/bin/env python3
"""V2c 角の作用素 A/B (plan boundary-node-farfield-characteristic §5.1 #3a、codex diagnose 2026-09-29)。
変えるのは領域の高さだけ: 低領域 = run_0063_v2c_A の格子 (H 0.3、上面 slip)、高領域 = run_0065_v2c_C の格子 (H 0.8、上面 slip)。
両方に C の最終場 (res_6000) を座標対応で与え、1 step の初回評価の帳簿 (FORGE_DUMP_LEDGER) を比べる。
対象 = 下面の角 x 0.195–0.230 (判定) と、再構成に要る近傍 x 0.17–0.26・y ≤ 0.05 (帳簿の印)。
判定 (事前固定): 各保存量の残差差 / その節点に接する内部面流束の絶対和 ≤ 1e-5 (規模 0 は自由流基準 ρ∞U∞ 等 × 面積)。
  python3 v2c_operator_ab.py setup      低・高の run を作って回す
  python3 v2c_operator_ab.py compare    幾何・状態・面・残差を照合して VERDICT を出す
"""
import csv, os, shutil, subprocess, sys
from collections import defaultdict
import h5py, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_LOW, SRC_HIGH = "run_0063_v2c_A", "run_0065_v2c_C"
LOW, HIGH = "run_0092_v2c_opab_low", "run_0093_v2c_opab_high"
SNAP = os.path.join(SRC_HIGH, "res_6000.h5")
BIN = os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff/forge")
CONS = ("ro", "roUx", "roUy", "roUz", "roe")


def coords(run):
    with h5py.File(os.path.join(run, "ramp.h5")) as f:
        return np.array(f["MESH/COORD"]).reshape(-1, 3)


def key(xyz):
    return [tuple(v) for v in np.round(xyz.astype(np.float64), 6)]


def marked(xyz):
    return np.where((xyz[:, 0] >= 0.17) & (xyz[:, 0] <= 0.26) & (xyz[:, 1] <= 0.05))[0]


def setup():
    for dst, src in ((LOW, SRC_LOW), (HIGH, SRC_HIGH)):
        os.makedirs(dst)
        for f in ("solverConfig.yaml", "bcondConfig.yaml", "ramp.h5", "probe.yaml"):
            shutil.copy(os.path.join(src, f), dst)
        s = open(os.path.join(dst, "solverConfig.yaml")).read()
        s = s.replace("nStepOuter: 6000", "nStepOuter: 1").replace("outStepInterval: 1000", "outStepInterval: 1")
        open(os.path.join(dst, "solverConfig.yaml"), "w").write(s)
    # 高領域: 同一格子なので保存量をそのまま、低領域: 座標対応で写す (共通領域の節点は座標が一致する)
    xh, xl = coords(HIGH), coords(LOW)
    idx = {k: i for i, k in enumerate(key(xh))}
    m = np.array([idx[k] for k in key(xl)])
    with h5py.File(SNAP) as s:
        V = {q: s["VALUE"][q][:] for q in CONS}
    for run, sel in ((HIGH, np.arange(len(xh))), (LOW, m)):
        with h5py.File(os.path.join(run, "ramp.h5"), "r+") as f:
            for q in CONS:
                f["VALUE"][q][...] = V[q][sel].astype(f["VALUE"][q].dtype)
        open(os.path.join(run, "IC_FROM.txt"), "w").write(f"V2c 作用素 A/B: {SNAP} の保存量を座標対応で写した ({len(sel)} 節点)\n")
    for run in (LOW, HIGH):
        ids = marked(coords(run))
        env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_DUMP_LEDGER="ledger.csv", FORGE_DUMP_LEDGER_CALLS="1",
                   FORGE_DUMP_LEDGER_NODES=",".join(map(str, ids)))
        r = subprocess.run([BIN], cwd=run, env=env, stdin=subprocess.DEVNULL, stdout=open(os.path.join(run, "forge_run.log"), "w"), stderr=subprocess.STDOUT)
        print(run, "exit", r.returncode, "印", len(ids))


def ledger(run):
    st, res = defaultdict(dict), defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(run, "ledger.csv"))):
        if r["call"] != "1":
            continue
        n, v = int(r["node"]), float(r["value"])
        if r["tag"] == "after_eos_bc":
            st[r["field"]][n] = v
        elif r["tag"] == "res_after_conv":
            res[r["field"]][n] = v
    faces = [r for r in csv.DictReader(open(os.path.join(run, "ledger.csv.faces"))) if r["call"] == "1"]
    return st, res, faces


def compare():
    xl, xh = coords(LOW), coords(HIGH)
    kl, kh = key(xl), key(xh)
    ih = {k: i for i, k in enumerate(kh)}
    lo_ids = marked(xl)
    to_h = {int(i): ih[kl[i]] for i in lo_ids}
    judge = [int(i) for i in lo_ids if 0.195 <= xl[i, 0] <= 0.230 and xl[i, 1] <= 1e-9 + (xl[i, 0] - 0.2) * np.tan(np.radians(10)) * (xl[i, 0] > 0.2)]
    # 1. 幾何: 双対体積
    with h5py.File(os.path.join(LOW, "ramp.h5")) as a, h5py.File(os.path.join(HIGH, "ramp.h5")) as b:
        vl, vh = a["CELLS/volume"][:], b["CELLS/volume"][:]
    dv = max(abs(vl[i] - vh[to_h[i]]) / vh[to_h[i]] for i in lo_ids)
    print(f"[幾何] 印 {len(lo_ids)} 節点の双対体積 最大相対差 {dv:.3e}")
    sl, rl, fl = ledger(LOW); sh, rh, fh = ledger(HIGH)
    # 2. 状態 (after_eos_bc)
    ds = {}
    for q in ("ro", "roUx", "roUy", "roUz", "roe", "P", "T"):
        if q in sl:
            ds[q] = max(abs(sl[q][i] - sh[q][to_h[i]]) / max(abs(sh[q][to_h[i]]), 1e-30) for i in lo_ids)
    print("[状態] after_eos_bc 最大相対差 " + "、".join(f"{q} {v:.2e}" for q, v in ds.items()))
    # 3. 面 (座標の組で照合)
    def fkey(r, xyz):
        return (tuple(np.round(xyz[int(r["ic0"])].astype(np.float64), 6)), tuple(np.round(xyz[int(r["ic1"])].astype(np.float64), 6)))
    FH = {fkey(r, xh): r for r in fh}
    cols = ["sx", "sy", "sz", "ss", "ro_L", "ro_R", "P_L", "P_R", "Ux_L", "Uy_L", "Ux_R", "Uy_R", "mdot", "F_ro", "F_roUx", "F_roUy", "F_roe"]
    worst = defaultdict(float); first = None; nmatch = 0; nmiss = 0
    for r in fl:
        k = fkey(r, xl)
        if k not in FH:
            k2 = (k[1], k[0])
            if k2 in FH:
                nmiss += 1   # 向きが逆 (比較は符号込みになるので別に数える)
            else:
                nmiss += 1
            continue
        h = FH[k]; nmatch += 1
        for c in cols:
            a, b = float(r[c]), float(h[c]); sc = max(abs(b), 1e-30)
            d = abs(a - b) / sc if abs(b) > 1e-12 else abs(a - b)
            if d > worst[c]:
                worst[c] = d
            if d > 1e-5 and first is None and c in ("sx", "sy", "sz", "ss", "ro_L", "ro_R", "P_L", "P_R"):
                first = (c, k, a, b)
    print(f"[面] 照合 {nmatch} 面 (向き違い・未対応 {nmiss})。最大相対差: " + "、".join(f"{c} {worst[c]:.1e}" for c in cols))
    if first:
        print(f"      最初に 1e-5 を超えた幾何/再構成: {first[0]} 面 {first[1]}: 低 {first[2]:.9g} / 高 {first[3]:.9g}")
    # 4. 残差 (判定)
    fabs = defaultdict(lambda: defaultdict(float))
    for r in fh:
        for c, q in (("F_ro", "res_ro"), ("F_roUx", "res_roUx"), ("F_roUy", "res_roUy"), ("F_roUz", "res_roUz"), ("F_roe", "res_roe")):
            for n in (int(r["ic0"]), int(r["ic1"])):
                fabs[n][q] += abs(float(r[c]))
    ok = True
    for q in ("res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"):
        w, wn = 0.0, None
        for i in judge:
            j = to_h[i]; sc = fabs[j][q] if fabs[j][q] > 0 else 1.0
            d = abs(rl[q][i] - rh[q][j]) / sc
            if d > w:
                w, wn = d, i
        ok &= w <= 1e-5
        print(f"[残差] {q:9s}: 角 {len(judge)} 節点の最大 |Δres|/Σ|F| {w:.3e}" + (f" (x {xl[wn,0]:.3f} y {xl[wn,1]:.4f} z {xl[wn,2]:.3f})" if wn is not None else ""))
    print(f"VERDICT: {'作用素差なし (≤ 1e-5) → 反復過程へ切り分けを移す' if ok else '作用素差あり (> 1e-5) → 最初に異なる幾何・再構成・流束へ絞る'}")


if __name__ == "__main__":
    if sys.argv[1] in ("setup", "compare"):
        {"setup": setup, "compare": compare}[sys.argv[1]]()


def fromC():
    """反復過程の切り分け: C の最終場から低領域 (A: slip / B: farfield) と C 自身を 6000 step 継続する"""
    xh = coords(SRC_HIGH)
    idx = {k: i for i, k in enumerate(key(xh))}
    with h5py.File(SNAP) as s:
        V = {q: s["VALUE"][q][:] for q in CONS}
    for dst, src in (("run_0094_v2c_B_fromC", "run_0064_v2c_B"), ("run_0095_v2c_A_fromC", "run_0063_v2c_A"), ("run_0096_v2c_C_cont", SRC_HIGH)):
        os.makedirs(dst)
        for f in ("solverConfig.yaml", "bcondConfig.yaml", "ramp.h5", "probe.yaml"):
            shutil.copy(os.path.join(src, f), dst)
        m = np.array([idx[k] for k in key(coords(dst))])
        with h5py.File(os.path.join(dst, "ramp.h5"), "r+") as f:
            for q in CONS:
                f["VALUE"][q][...] = V[q][m].astype(f["VALUE"][q].dtype)
        open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"{SNAP} の保存量を座標対応で写した ({len(m)} 節点)、設定は {src}\n")
    print("prepared")


if __name__ == "__main__" and sys.argv[1] == "fromC":
    fromC()
