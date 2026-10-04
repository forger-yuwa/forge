#!/usr/bin/env python3
"""case/50 の URANS (dual-time) run を、定常 run の収束場から作る。

定常 (局所 dt) の残差がリミットサイクルになる = 定常解が無い可能性 → 時間精度計算で平均±振幅を出す。
物理 dt は「開口せん断層 (f ~ St·U/w) を 20-40 分割」から決める。物理 CFL ≲ 12 を守る (memory)。

usage: python3 tools/make_urans.py --run run_0007_T1_urans --from run_0006_T1_wd0063_long [--dt 5e-8] [--steps 5000]
"""
import argparse, json, shutil, subprocess, os, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--from", dest="src", required=True)
    ap.add_argument("--dt", type=float, default=5e-8)
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--subiter", type=int, default=20)
    ap.add_argument("--out-int", type=int, default=50)
    ap.add_argument("--cfl", type=float, default=2.0)
    a = ap.parse_args()

    src, rd = CASE / a.src, CASE / a.run
    rd.mkdir(exist_ok=True)
    for f in ("bcondConfig.yaml", "species_db.yaml", "probe.yaml", "case_setup.json"):
        shutil.copy(src / f, rd / f)
    # 最終場を IC にする
    shutil.copy(src / "mesh.h5", rd / "mesh.h5")
    (rd / "CONTINUED_FROM").write_text(f"{a.src} (定常がリミットサイクル → URANS)\n")

    cfg = (src / "solverConfig.yaml").read_text()
    cfg = cfg.replace("  unsteady: 0", "  unsteady: 1").replace("  dualTime: 0", "  dualTime: 1")
    import re
    cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {a.steps}", cfg)
    cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {a.out_int}", cfg)
    cfg = re.sub(r"control: 1, dt: [^,]+,", f"control: 0, dt: {a.dt:.3e},", cfg)
    cfg = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {a.cfl}, cfl_pseudo: {a.cfl}", cfg)
    cfg = cfg.replace("  nStepInner: 4", f"  nStepInner: 4\n  bdfOrder: 2\n  nSubIterDualTime: {a.subiter}")
    (rd / "solverConfig.yaml").write_text(cfg)
    res = sorted(src.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if res:
        r = subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]), str(rd / "mesh.h5")],
                           env=ENV, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-1500:], r.stderr[-800:]); raise SystemExit("interp_field failed")
        print(f"  IC = {res[-1].name}")

    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    w = 1.270e-3
    f_sl = 0.4 * setup["U_inf"] / w
    print(f"URANS: dt = {a.dt:.2e} s, {a.steps} step → 物理時間 {a.dt*a.steps*1e6:.1f} µs")
    print(f"  開口せん断層 f ≈ {f_sl*1e-3:.0f} kHz (St 0.4) → 1 周期 {1/f_sl*1e6:.2f} µs = {1/f_sl/a.dt:.0f} step")
    print(f"  深さ方向 1/4 波長 f ≈ {340/(4*20.32e-3)*1e-3:.1f} kHz → 1 周期 {4*20.32e-3/340*1e6:.0f} µs")
    print(f"  物理 CFL ≈ {setup['U_inf']*a.dt/21.2e-6:.1f} (開口セル 21.2 µm)")
    print(f"  → {rd}")


if __name__ == "__main__":
    main()
