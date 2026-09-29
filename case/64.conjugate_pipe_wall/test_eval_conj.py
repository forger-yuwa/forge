#!/usr/bin/env python3
r"""eval_conj.py の負例試験 (A、N_r=32 の格子、Poiseuille の合成場)。参照解そのものを forge の出力として書いた合成 run が PASS し、
改変した run が FAIL / REFUSED になることを見る。合成 run は作業ディレクトリに実ファイルで作る (元の run には触らない)。

    python3 test_eval_conj.py [--work DIR]
"""
from __future__ import annotations
import argparse, shutil, subprocess, sys, tempfile
from pathlib import Path
import h5py
import numpy as np
import pipe_common as pc
from pipe_common import gc
import eval_conj as ev

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh" / "pipe_r32.h5"; DRY = HERE / "run_0002_dry_r32" / "res_wall_3_1.h5"
SOLID = HERE / "mesh" / "solid_run_0002_dry_r32_A1.h5"


def synth(d: Path, mod=None):
    d.mkdir(parents=True)
    shutil.copy(MESH, d / "mesh.h5"); shutil.copy(SOLID, d / "solid.h5")
    (d / "RUN_INPUTS.txt").write_text("kind cht A1\n")
    with h5py.File(MESH, "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    x, r = xyz[:, 0], xyz[:, 1]
    u = 2 * gc.U_M * np.clip(1 - (r / pc.R) ** 2, 0, None)
    p = gc.P_OUT + 8 * gc.MU * gc.U_M / pc.R ** 2 * (pc.L_HEAT + pc.L_DOWN - x)
    with h5py.File(d / "res_100.h5", "w") as h:
        for k, v in (("ro", p / (gc.RGAS * gc.T_IN)), ("Ux", u), ("Uy", 0 * u), ("P", p), ("T", np.full_like(u, gc.T_IN))):
            h.create_dataset(f"VALUE/{k}", data=v)
    pc.CASE = "A1"
    xs_f, ys_f, G, _ = ev.load_fields_only(d, 100)
    P, xs, ys, jw, k_s = ev.build("A", xs_f, ys_f, G, 2, pc=pc)
    T = ev._solve_rowwise(P, P.T_in_profile); q = ev.interface_q(P, T, jw, k_s, "A", 2)
    with h5py.File(DRY, "r") as w:
        c = np.asarray(w["MESH/COORD"][:], float); conne = np.asarray(w["MESH/CONNE"][:])
    cx = c.reshape(-1, 3)[:, 0]
    idx = [int(np.argmin(np.abs(xs - xv))) for xv in cx]
    Tw, qw, ok = T[idx, jw].copy(), q[idx].copy(), np.ones(len(cx))
    if mod == "signflip": qw = -qw
    if mod == "Tshift": Tw = Tw + 1.0
    if mod == "ok0": ok[len(ok) // 2] = 0
    with h5py.File(d / "res_wall_3_100.h5", "w") as w:
        w.create_dataset("MESH/COORD", data=c); w.create_dataset("MESH/CONNE", data=conne)
        w.create_dataset("VALUE/iface_q_eff", data=-qw); w.create_dataset("VALUE/iface_Tw_bc", data=Tw); w.create_dataset("VALUE/iface_ok", data=ok)


def run_eval(d):
    p = subprocess.run([sys.executable, str(HERE / "eval_conj.py"), "A", str(d)], capture_output=True, text=True)
    v = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("VERDICT")]
    return (v[-1].split(":")[1].strip().split()[0] if v else "(no VERDICT)"), p.returncode, p.stdout + p.stderr


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work"); a = ap.parse_args()
    work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix="evalconj_"))
    cases = [("正常 (参照解そのもの)", None, "PASS", 0), ("熱流束の符号反転", "signflip", "FAIL", 1),
             ("界面温度 +1 K", "Tshift", "FAIL", 1), ("iface_ok = 0 の節点", "ok0", "REFUSED", 2)]
    bad = 0
    for nm, mod, want, rc_want in cases:
        d = work / (mod or "ok"); synth(d, mod)
        w, rc, out = run_eval(d)
        good = (w == want and rc == rc_want); bad += 0 if good else 1
        print(f"  {'ok ' if good else 'NG '} {nm:<22} 判定 {w:<8} rc {rc} (期待 {want}/{rc_want})")
        if not good: print("     " + "\n     ".join(out.splitlines()[-9:]))
    print(f"VERDICT: {'PASS' if bad == 0 else 'FAIL'}  ({len(cases) - bad}/{len(cases)})")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
