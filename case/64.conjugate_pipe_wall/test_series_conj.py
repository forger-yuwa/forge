#!/usr/bin/env python3
r"""series_conj.py の負例試験 (2026-09-30 result レビュー M3)。test_eval_conj の合成 run を 8 スナップショットに複製し、
正常は PASS、NaN・壁節点の欠落・重複・尺度 0 は REFUSED になることを見る。作業ディレクトリに実ファイルで作る。"""
from __future__ import annotations
import argparse, shutil, subprocess, sys, tempfile
from pathlib import Path
import h5py
import numpy as np
import test_eval_conj as te

HERE = Path(__file__).resolve().parent


def multi(d: Path, mod=None):
    te.synth(d)
    for st in range(200, 900, 100):
        for stem in ("res", "res_wall_3", "res_solid_3"):
            shutil.copy(d / f"{stem}_100.h5", d / f"{stem}_{st}.h5")
    last = d / "res_wall_3_800.h5"
    if mod in ("nan", "missing", "dup", "zero"):
        with h5py.File(last, "a") as w:
            q = w["VALUE/iface_q_eff"][:]; T = w["VALUE/iface_Tw_bc"][:]; c = w["MESH/COORD"][:].reshape(-1, 3); ok = w["VALUE/iface_ok"][:]
            if mod == "nan": q[len(q) // 2] = np.nan
            if mod == "missing": q, T, c, ok = q[:-1], T[:-1], c[:-1], ok[:-1]
            if mod == "dup": c[1] = c[0]
            if mod == "zero": T[:] = 300.0
            for k in ("VALUE/iface_q_eff", "VALUE/iface_Tw_bc", "MESH/COORD", "VALUE/iface_ok"): del w[k]
            w.create_dataset("VALUE/iface_q_eff", data=q); w.create_dataset("VALUE/iface_Tw_bc", data=T)
            w.create_dataset("MESH/COORD", data=c.ravel()); w.create_dataset("VALUE/iface_ok", data=ok)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work"); a = ap.parse_args()
    work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix="seriesconj_"))
    cases = [("正常 (8 枚同一)", None, "PASS", 0), ("NaN", "nan", "REFUSED", 2), ("壁節点の欠落", "missing", "REFUSED", 2),
             ("壁節点の重複", "dup", "REFUSED", 2), ("尺度 0 (壁温 = T_in)", "zero", "REFUSED", 2)]
    bad = 0
    for nm, mod, want, rcw in cases:
        d = work / (mod or "ok"); multi(d, mod)
        p = subprocess.run([sys.executable, str(HERE / "series_conj.py"), "A", str(d)], capture_output=True, text=True)
        v = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("VERDICT")]
        w = v[-1].split(":")[1].strip().split()[0] if v else "(no VERDICT)"
        good = (w == want and p.returncode == rcw); bad += 0 if good else 1
        print(f"  {'ok ' if good else 'NG '} {nm:<22} 判定 {w:<8} rc {p.returncode} (期待 {want}/{rcw})")
        if not good: print("     " + "\n     ".join((p.stdout + p.stderr).splitlines()[-6:]))
    print(f"VERDICT: {'PASS' if bad == 0 else 'FAIL'}  ({len(cases) - bad}/{len(cases)})")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
