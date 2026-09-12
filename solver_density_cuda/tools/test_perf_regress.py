#!/usr/bin/env python3
"""perf_regress.py cmp のデータ検査の退行試験 (codex result-4 M1, 2026-09-12)。

小さな合成 res_1.h5 を作り、次を確認する:
  1. 正常データ: PASS (終了 0)
  2. ノイズ用ラベルから T が欠落: 数値判定に進まず終了 2 (旧実装は noise=inf で新版の T が 100 倍でも PASS した)
  3. 基準の ro の形状が新版・double 対照と違う: 終了 2 (旧実装は truth-dist base=inf で SHAPE を合格にした)
  4. 新版に NaN: 終了 2
使い方: python3 solver_density_cuda/tools/test_perf_regress.py
"""
import os, subprocess, sys, tempfile
import h5py, numpy as np

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "perf_regress.py")
FIELDS = ("ro", "P", "T", "Ux", "Uy", "Uz", "roe", "k", "omega")

def write(dst, label, n=8, scale=1.0, drop=(), nan=(), shape_override=None):
    d = os.path.join(dst, label); os.makedirs(d, exist_ok=True)
    with h5py.File(os.path.join(d, "res_1.h5"), "w") as f:
        g = f.create_group("VALUE")
        for i, k in enumerate(FIELDS):
            if k in drop: continue
            m = shape_override.get(k, n) if shape_override else n
            v = (np.linspace(1.0, 2.0, m) * (10.0 + i)).astype(np.float32) * scale
            if k in nan: v[0] = np.nan
            g.create_dataset(k, data=v)

def run(dst, *args):
    r = subprocess.run([sys.executable, TOOL, "cmp", dst] + list(args), capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr

def main():
    fails = 0
    with tempfile.TemporaryDirectory() as dst:
        open(os.path.join(dst, "solverConfig.yaml"), "w").write("time: {last: {nStepOuter: 1}}\n")
        write(dst, "base_r1"); write(dst, "base_r2"); write(dst, "dbl"); write(dst, "new")
        rc, out = run(dst, "base_r1", "new", "--noise", "base_r2", "--truth", "dbl")
        ok = (rc == 0 and "PASS" in out); print("1 normal           :", "ok" if ok else "FAIL rc=%d\n%s" % (rc, out)); fails += (not ok)
        # 2. ノイズラベルの T 欠落 + 新版の T が 100 倍
        write(dst, "noise_noT", drop=("T",)); write(dst, "new_bigT", scale=1.0); 
        with h5py.File(os.path.join(dst, "new_bigT", "res_1.h5"), "r+") as f: f["VALUE/T"][...] = f["VALUE/T"][...] * 100
        rc, out = run(dst, "base_r1", "new_bigT", "--noise", "noise_noT")
        ok = (rc == 2 and "DATA CHECK FAILED" in out and "MISSING field T" in out); print("2 noise missing T  :", "ok" if ok else "FAIL rc=%d\n%s" % (rc, out)); fails += (not ok)
        # 3. 基準の ro の形状が違う (新版・double 対照は同じ)
        write(dst, "ref_shape", shape_override={"ro": 9})
        rc, out = run(dst, "ref_shape", "new", "--noise", "base_r2", "--truth", "dbl")
        ok = (rc == 2 and "SHAPE ro" in out); print("3 ref shape ro     :", "ok" if ok else "FAIL rc=%d\n%s" % (rc, out)); fails += (not ok)
        # 4. 新版に NaN
        write(dst, "new_nan", nan=("P",))
        rc, out = run(dst, "base_r1", "new_nan", "--noise", "base_r2")
        ok = (rc == 2 and "NONFINITE P" in out); print("4 new NaN          :", "ok" if ok else "FAIL rc=%d\n%s" % (rc, out)); fails += (not ok)
    print("VERDICT:", "PASS" if fails == 0 else "FAIL (%d)" % fails); sys.exit(1 if fails else 0)

if __name__ == "__main__": main()
