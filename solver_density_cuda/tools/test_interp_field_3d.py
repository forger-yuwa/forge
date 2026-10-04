#!/usr/bin/env python3
"""interp_field.py の照合座標の回帰試験 (2026-09-27 修正: 3D で z を無視していた)。
合成 h5 (MESH/COORD + VALUE) で (1) 3D→3D は z まで見て正しい節点を拾う、(2) 2D (z 2 層) → 3D は x,y で押し出す、
(3) 3D→2D は拒否、を確かめる。化学種の照合は --force-species で外す。
  python3 solver_density_cuda/tools/test_interp_field_3d.py
"""
import os, subprocess, sys, tempfile
import h5py, numpy as np
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "interp_field.py")
CONS = ["ro", "roUx", "roUy", "roUz", "roe"]


def write(p, coords, vals):
    with h5py.File(p, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float32).reshape(-1)
        for k in CONS:
            f["VALUE/" + k] = vals[k].astype(np.float32)


def grid(nz):
    x, y, z = np.meshgrid(np.arange(3.0), np.arange(2.0), np.arange(nz, dtype=float), indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], 1)


def run(src, dst):
    return subprocess.run([sys.executable, T, src, dst, "--force-species"], capture_output=True, text=True)


fail = 0
def check(name, ok, detail=""):
    global fail
    print(("ok   " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    fail += 0 if ok else 1


with tempfile.TemporaryDirectory() as td:
    c3 = grid(4)
    v3 = {k: (1.0 + c3[:, 2] * 10.0 + c3[:, 0]) for k in CONS}   # z に強く依存
    # (1) 3D → 3D: DST は同じ節点を別順序で持つ
    perm = np.random.default_rng(0).permutation(len(c3))
    s, d = os.path.join(td, "s3.h5"), os.path.join(td, "d3.h5")
    write(s, c3, v3); write(d, c3[perm], {k: np.zeros(len(c3)) for k in CONS})
    r = run(s, d)
    with h5py.File(d, "r") as f:
        got = np.asarray(f["VALUE/ro"])
    check("3D→3D: z まで見て同じ節点の値を拾う", r.returncode == 0 and np.allclose(got, v3["ro"][perm]), r.stderr[-200:])
    # (2) 2D (z 2 層) → 3D: x,y で押し出し
    c2 = grid(2); v2 = {k: (1.0 + c2[:, 0] + 5.0 * c2[:, 1]) for k in CONS}
    s2 = os.path.join(td, "s2.h5"); d2 = os.path.join(td, "d23.h5")
    write(s2, c2, v2); write(d2, c3, {k: np.zeros(len(c3)) for k in CONS})
    r = run(s2, d2)
    with h5py.File(d2, "r") as f:
        got = np.asarray(f["VALUE/ro"])
    check("2D→3D: x,y で押し出す", r.returncode == 0 and np.allclose(got, 1.0 + c3[:, 0] + 5.0 * c3[:, 1]), r.stderr[-200:])
    # (3) 3D → 2D は拒否
    d32 = os.path.join(td, "d32.h5"); write(d32, c2, {k: np.zeros(len(c2)) for k in CONS})
    r = run(s, d32)
    check("3D→2D: 拒否する", r.returncode != 0 and "REFUSED" in (r.stdout + r.stderr))
print("ALL PASS" if fail == 0 else f"{fail} FAILED")
sys.exit(1 if fail else 0)
