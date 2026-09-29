#!/usr/bin/env python3
"""interp_field.py のトレーサ転送回帰 (codex result-5 M1): 元 res に原始量 Xi しか無いとき ρ·Xi を roXi として転送する。
usage: python3 tests/unit/test_interp_field_tracer.py  (tools/ の interp_field.py を subprocess で実行)"""
import os, subprocess, sys, tempfile
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parents[2] / "tools"
DB = '"N2":\n  MW: 0.0280134\n  nasa9_low: [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12, 710.846086, -10.76003744]\n  nasa9_high: [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86640027]\n'
CFG = 'physProp: {thermalMethod: 2, species: ["N2"], speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15, tracer: exhaust}\n'
fail = 0
def check(name, ok, info=""):
    global fail
    print(("[PASS] " if ok else "[FAIL] ") + name + (f"  [{info}]" if info else ""))
    fail += 0 if ok else 1
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    for sub in ("src", "dst"):
        (td / sub).mkdir(); (td / sub / "solverConfig.yaml").write_text(CFG); (td / sub / "species_db.yaml").write_text(DB)
    n = 8; xyz = np.column_stack([np.arange(n, dtype=float), np.zeros(n), np.zeros(n)])
    def mesh(f):
        f.create_dataset("MESH/COORD", data=xyz.ravel()); f.create_dataset("CELLS/centCoords", data=xyz.ravel())
    with h5py.File(td / "src" / "res_10.h5", "w") as f:      # res: 原始量のみ (Xi はあるが roXi は無い)
        mesh(f); V = f.create_group("VALUE"); ro = np.full(n, 0.8)
        V["ro"] = ro; V["P"] = np.full(n, 1e5); V["Ux"] = np.full(n, 10.0); V["Uy"] = np.zeros(n); V["Uz"] = np.zeros(n)
        V["roe"] = np.full(n, 2e5); V["Y0"] = np.ones(n); V["Xi"] = np.full(n, 0.3)
    with h5py.File(td / "dst" / "mesh.h5", "w") as f:        # 入力: roXi/ρ = 0.8 (転送で 0.3 に上書きされるべき)
        mesh(f); V = f.create_group("VALUE"); ro = np.ones(n)
        for k in ("ro", "roUx", "roUy", "roUz", "roe"): V[k] = ro.astype(np.float32)
        V["roY0"] = ro.astype(np.float32); V["roXi"] = (0.8 * ro).astype(np.float32)
    # 合成の res は化学種属性を持たない (未検証)。#3c 以降のツールは既定で停止するので、まずそれを確かめ、
    # 本題のトレーサ転送はその実行だけの許可 FORGE_ALLOW_UNVERIFIED_SPECIES=1 で回す (宛先は属性なしのまま)
    args = [sys.executable, str(HERE / "interp_field.py"), str(td / "src" / "res_10.h5"), str(td / "dst" / "mesh.h5")]
    env = dict(os.environ); env.pop("FORGE_ALLOW_UNVERIFIED_SPECIES", None)
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    with h5py.File(td / "dst" / "mesh.h5") as f:
        xi0 = f["VALUE/roXi"][:] / f["VALUE/ro"][:]
    check("interp_field (属性なしの SRC, 既定): 停止し宛先を書き換えない", r.returncode != 0 and "REFUSED" in r.stdout + r.stderr
          and np.allclose(xi0, 0.8, atol=1e-6), (r.stdout + r.stderr)[-200:].replace("\n", " | "))
    env["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    check("interp_field (res に Xi のみ, FORGE_ALLOW_UNVERIFIED_SPECIES=1): 終了コード 0", r.returncode == 0, (r.stdout + r.stderr)[-200:].replace("\n", " | "))
    with h5py.File(td / "dst" / "mesh.h5") as f:
        check("interp_field: 許可して通した宛先に化学種属性が無い", not any(k.startswith("species") for k in f.attrs), str(list(f.attrs)))
    with h5py.File(td / "dst" / "mesh.h5") as f:
        xi = f["VALUE/roXi"][:] / f["VALUE/ro"][:]
        check("interp_field: 宛先の roXi/ρ = 0.3 (ρ·Xi で転送)", np.allclose(xi, 0.3, atol=1e-6), f"{xi[0]:.6f}")
    check("interp_field: 転送一覧に roXi", "roXi" in r.stdout)
print("ALL PASS" if fail == 0 else f"{fail} FAILED"); sys.exit(1 if fail else 0)
