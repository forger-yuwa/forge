#!/usr/bin/env python3
"""同一メッシュの restart: res_*.h5 の**保存量をそのまま**入力 h5 の /VALUE へ写す。

    python3 solver_density_cuda/tools/restart_field.py SRC_res.h5 DST_input.h5 [--dry-run]

**なぜ専用ツールが要るか** (2026-09-23 に実際にやらかした):
`res_*.h5` は保存量 (`ro,roUx,roUy,roUz,roe,roK,roOmega,…`) と原始量 (`P,T,Ux,…`) の**両方**を持つ。
同一メッシュの restart で原始量から保存量を**組み直すと**、`roUx = ro * Ux` の丸めで元の値に戻らない。
`case/56.gap_tp1187` で `interp_field.py` (cross-mesh 用) を同一メッシュに使ったところ、
`roUy` が **62104/65194 セルで不一致・最大 1.6 %** になった。測る量が
$\\dot m = |\\int \\rho U_y dx|$ だったので、そのまま回せば実験ごと壊れていた
(既存の `case/*/restart_field.py` も原始量から組み直す実装で、同じ罠を持つ)。

本ツールは**保存量を index コピー**する。DST/VALUE にある各データセットについて、
SRC/VALUE に同名があればその値をそのまま書き、無ければ DST の値を残す。
`wall_dist` は**常に DST のものを残す** (メッシュ由来の量)。
最後に「SRC と DST が保存量でビット一致すること」を検査し、しなければ**失敗させる**。

メッシュが違うとき (解像度変更・quad↔tri) は本ツールでなく `interp_field.py` を使う。
"""
import argparse, sys
import h5py
import numpy as np

KEEP_FROM_DST = {"wall_dist"}

ap = argparse.ArgumentParser()
ap.add_argument("src", help="継続元の res_*.h5")
ap.add_argument("dst", help="forge 入力 h5 (mesh.h5 等)。/VALUE を上書きする")
ap.add_argument("--dry-run", action="store_true", help="書かずに何が起きるかだけ出す")
a = ap.parse_args()

with h5py.File(a.src, "r") as s, h5py.File(a.dst, "r" if a.dry_run else "r+") as d:
    if "VALUE" not in s or "VALUE" not in d:
        sys.exit("SRC / DST のどちらかに /VALUE が無い")
    sv, dv = s["VALUE"], d["VALUE"]

    # --- 同一メッシュであることを確かめる (違えば interp_field.py を使うべき) ---
    n_dst = dv[next(iter(dv))].shape[0]
    n_src = sv[next(iter(sv))].shape[0]
    if n_src != n_dst:
        sys.exit(f"セル数が違う (SRC {n_src} != DST {n_dst})。"
                 " 同一メッシュでないので interp_field.py を使うこと")
    if "MESH/COORD" in s and "MESH/COORD" in d:
        if not np.array_equal(np.asarray(s["MESH/COORD"]), np.asarray(d["MESH/COORD"])):
            sys.exit("MESH/COORD が一致しない。同一メッシュでないので interp_field.py を使うこと")
        print(f"同一メッシュを確認 (MESH/COORD 一致, {n_dst} cells)")
    else:
        print(f"セル数一致 {n_dst} (MESH/COORD がどちらかに無いので座標比較は省略)")

    moved, kept, missing = [], [], []
    for name in sorted(dv):
        if name in KEEP_FROM_DST:
            kept.append(name); continue
        if name not in sv:
            missing.append(name); continue
        val = np.asarray(sv[name])
        if val.shape != dv[name].shape:
            sys.exit(f"{name}: 形が違う (SRC {val.shape} != DST {dv[name].shape})")
        if not a.dry_run:
            dv[name][...] = val.astype(dv[name].dtype, copy=False)
        moved.append(name)

    print(f"移した保存量  : {moved}")
    print(f"DST を残す    : {kept}")
    if missing:
        print(f"SRC に無く据置: {missing}")

    if a.dry_run:
        print("--dry-run: 何も書いていない")
        sys.exit(0)

    # --- 検査: 写した量が SRC とビット一致すること ---
    # **例外は「SRC の方が広い型」のとき** (倍精度 run の res_*.h5 → float32 の入力 h5)。
    # forge の入力 h5 は float32 なので、倍精度の場を種にするときは丸めが必ず入る
    # (run_0020_double 自身も float32 の種から出発している)。この場合はビット一致を求めず、
    # **丸めで失われた大きさを報告**して続行する。
    bad, narrowed = [], []
    for n in moved:
        a = np.asarray(sv[n]); b = np.asarray(dv[n])
        if np.array_equal(a, b):
            continue
        if a.dtype.itemsize > b.dtype.itemsize:
            rel = np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))) / max(np.max(np.abs(a)), 1e-300)
            narrowed.append((n, rel))
        else:
            bad.append(n)
    if bad:
        sys.exit(f"検査 NG: 写したのに SRC と一致しない: {bad}")
    if narrowed:
        print(f"型の縮小 ({sv[moved[0]].dtype} -> {dv[moved[0]].dtype}) で丸めが入った "
              f"(入力 h5 が float32 のため不可避):")
        for n, rel in narrowed:
            print(f"    {n:<10} 相対 {rel:.3e}")
    print(f"VERDICT: OK ({len(moved)} 量を移した"
          f"{'、うち ' + str(len(narrowed)) + ' 量は型の縮小で丸めあり' if narrowed else '、SRC とビット一致'})")
