#!/usr/bin/env python3
"""forge の入力 HDF5 (変換器の出力) 2 つを、データセットごとにビット単位で比べる。

使い方:
    python3 compare_mesh_h5.py A.h5 B.h5 [--geometry-only] [--require-float64] [--verbose] [--show N]

plans/active/architecture-float-state-double-geometry.md §6.14 の比較器 (変換器の幾何の精度)。

- 対象: /MESH/* (/MESH/RENUMBER_PERM があればそれも)・/PLANES/*・/CELLS/*・/DUAL/*・/BCONDS/*・/VIZMESH/*・
  /VALUE/wall_dist。初期値 (/VALUE/ の他の量) は対象外。
- データセットごとに型 (dtype、バイト順を含む)・形・値をビット単位で比べる (符号付きゼロも区別する)。
  許容差は持たない。
- 非有限 (NaN/Inf) は、両側に同じ NaN があっても FAIL (浮動小数点のデータセットすべて)。
- 欠損 (片側にしか無い名前)・群とデータセットの取り違え・形の不一致は FAIL。
- 属性も、名前の集合・型・値を比べる (/MESH の nCells 等、/BCONDS/<id> の bcondKind 等)。
- `--geometry-only`: float の初期値・境界条件の実数値 (/BCONDS/<id>/VALUE/* の浮動小数点) の違いを許す
  (float のビルドと FP64 のビルドの変換器の比較、§6.14 の 2)。これらも名前の集合・形・有限性は確かめる。
  整数の表 (接続・番号・境界の所属・可視化の接続) と幾何は必ずビット単位で比べる。
- `--require-float64`: 幾何の浮動小数点のデータセット (/BCONDS/*/VALUE/* 以外) が両側とも binary64 で
  なければ FAIL (§6.14 の 2「幾何のデータセットが binary64」)。

変換器の HDF5 はバイト単位では再現しない (同じバイナリでも md5 が違う) ので、比べるのはデータセットの値。
終了コード: 0 = PASS、1 = FAIL、2 = 入力エラー (開けない等)。
"""
import argparse
import sys

import numpy as np
import h5py

SCOPE_GROUPS = ("MESH", "PLANES", "CELLS", "DUAL", "BCONDS", "VIZMESH")
SCOPE_DATASETS = ("VALUE/wall_dist",)
# 必要なデータセット (どちらの側にも無ければ欠損として FAIL。片側だけなら名前の集合の比較で FAIL)。
# /DUAL は node の置換後は書かれず、/VIZMESH・/MESH/RENUMBER_PERM は条件付きなので必須にしない。
REQUIRED = ("MESH/COORD", "MESH/CONNE", "PLANES/STRUCT", "PLANES/surfVect", "PLANES/surfArea", "PLANES/centCoords",
            "CELLS/STRUCT", "CELLS/volume", "CELLS/centCoords", "VALUE/wall_dist")


def in_scope(path):
    """path (先頭の / なし) が比較の対象か。"""
    if path in SCOPE_DATASETS:
        return True
    top = path.split("/", 1)[0]
    return top in SCOPE_GROUPS


def is_bc_value(path):
    """境界条件の実数値 (/BCONDS/<id>/VALUE/<name>) か。"""
    parts = path.split("/")
    return len(parts) == 4 and parts[0] == "BCONDS" and parts[2] == "VALUE"


def collect(h5):
    """対象の群とデータセットを {path: obj} で集める (path は先頭の / なし)。"""
    objs = {}

    def visit(name, obj):
        if in_scope(name):
            objs[name] = obj

    h5.visititems(visit)
    return objs


def kind_of(obj):
    if isinstance(obj, h5py.Dataset):
        return "dataset"
    if isinstance(obj, h5py.Group):
        return "group"
    return type(obj).__name__


def bits_view(a):
    """数値配列をビット列 (同じ幅の符号なし整数) として見る。比べられない型なら None。"""
    a = np.ascontiguousarray(np.atleast_1d(a))
    if a.dtype.kind in "fiub" and a.dtype.itemsize in (1, 2, 4, 8):
        return a.view(np.dtype(f"u{a.dtype.itemsize}"))
    return None


def fmt_elem(a, i):
    v = np.atleast_1d(a).ravel()[i]
    if np.asarray(v).dtype.kind == "f":
        bits = np.asarray(v).view(np.dtype(f"u{np.asarray(v).dtype.itemsize}"))
        return f"{float(v)!r} (0x{int(bits):0{2 * np.asarray(v).dtype.itemsize}x})"
    return repr(v.item() if hasattr(v, "item") else v)


def nonfinite_count(a):
    a = np.asarray(a)
    if a.dtype.kind in "fc":
        return int(np.count_nonzero(~np.isfinite(a)))
    return 0


def compare_attrs(path, oa, ob, fails):
    """属性の名前の集合・型・値を比べる。比べた属性の数を返す。"""
    na, nb = set(oa.attrs.keys()), set(ob.attrs.keys())
    for k in sorted(na - nb):
        fails.append(f"{path} @{k}: attribute only in A")
    for k in sorted(nb - na):
        fails.append(f"{path} @{k}: attribute only in B")
    n = 0
    for k in sorted(na & nb):
        n += 1
        ta = oa.attrs.get_id(k).dtype
        tb = ob.attrs.get_id(k).dtype
        va, vb = oa.attrs[k], ob.attrs[k]
        if ta != tb:
            fails.append(f"{path} @{k}: attribute dtype differs (A {ta}, B {tb})")
            continue
        aa, ab = np.asarray(va), np.asarray(vb)
        if aa.shape != ab.shape:
            fails.append(f"{path} @{k}: attribute shape differs (A {aa.shape}, B {ab.shape})")
            continue
        if nonfinite_count(aa) or nonfinite_count(ab):
            fails.append(f"{path} @{k}: attribute has non-finite values")
            continue
        ba, bb = bits_view(aa), bits_view(ab)
        if ba is not None and bb is not None:
            same = bool(np.array_equal(ba, bb))
        elif aa.dtype.kind == "O" or ab.dtype.kind == "O":
            same = bool(np.all(aa == ab))
        else:
            same = aa.dtype == ab.dtype and aa.tobytes() == ab.tobytes()
        if not same:
            fails.append(f"{path} @{k}: attribute value differs (A {va!r}, B {vb!r})")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--geometry-only", action="store_true",
                    help="境界条件の実数値 (/BCONDS/*/VALUE/* の浮動小数点) の型・値の違いを許す")
    ap.add_argument("--require-float64", action="store_true",
                    help="幾何の浮動小数点のデータセットが両側とも binary64 でなければ FAIL")
    ap.add_argument("--verbose", action="store_true", help="一致したデータセットも型・形つきで並べる")
    ap.add_argument("--show", type=int, default=3, help="不一致の要素を最初から何個まで出すか")
    args = ap.parse_args()

    try:
        fa = h5py.File(args.a, "r")
        fb = h5py.File(args.b, "r")
    except OSError as e:
        print(f"ERROR: cannot open: {e}", file=sys.stderr)
        return 2

    mode = "geometry-only" if args.geometry_only else "full"
    print(f"compare_mesh_h5: mode={mode}{' +require-float64' if args.require_float64 else ''}")
    print(f"  A = {args.a}")
    print(f"  B = {args.b}")

    oa, ob = collect(fa), collect(fb)
    fails = []
    n_bit = n_relaxed = n_attr = 0
    rows = []

    for path in sorted(set(oa) - set(ob)):
        fails.append(f"{path}: only in A ({kind_of(oa[path])})")
    for path in sorted(set(ob) - set(oa)):
        fails.append(f"{path}: only in B ({kind_of(ob[path])})")
    for path in REQUIRED:
        if path not in oa and path not in ob:
            fails.append(f"{path}: required dataset missing in both files")

    for path in sorted(set(oa) & set(ob)):
        xa, xb = oa[path], ob[path]
        ka, kb = kind_of(xa), kind_of(xb)
        if ka != kb:
            fails.append(f"{path}: kind differs (A {ka}, B {kb})")
            continue
        n_attr += compare_attrs(path, xa, xb, fails)
        if ka != "dataset":
            continue

        a, b = xa[()], xb[()]
        a, b = np.asarray(a), np.asarray(b)
        relaxed = args.geometry_only and is_bc_value(path) and a.dtype.kind == "f" and b.dtype.kind == "f"
        geomfloat = (not is_bc_value(path)) and (a.dtype.kind == "f" or b.dtype.kind == "f")

        nf_a, nf_b = nonfinite_count(a), nonfinite_count(b)
        if nf_a or nf_b:
            fails.append(f"{path}: non-finite values (A {nf_a}, B {nf_b})")
        if args.require_float64 and geomfloat and (a.dtype != np.dtype("<f8") or b.dtype != np.dtype("<f8")):
            fails.append(f"{path}: not binary64 (A {a.dtype.str}, B {b.dtype.str})")
        if a.shape != b.shape:
            fails.append(f"{path}: shape differs (A {a.shape}, B {b.shape})")
            continue
        if relaxed:
            n_relaxed += 1
            rows.append((path, f"A {a.dtype.str} / B {b.dtype.str}", a.shape, "shape+finite only (BC value)"))
            continue
        n_bit += 1
        if a.dtype != b.dtype:
            fails.append(f"{path}: dtype differs (A {a.dtype.str}, B {b.dtype.str})")
            rows.append((path, f"A {a.dtype.str} / B {b.dtype.str}", a.shape, "DTYPE DIFFERS"))
            continue
        ba, bb = bits_view(a), bits_view(b)
        if ba is None or bb is None:
            same = a.tobytes() == b.tobytes()
            if not same:
                fails.append(f"{path}: values differ (non-numeric bytes)")
            rows.append((path, a.dtype.str, a.shape, "bit-identical" if same else "DIFFERS"))
            continue
        diff = np.flatnonzero(ba.ravel() != bb.ravel())
        if diff.size:
            shown = ", ".join(f"[{i}] A={fmt_elem(a, i)} B={fmt_elem(b, i)}" for i in diff[: args.show])
            fails.append(f"{path}: {diff.size} / {ba.size} elements differ in bits; first: {shown}")
            rows.append((path, a.dtype.str, a.shape, f"DIFFERS {diff.size}/{ba.size}"))
        else:
            rows.append((path, a.dtype.str, a.shape, "bit-identical"))

    if args.verbose:
        print("datasets:")
        for path, dt, shape, st in rows:
            print(f"  {path:40s} {dt:16s} {str(shape):14s} {st}")
    print(f"summary: {n_bit} datasets compared bit-for-bit, {n_relaxed} BC-value datasets checked for "
          f"shape/finiteness only, {n_attr} attributes compared, {len(fails)} failures")
    for f in fails:
        print(f"  FAIL {f}")
    verdict = "PASS" if not fails else "FAIL"
    print(f"VERDICT: {verdict}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
