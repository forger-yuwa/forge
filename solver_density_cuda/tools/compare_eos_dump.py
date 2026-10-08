#!/usr/bin/env python3
"""EOS 1 回の前後の全 cell 配列のダンプ (FORGE_DUMP_EOS_STEP / FORGE_DUMP_EOS_FILE) を 2 つ比べ、ビット一致かを判定する。

plans/active/tooling-sern-te-wake-grid.md §5.1 #2 の「凍結入力に対する EOS 1 回の全出力のビット比較」の比較器。
ダンプの形式は solver_density_cuda/eosDump.hpp (`/pre/<名前>`・`/post/<名前>` = var.c_d の全配列、nCells_all 長、
`/db/*` = 物性 DB のバイト列、ルートの属性 = EOS が読む設定・経路・solverConfig.yaml の本文)。

    python3 compare_eos_dump.py A.h5 B.h5 [--json out.json]

**幅・閾値は持たない。差分バイト数 0 だけが一致。** 判定:
  INVALID (試験不成立): 読めない・設定/経路/格子の大きさ/物性 DB の不一致・配列構造 (名前の集合・shape・dtype) の不一致・
      EOS の引数配列 (下の EOS_ARGS) の非有限・**EOS の入力 (EOS 直前の EOS 引数配列) の不一致**・
      書き込み先の列挙 (EOS_WRITES) の取りこぼし (ファイル内で EOS の前後に変わった配列が列挙に無い)・未対応の経路。
  DIFFERENT: 入力が一致したうえで、EOS 直後の EOS 引数配列に 1 バイトでも差がある (最初に異なる配列・節点を表示)。
  IDENTICAL: それ以外。

配列の分類はダンプを見る前に dependentVariables_d.cu から決めたもの (反復結果から選ばない):
  EOS_ARGS   = カーネルに渡る cell 配列 (読む・書く・書かずに残す) と化学種 roY{s}。
  EOS_WRITES = カーネル本体が書く配列 (TP: ro roe Ux Uy Uz T P Ht sonic k omega gamma cp Rmix / CPG 単相: gamma cp を除く)。
EOS が読まない配列 (勾配・残差など) は全部比べて差分バイト数・非有限の数を表示するが、判定には入れない
(カーネルの引数に無いので EOS の出力に影響し得ない。ファイル内で EOS の前後に変わらないことは取りこぼしの検査で確かめる)。
"""
import argparse
import json
import sys

import h5py
import numpy as np

# dependentVariables_d (cuda_forge/dependentVariables_d.cu) の引数 ro..Rmix の順
EOS_ARGS_BASE = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega",
                 "P", "Ht", "sonic", "k", "omega", "T", "Ux", "Uy", "Uz",
                 "gamma", "cp", "Rmix"]
# 本体が書く配列 (単相・凝縮なし; dependentVariables_d.cu の行は HEAD 6717f7f9): 速度 :88-90、TP 分岐 :230-253、CPG 単相分岐 :328-338、共通 k・ω :342-343。
# gamma・cp・Rmix・P・Ht・sonic・k・omega の入力値はカーネルが読まない (書くだけ)。読む入力は ro・roU・roe・roK・roOmega・T (温度反転の初期推定)・roY
EOS_WRITES_TP = ["ro", "roe", "Ux", "Uy", "Uz", "T", "P", "Ht", "sonic", "k", "omega", "gamma", "cp", "Rmix"]
# カーネル内で値が決まる順 (DIFFERENT のとき「最初に異なる配列」をこの順で示す)
EOS_ORDER = ["Ux", "Uy", "Uz", "T", "P", "ro", "roe", "Ht", "sonic", "gamma", "cp", "Rmix", "k", "omega"]
EOS_WRITES_CPG = ["ro", "roe", "Ux", "Uy", "Uz", "T", "P", "Ht", "sonic", "k", "omega", "Rmix"]


def nbytes_diff(a, b):
    return int(np.count_nonzero(a.view(np.uint8) != b.view(np.uint8)))


def nonfinite(a):
    if a.dtype.kind != "f":
        return 0
    return int(np.count_nonzero(~np.isfinite(a)))


def attrs_of(f):
    out = {}
    for k, v in f.attrs.items():
        if isinstance(v, bytes):
            v = v.decode("utf-8", "replace")
        elif isinstance(v, np.generic):
            v = v.item()
        out[k] = v
    return out


def fmt_val(x):
    if x.dtype == np.float32:
        return f"{float(x):.9g} (0x{int(x.view(np.uint32)):08x})"
    return repr(x)


def compare(pa, pb, max_nodes=5):
    rep = {"A": pa, "B": pb, "invalid": [], "different": [], "notes": [], "datasets": {}}
    inv = rep["invalid"]
    try:
        fa, fb = h5py.File(pa, "r"), h5py.File(pb, "r")
    except Exception as e:  # noqa: BLE001
        inv.append(f"読めない: {e}")
        return finish(rep)
    with fa, fb:
        for tag, f in (("A", fa), ("B", fb)):
            for g in ("pre", "post"):
                if g not in f:
                    inv.append(f"{tag}: /{g} が無い (EOS の後まで書けていない)")
        if inv:
            return finish(rep)
        ma, mb = attrs_of(fa), attrs_of(fb)
        rep["meta"] = {k: v for k, v in ma.items() if k != "solverConfig_yaml"}
        for k in sorted(set(ma) | set(mb)):
            if ma.get(k) != mb.get(k):
                inv.append(f"設定・経路の属性 {k} が不一致 (A {str(ma.get(k))[:80]!r} / B {str(mb.get(k))[:80]!r})")
        if ma.get("isImplicit") != 1 or ma.get("unsteady") != 0:
            inv.append(f"経路が定常陰解法でない (isImplicit {ma.get('isImplicit')}, unsteady {ma.get('unsteady')})")
        if ma.get("condensation", 0) != 0:
            inv.append("凝縮 (condensation 1) は未対応 (rog の書き込み・二相分岐の列挙をしていない)")
        tm = ma.get("thermalMethod")
        if tm not in (0, 2):
            inv.append(f"thermalMethod {tm} は未対応")
        nsr = int(ma.get("nSpeciesRegistered", 1) or 1)
        eos_args = EOS_ARGS_BASE + ([f"roY{s}" for s in range(nsr)] if nsr >= 2 else [])
        writes = EOS_WRITES_TP if tm == 2 else EOS_WRITES_CPG
        rep["eos_args"], rep["eos_writes"] = eos_args, writes
        n_all = ma.get("nCells_all")
        n_real = ma.get("nCells")

        # 物性 DB (EOS の入力)
        dbn_a = sorted(fa["db"].keys()) if "db" in fa else []
        dbn_b = sorted(fb["db"].keys()) if "db" in fb else []
        if dbn_a != dbn_b:
            inv.append(f"物性 DB のデータセットが不一致 ({dbn_a} / {dbn_b})")
        for k in dbn_a:
            if k in dbn_b:
                a, b = fa["db"][k][()], fb["db"][k][()]
                d = nbytes_diff(a, b) if a.shape == b.shape and a.dtype == b.dtype else -1
                rep["datasets"][f"db/{k}"] = {"diff_bytes": d, "nbytes": int(a.nbytes)}
                if d != 0:
                    inv.append(f"物性 DB {k} が不一致 (差分バイト {d})")

        # 配列構造
        names = {}
        for tag, f in (("A", fa), ("B", fb)):
            for g in ("pre", "post"):
                names[(tag, g)] = set(f[g].keys())
        allnames = set().union(*names.values())
        for key, s in names.items():
            miss = sorted(allnames - s)
            if miss:
                inv.append(f"{key[0]}:/{key[1]} に無い配列 {miss[:8]}{' …' if len(miss) > 8 else ''}")
        for nm in eos_args:
            if nm not in allnames:
                inv.append(f"EOS の引数配列 {nm} がダンプに無い")
        if inv:
            return finish(rep)

        changed = {"A": [], "B": []}
        tot = {"eos_args_post_diff_bytes": 0, "eos_args_pre_diff_bytes": 0, "other_pre_diff_bytes": 0,
               "other_post_diff_bytes": 0, "eos_args_bytes": 0, "all_bytes": 0}
        for nm in sorted(allnames):
            arr = {}
            for tag, f in (("A", fa), ("B", fb)):
                for g in ("pre", "post"):
                    arr[(tag, g)] = f[g][nm][()]
            shapes = {k: (v.shape, v.dtype.str) for k, v in arr.items()}
            if len(set(shapes.values())) != 1 or arr[("A", "pre")].shape != (n_all,):
                inv.append(f"{nm}: shape/dtype の不一致 {shapes} (nCells_all {n_all})")
                continue
            is_arg = nm in eos_args
            d = {
                "class": "eos_arg" if is_arg else "other",
                "dtype": arr[("A", "pre")].dtype.str,
                "pre_diff_bytes": nbytes_diff(arr[("A", "pre")], arr[("B", "pre")]),
                "post_diff_bytes": nbytes_diff(arr[("A", "post")], arr[("B", "post")]),
                "nonfinite": {f"{k[0]}:{k[1]}": nonfinite(v) for k, v in arr.items()},
            }
            for tag in ("A", "B"):
                if nbytes_diff(arr[(tag, "pre")], arr[(tag, "post")]) != 0:
                    changed[tag].append(nm)
            rep["datasets"][nm] = d
            nb = int(arr[("A", "post")].nbytes)
            tot["all_bytes"] += 2 * nb
            if is_arg:
                tot["eos_args_bytes"] += nb
                tot["eos_args_pre_diff_bytes"] += d["pre_diff_bytes"]
                tot["eos_args_post_diff_bytes"] += d["post_diff_bytes"]
                nf = sum(d["nonfinite"].values())
                if nf:
                    inv.append(f"EOS の引数配列 {nm} に非有限 {d['nonfinite']}")
                if d["pre_diff_bytes"]:
                    inv.append(f"EOS の入力 {nm} (EOS 直前) が不一致: 差分バイト {d['pre_diff_bytes']}")
                if d["post_diff_bytes"]:
                    a, b = arr[("A", "post")], arr[("B", "post")]
                    ua = a.view(np.uint8).reshape(len(a), -1)
                    ub = b.view(np.uint8).reshape(len(b), -1)
                    neq = np.nonzero((ua != ub).any(axis=1))[0]
                    inp = {k2: fa["pre"][k2][()] for k2 in ("ro", "roe", "T") if k2 in fa["pre"]}
                    info = {"array": nm, "n_elements": int(len(neq)), "diff_bytes": d["post_diff_bytes"],
                            "first_nodes": [{"node": int(i), "real": bool(i < n_real),
                                             "A": fmt_val(a[i]), "B": fmt_val(b[i]),
                                             "input_pre": {k2: fmt_val(v2[i]) for k2, v2 in inp.items()}}
                                            for i in neq[:max_nodes]]}
                    rep["different"].append(info)
            else:
                tot["other_pre_diff_bytes"] += d["pre_diff_bytes"]
                tot["other_post_diff_bytes"] += d["post_diff_bytes"]
        rep["totals"] = tot
        # 「最初に異なる配列」はカーネル内の計算順 (速度 → 温度 → 圧力 → 保存量の再構成 → Ht・音速 → 物性 → k・ω) で並べる
        rep["different"].sort(key=lambda x: (EOS_ORDER.index(x["array"]) if x["array"] in EOS_ORDER else len(EOS_ORDER), x["array"]))
        rep["changed_by_eos"] = changed
        for tag in ("A", "B"):
            gap = sorted(set(changed[tag]) - set(writes))
            if gap:
                inv.append(f"{tag}: EOS の前後で変わった配列が書き込み先の列挙に無い (取りこぼし) {gap}")
        rep["other_pre_mismatch"] = sorted(nm for nm, d in rep["datasets"].items()
                                           if isinstance(d, dict) and d.get("class") == "other" and d["pre_diff_bytes"])
        rep["other_nonfinite"] = sorted(nm for nm, d in rep["datasets"].items()
                                        if isinstance(d, dict) and d.get("class") == "other" and sum(d["nonfinite"].values()))
    return finish(rep)


def finish(rep):
    if rep["invalid"]:
        rep["verdict"] = "INVALID"
    elif rep["different"]:
        rep["verdict"] = "DIFFERENT"
    else:
        rep["verdict"] = "IDENTICAL"
    return rep


def print_report(rep):
    print(f"A: {rep['A']}\nB: {rep['B']}")
    m = rep.get("meta")
    if m:
        print(f"step {m.get('step')}  nCells {m.get('nCells')}  nCells_all {m.get('nCells_all')}  thermalMethod {m.get('thermalMethod')}"
              f"  nSpeciesRegistered {m.get('nSpeciesRegistered')}  thermoFloat {m.get('thermoFloat')}  discretization {m.get('discretization')}")
    ds = {k: v for k, v in rep.get("datasets", {}).items() if not k.startswith("db/")}
    if ds:
        n_arg = sum(1 for v in ds.values() if v["class"] == "eos_arg")
        print(f"配列 {len(ds)} 本 (EOS の引数 {n_arg} 本・その他 {len(ds) - n_arg} 本) × {{pre, post}}、物性 DB "
              f"{[k for k in rep['datasets'] if k.startswith('db/')]}")
        print("EOS の引数配列 (pre 差分バイト / post 差分バイト / 非有限 A:pre,A:post,B:pre,B:post):")
        for nm in rep["eos_args"]:
            d = ds.get(nm)
            if d:
                nf = d["nonfinite"]
                print(f"  {nm:10s} {d['pre_diff_bytes']:>10d} {d['post_diff_bytes']:>10d}   "
                      f"{nf['A:pre']},{nf['A:post']},{nf['B:pre']},{nf['B:post']}")
        t = rep["totals"]
        print(f"合計: EOS の引数 post 差分バイト {t['eos_args_post_diff_bytes']} / {t['eos_args_bytes']} バイト、"
              f"pre 差分バイト {t['eos_args_pre_diff_bytes']}; その他の配列 pre 差分 {t['other_pre_diff_bytes']}・post 差分 {t['other_post_diff_bytes']}"
              f" (判定外、EOS が読まない)")
        print(f"その他の配列で A/B の pre が異なるもの: {rep['other_pre_mismatch'] or 'なし'}")
        print(f"その他の配列で非有限を含むもの: {rep['other_nonfinite'] or 'なし'}")
        for tag in ("A", "B"):
            print(f"{tag}: EOS の前後で変わった配列 {rep['changed_by_eos'][tag]}")
        print(f"書き込み先の列挙 (dependentVariables_d.cu から): {rep['eos_writes']}")
    for d in rep["different"]:
        print(f"差: {d['array']} {d['n_elements']} 要素 / {d['diff_bytes']} バイト、最初の節点:")
        for n in d["first_nodes"]:
            print(f"    node {n['node']} ({'実節点' if n['real'] else 'ghost'}): A {n['A']}  B {n['B']}  入力 {n['input_pre']}")
    for s in rep["invalid"]:
        print(f"不成立: {s}")
    print(f"VERDICT: {rep['verdict']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--json", default=None, help="結果を JSON でも書く")
    ap.add_argument("--max-nodes", type=int, default=5)
    a = ap.parse_args()
    rep = compare(a.a, a.b, a.max_nodes)
    print_report(rep)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1, ensure_ascii=False, default=str)
    return {"IDENTICAL": 0, "DIFFERENT": 1, "INVALID": 2}[rep["verdict"]]


if __name__ == "__main__":
    sys.exit(main())
