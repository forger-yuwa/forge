#!/usr/bin/env python3
"""EOS 1 回の前後の全 cell 配列のダンプ (FORGE_DUMP_EOS_STEP / FORGE_DUMP_EOS_FILE) を 2 つ比べ、ビット一致かを判定する。

plans/active/tooling-sern-te-wake-grid.md §5.1 #2 の「凍結入力に対する EOS 1 回の全出力のビット比較」の比較器。
ダンプの形式は solver_density_cuda/eosDump.hpp (`/pre/<名前>`・`/post/<名前>` = var.c_d の全配列、nCells_all 長、
`/db/*` = 物性 DB のバイト列、ルートの属性 = EOS が読む設定・経路・solverConfig.yaml の本文)。

    python3 compare_eos_dump.py A.h5 B.h5 [--json out.json]
    python3 compare_eos_dump.py --rejudge-json 保存済みの比較.json [--json out.json]   # ダンプを読まずに両規則で再判定

**幅・閾値は持たない。差分バイト数 0 だけが一致。** 判定は 2 本を並べて出す (終了コードと JSON の "verdict" は改訂規則 v2):

VERDICT_v1 (旧規則。2026-10-08 に登録した規則そのまま。記録として残す):
  INVALID: 読めない・設定/経路/格子の大きさ/物性 DB の不一致・配列構造 (名前の集合・shape・dtype) の不一致・
      **EOS の引数配列 (EOS_ARGS) の pre・post の非有限**・EOS の入力 (EOS 直前の EOS 引数配列) の不一致・
      書き込み先の列挙の取りこぼし (EOS の前後で変わった配列が列挙に無い)・未対応の経路。
  DIFFERENT: 入力が一致したうえで、EOS 直後の EOS 引数配列に 1 バイトでも差がある。
  IDENTICAL: それ以外。

VERDICT_v2 (改訂規則。codex diagnose notes/reviews/2026-10-08-eos-dump-rule-diagnose.md を採用):
  v1 から変えるのは**有限性の範囲だけ**で、「EOS が読む入力 (読む集合) の pre」と「EOS の全書き込み先 (書く集合) の post」に限る。
  書き込み専用の配列 (TP の gamma・cp など、カーネルが読まずに上書きする) の pre の非有限は判定に入れず、件数を情報として出す
  (無害という認定ではない: plan §5.1 #8 の別件)。維持するもの: 壁・ghost を含む全節点での EOS の引数配列の post の比較
  (差分 0 バイト)・EOS の引数配列の pre の同一性 (書き込み専用でも凍結入力の同一性として要求)・未変更入力の不変
  (EOS の前後で変わった配列 ⊂ 書く集合)・構造/設定/物性 DB の一致。

読み書きの集合は、ダンプを見る前に dependentVariables_d.cu (HEAD 30f86067) のカーネル本体から決めた定数 (反復結果から選ばない)。
カーネルは節点ごとに分岐する (TP / CPG、単相 / 二相) ので、設定で入り得る分岐の和集合を取る。行番号は dependentVariables_d.cu:
  全分岐に共通  読む ro roUx roUy roUz roe (:81 床事象の入力・:86 密度床・:88-90 速度・:93 内部エネルギー)・roK roOmega (:342-343)
                書く Ux Uy Uz (:88-90。書いた値を :92 の ek に読み直すので入力値は読まない)・ro (:232/:297/:309/:331)・k omega (:342-343)
  TP 単相       読む + roY{s} (組成 :108/:117、nSpeciesRegistered ≥ 2)・T (温度反転の初期推定 :126)
                書く + T P (:230-231) roe (:234) Ht (:236) sonic gamma cp Rmix (:250-253)
  TP 二相       読む + rog_{s} (:133)・roY[condGasSpecies] (:138/:156)。書くのは TP 単相と同じ (反転失敗の節点は roe を書かない :234)。
                condEquilibrium 2 は rog_0 も書き (:164)、wrapper が condensationPrimitive を続けて呼ぶ (:398-401)
  CPG 単相      読む 共通のみ (gamma・cp は設定のスカラー)。書く + T P (:328-329) roe (:332) Ht (:334) sonic (:336) Rmix (:338)。
                gamma・cp の配列は読みも書きもしない
  CPG 二相      読む + rog_{s} (:262)・T (初期推定 :278/:289)。書くのは CPG 単相と同じ (反転失敗の節点は ro・Rmix だけ :297-298)。
                condEquilibrium 2 は rog_0 も書く (:280)
凝縮 (condensation 1) は v1・v2 とも未対応 (INVALID): カーネルが読む凝縮の物性表と設定 (condTb・condOpts) がダンプに無く、
読む入力の同一性を確かめられない。二相の集合は定数として持つ (報告に出す) が、判定は未対応のまま。
EOS が読まない配列 (勾配・残差など) は全部比べて差分バイト数・非有限の数を表示するが、判定には入れない
(カーネルの引数に無いので EOS の出力に影響し得ない。EOS の前後で変わらないことは未変更入力の検査で確かめる)。
"""
import argparse
import json
import re
import sys

import h5py
import numpy as np

RULE_VERSION = 2

# dependentVariables_d (cuda_forge/dependentVariables_d.cu) の引数 ro..Rmix の順 (v1 の EOS_ARGS の基本部分)
EOS_ARGS_BASE = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega",
                 "P", "Ht", "sonic", "k", "omega", "T", "Ux", "Uy", "Uz",
                 "gamma", "cp", "Rmix"]
# v1 の書き込み先 (単相・凝縮なし)。v2 の単相の書く集合と同じ (test_compare_eos_dump.py で確かめる)
EOS_WRITES_TP = ["ro", "roe", "Ux", "Uy", "Uz", "T", "P", "Ht", "sonic", "k", "omega", "gamma", "cp", "Rmix"]
EOS_WRITES_CPG = ["ro", "roe", "Ux", "Uy", "Uz", "T", "P", "Ht", "sonic", "k", "omega", "Rmix"]
# カーネル内で値が決まる順 (DIFFERENT のとき「最初に異なる配列」をこの順で示す)
EOS_ORDER = ["Ux", "Uy", "Uz", "T", "P", "ro", "roe", "Ht", "sonic", "gamma", "cp", "Rmix", "k", "omega"]

# ---- v2 の読み書きの集合 (分岐ごと。根拠の行はモジュールの説明)。"roY*" = roY{s} (nSpeciesRegistered ≥ 2 のときの全種)、
#      "rog_*" = rog_{s} (nCondSpeciesRegistered 本) ----
EOS_READS_COMMON = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]
EOS_WRITES_COMMON = ["Ux", "Uy", "Uz", "ro", "k", "omega"]
EOS_BRANCHES = {
    "tp_1ph":  {"reads": EOS_READS_COMMON + ["roY*", "T"],
                "writes": EOS_WRITES_COMMON + ["T", "P", "roe", "Ht", "sonic", "gamma", "cp", "Rmix"]},
    "tp_2ph":  {"reads": EOS_READS_COMMON + ["roY*", "T", "rog_*"],
                "writes": EOS_WRITES_COMMON + ["T", "P", "roe", "Ht", "sonic", "gamma", "cp", "Rmix"]},
    "cpg_1ph": {"reads": list(EOS_READS_COMMON),
                "writes": EOS_WRITES_COMMON + ["T", "P", "roe", "Ht", "sonic", "Rmix"]},
    "cpg_2ph": {"reads": EOS_READS_COMMON + ["rog_*", "T"],
                "writes": EOS_WRITES_COMMON + ["T", "P", "roe", "Ht", "sonic", "Rmix"]},
}
# condEquilibrium 2 で二相の分岐が追加で書く配列 (wrapper の condensationPrimitive の書き込みは列挙していない = 未対応の理由)
EOS_WRITES_EQ2_EXTRA = ["rog_0"]

# 保存済みの比較 JSON (旧形式) の invalid の文言のうち、事実 (meta・datasets・changed_by_eos) から再計算できるもの。再判定でだけ使う
_RECOMPUTABLE = [r"^EOS の引数配列 \S+ に非有限", r"^EOS の入力 \S+ \(EOS 直前\) が不一致", r"取りこぼし",
                 r"^経路が定常陰解法でない", r"^凝縮 \(condensation 1\) は未対応", r"^thermalMethod \S+ は未対応"]


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


def _expand(names, nsr, ncond):
    out = []
    for nm in names:
        if nm == "roY*":
            out += [f"roY{s}" for s in range(nsr)] if nsr >= 2 else []
        elif nm == "rog_*":
            out += [f"rog_{s}" for s in range(ncond)]
        else:
            out.append(nm)
    return out


def eos_sets(meta):
    """設定から入り得る分岐と、その和集合の読む集合・書く集合を返す: (分岐のリスト, 読む, 書く, 未対応の理由のリスト)。"""
    tm = meta.get("thermalMethod")
    cond = int(meta.get("condensation", 0) or 0)
    nsr = int(meta.get("nSpeciesRegistered", 1) or 1)
    ncond = int(meta.get("nCondSpeciesRegistered", 0) or 0)
    unsup = []
    if tm == 2:
        br = ["tp_1ph"] + (["tp_2ph"] if cond == 1 else [])
    elif tm == 0:
        br = ["cpg_1ph"] + (["cpg_2ph"] if cond == 1 else [])
    else:
        br = []
        unsup.append(f"thermalMethod {tm} は未対応")
    reads, writes = [], []
    for b in br:
        reads += _expand(EOS_BRANCHES[b]["reads"], nsr, ncond)
        writes += _expand(EOS_BRANCHES[b]["writes"], nsr, ncond)
    if cond != 0:
        eq2 = (meta.get("condEquilibrium") == 2)
        if eq2:
            writes += EOS_WRITES_EQ2_EXTRA
        unsup.append("凝縮 (condensation 1) は未対応: カーネルが読む凝縮の物性表と設定 (condTb・condOpts) がダンプに無く、"
                     "読む入力の同一性を確かめられない" + ("。condEquilibrium 2 は wrapper が condensationPrimitive で rog・モーメント・"
                                                          "原始量も書く (列挙していない)" if eq2 else ""))
    return br, list(dict.fromkeys(reads)), list(dict.fromkeys(writes)), unsup


def v1_sets(meta):
    """旧規則の引数配列と書き込み先 (登録どおり)。"""
    tm = meta.get("thermalMethod")
    nsr = int(meta.get("nSpeciesRegistered", 1) or 1)
    args = EOS_ARGS_BASE + ([f"roY{s}" for s in range(nsr)] if nsr >= 2 else [])
    writes = EOS_WRITES_TP if tm == 2 else EOS_WRITES_CPG
    return args, writes


def judged_arrays(meta):
    """v2 で判定する配列 (= v1 の引数配列 + 読む集合のうち引数に無いもの)。"""
    args1 = v1_sets(meta)[0]
    return args1 + [r for r in eos_sets(meta)[1] if r not in args1]


def route_invalid(meta, rule):
    out = []
    if meta.get("isImplicit") != 1 or meta.get("unsteady") != 0:
        out.append(f"経路が定常陰解法でない (isImplicit {meta.get('isImplicit')}, unsteady {meta.get('unsteady')})")
    if rule == 1:
        if meta.get("condensation", 0) != 0:
            out.append("凝縮 (condensation 1) は未対応 (rog の書き込み・二相分岐の列挙をしていない)")
        tm = meta.get("thermalMethod")
        if tm not in (0, 2):
            out.append(f"thermalMethod {tm} は未対応")
    else:
        out += eos_sets(meta)[3]
    return out


def judge_v1(rep):
    """旧規則 (登録どおり)。rep の事実 (structural_invalid・meta・datasets・changed_by_eos・different) から判定する。"""
    inv = list(rep["structural_invalid"])
    meta = rep.get("meta") or {}
    diff = []
    if rep.get("facts_complete"):
        inv += route_invalid(meta, 1)
        args, writes = v1_sets(meta)
        ds = rep["datasets"]
        for nm in sorted(ds):
            if nm.startswith("db/") or nm not in args:
                continue
            d = ds[nm]
            if sum(d["nonfinite"].values()):
                inv.append(f"EOS の引数配列 {nm} に非有限 {d['nonfinite']}")
            if d["pre_diff_bytes"]:
                inv.append(f"EOS の入力 {nm} (EOS 直前) が不一致: 差分バイト {d['pre_diff_bytes']}")
        for tag in ("A", "B"):
            gap = sorted(set(rep["changed_by_eos"][tag]) - set(writes))
            if gap:
                inv.append(f"{tag}: EOS の前後で変わった配列が書き込み先の列挙に無い (取りこぼし) {gap}")
        diff = [d for d in rep.get("different", []) if d["array"] in args]
    return ("INVALID" if inv else ("DIFFERENT" if diff else "IDENTICAL")), inv


def judge_v2(rep):
    """改訂規則: 有限性は読む集合の pre と書く集合の post だけ。それ以外は v1 と同じ厳しさ。"""
    inv = list(rep["structural_invalid"])
    info = []
    meta = rep.get("meta") or {}
    diff = []
    if rep.get("facts_complete"):
        inv += route_invalid(meta, 2)
        _, reads, writes, _ = eos_sets(meta)
        args = judged_arrays(meta)
        ds = rep["datasets"]
        for nm in list(dict.fromkeys(reads + writes)):
            if nm not in ds:
                inv.append(f"EOS の読み書きの集合の配列 {nm} がダンプに無い")
        for nm in sorted(ds):
            if nm.startswith("db/") or nm not in args:
                continue
            d = ds[nm]
            nf = d["nonfinite"]
            if nm in reads and (nf["A:pre"] or nf["B:pre"]):
                inv.append(f"EOS が読む入力 {nm} の pre に非有限 (A {nf['A:pre']}, B {nf['B:pre']})")
            if nm in writes and (nf["A:post"] or nf["B:post"]):
                inv.append(f"EOS の書き込み先 {nm} の post に非有限 (A {nf['A:post']}, B {nf['B:post']})")
            if nm not in reads and (nf["A:pre"] or nf["B:pre"]):
                info.append(f"{'書き込み専用' if nm in writes else 'EOS が読みも書きもしない引数'} {nm} の pre に非有限 "
                            f"A {nf['A:pre']}・B {nf['B:pre']} (判定外。無害の認定ではない)")
            if nm not in writes and nm not in reads and (nf["A:post"] or nf["B:post"]):
                info.append(f"EOS が読みも書きもしない引数 {nm} の post に非有限 A {nf['A:post']}・B {nf['B:post']} (判定外)")
            if d["pre_diff_bytes"]:
                kind = "読む入力" if nm in reads else "引数配列 (書き込み専用・不使用でも凍結入力として同一を要求)"
                inv.append(f"EOS の{kind} {nm} の pre が不一致: 差分バイト {d['pre_diff_bytes']}")
        for tag in ("A", "B"):
            gap = sorted(set(rep["changed_by_eos"][tag]) - set(writes))
            if gap:
                inv.append(f"{tag}: EOS の前後で変わった配列が書く集合に無い (未変更入力の不変が破れた・取りこぼし) {gap}")
        diff = [d for d in rep.get("different", []) if d["array"] in args]
    return ("INVALID" if inv else ("DIFFERENT" if diff else "IDENTICAL")), inv, info


def compare(pa, pb, max_nodes=5):
    rep = {"A": pa, "B": pb, "rule_version": RULE_VERSION, "structural_invalid": [], "different": [], "notes": [],
           "datasets": {}, "facts_complete": False}
    inv = rep["structural_invalid"]
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
        args1, writes1 = v1_sets(ma)
        branches, reads2, writes2, _ = eos_sets(ma)
        rep["eos_args"], rep["eos_writes"] = args1, writes1
        rep["eos_branches"], rep["eos_reads_v2"], rep["eos_writes_v2"] = branches, reads2, writes2
        args = judged_arrays(ma)
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
        for nm in args1:
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
            is_arg = nm in args
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
        rep["other_pre_mismatch"] = sorted(nm for nm, d in rep["datasets"].items()
                                           if isinstance(d, dict) and d.get("class") == "other" and d["pre_diff_bytes"])
        rep["other_nonfinite"] = sorted(nm for nm, d in rep["datasets"].items()
                                        if isinstance(d, dict) and d.get("class") == "other" and sum(d["nonfinite"].values()))
        rep["facts_complete"] = not inv
    return finish(rep)


def finish(rep):
    rep["verdict_v1"], rep["invalid_v1"] = judge_v1(rep)
    rep["verdict_v2"], rep["invalid_v2"], rep["info_v2"] = judge_v2(rep)
    # 採用中の規則は v2 (終了コード・"verdict"・"invalid" も v2)。旧規則の判定は verdict_v1・invalid_v1 に残す
    rep["verdict"], rep["invalid"] = rep["verdict_v2"], rep["invalid_v2"]
    return rep


def rejudge_json(path):
    """保存済みの比較 JSON (rule_version の無い旧形式を含む) を、記録された事実から両規則で再判定する (ダンプは読まない)。

    使う事実: meta・datasets (配列ごとの非有限の数・pre/post の差分バイト)・changed_by_eos・different。
    旧形式の invalid の文言のうち事実から再計算できないもの (読めない・設定や構造の不一致など) は両規則にそのまま引き継ぐ。
    事実が欠けていれば両方「再判定不能」。v1 の再計算が保存された判定と食い違えば v2 も「再判定不能」にする (読み取りを疑う)。"""
    j = json.load(open(path))
    ver = j.get("rule_version", 1)
    out = {"source": path, "stored_verdict": j.get("verdict"), "stored_rule_version": ver}
    miss = [k for k in ("meta", "datasets", "changed_by_eos", "different") if k not in j]
    if "structural_invalid" in j:
        carried = list(j["structural_invalid"])
    else:
        carried = [s for s in j.get("invalid", []) if not any(re.search(p, s) for p in _RECOMPUTABLE)]
    rep = {"structural_invalid": carried, "meta": j.get("meta", {}), "datasets": j.get("datasets", {}),
           "changed_by_eos": j.get("changed_by_eos", {"A": [], "B": []}), "different": j.get("different", []),
           "facts_complete": not miss and not carried}
    if miss and not carried:
        out.update(verdict_v1="再判定不能", verdict_v2="再判定不能", reason=f"JSON に事実が無い {miss}")
        return out
    if rep["facts_complete"]:
        lack = [nm for nm in judged_arrays(rep["meta"])
                if nm not in rep["datasets"] or "nonfinite" not in rep["datasets"][nm]]
        if lack:
            out.update(verdict_v1="再判定不能", verdict_v2="再判定不能", reason=f"配列の事実が欠ける {lack[:8]}")
            return out
    v1, inv1 = judge_v1(rep)
    v2, inv2, info2 = judge_v2(rep)
    br, reads, writes, _ = eos_sets(rep["meta"])
    out.update(verdict_v1=v1, invalid_v1=inv1, verdict_v2=v2, invalid_v2=inv2, info_v2=info2, eos_branches=br)
    stored_v1 = j.get("verdict") if ver == 1 else j.get("verdict_v1")
    if stored_v1 is not None and v1 != stored_v1:
        out["verdict_v2"] = "再判定不能"
        out["reason"] = f"v1 の再計算 {v1} が保存された v1 の判定 {stored_v1} と一致しない (事実の読み取りを疑う)"
    t = j.get("totals", {})
    out["eos_args_post_diff_bytes"] = t.get("eos_args_post_diff_bytes")
    out["eos_args_pre_diff_bytes"] = t.get("eos_args_pre_diff_bytes")
    out["eos_args_bytes"] = t.get("eos_args_bytes")
    ds = rep["datasets"]
    out["reads_pre_nonfinite"] = {nm: ds[nm]["nonfinite"]["A:pre"] + ds[nm]["nonfinite"]["B:pre"] for nm in reads if nm in ds}
    out["writes_post_nonfinite"] = {nm: ds[nm]["nonfinite"]["A:post"] + ds[nm]["nonfinite"]["B:post"] for nm in writes if nm in ds}
    return out


def print_report(rep):
    print(f"A: {rep['A']}\nB: {rep['B']}")
    m = rep.get("meta")
    if m:
        print(f"step {m.get('step')}  nCells {m.get('nCells')}  nCells_all {m.get('nCells_all')}  thermalMethod {m.get('thermalMethod')}"
              f"  nSpeciesRegistered {m.get('nSpeciesRegistered')}  thermoFloat {m.get('thermoFloat')}  discretization {m.get('discretization')}"
              + (f"  replay {m.get('replay_file')}" if m.get("replay_file") else ""))
    ds = {k: v for k, v in rep.get("datasets", {}).items() if not k.startswith("db/")}
    if ds:
        n_arg = sum(1 for v in ds.values() if v["class"] == "eos_arg")
        print(f"配列 {len(ds)} 本 (判定する EOS の配列 {n_arg} 本・その他 {len(ds) - n_arg} 本) × {{pre, post}}、物性 DB "
              f"{[k for k in rep['datasets'] if k.startswith('db/')]}")
        print(f"v2 の分岐 {rep.get('eos_branches')}  読む {rep.get('eos_reads_v2')}")
        print(f"  書く {rep.get('eos_writes_v2')}")
        print("EOS の配列 (pre 差分バイト / post 差分バイト / 非有限 A:pre,A:post,B:pre,B:post / v2 の分類):")
        reads, writes = set(rep.get("eos_reads_v2", [])), set(rep.get("eos_writes_v2", []))
        for nm in judged_arrays(rep.get("meta") or {}):
            d = ds.get(nm)
            if d:
                nf = d["nonfinite"]
                cls = ("読む・書く" if nm in reads and nm in writes else "読む" if nm in reads else "書く" if nm in writes else "不使用")
                print(f"  {nm:10s} {d['pre_diff_bytes']:>10d} {d['post_diff_bytes']:>10d}   "
                      f"{nf['A:pre']},{nf['A:post']},{nf['B:pre']},{nf['B:post']}   {cls}")
        t = rep["totals"]
        print(f"合計: EOS の配列の post 差分バイト {t['eos_args_post_diff_bytes']} / {t['eos_args_bytes']} バイト、"
              f"pre 差分バイト {t['eos_args_pre_diff_bytes']}; その他の配列 pre 差分 {t['other_pre_diff_bytes']}・post 差分 {t['other_post_diff_bytes']}"
              f" (判定外、EOS が読まない)")
        print(f"その他の配列で A/B の pre が異なるもの: {rep['other_pre_mismatch'] or 'なし'}")
        print(f"その他の配列で非有限を含むもの: {rep['other_nonfinite'] or 'なし'}")
        for tag in ("A", "B"):
            print(f"{tag}: EOS の前後で変わった配列 {rep['changed_by_eos'][tag]}")
        print(f"v1 の書き込み先の列挙: {rep['eos_writes']}")
    for d in rep["different"]:
        print(f"差: {d['array']} {d['n_elements']} 要素 / {d['diff_bytes']} バイト、最初の節点:")
        for n in d["first_nodes"]:
            print(f"    node {n['node']} ({'実節点' if n['real'] else 'ghost'}): A {n['A']}  B {n['B']}  入力 {n['input_pre']}")
    for s in rep["invalid_v1"]:
        print(f"不成立 (v1): {s}")
    for s in rep["invalid_v2"]:
        print(f"不成立 (v2): {s}")
    for s in rep.get("info_v2", []):
        print(f"情報 (v2): {s}")
    print(f"VERDICT_v1: {rep['verdict_v1']}  (旧規則: EOS の引数配列の pre・post 全部に有限性。記録用)")
    print(f"VERDICT_v2: {rep['verdict_v2']}  (改訂規則: 有限性は読む入力の pre と書き込み先の post。採用中・終了コード)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", nargs="?")
    ap.add_argument("b", nargs="?")
    ap.add_argument("--json", default=None, help="結果を JSON でも書く")
    ap.add_argument("--max-nodes", type=int, default=5)
    ap.add_argument("--rejudge-json", default=None, help="保存済みの比較 JSON を両規則で再判定する (ダンプは読まない)")
    a = ap.parse_args()
    if a.rejudge_json:
        out = rejudge_json(a.rejudge_json)
        for k in ("source", "stored_rule_version", "stored_verdict", "eos_branches", "eos_args_pre_diff_bytes",
                  "eos_args_post_diff_bytes", "eos_args_bytes", "reads_pre_nonfinite", "writes_post_nonfinite", "reason"):
            if k in out:
                print(f"{k}: {out[k]}")
        for s in out.get("invalid_v1", []):
            print(f"不成立 (v1): {s}")
        for s in out.get("invalid_v2", []):
            print(f"不成立 (v2): {s}")
        for s in out.get("info_v2", []):
            print(f"情報 (v2): {s}")
        print(f"VERDICT_v1: {out['verdict_v1']}")
        print(f"VERDICT_v2: {out['verdict_v2']}")
        if a.json:
            with open(a.json, "w") as fh:
                json.dump(out, fh, indent=1, ensure_ascii=False, default=str)
        return {"IDENTICAL": 0, "DIFFERENT": 1}.get(out["verdict_v2"], 2)
    if not (a.a and a.b):
        ap.error("A.h5 と B.h5 (または --rejudge-json) が要る")
    rep = compare(a.a, a.b, a.max_nodes)
    print_report(rep)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1, ensure_ascii=False, default=str)
    return {"IDENTICAL": 0, "DIFFERENT": 1, "INVALID": 2}[rep["verdict_v2"]]


if __name__ == "__main__":
    sys.exit(main())
