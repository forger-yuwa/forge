#!/usr/bin/env python3
"""compare_eos_dump.py の単体試験 (合成ダンプ、GPU 不要)。

一致 → IDENTICAL、出力 1 ビット差 → DIFFERENT (配列・節点を特定)、入力不一致・非有限・構造不一致・設定不一致・
物性 DB 不一致・post 欠落・書き込み先の取りこぼし・凝縮/非定常の経路 → INVALID、EOS が読まない配列の差は判定外 (IDENTICAL のまま報告)。
"verdict" は改訂規則 (v2)。改訂規則の試験 (後半): 書き込み専用配列の pre の NaN は v2 で INVALID にしない (v1 は INVALID のまま)、
読む入力の pre の NaN と書き込み先の post の NaN は v2 でも INVALID、読み書きの集合の定数、保存済み JSON の再判定。

    python3 test_compare_eos_dump.py
"""
import json
import os
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_eos_dump as ced   # noqa: E402

N, NREAL = 50, 40
ok_all = True


def chk(name, cond, detail=""):
    global ok_all
    ok_all = ok_all and bool(cond)
    print("  [%s] %-60s %s" % ("OK " if cond else "NG ", name, detail))


def base_arrays(rng):
    pre = {}
    for nm in ced.EOS_ARGS_BASE + ["roY0", "roY1", "dUxdx", "volume"]:
        pre[nm] = rng.random(N, dtype=np.float32) + 0.5
    post = {k: v.copy() for k, v in pre.items()}
    for nm in ced.EOS_WRITES_TP:
        post[nm] = (pre[nm] * np.float32(1.25)).astype(np.float32)
    return pre, post


def write(path, pre, post, meta=None, db=b"\x01\x02\x03\x04", skip_post=False):
    m = {"step": 1, "nCells": NREAL, "nCells_all": N, "sizeof_flow_float": 4, "thermalMethod": 2, "nSpeciesRegistered": 2,
         "nSpecies": 2, "thermoFloat": 1, "condensation": 0, "isImplicit": 1, "unsteady": 0, "dualTime": 0,
         "discretization": "node", "solverConfig_yaml": "a: 1\n"}
    m.update(meta or {})
    with h5py.File(path, "w") as f:
        for k, v in m.items():
            f.attrs[k] = v
        f.create_dataset("db/species_thermo", data=np.frombuffer(db, dtype=np.uint8))
        for g, d in (("pre", pre), ("post", post)):
            if g == "post" and skip_post:
                continue
            for k, v in d.items():
                f.create_dataset(f"{g}/{k}", data=v)


def run_case(tmp, name, mutate_b=None, mutate_a=None, **kw_b):
    rng = np.random.default_rng(0)
    pre, post = base_arrays(rng)
    a = os.path.join(tmp, f"{name}_A.h5")
    b = os.path.join(tmp, f"{name}_B.h5")
    preA = {k: v.copy() for k, v in pre.items()}
    postA = {k: v.copy() for k, v in post.items()}
    preB = {k: v.copy() for k, v in pre.items()}
    postB = {k: v.copy() for k, v in post.items()}
    if mutate_a:
        mutate_a(preA, postA)
    if mutate_b:
        mutate_b(preB, postB)
    write(a, preA, postA)
    write(b, preB, postB, **kw_b)
    return ced.compare(a, b)


def flip(arr, i):
    v = arr.view(np.uint32)
    v[i] ^= np.uint32(1)


def pair(tmp, name, mut_a=None, mut_b=None, meta=None, cpg=False):
    """A・B に同じ改変 (mut_a) と B だけの改変 (mut_b) を入れたダンプを作って比べる。cpg=True は CPG 単相の形 (gamma・cp を書かない)。"""
    pre, post = base_arrays(np.random.default_rng(2))
    m = dict(meta or {})
    if cpg:
        for nm in ("gamma", "cp"):
            post[nm] = pre[nm].copy()
        m.update({"thermalMethod": 0, "nSpeciesRegistered": 1})
    preA, postA = {k: v.copy() for k, v in pre.items()}, {k: v.copy() for k, v in post.items()}
    preB, postB = {k: v.copy() for k, v in pre.items()}, {k: v.copy() for k, v in post.items()}
    if mut_a:
        mut_a(preA, postA); mut_a(preB, postB)
    if mut_b:
        mut_b(preB, postB)
    a, b = os.path.join(tmp, f"{name}_A.h5"), os.path.join(tmp, f"{name}_B.h5")
    write(a, preA, postA, meta=m)
    write(b, preB, postB, meta=m)
    return ced.compare(a, b), a, b


def v2_rule_tests():
    nan = np.float32("nan")
    print("--- 改訂規則 (v2) ---")
    # 読み書きの集合の定数
    br, rd, wr, un = ced.eos_sets({"thermalMethod": 2, "nSpeciesRegistered": 2, "condensation": 0})
    chk("集合: TP 単相の書く集合 = v1 の TP の書き込み先", sorted(wr) == sorted(ced.EOS_WRITES_TP) and br == ["tp_1ph"] and not un, str(wr))
    chk("集合: TP 単相の読む集合 = ρ・ρU・ρE・ρk・ρω・T・ρY_s", sorted(rd) == sorted(
        ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "T", "roY0", "roY1"]), str(rd))
    chk("集合: TP 単相は gamma・cp・Rmix・P・Ht・sonic・k・omega・速度を読まない",
        not (set(rd) & {"gamma", "cp", "Rmix", "P", "Ht", "sonic", "k", "omega", "Ux", "Uy", "Uz"}))
    br, rd, wr, un = ced.eos_sets({"thermalMethod": 0, "nSpeciesRegistered": 1, "condensation": 0})
    chk("集合: CPG 単相の書く集合 = v1 の CPG の書き込み先、読む集合 = 共通 (T を読まない)",
        sorted(wr) == sorted(ced.EOS_WRITES_CPG) and sorted(rd) == sorted(ced.EOS_READS_COMMON) and not un)
    br, rd, wr, un = ced.eos_sets({"thermalMethod": 2, "nSpeciesRegistered": 2, "condensation": 1, "nCondSpeciesRegistered": 1})
    chk("集合: TP 凝縮は単相と二相の和 (rog_0 を読む) で未対応", br == ["tp_1ph", "tp_2ph"] and "rog_0" in rd and un, str(un)[:60])
    br, rd, wr, un = ced.eos_sets({"thermalMethod": 0, "condensation": 1, "condEquilibrium": 2, "nCondSpeciesRegistered": 1})
    chk("集合: CPG 二相は T と rog_0 を読み、平衡 2 は rog_0 を書く", "T" in rd and "rog_0" in rd and "rog_0" in wr and un)
    with tempfile.TemporaryDirectory() as tmp:
        def wo_nan(pre, post):   # 書き込み専用 (TP の gamma・cp) の pre の ghost に NaN、post は有限
            pre["gamma"][NREAL:] = nan; pre["cp"][NREAL:] = nan
        r, a, b = pair(tmp, "v2_wonan", mut_a=wo_nan)
        chk("TP: 書き込み専用 gamma・cp の pre の NaN → v1 INVALID / v2 IDENTICAL",
            r["verdict_v1"] == "INVALID" and r["verdict_v2"] == "IDENTICAL" and r["verdict"] == "IDENTICAL",
            f"{r['verdict_v1']}/{r['verdict_v2']}")
        chk("TP: v2 は書き込み専用の pre の非有限を情報として出す", any("gamma" in s for s in r["info_v2"])
            and any("cp" in s for s in r["info_v2"]), "; ".join(r["info_v2"])[:100])
        for nm in ("roe", "T", "roY0", "roK"):
            def rd_nan(pre, post, nm=nm):   # 読む入力の pre に NaN (post は有限のまま: 比較器が pre を見ていることだけを試す)
                pre[nm][3] = nan
            r, _, _ = pair(tmp, f"v2_rdnan_{nm}", mut_a=rd_nan)
            chk(f"TP: 読む入力 {nm} の pre の NaN → v2 INVALID", r["verdict_v2"] == "INVALID"
                and any(nm in s and "読む入力" in s and "非有限" in s for s in r["invalid_v2"]), r["verdict_v2"])
        def rd_nan_ghost(pre, post):
            pre["roe"][N - 1] = nan
        r, _, _ = pair(tmp, "v2_rdnan_ghost", mut_a=rd_nan_ghost)
        chk("TP: 読む入力 roe の ghost の pre の NaN も v2 INVALID (ghost を除外しない)", r["verdict_v2"] == "INVALID")
        for nm in ("Ht", "gamma", "k"):
            def wr_nan(pre, post, nm=nm):
                post[nm][N - 2] = nan
            r, _, _ = pair(tmp, f"v2_wrnan_{nm}", mut_a=wr_nan)
            chk(f"TP: 書き込み先 {nm} の post の NaN (ghost) → v2 INVALID", r["verdict_v2"] == "INVALID"
                and any(nm in s and "post に非有限" in s for s in r["invalid_v2"]), r["verdict_v2"])
        def wo_mis(pre, post):
            flip(pre["gamma"], 1)
        r, _, _ = pair(tmp, "v2_womis", mut_b=wo_mis)
        chk("TP: 書き込み専用 gamma の pre が A/B で違う → v2 も INVALID (凍結入力の同一性)", r["verdict_v2"] == "INVALID"
            and any("gamma" in s and "pre が不一致" in s for s in r["invalid_v2"]))
        r, _, _ = pair(tmp, "v2_wonan_diff", mut_a=lambda p, q: (p["gamma"].__setitem__(slice(NREAL, None), nan)))
        r2, _, _ = pair(tmp, "v2_wonan_diff2", mut_a=lambda p, q: (p["gamma"].__setitem__(slice(NREAL, None), nan)),
                        mut_b=lambda p, q: flip(q["Ht"], 2))
        chk("TP: 書き込み専用の pre の NaN があっても出力差は v2 DIFFERENT", r["verdict_v2"] == "IDENTICAL" and r2["verdict_v2"] == "DIFFERENT",
            f"{r['verdict_v2']}/{r2['verdict_v2']}")
        r, _, _ = pair(tmp, "v2_unch", mut_b=lambda p, q: flip(q["roUy"], 4))
        chk("TP: 書かないはずの roUy が EOS の前後で変わる → v2 INVALID (未変更入力の不変)", r["verdict_v2"] == "INVALID"
            and any("roUy" in s and "不変" in s for s in r["invalid_v2"]))
        # CPG 単相: T は書き込み専用、gamma・cp は読みも書きもしない
        r, _, _ = pair(tmp, "v2_cpg_T", mut_a=lambda p, q: p["T"].__setitem__(slice(NREAL, None), nan), cpg=True)
        chk("CPG: T (書き込み専用) の pre の NaN → v1 INVALID / v2 IDENTICAL", r["verdict_v1"] == "INVALID" and r["verdict_v2"] == "IDENTICAL",
            f"{r['verdict_v1']}/{r['verdict_v2']}")
        def cpg_gc(pre, post):
            pre["cp"][NREAL:] = nan; post["cp"][NREAL:] = nan
        r, _, _ = pair(tmp, "v2_cpg_cp", mut_a=cpg_gc, cpg=True)
        chk("CPG: 読みも書きもしない cp の pre・post の NaN → v2 IDENTICAL・情報に出る", r["verdict_v2"] == "IDENTICAL"
            and any("cp" in s for s in r["info_v2"]), f"{r['verdict_v2']} {r['info_v2']}"[:100])
        r, _, _ = pair(tmp, "v2_cpg_roe", mut_a=lambda p, q: p["roe"].__setitem__(5, nan), cpg=True)
        chk("CPG: 読む入力 roe の pre の NaN → v2 INVALID", r["verdict_v2"] == "INVALID")
        r, _, _ = pair(tmp, "v2_cpg_Tpost", mut_a=lambda p, q: q["T"].__setitem__(5, nan), cpg=True)
        chk("CPG: 書き込み先 T の post の NaN → v2 INVALID", r["verdict_v2"] == "INVALID")
        r, _, _ = pair(tmp, "v2_cond", meta={"condensation": 1, "nCondSpeciesRegistered": 0})
        chk("凝縮 → v2 も INVALID (物性表がダンプに無い)", r["verdict_v2"] == "INVALID" and any("物性表" in s for s in r["invalid_v2"]))
        # 保存済み JSON の再判定 (旧形式 = rule_version 無し・invalid は v1 の文言)
        r, a, b = pair(tmp, "v2_rej", mut_a=wo_nan)
        legacy = {k: r[k] for k in ("A", "B", "different", "notes", "datasets", "meta", "eos_args", "eos_writes", "totals",
                                    "changed_by_eos", "other_pre_mismatch", "other_nonfinite")}
        legacy["invalid"], legacy["verdict"] = r["invalid_v1"], r["verdict_v1"]
        jp = os.path.join(tmp, "legacy.json")
        json.dump(legacy, open(jp, "w"), ensure_ascii=False)
        o = ced.rejudge_json(jp)
        chk("再判定: 旧形式 JSON (書き込み専用の NaN) → v1 INVALID (保存と一致)・v2 IDENTICAL",
            o["verdict_v1"] == "INVALID" and o["verdict_v2"] == "IDENTICAL", f"{o['verdict_v1']}/{o['verdict_v2']} {o.get('reason', '')}")
        r, a, b = pair(tmp, "v2_rej2", mut_a=lambda p, q: q["P"].__setitem__(2, nan))
        legacy2 = dict(legacy, datasets=r["datasets"], invalid=r["invalid_v1"], verdict=r["verdict_v1"])
        json.dump(legacy2, open(jp, "w"), ensure_ascii=False)
        o = ced.rejudge_json(jp)
        chk("再判定: 旧形式 JSON (書き込み先の post の NaN) → v2 INVALID", o["verdict_v2"] == "INVALID", o["verdict_v2"])
        bad = {k: v for k, v in legacy.items() if k != "datasets"}
        json.dump(bad, open(jp, "w"), ensure_ascii=False)
        o = ced.rejudge_json(jp)
        chk("再判定: 事実 (datasets) の無い JSON → 再判定不能", o["verdict_v2"] == "再判定不能", o.get("reason", ""))
        lie = dict(legacy, verdict="IDENTICAL")
        json.dump(lie, open(jp, "w"), ensure_ascii=False)
        o = ced.rejudge_json(jp)
        chk("再判定: v1 の再計算が保存と食い違う → 再判定不能", o["verdict_v2"] == "再判定不能", o.get("reason", ""))
        struct = dict(legacy, invalid=["設定・経路の属性 cp が不一致 (A 1 / B 2)"] + r["invalid_v1"])
        json.dump(struct, open(jp, "w"), ensure_ascii=False)
        o = ced.rejudge_json(jp)
        chk("再判定: 事実から再計算できない不成立 (設定の不一致) は v2 に引き継ぐ", o["verdict_v2"] == "INVALID"
            and any("設定" in s for s in o["invalid_v2"]))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        r = run_case(tmp, "same")
        chk("同一 → IDENTICAL", r["verdict"] == "IDENTICAL", r["verdict"])
        chk("同一: 変わった配列 = 書き込み先", sorted(r["changed_by_eos"]["A"]) == sorted(ced.EOS_WRITES_TP))

        r = run_case(tmp, "post1bit", mutate_b=lambda pre, post: flip(post["Ht"], 7))
        chk("出力 Ht の 1 ビット差 → DIFFERENT", r["verdict"] == "DIFFERENT", r["verdict"])
        d = r["different"][0] if r["different"] else {}
        chk("DIFFERENT: 最初の配列 Ht・節点 7・差分 1 バイト", d.get("array") == "Ht" and d["first_nodes"][0]["node"] == 7
            and d.get("diff_bytes") == 1, str(d.get("array")))

        def two(pre, post):
            flip(post["gamma"], 3); flip(post["T"], 45)
        r = run_case(tmp, "post2", mutate_b=two)
        chk("T と gamma の差 → 最初の配列は計算順で T、節点 45 は ghost", r["verdict"] == "DIFFERENT" and r["different"][0]["array"] == "T"
            and r["different"][0]["first_nodes"][0]["real"] is False, [x["array"] for x in r["different"]])

        r = run_case(tmp, "unchanged_input_post", mutate_b=lambda pre, post: flip(post["roUx"], 2))
        chk("書かないはずの入力 roUx の post 差 → INVALID (取りこぼし) か DIFFERENT", r["verdict"] in ("INVALID", "DIFFERENT")
            and any("roUx" in s for s in r["invalid"]), "; ".join(r["invalid"])[:120])

        def inmis(pre, post):
            flip(pre["roe"], 5); flip(post["roe"], 5)
        r = run_case(tmp, "input", mutate_b=inmis)
        chk("EOS 入力 roe (pre) の不一致 → INVALID", r["verdict"] == "INVALID" and any("roe" in s and "入力" in s for s in r["invalid"]))

        r = run_case(tmp, "roY", mutate_b=lambda pre, post: (flip(pre["roY1"], 0), flip(post["roY1"], 0)))
        chk("化学種 roY1 (pre) の不一致 → INVALID", r["verdict"] == "INVALID" and any("roY1" in s for s in r["invalid"]))

        def nan(pre, post):
            post["P"][4] = np.float32("nan")
        r = run_case(tmp, "nan", mutate_a=nan, mutate_b=nan)
        chk("EOS の出力に非有限 (両方同じ NaN) → INVALID", r["verdict"] == "INVALID" and any("非有限" in s for s in r["invalid"]))

        def other(pre, post):
            flip(pre["dUxdx"], 1); flip(post["dUxdx"], 1)
        r = run_case(tmp, "other", mutate_b=other)
        chk("EOS が読まない配列 (勾配) の差 → 判定外 (IDENTICAL)・報告に出る", r["verdict"] == "IDENTICAL"
            and r["other_pre_mismatch"] == ["dUxdx"], str(r["other_pre_mismatch"]))

        def other_written(pre, post):
            flip(post["volume"], 1)
        r = run_case(tmp, "gap", mutate_b=other_written)
        chk("EOS の前後で列挙外の配列が変わった → INVALID (取りこぼし)", r["verdict"] == "INVALID"
            and any("取りこぼし" in s and "volume" in s for s in r["invalid"]))

        r = run_case(tmp, "struct", mutate_b=lambda pre, post: (pre.pop("cp"), post.pop("cp")))
        chk("配列の欠落 → INVALID", r["verdict"] == "INVALID")

        r = run_case(tmp, "shape", mutate_b=lambda pre, post: (pre.__setitem__("Ux", pre["Ux"][:-1]), post.__setitem__("Ux", post["Ux"][:-1])))
        chk("shape の不一致 → INVALID", r["verdict"] == "INVALID")

        r = run_case(tmp, "dtype", mutate_b=lambda pre, post: (pre.__setitem__("k", pre["k"].astype(np.float64)),
                                                              post.__setitem__("k", post["k"].astype(np.float64))))
        chk("dtype の不一致 → INVALID", r["verdict"] == "INVALID")

        r = run_case(tmp, "cfg", meta={"solverConfig_yaml": "a: 2\n"})
        chk("設定 (solverConfig.yaml) の不一致 → INVALID", r["verdict"] == "INVALID")

        r = run_case(tmp, "db", db=b"\x01\x02\x03\x05")
        chk("物性 DB の不一致 → INVALID", r["verdict"] == "INVALID" and any("物性 DB" in s for s in r["invalid"]))

        r = run_case(tmp, "nopost", skip_post=True)
        chk("post 欠落 (EOS の後まで書けていない) → INVALID", r["verdict"] == "INVALID")

        rng = np.random.default_rng(0)
        pre, post = base_arrays(rng)
        a, b = os.path.join(tmp, "dual_A.h5"), os.path.join(tmp, "dual_B.h5")
        write(a, pre, post, meta={"unsteady": 1, "dualTime": 1})
        write(b, pre, post, meta={"unsteady": 1, "dualTime": 1})
        r = ced.compare(a, b)
        chk("非定常 (dual-time) の経路 → INVALID", r["verdict"] == "INVALID")
        write(a, pre, post, meta={"condensation": 1})
        write(b, pre, post, meta={"condensation": 1})
        r = ced.compare(a, b)
        chk("凝縮 → INVALID (未対応)", r["verdict"] == "INVALID")
        r = ced.compare(a, os.path.join(tmp, "missing.h5"))
        chk("ファイルが無い → INVALID", r["verdict"] == "INVALID")

        # CPG 単相: gamma/cp は書かない
        pre, post = base_arrays(np.random.default_rng(1))
        for nm in ("gamma", "cp"):
            post[nm] = pre[nm].copy()
        write(a, pre, post, meta={"thermalMethod": 0, "nSpeciesRegistered": 1})
        write(b, pre, post, meta={"thermalMethod": 0, "nSpeciesRegistered": 1})
        r = ced.compare(a, b)
        chk("CPG 単相: gamma/cp 不変で IDENTICAL", r["verdict"] == "IDENTICAL", r["verdict"])
        post2 = {k: v.copy() for k, v in post.items()}
        flip(post2["cp"], 0)
        write(b, pre, post2, meta={"thermalMethod": 0, "nSpeciesRegistered": 1})
        r = ced.compare(a, b)
        chk("CPG 単相: cp が変わった → INVALID (CPG の書き込み先に無い)", r["verdict"] == "INVALID"
            and any("取りこぼし" in s for s in r["invalid"]))

    v2_rule_tests()
    print("\nVERDICT: %s" % ("PASS" if ok_all else "FAIL"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
