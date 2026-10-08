#!/usr/bin/env python3
"""compare_eos_dump.py の単体試験 (合成ダンプ、GPU 不要)。

一致 → IDENTICAL、出力 1 ビット差 → DIFFERENT (配列・節点を特定)、入力不一致・非有限・構造不一致・設定不一致・
物性 DB 不一致・post 欠落・書き込み先の取りこぼし・凝縮/非定常の経路 → INVALID、EOS が読まない配列の差は判定外 (IDENTICAL のまま報告)。

    python3 test_compare_eos_dump.py
"""
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

    print("\nVERDICT: %s" % ("PASS" if ok_all else "FAIL"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
