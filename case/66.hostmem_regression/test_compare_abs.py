#!/usr/bin/env python3
"""追加診断 B (compare_runs.judge_abs) の単体試験 (plan architecture-solver-host-memory §6.2)。

    python3 test_compare_abs.py        # 全試験を回し、最後に "ALL PASS" を出す (pytest でも回る)

codex (diagnose 2026-10-07) の最小再現 2 つで、登録判定 A の尺度 m = max|A−B|/max|A| の欠陥 (並び順で判定が変わる・
比較不能を PASS にする) を確かめ、B の尺度 d = max|A−B| では起きないことを確かめる。
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_runs as cr  # noqa: E402


def runs(*vals):
    """1 run = 1 要素の float32 配列 (forge の出力と同じ型)。"""
    return [np.array([v], dtype=np.float32) for v in vals]


def judge_A(base, new):
    """登録判定 A の (b) と同じ組み方 (compare_runs.compare の eval_set を写したもの)。"""
    pb = [(base[i], base[j]) for i in range(len(base)) for j in range(i + 1, len(base))]
    pn = [(new[i], new[j]) for i in range(len(new)) for j in range(i + 1, len(new))]
    px = [(b, n) for b in base for n in new]
    sb = max(cr.metric(a, b) for a, b in pb)
    sn = max(cr.metric(a, b) for a, b in pn)
    dd = max(cr.metric(a, b) for a, b in px)
    s = max(sb, sn)
    ok = (dd == 0.0) if s == 0.0 else (dd <= 2.0 * s)
    return s, dd, ("PASS" if ok else "FAIL")


def test_codex_repro_order_dependence():
    # base [100,2,1] / new [100,2,1] と new を [1,2,100] に並べ替えたもの
    b = runs(100, 2, 1)
    n1 = runs(100, 2, 1)
    n2 = runs(1, 2, 100)
    # A の欠陥の再現: 同じ値の集合で FAIL → 並べ替えで PASS
    assert judge_A(b, n1)[2] == "FAIL"
    assert judge_A(b, n2)[2] == "PASS"
    # B: 両方 PASS で同じ値 (S_abs = 99 = |100−1|、D_abs = 99)
    r1 = cr.judge_abs(b, n1)
    r2 = cr.judge_abs(b, n2)
    assert r1["verdict"] == r2["verdict"] == "PASS", (r1, r2)
    assert r1["S"] == r2["S"] == 99.0 and r1["D"] == r2["D"] == 99.0, (r1, r2)


def test_codex_repro_inf_pass():
    # base [0,1,1] / new [100,100,100]: A は inf ≤ 2·inf で PASS、B は S_abs 1・D_abs 100 で FAIL
    b = runs(0, 1, 1)
    n = runs(100, 100, 100)
    assert judge_A(b, n)[2] == "PASS"
    r = cr.judge_abs(b, n)
    assert r["verdict"] == "FAIL" and r["S"] == 1.0 and r["D"] == 100.0, r


def test_incomparable_is_fail():
    b = runs(1, 2, 3)
    assert cr.judge_abs(b, runs(1, 2, np.nan))["verdict"] == "FAIL"        # 非有限
    assert cr.judge_abs(b, runs(1, 2, np.inf))["verdict"] == "FAIL"
    assert cr.judge_abs(b, [b[0], b[1], None])["verdict"] == "FAIL"        # 欠落
    bad_shape = [np.zeros(2, np.float32), np.zeros(2, np.float32), np.zeros(3, np.float32)]
    r = cr.judge_abs(b, bad_shape)
    assert r["verdict"] == "FAIL" and "shape" in r["reason"], r
    assert cr.judge_abs(b[:1], runs(1, 2, 3))["verdict"] == "FAIL"         # 反復が足りない


def test_zero_spread():
    z = runs(5, 5, 5)
    assert cr.judge_abs(z, runs(5, 5, 5))["verdict"] == "PASS"            # S = 0、D = 0
    r = cr.judge_abs(z, runs(5, 5, 5.5))                                   # new 内で割れると S > 0 になる
    assert r["S"] == 0.5 and r["D"] == 0.5 and r["verdict"] == "PASS", r
    r = cr.judge_abs(runs(5, 5), runs(6, 6, 6))                            # S = 0 で D > 0 は FAIL
    assert r["S"] == 0.0 and r["verdict"] == "FAIL", r


def test_bound_is_inclusive():
    # D_abs = 2·S_abs ちょうどは PASS、超えたら FAIL (base 内の差 1、new 内の差 0 → S_abs = 1)
    r = cr.judge_abs(runs(0, 1, 1), runs(2, 2, 2))            # D_abs = |0−2| = 2 = 2·S_abs
    assert r["S"] == 1.0 and r["D"] == 2.0 and r["verdict"] == "PASS", r
    r = cr.judge_abs(runs(0, 1, 1), runs(2.5, 2.5, 2.5))      # D_abs = 2.5 > 2
    assert r["D"] == 2.5 and r["verdict"] == "FAIL", r


def test_integer_exact():
    a = [np.arange(5, dtype=np.int32) for _ in range(3)]
    assert cr.judge_abs(a, [x.copy() for x in a])["verdict"] == "PASS"
    c = [x.copy() for x in a]
    c[2][0] = 9
    assert cr.judge_abs(a, c)["verdict"] == "FAIL"
    s = [np.array(["inner_iter"], dtype=object) for _ in range(3)]
    assert cr.judge_abs(s, [x.copy() for x in s])["verdict"] == "PASS"


def test_self_check_random():
    rng = np.random.default_rng(1)
    for _ in range(20):
        b = [rng.normal(size=50).astype(np.float32) for _ in range(3)]
        n = [rng.normal(size=50).astype(np.float32) * rng.uniform(0.5, 3) for _ in range(3)]
        ref = cr.judge_abs(b, n)
        n_perm, n_bad, swap_ok = cr.self_check_abs(b, n, ref)
        assert n_perm == 36 and n_bad == 0 and swap_ok, (n_perm, n_bad, swap_ok)
        # base/new を入れ替えても D_abs は同じ (S は同じ集合なので同じ)
        sw = cr.judge_abs(n, b)
        assert sw["D"] == ref["D"] and sw["S"] == ref["S"]
        # 並べ替えの全通りで同じ値
        vals = {(cr.judge_abs(list(pb), list(pn))["S"], cr.judge_abs(list(pb), list(pn))["D"])
                for pb in itertools.permutations(b) for pn in itertools.permutations(n)}
        assert len(vals) == 1


def test_dabs_symmetric():
    rng = np.random.default_rng(2)
    a = rng.normal(size=1000).astype(np.float32)
    b = rng.normal(size=1000).astype(np.float32)
    assert cr.dabs(a, b) == cr.dabs(b, a)
    assert cr.dabs(a, a) == 0.0
    assert cr.dabs(a, b[:10]) is None


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    print("ALL PASS" if fails == 0 else f"{fails} FAIL")
    sys.exit(1 if fails else 0)
