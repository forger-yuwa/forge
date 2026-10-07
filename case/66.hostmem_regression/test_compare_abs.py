#!/usr/bin/env python3
"""追加診断 B (compare_runs.judge_abs) の単体試験 (plan architecture-solver-host-memory §6.2)。

    python3 test_compare_abs.py        # 全試験を回し、最後に "ALL PASS" を出す (pytest でも回る)

codex (diagnose 2026-10-07) の最小再現 2 つで、登録判定 A の尺度 m = max|A−B|/max|A| の欠陥 (並び順で判定が変わる・
比較不能を PASS にする) を確かめ、B の尺度 d = max|A−B| では起きないことを確かめる。

後半 (test_complete_* / test_fw_*) は比較器の完全性検査 (plan §5.1 #6、2026-10-07 result レビュー M1) の試験。
一時ディレクトリに小さな run (h5・残差 CSV・.state・NANCHECK.txt・solverConfig.yaml) を合成し、存在しないディレクトリ・
反復不足・初期出力だけ・予定 N 未到達・必須列の欠落・比較量 0 が A・B・固定幅のどれでも合格にならないことを確かめる
(陽性対照として、そろった run は完全性 PASS になることも確かめる)。
"""
import contextlib
import csv
import io
import itertools
import json
import os
import subprocess
import sys
import tempfile
import traceback
from types import SimpleNamespace

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_runs as cr  # noqa: E402
import fixedwidth_eval as fw  # noqa: E402
import matrix_spec as ms  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


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


# ================================================================== 完全性検査 (plan §5.1 #6、result レビュー M1)
RES_COLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe"]


def make_run(root, name, n=4, out=None, rc=0, nan="PASS", state="done", res_steps=None, csv_n=None, csv_cols=None,
             cfg_n=None, sst=False, total_time=0.5, value=1.0, last_phase="outer_end", extra=()):
    """forge の run を 1 本合成する (run_matrix.py のワーカーが残すものと同じ形)。
    res_steps: 書く res_* の step (既定は 0 と予定の出力すべて)。csv_n: 残差 CSV の step 数 (既定 n)。"""
    d = os.path.join(root, name)
    os.makedirs(d)
    out = out or n
    with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
        f.write(f"time:\n  last: {{nStepOuter: {cfg_n or n}}}\n  outStepInterval: {out}\n"
                f"turbulence: {{model: {'sst' if sst else 'none'}}}\n")
    with open(os.path.join(d, "INPUT_FILES"), "w") as f:
        f.write("solverConfig.yaml\n")
    if state is not None:
        with open(os.path.join(d, ".state"), "w") as f:
            f.write(f"done rc={rc} nan={nan} wall=1s\n" if state == "done" else f"{state}\n")
    if nan is not None:
        with open(os.path.join(d, "NANCHECK.txt"), "w") as f:
            f.write(f"h5: 浮動小数データセット 2 個を検査\nNANCHECK: {nan}\n")
    steps = [0] + [k * out for k in range(1, n // out + 1)] if res_steps is None else res_steps
    for s in steps:
        with h5py.File(os.path.join(d, f"res_{s}.h5"), "w") as f:
            f.create_dataset("VALUE/ro", data=np.array([1.0, 2.0, 3.0], dtype=np.float32) * value)
            f.create_dataset("MESH/x", data=np.arange(3, dtype=np.float32))
            if s > 0 and total_time is not None:
                f.create_group("CHECKPOINT").attrs["totalTime"] = total_time
    cols = csv_cols if csv_cols is not None else RES_COLS + (["rms_roK", "rms_roOmega"] if sst else [])
    with open(os.path.join(d, "residual_history.csv"), "w") as f:
        w = csv.writer(f)
        w.writerow(["step", "inner", "phase"] + cols)
        for s in range(n if csv_n is None else csv_n):
            w.writerow([s, -1, last_phase] + [f"{1e-3 * value / (s + 1):.6e}"] * len(cols))
    for fn in extra:
        open(os.path.join(d, fn), "w").write("t , v\n0 , 1\n")
    return d


def make_convert(root, name, rc=1, written=1, datasets=True, mesh=True):
    """変換器の run を 1 本合成する。written=None は written の記録が無い .state (記録を足す前のワーカー)。"""
    d = os.path.join(root, name)
    os.makedirs(d)
    with open(os.path.join(d, ".state"), "w") as f:
        f.write(f"done rc={rc} nan=PASS" + ("" if written is None else f" written={written}") + " wall=1s\n")
    with open(os.path.join(d, "NANCHECK.txt"), "w") as f:
        f.write("NANCHECK: PASS\n")
    with h5py.File(os.path.join(d, "converted.h5"), "w") as f:
        f.create_group("VALUE")
        if mesh:
            f.create_group("MESH")
        if datasets:
            f.create_dataset("VALUE/ro", data=np.ones(3, dtype=np.float32))
    return d


def sub(t, name):
    d = os.path.join(t, name)
    os.makedirs(d)
    return d


def plan4(**kw):
    return ms.plan_from("t", "forge", 3, 4, None, kw.get("extra", ()))


def six(root, **kw):
    b = [make_run(root, f"b{i}", **kw) for i in range(1, 4)]
    n = [make_run(root, f"n{i}", **kw) for i in range(1, 4)]
    return b, n


def comp_fail(c, *needles):
    """完全性が FAIL で、理由のどれかに needles が全部入っていること。"""
    assert c["verdict"] == "FAIL", c
    txt = " | ".join(c["reasons"])
    for x in needles:
        assert x in txt, (x, txt)


def test_complete_positive_control():
    # 陽性対照: 予定どおりそろった 3 + 3 本は完全性 PASS、A・B とも量 > 0 で比較する
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t)
        c = cr.completeness(plan4(), b, n)
        assert c["verdict"] == "PASS", c["reasons"]
        lines, rows, sc = cr.compare_abs(b, n, "t", plan4())
        assert sc["quantities"] > 0 and all(r["verdict"] == "PASS" for r in rows), rows
        assert any(r["file"] == "res_4.h5" for r in rows) and not any(r["file"] == "res_0.h5" for r in rows)
        rep = cr.compare(b, n, cr.Report(), "t", plan4())
        assert rep.verdicts["NSTEP"] == "PASS" and rep.verdicts["INIT_OUT"] == "PASS", rep.verdicts


def test_complete_missing_directory():
    # codex の再現: 存在しないディレクトリを base/new 各 3 本で渡すと、以前は「量 0、FAIL 0」だった
    with tempfile.TemporaryDirectory() as t:
        b = [os.path.join(t, f"nope_b{i}") for i in range(3)]
        n = [os.path.join(t, f"nope_n{i}") for i in range(3)]
        comp_fail(cr.completeness(plan4(), b, n), "ディレクトリが無い")
        # 完全性検査を通さず compare_abs を直接呼んでも、比較量 0 は FAIL の行になる
        _, rows, sc = cr.compare_abs(b, n, "t", plan4())
        assert sc["quantities"] == 0 and any(r["verdict"] == "FAIL" and r["file"] == "(比較量)" for r in rows), rows
        # 1 本だけ無い場合も FAIL
        b2, n2 = six(sub(t, "x"))
        comp_fail(cr.completeness(plan4(), b2[:2] + [os.path.join(t, "gone")], n2), "gone: ディレクトリが無い")
        # CLI (B と A): 完全性 FAIL が summary・completeness.tsv・all_quantities.tsv に出る
        od = os.path.join(t, "outB")
        r = subprocess.run([sys.executable, os.path.join(HERE, "compare_runs.py"), "--metric", "abs", "--tag", "t",
                            "--base", *b, "--new", *n, "--steps", "4", "--out-dir", od], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        assert "完全性 FAIL" in open(os.path.join(od, "summary.txt")).read()
        rows = list(csv.DictReader(open(os.path.join(od, "completeness.tsv")), delimiter="\t"))
        assert rows[0]["verdict"] == "FAIL" and "ディレクトリが無い" in rows[0]["reasons"], rows
        aq = list(csv.DictReader(open(os.path.join(od, "all_quantities.tsv")), delimiter="\t"))
        assert len(aq) == 1 and aq[0]["verdict"] == "FAIL" and aq[0]["file"] == "(完全性)", aq
        od = os.path.join(t, "outA")
        r = subprocess.run([sys.executable, os.path.join(HERE, "compare_runs.py"), "--tag", "t",
                            "--base", *b, "--new", *n, "--steps", "4", "--out-dir", od], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        assert "COMPLETE FAIL" in open(os.path.join(od, "summary.txt")).read()


def test_complete_no_planned_n():
    # 直接指定で予定 N (--steps) が無いと、run がそろっていても FAIL (存在する最大 step を最終と見なさない)
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t)
        comp_fail(cr.completeness(ms.plan_from("t", "forge", 3, None), b, n), "予定 N が無い")


def test_complete_too_few_reps():
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t)
        comp_fail(cr.completeness(plan4(), b[:2], n), "base 2 本 (予定 3 本)")
        comp_fail(cr.completeness(plan4(), b, n[:1]), "new 1 本 (予定 3 本)")
        comp_fail(cr.completeness(plan4(), b, []), "new 0 本")
        comp_fail(cr.completeness(plan4(), [b[0], b[0], b[1]], n), "2 回以上")
        # 予定の反復数は構成表から (既定 3)
        assert ms.planned("c36node")["reps"] == ms.REPS == 3


def test_complete_registry_selection():
    # registry の rep 1..3 を 1 本ずつ選ぶ。足りない rep・除外されていない重複・終わっていない run は FAIL
    with tempfile.TemporaryDirectory() as t:
        reg = [("run_0001_x_base_r1", "base", 1), ("run_0002_x_base_r2", "base", 2), ("run_0003_x_base_r3", "base", 3),
               ("run_0004_x_new_r1", "new", 1), ("run_0005_x_new_r2", "new", 2),           # new r3 が無い
               ("run_0006_x_new_r3", "new", 3), ("run_0007_x_base_r4", "base", 4)]          # new r3 は .exclude、r4 は予定外
        for name, _b, _r in reg:
            make_run(t, name)
        open(os.path.join(t, "run_0006_x_new_r3", ".exclude"), "w").close()
        with open(os.path.join(t, "registry.tsv"), "w") as f:
            f.write("run\tcfg\tbuild\trep\n" + "".join(f"{a}\tx\t{b}\t{r}\n" for a, b, r in reg))
        p = plan4()
        b, wb, nb = cr.select_runs("x", "base", p, root=t)
        assert len(b) == 3 and not wb and any("予定外" in x for x in nb), (b, wb, nb)
        n, wn, _ = cr.select_runs("x", "new", p, root=t)
        assert len(n) == 2 and any("new r3" in x for x in wn), wn
        comp_fail(cr.completeness(p, b, n, wb + wn), "new r3 の run が registry に無い")
        # 同じ rep が除外されずに 2 本 → どれを使うか決まらない
        make_run(t, "run_0008_x_base_r1")
        with open(os.path.join(t, "registry.tsv"), "a") as f:
            f.write("run_0008_x_base_r1\tx\tbase\t1\n")
        _, wb2, _ = cr.select_runs("x", "base", p, root=t)
        assert any("base r1 が複数" in x for x in wb2), wb2
        # registry にあるがディレクトリが無い run
        with open(os.path.join(t, "registry.tsv"), "w") as f:
            f.write("run\tcfg\tbuild\trep\nrun_0099_y_base_r1\ty\tbase\t1\n")
        b3, w3, _ = cr.select_runs("y", "base", p, root=t)
        comp_fail(cr.completeness(p, b3, [], w3, need_new=False), "ディレクトリが無い", "base r2")


def test_complete_cli_registry():
    # registry から選ぶ経路 (--cfg、A と B) の通しの試験: そろえば完全性 PASS で量 > 0、new 1 本の最終出力を消すと構成ごと FAIL
    with tempfile.TemporaryDirectory() as t:
        reg = [(f"run_{i:04d}_c09ckpt100_{b}_r{r}", b, r) for i, (b, r) in enumerate(
            [(b, r) for b in ("base", "new") for r in (1, 2, 3)], start=1)]
        for name, _b, _r in reg:
            make_run(t, name, n=100)
        with open(os.path.join(t, "registry.tsv"), "w") as f:
            f.write("run\tcfg\tbuild\trep\n" + "".join(f"{a}\tc09ckpt100\t{b}\t{r}\n" for a, b, r in reg))

        def cli(metric, od):
            r = subprocess.run([sys.executable, os.path.join(HERE, "compare_runs.py"), "--metric", metric, "--cfg", "c09ckpt100",
                                "--new-build", "new", "--root", t, "--out-dir", os.path.join(t, od)], capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
            return (open(os.path.join(t, od, "summary.txt")).read(),
                    list(csv.DictReader(open(os.path.join(t, od, "completeness.tsv")), delimiter="\t")))
        sB, cB = cli("abs", "B1")
        assert cB[0]["verdict"] == "PASS" and "完全性 PASS" in sB and "FAIL   0" in sB, sB
        aq = list(csv.DictReader(open(os.path.join(t, "B1", "all_quantities.tsv")), delimiter="\t"))
        assert len(aq) > 0 and all(r["verdict"] == "PASS" for r in aq), aq
        sA, cA = cli("m", "A1")
        assert cA[0]["verdict"] == "PASS" and "COMPLETE PASS" in sA and "NSTEP PASS" in sA, sA
        os.remove(os.path.join(t, reg[4][0], "res_100.h5"))                 # new r2 の最終出力が無い
        sB, cB = cli("abs", "B2")
        assert cB[0]["verdict"] == "FAIL" and "res_100.h5" in cB[0]["reasons"] and "完全性 FAIL" in sB, sB
        sA, cA = cli("m", "A2")
        assert cA[0]["verdict"] == "FAIL" and "COMPLETE FAIL" in sA and "NSTEP" not in sA, sA


def test_complete_run_failed():
    # 実行成功でない run: 終了コード ≠ 0・終わっていない・.state が無い・NaN 検査 FAIL / 無し
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t)
        bad = [make_run(t, "rc1", rc=1), make_run(t, "running", state="running"), make_run(t, "nostate", state=None),
               make_run(t, "nanfail", nan="FAIL"), make_run(t, "nonan", nan=None)]
        for d, needle in zip(bad, ["終了コード 1", "終了していない", ".state が無い", "NaN 検査が合格でない", "NANCHECK.txt が無い"]):
            comp_fail(cr.completeness(plan4(), b[:2] + [d], n), needle)


def test_complete_initial_output_only():
    # 初期出力 res_0 だけで最終出力が無い: 旧い final_files は res_0 を「最終」と見なしていた
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t, res_steps=[0])
        assert cr.final_files(b[0])[0] == ["res_0.h5"]          # 旧い挙動 (存在する最大 step = 0 を最終に)
        assert cr.final_files(b[0], 4)[0] == []                 # 予定 N で選ぶと最終出力は無い
        comp_fail(cr.completeness(plan4(), b, n), "必要ファイルが無い: res_4.h5", "最終出力 res_4.h5 が無い")


def test_complete_n_not_reached():
    with tempfile.TemporaryDirectory() as t:
        # (i) 全反復が同じ所 (step 2) で止まった: 最後の出力 res_2・残差 CSV は step 0..1
        b, n = six(sub(t, "i"), res_steps=[0, 2], csv_n=2)
        assert cr.final_files(b[0])[0] == ["res_2.h5"]           # 旧い挙動は res_2 を最終と見なし、全反復一致で合格していた
        comp_fail(cr.completeness(plan4(), b, n), "最終出力 res_4.h5 が無い", "step が 0..3 でない")
        # (ii) 設定ごと違う (nStepOuter 2 で 2 step を完走): 構成表の予定 N 4 と違うので FAIL
        b, n = six(sub(t, "ii"), n=2, out=2)
        comp_fail(cr.completeness(plan4(), b, n), "nStepOuter 2 が予定 N 4 と違う", "res_4.h5 が無い")
        # (iii) 予定 N より後の出力がある (N を取り違えた run)
        b, n = six(sub(t, "iii"), res_steps=[0, 4, 8])
        comp_fail(cr.completeness(plan4(), b, n), "予定 N 4 より後の出力")
        # (iv) 最後の残差行が outer_end でない (最後の step の途中で止まった)
        b, n = six(sub(t, "iv"), last_phase="inner_iter")
        comp_fail(cr.completeness(plan4(), b, n), "outer_end でない")


def test_complete_missing_required_column():
    with tempfile.TemporaryDirectory() as t:
        b, n = six(sub(t, "a"), csv_cols=["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz"])   # rms_roe が無い
        comp_fail(cr.completeness(plan4(), b, n), "必須の列が無い: rms_roe")
        b, n = six(sub(t, "b"), sst=True, csv_cols=RES_COLS)                                 # SST なのに rms_roK 等が無い
        comp_fail(cr.completeness(plan4(), b, n), "rms_roK", "rms_roOmega")
        b, n = six(sub(t, "c"), sst=True)                                                     # 対照: SST で列あり
        assert cr.completeness(plan4(), b, n)["verdict"] == "PASS"
        # 構成ごとの追加の必要ファイル (境界出力・probe・診断 CSV) が無い
        b, n = six(sub(t, "d"))
        comp_fail(cr.completeness(plan4(extra=["res_wall_4_{N}.h5", "diag.csv"]), b, n), "res_wall_4_4.h5", "diag.csv")


def test_complete_zero_quantities():
    # 完全性は通るが比較できる量が 0 (空の converted.h5): A・B とも FAIL (量 0・FAIL 0 を合格にしない)
    with tempfile.TemporaryDirectory() as t:
        b = [make_convert(t, f"b{i}", datasets=False) for i in range(3)]
        n = [make_convert(t, f"n{i}", datasets=False) for i in range(3)]
        p = ms.plan_from("v", "convert", 3)
        assert cr.completeness(p, b, n)["verdict"] == "PASS"
        _, rows, sc = cr.compare_abs(b, n, "v", p)
        assert sc["quantities"] == 0 and [r["file"] for r in rows if r["verdict"] == "FAIL"] == ["(比較量)"], rows
        rep = cr.compare(b, n, cr.Report(), "v", p)
        assert rep.verdicts["NSTEP"].startswith("FAIL") and "量 0" in rep.verdicts["NSTEP"], rep.verdicts


def test_complete_converter_rc():
    # 変換器: 既知の終了時 rc 1 は written=1 なら可 (今の扱い)。written≠1・rc が 0/1 以外は FAIL
    with tempfile.TemporaryDirectory() as t:
        p = ms.plan_from("v", "convert", 3)
        b = [make_convert(t, f"b{i}", rc=1) for i in range(3)]
        n = [make_convert(t, f"n{i}", rc=0) for i in range(3)]
        assert cr.completeness(p, b, n)["verdict"] == "PASS"
        comp_fail(cr.completeness(p, b, n[:2] + [make_convert(t, "w0", rc=1, written=0)]), "変換結果が書かれていない")
        comp_fail(cr.completeness(p, b, n[:2] + [make_convert(t, "segv", rc=-11)]), "終了コード -11")
        # written の記録が無い .state は converted.h5 の /VALUE・/MESH で判定 (run_matrix.py の扱いのまま)
        assert cr.completeness(p, b, n[:2] + [make_convert(t, "old", written=None)])["verdict"] == "PASS"
        comp_fail(cr.completeness(p, b, n[:2] + [make_convert(t, "oldbad", written=None, mesh=False)]), "/VALUE か /MESH が無い")
        _, rows, sc = cr.compare_abs(b, n, "v", p)
        assert sc["quantities"] == 1 and all(r["verdict"] == "PASS" for r in rows), rows


def test_complete_physical_time_mismatch():
    # 最終出力の物理時刻 (CHECKPOINT の totalTime) が run で違う → FAIL
    with tempfile.TemporaryDirectory() as t:
        b, n = six(t)
        odd = make_run(t, "odd", total_time=0.75)
        comp_fail(cr.completeness(plan4(), b, n[:2] + [odd]), "物理時刻", "totalTime")


def test_complete_spec_planned():
    # 構成表の予定: 全構成で決まる。中間出力・境界出力・変換器の必要ファイル
    for c in ms.CONFIGS:
        p = ms.planned(c)
        assert p["reps"] == 3 and p["required"], (c, p)
        assert (p["N"] is None) == (p["kind"] == "convert"), (c, p)
    p = ms.planned("c09cont200")
    assert p["N"] == 200 and p["out"] == 100 and {"res_0.h5", "res_100.h5", "res_200.h5", "residual_history.csv"} <= set(p["required"])
    assert "res_wall_4_200.h5" in ms.planned("c36node")["required"] and "diag.csv" in ms.planned("c36node_impdiag")["required"]
    assert ms.planned("v44")["required"] == ["converted.h5"]


def _fw_setup(t, bad=None):
    """固定幅 eval の最小の計画 (構成 c44dual_ckpt100: N 100) を合成。bad = (run 名, make_run の引数) で 1 本を不完全にする。"""
    bnd = ["res_outlet_2_100.h5", "res_wall_3_100.h5"]       # 構成表 EXTRA_OUTPUTS の境界出力
    assert set(bnd) <= set(ms.planned("c44dual_ckpt100")["required"])
    names = ["base_e1", "base_e2", "base_e3"] + [f"fw_{g}{i}" for i in range(1, 7) for g in ("b", "n")]
    for nm in names:
        d = make_run(t, nm, n=100, **(bad[1] if bad and bad[0] == nm else {}))
        if os.path.exists(os.path.join(d, "res_100.h5")):      # 境界出力は最終出力と同じ時点にだけ書く
            for fn in bnd:
                with h5py.File(os.path.join(d, fn), "w") as f:
                    f.create_dataset("VALUE/q", data=np.ones(2, dtype=np.float32))
    pdir = os.path.join(t, "fwplan")
    os.makedirs(pdir)
    tp = os.path.join(pdir, "T_frozen.tsv")
    with open(tp, "w") as f:
        f.write("file\tquantity\tmode\tS0\tT\tS0_pair\nres_100.h5\tVALUE/ro\t幅\t1.0\t2.0\tB1-B2\n")
    plan = dict(input="c44dual_ckpt100", nStepOuter=100, outStepInterval=100, base_existing=names[:3],
                T_frozen_sha256=fw.sha256(tp),
                runs=[dict(run=nm, group=nm[3], build="x", rep=1) for nm in names[3:]])
    json.dump(plan, open(os.path.join(pdir, "plan.json"), "w"))
    return SimpleNamespace(plan=os.path.join(pdir, "plan.json"), out=os.path.join(pdir, "result"), root=t)


def _fw_verdict(a):
    with contextlib.redirect_stdout(io.StringIO()):
        fw.cmd_eval(a)
    txt = open(os.path.join(a.out, "fixedwidth_result.txt")).read()
    return txt.strip().splitlines()[-1], txt


def test_fw_positive_and_incomplete():
    # 固定幅 eval: そろった 15 本は全量幅内 (陽性対照)。新しい run が初期出力だけ・既存 base が N 未到達なら判定不能
    with tempfile.TemporaryDirectory() as t:
        v, _ = _fw_verdict(_fw_setup(t))
        assert "全量で幅内" in v, v
    with tempfile.TemporaryDirectory() as t:
        v, txt = _fw_verdict(_fw_setup(t, ("fw_n3", dict(res_steps=[0]))))
        assert "判定不能" in v and "不完全" in v and "最終出力 res_100.h5 が無い" in txt, v
    with tempfile.TemporaryDirectory() as t:
        v, txt = _fw_verdict(_fw_setup(t, ("base_e2", dict(res_steps=[0, 50], csv_n=50))))
        assert "判定不能" in v and "既存 base" in v, v
    with tempfile.TemporaryDirectory() as t:
        v, _ = _fw_verdict(_fw_setup(t, ("fw_b1", dict(rc=139))))
        assert "判定不能" in v, v


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
            except Exception:   # 試験自体の例外も合格にしない
                fails += 1
                print(f"FAIL {name}: 例外\n{traceback.format_exc()}")
    print("ALL PASS" if fails == 0 else f"{fails} FAIL")
    sys.exit(1 if fails else 0)
