#!/usr/bin/env python3
r"""V6′ 評価器と帯内 G-if の**負例**試験 (codex result 8 巡目 M1/M2 の再発防止)。

実 run (`run_0012_v6p_df5_i50_nlog`) の成果物を一時ディレクトリに写し (大きいファイルはリンク)、
1 か所だけ壊して、評価器・判定ツールが**合格させないこと**を確かめる。正例 (無改変) が PASS することも確かめる。

    python3 case/58.conjugate_slot/test_v6p_negative.py [--run <run_dir>]
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EVAL = os.path.join(HERE, "eval_v6p.py")
GIF = os.path.join(ROOT, "solver_density_cuda", "tools", "check_cht_interface.py")
BAND = ["-0.004166299", "-0.019779487"]


def clone(src, dst, copy=()):
    os.makedirs(dst)
    for f in os.listdir(src):
        s = os.path.join(os.path.abspath(src), f)
        if os.path.isdir(s):
            continue
        # **リンクは読むだけの h5 に限る**。評価器・判定ツールは run 直下に成果物 (v6p_band*.json/csv・
        # *_VERDICT.txt) を書くので、それらまでリンクすると元 run の成果物を上書きする (初版で実際に起きた)
        if f in copy or not f.endswith(".h5"):
            shutil.copy(s, os.path.join(dst, f))
        else:
            os.symlink(s, os.path.join(dst, f))


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = p.stdout + p.stderr
    v = [l for l in out.splitlines() if l.startswith("VERDICT")]
    return (v[-1] if v else "(no VERDICT)"), out


def edit_csv(path, fn):
    lines = open(path).read().splitlines()
    hdr, body = lines[0], [l.split(",") for l in lines[1:]]
    fn(body)
    open(path, "w").write("\n".join([hdr] + [",".join(r) for r in body]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.path.join(HERE, "run_0012_v6p_df5_i50_nlog"))
    a = ap.parse_args()
    step = 100000
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        def case(name, copy, mutate, tool, expect_pass):
            d = os.path.join(tmp, name)
            clone(a.run, d, copy)
            mutate(d)
            cmd = ([sys.executable, EVAL, d, "--fix-band", *BAND] if tool == "eval"
                   else [sys.executable, GIF, d, "--band-y", *BAND])
            v, _ = run(cmd)
            ok = ("PASS" in v) == expect_pass
            results.append(ok)
            print(f"[{'OK' if ok else 'NG'}] {name:<34} 期待 {'PASS' if expect_pass else '不合格'}  -> {v}")

        def flip_back(d):
            with h5py.File(os.path.join(d, f"res_slot_back_6_{step}.h5"), "r+") as f:
                f["VALUE/iface_q_eff"][...] = -f["VALUE/iface_q_eff"][:]

        def nan_coord(d):   # 帯上端の節点 212 の座標を NaN、その節点の残差を大きく
            edit_csv(os.path.join(d, "conjugate_iface_nodes_5.csv"), lambda b: b[212].__setitem__(3, "nan"))
            def f(b):
                for r in b[-321:]:
                    if r[2] == "212":
                        r[3] = "1000"
            edit_csv(os.path.join(d, "conjugate_iface_log_5.csv"), f)

        def stale_step(d):  # 最終更新の節点 212 だけ step=0
            def f(b):
                for r in b[-321:]:
                    if r[2] == "212":
                        r[1] = "0"
            edit_csv(os.path.join(d, "conjugate_iface_log_5.csv"), f)

        def dup_id(d):      # 最終更新で節点 212 の行を 211 に書き換え (重複 + 欠け)
            def f(b):
                for r in b[-321:]:
                    if r[2] == "212":
                        r[2] = "211"
            edit_csv(os.path.join(d, "conjugate_iface_log_5.csv"), f)

        wall_b = f"res_slot_back_6_{step}.h5"
        logs = ("conjugate_iface_log_5.csv", "conjugate_iface_nodes_5.csv", "v6p_band.json")
        case("正例 eval (無改変)", ("v6p_band.json",), lambda d: None, "eval", True)
        case("正例 G-if (無改変)", logs, lambda d: None, "gif", True)
        case("M1 後壁 q_eff の符号反転", (wall_b, "v6p_band.json"), flip_back, "eval", False)
        case("M2 節点 212 の座標 NaN + 残差 1000", logs, nan_coord, "gif", False)
        case("M2 最終更新で節点 212 だけ step=0", logs, stale_step, "gif", False)
        case("M2 最終更新で節点 ID 重複", logs, dup_id, "gif", False)
    n_ok = sum(results)
    print(f"\nVERDICT: {'PASS' if n_ok == len(results) else 'FAIL'} ({n_ok}/{len(results)})")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
