#!/usr/bin/env python3
"""#4f (plan gradient-scalar-lsq-unification §5.1 #4f、測る前に固定): tgv lsq のダンプ無し・有りを各 +5 本 (d–h)、ABBA 順で 1 step。
入力は既存 `pg_tgv_lsq_a` の入力ファイルの複製 (出力は持ち込まない)。バイナリは既存 12 本と同じ。
使い方 (AWS): python3 lsq8_run.py --src ~/sglsq/pg/pg_tgv_lsq_a --scratch ~/sglsq/pg8 --bin ~/sglsq/forge_hook4a"""
import argparse, os, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G  # noqa: E402

ORDER = [("d", "off"), ("d", "on"), ("e", "on"), ("e", "off"), ("f", "off"), ("f", "on"),
         ("g", "on"), ("g", "off"), ("h", "off"), ("h", "on")]
INPUTS = ("solverConfig.yaml", "bcondConfig.yaml", "mesh.h5", "probe.yaml")

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True); ap.add_argument("--scratch", required=True); ap.add_argument("--bin", required=True)
a = ap.parse_args()
os.makedirs(a.scratch, exist_ok=True)
for rep, side in ORDER:
    d = os.path.join(a.scratch, f"lsq_{side}_{rep}")
    if os.path.exists(os.path.join(d, "res_1.h5")):
        print("exists", d); continue
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in INPUTS:
        if os.path.exists(os.path.join(a.src, f)):
            shutil.copy2(os.path.join(a.src, f), d)
    if side == "on":
        os.environ["FORGE_DUMP_PREGATHER"] = os.path.join(d, "pregather")
    G.run_forge(d, expect_merge=True, forge_bin=a.bin)
    os.environ.pop("FORGE_DUMP_PREGATHER", None)
    print("ran", d, flush=True)
