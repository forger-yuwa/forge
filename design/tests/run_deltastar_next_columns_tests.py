#!/usr/bin/env python3
"""`delta_r_next.csv` の列の区別 (δ_target = 緩和後の壁入力 / δ_E = 未緩和の平滑化抽出) の試験。

plan tooling-nozzle-cfd-pinned-initial-line §5.1 #8 (codex result M1: solve_rt・c2pin_solve・c2pin_hand0 が第 2 列 [緩和後]
を δ_E として読んでいた取り違えの訂正)。
  1. 既知の CSV (新しい列名 / 旧列名) で `read_delta_r_next` が δ_E = 第 4 列、δ_target = 第 2 列を返す。想定外の列は拒否。
  2. `solve_rt(prev_run=...)` が δ_E を読み、どの列を読んだかを返り値に記録する。
  3. (run があれば) 既存 run の訂正値: x_F で δ_E/δ_C = run_0090 1.0072、run_0092 1.0038。

usage: design/.venv-opt/bin/python design/tests/run_deltastar_next_columns_tests.py
  既存 run の場所: $FORGE_CASE45 (既定 /home/sano/work/forge/case/45.isobutane_m6_d155)。無ければ 3 は省略 (表示する)。
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.feedback.deltastar_loop import (DELTA_R_NEXT_COLUMNS, read_delta_r_next,  # noqa: E402
                                                   solve_rt)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


tmp = Path(tempfile.mkdtemp(prefix="dnext_"))
x = np.linspace(-12.5, 100.0, 50)
d_target = 0.001 * (x + 13.0)          # 緩和後 (第 2 列)
d_prev = 0.0008 * (x + 13.0)
d_E = 0.0012 * (x + 13.0)              # 未緩和の抽出 (第 4 列)
held = np.zeros_like(x)
for kind, head in (("named", ",".join(DELTA_R_NEXT_COLUMNS) + " (omega=0.5; test)"),
                   ("legacy", "x_rt,delta_r,delta_in_prev,delta_r_use_smoothed,held (omega=0.5; test)")):
    f = tmp / f"{kind}.csv"
    np.savetxt(f, np.c_[x, d_target, d_prev, d_E, held], delimiter=",", comments="", header=head)
    t = read_delta_r_next(f)
    check(f"{kind}: delta_E = 第 4 列・delta_target = 第 2 列 (header_kind {t['header_kind']})",
          t["header_kind"] == kind and np.array_equal(t["delta_E"], d_E) and np.array_equal(t["delta_target"], d_target))
bad = tmp / "bad.csv"
np.savetxt(bad, np.c_[x, d_E, d_target, d_prev, held], delimiter=",", comments="", header="x_rt,delta_E,delta_target,a,b")
try:
    read_delta_r_next(bad)
    check("並びの違う CSV を拒否", False)
except ValueError:
    check("並びの違う CSV を拒否", True)

# --- solve_rt (prev_run あり) が δ_E を読む ----------------------------------------
PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_euler_c2final_n2400.yaml"
prev = tmp / "run_prev"
prev.mkdir()
S_prev = 0.0768
(prev / "prepare_info.json").write_text(json.dumps({"scale_m": S_prev}))
(prev / "delta_r_equiv.csv").write_text("x_rt,delta_r_raw\n0,0\n")
np.savetxt(prev / "delta_r_next.csv", np.c_[x, d_target, d_prev, d_E, held], delimiter=",", comments="",
           header=",".join(DELTA_R_NEXT_COLUMNS) + " (test)")
rs = solve_rt(PROB, 0.775, prev_run=prev)
dE_xF = float(np.interp(rs["x_F_rt"], x, d_E))
check(f"solve_rt: x_F の δ_E を読む ({rs['delta_measured_xF_rt']:.6f} vs δ_E {dE_xF:.6f}, δ_target "
      f"{float(np.interp(rs['x_F_rt'], x, d_target)):.6f})", abs(rs["delta_measured_xF_rt"] - dE_xF) < 1e-15)
check(f"solve_rt: 読んだ列を記録 ({rs['delta_column']})", str(rs["delta_column"]).startswith("delta_E"))
rt = rs["r_t_m"]
check("solve_rt: R = r_t (r_F + δ_E (r_t/S)^−0.2) を満たす",
      abs(rt * (rs["r_F_rt"] + dE_xF * (rt / S_prev) ** -0.2) - 0.775) < 1e-9)

# --- 既存 run の訂正値 ---------------------------------------------------------------
C45 = Path(os.environ.get("FORGE_CASE45", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
for run, ref in (("run_0090_ns_c2pin_final", 1.0072), ("run_0092_ns_c2pin_pass2", 1.0038)):
    rd = C45 / run
    if not (rd / "_extract_edge/delta_r_next.csv").exists():
        print(f"skip: {rd} が無い (訂正値の再現は省略)")
        continue
    nx = read_delta_r_next(rd / "_extract_edge/delta_r_next.csv")
    di = np.loadtxt(rd / "delta_r_initial.csv", delimiter=",", skiprows=1)
    xF = float(di[-1, 0])
    dC = float(np.interp(xF, di[:, 0], di[:, 1]))
    rE = float(np.interp(xF, nx["x_rt"], nx["delta_E"])) / dC
    rT = float(np.interp(xF, nx["x_rt"], nx["delta_target"])) / dC
    check(f"{run}: x_F で δ_E/δ_C = {rE:.4f} (訂正値 {ref}; 取り違えた δ_target/δ_C は {rT:.4f})", abs(rE - ref) < 5e-5)

print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
