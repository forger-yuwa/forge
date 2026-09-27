#!/usr/bin/env python3
r"""V-ax2 同心円環 (ソルバ内連成) の合否を当てる。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/accepted/boundary-cht-axisymmetric-fem2d.md) §6 **V-ax2**。
**閾値・解析解は登録値を引数の既定に固定してあり、結果を見て動かさない**。

## 解析解 (単位軸長・ラジアンあたり、登録値)

    R'_f = ln(r_b/r_a)/k_f = 28.761,  R'_s = ln(r_c/r_b)/k_s = 25.672,  R'_h = 1/(h r_c) = 5e-7   [m·K/(W·rad)]
    Q'   = (T_a − T_c)/ΣR' = 0.91855 W/(m·rad)
    T_w* = T_c + Q'(R'_s + R'_h) = 323.581 K   固体の温度降下 23.581 K (0.5 % = 0.118 K)
    q*   = Q'/r_b = 91.86 W/m²                 節点荷重の総和との比較は Q'L_x = 1.8371e-3 W/rad

## 合否 (本スクリプトが判定する項目。すべて全界面節点・最終 step)

    (a) max |T_w − T_w*|                              <= 0.5 % of 23.581 K    T_w = 壁ダンプ `Ts`
    (b) G-cons: 流体の内壁放熱・固体入熱 ΣQ_f・Robin 持ち去り Σq_hole の相互差 (最大) <= 0.5 % of Q'L_x
    (c) 連成保存: max_i |Q_sol,i − Q_f,i| / A_i^r       <= 0.1 % of q* (0.092 W/m²)
        Q_sol = (K u − b)_iface は**固体の物理作用素** `Fem2DOperator(axisym=True)` を固体ダンプの温度に当てたもの
        (case/58 `eval_v6p.py` と同じ作法)、Q_f = 共役壁ダンプの `iface_Qf_eff` [W/rad]、A_i^r = 固体側の集中量 ∫N_i r ds
    (d) `iface_q_eff` 全節点 max |q − q*|               <= 0.5 % of q*
    (e) `iface_q_eff_raw` 全節点 max |q − q*|           <= 0.5 % of q*
    (f) 準定常: 節点ログから全界面節点の (T_w − T_c) と q = Q_f,i/A_i^r の毎更新系列を `vax2_series.csv` に書き、
        `check_quasisteady.py --series-csv ... --tail 0.5 --drift 0.001 --osc 0.001` を呼ぶ (**絶対温度は渡さない**)

**流体の内壁放熱の取り方**: 内壁 (physID 3、等温 350 K、非連成) の壁ダンプ `iface_Qf_eff` は「流体 → 壁」を正とする
積分済み荷重 [W/rad] なので、熱い内壁が流体へ出す熱は $Q_{\rm in}=-\sum_i$ `iface_Qf_eff`。
共役壁 (physID 4) は $\sum_i$ `iface_Qf_eff` (流体 → 固体、正)、Robin は固体ダンプの $\sum$ `q_hole` (固体 → 冷却側、正)。

**本スクリプトが判定しない登録項目** (別ツール。総合判定には全部要る):

    G-if     python3 solver_density_cuda/tools/check_cht_interface.py <run>          (conjugate_gate.json の登録値)
    流体残差 python3 solver_density_cuda/tools/check_convergence.py <run>           (一様 IC からの単一区間)
    感度     python3 case/61.conjugate_annulus/sens_vax2.py ...

負例 (plan §6 V-ax2 / V-ax2b の M2): `--simulate-old-denominator` は `iface_q_eff`・`iface_q_eff_raw` を
「分母修正を外した出力」$q_{\rm bad}=q_{\rm eff}\,A^r_{\rm fluid}/A_{\rm planar}$ に置き換えて判定する (値は r_b 倍 = 0.9186 W/m²)。
(d)(e) が FAIL しなければ評価器に検出力が無い。試験は `test_eval_vax2.py`。

使い方:
  python3 case/61.conjugate_annulus/eval_vax2.py <run_dir> [--step N]
結果は <run>/vax2_eval.json にも書く (項目ごとの値・許容・判定)。
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

import axcht

# ---- 登録値 (plan §6 V-ax2) ----
RA, RB, RC, LX = 5.0e-3, 10.0e-3, 20.0e-3, 2.0e-3
KF, KS, H, TA, TC = 0.0241, 0.027, 1.0e8, 350.0, 300.0
PID_IN, PID_CJ = 3, 4
NAME_IN, NAME_CJ = "wall_inner", "wall_outer"


def analytic():
    Rf = math.log(RB / RA) / KF
    Rs = math.log(RC / RB) / KS
    Rh = 1.0 / (H * RC)
    Qp = (TA - TC) / (Rf + Rs + Rh)
    Tw = TC + Qp * (Rs + Rh)
    return {"Rf": Rf, "Rs": Rs, "Rh": Rh, "Qp": Qp, "Tw": Tw, "drop": Tw - TC, "q": Qp / RB, "QL": Qp * LX}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--tol-tw", type=float, default=0.5, help="(a) [%% of 降下] (登録値)")
    ap.add_argument("--tol-gcons", type=float, default=0.5, help="(b) [%% of Q'L_x] (登録値)")
    ap.add_argument("--tol-cons", type=float, default=0.1, help="(c) [%% of q*] (登録値)")
    ap.add_argument("--tol-q", type=float, default=0.5, help="(d)(e) [%% of q*] (登録値)")
    ap.add_argument("--no-qs", action="store_true", help="準定常 (f) を呼ばない (1〜2 step の試験用)")
    ap.add_argument("--simulate-old-denominator", action="store_true",
                    help="負例: q_eff を分母修正前 (Q/A_planar) の値に置き換えて判定する")
    a = ap.parse_args(argv)
    run = Path(a.run.rstrip("/"))
    step = a.step if a.step is not None else axcht.last_step(run)
    an = analytic()
    rfloor = axcht.axis_r_floor(run)

    # ---- 流体の壁 (r 重み面積は mesh.h5 の PLANES から、ソルバと同じ定義)
    xyz_cj, Apl_cj, Ar_cj, _ = axcht.fluid_wall_areas(run, PID_CJ, rfloor)
    xyz_in, _, _, _ = axcht.fluid_wall_areas(run, PID_IN, rfloor)
    wcj = axcht.wall_dump(run, NAME_CJ, PID_CJ, step)
    win = axcht.wall_dump(run, NAME_IN, PID_IN, step)
    axcht.check_wall_order(wcj["xyz"], xyz_cj, NAME_CJ)
    axcht.check_wall_order(win["xyz"], xyz_in, NAME_IN)
    q_eff, q_raw = wcj["iface_q_eff"], wcj["iface_q_eff_raw"]
    if a.simulate_old_denominator:
        q_eff = q_eff * Ar_cj / Apl_cj
        q_raw = q_raw * Ar_cj / Apl_cj

    # ---- 固体
    sol = axcht.solid_state(run, PID_CJ, step)
    jw = axcht.match_iface_to_wall(sol["iface_xy"], wcj["xyz"])      # 界面節点 i -> 壁節点
    Qf_if = wcj["iface_Qf_eff"][jw]
    Tw = wcj["Ts"]

    rows = []
    ok = []
    print(f"=== V-ax2 評価: {run}  step {step}  (出力の精度 {wcj['dtype']}) ===")
    print(f"解析解: Q' = {an['Qp']:.6f} W/(m·rad), T_w* = {an['Tw']:.4f} K, 降下 {an['drop']:.4f} K, "
          f"q* = {an['q']:.4f} W/m², Q'L_x = {an['QL']:.6e} W/rad")
    # (a)
    ea = float(np.max(np.abs(Tw - an["Tw"]))) / an["drop"] * 100
    ok.append(axcht.gate(rows, "(a) max|T_w - T_w*| [% of 降下]", ea, a.tol_tw, "%"))
    # (b)
    Q_in = -float(np.sum(win["iface_Qf_eff"]))
    Q_f = float(np.sum(wcj["iface_Qf_eff"]))
    Q_h = float(np.sum(sol["q_hole"]))
    trio = [Q_in, Q_f, Q_h]
    eb = max(abs(x - y) for x in trio for y in trio) / an["QL"] * 100
    ok.append(axcht.gate(rows, "(b) G-cons 相互差の最大 [% of Q'L_x]", eb, a.tol_gcons, "%"))
    # (c)
    ec = float(np.max(np.abs(sol["Qsol_if"] - Qf_if) / sol["A_r"])) / an["q"] * 100
    ok.append(axcht.gate(rows, "(c) 連成保存 max|Q_sol-Q_f|/A^r [% of q*]", ec, a.tol_cons, "%"))
    # (d)(e)
    ed, _, _ = axcht.q_eff_check(q_eff, an["q"])
    ee, _, _ = axcht.q_eff_check(q_raw, an["q"])
    ok.append(axcht.gate(rows, "(d) iface_q_eff max|q - q*| [% of q*]", ed / an["q"] * 100, a.tol_q, "%"))
    ok.append(axcht.gate(rows, "(e) iface_q_eff_raw max|q - q*| [% of q*]", ee / an["q"] * 100, a.tol_q, "%"))

    # (f) 準定常 (全界面節点)
    qs = None
    if not a.no_qs:
        try:
            cols, nup = axcht.write_series_csv(run, PID_CJ, TC, run / "vax2_series.csv")
        except ValueError as e:
            print(f"節点ログが使えない: {e}")
            cols, nup = None, 0
        if cols is None:
            qs = "REFUSED"
        else:
            cmd = [sys.executable, str(axcht.TOOLS / "check_quasisteady.py"), "--series-csv",
                   str(run / "vax2_series.csv"), "--series-cols", ",".join(cols),
                   "--tail", "0.5", "--drift", "0.001", "--osc", "0.001"]
            p = subprocess.run(cmd, capture_output=True, text=True)
            txt = p.stdout + p.stderr
            (run / "QUASISTEADY_vax2.txt").write_text(" ".join(cmd) + "\n" + txt)
            head = [l for l in txt.splitlines() if l.startswith("===")]
            qs = head[-1].split("->")[-1].strip(" =") if head else "(判定なし)"
            print(f"\n準定常 ({nup} 更新 × {len(cols)//2} 節点、T_w−T_c と q): {qs}  -> {run}/QUASISTEADY_vax2.txt")
            print("  " + " ".join(cmd))
        ok.append(qs == "STEADY")
        rows.append(("(f) 準定常 (全界面節点, check_quasisteady)", float("nan"), float("nan"), qs, qs == "STEADY", "{}"))

    axcht.print_rows([r for r in rows if not r[0].startswith("(f)")])
    if qs is not None:
        print(f"{'(f) 準定常 (全界面節点)':<56}{qs:>26}  {'PASS' if qs == 'STEADY' else '**FAIL**'}")

    # ---- 参考 (判定に入れない)
    Tsol_if = sol["T"][sol["iface"]]
    print("\n参考:")
    print(f"  T_w 全節点 {Tw.min():.5f} .. {Tw.max():.5f} K (解析 {an['Tw']:.5f})")
    print(f"  固体界面温度 − 壁 Ts  max {np.max(np.abs(Tsol_if - Tw[jw])):.3e} K")
    print(f"  Q_in = {Q_in:.6e}, ΣQ_f = {Q_f:.6e}, Σq_hole = {Q_h:.6e} W/rad (Q'L_x {an['QL']:.6e}; "
          f"各々の解析差 {', '.join(f'{(x/an["QL"]-1)*100:+.4f} %' for x in trio)})")
    print(f"  iface_q_eff {q_eff.min():.5f} .. {q_eff.max():.5f}, raw {q_raw.min():.5f} .. {q_raw.max():.5f} W/m²"
          + ("  **(分母修正を外した模擬値)**" if a.simulate_old_denominator else ""))
    print(f"  流体側 A^r 合計 {Ar_cj.sum():.6e} (r_b L_x = {RB*LX:.6e}) / 固体側 A^r 合計 {sol['A_r'].sum():.6e} m²/rad")

    v = "PASS" if all(ok) else "FAIL"
    bad = [r[0] for r in rows if not r[4]]
    print(f"\nVERDICT: {v}" + (f"  (外れた項目: {', '.join(bad)})" if bad else "")
          + "   ※ G-if・流体残差・感度は別ツール (docstring)")
    res = {"run": str(run), "step": step, "simulated_old_denominator": a.simulate_old_denominator,
           "analytic": an, "items": [{"name": r[0], "value": r[1], "tol": r[2], "pass": r[4]} for r in rows],
           "verdict": v}
    out = run / ("vax2_eval_simold.json" if a.simulate_old_denominator else "vax2_eval.json")
    out.write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
