#!/usr/bin/env python3
r"""V-ax2b 同軸円板 (ソルバ内連成、半径が変わる界面の局所整合) の合否を当てる。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/accepted/boundary-cht-axisymmetric-fem2d.md) §6 **V-ax2b**。
**閾値・解析解は登録値を引数の既定に固定してあり、結果を見て動かさない**。

## 解析解 (解は $x$ の 1 次関数で $r$ に依らない)

    q    = (350 − 300) / (H/k_f + t/k_s + 1/h)     H = t = 5 mm, k_f 0.0241, k_s 0.027, h 1e8
    T_w* = T_c + q (t/k_s + 1/h)                   固体の温度降下 = T_w* − T_c
    Q_tot = q (r_2² − r_1²)/2                      [W/rad] (G-cons の規格化。V-ax2 の Q'L_x に当たる)
(登録文の式は 1/h を省いているが、h = 1e8 の寄与は T_w* で 1e-6 K 程度。ここでは含める)

## 合否 (本スクリプトが判定する項目。すべて全界面節点 = 端点を含む・最終 step)

    (a) max |T_w − T_w*|                                <= 0.5 % of 降下
        (登録の合格は N_r=32 一様と非一様に対して。8 / 16 は系列判定 `series_vax2b.py` の素材)
    (b) G-cons: 熱い壁 (x=0) の放熱・固体入熱 ΣQ_f・Robin 持ち去り Σq_hole の相互差 (最大) <= 0.5 % of Q_tot
    (c) 連成保存: max_i |Q_sol,i − Q_f,i| / A_i^r         <= 0.1 % of q
    (d) 恒等式: max_i |q_eff,i A^r_fluid,i − Q_f,eff,i| / |Q_f,eff,i| <= 1e-12 (**FP64 run のみ**。FP32 は判定不能)
        あわせて iface_q_eff_raw A^r = R_raw − F_w (ソルバの定義) も同じ許容で照合する
    (e) 準定常: 全界面節点の (T_w − T_c) と q = Q_f,i/A_i^r の毎更新系列 → `check_quasisteady.py --tail 0.5 --drift 0.001 --osc 0.001`

参考 (判定に入れない): `iface_q_eff` と q の差 (流体面の厳密な面積分なので一様熱流束なら q に一致するはず)、
固体側 A_i^r で割った Q_f/A_i^r の端点の偏り (plan §3: 一様熱流束でも両端で 0.9375 q / 1.05 q の型)。

**本スクリプトが判定しない登録項目**: G-if (`check_cht_interface.py <run>`)、流体残差 (`check_convergence.py <run>`)、
系列の収束 (`series_vax2b.py`)。

`--identity-only` は軸に触れる負例 run (template/*_axisprobe.yaml、共役なし 1 step) 用: 両壁の全節点で恒等式を照合し、
`axisRFloor` の床に当たった節点が 1 つ以上あること、**床を外した面積では恒等式が崩れること** (検出力) を見る。

使い方:
  python3 case/62.conjugate_disk/eval_vax2b.py <run_dir> [--step N]
  python3 case/62.conjugate_disk/eval_vax2b.py <axisprobe_run> --identity-only
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "61.conjugate_annulus"))
import axcht  # noqa: E402

# ---- 登録値 (plan §6 V-ax2b) ----
H, T, R1, R2 = 5.0e-3, 5.0e-3, 5.0e-3, 20.0e-3
KF, KS, HR, THOT, TC = 0.0241, 0.027, 1.0e8, 350.0, 300.0
PID_HOT, PID_CJ = 3, 4
NAME_HOT, NAME_CJ = "wall_hot", "wall_cj"
TOL_ID = 1.0e-12


def analytic():
    Rs = T / KS + 1.0 / HR
    q = (THOT - TC) / (H / KF + Rs)
    Tw = TC + q * Rs
    return {"q": q, "Tw": Tw, "drop": Tw - TC, "Qtot": q * (R2 ** 2 - R1 ** 2) / 2.0}


def identity_items(run, step, rfloor, pids):
    rows, ok = [], True
    for pid, name in pids:
        r = axcht.identity_check(run, pid, name, step, rfloor)
        fp64 = r["dtype"] == "float64"
        for key, lab in (("rel_q_eff", "q_eff·A^r = Q_f,eff"), ("rel_q_raw", "q_eff_raw·A^r = R_raw−F_w")):
            v = r[key]
            if v is None:
                rows.append((f"恒等式 {name} {lab}", "列なし", False)); ok = False
                continue
            p = fp64 and v <= TOL_ID
            ok &= p
            rows.append((f"恒等式 {name} {lab} [相対 max]", f"{v:.3e}" + ("" if fp64 else f" ({r['dtype']}: 判定不能)"), p))
        rows.append((f"  ({name}: {r['n']} 節点, 床 {r['rfloor']:g} m に当たった面 {r['n_floor']}, "
                     f"面重心 r 最小 {r['r_face_min']:.4e} m)", "", None))
    return rows, ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--identity-only", action="store_true", help="軸に触れる負例 run (共役なし) の恒等式だけ")
    ap.add_argument("--tol-tw", type=float, default=0.5, help="(a) [%% of 降下] (登録値)")
    ap.add_argument("--tol-gcons", type=float, default=0.5, help="(b) [%% of Q_tot] (登録値)")
    ap.add_argument("--tol-cons", type=float, default=0.1, help="(c) [%% of q] (登録値)")
    ap.add_argument("--no-qs", action="store_true", help="準定常 (e) を呼ばない (試験用)")
    a = ap.parse_args(argv)
    run = Path(a.run.rstrip("/"))
    step = a.step if a.step is not None else axcht.last_step(run)
    rfloor = axcht.axis_r_floor(run)
    an = analytic()

    if a.identity_only:
        print(f"=== V-ax2b 軸負例 (恒等式のみ): {run}  step {step}  axisRFloor {rfloor:g} m ===")
        rows, ok = identity_items(run, step, rfloor, [(PID_HOT, NAME_HOT), (PID_CJ, NAME_CJ)])
        nfloor = sum(axcht.identity_check(run, p, n, step, rfloor)["n_floor"] for p, n in
                     [(PID_HOT, NAME_HOT), (PID_CJ, NAME_CJ)])
        # 検出力: 床を外した面積 |S_f| r_f では床に当たった節点で恒等式が崩れるはず
        nofloor = [axcht.identity_check(run, p, n, step, rfloor, use_floor=False)["rel_q_eff"]
                   for p, n in [(PID_HOT, NAME_HOT), (PID_CJ, NAME_CJ)]]
        for nm, v, p in rows:
            print(f"{nm:<70}{v:>28}  " + ("" if p is None else ("PASS" if p else "**FAIL**")))
        print(f"床に当たった面 {nfloor} (1 以上であること)")
        print(f"対照: 床を外した面積での恒等式 相対 max = {', '.join(f'{v:.3e}' for v in nofloor)} "
              f"(床に当たった節点があれば 1e-12 を大きく超えるはず)")
        power = nfloor > 0 and min(nofloor) > 1e-6
        v = "PASS" if (ok and power) else "FAIL"
        dts = {axcht.wall_dump(run, n, p, step)["dtype"] for p, n in [(PID_HOT, NAME_HOT), (PID_CJ, NAME_CJ)]}
        if dts != {"float64"}:
            v = "REFUSED"
            print(f"出力が {sorted(dts)} — 相対 1e-12 の恒等式は FP64 run でしか判定できない")
        if nfloor == 0:
            v = "REFUSED"
            print("床に当たった面が無い — この run は床の恒等式を試していない")
        print(f"\nVERDICT: {v}")
        return 0 if v == "PASS" else 1

    xyz_cj, _, Ar_cj, _ = axcht.fluid_wall_areas(run, PID_CJ, rfloor)
    xyz_h, _, _, _ = axcht.fluid_wall_areas(run, PID_HOT, rfloor)
    wcj = axcht.wall_dump(run, NAME_CJ, PID_CJ, step)
    wh = axcht.wall_dump(run, NAME_HOT, PID_HOT, step)
    axcht.check_wall_order(wcj["xyz"], xyz_cj, NAME_CJ)
    axcht.check_wall_order(wh["xyz"], xyz_h, NAME_HOT)
    sol = axcht.solid_state(run, PID_CJ, step)
    jw = axcht.match_iface_to_wall(sol["iface_xy"], wcj["xyz"])
    Qf_if = wcj["iface_Qf_eff"][jw]
    Tw = wcj["Ts"]
    r_w = wcj["xyz"][:, 1]

    print(f"=== V-ax2b 評価: {run}  step {step}  (出力の精度 {wcj['dtype']}, 界面 {len(Tw)} 節点) ===")
    print(f"解析解: q = {an['q']:.5f} W/m², T_w* = {an['Tw']:.5f} K, 降下 {an['drop']:.5f} K, Q_tot = {an['Qtot']:.6e} W/rad")
    rows, ok = [], []
    ea = np.abs(Tw - an["Tw"]) / an["drop"] * 100
    ok.append(axcht.gate(rows, "(a) max|T_w - T_w*| [% of 降下] (端点を含む)", float(ea.max()), a.tol_tw, "%"))
    Q_hot = -float(np.sum(wh["iface_Qf_eff"]))
    Q_f = float(np.sum(wcj["iface_Qf_eff"]))
    Q_h = float(np.sum(sol["q_hole"]))
    trio = [Q_hot, Q_f, Q_h]
    eb = max(abs(x - y) for x in trio for y in trio) / an["Qtot"] * 100
    ok.append(axcht.gate(rows, "(b) G-cons 相互差の最大 [% of Q_tot]", eb, a.tol_gcons, "%"))
    ec = float(np.max(np.abs(sol["Qsol_if"] - Qf_if) / sol["A_r"])) / an["q"] * 100
    ok.append(axcht.gate(rows, "(c) 連成保存 max|Q_sol-Q_f|/A^r [% of q]", ec, a.tol_cons, "%"))
    axcht.print_rows(rows)
    idr, idok = identity_items(run, step, rfloor, [(PID_CJ, NAME_CJ)])
    for nm, v, p in idr:
        print(f"{'(d) ' + nm if p is not None else nm:<70}{v:>22}  " + ("" if p is None else ("PASS" if p else "**FAIL**")))
    ok.append(idok)

    qs = None
    if not a.no_qs:
        try:
            cols, nup = axcht.write_series_csv(run, PID_CJ, TC, run / "vax2b_series.csv")
        except ValueError as e:
            print(f"節点ログが使えない: {e}")
            cols = None
        if cols is None:
            qs = "REFUSED"
        else:
            cmd = [sys.executable, str(axcht.TOOLS / "check_quasisteady.py"), "--series-csv",
                   str(run / "vax2b_series.csv"), "--series-cols", ",".join(cols),
                   "--tail", "0.5", "--drift", "0.001", "--osc", "0.001"]
            p = subprocess.run(cmd, capture_output=True, text=True)
            txt = p.stdout + p.stderr
            (run / "QUASISTEADY_vax2b.txt").write_text(" ".join(cmd) + "\n" + txt)
            head = [l for l in txt.splitlines() if l.startswith("===")]
            qs = head[-1].split("->")[-1].strip(" =") if head else "(判定なし)"
        print(f"{'(e) 準定常 (全界面節点、T_w−T_c と q)':<70}{qs:>22}  {'PASS' if qs == 'STEADY' else '**FAIL**'}")
        ok.append(qs == "STEADY")

    # ---- 参考
    o = np.argsort(r_w)
    print("\n参考 (判定に入れない):")
    print(f"  T_w 誤差 [% of 降下] r 昇順: " + " ".join(f"{v:.4f}" for v in ea[o]))
    print(f"  iface_q_eff − q: max {np.max(np.abs(wcj['iface_q_eff'] - an['q'])):.4e} W/m² "
          f"({np.max(np.abs(wcj['iface_q_eff'] - an['q'])) / an['q'] * 100:.4f} % of q)、"
          f"raw {np.max(np.abs(wcj['iface_q_eff_raw'] - an['q'])) / an['q'] * 100:.4f} %")
    qs_solid = Qf_if / sol["A_r"]
    os_ = np.argsort(sol["iface_xy"][:, 1])
    print(f"  Q_f/A_i^r (固体側の集中量で割った値) / q: r_1 端 {qs_solid[os_[0]]/an['q']:.5f}, "
          f"中央 {qs_solid[os_[len(os_)//2]]/an['q']:.5f}, r_2 端 {qs_solid[os_[-1]]/an['q']:.5f}")
    print(f"  固体界面温度 − 壁 Ts  max {np.max(np.abs(sol['T'][sol['iface']] - Tw[jw])):.3e} K")
    print(f"  Q_hot = {Q_hot:.6e}, ΣQ_f = {Q_f:.6e}, Σq_hole = {Q_h:.6e} W/rad (Q_tot {an['Qtot']:.6e})")

    v = "PASS" if all(ok) else "FAIL"
    print(f"\nVERDICT: {v}   ※ G-if・流体残差・系列は別ツール (docstring)")
    res = {"run": str(run), "step": step, "analytic": an, "dtype": wcj["dtype"],
           "Tw_err_pct_max": float(ea.max()), "Tw_err_pct_by_r": [[float(r_w[i]), float(ea[i])] for i in o],
           "items": [{"name": r[0], "value": r[1], "tol": r[2], "pass": r[4]} for r in rows],
           "identity_pass": bool(idok), "quasisteady": qs, "verdict": v}
    (run / "vax2b_eval.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
