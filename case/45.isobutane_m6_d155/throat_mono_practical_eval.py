"""plan tooling-nozzle-throat-monotone-r2 §6 E′ の実務判定 (ユーザ決定 2026-10-06「まずは 1 で」)。
腕 A (現行壁) run_0140〜0142 と腕 B (単調壁) run_0143〜0145 の評価量を、本段 step 6000〜18000 の 13 枚の時間平均で比べる。
入力: 各 run の wallfit_series_e3.csv (eval_wallfit_euler.py --e3 が書く)。予備 A/B の α−β の差は wallfit_series_icab.csv から併記する。
統計 (§6 E′):
  run i の窓平均 m_i、時間変動の標準誤差 s_i = sd_i/√n (自己相関は補正しない)。
  腕の平均 M = mean(m_i)、腕の標準誤差 SE = max(sd(m_i)/√3, √(Σ s_i²)/3)。
  D = M_B − M_A、SE_D = √(SE_A² + SE_B²)。
判定 (大きいほど悪い量の片側): D + 2·SE_D ≤ Δq → 「許容幅未満」/ D − 2·SE_D ≥ Δq → 「悪化」/ それ以外 → 保留 (ユーザ判断)。
  別欄で差の検出の有無 (|D| > 2·SE_D なら「検出」) を記録する。判定は許容幅で行い、検出の有無は判定に使わない。
usage: python3 throat_mono_practical_eval.py [case_dir]   → _band_ab/throat_mono_practical_eval.json と標準出力の表
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

C = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
ARMS = {"A": [f"run_0{n}_euler_wallfit_pinG1_r{k}" for n, k in ((140, 1), (141, 2), (142, 3))],
        "B": [f"run_0{n}_euler_wallfit_monoG1_r{k}" for n, k in ((143, 1), (144, 2), (145, 3))]}
ICAB = ("run_0143_euler_wallfit_monoG1_r1", "run_0146_euler_icab_monoG1_nn")   # β (番号写像), α (最近傍)
WIN = (6000, 18000)
DQ = {"M_wave_eta0.1": 0.001, "P_wave_eta0.1": 0.010, "overshoot_eta0.1": 0.003, "overshoot_exitnorm_eta0.1": 0.003,
      "P_slope_abs_eta0.1": 0.03, "exit_M_dev": 0.00018}
RECORD_ONLY = ("exit_core_M",)


def series(run: str, name: str) -> dict:
    f = C / run / name
    if not f.is_file():
        raise FileNotFoundError(f"{f} が無い (eval_wallfit_euler.py --e3 / --icab を先に回すこと)")
    rows = list(csv.DictReader(open(f)))
    win = [r for r in rows if WIN[0] <= int(float(r["step"])) <= WIN[1]]
    steps = [int(float(r["step"])) for r in win]
    if steps != list(range(WIN[0], WIN[1] + 1, 1000)):
        raise ValueError(f"{run}: 窓 {WIN} の 1000 step ごとの 13 枚がそろっていない (実際 {steps})")
    out = {}
    for k in list(DQ) + list(RECORD_ONLY):
        v = np.array([float(r[k]) for r in win])
        if not np.all(np.isfinite(v)):
            raise ValueError(f"{run}: {k} に非有限値")
        out[k] = v
    return out


def main():
    S = {arm: {r: series(r, "wallfit_series_e3.csv") for r in runs} for arm, runs in ARMS.items()}
    ic = {r: series(r, "wallfit_series_icab.csv") for r in ICAB}
    rows, verdicts = [], {}
    for k in list(DQ) + list(RECORD_ONLY):
        st = {}
        for arm, runs in S.items():
            m = np.array([runs[r][k].mean() for r in runs])
            s = np.array([runs[r][k].std(ddof=1) / np.sqrt(len(runs[r][k])) for r in runs])
            se = max(m.std(ddof=1) / np.sqrt(len(m)), float(np.sqrt((s ** 2).sum())) / len(m))
            st[arm] = dict(mean=float(m.mean()), run_means=[float(x) for x in m], se=float(se), temporal_se=[float(x) for x in s])
        D = st["B"]["mean"] - st["A"]["mean"]
        SE = float(np.hypot(st["A"]["se"], st["B"]["se"]))
        b, a = ic[ICAB[0]][k], ic[ICAB[1]][k]
        ic_diff = float(b.mean() - a.mean())
        ic_se = float(np.sqrt(b.var(ddof=1) / len(b) + a.var(ddof=1) / len(a)))
        row = dict(qty=k, A=st["A"], B=st["B"], D=float(D), SE_D=SE, ic_beta_minus_alpha=ic_diff, ic_se=ic_se)
        if k in DQ:
            dq = DQ[k]
            if D + 2 * SE <= dq:
                v = "許容幅未満"
            elif D - 2 * SE >= dq:
                v = "悪化"
            else:
                v = "保留"
            row.update(dq=dq, verdict=v, D_plus_2SE_over_dq=float((D + 2 * SE) / dq), detected=bool(abs(D) > 2 * SE))
            verdicts[k] = v
        else:
            row.update(verdict="記録のみ")
        rows.append(row)
    if all(v == "許容幅未満" for v in verdicts.values()):
        det = [r["qty"] for r in rows if r.get("detected")]
        overall = ("単調壁を候補形状として採用 (Euler で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし"))
    elif any(v == "悪化" for v in verdicts.values()):
        overall = "不採用 (悪化した量がある)"
    else:
        overall = "保留 (ユーザ判断)"
    out = dict(plan="plans/active/tooling-nozzle-throat-monotone-r2.md §6 E′", window_steps=WIN, n_per_run=13,
               arms=ARMS, icab=ICAB, rows=rows, overall=overall,
               limits=["自己相関は補正していない", "窓の開始 6000 は予備 A/B の時系列を見た後に決めた",
                       "IC 写像による差 (α−β) は Δq/10 の精度では除外できていない (表の ic 列)", "軸 (η0) の量は判定対象外"])
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / "_band_ab/throat_mono_practical_eval.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(f"{'量':28s} {'Δq':>8s} {'A 平均':>11s} {'B 平均':>11s} {'D=B−A':>10s} {'2·SE_D':>9s} {'(D+2SE)/Δq':>10s} {'IC β−α':>10s} {'検出':>4s}  判定")
    for r in rows:
        print(f"{r['qty']:28s} {r.get('dq', float('nan')):8.2e} {r['A']['mean']:11.6g} {r['B']['mean']:11.6g} {r['D']:+10.2e} "
              f"{2 * r['SE_D']:9.2e} {r.get('D_plus_2SE_over_dq', float('nan')):10.2f} {r['ic_beta_minus_alpha']:+10.2e} {('あり' if r.get('detected') else 'なし') if 'detected' in r else '-':>4s}  {r['verdict']}")
    print("総合:", overall)


if __name__ == "__main__":
    main()
