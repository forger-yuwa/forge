"""plan tooling-nozzle-throat-monotone-r2 §6 E′ の実務判定 (ユーザ決定 2026-10-06「まずは 1 で」)。
腕 A (現行壁) run_0140〜0142 と腕 B (単調壁) run_0143〜0145 の評価量を、本段 step 6000〜18000 の 13 枚の時間平均で比べる。
入力: 各 run の wallfit_series_e3.csv (eval_wallfit_euler.py --e3 が書く)。予備 A/B の β−α の差 (番号写像 − 最近傍)は wallfit_series_icab.csv から併記する。
統計 (§6 E′):
  run i の窓平均 m_i、時間変動の標準誤差 s_i = sd_i/√n (自己相関は補正しない)。
  腕の平均 M = mean(m_i)、腕の標準誤差 SE = max(sd(m_i)/√3, √(Σ s_i²)/3)。
  D = M_B − M_A、SE_D = √(SE_A² + SE_B²)。
前提検査 (codex result 段レビュー 2026-10-07 M2。どれか不成立なら総合は「保留 (前提不成立)」):
  残差: 判定区間 (--segment) の判定が pass か既知の plateau。DIVERGED・RISING・欠損・判定不能は不成立。
  壁: prepare_info の wall_fit.mono_r2 が腕 A は None、腕 B は [0.0, 1.5] と完全一致。メッシュの壁節点が自腕の当てはめ後 spline と一致し
      (throat_mono_ab.wall_evidence の consistent)、対の他腕 (A rk ↔ B rk) との照合が ok。
  初期値: 腕 B は IC_MAP.json の VERDICT OK・mode index で、3 本の写像後 sha256 が同じ (同じ準備から作った)。
          腕 A は新形式の IC 記録より前に準備したので、旧形式の例外として restart_field.log の「SRC とビット一致」を要する (例外として記録)。
  乾式確認 (DRY) の準備から作った run は不成立。
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


def spline_shape_checks(info: dict, arm: str) -> list:
    """保存された当てはめ後 spline (prepare_info の wall_fit.spline) 自体の形状検査 (2 回目の result 段レビュー M1)。
    腕 B: 形状ゲート S1 と同じ形状用許容差で、[0, 1.5] の r‴ 最大 ≤ 1e-6 かつ r″ の最大増加 ≤ 1e-7。
    腕 A: 単調拘束が掛かっていない旧壁であること (r″ の最大増加 > 1e-4; 旧壁は 0.0102)。"""
    from scipy.interpolate import BSpline
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
    from forge_design.geometry.wall_axismach import r3_piecewise_exact
    sp = (info.get("wall_fit") or {}).get("spline")
    if not sp:
        return ["prepare_info の wall_fit.spline が無い (形状を検査できない)"]
    spl = BSpline(np.asarray(sp["t"], dtype=float), np.asarray(sp["c"], dtype=float), int(sp["k"]))
    sh = r3_piecewise_exact(spl, None, 0.0, 1.5)
    if arm == "B" and not (sh["r3_max"] <= 1e-6 and sh["r2_max_increase"] <= 1e-7):
        return [f"腕 B の保存 spline が単調でない (r‴ 最大 {sh['r3_max']:.3g}、r″ の最大増加 {sh['r2_max_increase']:.3g})"]
    if arm == "A" and not sh["r2_max_increase"] > 1e-4:
        return [f"腕 A の保存 spline が旧壁に見えない (r″ の最大増加 {sh['r2_max_increase']:.3g} ≤ 1e-4)"]
    return []


def preconditions() -> tuple:
    """§冒頭の前提検査。戻り値 (不成立の一覧 {run: [理由]}, 例外の記録 [文])。証拠が読めないこと自体を不成立にする。"""
    from throat_mono_judge import SEGMENT_VERDICT_FILE, parse_segment_verdict, mono_r2_matches
    bad, exceptions, sha_b = {}, [], {}
    for arm, runs in ARMS.items():
        for k, run in enumerate(runs):
            rd = C / run; why = []
            try:
                seg = parse_segment_verdict((rd / SEGMENT_VERDICT_FILE).read_text() if (rd / SEGMENT_VERDICT_FILE).is_file() else None)
                if seg["status"] not in ("pass", "plateau"):
                    why.append(f"判定区間の残差判定が {seg['status']} ({seg.get('reason') or seg.get('line')})")
                info = json.loads((rd / "prepare_info.json").read_text())
                if info.get("DRY") or info.get("DRY_NO_IC"):
                    why.append("乾式確認 (DRY) の準備から作った run")
                wf = info.get("wall_fit") or {}
                exp = None if arm == "A" else [0.0, 1.5]
                if "mono_r2" not in wf:
                    why.append("prepare_info の wall_fit に mono_r2 の記録が無い")
                elif not mono_r2_matches(wf.get("mono_r2"), exp):
                    why.append(f"mono_r2 {wf.get('mono_r2')!r} が期待 {exp!r} と一致しない")
                from throat_mono_ab import wall_evidence
                ev = wall_evidence(rd, C / ARMS["B" if arm == "A" else "A"][k])
                if ev.get("status") != "consistent":
                    why.append(f"壁の証拠が {ev.get('status')} ({ev.get('reason', '')})")
                vo = ev.get("vs_other") or {}
                if vo.get("status") != "ok":
                    why.append(f"他腕との壁の照合が {vo.get('status')}")
                else:
                    # 2 回目の result 段レビュー M1: 同じ壁を両腕に使っても照合は ok になるので、見分けられる節点があり、
                    # それがすべて自腕の当てはめに一致することを必須にする
                    nd, nown = int(vo.get("n_discriminable", 0)), int(vo.get("n_discriminable_matching_own", -1))
                    if nd <= 0:
                        why.append("他腕と見分けられる壁節点が無い (両腕が同じ壁の可能性)")
                    elif nown != nd:
                        why.append(f"見分けられる壁節点 {nd} のうち自腕の当てはめに一致するのは {nown} (壁の取り違えの可能性)")
                why += spline_shape_checks(info, arm)
                if arm == "B":
                    rec = json.loads((rd / "IC_MAP.json").read_text())
                    if rec.get("VERDICT") != "OK" or rec.get("mode") != "index":
                        why.append(f"IC_MAP.json の VERDICT {rec.get('VERDICT')!r}・mode {rec.get('mode')!r} (OK・index であること)")
                    sha_b[run] = rec.get("dst_sha256_after")
                else:
                    log = (rd / "restart_field.log").read_text() if (rd / "restart_field.log").is_file() else ""
                    if "SRC とビット一致" not in log:
                        why.append("restart_field.log に「SRC とビット一致」が無い (腕 A の IC の証拠)")
                    else:
                        exceptions.append(f"{run}: 新形式の IC 記録より前に準備したため、旧形式 (restart_field.log のビット一致) で受けた")
            except Exception as e:  # noqa: BLE001  証拠が読めないことは不成立
                why.append(f"証拠を読めない: {type(e).__name__}: {e}")
            if why:
                bad[run] = why
    if len(set(sha_b.values())) > 1 or any(v is None for v in sha_b.values()):
        for r in sha_b:
            bad.setdefault(r, []).append(f"腕 B の 3 本の写像後 sha256 がそろわない ({sorted(set(map(str, sha_b.values())))})")
    return bad, exceptions


def main():
    pre_bad, pre_exc = preconditions()
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
    if pre_bad:
        overall = "保留 (前提不成立): " + "; ".join(f"{r}: {', '.join(w)}" for r, w in pre_bad.items())
    elif all(v == "許容幅未満" for v in verdicts.values()):
        det = [r["qty"] for r in rows if r.get("detected")]
        overall = ("単調壁を候補形状として採用 (Euler で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし"))
    elif any(v == "悪化" for v in verdicts.values()):
        overall = "不採用 (悪化した量がある)"
    else:
        overall = "保留 (ユーザ判断)"
    out = dict(plan="plans/accepted/tooling-nozzle-throat-monotone-r2.md §6 E′", window_steps=WIN, n_per_run=13,
               arms=ARMS, icab=ICAB, rows=rows, overall=overall, preconditions_failed=pre_bad, precondition_exceptions=pre_exc,
               limits=["自己相関は補正していない", "窓の開始 6000 は予備 A/B の時系列を見た後に決めた",
                       "IC 写像による差 (β−α = 番号写像 − 最近傍) は Δq/10 の精度では除外できていない (表の ic 列)", "軸 (η0) の量は判定対象外"])
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / "_band_ab/throat_mono_practical_eval.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(f"{'量':28s} {'Δq':>8s} {'A 平均':>11s} {'B 平均':>11s} {'D=B−A':>10s} {'2·SE_D':>9s} {'(D+2SE)/Δq':>10s} {'IC β−α':>10s} {'検出':>4s}  判定")
    for r in rows:
        print(f"{r['qty']:28s} {r.get('dq', float('nan')):8.2e} {r['A']['mean']:11.6g} {r['B']['mean']:11.6g} {r['D']:+10.2e} "
              f"{2 * r['SE_D']:9.2e} {r.get('D_plus_2SE_over_dq', float('nan')):10.2f} {r['ic_beta_minus_alpha']:+10.2e} {('あり' if r.get('detected') else 'なし') if 'detected' in r else '-':>4s}  {r['verdict']}")
    for e in pre_exc:
        print("前提の例外:", e)
    print("総合:", overall)


if __name__ == "__main__":
    main()
