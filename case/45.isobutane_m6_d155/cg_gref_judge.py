"""plan tooling-nozzle-core-grid §4.18 (事前登録、codex 諮問 2026-10-11): 参照の格子 Gref との比較の集計と判定。

2 段に分ける (場のある場所と判定する場所が違うため):
  series: 1 本の run の出力 (通算 step の範囲) から登録した指標の系列を作る (場のあるインスタンスで回す。run には何も書かない)
      python3 cg_gref_judge.py series <名前> <run> [--offset N] [--from 160000] [--to 200000]
      → _band_ab/core_grid/gref_series_<名前>.csv (step は run_0183 の res_100000 からの通算 = run の step + offset)
  judge: 系列の CSV から §4.18 の式で判定する (どこでもよい)
      python3 cg_gref_judge.py judge [--tau-theta 1 ...] → _band_ab/core_grid/gref_judge.json

指標 (§4.5): θ_r・δ_loc = PCHIP でつないだ x の窓 [40,50]・[65,75]・[84,94] の平均、Q_w、出口の主流 M の線平均、
軸と η 0.1 の試験部の平均 100(M/6−1) [%pt]。
式 (§4.18): μ = 180k〜200k の時間窓平均 (台形積分 / 窓長)、d = |平均(160k〜180k) − μ|、U = max_{160k〜200k} |q − μ|、
E = |μ_G − μ_ref|、B = E + U_G + U_ref + s_G + s_ref + p_G + p_ref。相対量は共通の分母 |μ_ref| で % にする。
M の %pt の量はそのまま %pt。判定資格: 全量で check_quasisteady (系列) が STEADY、d ≤ τ/8、U ≤ τ/4。
s が未測定なら合格にしない (保留)。160k〜180k の系列が無い腕 (G1・G1x は 180k〜200k の 9 出力しか残っていない) は暫定の比較に限る。"""
import argparse, json, re, subprocess, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUTD = HERE / "_band_ab" / "core_grid"
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
WINDOWS = ((40.0, 50.0), (65.0, 75.0), (84.0, 94.0))
W1, W2 = (160000, 180000), (180000, 200000)
# 後処理の予算 p (§4.10・§4.18、真の分布 2 通りの最大、Gpost (nj 560) 比)。M の量は全格子で同じ値 (§4.10)
P_BUDGET = {"G0": {"theta": 0.006, "delta": 0.026, "Q_w": 0.009}, "G1": {"theta": 0.147, "delta": 0.225, "Q_w": 0.028},
            "G1x": {"theta": 0.147, "delta": 0.225, "Q_w": 0.06}, "Gref": {"theta": 0.010, "delta": 0.023, "Q_w": 0.009}}
P_M = {"exit_M_line": 0.0006, "eta_pctpt": 0.0014}
# check_quasisteady (系列、--tail 1.0) の drift/osc の閾値。θ・δ・Q_w は §4.15 と同じ、M の量は許容差 0.02 % の 1/10 (eta は M/6 の比の系列で判定)
QS_TOL = {"theta": 1e-4, "delta": 1e-4, "Q_w": 2e-4, "exit_M_line": 2e-5, "eta_pctpt": 2e-5}


def cmd_series(a):
    import h5py
    sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "design"))
    import cold_xcheck as XC
    from forge_design.report import nozzle_report as NR
    yb_x, yb = XC.common_yb(); run = Path(a.run).resolve()
    res = sorted((int(m.group(1)), f.name) for f in run.iterdir() if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)))
    res = [(n, f) for n, f in res if a.lo <= n + a.offset <= a.hi]
    rows = []
    for n, f in res:
        with h5py.File(run / f, "r") as h:
            xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
            v = [np.asarray(h["VALUE/" + k][:], float) for k in ("ro", "Ux", "Uy", "T", "k")]
        bad = int(sum(np.count_nonzero(~np.isfinite(x)) for x in v))
        ni, nj, S = XC.mesh_info(len(v[0]), run)
        o = XC.reduce_fields(xy, *v, False, ni, nj, S, yb_x, yb, profile="pchip")
        rec = {"step": n + a.offset, "nonfinite": bad}
        for lo, hi in WINDOWS:
            w = (o["x"] >= lo) & (o["x"] <= hi)
            for q in ("theta_r", "delta_loc"):
                rec[f"{q}_w{int(lo)}_{int(hi)}"] = float(np.trapezoid(o[q][w], o["x"][w]) / (o["x"][w][-1] - o["x"][w][0]))
        rec["Q_w"] = float(o["Q_w"])
        F = NR.load_field(run, f); info = F["info"]; Md = float(info.get("Md", 6.0)); xE = float(info["x_E"]); xF = float(F["X"][-1, 0])
        eta_last = F["R"][-1] / F["R"][-1, -1]; rec["exit_M_line"] = NR.exit_core_line_mean(eta_last, F["V"]["M"][-1])
        xq = np.linspace(float(F["X"][0, 0]), xF, 2401); w = (xq >= xE + 2) & (xq <= xF - 1)
        for et in (0.0, 0.1):
            rec[f"eta{et}_mean_pct"] = float((100 * (NR.eta_line(F, "M", et, xq) / Md - 1))[w].mean())
        rows.append(rec); print(f"[series {a.name}] {f} → 通算 {rec['step']}", flush=True)
    if not rows:
        raise SystemExit(f"{run}: 通算 {a.lo}〜{a.hi} の出力が無い")
    keys = [k for k in rows[0] if k != "step"]
    p = OUTD / f"gref_series_{a.name}.csv"
    with open(p, "w") as fh:
        fh.write("step," + ",".join(keys) + "\n")
        for r in rows:
            fh.write(f"{r['step']}," + ",".join(f"{r[k]:.12g}" for k in keys) + "\n")
    (OUTD / f"gref_series_{a.name}.json").write_text(json.dumps({"run": str(run), "offset": a.offset, "from": a.lo, "to": a.hi, "n": len(rows)}, ensure_ascii=False))
    print(f"[series] {p} ({len(rows)} 出力)")


def tmean(s, v, lo, hi):
    m = (s >= lo) & (s <= hi)
    if m.sum() < 2 or s[m][0] != lo or s[m][-1] != hi:
        return None
    return float(np.trapezoid(v[m], s[m]) / (hi - lo))


def cmd_judge(a):
    import csv
    tau = {"theta": a.tau_theta, "delta": a.tau_delta, "Q_w": a.tau_qw, "exit_M_line": a.tau_exit, "eta_pctpt": a.tau_eta}
    def kind(k):
        return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"
    S = {}
    for f in sorted(OUTD.glob("gref_series_*.csv")):
        rows = list(csv.DictReader(open(f)))
        S[f.stem.replace("gref_series_", "")] = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    if "Gref" not in S:
        raise SystemExit("Gref の系列が無い")
    keys = [k for k in S["Gref"] if k not in ("step", "nonfinite")]
    out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.18", "tau": tau, "windows": {"mu": W2, "d_ref": W1}, "arms": {}}
    # 判定資格: check_quasisteady の系列モード (CSV だけを渡す。run は読まない)。M の %pt の量は M/6 の比の系列に直して判定する
    for name, D in S.items():
        s = D["step"]; m = (s >= W1[0]) & (s <= W2[1]); tmp = OUTD / f"_qs_{name}.csv"
        cols = {k: D[k][m] for k in keys}
        for k in keys:
            if k.startswith("eta"):
                cols[k] = 1 + cols[k] / 100
        with open(tmp, "w") as fh:
            fh.write("step," + ",".join(keys) + "\n")
            for i in range(int(m.sum())):
                fh.write(f"{s[m][i]:.0f}," + ",".join(f"{cols[k][i]:.12g}" for k in keys) + "\n")
        qs = {}
        for k in keys:
            tol = QS_TOL[kind(k)]
            p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                                "--tail", "1.0", "--drift", f"{tol:g}", "--osc", f"{tol:g}"], capture_output=True, text=True)
            ln = [l for l in p.stdout.splitlines() if l.strip().startswith(k)]
            qs[k] = {"line": ln[-1].strip() if ln else f"(行なし) rc={p.returncode}", "steady": p.returncode == 0, "tol": tol}
        tmp.unlink()
        out["arms"][name] = {"n": int(m.sum()), "steps": [float(s[m][0]), float(s[m][-1])] if m.any() else None,
                             "nonfinite": int(D["nonfinite"][m].sum()), "check_quasisteady": qs}
    ref = S["Gref"]; sr = ref["step"]
    for name, D in S.items():
        s = D["step"]; arm = out["arms"][name]; arm["q"] = {}; full = tmean(s, s, *W1) is not None
        arm["window_160_200"] = full
        for k in keys:
            mu = tmean(s, D[k], *W2); mu_ref = tmean(sr, ref[k], *W2)
            if mu is None or mu_ref is None:
                arm["q"][k] = {"note": "180k〜200k の窓が揃わない"}; continue
            pct = kind(k) == "eta_pctpt"
            scale = (lambda x: x) if pct else (lambda x: 100 * x / abs(mu_ref))
            m = (s >= (W1[0] if full else W2[0])) & (s <= W2[1])
            U = scale(float(np.abs(D[k][m] - mu).max()))
            d = scale(abs(tmean(s, D[k], *W1) - mu)) if full else None
            r = {"mu": mu, "U": U, "d": d, "tau": tau[kind(k)], "steady": arm["check_quasisteady"][k]["steady"]}
            if name != "Gref":
                mr = (sr >= W1[0]) & (sr <= W2[1]); U_ref = scale(float(np.abs(ref[k][mr] - mu_ref).max()))
                pk = P_M[kind(k)] if kind(k) in P_M else P_BUDGET.get(name, {}).get(kind(k))
                pr = P_M[kind(k)] if kind(k) in P_M else P_BUDGET["Gref"][kind(k)]
                E = scale(abs(mu - mu_ref)); r.update(E=E, signed=scale(mu - mu_ref), U_ref=U_ref, p=pk, p_ref=pr,
                                                      B_without_s=(E + U + U_ref + (pk or 0) + pr))
            arm["q"][k] = r
        qual = [arm["nonfinite"] == 0, full] + [v.get("steady", False) and v.get("d") is not None and v["d"] <= v["tau"] / 8 and v["U"] <= v["tau"] / 4
                                               for v in arm["q"].values()]
        arm["qualified"] = all(qual)
    # Gref の窓の A/B (§4.18): 160k〜180k と 180k〜200k の平均の差 d ≤ τ/8 と両窓の資格 (= 上の qualified)
    g = out["arms"]["Gref"]; g["verdict"] = "資格あり (次へ)" if g.get("qualified") else "判定不能 (Gref の資格なし)"
    for name, arm in out["arms"].items():           # 資格を全腕で出してから判定する (Gref の資格を先に要る)
        if name == "Gref":
            continue
        within = all(v.get("B_without_s", np.inf) <= v.get("tau", 0) for v in arm["q"].values())
        arm["verdict"] = ("暫定の比較のみ (160k〜180k の系列が無い)" if not arm["window_160_200"] else
                          "判定不能 (資格なし)" if not (arm["qualified"] and g["qualified"]) else
                          "保留 (s 未測定。s を除いた予算は τ 以内)" if within else "予算超過 (s を除いても τ を超える。真の格子誤差が τ 外とは断定しない)")
    (OUTD / "gref_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for name, arm in out["arms"].items():
        print(f"== {name}: n {arm['n']} 範囲 {arm['steps']} 非有限 {arm['nonfinite']} 資格 {arm.get('qualified')} → {arm.get('verdict')}")
        for k, v in arm["q"].items():
            print(f"   {k:22s} " + " ".join(f"{x} {v[x]:.4g}" if isinstance(v.get(x), float) else f"{x} {v.get(x)}" for x in ("signed", "U", "d", "U_ref", "p", "B_without_s", "tau") if x in v)
                  + f" | {arm['check_quasisteady'][k]['line'][-40:]}")
    print(f"[judge] {OUTD / 'gref_judge.json'}")


# ---- §4.21 の A/B (Gref の継続 A 対 Grefx B): 新しい資格の基準 (codex 諮問 §4.20、ユーザ決定で M の τ 0.05) ----
W1B, W2B = (20000, 40000), (40000, 60000)     # run の step (各 6 万 step、最初の 2 万は外す)
P_GREFX = None                                 # Grefx の後処理の予算 (§4.21、post-ab で測った値を judge2 の引数で渡す)


def cmd_judge2(a):
    """A = gref_series_<a.arm_a>.csv、B = gref_series_<a.arm_b>.csv (run の step)。量ごとに:
    μ = 40k〜60k の時間窓平均、d = |平均(20k〜40k) − μ|、U = max_{20k〜60k} |q − μ|、相対量の分母は D = |μ_A|、M の %pt の量は差。
    資格 (両腕): 非有限 0、check_quasisteady (系列、20k〜60k、--tail 1.0) が STEADY で、drift の閾値 = τ/8・osc (全振幅) の閾値 = τ/2 を
    系列の平均を分母とする比に換算して渡す (M の %pt の量は M/6 の比の系列で τ/100)。d ≤ τ/8、U ≤ τ/4。
    単調な量は check_quasisteady の漸近値の「最終比」が τ/8 を超えたら保留。
    結果: 両腕が資格を満たし全量で E + U_A + U_B + p_A + p_B ≤ τ → A (x の粗さの大きな影響を退ける → s の測定へ)。
    両腕が資格を満たし、どれかの量で E − U_A − U_B − p_A − p_B > τ → B (この粗さの幅での候補化を見送る)。それ以外 → 判定不能。"""
    import csv
    tau = {"theta": a.tau_theta, "delta": a.tau_delta, "Q_w": a.tau_qw, "exit_M_line": a.tau_exit, "eta_pctpt": a.tau_eta}
    pB = json.loads(a.p_b)
    P = {"A": {"theta": P_BUDGET["Gref"]["theta"], "delta": P_BUDGET["Gref"]["delta"], "Q_w": P_BUDGET["Gref"]["Q_w"], **P_M},
         "B": {**{k: pB[k] for k in ("theta", "delta", "Q_w")}, **P_M}}
    def kind(k):
        return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"
    S = {}
    for tag, nm in (("A", a.arm_a), ("B", a.arm_b)):
        rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv")))
        S[tag] = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    keys = [k for k in S["A"] if k not in ("step", "nonfinite")]
    def tm(D, k, lo, hi):
        s = D["step"]; m = (s >= lo) & (s <= hi)
        if m.sum() < 2 or s[m][0] != lo or s[m][-1] != hi:
            raise SystemExit(f"窓 {lo}〜{hi} の出力が揃わない")
        return float(np.trapezoid(D[k][m], s[m]) / (hi - lo))
    out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.21", "tau": tau, "p": P, "arms": {"A": a.arm_a, "B": a.arm_b}, "q": {}, "qual": {}}
    muA = {k: tm(S["A"], k, *W2B) for k in keys}
    for tag in ("A", "B"):
        D = S[tag]; s = D["step"]; m = (s >= W1B[0]) & (s <= W2B[1]); qual = {"nonfinite": int(D["nonfinite"][m].sum())}
        tmp = OUTD / f"_qs2_{tag}.csv"
        with open(tmp, "w") as fh:
            fh.write("step," + ",".join(keys) + "\n")
            for i in np.where(m)[0]:
                fh.write(f"{s[i]:.0f}," + ",".join(f"{(1 + D[k][i] / 100) if k.startswith('eta') else D[k][i]:.12g}" for k in keys) + "\n")
        for k in keys:
            t = tau[kind(k)]; mu = tm(D, k, *W2B)
            if kind(k) == "eta_pctpt":
                ratio_mean = 1 + mu / 100; rel_tau = (t / 100) / abs(ratio_mean)          # %pt の τ を M/6 の比の系列の相対に
                scale = (lambda x: x)
            else:
                rel_tau = (t / 100) * abs(muA[k]) / abs(mu)                                 # 共通の分母 D = |μ_A| の τ を、この系列の平均の相対に
                scale = (lambda x, D_=abs(muA[k]): 100 * x / D_)
            p_ = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                                 "--tail", "1.0", "--drift", f"{rel_tau / 8:.6g}", "--osc", f"{rel_tau / 2:.6g}"], capture_output=True, text=True)
            ln = [l for l in p_.stdout.splitlines() if l.strip().startswith(k)]; line = ln[-1].strip() if ln else f"(行なし) rc={p_.returncode}"
            mono = re.search(r"最終比 ([+\-−]?[\d.]+) ?%", line); asym_ok = True
            if mono:
                y = abs(float(mono.group(1).replace("−", "-")))                           # 漸近値と最終値の差 [最終値の %]
                asym_ok = (y * (1 + mu / 100) if kind(k) == "eta_pctpt" else y * abs(mu) / abs(muA[k])) <= t / 8   # %pt / 共通分母の % に直す
            d = scale(abs(tm(D, k, *W1B) - mu)); U = scale(float(np.abs(D[k][m] - mu).max()))
            qual[k] = {"mu": mu, "d": d, "U": U, "tau": t, "qs": line, "steady": p_.returncode == 0, "asym_ok": asym_ok,
                       "drift_tol": rel_tau / 8, "osc_tol": rel_tau / 2,
                       "ok": p_.returncode == 0 and asym_ok and d <= t / 8 and U <= t / 4}
        tmp.unlink()
        qual["qualified"] = qual["nonfinite"] == 0 and all(v["ok"] for k, v in qual.items() if isinstance(v, dict))
        out["qual"][tag] = qual
    within, beyond = True, False
    for k in keys:
        qa, qb = out["qual"]["A"][k], out["qual"]["B"][k]; t = tau[kind(k)]; kk = kind(k)
        E = (qb["mu"] - qa["mu"]) if kk == "eta_pctpt" else 100 * (qb["mu"] - qa["mu"]) / abs(muA[k])
        Ub = qa["U"] + qb["U"] + P["A"][kk] + P["B"][kk]
        out["q"][k] = {"E_signed": E, "budget_U_p": Ub, "upper": abs(E) + Ub, "lower": abs(E) - Ub, "tau": t}
        within &= abs(E) + Ub <= t; beyond |= abs(E) - Ub > t
    both = out["qual"]["A"]["qualified"] and out["qual"]["B"]["qualified"]
    out["result"] = ("判定不能 (資格なし: " + ", ".join(tag for tag in ("A", "B") if not out["qual"][tag]["qualified"]) + ")") if not both else \
                    "A: x の粗さの大きな影響を退ける (s 抜きの予算は τ 以内) → Grefx と Gref の s の測定へ" if within else \
                    "B: この粗さの幅での候補化を見送る (どれかの量で E − U − p > τ。真の格子誤差が τ を超えたとは断定しない)" if beyond else "判定不能 (境界域)"
    (OUTD / "grefx_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for tag in ("A", "B"):
        q = out["qual"][tag]; print(f"== {tag} ({out['arms'][tag]}): 非有限 {q['nonfinite']} 資格 {q['qualified']}")
        for k in keys:
            v = q[k]; print(f"   {k:18s} d {v['d']:.4f} U {v['U']:.4f} τ {v['tau']:g} STEADY {v['steady']} 漸近 {v['asym_ok']} | {v['qs'][-48:]}")
    for k in keys:
        v = out["q"][k]; print(f"   {k:18s} E {v['E_signed']:+.4f}  U+p {v['budget_U_p']:.4f}  上 {v['upper']:.4f}  下 {v['lower']:.4f}  τ {v['tau']:g}")
    print("結果:", out["result"]); print(f"[judge2] {OUTD / 'grefx_judge.json'}")


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("series"); s.add_argument("name"); s.add_argument("run")
    s.add_argument("--offset", type=int, default=0, help="run の step に足して通算にする (run_0485 は 140000)")
    s.add_argument("--from", dest="lo", type=int, default=W1[0]); s.add_argument("--to", dest="hi", type=int, default=W2[1])
    j = sp.add_parser("judge")
    j.add_argument("--tau-theta", type=float, default=1.0); j.add_argument("--tau-delta", type=float, default=1.0); j.add_argument("--tau-qw", type=float, default=1.0)
    j.add_argument("--tau-exit", type=float, default=0.05); j.add_argument("--tau-eta", type=float, default=0.05)   # M の許容差: 2026-10-11 ユーザ決定で 0.02 → 0.05 (§4.21 から。§4.19 の gref_judge.json は 0.02 で出した)
    j2 = sp.add_parser("judge2", help="§4.21 の A/B (Gref の継続 対 Grefx)")
    j2.add_argument("arm_a"); j2.add_argument("arm_b"); j2.add_argument("--p-b", required=True, help='Grefx の p (JSON: {"theta":..,"delta":..,"Q_w":..} [%%])')
    j2.add_argument("--tau-theta", type=float, default=1.0); j2.add_argument("--tau-delta", type=float, default=1.0); j2.add_argument("--tau-qw", type=float, default=1.0)
    j2.add_argument("--tau-exit", type=float, default=0.05); j2.add_argument("--tau-eta", type=float, default=0.05)
    a = ap.parse_args()
    {"series": cmd_series, "judge": cmd_judge, "judge2": cmd_judge2}[a.cmd](a)


if __name__ == "__main__":
    main()
