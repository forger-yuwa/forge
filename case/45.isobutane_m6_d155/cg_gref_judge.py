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


# §4.25 (ユーザ決定 2026-10-11): 判定の量は試験部の入口側の窓 [40,50] を外した 8 量 ([40,50] は記録だけ)
REQUIRED8 = ("theta_r_w65_75", "theta_r_w84_94", "delta_loc_w65_75", "delta_loc_w84_94", "Q_w", "exit_M_line", "eta0.0_mean_pct", "eta0.1_mean_pct")


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
    WIN = {"A": (W1B, W2B), "B": ((W1B[0] + a.b_shift, W1B[1] + a.b_shift), (W2B[0] + a.b_shift, W2B[1] + a.b_shift))}   # §4.23: B の延長は通算 80k〜120k
    for tag, nm in (("A", a.arm_a), ("B", a.arm_b)):
        rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv")))
        S[tag] = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    keys = [k for k in S["A"] if k not in ("step", "nonfinite")]
    if a.required8:
        miss = [k for k in REQUIRED8 if k not in keys]
        if miss:
            raise SystemExit(f"必須の量が系列に無い (判定不能): {miss}")
        keys = list(REQUIRED8)
    def tm(D, k, lo, hi):
        s = D["step"]; m = (s >= lo) & (s <= hi)
        if m.sum() < 2 or s[m][0] != lo or s[m][-1] != hi:
            raise SystemExit(f"窓 {lo}〜{hi} の出力が揃わない")
        return float(np.trapezoid(D[k][m], s[m]) / (hi - lo))
    out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.21", "tau": tau, "p": P, "arms": {"A": a.arm_a, "B": a.arm_b}, "q": {}, "qual": {}}
    muA = {k: tm(S["A"], k, *W2B) for k in keys}
    for tag in ("A", "B"):
        w1, w2 = WIN[tag]
        D = S[tag]; s = D["step"]; m = (s >= w1[0]) & (s <= w2[1]); qual = {"nonfinite": int(D["nonfinite"][m].sum()), "window": [w1, w2]}
        tmp = OUTD / f"_qs2_{tag}.csv"
        with open(tmp, "w") as fh:
            fh.write("step," + ",".join(keys) + "\n")
            for i in np.where(m)[0]:
                fh.write(f"{s[i]:.0f}," + ",".join(f"{(1 + D[k][i] / 100) if k.startswith('eta') else D[k][i]:.12g}" for k in keys) + "\n")
        for k in keys:
            t = tau[kind(k)]; mu = tm(D, k, *w2); cmean = float(np.mean(D[k][m])); last = float(D[k][m][-1])   # 分類の窓の標本平均・最終値
            if kind(k) == "eta_pctpt":
                rel_tau = (t / 100) / abs(1 + cmean / 100)                                  # %pt の τ を M/6 の比の系列 (分類の窓の平均) の相対に
                scale = (lambda x: x)
            else:
                rel_tau = (t / 100) * abs(muA[k]) / abs(cmean)                              # 共通の分母 |μ_A| の τ を、分類の窓の標本平均の相対に
                scale = (lambda x, D_=abs(muA[k]): 100 * x / D_)
            p_ = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                                 "--tail", "1.0", "--drift", f"{rel_tau / 8:.6g}", "--osc", f"{rel_tau / 2:.6g}"], capture_output=True, text=True)
            ln = [l for l in p_.stdout.splitlines() if l.strip().startswith(k)]; line = ln[-1].strip() if ln else f"(行なし) rc={p_.returncode}"
            mono = re.search(r"最終比 ([+\-−]?[\d.]+) ?%", line); asym_ok = True
            if mono:
                y = abs(float(mono.group(1).replace("−", "-")))                           # 漸近値と最終値の差 [最終値の %]
                asym_ok = (y * (1 + last / 100) if kind(k) == "eta_pctpt" else y * abs(last) / abs(muA[k])) <= t / 8   # 最終値の % を %pt / 共通分母の % に直す
            d = scale(abs(tm(D, k, *w1) - mu)); U = scale(float(np.abs(D[k][m] - mu).max()))
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
    (OUTD / a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for tag in ("A", "B"):
        q = out["qual"][tag]; print(f"== {tag} ({out['arms'][tag]}): 非有限 {q['nonfinite']} 資格 {q['qualified']}")
        for k in keys:
            v = q[k]; print(f"   {k:18s} d {v['d']:.4f} U {v['U']:.4f} τ {v['tau']:g} STEADY {v['steady']} 漸近 {v['asym_ok']} | {v['qs'][-48:]}")
    for k in keys:
        v = out["q"][k]; print(f"   {k:18s} E {v['E_signed']:+.4f}  U+p {v['budget_U_p']:.4f}  上 {v['upper']:.4f}  下 {v['lower']:.4f}  τ {v['tau']:g}")
    print("結果:", out["result"]); print(f"[judge2] {OUTD / a.out}")


def cmd_judge3(a):
    """§4.26 (事前登録、codex 諮問 2026-10-11): 反復条件の感度 s を含めた判定。基準 (B0) の窓は §4.24 のまま固定し、
    point の継続 (各 6 万 step) を両格子で 1 回だけ回す。量は REQUIRED8。相対量の分母は D = |μ_ref,base| (固定)、M の量は %pt の差。
      基準: μ_base・U_base は §4.24 と同じ (Gref = GrefA の 40k〜60k / 20k〜60k、Grefx = Grefx2 の 100k〜120k / 80k〜120k)
      point: μ_point = 40k〜60k の時間窓平均、d = |平均(20k〜40k) − μ_point|、U_point = max_{20k〜60k} |q − μ_point|、
             資格 = 非有限 0・17 出力・check_quasisteady (20k〜60k、--tail 1.0、drift τ/8・osc τ/2 を換算) STEADY・d ≤ τ/8・U ≤ τ/4・漸近値の差 ≤ τ/8
      s = |μ_point − μ_base|、B_total = |E| + U_ref,base + U_x,base + p_ref + p_x + s_ref + s_x + U_ref,point + U_x,point ≤ τ
    結果: 両 point が資格を満たし全 8 量で B_total ≤ τ → A (精度・反復の頑健性の条件で候補としてユーザに諮る)。資格ありでどれかが超過 → B。資格なし → 判定不能。延長しない。"""
    import csv
    tau = {"theta": a.tau_theta, "delta": a.tau_delta, "Q_w": a.tau_qw, "exit_M_line": a.tau_exit, "eta_pctpt": a.tau_eta}
    pX = json.loads(a.p_x)
    P = {"ref": {"theta": P_BUDGET["Gref"]["theta"], "delta": P_BUDGET["Gref"]["delta"], "Q_w": P_BUDGET["Gref"]["Q_w"], **P_M},
         "x": {**{k: pX[k] for k in ("theta", "delta", "Q_w")}, **P_M}}
    def kind(k):
        return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"
    def load(nm):
        rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv")))
        return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    S = {"ref_base": load(a.ref_base), "x_base": load(a.x_base), "ref_point": load(a.ref_point), "x_point": load(a.x_point)}
    for nm, D in S.items():
        miss = [k for k in REQUIRED8 if k not in D]
        if miss:
            raise SystemExit(f"{nm}: 必須の量が系列に無い (判定不能): {miss}")
    sh = a.x_base_shift
    WB = {"ref_base": ((20000, 40000), (40000, 60000)), "x_base": ((20000 + sh, 40000 + sh), (40000 + sh, 60000 + sh)),
          "ref_point": ((20000 + a.ref_point_shift, 40000 + a.ref_point_shift), (40000 + a.ref_point_shift, 60000 + a.ref_point_shift)),
          "x_point": ((20000, 40000), (40000, 60000))}                                                   # §4.28: Gref の point の延長は通算 80k〜120k
    def tm(D, k, lo, hi):
        s_ = D["step"]; m = (s_ >= lo) & (s_ <= hi)
        if m.sum() < 2 or s_[m][0] != lo or s_[m][-1] != hi:
            raise SystemExit(f"窓 {lo}〜{hi} の出力が揃わない (判定不能)")
        return float(np.trapezoid(D[k][m], s_[m]) / (hi - lo))
    D0 = {k: abs(tm(S["ref_base"], k, *WB["ref_base"][1])) for k in REQUIRED8}
    def sc(k, x):
        return x if kind(k) == "eta_pctpt" else 100 * x / D0[k]
    st = {}
    for nm, D in S.items():
        w1, w2 = WB[nm]; s_ = D["step"]; m = (s_ >= w1[0]) & (s_ <= w2[1]); r = {"n": int(m.sum()), "nonfinite": int(D["nonfinite"][m].sum()), "q": {}}
        if nm.endswith("point"):
            tmp = OUTD / f"_qs3_{nm}.csv"
            with open(tmp, "w") as fh:
                fh.write("step," + ",".join(REQUIRED8) + "\n")
                for i in np.where(m)[0]:
                    fh.write(f"{s_[i]:.0f}," + ",".join(f"{(1 + D[k][i] / 100) if k.startswith('eta') else D[k][i]:.12g}" for k in REQUIRED8) + "\n")
        for k in REQUIRED8:
            t = tau[kind(k)]; mu = tm(D, k, *w2); U = sc(k, float(np.abs(D[k][m] - mu).max()))
            q = {"mu": mu, "U": U}
            if nm.endswith("point"):
                cmean = float(np.mean(D[k][m])); last = float(D[k][m][-1])
                rel = (t / 100) / abs(1 + cmean / 100) if kind(k) == "eta_pctpt" else (t / 100) * D0[k] / abs(cmean)
                pr = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                                     "--tail", "1.0", "--drift", f"{rel / 8:.6g}", "--osc", f"{rel / 2:.6g}"], capture_output=True, text=True)
                ln = [l for l in pr.stdout.splitlines() if l.strip().startswith(k)]; line = ln[-1].strip() if ln else f"(行なし) rc={pr.returncode}"
                mono = re.search(r"最終比 ([+\-−]?[\d.]+) ?%", line); asym_ok = True
                if mono:
                    y = abs(float(mono.group(1).replace("−", "-")))
                    asym_ok = (y * (1 + last / 100) if kind(k) == "eta_pctpt" else y * abs(last) / D0[k]) <= t / 8
                d = sc(k, abs(tm(D, k, *w1) - mu))
                q.update(d=d, qs=line, steady=pr.returncode == 0, asym_ok=asym_ok, ok=pr.returncode == 0 and asym_ok and d <= t / 8 and U <= t / 4)
            r["q"][k] = q
        if nm.endswith("point"):
            tmp.unlink()
            r["qualified"] = r["nonfinite"] == 0 and r["n"] == 17 and all(v["ok"] for v in r["q"].values())
        st[nm] = r
    out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.26", "tau": tau, "p": P, "series": {k: getattr(a, k) for k in ("ref_base", "x_base", "ref_point", "x_point")},
           "stats": st, "q": {}}
    within, beyond = True, False
    for k in REQUIRED8:
        t = tau[kind(k)]; kk = kind(k)
        E = sc(k, st["x_base"]["q"][k]["mu"] - st["ref_base"]["q"][k]["mu"])
        s_ref = sc(k, abs(st["ref_point"]["q"][k]["mu"] - st["ref_base"]["q"][k]["mu"])); s_x = sc(k, abs(st["x_point"]["q"][k]["mu"] - st["x_base"]["q"][k]["mu"]))
        B = (abs(E) + st["ref_base"]["q"][k]["U"] + st["x_base"]["q"][k]["U"] + P["ref"][kk] + P["x"][kk]
             + s_ref + s_x + st["ref_point"]["q"][k]["U"] + st["x_point"]["q"][k]["U"])
        out["q"][k] = {"E_signed": E, "s_ref": s_ref, "s_x": s_x, "B_total": B, "tau": t}
        within &= B <= t
    both = st["ref_point"]["qualified"] and st["x_point"]["qualified"]
    out["result"] = ("判定不能 (point の資格なし: " + ", ".join(n for n in ("ref_point", "x_point") if not st[n]["qualified"]) + ")") if not both else \
                    "A: 精度・反復の頑健性の条件で候補 (全 8 量で B_total ≤ τ) → ユーザに諮る (総時間の条件は別)" if within else \
                    "B: この予算では候補化を見送る (どれかの量で B_total > τ。真の格子誤差が τ を超えたとは断定しない)"
    (OUTD / a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for nm in ("ref_point", "x_point"):
        r = st[nm]; print(f"== {nm} ({getattr(a, nm)}): n {r['n']} 非有限 {r['nonfinite']} 資格 {r['qualified']}")
        for k in REQUIRED8:
            v = r["q"][k]; print(f"   {k:18s} d {v['d']:.4f} U {v['U']:.4f} STEADY {v['steady']} 漸近 {v['asym_ok']} | {v['qs'][-44:]}")
    for k in REQUIRED8:
        v = out["q"][k]; print(f"   {k:18s} E {v['E_signed']:+.4f}  s_ref {v['s_ref']:.4f}  s_x {v['s_x']:.4f}  B_total {v['B_total']:.4f}  τ {v['tau']:g}")
    print("結果:", out["result"]); print(f"[judge3] {OUTD / a.out}")


def cmd_judge4(a):
    """§4.30 (事前登録、codex 諮問 2026-10-11): 物差しを「Gref の point の資格を満たした値」に替えた新しい比較。旧い試験 (§4.26・§4.28) の判定不能は残す。
      参照: Gref の point (通算 160k〜200k、μ_ref,point = 180k〜200k、d = |平均(160k〜180k) − μ|、U_ref,point = 160k〜200k)。資格は §4.26 と同じ (8 量)。
      候補: Grefx の B0 (μ_x,base = 通算 100k〜120k、U_x,base = 80k〜120k) と point (μ = 40k〜60k、U_x,point = 20k〜60k)、s_x = |μ_x,point − μ_x,base|。
      E′ = μ_x,base − μ_ref,point、B′ = |E′| + U_x,base + U_ref,point + p_x + p_ref + s_x + U_x,point ≤ τ。相対量の分母は D = |μ_ref,base| (Gref の B0 の 40k〜60k、従来と同じ)。
      s_ref (Gref の B0 ↔ point) と U_ref,base は予算に入れず記録する。
    結果: 参照の資格あり (かつ Grefx の point の資格あり) で全 8 量 B′ ≤ τ → 指定の参照・許容差に対する精度面の候補 (総時間の条件は別)。資格ありで超過 → 見送り。資格なし → 判定不能 (延長しない)。"""
    import csv
    tau = {"theta": a.tau_theta, "delta": a.tau_delta, "Q_w": a.tau_qw, "exit_M_line": a.tau_exit, "eta_pctpt": a.tau_eta}
    pX = json.loads(a.p_x)
    P = {"ref": {"theta": P_BUDGET["Gref"]["theta"], "delta": P_BUDGET["Gref"]["delta"], "Q_w": P_BUDGET["Gref"]["Q_w"], **P_M},
         "x": {**{k: pX[k] for k in ("theta", "delta", "Q_w")}, **P_M}}
    def kind(k):
        return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"
    def load(nm):
        rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv")))
        return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    S = {"ref_base": load(a.ref_base), "ref_point": load(a.ref_point), "x_base": load(a.x_base), "x_point": load(a.x_point)}
    for nm, D in S.items():
        miss = [k for k in REQUIRED8 if k not in D]
        if miss:
            raise SystemExit(f"{nm}: 必須の量が系列に無い (判定不能): {miss}")
    WB = {"ref_base": ((20000, 40000), (40000, 60000)), "ref_point": ((160000, 180000), (180000, 200000)),
          "x_base": ((80000, 100000), (100000, 120000)), "x_point": ((20000, 40000), (40000, 60000))}
    def tm(D, k, lo, hi):
        s_ = D["step"]; m = (s_ >= lo) & (s_ <= hi)
        if m.sum() < 2 or s_[m][0] != lo or s_[m][-1] != hi:
            raise SystemExit(f"窓 {lo}〜{hi} の出力が揃わない (判定不能)")
        return float(np.trapezoid(D[k][m], s_[m]) / (hi - lo))
    D0 = {k: abs(tm(S["ref_base"], k, *WB["ref_base"][1])) for k in REQUIRED8}
    def sc(k, x):
        return x if kind(k) == "eta_pctpt" else 100 * x / D0[k]
    st = {}
    for nm, D in S.items():
        w1, w2 = WB[nm]; s_ = D["step"]; m = (s_ >= w1[0]) & (s_ <= w2[1]); r = {"n": int(m.sum()), "nonfinite": int(D["nonfinite"][m].sum()), "q": {}}
        qual_needed = nm in ("ref_point", "x_point")
        if qual_needed:
            tmp = OUTD / f"_qs4_{nm}.csv"
            with open(tmp, "w") as fh:
                fh.write("step," + ",".join(REQUIRED8) + "\n")
                for i in np.where(m)[0]:
                    fh.write(f"{s_[i]:.0f}," + ",".join(f"{(1 + D[k][i] / 100) if k.startswith('eta') else D[k][i]:.12g}" for k in REQUIRED8) + "\n")
        for k in REQUIRED8:
            t = tau[kind(k)]; mu = tm(D, k, *w2); U = sc(k, float(np.abs(D[k][m] - mu).max())); q = {"mu": mu, "U": U}
            if qual_needed:
                cmean = float(np.mean(D[k][m])); last = float(D[k][m][-1])
                rel = (t / 100) / abs(1 + cmean / 100) if kind(k) == "eta_pctpt" else (t / 100) * D0[k] / abs(cmean)
                pr = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                                     "--tail", "1.0", "--drift", f"{rel / 8:.6g}", "--osc", f"{rel / 2:.6g}"], capture_output=True, text=True)
                ln = [l for l in pr.stdout.splitlines() if l.strip().startswith(k)]; line = ln[-1].strip() if ln else f"(行なし) rc={pr.returncode}"
                mono = re.search(r"最終比 ([+\-−]?[\d.]+) ?%", line); asym_ok = True
                if mono:
                    y = abs(float(mono.group(1).replace("−", "-")))
                    asym_ok = (y * (1 + last / 100) if kind(k) == "eta_pctpt" else y * abs(last) / D0[k]) <= t / 8
                d = sc(k, abs(tm(D, k, *w1) - mu))
                q.update(d=d, qs=line, steady=pr.returncode == 0, asym_ok=asym_ok, ok=pr.returncode == 0 and asym_ok and d <= t / 8 and U <= t / 4)
            r["q"][k] = q
        if qual_needed:
            tmp.unlink()
            r["qualified"] = r["nonfinite"] == 0 and r["n"] == 17 and all(v["ok"] for v in r["q"].values())
        st[nm] = r
    out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.30", "tau": tau, "p": P, "windows": WB,
           "series": {k: getattr(a, k) for k in ("ref_base", "ref_point", "x_base", "x_point")}, "stats": st, "q": {}}
    within = True
    for k in REQUIRED8:
        t = tau[kind(k)]; kk = kind(k)
        E = sc(k, st["x_base"]["q"][k]["mu"] - st["ref_point"]["q"][k]["mu"])
        s_x = sc(k, abs(st["x_point"]["q"][k]["mu"] - st["x_base"]["q"][k]["mu"]))
        B = abs(E) + st["x_base"]["q"][k]["U"] + st["ref_point"]["q"][k]["U"] + P["x"][kk] + P["ref"][kk] + s_x + st["x_point"]["q"][k]["U"]
        rec_s_ref = sc(k, abs(st["ref_point"]["q"][k]["mu"] - st["ref_base"]["q"][k]["mu"]))
        out["q"][k] = {"E_signed": E, "s_x": s_x, "B_prime": B, "tau": t, "record_s_ref": rec_s_ref, "record_U_ref_base": st["ref_base"]["q"][k]["U"]}
        within &= B <= t
    both = st["ref_point"]["qualified"] and st["x_point"]["qualified"]
    out["result"] = ("判定不能 (資格なし: " + ", ".join(n for n in ("ref_point", "x_point") if not st[n]["qualified"]) + ")") if not both else \
                    "A: 指定の参照・許容差に対する精度面の候補 (全 8 量で B′ ≤ τ。総時間の条件は別)" if within else \
                    "B: この予算では候補化を見送る (どれかの量で B′ > τ)"
    (OUTD / a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1))
    for nm in ("ref_point", "x_point"):
        r = st[nm]; print(f"== {nm} ({getattr(a, nm)}): n {r['n']} 非有限 {r['nonfinite']} 資格 {r['qualified']}")
        for k in REQUIRED8:
            v = r["q"][k]; print(f"   {k:18s} d {v['d']:.4f} U {v['U']:.4f} STEADY {v['steady']} 漸近 {v['asym_ok']} | {v['qs'][-44:]}")
    for k in REQUIRED8:
        v = out["q"][k]; print(f"   {k:18s} E′ {v['E_signed']:+.4f}  s_x {v['s_x']:.4f}  B′ {v['B_prime']:.4f}  τ {v['tau']:g}  (記録 s_ref {v['record_s_ref']:.4f})")
    print("結果:", out["result"]); print(f"[judge4] {OUTD / a.out}")


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
    j2.add_argument("--b-shift", type=int, default=0, help="B の窓を通算でずらす step (§4.23 の延長は 60000、系列は --offset 60000 で作る)")
    j2.add_argument("--out", default="grefx_judge.json", help="出力 (§4.22 の判定は grefx_judge.json、§4.23 は grefx_ext_judge.json)")
    j2.add_argument("--required8", action="store_true", help="判定の量を §4.25 の 8 量に限る ([40,50] の窓は記録だけ。欠けていれば判定不能)")
    j3 = sp.add_parser("judge3", help="§4.26 の s を含めた判定 (Gref・Grefx の point の継続)")
    j3.add_argument("ref_base"); j3.add_argument("x_base"); j3.add_argument("ref_point"); j3.add_argument("x_point")
    j3.add_argument("--x-base-shift", type=int, default=60000, help="Grefx の基準の系列 (Grefx2) は通算 step なので窓をずらす")
    j3.add_argument("--p-x", required=True, help='Grefx の p (JSON: {"theta":..,"delta":..,"Q_w":..} [%%])')
    j3.add_argument("--tau-theta", type=float, default=1.5); j3.add_argument("--tau-delta", type=float, default=1.5); j3.add_argument("--tau-qw", type=float, default=1.5)
    j3.add_argument("--tau-exit", type=float, default=0.05); j3.add_argument("--tau-eta", type=float, default=0.05)
    j3.add_argument("--out", default="grefx_s_judge.json")
    j3.add_argument("--ref-point-shift", type=int, default=0, help="Gref の point の系列の窓をずらす (§4.28 の延長は 60000、系列は --offset 60000 で作る)")
    j4 = sp.add_parser("judge4", help="§4.30 の比較 (物差し = Gref の point の資格を満たした値)")
    j4.add_argument("ref_base"); j4.add_argument("ref_point"); j4.add_argument("x_base"); j4.add_argument("x_point")
    j4.add_argument("--p-x", required=True, help='Grefx の p (JSON: {"theta":..,"delta":..,"Q_w":..} [%%])')
    j4.add_argument("--tau-theta", type=float, default=1.5); j4.add_argument("--tau-delta", type=float, default=1.5); j4.add_argument("--tau-qw", type=float, default=1.5)
    j4.add_argument("--tau-exit", type=float, default=0.05); j4.add_argument("--tau-eta", type=float, default=0.05)
    j4.add_argument("--out", default="grefx_refpoint_judge.json")
    a = ap.parse_args()
    {"series": cmd_series, "judge": cmd_judge, "judge2": cmd_judge2, "judge3": cmd_judge3, "judge4": cmd_judge4}[a.cmd](a)


if __name__ == "__main__":
    main()
