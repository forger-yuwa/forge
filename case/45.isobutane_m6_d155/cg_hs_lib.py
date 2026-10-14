"""plan tooling-nozzle-core-grid §4.104 P-hs の共通の計算 (見張り `cg_hs_check.py` と判定 `cg_hs_judge.py` が使う)。
- 系列: 出力ごとの `gref_series_<run>_<k>.csv` (8 量、cg_gref_judge.py series) と `stoprows_<run>.jsonl` (停止規則の量) をつなぐ。step は継続の中の番号。
- 停止規則: `stop_rule.judge_rows` (719dfbf4) をそのまま当てる (r の窓平均・振幅、x/x_E の θ_r の D・振幅、窓 2 万 step・9 点以上)。
- 8 量の資格 (§4.30 と同じ式): 窓 [s − 40k, s] の 17 出力、μ = (s − 20k, s] の台形平均、d = |[s − 40k, s − 20k] の台形平均 − μ|、U = 窓の max|q − μ|、
  check_quasisteady (--tail 1.0、drift τ/8・osc τ/2 を系列の値に換算) STEADY、単調なら漸近値と最終値の差 ≤ τ/8。d ≤ τ/8・U ≤ τ/4。
  相対量の分母は呼び出し側が渡す D0 (既定はその腕の μ)。τ は θ_r・δ_loc・Q_w 1.5 %、出口の主流 M 0.05 %、軸・η 0.1 の M 0.05 %pt。"""
import csv, json, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
TOOLS = HERE.parents[1] / "solver_density_cuda" / "tools"
sys.path.insert(0, str(Path.home() / "forge-coregrid" / "_ext" / "f5design_719dfbf4" / "design"))
REQUIRED8 = ("theta_r_w65_75", "theta_r_w84_94", "delta_loc_w65_75", "delta_loc_w84_94", "Q_w", "exit_M_line", "eta0.0_mean_pct", "eta0.1_mean_pct")
TAU = {"theta": 1.5, "delta": 1.5, "Q_w": 1.5, "exit_M_line": 0.05, "eta_pctpt": 0.05}
CONF = 40000
_TRAP = getattr(np, "trapezoid", None) or np.trapz


def kind(k):
    return "theta" if k.startswith("theta_r") else "delta" if k.startswith("delta_loc") else "Q_w" if k == "Q_w" else "exit_M_line" if k == "exit_M_line" else "eta_pctpt"


def series(run, out=2500, upto=None):
    rows = []
    k = out
    while True:
        p = OUTD / f"gref_series_{run}_{k}.csv"
        if not p.exists() or (upto and k > upto):
            break
        r = list(csv.DictReader(open(p)))
        if not r:
            break
        rows.append({"step": k, **{q: float(v) for q, v in r[-1].items() if q != "step"}})
        k += out
    return rows


def stoprows(run):
    p = OUTD / f"stoprows_{run}.jsonl"
    if not p.exists():
        return []
    return sorted({r["step"]: r for r in (json.loads(l) for l in open(p))}.values(), key=lambda r: r["step"])


def stop_judgements(run):
    from forge_design.evaluate import stop_rule as SR
    rows = stoprows(run)
    return rows, [SR.judge_rows(rows, i) for i in range(len(rows))]


def tm(S, V, lo, hi):
    m = (S >= lo) & (S <= hi)
    if m.sum() < 2 or S[m][0] != lo or S[m][-1] != hi:
        return None
    return float(_TRAP(V[m], S[m]) / (hi - lo))


def qualify(rows, s_end, D0=None):
    """窓 [s_end − 40k, s_end] の 8 量の資格。戻り {"ok", "n", "nonfinite", "q": {k: {...}}}。"""
    S = np.array([r["step"] for r in rows], float)
    m = (S >= s_end - CONF) & (S <= s_end)
    out = {"s_end": s_end, "n": int(m.sum()), "q": {}}
    if out["n"] != 17:
        out.update(ok=False, why=f"窓の出力 {out['n']} (17 であること)")
        return out
    out["nonfinite"] = int(sum(r.get("nonfinite", 0) for r, mm in zip(rows, m) if mm))
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, dir=OUTD) as fh:
        tmp = Path(fh.name)
        fh.write("step," + ",".join(REQUIRED8) + "\n")
        for r, mm in zip(rows, m):
            if mm:
                fh.write(f"{r['step']:.0f}," + ",".join(f"{(1 + r[k] / 100) if k.startswith('eta') else r[k]:.12g}" for k in REQUIRED8) + "\n")
    allok = out["nonfinite"] == 0
    for k in REQUIRED8:
        V = np.array([r[k] for r in rows], float)
        t = TAU[kind(k)]
        mu = tm(S, V, s_end - CONF / 2, s_end); w1 = tm(S, V, s_end - CONF, s_end - CONF / 2)
        den = abs(D0[k]) if D0 else abs(mu)
        sc = (lambda x: x) if kind(k) == "eta_pctpt" else (lambda x: 100 * x / den)
        U = sc(float(np.abs(V[m] - mu).max())); d = sc(abs(w1 - mu))
        cmean = float(V[m].mean()); last = float(V[m][-1])
        rel = (t / 100) / abs(1 + cmean / 100) if kind(k) == "eta_pctpt" else (t / 100) * den / abs(cmean)
        pr = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(tmp), "--series-cols", k,
                             "--tail", "1.0", "--drift", f"{rel / 8:.6g}", "--osc", f"{rel / 2:.6g}"], capture_output=True, text=True)
        ln = [l for l in pr.stdout.splitlines() if l.strip().startswith(k)]; line = ln[-1].strip() if ln else f"(行なし) rc={pr.returncode}"
        mono = re.search(r"最終比 ([+\-−]?[\d.]+) ?%", line); asym_ok = True
        if mono:
            y = abs(float(mono.group(1).replace("−", "-")))
            asym_ok = (y * (1 + last / 100) if kind(k) == "eta_pctpt" else y * abs(last) / den) <= t / 8
        ok = pr.returncode == 0 and asym_ok and d <= t / 8 and U <= t / 4
        out["q"][k] = {"mu": mu, "d": d, "U": U, "tau": t, "steady": pr.returncode == 0, "asym_ok": asym_ok, "qs": line, "ok": ok}
        allok = allok and ok
    tmp.unlink()
    out["ok"] = bool(allok)
    return out
