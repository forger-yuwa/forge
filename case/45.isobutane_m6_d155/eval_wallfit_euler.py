"""補間壁 (A) vs 同時当てはめ壁 (B) の Euler A/B を固定評価器で判定する。plan verification-m6-axis-wave-mesh-su2 §5.1 #15 (事前登録)。
A: run_0053〜0055 / B: run_0056〜0058 (各腕 3 回の独立再実行)。
延長 (準定常未達時の 6000 step × 1 回, run_0059〜0064) があれば、本段 (step 1000〜12000) に延長 (12000+1000〜) を連結して判定する。
v2 (2026-10-05, codex 6 回目): 分布差の U に 2T・2E を入れ、E ≤ Δ/10 を前提ゲートに、残差 VERDICT は全文の行を保存。登録時の出力 (wallfit_euler_ab.json) は上書きしない。
usage: python3 eval_wallfit_euler.py [case_dir] [--fixed-coef] [--pair=A,B (既定 interp,fit; 例 fit,v4)]   (→ _band_ab/wallfit_euler_ab_v2.json / _band_ab/wallfit_euler_ab_diag_fixedcoef.json) → _band_ab/wallfit_euler_ab.json, 各 run の wallfit_series.csv と QUASISTEADY_wallfit.txt
--e3 (plan tooling-nozzle-throat-monotone-r2 §6 E1〜E4、§5.1 #5; 2026-10-06): 単調壁の A/B 用の判定。既存の呼び方 (--e3 なし) の挙動は変えない。
  X_E・X_F を各腕の prepare_info.json (x_E) と wall_design.csv (最後の x) から読む (腕・run で一致しなければ止める) /
  準定常 (check_quasisteady --tail) と T を同じ末尾 5 枚に / 出口評価誤差 E_exit (出口帯の端 ±0.05・定義 = 共通標本の節点平均) を
  |出口 M − 6| と出口規格化オーバーシュートの U に入れる / 判定は throat_mono_judge (片側の非劣化・不採用・保留、η0 除外、未定常・欠損) /
  壁の証拠 (M4、throat_mono_ab.wall_evidence) / 腕の run 数が足りない・snapshot が足りないときは欠損として保留。
  例: python3 eval_wallfit_euler.py . --fixed-coef --e3 --pair=pinG1,monoG1  → _band_ab/throat_mono_euler_ab_pinG1_vs_monoG1_fixedcoef.json
  (各 run には wallfit_series_e3.csv・QUASISTEADY_wallfit_e3.txt を書く。既存の評価の出力は上書きしない。--pair は必須)
"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
from scipy.interpolate import BSpline
E3 = "--e3" in sys.argv   # plan tooling-nozzle-throat-monotone-r2 §6 E3 の判定 (下の run_e3)
FIXED = "--fixed-coef" in sys.argv   # codex 6 回目の診断 B: 刻み感度 E を、0.05 で当てた平滑化曲線を固定して 0.025 で評価する
_args = [a for a in sys.argv[1:] if not a.startswith("--")]
C = Path(_args[0]).resolve() if _args else Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.report.nozzle_report import load_field, eta_line  # noqa: E402

X_E, X_F = 39.82004263, 95.22667765
WIN_T = (41.82004263, 94.22667765); WIN_O = (24.82004263, 95.22667765)
P_REF, MD = 2237.0, 6.0
DELTA = dict(M_wave=0.001, M_resid_diff=0.001, overshoot=0.003, exit_core_M=0.00018, P_wave=0.010, P_resid_diff=0.010,
             P_slope_abs=0.03, overshoot_exitnorm=0.003, exit_M_dev=0.00018)   # 2026-10-05 codex plan M1 (CFD ピン): |傾き|・出口 M 規格化オーバーシュート・|出口 M − 6|
_pair = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--pair=")), "interp,fit").split(",")   # A 腕, B 腕 (interp / fit / v4)


def _arm_runs(arm):
    rs = sorted(q.name for q in C.glob(f"run_[0-9][0-9][0-9][0-9]_euler_wallfit_{arm}_r[0-9]") if q.is_dir())   # run 番号 4 桁 (旧: 00xx だけ)
    return rs


ARMS = {"A": _arm_runs(_pair[0]), "B": _arm_runs(_pair[1])}


def grid(h):
    g = np.arange(0.0, X_F + 1e-9, h)
    return np.unique(np.r_[g, WIN_T, WIN_O, X_E, X_F])


def pspline_fixed(x, v, lam=1.0, k=3, coef=False):
    """3 次 B-spline、内部ノット x = 10k (固定)、2 階差分罰則。フィット区間 [X_E, X_F] 固定。coef=True なら BSpline も返す。"""
    m = (x >= X_E) & (x <= X_F)
    ki = np.arange(10.0, X_F, 10.0); ki = ki[(ki > X_E) & (ki < X_F)]
    t = np.r_[[X_E] * (k + 1), ki, [X_F] * (k + 1)]; nc = len(t) - k - 1
    B = BSpline.design_matrix(x[m], t, k).toarray()
    D = np.diff(np.eye(nc), 2, axis=0)
    c = np.linalg.solve(B.T @ B + lam * D.T @ D, B.T @ v[m])
    out = np.full_like(v, np.nan); out[m] = B @ c
    return (out, BSpline(t, c, k)) if coef else out


def quantities(F, h, fixed_from=None):
    """fixed_from=h0 なら、平滑化曲線を h0 格子で当てはめて固定し、h 格子ではその曲線を評価するだけ (再フィットしない; codex 6 回目の診断 B)。"""
    xq = grid(h); wt = (xq >= WIN_T[0]) & (xq <= WIN_T[1]); wo = (xq >= WIN_O[0]) & (xq <= WIN_O[1])
    q, dist = {}, {}
    for eta in (0.0, 0.1):
        dm = 100 * (eta_line(F, "M", eta, xq) / MD - 1); dp = 100 * (eta_line(F, "P", eta, xq) / P_REF - 1)
        if fixed_from is None:
            rm = dm - pspline_fixed(xq, dm); rp = dp - pspline_fixed(xq, dp)
        else:
            x0 = grid(fixed_from)
            _, sm = pspline_fixed(x0, 100 * (eta_line(F, "M", eta, x0) / MD - 1), coef=True)
            _, sp = pspline_fixed(x0, 100 * (eta_line(F, "P", eta, x0) / P_REF - 1), coef=True)
            rm = dm - sm(xq); rp = dp - sp(xq)
        q[f"M_wave_eta{eta}"] = float(np.abs(rm[wt]).max()); q[f"P_wave_eta{eta}"] = float(np.abs(rp[wt]).max())
        q[f"overshoot_eta{eta}"] = float(dm[wo].max())
        cf = np.polyfit(xq[wt], dp[wt], 1); q[f"P_slope_eta{eta}"] = float(cf[0] * (WIN_T[1] - WIN_T[0]))
        dist[f"M_resid_eta{eta}"] = (xq[wt], rm[wt]); dist[f"P_resid_eta{eta}"] = (xq[wt], rp[wt]); dist[f"P_raw_eta{eta}"] = (xq[wt], dp[wt])
    eta_last = F["R"][-1] / F["R"][-1, -1]; eg = np.linspace(0.05, 0.7, 131)
    Me = np.interp(eg, eta_last, F["V"]["M"][-1])
    q["exit_core_M"] = float(np.trapezoid(Me * eg, eg) / np.trapezoid(eg, eg))
    q["exit_M_dev"] = abs(q["exit_core_M"] - MD)
    for eta in (0.0, 0.1):
        q[f"P_slope_abs_eta{eta}"] = abs(q[f"P_slope_eta{eta}"])
        dmx = 100 * (eta_line(F, "M", eta, xq) / q["exit_core_M"] - 1)          # 出口コア M で規格化 (m* によるレベルシフトを除く)
        q[f"overshoot_exitnorm_eta{eta}"] = float(dmx[wo].max())
    return q, dist


def snaps(run):
    fs = [f.name for f in (C / run).glob("res_[0-9]*.h5")]
    return sorted(fs, key=lambda f: int(re.findall(r"\d+", f)[0]))


def ext_of(run):
    m = sorted(C.glob("run_[0-9][0-9][0-9][0-9]_" + run[9:] + "_ext6k"))
    return m[0] if m else None


def series_files(run):
    """[(run_dir, res 名, 通算 step)]。延長があれば本段 12000 の後に連結 (延長の step 0 は本段 12000 と同じ場なので除く)。"""
    rd = C / run; out = [(rd, f, int(re.findall(r"\d+", f)[0])) for f in snaps(run) if int(re.findall(r"\d+", f)[0]) > 0]
    e = ext_of(run)
    if e is not None:
        out += [(e, f, 12000 + int(re.findall(r"\d+", f)[0])) for f in snaps(e.name) if int(re.findall(r"\d+", f)[0]) > 0]
    return out


# --- --e3: plan tooling-nozzle-throat-monotone-r2 §6 E1〜E4 (単調壁の A/B) ---------------------------------------
E3_TAIL = 5                                          # 準定常と T を同じ末尾 5 枚で
EXIT_BAND = (0.05, 0.7)                              # 出口コア M の帯 (η 重み付き面積平均、既定)
EXIT_VARIANTS = {"band_lo-0.05": ("area", 0.0, 0.7), "band_lo+0.05": ("area", 0.1, 0.7),
                 "band_hi-0.05": ("area", 0.05, 0.65), "band_hi+0.05": ("area", 0.05, 0.75),
                 "M_common": ("common", 0.05, 0.7)}  # 出口の評価誤差 (E2: 帯の端 ±0.05、定義 = 共通標本の節点平均 [euler_grid_ab.py の M_common])


def e3_geometry(runs):
    """各 run の prepare_info.json の x_E と wall_design.csv の最後の x (/ scale_m) → (X_E, X_F)。run 間で 1e-9 を超えて違えば止める。"""
    g = []
    for run in runs:
        info = json.loads((C / run / "prepare_info.json").read_text())
        wd = np.loadtxt(C / run / "wall_design.csv", delimiter=",", skiprows=1)
        g.append((float(info["x_E"]), float(wd[-1, 0]) / float(info["scale_m"])))
    if not g:
        raise SystemExit("--e3: 腕の run が 1 本も無い (X_E・X_F を決められない)")
    g = np.array(g)
    if np.ptp(g[:, 0]) > 1e-9 or np.ptp(g[:, 1]) > 1e-9:
        raise SystemExit(f"--e3: run ごとの X_E / X_F が一致しない (A/B の MOC 点群が違う): {g.tolist()}")
    return float(g[0, 0]), float(g[0, 1])


def exit_core(F, kind, lo, hi, etaS=None):
    """出口断面 (最後の i 列) の M。area: η 重み付き面積平均 (刻み 0.005; 既定帯 [0.05, 0.7] は quantities と同一) / common: 共通標本 etaS の節点平均。"""
    eta_last = F["R"][-1] / F["R"][-1, -1]
    if kind == "common":
        return float(np.interp(etaS, eta_last, F["V"]["M"][-1]).mean())
    eg = np.linspace(0.05, 0.7, 131) if (lo, hi) == EXIT_BAND else np.linspace(lo, hi, int(round((hi - lo) / 0.005)) + 1)
    Me = np.interp(eg, eta_last, F["V"]["M"][-1])
    return float(np.trapezoid(Me * eg, eg) / np.trapezoid(eg, eg))


def exit_variant_quantities(F, q, etaS):
    """出口評価の各変種での |出口 M − 6| と出口規格化オーバーシュート (η 0, 0.1)。base は既定の帯。
    規格化オーバーシュートは max_x 100 (M_η/M_exit − 1) = 100 (max M_η / M_exit − 1) なので、quantities のオーバーシュートから戻す。"""
    out = {}
    for name, (kind, lo, hi) in [("base", ("area",) + EXIT_BAND)] + list(EXIT_VARIANTS.items()):
        Mx = exit_core(F, kind, lo, hi, etaS)
        out[name] = {"exit_core_M": Mx, "exit_M_dev": abs(Mx - MD),
                     **{f"overshoot_exitnorm_eta{e}": 100.0 * (MD * (1.0 + q[f"overshoot_eta{e}"] / 100.0) / Mx - 1.0) for e in (0.0, 0.1)}}
    return out


def _verdict_lines(qs, cols):
    v = {}
    for c in cols:
        m = re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", qs, re.M)
        v[c] = m.group(1) if m else "UNKNOWN"
    return v


def e3_run(run, arm, etaS, other_run):
    """1 run の評価: 時系列 CSV・末尾 5 枚の準定常・代表値・T・E・出口変種・分布・壁の証拠。足りなければ status=missing。"""
    from throat_mono_ab import wall_evidence
    rd = C / run
    sf = series_files(run)
    rec = {"arm": arm, "run": run}
    if len(sf) < E3_TAIL:
        rec.update(status="missing", reason=f"snapshot {len(sf)} 枚 (< {E3_TAIL})")
        return rec
    series, exv, dists = [], [], []
    for k, (d_, f, st_) in enumerate(sf):
        F = load_field(d_, f)
        q, dist = quantities(F, 0.05)
        q["step"] = st_
        series.append(q)
        if k >= len(sf) - E3_TAIL:
            exv.append(exit_variant_quantities(F, q, etaS))
            dists.append(dist)
    cols = [k for k in series[0] if k != "step"]
    # 既存の評価の成果物 (wallfit_series.csv・QUASISTEADY_wallfit.txt) を上書きしないよう --e3 は別名で書く
    with open(rd / "wallfit_series_e3.csv", "w") as fh:
        fh.write("step," + ",".join(cols) + "\n")
        for q in series:
            fh.write(f"{q['step']}," + ",".join(f"{q[c]:.10g}" for c in cols) + "\n")
    tail_frac = (E3_TAIL - 0.5) / len(series)       # check_quasisteady は k = max(3, ceil(tail·n)) 枚 → ちょうど 5 枚
    qs = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(rd / "wallfit_series_e3.csv"),
                         "--series-cols", ",".join(cols), "--tail", f"{tail_frac:.6f}"], capture_output=True, text=True).stdout
    (rd / "QUASISTEADY_wallfit_e3.txt").write_text(qs)
    tail = series[-E3_TAIL:]
    rep = {c: float(np.mean([q[c] for q in tail])) for c in cols}
    T = {}
    for c in cols:
        v = np.array([q[c] for q in tail]); st = np.array([q["step"] for q in tail], float)
        T[c] = max(float(np.ptp(v)), float(abs(np.polyfit(st, v, 1)[0]) * (st[-1] - st[0])))
    qE, dE = quantities(load_field(sf[-1][0], sf[-1][1]), 0.025, fixed_from=0.05 if FIXED else None)
    E = {c: abs(qE[c] - series[-1][c]) for c in cols if c in qE}
    exrep = {v: {k: float(np.mean([e[v][k] for e in exv])) for k in exv[0][v]} for v in exv[0]}
    base_err = max(abs(exrep["base"][k] - rep[k]) for k in ("exit_core_M", "exit_M_dev"))
    cv = sf[-1][0] / "CONVERGENCE_VERDICT.txt"
    conv = cv.read_text() if cv.exists() else "missing"
    conv_overall = [l.strip() for l in conv.splitlines() if "-> " in l and l.startswith("===")]
    info = json.loads((rd / "prepare_info.json").read_text())
    try:
        ev = wall_evidence(rd, (C / other_run) if other_run else None)
    except Exception as e:  # noqa: BLE001 — 証拠が取れないことも記録 (前提未達として総合を保留にする)
        ev = {"status": "error", "error": str(e)}
    rec.update(status="ok", rep=rep, T=T, E=E, quasisteady=_verdict_lines(qs, cols), n_snaps=len(series), last_step=series[-1]["step"],
               tail_steps=[q["step"] for q in tail], tail_frac=tail_frac, exit_variants_tailmean=exrep, exit_base_check_maxabs=base_err,
               convergence=(conv_overall[-1] if conv_overall else conv.strip().splitlines()[-1] if conv != "missing" else "missing"),
               wall_evidence=ev, ic=info.get("ic"), wall_fit_mono_r2=(info.get("wall_fit") or {}).get("mono_r2"),
               _dists=dists)
    return rec


def run_e3():
    global X_E, X_F, WIN_T, WIN_O
    from throat_mono_judge import EXIT_QUANTITIES, delta_key, judge_quantity, overall
    runs_all = ARMS["A"] + ARMS["B"]
    X_E, X_F = e3_geometry(runs_all)
    WIN_T = (X_E + 2.0, X_F - 1.0); WIN_O = (X_E - 15.0, X_F)          # c2final n2400 の定数と同じ取り方 (X_E+2, X_F−1) / (X_E−15, X_F)
    # 出口の共通標本: 腕 A の最初の run の最終場の出口節点 (η 0.05〜0.7)
    r0 = ARMS["A"][0] if ARMS["A"] else ARMS["B"][0]
    F0 = load_field(C / r0, snaps(r0)[-1])
    eta0 = F0["R"][-1] / F0["R"][-1, -1]
    etaS = eta0[(eta0 >= EXIT_BAND[0]) & (eta0 <= EXIT_BAND[1])]
    R_ = {}
    for arm, runs in ARMS.items():
        other = ARMS["B" if arm == "A" else "A"]
        for run in runs:
            R_[run] = e3_run(run, arm, etaS, other[0] if other else None)
    ok = {arm: [r for r in runs if R_[r]["status"] == "ok"] for arm, runs in ARMS.items()}
    cols = list(next(R_[r] for r in R_ if R_[r]["status"] == "ok")["rep"].keys()) if any(ok.values()) else []
    # 出口評価誤差 E_exit,q: 各変種での (B − A) と既定での (B − A) の差の最大
    E_exit, exit_rows = {}, {}
    if ok["A"] and ok["B"]:
        variants = list(R_[ok["A"][0]]["exit_variants_tailmean"].keys())
        for qn in EXIT_QUANTITIES + ("exit_core_M",):
            ba = {}
            for v in variants:
                a = np.mean([R_[r]["exit_variants_tailmean"][v][qn] for r in ok["A"]])
                b = np.mean([R_[r]["exit_variants_tailmean"][v][qn] for r in ok["B"]])
                ba[v] = {"A": float(a), "B": float(b), "B_minus_A": float(b - a)}
            exit_rows[qn] = ba
            if qn in EXIT_QUANTITIES:
                E_exit[qn] = float(max(abs(ba[v]["B_minus_A"] - ba["base"]["B_minus_A"]) for v in variants))
    judge = []
    for c in cols:
        if delta_key(c) is None:
            continue
        a = [R_[r]["rep"][c] for r in ok["A"]]; b = [R_[r]["rep"][c] for r in ok["B"]]
        Tm = max([R_[r]["T"][c] for r in ok["A"] + ok["B"]], default=float("nan"))
        Em = max([R_[r]["E"].get(c, float("nan")) for r in ok["A"] + ok["B"]], default=float("nan"))
        steady = all(R_[r]["quasisteady"].get(c) == "STEADY" for r in ok["A"] + ok["B"])
        judge.append(judge_quantity(c, a, b, Tm, Em, steady, E_exit=E_exit.get(c, 0.0) if c in EXIT_QUANTITIES else 0.0))
    ov = overall(judge)
    # 前提: 全 run の壁の証拠が一致し、腕 B の壁に mono_r2、出口の既定帯が quantities と一致
    pre = {r: {"wall_evidence": R_[r].get("wall_evidence", {}).get("status"), "mono_r2": R_[r].get("wall_fit_mono_r2"),
               "exit_base_check_maxabs": R_[r].get("exit_base_check_maxabs")} for r in R_ if R_[r]["status"] == "ok"}
    pre_ok = (all(v["wall_evidence"] == "consistent" for v in pre.values())
              and all(R_[r].get("wall_fit_mono_r2") is not None for r in ok["B"])
              and all((v["exit_base_check_maxabs"] or 0.0) <= 1e-12 for v in pre.values()))
    verdict = ov["verdict"] if pre_ok else f"保留 (前提未達: 壁の証拠・腕 B の mono_r2・出口既定帯の照合) / 量の判定は {ov['verdict']}"
    # 分布 (残差分布の差) は記録のみ (E3 の判定量ではない)
    dist_rec = {}
    if ok["A"] and ok["B"]:
        for k in R_[ok["A"][0]]["_dists"][0]:
            A_ = np.array([np.mean([d[k][1] for d in R_[r]["_dists"]], axis=0) for r in ok["A"]])
            B_ = np.array([np.mean([d[k][1] for d in R_[r]["_dists"]], axis=0) for r in ok["B"]])
            if A_.shape[1:] == B_.shape[1:]:
                dist_rec[k] = {"B_minus_A_maxabs": float(np.abs(B_.mean(0) - A_.mean(0)).max()), "R_A": float(np.ptp(A_, 0).max()),
                               "R_B": float(np.ptp(B_, 0).max())}
    for r in R_:
        R_[r].pop("_dists", None)
    out = {"plan": "plans/active/tooling-nozzle-throat-monotone-r2.md §6 E1〜E4", "pair": _pair, "arms": ARMS, "fixed_coef": FIXED,
           "geometry": {"X_E": X_E, "X_F": X_F, "WIN_T": list(WIN_T), "WIN_O": list(WIN_O)}, "exit_common_sample_from": r0,
           "runs": R_, "E_exit": E_exit, "exit_variants_B_minus_A": exit_rows, "judge": judge, "overall": ov, "preconditions": pre,
           "preconditions_ok": pre_ok, "VERDICT": verdict, "distributions_record_only": dist_rec}
    (C / "_band_ab").mkdir(exist_ok=True)
    dst = C / f"_band_ab/throat_mono_euler_ab_{_pair[0]}_vs_{_pair[1]}{'_fixedcoef' if FIXED else ''}.json"
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    for row in judge:
        print(f"{row['col']:28s} {row['status']:9s} " + ("" if "U" not in row else
              f"B−A {row['B_minus_A']:+.3e}  U {row['U']:.3e}  Δq {row['Delta']:.3g}") + (f"  ({row['reason']})" if "reason" in row else ""))
    print(f"VERDICT: {verdict}  -> {dst}")


if E3:
    if not any(a.startswith("--pair=") for a in sys.argv):
        raise SystemExit("--e3 は --pair=A,B を明示する (既定の interp,fit の run に e3 の成果物を書かない)")
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    run_e3()
    sys.exit(0)

R = {}
for arm, runs in ARMS.items():
    for run in runs:
        rd = C / run
        sf = series_files(run)
        series = []
        for d_, f, st_ in sf:
            F = load_field(d_, f); q, _ = quantities(F, 0.05); q["step"] = st_; series.append(q)
        cols = [k for k in series[0] if k != "step"]
        with open(rd / "wallfit_series.csv", "w") as fh:
            fh.write("step," + ",".join(cols) + "\n")
            for q in series:
                fh.write(f"{q['step']}," + ",".join(f"{q[c]:.10g}" for c in cols) + "\n")
        qs = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(rd / "wallfit_series.csv"),
                             "--series-cols", ",".join(cols)], capture_output=True, text=True).stdout
        (rd / "QUASISTEADY_wallfit.txt").write_text(qs)
        verd = {c: re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", qs, re.M).group(1) for c in cols}
        tail = series[-5:]
        rep = {c: float(np.mean([q[c] for q in tail])) for c in cols}
        T = {}
        for c in cols:
            v = np.array([q[c] for q in tail]); st = np.array([q["step"] for q in tail], float)
            T[c] = max(float(np.ptp(v)), float(abs(np.polyfit(st, v, 1)[0]) * (st[-1] - st[0])))
        last = sf[-5:]
        dists = [quantities(load_field(d_, f), 0.05)[1] for d_, f, _ in last]
        dmean = {k: (dists[0][k][0], np.mean([d[k][1] for d in dists], axis=0)) for k in dists[0]}
        qE, _ = quantities(load_field(sf[-1][0], sf[-1][1]), 0.025, fixed_from=0.05 if FIXED else None); q5 = series[-1]
        E = {c: abs(qE[c] - q5[c]) for c in cols if c in qE}
        cv = (sf[-1][0] / "CONVERGENCE_VERDICT.txt")
        conv = cv.read_text() if cv.exists() else "missing"
        conv_overall = [l.strip() for l in conv.splitlines() if "-> " in l and l.startswith("===")]
        R[run] = dict(arm=arm, rep=rep, T=T, E=E, quasisteady=verd, dist=dmean, n_snaps=len(series), last_step=series[-1]["step"],
                      convergence=conv_overall[-1] if conv_overall else conv.strip().splitlines()[-1] if conv != "missing" else "missing",
                      convergence_file=str(cv.relative_to(C)) if cv.exists() else None)
        # 分布の時間変動 T_dist (同じ x ごとの末尾 5 枚の ptp と線形ドリフトの最大) と刻み感度 E_dist (0.025 で同じ曲線を評価し 0.05 の点へ戻した差)
        st5 = np.array([st_ for _, _, st_ in last], float)
        R[run]["T_dist"] = {k: float(max(np.ptp(np.array([d[k][1] for d in dists]), axis=0).max(),
                                         (np.abs(np.polyfit(st5, np.array([d[k][1] for d in dists]), 1)[0]) * (st5[-1] - st5[0])).max()))
                            for k in dists[0]}
        _, dE = quantities(load_field(sf[-1][0], sf[-1][1]), 0.025, fixed_from=0.05 if FIXED else None)
        R[run]["E_dist"] = {k: float(np.abs(np.interp(dists[-1][k][0], dE[k][0], dE[k][1]) - dists[-1][k][1]).max()) for k in dE}

# メッシュ照合: B が使われた証拠 (準備ディレクトリの壁節点)
mesh_check = {}
try:
    W = {}
    for arm in _pair:
        pd = C / f"_prep_wallfit_{arm}"; info = json.loads((pd / "prepare_info.json").read_text())
        with h5py.File(pd / "nozzle.h5") as f:
            nc = f["/MESH/COORD"][:].reshape(-1, 3)
        ni = int(info["mesh"]["ni"]); S = float(info["scale_m"]); nj = nc.shape[0] // ni
        W[arm] = ((nc[:, 0] / S).reshape(ni, nj)[:, -1], (nc[:, 1] / S).reshape(ni, nj)[:, -1])
    xa, ra = W[_pair[0]]; xb, rb = W[_pair[1]]
    mesh_check = dict(same_x=bool(np.allclose(xa, xb, atol=1e-9)), wall_dr_max=float(np.abs(rb - ra).max()),
                      x_at_max=float(xa[np.argmax(np.abs(rb - ra))]))
except Exception as e:  # noqa: BLE001
    mesh_check = {"error": str(e)}

scal = list(next(iter(R.values()))["rep"].keys())
out = {"runs": {k: {kk: vv for kk, vv in v.items() if kk != "dist"} for k, v in R.items()}, "mesh_check": mesh_check, "delta": DELTA}
dkey = lambda c: ("M_wave" if c.startswith("M_wave") else "P_wave" if c.startswith("P_wave") else "overshoot_exitnorm" if c.startswith("overshoot_exitnorm")
                  else "overshoot" if c.startswith("overshoot") else "P_slope_abs" if c.startswith("P_slope_abs") else "exit_M_dev" if c == "exit_M_dev"
                  else "exit_core_M" if c == "exit_core_M" else None)
judge = {}
for c in scal:
    dk = dkey(c)
    a = np.array([R[r]["rep"][c] for r in ARMS["A"]]); b = np.array([R[r]["rep"][c] for r in ARMS["B"]])
    T = max(R[r]["T"][c] for r in ARMS["A"] + ARMS["B"]); E = max(R[r]["E"].get(c, 0.0) for r in ARMS["A"] + ARMS["B"])
    row = dict(A_mean=float(a.mean()), B_mean=float(b.mean()), B_minus_A=float(b.mean() - a.mean()), R_A=float(np.ptp(a)), R_B=float(np.ptp(b)), T=T, E=E)
    if dk:
        D = DELTA[dk]; U = max(3 * row["R_A"], 3 * row["R_B"], 2 * T, 2 * E, D / 10)
        row.update(Delta=D, U=U, tail_ok=bool(T <= D / 4), step_ok=bool(E <= D / 10))
        row["judge"] = ("small" if (U <= D / 2 and abs(row["B_minus_A"]) + U <= D) else "different" if abs(row["B_minus_A"]) - U > D else "hold")
        if dk != "exit_core_M":   # 大きいほど悪い量: 非劣化 (B−A+U ≤ Δ)・改善 (B−A ≤ −U)・悪化 (B−A ≥ +U) (codex plan M1)
            d_ = row["B_minus_A"]
            row["noninferior"] = bool(d_ + U <= D); row["direction"] = "improved" if d_ <= -U else "worse" if d_ >= U else "within_U"
    judge[c] = row
for k in ("M_resid_eta0.0", "M_resid_eta0.1", "P_resid_eta0.0", "P_resid_eta0.1"):
    dk = "M_resid_diff" if k.startswith("M") else "P_resid_diff"; D = DELTA[dk]
    A_ = np.array([R[r]["dist"][k][1] for r in ARMS["A"]]); B_ = np.array([R[r]["dist"][k][1] for r in ARMS["B"]])
    diff = float(np.abs(B_.mean(0) - A_.mean(0)).max()); RA_ = float(np.ptp(A_, 0).max()); RB_ = float(np.ptp(B_, 0).max())
    Td = max(R[r]["T_dist"][k] for r in ARMS["A"] + ARMS["B"]); Ed = max(R[r]["E_dist"][k] for r in ARMS["A"] + ARMS["B"])
    U = max(3 * RA_, 3 * RB_, 2 * Td, 2 * Ed, D / 10)
    judge[k + "_dist"] = dict(B_minus_A_maxabs=diff, R_A=RA_, R_B=RB_, T=Td, E=Ed, Delta=D, U=U, tail_ok=bool(Td <= D / 4), step_ok=bool(Ed <= D / 10),
                              judge=("small" if (U <= D / 2 and diff + U <= D) else "different" if diff - U > D else "hold"))
for k in ("P_raw_eta0.0", "P_raw_eta0.1"):
    A_ = np.array([R[r]["dist"][k][1] for r in ARMS["A"]]); B_ = np.array([R[r]["dist"][k][1] for r in ARMS["B"]])
    judge[k + "_dist"] = dict(B_minus_A_maxabs=float(np.abs(B_.mean(0) - A_.mean(0)).max()), note="併記のみ")
qs_all = {r: all(v == "STEADY" for v in R[r]["quasisteady"].values()) for r in R}
out["judge"] = judge; out["all_steady"] = qs_all
gated = [v["judge"] for v in judge.values() if "judge" in v]
tail_ok = all(v.get("tail_ok", True) and v.get("step_ok", True) for v in judge.values())
out["overall"] = ("保留 (準定常・刻み精度の前提未達)" if not (all(qs_all.values()) and tail_ok) else
                  "支持: 両壁の差は小さい" if all(g == "small" for g in gated) else
                  "棄却: 差あり (向きを確認)" if any(g == "different" for g in gated) else "保留")
(C / "_band_ab").mkdir(exist_ok=True)
(C / ("_band_ab/wallfit_euler_ab" + ("" if _pair == ["interp", "fit"] else f"_{_pair[0]}_vs_{_pair[1]}") + ("_diag_fixedcoef.json" if FIXED else "_v2.json"))).write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps({"mesh_check": mesh_check, "all_steady": qs_all, "overall": out["overall"],
                  "judge": {k: {kk: (round(vv, 7) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in judge.items()}}, indent=1, ensure_ascii=False))
