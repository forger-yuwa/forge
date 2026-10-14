"""plan architecture-float-state-double-geometry §6.61 (事前登録、codex diagnose `2026-10-13-testg-result-diagnose.md`) の判定: 速さの比較 S。
上限 50 (生産の B0)・FP64 に固定し、格子 (と登録した初期場の移植) だけを変える。S_g0c = G0c `_g0c/q32_init.h5` (run_0804_sp_g0c)、S_g0 = G0 `_tr/q32_init.h5` (run_0805_sp_g0)。
B で GPU を占有して逐次 (G0c → G0)。停止規則 §6.11 (m9_watch.py、上限 20 万 step)。
出力要件 (各腕): 停止規則に到達し、停止の出力に非有限が無く、T2 の 2 本 (G0・上限 1・FP64 の停止の出力) のそれぞれと M1〜M7 (§6.58 と同じ式・許容) が許容内。
時間: T_reach = 起動 (START_EPOCH) から残差の履歴が到達の step に達した時刻まで (timing.tsv の 30 秒ごとの記録を線形に補間、起動・出力・見張りを含む)。
  1 step の時間 = step 2,500 から到達までの傾き (起動を除く)。占有 = 起動から到達まで timing.tsv の forge の数がすべて 1 以下。
判定:
  S_g0c が出力要件を満たさない → 「試験 G の G2 を占有で再現しない」(速さは言わない)。
  S_g0 が上限までに到達しない、または出力要件を満たさない → 「G0c だけが今回の出力要件を満たした」(速度の倍率は出さない。G0 の所要の下限は記録)。
  両腕が満たす: 占有でない腕があれば速さは判定不能。比 ρ = T_reach(G0c) / T_reach(G0)。ρ ≤ 0.9 → 結果 A (この起点・停止の手順での速さの改善を支持)、
    ρ ≥ 1.0 → 結果 B (速さの改善を退ける)、その間 → 未確定 (計時の揺れと節点数の差 5 % の範囲)。
記録: 到達の step、全残差列の check_convergence、θ_r・Q_w の check_quasisteady (末尾 2 万 step)、各組の全指標。冷えた起点からの所要・収束解の一致は言わない。
→ _band_ab/cold_pair/sp_judge.json"""
import csv, json, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
CAND = {"S_g0c": "run_0804_sp_g0c", "S_g0": "run_0805_sp_g0"}
REF = {"T2a": "run_0752_t2_d1", "T2b": "run_0754_t2_d2"}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
OUT = {"undecidable": [], "arms": {}, "pairs": {}, "timing": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
yb_x, yb = XC.common_yb()
def reach_of(run):
    p = HERE / run / "m9_watch.json"
    if not p.exists(): return None, None
    d = json.loads(p.read_text()); return d.get("status"), d.get("reach_step")
def fields(run, step):   # gg_judge.py (§6.58) と同じ
    xy, ro, ux, uy, T, k, h0 = XC.load_forge(HERE / run, step, with_h0=True)
    ni, nj, S = XC.mesh_info(len(ro))
    o = XC.reduce_fields(xy, ro, ux, uy, T, k, False, ni, nj, S, yb_x, yb, h0=h0)
    with h5py.File(HERE / run / f"res_{step}.h5", "r") as h:
        M = (np.hypot(np.asarray(h["VALUE/Ux"][:], float), np.asarray(h["VALUE/Uy"][:], float)) / np.asarray(h["VALUE/sonic"][:], float)).reshape(ni, nj)
        nonfin = sum(int(np.count_nonzero(~np.isfinite(h["VALUE"][q][:]))) for q in h["VALUE"] if h["VALUE"][q].dtype.kind == "f")
    r = (xy[:, 1] / S).reshape(ni, nj)
    o["M_exit"] = M[-1]; o["r_exit"] = r[-1]; o["nonfinite"] = nonfin
    for xs in XC.XS:
        i = int(np.argmin(np.abs(o["x"] - xs))); o[f"theta_r_{int(xs)}"] = float(o["theta_r"][i]); o[f"delta_loc_{int(xs)}"] = float(o["delta_loc"][i])
    return o
F = {}
def compare(c, r):   # gg_judge.py (§6.58、2026-10-13 に平均 |ΔM| を直した版) と同じ
    a, b = F[c], F[r]; res = {}; ok = True
    for k in KEYS: v = 100 * (a[k] - b[k]) / abs(b[k]); res[k] = v; ok &= abs(v) <= 0.1
    for k in ("delta_loc_40", "delta_loc_70", "delta_loc_94"): v = 100 * (a[k] - b[k]) / abs(b[k]); res[k] = v; ok &= abs(v) <= 0.1
    v = 100 * (float(a["Q_w"]) - float(b["Q_w"])) / abs(float(b["Q_w"])); res["Q_w_proxy"] = v; ok &= abs(v) <= 0.1
    for key, nm in (("q_w", "q_w_proxy"), ("tau_w_signed", "tau_w_signed")):
        qr = np.interp(a["x_w"], b["x_w"], b[key]); dq = np.abs(a[key] - qr); tol = 0.01 * np.abs(qr) + 0.01 * np.max(np.abs(b[key]))
        res[f"{nm} 超過の節点"] = int(np.count_nonzero(dq > tol)); res[f"{nm} 最大の相対差 [%]"] = float(100 * np.max(dq / np.maximum(np.abs(qr), 1e-300)))
        ok &= res[f"{nm} 超過の節点"] == 0
    flips = lambda s: int(np.count_nonzero(np.diff(np.sign(s)) != 0))  # noqa: E731
    res["τ_w の符号の反転の数 (候補 / 参照)"] = [flips(a["tau_w_signed"]), flips(b["tau_w_signed"])]; ok &= flips(a["tau_w_signed"]) == flips(b["tau_w_signed"])
    re_ = b["r_exit"][-1]; rr = a["r_exit"] / re_
    m95 = (rr <= 0.95); m95[0] = False; m80 = (rr <= 0.8); m80[0] = False
    dM = a["M_exit"] - b["M_exit"]
    res["M_exit 最大 |ΔM| (r/r_e ≤ 0.95、軸の行を除く)"] = float(np.max(np.abs(dM[m95]))); res["M_exit 平均 |ΔM| (r/r_e ≤ 0.8)"] = float(np.mean(np.abs(dM[m80]))); res["M_exit 軸の行の ΔM (記録)"] = float(dM[0])
    ok &= res["M_exit 最大 |ΔM| (r/r_e ≤ 0.95、軸の行を除く)"] <= 0.002 and res["M_exit 平均 |ΔM| (r/r_e ≤ 0.8)"] <= 0.001
    for k in ("mdot_in", "mdot_out", "Hdot_in", "Hdot_out"):
        v = 100 * (float(a[k]) - float(b[k])) / abs(float(b[k])); res[k] = v; ok &= abs(v) <= 0.1
    res["許容内"] = bool(ok); return res
def timing(run, n_end):
    d = HERE / run
    try:
        t0 = float((d / "START_EPOCH").read_text().strip())
        rows = [r for r in csv.reader(open(d / "timing.tsv"), delimiter="\t") if len(r) >= 3 and r[1].strip()]
    except (OSError, ValueError) as e:
        return None, f"時刻の記録を読めない ({e})"
    t = np.array([float(r[0]) for r in rows]); s = np.array([int(r[1]) for r in rows]); nf = np.array([int(r[2]) for r in rows])
    if len(t) < 3 or s.max() < n_end: return None, f"時刻の記録が到達の step {n_end} まで無い (最大 {s.max() if len(s) else None})"
    k = int(np.argmax(s >= n_end)); tr = t[k] if k == 0 or s[k] == s[k - 1] else t[k - 1] + (t[k] - t[k - 1]) * (n_end - s[k - 1]) / (s[k] - s[k - 1])
    j = int(np.argmax(s >= 2500)); per = (tr - t[j]) / max(n_end - s[j], 1)
    excl = bool(np.all(nf[t <= tr] <= 1))
    others = sorted({r[3].strip() for r in rows if len(r) > 3 and r[3].strip() and float(r[0]) <= tr})
    return {"T_reach_s": float(tr - t0), "per_step_ms": float(1e3 * per), "exclusive": excl, "others": others, "n_end": int(n_end), "samples": int(len(t))}, None
for nm, run in list(CAND.items()) + list(REF.items()):
    st, n = reach_of(run)
    OUT["arms"][nm] = {"run": run, "status": st, "reach_step": n}
    if st != "REACHED":
        if nm in CAND and st == "CENSORED": continue
        und(f"{nm} ({run}): 見張りの状態 {st}"); continue
    try:
        F[nm] = fields(run, n)
        if F[nm]["nonfinite"]: und(f"{nm}: 停止の出力に非有限 {F[nm]['nonfinite']}")
    except (OSError, KeyError, SystemExit) as e:
        und(f"{nm}: 停止の出力を読めない ({e})")
okout = {}
for c in CAND:
    if c not in F: okout[c] = False; continue
    oks = []
    for r in REF:
        if r in F: OUT["pairs"][f"{c}/{r}"] = compare(c, r); oks.append(OUT["pairs"][f"{c}/{r}"]["許容内"])
    okout[c] = bool(oks) and len(oks) == len(REF) and all(oks)
for c, run in CAND.items():
    st, n = OUT["arms"][c]["status"], OUT["arms"][c]["reach_step"]
    end = n if st == "REACHED" else (200000 if st == "CENSORED" else None)
    if end:
        tm, e = timing(run, end); OUT["timing"][c] = tm if tm else {"error": e}
for k, v in OUT["pairs"].items():
    print(k, "許容内" if v["許容内"] else "超過", {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "許容内"})
print("時間", json.dumps(OUT["timing"], ensure_ascii=False))
if OUT["undecidable"]:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
elif not okout["S_g0c"]:
    verdict = "試験 G の G2 を占有で再現しない (S_g0c が到達しないか、T2 との出力要件を外す)。速さは言わない"
elif not okout["S_g0"]:
    tg = OUT["timing"].get("S_g0", {})
    verdict = ("G0c だけが今回の出力要件を満たした (G0 は " + ("上限 20 万 step までに到達しない" if OUT["arms"]["S_g0"]["status"] == "CENSORED" else "停止したが出力要件を外す")
               + f")。速度の倍率は出さない。G0 の記録の時間 {tg.get('T_reach_s')} s")
else:
    a, b = OUT["timing"].get("S_g0c", {}), OUT["timing"].get("S_g0", {})
    if "T_reach_s" not in a or "T_reach_s" not in b: verdict = "速さは判定不能 (時刻の記録が欠ける)"
    elif not (a["exclusive"] and b["exclusive"]): verdict = "速さは判定不能 (GPU を占有できなかった腕がある)"
    else:
        rho = a["T_reach_s"] / b["T_reach_s"]; OUT["rho"] = rho
        verdict = (f"結果 A (ρ = {rho:.3f} ≤ 0.9: この起点・停止の手順で G0c の方が速い)" if rho <= 0.9 else
                   f"結果 B (ρ = {rho:.3f} ≥ 1.0: G0c は速くない)" if rho >= 1.0 else f"未確定 (ρ = {rho:.3f}、0.9〜1.0)")
for nm, run in CAND.items():
    if not (HERE / run / "residual_history.csv").exists(): continue
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / run)], capture_output=True, text=True)
    lines = (p.stdout + p.stderr).splitlines()
    OUT["records"][f"{nm} check_convergence"] = {"overall": next((l for l in lines if "OVERALL" in l), "(行なし)"), "failed": [l.strip()[:120] for l in lines if "<--" in l]}
    n = OUT["arms"][nm]["reach_step"]
    if n:
        st = json.loads((HERE / run / "m9_watch.json").read_text()); rows = sorted([x for x in st["rows"] if n - 20000 <= x["step"] <= n], key=lambda x: x["step"])
        sc = HERE / run / "sp_series_tail.csv"
        with open(sc, "w") as fh:
            fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
            for x in rows: fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
        q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / run), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + ("Q_w",)), "--tail", "1.0", "--drift", "0.0001", "--osc", "0.0001"], capture_output=True, text=True)
        OUT["records"][f"{nm} check_quasisteady (末尾 2 万 step)"] = [l.strip() for l in q.stdout.splitlines() if ("theta_r_" in l or "Q_w" in l) and ":" in l]
OUT["verdict"] = verdict
print(f"== VERDICT §6.61: {verdict}")
(HERE / "_band_ab/cold_pair/sp_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/sp_judge.json")
