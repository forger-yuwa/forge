"""plan architecture-float-state-double-geometry §6.58 (事前登録、codex diagnose `2026-10-12-g0c-stop-compare-2-diagnose.md`) の判定:
停止規則 (§6.11) での出力の互換性の試験。G0c・FP64・G0c の Q32 の起点から、G1 = 上限 1 (run_0760_gg_c1)、G2a・G2b = 上限 50 (run_0761_gg_c50a・run_0762_gg_c50b)。
参照 (G0c の採用の判断の追加の比較): G0・上限 1・FP64 の T2 の 2 本 (run_0752_t2_d1・run_0754_t2_d2、停止の出力)。収束した解の比較ではない。
指標 (停止の出力どうし、候補 c と参照 r の組ごと、分母は |r|):
  M1 θ_r(40/70/94): 相対差 ≤ 0.1 %。M2 δ_loc(40/70/94): ≤ 0.1 %。M3 代理 Q_w (CP 1360・Pr 0.72・Sutherland の代理の物性): ≤ 0.1 %。
  M4 代理 q_w(x) の局所の相対差: 壁の全節点で |Δq_w| ≤ 0.01·|q_w,r| + 0.01·max|q_w,r| (ゼロの近くの絶対の許容 = 最大の 1 %)。
  M5 符号つきのせん断 τ_w(x): 同じ式で ≤ 1 % + 最大の 1 %、符号の反転の数が同じ。
  M6 出口の列のマッハ数 (軸の行を除く): r/r_e ≤ 0.95 の行の最大 |ΔM| ≤ 0.002、r/r_e ≤ 0.8 の行の平均の |ΔM| ≤ 0.001。軸の行の差は別に記録。
  M7 収支: 入口・出口の質量流量・全エンタルピー流量 (VALUE/h0) の相対差 ≤ 0.1 %。停止の出力の欠損は停止規則で ≤ 0.1 kg/s。
  格子が違う組 (G0c 対 G0) は、参照の壁の分布を候補の壁の x へ線形に補間し、出口の列は同じ j の行で比べる (出口の列の節点は区間の外でビット一致)。
判定 (主の比較 G2a/G1・G2b/G1):
  結果 A: 3 本とも停止規則に到達し、両方の組で M1〜M7 がすべて許容内 → 「上限の変更による許容超過」を支持しない (この起点・停止の手順に限る)。
  結果 B: 到達した腕で許容超過がある、または上限 20 万 step までに到達しない腕がある → 「上限 50 がこの運用の要件を満たす」を退ける。
  G2a と G2b で採否が割れたら再現性未確認。欠損・非有限は判別不能。
G0c の採用の判断には、G1・G2a・G2b のそれぞれを T2 の 2 本と直接比べる (6 組)。全組で M1〜M7 が許容内のときだけ「G0 との出力の互換性: 許容内」。
記録: 到達の step、全残差列の check_convergence の VERDICT と落ちた列、θ_r・Q_w の check_quasisteady (末尾 2 万 step)、各組の全指標の値。→ _band_ab/cold_pair/gg_judge.json"""
import csv, json, math, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
CAND = {"G1": "run_0760_gg_c1", "G2a": "run_0761_gg_c50a", "G2b": "run_0762_gg_c50b"}
REF = {"T2a": "run_0752_t2_d1", "T2b": "run_0754_t2_d2"}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94")
OUT = {"undecidable": [], "arms": {}, "pairs": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
yb_x, yb = XC.common_yb()
def reach_of(run):
    p = HERE / run / "m9_watch.json"
    if not p.exists(): return None, None
    d = json.loads(p.read_text()); return d.get("status"), d.get("reach_step")
def fields(run, step):
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
for nm, run in list(CAND.items()) + list(REF.items()):
    st, n = reach_of(run)
    OUT["arms"][nm] = {"run": run, "status": st, "reach_step": n}
    if st != "REACHED":
        if nm in CAND and st == "CENSORED": continue      # 未到達は結果 B の材料 (下で扱う)
        und(f"{nm} ({run}): 見張りの状態 {st}"); continue
    try:
        F[nm] = fields(run, n)
        if F[nm]["nonfinite"]: und(f"{nm}: 停止の出力に非有限 {F[nm]['nonfinite']}")
    except (OSError, KeyError, SystemExit) as e:
        und(f"{nm}: 停止の出力を読めない ({e})")
def compare(c, r):
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
    res["M_exit 最大 |ΔM| (r/r_e ≤ 0.95、軸の行を除く)"] = float(np.max(np.abs(dM[m95]))); res["M_exit 平均 |ΔM| (r/r_e ≤ 0.8)"] = float(np.mean(np.abs(dM[m80])))  # 2026-10-13 codex: |平均| でなく平均 |ΔM| (登録どおり); res["M_exit 軸の行の ΔM (記録)"] = float(dM[0])
    ok &= res["M_exit 最大 |ΔM| (r/r_e ≤ 0.95、軸の行を除く)"] <= 0.002 and res["M_exit 平均 |ΔM| (r/r_e ≤ 0.8)"] <= 0.001
    for k in ("mdot_in", "mdot_out", "Hdot_in", "Hdot_out"):
        v = 100 * (float(a[k]) - float(b[k])) / abs(float(b[k])); res[k] = v; ok &= abs(v) <= 0.1
    res["許容内"] = bool(ok); return res
verdict, adopt = None, None
if not OUT["undecidable"]:
    unreached = [nm for nm in CAND if OUT["arms"][nm]["status"] != "REACHED"]
    if "G1" in unreached: verdict = "結果 B (上限 1 の G1 が上限 20 万 step までに到達しない)"
    else:
        main = {}
        for c in ("G2a", "G2b"):
            if c in F: main[c] = compare(c, "G1"); OUT["pairs"][f"{c}/G1"] = main[c]
        oks = [main[c]["許容内"] for c in main]
        if unreached: verdict = f"結果 B (上限 50 の {', '.join(unreached)} が到達しない)"
        elif len(set(oks)) > 1: verdict = "再現性未確認 (G2a と G2b で採否が割れた)"
        elif all(oks): verdict = "結果 A (G0c で上限 50 と 1 の停止の出力が全指標で許容内)"
        else: verdict = "結果 B (上限 50 と 1 の停止の出力に許容超過がある)"
    ad = []
    for c in CAND:
        if c not in F: continue
        for r in REF:
            OUT["pairs"][f"{c}/{r}"] = compare(c, r); ad.append(OUT["pairs"][f"{c}/{r}"]["許容内"])
    adopt = "G0 との出力の互換性: 許容内 (全組)" if ad and all(ad) else "G0 との出力の互換性: 許容超過あり" if ad else "比較なし"
    for k, v in OUT["pairs"].items():
        print(k, "許容内" if v["許容内"] else "超過", {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "許容内"})
else:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
for nm, run in list(CAND.items()) + list(REF.items()):
    if not (HERE / run / "residual_history.csv").exists(): continue
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(HERE / run)], capture_output=True, text=True)
    lines = (p.stdout + p.stderr).splitlines()
    OUT["records"][f"{nm} check_convergence"] = {"overall": next((l for l in lines if "OVERALL" in l), "(行なし)"), "failed": [l.strip()[:120] for l in lines if "<--" in l]}
    n = OUT["arms"][nm]["reach_step"]
    if n:
        st = json.loads((HERE / run / "m9_watch.json").read_text()); rows = sorted([x for x in st["rows"] if n - 20000 <= x["step"] <= n], key=lambda x: x["step"])
        sc = HERE / run / "gg_series_tail.csv"
        with open(sc, "w") as fh:
            fh.write("step," + ",".join(KEYS + ("Q_w",)) + "\n")
            for x in rows: fh.write(f"{x['step']}," + ",".join(f"{x[k]:.10g}" for k in KEYS + ("Q_w",)) + "\n")
        q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(HERE / run), "--series-csv", str(sc), "--series-cols", ",".join(KEYS + ("Q_w",)), "--tail", "1.0", "--drift", "0.0001", "--osc", "0.0001"], capture_output=True, text=True)
        OUT["records"][f"{nm} check_quasisteady (末尾 2 万 step)"] = [l.strip() for l in q.stdout.splitlines() if ("theta_r_" in l or "Q_w" in l) and ":" in l]
OUT["verdict"] = verdict; OUT["adoption"] = adopt
print(f"== VERDICT §6.58: {verdict}\n== G0c の採用の判断の比較: {adopt}")
(HERE / "_band_ab/cold_pair/gg_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/gg_judge.json")
