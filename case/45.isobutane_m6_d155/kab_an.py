"""plan axisymmetric-freestream-hoop-gauge §4.13 (事前登録) の判定: キーだけを変えた 60,000 step の A/B の出力ごとの揺れ。kab.sh の後に回す。
A = run_0485_kab_k1 (キー 1)、B = run_0486_kab_k0 (キー 0)。各出力 (2,500〜60,000 step) から cold_series.snapshot で Q_w と θ_r(40/70/94) を作る。
窓: 主 20,000〜60,000 (17 点)、前半 20,000〜40,000・後半 40,000〜60,000 (各 9 点)。物差し: J = RMS(qᵢ − qᵢ₋₁)/|平均|、傾きを除いた RMS/|平均|、
(最大 − 最小)/|平均|、ドリフト (窓の両端の差/平均)。分岐は Q_w で決める (θ_r は記録):
- 支持 (暫定): 主窓の J_A ≥ 0.0385 % (run_0483 の末尾 0.077 % の半分) かつ、B/A が J と傾きを除いた RMS の両方で 3 窓すべて ≤ 0.5。
- 支持しない: J_A が維持され、上の半減が 3 窓で揃わない。
- 判別不能: J_A < 0.0385 %、または窓によって B/A ≤ 0.5 と B/A > 1 が入れ替わる、または非有限・欠損。
結果は標準出力と _band_ab/cold_pair/kab_judge.json。"""
import json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ARMS = {"A_k1": HERE / "run_0485_kab_k1", "B_k0": HERE / "run_0486_kab_k0"}
STEPS = list(range(2500, 60001, 2500))
WIN = {"主 20k〜60k": (20000, 60000), "前半 20k〜40k": (20000, 40000), "後半 40k〜60k": (40000, 60000)}
QTY = ("Q_w", "theta_r_40", "theta_r_70", "theta_r_94")
J_REF_HALF = 0.0385
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
yb_x, yb = XC.common_yb()
OUT = {"undecidable": [], "series": {}, "metrics": {}}
for a, run in ARMS.items():
    ser = {}
    for s in STEPS:
        if not (run / f"res_{s}.h5").exists(): OUT["undecidable"].append(f"{a}: res_{s} が無い"); break
        r = CS.snapshot(run, s, yb_x, yb)
        if r.get("nonfinite") or any(not np.isfinite(r[k]) for k in QTY): OUT["undecidable"].append(f"{a}: step {s} に非有限"); break
        ser[s] = {k: float(r[k]) for k in QTY}
    OUT["series"][a] = ser
def metr(ser, w):
    st = [s for s in STEPS if w[0] <= s <= w[1] and s in ser]
    out = {}
    for k in QTY:
        v = np.array([ser[s][k] for s in st]); m = abs(v.mean())
        if len(v) < 3 or m == 0 or not np.isfinite(m): out[k] = None; continue
        d = np.diff(v); x = np.arange(len(v)); res = v - np.polyval(np.polyfit(x, v, 1), x)
        out[k] = dict(n=len(v), J=100 * np.sqrt(np.mean(d ** 2)) / m, detr=100 * np.sqrt(np.mean(res ** 2)) / m,
                      ptp=100 * np.ptp(v) / m, drift=100 * (v[-1] - v[0]) / m)
    return out
if not OUT["undecidable"]:
    for a in ARMS: OUT["metrics"][a] = {wn: metr(OUT["series"][a], w) for wn, w in WIN.items()}
    print("窓            量          A(キー1) J / 傾き除き / 幅 / ドリフト [%]      B(キー0) J / 傾き除き / 幅 / ドリフト [%]   B/A (J, 傾き除き)")
    ratios = {}
    for wn in WIN:
        for k in QTY:
            ma, mb = OUT["metrics"]["A_k1"][wn][k], OUT["metrics"]["B_k0"][wn][k]
            if ma is None or mb is None: OUT["undecidable"].append(f"{wn} {k}: 物差しが作れない"); continue
            rj = mb["J"] / ma["J"] if ma["J"] > 0 else np.inf; rd = mb["detr"] / ma["detr"] if ma["detr"] > 0 else np.inf
            ratios[(wn, k)] = (rj, rd)
            print(f"{wn:12s}  {k:10s}  {ma['J']:.4f} / {ma['detr']:.4f} / {ma['ptp']:.4f} / {ma['drift']:+.4f}    {mb['J']:.4f} / {mb['detr']:.4f} / {mb['ptp']:.4f} / {mb['drift']:+.4f}    {rj:.2f}, {rd:.2f}")
    OUT["ratios"] = {f"{wn}|{k}": v for (wn, k), v in ratios.items()}
if OUT["undecidable"]:
    verdict = "判別不能: " + "; ".join(OUT["undecidable"][:5])
else:
    jA = OUT["metrics"]["A_k1"]["主 20k〜60k"]["Q_w"]["J"]
    q = [ratios[(wn, "Q_w")] for wn in WIN]
    if jA < J_REF_HALF:
        verdict = f"判別不能 (A が静まった: 主窓の J_A {jA:.4f} % < {J_REF_HALF} %)"
    elif any(max(r) <= 0.5 for r in q) and any(min(r) > 1.0 for r in q):
        verdict = "判別不能 (窓によって B/A ≤ 0.5 と B/A > 1 が入れ替わる)"
    elif all(r[0] <= 0.5 and r[1] <= 0.5 for r in q):
        verdict = "暫定支持 (G0 のこの状態・期間では、キー 1 が Q_w の揺れを 2 倍以上に維持している)"
    else:
        verdict = "支持しない (キーで揺れが半減するとは言えない。キー 0・1 に共通の変更と履歴の候補が残る)"
OUT["verdict"] = verdict
print(f"== VERDICT §4.13: {verdict}")
print("== 記録: 判定ツール")
for a, run in ARMS.items():
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT.setdefault("tools", {})[f"{a} check_convergence"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
    if a in OUT["series"] and OUT["series"][a]:
        sc = run / "kab_series_main.csv"
        with open(sc, "w") as fh:
            fh.write("step," + ",".join(QTY) + "\n")
            for s in STEPS:
                if 20000 <= s <= 60000 and s in OUT["series"][a]: fh.write(f"{s}," + ",".join(f"{OUT['series'][a][s][k]:.10g}" for k in QTY) + "\n")
        p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(sc), "--series-cols", ",".join(QTY), "--tail", "1.0"], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("OVERALL")]
        OUT["tools"][f"{a} check_quasisteady"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"
for k, v in OUT.get("tools", {}).items(): print(f"  {k}: {v}")
o = HERE / "_band_ab/cold_pair/kab_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
