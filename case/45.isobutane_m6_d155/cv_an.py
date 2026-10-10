"""plan axisymmetric-freestream-hoop-gauge §4.14 (事前登録) の判定: 旧い変換と新しい変換の格子の A/B (キー 0、同じビルド・同じ保存量)。
A の上で回す。旧 = run_0490_cv_old_k0 (A、ここで系列を作る)、新 = run_0486_kab_k0 (B、kab_judge.json の系列 B_k0 を読む。B から _band_ab/cold_pair/kab_judge.json を写しておく)。
窓: 20,000〜40,000・40,000〜60,000 step。主指標: Q_w の J = RMS(出力ごとの差分)/|窓の平均|。
両区間で J_新/J_旧 ≥ 2 → 「変換器の変更だけで揺れが大きく増える」を支持。両区間で ≤ 1.25 → 「変換器の変更が主因」を退ける。それ以外・欠損・非有限 → 判定不能。
θ_r の J・傾きを除いた RMS・窓の平均を併記する。結果は _band_ab/cold_pair/cv_judge.json。"""
import json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
OLD = HERE / "run_0490_cv_old_k0"
STEPS = list(range(2500, 60001, 2500)); QTY = ("Q_w", "theta_r_40", "theta_r_70", "theta_r_94")
WIN = {"20k〜40k": (20000, 40000), "40k〜60k": (40000, 60000)}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": []}
yb_x, yb = XC.common_yb()
old = {}
for s in STEPS:
    if not (OLD / f"res_{s}.h5").exists(): OUT["undecidable"].append(f"旧: res_{s} が無い"); break
    r = CS.snapshot(OLD, s, yb_x, yb)
    if r.get("nonfinite") or any(not np.isfinite(r[k]) for k in QTY): OUT["undecidable"].append(f"旧: step {s} に非有限"); break
    old[s] = {k: float(r[k]) for k in QTY}
kj = HERE / "_band_ab/cold_pair/kab_judge.json"
new = {int(s): v for s, v in json.loads(kj.read_text())["series"]["B_k0"].items()} if kj.exists() else {}
if not new: OUT["undecidable"].append("新 (run_0486) の系列が無い (kab_judge.json を写す)")
OUT["series"] = {"old_run_0490": old, "new_run_0486": new}
def metr(ser, w, k):
    v = np.array([ser[s][k] for s in STEPS if w[0] <= s <= w[1] and s in ser]); m = abs(v.mean())
    if len(v) != 9 or m == 0: return None
    x = np.arange(len(v)); res = v - np.polyval(np.polyfit(x, v, 1), x)
    return dict(J=100 * np.sqrt(np.mean(np.diff(v) ** 2)) / m, detr=100 * np.sqrt(np.mean(res ** 2)) / m, mean=float(v.mean()))
ratios = []
if not OUT["undecidable"]:
    print("区間        量          旧 J / 傾き除き / 平均            新 J / 傾き除き / 平均            新/旧 J")
    for wn, w in WIN.items():
        for k in QTY:
            mo, mn = metr(old, w, k), metr(new, w, k)
            if mo is None or mn is None: OUT["undecidable"].append(f"{wn} {k}: 9 点そろわない"); continue
            r = mn["J"] / mo["J"] if mo["J"] > 0 else np.inf
            OUT.setdefault("metrics", {})[f"{wn}|{k}"] = {"old": mo, "new": mn, "ratio_J": r}
            print(f"{wn:10s}  {k:10s}  {mo['J']:.4f} / {mo['detr']:.4f} / {mo['mean']:.6g}    {mn['J']:.4f} / {mn['detr']:.4f} / {mn['mean']:.6g}    {r:.2f}")
            if k == "Q_w": ratios.append(r)
if OUT["undecidable"] or len(ratios) != 2:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
elif all(r >= 2 for r in ratios):
    verdict = "支持 (変換器の変更だけで揺れが大きく増える)"
elif all(r <= 1.25 for r in ratios):
    verdict = "退ける (変換器の変更は主因でない)"
else:
    verdict = "判定不能 (中間、または区間で逆転)"
OUT["verdict"] = verdict
print(f"== VERDICT §4.14: {verdict}  (Q_w の新/旧 J: {', '.join(f'{r:.2f}' for r in ratios)})")
p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(OLD)], capture_output=True, text=True)
ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
OUT["old_check_convergence"] = ov[-1] if ov else f"(行なし) rc={p.returncode}"; print("  旧 check_convergence:", OUT["old_check_convergence"])
o = HERE / "_band_ab/cold_pair/cv_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
