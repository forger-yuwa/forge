"""plan architecture-float-state-double-geometry §6.28 (事前登録) の判定: SST の輸送の更新を止めた (FORGE_FREEZE_TURB=1) 軌跡の A/B。fz.sh の後に B で回す。
A = run_0491_fz_fp64、B = run_0492_fz_f32、B' = run_0493_fz_f32b (共通の Q32 の起点、各 30,000 step・500 ごと)。
前提 (欠ければ判定不能): 全期間の残差の全列・全出力の全データセットが有限、正式ツールが実行できる、内部の ρk・ρω が起点と step 30,000 で変わらない (凍結が効いている)。
主判定 (θ_r の 3 断面、W1 = 20,000〜25,000・W2 = 25,000〜30,000): Δ = (窓の平均 − 自分の起点)/A の起点 ×100、δB = ΔB − ΔA、δB' = ΔB' − ΔA、D = (+1.892, +1.647, +1.431) %。
  差が残る: 両窓・3 断面とも D の向きに δB・δB' ≥ 0.5|D|、かつ |mean(δB,δB')| > 10|δB − δB'| → 第 2 仮説を支持。
  差が抑えられる: 両窓・3 断面とも max(|δB|,|δB'|) + 10|δB − δB'| ≤ 0.1|D| → 「SST の更新を止めても従来の規模の差が育つ」説を退ける。
  それ以外は判別不能。Q_w は副判定として同じ式で記録。結果は _band_ab/cold_pair/fz_judge.json。"""
import csv, json, math, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
ARMS = {"A": HERE / "run_0491_fz_fp64", "B": HERE / "run_0492_fz_f32", "Bp": HERE / "run_0493_fz_f32b"}
STEPS = list(range(0, 30001, 500))
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
NJ = 121
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "series": {}, "delta": {}, "freeze_check": {}, "tools": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
# ---- 前提: 有限性・ツール ----
for a, run in ARMS.items():
    rows = [x for x in csv.DictReader(open(run / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{a}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 30000: und(f"{a}: 30,000 step に届かない")
    for s in STEPS:
        p = run / f"res_{s}.h5"
        if not p.exists():
            if s == 0: continue
            und(f"{a}: res_{s} が無い"); break
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][k].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][k][:])) for k in h["VALUE"]): und(f"{a}: res_{s} に非有限"); break
    for tool, args in (("check_convergence", [str(run)]),):
        p = subprocess.run([sys.executable, str(TOOLS / f"{tool}.py")] + args, capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
        if not ov: und(f"{a}: {tool} が実行できない (rc {p.returncode})")
        OUT["tools"][f"{a} {tool}"] = ov[-1] if ov else None
# ---- 凍結が効いているか (内部の ρk・ρω) ----
with h5py.File(ARMS["A"] / "nozzle.h5", "r") as h: n = h["MESH/COORD"].shape[0] // 3
jj = np.arange(n) % NJ; ii = np.arange(n) // NJ; interior = (jj > 0) & (jj < NJ - 1) & (ii > 0) & (ii < ii.max())
for a, run in ARMS.items():
    s0 = 0 if (run / "res_0.h5").exists() else 500
    with h5py.File(run / f"res_{s0}.h5", "r") as h0, h5py.File(run / "res_30000.h5", "r") as h1:
        rec = {}
        for k in ("roK", "roOmega"):
            x0, x1 = np.asarray(h0["VALUE/" + k][:], np.float64), np.asarray(h1["VALUE/" + k][:], np.float64)
            d = np.abs(x1 - x0)[interior] / np.maximum(np.abs(x0[interior]), 1e-300)
            rec[k] = {"changed_nodes": int(np.count_nonzero(d > 0)), "max_rel": float(d.max())}
        for k in ("k", "omega", "vis_turb"):
            x0, x1 = np.asarray(h0["VALUE/" + k][:], np.float64), np.asarray(h1["VALUE/" + k][:], np.float64)
            rec[k + "_rel_l2"] = float(np.sqrt(np.mean((x1 - x0)[interior] ** 2)) / max(np.sqrt(np.mean(x0[interior] ** 2)), 1e-300))
    OUT["freeze_check"][a] = rec
    print(f"  凍結の確認 {a}: 内部の ρk 変化 {rec['roK']['changed_nodes']} 節点 (最大相対 {rec['roK']['max_rel']:.2e})、ρω {rec['roOmega']['changed_nodes']} 節点 ({rec['roOmega']['max_rel']:.2e})、k・ω・μt の相対 L2 {rec['k_rel_l2']:.2e}・{rec['omega_rel_l2']:.2e}・{rec['vis_turb_rel_l2']:.2e}")
    if rec["roK"]["changed_nodes"] or rec["roOmega"]["changed_nodes"]: und(f"{a}: 内部の ρk・ρω が変わった (凍結が効いていない)")
# ---- 系列と判定 ----
yb_x, yb = XC.common_yb()
if not OUT["undecidable"]:
    for a, run in ARMS.items():
        ser = {}
        for s in STEPS:
            if not (run / f"res_{s}.h5").exists(): continue
            r = CS.snapshot(run, s, yb_x, yb); ser[s] = {k: float(r[k]) for k in KEYS + (QW,)}
        OUT["series"][a] = ser
start = lambda a: OUT["series"][a][min(OUT["series"][a])]
verdict = None
if not OUT["undecidable"]:
    A0 = start("A"); keep, supp = [], []
    print("窓  断面         ΔA       ΔB       ΔB'      δB       δB'      |mean δ|  10|δB−δB'|  0.5|D|  0.1|D|")
    for wn, (w0, w1) in WIN.items():
        for k in KEYS + (QW,):
            dl = {a: 100 * (np.mean([OUT["series"][a][s][k] for s in OUT["series"][a] if w0 <= s <= w1]) - start(a)[k]) / abs(A0[k]) for a in ARMS}
            dB, dBp = dl["B"] - dl["A"], dl["Bp"] - dl["A"]; sg = np.sign(D[k]); Dk = abs(D[k]); md = 0.5 * (dB + dBp)
            k_ = sg * dB >= 0.5 * Dk and sg * dBp >= 0.5 * Dk and abs(md) > 10 * abs(dB - dBp)
            s_ = max(abs(dB), abs(dBp)) + 10 * abs(dB - dBp) <= 0.1 * Dk
            OUT["delta"][f"{wn}|{k}"] = dict(dA=dl["A"], dB=dl["B"], dBp=dl["Bp"], deltaB=dB, deltaBp=dBp, keep=bool(k_), suppressed=bool(s_))
            print(f"{wn}  {k:11s} {dl['A']:+7.3f}  {dl['B']:+7.3f}  {dl['Bp']:+7.3f}  {dB:+7.3f}  {dBp:+7.3f}  {abs(md):7.3f}  {10 * abs(dB - dBp):8.3f}  {0.5 * Dk:6.3f}  {0.1 * Dk:6.3f}  {'残る' if k_ else ''}{'抑えられる' if s_ else ''}{' (副)' if k == QW else ''}")
            if k != QW: keep.append(k_); supp.append(s_)
    verdict = ("差が残る (SST の更新を止めても残る経路だけで偏りが育つ: 第 2 仮説を支持)" if all(keep)
               else "差が抑えられる (この精度差の発達には SST の更新・フィードバックが要った。SST の中の演算が発生源とは言わない)" if all(supp) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:5])
OUT["verdict"] = verdict
print(f"== VERDICT §6.28: {verdict}")
for k, v in OUT["tools"].items(): print(f"  {k}: {v}")
o = HERE / "_band_ab/cold_pair/fz_judge.json"; o.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", o)
