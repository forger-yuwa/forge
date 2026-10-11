"""plan architecture-float-state-double-geometry §6.36 (事前登録) の判定: FP64 のビルドで step の末尾に状態を float32 に丸める十分性の試験。
入力: _band_ab/cold_pair/sq_series_A.json・sq_series_B.json (sq_series.py が各機械で書く)。出力: _band_ab/cold_pair/sq_judge.json。
腕: A* = run_0595〜0599 (量子化 OFF)、G = run_0620〜0624 (量子化 ON)。
前提 (一つでも欠ければ判定不能):
  1. 10 本とも残差・使う出力 (step 0・20,000〜30,000) が有限、RUN_RC 0、check_convergence が実行できる。
  2. res_0 の 9 保存量が 10 本でビット一致 (旧 A の res_0 との一致は記録)。設定は 10 本とも run_0487 と差 0 行。
  3. G の 5 本は使う全出力で 9 保存量の全 CV が float32 で表せる。自己点検 (step 1 と 30,000、ρY の履歴を含む) の表せない値が 0。
     A* は step 30,000 の表せない値の割合を記録する。
  4. A* の 5 本の Δ がすべての組で |Δ| ≤ 0.1|D|。
判定 (θ_r の 3 断面 × W1 = 20,000〜25,000・W2 = 25,000〜30,000、Δ = (窓の平均 − 自分の step 0) / |A*1 の step 0| × 100、D = (+1.892, +1.647, +1.431) %):
  標本は run ごとの窓の平均。e = sign(D)·(mean(ΔG) − mean(ΔA*))。腕ごとの標本分散から Welch の区間、6 組に Bonferroni (各組 t の分位 1 − 0.05/12、
  自由度 Welch–Satterthwaite、両腕の分散が 0 なら区間は点)。
  結果 A: 6 組すべてで L ≥ 0.5|D|。結果 B: 6 組すべてで [L, U] ⊂ [−0.1|D|, +0.1|D|]。それ以外は判別不能 (有意差が無いだけでは B にしない)。
  Q_w は同じ式の副判定 (D = 2.114 %)。"""
import json, re, sys
from pathlib import Path
import numpy as np
from scipy import stats
HERE = Path(__file__).resolve().parent
CP = HERE / "_band_ab/cold_pair"
ASTAR = [f"run_059{i}_sq_off{i - 4}" for i in range(5, 10)]
G = [f"run_062{i}_sq_on{i + 1}" for i in range(5)]
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
D = {"theta_r_40": 1.892, "theta_r_70": 1.647, "theta_r_94": 1.431, "Q_w": 2.114}
WIN = {"W1": (20000, 25000), "W2": (25000, 30000)}
ALPHA, NPAIR = 0.05, 6
OUT = {"undecidable": [], "pre": {}, "pairs": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
R = {}
for m in ("A", "B"):
    f = CP / f"sq_series_{m}.json"
    if not f.exists(): und(f"{f.name} が無い"); continue
    d = json.loads(f.read_text())
    for k, v in d["runs"].items(): v["machine"] = m; R[k] = v
    if d.get("ref_run_0487_res0_sha256"): OUT["pre"][f"旧 A の res_0 のハッシュ ({m})"] = d["ref_run_0487_res0_sha256"]
miss = [r for r in ASTAR + G if r not in R]
if miss: und(f"run の記録が無い: {miss}")
if not OUT["undecidable"]:
    # 1
    for r in ASTAR + G:
        if R[r]["problems"]: und(f"{r}: {R[r]['problems']}")
        OUT["pre"][f"{r} ({R[r]['machine']}) check_convergence"] = R[r]["check_convergence"]
    # 2
    hs = {R[r]["res0_q9_sha256"] for r in ASTAR + G}
    if len(hs) != 1: und(f"res_0 の 9 保存量が 10 本で一致しない ({len(hs)} 通り)")
    OUT["pre"]["res_0 のハッシュ (10 本)"] = sorted(hs)
    ref = [v for k, v in OUT["pre"].items() if k.startswith("旧 A の res_0")]
    OUT["pre"]["旧 A の res_0 と一致 (記録)"] = (ref[0] in hs) if ref else None
    for r in ASTAR + G:
        if R[r]["config_diff_lines"] != 0: und(f"{r}: run_0487 と設定が {R[r]['config_diff_lines']} 行違う")
    # 3
    for r in G:
        bad = {s: v["total"] for s, v in R[r]["not_float32_representable"].items() if v["total"] != 0}
        if bad: und(f"{r}: 量子化したはずの出力に float32 で表せない値 {bad}")
        lines = R[r]["quantize_log"]
        chk = [l for l in lines if re.search(r"step \d+:", l) and "表せない" in l]
        steps = {int(m.group(1)) for l in chk for m in [re.search(r"step (\d+):", l)] if m}
        nums = [int(x) for l in chk for x in re.findall(r"(?<![\w.])(\d+)(?![\w.])", l.split("表せない", 1)[1])]
        if not ({1, 30000} <= steps): und(f"{r}: 自己点検の行が step 1・30,000 の両方に無い ({sorted(steps)})")
        if any(nums): und(f"{r}: 自己点検で表せない値がある ({chk[:2]})")
        OUT["pre"][f"{r} 自己点検"] = chk[:6]
    for r in ASTAR:
        v = R[r]["not_float32_representable"].get("30000", {})
        OUT["pre"][f"{r} step 30,000 の表せない値の割合 (記録)"] = (v.get("total", 0) / (9 * v["n"])) if v.get("n") else None
# 系列と Δ
def delta(r, k, w, base):
    ser = {int(s): v for s, v in R[r]["series"].items()}
    w0, w1 = WIN[w]
    return 100 * (np.mean([ser[s][k] for s in ser if w0 <= s <= w1]) - ser[0][k]) / abs(base)
verdict = None
if not OUT["undecidable"]:
    base = {k: {int(s): v for s, v in R[ASTAR[0]]["series"].items()}[0][k] for k in KEYS + (QW,)}
    for w in WIN:
        for k in KEYS:
            bad = [r for r in ASTAR if abs(delta(r, k, w, base[k])) > 0.1 * D[k]]
            if bad: und(f"前提 4: A* の {bad} が {w} {k} で |Δ| > 0.1|D|")
if not OUT["undecidable"]:
    tq = 1 - ALPHA / (2 * NPAIR); ra, rb = [], []
    print("窓  断面         mean ΔA*  sd ΔA*   mean ΔG   sd ΔG    e        df     t       L        U        0.5|D|  0.1|D|  判定")
    for w in WIN:
        for k in KEYS + (QW,):
            a = np.array([delta(r, k, w, base[k]) for r in ASTAR]); g = np.array([delta(r, k, w, base[k]) for r in G])
            sg = np.sign(D[k]); Dk = abs(D[k])
            e = sg * (g.mean() - a.mean()); va, vg = a.var(ddof=1) / len(a), g.var(ddof=1) / len(g); se = np.sqrt(va + vg)
            if se == 0: df, t = float("inf"), 0.0
            else:
                df = (va + vg) ** 2 / ((va ** 2 / (len(a) - 1) if va else 0) + (vg ** 2 / (len(g) - 1) if vg else 0)); t = float(stats.t.ppf(tq, df))
            L, U = e - t * se, e + t * se
            A_ = L >= 0.5 * Dk; B_ = (L >= -0.1 * Dk) and (U <= 0.1 * Dk)
            OUT["pairs"][f"{w}|{k}"] = dict(dA=a.tolist(), dG=g.tolist(), e=e, se=se, df=df, t=t, L=L, U=U, resultA=bool(A_), resultB=bool(B_))
            print(f"{w}  {k:11s} {a.mean():+8.3f} {a.std(ddof=1):7.3f}  {g.mean():+8.3f} {g.std(ddof=1):7.3f}  {e:+7.3f}  {df:6.1f} {t:6.2f}  {L:+7.3f}  {U:+7.3f}  {0.5 * Dk:.3f}  {0.1 * Dk:.3f}  {'A' if A_ else ''}{'B' if B_ else ''}{' (副)' if k == QW else ''}")
            if k != QW: ra.append(A_); rb.append(B_)
    verdict = ("結果 A: この量子化の介入 (step 末尾の 9 保存量の float32 化) だけで、観測の差の少なくとも半分の規模を生める (十分性。必要性・float の真因は確定しない)" if all(ra)
               else "結果 B: この位置・頻度の量子化だけでは偏りを生まない (次の対象は途中の丸め・演算との組み合わせ)" if all(rb) else "判別不能")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.36: {verdict}")
(CP / "sq_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→", CP / "sq_judge.json")
