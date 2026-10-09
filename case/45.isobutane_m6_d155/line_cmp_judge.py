"""ライン解の比較の判定 (plan time_integration-line-implicit-speed §6.4 (1)・§6.7 (1)、codex 2026-10-09 plan-3・2026-10-10 plan の反映)。
FORGE_LINE_COMPARE=1 の forge_run.log から、番号ごとに行がちょうど 1 本ずつそろっていることを確かめてから判定する。
  factor N (1..F): [lineCompare] factor N … 失敗の不一致 M / [lineNonfinite] factor N: … / (--bitwise かつ LAYOUT) [lineLayoutCmp] factor N: …
  solve N (1..S): [lineCompare] solve N: dq 最大差 … (不一致 …); 成分ごと … (非有限の累計 …) / [lineNonfinite] solve N: … / [lineEta] solve N <腕> と LU
FAIL: 非有限 (どの行でも) 1 件以上、分解の失敗の不一致・どちらかの腕の失敗、腕の η > 上限で従来の η ≤ 1e-11、(--dq-limit) 補正の差、
      (--bitwise) dq の差 ≠ 0 か不一致 (ビット列) ≠ 0、並べ替えた因子の不一致 ≠ 0。
INDETERMINATE: 行の欠け・重複・読めない数値、評価本数 ≠ ライン数、腕の η > 上限で従来の η も 1e-11 超。
終了コード: 0 = PASS、1 = FAIL、2 = INDETERMINATE。結果は <run>/cmp_judge_<腕>.json。"""
import argparse, json, math, re, sys
from collections import defaultdict
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--arm", required=True); ap.add_argument("--eta-limit", type=float, required=True)
ap.add_argument("--dq-limit", type=float, default=-1.0); ap.add_argument("--bitwise", action="store_true")
ap.add_argument("--factors", type=int, default=20); ap.add_argument("--solves", type=int, default=100); ap.add_argument("--lines", type=int, default=4719)
ap.add_argument("--ro-ref", type=float, default=0.8739869154); ap.add_argument("--a-ref", type=float, default=359.7712185)
a = ap.parse_args()
run = Path(a.run); log = (run / "forge_run.log").read_text(errors="replace")
NUM = r"(\S+)"
F, I = [], []                        # FAIL / INDETERMINATE の理由
def fnum(x, what):
    try:
        v = float(x)
    except ValueError:
        I.append(f"{what}: 数値が読めない ({x})"); return math.nan
    if not math.isfinite(v): F.append(f"{what}: 非有限 ({x})")
    return v
def one(pattern, n, what):
    """番号 1..n の行がちょうど 1 本ずつあるか。{番号: match} を返す。"""
    got = defaultdict(list)
    for m in re.finditer(pattern, log, re.M): got[int(m.group(1))].append(m)
    miss = [k for k in range(1, n + 1) if len(got[k]) != 1]
    extra = [k for k in got if not (1 <= k <= n)]
    if miss or extra: I.append(f"{what}: 番号ごとに 1 本でない (欠け・重複 {miss[:5]}…、範囲外 {extra[:5]})")
    return {k: v[0] for k, v in got.items() if len(v) == 1}
fac = one(r"^\[lineCompare\] factor (\d+)\b.*失敗の不一致 (\d+)$", a.factors, "factor の比較の行")
fnf = one(r"^\[lineNonfinite\] factor (\d+): 腕 LU (\d+)・W (\d+)、従来 LU (\d+)・W (\d+)$", a.factors, "factor の非有限の行")
sol = one(r"^\[lineCompare\] solve (\d+): dq 最大差 " + NUM + " / 最大 " + NUM + r" \(不一致 (\d+)\); 成分ごと 差/最大 = " + " ".join([NUM + "/" + NUM] * 5) + r" \(非有限の累計 (\d+)\)$", a.solves, "solve の比較の行")
snf = one(r"^\[lineNonfinite\] solve (\d+): 腕 y (\d+)、従来 y (\d+)$", a.solves, "solve の非有限の行")
eta = {arm: one(r"^\[lineEta\] solve (\d+) " + arm + r": η 最大 " + NUM + r" \(ライン (-?\d+)\)、評価 (\d+) 本、1e-11 超 (\d+) 本、非有限 (\d+) 本、分解の失敗 (\d+) 本", a.solves, f"η の行 ({arm})")
       for arm in (a.arm, "LU")}
nf = sum(int(x) for x in re.findall(r"^\[lineCompare\] 非有限 (\d+) 件", log, re.M))
nf += sum(sum(int(g) for g in m.groups()[1:]) for m in fnf.values()) + sum(sum(int(g) for g in m.groups()[1:]) for m in snf.values())
fmis = sum(int(m.group(2)) for m in fac.values())
if fmis: F.append(f"分解の失敗の不一致 {fmis}")
sc = [a.ro_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref ** 2]
wp = [0.0] * 5; wr = [0.0] * 5; bw_bad = []
for k, m in sol.items():
    d0, mx0, nm = fnum(m.group(2), f"solve {k} の最大差"), fnum(m.group(3), f"solve {k} の最大"), int(m.group(4))
    v = [fnum(x, f"solve {k} の成分") for x in m.groups()[4:14]]
    nf += int(m.group(15))
    for c in range(5):
        wp[c] = max(wp[c], v[2 * c] / sc[c]); wr[c] = max(wr[c], v[2 * c] / v[2 * c + 1] if v[2 * c + 1] > 0 else 0.0)
    if a.bitwise and (d0 != 0.0 or nm != 0): bw_bad.append((k, d0, nm))
if bw_bad: F.append(f"dq がビットで一致しない solve {len(bw_bad)} 回 (最初: solve {bw_bad[0][0]}・最大差 {bw_bad[0][1]:.2e}・不一致 {bw_bad[0][2]})")
if a.dq_limit >= 0 and max(wp) > a.dq_limit: F.append(f"補正の差 {max(wp):.2e} > {a.dq_limit:.0e}")
E = {}
for arm in (a.arm, "LU"):
    E[arm] = {}
    for k, m in eta[arm].items():
        r = {"max": fnum(m.group(2), f"η ({arm}, solve {k})"), "line": int(m.group(3)), "n_eval": int(m.group(4)), "n_over": int(m.group(5)), "n_nonfinite": int(m.group(6)), "n_fail": int(m.group(7))}
        E[arm][k] = r; nf += r["n_nonfinite"]
        if r["n_fail"]: F.append(f"{arm} の分解の失敗 {r['n_fail']} 本 (solve {k})")
        elif r["n_eval"] != a.lines: I.append(f"{arm} の評価本数 {r['n_eval']} ≠ {a.lines} (solve {k})")
for k, r in E[a.arm].items():
    if r["max"] > a.eta_limit:
        lu = E["LU"].get(k)
        if lu is not None and lu["max"] > 1e-11: I.append(f"solve {k}: 腕の η {r['max']:.3e} > 上限、従来も {lu['max']:.3e} > 1e-11")
        else: F.append(f"solve {k}: 腕の η {r['max']:.3e} > {a.eta_limit:.0e}")
if a.bitwise and a.arm == "LAYOUT":
    lay = one(r"^\[lineLayoutCmp\] factor (\d+): LU 不一致 (\d+)・W 不一致 (\d+)・ピボット 不一致 (\d+)・失敗 不一致 (\d+)$", a.factors, "並べ替えた因子の比較の行")
    bad = {k: [int(g) for g in m.groups()[1:]] for k, m in lay.items() if any(int(g) for g in m.groups()[1:])}
    if bad: F.append(f"並べ替えた因子が従来の因子とビットで一致しない factor {sorted(bad)[:5]} (最初: {bad[min(bad)]})")
if nf: F.append(f"非有限 {nf} 件")
verdict = "FAIL" if F else ("INDETERMINATE" if I else "PASS")
out = {"run": run.name, "arm": a.arm, "verdict": verdict, "fail": F, "indeterminate": I, "nonfinite": nf,
       "eta_arm_max": max((r["max"] for r in E[a.arm].values()), default=None), "eta_LU_max": max((r["max"] for r in E["LU"].values()), default=None),
       "max_dq_diff_over_scale": wp, "max_dq_diff_over_max": wr, "n_factor": len(fac), "n_solve": len(sol)}
(run / f"cmp_judge_{a.arm}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
sys.exit({"PASS": 0, "FAIL": 1, "INDETERMINATE": 2}[verdict])
