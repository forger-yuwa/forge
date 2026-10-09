"""ライン解の比較の判定 (plan time_integration-line-implicit-speed §6.4 (1)、codex 2026-10-09 plan-3 の反映)。FORGE_LINE_COMPARE=1 の forge_run.log から:
  - そろい: factor の番号がちょうど 1..F、solve の番号がちょうど 1..S、各 solve に腕と従来 (LU) の [lineEta] の行が 1 本ずつ、どちらも評価本数 = ライン数
  - 非有限: [lineNonfinite] (実際に使った因子・y)、[lineCompare] 非有限、[lineEta] の非有限 を腕・従来とも合計し、1 件でもあれば FAIL
  - 分解の失敗: 両腕とも 0 本、不一致 0
  - 後退誤差: 腕の η > 上限の solve で、従来の η ≤ 1e-11 なら FAIL、従来も 1e-11 超なら INDETERMINATE
  - (任意) 成分ごとの max|Δdq_c|/尺度_c ≤ --dq-limit
終了コード: 0 = PASS、1 = FAIL、2 = INDETERMINATE (行の欠け・重複を含む)。結果は <run>/cmp_judge_<腕>.json。
usage: python3 line_cmp_judge.py <run> --arm F32cs --eta-limit 1e-5 [--dq-limit -1] [--factors 20 --solves 100 --lines 4719]"""
import argparse, json, re, sys
from collections import Counter
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--arm", required=True); ap.add_argument("--eta-limit", type=float, required=True)
ap.add_argument("--dq-limit", type=float, default=-1.0); ap.add_argument("--factors", type=int, default=20); ap.add_argument("--solves", type=int, default=100)
ap.add_argument("--lines", type=int, default=4719); ap.add_argument("--ro-ref", type=float, default=0.8739869154); ap.add_argument("--a-ref", type=float, default=359.7712185)
a = ap.parse_args()
run = Path(a.run); log = (run / "forge_run.log").read_text(errors="replace")
num = r"([0-9.eE+-]+|nan|-?inf)"
reasons_fail, reasons_ind = [], []
fac = [int(m.group(1)) for m in re.finditer(r"\[lineCompare\] factor (\d+)", log)]
sol = re.findall(r"\[lineCompare\] solve (\d+): .*成分ごと 差/最大 = " + " ".join([num + "/" + num] * 5), log)
if Counter(fac) != Counter(range(1, a.factors + 1)): reasons_ind.append(f"factor の番号が 1..{a.factors} ちょうどでない ({len(fac)} 行)")
if Counter(int(s[0]) for s in sol) != Counter(range(1, a.solves + 1)): reasons_ind.append(f"solve の番号が 1..{a.solves} ちょうどでない ({len(sol)} 行)")
facts = re.findall(r"\[lineCompare\] factor (\d+).*?失敗の不一致 (\d+)", log)
fmis = sum(int(f[1]) for f in facts)
if fmis: reasons_fail.append(f"分解の失敗の不一致 {fmis}")
nf = sum(int(x) for x in re.findall(r"\[lineCompare\] 非有限 (\d+) 件", log))
for m in re.finditer(r"\[lineNonfinite\] (?:factor|solve) \d+: (.*)", log):
    nf += sum(int(x) for x in re.findall(r"(\d+)(?:・|、|$)", m.group(1)))
etas = re.findall(r"\[lineEta\] solve (\d+) (\w+): η 最大 " + num + r" \(ライン (-?\d+)\)、評価 (\d+) 本、1e-11 超 (\d+) 本、非有限 (\d+) 本、分解の失敗 (\d+) 本", log)
E = {a.arm: {}, "LU": {}}
for e in etas:
    if e[1] not in E: continue
    k = int(e[0])
    if k in E[e[1]]: reasons_ind.append(f"η の行の重複 (solve {k}・{e[1]})")
    E[e[1]][k] = {"max": float(e[2]), "line": int(e[3]), "n_eval": int(e[4]), "n_over": int(e[5]), "n_nonfinite": int(e[6]), "n_fail": int(e[7])}
for arm in (a.arm, "LU"):
    if set(E[arm]) != set(range(1, a.solves + 1)): reasons_ind.append(f"{arm} の η の行が solve 1..{a.solves} にそろわない ({len(E[arm])} 本)")
    for k, r in E[arm].items():
        nf += r["n_nonfinite"]
        if r["n_fail"]: reasons_fail.append(f"{arm} の分解の失敗 {r['n_fail']} 本 (solve {k})")
        elif r["n_eval"] != a.lines: reasons_ind.append(f"{arm} の評価本数 {r['n_eval']} ≠ {a.lines} (solve {k})")
if nf: reasons_fail.append(f"非有限 {nf} 件")
for k, r in E[a.arm].items():
    if r["max"] > a.eta_limit:
        lu = E["LU"].get(k)
        if lu is not None and lu["max"] > 1e-11: reasons_ind.append(f"solve {k}: 腕の η {r['max']:.2e} > 上限、従来も {lu['max']:.2e} > 1e-11")
        else: reasons_fail.append(f"solve {k}: 腕の η {r['max']:.2e} > {a.eta_limit:.0e}")
sc = [a.ro_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref ** 2]
wp = [0.0] * 5; wr = [0.0] * 5
for s in sol:
    v = [float(x) for x in s[1:11]]
    for c in range(5):
        wp[c] = max(wp[c], v[2 * c] / sc[c]); wr[c] = max(wr[c], v[2 * c] / v[2 * c + 1] if v[2 * c + 1] > 0 else 0.0)
if a.dq_limit >= 0 and max(wp) > a.dq_limit: reasons_fail.append(f"補正の差 {max(wp):.2e} > {a.dq_limit:.0e}")
verdict = "FAIL" if reasons_fail else ("INDETERMINATE" if reasons_ind else "PASS")
out = {"run": run.name, "arm": a.arm, "verdict": verdict, "fail": reasons_fail, "indeterminate": reasons_ind, "nonfinite": nf,
       "eta_arm_max": max((r["max"] for r in E[a.arm].values()), default=None), "eta_LU_max": max((r["max"] for r in E["LU"].values()), default=None),
       "eta_arm_worst": max(E[a.arm].items(), key=lambda kv: kv[1]["max"], default=None), "max_dq_diff_over_scale": wp, "max_dq_diff_over_max": wr,
       "n_factor_lines": len(fac), "n_solve_lines": len(sol)}
(run / f"cmp_judge_{a.arm}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
sys.exit({"PASS": 0, "FAIL": 1, "INDETERMINATE": 2}[verdict])
