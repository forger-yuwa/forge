"""§6.2 (1) の判定 (plan time_integration-line-implicit-speed、codex 2026-10-09 反映): FORGE_LINE_COMPARE=1 + FORGE_LINE_INV=1 の forge_run.log から
  - 全 solve の成分ごとの max|Δdq_c|/尺度_c (≤ 1e-8)、W の差 0・不一致 0、分解の失敗の不一致 0、非有限 0 件
  - 全ラインの後退誤差 ([lineEta] の行): 逆行列の η ≤ 1e-11 (評価した全ライン・全 solve)。従来の η も 1e-11 を超えた solve は判別不能として別に数える
を読む。usage: python3 inv_cmp_judge.py <run> [--ro-ref .. --a-ref ..] → 標準出力と <run>/inv_cmp_judge.json、不合格・判定不能なら終了コード 1。"""
import argparse, json, re, sys
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--ro-ref", type=float, default=0.8739869154); ap.add_argument("--a-ref", type=float, default=359.7712185)
ap.add_argument("--arm", default="INV", help="比べる腕の名前 (INV / F32c / F32cs)"); ap.add_argument("--eta-limit", type=float, default=1e-11)
ap.add_argument("--dq-limit", type=float, default=1e-8, help="成分ごとの max|Δdq|/尺度 の上限 (負なら判定に使わない)")
a = ap.parse_args()
run = Path(a.run); log = (run / "forge_run.log").read_text(errors="replace")
sc = [a.ro_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref, a.ro_ref * a.a_ref ** 2]
num = r"([0-9.eE+-]+|nan|inf)"
solves = re.findall(r"\[lineCompare\] solve (\d+): .*成分ごと 差/最大 = " + " ".join([num + "/" + num] * 5) + r" \(非有限の累計 (\d+)\)", log)
facts = re.findall(r"\[lineCompare\] factor (\d+).*?W " + num + " / " + num + r" \(不一致 (\d+)\)、ピボットの不一致 (\d+)、失敗の不一致 (\d+)", log)
etas = re.findall(r"\[lineEta\] solve (\d+) (\w+): η 最大 " + num + r" \(ライン (-?\d+)\)、評価 (\d+) 本、1e-11 超 (\d+) 本、非有限 (\d+) 本、分解の失敗 (\d+) 本", log)
nonfin_lines = re.findall(r"\[lineCompare\] 非有限 (\d+) 件", log)
if not solves or not facts or not etas:
    print(f"比較の行が足りない (solve {len(solves)}・factor {len(facts)}・η {len(etas)}) — 判定不能"); sys.exit(1)
worst_phys = [0.0] * 5; worst_rel = [0.0] * 5
for s in solves:
    v = [float(x) for x in s[1:11]]
    for c in range(5):
        d, m = v[2 * c], v[2 * c + 1]
        worst_phys[c] = max(worst_phys[c], d / sc[c]); worst_rel[c] = max(worst_rel[c], d / m if m > 0 else 0.0)
nonfinite = max(int(s[11]) for s in solves) + sum(int(x) for x in nonfin_lines) * 0   # 累計は最後の行が全体
wdiff = max(float(f[1]) for f in facts); wmis = sum(int(f[3]) for f in facts); fmis = sum(int(f[5]) for f in facts)
eta = {a.arm: [], "LU": []}
for e in etas:
    if e[1] not in eta: continue
    eta[e[1]].append({"solve": int(e[0]), "max": float(e[2]), "line": int(e[3]), "n_eval": int(e[4]), "n_over": int(e[5]), "n_nonfinite": int(e[6]), "n_fail": int(e[7])})
lu_over = {r["solve"] for r in eta["LU"] if r["n_over"] > 0 or r["n_nonfinite"] > 0}
inv_over = [r for r in eta[a.arm] if (r["max"] > a.eta_limit or r["n_nonfinite"] > 0)]
inv_over_judgeable = [r for r in inv_over if r["solve"] not in lu_over]
n_fail_arm = max(r["n_fail"] for r in eta[a.arm]) if eta[a.arm] else -1
ok = ((a.dq_limit < 0 or max(worst_phys) <= a.dq_limit) and (a.arm != "INV" or (wdiff == 0.0 and wmis == 0)) and fmis == 0 and nonfinite == 0
      and not inv_over_judgeable and bool(eta[a.arm]) and n_fail_arm == 0)
indeterminate = sorted(lu_over)
out = {"run": run.name, "n_solve": len(solves), "n_factor": len(facts), "max_dq_diff_over_scale": worst_phys, "max_dq_diff_over_max": worst_rel,
       "W_max_diff": wdiff, "W_mismatch": wmis, "fail_mismatch": fmis, "nonfinite_total": nonfinite,
       "arm": a.arm, "eta_limit": a.eta_limit, "eta_arm_max": max(r["max"] for r in eta[a.arm]), "eta_arm_worst": max(eta[a.arm], key=lambda r: r["max"]),
       "eta_LU_max": max(r["max"] for r in eta["LU"]), "eta_LU_worst": max(eta["LU"], key=lambda r: r["max"]),
       "n_fail_lines_arm_max": n_fail_arm, "n_fail_lines_LU_max": max(r["n_fail"] for r in eta["LU"]),
       "solves_indeterminate_LU_also_over": indeterminate, "inv_over_judgeable": inv_over_judgeable, "pass": ok}
(run / f"cmp_judge_{a.arm}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False)); sys.exit(0 if ok else 1)
