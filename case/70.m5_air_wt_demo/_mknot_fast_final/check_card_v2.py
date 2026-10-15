#!/usr/bin/env python3
"""要求カードの生成の出力 (規則 coarse_min_theta_v2) を照合する (plan tooling-nozzle-mknot-fast §5.1 #10、result 段 NO-GO の M1)。
usage: python check_card_v2.py <design のパス> <problem.yaml> [gen.log]
見るもの: 規則の名前・自動、YAML の M_knot = 記録の選択、最終の 3 点が生産の網で成立、寸法の判定、YAML の要約と card_log の照合
(`wt_card.mknot_card_log_check`)、M_knot の所要の内訳とカード全体の時間 (と gen.log の壁時計)、実効の入力、並列。"""
import json
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, sys.argv[1])
from forge_design.intake import wt_card as WC  # noqa: E402

p = Path(sys.argv[2])
prob = yaml.safe_load(p.read_text())
cl = json.loads(p.with_name(p.name + ".card_log.json").read_text())
mk = cl["M_knot_rule"]
sel = mk["selection"]
fails = []


def check(name, cond, info=""):
    print(("ok   " if cond else "FAIL ") + name + (f"  [{info}]" if info else ""))
    if not cond:
        fails.append(name)


check("規則 coarse_min_theta_v2・自動", mk["rule"] == "coarse_min_theta_v2" and mk["mode"] == "auto", f"{mk['rule']} {mk['mode']}")
check("YAML の geometry.M_knot = 記録の選択", float(prob["geometry"]["M_knot"]) == float(mk["chosen"]), f"{prob['geometry']['M_knot']} / {mk['chosen']}")
f3 = mk["final_three"]
rows = [(r.get("M_K"), r.get("ok"), r.get("theta_max_deg"), r.get("reasons")) for r in f3]
check("最終の 3 点が生産の網で成立", len(f3) == 3 and all(r.get("ok") for r in f3), "; ".join(f"M_K {a} ok {b} θ_max {c}" for a, b, c, _ in rows))
dims = cl.get("dimensions") or {}
check("寸法の判定 PASS", bool(dims.get("pass")), "; ".join(f"{r['name']} {r['diff_m']:+.3e} m (公差 {r['tol_m']:g})" for r in dims.get("rows", [])))
c = WC.mknot_card_log_check(p)
check("YAML の要約と card_log の照合 (sha256・作り直した要約)", c.get("ok") is True, str({k: c.get(k) for k in ("ok", "summary_matches")}))
eb = mk.get("elapsed_breakdown", {})
wall = None
if len(sys.argv) > 3:
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", Path(sys.argv[3]).read_text())
    wall = m.group(1) if m else None
check("M_knot の所要 ≤ 600 s", float(mk["elapsed_s"]) <= 600.0,
      f"M_knot {mk['elapsed_s']:.1f} s (選択 {eb.get('select_s', float('nan')):.1f} s・r_t の固定点のやり直し {eb.get('refix_s')} s)、"
      f"記録の時点のカード全体 {mk.get('card_elapsed_s', float('nan')):.1f} s、壁時計 (gen.log) {wall}")
print("info 選択:", mk["reason"][:200])
print("info 実効の入力 (最終の組):", json.dumps(sel["inputs"], ensure_ascii=False))
print("info 選択の段:", [(s["stage"], s["chosen"], round(s["elapsed_s"], 1), s.get("discarded")) for s in mk["selections"]], "redo", mk.get("redo"),
      "n_restart", mk.get("n_restart"))
print("info 並列:", {k: mk["workers"].get(k) for k in ("cpu", "coarse", "prod", "per_eval_gb_coarse", "per_eval_gb_prod", "available_gb")})
print("info r_t:", mk["r_t_m"], "θ_max:", mk["theta_max_deg"], "余裕:", mk["margin_deg"], "最小 sin 角:", mk["min_sin_angle"])
print("VERDICT:", "PASS" if not fails else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
