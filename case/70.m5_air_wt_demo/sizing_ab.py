"""case/70: 寸法の A/B (plan tooling-nozzle-deltastar-band-root-zone §5.1 #2、codex 諮問 2026-10-13 の推奨、forge 0 step)。

同じ問題 YAML・同じ物理スロート径・同じ壁の構築で、補正 δ_r だけを替えて r_t を解く:
  腕 A: run_0010 の壁と同じ補正 (CONTUR 積分法 k_f 0.8246; solve_rt_throat の CFD 前の経路)
  腕 B: run_0010 の修正後の抽出 (zone の腕) の δ_E 全分布 (solve_rt_throat --prev の経路)
腕 B は #18 (表の終点が x_e より float32 の丸め分だけ手前) を、計測のために終点の x を x_e に置き換えた写しで避ける
(置き換える前に差が ε_f32 = ulp_f32(x_e·S)/S + ulp_f32(x_e) 以内であることを確かめる)。

usage (case dir, design/.venv-opt):  python sizing_ab.py
"""
import json
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.probdef import load_problem  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.feedback.deltastar_loop import solve_rt_throat  # noqa: E402

PROB = C / "problem_m5_air_c2b.yaml"
R_THROAT = 0.025                          # カード: 物理スロート径 50 mm
CARD = {"pipe": 0.20, "contraction": 0.30, "divergent": 1.20}
TOL_L = 1.0e-3
NS, EU = C / "run_0010_ns_c2b_final", C / "run_0003_euler_36k"
ZONE_NEXT = C / "_ab_root_zone/run_0010/zone/delta_r_next.csv"
OUT = C / "_sizing_ab"
OUT.mkdir(exist_ok=True)

p = load_problem(PROB); d = design_chain(p)
x_e = float(d["wall"].x_e)
S_ns = float(json.loads((NS / "prepare_info.json").read_text())["scale_m"])
S_eu = float(json.loads((EU / "prepare_info.json").read_text())["scale_m"])

# 腕 B の表: 終点を x_e に揃えた写し
lines = ZONE_NEXT.read_text().splitlines()
hdr, rows = lines[0], [l for l in lines[1:] if l.strip()]
x_last = float(rows[-1].split(",")[0])
eps = float(np.spacing(np.float32(x_e * S_ns))) / S_ns + float(np.spacing(np.float32(x_e)))
gap = x_e - x_last
assert 0.0 < gap <= eps, f"表の終点の差 {gap!r} が ε_f32 {eps!r} の外 — 揃えない"
f = rows[-1].split(","); f[0] = repr(x_e); rows[-1] = ",".join(f)
snapped = OUT / "delta_r_next_snapped.csv"
snapped.write_text("\n".join([hdr + f" [sizing_ab: 終点 x {x_last!r} → x_e {x_e!r}、差 {gap:.3e} ≤ ε_f32 {eps:.3e}]"] + rows) + "\n")


out = {"problem": PROB.name, "R_throat_m": R_THROAT, "S_ns_run_0010_m": S_ns, "S_euler_run_0003_m": S_eu,
       "x_e_rt": x_e, "table_end_gap_rt": gap, "eps_f32_rt": eps}
arms = {"A": dict(), "B": dict(prev_run=NS, delta_r_out=OUT / "delta_r_B.csv", delta_next=snapped)}
for arm, kw in arms.items():
    r = solve_rt_throat(PROB, R_THROAT, **kw)
    rt = r["r_t_m"]
    L = {"pipe_m": float(p.geometry["L_pipe"]) * rt, "contraction_m": float(p.geometry["L_U"]) * rt, "divergent_m": r["length_m"]}
    err = {"pipe": L["pipe_m"] - CARD["pipe"], "contraction": L["contraction_m"] - CARD["contraction"], "divergent": L["divergent_m"] - CARD["divergent"]}
    out[arm] = {"r_t_m": rt, "rt_over_run0010": rt / S_ns, "residual_m": r["residual_m"], "n_iter": r["n_iter"],
                "physical_throat": r["physical_throat"], "exit_radius_m": r["exit_radius_m"], "lengths": L,
                "length_err_m": err, "within_tol": all(abs(v) <= TOL_L for v in err.values()),
                "area_ratio_vs_euler_scale": (rt / S_eu) ** 2, "delta_r_source": str(r.get("delta_r_source"))[:200]}
A, B = out["A"], out["B"]
out["A_reproduces_run0010"] = abs(A["r_t_m"] - S_ns) <= 1e-9
if not out["A_reproduces_run0010"]:
    verdict = "比較不成立 (A が run_0010 の r_t を再現しない)"
elif A["within_tol"] and not B["within_tol"]:
    verdict = "r_t 単位の設計を据え置く経路を棄却 (B が長さの公差外) → カードから実寸の長さを保って作り直す"
elif B["within_tol"]:
    verdict = "寸法の公差破綻の仮説を棄却 (B も公差内)"
else:
    verdict = "A も公差外 (判定条件の想定外)"
out["verdict"] = verdict
(OUT / "sizing_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
print(json.dumps({k: out[k] for k in ("A_reproduces_run0010", "verdict")}, ensure_ascii=False))
for arm in ("A", "B"):
    a = out[arm]
    print(arm, f"r_t {a['r_t_m']:.9f} (×{a['rt_over_run0010']:.6f})", "throat x/r [rt]", round(a["physical_throat"]["x_rt"], 5), round(a["physical_throat"]["r_rt"], 6),
          "len err [mm]", {k: round(v * 1e3, 3) for k, v in a["length_err_m"].items()}, "area/Euler", round(a["area_ratio_vs_euler_scale"], 6))
