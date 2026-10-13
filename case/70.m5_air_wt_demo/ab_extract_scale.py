"""case/70: 抽出の尺度の A/B (plan tooling-nozzle-wt-deliverables §5.1 #19、codex 諮問 2026-10-13 m5-demo-throat-delta)。

固定した NS 場と Euler 場を同じ抽出 (`extract_and_merge(.., band_select="edge")`) に 2 通りで通す (forge は回さない):
  腕 A: 現行 (各 run を自分の r_t = S で無次元化したまま照合)
  腕 B: Euler の x・r を S_E/S_NS 倍して NS の尺度に揃える (ρU_x は実寸なのでそのまま)
主判定量 D = max_{0 ≤ x/r_t ≤ 3} |δ_E| (delta_r_next.csv の delta_E 列)。判定は plan §5.1 #19 の事前登録どおり。

usage (case dir, design/.venv-opt):  python ab_extract_scale.py NS_RUN EULER_RUN
"""
import json
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next  # noqa: E402

ns, eu = C / sys.argv[1], C / sys.argv[2]
S_ns = float(json.loads((ns / "prepare_info.json").read_text())["scale_m"])
S_eu = float(json.loads((eu / "prepare_info.json").read_text())["scale_m"])
out = C / "_ab_extract_scale"
res = {"ns_run": ns.name, "euler_run": eu.name, "S_ns_m": S_ns, "S_euler_m": S_eu, "scale_ratio": S_eu / S_ns}
orig_load = DS._load_structured


def scaled_load(run_dir):
    d = orig_load(run_dir)
    if Path(run_dir).resolve() == eu.resolve():
        f = S_eu / S_ns                                   # Euler の r_t 単位 → NS の r_t 単位
        d = dict(d, x=d["x"] * f, r=d["r"] * f)
    return d


for arm, loader in (("A", orig_load), ("B", scaled_load)):
    DS._load_structured = loader
    od = out / arm
    od.mkdir(parents=True, exist_ok=True)
    summ = extract_and_merge(ns, eu, band_select="edge", out_dir=od)
    nx = read_delta_r_next(od / "delta_r_next.csv")
    x, dE = np.asarray(nx["x_rt"]), np.asarray(nx["delta_E"])
    m = (x >= 0.0) & (x <= 3.0)
    hdr = (od / "delta_r_next.csv").read_text().splitlines()[0]
    resid = hdr.split("resid_rel_rms=")[1].split(")")[0] if "resid_rel_rms=" in hdr else None
    res[arm] = {"D_max_abs_deltaE_0_3_rt": float(np.max(np.abs(dE[m]))),
                "deltaE_at": {str(xs): float(np.interp(xs, x, dE)) for xs in (0.0, 0.5, 1.0, 2.0, 3.0, 10.0, 30.0, 48.0)},
                "resid_rel_rms": resid, "delta_r_throat_use": summ.get("delta_r_throat_use"),
                "massflow": {k: v for k, v in (summ.get("massflow") or {}).items() if not str(k).startswith("res_")}}
DS._load_structured = orig_load
DA, DB = res["A"]["D_max_abs_deltaE_0_3_rt"], res["B"]["D_max_abs_deltaE_0_3_rt"]
if DB <= 0.1 * DA and DB <= 0.01:
    verdict = "H1 支持 (尺度の混在が主因; 抽出の物理的妥当性の合格ではない)"
elif DB >= 0.9 * DA:
    verdict = "H1 棄却 (符号つきの帯選択を優先)"
else:
    verdict = "判別未了 (中間)"
res["verdict"] = verdict
res["criteria"] = "plan §5.1 #19: B の D ≤ 0.1·A の D かつ ≤ 0.01 r_t → H1 支持; B の D ≥ 0.9·A の D → H1 棄却; 他は判別未了"
(out / "ab_extract_scale.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float))
print(json.dumps({"D_A": DA, "D_B": DB, "verdict": verdict, "A": res["A"]["deltaE_at"], "B": res["B"]["deltaE_at"],
                  "resid": (res["A"]["resid_rel_rms"], res["B"]["resid_rel_rms"])}, indent=1, ensure_ascii=False))
