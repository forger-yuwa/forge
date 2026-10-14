"""M6 (case/45 の H) の第一層の表 ÷ 見積もりの要求 の断面ごとの比を CSV と区間ごとの要約に書く (plan tooling-nozzle-card-local-mesh
§5.1 #3a、codex plan 段レビュー M3「比の分布を固定して commit」への対応)。

入力: step1_c45.py の出力 `c45/sections.npz` (見積もりの式で作った y1,req/r_w、Tw 370 K 一様、設計壁 + H の上流; AWS A `~/m55_mesh_est/c45/`) と
`c45/h_frac_table.json` (problem_h_hall.yaml の mesh.wall_first_frac_table。300 K 用に断熱 NS の y1+ × 冷却倍率 8.5〜9.4 で作った表)。
比の定義と区間は c45_hcheck.py と同じ (表の端は設計壁の範囲まで端値で延ばす、log 線形の補間)。
出力: `c45/h_over_est_ratio.csv` (断面 8001 点: x_rt, 見積もりの y1,req/r_w, H の表の値, 比) と `c45/h_over_est_ratio_summary.json` (区間ごとの最小・最大・中央と位置)。
使い方: python3 c45_ratio.py --dir <sections.npz と h_frac_table.json のあるディレクトリ> [--out c45]"""
import argparse
import json
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "c45"))
a = ap.parse_args()
D, O = Path(a.dir), Path(a.out)
sec = np.load(D / "sections.npz")
x, f = sec["x_rt"], sec["y1_over_rw"]
ft = np.asarray(json.load(open(D / "h_frac_table.json")), float)
ft[0, 0] = min(ft[0, 0], float(x[0]))
ft[-1, 0] = max(ft[-1, 0], float(x[-1]))
fh = np.exp(np.interp(x, ft[:, 0], np.log(ft[:, 1])))
with np.errstate(divide="ignore", invalid="ignore"):
    rat = fh / np.where(f > 0, f, np.nan)
with open(O / "h_over_est_ratio.csv", "w") as fh_:
    fh_.write("x_rt,f_req_est,f_H_table,ratio\n")
    for i in range(len(x)):
        fh_.write(f"{x[i]:.6f},{f[i]:.6e},{fh[i]:.6e},{rat[i]:.6f}\n")
x1 = float(x[-1])
segs = (("入口の配管", -12.5, -9.0), ("縮流部 x<-2", -9.0, -2.0), ("スロート付近 |x|<=2", -2.0, 2.0),
        ("拡大部 2<x<=x_e/2", 2.0, x1 / 2), ("拡大部 x>x_e/2", x1 / 2, x1))
out = {"source": "c45/sections.npz (step1_c45.py) と c45/h_frac_table.json; 比 = H の表 / 見積もりの要求 (c45_hcheck.py と同じ定義)",
       "Tw_K_estimate": float(np.median(sec["Tw_K"])), "segments": []}
for nm, lo, hi in segs:
    m = (x >= lo) & (x <= hi) & np.isfinite(rat)
    k_lo, k_hi = int(np.nanargmin(np.where(m, rat, np.nan))), int(np.nanargmax(np.where(m, rat, np.nan)))
    out["segments"].append({"seg": nm, "x_lo": lo, "x_hi": hi, "ratio_min": float(rat[k_lo]), "x_at_min": float(x[k_lo]),
                            "ratio_max": float(rat[k_hi]), "x_at_max": float(x[k_hi]), "ratio_median": float(np.nanmedian(rat[m]))})
    print(f"{nm:<20s} 比 {rat[k_lo]:.5f} (x {x[k_lo]:.3f}) – {rat[k_hi]:.5f} (x {x[k_hi]:.3f})  中央 {np.nanmedian(rat[m]):.4f}")
(O / "h_over_est_ratio_summary.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
