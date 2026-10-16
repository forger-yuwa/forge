#!/usr/bin/env python3
r"""case/71 run_0002_ns の θ_r の測定帯の A/B (forge 0 step)。plan tooling-nozzle-card-local-mesh §5.1 #3b、
codex 諮問 `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md` の「判別 A/B」。

保持済みの保存場の列 (既定 100k〜183k)・同じ Euler の参照 (run_0001_euler の res_12000.h5、sha256 で固定)・同じ 3 断面
(停止規則の x_fracs × 設計壁の出口) で、変えるのは積分帯の外縁 y_b だけ:
  - A: `stop_rule.json` の固定帯 (本段の最初の出力 res_1000 から見張りが抽出した帯)
  - B: 後期の場 (既定 res_183000) から同じ抽出 (`stop_rule.fix_band`、同じ抽出の条件) で作った帯を全時点に固定
θ_r は `stop_rule.snapshot` (= `theta_r_sections`) で計算する (見張りと同じ関数)。

D_θ (2 万 step あたり) = |後半の平均 − 前半の平均| / |全体の平均| × 20,000 / (後半の平均 step − 前半の平均 step)。
前半・後半は窓の中点で分け、中点ちょうどは後半 (`stop_rule.judge_rows` と同じ割り方; D_raw と step の差は judge_rows の戻りを使う)。
主の窓 = 保持済みの全時点 (100k〜183k)。副 (判定に使わない): 10k 等間隔の 100k〜180k、最小二乗の傾き、腕 A の 1,000 step ごとの
見張りの系列 (stop_rule.json の行) の同じ窓。

判定 (呼び出し側が事前に固定、閾値は停止規則の θ_drift_max = 0.05 %): 断面ごとに
  - 両腕とも D_θ > 0.05 % → 「帯の選び方だけで未到達を説明できる」を棄却
  - A だけ超え B が閾値以下 → 測定帯への依存
  - それ以外 → そのまま記録 (判定は上位に渡す)
これは測定の感度の診断で、間引いた系列による正式な停止の判定ではない。

副産物: 停止規則の系列の CSV (`check_quasisteady.py --series-csv` 用; stop_rule.json の行 1,000〜182,000 + 腕 A で計算し直した
183,000 の行。保存に失敗して stop_rule.json に 183,000 の行が無いため)。

使い方 (AWS A、design/ を import できる python3):
  python3 theta_band_ab.py --design ~/forge-cardmesh/design --run <run_0002_ns> --euler-run <run_0001_euler> \
      --band-b-dir <出力先>/band_B_183000 --out-json <出力先>/theta_band_ab.json --series-csv <出力先>/stop_rule_series.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

PER = 20000          # D_θ の換算の幅 (停止規則の窓と同じ)
RULE = ("断面ごとに: 両腕とも D_θ > 0.05 % → 「帯の選び方だけで未到達を説明できる」を棄却; A だけ超え B が閾値以下 → 測定帯への依存; "
        "それ以外はそのまま記録 (呼び出し側が事前に固定、2026-10-15)")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def drift_per(rows: list, keys: list, SR, window: int) -> dict:
    """rows (step 昇順) の最後の時点で window の窓の D_raw と前半・後半の step の差を judge_rows から取り、2 万 step あたりに換算する。"""
    th = SR.thresholds()
    th.update(window=int(window), min_points=min(len(rows), th["min_points"]))
    fr = [float(k.split("_")[-1]) for k in keys]
    th["x_fracs"] = fr
    j = SR.judge_rows(rows, **th)
    out = {"window": int(window), "n_window": j.get("n_window"), "half_step_gap": j.get("half_step_gap"), "why_not_judged": j.get("why")
           if "theta" not in j else None, "sections": {}}
    if "theta" not in j:
        return out
    for k in keys:
        t = j["theta"][k]
        out["sections"][k] = {"mean": t["mean"], "D_raw": t["D_raw"], "D_per20k": t["D_raw"] * PER / j["half_step_gap"],
                              "amp_rel": t["amp"]}
    return out


def ols_per(rows: list, keys: list) -> dict:
    s = np.array([r["step"] for r in rows], float)
    out = {}
    for k in keys:
        v = np.array([r[k] for r in rows], float)
        slope = float(np.polyfit(s, v, 1)[0])
        out[k] = {"slope_per_step": slope, "rel_per20k": abs(slope) * PER / abs(float(v.mean()))}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--design", required=True, help="design/ (forge_design を import する)")
    ap.add_argument("--run", required=True)
    ap.add_argument("--euler-run", required=True)
    ap.add_argument("--euler-res", default="res_12000.h5")
    ap.add_argument("--late-step", type=int, default=183000)
    ap.add_argument("--steps", default="100000,110000,120000,130000,140000,150000,160000,170000,180000,181000,182000,183000")
    ap.add_argument("--band-b-dir", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--series-csv", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.design).resolve()))
    from forge_design.evaluate import stop_rule as SR

    t0 = time.time()
    run = Path(a.run).resolve(); eu = Path(a.euler_run).resolve()
    rec = json.loads((run / SR.RECORD).read_text())
    sr = rec["settings"]
    band_A = rec["band"]
    fr = [float(f) for f in sr["x_fracs"]]
    keys = [SR.theta_key(f) for f in fr]
    steps = [int(s) for s in a.steps.split(",")]
    # Euler の参照: 腕 A の帯は見張りが run dir だけで指定 (= その時の最新の出力) した。その時の最新が res_12000 だったことは
    # stop_rule_band/delta_r_extract_summary.json の massflow.res_euler に残っている。腕 B はそのファイルを sha256 で固定して使う
    eu_steps = sorted(int(m.group(1)) for f in eu.iterdir() if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)))
    summ_A = json.loads((run / SR.BAND_DIR / "delta_r_extract_summary.json").read_text())
    eu_spec = {"run": str(eu), "res": a.euler_res, "sha256": sha256(eu / a.euler_res)}
    # 腕 B の帯 (後期の場から、見張りと同じ関数・同じ抽出の条件)
    band_B = SR.fix_band(run, a.late_step, eu_spec, Path(a.band_b_dir), **SR.band_extract_kw(sr))
    summ_B = json.loads((Path(a.band_b_dir) / "delta_r_extract_summary.json").read_text())

    geo = SR.RunGeometry(run)
    rows_A, rows_B, res_sha = [], [], {}
    for n in steps:
        p = run / f"res_{n}.h5"
        V = SR.read_values(p)
        res_sha[str(n)] = sha256(p)
        rows_A.append(SR.snapshot(geo, n, band_A, fr, values=V))
        rows_B.append(SR.snapshot(geo, n, band_B, fr, values=V))
        del V
    # 腕 A の再計算が見張りの記録 (stop_rule.json の行) と一致するか
    rec_by = {r["step"]: r for r in rec["rows"]}
    repro = {}
    for r in rows_A:
        q = rec_by.get(r["step"])
        if q is None:
            repro[str(r["step"])] = "stop_rule.json に行が無い (保存の失敗で 182,000 までしか無い)"
            continue
        repro[str(r["step"])] = {k: (r[k] - q[k]) for k in keys + ["r", "mdot_in", "mdot_out", "exit_core_M"]}
    max_abs = max((abs(v) for d in repro.values() if isinstance(d, dict) for v in d.values()), default=None)

    W = steps[-1] - steps[0]
    main_A = drift_per(rows_A, keys, SR, W)
    main_B = drift_per(rows_B, keys, SR, W)
    even = [n for n in steps if n % 10000 == 0]
    ev_A = drift_per([r for r in rows_A if r["step"] in even], keys, SR, even[-1] - even[0])
    ev_B = drift_per([r for r in rows_B if r["step"] in even], keys, SR, even[-1] - even[0])
    watch_rows = [r for r in rec["rows"] if steps[0] <= r["step"] <= steps[-1]]
    watch_A = drift_per(watch_rows, keys, SR, watch_rows[-1]["step"] - watch_rows[0]["step"])

    thr = float(sr["theta_drift_max"])
    verdict = {}
    for k in keys:
        dA = main_A["sections"][k]["D_per20k"]; dB = main_B["sections"][k]["D_per20k"]
        if dA > thr and dB > thr:
            v = "両腕とも超過 → 帯の選び方だけで未到達を説明できる、を棄却"
        elif dA > thr and dB <= thr:
            v = "A だけ超過・B は閾値以下 → 測定帯への依存"
        else:
            v = f"どちらにも当たらない (A {'超過' if dA > thr else '以下'}、B {'超過' if dB > thr else '以下'})"
        verdict[k] = {"D_A_per20k": dA, "D_B_per20k": dB, "threshold": thr, "verdict": v}
    all_reject = all(v["verdict"].startswith("両腕とも") for v in verdict.values())
    all_band = all(v["verdict"].startswith("A だけ") for v in verdict.values())
    overall = ("棄却 (3 断面とも両腕で D_θ > 0.05 %): 帯の選び方だけでは未到達を説明できない" if all_reject else
               "測定帯への依存 (3 断面とも A だけ超過)" if all_band else "断面で結果が分かれた (断面ごとの欄を見る)")

    # 帯の比較 (断面の位置と使用域 [0, x_w])
    xA, yA = np.asarray(band_A["x"]), np.asarray(band_A["y_b"])
    xB, yB = np.asarray(band_B["x"]), np.asarray(band_B["y_b"])
    same_x = bool(len(xA) == len(xB) and np.array_equal(xA, xB))
    sec = {}
    for k in keys:
        x0 = rows_A[-1][f"x_{k}"]
        ya = float(np.interp(x0, xA, yA)); yb = float(np.interp(x0, xB, yB))
        sec[k] = {"x": x0, "y_b_A": ya, "y_b_B": yb, "ratio_B_over_A": yb / ya,
                  "late_edge_ystar": yb / 1.25, "late_edge_inside_A": bool(yb / 1.25 <= ya)}
    use = (xA >= 0.0) & (xA <= geo.x_wall_exit + 1e-9)
    ratio = np.interp(xA[use], xB, yB) / yA[use]
    band_cmp = {"same_x": same_x, "sections": sec,
                "use_domain": [0.0, geo.x_wall_exit],
                "ratio_B_over_A": {"min": float(ratio.min()), "p10": float(np.percentile(ratio, 10)), "median": float(np.median(ratio)),
                                   "p90": float(np.percentile(ratio, 90)), "max": float(ratio.max())},
                "frac_late_edge_outside_A": float(np.mean(np.interp(xA[use], xB, yB) / 1.25 > yA[use]))}

    out = {
        "what": "θ_r の測定帯の A/B (forge 0 step)",
        "plan": "plans/active/tooling-nozzle-card-local-mesh.md §5.1 #3b",
        "consult": "notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md (判別 A/B)",
        "rule_prefixed": RULE,
        "created": datetime.now().astimezone().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "run": str(run), "steps": steps, "res_sha256": res_sha,
        "euler_ref": {**eu_spec, "euler_res_steps_present": eu_steps, "band_A_res_euler": summ_A.get("massflow", {}).get("res_euler"),
                      "band_B_res_euler": summ_B.get("massflow", {}).get("res_euler")},
        "sections_requested": rec.get("sections"),
        "band_A": {"source": "stop_rule.json の band (見張りが本段の最初の出力から抽出)", "step": band_A["step"],
                   "res_sha256": band_A["res_sha256"], "euler_ref": band_A["euler_ref"], "extract": band_A["extract"]},
        "band_B": {"source": "fix_band (後期の場)", "step": band_B["step"], "res_sha256": band_B["res_sha256"],
                   "euler_ref": band_B["euler_ref"], "extract": band_B["extract"], "out_dir": band_B["out_dir"], "seconds": band_B["seconds"]},
        "band_compare": band_cmp,
        "arm_A_reproduces_watch_rows": {"max_abs_diff": max_abs, "per_step": repro},
        "rows": {"A": rows_A, "B": rows_B},
        "D_theta_main": {"note": "主の窓: 保持済みの全時点。判定はこれ", "A": main_A, "B": main_B},
        "D_theta_secondary": {"note": "判定に使わない", "even_10k_A": ev_A, "even_10k_B": ev_B,
                              "ols_A": ols_per(rows_A, keys), "ols_B": ols_per(rows_B, keys),
                              "watch_rows_1000_A": watch_A,
                              "watch_last_window_20k_A": rec.get("judge_last")},
        "verdict": {"per_section": verdict, "overall": overall},
        "seconds": round(time.time() - t0, 1),
    }
    Path(a.out_json).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=SR._json_default))

    # 停止規則の系列の CSV (check_quasisteady --series-csv 用)
    cols = ["step", "r", "mdot_in", "mdot_out"] + keys + ["exit_core_M", "mdot_rt", "source"]
    with open(a.series_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        last = 0
        for r in rec["rows"]:
            w.writerow([r["step"]] + [repr(float(r[c])) for c in cols[1:-1]] + ["stop_rule.json"])
            last = r["step"]
        for r in rows_A:
            if r["step"] > last:
                w.writerow([r["step"]] + [repr(float(r[c])) for c in cols[1:-1]] + ["recomputed_band_A"])
    print(json.dumps({"overall": overall, "per_section": verdict, "band_compare": band_cmp, "repro_max_abs": max_abs,
                      "seconds": out["seconds"]}, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
