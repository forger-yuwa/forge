#!/usr/bin/env python3
"""V-h370 の Euler 参照の擬似 CFL の A/B の判定 (plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36j・§6、
2026-10-14 事前登録 b1470a4a、codex diagnose notes/reviews/2026-10-14-h370-euler-ref-plateau-diagnose.md)。

AWS B の case dir で h370_euler_cfl_ab.sh が最後に呼ぶ。出力 `_band_ab/half_inlet/h370_euler_cfl_ab.json`。
- 収束: 元の run_0950 の mid → main と各腕の継続を連結した履歴 (step を 0 / 3000 / 63000 だけずらす) に
  `check_convergence.py --drop 3` (段の設定の差は cfl・cfl_pseudo と nStepOuter だけ = 連結可、h370_euler_cfl_ab.sh が確かめる)。
- 残差の統計: 継続の前半 [0, 20000)・後半 [20000, 40000) の outer_end の行で、列ごとの中央値・最大/最小・log10 の傾き (dec/1e4)。
  対照として run_0950 の本段の [40000, 60000) も同じ量を出す。
- 量: 出力 1,000 ごとの 8 量 (h370_judge.euler_series と同じ定義) の `check_quasisteady.py` の VERDICT と、末尾 40 % の窓の生の数値
  (平均・(最大−最小)/|平均|・(後半の平均 − 前半の平均)/|平均|; 表示の桁で丸めない)。
- 分類 (登録の文言を機械的に当てたもの。解釈の確定は上位に諮る): A が停滞 (NOT CONVERGED)・B が全残差 PASS → "B_PASS_A_STALL"、
  両方 NOT CONVERGED → "BOTH_STALL"、非有限・restart のビット一致の不成立・A が PASS → "INDETERMINATE"。
"""
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WT = Path.home() / "forge-g0cchain"
sys.path.insert(0, str(WT / "design"))
sys.path.insert(0, str(HERE))
import h370_judge as HJ  # noqa: E402  (euler_series・run_qs・QUAL_COLS を同じ定義で使う)

TOOLS = WT / "solver_density_cuda" / "tools"
SRC = HERE / "run_0950_h370_euler_ref"
ARMS = {"A": ("run_0952_h370_euler_cfl2", 2.0), "B": ("run_0953_h370_euler_cfl05", 0.5)}
OUT = HERE / "_band_ab" / "half_inlet" / "h370_euler_cfl_ab.json"
OFFSETS = (("residual_history_mid.csv", 0), ("residual_history_main.csv", 3000))
ARM_OFFSET = 63000
HALVES = ((0, 20000), (20000, 40000))


def read_csv(p: Path):
    with open(p) as fh:
        rows = list(csv.reader(fh))
    return rows[0], rows[1:]


def outer(p: Path):
    h, rows = read_csv(p)
    rows = [r for r in rows if len(r) > 2 and r[2] == "outer_end"]
    step = np.array([float(r[0]) for r in rows])
    cols = {c: np.array([float(r[i]) for r in rows]) for i, c in enumerate(h) if c.startswith("rms_")}
    return step, cols


def stats(step, cols, a, b):
    m = (step >= a) & (step < b)
    out = {}
    for c, v in cols.items():
        x = v[m]
        if x.size < 3 or not np.all(np.isfinite(x)):
            out[c] = {"n": int(x.size), "nonfinite": int(np.count_nonzero(~np.isfinite(x)))}
            continue
        if np.all(x == 0):
            out[c] = {"n": int(x.size), "all_zero": True}
            continue
        y = np.log10(np.maximum(x, 1e-300))
        out[c] = {"n": int(x.size), "median": float(np.median(x)), "min": float(x.min()), "max": float(x.max()),
                  "max_over_min": float(x.max() / x.min()) if x.min() > 0 else None,
                  "slope_dec_per_1e4": float(np.polyfit(step[m], y, 1)[0] * 1e4)}
    return out


def concat_convergence(arm_dir: Path) -> dict:
    hdr, out = None, []
    parts = [(SRC / f, off) for f, off in OFFSETS] + [(arm_dir / "residual_history.csv", ARM_OFFSET)]
    for p, off in parts:
        h, rows = read_csv(p)
        if hdr is None:
            hdr = h
        elif h != hdr:
            return {"error": f"列が違う: {p}"}
        for r in rows:
            r = list(r)
            r[0] = str(int(float(r[0])) + off)
            out.append(r)
    with tempfile.TemporaryDirectory() as td:
        with open(Path(td) / "residual_history.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(hdr)
            w.writerows(out)
        r = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), td, "--drop", "3"], capture_output=True, text=True)
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    verdict = next((ln for ln in lines if ln.startswith("===")), None)
    return {"cmd": "check_convergence.py <mid+main+arm> --drop 3", "rc": r.returncode, "verdict": verdict,
            "columns": [ln.strip() for ln in lines if ln.strip().startswith("rms_")],
            "overall": next((ln for ln in lines if ln.startswith("OVERALL")), None),
            "pass": r.returncode == 0, "rows": len(out)}


def raw_tail(rows, cols, tail=0.4):
    n = len(rows)
    k0 = int(math.floor(n * (1 - tail)))
    out = {}
    for c in cols:
        v = np.array([float(r.get(c, float("nan"))) for r in rows[k0:]])
        if v.size < 2 or not np.all(np.isfinite(v)):
            out[c] = {"n": int(v.size), "nonfinite": True}
            continue
        mu = float(np.mean(v))
        h = v.size // 2
        out[c] = {"n": int(v.size), "mean": mu, "range_rel": float((v.max() - v.min()) / abs(mu)) if mu else None,
                  "drift_rel": float((np.mean(v[h:]) - np.mean(v[:h])) / abs(mu)) if mu else None}
    return out


def arm_report(key: str) -> dict:
    name, cfl = ARMS[key]
    d = HERE / name
    rep = {"run": name, "cfl_pseudo": cfl}
    rlog = (d / "restart_field.log").read_text() if (d / "restart_field.log").exists() else ""
    rep["restart_verdict"] = [ln.strip() for ln in rlog.splitlines() if "VERDICT" in ln.upper()][-2:]
    rep["restart_bit_identical"] = any(ln.startswith("VERDICT: OK") and "ビット一致" in ln and "丸めあり" not in ln for ln in rep["restart_verdict"])
    rep["config_diff"] = (d / "config_diff.txt").read_text() if (d / "config_diff.txt").exists() else None
    step, cols = outer(d / "residual_history.csv")
    rep["last_step"] = float(step.max()) if step.size else None
    rep["residual_nonfinite"] = int(sum(np.count_nonzero(~np.isfinite(v)) for v in cols.values()))
    rep["residual_halves"] = {f"[{a},{b})": stats(step, cols, a, b) for a, b in HALVES}
    rep["convergence_concat"] = concat_convergence(d)
    steps = HJ.res_steps(d)
    steps = [s for s in steps if s > 0]
    rows, probs = HJ.euler_series(d, steps)
    rep["outputs"] = {"n": len(steps), "first": steps[0] if steps else None, "last": steps[-1] if steps else None}
    rep["field_problems"] = probs
    sc = d / "h370_euler_cfl_series.csv"
    HJ.write_series(sc, rows, HJ.QUAL_COLS + HJ.REC_COLS)
    rep["quasisteady"] = HJ.run_qs(TOOLS, sc, HJ.QUAL_COLS)
    rep["raw_tail"] = raw_tail(rows, HJ.QUAL_COLS)
    return rep


def main() -> int:
    res = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36j・§6 (事前登録 b1470a4a)",
           "source": {"run": SRC.name, "res": "res_60000.h5", "sha256": HJ.sha256_file(SRC / "res_60000.h5")}}
    s0, c0 = outer(SRC / "residual_history_main.csv")
    res["run_0950_main_tail"] = stats(s0, c0, 40000, 60000)
    res["arms"] = {}
    for k in ARMS:
        try:
            res["arms"][k] = arm_report(k)
        except Exception as e:  # noqa: BLE001
            res["arms"][k] = {"error": repr(e)}
    A, B = res["arms"].get("A", {}), res["arms"].get("B", {})
    bad = [k for k, r in res["arms"].items() if "error" in r or r.get("residual_nonfinite") or r.get("field_problems")
           or not r.get("restart_bit_identical")]
    pa = A.get("convergence_concat", {}).get("pass")
    pb = B.get("convergence_concat", {}).get("pass")
    if bad:
        cls, why = "INDETERMINATE", f"腕 {bad} に誤り・非有限・場の問題・restart のビット一致の不成立"
    elif pa:
        cls, why = "INDETERMINATE", "対照 A が全残差 PASS (元の停滞を再現しない)"
    elif pb:
        cls, why = "B_PASS_A_STALL", "A が停滞・B が全残差 PASS (CFL に依る反復の問題を支持。リミッタが原因の確定でも元の res_60000 の適格化でもない)"
    else:
        cls, why = "BOTH_STALL", "A・B とも不合格 (CFL を下げるだけで解ける、を退ける候補。B がまだ減衰中なら判別不能 — 後半の傾きで確かめる)"
    res["classification"] = {"mechanical": cls, "reason": why,
                             "note": "登録の文言を機械的に当てたもの。B の後半の傾き・A の床と run_0950 の末尾の比較を見て、解釈の確定は上位に諮る (条件 7)"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
    print(f"== Euler の CFL の A/B: {cls} — {why}")
    for k, r in res["arms"].items():
        cc = r.get("convergence_concat", {})
        print(f"  腕 {k} {r.get('run')}: {cc.get('verdict')}")
        for ln in cc.get("columns", []):
            print(f"     {ln}")
        qs = r.get("quasisteady", {}).get("verdicts", {})
        print(f"     量: " + ", ".join(f"{c}={v['verdict']}" for c, v in qs.items()))
    print(f"→ {OUT}")
    return 0 if cls != "INDETERMINATE" else 2


if __name__ == "__main__":
    sys.exit(main())
