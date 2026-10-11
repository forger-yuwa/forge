#!/usr/bin/env python3
"""全部入り (L7) の再実行 1 本の判定 (plan time_integration-line-viscous-jacobian-dt §6.5、2026-10-11 事前登録、diagnostician 諮問を採用)。

L7r = run_0567_lvcmp_l7r を、§6.3 の判定器 (lvcmp_judge.py) と同じゲート・腕の照合 (見張りの記録のスキーマ・系列・記録係・到達の出力・区間付きの収束) で調べ、
L0 (run_0554)・L5 (run_0555)・L7 (run_0553) の記録と比べる。
  量: 到達の θ_r(40/70/94)・Q_w の L0 からの相対差 [%]。床 = 到達の窓 [N − 20000, N] の rms_ro の中央値 (residual_history.csv の outer_begin の行、L7r・L0・L5・L7 で同じ規則)。
  判定: 再現 = 4 量すべて + の向きで ≥ 0.3 % かつ 床の比 (L7r/L0) ≥ 2.0。再現せず = 4 量すべて |ずれ| ≤ 0.1 % かつ 床の比 < 1.5。
        L7r が DIVERGED・CENSORED は別の判定 (CENSORED は値を比べず、床は最後の 2 万 step の窓)。それ以外は判別不能 (値と床を別々に記録)。閾値は事後に変えない。
  I の追加: L7r の設定が run_0553 と全項目同じ。
出力は _band_ab/cold_pair/lvcmp_rerun_judge.json。ゲートの不合格か INVALID があれば判定せず、終了コード 1。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.argv = sys.argv[:1]
import lvcmp_judge as J  # noqa: E402

J.ARMS["L7r"] = ("run_0567_lvcmp_l7r", "7", "m7")
KEYS = J.KEYS
REP_DEV, NOREP_DEV, REP_FLOOR, NOREP_FLOOR = 0.3, 0.1, 2.0, 1.5


def floor_median(run, n_end):
    """到達の窓 [n_end − 20000, n_end] の rms_ro の中央値 (outer_begin の行)"""
    with open(HERE / run / "residual_history.csv", newline="") as f:
        rows = [r for r in csv.DictReader(f) if r[list(r.keys())[2]] == "outer_begin"]
    c0 = list(rows[0].keys())[0]
    vals = [float(r["rms_ro"]) for r in rows if n_end - 20000 <= int(r[c0]) <= n_end]
    if len(vals) < 1000 or not all(np.isfinite(vals)):
        raise ValueError(f"{run}: 床の窓の行が足りないか非有限 ({len(vals)} 行)")
    return float(np.median(vals))


def reach_row(name):
    w = json.loads((HERE / J.ARMS[name][0] / "m9_watch.json").read_text())
    rows = sorted([x for x in w["rows"] if x["step"] > 0], key=lambda x: x["step"])
    return w, rows


def main():
    out = J.D / "lvcmp_rerun_judge.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6.5"}
    try:
        gates = J.input_gates()
    except Exception as e:
        gates = [("入力のゲートを評価できない", False, f"{type(e).__name__}: {e}")]
    try:
        c7, cr = J.cfg(J.ARMS["L7"][0]), J.cfg(J.ARMS["L7r"][0])
        dd = sorted(k for k in set(c7) | set(cr) if c7.get(k) != cr.get(k))
        gates.append(("L7r の設定が run_0553 と全項目同じ", dd == [], dd[:10]))
    except Exception as e:
        gates.append(("L7r の設定を照合できない", False, f"{type(e).__name__}: {e}"))
    cls = {}
    for n in ("L7", "L0", "L5", "L7r"):
        try:
            cls[n] = J.assess(n)
        except Exception as e:
            cls[n] = ("INVALID", {"why": f"例外 {type(e).__name__}: {e}"})
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = all(g[1] for g in gates)
    rec["arms"] = {n: {"class": c[0], **c[1]} for n, c in cls.items()}
    invalid = [n for n in cls if cls[n][0] == "INVALID"] + [n for n in ("L7", "L0", "L5") if cls[n][0] != "REACHED"]
    if not rec["gates_ok"] or invalid:
        rec["VERDICT"] = f"INVALID (ゲート不合格 {sum(not g[1] for g in gates)} 件、INVALID・未到達の基準の腕 {invalid}): 判定しない"
    else:
        try:
            n0 = cls["L0"][1]["reach_step"]
            ref = reach_row("L0")[1][-1]
            floors = {n: floor_median(J.ARMS[n][0], cls[n][1]["reach_step"]) for n in ("L7", "L0", "L5")}
            st = cls["L7r"][0]
            if st == "REACHED":
                nr = cls["L7r"][1]["reach_step"]
                row = reach_row("L7r")[1][-1]
                dev = {k: 100.0 * (row[k] / ref[k] - 1.0) for k in KEYS}
                floors["L7r"] = floor_median(J.ARMS["L7r"][0], nr)
                ratio = floors["L7r"] / floors["L0"]
                val_rep = all(d >= REP_DEV for d in dev.values())
                val_norep = all(abs(d) <= NOREP_DEV for d in dev.values())
                fl_rep, fl_norep = ratio >= REP_FLOOR, ratio < NOREP_FLOOR
                if val_rep and fl_rep:
                    v = "再現 (値のずれと高い床が再実行でも出た)"
                elif val_norep and fl_norep:
                    v = "再現せず (値のずれ・床とも生産と区別できない)"
                else:
                    v = "判別不能"
                rec.update(dev_vs_L0_pct=dev, floor_ratio_vs_L0=ratio, value_reproduced=val_rep, value_not_reproduced=val_norep,
                           floor_reproduced=fl_rep, floor_not_reproduced=fl_norep, reach_step_L7r=nr, reach_step_diff_vs_L7=nr - cls["L7"][1]["reach_step"])
                if v == "判別不能":
                    v += f" (値: {'再現' if val_rep else '再現せず' if val_norep else '中間'}、床: {'再現' if fl_rep else '再現せず' if fl_norep else '中間'})"
            elif st == "DIVERGED":
                v = "L7r が DIVERGED (全部入りは上限 50 で破綻しうる。値の再現性は判別不能)"
            elif st == "CENSORED":
                floors["L7r"] = floor_median(J.ARMS["L7r"][0], J.BUDGET)
                rec["floor_ratio_vs_L0"] = floors["L7r"] / floors["L0"]
                v = "L7r が CENSORED (15 万 step で水準に入らない。値は比べない、床は最後の 2 万 step の窓で記録)"
            else:
                v = f"判別不能 (L7r の状態 {st})"
            rec["floors_median_rms_ro"] = floors
            rec["VERDICT"] = v
        except Exception as e:
            rec["VERDICT"] = f"INVALID (判定で例外 {type(e).__name__}: {e})"
            invalid = ["judge"]
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for n, c in cls.items():
        print(n, c[0], {k: c[1].get(k) for k in ("run", "reach_step", "segment_verdict", "why")})
    for k in ("dev_vs_L0_pct", "floor_ratio_vs_L0", "floors_median_rms_ro", "reach_step_diff_vs_L7"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", rec["VERDICT"])
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
