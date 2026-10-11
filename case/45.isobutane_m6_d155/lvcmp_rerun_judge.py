#!/usr/bin/env python3
"""全部入り (L7) の再実行 1 本の判定 (plan time_integration-line-viscous-jacobian-dt §6.5、2026-10-11 事前登録、diagnostician 諮問を採用)。

L7r = run_0567_lvcmp_l7r を、§6.3 の判定器 (lvcmp_judge.py) と同じゲート・腕の照合 (見張りの記録のスキーマ・系列・記録係・到達の出力・区間付きの収束) で調べ、
L0 (run_0554)・L5 (run_0555)・L7 (run_0553) の記録と比べる。
  量: 到達の θ_r(40/70/94)・Q_w の L0 からの相対差 [%]。床 = 到達の窓 [N − 20000, N] の rms_ro の中央値 (residual_history.csv の outer_begin の行、L7r・L0・L5・L7 で同じ規則)。
  判定: 再現 = 4 量すべて + の向きで ≥ 0.3 % かつ 床の比 (L7r/L0) ≥ 2.0、かつ旧 L7 も同じ規則でパターンあり。再現せず = 4 量すべて |ずれ| ≤ 0.1 % かつ 床の比 < 1.5。
        言えるのは「同じ停止規則による到達の時点の差のパターンが再び出たか」まで (止まる状態の再現・統計的な同等性ではない。codex plan-3 M4)。
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
J.KIND["L7r"] = "L7"                          # L7r は L7 の構成 (codex plan-3 M1)
KEYS = J.KEYS
REP_DEV, NOREP_DEV, REP_FLOOR, NOREP_FLOOR = 0.3, 0.1, 2.0, 1.5


import hashlib  # noqa: E402
import math  # noqa: E402
RES_COLS = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega", "rms_roY0", "rms_roY1")


def floor_median(run, n_end):
    """到達の窓 [n_end − 20000, n_end] の outer_begin の行の残差の中央値 (列ごと)。窓の全 step が重複なくそろい、
    残差が有限・非負であることを要求する (codex plan-3 M2・m7)。戻り値: {"rms_ro": 中央値, "all": {列: 中央値}, "rows": 行数, "csv_sha256": ...}"""
    p = HERE / run / "residual_history.csv"
    want = set(range(n_end - 20000, n_end + 1))
    seen, vals = set(), {c: [] for c in RES_COLS}
    with open(p, newline="") as f:
        for r in csv.DictReader(f):
            if r["phase"] != "outer_begin":
                continue
            st = int(r["step"])
            if st not in want:
                continue
            if st in seen:
                raise ValueError(f"{run}: 床の窓の step {st} が重複")
            seen.add(st)
            for c in RES_COLS:
                v = float(r[c])
                if not (math.isfinite(v) and v >= 0):
                    raise ValueError(f"{run}: 床の窓の {c} (step {st}) が非有限か負")
                vals[c].append(v)
    if seen != want:
        raise ValueError(f"{run}: 床の窓の step がそろっていない ({len(seen)} / {len(want)})")
    med = {c: float(np.median(v)) for c, v in vals.items()}
    if not med["rms_ro"] > 0:
        raise ValueError(f"{run}: 床の中央値が正でない")
    return {"rms_ro": med["rms_ro"], "all": med, "rows": len(seen), "window": [n_end - 20000, n_end],
            "csv_sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


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
            ref = reach_row("L0")[1][-1]
            floors = {n: floor_median(J.ARMS[n][0], cls[n][1]["reach_step"]) for n in ("L7", "L0", "L5")}

            def devs(name):
                row = reach_row(name)[1][-1]
                d = {k: 100.0 * (float(row[k]) / float(ref[k]) - 1.0) for k in KEYS}
                if not all(math.isfinite(float(ref[k])) and float(ref[k]) > 0 for k in KEYS) or not all(math.isfinite(v) for v in d.values()):
                    raise ValueError(f"{name}: 到達の量か相対差が非有限 (codex plan-3 M3)")
                return d

            def classify(d, ratio):
                if not math.isfinite(ratio):
                    raise ValueError("床の比が非有限")
                vr, vn = all(x >= REP_DEV for x in d.values()), all(abs(x) <= NOREP_DEV for x in d.values())
                fr, fn = ratio >= REP_FLOOR, ratio < NOREP_FLOOR
                return vr, vn, fr, fn
            # 旧 L7 にも同じ規則を当てる (「2 本で出た」と言う前提、codex plan-3 M2)
            d7 = devs("L7"); r7 = floors["L7"]["rms_ro"] / floors["L0"]["rms_ro"]
            c7 = classify(d7, r7)
            rec["L7_under_same_rule"] = {"dev_vs_L0_pct": d7, "floor_ratio_vs_L0": r7, "value_pattern": c7[0], "floor_pattern": c7[2]}
            rec["floor_ratio_L5_vs_L0"] = floors["L5"]["rms_ro"] / floors["L0"]["rms_ro"]
            st = cls["L7r"][0]
            if st == "REACHED":
                nr = cls["L7r"][1]["reach_step"]
                dev = devs("L7r")
                floors["L7r"] = floor_median(J.ARMS["L7r"][0], nr)
                ratio = floors["L7r"]["rms_ro"] / floors["L0"]["rms_ro"]
                vr, vn, fr, fn = classify(dev, ratio)
                if vr and fr and c7[0] and c7[2]:
                    v = "再現: 同じ停止規則による到達の時点で、登録幅以上の正の差と高い残差の中央値が再び観測された (2 本、各 1 本の L0 に対して)"
                elif vn and fn:
                    v = "再現せず: 今回は、登録した差のパターン (値のずれ・高い床) を再現しなかった"
                else:
                    v = f"判別不能 (値: {'パターンあり' if vr else 'パターンなし' if vn else '中間'}、床: {'パターンあり' if fr else 'パターンなし' if fn else '中間'}、旧 L7 の値・床のパターン {c7[0]}・{c7[2]})"
                rec.update(dev_vs_L0_pct=dev, floor_ratio_vs_L0=ratio, value_pattern=vr, value_no_pattern=vn, floor_pattern=fr, floor_no_pattern=fn,
                           reach_step_L7r=nr, reach_step_diff_vs_L7=nr - cls["L7"][1]["reach_step"])
            elif st == "DIVERGED":
                v = "L7r が DIVERGED (全部入りは上限 50 で破綻しうる。値のパターンの再現は判別不能)"
            elif st == "CENSORED":
                floors["L7r"] = floor_median(J.ARMS["L7r"][0], J.BUDGET)
                rec["floor_ratio_vs_L0"] = floors["L7r"]["rms_ro"] / floors["L0"]["rms_ro"]
                v = "L7r が CENSORED (15 万 step で水準に入らない。値は比べない、床は最後の 2 万 step の窓で記録)"
            else:
                v = f"判別不能 (L7r の状態 {st})"
            rec["floors"] = floors
            rec["quasisteady"] = {n: (cls[n][1].get("quasisteady") or {}).get("overall") for n in cls}
            rec["VERDICT"] = v
        except Exception as e:
            rec["VERDICT"] = f"INVALID (判定で例外 {type(e).__name__}: {e})"
            invalid = ["judge"]
    try:
        out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str, allow_nan=False))
    except ValueError:
        rec = {"plan": rec.get("plan"), "VERDICT": "INVALID (記録に非有限の値)"}; invalid = ["nan"]
        out.write_text(json.dumps(rec, ensure_ascii=False, indent=1))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for n, c in cls.items():
        print(n, c[0], {k: c[1].get(k) for k in ("run", "reach_step", "segment_verdict", "why")})
    for k in ("dev_vs_L0_pct", "floor_ratio_vs_L0", "L7_under_same_rule", "floor_ratio_L5_vs_L0", "reach_step_diff_vs_L7", "quasisteady"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", rec["VERDICT"])
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
