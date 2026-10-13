"""plan tooling-nozzle-core-grid §4.104 P-hs の判定 (事前登録、ユーザ指示・codex 諮問 2026-10-14 hi-settle)。
腕 L = run_0917_hs_b0 (H・B0)、腕 P = run_0918_hs_point (H・point)。共通の起点は run_0916_hi_h_cont の res_40000 (本段の通算 10 万 step)。step は継続の中の番号。
判定の順序:
 1. 停止規則への初回到達 (stop_rule.judge_rows、719dfbf4、帯は起点で固定): L の s₁。N_reached = 100,000 + s₁ (本段の通算。起動段 4,000 step は別、
    移送元 C の計算は含まない)。探索の上限は本段の通算 20 万 (継続 10 万)。
 2. 確認と資格: L は s₁ の後の連続 4 万 step ([s − 40k, s]、s − 40k ≥ s₁) の全出力で停止規則を保ち、同じ窓で 8 量の資格 (§4.30 の式、cg_hs_lib.qualify) を満たす最初の s
    (= 見張りが止めた STOP_RULE の出力)。N_confirmed = 100,000 + s。P は [s − 40k, s] で 8 量の資格を満たし、s で停止規則を満たす (point 切替後の step で数える)。
    相対量の分母は |μ_H,point| (P の資格の窓の後半の平均)。残差の定常は求めない (全保存量の check_convergence は別に記録)。
 3. 方式の差: 8 量で |μ_L − μ_P| + U_L + U_P + p_diff ≤ τ (θ_r・δ_loc・Q_w 1.5 %、出口の主流 M 0.05 %、軸・η 0.1 の M 0.05 %pt)。p_diff = Gref の後処理の予算 × 2。
    形の差 (記録): μ_P − μ_C,point (C の point の参照 §4.31、grefx_refpoint_judge.json の ref_point) を |μ_C,point| で割る。U と p を添える。
結果: A = L が通算 20 万までに到達・確認に合格・P が資格を取得・8 量すべて予算 ≤ τ。
      B = 正常な計算・有効なデータで、L が通算 20 万までに未到達か確認が上限内に成立しない (B-予算不足)、または両腕が資格を得て予算超過 (B-許容超過)。
      判定不能 = P が上限までに資格を得ない・抽出や履歴の欠け・実行やデータの異常 (L の到達の実績は報告してよい)。発散は DIVERGED と書く。
参考 (判定に使わない): M6 の見張りの断面 x = 40・70・94 r_t の θ_r に同じ D の式と |欠損| ≤ 0.1 kg/s を当てた最初の到達。
usage (AWS の case/45 で): python3 cg_hs_judge.py → _band_ab/core_grid/hs_judge.json"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cg_hs_lib as L  # noqa: E402
RL, RP = "run_0917_hs_b0", "run_0918_hs_point"
CAP = {RL: 140000, RP: 200000}; REACH_MAX = 100000; OFF = 100000
P_ONE = {"theta": 0.010, "delta": 0.023, "Q_w": 0.009, "exit_M_line": 0.0006, "eta_pctpt": 0.0014}
out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.104", "runs": {"L": RL, "P": RP}, "offset_main_steps": OFF}
qual, info = {}, {}
for run in (RL, RP):
    d = HERE / run; i = {}
    rl = d / "restart.log"; i["restart_bitwise_9"] = rl.exists() and "VERDICT: OK (9 量を移した、SRC とビット一致)" in rl.read_text()
    lg = (HERE / "cg_hs.log").read_text() if (HERE / "cg_hs.log").exists() else ""
    i["prepared"] = f"{run}: 準備 (" in lg
    i["RUN_RC"] = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    i["STOP_RULE"] = int((d / "STOP_RULE").read_text()) if (d / "STOP_RULE").exists() else None
    ser = L.series(run); rows, jd = L.stop_judgements(run)
    phys = [json.loads(l) for l in open(L.OUTD / f"phys_{run}.jsonl")] if (L.OUTD / f"phys_{run}.jsonl").exists() else []
    steps = [r["step"] for r in ser]
    i["n_series"] = len(ser); i["last_step"] = steps[-1] if steps else None
    i["contiguous"] = steps == list(range(2500, (steps[-1] if steps else 0) + 1, 2500)) and [r["step"] for r in rows] == steps
    i["finite_physical"] = not any(r.get("nonfinite") for r in ser) and not any(r.get("nonfinite_all") for r in rows) \
        and all(p.get("min_ro", 1) > 0 and p.get("min_P", 1) > 0 and p.get("min_T", 1) > 0 for p in phys)
    i["diverged"] = any(r.get("nonfinite_all") for r in rows)
    i["first_reach"] = next((j["step"] for j in jd if j.get("ok")), None)
    i["normal_end"] = (i["RUN_RC"] == "0" and i["last_step"] == CAP[run]) or (i["STOP_RULE"] is not None)
    info[run] = i; info[run]["_ser"] = ser; info[run]["_okmap"] = {j["step"]: bool(j.get("ok")) for j in jd}; info[run]["_rows"] = rows
valid = all(info[r]["restart_bitwise_9"] and info[r]["prepared"] and info[r]["contiguous"] and info[r]["finite_physical"] and info[r]["normal_end"] for r in (RL, RP))
# P の資格 (分母は自分の μ で一度出し、その μ を全体の分母にする)
sP = info[RP]["STOP_RULE"]
qP = L.qualify(info[RP]["_ser"], sP) if sP else None
D0 = {k: v["mu"] for k, v in qP["q"].items()} if qP else None
if qP:
    qP = L.qualify(info[RP]["_ser"], sP, D0=D0)
P_ok = bool(qP and qP["ok"] and info[RP]["_okmap"].get(sP))
sL = info[RL]["STOP_RULE"]; frL = info[RL]["first_reach"]
qL = L.qualify(info[RL]["_ser"], sL, D0=D0) if (sL and D0) else None
L_reached = frL is not None and frL <= REACH_MAX
L_confirmed = bool(sL and L_reached and sL - L.CONF >= frL and all(info[RL]["_okmap"].get(k) for k in range(sL - L.CONF, sL + 1, 2500)) and qL and qL["ok"])
budget = {}
if qP and qL:
    for k in L.REQUIRED8:
        kd = L.kind(k); t = L.TAU[kd]
        dmu = abs(qL["q"][k]["mu"] - qP["q"][k]["mu"]); dmu = dmu if kd == "eta_pctpt" else 100 * dmu / abs(D0[k])
        B = dmu + qL["q"][k]["U"] + qP["q"][k]["U"] + 2 * P_ONE[kd]
        budget[k] = {"dmu": dmu, "U_L": qL["q"][k]["U"], "U_P": qP["q"][k]["U"], "p_diff": 2 * P_ONE[kd], "B": B, "tau": t, "ok": B <= t}
if any(info[r]["diverged"] for r in (RL, RP)):
    res = "判定不能 (DIVERGED: " + ", ".join(r for r in (RL, RP) if info[r]["diverged"]) + ")"
elif not valid:
    res = "判定不能 (実行・データの異常: " + json.dumps({r: {k: v for k, v in info[r].items() if not k.startswith("_")} for r in (RL, RP)}, ensure_ascii=False) + ")"
elif not P_ok:
    res = "判定不能 (P が上限までに参照の資格を得ない)"
elif not L_reached or not L_confirmed:
    res = "B-予算不足 (" + ("L が本段の通算 20 万 step までに停止規則に未到達" if not L_reached else "L の到達後の 4 万 step の確認が上限内に成立しない") + ")"
elif all(v["ok"] for v in budget.values()):
    res = "A: 今回の起点・条件で、B0 は停止規則に到達し、point の参照に対して 8 量とも許容内"
else:
    res = "B-許容超過 (両腕が資格を得たが、方式の差の予算が τ を超える量がある: " + ", ".join(k for k, v in budget.items() if not v["ok"]) + ")"
# 形の差 (記録)
ref = json.loads((L.OUTD / "grefx_refpoint_judge.json").read_text())["stats"]["ref_point"]["q"]
shape = {}
if qP:
    for k in L.REQUIRED8:
        kd = L.kind(k); muC = ref[k]["mu"]; muH = qP["q"][k]["mu"]
        dd = (muH - muC) if kd == "eta_pctpt" else 100 * (muH - muC) / abs(muC)
        shape[k] = {"mu_H_point": muH, "mu_C_point": muC, "diff": dd, "U_H_point_on_C_scale": qP["q"][k]["U"] if kd == "eta_pctpt" else qP["q"][k]["U"] * abs(D0[k]) / abs(muC),
                    "U_C_point": ref[k]["U"], "p": 2 * P_ONE[kd]}
# 参考: M6 の断面
def m9_first(run):
    rows = info[run]["_rows"]; S = np.array([r["step"] for r in rows], float)
    for idx, s in enumerate(S):
        m = (S >= s - 20000) & (S <= s)
        if S[0] > s - 20000 or m.sum() < 9:
            continue
        ok = abs(rows[idx].get("deficit_kg_s", 1e9)) <= 0.1
        for x in (40, 70, 94):
            v = np.array([r.get(f"m9_theta_r_{x}", np.nan) for r in rows])[m]; xs = S[m]; mid = s - 10000
            a, b = v[xs < mid], v[xs >= mid]
            Dm = abs(b.mean() - a.mean()) / abs(v.mean()) * 20000 / (xs[xs >= mid].mean() - xs[xs < mid].mean())
            ok = ok and Dm <= 5e-4
        if ok:
            return int(s)
    return None
conv = {}
for run in (RL, RP):
    o = subprocess.run([sys.executable, str(L.TOOLS / "check_convergence.py"), str(HERE / run)], capture_output=True, text=True)
    conv[run] = o.stdout.strip().splitlines()[-14:]
out.update(info={r: {k: v for k, v in info[r].items() if not k.startswith("_")} for r in (RL, RP)}, valid=valid,
           L={"first_reach_cont": frL, "N_reached_main": (OFF + frL) if frL else None, "stop_cont": sL, "N_confirmed_main": (OFF + sL) if sL else None,
              "reached_within_budget": L_reached, "confirmed": L_confirmed, "qual": qL},
           P={"stop_cont_after_switch": sP, "first_reach_after_switch": info[RP]["first_reach"], "qualified": P_ok, "qual": qP},
           budget=budget, result=res, shape_H_point_minus_C_point_pct=shape,
           m9_like_first_reach={r: m9_first(r) for r in (RL, RP)}, convergence=conv,
           note="N は H の本段の通算 (起点は C の B0 の場を移した後の起動段 4,000 step の後)。他のセッションの開始の場・壁温・形・ガスでの必要 step 数として使わない (§4.104)。")
(L.OUTD / "hs_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=str))
print("結果:", res)
print("L: 初回到達 (継続)", frL, "→ 本段の通算", out["L"]["N_reached_main"], "/ 確認", sL, "→", out["L"]["N_confirmed_main"])
print("P: 資格", P_ok, "停止", sP, "初回到達", info[RP]["first_reach"])
for k, v in budget.items(): print(f"  {k:18s} |Δμ| {v['dmu']:.4f} U_L {v['U_L']:.4f} U_P {v['U_P']:.4f} B {v['B']:.4f} τ {v['tau']} {'OK' if v['ok'] else 'NG'}")
for k, v in shape.items(): print(f"  形の差 {k:18s} H {v['mu_H_point']:.6g} C {v['mu_C_point']:.6g} 差 {v['diff']:+.3f}")
print("M6 の断面の参考の到達:", out["m9_like_first_reach"])
