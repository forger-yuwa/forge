"""§6.18 の総時間の判定 (plan time_integration-line-implicit-speed、diagnostician 2026-10-10 の推奨)。
usage (AWS の case dir):
  python3 tt_judge.py budgets  → 各腕の単価と損益分岐の予算 (step、2500 の倍数に切り上げ) を _band_ab/cold_pair/tt_budgets.json と標準出力 ("S2 162500" の行) に
  python3 tt_judge.py final    → _band_ab/cold_pair/tt_judge.json
単価: B0・S2・M64 は §6.16 のふるい (scr_check.json) と run_0376_ttu_* の B0、M64S2 は run_0376_ttu_* (B0 で挟んだ 3 本)。確認 (scr_check.json の ok) が欠けた run があれば判定不能 (終了コード 2)。
B0 の総時間 = 125000 step (系列を 2 出力連続で数え直した到達、§6.18) × B0 の単価 + 25 回の出力 + 起動 (m9_unit.json の L0 の値)。
腕の総時間 = 到達 step × 単価 + (到達/2500) 回の出力 + 起動。分解能 = 出力の間隔 × 単価 + step 数 × 単価の 3 本の幅 (腕と B0 の和)。
R1: 総時間が B0 − 分解能 より小さい → 速い、B0 + 分解能 より大きい・予算で打ち切り → 遅い、ほか → 判別不能。発散・異常 → 比較不可。
R2: 速い腕のうち総時間が最小のものを採用の候補にする (R3 のフラグがあれば候補にしない)。
R3: 到達時の状態を B0 と比べ、|Δθ_r| > 0.05 % (3 断面のどれか)・|ΔQ_w| > 0.1 %・Σ|res_ro|·2π が B0 の 2 倍超・末尾 2 万 step の rms_* の中央値が B0 の 2 倍超、のどれかでフラグ。
    B0 の θ_r・Q_w は系列の 125000、Σ|res_ro| は残っている run_0252 の res_60000 (通算 135000、到達の出力は削除済み)、末尾の残差は run_0252 の局所 30000〜50000 (通算 105000〜125000)。
R4: 単独の腕 (S2・M64) の到達 step が B0 の 0.9〜1.1 倍なら、step 数への効果は断定しない (記録)。"""
import csv, collections, json, math, statistics, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
NI, NJ, S = 4719, 121, 0.07667
B0_STEPS, B0_OUT_INT, ARM_OUT_INT = 125000, 5000, 2500
def bad(msg): print(msg + " — 判定不能"); sys.exit(2)
def ms(run):
    p = HERE / run / "scr_check.json"
    if not p.exists(): bad(f"{run}: scr_check.json が無い")
    j = json.loads(p.read_text())
    if not j.get("ok"): bad(f"{run}: 確認が ok でない ({j.get('why')})")
    return j["ms"]
def units():
    scr = [p.name for p in HERE.glob("run_037[345]_scr_r*")]
    ttu = [p.name for p in HERE.glob("run_0376_ttu_*")]
    if len(scr) != 33 or len(ttu) != 7: bad(f"単価の run の数が違う (ふるい {len(scr)}・ttu {len(ttu)})")
    u = {"B0": [ms(r) for r in scr + ttu if r.endswith("_B0")], "S2": [ms(r) for r in scr if r.endswith("_S2")],
         "M64": [ms(r) for r in scr if r.endswith("_M64")], "M64S2": [ms(r) for r in ttu if r.endswith("_M64S2")]}
    if [len(u[k]) for k in ("S2", "M64", "M64S2")] != [3, 3, 3]: bad("案の単価が 3 本ずつそろわない")
    return {k: {"ms": v, "median": statistics.median(v), "spread": max(v) - min(v)} for k, v in u.items()}
U9 = json.loads((D / "m9_unit.json").read_text())["L0"]
OUT_S, ST_S = U9["median_output_ms"] / 1000, U9["median_startup_s"]
def total(n, u, out_int): return n * u / 1000 + (n // out_int) * OUT_S + ST_S
def budgets():
    u = units(); tb0 = total(B0_STEPS, u["B0"]["median"], B0_OUT_INT); res = {"B0": {"unit": u["B0"], "total_s": tb0}}
    for k in ("S2", "M64", "M64S2"):
        per = u[k]["median"] / 1000 + OUT_S / ARM_OUT_INT
        n = math.ceil((tb0 - ST_S) / per / ARM_OUT_INT) * ARM_OUT_INT
        res[k] = {"unit": u[k], "budget": n}
    (D / "tt_budgets.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    for k in ("S2", "M64", "M64S2"): print(k, res[k]["budget"])
def j_dist(h5):
    with h5py.File(h5, "r") as h:
        r = np.abs(np.asarray(h["VALUE/res_ro"][:])).reshape(NI, NJ); c = np.asarray(h["MESH/COORD"][:]).reshape(-1, 3)
    xw = c[:, 0].reshape(NI, NJ)[:, -1] / S; tot = float(r.sum())
    return {"sum_2pi": tot * 2 * math.pi, "frac_j20_59": float(r[:, 20:60].sum() / tot), "frac_x_m5_0": float(r[(xw >= -5) & (xw < 0)].sum() / tot)}
def tail_rms(run, lo, hi):
    with open(HERE / run / "residual_history.csv") as f:
        rd = csv.DictReader(f); cols = [c for c in rd.fieldnames if c.startswith("rms_") and not c.startswith("rms_dq")]
        v = collections.defaultdict(list)
        for r in rd:
            if r["phase"] != "outer_end": continue
            s = int(r["step"])
            if lo <= s < hi:
                for c in cols: v[c].append(float(r[c]))
    return {c: statistics.median(x) for c, x in v.items() if x}
def final():
    B = json.loads((D / "tt_budgets.json").read_text()); u = units()
    b0s = {r["step"]: r for r in json.loads((D / "series_run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2.json").read_text())["rows"]}
    b0 = b0s[B0_STEPS - 75000]
    b0_tail = tail_rms("run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2", 30000, 50000)
    b0_dist = j_dist(HERE / "run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2" / "res_60000.h5")
    tb0 = total(B0_STEPS, u["B0"]["median"], B0_OUT_INT); rb0 = B0_OUT_INT * u["B0"]["median"] / 1000 + B0_STEPS * u["B0"]["spread"] / 1000
    out = {"B0": {"steps": B0_STEPS, "total_s": tb0, "res_s": rb0, "theta_q": {k: b0[k] for k in ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")}, "res_ro_dist_135000": b0_dist, "tail_rms": b0_tail}, "arms": {}}
    for k, run in (("S2", "run_0377_tt_S2"), ("M64", "run_0378_tt_M64"), ("M64S2", "run_0379_tt_M64S2")):
        w = json.loads((HERE / run / "m9_watch.json").read_text()); a = {"run": run, "status": w["status"], "budget": B[k]["budget"]}
        conv = (HERE / run / "CONVERGENCE_VERDICT.txt"); a["convergence"] = ("DIVERGED" if "DIVERGED" in conv.read_text() else "NOT CONVERGED" if "NOT CONVERGED" in conv.read_text() else "PASS" if "-> PASS" in conv.read_text() else "記録なし") if conv.exists() else "記録なし"
        if w["status"] == "CENSORED": a["verdict"] = "遅い (損益分岐の予算で打ち切り)"
        elif w["status"] != "REACHED" or a["convergence"] in ("DIVERGED", "記録なし"): a["verdict"] = f"比較不可 ({w['status']}・{a['convergence']})"
        else:
            n = int(w["reach_step"]); um = u[k]["median"]; t = total(n, um, ARM_OUT_INT); r = ARM_OUT_INT * um / 1000 + n * u[k]["spread"] / 1000
            a.update(steps=n, total_s=t, res_s=r, step_ratio=n / B0_STEPS)
            a["verdict"] = "速い" if t < tb0 - (r + rb0) else "遅い" if t > tb0 + (r + rb0) else "判別不能"
            row = next(x for x in w["rows"] if x["step"] == n)
            dq = {kk: 100 * (row[kk] / b0[kk] - 1) for kk in ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")}
            dist = j_dist(HERE / run / f"res_{n}.h5"); tail = tail_rms(run, n - 20000, n + 1)
            flags = [f"θ_r {kk} {v:+.3f} %" for kk, v in dq.items() if kk != "Q_w" and abs(v) > 0.05]
            if abs(dq["Q_w"]) > 0.1: flags.append(f"Q_w {dq['Q_w']:+.3f} %")
            if dist["sum_2pi"] > 2 * b0_dist["sum_2pi"]: flags.append(f"Σ|res_ro| {dist['sum_2pi']:.3g} > 2 × B0 {b0_dist['sum_2pi']:.3g}")
            for c, v in tail.items():
                if c in b0_tail and b0_tail[c] > 0 and v > 2 * b0_tail[c]: flags.append(f"{c} 中央値 {v:.3g} > 2 × B0 {b0_tail[c]:.3g}")
            a.update(diff_vs_b0_pct=dq, res_ro_dist=dist, tail_rms=tail, r3_flags=flags)
            if k in ("S2", "M64") and 0.9 <= n / B0_STEPS <= 1.1: a["r4"] = "step 数が B0 の 0.9〜1.1 倍: step 数への効果は断定しない"
        out["arms"][k] = a
    cand = [(v["total_s"], k) for k, v in out["arms"].items() if v["verdict"] == "速い" and not v.get("r3_flags")]
    out["r2"] = f"採用の候補: {min(cand)[1]}" if cand else "速くフラグの無い腕なし → B0 のまま"
    (D / "tt_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float)); print(json.dumps(out, indent=1, ensure_ascii=False, default=float))
{"budgets": budgets, "final": final}[sys.argv[1]]()
