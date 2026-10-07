"""throat_mono_judge (plan tooling-nozzle-throat-monotone-r2 §6 E3・E4) の判定関数の試験。
改善・悪化・境界・未定常・欠損・η0 除外・総合判定・不正入力。値は 2 進で正確に表せるものを使い、境界の等号を固定する。
usage: python3 test_throat_mono_judge.py   (FAIL 0 で exit 0)
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from throat_mono_judge import DELTA_Q, delta_key, judge_quantity, overall  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def raises(name, fn):
    try:
        fn()
    except ValueError as e:
        check(f"{name} を拒否 ({e})", True)
        return
    check(f"{name} を拒否 (通ってしまった)", False)


D = 0.5                                 # 判定の境界を正確に作るための Δq (2 進で正確)
A3 = [1.0, 1.0, 1.0]
q = lambda d, **kw: judge_quantity("P_wave_eta0.1", A3, [1.0 + d] * 3, kw.pop("T", 0.0), kw.pop("E", 0.0),  # noqa: E731
                                   kw.pop("steady", True), delta=D, **kw)
# U = Δq/10 = 0.05 (幅 0 のとき)
check("改善 (B−A = −0.25) → 採用", q(-0.25)["status"] == "adopt")
check("同等 (B−A = 0) → 採用", q(0.0)["status"] == "adopt")
r = judge_quantity("P_wave_eta0.1", [1.0] * 3, [1.375] * 3, 0.0625, 0.0, True, delta=D)   # U = 2T = 0.125, d = 0.375 → d+U = 0.5 = Δq
check(f"境界: d + U = Δq (0.375 + 0.125) → 採用 (U {r['U']}, d {r['B_minus_A']})", r["status"] == "adopt" and r["U"] == 0.125)
r = judge_quantity("P_wave_eta0.1", [1.0] * 3, [1.625] * 3, 0.0625, 0.0, True, delta=D)   # d = 0.625 = Δq + U → 不採用 (≥)
check(f"境界: d = Δq + U (0.625) → 不採用", r["status"] == "reject")
r = judge_quantity("P_wave_eta0.1", [1.0] * 3, [1.5] * 3, 0.0625, 0.0, True, delta=D)     # 0.375 < d = 0.5 < 0.625 → 保留
check("中間 (Δq − U < d < Δq + U) → 保留", r["status"] == "hold")
check("悪化 (B−A = 1.0) → 不採用", q(1.0)["status"] == "reject")
# U の各項
r = judge_quantity("P_wave_eta0.1", [1.0, 1.0625, 1.0], [1.0] * 3, 0.0, 0.0, True, delta=D)
check(f"U に 3R_A が入る (R_A 0.0625 → U {r['U']})", r["U"] == 0.1875 and r["R_A"] == 0.0625)
r = judge_quantity("P_wave_eta0.1", [1.0] * 3, [1.0, 1.125, 1.0], 0.0, 0.0, True, delta=D)
check(f"U に 3R_B が入る (U {r['U']})", r["U"] == 0.375)
r = judge_quantity("P_wave_eta0.1", [1.0] * 3, [1.0] * 3, 0.0, 0.125, True, delta=D)
check(f"U に 2E が入る (U {r['U']})", r["U"] == 0.25)
r = judge_quantity("exit_M_dev", [1.0] * 3, [1.0] * 3, 0.0, 0.0, True, delta=D, E_exit=0.3125)
check(f"U に E_exit が入る (U {r['U']})", r["U"] == 0.3125)
r = judge_quantity("exit_M_dev", [1.0] * 3, [1.25] * 3, 0.0, 0.0, True, delta=D, E_exit=0.3125)
check("E_exit で U が広がり、採用できない差が保留になる (d 0.25 + U 0.3125 > Δq)", r["status"] == "hold")
# 未定常・欠損
check("末尾 5 枚が STEADY でない → unsteady (判定に入れず保留)", q(-0.25, steady=False)["status"] == "unsteady")
check("未定常は悪化でも reject にしない", q(1.0, steady=False)["status"] == "unsteady")
check("腕 B が 2 run → missing", judge_quantity("P_wave_eta0.1", A3, [1.0, 1.0], 0.0, 0.0, True, delta=D)["status"] == "missing")
check("腕 A が空 → missing", judge_quantity("P_wave_eta0.1", [], [1.0] * 3, 0.0, 0.0, True, delta=D)["status"] == "missing")
check("代表値に NaN → missing", judge_quantity("P_wave_eta0.1", A3, [1.0, math.nan, 1.0], 0.0, 0.0, True, delta=D)["status"] == "missing")
check("代表値に None → missing", judge_quantity("P_wave_eta0.1", A3, [1.0, None, 1.0], 0.0, 0.0, True, delta=D)["status"] == "missing")
check("T が inf → missing", judge_quantity("P_wave_eta0.1", A3, A3, math.inf, 0.0, True, delta=D)["status"] == "missing")
check("E が NaN → missing", judge_quantity("P_wave_eta0.1", A3, A3, 0.0, math.nan, True, delta=D)["status"] == "missing")
# η0 の除外 (T ≤ Δq/4 のときだけ判定)
r = judge_quantity("P_wave_eta0.0", A3, [3.0] * 3, 0.125, 0.0, True, delta=D)     # T = Δq/4 ちょうど → 判定に入る
check(f"η0: T = Δq/4 ちょうど → 判定に入る ({r['status']})", r["status"] == "reject")
r = judge_quantity("P_wave_eta0.0", A3, [3.0] * 3, 0.125 + 2 ** -10, 0.0, True, delta=D)
check(f"η0: T > Δq/4 → 除外 ({r['status']})", r["status"] == "excluded")
r = judge_quantity("P_wave_eta0.1", A3, [3.0] * 3, 0.125 + 2 ** -10, 0.0, True, delta=D)
check("η0.1 は T > Δq/4 でも除外しない", r["status"] == "reject")
r = judge_quantity("P_wave_eta0.0", A3, A3, 0.375, 0.0, False, delta=D)
check("η0 で除外に当たるなら未定常でも除外 (判定に入れない)", r["status"] == "excluded")
# 不正入力
raises("Δq = 0", lambda: judge_quantity("P_wave_eta0.1", A3, A3, 0.0, 0.0, True, delta=0.0))
raises("Δq < 0", lambda: judge_quantity("P_wave_eta0.1", A3, A3, 0.0, 0.0, True, delta=-1.0))
raises("Δq = NaN", lambda: judge_quantity("P_wave_eta0.1", A3, A3, 0.0, 0.0, True, delta=math.nan))
raises("T < 0", lambda: judge_quantity("P_wave_eta0.1", A3, A3, -1.0, 0.0, True, delta=D))
raises("判定対象でない列 (exit_core_M)", lambda: judge_quantity("exit_core_M", A3, A3, 0.0, 0.0, True))
# Δq の割り当て (E2)
want = {"P_slope_abs_eta0.0": 0.03, "P_slope_abs_eta0.1": 0.03, "overshoot_eta0.0": 0.003, "overshoot_eta0.1": 0.003,
        "overshoot_exitnorm_eta0.0": 0.003, "overshoot_exitnorm_eta0.1": 0.003, "M_wave_eta0.0": 0.001, "M_wave_eta0.1": 0.001,
        "P_wave_eta0.0": 0.010, "P_wave_eta0.1": 0.010, "exit_M_dev": 0.00018}
check("E2 の Δq の割り当て", all(DELTA_Q[delta_key(c)] == v for c, v in want.items()))
check("判定しない列は None (exit_core_M, P_slope_eta0.0 = 符号付き)", delta_key("exit_core_M") is None and delta_key("P_slope_eta0.0") is None)
# 総合判定
mk = lambda c, s: {"col": c, "status": s}  # noqa: E731
check("総合: 全量採用 → 採用", overall([mk("a", "adopt"), mk("b", "adopt")])["verdict"].startswith("採用"))
check("総合: 採用 + η0 除外 → 採用", overall([mk("a", "adopt"), mk("b", "excluded")])["verdict"].startswith("採用"))
check("総合: 1 つ不採用 → 不採用 (保留・欠損があっても)", overall([mk("a", "adopt"), mk("b", "reject"), mk("c", "missing")])["verdict"] == "不採用")
check("総合: 保留がある → 保留", overall([mk("a", "adopt"), mk("b", "hold")])["verdict"].startswith("保留"))
check("総合: 未定常がある → 保留", overall([mk("a", "adopt"), mk("b", "unsteady")])["verdict"].startswith("保留"))
check("総合: 欠損がある → 保留", overall([mk("a", "adopt"), mk("b", "missing")])["verdict"].startswith("保留"))
check("総合: 全量除外 → 保留 (採用にしない)", overall([mk("a", "excluded")])["verdict"].startswith("保留"))
check("総合: 空 → 保留", overall([])["verdict"].startswith("保留"))

# --- §5.1 #5b: 判定区間の収束判定 (E4、諮問 ③) ------------------------------------------------------------------------
from throat_mono_judge import (ABS_TOL_WITHIN, ICAB_U_FLOOR_FRAC, abs_tolerance_verdict, check_run_preconditions,  # noqa: E402
                               icab_overall, judge_icab_quantity, mono_r2_matches, parse_segment_verdict)
SEGHEAD = "  [segment] 判定区間 = main  (18000 行) -> run_x/residual_history_segment.csv\n\n"
PLATEAU = SEGHEAD + ("=== run_x  [last step 17999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===\n"
                     "  rms_ro      : init=1e-05 fin=4e-07 drop= 1.4dec flat     <-- STALLED (plateau)\n"
                     "  rms_roUz    : all-zero (inactive, skip)\n\nOVERALL: CHECK FAILURES ABOVE\n")
seg = lambda t: parse_segment_verdict(t)["status"]  # noqa: E731
check("区間判定: 既知の plateau → plateau", seg(PLATEAU) == "plateau")
check("区間判定: PASS → pass", seg(SEGHEAD + "=== run_x  [last step 17999]  -> PASS (converged) ===\n  rms_ro : ok\n") == "pass")
check("区間判定: DIVERGED → diverged (plateau と同列にしない)",
      seg(SEGHEAD + "=== run_x  [last step 99]  -> DIVERGED (NaN/Inf) ===\n  rms_ro      : NaN/Inf present  <-- DIVERGED\n") == "diverged")
check("区間判定: 全体行が plateau でも RISING の列がある → rising",
      seg(PLATEAU.replace("rms_roUz    : all-zero (inactive, skip)", "rms_roe : init=1 fin=9 drop=0.0dec rising  <-- RISING (divergent)")) == "rising")
check("区間判定: 判定不能の行 → undeterminable",
      seg(SEGHEAD + "=== run_x  [last step 1]  -> NOT CONVERGED ===\n  (入力)      : 必須の保存量残差列が無い: rms_roe  <-- 判定不能\n") == "undeterminable")
check("区間判定: --segment の区間の行が無い (通常判定の出力) → undeterminable", seg(PLATEAU.replace(SEGHEAD, "")) == "undeterminable")
check("区間判定: stage_manifest が無い (判定不能の出力のみ) → undeterminable",
      seg("[run_x] stage_manifest.json が無い -> --segment は使えない (判定区間を人が明示すること)  <-- 判定不能\n\nOVERALL: CHECK FAILURES ABOVE\n") == "undeterminable")
check("区間判定: ファイル無し (None) → missing", seg(None) == "missing")
check("区間判定: 空文字 → undeterminable", seg("") == "undeterminable")
check("区間判定: 全体行が 2 行 (別 run の出力が混在) → undeterminable", seg(PLATEAU + PLATEAU.replace(SEGHEAD, "")) == "undeterminable")
check("区間判定: still converging → converging (保留側)",
      seg(SEGHEAD + "=== run_x  [last step 17999]  -> NOT CONVERGED (still converging — run more steps) ===\n") == "converging")

# --- 近零量の「絶対許容内」(E4、諮問 ②) -----------------------------------------------------------------------------
TOL = 2.0 ** -16                                          # 2 進で正確な許容 (1.5e-5 相当)
st10 = [1000.0 * k for k in range(9, 19)]
av = lambda v, **kw: abs_tolerance_verdict(st10[-len(v):] if len(v) <= 10 else [1000.0 * k for k in range(len(v))], v, TOL, **kw)  # noqa: E731
osc = [TOL / 4 * (1 if k % 2 else -1) for k in range(10)]  # 振動 (単調でない)、幅 TOL/2
check(f"絶対許容内: 小さい振動 → {ABS_TOL_WITHIN}", av(osc)["status"] == ABS_TOL_WITHIN)
check("絶対許容内: 一定値 (増分 0) → 絶対許容内 (単調扱いしない)", av([TOL / 2] * 10)["status"] == ABS_TOL_WITHIN)
check("絶対許容内: 9 枚 → 適用外 (枚数不足)", av(osc[:9])["status"] == "適用外")
check("絶対許容内: 最大絶対値 > tol → 適用外 (元の VERDICT で扱う)", av([TOL * 1.5] + osc[1:])["status"] == "適用外")
check("絶対許容内: 最大絶対値 = tol ちょうど → 適用範囲 (≤)",
      av([TOL] + [TOL / 2] * 8 + [TOL])["status"] == ABS_TOL_WITHIN)
wide = [-TOL * 0.75, TOL * 0.75] * 5                      # |値| ≤ tol だが幅 1.5 tol
check("絶対許容内: 幅 > tol → 許容外", av(wide)["status"] == "許容外")
step = [-TOL * 0.6] * 5 + [TOL * 0.45] * 5                # 幅 1.05 tol > tol
check("絶対許容内: 隣接 5 枚の平均差 > tol → 許容外", av(step)["status"] == "許容外")
lin = [TOL * (-0.4 + 0.08 * k) for k in range(10)]        # 単調・増分一定 (減衰しない)、幅 0.72 tol
check(f"絶対許容内: 単調で減衰しない → 保留 ({av(lin)['status']})", av(lin)["status"] == "保留 (単調・減衰なし)")
dec = [TOL * (0.5 - 0.5 * 0.5 ** k) for k in range(10)]   # 単調・増分が半減 (減衰)
check(f"絶対許容内: 単調でも増分が減衰 → 絶対許容内 ({av(dec)['status']})", av(dec)["status"] == ABS_TOL_WITHIN)
slow = [TOL * (0.5 - 0.5 * 0.9 ** k) for k in range(10)]   # 単調・増分の比 0.9 (後半/前半 = 0.59 > 0.5): ゆっくりしか減衰しない
check(f"絶対許容内: 単調で減衰が遅い (比 0.9) → 保留 ({av(slow)['status']})", av(slow)["status"] == "保留 (単調・減衰なし)")
q8 = [TOL * (0.5 - 0.5 * 0.8 ** k) for k in range(10)]     # 比 0.8 (後半/前半 = 0.33 ≤ 0.5)
check(f"絶対許容内: 単調で比 0.8 の減衰 → 絶対許容内 ({av(q8)['status']})", av(q8)["status"] == ABS_TOL_WITHIN)
jit = [TOL * (-0.4 + 0.08 * k + (2 ** -30 if k in (6, 8) else 0.0)) for k in range(10)]   # 一定速度 + 丸め程度のジッタ
check(f"絶対許容内: 一定速度のドリフトは丸めのジッタがあっても保留 ({av(jit)['status']})", av(jit)["status"] == "保留 (単調・減衰なし)")
grow = [TOL * 0.01 * 1.5 ** k for k in range(10)]         # 単調・増分が増える
check("絶対許容内: 単調で増分が増える → 保留", av(grow)["status"] == "保留 (単調・減衰なし)")
check("絶対許容内: 非有限 → 非有限 (保留)", av(osc[:5] + [math.nan] + osc[6:])["status"] == "非有限")
check("絶対許容内: None → 非有限 (保留)", av(osc[:5] + [None] + osc[6:])["status"] == "非有限")
check("絶対許容内: step が増えない → 非有限 (保留)",
      abs_tolerance_verdict([0.0] * 10, osc, TOL)["status"] == "非有限")
check("絶対許容内: 末尾 10 枚だけを見る (前の大きな過渡は無視)", av([1.0] * 5 + osc)["status"] == ABS_TOL_WITHIN)
r = av(osc)
check("絶対許容内: 幅・平均差・最大絶対値を記録", all(k in r for k in ("width", "adjacent_mean_diff", "max_abs")))
raises("tol = 0", lambda: abs_tolerance_verdict(st10, osc, 0.0))
raises("tol = NaN", lambda: abs_tolerance_verdict(st10, osc, math.nan))
raises("tol < 0", lambda: abs_tolerance_verdict(st10, osc, -1.0))
raises("steps と値の数が違う", lambda: abs_tolerance_verdict(st10[:9], osc, TOL))
raises("n < 2·half", lambda: abs_tolerance_verdict(st10, osc, TOL, n=8, half=5))

# --- 前提 (諮問 ③) -----------------------------------------------------------------------------------------------------
check("mono_r2: None は None だけ", mono_r2_matches(None, None) and not mono_r2_matches([0.0, 1.5], None) and not mono_r2_matches(0, None))
check("mono_r2: [0.0, 1.5] と完全一致 (int 0 も数値として一致)", mono_r2_matches([0.0, 1.5], [0.0, 1.5]) and mono_r2_matches((0, 1.5), [0.0, 1.5]))
check("mono_r2: 値違い・長さ違い・文字列・bool・None は不一致",
      not any(mono_r2_matches(v, [0.0, 1.5]) for v in ([0.0, 1.0], [0.0, 1.5, 2.0], ["0", "1.5"], [False, 1.5], None, [0.0, math.nan])))
EV = {"status": "consistent", "resid_over_f32tol_max": 0.5, "vs_other": {"status": "ok", "actual_dr_max_m": 4.99e-7}}
GOOD = {"wall_fit_has_mono_r2": True, "mono_r2": [0.0, 1.5], "wall_evidence": EV, "exit_base_check_maxabs": 0.0,
        "segment": {"status": "plateau"}, "ic": {"tool": "ic_index_map", "mode": "index", "VERDICT": "OK"}, "dry": False}
ICB = {"tool": "ic_index_map", "mode": "index"}
pre = lambda **kw: check_run_preconditions({**GOOD, **kw}, [0.0, 1.5], ICB)  # noqa: E731
check("前提: そろっていれば成立", pre() == [])
check("前提: segment pass も成立", pre(segment={"status": "pass"}) == [])
check("前提: A (mono_r2 None・restart_field) も成立",
      check_run_preconditions({**GOOD, "mono_r2": None, "ic": {"tool": "restart_field", "VERDICT": "OK"}}, None, {"tool": "restart_field"}) == [])
check("前提: A の mono_r2 が [0, 1.5] → 不成立",
      check_run_preconditions({**GOOD, "ic": {"tool": "restart_field", "VERDICT": "OK"}}, None, {"tool": "restart_field"}) != [])
check("前提: wall_fit に mono_r2 の記録が無い → 不成立 (None と区別)", pre(wall_fit_has_mono_r2=False, mono_r2=None) != []
      and check_run_preconditions({**GOOD, "wall_fit_has_mono_r2": False, "mono_r2": None, "ic": {"tool": "restart_field", "VERDICT": "OK"}},
                                  None, {"tool": "restart_field"}) != [])
check("前提: B の mono_r2 が [0, 1.0] → 不成立", pre(mono_r2=[0.0, 1.0]) != [])
check("前提: 出口照合が None → 不成立 (or 0.0 で通さない)", pre(exit_base_check_maxabs=None) != [])
check("前提: 出口照合が NaN → 不成立", pre(exit_base_check_maxabs=math.nan) != [])
check("前提: 出口照合が 1e-9 → 不成立", pre(exit_base_check_maxabs=1e-9) != [])
check("前提: 壁の証拠が mismatch → 不成立", pre(wall_evidence={**EV, "status": "mismatch"}) != [])
check("前提: 壁の証拠が無い → 不成立", pre(wall_evidence=None) != [])
check("前提: 壁の証拠の残差比が NaN → 不成立", pre(wall_evidence={**EV, "resid_over_f32tol_max": math.nan}) != [])
check("前提: 他腕との照合が無い → 不成立 (require_vs_other)", pre(wall_evidence={k: v for k, v in EV.items() if k != "vs_other"}) != [])
check("前提: 他腕との照合が無くても require_vs_other=False なら成立 (予備 A/B)",
      check_run_preconditions({**GOOD, "wall_evidence": {k: v for k, v in EV.items() if k != "vs_other"}}, [0.0, 1.5], ICB,
                              require_vs_other=False) == [])
for st in ("diverged", "missing", "undeterminable", "rising", "converging", "other"):
    check(f"前提: 区間判定 {st} → 不成立 (保留)", pre(segment={"status": st, "reason": st}) != [])
check("前提: 区間判定が None → 不成立", pre(segment=None) != [])
check("前提: IC の道具違い (interp_field) → 不成立", pre(ic={"tool": "interp_field", "VERDICT": "OK"}) != [])
check("前提: IC のモード違い (nearest) → 不成立", pre(ic={"tool": "ic_index_map", "mode": "nearest", "VERDICT": "OK"}) != [])
check("前提: IC の VERDICT が無い → 不成立", pre(ic={"tool": "ic_index_map", "mode": "index"}) != [])
check("前提: IC の記録が無い → 不成立", pre(ic=None) != [])
check("前提: 乾式確認の prep → 不成立", pre(dry=True) != [])

# --- 予備 A/B (§6 E1) --------------------------------------------------------------------------------------------------
DQ = 0.5
ic_ = lambda d, fl=0.0, **kw: judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0 + d, kw.pop("T", 0.0), kw.pop("E", 0.0),  # noqa: E731
                                                   kw.pop("steady", True), delta=DQ, u_floor_frac=fl, **kw)
# Δq/10 = 0.05。下限項を外した (fl = 0) ときの 3 区分
check("予備 A/B: 差 0・U 0 → negligible", ic_(0.0)["status"] == "negligible")
r = ic_(0.03125, T=0.009765625)                            # |d| + U = 0.03125 + 0.01953125 = 0.05078 > 0.05
check(f"予備 A/B: |d| + U がわずかに Δq/10 を超える → indeterminate ({r['status']}, U {r['U']})", r["status"] == "indeterminate")
r = ic_(0.03125, T=0.0078125)                              # |d| + U = 0.03125 + 0.015625 = 0.046875 ≤ 0.05
check(f"予備 A/B: |d| + U ≤ Δq/10 → negligible ({r['status']})", r["status"] == "negligible")
r = ic_(-0.25, T=0.0625)                                   # |d| − U = 0.25 − 0.125 = 0.125 > 0.05 (符号によらない)
check(f"予備 A/B: |d| − U > Δq/10 (負の差) → dependent ({r['status']})", r["status"] == "dependent")
r = ic_(0.25, T=0.0625)
check("予備 A/B: |d| − U > Δq/10 (正の差) → dependent", r["status"] == "dependent")
r = ic_(0.125, T=0.0390625)                                # |d| − U = 0.125 − 0.078125 = 0.046875 ≤ 0.05 → 境界域
check(f"予備 A/B: |d| − U ≤ Δq/10 かつ |d| + U > Δq/10 → indeterminate ({r['status']})", r["status"] == "indeterminate")
r = ic_(0.0, E=0.03125)
check(f"予備 A/B: U に 2E (U {r['U']}) → |d| + U = 0.0625 > 0.05 → indeterminate", r["U"] == 0.0625 and r["status"] == "indeterminate")
r = ic_(0.0, E_exit=0.04)
check(f"予備 A/B: U に E_exit (U {r['U']})", r["U"] == 0.04 and r["status"] == "negligible")
r = ic_(0.0, T=0.0078125)
check("予備 A/B: U の項は 2T・2E・E_exit・下限で、3R_A・3R_B は無い (各 1 本)", set(r["U_terms"]) == {"2T", "2E", "E_exit", "floor"} and "3R" not in "".join(r["U_terms"]))
check("予備 A/B: 未定常 → unsteady (差が大きくても dependent にしない)", ic_(1.0, steady=False)["status"] == "unsteady")
check("予備 A/B: 代表値 NaN → missing", judge_icab_quantity("P_wave_eta0.1", math.nan, 1.0, 0.0, 0.0, True, delta=DQ)["status"] == "missing")
check("予備 A/B: T None → missing", judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0, None, 0.0, True, delta=DQ)["status"] == "missing")
check("予備 A/B: E_exit inf → missing", judge_icab_quantity("exit_M_dev", 1.0, 1.0, 0.0, 0.0, True, E_exit=math.inf, delta=DQ)["status"] == "missing")
raises("予備 A/B: T < 0", lambda: judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0, -1.0, 0.0, True, delta=DQ))
raises("予備 A/B: Δq = 0", lambda: judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0, 0.0, 0.0, True, delta=0.0))
raises("予備 A/B: 判定対象でない列", lambda: judge_icab_quantity("exit_core_M", 1.0, 1.0, 0.0, 0.0, True))
raises("予備 A/B: u_floor_frac < 0", lambda: judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0, 0.0, 0.0, True, delta=DQ, u_floor_frac=-0.1))
check("予備 A/B: Δq は E2 の割り当て (exit_M_dev 1.8e-4 → 閾値 1.8e-5)",
      abs(judge_icab_quantity("exit_M_dev", 1.0, 1.0, 0.0, 0.0, True)["threshold"] - 1.8e-5) < 1e-18)
# 既定の下限項は 0 (plan §6 E1、2026-10-06 決定: 予備 A/B の U は評価・時間変動の不確かさだけ)
check(f"予備 A/B 既定の下限項は 0 (実際 {ICAB_U_FLOOR_FRAC})", ICAB_U_FLOOR_FRAC == 0.0)
_thr = DQ / 10.0
check("予備 A/B 既定: 差 0.5·閾値・U 0 → negligible",
      judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0 + 0.5 * _thr, 0.0, 0.0, True, delta=DQ)["status"] == "negligible")
check("予備 A/B 既定: 差 0.5·閾値・2T 0.6·閾値 → indeterminate (|差| + U > 閾値、|差| − U ≤ 閾値)",
      judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0 + 0.5 * _thr, 0.3 * _thr, 0.0, True, delta=DQ)["status"] == "indeterminate")
check("予備 A/B 既定: 差 2·閾値・U 0.5·閾値 → dependent",
      judge_icab_quantity("P_wave_eta0.1", 1.0, 1.0 + 2.0 * _thr, 0.25 * _thr, 0.0, True, delta=DQ)["status"] == "dependent")
mk2 = lambda c, s: {"col": c, "status": s}  # noqa: E731
check("予備 A/B 総合: 全量 negligible → 残留 IC 依存説を棄却",
      icab_overall([mk2("a", "negligible"), mk2("b", "negligible")])["verdict"].startswith("残留 IC 依存説を棄却"))
check("予備 A/B 総合: 1 つ dependent → 前提を棄却 (未定常・欠損があっても)",
      icab_overall([mk2("a", "negligible"), mk2("b", "dependent"), mk2("c", "unsteady")])["verdict"].startswith("前提を棄却"))
for s_ in ("indeterminate", "unsteady", "missing"):
    check(f"予備 A/B 総合: {s_} が 1 つ → 判別不能", icab_overall([mk2("a", "negligible"), mk2("b", s_)])["verdict"].startswith("判別不能"))
check("予備 A/B 総合: 前提未達 → 判別不能 (dependent があっても)",
      icab_overall([mk2("a", "dependent")], preconditions_ok=False)["verdict"].startswith("判別不能 (前提未達"))
check("予備 A/B 総合: 空 → 判別不能", icab_overall([])["verdict"].startswith("判別不能"))
print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
