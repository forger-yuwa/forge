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
print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
