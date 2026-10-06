"""Euler A/B の判定関数 (plan tooling-nozzle-throat-monotone-r2 §6 E3・E4、§5.1 #5 plan レビュー M2)。
評価器 (`eval_wallfit_euler.py --e3`) から呼ぶ純粋関数だけを置く (import しても何も走らない)。試験は test_throat_mono_judge.py。

量はどれも「大きいほど悪い」(|P 傾き|・オーバーシュート・波・|出口 M − 6|)。腕 A = 現行壁、腕 B = 単調壁。
  U = max(3R_A, 3R_B, 2T, 2E_fixed, E_exit,q, Δq/10)   (R: 腕内の幅、T: 末尾 5 枚の幅の最大、E: 評価刻み感度の最大)
  採用 (候補として非劣化): (q_B − q_A) + U ≤ Δq
  不採用:                   q_B − q_A ≥ Δq + U
  それ以外:                 保留 (諮問)
  η0 の量は T ≤ Δq/4 のときだけ判定に入れる (超えたら除外 = 判定に入れない。V3 と同じ)。
  末尾 5 枚が STEADY でない量は判定に入れず保留 (E4)。腕の run 数が足りない・値が非有限なら欠損 (保留)。
総合: どれか不採用 → 不採用 / 保留・未定常・欠損が 1 つでもあれば保留 / 判定に入った量がすべて採用 → 採用。
"""
import math

# E2 の許容悪化幅 Δq (量の列名 → Δq)。η は 0.0 / 0.1 の 2 本
DELTA_Q = {"P_slope_abs": 0.03, "overshoot": 0.003, "overshoot_exitnorm": 0.003, "M_wave": 0.001, "P_wave": 0.010,
           "exit_M_dev": 0.00018}
EXIT_QUANTITIES = ("exit_M_dev", "overshoot_exitnorm_eta0.0", "overshoot_exitnorm_eta0.1")   # U に E_exit,q を入れる量


def delta_key(col: str):
    """列名 → DELTA_Q のキー (判定しない列は None)。"""
    if col == "exit_M_dev":
        return "exit_M_dev"
    for k in ("P_slope_abs", "overshoot_exitnorm", "M_wave", "P_wave"):
        if col.startswith(k + "_eta"):
            return k
    if col.startswith("overshoot_eta"):
        return "overshoot"
    return None


def is_eta0(col: str) -> bool:
    return col.endswith("_eta0.0")


def _finite(v) -> bool:
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def judge_quantity(col: str, a_vals, b_vals, T: float, E: float, steady: bool, delta: float | None = None,
                   E_exit: float = 0.0, n_required: int = 3) -> dict:
    """1 量の E3 判定。a_vals / b_vals: 各腕の run ごとの代表値 (末尾 5 枚平均)。T / E: 全 run の最大。
    steady: 全 run で末尾 5 枚 STEADY。戻り値の status: adopt / reject / hold / excluded / unsteady / missing。"""
    if delta is None:
        dk = delta_key(col)
        if dk is None:
            raise ValueError(f"{col} は判定対象の量でない")
        delta = DELTA_Q[dk]
    if not (_finite(delta) and float(delta) > 0):
        raise ValueError(f"{col}: Δq = {delta!r} は正の有限値であること")
    a = [float(v) if _finite(v) else math.nan for v in (a_vals or [])]
    b = [float(v) if _finite(v) else math.nan for v in (b_vals or [])]
    row = {"col": col, "Delta": float(delta), "n_A": len(a), "n_B": len(b)}
    bad = [nm for nm, v in (("T", T), ("E", E), ("E_exit", E_exit)) if not _finite(v)]
    if len(a) < n_required or len(b) < n_required or any(math.isnan(v) for v in a + b) or bad:
        row.update(status="missing", reason=(f"run 数 A {len(a)} / B {len(b)} (< {n_required})" if len(a) < n_required or len(b) < n_required
                                             else f"非有限値 ({', '.join(bad) or '代表値'})"))
        return row
    for nm, v in (("T", T), ("E", E), ("E_exit", E_exit)):
        if float(v) < 0:
            raise ValueError(f"{col}: {nm} = {v} が負 (幅・感度は非負)")
    A_ = sum(a) / len(a)
    B_ = sum(b) / len(b)
    R_A, R_B = max(a) - min(a), max(b) - min(b)
    U = max(3 * R_A, 3 * R_B, 2 * float(T), 2 * float(E), float(E_exit), delta / 10)
    d = B_ - A_
    row.update(A_mean=A_, B_mean=B_, B_minus_A=d, R_A=R_A, R_B=R_B, T=float(T), E=float(E), E_exit=float(E_exit), U=U,
               U_terms={"3R_A": 3 * R_A, "3R_B": 3 * R_B, "2T": 2 * float(T), "2E": 2 * float(E), "E_exit": float(E_exit),
                        "Delta/10": delta / 10})
    if is_eta0(col) and float(T) > delta / 4:
        row.update(status="excluded", reason=f"η0: T {float(T):.3g} > Δq/4 {delta / 4:.3g} (判定に入れない)")
        return row
    if not steady:
        row.update(status="unsteady", reason="末尾 5 枚が STEADY でない run がある (判定に入れず保留)")
        return row
    if d + U <= delta:
        row["status"] = "adopt"
    elif d >= delta + U:
        row["status"] = "reject"
    else:
        row["status"] = "hold"
    return row


def overall(rows) -> dict:
    """量ごとの判定 (judge_quantity の戻り値の列) → 総合判定。"""
    st = {r["col"]: r["status"] for r in rows}
    if not st:
        return {"verdict": "保留 (判定する量が無い)", "counts": {}}
    counts = {}
    for s in st.values():
        counts[s] = counts.get(s, 0) + 1
    rej = [c for c, s in st.items() if s == "reject"]
    pend = [c for c, s in st.items() if s in ("hold", "unsteady", "missing")]
    adopt = [c for c, s in st.items() if s == "adopt"]
    if rej:
        v = "不採用"
    elif pend:
        v = "保留 (諮問)"
    elif adopt:
        v = "採用 (候補として非劣化)"
    else:
        v = "保留 (判定に入った量が無い)"
    return {"verdict": v, "counts": counts, "reject": rej, "pending": pend,
            "excluded": [c for c, s in st.items() if s == "excluded"]}
