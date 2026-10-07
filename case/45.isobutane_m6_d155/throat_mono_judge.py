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
§5.1 #5b (諮問 2026-10-06 ②③、§6 E1・E4): 判定区間の収束判定の解釈 (parse_segment_verdict)、近零量の「絶対許容内」
(abs_tolerance_verdict; 元の VERDICT は残す)、run ごとの前提 (check_run_preconditions; 証拠の欠損・非有限値は既定値で通さない)、
IC 写像の予備 A/B の判定 (judge_icab_quantity・icab_overall)。
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


# --- §5.1 #5b (諮問 2026-10-06 ②③、§6 E1・E4) ------------------------------------------------------------------------
SEGMENT_VERDICT_FILE = "CONVERGENCE_VERDICT_segment.txt"   # 投入スクリプトが check_convergence --segment の出力を書き、評価器が読む (同じ区間)
ACCEPTED_SEGMENT = ("pass", "plateau")                      # plateau の NOT CONVERGED は既知として記録のみ (E4)。ほかはすべて保留
ABS_TOL_N = 10                                              # 近零量の「絶対許容内」: 同一設定区間の連続 10 枚以上
ABS_TOL_HALF = 5                                            # 隣接 5 枚の平均差
ABS_TOL_WITHIN = "絶対許容内"                                # STEADY とは表示しない (E4)
# 「減衰」の定義: 単調なとき、後半の増分 (step あたり) の絶対値の和が前半の ABS_TOL_DECAY_RATIO 倍以下。
# 厳密な「後半 < 前半」では、一定速度のドリフトが丸め・ジッタで半々に「減衰」と判定される (再現しない判定) ので余裕を持たせる。
# 0.5 は 1 枚あたりの比 q ≈ 0.87 (q^5 = 0.5) に当たる。plan に数値の定めが無いので実装上の定義 (§5.1 #5b 報告)。
ABS_TOL_DECAY_RATIO = 0.5
EXIT_BASE_TOL = 1e-12                                       # 出口の既定帯と quantities の照合 (評価定義の一致)


def parse_segment_verdict(text) -> dict:
    """check_convergence.py <run> --segment の出力 → {"status", "line", "reason"}。
    status: pass / plateau (既知、記録のみ) / diverged / rising / converging / undeterminable / missing / other。
    plateau は「全体行が stalled/plateau」かつ列に RISING・DIVERGED・判定不能が無いときだけ (check_convergence は stalled の列が
    1 つでもあれば RISING の列があっても全体行を plateau にするので、列の行も見る)。"""
    if text is None:
        return {"status": "missing", "line": None, "reason": f"{SEGMENT_VERDICT_FILE} が無い"}
    if "[segment] 判定区間" not in text:
        return {"status": "undeterminable", "line": None, "reason": "--segment の判定区間の行が無い (区間が不明)"}
    heads = [l.strip() for l in text.splitlines() if l.startswith("=== ") and "-> " in l]
    if len(heads) != 1:
        return {"status": "undeterminable", "line": None, "reason": f"判定の全体行が {len(heads)} 行 (1 run 1 行であること)"}
    line = heads[0]
    v = line.split("-> ", 1)[1]
    if "判定不能" in text:
        return {"status": "undeterminable", "line": line, "reason": "判定不能の行がある"}
    if "DIVERGED" in v or "NaN/Inf present" in text:
        return {"status": "diverged", "line": line, "reason": "DIVERGED (NaN/Inf)"}
    if "RISING" in text:
        return {"status": "rising", "line": line, "reason": "RISING の列がある"}
    if v.startswith("PASS"):
        return {"status": "pass", "line": line, "reason": None}
    if v.startswith("NOT CONVERGED (stalled/plateau"):
        return {"status": "plateau", "line": line, "reason": "既知の plateau (記録のみ)"}
    if "still converging" in v:
        return {"status": "converging", "line": line, "reason": "まだ下がっている (区間が足りない)"}
    return {"status": "other", "line": line, "reason": "解釈できない全体行"}


def abs_tolerance_verdict(steps, signed_vals, tol, n=ABS_TOL_N, half=ABS_TOL_HALF) -> dict:
    """近零量の「絶対許容内」判定 (§6 E4、諮問 ②)。元の check_quasisteady の VERDICT は別に残す (置き換えない)。
    signed_vals: 符号付きの量 (|出口 M − 6| なら M_exit − 6) の時系列。末尾 n 枚 (同一設定区間の連続 n 枚であることは呼び出し側が保証) で:
      適用範囲: n 枚そろい、最大絶対値 ≤ tol。外れたら「適用外」(元の VERDICT だけで扱う)。
      条件: n 枚の幅 (max − min) と隣接 half 枚の平均差 |mean(末尾 half) − mean(その前の half)| がともに ≤ tol。外れたら「許容外」。
      単調 (増分がすべて ≥ 0 か ≤ 0、0 でない増分あり) で、増分 (step あたり) の絶対値の和が後半で前半の ABS_TOL_DECAY_RATIO 倍を
      超える (= 減衰しない) なら「保留 (単調・減衰なし)」。
    すべて満たせば「絶対許容内」。非有限値・step の不正は「非有限」(保留)。"""
    if not (isinstance(tol, (int, float)) and _finite(tol) and tol > 0):
        raise ValueError(f"tol = {tol!r} は正の有限値であること")
    if n < 2 * half or half < 1:
        raise ValueError(f"n {n} は 2·half {2 * half} 以上")
    s = [float(x) if _finite(x) else math.nan for x in (steps or [])]
    v = [float(x) if _finite(x) else math.nan for x in (signed_vals or [])]
    out = {"tol": float(tol), "n": n, "n_available": len(v)}
    if len(v) != len(s):
        raise ValueError(f"steps {len(s)} と値 {len(v)} の数が違う")
    if len(v) < n:
        out.update(status="適用外", reason=f"枚数 {len(v)} < {n}")
        return out
    s, v = s[-n:], v[-n:]
    if any(math.isnan(x) for x in v + s):
        out.update(status="非有限", reason="非有限の値・step がある")
        return out
    ds = [b - a for a, b in zip(s, s[1:])]
    if any(d <= 0 for d in ds):
        out.update(status="非有限", reason="step が狭義単調増加でない")
        return out
    maxabs = max(abs(x) for x in v)
    width = max(v) - min(v)
    mdiff = abs(sum(v[-half:]) / half - sum(v[-2 * half:-half]) / half)
    out.update(max_abs=maxabs, width=width, adjacent_mean_diff=mdiff, steps=[s[0], s[-1]])
    if maxabs > tol:
        out.update(status="適用外", reason=f"最大絶対値 {maxabs:.3g} > {tol:.3g}")
        return out
    if width > tol or mdiff > tol:
        out.update(status="許容外", reason=f"幅 {width:.3g} / 隣接 {half} 枚の平均差 {mdiff:.3g} のどちらかが > {tol:.3g}")
        return out
    rate = [(b - a) / d for a, b, d in zip(v, v[1:], ds)]
    mono = any(r != 0 for r in rate) and (all(r >= 0 for r in rate) or all(r <= 0 for r in rate))
    m = len(rate) // 2
    early, late = sum(abs(r) for r in rate[:m]), sum(abs(r) for r in rate[-m:])
    out.update(monotone=mono, rate_abs_sum_first_half=early, rate_abs_sum_second_half=late, decay_ratio_limit=ABS_TOL_DECAY_RATIO)
    if mono and not late <= ABS_TOL_DECAY_RATIO * early:
        out.update(status="保留 (単調・減衰なし)", reason=f"単調で、後半の増分が前半の {ABS_TOL_DECAY_RATIO} 倍を超える")
        return out
    out.update(status=ABS_TOL_WITHIN, reason=None)
    return out


def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def mono_r2_matches(value, expected) -> bool:
    """mono_r2 の完全一致: expected None なら None だけ、[a, b] なら数 2 つ (bool・文字列は不可) で値が等しいこと。"""
    if expected is None:
        return value is None
    return (isinstance(value, (list, tuple)) and len(value) == len(expected) and all(_is_num(x) for x in value)
            and all(float(x) == float(y) for x, y in zip(value, expected)))


def check_run_preconditions(meta: dict, mono_r2_expected, ic_expected: dict, require_vs_other: bool = True) -> list:
    """E3・予備 A/B の 1 run の前提 (諮問 ③)。証拠の欠損・非有限値は既定値で通さない。戻り値 = 不成立の理由 (空なら成立)。
    meta: {"wall_fit_has_mono_r2": bool, "mono_r2", "wall_evidence": dict, "exit_base_check_maxabs", "segment": parse_segment_verdict の戻り値,
           "ic": prepare_info の ic, "dry": bool}
    ic_expected: {"tool": "restart_field" | "ic_index_map", "mode": "index" | "nearest" (ic_index_map のとき)}"""
    bad = []
    if not meta.get("wall_fit_has_mono_r2"):
        bad.append("prepare_info の wall_fit に mono_r2 の記録が無い (壁の当てはめ設定の証拠が欠損)")
    elif not mono_r2_matches(meta.get("mono_r2"), mono_r2_expected):
        bad.append(f"mono_r2 {meta.get('mono_r2')!r} が期待値 {mono_r2_expected!r} と完全一致しない")
    ev = meta.get("wall_evidence")
    if not isinstance(ev, dict):
        bad.append("壁の証拠 (M4) が無い")
    else:
        if ev.get("status") != "consistent":
            bad.append(f"壁の証拠 (M4) が consistent でない ({ev.get('status')}: {ev.get('reason') or ev.get('error')})")
        r = ev.get("resid_over_f32tol_max")
        if not (_finite(r) and float(r) <= 1.0):
            bad.append(f"壁の証拠の残差比 {r!r} が有限で ≤ 1 でない")
        if require_vs_other:
            vo = ev.get("vs_other")
            if not (isinstance(vo, dict) and vo.get("status") == "ok" and _finite(vo.get("actual_dr_max_m"))):
                bad.append(f"他腕との壁の照合 (vs_other) が無い・不成立 ({(vo or {}).get('status') if isinstance(vo, dict) else vo!r})")
    x = meta.get("exit_base_check_maxabs")
    if not (_finite(x) and float(x) <= EXIT_BASE_TOL):
        bad.append(f"出口の既定帯の照合 {x!r} が有限で ≤ {EXIT_BASE_TOL:g} でない")
    seg = meta.get("segment") or {"status": "missing"}
    if seg.get("status") not in ACCEPTED_SEGMENT:
        bad.append(f"判定区間の収束判定が {seg.get('status')} ({seg.get('reason')}) — DIVERGED・欠損・判定不能は plateau と同列に扱わない")
    ic = meta.get("ic")
    if not isinstance(ic, dict):
        bad.append("IC の記録 (prepare_info の ic) が無い")
    else:
        if ic.get("tool") != ic_expected.get("tool"):
            bad.append(f"IC の道具 {ic.get('tool')!r} が期待 {ic_expected.get('tool')!r} と違う")
        if ic_expected.get("mode") is not None and ic.get("mode") != ic_expected["mode"]:
            bad.append(f"IC の写像モード {ic.get('mode')!r} が期待 {ic_expected['mode']!r} と違う")
        if ic.get("VERDICT") != "OK":
            bad.append(f"IC 写像の VERDICT が OK でない ({ic.get('VERDICT')!r})")
        if ic.get("skipped"):
            bad.append("IC を入れていない (乾式確認の prep)")
    if meta.get("dry"):
        bad.append("乾式確認 (DRY) の prep から作られた run")
    return bad


# --- 予備 A/B (§6 E1): 同じ B 格子で IC の写像だけを変えた 2 本 (α 最近傍 run_0146 / β 番号写像 run_0143) ---------------
ICAB_WIN = 10                    # 本段の連続 10 枚 (代表値 = 窓平均、T も同じ窓)
# U の下限項。plan §6 E1 (2026-10-06 決定): 予備 A/B の U は評価・時間変動の不確かさ max(2T, 2E, E_exit) だけで、下限項を入れない。
# E3 の U の Δq/10 を入れると、閾値も Δq/10 なので「残留 IC 依存説を棄却」は差が厳密に 0 のときしか成立しない (plan の書き損じ)。
ICAB_U_FLOOR_FRAC = 0.0


def judge_icab_quantity(col: str, a, b, T, E, steady: bool, E_exit=0.0, delta: float | None = None,
                        u_floor_frac: float = ICAB_U_FLOOR_FRAC) -> dict:
    """予備 A/B の 1 量の判定。a = α (最近傍) の代表値、b = β (番号写像) の代表値 (各 1 本の窓平均)。
    U = max(2T, 2E, E_exit,q, u_floor_frac·Δq)。E3 と同じ項から、腕内 3 回の幅 3R_A・3R_B を除いた (各 1 本なので幅が無い)。
      negligible (残留 IC 依存なし): |b − a| + U ≤ Δq/10
      dependent  (IC 依存あり):      |b − a| − U > Δq/10
      それ以外 indeterminate (判別不能)。未定常は unsteady、非有限は missing (どちらも判別不能側)。"""
    if delta is None:
        dk = delta_key(col)
        if dk is None:
            raise ValueError(f"{col} は判定対象の量でない")
        delta = DELTA_Q[dk]
    if not (_finite(delta) and float(delta) > 0):
        raise ValueError(f"{col}: Δq = {delta!r} は正の有限値であること")
    if not (_finite(u_floor_frac) and float(u_floor_frac) >= 0):
        raise ValueError(f"u_floor_frac = {u_floor_frac!r} は非負の有限値であること")
    thr = float(delta) / 10.0
    row = {"col": col, "Delta": float(delta), "threshold": thr}
    vals = {"a": a, "b": b, "T": T, "E": E, "E_exit": E_exit}
    bad = [k for k, v in vals.items() if not _finite(v)]
    if bad:
        row.update(status="missing", reason=f"非有限・欠損 ({', '.join(bad)})")
        return row
    for nm in ("T", "E", "E_exit"):
        if float(vals[nm]) < 0:
            raise ValueError(f"{col}: {nm} = {vals[nm]} が負 (幅・感度は非負)")
    d = float(b) - float(a)
    terms = {"2T": 2 * float(T), "2E": 2 * float(E), "E_exit": float(E_exit), "floor": float(u_floor_frac) * float(delta)}
    U = max(terms.values())
    row.update(a=float(a), b=float(b), b_minus_a=d, T=float(T), E=float(E), E_exit=float(E_exit), U=U, U_terms=terms,
               note="U は E3 と同じ項から 3R_A・3R_B を除いたもの (α・β は各 1 本で腕内の幅が無い)")
    if not steady:
        row.update(status="unsteady", reason="窓の評価量が STEADY (または絶対許容内) でない (判別不能)")
        return row
    if abs(d) + U <= thr:
        row["status"] = "negligible"
    elif abs(d) - U > thr:
        row["status"] = "dependent"
    else:
        row["status"] = "indeterminate"
    return row


def icab_overall(rows, preconditions_ok: bool = True) -> dict:
    """予備 A/B の総合。前提不成立 → 判別不能。どれか dependent → 前提を棄却・形状採否を保留。全量 negligible → 残留 IC 依存説を棄却。
    それ以外 (未定常・境界域・欠損) は判別不能。いずれも「棄却」以外は諮問 (§6 E1)。"""
    st = {r["col"]: r["status"] for r in rows}
    counts = {}
    for s in st.values():
        counts[s] = counts.get(s, 0) + 1
    dep = [c for c, s in st.items() if s == "dependent"]
    if not preconditions_ok:
        v = "判別不能 (前提未達; 諮問)"
    elif not st:
        v = "判別不能 (判定する量が無い; 諮問)"
    elif dep:
        v = "前提を棄却・形状採否を保留 (IC の影響は記録だけでは不十分; 諮問)"
    elif all(s == "negligible" for s in st.values()):
        v = "残留 IC 依存説を棄却 (β で腕 B の r2・r3 へ進む)"
    else:
        v = "判別不能 (未定常・境界域・欠損; 諮問)"
    return {"verdict": v, "counts": counts, "dependent": dep,
            "not_negligible": [c for c, s in st.items() if s != "negligible"]}
