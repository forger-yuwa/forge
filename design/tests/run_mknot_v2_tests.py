#!/usr/bin/env python3
"""M_knot の規則 coarse_min_theta_v2 (`forge_design.evaluate.mknot_select.select_coarse_min_theta`) の試験 — plan
tooling-nozzle-mknot-fast §6 V1 (§5.1 #3)。

V1 (合成の評価、設計チェーンは回さない; (vii) だけ実際の問題):
  (i)   谷の取り違え (plan 段 M1): 粗い網 min(15 + 2(M−2)², 15.02 + 400(M−3.45)²)、生産の網で第 2 項の定数 14.5 → 3.45
        (粗い網だけで選ぶと 2.00)
  (ii)  ずらした後の比較 (plan 段 v2 M2): 粗い網は (i) と同じ、生産の網 min(15 + 400(M−2)², 15.02 + 400(M−3.45)²)、生産の網で 1.95 だけ
        不成立 → 3.45 (2.00 の谷は 2.05 にずれて 16.00°、それを返さない)
  (iii) 確認の分岐: 中心だけ不成立 → 谷を捨てる、両側不成立 → 捨てる、片側不成立 → 反対側へずらして再評価、移動先が区間の外 → 捨てる、
        移動 3 回の上限 → 捨てる (初回の確認 + 移動 3 回)、全部の谷を捨てる → 例外 (表つき)
  (iv)  候補の区間: 成立の端の点を候補にしない、M=[2.00, 2.05, 3.00, 3.05]・θ=[15, 16, 17, 18] (2 区間) で両方の区間の極小 (2.00 と
        3.00) を拾う、平らな列は 1 つの極小 (代表は推定の距離が最大、同じなら小さい M_K)、入力の順を逆にしても同じ結果
  (v)   同点: θ 15.000・15.009・15.018° で推定の距離 0.1・0.2・0.3 → 同点の集合は {15.000, 15.009} (15.018 は θ_min + 0.01 を超える)
        → 距離 0.2 の谷。比べる順 (6 通り) を変えても同じ。逐次比較なら順で結果が変わる (反例の確認)
  (vi)  粗い網で成立点なし・長さ合わせの全点未収束・候補なし → 例外 (表つき)
  (vii) --real のときだけ: 捕捉を 2 つに限っても、time_s を除く数値・判定の欄が全部を捕捉したとき (旧の挙動) と完全に一致 (実際の問題
        case/70 problem_m5_air_ns2_euler.yaml の 3 点、N 600 (L_total_tol 0.1 r_t = 粗い網) と N 2400 (カードの L_total_tol = 生産の網))。
        捕捉の外を採るように細工する (design_chain の後に逆 MOC を 2 回足す) と行の error・ok=False で、選択に使われない
  ほか: 格子 (整数 k、M_K = k/20、両端から 0.05 以上内側)。
  (6d)  長さ合わせの停止の判定 `v2_length_stall` (plan tooling-nozzle-card-local-mesh §5.1 #6d): 生産の網の未収束の行 (選べない例外・
        選べた結果の両方)・粗い網の全点未収束は停止、制約による不成立 (中心・成立点なし・候補なし) と全部成立は停止でない

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_mknot_v2_tests.py [--real]
所要: 数秒 (--real で 約 1.5 分: 実際の問題の評価 14 回; 1 プロセスで逐次、最大 RSS 約 1.3 GB)。
"""
import itertools
import math
import sys
import time
from pathlib import Path

import yaml

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import mknot_select as MS  # noqa: E402

FAIL = 0
REAL = "--real" in sys.argv[1:]


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""), flush=True)
    if not cond:
        FAIL += 1


def row(M, theta, *, margin=4.0, sin=0.05, mpp=1, error=None):
    """合成の評価の行 (evaluate_mknot と同じ鍵)。成立の判定は本物の feasibility_reasons で作る。"""
    if error:
        r = {"M_K": M, "error": error}
    else:
        r = {"M_K": M, "theta_max_deg": theta, "margin_deg": margin, "min_sin_angle": sin, "mpp_sign_changes": mpp,
             "mpp_quality_ok": mpp <= 1, "n_flip_interior": 0, "n_crossings_sampled": 0, "pileup_pairs_interior": 0,
             "qa_violations": [], "moc_gate_ok": True, "moc_gate_reason": None}
    r["reasons"] = MS.feasibility_reasons(r)
    r["ok"] = not r["reasons"]
    return r


def bad(M, theta=14.0):
    """不成立の行。格子の k が偶数なら θ_max は有限で小さいが余裕が床の下 (ok False; plan 段 v2 M1 の型)、奇数なら設計チェーンの例外
    (θ_max なし)。両方の型を混ぜる。"""
    if MS.v2_k(M) % 2:
        return row(M, None, error="ValueError: 壁 QA 不合格: スロート曲率不整合 (合成)")
    return row(M, theta, margin=0.5)


def near(a, b, tol=1e-9):
    return a is not None and b is not None and abs(a - b) <= tol


def sel(fc, fp, M_A=1.0, M_d=5.0, **kw):
    return MS.select_coarse_min_theta(fc, fp, M_A, M_d, **kw)


def ks(*Ms):
    return {MS.v2_k(m) for m in Ms}


def in_set(M, S):
    return MS.v2_k(M) in S


# ---------------------------------------------------------------------------------------------------------------------
print("== 格子 (整数 k、M_K = k/20、両端から 0.05 以上内側)")
g = MS.v2_grid(1.2293, 5.0)
check("M_A 1.2293・M_d 5.0 → 1.30〜4.95 (k 26〜99、74 点; 1.25 は M_A から 0.0207 で外、4.95 は M_d からちょうど 0.05 で内)",
      g[0] == 26 and g[-1] == 99 and len(g) == 74 and MS.v2_M(41) == 2.05 and MS.v2_k(2.05) == 41, f"{g[0]}〜{g[-1]} ({len(g)})")
try:
    MS.v2_k(2.051)
    check("格子の上に無い M_K → ValueError", False)
except ValueError:
    check("格子の上に無い M_K → ValueError", True)

print("== V1 (i) 谷の取り違え (plan 段 M1)")


def coarse_i(M):
    return row(M, min(15.0 + 2.0 * (M - 2.0) ** 2, 15.02 + 400.0 * (M - 3.45) ** 2))


def prod_i(M):
    return row(M, min(15.0 + 2.0 * (M - 2.0) ** 2, 14.5 + 400.0 * (M - 3.45) ** 2))


r = sel(coarse_i, prod_i)
coarse_best = min((x for x in r["coarse_table"] if MS._feasible(x)), key=lambda x: x["theta_max_deg"])
check("(i) 反例の確認: 粗い網だけの最小は 2.00 (谷を取り違える側)", near(coarse_best["M_K"], 2.00), f"{coarse_best['M_K']}")
check("(i) 谷を 2 つ拾い (2.00・3.45)、生産の網で両方確かめて 3.45 (θ 14.5°) を選ぶ",
      near(r["chosen"], 3.45) and near(r["theta_max_deg"], 14.5) and sorted(v["start"] for v in r["valleys"]) == [2.0, 3.45]
      and all(v["final"] is not None for v in r["valleys"]) and any("谷" in w for w in r["warnings"]),
      f"chosen {r['chosen']!r} 谷 {[(v['start'], v['final'], v['theta_max_deg']) for v in r['valleys']]}")
check("(i) 出す値は格子の値 (3.45 そのもの)・選んだ点と両隣の生産の網の行が成立",
      r["chosen"] == 3.45 and [x["M_K"] for x in r["final_three"]] == [3.4, 3.45, 3.5] and all(x["ok"] for x in r["final_three"]),
      str([x["M_K"] for x in r["final_three"]]))

print("== V1 (ii) ずらした後の比較 (plan 段 v2 M2)")


def prod_ii(M):
    if abs(M - 1.95) < 1e-9:
        return bad(M)
    return row(M, min(15.0 + 400.0 * (M - 2.0) ** 2, 15.02 + 400.0 * (M - 3.45) ** 2))


r = sel(coarse_i, prod_ii)
v200 = [v for v in r["valleys"] if near(v["start"], 2.00)][0]
check("(ii) 2.00 の谷は左隣 1.95 が不成立 → 2.05 へずらし、3 点成立で代表 2.05 (16.00°)",
      near(v200["final"], 2.05) and v200["n_moves"] == 1 and near(v200["theta_max_deg"], 16.0, 1e-9)
      and v200["steps"][0]["left_ok"] is False and "ずらす" in v200["steps"][0]["action"], str(v200["steps"]))
check("(ii) ずらした後に谷を比べて 3.45 (15.02°) を選ぶ (2.05・16.00° を返さない)",
      near(r["chosen"], 3.45) and near(r["theta_max_deg"], 15.02, 1e-9), f"chosen {r['chosen']!r} θ {r['theta_max_deg']!r}")

print("== V1 (iii) 確認の分岐")


def one_valley(lo=1.50, hi=2.50):
    """粗い網: [lo, hi] で成立、θ = 15 + (M − 2)² (谷 1 つ、代表 2.00)。"""
    def fn(M):
        if not (lo - 1e-9 <= M <= hi + 1e-9):
            return bad(M)
        return row(M, 15.0 + (M - 2.0) ** 2)
    return fn


def prod_with(bad_set, lo=1.0, hi=5.0):
    def fn(M):
        if in_set(M, bad_set) or not (lo - 1e-9 <= M <= hi + 1e-9):
            return bad(M)
        return row(M, 15.0 + (M - 2.0) ** 2)
    return fn


def expect_raise(name, fn, needle):
    try:
        fn()
        check(name, False, "例外なし")
        return None
    except MS.MknotSelectionError as e:
        check(name, needle in str(e) and len(e.table) > 0, str(e).splitlines()[0][:200])
        return e


e = expect_raise("(iii) 中心 2.00 だけ不成立 → 谷を捨てる → 谷が残らず例外 (表つき)",
                 lambda: sel(one_valley(), prod_with(ks(2.00))), "谷が無い")
check("(iii) 中心だけ不成立: 捨てた理由が中心の不成立", e is not None and "中心 2.00 が生産の網で不成立" in str(e), "")
e = expect_raise("(iii) 両隣 1.95・2.05 が不成立 → 谷を捨てる",
                 lambda: sel(one_valley(), prod_with(ks(1.95, 2.05))), "谷が無い")
check("(iii) 両側不成立: 捨てた理由", e is not None and "両隣が生産の網で不成立" in str(e), "")
calls = []


def prod_iii_c(M):
    calls.append(M)
    return prod_with(ks(1.95))(M)


r = sel(one_valley(), prod_iii_c)
v = r["valleys"][0]
check("(iii) 片側 (1.95) 不成立 → 反対側の 2.05 へずらして 2.10 も評価し、2.00・2.05・2.10 成立で代表 2.05",
      near(r["chosen"], 2.05) and v["n_moves"] == 1 and len(v["steps"]) == 2 and near(v["steps"][1]["center"], 2.05)
      and sorted(set(round(m, 9) for m in calls)) == [1.95, 2.0, 2.05, 2.1] and len(calls) == 4,
      f"chosen {r['chosen']!r} 評価 {sorted(set(calls))} ({len(calls)} 回)")
e = expect_raise("(iii) 移動先が区間の外: 粗い網の成立 [1.90, 2.05] (候補の区間 [1.95, 2.00]、代表 2.00) で 1.95 が不成立 → 2.05 は区間の外 → 捨てる",
                 lambda: sel(one_valley(1.90, 2.05), prod_with(ks(1.95))), "谷が無い")
check("(iii) 区間の外: 捨てた理由", e is not None and "候補の区間" in str(e) and "の外" in str(e), "")
# 上限: 生産の網の成立が {2.00, 2.05} だけ → 2.00 (左 ×) → 2.05 (右 ×) → 2.00 → 2.05 (移動 3 回) → 上限で捨てる
e = expect_raise("(iii) 移動 3 回の上限: 生産の網の成立が {2.00, 2.05} だけ → 2.00 → 2.05 → 2.00 → 2.05 で上限 → 捨てる",
                 lambda: sel(one_valley(), prod_with(set(), lo=2.00, hi=2.05)), "谷が無い")
check("(iii) 上限: 捨てた理由に上限 3 回", e is not None and "上限 3 回" in str(e), "")
vs_lim = MS.v2_confirm([{"interval_k": (31, 49), "plateau_k": [40], "rep_k": 40, "theta_coarse": 15.0}],
                       lambda Ms: [prod_with(set(), lo=2.00, hi=2.05)(m) | {"M_K": m} for m in Ms])
check("(iii) 上限: 確認は初回 + 移動 3 回の 4 回 (中心 2.00 → 2.05 → 2.00 → 2.05)",
      [s["center"] for s in vs_lim[0]["steps"]] == [2.0, 2.05, 2.0, 2.05] and vs_lim[0]["n_moves"] == 3 and vs_lim[0]["final"] is None,
      str([s["center"] for s in vs_lim[0]["steps"]]))


def two_valley_coarse(M):
    return row(M, min(15.0 + 2.0 * (M - 2.0) ** 2, 15.02 + 400.0 * (M - 3.45) ** 2))


e = expect_raise("(iii) 谷 2 つとも捨てる (両方の中心が不成立) → 例外 (表は粗い掃引と生産の網の行)",
                 lambda: sel(two_valley_coarse, lambda M: bad(M) if in_set(M, ks(2.00, 3.45)) else row(M, 15.0)), "谷が無い")
check("(iii) 全部を捨てた例外: 谷ごとの理由と生産の網の行 (confirm) を含む",
      e is not None and str(e).count("中心") >= 2 and any(x.get("role") == "confirm" for x in e.table), "")

print("== V1 (iv) 候補の区間")
feas_iv = ks(1.95, 2.00, 2.05, 2.10, 2.95, 3.00, 3.05, 3.10)
th_iv = {2.00: 15.0, 2.05: 16.0, 3.00: 17.0, 3.05: 18.0}


def coarse_iv(M):
    if not in_set(M, feas_iv):
        return bad(M)
    return row(M, th_iv.get(round(M, 2), 20.0))


rows_iv = [dict(coarse_iv(MS.v2_M(k)), M_K=MS.v2_M(k)) for k in range(30, 70)]
ana = MS.v2_coarse_analysis(rows_iv)
vals = MS.v2_valleys(ana)
check("(iv) 成立の端の点 (1.95・2.10・2.95・3.10) を候補にしない; 候補は 2.00・2.05・3.00・3.05 の 2 区間",
      [MS.v2_M(k) for k in ana["candidates_k"]] == [2.0, 2.05, 3.0, 3.05]
      and [(MS.v2_M(a), MS.v2_M(b)) for a, b in ana["intervals_k"]] == [(2.0, 2.05), (3.0, 3.05)],
      f"{[MS.v2_M(k) for k in ana['candidates_k']]} {ana['intervals_k']}")
check("(iv) 両方の区間の極小 (2.00 と 3.00; 区間の端も内側の隣と比べて極小なら含める) を拾う",
      [MS.v2_M(v["rep_k"]) for v in vals] == [2.0, 3.0], str([MS.v2_M(v["rep_k"]) for v in vals]))
for nm, perm in (("逆順", list(reversed(rows_iv))), ("入れ替え", rows_iv[1::2] + rows_iv[::2])):
    a2 = MS.v2_coarse_analysis(perm)
    v2 = MS.v2_valleys(a2)
    check(f"(iv) 入力の順を変えても同じ結果 ({nm})", a2["intervals_k"] == ana["intervals_k"] and v2 == vals, "")
# 平らな列: 成立 1.80〜2.25 (不成立の最寄り 1.75・2.30)、θ 1.95・2.00・2.05 が 15 (同じ値) → 極小は 1 つ。推定の距離は 1.95: 0.15、
# 2.00・2.05: 0.20 (同じ) → 小さい M_K の 2.00
feas_f = set(range(36, 46))
th_f = {39: 15.0, 40: 15.0, 41: 15.0}
rows_f = [dict(row(MS.v2_M(k), th_f.get(k, 16.0)) if k in feas_f else bad(MS.v2_M(k)), M_K=MS.v2_M(k)) for k in range(30, 52)]
af = MS.v2_coarse_analysis(rows_f)
vf = MS.v2_valleys(af)
check("(iv) 平らな列 (1.95・2.00・2.05 が同じ θ) は 1 つの極小、代表は推定の距離が最大の点で同じなら小さい M_K → 2.00",
      len(vf) == 1 and [MS.v2_M(k) for k in vf[0]["plateau_k"]] == [1.95, 2.0, 2.05] and MS.v2_M(vf[0]["rep_k"]) == 2.0
      and [MS.v2_distance(k, af["infeasible_k"]) for k in vf[0]["plateau_k"]] == [0.15, 0.2, 0.2],
      f"{[(MS.v2_M(v['rep_k']), [MS.v2_M(k) for k in v['plateau_k']]) for v in vf]} 距離 "
      f"{[MS.v2_distance(k, af['infeasible_k']) for k in (vf[0]['plateau_k'] if vf else [])]}")
check("(iv) 平らな列: 入力の順を逆にしても同じ", MS.v2_valleys(MS.v2_coarse_analysis(list(reversed(rows_f)))) == vf, "")
rows_e = [r if MS._feasible(r) else dict(row(r["M_K"], None, error="ValueError: 壁 QA 不合格"), M_K=r["M_K"]) for r in rows_iv]
check("(iv) 不成立の行が θ_max を持たない (設計チェーンの例外) ときも同じ区間・谷",
      MS.v2_valleys(MS.v2_coarse_analysis(rows_e)) == vals, "")

print("== V1 (v) 同点 (一度だけ作る同点の集合)")
# 谷 A (2.00、不成立の最寄りまで 3 格子 → 推定の距離 0.10)、B (2.50、5 格子 → 0.20)、C (3.20、7 格子 → 0.30)
feas_v = set(range(38, 43)) | set(range(46, 55)) | set(range(58, 71))
centers_v = {40: 15.000, 50: 15.009, 64: 15.018}


def nearest_c(k):
    return min(centers_v, key=lambda c: abs(c - k))


def coarse_v(M):
    k = MS.v2_k(M)
    if k not in feas_v:
        return bad(M)
    c = nearest_c(k)
    return row(M, 15.0 + 0.1 * (k - c) ** 2)


def prod_v(M):
    k = MS.v2_k(M)
    if k not in feas_v:
        return bad(M)
    c = nearest_c(k)
    return row(M, centers_v[c] + 0.1 * (k - c) ** 2)


r = sel(coarse_v, prod_v)
check("(v) 谷 3 つの代表 2.00・2.50・3.20、推定の距離 0.10・0.20・0.30、θ 15.000・15.009・15.018",
      [(v["final"], v["distance"], round(v["theta_max_deg"], 6)) for v in r["valleys"]]
      == [(2.0, 0.1, 15.0), (2.5, 0.2, 15.009), (3.2, 0.3, 15.018)], str([(v["final"], v["distance"], v["theta_max_deg"]) for v in r["valleys"]]))
check("(v) 同点の集合は {2.00, 2.50} (3.20 の 15.018 は θ_min + 0.01 を超える) → 距離 0.20 の 2.50",
      r["tie_set"] == [2.0, 2.5] and near(r["chosen"], 2.5) and near(r["distance"], 0.2), f"tie {r['tie_set']} chosen {r['chosen']!r}")
finals_v = [{"M_K": 2.0, "theta_max_deg": 15.000, "distance_steps": 2.0}, {"M_K": 2.5, "theta_max_deg": 15.009, "distance_steps": 4.0},
            {"M_K": 3.2, "theta_max_deg": 15.018, "distance_steps": 6.0}]
picks = {MS.v2_pick(list(pm))["chosen"]["M_K"] for pm in itertools.permutations(finals_v)}
check("(v) 比べる順 6 通りで全部 2.50", picks == {2.5}, str(picks))


def sequential(order):
    """逐次比較 (直す前の型): 今の勝者と次を比べ、差 ≤ 0.01° なら距離の大きい方、そうでなければ θ の小さい方。"""
    best = order[0]
    for f in order[1:]:
        if abs(f["theta_max_deg"] - best["theta_max_deg"]) <= 0.01:
            best = f if f["distance_steps"] > best["distance_steps"] else best
        elif f["theta_max_deg"] < best["theta_max_deg"]:
            best = f
    return best["M_K"]


seq = {tuple(f["M_K"] for f in pm): sequential(list(pm)) for pm in itertools.permutations(finals_v)}
check("(v) 反例の確認: 逐次比較なら比べる順で結果が変わる (2.00・2.50・3.20 のどれにもなる)", set(seq.values()) == {2.0, 2.5, 3.2}, str(seq))

print("== V1 (vi) 例外 (表つき)")
expect_raise("(vi) 粗い網で成立点なし → 例外", lambda: sel(lambda M: bad(M), lambda M: row(M, 15.0)), "成立点が無い")
UNC = "ValueError: Lc_mode: from_length 未収束 (|x_F − L_total| = 0.05 > tol 0.02)"
expect_raise("(vi) 粗い網の長さ合わせが全点で未収束 → 例外", lambda: sel(lambda M: row(M, None, error=UNC), lambda M: row(M, 15.0)),
             "全点で未収束")
expect_raise("(vi) 候補なし (成立点が 1 つおき = 両隣とも成立する点が無い) → 例外",
             lambda: sel(lambda M: row(M, 15.0) if MS.v2_k(M) % 2 == 0 else bad(M), lambda M: row(M, 15.0)), "候補が無い")
e = expect_raise("(vi) 一部だけ未収束で残りが不成立 → 「成立点が無い」(全点未収束ではない)",
                 lambda: sel(lambda M: row(M, None, error=UNC) if M < 2.0 else bad(M), lambda M: row(M, 15.0)), "成立点が無い")

print("== 長さ合わせの停止の判定 v2_length_stall (plan tooling-nozzle-card-local-mesh §5.1 #6d; §4.1 の拡張)")
UNCP = "ValueError: Lc_mode: from_length 未収束 (|x_F − L_total| = 0.0172 > tol 0.00476)"


def unc_row(M, res=-0.0172):
    r = row(M, None, error=UNCP)
    r.update(L_residual_rt=res, length_stop="from_length 未収束")
    return r


e = expect_raise("(6d) 谷 1 つの中心が生産の網で未収束 → 谷が残らず例外",
                 lambda: sel(one_valley(), lambda M: unc_row(M) if near(M, 2.00) else row(M, 15.0 + (M - 2.0) ** 2)), "谷が無い")
st = MS.v2_length_stall(error=e) if e is not None else None
check("(6d) その例外は長さの停止 (生産の網・M_K 2.00・谷 2.00 の確認で未収束・最良の残差 0.0172・選べた値なし)",
      st is not None and st["grid"] == "prod" and [r["M_K"] for r in st["rows"]] == [2.0] and st["best_abs_residual_rt"] == 0.0172
      and [(v["start"], v["unconverged"]) for v in st["valleys"]] == [(2.0, [2.0])] and st["selected"] is None
      and e.detail.get("stage") == "no_final_valley", str(st))
r2 = sel(two_valley_coarse, lambda M: unc_row(M) if near(M, 2.00) else row(M, 15.0 + min(2.0 * (M - 2.0) ** 2, 400.0 * (M - 3.45) ** 2)))
st = MS.v2_length_stall(result=r2)
check("(6d) 谷 2 つで 2.00 の中心だけ未収束: 選択は 3.45 で黙って成功する (警告だけ) が、判定は長さの停止を返す (選べた値 3.45)",
      near(r2["chosen"], 3.45) and st is not None and st["grid"] == "prod" and near(st["selected"], 3.45)
      and [v["unconverged"] for v in st["valleys"]] == [[2.0], []], str(st and {k: st[k] for k in ("selected", "valleys")}))
e = expect_raise("(6d) 中心が制約で不成立 (壁での余裕) → 谷が残らず例外", lambda: sel(one_valley(), prod_with(ks(2.00))), "谷が無い")
check("(6d) 制約による不成立は長さの停止でない (None)", e is not None and MS.v2_length_stall(error=e) is None, "")
r3 = sel(one_valley(), prod_with(set()))
check("(6d) 全部成立の選択は長さの停止でない (None)", MS.v2_length_stall(result=r3) is None, "")
e = expect_raise("(6d) 粗い網の長さ合わせが全点で未収束 → 例外", lambda: sel(lambda M: unc_row(M, 0.102), lambda M: row(M, 15.0)), "全点で未収束")
st = MS.v2_length_stall(error=e) if e is not None else None
check("(6d) 粗い網の全点未収束は長さの停止 (粗い網・全点の行・最良の残差 0.102)",
      st is not None and st["grid"] == "coarse" and len(st["rows"]) == len(e.table) and st["best_abs_residual_rt"] == 0.102
      and st["valleys"] == [] and e.detail.get("stage") == "coarse_length", str(st and st["grid"]))
e = expect_raise("(6d) 粗い網で成立点なし (制約) → 例外", lambda: sel(lambda M: bad(M), lambda M: row(M, 15.0)), "成立点が無い")
check("(6d) 粗い網で成立点なし・候補なしは長さの停止でない (None)", e is not None and MS.v2_length_stall(error=e) is None, "")
try:
    MS.v2_length_stall()
    check("(6d) result と error のどちらでもない → ValueError", False)
except ValueError:
    check("(6d) result と error のどちらでもない → ValueError", True)

# ---------------------------------------------------------------------------------------------------------------------
if REAL:
    print("== V1 (vii) 捕捉を 2 つに限る (実際の問題 case/70 problem_m5_air_ns2_euler.yaml、1 プロセスで逐次)")
    C = DESIGN.parent / "case" / "70.m5_air_wt_demo"
    raw = yaml.safe_load((C / "problem_m5_air_ns2_euler.yaml").read_text())
    PTS = (1.95, 2.05, 3.0)
    t0 = time.time()

    def ev(M, N, L, keep):
        old = MS._CAPTURE_KEEP
        MS._CAPTURE_KEEP = keep
        try:
            return MS.evaluate_mknot(raw, M, N, base_dir=C, L_total_tol=L)
        finally:
            MS._CAPTURE_KEEP = old

    def diff(a, b):
        ka = set(a) - {"time_s"}
        kb = set(b) - {"time_s"}
        out = sorted(k for k in ka | kb if (k not in a or k not in b or not _same(a[k], b[k])))
        return out

    def _same(x, y):
        if isinstance(x, float) and isinstance(y, float):
            return (math.isnan(x) and math.isnan(y)) or x == y
        return x == y

    for N, L, tag in ((600, MS.V2_L_COARSE, "N 600・L_total_tol 0.1 r_t (粗い網)"), (2400, None, "N 2400・カードの L_total_tol (生産の網)")):
        for M in PTS:
            a, b = ev(M, N, L, None), ev(M, N, L, 2)
            d = diff(a, b)
            check(f"(vii) {tag}・M_K {M}: 捕捉 2 と全部で time_s を除く全部の欄が一致 (成立 {b['ok']}、θ {b.get('theta_max_deg')})",
                  not d, f"違う欄 {d}; 停止 {b.get('length_stop')}")
    # 細工: design_chain の後に逆 MOC を 2 回足して、design_chain の壁を捕捉の外に押し出す
    import forge_design.evaluate.runner_axismach as RA
    real_chain = RA.design_chain

    def rigged(p):
        inner = RA.inverse_design                   # evaluate_mknot の捕捉の包み
        last = {}

        def rec(*a, **k):
            last["a"], last["k"] = a, k
            return inner(*a, **k)
        RA.inverse_design = rec
        try:
            d = real_chain(p)
        finally:
            RA.inverse_design = inner
        inner(*last["a"], **last["k"])
        inner(*last["a"], **last["k"])
        return d

    RA.design_chain = rigged
    try:
        rr = ev(2.05, 600, MS.V2_L_COARSE, 2)
        rr_all = ev(2.05, 600, MS.V2_L_COARSE, None)
    finally:
        RA.design_chain = real_chain
    check("(vii) 細工 (捕捉の外の壁を採る): 捕捉 2 では行の error・ok=False (全部の捕捉なら成立)",
          rr["ok"] is False and "捕捉できない" in (rr.get("error") or "") and rr_all["ok"] is True,
          f"error {rr.get('error')!r} / 全部の捕捉 ok {rr_all['ok']}")

    def coarse_rig(M):
        return dict(rr, M_K=M) if near(M, 2.05) else row(M, 15.0 + (M - 2.05) ** 2)

    rs = sel(coarse_rig, lambda M: row(M, 15.0 + (M - 2.10) ** 2), M_A=1.5, M_d=3.0)
    check("(vii) その行は選択に使われない (粗い網で不成立 → 2.05 は候補にも谷にもならない)",
          not near(rs["chosen"], 2.05) and 2.05 not in rs["candidates"] and all(not near(v["start"], 2.05) for v in rs["valleys"]),
          f"chosen {rs['chosen']!r} 候補に 2.05 {2.05 in rs['candidates']}")
    print(f"   (vii) 所要 {time.time() - t0:.0f} s")
else:
    print("== V1 (vii) は --real のときだけ (飛ばした)")

print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
