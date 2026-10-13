#!/usr/bin/env python3
"""M_knot の選び方 `forge_design.evaluate.mknot_select` の試験 — plan tooling-nozzle-mknot-auto §6 V2・V1b (§5.1 #2・#2b)。

V2 (合成の評価、設計チェーンは回さない):
  (a) 単峰 → 最小
  (b) 局所極小 2 つ (粗い格子では浅い谷が勝つ形; plan 段レビュー M3 の反例を θ_max の最小化に直したもの) → 全部詰めて大域の最小 + 警告
  (c) 最小との差 1e-3° 以内の並び → 境目から最も遠い点
  (d) 成立点なし → 例外 (表つき)
  (e) 制約のどれか 1 つを破る点 (M″ の符号反転・網の sin 角・余裕・反転・交差・重なり・壁 QA・例外) は選ばれない。
      制約は全部を数える (値が無いものは判定不能 = 不成立)
  (f) 最小が境目から 0.05 未満 → 境目から 0.05 以上の候補の中の最小 (同点の規則の後; この入力では境目 + 0.05 の点)。
      成立の範囲が 0.1 より狭い → 例外。境目 + 0.05 の点が不成立 → 選び直し
  (g) 解像度のゲート: 合格の例、選択が 0.01 を超えて動く・高い解像度だけ制約違反・最小が詰めの端 → 例外 (選択は未確定)
  (h) 選ぶ順序の反例 (諮問 mknot-result-v6): 境目 1.98・ΔM_safe 0.05・tie 1e-3°、θ_max {2.00: 15.0000 (A) / 15.0020 (B)、
      2.04: 15.0008、2.05: 15.0015} で A・B とも 2.05 (距離で除外される 2.00 の値が選択を変えない)
  (i) 2N の境目が保守側へ動く合成の入力: N と 2N の境目を全部残し (§5.1 #2c で対にしてまとめるのをやめた)、全部の境目で
      N と 2N を選び直して比べる (N の境目のままで判定しない)。逆向き (2N が保守的でない側へ動く) では効く境目は N の境目
  (j) MOC のゲート (require_moc_gate と同じ関数で判定) だけを破る点は選ばれない。fixed2 の pass = None は今の扱いのまま通す
  (k) 境目をまとめない (諮問 mknot-gate-output の反例): N の不成立の区間 (2.00, 2.02)・2N の不成立の区間 (2.08, 2.10) で、2.05 は
      両方の解像度で成立するが、全部の境目からの距離 0.03 < 0.05 で選ばれない (対にしてまとめると 0.05 で通ってしまう)。
      ゲートの row_chosen_N・row_chosen_2N は選び直した点の行 (入力の点の行は row_input_N・row_input_2N)
V1b (M6 の目的の回帰): case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv の L_c 45・50 (n_axis 600) の記録を評価に差し替え、
  境目の距離と同点の規則を切った設定 (ΔM_safe 0、tie 0) で L_c 45 → 2.6・50 → 2.7 を選ぶ (余裕最大の 3.0・3.1 を選んだら取り違え)。
  既定の設定 (tie 1e-3°・ΔM_safe 0.05) の選択値も記録する (合否なし)。CSV は 0.1 刻みなので格子は coarse = fine = 0.1、原点は
  CSV の 1.5 に合わせて M_A = 1.4 (M_d = 6)。CSV に無い点 (0.1 刻みの外・二分法の点・4.0 より上) は「記録なし」= 不成立 (保守側)。
  CSV に無い列: n_crossings_sampled・pileup_pairs_interior・壁 QA の違反は hard_gate (compare_axislaw_D.py の hard) に畳まれている
  → hard_gate True なら 0 / 空、False なら不成立の理由に「hard_gate False」を足す。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_mknot_select_tests.py
所要: 数秒 (設計チェーンは回さない)。
"""
import csv
import sys
from pathlib import Path

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import mknot_select as MS  # noqa: E402

FAIL = 0
M6_CSV = DESIGN.parent / "case" / "42.isobutane_wt" / "sweep_axislaw_A_MK_summary.csv"


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""), flush=True)
    if not cond:
        FAIL += 1


def row(M, theta, *, margin=4.0, sin=0.05, mpp=1, flip=0, cross=0, pile=0, qa=(), error=None, moc=None):
    """合成の評価の行 (evaluate_mknot と同じ鍵)。成立の判定は本物の feasibility_reasons で作る。
    moc: MOC のゲートの判定 (moc_gate_result の戻り)。None なら成立 ({moc_gate_ok: True})。"""
    if error:
        r = {"M_K": M, "error": error}
    else:
        r = {"M_K": M, "theta_max_deg": theta, "margin_deg": margin, "min_sin_angle": sin, "mpp_sign_changes": mpp,
             "mpp_quality_ok": mpp <= 1, "n_flip_interior": flip, "n_crossings_sampled": cross, "pileup_pairs_interior": pile,
             "qa_violations": list(qa), **(moc or {"moc_gate_ok": True, "moc_gate_reason": None})}
    r["reasons"] = MS.feasibility_reasons(r)
    r["ok"] = not r["reasons"]
    return r


def infeasible(M, why="余裕"):
    return row(M, 15.0, margin=0.5) if why == "余裕" else row(M, 15.0, sin=0.001)


def sel(fn, **kw):
    kw.setdefault("workers", 1)
    return MS.select_mknot(fn, kw.pop("M_A", 1.0), kw.pop("M_d", 5.0), **kw)


def near(a, b, tol=1e-9):
    return a is not None and b is not None and abs(a - b) <= tol


# ---------------------------------------------------------------------------------------------------------------------
print("== V2 (a) 単峰 → 最小")
r = sel(lambda M: row(M, 15.0 + (M - 2.53) ** 2))
check("(a) 単峰の最小 2.53 を選ぶ", near(r["chosen"], 2.53), f"chosen {r['chosen']!r}")
check("(a) 局所極小 1 個・警告なし・境目なし", len(r["local_minima"]) == 1 and not r["warnings"] and not r["boundaries"],
      f"{r['local_minima']} {r['warnings']}")
check("(a) 粗い格子の原点は M_A + 0.1 (1.1)、端は M_d − 0.1 (4.9)",
      near(min(x["M_K"] for x in r["table"] if x["role"] == "coarse"), 1.1)
      and near(max(x["M_K"] for x in r["table"] if x["role"] == "coarse"), 4.9))

print("== V2 (b) 局所極小 2 つ (粗い格子では浅い谷が勝つ)")


def two_valleys(M):
    return row(M, min(15.0 + 2.0 * (M - 2.0) ** 2, 14.5 + 400.0 * (M - 3.45) ** 2))


r = sel(two_valleys)
coarse = [x for x in r["table"] if x["role"] == "coarse"]
cw = min(coarse, key=lambda x: x["theta_max_deg"])
check("(b) 反例の確認: 粗い格子の最小は浅い谷 (2.0) の側", near(cw["M_K"], 2.0), f"{cw['M_K']} θ {cw['theta_max_deg']}")
check("(b) 全部の局所極小を詰めて大域の最小 3.45 を選ぶ", near(r["chosen"], 3.45) and near(r["chosen_row"]["theta_max_deg"], 14.5, 1e-9),
      f"chosen {r['chosen']!r}")
check("(b) 局所極小 2 つ以上の警告", len(r["local_minima"]) >= 2 and any("局所極小" in w for w in r["warnings"]), str(r["warnings"]))

print("== V2 (c) 同点 (差 ≤ 1e-3°) → 境目から最も遠い点")


def flat(M):
    if not (1.95 <= M <= 3.5):
        return infeasible(M)
    if abs(M - 2.5) <= 0.105:
        return row(M, 15.0 + 1e-4 * abs(M - 2.45))
    return row(M, 15.01 + (abs(M - 2.5) - 0.1))


r = sel(flat)
check("(c) 境目は下 1.95・上 3.5 (幅 ≤ 0.001)", len(r["boundaries"]) == 2
      and all(b["width"] <= 1e-3 + 1e-12 for b in r["boundaries"])
      and abs(r["boundaries"][0]["feasible_end"] - 1.95) <= 1e-3 and abs(r["boundaries"][1]["feasible_end"] - 3.5) <= 1e-3,
      str([(b["kind"], round(b["feasible_end"], 6), round(b["width"], 6)) for b in r["boundaries"]]))
check("(c) 素直な最小は 2.45", near(r["theta_min"]["M_K"], 2.45), str(r["theta_min"]))
check("(c) 同点の集合 (2.40〜2.60) から境目から最も遠い 2.60 を選ぶ", near(r["chosen"], 2.60) and "同点" in r["reason"],
      f"chosen {r['chosen']!r} tie {[round(m, 4) for m in r['tie_set']]}")
r0 = sel(flat, tie_deg=0.0)
check("(c) tie 0 なら素直な最小 2.45", near(r0["chosen"], 2.45), f"chosen {r0['chosen']!r}")

print("== V2 (d) 成立点なし → 例外 (表つき)")
try:
    sel(lambda M: infeasible(M))
    check("(d) 例外が出る", False)
except MS.MknotSelectionError as e:
    check("(d) 例外 (MknotSelectionError) に全評価点の表", len(e.table) == 39 and "成立点が無い" in str(e) and "余裕" in str(e),
          f"{len(e.table)} 行")

print("== V2 (e) 制約のどれか 1 つを破る点は選ばれない")
cases = {"M″ の符号反転": dict(mpp=2), "網の sin 角": dict(sin=0.019), "余裕": dict(margin=0.99), "網の反転": dict(flip=1),
         "交差": dict(cross=1), "重なり": dict(pile=1), "壁 QA": dict(qa=("スロート曲率不整合",)), "例外": dict(error="ValueError: x")}
needle = {"M″ の符号反転": "M″", "網の sin 角": "sin 角", "余裕": "余裕", "網の反転": "反転", "交差": "交差", "重なり": "重なり",
          "壁 QA": "壁 QA", "例外": "例外"}
for name, bad in cases.items():
    def fn(M, bad=bad):
        if abs(M - 2.5) < 1e-6:
            return row(M, 14.0, **bad)          # θ_max は一番小さいが制約を 1 つ破る
        return row(M, 15.0 + (M - 2.5) ** 2)
    r = sel(fn, tie_deg=0.0, dM_safe=0.0)
    badrow = [x for x in r["table"] if abs(x["M_K"] - 2.5) < 1e-9][0]
    r2 = sel(fn)
    check(f"(e) {name} を破る点 (2.5) は選ばれない", not near(r["chosen"], 2.5) and abs(r["chosen"] - 2.5) <= 0.0100001
          and badrow["ok"] is False and any(needle[name] in s for s in badrow["reasons"]),
          f"chosen {r['chosen']!r} 理由 {badrow['reasons']}")
    check(f"(e) {name}: 既定の設定では不成立の点を境目として 0.05 以上離れる",
          MS.boundary_distance(r2["chosen"], r2["boundaries"]) >= 0.05 - 1e-9 and abs(r2["chosen"] - 2.5) >= 0.05 - 1e-9,
          f"chosen {r2['chosen']!r}")
multi = row(2.0, 15.0, margin=0.5, sin=0.01, mpp=3)
check("(e) 制約は全部を数える (余裕・sin 角・M″ の 3 つを同時に挙げる)", len(multi["reasons"]) == 3, str(multi["reasons"]))
partial = {"M_K": 2.0, "theta_max_deg": 15.0, "margin_deg": 4.0, "mpp_quality_ok": True}
pr = MS.feasibility_reasons(partial)
check("(e) 値が無い制約は判定不能 = 不成立 (例外が出なかっただけでは成立にしない)",
      len(pr) >= 4 and all("判定不能" in s for s in pr if "sin" in s or "反転" in s or "交差" in s or "重なり" in s), str(pr))

print("== V2 (f) 最小が境目から 0.05 未満 → 境目から 0.05 以上の候補の中の最小 (同点の規則の後)")


def near_edge(M, lo=2.0037):
    if M < lo:
        return infeasible(M)
    return row(M, 15.0 + 10.0 * (M - 2.02) ** 2)


r = sel(near_edge)
b = [x for x in r["boundaries"] if x["kind"] == "lower"][0]
check("(f) 境目の区間は 2.0037 を挟み幅 ≤ 0.001", b["lo"] < 2.0037 <= b["hi"] and b["width"] <= 1e-3 + 1e-12,
      f"[{b['lo']:.6f}, {b['hi']:.6f}]")
check("(f) 素直な最小 2.02 は境目から 0.05 未満", near(r["theta_min"]["M_K"], 2.02) and r["theta_min"]["distance"] < 0.05,
      str(r["theta_min"]))
adm = [x for x in r["table"] if x.get("role") in MS.CANDIDATE_ROLES and MS._feasible(x)
       and MS.boundary_distance(x["M_K"], r["boundaries"]) >= 0.05 - 1e-9]
amin = min(adm, key=lambda x: x["theta_max_deg"])
check("(f) 境目から 0.05 以上の候補の中の最小 (同点の規則の後) を選ぶ — この入力では境目 (成立側の端) + 0.05 の点 (評価して成立)",
      near(r["chosen"], amin["M_K"]) and near(r["theta_min_admissible"]["M_K"], amin["M_K"])
      and near(r["chosen"], b["feasible_end"] + 0.05, 1e-12)
      and r["constrained"] and r["chosen_row"]["ok"] and r["chosen_row"]["role"] == "safe",
      f"chosen {r['chosen']!r} 境目 {b['feasible_end']!r} 理由 {r['reason']}")
try:
    sel(lambda M: near_edge(M, lo=1.97) if M <= 2.05 else infeasible(M))     # 成立 [1.97, 2.05] (幅 0.08、粗い格子の 2.0 を含む)
    check("(f) 成立の範囲が 2 ΔM_safe より狭い → 例外", False)
except MS.MknotSelectionError as e:
    check("(f) 成立の範囲が 2 ΔM_safe より狭い → 例外", "ΔM_safe" in str(e), str(e).splitlines()[0][:160])


def island(M):
    if 2.05 <= M <= 2.06:
        return infeasible(M, "sin")
    return near_edge(M)


r = sel(island)
dmin = MS.boundary_distance(r["chosen"], r["boundaries"])
check("(f) 境目 + 0.05 の点が不成立の島に落ちる → 島の境目からも 0.05 以上離れた成立点を選ぶ",
      r["chosen_row"]["ok"] and dmin >= 0.05 - 1e-9 and r["chosen"] > 2.06 + 0.05 - 1e-3,
      f"chosen {r['chosen']!r} 距離 {dmin:.4f} 境目 {[round(x['feasible_end'], 5) for x in r['boundaries']]}")

print("== V2 (g) 解像度のゲート")


def quad(c, lo=2.0, hi=3.0, k=100.0, sin_hole=None):
    def fn(M):
        if not (lo <= M <= hi):
            return infeasible(M)
        if sin_hole is not None and abs(M - sin_hole) < 0.015:
            return row(M, 15.0 + k * (M - c) ** 2, sin=0.015)
        return row(M, 15.0 + k * (M - c) ** 2)
    return fn


def gate(fN, f2N, **kw):
    s = sel(fN)
    return s, MS.resolution_gate(None, s["chosen"], 600, select_result=s, workers=(1, 1),
                                 eval_factory=lambda n: fN if n == 600 else f2N, **kw)


s, g = gate(quad(2.5), quad(2.503))
check("(g) 合格の例: 選択値 2.5 → 2N でも 2.5、全部の判定が成立", g["pass"] and near(g["chosen_2N"], 2.5) and near(s["chosen"], 2.5),
      f"{g['checks']} chosen_2N {g['chosen_2N']!r}")
for name, f2N, key, needle_ in (("選択が 0.01 を超えて動く", quad(2.53), "change", "動く"),
                                ("高い解像度だけ制約違反", quad(2.5, sin_hole=2.5), "ok_2N", "2N")):
    try:
        gate(quad(2.5), f2N)
        check(f"(g) {name} → 例外", False)
    except MS.MknotUndeterminedError as e:
        check(f"(g) {name} → 例外 (選択は未確定)", e.result["checks"][key] is False and needle_ in str(e) and "未確定" in str(e),
              f"{e.result['checks']}")


def drift(M):
    if not (2.0 <= M <= 3.0):
        return infeasible(M)
    return row(M, 15.0 + 1e-4 * (2.6 - M))       # 窓 (2.4〜2.6) の中で単調に下がる・幅 2e-5° (全部同点)


try:
    gate(quad(2.5), drift)
    check("(g) 最小が詰めの端 → 例外", False)
except MS.MknotUndeterminedError as e:
    c = e.result["checks"]
    check("(g) 最小が詰めの端 → 例外 (選択値は同点の規則で動かない; 端の判定だけが落ちる)",
          c["edge_2N"] is False and c["change"] and c["ok_N"] and c["ok_2N"] and c["distance"] and c["edge_N"]
          and e.result["min_offset_2N"] == 10, f"{c} offset {e.result['min_offset_2N']}")
s, g = gate(quad(2.5), quad(2.53), raise_on_fail=False)
check("(g) raise_on_fail=False なら結果の dict (pass False・理由)", g["pass"] is False and g["reasons"], str(g["reasons"])[:160])

print("== V2 (h) 選ぶ順序の反例 (諮問 mknot-result-v6): 距離で除外される点の θ_max が選択を変えない")
bnd_h = [{"lo": 1.979, "hi": 1.98, "width": 0.001, "kind": "lower", "feasible_end": 1.98, "reasons_infeasible": []}]
for tag, t200 in (("A", 15.0000), ("B", 15.0020)):
    rows_h = [dict(row(m, t), role="fine") for m, t in ((2.00, t200), (2.04, 15.0008), (2.05, 15.0015))]
    old_first = MS._pick([x for x in rows_h], bnd_h, 1e-3)[0]["M_K"]      # 先に全体で同点の集合を作る順序 (直す前)
    check(f"(h) {tag}: 反例の確認 — 先に全体で同点の集合を作る順序なら {'2.04' if tag == 'A' else '2.05'}",
          near(old_first, 2.04 if tag == "A" else 2.05), f"{old_first!r}")
    rh = MS.apply_rule(rows_h, bnd_h, tie_deg=1e-3, dM_safe=0.05)
    check(f"(h) {tag}: 2.00 の θ_max {t200:.4f} でも 2.05 を選ぶ (制約と距離 → 最小 → 同点)",
          near(rh["chosen"], 2.05) and near(rh["theta_min_admissible"]["M_K"], 2.04) and sorted(rh["tie_set"]) == [2.04, 2.05]
          and rh["excluded_by_distance"] == [2.00], f"chosen {rh['chosen']!r} 理由 {rh['reason']}")
    check(f"(h) {tag}: 記録の素の最小 (距離の前) は {'2.00 (距離で除外)' if tag == 'A' else '2.04'}",
          near(rh["theta_min"]["M_K"], 2.00 if tag == "A" else 2.04)
          and rh["theta_min"]["excluded_by_distance"] is (tag == "A") and rh["constrained"] is (tag == "A"), str(rh["theta_min"]))
    # 境目 + ΔM_safe の点 (2.03) を評価して足す: θ_max が同点の外なら選ばれない (足すが必ず選ぶわけではない)、最小なら選ぶ
    for t203, want in ((15.0030, 2.05), (14.9, 2.03)):
        rs_ = MS.apply_rule(rows_h, bnd_h, tie_deg=1e-3, dM_safe=0.05, safe_eval=lambda Ms, t=t203: [dict(row(m, t), role="safe") for m in Ms])
        check(f"(h) {tag}: 境目 + 0.05 の点 2.03 (θ_max {t203}) を評価して足す → {want}",
              [p["M_K"] for p in rs_["safe_points"]] == [2.03] and near(rs_["chosen"], want), f"chosen {rs_['chosen']!r}")

print("== V2 (i) 2N の境目が保守側へ動く → 全部の境目 (N と 2N、まとめない) で N と 2N を選び直して比べる")


def edge_fn(lo, hi=3.0):
    def fn(M):
        if not (lo <= M <= hi):
            return infeasible(M)
        return row(M, 15.0 + (M - 1.9) ** 2)            # 境目へ向かって下がる (最小は境目 + ΔM_safe の点)
    return fn


def lower_by_source(g):
    return {x["source"]: x for x in g["common_boundaries"] if x["kind"] == "lower"}


fN_i, f2N_i = edge_fn(2.0005), edge_fn(2.0205)        # 2N の下の境目が 0.02 上 (成立の範囲を狭める側) へ動く
s, g = gate(fN_i, f2N_i, raise_on_fail=False)
lowN = [x for x in s["boundaries"] if x["kind"] == "lower"][0]
cl = lower_by_source(g)
binding = max(x["feasible_end"] for x in cl.values())
check("(i) 下の境目は N (≈ 2.0005) と 2N (≈ 2.0205) の両方を残す (まとめない; N から見た 2N の同じ向き ≈ +0.02)",
      sorted(cl) == ["2N", "N"] and len([x for x in g["common_boundaries"] if x["kind"] == "lower"]) == 2
      and 2.0005 - 1e-12 <= cl["N"]["feasible_end"] <= 2.0015 and 2.0205 - 1e-12 <= cl["2N"]["feasible_end"] <= 2.0215
      and 0.0195 <= cl["N"]["nearest_other"]["offset"] <= 0.0205 and binding == cl["2N"]["feasible_end"],
      str([(x["kind"], round(x["feasible_end"], 6), x["source"], x["nearest_other"]) for x in g["common_boundaries"]]))
check("(i) N の選択は N の境目 + 0.05 (入力)", near(s["chosen"], lowN["feasible_end"] + 0.05, 1e-12), f"{s['chosen']!r}")
check("(i) 全部の境目で N と 2N をそれぞれ選び直すと両方 効く境目 (2N の端) + 0.05 (差 0) で合格",
      g["pass"] and near(g["chosen_N"], g["chosen_2N"]) and near(g["chosen_N"], binding + 0.05, 1e-9),
      f"N {g['chosen_N']!r} 2N {g['chosen_2N']!r} {g['checks']} {g['reasons']}")
check("(i) N の選び直しは入力から +0.02 動き、記録に残る (入力の点は全部の境目から 0.05 未満)",
      g["moved_N"] is not None and 0.0195 <= g["moved_N"] <= 0.0205 and g["distance_input"] < 0.05,
      f"moved {g['moved_N']!r} distance_input {g['distance_input']!r}")
check("(i) row_chosen_N・row_chosen_2N は選び直した点の行、row_input_N は入力の点の行",
      near(g["row_chosen_N"]["M_K"], g["chosen_N"]) and near(g["row_chosen_2N"]["M_K"], g["chosen_2N"])
      and near(g["row_input_N"]["M_K"], s["chosen"]) and near(g["row_input_2N"]["M_K"], s["chosen"]) and not near(s["chosen"], g["chosen_N"]),
      f"row_chosen_N {g['row_chosen_N']['M_K']!r} row_input_N {g['row_input_N']['M_K']!r}")
own2N = sel(f2N_i)["chosen"]
check("(i) N の境目のままで比べると (N 単独 vs 2N 単独の選択) 差 > 0.01 — 全部の境目で選び直さないと判定が変わる入力",
      abs(own2N - s["chosen"]) > 0.01, f"N {s['chosen']!r} 2N 単独 {own2N!r}")
s, g = gate(edge_fn(2.0205), edge_fn(2.0005), raise_on_fail=False)      # 逆向き: 2N が保守的でない側へ動く
cl = lower_by_source(g)
check("(i) 逆向き (2N の境目が下へ): 両方残し、効く境目は N の端、選び直しも入力と同じで合格",
      g["pass"] and sorted(cl) == ["2N", "N"] and max(cl.values(), key=lambda x: x["feasible_end"])["source"] == "N"
      and near(g["chosen_N"], s["chosen"]) and near(g["chosen_2N"], s["chosen"]) and g["moved_N"] == 0.0,
      f"{[(k, round(v['feasible_end'], 6)) for k, v in cl.items()]} N {g['chosen_N']!r} 2N {g['chosen_2N']!r}")

print("== V2 (j) MOC のゲート (require_moc_gate と同じ関数) だけを破る点は選ばれない")
from types import SimpleNamespace  # noqa: E402


def moc(corrector, gate_):
    p = SimpleNamespace(geometry={} if corrector is None else {"moc_corrector": corrector})
    return MS.moc_gate_result(p, {"moc": ({} if gate_ is None else {"gate": gate_})})


bad_gate = {"applicable": True, "pass": False, "reasons": ["単位過程の反復失敗 2 対"]}
cases_j = (("converge・pass True", moc("converge", {"applicable": True, "pass": True}), True, None),
           ("converge・pass False", moc("converge", bad_gate), False, "不合格"),
           ("converge・診断なし", moc("converge", None), False, "診断が無い"),
           ("キーなし (既定 converge)・pass None", moc(None, {"applicable": True, "pass": None}), False, "不合格"),
           ("fixed2・pass None (今の扱いのまま通す)", moc("fixed2", {"applicable": True, "pass": None}), True, None))
for name, res, want, needle_ in cases_j:
    check(f"(j) moc_gate_result: {name} → {'成立' if want else '不成立'}",
          res["moc_gate_ok"] is want and (needle_ is None or needle_ in (res["moc_gate_reason"] or "")), str(res)[:160])
for name, mres, chosen_ok in (("pass False", moc("converge", bad_gate), False), ("fixed2 の pass None", moc("fixed2", {"applicable": True, "pass": None}), True)):
    def fn_j(M, mres=mres):
        if abs(M - 2.5) < 1e-6:
            return row(M, 14.0, moc=mres)               # θ_max は一番小さく、MOC のゲート以外の制約は全部成立
        return row(M, 15.0 + (M - 2.5) ** 2)
    r = sel(fn_j, tie_deg=0.0, dM_safe=0.0)
    badrow = [x for x in r["table"] if abs(x["M_K"] - 2.5) < 1e-9][0]
    if chosen_ok:
        check(f"(j) {name} の点 (2.5) は成立して選ばれる (今の扱いのまま)", near(r["chosen"], 2.5) and badrow["ok"], f"chosen {r['chosen']!r}")
    else:
        check(f"(j) MOC のゲートだけを破る点 (2.5) は選ばれない (理由はゲートの 1 つだけ)",
              not near(r["chosen"], 2.5) and badrow["ok"] is False and len(badrow["reasons"]) == 1
              and "MOC" in badrow["reasons"][0] and "反復失敗" in badrow["reasons"][0],
              f"chosen {r['chosen']!r} 理由 {badrow['reasons']}")
nomoc = {k: v for k, v in row(2.0, 15.0).items() if not k.startswith("moc_gate")}
check("(j) MOC のゲートの値が無い行は判定不能 = 不成立", any("MOC" in x and "判定不能" in x for x in MS.feasibility_reasons(nomoc)),
      str(MS.feasibility_reasons(nomoc)))

print("== V2 (k) 境目をまとめない (諮問 mknot-gate-output の反例)")


def gap_fn(a, b):
    """[1.9, 2.3] で成立、ただし開区間 (a, b) は不成立。θ_max = 15 + (M − 2.05)² (最小は 2.05)。"""
    def fn(M):
        if not (1.9 <= M <= 2.3) or a < M < b:
            return infeasible(M)
        return row(M, 15.0 + (M - 2.05) ** 2)
    return fn


def bnd(kind, end):
    lo, hi = (end, end + 5e-4) if kind == "upper" else (end - 5e-4, end)
    return {"lo": lo, "hi": hi, "width": 5e-4, "kind": kind, "feasible_end": end, "reasons_infeasible": []}


bN_k, b2_k = [bnd("upper", 2.00), bnd("lower", 2.02)], [bnd("upper", 2.08), bnd("lower", 2.10)]
paired_k = [bnd("upper", min(2.00, 2.08)), bnd("lower", max(2.02, 2.10))]       # 同じ向きを対にして保守側を採る (直す前)
all_k = MS.conservative_boundaries(bN_k, b2_k)
check("(k) 反例の確認: 対にしてまとめた境目 (上 2.00・下 2.10) からは 2.05 の距離 0.05 (ΔM_safe を満たしてしまう)",
      near(MS.boundary_distance(2.05, paired_k), 0.05), f"{MS.boundary_distance(2.05, paired_k)!r}")
check("(k) N・2N の境目を全部 (4 個) 残し、2.05 の距離は 0.03 < 0.05",
      [(x["kind"], x["feasible_end"], x["source"]) for x in all_k]
      == [("upper", 2.00, "N"), ("lower", 2.02, "N"), ("upper", 2.08, "2N"), ("lower", 2.10, "2N")]
      and near(MS.boundary_distance(2.05, all_k), 0.03), f"{[(x['kind'], x['feasible_end'], x['source']) for x in all_k]}")
rows_k = [dict(row(m, 15.0 + (m - 2.05) ** 2), role="fine") for m in (1.95, 2.03, 2.04, 2.05, 2.06, 2.07, 2.15)]
rk = MS.apply_rule(rows_k, all_k, tie_deg=1e-3, dM_safe=0.05)
check("(k) apply_rule: θ_max の最小 2.05 は距離で除外され選ばれない (選ぶのは全部の境目から 0.05 以上の点)",
      2.05 in rk["excluded_by_distance"] and not near(rk["chosen"], 2.05) and rk["distance"] >= 0.05 - 1e-9,
      f"chosen {rk['chosen']!r} 距離 {rk['distance']!r} 除外 {rk['excluded_by_distance']}")
g = MS.resolution_gate(None, 2.05, 600, workers=(1, 1), eval_factory=lambda n: gap_fn(2.00, 2.02) if n == 600 else gap_fn(2.08, 2.10),
                       raise_on_fail=False)
sel205 = g["rows_selected"]["2.05"]
check("(k) ゲート (入力 2.05): 2.05 は両方の解像度で成立、全部の境目からは 0.03、選び直した N の値は 2.05 でなく全部の境目から ≥ 0.05",
      sel205["at_N"]["ok"] and sel205["at_2N"]["ok"] and near(g["distance_input"], 0.03, 1e-9) and not near(g["chosen_N"], 2.05)
      and MS.boundary_distance(g["chosen_N"], g["common_boundaries"]) >= 0.05 - 1e-9
      and sorted((x["kind"], x["source"]) for x in g["common_boundaries"])
      == [("lower", "2N"), ("lower", "N"), ("upper", "2N"), ("upper", "N")],
      f"chosen_N {g['chosen_N']!r} distance_input {g['distance_input']!r} 境目 "
      f"{[(x['kind'], round(x['feasible_end'], 6), x['source']) for x in g['common_boundaries']]}")
check("(k) ゲート: row_chosen_N は選び直した点 (chosen_N) の行、row_input_N は入力 2.05 の行",
      near(g["row_chosen_N"]["M_K"], g["chosen_N"]) and near(g["row_input_N"]["M_K"], 2.05) and near(g["row_input_2N"]["M_K"], 2.05),
      f"{g['row_chosen_N']['M_K']!r} / {g['row_input_N']['M_K']!r}")

# ---------------------------------------------------------------------------------------------------------------------
print("== V1b M6 の目的の回帰 (sweep_axislaw_A_MK_summary.csv、n_axis 600)")
recs = {}
with open(M6_CSV) as f:
    for rec in csv.DictReader(f):
        if rec["role"] == "sweep" and int(rec["n_axis"]) == 600 and float(rec["L_c"]) in (45.0, 50.0):
            recs.setdefault(float(rec["L_c"]), {})[round(float(rec["M_K"]), 6)] = rec


def csv_eval(table):
    def fn(M):
        rec = table.get(round(M, 6))
        if rec is None:
            return {"M_K": M, "ok": False, "reasons": ["CSV に記録なし (0.1 刻みの外・1.5〜4.0 の外)"]}
        if rec["hard_gate"] == "ERROR":
            return row(M, None, error=rec["error"])
        hard = rec["hard_gate"] == "True"
        r = {"M_K": M, "theta_max_deg": float(rec["theta_max_deg"]), "margin_deg": float(rec["min_mu_minus_theta_deg"]),
             "min_sin_angle": float(rec["topology_margin"]), "mpp_sign_changes": int(rec["Mpp_sign_changes"]),
             "mpp_quality_ok": rec["single_peak"] == "True",      # CSV の single_peak = axis_gates.mpp_quality_ok
             "n_flip_interior": int(rec["flip_interior"]),
             "n_crossings_sampled": 0 if hard else None, "pileup_pairs_interior": 0 if hard else None,
             "qa_violations": [], "L_c": float(rec["L_c"]), "x_K": float(rec["x_K"]),
             # MOC のゲート: CSV (2026-08-17) は corrector の導入 (2026-10-07) の前の旧 MOC (fixed2 相当)。require_moc_gate は
             # fixed2 を通す (pass = None) ので今の扱いのまま成立とする
             "moc_gate_ok": True, "moc_gate_reason": None}
        r["reasons"] = MS.feasibility_reasons(r) + ([] if hard else ["hard_gate False (CSV: 壁 QA・設計壁・単調・反転・交差・重なりのどれか)"])
        r["ok"] = not r["reasons"]
        return r
    return fn


v1b = {}
for L_c, want, want_margin in ((45.0, 2.6, 3.0), (50.0, 2.7, 3.1)):
    tab = recs.get(L_c, {})
    fn = csv_eval(tab)
    feas = [fn(m) for m in sorted(tab)]
    feas = [x for x in feas if x["ok"]]
    mx = max(feas, key=lambda x: x["margin_deg"])
    r_off = MS.select_mknot(fn, 1.4, 6.0, coarse=0.1, fine=0.1, tie_deg=0.0, dM_safe=0.0, workers=1)
    r_def = MS.select_mknot(fn, 1.4, 6.0, coarse=0.1, fine=0.1, workers=1)
    v1b[L_c] = (r_off["chosen"], r_def["chosen"])
    check(f"V1b L_c {L_c:g}: CSV の記録 26 点 (1.5〜4.0)", len(tab) == 26, f"{len(tab)} 点")
    check(f"V1b L_c {L_c:g}: 余裕最大は {want_margin} (取り違えの検出の前提)", near(mx["M_K"], want_margin, 1e-6), f"{mx['M_K']}")
    check(f"V1b L_c {L_c:g}: 規則を切った設定 (ΔM_safe 0・tie 0) で {want} を選ぶ", near(r_off["chosen"], want, 1e-6),
          f"chosen {r_off['chosen']!r} θ {r_off['chosen_row']['theta_max_deg']:.4f}° 余裕 {r_off['chosen_row']['margin_deg']:.3f}°")
    print(f"   L_c {L_c:g}: 既定の設定 (tie 1e-3°・ΔM_safe 0.05、格子 0.1) の選択 {r_def['chosen']!r} ({r_def['reason']}; "
          f"境目 {[(x['kind'], round(x['feasible_end'], 4)) for x in r_def['boundaries']]}) — 記録 (合否なし)")
print(f"   V1b の選択値 (規則を切った設定, 既定の設定): {v1b}")

print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
