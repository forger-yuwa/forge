#!/usr/bin/env python3
"""M_knot の選び方 `forge_design.evaluate.mknot_select` の試験 — plan tooling-nozzle-mknot-auto §6 V2・V1b (§5.1 #2・#2b・#2c・#3b)。

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
      ゲートの row_chosen_N・row_chosen_2N は選び直した点の行 (入力の点の行は row_input_N・row_input_2N)。この入力の窓 (半幅 0.1) では
      選び直した 1.95 が窓の端なので、§5.1 #3b の後は半幅 0.1 だけなら未確定・既定の列なら 0.2 に広げて内部で PASS
  (l) 窓の端の選択 (諮問 mknot-card-holes、§5.1 #3b): 諮問の合成の例 (N・2N 共通の境目 1.65 / 2.25、θ_max = 15 + 0.0005·((M_K − 2.05)/0.1)²)
      で、窓 1.95〜2.15 (半幅 0.1) では同点の規則の選択 1.95 が端 → 半幅 0.2 に広げて同じ 1.95 が内部に入り PASS、距離 0.30。
      半幅 0.1 だけなら例外。上限 0.3 まで広げても端 (さらに平らな谷) → 例外。θ_max の最小の位置の端 (同点の規則の前) だけでは
      広げない。広げた窓は N を ±(半幅 + ΔM_safe) まで評価して窓の外の境目を見る
  (m) カードの YAML の要約 (wt_card の `_mknot_record` → `_mknot_yaml_summary`、§4.4; 設計チェーンは回さない): 掃引の表は card_log だけ、
      YAML の `_card.M_knot_rule` の行 ≤ 4 kB、card_log の記録の sha256 で結び、`mknot_card_log_check` で照合できる
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
# 窓は半幅 0.1 だけ (広げない): (k) は「境目をまとめない」をこの窓で確かめる試験。選び直した 1.95 は窓の端なので、§5.1 #3b の後は
# 上限 (この呼び出しでは 0.1) で端 = 未確定 になる (広げたときの振る舞いは (l) で確かめる)
g = MS.resolution_gate(None, 2.05, 600, workers=(1, 1), eval_factory=lambda n: gap_fn(2.00, 2.02) if n == 600 else gap_fn(2.08, 2.10),
                       raise_on_fail=False, half_widths=(0.1,))
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
check("(k) 半幅 0.1 だけでは選び直した 1.95 が窓の端のまま上限 → 未確定 (以前はこの入力で PASS を出していた; plan §9 の穴 (1))",
      g["pass"] is False and near(g["chosen_N"], 1.95) and g["edge_values"] == {"N": 1.95, "2N": 1.95}
      and [k for k, v in g["checks"].items() if not v] == ["chosen_inside_N", "chosen_inside_2N"], f"{g['checks']} {g['edge_values']}")
gk = MS.resolution_gate(None, 2.05, 600, workers=(1, 1), eval_factory=lambda n: gap_fn(2.00, 2.02) if n == 600 else gap_fn(2.08, 2.10),
                        raise_on_fail=False)
check("(k) 既定の半幅の列では 0.2 に広げて選び直し、2.05 は全部の境目から 0.03 で選ばれないまま、選択は窓の内部で全部の境目から ≥ 0.05",
      gk["pass"] and gk["n_widen"] == 1 and gk["half_width"] == 0.2 and not near(gk["chosen_N"], 2.05)
      and 2.05 in gk["rule_N"]["excluded_by_distance"] and MS.boundary_distance(gk["chosen_N"], gk["common_boundaries"]) >= 0.05 - 1e-9
      and gk["checks"]["chosen_inside_N"] and gk["checks"]["chosen_inside_2N"],
      f"chosen_N {gk['chosen_N']!r} 半幅 {gk['half_width']} 広げた {gk['n_widen']} {gk['checks']}")

print("== V2 (l) 窓の端の選択 (諮問 mknot-card-holes): 選び直した値が窓の端なら半幅を 0.1 → 0.2 → 0.3 に広げて選び直す")


def holes_fn(lo, hi, c=2.05, a=0.0005):
    """諮問の合成の例: N・2N 共通で成立 [lo, hi]、θ_max = 15 + a·((M_K − c)/0.1)² ° (平らな谷)。"""
    def fn(M):
        if not (lo <= M <= hi):
            return infeasible(M)
        return row(M, 15.0 + a * ((M - c) / 0.1) ** 2)
    return fn


def gate_at(fn, chosen, **kw):
    """select_mknot の表 (境目を含む) を渡して、入力 chosen でゲートを回す (カードの最終の組のゲートと同じ形: 入力は探索の選択と違いうる)。"""
    return MS.resolution_gate(None, chosen, 600, select_result=sel(fn), workers=(1, 1), eval_factory=lambda n: fn, **kw)


f_l = holes_fn(1.65, 2.25)
g1 = gate_at(f_l, 2.05, raise_on_fail=False, half_widths=(0.1,))
check("(l) 諮問の例の A: 半幅 0.1 (窓 1.95〜2.15) では同点の規則の選択 1.95 が窓の端、θ_max の最小 2.05 は窓の中央",
      near(g1["chosen_N"], 1.95) and near(g1["chosen_2N"], 1.95) and g1["edge_values"] == {"N": 1.95, "2N": 1.95}
      and g1["min_offset_N"] == 0 and g1["min_offset_2N"] == 0 and g1["checks"]["edge_N"] and g1["checks"]["edge_2N"]
      and g1["window"] == [1.95, 2.15],
      f"chosen {g1['chosen_N']!r}/{g1['chosen_2N']!r} 端 {g1['edge_values']} min_offset {g1['min_offset_N']}/{g1['min_offset_2N']}")
g2 = gate_at(f_l, 2.05)
check("(l) 諮問の例の B: 既定の列で半幅 0.2 (窓 1.85〜2.25) に 1 回広げ、同じ 1.95 が内部に入り PASS、全部の境目から 0.30",
      g2["pass"] and near(g2["chosen_N"], 1.95) and near(g2["chosen_2N"], 1.95) and g2["n_widen"] == 1 and g2["half_width"] == 0.2
      and g2["window"] == [1.85, 2.25] and g2["edge_values"] == {} and near(g2["distance_to_conservative_boundary"], 0.30, 1e-9)
      and near(g2["distance_2N"], 0.30, 1e-9) and abs(g2["moved_N"] + 0.10) <= 1e-9,
      f"{g2['checks']} chosen {g2['chosen_N']!r}/{g2['chosen_2N']!r} 距離 {g2['distance_to_conservative_boundary']!r} 半幅 {g2['half_width']}")
check("(l) 記録: 広げた回数・最後の半幅・各半幅の端だった値と次の手 (θ_max の最小の位置の判定は各段で別に記録)",
      [s["half_width"] for s in g2["widen"]] == [0.1, 0.2] and g2["widen"][0]["edge_values"] == {"N": 1.95, "2N": 1.95}
      and "広げて" in g2["widen"][0]["action"] and g2["widen"][1]["edge_values"] == {} and "内部" in g2["widen"][1]["action"]
      and all(s["min_offset_N"] == 0 and s["min_offset_2N"] == 0 for s in g2["widen"]) and g2["settings"]["half_widths"] == [0.1, 0.2, 0.3],
      str(g2["widen"]))
try:
    gate_at(f_l, 2.05, half_widths=(0.1,))
    check("(l) 半幅 0.1 だけ (上限 0.1) なら端のまま → 例外", False)
except MS.MknotUndeterminedError as e:
    c = e.result["checks"]
    check("(l) 半幅 0.1 だけ (上限 0.1) なら端のまま → 例外 (落ちるのは選び直した値の位置の判定だけ)",
          [k for k, v in c.items() if not v] == ["chosen_inside_N", "chosen_inside_2N"] and "未確定" in str(e) and "端" in str(e), str(c))
f_top = holes_fn(1.5, 4.0, a=0.0001)          # 成立 [1.5, 4.0]、谷がさらに平ら (2.05 ± 0.316 が同点): 境目から最も遠い方へ窓の端が逃げ続ける
try:
    gate_at(f_top, 2.05)
    check("(l) 半幅 0.3 まで広げても端 → 例外", False)
except MS.MknotUndeterminedError as e:
    r_ = e.result
    check("(l) 半幅 0.3 まで広げても選び直した値が端 (2.15 → 2.25 → 2.35) → 例外 (選択は未確定)、θ_max の最小の位置は中央で成立",
          r_["n_widen"] == 2 and r_["half_width"] == 0.3 and r_["edge_values"] == {"N": 2.35, "2N": 2.35}
          and [s["edge_values"].get("N") for s in r_["widen"]] == [2.15, 2.25, 2.35]
          and [k for k, v in r_["checks"].items() if not v] == ["chosen_inside_N", "chosen_inside_2N"]
          and "上限 0.3" in str(e) and "上限の半幅 0.3 でも窓の端" in r_["widen"][-1]["action"],
          f"{r_['checks']} 端 {[s['edge_values'] for s in r_['widen']]}")
_, gd = gate(quad(2.5), drift, raise_on_fail=False)
check("(l) θ_max の最小が窓の端 (同点の規則の前) だけなら窓を広げない (別の判定; (g) の drift の入力)",
      gd["n_widen"] == 0 and gd["checks"]["edge_2N"] is False and gd["checks"]["chosen_inside_N"] and gd["checks"]["chosen_inside_2N"]
      and gd["min_offset_2N"] == 10, f"広げた {gd['n_widen']} {gd['checks']}")


def beyond_fn(M):
    """成立 [1.5, 2.27]、θ_max = 15 + 10 (M_K − 2.24)² (急な谷、最小 2.24 は上の境目から 0.03)。"""
    if not (1.5 <= M <= 2.27):
        return infeasible(M)
    return row(M, 15.0 + 10.0 * (M - 2.24) ** 2)


ctx_tab = [dict(beyond_fn(round(2.05 + 0.01 * j, 12)), role="context") for j in range(-15, 16)]     # カードの最終の組の文脈 ±0.15
gb = MS.resolution_gate(None, 2.05, 600, select_result={"settings": {}, "table": ctx_tab}, workers=(1, 1), eval_factory=lambda n: beyond_fn,
                        raise_on_fail=False)
upper = [b for b in gb["common_boundaries"] if b["kind"] == "upper"]
check("(l) 広げた窓は N を ±(半幅 + ΔM_safe) まで評価し、窓の外 0.02 の上の境目 2.27 を見て、最小 2.24 (境目から 0.03) を選ばない",
      gb["pass"] and gb["n_widen"] == 1 and near(gb["chosen_N"], 2.22) and near(gb["chosen_2N"], 2.22)
      and any(abs(b["feasible_end"] - 2.27) <= 1e-3 and b["source"] == "N" for b in upper)
      and gb["distance_to_conservative_boundary"] >= 0.05 - 1e-9
      and any(r.get("role") == "context" and r["M_K"] > 2.25 + 1e-9 for r in gb["rows_N"]),
      f"chosen {gb['chosen_N']!r}/{gb['chosen_2N']!r} 上の境目 {[(round(b['feasible_end'], 6), b['source']) for b in upper]} {gb['checks']}")
try:
    MS.resolution_gate(None, 2.05, 600, workers=(1, 1), eval_factory=lambda n: f_l, half_widths=(0.2, 0.1))
    check("(l) half_widths が増えない列 → ValueError", False)
except ValueError as e:
    check("(l) half_widths が増えない列 → ValueError", "half_widths" in str(e), str(e)[:120])

print("== (m) カードの YAML の要約 (wt_card、§4.4): 掃引の表は card_log だけ、YAML の _card.M_knot_rule の行 ≤ 4 kB、記録の sha256 で結ぶ")
import json  # noqa: E402
import tempfile  # noqa: E402

from forge_design.intake import wt_card as W  # noqa: E402

YAML_LINE_MAX = 4096                     # YAML の _card.M_knot_rule の 1 行の上限 [byte] (plan §5.1 #3b の例)
s_m = sel(f_l)                                                           # 探索の選択 1.95 (諮問の例)
# 最初のゲートの入力は 2.00 (合成; 探索の 1.95 を中心にすると窓 1.85〜2.05 の θ_max の最小 2.05 が端で別の判定が落ちる)
g1_m = MS.resolution_gate(None, 2.00, 600, select_result=s_m, workers=(1, 1), eval_factory=lambda n: f_l)
gF_m = MS.resolution_gate(None, 2.05, 600, select_result=s_m, workers=(1, 1), eval_factory=lambda n: f_l)   # 最終の組: 1 回広げて PASS
cN_m = gF_m["chosen_N"]
ev_m = [{"stage": "最初の r_t", "attempt": 0, "r_t_m": 0.025, "M_K_in": 2.00, "chosen_N": g1_m["chosen_N"],
         "gate_pass": True, "gate_n_widen": g1_m["n_widen"], "action": "select_mknot の選択 → ゲートで選び直した chosen_N を固定"},
        {"stage": "最終の組 (最終の組で詰め直し (N は ±(半幅 + ΔM_safe)))", "attempt": 0, "r_t_m": 0.0249, "M_K_in": 2.05,
         "M_K_used_in_fixed_point": 2.05, "chosen_N": cN_m, "diff_vs_used": cN_m - 2.05, "gate_pass": True,
         "gate_n_widen": gF_m["n_widen"], "action": "差 > 0.01 → chosen_N で r_t の固定点をやり直す (1 回目、始まりは前の r_t)"}]
mk_m = {"gate_final": gF_m, "gate_first": g1_m, "sel": s_m, "M_K": cN_m, "sel_rt": 0.025, "N": 600, "events": ev_m, "redo": 1,
        "workers": {"N": 1, "2N": 1, "available_gb": None}, "t": 1.0}
card_m = {"recipe": "axismach_final_v1", "M_exit": 5.0, "M_knot_auto_source": "カードに geometry.M_knot を書いていない (既定は自動)"}
full_m = W._mknot_record(card_m, {"geometry": {"M_knot": cN_m}}, {"qa": {"theta_max_deg": gF_m["row_chosen_N"]["theta_max_deg"]}},
                         0.0249, mk_m)
CL = "problem_m.yaml.card_log.json"
summ_m = W._mknot_yaml_summary(full_m, CL)
text_m = W.render_problem_yaml({"name": "m", "_card": {"M_knot_rule": summ_m}}, ["要約の試験"], {})
line_m = [ln for ln in text_m.splitlines() if ln.startswith("  M_knot_rule:")]
line_full = W.render_problem_yaml({"name": "m", "_card": {"M_knot_rule": full_m}}, [], {}).splitlines()
nfull = max(len(ln.encode()) for ln in line_full if ln.startswith("  M_knot_rule:"))
check(f"(m) YAML の _card.M_knot_rule の行 ≤ {YAML_LINE_MAX} byte (全部の記録なら {nfull} byte)、掃引の表・行の詳細を持たない",
      len(line_m) == 1 and len(line_m[0].encode()) <= YAML_LINE_MAX < nfull and "table" not in summ_m
      and "rows_selected" not in json.dumps(summ_m) and "constraints_at_output" not in summ_m,
      f"{len(line_m[0].encode()) if line_m else None} byte / 全部 {nfull} byte")
check("(m) 要約: 規則の名前・選んだ値・距離・全部の境目・ゲートの判定 (広げた回数・端だった値)・選び直しの要約・card_log の記録の sha256",
      summ_m["rule"] == "min_theta_max_v1" and summ_m["chosen"] == cN_m and summ_m["distance_to_boundary"] == full_m["distance_to_boundary"]
      and [b["feasible_end"] for b in summ_m["boundaries"]] == [b["feasible_end"] for b in full_m["boundaries"]]
      and summ_m["gate"]["pass"] and summ_m["gate"]["n_widen"] == 1 and summ_m["gate"]["half_width"] == 0.2
      and summ_m["gate"]["widen"][0]["edge_values"] == {"N": 1.95, "2N": 1.95} and summ_m["ok_at_output"] == {"N": True, "2N": True}
      and [e["gate_n_widen"] for e in summ_m["events"]] == [g1_m["n_widen"], 1] and summ_m["redo"] == 1
      and summ_m["card_log"]["file"] == CL and summ_m["card_log"]["key"] == "M_knot_rule"
      and summ_m["card_log"]["sha256"] == W._mknot_record_sha256(json.loads(json.dumps(full_m, ensure_ascii=False))),
      str({k: summ_m.get(k) for k in ("chosen", "distance_to_boundary", "redo")}))
check("(m) geometry.M_knot のコメントの文は要約からでも全部の記録からでも同じ (窓を広げた記録を含む)",
      W._mknot_rule_text(summ_m) == W._mknot_rule_text(full_m) and "窓を 1 回広げた" in W._mknot_rule_text(summ_m),
      W._mknot_rule_text(summ_m)[:200])
with tempfile.TemporaryDirectory() as td:
    yp = Path(td) / "problem_m.yaml"
    yp.write_text(text_m)
    (Path(td) / CL).write_text(json.dumps(W._json_safe({"out_sha256": "x", "M_knot_rule": full_m}), indent=1, ensure_ascii=False))
    ok_chk = W.mknot_card_log_check(yp)
    bad = json.loads((Path(td) / CL).read_text())
    bad["M_knot_rule"]["table"]["select"][0]["theta_max_deg"] += 1e-12     # card_log の表を 1 箇所だけ変える
    (Path(td) / CL).write_text(json.dumps(bad, indent=1, ensure_ascii=False))
    bad_chk = W.mknot_card_log_check(yp)
    yp.write_text(text_m.replace(f"chosen: {cN_m!r}", f"chosen: {cN_m + 0.01!r}", 1))     # YAML の要約の選んだ値だけを変える
    bad_y = W.mknot_card_log_check(yp)
check("(m) mknot_card_log_check: 書いた YAML の要約と card_log の記録が sha256・作り直した要約で一致。card_log の表を 1 箇所変えると "
      "sha256 だけが不一致、YAML の要約の値を変えると要約が不一致",
      ok_chk["ok"] is True and ok_chk["summary_matches"] and bad_chk["ok"] is False and bad_chk["summary_matches"]
      and bad_chk["sha256_card_log"] != bad_chk["sha256_yaml"] and bad_y["ok"] is False and bad_y["summary_matches"] is False,
      f"{ok_chk['ok']} / 表 {bad_chk['ok']} {bad_chk['summary_matches']} / 要約 {bad_y['ok']} {bad_y['summary_matches']}")

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
