"""case/70: M_knot の選び方 (forge_design.evaluate.mknot_select、規則 min_theta_max_v1) の V1a・V6 (設計計算だけ、forge 0 step)。

plan: plans/active/tooling-nozzle-mknot-auto.md §6 V1a・V6 (§5.1 #2・#2b・#2c)。問題は problem_m5_air_ns2_euler.yaml (r_t 0.0250733 m を固定、
長さ 1.2 m を from_length で固定)。出力は _mknot_select/ (JSON とログ)。

  v6   低い側の境目と最小の A/B: n_axis_inv 2400 と 4800 で M_K 1.90〜2.20 を 0.01 刻みに評価し、成立の切り替わりの区間を
       幅 ≤ 0.001 まで二分法で詰める。2400 と 4800 の境目を全部残し (対にしてまとめない; mknot_select.conservative_boundaries、
       §5.1 #2c)、全部の境目の成立側の端からの距離で両方の解像度の選択をそれぞれ選び直す (§4.2 の規則: 制約と距離 → 最小 →
       同点; 0.01 の格子なので詰めは格子そのもの)。各解像度の境目だけで選んだ値 (#2 の V6 と同じ選び方) も記録する (合否には
       使わない)。
       合格 (§6 V6): 選択値の差 ≤ 0.01・両方の選択点が両方の解像度で全部の制約を満たす・選択点から全部の境目まで ≥ 0.05。
       距離は評価表 (全部の評価済みの点) から求め直した境目で測り、選び直しの境目の探し直しが上限 (3 回) で収まらなければ
       不合格 (settled; §5.1 #3d、result 段レビュー M1)。
       不合格でも表と判定を JSON に残す (境目が動いたのか最小が動いたのかを表で区別)。
  v1a  固定の入力の再現: select_mknot を全範囲で回す (N 2400)。合格 (§6 V1a): (a) 選んだ点が 2.03〜2.10 で保守的な境目から ≥ 0.05、
       (b) 評価の値が _mknot_sweep/low_end_N2400.json・low_end2_N2400.json・res_N2400.json と同じ M_K で相対 1e-9。
       §4.2 の格子の原点は M_A + 0.1 (M_A は軸のアンカーの M、M5 では 1.2293) なので select_mknot の評価点は掃引の丸い M_K
       (1.98・2.05 …) と一致しない。そのため (b) は同じ評価の関数 (evaluate_mknot) を掃引の M_K で呼んで照合する (select_mknot の
       評価点と掃引の M_K の一致数も記録する)。--gate で解像度のゲート (§4.2、N 2400/2N 4800) も回して記録する (V1a の合否には入れない)。

usage (case dir か forge-deliv の直下、design/.venv-opt):
  python case/70.m5_air_wt_demo/mknot_select_ab.py v6  [--workers 4 --workers-2n 2]
  python case/70.m5_air_wt_demo/mknot_select_ab.py v1a [--gate] [--workers 4 --workers-2n 2]
一時の YAML は $TMPDIR に書く (case の中に置かない)。
"""
import json
import sys
import time
from pathlib import Path

import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))

from forge_design.evaluate import mknot_select as MS  # noqa: E402

PROB = C / "problem_m5_air_ns2_euler.yaml"
OUT = C / "_mknot_select"
SWEEP = C / "_mknot_sweep"
CMP_KEYS = ("theta_max_deg", "margin_deg", "x_at_min_margin", "min_sin_angle", "n_flip_interior", "n_crossings_sampled",
            "pileup_pairs_interior", "n_axis", "L_c", "x_K", "x_F", "max_kappa_down", "J_dkappa_down")


def _arg(argv, name, default):
    if name in argv:
        k = argv.index(name)
        v = argv[k + 1]
        del argv[k:k + 2]
        return type(default)(v)
    return default


class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        self.f.write(line + "\n")
        self.f.flush()


def _fmt(v, f):
    return "—" if v is None else format(v, f)


def _slim(r):
    keys = ("M_K", "ok", "reasons", "role", "theta_max_deg", "margin_deg", "min_sin_angle", "n_flip_interior",
            "n_crossings_sampled", "pileup_pairs_interior", "mpp_sign_changes", "mpp_quality_ok", "qa_violations",
            "law_violations", "L_c", "x_K", "x_F", "n_axis", "M_A", "M_d_moc", "moc_corrector", "moc_gate_pass", "moc_gate_ok",
            "moc_gate_reason", "error", "time_s")
    return {k: r.get(k) for k in keys if k in r}


# --- V6 ---------------------------------------------------------------------------------------------------------------
def run_v6(raw, workers, workers_2n, log):
    grid = [round(1.90 + 0.01 * i, 2) for i in range(31)]
    rng = (grid[0], grid[-1])
    res, evs = {}, {}
    for N, w in ((2400, workers), (4800, workers_2n)):
        log(f"V6: n_axis_inv {N} で M_K 1.90〜2.20 (31 点、並列 {w})")
        ev = MS.Evaluator(MS.make_eval_fn(raw, N, base_dir=C), w, log)
        with ev:
            ev.many(grid, "grid")
            bnd = MS.find_boundaries(ev, 1e-3)
        log(f"  境目: {[(b['kind'], round(b['lo'], 6), round(b['hi'], 6)) for b in bnd]}")
        # 記録: その解像度の境目だけで選んだ値 (#2 の V6 と同じ選び方の比較用; 合否には使わない)
        try:
            own = MS.apply_rule([r for r in ev.table() if r.get("role") == "grid"], bnd, tie_deg=1e-3, dM_safe=0.05,
                                safe_eval=lambda Ms, ev=ev: ev.many(Ms, "safe_own"), safe_range=rng)
            own = {k: v for k, v in own.items() if k != "chosen_row"} | {"chosen_row": _slim(own["chosen_row"])}
        except MS.MknotSelectionError as e:
            own = {"chosen": None, "error": str(e).splitlines()[0]}
        log(f"  その解像度の境目だけの選択 (記録): {own.get('chosen')!r} ({own.get('reason') or own.get('error')})")
        res[N] = {"boundaries": bnd, "rule_own_boundary": own}
        evs[N] = ev
    # 全部の境目 (2400 と 4800 の両方、まとめない) で両方を選び直す (§4.2; 候補は 0.01 の格子と境目 + ΔM_safe の点)
    grid_rows = {N: [r for r in evs[N].table() if r.get("role") == "grid"] for N in (2400, 4800)}
    # 境目の探し直しが上限で収まらなければ未確定 (settled False → 不合格; §5.1 #3d、result 段レビュー M1)
    rc = MS.reselect_common(evs[2400], evs[4800], grid_rows[2400], grid_rows[4800], tie_deg=1e-3, dM_safe=0.05,
                            bisect_tol=1e-3, safe_range=rng, candidate_roles=("grid", "safe"), raise_unsettled=False, log=log)
    for N, tag in ((2400, "N"), (4800, "2N")):
        rule = rc[tag]
        if rule.get("chosen") is not None:
            rule = {k: v for k, v in rule.items() if k != "chosen_row"} | {"chosen_row": _slim(rule["chosen_row"])}
        res[N]["rule"] = rule
        res[N]["boundaries_after_reselect"] = rc["boundaries_N" if tag == "N" else "boundaries_2N"]
        res[N]["n_eval"] = evs[N].n_eval
        log(f"  {N}: 全部の境目で選び直した選択 {rule.get('chosen')!r} ({rule.get('reason') or rule.get('error')})")
    c24, c48 = res[2400]["rule"].get("chosen"), res[4800]["rule"].get("chosen")
    # 交差の評価: 各解像度の選択点を相手の解像度でも評価する (格子の外の点なら新しく評価; 少数なので逐次)
    pts = [c for c in (c24, c48) if c is not None]
    cross = {}
    for N in (2400, 4800):
        at = {MS._key(r["M_K"]): r for r in evs[N].many(pts, "cross")}
        for c, src in ((c24, 2400), (c48, 4800)):
            if c is not None:
                cross[f"chosen_{src}_at_{N}"] = _slim(at[MS._key(c)])
    # 距離は評価表 (全部の評価済みの点; 交差の評価を含む) から求め直した境目で判定する (§5.1 #3d、レビュー M1)。
    # 選び直しに使った境目 (rc["common"]) と食い違えば settled も不合格
    b_fin = {N: MS.find_boundaries(evs[N], 1e-3) for N in (2400, 4800)}
    common = MS.conservative_boundaries(b_fin[2400], b_fin[4800])
    settled = (rc["settled"] and MS._same_boundaries(b_fin[2400], rc["boundaries_N"])
               and MS._same_boundaries(b_fin[4800], rc["boundaries_2N"]))
    evs = {N: evs[N].table() for N in (2400, 4800)}
    rows_at = {N: {MS._key(r["M_K"]): r for r in evs[N]} for N in (2400, 4800)}

    def ok_at(c, N):
        r = rows_at[N].get(MS._key(c)) if c is not None else None
        return bool(r and MS._feasible(r))

    checks = {
        "settled": settled,
        "diff_le_0p01": c24 is not None and c48 is not None and abs(c24 - c48) <= 0.01 + 1e-9,
        "constraints_both": all(ok_at(c, N) for c in (c24, c48) for N in (2400, 4800)),
        "distance_ge_0p05": all(c is not None and MS.boundary_distance(c, common) >= 0.05 - 1e-9 for c in (c24, c48)),
    }
    judg = {"checks": checks, "pass": all(checks.values()), "chosen_2400": c24, "chosen_4800": c48,
            "diff": (None if c24 is None or c48 is None else c48 - c24),
            "distance_2400": (None if c24 is None else MS.boundary_distance(c24, common)),
            "distance_4800": (None if c48 is None else MS.boundary_distance(c48, common)),
            "conservative_boundary": [max((b["feasible_end"] for b in common if b["kind"] == "lower"), default=None),
                                      min((b["feasible_end"] for b in common if b["kind"] == "upper"), default=None)],
            "common_boundaries": common, "common_boundaries_used": rc["common"], "reselect_warnings": rc["warnings"],
            "chosen_own_boundary": {str(N): res[N]["rule_own_boundary"].get("chosen") for N in (2400, 4800)}}
    # 表 (0.01 の格子と二分法の点)
    allM = sorted(set(rows_at[2400]) | set(rows_at[4800]))
    table = []
    for k in allM:
        a, b = rows_at[2400].get(k), rows_at[4800].get(k)
        table.append({"M_K": k, "N2400": (None if a is None else _slim(a)), "N4800": (None if b is None else _slim(b))})
    log("V6 の表 (○ = 全部の制約が成立):")
    log("  M_K       | 2400: 成立 θ_max      余裕    sin    M″ | 4800: 成立 θ_max      余裕    sin    M″ | 理由")
    for t in table:
        cells = []
        why = []
        for N in ("N2400", "N4800"):
            r = t[N]
            if r is None:
                cells.append("     (未評価)                           ")
                continue
            cells.append(f"{'○' if MS._feasible(r) else '×'}  {_fmt(r.get('theta_max_deg'), '10.6f')}  {_fmt(r.get('margin_deg'), '6.3f')}  "
                         f"{_fmt(r.get('min_sin_angle'), '.4f')}  {r.get('mpp_sign_changes', '—')}")
            if r.get("reasons"):
                why.append(f"{N}: " + "; ".join(r["reasons"])[:90])
        log(f"  {t['M_K']:<9.6f} | {cells[0]} | {cells[1]} | {' / '.join(why)}")
    log(f"V6 判定: {'PASS' if judg['pass'] else 'FAIL'} {checks} 選択 (全部の境目で選び直し) 2400 {c24!r}・4800 {c48!r}、距離 "
        f"{_fmt(judg['distance_2400'], '.4f')}・{_fmt(judg['distance_4800'], '.4f')}、効いている境目 (下の最大・上の最小) {judg['conservative_boundary']}"
        f"、各解像度の境目だけの選択 (記録) {judg['chosen_own_boundary']}")
    return {"problem": str(PROB), "grid": grid, "rule": MS.RULE_NAME, "N": {str(k): v for k, v in res.items()},
            "cross": cross, "judgment": judg, "table": table}


# --- V1a --------------------------------------------------------------------------------------------------------------
def _rel_ok(a, b):
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, int) and isinstance(b, int):
        return a == b
    if a is None or b is None:
        return a is b
    a, b = float(a), float(b)
    return abs(a - b) <= 1e-9 * max(abs(a), abs(b), 1e-300)


def run_v1a(raw, workers, workers_2n, gate, log):
    log("V1a: M_A と今の M_knot の行 (一回評価)")
    r0 = MS.evaluate_mknot(raw, float(raw["geometry"]["M_knot"]), base_dir=C)
    M_A, M_d = r0["M_A"], r0["M_d_moc"]
    log(f"  M_A {M_A!r}・M_d {M_d!r}・今の M_knot {raw['geometry']['M_knot']} の θ_max {r0.get('theta_max_deg')!r}")
    sel = MS.select_mknot(MS.make_eval_fn(raw, 2400, base_dir=C), M_A, M_d, workers=workers, log=log)
    log(f"  選択 {sel['chosen']!r} ({sel['reason']})、境目 {[(b['kind'], round(b['lo'], 6), round(b['hi'], 6)) for b in sel['boundaries']]}"
        f"、距離 {sel['distance_to_boundary']:.4f}、警告 {sel['warnings']}、評価 {sel['n_eval']} 点 {sel['elapsed_s']:.0f} s")
    log("select_mknot の表:\n" + MS.format_table(sel["table"]))
    # (b) 掃引の記録との照合
    ref = []
    for name in ("low_end_N2400.json", "low_end2_N2400.json", "res_N2400.json"):
        for r in json.loads((SWEEP / name).read_text()):
            ref.append((name, r))
    selM = {MS._key(r["M_K"]) for r in sel["table"]}
    overlap = sorted({MS._key(r["M_K"]) for _, r in ref} & selM)
    log(f"(b) select_mknot の評価点と掃引の M_K の一致: {len(overlap)} 点 {overlap} — 掃引の M_K を evaluate_mknot で評価して照合")
    with MS.Evaluator(MS.make_eval_fn(raw, 2400, base_dir=C), workers, log) as ev:
        ev.seed(sel["table"])
        mine = {MS._key(r["M_K"]): r for r in ev.many([r["M_K"] for _, r in ref], "sweep_ref")}
    cmp_rows, n_bad = [], 0
    for name, r in ref:
        m = mine[MS._key(r["M_K"])]
        if not r.get("ok"):
            same = (not m.get("ok")) and m.get("error") == r.get("error")
            diffs = {} if same else {"error": [r.get("error"), m.get("error")]}
        else:
            diffs = {k: [r.get(k), m.get(k)] for k in CMP_KEYS if not _rel_ok(r.get(k), m.get(k))}
            if list(r.get("qa_violations") or []) != list(m.get("qa_violations") or []):
                diffs["qa_violations"] = [r.get("qa_violations"), m.get("qa_violations")]
            same = not diffs
        n_bad += (not same)
        cmp_rows.append({"file": name, "M_K": r["M_K"], "same": same, "diffs": diffs, "ok_ref": r.get("ok"),
                         "ok_select_rule": m.get("ok"), "reasons": m.get("reasons")})
    b_lower = [b for b in sel["boundaries"] if b["kind"] == "lower"]
    checks = {"a_range_2.03_2.10": 2.03 - 1e-9 <= sel["chosen"] <= 2.10 + 1e-9,
              "a_distance_ge_0p05": sel["distance_to_boundary"] >= 0.05 - 1e-9,
              "b_same_values_rel1e-9": n_bad == 0 and len(cmp_rows) > 0}
    log(f"V1a 判定: {'PASS' if all(checks.values()) else 'FAIL'} {checks}; 照合 {len(cmp_rows) - n_bad}/{len(cmp_rows)} 点が相対 1e-9 で一致")
    out = {"problem": str(PROB), "rule": MS.RULE_NAME, "M_A": M_A, "M_d": M_d, "row_current_M_knot": _slim(r0),
           "select": {k: v for k, v in sel.items() if k not in ("table", "chosen_row")} | {"chosen_row": _slim(sel["chosen_row"])},
           "table": [_slim(r) for r in sel["table"]], "lower_boundaries": b_lower,
           "compare_sweep": {"overlap_select_points": overlap, "rows": cmp_rows, "n_mismatch": n_bad},
           "checks": checks, "pass": all(checks.values())}
    if gate:
        log("解像度のゲート (§4.2; V1a の合否には入れない)")
        g = MS.resolution_gate(raw, sel["chosen"], 2400, select_result=sel, base_dir=C, workers=(workers, workers_2n),
                               raise_on_fail=False, log=log)
        log(f"  ゲート: {'PASS' if g['pass'] else 'FAIL'} {g['checks']} 全部の境目で選び直した N {g['chosen_N']!r} "
            f"({g['chosen_N_reason']})・2N {g['chosen_2N']!r} ({g['chosen_2N_reason']})、入力からの動き {g['moved_N']!r}、"
            f"全部の境目まで N {_fmt(g['distance_to_conservative_boundary'], '.4f')}・2N {_fmt(g['distance_2N'], '.4f')}、全部の境目 "
            f"{[(b['kind'], round(b['feasible_end'], 6), b['source'], b['nearest_other']) for b in g['common_boundaries'] if abs(b['feasible_end'] - sel['chosen']) < 0.5]}"
            f"、2N の境目 {[(b['kind'], round(b['lo'], 6), round(b['hi'], 6)) for b in g['boundaries_2N']]}、理由 {g['reasons']}、"
            f"警告 {g['warnings']}")
        out["resolution_gate"] = {k: v for k, v in g.items()
                                  if k not in ("rows_N", "rows_2N", "row_chosen_N", "row_chosen_2N", "row_input_N", "row_input_2N",
                                               "rows_selected")} | {
            "rows_N": [_slim(r) for r in g["rows_N"]], "rows_2N": [_slim(r) for r in g["rows_2N"]],
            **{k: (None if g[k] is None else _slim(g[k])) for k in ("row_chosen_N", "row_chosen_2N", "row_input_N", "row_input_2N")},
            "rows_selected": {k: {kk: _slim(vv) for kk, vv in v.items()} for k, v in g["rows_selected"].items()}}
    return out


if __name__ == "__main__":
    argv = list(sys.argv[1:])
    workers = _arg(argv, "--workers", 4)
    workers_2n = _arg(argv, "--workers-2n", 2)
    gate = "--gate" in argv
    if gate:
        argv.remove("--gate")
    what = argv[0] if argv else ""
    if what not in ("v6", "v1a"):
        sys.exit(__doc__)
    OUT.mkdir(exist_ok=True)
    raw = yaml.safe_load(PROB.read_text())
    log = Log(OUT / f"{what}.log")
    t0 = time.time()
    res = run_v6(raw, workers, workers_2n, log) if what == "v6" else run_v1a(raw, workers, workers_2n, gate, log)
    res["elapsed_s"] = time.time() - t0
    (OUT / f"{what}.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))
    log(f"書いた: {OUT / (what + '.json')} ({res['elapsed_s']:.0f} s)")
