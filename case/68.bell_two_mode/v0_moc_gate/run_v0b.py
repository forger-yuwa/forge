#!/usr/bin/env python3
"""§6 V0-b: 前進 MOC ゲートの A/B の回し直し (CFD なし、事前登録の分類) — plan tooling-nozzle-bell-two-mode §6 V0・§5.1 #3c。

V0-a (`run_v0.py`、出力 `v0a_*`) は (D) の成立点が 8/135 に退化して判定不能だった (diagnostician 2026-10-11)。V0-b は
(D) の $x_m$ を設計変数から外して内部で決め (`wall_bell_c2.design_bell_c2`: 0.5〜4.0 r_t を 141 点で走査 → 成立区間で
$J=\\int(r''')^2$ 最小を黄金分割)、前進 MOC ゲートのガード 3 点 (§5.1 #3c) を入れた状態で回す。

対象:
- (D) の格子: θ_max {20, 25, 32.5, 40}° × θ_e {5, 10, 15}° × L_k (等面積比 15° 円錐の長さ (√ε − 1)/tan15° = 22.6575 r_t の
  60・70・80・90・100 %、ε 50) = 60 点。成立した点 (`design_bell_c2` が OK) だけゲートにかける。
- TOP (Phase 2 の最適形 θ_n 25.8°・θ_e 10°・R_u = R_d = R_t = 1.0、`geometry.wall.NozzleWall` の `bell_type: top`) を各 L_k で 5 点。
ゲート: `geometry.moc_forward.forward_moc` (①の interior_vec、analytic・converge)、ガス = 基準ケースの凍結 TP
(`RocketExhaust.moc_gas`、../exhaust_lox_lh2_pc100_of6.yaml)、初期線 = Hall のスロート特性線 (R = R_t = 1.0、γ = スロートの
局所 γ*)。主設定 n 321 (`hall_n321`)、確認 n 161 (`hall_n161`、合否の一致を記録)。ガス・ε・寸法は V0-a と同じ。

**事前登録の判定** (plan §6 V0、この順で出す。解釈はしない):
1. **有効条件**: 各 L_k に (D) の成立点 ≥ 3 かつ合計 ≥ 15。満たさなければ「判定不能」。
2. 有効なら通過率 ≥ 80 % → A、≤ 20 % → B、その間 → 中間。通過率の分母は成立点のうちゲートが判定できた点
   (verdict が pass か reject)。ゲートの判定不能 (初期線不整合・軸近傍・初期線の検査) は分母から外して別に数え、
   成立点を分母にした値も併記する。
3. C⁻ 間隔比の最小 < 0.2 の通過点は marginal と注記 (分類には入れる)。
4. **健全性の確認** (既知のパターン): (i) θ_max 40° は全 L で不受理、(ii) L60 は 32.5° でも不受理、(iii) 32.5°・L ≥ 70 % は
   通過。成立点についてパターンとの一致・不一致を記録し、崩れたら「ゲートが帯内の x_m に敏感」と記録する。

**走査範囲の確認** (1 回、§5.1 #3c): 60 点それぞれで x_m を (0, L) の 0.01 r_t 刻みで振って成立条件を見て、成立する x_m が
0.5〜4.0 r_t の外にあるかを記録する (`v0b_xm_scan_coverage.csv`)。

出力 (このディレクトリ): `v0b_d_design.csv` (60 点の x_m の内部決定)、`v0b_xm_scan_coverage.csv`、`v0b_gate.csv`
(ゲートの結果、設定ごとの縦持ち)、`v0b_summary.json` (判定・内訳・TOP・n 161 との一致・来歴)、`v0b_run.log`
(標準出力の写し)、図 `v0b_gate_map.png` (git 管理外)。
実行: design/.venv-opt の Python で `python case/68.bell_two_mode/v0_moc_gate/run_v0b.py` (数分)。
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry.moc_forward import (CallableWall, forward_moc, gate_record,  # noqa: E402
                                               hall_initial_line)
from forge_design.geometry.wall import NozzleWall, WallParams  # noqa: E402
from forge_design.geometry.wall_bell_c2 import (XM_SCAN_N, XM_SCAN_RANGE, BellC2Spec, BellC2Wall,  # noqa: E402
                                                cone15_length, design_bell_c2)

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
EPS = 50.0
R_T = 1.0
L_FRACS = (0.6, 0.7, 0.8, 0.9, 1.0)
THETA_MAX = (20.0, 25.0, 32.5, 40.0)
THETA_E = (5.0, 10.0, 15.0)
TOP = {"theta_n_deg": 25.8, "theta_e_deg": 10.0}
PRIMARY = "hall_n321"
CONFIGS = {"hall_n321": 321, "hall_n161": 161}
CLASS_A, CLASS_B = 0.80, 0.20
VALID_MIN_PER_L, VALID_MIN_TOTAL = 3, 15
MARGINAL_SEP = 0.2
COVER_DX = 0.01                      # 走査範囲の確認の刻み [r_t]
SOURCES = ("design/forge_design/geometry/wall_bell_c2.py", "design/forge_design/geometry/moc_forward.py",
           "design/forge_design/geometry/moc_kernel.py")


class Tee:
    def __init__(self, path):
        self.fh = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.fh.write(s + "\n")
        self.fh.flush()


def git_info() -> dict:
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                           capture_output=True, text=True).stdout
    sha = {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest()[:16] for s in SOURCES}
    return {"commit": head, "design_tree_dirty": [ln for ln in dirty.splitlines() if ln.strip()], "source_sha256_16": sha}


def pid_of(fr, tm, te):
    return f"D_L{int(round(fr * 100))}_tm{tm:g}_te{te:g}"


def classify(frac):
    return "A" if frac >= CLASS_A else ("B" if frac <= CLASS_B else "中間")


def coverage_scan(L, te, tm, spec) -> dict:
    """x_m を (0, L) の COVER_DX 刻みで振り、成立する x_m の範囲と、走査範囲 0.5〜4.0 r_t の外の成立点を数える。"""
    xs = np.arange(COVER_DX, L - 0.5 * COVER_DX, COVER_DX)
    oks = np.zeros(len(xs), dtype=bool)
    for i, x in enumerate(xs):
        try:
            oks[i] = BellC2Wall(L, te, tm, spec=spec, x_m=float(x)).feasibility()["ok"]
        except ValueError:
            oks[i] = False
    lo, hi = XM_SCAN_RANGE
    out_lo = xs[oks & (xs < lo)]
    out_hi = xs[oks & (xs > hi)]
    return {"n_cover": int(len(xs)), "n_cover_ok": int(np.sum(oks)),
            "xm_ok_min": float(xs[oks].min()) if np.any(oks) else None,
            "xm_ok_max": float(xs[oks].max()) if np.any(oks) else None,
            "n_ok_below_scan": int(len(out_lo)), "n_ok_above_scan": int(len(out_hi)),
            "ok_outside_scan": ";".join(f"{v:.2f}" for v in np.concatenate([out_lo, out_hi])[:20])}


def main() -> int:
    t_all = time.time()
    log = Tee(HERE / "v0b_run.log")
    gas = RocketExhaust.from_file(EXHAUST).moc_gas()
    spec = BellC2Spec(eps=EPS, R_t=R_T)
    Lc = cone15_length(EPS)
    # --- (D) の x_m の内部決定 (60 点) と走査範囲の確認 ---
    design_rows, cover_rows, d_walls = [], [], []
    for fr in L_FRACS:
        L = fr * Lc
        for tm in THETA_MAX:
            for te in THETA_E:
                pid = pid_of(fr, tm, te)
                d = design_bell_c2(L, te, tm, spec)
                rec = d.record()
                design_rows.append({"point": pid, "L_frac": fr, "L_rt": L, "theta_max_deg": tm, "theta_e_deg": te,
                                    **rec, "scan_violations": "; ".join(f"{k} {v}" for k, v in d.scan_violations.items()
                                                                        if not k.startswith("_"))})
                cv = coverage_scan(L, te, tm, spec)
                cover_rows.append({"point": pid, "L_frac": fr, "theta_max_deg": tm, "theta_e_deg": te,
                                   "design_status": d.status, "scan_intervals": rec["intervals"], **cv})
                if d.ok:
                    d_walls.append((pid, "D", fr, L, tm, te, d.x_m, d.wall))
                log(f"design {pid:20s} {d.status:10s} x_m {d.x_m if d.x_m is None else round(d.x_m, 5)}  "
                    f"区間 {rec['intervals'] or '-'}  端 {d.edge_touch}  走査 (0,L) の成立 x_m "
                    f"[{cv['xm_ok_min']}, {cv['xm_ok_max']}]・範囲外 {cv['n_ok_below_scan'] + cv['n_ok_above_scan']}")
    top_walls = []
    for fr in L_FRACS:
        L = fr * Lc
        nw = NozzleWall(WallParams(Ru=R_T, Rd=R_T, theta_a=np.deg2rad(TOP["theta_n_deg"]),
                                   theta_e=np.deg2rad(TOP["theta_e_deg"]), L_div=L, r_exit=np.sqrt(EPS), bell_type="top"))
        top_walls.append((f"TOP_L{int(round(fr * 100))}", "TOP", fr, L, TOP["theta_n_deg"], TOP["theta_e_deg"], None,
                          CallableWall(nw.r, nw.drdx)))
    # --- ゲート ---
    lines, gate_rows = {}, []
    for cfg, n in CONFIGS.items():
        lines[cfg] = hall_initial_line(R_T, gas, n)
        for (pid, fam, fr, L, a1, a2, xm, w) in d_walls + top_walls:
            res = forward_moc(w, L, gas, lines[cfg])
            rec = gate_record(res)
            row = {"config": cfg, "point": pid, "family": fam, "L_frac": fr, "L_rt": L,
                   "theta_max_or_n_deg": a1, "theta_e_deg": a2, "x_m": xm, **rec}
            gate_rows.append(row)
            log(f"{cfg:10s} {pid:20s} {rec['verdict']:11s} {rec['undecidable_kind'] or ''} {rec['codes'] or '-'}  "
                f"C⁻比 {rec['cminus_sep_ratio_min']}  交差開始 x {rec['reversal_onset_x']} (j {rec['reversal_onset_j']})  "
                f"根拠の j [{rec['evidence_j_min']}, {rec['evidence_j_max']}]  ({rec['time_s']:.1f} s)")
    # --- 集計 ---
    feas_by_L = {f"{int(round(fr * 100))}%": sum(1 for r in design_rows if r["design_status"] == "OK" and r["L_frac"] == fr)
                 for fr in L_FRACS}
    n_feas = sum(feas_by_L.values())
    valid = all(v >= VALID_MIN_PER_L for v in feas_by_L.values()) and n_feas >= VALID_MIN_TOTAL
    summary = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §6 V0 (V0-b、§5.1 #3c)", **git_info(),
               "gas": {"exhaust": str(EXHAUST.relative_to(ROOT)), "Tt_K": gas.Tt, "gamma_star": gas.gamma_star},
               "throat": {"R_t": R_T, "initial_line": "Hall throat characteristic (hall_initial_line)",
                          "gamma_hall": float(gas.gamma_throat())},
               "L_cone15_rt": Lc, "L_k_rt": [fr * Lc for fr in L_FRACS],
               "grid": {"theta_max_deg": THETA_MAX, "theta_e_deg": THETA_E, "L_frac": L_FRACS, "n_points": len(design_rows)},
               "xm_rule": {"scan_range_rt": XM_SCAN_RANGE, "scan_n": XM_SCAN_N,
                           "select": "成立区間のうち J 最小の走査点を含む区間 → 黄金分割で J 最小"},
               "rule": {"valid": f"各 L_k に成立 ≥ {VALID_MIN_PER_L} かつ合計 ≥ {VALID_MIN_TOTAL}",
                        "A": f"通過率 ≥ {CLASS_A:.0%}", "B": f"≤ {CLASS_B:.0%}", "中間": "その間",
                        "pass_fraction_denominator": "成立点のうちゲートが判定できた点 (pass + reject)。判定不能は別に数える",
                        "marginal": f"通過点で C⁻ 間隔比の最小 < {MARGINAL_SEP}"},
               "primary_config": PRIMARY}
    inf = {}
    for r in design_rows:
        if r["design_status"] == "OK":
            continue
        for part in r["scan_violations"].split("; "):
            if part:
                k, v = part.rsplit(" ", 1)
                inf[k] = inf.get(k, 0) + 1
    summary["d_design"] = {"n_feasible": n_feas, "feasible_by_L": feas_by_L,
                           "n_infeasible": len(design_rows) - n_feas,
                           "infeasible_points_by_violation_kind (走査点に 1 回でも現れた種類)": inf,
                           "n_edge_touch": sum(1 for r in design_rows if r["edge_touch"]),
                           "n_multi_interval": sum(1 for r in design_rows if (r["n_intervals"] or 0) >= 2),
                           "n_at_scan_point": sum(1 for r in design_rows if r["golden_at_scan_point"]),
                           "n_golden_fallback": sum(1 for r in design_rows if r["golden_fallback"]),
                           "x_m_range": [min(r["x_m"] for r in design_rows if r["x_m"] is not None),
                                         max(r["x_m"] for r in design_rows if r["x_m"] is not None)]}
    summary["xm_scan_coverage"] = {
        "dx_rt": COVER_DX, "n_points": len(cover_rows),
        "n_points_with_feasible_outside_scan": sum(1 for r in cover_rows if r["n_ok_below_scan"] + r["n_ok_above_scan"] > 0),
        "points_with_feasible_outside_scan": [r["point"] for r in cover_rows if r["n_ok_below_scan"] + r["n_ok_above_scan"] > 0],
        "feasible_xm_min_over_all": min((r["xm_ok_min"] for r in cover_rows if r["xm_ok_min"] is not None), default=None),
        "feasible_xm_max_over_all": max((r["xm_ok_max"] for r in cover_rows if r["xm_ok_max"] is not None), default=None),
        "n_points_feasible_in_fine_scan": sum(1 for r in cover_rows if r["n_cover_ok"] > 0)}
    summary["validity"] = {"valid": valid, "n_feasible": n_feas, "feasible_by_L": feas_by_L}
    summary["configs"] = {}
    for cfg in CONFIGS:
        rows = [r for r in gate_rows if r["config"] == cfg]
        d = [r for r in rows if r["family"] == "D"]
        npass = sum(r["verdict"] == "pass" for r in d)
        nrej = sum(r["verdict"] == "reject" for r in d)
        nund = sum(r["verdict"] == "undecidable" for r in d)
        frac = npass / (npass + nrej) if (npass + nrej) else None
        frac_all = npass / len(d) if d else None
        if not valid:
            cls = "判定不能 (有効条件を満たさない)"
        elif frac is None:
            cls = "判定不能 (判定できた点が 0)"
        else:
            cls = classify(frac)
        reasons = {}
        for r in d:
            for c in (r["codes"].split(";") if r["codes"] else []):
                reasons[c] = reasons.get(c, 0) + 1
        by_L = {}
        for fr in L_FRACS:
            dl = [r for r in d if r["L_frac"] == fr]
            by_L[f"{int(round(fr * 100))}%"] = {"feasible": len(dl), "pass": sum(r["verdict"] == "pass" for r in dl),
                                                "reject": sum(r["verdict"] == "reject" for r in dl),
                                                "undecidable": sum(r["verdict"] == "undecidable" for r in dl)}
        marginal = [r["point"] for r in d if r["verdict"] == "pass" and r["cminus_sep_ratio_min"] is not None
                    and r["cminus_sep_ratio_min"] < MARGINAL_SEP]
        # 健全性の確認 (既知のパターン)
        def pat(name, sel, want):
            pts = [r for r in d if sel(r)]
            bad = [f"{r['point']} ({r['verdict']})" for r in pts if r["verdict"] != want]
            # match: True = 一致、False = 崩れた、None = パターンに当たる成立点が無く検査できない
            return {"pattern": name, "expected": want, "n_points": len(pts), "match": (not bad) if pts else None,
                    "n_mismatch": len(bad), "mismatch": bad, "points": [r["point"] for r in pts]}
        health = [pat("θ_max 40° は全 L で不受理", lambda r: r["theta_max_or_n_deg"] == 40.0, "reject"),
                  pat("L60 は 32.5° でも不受理", lambda r: r["L_frac"] == 0.6 and r["theta_max_or_n_deg"] == 32.5, "reject"),
                  pat("32.5°・L ≥ 70 % は通過", lambda r: r["L_frac"] >= 0.7 and r["theta_max_or_n_deg"] == 32.5, "pass")]
        health_ok = all(h["match"] is True for h in health)
        summary["configs"][cfg] = {
            "n_line": CONFIGS[cfg], "d_feasible": len(d), "d_pass": npass, "d_reject": nrej, "d_undecidable": nund,
            "d_undecidable_kinds": {k: sum(1 for r in d if r["undecidable_kind"] == k)
                                    for k in sorted({r["undecidable_kind"] for r in d if r["undecidable_kind"]})},
            "pass_fraction": frac, "pass_fraction_over_feasible": frac_all, "classification": cls,
            "by_L": by_L, "marginal_pass": marginal, "d_fail_codes": reasons,
            "health_check": health, "health_all_match": health_ok,
            "health_note": (None if health_ok else
                            "既知のパターンが崩れた — ゲートが帯内の x_m に敏感" if any(h["match"] is False for h in health)
                            else "検査できる成立点が無いパターンがある"),
            "d_points": {r["point"]: {"verdict": r["verdict"], "undecidable_kind": r["undecidable_kind"], "codes": r["codes"],
                                      "x_m": r["x_m"], "cminus_sep_ratio_min": r["cminus_sep_ratio_min"],
                                      "cminus_sep_j": r["cminus_sep_j"], "reversal_onset_x": r["reversal_onset_x"],
                                      "reversal_onset_r": r["reversal_onset_r"], "reversal_onset_j": r["reversal_onset_j"],
                                      "evidence_j": [r["evidence_j_min"], r["evidence_j_max"]]} for r in d},
            "top": {r["point"]: {"verdict": r["verdict"], "undecidable_kind": r["undecidable_kind"], "codes": r["codes"],
                                 "cminus_sep_ratio_min": r["cminus_sep_ratio_min"],
                                 "reversal_onset_x": r["reversal_onset_x"], "reversal_onset_r": r["reversal_onset_r"],
                                 "reversal_onset_j": r["reversal_onset_j"],
                                 "evidence_j": [r["evidence_j_min"], r["evidence_j_max"]]}
                    for r in rows if r["family"] == "TOP"}}
        if cfg != PRIMARY:
            prim = {r["point"]: r["verdict"] for r in gate_rows if r["config"] == PRIMARY}
            diff = sorted(f"{r['point']} ({prim.get(r['point'])} → {r['verdict']})" for r in rows
                          if r["verdict"] != prim.get(r["point"]))
            summary["configs"][cfg]["verdict_differs_from_primary"] = diff
            summary["configs"][cfg]["agrees_with_primary"] = not diff
    summary["initial_lines"] = {cfg: {"n": len(ln["x"]), "x_axis_end": float(ln["x"][0]), "M_axis_end": float(ln["M"][0]),
                                      "M_wall": float(ln["M"][-1]), "thr_axis": float(ln["thr_axis"])}
                                for cfg, ln in lines.items()}
    summary["elapsed_s"] = time.time() - t_all
    # --- 書き出し ---
    for name, rows in (("v0b_d_design.csv", design_rows), ("v0b_xm_scan_coverage.csv", cover_rows),
                       ("v0b_gate.csv", gate_rows)):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
    (HERE / "v0b_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    try:
        plot_map(design_rows, gate_rows, Lc)
    except Exception as e:                      # 図は補助 (git 管理外)。失敗しても集計は残す
        log(f"図の作成に失敗: {e!r}")
    p = summary["configs"][PRIMARY]
    log(f"\nV0-b ({PRIMARY}): (D) 成立 {n_feas}/{len(design_rows)} {feas_by_L} → 有効条件 {'満たす' if valid else '満たさない'}・"
        f"通過 {p['d_pass']}・不受理 {p['d_reject']}・判定不能 {p['d_undecidable']}・通過率 {p['pass_fraction']}・"
        f"分類 {p['classification']}・marginal {p['marginal_pass']}・健全性 {'一致' if p['health_all_match'] else '不一致'}")
    for h in p["health_check"]:
        log(f"  健全性 {h['pattern']}: {({True: '一致', False: '不一致', None: '検査できない'})[h['match']]} "
            f"({h['n_points']} 点、不一致 {h['mismatch']})")
    log(f"  TOP: " + "  ".join(f"{k} {v['verdict']} (交差開始 x {v['reversal_onset_x']}, r {v['reversal_onset_r']})"
                               for k, v in p["top"].items()))
    c2 = summary["configs"]["hall_n161"]
    log(f"  n 161: 通過 {c2['d_pass']}・不受理 {c2['d_reject']}・判定不能 {c2['d_undecidable']} → {c2['classification']}・"
        f"主と違う点 {c2['verdict_differs_from_primary'] or 'なし'}")
    cv = summary["xm_scan_coverage"]
    log(f"  走査範囲の確認: 成立する x_m が 0.5〜4.0 r_t の外にある点 {cv['n_points_with_feasible_outside_scan']}/{cv['n_points']}・"
        f"成立する x_m の全体の範囲 [{cv['feasible_xm_min_over_all']}, {cv['feasible_xm_max_over_all']}]  "
        f"({summary['elapsed_s']:.0f} s)")
    return 0


def _font():
    import matplotlib
    from matplotlib import font_manager
    for f in (Path.home() / ".fonts").glob("NotoSansCJKjp-Regular.otf"):
        font_manager.fontManager.addfont(str(f))
        matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"


def plot_map(design_rows, gate_rows, Lc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _font()
    prim = {r["point"]: r for r in gate_rows if r["config"] == PRIMARY}
    des = {r["point"]: r for r in design_rows}
    fig, axs = plt.subplots(1, len(L_FRACS), figsize=(18, 4.6), squeeze=False)
    for i, fr in enumerate(L_FRACS):
        ax = axs[0, i]
        top = prim.get(f"TOP_L{int(round(fr * 100))}")
        for a, tm in enumerate(THETA_MAX):
            for b, te in enumerate(THETA_E):
                pid = pid_of(fr, tm, te)
                if des[pid]["design_status"] != "OK":
                    col, txt = "0.85", "不成立"
                else:
                    g = prim[pid]
                    col = {"pass": "#7fc97f", "reject": "#f08080", "undecidable": "#fdc086"}[g["verdict"]]
                    txt = (f"{ {'pass': '通過', 'reject': '不受理', 'undecidable': '判定不能'}[g['verdict']] }\n"
                           f"x_m {g['x_m']:.3f}\nC-比 {g['cminus_sep_ratio_min']:.2f}")
                ax.add_patch(plt.Rectangle((b, a), 1, 1, color=col, ec="w"))
                ax.text(b + 0.5, a + 0.5, txt, ha="center", va="center", fontsize=7)
        ax.set_xlim(0, len(THETA_E))
        ax.set_ylim(0, len(THETA_MAX))
        ax.set_xticks([k + 0.5 for k in range(len(THETA_E))], [f"θe {t:g}°" for t in THETA_E], fontsize=8)
        ax.set_yticks([k + 0.5 for k in range(len(THETA_MAX))], [f"θmax {t:g}°" for t in THETA_MAX], fontsize=8)
        ttl = f"L = {fr:.0%} ({fr * Lc:.2f} r_t)"
        if top is not None:
            ttl += f"\nTOP: {top['verdict']} (C-比 {top['cminus_sep_ratio_min']:.2f})"
        ax.set_title(ttl, fontsize=9)
    fig.suptitle(f"V0-b 前進 MOC ゲート ({PRIMARY}、x_m は内部決定): 灰 = 不成立、緑 = 通過、赤 = 不受理、橙 = 判定不能",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(HERE / "v0b_gate_map.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
