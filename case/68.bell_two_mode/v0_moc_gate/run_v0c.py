#!/usr/bin/env python3
"""§6 V0-c: (D) の成立条件の穴を塞いだ後の V0 の記録の差し替え (CFD なし、事前登録の判定) — plan tooling-nozzle-bell-two-mode
§6 V0-c・§5.1 #3d。

V0-b (`run_v0b.py`、出力 `v0b_*`) は「中間」で解釈が確定した (diagnostician 2026-10-11、
notes/reviews/2026-10-11-bell-v0b-diagnostician.md)。その後 #3d で (D) の壁 (`wall_bell_c2`) を 2 点直した:
(1) 成立条件 (ii) を「θ' の符号変化がちょうど 1 回」から「r'' ≥ −tol on [0, x_m] かつ r'' ≤ tol on [x_m, L]」へ
(x_m で r'' が 0 に接して符号を変えない壁が通っていた)、(2) 成立区間の端を走査点と隣の不成立点の間で二分法で詰めてから
黄金分割 (x_m の 0.025 r_t 刻みの量子化を消す。J 最小の規則は不変)。V0-c は族の定義が変わったための記録の更新で、
対象・ゲート・閾値は V0-b と同じ:

- (D) の格子: θ_max {20, 25, 32.5, 40}° × θ_e {5, 10, 15}° × L_k (等面積比 15° 円錐の長さの 60・70・80・90・100 %、ε 50)
  = 60 点。成立した点だけゲートにかける。TOP (θ_n 25.8°・θ_e 10°・R_u = R_d = R_t = 1.0) を各 L_k で 5 点。
- ゲート `geometry.moc_forward.forward_moc` (analytic・converge)、ガス = 基準ケースの凍結 TP、初期線 = Hall (R = R_t = 1.0)。
  主設定 n 321、確認 n 161。

**事前登録の判定** (plan §6 V0-c、数値で出すところまで。解釈はしない):
- V0 の分類 (V0-b と同じ規則): 有効条件 (各 L_k に成立 ≥ 3 かつ合計 ≥ 15) → 通過率 (分母 = pass + reject) ≥ 80 % A、
  ≤ 20 % B、その間「中間」。
- **A**: 分類が中間のまま、かつ n 161 との不一致 ≤ 2 → V0 を閉じる。
- **B**: 分類が変わる・有効条件が破れる・V0-b の既存 29 点 (V0-b で成立した 31 点から #3d で不成立になる接触型 2 点を除いたもの)
  のうち合否が 4 点以上反転 → 再諮問。
- 合否 = ゲートの `gate_pass` (pass なら True、reject と判定不能は False)。参考として pass↔reject だけの反転と、判定不能を含む
  verdict の変化も数える。n 161 との不一致は V0-b と同じく全点 (D + TOP) の verdict の違いで数え、D だけの数も併記する。
  読み方 (合否の定義・不一致の範囲) で判定が変わる場合は `judgement.reading_sensitive` に記録する。

出力 (このディレクトリ): `v0c_d_design.csv` (60 点の x_m の内部決定)、`v0c_xm_scan_coverage.csv`、`v0c_gate.csv`
(ゲートの結果、設定ごとの縦持ち)、`v0c_vs_v0b.csv` (V0-b で成立した 31 点の合否の対照表)、`v0c_summary.json`
(判定・内訳・TOP・n 161 との一致・V0-b との対照・来歴)、`v0c_run.log` (標準出力の写し)、図 `v0c_gate_map.png` (git 管理外)。
V0-b の出力 (`v0b_*`) は読むだけで上書きしない。
実行: design/.venv-opt の Python で `python case/68.bell_two_mode/v0_moc_gate/run_v0c.py` (数分)。
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
from forge_design.geometry.wall_bell_c2 import (XM_BISECT_TOL, XM_SCAN_N, XM_SCAN_RANGE, BellC2Spec,  # noqa: E402
                                                BellC2Wall, cone15_length, design_bell_c2)

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
# V0-c の事前登録 (plan §6 V0-c)
V0B_CLASS = "中間"                   # V0-b の分類 (diagnostician 2026-10-11 で確定)
N161_DISAGREE_MAX = 2                # A: n 161 との不一致 ≤ 2
FLIP_B_MIN = 4                       # B: 既存 29 点のうち合否の反転 ≥ 4
CONTACT_REMOVED = ("D_L60_tm20_te10", "D_L70_tm25_te5")   # V0-b で成立した接触型 2 点 (#3d で不成立)
V0B_GATE = HERE / "v0b_gate.csv"
V0B_DESIGN = HERE / "v0b_d_design.csv"
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


def _read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(v):
    return None if v in (None, "") else float(v)


def main() -> int:
    t_all = time.time()
    log = Tee(HERE / "v0c_run.log")
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
                log(f"design {pid:20s} {d.status:10s} x_m {d.x_m if d.x_m is None else round(d.x_m, 6)}  "
                    f"区間 {rec['intervals'] or '-'} → 詰めた区間 {rec['interval_refined'] or '-'}  "
                    f"端の違反 [{rec['binding_lo'] or '-'} | {rec['binding_hi'] or '-'}]  詰めた端で確定 {rec['golden_at_refined_end']}  "
                    f"端 {d.edge_touch}  走査 (0,L) の成立 x_m [{cv['xm_ok_min']}, {cv['xm_ok_max']}]・"
                    f"範囲外 {cv['n_ok_below_scan'] + cv['n_ok_above_scan']}")
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
    # --- 集計 (V0-b と同じ規則) ---
    feas_by_L = {f"{int(round(fr * 100))}%": sum(1 for r in design_rows if r["design_status"] == "OK" and r["L_frac"] == fr)
                 for fr in L_FRACS}
    n_feas = sum(feas_by_L.values())
    valid = all(v >= VALID_MIN_PER_L for v in feas_by_L.values()) and n_feas >= VALID_MIN_TOTAL
    summary = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §6 V0-c (§5.1 #3d)", **git_info(),
               "gas": {"exhaust": str(EXHAUST.relative_to(ROOT)), "Tt_K": gas.Tt, "gamma_star": gas.gamma_star},
               "throat": {"R_t": R_T, "initial_line": "Hall throat characteristic (hall_initial_line)",
                          "gamma_hall": float(gas.gamma_throat())},
               "L_cone15_rt": Lc, "L_k_rt": [fr * Lc for fr in L_FRACS],
               "grid": {"theta_max_deg": THETA_MAX, "theta_e_deg": THETA_E, "L_frac": L_FRACS, "n_points": len(design_rows)},
               "xm_rule": {"scan_range_rt": XM_SCAN_RANGE, "scan_n": XM_SCAN_N, "bisect_tol_rt": XM_BISECT_TOL,
                           "feasibility_ii": "r'' ≥ −tol on [0, x_m] かつ r'' ≤ tol on [x_m, L] (tol = RESID_TOL、#3d)",
                           "select": "成立区間のうち J 最小の走査点を含む区間 → 区間の端を二分法で詰める → 黄金分割で J 最小"},
               "rule": {"valid": f"各 L_k に成立 ≥ {VALID_MIN_PER_L} かつ合計 ≥ {VALID_MIN_TOTAL}",
                        "A": f"通過率 ≥ {CLASS_A:.0%}", "B": f"≤ {CLASS_B:.0%}", "中間": "その間",
                        "pass_fraction_denominator": "成立点のうちゲートが判定できた点 (pass + reject)。判定不能は別に数える",
                        "marginal": f"通過点で C⁻ 間隔比の最小 < {MARGINAL_SEP}"},
               "v0c_rule": {"A": f"分類が {V0B_CLASS} のまま、かつ n 161 との不一致 ≤ {N161_DISAGREE_MAX} → V0 を閉じる",
                            "B": f"分類が変わる・有効条件が破れる・既存 29 点のうち合否が {FLIP_B_MIN} 点以上反転 → 再諮問",
                            "pass_fail": "gate_pass (pass = True、reject・判定不能 = False)",
                            "n161_disagree": "全点 (D + TOP) で n 321 と verdict が違う点の数 (V0-b と同じ数え方)"},
               "primary_config": PRIMARY}
    inf = {}
    for r in design_rows:
        if r["design_status"] == "OK":
            continue
        for part in r["scan_violations"].split("; "):
            if part:
                k, v = part.rsplit(" ", 1)
                inf[k] = inf.get(k, 0) + 1
    ok_rows = [r for r in design_rows if r["design_status"] == "OK"]
    binding = {}
    for r in ok_rows:
        for side in ("binding_lo", "binding_hi"):
            if r[side]:
                binding[f"{side[-2:]}: {r[side]}"] = binding.get(f"{side[-2:]}: {r[side]}", 0) + 1
    summary["d_design"] = {"n_feasible": n_feas, "feasible_by_L": feas_by_L,
                           "n_infeasible": len(design_rows) - n_feas,
                           "infeasible_points_by_violation_kind (走査点に 1 回でも現れた種類)": inf,
                           "n_edge_touch": sum(1 for r in design_rows if r["edge_touch"]),
                           "n_multi_interval": sum(1 for r in design_rows if (r["n_intervals"] or 0) >= 2),
                           "n_at_scan_point": sum(1 for r in ok_rows if r["golden_at_scan_point"]),
                           "n_at_refined_end": sum(1 for r in ok_rows if r["golden_at_refined_end"]),
                           "n_golden_fallback": sum(1 for r in ok_rows if r["golden_fallback"]),
                           "refined_end_binding (端の外側の違反の種類、点の数)": binding,
                           "x_m_range": [min(r["x_m"] for r in ok_rows), max(r["x_m"] for r in ok_rows)] if ok_rows else None}
    summary["xm_scan_coverage"] = {
        "dx_rt": COVER_DX, "n_points": len(cover_rows),
        "n_points_with_feasible_outside_scan": sum(1 for r in cover_rows if r["n_ok_below_scan"] + r["n_ok_above_scan"] > 0),
        "points_with_feasible_outside_scan": [r["point"] for r in cover_rows if r["n_ok_below_scan"] + r["n_ok_above_scan"] > 0],
        "feasible_xm_min_over_all": min((r["xm_ok_min"] for r in cover_rows if r["xm_ok_min"] is not None), default=None),
        "feasible_xm_max_over_all": max((r["xm_ok_max"] for r in cover_rows if r["xm_ok_max"] is not None), default=None),
        "n_points_feasible_in_fine_scan": sum(1 for r in cover_rows if r["n_cover_ok"] > 0),
        "points_feasible_in_fine_scan_but_design_infeasible": [r["point"] for r in cover_rows
                                                               if r["n_cover_ok"] > 0 and r["design_status"] != "OK"]}
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

        # 健全性の確認 (既知のパターン、V0-b と同じ)
        def pat(name, sel, want):
            pts = [r for r in d if sel(r)]
            bad = [f"{r['point']} ({r['verdict']})" for r in pts if r["verdict"] != want]
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
            summary["configs"][cfg]["n_differs_from_primary"] = len(diff)
            summary["configs"][cfg]["n_differs_from_primary_D"] = sum(1 for s in diff if s.startswith("D_"))
            summary["configs"][cfg]["agrees_with_primary"] = not diff
    summary["initial_lines"] = {cfg: {"n": len(ln["x"]), "x_axis_end": float(ln["x"][0]), "M_axis_end": float(ln["M"][0]),
                                      "M_wall": float(ln["M"][-1]), "thr_axis": float(ln["thr_axis"])}
                                for cfg, ln in lines.items()}
    # --- V0-b との対照 (V0-b で成立した 31 点) ---
    vb_gate = _read_csv(V0B_GATE)
    vb_des = {r["point"]: r for r in _read_csv(V0B_DESIGN)}
    vb = {(r["config"], r["point"]): r for r in vb_gate}
    vc = {(r["config"], r["point"]): r for r in gate_rows}
    vc_des = {r["point"]: r for r in design_rows}
    v0b_ok = [p for p, r in vb_des.items() if r["design_status"] == "OK"]
    keep29 = [p for p in v0b_ok if p not in CONTACT_REMOVED]
    cmp_rows = []
    for p in v0b_ok:
        b3, c3 = vb.get((PRIMARY, p)), vc.get((PRIMARY, p))
        b1, c1 = vb.get(("hall_n161", p)), vc.get(("hall_n161", p))
        row = {"point": p, "in_existing_29": p in keep29, "v0c_design_status": vc_des[p]["design_status"],
               "x_m_v0b": _f(vb_des[p]["x_m"]), "x_m_v0c": vc_des[p]["x_m"],
               "dx_m": (vc_des[p]["x_m"] - _f(vb_des[p]["x_m"])) if vc_des[p]["x_m"] is not None else None,
               "v0c_at_refined_end": vc_des[p]["golden_at_refined_end"], "v0c_binding_hi": vc_des[p]["binding_hi"],
               "verdict_n321_v0b": b3["verdict"] if b3 else None, "verdict_n321_v0c": c3["verdict"] if c3 else None,
               "gate_pass_n321_v0b": (b3["gate_pass"] == "True") if b3 else None,
               "gate_pass_n321_v0c": bool(c3["gate_pass"]) if c3 else None,
               "verdict_n161_v0b": b1["verdict"] if b1 else None, "verdict_n161_v0c": c1["verdict"] if c1 else None,
               "cminus_sep_min_n321_v0b": _f(b3["cminus_sep_ratio_min"]) if b3 else None,
               "cminus_sep_min_n321_v0c": c3["cminus_sep_ratio_min"] if c3 else None,
               "reversal_onset_x_n321_v0b": _f(b3["reversal_onset_x"]) if b3 else None,
               "reversal_onset_x_n321_v0c": c3["reversal_onset_x"] if c3 else None,
               "codes_n321_v0b": b3["codes"] if b3 else None, "codes_n321_v0c": c3["codes"] if c3 else None}
        if c3 is not None:
            row["flip_gate_pass"] = row["gate_pass_n321_v0b"] != row["gate_pass_n321_v0c"]
            row["flip_pass_reject"] = {row["verdict_n321_v0b"], row["verdict_n321_v0c"]} == {"pass", "reject"}
            row["verdict_changed"] = row["verdict_n321_v0b"] != row["verdict_n321_v0c"]
        else:
            row["flip_gate_pass"] = row["flip_pass_reject"] = row["verdict_changed"] = None
        cmp_rows.append(row)
    c29 = [r for r in cmp_rows if r["in_existing_29"]]
    n_flip = sum(1 for r in c29 if r["flip_gate_pass"])
    n_flip_pr = sum(1 for r in c29 if r["flip_pass_reject"])
    n_vchg = sum(1 for r in c29 if r["verdict_changed"])
    n29_missing = [r["point"] for r in c29 if r["v0c_design_status"] != "OK"]
    top_same = all(vb.get((cfg, p)) is not None and vb[(cfg, p)]["verdict"] == vc[(cfg, p)]["verdict"]
                   for (cfg, p) in vc if p.startswith("TOP_"))
    summary["v0b_comparison"] = {
        "v0b_feasible": len(v0b_ok), "existing_29": keep29, "contact_removed": list(CONTACT_REMOVED),
        "contact_removed_v0c_status": {p: vc_des[p]["design_status"] for p in CONTACT_REMOVED},
        "existing_29_not_ok_in_v0c": n29_missing,
        "new_ok_not_in_v0b": [r["point"] for r in design_rows if r["design_status"] == "OK" and r["point"] not in v0b_ok],
        "n_flip_gate_pass_29": n_flip, "n_flip_pass_reject_29": n_flip_pr, "n_verdict_changed_29": n_vchg,
        "flips_gate_pass_29": [f"{r['point']} ({r['verdict_n321_v0b']} → {r['verdict_n321_v0c']})" for r in c29
                               if r["flip_gate_pass"]],
        "verdict_changes_29": [f"{r['point']} ({r['verdict_n321_v0b']} → {r['verdict_n321_v0c']})" for r in c29
                               if r["verdict_changed"]],
        "dx_m_29": {"min": min(r["dx_m"] for r in c29 if r["dx_m"] is not None),
                    "max": max(r["dx_m"] for r in c29 if r["dx_m"] is not None),
                    "abs_max": max(abs(r["dx_m"]) for r in c29 if r["dx_m"] is not None)} if c29 else None,
        "top_verdicts_same_as_v0b": top_same}
    # --- 事前登録の判定 (plan §6 V0-c) ---
    p = summary["configs"][PRIMARY]
    c2 = summary["configs"]["hall_n161"]
    n161_all, n161_D = c2["n_differs_from_primary"], c2["n_differs_from_primary_D"]

    def judge(n_disagree, n_flips):
        b_reasons = []
        if not valid:
            b_reasons.append("有効条件が破れた")
        if p["classification"] != V0B_CLASS:
            b_reasons.append(f"分類が変わった ({V0B_CLASS} → {p['classification']})")
        if n_flips >= FLIP_B_MIN:
            b_reasons.append(f"既存 29 点のうち合否の反転 {n_flips} ≥ {FLIP_B_MIN}")
        a_ok = p["classification"] == V0B_CLASS and n_disagree <= N161_DISAGREE_MAX
        if b_reasons:
            return "B", b_reasons
        if a_ok:
            return "A", []
        return "A/B いずれにも当たらない", [f"n 161 との不一致 {n_disagree} > {N161_DISAGREE_MAX}"]
    jv, jr = judge(n161_all, n_flip)
    alt = {"n161_D_only": judge(n161_D, n_flip)[0], "flip_pass_reject_only": judge(n161_all, n_flip_pr)[0],
           "flip_any_verdict_change": judge(n161_all, n_vchg)[0]}
    summary["judgement"] = {
        "verdict": jv, "reasons": jr,
        "inputs": {"classification": p["classification"], "valid": valid, "pass_fraction": p["pass_fraction"],
                   "n161_disagree_all": n161_all, "n161_disagree_D": n161_D,
                   "n_flip_gate_pass_29": n_flip, "n_flip_pass_reject_29": n_flip_pr, "n_verdict_changed_29": n_vchg},
        "alternative_readings": alt, "reading_sensitive": any(v != jv for v in alt.values())}
    summary["elapsed_s"] = time.time() - t_all
    # --- 書き出し ---
    for name, rows in (("v0c_d_design.csv", design_rows), ("v0c_xm_scan_coverage.csv", cover_rows),
                       ("v0c_gate.csv", gate_rows), ("v0c_vs_v0b.csv", cmp_rows)):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
    (HERE / "v0c_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    try:
        plot_map(design_rows, gate_rows, Lc)
    except Exception as e:                      # 図は補助 (git 管理外)。失敗しても集計は残す
        log(f"図の作成に失敗: {e!r}")
    log(f"\nV0-c ({PRIMARY}): (D) 成立 {n_feas}/{len(design_rows)} {feas_by_L} → 有効条件 {'満たす' if valid else '満たさない'}・"
        f"通過 {p['d_pass']}・不受理 {p['d_reject']}・判定不能 {p['d_undecidable']}・通過率 {p['pass_fraction']}・"
        f"分類 {p['classification']}・marginal {p['marginal_pass']}・健全性 {'一致' if p['health_all_match'] else '不一致'}")
    for h in p["health_check"]:
        log(f"  健全性 {h['pattern']}: {({True: '一致', False: '不一致', None: '検査できない'})[h['match']]} "
            f"({h['n_points']} 点、不一致 {h['mismatch']})")
    log("  TOP: " + "  ".join(f"{k} {v['verdict']} (交差開始 x {v['reversal_onset_x']}, r {v['reversal_onset_r']})"
                              for k, v in p["top"].items()) + f"  (V0-b と同じ verdict: {top_same})")
    log(f"  n 161: 通過 {c2['d_pass']}・不受理 {c2['d_reject']}・判定不能 {c2['d_undecidable']} → {c2['classification']}・"
        f"主と違う点 {c2['verdict_differs_from_primary'] or 'なし'}")
    cmpv = summary["v0b_comparison"]
    log(f"  V0-b との対照 (既存 29 点、n 321): 合否 (gate_pass) の反転 {n_flip} {cmpv['flips_gate_pass_29']}・"
        f"pass↔reject {n_flip_pr}・verdict の変化 {n_vchg} {cmpv['verdict_changes_29']}・x_m の差 {cmpv['dx_m_29']}・"
        f"V0-c で不成立 {n29_missing or 'なし'}・接触型 2 点 {cmpv['contact_removed_v0c_status']}")
    cv = summary["xm_scan_coverage"]
    log(f"  走査範囲の確認: 成立する x_m が 0.5〜4.0 r_t の外にある点 {cv['n_points_with_feasible_outside_scan']}/{cv['n_points']}・"
        f"成立する x_m の全体の範囲 [{cv['feasible_xm_min_over_all']}, {cv['feasible_xm_max_over_all']}]・"
        f"0.01 刻みで成立があるのに INFEASIBLE {cv['points_feasible_in_fine_scan_but_design_infeasible'] or 'なし'}")
    log(f"  事前登録の判定 (plan §6 V0-c): {jv} {jr or ''}・読み方を変えた場合 {alt}  ({summary['elapsed_s']:.0f} s)")
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
    fig.suptitle(f"V0-c 前進 MOC ゲート ({PRIMARY}、#3d 後の (D)): 灰 = 不成立、緑 = 通過、赤 = 不受理、橙 = 判定不能",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(HERE / "v0c_gate_map.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
