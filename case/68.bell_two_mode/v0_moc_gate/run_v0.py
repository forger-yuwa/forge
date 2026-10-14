#!/usr/bin/env python3
"""§6 V0: 前進 MOC ゲートの A/B (CFD なし、事前登録の分類) — plan tooling-nozzle-bell-two-mode §6 V0・§5.1 #3・#3b。

**V0-a の記録用 (2026-10-11、判定不能 = 標本の退化)**。出力は `v0a_*` に改名して残した (`v0a_gate.csv`・
`v0a_d_feasibility.csv`・`v0a_summary.json`・`v0a_run_v0.log`、図 `v0a_gate_map.png`・`v0a_nets.png`)。§5.1 #3c の
ガードの後は `hall_n641` の設定が `hall_initial_line` で拒否されるので、このスクリプトはそのままでは再実行できない
(再現は commit `c636e716` で)。回し直し (V0-b) は `run_v0b.py`。

対象:
- (D) の格子点: θ_max {25, 32.5, 40}° × θ_e {5, 10, 15}° × x_m/L {0.1, 0.2, 0.3} を各 L_k
  (等面積比 15° 円錐の長さ (√ε − 1)/tan15° = 22.6575 r_t の 60・70・80・90・100 %、ε 50)。成立条件 (§5.1 #3) を満たす点だけゲートにかける。
- Phase 2 の TOP 最適形 (θ_n 25.8°・θ_e 10°) を ε 50・R_u = R_d = R_t = 1.0 に置き直したもの (`geometry.wall.NozzleWall`、
  `bell_type: top`) を各 L_k で。
ゲート: `geometry.moc_forward.forward_moc` (①の interior_vec、analytic・converge)、ガス = 基準ケースの凍結 TP
(`RocketExhaust.moc_gas`、排気ファイル ../exhaust_lox_lh2_pc100_of6.yaml)、初期線 = Hall のスロート特性線 (R = R_t = 1.0、
γ = スロートの局所 γ*)。

**事前登録の分類** (plan §6 V0): (D) の成立点のうちゲートを通った割合が 80 % 以上 → A、20 % 以下 → B、その間 → 中間。
分類を出すまで (解釈はしない)。

設定 (`CONFIGS`):
- `hall_n321` (主): Hall の初期線 321 点。分類はこれで出す。
- `hall_n161`: 解像度の確認 (主と判定が一致するか)。
- `hall_n641`: 細かい解像度の確認。**Hall の初期線 (R 1) は C⁻ の適合式を 1.1e-2 rad 満たさない** ため、この解像度では
  スロート直後で前線 1 が初期線と交差し、壁によらず不受理になる (記録用; 試験 `design/tests/run_moc_forward_tests.py` の info)。
- `ref_consistent_n321/641/1281` (**参考、plan 未承認**): Hall の θ 場の中で C⁻ を壁点から軸まで、ν を C⁻ の適合式で・幾何を
  tan(θ−μ(ν)) で同時に積分した「適合式と整合する初期線」。初期線の作り方を変える提案の根拠として記録するだけで、分類には使わない。

出力 (このディレクトリ): v0_d_feasibility.csv (D の 135 点の成立判定)、v0_gate.csv (ゲートの結果、設定ごとの縦持ち)、
v0_summary.json (分類・理由の内訳・TOP・来歴)、v0_gate_map.png、v0_nets.png。
実行: design/.venv-opt の Python で `python case/68.bell_two_mode/v0_moc_gate/run_v0.py` (数分)。
"""
from __future__ import annotations

import csv
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
from forge_design.geometry.moc_kernel import axis_theta_r, pm_mach_vec, pm_nu  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.wall import NozzleWall, WallParams  # noqa: E402
from forge_design.geometry.wall_bell_c2 import BellC2Spec, BellC2Wall, cone15_length  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
EPS = 50.0
R_T = 1.0
L_FRACS = (0.6, 0.7, 0.8, 0.9, 1.0)
THETA_MAX = (25.0, 32.5, 40.0)
THETA_E = (5.0, 10.0, 15.0)
XM_FRAC = (0.1, 0.2, 0.3)
TOP = {"theta_n_deg": 25.8, "theta_e_deg": 10.0}
PRIMARY = "hall_n321"
CONFIGS = {"hall_n321": ("hall", 321), "hall_n161": ("hall", 161), "hall_n641": ("hall", 641),
           "ref_consistent_n321": ("consistent", 321), "ref_consistent_n641": ("consistent", 641),
           "ref_consistent_n1281": ("consistent", 1281)}
CLASS_A, CLASS_B = 0.80, 0.20


def consistent_hall_line(R: float, gas, n: int) -> dict:
    """**参考 (plan 未承認)**: Hall の θ 場の中で、壁点 (0, 1) から C⁻ を軸まで積分する。ν は C⁻ の適合式
    d(θ+ν) = sinμ (sinθ/r) dσ、幾何は dx/dr = 1/tan(θ−μ(ν)) を同時に解く (r を独立変数、DOP853)。Hall の M は壁点だけ使う。
    戻りは `hall_initial_line` と同じ形 (軸端の θ_r は Hall の軸上の M′ と新しい軸端の M から)。"""
    from scipy.integrate import solve_ivp
    ht = HallThroat(R=float(R), gamma=float(gas.gamma_throat()))
    nu_w = float(pm_nu(float(ht.mach(0.0, 1.0)), gas))

    def rhs(r, y):
        x, tpn = y
        th = float(ht.theta(x, r))
        M = float(pm_mach_vec(np.array([tpn - th]), gas)[0])
        mu = np.arcsin(1.0 / M)
        dxdr = 1.0 / np.tan(th - mu)
        s = np.sin(th) / r if r > 1e-9 else 0.0
        return [dxdr, np.sin(mu) * s / np.cos(th - mu) * dxdr]
    r_end = 1e-7
    sol = solve_ivp(rhs, (1.0, r_end), [0.0, nu_w], method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
    rq = np.linspace(0.0, 1.0, int(n))
    y = sol.sol(np.maximum(rq, r_end))
    x = y[0]
    th = np.asarray(ht.theta(x, np.maximum(rq, r_end)), dtype=float)
    th[0] = 0.0
    nu = y[1] - th
    M = np.asarray(pm_mach_vec(nu, gas), dtype=float)
    thr = float(axis_theta_r(M[0], ht.axis_anchor(float(x[0]))[1], gas))
    return {"x": x, "r": rq, "M": M, "theta": th, "thr_axis": thr, "source": "ref_consistent_hall_theta"}


def git_info() -> dict:
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                           capture_output=True, text=True).stdout
    return {"commit": head, "design_tree_dirty": [ln for ln in dirty.splitlines() if ln.strip()]}


def classify(frac: float | None) -> str:
    if frac is None:
        return "判定不能 (成立点 0)"
    return "A" if frac >= CLASS_A else ("B" if frac <= CLASS_B else "中間")


def main() -> int:
    t_all = time.time()
    gas = RocketExhaust.from_file(EXHAUST).moc_gas()
    spec = BellC2Spec(eps=EPS, R_t=R_T)
    Lc = cone15_length(EPS)
    # --- (D) の成立判定 (135 点) ---
    feas_rows, d_walls = [], []
    for fr in L_FRACS:
        L = fr * Lc
        for tm in THETA_MAX:
            for te in THETA_E:
                for xf in XM_FRAC:
                    w = BellC2Wall(L, te, tm, xf, spec)
                    f = w.feasibility()
                    pid = f"D_L{int(round(fr * 100))}_tm{tm:g}_te{te:g}_xm{xf:g}"
                    feas_rows.append({"point": pid, "L_frac": fr, "L_rt": L, "theta_max_deg": tm, "theta_e_deg": te,
                                      "xm_frac": xf, "feasible": f["ok"], "violations": " / ".join(f["violations"]),
                                      "rp_min": f["rp_min"], "theta_sign_changes": f["theta_sign_changes"],
                                      "kappa_Rt_max": f["kappa_Rt_max"],
                                      "constraint_resid_max": float(np.max(np.abs(w.constraint_residuals())))})
                    if f["ok"]:
                        d_walls.append((pid, "D", fr, L, tm, te, xf, w))
    top_walls = []
    for fr in L_FRACS:
        L = fr * Lc
        nw = NozzleWall(WallParams(Ru=R_T, Rd=R_T, theta_a=np.deg2rad(TOP["theta_n_deg"]),
                                   theta_e=np.deg2rad(TOP["theta_e_deg"]), L_div=L, r_exit=np.sqrt(EPS), bell_type="top"))
        top_walls.append((f"TOP_L{int(round(fr * 100))}", "TOP", fr, L, TOP["theta_n_deg"], TOP["theta_e_deg"], None,
                          CallableWall(nw.r, nw.drdx)))
    # --- ゲート ---
    lines = {}
    gate_rows = []
    nets = {}
    for cfg, (kind, n) in CONFIGS.items():
        lines[cfg] = hall_initial_line(R_T, gas, n) if kind == "hall" else consistent_hall_line(R_T, gas, n)
        for (pid, fam, fr, L, a1, a2, a3, w) in d_walls + top_walls:
            res = forward_moc(w, L, gas, lines[cfg])
            rec = gate_record(res)
            row = {"config": cfg, "point": pid, "family": fam, "L_frac": fr, "L_rt": L,
                   "theta_max_or_n_deg": a1, "theta_e_deg": a2, "xm_frac": a3, **rec}
            gate_rows.append(row)
            if cfg == PRIMARY:
                nets[pid] = res
            print(f"{cfg:22s} {pid:28s} {'PASS' if rec['gate_pass'] else 'FAIL'} {rec['codes'] or '-'}  "
                  f"C⁻比 {rec['cminus_sep_ratio_min']}  交差開始 x {rec['reversal_onset_x']}  ({rec['time_s']:.1f} s)", flush=True)
    # --- 集計 ---
    summary = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §6 V0 (§5.1 #3・#3b)", **git_info(),
               "gas": {"exhaust": str(EXHAUST.relative_to(ROOT)), "Tt_K": gas.Tt, "gamma_star": gas.gamma_star},
               "throat": {"R_t": R_T, "initial_line": "Hall throat characteristic (HallThroat.throat_characteristic)",
                          "gamma_hall": float(gas.gamma_throat())},
               "L_cone15_rt": Lc, "L_k_rt": [fr * Lc for fr in L_FRACS],
               "rule": {"A": f"成立点の {CLASS_A:.0%} 以上が通る", "B": f"{CLASS_B:.0%} 以下", "中間": "その間"},
               "primary_config": PRIMARY, "d_grid_points": len(feas_rows), "d_feasible": len(d_walls),
               "d_infeasible_reasons": {}, "configs": {}}
    inf = {}
    for r in feas_rows:
        if r["feasible"]:
            continue
        for v in r["violations"].split(" / "):
            key = "r' < 0" if v.startswith("r' < 0") else ("θ' の符号変化 ≠ 1" if "符号変化" in v else ("|κ| > 1/R_t" if "|κ|" in v else v))
            inf[key] = inf.get(key, 0) + 1
    summary["d_infeasible_reasons"] = inf
    summary["d_feasible_by_L"] = {f"{int(round(fr * 100))}%": sum(1 for r in feas_rows if r["feasible"] and r["L_frac"] == fr)
                                  for fr in L_FRACS}
    for cfg in CONFIGS:
        rows = [r for r in gate_rows if r["config"] == cfg]
        d = [r for r in rows if r["family"] == "D"]
        npass = sum(r["gate_pass"] for r in d)
        frac = npass / len(d) if d else None
        reasons = {}
        for r in d:
            for c in (r["codes"].split(";") if r["codes"] else []):
                reasons[c] = reasons.get(c, 0) + 1
        summary["configs"][cfg] = {
            "initial_line": CONFIGS[cfg][0], "n_line": CONFIGS[cfg][1],
            "d_feasible": len(d), "d_pass": npass, "d_pass_fraction": frac, "classification": classify(frac),
            "d_fail_codes": reasons,
            "d_points": {r["point"]: {"pass": r["gate_pass"], "codes": r["codes"],
                                      "cminus_sep_ratio_min": r["cminus_sep_ratio_min"],
                                      "reversal_onset_x": r["reversal_onset_x"], "reversal_onset_r": r["reversal_onset_r"]}
                         for r in d},
            "top": {r["point"]: {"pass": r["gate_pass"], "codes": r["codes"],
                                 "cminus_sep_ratio_min": r["cminus_sep_ratio_min"],
                                 "reversal_onset_x": r["reversal_onset_x"], "reversal_onset_r": r["reversal_onset_r"]}
                    for r in rows if r["family"] == "TOP"}}
        if cfg != PRIMARY:
            prim = {r["point"]: r["gate_pass"] for r in gate_rows if r["config"] == PRIMARY}
            summary["configs"][cfg]["verdict_differs_from_primary"] = sorted(
                r["point"] for r in rows if r["gate_pass"] != prim.get(r["point"]))
    line_info = {}
    for cfg, ln in lines.items():
        line_info[cfg] = {"x_axis_end": float(ln["x"][0]), "M_axis_end": float(ln["M"][0]), "M_wall": float(ln["M"][-1]),
                          "thr_axis": float(ln["thr_axis"])}
    summary["initial_lines"] = line_info
    summary["elapsed_s"] = time.time() - t_all
    # --- 書き出し ---
    with open(HERE / "v0_d_feasibility.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(feas_rows[0]))
        wr.writeheader()
        wr.writerows(feas_rows)
    with open(HERE / "v0_gate.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(gate_rows[0]))
        wr.writeheader()
        wr.writerows(gate_rows)
    (HERE / "v0_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    plot_map(feas_rows, gate_rows, Lc)
    plot_nets(nets, d_walls, top_walls)
    p = summary["configs"][PRIMARY]
    print(f"\nV0 ({PRIMARY}): (D) 成立 {summary['d_feasible']}/{summary['d_grid_points']}・通過 {p['d_pass']}・"
          f"割合 {p['d_pass_fraction']}・分類 {p['classification']}  ({summary['elapsed_s']:.0f} s)")
    for cfg, c in summary["configs"].items():
        print(f"  {cfg:22s} 通過 {c['d_pass']}/{c['d_feasible']} → {c['classification']}  TOP 通過 "
              f"{sum(v['pass'] for v in c['top'].values())}/{len(c['top'])}  主と違う点 {c.get('verdict_differs_from_primary', '-')}")
    return 0


def _font():
    import matplotlib
    from matplotlib import font_manager
    for f in (Path.home() / ".fonts").glob("NotoSansCJKjp-Regular.otf"):
        font_manager.fontManager.addfont(str(f))
        matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"


def plot_map(feas_rows, gate_rows, Lc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _font()
    prim = {r["point"]: r for r in gate_rows if r["config"] == PRIMARY}
    fig, axs = plt.subplots(len(L_FRACS), len(XM_FRAC), figsize=(11, 15), squeeze=False)
    for i, fr in enumerate(L_FRACS):
        top = prim.get(f"TOP_L{int(round(fr * 100))}")
        for jx, xf in enumerate(XM_FRAC):
            ax = axs[i, jx]
            for a, tm in enumerate(THETA_MAX):
                for b, te in enumerate(THETA_E):
                    pid = f"D_L{int(round(fr * 100))}_tm{tm:g}_te{te:g}_xm{xf:g}"
                    fe = next(r for r in feas_rows if r["point"] == pid)
                    if not fe["feasible"]:
                        col, txt = "0.85", "不成立"
                    else:
                        g = prim[pid]
                        col = "#7fc97f" if g["gate_pass"] else "#f08080"
                        txt = f"{'通過' if g['gate_pass'] else '不受理'}\nC-比 {g['cminus_sep_ratio_min']:.2f}"
                    ax.add_patch(plt.Rectangle((b, a), 1, 1, color=col, ec="w"))
                    ax.text(b + 0.5, a + 0.5, txt, ha="center", va="center", fontsize=7)
            ax.set_xlim(0, 3)
            ax.set_ylim(0, 3)
            ax.set_xticks([0.5, 1.5, 2.5], [f"θe {t:g}°" for t in THETA_E], fontsize=8)
            ax.set_yticks([0.5, 1.5, 2.5], [f"θmax {t:g}°" for t in THETA_MAX], fontsize=8)
            ttl = f"L = {fr:.0%} ({fr * Lc:.2f} r_t)・x_m/L {xf:g}"
            if jx == 0 and top is not None:
                ttl += f"\nTOP: {'通過' if top['gate_pass'] else '不受理'} (C-比 {top['cminus_sep_ratio_min']:.2f})"
            ax.set_title(ttl, fontsize=9)
    fig.suptitle(f"V0 前進 MOC ゲート ({PRIMARY}): 灰 = 成立条件で不成立、緑 = 通過、赤 = 不受理", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(HERE / "v0_gate_map.png", dpi=110)
    plt.close(fig)


def plot_nets(nets, d_walls, top_walls):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from forge_design.geometry.moc_forward import _signed_len
    from forge_design.geometry.moc_kernel import mu_vec
    _font()
    picks = []
    passed = [w for w in d_walls if nets[w[0]]["gate"]["pass"]]
    failed = [w for w in d_walls if not nets[w[0]]["gate"]["pass"]]
    if passed:
        picks.append(passed[0])
    if failed:
        picks.append(failed[0])
    picks.append(top_walls[2])
    fig, axs = plt.subplots(len(picks), 1, figsize=(12, 3.6 * len(picks)), squeeze=False)
    for ax, (pid, fam, fr, L, a1, a2, a3, w) in zip(axs[:, 0], picks):
        res = nets[pid]
        net = res["net"]
        X, R, TH, M = (net.a[f] for f in ("x", "r", "th", "M"))
        ok = (net.status == 0) & (X <= L)
        K, n1 = X.shape
        for k in range(0, K, max(1, K // 50)):
            ax.plot(X[k, ok[k]], R[k, ok[k]], "b-", lw=0.35)
        for a in range(-(n1 - 1), K, max(1, (K + n1) // 50)):
            kk = np.arange(max(a, 0), min(K, a + n1))
            jj = kk - a
            m = ok[kk, jj]
            ax.plot(X[kk, jj][m], R[kk, jj][m], "r-", lw=0.35)
        xs = np.linspace(0.0, L, 600)
        ax.plot(xs, np.asarray(w.r(xs, 0)), "k-", lw=1.4)
        MU = mu_vec(np.where(np.isfinite(M), M, 1.0))
        ang = 0.5 * (TH[:-1, :-1] + TH[1:, 1:]) + 0.5 * (MU[:-1, :-1] + MU[1:, 1:])
        both = ok[1:, 1:] & ok[:-1, :-1]
        rev = both & ~(_signed_len(X[1:, 1:], X[:-1, :-1], ang) > 0.0)
        ax.plot(X[1:, 1:][rev], R[1:, 1:][rev], "x", color="g", ms=3, label="網の反転 (C+ 沿い)")
        g = res["gate"]
        lab = (f"{pid}: {'通過' if g['pass'] else '不受理'} {g['codes'] or ''}  "
               f"C- 間隔比の最小 {g['indicators'].get('cminus_sep_ratio_min', float('nan')):.3f}")
        ax.set_title(lab, fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlim(0, L * 1.02)
        ax.set_ylim(0, np.sqrt(EPS) * 1.05)
        if np.any(rev):
            ax.legend(fontsize=8, loc="upper left")
    axs[-1, 0].set_xlabel("x / r_t (青 = C- 前線、赤 = C+、間引き表示)")
    fig.tight_layout()
    fig.savefig(HERE / "v0_nets.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
