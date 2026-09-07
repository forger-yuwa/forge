r"""排除厚さの固定点反復ドライバ (固定 Euler 基準・コア整合抽出)。

計画: plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5–4.6。

1 pass = 前 pass の NS run から $\delta_r(x)$ を抽出 (`metrics.deltastar.deltastar_from_core_matched_euler`)
→ 緩和して次の入力 $\delta_{in}^{k+1} = (1-\omega)\delta_{in}^k + \omega\,\delta_{use}^k$ (hard 不合格の断面は前回値保持)
→ 半径方向オフセットの物理壁で `prepare_ns` → 段階起動 NS → `collect` (質量流量帳簿込み)
→ 新 run からも抽出して固定点差を記録 (`fixed_point.json`)。

使い方 (リポジトリルートで):
  design/.venv-opt/bin/python -m forge_design.feedback.deltastar_loop \
      --problem case/45.isobutane_m6_d155/problem_d155_ns.yaml \
      --euler-ref case/45.isobutane_m6_d155/run_0001_euler_shortest_dry \
      --prev case/45.isobutane_m6_d155/run_0004_ns_v3 \
      --run-dir case/45.isobutane_m6_d155/run_0019_ns_cm_pass1 [--omega 0.5] [--prepare-only]

推奨 (2026-09-04, case/45 run_0026 で検証): 収束済み NS 場からの warm start では
  `--stages none --cfl 5 --implicit-relax 0.7 --steps 12000`
(soft/mid 6000 step は不要、cfl 5 + relax 0.7 で 5000〜9000 step で cfl1/24000 step と同じ残差水準・出口 M は 8000 step で凍結。
 NS 1 本 ≈ 95 s)。cold start (Euler → 中継 → 初回 NS) では `--stages full` (既定) か `--stages ramp --ramp 1,2,3.5`。
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..metrics.deltastar import deltastar_from_core_matched_euler, massflow_ratio


def extract_and_merge(prev_run, euler_run, omega: float = 0.5, smooth_lam: float = 1.0,
                      knot_spacing: float = 2.0, **kw) -> dict:
    """前 pass の NS run から抽出し、次 pass の入力 δ_r(x) を作る。
    出力 (prev_run 内): delta_r_equiv.csv / delta_r_equiv_diag.json / delta_r_next.csv。"""
    prev_run = Path(prev_run)
    d = deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=prev_run, **kw)
    d_in = d["delta_in"]
    d_use = d["delta_r_use"]
    held = ~np.isfinite(d_use)
    # 平滑化 (2026-09-04 ユーザ指摘: 3 次平滑化では壁曲率が凸凹): 抽出 raw を 5 次 P-spline (3 階差分ペナルティ,
    # ノット 2 r_t, λ=1) で「曲率が滑らか」な δ_ext(x) にしてから緩和する。hard 不合格の点は重み 0。
    from ..metrics.deltastar import smooth_delta_quintic
    x = d["x"]; r_inv = d["r_wall_euler"]
    raw = d["delta_r_raw"]; wts = d["hard_ok"].astype(float) & np.isfinite(raw) if False else (d["hard_ok"] & np.isfinite(raw)).astype(float)
    lam_used = None; diag = None
    for lam in (smooth_lam, smooth_lam * 10, smooth_lam * 100, smooth_lam * 1000):
        f_s, diag = smooth_delta_quintic(x, raw, weights=wts, knot_spacing=knot_spacing, lam=lam)
        d_ext = f_s(x)
        d_next = (1.0 - omega) * d_in + omega * d_ext
        if omega < 1.0:
            # 前回入力 δ_in (旧方式の壁など) の凸凹を引き継がないよう、緩和後も同じ P-spline を通す
            f_b, _ = smooth_delta_quintic(x, d_next, knot_spacing=knot_spacing, lam=lam)
            d_next = f_b(x)
        # 単調性ガード: 新しい物理壁 r_inv + δ_next がスロート下流で非単調なら λ を 10 倍して再平滑化
        r_new = r_inv + d_next
        m = x >= 0.5
        if np.all(np.diff(r_new[m]) > -1e-9):
            lam_used = lam; break
    else:
        lam_used = f"{lam} (still non-monotone)"
    d_use = d_ext                                   # 以降の帳簿は平滑化後の抽出値で取る
    held = ~np.isfinite(d["delta_r_use"])
    np.savetxt(prev_run / "delta_r_next.csv", np.c_[d["x"], d_next, d_in, d_use, held.astype(int)],
               delimiter=",", comments="",
               header=f"x_rt,delta_r,delta_in_prev,delta_r_use_smoothed,held (omega={omega}; quintic P-spline knot={knot_spacing} lam={lam_used}; resid_rel_rms={diag['resid_rel_rms']:.4f})")
    fin = np.isfinite(d_use) & (d_in > 1e-4)
    ratio = d_use[fin] / d_in[fin]
    summary = {"omega": omega, "smooth": {"kind": "quintic_pspline", "knot_spacing": knot_spacing, "lam": str(lam_used), **diag},
               "n_stations": int(len(d["x"])), "n_ok": int(d["ok"].sum()),
               "n_hard_ok": int(d["hard_ok"].sum()), "n_held": int(held.sum()),
               "use_over_in_median": float(np.median(ratio)) if fin.any() else None,
               "use_over_in_p10_p90": [float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))] if fin.any() else None,
               "max_abs_change": float(np.nanmax(np.abs(d_use - d_in))) if fin.any() else None,
               "massflow": d["massflow"],
               "delta_r_throat_use": float(np.interp(0.0, d["x"], np.where(np.isfinite(d_use), d_use, d_in))),
               "delta_r_throat_in": float(np.interp(0.0, d["x"], d_in))}
    (prev_run / "delta_r_extract_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def solve_rt(problem, R_exit_m: float, prev_run=None, euler_run=None, n_iter: int = 6) -> dict:
    r"""出口の物理半径 $R$ を仕様に合わせる **スロート半径 $r_t$ の 1 変数解**。

    設計は $r_t$ 無次元で不変なので $R = r_t\,[r_F/r_t + \delta_r(x_F)/r_t]$ の $r_t$ だけを解く。
    - prev_run なし: 積分法 (CONTUR) の $\delta_{r}(x_F; r_t)$ で Newton (CFD 前に使う)。
    - prev_run あり: その run の抽出 δ_r(x_F) を使い、$r_t$ 依存は $Re^{-0.2}$ で補正 (NS 後の最終補正; 再計算不要)。
    戻り: dict(r_t_m, delta_exit_rt, source, iters)。"""
    from ..probdef import load_problem
    from ..evaluate.runner_axismach import design_chain, _gam_or_gas
    from ..feedback.deltastar_integral import integral_bl
    p = load_problem(problem); d = design_chain(p)
    S0 = float(p.spec["r_throat"]); Pt = float(p.spec["Pt"]); Tt = float(p.spec["Tt"])
    rF = float(d["wall_inv"][-1, 1]); xF = float(d["wall_inv"][-1, 0])
    hist = []
    if prev_run is None:
        rt = S0
        for k in range(n_iter):
            de = float(integral_bl(d["wall"], d["wall_inv"], _gam_or_gas(p), p.cp, Pt, Tt, rt)["delta_r"][-1])
            rt_new = R_exit_m / (rF + de); hist.append((rt, de, rt_new))
            if abs(rt_new - rt) < 1e-7: rt = rt_new; break
            rt = rt_new
        src = "integral_bl (CONTUR)"
    else:
        prev_run = Path(prev_run)
        if not (prev_run / "delta_r_equiv.csv").exists():
            if euler_run is None: raise ValueError("prev_run に抽出結果が無く euler_run も未指定")
            deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=prev_run)
        # 次 pass の壁に実際に載る値 (P-spline 平滑化後 = delta_r_next.csv の x_F 端) を使う。無ければ生抽出の x_F−0.3
        if (prev_run / "delta_r_next.csv").exists():
            nx = np.loadtxt(prev_run / "delta_r_next.csv", delimiter=",", skiprows=1)
            d_meas = float(np.interp(xF, nx[:, 0], nx[:, 1]))
        else:
            e = np.genfromtxt(prev_run / "delta_r_equiv.csv", delimiter=",", names=True)
            d_meas = float(np.interp(xF - 0.3, e["x_rt"], e["delta_r_raw"]))
        S_prev = float(json.loads((prev_run / "prepare_info.json").read_text())["scale_m"])
        rt = S_prev
        for k in range(n_iter):
            de = d_meas * (rt / S_prev) ** -0.2
            rt_new = R_exit_m / (rF + de); hist.append((rt, de, rt_new)); rt = rt_new
        src = f"measured delta_r from {prev_run.name} (Re^-0.2 scaling)"
    return dict(r_t_m=float(rt), r_t_prev_m=S0, delta_exit_rt=float(hist[-1][1]), r_F_rt=rF, x_F_rt=xF,
                R_exit_m=R_exit_m, source=src, iters=hist, length_m=float(xF * rt))


def run_pass(problem, euler_ref, prev_run, run_dir, omega: float = 0.5, ic_from=None,
             prepare_only: bool = False, nsteps=None, cfl_main=None, implicit_relax=None,
             stages: str = "full", ramp=None, ramp_steps: int = 1000) -> dict:
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    prev_run = Path(prev_run); run_dir = Path(run_dir)
    summ = extract_and_merge(prev_run, euler_ref, omega=omega)
    print("extract(prev):", json.dumps({k: v for k, v in summ.items() if k != "massflow"}), flush=True)
    print("massflow(prev):", json.dumps(summ["massflow"]), flush=True)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from or prev_run,
                      delta_r_csv=prev_run / "delta_r_next.csv", offset="radial",
                      euler_ref=euler_ref, omega=omega, prev_run=prev_run,
                      cfl_main=cfl_main, implicit_relax=implicit_relax)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "mesh", "nStepOuter", "cfl_main")}, indent=1), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    # 固定点差: 新 run から抽出し、入力 (= 新 run の壁 − Euler 壁) と比較
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("fixed_point(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    return m


def run_pass0_integral(problem, euler_ref, run_dir, ic_from, initializer=None, prepare_only: bool = False,
                       nsteps=None, omega: float = 0.5, cfl_main=None, implicit_relax=None,
                       stages: str = "full", ramp=None, ramp_steps: int = 1000) -> dict:
    """pass 0: 積分法 (CONTUR) 初期壁で NS を立て、終了後に抽出して次 pass 用 delta_r_next.csv まで作る。
    initializer=None なら problem YAML の `deltastar_initializer` (無ければ断熱 contur)。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    run_dir = Path(run_dir)
    init = initializer if initializer is not None else {"model": "contur", "thermal_bc": {"mode": "adiabatic"}}
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from, initializer=init,
                      euler_ref=euler_ref, omega=omega, cfl_main=cfl_main, implicit_relax=implicit_relax)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "initializer", "mesh", "nStepOuter", "cfl_main")},
                     indent=1, default=str), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("extract(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    return m


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="δ_r 固定点反復 (固定 Euler 基準)")
    ap.add_argument("--problem", required=True)
    ap.add_argument("--euler-ref", required=True)
    ap.add_argument("--prev", default=None, help="前 pass の NS run (抽出元・IC 既定)。--init-integral では不要")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--omega", type=float, default=0.5)
    ap.add_argument("--ic-from", default=None)
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--extract-only", action="store_true", help="抽出と delta_r_next.csv だけ作る")
    ap.add_argument("--init-integral", action="store_true",
                    help="pass 0: 積分法 (CONTUR) 初期壁で NS を立てる (YAML の deltastar_initializer を使う)")
    ap.add_argument("--init-thermal", default=None, help="--init-integral の熱境界条件 JSON (例 '{\"mode\":\"adiabatic\"}')")
    ap.add_argument("--cfl", type=float, default=None, help="本段 cfl (YAML evaluate.cfl_main を上書き)")
    ap.add_argument("--implicit-relax", type=float, default=None, help="implicitRelax (YAML evaluate.implicit_relax を上書き)")
    ap.add_argument("--stages", default="full", choices=("full", "none", "ramp"), help="起動: full=soft/mid/本段, none=本段のみ, ramp=cfl を段階的に上げて本段")
    ap.add_argument("--ramp", default=None, help="--stages ramp の cfl 列 (例 '1,2,3.5')")
    ap.add_argument("--ramp-steps", type=int, default=1000)
    ap.add_argument("--solve-rt", type=float, default=None, metavar="R_EXIT_M",
                    help="出口物理半径 [m] を与えて r_t を解くだけ (--prev があれば実測 δ_r で、無ければ積分法で)")
    a = ap.parse_args(argv)
    if a.solve_rt is not None:
        r = solve_rt(a.problem, a.solve_rt, prev_run=a.prev, euler_run=a.euler_ref)
        print(json.dumps(r, indent=1, default=str))
        return 0
    if a.extract_only:
        s = extract_and_merge(a.prev, a.euler_ref, omega=a.omega)
        print(json.dumps(s, indent=1))
        return 0
    if a.init_integral:
        init = None
        if a.init_thermal:
            init = {"model": "contur", "thermal_bc": json.loads(a.init_thermal)}
        ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
        run_pass0_integral(a.problem, a.euler_ref, a.run_dir, a.ic_from, initializer=init,
                           prepare_only=a.prepare_only, nsteps=a.steps, omega=a.omega, cfl_main=a.cfl,
                           implicit_relax=a.implicit_relax, stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps)
        return 0
    if not a.prev:
        ap.error("--prev が必要 (--init-integral でなければ)")
    ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
    run_pass(a.problem, a.euler_ref, a.prev, a.run_dir, omega=a.omega, ic_from=a.ic_from,
             prepare_only=a.prepare_only, nsteps=a.steps, cfl_main=a.cfl, implicit_relax=a.implicit_relax,
             stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
