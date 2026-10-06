"""Euler A/B (plan tooling-nozzle-throat-monotone-r2 §6 E1、§5.1 #5) の run 準備・実行と、壁の証拠 (M4)・IC 写像の記録。
腕 A (pinG1) = 現行壁 problem_d155_euler_pin_G1_recal.yaml / 腕 B (monoG1) = 単調壁 problem_d155_euler_pin_G1_recal_mono.yaml
(G1_recal + geometry.wall_fit_mono_r2 [0, 1.5])。生産 Euler 格子 G1、IC はどちらも run_0114 の最終場:
  A は同一メッシュなので restart_field (保存量の index コピー、ビット一致検査つき)、
  B は壁節点が最大 0.5 µm 動くので interp_field (規則どおり) で、IC 写像 (対応節点の一致率・最大移動量・保存量の差) を記録する。
起動は euler_grid_ab.py prep と同じ RA.prepare + RA.run_staged: soft (1 次 cfl 0.5、3000) → 本段 2 次 cfl 2・implicitRelax 0.7・
18000 step・1000 ごと出力。
usage: [CASE_RUNS=<run_0062 / run_0114 のある case dir>] python3 throat_mono_ab.py prep <prep_dir> {pinG1|monoG1} [--no-ic (乾式確認)]
       python3 throat_mono_ab.py run <run_dir>               (準備済み run を soft → 本段で回す; forge は run_case.sh 経由)
       python3 throat_mono_ab.py evidence <run_dir> [<run_dir_other_arm>]   (壁の証拠だけ出す; JSON を標準出力へ)
壁の証拠 (M4): prepare_info.json の wall_fit.spline (当てはめ後のノット・係数・次数) を、変換後メッシュ (nozzle.h5, float32) の
壁節点の x で評価し、節点の r と照合する (許容 = float32 の丸め ulp(r) + |r′|·ulp(x))。wall_design.csv は当てはめ前の MOC 点群なので
証拠にならない。r″ の山の有無は壁節点からは認定しない (解析形で形状ゲート S1 が担う)。
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import BSpline

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
sys.path.insert(0, str(ROOT / "design"))
TOOLS = ROOT / "solver_density_cuda/tools"
RUNS = Path(os.environ.get("CASE_RUNS", C))           # 凍結源 run_0062・IC run_0114 の場所 (既定は case dir)
IC_RUN = "run_0114_euler_pin_G1_recal_ext6k"
ARMS = {"pinG1": "problem_d155_euler_pin_G1_recal.yaml", "monoG1": "problem_d155_euler_pin_G1_recal_mono.yaml"}
IC_MODE = {"pinG1": "restart_field", "monoG1": "interp_field"}
NSTEPS, CFL_MAIN, RELAX = 18000, 2.0, 0.7
CONS = ("ro", "roUx", "roUy", "roUz", "roe")


def _last_res(run: Path) -> Path:
    fs = sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not fs:
        raise FileNotFoundError(f"{run} に res_*.h5 が無い")
    return fs[-1]


def wall_spline(info: dict):
    """prepare_info の wall_fit.spline → BSpline (無ければ None)。"""
    sp = (info.get("wall_fit") or {}).get("spline")
    if not sp:
        return None
    return BSpline(np.asarray(sp["t"], dtype=float), np.asarray(sp["c"], dtype=float), int(sp["k"]))


def wall_nodes(run_dir: Path, info: dict):
    """nozzle.h5 の壁節点 (各 i 列の j = nj−1) の x, r [m] (float32 の値を float64 で)。並びと壁 BC を検査する。"""
    ni = int(info["mesh"]["ni"])
    with h5py.File(run_dir / "nozzle.h5") as f:
        raw = f["/MESH/COORD"][:]
        nc = raw.reshape(-1, 3)
        wall_ids = None
        for k in f["BCONDS"]:
            if str(f["BCONDS"][k].attrs.get("bcondKind", "")) in ("slip", "wall", "noslip", "wall_noSlip"):
                wall_ids = np.unique(np.asarray(f["BCONDS"][k]["vizBfaceNodes"]))
    if nc.shape[0] % ni:
        raise ValueError(f"{run_dir}: 節点数 {nc.shape[0]} が ni {ni} で割り切れない")
    nj = nc.shape[0] // ni
    X, Y = nc[:, 0].reshape(ni, nj), nc[:, 1].reshape(ni, nj)
    if not (np.all(X == X[:, :1]) and np.array_equal(Y[:, -1], Y.max(axis=1))):
        raise ValueError(f"{run_dir}: 節点の並びが (i, j) の構造格子でない (列の x が一定でない・壁が j = nj−1 でない)")
    if wall_ids is not None:
        idx = np.arange(nc.shape[0]).reshape(ni, nj)[:, -1]
        if not np.all(np.isin(idx, wall_ids)):
            raise ValueError(f"{run_dir}: j = nj−1 の節点が壁 BC の節点に含まれない")
    return X[:, -1].astype(float), Y[:, -1].astype(float), str(raw.dtype)


def wall_evidence(run_dir, other_dir=None) -> dict:
    """M4: 当てはめ後 spline を変換後メッシュの壁節点の x で評価し、節点の r と照合する。other_dir (他腕の run) があれば
    実座標差・解析差・float32 丸め後の差と、各節点が自腕の spline と他腕の spline のどちらに一致するかを記録する。"""
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    S = float(info["scale_m"])
    spl = wall_spline(info)
    out = {"run": run_dir.name, "mono_r2": (info.get("wall_fit") or {}).get("mono_r2")}
    if spl is None:
        out.update(status="missing", reason="prepare_info.json の wall_fit.spline が無い (当てはめ後の壁が保存されていない)")
        return out
    xw, rw, dtype = wall_nodes(run_dir, info)
    m = (xw / S >= spl.t[0]) & (xw / S <= spl.t[-1])
    xs = xw[m] / S
    r_an = S * spl(xs)
    tol = (np.spacing(np.abs(rw[m]).astype(np.float32)).astype(float)
           + np.abs(spl(xs, 1)) * np.spacing(np.abs(xw[m]).astype(np.float32)).astype(float))
    res = rw[m] - r_an
    ratio = np.abs(res) / tol
    out.update(coord_dtype=dtype, n_wall_nodes_design=int(m.sum()), resid_max_m=float(np.abs(res).max()),
               x_resid_max_rt=float(xs[np.argmax(np.abs(res))]), resid_over_f32tol_max=float(ratio.max()),
               f32_tol_max_m=float(tol.max()),
               status=("consistent" if ratio.max() <= 1.0 else "mismatch"))
    if other_dir is not None:
        other_dir = Path(other_dir)
        io = json.loads((other_dir / "prepare_info.json").read_text())
        so = wall_spline(io)
        xo, ro, _ = wall_nodes(other_dir, io)
        if so is None:
            out["vs_other"] = {"run": other_dir.name, "status": "missing", "reason": "他腕の wall_fit.spline が無い"}
        elif not np.array_equal(xo, xw):
            out["vs_other"] = {"run": other_dir.name, "status": "x_differs", "max_dx_m": float(np.abs(xo - xw).max())}
        else:
            ro_an = S * so(xs)
            d_an = r_an - ro_an
            d_f32 = r_an.astype(np.float32).astype(float) - ro_an.astype(np.float32).astype(float)
            d_act = rw[m] - ro[m]
            disc = np.abs(d_an) > 2.0 * tol                 # 2 本の解析壁が float32 で見分けられる節点
            own_better = np.abs(rw[m] - r_an) < np.abs(rw[m] - ro_an)
            out["vs_other"] = {"run": other_dir.name, "status": "ok",
                               "actual_dr_max_m": float(np.abs(d_act).max()), "x_actual_dr_max_rt": float(xs[np.argmax(np.abs(d_act))]),
                               "analytic_dr_max_m": float(np.abs(d_an).max()), "x_analytic_dr_max_rt": float(xs[np.argmax(np.abs(d_an))]),
                               "f32_rounded_analytic_dr_max_m": float(np.abs(d_f32).max()),
                               "actual_minus_f32_analytic_max_m": float(np.abs(d_act - d_f32).max()),
                               "n_discriminable": int(disc.sum()),
                               "n_discriminable_matching_own": int((disc & own_better).sum()),
                               "x_range_discriminable_rt": ([float(xs[disc].min()), float(xs[disc].max())] if disc.any() else None)}
    return out


def ic_map_record(src_res: Path, src_mesh: Path, dst_h5: Path) -> dict:
    """B の IC 写像 (interp_field の後): 対応節点の一致率・最大移動量、保存量の差 (最大・RMS、index 対応で比較)。"""
    from scipy.spatial import cKDTree
    with h5py.File(src_mesh) as f:
        cs = f["/MESH/COORD"][:].reshape(-1, 3).astype(float)
    with h5py.File(dst_h5) as f:
        cd = f["/MESH/COORD"][:].reshape(-1, 3).astype(float)
        vd = {k: f["/VALUE/" + k][:].astype(float) for k in CONS if "/VALUE/" + k in f}
    with h5py.File(src_res) as f:
        vs = {k: f["/VALUE/" + k][:].astype(float) for k in CONS if "/VALUE/" + k in f}
    out = {"src_res": str(src_res), "n_src": int(len(cs)), "n_dst": int(len(cd))}
    dist, idx = cKDTree(cs).query(cd)
    out["nn_match_rate"] = (float(np.mean(idx == np.arange(len(cd)))) if len(cs) == len(cd) else None)
    out["nn_dist_max_m"] = float(dist.max())
    if len(cs) == len(cd):
        disp = np.linalg.norm(cd - cs, axis=1)
        out["same_index_disp_max_m"] = float(disp.max())
        out["n_moved_nodes"] = int(np.count_nonzero(disp > 0))
        dv = {}
        for k in CONS:
            if k in vd and k in vs and len(vd[k]) == len(vs[k]):
                d = vd[k] - vs[k]
                sc = float(np.abs(vs[k]).max()) or 1.0
                dv[k] = {"max_abs": float(np.abs(d).max()), "rms": float(np.sqrt(np.mean(d * d))), "max_rel_to_maxabs": float(np.abs(d).max() / sc),
                         "n_differs": int(np.count_nonzero(d))}
        out["conserved_diff_vs_src"] = dv
    out["note"] = "IC の影響は消去済みとは扱わない (plan §6 E1、諮問 ④)。限界として記録する"
    return out


def _load_problem_with_runs():
    """design_chain が凍結源 run (geometry.initial_line_run、問題ファイルからの相対) を探す場所を CASE_RUNS に向ける。"""
    from forge_design.evaluate import runner_axismach as RA
    if os.environ.get("CASE_RUNS"):
        _lp = RA.load_problem

        def _load(path, _lp=_lp):
            p = _lp(path)
            il = p.geometry.get("initial_line_run")
            if il and not Path(str(il)).is_absolute():
                p.geometry["initial_line_run"] = str((RUNS / str(il)).resolve())
            return p
        RA.load_problem = _load
    return RA


def prep(prep_dir: Path, arm: str, no_ic: bool = False) -> dict:
    """腕の入力を作る。no_ic=True はローカルの乾式確認用 (restart_field / interp_field は forge --resolve-species を起動するので
    飛ばし、メッシュ・壁の証拠・座標一致だけを確かめる。その prep は投入に使えないので prepare_info に印を残す)。"""
    if arm not in ARMS:
        raise SystemExit(f"arm は {list(ARMS)}")
    RA = _load_problem_with_runs()
    ic_run = RUNS / IC_RUN
    src = _last_res(ic_run)
    info = RA.prepare(C / ARMS[arm], prep_dir, nsteps=NSTEPS, ic_from=None, cfl_main=CFL_MAIN, implicit_relax=RELAX)
    if arm == "monoG1" and (info.get("wall_fit") or {}).get("mono_r2") != [0.0, 1.5]:
        raise RuntimeError(f"腕 B の壁に mono_r2 [0, 1.5] が入っていない: {(info.get('wall_fit') or {}).get('mono_r2')}")
    if arm == "pinG1" and (info.get("wall_fit") or {}).get("mono_r2") is not None:
        raise RuntimeError("腕 A の壁に mono_r2 が入っている")
    tool = IC_MODE[arm]
    with h5py.File(prep_dir / "nozzle.h5") as f, h5py.File(ic_run / "nozzle.h5") as g:
        same_mesh = f["/MESH/COORD"].shape == g["/MESH/COORD"].shape and np.array_equal(f["/MESH/COORD"][:], g["/MESH/COORD"][:])
    if tool == "restart_field" and not same_mesh:
        # 同一メッシュでなければ index コピーの restart は誤り (規則: 座標が違えば interp_field)。座標検査を緩めない
        raise RuntimeError(f"腕 A のメッシュが IC run {IC_RUN} と座標一致しない — restart_field を使えない")
    if no_ic:
        ev = wall_evidence(prep_dir)
        info.update(wall_arm=arm, DRY_NO_IC=True, wall_evidence_prep=ev, ic={"tool": tool, "skipped": True,
                    "mesh_coords_identical_to_ic_run": bool(same_mesh)})
        (prep_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
        print(arm, "(dry, IC なし) | M4", ev["status"], f"resid/tol {ev.get('resid_over_f32tol_max', float('nan')):.3g}",
              "| 座標一致 (IC run)", same_mesh)
        return info
    cmd = [sys.executable, str(TOOLS / f"{tool}.py"), str(src), str(prep_dir / "nozzle.h5")]
    if tool == "restart_field":
        cmd += ["--dst-run", str(prep_dir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    (prep_dir / f"{tool}.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"{tool} が失敗 (rc {r.returncode}):\n{(r.stdout + r.stderr)[-2000:]}")
    ic = {"run": str(ic_run), "res": src.name, "tool": tool, "log": f"{tool}.log", "mesh_coords_identical_to_ic_run": bool(same_mesh),
          "tool_last_line": ((r.stdout.strip().splitlines() or [""])[-1])}
    if tool == "interp_field":
        ic["map"] = ic_map_record(src, ic_run / "nozzle.h5", prep_dir / "nozzle.h5")
    # MOC 点群 (当てはめ前) が IC run と同じこと (腕の違いは当てはめだけ)
    a = np.loadtxt(prep_dir / "wall_design.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(ic_run / "wall_design.csv", delimiter=",", skiprows=1)
    ic["wall_design_vs_ic_run_max_abs"] = (float(np.abs(a - b).max()) if a.shape == b.shape else f"shape {a.shape} vs {b.shape}")
    ev = wall_evidence(prep_dir)
    if ev.get("status") != "consistent":
        raise RuntimeError(f"壁の証拠 (M4) が不一致: {json.dumps(ev, ensure_ascii=False)}")
    info.update(stages="soft", wall_arm=arm, ic=ic, wall_evidence_prep=ev,
                plan="plans/active/tooling-nozzle-throat-monotone-r2.md §6 E1")
    (prep_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    (prep_dir / "IC_MAP.json").write_text(json.dumps(ic, indent=1, default=str, ensure_ascii=False))
    print(arm, "mesh", json.dumps(info.get("mesh")), "| IC", tool, src.name, "| M4", ev["status"],
          f"resid/tol {ev['resid_over_f32tol_max']:.3g}", "|", (prep_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
    return info


def main(argv) -> int:
    if len(argv) in (3, 4) and argv[0] == "prep" and (len(argv) == 3 or argv[3] == "--no-ic"):
        prep(Path(argv[1]).resolve(), argv[2], no_ic=(len(argv) == 4))
        return 0
    if len(argv) == 2 and argv[0] == "run":
        from forge_design.evaluate import runner_axismach as RA
        rd = Path(argv[1]).resolve()
        if json.loads((rd / "prepare_info.json").read_text()).get("DRY_NO_IC"):
            raise SystemExit(f"{rd} は乾式確認 (--no-ic) の prep から作られている — IC が無いので回さない")
        rc = RA.run_staged(rd, cfl_main=CFL_MAIN, mid_stage=False, stages="soft")
        print(f"forge exit={rc}")
        return rc
    if len(argv) in (2, 3) and argv[0] == "evidence":
        print(json.dumps(wall_evidence(Path(argv[1]), Path(argv[2]) if len(argv) == 3 else None), indent=1, ensure_ascii=False))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
