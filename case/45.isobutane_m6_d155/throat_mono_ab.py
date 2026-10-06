"""Euler A/B (plan tooling-nozzle-throat-monotone-r2 §6 E1、§5.1 #5・#5b) の run 準備・実行と、壁の証拠 (M4)・IC 写像の記録。
腕 A (pinG1) = 現行壁 problem_d155_euler_pin_G1_recal.yaml / 腕 B (monoG1) = 単調壁 problem_d155_euler_pin_G1_recal_mono.yaml
(G1_recal + geometry.wall_fit_mono_r2 [0, 1.5])。生産 Euler 格子 G1、IC はどちらも run_0114 の最終場の**保存済み保存量**:
  A は同一メッシュなので restart_field (保存量の index コピー、ビット一致検査つき)、
  B は壁と近傍の節点が最大 0.5 µm 動くので、検証付き番号写像 ic_index_map.py --mode index (節点数・接続・論理位置 (i, j)・境界種別・
  座標系と単位・化学種とエネルギー基準・要素の反転なし・移動 ≤ 1 µm を検査し、保存量を番号で直接コピー; 諮問 2026-10-06 ①)。
  予備 A/B (§6 E1): `--with-nn <dir>` で、B の prep を IC を入れる前に複製し (同じ問題・同じ格子)、最近傍対応 (interp_field と同じ規則) で
  保存済み保存量を直接転送した prep を作る → run_0146 euler_icab_monoG1_nn。run_0143 (B の r1) が (β) を兼ねる。
起動は euler_grid_ab.py prep と同じ RA.prepare + RA.run_staged: soft (1 次 cfl 0.5、3000) → 本段 2 次 cfl 2・implicitRelax 0.7・
18000 step・1000 ごと出力。
usage: [CASE_RUNS=<run_0062 / run_0114 のある case dir>] python3 throat_mono_ab.py prep <prep_dir> {pinG1|monoG1} [--with-nn <prep_dir_nn>] [--dry]
       python3 throat_mono_ab.py run <run_dir>               (準備済み run を soft → 本段で回す; forge は run_case.sh 経由)
       python3 throat_mono_ab.py verify-prep <prep_dir> [<run_dir> ...]   (prep の入力が IC 写像の直後のままで、run がその prep から作られたか)
       python3 throat_mono_ab.py evidence <run_dir> [<run_dir_other_arm>]   (壁の証拠だけ出す; JSON を標準出力へ)
--dry (旧名 --no-ic): ローカルの乾式確認 (forge を起動しない)。A は restart_field (forge --resolve-species を起動する) を飛ばし、
  B・最近傍は ic_index_map を --no-species-resolve で通す (検査・転送・記録は本番と同じで、化学種の属性だけ付けない)。
  その prep は prepare_info に DRY の印が付き、run・verify-prep が拒否する。
壁の証拠 (M4): prepare_info.json の wall_fit.spline (当てはめ後のノット・係数・次数) を、変換後メッシュ (nozzle.h5, float32) の
壁節点の x で評価し、節点の r と照合する (許容 = float32 の丸め ulp(r) + |r′|·ulp(x))。wall_design.csv は当てはめ前の MOC 点群なので
証拠にならない。r″ の山の有無は壁節点からは認定しない (解析形で形状ゲート S1 が担う)。
"""
import argparse
import json
import os
import shutil
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
MONO_R2 = {"pinG1": None, "monoG1": [0.0, 1.5]}       # prepare_info の wall_fit.mono_r2 の期待値 (完全一致)
# IC の入れ方 (variant → 道具): 腕 A は restart_field、腕 B は番号写像、予備 A/B の (α) は最近傍 (どちらも保存量の直接転送)
IC_MODE = {"pinG1": "restart_field", "monoG1": "index", "monoG1_nn": "nearest"}
NSTEPS, CFL_MAIN, RELAX = 18000, 2.0, 0.7


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




def _run_tool(cmd, log: Path) -> subprocess.CompletedProcess:
    r = subprocess.run([str(x) for x in cmd], capture_output=True, text=True)
    log.write_text(r.stdout + r.stderr)
    return r


def _map_ic(prep_dir: Path, variant: str, info: dict, src: Path, ic_run: Path, dry: bool, same_mesh: bool, extra=None) -> dict:
    """prep_dir の nozzle.h5 に IC を入れ、壁の証拠を照合して prepare_info.json・IC_MAP.json を書く。
    variant: pinG1 (restart_field) / monoG1 (ic_index_map --mode index) / monoG1_nn (ic_index_map --mode nearest)。"""
    from ic_index_map import _sha_file, mesh_digest
    mode = IC_MODE[variant]
    dst = prep_dir / "nozzle.h5"
    ic = {"run": str(ic_run), "res": src.name, "mesh_coords_identical_to_ic_run": bool(same_mesh), **(extra or {})}
    if mode == "restart_field":
        ic.update(tool="restart_field", log="restart_field.log")
        if dry:
            ic.update(skipped="restart_field は forge --resolve-species を起動するので乾式確認では飛ばす", VERDICT=None)
        else:
            r = _run_tool([sys.executable, TOOLS / "restart_field.py", src, dst, "--dst-run", prep_dir], prep_dir / "restart_field.log")
            last = (r.stdout.strip().splitlines() or [""])[-1]
            if r.returncode != 0 or not (last.startswith("VERDICT: OK") and "SRC とビット一致" in last):
                raise RuntimeError(f"restart_field が失敗・ビット一致でない (rc {r.returncode}):\n{(r.stdout + r.stderr)[-2000:]}")
            ic.update(tool_last_line=last, VERDICT="OK")
        ic.update(dst_sha256_after=_sha_file(dst), dst_mesh_digest=mesh_digest(dst))
    else:
        cmd = [sys.executable, C / "ic_index_map.py", src, dst, "--mode", mode, "--src-mesh", ic_run / "nozzle.h5",
               "--src-run", ic_run, "--dst-run", prep_dir, "--record", prep_dir / "IC_MAP.json"]
        if dry:
            cmd.append("--no-species-resolve")
        r = _run_tool(cmd, prep_dir / "ic_index_map.log")
        rec = json.loads((prep_dir / "IC_MAP.json").read_text()) if (prep_dir / "IC_MAP.json").exists() else {}
        if r.returncode != 0 or rec.get("VERDICT") != "OK":
            raise RuntimeError(f"ic_index_map ({mode}) が失敗 (rc {r.returncode}, VERDICT {rec.get('VERDICT')}):\n{(r.stdout + r.stderr)[-3000:]}")
        dr = rec["checks"]["displacement"]["detail"]
        nv = rec["nearest_vs_index"]
        ic.update(tool="ic_index_map", mode=mode, log="ic_index_map.log", record="IC_MAP.json", VERDICT=rec["VERDICT"],
                  species_resolved=rec["species_resolved"], transferred=rec["transferred"],
                  displacement_max_um=dr["max_um"], n_moved=dr["n_moved"], nearest_mismatch=nv["n_mismatch"],
                  nearest_mismatch_dj=nv["dj_counts"], dst_sha256_after=rec["dst_sha256_after"], dst_mesh_digest=rec["dst_mesh_digest"],
                  tool_last_line=(r.stdout.strip().splitlines() or [""])[-1])
    # MOC 点群 (当てはめ前) が IC run と同じこと (腕の違いは当てはめだけ)
    a = np.loadtxt(prep_dir / "wall_design.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(ic_run / "wall_design.csv", delimiter=",", skiprows=1)
    ic["wall_design_vs_ic_run_max_abs"] = (float(np.abs(a - b).max()) if a.shape == b.shape else f"shape {a.shape} vs {b.shape}")
    ev = wall_evidence(prep_dir)
    if ev.get("status") != "consistent":
        raise RuntimeError(f"壁の証拠 (M4) が不一致: {json.dumps(ev, ensure_ascii=False)}")
    info.update(stages="soft", wall_arm=variant, ic=ic, wall_evidence_prep=ev, plan="plans/active/tooling-nozzle-throat-monotone-r2.md §6 E1")
    info.pop("DRY_NO_IC", None)
    if dry:
        info["DRY"] = True
    else:
        info.pop("DRY", None)
    (prep_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    if mode == "restart_field":                      # B・最近傍の IC_MAP.json は ic_index_map の記録そのもの (上書きしない)
        (prep_dir / "IC_MAP.json").write_text(json.dumps(ic, indent=1, default=str, ensure_ascii=False))
    print(variant, "mesh", json.dumps(info.get("mesh")), "| IC", ic["tool"], ic.get("mode", ""), src.name, "VERDICT", ic.get("VERDICT"),
          ("(dry)" if dry else ""), "| M4", ev["status"], f"resid/tol {ev['resid_over_f32tol_max']:.3g}", "|",
          (prep_dir / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
    if mode != "restart_field":
        print(f"  移動 最大 {ic['displacement_max_um']:.4f} µm / 動いた節点 {ic['n_moved']} / 最近傍が番号と食い違う節点 {ic['nearest_mismatch']} "
              f"(j の差 {ic['nearest_mismatch_dj']})")
    return info


def prep(prep_dir: Path, arm: str, dry: bool = False, nn_dir: Path | None = None) -> dict:
    """腕の入力を作る。nn_dir (monoG1 のみ): 予備 A/B の最近傍の prep を、IC を入れる前の B の prep の複製から作る (同じ格子)。
    dry=True はローカルの乾式確認用 (forge を起動しない; その prep は投入に使えない)。"""
    if arm not in ARMS:
        raise SystemExit(f"arm は {list(ARMS)}")
    if nn_dir is not None and arm != "monoG1":
        raise SystemExit("--with-nn は monoG1 (腕 B) の prep にだけ付ける (予備 A/B は B の格子で行う)")
    if nn_dir is not None and Path(nn_dir).exists():
        raise SystemExit(f"{nn_dir} が既にある (消してから作る)")
    RA = _load_problem_with_runs()
    ic_run = RUNS / IC_RUN
    src = _last_res(ic_run)
    info = RA.prepare(C / ARMS[arm], prep_dir, nsteps=NSTEPS, ic_from=None, cfl_main=CFL_MAIN, implicit_relax=RELAX)
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf:
        raise RuntimeError("prepare_info の wall_fit に mono_r2 の記録が無い (壁の当てはめ設定の証拠が欠損)")
    from throat_mono_judge import mono_r2_matches
    if not mono_r2_matches(wf["mono_r2"], MONO_R2[arm]):
        raise RuntimeError(f"腕 {arm} の壁の mono_r2 {wf['mono_r2']!r} が期待 {MONO_R2[arm]!r} と完全一致しない")
    with h5py.File(prep_dir / "nozzle.h5") as f, h5py.File(ic_run / "nozzle.h5") as g:
        same_mesh = f["/MESH/COORD"].shape == g["/MESH/COORD"].shape and np.array_equal(f["/MESH/COORD"][:], g["/MESH/COORD"][:])
    if IC_MODE[arm] == "restart_field" and not same_mesh:
        # 同一メッシュでなければ index コピーの restart は誤り。restart_field の座標検査を緩めない
        raise RuntimeError(f"腕 A のメッシュが IC run {IC_RUN} と座標一致しない — restart_field を使えない")
    if nn_dir is not None:
        nn_dir = Path(nn_dir).resolve()
        shutil.copytree(prep_dir, nn_dir)            # IC を入れる前の B の prep (同じ問題・同じ格子・同じ設定) を複製
    out = _map_ic(prep_dir, arm, info, src, ic_run, dry, same_mesh)
    if nn_dir is not None:
        info_nn = json.loads((nn_dir / "prepare_info.json").read_text())
        _map_ic(nn_dir, "monoG1_nn", info_nn, src, ic_run, dry, same_mesh,
                extra={"grid_copied_from": str(prep_dir), "icab": "plan §6 E1 予備 A/B の (α)、(β) = 腕 B の r1"})
        a, b = (json.loads((d / "IC_MAP.json").read_text())["dst_mesh_digest"] for d in (prep_dir, nn_dir))
        if a != b:
            raise RuntimeError(f"予備 A/B の 2 つの prep の格子が同じでない ({a[:16]} / {b[:16]})")
    return out


def verify_prep(prep_dir: Path, run_dirs=()) -> list:
    """prep の入力が IC 写像の直後のまま (nozzle.h5 の sha256 = IC_MAP.json の記録) で、run_dirs がその prep から作られたか
    (run の IC_MAP.json の記録が同じ)。腕 B の r2・r3 を r1 と同じ入力から作るときに使う。戻り値 = 不成立の理由 (空なら成立)。"""
    from ic_index_map import _sha_file
    bad = []
    try:
        info = json.loads((prep_dir / "prepare_info.json").read_text())
        rec = json.loads((prep_dir / "IC_MAP.json").read_text())
    except (OSError, ValueError) as e:
        return [f"{prep_dir}: prepare_info.json / IC_MAP.json を読めない ({e})"]
    if info.get("DRY") or info.get("DRY_NO_IC"):
        bad.append(f"{prep_dir} は乾式確認 (DRY) の prep")
    if rec.get("VERDICT") != "OK":
        bad.append(f"{prep_dir}: IC 写像の VERDICT {rec.get('VERDICT')!r}")
    sha = rec.get("dst_sha256_after")
    if not sha or _sha_file(prep_dir / "nozzle.h5") != sha:
        bad.append(f"{prep_dir}/nozzle.h5 が IC 写像の記録 (sha256 {str(sha)[:16]}) と違う (写像の後に変わった)")
    for rd in run_dirs:
        try:
            r2 = json.loads((Path(rd) / "IC_MAP.json").read_text())
            i2 = json.loads((Path(rd) / "prepare_info.json").read_text())
        except (OSError, ValueError) as e:
            bad.append(f"{rd}: IC_MAP.json / prepare_info.json を読めない ({e})")
            continue
        if r2.get("dst_sha256_after") != sha:
            bad.append(f"{rd}: IC 写像の記録が prep と違う (別の prep から作られた)")
        if i2.get("wall_arm") != info.get("wall_arm"):
            bad.append(f"{rd}: wall_arm {i2.get('wall_arm')!r} が prep {info.get('wall_arm')!r} と違う")
    return bad


def main(argv) -> int:
    if argv and argv[0] == "prep":
        ap = argparse.ArgumentParser(prog="throat_mono_ab.py prep")
        ap.add_argument("prep_dir")
        ap.add_argument("arm", choices=list(ARMS))
        ap.add_argument("--with-nn", help="予備 A/B の最近傍の prep (monoG1 のみ)")
        ap.add_argument("--dry", "--no-ic", dest="dry", action="store_true", help="乾式確認 (forge を起動しない)")
        a = ap.parse_args(argv[1:])
        prep(Path(a.prep_dir).resolve(), a.arm, dry=a.dry, nn_dir=(Path(a.with_nn) if a.with_nn else None))
        return 0
    if len(argv) == 2 and argv[0] == "run":
        from forge_design.evaluate import runner_axismach as RA
        rd = Path(argv[1]).resolve()
        info = json.loads((rd / "prepare_info.json").read_text())
        if info.get("DRY") or info.get("DRY_NO_IC"):
            raise SystemExit(f"{rd} は乾式確認 (--dry) の prep から作られている — 化学種の属性・IC がそろっていないので回さない")
        if (info.get("ic") or {}).get("VERDICT") != "OK":
            raise SystemExit(f"{rd}: IC 写像の VERDICT が OK でない ({(info.get('ic') or {}).get('VERDICT')!r}) — 回さない")
        rc = RA.run_staged(rd, cfl_main=CFL_MAIN, mid_stage=False, stages="soft")
        print(f"forge exit={rc}")
        return rc
    if len(argv) >= 2 and argv[0] == "verify-prep":
        bad = verify_prep(Path(argv[1]).resolve(), [Path(x).resolve() for x in argv[2:]])
        print("VERIFY-PREP: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if len(argv) in (2, 3) and argv[0] == "evidence":
        print(json.dumps(wall_evidence(Path(argv[1]), Path(argv[2]) if len(argv) == 3 else None), indent=1, ensure_ascii=False))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
