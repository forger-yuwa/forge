"""plan verification-case45-euler-total-enthalpy §6 E2 (2026-10-07 登録、改訂 a2259768・f158bac2・再改訂 a59e61b5) の run 準備・実行。
半径方向の配点プロファイルへの感度: 変えるのは mesh.wall_first_frac と mesh.wall_first_frac_throat の 2 鍵だけ (連動して 1 要因)。
  A = problem_d155_euler_t0cluster_g1.yaml    (G1 のまま 1.3e-5・4.5e-6)  → run_0161_euler_t0cluster_g1
  B = problem_d155_euler_t0cluster_u5em3.yaml (両方 0.005)                 → run_0162_euler_t0cluster_u5em3
  どちらも E1 の run_0153 の問題 problem_d155_euler_pin_G1_recal_mono_moc.yaml の写し (name を除いて A は同一、B は上の 2 鍵だけが違う)。
準備 (prep): 問題の検査 (check-problems) → RA.prepare (等エントロピー IC。run_0153 と同じ生成手順) を A・B で → 各腕の検査
  (mono_r2・MOC のキーとゲート・壁の証拠 M4・メッシュ品質 PASS・解決済みの Mesh2DParams が登録の値) → A と B の照合
  (列の x と壁の r がビット単位で同一、BC・solverConfig・設計壁・目標軸分布が同一、prepare_info が mesh 欄を除いて同一、
  Mesh2DParams が 2 鍵を除いて同一) → A の格子が run_0153 の格子と同一 (座標・設計壁; 本番は必須) →
  起動前の IC の検査 (全節点で |T0 − 1600| ≤ 1 K、復元は E1 と同じ euler_t0_stage_ab.recon) → 実効のメッシュの記録 (E2_MESH.json・
  E2_MESH_sections.csv) → prepare_info.json に段・腕・IC の記録 → E2_PREP.json (入力の sha256・照合の結果)。
  不成立はすべて例外で止める (投入スクリプトは forge を 1 本も起動しない)。
実行 (run): runner_axismach.run_staged(stages="soft", mid_stage=False) — soft (1 次・cfl 0.5・3000 step) → 本段 (2 次・cfl 2・
  implicitRelax 0.7・54000 step・1000 ごと出力; 再改訂 a59e61b5)。本段の step 数・出力間隔は評価器の定数 (euler_t0_e2_eval.MAIN_NSTEPS・
  OUT_INTERVAL) を読み、prep の後と起動の直前に solverConfig.yaml の実効値と評価器の窓を照合する (食い違えば止める)。
  soft 段の出力の退避: runner は段の終わりに _restart_same_mesh を呼んでから
  `glob("res_*")` を消す (runner_axismach.py:885–887) ので、_restart_same_mesh を包み、呼ばれる前に run 直下の res_* (と段の
  forge_run.log・RUN_PROVENANCE.txt・run_case_stdout.log) を `_soft_stage/` に写す (サブ dir なので glob に当たらない)。
  写しの sha256 を元と照合し、`_soft_stage/SOFT_STAGE_SAVED.json` に記録する。design/ は変えない (import だけ)。
usage: [CASE_RUNS=<run_0062・run_0114 [・run_0153] のある case dir>] python3 euler_t0_e2.py check-problems
       python3 euler_t0_e2.py prep [--out-root DIR] [--dry]
       python3 euler_t0_e2.py verify-prep <prep_dir> <run_dir>
       python3 euler_t0_e2.py run <run_dir>
--dry: ローカルの乾式確認 (forge を起動しない)。FORGE_BIN を存在しない道に向け、FORGE_ALLOW_UNVERIFIED_SPECIES=1 にする (化学種の属性は
  付かない)。起動前の IC の検査の熱物性は CASE_RUNS の run_0114 の解決済み記録で代用する (記録に明記)。run_0153 が無ければ照合を
  「未確認 (DRY)」として続ける。その prep は prepare_info に DRY の印が付き、run・verify-prep が拒否する。
  記録は _band_ab/euler_t0_e2_mesh_dry.json (本番は _band_ab/euler_t0_e2_mesh.json)。
"""
from __future__ import annotations

import argparse
import copy
import csv
import dataclasses
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import throat_mono_ab as TM  # noqa: E402  (RUNS・問題の読み込み・壁の証拠・段の設定 CFL_MAIN/RELAX を共有する。step 数は EV から)
import euler_t0_e2_eval as EV  # noqa: E402  (領域の定義を評価器と共有する)

ROOT = C.parents[1]
PLAN = EV.PLAN
PLAN_REG_COMMIT = EV.PLAN_REG_COMMIT
BASE_PROBLEM = "problem_d155_euler_pin_G1_recal_mono_moc.yaml"
ARMS = {"A": {"problem": "problem_d155_euler_t0cluster_g1.yaml", "run": EV.RUNS["A"], "prep": "_prep_e2_g1",
              "wall_first_frac": 1.3e-5, "wall_first_frac_throat": 4.5e-6},
        "B": {"problem": "problem_d155_euler_t0cluster_u5em3.yaml", "run": EV.RUNS["B"], "prep": "_prep_e2_u5em3",
              "wall_first_frac": 0.005, "wall_first_frac_throat": 0.005}}
CHANGED = (("mesh", "wall_first_frac"), ("mesh", "wall_first_frac_throat"))
REF_RUN_WALL = "run_0153_euler_icdep_mocG1_isen"     # E1 の run (A はこれと同じ格子のはず)
REF_PREP_ISEN = "_prep_moc_v5_mocG1_isen"            # run_0153 の起動前の prep (等エントロピー IC; 記録用の照合)
DRY_THERMO_RUN = "run_0114_euler_pin_G1_recal_ext6k"  # 乾式確認で IC の検査に使う熱物性の記録 (CASE_RUNS)
MONO_R2 = [0.0, 1.5]
MOC_EXPECT = {"axis_limit": "analytic", "corrector": "converge"}
IC_TOL_K = 1.0
NI_NJ = (2000, 97)
SAME_FILES = ("bcondConfig.yaml", "solverConfig.yaml", "probe.yaml", "wall_design.csv", "target_axis_M.csv")
MESH_JSON, MESH_CSV, PREP_JSON = "E2_MESH.json", "E2_MESH_sections.csv", "E2_PREP.json"
SOFT_SAVE_EXTRA = ("forge_run.log", "RUN_PROVENANCE.txt", "run_case_stdout.log", "forge_launches.jsonl")
DESIGN_FILES = ("design/forge_design/evaluate/runner_axismach.py", "design/forge_design/evaluate/runner.py",
                "design/forge_design/meshing/mesh2d.py", "design/forge_design/evaluate/ic.py")
# 記録の帯 [r_t] (登録の 3 帯に加え、G1 のブレンドの区切り −9・−4・1・17 で分けた帯)
SUMMARY_BANDS_RT = (("x<-9", -np.inf, -9.0, "()"), ("-9<=x<-4", -9.0, -4.0, "[)"), ("-4<=x<=1", -4.0, 1.0, "[]"),
                    ("1<x<17", 1.0, 17.0, "()"), ("x>=17", 17.0, np.inf, "[)"))


def _sha(p) -> str | None:
    return EV._sha(p)


def _sha_bytes(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()


def _git(*args) -> str:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, timeout=30)
        return r.stdout.strip() if r.returncode == 0 else f"(git rc {r.returncode})"
    except (OSError, subprocess.SubprocessError) as e:
        return f"(git 不可: {type(e).__name__})"


# --- 問題の検査 ----------------------------------------------------------------------------------------------------------------
def _diff_paths(a, b, path=(), ignore=frozenset()) -> list:
    """2 つの YAML の木の違う位置 (ignore の位置は飛ばす)。型の違いも違いとする。"""
    if path in ignore:
        return []
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            if (path + (k,)) in ignore:
                continue
            if k not in a or k not in b:
                out.append(".".join(map(str, path + (k,))) + (" (A に無い)" if k not in a else " (B に無い)"))
            else:
                out += _diff_paths(a[k], b[k], path + (k,), ignore)
        return out
    if type(a) is not type(b) or a != b:
        return [f"{'.'.join(map(str, path))}: {a!r} / {b!r}"]
    return []


def check_problems(case: Path = C) -> list:
    """登録どおりか: A = 元の問題 (name を除いて同一)、B = A (name と 2 鍵を除いて同一)、2 鍵の値が登録の値 (float)。不成立の理由のリスト。"""
    load = lambda n: yaml.safe_load((case / n).read_text())  # noqa: E731
    base, a, b = load(BASE_PROBLEM), load(ARMS["A"]["problem"]), load(ARMS["B"]["problem"])
    bad = [f"A ({ARMS['A']['problem']}) と元の問題 ({BASE_PROBLEM}) の違い: {d}" for d in _diff_paths(base, a, ignore={("name",)})]
    bad += [f"B と A の違い (2 鍵と name 以外): {d}" for d in _diff_paths(a, b, ignore={("name",), *CHANGED})]
    for arm, doc in (("A", a), ("B", b)):
        for k in ("wall_first_frac", "wall_first_frac_throat"):
            v = (doc.get("mesh") or {}).get(k)
            if not isinstance(v, float) or v != ARMS[arm][k]:
                bad.append(f"{arm} の mesh.{k} = {v!r} (登録 {ARMS[arm][k]!r}、float であること)")
    names = {base.get("name"), a.get("name"), b.get("name")}
    if len(names) != 3:
        bad.append(f"name が 3 つの問題で別々でない ({base.get('name')}, {a.get('name')}, {b.get('name')})")
    return bad


# --- 実効のメッシュの記録 ----------------------------------------------------------------------------------------------------------
def _grid(d: Path):
    """nozzle.h5 → X, R [m] (ni, nj) float64 (float32 の値)、座標の dtype。並びは VIZMESH/CONNE で確かめる。"""
    from ic_index_map import structured_shape
    with h5py.File(d / "nozzle.h5", "r") as f:
        raw = f["/MESH/COORD"][:]
        nc = raw.reshape(-1, 3).astype(np.float64)
        ni, nj = structured_shape(f["VIZMESH/CONNE"], nc.shape[0])
    return nc[:, 0].reshape(ni, nj), nc[:, 1].reshape(ni, nj), nc[:, 2].reshape(ni, nj), str(raw.dtype)


def _topology_digest(h5path: Path) -> str:
    """接続 (座標以外の格子の位相) のハッシュ: ic_index_map の TOPOLOGY と境界の位相のデータセット。"""
    from ic_index_map import BCOND_TOPO, TOPOLOGY
    h = hashlib.sha256()
    with h5py.File(h5path, "r") as f:
        names = [t for t in TOPOLOGY if t in f]
        names += sorted(f"BCONDS/{b}/{k}" for b in f["BCONDS"] for k in BCOND_TOPO if k in f["BCONDS"][b])
        for n in names:
            a = np.asarray(f[n])
            h.update(n.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()


def mesh_record(d: Path, problem: Path, mp, S: float) -> tuple:
    """登録「記録する実効のメッシュ」: 解決済みの Mesh2DParams と、変換後の nozzle.h5 で測った断面ごとの量・領域内の節点数・ハッシュ。
    戻り (記録 dict, 断面の行のリスト)。"""
    from ic_index_map import mesh_digest
    X, R, Z, dtype = _grid(d)
    ni, nj = X.shape
    with h5py.File(d / "nozzle.h5", "r") as f:
        coord_sha = _sha_bytes(f["/MESH/COORD"][:])
    g = np.diff(R, axis=1)                       # 軸側 (j = 0) → 壁側
    rw = R[:, -1]
    ratio = np.maximum(g[:, 1:] / g[:, :-1], g[:, :-1] / g[:, 1:]).max(axis=1)
    eta = R / rw[:, None]
    eb = EV.eta_band(eta)
    rows = []
    for i in range(ni):
        rows.append({"i": i, "x_m": X[i, 0], "x_rt": X[i, 0] / S, "r_w_m": rw[i], "r_w_rt": rw[i] / S,
                     "wall_gap_m": g[i, -1], "wall_gap_rw": g[i, -1] / rw[i], "axis_gap_m": g[i, 0], "axis_gap_rw": g[i, 0] / rw[i],
                     "max_adjacent_ratio": ratio[i],
                     "n_axis": int(np.sum(eb[i] == 0)), "n_core": int(np.sum(eb[i] == 1)), "n_wall": int(np.sum(eb[i] == 2))})
    x_rt = X[:, 0] / S
    rid = EV.region_id(X / S, eta)

    def band(lo, hi, closed):
        """closed: 端の含み方 "[)"・"[]"・"()"・"(]"。"""
        m = ((x_rt >= lo) if closed[0] == "[" else (x_rt > lo)) & ((x_rt <= hi) if closed[1] == "]" else (x_rt < hi))
        if not m.any():
            return {"n_columns": 0}
        out = {"n_columns": int(m.sum()), "x_rt_range": [float(x_rt[m].min()), float(x_rt[m].max())]}
        for k, v in (("wall_gap_rw", g[m, -1] / rw[m]), ("axis_gap_rw", g[m, 0] / rw[m]), ("wall_gap_m", g[m, -1]),
                     ("axis_gap_m", g[m, 0]), ("max_adjacent_ratio", ratio[m]), ("r_w_rt", rw[m] / S)):
            out[k] = [float(v.min()), float(v.max())]
        return out
    rec = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "problem": problem.name, "problem_sha256": _sha(problem),
           "mesh2d_params": dataclasses.asdict(mp),
           "axis_cap": {"axis_cap_frac": mp.axis_cap_frac, "axis_gap_frac": mp.axis_gap_frac,
                        "capped": mp.axis_cap_frac is not None or mp.axis_gap_frac is not None},
           "ni": ni, "nj": nj, "scale_m": S, "coord_dtype": dtype, "z_all_zero": bool(np.all(Z == 0.0)),
           "summary_bands_rt": {name: band(lo, hi, cl) for name, lo, hi, cl in SUMMARY_BANDS_RT},
           "registration_bands_rt": {xb: band(lo, hi, cl) for xb, lo, hi, cl in
                                     (("up", -np.inf, EV.X_EDGES_RT[0], "()"), ("thr", EV.X_EDGES_RT[0], EV.X_EDGES_RT[1], "[]"),
                                      ("dn", EV.X_EDGES_RT[1], np.inf, "()"))},
           "nodes_by_region": {k: int(np.sum(rid == i)) for i, k in enumerate(EV.REGION_KEYS)},
           "hashes": {"coord_sha256": coord_sha, "topology_sha256": _topology_digest(d / "nozzle.h5"),
                      "mesh_digest_ic_index_map": mesh_digest(d / "nozzle.h5"), "sections_csv": None},
           "note": "prepare_info.json の mesh 欄 (ni・nj・wall_first_frac) はこれらを記録しないので別に保存する (登録)"}
    return rec, rows


def write_sections(path: Path, rows: list) -> None:
    keys = list(rows[0])
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([(f"{r[k]:.10g}" if isinstance(r[k], float) else r[k]) for k in keys])


# --- 起動前の IC の検査 --------------------------------------------------------------------------------------------------------------
def ic_check(d: Path, thermo_dir: Path, dry: bool) -> dict:
    """全節点で |T0 − 1600| ≤ 1 K (登録)。復元は E1 と同じ (euler_t0_stage_ab.recon)。thermo_dir はその解決済み熱物性を読む run
    (本番は prep 自身 = nozzle.h5 の化学種の属性、乾式は run_0114 の記録で代用)。"""
    import euler_t0_stage_ab as E1
    rec = {"file": f"{d.name}/nozzle.h5", "sha256": _sha(d / "nozzle.h5"), "thermo_dir": str(thermo_dir), "tol_K": IC_TOL_K,
           "Tt_reg": EV.TT_REG, "dry_thermo_substitute": bool(dry)}
    Tt = E1.V.read_Tt(d)
    rec["Tt_bcond"] = Tt
    r = E1.recon(thermo_dir, d / "nozzle.h5")
    T0 = np.asarray(r["T0"], dtype=np.float64)
    dev = T0 - EV.TT_REG
    X, R, _, _ = _grid(d)
    S = float(json.loads((d / "prepare_info.json").read_text())["scale_m"])
    k = int(np.nanargmax(np.abs(dev))) if np.any(np.isfinite(dev)) else 0
    i, j = divmod(k, X.shape[1])
    rec.update(thermo_source=r["src"], newton_ok=bool(r["newton_ok"]), n_nodes=int(T0.size), n_nonfinite=int(np.sum(~np.isfinite(T0))),
               T0_min=float(np.nanmin(T0)), T0_max=float(np.nanmax(T0)), max_abs_dev_K=float(np.nanmax(np.abs(dev))),
               argmax_abs={"i": i, "j": j, "j_from_wall": X.shape[1] - 1 - j, "x_rt": float(X[i, j] / S), "eta": float(R[i, j] / R[i, -1])},
               n_abs_dev_gt_tol=int(np.sum(np.abs(dev) > IC_TOL_K)))
    why = []
    if Tt != EV.TT_REG:
        why.append(f"入口の Tt {Tt!r} が {EV.TT_REG} でない")
    if not rec["newton_ok"]:
        why.append("復元の Newton が収束しない節点がある")
    if rec["n_nonfinite"]:
        why.append(f"非有限の全温 {rec['n_nonfinite']} 節点")
    if not rec["max_abs_dev_K"] <= IC_TOL_K:
        why.append(f"max|T0 − {EV.TT_REG:g}| = {rec['max_abs_dev_K']:.4g} K > {IC_TOL_K} K ({rec['n_abs_dev_gt_tol']} 節点)")
    rec["problems"] = why
    rec["VERDICT"] = ("OK" if not why else "FAIL") if not dry else ("DRY (熱物性は run_0114 の記録で代用)" if not why else "FAIL (DRY)")
    return rec


# --- A と B の照合・run_0153 との照合 ------------------------------------------------------------------------------------------------
def cross_check(dA: Path, dB: Path, mpA, mpB, infoA: dict, infoB: dict) -> dict:
    """A と B の照合 (登録「A と B の壁の座標と軸方向の節点が同一」と、変える要因以外の同一性)。不成立は failures に入れる。"""
    fails, rec = [], {}
    XA, RA_, ZA, dta = _grid(dA)
    XB, RB_, ZB, dtb = _grid(dB)
    rec["shape"] = {"A": list(XA.shape), "B": list(XB.shape)}
    if XA.shape != XB.shape or XA.shape != NI_NJ:
        fails.append(f"格子の形 A {XA.shape} / B {XB.shape} (登録 {NI_NJ})")
    else:
        rec["x_identical"] = bool(np.array_equal(XA, XB))
        rec["x_column_constant"] = bool(np.all(XA == XA[:, :1]) and np.all(XB == XB[:, :1]))
        rec["wall_r_identical"] = bool(np.array_equal(RA_[:, -1], RB_[:, -1]))
        rec["axis_r_zero"] = bool(np.all(RA_[:, 0] == 0.0) and np.all(RB_[:, 0] == 0.0))
        rec["z_identical_zero"] = bool(np.all(ZA == 0.0) and np.all(ZB == 0.0))
        rec["interior_r_max_abs_diff_m"] = float(np.abs(RA_ - RB_).max())
        for k in ("x_identical", "x_column_constant", "wall_r_identical", "axis_r_zero", "z_identical_zero"):
            if not rec[k]:
                fails.append(f"{k} が成り立たない")
        if dta != dtb:
            fails.append(f"座標の dtype が違う ({dta} / {dtb})")
    rec["topology_identical"] = _topology_digest(dA / "nozzle.h5") == _topology_digest(dB / "nozzle.h5")   # 記録 (同じ ni・nj なら同じはず)
    rec["same_files"] = {}
    for n in SAME_FILES:
        sa, sb = _sha(dA / n), _sha(dB / n)
        rec["same_files"][n] = {"A": sa, "B": sb, "identical": bool(sa and sa == sb)}
        if not (sa and sa == sb):
            fails.append(f"{n} が A と B で同一でない")
    ia, ib = copy.deepcopy(infoA), copy.deepcopy(infoB)
    rec["prepare_info_mesh"] = {"A": ia.pop("mesh", None), "B": ib.pop("mesh", None)}
    dpi = _diff_paths(ia, ib)
    rec["prepare_info_diff_except_mesh"] = dpi
    if dpi:
        fails.append(f"prepare_info.json が mesh 欄の外で違う: {dpi[:5]}")
    pa, pb = dataclasses.asdict(mpA), dataclasses.asdict(mpB)
    dmp = sorted(k for k in pa if pa[k] != pb[k])
    rec["mesh2d_params_diff_keys"] = dmp
    if dmp != sorted(k for _, k in CHANGED):
        fails.append(f"解決済みの Mesh2DParams の違いが 2 鍵だけでない: {dmp}")
    rec["failures"] = fails
    rec["ok"] = not fails
    return rec


def ref_check(dA: Path, dry: bool) -> dict:
    """A の格子が E1 の run_0153 と同一 (登録「E1 の run_0153 と同じ壁・2000 × 97・軸方向の節点」)。本番は run_0153 が無ければ不成立。
    run_0153 の起動前の prep (等エントロピー IC) と A の IC の保存量の一致は記録だけ (同じ生成手順の確認)。"""
    ref = TM.RUNS / REF_RUN_WALL
    rec = {"ref_run": str(ref), "failures": []}
    if not (ref / "nozzle.h5").is_file():
        rec["status"] = "未確認 (DRY: run_0153 が無い)" if dry else "run_0153 が無い"
        if not dry:
            rec["failures"].append(f"{ref}/nozzle.h5 が無い (同じ壁・格子の照合ができない)")
        rec["ok"] = not rec["failures"]
        return rec
    with h5py.File(dA / "nozzle.h5", "r") as fa, h5py.File(ref / "nozzle.h5", "r") as fr:
        ca, cr = fa["/MESH/COORD"][:], fr["/MESH/COORD"][:]
        rec["coord_identical"] = bool(ca.shape == cr.shape and ca.dtype == cr.dtype and np.array_equal(ca, cr))
        rec["vizmesh_conne_identical"] = bool(np.array_equal(fa["VIZMESH/CONNE"][:], fr["VIZMESH/CONNE"][:]))
    rec["wall_design_identical"] = _sha(dA / "wall_design.csv") == _sha(ref / "wall_design.csv")
    rec["topology_identical"] = _topology_digest(dA / "nozzle.h5") == _topology_digest(ref / "nozzle.h5")   # 記録
    for k in ("coord_identical", "vizmesh_conne_identical", "wall_design_identical"):
        if not rec[k]:
            rec["failures"].append(f"A と {REF_RUN_WALL} で {k} が成り立たない")
    pi = TM.RUNS / REF_PREP_ISEN / "nozzle.h5"
    if pi.is_file():
        with h5py.File(dA / "nozzle.h5", "r") as fa, h5py.File(pi, "r") as fp:
            names = sorted(k for k in fa["VALUE"] if k in fp["VALUE"])
            rec["ic_vs_run0153_prep"] = {k: bool(np.array_equal(fa[f"VALUE/{k}"][:], fp[f"VALUE/{k}"][:])) for k in names}
    else:
        rec["ic_vs_run0153_prep"] = f"{REF_PREP_ISEN} が無い (記録なし)"
    rec["status"] = "OK" if not rec["failures"] else "FAIL"
    rec["ok"] = not rec["failures"]
    return rec


def _arm_checks(d: Path, info: dict, mp, arm: str) -> tuple:
    bad = []
    from throat_mono_judge import mono_r2_matches
    wf = info.get("wall_fit") or {}
    if "mono_r2" not in wf or not mono_r2_matches(wf["mono_r2"], MONO_R2):
        bad.append(f"{arm}: 壁の mono_r2 {wf.get('mono_r2')!r} が {MONO_R2} と完全一致しない")
    moc = info.get("moc") or {}
    gate = moc.get("gate") or {}
    if {k: moc.get(k) for k in MOC_EXPECT} != MOC_EXPECT or gate.get("applicable") is not True or gate.get("pass") is not True:
        bad.append(f"{arm}: prepare_info の moc が {MOC_EXPECT}・ゲート合格でない ({ {k: moc.get(k) for k in MOC_EXPECT} }, {gate})")
    ev = TM.wall_evidence(d)
    if ev.get("status") != "consistent":
        bad.append(f"{arm}: 壁の証拠 (M4) が {ev.get('status')}")
    q = (d / "MESH_QUALITY.txt").read_text() if (d / "MESH_QUALITY.txt").is_file() else ""
    if not any(l.strip().startswith("VERDICT: PASS") for l in q.splitlines()):
        bad.append(f"{arm}: メッシュ品質が PASS でない")
    for k in ("wall_first_frac", "wall_first_frac_throat"):
        if getattr(mp, k) != ARMS[arm][k]:
            bad.append(f"{arm}: 解決済みの Mesh2DParams.{k} = {getattr(mp, k)!r} (登録 {ARMS[arm][k]!r})")
    if (mp.ni, mp.nj) != NI_NJ:
        bad.append(f"{arm}: 解決済みの (ni, nj) = {(mp.ni, mp.nj)} (登録 {NI_NJ})")
    if int((info.get("mesh") or {}).get("ni", -1)) != NI_NJ[0]:
        bad.append(f"{arm}: prepare_info の mesh.ni が {NI_NJ[0]} でない")
    return bad, ev


def _input_hashes(d: Path) -> dict:
    return {p.name: _sha(p) for p in sorted(d.iterdir()) if p.is_file() and p.name != PREP_JSON}


def _binaries() -> dict:
    out = {}
    for k in ("FORGE_BIN", "REAL_CONVERTER", "FORGE_CONVERTER"):
        v = os.environ.get(k)
        out[k] = {"path": v, "sha256": (_sha(v) if v else None)}
    return out


def prep(out_root: Path, dry: bool = False) -> dict:
    if dry:
        # 乾式確認: forge を起動しない (--resolve-species も含めて)。runner の _ENV は import 時に環境を写すので import より前に設定する
        os.environ["FORGE_BIN"] = str(Path("/nonexistent/forge-dry-run-guard"))
        os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"
    bad = check_problems()
    if bad:
        raise SystemExit("問題の検査が不成立 — 止める:\n  " + "\n  ".join(bad))
    out_root = Path(out_root).resolve()
    dirs = {arm: out_root / ARMS[arm]["prep"] for arm in ARMS}
    for d in dirs.values():
        if d.exists():
            raise SystemExit(f"{d} が既にある (投入スクリプトが消してから作る)")
    RA = TM._load_problem_with_runs()
    infos, mps, evs, steps = {}, {}, {}, {}
    fails = []
    for arm, d in dirs.items():
        prob = C / ARMS[arm]["problem"]
        infos[arm] = RA.prepare(prob, d, nsteps=EV.MAIN_NSTEPS, ic_from=None, cfl_main=TM.CFL_MAIN, implicit_relax=TM.RELAX)
        steps[arm] = EV.config_steps((d / "solverConfig.yaml").read_text())
        fails += [f"{arm}: step 数の同期: {x}" for x in EV.step_problems(steps[arm])]
        S = float(infos[arm]["scale_m"])
        mps[arm] = RA.mesh_params(RA.load_problem(prob), S, ni=321, nj=65, wall_first_frac=5.0e-3)   # prepare と同じ呼び方
        b, evs[arm] = _arm_checks(d, infos[arm], mps[arm], arm)
        fails += b
        print(arm, ARMS[arm]["problem"], "| mesh", json.dumps(infos[arm].get("mesh")), "| M4", evs[arm].get("status"), "|",
              (d / "MESH_QUALITY.txt").read_text().strip().splitlines()[-1])
    cross = cross_check(dirs["A"], dirs["B"], mps["A"], mps["B"], infos["A"], infos["B"])
    fails += cross["failures"]
    ref = ref_check(dirs["A"], dry)
    fails += ref["failures"]
    thermo = {arm: (TM.RUNS / DRY_THERMO_RUN if dry else d) for arm, d in dirs.items()}
    ics = {}
    for arm, d in dirs.items():
        try:
            ics[arm] = ic_check(d, thermo[arm], dry)
        except Exception as e:  # noqa: BLE001 — 検査できないことも不成立として記録する
            ics[arm] = {"VERDICT": "FAIL", "problems": [f"IC の検査を完了できない: {type(e).__name__}: {e}"]}
        if ics[arm]["problems"]:
            fails += [f"{arm}: 起動前の IC の検査: {p}" for p in ics[arm]["problems"]]
    meshes = {}
    for arm, d in dirs.items():
        rec, rows = mesh_record(d, C / ARMS[arm]["problem"], mps[arm], float(infos[arm]["scale_m"]))
        write_sections(d / MESH_CSV, rows)
        rec["hashes"]["sections_csv"] = _sha(d / MESH_CSV)
        rec["binaries"] = _binaries()
        rec["initial_line_source"] = {
            "run": (infos[arm].get("initial_line") or {}).get("run"), "sha256_16": (infos[arm].get("initial_line") or {}).get("sha256_16"),
            "res_sha256": _sha(TM.RUNS / str(yaml.safe_load((C / ARMS[arm]["problem"]).read_text())["geometry"]["initial_line_run"])
                               / str(yaml.safe_load((C / ARMS[arm]["problem"]).read_text())["geometry"]["initial_line_res"]))}
        rec["base_problem_sha256"] = _sha(C / BASE_PROBLEM)
        rec["design_sha256"] = {f: _sha(ROOT / f) for f in DESIGN_FILES}
        rec["git"] = {"head": _git("rev-parse", "HEAD"), "status_design": _git("status", "--porcelain", "--", "design")}
        (d / MESH_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
        meshes[arm] = rec
    summary = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "tool": "euler_t0_e2.prep", "dry": dry,
               "problems_check": "OK", "steps": {"evaluator": {"MAIN_NSTEPS": EV.MAIN_NSTEPS, "OUT_INTERVAL": EV.OUT_INTERVAL,
                                                                "windows": {k: list(v) for k, v in EV.WINDOWS.items()}},
                                                  "config": steps}, "cross_check": cross, "ref_run0153": ref, "ic_check": ics, "mesh": meshes,
               "failures": fails, "VERDICT": "OK" if not fails else "REFUSED"}
    (C / "_band_ab").mkdir(exist_ok=True)
    (C / f"_band_ab/euler_t0_e2_mesh{'_dry' if dry else ''}.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))
    if fails:
        raise RuntimeError("E2 の準備の検査が不成立 — 止める (forge を起動しない):\n  " + "\n  ".join(fails))
    for arm, d in dirs.items():
        info = json.loads((d / "prepare_info.json").read_text())
        info.update(stages="soft", plan=PLAN,
                    e2={"arm": arm, "run": ARMS[arm]["run"], "problem": ARMS[arm]["problem"], "plan_reg_commit": PLAN_REG_COMMIT},
                    ic={"mode": "isentropic", "tool": "paste_isentropic_ic (runner_axismach.prepare)", "VERDICT": ("OK" if not dry else "DRY"),
                        "check": {k: ics[arm].get(k) for k in ("max_abs_dev_K", "T0_min", "T0_max", "thermo_source", "VERDICT")}})
        if dry:
            info["DRY"] = True
        (d / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
        rec = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "arm": arm, "run": ARMS[arm]["run"], "problem": ARMS[arm]["problem"],
               "dry": dry, "ic_check": ics[arm], "cross_check_ok": cross["ok"], "ref_run0153": ref.get("status"),
               "steps": {"config": steps[arm], "MAIN_NSTEPS": EV.MAIN_NSTEPS, "OUT_INTERVAL": EV.OUT_INTERVAL},
               "files_sha256": _input_hashes(d)}
        (d / PREP_JSON).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
        mj = meshes[arm]
        tb = mj["summary_bands_rt"]
        print(f"{arm} 実効のメッシュ: 壁側 (/r_w) スロート {tb['-4<=x<=1'].get('wall_gap_rw')} 下流 x≥17 {tb['x>=17'].get('wall_gap_rw')} | "
              f"軸側 (/r_w) スロート {tb['-4<=x<=1'].get('axis_gap_rw')} 下流 {tb['x>=17'].get('axis_gap_rw')} | "
              f"IC max|T0−1600| {ics[arm]['max_abs_dev_K']:.4g} K ({ics[arm]['VERDICT']})")
    print("照合: A・B の列の x と壁の r が同一", cross["ok"], "| run_0153 との照合", ref.get("status"))
    return summary


def verify_prep(prep_dir: Path, run_dir: Path) -> list:
    """run_dir が prep_dir の複製のまま (E2_PREP.json の入力の sha256 が一致) か。戻り = 不成立の理由。"""
    bad = []
    try:
        rec = json.loads((prep_dir / PREP_JSON).read_text())
        info = json.loads((run_dir / "prepare_info.json").read_text())
    except (OSError, ValueError) as e:
        return [f"{PREP_JSON} / prepare_info.json を読めない ({e})"]
    if rec.get("dry") or info.get("DRY"):
        bad.append("乾式確認 (DRY) の prep")
    if (rec.get("ic_check") or {}).get("VERDICT") != "OK":
        bad.append(f"起動前の IC の検査が OK でない ({(rec.get('ic_check') or {}).get('VERDICT')!r})")
    if run_dir.name != rec.get("run"):
        bad.append(f"run の名前 {run_dir.name} が prep の記録 {rec.get('run')} と違う")
    for n, sha in (rec.get("files_sha256") or {}).items():
        if _sha(run_dir / n) != sha:
            bad.append(f"{run_dir.name}/{n} が prep の記録と違う")
    if not rec.get("files_sha256"):
        bad.append("prep の入力の sha256 の記録が無い")
    return bad


# --- 実行 (soft 段の出力の退避つき) ---------------------------------------------------------------------------------------------
def install_soft_saver(RA, rd: Path) -> list:
    """RA._restart_same_mesh を包む: 呼ばれる前に run 直下の res_* と段のログを rd/_soft_stage/ に写す (1 回だけ。soft 段の引き継ぎ)。
    戻り値は呼ばれた記録のリスト (呼び出し側が 1 回であることを確かめる)。"""
    orig = RA._restart_same_mesh
    calls = []

    def wrapped(res_h5, mesh_h5):
        res_h5 = Path(res_h5)
        if calls:
            raise RuntimeError("段の引き継ぎが 2 回目 — 登録の段は soft → 本段だけ (mid 段は無い)")
        if res_h5.parent.resolve() != rd.resolve() or res_h5.name != f"res_{EV.SOFT_STEPS}.h5":
            raise RuntimeError(f"soft 段の最終場が {res_h5} (期待 {rd}/res_{EV.SOFT_STEPS}.h5)")
        dst = rd / EV.SOFT_DIR
        dst.mkdir(exist_ok=False)
        rec = {"plan": PLAN, "tool": "euler_t0_e2.install_soft_saver", "stage": "soft", "restart_src": res_h5.name, "files": {}}
        srcs = sorted(p for p in rd.glob("res_*") if p.is_file()) + [rd / n for n in SOFT_SAVE_EXTRA if (rd / n).is_file()]
        for p in srcs:
            shutil.copy2(p, dst / p.name)
            a, b = _sha(p), _sha(dst / p.name)
            if a != b:
                raise RuntimeError(f"soft 段の出力の写し {p.name} が元と一致しない")
            rec["files"][p.name] = a
        rec["nozzle_h5_sha256_before_restart"] = _sha(mesh_h5)
        calls.append(str(res_h5))
        try:
            orig(res_h5, mesh_h5)
        finally:
            rec["nozzle_h5_sha256_after_restart"] = _sha(mesh_h5)
            (dst / "SOFT_STAGE_SAVED.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    RA._restart_same_mesh = wrapped
    return calls


def check_run_steps(rd: Path) -> list:
    """起動の直前: run の solverConfig.yaml (本段の config) の nStepOuter・outStepInterval が評価器の MAIN_NSTEPS・OUT_INTERVAL と同じで、
    評価器の窓がその定数から作られた形か。戻り = 食い違いの理由 (空なら同期)。"""
    try:
        return EV.step_problems(EV.config_steps((rd / "solverConfig.yaml").read_text()))
    except Exception as e:  # noqa: BLE001 — 読めないことも食い違いとして止める
        return [f"solverConfig.yaml を読めない: {type(e).__name__}: {e}"]


def run(rd: Path) -> int:
    from forge_design.evaluate import runner_axismach as RA
    info = json.loads((rd / "prepare_info.json").read_text())
    if info.get("DRY"):
        raise SystemExit(f"{rd} は乾式確認 (--dry) の prep から作られている — 回さない")
    e2 = info.get("e2") or {}
    arm = e2.get("arm")
    if arm not in ARMS or rd.name != ARMS[arm]["run"]:
        raise SystemExit(f"{rd}: E2 の腕 {arm!r} とその run 名でない — 回さない")
    if (info.get("ic") or {}).get("VERDICT") != "OK":
        raise SystemExit(f"{rd}: IC の VERDICT が OK でない ({(info.get('ic') or {}).get('VERDICT')!r}) — 回さない")
    if ((info.get("moc") or {}).get("gate") or {}).get("pass") is not True:
        raise SystemExit(f"{rd}: MOC のゲートが合格でない — 回さない")
    if (rd / EV.SOFT_DIR).exists() or any(rd.glob("res_*")):
        raise SystemExit(f"{rd}: 既に出力がある (_soft_stage または res_*) — 回さない")
    sp = check_run_steps(rd)
    if sp:
        raise SystemExit(f"{rd}: 本段の step 数が評価器と同期していない — 回さない:\n  " + "\n  ".join(sp))
    calls = install_soft_saver(RA, rd)
    rc = RA.run_staged(rd, cfl_main=TM.CFL_MAIN, mid_stage=False, stages="soft")
    if len(calls) != 1:
        raise RuntimeError(f"soft 段の出力の退避が {len(calls)} 回 (1 回のはず)")
    print(f"forge exit={rc}")
    return rc


def main(argv) -> int:
    if argv and argv[0] == "check-problems":
        bad = check_problems()
        print("CHECK-PROBLEMS: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if argv and argv[0] == "prep":
        ap = argparse.ArgumentParser(prog="euler_t0_e2.py prep")
        ap.add_argument("--out-root", default=str(C), help="_prep_e2_* を作る場所 (既定は case dir)")
        ap.add_argument("--dry", action="store_true", help="乾式確認 (forge を起動しない)")
        a = ap.parse_args(argv[1:])
        prep(Path(a.out_root), dry=a.dry)
        return 0
    if len(argv) == 3 and argv[0] == "verify-prep":
        bad = verify_prep(Path(argv[1]).resolve(), Path(argv[2]).resolve())
        print("VERIFY-PREP: " + ("OK" if not bad else "FAIL\n  " + "\n  ".join(bad)))
        return 0 if not bad else 1
    if len(argv) == 2 and argv[0] == "run":
        return run(Path(argv[1]).resolve())
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
