"""冷却壁の NS の対 (plan tooling-nozzle-isothermal-wall-chain §5.1 #14、§6 V-c45、2026-10-08 事前登録) の準備・実行。
生産の壁 (run_0167_ns_n012_N2 と同じ物理壁) を冷却壁用の格子 (cold_pair_mesh.py) に載せ、断熱と 300 K の 2 本を回す。
目的: CONTUR の熱閉包 (温度形 / エンタルピー形) のどちらが、冷却による δ_r の変化 δ_E(300 K)/δ_E(断熱) を当てるかを決める。

  壁   : prepare_ns の delta_r_csv = _band_ab/prod_confirm/prep/delta_r_initial.csv (生産の問題の準備、run_0167 の入力をビット同一で再現)。
         物理壁 (wall_repr.json) を同じ準備の壁と密な格子で比べ、|Δr| ≤ 1e-8 r_t を必須にする。
  IC   : 2 本とも生産の dry の最終場 (run_0179_ns_n012_N2_ext の res_20000、通算 100000) から prepare_ns の ic_from (interp_field の
         cross-mesh 移送)。段階起動 full (soft 1 次 cfl 0.5 3000 → mid 1 次 cfl 1 3000 → 本段)。
  本段 : 生産と同じ 2 次 convMethod 1・cfl 1・implicitRelax 0.7。100000 step・5000 ごと。判定窓は 80000〜100000 の 5 枚。
  検査 : メッシュ品質 PASS (--ar-max 5000)・壁距離の変換し直しと一致 (相対 1e-6)・壁の bcond (断熱 wall / 300 K wall_isothermal Ts 300)。

バイナリ: FP64 のビルド (~/forge-wallfit-bin-fp64; sha256 は FORGE_SHA・CONV_SHA、plan §5.1 #16)。

usage (AWS の case dir、FORGE_BIN 等は run_cold_pair.sh が設定する):
  python3 cold_pair.py prep <ad|tw300> <run>
  python3 cold_pair.py run <run>
  python3 cold_pair.py nan-scan <run>
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "design"))
import ns_n012 as NS  # noqa: E402  (小道具: set_config・mesh_quality_ok・reconvert_wall_dist・nan_scan・sha256_file・runner)

PROBLEMS = {"ad": "problem_d155_ns_prod_coldmesh.yaml", "tw300": "problem_d155_ns_prod_coldmesh_tw300.yaml"}
WALL_RUN = "run_0167_ns_n012_N2"                 # 生産の dry の run (記録)
# 生産の壁 (wall_repr.json) と delta_r の表: 生産の問題の準備 (prod_confirm.py、run_0167 の入力をビット同一で再現と確認済み)。
# run_0167 自体は single_bspline の採用前の準備なので wall_repr.json を持たない
WALL_REF = "_band_ab/prod_confirm/prep"
IC_SRC = "run_0179_ns_n012_N2_ext"               # 生産の dry の最終場 (res_20000 = 通算 100000)
IC_RES = "res_20000.h5"
EULER_REF = "run_0174_euler_v5d_M_r1"
MAIN_STEPS, OUT_INT, CFL_MAIN, RELAX = 100000, 5000, 1.0, 0.7
STAGES = "full"
WALL_TOL_RT = 1e-8
RECORD = "COLD_PAIR.json"
# FP64 のビルド (plan §5.1 #16、§6 V-c45): 生産のバイナリの commit e2696d8f0 + typedef 4 行 + 座標読み込み (stod)
FP64_TREE = Path.home() / "forge-wallfit-bin-fp64"
FORGE_SHA = "65be5e28ca1aed9f7683f0ae390adaa103488a39859ca522c2460e9cb71154ea"
CONV_SHA = "ac88861fa04ebdf41e3949fe7c871e12f7356ee59d239bd0c7c204b16d4df872"
FIRST_LAYER_TOL = 1e-6
RUN_RE = re.compile(r"^run_\d{4}_[A-Za-z0-9_.-]+$")


def geom_check(run: Path, scale: float) -> dict:
    """物理壁 (wall_repr.json から復元) を run_0167 の壁と密な x で比べる。"""
    from forge_design.geometry.wall_axismach import load_wall_file
    A = load_wall_file(run); B = load_wall_file(HERE / WALL_REF)
    lo = max(float(A["domain"][0]), float(B["domain"][0])); hi = min(float(A["domain"][1]), float(B["domain"][1]))
    x = np.linspace(lo, hi, 200001)
    d = np.abs(A["physical"].r(x) - B["physical"].r(x))
    out = {"domain": [lo, hi], "domains": [list(map(float, A["domain"])), list(map(float, B["domain"]))],
           "max_abs_dr_rt": float(d.max()), "x_at_max": float(x[np.argmax(d)]), "max_abs_dr_m": float(d.max() * scale),
           "tol_rt": WALL_TOL_RT, "ok": bool(d.max() <= WALL_TOL_RT and A["domain"] == B["domain"])}
    return out


def binary_record() -> dict:
    """FP64 の forge と変換器の sha256 を登録値と照合し、ソースの commit と全差分を記録する。"""
    import os
    import subprocess
    fb = Path(os.environ.get("FORGE_BIN", "")); cv = Path(os.environ.get("REAL_CONVERTER", ""))
    have = {"forge": NS.sha256_file(fb) if fb.is_file() else None, "converter": NS.sha256_file(cv) if cv.is_file() else None}
    if have["forge"] != FORGE_SHA or have["converter"] != CONV_SHA:
        raise SystemExit(f"FP64 のバイナリが登録と違う ({have}; forge {FORGE_SHA[:16]}…、変換器 {CONV_SHA[:16]}…) — 止める")
    head = subprocess.run(["git", "-C", str(FP64_TREE), "log", "--oneline", "-1"], capture_output=True, text=True).stdout.strip()
    diff = subprocess.run(["git", "-C", str(FP64_TREE), "diff"], capture_output=True, text=True).stdout
    return {"forge_bin": str(fb), "forge_sha256": have["forge"], "converter": str(cv), "converter_sha256": have["converter"],
            "source_head": head, "source_diff": diff}


def mesh_checks(run: Path, problem: Path) -> dict:
    """変換後の格子 (nozzle.h5) の検査: スキュー > 0.1 かつ AR > 1000 のセルが 0、第一層厚が生成時の倍精度座標と相対 1e-6 以内。"""
    import h5py
    import yaml
    from forge_design.geometry.wall_axismach import load_wall_file
    from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh
    W = load_wall_file(HERE / WALL_REF); ph = W["physical"]; d0, d1 = (float(v) for v in W["domain"]); rt = float(W["scale_m"])

    class Wall:
        x_in, x_e = d0, d1

        def r(self, x, d=0):
            return ph.r(np.asarray(x), d) if d else ph.r(np.asarray(x))
    m = yaml.safe_load(open(problem))["mesh"]
    prm = Mesh2DParams(ni=m["ni"], nj=m["nj"], wall_first_frac=m["wall_first_frac"], throat_refine=m["throat_refine"],
                       throat_width=m["throat_width"], wall_first_frac_table=m["wall_first_frac_table"],
                       x_density_table=m["x_density_table"], wall_normal_layer=m["wall_normal_layer"], scale=rt)
    coords, quads, _ = generate_axisym_mesh(Wall(), prm)
    with h5py.File(run / "nozzle.h5", "r") as h:
        C = np.array(h["MESH/COORD"]).reshape(-1, 3); dt = str(h["MESH/COORD"].dtype)
    out = {"coord_dtype": dt}
    if C.shape != coords.shape or not np.allclose(C[:, :2], coords[:, :2], rtol=0, atol=1e-6):
        raise SystemExit(f"{run}: nozzle.h5 の節点が生成時の格子と対応しない — 止める")
    ni, nj = prm.ni, prm.nj
    G = coords[:, :2].reshape(ni, nj, 2); H = C[:, :2].reshape(ni, nj, 2)
    rel = np.abs(np.linalg.norm(H[:, -1] - H[:, -2], axis=1) / np.linalg.norm(G[:, -1] - G[:, -2], axis=1) - 1.0)
    P = C[quads][:, :, :2]
    e = np.linalg.norm(np.roll(P, -1, axis=1) - P, axis=2); ar = e.max(1) / e.min(1)
    v1 = np.roll(P, -1, axis=1) - P; v0 = P - np.roll(P, 1, axis=1)
    c = np.sum(-v0 * v1, axis=2) / (np.linalg.norm(v0, axis=2) * np.linalg.norm(v1, axis=2))
    ang = np.degrees(np.arccos(np.clip(c, -1, 1))); sk = np.maximum((ang.max(1) - 90) / 90, (90 - ang.min(1)) / 90)
    bad = int(np.sum((sk > 0.1) & (ar > 1000)))
    out.update(first_layer_rel_err_max=float(rel.max()), ar_max=float(ar.max()), skew_max=float(sk.max()),
               n_skewed_high_ar=bad, ar_max_skewed=float(ar[sk > 0.1].max()))
    if rel.max() > FIRST_LAYER_TOL or bad or dt != "float64":
        raise SystemExit(f"{run}: 格子の検査が不成立 ({out}) — 止める")
    return out


def mesh_quality_strict(run: Path) -> str:
    """MESH_QUALITY.txt の VERDICT が厳密に PASS (SOFT-PASS は不可; ns_n012.mesh_quality_ok は部分文字列で SOFT-PASS も通す)。"""
    t = (run / "MESH_QUALITY.txt").read_text() if (run / "MESH_QUALITY.txt").is_file() else ""
    line = next((l.strip() for l in t.splitlines()[::-1] if "VERDICT" in l), "")
    if not line.startswith("VERDICT: PASS"):
        raise SystemExit(f"{run}/MESH_QUALITY.txt の VERDICT が厳密な PASS でない ({line}) — 止める")
    return line


def bcond_check(run: Path, kind: str) -> str:
    txt = (run / "bcondConfig.yaml").read_text()
    line = next((l for l in txt.splitlines() if "physID: 3" in l), "")
    ok = (re.search(r"kind:\s*wall_isothermal\s*,", line) and re.search(r"Ts:\s*300(\.0)?\b", line)) if kind == "tw300" \
        else re.search(r"kind:\s*wall\s*,", line)
    if not ok:
        raise SystemExit(f"{run}/bcondConfig.yaml の壁の行が {kind} と合わない: {line!r} — 止める")
    return line.strip()


def prep(kind: str, run: Path) -> dict:
    NS.check_dry_env(False)
    binrec = binary_record()
    if kind not in PROBLEMS:
        raise SystemExit(f"kind は {list(PROBLEMS)} のどれか")
    if not RUN_RE.match(run.name):
        raise SystemExit(f"run 名 {run.name!r} が run_NNNN_<slug> でない")
    if run.exists():
        raise SystemExit(f"{run} が既にある — 止める (既存 run は消さない)")
    problem = HERE / PROBLEMS[kind]
    dr_csv = HERE / WALL_REF / "delta_r_initial.csv"
    src = HERE / IC_SRC
    rs = NS.res_files(src)
    if not rs or rs[-1].name != IC_RES:
        raise SystemExit(f"IC のドナー {src.name} の最後の res が {rs[-1].name if rs else None} ({IC_RES} であること) — 止める")
    info = NS._prepare_ns(problem, run, nsteps=12000, ic_from=src, delta_r_csv=str(dr_csv), offset="radial",
                          euler_ref=str(HERE / EULER_REF), cfl_main=5.0, implicit_relax=RELAX)
    NS.set_config(run, {NS.CFL: "1.0", NS.CFLP: "1.0", NS.NSTEP: MAIN_STEPS, NS.OUTINT: OUT_INT},
                  {NS.CFL: CFL_MAIN, NS.CFLP: CFL_MAIN, NS.NSTEP: MAIN_STEPS, NS.OUTINT: OUT_INT, NS.RELAXP: RELAX, NS.CONVP: 1})
    info["stages"] = {"stages": STAGES, "ramp": None, "ramp_steps": 1000}
    NS.jdump(run / "prepare_info.json", info)
    mq = mesh_quality_strict(run)
    mck = mesh_checks(run, problem)
    scale = float(info["scale_m"])
    geo = geom_check(run, scale)
    if not geo["ok"]:
        raise SystemExit(f"{run}: 物理壁が {WALL_RUN} と合わない ({geo}) — 止める")
    wd_new, coord_new, conv = NS.reconvert_wall_dist(run)
    import h5py
    with h5py.File(run / "nozzle.h5", "r") as f:
        wd = np.array(f["VALUE/wall_dist"]); coord = np.array(f["MESH/COORD"])
    if wd.shape != wd_new.shape or not np.array_equal(coord, coord_new):
        raise SystemExit(f"{run}: 変換し直した格子が nozzle.h5 と合わない — 止める")
    wd_rel = float(np.max(np.abs(wd - wd_new) / np.maximum(np.abs(wd_new), 1e-30)))
    if wd_rel > NS.WALLDIST_RTOL:
        raise SystemExit(f"{run}: wall_dist が変換し直した値と相対 {wd_rel:.3e} ずれる — 止める")
    bl = bcond_check(run, kind)
    rec = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #14、§6 V-c45", "tool": "cold_pair.py prep",
           "created": NS.now(), "git_head": NS.git_head(), "kind": kind, "problem": problem.name, "problem_sha256": NS.sha256_file(problem),
           "delta_r_csv": f"{WALL_REF}/delta_r_initial.csv", "delta_r_csv_sha256": NS.sha256_file(dr_csv),
           "ic": {"src_run": src.name, "src_res": IC_RES, "src_res_sha256": NS.sha256_file(rs[-1]), "via": "prepare_ns(ic_from) → interp_field.py"},
           "euler_ref": EULER_REF, "stages": STAGES, "main_steps": MAIN_STEPS, "out_interval": OUT_INT, "cfl_main": CFL_MAIN,
           "implicit_relax": RELAX, "mesh_quality": mq, "mesh": info.get("mesh"), "wall_thermal": info.get("wall_thermal"),
           "binary": binrec, "mesh_checks": mck, "geometry_vs_production": geo, "wall_dist_rel_max": wd_rel, "converter": conv, "bcond_wall": bl,
           "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / RECORD, rec)
    print(f"[cold_pair prep] {kind} {run.name}: {mq}; 第一層 {mck['first_layer_rel_err_max']:.1e}; 高 AR のスキュー {mck['n_skewed_high_ar']}; "
          f"壁の差 {geo['max_abs_dr_rt']:.2e} r_t; wall_dist 相対 {wd_rel:.1e}; {bl}")
    return rec


def run_one(run: Path) -> int:
    rec = NS.jload(run / RECORD)
    NS.check_dry_env(False)
    if any((run / f).exists() for f in ("stage_manifest.json", "residual_history.csv")) or NS.res_files(run):
        raise SystemExit(f"{run}: 既に出力がある — 回さない")
    if NS.sha256_file(run / "nozzle.h5") != rec.get("nozzle_sha256_after_prep"):
        raise SystemExit(f"{run}: nozzle.h5 が準備の後に変わった — 回さない")
    binary_record()                                   # 回す直前にもバイナリを照合
    rc = NS.runner().run_staged_ns(run, stages=rec.get("stages", STAGES))
    last = NS.res_files(run)
    print(f"forge exit={rc} last_res={last[-1].name if last else None}")
    (run / "RUN_RC").write_text(f"{rc}\n")
    return rc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep"); p.add_argument("kind"); p.add_argument("run")
    p = sp.add_parser("run"); p.add_argument("run")
    p = sp.add_parser("nan-scan"); p.add_argument("run")
    a = ap.parse_args()
    if a.cmd == "prep":
        prep(a.kind, HERE / a.run)
    elif a.cmd == "run":
        sys.exit(run_one(HERE / a.run))
    elif a.cmd == "nan-scan":
        out = NS.nan_scan(HERE / a.run)
        NS.jdump(HERE / a.run / "NAN_SCAN.json", out)
        print(json.dumps({k: out.get(k) for k in ("first_nonfinite", "first_bad_field")}, ensure_ascii=False))
