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
  python3 cold_pair.py extract <run>            (AWS: 判定窓の δ_E・感度・壁温・Q_w → _band_ab/cold_pair/extract_<run>.npz)
  python3 cold_pair.py gates <ad_run> <tw_run>  (AWS: check_convergence・残差の床・壁解像 → gates_aws.json、y1p_<run>_wall.csv)
  python3 cold_pair.py judge <ad_run> <tw_run>  (手元: V-c45 の判定 → V_c45.json; NaN 検査は _band_ab/cold_pair/NAN_SCAN_<run>.json に写しておく)
  python3 cold_pair.py extract-b <run>          (AWS: 抽出の A/B の B 腕 = 一定 x 断面 → extractB_<run>.npz)
  python3 cold_pair.py ab-compare <ad> <tw>     (手元: 抽出の A/B の判定 → extract_ab.json)
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


# --- 判定 (plan §6 V-c45) --------------------------------------------------------------------------------------------
WINDOW = (80000, 85000, 90000, 95000, 100000)
XJ = np.arange(40.0, 94.0 + 1e-9, 1.0)                 # 定常性の評価点 (1 r_t ごと)
XE = np.arange(40.0, 94.0 + 1e-9, 0.25)                # 判定の評価点 (0.25 r_t ごと)
OUTD = HERE / "_band_ab" / "cold_pair"
RES_COLS = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roK", "rms_roOmega")


def _wall_dump(run: Path, step: int, scale: float, axisym: bool = True):
    """res_wall_3_<step>.h5 → x [r_t] 順の (x, T_w, q_w, 面積の重み)。"""
    import h5py
    sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda/tools"))
    from check_wall_resolution import point_weights
    with h5py.File(run / f"res_wall_3_{step}.h5", "r") as h:
        xyz = np.array(h["MESH/COORD"]).reshape(-1, 3).astype(float); conne = np.array(h["MESH/CONNE"])
        Ts = np.array(h["VALUE/Ts"]).astype(float); qw = np.array(h["VALUE/qwall"]).astype(float)
    w = point_weights(xyz, conne, len(Ts), False, axisym)
    o = np.argsort(xyz[:, 0])
    return xyz[o, 0] / scale, Ts[o], qw[o], (w[o] if w is not None else None)


def extract(run: Path) -> Path:
    """AWS: 判定窓の 5 枚から δ_E (帯 E、extract_and_merge の delta_E = 未緩和の平滑化抽出、300 K は符号付き)、抽出の感度、
    壁温、壁の熱流束の積分 Q_w を取り出して _band_ab/cold_pair/extract_<run>.npz に保存する。"""
    import os
    import shutil
    import tempfile
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    from forge_design.metrics.deltastar import deltastar_from_core_matched_euler
    info = NS.jload(run / "prepare_info.json"); scale = float(info["scale_m"])
    out = {}
    steps_all = sorted(NS.step_of(f) for f in NS.res_files(run))
    window = tuple(steps_all[-5:])                     # 判定窓 = その run の最後の 5 枚 (本段 80000〜100000、延長は延長の最後の 5 枚)
    for k, step in enumerate(window):
        with tempfile.TemporaryDirectory() as td:
            dd = Path(td) / "ns"; dd.mkdir()
            for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
                shutil.copy(run / f, dd / f)
            os.symlink((run / "nozzle.h5").resolve(), dd / "nozzle.h5")
            os.symlink((run / f"res_{step}.h5").resolve(), dd / f"res_{step}.h5")
            x = deltastar_from_core_matched_euler(dd, HERE / EULER_REF, band_select="edge")
            extract_and_merge(dd, HERE / EULER_REF, band_select="edge")
            nx = read_delta_r_next(dd / "delta_r_next.csv")
        sens = np.asarray(x["delta_r_sens"], dtype=float)
        out[f"s{k}_x"] = np.asarray(x["x"]); out[f"s{k}_raw"] = np.asarray(x["delta_r_raw"]); out[f"s{k}_sens"] = sens
        out[f"s{k}_ok"] = np.asarray(x["ok"]); out[f"s{k}_mx"] = np.asarray(nx["x_rt"]); out[f"s{k}_dE"] = np.asarray(nx["delta_E"])
        wx, Tw, qw, ww = _wall_dump(run, step, scale)
        out[f"s{k}_wx"] = wx; out[f"s{k}_Tw"] = Tw; out[f"s{k}_qw"] = qw
        out[f"s{k}_Qw"] = np.array(float(np.sum(qw * ww)) if ww is not None else np.nan)
        print(f"[extract] {run.name} res_{step}: δ 列 {len(x['x'])}、sens {sens.shape}、Q_w {float(out[f's{k}_Qw']):.4g} W", flush=True)
    OUTD.mkdir(parents=True, exist_ok=True)
    p = OUTD / f"extract_{run.name}.npz"
    np.savez(p, steps=np.array(window), scale=np.array(scale), **out)
    return p


def residual_floor(run: Path, last_steps: int = 5000) -> dict:
    import csv
    f = run / "residual_history.csv"
    rows = [r for r in csv.DictReader(open(f)) if r.get("phase", "").strip() == "outer_end"]
    smax = max(int(r["step"]) for r in rows)
    sel = [r for r in rows if int(r["step"]) > smax - last_steps]
    return {c: float(np.median([float(r[c]) for r in sel])) for c in RES_COLS if c in rows[0]}


def gates_aws(ad: Path, tw: Path) -> dict:
    """AWS: check_convergence --segment、残差の床 (生産の run_0179 との比較)、壁解像 (面積、領域別の分布 CSV)。"""
    import subprocess
    tools = HERE.parents[1] / "solver_density_cuda/tools"
    out = {"floor_ref": residual_floor(HERE / IC_SRC)}
    for run in (ad, tw):
        r = {}
        cc = subprocess.run([sys.executable, str(tools / "check_convergence.py"), str(run), "--segment"], capture_output=True, text=True)
        r["convergence"] = [l for l in cc.stdout.splitlines() if "->" in l or "VERDICT" in l or "RISING" in l][-6:]
        r["convergence_rc"] = cc.returncode
        r["floor"] = residual_floor(run)
        r["floor_ok"] = all(r["floor"][c] <= out["floor_ref"][c] for c in r["floor"] if c in out["floor_ref"])
        wr = subprocess.run([sys.executable, str(tools / "check_wall_resolution.py"), str(run), "--groups", "wall", "--weight", "area",
                             "--over-frac", "5", "--profile-csv", str(OUTD / f"y1p_{run.name}.csv")], capture_output=True, text=True)
        r["wall_resolution"] = [l for l in wr.stdout.splitlines() if l.strip()][-4:]
        out[run.name] = r
    NS.jdump(OUTD / "gates_aws.json", out)
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


def judge(ad_name: str, tw_name: str) -> dict:
    """手元: 抽出 (extract_*.npz)・ゲート (gates_aws.json、y1p_*.csv)・CONTUR の予測 (../delta_contur/cooling_ratio_predictions.json) から
    V-c45 を判定して _band_ab/cold_pair/V_c45.json に書く。"""
    sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda/tools"))
    from check_quasisteady import classify
    A = np.load(OUTD / f"extract_{ad_name}.npz"); B = np.load(OUTD / f"extract_{tw_name}.npz")
    rt_mm = float(A["scale"]) * 1e3
    n = len(A["steps"])
    dE = {}
    for lab, Z in (("ad", A), ("tw", B)):
        dE[lab] = np.array([np.interp(XE, Z[f"s{k}_mx"], Z[f"s{k}_dE"]) for k in range(n)])
    Rk = dE["tw"] / dE["ad"]                                       # 各枚の比
    R = dE["tw"].mean(0) / dE["ad"].mean(0)                        # 5 枚平均の比
    out = {"x": XE.tolist(), "R_NS": R.tolist(), "gates": {}}
    # --- ゲート 3: 準定常 (classify、5 枚全部、drift・osc ≤ 0.1 %) ---
    ji = [int(np.argmin(np.abs(XE - x))) for x in XJ]
    def qs(series, name):
        # 5 枚すべて (tail_frac 1.0)、drift・osc の許容 0.1 %、最少 5 枚 (plan §6 V-c45 ゲート 3)
        return classify(np.array(B["steps"], dtype=float), np.asarray(series, dtype=float), 1.0, 0.001, 0.001, 5)
    bad = []
    for lab in ("ad", "tw"):
        for i in ji:
            v = qs(dE[lab][:, i], f"dE_{lab}")
            if v[0] != "STEADY":
                bad.append((f"δ_E {lab} x={XE[i]:.0f}", v[0]))
    for i in ji:
        v = qs(Rk[:, i], "R")
        if v[0] != "STEADY":
            bad.append((f"R x={XE[i]:.0f}", v[0]))
    TwA = np.array([np.interp(XJ, A[f"s{k}_wx"], A[f"s{k}_Tw"]) for k in range(n)])
    for j, x in enumerate(XJ):
        v = qs(TwA[:, j], "Tw")
        if v[0] != "STEADY":
            bad.append((f"断熱壁温 x={x:.0f}", v[0]))
    Qw = np.array([float(B[f"s{k}_Qw"]) for k in range(n)])
    vq = qs(Qw, "Qw")
    if vq[0] != "STEADY":
        bad.append(("Q_w (300 K)", vq[0]))
    out["gates"]["quasisteady"] = {"ok": not bad, "not_steady": bad[:40], "n_not_steady": len(bad), "Qw_W": Qw.tolist(), "Qw_verdict": vq[0]}
    # --- ゲート 2・4 (AWS の結果) ---
    G = json.loads((OUTD / "gates_aws.json").read_text())
    out["gates"]["convergence"] = {k: {"lines": G[k]["convergence"], "floor": G[k]["floor"], "floor_ok": G[k]["floor_ok"]}
                                   for k in (ad_name, tw_name)}
    out["gates"]["convergence"]["floor_ref_run_0179"] = G["floor_ref"]
    def conv_one(k):
        lines = G[k].get("convergence") or []
        txt = " ".join(lines)
        if not lines or "->" not in txt:
            return False                                           # 判定の記録が無い = 不合格
        if any(w in txt for w in ("RISING", "DIVERGED", "INDETERMINATE", "判定不能")):
            return False
        allowed = ("PASS" in txt) or ("stalled/plateau" in txt)     # plateau は許す (VERDICT の文言はそのまま残す)
        return bool(allowed and G[k]["floor_ok"])
    conv_ok = all(conv_one(k) for k in (ad_name, tw_name))
    # ゲート 1: NaN 検査 (各 run の NAN_SCAN.json を手元に写したもの; 無ければ不合格)
    nan_ok = True; nan_rec = {}
    for k in (ad_name, tw_name):
        pth = OUTD / f"NAN_SCAN_{k}.json"
        if not pth.is_file():
            nan_ok = False; nan_rec[k] = "記録なし (不合格)"; continue
        z = json.loads(pth.read_text())
        ok = z.get("first_nonfinite") is None and z.get("first_bad_field") is None and bool(z.get("fields_scanned"))
        nan_ok &= ok; nan_rec[k] = {"ok": ok, "fields": len(z.get("fields_scanned") or [])}
    out["gates"]["nan"] = nan_rec
    wr = {}
    for run_name in (ad_name, tw_name):
        P = np.loadtxt(OUTD / f"y1p_{run_name}_wall.csv", delimiter=",", skiprows=1)
        x = P[:, 0] / (rt_mm * 1e-3); yp = P[:, 4]; w = P[:, 9]
        regs = {"all": np.ones_like(x, bool), "test[40,94]": (x >= 40) & (x <= 94), "[-1,40)": (x >= -1) & (x < 40),
                "contraction x<-1 (入口の角を除く)": (x < -1) & (x > x.min() + 0.05), "入口の角 0.05 r_t": x <= x.min() + 0.05}
        lim = {"all": 5.0, "test[40,94]": 1.0, "[-1,40)": 5.0, "contraction x<-1 (入口の角を除く)": 10.0, "入口の角 0.05 r_t": None}
        rr = {}
        for k, m in regs.items():
            good = m & np.isfinite(yp)
            pct = float(100 * w[good & (yp > 1)].sum() / w[good].sum()) if w[good].sum() > 0 else float("nan")
            rr[k] = {"over1_area_pct": pct, "limit_pct": lim[k], "ok": (lim[k] is None) or pct <= lim[k],
                     "y1p_max": float(np.nanmax(yp[m])) if m.any() else None,
                     "x_of_max": float(x[m][np.nanargmax(yp[m])]) if m.any() else None}
        wr[run_name] = rr
    out["gates"]["wall_resolution"] = wr
    wall_ok = all(v["ok"] for v in wr[tw_name].values())
    # --- 判定の量と不確かさ ---
    u_t = np.max(np.abs(Rk / R - 1.0), axis=0)
    def ext_rel(Z):
        d = []
        for k in range(n):
            raw = np.interp(XE, Z[f"s{k}_x"], Z[f"s{k}_raw"])
            sens = np.atleast_2d(Z[f"s{k}_sens"])
            if sens.shape[0] != len(Z[f"s{k}_x"]):
                sens = sens.T
            dev = np.max(np.abs(np.array([np.interp(XE, Z[f"s{k}_x"], sens[:, j]) for j in range(sens.shape[1])]) / raw - 1.0), axis=0)
            d.append(dev)
        return np.max(d, axis=0)
    u_ext_raw = ext_rel(A) + ext_rel(B)
    win = int(round(2.0 / 0.25))
    ker = np.ones(win + 1)
    u_ext = np.convolve(u_ext_raw, ker, mode="same") / np.convolve(np.ones_like(u_ext_raw), ker, mode="same")   # 端は有効点数で割る
    u_c = 0.001
    U = u_t + u_ext + u_c
    out.update(u_t=u_t.tolist(), u_ext=u_ext.tolist(), u_c=u_c, U=U.tolist())
    pred = json.loads((HERE / "_band_ab/delta_contur/cooling_ratio_predictions.json").read_text())
    px = np.array(pred["x"])
    dec = {}
    for kf in pred["k_f"]:
        e = {}
        for arm in ("A", "B"):
            Rm = np.interp(XE, px, np.array(pred["arms"][f"{arm}_kf{kf:.6f}"]["R"]))
            lo_R, hi_R = R * (1 - U), R * (1 + U)
            c1, c2 = Rm / hi_R - 1.0, Rm / lo_R - 1.0                    # R が区間を動くときの R_m/R − 1 の範囲
            lo = np.where(c1 * c2 <= 0, 0.0, np.minimum(np.abs(c1), np.abs(c2)))
            hi = np.maximum(np.abs(c1), np.abs(c2))
            e[arm] = {"lo": float(lo.max()), "hi": float(hi.max()), "point": float(np.max(np.abs(Rm / R - 1.0))),
                      "x_at_point_max": float(XE[np.argmax(np.abs(Rm / R - 1.0))])}
        if e["B"]["hi"] < e["A"]["lo"]:
            v = "エンタルピー形を支持"; win_ = "B"
        elif e["A"]["hi"] < e["B"]["lo"]:
            v = "温度形を支持"; win_ = "A"
        else:
            v = "判定保留 (区間が分離しない)"; win_ = None
        fit = None if win_ is None else ("この格子で冷却の効果を 1 % 以内で当てる" if e[win_]["hi"] <= 0.01 else "冷却の効果は当てきれない (1 % 超)")
        dec[f"kf{kf:.6f}"] = {"e": e, "verdict": v, "winner": win_, "fit": fit}
    out["decision_by_kf"] = dec
    wins = {d["winner"] for d in dec.values()}
    gates_ok = nan_ok and out["gates"]["quasisteady"]["ok"] and conv_ok and wall_ok
    if not gates_ok:
        final = "判定不能 (ゲート不成立)"
    elif len(wins) == 1 and None not in wins:
        final = next(iter(dec.values()))["verdict"]
    else:
        final = "判定保留 (k_f で判定が違う、または区間が分離しない)"
    out["gates_ok"] = {"nan": nan_ok, "quasisteady": out["gates"]["quasisteady"]["ok"], "convergence": conv_ok, "wall_resolution": wall_ok}
    out["VERDICT"] = final
    # 記録のみ: 格子・精度を替えた感度 (この格子の断熱 vs 生産の FP32)
    prod = HERE / "_band_ab/delta_contur/extract.npz"
    if prod.is_file():
        Z = np.load(prod); nps = len(json.loads((HERE / "_band_ab/delta_contur/extract.json").read_text())["snaps"])
        dprod = np.array([np.interp(XE, Z[f"snap{i}_merged_x"], Z[f"snap{i}_merged_dE"]) for i in range(nps)]).mean(0)
        out["grid_precision_sensitivity_ad"] = {"max_rel": float(np.max(np.abs(dE["ad"].mean(0) / dprod - 1.0))),
                                                "test_mean_rel": float(np.mean(dE["ad"].mean(0) / dprod - 1.0))}
    out["delta_E_mean_mm"] = {"ad": (dE["ad"].mean(0) * rt_mm).tolist(), "tw": (dE["tw"].mean(0) * rt_mm).tolist()}
    NS.jdump(OUTD / "V_c45.json", out)
    brief = {k: out[k] for k in ("VERDICT", "gates_ok")}
    brief["decision_by_kf"] = dec
    brief["n_not_steady"] = out["gates"]["quasisteady"]["n_not_steady"]
    brief["U_range"] = [float(U.min()), float(U.max())]
    print(json.dumps(brief, indent=1, ensure_ascii=False))
    return out


def prep_ext(src: Path, run: Path, steps: int) -> dict:
    """延長 (plan §6 V-c45「延長の決め方」2026-10-08): src の最後の res から restart_field (同一格子、ビット一致、FP64 の型のまま) で
    新しい run に継ぎ、設定は nStepOuter だけを変えて steps step (5000 ごと)、段なし。"""
    import shutil
    import subprocess
    NS.check_dry_env(False)
    binrec = binary_record()
    if not RUN_RE.match(run.name) or run.exists():
        raise SystemExit(f"{run} の名前が不正か既にある — 止める")
    srec = NS.jload(src / RECORD)
    rs = NS.res_files(src)
    if not rs:
        raise SystemExit(f"{src} に res が無い")
    src_h5 = rs[-1]
    ys = NS.yaml_strict()
    ptext = (src / "solverConfig.yaml").read_text()
    ctext = ys.replace_scalars(ptext, {NS.NSTEP: str(int(steps))})
    if not set(NS.MK.diff_paths(ys.load(ptext), ys.load(ctext))) <= {NS.NSTEP}:   # 同じ step 数なら差は空
        raise SystemExit("延長の solverConfig に nStepOuter 以外の差がある — 止める")
    run.mkdir(parents=True)
    for fn in NS.EXT_COPY + ("wall_repr.json", "bcondConfig.yaml", "species_meta.yaml"):
        if (src / fn).is_file():
            shutil.copy2(src / fn, run / fn)
    for p in sorted(src.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
    (run / "solverConfig.yaml").write_text(ctext)
    cmd = [sys.executable, str(NS.TOOLS / "restart_field.py"), str(src_h5), str(run / "nozzle.h5"), "--dst-run", str(run), "--keep-src-dtype"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=NS.runner()._ENV)
    (run / "restart_field.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or "ビット一致" not in (r.stdout + r.stderr):
        print((r.stdout + r.stderr)[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める")
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src.name, restart_from=f"{src.name}/{src_h5.name}")
    NS.jdump(run / "prepare_info.json", info)
    rec = {**{k: srec[k] for k in ("plan", "kind", "problem", "problem_sha256", "delta_r_csv", "euler_ref", "cfl_main", "implicit_relax",
                                   "out_interval", "mesh_checks", "geometry_vs_production", "wall_thermal", "bcond_wall")},
           "tool": "cold_pair.py prep-ext", "created": NS.now(), "git_head": NS.git_head(), "binary": binrec, "stages": "none",
           "parent": src.name, "parent_res": src_h5.name, "parent_res_sha256": NS.sha256_file(src_h5), "ext_steps": int(steps),
           "restart_field_tail": (r.stdout + r.stderr).strip().splitlines()[-1:], "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / RECORD, rec)
    print(f"[cold_pair prep-ext] {run.name} ← {src.name}/{src_h5.name}: {rec['restart_field_tail']}")
    return rec


def extract_b(run: Path) -> Path:
    """AWS: 抽出の A/B の B 腕 (plan §5.1 #21、codex diagnose 2026-10-08): 判定窓の 5 枚で、列を一定 x の断面に補間してから
    同じ抽出 (帯 E、extract_and_merge の delta_E) を行う。A 腕は extract() の結果 (列をそのまま使う)。
    → _band_ab/cold_pair/extractB_<run>.npz (s{k}_mx, s{k}_dE)。"""
    import os
    import shutil
    import tempfile
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    out = {}
    steps_all = sorted(NS.step_of(f) for f in NS.res_files(run))
    window = tuple(steps_all[-5:])
    for k, step in enumerate(window):
        with tempfile.TemporaryDirectory() as td:
            dd = Path(td) / "ns"; dd.mkdir()
            for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
                shutil.copy(run / f, dd / f)
            os.symlink((run / "nozzle.h5").resolve(), dd / "nozzle.h5")
            os.symlink((run / f"res_{step}.h5").resolve(), dd / f"res_{step}.h5")
            extract_and_merge(dd, HERE / EULER_REF, band_select="edge", remap_constant_x=True)
            nx = read_delta_r_next(dd / "delta_r_next.csv")
        out[f"s{k}_mx"] = np.asarray(nx["x_rt"]); out[f"s{k}_dE"] = np.asarray(nx["delta_E"])
        print(f"[extract-B] {run.name} res_{step}", flush=True)
    OUTD.mkdir(parents=True, exist_ok=True)
    p = OUTD / f"extractB_{run.name}.npz"
    np.savez(p, steps=np.array(window), **out)
    return p


def ab_compare(ad_name: str, tw_name: str) -> dict:
    """手元: 抽出の A/B の判定 (事前登録 §5.1 #21): 各 5 枚の対で R = δ_E(300 K)/δ_E(断熱) を A・B で作り、試験部 [40, 94] (0.25 刻み) の
    max |R_B/R_A − 1| が 1 % 以上なら「抽出の座標の不整合の寄与を支持」、1 % 未満なら棄却。各腕の δ_E の A/B の差も記録する。"""
    A = {lab: np.load(OUTD / f"extract_{n}.npz") for lab, n in (("ad", ad_name), ("tw", tw_name))}
    B = {lab: np.load(OUTD / f"extractB_{n}.npz") for lab, n in (("ad", ad_name), ("tw", tw_name))}
    nA = len(A["ad"]["steps"])
    d = lambda Z, k: np.interp(XE, Z[f"s{k}_mx"], Z[f"s{k}_dE"])  # noqa: E731
    rel = []; arm = {"ad": [], "tw": []}
    for k in range(nA):
        RA = d(A["tw"], k) / d(A["ad"], k); RB = d(B["tw"], k) / d(B["ad"], k)
        rel.append(np.abs(RB / RA - 1.0))
        for lab in ("ad", "tw"):
            arm[lab].append(np.abs(d(B[lab], k) / d(A[lab], k) - 1.0))
    rel = np.array(rel)
    mx = float(rel.max()); i = np.unravel_index(int(np.argmax(rel)), rel.shape)
    out = {"criterion": "max |R_B/R_A − 1| (5 枚・試験部) ≥ 1 % で支持", "max_rel_R": mx, "at": {"snapshot": int(i[0]), "x": float(XE[i[1]])},
           "mean_rel_R_by_snapshot": rel.mean(1).tolist(),
           "delta_E_AB_max_rel": {lab: float(np.max(arm[lab])) for lab in arm},
           "delta_E_AB_mean_rel": {lab: float(np.mean(arm[lab])) for lab in arm},
           "R_A_mean": (np.mean([d(A["tw"], k) / d(A["ad"], k) for k in range(nA)], axis=0)[[0, 80, 160, 216]]).tolist(),
           "R_B_mean": (np.mean([d(B["tw"], k) / d(B["ad"], k) for k in range(nA)], axis=0)[[0, 80, 160, 216]]).tolist(),
           "x_report": [float(XE[j]) for j in (0, 80, 160, 216)],
           "VERDICT": ("抽出の座標の不整合の寄与を支持 (≥ 1 %)" if mx >= 0.01 else "第 1 仮説を棄却 (< 1 %)")}
    NS.jdump(OUTD / "extract_ab.json", out)
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


def theta_diag(run: Path, step: int | None = None) -> Path:
    """AWS: NS の θ の抽出 (plan §5.1 #24、記録)。δ_E と同じ抽出 (帯 E、一定 x の断面) から帯の外縁 y_b を取り、断面の NS の分布で
    外縁の値 (ρ_e, u_e) を縁として、面積で等価な厚さを作る:
      質量: 2π∫_{r_w−y_b}^{r_w} (ρ_e u_e − ρu) r dr = π(r_w² − (r_w − δ_loc)²) ρ_e u_e
      運動量: 2π∫ ρu (u_e − u) r dr = π(r_w² − (r_w − θ_r)²) ρ_e u_e²
    比 (300 K/断熱) だけを CONTUR と比べる (δ_loc/θ_r は CONTUR の H と定義が違う)。→ _band_ab/cold_pair/theta_<run>_<step>.npz"""
    import os
    import shutil
    import tempfile
    import h5py
    from forge_design.metrics.deltastar import deltastar_from_core_matched_euler
    info = NS.jload(run / "prepare_info.json"); S = float(info["scale_m"]); ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    if step is None:
        step = sorted(NS.step_of(f) for f in NS.res_files(run))[-1]
    with tempfile.TemporaryDirectory() as td:
        dd = Path(td) / "ns"; dd.mkdir()
        for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
            shutil.copy(run / f, dd / f)
        os.symlink((run / "nozzle.h5").resolve(), dd / "nozzle.h5")
        os.symlink((run / f"res_{step}.h5").resolve(), dd / f"res_{step}.h5")
        ex = deltastar_from_core_matched_euler(dd, HERE / EULER_REF, band_select="edge", remap_constant_x=True)
    with h5py.File(run / "nozzle.h5") as h:
        nc = h["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(run / f"res_{step}.h5") as h:
        ro = h["/VALUE/ro"][:].astype(float); ux = h["/VALUE/Ux"][:].astype(float)
    x = nc[:, 0].reshape(ni, nj) / S; r = nc[:, 1].reshape(ni, nj) / S
    ro = ro.reshape(ni, nj); ux = ux.reshape(ni, nj)
    xt = x[:, 0].copy()
    R = np.empty_like(r); RO = np.empty_like(ro); UX = np.empty_like(ux)
    for j in range(nj):                                   # 一定 x の断面 (抽出の B 腕と同じ)
        R[:, j] = np.interp(xt, x[:, j], r[:, j]); RO[:, j] = np.interp(xt, x[:, j], ro[:, j]); UX[:, j] = np.interp(xt, x[:, j], ux[:, j])
    xs = np.asarray(ex["x"]); yb = np.asarray(ex["band_y_b"], dtype=float)
    th = np.full(len(xs), np.nan); dl = np.full(len(xs), np.nan)
    for i, xv in enumerate(xs):
        k = int(np.argmin(np.abs(xt - xv)))
        if not np.isfinite(yb[i]) or abs(xt[k] - xv) > 1e-9:
            continue
        rr = R[k]; rw = rr[-1]; rb = rw - yb[i]
        if rb <= rr[0]:
            continue
        rf = np.linspace(rb, rw, 4001)
        rho = np.interp(rf, rr, RO[k]); u = np.interp(rf, rr, UX[k])
        re_, ue_ = rho[0], u[0]
        qm = np.trapezoid((re_ * ue_ - rho * u) * rf, rf)
        qq = np.trapezoid(rho * u * (ue_ - u) * rf, rf)
        dl[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qm / (re_ * ue_), 0.0))
        th[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qq / (re_ * ue_ ** 2), 0.0))
    OUTD.mkdir(parents=True, exist_ok=True)
    p = OUTD / f"theta_{run.name}_{step}.npz"
    np.savez(p, x=xs, theta_r=th, delta_loc=dl, delta_r_raw=np.asarray(ex["delta_r_raw"]), band_y_b=yb, step=np.array(step), scale=np.array(S))
    print(f"[theta] {run.name} res_{step}: θ の有効点 {int(np.isfinite(th).sum())}/{len(th)}")
    return p


# --- SU2 との照合 (plan §5.1 #26): 壁・格子は TP の対と同じ、CFD だけ CPG -------------------------------------------
CPG_GAMMA, CPG_CP = 1.27354, 1360.0
CPG_PHYS = ("physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, "
            f"cp: {CPG_CP}, gamma: {CPG_GAMMA}}}")
TURB_PLAIN = 'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: 0, wallTreatmentSST: 0, turbulentPrandtl: 0.9}'


def prep_cpg(mesh_run: Path, field_res: Path, run: Path, variant: str, main_steps: int = 100000, out_int: int = 20000) -> dict:
    """準備済みの TP の run (mesh_run: 格子・壁・BC) を複製して CFD を CPG に書き換える。初期値は field_res (収束した TP の解) の
    ρ・U・P・k・ω を CPG の保存量 (roe = P/(γ−1) + ½ρ|U|²) に組み直したもの。variant: plain (dilat 0・KL 0) / dilat2 (生産と同じ)。"""
    import re as _re
    import shutil
    import h5py
    NS.check_dry_env(False)
    binrec = binary_record()
    if variant not in ("plain", "dilat2"):
        raise SystemExit("variant は plain | dilat2")
    if not RUN_RE.match(run.name) or run.exists():
        raise SystemExit(f"{run} の名前が不正か既にある — 止める")
    run.mkdir(parents=True)
    for fn in ("bcondConfig.yaml", "solverConfig.yaml", "nozzle.h5", "nozzle.msh", "prepare_info.json", "probe.yaml", "wall_repr.json",
               "wall_physical.csv", "wall_design.csv", "target_axis_M.csv", "MESH_QUALITY.txt"):
        if (mesh_run / fn).is_file():
            shutil.copy2(mesh_run / fn, run / fn)
    cfg = (run / "solverConfig.yaml").read_text()
    i0 = cfg.index("physProp: {"); depth = 0
    for k in range(i0, len(cfg)):
        if cfg[k] == "{": depth += 1
        elif cfg[k] == "}":
            depth -= 1
            if depth == 0:
                i1 = k + 1; break
    cfg = cfg[:i0] + CPG_PHYS + cfg[i1:]
    if variant == "plain":
        cfg = _re.sub(r"^turbulence: \{.*\}$", TURB_PLAIN, cfg, count=1, flags=_re.M)
    cfg = _re.sub(r"nStepOuter: *\d+", f"nStepOuter: {int(main_steps)}", cfg, count=1)
    cfg = _re.sub(r"outStepInterval: *\d+", f"outStepInterval: {int(out_int)}", cfg, count=1)
    if "thermalMethod: 0" not in cfg or ("dilatationCorrection: 0" in cfg) != (variant == "plain"):
        raise SystemExit("solverConfig の書き換えの検査が不成立 — 止める")
    # 物性の検査 (codex diagnose 2026-10-08 Major): 実効の physProp が Sutherland (1.716e-5/273/111) + 定 Pr 0.72 で、
    # 化学種・transport が無いこと。λ(T) = μ(T)·c_p/Pr を SU2 の CONSTANT_PRANDTL と同じ式で照合する
    import yaml as _yaml
    php = _yaml.safe_load(cfg)["physProp"]
    need = {"thermalMethod": 0, "viscMethod": 1, "thermCondMethod": 1, "prandtlLam": 0.72, "cp": CPG_CP, "gamma": CPG_GAMMA}
    bad = {k: php.get(k) for k, v in need.items() if php.get(k) != v}
    if bad or "species" in php or "transport" in php:
        raise SystemExit(f"physProp が CPG の照合の条件と違う ({bad}、species/transport の有無 {'species' in php}/{'transport' in php}) — 止める")
    mu = lambda T: 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)  # noqa: E731  (forge viscMethod 1 = SU2 の指定)
    lam = {T: mu(T) * CPG_CP / 0.72 for T in (300.0, 1500.0)}
    (run / "solverConfig.yaml").write_text(cfg)
    bc = (run / "bcondConfig.yaml").read_text()
    bc2 = _re.sub(r"Y0: *[0-9.eE+-]+, *Y1: *[0-9.eE+-]+, *", "", bc)
    if "Y0" in bc2 or bc2 == bc:
        raise SystemExit("bcond の入口から化学種を外せない — 止める")
    (run / "bcondConfig.yaml").write_text(bc2)
    # 初期値: TP の収束解 → CPG の保存量
    with h5py.File(field_res, "r") as h:
        V = {k: np.array(h["VALUE/" + k], dtype=float) for k in ("ro", "Ux", "Uy", "Uz", "P", "k", "omega")}
    roe = V["P"] / (CPG_GAMMA - 1.0) + 0.5 * V["ro"] * (V["Ux"] ** 2 + V["Uy"] ** 2 + V["Uz"] ** 2)
    new = {"ro": V["ro"], "roUx": V["ro"] * V["Ux"], "roUy": V["ro"] * V["Uy"], "roUz": V["ro"] * V["Uz"], "roe": roe,
           "roK": V["ro"] * V["k"], "roOmega": V["ro"] * V["omega"]}
    with h5py.File(run / "nozzle.h5", "r+") as h:
        for k, v in new.items():
            ds = h["VALUE/" + k]
            if ds.shape != v.shape:
                raise SystemExit(f"{k} の長さが違う ({ds.shape} vs {v.shape})")
            ds[...] = v.astype(ds.dtype)
        for k in ("roY0", "roY1"):
            if "VALUE/" + k in h:
                del h["VALUE/" + k]
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": STAGES, "ramp": None, "ramp_steps": 1000}, cfd_gas="cpg (#26: 設定を書き換え)", turbulence_variant=variant)
    NS.jdump(run / "prepare_info.json", info)
    rec = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #26", "tool": "cold_pair.py prep-cpg", "created": NS.now(),
           "git_head": NS.git_head(), "binary": binrec, "kind": "cpg_" + variant, "mesh_run": mesh_run.name, "field_res": str(field_res),
           "field_res_sha256": NS.sha256_file(field_res), "gas": {"gamma": CPG_GAMMA, "cp": CPG_CP, "R": CPG_CP * (CPG_GAMMA - 1) / CPG_GAMMA,
           "mu": "Sutherland 1.716e-5/273/111", "Pr": 0.72, "Prt": 0.9}, "variant": variant, "stages": STAGES, "main_steps": int(main_steps),
           "out_interval": int(out_int), "bcond_wall": next(l for l in bc2.splitlines() if "physID: 3" in l).strip(),
           "physProp_effective": php, "lambda_W_mK": {str(k): v for k, v in lam.items()},
           "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / RECORD, rec)
    print(f"[prep-cpg] {run.name}: {variant}、格子 {mesh_run.name}、初期値 {field_res.parent.name}/{field_res.name}、{rec['bcond_wall']}")
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep"); p.add_argument("kind"); p.add_argument("run")
    p = sp.add_parser("run"); p.add_argument("run")
    p = sp.add_parser("nan-scan"); p.add_argument("run")
    p = sp.add_parser("prep-ext"); p.add_argument("src"); p.add_argument("run"); p.add_argument("steps", type=int)
    p = sp.add_parser("extract"); p.add_argument("run")
    p = sp.add_parser("extract-b"); p.add_argument("run")
    p = sp.add_parser("theta"); p.add_argument("run"); p.add_argument("--step", type=int, default=None)
    p = sp.add_parser("prep-cpg"); p.add_argument("mesh_run"); p.add_argument("field_res"); p.add_argument("run"); p.add_argument("variant")
    p = sp.add_parser("ab-compare"); p.add_argument("ad"); p.add_argument("tw")
    p = sp.add_parser("gates"); p.add_argument("ad"); p.add_argument("tw")
    p = sp.add_parser("judge"); p.add_argument("ad"); p.add_argument("tw")
    a = ap.parse_args()
    if a.cmd == "prep":
        prep(a.kind, HERE / a.run)
    elif a.cmd == "run":
        sys.exit(run_one(HERE / a.run))
    elif a.cmd == "prep-ext":
        prep_ext(HERE / a.src, HERE / a.run, a.steps)
    elif a.cmd == "extract":
        print(extract(HERE / a.run))
    elif a.cmd == "prep-cpg":
        prep_cpg(HERE / a.mesh_run, HERE / a.field_res, HERE / a.run, a.variant)
    elif a.cmd == "theta":
        print(theta_diag(HERE / a.run, a.step))
    elif a.cmd == "extract-b":
        print(extract_b(HERE / a.run))
    elif a.cmd == "ab-compare":
        ab_compare(a.ad, a.tw)
    elif a.cmd == "gates":
        gates_aws(HERE / a.ad, HERE / a.tw)
    elif a.cmd == "judge":
        judge(a.ad, a.tw)
    elif a.cmd == "nan-scan":
        out = NS.nan_scan(HERE / a.run)
        NS.jdump(HERE / a.run / "NAN_SCAN.json", out)
        print(json.dumps({k: out.get(k) for k in ("first_nonfinite", "first_bad_field")}, ensure_ascii=False))
