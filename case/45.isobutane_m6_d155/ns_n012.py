"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「NS の 3 条件」(N0・N1・N2) の run の準備・検査・実行・記録の時系列。
順序と根拠: plan discretization-moc-axis-limit-and-corrector §6 V5d、諮問 notes/reviews/2026-10-07-euler-grid-switch-plan-diagnose.md。
起動スクリプト run_ns_n012.sh が呼ぶ (forge を直接回さない; 実行は runner_axismach.run_staged_ns → run_case.sh)。

登録 (U4) から決まること:
  - 3 条件の問題は make_ns_n012_problems.py が作る (較正値・k_f・r_t・NS の格子・実効の設定は共通)。
  - IC: 3 条件とも生産の dry NS の最終場 (run_0147 + 延長 run_0149 の通算 80000 step = run_0149 の res_20000.h5) から
    interp_field.py の cross-mesh 移送 (prepare_ns の ic_from = 生産の入口 prep_c2pin.py と同じ経路)。番号写像は使わない。
    移送の検査 (ic-check): ドナーの対応・壁の層の対応・化学種と EOS の基底・新しい壁距離。移送直後の場は定常の比較に使わない。
  - 段階起動 (soft 段から) → 本段は run_0147 と同じ (2 次 convMethod 1・cfl 1・implicitRelax 0.7)、80000 step、5000 ごと出力
    (run_0147 + 0149 と同じ間隔)。判定窓は 60000〜80000 の 5 枚 (評価器 ns_n012_eval.py)。
  - 延長 (主セッションが判断): 本段の最終場 res_80000.h5 から restart_field (同一メッシュ、ビット一致) で 20000 step、soft 段なし。
  - 凝縮 (K): dry の最終場から convert_species_field (run_0148 と同じ手順: cfl 1・implicitRelax なし・18000 step・1000 ごと)。
登録に無く、ここで決めた細部 (plan に書いていない。報告に列挙):
  - 段の構成は stages="full" (soft 1 次 cfl 0.5 3000 → mid 1 次 cfl 1 3000 → 本段)。run_0147 は段なし (番号写像だったため) で、
    cross-mesh の移送から段階起動した生産の実績は run_0116 (同じ "full")。
  - IC の検査の合否: 下の ic_check の「必須」は不成立なら forge を起動しない。「記録」は合否に使わない (壁の層の取り違えは
    壁が動く限り起こるので、場所と量を記録する)。
  - 新しい壁距離: 宛先の nozzle.msh を同じ変換器で変換し直した wall_dist と完全一致 (または相対 1e-6 以内) を必須にする。

usage (case dir、引数の run は case dir からの相対):
  python3 ns_n012.py prep-dry  <N0|N1|N2> <run> --md-offset V [--ic-src-run DIR (乾式のみ)] [--dry]
  python3 ns_n012.py ic-check  <run> [--dry]
  python3 ns_n012.py verify-set <run_N0> <run_N1> <run_N2> --md-offset V [--dry] [--out JSON]
  python3 ns_n012.py prep-cond <N0|N1|N2> <src_dry_run> <cond_run> --md-offset V [--dry]
  python3 ns_n012.py prep-ext  <src_dry_run> <ext_run> [--dry]
  python3 ns_n012.py run       <run>
  python3 ns_n012.py nan-scan  <run>
  python3 ns_n012.py record-series <run> [--euler RUN] [--out CSV] [--nproc 3]
--dry: ローカルの乾式確認。FORGE_BIN が存在しない道を指し FORGE_ALLOW_UNVERIFIED_SPECIES=1 であること (run_ns_n012.sh DRY=1 が設定)。
  乾式の準備は prepare_info.json と NS_N012.json に DRY の印が付き、run が拒否する。乾式の prep-cond は convert_species_field を
  回さない (両 run の解決済みの熱物性 = forge --resolve-species の記録が要るため。本番だけ)。
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TOOLS = REPO / "solver_density_cuda/tools"
DESIGN = REPO / "design"
sys.path.insert(0, str(HERE))
import make_ns_n012_problems as MK  # noqa: E402

PLAN = MK.PLAN_U4
IC_SRC_RUN = "run_0149_ns_mono_final_ext"      # 生産の dry NS の最終場 (run_0147 + 延長 run_0149、通算 80000 step)
IC_SRC_RES = "res_20000.h5"                    # run_0149 の最終 res (= 通算 80000)
REF_PROV_RUN = "run_0147_ns_mono_final"        # forge の sha256 と solverConfig の比較元
SOLVE_JSON = "c2pin_solve_recal.json"          # k_f (C2 方式の較正; 生産の YAML の deltastar_initializer と同じ値であること)
MAIN_STEPS, OUT_INT = 80000, 5000
EXT_STEPS = 20000
COND_STEPS, COND_OUT = 18000, 1000
STAGES = "full"
CFL_MAIN, RELAX = 1.0, 0.7
RECORD = "NS_N012.json"
IC_CHECK = "IC_CHECK.json"
NAN_SCAN = "NAN_SCAN.json"
SET_CHECK = "_band_ab/ns_n012_prep_check.json"
RUN_RE = re.compile(r"^run_\d{4}_[A-Za-z0-9_.-]+$")
CONS_RE = re.compile(r"^(ro|roU[xyz]|roe|roK|roOmega|roGamma|roReth|roXi|roY\d+|rog_.+|roQ[012]_.+)$")
NEAR_WALL_LAYERS = 10                          # 壁の層の対応を記録する層の数 (壁 = 層 0)
WALLDIST_RTOL = 1e-6                           # 変換し直した wall_dist との一致の許容 (相対; 0 が期待値)
NSTEP = ("time", "last", "nStepOuter")
OUTINT = ("time", "outStepInterval")
CFL = ("time", "deltaT", "cfl")
CFLP = ("time", "deltaT", "cfl_pseudo")
RELAXP = ("time", "deltaT", "implicitRelax")
CONVP = ("space", "convMethod")
# 延長の子へ複製する親の入力 (run_mono_ns_ext.sh と同じ集合。出力・判定ファイル・IC の記録は持ち込まない)
EXT_COPY = ("solverConfig.yaml", "bcondConfig.yaml", "nozzle.h5", "prepare_info.json", "probe.yaml", "species_meta.yaml",
            "MESH_QUALITY.txt", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "delta_r_initial.csv",
            "delta_r_initial.json", "nozzle.msh", RECORD)


# --- 小道具 ---------------------------------------------------------------------------------------------------------
def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def git_head() -> str:
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or "unknown"
    except OSError:
        return "unknown"


def res_files(run: Path) -> list:
    """res_<n>.h5 (n の昇順、res_0 を含む)。"""
    fs = [f for f in Path(run).glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)]
    return sorted(fs, key=lambda f: int(re.findall(r"\d+", f.name)[0]))


def step_of(f) -> int:
    return int(re.findall(r"\d+", Path(f).name)[0])


def jload(p) -> dict:
    return json.loads(Path(p).read_text())


def jdump(p, d) -> None:
    Path(p).write_text(json.dumps(d, indent=1, ensure_ascii=False, default=str))


def yaml_strict():
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import yaml_strict as ys
    return ys


def cfg_get(doc, path):
    v = doc
    for k in path:
        if not isinstance(v, dict) or k not in v:
            return None
        v = v[k]
    return v


def check_dry_env(dry: bool) -> None:
    """乾式では forge を起動しない (FORGE_BIN が存在しない道)。本番では化学種の未検証の許可を立てない。"""
    fb = os.environ.get("FORGE_BIN", "")
    if dry:
        if fb and Path(fb).exists():
            raise SystemExit(f"--dry なのに FORGE_BIN={fb} が存在する (乾式では forge を起動しない。run_ns_n012.sh DRY=1 を使う)")
        if os.environ.get("FORGE_ALLOW_UNVERIFIED_SPECIES") != "1":
            raise SystemExit("--dry には FORGE_ALLOW_UNVERIFIED_SPECIES=1 が要る (forge --resolve-species を呼ばないため)")
    else:
        if os.environ.get("FORGE_ALLOW_UNVERIFIED_SPECIES"):
            raise SystemExit("本番で FORGE_ALLOW_UNVERIFIED_SPECIES が立っている — 化学種の属性を検証せずに移すことになる。止める")
        if not fb or not Path(fb).is_file():
            raise SystemExit(f"FORGE_BIN ({fb or '未設定'}) が無い — 止める")


def runner():
    if str(DESIGN) not in sys.path:
        sys.path.insert(0, str(DESIGN))
    from forge_design.evaluate import runner_axismach as RA
    return RA


def k_f_checked(problem: Path) -> float:
    """k_f: c2pin_solve_recal.json と問題 YAML の deltastar_initializer.cf_scale が同じ float (生産の入口と同じ値)。"""
    kf = float(jload(HERE / SOLVE_JSON)["k_f"])
    ys = yaml_strict()
    y = ys.load(problem.read_text())
    cf = (y.get("deltastar_initializer") or {}).get("cf_scale")
    if not isinstance(cf, float) or cf != kf:
        raise SystemExit(f"k_f が一致しない: {SOLVE_JSON} {kf!r} / {problem.name} deltastar_initializer.cf_scale {cf!r} — 止める")
    return kf


def problems_ok(md_offset: str) -> dict:
    """make_ns_n012_problems.py --check と同じ照合 (6 本の中身と記録が --md-offset の値・今の元の YAML と一致)。"""
    rc = MK.main(["--md-offset", md_offset, "--check", "--out-dir", str(HERE)])
    if rc != 0:
        raise SystemExit(f"問題 YAML が --md-offset {md_offset} と一致しない — make_ns_n012_problems.py を回し直すか較正値を確かめる。止める")
    return jload(HERE / MK.RECORD)


def initializer(kf: float) -> dict:
    """prep_c2pin.py と同じ (生産の入口)。"""
    return {"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}


def set_config(run: Path, updates: dict, want: dict) -> None:
    """solverConfig.yaml の値トークンだけを置換し (yaml_strict.replace_scalars)、読み直して want の値を確かめる。"""
    ys = yaml_strict()
    p = run / "solverConfig.yaml"
    txt = ys.replace_scalars(p.read_text(), {k: str(v) for k, v in updates.items()})
    d = ys.load(txt)
    bad = [f"{'.'.join(k)} = {cfg_get(d, k)!r} ({v!r} であること)" for k, v in want.items() if cfg_get(d, k) != v]
    if bad:
        raise SystemExit(f"{p}: 設定の検査が不成立 — 止める: " + "; ".join(bad))
    p.write_text(txt)


def mesh_quality_ok(run: Path) -> str:
    t = (run / "MESH_QUALITY.txt").read_text() if (run / "MESH_QUALITY.txt").is_file() else ""
    line = next((l.strip() for l in t.splitlines()[::-1] if "VERDICT" in l), None)
    if not line or "PASS" not in line:
        raise SystemExit(f"{run}/MESH_QUALITY.txt の VERDICT が PASS でない ({line}) — 止める (U4 の投入の条件)")
    return line


def _prepare_ns(problem: Path, run: Path, **kw) -> dict:
    RA = runner()
    try:
        return RA.prepare_ns(problem, run, **kw)
    except subprocess.CalledProcessError as e:
        print(f"[ns_n012] 子プロセスが失敗: {e.cmd}\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}")
        raise SystemExit(f"prepare_ns が失敗 — 止める (作りかけの {run} は消さない)")


# --- prep-dry ------------------------------------------------------------------------------------------------------
def prep_dry(cond: str, run: Path, md_offset: str, ic_src_run: str | None, dry: bool) -> dict:
    check_dry_env(dry)
    if cond not in MK.CONDS:
        raise SystemExit(f"条件 {cond!r} は {MK.CONDS} のどれか")
    if not RUN_RE.match(run.name):
        raise SystemExit(f"run 名 {run.name!r} が run_NNNN_<slug> でない")
    if run.exists():
        raise SystemExit(f"{run} が既にある — 止める (既存 run は消さない)")
    if ic_src_run and not dry:
        raise SystemExit("--ic-src-run は乾式確認だけ (本番の IC は登録どおり run_0149 の最終場)")
    rec_p = problems_ok(md_offset)
    value = MK.parse_offset(md_offset)
    problem = HERE / MK.out_name(cond, "dry")
    kf = k_f_checked(problem)
    src_run = HERE / (ic_src_run or IC_SRC_RUN)
    rs = res_files(src_run)
    if not rs:
        raise SystemExit(f"IC のドナー {src_run} に res が無い — 止める")
    if not dry and rs[-1].name != IC_SRC_RES:
        raise SystemExit(f"IC のドナー {src_run.name} の最後の res が {rs[-1].name} ({IC_SRC_RES} であること) — 止める")
    src_res = rs[-1]
    src_sha = sha256_file(src_res)
    info = _prepare_ns(problem, run, nsteps=12000, ic_from=src_run, initializer=initializer(kf), cfl_main=5.0, implicit_relax=RELAX)
    # 生産の投入スクリプト (run_mono_ns_chain.sh) と同じ書き換え: cfl 1・本段 80000・5000 ごと (値トークンだけ置換)
    set_config(run, {CFL: "1.0", CFLP: "1.0", NSTEP: MAIN_STEPS, OUTINT: OUT_INT},
               {CFL: CFL_MAIN, CFLP: CFL_MAIN, NSTEP: MAIN_STEPS, OUTINT: OUT_INT, RELAXP: RELAX, CONVP: 1})
    info["stages"] = {"stages": STAGES, "ramp": None, "ramp_steps": 1000}
    if dry:
        info["DRY"] = True
    jdump(run / "prepare_info.json", info)
    mq = mesh_quality_ok(run)
    rec = {"plan": PLAN, "tool": "ns_n012.py prep-dry", "created": now(), "git_head": git_head(), "dry": bool(dry),
           "role": "dry", "condition": cond, "md_offset": value, "md_offset_token": rec_p["md_offset_token"],
           "problem": problem.name, "problem_sha256": sha256_file(problem), "problems_record_sha256": sha256_file(HERE / MK.RECORD),
           "k_f": kf, "solve_json": SOLVE_JSON,
           "ic": {"src_run": src_run.name, "src_dir": os.path.relpath(src_run, HERE), "src_res": src_res.name, "src_res_sha256": src_sha,
                  "via": "prepare_ns(ic_from) → interp_field.py",
                  "registered": (src_run.name == IC_SRC_RUN and src_res.name == IC_SRC_RES)},
           "stages": STAGES, "main_steps": MAIN_STEPS, "out_interval": OUT_INT, "cfl_main": CFL_MAIN, "implicit_relax": RELAX,
           "mesh_quality": mq, "effective": {"pw_upstream": info.get("pw_upstream"), "moc": info.get("moc"),
                                             "Md_moc_offset": info.get("Md_moc_offset")},
           "nozzle_sha256_after_prep": sha256_file(run / "nozzle.h5")}
    jdump(run / RECORD, rec)
    print(f"[prep-dry] {cond} {run.name}: 問題 {problem.name}、IC {src_run.name}/{src_res.name}、{mq}")
    return rec


# --- ic-check -------------------------------------------------------------------------------------------------------
def interp_like_fields(src_res: Path) -> dict:
    """interp_field.py の res 経路と同じ式で転送配列を作る (同じ numpy の演算; ドナーの対応を独立に再現するため)。
    res に保存量 roe が無いとき interp_field は CPG 式で組み直す (TP では誤り) ので、ここでは拒否する。"""
    import h5py
    with h5py.File(src_res, "r") as s:
        V = s["VALUE"]
        if not ("P" in V and "Ux" in V):
            raise ValueError(f"{src_res} は res (原始量を持つ出力) でない")
        if "roe" not in V:
            raise ValueError(f"{src_res} に保存量 roe が無い (interp_field は CPG 式で組み直す — TP の EOS の基底が崩れる)")
        ro = np.array(V["ro"])
        Ux, Uy, Uz = np.array(V["Ux"]), np.array(V["Uy"]), np.array(V["Uz"])
        fields = {"ro": ro, "roUx": ro * Ux, "roUy": ro * Uy, "roUz": ro * Uz, "roe": np.array(V["roe"])}
        if "k" in V and "omega" in V:
            fields["roK"] = ro * np.array(V["k"])
            fields["roOmega"] = ro * np.array(V["omega"])
        if "roGamma" in V and "roReth" in V:
            fields["roGamma"] = np.array(V["roGamma"]); fields["roReth"] = np.array(V["roReth"])
        elif "gammaTr" in V and "reTheta" in V:
            fields["roGamma"] = ro * np.array(V["gammaTr"]); fields["roReth"] = ro * np.array(V["reTheta"])
        for key in V:
            if key.startswith("Y") and key[1:].isdigit():
                fields["ro" + key] = ro * np.array(V[key])
            if key == "roXi":
                fields["roXi"] = np.array(V[key])
            if key.startswith(("g_", "Q0_", "Q1_", "Q2_")):
                fields["ro" + key] = ro * np.array(V[key])
        if "roXi" not in fields and "Xi" in V:
            fields["roXi"] = ro * np.clip(np.array(V["Xi"], dtype=np.float64), 0.0, 1.0)
    return fields


def _structured(h5, n_nodes):
    from ic_index_map import structured_shape
    return structured_shape(h5["VIZMESH/CONNE"], n_nodes)


def wall_layer_record(cd, cs, idx, ni, nj, nis, njs, scale) -> dict:
    """壁の層の対応 (記録)。宛先の各節点の壁からの層 L = nj−1−j と、最近傍のドナーの層を比べる。
    壁の移動 Δr_w (宛先の壁節点の r − 同じ x のドナーの壁の r) と第 1 セル厚 h1 も列ごとに記録する。"""
    jd = np.arange(ni * nj) % nj
    js = idx % njs
    Ld, Ls = nj - 1 - jd, njs - 1 - js
    Xd, Rd = cd[:, 0].reshape(ni, nj), cd[:, 1].reshape(ni, nj)
    Xs, Rs = cs[:, 0].reshape(nis, njs), cs[:, 1].reshape(nis, njs)
    xw_s, rw_s = Xs[:, -1], Rs[:, -1]
    o = np.argsort(xw_s)
    dRw = Rd[:, -1] - np.interp(Xd[:, -1], xw_s[o], rw_s[o])
    h1 = Rd[:, -1] - Rd[:, -2]
    col = np.arange(ni * nj) // nj
    out = {"layers": {}, "n_columns": int(ni), "n_layers": int(nj)}
    for L in range(min(NEAR_WALL_LAYERS, nj)):
        m = Ld == L
        bad = m & (Ls != L)
        r = {"n": int(m.sum()), "n_donor_other_layer": int(bad.sum()),
             "max_abs_layer_diff": int(np.abs(Ls[m] - L).max()) if m.any() else 0}
        if bad.any():
            xb = cd[bad, 0] / scale
            r.update(x_rt_range=[float(xb.min()), float(xb.max())],
                     layer_diff_counts={str(int(v)): int(np.count_nonzero((Ls[bad] - L) == v)) for v in np.unique(Ls[bad] - L)})
        out["layers"][str(L)] = r
    near = Ld < NEAR_WALL_LAYERS
    bad_near = near & (Ls != Ld)
    expl = np.abs(dRw[col]) > 0.5 * h1[col]
    out["near_wall"] = {"n": int(near.sum()), "n_donor_other_layer": int(bad_near.sum()),
                        "n_in_columns_wall_moved_gt_half_h1": int((bad_near & expl).sum()),
                        "n_in_columns_wall_moved_le_half_h1": int((bad_near & ~expl).sum())}
    k = int(np.argmax(np.abs(dRw)))
    out["wall_shift"] = {"max_abs_um": float(np.abs(dRw).max() * 1e6), "x_rt_at_max": float(Xd[k, -1] / scale),
                         "signed_um_at_max": float(dRw[k] * 1e6),
                         "n_columns_gt_half_h1": int(np.count_nonzero(np.abs(dRw) > 0.5 * h1)),
                         "x_rt_range_gt_half_h1": ([float(Xd[np.abs(dRw) > 0.5 * h1, -1].min() / scale),
                                                    float(Xd[np.abs(dRw) > 0.5 * h1, -1].max() / scale)]
                                                   if np.any(np.abs(dRw) > 0.5 * h1) else None),
                         "h1_um_min": float(h1.min() * 1e6), "h1_um_median": float(np.median(h1) * 1e6)}
    if (ni, nj) == (nis, njs):
        dx = (Xd[:, 0] - Xs[:, 0]) / scale
        out["column_x_shift_rt"] = {"max_abs": float(np.abs(dx).max()), "median_abs": float(np.median(np.abs(dx)))}
    out["note"] = ("最近傍は、壁の移動または列の x のずれが壁際のセル厚を超える所で別の層の値を拾う。cross-mesh の移送では避けられないので"
                   "記録だけにし、移送直後の場は定常の比較に使わない (段階起動の後の本段の窓だけで判定する)")
    return out


def reconvert_wall_dist(run: Path):
    """宛先の nozzle.msh を同じ変換器で変換し直した /VALUE/wall_dist と /MESH/COORD (乾式でも変換器は起動してよい)。"""
    import h5py
    RA = runner()
    conv = RA.converter_path()
    with tempfile.TemporaryDirectory(prefix="ns_n012_wd_") as td:
        td = Path(td)
        shutil.copy2(run / "nozzle.msh", td / "nozzle.msh")
        for f in run.glob("*.yaml"):
            shutil.copy2(f, td / f.name)
        r = subprocess.run([str(conv), "nozzle.msh", "wd.h5"], cwd=td, env=RA._ENV, capture_output=True, text=True)
        if r.returncode != 0 or not (td / "wd.h5").is_file():
            raise RuntimeError(f"変換器 {conv} が失敗 (rc {r.returncode}): {(r.stdout + r.stderr)[-1500:]}")
        with h5py.File(td / "wd.h5", "r") as f:
            return np.array(f["VALUE/wall_dist"]), np.array(f["MESH/COORD"]), str(conv)


def ic_check(run: Path, dry: bool) -> dict:
    """移送の検査 (U4 の登録)。必須の不成立があれば VERDICT NG (forge を起動しない)。"""
    import h5py
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import forge_species as fsp
    from ic_index_map import nearest_index, nearest_vs_index
    rec0 = jload(run / RECORD)
    src_run = HERE / rec0["ic"].get("src_dir", rec0["ic"]["src_run"])
    src_res = src_run / rec0["ic"]["src_res"]
    dst = run / "nozzle.h5"
    must, out = [], {"plan": PLAN, "tool": "ns_n012.py ic-check", "date": now(), "run": run.name, "dry": bool(dry),
                     "src_res": f"{src_run.name}/{src_res.name}",
                     "note": "移送直後の場は定常の比較に使わない (評価は本段の窓 60000〜80000 だけ)"}
    # ---- I1 ドナーの対応 (必須) ----
    i1 = {}
    if sha256_file(src_res) != rec0["ic"]["src_res_sha256"]:
        must.append("I1: ドナーの res の sha256 が準備の記録と違う")
    last = res_files(src_run)[-1]
    if last.name != src_res.name:
        must.append(f"I1: ドナーの最後の res が {last.name} (prepare_ns の ic_from が選ぶのは最後の res; 記録は {src_res.name})")
    try:
        fields = interp_like_fields(src_res)
    except ValueError as e:
        must.append(f"I1/I3: {e}")
        fields = {}
    idx, dist, nd = nearest_index(src_res, dst)
    i1["nn_dim"] = int(nd)
    i1["nn_dist_m"] = {"max": float(dist.max()), "p50": float(np.quantile(dist, 0.5)), "p99": float(np.quantile(dist, 0.99)),
                       "n_zero": int(np.count_nonzero(dist == 0)), "n": int(dist.size)}
    mism, floored, moved = [], None, []
    with h5py.File(dst, "r") as d, h5py.File(src_res, "r") as s:
        dnames = sorted(n for n in d["VALUE"] if n != "wall_dist")
        notmoved = [n for n in dnames if n not in fields]
        if notmoved:
            must.append(f"I1: 宛先の量 {notmoved} がドナーから移っていない (等エントロピーの IC のまま)")
        for n, arr in fields.items():
            if n not in d["VALUE"]:
                must.append(f"I1: ドナーの量 {n} が宛先に無い")
                continue
            got = np.array(d["VALUE"][n])
            exp = arr[idx].astype(got.dtype)
            if n == "roOmega":
                ok = np.all(got >= exp)
                floored = int(np.count_nonzero(got != exp))
                if not ok:
                    mism.append(f"{n} (ω の床より小さい値がある)")
            elif not np.array_equal(got, exp):
                mism.append(f"{n} (不一致 {int(np.count_nonzero(got != exp))} 節点)")
            moved.append(n)
        om = float(np.max(np.array(s["VALUE/omega"]))) if "omega" in s["VALUE"] else 0.0
        fin = [n for n in dnames if not np.all(np.isfinite(np.array(d["VALUE"][n])))]
        if fin:
            must.append(f"I1: 宛先に非有限値 {fin}")
        if not np.all(np.array(d["VALUE/ro"]) > 0):
            must.append("I1: 宛先の ro に 0 以下")
        cd = np.array(d["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        n_d = cd.shape[0]
        ni, nj = _structured(d, n_d)
        wd_dst = np.array(d["VALUE/wall_dist"])
    if mism:
        must.append("I1: 宛先の値が interp_field の最近傍の転送 (ドナー = 記録の res) と一致しない: " + "; ".join(mism))
    if om < 1.0:
        must.append(f"I1: ドナーの omega の最大 {om:.3g} < 1 (Euler の場とみなされ k・ω が入口値で上書きされる)")
    i1.update(transferred=moved, roOmega_floored_nodes=floored, donor_omega_max=om)
    out["I1_donor"] = i1
    # ---- I2 壁の層の対応 (記録) ----
    with h5py.File(src_run / "nozzle.h5", "r") as fs:
        cs = np.array(fs["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        nis, njs = _structured(fs, cs.shape[0])
        wd_src = np.array(fs["VALUE/wall_dist"])
    with h5py.File(src_res, "r") as fr:
        if not np.array_equal(np.array(fr["MESH/COORD"]).reshape(-1, 3).astype(np.float64), cs):
            must.append("I2: ドナーの res の座標がドナーの nozzle.h5 と違う")
    scale = float(jload(run / "prepare_info.json")["scale_m"])
    i2 = wall_layer_record(cd, cs, idx, ni, nj, nis, njs, scale)
    if (ni, nj) == (nis, njs):
        names = [n for n in ("ro", "roUx", "roe") if n in fields]
        sv = {n: fields[n] for n in names}
        cmp_ = nearest_vs_index(idx, dist, cd, ni, nj, scale, sv, names, list_max=50)
        i2["nearest_vs_same_index"] = {k: cmp_[k] for k in ("n_mismatch", "match_rate", "dj_counts", "x_range_rt", "j_range",
                                                            "conserved_diff_nearest_minus_index")}
    out["I2_wall_layers"] = i2
    # ---- I3 化学種と EOS の基底 (必須) ----
    i3 = {}
    ys = yaml_strict()
    cfg_s, cfg_d = ys.load((src_run / "solverConfig.yaml").read_text()), ys.load((run / "solverConfig.yaml").read_text())
    for key in ("physProp", "turbulence"):
        if (cfg_s.get(key) or {}) != (cfg_d.get(key) or {}):
            must.append(f"I3: solverConfig.yaml の {key} がドナーと違う")
    for f in ("species_meta.yaml",):
        a_, b_ = src_run / f, run / f
        if a_.is_file() != b_.is_file() or (a_.is_file() and sha256_file(a_) != sha256_file(b_)):
            must.append(f"I3: {f} がドナーと違う")
    try:
        sig_s, sig_d = fsp.species_signature(str(src_run)), fsp.species_signature(str(run))
        diffs = fsp.compare_signatures(sig_s, sig_d)
        hard = [x for x in diffs if "unverifiable" not in x]
        if hard:
            must.append(f"I3: 化学種の署名が違う: {hard}")
        i3.update(names=sig_s["names"], thermoHrefTemp=sig_s["thermoHrefTemp"], thermalMethod=sig_s["thermalMethod"],
                  unverifiable=[x for x in diffs if "unverifiable" in x])
    except Exception as e:  # noqa: BLE001
        must.append(f"I3: 化学種の署名を作れない: {type(e).__name__}: {e}")
    at_s, at_d = fsp.field_species_attrs(str(src_res)), fsp.field_species_attrs(str(dst))
    i3.update(src_attrs=at_s, dst_attrs=at_d)
    if not dry:
        if not at_d["species_hash"] or at_d["species_hash"] != at_s["species_hash"] or at_d["species_input_unverified"] not in (None, 0):
            must.append(f"I3: 宛先の化学種の属性がドナーから継承されていない (宛先 {at_d}, ドナー {at_s})")
    else:
        i3["dry_note"] = "乾式: forge --resolve-species を呼ばないので属性は付かない (本番では継承を必須にする)"
    with h5py.File(src_res, "r") as s:
        if "h0" in s["VALUE"]:
            a = s["VALUE/h0"].attrs.get("h0_includes_k")
            i3["src_h0_includes_k"] = None if a is None else int(a)
    out["I3_species_eos"] = i3
    # ---- I4 新しい壁距離 (必須) ----
    i4 = {}
    try:
        wd_new, coord_new, conv = reconvert_wall_dist(run)
        with h5py.File(dst, "r") as d:
            same_coord = np.array_equal(coord_new, np.array(d["MESH/COORD"]))
        if not same_coord:
            must.append("I4: 変換し直した座標が宛先の nozzle.h5 と違う")
        exact = np.array_equal(wd_new, wd_dst)
        rel = float(np.max(np.abs(wd_new.astype(np.float64) - wd_dst) / np.maximum(np.abs(wd_new.astype(np.float64)), 1e-30)))
        i4.update(converter=conv, equal_to_reconverted=bool(exact), max_rel_diff_vs_reconverted=rel)
        if not exact and not rel <= WALLDIST_RTOL:
            must.append(f"I4: 宛先の wall_dist が変換し直した値と違う (最大相対差 {rel:.3e} > {WALLDIST_RTOL:g})")
    except Exception as e:  # noqa: BLE001
        must.append(f"I4: 宛先の nozzle.msh を変換し直せない: {type(e).__name__}: {e}")
    dmapped = wd_src[idx].astype(np.float64)
    ad = np.abs(wd_dst.astype(np.float64) - dmapped)
    pos = wd_dst > 0
    rel = ad[pos] / wd_dst[pos].astype(np.float64)
    # 記録: 宛先の wall_dist がドナーの値 (最近傍で写したもの) と違うこと = 移送されていないことの傍証 (壁が動いた・列がずれた所で違う)
    i4["vs_donor_mapped"] = {"n_equal": int(np.count_nonzero(wd_dst == wd_src[idx].astype(wd_dst.dtype))), "n": int(wd_dst.size),
                             "max_abs_diff_um": float(ad.max() * 1e6), "p99_abs_diff_um": float(np.quantile(ad, 0.99) * 1e6),
                             "p99_rel_diff_wd_gt0": (float(np.quantile(rel, 0.99)) if rel.size else None),
                             "n_wall_dist_zero": int(np.count_nonzero(~pos))}
    out["I4_wall_dist"] = i4
    out["must_failures"] = must
    out["VERDICT"] = "OK" if not must else "NG"
    out["nozzle_sha256"] = sha256_file(dst)
    jdump(run / IC_CHECK, out)
    print(f"[ic-check] {run.name}: {out['VERDICT']}  最近傍距離の最大 {i1['nn_dist_m']['max']:.3e} m、"
          f"壁際 {NEAR_WALL_LAYERS} 層で別の層を拾った節点 {i2['near_wall']['n_donor_other_layer']} / {i2['near_wall']['n']} "
          f"(壁の移動の最大 {i2['wall_shift']['max_abs_um']:.3f} µm、x/r_t {i2['wall_shift']['x_rt_at_max']:.3f})")
    for m in must:
        print("  NG: " + m)
    return out


# --- verify-set (3 条件の固定) ----------------------------------------------------------------------------------------
def verify_set(runs: list, md_offset: str, dry: bool, out_path: Path | None = None) -> dict:
    """3 条件で較正値・k_f・r_t・NS の格子・実効の設定・IC・段が同じで、変えたのが条件のキーだけかを調べる。"""
    ys = yaml_strict()
    value = MK.parse_offset(md_offset)
    why, rows = [], {}
    recs = {}
    for cond, run in zip(MK.CONDS, runs):
        run = Path(run)
        r = jload(run / RECORD)
        info = jload(run / "prepare_info.json")
        recs[cond] = (run, r, info)
        if r.get("condition") != cond or r.get("role") != "dry":
            why.append(f"{run.name}: 記録の条件 {r.get('condition')}/{r.get('role')} が {cond}/dry でない")
        if bool(r.get("dry")) != bool(dry) or bool(info.get("DRY")) != bool(dry):
            why.append(f"{run.name}: 乾式の印が要求 ({dry}) と違う")
        ic = jload(run / IC_CHECK) if (run / IC_CHECK).is_file() else {}
        if ic.get("VERDICT") != "OK":
            why.append(f"{run.name}: IC の検査が OK でない ({ic.get('VERDICT')})")
        sp = MK.SPEC[cond]
        pwu = (info.get("pw_upstream") or {}).get("value")
        moc = info.get("moc") or {}
        eff = {"pw_upstream": pwu, "moc": moc, "Md_moc_offset": info.get("Md_moc_offset"), "scale_m": info.get("scale_m"),
               "cf_scale": (info.get("initializer") or {}).get("cf_scale")}
        eff["moc"] = {k: moc.get(k) for k in ("axis_limit", "corrector", "pairs", "gate")}
        rows[cond] = {"run": run.name, **eff, "ic": r.get("ic"), "mesh": info.get("mesh")}
        if pwu != sp["pw_upstream"]:
            why.append(f"{run.name}: 実効の pw_upstream {pwu!r} ({sp['pw_upstream']!r} であること)")
        for k, v in (("axis_limit", sp["moc_axis_limit"]), ("corrector", sp["moc_corrector"])):
            if moc.get(k) != v:
                why.append(f"{run.name}: prepare_info の moc.{k} = {moc.get(k)!r} ({v!r} であること)")
        if info.get("Md_moc_offset") != value:
            why.append(f"{run.name}: 実効の Md_moc_offset {info.get('Md_moc_offset')!r} ({value!r} であること)")
    base_run, base_r, base_info = recs["N0"]
    for cond in ("N1", "N2"):
        run, r, info = recs[cond]
        for k in ("scale_m",):
            if info.get(k) != base_info.get(k):
                why.append(f"{run.name}: {k} が N0 と違う ({info.get(k)!r} / {base_info.get(k)!r})")
        if (info.get("initializer") or {}).get("cf_scale") != (base_info.get("initializer") or {}).get("cf_scale"):
            why.append(f"{run.name}: k_f が N0 と違う")
        mi, mb = dict(info.get("mesh") or {}), dict(base_info.get("mesh") or {})
        for m in (mi, mb):
            m.pop("hashes", None)            # 座標は壁に依るので条件で違う (ハッシュは比べない)
        if mi != mb:
            why.append(f"{run.name}: 格子の設定 (prepare_info.mesh) が N0 と違う")
        for f in ("solverConfig.yaml", "bcondConfig.yaml", "species_meta.yaml", "probe.yaml"):
            if (run / f).read_bytes() != (base_run / f).read_bytes():
                why.append(f"{run.name}: {f} が N0 と違う (実効の設定を固定)")
        if r.get("ic") != base_r.get("ic"):
            why.append(f"{run.name}: IC のドナーが N0 と違う")
        for k in ("stages", "main_steps", "out_interval", "cfl_main", "implicit_relax", "k_f", "md_offset"):
            if r.get(k) != base_r.get(k):
                why.append(f"{run.name}: 記録の {k} が N0 と違う")
    # 生産 (run_0147) の設定との差は nStepOuter だけ (本番のみ; 乾式は比較元が無ければ省く)
    ref = HERE / REF_PROV_RUN
    cmp_ref = None
    if ref.is_dir():
        a = ys.load((ref / "solverConfig.yaml").read_text())
        b = ys.load((base_run / "solverConfig.yaml").read_text())
        dd = {".".join(map(str, p)): v for p, v in MK.diff_paths(a, b).items()}
        cmp_ref = {k: [repr(x) for x in v] for k, v in dd.items()}
        if set(dd) != {".".join(NSTEP)}:
            why.append(f"solverConfig.yaml の {REF_PROV_RUN} との差が nStepOuter だけでない: {sorted(dd)}")
        if ys.load((ref / "bcondConfig.yaml").read_text()) != ys.load((base_run / "bcondConfig.yaml").read_text()):
            why.append(f"bcondConfig.yaml が {REF_PROV_RUN} と違う")
    elif not dry:
        why.append(f"比較元 {REF_PROV_RUN} が無い")
    out = {"plan": PLAN, "tool": "ns_n012.py verify-set", "date": now(), "dry": bool(dry), "md_offset": value,
           "runs": rows, "solverConfig_vs_" + REF_PROV_RUN: cmp_ref, "failures": why, "VERDICT": "OK" if not why else "NG"}
    p = out_path or (HERE / (SET_CHECK.replace(".json", "_dry.json") if dry else SET_CHECK))
    p.parent.mkdir(parents=True, exist_ok=True)
    jdump(p, out)
    print(f"[verify-set] {out['VERDICT']} → {p.relative_to(HERE) if p.is_relative_to(HERE) else p}")
    for w in why:
        print("  NG: " + w)
    return out


# --- prep-cond ------------------------------------------------------------------------------------------------------
def same_mesh(a: Path, b: Path) -> list:
    """2 つの入力 h5 の座標・接続・境界が同一か (index コピーの前提)。"""
    import h5py
    from ic_index_map import TOPOLOGY, BCOND_TOPO
    why = []
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        for n in ("MESH/COORD",) + TOPOLOGY:
            if (n in fa) != (n in fb) or (n in fa and not np.array_equal(np.array(fa[n]), np.array(fb[n]))):
                why.append(f"{n} が違う")
        if sorted(fa["BCONDS"]) != sorted(fb["BCONDS"]):
            why.append("BCONDS の physID が違う")
        else:
            for g in fa["BCONDS"]:
                for k in BCOND_TOPO:
                    if (k in fa["BCONDS"][g]) != (k in fb["BCONDS"][g]) or (
                            k in fa["BCONDS"][g] and not np.array_equal(np.array(fa["BCONDS"][g][k]), np.array(fb["BCONDS"][g][k]))):
                        why.append(f"BCONDS/{g}/{k} が違う")
    return why


def prep_cond(cond: str, src_run: Path, run: Path, md_offset: str, dry: bool) -> dict:
    """run_mono_ns_chain.sh の K と同じ手順 (prepare_ns: cfl 1・implicitRelax なし・18000 step、IC は src の最終場を
    convert_species_field で)。src は同じ条件の dry (本段か延長) の run。乾式では src の nozzle.h5 (移送後の IC) を使う。"""
    check_dry_env(dry)
    if not RUN_RE.match(run.name):
        raise SystemExit(f"run 名 {run.name!r} が run_NNNN_<slug> でない")
    if run.exists():
        raise SystemExit(f"{run} が既にある — 止める (既存 run は消さない)")
    rec_p = problems_ok(md_offset)
    value = MK.parse_offset(md_offset)
    srec = jload(src_run / RECORD)
    if srec.get("condition") != cond or srec.get("role") not in ("dry", "ext") or srec.get("md_offset") != value:
        raise SystemExit(f"src {src_run.name} の記録 ({srec.get('condition')}/{srec.get('role')}/{srec.get('md_offset')}) が "
                         f"{cond}/dry|ext/{value} でない — 止める")
    if bool(srec.get("dry")) != bool(dry):
        raise SystemExit(f"src {src_run.name} の乾式の印 {srec.get('dry')} が要求 ({dry}) と違う — 止める")
    if dry:
        src_h5 = src_run / "nozzle.h5"
    else:
        rs = res_files(src_run)
        want = MAIN_STEPS if srec["role"] == "dry" else EXT_STEPS
        if not rs or step_of(rs[-1]) != want:
            raise SystemExit(f"src {src_run.name} の最後の res が {rs[-1].name if rs else None} (res_{want}.h5 であること) — 止める")
        src_h5 = rs[-1]
    problem = HERE / MK.out_name(cond, "cond")
    kf = k_f_checked(problem)
    if kf != srec.get("k_f"):
        raise SystemExit("k_f が dry の記録と違う — 止める")
    info = _prepare_ns(problem, run, nsteps=COND_STEPS, ic_from=None, initializer=initializer(kf), cfl_main=1.0, implicit_relax=None)
    info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}
    info["restart_from"] = f"{src_run.name}/{src_h5.name} (convert_species_field)"
    if dry:
        info["DRY"] = True
    jdump(run / "prepare_info.json", info)
    ys = yaml_strict()
    c = ys.load((run / "solverConfig.yaml").read_text())
    bad = [f"{'.'.join(k)} = {cfg_get(c, k)!r} ({v!r} であること)" for k, v in
           ((NSTEP, COND_STEPS), (OUTINT, COND_OUT), (CFL, 1.0), (CFLP, 1.0), (CONVP, 1), (RELAXP, None)) if cfg_get(c, k) != v]
    if "condensation" not in c:
        bad.append("solverConfig.yaml に condensation が無い")
    if bad:
        raise SystemExit(f"{run}/solverConfig.yaml の検査が不成立 (run_0148 と同じ手順でない) — 止める: " + "; ".join(bad))
    mq = mesh_quality_ok(run)
    mesh_why = same_mesh(src_run / "nozzle.h5", run / "nozzle.h5")
    if mesh_why:
        raise SystemExit(f"凝縮の格子が dry ({src_run.name}) と違う (convert_species_field は index コピー) — 止める: " + "; ".join(mesh_why))
    log = run / "convert_species_field.log"
    cmd = [sys.executable, str(TOOLS / "convert_species_field.py"), str(src_h5), str(run / "nozzle.h5"),
           "--meta", str(run / "species_meta.yaml"), "--src-run", str(src_run), "--dst-run", str(run)]
    if dry:
        # 乾式: 変換器は両 run の解決済みの熱物性 (forge --resolve-species の記録、凝縮では液相つき) を要るので回さない (本番だけ)
        tail = ["乾式: convert_species_field は回していない (forge --resolve-species の記録が要る)"]
        log.write_text(tail[0] + "\n" + " ".join(cmd) + "\n")
    else:
        r = subprocess.run(cmd, capture_output=True, text=True, env=runner()._ENV)
        log.write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            print((r.stdout + r.stderr)[-3000:])
            raise SystemExit(f"convert_species_field が失敗 (rc {r.returncode}) — 止める ({log})")
        tail = (r.stdout + r.stderr).strip().splitlines()[-1:]
    rec = {"plan": PLAN, "tool": "ns_n012.py prep-cond", "created": now(), "git_head": git_head(), "dry": bool(dry),
           "role": "cond", "condition": cond, "md_offset": value, "md_offset_token": rec_p["md_offset_token"],
           "problem": problem.name, "problem_sha256": sha256_file(problem), "k_f": kf,
           "parent": src_run.name, "parent_role": srec["role"], "src": src_h5.name, "src_sha256": sha256_file(src_h5),
           "stages": "none", "steps": COND_STEPS, "out_interval": COND_OUT, "mesh_quality": mq,
           "convert_species_field_tail": tail,
           "nozzle_sha256_after_prep": sha256_file(run / "nozzle.h5")}
    jdump(run / RECORD, rec)
    print(f"[prep-cond] {cond} {run.name}: 問題 {problem.name}、IC {src_run.name}/{src_h5.name} (convert_species_field)")
    return rec


# --- prep-ext -------------------------------------------------------------------------------------------------------
def prep_ext(src_run: Path, run: Path, dry: bool) -> dict:
    """延長 1 回 (登録「未達の量があれば 20000 step」): src (本段) の res_80000.h5 から restart_field (同一メッシュ、ビット一致) で
    新しい run に継ぎ、本段と同じ設定で 20000 step (5000 ごと)。soft 段なし。乾式では src の nozzle.h5 から継ぐ (検査の通しだけ)。"""
    check_dry_env(dry)
    if not RUN_RE.match(run.name):
        raise SystemExit(f"run 名 {run.name!r} が run_NNNN_<slug> でない")
    if run.exists():
        raise SystemExit(f"{run} が既にある — 止める (既存 run は消さない)")
    srec = jload(src_run / RECORD)
    if srec.get("role") != "dry" or bool(srec.get("dry")) != bool(dry):
        raise SystemExit(f"src {src_run.name} は本段の dry の run でない ({srec.get('role')}, 乾式 {srec.get('dry')}) — 止める")
    if dry:
        src_h5 = src_run / "nozzle.h5"
    else:
        rs = res_files(src_run)
        if not rs or rs[-1].name != f"res_{MAIN_STEPS}.h5":
            raise SystemExit(f"src {src_run.name} の最後の res が {rs[-1].name if rs else None} (res_{MAIN_STEPS}.h5 であること) — 止める")
        src_h5 = rs[-1]
    ys = yaml_strict()
    ptext = (src_run / "solverConfig.yaml").read_text()
    P = ys.load(ptext)
    if cfg_get(P, NSTEP) != MAIN_STEPS or cfg_get(P, OUTINT) != OUT_INT:
        raise SystemExit(f"src の solverConfig が本段 (nStepOuter {MAIN_STEPS}・{OUT_INT} ごと) でない — 止める")
    ctext = ys.replace_scalars(ptext, {NSTEP: str(EXT_STEPS)})
    dd = MK.diff_paths(P, ys.load(ctext))
    if set(dd) != {NSTEP}:
        raise SystemExit(f"延長の solverConfig の差が nStepOuter だけでない ({sorted(dd)}) — 止める")
    run.mkdir(parents=True)
    copied = []
    for f in EXT_COPY:
        if (src_run / f).is_file():
            shutil.copy2(src_run / f, run / f)
            copied.append(f)
    for p in sorted(src_run.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
        copied.append(p.name)
    (run / RECORD).unlink(missing_ok=True)
    (run / "solverConfig.yaml").write_text(ctext)
    log = run / "restart_field.log"
    cmd = [sys.executable, str(TOOLS / "restart_field.py"), str(src_h5), str(run / "nozzle.h5"), "--dst-run", str(run)]
    r = subprocess.run(cmd, capture_output=True, text=True, env=runner()._ENV)
    log.write_text(r.stdout + r.stderr)
    if r.returncode != 0 or "ビット一致" not in (r.stdout + r.stderr):
        print((r.stdout + r.stderr)[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める ({log}; 作りかけの {run} は消さない)")
    info = jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src_run.name,
                restart_from=f"{src_run.name}/{src_h5.name} (restart_field)")
    if dry:
        info["DRY"] = True
    jdump(run / "prepare_info.json", info)
    rec = {**{k: srec[k] for k in ("plan", "condition", "md_offset", "md_offset_token", "problem", "problem_sha256", "k_f",
                                   "cfl_main", "implicit_relax", "out_interval")},
           "tool": "ns_n012.py prep-ext", "created": now(), "git_head": git_head(), "dry": bool(dry), "role": "ext",
           "parent": src_run.name, "parent_end": MAIN_STEPS, "ext_steps": EXT_STEPS, "stages": "none",
           "src": src_h5.name, "src_sha256": sha256_file(src_h5), "copied": copied,
           "restart_field_tail": (r.stdout + r.stderr).strip().splitlines()[-1:],
           "nozzle_sha256_after_prep": sha256_file(run / "nozzle.h5")}
    jdump(run / RECORD, rec)
    print(f"[prep-ext] {run.name} ← {src_run.name}/{src_h5.name}: {EXT_STEPS} step、{rec['restart_field_tail']}")
    return rec


# --- run ------------------------------------------------------------------------------------------------------------
def run_one(run: Path) -> int:
    """forge を回す (runner_axismach.run_staged_ns → run_case.sh)。乾式の準備・検査が OK でない準備・準備の後に変わった
    nozzle.h5・既に出力がある run は回さない。"""
    info = jload(run / "prepare_info.json")
    rec = jload(run / RECORD)
    if info.get("DRY") or rec.get("dry"):
        raise SystemExit(f"{run} は乾式確認の準備 — 回さない")
    check_dry_env(False)
    if any((run / f).exists() for f in ("stage_manifest.json", "residual_history.csv")) or res_files(run):
        raise SystemExit(f"{run}: 既に出力がある — 回さない (既存 run は消さない)")
    if sha256_file(run / "nozzle.h5") != rec.get("nozzle_sha256_after_prep"):
        raise SystemExit(f"{run}: nozzle.h5 が準備の後に変わった — 回さない")
    if rec["role"] == "dry":
        ic = jload(run / IC_CHECK) if (run / IC_CHECK).is_file() else {}
        if ic.get("VERDICT") != "OK" or ic.get("nozzle_sha256") != rec.get("nozzle_sha256_after_prep"):
            raise SystemExit(f"{run}: IC の検査 ({IC_CHECK}) が OK でない・検査の後に nozzle.h5 が変わった — 回さない")
        sc = HERE / SET_CHECK
        s = jload(sc) if sc.is_file() else {}
        if s.get("VERDICT") != "OK" or run.name not in [v.get("run") for v in (s.get("runs") or {}).values()]:
            raise SystemExit(f"{run}: 3 条件の固定の検査 ({SET_CHECK}) が OK でない・この run を含まない — 回さない")
    RA = runner()
    stages = rec.get("stages", "none")
    rc = RA.run_staged_ns(run, stages=stages)
    last = res_files(run)
    print(f"forge exit={rc} stages={stages} last_res={last[-1].name if last else None}")
    return rc


# --- nan-scan -------------------------------------------------------------------------------------------------------
def nan_scan(run: Path) -> dict:
    """残差履歴 (段ごとの residual_history_<段>.csv と本段) の rms_* 列 (rms_dq_* を除く) で最初の非有限の行と、最終 res の
    /VALUE の非有限・ρ ≤ 0・T ≤ 0・P ≤ 0。"""
    import h5py
    out = {"run": run.name, "date": now(), "residual": {}, "final_res": None}
    files = sorted(run.glob("residual_history_*.csv"))
    if (run / "residual_history.csv").is_file():
        files.append(run / "residual_history.csv")
    first = None
    for f in files:
        if f.name.endswith("_segment.csv"):
            continue
        with open(f) as fh:
            rd = csv.reader(fh)
            head = [h.strip() for h in next(rd, [])]
            cols = [i for i, h in enumerate(head) if h.startswith("rms_") and not h.startswith("rms_dq_")]
            si = head.index("step") if "step" in head else (head.index("iter") if "iter" in head else 0)
            n, bad = 0, None
            for row in rd:
                n += 1
                for i in cols:
                    try:
                        v = float(row[i])
                    except (ValueError, IndexError):
                        v = float("nan")
                    if not math.isfinite(v):
                        bad = {"row": n, "step": row[si] if si < len(row) else None, "col": head[i]}
                        break
                if bad:
                    break
        out["residual"][f.name] = {"rows": n, "first_nonfinite": bad}
        if bad and first is None:
            first = {"file": f.name, **bad}
    out["first_nonfinite"] = first
    rs = res_files(run)
    if rs:
        fr = {"file": rs[-1].name, "nonfinite": [], "nonpositive": []}
        with h5py.File(rs[-1], "r") as h:
            for k in h["VALUE"]:
                a = np.array(h["VALUE"][k])
                if a.dtype.kind == "f" and not np.all(np.isfinite(a)):
                    fr["nonfinite"].append(k)
            for k in ("ro", "T", "P"):
                if k in h["VALUE"] and not np.all(np.array(h["VALUE"][k]) > 0):
                    fr["nonpositive"].append(k)
        out["final_res"] = fr
    nan_files = [p.name for p in run.glob("res_nan_*.h5")]
    out["res_nan_files"] = nan_files
    ok = first is None and not nan_files and (out["final_res"] is None or not (out["final_res"]["nonfinite"] or out["final_res"]["nonpositive"]))
    out["VERDICT"] = "CLEAN" if ok else "NAN"
    jdump(run / NAN_SCAN, out)
    print(f"[nan-scan] {run.name}: {out['VERDICT']}" + (f" 最初の非有限 {first}" if first else "")
          + (f" 最終 res {out['final_res']}" if out["final_res"] and not ok else ""))
    return out


# --- record-series (判定なしの記録) ----------------------------------------------------------------------------------
REC_X_PW = (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0)          # 壁圧を記録する x/r_t (スロート付近)
REC_X_DE = (0.0, 0.5, 1.0, 2.0, 5.0)                          # 排除厚の補正量 (抽出 δ_E) を記録する x/r_t
REC_ETA_SONIC = (0.0, 0.5, 0.9)                               # 音速線を探す r/r_w


def _trap(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def _record_one(args):
    run, res, euler, Pt = args
    sys.path.insert(0, str(DESIGN))
    from forge_design.report.nozzle_report import load_field, eta_line
    G = load_field(run, res)
    X, R, V, S = G["X"], G["R"], G["V"], G["S"]
    info = G["info"]
    xE = float(info.get("x_E", 40.0))
    row = {"step": step_of(res)}
    # 質量流量 [kg/s]: 列ごとに 2π ∫ ρ U_x r dr (列は x 一定の縦線; r は m)
    md = np.array([2 * np.pi * _trap(V["ro"][i] * V["Ux"][i] * R[i] * S, R[i] * S) for i in range(X.shape[0])])
    xa = X[:, 0]
    sel = (xa > -2.5) & (xa < xE)
    row["mdot_throat_kgs"] = float(md[int(np.argmin(np.abs(xa)))])
    row["mdot_median_kgs"] = float(np.median(md[sel]))
    row["mdot_spread_rel"] = float((md[sel].max() - md[sel].min()) / np.median(md[sel]))
    row["mdot_exit_kgs"] = float(md[-1])
    # 音速線: r/r_w = η の線で M = 1 を最初に横切る x (x ∈ [−3, 3]、線形補間)
    xq = np.linspace(-3.0, 3.0, 6001)
    for eta in REC_ETA_SONIC:
        m = eta_line(G, "M", eta, xq) - 1.0
        k = np.flatnonzero((m[:-1] < 0) & (m[1:] >= 0))
        row[f"x_sonic_eta{eta:g}"] = float(xq[k[0]] - m[k[0]] * (xq[k[0] + 1] - xq[k[0]]) / (m[k[0] + 1] - m[k[0]])) if k.size else float("nan")
    # スロート付近の壁圧 P_w/Pt (壁の節点 j = nj−1 を x で線形補間)
    for x in REC_X_PW:
        row[f"pw_over_pt_x{x:g}"] = float(np.interp(x, X[:, -1], V["P"][:, -1]) / Pt)
    # 排除厚の補正量: CFD から抽出した δ_E (未緩和・平滑化、r_t 単位; exitM_sampling_ab.py と同じ extract_and_merge)
    if euler is not None:
        from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
        with tempfile.TemporaryDirectory(prefix=f"ns_n012_dE_{Path(run).name[:8]}_{row['step']}_") as dd:
            dd = Path(dd)
            for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"):
                shutil.copy(Path(run) / f, dd / f)
            os.symlink((Path(run) / "nozzle.h5").resolve(), dd / "nozzle.h5")
            os.symlink((Path(run) / res).resolve(), dd / res)
            extract_and_merge(dd, euler, band_select="edge")
            nx = read_delta_r_next(dd / "delta_r_next.csv")
        xs, de = nx["x_rt"], nx["delta_E"]
        for x in REC_X_DE:
            row[f"deltaE_rt_x{x:g}"] = float(np.interp(x, xs, de)) if xs.min() <= x <= xs.max() else float("nan")
    else:
        for x in REC_X_DE:
            row[f"deltaE_rt_x{x:g}"] = float("nan")
    return row


def record_series(run: Path, euler: Path | None, out_csv: Path, nproc: int = 3) -> Path:
    """判定に使わない記録の時系列 (U4 の「記録する量」): 質量流量、音速線の位置、スロート付近の壁圧、排除厚の補正量 (抽出 δ_E)。
    1 行 = 1 スナップショット (res_0 を除く本段の res)。加えて、壁に足した δ_r (delta_r_initial.csv) を同じ x で JSON に記録する。"""
    from concurrent.futures import ProcessPoolExecutor
    rec = jload(run / RECORD) if (run / RECORD).is_file() else {}
    ys = yaml_strict()
    prob = rec.get("problem")
    if prob:
        Pt = float(ys.load((HERE / prob).read_text())["spec"]["Pt"])
    else:
        Pt = float(ys.load((HERE / MK.BASE["dry"]).read_text())["spec"]["Pt"])
    files = [f.name for f in res_files(run) if f.name != "res_0.h5"]
    if not files:
        raise SystemExit(f"{run}: 本段の res が無い")
    args = [(str(run), f, (str(euler) if euler else None), Pt) for f in files]
    if nproc > 1:
        with ProcessPoolExecutor(nproc) as ex:
            rows = list(ex.map(_record_one, args))
    else:
        rows = [_record_one(a) for a in args]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    side = {"run": run.name, "euler_ref": (str(euler) if euler else None), "Pt": Pt, "date": now(),
            "columns": list(rows[0]), "note": "判定に使わない記録 (plan U4「記録する量」)。δ_E は extract_and_merge (band_select edge)"}
    dr = run / "delta_r_initial.csv"
    if dr.is_file():
        a = np.loadtxt(dr, delimiter=",", skiprows=1)
        side["delta_r_applied_rt"] = {f"x{x:g}": float(np.interp(x, a[:, 0], a[:, 1])) for x in REC_X_DE}
    pi = jload(run / "prepare_info.json")
    side["throat_physical"] = pi.get("throat_physical")
    jdump(out_csv.with_suffix(".json"), side)
    print(f"[record-series] {run.name}: {len(rows)} 行 → {out_csv}")
    return out_csv


# --- main -----------------------------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("prep-dry"); a1.add_argument("cond"); a1.add_argument("run"); a1.add_argument("--md-offset", required=True)
    a1.add_argument("--ic-src-run"); a1.add_argument("--dry", action="store_true")
    a2 = sub.add_parser("ic-check"); a2.add_argument("run"); a2.add_argument("--dry", action="store_true")
    a3 = sub.add_parser("verify-set"); a3.add_argument("runs", nargs=3); a3.add_argument("--md-offset", required=True)
    a3.add_argument("--dry", action="store_true"); a3.add_argument("--out")
    a4 = sub.add_parser("prep-cond"); a4.add_argument("cond"); a4.add_argument("src"); a4.add_argument("run")
    a4.add_argument("--md-offset", required=True); a4.add_argument("--dry", action="store_true")
    a5 = sub.add_parser("prep-ext"); a5.add_argument("src"); a5.add_argument("run"); a5.add_argument("--dry", action="store_true")
    a6 = sub.add_parser("run"); a6.add_argument("run")
    a7 = sub.add_parser("nan-scan"); a7.add_argument("run")
    a8 = sub.add_parser("record-series"); a8.add_argument("run"); a8.add_argument("--euler"); a8.add_argument("--out")
    a8.add_argument("--nproc", type=int, default=3)
    a = ap.parse_args(argv)
    P = lambda s: (HERE / s) if not os.path.isabs(s) else Path(s)  # noqa: E731
    if a.cmd == "prep-dry":
        prep_dry(a.cond, P(a.run), a.md_offset, a.ic_src_run, a.dry)
    elif a.cmd == "ic-check":
        return 0 if ic_check(P(a.run), a.dry)["VERDICT"] == "OK" else 2
    elif a.cmd == "verify-set":
        return 0 if verify_set([P(r) for r in a.runs], a.md_offset, a.dry, P(a.out) if a.out else None)["VERDICT"] == "OK" else 2
    elif a.cmd == "prep-cond":
        prep_cond(a.cond, P(a.src), P(a.run), a.md_offset, a.dry)
    elif a.cmd == "prep-ext":
        prep_ext(P(a.src), P(a.run), a.dry)
    elif a.cmd == "run":
        return run_one(P(a.run))
    elif a.cmd == "nan-scan":
        return 0 if nan_scan(P(a.run))["VERDICT"] == "CLEAN" else 1
    elif a.cmd == "record-series":
        run = P(a.run)
        out = P(a.out) if a.out else run / "ns_n012_record_series.csv"
        record_series(run, P(a.euler) if a.euler else None, out, a.nproc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
