"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36r の前半 (解析と表の案、CFD 0 step): V-h370-diag (run_0951_h370_ns_b0) の保持の窓
112,500〜152,500 の 17 出力から壁の局所の τ_t・y₁⁺ の時間の包絡を作り、上流側 (x < 0) の第一層の表 `wall_first_frac_table` の案を決めて、
格子を生成器 (`mesh2d`) で作り直して品質と y₁⁺ の見積もりを出す (codex diagnose `notes/reviews/2026-10-15-h370-wallres-timeband-result-diagnose.md`
の Major 4・5: 直すのは上流の表の値、`wall_first_up_*` は併用できない)。新しい NS は回さない。結果の解釈はしない (事実だけを JSON に書く)。

envelope (B の case dir で): 各出力 s について `check_wall_resolution.py <run> --step s --groups wall --weight area --profile-csv` を回し
  (#36o の `h370_wallres_split.py` と同じ呼び方・同じ定義: y₁ = 第一内部節点までの壁法線の距離、τ_t = 接線の壁応力、μ_w = 体積出力の vis_lam)、
  壁の節点ごとに τ_t・ρ_w・μ_w・y₁⁺ = y₁ √(ρ_w |τ_t|)/μ_w の時間の最大・最小・中央を出す。節点の座標・y₁・面積の重みが全出力で同じことを確かめる。
  → `_band_ab/half_inlet/h370_wall_envelope.json`・`.png`。
design (B の case dir で、envelope の後): run の `wall_repr.json`・`prepare_info.json` の格子のブロック (問題 YAML の mesh 節と照合) から、
  元の表で格子を作り直して生成器の sha256 (prepare_info.mesh.hashes.generator_sha256) と壁の節点の座標が run と一致することを確かめ (対照)、
  案ごとに上流の表の値を決める:
  - 各壁の節点 i (x < 0) で、流れを包絡の最大に固定した一次の見積もり y₁⁺_new(i) = y₁⁺_max(i) · y₁_new(i) / y₁_old(i) ≤ 目標 になる
    第一層の比 fr_req(i) = fr_old(i) · 目標 / y₁⁺_max(i) を求める (y₁ = fr · r_w は生成器の定義; #36o で y₁/設計 d₁ = 1 ± 1e−9)。
  - 表の x < 0 の点の値 v_k (log) を変数に、生成器と同じ log 線形の補間が全節点で fr_req 以下という線形の制約の下で、Σ (台形の重み) log v_k
    を最大にする線形計画で決める (HiGHS)。既定は元の値より厚くしない (u_k ≤ log v_old)。x ≥ 0 の点は変えない。値は 6 桁に切り下げる (制約を保つ向き)。
  - 案: 目標 0.8・0.9 × 表の点 (今の点 / x < 0 を 0.125 r_t 刻みに足す) (T08・T09・T08f・T09f)。表だけでは上流の壁のセルの AR が 5000 を
    超えると分かったので、参考 (方針の外、採否は上位) として上流の流れ方向の細分 (x_density_table の上流の値と ni) を併せた案 (X08・X09・
    X08f・X09f、`x_refine_density`) も作る。細分しない区間の列は元の格子の列に丸めの範囲で一致する (列の番号が足した数だけずれる)。
  作り直した格子で: 生成の座標で AR・スキュー (check_mesh_quality の関数そのもの)、x ≥ 0 の列が元とビット同一か、径方向の等比の比の最大、
  第一層の厚さ y₁ (生成の座標) と y₁⁺ の見積もり (最大・1 超えと目標超えの面積)、節点数。--convert を付けると B のビルド済みの変換器で
  node 変換して `check_mesh_quality.py --ar-max 5000` の VERDICT、AR ≤ 5000 の例外の条件 (`runner.mesh_ar_relax_check`: AR > 1000 の全セルで
  スキュー ≤ 0.1・入口と出口が平面)、`check_wall_resolution.py --geometry-only --bcond` の VERDICT と節点ごとの y₁ (同じツールの
  `wall_first_distance`) を出す。対照 (元の表) の変換後の座標の sha256 が run の記録と同じことも確かめる。
  → `_band_ab/half_inlet/h370_wall_table_design.json`・`.png`。

使い方 (B の case dir = ~/forge-wallfit/case/45.isobutane_m6_d155 で、system の python3; 台本は ~/forge-g0cchain の版):
  python3 ~/forge-g0cchain/case/45.isobutane_m6_d155/h370_wall_envelope.py envelope [--run run_0951_h370_ns_b0]
  python3 ~/forge-g0cchain/case/45.isobutane_m6_d155/h370_wall_envelope.py design [--convert] [--workdir ~/forge-g0cchain/_tmp_36r]
  出力の既定は台本の置き場所 (~/forge-g0cchain/case/45.isobutane_m6_d155/_band_ab/half_inlet/)。
終了コード 0 書けた / 2 入力が足りない・照合の不一致 (JSON の errors)。
"""
import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36r の前半 (codex diagnose "
        "notes/reviews/2026-10-15-h370-wallres-timeband-result-diagnose.md の Major 4・5)")
RUN = "run_0951_h370_ns_b0"
STEPS = list(range(112500, 152501, 2500))       # 保持の窓の 17 出力 (§6 #36p と同じ)
PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
PROBLEM_SHA_RUN = "0d784f7c49922444f9e476579b3fee092f4188e7717354caf37184ba84250eea"   # run の H370_PRECHECK.json の N8_sha.problem
OUT_ENV = "_band_ab/half_inlet/h370_wall_envelope.json"
OUT_DES = "_band_ab/half_inlet/h370_wall_table_design.json"
CONVERTER = Path.home() / "forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge"
CONVERTER_SHA = "711f44389209e16d6ae3f8111f21faec31c9d3e0ba3382722e8d33a69579d7bf"    # run の H370_PRECHECK.json の converter_sha256
AR_MAX = 5000.0
BINS_RT = [-12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.0, -2.0,
           -1.0, 0.0, 1.0, 3.0, 5.0, 10.0, 40.0, 94.0, None]   # #36o と同じ区切り [r_t] (最後は出口まで)
# 案 (結果を見る前に決めた並び): 目標 0.8・0.9 × 表の点 (今の点 / x < 0 を 0.125 r_t 刻み)。どれも元より厚くしない・x ≥ 0 は変えない
VARIANTS = [
    {"name": "T08", "target": 0.8, "spacing": None, "thicken": False},
    {"name": "T09", "target": 0.9, "spacing": None, "thicken": False},
    {"name": "T08f", "target": 0.8, "spacing": 0.125, "thicken": False},
    {"name": "T09f", "target": 0.9, "spacing": 0.125, "thicken": False},
    # 参考 (方針の外、表だけでは AR ≤ 5000 を満たせないと分かった後に足した): 上流の流れ方向の細分 (x_density_table・ni) を併せる。採否は上位
    # ar_goal 4500 は元の格子の設計の見込み (問題 YAML の mesh 節のコメント「AR ≤ 4500 の見込み」) と同じ値
    {"name": "X08", "target": 0.8, "spacing": None, "thicken": False, "x_refine": {"base": "T08", "ar_goal": 4500.0}},
    {"name": "X09", "target": 0.9, "spacing": None, "thicken": False, "x_refine": {"base": "T09", "ar_goal": 4500.0}},
    {"name": "X08f", "target": 0.8, "spacing": 0.125, "thicken": False, "x_refine": {"base": "T08f", "ar_goal": 4500.0}},
    {"name": "X09f", "target": 0.9, "spacing": 0.125, "thicken": False, "x_refine": {"base": "T09f", "ar_goal": 4500.0}},
]
VERDICT_RE = re.compile(r"^VERDICT")


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fnum(v, sig: int = 10):
    """JSON に書ける数 (非有限は None)。"""
    if v is None:
        return None
    v = float(v)
    return float(f"{v:.{sig}g}") if math.isfinite(v) else None


def flist(a, sig: int = 7) -> list:
    return [fnum(v, sig) for v in np.asarray(a, dtype=float)]


def jwrite(p: Path, d) -> None:
    """JSON を書く。節点ごとの配列 (キー nodes) だけは詰めて書く (大きさを抑える)。"""
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    dflt = lambda o: o.tolist() if hasattr(o, "tolist") else str(o)  # noqa: E731
    d2 = dict(d)
    nodes = d2.pop("nodes", None)
    txt = json.dumps(d2, ensure_ascii=False, indent=1, default=dflt)
    if nodes is not None:
        txt = txt[:-2] + ',\n "nodes": ' + json.dumps(nodes, ensure_ascii=False, separators=(",", ":"), default=dflt) + "\n}"
    tmp.write_text(txt)
    tmp.replace(p)


def git_head() -> dict:
    def g(*a):
        r = subprocess.run(["git", *a], cwd=str(ROOT), capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    st = g("status", "--short", "--", "case/45.isobutane_m6_d155/h370_wall_envelope.py")
    return {"head": g("rev-parse", "HEAD"), "script_dirty": bool(st), "script_sha256": sha256_file(__file__)}


def bins_of(x_rt: np.ndarray):
    for lo, hi in zip(BINS_RT[:-1], BINS_RT[1:]):
        m = (x_rt >= lo) & ((x_rt < hi) if hi is not None else True)
        if m.any():
            yield lo, hi, m


def runs_of(mask) -> list:
    idx = np.flatnonzero(mask)
    if not len(idx):
        return []
    out, s, p = [], int(idx[0]), int(idx[0])
    for i in idx[1:]:
        if i != p + 1:
            out.append((s, p))
            s = int(i)
        p = int(i)
    out.append((s, p))
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
# envelope
# ---------------------------------------------------------------------------------------------------------------------------------
def run_tool(tools: Path, case: Path, run: str, step: int, prof: Path) -> dict:
    cmd = [sys.executable, str(tools / "check_wall_resolution.py"), run, "--step", str(step), "--groups", "wall", "--weight", "area",
           "--profile-csv", str(prof)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(case))
    lines = [ln.rstrip() for ln in (r.stdout + r.stderr).splitlines() if ln.strip()]
    over = None
    for ln in lines:
        m = re.search(r"y1\+ > [0-9.]+ が面積割合 ([0-9.]+) %", ln)
        if m:
            over = float(m.group(1))
    return {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdict": [ln.strip() for ln in lines if VERDICT_RE.match(ln.strip())],
            "area_over_pct_printed": over, "lines": [ln for ln in lines if "分布 ->" not in ln][:12]}


def read_profile(p: Path) -> dict:
    d = np.genfromtxt(p, delimiter=",", names=True)
    o = np.argsort(d["x"], kind="stable")
    return {k: np.asarray(d[k], dtype=float)[o] for k in d.dtype.names}


def envelope(a) -> int:
    case = Path(a.case_dir).resolve()
    tools = Path(a.tools)
    run = case / a.run
    out_p = Path(a.out) if a.out else HERE / OUT_ENV
    steps = [int(s) for s in a.steps.split(",")] if a.steps else STEPS
    R = {"plan": PLAN, "created": now(), "kind": "包絡 (§5.1 #36r の前半)。CFD 0 step・結果の解釈はしない", "run": a.run,
         "run_dir": str(run), "steps": steps, "script": git_head(), "errors": []}
    err = R["errors"]
    need = [run / "nozzle.h5", run / "solverConfig.yaml", run / "prepare_info.json"]
    need += [run / f"res_wall_3_{s}.h5" for s in steps] + [run / f"res_{s}.h5" for s in steps]
    for p in need:
        if not p.is_file():
            err.append(f"{p} が無い")
    if err:
        jwrite(out_p, R)
        return 2
    info = json.loads((run / "prepare_info.json").read_text())
    rt = float(info["scale_m"])
    R["inputs"] = {"r_t_m": rt, "tools_dir": str(tools), "tool_sha256": sha256_file(tools / "check_wall_resolution.py"),
                   "sha256": {n: sha256_file(run / n) for n in ("nozzle.h5", "solverConfig.yaml", "bcondConfig.yaml", "prepare_info.json")},
                   "wall_dump_sha256": {str(s): sha256_file(run / f"res_wall_3_{s}.h5") for s in steps}}
    P, per = [], []
    with tempfile.TemporaryDirectory(prefix="h370_env_") as td:
        for s in steps:
            t = run_tool(tools, case, a.run, s, Path(td) / f"p{s}.csv")
            pp = Path(td) / f"p{s}_wall.csv"
            if t["rc"] not in (0, 1) or not pp.is_file():
                err.append(f"step {s}: check_wall_resolution.py が判定できない (rc {t['rc']}) {t['verdict']}")
                per.append({"step": s, "tool": t})
                continue
            q = read_profile(pp)
            P.append(q)
            g = np.isfinite(q["y1plus"])
            w = q["area_weight"]
            im = int(np.nanargmax(np.where(g, q["y1plus"], -np.inf)))
            per.append({"step": s, "tool": t, "n_eval": int(g.sum()),
                        "area_over1_pct": fnum(100 * w[g & (q["y1plus"] > 1.0)].sum() / w[g].sum()),
                        "area_over08_pct": fnum(100 * w[g & (q["y1plus"] > 0.8)].sum() / w[g].sum()),
                        "y1plus_max": fnum(q["y1plus"][im]), "x_rt_max": fnum(q["x"][im] / rt),
                        "tau_t_max": fnum(np.nanmax(q["tau_t"])), "x_rt_tau_max": fnum(q["x"][int(np.nanargmax(q["tau_t"]))] / rt)})
            print(f"step {s}: rc {t['rc']} {t['verdict']} 1 超えの面積 {per[-1]['area_over1_pct']} % 最大 {per[-1]['y1plus_max']} "
                  f"(x {per[-1]['x_rt_max']} r_t)", flush=True)
    R["per_step"] = per
    if err:
        jwrite(out_p, R)
        return 2
    # 節点の座標・y₁・面積の重みが全出力で同じか (格子は同じなので一致するはず)
    base = P[0]
    same = {k: fnum(max(float(np.max(np.abs(q[k] - base[k]))) for q in P)) for k in ("x", "y", "y1_m", "area_weight", "mu_w", "T_w")}
    R["static_check_max_abs_diff"] = same
    for k in ("x", "y", "y1_m", "area_weight"):
        if same[k] != 0.0:
            err.append(f"{k} が出力の間で一致しない (最大差 {same[k]})")
    S = {k: np.stack([q[k] for q in P]) for k in ("tau_t", "rho_w", "mu_w", "y1plus")}
    x, r, y1, w = base["x"], base["y"], base["y1_m"], base["area_weight"]
    xr = x / rt
    env = {}
    for k in ("tau_t", "rho_w", "mu_w", "y1plus"):
        env[k] = {"max": np.max(S[k], 0), "min": np.min(S[k], 0), "median": np.median(S[k], 0)}
    imax_t = np.argmax(S["y1plus"], 0)
    # y₁⁺ の時間の最大は「流れを各時刻に固定した y₁⁺」の最大。τ_t の最大と ρ_w の組を固定した値 (τ_max と同じ時刻の ρ) との差も記録
    it_tau = np.argmax(S["tau_t"], 0)
    jj = np.arange(len(x))
    y1p_at_taumax = S["y1plus"][it_tau, jj]
    rel_w = lambda k: (env[k]["max"] - env[k]["min"]) / env[k]["median"]  # noqa: E731
    R["wall"] = {"n_nodes": int(len(x)), "area_total_m2": fnum(w.sum()),
                 "mu_w": {"min": fnum(S["mu_w"].min()), "max": fnum(S["mu_w"].max())},
                 "rho_w_rel_time_width_max": fnum(np.max(rel_w("rho_w"))),
                 "y1plus_env_max_vs_at_tau_max_max_rel_diff": fnum(np.max(np.abs(env["y1plus"]["max"] / y1p_at_taumax - 1.0))),
                 "note": ("包絡は 17 出力の節点ごとの最大・最小・中央。y₁ は時間によらない (static_check)。y₁⁺ の時間の最大と、τ_t が最大の時刻の y₁⁺ の差は "
                          "ρ_w の時間の変化だけ (相対の最大を記録)。")}
    ymax = env["y1plus"]["max"]
    over = ymax > 1.0
    R["envelope_over1"] = {"area_pct": fnum(100 * w[over].sum() / w.sum()), "n_nodes": int(over.sum()),
                           "regions": [{"x_rt": [fnum(xr[a_]), fnum(xr[b_])], "n": int(b_ - a_ + 1),
                                        "area_pct": fnum(100 * w[a_:b_ + 1].sum() / w.sum()),
                                        "y1plus_max": fnum(ymax[a_:b_ + 1].max()),
                                        "x_rt_max": fnum(xr[a_ + int(np.argmax(ymax[a_:b_ + 1]))])} for a_, b_ in runs_of(over)]}
    R["envelope_over08"] = {"area_pct": fnum(100 * w[ymax > 0.8].sum() / w.sum()), "n_nodes": int((ymax > 0.8).sum()),
                            "regions": [{"x_rt": [fnum(xr[a_]), fnum(xr[b_])], "n": int(b_ - a_ + 1)} for a_, b_ in runs_of(ymax > 0.8)]}
    bins = []
    for lo, hi, m in bins_of(xr):
        i1 = np.flatnonzero(m)[int(np.argmax(ymax[m]))]
        i2 = np.flatnonzero(m)[int(np.argmax(env["tau_t"]["max"][m]))]
        bins.append({"x_rt": [lo, hi if hi is not None else fnum(xr[m].max())], "n": int(m.sum()), "area_m2": fnum(w[m].sum()),
                     "y1_m": [fnum(y1[m].min()), fnum(y1[m].max())],
                     "y1plus_env_max": fnum(ymax[i1]), "x_rt_at_y1plus_max": fnum(xr[i1]),
                     "y1plus_at_that_node_min_median_max": [fnum(env["y1plus"]["min"][i1]), fnum(env["y1plus"]["median"][i1]), fnum(ymax[i1])],
                     "y1plus_rel_time_width_max": fnum(np.max(rel_w("y1plus")[m])),
                     "tau_t_env_max": fnum(env["tau_t"]["max"][i2]), "x_rt_at_tau_max": fnum(xr[i2]),
                     "tau_t_at_that_node_min_median_max": [fnum(env["tau_t"]["min"][i2]), fnum(env["tau_t"]["median"][i2]),
                                                           fnum(env["tau_t"]["max"][i2])],
                     "tau_t_rel_time_width_max": fnum(np.max(rel_w("tau_t")[m])),
                     "area_over1_env_pct_of_wall": fnum(100 * w[m & over].sum() / w.sum()),
                     "area_over08_env_pct_of_wall": fnum(100 * w[m & (ymax > 0.8)].sum() / w.sum())})
    R["bins"] = bins
    i = int(np.argmax(ymax))
    R["max_node"] = {"i": i, "x_rt": fnum(xr[i]), "r_rt": fnum(r[i] / rt), "y1_m": fnum(y1[i]),
                     "y1plus_min_median_max": [fnum(env["y1plus"]["min"][i]), fnum(env["y1plus"]["median"][i]), fnum(ymax[i])],
                     "tau_t_min_median_max": [fnum(env["tau_t"]["min"][i]), fnum(env["tau_t"]["median"][i]), fnum(env["tau_t"]["max"][i])],
                     "rho_w_min_max": [fnum(env["rho_w"]["min"][i]), fnum(env["rho_w"]["max"][i])], "mu_w": fnum(env["mu_w"]["max"][i]),
                     "step_of_max": steps[int(imax_t[i])],
                     "y1_for_target_0p8_m": fnum(0.8 * y1[i] / ymax[i]), "y1_for_target_0p9_m": fnum(0.9 * y1[i] / ymax[i])}
    # 節点ごとの包絡 (x 昇順)。全出力の値そのもの (17 × 節点) は npz (B だけ) に残し、sha256 を書く
    npz = out_p.with_suffix(".npz")
    np.savez_compressed(npz, steps=np.asarray(steps), x=x, r=r, y1=y1, area_weight=w, **{k: S[k] for k in S})
    R["npz"] = {"file": npz.name, "sha256": sha256_file(npz), "note": "17 出力の節点ごとの tau_t・rho_w・mu_w・y1plus (B だけ、commit しない)"}
    R["nodes"] = {"note": ("x 昇順の壁の節点 (check_wall_resolution.py の profile CSV の並び; 全出力で同じ)。座標は 13 桁・y₁・面積・y₁⁺ の最大は 10 桁・"
                           "ほかは 7 桁に丸め。μ_w は全節点・全出力で同じ (wall.mu_w)"),
                  "x_m": flist(x, 13), "r_m": flist(r, 13), "y1_m": flist(y1, 10), "area_weight_m2": flist(w, 10),
                  "tau_t_max": flist(env["tau_t"]["max"]), "tau_t_min": flist(env["tau_t"]["min"]), "tau_t_median": flist(env["tau_t"]["median"]),
                  "y1plus_max": flist(ymax, 10), "y1plus_min": flist(env["y1plus"]["min"]), "y1plus_median": flist(env["y1plus"]["median"])}
    R["facts_summary"] = [
        f"17 出力の包絡: 1 超えの面積 {R['envelope_over1']['area_pct']} % ({R['envelope_over1']['n_nodes']} 節点)、0.8 超え "
        f"{R['envelope_over08']['area_pct']} %",
        f"最大 y₁⁺ {R['max_node']['y1plus_min_median_max']} (最小・中央・最大) at x {R['max_node']['x_rt']} r_t、y₁ {R['max_node']['y1_m']} m、"
        f"τ_t {R['max_node']['tau_t_min_median_max']} Pa、目標 0.8 の y₁ {R['max_node']['y1_for_target_0p8_m']} m",
        "出力ごとの 1 超えの面積 [%]: " + ", ".join(f"{p['step']}: {p['area_over1_pct']}" for p in per),
    ]
    png = out_p.with_suffix(".png")
    try:
        plot_envelope(png, xr, y1, env, rt, per, steps)
        R["png"] = png.name
    except Exception as e:                      # noqa: BLE001
        R["png_error"] = repr(e)
    jwrite(out_p, R)
    print("\n".join(R["facts_summary"]))
    print(f"→ {out_p}")
    return 0 if not err else 2


def plot_envelope(png: Path, xr, y1, env, rt, per, steps) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 1, figsize=(11, 13))
    for a_, sl, ttl in ((ax[0], slice(None), "whole wall"), (ax[1], xr < 0.0, "x < 0 (pipe + contraction)")):
        a_.fill_between(xr[sl], env["y1plus"]["min"][sl], env["y1plus"]["max"][sl], color="C0", alpha=0.35, lw=0, label="y1+ min-max (17 outputs)")
        a_.plot(xr[sl], env["y1plus"]["median"][sl], "C0-", lw=0.8, label="y1+ median")
        a_.axhline(1.0, color="r", lw=0.8, ls="--", label="1")
        a_.axhline(0.8, color="k", lw=0.8, ls=":", label="0.8")
        a_.set_ylabel("local y1+")
        a_.set_title(f"run_0951 wall y1+ envelope, steps {steps[0]}-{steps[-1]} ({ttl})")
        a_.legend(fontsize=8, loc="upper right")
        a_.grid(alpha=0.3)
    sl = xr < 0.0
    ax[2].fill_between(xr[sl], env["tau_t"]["min"][sl], env["tau_t"]["max"][sl], color="C1", alpha=0.35, lw=0, label="tau_t min-max")
    ax[2].set_yscale("log")
    ax[2].plot(xr[sl], env["tau_t"]["median"][sl], "C1-", lw=0.8, label="tau_t median")
    ax[2].set_ylabel("tangential wall stress [Pa]")
    ax2 = ax[2].twinx()
    ax2.plot(xr[sl], y1[sl] * 1e6, "k-", lw=0.8, label="y1 [um]")
    ax2.set_ylabel("y1 [um]")
    ax[2].legend(fontsize=8, loc="upper left")
    ax2.legend(fontsize=8, loc="upper right")
    ax[2].grid(alpha=0.3)
    ax[3].plot([p["step"] for p in per], [p["area_over1_pct"] for p in per], "o-", label="area y1+ > 1 [%]")
    ax[3].plot([p["step"] for p in per], [p["area_over08_pct"] for p in per], "s-", label="area y1+ > 0.8 [%]")
    ax[3].set_xlabel("step")
    ax[3].set_ylabel("area fraction [%]")
    ax[3].legend(fontsize=8)
    ax[3].grid(alpha=0.3)
    for a_ in ax[:3]:
        a_.set_xlabel("x / r_t")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------------------------
# design
# ---------------------------------------------------------------------------------------------------------------------------------
def floor_sig(v: float, sig: int = 6) -> float:
    """有効数字 sig 桁に切り下げる (表の値を薄い側へ丸めて制約を保つ)。"""
    e = math.floor(math.log10(v)) - (sig - 1)
    q = math.floor(v / 10.0 ** e)
    out = float(f"{q}e{e}")
    while out > v:                               # 10 進の表記の丸めで上にずれたら 1 単位下げる
        q -= 1
        out = float(f"{q}e{e}")
    return out


def design_table(old_tbl: np.ndarray, xs_w: np.ndarray, fr_old: np.ndarray, y1p_max: np.ndarray, target: float, spacing, thicken: bool) -> dict:
    """上流 (x < 0) の表の値を線形計画で決める。戻り {table (list), lp, ...}。xs_w・fr_old・y1p_max は壁の各 station (r_t 単位の x)。"""
    from scipy.optimize import linprog
    tx, tv = old_tbl[:, 0], old_tbl[:, 1]
    if spacing is None:
        nx = tx.copy()
    else:
        add = np.arange(-12.5 + spacing, 0.0 - 1e-9, spacing)
        nx = np.unique(np.concatenate([tx, add[(add > tx[0]) & (add < 0.0)]]))
    old_at = np.exp(np.interp(nx, tx, np.log(tv)))       # 元の表 (log 線形) の値を新しい点で
    old_at[np.isin(nx, tx)] = tv[np.isin(tx, nx)]         # 既存の点は元の値そのもの
    var = np.flatnonzero(nx < 0.0)                        # 変数 (x < 0)。x ≥ 0 の点は元のまま
    fixed = np.log(old_at)
    nv = len(var)
    col = {int(k): j for j, k in enumerate(var)}
    sel = (xs_w < 0.0) & np.isfinite(y1p_max) & (y1p_max > 0)
    lo_rows, rhs = [], []
    req = np.log(fr_old * target / np.where(y1p_max > 0, y1p_max, np.nan))
    for i in np.flatnonzero(sel):
        k = int(np.clip(np.searchsorted(nx, xs_w[i], side="right") - 1, 0, len(nx) - 2))
        t = (xs_w[i] - nx[k]) / (nx[k + 1] - nx[k])
        row = np.zeros(nv)
        b = req[i]
        for kk, c in ((k, 1.0 - t), (k + 1, t)):
            if kk in col:
                row[col[kk]] += c
            else:
                b -= c * fixed[kk]
        lo_rows.append(row)
        rhs.append(b)
    A = np.array(lo_rows)
    bnd = [(None, (np.log(old_at[k]) if not thicken else None)) for k in var]
    wq = np.zeros(len(nx))
    dx = np.diff(nx)
    wq[:-1] += 0.5 * dx
    wq[1:] += 0.5 * dx
    c = -wq[var]
    res = linprog(c, A_ub=A, b_ub=np.array(rhs), bounds=bnd, method="highs")
    if res.status != 0:
        return {"ok": False, "lp": {"status": int(res.status), "message": res.message}}
    newv = old_at.copy()
    for j, k in enumerate(var):
        v = float(np.exp(res.x[j]))
        newv[k] = old_at[k] if (not thicken and v >= old_at[k] * (1 - 1e-12)) else floor_sig(v)
    tbl = [[float(a_), float(b_)] for a_, b_ in zip(nx, newv)]
    return {"ok": True, "table": tbl, "n_points": len(tbl), "n_var": nv, "n_constraints": int(len(rhs)),
            "n_changed": int(np.sum(newv[var] != old_at[var])),
            "lp": {"status": int(res.status), "message": res.message, "objective": fnum(res.fun)},
            "binding": int(np.sum(A @ res.x >= np.array(rhs) - 1e-9))}


def x_refine_density(block: dict, wall, base: dict, ar_goal: float, ar_ctrl: np.ndarray) -> dict:
    """参考の案 (表の案の基 base の格子で、表を薄くしたために AR が ar_goal と元の格子の AR (ar_ctrl、同じ列) の大きい方を超えた壁のセルが
    ある区間だけ流れ方向を細かくする; 方針の外、採否は上位)。表を変えていない区間 (元の格子で AR 4574 のスロート付近など) は細かくしない。
    x_density_table の各点 x_k の値を d_k (1 + α (F_k − 1)) にする。F_k = max(1, AR / max(ar_goal, AR_元)) の、点の両隣の区間 [x_{k−1}, x_{k+1}] に
    掛かる壁のセルでの最大 (線形補間した密度が全区間で要る倍率以上になる)。α ≥ 1 は、足した密度の積分 (mesh2d._x_stations と同じ 8001 点の
    台形) が元の列の間隔 h = ∫d / (ni − 1) のちょうど m 倍になるように決め、ni を m 増やす (細分しない区間の station は元の列そのものに
    丸めの範囲で一致する: 累積の差が定数 m·h なので列の番号が m ずれるだけ)。表の値は丸めない (積分の整数倍を保つため)。"""
    from forge_design.meshing import mesh2d
    M = base["M"]
    dt = np.asarray(block["x_density_table"], dtype=float)
    x_in, x_e = float(wall.x_in), float(wall.x_e)
    ar_w = M["ar2"][:, -1]                                   # 壁に接するセル (列 i〜i+1) の AR
    xa, xb = M["xw"][:-1], M["xw"][1:]
    f = ar_w / np.maximum(ar_goal, ar_ctrl)
    f = np.where(f > 1.0 + 1e-6, f, 1.0)                    # 表を変えていない列の丸めの差 (等比の解き直しで相対 1e−9 程度) では細かくしない
    F = np.ones(len(dt))
    for k in range(len(dt)):
        lo = dt[max(k - 1, 0), 0]
        hi = dt[min(k + 1, len(dt) - 1), 0]
        sel = (xb >= lo) & (xa <= hi)
        if sel.any():
            F[k] = float(f[sel].max())
    grid = np.linspace(x_in, x_e, 8001)
    prm = dict(refine=float(block.get("throat_refine", 3.0)), width=float(block.get("throat_width", 1.5)))

    def integ(vals):
        t = np.c_[dt[:, 0], vals]
        dens = mesh2d._x_density(grid, prm["refine"], prm["width"], table=mesh2d._check_table("x_density_table", t, x_in, x_e))
        return float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(grid)))
    ni = int(block["ni"])
    tot0 = integ(dt[:, 1])
    h = tot0 / (ni - 1)
    inc1 = integ(dt[:, 1] * F) - tot0
    if inc1 <= 0:
        return {"table": dt.tolist(), "ni": ni, "m": 0, "alpha": 1.0, "F_max": 1.0}
    m = int(math.ceil(inc1 / h))
    alpha = m * h / inc1
    vals = dt[:, 1] * (1.0 + alpha * (F - 1.0))
    tbl = [[float(x), float(v)] for x, v in zip(dt[:, 0], vals)]
    ch = F > 1.0
    return {"table": tbl, "ni": ni + m, "m": m, "alpha": fnum(alpha), "ar_goal": ar_goal, "F_max": fnum(F.max()),
            "x_rt_refined": [fnum(dt[ch, 0].min()), fnum(dt[ch, 0].max())] if ch.any() else None,
            "added_integral_over_h": fnum((integ(vals) - tot0) / h),
            "note": "参考 (方針の外): x_density_table の上流の値と ni も変える。採否は上位"}


def mesh_from(block: dict, scale: float, wall, table=None):
    from forge_design.evaluate.runner_axismach import _mesh_params_from
    from forge_design.meshing import mesh2d
    b = dict(block)
    if table is not None:
        b["wall_first_frac_table"] = table
    prm = _mesh_params_from(b, scale, int(b.get("ni", 0)), int(b.get("nj", 0)), float(b.get("wall_first_frac", 0.0)))
    coords, quads, bedges = mesh2d.generate_axisym_mesh(wall, prm)
    return prm, coords, quads, bedges


def generator_sha(coords, quads, bedges) -> str:
    """runner_axismach._mesh_hashes の generator_sha256 と同じ式。"""
    h = hashlib.sha256()

    def upd(name, a):
        a = np.ascontiguousarray(a)
        h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    upd("coords", np.asarray(coords, dtype=np.float64))
    upd("quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        upd(g, np.asarray(bedges[g], dtype=np.int64))
    return h.hexdigest()


def mesh_metrics(prm, coords, quads, scale: float) -> dict:
    """生成の座標の primal セルで AR・スキュー (check_mesh_quality の関数そのもの)、壁の列ごとの量。"""
    import check_mesh_quality as cmq
    from forge_design.meshing import mesh2d
    ni, nj = prm.ni, len(coords) // prm.ni
    conn = np.asarray(quads, dtype=np.int64).ravel()
    offs = np.arange(4, 4 * len(quads) + 1, 4)
    ars, sk = cmq.metrics_vectorized(coords[:, :2], conn, offs, np.full(len(quads), 9), False)
    X = coords[:, 0].reshape(ni, nj)
    Y = coords[:, 1].reshape(ni, nj)
    xw = X[:, -1] / scale
    y1 = np.hypot(X[:, -1] - X[:, -2], Y[:, -1] - Y[:, -2])       # 第一内部節点までの距離 (壁法線の層では法線の距離)
    ar2 = ars.reshape(ni - 1, nj - 1)
    sk2 = sk.reshape(ni - 1, nj - 1)
    fr = mesh2d._first_frac_profile(xw, prm)
    ratio = np.array([mesh2d._radial_fracs_capfixed(nj, float(f), float(prm.axis_cap_frac))[1] for f in fr]) \
        if prm.axis_cap_frac is not None else None
    hi = ars > 1000.0
    return {"ars": ars, "sk": sk, "ar2": ar2, "sk2": sk2, "xw": xw, "y1": y1, "fr": fr, "ratio": ratio, "X": X, "Y": Y,
            "summary": {"cells": int(len(quads)), "nodes": int(len(coords)), "ni": int(ni), "nj": int(nj),
                        "ar_max": fnum(ars.max()), "skew_max": fnum(sk.max()), "n_ar_gt_5000": int((ars > AR_MAX).sum()),
                        "n_ar_gt_1000": int(hi.sum()), "skew_max_ar_gt_1000": fnum(sk[hi].max()) if hi.any() else None,
                        "n_ar_gt_1000_skew_gt_0p1": int((hi & (sk > 0.1)).sum())}}


def wall_area_weights(xw_m, rw_m) -> np.ndarray:
    """check_wall_resolution.element_areas と同じ: 壁の線分の長さ × 2π r_mid を両端の節点に等分。"""
    seg = np.hypot(np.diff(xw_m), np.diff(rw_m))
    a = seg * 2.0 * np.pi * 0.5 * (np.abs(rw_m[1:]) + np.abs(rw_m[:-1]))
    w = np.zeros(len(xw_m))
    w[:-1] += 0.5 * a
    w[1:] += 0.5 * a
    return w


def estimate(xr, y1_new, y1_old, y1p_max, w, target) -> dict:
    est = y1p_max * y1_new / y1_old
    g = np.isfinite(est)
    i = int(np.nanargmax(np.where(g, est, -np.inf)))
    out = {"y1plus_est_max": fnum(est[i]), "x_rt_at_max": fnum(xr[i]),
           "area_over1_pct": fnum(100 * w[g & (est > 1.0)].sum() / w[g].sum()),
           "area_over_target_pct": fnum(100 * w[g & (est > target * (1 + 1e-9))].sum() / w[g].sum()),
           "n_over_target": int((g & (est > target * (1 + 1e-9))).sum()),
           "x_lt0": {"y1plus_est_max": fnum(np.nanmax(est[xr < 0])),
                     "area_over_target_pct": fnum(100 * w[g & (xr < 0) & (est > target * (1 + 1e-9))].sum() / w[g].sum())},
           "x_ge0": {"y1plus_est_max": fnum(np.nanmax(est[xr >= 0])), "x_rt_at_max": fnum(xr[xr >= 0][int(np.nanargmax(est[xr >= 0]))]),
                     "area_over1_pct": fnum(100 * w[g & (xr >= 0) & (est > 1.0)].sum() / w[g].sum())}}
    return {"summary": out, "est": est}


def convert_and_check(work: Path, run: Path, coords, quads, bedges, scale: float, digits: int, tools: Path, keep: bool) -> dict:
    """B のビルド済みの変換器で node 変換し、品質の VERDICT・AR の例外の条件・壁の第一層の距離 (ツールの関数) を出す。"""
    from forge_design.evaluate.runner import mesh_ar_relax_check, run_converter
    from forge_design.meshing import mesh2d
    import check_wall_resolution as cwr
    import h5py
    work.mkdir(parents=True, exist_ok=True)
    for n in ("solverConfig.yaml", "bcondConfig.yaml"):
        shutil.copy2(run / n, work / n)
    mesh2d.write_msh41_2d(work / "nozzle.msh", coords, quads, bedges, digits=digits)
    out = {"workdir": str(work), "msh_sha256": sha256_file(work / "nozzle.msh"), "converter": str(CONVERTER),
           "converter_sha256": sha256_file(CONVERTER)}
    try:
        # 生産の準備 (prepare_ns) と同じ関数: 終了時の既知の GPU エラー (HDF5 は書き終えている) だけは通して隣に記録を書く
        r = run_converter([str(CONVERTER), "nozzle.msh", "nozzle.h5"], cwd=work)
        (work / "converter.log").write_text(r.stdout + r.stderr)
        out["converter_rc"] = r.returncode
        out["converter_exit_log"] = (work / "nozzle_converter_exit.log").read_text().splitlines()[0] \
            if (work / "nozzle_converter_exit.log").is_file() else None
    except subprocess.CalledProcessError as e:
        (work / "converter.log").write_text((e.stdout or "") + (e.stderr or ""))
        out["converter_rc"] = e.returncode
        out["error"] = "変換に失敗"
        cleanup(work, keep)
        return out
    with h5py.File(work / "nozzle.h5", "r") as f:
        c = f["/MESH/COORD"][:]
    hc = hashlib.sha256()
    hc.update(str(c.dtype).encode() + str(c.shape).encode() + np.ascontiguousarray(c).tobytes())
    out["coord_sha256"] = hc.hexdigest()
    q = subprocess.run([sys.executable, str(tools / "check_mesh_quality.py"), "nozzle.h5", "--ar-max", str(int(AR_MAX))],
                       cwd=str(work), capture_output=True, text=True)
    out["quality"] = {"rc": q.returncode, "lines": [ln for ln in (q.stdout + q.stderr).splitlines() if ln.strip()]}
    relax = mesh_ar_relax_check(work / "nozzle.h5", coords, bedges, scale, AR_MAX)
    out["ar_relax"] = relax
    g = subprocess.run([sys.executable, str(tools / "check_wall_resolution.py"), "nozzle.h5", "--geometry-only", "--bcond", "bcondConfig.yaml"],
                       cwd=str(work), capture_output=True, text=True)
    out["geometry_only"] = {"rc": g.returncode, "lines": [ln for ln in (g.stdout + g.stderr).splitlines() if ln.strip()]}
    # 節点ごとの y₁ (ツールの wall_first_distance、第一内部点 = 法線に最も沿った隣接 DOF)
    mesh = str(work / "nozzle.h5")
    nb = cwr.dof_neighbors(mesh)
    with h5py.File(mesh, "r") as f:
        nod = np.array(f["MESH/COORD"]).reshape(-1, 3)
    wd = cwr.wall_first_distance(mesh, 3, nb, 0.5, coords=nod)
    ids = np.array(sorted(wd))
    xy = nod[ids]
    o = np.argsort(xy[:, 0], kind="stable")
    out["_y1_tool"] = {"x": xy[o, 0], "y": xy[o, 1], "y1": np.array([wd[int(k)][0] for k in ids[o]]),
                       "ok": np.array([wd[int(k)][2] for k in ids[o]])}
    cleanup(work, keep)
    return out


def cleanup(work: Path, keep: bool) -> None:
    """変換の大きい生成物 (msh・h5) を消す (判定の記録は JSON に残す)。"""
    if keep:
        return
    for n in ("nozzle.h5", "nozzle.msh", "nozzle.xmf"):
        p = work / n
        if p.exists():
            p.unlink()


def design(a) -> int:
    sys.path.insert(0, str(ROOT / "design"))
    sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
    from forge_design.geometry.wall_axismach import load_wall_file
    case = Path(a.case_dir).resolve()
    run = case / a.run
    tools = Path(a.tools)
    out_p = Path(a.out) if a.out else HERE / OUT_DES
    env_p = Path(a.envelope) if a.envelope else HERE / OUT_ENV
    R = {"plan": PLAN, "created": now(), "kind": "表の案と作り直した格子 (§5.1 #36r の前半)。CFD 0 step・NS は回さない・結果の解釈はしない",
         "run": a.run, "run_dir": str(run), "envelope": str(env_p), "script": git_head(), "errors": []}
    err = R["errors"]
    for p in (env_p, run / "wall_repr.json", run / "prepare_info.json", run / "bcondConfig.yaml", run / "solverConfig.yaml", case / a.problem):
        if not p.is_file():
            err.append(f"{p} が無い")
    if err:
        jwrite(out_p, R)
        return 2
    E = json.loads(env_p.read_text())
    info = json.loads((run / "prepare_info.json").read_text())
    scale = float(info["scale_m"])
    block = info["mesh"]["block"]
    prob = yaml.safe_load((case / a.problem).read_text())
    pm = prob.get("mesh") or {}
    psha = sha256_file(case / a.problem)
    R["inputs"] = {"envelope_sha256": sha256_file(env_p), "wall_repr_sha256": sha256_file(run / "wall_repr.json"),
                   "prepare_info_sha256": sha256_file(run / "prepare_info.json"), "problem_yaml": a.problem, "problem_sha256": psha,
                   "problem_is_run_input": psha == PROBLEM_SHA_RUN,
                   "problem_mesh_equals_prepare_info_block": all(pm.get(k) == block.get(k) for k in set(pm) | set(block)),
                   "r_t_m": scale}
    if not R["inputs"]["problem_mesh_equals_prepare_info_block"]:
        err.append("問題 YAML の mesh 節が run の prepare_info.mesh.block と違う")
    W = load_wall_file(run / "wall_repr.json")
    wall = W["physical"]
    # --- 対照: 元の表で作り直す ---
    prm0, c0, q0, b0 = mesh_from(block, scale, wall)
    gsha = generator_sha(c0, q0, b0)
    M0 = mesh_metrics(prm0, c0, q0, scale)
    xr = np.asarray(E["nodes"]["x_m"]) / scale
    xw_m, rw_m = M0["X"][:, -1], M0["Y"][:, -1]
    y1_old = np.asarray(E["nodes"]["y1_m"], dtype=float)
    y1p_max = np.asarray(E["nodes"]["y1plus_max"], dtype=float)
    w = np.asarray(E["nodes"]["area_weight_m2"], dtype=float)
    ctrl = {"generator_sha256": gsha, "expected": info["mesh"]["hashes"]["generator_sha256"],
            "generator_sha_match": gsha == info["mesh"]["hashes"]["generator_sha256"],
            "n_wall_nodes_mesh": int(len(xw_m)), "n_wall_nodes_envelope": int(len(xr))}
    if len(xw_m) == len(xr):
        ctrl["wall_node_x_max_abs_diff_m"] = fnum(np.max(np.abs(xw_m - np.asarray(E["nodes"]["x_m"]))))
        ctrl["wall_node_r_max_abs_diff_m"] = fnum(np.max(np.abs(rw_m - np.asarray(E["nodes"]["r_m"]))))
        ctrl["y1_generator_over_tool"] = [fnum(np.min(M0["y1"] / y1_old)), fnum(np.max(M0["y1"] / y1_old))]
        ctrl["area_weight_rel_diff_max"] = fnum(np.max(np.abs(wall_area_weights(xw_m, rw_m) / w - 1.0)))
    ctrl["metrics"] = M0["summary"]
    R["control"] = ctrl
    if not ctrl["generator_sha_match"]:
        err.append("対照: 生成器の sha256 が run の記録と違う (作り直しが run の格子を再現しない)")
    if len(xw_m) != len(xr) or ctrl.get("wall_node_x_max_abs_diff_m", 1.0) > 1e-9:
        err.append("対照: 壁の節点が包絡の節点と対応しない")
    if err:
        jwrite(out_p, R)
        return 2
    old_tbl = np.asarray(block["wall_first_frac_table"], dtype=float)
    fr_old = M0["fr"]
    xs_w = M0["xw"]
    ctrl["envelope_estimate"] = estimate(xr, y1_old, y1_old, y1p_max, w, 0.8)["summary"]
    if a.convert:
        cc = convert_and_check(Path(a.workdir) / "control", run, c0, q0, b0, scale, int(block.get("msh_digits", 10)), tools, keep=False)
        cc["coord_sha256_expected"] = info["mesh"]["hashes"]["coord_sha256"]
        cc["coord_sha_match"] = cc.get("coord_sha256") == cc["coord_sha256_expected"]
        y1t = cc.pop("_y1_tool", None)
        if y1t is not None:
            cc["y1_tool_over_envelope_y1"] = [fnum(np.min(y1t["y1"] / y1_old)), fnum(np.max(y1t["y1"] / y1_old))]
        R["control"]["convert"] = cc
    plots = {"control": M0}
    R["variants"] = []
    # 一次の見積もりの係数: 逆粘性長さ L = y₁⁺_max / y₁_old [1/m] (節点の動く案は x で線形補間)
    Linv = y1p_max / y1_old
    done = {}
    for v in VARIANTS if not a.only else [v for v in VARIANTS if v["name"] in a.only.split(",")]:
        V = {"spec": v}
        blk = block
        if v.get("x_refine"):
            base = done.get(v["x_refine"]["base"])
            if base is None:
                V["error"] = f"基の案 {v['x_refine']['base']} が無い"
                R["variants"].append(V)
                continue
            xinfo = x_refine_density(block, wall, base, float(v["x_refine"]["ar_goal"]), M0["ar2"][:, -1])
            V["x_refine"] = {k: xinfo[k] for k in xinfo if k != "table"}
            blk = {**block, "x_density_table": xinfo["table"], "ni": xinfo["ni"]}
            prm_s, cs, qs, bs = mesh_from(blk, scale, wall, base["table"])          # station の位置だけを取る (表は station に効かない)
            xs_v = cs[:, 0].reshape(prm_s.ni, -1)[:, -1] / scale
            from forge_design.meshing import mesh2d
            fr_old_v = mesh2d._first_frac_profile(xs_v, prm0)
            # 新しい station の要求: 見積もり (y₁_new · L(x)、L は x の線形補間) と同じ係数で書く: y₁⁺_max 相当 = L(x) · y₁_old 相当 (fr_old · r_w)
            y1p_v = np.interp(xs_v, xr, Linv) * fr_old_v * np.asarray(wall.r(xs_v), dtype=float) * scale
            T = design_table(old_tbl, xs_v, fr_old_v, y1p_v, v["target"], v["spacing"], v["thicken"])
        else:
            T = design_table(old_tbl, xs_w, fr_old, y1p_max, v["target"], v["spacing"], v["thicken"])
        V["table"] = {k: T[k] for k in T if k != "table"}
        if not T["ok"]:
            V["error"] = "線形計画が解けない"
            R["variants"].append(V)
            continue
        tbl = T["table"]
        V["table"]["x_lt0"] = [p for p in tbl if p[0] < 0.0]
        V["table"]["x_ge0_unchanged"] = [p for p in tbl if p[0] >= 0.0] == [list(map(float, p)) for p in old_tbl if p[0] >= 0.0]
        V["table"]["old_x_lt0"] = [[float(p[0]), float(p[1])] for p in old_tbl if p[0] < 0.0]
        prm, c1, q1, b1 = mesh_from(blk, scale, wall, tbl)
        M1 = mesh_metrics(prm, c1, q1, scale)
        xr1 = M1["xw"]
        same_nodes = bool(len(xr1) == len(xr) and np.array_equal(M1["X"][:, -1], M0["X"][:, -1]) and np.array_equal(M1["Y"][:, -1], M0["Y"][:, -1]))
        V["wall_nodes_bit_identical"] = same_nodes
        if same_nodes:
            L1, w1, y1o = Linv, w, y1_old
            cols = M1["xw"] >= 0.0
            V["downstream_bit_identical"] = bool(np.array_equal(M1["X"][cols], M0["X"][cols]) and np.array_equal(M1["Y"][cols], M0["Y"][cols]))
        else:
            L1 = np.interp(xr1, xr, Linv)
            w1 = wall_area_weights(M1["X"][:, -1], M1["Y"][:, -1])
            y1o = np.interp(xr1, xr, y1_old)
            sh = len(xr1) - len(xr)                       # 列の数の増分 (上流に足した列)
            cols = M0["xw"] >= 0.0
            V["downstream_shifted_columns"] = {"shift": int(sh),
                                               "x_max_abs_diff_m": fnum(np.max(np.abs(M1["X"][sh:][cols] - M0["X"][cols]))),
                                               "r_max_abs_diff_m": fnum(np.max(np.abs(M1["Y"][sh:][cols] - M0["Y"][cols]))),
                                               "bit_identical": bool(np.array_equal(M1["X"][sh:][cols], M0["X"][cols])
                                                                     and np.array_equal(M1["Y"][sh:][cols], M0["Y"][cols])),
                                               "note": "x ≥ 0 の列を、上流に足した列の数だけずらして元の格子と比べた"}
        V["metrics"] = M1["summary"]
        V["metrics"]["nodes_change"] = int(M1["summary"]["nodes"] - M0["summary"]["nodes"])
        up = M1["xw"][:-1] < 0.0
        V["metrics"]["x_lt0"] = {"ar_max": fnum(M1["ar2"][up].max()), "n_ar_gt_5000": int((M1["ar2"][up] > AR_MAX).sum()),
                                 "x_rt_ar_gt_5000": ([fnum(M1["xw"][:-1][up][(M1["ar2"][up] > AR_MAX).any(1)].min()),
                                                      fnum(M1["xw"][:-1][up][(M1["ar2"][up] > AR_MAX).any(1)].max())]
                                                     if (M1["ar2"][up] > AR_MAX).any() else None),
                                 "ar_max_control": fnum(M0["ar2"][M0["xw"][:-1] < 0.0].max())}
        if M1["ratio"] is not None:
            V["metrics"]["radial_ratio_max"] = {"x_lt0": fnum(M1["ratio"][M1["xw"] < 0].max()), "control_x_lt0": fnum(M0["ratio"][M0["xw"] < 0].max()),
                                                "all": fnum(M1["ratio"].max()), "control_all": fnum(M0["ratio"].max())}
        est = estimate(xr1, M1["y1"], 1.0, L1, w1, v["target"])
        V["estimate"] = est["summary"]
        # 位置の区切りごとの y₁ (生成の座標) と見積もり
        bins = []
        for lo, hi, m in bins_of(xr1):
            m0 = next(mm for l0, h0, mm in bins_of(xr) if l0 == lo)
            bins.append({"x_rt": [lo, hi], "n_nodes": int(m.sum()), "n_nodes_control": int(m0.sum()),
                         "y1_new_m": [fnum(M1["y1"][m].min()), fnum(M1["y1"][m].max())],
                         "y1_old_m": [fnum(y1o[m].min()), fnum(y1o[m].max())],
                         "y1plus_est_max": fnum(est["est"][m].max()), "y1plus_env_max": fnum(y1p_max[m0].max()),
                         "ar_wall_cell_max": fnum(M1["ar2"][:, -1][m[:-1]].max()) if m[:-1].any() else None,
                         "ar_wall_cell_max_control": fnum(M0["ar2"][:, -1][m0[:-1]].max()) if m0[:-1].any() else None})
        V["bins"] = bins
        # AR ≤ 5000 に要る流れ方向の間隔 (壁のセル: 最長辺 = 流れ方向の辺、最短辺 = 第一層)
        dxw = np.hypot(np.diff(M1["X"][:, -1]), np.diff(M1["Y"][:, -1]))
        need = dxw / (AR_MAX * np.minimum(M1["y1"][:-1], M1["y1"][1:]))
        V["x_refine_needed"] = {"max_factor": fnum(need.max()), "x_rt_at_max": fnum(M1["xw"][:-1][int(np.argmax(need))]),
                                "extra_columns_estimate": int(np.ceil(np.sum(np.maximum(need, 1.0) - 1.0))),
                                "x_rt_range_factor_gt_1": ([fnum(M1["xw"][:-1][need > 1].min()), fnum(M1["xw"][:-1][need > 1].max())]
                                                           if (need > 1).any() else None),
                                "note": "壁のセルの AR ≒ 流れ方向の辺 / 第一層。流れ方向の間隔をこの倍率で細かくすれば AR ≤ 5000 (列の数は Σ (倍率 − 1) の見積もり)"}
        if a.convert:
            cc = convert_and_check(Path(a.workdir) / v["name"], run, c1, q1, b1, scale, int(block.get("msh_digits", 10)), tools, keep=False)
            y1t = cc.pop("_y1_tool", None)
            if y1t is not None and len(y1t["y1"]) == len(xr1):
                cc["y1_tool_vs_generator_rel_max"] = fnum(np.max(np.abs(y1t["y1"] / M1["y1"] - 1.0)))
                cc["y1_tool_eval_frac"] = fnum(np.mean(y1t["ok"]))
                e2 = estimate(xr1, y1t["y1"], 1.0, L1, w1, v["target"])
                cc["estimate_tool_y1"] = e2["summary"]
            V["convert"] = cc
        if v.get("x_refine"):
            V["x_density_table_x_lt0"] = [p for p in blk["x_density_table"] if p[0] < 0.0]
        R["variants"].append(V)
        done[v["name"]] = {"table": tbl, "M": M1}
        plots[v["name"]] = {**M1, "est": est["est"], "table": tbl}
        print(f"{v['name']}: 変えた点 {T['n_changed']}/{T['n_var']}、節点 {M1['summary']['nodes']} ({V['metrics']['nodes_change']:+d})、"
              f"AR 最大 {M1['summary']['ar_max']} (x<0 で 5000 超え {V['metrics']['x_lt0']['n_ar_gt_5000']} セル)、"
              f"見積もり y₁⁺ 最大 {est['summary']['y1plus_est_max']} 1 超え {est['summary']['area_over1_pct']} %、"
              f"AR ≤ 5000 に要る流れ方向の細分 ×{V['x_refine_needed']['max_factor']}"
              + (f"、変換後 {[ln for ln in V['convert'].get('quality', {}).get('lines', []) if ln.startswith('VERDICT')]}" if a.convert else ""),
              flush=True)
    png = out_p.with_suffix(".png")
    try:
        plot_design(png, xr, y1_old, y1p_max, plots, old_tbl)
        R["png"] = png.name
    except Exception as e:                      # noqa: BLE001
        R["png_error"] = repr(e)
    jwrite(out_p, R)
    print(f"→ {out_p}")
    return 0 if not err else 2


def plot_design(png: Path, xr, y1_old, y1p_max, plots: dict, old_tbl) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 1, figsize=(11, 14))
    sl = xr < 0.5
    tsl = old_tbl[:, 0] < 0.5
    ax[0].semilogy(old_tbl[tsl, 0], old_tbl[tsl, 1], "ko-", ms=3, lw=0.8, label="table (run_0951)")
    for k, P in plots.items():
        if k == "control":
            continue
        t = np.asarray(P["table"])
        ts = t[:, 0] < 0.5
        ax[0].semilogy(t[ts, 0], t[ts, 1], ".-", ms=3, lw=0.8, label=f"table {k}")
    ax[0].set_ylabel("wall_first_frac_table (/ r_w)")
    ax[1].semilogy(xr[sl], y1_old[sl] * 1e6, "k-", lw=1.0, label="y1 run_0951")
    for k, P in plots.items():
        if k != "control":
            s1 = P["xw"] < 0.5
            ax[1].semilogy(P["xw"][s1], P["y1"][s1] * 1e6, lw=0.8, label=f"y1 {k}")
    ax[1].set_ylabel("y1 [um]")
    ax[2].plot(xr[sl], y1p_max[sl], "k-", lw=1.0, label="y1+ envelope max (run_0951)")
    for k, P in plots.items():
        if k != "control":
            s1 = P["xw"] < 0.5
            ax[2].plot(P["xw"][s1], P["est"][s1], lw=0.8, label=f"estimate {k}")
    ax[2].axhline(1.0, color="r", lw=0.8, ls="--")
    ax[2].axhline(0.8, color="k", lw=0.8, ls=":")
    ax[2].set_ylabel("y1+ (flow fixed at envelope max)")
    for k, P in plots.items():
        xc = 0.5 * (P["xw"][1:] + P["xw"][:-1])
        s2 = xc < 0.5
        ax[3].semilogy(xc[s2], P["ar2"][s2, -1], lw=0.8 if k != "control" else 1.2, color="k" if k == "control" else None, label=f"AR wall cell {k}")
    ax[3].axhline(AR_MAX, color="r", lw=0.8, ls="--", label="5000")
    ax[3].set_ylabel("AR (wall cell)")
    for a_ in ax:
        a_.set_xlabel("x / r_t")
        a_.grid(alpha=0.3, which="both")
        a_.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=("envelope", "design"))
    ap.add_argument("--case-dir", default=str(Path.home() / "forge-wallfit/case/45.isobutane_m6_d155"))
    ap.add_argument("--run", default=RUN)
    ap.add_argument("--steps", default=None, help="カンマ区切り (既定は保持の窓の 17 出力)")
    ap.add_argument("--tools", default=str(ROOT / "solver_density_cuda" / "tools"))
    ap.add_argument("--problem", default=PROBLEM)
    ap.add_argument("--envelope", default=None)
    ap.add_argument("--convert", action="store_true", help="B のビルド済みの変換器で変換して品質・第一層の距離を確かめる")
    ap.add_argument("--workdir", default=str(Path.home() / "forge-g0cchain/_tmp_36r"))
    ap.add_argument("--only", default=None, help="案の名前 (カンマ区切り)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    return envelope(a) if a.phase == "envelope" else design(a)


if __name__ == "__main__":
    sys.exit(main())
