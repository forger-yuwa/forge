"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36o: V-h370-diag (run_0951_h370_ns_b0) の停止の出力 152,500 の壁の解像の切り分け
(codex diagnose `notes/reviews/2026-10-14-h370-diag-result-diagnose.md` の Major 3、`notes/reviews/2026-10-14-h370-fixedband-result-diagnose.md`
の Major「新しい計算の前に既存の出力で壁の解像を点検」)。CFD 0 step。格子は直さない。結果の解釈はしない (事実だけを JSON に書く)。

1. 境界の分類 (出口の群がなぜ壁の解像の対象に入るか): run の bcondConfig.yaml の各境界の種類 (kind)・outputHDFflg、停止の出力の境界の出力
   (`res_<群>_<physID>_<step>.h5`) の有無と VALUE の量 (twall_* の有無)、`check_wall_resolution.py` を --groups なし (受け渡しの関門
   `deltastar_loop.wall_resolution_verdict` と同じ呼び方) と --groups wall で回した各群の行・VERDICT・終了コード、それを決めている
   ソースの行 (ファイル:行、実行時に文字列で探す; 見つからなければ None)。
2. 壁の群 (--groups wall --weight area --profile-csv) の位置ごとの y₁ (第一内部節点までの壁法線の距離)・接線の壁応力 τ_t・ρ_w・μ_w・y₁⁺・
   面積の重み (ツールの出力そのもの; 式を書き写さない)。y₁⁺ > 1 の連続した区間 (x の昇順) ごとに x の範囲・節点の数・面積と割合・最大と位置・
   形 (直管 x ≤ x_contraction / 縮流部 / x ≥ 0) の内訳。**近くの y₁⁺ ≤ 1 の場所** = 区間の端に接する y₁⁺ ≤ 1 の節点のうち、端から x 方向に
   FLANK_RT [r_t] 以内で、y₁⁺ > 1 の節点に当たるまでの連続した節点 (上流側・下流側。入口の端に接する側は空)。区間と近くの場所の比を、
   面積の重みつきの幾何平均 (log の重みつき平均の exp) で y₁・τ_t・ρ_w・μ_w・y₁⁺ について出す。y₁⁺ = y₁·√(ρ_w τ_t)/μ_w なので log の和に
   厳密に分かれる: log10(y₁⁺ の比) = log10(y₁ の比) + ½ log10(ρ の比) + ½ log10(τ の比) − log10(μ の比) (各項を桁で記録)。最大の点と
   近くの場所の最初の節点 (端に接する節点) の比も同じ形で記録する。区間の外も含め、位置の区切り (bins) ごとの代表値も残す。
3. 設計の第一層: run の prepare_info.json の mesh.block (問題 YAML の mesh 節と照合) から格子の生成と同じ関数 (`runner_axismach._mesh_params_from`・
   `mesh2d._first_frac_profile`・`mesh2d._radial_fracs_capfixed` / `_radial_fracs`) で各壁の節点の station の第一セル厚 d₁ = r_w (1 − s_{nj−2})
   [m] を出す (wall_normal_layer があれば第一内部節点は壁法線上の距離 d₁ に置かれる = 壁法線の距離の設計値)。実際の y₁ との比を全域と区間で。
4. 補助 (記録だけ): 形 C・370 K の run (既定 run_0813_tw370_g0c_c50 の最後の壁の出力) の同じ x での y₁・τ_t・ρ_w・μ_w・y₁⁺ (線形補間) と
   H/C の比、両方の設計の第一層の実寸 (表の比 × 局所半径) の比。C の run の格子・収束の状態は違う (g0c、nj 121) ので比較は補助の記録。

使い方 (B の case dir で、system の python3):
  python3 h370_wallres_split.py [--run run_0951_h370_ns_b0] [--step 152500] [--aux-run run_0813_tw370_g0c_c50] [--aux-step 97500]
         [--design ~/forge-g0cchain/design] [--tools ~/forge-g0cchain/solver_density_cuda/tools] [--out _band_ab/half_inlet/h370_wallres_split.json]
終了コード 0 書けた / 2 入力が足りない (足りないものは JSON の errors)。
"""
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
PLAN = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36o (codex diagnose notes/reviews/2026-10-14-h370-diag-result-diagnose.md "
        "の Major 3・notes/reviews/2026-10-14-h370-fixedband-result-diagnose.md の Major)")
RUN, STEP = "run_0951_h370_ns_b0", 152500
AUX_RUN, AUX_STEP = "run_0813_tw370_g0c_c50", 97500
PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
OUT = "_band_ab/half_inlet/h370_wallres_split.json"
TARGET = 1.0            # 局所 y₁⁺ の目標 (check_wall_resolution.py の既定と同じ)
FLANK_RT = 1.0          # 近くの y₁⁺ ≤ 1 の場所の窓 [r_t] (結果を見る前に決めた値)
BINS_RT = [-12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.0, -2.0,
           -1.0, 0.0, 1.0, 3.0, 5.0, 10.0, 40.0, 94.0, None]   # 位置の区切り [r_t] (最後は出口まで)
DUMP_RE = re.compile(r"res_(.+)_(\d+)_(\d+)\.h5$")
# 出口の群が判定に入る理由を決めているソースの行 (実行時に文字列で探す)
SOURCE_PATTERNS = [
    ("tools", "check_wall_resolution.py", 'dumps = sorted(glob.glob(os.path.join(a.run, "res_*_*_*.h5")))',
     "main: 境界の出力 res_<群>_<physID>_<step>.h5 を種類によらずすべて拾う"),
    ("tools", "check_wall_resolution.py", "if want and name not in want:", "main: 群を絞るのは --groups を渡したときだけ"),
    ("tools", "check_wall_resolution.py", 'need = ("twall_x", "twall_y", "twall_z", "ro", "Ts")', "main: twall_* が無い群は判定不能にする"),
    ("tools", "check_wall_resolution.py", 'str(d.get("kind", "")).startswith("wall")', "--geometry-only だけは bcondConfig の kind が wall* の境界に絞る"),
    ("design", "forge_design/feedback/deltastar_loop.py", 'cmd = [sys.executable, str(FORGE_TOOLS / "check_wall_resolution.py"), str(run_dir)]',
     "受け渡しの関門 wall_resolution_verdict は --groups を渡さない"),
    ("design", "forge_design/report/nozzle_report.py", 'cmd = [sys.executable, str(TOOLS / "check_wall_resolution.py"), str(run), "--groups", ",".join(walls)]',
     "標準出力の報告 (nozzle_report.wall_resolution) は kind が wall* の境界だけを --groups で渡す"),
    ("design", "forge_design/evaluate/runner_wt.py", "outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1",
     "bcondConfig の出口の行は outputHDFflg 1 (境界の出力を書く)"),
    ("design", "forge_design/evaluate/runner_axismach.py", '_bcond(p, euler=False, profile_ok=True)', "NS の準備 (prepare_ns) は runner_wt._bcond で bcondConfig を書く"),
]


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def fnum(v):
    """JSON に書ける数 (非有限は None)。"""
    if v is None:
        return None
    v = float(v)
    return v if math.isfinite(v) else None


def source_facts(tools: Path, design: Path) -> list:
    out = []
    for root, rel, pat, what in SOURCE_PATTERNS:
        p = (tools if root == "tools" else design) / rel
        hit = None
        if p.is_file():
            for i, ln in enumerate(p.read_text().splitlines(), 1):
                if pat in ln:
                    hit = {"line": i, "text": ln.strip()[:240]}
                    break
        out.append({"file": f"{'solver_density_cuda/tools' if root == 'tools' else 'design'}/{rel}", "what": what,
                    **(hit or {"line": None, "text": None, "note": "文字列が見つからない"})})
    return out


def bcond_table(run: Path) -> list:
    d = yaml.safe_load((run / "bcondConfig.yaml").read_text()) or {}
    return [{"name": k, "physID": int(v.get("physID")), "kind": v.get("kind"), "outputHDFflg": int(v.get("outputHDFflg") or 0)}
            for k, v in d.items() if isinstance(v, dict)]


def boundary_dumps(run: Path, step: int) -> list:
    import h5py
    out = []
    for p in sorted(run.glob(f"res_*_*_{step}.h5")):
        m = DUMP_RE.search(p.name)
        if not m or int(m.group(3)) != step:
            continue
        with h5py.File(p, "r") as f:
            keys = sorted(f["VALUE"].keys()) if "VALUE" in f else []
            n = int(np.asarray(f["MESH/COORD"]).size // 3) if "MESH/COORD" in f else None
        out.append({"file": p.name, "group": m.group(1), "physID": int(m.group(2)), "n_points": n, "value_keys": keys,
                    "has_twall": all(k in keys for k in ("twall_x", "twall_y", "twall_z")), "sha256": sha256_file(p)})
    return out


def run_tool(tools: Path, case: Path, run: str, step: int, groups=None, weight="area", profile=None) -> dict:
    cmd = [sys.executable, str(tools / "check_wall_resolution.py"), run, "--step", str(step), "--weight", weight]
    if groups:
        cmd += ["--groups", groups]
    if profile:
        cmd += ["--profile-csv", str(profile)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(case))
    out = r.stdout + r.stderr
    lines = [ln.rstrip() for ln in out.splitlines() if ln.strip()]
    return {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdict": [ln.strip() for ln in lines if ln.startswith("VERDICT")],
            "lines": [ln for ln in lines if "分布 ->" not in ln][:40]}


def read_profile(p: Path) -> dict:
    d = np.genfromtxt(p, delimiter=",", names=True)
    o = np.argsort(d["x"], kind="stable")
    return {k: np.asarray(d[k], dtype=float)[o] for k in d.dtype.names}


def design_first_layer(block: dict, scale: float, x_m, r_m) -> dict:
    """格子の生成と同じ関数で、壁の節点 (x_m, r_m) [m] の station の第一セル厚 d₁ [m] (壁法線の層のとき壁法線の距離) を出す。"""
    from forge_design.evaluate.runner_axismach import _mesh_params_from
    from forge_design.meshing import mesh2d
    prm = _mesh_params_from(block, scale, int(block.get("ni", 0)), int(block.get("nj", 0)), float(block.get("wall_first_frac", 0.0)))
    xs = np.asarray(x_m, dtype=float) / scale
    rw = np.asarray(r_m, dtype=float) / scale
    nj = prm.nj
    if prm.axis_gap_frac is not None:
        raise ValueError("axis_gap_frac の格子は対象外 (この台本は axis_cap_frac か表・一様の格子だけ)")
    fr = mesh2d._first_frac_profile(xs, prm)
    g = np.empty(len(xs))
    for i in range(len(xs)):
        if prm.axis_cap_frac is not None:
            s, _ = mesh2d._radial_fracs_capfixed(nj, float(fr[i]), float(prm.axis_cap_frac))
        else:
            s = mesh2d._radial_fracs(nj, float(fr[i]))
        g[i] = 1.0 - s[-2]
    d1 = rw * g * scale
    return {"frac_table": fr, "gap_frac": g, "d1_m": d1, "normal_layer": prm.wall_normal_layer is not None,
            "branch": "axis_cap_frac (_radial_fracs_capfixed)" if prm.axis_cap_frac is not None else "_radial_fracs",
            "abs_table_m": fr * rw * scale}


def gmean(v, w, m):
    """面積の重みつきの幾何平均 (m の中、正の有限値)。"""
    m = m & np.isfinite(v) & (v > 0)
    if not m.any():
        return None
    return float(np.exp(np.sum(w[m] * np.log(v[m])) / np.sum(w[m])))


def ratio_block(P: dict, mreg, mref) -> dict:
    """区間 mreg と近くの場所 mref の幾何平均の比と、log10(y₁⁺ の比) の厳密な分け方。"""
    w = P["area_weight"]
    q = {"y1_m": "y1", "tau_t": "tau_t", "rho_w": "rho_w", "mu_w": "mu_w", "y1plus": "y1plus"}
    gr = {name: gmean(P[k], w, mreg) for k, name in q.items()}
    gf = {name: gmean(P[k], w, mref) for k, name in q.items()}
    if not mref.any() or any(v is None for v in list(gr.values()) + list(gf.values())):
        return {"region_gmean": gr, "flank_gmean": gf, "ratio": None, "note": "近くの場所が空か非有限"}
    rt = {k: gr[k] / gf[k] for k in gr}
    dec = {"y1": math.log10(rt["y1"]), "tau_t (½)": 0.5 * math.log10(rt["tau_t"]), "rho_w (½)": 0.5 * math.log10(rt["rho_w"]),
           "mu_w (−1)": -math.log10(rt["mu_w"])}
    return {"region_gmean": gr, "flank_gmean": gf, "ratio_region_over_flank": rt, "log10_y1plus_ratio": math.log10(rt["y1plus"]),
            "log10_terms": dec, "log10_terms_sum": sum(dec.values()),
            "larger_contribution": max(dec, key=lambda k: abs(dec[k]))}


def point_ratio(P: dict, i: int, j: int) -> dict:
    """節点 i と節点 j の比 (i / j) と log10 の分け方。"""
    r = {k: P[k][i] / P[k][j] for k in ("y1_m", "tau_t", "rho_w", "mu_w", "y1plus")}
    dec = {"y1": math.log10(r["y1_m"]), "tau_t (½)": 0.5 * math.log10(r["tau_t"]), "rho_w (½)": 0.5 * math.log10(r["rho_w"]),
           "mu_w (−1)": -math.log10(r["mu_w"])}
    return {"ratio": {k: fnum(v) for k, v in r.items()}, "log10_terms": dec, "log10_y1plus_ratio": math.log10(r["y1plus"])}


def node_rec(P: dict, i: int, rt: float, D=None) -> dict:
    r = {"i": int(i), "x_m": fnum(P["x"][i]), "x_rt": fnum(P["x"][i] / rt), "r_rt": fnum(P["y"][i] / rt), "y1_m": fnum(P["y1_m"][i]),
         "y1plus": fnum(P["y1plus"][i]), "rho_w": fnum(P["rho_w"][i]), "tau_t": fnum(P["tau_t"][i]), "mu_w": fnum(P["mu_w"][i]),
         "T_w": fnum(P["T_w"][i]), "area_weight_m2": fnum(P["area_weight"][i])}
    if D is not None:
        r.update(design_d1_m=fnum(D["d1_m"][i]), y1_over_design=fnum(P["y1_m"][i] / D["d1_m"][i]), frac_table=fnum(D["frac_table"][i]))
    return r


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


def flank(P: dict, a: int, b: int, rt: float, side: str):
    """区間 [a, b] (x の昇順の添字) の端に接する y₁⁺ ≤ 1 の節点のうち、端から FLANK_RT 以内で y₁⁺ > 1 に当たるまでの連続した節点。"""
    x, yp, n = P["x"] / rt, P["y1plus"], len(P["x"])
    out = []
    if side == "down":
        k = b + 1
        while k < n and np.isfinite(yp[k]) and yp[k] <= TARGET and x[k] - x[b] <= FLANK_RT:
            out.append(k)
            k += 1
    else:
        k = a - 1
        while k >= 0 and np.isfinite(yp[k]) and yp[k] <= TARGET and x[a] - x[k] <= FLANK_RT:
            out.append(k)
            k -= 1
    return out


def bins_table(P: dict, rt: float, D=None, aux=None) -> list:
    x = P["x"] / rt
    w = P["area_weight"]
    out = []
    for lo, hi in zip(BINS_RT[:-1], BINS_RT[1:]):
        m = (x >= lo) & ((x < hi) if hi is not None else True)
        if not m.any():
            continue
        r = {"x_rt": [lo, hi if hi is not None else fnum(x[m].max())], "n": int(m.sum()), "area_m2": fnum(w[m].sum()),
             "y1plus_max": fnum(np.nanmax(P["y1plus"][m])), "n_over": int((P["y1plus"][m] > TARGET).sum()),
             "area_over_m2": fnum(w[m & (P["y1plus"] > TARGET)].sum()),
             "gmean": {k: gmean(P[k], w, m) for k in ("y1_m", "tau_t", "rho_w", "mu_w", "y1plus")}}
        if D is not None:
            r["gmean"]["design_d1_m"] = gmean(D["d1_m"], w, m)
            r["y1_over_design_minmax"] = [fnum(np.min(P["y1_m"][m] / D["d1_m"][m])), fnum(np.max(P["y1_m"][m] / D["d1_m"][m]))]
        if aux is not None:
            r["aux_C_gmean"] = {k: gmean(aux[k], w, m) for k in ("y1_m", "tau_t", "rho_w", "mu_w", "y1plus")}
        out.append(r)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", default=str(HERE))
    ap.add_argument("--run", default=RUN)
    ap.add_argument("--step", type=int, default=STEP)
    ap.add_argument("--aux-run", default=AUX_RUN, help="補助の形 C・370 K の run (空文字で回さない)")
    ap.add_argument("--aux-step", type=int, default=AUX_STEP)
    ap.add_argument("--problem", default=PROBLEM)
    ap.add_argument("--design", default=str(Path.home() / "forge-g0cchain/design"))
    ap.add_argument("--tools", default=str(Path.home() / "forge-g0cchain/solver_density_cuda/tools"))
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case = Path(a.case_dir).resolve()
    design, tools = Path(a.design), Path(a.tools)
    if not design.is_dir():
        design = HERE.parents[1] / "design"
    if not tools.is_dir():
        tools = HERE.parents[1] / "solver_density_cuda/tools"
    sys.path.insert(0, str(design))
    out_p = Path(a.out) if a.out else case / OUT
    run = case / a.run
    R = {"plan": PLAN, "created": now(), "kind": "診断 (§5.1 #36o)。格子は直さない・結果の解釈はしない (事実だけ)", "run": a.run, "step": a.step,
         "target_y1plus": TARGET, "flank_rt": FLANK_RT, "errors": []}
    err = R["errors"]
    wall_dump = run / f"res_wall_3_{a.step}.h5"
    for p in (run / "bcondConfig.yaml", run / "prepare_info.json", run / "nozzle.h5", run / f"res_{a.step}.h5", wall_dump):
        if not p.is_file():
            err.append(f"{p.relative_to(case)} が無い")
    if err:
        jwrite(out_p, R)
        return 2
    info = jload(run / "prepare_info.json") or {}
    rt = float(info["scale_m"])
    R["inputs"] = {"r_t_m": rt, "sha256": {n: sha256_file(run / n) for n in ("bcondConfig.yaml", "nozzle.h5", f"res_{a.step}.h5", wall_dump.name,
                                                                             "solverConfig.yaml", "prepare_info.json")},
                   "tools_dir": str(tools), "design_dir": str(design),
                   "tool_sha256": sha256_file(tools / "check_wall_resolution.py")}
    # --- 1. 境界の分類 ---
    bc = bcond_table(run)
    dumps = boundary_dumps(run, a.step)
    t_all = run_tool(tools, case, a.run, a.step)
    with tempfile.TemporaryDirectory(prefix="h370_wallres_") as td:
        t_wall = run_tool(tools, case, a.run, a.step, groups="wall", profile=Path(td) / "prof.csv")
        prof_p = Path(td) / "prof_wall.csv"
        P = read_profile(prof_p) if prof_p.is_file() else None
        aux_P, aux_tool = None, None
        if a.aux_run and (case / a.aux_run).is_dir():
            aux_tool = run_tool(tools, case, a.aux_run, a.aux_step, groups="wall", profile=Path(td) / "aux.csv")
            ap_ = Path(td) / "aux_wall.csv"
            aux_P = read_profile(ap_) if ap_.is_file() else None
    kinds = {b["physID"]: b for b in bc}
    facts = []
    for d in dumps:
        b = kinds.get(d["physID"]) or {}
        facts.append(f"境界の出力 {d['file']}: bcondConfig の {b.get('name')} (kind {b.get('kind')}, outputHDFflg {b.get('outputHDFflg')})、"
                     f"twall_* {'あり' if d['has_twall'] else 'なし'}")
    no_dump = [b for b in bc if b["outputHDFflg"] == 0]
    facts.append("outputHDFflg 0 で境界の出力の無い境界: " + ", ".join(f"{b['name']} ({b['kind']})" for b in no_dump))
    R["boundary_classification"] = {"bcond": bc, "dumps_at_step": dumps, "source": source_facts(tools, design),
                                    "tool_all_groups": t_all, "tool_wall_only": t_wall, "facts": facts}
    if P is None:
        err.append("--groups wall の分布 (profile CSV) を作れない")
        jwrite(out_p, R)
        return 2
    # --- 3. 設計の第一層 (2 の表にも使う) ---
    block = (info.get("mesh") or {}).get("block") or {}
    prob = yaml.safe_load((case / a.problem).read_text()) if (case / a.problem).is_file() else None
    pm = (prob or {}).get("mesh") or {}
    same_block = bool(prob) and all(pm.get(k) == block.get(k) for k in set(pm) | set(block))
    D = design_first_layer(block, rt, P["x"], P["y"])
    pu = ((info.get("physical_upstream") or {}).get("wall") or {})
    x_con = pu.get("x_contraction")
    # --- 2. 壁の群 ---
    x = P["x"] / rt
    w = P["area_weight"]
    yp = P["y1plus"]
    good = np.isfinite(yp)
    over = good & (yp > TARGET)
    imax = int(np.nanargmax(np.where(good, yp, -np.inf)))
    regions = []
    for a_, b_ in runs_of(over):
        m = np.zeros(len(x), dtype=bool)
        m[a_:b_ + 1] = True
        fl_d, fl_u = flank(P, a_, b_, rt, "down"), flank(P, a_, b_, rt, "up")
        mf = np.zeros(len(x), dtype=bool)
        mf[fl_d + fl_u] = True
        im = a_ + int(np.argmax(yp[a_:b_ + 1]))
        geo = {"pipe (x ≤ x_contraction)": int((m & (x <= x_con)).sum()) if x_con is not None else None,
               "contraction (x_contraction < x < 0)": int((m & (x > x_con) & (x < 0)).sum()) if x_con is not None else None,
               "x ≥ 0": int((m & (x >= 0)).sum())}
        ref_pts = {}
        for side, fl in (("down", fl_d), ("up", fl_u)):
            if fl:
                ref_pts[f"max_vs_first_{side}stream_flank_node"] = {"flank_node": node_rec(P, fl[0], rt, D), **point_ratio(P, im, fl[0])}
        regions.append({
            "i_range": [int(a_), int(b_)], "x_rt": [fnum(x[a_]), fnum(x[b_])], "x_m": [fnum(P["x"][a_]), fnum(P["x"][b_])],
            "n_nodes": int(b_ - a_ + 1),
            "area_m2": fnum(w[m].sum()), "area_frac_pct": fnum(100 * w[m].sum() / w[good].sum()), "geometry_nodes": geo,
            "max": node_rec(P, im, rt, D),
            "flank": {"down": {"n": len(fl_d), "x_rt": [fnum(x[fl_d[0]]), fnum(x[fl_d[-1]])] if fl_d else None},
                      "up": {"n": len(fl_u), "x_rt": [fnum(x[fl_u[-1]]), fnum(x[fl_u[0]])] if fl_u else None,
                             "note": "入口の端に接する (上流の y₁⁺ ≤ 1 の節点が無い)" if not fl_u and a_ == 0 else None}},
            "region_vs_flank": ratio_block(P, m, mf),
            "region_vs_flank_design_d1": {"region_gmean": gmean(D["d1_m"], w, m), "flank_gmean": gmean(D["d1_m"], w, mf),
                                          "ratio": (gmean(D["d1_m"], w, m) / gmean(D["d1_m"], w, mf)) if mf.any() else None},
            **ref_pts,
            "y1_over_design": {"min": fnum(np.min(P["y1_m"][m] / D["d1_m"][m])), "max": fnum(np.max(P["y1_m"][m] / D["d1_m"][m]))},
            "nodes": [node_rec(P, i, rt, D) for i in sorted(set(range(a_, b_ + 1)) | set(fl_d) | set(fl_u))],
        })
    R["wall"] = {
        "tool": t_wall, "n_nodes": int(len(x)), "n_evaluated": int(good.sum()), "area_total_m2": fnum(w.sum()),
        "y1plus": {"area_mean": fnum(np.sum(w[good] * yp[good]) / np.sum(w[good])), "node_mean": fnum(np.mean(yp[good])),
                   "max": node_rec(P, imax, rt, D), "n_over": int(over.sum()), "area_over_m2": fnum(w[over].sum()),
                   "area_over_pct": fnum(100 * w[over].sum() / w[good].sum())},
        "mu_w": {"min": fnum(np.min(P["mu_w"])), "max": fnum(np.max(P["mu_w"])),
                 "source": "check_wall_resolution.py: 体積出力 res_<step>.h5 の VALUE/vis_lam があればそれ、無ければ Sutherland (viscMethod 1)"},
        "T_w": {"min": fnum(np.min(P["T_w"])), "max": fnum(np.max(P["T_w"]))},
        "geometry": {"x_in_rt": pu.get("x_in"), "x_contraction_rt": x_con, "r_pipe_rt": pu.get("r_U"), "source": "prepare_info.json physical_upstream.wall"},
        "over_regions": regions,
        "bins": bins_table(P, rt, D),
        "note": ("y₁ は check_wall_resolution.py の第一内部節点までの壁法線の距離、τ_t は接線の壁応力、μ_w は vis_lam。区間と近くの場所の比は"
                 "面積の重みつきの幾何平均。log10_terms は log10(y₁⁺ の比) の厳密な分け方 (和が log10_y1plus_ratio)。"),
    }
    # --- 3. 設計の第一層の記録 ---
    rr = P["y1_m"] / D["d1_m"]
    R["design_first_layer"] = {
        "mesh_block_source": "prepare_info.json mesh.block", "problem_yaml": a.problem,
        "problem_yaml_sha256": sha256_file(case / a.problem) if (case / a.problem).is_file() else None,
        "same_as_problem_yaml_mesh": same_block, "branch": D["branch"], "wall_normal_layer": block.get("wall_normal_layer"),
        "normal_distance_is_d1": D["normal_layer"], "nj": block.get("nj"), "axis_cap_frac": block.get("axis_cap_frac"),
        "y1_over_design_all": {"min": fnum(np.min(rr)), "max": fnum(np.max(rr)), "median": fnum(np.median(rr)),
                               "i_max_dev": int(np.argmax(np.abs(rr - 1))), "x_rt_max_dev": fnum(x[int(np.argmax(np.abs(rr - 1)))])},
        "note": ("d₁ = r_w (1 − s_{nj−2}) [m]: 格子の生成と同じ関数で各 station の第一セル厚 (局所半径比の表 × r_w を axis_cap_frac の等比で "
                 "正規化した値)。wall_normal_layer があれば第一内部節点は壁法線上の距離 d₁ に置かれる (d₁/r_w ≪ d_n)。"),
    }
    # --- 4. 補助: 形 C・370 K の run ---
    if aux_P is not None:
        ainfo = jload(case / a.aux_run / "prepare_info.json") or {}
        ablock = (ainfo.get("mesh") or {}).get("block") or {}
        Da = design_first_layer(ablock, float(ainfo.get("scale_m", rt)), aux_P["x"], aux_P["y"])
        A = {k: np.interp(P["x"], aux_P["x"], aux_P[k]) for k in ("y", "y1_m", "tau_t", "rho_w", "mu_w", "y1plus")}
        A["abs_table_m"] = np.interp(P["x"], aux_P["x"], Da["abs_table_m"])
        A["d1_m"] = np.interp(P["x"], aux_P["x"], Da["d1_m"])
        cv = (case / a.aux_run / "CONVERGENCE_VERDICT.txt")
        cvl = [ln.strip() for ln in cv.read_text().splitlines() if ln.startswith(("OVERALL", "CONVERGENCE_EXIT"))] if cv.is_file() else None
        reg = []
        for r_ in regions:
            i0, i1 = r_["i_range"]
            m = np.zeros(len(x), dtype=bool)
            m[i0:i1 + 1] = True
            im = int(r_["max"]["i"])
            reg.append({"x_rt": r_["x_rt"],
                        "H_over_C_gmean": {k: (gmean(P[k], w, m) / gmean(A[k], w, m)) for k in ("y1_m", "tau_t", "rho_w", "mu_w", "y1plus")},
                        "H_over_C_design_d1_gmean": gmean(D["d1_m"], w, m) / gmean(A["d1_m"], w, m),
                        "H_over_C_abs_table_minmax": [fnum(np.min(D["abs_table_m"][m] / A["abs_table_m"][m])),
                                                      fnum(np.max(D["abs_table_m"][m] / A["abs_table_m"][m]))],
                        "C_y1plus_max_in_range": fnum(np.max(A["y1plus"][m])),
                        "at_H_max": {"x_rt": r_["max"]["x_rt"], "C": {k: fnum(A[k][im]) for k in ("y", "y1_m", "tau_t", "rho_w", "mu_w", "y1plus")},
                                     "H_over_C": {k: fnum(P[k][im] / A[k][im]) for k in ("y", "y1_m", "tau_t", "rho_w", "mu_w", "y1plus")}}})
        R["aux_C"] = {"run": a.aux_run, "step": a.aux_step, "tool": aux_tool, "convergence_verdict_file": cvl,
                      "run_end_m9_watch": (jload(case / a.aux_run / "RUN_END.json") or {}).get("m9_watch"),
                      "mesh": {k: ablock.get(k) for k in ("ni", "nj", "axis_cap_frac", "x_respace", "wall_normal_layer")},
                      "physical_upstream": ainfo.get("physical_upstream"),
                      "sha256": {"wall_dump": sha256_file(case / a.aux_run / f"res_wall_3_{a.aux_step}.h5")},
                      "regions": reg, "bins_C": bins_table(aux_P, rt),
                      "note": ("補助の記録だけ (主の結論に使わない): 形 C・370 K・g0c の格子 (nj 121、x_respace) の run の最後の壁の出力を H の壁の節点の x に"
                               "線形補間した値。収束の状態・格子・段の構成は run_0951 と違う。abs_table = 第一層の表の比 × 局所半径 (格子の生成の前の値)")}
    R["facts_summary"] = summary_lines(R)
    jwrite(out_p, R)
    print("\n".join(R["facts_summary"]))
    print(f"→ {out_p}")
    return 0 if not err else 2


def summary_lines(R: dict) -> list:
    out = []
    bc = R.get("boundary_classification") or {}
    out += bc.get("facts", [])
    for k in ("tool_all_groups", "tool_wall_only"):
        t = bc.get(k) or {}
        out.append(f"{k}: rc {t.get('rc')} {t.get('verdict')}")
    wl = R.get("wall") or {}
    yp = wl.get("y1plus") or {}
    out.append(f"壁: 面積平均 y₁⁺ {yp.get('area_mean')}, 1 超えの面積 {yp.get('area_over_pct')} %, 最大 {(yp.get('max') or {}).get('y1plus')} "
               f"(x {(yp.get('max') or {}).get('x_rt')} r_t)")
    for r_ in wl.get("over_regions", []):
        rv = r_["region_vs_flank"]
        out.append(f"y₁⁺ > 1 の区間 x {r_['x_rt']} r_t: 節点 {r_['n_nodes']}, 面積 {r_['area_frac_pct']} %, 形 {r_['geometry_nodes']}, "
                   f"区間/近くの場所 {rv.get('ratio_region_over_flank')}, log10 の項 {rv.get('log10_terms')}")
    df = R.get("design_first_layer") or {}
    out.append(f"y₁/設計 d₁ (全域): {df.get('y1_over_design_all')}")
    for r_ in (R.get("aux_C") or {}).get("regions", []):
        out.append(f"補助 H/C (x {r_['x_rt']}): {r_['H_over_C_gmean']}, 設計 d₁ の比 {r_['H_over_C_design_d1_gmean']}, C の y₁⁺ 最大 {r_['C_y1plus_max_in_range']}")
    return out


def jwrite(p: Path, d) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(d, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    tmp.replace(p)


if __name__ == "__main__":
    sys.exit(main())
