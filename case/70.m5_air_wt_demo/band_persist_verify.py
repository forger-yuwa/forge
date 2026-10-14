"""case/70: edge の縁 y* の推定の持続 (`edge_band_positions(..., persist)`、既定 4) の検証 — plan tooling-nozzle-cooled-wall-ns-settling
§6 S2b (a)(b)(d) (§5.1 #3、forge 0 step、NS・Euler の場は読むだけ)。(c) は単体の試験 design/tests/run_band_persist_tests.py。

(a) 窓をまたぐ再現性: NS run_0010b_ns_c2b_cont60k (ファイル名の step は継続の中の step) の非重複の窓
    W1 = res_10000・15000・20000・25000・30000 (通算 70k〜90k)、W2 = res_40000・45000・50000・55000・60000 (通算 100k〜120k)、
    Euler run_0012b_euler_c2b_cont12k (最新の res)。各窓で帯位置を独立に選ぶ抽出を 5 枚ずつ行い δ_E を窓で平均 (独立)。
    W2 の各枚の帯位置 y_b(x) を窓の中で (算術) 平均したものに固定して W1・W2 の各枚を抽出し窓で平均 (共通帯)。
    帯の再選択による差 = (独立の W1 平均 − 独立の W2 平均) − (共通帯の W1 平均 − 共通帯の W2 平均)、使用域 [0, x_e]
    (x_e = run_0010b の wall_repr.json の domain[1]) の最大 ≤ 1e-4 r_t で合格。persist 1 でも同じ量を出す (記録)。
(b) 非劣化: case/45 の 4 対 (NS × Euler、最新の res、/home/sano/work/forge の run を読むだけ) で persist 4 と 1 の δ_E の差。
    x ∈ [0, x_e − 3] の最大 ≤ 1e-4 r_t で合格 (x_e = 抽出の表の終点)、出口の端 (x > x_e − 3) は値と位置を記録。y* の NaN の割合も記録。
(d) NaN の区間の埋め方の A/B (persist 4): A = 今の P-spline (= (a)(b) の persist 4 の独立の抽出)、B = y* の連続した NaN の区間
    (使用域 [0, x_e] と hold の基準位置にかかるもの) を、その両隣の有効点での A の kappa_s (hold の前の P-spline の値) を結ぶ直線に
    置き換えてから同じ hold を当てる (有効点の値は変えない)。W1・W2 の 10 枚と case/45 の 4 対で、使用域 [0, x_e] の
    max|δ_E,B − δ_E,A| ≤ 1e-4 r_t で合格。最大の位置と NaN の区間の長さを記録。対象の区間で片側の支持点が無ければ「確認不能」(合格にしない)。

帯位置の固定 (共通帯) と B の埋め方は `edge_band_positions` の戻り値を差し替える測定用の方法 (実装ではない)。抽出 1 本 約 45 s。

usage (design/.venv-opt):
  python band_persist_verify.py [--jobs 4]      # 抽出 (済んだものは飛ばす) → 要約。出力 _band_persist_verify/verify.json と標準出力
  python band_persist_verify.py summary         # 要約だけ
  python band_persist_verify.py one JOB_JSON    # 抽出 1 本 (内部用)
"""
from __future__ import annotations

import inspect
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

OUT = C / "_band_persist_verify"
NS_B = C / "run_0010b_ns_c2b_cont60k"
EU_B = C / "run_0012b_euler_c2b_cont12k"
W1 = (10000, 15000, 20000, 25000, 30000)            # 通算 70k〜90k
W2 = (40000, 45000, 50000, 55000, 60000)            # 通算 100k〜120k
C45 = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
PAIRS45 = (("run_0038_ns_final_rt77p02", "run_0037_euler_rt77p02"),
           ("run_0036_ns_final_rt77p05", "run_0035_euler_rt77p05"),
           ("run_0034_ns_ib_pass0_Lpipe10", "run_0033_euler_Lpipe10"),
           ("run_0022_ns_ib_pass0", "run_0001_euler_shortest_dry"))
TOL = 1e-4                                          # [r_t]
EXIT_BAND = 3.0                                     # (b) 出口の端 x > x_e − 3 (hold の区間)


def x_e_case70() -> float:
    return float(json.loads((NS_B / "wall_repr.json").read_text())["domain"][1])


# ================================================================================================ 抽出 1 本
def patch_step(ns_run: Path, step: int):
    """NS run の res_<step>.h5 を読む (band_persist_ab.py と同じ差し替え)。"""
    import h5py
    orig = DS._load_structured
    ns_run = ns_run.resolve()

    def load(run_dir):
        run_dir = Path(run_dir)
        if run_dir.resolve() != ns_run:
            return orig(run_dir)
        info = json.loads((run_dir / "prepare_info.json").read_text())
        ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
        S = float(info["scale_m"])
        with h5py.File(run_dir / "nozzle.h5") as f:
            nc = f["/MESH/COORD"][:].reshape(-1, 3)
        res = run_dir / f"res_{step}.h5"
        with h5py.File(res) as f:
            q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
        return dict(info=info, S=S, res=res.name, x=nc[:, 0].reshape(ni, nj) / S, r=nc[:, 1].reshape(ni, nj) / S, q=q.reshape(ni, nj))

    DS._load_structured = load


def patch_fix_yb(csv: Path, rec: dict):
    """帯位置 y_b(x) を表 (x_rt, y_b) の値に固定する (梯子・縁の推定は計算して記録するが、y_b だけ差し替える)。"""
    tb = np.loadtxt(csv, delimiter=",", skiprows=1)
    xf, ybf = tb[:, 0], tb[:, 1]
    orig = DS.edge_band_positions

    def fixed(x, *a, **k):
        out = dict(orig(x, *a, **k))
        xx = np.asarray(x, dtype=float)
        if not (xx.min() >= xf.min() - 1e-6 and xx.max() <= xf.max() + 1e-6):
            raise RuntimeError("固定する帯位置の表が抽出の範囲を覆わない")
        rec["same_x_as_table"] = bool(len(xx) == len(xf) and np.array_equal(xx, xf))
        out["y_b_reselected"] = out["y_b"]
        out["y_b"] = np.interp(xx, xf, ybf)
        return out

    DS.edge_band_positions = fixed


def nan_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """True の連続区間 (s, e) (両端を含む添字) の一覧。"""
    runs, i, n = [], 0, len(mask)
    while i < n:
        if mask[i]:
            s = i
            while i + 1 < n and mask[i + 1]:
                i += 1
            runs.append((s, i))
        i += 1
    return runs


def apply_hold(x, ks, hold_lo, hold_hi_from_end):
    """edge_band_positions の端の値保持と同じ式。"""
    if hold_lo is not None and np.any(x >= hold_lo):
        ks = np.where(x < hold_lo, float(np.interp(hold_lo, x, ks)), ks)
    if hold_hi_from_end is not None:
        xh = float(x.max()) - hold_hi_from_end
        ks = np.where(x > xh, float(np.interp(xh, x, ks)), ks)
    return ks


def patch_fill_b(use_hi: float, rec: dict):
    """(d) の B: y* の連続した NaN の区間のうち使用域 [0, use_hi] か hold の基準位置にかかるものを、両隣の有効点での
    A の kappa_s (hold の前の P-spline の値) を結ぶ直線に置き換えてから同じ hold を当てる。有効点の値は A のまま。"""
    orig = DS.edge_band_positions
    dflt = {k: v.default for k, v in inspect.signature(orig).parameters.items() if v.default is not inspect.Parameter.empty}

    def fill_b(x, y_ladders, d_ladders, d_in, r_w, **k):
        hold_lo = k.get("hold_lo", dflt["hold_lo"]); hold_hi = k.get("hold_hi_from_end", dflt["hold_hi_from_end"])
        c = k.get("c", dflt["c"])
        A = orig(x, y_ladders, d_ladders, d_in, r_w, **k)
        pre = orig(x, y_ladders, d_ladders, d_in, r_w, **dict(k, hold_lo=None, hold_hi_from_end=None))
        xx = np.asarray(x, dtype=float)
        ks_pre = np.asarray(pre["kappa_s"], dtype=float)
        rec["A_reproduced_from_prehold"] = bool(np.array_equal(apply_hold(xx, ks_pre, hold_lo, hold_hi), A["kappa_s"]))
        rec["y_star_same"] = bool(np.array_equal(pre["y_star"], A["y_star"], equal_nan=True))
        refs = []                                       # hold の基準位置 (その点の kappa_s が端の区間に写される)
        if hold_lo is not None and np.any(xx >= hold_lo):
            refs.append(float(hold_lo))
        if hold_hi is not None:
            refs.append(float(xx.max()) - hold_hi)
        ks_b = ks_pre.copy()
        ivs = []
        n = len(xx)
        for s, e in nan_runs(np.isnan(A["y_star"])):
            lo_x = xx[s - 1] if s > 0 else xx[s]
            hi_x = xx[e + 1] if e + 1 < n else xx[e]
            in_use = bool(xx[e] >= 0.0 and xx[s] <= use_hi)
            at_ref = [float(r) for r in refs if lo_x < r < hi_x or xx[s] <= r <= xx[e]]
            if not (in_use or at_ref):
                continue
            two_sided = s > 0 and e + 1 < n
            iv = dict(i0=int(s), i1=int(e), n=int(e - s + 1), x0=float(xx[s]), x1=float(xx[e]), dx=float(hi_x - lo_x),
                      in_use=in_use, at_hold_ref=at_ref, two_sided=bool(two_sided))
            if two_sided:
                xl, xr = xx[s - 1], xx[e + 1]
                ks_b[s:e + 1] = ks_pre[s - 1] + (ks_pre[e + 1] - ks_pre[s - 1]) * (xx[s:e + 1] - xl) / (xr - xl)
                iv["max_abs_dkappa_prehold"] = float(np.max(np.abs(ks_b[s:e + 1] - ks_pre[s:e + 1])))
            ivs.append(iv)
        ks_b = apply_hold(xx, ks_b, hold_lo, hold_hi)
        rec["intervals"] = ivs
        rec["hold_refs"] = refs
        rec["undeterminable"] = [iv for iv in ivs if not iv["two_sided"]]
        out = dict(A)
        out["kappa_s_A"] = A["kappa_s"]
        out["kappa_s"] = ks_b
        out["y_b"] = c * np.asarray(A["y_ref"]) * np.exp(ks_b)
        return out

    DS.edge_band_positions = fill_b


def one(job: dict):
    od = Path(job["out"]); od.mkdir(parents=True, exist_ok=True)
    ns, eu = Path(job["ns"]), Path(job["eu"])
    if job.get("step"):
        patch_step(ns, int(job["step"]))
    rec: dict = {}
    if job.get("fix_yb"):
        patch_fix_yb(Path(job["fix_yb"]), rec)
    if job.get("fill_b"):
        patch_fill_b(float(job["use_hi"]), rec)
    summ = DL.extract_and_merge(ns, eu, band_select="edge", edge_persist=int(job["persist"]), out_dir=od)
    (od / "done.json").write_text(json.dumps(dict(job=job, rec=rec, res_ns=summ["massflow"]["res_ns"],
                                                  res_euler=summ["massflow"]["res_euler"]), indent=1, ensure_ascii=False))


# ================================================================================================ 段取り
def jd(**k):
    return {kk: (str(v) if isinstance(v, Path) else v) for kk, v in k.items()}


def jobs_phase1(x_e70: float) -> list[dict]:
    J = []
    for K in (4, 1):
        for st in W1 + W2:
            J.append(jd(ns=NS_B, eu=EU_B, step=st, persist=K, out=OUT / "r0010b" / f"K{K}" / "indep" / str(st)))
    for ns, eu in PAIRS45:
        for K in (4, 1):
            J.append(jd(ns=C45 / ns, eu=C45 / eu, step=0, persist=K, out=OUT / "c45" / ns / f"K{K}"))
    for st in W1 + W2:
        J.append(jd(ns=NS_B, eu=EU_B, step=st, persist=4, fill_b=True, use_hi=x_e70, out=OUT / "r0010b" / "K4" / "fillB" / str(st)))
    return J


def jobs_fillb45() -> list[dict]:
    """case/45 の B: 使用域の上端 = persist 4 の抽出の表の終点 (A の抽出の後に決まる)。"""
    J = []
    for ns, eu in PAIRS45:
        a = DL.read_delta_r_next(OUT / "c45" / ns / "K4" / "delta_r_next.csv")
        J.append(jd(ns=C45 / ns, eu=C45 / eu, step=0, persist=4, fill_b=True, use_hi=float(np.asarray(a["x_rt"])[-1]),
                    out=OUT / "c45" / ns / "K4_fillB"))
    return J


def common_band(K: int) -> Path:
    """W2 の独立の抽出の帯位置 y_b(x) を窓の中で算術平均した表。"""
    tabs = [np.loadtxt(OUT / "r0010b" / f"K{K}" / "indep" / str(st) / "delta_r_edge.csv", delimiter=",", skiprows=1) for st in W2]
    x0 = tabs[0][:, 0]
    assert all(np.array_equal(t[:, 0], x0) for t in tabs), "W2 の各枚の帯位置の x が揃わない"
    yb = np.mean([t[:, 5] for t in tabs], axis=0)
    p = OUT / "r0010b" / f"K{K}" / "common_yb_W2mean.csv"
    np.savetxt(p, np.c_[x0, yb], delimiter=",", comments="", header="x_rt,y_b_W2mean")
    return p


def jobs_phase2() -> list[dict]:
    J = []
    for K in (4, 1):
        p = common_band(K)
        for st in W1 + W2:
            J.append(jd(ns=NS_B, eu=EU_B, step=st, persist=K, fix_yb=p, out=OUT / "r0010b" / f"K{K}" / "common" / str(st)))
    return J


def run_jobs(J: list[dict], n: int):
    todo = [j for j in J if not (Path(j["out"]) / "done.json").exists()]
    print(f"抽出 {len(todo)} 本 (済み {len(J) - len(todo)} 本) を {n} 本ずつ", flush=True)

    def go(j):
        Path(j["out"]).mkdir(parents=True, exist_ok=True)
        with open(Path(j["out"]) / "log.txt", "w") as lf:
            r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "one", json.dumps(j)], stdout=lf, stderr=subprocess.STDOUT)
        print(f"  {'ok ' if r.returncode == 0 else 'ERR'} {Path(j['out']).relative_to(OUT)}", flush=True)
        return r.returncode

    with ThreadPoolExecutor(max_workers=n) as ex:
        rcs = list(ex.map(go, todo))
    if any(rcs):
        raise SystemExit(f"抽出の失敗 {sum(1 for r in rcs if r)} 本 (各 out の log.txt)")


# ================================================================================================ 要約
def dE(d: Path):
    t = DL.read_delta_r_next(d / "delta_r_next.csv")
    return np.asarray(t["x_rt"], dtype=float), np.asarray(t["delta_E"], dtype=float)


def edge_tab(d: Path):
    return np.loadtxt(d / "delta_r_edge.csv", delimiter=",", skiprows=1)


def argmax_abs(x, v, m):
    i = int(np.argmax(np.where(m, np.abs(v), -np.inf)))
    return float(np.abs(v[i])), float(x[i])


def nan_stats(tab, lo=None, hi=None):
    x, ys = tab[:, 0], tab[:, 2]
    m = np.ones(len(x), bool) if lo is None else (x >= lo) & (x <= hi)
    return dict(n=int(m.sum()), n_nan=int(np.isnan(ys[m]).sum()), frac=float(np.isnan(ys[m]).mean()))


def summary_a(x_e: float) -> dict:
    res = {"x_e": x_e, "use": [0.0, x_e], "W1_files": [f"res_{s}" for s in W1], "W2_files": [f"res_{s}" for s in W2]}
    for K in (4, 1):
        base = OUT / "r0010b" / f"K{K}"
        mean = {}
        x0 = None
        for arm in ("indep", "common"):
            for wn, W in (("W1", W1), ("W2", W2)):
                arr = []
                for st in W:
                    x, d = dE(base / arm / str(st))
                    if x0 is None:
                        x0 = x
                    assert np.array_equal(x, x0), "抽出の x が揃わない"
                    arr.append(d)
                arr = np.array(arr)
                mean[(arm, wn)] = arr.mean(axis=0)
                mean[(arm, wn, "spread")] = arr.max(axis=0) - arr.min(axis=0)
        use = (x0 >= 0.0) & (x0 <= x_e)
        D_ind = mean[("indep", "W1")] - mean[("indep", "W2")]
        D_com = mean[("common", "W1")] - mean[("common", "W2")]
        resel = D_ind - D_com
        r_max, r_x = argmax_abs(x0, resel, use)
        i_max, i_x = argmax_abs(x0, D_ind, use)
        c_max, c_x = argmax_abs(x0, D_com, use)
        nan_w = {}
        for st in W1 + W2:
            nan_w[str(st)] = nan_stats(edge_tab(base / "indep" / str(st)), 0.0, x_e)["n_nan"]
        fixed_x = all(json.loads((base / "common" / str(st) / "done.json").read_text())["rec"].get("same_x_as_table") for st in W1 + W2)
        res[f"K{K}"] = dict(
            reselection_max_rt=r_max, reselection_x=r_x,
            indep_W1_minus_W2_max_rt=i_max, indep_x=i_x,
            common_W1_minus_W2_max_rt=c_max, common_x=c_x,
            indep_within_window_spread_max_rt={w: float(np.max(mean[("indep", w, "spread")][use])) for w in ("W1", "W2")},
            common_within_window_spread_max_rt={w: float(np.max(mean[("common", w, "spread")][use])) for w in ("W1", "W2")},
            y_star_nan_in_use_per_file=nan_w, n_stations_use=int(use.sum()),
            common_band_same_x=bool(fixed_x),
            pass_=bool(r_max <= TOL))
    # 記録: 同じ枚の persist 4 と 1 の独立の抽出の δ_E の差 (使用域; 判定には使わない)
    k41 = {}
    for st in W1 + W2:
        x4, e4 = dE(OUT / "r0010b" / "K4" / "indep" / str(st)); x1, e1 = dE(OUT / "r0010b" / "K1" / "indep" / str(st))
        assert np.array_equal(x4, x1)
        mx, xx_ = argmax_abs(x4, e4 - e1, (x4 >= 0.0) & (x4 <= x_e))
        k41[str(st)] = dict(max_abs_rt=mx, x=xx_)
    res["K4_minus_K1_same_file_use"] = k41
    res["verdict"] = "PASS" if res["K4"]["pass_"] else "FAIL"
    return res


def summary_b() -> dict:
    res = {"pairs": {}}
    ok = True
    for ns, eu in PAIRS45:
        d4, d1 = OUT / "c45" / ns / "K4", OUT / "c45" / ns / "K1"
        x4, e4 = dE(d4); x1, e1 = dE(d1)
        assert np.array_equal(x4, x1)
        x_e = float(x4[-1])
        diff = e4 - e1
        main = (x4 >= 0.0) & (x4 <= x_e - EXIT_BAND)
        ex = x4 > x_e - EXIT_BAND
        m_max, m_x = argmax_abs(x4, diff, main)
        e_max, e_x = argmax_abs(x4, diff, ex)
        e_val = float(diff[int(np.argmin(np.abs(x4 - e_x)))])
        dn = json.loads((d4 / "done.json").read_text())
        rec = dict(ns=ns, euler=eu, res_ns=dn["res_ns"], res_euler=dn["res_euler"], x_e=x_e,
                   main_max_rt=m_max, main_x=m_x, exit_max_rt=e_max, exit_x=e_x, exit_signed_rt=e_val,
                   exit_rel_to_dE=float(e_val / e1[int(np.argmin(np.abs(x4 - e_x)))]),
                   y_star_nan_K4=nan_stats(edge_tab(d4)), y_star_nan_K1=nan_stats(edge_tab(d1)),
                   y_star_nan_K4_use=nan_stats(edge_tab(d4), 0.0, x_e), y_star_nan_K1_use=nan_stats(edge_tab(d1), 0.0, x_e),
                   pass_=bool(m_max <= TOL))
        ok &= rec["pass_"]
        res["pairs"][ns] = rec
    res["verdict"] = "PASS" if ok else "FAIL"
    return res


def summary_d(x_e70: float) -> dict:
    items = []
    for st in W1 + W2:
        items.append((f"r0010b/res_{st}", OUT / "r0010b" / "K4" / "indep" / str(st), OUT / "r0010b" / "K4" / "fillB" / str(st), x_e70))
    for ns, _ in PAIRS45:
        a = OUT / "c45" / ns / "K4"
        items.append((f"c45/{ns}", a, OUT / "c45" / ns / "K4_fillB", None))
    out, ok, undet = {}, True, False
    for tag, da, db, x_e in items:
        xa, ea = dE(da); xb, eb = dE(db)
        assert np.array_equal(xa, xb)
        if x_e is None:
            x_e = float(xa[-1])
        use = (xa >= 0.0) & (xa <= x_e)
        diff = eb - ea
        mx, mxx = argmax_abs(xa, diff, use)
        m3, m3x = argmax_abs(xa, diff, (xa >= 0.0) & (xa <= x_e - EXIT_BAND))
        rec = json.loads((db / "done.json").read_text())["rec"]
        ivs = rec["intervals"]
        near = min(ivs, key=lambda iv: 0.0 if iv["x0"] <= mxx <= iv["x1"] else min(abs(mxx - iv["x0"]), abs(mxx - iv["x1"]))) if ivs else None
        ya, yb = edge_tab(da)[:, 2], edge_tab(db)[:, 2]
        r = dict(x_e=x_e, max_abs_rt=mx, x_at_max=mxx, max_abs_excl_exit_end_rt=m3, x_at_max_excl_exit_end=m3x,
                 n_intervals=len(ivs), interval_n_stations=[iv["n"] for iv in ivs], interval_dx_rt=[round(iv["dx"], 4) for iv in ivs],
                 interval_nearest_max=near, longest_interval=(max(ivs, key=lambda iv: iv["dx"]) if ivs else None),
                 undeterminable=rec["undeterminable"], A_reproduced_from_prehold=rec["A_reproduced_from_prehold"],
                 y_star_same_as_A=bool(np.array_equal(ya, yb, equal_nan=True)), hold_refs=rec["hold_refs"],
                 pass_=bool(mx <= TOL and not rec["undeterminable"]))
        ok &= r["pass_"]
        undet |= bool(rec["undeterminable"])
        out[tag] = r
    return dict(items=out, verdict=("PASS" if ok else ("UNDETERMINABLE" if undet else "FAIL")))


def summary(x_e70: float):
    a = summary_a(x_e70); b = summary_b(); d = summary_d(x_e70)
    res = dict(plan="tooling-nozzle-cooled-wall-ns-settling §6 S2b", tol_rt=TOL, a=a, b=b, d=d)
    (OUT / "verify.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"\n=== (a) 窓をまたぐ再現性 (run_0010b W1 = res_10000〜30000 / W2 = res_40000〜60000, 使用域 [0, {x_e70:.4f}]) ===")
    for K in (4, 1):
        k = a[f"K{K}"]
        print(f"  persist {K}: 帯の再選択による差 {k['reselection_max_rt']:.3e} r_t (x {k['reselection_x']:.3f})  "
              f"独立 W1−W2 {k['indep_W1_minus_W2_max_rt']:.3e} (x {k['indep_x']:.3f})  共通帯 W1−W2 {k['common_W1_minus_W2_max_rt']:.3e} "
              f"(x {k['common_x']:.3f})  窓の中の幅 (独立) {k['indep_within_window_spread_max_rt']}  {'PASS' if k['pass_'] else 'FAIL'}")
    print("  記録: 同じ枚の persist 4 − 1 の δ_E (使用域の最大): " + ", ".join(
        f"res_{st} {v['max_abs_rt']:.2e} (x {v['x']:.2f})" for st, v in a["K4_minus_K1_same_file_use"].items()))
    print(f"  VERDICT (a, persist 4): {a['verdict']}")
    print("\n=== (b) 非劣化 (case/45, persist 4 − persist 1) ===")
    for ns, r in b["pairs"].items():
        print(f"  {ns} ({r['res_ns']}): [0, x_e−3] 最大 {r['main_max_rt']:.3e} r_t (x {r['main_x']:.3f})  出口の端 {r['exit_signed_rt']:+.3e} "
              f"(x {r['exit_x']:.3f}, δ_E 比 {r['exit_rel_to_dE']:+.4%})  y* NaN K4 {r['y_star_nan_K4']['frac']:.3f} / K1 "
              f"{r['y_star_nan_K1']['frac']:.3f}  {'PASS' if r['pass_'] else 'FAIL'}")
    print(f"  VERDICT (b): {b['verdict']}")
    print("\n=== (d) NaN の区間の埋め方の A/B (persist 4, 使用域 [0, x_e]) ===")
    for tag, r in d["items"].items():
        nr = r["interval_nearest_max"]
        near = "-" if nr is None else "x {:.3f}〜{:.3f} ({} 断面)".format(nr["x0"], nr["x1"], nr["n"])
        longest = float("nan") if r["longest_interval"] is None else r["longest_interval"]["dx"]
        print(f"  {tag}: max|B−A| {r['max_abs_rt']:.3e} r_t (x {r['x_at_max']:.3f})  区間 {r['n_intervals']} (最長 {longest:.3f} r_t)  "
              f"最大に近い区間 {near}  確認不能 {len(r['undeterminable'])}  {'PASS' if r['pass_'] else 'FAIL'}")
    print(f"  VERDICT (d): {d['verdict']}")
    print(f"\n{OUT / 'verify.json'}")


def main(argv):
    if argv and argv[0] == "one":
        one(json.loads(argv[1]))
        return
    x_e70 = x_e_case70()
    if argv and argv[0] == "summary":
        summary(x_e70)
        return
    n = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 4
    OUT.mkdir(exist_ok=True)
    run_jobs(jobs_phase1(x_e70), n)
    run_jobs(jobs_fillb45() + jobs_phase2(), n)
    summary(x_e70)


if __name__ == "__main__":
    main(sys.argv[1:])
