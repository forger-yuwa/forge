"""plan tooling-nozzle-wall-single-bspline §6.0 W0〜W5 の判定 (CFD 0 step、forge は起動しない)。判定条件は §6.0 の事前登録どおり。

前提: wsb_prepare.py で次の成果物を作ってあること (_band_ab/wsb/ 以下):
  base_legacy (#2b 基準 = HEAD 7a505415 の design/ の写し、キー無し) / base_legacy_rep (同じ基準の再実行、自己再現の対照) /
  after_nokey (変更後、キー無し) / w3_A_legacy (変更後、physical_wall_repr: legacy) / w3_B_single_bspline (変更後、single_bspline)

usage: [CASE_RUNS=...] FORGE_CONVERTER=... python wsb_verify.py {w0|w3|w1w2|w4|w5|neg|all} → _band_ab/wsb/W*.json
  w0 は design/tests の変更後の結果 (_band_ab/wsb/tests_after.txt) も読む。
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
WSB = C / "_band_ab" / "wsb"
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
PN = "problem_d155_ns_finemesh_recal_final_mono.yaml"
SOLVER_FILES = ("solverConfig.yaml", "bcondConfig.yaml")
# §6.0 W3 の比較項目 → nozzle.h5 のデータセット (前方一致)。どれにも当たらないデータセットも「その他」として全部比べる
W3_ITEMS = {
    "座標": ("MESH/COORD",),
    "接続": ("MESH/CONNE", "CELLS/STRUCT", "PLANES/STRUCT", "VIZMESH/CONNE", "CELLS/regionId"),
    "境界対応": ("BCONDS/*/iBPlanes", "BCONDS/*/iCells", "BCONDS/*/iPlanes", "BCONDS/*/vizBfaceNodes", "BCONDS/*/vizBfaceSizes"),
    "体積": ("CELLS/volume",),
    "面の幾何": ("PLANES/surfArea", "PLANES/surfVect", "PLANES/centCoords"),
    "centCoords": ("CELLS/centCoords",),
    "wall_dist": ("VALUE/wall_dist",),
    "全初期保存量": ("VALUE/ro", "VALUE/roUx", "VALUE/roUy", "VALUE/roUz", "VALUE/roe", "VALUE/roK", "VALUE/roOmega",
                 "VALUE/roY*", "BCONDS/*/VALUE/*"),
}


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _match(name: str, pat: str) -> bool:
    import fnmatch
    return fnmatch.fnmatchcase(name, pat)


def item_of(name: str) -> str:
    for item, pats in W3_ITEMS.items():
        if any(_match(name, pt) for pt in pats):
            return item
    return "その他"


def h5_compare(pa: Path, pb: Path) -> dict:
    """2 つの HDF5 をデータセット単位で比べる (ファイルのバイト列は同じ入力でも再現しないので使わない)。
    データセットごとに: 一致 / 不一致 (要素の不一致数・最大絶対差・float32 の ULP 差の最大)。属性も全部比べる。"""
    import h5py

    def walk(path):
        ds, at = {}, {}
        with h5py.File(path, "r") as f:
            def attrs(name, obj):
                if len(obj.attrs):
                    at[name or "/"] = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in obj.attrs.items()}
            attrs("/", f)

            def visit(name, obj):
                attrs(name, obj)
                if isinstance(obj, h5py.Dataset):
                    ds[name] = np.asarray(obj[()])
            f.visititems(visit)
        return ds, at

    da, aa = walk(pa)
    db, ab = walk(pb)
    rows = {}
    for name in sorted(set(da) | set(db)):
        if name not in da or name not in db:
            rows[name] = {"item": item_of(name), "identical": False, "why": f"片方にだけある (A: {name in da}, B: {name in db})"}
            continue
        a, b = da[name], db[name]
        if a.dtype != b.dtype or a.shape != b.shape:
            rows[name] = {"item": item_of(name), "identical": False, "why": f"dtype/shape {a.dtype}{a.shape} ≠ {b.dtype}{b.shape}"}
            continue
        if a.tobytes() == b.tobytes():
            rows[name] = {"item": item_of(name), "identical": True}
            continue
        ne = a != b
        r = {"item": item_of(name), "identical": False, "n_diff": int(ne.sum()), "n": int(a.size)}
        if a.dtype.kind == "f":
            af, bf = a.astype(float), b.astype(float)
            r["max_abs_diff"] = float(np.abs(af - bf).max())
            r["max_rel_diff"] = float((np.abs(af - bf) / np.maximum(np.abs(af), 1e-300)).max())
            if a.dtype == np.float32:
                ia, ib = a.view(np.int32).astype(np.int64), b.view(np.int32).astype(np.int64)
                r["max_ulp"] = int(np.abs(ia - ib).max())
            idx = np.where(ne.ravel())[0]
            r["first_diff_index"] = [int(i) for i in idx[:10]]
        rows[name] = r
    by_item = {}
    for name, r in rows.items():
        it = by_item.setdefault(r["item"], {"n_datasets": 0, "n_diff_datasets": 0, "diff": []})
        it["n_datasets"] += 1
        if not r["identical"]:
            it["n_diff_datasets"] += 1
            it["diff"].append(name)
    return {"datasets": rows, "by_item": by_item, "attrs_identical": aa == ab,
            "attrs_diff": (None if aa == ab else {k: [aa.get(k), ab.get(k)] for k in set(aa) | set(ab) if aa.get(k) != ab.get(k)}),
            "all_identical": bool(all(r["identical"] for r in rows.values()) and aa == ab)}


def dict_diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: {'(無)' if k not in a else '(有)'} ≠ {'(無)' if k not in b else '(有)'}")
            else:
                out += dict_diff(a[k], b[k], f"{path}/{k}")
    elif a != b:
        out.append(f"{path}: {str(a)[:80]} ≠ {str(b)[:80]}")
    return out


def compare_dirs(A: Path, B: Path) -> dict:
    """2 つの prepare_ns の成果物を比べる: 全ファイルの sha256、nozzle.h5 はデータセット単位、設定は辞書でも、
    prepare_info.json はキーの差、nozzle.msh は文字列の差 (中間ファイル。ソルバは読まない)。"""
    import yaml
    fa = {p.name for p in A.iterdir() if p.is_file()}
    fb = {p.name for p in B.iterdir() if p.is_file()}
    files = {}
    for n in sorted(fa | fb):
        if n in fa and n in fb:
            files[n] = {"identical_bytes": sha_file(A / n) == sha_file(B / n)}
        else:
            files[n] = {"only_in": "A" if n in fa else "B"}
    cfg = {}
    for n in SOLVER_FILES:
        ya, yb = yaml.safe_load((A / n).read_text()), yaml.safe_load((B / n).read_text())
        cfg[n] = {"identical_bytes": sha_file(A / n) == sha_file(B / n), "identical_parsed": ya == yb, "diff": dict_diff(ya, yb)}
    ia, ib = json.loads((A / "prepare_info.json").read_text()), json.loads((B / "prepare_info.json").read_text())
    msh = None
    if (A / "nozzle.msh").exists() and (B / "nozzle.msh").exists():
        la, lb = (A / "nozzle.msh").read_text().splitlines(), (B / "nozzle.msh").read_text().splitlines()
        nd = [i for i, (x, y) in enumerate(zip(la, lb)) if x != y]
        msh = {"identical": la == lb, "n_lines": [len(la), len(lb)], "n_diff_lines": len(nd) + abs(len(la) - len(lb)),
               "diff_lines": [{"line": i + 1, "A": la[i][:120], "B": lb[i][:120]} for i in nd[:10]]}
    return {"A": str(A), "B": str(B), "files": files, "solver_config": cfg, "h5": h5_compare(A / "nozzle.h5", B / "nozzle.h5"),
            "prepare_info_diff": dict_diff(ia, ib)[:60], "nozzle_msh": msh}


def _tests_summary(path: Path, logs: Path | None = None) -> dict:
    """run_all_tests.sh の summary (各テストの rc) と各ログの FAIL 行 (「FAIL 件数: 0」の集計行は数えない)。
    変更前の 7a505415 でも失敗する 3 本 (tests_before_7a505415.txt で確認) は除いて判定する。"""
    import re
    if not path.exists():
        return {"status": "未実行", "path": str(path)}
    excl = ("run_sern_gates_tests.py", "run_sern_moc_tests.py", "run_species_attrs_ic_tests.py")
    rows = {}
    for line in path.read_text().splitlines():
        if line.startswith("run_") and " rc=" in line:
            nm, rest = line.split(" ", 1)
            kv = dict(x.split("=") for x in rest.split())
            nf = None
            if logs is not None and (logs / nm.replace(".py", ".log")).exists():
                nf = sum(1 for l in (logs / nm.replace(".py", ".log")).read_text().splitlines() if re.match(r"^FAIL(?! 件数)", l))
            rows[nm] = {"rc": int(kv["rc"]), "fail_lines": nf}
    bad = {k: v for k, v in rows.items() if k not in excl and (v["rc"] != 0 or (v["fail_lines"] or 0) > 0)}
    return {"path": str(path), "logs": (str(logs) if logs else None), "n": len(rows), "excluded_pre_existing": list(excl),
            "excluded_rows": {k: rows.get(k) for k in excl}, "rows": rows, "fail_tests": bad,
            "pass": (not bad) and len(rows) > 0 and "run_wall_single_bspline_tests.py" in rows}


# ----------------------------------------------------------------------------------------------- W0
def judge_w0(cmp_main: dict, cmp_ctrl: dict | None = None) -> dict:
    """W0: 基準 (変更前のコード) と変更後 (キー無し) の成果物がビット同一。nozzle.h5 はデータセット単位 + 属性、その他のファイルはバイト列。
    ファイルの集合も同じであること (キー無しで壁ファイルが増えない)。"""
    files = cmp_main["files"]
    only = [n for n, v in files.items() if "only_in" in v]
    bytes_diff = [n for n, v in files.items() if n != "nozzle.h5" and not v.get("identical_bytes", False)]
    h5_ok = cmp_main["h5"]["all_identical"]
    ok = (not only) and (not bytes_diff) and h5_ok
    out = {"pass": bool(ok), "files_only_in_one": only, "files_bytes_diff_except_h5": bytes_diff,
           "nozzle_h5_datasets_and_attrs_identical": h5_ok,
           "nozzle_h5_file_bytes_identical": files.get("nozzle.h5", {}).get("identical_bytes")}
    if cmp_ctrl is not None:
        out["control_self_reproducibility"] = {
            "nozzle_h5_datasets_identical": cmp_ctrl["h5"]["all_identical"],
            "nozzle_h5_file_bytes_identical": cmp_ctrl["files"].get("nozzle.h5", {}).get("identical_bytes"),
            "other_files_bytes_diff": [n for n, v in cmp_ctrl["files"].items() if n != "nozzle.h5" and not v.get("identical_bytes", False)]}
    return out


def run_w0():
    base, rep, after = WSB / "base_legacy", WSB / "base_legacy_rep", WSB / "after_nokey"
    cm, cc = compare_dirs(base, after), compare_dirs(base, rep)
    j = judge_w0(cm, cc)
    tests = _tests_summary(WSB / "tests_after.txt", WSB / "tests_after_logs")
    hb = json.loads((WSB / "base_legacy.hashes.json").read_text())
    out = {"item": "W0 既定のビット同一 (plan §6.0)",
           "criterion": "キー無しで、#2b の基準成果物 (変更前のコード・同じ固定入力) と変更後のコードの成果物がビット単位で一致 "
                        "(nozzle.h5 はデータセット単位と属性; 変換器の HDF5 はバイト列が同じ入力でも再現しないため)。design/tests の既存テストが FAIL 0 "
                        "(変更前でも失敗する 3 本を除く)",
           "verdict": "PASS" if (j["pass"] and tests.get("pass")) else "FAIL",
           "artifacts": j, "tests_after": tests, "compare_base_vs_after": cm, "compare_base_vs_base_rep": cc,
           "base_hashes": str(WSB / "base_legacy.hashes.json"),
           "past_run_reproduction_run_0147": {"note": "別項目 (§6.0 W0)。ローカルで作り直した δ_r と物理壁は run_0147 (AWS で作成) の保存物と違う",
                                              "delta_r_initial": hb.get("delta_r_initial_vs_run_0147"),
                                              "wall_physical": hb.get("wall_physical_vs_run_0147"),
                                              "nozzle_h5": hb.get("run_0147_results")}}
    (WSB / "W0.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


# ----------------------------------------------------------------------------------------------- W3
def judge_w3(cmp: dict) -> dict:
    """W3: ソルバ入力 (nozzle.h5 の全データセットと属性 + solverConfig.yaml・bcondConfig.yaml) が全項目ビット同一か。"""
    h5 = cmp["h5"]
    cfg_ok = all(v["identical_bytes"] for v in cmp["solver_config"].values())
    items = {it: {"identical": v["n_diff_datasets"] == 0, "n_datasets": v["n_datasets"], "diff": v["diff"]}
             for it, v in h5["by_item"].items()}
    items["実効設定 (solverConfig.yaml・bcondConfig.yaml)"] = {"identical": cfg_ok,
                                                           "diff": {n: v["diff"] for n, v in cmp["solver_config"].items() if not v["identical_bytes"]}}
    items["HDF5 の属性"] = {"identical": h5["attrs_identical"], "diff": h5["attrs_diff"]}
    ok = all(v["identical"] for v in items.values())
    return {"pass": bool(ok), "items": items}


def run_w3():
    A, B = WSB / "w3_A_legacy", WSB / "w3_B_single_bspline"
    cmp = compare_dirs(A, B)
    j = judge_w3(cmp)
    ha, hb = (json.loads((WSB / f"{n}.hashes.json").read_text()) for n in ("w3_A_legacy", "w3_B_single_bspline"))
    env_same = {k: ha.get(k) == hb.get(k) for k in ("design_dir", "python", "numpy", "scipy", "converter_sha256", "FORGE_ALLOW_UNVERIFIED_SPECIES")}
    extra = compare_dirs(WSB / "base_legacy", A)
    out = {"item": "W3 ソルバ入力の同一性 (判別 A/B、plan §6.0)",
           "arms": {"A": "physical_wall_repr: legacy", "B": "physical_wall_repr: single_bspline"},
           "criterion": "変えるのは geometry.physical_wall_repr だけ。変換後の座標・接続・境界対応・体積・面の幾何・centCoords・wall_dist・"
                        "全初期保存量・実効設定が全項目ビット同一なら支持、1 つでも違えば「入力は変わらない」を棄却して差の経路を特定して止める",
           "verdict": "PASS (全項目ビット同一)" if j["pass"] else "FAIL (入力が違う — 差の経路を特定して止める)",
           "items": j["items"], "same_environment": env_same,
           "problem_copy_diff": _problem_diff(ha["problem_copy"], hb["problem_copy"]),
           "nozzle_msh_text": cmp["nozzle_msh"], "files": cmp["files"], "prepare_info_diff": cmp["prepare_info_diff"],
           "h5_datasets": {k: v for k, v in cmp["h5"]["datasets"].items() if not v["identical"]},
           "extra_base_vs_A": {"note": "参考: #2b の基準 (キー無し) と腕 A (legacy を明示) のソルバ入力",
                               "h5_all_identical": extra["h5"]["all_identical"],
                               "solver_config_identical": all(v["identical_bytes"] for v in extra["solver_config"].values()),
                               "files": extra["files"], "prepare_info_diff_keys": extra["prepare_info_diff"][:20]}}
    (WSB / "W3.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


def _problem_diff(pa, pb):
    la, lb = Path(pa).read_text().splitlines(), Path(pb).read_text().splitlines()
    import difflib
    return [l for l in difflib.unified_diff(la, lb, lineterm="", n=0) if not l.startswith(("---", "+++", "@@"))]


# ----------------------------------------------------------------------------------------------- 壁の作り直し (W1・W2・W4・W5)
def build_walls(problem_yaml: Path):
    """固定入力の問題 (wsb_prepare.py の写し) から prepare_ns と同じ経路で今の物理壁 (legacy) と 1 本の B-spline を作る。"""
    from forge_design.evaluate.runner_axismach import design_chain, integral_delta_r, load_problem, _gam_or_gas
    from forge_design.geometry.wall_axismach import PhysicalNozzleWall, SingleBSplinePhysicalWall
    p = load_problem(problem_yaml)
    d = design_chain(p)
    _, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    pw = p.geometry.get("pw_ramp")
    PW = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]),
                            _gam_or_gas(p), p.cp, offset="radial", delta_r_x=drx,
                            ramp=(None if pw is None else tuple(float(v) for v in pw)))
    return p, d, drx, PW, SingleBSplinePhysicalWall(PW)


W1_TOL = {"r": 1.3e-7, "r1": 1e-7, "r2": 1e-5, "joint_jump_d1_d2": 1e-8, "S6_r2_max_increase": 0.002}
W2_TOL = {"x": 1e-6, "r": 1.3e-7, "kappa": 1e-5}


def judge_w1(B, PW, mono) -> dict:
    """W1: 許容誤差 (密な点 + 区間多項式の極値; B.fit_diag でなく判定側で測り直す)・継ぎ目の跳び・区間内 C⁴ (ノットの重複度)・形のゲート。"""
    from forge_design.geometry.wall_axismach import bspline_piece_limits, bspline_wall_errors, r3_piecewise_exact
    spl = B.spline
    t = np.asarray(spl.t)
    dist, mult = np.unique(t, return_counts=True)
    iv = np.c_[dist[:-1], dist[1:]]
    err = bspline_wall_errors(spl, PW.r, iv)
    emax = {n: float(err[n][0].max()) for n in range(4)}
    xat = {n: float(err[n][1][int(np.argmax(err[n][0]))]) for n in range(4)}
    seg = []
    brk = [B.x_in] + list(B.joints) + [B.x_e]
    for a, b in zip(brk[:-1], brk[1:]):
        m = (iv[:, 0] >= a) & (iv[:, 1] <= b)
        seg.append({"range": [a, b], "n_intervals": int(m.sum()),
                    **{f"max_d{n}": float(err[n][0][m].max()) for n in range(4)}})
    tol_ok = emax[0] <= W1_TOL["r"] and emax[1] <= W1_TOL["r1"] and emax[2] <= W1_TOL["r2"]
    jl, jr = bspline_piece_limits(spl, B.joints, nmax=3)
    jumps = [{"x": float(x), **{f"d{n}": float(abs(jr[n, q] - jl[n, q])) for n in range(4)}} for q, x in enumerate(B.joints)]
    jump_ok = all(j["d1"] <= W1_TOL["joint_jump_d1_d2"] and j["d2"] <= W1_TOL["joint_jump_d1_d2"] for j in jumps)
    # 区間の中は 4 階微分まで連続 = 継ぎ目以外の内部ノットの重複度がすべて 1 (5 次の B-spline は単純ノットで C⁴)。構造で判定し、数値の跳びは記録
    inner = (dist > B.x_in) & (dist < B.x_e)
    is_joint = np.isin(dist, np.asarray(B.joints))
    simple = inner & ~is_joint
    struct_ok = bool(np.all(mult[simple] == 1) and np.all(mult[is_joint] == 3) and mult[0] == 6 and mult[-1] == 6)
    xs_simple = dist[simple]
    sl, sr = bspline_piece_limits(spl, xs_simple, nmax=4)
    rel = np.abs(sr - sl) / np.maximum(1.0, np.maximum(np.abs(sl), np.abs(sr)))
    c4 = {f"max_rel_jump_d{n}": float(np.nanmax(rel[n])) for n in range(5)}
    # 形のゲート: S6 (物理壁 r″ の [mono] での最大増加 ≤ 0.002・validate 空)、ランプのゲート (1 本の B-spline で評価し直したもの)
    m6 = r3_piecewise_exact(spl, None, *mono)
    m6_legacy = r3_piecewise_exact(lambda xq, k: PW.r(xq, k), np.r_[PW.design._spl.t, PW._dr.spline.t], *mono)
    v = B.validate()
    s6_ok = m6["r2_max_increase"] <= W1_TOL["S6_r2_max_increase"] and v == []
    # ランプのゲートは判定側で 1 本の B-spline から評価し直す (壁の属性 ramp_gate を信用しない; 式と閾値は元の壁と同じ)
    lo, hi = B._ramp
    xg = np.linspace(lo, hi, 6001)
    rg = {"ramp": [lo, hi], "max_abs_d2_change": float(np.max(np.abs(spl(xg, 2) - B.design.r(xg, 2)))),
          "max_r1": float(np.max(spl(xg, 1))), "limit_d2": float(PW.ramp_gate["limit_d2"]), "source_wall": PW.ramp_gate}
    ramp_ok = bool(rg["max_abs_d2_change"] <= rg["limit_d2"] and rg["max_r1"] < 0.0)
    ok = tol_ok and jump_ok and struct_ok and s6_ok and ramp_ok
    return {"pass": bool(ok), "tol": dict(W1_TOL),
            "max_err": {f"d{n}": emax[n] for n in range(4)}, "x_at_max_err": {f"d{n}": xat[n] for n in range(4)},
            "tol_pass": bool(tol_ok), "segments": seg,
            "joint_jumps": jumps, "joint_jump_pass": bool(jump_ok),
            "c4_structure_pass": struct_ok, "n_simple_interior_knots": int(simple.sum()), "joint_mults": mult[is_joint].tolist(),
            "c4_numeric_record": c4,
            "S6": {"pass": bool(s6_ok), "interval": list(mono), "r2_max_increase": m6["r2_max_increase"], "r2_max": m6["r2_max"],
                   "x_r2_max": m6["x_r2_max"], "limit": W1_TOL["S6_r2_max_increase"], "validate": v,
                   "legacy_wall": {"r2_max_increase": m6_legacy["r2_max_increase"], "r2_max": m6_legacy["r2_max"]}},
            "ramp_gate": {"pass": ramp_ok, **{k: rg[k] for k in ("ramp", "max_abs_d2_change", "max_r1", "limit_d2")},
                          "legacy_wall": {k: rg["source_wall"][k] for k in ("max_abs_d2_change", "max_r1")}},
            "counts": {"n_coef": int(len(spl.c)), "n_knots": int(len(t)), "n_distinct_knots": int(len(dist)),
                       "min_knot_gap_rt": float(np.diff(dist).min()), "joints": list(B.joints)},
            "fit_record": {k: B.fit_diag[k] for k in ("lstsq_rank", "lstsq_cond", "n_gauss", "n_rows", "delta_r_table_range")}}


def judge_w2(B, PW) -> dict:
    """W2: スロート (r′ = 0 の根とそこでの r・r″) を判定側で 1 本の B-spline から求め直し (囲い込みは今の (−0.3, 0.2))、
    元の壁の値との差を許容差と比べる。壁の属性 (x_throat ほか) が判定側の値と一致することも確かめる。"""
    from scipy.optimize import brentq
    from forge_design.geometry.wall_axismach import THROAT_BRACKET
    spl = B.spline
    xl, xr = THROAT_BRACKET
    xt = float(brentq(lambda x: float(spl(x, 1)), xl, xr, xtol=1e-14))
    th = {"x": xt, "r": float(spl(xt)), "kappa": float(spl(xt, 2))}
    d = {"x": th["x"] - PW.x_throat, "r": th["r"] - PW.r_throat, "kappa": th["kappa"] - PW.kappa_throat}
    attr = {"x": B.x_throat - th["x"], "r": B.r_throat - th["r"], "kappa": B.kappa_throat - th["kappa"]}
    attr_ok = abs(attr["x"]) <= 1e-12 and abs(attr["r"]) <= 1e-14 and abs(attr["kappa"]) <= 1e-10
    ok = all(abs(d[k]) <= W2_TOL[k] for k in d) and attr_ok
    return {"pass": bool(ok), "tol": dict(W2_TOL), "diff": d, "single_bspline": th,
            "legacy": {"x": PW.x_throat, "r": PW.r_throat, "kappa": PW.kappa_throat},
            "wall_attr_minus_recomputed": attr, "wall_attr_consistent": bool(attr_ok)}


def run_w1w2():
    ha = json.loads((WSB / "w3_B_single_bspline.hashes.json").read_text())
    p, d, drx, PW, B = build_walls(Path(ha["problem_copy"]))
    mono = tuple(float(v) for v in p.geometry["wall_fit_mono_r2"])
    j1 = judge_w1(B, PW, mono)
    # 判定した壁が W3 腕 B の run が保存した壁と同じこと (係数・ノットのビット一致)
    from forge_design.geometry.wall_axismach import load_wall_file
    Wf = load_wall_file(WSB / "w3_B_single_bspline")
    same = bool(np.array_equal(Wf["physical"].spline.t, B.spline.t) and np.array_equal(Wf["physical"].spline.c, B.spline.c))
    out1 = {"item": "W1 作り直しの精度 (plan §6.0、CFD 0 step)", "verdict": "PASS" if j1["pass"] else "FAIL", **j1,
            "same_as_run_wall_file": same, "run_wall_file": Wf["path"],
            "achieved_vs_limit_orders": {f"d{n}": float(np.log10(lim / max(j1["max_err"][f"d{n}"], 1e-300)))
                                         for n, lim in ((0, W1_TOL["r"]), (1, W1_TOL["r1"]), (2, W1_TOL["r2"]))}}
    (WSB / "W1.json").write_text(json.dumps(out1, indent=1, ensure_ascii=False, default=str))
    j2 = judge_w2(B, PW)
    out2 = {"item": "W2 スロート量 (plan §6.0)", "verdict": "PASS" if j2["pass"] else "FAIL", **j2,
            "scale_m": float(p.spec["r_throat"]),
            "diff_um": {"x": j2["diff"]["x"] * float(p.spec["r_throat"]) * 1e6, "r": j2["diff"]["r"] * float(p.spec["r_throat"]) * 1e6}}
    (WSB / "W2.json").write_text(json.dumps(out2, indent=1, ensure_ascii=False, default=str))
    return out1, out2


# ----------------------------------------------------------------------------------------------- W4
# §4.2 の属性の表: 物理壁の属性 → 下流の読み出し箇所 (ファイルと探す文字列; 行番号は実行時に探す)
ATTR_READERS = {
    "x_in": [("design/forge_design/meshing/mesh2d.py", "_x_stations(wall.x_in")],
    "x_e": [("design/forge_design/meshing/mesh2d.py", "wall.x_e"), ("design/forge_design/evaluate/runner_axismach.py", "np.linspace(wall.x_throat, wall.x_e")],
    "r(x, deriv)": [("design/forge_design/meshing/mesh2d.py", "rw = wall.r(xs)"), ("design/forge_design/evaluate/ic.py", "wall.r(xn)"),
                    ("design/forge_design/evaluate/runner_axismach.py", "wall.r(xs_p)")],
    "validate()": [("design/forge_design/evaluate/runner_axismach.py", "msgs = wall.validate()")],
    "x_throat": [("design/forge_design/evaluate/ic.py", "float(wall.x_throat)"), ("design/forge_design/evaluate/runner_axismach.py", '"x": wall.x_throat')],
    "r_throat": [("design/forge_design/evaluate/ic.py", "float(wall.r_throat)"), ("design/forge_design/evaluate/runner_axismach.py", '"r": wall.r_throat')],
    "kappa_throat": [("design/forge_design/evaluate/runner_axismach.py", "wall.kappa_throat")],
    "_dstar_hist": [("design/forge_design/evaluate/runner_axismach.py", "wall._dstar_hist(0.0)")],
    "offset_mode": [("design/forge_design/evaluate/runner_axismach.py", "wall.offset_mode")],
    "ramp_gate": [("design/forge_design/evaluate/runner_axismach.py", 'getattr(wall, "ramp_gate", None)')],
    "theta(x) [互換の API; 下流の読み出しなし、定義の位置]": [("design/forge_design/geometry/wall_axismach.py", "def theta(self, x)")],
    "design / r_U / L_U / _herm_x0 / _ramp / _throat_diag / _delta_r_applied [壁の内部・validate・診断; 定義の位置]": [
        ("design/forge_design/geometry/wall_axismach.py", "self._herm_x0 = -self.L_U")],
}


def _find_lines(rel: str, pat: str) -> list:
    path = C.parents[1] / rel
    return [f"{rel}:{i + 1}" for i, l in enumerate(path.read_text().splitlines()) if pat in l]


def run_w4():
    import copy
    from forge_design.evaluate import ic as ic_mod
    from forge_design.evaluate.ic import paste_isentropic_ic
    from forge_design.geometry.wall_axismach import (WALL_FILE, SingleBSplinePhysicalWall, check_required_attrs, load_wall_file)
    from forge_design.report.nozzle_report import fig_wall_shape
    import h5py
    ha = json.loads((WSB / "w3_B_single_bspline.hashes.json").read_text())
    p, d, drx, PW, B = build_walls(Path(ha["problem_copy"]))
    S = float(p.spec["r_throat"])
    RB, RA = WSB / "w3_B_single_bspline", WSB / "w3_A_legacy"
    out = {"item": "W4 下流の道具 (plan §6.0)", "checks": {}}
    ck = out["checks"]
    # (1) 属性の表
    table = {}
    for attr, readers in ATTR_READERS.items():
        names = [a.split("(")[0].split("[")[0].strip() for a in attr.split(" / ")]
        present = {n: (getattr(B, n, None) is not None) for n in names}
        table[attr] = {"readers": sum((_find_lines(r, pt) for r, pt in readers), []), "present_non_none": present}
    miss = [a for a in SingleBSplinePhysicalWall.REQUIRED_ATTRS if getattr(B, a, None) is None]
    ck["attribute_table"] = {"pass": (not miss) and all(all(v["present_non_none"].values()) and v["readers"] for v in table.values()),
                             "table": table, "required_attrs": list(SingleBSplinePhysicalWall.REQUIRED_ATTRS), "missing": miss}
    # (2) prepare_ns が新しい壁を使った (腕 B の prepare_info・壁ファイル・wall_physical.csv)
    ib, ia = json.loads((RB / "prepare_info.json").read_text()), json.loads((RA / "prepare_info.json").read_text())
    Wf = load_wall_file(RB)
    pwi = ib.get("physical_wall") or {}
    wp = np.loadtxt(RB / "wall_physical.csv", delimiter=",", skiprows=1)
    ck["prepare_ns"] = {
        "physical_wall_repr": pwi.get("repr"), "wall_file": pwi.get("file"),
        "wall_file_sha256_matches": pwi.get("sha256") == sha_file(RB / WALL_FILE),
        "throat_physical_is_new_wall": ib["throat_physical"]["x"] == Wf["throat"]["x"] == B.x_throat,
        "throat_physical_A_minus_B": {k: ia["throat_physical"][k] - ib["throat_physical"][k] for k in ("x", "r", "kappa")},
        "pw_ramp_gate_repr": (ib.get("pw_ramp_gate") or {}).get("repr"),
        "wall_physical_csv_from_new_wall_max_abs_m": float(np.abs(Wf["physical"].r(np.clip(wp[:, 0] / S, *Wf["domain"])) * S - wp[:, 1]).max())}
    ck["prepare_ns"]["pass"] = bool(ck["prepare_ns"]["physical_wall_repr"] == "single_bspline" and ck["prepare_ns"]["wall_file_sha256_matches"]
                                    and ck["prepare_ns"]["throat_physical_is_new_wall"] and ck["prepare_ns"]["pw_ramp_gate_repr"] == "single_bspline"
                                    and ck["prepare_ns"]["wall_physical_csv_from_new_wall_max_abs_m"] <= 1e-15)
    # (3) 初期値 (ic.py) が新しい壁で動く: 腕 B の nozzle.h5 の写しに作り直した 1 本の壁で IC を貼り直し、保存された IC と比べる
    with tempfile.TemporaryDirectory() as td:
        h5c = Path(td) / "nozzle.h5"
        shutil.copy(RB / "nozzle.h5", h5c)
        from forge_design.evaluate.runner_axismach import _tp_species_Y
        paste_isentropic_ic(h5c, B, S, float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                            gas=(None if str(p.evaluate.get("cfd_gas", "same")) == "cpg" else p.gas_model),
                            h_ref_T=float(p.evaluate.get("thermo_href_temp", 298.15)), species_Y=_tp_species_Y(p))
        same = {}
        with h5py.File(h5c) as f1, h5py.File(RB / "nozzle.h5") as f2:
            for k in f2["/VALUE"]:
                same[k] = bool(np.array_equal(f1[f"/VALUE/{k}"][:], f2[f"/VALUE/{k}"][:]))
    Bx = copy.copy(B)
    del Bx.x_throat
    try:
        ic_mod._throat_of(Bx)
        ic_raise = False
    except ValueError:
        ic_raise = True
    ck["initial_field"] = {"pass": bool(all(same.values()) and ic_raise and ic_mod._throat_of(B) == (B.x_throat, B.r_throat)),
                           "repaste_identical": same, "missing_throat_attr_raises": ic_raise}
    # (4) 報告の単体検査 4 種 (fig_wall_shape は F["S"] だけを読む; 結果ファイル不要)
    F = {"S": S}
    rep = {}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        o_old = fig_wall_shape(WSB / "base_legacy", F, td / "old.png")
        rep["old_run_no_wall_file"] = {"pass": o_old.get("wall_source") == "legacy_csv" and "旧経路" in o_old.get("wall_source_note", ""),
                                       "run": str(WSB / "base_legacy"), "wall_source": o_old.get("wall_source"), "note": o_old.get("wall_source_note")}
        o_new = fig_wall_shape(RB, F, td / "new.png")
        rep["new_format_saved_coefficients"] = {
            "pass": o_new.get("wall_source") == "saved_coefficients" and o_new["exit_radius_m"] == float(Wf["physical"].r(np.array([Wf["domain"][1]]))[0] * S),
            "run": str(RB), "wall_source": o_new.get("wall_source"), "exit_radius_m": o_new.get("exit_radius_m"),
            "r2_highfreq_max_x_gt2": {"saved": o_new["r2_highfreq_max_x_gt2"], "legacy_csv_same_run": None}}
        # 同じ run の CSV 経路の値 (壁ファイルを外した写しで) と並べる
        cp = td / "csv_only"; cp.mkdir()
        for n in ("wall_physical.csv", "wall_design.csv"):
            shutil.copy(RB / n, cp / n)
        rep["new_format_saved_coefficients"]["r2_highfreq_max_x_gt2"]["legacy_csv_same_run"] = fig_wall_shape(cp, F, td / "c.png")["r2_highfreq_max_x_gt2"]
        bad = td / "bad"; bad.mkdir()
        r_ = json.loads((RB / WALL_FILE).read_text()); r_["physical_wall"].pop("c")
        (bad / WALL_FILE).write_text(json.dumps(r_))
        try:
            fig_wall_shape(bad, F, td / "bad.png"); raised = None
        except ValueError as e:
            raised = str(e)
        rep["missing_coefficients"] = {"pass": raised is not None, "error": raised}
        from scipy.interpolate import BSpline
        spl_info = ib["wall_fit"]["spline"]
        ext = float(BSpline(np.asarray(spl_info["t"]), np.asarray(spl_info["c"]), int(spl_info["k"]))(-6.0))
        h6 = float(Wf["design"].r(np.array([-6.0]))[0])
        try:
            Wf["design"].S(np.array([-6.0])); s_raise = False
        except ValueError:
            s_raise = True
        rep["upstream_difference_plot"] = {"pass": h6 == float(d["wall"].up.r(np.array([-6.0]))[0]) and s_raise and abs(h6 - ext) > 1.0,
                                           "design_at_x_minus6_restored": h6, "hermite_at_x_minus6": float(d["wall"].up.r(np.array([-6.0]))[0]),
                                           "prepare_info_wall_fit_spline_extrapolated_at_x_minus6": ext, "S_outside_domain_raises": s_raise}
        shutil.copy(td / "new.png", WSB / "W4_fig_wall_shape_saved.png")
        shutil.copy(td / "old.png", WSB / "W4_fig_wall_shape_legacy_csv.png")
    ck["report_unit"] = {"pass": all(v["pass"] for v in rep.values()), **rep}
    # (5) 統合検査: 既存の結果 (run_0147 の nozzle.h5・res_*.h5) はローカルに無い
    hb = json.loads((WSB / "base_legacy.hashes.json").read_text())
    ck["report_integration"] = {"pass": None, "status": "未実施 (判定不能)", "reason": "run_0147 の nozzle.h5・res_*.h5 がローカルに無い (AWS、未取得)",
                                "source": hb.get("run_0147_results")}
    # (6) 一般性: 違う scale_m・L_U・L_pipe・既定ランプ (design/tests/run_wall_single_bspline_tests.py §2 と同じ設定)
    from forge_design.evaluate.runner_axismach import design_chain, integral_delta_r, load_problem, _gam_or_gas
    from forge_design.geometry.wall_axismach import PhysicalNozzleWall, default_pw_ramp, save_wall_file
    from forge_design.export.wall_step import step_curve_data
    pg = load_problem(C / PN)
    g = pg.geometry
    g["initial_line"] = "hall"; g.pop("initial_line_run"); g.pop("initial_line_res"); g.pop("pw_ramp")
    g["L_U"], g["L_pipe"], g["n_axis_inv"] = 10.0, 1.0, 800
    pg.spec["r_throat"] = 0.05
    dg = design_chain(pg)
    _, drg, _ = integral_delta_r(pg, dg, pg.raw["deltastar_initializer"])
    PWg = PhysicalNozzleWall(dg["wall"], dg["wall_inv"], 0.05, float(pg.spec["Pt"]), float(pg.spec["Tt"]), _gam_or_gas(pg), pg.cp,
                             offset="radial", delta_r_x=drg, ramp=None)
    Bg = SingleBSplinePhysicalWall(PWg)
    lo_d, _ = default_pw_ramp(dg["wall"])
    with tempfile.TemporaryDirectory() as td:
        save_wall_file(td, Bg, 0.05, "single_bspline")
        dat = step_curve_data(load_wall_file(td)["record"])
    ck["generality"] = {"pass": bool(Bg.joints == [-10.0, lo_d, -5.0, 0.0] and Bg.x_in == -11.0 and dat["scale_mm_per_rt"] == 50.0
                                     and np.array_equal(dat["t_mm"], np.asarray(Bg.spline.t) * 50.0)),
                        "setting": "case/45 の問題の初期線を Hall に、r_throat 0.05 m・L_U 10・L_pipe 1・pw_ramp 無し (既定ランプ)・n_axis_inv 800",
                        "joints": Bg.joints, "expected_joints": [-10.0, lo_d, -5.0, 0.0], "domain": [Bg.x_in, Bg.x_e],
                        "ramp_source": PWg._ramp_source, "scale_mm_per_rt": dat["scale_mm_per_rt"], "max_err": Bg.fit_diag["max_err"]}
    unit_ok = all(ck[k]["pass"] for k in ("attribute_table", "prepare_ns", "initial_field", "report_unit", "generality"))
    out["verdict"] = ("PASS (単体・一般性) / 統合検査は未実施 (判定不能: run_0147 の結果がローカルに無い)" if unit_ok
                      else "FAIL")
    (WSB / "W4.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


# ----------------------------------------------------------------------------------------------- W5
def judge_w5_orig(data, rd, pts, B, PW, s_mm):
    """W5「元の壁との差」: STEP (OpenCascade の評価) と元の物理壁の関数を同じ物理 x で比べる。予算 = 表現の誤差 + 転送誤差。
    表現の誤差は保存した 1 本の B-spline と元の壁の接線角・曲率を同じ点で直接比べた最大 (r″ の許容差の換算ではない)。
    半径の予算は W1 の上限 (1.3e-7 r_t を mm に換算) + 1e-6 mm。接線・曲率は OCC の u 微分から (u = x [mm])。"""
    from forge_design.export.wall_step import TRANSFER_TOL, angle_diff, tangent_angle_curvature
    x_in, x_e = B.x_in, B.x_e
    rows = []
    for (u, side, lab), (p0, d1, d2) in zip(pts, rd["eval"]):
        x = min(max(p0[0] / s_mm, x_in), x_e)
        a_s, k_s = tangent_angle_curvature(np.array(d1[:2]), np.array(d2[:2]))
        r0, r1, r2 = (float(PW.r(np.array([x]), n)[0]) for n in range(3))
        b0, b1, b2 = (float(B.spline(x, n)) for n in range(3))
        a_o, a_b = np.arctan(r1), np.arctan(b1)
        k_o, k_b = (r2 / s_mm) / (1 + r1 * r1) ** 1.5, (b2 / s_mm) / (1 + b1 * b1) ** 1.5
        rows.append({"u": u, "x": x, "dr_step_orig_mm": p0[1] - r0 * s_mm,
                     "dang_step_orig": float(abs(angle_diff(a_s, a_o))), "dkap_step_orig": float(abs(k_s - k_o)),
                     "dang_repr": float(abs(angle_diff(a_b, a_o))), "dkap_repr": float(abs(k_b - k_o)), "kap_orig": k_o,
                     # 単位の確認: OCC の dy/du (無次元) = r′(x)、d²y/du² [1/mm] = r″(x) [1/r_t] / s
                     "d1_occ_minus_r1": d1[1] - b1, "d2_occ_minus_r2_over_s": d2[1] - b2 / s_mm, "dx_du_minus_1": d1[0] - 1.0, "d2x_du2": d2[0]})
    dr = np.array([abs(r["dr_step_orig_mm"]) for r in rows])
    da, dk = np.array([r["dang_step_orig"] for r in rows]), np.array([r["dkap_step_orig"] for r in rows])
    ra, rk = np.array([r["dang_repr"] for r in rows]), np.array([r["dkap_repr"] for r in rows])
    budget = {"radius_mm": 1.3e-7 * s_mm + TRANSFER_TOL["pos_mm"], "angle_rad": float(ra.max()) + TRANSFER_TOL["angle_rad"],
              "kappa_per_mm": float(rk.max()) + TRANSFER_TOL["kappa_abs_per_mm"]}
    kap_o = np.array([abs(r["kap_orig"]) for r in rows])
    kap_ok = (dk <= budget["kappa_per_mm"]) | (dk <= float(rk.max()) + TRANSFER_TOL["kappa_rel"] * kap_o)
    units = {"max_abs_d1_occ_minus_r1": float(max(abs(r["d1_occ_minus_r1"]) for r in rows)),
             "max_abs_d2_occ_minus_r2_over_s_per_mm": float(max(abs(r["d2_occ_minus_r2_over_s"]) for r in rows)),
             "max_abs_dx_du_minus_1": float(max(abs(r["dx_du_minus_1"]) for r in rows)),
             "max_abs_d2x_du2": float(max(abs(r["d2x_du2"]) for r in rows))}
    ok = bool(dr.max() <= budget["radius_mm"] and da.max() <= budget["angle_rad"] and kap_ok.all())
    return {"pass": ok, "budget": budget, "n_points": len(rows),
            "radius_max_mm": float(dr.max()), "angle_max_rad": float(da.max()), "kappa_max_per_mm": float(dk.max()),
            "repr_angle_max_rad": float(ra.max()), "repr_kappa_max_per_mm": float(rk.max()),
            "kappa_points_failing": int((~kap_ok).sum()), "derivative_units_check": units,
            "derivative_units_pass": bool(units["max_abs_d1_occ_minus_r1"] <= 1e-9 and units["max_abs_d2_occ_minus_r2_over_s_per_mm"] <= 1e-9
                                          and units["max_abs_dx_du_minus_1"] <= 1e-9)}


def run_w5():
    import re
    from forge_design.export.wall_step import export_run, read_step, sample_params, step_curve_data, transfer_check
    from forge_design.geometry.wall_axismach import load_wall_file
    RB = WSB / "w3_B_single_bspline"
    outd = WSB / "step"
    if outd.exists():
        shutil.rmtree(outd)
    sc = export_run(RB, outd, n_per_interval=3)
    W = load_wall_file(RB)
    data = step_curve_data(W["record"])
    pts = sample_params(data, 3)
    rd = read_step(outd / "wall_physical.step", [u for u, _, _ in pts], revolve=False)
    ha = json.loads((WSB / "w3_B_single_bspline.hashes.json").read_text())
    _, _, _, PW, B = build_walls(Path(ha["problem_copy"]))
    same = bool(np.array_equal(B.spline.c, W["physical"].spline.c) and np.array_equal(B.spline.t, W["physical"].spline.t))
    orig = judge_w5_orig(data, rd, pts, B, PW, data["scale_mm_per_rt"])
    tr = sc["readback"]["transfer"]
    txt = (outd / "wall_physical.step").read_text()
    reals = re.findall(r"-?\d+\.\d*(?:E[+-]\d+)?", txt.split("B_SPLINE_CURVE_WITH_KNOTS")[0])
    sig = max(len(re.sub(r"[-.]|E.*", "", v).lstrip("0")) for v in reals) if reals else None
    # 回転面の面積と解析値 2π∫ r √(1 + r′²) dx の比較 (記録)
    gx, gw = np.polynomial.legendre.leggauss(8)
    kn = data["knots_mm"]
    a, b = kn[:-1], kn[1:]
    xq = (0.5 * (a + b))[:, None] + (0.5 * (b - a))[:, None] * gx[None, :]
    s = data["scale_mm_per_rt"]
    rq, r1q = B.spline(xq / s) * s, B.spline(xq / s, 1)
    area = float((2 * np.pi * rq * np.sqrt(1 + r1q ** 2) * (0.5 * (b - a))[:, None] * gw[None, :]).sum())
    rv = sc["readback"]["revolve"] or {}
    out = {"item": "W5 STEP (plan §6.0)", "step": str(outd / "wall_physical.step"), "sidecar": str(outd / "wall_physical_step.json"),
           "same_as_rebuilt_wall": same,
           "structure": tr["structure"], "structure_pass": all(tr["structure"].values()),
           "transfer": {k: tr[k] for k in tr if k != "structure"},
           "orig": orig,
           "step_real_significant_digits": sig,
           "revolve": {**rv, "pass": bool(rv.get("is_valid") is True and rv.get("area_mm2", 0) > 0),
                       "area_analytic_mm2": area, "area_rel_diff": (rv.get("area_mm2", 0) - area) / area},
           "receiving_cad": "未確認 (受け取り側の CAD の種類・版での読み込みはユーザ側で確認; plan §6 W5)",
           "freecad_version": sc["readback"].get("freecad_version"),
           "chord_minus_curve_um": sc["cad_vs_cfd"].get("chord_minus_curve")}
    ok = out["structure_pass"] and tr["pass"] and orig["pass"] and orig["derivative_units_pass"] and out["revolve"]["pass"]
    out["verdict"] = "PASS" if ok else "FAIL"
    (WSB / "W5.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    return out


# ----------------------------------------------------------------------------------------------- 負例 (判定が FAIL を FAIL と出すか)
def run_neg():
    import copy
    import h5py
    from scipy.interpolate import BSpline, insert
    from forge_design.export.wall_step import read_step, sample_params, step_curve_data, transfer_check, write_step
    from forge_design.geometry.wall_axismach import WALL_FILE, PhysicalNozzleWall, SingleBSplinePhysicalWall, load_wall_file
    from forge_design.evaluate.runner_axismach import _gam_or_gas, delta_r_from_table
    from forge_design.export.wall_step import greville
    from forge_design.report.nozzle_report import fig_wall_shape
    res = {}
    ha = json.loads((WSB / "w3_B_single_bspline.hashes.json").read_text())
    p, d, drx, PW, B = build_walls(Path(ha["problem_copy"]))
    mono = tuple(float(v) for v in p.geometry["wall_fit_mono_r2"])
    with tempfile.TemporaryDirectory(dir=str(Path(os.environ.get("TMPDIR", "/tmp")))) as td:
        td = Path(td)
        # W3: 腕 B の写しの VALUE/ro を 1 要素だけ 1 ulp ずらす / solverConfig.yaml に空白を 1 つ足す
        A, Bd = WSB / "w3_A_legacy", WSB / "w3_B_single_bspline"
        b1 = td / "B_ulp"; shutil.copytree(Bd, b1)
        with h5py.File(b1 / "nozzle.h5", "r+") as f:
            v = f["/VALUE/ro"][:]; v[1234] = np.nextafter(v[1234], np.float32(np.inf)); f["/VALUE/ro"][:] = v
        j = judge_w3(compare_dirs(A, b1))
        res["W3_ro_one_ulp"] = {"detected": not j["pass"], "why": [k for k, v in j["items"].items() if not v["identical"]]}
        b2 = td / "B_cfg"; shutil.copytree(Bd, b2)
        (b2 / "solverConfig.yaml").write_text((b2 / "solverConfig.yaml").read_text() + " ")
        j = judge_w3(compare_dirs(A, b2))
        res["W3_solverConfig_one_byte"] = {"detected": not j["pass"], "why": [k for k, v in j["items"].items() if not v["identical"]]}
        b3 = td / "B_wd"; shutil.copytree(Bd, b3)
        with h5py.File(b3 / "nozzle.h5", "r+") as f:
            v = f["/VALUE/wall_dist"][:]; v[-1] = np.nextafter(v[-1], np.float32(np.inf)); f["/VALUE/wall_dist"][:] = v
        j = judge_w3(compare_dirs(A, b3))
        res["W3_wall_dist_one_ulp"] = {"detected": not j["pass"], "why": [k for k, v in j["items"].items() if not v["identical"]]}
        # W0: キー無しの写しに壁ファイルが増えた / prepare_info.json が 1 バイト違う
        n1 = td / "nokey_extra"; shutil.copytree(WSB / "after_nokey", n1)
        shutil.copy(Bd / WALL_FILE, n1 / WALL_FILE)
        j = judge_w0(compare_dirs(WSB / "base_legacy", n1))
        res["W0_extra_wall_file"] = {"detected": not j["pass"], "why": j["files_only_in_one"]}
        n2 = td / "nokey_info"; shutil.copytree(WSB / "after_nokey", n2)
        (n2 / "prepare_info.json").write_text((n2 / "prepare_info.json").read_text() + "\n")
        j = judge_w0(compare_dirs(WSB / "base_legacy", n2))
        res["W0_prepare_info_one_byte"] = {"detected": not j["pass"], "why": j["files_bytes_diff_except_h5"]}
        # W1: 係数を 1 個ずらす (値に 1e-6 を足す / 並びを 1 つずらす)・単純ノットを重複させる (同じ関数で重複度だけ 2)
        def with_spline(spl):
            Bx = copy.copy(B); Bx.spline = spl
            return Bx
        c = np.array(B.spline.c); c[900] += 1e-6
        def w1why(j):
            return {"tol_pass": j["tol_pass"], "max_err": j["max_err"], "joint_jump_pass": j["joint_jump_pass"],
                    "c4_structure_pass": j["c4_structure_pass"], "S6": j["S6"]["pass"], "ramp_gate": j["ramp_gate"]["pass"]}
        j = judge_w1(with_spline(BSpline(B.spline.t, c, 5)), PW, mono)
        res["W1_coef_plus_1e-6"] = {"detected": (not j["pass"]) and (not j["tol_pass"]), "why": w1why(j)}
        j = judge_w1(with_spline(BSpline(B.spline.t, np.roll(B.spline.c, 1), 5)), PW, mono)
        res["W1_coef_shift_by_one"] = {"detected": (not j["pass"]) and (not j["tol_pass"]), "why": w1why(j)}
        tk = np.asarray(B.spline.t); xk = float(np.unique(tk)[700])
        t2, c2, _ = insert(xk, (tk, np.asarray(B.spline.c), 5))
        j_ins = judge_w1(with_spline(BSpline(t2, c2[:len(t2) - 6], 5)), PW, mono)
        res["W1_duplicate_simple_knot"] = {"detected": (not j_ins["pass"]) and (not j_ins["c4_structure_pass"]) and j_ins["tol_pass"],
                                           "why": w1why(j_ins), "knot_x": xk}
        # W2: スロートの属性を 2e-6 ずらす / スロート近傍の係数を 1e-5 ずらす
        Bt = copy.copy(B); Bt.x_throat = B.x_throat + 2e-6
        j = judge_w2(Bt, PW)
        res["W2_attr_shift"] = {"detected": (not j["pass"]) and (not j["wall_attr_consistent"]), "why": j["wall_attr_minus_recomputed"]}
        i_th = int(np.argmin(np.abs(greville(B.spline.t, 5) - B.x_throat)))     # スロートに最も近いグレビル点の係数
        c3 = np.array(B.spline.c); c3[i_th] += 1e-5
        Bc = with_spline(BSpline(B.spline.t, c3, 5))
        j = judge_w2(Bc, PW)
        res["W2_coef_near_throat"] = {"detected": (not j["pass"]) and any(abs(j["diff"][k]) > W2_TOL[k] for k in j["diff"]),
                                      "why": j["diff"], "coef_index": i_th}
        # δ_r の表を短くする → 1 本の B-spline を作らない (例外)
        x_tab = np.linspace(-12.499999, B.x_e - 2.0, 1400)
        short = delta_r_from_table(x_tab, drx(x_tab))
        PWs = PhysicalNozzleWall(PW.design, d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]),
                                 _gam_or_gas(p), p.cp, offset="radial", delta_r_x=short, ramp=PW._ramp)
        try:
            SingleBSplinePhysicalWall(PWs); res["delta_r_table_short_raises"] = False
        except ValueError:
            res["delta_r_table_short_raises"] = True
        # 壁ファイルの要素を欠く → 読み込み・報告が例外
        for key in ("c", "t", "throat"):
            r_ = json.loads((Bd / WALL_FILE).read_text()); r_["physical_wall"].pop(key)
            q = td / f"miss_{key}"; q.mkdir(); (q / WALL_FILE).write_text(json.dumps(r_))
            try:
                load_wall_file(q); res[f"wall_file_missing_{key}_raises"] = False
            except ValueError:
                res[f"wall_file_missing_{key}_raises"] = True
        try:
            fig_wall_shape(td / "miss_c", {"S": float(p.spec["r_throat"])}, td / "x.png"); res["report_missing_c_raises"] = False
        except ValueError:
            res["report_missing_c_raises"] = True
        # W5: 制御点 1 個を 1e-5 mm ずらした STEP / ノット 1 個を 1e-5 mm ずらした STEP を、ずらす前の曲線と比べる
        W = load_wall_file(Bd)
        data = step_curve_data(W["record"])
        pts = sample_params(data, 1)
        for lab, mut in (("pole", lambda dd: dd["poles_mm"].__setitem__((800, 1), dd["poles_mm"][800, 1] + 1e-5)),
                         ("knot", lambda dd: dd["knots_mm"].__setitem__(900, dd["knots_mm"][900] + 1e-5))):
            dd = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v) for k, v in data.items()}
            mut(dd)
            stp = td / f"neg_{lab}.step"
            write_step(dd, stp)
            rd = read_step(stp, [u for u, _, _ in pts])
            j = transfer_check(data, rd, pts)
            res[f"W5_{lab}_shift_1e-5mm"] = {"detected": not j["pass"], "why": {"structure": j["structure"], "pos_max_mm": j.get("pos_max_mm"),
                                                                              "angle_max_rad": j.get("angle_max_rad"),
                                                                              "kappa_abs_max_per_mm": j.get("kappa_abs_max_per_mm")}}
    det = {k: (v["detected"] if isinstance(v, dict) else bool(v)) for k, v in res.items()}
    out = {"item": "判定スクリプトの負例 (FAIL を FAIL と出すか; detected = 意図した検査で FAIL・例外になった)", "cases": res,
           "verdict": "PASS (全負例を FAIL・例外として検出)" if all(det.values()) else "FAIL (検出できない負例がある)"}
    (WSB / "negatives.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    todo = ("w0", "w3", "w1w2", "w4", "w5", "neg") if cmd == "all" else (cmd,)
    for c in todo:
        r = globals()[f"run_{c}"]()
        if isinstance(r, tuple):
            for x in r:
                print(c, x["verdict"])
        else:
            print(c, r["verdict"])
