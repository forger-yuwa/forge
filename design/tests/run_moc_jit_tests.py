#!/usr/bin/env python3
"""逆 MOC の単位過程の numba の版 (`moc_fill_nb.pairs_nb`、`FORGE_MOC_KERNEL`) の試験 — plan tooling-nozzle-mknot-fast
§6 V5 (§5.1 #4・#4b、§4.4・§4.5)。numba の版は numpy の版と**ビット同一を求めない** (§4.5、ユーザ決定 2026-10-15) ので、許容差で照合する。

既定 (合成の入力、数十秒):
  1. backend の選び方: 未設定 = numba、`numpy`、他の値 (空を含む) は ValueError。numba が import できなければ ImportError
     (黙って numpy の版に落ちない、別プロセスで numba を隠して確かめる)。単位過程の差し替えを検出して止める。扱わないガスは ValueError
  2. 部品: ν→M (CPG の Newton、前の M から始める) が `pm_mach_vec` と |ΔM| ≤ 1e−12、表の補間が `np.interp` と |ΔM| ≤ 1e−12
     (NaN は NaN、範囲外は端の値、表の節点はそのままの値)
  3. 放射源流の網 (CPG γ 1.4、2 解像度) と semi-perfect の表 (case/71 の燃焼ガス) で、4 つの組 (legacy/analytic × fixed2/converge) と
     `fill_levels`・`fill_arrays` の両方で、網の NaN の位置が同じ、位置 ≤ 1e−10 r_t・θ/ν ≤ 1e−10 rad、集計 (`last_diag`) の整数
     (分類の数・分岐の数・位置の添字・反復回数の分布) が同じ、実数 (位置・残差) は許容差の中。`AXIS_LIMIT_FRAC` を実行時に変えても同じ
  4. 異常系 (plan 段 jit の M2): もともとの欠損・反復中の非有限・平行な特性線・軸より下・反復の上限・θ_r の欠如 (ValueError の文) を
     両方の backend の `fill_levels` に通し、段・対ごとの状態 (`PAIR_*`)・反復回数・源項の分岐が同じ、最終残差が許容差の中、
     棄却の位置・`moc_gate` の判定が同じ (段・対ごとの値は `MocFillDiag.add` の入力を両方の backend で記録して比べる)
--real (V5 の実データ、所要 10〜20 分、最大 RSS 約 3 GB):
  5. M5 (case/70 c2b)・M5.5 (case/71 のカード)・M6 (case/45 の生産の問題、CFD の初期線) を N 600 (L_total_tol 0.1 r_t) と
     2400 (問題の値) で `evaluate_mknot` に通し、逆 MOC の呼び出しごとに網・壁・集計を比べる。合格: 網の有限の点の位置の差 ≤ 1e−10 r_t、
     θ・ν の差 ≤ 1e−10 rad、NaN の位置が同じ、壁の位置の差 ≤ 1e−10 r_t、θ_max の差 ≤ 1e−8°、判定の欄 (ok・reasons・error) が同じ、
     集計 (`moc`) の分類の数が同じ。差の最大を表で出す
  6. 所要: 放射源流 n_ax 2400 の充填と N 2400 の 1 評価 (numpy と numba、コンパイル済み) と、空のキャッシュでのコンパイル込みの初回
     (別プロセス)。判定 (3 倍) は AWS A の値で行う (plan §6 V5)、手元の値は記録

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_moc_jit_tests.py [--real]
凍結源 (case/45 の run_0062) や case/71 のカードが無ければその項目を飛ばして exit 2 (合格扱いにしない)。
環境変数: CASE_RUNS (case/45 の run の場所、既定 /home/sano/work/forge/case/45.isobutane_m6_d155)、
CARD71 (case/71 のカード、既定 /home/sano/work/forge-m55/case/71.m55_combustion_wt/card_m55.yaml)。
"""
import copy
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from forge_design.geometry import moc_fill_nb as NB  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from forge_design.geometry import moc_kernel as MK  # noqa: E402
from forge_design.geometry.moc_inverse import InverseMOC, MocFillDiag, axis_theta_r_init, moc_gate  # noqa: E402
from forge_design.geometry.moc_kernel import (PAIR_BELOW_AXIS, PAIR_CONVERGED, PAIR_DONE, PAIR_GEOM,  # noqa: E402
                                              PAIR_MAXITER, PAIR_MISSING, PAIR_NONFINITE, _Pt, pm_mach_vec, pm_nu)
import moc_axis_limit_radial as RAD  # noqa: E402

FAIL = 0
SKIPPED = []
REAL = "--real" in sys.argv[1:]
G = 1.4
MODES = (("legacy", "fixed2"), ("analytic", "fixed2"), ("legacy", "converge"), ("analytic", "converge"))
C45 = ROOT / "case/45.isobutane_m6_d155"
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
CARD71 = Path(os.environ.get("CARD71", "/home/sano/work/forge-m55/case/71.m55_combustion_wt/card_m55.yaml"))
PY = sys.executable


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def backend(name):
    if name is None:
        os.environ.pop(MI.MOC_KERNEL_ENV, None)
    else:
        os.environ[MI.MOC_KERNEL_ENV] = name


def same_obj(a, b, path=""):
    """入れ子の dict・list・tuple・数の完全一致 (NaN どうしは一致、辞書の順と型も見る)。戻り: 違いの一覧。"""
    if type(a) is not type(b):
        return [f"{path}: 型 {type(a).__name__} vs {type(b).__name__}"]
    if isinstance(a, dict):
        if list(a) != list(b):
            return [f"{path}: キー {list(a)} vs {list(b)}"]
        return [e for k in a for e in same_obj(a[k], b[k], f"{path}.{k}")]
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return [f"{path}: 長さ {len(a)} vs {len(b)}"]
        return [e for i, (x, y) in enumerate(zip(a, b)) for e in same_obj(x, y, f"{path}[{i}]")]
    if isinstance(a, float):
        return [] if (a == b or (a != a and b != b)) else [f"{path}: {a!r} vs {b!r}"]
    if isinstance(a, np.ndarray):
        return [] if (a.shape == b.shape and np.array_equal(a, b, equal_nan=True)) else [f"{path}: 配列が違う"]
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def bits_equal(a, b):
    return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)


TOL_X, TOL_A, TOL_RES = 1e-10, 1e-10, 1e-12      # 位置 [r_t]・θ/ν [rad]・最終残差 (V5)


def close_obj(a, b, atol, path=""):
    """入れ子の dict・list・tuple の照合: 整数・文字列・真偽は完全一致、実数は |差| ≤ atol (NaN どうし・同じ無限は一致)、
    辞書はキーの順まで同じ。戻り: 違いの一覧。"""
    if isinstance(a, (float, np.floating)) and isinstance(b, (float, np.floating)):
        if a == b or (a != a and b != b) or abs(a - b) <= atol:
            return []
        return [f"{path}: {a!r} vs {b!r}"]
    if type(a) is not type(b):
        return [f"{path}: 型 {type(a).__name__} vs {type(b).__name__}"]
    if isinstance(a, dict):
        if list(a) != list(b):
            return [f"{path}: キー {list(a)} vs {list(b)}"]
        return [e for k in a for e in close_obj(a[k], b[k], atol, f"{path}.{k}")]
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return [f"{path}: 長さ {len(a)} vs {len(b)}"]
        return [e for i, (x, y) in enumerate(zip(a, b)) for e in close_obj(x, y, atol, f"{path}[{i}]")]
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def lev_diff(a, b):
    """網 (…, 5) の照合: (形と NaN の位置が同じか, 位置の差の最大, θ・ν の差の最大)。"""
    if a.shape != b.shape:
        return False, np.inf, np.inf
    fa, fb = np.isfinite(a), np.isfinite(b)
    both = fa & fb
    dx = float(np.max(np.abs(a[..., :2] - b[..., :2])[both[..., :2]], initial=0.0))
    da = float(np.max(np.abs(a[..., 2:4] - b[..., 2:4])[both[..., 2:4]], initial=0.0))
    return bool(np.array_equal(fa, fb)), dx, da


class RecDiag(MocFillDiag):
    """段・対ごとの状態 (`add` の中の軸より下の上書きの後)・反復回数・最終残差・分岐を記録する (両方の backend で同じ `add` を通る)。"""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.rec = {}

    def add(self, level, Ax, Ar, Bx, Br, xP, rP, st, below):
        status = np.array(st["status"], dtype=np.int8)
        status[(status == PAIR_CONVERGED) & below] = PAIR_BELOW_AXIS
        status[(status == PAIR_DONE) & below] = PAIR_BELOW_AXIS
        b = st["branch"]
        self.rec[int(level)] = {"status": status, "n_iter": np.asarray(st["n_iter"], dtype=np.int64),
                                "rg": np.asarray(st["resid_geom"], dtype=float), "rc": np.asarray(st["resid_comp"], dtype=float),
                                "branch": np.array([np.asarray(b[k]).astype(np.int8) for k in NB.BRANCH_KEYS])}
        super().add(level, Ax, Ar, Bx, Br, xP, rP, st, below)


class RecMOC(InverseMOC):
    def _new_diag(self):
        self._rec = RecDiag(self.axis_limit, self.corrector, self.n_corr, self.tol, self.max_corr)
        return self._rec


def run_both(init, al, co, thr=None, gas=G, how="fill_levels", **kw):
    """両方の backend で充填 (`fill_levels` か `fill_arrays`)。戻り: {backend: (網, last_diag, 段ごとの記録)}。"""
    out = {}
    for be in ("numpy", "numba"):
        backend(be)
        inv = RecMOC(gamma=gas, axis_limit=al, corrector=co, **kw)
        lev = getattr(inv, how)(init, axis_thr=thr if al == "analytic" else None)
        out[be] = (lev, inv.last_diag, inv._rec.rec)
    backend(None)
    return out


def rec_matches(ra, rb):
    """段ごとの記録の照合: 状態・反復回数・分岐は完全一致、最終残差は NaN の位置が同じで差 ≤ TOL_RES。戻り: 違いの一覧。"""
    errs = []
    if list(ra) != list(rb):
        return [f"段の並び numpy {list(ra)[:5]}… ({len(ra)}) vs numba {list(rb)[:5]}… ({len(rb)})"]
    for k in ra:
        a, b = ra[k], rb[k]
        for key in ("status", "n_iter", "branch"):
            if not np.array_equal(a[key], b[key]):
                errs.append(f"段 {k} の {key}")
        for key in ("rg", "rc"):
            fa, fb = np.isfinite(a[key]), np.isfinite(b[key])
            if not np.array_equal(fa, fb) or float(np.max(np.abs(a[key] - b[key])[fa], initial=0.0)) > TOL_RES:
                errs.append(f"段 {k} の {key}")
    return errs


def compare(label, res, gate_wall=None, thr_info=None):
    a, b = res["numpy"], res["numba"]
    nan_same, dx, da = lev_diff(a[0], b[0])
    da_ = dict(a[1]); db_ = dict(b[1])
    ra, rb = da_.pop("resid"), db_.pop("resid")
    # 最終残差は丸めの大きさ (~1e-15) なので、最大の値は許容差で比べ、最大になる段・対 (at_*) は丸めで入れ替わるので比べない
    dd = close_obj(da_, db_, TOL_X) + (close_obj({k: v for k, v in ra.items() if not k.startswith("at_")},
                                                 {k: v for k, v in rb.items() if not k.startswith("at_")}, TOL_RES)
                                       if ra is not None and rb is not None else close_obj(ra, rb, TOL_RES))
    te = rec_matches(a[2], b[2])
    ga, gb = moc_gate(dict(a[1]), thr_info, gate_wall), moc_gate(dict(b[1]), thr_info, gate_wall)
    check(f"{label}: 網 {a[0].shape} の NaN の位置が同じ・位置の差 {dx:.1e} ≤ {TOL_X:g}・θ/ν の差 {da:.1e} ≤ {TOL_A:g}、"
          f"集計の整数が同じで実数が許容差の中 (ビット同一 {bits_equal(a[0], b[0])})" + (f" (違い: {dd[:3]})" if dd else ""),
          nan_same and dx <= TOL_X and da <= TOL_A and not dd)
    check(f"{label}: 段・対ごとの状態・反復回数・源項の分岐が同じ、最終残差が許容差の中" + (f" ({te[:3]})" if te else ""), not te)
    check(f"{label}: moc_gate の判定が同じ (pass {ga['pass']}、理由 {len(ga['reasons'])} 件)", not close_obj(ga, gb, TOL_RES))
    return a[1]


# --- 1. backend の選び方 ----------------------------------------------------------------------------------------------
print("== 1. backend の選び方")
backend(None)
check("未設定は numba", MI.moc_kernel_backend() == "numba")
for v in ("numpy", "numba"):
    backend(v)
    check(f"{v!r} はそのまま", MI.moc_kernel_backend() == v)
for bad in ("", "NUMBA", "numpy ", "cuda"):
    backend(bad)
    try:
        MI.moc_kernel_backend()
        check(f"不正値 {bad!r} を拒否 (通ってしまった)", False)
    except ValueError as e:
        check(f"不正値 {bad!r} を拒否 ({str(e)[:60]})", MI.MOC_KERNEL_ENV in str(e))
    try:
        InverseMOC(gamma=G).fill_levels(RAD.build_init(10, 5))
        check(f"不正値 {bad!r} で fill_levels が止まる (通ってしまった)", False)
    except ValueError:
        check(f"不正値 {bad!r} で fill_levels が止まる", True)
backend(None)
rec = MI.moc_kernel_record()
check(f"記録: backend と numba・llvmlite の版 ({rec})", rec["backend"] == "numba" and rec["numba"] and rec["llvmlite"])
# numba を隠すと ImportError (別プロセス: numpy の版には落ちない、numpy を選べば動く)
_snip = (f"import sys; sys.modules['numba'] = None; sys.path.insert(0, {str(ROOT / 'design')!r}); "
         f"sys.path.insert(0, {str(Path(__file__).resolve().parent)!r}); import os\n"
         "import moc_axis_limit_radial as RAD\nfrom forge_design.geometry.moc_inverse import InverseMOC\n"
         "init = RAD.build_init(10, 5)\n"
         "try:\n    InverseMOC().fill_levels(init); print('NO_ERROR')\n"
         "except ImportError as e:\n    print('IMPORT_ERROR', 'FORGE_MOC_KERNEL=numpy' in str(e))\n"
         "os.environ['FORGE_MOC_KERNEL'] = 'numpy'; print('NUMPY_OK', InverseMOC().fill_levels(init).shape)\n")
_env = {k: v for k, v in os.environ.items() if k != MI.MOC_KERNEL_ENV}
_r = subprocess.run([PY, "-c", _snip], capture_output=True, text=True, env=_env)
check(f"numba が無い環境: 既定は ImportError (案内つき)、numpy を選べば動く ({_r.stdout.strip()!r} {_r.stderr[-200:]!r})",
      "IMPORT_ERROR True" in _r.stdout and "NUMPY_OK" in _r.stdout and "NO_ERROR" not in _r.stdout)
# 単位過程の差し替えを検出 (numba の版には効かない)
for mod, name in ((MK, "_sin_over_r_vec"), (MK, "pm_mach_vec"), (MI, "interior_vec")):
    orig = getattr(mod, name)
    setattr(mod, name, lambda *a, **k: orig(*a, **k))
    try:
        InverseMOC(gamma=G).fill_levels(RAD.build_init(10, 5))
        check(f"差し替え ({mod.__name__}.{name}) を検出 (通ってしまった)", False)
    except RuntimeError as e:
        check(f"差し替え ({mod.__name__}.{name}) を検出して止める", "FORGE_MOC_KERNEL=numpy" in str(e))
    finally:
        setattr(mod, name, orig)
# ガスの振り分け
from forge_design.gas import GasCPG  # noqa: E402
from forge_design.gas.isentrope import Isentrope  # noqa: E402
init_s = RAD.build_init(40, 9)
lev_f = InverseMOC(gamma=1.3).fill_levels(init_s)
for lab, gg in (("GasCPG(1.3)", GasCPG(1.3)), ("Isentrope(1.3) (CPG)", Isentrope(1.3)), ("np.float64(1.3)", np.float64(1.3))):
    backend("numpy")
    ref = InverseMOC(gamma=gg).fill_levels(init_s)
    backend(None)
    got = InverseMOC(gamma=gg).fill_levels(init_s)
    ok_, dx, da = lev_diff(ref, got)
    check(f"ガス {lab}: numba の版が float の γ と同一、numpy の版と許容差の中 (位置 {dx:.1e}・θ/ν {da:.1e})",
          bits_equal(got, lev_f) and ok_ and dx <= TOL_X and da <= TOL_A)
for lab, gg in (("bool", True), ("TP の Isentrope 以外の物 (dict)", {"gamma": 1.4})):
    try:
        InverseMOC(gamma=gg).fill_levels(init_s)
        check(f"扱わないガス {lab} を拒否 (通ってしまった)", False)
    except (ValueError, TypeError) as e:
        check(f"扱わないガス {lab} を拒否 ({str(e)[:60]})", True)

# --- 2. 部品 -----------------------------------------------------------------------------------------------------------
print("== 2. 部品 (ν→M の Newton と表の補間)")
rng = np.random.default_rng(20261015)
nus = np.concatenate([rng.uniform(0.0, 1.6, 20000), [0.0, -0.1, 1e-7, 2.2]])
E = np.zeros(1)
for gam in (1.4, 1.3, 1.15):
    ref = pm_mach_vec(nus, gam)
    for lab, M0 in (("冷えた初期値 (NaN)", np.full(nus.size, np.nan)), ("前の M ± 1 %", ref * (1.0 + rng.uniform(-0.01, 0.01, nus.size))),
                    ("前の M = 正解", ref)):
        got = np.array([NB._mach(v, m0, NB.GAS_CPG, gam, E, E) for v, m0 in zip(nus, M0)])
        d = float(np.max(np.abs(got - ref)))
        check(f"Newton (γ {gam}、{lab}、{nus.size} 点) が pm_mach_vec と |ΔM| {d:.1e} ≤ 1e-12", d <= 1e-12)
check("Newton: NaN は NaN", np.isnan(NB._mach(np.nan, 2.0, NB.GAS_CPG, 1.4, E, E)))
for n_t in (1, 2, 3, 400):
    tx = np.cumsum(rng.uniform(0.001, 0.01, n_t))
    ty = np.sort(rng.uniform(1.0, 6.0, n_t))
    xs = np.concatenate([rng.uniform(tx[0] - 0.1, tx[-1] + 0.1, 5000), tx])
    got = np.array([NB._mach(v, np.nan, NB.GAS_TABLE, 0.0, tx, ty) for v in xs])
    ref = np.interp(xs, tx, ty)
    d = float(np.max(np.abs(got - ref)))
    check(f"表の補間 (表 {n_t} 点、{xs.size} 点): np.interp と |ΔM| {d:.1e} ≤ 1e-12、表の節点・範囲外はそのままの値",
          d <= 1e-12 and np.array_equal(got[-n_t:], ty) and got[np.argmin(xs)] == ty[0] and got[np.argmax(xs)] == ty[-1])
check("表の補間: NaN は NaN", np.isnan(NB._mach(np.nan, np.nan, NB.GAS_TABLE, 0.0, tx, ty)))

# --- 3. 放射源流の網と semi-perfect の表 ------------------------------------------------------------------------------------
print("== 3. 網と集計 (4 つの組)")
for n_ax, n_st in ((60, 13), (300, 41)):
    init = RAD.build_init(n_ax, n_st)
    thr, info = axis_theta_r_init(init, n_ax, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                                  M_line_axis=RAD.src_M(RAD.X1, 0.0))
    for al, co in MODES:
        for how in ("fill_levels", "fill_arrays"):
            d = compare(f"放射源流 n_ax {n_ax} {al}+{co} ({how})", run_both(init, al, co, thr, how=how),
                        thr_info=info if al == "analytic" else None)
    # AXIS_LIMIT_FRAC を実行時に変える (numba の版は呼ぶたびに読む)
    frac0 = MK.AXIS_LIMIT_FRAC
    MK.AXIS_LIMIT_FRAC = 0.4
    try:
        d = compare(f"放射源流 n_ax {n_ax} AXIS_LIMIT_FRAC 0.4 (legacy+converge)", run_both(init, "legacy", "converge"))
        check(f"  AXIS_LIMIT_FRAC 0.4 で frac の分岐が発火している ({sum(d['branch'][k] for k in ("frac_A", "frac_B", "frac_PA", "frac_PB"))} 回)",
              sum(d["branch"][k] for k in ("frac_A", "frac_B", "frac_PA", "frac_PB")) > 0)
    finally:
        MK.AXIS_LIMIT_FRAC = frac0
# semi-perfect の表 (case/71 の燃焼ガス)
from forge_design.gas import GasSemiPerfect  # noqa: E402
from forge_design.gas.composition import ResolvedSpeciesDB  # noqa: E402
gas71 = GasSemiPerfect({"CO2": 0.1677, "H2O": 0.0858, "O2": 0.0220, "N2": 0.7245}, Tt=1250.0, db=ResolvedSpeciesDB.builtin())
init0 = RAD.build_init(120, 25)
init_g = [_Pt(p.x, p.r, p.th, float(gas71.nu(p.M)), gas71, float(gas71.mach_of_nu(float(gas71.nu(p.M))))) for p in init0]
thr_g = np.full(len(init_g), np.nan)
for i, p in enumerate(init_g):            # 軸上の既知点 (初期線の軸端を含む) に θ_r
    if p.r <= MK.AXIS_R_EPS:
        thr_g[i] = float(MK.axis_theta_r(p.M, RAD.src_dMdx_axis(p.x), gas71))
for al, co in MODES:
    compare(f"semi-perfect (case/71 の燃焼ガス) {al}+{co}", run_both(init_g, al, co, thr_g, gas=gas71))

# --- 4. 異常系 -------------------------------------------------------------------------------------------------------------
print("== 4. 異常系 (両方の backend)")


def P(x, r, th, M):
    return _Pt(float(x), float(r), float(th), float(pm_nu(M, G)), G)


abn = {
    # 平行な対だけの前線 (隣の対を置くと θ+μ = 90° の向きになり、丸めの差が大きく増幅される — 照合の入力に向かない)
    "平行な特性線 (θ ∓60°、M 2)": ([P(1.2, 0.6, np.pi / 3, 2.0), P(1.0, 0.5, -np.pi / 3, 2.0)], None, {"geom_parallel"}),
    "反復中の非有限 (座標 1e308)": ([P(3.0, 0.0, 0.0, 2.0), P(2.0, 0.1, 0.0, 2.0), P(1e308, 0.2, 0.1, 2.0), P(1.0, 0.3, 0.1, 2.0),
                                     P(0.5, 0.4, 0.1, 2.0)], None, {"geom_parallel"}),
    "軸より下 (軸の点が上流から下流の順)": ([P(1.0, 0.0, 0.0, 2.0), P(1.2, 0.0, 0.0, 2.1), P(1.4, 0.0, 0.0, 2.2), P(1.3, 0.3, 0.1, 2.0),
                                           P(1.2, 0.5, 0.15, 1.9)], np.array([0.3, 0.3, 0.3, np.nan, np.nan]), {"geom_below_axis", "missing"}),
    "反復の上限 (converge で 50 回)": ([P(1.2, 0.001, -0.6, 2.0), P(1.0, 0.002, -0.6, 2.0), P(0.9, 0.05, -0.2, 2.0),
                                       P(0.8, 0.1, 0.0, 2.0)], None, set()),
}
init_m = RAD.build_init(60, 13)
init_m[30].th = float("nan")
thr_m, info_m = axis_theta_r_init(RAD.build_init(60, 13), 60, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                                  M_line_axis=RAD.src_M(RAD.X1, 0.0))
abn["もともとの欠損 (初期前線の 1 点が NaN)"] = (init_m, thr_m, {"missing"})
seen = {k: 0 for k in ("missing", "geom_parallel", "geom_below_axis", "iter_nonfinite", "iter_maxiter")}
for lab, (init, thr, need) in abn.items():
    for al, co in MODES:
        if al == "analytic" and thr is None:
            continue
        d = compare(f"{lab} {al}+{co}", run_both(init, al, co, thr))
        for k in seen:
            seen[k] += d["pairs"][k]
        miss = [k for k in need if d["pairs"][k] == 0]
        check(f"  {lab} {al}+{co}: 期待の分類が出ている ({d['pairs']})", not miss)
# 反復の上限を小さくした充填 (放射源流)
init = RAD.build_init(60, 13)
d = compare("放射源流 converge max_corr 1 (上限到達)", run_both(init, "legacy", "converge", max_corr=1))
seen["iter_maxiter"] += d["pairs"]["iter_maxiter"]
check(f"異常系の分類が全部出た ({seen})", all(v > 0 for v in seen.values()))
# moc_gate (壁あり): 幾何的棄却の位置で判定が変わる入力
res = run_both(abn["平行な特性線 (θ ∓60°、M 2)"][0], "legacy", "converge")
for r_w, want in ((0.45, True), (0.8, False)):           # 棄却の対 (r 0.5・0.6) が壁の外 → 合格、壁の内側 → 不合格
    wall = np.column_stack([np.linspace(0.5, 3.0, 20), np.full(20, r_w)])
    ga, gb = moc_gate(dict(res["numpy"][1]), None, wall), moc_gate(dict(res["numba"][1]), None, wall)
    check(f"moc_gate (壁 r {r_w}) が同じで pass {want} ({ga['reasons']})", not same_obj(ga, gb) and ga["pass"] is want)
# θ_r の欠如 (段 1 の軸上の既知点、ValueError の文が同じ)
init_l, _, _ = abn["軸より下 (軸の点が上流から下流の順)"]
msgs = {}
for be in ("numpy", "numba"):
    backend(be)
    try:
        InverseMOC(gamma=G, axis_limit="analytic", corrector="converge").fill_levels(init_l, axis_thr=np.full(len(init_l), np.nan))
        msgs[be] = None
    except ValueError as e:
        msgs[be] = str(e)
backend(None)
check(f"θ_r の欠如: 両方とも同じ文の ValueError ({msgs['numba']})", msgs["numpy"] is not None and msgs["numpy"] == msgs["numba"])

# --- 5・6. 実データ (--real) ----------------------------------------------------------------------------------------------
if REAL:
    print("== 5. 実データ (V5)")
    import forge_design.evaluate.runner_axismach as RA  # noqa: E402
    from forge_design.evaluate.mknot_select import evaluate_mknot  # noqa: E402
    from forge_design.intake import wt_card as WC  # noqa: E402

    probs = []
    p70 = ROOT / "case/70.m5_air_wt_demo/_mknot_card_v4_3d/problem_c2b_auto_final.yaml"
    if p70.exists():
        probs.append(("M5 case/70 c2b", yaml.safe_load(p70.read_text()), p70.parent, 2.05))
    else:
        SKIPPED.append(f"5 M5 ({p70} が無い)")
    if CARD71.exists():
        card = WC.load_card(CARD71)
        probs.append(("M5.5 case/71", WC.card_to_problem(card, 0.1, M_knot=2.80), CARD71.parent, 2.80))
    else:
        SKIPPED.append(f"5 M5.5 ({CARD71} が無い)")
    p45 = C45 / "problem_d155_ns_prod.yaml"
    run45 = RUNS / "run_0062_euler_wallfit_fit_r1_ext6k"
    if p45.exists() and (run45 / "res_6000.h5").exists():
        raw45 = yaml.safe_load(p45.read_text())
        raw45["geometry"]["initial_line_run"] = str(run45)
        probs.append(("M6 case/45 prod", raw45, C45, float(raw45["geometry"]["M_knot"])))
    else:
        SKIPPED.append(f"5 M6 ({p45} か {run45} が無い)")

    TOL_X, TOL_A, TOL_TH = 1e-10, 1e-10, 1e-8
    table = []
    orig_inv = RA.inverse_design
    for lab, raw, base, MK_ in probs:
        for N in (600, 2400):
            L_tol = 0.1 if N == 600 else None
            calls = {"numpy": [], "numba": []}
            rows, times = {}, {}
            for be in ("numpy", "numba"):
                backend(be)

                def cap(*a, _be=be, **k):
                    r = orig_inv(*a, **k)
                    calls[_be].append({"levels": r["levels"], "wall": r["wall"], "wall_full": r["wall_full"], "moc": r.get("moc")})
                    if _be == "numba":          # numpy の版の同じ番の呼び出しとその場で比べて捨てる (メモリ)
                        i = len(calls["numba"]) - 1
                        c0, c1 = calls["numpy"][i], calls["numba"][i]
                        la, lb = c0.pop("levels"), c1.pop("levels")
                        fa, fb = np.isfinite(la), np.isfinite(lb)
                        both = fa & fb
                        dx = float(np.max(np.abs(la[..., :2] - lb[..., :2])[both[..., :2]], initial=0.0))
                        da = float(np.max(np.abs(la[..., 2:4] - lb[..., 2:4])[both[..., 2:4]], initial=0.0))
                        c1["cmp"] = {"shape": la.shape == lb.shape, "nan_same": bool(np.array_equal(fa, fb)), "dx": dx, "da": da,
                                     "bits": bits_equal(la, lb)}
                    return r

                RA.inverse_design = cap
                t = time.perf_counter()
                try:
                    rows[be] = evaluate_mknot(copy.deepcopy(raw), MK_, n_axis_inv=N, base_dir=base, L_total_tol=L_tol)
                finally:
                    RA.inverse_design = orig_inv
                times[be] = time.perf_counter() - t
            backend(None)
            ra, rb = rows["numpy"], rows["numba"]
            n_calls = (len(calls["numpy"]), len(calls["numba"]))
            ok_calls = n_calls[0] == n_calls[1] and n_calls[0] > 0
            worst = {"dx": 0.0, "da": 0.0, "dw": 0.0}
            bits = True
            n_moc_diff = 0
            for c0, c1 in zip(calls["numpy"], calls["numba"]):
                cm = c1["cmp"]
                ok_calls &= cm["shape"] and cm["nan_same"]
                bits &= cm["bits"]
                worst["dx"] = max(worst["dx"], cm["dx"])
                worst["da"] = max(worst["da"], cm["da"])
                for key in ("wall", "wall_full"):
                    wa, wb = np.asarray(c0[key]), np.asarray(c1[key])
                    if wa.shape != wb.shape:
                        ok_calls = False
                        continue
                    worst["dw"] = max(worst["dw"], float(np.max(np.abs(wa[:, :2] - wb[:, :2]), initial=0.0)))
                    bits &= bits_equal(wa, wb)
                ma, mb = c0["moc"] or {}, c1["moc"] or {}
                ok_calls &= all(ma.get(k) == mb.get(k) for k in ("pairs", "branch", "n_levels"))
                ok_calls &= (ma.get("gate") or {}).get("pass") == (mb.get("gate") or {}).get("pass")
                n_moc_diff += len(close_obj({k: v for k, v in ma.items() if k != "resid"}, {k: v for k, v in mb.items() if k != "resid"}, TOL_X))
            dth = (abs(ra["theta_max_deg"] - rb["theta_max_deg"]) if ("theta_max_deg" in ra and "theta_max_deg" in rb)
                   else (0.0 if ("theta_max_deg" not in ra and "theta_max_deg" not in rb) else float("inf")))
            judge = all(ra.get(k) == rb.get(k) for k in ("ok", "reasons", "error"))
            num_keys = [k for k in ra if k != "time_s" and isinstance(ra[k], float)]
            num_same = all((ra[k] == rb.get(k)) or (ra[k] != ra[k] and rb.get(k) != rb.get(k)) for k in num_keys)
            check(f"{lab} N {N}: 呼び出し {n_calls} の網の形・NaN の位置・集計が同じ、位置の差 {worst['dx']:.1e} ≤ {TOL_X:g} r_t・"
                  f"θ/ν の差 {worst['da']:.1e} ≤ {TOL_A:g} rad・壁の差 {worst['dw']:.1e} ≤ {TOL_X:g} r_t",
                  ok_calls and worst["dx"] <= TOL_X and worst["da"] <= TOL_A and worst["dw"] <= TOL_X)
            check(f"{lab} N {N}: θ_max の差 {dth:.1e} ≤ {TOL_TH:g}°、判定の欄 (ok {ra.get('ok')}・reasons・error) が同じ", dth <= TOL_TH and judge)
            print(f"info {lab} N {N}: ビット同一 {bits}、行の数値の欄 {len(num_keys)} 個が完全一致 {num_same}、"
                  f"集計の分類・分岐の数以外で許容差を超えた欄 {n_moc_diff} 個 (反復回数の分布など、記録)、"
                  f"評価 numpy {times['numpy']:.1f} s / numba {times['numba']:.1f} s ({times['numpy'] / times['numba']:.2f} 倍)")
            table.append({"problem": lab, "N": N, "calls": n_calls[0], "bit_identical": bool(bits), "row_numeric_identical": bool(num_same),
                          "max_dx_rt": worst["dx"], "max_dangle_rad": worst["da"], "max_dwall_rt": worst["dw"], "dtheta_max_deg": dth,
                          "moc_other_diffs": n_moc_diff, "ok": ra.get("ok"), "t_numpy_s": times["numpy"], "t_numba_s": times["numba"]})
            del calls
    print("== 6. 所要")
    init24 = RAD.build_init(2400, 41)
    thr24, _ = axis_theta_r_init(init24, 2400, lambda x: RAD.src_M(x, 0.0), G, target_dM=RAD.src_dMdx_axis,
                                 M_line_axis=RAD.src_M(RAD.X1, 0.0))
    t_fill = {}
    for be in ("numpy", "numba", "numba"):
        backend(be)
        t = time.perf_counter()
        InverseMOC(axis_limit="analytic", corrector="converge").fill_levels(init24, axis_thr=thr24)
        t_fill[be] = time.perf_counter() - t
    backend(None)
    print(f"info 放射源流 n_ax 2400 の fill_levels (analytic + converge): numpy {t_fill['numpy']:.2f} s / numba {t_fill['numba']:.2f} s "
          f"({t_fill['numpy'] / t_fill['numba']:.1f} 倍)")
    _snip6 = (f"import sys, time; sys.path.insert(0, {str(ROOT / 'design')!r}); sys.path.insert(0, {str(Path(__file__).resolve().parent)!r})\n"
              "import moc_axis_limit_radial as RAD\nfrom forge_design.geometry.moc_inverse import InverseMOC\n"
              "init = RAD.build_init(2400, 41)\nfor co in ('converge', 'converge'):\n"
              "    t = time.perf_counter(); InverseMOC(axis_limit='legacy', corrector=co).fill_levels(init); print(time.perf_counter() - t)\n")
    with tempfile.TemporaryDirectory() as td:
        env = {**{k: v for k, v in os.environ.items() if k != MI.MOC_KERNEL_ENV}, "NUMBA_CACHE_DIR": td}
        r = subprocess.run([PY, "-c", _snip6], capture_output=True, text=True, env=env)
        ts = [float(v) for v in r.stdout.split()] if r.returncode == 0 else []
    check(f"コンパイル込みの初回を測れた (放射源流 n_ax 2400 の converge: 初回 {ts[:1]} s、2 回目 {ts[1:2]} s; {r.stderr[-200:]!r})", len(ts) == 2)
    out_json = Path(tempfile.gettempdir()) / "run_moc_jit_tests_real.json"
    out_json.write_text(json.dumps({"table": table, "fill_2400_s": t_fill, "compile_first_s": ts[:1], "compiled_s": ts[1:2],
                                    "env": MI.moc_kernel_record(),
                                    "host": os.uname().nodename}, ensure_ascii=False, indent=1))
    print(f"info 表を {out_json} に書いた")
else:
    print("== 5・6 (実データ・所要) は --real のときだけ (飛ばした)")

print(f"\n{'FAIL ' + str(FAIL) if FAIL else 'ALL OK'}" + (f" (飛ばした: {SKIPPED})" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
