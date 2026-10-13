#!/usr/bin/env python3
"""edge の縁 y* の推定の持続 (`edge_band_positions(..., persist)`、既定 1 = 旧の抽出; K = 4 は検証 S2b で不採用 2026-10-13) の試験 — plan tooling-nozzle-cooled-wall-ns-settling
§4.2b・§5.1 #3・§6 S2b (c)。

  合成の梯子 (各段の R を直接与える: d_0..d_7 = 1、d_{k+8} = d_k (1 + R_k)。梯子の比 1.25^(1/8) なので 8 段先が 1.25 倍):
  (a) 一番内側の段だけ |R| < ε で次の段が外れる梯子 → persist 4 で NaN、persist 1 でその段 (段 0 なら y* = y_0、途中の段なら
      その段と 1 つ内側の段の間)。1 段だけの偶然の通過の後に 4 段続く平らな部分があれば persist 4 はそちらを拾う。梯子の終わりで
      persist に足りない平らな部分は拾わない。
  (b) 4 段以上続く平らな部分 → persist 1 と 4 で同じ y*・κ・y_b (ビット同一)。
  ほか: persist < 1 は拒否。既定は 1 (`edge_band_positions` と `deltastar_from_core_matched_euler` の両方; K = 4 は plan §6 S2b で FAIL のため不採用)。
  乱数の梯子 (各段の |R| を ε の内外から無作為に) で、persist K (2・4・8) が測定用の差し替え (case/70 band_persist_ab.py の
  文字列の差し替えを基準の版のソースに当てたもの) と全出力ビット同一。persist 1 は基準の版そのものとビット同一。
  実データ (case/70 run_0010 × run_0012、最新の場 res_60000; run が無ければ skip):
  (c) `edge_persist=1` の抽出が基準の版 (下の REF_COMMIT; 持続を入れる前の HEAD) の既定の抽出と `delta_r_next.csv`・`delta_r_equiv.csv`・
      `delta_r_edge.csv` がバイト同一 (`delta_r_equiv_diag.json` は settings.edge.persist の記録だけ違う)。
  (d) `edge_persist=4` の抽出が、基準の版に band_persist_ab.py の差し替え (K = 4) を当てた抽出と同じ 3 ファイルでバイト同一。

基準の版 REF_COMMIT は持続を入れる前の HEAD (32454f5f; deltastar.py は 1ca7160b から不変) に固定する
(HEAD のままだと commit の後に新しい実装どうしの比較になる)。`git show REF_COMMIT:design/forge_design/metrics/deltastar.py` を一時ファイルに
取り出して import する (相対 import の無いモジュール)。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_band_persist_tests.py [--quick]
所要: 合成 < 5 s、実データ (抽出 4 回) 約 3 分。--quick は実データを省く。
"""
from __future__ import annotations

import filecmp
import importlib.util
import inspect
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
sys.path.insert(0, str(DESIGN))

import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

REF_COMMIT = "32454f5f5660f9689dd60be54b5c791ae8f8af57"   # 持続を入れる前 (K = 1 の式) の HEAD
C70 = ROOT / "case" / "70.m5_air_wt_demo"
NS10 = C70 / "run_0010_ns_c2b_final"
EU12 = C70 / "run_0012_euler_c2b"
QUICK = "--quick" in sys.argv[1:]
G = 1.25 ** 0.125            # 抽出の梯子の比 (deltastar_from_core_matched_euler の edge_step 既定)
EPS = 0.003
FAIL = 0
FILES = ("delta_r_next.csv", "delta_r_equiv.csv", "delta_r_edge.csv")


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def skip(name, why):
    print(f"skip {name}  [{why}]", flush=True)


def load_ref_module(tmp: Path):
    """基準の版の deltastar.py を一時ファイルに取り出して別名で import する。"""
    src = subprocess.run(["git", "-C", str(ROOT), "show", f"{REF_COMMIT}:design/forge_design/metrics/deltastar.py"],
                         capture_output=True, text=True, check=True).stdout
    p = tmp / "deltastar_ref.py"
    p.write_text(src)
    spec = importlib.util.spec_from_file_location("deltastar_ref", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ab_patched(ref_mod, K: int):
    """case/70 band_persist_ab.py と同じ文字列の差し替え (交点の行を「K 段続けて |R| < ε の最初の段」に) を基準の版に当てた
    edge_band_positions。戻り: 差し替えた関数を持つモジュールの名前空間 (deltastar_from_core_matched_euler もその関数を呼ぶ)。"""
    src = inspect.getsource(ref_mod.edge_band_positions)
    old = "        j = int(np.argmax(Rabs < 0.0)) if np.any(Rabs < 0.0) else -1\n"
    assert src.count(old) == 1, "基準の版の edge_band_positions に交点の行が無い"
    new = ("        _ok = Rabs < 0.0\n"
           f"        _K = {K}\n"
           "        _c = [jj for jj in range(len(_ok) - _K + 1) if bool(np.all(_ok[jj:jj + _K]))]\n"
           "        j = _c[0] if _c else -1\n")
    g = dict(ref_mod.__dict__)
    exec(compile(src.replace(old, new), "patched", "exec"), g)
    return g["edge_band_positions"]


# ================================================================================================ 合成の梯子
def ladder_from_R(R, y0=0.02):
    """段 k の R が R[k] になる梯子 (y_k = y0 G^k、d_0..d_7 = 1、d_{k+8} = d_k (1 + R_k))。"""
    R = np.asarray(R, dtype=float)
    n = len(R) + 8
    y = y0 * G ** np.arange(n)
    d = np.ones(n)
    for k, r in enumerate(R):
        d[k + 8] = d[k] * (1.0 + r)
    return y, d


def ladder_R(y, d):
    """edge_band_positions と同じ式で各段の R を出す (合成の梯子の前提の確認用)。"""
    ly = np.log(y); lq = ly + np.log(1.25); m = lq <= ly[-1] + 1e-12
    return np.interp(lq[m], ly, d) / d[m] - 1.0


GOOD = [0.2, 0.15, 0.1, 0.05, 0.02, 0.002, 0.001, 0.0008, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001, 0.0001, 0.0001, 0.0001]


def stations(test_R, n=60, every=2):
    """n 断面 (x −10〜50)。every 番目ごとに試験の梯子、ほかは平らな部分が長い梯子 (P-spline に支持点を与える)。"""
    x = np.linspace(-10.0, 50.0, n)
    yl, dl, idx = [], [], []
    for i in range(n):
        if i % every == 1:
            y, d = ladder_from_R(test_R); idx.append(i)
        else:
            y, d = ladder_from_R(GOOD)
        yl.append(y); dl.append(d)
    return x, yl, dl, np.full(n, 0.01), np.ones(n), idx


def run_edge(fn, x, yl, dl, din, rw, **kw):
    return fn(x, yl, dl, din, rw, eps=EPS, **kw)


def test_synthetic(ref_mod):
    print("\n--- 合成の梯子 ---")
    E = DS.edge_band_positions
    sig = inspect.signature(E).parameters["persist"].default
    sig2 = inspect.signature(DS.deltastar_from_core_matched_euler).parameters["edge_persist"].default
    check("既定: edge_band_positions(persist=1)・deltastar_from_core_matched_euler(edge_persist=1) (K = 4 は不採用)", sig == 1 and sig2 == 1, f"{sig}, {sig2}")
    try:
        E(*stations(GOOD)[:5], persist=0)
        check("persist 0 を拒否", False)
    except ValueError:
        check("persist 0 を拒否", True)

    # (a) 一番内側の段 (段 0) だけ |R| < ε、次の段が外れる
    Ra = [0.002] + [0.1] * 15
    y, d = ladder_from_R(Ra)
    Rchk = ladder_R(y, d)
    pre = abs(Rchk[0]) < EPS and np.all(np.abs(Rchk[1:]) > EPS)
    x, yl, dl, din, rw, idx = stations(Ra)
    p1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    p4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    y0 = float(np.exp(np.log(yl[idx[0]][0])))
    check("(a) 段 0 だけ |R| < ε・次の段は外れる: persist 1 は y* = 段 0、persist 4 は NaN",
          pre and np.all(p1["y_star"][idx] == y0) and np.all(np.isnan(p4["y_star"][idx])),
          f"R[0:3] {Rchk[:3]}, persist 1 の y* {float(p1['y_star'][idx[0]])!r} (段 0 {y0!r}), persist 4 {float(p4['y_star'][idx[0]])!r}")
    good_idx = [i for i in range(len(x)) if i not in idx]
    check("(a) 平らな梯子の断面は persist 1 と 4 で同じ y*", np.array_equal(p1["y_star"][good_idx], p4["y_star"][good_idx]))

    # (a') 途中の段 (段 3) だけ |R| < ε → persist 1 は段 2 と 3 の間、persist 4 は NaN
    Rb = [0.1, 0.08, 0.05, 0.002, 0.05, 0.1] + [0.1] * 10
    x, yl, dl, din, rw, idx = stations(Rb)
    p1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    p4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    yy = yl[idx[0]]
    check("(a') 段 3 だけ |R| < ε: persist 1 の y* は段 2〜3 の間、persist 4 は NaN",
          np.all((p1["y_star"][idx] > yy[2]) & (p1["y_star"][idx] < yy[3])) and np.all(np.isnan(p4["y_star"][idx])),
          f"y* {p1['y_star'][idx[0]]:.6e} (段 2 {yy[2]:.6e}・段 3 {yy[3]:.6e})")

    # (a'') 段 1 の偶然の通過の後、段 4〜7 の 4 段が |R| < ε → persist 4 は段 3〜4 の間、persist 1 は段 0〜1 の間
    Rc = [0.1, 0.002, 0.1, 0.05, 0.001, 0.001, 0.0015, 0.002, 0.05, 0.1] + [0.1] * 6
    x, yl, dl, din, rw, idx = stations(Rc)
    p1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    p4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    yy = yl[idx[0]]
    check("(a'') 1 段の通過の後に 4 段の平らな部分: persist 4 は段 3〜4 の間、persist 1 は段 0〜1 の間",
          np.all((p4["y_star"][idx] > yy[3]) & (p4["y_star"][idx] < yy[4]))
          and np.all((p1["y_star"][idx] > yy[0]) & (p1["y_star"][idx] < yy[1])),
          f"persist 4 {p4['y_star'][idx[0]]:.6e}, persist 1 {p1['y_star'][idx[0]]:.6e}")

    # (a''') 梯子の終わりの 3 段だけ |R| < ε → persist 4 は NaN、persist 3 は拾う
    Rd = [0.1] * 13 + [0.001, 0.001, 0.001]
    x, yl, dl, din, rw, idx = stations(Rd)
    p3 = run_edge(E, x, yl, dl, din, rw, persist=3)
    p4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    yy = yl[idx[0]]
    check("(a''') 梯子の終わりで persist に足りない平らな部分は拾わない (persist 4 は NaN、persist 3 は段 12〜13 の間)",
          np.all(np.isnan(p4["y_star"][idx])) and np.all((p3["y_star"][idx] > yy[12]) & (p3["y_star"][idx] < yy[13])))

    # (b) 4 段以上続く平らな部分 → persist 1 と 4 でビット同一 (全断面がこの形)
    x, yl, dl, din, rw, idx = stations(GOOD, every=1)
    p1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    p4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    same = all(np.array_equal(p1[k], p4[k], equal_nan=True) for k in ("y_star", "kappa", "kappa_med", "kappa_s", "y_b"))
    check("(b) 段 5 から平らな梯子: persist 1 と 4 で y*・κ・κ_s・y_b がビット同一", same and np.all(np.isfinite(p1["y_star"])),
          f"y* {p1['y_star'][0]:.6e}")
    # 交点の段の 1 つ内側も平ら (段 4 が |R| < ε でも 1 段手前) の場合: 段 4 から平らなら同じく同一
    Rf = [0.2, 0.1, 0.05, 0.01, 0.0029, 0.002, 0.001, 0.001, 0.0005] + [0.0001] * 7
    x, yl, dl, din, rw, idx = stations(Rf, every=1)
    q1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    q4 = run_edge(E, x, yl, dl, din, rw, persist=4)
    check("(b) 段 4 から平らな梯子 (交点の段が |R| ≈ ε のすぐ内側): persist 1 と 4 でビット同一",
          all(np.array_equal(q1[k], q4[k], equal_nan=True) for k in ("y_star", "kappa_s", "y_b")))

    # 乱数の梯子: persist K = 測定用の差し替え (K)、persist 1 = 基準の版
    rng = np.random.default_rng(20261013)
    n_st = 60
    x = np.linspace(-12.0, 60.0, n_st)
    yl, dl = [], []
    for i in range(n_st):
        if i % 3 == 0:
            y, d = ladder_from_R(GOOD)
        else:
            L = int(rng.integers(6, 40))
            R = np.where(rng.random(L) < 0.45, rng.uniform(-0.0029, 0.0029, L), rng.choice([-1.0, 1.0], L) * rng.uniform(0.0031, 0.05, L))
            y, d = ladder_from_R(R, y0=float(rng.uniform(0.01, 0.03)))
        yl.append(y); dl.append(d)
    din = rng.uniform(0.0, 0.02, n_st); rw = np.linspace(2.0, 1.0, n_st)
    keys = ("y_star", "kappa", "kappa_med", "kappa_s", "y_b", "y_ref")
    r1 = run_edge(ref_mod.edge_band_positions, x, yl, dl, din, rw)
    n1 = run_edge(E, x, yl, dl, din, rw, persist=1)
    check("乱数の梯子: persist 1 = 基準の版 (全出力ビット同一)", all(np.array_equal(r1[k], n1[k], equal_nan=True) for k in keys),
          f"y* の NaN {int(np.isnan(n1['y_star']).sum())}/{n_st}")
    for K in (2, 4, 8):
        rk = run_edge(ab_patched(ref_mod, K), x, yl, dl, din, rw)
        nk = run_edge(E, x, yl, dl, din, rw, persist=K)
        check(f"乱数の梯子: persist {K} = band_persist_ab.py の差し替え (K = {K}) (全出力ビット同一)",
              all(np.array_equal(rk[k], nk[k], equal_nan=True) for k in keys),
              f"y* の NaN {int(np.isnan(nk['y_star']).sum())}/{n_st}、persist 1 と違う断面 "
              f"{int(np.sum(~((n1['y_star'] == nk['y_star']) | (np.isnan(n1['y_star']) & np.isnan(nk['y_star'])))))}")


# ================================================================================================ 実データ
def diag_without_persist(p: Path) -> dict:
    d = json.loads(p.read_text())
    e = d.get("settings", {}).get("edge")
    if isinstance(e, dict):
        e.pop("persist", None)
    return d


def compare_dirs(a: Path, b: Path) -> dict:
    return {f: filecmp.cmp(a / f, b / f, shallow=False) for f in FILES}


def extract_with(mod_fn, out, **kw):
    """extract_and_merge が呼ぶ抽出器を mod_fn に差し替えて抽出する (後処理の 5 次 P-spline 平滑化は共通)。"""
    orig = DL.deltastar_from_core_matched_euler
    DL.deltastar_from_core_matched_euler = mod_fn
    try:
        DL.extract_and_merge(NS10, EU12, band_select="edge", out_dir=out, **kw)
    finally:
        DL.deltastar_from_core_matched_euler = orig


def test_real(ref_mod, tmp: Path):
    print("\n--- 実データ: case/70 run_0010 × run_0012 (最新の場) ---")
    if not (NS10.exists() and EU12.exists()):
        skip("(c)(d) 実データ", "run_0010・run_0012 が無い")
        return
    res = sorted(NS10.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1].name
    # (c) persist 1 = 基準の版の既定
    DL.extract_and_merge(NS10, EU12, band_select="edge", edge_persist=1, out_dir=tmp / "new_p1")
    extract_with(ref_mod.deltastar_from_core_matched_euler, tmp / "ref")
    same = compare_dirs(tmp / "new_p1", tmp / "ref")
    dn, dr = diag_without_persist(tmp / "new_p1" / "delta_r_equiv_diag.json"), json.loads((tmp / "ref" / "delta_r_equiv_diag.json").read_text())
    p1rec = json.loads((tmp / "new_p1" / "delta_r_equiv_diag.json").read_text())["settings"]["edge"].get("persist")
    check(f"(c) edge_persist=1 の抽出 = 基準の版 {REF_COMMIT[:8]} の既定の抽出 ({res}; 3 ファイルがバイト同一)", all(same.values()), f"{same}")
    check("(c) delta_r_equiv_diag.json は settings.edge.persist (= 1) の記録だけ違う", dn == dr and p1rec == 1, f"persist の記録 {p1rec!r}")
    # (d) 既定 (persist 4) = band_persist_ab.py の差し替え (K = 4)
    DL.extract_and_merge(NS10, EU12, band_select="edge", edge_persist=4, out_dir=tmp / "new_default")
    g = dict(ref_mod.__dict__)
    exec(compile(inspect.getsource(ref_mod.deltastar_from_core_matched_euler), "ref_patched", "exec"), g)
    g["edge_band_positions"] = ab_patched(ref_mod, 4)
    extract_with(g["deltastar_from_core_matched_euler"], tmp / "ab_k4")
    same4 = compare_dirs(tmp / "new_default", tmp / "ab_k4")
    p4rec = json.loads((tmp / "new_default" / "delta_r_equiv_diag.json").read_text())["settings"]["edge"].get("persist")
    ys_new = np.loadtxt(tmp / "new_default" / "delta_r_edge.csv", delimiter=",", skiprows=1)[:, 2]
    ys_p1 = np.loadtxt(tmp / "new_p1" / "delta_r_edge.csv", delimiter=",", skiprows=1)[:, 2]
    check(f"(d) edge_persist=4 の抽出 = band_persist_ab.py の差し替え (K = 4) の抽出 ({res}; 3 ファイルがバイト同一)", all(same4.values()),
          f"{same4}、y* の NaN: persist 4 {int(np.isnan(ys_new).sum())} / persist 1 {int(np.isnan(ys_p1).sum())} (全 {len(ys_new)})")
    check("(d) edge_persist=4 の抽出の settings.edge.persist = 4", p4rec == 4, f"{p4rec!r}")
    d_new = DL.read_delta_r_next(tmp / "new_default" / "delta_r_next.csv")["delta_E"]
    d_p1 = DL.read_delta_r_next(tmp / "new_p1" / "delta_r_next.csv")["delta_E"]
    print(f"     参考: persist 4 と 1 の δ_E の差の最大 {float(np.max(np.abs(d_new - d_p1))):.3e} r_t (全域)", flush=True)


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ref_mod = load_ref_module(tmp)
        test_synthetic(ref_mod)
        if QUICK:
            skip("(c)(d) 実データ", "--quick")
        else:
            test_real(ref_mod, tmp)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
