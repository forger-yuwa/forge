#!/usr/bin/env python3
"""③ベル (M) の δ_r 抽出と contur_v2 の接続 (plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)・(7b-3) 前半、codex diagnose
2026-10-12 M1〜M5・m1) の試験 (pytest 非依存・素の assert)。

  (1) ③の抽出が一定 x の断面で回る: 列の中で x が最大 0.31 r_t ずれる合成の③の run (近壁が壁法線の層に相当) で、
      `extract_bell` (一定 x) の δ_r が真値 (解析場の一定 x の質量欠損) に 2 % 以内、列のままの抽出は 10 % 以上ずれる
  (2) x_E 欠け: ③の prepare_info (x_E 無し・bell_mode・wall_design.L) で流量の帳簿が止まらず L を使い出所を書く。
      x_E も L も無い run は例外にせず `skipped` と理由を返し、抽出は最後まで回る。一定 x の帳簿は列のままより正しい
  (3) 帯位置の感度: k_band の感度系列は主系列とビット一致 (codex M2 の前提)、帯位置の倍率 1.0 は主系列とビット一致、
      境界層に入る倍率は主系列と違う値を出す。受理の判定 (外の倍率だけなら受理、内の倍率を入れると不受理)
  (4) ①の既定の抽出がビット同一: 基準コミット 06089ab5 の `metrics/deltastar.py` と、合成の①の run (x_E あり) で
      既定・edge+一定 x・band_smooth_x・core 法・流量の帳簿の全出力が一致 (手元に case/45 の run があれば実データでも)
  (5) δ_r の時系列: 3 枚のスナップショットの合成 run で CSV (step 列・代表 x・L2) を書き、値が step で変わる。
      `check_quasisteady.py --series-csv` が列を読める
  (6) contur_v2 の入口: edge_M_full は表の範囲外・contur_v1 を拒む、表の M をそのまま使う。環状の定義 (M5) の逆変換が厳密。
      手元に③の run があれば: 既定の経路 (edge_M_full 無し) が基準コミットの deltastar_integral.py とビット同一、
      ③の抽出が codex の再計算 (x 1.0: 列 0.0031796 → 一定 x 0.0041637) を再現
usage: design/.venv-opt/bin/python design/tests/run_bell_deltastar_tests.py
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.metrics import bell_deltastar as bd  # noqa: E402
from forge_design.metrics import deltastar as ds  # noqa: E402

BASE_COMMIT = "06089ab5"          # 変更前 (③の拡張の前) のコミット
MAIN_BASE = "81e92491"            # (4) の基準: 取り込んだ main (2026-10-15)。③の追加が main の①の既定の抽出を変えないことを見る
                                  # (06089ab5 と比べると main 側の正当な変更 [not_extracted の鍵など] まで差に数えてしまう)
fails = []
_BASE_TMP = []


def check(name, ok, detail=""):
    print(f"[{'OK ' if ok else 'NG '}] {name}  {detail}")
    if not ok:
        fails.append(name)


def load_base(relpath, modname, commit=None):
    """基準コミット (既定 BASE_COMMIT) のモジュールを読み込む (パッケージ内の名前で — 相対 import のため)。取れなければ None。"""
    commit = commit or BASE_COMMIT
    try:
        src = subprocess.run(["git", "show", f"{commit}:{relpath}"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    except Exception as e:                                           # noqa: BLE001
        print(f"基準コミット {commit} の {relpath} が取れない: {e}")
        return None
    d = Path(tempfile.mkdtemp(prefix="bell_dstar_base_")); _BASE_TMP.append(d)
    tmp = d / Path(relpath).name
    tmp.write_text(src)
    spec = importlib.util.spec_from_file_location(modname, tmp)
    mod = importlib.util.module_from_spec(spec); sys.modules[modname] = mod; spec.loader.exec_module(mod)
    return mod


def same(a, b, path=""):
    """dict / 配列 / 値のビット一致 (鍵の順も)。違う箇所を返す (一致なら None)。"""
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)) or list(a) != list(b):
            return f"{path}: 鍵が違う {set(a) ^ set(b) if isinstance(a, dict) and isinstance(b, dict) else type(b)}"
        for k in a:
            r = same(a[k], b[k], f"{path}/{k}")
            if r:
                return r
        return None
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a_, b_ = np.asarray(a), np.asarray(b)
        ok = a_.dtype == b_.dtype and a_.shape == b_.shape and (
            np.array_equal(a_, b_, equal_nan=True) if a_.dtype.kind == "f" else np.array_equal(a_, b_))
        return None if ok else f"{path}: 配列が違う"
    if isinstance(a, float) and isinstance(b, float) and np.isnan(a) and np.isnan(b):
        return None
    return None if a == b else f"{path}: {a!r} != {b!r}"


# ------------------------------------------------------------------ 合成の run
NI, NJ, S, P = 81, 81, 0.1, 2.5
XI = np.linspace(0.0, 10.0, NI)
SIN = float(np.sin(np.arctan(0.25)))


def rw(x):
    return 1.0 + 0.25 * x


def Q(x):
    return 1.0 + 0.3 * x                              # 核の流量密度が x で変わる (列の x のずれが比に出る)


def dl(x, k=1.0):
    return k * (0.01 + 0.004 * x)                     # 境界層の厚さ (帯の下限 0.02 r_w より薄い)


def q_ns(x, r, k=1.0):
    s = np.clip((rw(x) - r) / dl(x, k), 0.0, None)
    return Q(x) * np.where(s < 1.0, s ** (1.0 / 7.0), 1.0)


def true_dr(x, k=1.0):
    """一定 x の断面の質量欠損 → 環状の等価厚さ (q_ref = Q(x) 一様)。"""
    r = np.linspace(rw(x) - dl(x, k), rw(x), 200001)
    a = np.trapezoid((Q(x) - q_ns(x, r, k)) * r, r) / Q(x)
    return rw(x) - np.sqrt(rw(x) ** 2 - 2.0 * a)


def write_run(rd: Path, slant: float, info_extra: dict, field, steps=(100,)):
    """node 構造格子 (index = i*nj + j、j=0 軸 → j=nj-1 壁) の run。slant > 0 で列の中の x を壁から内側へずらす。"""
    rd.mkdir(parents=True, exist_ok=True)
    t = 1.0 - np.arange(NJ) / (NJ - 1)
    X = np.empty((NI, NJ)); R = np.empty((NI, NJ))
    for i, x0 in enumerate(XI):
        y = rw(x0) * t ** P
        R[i] = rw(x0) - y
        X[i] = x0 + slant * y * SIN * (1.0 - t)
    co = np.zeros((NI * NJ, 3)); co[:, 0] = X.ravel() * S; co[:, 1] = R.ravel() * S
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = co.ravel()
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": S, **info_extra}))
    for st in steps:
        with h5py.File(rd / f"res_{st}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(NI * NJ); f["/VALUE/Ux"] = np.asarray(field(X, R, st), dtype=float).ravel()
    return rd


TMP = Path(tempfile.mkdtemp(prefix="bell_dstar_test_"))
BELL = {"bell_mode": "moc_rao", "wall_design": {"L": 10.0}}
ns3 = write_run(TMP / "bell_ns", 3.0, BELL, lambda X, R, st: q_ns(X, R))
eu3 = write_run(TMP / "bell_eu", 0.0, BELL, lambda X, R, st: Q(X))

# (1) 一定 x の断面
col_dx = float(np.max(np.ptp(ds._load_structured(ns3)["x"], axis=1)))
o_new = bd.extract_bell(ns3, eu3, yb_sens_factors=(1.5, 2.0))
o_old = ds.deltastar_from_core_matched_euler(ns3, eu3)
xs = [1.0, 2.0, 4.0, 6.0, 8.0, 9.5]
idx = [int(np.argmin(abs(o_new["x"] - v))) for v in xs]
err_new = [o_new["delta_r_raw"][i] / true_dr(o_new["x"][i]) - 1.0 for i in idx]
err_old = [o_old["delta_r_raw"][i] / true_dr(o_old["x"][i]) - 1.0 for i in idx]
check("(1) 合成の③の run の列は x が大きくずれる", col_dx > 0.25, f"列内の x の幅 最大 {col_dx:.3f} r_t")
check("(1) 一定 x の抽出が真値に 2 % 以内", max(abs(e) for e in err_new) < 0.02, "相対差 " + " ".join(f"{e:+.4f}" for e in err_new))
check("(1) 列のままの抽出は 10 % 以上ずれる (直した意味がある)", min(abs(e) for e in err_old) > 0.10,
      "相対差 " + " ".join(f"{e:+.4f}" for e in err_old))
check("(1) ③の呼び出しは一定 x が既定", o_new["settings"].get("remap_constant_x") is True and o_old["settings"].get("remap_constant_x") is None)

# (2) x_E 欠け
mf = o_new["massflow"]
check("(2) ③ (x_E 無し) の流量の帳簿は wall_design.L を使い出所を書く",
      not mf.get("skipped") and mf.get("x_E_ns") == 10.0 and "wall_design.L" in str(mf.get("x_E_ns_source")),
      f"x_E {mf.get('x_E_ns')} ({mf.get('x_E_ns_source')})")
ns_bare = write_run(TMP / "bare_ns", 3.0, {}, lambda X, R, st: q_ns(X, R))
eu_bare = write_run(TMP / "bare_eu", 0.0, {}, lambda X, R, st: Q(X))
try:
    mf_b = ds.massflow_ratio(ns_bare, eu_bare)
    o_b = ds.deltastar_from_core_matched_euler(ns_bare, eu_bare)
    check("(2) x_E も L も無い run は止めずに理由を返す", mf_b.get("skipped") is True and "x_E" in mf_b.get("reason", "")
          and o_b["massflow"].get("skipped") is True and len(o_b["x"]) == NI, mf_b.get("reason", "")[:90])
except Exception as e:                                               # noqa: BLE001
    check("(2) x_E も L も無い run は止めずに理由を返す", False, f"例外 {e}")
# 一定 x の帳簿: 列のまま (−ρU_r dx の項が落ちる) より流量の比が真値に近い (真値 = 一定 x の断面の厳密な比)
xs_all = XI[(XI > -2.5) & (XI < 10.0)]
ratio_true = []
for x in xs_all:
    r = np.linspace(0.0, rw(x), 400001)
    ratio_true.append(np.trapezoid(q_ns(x, r) * r, r) / np.trapezoid(Q(x) * r, r))
m_true = float(np.median(ratio_true))
rel_mx = abs(mf["mdot_ratio"] / m_true - 1.0); rel_col = abs(o_new["massflow_column"]["mdot_ratio"] / m_true - 1.0)
check("(2) 一定 x の流量の帳簿は列のままより正しい", rel_mx < 0.2 * rel_col,
      f"真値に対し 一定 x {rel_mx:.2e} / 列のまま {rel_col:.2e} (中央値の比どうし)")

# (3) 帯位置の感度
o3 = bd.extract_bell(ns3, eu3, yb_sens_factors=(0.3, 1.0, 2.0))
check("(3) k_band の感度系列は主系列とビット一致 (同じ壁、delta_in = 0)", np.array_equal(o3["delta_r_sens"][:, 0], o3["delta_r_raw"])
      and np.array_equal(o3["delta_r_sens"][:, 1], o3["delta_r_raw"]))
check("(3) 帯位置の倍率 1.0 は主系列とビット一致", np.array_equal(o3["delta_r_sens_yb"][:, 1], o3["delta_r_raw"]))
d03 = o3["delta_r_sens_yb"][:, 0] / o3["delta_r_raw"] - 1.0
check("(3) 境界層に入る帯位置 (0.3 倍) は主系列と違う値を出す", np.all(np.isfinite(d03)) and np.min(np.abs(d03)) > 1e-3,
      f"相対差 {np.min(d03):+.3f}〜{np.max(d03):+.3f}")
u_ok, u_ng = o_new["uncertainty"], o3["uncertainty"]
check("(3) 受理の判定: 外の倍率 (1.5, 2.0) だけなら受理、0.3 を入れると不受理と理由 yb_sens",
      u_ok["accepted"] and not u_ng["accepted"] and u_ng["n_over"] > 0 and any("yb_sens" in r for r in o3["reason"]),
      f"外のみ u_max {u_ok['u_max']:.4f} / 0.3 込み u_max {u_ng['u_max']:.3f}")
check("(3) 既定の抽出器は帯位置の感度の鍵を足さない", "delta_r_sens_yb" not in o_old and "yb_sens_factors" not in o_old["settings"]
      and not any(k.startswith("yb_sens_") or k == "band_ratio_a" for k in o_old))

# (3b) 帯の診断と A/B の判定 (plan §5.1 #7「次の A/B (事前登録)」)
same_diag = all(np.array_equal(o3["yb_sens_" + k][:, 1], o3[k]) for k in ("band_slope", "band_rms", "band_deficit_share", "band_ratio_a", "alpha"))
check("(3b) 感度帯の診断: 倍率 1.0 の列は主系列の帯の診断とビット一致", same_diag and o3["yb_sens_band_y_b"].shape == (len(o3["x"]), 3)
      and np.array_equal(o3["yb_sens_band_y_b"][:, 1], o3["band_y_b"]))
acc_ok = bd.band_acceptance(o_new)
acc_ng = bd.band_acceptance(o3)
check("(3b) A/B の判定: 外の倍率だけ (1.5, 2.0) は全条件を満たし、0.3 倍は感度と感度帯の停止条件で落ちる",
      acc_ok["all_pass"] and not acc_ng["all_pass"] and acc_ng["summary"]["sens"]["n_fail"] > 0
      and acc_ng["summary"]["sens_stop_parts"]["f0.3"]["slope"]["n_fail"] > 0 and acc_ng["summary"]["main_stop"]["n_fail"] == 0
      and "edge_found" not in acc_ok["criteria"],
      f"0.3 倍の帯の傾き不成立 {acc_ng['summary']['sens_stop_parts']['f0.3']['slope']['n_fail']} 断面")
o_e = bd.extract_bell(ns3, eu3, yb_sens_factors=(0.75, 1.5, 2.0), band_select="edge", return_ladders=True)
acc_e = bd.band_acceptance(o_e)
E = acc_e["edge"]
check("(3b) edge の判定: 局所縁・帯が縁の外・平滑化の前後の縁の位置・端の保持を断面ごとに返す",
      {"edge_found", "band_outside_edge"} <= set(acc_e["criteria"]) and E is not None
      and np.array_equal(E["y_b"], o_e["band_y_b"]) and np.allclose(E["y_b"], E["c"] * E["y_edge_smooth"], rtol=1e-12)
      and np.array_equal(E["held"], (o_e["x"] < 3.0) | (o_e["x"] > o_e["x"].max() - 3.0)) and np.all(np.isfinite(E["y_star"])),
      f"all_pass {acc_e['all_pass']}、不成立 " + ", ".join(f"{k} {acc_e['summary'][k]['n_fail']}" for k in acc_e["criteria"]))

# (4) ①の既定がビット同一 (基準 = 取り込んだ main の deltastar.py)
base_ds = load_base("design/forge_design/metrics/deltastar.py", "forge_design.metrics._deltastar_base", commit=MAIN_BASE)
WT = {"x_E": 8.0}
ns1 = write_run(TMP / "wt_ns", 3.0, WT, lambda X, R, st: q_ns(X, R))
eu1 = write_run(TMP / "wt_eu", 0.0, WT, lambda X, R, st: Q(X))
if base_ds is None:
    check("(4) 基準コミットのモジュールが読める", False, MAIN_BASE)
else:
    cases = [{}, {"band_select": "edge", "remap_constant_x": True}, {"band_smooth_x": 4.0}, {"method": "core"}]
    for kw in cases:
        r_ = same(base_ds.deltastar_from_core_matched_euler(ns1, eu1, **kw), ds.deltastar_from_core_matched_euler(ns1, eu1, **kw))
        check(f"(4) ①の抽出がビット同一 (合成、{kw or '既定'})", r_ is None, r_ or "")
    r_ = same(base_ds.massflow_ratio(ns1, eu1), ds.massflow_ratio(ns1, eu1))
    check("(4) ①の流量の帳簿がビット同一 (合成)", r_ is None, r_ or "")
    real = [(c / "run_0038_ns_final_rt77p02", c / "run_0037_euler_rt77p02")
            for c in (ROOT / "case/45.isobutane_m6_d155", Path("/home/sano/work/forge/case/45.isobutane_m6_d155"))]
    real = [(a, b) for a, b in real if (a / "nozzle.h5").is_file() and (b / "nozzle.h5").is_file()]
    if real:
        a, b = real[0]
        for kw in ({}, {"band_smooth_x": 4.0}):
            r_ = same(base_ds.deltastar_from_core_matched_euler(a, b, **kw), ds.deltastar_from_core_matched_euler(a, b, **kw))
            check(f"(4) ①の抽出がビット同一 (実データ {a.name}、{kw or '既定'})", r_ is None, r_ or "")
        r_ = same(base_ds.massflow_ratio(a, b), ds.massflow_ratio(a, b))
        check("(4) ①の流量の帳簿がビット同一 (実データ)", r_ is None, r_ or "")
    else:
        print("[-- ] (4) 実データ (case/45 run_0038/0037) が無い — 合成だけ")

# (5) 時系列
ns_t = write_run(TMP / "series_ns", 3.0, BELL, lambda X, R, st: q_ns(X, R, k=1.0 + 0.002 * (st - 200) / 100.0),
                 steps=(100, 200, 300))
(ns_t / "res_wall_3_300.h5").write_bytes(b"")                       # res_<n>.h5 以外は拾わない
csv = TMP / "series.csv"
sr = bd.deltastar_series(ns_t, eu3, csv, out_npz=TMP / "series.npz")
rows = csv.read_text().splitlines()
hdr = rows[0].split(",")
vals = np.array([[float(v) for v in r.split(",")[2:-1]] for r in rows[1:]])
check("(5) 時系列の CSV: step 列・代表 x・出口・L2、3 行", hdr[:2] == ["step", "res"] and "dr_x1" in hdr and "dr_exit" in hdr
      and "dr_l2" in hdr and len(rows) == 4 and [r.split(",")[0] for r in rows[1:]] == ["100", "200", "300"], ",".join(hdr))
check("(5) δ_r が step で変わる (境界層を厚くした)", np.all(np.diff(vals[:, hdr.index("dr_l2") - 2]) > 0))
cq = subprocess.run([sys.executable, str(ROOT / "solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(csv),
                     "--series-cols", ",".join(sr["cols"]), "--min-snaps", "3"], capture_output=True, text=True)
check("(5) check_quasisteady --series-csv が列を読める", "ERROR" not in cq.stdout and "dr_l2" in cq.stdout, cq.stdout.strip().splitlines()[0][:100]
      if cq.stdout.strip() else cq.stderr[-200:])

# (6) contur_v2 の入口 + M5 の定義
from forge_design.feedback import deltastar_integral as DI  # noqa: E402
rwv, cphi, dr0 = 3.0, np.cos(np.deg2rad(20.0)), 0.12
dstar_n = (rwv * dr0 - 0.5 * dr0 ** 2) * cphi / rwv       # 2π q_w (r_w δ_r − δ_r²/2) = 2π ρu r_w δ*_n、q_w = ρu cos φ
check("(6) 環状の定義 (M5) の逆変換", abs(float(bd.delta_r_annular(dstar_n, cphi, rwv)) - dr0) < 1e-14,
      f"{float(bd.delta_r_annular(dstar_n, cphi, rwv))!r}")

C68 = ROOT / "case/68.bell_two_mode"
NS68, EU68, PB68 = C68 / "run_0019_tp_m_ns_l100", C68 / "run_0018_tp_m_euler_l100_pin", C68 / "problem_bell_tp_m_ns_l100.yaml"
if not ((NS68 / "res_60000.h5").is_file() and (EU68 / "res_30000.h5").is_file()):
    print("[-- ] (6) ③の run (case/68 run_0019/0018) が手元に無い — 実データの試験を飛ばす")
else:
    from forge_design.evaluate.runner import build_wall
    from forge_design.evaluate.runner_axismach import contur_mu_fn
    from forge_design.probdef import load_problem
    p = load_problem(PB68); wall, _ = build_wall(p); gas = p.gas_model; rt = float(p.spec["r_throat"])
    ed = bd.euler_wall_edge(EU68, gas)
    m = ed["x"] <= wall.x_e
    wtbl = np.c_[ed["x"][m], ed["r"][m], np.zeros(int(m.sum())), ed["M"][m]]
    tbl = np.c_[ed["x"][m], ed["M"][m]]
    args = (wall, wtbl, gas, p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]), rt)
    kw = dict(closure_version="contur_v2", mu_fn=contur_mu_fn(p), x_virtual_m=rt, x_out=np.linspace(-3.9, 22.6, 200))
    ec = DI.EdgeConditions(*args, version="contur_v2", mu_fn=kw["mu_fn"], edge_M_full=tbl)
    check("(6) edge_M_full: 縁の M は表の線形補間そのもの", np.array_equal(ec._M, np.maximum(np.interp(ec._xg, tbl[:, 0], tbl[:, 1]), 1e-3)))
    for name, bad in (("範囲外", np.c_[tbl[5:, 0], tbl[5:, 1]]), ("contur_v1", tbl)):
        try:
            DI.EdgeConditions(*args, version=("contur_v1" if name == "contur_v1" else "contur_v2"), mu_fn=kw["mu_fn"], edge_M_full=bad)
            check(f"(6) edge_M_full は{name}を拒む", False, "例外が出なかった")
        except ValueError as e:
            check(f"(6) edge_M_full は{name}を拒む", True, str(e)[:60])
    base_di = load_base("design/forge_design/feedback/deltastar_integral.py", "forge_design.feedback._di_base")
    if base_di is None:
        check("(6) 基準コミットの deltastar_integral.py が読める", False, BASE_COMMIT)
    else:
        a = base_di.integral_bl(*args, **kw); b = DI.integral_bl(*args, **kw)
        r_ = same(a, b)
        check("(6) contur_v2 の既定の経路 (edge_M_full 無し) がビット同一 (③の壁・Euler の M)", r_ is None, r_ or "")
        a = base_di.integral_bl(*args, x_out=kw["x_out"]); b = DI.integral_bl(*args, x_out=kw["x_out"])
        r_ = same(a, b)
        check("(6) contur_v1 の既定の経路がビット同一 (③の壁)", r_ is None, r_ or "")
    c = bd.contur_bell(PB68, EU68, x_out=np.linspace(0.5, 22.6, 50))
    check("(6) contur_v2 が③の条件で回る (全域で有限・正、下流で単調増加)",
          np.all(np.isfinite(c["delta_r_ann"])) and np.all(c["delta_r_ann"] > 0) and np.all(np.diff(c["delta_r_ann"][c["x"] > 1.0]) > 0)
          and c["settings"].get("edge", "").startswith("table_full"),
          f"x 0.5 {c['delta_r_ann'][0]:.5f}・出口 {c['delta_r_ann'][-1]:.5f}")
    cmp_ = bd.remap_comparison(NS68, EU68, stations=(1.0,))
    r1 = cmp_["rows"][0]
    check("(6) ③の抽出が codex の再計算を再現 (x 1.0: 列 0.0031796 → 一定 x 0.0041637)",
          abs(r1["raw_column"] - 0.0031796) < 5e-8 and abs(r1["raw_constx"] - 0.0041637) < 5e-8,
          f"列 {r1['raw_column']:.7f}・一定 x {r1['raw_constx']:.7f}")

import shutil  # noqa: E402
shutil.rmtree(TMP, ignore_errors=True)
for d in _BASE_TMP:
    shutil.rmtree(d, ignore_errors=True)
print(f"FAIL 件数: {len(fails)}" + (f" ({', '.join(fails)})" if fails else ""))
sys.exit(1 if fails else 0)
