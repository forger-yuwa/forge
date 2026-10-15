"""CONTUR の較正の道具 `feedback/contur_calibrate.py` の試験 (plan tooling-nozzle-card-local-mesh §4.3・§5.1 #6a)。

合成のデータ:
 (1) 閉形式の 3 係数: 形の残差の無い δ_E から目標の係数を取り戻す (相対 1e−10)。残差があっても端点 (スロート・出口) は等式、
     c₁ は台形の重みの拘束つき最小二乗の解 (KKT の連立を一般の解法で解いた値と相対 1e−9、目的関数が c₁ ± ε で増える)。
 (2) Δx 重み: 同じ関数を列の寄せ方だけ変えて標本化しても c₁ がほぼ同じ (等重みより近い)。
 (3) 合成の run: 2 係数・3 係数・pass 0 (係数なし) の壁から目標の係数を取り戻す。停止規則 v2 の抽出の表を探す・sha256 の照合・
     状態が ok でない・held・delta_r_csv の壁・係数の不一致で止める。overrides の断片が YAML で同じ値に戻り、生成器の検査を通る。
 (4) CLI: 書く・--no-write・止めたら終了コード 1。
実データ (手元に記録の写しがあるとき):
 (5) case/71 run_0008 (2 係数の壁) → §4.3 の 3 係数 (c₀ 0.6779482479245829・c₁ 0.021036335310539667 /m・δ₀ −3.057205936466481e−4 m)
     を相対 1e−9 以内で再現。誤差の要約は §4.3 (旧 NS への射影: 最大 0.612 mm @ x 1.07 m、x ≥ 3 m 0.198 mm、2 係数は最大 2.14 mm)。
 (6) case/71 run_0010 (3 係数の壁) で較正し直し: 最大 0.624 → 0.612 mm・x ≥ 3 m 0.257 → 0.214 mm (plan §5.1 #3t)。
使い方: design/.venv-opt/bin/python design/tests/run_contur_calibrate_tests.py  (終了コード 0 全合格 / 1 不合格 / 2 実データを飛ばした)
"""
import io
import json
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from forge_design.feedback import contur_calibrate as CC  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
import contur_calibrate_fixture as FX  # noqa: E402

FAIL, SKIPPED = [], []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------ (1) 閉形式
print("== (1) 閉形式の 3 係数")
S = FX.S_DEF
x_m = FX.x_columns(with_zero=True) * S
C = FX.contur_m(x_m)
tgt = (0.68, 0.021, -3.06e-4)
E = (tgt[0] + tgt[1] * x_m) * C + tgt[2]
r = CC.calibrate_coefficients(x_m, C, E)
check("(1) 形の残差の無い δ_E から目標の係数を取り戻す (x = 0 の列のある格子、相対 1e−12)",
      rel(r["c0"], tgt[0]) < 1e-12 and rel(r["c1_per_m"], tgt[1]) < 1e-12 and rel(r["d0_m"], tgt[2]) < 1e-12,
      f"{r['c0']!r}, {r['c1_per_m']!r}, {r['d0_m']!r}")
Eb = E + FX.bump_m(x_m, 3e-4)
rb = CC.calibrate_coefficients(x_m, C, Eb)
db = CC.model_delta(x_m, C, rb["c0"], rb["c1_per_m"], rb["d0_m"])
et = abs(rb["c0"] * np.interp(0.0, x_m, C) + rb["d0_m"] - np.interp(0.0, x_m, Eb)); ee = abs(db[-1] - Eb[-1])
check("(1) 残差があってもスロート (c₀·C_t + δ₀ = E_t、C_t・E_t は x = 0 の補間)・出口 (最後の列) は等式 (≤ 1e−15 m)",
      et <= 1e-15 and ee <= 1e-15, f"{et:.1e}, {ee:.1e} m")
xq = FX.x_columns() * S
rq = CC.calibrate_coefficients(xq, FX.contur_m(xq), (tgt[0] + tgt[1] * xq) * FX.contur_m(xq) + tgt[2])
check("(1) 0 の列の無い格子 (実際の表) でも相対 1e−8 (スロートの値は両隣の列の直線補間なので厳密ではない)",
      max(rel(rq["c0"], tgt[0]), rel(rq["c1_per_m"], tgt[1]), rel(rq["d0_m"], tgt[2])) < 1e-8,
      f"相対 {max(rel(rq['c0'], tgt[0]), rel(rq['c1_per_m'], tgt[1]), rel(rq['d0_m'], tgt[2])):.1e}")
# KKT: min Σ w (A p − E)² s.t. G p = h、p = (c0, c1, d0)、A = [C, x C, 1]
m = x_m >= 0
w = CC.trapezoid_weights(x_m[m])
A = np.c_[C[m], x_m[m] * C[m], np.ones(m.sum())]
Ct, Ce, L = np.interp(0.0, x_m, C), C[-1], x_m[-1]
G = np.array([[Ct, 0.0, 1.0], [Ce, L * Ce, 1.0]])
h = np.array([np.interp(0.0, x_m, Eb), Eb[-1]])
K = np.block([[2 * A.T @ (w[:, None] * A), G.T], [G, np.zeros((2, 2))]])
sol = np.linalg.solve(K, np.r_[2 * A.T @ (w * Eb[m]), h])
check("(1) 閉形式 = 拘束つき重み付き最小二乗の KKT の解 (相対 1e−9)",
      rel(rb["c0"], sol[0]) < 1e-9 and rel(rb["c1_per_m"], sol[1]) < 1e-9 and rel(rb["d0_m"], sol[2]) < 1e-9,
      f"c₁ {rb['c1_per_m']!r} vs {sol[1]!r}")


def obj(c1):
    a = (h[1] - h[0]) / (Ce - Ct)
    c0 = a - c1 * L * Ce / (Ce - Ct)
    d0 = h[0] - c0 * Ct
    return float(np.sum(w * (CC.model_delta(x_m[m], C[m], c0, c1, d0) - Eb[m]) ** 2))


f0 = obj(rb["c1_per_m"])
check("(1) 目的関数は c₁ ± 1e−4 で増える (最小)", obj(rb["c1_per_m"] + 1e-4) > f0 and obj(rb["c1_per_m"] - 1e-4) > f0)
for bad, why in (((x_m[::-1], C[::-1], E[::-1]), "昇順でない"), ((x_m[x_m > 0], C[x_m > 0], E[x_m > 0]), "x = 0 を挟まない"),
                 ((x_m, np.full_like(C, 1e-3), E), "C_e = C_t")):
    try:
        CC.calibrate_coefficients(*bad)
        check(f"(1) 不正な入力 ({why}) は ValueError", False, "例外が出ない")
    except ValueError:
        check(f"(1) 不正な入力 ({why}) は ValueError", True)

# ------------------------------------------------------------------ (2) Δx 重み
print("== (2) Δx 重み")
xu = FX.x_columns(clustered=False, n_down=6000) * S
xc = FX.x_columns(clustered=True) * S
res = {}
for nm, xx in (("uniform", xu), ("clustered", xc)):
    Cx = FX.contur_m(xx)
    Ex = (tgt[0] + tgt[1] * xx) * Cx + tgt[2] + FX.bump_m(xx, 3e-4)
    res[nm] = CC.calibrate_coefficients(xx, Cx, Ex)["c1_per_m"]
    mm = xx >= 0
    # 等重みの c₁ (同じ式で w = 1)
    Ct_, Et_, Ce_, Ee_, L_ = np.interp(0, xx, Cx), np.interp(0, xx, Ex), Cx[-1], Ex[-1], xx[-1]
    a_ = (Ee_ - Et_) / (Ce_ - Ct_)
    b_ = xx[mm] * Cx[mm] - L_ * Ce_ * (Cx[mm] - Ct_) / (Ce_ - Ct_)
    r_ = Ex[mm] - Et_ - a_ * (Cx[mm] - Ct_)
    res[nm + "_eq"] = float(np.sum(b_ * r_) / np.sum(b_ * b_))
d_trap = rel(res["clustered"], res["uniform"]); d_eq = rel(res["clustered_eq"], res["uniform_eq"])
check("(2) 列の寄せ方を変えても台形の重みの c₁ はほぼ同じ (相対 < 1e−3) で、等重みより近い", d_trap < 1e-3 and d_trap < d_eq,
      f"台形 {d_trap:.2e}、等重み {d_eq:.2e}")

# ------------------------------------------------------------------ (3) 合成の run
print("== (3) 合成の run")
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    for nm, coefs in (("2 係数", (0.85, 0.0, -3.2e-4)), ("3 係数", (0.6779, 0.0210, -3.06e-4)), ("pass 0", None)):
        run = FX.make_run(td / f"run_{nm.replace(' ', '')}", coefs=coefs, target=tgt)
        rec = CC.calibrate_run(run)
        co = rec["coefficients"]
        ok = rel(co["c0"], tgt[0]) < 1e-8 and rel(co["c1_per_m"], tgt[1]) < 1e-8 and rel(co["d0_m"], tgt[2]) < 1e-8
        check(f"(3) {nm} の壁から目標の係数を取り戻す (相対 1e−8)", ok, f"{co}")
        check(f"(3) {nm}: 今の係数を記録 (pass 0 は (1, 0, 0))", (rec["current"]["c0"], rec["current"]["c1_per_m"], rec["current"]["d0_m"])
              == (coefs if coefs is not None else (1.0, 0.0, 0.0)) and rec["current"]["pass0"] == (coefs is None))
        pj = rec["errors"]["projected_after_calibration"]
        check(f"(3) {nm}: 射影の誤差 ≈ 0 (形の残差が無い; ≤ 1e−6 mm)", pj["max_abs_mm"] < 1e-6, f"{pj['max_abs_mm']:.2e} mm")
        dc = rec["checks"]["delta_C_vs_recorded_contur"]
        check(f"(3) {nm}: δ_in から戻した δ_C が記録の CONTUR のスロート・出口と合う (≤ 1e−9 mm)",
              abs(dc["throat"]["diff_mm"]) < 1e-9 and abs(dc["exit"]["diff_mm"]) < 1e-9, f"{dc['throat']['diff_mm']:.1e}, {dc['exit']['diff_mm']:.1e}")
    # 残差つき: 誤差の要約
    run = FX.make_run(td / "run_bump", coefs=(0.85, 0.0, -3.2e-4), target=tgt, bump=3e-4)
    rec = CC.calibrate_run(run)
    cur, pj = rec["errors"]["current_wall"], rec["errors"]["projected_after_calibration"]
    x_m2 = FX.x_columns() * S
    din2 = 0.85 * FX.contur_m(x_m2) - 3.2e-4
    E2 = (tgt[0] + tgt[1] * x_m2) * FX.contur_m(x_m2) + tgt[2] + FX.bump_m(x_m2, 3e-4)
    want = np.max(np.abs((din2 - E2)[x_m2 >= 0])) * 1e3
    check("(3) 今の壁の誤差の最大は |δ_in − δ_E| (x ≥ 0) の最大", abs(cur["max_abs_mm"] - want) < 1e-9, f"{cur['max_abs_mm']:.6f} vs {want:.6f} mm")
    check("(3) 射影の誤差はスロート・出口で 0、最大は形の残差の大きさ以下",
          abs(pj["throat_mm"]) < 1e-6 and abs(pj["exit_mm"]) < 1e-9 and pj["max_abs_mm"] <= 0.3 + 1e-9, f"{pj}")
    # overrides の断片
    ov = yaml.safe_load(rec["overrides_yaml"])["overrides"]
    check("(3) overrides の断片が YAML で同じ値 (float の往復) に戻り、model: contur を含む",
          ov == rec["overrides"] and ov["deltastar_initializer.model"] == "contur" and set(ov) == set(CC.OVERRIDE_KEYS))
    init_cfg = {k.split(".", 1)[1]: v for k, v in ov.items()}
    two = RA._two_coef(init_cfg)
    check("(3) 生成器の検査 (runner_axismach._two_coef・_slope_coef) を通る",
          two == (ov["deltastar_initializer.delta_scale"], ov["deltastar_initializer.delta_offset_m"])
          and RA._slope_coef(init_cfg) == ov["deltastar_initializer.delta_scale_slope_per_m"])
    # 停止規則 v2 の抽出の表
    run = FX.make_run(td / "run_v2", coefs=(0.85, 0.0, -3.2e-4), target=tgt, v2=True)
    loc = CC.locate_delta_e(run)
    check("(3) 停止規則 v2 の抽出の表を探す (状態の dir が別の機械のパスでも run の stop_rule_v2/delta_E_<k>/ に戻る)",
          loc["path"].endswith("stop_rule_v2/delta_E_46202/delta_r_next.csv") and "停止規則 v2" in loc["source"], loc["source"])
    rec = CC.calibrate_run(run)
    check("(3) v2 の表から目標の係数を取り戻す", rel(rec["coefficients"]["c1_per_m"], tgt[1]) < 1e-8)
    (run / "delta_r_next.csv").write_text("x_rt,delta_target,delta_in_prev,delta_E,held\n0,0,0,0,0\n")
    check("(3) 状態に δ_E の記録があれば run 直下の表より優先", "停止規則 v2" in CC.locate_delta_e(run)["source"])
    st = json.loads((run / "STOP_RULE_STATE.json").read_text())
    st["delta_E"]["delta_r_next_sha256"] = "0" * 64
    (run / "STOP_RULE_STATE.json").write_text(json.dumps(st))
    for nm, mut in (("sha256 が違う", None), ("状態が ok でない", "failed")):
        if mut is not None:
            st = json.loads((run / "STOP_RULE_STATE.json").read_text())
            st["delta_E"]["status"] = mut
            (run / "STOP_RULE_STATE.json").write_text(json.dumps(st))
        try:
            CC.locate_delta_e(run)
            check(f"(3) v2 の表の {nm} → ValueError", False, "例外が出ない")
        except ValueError as e:
            check(f"(3) v2 の表の {nm} → ValueError", True, str(e)[:80])
    check("(3) --delta-next の明示は状態より優先", CC.locate_delta_e(run, run / "delta_r_next.csv")["source"].startswith("explicit"))
    for nm, kw in (("x ≥ 0 に held", {"held_at": 30.0}), ("delta_r_csv の壁", {"delta_r_csv": True}),
                   ("prepare_info と delta_r_initial.json の係数が違う", {"mismatch_json": True})):
        run = FX.make_run(td / f"run_bad_{len(nm)}", target=tgt, **kw)
        try:
            CC.calibrate_run(run)
            check(f"(3) {nm} → ValueError", False, "例外が出ない")
        except ValueError as e:
            check(f"(3) {nm} → ValueError", True, str(e)[:80])
    run = FX.make_run(td / "run_noinit", target=tgt, initializer_in_prepare=False)
    (run / "delta_r_initial.json").unlink()
    try:
        CC.calibrate_run(run)
        check("(3) 積分法の初期壁でない (initializer が無い) → ValueError", False, "例外が出ない")
    except ValueError:
        check("(3) 積分法の初期壁でない (initializer が無い) → ValueError", True)
    run = FX.make_run(td / "run_held_up", target=tgt, held_at=-5.0)
    check("(3) 上流 (x < 0) の held は較正に使わないので止めない", rel(CC.calibrate_run(run)["coefficients"]["c0"], tgt[0]) < 1e-8)
    # --delta-in function
    run = FX.make_run(td / "run_fun", coefs=(0.85, 0.0, -3.2e-4), target=tgt)
    rf = CC.calibrate_run(run, delta_in="function")
    check("(3) --delta-in function (関数の値の 5 次補間) でも目標の係数 (相対 1e−6)", rel(rf["coefficients"]["c1_per_m"], tgt[1]) < 1e-6
          and rf["inputs"]["delta_in"]["source"] == "function", f"{rf['coefficients']}")
    chk = rf["checks"]["delta_in"]
    check("(3) 格子と関数の δ_in の差を記録する", chk.get("max_abs_grid_minus_function_x_ge_0_mm") is not None and not chk["warn"],
          f"{chk.get('max_abs_grid_minus_function_x_ge_0_mm')}")

    # ------------------------------------------------------------------ (4) CLI
    print("== (4) CLI")
    run = FX.make_run(td / "run_cli", coefs=(0.85, 0.0, -3.2e-4), target=tgt, bump=1e-4)
    out = io.StringIO()
    with redirect_stdout(out):
        rc = CC.main([str(run)])
    js = run / CC.OUT_DIRNAME / CC.OUT_JSON
    check("(4) 既定で <run>/contur_calibrate/ に JSON と overrides の YAML を書き終了コード 0",
          rc == 0 and js.is_file() and (run / CC.OUT_DIRNAME / CC.OUT_YAML).is_file()
          and json.loads(js.read_text())["format"] == CC.FORMAT)
    od = td / "elsewhere"
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        rc2 = CC.main([str(run), "--out-dir", str(od), "--json"])
    check("(4) --out-dir に書く", rc2 == 0 and (od / CC.OUT_JSON).is_file())
    run2 = FX.make_run(td / "run_cli_nowrite", target=tgt)
    with redirect_stdout(io.StringIO()):
        rc3 = CC.main([str(run2), "--no-write"])
    check("(4) --no-write は書かない", rc3 == 0 and not (run2 / CC.OUT_DIRNAME).exists())
    err = io.StringIO()
    with redirect_stderr(err), redirect_stdout(io.StringIO()):
        rc4 = CC.main([str(td / "no_such_run")])
    check("(4) 入力の不備は終了コード 1 と理由", rc4 == 1 and "止めた" in err.getvalue(), err.getvalue().strip()[:80])

# ------------------------------------------------------------------ (5)(6) 実データ
print("== (5)(6) case/71 の実データ")
R = ROOT / "case/71.m55_combustion_wt/_runs"
r08, r10 = R / "run_0008_ns_b0", R / "run_0010_ns_b0"
if (r08 / "prepare_info.json").is_file() and (r08 / "judge/extract/delta_r_next.csv.gz").is_file():
    rec = CC.calibrate_run(r08, delta_next=r08 / "judge/extract/delta_r_next.csv.gz")
    co = rec["coefficients"]
    T = (0.6779482479245829, 0.021036335310539667, -3.057205936466481e-04)
    rr = (rel(co["c0"], T[0]), rel(co["c1_per_m"], T[1]), rel(co["d0_m"], T[2]))
    check("(5) run_0008 → §4.3 の 3 係数を相対 1e−9 以内で再現", max(rr) <= 1e-9, f"相対差 {rr[0]:.1e}・{rr[1]:.1e}・{rr[2]:.1e}")
    cur, pj = rec["errors"]["current_wall"], rec["errors"]["projected_after_calibration"]
    check("(5) run_0008 の 2 係数の壁の誤差 最大 2.14 mm", round(cur["max_abs_mm"], 2) == 2.14, f"{cur['max_abs_mm']:.4f} mm @ {cur['x_m_at_max']:.2f} m")
    check("(5) 旧 NS への射影 最大 0.612 mm @ x 1.07 m・x ≥ 3 m 0.198 mm (§4.3)",
          round(pj["max_abs_mm"], 3) == 0.612 and round(pj["x_m_at_max"], 2) == 1.07 and round(pj["max_abs_x_ge_3m_mm"], 3) == 0.198,
          f"{pj['max_abs_mm']:.4f} mm @ {pj['x_m_at_max']:.3f} m、x ≥ 3 m {pj['max_abs_x_ge_3m_mm']:.4f} mm")
    check("(5) run_0008 の今の係数は 2 係数のカードの値 (c 0.852087987522164・δ₀ −3.220418618084591e−4 m)",
          rec["current"]["c0"] == 0.852087987522164 and rec["current"]["d0_m"] == -3.220418618084591e-04 and rec["current"]["c1_per_m"] == 0.0)
    rf = CC.calibrate_run(r08, delta_next=r08 / "judge/extract/delta_r_next.csv.gz", delta_in="function")
    rf_rel = max(rel(rf["coefficients"][k], v) for k, v in zip(("c0", "c1_per_m", "d0_m"), T))
    check("(5) 参考: 関数の δ_in で求めると §4.3 の値から相対 ≤ 1e−3 (格子と関数の差 ~1e−4 mm の分だけ違う)", rf_rel <= 1e-3, f"相対 {rf_rel:.1e}")
else:
    SKIPPED.append(f"(5) {r08.relative_to(ROOT)} の記録の写しが無い")
if (r10 / "prepare_info.json").is_file() and (r10 / "judge/extract/delta_r_next.csv.gz").is_file():
    rec = CC.calibrate_run(r10, delta_next=r10 / "judge/extract/delta_r_next.csv.gz")
    cur, pj = rec["errors"]["current_wall"], rec["errors"]["projected_after_calibration"]
    check("(6) run_0010 の較正し直し: 最大 0.624 → 0.612 mm・x ≥ 3 m 0.257 → 0.214 mm (plan #3t)",
          round(cur["max_abs_mm"], 3) == 0.624 and round(pj["max_abs_mm"], 3) == 0.612
          and round(cur["max_abs_x_ge_3m_mm"], 3) == 0.257 and round(pj["max_abs_x_ge_3m_mm"], 3) == 0.214,
          f"{cur['max_abs_mm']:.4f} → {pj['max_abs_mm']:.4f} mm、x ≥ 3 m {cur['max_abs_x_ge_3m_mm']:.4f} → {pj['max_abs_x_ge_3m_mm']:.4f} mm")
    check("(6) run_0010 の今の係数は 3 係数のカードの値", (rec["current"]["c0"], rec["current"]["c1_per_m"], rec["current"]["d0_m"])
          == (0.6779482479245829, 0.021036335310539667, -3.057205936466481e-04))
    check("(6) 最大の位置 x 1.06 m (plan #3t)", round(cur["x_m_at_max"], 2) == 1.06, f"{cur['x_m_at_max']:.3f} m")
else:
    SKIPPED.append(f"(6) {r10.relative_to(ROOT)} の記録の写しが無い")

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if not FAIL else f'{len(FAIL)} FAILURES: {FAIL}'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
