"""e4_recal_eval.py と e4_recal.py (問題の検査・段 2 の問題の作成・step 数の同期) の試験 (plan verification-case45-euler-total-enthalpy §6 E4、
登録 99431498)。forge・AWS は使わない (合成した出口 M・全温の主指標・残差 CSV・格子だけ。E4 の run の結果は読まない)。
check_convergence.py と check_quasisteady.py は本物を呼ぶ (出力の書式を合成で真似ない)。
確かめること:
  段 1 の判定: 13 枚すべて |M_common − 6| ≤ 1e-4 → 据え置き / 1 枚でも外れる → 更新 (δ₁ = δ₀ − (13 枚の平均 − 6)、末尾 5 枚の平均による値も記録) /
    前提未達・窓の欠け → 判別不能。境界は float の ≤ そのまま (6.0001 − 6 = 9.999999999976694e-05 は内、6.00010001 は外)。窓の外の時点は判定に入らない。
  段 2 の判定: 段 1 が更新・Md_moc_offset が δ₁ と一致・前提成立で 13 枚すべて内 → 合格 / 外 → 保留 (補正を重ねない) / 一致しない・段 1 が更新でない → 判別不能。
  出口の評価量: 線形の M(η) で M_common = 共通の η の列での値の平均 (線形補間は厳密)、M_own は帯内 (端を含む) の節点の平均、η の不正は例外。
  共通の η の列: NS の問題から 12 点・登録の sha256、すべて [0.05, 0.7]。NS の基準 run の nozzle.h5 (合成) と照合: 一致 → OK、
    帯内の点の数が違う・|Δη| > 1e-6 → 不成立、run が無い → 未確認 (不成立にしない)。
  出口 M の前提: 一定 + 微小な揺れ → 成立 / 13 枚の幅 6e-5 → 不成立 / ドリフト (STEADY でない) → 不成立 / 54000 の欠け → 不成立。
  全温の健全性: 13 枚すべて ±1 K (境界 1.0 K は内) → 成立、1 枚 1.01 K → 不成立。
  結合 (evaluate を合成の run で): 段 1 更新 → 段 2 待ち → 段 2 合格 (最終 = 更新して合格)、段 2 の Md_moc_offset の食い違い → 判別不能、
    段 1 据え置き、段 1 の残差に RISING → 判別不能、forge の sha256 の食い違い → 判別不能。
  e4_recal: check-problem (d0 の実物は OK、MOC のキー・mesh_euler・Md_moc_offset・mono_r2 の改変は止まる)、make-d1 (段 1 が更新で δ が一致したときだけ
    作る、YAML の float の書き方)、step 数の食い違いで run が forge の前に止まる。
usage: python3 test_e4_recal_eval.py
"""
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import e4_recal_eval as EV  # noqa: E402
import euler_t0_e2_eval as EV2  # noqa: E402

sys.path.insert(0, str(EV.TOOLS))
from stage_manifest import StageManifest  # noqa: E402

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def series(fn):
    return {s: fn(s) for s in EV.MAIN_STEPS}


W = EV.WIN13
check("窓: 判定窓 42000〜54000 (13 枚)・末尾 50000〜54000・直前 45000〜49000・本段 55 枚 (54000)",
      W == tuple(range(42000, 54001, 1000)) and EV.TAIL5 == tuple(range(50000, 54001, 1000)) and EV.PREV5 == tuple(range(45000, 49001, 1000))
      and len(EV.MAIN_STEPS) == 55 and EV.MAIN_NSTEPS == 54000 and EV.window_problems() == [])

# --- 段 1 の判定 -------------------------------------------------------------------------------------------------------------
r = EV.judge_stage1(series(lambda s: 6.00005), True)
check("段 1: 13 枚すべて +5e-5 → 据え置き (δ₀ を採用)", r["verdict"] == EV.LBL_KEEP and r["adopted_delta"] == EV.DELTA0 == 3.770e-4)
r = EV.judge_stage1(series(lambda s: 6.0 + 0.99e-4), True)
check("段 1: 境界の内側 (6 + 0.99e-4) → 据え置き", r["verdict"] == EV.LBL_KEEP)
r = EV.judge_stage1(series(lambda s: 6.0001), True)
check("段 1: 境界 6.0001 は float で |·| = 9.999999999976694e-05 ≤ 1e-4 (≤ の比較そのまま、許容を足さない) → 据え置き",
      r["verdict"] == EV.LBL_KEEP and r["facts"]["max_abs_dev"] == 6.0001 - 6.0)
r = EV.judge_stage1(series(lambda s: 6.00010001), True)
check("段 1: 6.00010001 (|·| > 1e-4) → 更新", r["verdict"] == EV.LBL_UPDATE and r["facts"]["max_abs_dev"] > 1e-4)
s1 = series(lambda s: 6.00005 + (0.6e-4 if s == W[5] else 0.0))
r = EV.judge_stage1(s1, True)
m13 = float(np.mean([s1[s] for s in W]))
check("段 1: 13 枚のうち 1 枚だけ +1.1e-4 → 更新、δ₁ = δ₀ − (13 枚の平均 − 6)",
      r["verdict"] == EV.LBL_UPDATE and r["delta1"] == EV.DELTA0 - (m13 - 6.0) and r["delta1_repr"] == repr(r["delta1"])
      and r["facts"]["basis"] == "win13")
check("段 1: 末尾 5 枚の平均による δ₁ も記録", r["facts"]["delta1_if_tail5"] == EV.DELTA0 - (float(np.mean([s1[s] for s in EV.TAIL5])) - 6.0))
r = EV.judge_stage1(series(lambda s: 5.9998), True)
check("段 1: 負側 (5.9998) → 更新、δ₁ = δ₀ + 2e-4 (float の丸め 1e-14 以内)", r["verdict"] == EV.LBL_UPDATE and abs(r["delta1"] - (3.770e-4 + 2e-4)) < 1e-14)
r = EV.judge_stage1(series(lambda s: 6.00005 + (5e-4 if s == W[0] - 1000 else 0.0)), True)
check("段 1: 窓の直前 (41000) の外れは判定に入らない (据え置き)", r["verdict"] == EV.LBL_KEEP)
r = EV.judge_stage1(series(lambda s: 6.00005), False, ["[P3] 合成"])
check("段 1: 前提未達なら目標内でも判別不能 (理由と事実を残す)", r["verdict"] == EV.LBL_UNDET and "[P3] 合成" in r["reasons"]
      and r["facts"]["all_within_1e-4"])
sm = series(lambda s: 6.00005)
del sm[54000]
check("段 1: 判定窓の 54000 が欠ける → 判別不能", EV.judge_stage1(sm, True)["verdict"] == EV.LBL_UNDET)
sn = series(lambda s: 6.00005)
sn[W[3]] = math.nan
check("段 1: 判定窓に非有限 → 判別不能", EV.judge_stage1(sn, True)["verdict"] == EV.LBL_UNDET)

# --- 段 2 の判定 -------------------------------------------------------------------------------------------------------------
j1 = EV.judge_stage1(series(lambda s: 5.9998), True)
d1 = j1["delta1"]
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], j1, d1)
check("段 2: 段 1 更新・δ₁ 一致・13 枚すべて内 → 合格 (δ₁ を採用)", r["verdict"] == EV.LBL_PASS2 and r["adopted_delta"] == d1)
r = EV.judge_stage2(series(lambda s: 6.00002 + (2e-4 if s == W[-1] else 0.0)), True, [], j1, d1)
check("段 2: 1 枚でも外 → 保留 (補正を重ねない)", r["verdict"] == EV.LBL_HOLD2 and "重ねない" in EV.LBL_HOLD2)
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], j1, d1 + 1e-12)
check("段 2: run の Md_moc_offset が δ₁ と一致しない → 判別不能", r["verdict"] == EV.LBL_UNDET and any("一致しない" in x for x in r["reasons"]))
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], EV.judge_stage1(series(lambda s: 6.00005), True), EV.DELTA0)
check("段 2: 段 1 が据え置きなら段 2 は登録外 → 判別不能", r["verdict"] == EV.LBL_UNDET and any("登録外" in x for x in r["reasons"]))
r = EV.judge_stage2(series(lambda s: 6.00002), False, ["[P1] 合成"], j1, d1)
check("段 2: 前提未達 → 判別不能", r["verdict"] == EV.LBL_UNDET)

# E4V (plan §6 E4V、登録 6b5b2cdb): 段 1 が判別不能でも、facts の δ₁ の repr が登録の候補と一致し、run の値も一致するときだけ判定する
cand = float(EV.E4V_DELTA_CAND_REPR)
j1u = {"verdict": EV.LBL_UNDET, "reasons": ["前提の未達"], "facts": {"delta1": cand, "delta1_repr": EV.E4V_DELTA_CAND_REPR}}
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], j1u, cand)
check("E4V: 段 1 が判別不能・候補一致・全前提成立・目標内 → E4V 合格", r["verdict"] == EV.LBL_PASS_E4V and r["mode"] == "E4V")
r = EV.judge_stage2(series(lambda s: 6.0003), True, [], j1u, cand)
check("E4V: 全前提成立で目標外 → E4V 保留", r["verdict"] == EV.LBL_HOLD_E4V)
r = EV.judge_stage2(series(lambda s: 6.00002), False, ["[P2] 合成"], j1u, cand)
check("E4V: 前提の未達 → 判別不能", r["verdict"] == EV.LBL_UNDET)
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], j1u, cand + 1e-12)
check("E4V: run の値が候補と違う → 判別不能", r["verdict"] == EV.LBL_UNDET)
j1x = {"verdict": EV.LBL_UNDET, "reasons": ["前提の未達"], "facts": {"delta1": 1e-5, "delta1_repr": "1e-05"}}
r = EV.judge_stage2(series(lambda s: 6.00002), True, [], j1x, 1e-5)
check("E4V: 段 1 の δ₁ が登録の候補でない判別不能 → 段 2 は登録外 (判別不能)", r["verdict"] == EV.LBL_UNDET and "登録外" in " ".join(r["reasons"]))

# --- 出口の評価量 ---------------------------------------------------------------------------------------------------------------
eta_own = np.r_[0.0, np.linspace(0.05, 0.7, 14), 0.8, 0.9, 1.0]
eta_c = np.array([0.0887581, 0.2, 0.45, 0.6722])
M = 5.9 + 0.2 * eta_own
v = EV.exit_values(M, eta_own, eta_c)
check("出口: 線形の M(η) で M_common = 共通の η での値の平均 (線形補間は厳密)", abs(v["exitM_common"] - float(np.mean(5.9 + 0.2 * eta_c))) < 1e-14)
check("出口: M_own は帯内 (0.05・0.7 の端を含む) の節点の平均", v["n_own"] == 14 and abs(v["exitM_own"] - float(np.mean(5.9 + 0.2 * np.linspace(0.05, 0.7, 14)))) < 1e-14)
for bad_eta, why in ((eta_own + 0.01, "0 から始まらない"), (eta_own[::-1], "単調でない")):
    try:
        EV.exit_values(M, bad_eta, eta_c)
        ok_ = False
    except ValueError:
        ok_ = True
    check(f"出口: η が {why}と例外", ok_)

# --- 共通の η の列 ---------------------------------------------------------------------------------------------------------------
eta_reg = EV.eta_common_from_problem(C / EV.NS_REF_PROBLEM)
check("η の列: NS の問題から 12 点・登録の sha256・すべて [0.05, 0.7]・0.08876〜0.67224",
      eta_reg.size == EV.ETA_COMMON_N == 12 and EV.eta_digest(eta_reg) == EV.ETA_COMMON_SHA256
      and eta_reg.min() >= 0.05 and eta_reg.max() <= 0.7 and abs(eta_reg[0] - 0.0887581) < 1e-6 and abs(eta_reg[-1] - 0.6722423) < 1e-6)


def write_struct_h5(path: Path, X, R, extra_values=None):
    """構造格子 (ni, nj) の合成 nozzle.h5 (/MESH/COORD float32、VIZMESH/CONNE の構造格子の並び)。"""
    import h5py
    ni, nj = X.shape
    c = np.zeros((ni * nj, 3), dtype=np.float32)
    c[:, 0], c[:, 1] = X.ravel(), R.ravel()
    i, j = np.meshgrid(np.arange(ni - 1), np.arange(nj - 1), indexing="ij")
    base = (i * nj + j).ravel()
    q = np.stack([np.full(base.size, 9), base, base + nj, base + nj + 1, base + 1], axis=1)
    with h5py.File(path, "w") as f:
        f["/MESH/COORD"] = c.ravel()
        f["VIZMESH/CONNE"] = q.ravel().astype(np.int64)
        for k, v_ in (extra_values or {}).items():
            f[f"/VALUE/{k}"] = v_


def flat_ns_grid(scale_eta=1.0):
    sys.path.insert(0, str(EV.ROOT / "design"))
    from forge_design.evaluate import runner_axismach as RA
    from forge_design.meshing.mesh2d import generate_axisym_mesh
    mp = RA.mesh_params(RA.load_problem(C / EV.NS_REF_PROBLEM), 1.0, 561, 97, 4.5e-5)
    c, _, _ = generate_axisym_mesh(EV._FlatWall(), mp)
    X, R = c[:, 0].reshape(mp.ni, mp.nj), c[:, 1].reshape(mp.ni, mp.nj) * 0.0767531 * 9.4
    if scale_eta != 1.0:
        R[-1, 1:-1] *= scale_eta
    return X * 0.0767531, R


with tempfile.TemporaryDirectory() as td:
    case = Path(td)
    shutil.copy(C / EV.NS_REF_PROBLEM, case / EV.NS_REF_PROBLEM)
    eta, rec = EV.eta_common(case)
    check("η の列: NS の基準 run が無ければ照合は未確認 (不成立にしない)", not rec["problems"] and "未確認" in rec["xcheck_ns_run"]["status"])
    (case / EV.NS_REF_RUN).mkdir()
    X, R = flat_ns_grid()
    write_struct_h5(case / EV.NS_REF_RUN / "nozzle.h5", X, R)
    eta, rec = EV.eta_common(case)
    check(f"η の列: 合成の NS の基準 run (float32 の座標) と照合 OK (|Δη| {rec['xcheck_ns_run'].get('max_abs_deta')})",
          not rec["problems"] and rec["xcheck_ns_run"]["ok"] and rec["xcheck_ns_run"]["max_abs_deta"] <= 1e-6)
    X, R = flat_ns_grid(scale_eta=1.00001)
    write_struct_h5(case / EV.NS_REF_RUN / "nozzle.h5", X, R)
    eta, rec = EV.eta_common(case)
    check("η の列: NS の基準 run の η が 1e-5 相対でずれる → 照合が不成立 (固定の条件)", bool(rec["problems"]) and not rec["xcheck_ns_run"]["ok"])

# --- 出口 M の前提 (本物の check_quasisteady) -------------------------------------------------------------------------------------


def exit_case(root: Path, name: str, fn, steps=None):
    steps = EV.MAIN_STEPS if steps is None else steps
    ser = {s: fn(s) for s in steps}
    qs = {}
    for wn, ws in EV.WINDOWS.items():
        rows = [{"stage": "main", "step": s, "exitM_common": ser[s], "exitM_own": ser[s] + 1e-5} for s in ws if s in ser]
        p = root / f"{name}_{wn}.csv"
        EV.series_csv(p, rows, list(EV.EXIT_COLS))
        qs[wn] = EV2.run_quasisteady(p, list(EV.EXIT_COLS), root / f"{name}_{wn}.txt")
    return EV.exit_window_check(ser, qs)


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    rc = exit_case(root, "flat", lambda s: 6.00003 + 2e-6 * math.sin(1.3 * s / 1000.0))
    check(f"出口 M の前提: 一定 + 2e-6 の揺れ → 両窓 STEADY・幅・平均差で成立 ({rc['verdicts']})", rc["ok"])
    rc = exit_case(root, "range", lambda s: 6.00003 + (6e-5 if s == W[4] else 0.0))
    check("出口 M の前提: 13 枚の幅 6e-5 > 5e-5 → 不成立", (not rc["ok"]) and any("13 枚の幅" in p for p in rc["problems"]))
    rc = exit_case(root, "drift", lambda s: 6.0 + 2e-3 * s / 54000.0)
    check(f"出口 M の前提: ドリフト (VERDICT {rc['verdicts']}) → 不成立", (not rc["ok"]) and rc["range13"] > 5e-5)
    rc = exit_case(root, "edge", lambda s: 6.00003 + (5e-5 if s == W[4] else 0.0))
    check(f"出口 M の前提: 幅ちょうど 5e-5 (float {rc.get('range13')!r}) は ≤ の比較そのまま", rc["ok"] == (rc["range13"] <= 5e-5))
    rc = exit_case(root, "short", lambda s: 6.00003, steps=EV.MAIN_STEPS[:-1])
    check("出口 M の前提: 54000 が欠ける → 不成立", (not rc["ok"]) and any("欠ける" in p for p in rc["problems"]))
    ser = {s: 6.00003 for s in EV.MAIN_STEPS}
    rc = EV.exit_window_check(ser, {"win13": "", "tail5": ""})
    check("出口 M の前提: check_quasisteady の出力が読めない → 不成立", not rc["ok"])

# --- 全温の健全性 -----------------------------------------------------------------------------------------------------------------


def ind_const(base=0.0, spread=0.2):
    return {r_: {"n": 100, "max": base + spread, "min": base - spread, "q99": abs(base) + spread / 2, **{s: 0.0 for s in EV2.FRAC_STATS}}
            for r_ in EV2.REGION_KEYS}


check("全温: 13 枚すべて ±1.0 K (境界) → 成立", EV.t0_within_check(series(lambda s: ind_const(0.0, 1.0)))["ok"])
r = EV.t0_within_check(series(lambda s: ind_const(0.0, 1.0 + (0.01 if s == W[2] else 0.0))))
check("全温: 1 枚だけ 1.01 K → 不成立", not r["ok"] and abs(r["max_over_steps_K"] - 1.01) < 1e-12)
ib = series(lambda s: ind_const(0.0, 0.1))
del ib[W[7]]
check("全温: 判定窓の主指標が欠ける → 不成立", not EV.t0_within_check(ib)["ok"])

# --- 結合 (evaluate を合成の run で; 全温の復元と出口 M の読みは差し替え) ---------------------------------------------------------
RCOLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roY0", "rms_roY1"]
BC = ("inlet: {physID: 1, kind: inlet_Pressure, floats: {Pt: 5500000.0, Tt: 1600.0}}\n")
FSHA = "6b47811bad9daa6807fb3f1e0808a1472c8791b572c2b7c6053fd9e225fa57d8"


def cfg(cfl, conv, nsteps, out, relax=0.7):
    return ("solver: \"SLAU\"\nphysProp: {thermalMethod: 2}\ntime:\n  last: {nStepOuter: %d}\n  outStepInterval: %d\n"
            "  deltaT: {cfl: %s, cfl_pseudo: %s, blockDPLUR: 1, implicitRelax: %s}\n  nStepInner: 5\n"
            "space: {convMethod: %d, limiter: 2}\nturbulence: {model: \"none\"}\n" % (nsteps, out, cfl, cfl, relax, conv))


def resid(kind, n=300):
    i = np.arange(n)
    if kind == "stall":
        return 1e-3 * (1.0 + 0.02 * np.sin(1.7 * i))
    if kind == "rising":
        v = resid("stall", n)
        k = int(0.8 * n)
        v[k:] = 1e-3 * 10.0 ** (2.0 * (i[k:] - k) / (n - 1 - k))
        return v
    return 10.0 ** (-4.0 * i / (n - 1))


def make_run(case: Path, stage: str, delta, kinds=None, fsha=FSHA, soft_cfl="0.5"):
    rd = case / EV.RUNS[stage]
    rd.mkdir()
    sm = StageManifest(rd)
    sm.add("soft", cfg(soft_cfl, 0, EV.SOFT_STEPS, EV.SOFT_STEPS), BC, history="residual_history_soft.csv")
    main_cfg = cfg("2.0", 1, EV.MAIN_NSTEPS, EV.OUT_INTERVAL)
    sm.add("main", main_cfg, BC, history="residual_history_main.csv")
    sm.write()
    (rd / "solverConfig.yaml").write_text(main_cfg)
    (rd / "bcondConfig.yaml").write_text(BC)
    for tag, ns, kk in (("soft", EV.SOFT_STEPS, {}), ("main", EV.MAIN_NSTEPS, kinds or {})):
        cols = {c: resid(kk.get(c, "stall" if tag == "main" else "pass")) for c in RCOLS}
        n = len(cols[RCOLS[0]])
        stp = np.round(np.linspace(0, ns - 1, n)).astype(int)
        with open(rd / f"residual_history_{tag}.csv", "w") as f:
            f.write("step," + ",".join(RCOLS) + "\n")
            for k in range(n):
                f.write(f"{stp[k]}," + ",".join(f"{cols[c][k]:.6e}" for c in RCOLS) + "\n")
    r = subprocess.run([sys.executable, str(EV.TOOLS / "check_convergence.py"), str(rd), "--segment"], capture_output=True, text=True)
    (rd / EV.SEGMENT_VERDICT_FILE).write_text(r.stdout + r.stderr)
    (rd / EV.SOFT_DIR).mkdir()
    for p in (rd / "RUN_PROVENANCE.txt", rd / EV.SOFT_DIR / "RUN_PROVENANCE.txt"):
        p.write_text(f"forge_sha256: {fsha}\n")
    info = {"moc": {"axis_limit": "legacy", "corrector": "fixed2", "gate": {"applicable": False, "pass": None}},
            "wall_fit": {"mono_r2": [0.0, 1.5]}, "Md_moc_offset": delta, "scale_m": EV.SCALE_M,
            "initial_line": {"run": f"/x/{EV.INITIAL_LINE[0]}", "res": [EV.INITIAL_LINE[1]]},
            "mesh": {"ni": 4, "nj": 5, "source": "mesh_euler",
                     "params": {**EV.MESH_EULER_EXPECT, "scale": EV.SCALE_M, "local_center": 0.0, "local_refine": 1.0, "local_width": 0.75,
                                "wall_first_blend_x0": 0.5, "wall_first_blend_x1": 6.0, "wall_first_up_x0": None, "wall_first_up_x1": None}}}
    (rd / "prepare_info.json").write_text(json.dumps(info))
    files = {n: "same" for n in EV.SAME_AS_D0}
    (rd / EV.PREP_JSON).write_text(json.dumps({"stage": stage, "ic_check": {"VERDICT": "OK"}, "files_sha256": files}))
    return rd


x4 = np.array([-5.0, -1.0, 2.0, 90.0])
s5 = np.array([0.0, 0.1, 0.4, 0.75, 1.0])
GX = np.repeat(x4[:, None], s5.size, axis=1)
GR = np.array([3.0, 1.2, 1.5, 9.4])[:, None] * s5[None, :]


def run_evaluate(exit_by_stage: dict, stages=("d0",), kinds=None, delta_d1=None, fsha=FSHA, t0_spread=0.05):
    with tempfile.TemporaryDirectory() as td:
        case = Path(td)
        shutil.copy(C / EV.NS_REF_PROBLEM, case / EV.NS_REF_PROBLEM)
        (case / EV.REF_RUN_BIN).mkdir()
        (case / EV.REF_RUN_BIN / "RUN_PROVENANCE.txt").write_text(f"forge_sha256: {FSHA}\n")
        for st in stages:
            make_run(case, st, EV.DELTA0 if st == "d0" else delta_d1, (kinds or {}).get(st), fsha)
        sv = (EV2.load_geometry, EV2.eval_snapshot, EV.exit_snapshot)
        EV2.load_geometry = lambda rd: (GX, GR, 0.0768075)
        stage_of = {v_: k for k, v_ in EV.RUNS.items()}

        def fake_t0(rd, h5, X, R, eta, ep):
            return {"file": h5.name, "problems": [], "indicators": ind_const(0.0, t0_spread), "sha256": None}

        def fake_exit(h5, ni, nj, eta_last, eta_c):
            st = stage_of[h5.parent.name]
            step = int(h5.stem.split("_")[1])
            m = exit_by_stage[st](step)
            return {"file": h5.name, "problems": [], "exitM_common": m, "exitM_own": m + 3e-5, "n_own": 2}
        EV2.eval_snapshot, EV.exit_snapshot = fake_t0, fake_exit
        try:
            out = EV.evaluate(case)
        finally:
            EV2.load_geometry, EV2.eval_snapshot, EV.exit_snapshot = sv
        return out


flat = lambda m: (lambda s: m + 2e-6 * math.sin(1.3 * s / 1000.0))  # noqa: E731
out = run_evaluate({"d0": flat(5.99980)})
j = out["stages"]["d0"]["judgment"]
check(f"結合: 段 1 が -2e-4 で整定 → 前提成立・更新 (δ₁ {j.get('delta1_repr')})、最終は段 2 待ち",
      out["stages"]["d0"]["precondition"]["ok"] and j["verdict"] == EV.LBL_UPDATE and out["final"]["status"] == "段 2 待ち"
      and j["delta1"] == EV.DELTA0 - (float(np.mean([flat(5.99980)(s) for s in W])) - 6.0)
      and out["stages"]["d1"]["judgment"]["verdict"] == EV.LBL_NOTRUN)
check("結合: 出力に評価器の sha256・登録の commit 99431498・η の列・定義・停滞のみの表示",
      out["evaluator_sha256"] == EV._sha(Path(EV.__file__)) and out["plan_reg_commit"] == "6b5b2cdb" and out["plan_reg_commit_e4"] == "99431498"
      and out["eta_common"]["n"] == 12 and out["eta_common"]["sha256"] == EV.ETA_COMMON_SHA256
      and "PASS ではない" in out["stages"]["d0"]["convergence"]["label"])
d1v = j["delta1"]
out = run_evaluate({"d0": flat(5.99980), "d1": flat(6.00001)}, stages=("d0", "d1"), delta_d1=d1v)
check("結合: 段 2 (Md_moc_offset = δ₁) が +1e-5 → 合格、最終 = 更新して合格・Md_moc_offset δ₁・参照 run_0164",
      out["stages"]["d1"]["judgment"]["verdict"] == EV.LBL_PASS2 and out["final"]["status"] == "更新して合格"
      and out["final"]["Md_moc_offset"] == d1v and out["final"]["reference_run"] == EV.RUNS["d1"])
out = run_evaluate({"d0": flat(5.99980), "d1": flat(6.00001)}, stages=("d0", "d1"), delta_d1=d1v + 1e-9)
check("結合: 段 2 の Md_moc_offset が δ₁ と違う → 判別不能・最終は保留",
      out["stages"]["d1"]["judgment"]["verdict"] == EV.LBL_UNDET and out["final"]["status"] == "保留")
out = run_evaluate({"d0": flat(5.99980), "d1": flat(6.00030)}, stages=("d0", "d1"), delta_d1=d1v)
check("結合: 段 2 でも外れる → 保留 (補正を重ねない)", out["stages"]["d1"]["judgment"]["verdict"] == EV.LBL_HOLD2 and out["final"]["status"] == "保留")
out = run_evaluate({"d0": flat(6.00004)})
check("結合: 段 1 が +4e-5 → 据え置き (最終 Md_moc_offset = δ₀、参照 run_0163)",
      out["stages"]["d0"]["judgment"]["verdict"] == EV.LBL_KEEP and out["final"] == {"status": "据え置き", "Md_moc_offset": EV.DELTA0,
                                                                                      "reference_run": EV.RUNS["d0"]})
out = run_evaluate({"d0": flat(6.00004)}, kinds={"d0": {"rms_roUy": "rising"}})
check("結合: 段 1 の残差に RISING の列 → 前提未達・判別不能", out["stages"]["d0"]["judgment"]["verdict"] == EV.LBL_UNDET
      and any(x.startswith("[P1]") for x in out["stages"]["d0"]["precondition"]["reasons"]))
out = run_evaluate({"d0": flat(6.00004)}, fsha="0" * 64)
check("結合: forge の sha256 が run_0143 と違う → 判別不能", out["stages"]["d0"]["judgment"]["verdict"] == EV.LBL_UNDET
      and any("forge の sha256" in x for x in out["stages"]["d0"]["precondition"]["reasons"]))
out = run_evaluate({"d0": flat(6.00004)}, t0_spread=1.5)
check("結合: 全温が ±1.5 K → 前提未達・判別不能", out["stages"]["d0"]["judgment"]["verdict"] == EV.LBL_UNDET
      and any(x.startswith("[P2a]") for x in out["stages"]["d0"]["precondition"]["reasons"]))
out = run_evaluate({"d0": lambda s: 6.00004 + (7e-5 if s == W[6] else 0.0)})
check("結合: 出口 M の 13 枚の幅 7e-5 → 前提未達・判別不能", out["stages"]["d0"]["judgment"]["verdict"] == EV.LBL_UNDET
      and any(x.startswith("[P3]") for x in out["stages"]["d0"]["precondition"]["reasons"]))

# --- step 数の同期 -------------------------------------------------------------------------------------------------------------
ok_eff = {"nStepOuter": 54000, "outStepInterval": 1000, "manifest_main_nStepOuter": 54000, "residual_last_step": 53999}
check("同期: 実効値 54000・1000・manifest 54000・最終 step 53999 は食い違いなし", EV.step_problems(ok_eff) == [])
for k_, v_ in (("nStepOuter", 18000), ("outStepInterval", 500), ("manifest_main_nStepOuter", 18000), ("residual_last_step", 17999)):
    check(f"同期: {k_} = {v_} は食い違い", len(EV.step_problems({**ok_eff, k_: v_})) == 1)
with tempfile.TemporaryDirectory() as td:
    rd = Path(td) / "x"
    rd.mkdir()
    sm = StageManifest(rd)
    sm.add("soft", cfg("1.0", 0, EV.SOFT_STEPS, EV.SOFT_STEPS), BC, history="h_soft.csv")
    sm.add("main", cfg("2.0", 1, EV.MAIN_NSTEPS, EV.OUT_INTERVAL, relax=0.5), BC, history="h_main.csv")
    sm.write()
    _, mp_ = EV.manifest_problems(rd)
    check("同期: stage_manifest の soft の cfl 1.0・本段の implicitRelax 0.5 は不成立", len(mp_) == 3)

# --- e4_recal: 問題の検査・段 2 の問題の作成・run の前の停止 ----------------------------------------------------------------------
import e4_recal as RUN  # noqa: E402

check("check-problem: 実物の d0 は OK", RUN.check_problem("d0") == [])
with tempfile.TemporaryDirectory() as td:
    case = Path(td)
    shutil.copy(C / RUN.BASE_PROBLEM, case / RUN.BASE_PROBLEM)
    src = (C / EV.PROBLEMS["d0"]).read_text()
    for name, edit in (("MOC のキー", lambda t: t.replace("  wall_fit_mono_r2:", "  moc_corrector: converge\n  wall_fit_mono_r2:", 1)),
                       ("mesh_euler", lambda t: t.replace("  wall_first_frac: 0.005\n", "  wall_first_frac: 0.004\n", 1)),
                       ("Md_moc_offset", lambda t: t.replace("Md_moc_offset: 3.7700e-04", "Md_moc_offset: 3.7600e-04", 1)),
                       ("mono_r2", lambda t: t.replace("  wall_fit_mono_r2: [0.0, 1.5]", "  wall_fit_mono_r2: [0.0, 1.4]", 1))):
        t = edit(src)
        assert t != src, name
        (case / EV.PROBLEMS["d0"]).write_text(t)
        check(f"check-problem: d0 の {name} を変えると止まる", bool(RUN.check_problem("d0", case=case)))
    (case / EV.PROBLEMS["d0"]).write_text(src)
    check("check-problem: 写しの d0 は OK", RUN.check_problem("d0", case=case) == [])
    d = 2.1234567890123e-05
    ej = case / "eval.json"
    sha = EV._sha(C / "e4_recal_eval.py")
    ej.write_text(json.dumps({"evaluator_sha256": sha, "stages": {"d0": {"judgment": {"verdict": EV.LBL_KEEP}}}}))
    try:
        RUN.make_d1(repr(d), ej, case=case)
        ok_ = False
    except SystemExit:
        ok_ = True
    check("make-d1: 段 1 が据え置きなら作らない", ok_ and not (case / EV.PROBLEMS["d1"]).exists())
    ej.write_text(json.dumps({"evaluator_sha256": sha, "stages": {"d0": {"judgment": {"verdict": EV.LBL_UPDATE, "delta1": d,
                                                                                         "delta1_repr": repr(d)}}}}))
    try:
        RUN.make_d1(repr(d + 1e-12), ej, case=case)
        ok_ = False
    except SystemExit:
        ok_ = True
    check("make-d1: --delta が段 1 の δ₁ と違えば作らない", ok_ and not (case / EV.PROBLEMS["d1"]).exists())
    out_p = RUN.make_d1(repr(d), ej, case=case)
    doc = yaml.safe_load(out_p.read_text())
    check("make-d1: δ₁ の問題を作り、d0 と name・Md_moc_offset だけが違い、値は float で完全一致",
          doc["geometry"]["Md_moc_offset"] == d and isinstance(doc["geometry"]["Md_moc_offset"], float)
          and RUN.check_problem("d1", d, case=case) == [] and doc["name"].endswith("_d1"))
    check("make-d1: check-problem d1 は違う δ で止まる", bool(RUN.check_problem("d1", d * 2, case=case)))
check("yaml_float: 1e-05 → 1.0e-05 (PyYAML の float)・負の値・普通の値", RUN.yaml_float("1e-05") == "1.0e-05"
      and RUN.yaml_float(repr(-0.00012345)) == repr(-0.00012345) and RUN.yaml_float("0.000377") == "0.000377")
with tempfile.TemporaryDirectory() as td:
    rd = Path(td) / EV.RUNS["d0"]
    rd.mkdir()
    info = {"e4": {"stage": "d0", "Md_moc_offset": EV.DELTA0}, "ic": {"VERDICT": "OK"},
            "moc": {"axis_limit": "legacy", "corrector": "fixed2"}, "wall_fit": {"mono_r2": [0.0, 1.5]}, "Md_moc_offset": EV.DELTA0,
            "scale_m": EV.SCALE_M, "initial_line": {"run": EV.INITIAL_LINE[0], "res": [EV.INITIAL_LINE[1]]},
            "mesh": {"ni": 2000, "source": "mesh_euler", "params": dict(EV.MESH_EULER_EXPECT)}}
    (rd / "prepare_info.json").write_text(json.dumps(info))
    (rd / "solverConfig.yaml").write_text(cfg("2.0", 1, 18000, 1000))
    check("実行側: run の config が 18000 なら check_run_steps が食い違いを返す", bool(RUN.check_run_steps(rd)))
    try:
        RUN.run(rd)
        ok_, msg = False, "起動した"
    except SystemExit as e:
        ok_, msg = "同期していない" in str(e), str(e)
    check(f"実行側: run は step 数の食い違いで forge の前に止まる ({msg[:40]})", ok_)
    (rd / "solverConfig.yaml").write_text(cfg("2.0", 1, 54000, 1000))
    check("実行側: config が 54000・1000 なら食い違いなし", RUN.check_run_steps(rd) == [])
    info["mesh"]["source"] = "mesh"
    (rd / "prepare_info.json").write_text(json.dumps(info))
    try:
        RUN.run(rd)
        ok_ = False
    except SystemExit as e:
        ok_ = "固定の条件" in str(e)
    check("実行側: 格子の採用元が mesh (mesh_euler でない) なら forge の前に止まる", ok_)

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
