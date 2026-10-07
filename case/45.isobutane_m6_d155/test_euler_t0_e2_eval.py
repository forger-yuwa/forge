"""euler_t0_e2_eval.py と euler_t0_e2.py (step 数の同期) の試験 (plan verification-case45-euler-total-enthalpy §6 E2、再改訂 a59e61b5)。
forge・AWS・HDF5 は使わない
(合成した主指標・残差 CSV・格子だけ。E2 の run の結果は読まない)。check_convergence.py と check_quasisteady.py は本物を呼ぶ
(出力の書式を合成で真似ない)。
確かめること:
  判定の 3 区分: 支持 (A のスロート > 100 K・B が全領域 ±1 K) / 退ける (両側の異常・差 ≤ 1 K) / 判別不能 (A が再現しない・中間・前提未達)。
    「評価窓の全保存時点」: 13 枚のうち 1 枚だけ外れても判別不能。境界 (B の最大 = 1.0 K、差 = 1.0 K) は含む。
  前提 P1 (列ごとの内訳): PASS / 全列が停滞 → 成立。総合は plateau だが 1 列が RISING (b) / still converging / NaN / 区間が main でない /
    必須列の欠け / 判定ファイルなし → 不成立。停滞のみの表示を PASS と書かない。
  前提 P2 (a59e61b5: 絶対の幅だけ。VERDICT は併記): 300→304 K の 13 点は check_quasisteady (--tail 1) で STEADY だが、幅 4 K > 0.1 K で
    不成立 (a)。近零の偏差で STEADY でなくても (DRIFTING 等) 幅 ≤ 0.1 K なら成立。幅 0.2 K の近零の上昇は不成立。末尾 5 枚の窓も同じ。
    54000 を欠く窓 (13 枚・5 枚) は不成立。
  step 数の同期: 評価器の窓は MAIN_NSTEPS (54000) の最後の 13 枚・5 枚。solverConfig の nStepOuter・outStepInterval、stage_manifest の main、
    本段の残差の最終 step が食い違えば、評価は固定の条件の不成立 (判別不能)、実行側 (euler_t0_e2.run) は起動前に止まる。
    評価器の定数だけを書き換えた場合 (窓が追従しない) も食い違いとして拾う。
  評価点: 列の x・壁の r の不一致は例外。自分の節点では補間値 = 節点値、評価点の最大・最小 = 節点の最大・最小。領域の境界 (x = −4・1、η = 0.1・0.9)。
  結合 (evaluate を合成の格子・主指標で): 支持の場合・(a) の幅で前提未達・(b) の RISING で前提未達・固定の条件の不成立 → 判別不能。
usage: python3 test_euler_t0_e2_eval.py
"""
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import euler_t0_e2_eval as EV  # noqa: E402

sys.path.insert(0, str(EV.TOOLS))
from stage_manifest import StageManifest  # noqa: E402

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


# --- 合成の主指標 ---------------------------------------------------------------------------------------------------------------
def ind_const(base=0.0, throat=None, spread=0.2):
    """全領域で max = base + spread、min = base − spread、q99 = |base| + spread/2。throat (スロートの 3 領域の max) を上書きできる。"""
    out = {}
    for r in EV.REGION_KEYS:
        mx, mn = base + spread, base - spread
        if throat is not None and r in EV.THROAT_REGIONS:
            mx = throat
        out[r] = {"n": 100, "max": mx, "min": mn, "q99": abs(base) + spread / 2,
                  **{s: 0.0 for s in EV.FRAC_STATS}}
    return out


def series(fn):
    """fn(step) → 主指標 の {step: 主指標} (本段の 19 枚)。"""
    return {s: fn(s) for s in EV.MAIN_STEPS}


# --- 判定 -------------------------------------------------------------------------------------------------------------------------
A300 = series(lambda s: ind_const(0.0, throat=300.0))
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=0.5))}, True)
check("支持: A スロート 300 K・B 全領域 ±0.5 K・前提成立", r["verdict"] == EV.LBL_SUPPORT)
check("支持の文言は観測区間に限定し「収束解」と言わない", "この観測区間" in EV.LBL_SUPPORT and "収束解が一致" not in EV.LBL_SUPPORT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=1.0))}, True)
check("境界: B の最大 = +1.0 K・最小 = −1.0 K は ±1 K 以内 (支持)", r["verdict"] == EV.LBL_SUPPORT)
W = EV.WIN13
check("評価窓は本段 42000〜54000 の 13 枚、末尾は 50000〜54000 の 5 枚 (本段 54000、55 枚)",
      W == tuple(range(42000, 54001, 1000)) and EV.TAIL5 == tuple(range(50000, 54001, 1000)) and len(EV.MAIN_STEPS) == 55
      and EV.MAIN_NSTEPS == 54000)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=1.0 + (0.01 if s == W[6] else 0.0)))}, True)
check(f"全保存時点: B が 13 枚のうち 1 枚 ({W[6]}) だけ 1.01 K → 判別不能", r["verdict"] == EV.LBL_UNDET)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=1.0 + (5.0 if s == W[0] - 1000 else 0.0)))}, True)
check(f"評価窓の直前 ({W[0] - 1000}) の B の 6 K は判定に入らない (支持)", r["verdict"] == EV.LBL_SUPPORT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=1.0 + (5.0 if s == 0 else 0.0)))}, True)
check("評価窓の外 (本段 0) の B の 6 K は判定に入らない (支持)", r["verdict"] == EV.LBL_SUPPORT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, spread=0.5))}, False, ["[P1] A: 不成立"])
check("前提未達なら支持の形でも判別不能 (理由を残す)", r["verdict"] == EV.LBL_UNDET and "[P1] A: 不成立" in r["reasons"]
      and r["facts"]["A_throat_anomaly_all_steps"] and r["facts"]["B_within_1K_all_steps"])
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, throat=300.5))}, True)
check("退ける: 両側にスロート 300 K の異常・差 0.5 K", r["verdict"] == EV.LBL_REJECT)
check("退ける の文言は収束解の依存を否定しない", "収束解に対する配点の依存までは否定しない" in EV.LBL_REJECT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, throat=301.0))}, True)
check("境界: 差 = 1.0 K は 1 K 以内 (退ける)", r["verdict"] == EV.LBL_REJECT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, throat=300.5 + (0.7 if s == W[-1] else 0.0)))}, True)
check(f"全保存時点: 差が最後の 1 枚 ({W[-1]}) だけ 1.2 K → 判別不能", r["verdict"] == EV.LBL_UNDET)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.3, throat=300.5))}, True)
check("全領域の差: スロートの max は同じでも他の領域の min・q99 の差 0.3 K は 1 K 以内 (退ける)", r["verdict"] == EV.LBL_REJECT)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(1.5, throat=300.5))}, True)
check("全領域の差: 他の領域の差 1.5 K → 判別不能", r["verdict"] == EV.LBL_UNDET)
r = EV.judge({"A": A300, "B": series(lambda s: ind_const(0.0, throat=50.0))}, True)
check("中間的な改善 (B のスロート 50 K) → 判別不能", r["verdict"] == EV.LBL_UNDET and any("中間" in x for x in r["reasons"]))
r = EV.judge({"A": series(lambda s: ind_const(0.0, throat=60.0)), "B": series(lambda s: ind_const(0.0, spread=0.5))}, True)
check("A がスロートで > 100 K を再現しない → 判別不能", r["verdict"] == EV.LBL_UNDET and any("再現しない" in x for x in r["reasons"]))
r = EV.judge({"A": series(lambda s: ind_const(0.0, throat=300.0 if s != W[3] else 99.0)), "B": series(lambda s: ind_const(0.0, spread=0.5))}, True)
check(f"A の異常が 13 枚のうち 1 枚 ({W[3]}) で 99 K → 判別不能", r["verdict"] == EV.LBL_UNDET)
ia2 = ind_const(0.0, throat=None)
ia2["thr_wall"]["max"] = 150.0                  # 3 領域のうち壁域だけが > 100 K
r = EV.judge({"A": series(lambda s: ia2), "B": series(lambda s: ind_const(0.0, spread=0.5))}, True)
check("スロートの 3 領域を合わせた max(|最大|, |最小|) > 100 K (壁域だけで超える) を異常と読む (各領域すべては要求しない; 支持)",
      r["verdict"] == EV.LBL_SUPPORT and EV.anomaly_K(ia2) == 150.0)
b = series(lambda s: ind_const(0.0, spread=0.5))
del b[W[9]]
r = EV.judge({"A": A300, "B": b}, True)
check(f"評価窓の 1 枚 (B の {W[9]}) が欠ける → 判別不能", r["verdict"] == EV.LBL_UNDET and any("欠ける" in x for x in r["reasons"]))
b = series(lambda s: ind_const(0.0, spread=0.5))
del b[54000]
r = EV.judge({"A": A300, "B": b}, True)
check("評価窓の最後 (B の 54000) が欠ける → 判別不能", r["verdict"] == EV.LBL_UNDET and any("B:54000" in x for x in r["reasons"]))
bn = series(lambda s: ind_const(0.0, spread=0.5))
bn[W[10]]["dn_core"]["q99"] = math.nan
r = EV.judge({"A": A300, "B": bn}, True)
check("評価窓の主指標に非有限 → 判別不能", r["verdict"] == EV.LBL_UNDET)
ia = ind_const(0.0, throat=None)
ia["thr_core"]["min"] = -250.0
r = EV.judge({"A": series(lambda s: ia), "B": series(lambda s: ind_const(0.0, spread=0.5))}, True)
check("負側の異常 (A のスロートの最小 −250 K) も > 100 K の異常と読む (支持)", r["verdict"] == EV.LBL_SUPPORT)


# --- 前提 P1: 本物の check_convergence の出力 -------------------------------------------------------------------------------------
def cfg(cfl, conv, nsteps=None, out=None):
    """合成の solverConfig (本段の step 数と出力間隔は既定で評価器の定数)。"""
    nsteps = EV.MAIN_NSTEPS if nsteps is None else nsteps
    out = EV.OUT_INTERVAL if out is None else out
    return ("solver: \"SLAU\"\nphysProp: {thermalMethod: 2}\ntime:\n  last: {nStepOuter: %d}\n  outStepInterval: %d\n"
            "  deltaT: {cfl: %s, cfl_pseudo: %s}\nspace: {convMethod: %d, limiter: 2}\nturbulence: {model: \"none\"}\n"
            % (nsteps, out, cfl, cfl, conv))
BC = "inlet: {physID: 1, kind: inlet_Pressure, floats: {Pt: 5500000.0, Tt: 1600.0}}\n"
RCOLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roY0", "rms_roY1"]


def resid(kind, n=300):
    i = np.arange(n)
    if kind == "pass":
        return 10.0 ** (-4.0 * i / (n - 1))
    if kind == "stall":
        return 1e-3 * (1.0 + 0.02 * np.sin(1.7 * i))
    if kind == "converging":
        return 1e-3 * 10.0 ** (-2.0 * i / (n - 1))
    if kind == "rising":
        v = 1e-3 * (1.0 + 0.02 * np.sin(1.7 * i))
        k = int(0.8 * n)
        v[k:] = 1e-3 * 10.0 ** (2.0 * (i[k:] - k) / (n - 1 - k))
        return v
    if kind == "nan":
        v = resid("stall", n)
        v[-3:] = np.nan
        return v
    raise ValueError(kind)


def make_run(root: Path, name: str, kinds: dict, joined=False, cfg_steps=None, resid_steps=None) -> Path:
    """soft (1 次) → main (2 次) の段の記録・残差履歴・本段の solverConfig.yaml を持つ合成の run。joined=True は soft も 2 次にして、
    check_convergence --segment が soft と main を 1 区間に連結する (区間が main だけでない) 場合を作る。cfg_steps: config・manifest の
    本段の nStepOuter (既定 MAIN_NSTEPS)。resid_steps: 本段の残差の step 数 (最終 step = resid_steps − 1; 既定 MAIN_NSTEPS)。"""
    rd = root / name
    rd.mkdir()
    cfg_steps = EV.MAIN_NSTEPS if cfg_steps is None else cfg_steps
    resid_steps = EV.MAIN_NSTEPS if resid_steps is None else resid_steps
    sm = StageManifest(rd)
    sm.add("soft", cfg("0.5", 1 if joined else 0, nsteps=EV.SOFT_STEPS, out=EV.SOFT_STEPS), BC, history="residual_history_soft.csv")
    main_cfg = cfg("2.0", 1, nsteps=cfg_steps)
    sm.add("main", main_cfg, BC, history="residual_history_main.csv")
    sm.write()
    (rd / "solverConfig.yaml").write_text(main_cfg)
    for tag, kk, ns in (("soft", {c: "pass" for c in RCOLS}, EV.SOFT_STEPS), ("main", kinds, resid_steps)):
        cols = {c: resid(kk.get(c, "stall")) for c in RCOLS}
        n = len(cols[RCOLS[0]])
        stp = np.round(np.linspace(0, ns - 1, n)).astype(int)
        with open(rd / f"residual_history_{tag}.csv", "w") as f:
            f.write("step," + ",".join(RCOLS) + "\n")
            for s in range(n):
                f.write(f"{stp[s]}," + ",".join(f"{cols[c][s]:.6e}" for c in RCOLS) + "\n")
    return rd


def conv_text(rd: Path) -> str:
    r = subprocess.run([sys.executable, str(EV.TOOLS / "check_convergence.py"), str(rd), "--segment"], capture_output=True, text=True)
    return r.stdout + r.stderr


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    cases = {
        "pass": ({c: "pass" for c in RCOLS}, "pass", True),
        "stall_all": ({c: "stall" for c in RCOLS}, "stall_only", True),
        "stall_one_pass_rest": ({**{c: "pass" for c in RCOLS}, "rms_roY1": "stall"}, "stall_only", True),
        "plateau_but_rising": ({**{c: "stall" for c in RCOLS}, "rms_roUy": "rising"}, "fail", False),
        "plateau_but_converging": ({**{c: "stall" for c in RCOLS}, "rms_roe": "converging"}, "fail", False),
        "plateau_but_nan": ({**{c: "stall" for c in RCOLS}, "rms_ro": "nan"}, "fail", False),
    }
    texts = {}
    for nm, (kinds, want, ok) in cases.items():
        t = conv_text(make_run(root, nm, kinds))
        texts[nm] = t
        p = EV.parse_convergence(t)
        check(f"P1 {nm}: status {p['status']} (期待 {want})、前提 {'成立' if ok else '不成立'}、最終 step {p.get('last_step')}",
              p["status"] == want and EV.convergence_ok(p) == ok and p["segment"] == "main" and p.get("last_step") == EV.MAIN_NSTEPS - 1)
    head_b = [l for l in texts["plateau_but_rising"].splitlines() if l.startswith("=== ")]
    check("(b) の前提の確認: check_convergence の総合表示は plateau (RISING の列があっても)",
          len(head_b) == 1 and "stalled/plateau" in head_b[0] and "RISING" in texts["plateau_but_rising"])
    check("(b) 総合は plateau だが 1 列が RISING → 不成立で、理由に RISING の列名",
          any("rms_roUy" in x and "RISING" in x for x in EV.parse_convergence(texts["plateau_but_rising"])["reasons"]))
    lab = EV.convergence_label(EV.parse_convergence(texts["stall_all"]))
    check("停滞のみの表示は PASS と書かない (「PASS ではない」と明記)", lab != "PASS" and "PASS ではない" in lab)
    t = conv_text(make_run(root, "joined", {c: "stall" for c in RCOLS}, joined=True))
    p = EV.parse_convergence(t)
    check("区間が main だけでない (soft → main を連結) → 判定不能", p["status"] == "undeterminable" and EV.convergence_ok(p) is False)
    t = texts["stall_all"].replace("  rms_roe     :", "  rms_roeX    :")
    check("必須の残差列 (rms_roe) の行が無い → 判定不能", EV.parse_convergence(t)["status"] == "undeterminable")
    check("判定ファイルなし → missing・不成立", EV.parse_convergence(None)["status"] == "missing" and not EV.convergence_ok(EV.parse_convergence(None)))
    t = texts["stall_all"].replace("[segment] 判定区間", "[segmentX] 区間")
    check("--segment の区間の行が無い → 判定不能", EV.parse_convergence(t)["status"] == "undeterminable")
    t = texts["pass"].replace("PASS (converged)", "NOT CONVERGED (stalled/plateau — x)")
    check("全体行と列の内訳の食い違い (全列合格なのに plateau) → 判定不能", EV.parse_convergence(t)["status"] == "undeterminable")


# --- 前提 P2: 本物の check_quasisteady と幅の条件 ------------------------------------------------------------------------------------
QC = EV.qs_columns()


def window_case(root: Path, name: str, steps, fn, want=None) -> dict:
    """fn(列, k (窓内の番号)) → 値。窓の CSV を書き、本物の check_quasisteady を回し、window_check の結果を返す (want = 登録の時点)。"""
    rows = [{"stage": "main", "step": s, **{c: fn(c, k) for c in QC}} for k, s in enumerate(steps)]
    p = root / f"{name}.csv"
    EV.write_series_csv(p, rows, QC)
    txt = EV.run_quasisteady(p, QC, root / f"{name}.txt")
    st, vals = EV.read_csv_cols(p, QC)
    return EV.window_check(st, vals, txt, QC, steps if want is None else want), txt


def const_fn(c, k):
    return 300.0 if c.startswith("max_thr") else (0.4 if c.startswith(("max", "q99")) else -0.4)


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    wc, _ = window_case(root, "const13", EV.WIN13, const_fn)
    check("P2: 一定の 27 列 (13 枚) は STEADY・幅 0 で成立", wc["ok"] and wc["n_bad"] == 0)
    ramp = lambda c, k: 300.0 + 4.0 * k / 12.0 if c == "max_thr_wall" else const_fn(c, k)  # noqa: E731
    wc, txt = window_case(root, "ramp13", EV.WIN13, ramp)
    col = wc["columns"]["max_thr_wall"]
    check("(a) 300→304 K の 13 点: check_quasisteady (--tail 1) の VERDICT は STEADY (併記)", col["verdict_quasisteady_recorded"] == "STEADY")
    check("(a) 同じ列の幅 4 K > 0.1 K で前提不成立 (窓全体も不成立)", (not col["ok"]) and abs(col["range_K"] - 4.0) < 1e-9
          and not wc["ok"] and wc["bad"] == ["max_thr_wall"])
    check("(a) --tail 1 で回している (出力に --tail 1)", "--tail 1" in txt)
    wc, _ = window_case(root, "ramp5", EV.TAIL5, lambda c, k: 300.0 + 4.0 * k / 4.0 if c == "max_thr_wall" else const_fn(c, k))
    check("(a) 末尾 5 枚の 300→304 K も STEADY だが幅で不成立",
          wc["columns"]["max_thr_wall"]["verdict_quasisteady_recorded"] == "STEADY" and not wc["columns"]["max_thr_wall"]["ok"])
    small = lambda c, k: 300.0 + 0.05 * k / 12.0 if c == "max_thr_wall" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "small13", EV.WIN13, small)
    check("P2: 幅 0.05 K の緩い上昇は幅 ≤ 0.1 K で成立", wc["ok"])
    edge = lambda c, k: (0.1 if k == 6 else 0.0) if c == "max_up_core" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "edge13", EV.WIN13, edge)
    check("P2: 幅 = 0.1 K ちょうど (0.0 と 0.1、float でも 0.1) は成立 (≤)", wc["columns"]["max_up_core"]["range_K"] == 0.1
          and wc["columns"]["max_up_core"]["ok"])
    edge3 = lambda c, k: 300.0 + (0.1 if k == 6 else 0.0) if c == "max_thr_wall" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "edge300", EV.WIN13, edge3)
    check("P2: 比較は float の ≤ そのまま (許容を足さない): 300.1 − 300 = 0.1 + 2.3e-14 は不成立",
          wc["columns"]["max_thr_wall"]["range_K"] > 0.1 and not wc["columns"]["max_thr_wall"]["ok"])
    drift = lambda c, k: 300.0 * (1.0 + 0.2 * k / 12.0) if c == "max_thr_wall" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "drift13", EV.WIN13, drift)
    check("P2: 相対 20 % のドリフト (幅 60 K) → 不成立 (VERDICT は STEADY 以外を併記)",
          wc["columns"]["max_thr_wall"]["verdict_quasisteady_recorded"] != "STEADY" and not wc["ok"])
    nz = lambda c, k: 0.07 + 0.001 * k + 0.005 * math.sin(2.3 * k) if c == "max_dn_core" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "nearzero13", EV.WIN13, nz)
    col = wc["columns"]["max_dn_core"]
    check(f"P2 (a59e61b5): 近零の偏差 (平均 0.07 K) は VERDICT {col['verdict_quasisteady_recorded']} (STEADY でない) でも幅 "
          f"{col['range_K']:.3f} K ≤ 0.1 K なら成立", col["verdict_quasisteady_recorded"] != "STEADY" and col["ok"] and wc["ok"])
    nz5 = lambda c, k: 0.01 + 0.004 * k if c == "q99_up_axis" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "nearzero5", EV.TAIL5, nz5)
    col = wc["columns"]["q99_up_axis"]
    check(f"P2: 末尾 5 枚の近零の単調な上昇 (0.010→0.026 K、VERDICT {col['verdict_quasisteady_recorded']}) も幅で成立",
          col["verdict_quasisteady_recorded"] != "STEADY" and col["ok"])
    nzbig = lambda c, k: 0.02 + 0.2 * k / 12.0 if c == "min_thr_axis" else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "nearzero_wide", EV.WIN13, nzbig)
    check("P2: 近零でも幅 0.2 K の上昇は不成立", not wc["columns"]["min_thr_axis"]["ok"] and not wc["ok"])
    wc, _ = window_case(root, "short", EV.WIN13[:-1], const_fn, want=EV.WIN13)
    check("P2: 54000 を欠く 13 枚の窓 (12 枚) → 不成立", (not wc["ok"]) and wc["problems"])
    wc, _ = window_case(root, "short5", EV.TAIL5[:-1], const_fn, want=EV.TAIL5)
    check("P2: 54000 を欠く末尾 5 枚の窓 (4 枚) → 不成立", (not wc["ok"]) and wc["problems"])
    wc, _ = window_case(root, "old18000", tuple(range(6000, 18001, 1000)), const_fn, want=EV.WIN13)
    check("P2: 旧登録の窓 (6000〜18000) の時点は登録の窓と違い不成立", (not wc["ok"]) and wc["problems"])
    nanf = lambda c, k: math.nan if (c == "q99_dn_core" and k == 3) else const_fn(c, k)  # noqa: E731
    wc, _ = window_case(root, "nan13", EV.WIN13, nanf)
    check("P2: 非有限の値 → 不成立", not wc["columns"]["q99_dn_core"]["ok"] and not wc["ok"])


# --- 評価点と領域 -------------------------------------------------------------------------------------------------------------------
x = np.array([-5.0, -4.0, 0.0, 1.0, 2.0])
rw = np.array([3.0, 2.0, 1.0, 1.5, 2.5])
sA = np.array([0.0, 0.3, 0.7, 0.95, 0.999, 1.0])
sB = np.array([0.0, 0.2, 0.5, 0.8, 0.95, 1.0])
XA = np.repeat(x[:, None], sA.size, axis=1)
RA_, RB_ = rw[:, None] * sA[None, :], rw[:, None] * sB[None, :]
ep = EV.eval_points(XA, RA_, XA.copy(), RB_)
check("評価点: 列ごとに A と B の η の和集合 (重複は 1 つ)", ep["eta"].size == x.size * np.unique(np.r_[sA, sB]).size)
devA = np.sin(3.0 * RA_) + XA
vals = EV.to_points(devA, EV.column_eta(RA_), ep)
own = np.isin(ep["eta"], sA)
check("評価点: 自分の節点では補間値 = 節点値", np.allclose(vals[own], (np.sin(3.0 * rw[ep["col"]] * ep["eta"]) + ep["x_rt"])[own], atol=1e-12))
check("評価点: 評価点の最大・最小 = 節点の最大・最小", vals.max() == devA.max() and vals.min() == devA.min())
try:
    EV.eval_points(XA, RA_, XA + 1e-9, RB_)
    ok_ = False
except ValueError:
    ok_ = True
check("評価点: 列の x が A と B で違えば例外", ok_)
try:
    RB2 = RB_.copy()
    RB2[2, -1] *= 1.0 + 1e-7
    EV.eval_points(XA, RA_, XA, RB2)
    ok_ = False
except ValueError:
    ok_ = True
check("評価点: 壁の r が A と B で違えば例外", ok_)
check("領域の境界: x = −4・1 は thr、−4.0001 は up、1.0001 は dn",
      list(EV.x_band([-4.0001, -4.0, 1.0, 1.0001])) == [0, 1, 1, 2])
check("領域の境界: η = 0.1・0.9 は core、0.0999 は axis、0.9001 は wall",
      list(EV.eta_band([0.0999, 0.1, 0.9, 0.9001])) == [0, 1, 1, 2])
iv = EV.indicators(np.array([150.0, -2.0, 0.5, 0.0]), np.array([4, 4, 4, 4]))
k = EV.REGION_KEYS[4]
check(f"主指標: 領域 {k} の最大・最小・割合 (+1 K 超 1/4、−1 K 未満 1/4、|·| > 100 K 1/4)",
      k == "thr_core" and iv[k]["max"] == 150.0 and iv[k]["min"] == -2.0 and iv[k]["fpos1"] == 0.25 and iv[k]["fneg1"] == 0.25
      and iv[k]["fabs100"] == 0.25 and iv["up_axis"]["n"] == 0 and iv["up_axis"]["max"] is None)


# --- 結合 (evaluate を合成の格子・主指標で; 復元・h5 は差し替え) ---------------------------------------------------------------------
def run_evaluate(ind_by_arm: dict, conv_kinds: dict, fixed_ok=True, run_kw=None):
    with tempfile.TemporaryDirectory() as td:
        case = Path(td)
        for arm, rn in EV.RUNS.items():
            rd = make_run(case, rn, conv_kinds[arm], **((run_kw or {}).get(arm, {})))
            (rd / EV.SEGMENT_VERDICT_FILE).write_text(conv_text(rd))
        geo = (XA, RA_, 0.0768)
        geoB = (XA, RB_, 0.0768)
        sv = (EV.load_geometry, EV.eval_snapshot, EV.fixed_conditions)
        EV.load_geometry = lambda rd: geo if rd.name == EV.RUNS["A"] else geoB
        arm_of = {v: k for k, v in EV.RUNS.items()}

        def fake_snap(rd, h5, X, R, eta, ep_):
            arm = arm_of[rd.name]
            step = int(h5.stem.split("_")[1])
            if h5.parent.name == EV.SOFT_DIR:
                return {"file": str(h5), "problems": [], "indicators": ind_const(0.0, spread=50.0)}
            return {"file": str(h5), "problems": [], "indicators": ind_by_arm[arm][step]}
        EV.eval_snapshot = fake_snap
        EV.fixed_conditions = lambda case_, ok, why: {"ok": fixed_ok and ok, "problems": ([] if (fixed_ok and ok) else ["合成: 固定の条件の不成立"])}
        try:
            out = EV.evaluate(case)
        finally:
            EV.load_geometry, EV.eval_snapshot, EV.fixed_conditions = sv
        csv_ok = all((case / rn / f"euler_t0_e2_series_{w}.csv").is_file() for rn in EV.RUNS.values() for w in EV.WINDOWS)
        return out, csv_ok


STALL = {c: "stall" for c in RCOLS}
B05 = series(lambda s: ind_const(0.0, spread=0.5))
out, csv_ok = run_evaluate({"A": A300, "B": B05}, {"A": STALL, "B": STALL})
check("結合: 両側停滞のみ・主指標一定 → 前提成立・支持", out["precondition"]["ok"] and out["VERDICT"] == EV.LBL_SUPPORT and csv_ok)
check("結合: 収束の表示は停滞のみ (PASS と書かない)", all("PASS ではない" in v for v in out["precondition"]["convergence_labels"].values()))
check("結合: 出力に評価器の sha256・登録の commit a59e61b5・前提の定義・step 数", out["evaluator_sha256"] == EV._sha(Path(EV.__file__))
      and out["plan_reg_commit"] == "a59e61b5" and out["precondition_definition"]["range_tol_K"] == 0.1
      and out["precondition_definition"]["qs_tail"] == 1.0 and "定常化・収束の認定ではない" in out["precondition_definition"]["p2_note"]
      and out["steps"]["evaluator"]["MAIN_NSTEPS"] == 54000
      and all(v["effective"]["residual_last_step"] == 53999 and v["effective"]["nStepOuter"] == 54000
              and v["effective"]["manifest_main_nStepOuter"] == 54000 and not v["problems"] for v in out["steps"]["runs"].values()))
Bnz = series(lambda s: ind_const(0.07 + 0.001 * (s / 1000.0) + 0.005 * math.sin(2.3 * s / 1000.0), spread=0.3))
out, _ = run_evaluate({"A": A300, "B": Bnz}, {"A": STALL, "B": STALL})
qb = out["per_run"]["B"]["quasisteady"]["win13"]
check("結合 (a59e61b5): B の近零の主指標が STEADY でなくても (VERDICT 併記) 幅 ≤ 0.1 K なら前提成立・支持",
      out["precondition"]["ok"] and out["VERDICT"] == EV.LBL_SUPPORT
      and any(v["verdict_quasisteady_recorded"] != "STEADY" for v in qb["columns"].values()) and qb["ok"])
Aramp = series(lambda s: ind_const(0.0, throat=300.0 + (4.0 * (s - W[0]) / 12000.0 if s >= W[0] else 0.0)))
out, _ = run_evaluate({"A": Aramp, "B": B05}, {"A": STALL, "B": STALL})
check("結合 (a): A のスロートの最大が 300→304 K (STEADY) → 幅で前提未達・判別不能",
      (not out["precondition"]["ok"]) and out["VERDICT"] == EV.LBL_UNDET
      and out["per_run"]["A"]["quasisteady"]["win13"]["columns"]["max_thr_wall"]["verdict_quasisteady_recorded"] == "STEADY"
      and not out["per_run"]["A"]["quasisteady"]["win13"]["columns"]["max_thr_wall"]["ok"])
out, _ = run_evaluate({"A": A300, "B": B05}, {"A": STALL, "B": {**STALL, "rms_roUy": "rising"}})
check("結合 (b): B の総合は plateau だが 1 列が RISING → 前提未達・判別不能",
      (not out["precondition"]["ok"]) and out["VERDICT"] == EV.LBL_UNDET and out["per_run"]["B"]["convergence"]["status"] == "fail")
out, _ = run_evaluate({"A": A300, "B": B05}, {"A": STALL, "B": STALL}, fixed_ok=False)
check("結合: 固定の条件の不成立 → 判別不能", out["VERDICT"] == EV.LBL_UNDET and any(x.startswith("[固定]") for x in out["precondition"]["reasons"]))
out, _ = run_evaluate({"A": A300, "B": series(lambda s: ind_const(0.0, throat=300.5))}, {"A": {c: "pass" for c in RCOLS}, "B": STALL})
check("結合: A は PASS・B は停滞のみ、両側の異常・差 0.5 K → 退ける", out["VERDICT"] == EV.LBL_REJECT)

out, _ = run_evaluate({"A": A300, "B": B05}, {"A": STALL, "B": STALL}, run_kw={"B": {"cfg_steps": 18000}})
check("結合 (同期): B の本段の config・manifest が 18000 → [固定] の不成立・判別不能",
      out["VERDICT"] == EV.LBL_UNDET and any(x.startswith("[固定] B:") and "nStepOuter" in x for x in out["precondition"]["reasons"]))
out, _ = run_evaluate({"A": A300, "B": B05}, {"A": STALL, "B": STALL}, run_kw={"A": {"resid_steps": 18000}})
check("結合 (同期): A の本段の残差が 17999 で終わる (config は 54000) → 判別不能",
      out["VERDICT"] == EV.LBL_UNDET and any("最終 step 17999" in x for x in out["precondition"]["reasons"]))


# --- step 数の同期 (評価器と実行側) ---------------------------------------------------------------------------------------------
ok_eff = {"nStepOuter": 54000, "outStepInterval": 1000, "manifest_main_nStepOuter": 54000, "residual_last_step": 53999}
check("同期: 評価器の窓は MAIN_NSTEPS から作られた形 (食い違いなし)", EV.window_problems() == [])
check("同期: 実効値 54000・1000・manifest 54000・最終 step 53999 は食い違いなし", EV.step_problems(ok_eff) == [])
for k_, v_ in (("nStepOuter", 18000), ("outStepInterval", 500), ("manifest_main_nStepOuter", 18000), ("residual_last_step", 17999)):
    check(f"同期: {k_} = {v_} は食い違い", len(EV.step_problems({**ok_eff, k_: v_})) == 1)
check("同期: config の読み (yaml_strict)", EV.config_steps(cfg("2.0", 1)) == {"nStepOuter": 54000, "outStepInterval": 1000})
try:
    EV.config_steps(cfg("2.0", 1) + "time:\n  last: {nStepOuter: 3}\n")
    ok_ = False
except Exception:  # noqa: BLE001
    ok_ = True
check("同期: config の重複キーは例外 (solver の先勝ちと食い違う値を照合しない)", ok_)
sv_n = EV.MAIN_NSTEPS
EV.MAIN_NSTEPS = 60000              # 評価器の定数だけを書き換え、窓を作り直さなかった場合
try:
    wp = EV.window_problems()
    sp = EV.step_problems({**ok_eff, "nStepOuter": 60000, "manifest_main_nStepOuter": 60000, "residual_last_step": 59999})
finally:
    EV.MAIN_NSTEPS = sv_n
check("同期: 評価器の MAIN_NSTEPS だけを変えて窓が追従しない → 食い違い (実行側も同じ関数で止まる)", bool(wp) and bool(sp))
import euler_t0_e2 as RUN  # noqa: E402
with tempfile.TemporaryDirectory() as td:
    rd = Path(td) / EV.RUNS["A"]
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(cfg("2.0", 1, nsteps=18000))
    (rd / "prepare_info.json").write_text(
        '{"e2": {"arm": "A"}, "ic": {"VERDICT": "OK"}, "moc": {"gate": {"pass": true}}, "mesh": {"ni": 2000}}')
    cp = RUN.check_run_steps(rd)
    check("実行側: run の config が 18000 なら check_run_steps が食い違いを返す", bool(cp) and "18000" in cp[0])
    try:
        RUN.run(rd)
        ok_, msg = False, "起動した"
    except SystemExit as e:
        ok_, msg = "同期していない" in str(e), str(e)
    check("実行側: run は step 数の食い違いで forge の前に止まる (SystemExit)", ok_)
    (rd / "solverConfig.yaml").write_text(cfg("2.0", 1))
    check("実行側: config が 54000・1000 なら食い違いなし", RUN.check_run_steps(rd) == [])
    (rd / "solverConfig.yaml").write_text(cfg("2.0", 1, out=500))
    check("実行側: outStepInterval 500 は食い違い", bool(RUN.check_run_steps(rd)))

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
