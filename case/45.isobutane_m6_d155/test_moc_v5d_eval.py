"""moc_v5d_eval.py と moc_v5d.py (問題の生成と検査・run の前の停止・起動の順) の試験 (plan discretization-moc-axis-limit-and-corrector §6 V5d、
登録 99431498)。forge・AWS は使わない (V5d の run の結果は読まない。合成した評価量の系列・全温の主指標・M_common・残差 CSV だけ)。
check_convergence.py と check_quasisteady.py は本物を呼ぶ (出力の書式を合成で真似ない)。保存場を読む部分 (格子・全温の復元・M_common・
壁の照合・評価座標の感度) は差し替える。
確かめること:
  定数: 窓 42000〜54000 (13 枚)・50000〜54000 (5 枚) が E4 と同じ / 腕の構成 (既定 6 本、共用で run_0164 を腕 B に・run_0173 を作らない) /
    起動の順は腕を交互に / --md-offset の正規化と範囲。
  部品: exit_M_dev の例外 (V5 と同じ条件) / 幅の条件 (Δq/10・符号付き出口 M 1e-5、境界は float の ≤) / 準定常の判定 (両窓 STEADY、
    exit_M_dev だけ絶対許容内) / 出口較正の 3 区分 / M = 6 の判定 (達成・未達・判別不能)。
  結合 (evaluate、6 本の合成 run): 陽性の対照 (許容幅内・据え置き・達成) / 悪化 / 保留 (ユーザ判断) / 前提の不成立で保留:
    残差に RISING・全温 ±1.5 K・全温の幅 0.15 K・腕 B だけ非 STEADY (対称)・幅だけ超える (STEADY でも)・符号付き出口 M の幅 1.2e-5
    (例外の範囲内でも)・forge の sha256・Md_moc_offset の食い違い・E4 の結果が採用でない・E4 の値の食い違い・起動の記録の食い違い・
    X_F のずれ・古い CSV・旧 G1 の格子 (mesh)・腕 M のゲート不合格・壁の取り違え・腕の中の壁の不一致・prep の設定ファイルの不一致・
    RUN_RC・EARLY_STOP・乾式 prep・run の欠損 / 出口較正のやり直し / 腕 M の M_common が 1e-4 を外れる → 未達 (NS を始めない) /
    腕 M の M_common の幅 6e-5 → 判別不能 / 共用 (run_0164 を腕 B に): 陽性・E4 の参照 run の食い違い・共用しない構成に run_0164。
  moc_v5d: 問題の生成 (B = E4 の d0 の name・Md_moc_offset だけ、M = B + MOC の 2 キー、M = V5 の腕 M の問題の Md_moc_offset だけ違う) と改変の拒否、
    E4 の採用値の要求 (本番)・乾式は記録だけ、run dir があれば別の値で書き直さない、共用の照合、run は forge の前に止まる (乾式 prep・
    run 名・step 数・旧格子)、verify-prep、共用の照合 (腕 B の prep と run_0164: 初期線の絶対パスは比べない、Md_moc_offset・格子のハッシュ・
    設定ファイルの違いで不成立、run_0164 が無ければ本番は不成立・乾式は未確認)、共用の前提 (完走・窓の res・残差判定・較正値)。
usage: python3 test_moc_v5d_eval.py
"""
import os

os.environ["FORGE_BIN"] = "/nonexistent/forge-test-guard"     # 念のため (この試験は forge を起動しない)
os.environ["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"

import copy  # noqa: E402
import csv  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import yaml  # noqa: E402
from scipy.interpolate import make_interp_spline  # noqa: E402

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import moc_v5d_eval as EV  # noqa: E402
import euler_t0_e2_eval as EV2  # noqa: E402
import e4_recal_eval as EV4  # noqa: E402

sys.path.insert(0, str(EV.TOOLS))
from stage_manifest import StageManifest  # noqa: E402

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


W = EV.WIN13
MD = 6.8825162455159465e-06
MDR = repr(MD)

# --- 定数 ---------------------------------------------------------------------------------------------------------------------
check("窓: 42000〜54000 の 13 枚・50000〜54000 の 5 枚、E4 の定数と同じ (window_problems 空)",
      W == tuple(range(42000, 54001, 1000)) and EV.TAIL5 == tuple(range(50000, 54001, 1000)) and EV.window_problems() == []
      and EV.MAIN_NSTEPS == EV4.MAIN_NSTEPS == 54000 and EV.SOFT_STEPS == 3000)
check("量と Δq: V5 の 6 量 (0.001 / 0.010 / 0.003 / 0.003 / 0.03 / 1.8e-4)",
      EV.DQ == {"M_wave_eta0.1": 0.001, "P_wave_eta0.1": 0.010, "overshoot_eta0.1": 0.003, "overshoot_exitnorm_eta0.1": 0.003,
                "P_slope_abs_eta0.1": 0.03, "exit_M_dev": 0.00018})
a0, a1 = EV.arms(False), EV.arms(True)
check("腕: 既定は run_0171〜0173 (B)・run_0174〜0176 (M)、共用は B = run_0171・run_0172・run_0164 (run_0173 は作らない)",
      a0["B"] == ["run_0171_euler_v5d_B_r1", "run_0172_euler_v5d_B_r2", "run_0173_euler_v5d_B_r3"]
      and a0["M"] == ["run_0174_euler_v5d_M_r1", "run_0175_euler_v5d_M_r2", "run_0176_euler_v5d_M_r3"]
      and a1["B"] == ["run_0171_euler_v5d_B_r1", "run_0172_euler_v5d_B_r2", "run_0164_euler_e4_recal_d1"] and a1["M"] == a0["M"]
      and EV.new_runs(True)["B"] == a1["B"][:2] and EV.GEOM_REF == "run_0171_euler_v5d_B_r1")
check("--md-offset: repr に正規化 (6.8825162455159465e-6 → …e-06)", EV.parse_md("6.8825162455159465e-6") == (MD, MDR))
for bad_md in ("nan", "inf", "0.02", "abc", "-0.0101"):
    try:
        EV.parse_md(bad_md)
        ok_ = False
    except ValueError:
        ok_ = True
    check(f"--md-offset {bad_md!r} は拒否", ok_)

# --- 部品 ---------------------------------------------------------------------------------------------------------------------
sig = lambda amp, ph=0.0: np.array([amp * math.sin(1.3 * k + ph) for k in range(13)])  # noqa: E731
r = EV.exit_exception(-1.0e-5 + sig(2e-6))
check(f"exit_M_dev の例外: 符号付き −1e-5 ± 2e-6 → 絶対許容内 ({r.get('reason')})", r["applies"] and r["label"] == EV.LBL_ABS)
r = EV.exit_exception(-1.7e-5 + sig(2e-6))
check("exit_M_dev の例外: 最大絶対値 1.9e-5 > 1.8e-5 → 当たらない", not r["applies"])
r = EV.exit_exception(np.r_[np.full(12, -1e-5), np.nan])
check("exit_M_dev の例外: 非有限 → 当たらない", not r["applies"])

base_vals = {"M_wave_eta0.1": 0.0064 + sig(2e-5), "P_wave_eta0.1": 0.041 + sig(2e-4), "overshoot_eta0.1": 0.038 + sig(6e-5),
             "overshoot_exitnorm_eta0.1": 0.038 + sig(6e-5), "P_slope_abs_eta0.1": 0.19 + sig(6e-4),
             "exit_core_M": 5.99999 + sig(2e-6)}
base_vals["exit_M_dev"] = np.abs(base_vals["exit_core_M"] - 6.0)
bad, rec = EV.width_checks(base_vals)
check("幅: 合成の基準は全量 Δq/10 以内・符号付き出口 M 1e-5 以内", bad == [] and all(v["ok"] for v in rec.values()))
v2 = dict(base_vals)
v2["M_wave_eta0.1"] = np.r_[np.full(12, 0.0064), 0.0064 + 1e-4]
bad, rec = EV.width_checks(v2)
check(f"幅: M 波の幅 = Δq/10 (float {rec['M_wave_eta0.1']['width']!r}) は ≤ の比較そのまま",
      (rec["M_wave_eta0.1"]["ok"]) == (rec["M_wave_eta0.1"]["width"] <= 1e-4))
v2["M_wave_eta0.1"] = np.r_[np.full(12, 0.0064), 0.0064 + 1.05e-4]
bad, rec = EV.width_checks(v2)
check("幅: M 波の幅 1.05e-4 > Δq/10 → 不成立", any("M_wave" in b for b in bad))
v2 = dict(base_vals)
v2["exit_core_M"] = 5.99999 + np.r_[np.zeros(12), 1.2e-5]
bad, rec = EV.width_checks(v2)
check("幅: 符号付き出口 M の幅 1.2e-5 > 1e-5 → 不成立", any("符号付き" in b for b in bad))
st_ = {c: "STEADY" for c in EV.QS_COLS}
why, status = EV.judge_quasisteady({"win13": dict(st_), "tail5": dict(st_)}, {"applies": False})
check("準定常: 両窓 STEADY → 全量 STEADY", why == [] and all(v == "STEADY" for v in status.values()))
why, status = EV.judge_quasisteady({"win13": dict(st_, exit_M_dev="OSCILLATING"), "tail5": dict(st_)}, {"applies": True})
check("準定常: exit_M_dev だけ OSCILLATING で例外が成り立つ → 絶対許容内 (STEADY と表示しない)", why == [] and status["exit_M_dev"] == EV.LBL_ABS)
why, status = EV.judge_quasisteady({"win13": dict(st_), "tail5": dict(st_, **{"overshoot_eta0.1": "DRIFTING"})}, {"applies": True})
check("準定常: 他の量は末尾 5 枚だけ DRIFTING でも保留 (例外は exit_M_dev だけ)", why != [] and status["overshoot_eta0.1"] == EV.LBL_HOLD)
why, status = EV.judge_quasisteady({"win13": dict(st_, exit_core_M="OSCILLATING"), "tail5": dict(st_)}, {"applies": True})
check("準定常: 出口コア M (記録用の量) も STEADY を要求", status["exit_core_M"] == EV.LBL_HOLD)
check("出口較正: |D| + 2SE ≤ 1e-4 → 据え置き / |D| − 2SE > 1e-4 → やり直す / 間 → 保留",
      EV.judge_exit_calibration(4e-5, 2e-5) == EV.LBL_CAL_KEEP and EV.judge_exit_calibration(2e-4, 2e-5) == EV.LBL_CAL_REDO
      and EV.judge_exit_calibration(9e-5, 2e-5) == EV.LBL_CAL_HOLD and "M = 6 の達成の証明ではない" in EV.LBL_CAL_KEEP)
ser_ok = {s: 6.00002 for s in W}
check("M = 6: 13 枚すべて +2e-5 → within", EV.m6_target(ser_ok)["within"])
check("M = 6: 境界 6.0001 (float で |·| ≤ 1e-4) は内", EV.m6_target({s: 6.0001 for s in W})["within"])
check("M = 6: 1 枚 6.00010001 → 外", not EV.m6_target({**ser_ok, W[4]: 6.00010001})["within"])
check("M = 6: 1 枚欠ける → 不完全", not EV.m6_target({s: 6.0 for s in W[1:]})["complete"])
runsM = a0["M"]
p3ok = {r_: {"ok": True} for r_ in runsM}
mm = {r_: EV.m6_target(ser_ok) for r_ in runsM}
check("M = 6 の判定: 3 本とも内・前提成立 → 達成", EV.judge_m6(runsM, mm, p3ok, True, [])["verdict"] == EV.LBL_M6_OK)
mm2 = dict(mm, **{runsM[2]: EV.m6_target({s: 6.00012 for s in W})})
j = EV.judge_m6(runsM, mm2, p3ok, True, [])
check("M = 6 の判定: 1 本が外 → 未達 (NS を始めない・共通の較正値で両腕を確かめ直す)", j["verdict"] == EV.LBL_M6_NG and "NS を始めない" in j["verdict"])
check("M = 6 の判定: E4 の P3 が 1 本不成立 → 判別不能", EV.judge_m6(runsM, mm2, dict(p3ok, **{runsM[0]: {"ok": False}}), True, [])["verdict"] == EV.LBL_M6_UNDET)
check("M = 6 の判定: 腕 M の前提の不成立 → 判別不能", EV.judge_m6(runsM, mm, p3ok, False, ["x"])["verdict"] == EV.LBL_M6_UNDET)

# --- 結合 (evaluate を 6 本の合成 run で) ------------------------------------------------------------------------------------
RCOLS = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roY0", "rms_roY1"]
BC = "inlet: {physID: 1, kind: inlet_Pressure, floats: {Pt: 5500000.0, Tt: 1600.0}}\n"
FSHA = "6b47811bad9daa6807fb3f1e0808a1472c8791b572c2b7c6053fd9e225fa57d8"
PROV = f"forge_sha256: {FSHA}\n'slauWallNormalChi' effective: 1 (auto)\n'scalarGradient' effective: lsq (default)\n"
X_E0, X_F0 = 39.82004263, 95.22667765
GEOM = {"X_E": X_E0, "X_F": X_F0, "WIN_T": [X_E0 + 2.0, X_F0 - 1.0], "WIN_O": [X_E0 - 15.0, X_F0]}
SER_COLS = ["step", "M_wave_eta0.1", "P_wave_eta0.1", "overshoot_eta0.1", "overshoot_exitnorm_eta0.1", "P_slope_abs_eta0.1",
            "exit_M_dev", "exit_core_M", "M_wave_eta0.0"]
BASEQ = {"M_wave_eta0.1": 0.0064, "P_wave_eta0.1": 0.041, "overshoot_eta0.1": 0.038, "overshoot_exitnorm_eta0.1": 0.038,
         "P_slope_abs_eta0.1": 0.19, "exit_core_M": 5.99999, "M_wave_eta0.0": 0.02}
AMP = {"M_wave_eta0.1": 2e-5, "P_wave_eta0.1": 2e-4, "overshoot_eta0.1": 6e-5, "overshoot_exitnorm_eta0.1": 6e-5,
       "P_slope_abs_eta0.1": 6e-4, "exit_core_M": 2e-6, "M_wave_eta0.0": 1e-4}


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


def write_resid(rd: Path, kinds=None):
    for tag, ns, kk in (("soft", EV.SOFT_STEPS, {}), ("main", EV.MAIN_NSTEPS, kinds or {})):
        cols = {c: resid(kk.get(c, "stall" if tag == "main" else "pass")) for c in RCOLS}
        n = len(cols[RCOLS[0]])
        stp = np.round(np.linspace(0, ns - 1, n)).astype(int)
        with open(rd / f"residual_history_{tag}.csv", "w") as f:
            f.write("step," + ",".join(RCOLS) + "\n")
            for k in range(n):
                f.write(f"{stp[k]}," + ",".join(f"{cols[c][k]:.6e}" for c in RCOLS) + "\n")
    r_ = subprocess.run([sys.executable, str(EV.TOOLS / "check_convergence.py"), str(rd), "--segment"], capture_output=True, text=True)
    (rd / EV4.SEGMENT_VERDICT_FILE).write_text(r_.stdout + r_.stderr)


def spline_dict(shift=0.0):
    x = np.linspace(-1.0, 3.0, 801)
    s_ = make_interp_spline(x, 1 + x ** 2 / 4 - x ** 3 / 30 + shift, k=5)
    return {"t": [float(v) for v in s_.t], "c": [float(v) for v in s_.c], "k": 5}


def write_series(rd: Path, k: int, offset=None, fn=None):
    """評価量の系列 (本段 1000〜54000、54 枚)。offset: {量: 加算}、fn: {量: step → 値} (上書き)。exit_M_dev は |exit_core_M − 6| で作る。"""
    offset, fn = offset or {}, fn or {}
    with open(rd / EV.SERIES_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(SER_COLS)
        for st in range(1000, EV.MAIN_NSTEPS + 1, 1000):
            row = {c: BASEQ[c] + offset.get(c, 0.0) + AMP[c] * math.sin(1.3 * st / 1000.0 + 0.7 * k) for c in BASEQ}
            row.update({c: g(st) for c, g in fn.items()})
            row["exit_M_dev"] = abs(row["exit_core_M"] - 6.0)
            w.writerow([st] + [f"{row[c]:.10g}" for c in SER_COLS[1:]])
    (rd / EV.SERIES_QS).write_text("=== series -> STEADY ===\n" + "".join(f"  {c:16s}: tail mean=0  drift=0.0%/tail  fluct=0.1%   STEADY\n"
                                                                         for c in SER_COLS[1:]))


def write_run(case: Path, run: str, arm: str, k: int, share: bool):
    rd = case / run
    rd.mkdir()
    sm = StageManifest(rd)
    sm.add("soft", cfg("0.5", 0, EV.SOFT_STEPS, EV.SOFT_STEPS), BC, history="residual_history_soft.csv")
    main_cfg = cfg("2.0", 1, EV.MAIN_NSTEPS, EV.OUT_INTERVAL)
    sm.add("main", main_cfg, BC, history="residual_history_main.csv")
    sm.write()
    (rd / "solverConfig.yaml").write_text(main_cfg)
    (rd / "bcondConfig.yaml").write_text(BC)
    write_resid(rd)
    (rd / EV4.SOFT_DIR).mkdir()
    (rd / "RUN_PROVENANCE.txt").write_text(PROV)
    (rd / EV4.SOFT_DIR / "RUN_PROVENANCE.txt").write_text(PROV)
    (rd / "RUN_RC").write_text("0\n")
    moc = {"axis_limit": "legacy", "corrector": "fixed2", "gate": {"applicable": False, "pass": None}} if arm == "B" else \
          {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
    info = {"moc": moc, "wall_fit": {"mono_r2": [0.0, 1.5], "spline": spline_dict(0.0 if arm == "B" else 1e-5)}, "Md_moc_offset": MD,
            "scale_m": EV4.SCALE_M, "x_E": X_E0, "initial_line": {"run": f"/x/{EV4.INITIAL_LINE[0]}", "res": [EV4.INITIAL_LINE[1]]},
            "mesh": {"ni": 4, "nj": 5, "source": "mesh_euler", "params": {**EV4.MESH_EULER_EXPECT, "scale": EV4.SCALE_M},
                     "hashes": {"coord_sha256": f"coord{arm}", "topology_sha256": "topo"}},
            "ic": {"mode": "isentropic", "VERDICT": "OK"}}
    files = {n: f"same-{n}" for n in EV.SAME_FILES}
    if run == EV.E4V_RUN:
        info["e4"] = {"stage": "d1", "run": run}
        (rd / EV4.PREP_JSON).write_text(json.dumps({"stage": "d1", "run": run, "dry": False, "ic_check": {"VERDICT": "OK"},
                                                    "Md_moc_offset": MD, "files_sha256": files}))
    else:
        info["v5d"] = {"arm": arm, "md_offset_repr": MDR, "share_e4v": share}
        (rd / EV.PREP_JSON).write_text(json.dumps({"arm": arm, "runs": EV.new_runs(share)[arm], "md_offset_repr": MDR, "share_e4v": share,
                                                   "dry": False, "ic_check": {"VERDICT": "OK"},
                                                   "same_as_e4v": ({"status": "OK"} if share and arm == "B" else None), "files_sha256": files}))
    (rd / "prepare_info.json").write_text(json.dumps(info))
    write_series(rd, k)
    return rd


def refresh_series_record(case: Path):
    p = case / EV.SERIES_JSON
    r_ = json.loads(p.read_text())
    for run in r_["runs"]:
        if (case / run).is_dir():
            r_["runs"][run]["csv_sha256"] = sha(case / run / EV.SERIES_CSV)
            r_["runs"][run]["quasisteady_sha256"] = sha(case / run / EV.SERIES_QS)
    p.write_text(json.dumps(r_))


def build_base(root: Path, share: bool) -> Path:
    case = root / ("base_share" if share else "base")
    case.mkdir()
    shutil.copy(C / EV4.NS_REF_PROBLEM, case / EV4.NS_REF_PROBLEM)
    (case / EV4.REF_RUN_BIN).mkdir()
    (case / EV4.REF_RUN_BIN / "RUN_PROVENANCE.txt").write_text(PROV)
    a = EV.arms(share)
    runs = {}
    for k, run in enumerate(a["B"] + a["M"]):
        arm = "B" if run in a["B"] else "M"
        write_run(case, run, arm, k, share)
        runs[run] = {"status": "ok", "reason": None, "dX_E": 0.0, "dX_F": (0.0 if arm == "B" else 9.9e-9), "n_snaps": 54, "last_step": 54000}
    (case / "_band_ab").mkdir()
    (case / EV.SERIES_JSON).write_text(json.dumps({"tag": EV.TAG, "geom_ref": EV.GEOM_REF, "geometry": GEOM, "runs": runs}))
    refresh_series_record(case)
    (case / EV.E4_EVAL_JSON).write_text(json.dumps({"evaluator_sha256": sha(C / "e4_recal_eval.py"),
                                                    "final": {"status": EV.E4_ADOPTED[2], "Md_moc_offset": MD, "Md_moc_offset_repr": MDR,
                                                              "reference_run": EV.E4V_RUN}}))
    (case / EV.LAUNCH_JSON).write_text(json.dumps({"dry": False, "md_offset_repr": MDR, "share_e4v": share, "arms": a}))
    return case


x4 = np.array([-5.0, -1.0, 2.0, 90.0])
s5 = np.array([0.0, 0.1, 0.4, 0.75, 1.0])
GX = np.repeat(x4[:, None], s5.size, axis=1)
GR = np.array([3.0, 1.2, 1.5, 9.4])[:, None] * s5[None, :]


def ind_const(base=0.0, spread=0.2):
    return {r_: {"n": 100, "max": base + spread, "min": base - spread, "q99": abs(base) + spread / 2, **{s: 0.0 for s in EV2.FRAC_STATS}}
            for r_ in EV2.REGION_KEYS}


def stub_wall(case, run, other):
    out = {"own": {"status": "consistent"}}
    if other is not None:
        out["vs_other"] = {"status": "ok", "n_discriminable": 12, "n_discriminable_matching_own": 12}
    return out


def stub_coord(case, run, geom, shifts):
    return {"step": W[-1], "base": {c: 0.0 for c in EV.QS_COLS}, "diff": {n: {c: 0.0 for c in EV.QS_COLS} for n in shifts}}


TMP = Path(tempfile.mkdtemp(prefix="test_moc_v5d_"))
BASES = {False: build_base(TMP, False), True: build_base(TMP, True)}
_n = [0]


def scenario(mod=None, share=False, md=MDR, m_common=None, t0=None, wall=stub_wall):
    """基準の合成 case を複製し、mod(case) で変えて評価する。m_common: (run, step) → M_common、t0: (run, step) → 主指標。"""
    _n[0] += 1
    case = TMP / f"s{_n[0]}"
    shutil.copytree(BASES[share], case)
    if mod:
        mod(case)
    sv = (EV2.load_geometry, EV2.eval_snapshot, EV4.exit_snapshot)
    EV2.load_geometry = lambda rd: (GX, GR, 0.0768075)
    mc = m_common or (lambda run, s: 6.00002 + 2e-6 * math.sin(1.3 * s / 1000.0))
    tf = t0 or (lambda run, s: ind_const(0.0, 0.05))

    def fake_t0(rd, h5, X, R, eta, ep):
        return {"file": h5.name, "problems": [], "indicators": tf(rd.name, int(h5.stem.split("_")[1])), "sha256": None}

    def fake_exit(h5, ni, nj, eta_last, eta_c):
        m = mc(h5.parent.name, int(h5.stem.split("_")[1]))
        return {"file": h5.name, "problems": [], "exitM_common": m, "exitM_own": m + 3e-5, "n_own": 2}
    EV2.eval_snapshot, EV4.exit_snapshot = fake_t0, fake_exit
    try:
        return EV.evaluate(case, md, share, wall_check=wall, coord_sensitivity=stub_coord)
    finally:
        EV2.load_geometry, EV2.eval_snapshot, EV4.exit_snapshot = sv
        shutil.rmtree(case, ignore_errors=True)


def held(out, needle=None, run=None):
    """総合が「保留 (前提不成立)」で、(あれば) run の理由に needle を含む。"""
    if not out["overall"].startswith("保留 (前提不成立)"):
        return False
    if needle is None:
        return True
    pool = out["preconditions_failed"].get(run, []) if run else [x for v in out["preconditions_failed"].values() for x in v]
    return any(needle in x for x in pool)


def edit_json(p: Path, fn):
    d = json.loads(p.read_text())
    fn(d)
    p.write_text(json.dumps(d))


B0, B1, B2 = a0["B"]
M0, M1, M2 = a0["M"]

out = scenario()
check(f"結合 陽性: 許容幅内・据え置き・腕 M の M = 6 達成 (総合 {out['overall'][:60]})",
      out["overall"].startswith("許容幅内") and out["exit_calibration"]["verdict"] == EV.LBL_CAL_KEEP
      and out["armM_M6"]["verdict"] == EV.LBL_M6_OK and not out["preconditions_failed"])
check("結合 陽性: 準定常は exit_M_dev だけ絶対許容内、他は両窓 STEADY",
      all(m["window"]["status"]["exit_M_dev"] in ("STEADY", EV.LBL_ABS) and
          all(v == "STEADY" for c, v in m["window"]["status"].items() if c != "exit_M_dev") for m in out["run_meta"].values()))
check("結合 陽性: 出力に評価器の sha256・登録 99431498・較正値・窓・範囲の注記・M_common の参考の 3 区分",
      out["evaluator_sha256"] == sha(C / "moc_v5d_eval.py") and out["plan_registration_commit"] == "99431498"
      and out["md_offset_repr"] == MDR and out["windows"]["win13"] == list(W) and "等エントロピー" in out["scope_note"]
      and out["exit_calibration_M_common_record_only"] is not None and out["eta_common"]["n"] == 12)
check("結合 陽性: 残差は停滞のみ (PASS と書かない)", all("PASS ではない" in m["convergence"]["label"] for m in out["run_meta"].values()))
_js = json.loads(json.dumps(out, ensure_ascii=False, default=str))
check("結合 陽性: 出力が JSON に書ける (main と同じ json.dumps)", _js["overall"] == out["overall"] and len(_js["rows"]) == len(EV.QS_COLS))

def off_M(d):
    def mod(case):
        for k, r_ in enumerate(a0["M"]):
            write_series(case / r_, 3 + k, offset=d)
        refresh_series_record(case)
    return mod


out = scenario(off_M({"P_slope_abs_eta0.1": 0.05}))
check(f"結合: 腕 M の |P 傾き| +0.05 → 悪化 ({out['verdicts'].get('P_slope_abs_eta0.1')})", out["overall"].startswith("悪化")
      and out["verdicts"]["P_slope_abs_eta0.1"] == EV.LBL_WORSE)
out = scenario(off_M({"P_slope_abs_eta0.1": 0.0299}))
check(f"結合: 腕 M の |P 傾き| +0.0299 → 保留 (ユーザ判断) ({out['overall'][:40]})", out["overall"].startswith("保留 (ユーザ判断")
      and out["verdicts"]["P_slope_abs_eta0.1"] == EV.LBL_HOLD)
out = scenario(off_M({"P_slope_abs_eta0.1": -0.04}))
check("結合: 腕 M の |P 傾き| −0.04 (改善) → 許容幅内 (片側)", out["overall"].startswith("許容幅内")
      and "P_slope_abs_eta0.1" in out["overall"])
out = scenario(off_M({"exit_core_M": 2.0e-4}))
check(f"結合: 腕 M の出口コア M +2e-4 (整定、前提成立) → 出口較正をやり直す ({out['exit_calibration']['verdict'][:12]})",
      not out["preconditions_failed"] and out["exit_calibration"]["verdict"] == EV.LBL_CAL_REDO and out["exit_calibration"]["decision_complete"])


def rising(case):
    write_resid(case / B1, {"rms_roUy": "rising"})


out = scenario(rising)
check("結合: 腕 B の 1 本の残差に RISING → 保留 (前提不成立 [P1])", held(out, "[P1]", B1) and out["exit_calibration"]["verdict"].startswith("保留"))
out = scenario(t0=lambda run, s: ind_const(0.0, 1.5 if run == M2 else 0.05))
check("結合: 腕 M の 1 本の全温 ±1.5 K → 保留 ([P2a])", held(out, "[P2a]", M2))
out = scenario(t0=lambda run, s: ind_const(0.0, 0.05 + (0.15 if (run == B0 and s == W[3]) else 0.0)))
check("結合: 腕 B の 1 本の全温の時間の幅 0.15 K > 0.1 K → 保留 ([P2b])", held(out, "[P2b]", B0))


def drift_B(case):
    write_series(case / B2, 2, fn={"overshoot_eta0.1": lambda st: 0.030 + 0.012 * st / 54000.0})
    refresh_series_record(case)


out = scenario(drift_B)
check("結合: 腕 B の 1 本だけ非 STEADY (対称に要求) → 保留", held(out, "準定常でない量", B2))


def width_only(case):
    write_series(case / M0, 3, fn={"M_wave_eta0.1": lambda st: 0.0064 + (1.5e-4 if st == W[6] else 0.0)})
    refresh_series_record(case)


out = scenario(width_only)
st_w = out["run_meta"][M0]["window"]["status"]["M_wave_eta0.1"]
check(f"結合: M 波の 13 枚の幅 1.5e-4 > Δq/10 (準定常は {st_w}) → 保留 (幅)", held(out, "M_wave_eta0.1 の 13 枚の幅", M0) and st_w == "STEADY")


def exit_width(case):
    write_series(case / M1, 4, fn={"exit_core_M": lambda st: 5.99999 + (1.2e-5 if st == W[2] else 0.0)})
    refresh_series_record(case)


out = scenario(exit_width)
check("結合: 符号付き出口 M の幅 1.2e-5 > 1e-5 (exit_M_dev は例外の 1.8e-5 以内) → 保留 (幅)", held(out, "符号付きの出口 M", M1))


def forge_sha(case):
    (case / M1 / "RUN_PROVENANCE.txt").write_text(PROV.replace(FSHA, "0" * 64))


out = scenario(forge_sha)
check("結合: 腕 M の 1 本の forge の sha256 が run_0143 と違う → 保留", held(out, "forge の sha256", M1))


def md_run(case):
    edit_json(case / B2 / "prepare_info.json", lambda d: d.update(Md_moc_offset=MD * (1 + 1e-12)))


out = scenario(md_run)
check("結合: 腕 B の 1 本の Md_moc_offset が --md-offset と違う → 保留 (両腕で同じ値の検査)", held(out, "Md_moc_offset", B2))
out = scenario(md="6.88e-06")
check("結合: 評価の --md-offset が E4・起動の記録・run と違う → 保留", held(out, "[E4]") and held(out, "[起動の記録]"))


def e4_hold(case):
    edit_json(case / EV.E4_EVAL_JSON, lambda d: d["final"].update(status="保留"))


out = scenario(e4_hold)
check("結合: E4 の最終が採用でない (保留) → 保留", held(out, "採用でない"))


def e4_sha(case):
    edit_json(case / EV.E4_EVAL_JSON, lambda d: d.update(evaluator_sha256="0" * 64))


out = scenario(e4_sha)
check("結合: E4 の評価器の sha256 が今の e4_recal_eval.py と違う → 保留", held(out, "E4 の評価器の sha256"))


def launch_share(case):
    edit_json(case / EV.LAUNCH_JSON, lambda d: d.update(share_e4v=True))


out = scenario(launch_share)
check("結合: 起動の記録の共用が評価の引数と違う → 保留", held(out, "起動の記録の共用"))


def xf_shift(case):
    edit_json(case / EV.SERIES_JSON, lambda d: d["runs"][M2].update(dX_F=2e-6))


out = scenario(xf_shift)
check("結合: 腕 M の 1 本の X_F のずれ 2e-6 r_t > 1e-6 → 保留", held(out, "座標の整合の許容差", M2))


def stale_csv(case):
    with open(case / B0 / EV.SERIES_CSV, "a") as f:
        f.write("")
    (case / B0 / EV.SERIES_QS).write_text("別の呼び出し\n")


out = scenario(stale_csv)
check("結合: 準定常の記録ファイルが時系列の記録と違う (古い成果物) → 保留", held(out, "時系列の記録と違う", B0))


def g1_mesh(case):
    edit_json(case / B1 / "prepare_info.json", lambda d: d["mesh"].update(source="mesh"))


out = scenario(g1_mesh)
check("結合: 旧 G1 の格子 (採用元 mesh) の run を混ぜる → 保留", held(out, "mesh_euler でない", B1))


def gate_fail(case):
    edit_json(case / M0 / "prepare_info.json", lambda d: d["moc"]["gate"].update(**{"pass": False}))


out = scenario(gate_fail)
check("結合: 腕 M の MOC のゲート不合格 → 保留", held(out, "ゲートが合格でない", M0))
out = scenario(wall=lambda case, run, other: {"own": {"status": "consistent"},
                                              **({"vs_other": {"status": "ok", "n_discriminable": 12, "n_discriminable_matching_own": 7}}
                                                 if other else {})})
check("結合: 腕 M の壁節点が腕 B の当てはめに一致 (取り違え) → 保留", held(out, "壁の取り違え"))


def spline_mismatch(case):
    edit_json(case / B2 / "prepare_info.json", lambda d: d["wall_fit"].update(spline=spline_dict(3e-6)))


out = scenario(spline_mismatch)
check("結合: 腕 B の中で保存 spline がそろわない → 保留", held(out, "腕 B の run の保存 spline"))


def same_wall(case):
    for r_ in a0["M"]:
        edit_json(case / r_ / "prepare_info.json", lambda d: d["wall_fit"].update(spline=spline_dict(0.0)))


out = scenario(same_wall)
check("結合: 腕 B と腕 M の保存 spline が同じ (MOC の変更が壁に入っていない) → 保留", held(out, "MOC の変更が壁に入っていない"))


def cfg_mismatch(case):
    edit_json(case / M2 / EV.PREP_JSON, lambda d: d["files_sha256"].update(**{"solverConfig.yaml": "other"}))


out = scenario(cfg_mismatch)
check("結合: prep 時の solverConfig が腕 B と腕 M で違う → 保留", held(out, "prep 時の solverConfig.yaml が run の間で違う"))


def rc1(case):
    (case / M1 / "RUN_RC").write_text("1\n")
    (case / B0 / "EARLY_STOP.txt").write_text("x\n")


out = scenario(rc1)
check("結合: RUN_RC 1・EARLY_STOP.txt → 保留", held(out, "RUN_RC", M1) and held(out, "EARLY_STOP", B0))


def dry_prep(case):
    edit_json(case / B1 / "prepare_info.json", lambda d: d.update(DRY=True))


out = scenario(dry_prep)
check("結合: 乾式確認の prep から作った run → 保留", held(out, "乾式確認", B1))


def drop_run(case):
    shutil.rmtree(case / M2)


out = scenario(drop_run)
check("結合: 腕 M の 1 本が無い → 保留 (欠損を既定値で埋めない)", held(out, "欠損", M2) and out["rows"] == [] and M2 in out["missing"])
out = scenario(m_common=lambda run, s: 6.00012 if run == M1 else 6.00002 + 2e-6 * math.sin(1.3 * s / 1000.0))
check(f"結合: 腕 M の 1 本の M_common が +1.2e-4 (整定) → M = 6 未達 (NS を始めない) ({out['armM_M6']['verdict'][:10]})",
      out["armM_M6"]["verdict"] == EV.LBL_M6_NG and out["overall"].startswith("許容幅内"))
out = scenario(m_common=lambda run, s: 6.00002 + (6e-5 if (run == M0 and s == W[5]) else 0.0))
check("結合: 腕 M の 1 本の M_common の幅 6e-5 (E4 の P3 不成立) → M = 6 は判別不能 (比較の判定とは別)",
      out["armM_M6"]["verdict"] == EV.LBL_M6_UNDET and out["overall"].startswith("許容幅内"))
out = scenario(m_common=lambda run, s: 6.00012 if run == B0 else 6.00002)
check("結合: 腕 B の M_common が外れても腕 M の M = 6 の判定には入らない (腕 B は記録だけ)",
      out["armM_M6"]["verdict"] == EV.LBL_M6_OK and out["armB_M6_record_only"][B0]["within"] is False)

# --- 共用 (run_0164 を腕 B に) --------------------------------------------------------------------------------------------------
out = scenario(share=True)
check(f"共用: 陽性 (腕 B = run_0171・run_0172・run_0164) → 許容幅内 ({out['overall'][:30]})",
      out["overall"].startswith("許容幅内") and out["arms"]["B"][2] == EV.E4V_RUN and not out["preconditions_failed"])


def e4_ref(case):
    edit_json(case / EV.E4_EVAL_JSON, lambda d: d["final"].update(reference_run="run_0163_euler_e4_recal_d0"))


out = scenario(e4_ref, share=True)
check("共用: E4 の参照 run が run_0164 でない → 保留", held(out, "参照 run"))


def e4v_wall(case):
    edit_json(case / EV.E4V_RUN / "prepare_info.json", lambda d: d["wall_fit"].update(spline=spline_dict(2e-6)))


out = scenario(e4v_wall, share=True)
check("共用: run_0164 の壁が腕 B の新しい run と違う → 保留", held(out, "腕 B の run の保存 spline"))


def e4v_prep_stage(case):
    edit_json(case / EV.E4V_RUN / EV4.PREP_JSON, lambda d: d.update(stage="d0"))


out = scenario(e4v_prep_stage, share=True)
check("共用: run_0164 の E4_PREP.json が段 d1 でない → 保留", held(out, "d1", EV.E4V_RUN))

# --- moc_v5d: 問題・run の前の停止 ----------------------------------------------------------------------------------------------
import moc_v5d as RUN  # noqa: E402

with tempfile.TemporaryDirectory() as td:
    case = Path(td)
    for n in (RUN.BASE_B, RUN.V5_M_PROBLEM):
        shutil.copy(C / n, case / n)
    texts = RUN.problem_texts(MDR, case)
    for arm, t in texts.items():
        (case / EV.PROBLEMS[arm]).write_text(t)
    bad, _ = RUN.check_problems(MD, False, case=case)
    check(f"問題: 生成した B・M は検査 OK ({bad[:1]})", bad == [])
    dB, dM = (yaml.safe_load(texts[x]) for x in ("B", "M"))
    check("問題: B は MOC のキー無し・M は analytic・converge・Md_moc_offset は float で完全一致・name が腕ごと",
          "moc_axis_limit" not in dB["geometry"] and dM["geometry"]["moc_axis_limit"] == "analytic" and dM["geometry"]["moc_corrector"] == "converge"
          and dB["geometry"]["Md_moc_offset"] == MD == dM["geometry"]["Md_moc_offset"] and isinstance(dB["geometry"]["Md_moc_offset"], float)
          and dB["name"] != dM["name"])
    for name, arm, edit in (("腕 M の corrector を fixed2", "M", lambda t: t.replace("moc_corrector: converge", "moc_corrector: fixed2")),
                            ("腕 B の mesh_euler", "B", lambda t: t.replace("  wall_first_frac: 0.005\n", "  wall_first_frac: 0.004\n", 1)),
                            ("腕 M の Md_moc_offset", "M", lambda t: t.replace(f"Md_moc_offset: {MDR}", "Md_moc_offset: 6.9e-06")),
                            ("腕 B の mono_r2", "B", lambda t: t.replace("wall_fit_mono_r2: [0.0, 1.5]", "wall_fit_mono_r2: [0.0, 1.4]"))):
        t2 = edit(texts[arm])
        assert t2 != texts[arm], name
        (case / EV.PROBLEMS[arm]).write_text(t2)
        check(f"問題: {name} を変えると止まる", bool(RUN.check_problems(MD, False, case=case)[0]))
        (case / EV.PROBLEMS[arm]).write_text(texts[arm])
    check("問題: 共用で E4V の問題が無い → 本番は止まる・乾式は未確認として通す",
          bool(RUN.check_problems(MD, True, case=case)[0]) and RUN.check_problems(MD, True, dry=True, case=case)[0] == [])
    d0t = (case / RUN.BASE_B).read_text()
    (case / RUN.E4V_PROBLEM).write_text(d0t.replace("name: isobutane_m6_d155_euler_e4_recal_d0", "name: isobutane_m6_d155_euler_e4_recal_d1")
                                        .replace("Md_moc_offset: 3.7700e-04", f"Md_moc_offset: {MDR}"))
    check("問題: 共用で E4V の問題 (name・Md_moc_offset を変えた d0) と腕 B の差が name だけ → OK", RUN.check_problems(MD, True, case=case)[0] == [])
    (case / RUN.E4V_PROBLEM).write_text(d0t.replace("name: isobutane_m6_d155_euler_e4_recal_d0", "name: x"))
    check("問題: 共用で E4V の問題の Md_moc_offset が違う → 止まる", bool(RUN.check_problems(MD, True, case=case)[0]))
    (case / "_band_ab").mkdir()
    (case / EV.E4_EVAL_JSON).write_text(json.dumps({"evaluator_sha256": sha(C / "e4_recal_eval.py"), "final": {"status": "保留"}}))
    try:
        RUN.make_problems(MDR, False, dry=False, case=case)
        ok_ = False
    except SystemExit as e:
        ok_ = "E4 の結果と合わない" in str(e)
    check("make-problems: 本番で E4 の最終が採用でない → 作らない", ok_)
    rec = RUN.make_problems(MDR, False, dry=True, case=case)
    check("make-problems: 乾式は E4 の不一致を記録だけして作る", rec["e4_adoption"].get("status") == "記録だけ (DRY)" and rec["check"]["problems"] == [])
    (case / EV.E4_EVAL_JSON).write_text(json.dumps({"evaluator_sha256": sha(C / "e4_recal_eval.py"),
                                                    "final": {"status": EV.E4_ADOPTED[2], "Md_moc_offset": MD, "reference_run": EV.E4V_RUN}}))
    rec = RUN.make_problems(MDR, False, dry=False, case=case)
    check("make-problems: 本番で E4 の採用値と同じ → 作る (中身が同じなら書き直さない)", rec["e4_adoption"]["ok"] and rec["written"] == [])
    (case / EV.E4_EVAL_JSON).write_text(json.dumps({"evaluator_sha256": sha(C / "e4_recal_eval.py"),
                                                    "final": {"status": EV.E4_ADOPTED[2], "Md_moc_offset": 7e-6, "reference_run": EV.E4V_RUN}}))
    rec = RUN.make_problems("7e-06", False, dry=False, case=case)
    check("make-problems: run が無い間は別の値で書き直す", sorted(rec["written"]) == sorted(EV.PROBLEMS.values())
          and yaml.safe_load((case / EV.PROBLEMS["B"]).read_text())["geometry"]["Md_moc_offset"] == 7e-6)
    (case / M0).mkdir()
    (case / EV.E4_EVAL_JSON).write_text(json.dumps({"evaluator_sha256": sha(C / "e4_recal_eval.py"),
                                                    "final": {"status": EV.E4_ADOPTED[2], "Md_moc_offset": MD, "reference_run": EV.E4V_RUN}}))
    try:
        RUN.make_problems(MDR, False, dry=False, case=case)
        ok_ = False
    except SystemExit as e:
        ok_ = "書き直さない" in str(e)
    check("make-problems: V5d の run dir があるときは別の値で書き直さない", ok_)
    tok = __import__("e4_recal").yaml_float("1e-05")
    check("問題: YAML の float の書き方 (1e-05 → 1.0e-05)", tok == "1.0e-05")

check("起動の順: 腕を交互に (B1 M1 B2 M2 B3 M3 / 共用 B1 M1 B2 M2 M3)",
      RUN.launch_order(False) == [B0, M0, B1, M1, B2, M2] and RUN.launch_order(True) == [B0, M0, B1, M1, M2])
check("plan-runs: 共用の印・全 run・基準 run", "SHARED run_0164_euler_e4_recal_d1" in RUN.plan_runs(True)
      and f"GEOM_REF {EV.GEOM_REF}" in RUN.plan_runs(False) and "SHARED -" in RUN.plan_runs(False))

with tempfile.TemporaryDirectory() as td:
    rd = Path(td) / M0
    rd.mkdir()
    info = {"v5d": {"arm": "M", "md_offset_repr": MDR, "share_e4v": False}, "ic": {"mode": "isentropic", "VERDICT": "OK"},
            "moc": {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}},
            "wall_fit": {"mono_r2": [0.0, 1.5]}, "Md_moc_offset": MD, "scale_m": EV4.SCALE_M,
            "initial_line": {"run": EV4.INITIAL_LINE[0], "res": [EV4.INITIAL_LINE[1]]},
            "mesh": {"ni": 2000, "source": "mesh_euler", "params": dict(EV4.MESH_EULER_EXPECT)}}

    def run_refused(inf, cfg_text, name=None):
        d = rd if name is None else rd.parent / name
        d.mkdir(exist_ok=True)
        (d / "prepare_info.json").write_text(json.dumps(inf))
        (d / "solverConfig.yaml").write_text(cfg_text)
        try:
            RUN.run(d)
            return False, "起動した"
        except SystemExit as e:
            return True, str(e)
    ok_, msg = run_refused(info, cfg("2.0", 1, 18000, 1000))
    check(f"run: step 数 18000 は forge の前に止まる ({msg[-30:]})", ok_ and "同期していない" in msg)
    ok_, msg = run_refused(dict(info, DRY=True), cfg("2.0", 1, 54000, 1000))
    check("run: 乾式確認の prep は止まる", ok_ and "乾式" in msg)
    ok_, msg = run_refused(dict(info, mesh=dict(info["mesh"], source="mesh")), cfg("2.0", 1, 54000, 1000))
    check("run: 格子の採用元が mesh (旧 G1) は止まる", ok_ and "固定の条件" in msg)
    ok_, msg = run_refused(dict(info, v5d={"arm": "B", "md_offset_repr": MDR, "share_e4v": True}), cfg("2.0", 1, 54000, 1000), name=B2)
    check("run: 共用の構成に run_0173 (作らない run) は止まる", ok_ and "run 名でない" in msg)
    ok_, msg = run_refused(dict(info, Md_moc_offset=MD * 2), cfg("2.0", 1, 54000, 1000))
    check("run: prepare_info の Md_moc_offset が v5d の記録と違う → 止まる", ok_ and "固定の条件" in msg)
    (rd / "res_0.h5").write_text("x")
    ok_, msg = run_refused(info, cfg("2.0", 1, 54000, 1000))
    check("run: 既に出力がある → 止まる", ok_ and "既に出力" in msg)
    prep = Path(td) / "_prep_v5d_M"
    prep.mkdir()
    (prep / "a.txt").write_text("1")
    (prep / EV.PREP_JSON).write_text(json.dumps({"runs": EV.new_runs(False)["M"], "dry": False, "ic_check": {"VERDICT": "OK"},
                                                 "files_sha256": {"a.txt": sha(prep / "a.txt")}}))
    r2 = Path(td) / M1
    r2.mkdir()
    (r2 / "a.txt").write_text("1")
    (r2 / "prepare_info.json").write_text("{}")
    check("verify-prep: 同じ入力・腕の run → OK", RUN.verify_prep(prep, [r2]) == [])
    (r2 / "a.txt").write_text("2")
    check("verify-prep: 入力が変わった → 不成立", bool(RUN.verify_prep(prep, [r2])))
    r3 = Path(td) / B0
    r3.mkdir()
    (r3 / "a.txt").write_text("1")
    (r3 / "prepare_info.json").write_text("{}")
    check("verify-prep: 別の腕の run 名 → 不成立", bool(RUN.verify_prep(prep, [r3])))

# --- moc_v5d: 共用の照合 (腕 B の prep と run_0164)・共用の前提 ----------------------------------------------------------------
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    d = td / "_prep_v5d_B"
    d.mkdir()
    for n in EV.SAME_FILES:
        (d / n).write_text(n)
    info = {"wall_fit": {"spline": spline_dict()}, "mesh": {"params": {"ni": 2000}, "hashes": {"coord_sha256": "c"}}, "Md_moc_offset": MD,
            "x_E": 39.8, "moc": {"axis_limit": "legacy", "corrector": "fixed2"}, "scale_m": EV4.SCALE_M,
            "initial_line": {"run": "/aws/x/run_0062_euler_wallfit_fit_r1_ext6k", "res": ["res_6000.h5"], "sha256_16": {"line": "a"}}}
    e = td / "runs" / EV.E4V_RUN
    e.mkdir(parents=True)
    (e / EV4.PREP_JSON).write_text(json.dumps({"files_sha256": {n: sha(d / n) for n in EV.SAME_FILES}}))
    sv = RUN.TM.RUNS
    RUN.TM.RUNS = td / "runs"
    try:
        (e / "prepare_info.json").write_text(json.dumps(dict(info, initial_line=dict(info["initial_line"], run="/other/run_0062_euler_wallfit_fit_r1_ext6k"))))
        r_ = RUN.compare_with_e4v(d, info, dry=False)
        check(f"共用の照合: 初期線の絶対パスだけが違う → OK ({r_['failures']})", r_["status"] == "OK")
        (e / "prepare_info.json").write_text(json.dumps(dict(info, Md_moc_offset=MD * 2)))
        check("共用の照合: run_0164 の Md_moc_offset が違う → 不成立", RUN.compare_with_e4v(d, info, dry=False)["status"] == "FAIL")
        (e / "prepare_info.json").write_text(json.dumps(dict(info, mesh={"params": {"ni": 2000}, "hashes": {"coord_sha256": "other"}})))
        check("共用の照合: 格子のハッシュが違う → 不成立", RUN.compare_with_e4v(d, info, dry=False)["status"] == "FAIL")
        (e / "prepare_info.json").write_text(json.dumps(info))
        (d / "solverConfig.yaml").write_text("changed")
        check("共用の照合: prep 時の solverConfig が違う → 不成立", RUN.compare_with_e4v(d, info, dry=False)["status"] == "FAIL")
        check("共用の前提: run_0164 に RUN_RC・判定窓の res・残差判定が無い → 不成立 (本番)", len(RUN.check_share(MD, dry=False)) >= 3)
        (e / "RUN_RC").write_text("0\n")
        (e / EV4.SEGMENT_VERDICT_FILE).write_text("x")
        for s in W:
            (e / f"res_{s}.h5").write_text("x")
        check("共用の前提: 完走・窓の res・残差判定・Md_moc_offset がそろう → OK", RUN.check_share(MD, dry=False) == [])
        check("共用の前提: --md-offset が run_0164 と違う → 不成立", bool(RUN.check_share(MD * 2, dry=False)))
        shutil.rmtree(e)
        check("共用の照合: run_0164 が無い → 本番は不成立・乾式は未確認",
              bool(RUN.compare_with_e4v(d, info, dry=False)["failures"]) and "未確認" in RUN.compare_with_e4v(d, info, dry=True)["status"]
              and bool(RUN.check_share(MD, dry=False)) and RUN.check_share(MD, dry=True) == [])
    finally:
        RUN.TM.RUNS = sv

shutil.rmtree(TMP, ignore_errors=True)
print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
