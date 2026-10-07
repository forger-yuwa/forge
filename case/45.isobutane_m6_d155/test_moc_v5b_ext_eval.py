"""moc_v5b_ext_eval.py と run_moc_v5b_ext.sh の試験 (plan discretization-moc-axis-limit-and-corrector §6 V5b)。forge・AWS は使わない。
合成した系列と、一時ディレクトリの模擬 case (V5 の親 7 本 + 延長の子 7 本の CSV・JSON・設定) だけを使う。V5・V5b の run の結果は読まない。
  接続: 重複・欠落・余分・逆順・子の local 0・非有限・列の欠損を拒否、正しい接続は通算 1000〜54000。
  分類: 減衰する系列 → 減衰側 / 3 周期以上の定振幅振動 → 持続振動側 / 周期不足・振幅の減衰が中途半端・折り返し (1 間隔ごと)・
    窓の平均の移動だけ・両方に当たる → 判別不能 / 符号付き出口 M の許容 1e-5。
  総合: 7 本すべてが減衰側 → 窓 B で V5 と同じ比較 / 1 本でも減衰側でない → 窓 B の差は参考値 / 早期停止 (NaN) → 保留 /
    設定・バイナリ・環境変数・段・評価座標・評価器・RUN_RC の欠損 → 保留 (既定で通さない) / V5 の静的な前提 (moc_v5_euler_eval.preconditions)
    の親の判定区間の理由だけを外す。
  延長の準備: prep-child・restart_field (本物)・verify-child を合成の h5 で (保存量の取りこぼし・設定の違いは不成立)。
  起動スクリプト: DRY=1 で合成の親 7 本から準備・restart・検査まで通り forge を起動しない / 本番で子が既にある・バイナリ違いは止まる。
usage: python3 test_moc_v5b_ext_eval.py
"""
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import make_interp_spline

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import moc_v5b_ext_eval as E  # noqa: E402
import moc_v5_euler_eval as EV  # noqa: E402

sys.path.insert(0, str(C.parents[1] / "solver_density_cuda/tools"))
from check_quasisteady import classify_series  # noqa: E402

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def edit_json(p: Path, fn):
    d = json.loads(p.read_text()); fn(d); p.write_text(json.dumps(d))


V5_SHA_BEFORE = sha(C / "_band_ab/moc_v5_euler_eval.json") if (C / "_band_ab/moc_v5_euler_eval.json").is_file() else None
COLS = ["step", "M_wave_eta0.0", "P_wave_eta0.0", "overshoot_eta0.0", "P_slope_eta0.0", "M_wave_eta0.1", "P_wave_eta0.1",
        "overshoot_eta0.1", "P_slope_eta0.1", "exit_core_M", "exit_M_dev", "P_slope_abs_eta0.0", "overshoot_exitnorm_eta0.0",
        "P_slope_abs_eta0.1", "overshoot_exitnorm_eta0.1"]
BASE = {"M_wave_eta0.1": 0.0064, "P_wave_eta0.1": 0.041, "overshoot_eta0.1": 0.038, "overshoot_exitnorm_eta0.1": 0.038,
        "P_slope_abs_eta0.1": 0.19, "exit_core_M": 5.99999,
        "M_wave_eta0.0": 0.02, "P_wave_eta0.0": 0.05, "overshoot_eta0.0": 0.04, "P_slope_eta0.0": 0.3, "P_slope_eta0.1": -0.19,
        "P_slope_abs_eta0.0": 0.3, "overshoot_exitnorm_eta0.0": 0.04}
TOLS = E.TOL


def decay_fn(q):
    """減衰する過渡: q(s) = 基準 + 10·tol·exp(−s/4000) (通算 s)。"""
    return lambda s: BASE[q] + 10.0 * TOLS[q] * math.exp(-s / 4000.0)


def series_values(spec: dict, s: int, rng) -> dict:
    """通算 step s の 1 行。spec {量: s → 値} が無い量は減衰する過渡 + 小さな雑音 (振れ 0.02·tol、転回点のしきい値より十分小さい)。"""
    row = {}
    for q in BASE:
        if q in spec:
            row[q] = spec[q](s)
        elif q in TOLS:
            row[q] = decay_fn(q)(s) + 0.02 * TOLS[q] * (rng.random() - 0.5)
        else:
            row[q] = BASE[q]
    row["exit_M_dev"] = abs(row["exit_core_M"] - 6.0)
    return row


def write_csv(path: Path, steps, spec: dict, offset: int, seed: int, drop=(), extra_rows=(), nan_at=None):
    rng = np.random.default_rng(seed)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for st in steps:
            row = series_values(spec, st + offset, rng)
            if st in drop:
                continue
            vals = [st] + [row[c] for c in COLS[1:]]
            if nan_at is not None and st == nan_at:
                vals[COLS.index("overshoot_eta0.1")] = float("nan")
            w.writerow(vals)
        for r in extra_rows:
            w.writerow(r)


def write_qs(csv_path: Path, qs_path: Path, tail: int = 5):
    """eval_wallfit_euler (check_quasisteady --series-csv) と同じ形の準定常ファイル (末尾 tail 枚)。"""
    rows = list(csv.DictReader(open(csv_path)))
    steps = [float(r["step"]) for r in rows]
    lines = [f"\n=== {csv_path}  [{len(rows)} rows]  -> X ==="]
    for c in COLS[1:]:
        v, d, _ = classify_series(steps, [float(r[c]) for r in rows], (tail - 0.5) / len(rows), 0.05, 0.10, 4)
        lines.append(f"  {c:16s}: {d:55s} {v}")
    qs_path.write_text("\n".join(lines) + "\n")


# --- 模擬 case (V5 の親 7 本 + 子 7 本) -------------------------------------------------------------------------------
X_E0, X_F0 = 39.80895834812345, 95.24486827926533
GEOM = {"X_E": X_E0, "X_F": X_F0, "WIN_T": [X_E0 + 2.0, X_F0 - 1.0], "WIN_O": [X_E0 - 15.0, X_F0]}
SOLVER_P = ('mesh: {discretization: "node", isAxisymmetric: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}\n'
            'gpu: 1\nsolver: "SLAU"\nphysProp: {thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}\n'
            "time:\n  unsteady: 0\n  last: {nStepOuter: 18000}\n"
            "  deltaT: {control: 1, cfl: 2.0, cfl_pseudo: 2.0, blockDPLUR: 1, implicitRelax: 0.7, detectNaN: 1}\n"
            "  outStepStart: 0\n  outStepInterval: 1000\n  nStepInner: 5\nspace: {convMethod: 1, limiter: 2}\n")
SOLVER_C = SOLVER_P.replace("nStepOuter: 18000", "nStepOuter: 36000")
BCOND = "inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 5500000.0, Tt: 1600.0}}\n"
PROV = ("forge_sha256: " + "a" * 64 + "\nenv         : FORGE_CUDA_BLOCKSIZE=128\nenv         : FORGE_BIN=/x/forge\n"
        "'slauWallNormalChi' effective: 1 (auto)\n'scalarGradient' effective: lsq (default)\n")
SEG_PLATEAU = "[segment] 判定区間 main\n=== x  [last step 35999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===\n"
SEG_CONVERGING = "[segment] 判定区間 main\n=== x  [last step 17999]  -> NOT CONVERGED (still converging) ===\n"
SEG_DIVERGED = "[segment] 判定区間 main\n=== x  [last step 20999]  -> DIVERGED (NaN/Inf) ===\n"
MAIN_KEY = {"convMethod": "1", "limiter": "2", "bcond_sha1": "abc"}
PARENTS = list(E.CHILDREN.values())
CHILDREN = list(E.CHILDREN)
SER_SHA = "e" * 64


def spline_dict(shift=0.0):
    x = np.linspace(-1.0, 3.0, 801)
    s_ = make_interp_spline(x, 1 + x ** 2 / 4 - x ** 3 / 30 + shift, k=5)
    return {"t": [float(v) for v in s_.t], "c": [float(v) for v in s_.c], "k": 5}


def mock_parent(root: Path, run: str, k: int):
    """V5 の静的な前提を満たす親 (moc_v5_euler_eval の試験の模擬と同じ形) と V5 の時系列。"""
    d = root / run
    d.mkdir(parents=True)
    write_csv(d / "wallfit_series_v5.csv", range(1000, 18001, 1000), {}, 0, seed=100 + k)
    write_qs(d / "wallfit_series_v5.csv", d / "QUASISTEADY_wallfit_v5.txt")
    (d / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_PLATEAU)
    (d / "solverConfig.yaml").write_text(SOLVER_P)
    (d / "bcondConfig.yaml").write_text(BCOND)
    (d / "RUN_PROVENANCE.txt").write_text(PROV)
    (d / "stage_manifest.json").write_text(json.dumps({"manifest_version": 2, "stages": [
        {"tag": "soft", "key": {**MAIN_KEY, "convMethod": "0"}}, {"tag": "main", "key": MAIN_KEY}]}))
    role = "B" if run in EV.ARMS["B"] else ("M" if run in EV.ARMS["M"] else "ISEN")
    info = {"scale_m": 0.0768075, "x_E": X_E0, "wall_fit": {"mono_r2": [0.0, 1.5], "spline": spline_dict(0.0 if role == "B" else 1e-5)}}
    src = f"/aws/case/{EV.IC_RUN_NAME}/{EV.IC_RES}"
    if role == "B":
        ic = {"VERDICT": "OK", "mode": "index", "src_res": src, "dst_sha256_after": "shaB", "dst_mesh_digest": "gridB"}
    elif role == "M":
        info["moc"] = {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
        info["ic"] = {"VERDICT": "OK", "inspection": {"limit_m": 6.2e-6}}
        ic = {"VERDICT": "OK", "mode": "index", "src_res": src, "dst_sha256_after": "shaM", "dst_mesh_digest": "gridM",
              "checks": {"displacement": {"detail": {"limit_m": 6.2e-6, "max_m": 6.19e-6}}}}
        (d / "IC_INSPECTION.json").write_text(json.dumps({"VERDICT": "OK", "limit_m": 6.2e-6,
                                                          "conditions": {kk: {"ok": True} for kk in EV.IC_INSPECTION_KEYS}}))
    else:
        info["moc"] = {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
        ic = {"VERDICT": "OK", "mode": "isentropic", "dst_sha256_after": "shaI", "dst_mesh_digest": "gridM"}
    info.setdefault("ic", {}).update({"VERDICT": ic["VERDICT"]})
    (d / "prepare_info.json").write_text(json.dumps(info))
    (d / "IC_MAP.json").write_text(json.dumps(ic))
    return {"status": "ok", "reason": None, "dX_E": 0.0, "dX_F": (0.0 if role == "B" else 9.9e-9), "n_snaps": 18, "last_step": 18000,
            "csv_sha256": sha(d / "wallfit_series_v5.csv"), "quasisteady_sha256": sha(d / "QUASISTEADY_wallfit_v5.txt")}


def write_child_series(root: Path, child: str, spec=None, seed=0, last=36000, drop=(), extra_rows=(), nan_at=None):
    d = root / child
    write_csv(d / E.SERIES_CSV, range(1000, last + 1, 1000), spec or {}, E.PARENT_END, seed=seed, drop=drop, extra_rows=extra_rows,
              nan_at=nan_at)
    write_qs(d / E.SERIES_CSV, d / E.QS_FILE)


def refresh_record(root: Path, child: str):
    p = root / E.SERIES_REC
    r = json.loads(p.read_text())
    r["runs"][child]["csv_sha256"] = sha(root / child / E.SERIES_CSV)
    r["runs"][child]["quasisteady_sha256"] = sha(root / child / E.QS_FILE)
    p.write_text(json.dumps(r))


def mock_child(root: Path, child: str, parent: str, k: int, rec5: dict):
    d = root / child
    d.mkdir(parents=True)
    write_child_series(root, child, seed=200 + k)
    (d / "nozzle.h5").write_bytes(b"mock nozzle " + child.encode())
    (d / "solverConfig.yaml").write_text(SOLVER_C)
    (d / "bcondConfig.yaml").write_text(BCOND)
    (d / "RUN_PROVENANCE.txt").write_text(PROV)
    (d / "stage_manifest.json").write_text(json.dumps({"manifest_version": 2, "stages": [{"tag": "main", "key": MAIN_KEY}]}))
    info = json.loads((root / parent / "prepare_info.json").read_text())
    info.update(stages="none", extends=parent)
    (d / "prepare_info.json").write_text(json.dumps(info))
    (d / E.EXT_RECORD).write_text(json.dumps({"VERDICT": "OK", "dry": False, "parent": parent, "child": child,
                                              "src_res": f"{parent}/res_18000.h5", "nozzle_sha256_after": sha(d / "nozzle.h5"),
                                              "failures": []}))
    (d / E.RUN_RC_FILE).write_text("0\n")
    (d / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_PLATEAU)
    return {"status": "ok", "reason": None, "dX_E": rec5["dX_E"], "dX_F": rec5["dX_F"], "n_snaps": 36, "last_step": 36000,
            "csv_sha256": sha(d / E.SERIES_CSV), "quasisteady_sha256": sha(d / E.QS_FILE)}


def mock_case(root: Path):
    runs5, runs5b = {}, {}
    for k, p in enumerate(PARENTS):
        runs5[p] = mock_parent(root, p, k)
    for k, (c, p) in enumerate(E.CHILDREN.items()):
        runs5b[c] = mock_child(root, c, p, k, runs5[p])
    (root / "_band_ab").mkdir()
    base = {"geom_ref": EV.GEOM_REF, "geometry": GEOM, "fixed_coef": True, "evaluator": "eval_wallfit_euler.py", "evaluator_sha256": SER_SHA}
    (root / E.SERIES_REC_V5).write_text(json.dumps({**base, "tag": "v5", "runs": runs5}))
    (root / E.SERIES_REC).write_text(json.dumps({**base, "tag": "v5b", "runs": runs5b}))


def stub_wall(case, run, other):
    out = {"own": {"status": "consistent"}}
    if other is not None:
        out["vs_other"] = {"status": "ok", "n_discriminable": 12, "n_discriminable_matching_own": 12}
    return out


def run_eval(mutate=None, v5_static=None):
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        mock_case(root)
        if mutate:
            mutate(root)
        out = E.evaluate(root, wall_check=stub_wall, v5_static=v5_static)
        out["_v5_json_written"] = (root / "_band_ab/moc_v5_euler_eval.json").exists()
        return out


def reasons(out):
    return " ".join(" ".join(v) for v in out.get("preconditions_failed", {}).values())


# === 定数・登録との対応 ===============================================================================================
check("窓 A = 通算 24000〜36000 の 13 枚、窓 B = 42000〜54000 の 13 枚", E.WIN_A_STEPS == list(range(24000, 36001, 1000))
      and E.WIN_B_STEPS == list(range(42000, 54001, 1000)) and len(E.WIN_A_STEPS) == len(E.WIN_B_STEPS) == 13)
check("周期の区間 = 通算 24000〜54000 の 31 枚 (再開後 6000 step を除く)", E.OSC_STEPS == list(range(24000, 54001, 1000))
      and E.PARENT_END + E.SKIP_AFTER_RESTART == 24000)
check("追加 36000 step・通算 54000・親の終了 step 18000", (E.EXT_STEPS, E.TOTAL_END, E.PARENT_END) == (36000, 54000, 18000))
check("許容: 6 量は Δq/10、符号付き出口 M は 1e-5", all(abs(E.TOL[q] - EV.DQ[q] / 10) < 1e-18 for q in EV.DQ)
      and E.TOL["exit_core_M"] == 1e-5 and set(E.TOL) == set(EV.QS_COLS))
check("3 周期・半周期 ≥ 2 出力間隔", E.MIN_PERIODS == 3 and E.MIN_HALF_INTERVALS == 2)
check("run 名: run_0154〜0160 が親 run_0143〜0145・0150〜0153 に対応", CHILDREN == [
    "run_0154_euler_wallfit_monoG1_r1_ext36k", "run_0155_euler_wallfit_monoG1_r2_ext36k", "run_0156_euler_wallfit_monoG1_r3_ext36k",
    "run_0157_euler_wallfit_mocG1_r1_ext36k", "run_0158_euler_wallfit_mocG1_r2_ext36k", "run_0159_euler_wallfit_mocG1_r3_ext36k",
    "run_0160_euler_icdep_mocG1_isen_ext36k"] and PARENTS == EV.ARMS["B"] + EV.ARMS["M"] + [EV.ISEN])
check("登録 commit 4dd94af2・V5 の評価器の登録版 sha256", E.PLAN_REG_COMMIT == "4dd94af2"
      and E.V5_EVAL_SHA256 == sha(Path(EV.__file__).resolve()))

# === 接続 ==============================================================================================================
P_ROWS = [{"step": str(s), "a": "1.0"} for s in range(1000, 18001, 1000)]
C_ROWS = [{"step": str(s), "a": "2.0"} for s in range(1000, 36001, 1000)]


def join_err(p, c):
    try:
        E.join_series(p, c, ["a"])
        return None
    except E.JoinError as e:
        return str(e)


g, v = E.join_series(P_ROWS, C_ROWS, ["a"])
check("接続: 通算 1000〜54000 の 54 枚、子は local + 18000", list(g) == list(range(1000, 54001, 1000))
      and v["a"][17] == 1.0 and v["a"][18] == 2.0 and int(g[18]) == 19000)
e = join_err(P_ROWS, C_ROWS[:5] + [C_ROWS[4]] + C_ROWS[5:])
check("接続: 子の重複した step を拒否", e is not None and "重複" in e)
e = join_err(P_ROWS, [{"step": "0", "a": "2.0"}] + C_ROWS)
check("接続: 子の local 0 (親の 18000 と同じ場) を重複として拒否", e is not None and "重複" in e and "18000" in e)
e = join_err(P_ROWS, C_ROWS[:10] + C_ROWS[11:])
check("接続: 子の欠落 (11000) を拒否", e is not None and "欠落" in e and "11000" in e)
e = join_err(P_ROWS[:-1], C_ROWS)
check("接続: 親が 18000 で終わらない (欠落) を拒否", e is not None and "親" in e and "18000" in e)
e = join_err(P_ROWS, C_ROWS[:3] + [C_ROWS[4], C_ROWS[3]] + C_ROWS[5:])
check("接続: 逆順を拒否", e is not None and "逆順" in e)
e = join_err(P_ROWS[:3] + [P_ROWS[4], P_ROWS[3]] + P_ROWS[5:], C_ROWS)
check("接続: 親の逆順を拒否", e is not None and "親" in e and "逆順" in e)
e = join_err(P_ROWS, C_ROWS + [{"step": "37000", "a": "2.0"}])
check("接続: 子の余分な step (37000) を拒否 (固定加算の延長の取り込み)", e is not None and "余分" in e)
e = join_err(P_ROWS, C_ROWS[:-1] + [{"step": "36000", "a": "nan"}])
check("接続: 非有限値を拒否 (0 で埋めない)", e is not None and "非有限" in e)
e = join_err(P_ROWS, C_ROWS[:-1] + [{"step": "36000", "a": ""}])
check("接続: 空欄を拒否", e is not None and "非有限" in e)
e = join_err(P_ROWS, [{"step": r["step"]} for r in C_ROWS])
check("接続: 列の欠損を拒否", e is not None and "列 a が無い" in e)
e = join_err(P_ROWS, [{"step": str(s + 0.5), "a": "2.0"} for s in range(1000, 36001, 1000)])
check("接続: 整数でない step を拒否", e is not None and "整数" in e)
e = join_err([{"step": str(s), "a": "1.0"} for s in range(1000, 30001, 1000)], [])
check("接続: 12000 を足した旧延長の形 (親 1000〜30000) を拒否", e is not None and "親" in e)

# === 周期・転回点 (単体) ================================================================================================
S31 = E.OSC_STEPS
h = 1e-3
sine = [0.1 + 5 * h * math.sin(2 * math.pi * s / 8000.0) for s in S31]
pa = E.period_analysis(S31, sine, h)
check("周期: 周期 8000 の正弦 (31 枚) は転回点 7・完全な周期 3・持続", pa["n_periods"] == 3 and len(pa["turning_points"]) == 7
      and pa["resolved"] and pa["persistent"])
check("周期: 区間の最初の点 (24000) は転回点にしない", all(t["step"] != 24000 for t in pa["turning_points"]))
check("周期: 周期は同じ種類の転回点の間 (26000→34000→42000→50000)",
      [(p_["start_step"], p_["end_step"]) for p_ in pa["periods"]] == [(26000, 34000), (34000, 42000), (42000, 50000)])
check("周期: 周期の平均・半振幅 (両端を含む)", all(abs(p_["half_amp"] - 5 * h) < 1e-12 and abs(p_["mean"] - 0.1) < 1e-12 for p_ in pa["periods"]))
pa = E.period_analysis(S31, [0.1 + 5 * h * math.sin(2 * math.pi * s / 15000.0) for s in S31], h)
check("周期: 周期 15000 (区間に 2 周期) は周期不足で持続にしない", not pa["persistent"] and pa["n_periods"] < 3 and "周期" in pa["reasons"][0])
pa = E.period_analysis(S31, [0.1 + 5 * h * (-1) ** (s // 1000) for s in S31], h)
check("周期: 1 間隔ごとの上下 (出力間隔の折り返し) は解像していない", not pa["persistent"] and not pa["resolved"]
      and any("解像していない" in r for r in pa["reasons"]))
pa = E.period_analysis(S31, [0.1 + 10 * h * math.exp(-(s - 24000) / 20000.0) * math.sin(2 * math.pi * s / 8000.0) for s in S31], h)
check("周期: 振幅が減衰する振動は持続にしない (3 周期あっても半振幅の差 > tol・減衰)", not pa["persistent"] and pa["n_periods"] >= 3
      and any("半振幅の差" in r for r in pa["reasons"]) and any("減衰" in r for r in pa["reasons"]))
pa = E.period_analysis(S31, [0.1 + 5 * h * math.sin(2 * math.pi * s / 8000.0) + 0.3 * h * (s - 24000) / 1000.0 for s in S31], h)
check("周期: 平均が周期ごとに動く振動は持続にしない (平均の差 > tol)", not pa["persistent"] and any("平均の差" in r for r in pa["reasons"]))
pa = E.period_analysis(S31, [0.1 + 0.4 * h * math.sin(2 * math.pi * s / 8000.0) for s in S31], h)
check("周期: 振れが h 未満の振動は転回点にならない (振幅の分解能)", pa["turning_points"] == [] and not pa["persistent"])
pa = E.period_analysis(S31, [0.1] * 31, h)
check("周期: 一定値は転回点なし", pa["turning_points"] == [] and pa["n_periods"] == 0)
amp_g = [0.1 + (5 * h + 0.5 * h * (s - 24000) / 8000.0) * math.sin(2 * math.pi * s / 8000.0) for s in S31]
pa = E.period_analysis(S31, amp_g, h)
check("周期: 振幅が周期あたり 0.5·tol 増えるだけの振動は持続 (減衰しない)", pa["persistent"])
try:
    E.turning_points([0.0, float("nan"), 1.0], h)
    ok_ = False
except ValueError:
    ok_ = True
check("転回点: 非有限の値は例外 (黙って落とさない)", ok_)
ex = E.extrema_record(list(range(1000, 54001, 1000)), [math.sin(2 * math.pi * s / 8000.0) for s in range(1000, 54001, 1000)], 0.1)
check("連続系列の極大・極小を記録 (最大・最小と転回点)", abs(ex["max"]["value"] - 1.0) < 1e-12 and abs(ex["min"]["value"] + 1.0) < 1e-12
      and len(ex["turning_points"]) > 10)

# --- 窓 B の準定常は 13 枚の並びだけで決まる (V5 の関数が step 6000〜18000 の札で判定しても同じ) ---
B13 = {**BASE, "exit_M_dev": 1e-5}
vals13 = {q: np.array([B13[q] + 0.07 * B13[q] * math.sin(0.9 * i) + 0.01 * B13[q] * i for i in range(13)]) for q in EV.QS_COLS}
qv5 = EV.quasisteady_window(vals13)
same = all(qv5[q]["full13"] == classify_series(E.WIN_B_STEPS, list(vals13[q]), 1.0, 0.05, 0.10, 4)[0]
           and qv5[q]["tail5"] == classify_series(E.WIN_B_STEPS, list(vals13[q]), 4.5 / 13, 0.05, 0.10, 4)[0] for q in EV.QS_COLS)
check("窓 B の準定常: V5 の関数の判定 = 実際の step (42000〜54000) での判定", same)


# === 量ごとの分類 (classify_run、合成の連続系列) =========================================================================
def classify(spec, seed=1):
    rng = np.random.default_rng(seed)
    steps = list(range(1000, 54001, 1000))
    rows = [series_values(spec, s, rng) for s in steps]
    vals = {q: np.array([r[q] for r in rows]) for q in EV.QS_COLS}
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "x.csv"
        with open(p, "w", newline="") as f:
            w = csv.writer(f); w.writerow(COLS)
            for s, r in zip(steps[18:], rows[18:]):
                w.writerow([s - 18000] + [r[c] for c in COLS[1:]])
        write_qs(p, Path(td) / "qs.txt")
        return E.classify_run(np.array(steps), vals, (Path(td) / "qs.txt").read_text())


cl = classify({})
check("分類: 減衰する系列は 7 量とも減衰側、run は減衰側", cl["label"] == E.DECAY and all(v == E.DECAY for v in cl["by_label"].values()))
x = cl["quantities"]["overshoot_eta0.1"]
check("分類: 記録 (窓の平均・幅・傾き・末尾 5 枚と全 13 枚の VERDICT・連続系列の極値)",
      all(k in x["window_A"] for k in ("mean", "width", "slope_per_step", "change_over_window", "tail5", "full13"))
      and "file_tail5" in x["window_B"] and "max" in x["extrema_full"] and x["oscillation"]["range"] == [24000, 54000])
check("分類: 符号付き出口 M は M − 6 で記録", abs(cl["quantities"]["exit_core_M"]["window_B"]["mean"] - (-1e-5)) < 1e-7)
osc_q = "overshoot_eta0.1"
cl = classify({osc_q: lambda s: BASE[osc_q] + 0.2 * BASE[osc_q] * math.sin(2 * math.pi * s / 8000.0)})
check("分類: 3 周期以上の定振幅振動 (周期 8000、振幅 20 %) は持続振動側、run は持続振動側を含む",
      cl["by_label"][osc_q] == E.OSC and cl["label"] == E.RUN_OSC
      and sum(v == E.DECAY for v in cl["by_label"].values()) == 6)
cl = classify({osc_q: lambda s: BASE[osc_q] + 0.2 * BASE[osc_q] * math.sin(2 * math.pi * s / 15000.0)})
check("分類: 周期不足 (周期 15000) は判別不能", cl["by_label"][osc_q] == E.UNDET and cl["label"] == E.UNDET)
cl = classify({osc_q: lambda s: BASE[osc_q] + 0.01 * math.exp(-(s - 24000) / 15000.0) * math.sin(2 * math.pi * s / 8000.0)})
xo = cl["quantities"][osc_q]
check("分類: 振幅の減衰が中途半端 (窓 B は未定常、周期の半振幅は tol を超えて減る) は判別不能",
      xo["label"] == E.UNDET and not xo["decaying"] and not xo["persistent"])
cl = classify({osc_q: lambda s: BASE[osc_q] + 0.2 * BASE[osc_q] * (-1) ** (s // 1000)})
check("分類: 1 間隔ごとの上下 (折り返し) は判別不能", cl["by_label"][osc_q] == E.UNDET)
cl = classify({osc_q: lambda s: BASE[osc_q] + (0.0 if s < 39000 else 2 * TOLS[osc_q])})
check("分類: 窓 B は定常でも窓 A → B の平均の移動 > Δq/10 なら判別不能",
      cl["by_label"][osc_q] == E.UNDET and any("平均の移動" in r for r in cl["quantities"][osc_q]["decaying_reasons"]))
cl = classify({"exit_core_M": lambda s: 5.99999 + (0.0 if s < 39000 else 1.5e-5)})
check("分類: 符号付き出口 M の移動 1.5e-5 (> 1e-5) は判別不能 (|M − 6| は 1.8e-5 以内でも)",
      cl["by_label"]["exit_core_M"] == E.UNDET and abs(cl["quantities"]["exit_core_M"]["mean_shift_A_to_B"] - 1.5e-5) < 1e-9)
# 周期 6000 (窓 A と B は 3 周期ずれで同じ位相の並び → 平均の移動 0)、振幅 4 % (窓 B は STEADY、末尾 2 枚は極値でない)
cl = classify({osc_q: lambda s: BASE[osc_q] + 0.04 * BASE[osc_q] * math.sin(2 * math.pi * (s + 500) / 6000.0)})
xb = cl["quantities"][osc_q]
check("分類: 減衰側と持続振動側の両方に当たる量は判別不能 (相対 4 % の振動、窓の平均は動かない)",
      xb["decaying"] and xb["persistent"] and xb["both"] and xb["label"] == E.UNDET)
check("分類: 記録の準定常ファイルが無ければ減衰側にしない (UNKNOWN を通さない)",
      E.classify_run(np.arange(1000, 54001, 1000), {q: np.array([B13[q]] * 54) for q in EV.QS_COLS}, None)["label"] == E.UNDET)

# === 総合 (模擬 case) ==================================================================================================
out = run_eval()
check("陽性の対照: 7 本すべて減衰側 → 窓 B で V5 と同じ比較 (V5b の判定)", out["all_decaying"] and out["window_B_comparison_is_v5b_result"]
      and out["comparison"]["mode"].startswith("V5b の判定") and str(out["overall"]).startswith("許容幅内"))
check("陽性の対照: 前提不成立なし", out["preconditions_failed"] == {})
check("陽性の対照: 6 量が許容幅内、出口コア M は記録のみ、出口較正は据え置き、IC 依存は許容幅内",
      sum(r["verdict"] == "許容幅内" for r in out["comparison"]["rows"]) == 6
      and [r["qty"] for r in out["comparison"]["rows"] if r["verdict"] == "記録のみ"] == ["exit_core_M"]
      and out["comparison"]["exit_calibration"]["verdict"].startswith("据え置き") and out["comparison"]["ic_overall"] == EV.IC_WIN_LABEL)
check("出力に評価器の sha256・登録 commit 4dd94af2・V5 の評価器の sha256", out["evaluator_sha256"] == sha(Path(E.__file__).resolve())
      and out["plan_registration_commit"] == "4dd94af2" and out["v5_evaluator_sha256"] == E.V5_EVAL_SHA256)
check("V5 の結果ファイルは書かない", out["_v5_json_written"] is False)
check("窓 B の比較の値は窓 B の 13 枚 (腕 B の run の平均)", len(out["runs"][CHILDREN[0]]["classification"]["window_B_values"]["overshoot_eta0.1"]) == 13)


def osc_child(i, spec=None):
    def f(root):
        write_child_series(root, CHILDREN[i], spec=spec or {osc_q: lambda s: BASE[osc_q] + 0.2 * BASE[osc_q] * math.sin(2 * math.pi * s / 8000.0)},
                           seed=300 + i)
        refresh_record(root, CHILDREN[i])
    return f


out = run_eval(osc_child(4))
cmp_ = out["comparison"]
check("1 本 (run_0158) が持続振動側 → 窓 B の差は参考値、総合は保留", not out["all_decaying"] and cmp_["mode"].startswith("参考値")
      and all(r["verdict"] == E.REF_LABEL for r in cmp_["rows"] + cmp_["ic_rows"]) and str(out["overall"]).startswith("保留 (V5b")
      and out["labels"][CHILDREN[4]] == E.RUN_OSC and out["window_B_comparison_is_v5b_result"] is False)
check("参考値: 前提成立時の判定を別欄に残す", all("verdict_if_registered" in r for r in cmp_["rows"])
      and cmp_["exit_calibration"]["verdict"] == E.REF_LABEL and cmp_["exit_calibration"]["decision_complete"] is False)
out = run_eval(osc_child(6, {osc_q: lambda s: BASE[osc_q] + 0.2 * BASE[osc_q] * math.sin(2 * math.pi * s / 15000.0)}))
check("ISEN だけ判別不能 → 参考値・保留", out["labels"][E.ISEN] == E.UNDET and out["comparison"]["mode"].startswith("参考値")
      and str(out["overall"]).startswith("保留 (V5b"))


def early_stop(root):
    c = CHILDREN[2]
    write_child_series(root, c, seed=9, last=21000)
    refresh_record(root, c)
    (root / c / E.EARLY_STOP_FILE).write_text("date: x\nreason: 残差に NaN\n")
    (root / c / E.RUN_RC_FILE).write_text("143\n")
    (root / c / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_DIVERGED)


out = run_eval(early_stop)
check("早期停止 (NaN) の run → 保留、総合も保留 (窓 B の値がそろわず参考値も無い)",
      out["labels"][CHILDREN[2]] == E.HOLD and "早期停止" in reasons(out) and str(out["overall"]).startswith("保留")
      and out["comparison"] is None and "参考値も無い" in out["overall"])
check("早期停止: 残差判定 DIVERGED・RUN_RC・接続の欠落を理由に残す", "diverged" in reasons(out) and "RUN_RC" in reasons(out)
      and "欠落" in reasons(out))


def nan_in_series(root):
    write_child_series(root, CHILDREN[5], seed=11, nan_at=30000)
    refresh_record(root, CHILDREN[5])


out = run_eval(nan_in_series)
check("時系列に NaN → 保留 (接続を拒否)", out["labels"][CHILDREN[5]] == E.HOLD and "非有限" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[1] / E.RUN_RC_FILE).unlink())
check("RUN_RC が無い → 保留 (既定で通さない)", out["labels"][CHILDREN[1]] == E.HOLD and "RUN_RC が 0 でない (None)" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[3] / "CONVERGENCE_VERDICT_segment.txt").unlink())
check("子の残差判定が無い → 保留", out["labels"][CHILDREN[3]] == E.HOLD and "missing" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[0] / "solverConfig.yaml").write_text(SOLVER_C.replace("cfl: 2.0,", "cfl: 1.0,")))
check("子の設定が nStepOuter 以外も違う (cfl) → 保留", out["labels"][CHILDREN[0]] == E.HOLD and "nStepOuter" in reasons(out)
      and "deltaT.cfl" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[0] / "solverConfig.yaml").write_text(SOLVER_P))
check("子の nStepOuter が 18000 のまま → 保留", out["labels"][CHILDREN[0]] == E.HOLD and "nStepOuter" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[0] / "solverConfig.yaml").write_text(SOLVER_C.replace("outStepInterval: 1000", "outStepInterval: 2000")))
check("子の出力間隔が 1000 でない → 保留", out["labels"][CHILDREN[0]] == E.HOLD and "outStepInterval" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[0] / "solverConfig.yaml").write_text("# 注記だけ違う\n" + SOLVER_C))
check("子の設定がコメントだけ違う → 減衰側のまま (YAML として同じ)", out["labels"][CHILDREN[0]] == E.DECAY)
out = run_eval(lambda r: (r / CHILDREN[4] / "solverConfig.yaml").write_text(SOLVER_C + "space: {convMethod: 0}\n"))
check("子の設定に重複キー → 保留 (先勝ち・後勝ちで割れる)", out["labels"][CHILDREN[4]] == E.HOLD and "YAML として読めない" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[5] / "bcondConfig.yaml").write_text(BCOND.replace("1600.0", "1500.0")))
check("子の bcond が違う → 保留", out["labels"][CHILDREN[5]] == E.HOLD and "bcondConfig" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[6] / "RUN_PROVENANCE.txt").write_text(PROV.replace("a" * 64, "b" * 64)))
check("子の forge の sha256 が run_0143 と違う → 保留", out["labels"][E.ISEN] == E.HOLD and "forge の sha256" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[6] / "RUN_PROVENANCE.txt").write_text(PROV.replace("BLOCKSIZE=128", "BLOCKSIZE=512")))
check("子の FORGE_CUDA_BLOCKSIZE が親と違う → 保留", out["labels"][E.ISEN] == E.HOLD and "環境変数" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[6] / "RUN_PROVENANCE.txt").write_text(PROV.replace("FORGE_BIN=/x/forge", "FORGE_BIN=/y/forge")))
check("FORGE_BIN の道だけ違う → 減衰側のまま (バイナリは sha256 で照合)", out["labels"][E.ISEN] == E.DECAY)
out = run_eval(lambda r: (r / CHILDREN[6] / "RUN_PROVENANCE.txt").write_text(PROV.replace("'scalarGradient' effective: lsq", "'scalarGradient' effective: gg")))
check("子の起動時の実効値が親と違う → 保留", out["labels"][E.ISEN] == E.HOLD and "実効値" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[1] / "stage_manifest.json").write_text(json.dumps({"stages": [
    {"tag": "soft", "key": MAIN_KEY}, {"tag": "main", "key": MAIN_KEY}]})))
check("子に soft 段がある → 保留", out["labels"][CHILDREN[1]] == E.HOLD and "段の並び" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[1] / "stage_manifest.json").write_text(json.dumps({"stages": [
    {"tag": "main", "key": {**MAIN_KEY, "limiter": "1"}}]})))
check("子の main の hard キーが親と違う → 保留", out["labels"][CHILDREN[1]] == E.HOLD and "hard キー" in reasons(out))
out = run_eval(lambda r: edit_json(r / CHILDREN[2] / E.EXT_RECORD, lambda d: d.update(VERDICT="REFUSED")))
check("延長の記録が OK でない → 保留", out["labels"][CHILDREN[2]] == E.HOLD and "延長の記録" in reasons(out))
out = run_eval(lambda r: edit_json(r / CHILDREN[2] / E.EXT_RECORD, lambda d: d.pop("dry")))
check("延長の記録の dry 欠損 → 保留 (既定で通さない)", out["labels"][CHILDREN[2]] == E.HOLD and "延長の記録" in reasons(out))
out = run_eval(lambda r: edit_json(r / CHILDREN[2] / E.EXT_RECORD, lambda d: d.update(src_res=f"{PARENTS[2]}/res_17000.h5")))
check("延長の restart 元が res_18000 でない → 保留", out["labels"][CHILDREN[2]] == E.HOLD and "res_17000" in reasons(out))
out = run_eval(lambda r: edit_json(r / CHILDREN[3] / "prepare_info.json", lambda d: d.update(DRY=True)))
check("乾式の準備から作った子 → 保留", out["labels"][CHILDREN[3]] == E.HOLD and "乾式" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[3] / E.EXT_RECORD).unlink())
check("延長の記録が無い → 保留", out["labels"][CHILDREN[3]] == E.HOLD and "読めない" in reasons(out))
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d["runs"][CHILDREN[4]].update(dX_F=0.0)))
check("子の評価座標のずれが親と違う → 保留", out["labels"][CHILDREN[4]] == E.HOLD and "親" in reasons(out) and "X_F" in reasons(out))
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d["runs"][CHILDREN[4]].pop("dX_E")))
check("子の評価座標のずれの記録が欠損 → 保留 (0 とみなさない)", out["labels"][CHILDREN[4]] == E.HOLD and "X_E と評価の基準の差 None" in reasons(out))
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d["runs"].pop(CHILDREN[0])))
check("子が時系列の記録に無い → 保留", out["labels"][CHILDREN[0]] == E.HOLD and "時系列の記録" in reasons(out))


def stale_child_csv(root):
    write_child_series(root, CHILDREN[0], seed=999)     # 記録の後に CSV が書き換わった


out = run_eval(stale_child_csv)
check("子の CSV が記録と違う (古い成果物) → 保留", out["labels"][CHILDREN[0]] == E.HOLD and "時系列の記録と違う" in reasons(out))
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d["geometry"].update(X_F=X_F0 + 1e-3)))
check("子の時系列の評価座標が V5 と違う → 共通の前提不成立で保留 (窓 B は参考値)", str(out["overall"]).startswith("保留 (前提不成立: 共通")
      and out["comparison"]["mode"].startswith("参考値 (共通") and out["window_B_comparison_is_v5b_result"] is False)
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d.update(evaluator_sha256="f" * 64)))
check("時系列の評価器が V5 と違う → 共通の前提不成立で保留", str(out["overall"]).startswith("保留 (前提不成立: 共通") and "eval_wallfit_euler" in out["overall"])
out = run_eval(lambda r: edit_json(r / E.SERIES_REC, lambda d: d.update(geom_ref=PARENTS[3])))
check("子の時系列の基準が腕 B の r1 でない → 共通の前提不成立", str(out["overall"]).startswith("保留 (前提不成立: 共通"))
out = run_eval(lambda r: (r / E.SERIES_REC).unlink())
check("子の時系列の記録が無い → 保留", str(out["overall"]).startswith("保留"))
out = run_eval(lambda r: (r / PARENTS[0] / "RUN_PROVENANCE.txt").write_text("forge_bin: x\n"))
check("基準 (run_0143) の forge の sha256 を読めない → 保留", str(out["overall"]).startswith("保留"))


def parent_csv_changed(root):
    write_csv(root / PARENTS[5] / "wallfit_series_v5.csv", range(1000, 18001, 1000), {}, 0, seed=77)


out = run_eval(parent_csv_changed)
check("親の V5 の CSV が V5 の記録と違う → 保留", out["labels"][CHILDREN[5]] == E.HOLD and "V5 の記録と違う" in reasons(out))

# --- V5 の静的な前提 (moc_v5_euler_eval.preconditions をそのまま呼ぶ) ---
out = run_eval(lambda r: (r / PARENTS[6] / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_CONVERGING))
vs = out["v5_static"]
check("V5 の静的な前提: 親の判定区間 (still converging) の理由は使わず記録だけ → V5b の判定は許容幅内",
      str(out["overall"]).startswith("許容幅内") and PARENTS[6] in vs["segment_reasons_not_used"] and vs["failed"] == {})
out = run_eval(lambda r: edit_json(r / PARENTS[4] / "IC_INSPECTION.json", lambda d: d["conditions"].pop(EV.IC_INSPECTION_KEYS[2])))
check("V5 の静的な前提: 親 (腕 M) の IC の検査 C3 の欠け → 比較は保留 (前提不成立)", str(out["overall"]).startswith("保留 (前提不成立)")
      and f"親 {PARENTS[4]}" in out["preconditions_failed"] and EV.IC_INSPECTION_KEYS[2] in reasons(out))
check("V5 の静的な前提の不成立は分類 (減衰側) を変えない", out["all_decaying"] is True)
out = run_eval(lambda r: edit_json(r / PARENTS[3] / "prepare_info.json", lambda d: d["moc"]["gate"].update({"pass": False})))
check("V5 の静的な前提: 親の MOC のゲート不合格 → 保留", str(out["overall"]).startswith("保留 (前提不成立)") and "MOC のゲート" in reasons(out))
_reg = E.V5_EVAL_SHA256
E.V5_EVAL_SHA256 = "0" * 64
try:
    out = run_eval()
finally:
    E.V5_EVAL_SHA256 = _reg
check("V5 の評価器が登録版でない → 保留 (理由の読み分けが保証されない)", str(out["overall"]).startswith("保留 (前提不成立)")
      and "登録版" in reasons(out))
out = run_eval(lambda r: (r / CHILDREN[2] / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_CONVERGING))
check("子の判定区間が still converging → 分類は減衰側のまま、比較は保留 (V5 と同じく pass か plateau)",
      out["labels"][CHILDREN[2]] == E.DECAY and str(out["overall"]).startswith("保留 (前提不成立)") and "converging" in reasons(out))

# --- 判定の向き (窓 B) ---
out = run_eval(lambda r: [osc_child(i, {osc_q: lambda s: BASE[osc_q] + 0.01 + 10.0 * TOLS[osc_q] * math.exp(-s / 4000.0)})(r)
                          for i in (3, 4, 5, 6)])
check("腕 M・ISEN のオーバーシュート +0.01 %pt (減衰側) → 悪化", out["all_decaying"] and str(out["overall"]).startswith("悪化")
      and "overshoot_eta0.1" in out["overall"])
out = run_eval(lambda r: [osc_child(i, {"exit_core_M": lambda s: 5.99999 + 3e-4 + 1e-4 * math.exp(-s / 4000.0)})(r) for i in (3, 4, 5, 6)])
check("腕 M・ISEN の出口コア M +3e-4 → 出口較正はやり直し、|出口 M − 6| は悪化",
      out["comparison"]["exit_calibration"]["verdict"] == "較正をやり直す" and str(out["overall"]).startswith("悪化"))
out = run_eval(lambda r: osc_child(6, {"P_slope_abs_eta0.1": lambda s: BASE["P_slope_abs_eta0.1"] - 0.1 + 10.0 * TOLS["P_slope_abs_eta0.1"] * math.exp(-s / 4000.0)})(r))
check("ISEN が良い向きに大差 → IC 依存あり (両側) で保留", out["comparison"]["ic_overall"].startswith("IC 依存あり")
      and str(out["overall"]).startswith("保留 (IC 依存"))


# === 延長の準備 (prep-child・restart_field・verify-child) を合成の h5 で ===============================================
def make_h5(path: Path, vals: dict, coord):
    import h5py
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=coord)
        f.create_dataset("VIZMESH/CONNE", data=np.arange(12, dtype=np.int32))
        for k, v in vals.items():
            f.create_dataset(f"VALUE/{k}", data=v)


def make_parent_h5(d: Path, n=40, seed=0, drop_in_nozzle=()):
    rng = np.random.default_rng(seed)
    coord = rng.random((n, 3)).astype(np.float32)
    names = ["ro", "roUx", "roUy", "roUz", "roe", "roY0", "roY1", "P", "T", "wall_dist"]
    noz = {k: rng.random(n).astype(np.float32) for k in names if k not in drop_in_nozzle}
    res = {k: rng.random(n).astype(np.float32) for k in names if k != "wall_dist"}
    res["h0"] = rng.random(n).astype(np.float32)
    make_h5(d / "nozzle.h5", noz, coord)
    make_h5(d / "res_18000.h5", res, coord)
    for s in (17000,):
        make_h5(d / f"res_{s}.h5", res, coord)


def mock_parent_full(root: Path, run: str, k: int, drop_in_nozzle=()):
    rec = mock_parent(root, run, k)
    d = root / run
    make_parent_h5(d, seed=k, drop_in_nozzle=drop_in_nozzle)
    (d / "wall_design.csv").write_text("x_m,r_m,theta_rad,M_wall\n0.0,0.1,0.0,1.0\n7.3,0.5,0.0,6.0\n")
    (d / "probe.yaml").write_text("probes: []\n")
    (d / "MESH_QUALITY.txt").write_text("VERDICT: PASS\n")
    (d / "forge_run.log").write_text("log\n")
    (d / "residual_history.csv").write_text("step,rms_ro\n0,1\n")
    return rec


def restart(parent: Path, child: Path) -> int:
    env = {**os.environ, "FORGE_BIN": "/nonexistent/forge", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}
    r = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/restart_field.py"), str(parent / "res_18000.h5"),
                        str(child / "nozzle.h5"), "--dst-run", str(child)], capture_output=True, text=True, env=env)
    (child / "restart_field.log").write_text(r.stdout + r.stderr)
    return r.returncode


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    p0, c0 = root / PARENTS[0], root / CHILDREN[0]
    mock_parent_full(root, PARENTS[0], 0)
    r = E.prep_child(p0, c0, dry=True)
    check("prep-child: 入力だけを複製 (出力・判定ファイル・親の IC の記録は持ち込まない)",
          set(r["copied"]) == {"solverConfig.yaml", "bcondConfig.yaml", "nozzle.h5", "prepare_info.json", "wall_design.csv", "probe.yaml",
                               "MESH_QUALITY.txt"} and not (c0 / "forge_run.log").exists() and not (c0 / "IC_MAP.json").exists()
          and not (c0 / "residual_history.csv").exists() and not list(c0.glob("res_*.h5")))
    why, diff = E.config_failures(SOLVER_P, (c0 / "solverConfig.yaml").read_text())
    check("prep-child: solverConfig は nStepOuter 18000 → 36000 だけ", why == [] and diff == [["time.last.nStepOuter", 18000, 36000]])
    info = json.loads((c0 / "prepare_info.json").read_text())
    check("prep-child: prepare_info に extends・restart 元・DRY の印", info["extends"] == PARENTS[0] and info["DRY"] is True
          and info["restart_from"].startswith(f"{PARENTS[0]}/res_18000.h5") and info["stages"] == "none")
    try:
        E.prep_child(p0, c0, dry=True)
        ok_ = False
    except SystemExit as e_:
        ok_ = "既にある" in str(e_)
    check("prep-child: 子が既にあれば止める", ok_)
    try:
        E.prep_child(p0, root / CHILDREN[1], dry=True)
        ok_ = False
    except SystemExit as e_:
        ok_ = "登録の対応と違う" in str(e_)
    check("prep-child: 親と子の対応が登録と違えば止める", ok_)
    rc = restart(p0, c0)
    rec = E.verify_child(p0, c0, dry=True)
    check("restart_field (本物) → verify-child: OK、移した量の数と保存量 7 個 (roY0・roY1 を含む) を記録",
          rc == 0 and rec["VERDICT"] == "OK" and rec["n_src_conserved"] == 7 and set(rec["src_conserved"]) >= {"roY0", "roY1"}
          and rec["restart_field"]["n_moved"] == 9 and (c0 / E.EXT_RECORD).is_file())
    check("verify-child: 本番 (dry でない) では化学種の属性の継承を要求", E.verify_child(p0, c0, dry=False)["VERDICT"] == "REFUSED")
    (c0 / "solverConfig.yaml").write_text((c0 / "solverConfig.yaml").read_text().replace("cfl: 2.0,", "cfl: 1.0,"))
    rec = E.verify_child(p0, c0, dry=True)
    check("verify-child: 設定が nStepOuter 以外も違えば不成立", rec["VERDICT"] == "REFUSED" and any("cfl" in x for x in rec["failures"]))
    # 保存量の取りこぼし (親の nozzle.h5 に roY1 が無い → restart_field は移さない)
    p1, c1 = root / PARENTS[1], root / CHILDREN[1]
    mock_parent_full(root, PARENTS[1], 1, drop_in_nozzle=("roY1",))
    E.prep_child(p1, c1, dry=True)
    rc = restart(p1, c1)
    rec = E.verify_child(p1, c1, dry=True)
    check("verify-child: restart 元の保存量 roY1 が子に移っていなければ不成立 (skill §3 の化学種の罠)",
          rc == 0 and rec["VERDICT"] == "REFUSED" and rec["conserved_missing_in_child"] == ["roY1"])
    # restart をしていない子 (保存量が親の nozzle.h5 のまま)
    p2, c2 = root / PARENTS[2], root / CHILDREN[2]
    mock_parent_full(root, PARENTS[2], 2)
    E.prep_child(p2, c2, dry=True)
    rec = E.verify_child(p2, c2, dry=True)
    check("verify-child: restart_field を回していない子は不成立 (ビット一致・記録なし)", rec["VERDICT"] == "REFUSED"
          and rec["conserved_not_identical"] and any("restart_field の VERDICT" in x for x in rec["failures"]))
    # 親に 18000 より後の res がある
    p3, c3 = root / PARENTS[3], root / CHILDREN[3]
    mock_parent_full(root, PARENTS[3], 3)
    shutil.copy(p3 / "res_18000.h5", p3 / "res_19000.h5")
    E.prep_child(p3, c3, dry=True)
    restart(p3, c3)
    rec = E.verify_child(p3, c3, dry=True)
    check("verify-child: 親の最後の res が res_18000 でなければ不成立", rec["VERDICT"] == "REFUSED" and rec["parent_last_res"] == "res_19000.h5")
    # 乾式の子は run が拒否する
    try:
        E.run_child(c0)
        ok_ = False
    except SystemExit as e_:
        ok_ = "乾式" in str(e_)
    check("run: 乾式の子は回さない", ok_)
    # 親の nStepOuter が 18000 でない
    p4 = root / PARENTS[4]
    mock_parent_full(root, PARENTS[4], 4)
    (p4 / "solverConfig.yaml").write_text(SOLVER_P.replace("18000", "12000"))
    try:
        E.prep_child(p4, root / CHILDREN[4], dry=True)
        ok_ = False
    except SystemExit as e_:
        ok_ = "nStepOuter" in str(e_) and not (root / CHILDREN[4]).exists()
    check("prep-child: 親の nStepOuter が 18000 でなければ止める (子を作らない)", ok_)


# === 起動スクリプト (DRY と本番の前提で止まること) =========================================================================
def script_case(root: Path) -> Path:
    """合成のリポジトリ配置: <root>/case/x/ に起動スクリプトと評価器 (リンク)、<root>/solver_density_cuda/tools (リンク)、親 7 本。"""
    cd = root / "case/x"
    cd.mkdir(parents=True)
    (root / "solver_density_cuda").mkdir()
    os.symlink(C.parents[1] / "solver_density_cuda/tools", root / "solver_density_cuda/tools")
    for f in ("run_moc_v5b_ext.sh", "moc_v5b_ext_eval.py", "throat_mono_judge.py", "eval_wallfit_euler.py", "moc_v5_euler_eval.py"):
        os.symlink(C / f, cd / f)
    runs5 = {p: mock_parent_full(cd, p, k) for k, p in enumerate(PARENTS)}
    (cd / "_band_ab").mkdir()
    (cd / E.SERIES_REC_V5).write_text(json.dumps({"tag": "v5", "geom_ref": EV.GEOM_REF, "geometry": GEOM, "runs": runs5}))
    return cd


def run_script(cd: Path, env_extra: dict) -> subprocess.CompletedProcess:
    env = {**os.environ, **env_extra}
    return subprocess.run(["bash", str(cd / "run_moc_v5b_ext.sh"), "2"], capture_output=True, text=True, env=env, timeout=600)


with tempfile.TemporaryDirectory() as td:
    cd = script_case(Path(td))
    r = run_script(cd, {"DRY": "1", "FORGE_BIN": "/nonexistent/forge"})
    dry_ok = [json.loads((cd / "_dry_moc_v5b" / c / E.EXT_RECORD).read_text()) for c in CHILDREN
              if (cd / "_dry_moc_v5b" / c / E.EXT_RECORD).is_file()]
    check("起動スクリプト DRY: 合成の親 7 本から準備・restart・検査まで通り、終了コード 0",
          r.returncode == 0 and len(dry_ok) == 7 and all(x["VERDICT"] == "OK" and x["dry"] is True for x in dry_ok)
          and "DRY: 準備・restart・検査まで" in r.stdout)
    check("起動スクリプト DRY: forge を起動しない (RUN_RC・残差・res を作らない)、本番の run 名の dir を作らない",
          not any((cd / "_dry_moc_v5b" / c / f).exists() for c in CHILDREN for f in ("RUN_RC", "residual_history.csv", "stage_manifest.json"))
          and not any((cd / c).exists() for c in CHILDREN) and not list((cd / "_dry_moc_v5b").glob("*/res_*.h5")))
    check("起動スクリプト DRY: 移した量の数を記録 (7 本とも保存量 7・移した量 9)", all(x["n_src_conserved"] == 7 and x["restart_field"]["n_moved"] == 9
                                                                for x in dry_ok))
    check("起動スクリプト DRY: 起動の記録 (評価器の sha256)", (cd / "_band_ab/moc_v5b_launch_dry.txt").is_file()
          and sha(C / "moc_v5b_ext_eval.py") in (cd / "_band_ab/moc_v5b_launch_dry.txt").read_text())
    (cd / CHILDREN[3]).mkdir()
    r = run_script(cd, {"FORGE_BIN": "/nonexistent/forge"})
    check("起動スクリプト 本番: 子が既にあれば止める (終了コード 2、他の子を作らない)", r.returncode == 2 and "既にある" in r.stdout
          and not any((cd / c).exists() for c in CHILDREN if c != CHILDREN[3]))
    (cd / CHILDREN[3]).rmdir()
    fake = Path(td) / "fake_forge"
    fake.write_text("#!/bin/sh\nexit 0\n")
    r = run_script(cd, {"FORGE_BIN": str(fake)})
    check("起動スクリプト 本番: forge の sha256 が run_0143 と違えば止める (子を作らない)", r.returncode == 2
          and "同じバイナリでない" in r.stdout and not any((cd / c).exists() for c in CHILDREN))
    (cd / PARENTS[5] / "RUN_PROVENANCE.txt").write_text(PROV.replace("a" * 64, "c" * 64))
    r = run_script(cd, {"DRY": "1", "FORGE_BIN": "/nonexistent/forge"})
    check("起動スクリプト: 親の forge の sha256 が run_0143 と違えば止める", r.returncode == 2 and PARENTS[5] in r.stdout)


# === 起動スクリプトの早期停止 (nan_in・forge_desc・watch_nan を取り出して、偽の forge で確かめる) ===========================
SH = (C / "run_moc_v5b_ext.sh").read_text()
FUNCS = SH[SH.index("nan_in() {"):SH.index("run_one() {")]


def bash_funcs(body: str, timeout=60) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-c", "set -uo pipefail\nWATCH_SEC=1\n" + FUNCS + "\n" + body], capture_output=True, text=True,
                          timeout=timeout)


HDR = "step,inner,phase,rms_ro,rms_roUx,rms_dq_ro\n"
with tempfile.TemporaryDirectory() as td:
    t = Path(td)
    cases = {"clean": HDR + "0,-1,outer_end,1e-5,2e-3,0\n1,-1,outer_end,9e-6,1e-3,0\n",
             "nan": HDR + "0,-1,outer_end,1e-5,2e-3,0\n1,-1,outer_end,nan,1e-3,0\n",
             "negnan": HDR + "0,-1,outer_end,1e-5,-nan,0\n",
             "inf": HDR + "0,0,inner_end,1e-5,inf,0\n",
             "dqnan": HDR + "0,-1,outer_end,1e-5,2e-3,nan\n"}
    for k, txt in cases.items():
        (t / k).mkdir()
        (t / k / "residual_history.csv").write_text(txt)
    (t / "resnan").mkdir()
    (t / "resnan" / "residual_history.csv").write_text(cases["clean"])
    (t / "resnan" / "res_nan_1200.h5").write_bytes(b"x")
    (t / "none").mkdir()
    r = bash_funcs(" ".join(f'if nan_in {t / k}; then echo "{k}=1"; else echo "{k}=0"; fi;' for k in list(cases) + ["resnan", "none"]))
    got = dict(x.split("=") for x in r.stdout.split())
    check("早期停止の検出: rms_* 列の nan・-nan・inf と res_nan_*.h5 は真、正常・rms_dq_* だけの nan・残差なしは偽",
          got == {"clean": "0", "nan": "1", "negnan": "1", "inf": "1", "dqnan": "0", "resnan": "1", "none": "0"})
    # 偽の forge: 自分の run (nan_flag あり) は 1 秒後に残差へ nan を書き、他の dir の forge は正常のまま。watch_nan は自分の run の forge だけを止める
    bind = t / "bin"
    bind.mkdir()
    fake = bind / "forge"
    fake.write_text("#!/bin/bash\nprintf '" + HDR.replace("\n", "\\n") + "0,-1,outer_end,1e-5,2e-3,0\\n' > residual_history.csv\n"
                    "sleep 1\nif [ -f nan_flag ]; then printf '1,-1,outer_end,nan,1e-3,0\\n' >> residual_history.csv; fi\nsleep 30\n")
    fake.chmod(0o755)
    (t / "mine").mkdir()
    (t / "mine" / "nan_flag").write_text("")
    (t / "other").mkdir()
    body = (f'bash -c \'(cd "{t / "other"}" && exec "{fake}") & cd "{t / "mine"}"; "{fake}"; wait\' > /dev/null 2>&1 &\n'
            "P=$!\nsleep 0.3\n"
            f'watch_nan "{t / "mine"}" "$P"\n'
            'sleep 0.3\nfor k in $(pgrep -P "$P"); do echo "alive $k $(readlink /proc/$k/cwd)"; done\n'
            'for k in $(pgrep -P "$P"); do kill "$k" 2>/dev/null; done; kill "$P" 2>/dev/null; wait 2>/dev/null; echo end\n')
    r = bash_funcs(body, timeout=90)
    es = t / "mine" / "EARLY_STOP.txt"
    alive = [ln for ln in r.stdout.splitlines() if ln.startswith("alive")]
    check("早期停止: 自分の run の forge (cwd 一致・自分の子孫) を止め、EARLY_STOP.txt に記録",
          es.is_file() and "killed: forge pid" in es.read_text() and "実行中に検出" in es.read_text())
    check("早期停止: 同じ親の下でも cwd が違う forge (他の run) は止めない",
          len(alive) == 1 and alive[0].endswith(str(t / "other")) and not (t / "other" / "EARLY_STOP.txt").exists())
    if not (es.is_file() and len(alive) == 1):
        print(r.stdout, r.stderr)

check("V5 の結果ファイル (case の _band_ab/moc_v5_euler_eval.json) は試験の前後で変わらない",
      V5_SHA_BEFORE == (sha(C / "_band_ab/moc_v5_euler_eval.json") if (C / "_band_ab/moc_v5_euler_eval.json").is_file() else None))
check("試験は case の _band_ab に V5b の出力を書かない", not (C / E.OUT_JSON).exists())
print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
