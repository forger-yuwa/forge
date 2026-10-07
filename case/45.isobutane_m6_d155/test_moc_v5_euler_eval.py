"""moc_v5_euler_eval.py の試験 (plan discretization-moc-axis-limit-and-corrector §6 V5)。forge・AWS・HDF5 は使わない
(一時ディレクトリの模擬 case に CSV と JSON だけを置き、壁の照合 (nozzle.h5 を読む部分) と評価座標の感度 (res_*.h5 を読む部分) は
差し替える)。V5 の run の結果は読まない (合成した系列だけ)。
負例: 証拠の欠損・腕の欠損・DIVERGED・乾式 prep・moc のゲート不合格・腕 B の moc が legacy でない・IC 記録の不一致・
古い CSV・未定常・窓の欠け・非有限値・壁の取り違え・ISEN の壁違い・X_F のずれ・基準の取り違え で「保留」になること。
陽性の対照 (許容幅内・据え置き・IC 依存は許容幅内) と、悪化・IC 依存あり (両側)・較正のやり直しの判定も確かめる。
§6 V5 の追加登録 (2026-10-07、諮問 notes/reviews/2026-10-07-moc-v5-eval-interpretation-diagnose.md) の反例:
  (a) 前半 8 枚 0.03・後半 5 枚 0.04 → 末尾 5 枚 STEADY・全 13 枚 DRIFTING で保留 / (b) 近零の出口誤差が DRIFTING でも絶対許容内で
  比較に使う / (c) 出口 M 5.9998 と 6.0002 は出口誤差の差 0 でも IC の出口 M 条件で不合格 / (d) C1〜C5 の欠け → 保留 /
  (e) 腕 B だけ非 STEADY でも保留 (対称)。
usage: python3 test_moc_v5_euler_eval.py
"""
import copy
import csv
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import make_interp_spline

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import moc_v5_euler_eval as EV  # noqa: E402

B, M, I = EV.ARMS["B"], EV.ARMS["M"], EV.ISEN
COLS = ["step", "M_wave_eta0.1", "P_wave_eta0.1", "overshoot_eta0.1", "overshoot_exitnorm_eta0.1", "P_slope_abs_eta0.1",
        "exit_M_dev", "exit_core_M", "M_wave_eta0.0"]
BASE = {"M_wave_eta0.1": 0.0064, "P_wave_eta0.1": 0.041, "overshoot_eta0.1": 0.038, "overshoot_exitnorm_eta0.1": 0.038,
        "P_slope_abs_eta0.1": 0.19, "exit_core_M": 5.99999, "M_wave_eta0.0": 0.02}
# exit_core_M の雑音は、近零の出口誤差 (M − 6 ≈ −1e-5) が窓で OSCILLATING になりつつ、絶対許容内 (1.8e-5) に余裕を持って入る大きさ
NOISE = {"M_wave_eta0.1": 2e-5, "P_wave_eta0.1": 2e-4, "overshoot_eta0.1": 3e-4, "overshoot_exitnorm_eta0.1": 3e-4,
         "P_slope_abs_eta0.1": 1e-3, "exit_core_M": 2e-6, "M_wave_eta0.0": 1e-4}
X_E0, X_F0 = 39.82004263, 95.22667765            # 時系列の記録の共通の評価座標 (run_series と同じ形: WIN_T = (X_E + 2, X_F − 1) ほか)
GEOM = {"X_E": X_E0, "X_F": X_F0, "WIN_T": [X_E0 + 2.0, X_F0 - 1.0], "WIN_O": [X_E0 - 15.0, X_F0]}
SOLVER_CFG = ("mesh: {discretization: \"node\", isAxisymmetric: 1, meshFileName: \"nozzle.h5\"}\nsolver: \"SLAU\"\n"
              "time:\n  last: {nStepOuter: 18000}\n  deltaT: {cfl: 2.0, cfl_pseudo: 2.0, blockDPLUR: 1, implicitRelax: 0.7}\n"
              "space: {convMethod: 1, limiter: 2}\n")
BCOND_CFG = "inlet:  {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 5500000.0, Tt: 1600.0}}\n"
PROV = ("forge_sha256: " + "a" * 64 + "\n'slauWallNormalChi' effective: 1 (auto)\n'scalarGradient' effective: lsq (default)\n")
SEG_PLATEAU = "[segment] 判定区間 main\n=== x  [last step 17999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===\n"
SEG_DIVERGED = "[segment] 判定区間 main\n=== x  [last step 17999]  -> DIVERGED (nan) ===\n"
fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def spline_dict(kind="mono", shift=0.0):
    x = np.linspace(-1.0, 3.0, 801)
    r = (1 + x ** 2 / 4 - x ** 3 / 30 + shift) if kind == "mono" else 1 + x ** 2 / 4 + 1e-3 * np.exp(-((x - 0.1) / 0.05) ** 2)
    s_ = make_interp_spline(x, r, k=5)
    return {"t": [float(v) for v in s_.t], "c": [float(v) for v in s_.c], "k": 5}


def write_series(d: Path, run: str, offset: dict, seed: int, drop_step=None, nan_at=None, fn=None, last_step=18000):
    """fn: {列: step → 値} (雑音なしで上書き。exit_core_M を上書きすれば exit_M_dev も |M − 6| で作り直す)。"""
    rng = np.random.default_rng(seed)
    fn = fn or {}
    with open(d / "wallfit_series_v5.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(COLS)
        for st in range(1000, last_step + 1, 1000):
            if st == drop_step:
                continue
            row = {k: BASE[k] + offset.get(k, 0.0) + NOISE[k] * rng.standard_normal() for k in BASE}
            row.update({k: g(st) for k, g in fn.items()})
            row["exit_M_dev"] = abs(row["exit_core_M"] - 6.0)
            vals = [st] + [row[c] for c in COLS[1:]]
            if nan_at is not None and st == nan_at:
                vals[1] = float("nan")
            w.writerow(vals)
    (d / "QUASISTEADY_wallfit_v5.txt").write_text(
        "=== series  -> STEADY ===\n" + "".join(f"  {c:16s}: tail mean=0  drift=0.0%/tail  fluct=0.1%   STEADY\n" for c in COLS[1:]))


def mock_case(root: Path, full=True, offsets=None):
    """7 本の模擬 run。full=False なら評価量の CSV (と時系列の記録) だけ。offsets: {run: {量: 加算}}。"""
    offsets = offsets or {}
    runs = {}
    for k, run in enumerate(B + M + [I]):
        d = root / run; d.mkdir(parents=True)
        write_series(d, run, offsets.get(run, {}), seed=10 + k)
        runs[run] = {"status": "ok", "reason": None, "dX_E": 0.0, "dX_F": (0.0 if run in B else 9.9e-9), "n_snaps": 18,
                     "last_step": 18000, "csv_sha256": sha(d / "wallfit_series_v5.csv"),
                     "quasisteady_sha256": sha(d / "QUASISTEADY_wallfit_v5.txt")}
        if not full:
            continue
        (d / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_PLATEAU)
        (d / "solverConfig.yaml").write_text(SOLVER_CFG)
        (d / "bcondConfig.yaml").write_text(BCOND_CFG)
        (d / "RUN_PROVENANCE.txt").write_text(PROV)
        (d / "stage_manifest.json").write_text(json.dumps({"manifest_version": 2, "stages": [
            {"tag": "soft", "key": {"convMethod": "0", "bcond_sha1": "abc"}}, {"tag": "main", "key": {"convMethod": "1", "bcond_sha1": "abc"}}]}))
        role = "B" if run in B else ("M" if run in M else "ISEN")
        info = {"scale_m": 0.0768075, "wall_fit": {"mono_r2": [0.0, 1.5], "spline": spline_dict("mono", 0.0 if role == "B" else 1e-5)}}
        src = f"/aws/case/{EV.IC_RUN_NAME}/{EV.IC_RES}"
        if role == "B":
            ic = {"VERDICT": "OK", "mode": "index", "src_res": src, "dst_sha256_after": "shaB", "dst_mesh_digest": "gridB"}
        elif role == "M":
            info["moc"] = {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
            info["ic"] = {"VERDICT": "OK", "inspection": {"limit_m": 6.2e-6}}
            ic = {"VERDICT": "OK", "mode": "index", "src_res": src, "dst_sha256_after": "shaM", "dst_mesh_digest": "gridM",
                  "checks": {"displacement": {"detail": {"limit_m": 6.2e-6, "max_m": 6.19e-6}}}}
            (d / "IC_INSPECTION.json").write_text(json.dumps({"VERDICT": "OK", "limit_m": 6.2e-6,      # 条件名は moc_v5_euler.inspect_ic と同じ
                                                               "conditions": {k: {"ok": True} for k in EV.IC_INSPECTION_KEYS}}))
        else:
            info["moc"] = {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
            ic = {"VERDICT": "OK", "mode": "isentropic", "dst_sha256_after": "shaI", "dst_mesh_digest": "gridM"}
        info.setdefault("ic", {}).update({"VERDICT": ic["VERDICT"]})
        (d / "prepare_info.json").write_text(json.dumps(info))
        (d / "IC_MAP.json").write_text(json.dumps(ic))
    (root / "_band_ab").mkdir()
    (root / "_band_ab/wallfit_series_v5.json").write_text(json.dumps(
        {"tag": "v5", "geom_ref": EV.GEOM_REF, "geometry": GEOM, "runs": runs}))


def refresh_record(root: Path, run: str):
    """CSV・準定常ファイルを書き換えた後に時系列の記録の sha256 を合わせる (古い成果物の検査とは別の負例を作るため)。"""
    p = root / "_band_ab/wallfit_series_v5.json"
    r = json.loads(p.read_text())
    r["runs"][run]["csv_sha256"] = sha(root / run / "wallfit_series_v5.csv")
    r["runs"][run]["quasisteady_sha256"] = sha(root / run / "QUASISTEADY_wallfit_v5.txt")
    p.write_text(json.dumps(r))


def edit_json(p: Path, fn):
    d = json.loads(p.read_text()); fn(d); p.write_text(json.dumps(d))


def stub_wall(n_disc=12, n_own=12, own="consistent"):
    def f(case, run, other):
        out = {"own": {"status": own}}
        if other is not None:
            out["vs_other"] = {"status": "ok", "n_discriminable": n_disc, "n_discriminable_matching_own": n_own}
        return out
    return f


def stub_coord(diff=1e-7):
    """評価座標の感度の差し替え (res_*.h5 を読まない)。各ずれで M_wave だけ diff 動く。"""
    def f(case, run, geom, shifts):
        return {"step": 18000, "base": {c: 0.0 for c in EV.QS_COLS},
                "diff": {k: {c: (diff if c == "M_wave_eta0.1" else 0.0) for c in EV.QS_COLS} for k in shifts}}
    return f


def run_eval(mutate=None, wall=None, full=True, offsets=None, coord=None):
    with tempfile.TemporaryDirectory() as td:
        root = Path(td); mock_case(root, full=full, offsets=offsets)
        if mutate:
            mutate(root)
        return EV.evaluate(root, wall_check=wall or stub_wall(), coord_sensitivity=coord or stub_coord())


def reasons(out):
    return " ".join(" ".join(v) for v in out.get("preconditions_failed", {}).values())


def held(out):
    return str(out.get("overall", "")).startswith("保留")


# --- 判定関数の単体 ---
check("片側: D + 2SE ≤ Δq → 許容幅内", EV.judge_one_sided(0.0005, 0.0002, 0.001) == "許容幅内")
check("片側: D − 2SE ≥ Δq → 悪化", EV.judge_one_sided(0.002, 0.0004, 0.001) == "悪化")
check("片側: 境界域 → 保留", EV.judge_one_sided(0.0009, 0.0002, 0.001) == "保留")
check("片側: 大きな負の D (改善) は許容幅内", EV.judge_one_sided(-0.01, 0.0001, 0.001) == "許容幅内")
check("両側: 大きな負の D は IC 依存あり", EV.judge_two_sided(-0.01, 0.0001, 0.001) == "IC 依存あり")
check("両側: |D| + 2SE ≤ Δq → 許容幅内", EV.judge_two_sided(-0.0003, 0.0002, 0.001) == "IC 依存は許容幅内")
check("両側: 境界域 → 保留", EV.judge_two_sided(-0.0009, 0.0002, 0.001) == "保留")
check("較正: |D| + 2SE ≤ 1e-4 → 据え置き", EV.judge_exit_calibration(2e-5, 3e-5).startswith("据え置き"))
check("較正: |D| − 2SE > 1e-4 → やり直す", EV.judge_exit_calibration(-3e-4, 5e-5) == "較正をやり直す")
check("較正: 判別不能 → 保留", EV.judge_exit_calibration(9e-5, 2.5e-5).startswith("保留"))
st1 = EV.arm_stats([np.array([1.0, 2.0, 3.0])])
check("腕の統計: 1 本なら SE は時間変動の SE だけ", st1["se_between_runs"] is None and abs(st1["se"] - 1.0 / np.sqrt(3)) < 1e-12)
st3 = EV.arm_stats([np.full(13, 1.0), np.full(13, 2.0), np.full(13, 3.0)])
check("腕の統計: 3 本で時間変動 0 なら SE = sd(m_i)/√3", abs(st3["se"] - 1.0 / np.sqrt(3)) < 1e-12)

# --- 陽性の対照 ---
out = run_eval()
check("陽性の対照: 総合が許容幅内", str(out["overall"]).startswith("許容幅内"))
check("陽性の対照: 前提不成立なし", out["preconditions_failed"] == {})
check("陽性の対照: 腕 B の moc 無しを例外として記録 (3 本)", sum("legacy・fixed2 とみなした" in e for e in out["precondition_exceptions"]) == 3)
check("陽性の対照: 出口較正は据え置き", out["exit_calibration"]["verdict"].startswith("据え置き"))
check("陽性の対照: IC 依存は許容幅内", out["ic_dependence"]["overall"] == "IC 依存は許容幅内")
check("陽性の対照: 6 量が判定、1 量が記録のみ", sum(r["verdict"] == "許容幅内" for r in out["rows"]) == 6 and
      [r["qty"] for r in out["rows"] if r["verdict"] == "記録のみ"] == ["exit_core_M"])

# --- 負例: 保留になること ---
out = run_eval(full=False)
check("CSV だけ: 保留 (前提不成立)", str(out["overall"]).startswith("保留 (前提不成立)"))
check("CSV だけ: 7 本すべてが前提不成立に挙がる", all(r in out["preconditions_failed"] for r in B + M + [I]))
check("CSV だけ: 採用・許容幅内を返さない", "許容幅内" not in str(out["overall"]).split(":")[0])
out = run_eval(lambda r: shutil.rmtree(r / M[2]))
check("腕の欠損 (run_0152 なし): 保留", held(out) and "run dir が無い" in reasons(out) and M[2] in out["missing"])
out = run_eval(lambda r: shutil.rmtree(r / I))
check("ISEN の欠損: 保留", held(out) and I in out["preconditions_failed"])
out = run_eval(lambda r: (r / M[1] / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_DIVERGED))
check("残差 DIVERGED: 保留、理由に残差判定", held(out) and "残差判定が diverged" in reasons(out))
out = run_eval(lambda r: (r / B[0] / "CONVERGENCE_VERDICT_segment.txt").unlink())
check("残差判定ファイルの欠損 (腕 B): 保留", held(out) and "残差判定が missing" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "prepare_info.json", lambda d: d.update(DRY=True)))
check("乾式 prep: 保留、理由に DRY", held(out) and "乾式確認" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[1] / "prepare_info.json", lambda d: d["moc"]["gate"].update({"pass": False})))
check("moc のゲート不合格: 保留", held(out) and "MOC のゲートが合格でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[2] / "prepare_info.json", lambda d: d["moc"]["gate"].pop("pass")))
check("moc のゲートの pass 欠損: 保留 (既定で通さない)", held(out) and "MOC のゲートが合格でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "prepare_info.json", lambda d: d["moc"].update(corrector="fixed2")))
check("腕 M の moc が converge でない: 保留", held(out) and "analytic・converge でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "prepare_info.json", lambda d: d.pop("moc")))
check("腕 M の moc 欠損: 保留", held(out) and "analytic・converge でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / B[1] / "prepare_info.json",
                                   lambda d: d.update(moc={"axis_limit": "analytic", "corrector": "converge"})))
check("腕 B の moc が legacy でない: 保留", held(out) and "legacy・fixed2 でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[2] / "IC_MAP.json", lambda d: d.update(dst_sha256_after="other")))
check("IC 記録の不一致 (腕 M の写像後 sha256): 保留", held(out) and "sha256 がそろわない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[1] / "IC_INSPECTION.json", lambda d: d.update(limit_m=1e-5)))
check("IC 記録の不一致 (検査の上限と写像の上限): 保留", held(out) and "一致しない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "IC_INSPECTION.json", lambda d: d["conditions"][EV.IC_INSPECTION_KEYS[2]].update(ok=False)))
check("IC の検査の条件が不合格: 保留", held(out) and "IC_INSPECTION.json" in reasons(out))
out = run_eval(lambda r: (r / M[0] / "IC_INSPECTION.json").unlink())
check("IC の検査の記録が無い: 保留", held(out) and "証拠を読めない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[1] / "IC_MAP.json", lambda d: d["checks"]["displacement"]["detail"].update(max_m=7e-6)))
check("写像の移動が上限を超える記録: 保留", held(out) and "上限" in reasons(out))
out = run_eval(lambda r: edit_json(r / B[2] / "IC_MAP.json", lambda d: d.update(src_res="/x/run_0062_euler_wallfit_fit_r1_ext6k/res_6000.h5")))
check("IC の保存場の取り違え (腕 B): 保留", held(out) and "でない (/x/run_0062" in reasons(out))
out = run_eval(lambda r: edit_json(r / I / "IC_MAP.json", lambda d: d.update(VERDICT="DRY (化学種の属性なし)")))
check("ISEN の IC が OK でない: 保留", held(out) and "isentropic であること" in reasons(out))
out = run_eval(lambda r: edit_json(r / I / "IC_MAP.json", lambda d: d.update(dst_mesh_digest="gridX")))
check("ISEN の格子が腕 M と違う: 保留", held(out) and "格子ハッシュ" in reasons(out))
out = run_eval(lambda r: edit_json(r / I / "prepare_info.json", lambda d: d["wall_fit"].update(spline=spline_dict("mono", 2e-5))))
check("ISEN の壁が腕 M と違う: 保留", held(out) and "完全一致しない" in reasons(out))


def stale_csv(r):
    write_series(r / M[0], M[0], {}, seed=77)       # 記録の後に CSV が書き換わった (別の呼び出し・古い成果物)


out = run_eval(stale_csv)
check("古い CSV (記録の sha256 と違う): 保留", held(out) and "時系列の記録と違う" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"].pop(M[1])))
check("時系列の記録に無い run: 保留", held(out) and "記録 (wallfit_series_v5.json) に無い" in reasons(out))
out = run_eval(lambda r: (r / "_band_ab/wallfit_series_v5.json").unlink())
check("時系列の記録が無い: 保留", held(out) and "時系列の記録" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"][M[2]].update(dX_F=1e-5)))
check("X_F のずれが 1e-6 r_t を超える: 保留", held(out) and "X_F と評価の基準の差" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"][M[2]].pop("dX_E")))
check("X_E のずれの記録が欠損: 保留 (0 とみなさない)", held(out) and "X_E と評価の基準の差 None" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d.update(geom_ref=M[0])))
check("評価の基準が腕 B の r1 でない: 保留", held(out) and "腕 B の r1" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"][B[0]].update(status="missing")))
check("時系列の評価が missing: 保留", held(out) and "時系列の評価が 'missing'" in reasons(out))


def not_steady(r):
    # 追加登録 (2026-10-07) で exit_M_dev には絶対許容内の例外ができたので、出口誤差を許容 1.8e-5 の外 (M − 6 ≈ +2.9e-4) にして
    # 例外に当たらない形で確かめる (旧版は例外を適用しない規則だった)
    write_series(r / M[0], M[0], {"exit_core_M": 3e-4}, seed=13)
    p = r / M[0] / "QUASISTEADY_wallfit_v5.txt"
    p.write_text(p.read_text().replace("exit_M_dev      : tail mean=0  drift=0.0%/tail  fluct=0.1%   STEADY",
                                       "exit_M_dev      : tail mean=0  drift=9.0%/tail  fluct=9.1%   DRIFTING"))
    refresh_record(r, M[0])


out = run_eval(not_steady)
check("未定常 (exit_M_dev DRIFTING、絶対許容内の例外に当たらない): 保留、理由に STEADY でない量",
      held(out) and "exit_M_dev DRIFTING" in reasons(out))


def missing_qs_line(r):
    p = r / I / "QUASISTEADY_wallfit_v5.txt"
    p.write_text("\n".join(l for l in p.read_text().splitlines() if "P_wave_eta0.1" not in l) + "\n")
    refresh_record(r, I)


out = run_eval(missing_qs_line)
check("準定常の行の欠損: 保留 (UNKNOWN を通さない)", held(out) and "P_wave_eta0.1 UNKNOWN" in reasons(out))


def drop_step(r):
    write_series(r / M[1], M[1], {}, seed=99, drop_step=12000)
    refresh_record(r, M[1])


out = run_eval(drop_step)
check("窓の欠け (12000 なし): 保留 (欠損)", str(out["overall"]).startswith("保留 (欠損)") and M[1] in out["missing"])


def nan_value(r):
    write_series(r / B[1], B[1], {}, seed=98, nan_at=9000)
    refresh_record(r, B[1])


out = run_eval(nan_value)
check("非有限値: 保留 (欠損)", str(out["overall"]).startswith("保留 (欠損)") and B[1] in out["missing"])
out = run_eval(wall=stub_wall(12, 7))
check("壁の取り違え (見分けられる節点の一部が腕 B に一致): 保留", held(out) and "壁の取り違え" in reasons(out))
out = run_eval(wall=stub_wall(0, 0))
check("腕 B と見分けられない壁: 保留", held(out) and "見分けられる壁節点が無い" in reasons(out))
out = run_eval(wall=stub_wall(own="mismatch"))
check("壁節点が自腕の spline と不一致: 保留", held(out) and "自腕の当てはめ後 spline に一致しない" in reasons(out))


def bad_wall(case, run, other):
    raise OSError("nozzle.h5 が読めない")


out = run_eval(wall=bad_wall)
check("壁の証拠が読めない: 保留", held(out) and "nozzle.h5 が読めない" in reasons(out))
out = run_eval(lambda r: edit_json(r / B[0] / "prepare_info.json", lambda d: d["wall_fit"].update(spline=spline_dict("bump"))))
check("腕 B の保存 spline が単調でない: 保留", held(out) and "単調でない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "prepare_info.json", lambda d: d["wall_fit"].update(mono_r2=None)))
check("腕 M の mono_r2 が無効: 保留", held(out) and "mono_r2" in reasons(out))

out = run_eval(lambda r: (r / M[1] / "RUN_PROVENANCE.txt").write_text(PROV.replace("a" * 64, "b" * 64)))
check("forge のバイナリが腕 B と違う: 保留", held(out) and "実効設定が" in reasons(out) and "forge_sha256" in reasons(out))
out = run_eval(lambda r: (r / I / "solverConfig.yaml").write_text(SOLVER_CFG.replace("cfl: 2.0", "cfl: 1.0")))
check("本段の設定が腕 B と違う (cfl): 保留", held(out) and "solverConfig" in reasons(out))
out = run_eval(lambda r: (r / M[0] / "solverConfig.yaml").write_text(SOLVER_CFG.replace("{cfl: 2.0", "{cfl: 2e+0")))
check("数値の書式違い (2e+0 は PyYAML で文字列): 保留 (安全側)", held(out) and "solverConfig" in reasons(out))
out = run_eval(lambda r: (r / M[0] / "solverConfig.yaml").write_text("# コメントだけ違う\n" + SOLVER_CFG))
check("コメントだけ違う: 許容幅内 (YAML として同じ)", str(out["overall"]).startswith("許容幅内"))
out = run_eval(lambda r: (r / B[2] / "stage_manifest.json").unlink())
check("段の記録が無い: 保留", held(out) and "実効設定を読めない" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[2] / "stage_manifest.json", lambda d: d["stages"].insert(1, {"tag": "mid", "key": {}})))
check("段の並びが違う (mid あり): 保留", held(out) and "soft → main でない" in reasons(out))
out = run_eval(lambda r: (r / M[2] / "RUN_PROVENANCE.txt").write_text("forge_bin: x\n"))
check("RUN_PROVENANCE に sha256 が無い: 保留", held(out) and "forge_sha256 を読めない" in reasons(out))

# --- 判定の向き ---
off_M = {r: {"overshoot_eta0.1": 0.01} for r in M + [I]}
out = run_eval(offsets=off_M)
check("腕 M のオーバーシュート +0.01 %pt: 悪化", str(out["overall"]).startswith("悪化") and "overshoot_eta0.1" in out["overall"])
off_I = {I: {"P_slope_abs_eta0.1": -0.1}}
out = run_eval(offsets=off_I)
check("ISEN が良い向きに大差: IC 依存あり (両側) → 総合は保留", out["ic_dependence"]["overall"].startswith("IC 依存あり")
      and str(out["overall"]).startswith("保留 (IC 依存"))
check("IC 依存の前提が不成立なら出口較正も保留", out["exit_calibration"]["verdict"].startswith("保留"))
off_X = {r: {"exit_core_M": 3e-4} for r in M + [I]}
out = run_eval(offsets=off_X)
check("出口コア M +3e-4: 較正をやり直す", out["exit_calibration"]["verdict"] == "較正をやり直す")
check("出口コア M +3e-4: |出口 M − 6| が悪化", str(out["overall"]).startswith("悪化") and "exit_M_dev" in out["overall"])

# === §6 V5 の追加登録 (2026-10-07) ===================================================================================
WS = EV.WIN_STEPS


def qsw_all(v="STEADY"):
    return {c: {"tail5": v, "full13": v} for c in EV.QS_COLS}


def rewrite(r, run, seed, fn=None, offset=None, file_verdicts=None, last_step=18000):
    """run の CSV を書き直し (fn: 列 → step の関数)、記録の準定常ファイルの行を file_verdicts で置き換え、時系列の記録を合わせる。"""
    write_series(r / run, run, offset or {}, seed=seed, fn=fn, last_step=last_step)
    p = r / run / "QUASISTEADY_wallfit_v5.txt"
    t = p.read_text()
    for c, v in (file_verdicts or {}).items():
        t = t.replace(f"  {c:16s}: tail mean=0  drift=0.0%/tail  fluct=0.1%   STEADY", f"  {c:16s}: tail mean=0  drift=9.0%/tail  fluct=1.0%   {v}")
    p.write_text(t)
    refresh_record(r, run)


SEED = {run: 10 + k for k, run in enumerate(B + M + [I])}

# --- 準定常の判定 (単体) ---
why, stt = EV.judge_window_quasisteady(qsw_all(), {c: "STEADY" for c in EV.QS_COLS}, {"applies": False})
check("準定常: 3 つとも STEADY → 理由なし・全量 STEADY", why == [] and all(v == "STEADY" for v in stt.values()))
q = qsw_all(); q["exit_M_dev"] = {"tail5": "DRIFTING", "full13": "OSCILLATING"}
why, stt = EV.judge_window_quasisteady(q, {c: "STEADY" for c in EV.QS_COLS}, {"applies": True})
check("準定常: exit_M_dev が非 STEADY でも例外が成り立てば「絶対許容内」(STEADY と表示しない)",
      why == [] and stt["exit_M_dev"] == "絶対許容内" and stt["exit_M_dev"] != "STEADY")
why, _ = EV.judge_window_quasisteady(q, {c: "STEADY" for c in EV.QS_COLS}, {"applies": False, "reason": "x"})
check("準定常: exit_M_dev が非 STEADY で例外に当たらなければ保留", why != [] and "exit_M_dev DRIFTING (窓の末尾 5 枚)" in why[0])
q = qsw_all(); q["exit_core_M"] = {"tail5": "STEADY", "full13": "OSCILLATING"}
why, stt = EV.judge_window_quasisteady(q, {c: "STEADY" for c in EV.QS_COLS}, {"applies": True})
check("準定常: 例外は exit_M_dev だけ (exit_core_M の OSCILLATING は例外が成り立っても保留)",
      why != [] and "exit_core_M OSCILLATING (窓の全 13 枚)" in why[0] and stt["exit_core_M"] == "保留")
q = qsw_all(); q["P_wave_eta0.1"] = {"tail5": "STEADY", "full13": "DRIFTING"}
why, _ = EV.judge_window_quasisteady(q, {c: "STEADY" for c in EV.QS_COLS}, {"applies": True})
check("準定常: 末尾 5 枚 STEADY でも全 13 枚 DRIFTING なら保留", why != [] and "P_wave_eta0.1 DRIFTING (窓の全 13 枚)" in why[0])
why, _ = EV.judge_window_quasisteady(qsw_all(), {c: "STEADY" for c in EV.QS_COLS if c != "M_wave_eta0.1"}, {"applies": True})
check("準定常: 記録の準定常ファイルの行の欠損 (UNKNOWN) は保留", why != [] and "M_wave_eta0.1 UNKNOWN (記録の準定常ファイル)" in why[0])

# --- (a) 反例: 前半 8 枚 0.03・後半 5 枚 0.04 ---
step_series = [0.03] * 8 + [0.04] * 5
qa = EV.quasisteady_window({c: np.array(step_series) for c in EV.QS_COLS})["overshoot_eta0.1"]
check("(a) 単体: 前半 8 枚 0.03・後半 5 枚 0.04 は末尾 5 枚 STEADY、全 13 枚 DRIFTING",
      qa["tail5"] == "STEADY" and qa["full13"] == "DRIFTING")
sa = EV.arm_stats([np.array(step_series)])
check("(a) 単体: 比較相手を全窓平均にすると 2SE = 0.00280883 < Δq 0.003 (末尾だけ見ると IC 合格に化ける)",
      abs(2 * sa["se"] - 0.00280883) < 1e-8 and abs(sa["mean"] - 0.44 / 13) < 1e-15)


def isen_step(r):
    rewrite(r, I, SEED[I], fn={"overshoot_eta0.1": lambda s: 0.03 if s <= 13000 else 0.04})
    for run in M:
        rewrite(r, run, SEED[run], fn={"overshoot_eta0.1": lambda s: 0.44 / 13})


out = run_eval(isen_step)
wq = out["run_meta"][I]["window_quasisteady"]["window"]["overshoot_eta0.1"]
irow = next(x for x in out["ic_dependence"]["rows"] if x["qty"] == "overshoot_eta0.1")
check("(a) ISEN: 末尾 5 枚 STEADY・全 13 枚 DRIFTING", wq["tail5"] == "STEADY" and wq["full13"] == "DRIFTING")
check("(a) 保留 (前提不成立)、理由は ISEN の全 13 枚だけ", held(out) and list(out["preconditions_failed"]) == [I]
      and "overshoot_eta0.1 DRIFTING (窓の全 13 枚)" in reasons(out) and "(窓の末尾 5 枚)" not in reasons(out))
check("(a) IC の差は参考値 (末尾だけの検査なら「許容幅内」だった)", irow["verdict"] == EV.IC_REF_LABEL
      and irow["verdict_if_preconditions_held"] == "IC 依存は許容幅内" and out["ic_dependence"]["overall"].startswith("保留 (前提未達"))
check("(a) 出口較正も保留 (V5 の較正判断は未完了)", out["exit_calibration"]["verdict"].startswith("保留")
      and out["exit_calibration"]["decision_complete"] is False)

# --- (b) 近零の出口誤差: 末尾 10 枚の最大絶対値 5.1e-6・幅 3.5e-6 で減衰、DRIFTING でも絶対許容内 ---
def sig_b(s):
    return 1.4e-6 + 3.696e-6 * 0.7 ** ((s - 9000) / 1000.0)


sb = np.array([sig_b(s) for s in WS])
eb = EV.exit_exception(sb)
check("(b) 単体: 末尾 10 枚の最大絶対値 5.1e-6・幅 3.5e-6",
      abs(eb["e4_tail"]["max_abs"] - 5.096e-6) < 1e-9 and abs(eb["e4_tail"]["width"] - 3.546e-6) < 5e-9)
check("(b) 単体: E4 は絶対許容内、全 13 枚の最大絶対値・幅 ≤ 1.8e-5 → 例外が成り立つ",
      eb["e4_tail"]["status"] == "絶対許容内" and eb["full_max_abs"] <= 1.8e-5 and eb["full_width"] <= 1.8e-5 and eb["applies"] is True)
qb = EV.quasisteady_window({**{c: np.full(13, 1.0) for c in EV.QS_COLS}, "exit_M_dev": np.abs(sb)})["exit_M_dev"]
check("(b) 単体: |M − 6| は窓で DRIFTING (相対正規化のため)", qb["tail5"] == "DRIFTING" and qb["full13"] == "DRIFTING")


def near_zero_all(r):
    for run in B + M + [I]:                  # 全 7 本に対称に (腕 B・ISEN も同じ扱い)
        rewrite(r, run, SEED[run], fn={"exit_core_M": lambda s: 6.0 + sig_b(s)}, file_verdicts={"exit_M_dev": "DRIFTING"})


out = run_eval(near_zero_all)
st_b = out["quasisteady"]["status"]
check("(b) 全 7 本で exit_M_dev が DRIFTING でも絶対許容内 → 前提成立・許容幅内",
      out["preconditions_failed"] == {} and str(out["overall"]).startswith("許容幅内"))
check("(b) 表示は「絶対許容内」(STEADY と書かない)", all(st_b[r]["exit_M_dev"] == "絶対許容内" for r in B + M + [I]))
wb = out["run_meta"][B[0]]["window_quasisteady"]
check("(b) 元の VERDICT を残す (窓の末尾 5 枚・全 13 枚 DRIFTING、記録のファイル DRIFTING)",
      wb["window"]["exit_M_dev"]["tail5"] == "DRIFTING" and wb["window"]["exit_M_dev"]["full13"] == "DRIFTING"
      and wb["file_tail5"]["exit_M_dev"] == "DRIFTING")
check("(b) exit_M_dev は比較に使われる (行に判定がある)", next(x for x in out["rows"] if x["qty"] == "exit_M_dev")["verdict"] == "許容幅内")


def near_zero_big_head(r):
    near_zero_all(r)
    rewrite(r, B[2], SEED[B[2]], fn={"exit_core_M": lambda s: 6.0 + (2.5e-5 if s == 6000 else sig_b(s))},
            file_verdicts={"exit_M_dev": "DRIFTING"})


out = run_eval(near_zero_big_head)
eb2 = out["run_meta"][B[2]]["window_quasisteady"]["exit_M_dev_exception"]
check("(b′) 末尾 10 枚は絶対許容内でも、全 13 枚の最大絶対値 2.5e-5 > 1.8e-5 なら例外に当たらず保留",
      eb2["e4_tail"]["status"] == "絶対許容内" and eb2["applies"] is False and held(out) and list(out["preconditions_failed"]) == [B[2]]
      and "全 13 枚の最大絶対値" in reasons(out))

# --- (c) 出口 M 5.9998 と 6.0002: 出口誤差の差 0 でも IC の出口 M 条件で不合格 ---
out = run_eval(offsets={**{r: {"exit_core_M": -0.00019} for r in M}, I: {"exit_core_M": 0.00021}})
ic = {x["qty"]: x for x in out["ic_dependence"]["rows"]}
check("(c) 前提は成立 (腕 M・ISEN の出口誤差 2e-4 は STEADY)", out["preconditions_failed"] == {})
check("(c) |出口 M − 6| の IC の差は 0 付近で許容幅内", ic["exit_M_dev"]["verdict"] == "IC 依存は許容幅内"
      and abs(ic["exit_M_dev"]["D_ISEN_minus_M"]) < 1e-5)
check("(c) 出口コア M そのもの (符号つき) の差 4e-4 で IC 依存あり (許容 1e-4)", ic["exit_core_M"]["verdict"] == "IC 依存あり"
      and abs(ic["exit_core_M"]["D_ISEN_minus_M"] - 4e-4) < 1e-5 and ic["exit_core_M"]["dq"] == 1e-4)
check("(c) 総合は保留 (IC 依存の前提)", out["ic_dependence"]["overall"].startswith("IC 依存あり") and "exit_core_M" in out["ic_dependence"]["overall"]
      and str(out["overall"]).startswith("保留 (IC 依存"))
check("(c) 両側: 出口 M の差が逆向き (ISEN が小さい) でも同じく IC 依存あり",
      run_eval(offsets={**{r: {"exit_core_M": 0.00021} for r in M}, I: {"exit_core_M": -0.00019}})["ic_dependence"]["overall"].startswith("IC 依存あり"))
out = run_eval()
check("陽性の対照: IC 依存の行に出口コア M (Δ 1e-4) が入り、許容幅内",
      next(x for x in out["ic_dependence"]["rows"] if x["qty"] == "exit_core_M")["verdict"] == "IC 依存は許容幅内")

# --- (d) C1〜C5 の 5 キーを個別に要求 ---
C3K, C5K = EV.IC_INSPECTION_KEYS[2], EV.IC_INSPECTION_KEYS[4]
out = run_eval(lambda r: edit_json(r / M[1] / "IC_INSPECTION.json", lambda d: d["conditions"].pop(C3K)))
check("(d) C3 のキーが欠けたら保留 (残りが全部合格でも)", held(out) and f"{C3K} が無い" in reasons(out)
      and "OK・全条件合格でない" not in reasons(out))
out = run_eval(lambda r: edit_json(r / M[2] / "IC_INSPECTION.json", lambda d: d["conditions"][C5K].update(ok="true")))
check("(d) C5 の ok が true でない (文字列): 保留", held(out) and f"{C5K} が合格" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "IC_INSPECTION.json", lambda d: d.pop("conditions")))
check("(d) 条件の記録が丸ごと無い: 保留", held(out) and "C1_ic_index_map_checks" in reasons(out))
out = run_eval(lambda r: edit_json(r / M[0] / "IC_INSPECTION.json",
                                   lambda d: d.update(conditions={f"C{n}": {"ok": True} for n in range(1, 6)})))
check("(d) 条件名が inspect_ic と違う (C1〜C5 の短い名前) は保留", held(out) and "C2_columns_x_unchanged" in reasons(out))

# --- (e) 腕 B だけ非 STEADY でも保留 (対称) ---
out = run_eval(lambda r: rewrite(r, B[1], SEED[B[1]], fn={"P_wave_eta0.1": lambda s: 0.041 * (1 + 0.08 * (-1) ** (s // 1000))}))
check("(e) 腕 B の 1 本の P 波が OSCILLATING (窓の末尾 5 枚・全 13 枚): 保留", held(out) and list(out["preconditions_failed"]) == [B[1]]
      and "P_wave_eta0.1 OSCILLATING (窓の末尾 5 枚)" in reasons(out) and "P_wave_eta0.1 OSCILLATING (窓の全 13 枚)" in reasons(out))
out = run_eval(lambda r: rewrite(r, B[0], SEED[B[0]], fn={"exit_core_M": lambda s: 6.0 + 2.5e-5 * (1 + 0.3 * (-1) ** (s // 1000))}))
check("(e) 腕 B の出口誤差が許容 1.8e-5 の外で OSCILLATING: 例外に当たらず保留", held(out) and list(out["preconditions_failed"]) == [B[0]]
      and "exit_M_dev OSCILLATING" in reasons(out))

# --- 窓は固定 (延長分・窓の外を使わない) ---
out = run_eval(lambda r: rewrite(r, I, SEED[I], fn={"overshoot_eta0.1": lambda s: 0.5 if s > 18000 else 0.038}, last_step=19000))
ic = {x["qty"]: x for x in out["ic_dependence"]["rows"]}
check("窓は 6000〜18000 に固定: 19000 の枚は使わない (ISEN の窓平均 0.038、許容幅内)",
      abs(ic["overshoot_eta0.1"]["ISEN"]["mean"] - 0.038) < 1e-15 and str(out["overall"]).startswith("許容幅内")
      and out["window_fixed"] is True)

# --- 前提未達のとき IC の差は参考値 ---
out = run_eval(lambda r: (r / M[1] / "CONVERGENCE_VERDICT_segment.txt").write_text(SEG_DIVERGED))
check("前提未達: IC の差は参考値 (行に前提成立時の判定を残す)、総合の IC は保留",
      out["ic_dependence"]["overall"].startswith("保留 (前提未達") and all(x["verdict"] == EV.IC_REF_LABEL for x in out["ic_dependence"]["rows"])
      and all("verdict_if_preconditions_held" in x for x in out["ic_dependence"]["rows"]))

# --- 出口較正の文言 ---
out = run_eval()
ec = out["exit_calibration"]
check("出口較正: 「据え置き」に出口 M = 6 の達成の証明でないことを明記", ec["verdict"].startswith("据え置き")
      and "出口 M = 6 の達成の証明ではない" in ec["verdict"] and ec["note"] == EV.EXIT_CAL_NOTE and ec["decision_complete"] is True)
check("出口較正: 判別不能の文言に「較正判断は未完了」", "未完了" in EV.judge_exit_calibration(9e-5, 2.5e-5))

# --- 評価座標 ---
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"][M[2]].update(dX_F=1e-5)))
check("評価座標: 1e-6 r_t は「座標の整合の許容差」と表示", held(out) and "座標の整合の許容差" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["geometry"].pop("WIN_T")))
check("評価座標: 共通の選択範囲 (WIN_T) の記録が無い: 保留", held(out) and "WIN_T" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d.pop("geometry")))
check("評価座標: 共通の評価座標の記録が無い: 保留", held(out) and "共通の評価座標 (geometry) が無い" in reasons(out))
out = run_eval(lambda r: edit_json(r / "_band_ab/wallfit_series_v5.json", lambda d: d["runs"][B[0]].update(dX_E=1e-9)))
check("評価座標: 基準 (腕 B の r1) 自身の座標でない: 保留", held(out) and "自身の座標でない" in reasons(out))
out = run_eval()
co = out["evaluation_coordinates"]
sv = co["sensitivity"][M[0]]
check("評価座標の感度: 記録だけ (判定に使わない)、表示は座標の整合の許容差",
      co["used_in_judgement"] is False and co["tolerance_label"] == "座標の整合の許容差" and co["tolerance_r_t"] == 1e-6)
check("評価座標の感度: 自身のずれと ±1e-6 r_t (X_E・X_F) の 5 通りで差を取る",
      set(sv["diff"]) == {"own", "+tol_X_E", "-tol_X_E", "+tol_X_F", "-tol_X_F"})
check("評価座標の感度: 量ごとの最大絶対値と Δq に対する比", abs(sv["max_abs_diff"]["M_wave_eta0.1"] - 1e-7) < 1e-20
      and abs(sv["max_abs_diff_over_dq"]["M_wave_eta0.1"] - 1e-4) < 1e-15)
smp = co["samples"]
check("評価座標の標本数: 共通の座標で数え、腕 B (ずれ 0) は共通と同じ",
      smp["common"]["n_fit"] > 0 and smp["common"]["n_WIN_T"] > 0 and all(smp["own"][r]["same_as_common"] for r in B)
      and smp["tol_shifts_same_as_common"] is True)


def coord_boom(case, run, geom, shifts):
    raise OSError("res_18000.h5 が読めない")


out = run_eval(coord=coord_boom)
check("評価座標の感度が取れなくても判定は変わらない (記録にエラー)", str(out["overall"]).startswith("許容幅内")
      and "res_18000.h5 が読めない" in out["evaluation_coordinates"]["sensitivity"][I]["error"])
out = run_eval(coord=stub_coord(diff=1.0))
check("評価座標の感度が大きくても判定には使わない (記録だけ)", str(out["overall"]).startswith("許容幅内"))



# 既定の感度 (default_coord_sensitivity) を合成した場で: eval_wallfit_euler.quantities をソースから取り出して呼べること
# (quantities は np.trapezoid を使うので numpy 2 が要る。numpy 1 では例外になり、記録用なので判定には影響しない)
def synth_field(case, res=None):
    x = np.linspace(0.0, 96.0, 1921); eta = np.linspace(0.0, 1.0, 41)
    X = np.repeat(x[:, None], len(eta), axis=1); R = (1.0 + 0.02 * X) * eta[None, :]
    Mf = 6.0 * (1.0 + 1e-3 * np.sin(0.7 * X) * (1.0 - 0.5 * eta[None, :]))
    Pf = 2237.0 * (1.0 + 1e-2 * np.cos(0.3 * X) + 1e-4 * X)
    return {"X": X, "R": R, "V": {"M": Mf, "P": Pf}}


sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.report.nozzle_report as NR  # noqa: E402
_lf = NR.load_field
NR.load_field = synth_field
try:
    shifts = {"zero": (0.0, 0.0), "+tol_X_F": (0.0, 1e-6), "own_like": (0.0, 9.9e-9)}
    if hasattr(np, "trapezoid"):
        sd = EV.default_coord_sensitivity(Path("/nonexistent"), "run_x", GEOM, shifts)
        check("既定の感度 (合成した場、numpy 2): ずれ 0 の差は厳密に 0、有限、7 量そろう",
              all(v == 0.0 for v in sd["diff"]["zero"].values()) and set(sd["base"]) == set(EV.QS_COLS)
              and all(np.isfinite(v) for d in sd["diff"].values() for v in d.values()))
        check("既定の感度: 出口コア M は評価座標に依存しない (差 0)", sd["diff"]["+tol_X_F"]["exit_core_M"] == 0.0)
    else:
        try:
            EV.default_coord_sensitivity(Path("/nonexistent"), "run_x", GEOM, shifts)
            ok_ = False
        except AttributeError:
            ok_ = True
        check("既定の感度 (numpy 1): np.trapezoid が無く例外 (記録用なので判定に影響しないことは上の試験で確認)", ok_)
        check("既定の感度 (numpy 1): 取り出した quantities は eval_wallfit_euler と同じ定義 (grid が呼べる)",
              EV.sample_sets(EV._wallfit_defs(), GEOM)["n_fit"] > 0)
finally:
    NR.load_field = _lf

# --- 出力の記録 ---
out = run_eval()
check("出力に評価器自身の sha256", out["evaluator_sha256"] == sha(Path(EV.__file__).resolve()) and out["evaluator"] == "moc_v5_euler_eval.py")
check("出力に plan の追加登録の commit (9a09c0a5)", out["plan_registration_commit"] == "9a09c0a5")
check("出力に準定常の条件 (drift 0.05・osc 0.10・末尾 5 枚と全 13 枚)", out["quasisteady"]["drift"] == 0.05
      and out["quasisteady"]["osc"] == 0.10 and out["quasisteady"]["tail"] == 5 and out["quasisteady"]["full"] == 13)

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
