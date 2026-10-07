"""moc_v5_euler_eval.py の試験 (plan discretization-moc-axis-limit-and-corrector §6 V5)。forge・AWS・HDF5 は使わない
(一時ディレクトリの模擬 case に CSV と JSON だけを置き、壁の照合 (nozzle.h5 を読む部分) は差し替える)。
負例: 証拠の欠損・腕の欠損・DIVERGED・乾式 prep・moc のゲート不合格・腕 B の moc が legacy でない・IC 記録の不一致・
古い CSV・未定常・窓の欠け・非有限値・壁の取り違え・ISEN の壁違い・X_F のずれ・基準の取り違え で「保留」になること。
陽性の対照 (許容幅内・据え置き・IC 依存は許容幅内) と、悪化・IC 依存あり (両側)・較正のやり直しの判定も確かめる。
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
NOISE = {"M_wave_eta0.1": 2e-5, "P_wave_eta0.1": 2e-4, "overshoot_eta0.1": 3e-4, "overshoot_exitnorm_eta0.1": 3e-4,
         "P_slope_abs_eta0.1": 1e-3, "exit_core_M": 3e-6, "M_wave_eta0.0": 1e-4}
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


def write_series(d: Path, run: str, offset: dict, seed: int, drop_step=None, nan_at=None):
    rng = np.random.default_rng(seed)
    with open(d / "wallfit_series_v5.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(COLS)
        for st in range(1000, 18001, 1000):
            if st == drop_step:
                continue
            row = {k: BASE[k] + offset.get(k, 0.0) + NOISE[k] * rng.standard_normal() for k in BASE}
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
            (d / "IC_INSPECTION.json").write_text(json.dumps({"VERDICT": "OK", "limit_m": 6.2e-6,
                                                               "conditions": {f"C{n}": {"ok": True} for n in range(1, 6)}}))
        else:
            info["moc"] = {"axis_limit": "analytic", "corrector": "converge", "gate": {"applicable": True, "pass": True}}
            ic = {"VERDICT": "OK", "mode": "isentropic", "dst_sha256_after": "shaI", "dst_mesh_digest": "gridM"}
        info.setdefault("ic", {}).update({"VERDICT": ic["VERDICT"]})
        (d / "prepare_info.json").write_text(json.dumps(info))
        (d / "IC_MAP.json").write_text(json.dumps(ic))
    (root / "_band_ab").mkdir()
    (root / "_band_ab/wallfit_series_v5.json").write_text(json.dumps(
        {"tag": "v5", "geom_ref": EV.GEOM_REF, "geometry": {"X_E": 39.8, "X_F": 95.2}, "runs": runs}))


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


def run_eval(mutate=None, wall=None, full=True, offsets=None):
    with tempfile.TemporaryDirectory() as td:
        root = Path(td); mock_case(root, full=full, offsets=offsets)
        if mutate:
            mutate(root)
        return EV.evaluate(root, wall_check=wall or stub_wall())


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
out = run_eval(lambda r: edit_json(r / M[0] / "IC_INSPECTION.json", lambda d: d["conditions"]["C3"].update(ok=False)))
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
    p = r / M[0] / "QUASISTEADY_wallfit_v5.txt"
    p.write_text(p.read_text().replace("exit_M_dev      : tail mean=0  drift=0.0%/tail  fluct=0.1%   STEADY",
                                       "exit_M_dev      : tail mean=0  drift=9.0%/tail  fluct=9.1%   DRIFTING"))
    refresh_record(r, M[0])


out = run_eval(not_steady)
check("未定常 (exit_M_dev DRIFTING): 保留、理由に STEADY でない量", held(out) and "exit_M_dev DRIFTING" in reasons(out))


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

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
