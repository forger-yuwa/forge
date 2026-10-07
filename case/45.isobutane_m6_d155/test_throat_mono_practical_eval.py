"""throat_mono_practical_eval.py の前提検査の負例テスト (codex result 段レビュー 2026-10-07 M2)。
評価量の CSV だけがあり、残差判定・壁の証拠・IC 記録が無い入力で「採用」を返さず「保留 (前提不成立)」になること。
加えて、残差判定が DIVERGED の run が 1 本でもあれば保留になること。forge・AWS は使わない (一時ディレクトリの模擬 case)。
usage: python3 test_throat_mono_practical_eval.py
"""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

C = Path(__file__).resolve().parent
SCRIPT = C / "throat_mono_practical_eval.py"
COLS = ["step", "M_wave_eta0.1", "P_wave_eta0.1", "overshoot_eta0.1", "overshoot_exitnorm_eta0.1", "P_slope_abs_eta0.1", "exit_M_dev", "exit_core_M"]
RUNS = [f"run_0{n}_euler_wallfit_pinG1_r{k}" for n, k in ((140, 1), (141, 2), (142, 3))] + \
       [f"run_0{n}_euler_wallfit_monoG1_r{k}" for n, k in ((143, 1), (144, 2), (145, 3))]
fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def mock_case(root: Path, segment_text=None):
    for r in RUNS + ["run_0146_euler_icab_monoG1_nn"]:
        d = root / r; d.mkdir(parents=True)
        for fname in ("wallfit_series_e3.csv", "wallfit_series_icab.csv"):
            with open(d / fname, "w", newline="") as f:
                w = csv.writer(f); w.writerow(COLS)
                for st in range(1000, 18001, 1000):
                    w.writerow([st, 0.0064, 0.041, 0.038, 0.038, 0.19, 2e-5, 5.99999])
        if segment_text is not None:
            (d / "CONVERGENCE_VERDICT_segment.txt").write_text(segment_text)
    shutil.copy(C / "throat_mono_practical_eval.py", root / "throat_mono_practical_eval.py")


def run(root: Path) -> dict:
    q = subprocess.run([sys.executable, str(SCRIPT), str(root)], capture_output=True, text=True, cwd=str(C))
    f = root / "_band_ab/throat_mono_practical_eval.json"
    return json.loads(f.read_text()) if f.is_file() else {"rc": q.returncode, "stderr": q.stderr[-500:]}


with tempfile.TemporaryDirectory() as td:
    root = Path(td); mock_case(root)
    out = run(root)
    check("CSV だけ: 総合が「保留 (前提不成立)」", str(out.get("overall", "")).startswith("保留 (前提不成立)"))
    check("CSV だけ: 「採用」を返さない", "採用" not in str(out.get("overall", "")) or str(out.get("overall", "")).startswith("保留"))
    check("CSV だけ: 6 本すべてが前提不成立に挙がる", len(out.get("preconditions_failed", {})) == 6)
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    mock_case(root, segment_text="[segment] 判定区間 main\n=== x  [last step 17999]  -> DIVERGED (nan) ===\n")
    out = run(root)
    reasons = " ".join(" ".join(v) for v in out.get("preconditions_failed", {}).values())
    check("残差 DIVERGED: 保留", str(out.get("overall", "")).startswith("保留 (前提不成立)"))
    check("残差 DIVERGED: 理由に残差判定が入る", "残差判定" in reasons)


# --- 2 回目の result 段レビュー M1: 同じ非単調壁を両腕に使った負例を最終判定まで通す (壁の証拠は差し替え)。陽性の対照も置く ---
import importlib
import numpy as np
from scipy.interpolate import make_interp_spline


def spline_dict(kind):
    x = np.linspace(-1.0, 3.0, 801)
    r = 1 + x ** 2 / 4 - x ** 3 / 30 if kind == "mono" else 1 + x ** 2 / 4 + 1e-3 * np.exp(-((x - 0.1) / 0.05) ** 2)
    s_ = make_interp_spline(x, r, k=5)
    return {"t": [float(v) for v in s_.t], "c": [float(v) for v in s_.c], "k": 5}


def full_mock(root: Path, b_kind: str):
    mock_case(root, segment_text="[segment] 判定区間 main\n=== x  [last step 17999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===\n")
    for r in RUNS:
        d = root / r; arm = "A" if "pinG1" in r else "B"
        wf = {"mono_r2": None if arm == "A" else [0.0, 1.5], "spline": spline_dict("bump" if arm == "A" else b_kind)}
        (d / "prepare_info.json").write_text(json.dumps({"wall_fit": wf, "scale_m": 0.0766539}))
        if arm == "A":
            (d / "restart_field.log").write_text("VERDICT: OK (9 量を移した、SRC とビット一致)\n")
        else:
            (d / "IC_MAP.json").write_text(json.dumps({"VERDICT": "OK", "mode": "index", "dst_sha256_after": "same"}))


def run_inproc(root: Path, n_disc: int, n_own: int):
    sys.argv = [str(SCRIPT), str(root)]
    sys.path.insert(0, str(C))
    import throat_mono_ab
    orig = throat_mono_ab.wall_evidence
    throat_mono_ab.wall_evidence = lambda rd, other=None: {"status": "consistent", "vs_other": {"status": "ok", "n_discriminable": n_disc,
                                                                                               "n_discriminable_matching_own": n_own}}
    try:
        mod = importlib.import_module("throat_mono_practical_eval"); mod = importlib.reload(mod)
        mod.main()
    finally:
        throat_mono_ab.wall_evidence = orig
    return json.loads((root / "_band_ab/throat_mono_practical_eval.json").read_text())


with tempfile.TemporaryDirectory() as td:
    root = Path(td); full_mock(root, "bump")
    out = run_inproc(root, 0, 0)
    reasons = " ".join(" ".join(v) for v in out.get("preconditions_failed", {}).values())
    check("同じ非単調壁を両腕に: 総合が保留", str(out.get("overall", "")).startswith("保留 (前提不成立)"))
    check("同じ非単調壁を両腕に: 理由に「見分けられる壁節点が無い」", "見分けられる壁節点が無い" in reasons)
    check("同じ非単調壁を両腕に: 理由に「腕 B の保存 spline が単調でない」", "単調でない" in reasons)
with tempfile.TemporaryDirectory() as td:
    root = Path(td); full_mock(root, "mono")
    out = run_inproc(root, 12, 7)
    check("壁の取り違え (見分けられる節点の一部が他腕に一致): 保留", str(out.get("overall", "")).startswith("保留 (前提不成立)"))
with tempfile.TemporaryDirectory() as td:
    root = Path(td); full_mock(root, "mono")
    out = run_inproc(root, 12, 12)
    check("陽性の対照 (A 山あり・B 単調・見分けられる節点がすべて自腕): 採用", str(out.get("overall", "")).startswith("単調壁を候補形状として採用"))
print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
