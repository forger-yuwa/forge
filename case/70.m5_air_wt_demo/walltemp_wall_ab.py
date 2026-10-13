"""case/70: 壁温分布の δ_r が物理壁まで伝わることの A/B (plan tooling-nozzle-wall-temperature-profile §6 V1e、
codex result 段の諮問 2026-10-13 の推奨、forge 0 step)。

腕 A = problem_m5_air_walltemp.yaml (s の表 320/330/420/500/450/380/330 K)、腕 B = 同じでスロートより下流 (s = 0.2・0.5・1) だけ +100 K。
設計壁・scale (spec.r_throat)・補正の設定は同じ。各腕を `deltastar_loop --init-integral --prepare-only` で 1 回ずつ準備し (寸法の逆算は介さない)、
  Δr_ref(x) = δ_r^B(x) − δ_r^A(x)      (積分法の平滑化済み δ_r、delta_r_initial.csv; 5 次スプラインで補間)
  Δr_out(x) = r_phys^B(x) − r_phys^A(x)  (run に保存された物理壁 wall_repr.json の 1 本の B-spline)
を拡大部 x ∈ [1, x_e − 1] r_t の 4001 点で比べる。事前登録: ε = max(1e-8 m, Σ_腕 (物理壁の許容 tol.r + 旧表現との差 d0) × scale)、
判定 A (反映漏れの疑いを除外): max|Δr_ref| > 10 ε かつ max|Δr_out − Δr_ref| ≤ ε。期待差が 10 ε 以下なら判別不能 (合格にしない)。

usage (case dir、design/.venv-opt、FORGE_CONVERTER・FORGE_BIN):  python walltemp_wall_ab.py
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import BSpline, make_interp_spline

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
OUT = C / "_walltemp_wall_ab"
PY = sys.executable
SRC = C / "problem_m5_air_walltemp.yaml"

OUT.mkdir(exist_ok=True)
probs = {"A": SRC, "B": OUT / "problem_m5_air_walltemp_B.yaml"}
txt = SRC.read_text()
old = "Tw_K: [320.0, 330.0, 420.0, 500.0, 450.0, 380.0, 330.0]"
lines = txt.splitlines(keepends=True)
hit = [i for i, l in enumerate(lines) if l.startswith("  wall_thermal:") and old in l]
assert len(hit) == 1, "spec.wall_thermal の行に元の表が見つからない"           # _card の記録の行は書き換えない
lines[hit[0]] = lines[hit[0]].replace(old, "Tw_K: [320.0, 330.0, 420.0, 500.0, 550.0, 480.0, 430.0]")
probs["B"].write_text("".join(lines))

env = dict(os.environ)
env.setdefault("FORGE_CONVERTER", str(ROOT / "solver_density_cuda/.build-native/relwithdebinfo/convertGmshToForge"))
env.setdefault("FORGE_BIN", str(ROOT / "solver_density_cuda/.build-native/relwithdebinfo/forge"))
runs = {}
for arm, pr in probs.items():
    rd = OUT / arm
    if not (rd / "wall_repr.json").exists():
        subprocess.run([PY, "-m", "forge_design.feedback.deltastar_loop", "--problem", str(pr), "--euler-ref", str(C / "run_0003_euler_36k"),
                        "--init-integral", "--run-dir", str(rd), "--prepare-only"], cwd=ROOT / "design", env=env, check=True,
                       stdout=open(OUT / f"prepare_{arm}.log", "w"), stderr=subprocess.STDOUT)
    runs[arm] = rd


def load(rd):
    rep = json.loads((rd / "wall_repr.json").read_text())
    pw = rep["physical_wall"]
    spl = BSpline(np.asarray(pw["t"]), np.asarray(pw["c"]), int(pw["k"]))
    tb = np.loadtxt(rd / "delta_r_initial.csv", delimiter=",", skiprows=1)
    dr = make_interp_spline(tb[:, 0], tb[:, 1], k=5)
    tbc = json.loads((rd / "delta_r_initial.json").read_text())["thermal_bc"]
    return dict(rep=rep, r=spl, dr=dr, scale=float(rep["units"]["scale_m"]), x_e=float(rep["design_wall"]["x_e"]),
                tol_r=float(pw["tol"]["r"]), d0=float(pw["max_err_vs_legacy"]["d0"]), S_c=rep["design_wall"]["S"]["c"], tbc=tbc)


A, B = load(runs["A"]), load(runs["B"])
same_design = (A["S_c"] == B["S_c"]) and (A["scale"] == B["scale"]) and (A["x_e"] == B["x_e"])
x = np.linspace(1.0, A["x_e"] - 1.0, 4001)
d_ref = (B["dr"](x) - A["dr"](x)) * A["scale"]
d_out = (B["r"](x) - A["r"](x)) * A["scale"]
eps = max(1e-8, (A["tol_r"] + A["d0"] + B["tol_r"] + B["d0"]) * A["scale"])
err = np.abs(d_out - d_ref)
k = int(np.argmax(err))
res = {"same_design_wall_and_scale": same_design, "scale_m": A["scale"], "x_range_rt": [float(x[0]), float(x[-1])],
       "Tw_table_A": A["tbc"].get("Tw_table"), "Tw_table_B": B["tbc"].get("Tw_table"),
       "max_abs_dr_ref_m": float(np.max(np.abs(d_ref))), "dr_ref_at_exit_m": float(d_ref[-1]),
       "max_abs_out_minus_ref_m": float(err.max()), "at_x_rt": float(x[k]), "eps_m": eps,
       "ratio_signal_to_eps": float(np.max(np.abs(d_ref)) / eps)}
if not same_design:
    v = "比較不成立 (設計壁か scale が違う)"
elif res["max_abs_dr_ref_m"] <= 10 * eps:
    v = "判別不能 (期待差が 10 ε 以下)"
elif res["max_abs_out_minus_ref_m"] <= eps:
    v = "判定 A: 反映漏れの疑いを除外 (δ_r の差が物理壁に ε 以内で伝わる)"
else:
    v = "判定 B: 期待差があるのに物理壁の差が合わない (完了を保留)"
res["verdict"] = v
(OUT / "walltemp_wall_ab.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float))
print(json.dumps({kk: res[kk] for kk in ("same_design_wall_and_scale", "max_abs_dr_ref_m", "dr_ref_at_exit_m", "max_abs_out_minus_ref_m",
                                         "at_x_rt", "eps_m", "ratio_signal_to_eps", "verdict")}, indent=1, ensure_ascii=False))
