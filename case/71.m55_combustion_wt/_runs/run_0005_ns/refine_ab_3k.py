# plan tooling-nozzle-card-local-mesh §5.1 #3k: 抽出の径方向の標本化の A/B (forge 0 step、事前登録 2026-10-15、codex diagnose
# notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md)。変えるのは band_local_deficit の refine だけ (A = 4、B = 32)。
# 固定: run_0003_ns_cont/res_200000.h5 (pass 1 の壁の元)、Euler 参照 run_0001_euler/res_12000.h5 (sha256 93cb88b0…)、一定 x への補間
# (remap_constant_x、remap_edge missing)、各列の保存済み band_y_b (run_0003_ns_cont/delta_r_equiv.csv; 帯の梯子は選び直さない)、
# Euler の補間 (deltastar_from_core_matched_euler の euler_q_at と同じ式)、平滑化 (quintic P-spline ノット 2 r_t・λ 1・重み hard_ok)。
# 判定 (x 1〜3.5 m): max|Δδ_E| ≤ 0.01 mm かつ max|Δθ| ≤ 0.01° (θ = atan(r′_design + δ′_E)) → 結果 B、どちらかを超える → 結果 A。
# 拡張 (2 係数の前提; codex diagnose 2026-10-15 m55-two-coef-contur ⑤): 端点の |Δδ_E| ≤ 0.01 mm・x ≥ 0 の |Δ(c·δ_C + δ₀)| ≤ 0.01 mm・2 係数の壁角 |Δθ| ≤ 0.01°。
# 使い方 (AWS A): PYTHONPATH=~/forge-cardmesh/design python refine_ab_3k.py <case71 dir> <out dir>
import json
import sys
from pathlib import Path

import numpy as np

from forge_design.geometry.wall_axismach import load_wall_file
from forge_design.metrics.deltastar import (_load_structured, band_local_deficit, euler_cover_range, resolve_euler_ref,
                                            smooth_delta_quintic)

C = Path(sys.argv[1]); O = Path(sys.argv[2]); O.mkdir(parents=True, exist_ok=True)
NS = C / "run_0003_ns_cont"
info = json.loads((NS / "prepare_info.json").read_text())
er = info["euler_ref"] if isinstance(info.get("euler_ref"), dict) else None
if er is None:                                                 # prepare_info が run dir の文字列のときは停止規則の記録の固定の参照
    er = json.loads((NS / "stop_rule.json").read_text())["settings"].get("euler_ref")
if er is None or not str(er.get("sha256", "")).startswith("93cb88b0"):
    raise SystemExit(f"run_0003_ns_cont の euler_ref が想定と違う: {info.get('euler_ref')}")
EREF = {"run": str(C / "run_0001_euler"), "res": er["res"], "sha256": er["sha256"]}   # 固定の参照 (sha256 は prepare_info の完全な値)

saved = np.genfromtxt(NS / "delta_r_equiv.csv", delimiter=",", names=True)
N = _load_structured(NS, remap_constant_x=True, remap_edge="missing", res="res_200000.h5")
eref = resolve_euler_ref(EREF)
E = _load_structured(eref["run"], res=eref["res"])
xE = E["x"][:, 0]; rwE_col = E["r"][:, -1]
cover = euler_cover_range(xE, edge_tol=True)
miss = set(N.get("remap_missing_cols", {}).get("i", []))


def euler_q_at(x, r):
    k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
    w = float(np.clip((x - xE[k]) / max(xE[k + 1] - xE[k], 1e-30), 0.0, 1.0))
    q0 = np.interp(r, E["r"][k], E["q"][k]); q1 = np.interp(r, E["r"][k + 1], E["q"][k + 1])
    return (1 - w) * q0 + w * q1, (1 - w) * rwE_col[k] + w * rwE_col[k + 1]


xs_saved = np.asarray(saved["x_rt"]); yb_saved = np.asarray(saved["band_y_b"])
cols = []
for i in range(N["x"].shape[0]):
    x = float(N["x"][i, 0])
    if not (cover["lo"] <= x <= cover["hi"]) or i in miss:
        continue
    cols.append(i)
xN = np.array([float(N["x"][i, 0]) for i in cols])
if len(xN) != len(xs_saved) or np.max(np.abs(xN - xs_saved)) > 1e-12:
    raise SystemExit(f"保存済みの抽出と列が合わない (列 {len(xN)} / 保存 {len(xs_saved)})")

out = {}
for R in (4, 32):
    raw = np.full(len(cols), np.nan); hard = np.zeros(len(cols), bool); ybu = np.full(len(cols), np.nan)
    for k, i in enumerate(cols):
        x = float(N["x"][i, 0]); r0 = N["r"][i]; q0 = N["q"][i]
        _, rwE = euler_q_at(x, r0[:1])
        res = band_local_deficit(r0, q0, lambda rr, x=x: euler_q_at(x, rr)[0], float(rwE), k_band=1.5, refine=R,
                                 y_b_fixed=float(yb_saved[k]))
        if res is None:
            continue
        raw[k] = res["delta_r"]; ybu[k] = res.get("band_y_b", np.nan)
        hard[k] = np.isfinite(res["delta_r"]) and not any(h in res["reason"] for h in ("no_root", "nan", "nonpos_ref", "nonmono_zone"))
    f_s, diag = smooth_delta_quintic(xN, raw, weights=(hard & np.isfinite(raw)).astype(float), knot_spacing=2.0, lam=1.0, positive=False)
    out[R] = {"raw": raw, "hard": hard, "dE": f_s(xN), "yb": ybu, "diag": diag}

S = float(info["scale_m"]); mm = S * 1e3
W = load_wall_file(C / "run_0002_ns")                          # pass 0 の設計壁 (run_0003_ns_cont は run_0002_ns と同じ壁)
rp = np.asarray(W["design"].r(xN, 1), dtype=float)
th = {R: np.degrees(np.arctan(rp + np.gradient(out[R]["dE"], xN))) for R in out}
rep = {"saved_raw_vs_A": {"max_abs_mm": float(np.nanmax(np.abs(out[4]["raw"] - saved["delta_r_raw"])) * mm),
                          "n_nan_mismatch": int(np.sum(np.isfinite(out[4]["raw"]) != np.isfinite(saved["delta_r_raw"]))),
                          "bit_identical": bool(np.array_equal(out[4]["raw"], np.asarray(saved["delta_r_raw"]), equal_nan=True))},
       "yb_used_equals_saved": {R: bool(np.nanmax(np.abs(out[R]["yb"] - yb_saved)) < 1e-12) for R in out}}
xm = xN * S
for nm, (a, b) in {"judge_1_3.5m": (1.0, 3.5), "all_0_8m": (0.0, 8.01)}.items():
    m = (xm >= a) & (xm <= b)
    d_raw = (out[32]["raw"] - out[4]["raw"]) * mm; d_E = (out[32]["dE"] - out[4]["dE"]) * mm; d_th = th[32] - th[4]
    st4 = np.abs(np.diff(out[4]["raw"][m])) * mm; st32 = np.abs(np.diff(out[32]["raw"][m])) * mm
    i1 = int(np.nanargmax(np.abs(d_E[m]))); i2 = int(np.nanargmax(np.abs(d_th[m])))
    rep[nm] = {"max_abs_dE_mm": float(np.abs(d_E[m][i1])), "x_m_at": float(xm[m][i1]),
               "max_abs_dtheta_deg": float(np.abs(d_th[m][i2])), "x_m_at_theta": float(xm[m][i2]),
               "max_abs_draw_mm": float(np.nanmax(np.abs(d_raw[m]))),
               "raw_adjacent_step_max_mm": {"A4": float(np.nanmax(st4)), "B32": float(np.nanmax(st32))},
               "raw_adjacent_step_p99_mm": {"A4": float(np.nanpercentile(st4, 99)), "B32": float(np.nanpercentile(st32, 99))},
               "n_hard_fail": {"A4": int((~out[4]["hard"][m]).sum()), "B32": int((~out[32]["hard"][m]).sum())}}
j = rep["judge_1_3.5m"]
rep["tolerance"] = {"max_abs_dE_mm": 0.01, "max_abs_dtheta_deg": 0.01}
rep["result"] = "B (両方以内)" if (j["max_abs_dE_mm"] <= 0.01 and j["max_abs_dtheta_deg"] <= 0.01) else "A (どちらかを超える)"
# --- 拡張 (codex diagnose 2026-10-15 m55-two-coef-contur ⑤): 2 係数の較正 (plan §4.3) の前提 ---------------------------------
# δ_C = pass 0 の物理壁の δ (run_0003_ns_cont/delta_r_next.csv の delta_in_prev、r_t 単位)、設計壁は固定。各腕の平滑化後の δ_E の
# x = 0 と出口 (最後の列) で c・δ₀ を解き、2 係数の補正 c·δ_C + δ₀ と壁角 atan(r′_design + c·δ_C′) を比べる
from forge_design.feedback.deltastar_loop import read_delta_r_next
nx = read_delta_r_next(NS / "delta_r_next.csv")
if len(nx["x_rt"]) != len(xN) or np.max(np.abs(np.asarray(nx["x_rt"]) - xN)) > 1e-12:
    raise SystemExit("delta_r_next.csv の列が抽出と合わない")
dC = np.asarray(nx["delta_in_prev"], dtype=float)
rep["saved_dE_vs_A"] = {"max_abs_mm": float(np.max(np.abs(out[4]["dE"] - np.asarray(nx["delta_E"]))) * mm)}
rep["saved_hard_vs_A"] = {"n_diff": int(np.sum(out[4]["hard"] != (np.asarray(saved["hard_ok"]) != 0)))}
hard_same = bool(np.array_equal(out[4]["hard"], out[32]["hard"]))
cal = {}
for R in out:
    Et = float(np.interp(0.0, xN, out[R]["dE"])); Ee = float(out[R]["dE"][-1])
    Ct = float(np.interp(0.0, xN, dC)); Ce = float(dC[-1])
    c = (Ee - Et) / (Ce - Ct); d0 = Et - c * Ct
    corr = c * dC + d0
    th2 = np.degrees(np.arctan(rp + c * np.gradient(dC, xN)))
    cal[R] = {"c": c, "delta0_mm": d0 * mm, "Et_mm": Et * mm, "Ee_mm": Ee * mm, "corr": corr, "theta": th2}
dn = xN >= 0.0
d_corr = (cal[32]["corr"] - cal[4]["corr"]) * mm; d_th2 = cal[32]["theta"] - cal[4]["theta"]
i3 = int(np.argmax(np.abs(d_corr[dn]))); i4 = int(np.argmax(np.abs(d_th2[dn])))
ext = {"calibration": {R: {k: v for k, v in cal[R].items() if k not in ("corr", "theta")} for R in cal},
       "endpoint_abs_dE_mm": {"throat": abs(cal[32]["Et_mm"] - cal[4]["Et_mm"]), "exit": abs(cal[32]["Ee_mm"] - cal[4]["Ee_mm"])},
       "max_abs_dcorr_mm_x_ge_0": float(np.abs(d_corr[dn][i3])), "x_m_at_dcorr": float(xm[dn][i3]),
       "max_abs_dtheta_2coef_deg_x_ge_0": float(np.abs(d_th2[dn][i4])), "x_m_at_dtheta": float(xm[dn][i4]),
       "hard_same_between_arms": hard_same,
       "tolerance": {"endpoint_abs_dE_mm": 0.01, "max_abs_dcorr_mm": 0.01, "max_abs_dtheta_deg": 0.01}}
within = (max(ext["endpoint_abs_dE_mm"].values()) <= 0.01 and ext["max_abs_dcorr_mm_x_ge_0"] <= 0.01
          and ext["max_abs_dtheta_2coef_deg_x_ge_0"] <= 0.01)
reproduced = rep["saved_raw_vs_A"]["n_nan_mismatch"] == 0 and rep["saved_raw_vs_A"]["max_abs_mm"] <= 1e-9 and rep["saved_hard_vs_A"]["n_diff"] == 0
if not reproduced:
    ext["result"] = "判定不能 (A が保存の抽出を再現しない)"
elif not hard_same:
    ext["result"] = "判定不能 (hard の判定が腕で変わった — 単純な比較を合格にしない)"
else:
    ext["result"] = "すべて以内 (標本化の依存を 2 係数の阻害要因から外す)" if within else "超過 (2 点の較正でも抽出の感度は無視できない → 抽出を先に直す)"
rep["two_coef"] = ext
rep["smooth_diag"] = {R: {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v)) for k, v in out[R]["diag"].items()} for R in out}
np.savetxt(O / "refine_ab_3k.csv", np.c_[xN, xm, out[4]["raw"], out[32]["raw"], out[4]["dE"], out[32]["dE"], th[4], th[32],
                                         out[4]["hard"], out[32]["hard"], yb_saved, dC, cal[4]["corr"], cal[32]["corr"], cal[4]["theta"], cal[32]["theta"]],
           delimiter=",", comments="", header="x_rt,x_m,raw_A4,raw_B32,dE_A4,dE_B32,theta_A4_deg,theta_B32_deg,hard_A4,hard_B32,band_y_b_saved,dC,corr2_A4,corr2_B32,theta2_A4_deg,theta2_B32_deg")
(O / "refine_ab_3k.json").write_text(json.dumps(rep, indent=1, ensure_ascii=False, default=str))
print(json.dumps(rep, indent=1, ensure_ascii=False, default=str))
