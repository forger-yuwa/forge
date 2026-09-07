#!/usr/bin/env python3
"""排除厚さ更新計画 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) の単体テスト
(pytest 非依存・素の assert)。

対象:
  1. `core_matched_deficit` (固定 Euler 基準・コア整合質量欠損 → 半径方向等価排除厚)
     - q_NS = c·q_E なら δ_r = 0
     - 既知の壁側欠損から円環厚さを復元 (一様コア / 非一様コア)
     - Euler 壁と NS 壁が異なっても復元 (Euler 壁の外は壁値外挿)
     - 点数と壁近傍伸長率への収束
     - コア振幅差は α で除去、コア形状差はゲート (core_rms) で検出
     - コア範囲 25/30/35 % 感度
  2. 積分法初期推定器 (`feedback/deltastar_integral.py`, CONTUR 運動量積分)
     - CPG 断熱: 全域有限・正、スロートで薄い (平板相関の 1/10 程度)、下流で単調増加
     - CPG 等温 (定数 Tw / Tw(x) テーブル): 走る・有限
     - 指定壁温 → 回復壁温に近づけると断熱解へ連続的に一致
     - 入口 θ0 感度: スロート以降はほぼ無記憶 (×0.5/×2 で 5 % 以内)
     - 半径方向補正 δ_r = δ*_n / cos φ_w
     (設計は case/45 の problem YAML から design_chain で作る: ~1 s)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.metrics.deltastar import core_matched_deficit, band_local_deficit  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def wall_clustered_grid(rw, n=97, first_frac=4.5e-5, stretch=None):
    """壁側に幾何級数で集中した半径格子 (mesh2d._radial_fracs と同趣旨)。"""
    # 壁からの距離 y を等比級数で
    if stretch is None:
        # first_frac·rw が最初の区間になる伸長率を二分法で
        lo, hi = 1.0001, 1.5
        for _ in range(100):
            g = 0.5 * (lo + hi)
            tot = first_frac * (g ** (n - 1) - 1) / (g - 1)
            lo, hi = (g, hi) if tot < 1.0 else (lo, g)
        stretch = g
    dy = first_frac * stretch ** np.arange(n - 1)
    y = np.concatenate([[0.0], np.cumsum(dy)])
    y = y / y[-1]                       # 0 (壁) → 1 (軸)
    return (rw * (1.0 - y))[::-1]       # 軸 → 壁


def profile_with_annulus_deficit(r, q_core_fn, rw_ns, delta_r):
    """壁側の円環 [rw_ns-δ_r, rw_ns] を完全欠損 (q=0)、内側は q_core_fn(r)。
    → 等価排除厚は厳密に δ_r。"""
    q = q_core_fn(r).copy()
    q[r > rw_ns - delta_r] = 0.0
    return q


def profile_smooth_bl(r, q_core_fn, rw_ns, delta99, n_pow=7.0):
    """べき乗則境界層 (q/q_e = (y/δ)^{1/n})。等価排除厚は数値積分で別途求める。"""
    y = rw_ns - r
    f = np.where(y < delta99, (np.maximum(y, 0) / delta99) ** (1.0 / n_pow), 1.0)
    return q_core_fn(r) * f


def equiv_delta_exact(r_fine, q_ref, q_ns, rw):
    """参照解: 十分細かい格子で D と r_eff を直接計算。"""
    D = 2 * np.pi * np.trapezoid((q_ref - q_ns) * r_fine, r_fine)
    seg = 0.5 * (q_ref[1:] * r_fine[1:] + q_ref[:-1] * r_fine[:-1]) * np.diff(r_fine)
    F = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    return rw - float(np.interp(D, F, r_fine[::-1]))


# --- 1. q_NS = c·q_E → δ_r = 0 ------------------------------------------------------
rw = 5.0
r = wall_clustered_grid(rw)
qE = lambda r_: 1.0 + 0.2 * (r_ / rw) ** 2          # 非一様コア (壁側が高い)
res = core_matched_deficit(r, 1.7 * qE(r), qE(r), rw)
check("比例プロファイル: δ_r = 0", abs(res["delta_r"]) < 1e-12, f"{res['delta_r']:.2e}")
check("比例プロファイル: α = 1.7", abs(res["alpha"] - 1.7) < 1e-12)
check("比例プロファイル: コア RMS = 0", res["core_rms_noaxis"] < 1e-12)

# --- 2. 既知の円環欠損 (一様コア) ---------------------------------------------------
qU = lambda r_: np.ones_like(r_)
for d_true in (0.01, 0.05, 0.2):
    qN = profile_with_annulus_deficit(r, qU, rw, d_true)
    res = core_matched_deficit(r, qN, qU(r), rw)
    check(f"一様コア 円環欠損 δ_r={d_true}", abs(res["delta_r"] - d_true) < 2e-3 * rw,
          f"got {res['delta_r']:.5f}")

# --- 3. 非一様コア + 円環欠損 (α 一致・形状一致) --------------------------------------
for d_true in (0.02, 0.1):
    qN = profile_with_annulus_deficit(r, lambda r_: 1.3 * qE(r_), rw, d_true)
    res = core_matched_deficit(r, qN, qE(r), rw)
    check(f"非一様コア (×1.3) 円環欠損 δ_r={d_true}", abs(res["delta_r"] - d_true) < 2e-3 * rw,
          f"got {res['delta_r']:.5f}, α {res['alpha']:.4f}")
    check(f"  α が振幅差 1.3 を吸収", abs(res["alpha"] - 1.3) < 1e-9)
    check(f"  欠損は壁側に集中 (outer_share>0.8)", res["outer_share"] > 0.8, f"{res['outer_share']:.3f}")

# --- 4. Euler 壁 ≠ NS 壁 (NS が δ_in だけ太い; Euler 壁の外は壁値外挿) -----------------
rw_e = 5.0
for d_in in (0.05, 0.3):
    rw_n = rw_e + d_in
    rn = wall_clustered_grid(rw_n)
    qE_map = qE(np.minimum(rn, rw_e))                # Euler 壁の外は壁値
    d_true = 0.12
    qN = profile_with_annulus_deficit(rn, qE, rw_n, d_true)
    res = core_matched_deficit(rn, qN, qE_map, rw_e)
    check(f"壁半径差 δ_in={d_in}: 円環欠損 {d_true} を復元", abs(res["delta_r"] - d_true) < 3e-3 * rw_n,
          f"got {res['delta_r']:.5f}")
    check(f"  delta_in = {d_in}", abs(res["delta_in"] - d_in) < 1e-12)

# --- 5. 滑らかな境界層: 点数・伸長率への収束 ----------------------------------------
rw = 5.0; d99 = 0.4
r_fine = np.linspace(0, rw, 200001)
q_ref_f = qE(r_fine); q_ns_f = profile_smooth_bl(r_fine, qE, rw, d99)
d_exact = equiv_delta_exact(r_fine, q_ref_f, q_ns_f, rw)
errs = []
for n in (65, 97, 193, 385):
    rn = wall_clustered_grid(rw, n=n)
    res = core_matched_deficit(rn, profile_smooth_bl(rn, qE, rw, d99), qE, rw)
    errs.append(abs(res["delta_r"] - d_exact) / d_exact)
check("べき乗則 BL: n=97 で参照解 ±3 %", errs[1] < 0.03, f"exact {d_exact:.5f}, rel errs {['%.2e' % e for e in errs]}")
check("べき乗則 BL: 点数で単調収束", errs[0] > errs[1] > errs[2] > errs[3] * 0.999,
      f"{['%.2e' % e for e in errs]}")
errs_s = []
for ff in (2e-4, 4.5e-5, 1e-5):
    rn = wall_clustered_grid(rw, n=97, first_frac=ff)
    res = core_matched_deficit(rn, profile_smooth_bl(rn, qE, rw, d99), qE, rw)
    errs_s.append(abs(res["delta_r"] - d_exact) / d_exact)
check("べき乗則 BL: 壁近傍伸長率を変えても ±3 %", max(errs_s) < 0.03, f"{['%.2e' % e for e in errs_s]}")

# --- 6. コア形状差はゲートで検出 (振幅差は検出しない) --------------------------------
r = wall_clustered_grid(5.0)
qN_amp = 1.02 * qE(r)                                          # 振幅のみ 2 %
res_amp = core_matched_deficit(r, qN_amp, qE(r), 5.0)
qN_shape = qE(r) * (1.0 + 0.03 * (r / 5.0))                    # 形状 (勾配) 差 3 %
res_shape = core_matched_deficit(r, qN_shape, qE(r), 5.0)
check("振幅差 2 %: コア RMS ≈ 0 (α が吸収)", res_amp["core_rms_noaxis"] < 1e-12)
check("形状差 3 %: コア RMS > 0.1 % で検出", res_shape["core_rms_noaxis"] > 1e-3,
      f"{res_shape['core_rms_noaxis']:.4f}")

# --- 7. コア範囲 25/30/35 % 感度 (円環欠損なら不変) --------------------------------------
qN = profile_with_annulus_deficit(r, lambda r_: 1.1 * qE(r_), 5.0, 0.08)
ds = [core_matched_deficit(r, qN, qE(r), 5.0, core_frac=fc)["delta_r"] for fc in (0.25, 0.30, 0.35)]
check("コア範囲 25/30/35 %: δ_r 不変 (<0.5 %)", max(abs(d - ds[1]) for d in ds) / ds[1] < 5e-3,
      f"{['%.5f' % d for d in ds]}")

# --- 8. 負の欠損 (NS がコアより多く流す) は negative_deficit を返す -------------------
res_neg = core_matched_deficit(r, qE(r) * (1.0 + 0.05 * (r / 5.0) ** 4), qE(r), 5.0)
check("負の欠損: reason=negative_deficit, δ_r<0", "negative_deficit" in res_neg["reason"] and res_neg["delta_r"] < 0,
      f"{res_neg['delta_r']:.4f}")

# =============================== 1b. 帯局所参照 (生産抽出) =============================
# 同じ合成プロファイルで band_local_deficit を検証。加えて「コアに波 (形状差) がある」場合に
# コア全体方式が汚染され、帯局所方式は汚染されないことを確認する (V1 pass1 の知見)。
rw = 5.0; r = wall_clustered_grid(rw)
qE = lambda r_: 1.0 + 0.2 * (r_ / rw) ** 2
resb = band_local_deficit(r, 1.7 * qE(r), qE, rw, delta_in=0.05)
check("帯局所: 比例プロファイルで δ_r = 0", abs(resb["delta_r"]) < 1e-6, f"{resb['delta_r']:.2e}")
for d_true in (0.01, 0.05, 0.2):
    qN = profile_with_annulus_deficit(r, lambda r_: 1.3 * qE(r_), rw, d_true)
    resb = band_local_deficit(r, qN, qE, rw, delta_in=d_true)          # δ_in = 真値 (固定点近傍)
    check(f"帯局所: 円環欠損 δ_r={d_true} を復元 (δ_in=真値)", abs(resb["delta_r"] - d_true) < 2e-3 * rw,
          f"got {resb['delta_r']:.5f}, y_b {resb['band_y_b']:.3f}")
    resb2 = band_local_deficit(r, qN, qE, rw, delta_in=0.3 * d_true)   # δ_in が 3 倍過小でも (帯下限/再試行)
    check(f"帯局所: δ_in 過小 (×0.3) でも復元", abs(resb2["delta_r"] - d_true) < 4e-3 * rw,
          f"got {resb2['delta_r']:.5f}, retry {resb2['band_retry']}")
# 滑らかな BL: 参照解と比較
d99 = 0.4; r_fine = np.linspace(0, rw, 200001)
d_exact = equiv_delta_exact(r_fine, qE(r_fine), profile_smooth_bl(r_fine, qE, rw, d99), rw)
resb = band_local_deficit(r, profile_smooth_bl(r, qE, rw, d99), qE, rw, delta_in=d_exact)
# 適応帯は縁 (帯内の比の変化 <1 %) より外の尾部を取りこぼすため −3 % 程度の系統バイアスを持つ (既知・許容: 出口 M で ≤0.05 %)
check("帯局所: べき乗則 BL を参照解 ±4 % (尾部バイアス −3 % は既知)", abs(resb["delta_r"] - d_exact) / d_exact < 0.04,
      f"exact {d_exact:.5f} got {resb['delta_r']:.5f}")
# コアの波: q_NS = q_E (1 + ε sin(2π r/r_w)) × BL。全域コア整合は汚染、帯局所は復元
eps = 0.01
wave = lambda r_: qE(r_) * (1.0 + eps * np.sin(2 * np.pi * r_ / rw))
qN_w = profile_smooth_bl(r, wave, rw, d99)
res_core = core_matched_deficit(r, qN_w, qE, rw)
res_band = band_local_deficit(r, qN_w, qE, rw, delta_in=d_exact)
err_core = abs(res_core["delta_r"] - d_exact) / d_exact; err_band = abs(res_band["delta_r"] - d_exact) / d_exact
check("コア波 1 %: 帯局所は ±4 % で復元", err_band < 0.04, f"band {err_band:.3f}, core-matched {err_core:.3f}")
check("コア波 1 %: 帯局所の誤差 < コア全体方式の誤差", err_band < err_core, f"band {err_band:.3f} vs core {err_core:.3f}")
# Euler 壁 ≠ NS 壁
rw_n = rw + 0.3; rn = wall_clustered_grid(rw_n)
qN = profile_with_annulus_deficit(rn, qE, rw_n, 0.12)
resb = band_local_deficit(rn, qN, lambda r_: qE(np.minimum(r_, rw)), rw)
check("帯局所: 壁半径差 0.3 でも円環欠損 0.12 を復元", abs(resb["delta_r"] - 0.12) < 3e-3 * rw_n, f"got {resb['delta_r']:.5f}")

# =============================== 2. 積分法初期推定器 ==================================
from forge_design.probdef import load_problem  # noqa: E402
from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.feedback.deltastar_integral import integral_bl, delta_r_function  # noqa: E402
from forge_design.feedback.deltastar import dstar_flatplate  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
prob = load_problem(ROOT / "case/45.isobutane_m6_d155/problem_d155_ns.yaml")
dsn = design_chain(prob)
S = float(prob.spec["r_throat"]); Pt = float(prob.spec["Pt"]); Tt = float(prob.spec["Tt"])
ad = integral_bl(dsn["wall"], dsn["wall_inv"], prob.gamma, prob.cp, Pt, Tt, S)
k0 = int(np.argmin(np.abs(ad["x"])))
check("積分法 CPG 断熱: 全域有限", bool(np.all(np.isfinite(ad["delta_r"]))))
check("積分法 CPG 断熱: 全域正", bool(np.all(ad["delta_r"] > 0)), f"min {ad['delta_r'].min():.2e}")
m_dn = ad["x"] > 0.5
check("積分法 CPG 断熱: スロート下流で単調増加", bool(np.all(np.diff(ad["delta_r"][m_dn]) > -1e-9)))
# 平板相関 (入口起点 s) との比較: スロートで一桁薄い (加速による薄化を表現)
s_arc = (0.0 - dsn["wall"].x_in) * S
fp_t = float(dstar_flatplate(np.array([s_arc]), np.array([1.0]), Pt, Tt, prob.gamma, prob.cp)[0]) / S
check("積分法 CPG 断熱: スロート δ*_n が平板相関の 1/4 以下", ad["dstar_n"][k0] < 0.25 * fp_t,
      f"integral {ad['dstar_n'][k0]:.5f} vs flat-plate {fp_t:.5f}")
check("積分法: δ_r = δ*_n / cos φ_w", np.allclose(ad["delta_r"], ad["dstar_n"] / ad["cos_phi"]))
check("積分法: 断熱では Tw = Taw", np.allclose(ad["Tw"], ad["Taw"]))
iso = integral_bl(dsn["wall"], dsn["wall_inv"], prob.gamma, prob.cp, Pt, Tt, S,
                  thermal_bc={"mode": "prescribed_temperature", "Tw": 600.0})
check("積分法 CPG 等温 600 K: 全域有限・正", bool(np.all(np.isfinite(iso["delta_r"])) and np.all(iso["delta_r"] > 0)))
check("積分法 CPG 等温: Tw = 600 K が使われる", np.allclose(iso["Tw"], 600.0))
tab = integral_bl(dsn["wall"], dsn["wall_inv"], prob.gamma, prob.cp, Pt, Tt, S,
                  thermal_bc={"mode": "prescribed_temperature", "Tw_table": np.c_[ad["x"], ad["Taw"]].tolist()})
check("積分法 等温 (Tw=Taw テーブル) → 断熱解と一致 (<0.5 %)",
      float(np.max(np.abs(tab["delta_r"] - ad["delta_r"]) / np.maximum(ad["delta_r"], 1e-6))) < 5e-3)
# 連続性: Tw を Taw の 0.9, 0.99, 1.0 倍にすると断熱へ近づく
devs = []
for f in (0.9, 0.99):
    r_ = integral_bl(dsn["wall"], dsn["wall_inv"], prob.gamma, prob.cp, Pt, Tt, S,
                     thermal_bc={"mode": "prescribed_temperature", "Tw_table": np.c_[ad["x"], f * ad["Taw"]].tolist()})
    devs.append(float(np.max(np.abs(r_["delta_r"] - ad["delta_r"]) / np.maximum(ad["delta_r"], 1e-6))))
check("積分法 等温→断熱の連続性 (偏差が単調減少)", devs[0] > devs[1], f"{devs}")
# θ0 感度: スロート以降はほぼ無記憶
th = [integral_bl(dsn["wall"], dsn["wall_inv"], prob.gamma, prob.cp, Pt, Tt, S, theta0_m=f * ad["settings"]["theta0_m"])
      for f in (0.5, 2.0)]
dev_t = max(abs(t["dstar_n"][k0] - ad["dstar_n"][k0]) / ad["dstar_n"][k0] for t in th)
check("積分法 θ0 ×0.5/×2: スロート δ*_n の変化 < 5 %", dev_t < 0.05, f"{dev_t:.3e}")
fn = delta_r_function(ad)
check("積分法 delta_r_function: 補間が一致", abs(float(fn(ad["x"][k0])) - ad["delta_r"][k0]) < 1e-12)

# =============================== 3. δ_r の 5 次 P-spline 平滑化 ==========================
from forge_design.metrics.deltastar import smooth_delta_quintic  # noqa: E402
xg = np.linspace(-1.0, 95.0, 1250)
d_true = 0.0015 + 0.0075 * (np.maximum(xg, 0) / 95.0) ** 0.5 * 95.0 / 8.0 * 0.6   # 滑らかな成長 (スロート 0.0015 → 出口 ~0.7)
rng = np.random.default_rng(0)
d_noisy = d_true * (1.0 + 0.03 * rng.standard_normal(len(xg)))              # 抽出ノイズ 3 %
wts = np.ones_like(xg); wts[(xg < -0.5)] = 0.0                              # hard 不合格を混ぜる
f_s, dg = smooth_delta_quintic(xg, d_noisy, weights=wts, knot_spacing=2.0, lam=1.0)
d_s = f_s(xg)
err = np.abs(d_s - d_true) / d_true
check("P-spline: 3 % ノイズから真値を ±1.5 % (中央値) で復元", float(np.median(err[xg > 0])) < 0.015, f"median {np.median(err[xg>0]):.4f}")
check("P-spline: 全域正・有限", bool(np.all(np.isfinite(d_s)) and np.all(d_s >= 0)))
# 曲率の滑らかさ: 2 階差分の高周波成分がノイズ入り生値より 2 桁小さい
def hf2(v):
    d2 = np.gradient(np.gradient(v, xg), xg); w = 25; return float(np.std(d2 - np.convolve(d2, np.ones(w) / w, mode="same")))
check("P-spline: δ'' の高周波が生値の 1/50 以下", hf2(d_s) < hf2(d_noisy) / 50, f"{hf2(d_s):.2e} vs {hf2(d_noisy):.2e}")
check("P-spline: 重み 0 の点を無視 (x<-0.5 の外挿値が有限)", bool(np.all(np.isfinite(f_s(np.array([-0.9, -0.7]))))))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
