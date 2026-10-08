"""冷却壁 (300 K) の NS の対 (plan tooling-nozzle-isothermal-wall-chain §5.1 #14、§6 V-c45) のメッシュ設計と問題 YAML の生成。
2026-10-08 ユーザ「冷却壁は y+ 上がっちゃうから気をつけてね」。

入力 (どちらも _band_ab/cold_pair/ に置く):
  y1p_ad_wall.csv  : 生産の断熱 NS (run_0179_ns_n012_N2_ext、step 20000) の壁の各点の y1・y1+・ρ_w・接線せん断・μ_w・T_w
                     (`check_wall_resolution.py run_0179_ns_n012_N2_ext --groups wall --profile-csv y1p_ad.csv`、AWS で作成)
  ../delta_contur/hform_ab_profiles.json : CONTUR の C_f (断熱 / 300 K、エンタルピー形) — 冷却での摩擦の増え方に使う
冷却の y1+ の見積もり: y1+ ∝ y1·√(ρ_w τ_w)/μ_w。同じ圧力で ρ_w ∝ 1/T_w、τ_w は CONTUR の C_f 比、μ_w は混合気の μ(T)
  (NS と同じ CEA 輸送物性の値; MU_T・MU)。局所の倍率は 8.5〜9.4 (case/44 の空気 1060 K では 5〜5.4)。
設計: 第一セル厚 / 局所半径 = 0.8 / (局所の冷却 y1+ / 第一セル比) を表にする (入口の角 0.05 r_t は幾何的特異点として除く)。
  x 方向の密度は h ≤ min(0.06 r_t, 4500 · 第一セル厚) (AR ≤ 4500 の見込み、判定は check_mesh_quality --ar-max 5000)。
  ni はその密度の積分の 1.05 倍、nj 121 (半径方向の等比 ≤ 1.11、今の生産と同程度)。
出力: _band_ab/cold_pair/mesh_cold_tables.json、problem_d155_ns_prod_coldmesh.yaml (断熱)、problem_d155_ns_prod_coldmesh_tw300.yaml (300 K)。
usage: /home/sano/work/forge/design/.venv-opt/bin/python cold_pair_mesh.py
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design"))
from forge_design.meshing.mesh2d import Mesh2DParams, _first_frac_profile, _radial_fracs, _x_stations  # noqa: E402

OUT = HERE / "_band_ab" / "cold_pair"
SRC = "problem_d155_ns_prod.yaml"
RT = 0.07666536551630307
TARGET, AR_T, H_FAR, NJ, TW = 0.8, 4500.0, 0.06, 121, 300.0
MU_T = [250, 300, 350, 400, 600, 800, 1000, 1200, 1400, 1470, 1600]
MU = [15.281, 17.804, 20.166, 22.398, 30.376, 37.332, 43.613, 49.387, 54.804, 56.638, 59.975]   # μ_mix [µPa·s] (transport_reference.py)


def design():
    A = np.loadtxt(OUT / "y1p_ad_wall.csv", delimiter=",", skiprows=1)
    A = A[np.argsort(A[:, 0])]
    x = A[:, 0] / RT; rw = A[:, 1] / RT; y1 = A[:, 3]; yp = A[:, 4]; Tw = A[:, 8]
    HP = json.loads((HERE / "_band_ab/delta_contur/hform_ab_profiles.json").read_text())
    xh = np.array(HP["x"])
    cf = np.array(HP["B_Tw300"]["Cf"]) / np.array(HP["B_adiabatic"]["Cf"])
    cfr = np.interp(x, xh, cf)                       # 範囲外 (x < xh[0]) は端値 (入口の直管: 下で確認)
    fac = np.sqrt((Tw / TW) * cfr) * np.interp(Tw, MU_T, MU) / np.interp(TW, MU_T, MU)
    frac_now = y1 / (rw * RT)
    k = yp * fac / frac_now                          # 冷却の y1+ / 第一セル比
    x0, x1 = float(x[0]), float(x[-1])
    corner = x < x0 + 0.05
    k[corner] = k[~corner][0]
    xt = np.unique(np.r_[np.arange(x0, -4.0, 0.5), np.arange(-4.0, 6.0, 0.25), np.arange(6.0, x1, 1.0), x1])
    xt[0], xt[-1] = x0, x1
    req = TARGET / k
    ft = [float(req[(x >= (xt[i - 1] if i else x0)) & (x <= (xt[i + 1] if i + 1 < len(xt) else x1))].min()) for i in range(len(xt))]
    fr_tab = [[round(float(a), 6), float(f"{b:.4e}")] for a, b in zip(xt, ft)]
    fr_tab[0][0], fr_tab[-1][0] = x0 - 1e-6, x1 + 1e-6              # 丸めで範囲を外さない
    xd = np.linspace(x0, x1, 4001)
    tf = np.array(fr_tab)
    h = np.minimum(H_FAR, AR_T * np.exp(np.interp(xd, tf[:, 0], np.log(tf[:, 1]))) * np.interp(xd, x, rw))
    dens = 1.0 / h
    ni = int(np.ceil(1.05 * np.trapezoid(dens, xd))) + 1
    xdt = np.unique(np.r_[np.arange(x0, x1, 0.25), x1])
    dt = [float((dens[(xd >= (xdt[i - 1] if i else x0)) & (xd <= (xdt[i + 1] if i + 1 < len(xdt) else x1))].max()) / dens.min())
          for i in range(len(xdt))]
    den_tab = [[round(float(a), 6), float(f"{b:.4f}")] for a, b in zip(xdt, dt)]
    den_tab[0][0], den_tab[-1][0] = x0 - 1e-6, x1 + 1e-6
    # 見積もりの検査 (格子の節点で)
    prm = Mesh2DParams(ni=ni, nj=NJ, wall_first_frac_table=fr_tab, x_density_table=den_tab)
    xs = _x_stations(x0, x1, ni, 1.0, 1.0, density_table=den_tab)
    fr = _first_frac_profile(xs, prm)
    rws = np.interp(xs, x, rw)
    ar = np.gradient(xs) / (fr * rws)
    ypc = np.interp(xs, x, k) * fr
    area = np.gradient(xs) * rws
    g = np.diff(_radial_fracs(NJ, float(fr.min())))[::-1]
    est = {"ni": ni, "nj": NJ, "nodes": ni * NJ, "ar_est_max": float(ar.max()), "x_ar_est_max": float(xs[np.argmax(ar)]),
           "y1p_cold_est_max": float(ypc.max()), "y1p_cold_est_over1_area_pct": float(100 * area[ypc > 1].sum() / area.sum()),
           "cooling_factor_range": [float(fac.min()), float(fac.max())], "radial_ratio_max": float(g[1] / g[0]),
           "frac_range": [float(fr.min()), float(fr.max())], "target_y1p": TARGET, "ar_target": AR_T, "h_far": H_FAR,
           "inputs": {"y1p_ad_wall.csv": "run_0179_ns_n012_N2_ext step 20000 (check_wall_resolution --profile-csv)",
                      "cf_ratio": "CONTUR enthalpy form (hform_ab_profiles.json B_Tw300/B_adiabatic)"}}
    return fr_tab, den_tab, est


def flow(tbl) -> str:
    return "[" + ", ".join(f"[{a}, {b}]" for a, b in tbl) + "]"


def write_problems(fr_tab, den_tab, est):
    src = (HERE / SRC).read_text()
    m = re.search(r"^mesh:\n(?:  .*\n|\s*\n)+?(?=^\S)", src, flags=re.M)
    if not m:
        raise SystemExit("mesh ブロックが見つからない")
    block = ("mesh:\n  discretization: node\n"
             f"  ni: {est['ni']}\n  nj: {NJ}\n"
             "  # 冷却壁 (300 K) 用 (cold_pair_mesh.py、plan tooling-nozzle-isothermal-wall-chain §5.1 #13・#14): 第一セル厚の表 (局所半径比、log 線形) と\n"
             "  # x 方向の相対密度の表。断熱 NS の y1+ × 局所の冷却倍率 8.5〜9.4 から y1+_cold ≤ 0.8、AR ≤ 4500 の見込みで決めた。\n"
             "  # 断熱 (対の片方) も同じ格子で回す。壁法線の構造層なので AR ≤ 5000 (check_mesh_quality --ar-max 5000)\n"
             f"  wall_first_frac: {fr_tab[-1][1]}\n"
             f"  wall_first_frac_table: {flow(fr_tab)}\n"
             f"  x_density_table: {flow(den_tab)}\n"
             "  throat_refine: 1.0\n  throat_width: 1.0\n  ar_max: 5000\n\n")
    base = src[:m.start()] + block + src[m.end():]
    base = base.replace("name: problem_d155_ns_prod\n", "name: problem_d155_ns_prod_coldmesh\n", 1)
    head = ("# 冷却壁の NS の対 (plan tooling-nozzle-isothermal-wall-chain §5.1 #14、§6 V-c45) の断熱側。cold_pair_mesh.py が生産の\n"
            f"# {SRC} から作る: 違いは name と mesh ブロック (冷却壁用の格子) だけ。壁は生産と同じ (delta_r_csv = run_0167 の delta_r_initial.csv)。\n")
    ad = head + base
    tw = ad.replace("name: problem_d155_ns_prod_coldmesh\n", "name: problem_d155_ns_prod_coldmesh_tw300\n", 1)
    tw = tw.replace("の断熱側。", "の 300 K 側 (spec.wall_thermal だけが断熱側と違う)。", 1)
    tw = tw.replace("spec:\n", f"spec:\n  wall_thermal: {{mode: isothermal, Tw: {TW}}}\n", 1)
    for name, txt in (("problem_d155_ns_prod_coldmesh.yaml", ad), ("problem_d155_ns_prod_coldmesh_tw300.yaml", tw)):
        (HERE / name).write_text(txt)
    return ad, tw


if __name__ == "__main__":
    fr_tab, den_tab, est = design()
    (OUT / "mesh_cold_tables.json").write_text(json.dumps({"wall_first_frac_table": fr_tab, "x_density_table": den_tab, "estimate": est},
                                                          indent=1, ensure_ascii=False))
    write_problems(fr_tab, den_tab, est)
    print(json.dumps(est, indent=1, ensure_ascii=False))
