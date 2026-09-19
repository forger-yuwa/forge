#!/usr/bin/env python3
"""case/49 の作動条件を conditions.json から導出する単一ソース (plan §3.1)。

usage:
  python3 setup.py                 # 導出結果を表示して derived.json に保存
  python3 setup.py --json          # JSON だけ出す

他スクリプト (gen_runs.py / precursor / tools) は `from setup import load` で使う。

物性は **forge の実装値に合わせる**:
  - Sutherland (`viscMethod: 1`): mu0 1.716e-5, T0 273.0 K, S 111.0 K
    (solver_density_cuda/cuda_forge/gasProperties_d.cu)。`physProp.visc` は viscMethod 0 専用なので書かない
  - R = cp (gamma-1)/gamma   (case/48 と同じ。cp 1004.5, gamma 1.4 -> 287.0)
標準大気は P, T を与えるためだけに使い (R_ISA=287.05)、rho は forge の R で作る。
"""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

MU0, T0_SUTH, S_SUTH = 1.716e-5, 273.0, 111.0   # forge 実装値
_MW = {"N2": 0.0280134, "O2": 0.0319988, "AR": 0.039948, "CO2": 0.0440095,
       "H2O": 0.01801528, "CO": 0.0280104, "NO": 0.0300061}   # kg/mol
R_ISA, G0 = 287.05, 9.80665


def mu_sutherland(T):
    return MU0 * (T / T0_SUTH) ** 1.5 * (T0_SUTH + S_SUTH) / (T + S_SUTH)


def isa(z_km):
    """標準大気 (0-32 km) の (T [K], P [Pa])"""
    z = z_km * 1000.0
    if z <= 11000.0:
        T = 288.15 - 6.5e-3 * z
        P = 101325.0 * (T / 288.15) ** (G0 / (R_ISA * 6.5e-3))
    elif z <= 20000.0:
        T = 216.65
        P11 = 101325.0 * (216.65 / 288.15) ** (G0 / (R_ISA * 6.5e-3))
        P = P11 * math.exp(-G0 * (z - 11000.0) / (R_ISA * 216.65))
    elif z <= 32000.0:
        T = 216.65 + 1.0e-3 * (z - 20000.0)
        P11 = 101325.0 * (216.65 / 288.15) ** (G0 / (R_ISA * 6.5e-3))
        P20 = P11 * math.exp(-G0 * 9000.0 / (R_ISA * 216.65))
        P = P20 * (T / 216.65) ** (-G0 / (R_ISA * 1.0e-3))
    else:
        raise ValueError("altitude_km > 32 は未対応")
    return T, P


def tp_gas(Y):
    """semi-perfect 乾燥空気 (NASA-9)。design/forge_design が無ければ None"""
    import sys
    sys.path.insert(0, str(ROOT / "design"))
    try:
        from forge_design.gas.semiperfect import GasSemiPerfect
    except Exception as e:                       # noqa: BLE001
        print("[setup] WARNING: semi-perfect gas 不可 (%s)" % e)
        return None
    return GasSemiPerfect(Y, Tt=2000.0)


def _solve_T_of_h(gas, h_target, lo=150.0, hi=3000.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(gas.h_mass(mid)) < h_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def load_case(path=None):
    """最上位入力 case.json (旧 conditions.json 互換)。"""
    p = Path(path or HERE / "case.json")
    if not p.exists():
        p = HERE / "conditions.json"
    c = json.loads(p.read_text())
    return c if "conditions" in c else {"conditions": c}


def load(cfg_path=None):
    C = load_case(cfg_path)["conditions"]
    M, cp, gam = C["mach"], C["cp"], C["gamma"]
    R = cp * (gam - 1.0) / gam

    if C.get("P_inf") and C.get("T_inf"):
        T_inf, P_inf, src = C["T_inf"], C["P_inf"], "explicit"
    else:
        T_inf, P_inf = isa(C["altitude_km"])
        src = "ISA %.1f km" % C["altitude_km"]

    ro = P_inf / (R * T_inf)
    a = math.sqrt(gam * R * T_inf)
    U = M * a
    mu = mu_sutherland(T_inf)
    Re_m = ro * U / mu
    r_rec = C["prandtl_lam"] ** (1.0 / 3.0)

    d = dict(
        source=src, mach=M, R=R, cp=cp, gamma=gam,
        T_inf=T_inf, P_inf=P_inf, ro_inf=ro, a_inf=a, U_inf=U, mu_inf=mu, Re_per_m=Re_m,
        recovery_factor=r_rec,
        Tt_cpg=T_inf * (1.0 + 0.5 * (gam - 1.0) * M * M),
        Taw_cpg=T_inf * (1.0 + r_rec * 0.5 * (gam - 1.0) * M * M),
        Pt_cpg=P_inf * (1.0 + 0.5 * (gam - 1.0) * M * M) ** (gam / (gam - 1.0)),
        wall_T=C["wall_T"], cyl_top_thermal=C["cyl_top_thermal"],
        prandtl_lam=C["prandtl_lam"], prandtl_turb=C["prandtl_turb"],
        delta_target_m=C["delta_target_mm"] * 1e-3, gas=C["gas"], dry_air_Y=C["dry_air_Y"],
    )
    # 自由流乱流 (k, omega): k=1.5(TI U)^2, omega = ro k / (mu (mut/mu))
    TI = C["turb_intensity"]
    d["k_inf"] = 1.5 * (TI * U) ** 2
    d["omega_inf"] = ro * d["k_inf"] / (mu * C["mut_over_mu"])
    d["turb_intensity"] = TI
    d["mut_over_mu"] = C["mut_over_mu"]

    # semi-perfect (TP) の全温・回復温度・cp 比 (CPG の偏りを測るため常に出す)
    g = tp_gas(C["dry_air_Y"])
    if g is not None:
        h_inf = float(g.h_mass(T_inf))
        d["Tt_tp"] = _solve_T_of_h(g, h_inf + 0.5 * U * U)
        d["Taw_tp"] = _solve_T_of_h(g, h_inf + r_rec * 0.5 * U * U)
        d["cp_tp_298"] = float(g.cp_mass(298.15))
        d["cp_tp_at_Tt"] = float(g.cp_mass(d["Tt_tp"]))
        d["gamma_tp_at_Tt"] = float(g.gamma(d["Tt_tp"]))
        d["cp_rise_pct"] = 100.0 * (d["cp_tp_at_Tt"] / d["cp_tp_298"] - 1.0)
        d["Taw_cpg_minus_tp"] = d["Taw_cpg"] - d["Taw_tp"]
    d["Tw_over_Taw"] = d["wall_T"] / d.get("Taw_tp", d["Taw_cpg"])
    # TP (semi-perfect 1 擬似種) の気体定数と、それに整合する自由流密度。
    # gas == "TP" のときは BC/IC をこちらで組む (CPG の R 287.0 と 0.017 % 違う)。
    MW_mix = 1.0 / sum(y / _MW[k] for k, y in C["dry_air_Y"].items())
    d["R_tp"] = 8.314462618 / MW_mix
    d["MW_tp"] = MW_mix
    if C["gas"].upper() == "TP":
        d["ro_inf"] = P_inf / (d["R_tp"] * T_inf)
        d["Re_per_m"] = d["ro_inf"] * U / mu
        d["k_inf"] = 1.5 * (TI * U) ** 2
        d["omega_inf"] = d["ro_inf"] * d["k_inf"] / (mu * C["mut_over_mu"])
    return d


def geometry_manifest(G):
    """幾何 (mm 入力) を SI に直し、すきま・グループ別解析面積を導出する。
    **面積式の正本はここ 1 箇所** (build_geom.py / mesh_salome.py / 評価はこれを読む)。"""
    m = {k: (v * 1e-3 if k in ("Ro", "Ri", "x_off", "depth", "x_in", "x_plate", "x_out",
                               "y_max", "z_top", "r_patch", "protrude") else v)
         for k, v in G.items() if not k.startswith("_")}
    Ro, Ri, off, d = m["Ro"], m["Ri"], m["x_off"], m["depth"]
    if abs(off) >= Ro - Ri:
        raise SystemExit("|x_off| が すきま Ro-Ri 以上 (内側円柱が外筒に接触)")
    Lx = m["x_out"] - m["x_in"]
    half = bool(m.get("half_model", True))
    fy = 1.0 if half else 2.0          # y 方向の広がり倍率 (全周は -y_max..+y_max)
    has_runup = m["x_plate"] > m["x_in"] + 1e-12
    has_patch = m["r_patch"] > Ro
    plug = bool(m.get("plug_cavity", False))
    ring = math.pi * (Ro ** 2 - Ri ** 2) / 2.0 * fy            # 環 (偏心でも不変)
    a = {"inlet": fy * m["y_max"] * m["z_top"], "outlet": fy * m["y_max"] * m["z_top"],
         "top": Lx * m["y_max"] * fy, "side": Lx * m["z_top"] * fy}
    if half:
        a["sym"] = Lx * m["z_top"] + (0.0 if plug else 2.0 * (Ro - Ri) * d)
    if not plug:
        a.update(cav_outer=math.pi * Ro * d * fy, cyl_side=math.pi * Ri * d * fy,
                 cav_floor=ring, cyl_top=math.pi * Ri ** 2 / 2.0 * fy)
    rp = m["r_patch"] if has_patch else Ro
    a_plate_all = ((m["x_out"] - m["x_plate"]) if has_runup else Lx) * m["y_max"] * fy
    if has_runup:
        a["runup"] = (m["x_plate"] - m["x_in"]) * m["y_max"] * fy
    if has_patch:
        a["plate_in"] = math.pi * (rp ** 2 - (Ro ** 2 if not plug else 0.0)) / 2.0 * fy
        a["plate"] = a_plate_all - math.pi * rp ** 2 / 2.0 * fy
    else:
        a["plate"] = a_plate_all - (ring if not plug else 0.0) \
            - (math.pi * Ri ** 2 / 2.0 * fy if not plug else 0.0)
    m.update(gap_nom=Ro - Ri, gap_min=Ro - Ri - abs(off), gap_max=Ro - Ri + abs(off),
             has_runup=has_runup, has_patch=has_patch, plug_cavity=plug, half_model=half,
             group_area_m2=a, total_area_m2=sum(a.values()),
             group_area_mm2={k: v * 1e6 for k, v in a.items()})
    return m


def mesh_manifest(M, geo):
    """メッシュ設定を stage 解決し、VL (第一層/stretch/層数/総厚) を gap_min から引く。"""
    st = M["stage"]
    pick = lambda k: (M[k][st] if isinstance(M[k], dict) else M[k])          # noqa: E731
    f = pick("vl_first_um") * 1e-6
    s_ = pick("vl_stretch")
    t_max = M["vl_total_frac_gap"] * geo["gap_min"]
    n = 1
    while True:
        tot = f * (s_ ** (n + 1) - 1.0) / (s_ - 1.0)
        if tot > t_max:
            break
        n += 1
    total = f * (s_ ** n - 1.0) / (s_ - 1.0)
    out = {"stage": st, "vl_first_m": f, "vl_stretch": s_, "vl_nlayers": n,
           "vl_total_m": total, "vl_last_m": f * s_ ** (n - 1),
           "vl_total_max_m": t_max, "minh_m": M["minh"] * 1e-3, "fineness": M["fineness"],
           "node_budget": pick("node_budget")}
    for k in ("maxh", "size_plate", "size_plate_in", "size_cav", "size_floor",
              "size_cyl_top", "size_edge_opening"):
        out[k + "_m"] = pick(k) * 1e-3
    # ---- 整合ガード (2026-09-19 実測) ----
    # VL 総厚が**壁面の接線セルサイズを超える**と、開口リップのような凸角で層が自己交差し、
    # NETGEN の tet 充填が**無言で失敗**する (Compute: True なのに tet 0 / 体積の 2 % しか
    # 埋まらない)。実測: VL 0.8695 mm × リップ 0.3 mm -> tet 0、VL 0.426 mm × リップ 0.5 mm
    # と VL 0.8695 mm × リップ 1.0 mm はいずれも成功。よって全壁面サイズに下限を課す。
    lim = out["vl_total_m"]
    out["vl_size_floor_m"] = lim
    # 下限は **壁面サイズ全部** に課す (2026-09-19 実測 3 例):
    #   リップ 0.30 / 面 0.8-1.2 × VL 0.87 -> tet 0 (失敗)
    #   リップ 1.00 / 面 0.8-1.2 × VL 0.87 -> 成功
    #   リップ 0.90 / 面 0.5     × VL 0.90 -> tet 0 (失敗)   ← エッジだけの下限では不十分
    # すなわち「凸角まわりの**接線セルサイズ**が VL 総厚を下回ると層が自己交差する」。
    # よって**キャビティ内面を細かくしたいときは VL 総厚も一緒に薄くする**必要がある
    # (第一層 y1 は独立に決められるので y+ は保てる)。
    bumped = {}
    if M.get("guard_edge_by_vl", True):
        for k in ("size_plate_in", "size_cav", "size_floor", "size_cyl_top", "size_edge_opening"):
            if out[k + "_m"] < lim:
                bumped[k] = (out[k + "_m"], lim)
                out[k + "_m"] = lim
    out["size_bumped_by_vl"] = {k: [a * 1e3, b * 1e3] for k, (a, b) in bumped.items()}
    if "hex" in M:                       # 全ヘキサ版 (cad/build_hex_mesh.py) の分割数
        out["hex"] = {k: v for k, v in M["hex"].items() if not k.startswith("_")}
    # 継ぎ目の段差: 最終プリズム層厚 / 隣接する tet の代表サイズ。1 に近いほど滑らか。
    wall_min = min(out["size_cav_m"], out["size_plate_in_m"], out["size_floor_m"], out["size_cyl_top_m"])
    out["junction_ratio"] = out["vl_last_m"] / wall_min
    out["growth_rate"] = M.get("growth_rate", 0.3)
    return out


def resolve(path=None):
    C = load_case(path)
    man = {"conditions": load(path), "geometry": geometry_manifest(C["geometry"]),
           "eval": {k: v for k, v in C["eval"].items() if not k.startswith("_")}}
    man["mesh"] = mesh_manifest(C["mesh"], man["geometry"])
    man["eval"]["shrink_m"] = man["eval"].pop("shrink_mm") * 1e-3
    # physID の割り当て (正本)。面グループ -> physID、流体ボリューム -> fluid_id。
    order = ["inlet", "outlet", "top", "side", "sym", "runup",
             "plate", "plate_in", "cav_outer", "cav_floor", "cyl_side", "cyl_top"]
    present = [g for g in order if g in man["geometry"]["group_area_m2"]]
    man["phys_id"] = {g: i + 1 for i, g in enumerate(present)}
    man["fluid_id"] = 100
    man["_provenance"] = {"source": "case.json", "tool": "setup.py --resolve"}
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--resolve", action="store_true", help="manifest.json を書く")
    ap.add_argument("--config", default=None)
    a = ap.parse_args()
    if a.resolve:
        man = resolve(a.config)
        (HERE / "manifest.json").write_text(json.dumps(man, indent=2, ensure_ascii=False))
        g, me = man["geometry"], man["mesh"]
        print("gap nom/min/max [mm]: %.3f / %.3f / %.3f" % (g["gap_nom"] * 1e3, g["gap_min"] * 1e3, g["gap_max"] * 1e3))
        print("VL stage %s: first %.1f um x %.3f x %d layers -> total %.4f mm (<= %.4f), last %.4f mm"
              % (me["stage"], me["vl_first_m"] * 1e6, me["vl_stretch"], me["vl_nlayers"],
                 me["vl_total_m"] * 1e3, me["vl_total_max_m"] * 1e3, me["vl_last_m"] * 1e3))
        print("groups:", ", ".join("%s %.1f mm2" % (k, v) for k, v in sorted(g["group_area_mm2"].items())))
        print("wrote", HERE / "manifest.json")
        return
    d = load(a.config)
    (HERE / "derived.json").write_text(json.dumps(d, indent=2, ensure_ascii=False))
    if a.json:
        print(json.dumps(d, indent=2, ensure_ascii=False))
        return
    print("=== case/49 作動条件 (%s, M=%.3g) ===" % (d["source"], d["mach"]))
    rows = [
        ("T_inf [K]", d["T_inf"]), ("P_inf [Pa]", d["P_inf"]), ("ro_inf [kg/m3]", d["ro_inf"]),
        ("a_inf [m/s]", d["a_inf"]), ("U_inf [m/s]", d["U_inf"]), ("mu_inf [Pa s]", d["mu_inf"]),
        ("Re/m [1/m]", d["Re_per_m"]), ("Pt (CPG) [Pa]", d["Pt_cpg"]),
        ("Tt  (CPG) [K]", d["Tt_cpg"]), ("Taw (CPG) [K]", d["Taw_cpg"]),
    ]
    if "Tt_tp" in d:
        rows += [("Tt  (TP)  [K]", d["Tt_tp"]), ("Taw (TP)  [K]", d["Taw_tp"]),
                 ("cp(Tt)/cp(298)-1 [%]", d["cp_rise_pct"]),
                 ("Taw CPG-TP [K]", d["Taw_cpg_minus_tp"])]
    rows += [("wall_T [K]", d["wall_T"]), ("Tw/Taw", d["Tw_over_Taw"]),
             ("k_inf [m2/s2]", d["k_inf"]), ("omega_inf [1/s]", d["omega_inf"]),
             ("delta target [m]", d["delta_target_m"])]
    for k, v in rows:
        print("  %-22s %12.6g" % (k, v))
    print("wrote", HERE / "derived.json")


if __name__ == "__main__":
    main()
