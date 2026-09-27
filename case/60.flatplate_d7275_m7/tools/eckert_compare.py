#!/usr/bin/env python3
"""D-7275 試験 26: forge・実験・Eckert 参照温度法を同じ座標 (St*_l vs R*_l) で並べる (case/60)。

    python3 tools/eckert_compare.py PLATE_H5 [--json out.json]

Eckert の式は D-7275 p.14 (ref. 12 = Kays) の乱流の関係:
    St*_l = 0.0296 · Pr*^(−2/3) · (R*_l)^(−1/5)
座標は論文の Table III の換算係数 (R*_l/R_l,∞ = 0.180、(ρVcp)*_l/(ρVcp)_∞ = 0.413、T* = 600 K) をそのまま使う
(実験点の還元と同じ)。Pr* は forge と同じ燃焼ガスモデル (case/56 gas_htst、R = 296.2) の T* での値。

forge の St*_l は q_w / [(ρVcp)*_l (Taw − Tw)]、Taw = 1728 K (Table III)、Tw = 300 K (run)。
感度として (a) T* を Tw = 300 K から自前で組み直した場合、(b) 乱流起点をトリップ位置 (12.7 cm) にした場合も出す。
(c) --forge-props (forge の run の場から取った T* での μ・cp・Pr、`{"600": {"mu":..,"cp":..,"Pr":..}}`) を渡すと、
Eckert の q を forge 自身の輸送物性 (kinetic theory、Pr ≈ 0.745) で組み、forge の q_w と直接比べる (2026-09-27 追加:
Python 側の気体モデルの Pr* 0.6905 は forge の Pr 0.745 と違い、Eckert の値が約 5 % 変わる)。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
T56 = CASE.parent / "56.gap_tp1187" / "tools"
sys.path.insert(0, str(T56))
from gas_htst import htst                                        # noqa: E402
import importlib.util as _u                                      # noqa: E402
_sp = _u.spec_from_file_location("t4c", T56 / "conditions.py")
_m = _u.module_from_spec(_sp); _sp.loader.exec_module(_m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plate_h5")
    ap.add_argument("--json", default=None)
    ap.add_argument("--tw", type=float, default=300.0)
    ap.add_argument("--forge-props", default=None, help="forge の場から取った T* での μ・cp・Pr (JSON)")
    ap.add_argument("--case-setup", default=None, help="run の case_setup.json (p∞・U∞・R)。--forge-props と併用")
    a = ap.parse_args()
    d = json.loads((CASE / "digitize_d7275_fig20_test26.json").read_text())
    t2, t3 = d["test26_table"]["TableII"], d["test26_table"]["TableIII"]
    gas = _m.Gas(htst(float(t2["Tt_K"])))
    Rfac = t3["RstarLl_over_RLinf"] * t2["Re_inf_per_m"]           # R*_l = Rfac · x
    rvc = t3["rhoVcp_star_l_over_inf"] * t2["rhoVcp_inf_kW_m2K"] * 1e3
    Taw, Tstar, Tl = t3["Taw_K"], t3["T_star_K"], t3["T_l_K"]
    Pr_s = gas.Pr(Tstar)
    eck = lambda R: 0.0296 * Pr_s ** (-2.0 / 3.0) * R ** -0.2

    # 感度 (a): 自前の T* (Tw = 300 K) で換算係数を組み直す
    Tst2 = Tl + 0.5 * (a.tw - Tl) + 0.22 * (Taw - Tl)
    Rfac2 = t2["Re_inf_per_m"] * (Tl / Tst2) * (gas.mu(Tl) / gas.mu(Tst2))
    rvc2 = t2["rhoVcp_inf_kW_m2K"] * 1e3 * (Tl / Tst2) * (gas.cp(Tst2) / gas.cp(Tl))
    Pr_s2 = gas.Pr(Tst2)

    with h5py.File(a.plate_h5, "r") as h:
        c = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        q = -np.asarray(h["VALUE/qwall"], float)
    m = (c[:, 0] > 0.0) & (np.abs(c[:, 1]) < 1e-9)
    x = np.round(c[m, 0], 9); ux = np.unique(x)
    qx = np.array([q[m][x == u].mean() for u in ux])
    St_f = qx / (rvc * (Taw - a.tw))
    R_f = Rfac * ux

    pts = [(p["R_star_l"], p["St_star_l"]) for p in d["points"]]
    rows = []
    print(f"Pr*({Tstar:.0f} K) = {Pr_s:.4f}   (ρVcp)*_l = {rvc:.0f} W/m²K   R*_l = {Rfac:.4e}·x   Taw {Taw} K  Tw {a.tw} K")
    print(f"{'x[m]':>6} {'R*_l':>10} {'St_exp':>9} {'St_Eck':>9} {'St_forge':>9} {'exp/Eck':>8} {'forge/Eck':>9} {'forge/exp':>9}")
    for R, se in pts:
        xx = R / Rfac
        sf = float(np.interp(xx, ux, St_f))
        sE = eck(R)
        # 感度 (b): 乱流起点をトリップ (x − 0.127) にした Eckert
        sEb = eck(Rfac * (xx - 0.127))
        rows.append(dict(x=xx, R=R, St_exp=se, St_eck=sE, St_eck_trip=sEb, St_forge=sf))
        print(f"{xx:6.3f} {R:10.3e} {se:9.3e} {sE:9.3e} {sf:9.3e} {se/sE:8.3f} {sf/sE:9.3f} {sf/se:9.3f}")
    r = lambda k1, k2: np.array([row[k1] / row[k2] for row in rows])
    summ = dict(
        exp_over_eck=dict(mean=float(r("St_exp", "St_eck").mean()), min=float(r("St_exp", "St_eck").min()), max=float(r("St_exp", "St_eck").max())),
        forge_over_eck=dict(mean=float(r("St_forge", "St_eck").mean()), min=float(r("St_forge", "St_eck").min()), max=float(r("St_forge", "St_eck").max())),
        forge_over_exp=dict(mean=float(r("St_forge", "St_exp").mean()), min=float(r("St_forge", "St_exp").min()), max=float(r("St_forge", "St_exp").max())),
        forge_over_eck_trip=dict(mean=float(r("St_forge", "St_eck_trip").mean())),
        exp_over_eck_trip=dict(mean=float(r("St_exp", "St_eck_trip").mean())),
    )
    # 感度 (a): 自前 T* の座標で forge と Eckert を比べる (実験の還元は Table III のまま)
    St_f2 = qx / (rvc2 * (Taw - a.tw)); R_f2 = Rfac2 * ux
    eck2 = lambda R: 0.0296 * Pr_s2 ** (-2.0 / 3.0) * R ** -0.2
    xs = np.array([row["x"] for row in rows])
    fe2 = np.interp(xs, ux, St_f2) / eck2(Rfac2 * xs)
    summ["forge_over_eck_ownTstar"] = dict(Tstar=Tst2, Pr=Pr_s2, mean=float(fe2.mean()), min=float(fe2.min()), max=float(fe2.max()))
    # x 全域の forge/Eckert
    sel = (ux > 0.05)
    fe_curve = St_f[sel] / eck(R_f[sel])
    summ["forge_over_eck_x"] = {f"{xx:.1f}": float(np.interp(xx, ux[sel], fe_curve)) for xx in (0.2, 0.5, 1.0, 1.5, 2.0, 2.5)}
    if a.forge_props and a.case_setup:
        FP = json.loads(Path(a.forge_props).read_text()); cs = json.loads(Path(a.case_setup).read_text())
        for key, fp in FP.items():
            Ts = float(key); rho = cs["p_inf"] / (cs["R"] * Ts); U = cs["U_inf"]
            qE = lambda xx: 0.0296 * fp["Pr"] ** (-2 / 3) * (rho * U * xx / fp["mu"]) ** -0.2 * rho * U * fp["cp"] * (Taw - a.tw)
            rr = np.interp(xs, ux, qx) / qE(xs)
            summ[f"forge_over_eck_forgeprops_T{key}"] = dict(Pr=fp["Pr"], mu=fp["mu"], cp=fp["cp"], mean=float(rr.mean()),
                                                             min=float(rr.min()), max=float(rr.max()))
    print("\n" + json.dumps(summ, ensure_ascii=False, indent=1))
    if a.json:
        step = max(1, len(ux) // 400)
        out = dict(meta=dict(Pr_star=Pr_s, Tstar=Tstar, Taw=Taw, Tw=a.tw, rhoVcp_star=rvc, Rfac=Rfac,
                             eckert="St*_l = 0.0296 Pr*^(-2/3) R*_l^(-1/5) (D-7275 p.14, ref. 12 Kays)"),
                   points=rows, summary=summ,
                   forge_curve=dict(x=ux[::step].tolist(), R=R_f[::step].tolist(), St=St_f[::step].tolist(),
                                    q_kW=(qx[::step] / 1e3).tolist()))
        Path(a.json).write_text(json.dumps(out, ensure_ascii=False, indent=1))
        print(f"-> {a.json}")


if __name__ == "__main__":
    main()
