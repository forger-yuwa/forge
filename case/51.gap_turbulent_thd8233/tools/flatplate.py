#!/usr/bin/env python3
"""T2 Gate B — 分母 h_fp = C·Re'^0.69 を独立に検算する。

TN D-8233 の分母は**較正パネルの実測から作った相関**で、CFD からは直接出せない。
ここでは Van Driest II (Karman-Schoenherr + 圧縮性変換) と Reynolds 相似で
同じ条件の平板 h を独立に組み、相関値と桁・傾きが合うかを見る。合えば分母の定義
(r=0.89 の T_aw、Re' べき 0.69) を CFD 側で再現してよい、と言える。

すきまに近づく BL は**トンネル壁の厚い乱流 BL** (δ* = 10.5-12.8 cm = すきま幅の 46-56 倍)
なので、前縁から発達させるのではなく Fig 4 の θ を Re_θ として直接与える。
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_air import dry_air                                      # noqa: E402
# case/50 の tools も sys.path に入る (gas_air 経由) ので conditions は明示パスで読む
import importlib.util as _ilu                                    # noqa: E402
_sp = _ilu.spec_from_file_location("t2_conditions", HERE / "conditions.py")
_m = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_m)
Air = _m.Air

R_FACTOR = 0.89          # 乱流回復係数 (報告の定義)


def karman_schoenherr(Re_theta_i):
    L = np.log10(Re_theta_i)
    return 1.0 / (17.08 * L ** 2 + 25.11 * L + 6.012)


def van_driest_ii(M, T_inf, T_w, air, Re_theta, T_aw=None):
    """圧縮性 Cf を返す。Fc / F_theta は Hopkins-Inouye の表記。

    **T_aw は与えられたものを使う**。半完全気体では T_t は一定 gamma の式より 5 % 低いので、
    T_aw = T_inf(1 + r (gamma-1)/2 M^2) で作り直すと Crocco 二次式が端点で破綻する
    (T(xi=1) が T_inf に戻らない)。A^2, B は端点整合の形で書く。
    """
    if T_aw is None:
        gam = air.cp(T_inf) / (air.cp(T_inf) - air.R)
        T_aw = T_inf * (1.0 + R_FACTOR * 0.5 * (gam - 1.0) * M ** 2)
    A2 = (T_aw - T_inf) / T_w          # = (T_aw/T_inf - 1)(T_inf/T_w)
    B = (T_aw - T_w) / T_w
    A = np.sqrt(A2)
    den = np.sqrt(B ** 2 + 4.0 * A2)
    Fc = (T_aw / T_inf - 1.0) / (np.arcsin((2.0 * A2 - B) / den) + np.arcsin(B / den)) ** 2
    F_theta = air.mu(T_inf) / air.mu(T_w)
    Cf = karman_schoenherr(F_theta * Re_theta) / Fc
    return Cf, T_aw, Fc, F_theta


def h_vandriest(M, fs, T_w, air, theta):
    Re_theta = fs["ro"] * fs["U"] * theta / fs["mu"]
    Cf, T_aw, Fc, F_th = van_driest_ii(M, fs["T"], T_w, air, Re_theta, T_aw=fs.get("T_aw"))
    # 参照温度で Pr と cp を取り、Reynolds 相似係数 s = Pr^{-2/3}
    T_ref = fs["T"] + 0.5 * (T_w - fs["T"]) + 0.22 * (T_aw - fs["T"])
    Pr = air.Pr(T_ref)
    cp = air.cp(T_ref)
    St = 0.5 * Cf * Pr ** (-2.0 / 3.0)
    h = fs["ro"] * fs["U"] * cp * St
    return dict(h=h, Cf=Cf, St=St, Re_theta=Re_theta, T_aw=T_aw, T_ref=T_ref,
                Pr=Pr, cp=cp, Fc=Fc, F_theta=F_th)


def main():
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    der = json.loads((HERE.parent / "derived.json").read_text(encoding="utf-8"))
    M = cond["M"]
    C = cond["reference"]["C_digitized"]
    air = Air(dry_air(1000.0))
    # theta は **digitize 値ではなく bl_derived.json の導出値**を使う。
    # Fig 4 から読んだ theta は delta* と両立しない (delta が試験部を超える) ため棄却した。
    bld = HERE.parent / "bl_derived.json"
    if bld.exists():
        rows = json.loads(bld.read_text(encoding="utf-8"))["rows"]
        theta_of = {round(r["Re_m"] / 1e5): r["theta_cm_derived"] * 1e-2 for r in rows}
        print("  (theta は bl_derived.json の導出値を使用)\n")
    else:
        theta_of = {round(s["Re_m"] / 1e5): s["theta_cm"] * 1e-2
                    for s in cond["boundary_layer"]["series"]}

    print(f"=== T2 Gate B : 分母 h_fp = {C:.3e}·Re'^0.69 の独立検算 (Van Driest II) ===\n")
    rows = []
    for T_w in (300.0, 400.0, 480.0):
        print(f"--- 壁温 T_w = {T_w:.0f} K (T_w/T_aw = {T_w / 990:.2f}) ---")
        print(f"{'Re_inf':>11} {'Re_theta':>9} {'Cf':>9} {'h_VDII':>9} {'h_corr':>9} {'比':>7}")
        for s in der["series"]:
            Re = s["Re_m"]
            # θ は BL 系列から最も近い Re のものを使う
            key = min(theta_of, key=lambda k: abs(k - round(Re / 1e5)))
            theta = theta_of[key]
            fs = dict(ro=s["ro_inf"], U=s["U_inf"], mu=s["mu_inf"], T=s["T_inf"], T_aw=s["T_aw"])
            r = h_vandriest(M, fs, T_w, air, theta)
            h_corr = C * Re ** 0.69
            ratio = r["h"] / h_corr
            print(f"{Re:11.3e} {r['Re_theta']:9.3e} {r['Cf']:9.3e} "
                  f"{r['h']:9.3f} {h_corr:9.3f} {ratio:7.3f}")
            rows.append(dict(T_w=T_w, Re_m=Re, h_vdii=r["h"], h_corr=h_corr,
                             ratio=ratio, Cf=r["Cf"], Re_theta=r["Re_theta"],
                             T_aw=r["T_aw"], T_ref=r["T_ref"]))
        print()

    # Re' に対する冪指数も比べる (相関は 0.69)
    for T_w in (300.0, 400.0, 480.0):
        sub = [r for r in rows if r["T_w"] == T_w]
        p = np.polyfit(np.log([r["Re_m"] for r in sub]),
                       np.log([r["h_vdii"] for r in sub]), 1)[0]
        print(f"T_w={T_w:.0f} K : Van Driest II の実効冪 d(ln h)/d(ln Re') = {p:.3f}  "
              f"(相関 0.69)")

    rs = np.array([r["ratio"] for r in rows])
    print(f"\nh_VDII / h_corr : {rs.min():.3f} - {rs.max():.3f}  (中央 {np.median(rs):.3f})")
    verdict = "PASS" if 0.5 < np.median(rs) < 2.0 else "FAIL"
    print(f"GATE B VERDICT: {verdict}  (独立理論が相関の 0.5-2.0 倍に入るか)")
    print("  注: 相関は較正パネルの**実測**であり理論ではない。桁と傾きが合えば分母の"
          "定義 (r=0.89 の T_aw, Re' べき) を CFD 側で再現してよい、という確認。")

    (HERE.parent / "gateB.json").write_text(json.dumps(
        dict(_gate="B", C=C, verdict=verdict, rows=rows), indent=2), encoding="utf-8")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
