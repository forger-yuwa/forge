#!/usr/bin/env python3
"""TN D-8233 のすきま薄肉 (304SS 0.406 mm) の横方向伝導を **順方向に解いて**、
実測 $h/h_{fp}$ がどれだけ嵩上げされているかを係数で出す。

[TH76] はすきま壁に 0.406 mm の 304SS 薄板を貼り、薄肉過渡法

    q_meas(z) = rho c tau (dT/dt)_z

で局所熱流束を出している (`conditions.json` の `model.skin`)。薄板の実際の式は

    dT/dt = q_conv(z)/(rho c tau) + alpha_s d2T/dz2,      alpha_s = lam_s/(rho c)

で、304SS の横方向拡散長は 1 s で **2.01 mm = すきま幅 W の 0.88 倍**。すきま内の
$q(z)$ は W の数分の一で減衰するので、**横方向伝導は上端の熱を深部へ再配分し、深部の
実測を嵩上げする**。生の $q_{conv}$ と実測を直接比べると forge が不当に低く見える
(plan §5.1 #24 の規定)。

ここでは case/50 と同じく **逆問題 (後退拡散) を解かない**。真の対流分布を仮定して
薄板過渡を解き、実測と同じ演算 $\\rho c\\tau\\,dT/dt$ を掛けて「見かけの熱流束」を作り、

    嵩上げ係数 f(z) = (真の対流 h) / (見かけの実測 h)

を出す。$f<1$ なら実測は嵩上げされている = **実測を下げる側**
(`tolerances.json` の `skin_smear.direction`)。$f$ は模型が線形なので $q_{fp}$ の
スケールに依らない。

**本ツールは係数を出すだけで、分類 (alpha/beta/gamma) は決めない** (親が `gap_eval.py` 側で行う)。

case/50 からの継承と相違:
  * 継承: 順方向モデル、リップ端を両極 (断熱 / 連続) で挟む、**時間更新と観測量で同じ
    境界演算子**を使う (2026-09-20 codex result-2 M8。異なる演算子だと全経路積分が
    1.017 倍ずれる)、リップ Dirichlet は表面が $q_{fp}$ を受けて温まる温度のランプ。
  * 相違: (a) case/50 は「後壁 + 床」を 1 本の経路にしたが (床座標が $y/d$ 基準という
    M8 の訂正もここ)、[TH76] は $D/W$=20.0 で床が拡散長の 10 倍以上先にあるため、
    経路は**下流壁 1 本 (リップ → 床隅)** で閉じる (床隅の温度上昇が 0 であることを
    毎回印字して確かめる)。(b) 入力は CFD の run ではなく、**実測プロファイル**または
    解析形 (`--profile forge`) を「真の分布」として与える。
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

# numpy 2.x には np.trapz が無い
_trapz = getattr(np, "trapezoid", None) or np.trapz

HERE = Path(__file__).resolve().parent
CASE = HERE.parent

# 304SS (既定。CLI で変更可)
LAM_S_DEF = 16.0        # 熱伝導率 [W/(m K)]
RHO_DEF = 7900.0        # 密度 [kg/m^3]
CP_DEF = 500.0          # 比熱 [J/(kg K)]
TAU_DEF = 0.406e-3      # 薄板厚 [m] (conditions.json の skin.thickness_cm 0.0406)

# acceptance.json の T2p-C が判定に使う深さ
JUDGE_ZW = (1.4, 2.2, 3.0)


# ---------------------------------------------------------------- 入力
def read_measured(path, t_ratio=1.00, wall="down"):
    """ref CSV から (z/W, h/h_fp) を取る。solid=1 (黒塗り記号) は除外する。"""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("T_ratio"):
            continue
        tr, w, zw, hr, solid = line.split(",")
        if abs(float(tr) - t_ratio) > 1e-9 or w != wall or int(solid) != 0:
            continue
        rows.append([float(zw), float(hr)])
    if not rows:
        raise SystemExit(f"REFUSED: {path} に T_ratio={t_ratio} / wall={wall} / solid=0 の行が無い")
    arr = np.array(sorted(rows))
    return arr


def profile_meas(zw, ref, tail_decay, interp="linear"):
    """実測点を補間し、最深点より先は exp(-tail_decay * dz/W) で外挿する。

    interp="linear"    : case/50 と同じ線形補間。**実測 6 点では折れ点に曲率が集中する**
                         ので、その点の係数は補間の作り方に敏感になる。
    interp="loglinear" : log q を線形補間 (指数減衰を区分的に繋ぐ)。折れの大きさが減る。
    どちらが真かは決められないので、両方を回して差を見る。
    """
    z0, q0 = ref[:, 0], ref[:, 1]
    if interp == "loglinear" and np.all(q0 > 0.0):
        q = np.exp(np.interp(zw, z0, np.log(q0)))
    else:
        q = np.interp(zw, z0, q0)
    deep = zw > z0[-1]
    q[deep] = q0[-1] * np.exp(-tail_decay * (zw[deep] - z0[-1]))
    return q


def profile_forge(zw, anchor, decay):
    """forge のような急峻な分布 (q ∝ exp(-decay * z/W))。"""
    return anchor * np.exp(-decay * zw)


# ---------------------------------------------------------------- 薄板過渡
def solve(s, qhat, t_end, alpha, rho_c_tau, lip="adiabatic", n_sub=None):
    """dT/dt = q/(rho c tau) + alpha T'' を陽解法で t_end まで解く。

    s[0] = リップ (開口)、s[-1] = 床隅。熱源 qhat は q/q_fp (無次元) で、
    温度も T/q_fp で持つ (模型は線形なので係数はスケールに依らない)。

      lip="adiabatic" : リップ端 断熱 (すきま薄板が断熱材 + 水冷通路で絶縁されている側の極)
      lip="continuous": リップ端が模型表面と連続。表面は q_fp を受けて温まるので、
                        同じ厚さの薄板として T_lip(t) = t/(rho c tau) * q_fp を Dirichlet で与える
                        (= リップからの伝導寄与を過小評価しない側の極)

    戻り値は (T/q_fp, q_apparent/q_fp)。**観測量は時間更新と同じ境界演算子**で作る
    (case/50 codex result-2 M8。違う演算子を使うと全経路積分が保存しない)。
    """
    ds = s[1] - s[0]
    nst = int(math.ceil(t_end / (0.4 * ds ** 2 / alpha))) if n_sub is None else n_sub
    dt = t_end / nst
    T = np.zeros_like(qhat)
    src = qhat / rho_c_tau
    lip_rate = 1.0 / rho_c_tau          # 表面 (q/q_fp = 1) の温度上昇率

    def laplace(T):
        lap = np.zeros_like(T)
        lap[1:-1] = (T[2:] - 2.0 * T[1:-1] + T[:-2]) / ds ** 2
        lap[-1] = 2.0 * (T[-2] - T[-1]) / ds ** 2       # 床隅側は断熱 (床は拡散長の外)
        if lip == "adiabatic":
            lap[0] = 2.0 * (T[1] - T[0]) / ds ** 2
        return lap

    for k in range(nst):
        T += dt * (src + alpha * laplace(T))
        if lip == "continuous":
            T[0] = lip_rate * (k + 1) * dt          # 表面と同じ速さで温まる

    lap = laplace(T)
    if lip == "continuous":
        lap[0] = lap[1]                 # Dirichlet 端では壁側の値を使わない
    return T, rho_c_tau * (src + alpha * lap)


# ---------------------------------------------------------------- 出力
def fmt(x, w=9, p=4):
    return f"{x:{w}.{p}f}" if np.isfinite(x) else f"{'n/a':>{w}}"


def main():
    ap = argparse.ArgumentParser(
        description="TN D-8233 すきま薄肉の横方向伝導による実測の嵩上げ係数 (順方向モデル)")
    ap.add_argument("--profile", choices=["meas", "forge"], default="meas",
                    help="真の対流分布。meas=実測プロファイルをそのまま真値とみなす / "
                         "forge=急峻な解析形 q ∝ exp(-decay z/W)")
    ap.add_argument("--ref", default="ref/th76_fig5a_re147.csv")
    ap.add_argument("--t-ratio", type=float, default=1.00, help="ref CSV の壁温比 (既定 1.00)")
    ap.add_argument("--wall", default="down", choices=["down", "up"])
    ap.add_argument("--times", default="0.25,1,2,5",
                    help="薄肉過渡法の評価時刻 [s] (既定 0.25-5 s。原典の還元手順は未確認)")
    ap.add_argument("--interp", choices=["linear", "loglinear"], default="linear",
                    help="--profile meas の測点間の繋ぎ方 (既定 linear = case/50 踏襲)")
    ap.add_argument("--tail-decay", type=float, default=math.pi,
                    help="実測最深点より深い側の外挿減衰率 (per z/W, 既定 pi)")
    ap.add_argument("--forge-anchor", type=float, default=None,
                    help="--profile forge の開口値 (既定は実測の z/W=0)")
    ap.add_argument("--forge-decay", type=float, default=math.pi,
                    help="--profile forge の減衰率 (per z/W, 既定 pi)")
    ap.add_argument("--lam", type=float, default=LAM_S_DEF, help="薄板の熱伝導率 [W/(m K)]")
    ap.add_argument("--rho", type=float, default=RHO_DEF, help="薄板の密度 [kg/m^3]")
    ap.add_argument("--cp", type=float, default=CP_DEF, help="薄板の比熱 [J/(kg K)]")
    ap.add_argument("--tau", type=float, default=TAU_DEF, help="薄板の板厚 [m]")
    ap.add_argument("--width", type=float, default=None, help="すきま幅 [m] (既定 geometry.json)")
    ap.add_argument("--depth", type=float, default=None, help="すきま深さ [m] (既定 geometry.json)")
    ap.add_argument("--n", type=int, default=1201, help="深さ方向の格子点数")
    a = ap.parse_args()

    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    W = a.width if a.width else geom["gap"]["width"] * 1e-3
    D = a.depth if a.depth else geom["gap"]["depth"] * 1e-3
    rho_c = a.rho * a.cp
    alpha = a.lam / rho_c
    rho_c_tau = rho_c * a.tau

    ref = read_measured(CASE / a.ref, a.t_ratio, a.wall)

    # 物理温度の印字用。係数そのものはスケールに依らない
    acc = json.loads((CASE / "acceptance.json").read_text(encoding="utf-8"))
    hfp = float(acc["case"]["h_fp"]["value_W_m2K"])
    t_aw = float(acc["case"]["freestream"]["T_aw_K"])
    t_w = float(acc["case"]["wall_K"])
    qfp = hfp * (t_aw - t_w)

    s = np.linspace(0.0, D, a.n)
    zw = s / W
    if a.profile == "meas":
        qhat = profile_meas(zw, ref, a.tail_decay, a.interp)
        pdesc = (f"meas ({a.ref}, T_ratio={a.t_ratio:.2f} / wall={a.wall} / solid=0 の "
                 f"{len(ref)} 点を {a.interp} 補間、z/W>{ref[-1,0]:.1f} は "
                 f"exp(-{a.tail_decay:.2f} dz/W) で外挿)")
    else:
        anchor = a.forge_anchor if a.forge_anchor is not None else float(ref[0, 1])
        qhat = profile_forge(zw, anchor, a.forge_decay)
        pdesc = (f"forge (解析形 q/q_fp = {anchor:.3f} exp(-{a.forge_decay:.3f} z/W)"
                 + ("" if a.forge_anchor is not None else "、開口値は実測 z/W=0 に合わせた") + ")")

    times = [float(x) for x in a.times.split(",")]
    lips = ["adiabatic", "continuous"]
    lip_ja = {"adiabatic": "断熱", "continuous": "連続"}

    print("=" * 96)
    print("TN D-8233 すきま薄肉の横方向伝導 — 順方向観測モデル (逆問題は解かない)")
    print("=" * 96)
    print(f"  すきま      W = {W*1e3:.3f} mm,  D = {D*1e3:.3f} mm,  D/W = {D/W:.2f}")
    print(f"  薄板 304SS  tau = {a.tau*1e3:.3f} mm,  lam = {a.lam:.1f} W/mK,  "
          f"rho = {a.rho:.0f},  c = {a.cp:.0f}  →  alpha = {alpha:.3e} m²/s")
    print(f"  真の分布    {pdesc}")
    print(f"  参考        q_fp = {qfp:.1f} W/m² (h_fp {hfp:.4f} × (T_aw {t_aw:.2f} - T_w {t_w:.1f}))"
          "  ※ 係数は線形模型なのでスケールに依らない")
    print("  横方向拡散長 sqrt(alpha t): " + ", ".join(
        f"t={t}s → {math.sqrt(alpha*t)*1e3:.2f} mm ({math.sqrt(alpha*t)/W:.2f} W)" for t in times))
    print("  **原典 (TN D-8233) の還元手順 (初期勾配か時間窓平均か・時間窓の長さ・伝導補正の有無) は未確認**。")
    print("  時刻は決め打ちにせず CLI (--times) とし、既定は case/50 と同じ 0.25-5 s を仮に置いている")
    print("  (tolerances.json の skin_smear._missing の規定どおり。取得できたら --times で差し替える)。")

    # ---- 8 隅を解く ---------------------------------------------------
    res = {}
    for t in times:
        for lip in lips:
            T, qapp = solve(s, qhat, t, alpha, rho_c_tau, lip=lip)
            res[(t, lip)] = (T, qapp)

    zs = sorted(set(list(ref[:, 0]) + list(JUDGE_ZW)))
    idx = {z: int(np.argmin(np.abs(zw - z))) for z in zs}

    print("\n嵩上げ係数 f(z) = (真の対流 h) / (見かけの実測 h)。")
    print("  f < 1 : 実測が嵩上げされている → 実測を下げる側 (比 forge/実測 の上端を上げる)")
    print("  f > 1 : 実測が削られている側")
    print("  [*] = acceptance.json T2p-C の判定深さ")
    head = f"{'t[s]':>6} {'リップ端':>8} " + "".join(
        f"{('z/W='+format(z,'.1f')+('*' if z in JUDGE_ZW else '')):>11}" for z in zs)
    print("\n" + head)
    print("-" * len(head))
    for t in times:
        for lip in lips:
            qapp = res[(t, lip)][1]
            print(f"{t:6.2f} {lip_ja[lip]:>8} " + "".join(
                f"{qhat[idx[z]]/qapp[idx[z]]:11.4f}" for z in zs))

    print("\n8 隅の最小・最大 (最小 = 実測を最も下げる側 = 比の上端に最も効く):")
    print(f"{'z/W':>7} {'f_min':>9} {'(時刻/端)':>18} {'f_max':>9} {'(時刻/端)':>18} "
          f"{'実測':>8} {'補正後 (f_min/f_max)':>24}")
    summary = {}
    for z in zs:
        vals = {(t, lip): qhat[idx[z]] / res[(t, lip)][1][idx[z]] for t in times for lip in lips}
        kmin = min(vals, key=vals.get)
        kmax = max(vals, key=vals.get)
        summary[z] = (vals[kmin], vals[kmax])
        meas = np.interp(z, ref[:, 0], ref[:, 1]) if ref[0, 0] <= z <= ref[-1, 0] else float("nan")
        tag = "*" if z in JUDGE_ZW else " "
        print(f"{z:6.1f}{tag} {vals[kmin]:9.4f} {f'{kmin[0]}s/{lip_ja[kmin[1]]}':>18} "
              f"{vals[kmax]:9.4f} {f'{kmax[0]}s/{lip_ja[kmax[1]]}':>18} "
              f"{fmt(meas, 8, 3)} {fmt(meas*vals[kmin], 11, 4)} {fmt(meas*vals[kmax], 11, 4)}")

    print("\n判定深さ (z/W = " + " / ".join(f"{z:.1f}" for z in JUDGE_ZW) + ") の区間:")
    for z in JUDGE_ZW:
        lo, hi = summary[z]
        print(f"  z/W {z:.1f}: f ∈ [{lo:.4f}, {hi:.4f}]   "
              f"→ 比 forge/実測 に掛かる係数 1/f ∈ [{1/hi:.3f}, {1/lo:.3f}]")

    # ---- 保存の自己チェック -------------------------------------------
    print("\n保存の自己チェック (両端断熱なら全経路の積分が厳密保存する — case/50 M8 と同じ性質):")
    I_true = _trapz(qhat, s)
    for t in times:
        qapp = res[(t, "adiabatic")][1]
        I_app = _trapz(qapp, s)
        print(f"  t={t:5.2f}s 断熱   ∫q_app ds / ∫q_true ds = {I_app/I_true:.12f}  "
              f"(相対差 {I_app/I_true-1:+.3e})")
    for t in times:
        qapp = res[(t, "continuous")][1]
        I_app = _trapz(qapp, s)
        print(f"  t={t:5.2f}s 連続   ∫q_app ds / ∫q_true ds = {I_app/I_true:.6f}  "
              f"(差はリップ端から流入する熱。保存しないのが正しい)")

    # ---- 経路の前提が成立しているか ------------------------------------
    print("\n経路の前提 (床隅 s=D が拡散長の外にあること):")
    for t in times:
        T = res[(t, "continuous")][0]
        print(f"  t={t:5.2f}s: リップ温度上昇 {T[0]*qfp:8.3f} K, "
              f"z/W=2.2 {T[idx[2.2]]*qfp:7.3f} K, 床隅 {T[-1]*qfp:.3e} K")
    print("  床隅の温度上昇が 0 なら、下流壁 1 本の経路で閉じてよい (床・上流壁と熱的に切れている)。")

    print("\n注意: 本ツールは係数を出すだけで、分類 (alpha/beta/gamma) は決めない。")
    print("      --profile meas は**既に嵩上げされた実測を真値とみなす**ので、嵩上げを過小評価する側。")
    print("      --profile forge は真の分布が急峻な場合で、嵩上げの上限側を見るためのもの。")
    print("      熱源 q_conv は計測窓中で一定とし、壁温上昇による q の低下は入れていない。")
    if a.profile == "forge":
        print("      --profile forge の深部は f が 1e-2 以下になりうる。これは『深部の見かけ熱流束は"
              "ほぼ全部が横方向伝導』という意味であって、補正係数として実測に掛けて使う値ではない"
              "(校正には使えない。嵩上げが上限どれだけ効きうるかの提示)。")
        print("      開口値を実測 z/W=0 に合わせているが、その実測自体が既に滲んだ値である"
              "(自己整合ではない)。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
