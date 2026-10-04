#!/usr/bin/env python3
"""幅依存性の比較 — [W70] Fig 12 (開口面積平均 q̄_c/q_fp vs w/d) と forge の突き合わせ。

**なぜ比でなく絶対量を見るか**: 比 q̄_c/q_fp = Q_c/(w·ℓ·q_fp) は幅 w で割っている。
薄板の横方向伝導がリップから持ち込む熱は**幅に依存しない絶対量**なので、比で見ると
狭いすきまほど拡大される。`Q_c` [W/m] で並べれば、比例誤差 (ソルバのモデル誤差) と
加算誤差 (観測の混入) が分離できる。

前壁は [W70] と同じく Burggraf 理論で置換した (b) を一次比較量に使う (§4.4)。
"""
import json, subprocess, sys
from pathlib import Path
import numpy as np

CASE = Path(__file__).resolve().parents[1]
TOOLS = CASE.parents[1] / "solver_density_cuda" / "tools"

RUNS = {0.063: "run_0016_T1_wd0063_fine_settle", 0.211: "run_0008_T1_wd0211_long",
        0.383: "run_0009_T1_wd0383_long",        0.524: "run_0011_T1_wd0524_long"}

# 薄板 (304SS, W70 の薄肉過渡法) — case/50 tools/skin_smear.py と同じ値
ALPHA_S = 3.797e-6      # m²/s
T_MEAS = (0.4, 0.8, 1.6)  # s (報告から確定できないので幅で持つ)


def read_fig12():
    out = {}
    for ln in (CASE / "ref/w70_fig12_avg.csv").read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#") or ln[0].isalpha():
            continue
        k, v = ln.split(",")
        out[round(float(k), 3)] = float(v)
    return out


def verdict(run, quantity="Qc_front_theory_W_m"):
    """報告量の準定常判定。**一次比較量 (b) = 前壁置換後の熱量**の列を見る\n    (codex result-3 m1: 置換前の Qc で判定して置換後の量を報告していた)。"""
    csv = CASE / run / "cavity_series.csv"
    if not csv.exists():
        return "(cavity_series.csv なし — cavity_eval.py --series を先に回す)"
    r = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"),
                        "--series-csv", str(csv), "--series-cols", quantity],
                       capture_output=True, text=True)
    for l in r.stdout.splitlines():
        if quantity in l:
            return l.strip()
    return "(判定を取得できず)"


def main():
    lit = read_fig12()
    rows = []
    for wd, run in RUNS.items():
        j = json.loads((CASE / run / "cavity_eval.json").read_text(encoding="utf-8"))
        w, qfp = j["w"], j["q_lit"]
        Qb = j["Qc_per_span_front_theory"]
        Qm = lit[wd] * w * qfp
        rows.append(dict(wd=wd, run=run, w=w, qfp=qfp, Qb=Qb, Qm=Qm,
                         ratio_b=j["qbar_over_lit_front_theory"], ratio_lit=lit[wd]))

    print("[比で見たとき] — 狭いすきまほど食い違う")
    print(f"{'w/d':>6} {'w[mm]':>7} {'文献 q̄/q_fp':>12} {'forge (b)':>10} {'forge/文献':>11}")
    for r in rows:
        print(f"{r['wd']:6.3f} {r['w']*1e3:7.3f} {r['ratio_lit']:12.2f} {r['ratio_b']:10.3f} "
              f"{r['ratio_b']/r['ratio_lit']:11.3f}")

    Qf = np.array([r["Qb"] for r in rows]); Qm = np.array([r["Qm"] for r in rows])
    print("\n[絶対量で見たとき] — 食い違いは幅に依らない加算量")
    print(f"{'w/d':>6} {'Q_meas[W/m]':>12} {'Q_forge[W/m]':>13} {'差[W/m]':>9}")
    for r, a, b in zip(rows, Qm, Qf):
        print(f"{r['wd']:6.3f} {a:12.2f} {b:13.2f} {a-b:9.2f}")

    rms = lambda x: float(np.sqrt(np.mean(x ** 2)))
    k = float(np.sum(Qf * Qm) / np.sum(Qf * Qf))
    c = float(np.mean(Qm - Qf))
    slope, icpt = np.linalg.lstsq(np.c_[Qf, np.ones_like(Qf)], Qm, rcond=None)[0]
    print(f"\n  (1) 比例のみ   Q_meas = {k:.3f}·Q_forge            残差 RMS {rms(Qm-k*Qf):6.2f} W/m")
    print(f"  (2) 加算のみ   Q_meas = Q_forge + {c:.1f}          残差 RMS {rms(Qm-Qf-c):6.2f} W/m")
    print(f"  (3) 1 次       Q_meas = {slope:.3f}·Q_forge + {icpt:.1f}   残差 RMS {rms(Qm-(slope*Qf+icpt)):6.2f} W/m")
    print(f"  → 加算モデルが比例モデルより残差 {rms(Qm-k*Qf)/rms(Qm-Qf-c):.1f} 倍小さい。"
          f"傾きは {slope:.3f} (1 から {abs(1-slope)*100:.0f} %)。")

    # --- 4 点の回帰をどこまで信じてよいか (codex result-3 M1) ---
    # 「加算項は幅に依らない」と読むには、隣接幅の増分比が 1 付近で一定でなければならない。
    # 実際には単調に落ちるので、**差は幅に依存する**。傾きも 1 点除くだけで大きく動く。
    inc = (Qm[1:] - Qm[:-1]) / (Qf[1:] - Qf[:-1])
    print("\n[この回帰をどこまで信じてよいか] — 4 点しかないことの帰結")
    print("  隣接幅の増分比 ΔQ_meas/ΔQ_forge:")
    for i, v in enumerate(inc):
        print(f"     w/d {rows[i]['wd']:.3f} → {rows[i+1]['wd']:.3f} : {v:.3f}")
    print(f"  → {inc[0]:.3f} → {inc[-1]:.3f} と**単調に落ちる**。加算項が幅に依らないなら"
          f" 1 付近で一定のはずなので、**差は幅に依存する**。")
    print("  1 点除外の傾き:")
    los = []
    for kk in range(len(Qf)):
        m = np.ones(len(Qf), bool); m[kk] = False
        sl, ic = np.linalg.lstsq(np.c_[Qf[m], np.ones(int(m.sum()))], Qm[m], rcond=None)[0]
        los.append(float(sl))
        print(f"     w/d={rows[kk]['wd']:.3f} を除く: 傾き {sl:.3f}, 切片 {ic:+.1f} W/m")
    print(f"  → 傾きは {min(los):.3f}–{max(los):.3f} に動く。**「実測と 12 % 以内で一致」とは言えない**"
          f" (最も狭い点に依存している)。")
    print(f"  差そのもの {np.round(Qm-Qf,2)} W/m も単調でなく、最大は w/d={rows[int(np.argmax(Qm-Qf))]['wd']:.3f}。")
    print("  ** 結論: これは 4 点への事後的な回帰であって、検証精度ではない。")
    print("     測定・digitize・前壁理論補完・格子の誤差を傾きに伝播させ、独立条件で確かめるまで")
    print("     case/49 へ定量値として渡さない。")

    print("\n[加算量の出所の上界] — 薄板リップの横方向伝導 (両リップ) が供給しうる量")
    qfp_typ = float(np.mean([r["qfp"] for r in rows]))
    for t in T_MEAS:
        L = np.sqrt(ALPHA_S * t)
        print(f"   t = {t:.1f} s   L_diff = {L*1e3:.2f} mm   2·q_fp·L_diff = {2*qfp_typ*L:6.1f} W/m")
    L08 = np.sqrt(ALPHA_S * 0.8)
    print(f"  → 観測された加算量 {c:.1f} W/m は上界の {c/(2*qfp_typ*L08)*100:.0f} % (t=0.8 s)。"
          f" **供給可能な範囲に収まる**が、これは整合であって証明ではない"
          f" (リップの熱的接続は原報に書かれていない — §5.1 #32)。")

    print("\n[準定常判定] — 上の Q_forge を『結果』として使ってよいかの根拠")
    for r in rows:
        print(f"   w/d={r['wd']:.3f}  {r['run']}")
        print(f"      {verdict(r['run'])}")
    print("\n  注: 残差はいずれも NOT CONVERGED (プラトー)。報告量の準定常だけで使っている。")

    out = dict(rows=[{k2: v for k2, v in r.items()} for r in rows],
               fit_proportional=k, fit_offset=c, fit_slope=float(slope), fit_intercept=float(icpt),
               rms_proportional=rms(Qm - k * Qf), rms_offset=rms(Qm - Qf - c),
               increment_ratios=[float(x) for x in (Qm[1:]-Qm[:-1])/(Qf[1:]-Qf[:-1])],
               leave_one_out_slopes=[float(np.linalg.lstsq(
                   np.c_[Qf[np.arange(len(Qf)) != kk], np.ones(len(Qf)-1)],
                   Qm[np.arange(len(Qf)) != kk], rcond=None)[0][0]) for kk in range(len(Qf))],
               caveat="4 点への事後回帰。増分比が単調に落ちるので加算項は幅に依らない訳ではない。"
                      "検証精度として引用しないこと (codex result-3 M1)。",
               lip_conduction_bound_W_per_m={str(t): 2 * qfp_typ * float(np.sqrt(ALPHA_S * t))
                                             for t in T_MEAS})
    (CASE / "width_scaling.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                             encoding="utf-8")
    print(f"\n  → {CASE/'width_scaling.json'}")


if __name__ == "__main__":
    main()
