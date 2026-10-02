#!/usr/bin/env python3
"""二相拡散 A/B (plan condensation-two-phase-transport §5.1 #1b) の差を出す (記録用; 合否は付けない)。

    python3 twophase_ab_compare.py RUN_A RUN_B [--step N] [--gmin 1e-6] [--out PREFIX]

出すもの (同じ step のスナップショット同士。既定は両 run に共通の最大 step):
  1. 共通凝縮域 C = {g_A > gmin かつ g_B > gmin の節点} (既定 gmin 1e-6) での |ΔT| = |T_B − T_A| の体積重み分位点
     (p50 / p95 / p99 / max) と体積重み平均の符号付き ΔT。体積は VALUE/volume (A の値; メッシュは同一であることを検査する)。
     片方だけ g > gmin の体積 (onset 前縁・後縁のずれ) も出す。
     p95|ΔT| と 1 K の比較は**影響の大小の分類**であって合否ではない (plan §5.1 #1b)。
  2. 報告量 (twophase_ab_series.py の metrics(); 定義はそちらの docstring) の A・B・差・相対差の表。
     M0 (共通初期場の凝縮域) は A の valueFileName から作る (A・B は同じ初期場なので、B の値ファイルと一致するかも検査する)。

**共通初期場での EOS 投影 δṪ = (∂T/∂U)(R_B − R_A)/V は本スクリプトでは出さない (省略)**。理由:
  - res_*.h5 の step 0 出力は残差の組立て前に書かれる (main.cpp の writeInitialOutputs → twoPhaseAudit(0) → 時間ループの順) ので、
    0 step 出力に残差は無い。定常陰解法は 1 step に 1 回だけ残差を組み (implicitNonlinearUpdate → assembleResidual)、
    DPLUR は残差を固定して掃くので、1 step 回して res_1.h5 に output.extraFields で res_ro / res_roe / res_roY* / res_rog_0 を出せば
    「共通初期場での R」になる見込みはある。ただし (a) step 内の更新 (化学種・二相の非分割更新・ピン・周期集約) が res_* を
    上書きしないことは未確認、(b) ∂T/∂U は二相 EOS (e = e_v(Y_w 込み全蒸気) + g(R_w T − L(T)), L は液相多項式と種 DB の datum から)
    の接線であり、Python 側にソルバと同じ L(T)・datum・c_v,eff を持つ実装が無い (近似式では潜熱項の打ち消し — codex 2026-09-27 — を
    再現できず、量そのものが誤る)。(b) をソルバの host 関数で組む C++ ツール (#4c ハーネスの eos_T と同じ方式) が要るので、#1b の
    範囲では重い。定常の場の差 (1.) で影響を記録する。
"""
import argparse, importlib.util, os, sys
import h5py, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("twophase_ab_series", os.path.join(HERE, "twophase_ab_series.py"))
S = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(S)

QUANT = [("onset_c_g1e4_mm", "mm"), ("onset_c_g1e3_mm", "mm"), ("g_exit_mw", "-"), ("g_exit_c", "-"),
         ("pw_mean_x10", "p/p0"), ("pw21", "p/p0"), ("pw42", "p/p0"), ("pw52", "p/p0"), ("dev_pct", "%"),
         ("Tw_mean_x10_K", "K"), ("T_cond0_vmean_K", "K"), ("g_cond0_vmean", "-"), ("gmax", "-"), ("nonfinite", "count")]


def wquantile(x, w, q):
    """体積重み分位点: 昇順に並べた累積重みが q·Σw に初めて達する値。"""
    o = np.argsort(x); x, w = x[o], w[o]
    c = np.cumsum(w)
    return float(x[np.searchsorted(c, q*c[-1], side="left")])


def same_values(f1, f2):
    with h5py.File(f1, "r") as a, h5py.File(f2, "r") as b:
        va, vb = a["VALUE"], b["VALUE"]
        if set(va) != set(vb):
            return False, f"データセットの集合が違う: {sorted(set(va) ^ set(vb))}"
        bad = [k for k in va if not np.array_equal(np.asarray(va[k]), np.asarray(vb[k]))]
        return (not bad), (f"値が違う: {bad}" if bad else "全データセット一致")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_a"); ap.add_argument("run_b")
    ap.add_argument("--step", type=int, default=None, help="比べる step (既定: 両 run に共通の最大 step)")
    ap.add_argument("--gmin", type=float, default=1e-6, help="共通凝縮域の g 閾値 (既定 1e-6)")
    ap.add_argument("--out", default=None, help="表を PREFIX.md にも書く")
    a = ap.parse_args()

    fa = dict(S.res_files(a.run_a)); fb = dict(S.res_files(a.run_b))
    common = sorted(set(fa) & set(fb))
    if not common:
        sys.exit("両 run に共通の step の res_*.h5 が無い")
    st = a.step if a.step is not None else common[-1]
    if st not in fa or st not in fb:
        sys.exit(f"step {st} が両方に無い (共通: {common[:3]}..{common[-3:]})")

    lines = []
    def out(s=""):
        print(s); lines.append(s)

    # 初期場の同一性 (A・B の valueFileName)
    ok_ic, why_ic = same_values(S.value_file(a.run_a), S.value_file(a.run_b))
    _, geo, mask0, coord0 = S.setup(a.run_a, None, 1e-6, first=fa[st])
    da, db = S.read(fa[st]), S.read(fb[st])
    if not (np.array_equal(da["coord"], coord0) and np.array_equal(db["coord"], coord0)):
        sys.exit("A と B の MESH/COORD が違う (同一メッシュの A/B ではない)")

    out(f"# 二相拡散 A/B の差 (記録用; 合否なし)  step {st}")
    out(f"- A = `{a.run_a.rstrip('/')}` ({os.path.basename(fa[st])}), B = `{a.run_b.rstrip('/')}` ({os.path.basename(fb[st])})")
    out(f"- 初期場 (valueFileName の VALUE) の一致: {'OK' if ok_ic else 'NG'} — {why_ic}")
    if not ok_ic:
        out("  **初期場が違うので共通初期場からの A/B として読めない**")

    # 1. 共通凝縮域の |ΔT|
    vol = da["volume"]
    ga, gb = da["g_0"], db["g_0"]
    cm = (ga > a.gmin) & (gb > a.gmin)
    onlyA = (ga > a.gmin) & ~(gb > a.gmin); onlyB = (gb > a.gmin) & ~(ga > a.gmin)
    dT = db["T"] - da["T"]
    out("")
    out(f"## 1. 共通凝縮域 C (g_A > {a.gmin:g} かつ g_B > {a.gmin:g}) の温度差 ΔT = T_B − T_A")
    out(f"- C: {int(cm.sum())} 節点, 体積 {vol[cm].sum():.6e} (全体の {vol[cm].sum()/vol.sum()*100:.2f} %)")
    out(f"- A のみ凝縮: {int(onlyA.sum())} 節点・体積 {vol[onlyA].sum():.3e} / B のみ凝縮: {int(onlyB.sum())} 節点・体積 {vol[onlyB].sum():.3e}")
    if cm.sum() == 0:
        out("- C が空 (比較不能)")
    else:
        x, w = np.abs(dT[cm]), vol[cm]
        p50, p95, p99 = (wquantile(x, w, q) for q in (0.50, 0.95, 0.99))
        mean_signed = float(np.sum(dT[cm]*w)/np.sum(w))
        imax = np.where(cm)[0][np.argmax(x)]
        out("")
        out("| 量 | 値 [K] |")
        out("| --- | --- |")
        out(f"| 体積重み p50 \\|ΔT\\| | {p50:.4g} |")
        out(f"| **体積重み p95 \\|ΔT\\|** | **{p95:.4g}** |")
        out(f"| 体積重み p99 \\|ΔT\\| | {p99:.4g} |")
        out(f"| max \\|ΔT\\| | {x.max():.4g} (x = {geo.xy[imax,0]*1e3:.2f} mm, y = {geo.xy[imax,1]*1e3:.3f} mm) |")
        out(f"| 体積重み平均 ΔT (符号付き) | {mean_signed:+.4g} |")
        out("")
        out(f"- 分類 (合否ではない): p95|ΔT| {'≥' if p95 >= 1.0 else '<'} 1 K → 影響 {'大' if p95 >= 1.0 else '小'} の側")
        dg = (gb - ga)[cm]
        out(f"- 参考: C 上の体積重み平均 Δg = {np.sum(dg*w)/np.sum(w):+.4e}, 体積重み p95 |Δg| = {wquantile(np.abs(dg), w, 0.95):.4e}")

    # 2. 報告量の差
    ma, mb = S.metrics(da, geo, mask0), S.metrics(db, geo, mask0)
    out("")
    out("## 2. 報告量 (定義は twophase_ab_series.py)")
    out("")
    out("| 量 | 単位 | A | B | B − A | (B − A)/\\|A\\| |")
    out("| --- | --- | --- | --- | --- | --- |")
    for q, u in QUANT:
        va, vb = ma[q], mb[q]
        d = vb - va
        rel = d/abs(va) if (np.isfinite(va) and va != 0) else np.nan
        out(f"| {q} | {u} | {va:.7g} | {vb:.7g} | {d:+.4g} | {rel:+.3e} |")
    out("")
    out("- 壁熱流束: 断熱壁で出力なし (列を作らない)。拡散の差は Tw_mean_x10_K (壁の回復温度) に出る。")
    out("- 共通初期場での EOS 投影 δṪ: 省略 (理由は本スクリプトの docstring)。")
    out("- 前提: 両 run の check_convergence PASS・報告量 STEADY (check_quasisteady --series-csv)・B の `[twophase-audit] VERDICT: PASS` は"
        " 本スクリプトでは確かめない。")
    if a.out:
        open(a.out + ".md", "w").write("\n".join(lines) + "\n")
        print(f"\n-> {a.out}.md")


if __name__ == "__main__":
    main()
