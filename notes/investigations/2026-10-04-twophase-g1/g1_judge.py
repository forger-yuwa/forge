#!/usr/bin/env python3
"""G1 (0 step の面作用素 A/B) の判定 — plans/active/condensation-two-phase-default.md §5.1 #4 の事前登録をそのまま適用する。

入力は forge の診断 D1 (環境変数 FORGE_DIAG_TP_FACES=<h5>) が書いた h5。
  - A = 二相拡散 OFF の作用素 (species_diffusion_d の診断用の写し; /off/*)
  - B = 二相拡散 ON の作用素 (本番の tp_build_face_in + tp_face_flux<float>; /on_f/*)、double 参照 (/on_d/*)
  符号: J は「セル ic0 へ入る向き」が正の面積分値 (h5 属性 sign_convention)。

判定 (事前登録; 結果を見て変えない。定数はこのファイルの先頭に固定):
  (i)  乱流域 (面の μt,f/μ_f ≥ 1) の全面で |J_w^A| ≤ 0.05·max|J_l^B|
  (ii) |Jl_float − Jl_double| ≤ 8ε₃₂·A_l、|j_v,float − j_v,double| ≤ 8ε₃₂·A_v (全評価面)。外れたら診断未成立
  (iii) 第 1 仮説の判別 (記録; 既定化のゲートではない): x 35〜70 mm、壁距離帯 0〜0.1・0.1〜0.4・0.4〜1.6 mm、上下壁を壁距離差と
        面の向きで分けて、壁法線方向の面積分和が 分子蒸気 = 壁から外向き・液乱流 = 壁向き、かつ誤差尺度の和を超える → 支持
前提 (判定の成立条件; 外れたら (i) は UNDETERMINED):
  P0 OFF の写しが本番 species_diffusion_d の拡散寄与と節点で一致 (float 和の丸め上界 nface·ε₃₂·Σ|J| 以内)
  P1 後処理が面計算の入力 (vis_turb を除く) を変えていない、ON の中間量の自己検査の不一致 0、OFF/ON の skip の不一致 0

本ツールが決めた実装上の選択 (plan に明記が無いもの; 呼び出し側の確認待ち):
  - (i) の max|J_l^B| は乱流域の面での最大
  - (iii) の面の壁距離 = 両端節点の wall_dist の算術平均、壁法線成分 = J·Δwd/|Δx_cc| (Δwd = wd1 − wd0, |∇wd| = 1 として
    2 点流束の向きの壁法線への射影)、上下壁 = 壁距離が増える向き (Δx_cc·sign(Δwd)) の y 成分の符号 (負 = 上壁、正 = 下壁)、
    誤差尺度の和 = Σ 8ε₃₂·A·|Δwd|/|Δx_cc|、総合は面のある全 (帯 × 壁) で支持のときだけ「支持」

使い方:
  python3 g1_judge.py <diag.h5> [--coord-scale 1.0] [--wall-dist-h5 IC.h5 --wall-dist-dataset VALUE/wall_dist]
終了コード: 0 = (i)(ii) とも PASS、1 = どちらかが FAIL、2 = どちらかが UNDETERMINED (FAIL が無いとき)。
"""
import argparse
import sys

import h5py
import numpy as np

# ---- 事前登録の定数 (plan §5.1 #4; 変えない) ----
EPS32 = 1.1920928955078125e-7
K_ERR = 8.0               # 8ε₃₂
TURB_RATIO = 1.0          # 乱流域: μt,f/μ_f ≥ 1
FRAC_I = 0.05             # (i) の 5 %
X_RANGE_MM = (35.0, 70.0)
WD_BANDS_MM = [(0.0, 0.1), (0.1, 0.4), (0.4, 1.6)]


def load(path):
    f = h5py.File(path, "r")
    d = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            d[name] = obj[()]
    f.visititems(visit)
    attrs = {k: f.attrs[k] for k in f.attrs}
    f.close()
    return d, attrs


def col(d, name, k=None):
    a = d[name]
    if k is not None:
        return a[:, k]
    return a


def verdict_line(tag, v, detail):
    print(f"VERDICT {tag}: {v}  — {detail}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5")
    ap.add_argument("--revised-4sr", action="store_true", help="plan #4sr の改訂後の誤差尺度でも (ii) を判定する (旧判定も併記)")
    ap.add_argument("--coord-scale", type=float, default=1.0, help="座標を m に直す倍率 (メッシュが m なら 1)")
    ap.add_argument("--wall-dist-h5", default=None, help="diag h5 に wall_dist が無いときの節点 wall_dist の出典 (IC h5 等)")
    ap.add_argument("--wall-dist-dataset", default="VALUE/wall_dist")
    args = ap.parse_args()

    d, at = load(args.h5)
    n = int(at["nSpecies"]); iw = int(at["iw"])
    names = [s.decode() if isinstance(s, bytes) else str(s) for s in at.get("species_names", [])]
    print(f"input {args.h5}: faces {int(at['nFaces'])}, nodes {int(at['nNodes'])}, nSpecies {n}, iw {iw} "
          f"({names[iw] if 0 <= iw < len(names) else '?'}), state {at.get('twophase_diffusion_state')}, "
          f"requested {at.get('condTwoPhaseDiffusion_requested')}, IC {at.get('input_value_file')}")
    print(f"sign: {at.get('sign_convention')}")

    # ---- 組立の前処理の変化 (記録) ----
    pn = [s.decode() if isinstance(s, bytes) else str(s) for s in d["pre/names"]]
    print("pre-part changes of conserved arrays (count real / max|d| real / count ghost / max|d| ghost):")
    for i, nm in enumerate(pn):
        print(f"  {nm:10s} {int(d['pre/count_real'][i]):8d} {d['pre/maxabs_real'][i]:.3e} "
              f"{int(d['pre/count_ghost'][i]):8d} {d['pre/maxabs_ghost'][i]:.3e}")

    undet = []   # (i) を UNDETERMINED にする理由

    # ---- P1 ----
    qn = [s.decode() if isinstance(s, bytes) else str(s) for s in d["post_inputs/names"]]
    bad_inp = [(nm, int(d["post_inputs/count_real"][i]), int(d["post_inputs/count_ghost"][i]))
               for i, nm in enumerate(qn) if nm != "vis_turb"
               and (d["post_inputs/count_real"][i] + d["post_inputs/count_ghost"][i]) > 0]
    if bad_inp:
        print(f"P1: post-part changed face inputs: {bad_inp}")
        undet.append("post-part changed face inputs")
    skip = col(d, "face/skip")
    ev = skip == 0
    n_odd = int(np.sum((skip == 2) | (skip == 3)))
    n_self = int(np.sum(col(d, "on_f/selfcheck_mismatch")[ev] != 0))
    print(f"faces: evaluated {int(ev.sum())}, node boundary half-faces {int(np.sum(skip == 1))}, "
          f"OFF/ON skip disagreement {n_odd}, intermediate self-check mismatch {n_self}")
    if n_odd or n_self:
        undet.append("skip disagreement or intermediate self-check mismatch")

    # ---- P0: OFF の写しと本番の拡散寄与 ----
    nf = np.maximum(col(d, "off_check/nface"), 1).astype(float)
    p0_ok = True
    print("P0: OFF copy vs production species_diffusion_d (node; bound = nface·eps32·Σ|J|):")
    checks = [(f"roY{s}" + ("(w)" if s == iw else ""), d["off_check/R_copy"][:, s], d["off_check/R_prod"][:, s].astype(float),
               d["off_check/R_absum"][:, s]) for s in range(n)]
    checks.append(("roe", d["off_check/roe_copy"], d["off_check/roe_prod"].astype(float), d["off_check/roe_absum"]))
    checks += [(f"diag{s}", d["off_check/diag_copy"][:, s], d["off_check/diag_prod"][:, s].astype(float),
                d["off_check/diag_copy"][:, s]) for s in range(n)]
    for nm, cp, pr, ab in checks:
        diff = np.abs(cp - pr)
        bound = nf * EPS32 * np.abs(ab)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(bound > 0, diff / bound, np.where(diff == 0, 0.0, np.inf))
        ok = bool(np.all(np.isfinite(cp)) and np.all(np.isfinite(pr)) and np.all(r <= 1.0))
        p0_ok &= ok
        print(f"  {nm:10s} max|copy-prod| {diff.max():.3e}  max|prod| {np.abs(pr).max():.3e}  max ratio {np.nanmax(r):.3f}  "
              f"over {int(np.sum(~(r <= 1.0)))}  {'ok' if ok else 'NG'}")
    if not p0_ok:
        undet.append("OFF copy does not match the production species_diffusion_d (P0)")

    # ---- (ii) float vs double ----
    Jl_f = col(d, "on_f/Jl")[ev].astype(float); Jl_d = col(d, "on_d/Jl")[ev]
    jv_f = col(d, "on_f/jv_mol")[ev].astype(float); jv_d = col(d, "on_d/jv_mol")[ev]
    Al = col(d, "scale/A_l")[ev]; Av = col(d, "scale/A_v")[ev]
    v2 = "UNDETERMINED"
    if ev.sum() == 0:
        detail2 = "no evaluated faces"
    else:
        fin = np.isfinite(Jl_f) & np.isfinite(Jl_d) & np.isfinite(jv_f) & np.isfinite(jv_d) & np.isfinite(Al) & np.isfinite(Av)
        el = np.abs(Jl_f - Jl_d); erv = np.abs(jv_f - jv_d)
        tl = K_ERR * EPS32 * Al; tv = K_ERR * EPS32 * Av
        okl = el <= tl; okv = erv <= tv
        with np.errstate(divide="ignore", invalid="ignore"):
            rl = np.where(tl > 0, el / tl, np.where(el == 0, 0.0, np.inf))
            rv = np.where(tv > 0, erv / tv, np.where(erv == 0, 0.0, np.inf))
        nbad_l = int(np.sum(~okl)); nbad_v = int(np.sum(~okv)); nnf = int(np.sum(~fin))
        v2 = "PASS" if (nbad_l == 0 and nbad_v == 0 and nnf == 0) else "FAIL"
        detail2 = (f"liquid: max |Jl_f−Jl_d|/(8ε A_l) {np.nanmax(rl):.3f}, faces over {nbad_l}; "
                   f"molecular vapour: max |jv_f−jv_d|/(8ε A_v) {np.nanmax(rv):.3f}, faces over {nbad_v}; non-finite {nnf}")
    verdict_line("(ii) float vs double [pre-registered]", v2, detail2)
    v2_used = v2
    if args.revised_4sr and ev.sum() > 0:
        # plan condensation-two-phase-default §5.1 #4sr: 結果を見た後の改訂 (旧判定は上の行に残す)。
        # 液流束の非正規化数の絶対項 E_abs = (3·|ctg| + 1)·2^-150 (FTZ なしの段階的アンダーフロー)
        ctg = (col(d, "on_in/ct")[ev].astype(float) * col(d, "on_in/geo")[ev].astype(float))
        eabs = (3.0 * np.abs(ctg) + 1.0) * 2.0 ** -150
        tl2 = tl + eabs
        okl2 = el <= tl2
        with np.errstate(divide="ignore", invalid="ignore"):
            rl2 = np.where(tl2 > 0, el / tl2, np.where(el == 0, 0.0, np.inf))
        nbad_l2 = int(np.sum(~okl2))
        v2r = "PASS" if (nbad_l2 == 0 and nbad_v == 0 and nnf == 0) else "FAIL"
        verdict_line("(ii') float vs double [revised #4sr, after seeing results]", v2r,
                     f"liquid: max |Jl_f−Jl_d|/(8ε A_l + E_abs) {np.nanmax(rl2):.3f}, faces over {nbad_l2}; "
                     f"molecular vapour unchanged (faces over {nbad_v}); non-finite {nnf}")
        v2_used = v2r
    if v2_used != "PASS":
        undet.append("(ii) double reference check did not pass (diagnostic not established)")

    # ---- (i) ----
    mu = col(d, "face/mu_f").astype(float); mut = col(d, "face/mut_f").astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        turb = ev & (mu > 0) & (mut / mu >= TURB_RATIO)
    JwA = col(d, "off/J", iw).astype(float)
    JlB = col(d, "on_f/Jl").astype(float)
    if turb.sum() == 0:
        v1, detail1 = "UNDETERMINED", "no turbulent faces (μt,f/μ_f ≥ 1)"
    else:
        mJl = float(np.max(np.abs(JlB[turb])))
        aw = np.abs(JwA[turb])
        lim = FRAC_I * mJl
        nover = int(np.sum(aw > lim))
        imax = int(np.argmax(aw))
        idx = np.nonzero(turb)[0][imax]
        sc = args.coord_scale
        detail1 = (f"turbulent faces {int(turb.sum())}; max|J_w^A| {aw.max():.3e} vs 0.05·max|J_l^B| {lim:.3e} "
                   f"(ratio {aw.max()/lim if lim > 0 else np.inf:.3f}); faces over {nover} ({nover/turb.sum():.3%}); "
                   f"max at face {int(col(d, 'face/ip')[idx])} x={col(d, 'face/pcx')[idx]*sc*1e3:.3f} mm "
                   f"y={col(d, 'face/pcy')[idx]*sc*1e3:.3f} mm; max|J_w^A| over all evaluated faces {np.abs(JwA[ev]).max():.3e}, "
                   f"max|J_l^B| over all evaluated faces {np.abs(JlB[ev]).max():.3e}")
        v1 = "PASS" if nover == 0 else "FAIL"
        if v1 == "PASS" and not np.all(np.isfinite(aw)):
            v1 = "FAIL"
        if undet:
            detail1 += "  [UNDETERMINED because: " + "; ".join(undet) + f"; the raw comparison was {v1}]"
            v1 = "UNDETERMINED"
    verdict_line("(i) |J_w^A| <= 0.05 max|J_l^B| on turbulent faces", v1, detail1)

    # ---- (iii) 第 1 仮説 (記録) ----
    print("RECORD (iii) hypothesis 1 (molecular vapour away from the wall, turbulent liquid toward the wall):")
    wd0 = col(d, "state0/wall_dist").astype(float); wd1 = col(d, "state1/wall_dist").astype(float)
    if not (np.all(np.isfinite(wd0[ev])) and np.all(np.isfinite(wd1[ev]))):
        if args.wall_dist_h5 is None:
            print("  wall_dist missing in the diag h5; pass --wall-dist-h5 (node field) — (iii) not evaluated")
            return finish(v1, v2)
        with h5py.File(args.wall_dist_h5, "r") as g:
            wdn = g[args.wall_dist_dataset][()].astype(float).ravel()
        ic0 = col(d, "face/ic0"); ic1 = col(d, "face/ic1")
        nn = int(at["nNodes"])
        wd0 = np.where(ic0 < nn, wdn[np.minimum(ic0, len(wdn) - 1)], np.nan)
        wd1 = np.where(ic1 < nn, wdn[np.minimum(ic1, len(wdn) - 1)], np.nan)
        print(f"  wall_dist from {args.wall_dist_h5}:{args.wall_dist_dataset}")
    sc = args.coord_scale
    x_mm = col(d, "face/pcx").astype(float) * sc * 1e3
    dx = (col(d, "state1/x") - col(d, "state0/x")).astype(float) * sc
    dy = (col(d, "state1/y") - col(d, "state0/y")).astype(float) * sc
    dz = (col(d, "state1/z") - col(d, "state0/z")).astype(float) * sc
    dl = np.sqrt(dx * dx + dy * dy + dz * dz)
    dwd = (wd1 - wd0) * sc
    wdm_mm = 0.5 * (wd0 + wd1) * sc * 1e3
    with np.errstate(divide="ignore", invalid="ignore"):
        cosn = np.where(dl > 0, dwd / dl, 0.0)       # 2 点流束の向き (ic0→ic1) の壁法線 (壁距離増加) 方向への射影
        ey = np.where(dwd != 0, dy * np.sign(dwd), 0.0)   # 壁距離が増える向きの y 成分
    side = np.where(ey < 0, "upper", np.where(ey > 0, "lower", "none"))
    # 壁から外向き (壁距離増加) を正にした壁法線成分: J はセル 0 へ入る向き = ic1→ic0 なので −J·cos
    jv_away_f = -col(d, "on_f/jv_mol").astype(float) * cosn
    jv_away_d = -col(d, "on_d/jv_mol") * cosn
    jl_away_f = -col(d, "on_f/Jl").astype(float) * cosn
    jl_away_d = -col(d, "on_d/Jl") * cosn
    sv = K_ERR * EPS32 * col(d, "scale/A_v") * np.abs(cosn)
    sl = K_ERR * EPS32 * col(d, "scale/A_l") * np.abs(cosn)
    inx = ev & (x_mm >= X_RANGE_MM[0]) & (x_mm <= X_RANGE_MM[1]) & (dwd != 0) & np.isfinite(wdm_mm)
    print(f"  faces in x {X_RANGE_MM[0]}–{X_RANGE_MM[1]} mm with Δwd ≠ 0: {int(inx.sum())} (wall side 'none': {int(np.sum(inx & (side == 'none')))})")
    print(f"  {'band [mm]':12s} {'wall':6s} {'faces':>6s} {'Σ j_v,mol away (f)':>19s} {'(d)':>11s} {'Σ scale_v':>10s} "
          f"{'Σ J_l away (f)':>15s} {'(d)':>11s} {'Σ scale_l':>10s}  result")
    allsup = True; any_cell = False
    for lo, hi in WD_BANDS_MM:
        last = (hi == WD_BANDS_MM[-1][1])
        inb = inx & (wdm_mm >= lo) & ((wdm_mm < hi) if not last else (wdm_mm <= hi))
        for sd in ("upper", "lower"):
            m = inb & (side == sd)
            if m.sum() == 0:
                print(f"  {lo:.1f}–{hi:.1f}{'':6s} {sd:6s} {0:6d}  (no faces)")
                continue
            any_cell = True
            Sv_f, Sv_d, Ssv = jv_away_f[m].sum(), jv_away_d[m].sum(), sv[m].sum()
            Sl_f, Sl_d, Ssl = jl_away_f[m].sum(), jl_away_d[m].sum(), sl[m].sum()
            sup = (Sv_f > Ssv) and (Sl_f < -Ssl)
            allsup &= sup
            print(f"  {lo:.1f}–{hi:.1f}{'':6s} {sd:6s} {int(m.sum()):6d} {Sv_f:19.4e} {Sv_d:11.4e} {Ssv:10.3e} "
                  f"{Sl_f:15.4e} {Sl_d:11.4e} {Ssl:10.3e}  {'支持' if sup else '棄却'}")
    if not any_cell:
        print("RECORD (iii): 判定不能 (対象の面が無い)")
    else:
        note = "" if v2 == "PASS" else " (ただし (ii) が PASS でないので診断未成立)"
        print(f"RECORD (iii): {'支持' if allsup else '棄却'}{note}")
    return finish(v1, v2)


def finish(v1, v2):
    vs = (v1, v2)
    if "FAIL" in vs:
        return 1
    if "UNDETERMINED" in vs:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
