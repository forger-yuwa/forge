"""plan condensation-two-phase-default §5.1 #4pj: 射影の打ち消しの原因の判別 A/B (液・Q の DPLUR 分母の共通化) の事前登録判定。
  python3 pj_ab_judge.py <A tp_update.h5> <B tp_update.h5> [--length-scale 1000] [--ref-q1 0.913] [--ref-q2 0.974] [--ref-tol 0.05]

入力 (どちらも FORGE_DIAG_TP_UPDATE=<h5> で forge が書いた 1 更新の記録; main.cpp runTpUpdateDiag):
  A = 本番の成分別分母 (FORGE_DIAG_TP_UPDATE だけ; h5 属性 common_diag = 0)
  B = 液・Q2・Q1・Q0 の分母を節点ごとの max(D_g, D_Q2, D_Q1, D_Q0) に共通化 (加えて FORGE_DIAG_TP_COMMON_DIAG=1; common_diag = 1)
  どちらも同じ再開場 (run_0561 の res_48000 を restart_field --force-species で写したもの) と同じ config から。
判断の出典: notes/reviews/2026-10-04-twophase-projection-fight-diagnose.md (判別 A/B)。

集合 (更新のバッファ、二相 ON の経路 TPC → RNS → PFL → RZL → RZU → RZP → RZR):
  射影集合 P = 射影の分岐に入った (RZP_input/kind ≥ 0 ⇔ ρQ0 > 0 かつ ρg > 0) かつ射影が Q1 または Q2 を変えた節点
              (condensationRealizability_d.cuh の float 実体; kind 3 [Q3 が表現できず Q1 = Q2 = 0] も P に入る)。
  塵集合 G0 = 下限・上限の後に ρg = 0 で ρQ0・ρQ1・ρQ2 のどれかが正の節点 (射影は通らず、後段の液滴消滅 RZR が扱う)。
  両者は排他 (P は ρg > 0)。G0 の補正は RZR の after − before (物理の消滅処理と同じ操作) として別に出す。

指標 1 (射影比、Ω ごと・成分ごと): ΣV|C_projection| / ΣV|δq_limited|
  C_projection = RZP の after − before (Q1・Q2 だけが動く; g・Q0 は射影で不変)、δq_limited = 二相 commit (TPC) の制限後増分 (delta_d)。
  分母は Ω の全節点 (g3b_judge.py の G3B_RESULT と同じ定義; 判定に使う)。参考に分母を P の節点に限った比も出す。
指標 2 (外向き成分):
  尺度付きモーメント空間: 更新開始状態 s (TPC の before) の Q0_s と Q3_s = ρg_s/((4/3)π ρ_l) から r_s = (Q3_s/Q0_s)^{1/3} を作り、
    m_k = Q_k/(Q0_s r_s^k) (k = 0..3)。開始状態は m_s = (1, x_s, y_s, 1) (x, y はクランプの cond_moment_xy と同じ無次元量)。
    ρ_l はクランプが射影に使った値 (RZP_input/rho_l = cond_rho_cond(T)) を開始状態にも使う (同じ節点の 1 更新内)。
  Hankel 条件 (クランプの許容領域 {x² ≤ y ≤ √x} と同値; Q0・Q3 固定なら H1 = y − x²、H2 = x − y²):
    H1 = m0 m2 − m1² ≥ 0、∇H1 = (m2, −2 m1, m0, 0)
    H2 = m1 m3 − m2² ≥ 0、∇H2 = (0, m3, −2 m2, m1)
  増分の方向: δm = (δQ0/Q0_s, δQ1/(Q0_s r_s), δQ2/(Q0_s r_s²), δg/ρg_s) (δ は TPC の delta_d)、s = |δm|₂。
  方向微分 dH_k/ds = ∇H_k(m_s)·δm / |δm|₂ (開始状態での 1 次変化を増分の大きさで正規化した量; 無次元)。
  境界の同定 (どの条件の外へ出たか): 射影の入力 (RZP の before; Q0 と g は RZL/RZU の後) を同じ ρ_l で無次元化し、クランプと同じ
    相対許容 eps = 1e-6 で H1 違反 y < x²(1 − eps)、H2 違反 y² > x(1 + eps)。A と B は同じ開始状態なので、節点ごとの
    「該当条件」= A で違反 ∪ B で違反 (両腕で同じ集合を使う)。
  評価集合 S = A の射影集合 P_A (A で射影が働いた節点; 両腕で同じ節点・同じ開始状態)。B で射影が起きなくなった節点でも
    B の増分の向きを同じ節点で評価する (B 自身の射影集合で数えると「射影が起きた節点は外へ出た節点」なので割合が常に ~1 になり、
    解消を測れない)。
  外向き (腕 X): 該当条件のどれかで dH_k/ds(δ_X) < −tol (tol = --tangent-tol、既定 1e-6) の節点。|dH_k/ds| ≤ tol (接線) の節点は
    有限増分での厳密な変化 ΔH_k = H_k(m_s + δm) − H_k(m_s) (二次まで厳密) の符号で判定する (単分散の角 x = y = 1 では両 Hankel 条件とも
    一様 ṙ の方向に 1 次の変化が 0 で、−h² の 2 次で外へ出る [diagnose の第 2 仮説の反例] ので、1 次だけでは丸め誤差の符号で決まってしまう)。
    接線で 2 次判定した節点数は別に出す。判定不能 (開始状態に液滴が無い [Q0_s ≤ 0 または ρg_s ≤ 0]、
    腕 X の増分が 0、該当条件が同定できない [kind 2 で x, y がともに 0 / kind 3 の Q3 表現不能]) は別に数え、
    外向き割合では**外向き側に数える** (B で「解消」を過大に言わない側の保守的な扱い)。参考に判定不能を除いた割合も出す。
  外向き割合_X = |外向き_X ∪ 判定不能_X| / |S|。体積加重平均 = Σ_{判定可能} V·min_{該当 k} dH_k/ds / Σ V (負 = 外向き)。
  参考: 各腕自身の射影集合 P_X での同じ量も出す (判定には使わない)。

事前登録の判定 (Ω = x 35–70 mm・壁距離 < 0.4 mm):
  前提 (A の再現): A の射影比 Q1・Q2 が参照値 (run_0580 の G3B_RESULT: Q1 0.913、Q2 0.974) から ±ref-tol 以内。
             再現しなければ UNDECIDED (A が前回と違う更新をしているので B との比が意味を持たない)。
  SUPPORT (第 1 仮説 = 液と Q の分母の差が主因): B の Q1・Q2 の射影比がともに A の 0.1 倍以下 かつ
             外向き成分の解消 (B の外向き割合 ≤ 0.1 × A の外向き割合)。
  REJECT  (対角差の主因説を棄却 → 第 2 仮説を優先): B の Q1・Q2 の射影比と外向き割合がすべて A の 0.5 倍以上。
  それ以外 UNDECIDED。射影の絶対量の減少だけでは支持としない (比で判定する)。
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g3b_judge as g3   # noqa: E402  (Rec・regions・ba・commit_x を共用)

C_G, C_Q2, C_Q1, C_Q0 = 1, 2, 3, 4
QCOMP = [("g", C_G), ("Q2", C_Q2), ("Q1", C_Q1), ("Q0", C_Q0)]
FOUR_THIRDS_PI = 4.0 / 3.0 * np.pi
EPS_PROJ = 1.0e-6          # cond_realizability_project の相対許容
VERDICT_REGION = "x35-70 wd<0.4"


def slot(rec, name):
    if name not in rec.idx:
        raise SystemExit(f"{rec.path}: slot '{name}' not found (h5 written by a forge without the #4pj slots?)")
    return rec.s("upd", name)


def ba_req(rec, short, c):
    r = g3.ba(rec, "upd", short, c)
    if r is None:
        raise SystemExit(f"{rec.path}: operation {short} component {g3.COMPS[c]} not on the recorded path")
    return r


def analyse(rec):
    """節点ごとの量 (dict)。"""
    if not rec.on:
        raise SystemExit(f"{rec.path}: two-phase OFF record; #4pj needs the ON path")
    out = {}
    # 開始状態と制限後増分 (二相 commit)
    for nm, c in QCOMP:
        b, _a = ba_req(rec, "TPC", c)
        out["s_" + nm] = b
        out["d_" + nm] = g3.z(g3.commit_x(rec, "TPC", c, "delta_d"))
    # 射影 (Q1・Q2) と射影の入力
    for nm, c in (("Q2", C_Q2), ("Q1", C_Q1)):
        b, a = ba_req(rec, "RZP", c)
        out["pin_" + nm] = b
        out["C_" + nm] = g3.z(a - b)
    out["C_g"] = np.zeros(rec.n); out["C_Q0"] = np.zeros(rec.n)
    _b, gu = ba_req(rec, "RZU", C_G)          # 下限・上限の後の液
    out["pin_g"] = gu
    _b, q0l = ba_req(rec, "RZL", C_Q0)        # 下限の後の Q0 (射影は Q0 を変えない)
    out["pin_Q0"] = q0l
    kind = slot(rec, "RZP_input/kind")
    rhol = slot(rec, "RZP_input/rho_l")
    out["kind"] = kind; out["rho_l"] = rhol
    if not np.isfinite(kind).any():
        raise SystemExit(f"{rec.path}: RZP_input/kind is all NaN (the float realizability clamp did not run in the update)")
    changed = (out["C_Q1"] != 0.0) | (out["C_Q2"] != 0.0)
    P = (kind >= 0) & changed
    out["P"] = P
    # 塵集合 (RZR の before で判定; 射影は通らない)
    rb = {nm: ba_req(rec, "RZR", c) for nm, c in QCOMP}
    G0 = (rb["g"][0] == 0.0) & ((rb["Q0"][0] > 0) | (rb["Q1"][0] > 0) | (rb["Q2"][0] > 0))
    out["G0"] = G0
    for nm, _c in QCOMP:
        out["R_" + nm] = np.where(G0, g3.z(rb[nm][1] - rb[nm][0]), 0.0)
    out["P_nochange_branch"] = (kind >= 0) & ~changed

    # 外向き成分の材料 (全節点; 集合はあとで選ぶ)
    with np.errstate(all="ignore"):
        g_s, q0_s, q1_s, q2_s = out["s_g"], out["s_Q0"], out["s_Q1"], out["s_Q2"]
        okS = (g_s > 0) & (q0_s > 0) & np.isfinite(rhol) & (rhol > 0)
        q3_s = g_s / (FOUR_THIRDS_PI * rhol)
        lr = (np.log(q3_s) - np.log(q0_s)) / 3.0                  # ln r_s (指数分離; クランプの cond_moment_xy と同じ組み方)
        x_s = np.where(q1_s > 0, np.exp(np.log(np.maximum(q1_s, 1e-300)) - np.log(q0_s) - lr), 0.0)
        y_s = np.where(q2_s > 0, np.exp(np.log(np.maximum(q2_s, 1e-300)) - np.log(q0_s) - 2 * lr), 0.0)
        dm0 = out["d_Q0"] / q0_s
        dm1 = np.exp(np.log(np.abs(out["d_Q1"]) + 1e-300) - np.log(q0_s) - lr) * np.sign(out["d_Q1"])
        dm2 = np.exp(np.log(np.abs(out["d_Q2"]) + 1e-300) - np.log(q0_s) - 2 * lr) * np.sign(out["d_Q2"])
        dm3 = out["d_g"] / g_s
        sn = np.sqrt(dm0**2 + dm1**2 + dm2**2 + dm3**2)
        okD = okS & np.isfinite(sn) & (sn > 0)
        dH1 = (y_s * dm0 - 2 * x_s * dm1 + 1.0 * dm2) / sn
        dH2 = (1.0 * dm1 - 2 * y_s * dm2 + x_s * dm3) / sn
        H1_s = y_s - x_s**2; H2_s = x_s - y_s**2
        # 有限増分での厳密な変化 (二次まで; m0 = m3 = 1 の展開形で桁落ちを避ける)。1 次の項が消える接線方向 (単分散の角 x = y = 1 など) の判定に使う
        fH1 = (y_s * dm0 - 2 * x_s * dm1 + dm2) + dm0 * dm2 - dm1**2
        fH2 = (dm1 - 2 * y_s * dm2 + x_s * dm3) + dm1 * dm3 - dm2**2
        # 射影の入力での違反条件 (射影の分岐に入った節点だけ; kind 3 は x, y を作れないので同定しない)
        q3_p = out["pin_g"] / (FOUR_THIRDS_PI * rhol)
        lrp = (np.log(q3_p) - np.log(out["pin_Q0"])) / 3.0
        x_p = np.where(out["pin_Q1"] > 0, np.exp(np.log(np.maximum(out["pin_Q1"], 1e-300)) - np.log(out["pin_Q0"]) - lrp), 0.0)
        y_p = np.where(out["pin_Q2"] > 0, np.exp(np.log(np.maximum(out["pin_Q2"], 1e-300)) - np.log(out["pin_Q0"]) - 2 * lrp), 0.0)
        br = (kind >= 0) & (kind != 3) & np.isfinite(x_p) & np.isfinite(y_p)
        v1 = br & (y_p < x_p**2 * (1 - EPS_PROJ))
        v2 = br & (y_p**2 > x_p * (1 + EPS_PROJ))
    out.update(dict(okS=okS, okD=okD, dH1=dH1, dH2=dH2, fH1=fH1, fH2=fH2, H1_s=H1_s, H2_s=H2_s, v1=v1, v2=v2))
    return out


TANGENT_TOL = 1.0e-6   # |dH/ds| がこれ以下なら 1 次の変化なし (接線) とみなし、有限増分の厳密な変化の符号で判定する (main で上書き)


def outward_k(X, k):
    """条件 k の外向き: 1 次で dH_k/ds < −tol、または接線 (|dH_k/ds| ≤ tol) で有限増分の変化 < 0。(外向き, 接線で 2 次判定した外向き)"""
    d = X[f"dH{k}"]; f = X[f"fH{k}"]
    first = d < -TANGENT_TOL
    second = (np.abs(d) <= TANGENT_TOL) & (f < 0)
    return first | second, second


def outward_stats(rec, X, S, act1, act2, mask):
    """腕 X の増分の向きを集合 S (∩ mask) で評価する。act1/act2: 節点ごとの該当条件 (H1/H2)。"""
    V = rec.V
    Sm = S & mask
    det = Sm & X["okD"] & (act1 | act2)
    w1, s1 = outward_k(X, 1); w2, s2 = outward_k(X, 2)
    o1 = det & act1 & w1; o2 = det & act2 & w2
    outward = o1 | o2
    undet = Sm & ~det
    with np.errstate(all="ignore"):
        dmin = np.where(act1 & act2, np.minimum(X["dH1"], X["dH2"]), np.where(act1, X["dH1"], X["dH2"]))
    n = int(Sm.sum()); no = int(outward.sum()); nu = int(undet.sum()); nd = int(det.sum())
    st = dict(n=n, n_out=no, n_undet=nu, n_det=nd,
              frac=((no + nu) / n if n > 0 else 0.0), frac_det=(no / nd if nd > 0 else np.nan),
              vwmean=(float(np.sum(V[det] * dmin[det]) / np.sum(V[det])) if nd > 0 else np.nan),
              n_a1=int((Sm & act1).sum()), n_a2=int((Sm & act2).sum()), n_o1=int(o1.sum()), n_o2=int(o2.sum()),
              n_tan=int((o1 & s1 | o2 & s2).sum()),
              n_tangent_any=int((det & ((act1 & (np.abs(X["dH1"]) <= TANGENT_TOL)) | (act2 & (np.abs(X["dH2"]) <= TANGENT_TOL)))).sum()),
              n_noS=int((Sm & ~X["okS"]).sum()), n_zero=int((Sm & X["okS"] & ~X["okD"]).sum()),
              n_noact=int((Sm & X["okD"] & ~(act1 | act2)).sum()))
    for k, act in ((1, act1), (2, act2)):
        mk = det & act
        st[f"vw_dH{k}"] = float(np.sum(V[mk] * X[f"dH{k}"][mk]) / np.sum(V[mk])) if mk.any() else np.nan
    ms = Sm & np.isfinite(X["H1_s"]) & X["okS"]
    st["H1s_med"] = float(np.median(X["H1_s"][ms])) if ms.any() else np.nan
    st["H2s_med"] = float(np.median(X["H2_s"][ms])) if ms.any() else np.nan
    return st


def print_outward(label, st, indent="    "):
    print(f"{indent}outward [{label}]: fraction (outward + undetermined)/|S| = {st['frac']:.4f}  (|S| {st['n']}; outward {st['n_out']}, "
          f"undetermined {st['n_undet']} = no droplets at start {st['n_noS']} + zero increment {st['n_zero']} + no identified condition {st['n_noact']})")
    print(f"{indent}  fraction among determined = {st['frac_det']:.4f} ({st['n_det']} nodes); vol-weighted mean of min dH/ds = {st['vwmean']:+.4e}")
    print(f"{indent}  tangent (|dH/ds| <= {TANGENT_TOL:g} for an active condition) {st['n_tangent_any']} nodes; "
          f"outward decided at second order (finite increment) {st['n_tan']} nodes")
    print(f"{indent}  H1 active {st['n_a1']} (outward {st['n_o1']}, vw mean dH1/ds {st['vw_dH1']:+.4e}); "
          f"H2 active {st['n_a2']} (outward {st['n_o2']}, vw mean dH2/ds {st['vw_dH2']:+.4e})")
    print(f"{indent}  start state on S: median H1_s = {st['H1s_med']:+.3e}, median H2_s = {st['H2s_med']:+.3e} (0 = on the boundary)")


def region_metrics(rec, A, mask):
    V = rec.V
    m = {}
    for nm, _c in QCOMP:
        den = float(np.sum(V[mask] * np.abs(A["d_" + nm][mask])))
        mP = mask & A["P"]
        num = float(np.sum(V[mP] * np.abs(A["C_" + nm][mP])))
        denP = float(np.sum(V[mP] * np.abs(A["d_" + nm][mP])))
        mG = mask & A["G0"]
        numG = float(np.sum(V[mG] * np.abs(A["R_" + nm][mG])))
        m[nm] = dict(num=num, den=den, ratio=(num / den if den > 0 else (0.0 if num == 0 else np.inf)),
                     ratioP=(num / denP if denP > 0 else (0.0 if num == 0 else np.inf)),
                     numG=numG, ratioG=(numG / den if den > 0 else (0.0 if numG == 0 else np.inf)))
    m["n_P"] = int(np.count_nonzero(mask & A["P"])); m["n_G0"] = int(np.count_nonzero(mask & A["G0"]))
    m["n_nochange"] = int(np.count_nonzero(mask & A["P_nochange_branch"]))
    m["own"] = outward_stats(rec, A, A["P"], A["v1"], A["v2"], mask)   # 参考: 自身の射影集合・自身の違反条件
    return m


def diag_check(rec):
    """DPLUR 分母の記録: B は液・Q の行が節点ごとに等しく本番の max、A は used == production。"""
    prod = np.stack([slot(rec, f"DPLUR_denom/{r}/production") for r in ("v", "g", "Q2", "Q1", "Q0")])
    used = np.stack([slot(rec, f"DPLUR_denom/{r}/used") for r in ("v", "g", "Q2", "Q1", "Q0")])
    fin = np.all(np.isfinite(prod), axis=0)
    cd = int(rec.attrs.get("common_diag", -1))
    msgs = [f"common_diag attribute = {cd}; DPLUR denominators recorded at {int(fin.sum())} / {rec.n} nodes"]
    ok = True
    if not fin.any():
        return False, msgs + ["  no DPLUR denominator recorded (condTwoPhaseSolver 1 path not taken?)"]
    if np.any(used[0, fin] != prod[0, fin]):
        ok = False; msgs.append("  vapour denominator changed (must not)")
    if cd == 1:
        mx = np.max(prod[1:, fin], axis=0)
        if not np.all(used[1:, fin] == mx[None, :]):
            ok = False; msgs.append("  B: liquid/moment rows are not the per-node max of production")
    elif cd == 0:
        if np.any(used[:, fin] != prod[:, fin]):
            ok = False; msgs.append("  A: used differs from production")
    else:
        ok = False; msgs.append("  common_diag attribute missing")
    with np.errstate(all="ignore"):
        spread = np.max(prod[1:, fin], axis=0) / np.min(prod[1:, fin], axis=0)
    msgs.append(f"  production spread max/min over (g,Q2,Q1,Q0): median {np.nanmedian(spread):.3e}, 90% {np.nanpercentile(spread, 90):.3e}, max {np.nanmax(spread):.3e}")
    return ok, msgs


def print_rec(tag, rec, A, rows):
    print(f"\n=== {tag}: {rec.path} ===")
    for k in ("common_diag", "input_value_file", "cfl_pseudo", "implicitRelax", "nStepInner", "condTwoPhaseSolver", "condTwoPhaseRelax"):
        if k in rec.attrs:
            print(f"  {k} = {rec.attrs[k]}")
    ok, msgs = diag_check(rec)
    for l in msgs:
        print("  " + l)
    for rname, mask in g3.regions(rec):
        if not mask.any():
            continue
        m = region_metrics(rec, A, mask)
        rows[(tag, rname)] = m
        print(f"\n  -- Omega = {rname} ({int(mask.sum())} nodes); projection set P {m['n_P']} nodes "
              f"(branch entered but unchanged {m['n_nochange']}), g=0 dust set G0 {m['n_G0']} nodes --")
        print(f"    {'comp':4s} {'sumV|C_proj| (P)':>18s} {'sumV|delta_lim|':>18s} {'ratio':>10s} {'ratio(den on P)':>16s} {'sumV|C_rm| (G0)':>18s} {'G0 ratio':>10s}")
        for nm, _c in QCOMP:
            q = m[nm]
            print(f"    {nm:4s} {q['num']:18.6e} {q['den']:18.6e} {q['ratio']:10.4f} {q['ratioP']:16.4f} {q['numG']:18.6e} {q['ratioG']:10.4f}")
        print_outward(f"reference, own projection set P_{tag}, own violated conditions", m["own"])
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--length-scale", type=float, default=1000.0)
    ap.add_argument("--ref-q1", type=float, default=0.913, help="A の再現の参照 (run_0580 G3B_RESULT, x35-70 wd<0.4 の Q1 射影比)")
    ap.add_argument("--ref-q2", type=float, default=0.974, help="同 Q2")
    ap.add_argument("--ref-tol", type=float, default=0.05, help="A の再現の許容 (絶対差)")
    ap.add_argument("--tangent-tol", type=float, default=1.0e-6, help="1 次の変化なし (接線) とみなす |dH/ds| の上限 (無次元)")
    a = ap.parse_args()
    global TANGENT_TOL
    TANGENT_TOL = a.tangent_tol
    recA = g3.Rec(a.A, a.length_scale); recB = g3.Rec(a.B, a.length_scale)
    print(f"# A {a.A}: {recA.n} nodes, x [{recA.x.min():.3f}, {recA.x.max():.3f}]; B {a.B}: {recB.n} nodes")
    invalid = []
    if recA.n != recB.n:
        invalid.append("node counts differ")
    if int(recA.attrs.get("common_diag", -1)) != 0:
        invalid.append("A must have common_diag = 0")
    if int(recB.attrs.get("common_diag", -1)) != 1:
        invalid.append("B must have common_diag = 1")
    if not invalid:
        sa = recA.snapshot("upd_start", 0); sb = recB.snapshot("upd_start", 0)
        nd = sum(int(np.count_nonzero(recA.snapshot("upd_start", c) != recB.snapshot("upd_start", c))) for c in range(5)) if sa is not None and sb is not None else -1
        print(f"# update start state (5 components) differs at {nd} values between A and B (expected 0: same restart and assembly)")
        if nd != 0:
            invalid.append("A and B do not start from the same stored state")
        for k in ("input_value_file", "cfl_pseudo", "implicitRelax", "nStepInner", "condTwoPhaseRelax", "condTwoPhaseSolver"):
            if k in recA.attrs and k in recB.attrs and str(recA.attrs[k]) != str(recB.attrs[k]):
                invalid.append(f"attribute {k} differs ({recA.attrs[k]} vs {recB.attrs[k]})")
    AA = analyse(recA); AB = analyse(recB)
    rows = {}
    okA = print_rec("A", recA, AA, rows); okB = print_rec("B", recB, AB, rows)
    if not okA: invalid.append("A denominator record inconsistent")
    if not okB: invalid.append("B denominator record inconsistent")

    # 共通の評価集合 S = P_A、該当条件 = A で違反 ∪ B で違反 (両腕で同じ)
    act1 = AA["v1"] | AB["v1"]; act2 = AA["v2"] | AB["v2"]
    common = {}
    print(f"\n=== outward component on the common set S = P_A (A's projection set; same nodes and start state in both arms) ===")
    for rname, mask in g3.regions(recA):
        if not mask.any():
            continue
        stA = outward_stats(recA, AA, AA["P"], act1, act2, mask); stB = outward_stats(recB, AB, AA["P"], act1, act2, mask)
        common[rname] = (stA, stB)
        print(f"\n  -- Omega = {rname} --")
        print_outward("A increment on S", stA); print_outward("B increment on S", stB)

    print(f"\n=== pre-registered verdict (Omega = {VERDICT_REGION}) ===")
    if (("A", VERDICT_REGION) not in rows) or (("B", VERDICT_REGION) not in rows):
        print("VERDICT: UNDECIDED (verdict region empty: wall_dist missing or x range not in mm; check --length-scale)"); return 2
    mA = rows[("A", VERDICT_REGION)]; mB = rows[("B", VERDICT_REGION)]
    rA1, rA2, rB1, rB2 = mA["Q1"]["ratio"], mA["Q2"]["ratio"], mB["Q1"]["ratio"], mB["Q2"]["ratio"]
    oA, oB = common[VERDICT_REGION][0]["frac"], common[VERDICT_REGION][1]["frac"]
    print(f"  projection ratio Q1: A {rA1:.4f}  B {rB1:.4f}  B/A {rB1/rA1 if rA1 > 0 else float('nan'):.4f}")
    print(f"  projection ratio Q2: A {rA2:.4f}  B {rB2:.4f}  B/A {rB2/rA2 if rA2 > 0 else float('nan'):.4f}")
    print(f"  outward fraction   : A {oA:.4f}  B {oB:.4f}  B/A {oB/oA if oA > 0 else float('nan'):.4f}"
          f"   (on S = P_A; outward = increment decreases a Hankel condition violated at the projection input of A or B; undetermined counted as outward)")
    print(f"  A reproduction: |A_Q1 - {a.ref_q1}| <= {a.ref_tol} and |A_Q2 - {a.ref_q2}| <= {a.ref_tol}")
    print("  SUPPORT: A reproduces, B_Q1 <= 0.1 A_Q1, B_Q2 <= 0.1 A_Q2, and outward resolved := B_outward <= 0.1 A_outward")
    print("  REJECT : A reproduces, B_Q1 >= 0.5 A_Q1, B_Q2 >= 0.5 A_Q2 and B_outward >= 0.5 A_outward;  otherwise UNDECIDED")
    if invalid:
        print("VERDICT: INVALID (" + "; ".join(invalid) + ")"); return 3
    repro = abs(rA1 - a.ref_q1) <= a.ref_tol and abs(rA2 - a.ref_q2) <= a.ref_tol
    if not repro:
        print(f"VERDICT: UNDECIDED (A did not reproduce the recorded projection ratios: Q1 {rA1:.4f}, Q2 {rA2:.4f})"); return 0
    if rB1 <= 0.1 * rA1 and rB2 <= 0.1 * rA2 and oB <= 0.1 * oA:
        print("VERDICT: SUPPORT (hypothesis 1: the different liquid/moment preconditioner denominators drive the update out of the realizable set)")
    elif rB1 >= 0.5 * rA1 and rB2 >= 0.5 * rA2 and oB >= 0.5 * oA:
        print("VERDICT: REJECT (diagonal mismatch is not the main cause; prioritise hypothesis 2)")
    else:
        print("VERDICT: UNDECIDED (intermediate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
