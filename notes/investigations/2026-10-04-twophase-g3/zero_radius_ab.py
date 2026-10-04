"""#4g3c: 半径 0 の重みを蒸発の Q1 ソースから除く反実仮想 (0 step、新規 run なし)。

codex diagnose (notes/reviews/2026-10-04-twophase-g3a-result-diagnose.md) の判別 A/B。本番修正案ではない。
  A: 記録された現行ソース (tp_operator.h5 の /source/slots term_*)。
  B: 蒸発分岐の節点で S_Q1 だけを S_Q1·q1²/(q0 q2) に置換 (同じ ṙ のまま q0 を正半径側の重み w1 = q1²/q2 に替える)。
     S_Q0・S_Q2・S_g・輸送項は固定。
余裕 H_up = ln q1 + ln q3 − 2 ln q2 (q3 = ρg/(4/3 π ρ_l)) の物性固定の方向微分 Ḣ = Σ a_c R_c/(V q_c)
(温度変化による −d ln ρ_l/dt は含めない)。外向き寄与 N = Σ V q1 · max(−Ḣ, 0)。
対象集合 T = 射影集合 (run_0592、extract_projection_mask.py) ∩ 全成分正 ∩ 蒸発分岐 ∩ abs(D)/(q1 q3 + q2²) ≤ 1e-6
  (D = q1 q3 − q2²、ρ_l = h2o_rho_cond(T) 倍精度の式。本番 float 経路は物性表なので表と式の差は記録のみ)。
事前登録 (plan #4g3c): 覆い N_A(T)/N_A(射影集合の全成分正) ≥ 0.9 かつ 3 本とも N_B/N_A ≤ 0.1 → 支持、
  N_B/N_A ≥ 0.5 → 棄却、それ以外・覆い不足・記録不整合 → 判別不能。ソース単独と全残差の両方を出す。
使い方: python3 zero_radius_ab.py <proj.npy> <ON tp_operator.h5> [...]
"""
import sys
import h5py
import numpy as np

PI = np.pi


def rho_l_h2o(T):
    r = 1000.0 - 0.12 * (277.0 - T)
    return np.maximum(r, 920.0)


def names(ds):
    return [x.decode() if isinstance(x, bytes) else x for x in ds[:]]


def evaluate(path, mask):
    f = h5py.File(path, "r")
    comps = names(f.attrs["components"]) if not isinstance(f.attrs["components"][0], str) else list(f.attrs["components"])
    ci = {c: comps.index(c) for c in ["g", "Q2", "Q1", "Q0"]}
    st = f["node/state"][:].astype(np.float64)
    V = f["node/volume"][:].astype(np.float64)
    T = f["node/T"][:].astype(np.float64)
    labels = names(f["res/labels"])
    snap = f["res/snap"][:].astype(np.float64)
    slots = names(f["source/slot_names"])
    sl = f["source/slots"][:].astype(np.float64)
    idx = np.asarray(mask, dtype=np.int64)
    q = {c: st[ci[c], idx] for c in ci}
    pos = np.all([q[c] > 0 for c in ci], axis=0)
    idx = idx[pos]
    q = {c: st[ci[c], idx] for c in ci}
    v = V[idx]
    rl = rho_l_h2o(T[idx])
    q3 = q["g"] / (4.0 / 3.0 * PI * rl)
    D = q["Q1"] * q3 - q["Q2"] ** 2
    Drel = np.abs(D) / (q["Q1"] * q3 + q["Q2"] ** 2)
    branch = sl[slots.index("branch"), idx]
    evap = (branch == 2) | (branch == 12)
    r30 = np.cbrt(q3 / q["Q0"])
    cap1 = q["Q1"] > q["Q0"] * r30
    cap2 = q["Q2"] > q["Q0"] * r30 * r30
    term = {c: sl[slots.index("term_" + c), idx] for c in ci}
    # 記録の整合: ソース段の残差差分 (cond_src − tp_diff) と term_* の一致
    L = {l: labels.index(l) for l in labels}
    srcR = {c: snap[L["cond_src"], ci[c], idx] - snap[L["tp_diff"], ci[c], idx] for c in ci}
    base = {c: np.abs(snap[L["tp_diff"], ci[c], idx]) + np.abs(term[c]) + 1e-300 for c in ci}
    mism = max(float(np.max(np.abs(srcR[c] - term[c]) / base[c])) for c in ci)  # float の加算丸め (~6e-8) 程度なら整合
    # 改訂 (2026-10-04、結果を見た後): 不一致はすべて float の非正規化数の範囲 (abs < FLT_MIN) の丸め (刻み 2^-149)。
    # 丸め回数 nround ごとに絶対 2^-149 の床を足した整合 (G1 #4sr と同じ型)。旧基準の値も併記する。
    nround = sl[slots.index("nround"), idx]
    floor = np.maximum(nround, 1.0) * 2.0 ** -149
    mism_rev = max(float(np.max(np.maximum(np.abs(srcR[c] - term[c]) - floor, 0.0) / base[c])) for c in ci)
    finR = {c: snap[L["final"], ci[c], idx] for c in ci}

    def hdot(R):
        return R["Q1"] / (v * q["Q1"]) + R["g"] / (v * q["g"]) - 2.0 * R["Q2"] / (v * q["Q2"])

    termB = dict(term)
    termB["Q1"] = np.where(evap, term["Q1"] * q["Q1"] ** 2 / (q["Q0"] * q["Q2"]), term["Q1"])
    finB = dict(finR)
    finB["Q1"] = finR["Q1"] - term["Q1"] + termB["Q1"]
    H = {"src A": hdot(term), "src B": hdot(termB), "final A": hdot(finR), "final B": hdot(finB)}
    w = v * q["Q1"]
    N = {k: w * np.maximum(-h, 0.0) for k, h in H.items()}
    tgt = evap & (Drel <= 1e-6)
    out = {"n_pos": int(idx.size), "n_evap": int(evap.sum()), "n_tgt": int(tgt.sum()),
           "n_cap1_tgt": int((cap1 & tgt).sum()), "n_cap2_tgt": int((cap2 & tgt).sum()),
           "w0_frac_med": float(np.median(1.0 - (q["Q1"] ** 2 / q["Q2"])[tgt] / q["Q0"][tgt])) if tgt.any() else np.nan,
           "rec_mismatch_max_rel": mism, "rec_mismatch_max_rel_rev": mism_rev}
    for k in H:
        out[f"N {k} (T)"] = float(N[k][tgt].sum())
        out[f"N {k} (all pos)"] = float(N[k].sum())
        out[f"N {k} (outside T)"] = float(N[k][~tgt].sum())
        out[f"frac<0 {k} (T)"] = float(np.mean(H[k][tgt] < 0)) if tgt.any() else np.nan
        out[f"median Hdot {k} (T)"] = float(np.median(H[k][tgt])) if tgt.any() else np.nan
    out["coverage"] = out["N src A (T)"] / out["N src A (all pos)"]
    out["ratio src"] = out["N src B (T)"] / out["N src A (T)"]
    out["ratio final"] = out["N final B (T)"] / out["N final A (T)"]
    return out


REVISED = False


def main():
    global REVISED
    args = [a for a in sys.argv[1:] if a != "--revised-floor"]
    REVISED = len(args) != len(sys.argv) - 1
    print("記録の整合: " + ("改訂 (非正規化数の絶対床 nround·2^-149 を許す; 結果を見た後の改訂)" if REVISED else "事前登録 (相対 1e-5)"))
    sys.argv = [sys.argv[0]] + args
    mask = np.load(sys.argv[1])
    verdicts = []
    for p in sys.argv[2:]:
        o = evaluate(p, mask)
        print(f"\n## {p}")
        for k, val in o.items():
            print(f"  {k:28s} {val:.6g}" if isinstance(val, float) else f"  {k:28s} {val}")
        key = "rec_mismatch_max_rel_rev" if REVISED else "rec_mismatch_max_rel"
        if o[key] > 1e-5 or o["coverage"] < 0.9:
            vd = "判別不能"
        elif o["ratio src"] <= 0.1:
            vd = "支持"
        elif o["ratio src"] >= 0.5:
            vd = "棄却"
        else:
            vd = "判別不能"
        print(f"  -> {vd} (coverage {o['coverage']:.3f}, N_B/N_A src {o['ratio src']:.3g}, final {o['ratio final']:.3g})")
        verdicts.append(vd)
    allv = set(verdicts)
    v = verdicts[0] if len(allv) == 1 else "判別不能"
    print(f"\nVERDICT #4g3c: {v}  ({', '.join(verdicts)}); 全残差の外向きの残りは 'ratio final' で別記")


if __name__ == "__main__":
    main()
