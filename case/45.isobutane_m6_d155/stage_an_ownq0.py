"""(登録外の記録) §6.30 の分解を、各 run 自身の Q0 を起点にして計算する変種。前提の「Q0 の一致」は ρE で満たさない (初期化の ρE の組み直しが精度ごとに違う) ので、判定ではなく記録。元の説明:
plan architecture-float-state-double-geometry §6.30 (事前登録) の判定: 共通の Q32 の状態からの最初の 1 回の更新の段ごとの差。stage.sh の後に B で回す。
FP64 = run_0494〜0496、float = run_0497〜0499 (各 stage.h5)。保存量 ρ・ρu・ρv・ρE (ρw は 2D で全て 0 なら外す)。
  a = q − Q0、bb = b − Q0、e = q − b − d (double)。精度の間の差 Δx = mean_float(x) − mean_FP64(x) で Δa = Δbb + Δd + Δe を照合する。
  領域: 収縮部の内部 (x/r_t ∈ [−5, −1)・j 20〜60)、軸の近く (j 0〜8)、壁の近く (j 117〜120)。L2 は領域の節点で取る。
  再実行の差 = 各精度の中の 3 本の組の ‖a_i − a_j‖ の最大。‖Δa‖ ≤ 10·(再実行の差) の組は「差が小さい」。
  結果 A: ‖Δe‖ ≤ 0.1‖Δa‖ (commit の前で生じる)。結果 B: ‖Δbb + Δd‖ ≤ 0.1‖Δa‖ (commit の丸め)。それ以外は中間。
  総合: 判定対象 (差が小さくない組) がすべて A → A、すべて B → B、それ以外・判定対象なし → 判別不能。
前提: 6 本の Q0 が同じ (double で一致)、記録の経路が lineImplicit、commit の自己照合 (q == fl(b + d)) の不一致 0。
診断の有無の確認: _stage/off_* の res_1 と本番の 1 本目の res_1 の保存量の差が、診断なしどうしの差以下。
結果は _band_ab/cold_pair/stage_judge.json。"""
import json
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
F64 = ("run_0494_st_f64a", "run_0495_st_f64b", "run_0496_st_f64c"); F32 = ("run_0497_st_f32a", "run_0498_st_f32b", "run_0499_st_f32c")
Q = ("ro", "roUx", "roUy", "roUz", "roe")
NJ, RT = 121, 0.076807
OUT = {"undecidable": [], "pre": {}, "table": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def load(r):
    with h5py.File(HERE / r / "stage.h5", "r") as h:
        d = {g: {q: np.asarray(h[f"/{g}/{q}"][:], np.float64) for q in Q} for g in ("Q0", "Q_asm", "R", "b", "d", "q", "q_final")}
        def sc(v):   # 属性: 値が 1 つなら数値に、配列ならリストに
            if isinstance(v, (bytes, str)): return v.decode() if isinstance(v, bytes) else v
            a = np.asarray(v); return a.item() if a.size == 1 else a.tolist()
        d["attrs"] = {k: sc(h.attrs[k]) for k in h.attrs}
        d["cc"] = np.asarray(h["/cells/cc64"][:]).reshape(-1, 3) if "/cells/cc64" in h else None
    return d
S = {}
for r in F64 + F32:
    if not (HERE / r / "stage.h5").exists(): und(f"{r}: stage.h5 が無い"); continue
    S[r] = load(r)
if not OUT["undecidable"]:
    for r, d in S.items():
        at = d["attrs"]; OUT["pre"][r] = {k: at.get(k) for k in ("path", "flow_float_bytes", "n_commit_mismatch_q_vs_fl_b_plus_d", "FORGE_FREEZE_TURB", "step")}
        if at.get("path") != "lineImplicit": und(f"{r}: 経路が lineImplicit でない ({at.get('path')})")
        if int(at.get("n_commit_mismatch_q_vs_fl_b_plus_d", -1)) != 0: und(f"{r}: commit の自己照合の不一致 {at.get('n_commit_mismatch_q_vs_fl_b_plus_d')}")
        if int(at.get("FORGE_FREEZE_TURB", 0)) != 1: und(f"{r}: SST を止めていない")
    r0 = F64[0]
    for r in S:
        pass   # 変種: Q0 の一致は求めない (各 run の Q0 を起点にする)
    print("前提:", json.dumps(OUT["pre"], ensure_ascii=False))
if not OUT["undecidable"]:
    cc = S[F64[0]]["cc"]; n = len(S[F64[0]]["Q0"]["ro"]); jj = np.arange(n) % NJ; xr = cc[:n, 0] / RT
    REG = {"収縮部の内部": (xr >= -5) & (xr < -1) & (jj >= 20) & (jj <= 60), "軸の近く": jj <= 8, "壁の近く": jj >= 117}
    qs = [q for q in Q if not (q == "roUz" and all(np.all(S[r]["q"][q] == 0) for r in S))]
    Q0 = S[F64[0]]["Q0"]
    def parts(r, q):
        d = S[r]; return {"a": d["q"][q] - d["Q0"][q], "bb": d["b"][q] - d["Q0"][q], "d": d["d"][q], "e": d["q"][q] - d["b"][q] - d["d"][q],
                          "R": d["R"][q], "asm": d["Q_asm"][q] - d["b"][q]}
    P = {r: {q: parts(r, q) for q in qs} for r in S}
    l2 = lambda v, m: float(np.sqrt(np.mean(v[m] ** 2)))
    cls_all = []
    print("領域        量      ‖Δa‖       ‖Δe‖/‖Δa‖  ‖Δbb+Δd‖/‖Δa‖  再実行の差  10×   照合の残り   ‖ΔR‖/‖R‖   分類   (符号つきの和 Δa・Δd・Δe)")
    for rn, m in REG.items():
        for q in qs:
            mean = lambda runs, k: np.mean([P[r][q][k] for r in runs], axis=0)
            D = {k: mean(F32, k) - mean(F64, k) for k in ("a", "bb", "d", "e", "R")}
            na, ne, nbd = l2(D["a"], m), l2(D["e"], m), l2(D["bb"] + D["d"], m)
            ident = l2(D["a"] - (D["bb"] + D["d"] + D["e"]), m)
            noise = max(l2(P[x][q]["a"] - P[y][q]["a"], m) for runs in (F64, F32) for i, x in enumerate(runs) for y in runs[i + 1:])
            nR = l2(mean(F64, "R"), m)
            if na <= 10 * noise or na == 0: c = "小さい"
            elif ne <= 0.1 * na: c = "A"
            elif nbd <= 0.1 * na: c = "B"
            else: c = "中間"
            if c != "小さい": cls_all.append(c)
            OUT["table"][f"{rn}|{q}"] = dict(norm_da=na, ratio_e=ne / na if na else None, ratio_bd=nbd / na if na else None, noise=noise, identity_residual=ident,
                                             dR_over_R=l2(D["R"], m) / nR if nR else None, cls=c,
                                             sum_da=float(D["a"][m].sum()), sum_dd=float(D["d"][m].sum()), sum_de=float(D["e"][m].sum()), sum_dbb=float(D["bb"][m].sum()))
            t = OUT["table"][f"{rn}|{q}"]
            print(f"{rn:10s}  {q:6s}  {na:.3e}  {t['ratio_e'] if t['ratio_e'] is not None else float('nan'):9.3e}  {t['ratio_bd'] if t['ratio_bd'] is not None else float('nan'):12.3e}  {noise:.2e}  {10 * noise:.2e}  {ident:.1e}  {t['dR_over_R'] if t['dR_over_R'] is not None else float('nan'):.2e}  {c:4s}  ({t['sum_da']:+.2e}, {t['sum_dd']:+.2e}, {t['sum_de']:+.2e})")
    if not cls_all: verdict = "判別不能 (どの組も差が小さい)"
    elif all(c == "A" for c in cls_all): verdict = "結果 A: 最初の更新の差は主に commit の前で生じる (第 1 仮説を支持、commit の丸めを主因とする説を退ける)"
    elif all(c == "B" for c in cls_all): verdict = "結果 B: 最初の更新の差は主に commit の丸めで生じる"
    else: verdict = f"判別不能 (領域・量で分類が一致しない: {sorted(set(cls_all))})"
else:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== (記録・登録外、各 run の Q0 起点) 分類: {verdict}")
# 診断の有無の確認
chk = {}
for prec, on, offs in (("FP64", F64[0], ("_stage/off_f64_1", "_stage/off_f64_2")), ("float", F32[0], ("_stage/off_f32_1", "_stage/off_f32_2"))):
    try:
        rd = lambda r: {q: np.asarray(h5py.File(HERE / r / "res_1.h5", "r")["VALUE/" + q][:], np.float64) for q in Q}
        a, b, c = rd(on), rd(offs[0]), rd(offs[1])
        d_on = max(float(np.max(np.abs(a[q] - b[q]))) for q in Q); d_off = max(float(np.max(np.abs(b[q] - c[q]))) for q in Q)
        chk[prec] = dict(on_vs_off=d_on, off_vs_off=d_off, ok=bool(d_on <= d_off or d_on == 0))
        print(f"  診断の有無 {prec}: 診断あり−なし の最大 {d_on:.3e}、なしどうし {d_off:.3e} → {'OK' if chk[prec]['ok'] else '要確認'}")
    except Exception as ex:
        chk[prec] = dict(error=str(ex)); print(f"  診断の有無 {prec}: 確認できない ({ex})")
OUT["diag_on_off"] = chk
(HERE / "_band_ab/cold_pair/stage_ownq0_record.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/stage_ownq0_record.json")
