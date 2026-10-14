#!/usr/bin/env python3
r"""case/71 run_0002_ns の全圧の超過の A/B (forge 0 step)。plan tooling-nozzle-card-local-mesh §5.1 #3c、
codex 諮問 `notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md` の「判別 A/B」。

同じ保存場 (既定 res_100000〜180000 の 1 万 step ごと + 181,000〜183,000、12 枚)・同じ節点・同じ物性の記録
(res の属性が指す `resolved_species_*.yaml`; `forge_species.run_thermo` → `thermo_gas`)・同じ datum (thermoHrefTemp) で、
変えるのは全温・全圧の逆算に渡す状態の取得元だけ:

  - A: 今の `total_quantities.total_state` (出力の原始量 T・P・Y と `VALUE/h0` から) の T0・P0。
  - B: 保存量 ro・roUx・roUy・roUz・roe・roY* から、ソルバ (`cuda_forge/dependentVariables_d.cu`、thermalMethod 2・
       thermoFloat 1・凝縮 off の経路) と同じ規約で状態を作り直す:
         ρ = max(ro, roMin)、U = ρU/ρ、e = roe/ρ − |U|²/2 (h0_includes_k 0 なので roe に k は入っていない)、
         組成は float32 で Y_s = max(ρY_s·(float)(1/ρ), 0) → Σ で正規化 (ハイブリッド経路の非負化・正規化と同じ演算)、
         T は e_mix(T) = h_mix(T) − R_mix T = e の Newton (thermo_gas、範囲 [50, 6000] K)、P = ρ R_mix T (下限 pMin)、
         h0 = e_mix(T) + P/ρ + |U|²/2 (= ソルバの Ht; `output/output.cpp` の h0 は includes_k 0 のとき Ht そのもの)。
       全温・全圧は A と同じ式 (h_mix(T0) = h0 の Newton、P0 = P exp((s°(T0) − s°(T))/R_mix)、凍結組成)。定数 cp の近似は使わない。
  - 補助 B_dbl (判定に使わない): B の組成の正規化を double で行ったもの (thermoFloat の float の組成の丸めの寄与を分ける)。

ソルバの手順の上の位置 (`main.cpp` advanceImplicitSteady): 出力 (writeStepOutputs) は陰解法の更新と等温壁のピンの**後**、
次の step の dependentVariables の**前**に書かれる。したがって出力の原始量・Ht はその step の頭 (更新前、境界条件の上書きの後) の状態、
保存量は更新後の状態である (A − B はこの 1 step のずれ + 境界の上書き + 組成の丸めを含む)。

節点ごとに: ΔP0/Pt (A・B・A−B)、超過 (P0 > Pt (1 + 1e−6)) の点数・分位点・位置 (入口の列 i = 0・出口の列 i = ni−1・壁の行
j = nj−1・軸の行 j = 0・内部の別; index = i·nj + j)、T0 > Tt (1 + 1e−6) と T0 > Tt + 1 K の点、エネルギーの逆算の残差
(B の T で e を作り直したときの e_mix(T_B) − e、および出力の T を保存量の組成・e に当てたときの e_mix(T_A) − e)。

判定 (呼び出し側が事前に固定、2026-10-15; 条件・閾値は変えない):
  - 両経路の差 |P0_A − P0_B| が全時点で ≤ 1e−6 Pt (全節点) かつ B で超過が残る → 「出力の状態の不整合が +0.24 % の主因」を棄却
    (保存場の側に偏差)
  - B で超過が消え (全時点で B の超過 0)、A/B の差が超過を説明 (A の超過の点で B ≤ Pt (1 + 1e−6)) → 「保存場の側」を棄却し、
    出力・後処理の経路を直す対象として記録
  - それ以外 (一部だけ) → 複合要因として記録 (過渡の原因の確定には進まない)
補助 (判定に使わない): 同じ差の最大を A の超過の点の集合・内部の節点だけで取ったもの。

時系列 (`--series-csv`、`check_quasisteady.py --series-csv` 用): 最後の時点 (既定 183,000) の A の超過の点の集合 S を固定し、
S 全体・S の軸の近く (j ≤ 2、入口・出口の列を除く)・S の内部 (3 ≤ j ≤ nj − 2、同) の平均の ΔP0/Pt と h0 − h0_in
(h0_in = h_mix(Tt, 入口の Y))、各部分の最大の超過の節点 (代表点) の同じ量、超過の点数 (A・B)。

使い方 (AWS A):
  python3 p0_ab.py --tools ~/forge-cardmesh/solver_density_cuda/tools --run <run_0002_ns> \
      --out-json <出力先>/p0_ab.json --series-csv <出力先>/p0_ab_series.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import yaml

THR = 1.0e-6                     # 超過・差の閾値 (Pt に対する相対; 呼び出し側が事前に固定)
TMIN, TMAX = 50.0, 6000.0        # dependentVariables_d.cu の DEPVAR_TMIN / DEPVAR_TMAX
STEPS_DEFAULT = [100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 181000, 182000, 183000]
QS = (0, 1, 10, 50, 90, 99, 99.9, 100)
RULE = ("両経路の差 |P0_A − P0_B| が全時点で ≤ 1e−6 Pt かつ B で超過が残る → 「出力の状態の不整合が +0.24 % の主因」を棄却 (保存場の側に偏差); "
        "B で超過が消え A/B の差が超過を説明 → 「保存場の側」を棄却し、出力・後処理の経路を直す対象として記録; 一部だけ → 複合要因として記録 "
        "(呼び出し側が事前に固定、2026-10-15)")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def quant(v) -> dict:
    v = np.asarray(v, float)
    if v.size == 0:
        return {"n": 0}
    return {"n": int(v.size), **{str(q): float(np.percentile(v, q)) for q in QS}}


class Where:
    """節点の位置の分類 (node の構造格子: index = i·nj + j、j = 0 が軸・j = nj − 1 が壁)。
    優先順: 入口の列 (i = 0) → 出口の列 (i = ni − 1) → 壁の行 (j = nj − 1) → 軸の行 (j = 0) → 内部。"""

    def __init__(self, ni, nj, x):
        self.ni, self.nj, self.x = ni, nj, x
        idx = np.arange(ni * nj)
        self.i, self.j = idx // nj, idx % nj
        cat = np.full(ni * nj, "interior", dtype=object)
        cat[self.j == 0] = "axis_row"
        cat[self.j == nj - 1] = "wall_row"
        cat[self.i == ni - 1] = "outlet_col"
        cat[self.i == 0] = "inlet_col"
        self.cat = cat
        self.names = ("inlet_col", "outlet_col", "wall_row", "axis_row", "interior")

    def of(self, mask) -> dict:
        idx = np.nonzero(mask)[0]
        out = {"n": int(idx.size)}
        if idx.size == 0:
            return out
        out["by_cat"] = {c: int(np.count_nonzero(self.cat[idx] == c)) for c in self.names}
        jw = self.nj - 1 - self.j[idx]
        out["x_m"] = [float(self.x[idx].min()), float(self.x[idx].max())]
        out["i"] = [int(self.i[idx].min()), int(self.i[idx].max())]
        out["j_from_wall"] = [int(jw.min()), int(jw.max())]
        xb = np.floor(self.x[idx]).astype(int)
        out["by_x_1m"] = {f"[{b},{b + 1})": int(c) for b, c in zip(*np.unique(xb, return_counts=True))}
        return out


def newton_e(gas, Y, R, e, T_start):
    """e_mix(T) = h_mix(T) − R_mix T = e (thermo_T_from_e と同じ f・f' = cv、毎反復 [TMIN, TMAX] にクランプ)。厳密な根まで回す。"""
    T = np.clip(np.where(T_start > TMIN, T_start, 300.0), TMIN, TMAX)
    it = 0
    for it in range(1, 61):
        cv = gas.cp(Y, T) - R
        cvf = np.maximum(cv, 1.0e-2 * R)
        dT = (gas.h(Y, T) - R * T - e) / cvf
        dT = np.clip(dT, -0.5 * T, 0.5 * T)
        T = np.clip(T - dT, TMIN, TMAX)
        if np.max(np.abs(dT)) < 1.0e-10 * TMAX:
            break
    return T, it, float(np.max(np.abs(dT)))


def newton_h(gas, Y, h0, T_start):
    """h_mix(T0) = h0 (total_state と同じ式; 収束はより厳しく)。"""
    T0 = np.asarray(T_start, float).copy()
    it = 0
    for it in range(1, 61):
        r = gas.h(Y, T0) - h0
        dT = -r / np.maximum(gas.cp(Y, T0), 1.0)
        T0 = np.maximum(T0 + dT, 1.0)
        if np.max(np.abs(dT)) < 1.0e-10 * TMAX:
            break
    return T0, it, float(np.max(np.abs(dT)))


def compose_float(ro_t, roY):
    """dependentVariables_d.cu の useHybrid 経路の組成: float inv_ro = 1.0f/ro (double 除算 → float)、
    float y = roY·inv_ro (double 積 → float)、y < 0 → 0、float の和、float inv = 1.0f/ysum、Yf *= inv (float)、Y = (double)Yf。"""
    inv_ro = (1.0 / ro_t).astype(np.float32)
    ys = []
    nneg = 0
    for r in roY:
        y = (r * inv_ro.astype(np.float64)).astype(np.float32)
        nneg += int(np.count_nonzero(y < 0))
        ys.append(np.where(y < np.float32(0), np.float32(0), y))
    ysum = np.zeros_like(ys[0])
    for y in ys:
        ysum = (ysum + y).astype(np.float32)
    inv = (np.float32(1.0) / np.maximum(ysum, np.float32(1.0e-30))).astype(np.float32)
    return [(y * inv).astype(np.float32).astype(np.float64) for y in ys], nneg


def compose_double(ro_t, roY):
    """非ハイブリッド経路の組成 (double; 補助 B_dbl 用)。"""
    ys = [np.maximum(r / ro_t, 0.0) for r in roY]
    s = np.maximum(sum(ys), 1.0e-30)
    return [y / s for y in ys]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--steps", default=",".join(str(s) for s in STEPS_DEFAULT))
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--series-csv", required=True)
    a = ap.parse_args()
    t_start = time.time()
    sys.path.insert(0, str(Path(a.tools).resolve()))
    import forge_species as fsp
    from total_quantities import total_state

    run = Path(a.run).resolve()
    steps = [int(s) for s in a.steps.split(",")]
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    pp = cfg["physProp"]
    assert int(pp.get("thermalMethod", 0)) == 2, "thermalMethod 2 (TP) 専用"
    roMin = float(pp.get("roMin", 1.0e-4))
    pMin = float(pp.get("pMin", 1.0))
    thermoFloat = int(pp.get("thermoFloat", 1))
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    inl = bc["inlet"]["floats"]
    Pt, Tt = float(inl["Pt"]), float(inl["Tt"])
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.array(h["/MESH/COORD"][:]).reshape(-1, 3)
    W = Where(ni, nj, xy[:, 0])

    per = {}
    keep = {}
    th_ref = None
    for st in steps:
        t0 = time.time()
        res = run / f"res_{st}.h5"
        with h5py.File(res, "r") as h:
            V = {k: np.array(h["VALUE"][k][:], dtype=np.float64) for k in
                 ("ro", "roUx", "roUy", "roUz", "roe", "P", "T", "h0", "Y0", "Y1", "Ux", "Uy", "Uz")}
            nsp = len([k for k in h["VALUE"].keys() if k.startswith("roY")])
            roY = [np.array(h["VALUE"][f"roY{s}"][:], dtype=np.float64) for s in range(nsp)]
            inc_k = int(h["VALUE/h0"].attrs.get("h0_includes_k", 0))
            if any(k.startswith("g_") or k.startswith("rog") for k in h["VALUE"].keys()):
                raise SystemExit(f"{res}: 凝縮の量がある (この台本は凝縮 off の経路だけを写している)")
            rec_attr = h.attrs.get("species_record_file")
            rec_attr = rec_attr.decode() if isinstance(rec_attr, bytes) else str(rec_attr)
        if inc_k != 0:
            raise SystemExit(f"{res}: h0_includes_k={inc_k} (この台本は k を含まない h0 だけを扱う)")
        th = fsp.run_thermo(str(run), res_path=str(res), source="auto")
        if th_ref is None:
            th_ref = th
        elif th["path"] != th_ref["path"]:
            raise SystemExit(f"物性の記録が時点で違う: {th['path']} vs {th_ref['path']}")
        gas = fsp.thermo_gas(th)

        # ---- A: total_state (出力の原始量・h0) ----
        sA = total_state(str(run), str(res))
        T0A, P0A = np.asarray(sA["T0"], float), np.asarray(sA["P0"], float)
        TA, PA = V["T"], V["P"]
        YA = [V["Y0"], V["Y1"]][:nsp]

        # ---- B: 保存量からソルバの規約で ----
        ro_t = np.maximum(V["ro"], roMin)
        n_ro_clamp = int(np.count_nonzero(V["ro"] < roMin))
        Ux, Uy, Uz = V["roUx"] / ro_t, V["roUy"] / ro_t, V["roUz"] / ro_t
        ek = 0.5 * (Ux * Ux + Uy * Uy + Uz * Uz)
        e_in = V["roe"] / ro_t - ek
        out_st = {"res": res.name, "sha256": sha256(res), "species_record_attr": rec_attr}
        variants = {}
        for name, (Y, nneg) in (("B", compose_float(ro_t, roY) if thermoFloat else (compose_double(ro_t, roY), 0)),
                                ("B_dbl", (compose_double(ro_t, roY), 0))):
            R = gas.Rmix(Y)
            T, itT, dTlast = newton_e(gas, Y, R, e_in, TA)
            e_rec = gas.h(Y, T) - R * T
            P = np.maximum(ro_t * R * T, pMin)
            h0 = e_rec + P / ro_t + ek
            T0, it0, dT0last = newton_h(gas, Y, h0, T0A)
            P0 = P * np.exp((gas.s0(Y, T0) - gas.s0(Y, T)) / R)
            variants[name] = dict(Y=Y, R=R, T=T, P=P, h0=h0, T0=T0, P0=P0, e_rec=e_rec, nneg=nneg, itT=itT, dTlast=dTlast,
                                  it0=it0, dT0last=dT0last, n_T_clamp=int(np.count_nonzero((T <= TMIN) | (T >= TMAX))),
                                  n_P_floor=int(np.count_nonzero(ro_t * R * T < pMin)))
        B = variants["B"]
        Bd = variants["B_dbl"]
        # エネルギーの逆算の残差
        resB = B["e_rec"] - e_in
        RB = B["R"]
        resA_on_cons = (gas.h(B["Y"], TA) - RB * TA) - e_in      # 出力の T を保存量の組成・e に当てたときの差
        # 比較
        exA = P0A > Pt * (1 + THR)
        exB = B["P0"] > Pt * (1 + THR)
        exBd = Bd["P0"] > Pt * (1 + THR)
        dAB = (P0A - B["P0"]) / Pt
        dABd = (Bd["P0"] - B["P0"]) / Pt
        exA_amt = P0A - Pt
        expl_num = np.minimum(np.maximum(P0A - B["P0"], 0.0), exA_amt)[exA].sum()
        expl_den = exA_amt[exA].sum()
        interior = W.cat == "interior"
        out_st.update({
            "species_source": sA.get("species_source"), "method_A": sA["method"],
            "B_newton": {"T_iter": B["itT"], "T_last_dT_max": B["dTlast"], "T0_iter": B["it0"], "T0_last_dT_max": B["dT0last"],
                         "n_T_clamped": B["n_T_clamp"], "n_P_floor": B["n_P_floor"], "n_ro_floor": n_ro_clamp,
                         "n_Y_negative_clamped": B["nneg"]},
            "P0": {
                "A_over_Pt_m1": quant(P0A / Pt - 1), "B_over_Pt_m1": quant(B["P0"] / Pt - 1),
                "A_minus_B_over_Pt": quant(dAB), "abs_A_minus_B_over_Pt": quant(np.abs(dAB)),
                "abs_A_minus_B_over_Pt_max_on_A_exceed": float(np.abs(dAB[exA]).max()) if exA.any() else None,
                "abs_A_minus_B_over_Pt_max_interior": float(np.abs(dAB[interior]).max()),
                "abs_A_minus_B_over_Pt_max_by_cat": {c: float(np.abs(dAB[W.cat == c]).max()) for c in W.names},
                "n_abs_A_minus_B_gt_thr": W.of(np.abs(dAB) > THR),
                "exceed_A": W.of(exA), "exceed_B": W.of(exB),
                "exceed_A_and_B": int(np.count_nonzero(exA & exB)), "exceed_A_only": int(np.count_nonzero(exA & ~exB)),
                "exceed_B_only": int(np.count_nonzero(~exA & exB)),
                "on_A_exceed": {"A_over_Pt_m1": quant((P0A / Pt - 1)[exA]), "B_over_Pt_m1": quant((B["P0"] / Pt - 1)[exA]),
                                "A_minus_B_over_Pt": quant(dAB[exA])},
                "explained_fraction": float(expl_num / expl_den) if expl_den > 0 else None,
                "B_dbl": {"exceed": int(np.count_nonzero(exBd)), "abs_Bdbl_minus_B_over_Pt": quant(np.abs(dABd)),
                          "over_Pt_m1_max": float((Bd["P0"] / Pt - 1).max())},
            },
            "T0": {
                "A_minus_Tt": quant(T0A - Tt), "B_minus_Tt": quant(B["T0"] - Tt), "A_minus_B": quant(T0A - B["T0"]),
                "exceed_rel_A": W.of(T0A > Tt * (1 + THR)), "exceed_rel_B": W.of(B["T0"] > Tt * (1 + THR)),
                "exceed_1K_A": W.of(T0A > Tt + 1.0), "exceed_1K_B": W.of(B["T0"] > Tt + 1.0),
            },
            "static": {
                "P_A_minus_B_rel": quant((PA - B["P"]) / B["P"]), "T_A_minus_B": quant(TA - B["T"]),
                "h0_A_minus_B": quant(V["h0"] - B["h0"]), "Y0_A_minus_B": quant(YA[0] - B["Y"][0]),
                "T_B_minus_Tt": quant(B["T"] - Tt), "n_T_B_gt_Tt": W.of(B["T"] > Tt),
                "n_P_B_gt_P0_B": int(np.count_nonzero(B["P"] > B["P0"])),
                "max_P_B_over_P0_B": float((B["P"] / B["P0"]).max()),
            },
            "energy_residual": {
                "B_e_rec_minus_e": quant(resB), "B_abs_rel": quant(np.abs(resB) / np.maximum(np.abs(e_in), 1.0)),
                "outputT_on_conserved_e_minus_e": quant(resA_on_cons),
                "outputT_on_conserved_abs_rel": quant(np.abs(resA_on_cons) / np.maximum(np.abs(e_in), 1.0)),
            },
            "wall_s": time.time() - t0,
        })
        per[st] = out_st
        keep[st] = dict(P0A=P0A, P0B=B["P0"], h0A=V["h0"], h0B=B["h0"], YB=B["Y"], gas=gas,
                        dmax_all=float(np.abs(dAB).max()), nexA=int(exA.sum()), nexB=int(exB.sum()),
                        exB_on_A=int(np.count_nonzero(exA & exB)))
        print(f"[p0_ab] {res.name}: |A−B|/Pt max {np.abs(dAB).max():.3e} (interior {np.abs(dAB[interior]).max():.3e}), "
              f"exceed A {int(exA.sum())} B {int(exB.sum())} (A∩B {int(np.count_nonzero(exA & exB))}), "
              f"P0/Pt−1 max A {(P0A / Pt - 1).max():.6e} B {(B['P0'] / Pt - 1).max():.6e}, "
              f"|e_rec−e| max {np.abs(resB).max():.3e} J/kg, {time.time() - t0:.0f} s", flush=True)

    # ---- 判定 (事前に固定) ----
    dmax = {st: keep[st]["dmax_all"] for st in steps}
    cond_small = all(v <= THR for v in dmax.values())
    remains = all(keep[st]["nexB"] > 0 for st in steps)
    vanish = all(keep[st]["nexB"] == 0 for st in steps)
    explains = all(keep[st]["exB_on_A"] == 0 for st in steps)
    if cond_small and remains:
        verdict = "REJECT_OUTPUT_INCONSISTENCY (出力の状態の不整合を +0.24 % の主因とする仮説を棄却: 保存場の側に偏差)"
    elif vanish and explains:
        verdict = "REJECT_CONSERVED_SIDE (保存場の側という仮説を棄却: 出力・後処理の経路を直す対象として記録)"
    else:
        verdict = "COMPOUND (一部だけ: 複合要因として記録)"
    # 補助 (判定に使わない): 差の最大を A の超過の点・内部だけで
    aux = {st: {"on_A_exceed": per[st]["P0"]["abs_A_minus_B_over_Pt_max_on_A_exceed"],
                "interior": per[st]["P0"]["abs_A_minus_B_over_Pt_max_interior"]} for st in steps}

    # ---- 時系列 (最後の時点の A の超過の集合 S で固定) ----
    last = steps[-1]
    S = keep[last]["P0A"] > Pt * (1 + THR)
    cols_ok = (W.cat != "inlet_col") & (W.cat != "outlet_col")
    S_ax = S & cols_ok & (W.j <= 2)
    S_in = S & cols_ok & (W.j >= 3) & (W.j <= nj - 2)
    Yin = [float(inl[f"Y{s}"]) for s in range(len(th_ref["names"]))]
    gas_last = keep[last]["gas"]
    h0_in = float(gas_last.h([np.array([y]) for y in Yin], np.array([Tt]))[0])

    def rep(mask):
        idx = np.nonzero(mask)[0]
        return int(idx[np.argmax(keep[last]["P0A"][idx])]) if idx.size else None
    p_ax, p_in = rep(S_ax), rep(S_in)
    rows = []
    for st in steps:
        k = keep[st]
        dA, dB = k["P0A"] / Pt - 1, k["P0B"] / Pt - 1
        hA, hB = k["h0A"] - h0_in, k["h0B"] - h0_in
        r = {"step": st, "nexc_A": k["nexA"], "nexc_B": k["nexB"],
             "S_dP0_A_mean": dA[S].mean(), "S_dP0_B_mean": dB[S].mean(), "S_dP0_A_max": dA[S].max(), "S_dP0_B_max": dB[S].max(),
             "S_dh0_A_mean": hA[S].mean(), "S_dh0_B_mean": hB[S].mean(),
             "Sax_dP0_A_mean": dA[S_ax].mean(), "Sax_dP0_B_mean": dB[S_ax].mean(),
             "Sax_dh0_A_mean": hA[S_ax].mean(), "Sax_dh0_B_mean": hB[S_ax].mean(),
             "Sin_dP0_A_mean": dA[S_in].mean(), "Sin_dP0_B_mean": dB[S_in].mean(),
             "Sin_dh0_A_mean": hA[S_in].mean(), "Sin_dh0_B_mean": hB[S_in].mean(),
             "pax_dP0_A": dA[p_ax], "pax_dP0_B": dB[p_ax], "pax_dh0_A": hA[p_ax], "pax_dh0_B": hB[p_ax],
             "pin_dP0_A": dA[p_in], "pin_dP0_B": dB[p_in], "pin_dh0_A": hA[p_in], "pin_dh0_B": hB[p_in]}
        rows.append(r)
    with open(a.series_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({kk: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for kk, v in r.items()})

    def pinfo(p):
        return None if p is None else {"index": p, "i": int(W.i[p]), "j": int(W.j[p]), "j_from_wall": int(nj - 1 - W.j[p]),
                                       "x_m": float(xy[p, 0]), "y_m": float(xy[p, 1])}
    out = {
        "plan": "plans/active/tooling-nozzle-card-local-mesh.md §5.1 #3c",
        "consult": "notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md (判別 A/B)",
        "date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "run": str(run), "steps": steps, "ni": ni, "nj": nj, "Pt": Pt, "Tt": Tt, "thr_rel": THR,
        "species_record": {"source": th_ref["source"], "how": th_ref["how"], "path": th_ref["path"],
                           "sha256": sha256(th_ref["path"]) if th_ref["path"] else None,
                           "thermoHrefTemp": th_ref["thermoHrefTemp"], "names": th_ref["names"]},
        "solver_convention": {"roMin": roMin, "pMin": pMin, "thermoFloat": thermoFloat, "TMIN": TMIN, "TMAX": TMAX,
                              "h0_includes_k": 0, "composition_B": "float32 (dependentVariables_d.cu useHybrid)",
                              "output_timing": "原始量・Ht は step の頭 (更新前・BC 上書き後)、保存量は更新後 (main.cpp advanceImplicitSteady)"},
        "location_note": "index = i·nj + j (j = 0 軸、j = nj − 1 壁); 分類の優先順 inlet_col (i=0) → outlet_col (i=ni−1) → wall_row → axis_row → interior",
        "rule": RULE,
        "judge": {"max_abs_A_minus_B_over_Pt_all_nodes": dmax, "all_le_thr": cond_small,
                  "n_exceed_A": {st: keep[st]["nexA"] for st in steps}, "n_exceed_B": {st: keep[st]["nexB"] for st in steps},
                  "n_exceed_A_and_B": {st: keep[st]["exB_on_A"] for st in steps},
                  "B_exceed_remains_all_t": remains, "B_exceed_vanishes_all_t": vanish, "AB_diff_explains_all_t": explains,
                  "explained_fraction": {st: per[st]["P0"]["explained_fraction"] for st in steps},
                  "verdict": verdict,
                  "aux_not_for_verdict_max_abs_diff": aux},
        "series": {"csv": Path(a.series_csv).name, "S_defined_at": last, "S_def": "P0_A > Pt (1 + 1e−6) at the last step",
                   "n_S": int(S.sum()), "n_S_axis_near": int(S_ax.sum()), "n_S_interior": int(S_in.sum()),
                   "S_axis_near_def": "S ∩ j ≤ 2 (入口・出口の列を除く)", "S_interior_def": "S ∩ 3 ≤ j ≤ nj − 2 (同)",
                   "rep_axis_near": pinfo(p_ax), "rep_interior": pinfo(p_in), "rep_def": "各部分で最後の時点の P0_A が最大の節点",
                   "h0_in": h0_in, "h0_in_def": "h_mix(Tt, 入口の Y) (同じ物性の記録・datum)"},
        "per_step": per,
        "elapsed_s": time.time() - t_start,
    }
    Path(a.out_json).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"[p0_ab] verdict: {verdict}")
    print(f"[p0_ab] max|A−B|/Pt all nodes per step: " + ", ".join(f"{st}:{dmax[st]:.3e}" for st in steps))
    print(f"[p0_ab] exceed A/B per step: " + ", ".join(f"{st}:{keep[st]['nexA']}/{keep[st]['nexB']}" for st in steps))
    return 0


if __name__ == "__main__":
    sys.exit(main())
