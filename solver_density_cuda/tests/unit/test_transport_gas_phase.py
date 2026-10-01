#!/usr/bin/env python3
"""凝縮 carrier の気相組成の輸送物性 (plans/active/condensation-two-phase-transport.md §4.1, §5.1 #3, §6 物性) の試験。
物性だけを評価する CFD 0 step のハーネス (FORGE_TRANSPORT_PROBE; main.cpp runTransportProbe) を使う。

  python3 solver_density_cuda/tests/unit/test_transport_gas_phase.py --forge BIN --base-forge OLD_BIN [--seed-run DIR] [--keep]

合格条件 (plan §6 物性; 実装前に固定):
  (G0) g = 0 で現行とビット一致: 凝縮 OFF / 凝縮 ON (rog = 0) × 表引き ON / OFF (FORGE_TRANSPORT_TABLE=0) の 4 構成で、
       同じ状態表の probe 出力 (.bin: 読み戻した入力・セル/ghost の vis_lam・thermCond・double の μ・λ・X・壁の μ_w・λ_w・状態表 (5)(6))
       が旧バイナリとバイト一致。
  (G1) g > 0 の湿潤セル・湿潤壁: 気相組成 (Y_s/(1−g), (Y_w−g)/(1−g), 負は 0) を**独立参照** transport_reference.py に渡した値と
       - double 評価 (セル・ghost・壁・状態表 (5)) 相対 ≤1e-12、
       - float 格納 (vis_lam・thermCond のセル・ghost、壁の μ_w・λ_w、表引き状態表 (6)) 相対 ≤1e-5 (表範囲外 T の double 委譲を含む)。
       未更新・NaN/Inf・非正 0 件、どの状態もセル・ghost・壁の解点に乗る。
       判別: 同じ入力を総組成 (液を蒸気として数える現行の作り方) で参照した値との差が湿潤状態で 1e-5 を超えること (試験が組成の違いを見分ける)。
  (G2) 化学種の分子拡散係数 D_s (Fick 拡散と同じ組成・関数; speciesDmixProbe_d_wrapper) が、気相組成の X を渡した独立 double 参照
       (Chapman–Enskog 二元 + Neufeld Ω(1,1) + 混合平均 (1−X_i)/Σ X_j/D_ij) と相対 ≤1e-5 (セル・ghost)。判別は G1 と同じ。
構成: seed (case/44 run_0509, node 軸対称) の nozzle.h5、種 [MIXDRY, H2O]、TP carrier 凝縮 (condGasSpecies 1)、viscMethod 2。
  (a) seed の外部 DB (MIXDRY・H2O とも LJ あり)、transport {MIXDRY: kinetic, H2O: custom:h2o_iapws_cea_v1} — G1・G2
  (b) MIXDRY を lump {N2, O2, AR, CO2} (mole) にした va3 相当、transport 全 cea + H2O custom — G1
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, json, math, os, shutil, subprocess, sys, tempfile

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_transport_gpu as tg  # noqa: E402  (Runner・read_probe・seed の探索を共用)
from transport_reference import Reference  # noqa: E402

FAIL = 0
COND = {"condensation": 1, "nCondSpecies": 1, "condModel": 1, "condGasSpecies": 1}
T_LIST = [60.0, 120.0, 149.9, 180.0, 207.6, 253.15, 300.0, 600.0, 1000.0, 2000.0, 16000.0]   # 表範囲 [150, 15000] K の外を含む


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def rel(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


_H2O_MW = {}


def builtin_h2o_mw(R, exe):
    """exe の内蔵気相 H2O の MW (--resolve-species の記録から)。段 3 (#13-3) で CEA の 0.01801528 に変わった (旧 0.0180153)。"""
    if exe not in _H2O_MW:
        d = R.make("h2o_mw", ["H2O"], None, visc=0)
        r = subprocess.run([exe, "--resolve-species"], cwd=d, env=R.env(), capture_output=True, text=True, timeout=600)
        rec = [f for f in os.listdir(d) if f.startswith("resolved_species_")]
        assert r.returncode == 0 and rec, r.stdout + r.stderr
        _H2O_MW[exe] = float(yaml.safe_load(open(os.path.join(d, rec[0])))["species"][0]["MW"])
    return _H2O_MW[exe]


def seed_db_paired(seed_db, mw):
    """seed の外部 DB の H2O の MW を内蔵の値にした写し。凝縮 ON では気液ペアの契約 (気相 H2O = 内蔵とビット一致) が要る。
    seed (case/44 run_0509) の外部 H2O は段 3 (#13-3) 前の内蔵 H2O と同じ係数・MW だったので、段 3 後のバイナリでは MW だけ差し替える
    (係数は段 3 でも不変; 差し替え後に内蔵と一致しなければソルバが起動時に拒否する)。"""
    db = json.loads(json.dumps(seed_db))
    db["H2O"]["MW"] = float(mw)
    return db


def make(R, tag, species, transport, cond, keep_db=False, db=None, h2o_pair_mw=None):
    d = R.make(tag, species, transport, db=db, keep_db_of_seed=keep_db)
    if cond and keep_db and h2o_pair_mw is not None:
        dbp = os.path.join(d, "species_db.yaml")
        paired = seed_db_paired(yaml.safe_load(open(dbp)), h2o_pair_mw)
        with open(dbp, "w") as f:
            yaml.safe_dump(paired, f, sort_keys=False)
    p = os.path.join(d, "solverConfig.yaml")
    cfg = yaml.safe_load(open(p))
    if cond:
        cfg["condensation"] = dict(COND)
    else:
        cfg.pop("condensation", None)
    with open(p, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    return d


def probe(R, exe, d, states, nS, liq, table=True, name="states.txt"):
    """states: [(T, ro, Y[, g])]。liq=True なら 1 行目に 1 を付け各行の末尾に g。戻り値 (layout+data, bin bytes, err)。"""
    p = os.path.join(d, name)
    with open(p, "w") as f:
        f.write(f"{len(states)} {nS}" + (" 1" if liq else "") + "\n")
        for st in states:
            vals = [st[0], st[1]] + list(st[2]) + ([st[3]] if liq else [])
            f.write(" ".join("%.17g" % v for v in vals) + "\n")
    e = R.env()
    e["FORGE_TRANSPORT_PROBE"] = p
    if not table:
        e["FORGE_TRANSPORT_TABLE"] = "0"
    r = subprocess.run([exe], cwd=d, env=e, capture_output=True, text=True, timeout=1800)
    open(os.path.join(d, name + ".log"), "w").write(r.stdout + r.stderr)
    if r.returncode != 0:
        return None, None, (r.stdout + r.stderr)[-3000:]
    return tg.read_probe(p), open(p + ".bin", "rb").read(), ""


def y_of_x(ref, X):
    w = [x * m for x, m in zip(X, ref.mw)]
    return [v / sum(w) for v in w]


# ------------------------------------------------------------------ 独立参照: 分子拡散係数 (double)
def omega11(Ts):
    Ts = min(max(Ts, 0.3), 100.0)
    return 1.06036 * Ts**-0.15610 + 0.19300 * math.exp(-0.47635 * Ts) + 1.03587 * math.exp(-1.52996 * Ts) + 1.76474 * math.exp(-3.89411 * Ts)


def d_binary(a, b, T, P):
    Mi, Mj = a["MW"] * 1e3, b["MW"] * 1e3
    sig = 0.5 * (a["sigma"] + b["sigma"])
    eps = math.sqrt(a["eps"] * b["eps"])
    return 1.8583e-3 * math.sqrt(T**3 * (1 / Mi + 1 / Mj)) / ((P / 101325.0) * sig * sig * omega11(T / eps)) * 1e-4


def d_mix(sps, X, T, P):
    n = len(sps)
    out = []
    for i in range(n):
        den = sum(X[j] / max(d_binary(sps[i], sps[j], T, P), 1e-30) for j in range(n) if j != i)
        out.append(d_binary(sps[i], sps[i], T, P) if den < 1e-30 else (1 - X[i]) / den)
    return out


# ------------------------------------------------------------------ 組成
def gas_Y(roY, rog, ro, iw):
    """読み戻した float 入力から気相組成 (独立に作る): Y = max(ρY/ρ, 0)、Y_w ← max(Y_w − ρg/ρ, 0)、和で正規化。"""
    Y = [max(float(v) / float(ro), 0.0) for v in roY]
    Y[iw] = max(Y[iw] - float(rog) / float(ro), 0.0)
    t = sum(Y)
    return [v / t for v in Y]


def total_Y(roY, ro):
    Y = [max(float(v) / float(ro), 0.0) for v in roY]
    t = sum(Y)
    return [v / t for v in Y]


# ------------------------------------------------------------------ G0
def g0_bit_identity(R, a, species, transport, ref, keep_db):
    comps = [[0.99, 0.01], [0.94, 0.06], [0.6, 0.4], [1.0, 0.0]]
    states = [(T, 0.37, y_of_x(ref, X)) for X in comps for T in T_LIST]
    for cond in (False, True):
        for table in (True, False):
            tag = f"G0 cond {'ON (rog=0)' if cond else 'OFF'}, table {'ON' if table else 'OFF'}"
            outs, mws = [], []
            for which, exe in (("new", a.forge), ("old", a.base_forge)):
                mws.append(builtin_h2o_mw(R, exe))
                d = make(R, f"g0_{which}_c{int(cond)}_t{int(table)}", species, transport, cond, keep_db=keep_db,
                         h2o_pair_mw=mws[-1])
                o, b, err = probe(R, exe, d, states, len(species), False, table=table)
                if o is None:
                    check(False, f"{tag}: {which} probe failed: {err}")
                    outs = None
                    break
                outs.append((o, b))
            if outs is None:
                continue
            (on, bn), (oo, bo) = outs
            nd = int(np.count_nonzero(np.frombuffer(bn, np.uint8) != np.frombuffer(bo, np.uint8))) if len(bn) == len(bo) else -1
            h = on["h"]
            if cond and keep_db and mws[0] != mws[1]:
                # 凝縮 ON は外部 H2O を各バイナリの内蔵 H2O に合わせる (気液ペアの契約) ので、内蔵 H2O の MW が違うバイナリ同士
                # (段 3 #13-3 の前後) は入力が違う → ビット一致は求めず差の件数を記録する
                check(len(bn) == len(bo) and h.get("hasLiq", 0) == 0,
                      f"{tag}: built-in H2O MW differs between binaries ({mws[1]!r} -> {mws[0]!r}, #13-3); no bit identity required, "
                      f"differing bytes {nd} of {len(bn)}")
                continue
            check(len(bn) == len(bo) and nd == 0 and h.get("hasLiq", 0) == 0,
                  f"{tag}: probe output bytes new vs old {len(bn)}/{len(bo)}, differing {nd} "
                  f"(cells+ghosts {h['nCells_all']}, passes {h['passes']}, table {h['table']}, iw {h.get('iw')})")


# ------------------------------------------------------------------ G1 / G2
def g1(R, a, tag, species, transport, ref, keep_db, dsp=None, db=None):
    iw = 1
    comps = [[0.99, 0.01], [0.94, 0.06], [0.6, 0.4]]
    fr = [0.0, 0.3, 0.9, 0.99, 1.0, 1.2]    # g / Y_w (1.2 = 蒸気上限違反 → 気相の水 0)
    states = []
    for X in comps:
        Y = y_of_x(ref, X)
        for f in fr:
            for T in T_LIST:
                states.append((T, 0.37, Y, f * Y[iw]))
    d = make(R, "g1_" + tag, species, transport, True, keep_db=keep_db, db=db, h2o_pair_mw=builtin_h2o_mw(R, a.forge))
    out, _, err = probe(R, a.forge, d, states, len(species), True)
    if out is None:
        check(False, f"G1 {tag}: probe failed: {err}")
        return
    h = out["h"]
    nA, nC, nS, K = h["nCells_all"], h["nCells"], h["nSpecies"], h["K"]
    check(h.get("hasLiq") == 1 and h.get("iw") == iw and K == len(states),
          f"G1 {tag}: probe ran with liquid (iw {h.get('iw')}, states {K}, passes {h['passes']}, cells {nC}, ghosts {nA - nC})")
    cache = {}

    def refs(T, ro, roY, rog):
        key = (np.float32(T).tobytes(), np.float32(ro).tobytes(), np.asarray(roY, np.float32).tobytes(), np.float32(rog).tobytes())
        if key not in cache:
            rg = ref.state_Y(float(T), gas_Y(roY, rog, ro, iw))
            rt = ref.state_Y(float(T), total_Y(roY, ro))
            cache[key] = (rg, rt)
        return cache[key]

    W = dict(dmu=0.0, dlam=0.0, dX=0.0, wmu=0.0, wlam=0.0, fmu=0.0, flam=0.0, fwmu=0.0, fwlam=0.0, D=0.0, disc=0.0, discD=0.0,
             Dfar=0.0, Dnear=0.0, Dnear_ratio=0.0)
    bad = 0
    cov = {"cell": set(), "ghost": set(), "wall": set()}
    wet_wall = 0
    for q in out["passes"]:
        roY, rog = q["roY"], q["rog"]
        rmu, rlam, rX, tmu = np.empty(nA), np.empty(nA), np.empty((nA, h["nReal"])), np.empty(nA)
        Dref = np.full((nS, nA), np.nan)
        Dtot = np.full((nS, nA), np.nan)
        Xg = np.full((nS, nA), np.nan)
        for i in range(nA):
            rg, rt = refs(q["T"][i], q["ro"][i], roY[:, i], rog[i])
            rmu[i], rlam[i], rX[i], tmu[i] = rg["mu"], rg["lam"], rg["Xreal"], rt["mu"]
            if dsp is not None:
                T, P = float(q["T"][i]), float(q["P"][i])
                Xg[:, i] = ref.X_from_Y(gas_Y(roY[:, i], rog[i], q["ro"][i], iw))
                Dref[:, i] = d_mix(dsp, Xg[:, i], T, P)
                Dtot[:, i] = d_mix(dsp, ref.X_from_Y(total_Y(roY[:, i], q["ro"][i])), T, P)
        for arr in (q["vis"], q["cond"], q["mu"], q["lam"], q["wmu_f"], q["wlam_f"], q["wmu"], q["wlam"]):
            bad += int(np.sum(~np.isfinite(arr)) + np.sum(arr <= 0))
        W["dmu"] = max(W["dmu"], float(rel(q["mu"], rmu).max()))
        W["dlam"] = max(W["dlam"], float(rel(q["lam"], rlam).max()))
        W["dX"] = max(W["dX"], float(np.abs(q["X"] - rX).max()))
        W["fmu"] = max(W["fmu"], float(rel(q["vis"], rmu).max()))
        W["flam"] = max(W["flam"], float(rel(q["cond"], rlam).max()))
        wet = (rog > 0) & (roY[iw] > 0)
        W["disc"] = max(W["disc"], float(rel(rmu[wet], tmu[wet]).max()) if wet.any() else 0.0)
        wc = q["wc"]
        if len(wc):
            W["wmu"] = max(W["wmu"], float(rel(q["wmu"], rmu[wc]).max()))
            W["wlam"] = max(W["wlam"], float(rel(q["wlam"], rlam[wc]).max()))
            W["fwmu"] = max(W["fwmu"], float(rel(q["wmu_f"], rmu[wc]).max()))
            W["fwlam"] = max(W["fwlam"], float(rel(q["wlam_f"], rlam[wc]).max()))
            wet_wall += int(np.sum(wet[wc]))
        if dsp is not None:
            bad += int(np.sum(~np.isfinite(q["D"])) + np.sum(q["D"] <= 0))
            W["D"] = max(W["D"], float(rel(q["D"], Dref).max()))
            W["discD"] = max(W["discD"], float(rel(Dref[:, wet], Dtot[:, wet]).max()) if wet.any() else 0.0)
            # 内訳 (情報; 合否には使わない): 混合平均 (1−X_i)/Σ X_j/D_ij は X_i → 1 で float の桁落ち (|δX| ~ ε₃₂) が
            # 相対 ε₃₂/(1−X_i) の誤差になる。1−X_i ≥ 1e-2 の点と、それ未満の点の誤差/桁落ち限界 (6e-8/(1−X_i)) の比を出す。
            eD = rel(q["D"], Dref)
            far = (1.0 - Xg) >= 1e-2
            if far.any():
                W["Dfar"] = max(W["Dfar"], float(eD[far].max()))
            if (~far).any():
                W["Dnear"] = max(W["Dnear"], float(eD[~far].max()))
                W["Dnear_ratio"] = max(W["Dnear_ratio"], float((eD[~far] / (6e-8 / np.maximum(1.0 - Xg[~far], 1e-30))).max()))
        cov["cell"] |= set(q["st"][:nC].tolist())
        cov["ghost"] |= set(q["st"][nC:].tolist())
        cov["wall"] |= set(q["st"][wc].tolist())
    # 状態表 (5) double・(6) 表引き: 入力は double の Y と g (表引きは float に丸めたもの)
    s5m = s5l = s6m = s6l = 0.0
    for k, (T, _, Y, g) in enumerate(states):
        Yg = list(Y)
        Yg[iw] = max(Yg[iw] - g, 0.0)
        r = ref.state_Y(T, Yg)
        s5m = max(s5m, float(rel(out["states"]["mu"][k], r["mu"])))
        s5l = max(s5l, float(rel(out["states"]["lam"][k], r["lam"])))
        if h.get("table", 0):
            Yf = [float(np.float32(v)) for v in Y]
            Yf[iw] = float(np.float32(np.float32(Yf[iw]) - np.float32(g))) if np.float32(g) > 0 else Yf[iw]
            Yf[iw] = max(Yf[iw], 0.0)
            r6 = ref.state_Y(float(np.float32(T)), Yf)
            s6m = max(s6m, float(rel(out["states"]["tmu"][k], r6["mu"])))
            s6l = max(s6l, float(rel(out["states"]["tlam"][k], r6["lam"])))
    check(max(W["dmu"], W["dlam"], W["wmu"], W["wlam"], s5m, s5l) <= 1e-12 and W["dX"] <= 1e-12,
          f"G1 {tag}: double vs independent gas-phase reference: cells+ghosts mu {W['dmu']:.2e} lam {W['dlam']:.2e} |dX| {W['dX']:.2e}, "
          f"wall mu {W['wmu']:.2e} lam {W['wlam']:.2e}, states (5) mu {s5m:.2e} lam {s5l:.2e} (<= 1e-12)")
    check(max(W["fmu"], W["flam"], W["fwmu"], W["fwlam"], s6m, s6l) <= 1e-5,
          f"G1 {tag}: float storage ({'table' if h['table'] else 'double'} path; T outside [150, 15000] K delegates to double) "
          f"vs gas-phase reference: vis_lam {W['fmu']:.2e} thermCond {W['flam']:.2e}, wall mu_w {W['fwmu']:.2e} lam_w {W['fwlam']:.2e}, "
          f"table states (6) mu {s6m:.2e} lam {s6l:.2e} (<= 1e-5); wet wall points {wet_wall}")
    check(W["disc"] > 1e-5,
          f"G1 {tag}: discrimination — gas-phase vs total-composition reference on wet states max |d mu| {W['disc']:.2e} (> 1e-5)")
    allk = set(range(K))
    check(bad == 0 and cov["cell"] == allk and cov["ghost"] == allk and cov["wall"] == allk and wet_wall > 0,
          f"G1 {tag}: NaN/Inf/non-positive/non-updated {bad}; every state on cells/ghosts/boundary points "
          f"({len(cov['cell'])}/{len(cov['ghost'])}/{len(cov['wall'])} of {K})")
    if dsp is not None:
        print(f"[INFO] G2 {tag}: D_s error breakdown — points with 1-X_s >= 1e-2: {W['Dfar']:.2e}; points with 1-X_s < 1e-2 "
              f"(near-pure; float cancellation in 1-X_s of thermo_Dmix_species_f): {W['Dnear']:.2e} = "
              f"{W['Dnear_ratio']:.2f} x (6e-8/(1-X_s))", flush=True)
        check(W["D"] <= 1e-5 and W["discD"] > 1e-5,
              f"G2 {tag}: molecular D_s (Fick path composition) vs independent double gas-phase reference {W['D']:.2e} (<= 1e-5); "
              f"discrimination gas vs total on wet states {W['discD']:.2e} (> 1e-5)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", required=True)
    ap.add_argument("--base-forge", default="")
    ap.add_argument("--seed-run", default="")
    ap.add_argument("--blocksize", type=int, default=256)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    a.forge = os.path.abspath(a.forge)
    if a.base_forge:
        a.base_forge = os.path.abspath(a.base_forge)
    if not a.seed_run:
        for c in (os.path.join(tg.REPO, tg.SEED_REL), os.path.join(os.path.dirname(tg.REPO), "forge", tg.SEED_REL)):
            if os.path.exists(os.path.join(c, "nozzle.h5")):
                a.seed_run = c
                break
    if not a.seed_run:
        print("[FAIL] seed run with nozzle.h5 not found (--seed-run)")
        return 1
    root = tempfile.mkdtemp(prefix="transport_gas_phase_")
    print(f"work dir: {root}\nseed: {a.seed_run}\nforge: {a.forge}\nbase: {a.base_forge or '(none)'}", flush=True)
    R = tg.Runner(a, root)
    # 外部 H2O の MW は --forge の内蔵 H2O に揃える (凝縮 ON の気液ペア; G1 とその参照が同じ値を使う)
    seed_db = seed_db_paired(yaml.safe_load(open(os.path.join(a.seed_run, "species_db.yaml"))), builtin_h2o_mw(R, a.forge))

    # (a) 外部 DB (LJ 既知) — G0・G1・G2
    spA = ["MIXDRY", "H2O"]
    trA = {"MIXDRY": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}
    refA = Reference(spA, trA, seed_db)
    dsp = [{"MW": float(seed_db[n]["MW"]), "sigma": float(seed_db[n]["LJ_sigma"]), "eps": float(seed_db[n]["LJ_eps_kB"])} for n in spA]
    if a.base_forge:
        g0_bit_identity(R, a, spA, trA, refA, keep_db=True)
    else:
        print("[SKIP] G0 bit identity (no --base-forge)")
    g1(R, a, "a_extdb_kinetic_h2ocustom", spA, trA, refA, keep_db=True, dsp=dsp)

    # (b) va3 相当: MIXDRY lump (mole) + builtin H2O、全 cea + H2O custom — G1
    lump = {"name": "MIXDRY", "lump": {"N2": 0.708873, "O2": 0.230376, "AR": 0.00850387, "CO2": 0.0522474}, "basis": "mole"}
    spB = [lump, "H2O"]
    trB = {"N2": "cea", "O2": "cea", "AR": "cea", "CO2": "cea", "H2O": "custom:h2o_iapws_cea_v1"}
    refB = Reference(spB, trB)
    g1(R, a, "b_lump_cea", spB, trB, refB, keep_db=False)

    if not a.keep:
        shutil.rmtree(root, ignore_errors=True)
    print("ALL PASS" if FAIL == 0 else f"FAILED: {FAIL}")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
