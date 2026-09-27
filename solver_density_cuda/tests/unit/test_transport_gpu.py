#!/usr/bin/env python3
"""種ごとの輸送物性の GPU 接続 (physProp.transport, 段 2) の試験 — 物性だけを評価する CFD 0 step のハーネス。
plans/active/thermophysics-solver-owned-species-db.md §5.1 #5t2-2 (合格条件は実装前に固定したもの; codex diagnose
notes/reviews/2026-09-27-transport-stage2-gpu-diagnose.md)。

  python3 solver_density_cuda/tests/unit/test_transport_gpu.py --forge BIN [--base-forge OLD_BIN] [--seed-run DIR] [--keep]

forge を FORGE_TRANSPORT_PROBE=<states.txt> で起動する (main.cpp runTransportProbe)。初期化後、時間更新をせずに
状態表の (T, ρ, Y) を全セル (ghost 込み nCells_all) と全境界面の解点へ float で書き、実際のセル経路 (gasProperties_d)・
同じ評価関数の double 値・壁モデルの物性関数 (wmles_wall_props) を各 1 回呼び、**device から読み戻した実際の入力**
(float の T・ρ・roY) と出力を書いて終了する。参照は独立実装 transport_reference.py (実装とコードを共有しない) に
その実際の入力を渡して作る。メッシュは seed run (既定 case/44 run_0509, node 軸対称 23725 CV・境界面 860) の nozzle.h5。

  (D)  double 評価: セル (ghost 末尾まで)・壁・状態表 (float を経由しない double の Y) の μ・λ が独立参照と相対 ≤1e-12、
       展開後の実種モル分率 X が絶対 ≤1e-12。
  (S)  float 格納: vis_lam・thermCond (セル・ghost) と壁の μ_w・λ_w が独立 double 参照と相対 ≤1e-5、
       かつ参照を float に丸めた値と ≤2 ULP。
  (U)  未更新 (事前に NaN で埋めた値が残る)・NaN/Inf・非正の μ/λ が 0 件。どの状態もセル・ghost・境界面の解点に乗った。
  (AB) 判別 (T 400 K, ρ 1, 両実種 cea): A = full X_N2 0.4・X_He 0.6 / B = lump {N2 0.5, He 0.5} + 独立 He, X_L 0.8・X_He 0.2。
       B の展開後 X_He は 0.6 (基底取り違えなら 0.887)。A と B、列挙順を変えた A'・B' は double 入力 (状態表) で ≤1e-12、
       別々に float 化した roY 経由のセル値どうしは ≤1e-5 (それぞれは (D)(S) で自分の独立参照と照合済み)。
  範囲: 単成分 (roY 無し)・重複 lump・実種 12 と上限 32・ゼロ分率、T = 200/253.15/400/500/600/700/1000/2000 K と
       各フィット・接続・NASA 区間の境界 (境界そのものと両隣の float)。
  (B)  (--base-forge) physProp.transport なしでは同一の固定入力 (seed の保存量; res_0 の T・ρ・roY がビット一致することを
       確認してから) で vis_lam・thermCond が旧バイナリとビット一致 (viscMethod 2 の多成分・単成分、viscMethod 1)。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, json, math, os, shutil, subprocess, sys, tempfile

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, "..", ".."))
REPO = os.path.normpath(os.path.join(SOLVER, ".."))
sys.path.insert(0, HERE)
from transport_reference import Reference, TR  # noqa: E402

SEED_REL = os.path.join("case", "44.vitiated_air_wt", "run_0509_va3_M4.19_Lc8_dry_lumpX")
BASE_T = [200.0, 253.15, 400.0, 500.0, 600.0, 700.0, 1000.0, 2000.0]
FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def f32(x):
    return float(np.float32(x))


def ulp_dist(a, b):
    """正の float32 どうしの ULP 距離。"""
    ia = np.asarray(a, dtype=np.float32).view(np.int32).astype(np.int64)
    ib = np.asarray(b, dtype=np.float32).view(np.int32).astype(np.int64)
    return np.abs(ia - ib)


def rel(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


# ------------------------------------------------------------------ 試験用の外部 DB 種
N2_DATA = None


def n2_thermo():
    global N2_DATA
    if N2_DATA is None:
        n2 = [e for e in yaml.safe_load(open(os.path.join(SOLVER, "data", "species", "forge_species_v1.yaml")))["species"]
              if e["id"] == "N2"][0]
        N2_DATA = n2
    return N2_DATA


def ext_species(i, model):
    """外部 DB の試験種 XS<i> (熱物性は N2 の係数、MW・LJ・フィットは i でずらす)。"""
    n2 = n2_thermo()
    e = {"MW": 0.012 + 0.0027 * i, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
         "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"],
         "LJ_sigma": 3.0 + 0.071 * i, "LJ_eps_kB": 60.0 + 7.3 * i}
    if model == "kinetic" and i % 3 == 0:
        e["LJ_dipole"] = 0.3 + 0.1 * i
    if model == "fit":
        d = 0.02 * i
        e["transport_fit"] = {"V": [[150.0, 1000.0, 0.62526577, -31.779652, -1640.7983, 1.7454992 + d],
                                    [1000.0, 7000.0, 0.87395209, 561.52222, -173948.09, -0.39335958 + d]],
                              "C": [[150.0, 1000.0, 0.85439436, 105.73224, -12347.848, 0.47793128 - d],
                                    [1000.0, 7000.0, 0.88407146, 133.57293, -11429.64, 0.24417019 - d]]}
    return e


def fit_edges(ref):
    """実種のフィット区間・接続点・NASA 区間 (kinetic の c_p) と CEA 相互作用フィットの境界温度。"""
    E = set()
    for e in ref.reals:
        m, sp = e["model"], e["sp"]
        if m in ("cea", "custom:h2o_iapws_cea_v1"):
            for kind in ("V", "C"):
                for lo, hi, _ in TR[(e["trans"], "")][kind]:
                    E |= {lo, hi}
        if m == "custom:h2o_iapws_cea_v1":
            E |= {253.15, 500.0, 700.0}
        if m == "fit":
            for kind in ("V", "C"):
                for r in sp["fit"][kind]:
                    E |= {float(r[0]), float(r[1])}
        if m == "kinetic":
            E |= {sp["Tlo"], sp["Tmid"], sp["Thi"]}
    n = len(ref.reals)
    for a in range(n):
        for b in range(a + 1, n):
            A, B = ref.reals[a], ref.reals[b]
            if A["model"] == "kinetic" and B["model"] == "kinetic":
                continue
            if A["trans"] and B["trans"]:
                for k in ((A["trans"], B["trans"]), (B["trans"], A["trans"])):
                    if k in TR and TR[k]["V"]:
                        for lo, hi, _ in TR[k]["V"]:
                            E |= {lo, hi}
                        break
    return sorted(E)


def t_list(ref):
    T = {f32(t) for t in BASE_T}
    for b in fit_edges(ref):
        fb = np.float32(b)
        T |= {float(fb), float(np.nextafter(fb, np.float32(-np.inf))), float(np.nextafter(fb, np.float32(np.inf)))}
    return sorted(T)


# ------------------------------------------------------------------ run の準備と実行
class Runner:
    def __init__(self, a, root):
        self.a, self.root, self.n = a, root, 0
        self.cfg0 = yaml.safe_load(open(os.path.join(a.seed_run, "solverConfig.yaml")))
        self.bc0 = yaml.safe_load(open(os.path.join(a.seed_run, "bcondConfig.yaml")))

    def make(self, tag, species, transport, db=None, visc=2, keep_db_of_seed=False, nstep=0):
        self.n += 1
        d = os.path.join(self.root, f"{self.n:02d}_{tag}")
        os.makedirs(d)
        for fn in ("nozzle.h5", "probe.yaml"):
            shutil.copy(os.path.join(self.a.seed_run, fn), d)
        cfg = json.loads(json.dumps(self.cfg0))
        ph = cfg["physProp"]
        ph["viscMethod"] = visc
        ph["visc"], ph["thermCond"] = 1.8e-5, 0.026
        ph["species"] = species
        ph.pop("speciesDBFile", None)
        ph.pop("transport", None)
        if keep_db_of_seed:
            shutil.copy(os.path.join(self.a.seed_run, "species_db.yaml"), d)
            ph["speciesDBFile"] = "species_db.yaml"
        if db:
            with open(os.path.join(d, "db.yaml"), "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
            ph["speciesDBFile"] = "db.yaml"
        if transport:
            ph["transport"] = transport
        cfg["time"]["last"]["nStepOuter"] = nstep
        cfg["time"]["outStepInterval"] = 1000000
        cfg["output"] = {"level": 1, "extraFields": ["vis_lam", "thermCond", "T"]}
        with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False)
        bc = json.loads(json.dumps(self.bc0))
        fl = bc["inlet"]["floats"]
        for k in [k for k in fl if k.startswith("X")]:
            del fl[k]
        for s in range(len(species)):
            fl[f"X{s}"] = 1.0 if s == 0 else 0.0
        with open(os.path.join(d, "bcondConfig.yaml"), "w") as f:
            yaml.safe_dump(bc, f, sort_keys=False)
        return d

    def env(self):
        e = dict(os.environ)
        if self.a.blocksize:
            e["FORGE_CUDA_BLOCKSIZE"] = str(self.a.blocksize)
        return e

    def probe(self, d, states, nS):
        p = os.path.join(d, "states.txt")
        with open(p, "w") as f:
            f.write(f"{len(states)} {nS}\n")
            for T, ro, Y in states:
                f.write(" ".join("%.17g" % v for v in [T, ro] + list(Y)) + "\n")
        e = self.env()
        e["FORGE_TRANSPORT_PROBE"] = p
        r = subprocess.run([self.a.forge], cwd=d, env=e, capture_output=True, text=True, timeout=600)
        open(os.path.join(d, "forge_probe.log"), "w").write(r.stdout + r.stderr)
        if r.returncode != 0:
            return None, (r.stdout + r.stderr)[-3000:]
        return read_probe(p), ""


def read_probe(path):
    h = json.load(open(path + ".json"))
    assert h["sizeof_flow_float"] == 4
    nA, nS, nR, K, hasY = h["nCells_all"], h["nSpecies"], h["nReal"], h["K"], h["hasRoY"]
    buf = open(path + ".bin", "rb").read()
    off = [0]

    def take(dt, n):
        a = np.frombuffer(buf, dtype=dt, count=n, offset=off[0])
        off[0] += a.nbytes
        return a

    passes = []
    for p in range(h["passes"]):
        q = {}
        q["st"] = take(np.int32, nA)
        q["T"] = take(np.float32, nA)
        q["ro"] = take(np.float32, nA)
        q["roY"] = take(np.float32, nA * nS if hasY else 0).reshape((nS if hasY else 0), nA)
        q["vis"] = take(np.float32, nA)
        q["cond"] = take(np.float32, nA)
        q["mu"] = take(np.float64, nA)
        q["lam"] = take(np.float64, nA)
        q["X"] = take(np.float64, nA * nR).reshape(nA, nR)
        nw = h["wall_counts"][p]
        q["wb"] = take(np.int32, nw)
        q["wc"] = take(np.int64, nw)
        q["wmu_f"] = take(np.float32, nw)
        q["wlam_f"] = take(np.float32, nw)
        q["wmu"] = take(np.float64, nw)
        q["wlam"] = take(np.float64, nw)
        passes.append(q)
    s = {"mu": take(np.float64, K), "lam": take(np.float64, K), "X": take(np.float64, K * nR).reshape(K, nR)}
    assert off[0] == len(buf), (off[0], len(buf))
    return {"h": h, "passes": passes, "states": s}


# ------------------------------------------------------------------ 1 構成の検査
class Result:
    pass


def run_config(R, tag, species, transport, comps, db=None, ro=0.37, T_override=None):
    """comps: 輸送種のモル分率の組の列。戻り値: 状態ごとの double (状態表経由) とセル経路 (float roY 経由) の値。"""
    ref = Reference(species, transport, db)
    nS = len(species)
    Ts = T_override if T_override is not None else t_list(ref)
    states, keys = [], []
    for ci, X in enumerate(comps):
        # Y = X M / Σ X M (輸送種の分子量は参照側の独自値)
        w = [x * m for x, m in zip(X, ref.mw)]
        Y = [v / sum(w) for v in w]
        for T in Ts:
            states.append((T, ro, Y))
            keys.append((ci, T))
    d = R.make(tag, species, transport, db)
    out, err = R.probe(d, states, nS)
    if out is None:
        check(False, f"{tag}: forge probe failed: {err}")
        return None
    h = out["h"]
    nA, nC, nR, K = h["nCells_all"], h["nCells"], h["nReal"], h["K"]
    check(nR == len(ref.reals) and K == len(states),
          f"{tag}: probe ran (nCells {nC}, nCells_all {nA} = {nA - nC} ghosts, boundary planes {sum(h['wall_counts'][:1])}, "
          f"real species {nR}, states {K}, passes {h['passes']}, T values {len(Ts)} in [{min(Ts):.6g}, {max(Ts):.6g}])")

    cache = {}

    def ref_of(T, ro_, roY):
        key = (np.float32(T).tobytes(), np.float32(ro_).tobytes(), np.asarray(roY, np.float32).tobytes())
        if key not in cache:
            if nS >= 2 and h["hasRoY"]:
                Y = [float(v) / float(ro_) for v in roY]
            else:
                Y = [1.0]
            cache[key] = ref.state_Y(float(T), Y)
        return cache[key]

    worst = {"cell_mu": 0.0, "cell_lam": 0.0, "cell_X": 0.0, "ghost_mu": 0.0, "ghost_lam": 0.0,
             "wall_mu": 0.0, "wall_lam": 0.0, "fcell_mu": 0.0, "fcell_lam": 0.0, "fwall_mu": 0.0, "fwall_lam": 0.0}
    ulp_worst = 0
    bad = 0
    cov = {"cell": set(), "ghost": set(), "wall": set()}
    cell_by_state = {}
    for q in out["passes"]:
        roY = q["roY"]
        refmu = np.empty(nA)
        reflam = np.empty(nA)
        refX = np.empty((nA, nR))
        for i in range(nA):
            r = ref_of(q["T"][i], q["ro"][i], roY[:, i] if roY.shape[0] else [])
            refmu[i], reflam[i] = r["mu"], r["lam"]
            refX[i] = r["Xreal"]
        for arr in (q["vis"], q["cond"], q["mu"], q["lam"], q["wmu_f"], q["wlam_f"], q["wmu"], q["wlam"]):
            bad += int(np.sum(~np.isfinite(arr)) + np.sum(arr <= 0))
        bad += int(np.sum(~np.isfinite(q["X"])))
        em, el = rel(q["mu"], refmu), rel(q["lam"], reflam)
        worst["cell_mu"] = max(worst["cell_mu"], float(np.nanmax(em[:nC])))
        worst["cell_lam"] = max(worst["cell_lam"], float(np.nanmax(el[:nC])))
        if nA > nC:
            worst["ghost_mu"] = max(worst["ghost_mu"], float(np.nanmax(em[nC:])))
            worst["ghost_lam"] = max(worst["ghost_lam"], float(np.nanmax(el[nC:])))
        worst["cell_X"] = max(worst["cell_X"], float(np.nanmax(np.abs(q["X"] - refX))))
        fm, fl = rel(q["vis"], refmu), rel(q["cond"], reflam)
        worst["fcell_mu"] = max(worst["fcell_mu"], float(np.nanmax(fm)))
        worst["fcell_lam"] = max(worst["fcell_lam"], float(np.nanmax(fl)))
        ulp_worst = max(ulp_worst, int(ulp_dist(q["vis"], refmu.astype(np.float32)).max()),
                        int(ulp_dist(q["cond"], reflam.astype(np.float32)).max()))
        # 壁: 入力は同じ pass のセル配列の解点 ic
        wc = q["wc"]
        wmu = refmu[wc]
        wlam = reflam[wc]
        worst["wall_mu"] = max(worst["wall_mu"], float(np.max(rel(q["wmu"], wmu))) if len(wc) else 0.0)
        worst["wall_lam"] = max(worst["wall_lam"], float(np.max(rel(q["wlam"], wlam))) if len(wc) else 0.0)
        worst["fwall_mu"] = max(worst["fwall_mu"], float(np.max(rel(q["wmu_f"], wmu))) if len(wc) else 0.0)
        worst["fwall_lam"] = max(worst["fwall_lam"], float(np.max(rel(q["wlam_f"], wlam))) if len(wc) else 0.0)
        if len(wc):
            ulp_worst = max(ulp_worst, int(ulp_dist(q["wmu_f"], wmu.astype(np.float32)).max()),
                            int(ulp_dist(q["wlam_f"], wlam.astype(np.float32)).max()))
        cov["cell"] |= set(q["st"][:nC].tolist())
        cov["ghost"] |= set(q["st"][nC:].tolist())
        cov["wall"] |= set(q["st"][wc].tolist())
        for i in range(nC):
            cell_by_state.setdefault(int(q["st"][i]), (float(q["mu"][i]), float(q["lam"][i]), q["X"][i].copy()))
    # 状態表 (double の Y を直接)
    smu, slam, sX = out["states"]["mu"], out["states"]["lam"], out["states"]["X"]
    s_em = s_el = s_ex = 0.0
    for k, (T, _, Y) in enumerate(states):
        r = ref.state_Y(T, Y)
        s_em = max(s_em, float(rel(smu[k], r["mu"])))
        s_el = max(s_el, float(rel(slam[k], r["lam"])))
        s_ex = max(s_ex, float(np.max(np.abs(sX[k] - np.asarray(r["Xreal"])))))
        bad += int(not (np.isfinite(smu[k]) and smu[k] > 0 and np.isfinite(slam[k]) and slam[k] > 0))

    check(max(worst["cell_mu"], worst["cell_lam"], worst["ghost_mu"], worst["ghost_lam"], worst["wall_mu"], worst["wall_lam"],
              s_em, s_el) <= 1e-12 and max(worst["cell_X"], s_ex) <= 1e-12,
          f"{tag}: D double vs independent reference: cell mu {worst['cell_mu']:.2e} lam {worst['cell_lam']:.2e}, "
          f"ghost mu {worst['ghost_mu']:.2e} lam {worst['ghost_lam']:.2e}, wall mu {worst['wall_mu']:.2e} lam {worst['wall_lam']:.2e}, "
          f"double-Y states mu {s_em:.2e} lam {s_el:.2e}; |X - Xref| cell {worst['cell_X']:.2e} states {s_ex:.2e} (<= 1e-12)")
    check(max(worst["fcell_mu"], worst["fcell_lam"], worst["fwall_mu"], worst["fwall_lam"]) <= 1e-5 and ulp_worst <= 2,
          f"{tag}: S float storage vs independent double reference: vis_lam {worst['fcell_mu']:.2e} thermCond {worst['fcell_lam']:.2e} "
          f"(cells+ghosts), wall mu_w {worst['fwall_mu']:.2e} lam_w {worst['fwall_lam']:.2e} (<= 1e-5); max ULP vs float(ref) {ulp_worst} (<= 2)")
    allk = set(range(K))
    check(bad == 0 and cov["cell"] == allk and cov["ghost"] == allk and cov["wall"] == allk,
          f"{tag}: U non-updated/NaN/Inf/non-positive {bad}; every state on cells/ghosts/boundary points "
          f"({len(cov['cell'])}/{len(cov['ghost'])}/{len(cov['wall'])} of {K})")
    res = Result()
    res.ref, res.keys, res.states = ref, keys, states
    res.smu, res.slam, res.sX = smu, slam, sX
    res.cell = {keys[k]: v for k, v in cell_by_state.items()}
    res.real = [e["key"] for e in ref.reals]
    return res


def compare(tagA, A, tagB, B, pairs, tol_double=1e-12, tol_float=1e-5):
    """同じ物理組成の 2 構成を突き合わせる。pairs: [(A の組成番号, B の組成番号)]、温度は同じ値どうし。
    X は実種名で並べ替えて比べる。double 入力 (状態表) は ≤tol_double、別々に float 化した roY 経由のセル値は ≤tol_float。"""
    if A is None or B is None:
        check(False, f"AB {tagA} vs {tagB}: missing run")
        return
    ia = {k: i for i, k in enumerate(A.keys)}
    ib = {k: i for i, k in enumerate(B.keys)}
    em = el = ex = fm = fl = 0.0
    n = 0
    for ca, cb in pairs:
        for (c, T), i in ia.items():
            if c != ca or (cb, T) not in ib:
                continue
            j = ib[(cb, T)]
            em = max(em, float(rel(A.smu[i], B.smu[j])))
            el = max(el, float(rel(A.slam[i], B.slam[j])))
            xa = dict(zip(A.real, A.sX[i]))
            xb = dict(zip(B.real, B.sX[j]))
            ex = max(ex, max(abs(xa[r] - xb.get(r, 0.0)) for r in set(xa) | set(xb) if r in xa))
            pa, pb = A.cell[(ca, T)], B.cell[(cb, T)]
            fm = max(fm, float(rel(pa[0], pb[0])))
            fl = max(fl, float(rel(pa[1], pb[1])))
            n += 1
    check(n > 0 and em <= tol_double and el <= tol_double and ex <= tol_double and fm <= tol_float and fl <= tol_float,
          f"AB {tagA} vs {tagB} ({n} states): double input mu {em:.2e} lam {el:.2e} |dX| {ex:.2e} (<= {tol_double:g}); "
          f"via separately float-rounded roY mu {fm:.2e} lam {fl:.2e} (<= {tol_float:g})")


# ------------------------------------------------------------------ (B) 既定経路のビット一致
def bit_identity(R, a):
    cases = [
        ("seed MIXDRY/H2O ext DB, viscMethod 2", dict(species=["MIXDRY", "H2O"], transport=None, keep_db_of_seed=True, visc=2)),
        ("builtin N2/H2O/O2/AR/CO2, viscMethod 2", dict(species=["N2", "H2O", "O2", "AR", "CO2"], transport=None, visc=2)),
        ("builtin N2 single, viscMethod 2", dict(species=["N2"], transport=None, visc=2)),
        ("seed MIXDRY/H2O ext DB, viscMethod 1", dict(species=["MIXDRY", "H2O"], transport=None, keep_db_of_seed=True, visc=1)),
    ]
    for tag, kw in cases:
        outs = []
        for which, binp in (("new", a.forge), ("old", a.base_forge)):
            d = R.make("bit_" + which + "_" + tag.split(",")[0].replace(" ", "_").replace("/", "-"), kw["species"], None,
                       visc=kw["visc"], keep_db_of_seed=kw.get("keep_db_of_seed", False), nstep=0)
            r = subprocess.run([binp], cwd=d, env=R.env(), capture_output=True, text=True, timeout=600)
            open(os.path.join(d, "forge.log"), "w").write(r.stdout + r.stderr)
            f = os.path.join(d, "res_0.h5")
            if r.returncode != 0 or not os.path.exists(f):
                outs.append(None)
                continue
            with h5py.File(f, "r") as h:
                v = {k: np.asarray(h["VALUE"][k]) for k in h["VALUE"].keys()}
            outs.append(v)
        if None in outs:
            check(False, f"B {tag}: a run failed (see logs under the work dir)")
            continue
        n, o = outs
        inputs = ["T", "ro"] + [k for k in n if k.startswith("roY")]
        same_in = all(k in o and np.array_equal(n[k].view(np.uint32) if n[k].dtype == np.float32 else n[k],
                                                 o[k].view(np.uint32) if o[k].dtype == np.float32 else o[k]) for k in inputs)
        outk = ["vis_lam", "thermCond"]
        have = all(k in n and k in o for k in outk)
        nd = {k: int(np.count_nonzero(np.frombuffer(n[k].tobytes(), np.uint8) != np.frombuffer(o[k].tobytes(), np.uint8)))
              for k in outk} if have else {}
        rng = f"vis_lam [{n['vis_lam'].min():.4g}, {n['vis_lam'].max():.4g}]" if have else ""
        check(same_in and have and all(v == 0 for v in nd.values()),
              f"B {tag}: inputs {inputs} bit-identical {same_in}; vis_lam/thermCond differing bytes {nd} over {len(n['T'])} CVs ({rng})")


# ------------------------------------------------------------------ main
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
        for c in (os.path.join(REPO, SEED_REL), os.path.join(os.path.dirname(REPO), "forge", SEED_REL)):
            if os.path.exists(os.path.join(c, "nozzle.h5")):
                a.seed_run = c
                break
    if not a.seed_run:
        print("[FAIL] seed run with nozzle.h5 not found (--seed-run)")
        return 1
    root = tempfile.mkdtemp(prefix="transport_gpu_")
    print(f"work dir: {root}\nseed: {a.seed_run}\nforge: {a.forge}", flush=True)
    R = Runner(a, root)

    # ---- 判別 A/B (T 400 K, ρ 1, 両実種 cea) ----
    L = {"name": "L", "lump": {"N2": 0.5, "He": 0.5}, "basis": "mole"}
    Lr = {"name": "L", "lump": {"He": 0.5, "N2": 0.5}, "basis": "mole"}
    trNHe = {"N2": "cea", "He": "cea"}
    compsA = [[0.4, 0.6], [1.0, 0.0], [0.0, 1.0]]
    compsB = [[0.8, 0.2], [0.0, 1.0]]
    A = run_config(R, "A_full_N2_He", ["N2", "He"], trNHe, compsA, ro=1.0)
    B = run_config(R, "B_lump_N2He_plus_He", [L, "He"], trNHe, compsB, ro=1.0)
    Ar = run_config(R, "A_rev_He_N2", ["He", "N2"], {"He": "cea", "N2": "cea"}, [[0.6, 0.4], [0.0, 1.0]], ro=1.0)
    Br = run_config(R, "B_rev_He_plus_lump", ["He", Lr], {"He": "cea", "N2": "cea"}, [[0.2, 0.8]], ro=1.0)
    if B is not None:
        k = B.keys.index((0, f32(400.0)))
        xhe = dict(zip(B.real, B.sX[k]))["He"]
        cA = A.cell.get((0, f32(400.0))) if A else None
        cB = B.cell.get((0, f32(400.0)))
        check(abs(xhe - 0.6) <= 1e-12,
              f"AB B at 400 K: expanded X_He {xhe:.15f} (mole basis 0.6; mass-basis mix-up would give 0.887308)")
        if A is not None and cA and cB:
            print(f"       400 K: A mu {cA[0]:.12e} lam {cA[1]:.12e} / B mu {cB[0]:.12e} lam {cB[1]:.12e}")
    # A (X_N2 0.4, X_He 0.6) ↔ B (X_L 0.8, X_He 0.2) は組成 0 どうし、純 He は A の 2 と B の 1
    compare("A", A, "B", B, [(0, 0)])
    compare("A (pure He)", A, "B (pure He via the lone He)", B, [(2, 1)])
    # 列挙順: A' = [He, N2] (組成 0 = A の 0、1 = 純 N2 = A の 1)、B' = [He, L(He, N2)] (組成 0 = B の 0)
    compare("A", A, "A' (species order reversed)", Ar, [(0, 0), (1, 1)])
    compare("B", B, "B' (species and lump member order reversed)", Br, [(0, 0)])

    # ---- 単成分 (roY なし) ----
    run_config(R, "single_N2_cea", ["N2"], {"N2": "cea"}, [[1.0]])
    run_config(R, "single_H2O_custom", ["H2O"], {"H2O": "custom:h2o_iapws_cea_v1"}, [[1.0]])

    # ---- 実種 12 (2 lump + 外部 DB 6 種; N2 が 2 つの lump に重複) ----
    db12 = {f"XS{i:02d}": ext_species(i, "fit" if i % 2 else "kinetic") for i in range(1, 7)}
    sp12 = [{"name": "MIXA", "lump": {"N2": 0.70, "O2": 0.22, "AR": 0.01, "CO2": 0.07}, "basis": "mole"},
            {"name": "MIXB", "lump": {"N2": 0.4, "H2O": 0.35, "He": 0.25}, "basis": "mass"}] + list(db12)
    tr12 = {"N2": "cea", "O2": "kinetic", "AR": "cea", "CO2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1", "He": "cea"}
    tr12.update({k: ("fit" if int(k[2:]) % 2 else "kinetic") for k in db12})
    comps12 = [[0.30, 0.20, 0.10, 0.05, 0.15, 0.08, 0.07, 0.05],
               [0.50, 0.00, 0.10, 0.00, 0.15, 0.00, 0.25, 0.00]]   # ゼロ分率を含む
    r12 = run_config(R, "real12", sp12, tr12, comps12, db=db12)
    if r12:
        check(len(r12.real) == 12, f"real12: real species {len(r12.real)} == 12 ({r12.real})")

    # ---- 実種 32 (上限; 2 lump + 単独 9 種、N2 は lump と単独に重複) ----
    db32 = {f"XS{i:02d}": ext_species(i, "fit" if i % 2 else "kinetic") for i in range(1, 27)}
    lumpA = {"N2": 0.5, "O2": 0.2, "AR": 0.05, "CO2": 0.05}
    lumpA.update({f"XS{i:02d}": 0.025 for i in range(1, 9)})
    lumpB = {"H2O": 0.3, "He": 0.2}
    lumpB.update({f"XS{i:02d}": 0.05 for i in range(9, 19)})
    sp32 = [{"name": "LA", "lump": lumpA, "basis": "mole"}, {"name": "LB", "lump": lumpB, "basis": "mass"}, "N2"] + \
           [f"XS{i:02d}" for i in range(19, 27)]
    tr32 = {"N2": "cea", "O2": "cea", "AR": "kinetic", "CO2": "cea", "H2O": "custom:h2o_iapws_cea_v1", "He": "kinetic"}
    tr32.update({k: ("fit" if int(k[2:]) % 2 else "kinetic") for k in db32})
    nsp = len(sp32)
    c1 = [1.0 / nsp] * nsp
    c2 = [0.0] * nsp
    c2[0], c2[2], c2[5], c2[nsp - 1] = 0.4, 0.3, 0.2, 0.1        # LB と多くの単独種がゼロ
    r32 = run_config(R, "real32", sp32, tr32, [c1, c2], db=db32)
    if r32:
        check(len(r32.real) == 32, f"real32: real species {len(r32.real)} == 32 (TRANSPORT_MAX_REAL_SPECIES)")

    # ---- 既定経路のビット一致 ----
    if a.base_forge:
        bit_identity(R, a)
    else:
        print("[SKIP] B bit identity (no --base-forge)")

    if not a.keep:
        shutil.rmtree(root, ignore_errors=True)
    print("ALL PASS" if FAIL == 0 else f"FAILED: {FAIL}")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
