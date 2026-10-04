#!/usr/bin/env python3
"""種ごとの輸送物性の GPU 接続 (physProp.transport, 段 2) の試験 — 物性だけを評価する CFD 0 step のハーネス。
plans/active/thermophysics-solver-owned-species-db.md §5.1 #5t2-2 (合格条件は実装前に固定したもの; codex diagnose
notes/reviews/2026-09-27-transport-stage2-gpu-diagnose.md)。

  python3 solver_density_cuda/tests/unit/test_transport_gpu.py --forge BIN [--base-forge OLD_BIN] [--seed-run DIR] [--keep]
      [--lj-source gri30,svehla1962]

  --lj-source: 内蔵種の LJ の集合リスト (physProp.ljSource; plan §4.10 #14)。全 run の config と独立参照に同じリストを渡す
       (省略時はソルバの既定 [gri30, svehla1962] で config に書かない)。#14-L1 (iv) は集合ごとに回す。

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
       確認してから) で vis_lam・thermCond が旧バイナリとビット一致 (viscMethod 1 の多成分 3 構成、viscMethod 0)。
       旧バイナリと互換性ハッシュが違う構成 (段 3 #13-3 前のバイナリと内蔵 H2O/AR を含む構成) はビット一致を求めず差を記録し、
       経路の同一性は熱物性が不変な構成 (seed の外部 DB、内蔵 N2/O2/CO2) で見る。
  (B2) viscMethod 2 + physProp.transport なし (多成分・単成分) は起動時エラー (res_0 を書かない)。
       2026-09-27 に viscMethod 2 を種ごとの輸送物性へ置き換え (plan §4.3c 案 C)、旧 kinetic 経路 (Wilke 共用 φ) を計算から
       外したので、以前の「viscMethod 2 の transport なしが旧バイナリとビット一致」はこの期待に変えた。
  seed の nozzle.h5 は species 属性を持たない (記録導入前の場) ので、すべての実行に FORGE_ALLOW_UNVERIFIED_SPECIES=1 を付ける
  (#3c 以降、属性なしの場は既定で停止する)。

表引き (#5t2-3, 既定で有効; codex diagnose notes/reviews/2026-09-27-transport-tables-diagnose.md の事前固定条件):
  セル・ghost・壁の格納値 (S) は表引きの float 経路になる。ULP 条件は外し (新 float 経路に旧 double の基準を課さない)、
  独立 double 参照と ≤1e-5 で判定する。(D) の double 評価は範囲外の委譲先として残るので同じ基準で照合を続ける。
  (T1) 単体: 各表 (実種の μ・λ、CE・CEA 相互作用の組の η、剛体球の組は種別表の μ から実行時) を FORGE_TRANSPORT_TABLE_PROBE で
       GPU 評価し、実際の float T を独立参照へ渡して相対 ≤2e-6。検査点は各小区間の 17 等分点 (ln T で等分、float に丸め)、
       全境界 (参照側で独立に列挙: フィット区間・H2O 接続点・T* クランプ点 0.3ε/100ε・NASA Tlo/Tmid/Thi) の直前・直後 8 個の float、
       表範囲外。表の分割区間の端が参照側の境界をすべて含むことも確かめる。NaN/Inf・非正・未更新 0 件。
  (T2) 混合: 同じ Y・T を float に丸めて表引き経路へ渡した μ・λ (状態表 (6)) が独立 double 参照と ≤1e-5。
       T は全境界 ±8 float・表範囲外・小区間の 17 等分点 (小区間は間引き: 各分割区間の最初と最後は必ず、他は stride ごと)。
  (F)  全 CEA 指定の N2–H2O 16 状態 (T 400/600/1000/2000 K × X_H2O 0/0.1/0.5/1) を表引き経路で FCEA2 (frozen) と ≤0.1 %。
  (AB) 区間選択の判別 (単成分 fit、1000 K 以下 μ 1e-5・超過 1.01e-5 Pa s、λ 0.02/0.0202): 1000 K と両隣の float で、
       A = float の ln T で選ぶ (host ハーネス transport_table_ab_host.cpp にだけある) は不合格、
       B = 元の T で選ぶ (ソルバ; host 同関数と GPU 0 step の両方) は合格が期待。
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
from transport_reference import Reference, TR, set_lj_source  # noqa: E402

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
        if m == "kinetic":   # NASA-9 の全区間境界 (transport_reference の bounds; 段 3 #13-3 から 3 区間の種は 6000/20000 K も)
            E |= set(sp["bounds"]) | {0.3 * sp["eps"], 100.0 * sp["eps"]}
    return sorted(E | set(pair_edges(ref).get("all", [])))


def pair_edges(ref):
    """組ごとの境界 (参照側で独立に列挙)。戻り値 {(a, b): [...], "all": [...]}。剛体球は境界なし。"""
    out, allE = {}, set()
    n = len(ref.reals)
    for a in range(n):
        for b in range(a + 1, n):
            A, B = ref.reals[a], ref.reals[b]
            E = set()
            if A["model"] == "kinetic" and B["model"] == "kinetic":
                eps = math.sqrt(A["sp"]["eps"] * B["sp"]["eps"])
                E |= {0.3 * eps, 100.0 * eps}
                out[(a, b)] = sorted(E)
                allE |= E
                continue
            if A["trans"] and B["trans"]:
                for k in ((A["trans"], B["trans"]), (B["trans"], A["trans"])):
                    if k in TR and TR[k]["V"]:
                        for lo, hi, _ in TR[k]["V"]:
                            E |= {lo, hi}
                        break
            out[(a, b)] = sorted(E)
            allE |= E
    out["all"] = sorted(allE)
    return out


def species_edges(ref, r):
    """実種 r の境界 (参照側で独立に列挙)。"""
    e = ref.reals[r]
    m, sp = e["model"], e["sp"]
    E = set()
    if m in ("cea", "custom:h2o_iapws_cea_v1"):
        for kind in ("V", "C"):
            for lo, hi, _ in TR[(e["trans"], "")][kind]:
                E |= {lo, hi}
    if m == "custom:h2o_iapws_cea_v1":
        E |= {253.15, 500.0, 700.0}
    if m == "fit":
        for kind in ("V", "C"):
            for row in sp["fit"][kind]:
                E |= {float(row[0]), float(row[1])}
    if m == "kinetic":
        E |= set(sp["bounds"]) | {0.3 * sp["eps"], 100.0 * sp["eps"]}
    return sorted(E)


OUT_OF_RANGE = [60.0, 120.0, 149.9, 15000.5, 16000.0, 20000.0]   # 表の範囲 [150, 15000] K の外 (double 評価へ委譲)
NEAR = 8   # 境界の直前・直後に取る float の個数
T_LO, T_HI = 50.0, 30000.0   # 検査する T の範囲 (表の範囲外を含む)


def around(b, k=NEAR):
    """b を float に丸めた値と、その直前・直後 k 個の float。"""
    fb = np.float32(b)
    out = [float(fb)]
    lo = hi = fb
    for _ in range(k):
        lo = np.nextafter(lo, np.float32(-np.inf))
        hi = np.nextafter(hi, np.float32(np.inf))
        out += [float(lo), float(hi)]
    return out


def t_list(ref):
    T = {f32(t) for t in BASE_T} | {f32(t) for t in OUT_OF_RANGE} | {f32(150.0), f32(15000.0)}
    for b in fit_edges(ref) + [150.0, 15000.0]:
        if b > 0:
            T |= set(around(b))
    # 50 K 未満は物理的に使わない (T* クランプ点 0.3ε は数 K になりうるが、CEA フィットの C/T² が参照側で桁あふれする)
    return sorted(t for t in T if T_LO <= t <= T_HI)


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
        ph.pop("ljSource", None)
        if getattr(self.a, "lj_source", None):
            ph["ljSource"] = list(self.a.lj_source)
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
        e["FORGE_ALLOW_UNVERIFIED_SPECIES"] = "1"   # seed の場は species 属性なし (#3c 以降は既定で停止)
        if self.a.blocksize:
            e["FORGE_CUDA_BLOCKSIZE"] = str(self.a.blocksize)
        return e

    def probe(self, d, states, nS, nocells=False, name="states.txt"):
        p = os.path.join(d, name)
        with open(p, "w") as f:
            f.write(f"{len(states)} {nS}\n")
            for T, ro, Y in states:
                f.write(" ".join("%.17g" % v for v in [T, ro] + list(Y)) + "\n")
        e = self.env()
        e["FORGE_TRANSPORT_PROBE"] = p
        if nocells:
            e["FORGE_TRANSPORT_PROBE_NOCELLS"] = "1"
        r = subprocess.run([self.a.forge], cwd=d, env=e, capture_output=True, text=True, timeout=1800)
        open(os.path.join(d, "forge_probe.log"), "w").write(r.stdout + r.stderr)
        if r.returncode != 0:
            return None, (r.stdout + r.stderr)[-3000:]
        return read_probe(p), ""

    def table_probe(self, d, points, name="points.txt"):
        """FORGE_TRANSPORT_TABLE_PROBE: points = [(kind, idx, T)]。戻り値 (layout, v0, v1)。"""
        p = os.path.join(d, name)
        with open(p, "w") as f:
            f.write(f"{len(points)}\n")
            for k, i, T in points:
                f.write(f"{k} {i} {T:.9g}\n")
        e = self.env()
        e["FORGE_TRANSPORT_TABLE_PROBE"] = p
        r = subprocess.run([self.a.forge], cwd=d, env=e, capture_output=True, text=True, timeout=1800)
        open(os.path.join(d, "forge_table_probe.log"), "a").write(r.stdout + r.stderr)
        if r.returncode != 0:
            return None, None, None, (r.stdout + r.stderr)[-3000:]
        lay = json.load(open(p + ".json"))
        buf = np.fromfile(p + ".bin", dtype=np.float32)
        N = lay["N"]
        assert buf.size == 2 * N, (buf.size, N)
        return lay, buf[:N], buf[N:], ""


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
        if h.get("hasLiq", 0):   # 凝縮 carrier の気相組成の試験 (test_transport_gas_phase.py): 液・P・分子拡散係数
            q["rog"] = take(np.float32, nA)
            q["P"] = take(np.float32, nA)
            q["D"] = take(np.float32, nA * nS).reshape(nS, nA)
        passes.append(q)
    s = {"mu": take(np.float64, K), "lam": take(np.float64, K), "X": take(np.float64, K * nR).reshape(K, nR)}
    if h.get("table", 0):
        s["tmu"] = take(np.float32, K)
        s["tlam"] = take(np.float32, K)
    assert off[0] == len(buf), (off[0], len(buf))
    return {"h": h, "passes": passes, "states": s}


# ------------------------------------------------------------------ 1 構成の検査
class Result:
    pass


REFUSED_NO_LJ = object()   # run_config の戻り値: LJ が無く起動拒否を確かめた構成 (#14-L1; 比較は SKIP)


def refused_without_lj(R, tag, species, transport, db, ref):
    """--lj-source の集合に LJ の無い内蔵種がある構成で、LJ を読む使い方になるもの: forge が起動時に拒否することを確かめて True を返す
    (参照側も LJ なし。#14-L1)。R.make の既定 viscMethod 2 では化学種 2 以上なら LJ の混合平均拡散 (speciesDiffusionMethod 1) が
    全実種の LJ を読む (先に検査される)、単成分なら kinetic の種だけ。該当しなければ False (通常どおり試験する)。"""
    nolj = [e["name"] for e in ref.reals if e["sp"]["sigma"] is None]
    nolj_kin = [e["name"] for e in ref.reals if e["model"] == "kinetic" and e["sp"]["sigma"] is None]
    if len(species) >= 2 and nolj:
        needle, what = "have no Lennard-Jones data in any of the LJ sets searched", f"species {nolj} (mixture-averaged diffusion)"
    elif nolj_kin:
        needle, what = f"'{nolj_kin[0]}' has no Lennard-Jones data", f"kinetic species {nolj_kin}"
    else:
        return False
    d = R.make(tag + "_nolj", species, transport, db)
    r = subprocess.run([R.a.forge, "--resolve-species"], cwd=d, env=R.env(), capture_output=True, text=True, timeout=600)
    msg = r.stderr + r.stdout
    check(r.returncode != 0 and needle in msg,
          f"{tag}: {what} have no LJ in ljSource {list(R.a.lj_source or [])}: forge refuses at startup (rc={r.returncode})")
    return True


def run_config(R, tag, species, transport, comps, db=None, ro=0.37, T_override=None):
    """comps: 輸送種のモル分率の組の列。戻り値: 状態ごとの double (状態表経由) とセル経路 (float roY 経由) の値。"""
    ref = Reference(species, transport, db)
    if refused_without_lj(R, tag, species, transport, db, ref):
        return REFUSED_NO_LJ
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
    tab = bool(h.get("table", 0))
    check(max(worst["fcell_mu"], worst["fcell_lam"], worst["fwall_mu"], worst["fwall_lam"]) <= 1e-5 and (tab or ulp_worst <= 2),
          f"{tag}: S float storage ({'table' if tab else 'double'} path) vs independent double reference: vis_lam {worst['fcell_mu']:.2e} "
          f"thermCond {worst['fcell_lam']:.2e} (cells+ghosts), wall mu_w {worst['fwall_mu']:.2e} lam_w {worst['fwall_lam']:.2e} (<= 1e-5); "
          f"max ULP vs float(ref) {ulp_worst}" + (" (recorded; the table path is judged by 1e-5 only)" if tab else " (<= 2)"))
    if tab:
        tm, tl = out["states"]["tmu"], out["states"]["tlam"]
        t_em = t_el = 0.0
        tbad = 0
        for k, (T, _, Y) in enumerate(states):
            # 表引き経路の実際の入力: T と Y を float に丸めたもの
            r = ref.state_Y(f32(T), [f32(y) for y in Y])
            t_em = max(t_em, float(rel(tm[k], r["mu"])))
            t_el = max(t_el, float(rel(tl[k], r["lam"])))
            tbad += int(not (np.isfinite(tm[k]) and tm[k] > 0 and np.isfinite(tl[k]) and tl[k] > 0))
        check(t_em <= 1e-5 and t_el <= 1e-5 and tbad == 0,
              f"{tag}: T2 table path, float-rounded state inputs ({K} states): mu {t_em:.2e} lam {t_el:.2e} (<= 1e-5); "
              f"NaN/Inf/non-positive/non-updated {tbad}")
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
    if A is REFUSED_NO_LJ and B is REFUSED_NO_LJ:
        print(f"[SKIP] AB {tagA} vs {tagB}: both configurations refused for missing LJ under the --lj-source (checked above)")
        return
    if A is None or B is None or A is REFUSED_NO_LJ or B is REFUSED_NO_LJ:
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


# ------------------------------------------------------------------ 表引き (#5t2-3)
def fit_breaks(ivs):
    """区間の列 [(lo, hi, ...)] の内側の継ぎ目 (最後の区間の hi は外挿なので含めない)。"""
    his = sorted(float(v[1]) for v in ivs)
    return his[:-1]


def species_breaks(ref, r):
    """実種 r の式が切り替わる温度 (参照側で独立に列挙; 表の分割区間の端が含むべきもの)。"""
    e = ref.reals[r]
    m, sp = e["model"], e["sp"]
    B = set()
    if m in ("cea", "custom:h2o_iapws_cea_v1"):
        for kind in ("V", "C"):
            B |= set(fit_breaks(TR[(e["trans"], "")][kind]))
    if m == "custom:h2o_iapws_cea_v1":
        B |= {253.15, 500.0, 700.0}
    if m == "fit":
        for kind in ("V", "C"):
            B |= set(fit_breaks(sp["fit"][kind]))
    if m == "kinetic":
        B |= set(sp["bounds"]) | {0.3 * sp["eps"], 100.0 * sp["eps"]}
    return sorted(B)


def pair_breaks(ref, a, b):
    A, B = ref.reals[a], ref.reals[b]
    kind = ref.eta_pair(a, b, 1000.0, [1e-5] * len(ref.reals))[1]
    if kind == "ce":
        eps = math.sqrt(A["sp"]["eps"] * B["sp"]["eps"])
        return kind, [0.3 * eps, 100.0 * eps]
    if kind == "cea":
        for k in ((A["trans"], B["trans"]), (B["trans"], A["trans"])):
            if k in TR and TR[k]["V"]:
                return kind, fit_breaks(TR[k]["V"])
    return kind, []


def grid_points(segs, stride=1):
    """分割区間 [[Ta, Tb, m, Tupper], ...] の小区間ごとの 17 等分点 (ln T で等分、float に丸め)。
    stride > 1 のときは各分割区間の最初・最後の小区間と stride ごとの小区間だけ。"""
    T = []
    for Ta, Tb, m, _ in segs:
        la, lb = math.log(Ta), math.log(Tb)
        h = (lb - la) / m
        js = range(m) if stride <= 1 else sorted({0, m - 1} | set(range(0, m, stride)))
        for j in js:
            for q in range(17):
                T.append(f32(math.exp(la + (j + q / 16.0) * h)))
    return T


PAIR_KIND = {1: "ce", 2: "cea", 3: "rigid"}


def table_singles(R, tag, species, transport, db=None):
    """T1: 各表の単体値 (GPU, 表引き) を独立参照と比べる。戻り値: 表の配置 (layout) と参照。"""
    ref = Reference(species, transport, db)
    if refused_without_lj(R, tag + "_tab", species, transport, db, ref):
        return None, ref
    d = R.make(tag + "_tab", species, transport, db)
    lay, _, _, err = R.table_probe(d, [], name="layout.txt")
    if lay is None:
        check(False, f"{tag}: T1 table probe failed: {err}")
        return None, ref
    n = lay["nReal"]
    tabs = lay["tables"]
    npair = n * (n - 1) // 2
    ok_layout = (n == len(ref.reals) and len(tabs) == n + npair)
    # 分割区間の端が参照側の境界をすべて含むか (150 < b < 15000)
    missing = []
    kind_mismatch = []
    ends_of = {}
    for t in tabs:
        ends = set()
        for Ta, Tb, _, _ in t["segs"]:
            ends |= {Ta, Tb}
        ends_of[(t["kind"], t["idx"])] = ends
        if t["kind"] == 0:
            br = species_breaks(ref, t["idx"])
        else:
            rk, br = pair_breaks(ref, t["a"], t["b"])
            if PAIR_KIND[t["pair_kind"]] != rk:
                kind_mismatch.append((t["a"], t["b"], PAIR_KIND[t["pair_kind"]], rk))
            if rk == "rigid":
                br = []
        for b in br:
            if 150.0 < b < 15000.0 and not any(abs(b - e) <= 1e-12 * b for e in ends):
                missing.append((t["kind"], t["idx"], b))
    check(ok_layout and not missing and not kind_mismatch,
          f"{tag}: T1 layout: {n} species tables + {sum(1 for t in tabs if t['kind'] == 1 and t['segs'])} tabulated pairs "
          f"({sum(1 for t in tabs if t['kind'] == 1 and not t['segs'])} rigid-sphere), T in [{lay['Tmin']:.9g}, {lay['Tmax']:.9g}], "
          f"{lay['bytes']} bytes; independent breakpoints missing from segment ends {missing[:4]} (0), pair kind mismatch {kind_mismatch[:3]} (0)")
    points = []
    npts = {"grid": 0, "edge": 0, "out": 0}
    for t in tabs:
        k, i = t["kind"], t["idx"]
        if k == 0:
            segs, br = t["segs"], species_breaks(ref, i)
        else:
            rk, br = pair_breaks(ref, t["a"], t["b"])
            segs = t["segs"] if t["segs"] else tabs[t["a"]]["segs"]   # 剛体球は種 a の格子で検査
            if rk == "rigid":
                br = species_breaks(ref, t["a"]) + species_breaks(ref, t["b"])
        g = grid_points(segs)
        e = [x for b in br + [150.0, 15000.0] if b > 0 for x in around(b) if T_LO <= x <= T_HI]
        o = [f32(x) for x in OUT_OF_RANGE]
        npts["grid"] += len(g)
        npts["edge"] += len(e)
        npts["out"] += len(o)
        points += [(k, i, T) for T in g + e + o]
    lay2, v0, v1, err = R.table_probe(d, points)
    if lay2 is None:
        check(False, f"{tag}: T1 table probe failed: {err}")
        return lay, ref
    worst = {"mu": 0.0, "lam": 0.0, "eta": 0.0, "eta_rigid": 0.0}
    where = {}
    bad = 0
    cache = {}

    def sp_ref(r, T):
        key = (r, T)
        if key not in cache:
            cache[key] = ref.species(r, T)
        return cache[key]

    for (k, i, T), a, b in zip(points, v0, v1):
        if k == 0:
            m, l = sp_ref(i, T)
            bad += int(not (np.isfinite(a) and a > 0 and np.isfinite(b) and b > 0))
            for nm, got, want in (("mu", a, m), ("lam", b, l)):
                e = float(rel(got, want))
                if e > worst[nm]:
                    worst[nm], where[nm] = e, (i, T)
        else:
            t = tabs[n + i]
            eta = [0.0] * n
            eta[t["a"]] = sp_ref(t["a"], T)[0]
            eta[t["b"]] = sp_ref(t["b"], T)[0]
            want, rk = ref.eta_pair(t["a"], t["b"], T, eta)
            bad += int(not (np.isfinite(a) and a > 0))
            nm = "eta_rigid" if rk == "rigid" else "eta"
            e = float(rel(a, want))
            if e > worst[nm]:
                worst[nm], where[nm] = e, (t["a"], t["b"], T)
    check(max(worst.values()) <= 2e-6 and bad == 0,
          f"{tag}: T1 single values vs independent reference ({len(points)} points: 17/sub-interval {npts['grid']}, "
          f"boundaries +-{NEAR} floats {npts['edge']}, out of range {npts['out']}): mu {worst['mu']:.2e} lam {worst['lam']:.2e} "
          f"eta(CE/CEA) {worst['eta']:.2e} eta(rigid, from table mu) {worst['eta_rigid']:.2e} (<= 2e-6) at {where}; "
          f"NaN/Inf/non-positive/non-updated {bad}")
    return lay, ref


def table_mix_grid(R, tag, species, transport, comps, lay, db=None, stride=1, ro=0.37):
    """T2: 小区間の 17 等分点 (間引き stride)・境界 ±8 float・範囲外で、表引き経路の混合値を独立参照と比べる。"""
    ref = Reference(species, transport, db)
    Ts = set(t_list(ref))
    for t in lay["tables"]:
        if t["kind"] == 0:
            Ts |= set(grid_points(t["segs"], stride))
    Ts = sorted(Ts)
    states = []
    for X in comps:
        w = [x * m for x, m in zip(X, ref.mw)]
        Y = [v / sum(w) for v in w]
        states += [(T, ro, Y) for T in Ts]
    d = R.make(tag + "_mixgrid", species, transport, db)
    out, err = R.probe(d, states, len(species), nocells=True)
    if out is None:
        check(False, f"{tag}: T2 grid probe failed: {err}")
        return
    tm, tl = out["states"]["tmu"], out["states"]["tlam"]
    em = el = 0.0
    wm = None
    bad = 0
    for k, (T, _, Y) in enumerate(states):
        r = ref.state_Y(f32(T), [f32(y) for y in Y])
        e1, e2 = float(rel(tm[k], r["mu"])), float(rel(tl[k], r["lam"]))
        if max(e1, e2) > max(em, el):
            wm = (T, [round(y, 4) for y in Y[:4]])
        em, el = max(em, e1), max(el, e2)
        bad += int(not (np.isfinite(tm[k]) and tm[k] > 0 and np.isfinite(tl[k]) and tl[k] > 0))
    check(em <= 1e-5 and el <= 1e-5 and bad == 0,
          f"{tag}: T2 mixture on the grid ({len(comps)} compositions x {len(Ts)} T = {len(states)} states; sub-interval stride {stride}): "
          f"mu {em:.2e} lam {el:.2e} (<= 1e-5) worst at {wm}; NaN/Inf/non-positive/non-updated {bad}")


def fcea_via_table(R, work):
    """F: 全 CEA 指定の N2–H2O 16 状態を表引き経路で FCEA2 (frozen) と比べる。"""
    try:
        from test_species_transport import run_fcea, CEA
    except Exception as ex:  # noqa: BLE001
        check(False, f"F cannot import run_fcea: {ex}")
        return
    if not os.path.exists(os.path.join(CEA, "FCEA2")):
        check(False, "F FCEA2 not found at .venv-cea/nasa_cea/FCEA2")
        return
    fc = run_fcea(os.path.join(work, "fcea_table"))
    ref = Reference(["N2", "H2O"], {"N2": "cea", "H2O": "cea"})
    states, rows = [], []
    for xw, pts in fc.items():
        X = [1.0 - xw, xw]
        w = [x * m for x, m in zip(X, ref.mw)]
        Y = [v / sum(w) for v in w]
        for T, v, c in pts:
            states.append((T, 1.0, Y))
            rows.append((xw, T, v * 1e-4, c * 1e-1))
    d = R.make("fcea_table", ["N2", "H2O"], {"N2": "cea", "H2O": "cea"})
    out, err = R.probe(d, states, 2, nocells=True)
    if out is None:
        check(False, f"F probe failed: {err}")
        return
    wm = wl = 0.0
    for k, (xw, T, mu, lam) in enumerate(rows):
        em, el = out["states"]["tmu"][k] / mu - 1.0, out["states"]["tlam"][k] / lam - 1.0
        wm, wl = max(wm, abs(em)), max(wl, abs(el))
        print(f"       X_H2O {xw:4.2f} T {T:6.0f}  mu {100 * em:+.4f} %  lam {100 * el:+.4f} %  (table path)")
    check(len(rows) == 16 and wm <= 1e-3 and wl <= 1e-3,
          f"F all-CEA N2-H2O 16 states, table path vs FCEA2: max |mu| {100 * wm:.4f} %, max |lam| {100 * wl:.4f} % (<= 0.1 %)")


STEP_DB = None


def step_db():
    n2 = n2_thermo()
    ln = math.log
    return {"XSTEP": {"MW": 0.028, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
                      "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"],
                      "LJ_sigma": 3.6, "LJ_eps_kB": 97.0,
                      "transport_fit": {"V": [[150.0, 1000.0, 0.0, 0.0, 0.0, ln(100.0)], [1000.0, 20000.0, 0.0, 0.0, 0.0, ln(101.0)]],
                                        "C": [[150.0, 1000.0, 0.0, 0.0, 0.0, ln(200.0)], [1000.0, 20000.0, 0.0, 0.0, 0.0, ln(202.0)]]}}}


def ab_select(R, work):
    """AB: 区間選択だけを変える判別 (A = float ln T / B = 元の T)。"""
    db = step_db()
    ref = Reference(["XSTEP"], {"XSTEP": "fit"}, db)
    exe = os.path.join(work, "transport_table_ab_host")
    r = subprocess.run(["g++", "-O1", "-std=c++17", "-Wno-unknown-pragmas", "-I", SOLVER,
                        os.path.join(HERE, "transport_table_ab_host.cpp"), "-o", exe], capture_output=True, text=True)
    if r.returncode != 0:
        check(False, f"AB host harness build failed: {r.stderr[-2000:]}")
        return
    o = json.loads(subprocess.run([exe], capture_output=True, text=True).stdout)
    eA = eB = 0.0
    for p in o["points"]:
        T = p["T"]
        m, l = ref.species(0, T)
        a = max(float(rel(p["muA"], m)), float(rel(p["lamA"], l)))
        b = max(float(rel(p["muB"], m)), float(rel(p["lamB"], l)))
        eA, eB = max(eA, a), max(eB, b)
        print(f"       T {T:.9g} K (float ln T {p['lnT_float']:.9g}): ref mu {m:.6e}; A mu {p['muA']:.6e} ({a:.2e}), B mu {p['muB']:.6e} ({b:.2e})")
    # GPU (0 step): ソルバの表引き (B) を同じ点で
    d = R.make("ab_step_fit", ["XSTEP"], {"XSTEP": "fit"}, db)
    Ts = [p["T"] for p in o["points"]]
    lay, v0, v1, err = R.table_probe(d, [(0, 0, T) for T in Ts])
    eG = float("inf")
    if lay is not None:
        eG = max(max(float(rel(v0[k], ref.species(0, T)[0])), float(rel(v1[k], ref.species(0, T)[1]))) for k, T in enumerate(Ts))
    check(eA > 2e-6 and eB <= 2e-6 and eG <= 2e-6,
          f"AB interval selection (single fit, 1 % step at 1000 K; 1000 K and neighbouring floats): "
          f"A (float ln T) max {eA:.2e} -> {'FAIL as expected' if eA > 2e-6 else 'passes (hypothesis refuted)'}; "
          f"B (original T) host {eB:.2e}, GPU 0-step {eG:.2e} (<= 2e-6)")


# ------------------------------------------------------------------ (B) 既定経路のビット一致
def legacy_visc2_refused(R, a):
    """(B2) viscMethod 2 + transport なしは起動時エラー (旧 kinetic 経路は撤去; plan §4.3c 案 C)。"""
    for tag, species in (("builtin N2/H2O/O2/AR/CO2", ["N2", "H2O", "O2", "AR", "CO2"]), ("builtin N2 single", ["N2"])):
        d = R.make("visc2_notransport_" + tag.split(",")[0].replace(" ", "_").replace("/", "-"), species, None, visc=2, nstep=0)
        r = subprocess.run([a.forge], cwd=d, env=R.env(), capture_output=True, text=True, timeout=600)
        out = r.stdout + r.stderr
        open(os.path.join(d, "forge.log"), "w").write(out)
        check(r.returncode != 0 and "viscMethod: 2 requires physProp.transport" in out and "viscMethod: 1" in out
              and not os.path.exists(os.path.join(d, "res_0.h5")),
              f"B2 {tag}, viscMethod 2 without physProp.transport: refused at startup (rc={r.returncode})")


def bit_identity(R, a):
    # viscMethod 2 (transport なし) は起動時エラーになったので比較から外し (B2)、viscMethod 1 を 2 構成と viscMethod 0 を比べる
    cases = [
        ("seed MIXDRY/H2O ext DB, viscMethod 1", dict(species=["MIXDRY", "H2O"], transport=None, keep_db_of_seed=True, visc=1)),
        ("builtin N2/H2O/O2/AR/CO2, viscMethod 1", dict(species=["N2", "H2O", "O2", "AR", "CO2"], transport=None, visc=1)),
        # 段 3 (#13-3) で熱物性が不変な内蔵種だけの構成 (N2/O2/CO2 は 200–6000 K で Δ 0; 第 3 区間が増えただけ)
        # 互換性ハッシュは第 3 区間の分だけ変わるが、値 (T・μ・λ) はビット一致を求める
        ("builtin N2/O2/CO2, viscMethod 1", dict(species=["N2", "O2", "CO2"], transport=None, visc=1, values_unchanged=True)),
        ("seed MIXDRY/H2O ext DB, viscMethod 0", dict(species=["MIXDRY", "H2O"], transport=None, keep_db_of_seed=True, visc=0)),
    ]
    for tag, kw in cases:
        outs, hashes = [], []
        for which, binp in (("new", a.forge), ("old", a.base_forge)):
            d = R.make("bit_" + which + "_" + tag.split(",")[0].replace(" ", "_").replace("/", "-"), kw["species"], None,
                       visc=kw["visc"], keep_db_of_seed=kw.get("keep_db_of_seed", False), nstep=0)
            rr = subprocess.run([binp, "--resolve-species"], cwd=d, env=R.env(), capture_output=True, text=True, timeout=600)
            hashes.append((rr.stdout.strip().splitlines() or [""])[-1])
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
        if hashes[0] != hashes[1] and not kw.get("values_unchanged", False):
            # 旧バイナリと熱物性データが違う (段 3 #13-3 の前後: 内蔵 H2O の MW・AR の第 2 区間)。同じ保存量から温度反転した T が変わるので
            # ビット一致は求めず、差を記録する (経路の同一性は熱物性が不変な構成でビット一致を見る)
            drel = {k: float(np.max(np.abs(n[k].astype(np.float64) - o[k].astype(np.float64)) / np.abs(o[k].astype(np.float64))))
                    for k in ["T"] + outk} if have else {}
            check(have and set(n) == set(o),
                  f"B {tag}: species data differ from the base binary (compat {hashes[1][:16]} -> {hashes[0][:16]}, #13-3); "
                  f"no bit identity required — max rel diff {', '.join(f'{k} {v:.2e}' for k, v in drel.items())}; "
                  f"differing bytes {nd} over {len(n['T'])} CVs")
            continue
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
    ap.add_argument("--skip-table", action="store_true", help="表引きの試験 (T1/T2/F/AB) を省く")
    ap.add_argument("--stride12", type=int, default=4, help="T2 の実種 12 で 17 等分点を取る小区間の間引き")
    ap.add_argument("--stride32", type=int, default=32, help="T2 の実種 32 で 17 等分点を取る小区間の間引き")
    ap.add_argument("--lj-source", default=None, help="内蔵種の LJ の集合リスト (カンマ区切り; 省略はソルバの既定で config に書かない)")
    a = ap.parse_args()
    if a.lj_source:
        a.lj_source = tuple(s for s in a.lj_source.split(",") if s)
        set_lj_source(a.lj_source)
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
    print(f"work dir: {root}\nseed: {a.seed_run}\nforge: {a.forge}\nljSource: {list(a.lj_source) if a.lj_source else 'default'}", flush=True)
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
    if B is not None and B is not REFUSED_NO_LJ:
        k = B.keys.index((0, f32(400.0)))
        xhe = dict(zip(B.real, B.sX[k]))["He"]
        cA = A.cell.get((0, f32(400.0))) if A else None
        cB = B.cell.get((0, f32(400.0)))
        check(abs(xhe - 0.6) <= 1e-12,
              f"AB B at 400 K: expanded X_He {xhe:.15f} (mole basis 0.6; mass-basis mix-up would give 0.887308)")
        if A is not None and A is not REFUSED_NO_LJ and cA and cB:
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
    if r12 and r12 is not REFUSED_NO_LJ:
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
    if r32 and r32 is not REFUSED_NO_LJ:
        check(len(r32.real) == 32, f"real32: real species {len(r32.real)} == 32 (TRANSPORT_MAX_REAL_SPECIES)")

    # ---- 表引き (#5t2-3) ----
    if not a.skip_table:
        ab_select(R, root)
        fcea_via_table(R, root)
        layA, _ = table_singles(R, "A_full_N2_He", ["N2", "He"], trNHe)
        if layA:
            table_mix_grid(R, "A_full_N2_He", ["N2", "He"], trNHe, compsA, layA, ro=1.0)
        table_singles(R, "single_N2_cea", ["N2"], {"N2": "cea"})
        layH, _ = table_singles(R, "single_H2O_custom", ["H2O"], {"H2O": "custom:h2o_iapws_cea_v1"})
        if layH:
            table_mix_grid(R, "single_H2O_custom", ["H2O"], {"H2O": "custom:h2o_iapws_cea_v1"}, [[1.0]], layH)
        lay12, _ = table_singles(R, "real12", sp12, tr12, db12)
        if lay12:
            table_mix_grid(R, "real12", sp12, tr12, comps12, lay12, db12, stride=a.stride12)
        lay32, _ = table_singles(R, "real32", sp32, tr32, db32)
        if lay32:
            table_mix_grid(R, "real32", sp32, tr32, [c1, c2], lay32, db32, stride=a.stride32)

    # ---- viscMethod 2 + transport なしは起動時エラー ----
    legacy_visc2_refused(R, a)

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
