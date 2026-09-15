#!/usr/bin/env python3
"""化学種の集合/順序が違う run へ場を移す (種変換 restart)。同一メッシュ・index コピー。

旧 `[MIXDRY, H2O]` (擬似種) の収束場を新順序 (例 `[H2O, N2, O2, AR, CO2]` の full) へ移す、SERN の `[EXH, AIR]` (流れ lump) と
full (11 種 + トレーサ `roXi`) を相互に移す、といった「輸送種の配置が違う restart」を `species_meta.yaml` の情報で行う。
plans/active/thermophysics-cea-mole-fraction-species.md §2 (forge 本体) / §4.5 M4、codex 2026-09-16 result M2/M4。

  convert_species_field.py SRC_res.h5 DST_input.h5 --meta DST/species_meta.yaml [--src-meta SRC/species_meta.yaml]
      [--mode conserve|reinit] [--src-run SRC_DIR] [--dst-run DST_DIR] [--reconstruct-roe auto|always|never]
      [--T-tol 0.05] [--drop-moments] [--dry-run]

2 つの明示的な操作 (`--mode`):
  conserve (既定) — **保存的な展開/縮約**。source の各輸送種を `expansion` (lump 内質量分率) で実種に展開し実種ごとに合算、
      destination の輸送種へ「実種の行き先が一意」なときだけ縮約する (同名・同 expansion の輸送種は 1:1 コピー)。
      実種が複数の destination lump に属し得る (SERN の EXH/AIR はともに N2 を含む) 場合は拒否し、`--mode reinit` を要求する。
  reinit — **流れによる再初期化**。source の流入元分率 ξ (source の `roXi` があればそれ、無ければ source meta で純流入ラベルに
      なっている輸送種、例 Y_EXH) を使い、実種ごとの質量を destination の流れ lump へ `streams[...].Y_transport` の重み
      w_j(r) = (ξ Yt_in[j] + (1−ξ) Yt_ext[j]) · expansion_j[r] で分配する (ΣY=1)。**lump の内部組成は固定なので、source の
      組成が二流線形部分空間 ξY_in+(1−ξ)Y_ext から外れている分 (差動拡散・full/lumped の輸送差) は実種ごとの質量に射影損失が
      出る = 情報を落とす操作** (plan §4.5 M4)。損失は実種ごとに報告し `--reinit-tol` (既定 1e-3·max ρ) を超えたら失敗。
  どちらでも destination が `tracer.enabled` で source に `roXi` が無ければ **roXi = ρ ξ を生成**する (ξ が導けなければ拒否)。

エネルギーと温度:
  - source が res (T あり) なら T はソルバの値。input h5 (T なし) なら **ソルバと同じ二相 EOS** で `roe` から反転する:
    e = e_gas(Y_total, T) + g (R_w T − L(T)) (carrier 形, condensationEOS_d.cuh `cond_T_from_e_carrier`)、pure TP は
    e = e_v(T) + g R_mix T − g L(T) (`cond_T_from_e_onetemp`)。L(T) は condensationProperties_d.cuh の `h2o_latent` / `n2_latent` を移植。
  - DB (`species_db.yaml`) / datum (`thermoHrefTemp`) / 種集合が変わるときは `roe += ρ [e_gas,dst(Y_dst,T) − e_gas,src(Y_src,T)]`
    (差分形; 液相項 g(R_w T−L) は不変なので湿潤セルでも正しい)。
検査 (**1 つでも破れば書き込まず失敗終了**): 各セル ΣY_dst=1 (1e-6)、実種ごとの ρY 保存 (conserve: 1e-6·max ρ / reinit: 射影損失を
  `--reinit-tol`)、総水量 ρ·Y_w (Y_w は液相込みの総水分率なので ρ(Y_w+g) ではない) の保存 (conserve: rel 1e-9 / reinit: `--reinit-tol` 相対)、destination DB + 二相 EOS で `roe` を反転した T と source T の差
  (乾き・湿潤の全セル, `--T-tol` 既定 0.05 K)、roXi/ρ ∈ [0,1]。

- SRC: res_*.h5 (原始量 P,T,Ux,.. + Y{s}) か input h5 (保存量 roY{s})。DST: 同一メッシュ・同一 CV 数の input h5。
  ro/roU/roe/roK/roOmega・凝縮モーメント `rog_*/roQ*_*` (凝縮種が同名のとき) も index コピーする。
- 種名は `species_meta.yaml` (`species`, `expansion`, `streams`, `condensing_species`, `tracer`, `MW`) を正とし、無ければ run dir の
  `solverConfig.yaml` から取る (`--src-run/--dst-run` 省略時は h5 の隣)。condModel / condGasSpecies / thermoHrefTemp / species_db.yaml は
  run dir の solverConfig.yaml から読む。
"""
import argparse, os, sys
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from total_quantities import _TPGas  # noqa: E402
from forge_species import species_info, load_yaml_str  # noqa: E402


def _up(s):
    return str(s).upper()


# ----------------------------------------------------------------------------- 二相 EOS (ソルバ移植, numpy)
COND_T_PROP_FLOOR = 45.0
COND_N2_LATENT_TA = 70.0
COND_N2_CPV = 1038.8


def _h2o_nasa9_h_mass(a, T):
    Rw = 8.314462618 / 0.0180153
    hRT = (-a[0] / (T * T) + a[1] * np.log(T) / T + a[2] + a[3] * T / 2.0 + a[4] * T * T / 3.0
           + a[5] * T ** 3 / 4.0 + a[6] * T ** 4 / 5.0 + a[7] / T)
    return hRT * Rw * T


def h2o_latent(T):
    """condensationProperties_d.cuh h2o_latent の移植 (CEA 気相 H2O − 液相 H2O(L), 273.15 K 未満は cp_l 一定外挿, [1.5e6, 3.5e6] クランプ)。"""
    ag = np.array([-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,
                   -7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04])
    al = np.array([1.326371304e+09, -2.448295388e+07, 1.879428776e+05, -7.678995050e+02,
                   1.761556813e+00, -2.151167128e-03, 1.092570813e-06, 1.101760476e+08])
    Tf = 273.15
    Tg = np.clip(np.asarray(T, float), COND_T_PROP_FLOOR, 1000.0)
    hv = _h2o_nasa9_h_mass(ag, Tg)
    Tl = np.minimum(Tg, 373.15)
    hl_hi = _h2o_nasa9_h_mass(al, np.maximum(Tl, Tf))
    h0 = _h2o_nasa9_h_mass(al, np.array([Tf]))[0]
    cpl = (_h2o_nasa9_h_mass(al, np.array([Tf + 0.5]))[0] - _h2o_nasa9_h_mass(al, np.array([Tf - 0.5 + 1.0e-9]))[0]) / 1.0
    hl_lo = h0 - cpl * (Tf - Tg)
    hl = np.where(Tg >= Tf, hl_hi, hl_lo)
    return np.clip(hv - hl, 1.5e6, 3.5e6)


def n2_latent_poly(T):
    Tcl = np.clip(np.asarray(T, float), COND_T_PROP_FLOOR, 126.192 - 0.5)
    p1, p2, p3, p4, p5 = -2.137e-8, 7.18e-6, -9.142e-4, 0.05069, -0.809
    L = (p1 * Tcl ** 4 + p2 * Tcl ** 3 + p3 * Tcl ** 2 + p4 * Tcl + p5) * 1.0e6
    return np.maximum(L, 0.0)


def n2_latent_ex(T, lowT=1, cl=2000.0):
    T = np.asarray(T, float)
    lo = n2_latent_poly(np.array([COND_N2_LATENT_TA]))[0] + (COND_N2_CPV - cl) * (T - COND_N2_LATENT_TA)
    return np.where((lowT == 0) | (T >= COND_N2_LATENT_TA), n2_latent_poly(T), lo)


class CondEOS:
    """凝縮種の二相 EOS 定数 (condProps_H2O / condProps_N2 と同じ R) と潜熱。carrier=True で cond_T_from_e_carrier 形。"""
    def __init__(self, condModel, carrier, latentLowT=1, liquidCp=2000.0):
        self.model = int(condModel); self.carrier = bool(carrier)
        self.Rw = 461.5 if self.model == 1 else 296.8
        self.latentLowT = latentLowT; self.liquidCp = liquidCp

    def latent(self, T):
        if self.model == 1:
            return h2o_latent(T)
        # carrier N2 は n2_latent_ex、pure onetemp は旧多項式 n2_latent (ソルバと同じ)
        return n2_latent_ex(T, self.latentLowT, self.liquidCp) if self.carrier else n2_latent_poly(T)

    def e_liquid_term(self, T, g, Rmix):
        """e_mix − e_gas: carrier は g(R_w T − L), pure onetemp は g(R_mix T − L)。"""
        R = self.Rw if self.carrier else Rmix
        return g * (R * T - self.latent(T))


def T_from_e(gas, Y, e, T0, g=None, eos=None, Tmin=50.0, Tmax=6000.0):
    """e = e_gas(T) [+ 液相項] を Newton で反転 (ベクトル)。g=None/eos=None なら乾き気相。"""
    T = np.clip(np.asarray(T0, float).copy(), Tmin, Tmax); R = gas.Rmix(Y)
    g = np.zeros_like(T) if g is None else np.asarray(g, float)
    for _ in range(80):
        f = gas.h(Y, T) - R * T - e
        dfdT = gas.cp(Y, T) - R
        if eos is not None:
            f = f + eos.e_liquid_term(T, g, R)
            dfdT = dfdT + (eos.e_liquid_term(T + 0.1, g, R) - eos.e_liquid_term(T - 0.1, g, R)) / 0.2
        dT = np.clip(f / np.maximum(dfdT, 1.0e-2 * np.maximum(R, 1.0)), -0.5 * T, 0.5 * T)
        T = np.clip(T - dT, Tmin, Tmax)
        if np.max(np.abs(dT)) < 1e-10 * np.max(T):
            break
    return T


# ----------------------------------------------------------------------------- 配置の読込
def load_layout(meta_path, run_dir, label):
    """{names, expansion, streams, condensing, tracer, MW, db, Tref, condModel, condGasIndex, condensation, has_cfg}。"""
    meta = load_yaml_str(meta_path) if meta_path else None
    has_cfg = bool(run_dir) and os.path.exists(os.path.join(run_dir, "solverConfig.yaml"))
    info = species_info(run_dir) if has_cfg else None
    if meta is None and info is None:
        raise SystemExit(f"{label}: species_meta.yaml も solverConfig.yaml も無い (--meta / --src-meta / --src-run / --dst-run)")
    names = [_up(s) for s in (meta["species"] if meta else info["names"])]
    if info and [_up(s) for s in info["names"]] != names:
        raise SystemExit(f"{label}: species_meta.yaml の species {names} と solverConfig の {info['names']} が違う")
    exp = {}
    for s in names:
        row = (meta or {}).get("expansion", {}).get(s) if meta else None
        if row is None:
            row = {s: 1.0}
        exp[s] = {_up(k): float(v) for k, v in row.items()}
        tot = sum(exp[s].values())
        if abs(tot - 1.0) > 1e-9:
            raise SystemExit(f"{label}: expansion[{s}] の重み和 {tot} が 1 でない")
    streams = {}
    for st, v in ((meta or {}).get("streams") or {}).items():
        yt = [float(x) for x in v.get("Y_transport", [])]
        if len(yt) != len(names):
            raise SystemExit(f"{label}: streams[{st}].Y_transport の長さ {len(yt)} が species {len(names)} と違う")
        streams[str(st)] = np.array(yt)
    cond = None
    if meta and meta.get("condensing_species"):
        cond = _up(meta["condensing_species"])
    elif info and (info["condensing"] or info["h2o_index"] is not None):
        cond = _up(info["condensing"] or info["names"][info["h2o_index"]])
    tracer = bool(((meta or {}).get("tracer") or {}).get("enabled")) if meta else bool(info and info["tracer"])
    MW = {_up(k): float(v) for k, v in ((meta or {}).get("MW") or (info["MW"] if info else {})).items()}
    db = None; Tref = 0.0; condModel = 1; condGasIndex = None; condensation = False
    if has_cfg:
        cfg = load_yaml_str(os.path.join(run_dir, "solverConfig.yaml"))
        pp = cfg.get("physProp") or {}
        Tref = float(pp.get("thermoHrefTemp", 0.0))
        db_file = pp.get("speciesDBFile")
        if db_file:
            p = db_file if os.path.isabs(db_file) else os.path.join(run_dir, db_file)
            if os.path.exists(p):
                db = {_up(k): v for k, v in (load_yaml_str(p) or {}).items()}
        condensation = bool(info["condensation"]); condModel = int(info["condModel"])
        condGasIndex = info["condensing_index"]
    return {"names": names, "expansion": exp, "streams": streams, "condensing": cond, "tracer": tracer, "MW": MW,
            "db": db, "Tref": Tref, "run_dir": run_dir, "condModel": condModel, "condGasIndex": condGasIndex,
            "condensation": condensation, "has_cfg": has_cfg}


def eos_for(layout):
    """layout の凝縮設定から CondEOS (凝縮 OFF なら None)。carrier = condGasSpecies>=0。"""
    if not layout["condensation"]:
        return None
    return CondEOS(layout["condModel"], carrier=(layout["condGasIndex"] is not None))


def gas_for(layout):
    if layout["db"] is None:
        return None
    return _TPGas(dict(layout["db"]), layout["names"], layout["Tref"])


def inflow_label_index(layout):
    """純流入ラベルになっている輸送種 (inflow Y_transport=1, external=0) の index。無ければ None。"""
    st = layout["streams"]
    if "inflow" not in st:
        return None
    yin = st["inflow"]; yext = st.get("external", np.zeros_like(yin))
    cand = [j for j in range(len(yin)) if abs(yin[j] - 1.0) < 1e-12 and abs(yext[j]) < 1e-12]
    return cand[0] if len(cand) == 1 else None


# ----------------------------------------------------------------------------- 実種の質量と移送
def real_masses(layout, Y):
    """輸送種の質量分率 Y[ns, n] → 実種ごとの質量分率 {real: array}。"""
    m = {}
    for s, ss in enumerate(layout["names"]):
        for r, w in layout["expansion"][ss].items():
            m[r] = m.get(r, 0.0) + w * Y[s]
    return m


def transfer_conserve(src, dst):
    """保存的展開/縮約の行列 T[nd, ns]。同名・同 expansion は 1:1、他は実種ごとに一意の行き先を要求する。"""
    homes = {}
    for j, dj in enumerate(dst["names"]):
        for r in dst["expansion"][dj]:
            homes.setdefault(r, []).append(j)
    ns, nd = len(src["names"]), len(dst["names"])
    T = np.zeros((nd, ns)); missing = []; ambiguous = []
    for s, ss in enumerate(src["names"]):
        if ss in dst["names"]:
            j = dst["names"].index(ss)
            ea, eb = src["expansion"][ss], dst["expansion"][ss]
            if set(ea) == set(eb) and all(abs(ea[k] - eb[k]) < 1e-9 for k in ea):
                T[j, s] = 1.0
                continue
        for r, w in src["expansion"][ss].items():
            hj = homes.get(r, [])
            if len(hj) == 0:
                missing.append((ss, r))
            elif len(hj) > 1:
                ambiguous.append((ss, r, [dst["names"][j] for j in hj]))
            else:
                T[hj[0], s] += w
    if missing:
        raise SystemExit("REFUSED: source の実種に destination の行き先が無い: "
                         + ", ".join(f"{r} (from {ss})" for ss, r in missing)
                         + f". destination species {dst['names']} に入れるか lump に含めること")
    if ambiguous:
        raise SystemExit("REFUSED (--mode conserve): 実種の行き先が一意でない (流れ lump): "
                         + ", ".join(f"{r} (from {ss}) -> {hj}" for ss, r, hj in ambiguous)
                         + ". 流れによる再初期化 --mode reinit を使うこと")
    return T


def stream_fraction(src, Ysrc, roXi_src, ro):
    """流入元分率 ξ [n]: source roXi があればそれ、無ければ純流入ラベル種の Y。導けなければ None。"""
    if roXi_src is not None:
        return np.clip(roXi_src / ro, 0.0, 1.0), "source roXi"
    j = inflow_label_index(src)
    if j is not None:
        return np.clip(Ysrc[j], 0.0, 1.0), f"Y_{src['names'][j]} (pure inflow label in source species_meta)"
    return None, None


def transfer_reinit(src, dst, Ysrc, xi):
    """流れによる再初期化: 実種質量 m_r を w_j(r) = (ξ Yt_in[j] + (1−ξ) Yt_ext[j]) exp_j[r] で destination へ分配。返り値 Y_dst[nd, n]。"""
    if "inflow" not in dst["streams"] or "external" not in dst["streams"]:
        raise SystemExit("REFUSED (--mode reinit): destination species_meta.yaml に streams.inflow/external の Y_transport が無い")
    yin, yext = dst["streams"]["inflow"], dst["streams"]["external"]
    n = Ysrc.shape[1]; nd = len(dst["names"])
    m = real_masses(src, Ysrc)
    Yd = np.zeros((nd, n))
    for r, mr in m.items():
        w = np.zeros((nd, n))
        for j, dj in enumerate(dst["names"]):
            e = dst["expansion"][dj].get(r)
            if e is None:
                continue
            w[j] = (xi * yin[j] + (1.0 - xi) * yext[j]) * e
        wsum = w.sum(axis=0)
        bad = (wsum <= 0.0) & (mr > 1e-12)
        if np.any(bad):
            raise SystemExit(f"REFUSED (--mode reinit): 実種 {r} の質量があるのに、その ξ で重みが 0 の destination 種しか無いセルが {int(bad.sum())} 個")
        Yd += mr * np.where(wsum > 0.0, w / np.maximum(wsum, 1e-300), 0.0)
    return Yd


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--meta", required=True, help="destination の species_meta.yaml")
    ap.add_argument("--src-meta", help="source の species_meta.yaml (無ければ SRC run dir の solverConfig.yaml)")
    ap.add_argument("--src-run", help="source run dir (既定: SRC h5 の隣)")
    ap.add_argument("--dst-run", help="destination run dir (既定: DST h5 の隣)")
    ap.add_argument("--mode", choices=["conserve", "reinit"], default="conserve",
                    help="conserve: 保存的展開/縮約 (行き先一意), reinit: 流れ (ξ) による再初期化")
    ap.add_argument("--reconstruct-roe", choices=["auto", "always", "never"], default="auto",
                    help="roe の差分再構成 (auto: DB/datum/種集合が変わるときだけ)")
    ap.add_argument("--T-tol", type=float, default=0.05, help="T 保存検査の許容 [K] (全セル)")
    ap.add_argument("--reinit-tol", type=float, default=1e-3, help="reinit の射影損失の許容 (実種 ρY は ×max ρ [kg/m³], 総水量は相対)")
    ap.add_argument("--drop-moments", action="store_true", help="destination が凝縮 OFF のとき source の液相モーメントを捨てる (既定は拒否)")
    ap.add_argument("--dry-run", action="store_true", help="書き込まず検査だけ")
    a = ap.parse_args()

    src_run = a.src_run or os.path.dirname(os.path.abspath(a.src))
    dst_run = a.dst_run or os.path.dirname(os.path.abspath(a.dst))
    src = load_layout(a.src_meta, src_run, "source")
    dst = load_layout(a.meta, dst_run, "destination")
    print(f"[convert] mode={a.mode}  source {src['names']}  ->  destination {dst['names']}")
    gs, gd = gas_for(src), gas_for(dst)
    eos_s, eos_d = eos_for(src), eos_for(dst)

    # ---- source 読込 ----
    with h5py.File(a.src, "r") as f:
        V = f["VALUE"]; ro = np.array(V["ro"], np.float64); n = ro.shape[0]
        is_res = "P" in V and "Ux" in V
        ns = len(src["names"])
        if is_res:
            Ux, Uy, Uz = (np.array(V[k], np.float64) for k in ("Ux", "Uy", "Uz"))
            Tsrc = np.array(V["T"], np.float64)
            Ysrc = [np.array(V[f"Y{s}"], np.float64) if f"Y{s}" in V else None for s in range(ns)]
            roe = np.array(V["roe"], np.float64) if "roe" in V else None
            roK = ro * np.array(V["k"], np.float64) if "k" in V else None
            roOm = ro * np.array(V["omega"], np.float64) if "omega" in V else None
            moments = {"ro" + k: ro * np.array(V[k], np.float64) for k in V if k.startswith(("g_", "Q0_", "Q1_", "Q2_"))}
            roXi = np.array(V["roXi"], np.float64) if "roXi" in V else (ro * np.array(V["Xi"], np.float64) if "Xi" in V else None)
        else:
            roUx, roUy, roUz = (np.array(V[k], np.float64) for k in ("roUx", "roUy", "roUz"))
            Ux, Uy, Uz = roUx / ro, roUy / ro, roUz / ro
            roe = np.array(V["roe"], np.float64)
            Ysrc = [np.array(V[f"roY{s}"], np.float64) / ro if f"roY{s}" in V else None for s in range(ns)]
            roK = np.array(V["roK"], np.float64) if "roK" in V else None
            roOm = np.array(V["roOmega"], np.float64) if "roOmega" in V else None
            moments = {k: np.array(V[k], np.float64) for k in V if k.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_"))}
            roXi = np.array(V["roXi"], np.float64) if "roXi" in V else None
            Tsrc = None
        vol = np.array(V["volume"], np.float64) if "volume" in V else None
    if ns == 1 and Ysrc[0] is None:
        Ysrc[0] = np.ones(n)
    for s, y in enumerate(Ysrc):
        if y is None:
            raise SystemExit(f"source に Y{s}/roY{s} ({src['names'][s]}) が無い")
    Ysrc = np.array(Ysrc)
    ssum = Ysrc.sum(axis=0)
    print(f"[convert] source ΣY: min {ssum.min():.9f} max {ssum.max():.9f} (正規化して使う)")
    Ysrc = Ysrc / np.maximum(ssum, 1e-30)
    ke = 0.5 * (Ux ** 2 + Uy ** 2 + Uz ** 2)
    g_src = sum(v for k, v in moments.items() if k.startswith("rog_")) / ro if moments else np.zeros(n)

    # source T: source DB + 二相 EOS で roe から反転する (codex M2)。res の T はソルバが陰解法更新の前に評価した値で
    # 保存量 roe より 1 更新ぶん遅れる (未収束の過渡では ~1 K 違う) ので、参照は必ず roe と整合する反転値にし、
    # res の T との差は情報として出す。DB が無いときだけ res の T を使う。
    if gs is not None and roe is not None:
        Tinv = T_from_e(gs, list(Ysrc), roe / ro - ke, (Tsrc if Tsrc is not None else np.full(n, 300.0)), g=g_src, eos=eos_s)
        if Tsrc is not None:
            print(f"[convert] source T (res) と source roe の反転値の差: max {np.max(np.abs(Tinv - Tsrc)):.3e} K (情報; 未収束の過渡では非零)")
        Tsrc = Tinv
        print(f"[convert] source T を SRC DB {'+ 二相 EOS' if eos_s is not None else '(乾き)'} で反転: "
              f"{Tsrc.min():.2f}..{Tsrc.max():.2f} K (湿潤セル {int((g_src > 0).sum())})")
    elif Tsrc is None:
        raise SystemExit("source が input h5 で T が無く、source の species_db.yaml も読めない (--src-run)")

    # ---- 移送 ----
    xi, xi_how = stream_fraction(src, Ysrc, roXi, ro)
    if a.mode == "conserve":
        T = transfer_conserve(src, dst)
        for j, dj in enumerate(dst["names"]):
            parts = [f"{T[j, s]:.6g}*{src['names'][s]}" for s in range(ns) if T[j, s] > 0]
            print(f"[convert]   {dj:8s} <- " + (" + ".join(parts) if parts else "0 (not present in source)"))
        Ydst = T @ Ysrc
    else:
        if xi is None:
            raise SystemExit("REFUSED (--mode reinit): 流入元分率 ξ が導けない (source に roXi も純流入ラベル種も無い)")
        print(f"[convert]   ξ = {xi_how}: min {xi.min():.6g} max {xi.max():.6g} mean {xi.mean():.6g}")
        Ydst = transfer_reinit(src, dst, Ysrc, xi)
    nd = len(dst["names"])

    # ---- トレーサ ----
    roXi_out = None
    if dst["tracer"]:
        if roXi is not None:
            roXi_out = roXi.copy(); how = "copied from source roXi"
        elif xi is not None:
            roXi_out = ro * xi; how = f"generated roXi = ρ·ξ with ξ = {xi_how}"
        else:
            raise SystemExit("REFUSED: destination は tracer (roXi) を要るが source に roXi も純流入ラベル種も無い")
        print(f"[convert]   tracer: {how}")
    elif roXi is not None:
        print("[convert]   note: source roXi is dropped (destination has no tracer)")

    # ---- 凝縮モーメント ----
    if moments and not dst["condensation"]:
        if not a.drop_moments:
            raise SystemExit("REFUSED: source に液相モーメントがあるが destination は凝縮 OFF (--drop-moments で捨てる; 液相は蒸気に戻り T は保たれるが energy は変わる)")
        print("[convert]   WARNING: 液相モーメントを捨てる (--drop-moments)")
        moments_out = {}; g_dst = np.zeros(n)
    else:
        if moments and src["condensing"] and dst["condensing"] and src["condensing"] != dst["condensing"]:
            raise SystemExit(f"REFUSED: 凝縮種が違う (source {src['condensing']} / destination {dst['condensing']}); 凝縮モーメントを移せない")
        moments_out = dict(moments); g_dst = g_src

    # ---- 検査 1: ΣY, 実種保存, 総水量 ----
    fails = []
    dsum = Ydst.sum(axis=0)
    err_sum = float(np.max(np.abs(dsum - 1.0)))
    print(f"[convert] check ΣY_dst−1: max |{err_sum:.3e}| (tol 1e-6)")
    if err_sum > 1e-6:
        fails.append(f"ΣρY=ρ が破れる (max |ΣY−1| {err_sum:.3e})")
    ms, md = real_masses(src, Ysrc), real_masses(dst, Ydst)
    worst = 0.0
    for r in set(ms) | set(md):
        worst = max(worst, float(np.max(np.abs(ms.get(r, 0.0) - md.get(r, 0.0)) * ro)))
    lossy = (a.mode == "reinit")
    tol_m = (a.reinit_tol if lossy else 1e-6) * float(np.max(ro))
    print(f"[convert] check 実種ごとの ρY {'射影損失 (reinit は lump 内組成固定のため非保存)' if lossy else '保存'}: max |Δ| {worst:.3e} kg/m³ (tol {tol_m:.1e})")
    if worst > tol_m:
        fails.append(f"実種の質量{'の射影損失が大きい' if lossy else 'が保存されない'} (max |Δ ρY| {worst:.3e})")
    cond_name = dst["condensing"] or src["condensing"]
    if cond_name and cond_name in ms and cond_name in md:
        wt = vol if vol is not None else np.ones(n)
        tot_s = float(np.sum(ms[cond_name] * ro * wt)); tot_d = float(np.sum(md[cond_name] * ro * wt))
        rel = abs(tot_d - tot_s) / max(abs(tot_s), 1e-300)
        tol_w = a.reinit_tol if lossy else 1e-9
        print(f"[convert] check 総水量 ρ·Y_{cond_name} (液相込みの総水分率; {'体積重み' if vol is not None else 'CV 単純和'}): "
              f"source {tot_s:.9e} destination {tot_d:.9e} rel diff {rel:.3e} (tol {tol_w:.0e})")
        if rel > tol_w:
            fails.append(f"総水量が保存されない (rel {rel:.3e})")
    if roXi_out is not None:
        xr = roXi_out / ro
        print(f"[convert] check roXi/ρ ∈ [0,1]: min {xr.min():.6g} max {xr.max():.6g}")
        if xr.min() < -1e-6 or xr.max() > 1.0 + 1e-6:
            fails.append("roXi/ρ が [0,1] を外れる")
        roXi_out = np.clip(roXi_out, 0.0, ro)

    # ---- roe ----
    differs, why = _db_differs(src, dst)
    do_rec = (a.reconstruct_roe == "always") or (a.reconstruct_roe == "auto" and (differs or a.mode == "reinit" or (moments and not moments_out)))
    if do_rec:
        if gd is None:
            raise SystemExit("roe 再構成に destination の species_db.yaml が要る (--dst-run)")
        Yl = list(Ydst); Rd = gd.Rmix(Yl)
        e_dst = gd.h(Yl, Tsrc) - Rd * Tsrc
        if gs is not None and roe is not None and (moments_out or not moments):
            Ys = list(Ysrc); e_src = gs.h(Ys, Tsrc) - gs.Rmix(Ys) * Tsrc
            roe_new = roe + ro * (e_dst - e_src); how = "差分形 roe += ρ[e_gas,dst(T) − e_gas,src(T)] (液相項不変)"
        else:
            liq = eos_d.e_liquid_term(Tsrc, g_dst, Rd) if eos_d is not None else 0.0
            roe_new = ro * (e_dst + liq + ke); how = "完全再構成 ρ(e_gas,dst(T) + 液相項 + u²/2)"
        de = (roe_new - roe) / ro if roe is not None else np.full(n, np.nan)
        print(f"[convert] roe 再構成 ({why}; {how}): Δe max {np.nanmax(np.abs(de)):.3e} J/kg, mean {np.nanmean(np.abs(de)):.3e} J/kg")
        roe_out = roe_new
    else:
        if roe is None:
            raise SystemExit("source に roe が無く再構成も指定されていない (--reconstruct-roe always)")
        roe_out = roe
        print(f"[convert] roe はそのまま ({why})")

    # ---- 検査 2: T 保存 (destination DB + 二相 EOS で反転, 全セル) ----
    if gd is None:
        fails.append("destination の species_db.yaml が読めず T 保存を検査できない")
    else:
        Tchk = T_from_e(gd, list(Ydst), roe_out / ro - ke, Tsrc, g=g_dst, eos=eos_d)
        dT = np.abs(Tchk - Tsrc); wet = g_dst > 0
        print(f"[convert] check T 保存 (destination DB{' + 二相 EOS' if eos_d is not None else ''} で roe を反転): "
              f"max |ΔT| 全セル {dT.max():.3e} K, 乾き {dT[~wet].max() if (~wet).any() else 0:.3e} K, "
              f"湿潤 {dT[wet].max() if wet.any() else 0:.3e} K ({int(wet.sum())} セル) (tol {a.T_tol} K)")
        if dT.max() > a.T_tol:
            i = int(np.argmax(dT))
            fails.append(f"T が保存されない (max |ΔT| {dT.max():.3e} K at cell {i}: T_src {Tsrc[i]:.3f}, T_chk {Tchk[i]:.3f}, g {g_dst[i]:.3e})")

    if fails:
        print("[convert] FAILED (書き込みなし):"); [print("   - " + m) for m in fails]
        sys.exit(1)
    if a.dry_run:
        print("[convert] --dry-run: 検査 OK, 書き込みなし"); return

    # ---- 書き込み (同一メッシュ index コピー) ----
    with h5py.File(a.dst, "r+") as d:
        nd_ = d["VALUE/ro"].shape[0]
        if nd_ != n:
            raise SystemExit(f"REFUSED: CV 数が違う (source {n}, destination {nd_}); 同一メッシュの input h5 を指定する")
        dt = d["VALUE/ro"].dtype

        def put(name, arr):
            ds = "VALUE/" + name
            if ds in d:
                d[ds][...] = arr.astype(dt)
            else:
                d.create_dataset(ds, data=arr.astype(dt))

        put("ro", ro); put("roUx", ro * Ux); put("roUy", ro * Uy); put("roUz", ro * Uz); put("roe", roe_out)
        if roK is not None: put("roK", roK)
        if roOm is not None: put("roOmega", roOm)
        for k in list(d["VALUE"].keys()):
            if (k.startswith("roY") and k[3:].isdigit()) or k.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or k == "roXi":
                del d["VALUE/" + k]
        for j in range(nd):
            put(f"roY{j}", ro * Ydst[j])
        for k, v in moments_out.items():
            put(k, v)
        if roXi_out is not None:
            put("roXi", roXi_out)
        moved = ["ro", "roUx", "roUy", "roUz", "roe"] + [f"roY{j}" for j in range(nd)] + list(moments_out) \
            + (["roK", "roOmega"] if roK is not None else []) + (["roXi"] if roXi_out is not None else [])
    print(f"[convert] wrote {a.dst}: {moved}")
    print("[convert] SUMMARY: all checks passed (ΣY, real-species mass, total water, T, roXi range)")


def _db_differs(src, dst):
    if src["Tref"] != dst["Tref"]:
        return True, f"thermoHrefTemp {src['Tref']} -> {dst['Tref']}"
    if src["db"] is None or dst["db"] is None:
        return True, "species_db.yaml が片方で読めない"
    if src["names"] != dst["names"]:
        return True, "species set/order changed"
    for n in dst["names"]:
        a, b = src["db"].get(n), dst["db"].get(n)
        if a is None or b is None:
            return True, f"{n} not in one DB"
        for k in ("MW", "nasa9_low", "nasa9_high", "Tmid"):
            if not np.array_equal(np.asarray(a[k], float), np.asarray(b[k], float)):
                return True, f"{n}.{k} differs"
    return False, "same DB and datum"


if __name__ == "__main__":
    main()
