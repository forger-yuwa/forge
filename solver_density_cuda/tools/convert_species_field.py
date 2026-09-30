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
  reinit — **組成の再初期化 (作動点変更)**。source から持ち越すのは各セルの流入元分率 ξ **だけ** (source の `roXi/ρ`、無ければ
      source meta の `exhaust_fraction` が指す種ラベル配列 (例 Y_EXH)、それも無ければ streams.Y_transport の純流入ラベル種) で、
      **source の組成そのものは捨てる**。destination meta の流れ組成 (`streams.inflow.Y` / `streams.external.Y` = 実種質量分率、
      輸送種では `Y_transport`) から各セルの輸送種を Y_t[j] = ξ·Yt_in^dst[j] + (1−ξ)·Yt_ext^dst[j] と作る (実種でも
      Y_real(r) = ξ·Y_in^dst(r) + (1−ξ)·Y_ext^dst(r) と厳密に同じ)。**情報を落とす操作** (plan §4.5 M4, codex result-2 M1):
      新作動点の入口組成に置き換わるので実種の質量・総水量は保存しない (変化量は情報として表示)。ΣY=1・T 保存・有限性は検査する。
  どちらでも destination が `tracer.enabled` で source に `roXi` が無ければ **roXi = ρ ξ を生成**する (ξ が導けなければ拒否)。

エネルギーと温度:
  - source が res (T あり) なら T はソルバの値。input h5 (T なし) なら **ソルバと同じ二相 EOS** で `roe` から反転する:
    e = e_gas(Y_total, T) + g (R_w T − L(T)) (carrier 形, condensationEOS_d.cuh `cond_T_from_e_carrier`)、pure TP は
    e = e_v(T) + g R_mix T − g L(T) (`cond_T_from_e_onetemp`)。N2 の L(T) は condensationProperties_d.cuh の `n2_latent` を移植、
    H2O の L(T) は **run の解決済み記録の気液ペア** (condensed: 液相 H2O(L) とペアの気相種; plan thermophysics-solver-owned-species-db #10)
    から `h2o_latent(CondLatentPair)` と同式で作る。記録に液相が無い湿り場 (#10 以前の H2O 凝縮 run) は潜熱モデルが分からないので
    **既定で拒否**し、旧モデル (#10 以前の `h2o_latent`: 気相 H2O 200–1000 K 係数の多項式外挿・1000 K 頭打ち) で作った場と
    確かめたときだけ `--src-latent legacy-v0` で明示する (移行手順)。
  - DB (熱物性; `forge_species.run_thermo` = ソルバの解決済み記録か従来の speciesDBFile) / datum (`thermoHrefTemp`) / 種集合 /
    **液相 (潜熱モデル)** が変わるときは差分形で
    `roe += ρ {e_gas,dst(Y_dst,T) − e_gas,src(Y_src,T) + g_dst(R_dst T − L_dst(T)) − g_src(R_src T − L_src(T))}`
    (R は carrier 形で R_w、pure 形で R_mix; codex diagnose 2026-09-27: 液相項を入れないと L だけ変わる変換で補正が欠落する)。
検査 (**1 つでも破れば書き込まず失敗終了**; NaN は必ず失敗になるよう有限性を先に見る, codex result-2 M3):
  source の必須データセット (ro, roUx/Ux, roUy, roUz, roe, roY{s}/Y{s} 全種, tracer なら roXi/Xi) の存在、ρ>0 と全保存量・組成・
  T (source/destination) の有限性、組成の非負 (Y < −1e-9 は拒否、|Y| < 1e-9 は 0 にクリップして件数を報告)、|ΣY_src−1| ≤ 1e-4、
  EOS 反転の残差 (|e(T)−e| ≤ 1e-6|e| + 1 J/kg) と括弧端 (T_min=50 K / T_max=6000 K に張り付いたら拒否)、各セル ΣY_dst=1 (1e-6)、
  conserve では実種ごとの ρY 保存 (1e-6·max ρ) と総水量 ρ·Y_w (Y_w は液相込みの総水分率なので ρ(Y_w+g) ではない) の保存 (rel 1e-9)、
  destination DB + 二相 EOS で `roe` を反転した T と source T の差 (乾き・湿潤の全セル, `--T-tol` 既定 0.05 K)、roXi/ρ ∈ [0,1]。
失敗系の試験: tests/unit/test_convert_species_field_fail.py。

- **化学種の属性** (plans/active/thermophysics-solver-owned-species-db.md §4.3, #3b): 入力が属性と検証できる記録を持つときは、
  記録の完全性・SRC run の設定を `forge --resolve-species` で解決したハッシュ = 場の属性、を確かめ、宛先 run を解決して
  変換器が使う宛先の物性 (forge_species.run_thermo の熱物性 + datum) が宛先の記録と一致することを確かめてから、**変換の成功後に宛先のハッシュを付ける**。
  入力が未検証 (属性なし / `species_input_unverified=1`)・宛先を解決できないときは**既定で書き込まずに停止** (ソルバと同じ規約, #3c)。
  許可はその実行だけの `FORGE_ALLOW_UNVERIFIED_SPECIES=1` か `--force-species` で、そのとき変換後も未検証 (属性なし)。
- SRC: res_*.h5 (原始量 P,T,Ux,.. + Y{s}) か input h5 (保存量 roY{s})。DST: 同一メッシュ・同一 CV 数の input h5。
  ro/roU/roe/roK/roOmega・凝縮モーメント `rog_*/roQ*_*` (凝縮種が同名のとき) も index コピーする。
- 両 run dir (`--src-run/--dst-run` 省略時は h5 の隣) の `solverConfig.yaml` と熱物性 (ソルバの記録 `resolved_species_*.yaml`、
  無ければ従来の `speciesDBFile`、それも無ければ `forge --resolve-species`) が必須。lump 記法の種は `species_meta.yaml` が無くても
  記録の lump 構成から展開する: **トレーサの有無と必須保存量
  (`forge_species.required_conserved`) は config から決める**。`species_meta.yaml` は lump の展開・流れ組成・exhaust_fraction に使い、
  config と species の名前/順序または tracer.enabled が矛盾すれば書き込み前に拒否する (codex result-3 M1)。destination config が
  `tracer: exhaust` なら source の `roXi` を必ず持ち越す (conserve) か ρ·ξ で再生成する (reinit)。書き込む配列は 1 つの dict にまとめ、
  全配列 (roK/roOmega/roXi/凝縮モーメント込み) の有限性・ρ>0・非負を最後に検査する (result-3 M3)。
"""
import argparse, os, sys
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from total_quantities import _TPGas, RU  # noqa: E402
from forge_species import species_info, load_yaml_str, species_signature, required_conserved  # noqa: E402
import forge_species as fsp  # noqa: E402


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


def h2o_latent_legacy_v0(T):
    """**旧モデル** (plan thermophysics-solver-owned-species-db #10 以前のソルバの h2o_latent) の移植。移行専用 (--src-latent legacy-v0)。
    気相 H2O の 200–1000 K 係数を再ハードコードし、[45, 1000] K にクランプして多項式のまま評価 (200 K 未満も多項式外挿)、
    液相 H2O(L) は 273.15 K 未満 cp_l 一定外挿・373.15 K 超頭打ち、[1.5e6, 3.5e6] クランプ。"""
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


LEGACY_LATENT_KEY = ("legacy-v0: h2o_latent before plan thermophysics-solver-owned-species-db #10",)


class H2OLatentPair:
    """記録の気液ペア (plan thermophysics-solver-owned-species-db #10) による H2O の潜熱。
    condensationProperties_d.cuh h2o_latent_pair と同式: L = h_v − h_l を [1.5e6, 3.5e6] にクランプ、評価温度の下限 45 K。
    h_v はペアの気相種 (記録の係数; _TPGas と同じ区間・外挿規約)、h_l は記録の液相 (condensed[0]) を気相と同じ MW で質量換算し、
    273.15 K 未満は h(Tlo) から c_p,l = h(Tlo+0.5) − h(Tlo−0.5+1e-9) 一定、373.15 K 超は h(Thi) で頭打ち。
    datum (thermoHrefTemp) は気液に同じ定数を掛けるので差 L には効かない (ここでは両方とも絶対基準で評価する)。"""
    def __init__(self, gas, cond):
        import forge_species as _fsp
        if (cond["extension"]["below"], cond["extension"]["above"]) != ("linear_cp_fd_at_Tlo", "hold_at_Thi") \
                or cond["rule"] != _fsp.SPECIES_CONDENSED_EXTENSION or cond["latent"] != _fsp.SPECIES_CONDENSED_LATENT:
            raise SystemExit(f"REFUSED: 記録の液相 {cond['name']} の延長規約/潜熱規約がこの変換器の実装と違う "
                             f"({cond['extension']}, rule {cond['rule']!r}); 同じ版の tools を使うこと")
        if float(cond["MW"]) != float(gas["MW"]):
            raise SystemExit(f"REFUSED: 記録の液相 {cond['name']} の MW {cond['MW']} がペアの気相 MW {gas['MW']} と違う")
        self.gas = gas; self.MW = float(gas["MW"])
        self.a = np.asarray(cond["nasa9"][:8], float); self.Tlo = float(cond["Tlo"]); self.Thi = float(cond["Thi"])
        self._g = _TPGas({"G": gas}, ["G"], 0.0)
        self.hlLo = self._hl_poly(self.Tlo)
        self.cpl = (self._hl_poly(self.Tlo + 0.5) - self._hl_poly(self.Tlo - 0.5 + 1.0e-9)) / 1.0
        self.hlHi = self._hl_poly(self.Thi)
        _Tb, _co = _fsp.nasa9_intervals(gas)
        if len(_co) == 2:   # 2 区間は #13-1 前と同じキー
            self.key = tuple(_fsp.condensed_compat_lines(cond)) + (repr([gas["MW"], gas["Tlo"], gas["Tmid"], gas["Thi"],
                                                                         list(gas["nasa9_low"]), list(gas["nasa9_high"])]),)
        else:
            self.key = tuple(_fsp.condensed_compat_lines(cond)) + (repr([gas["MW"], list(_Tb), [list(a) for a in _co]]),)

    def _hl_poly(self, T):
        a = self.a; T = np.asarray(T, float)
        hRT = (-a[0] / (T * T) + a[1] * np.log(T) / T + a[2] + a[3] * T / 2.0 + a[4] * T * T / 3.0
               + a[5] * T * T * T / 4.0 + a[6] * T * T * T * T / 5.0 + a[7] / T)
        return hRT * (RU / self.MW) * T   # ソルバ cond_liquid_h_abs_poly と同じ演算順 (液相多項式は打ち消しが強く、順序で ~1e-11 相対ずれる)

    def h_liquid(self, T):
        T = np.asarray(T, float)
        return np.where(T < self.Tlo, self.hlLo - self.cpl * (self.Tlo - T),
                        np.where(T > self.Thi, self.hlHi, self._hl_poly(np.clip(T, self.Tlo, self.Thi))))

    def h_vapor(self, T):
        T = np.asarray(T, float)
        return self._g._props(self.gas, np.atleast_1d(T).ravel())[1].reshape(T.shape)   # _TPGas は 1 次元配列を取る

    def __call__(self, T):
        Tg = np.maximum(np.asarray(T, float), COND_T_PROP_FLOOR)
        L = np.clip(self.h_vapor(Tg) - self.h_liquid(Tg), 1.5e6, 3.5e6)
        return float(L) if L.shape == () else L


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
    """凝縮種の二相 EOS 定数 (condProps_H2O / condProps_N2 と同じ R) と潜熱。carrier=True で cond_T_from_e_carrier 形。
    H2O の潜熱は h2o_latent (H2OLatentPair = 記録の気液ペア, または移行用の h2o_latent_legacy_v0) を渡す。None は「潜熱モデル不明」で、
    液相 g>0 のセルで評価しようとすると拒否する (#10 以前の記録の湿り場; --src-latent legacy-v0 で明示)。"""
    def __init__(self, condModel, carrier, latentLowT=1, liquidCp=2000.0, h2o_latent=None, h2o_latent_key=None, why_unknown=None):
        self.model = int(condModel); self.carrier = bool(carrier)
        self.Rw = 461.5 if self.model == 1 else 296.8
        self.latentLowT = latentLowT; self.liquidCp = liquidCp
        self.h2o_latent = h2o_latent; self.why_unknown = why_unknown
        # 潜熱モデルの同一性 (変換で液相モデルだけ違うかの判定; N2 はソルバ内蔵の固定式)
        self.latent_key = ("n2", self.carrier, latentLowT, liquidCp) if self.model != 1 else h2o_latent_key

    def latent(self, T):
        if self.model == 1:
            if self.h2o_latent is None:
                raise SystemExit(f"REFUSED: H2O の潜熱モデルが分からない湿り場 ({self.why_unknown})。"
                                 "#10 (種 DB の気液ペア) 以前の run の場なら、旧モデルで作った場と確かめたうえで "
                                 "--src-latent legacy-v0 を付けて変換する (移行手順)")
            return self.h2o_latent(T)
        # carrier N2 は n2_latent_ex、pure onetemp は旧多項式 n2_latent (ソルバと同じ)
        return n2_latent_ex(T, self.latentLowT, self.liquidCp) if self.carrier else n2_latent_poly(T)

    def e_liquid_term(self, T, g, Rmix):
        """e_mix − e_gas: carrier は g(R_w T − L), pure onetemp は g(R_mix T − L)。g が全セル 0 なら潜熱を評価しない。"""
        g = np.asarray(g, float)
        if not np.any(g != 0.0):
            return np.zeros(np.broadcast(np.asarray(T, float), g).shape)
        R = self.Rw if self.carrier else Rmix
        return g * (R * T - self.latent(T))


def roe_delta_per_mass(gs, gd, Ys, Yd, T, g_src, g_dst, eos_s, eos_d):
    """差分形の保存エネルギー補正 Δ(roe)/ρ = e_gas,dst(Y_dst,T) − e_gas,src(Y_src,T) + 液相項_dst − 液相項_src
    (液相項 = g(R T − L(T)); codex diagnose 2026-09-27)。気相と液相の両方を SRC/DST それぞれのモデルで評価するので、
    DB/datum/種集合の変更と液相 (潜熱) モデルだけの変更の両方を正しく移す。"""
    Rs, Rd = gs.Rmix(Ys), gd.Rmix(Yd)
    de = (gd.h(Yd, T) - Rd * T) - (gs.h(Ys, T) - Rs * T)
    ls = eos_s.e_liquid_term(T, g_src, Rs) if eos_s is not None else 0.0
    ld = eos_d.e_liquid_term(T, g_dst, Rd) if eos_d is not None else 0.0
    return de + ld - ls


T_MIN, T_MAX = 50.0, 6000.0   # dependentVariables_d.cu DEPVAR_TMIN/TMAX


def T_from_e(gas, Y, e, T0, g=None, eos=None, Tmin=T_MIN, Tmax=T_MAX, diag=None):
    """e = e_gas(T) [+ 液相項] を Newton で反転 (ベクトル)。g=None/eos=None なら乾き気相。
    diag (dict) を渡すと残差 |e(T)−e| と括弧端張り付きの情報を入れる (呼び手が拒否判定に使う)。NaN 入力は NaN を返す。"""
    e = np.asarray(e, float)
    T = np.clip(np.asarray(T0, float).copy(), Tmin, Tmax); R = gas.Rmix(Y)
    g = np.zeros_like(T) if g is None else np.asarray(g, float)

    def resid(T):
        f = gas.h(Y, T) - R * T - e
        return f + (eos.e_liquid_term(T, g, R) if eos is not None else 0.0)

    for _ in range(80):
        f = resid(T)
        dfdT = gas.cp(Y, T) - R
        if eos is not None:
            dfdT = dfdT + (eos.e_liquid_term(T + 0.1, g, R) - eos.e_liquid_term(T - 0.1, g, R)) / 0.2
        dT = np.clip(f / np.maximum(dfdT, 1.0e-2 * np.maximum(R, 1.0)), -0.5 * T, 0.5 * T)
        T = np.clip(T - dT, Tmin, Tmax)
        fin = np.isfinite(dT)
        if not fin.any() or np.max(np.abs(dT[fin])) < 1e-10 * np.max(T[fin]):
            break
    if diag is not None:
        f = resid(T)
        tol = 1e-6 * np.abs(e) + 1.0
        bad_res = ~(np.abs(f) <= tol)                      # NaN も True
        at_end = (T <= Tmin * (1 + 1e-9)) | (T >= Tmax * (1 - 1e-9)) | ~np.isfinite(T)
        diag.update({"resid_max": float(np.max(np.abs(f))) if np.isfinite(f).all() else float("inf"),
                     "n_bad_resid": int(bad_res.sum()), "n_at_bracket": int(at_end.sum()),
                     "i_bad": int(np.argmax(bad_res | at_end)) if (bad_res | at_end).any() else -1})
    return T


def _amax(x):
    """NaN/Inf が 1 つでもあれば inf (比較で必ず失敗させる)。"""
    x = np.asarray(x, float)
    return float(np.max(np.abs(x))) if np.isfinite(x).all() else float("inf")


def _check_finite(fails, name, arr):
    arr = np.asarray(arr, float)
    n = int((~np.isfinite(arr)).sum())
    if n:
        fails.append(f"{name}: 非有限値が {n} 個 (最初 index {int(np.argmax(~np.isfinite(arr)))})")
    return n == 0


def _invert_checked(fails, label, gas, Y, e, T0, g, eos):
    """反転 + 残差/括弧端の検査 (失敗は fails に積む)。"""
    d = {}
    T = T_from_e(gas, Y, e, T0, g=g, eos=eos, diag=d)
    if d["n_bad_resid"] or d["n_at_bracket"]:
        i = d["i_bad"]
        fails.append(f"{label}: EOS 反転が収束しないか括弧端に張り付く (残差超過 {d['n_bad_resid']} セル, 端 {d['n_at_bracket']} セル; "
                     f"例 cell {i}: T {T[i]:.3f} K, e {float(np.asarray(e)[i]):.6g} J/kg, g {float(np.asarray(g)[i]) if g is not None else 0:.3e})")
    return T


# ----------------------------------------------------------------------------- 配置の読込
def load_layout(meta_path, run_dir, label, h5=None, forge=None, resolve_latent=False):
    """{names, expansion, streams, condensing, tracer, MW, db, Tref, condModel, condGasIndex, condensation, has_cfg, h2o_latent, ...}。
    熱物性 (db, Tref) は forge_species.run_thermo (h5 の属性が指す記録 > run dir の記録 > speciesDBFile > --resolve-species)。
    H2O の潜熱 (plan #10) は熱物性の記録の液相 (condensed) から作る。resolve_latent=True (宛先) で記録に液相が無い (speciesDBFile 経由) ときは
    `forge --resolve-species` の記録から液相を取り、ペアの気相が熱物性と同一なことを確かめて使う (宛先はこれから今のソルバで回すため)。"""
    meta = load_yaml_str(meta_path) if meta_path else None
    has_cfg = bool(run_dir) and os.path.exists(os.path.join(run_dir, "solverConfig.yaml"))
    if not has_cfg:
        raise SystemExit(f"{label}: solverConfig.yaml が無い ({run_dir}); 必須保存量・トレーサ・DB は config から決めるので run dir が要る (--src-run / --dst-run)")
    info = species_info(run_dir)
    # 署名 (config + 解決済み DB): tracer の有無と必須保存量はここから決める。species_meta.yaml との矛盾 (名前/順序/tracer) は拒否
    # (codex 2026-09-16 result-3 M1)。
    try:
        sig = species_signature(run_dir)
    except Exception as e:  # noqa: BLE001
        raise SystemExit(f"{label}: 化学種署名が解決できない ({run_dir}): {e}")
    if meta is not None and meta.get("species") is not None and [_up(x) for x in meta["species"]] != [_up(x) for x in sig["names"]]:
        raise SystemExit(f"{label}: species_meta.yaml ({meta_path}) の species {meta['species']} と solverConfig.yaml の physProp.species {sig['names']} が矛盾する (REFUSED)")
    if meta is not None and (meta.get("tracer") or {}).get("enabled") is not None \
            and bool((meta.get("tracer") or {}).get("enabled")) != bool(sig["tracer"]):
        raise SystemExit(f"{label}: species_meta.yaml ({meta_path}) の tracer.enabled={(meta.get('tracer') or {}).get('enabled')} と solverConfig.yaml の physProp.tracer={sig['tracer']} が矛盾する (REFUSED)")
    names = [_up(s) for s in (meta["species"] if meta else info["names"])]
    if info and [_up(s) for s in info["names"]] != names:
        raise SystemExit(f"{label}: species_meta.yaml の species {names} と solverConfig.yaml の physProp.species {info['names']} が矛盾する (REFUSED)")
    # 熱物性 (plan thermophysics-solver-owned-species-db #8): species_db.yaml を前提にしない
    th = None; th_why = None
    try:
        th = fsp.run_thermo(run_dir, res_path=h5)
    except ValueError as e:
        th_why = str(e)
    exp = {}
    for s in names:
        row = (meta or {}).get("expansion", {}).get(s) if meta else None
        if row is None and th is not None:
            row = fsp.lump_mass_expansion(th, s)     # config の lump 記法 (記録の lump 内モル分率と構成種 MW から)
        if row is None:
            row = {s: 1.0}
        exp[s] = {_up(k): float(v) for k, v in row.items()}
        tot = sum(exp[s].values())
        if abs(tot - 1.0) > 1e-9:
            raise SystemExit(f"{label}: expansion[{s}] の重み和 {tot} が 1 でない")
    streams = {}; stream_Y = {}
    for st, v in ((meta or {}).get("streams") or {}).items():
        yt = [float(x) for x in v.get("Y_transport", [])]
        if len(yt) != len(names):
            raise SystemExit(f"{label}: streams[{st}].Y_transport の長さ {len(yt)} が species {len(names)} と違う")
        if not np.isfinite(yt).all() or abs(sum(yt) - 1.0) > 1e-6 or min(yt) < 0.0:
            raise SystemExit(f"{label}: streams[{st}].Y_transport が非有限/負/和≠1 ({sum(yt)!r})")
        streams[str(st)] = np.array(yt)
        stream_Y[str(st)] = {_up(k): float(x) for k, x in (v.get("Y") or {}).items()}
    xi_spec = (meta or {}).get("exhaust_fraction") if meta else None
    cond = None
    if meta and meta.get("condensing_species"):
        cond = _up(meta["condensing_species"])
    elif info and (info["condensing"] or info["h2o_index"] is not None):
        cond = _up(info["condensing"] or info["names"][info["h2o_index"]])
    tracer = bool(sig["tracer"])   # config が正 (meta との矛盾は上で拒否済み)
    MW = {_up(k): float(v) for k, v in ((meta or {}).get("MW") or (info["MW"] if info else {})).items()}
    db = None; Tref = 0.0; condModel = 1; condGasIndex = None; condensation = False
    if has_cfg:
        cfg = load_yaml_str(os.path.join(run_dir, "solverConfig.yaml"))
        pp = cfg.get("physProp") or {}
        Tref = float(pp.get("thermoHrefTemp", 0.0))
        if th is not None:
            db = {_up(k): v for k, v in th["species"].items()}
            Tref = float(th["thermoHrefTemp"])
            print(f"[convert] {label} thermophysics: {th['source']} ({th['how']})")
        elif int(pp.get("thermalMethod", 0)) == 2:
            print(f"[convert] {label} thermophysics: unavailable ({th_why})")
        condensation = bool(info["condensation"]); condModel = int(info["condModel"])
        condGasIndex = info["condensing_index"]
    # H2O の潜熱 (plan #10): 解決済み記録の気液ペア (condensed) から作る。記録に液相が無ければ潜熱モデルは不明 (None)。
    h2o_lat = None; h2o_key = None; lat_why = None
    if condensation and condModel == 1:
        c = th.get("condensed") if th is not None else None
        if c is None and resolve_latent and db is not None:
            try:
                r = fsp.resolve_species(run_dir, forge, inplace=False)
            except (fsp.SpeciesResolveUnavailable, fsp.SpeciesCheckError) as e:
                r = None; lat_why = f"{label}: 記録に液相が無く --resolve-species もできない ({e})"
            rc = (r or {}).get("record") or {}
            if rc.get("condensed") is not None:
                gi = int(rc["condensed"]["gas_index"])
                ge = rc["species"][gi] if 0 <= gi < len(rc["species"]) else None
                ga = db.get(names[gi]) if 0 <= gi < len(names) else None
                same = ge is not None and ga is not None and float(ge["MW"]) == float(ga["MW"]) \
                    and fsp.nasa9_intervals(ge) == fsp.nasa9_intervals(ga)   # 区間 (数・境界) と全係数 (#13-1)
                if not same:
                    raise SystemExit(f"REFUSED: {label}: 熱物性 ({th['how'] if th else '?'}) のペアの気相が --resolve-species の記録と違う; 潜熱の気液ペアを組めない")
                c = rc["condensed"]
                print(f"[convert] {label} liquid phase taken from forge --resolve-species (the thermophysics source {th['source']} has none)")
        if c is not None:
            gi = int(c["gas_index"])
            if not (0 <= gi < len(names)) or db is None or names[gi] not in db:
                raise SystemExit(f"{label}: 記録の液相 {c['name']} のペアの気相 (gas_index {gi}) が種リスト {names} に無い")
            h2o_lat = H2OLatentPair(db[names[gi]], c); h2o_key = h2o_lat.key
            print(f"[convert] {label} H2O latent heat: pair {c['name']} + gas {names[gi]} from the record "
                  f"(L(150/250/300 K) = {h2o_lat(150.0):.1f} / {h2o_lat(250.0):.1f} / {h2o_lat(300.0):.1f} J/kg)")
        else:
            lat_why = lat_why or (f"{label}: 熱物性 {th['source'] + ' (' + th['how'] + ')' if th is not None else '(解決できない)'} に液相 (condensed) が無い"
                       " = plan #10 以前の記録か speciesDBFile")
    return {"names": names, "expansion": exp, "streams": streams, "stream_Y": stream_Y, "xi_spec": xi_spec,
            "condensing": cond, "tracer": tracer, "MW": MW, "sig": sig, "required": required_conserved(sig),
            "db": db, "Tref": Tref, "run_dir": run_dir, "condModel": condModel, "condGasIndex": condGasIndex,
            "condensation": condensation, "has_cfg": has_cfg, "label": label,
            "h2o_latent": h2o_lat, "h2o_latent_key": h2o_key, "latent_why": lat_why}


def eos_for(layout):
    """layout の凝縮設定から CondEOS (凝縮 OFF なら None)。carrier = condGasSpecies>=0。"""
    if not layout["condensation"]:
        return None
    return CondEOS(layout["condModel"], carrier=(layout["condGasIndex"] is not None), h2o_latent=layout["h2o_latent"],
                   h2o_latent_key=layout["h2o_latent_key"], why_unknown=layout["latent_why"])


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
    """流入元分率 ξ [n]: source roXi/ρ → source meta の exhaust_fraction (kind species の Y{i}) → 純流入ラベル種の Y。導けなければ None。"""
    if roXi_src is not None:
        return np.clip(roXi_src / ro, 0.0, 1.0), "source roXi/ρ"
    spec = src.get("xi_spec")
    if spec and spec.get("kind") == "species":
        nm = _up(spec.get("species", ""))
        if nm in src["names"]:
            j = src["names"].index(nm)
            return np.clip(Ysrc[j], 0.0, 1.0), f"Y_{nm} (source species_meta exhaust_fraction)"
    j = inflow_label_index(src)
    if j is not None:
        return np.clip(Ysrc[j], 0.0, 1.0), f"Y_{src['names'][j]} (pure inflow label in source species_meta)"
    return None, None


def transfer_reinit(dst, xi):
    """組成の再初期化: Y_t[j] = ξ·Yt_in^dst[j] + (1−ξ)·Yt_ext^dst[j] (destination meta の流れ組成; source の組成は捨てる)。
    Y_transport は実種組成 streams.<st>.Y を destination の輸送種 (lump 込み) に縮約した既定ベクトルなので、実種で
    Y_real(r) = ξ·Y_in(r) + (1−ξ)·Y_ext(r) を作ってから縮約するのと厳密に同じ。返り値 Y_dst[nd, n]。"""
    if "inflow" not in dst["streams"]:
        raise SystemExit("REFUSED (--mode reinit): destination species_meta.yaml に streams.inflow の Y_transport が無い")
    yin = dst["streams"]["inflow"]
    yext = dst["streams"].get("external")
    if yext is None:
        yext = yin
        print("[convert]   note: destination has no external stream; ξ<1 cells also get the inflow composition")
    return np.outer(yin, xi) + np.outer(yext, 1.0 - xi)


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
    ap.add_argument("--src-sum-tol", type=float, default=1e-4, help="source の |ΣY−1| の許容 (これを超える source は壊れているとみなす)")
    ap.add_argument("--drop-moments", action="store_true", help="destination が凝縮 OFF のとき source の液相モーメントを捨てる (既定は拒否)")
    ap.add_argument("--dry-run", action="store_true", help="書き込まず検査だけ")
    ap.add_argument("--forge", help="--resolve-species を持つ forge (既定: FORGE_BIN, solver_density_cuda/build/forge)")
    ap.add_argument("--force-species", action="store_true", help="入力の化学種記録の検証失敗を無視する (変換後は属性なし = 未検証)")
    ap.add_argument("--src-latent", choices=["record", "legacy-v0"], default="record",
                    help="source の H2O 潜熱モデル。record (既定): source の解決済み記録の気液ペア (plan #10)。"
                         "legacy-v0: #10 以前のソルバの h2o_latent (記録に液相が無い湿り場の移行専用; 記録に液相があれば拒否)")
    a = ap.parse_args()

    src_run = a.src_run or os.path.dirname(os.path.abspath(a.src))
    dst_run = a.dst_run or os.path.dirname(os.path.abspath(a.dst))
    src = load_layout(a.src_meta, src_run, "source", h5=a.src)
    dst = load_layout(a.meta, dst_run, "destination", h5=a.dst, forge=a.forge, resolve_latent=True)
    print(f"[convert] mode={a.mode}  source {src['names']}  ->  destination {dst['names']}")
    gs, gd = gas_for(src), gas_for(dst)
    if a.src_latent == "legacy-v0":
        # 移行手順 (plan #10): 記録に液相が無い H2O 湿り場を、#10 以前の潜熱モデルで作った場として読む (明示指定のときだけ)
        if not (src["condensation"] and src["condModel"] == 1):
            raise SystemExit("REFUSED: --src-latent legacy-v0 は source が H2O 凝縮 (condModel 1) の run のときだけ使う")
        if src["h2o_latent"] is not None:
            raise SystemExit("REFUSED: --src-latent legacy-v0 だが source の記録は液相 (気液ペア) を持つ; 旧モデルの指定は記録と矛盾する")
        src["h2o_latent"] = h2o_latent_legacy_v0; src["h2o_latent_key"] = LEGACY_LATENT_KEY; src["latent_why"] = None
        print("[convert] source H2O latent heat: legacy-v0 (h2o_latent before plan #10; explicit --src-latent migration)")
    if dst["condensation"] and dst["condModel"] == 1 and dst["h2o_latent"] is None:
        raise SystemExit(f"REFUSED: destination の H2O 潜熱モデルが無い ({dst['latent_why']}); 宛先は #10 以降の forge "
                         "(--resolve-species) で解決した記録が要る")
    eos_s, eos_d = eos_for(src), eos_for(dst)
    # 化学種の属性 (§4.3): 入力の検証と宛先の解決 (書き込み前)。付ける属性は変換の成功後に書く
    try:
        species_plan = fsp.plan_convert(a.src, dst_run, src_run_dir=src_run, forge=a.forge, force=a.force_species,
                                        tool="convert", inplace=not a.dry_run, legacy_latent=(a.src_latent == "legacy-v0"))
    except fsp.SpeciesCheckError as e:
        raise SystemExit(f"[convert] REFUSED (nothing written): {e}")
    if species_plan["attrs"] is not None:
        # 変換器が宛先の roe を作る/検査する物性 (run_thermo の熱物性 + thermoHrefTemp) がソルバの解決結果と同じか
        drec = species_plan["dst"]["record"]
        if gd is None:
            probs = ["destination thermophysics could not be resolved by the converter"]
        else:
            nd_ = len(dst["names"])
            mixes = [(f"species {nm}", [1.0 if j == i else 0.0 for j in range(nd_)], gd.R[i],
                      (lambda T, _i=i: gd.h([1.0 if j == _i else 0.0 for j in range(nd_)], T) - gd.R[_i] * T))
                     for i, nm in enumerate(dst["names"])]
            probs = fsp.check_ic_against_record(drec, dst["Tref"], dst["names"], [float(sp_["MW"]) for sp_ in gd.sp], mixes)
        if probs:
            raise SystemExit("[convert] REFUSED (nothing written): the converter's destination thermophysics differ from "
                             f"the destination record (forge --resolve-species):\n" + "".join(f"    {x}\n" for x in probs))

    # ---- source 読込 (必須データセットの存在を先に検査; codex result-2 M2) ----
    ns = len(src["names"])
    with h5py.File(a.src, "r") as f:
        V = f["VALUE"]
        is_res = "P" in V and "Ux" in V
        # 必須データセットは source の config 署名 (required_conserved) から決める。res では原始量 (Ux, Y{s}, Xi) で代替可。
        alt = {"roUx": "Ux", "roUy": "Uy", "roUz": "Uz", "roXi": "Xi"}
        need = []
        for k in src["required"]:
            if k in V:
                need.append(k)
            elif is_res and (alt.get(k) or (k[2:] if k.startswith("roY") else None)) in V:
                need.append(alt.get(k) or k[2:])
            else:
                need.append(k)
        missing = [k for k in need if k not in V]
        if missing:
            raise SystemExit(f"REFUSED: source {a.src} に必須データセットが無い: {missing} (species {src['names']}, tracer {src['tracer']})")
        ro = np.array(V["ro"], np.float64); n = ro.shape[0]
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

    # ---- 検査 0: 有限性・ρ>0・組成の非負 (書き込み前 hard fail; codex result-2 M3) ----
    fails = []
    _check_finite(fails, "ro", ro)
    if not (np.isfinite(ro).all() and (ro > 0.0).all()):
        fails.append(f"ρ>0 が破れる (min ρ {np.nanmin(ro):.6g}, ρ<=0 が {int((~(ro > 0.0)).sum())} セル)")
    for nm, arr in (("Ux", Ux), ("Uy", Uy), ("Uz", Uz)):
        _check_finite(fails, nm, arr)
    if roe is not None:
        _check_finite(fails, "roe", roe)
    for s in range(ns):
        _check_finite(fails, f"Y{s} ({src['names'][s]})", Ysrc[s])
    for k, v in moments.items():
        _check_finite(fails, k, v)
        if np.isfinite(v).all() and (v < -1e-9 * np.max(ro)).any():
            fails.append(f"{k}: 負の液相モーメントがある (min {v.min():.3e})")
    if roXi is not None:
        _check_finite(fails, "roXi", roXi)
    for nm, arr in (("roK", roK), ("roOmega", roOm)):
        if arr is not None:
            _check_finite(fails, nm, arr)
            if np.isfinite(arr).all() and (arr < -1e-9 * max(float(np.max(np.abs(arr))), 1e-300)).any():
                fails.append(f"{nm}: 負値がある (min {arr.min():.3e})")
    if Tsrc is not None:
        _check_finite(fails, "T (source res)", Tsrc)
    if fails:
        print("[convert] FAILED (書き込みなし; 入力の有限性/正値):"); [print("   - " + m) for m in fails]
        sys.exit(1)
    neg = Ysrc < -1e-9
    if neg.any():
        i = np.argwhere(neg)[0]
        raise SystemExit(f"REFUSED: source の組成に負値 (Y{i[0]}[{i[1]}] = {Ysrc[i[0], i[1]]:.3e} < -1e-9)")
    tiny = (np.abs(Ysrc) < 1e-9) & (Ysrc != 0.0)
    if tiny.any():
        print(f"[convert] |Y| < 1e-9 の {int(tiny.sum())} 値を 0 にクリップ")
        Ysrc = np.where(tiny, 0.0, Ysrc)
    Ysrc = np.maximum(Ysrc, 0.0)
    ssum = Ysrc.sum(axis=0)
    print(f"[convert] source ΣY: min {ssum.min():.9f} max {ssum.max():.9f} (正規化して使う; tol {a.src_sum_tol:.0e})")
    if _amax(ssum - 1.0) > a.src_sum_tol:
        raise SystemExit(f"REFUSED: source の |ΣY−1| が {_amax(ssum - 1.0):.3e} > {a.src_sum_tol:.0e} (組成が壊れている)")
    Ysrc = Ysrc / ssum
    ke = 0.5 * (Ux ** 2 + Uy ** 2 + Uz ** 2)
    g_src = sum(v for k, v in moments.items() if k.startswith("rog_")) / ro if moments else np.zeros(n)
    g_src = np.maximum(g_src, 0.0)

    # source T: source DB + 二相 EOS で roe から反転する (codex M2)。res の T はソルバが陰解法更新の前に評価した値で
    # 保存量 roe より 1 更新ぶん遅れる (未収束の過渡では ~1 K 違う) ので、参照は必ず roe と整合する反転値にし、
    # res の T との差は情報として出す。DB が無いときだけ res の T を使う。
    if gs is not None and roe is not None:
        Tinv = _invert_checked(fails, "source roe の反転", gs, list(Ysrc), roe / ro - ke,
                               (Tsrc if Tsrc is not None else np.full(n, 300.0)), g_src, eos_s)
        if fails:
            print("[convert] FAILED (書き込みなし):"); [print("   - " + m) for m in fails]
            sys.exit(1)
        if Tsrc is not None:
            print(f"[convert] source T (res) と source roe の反転値の差: max {np.max(np.abs(Tinv - Tsrc)):.3e} K (情報; 未収束の過渡では非零)")
        Tsrc = Tinv
        print(f"[convert] source T を SRC DB {'+ 二相 EOS' if eos_s is not None else '(乾き)'} で反転: "
              f"{Tsrc.min():.2f}..{Tsrc.max():.2f} K (湿潤セル {int((g_src > 0).sum())})")
    elif Tsrc is None:
        raise SystemExit("source が input h5 で T が無く、source の熱物性 (記録 / speciesDBFile) も解決できない (--src-run)")

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
            raise SystemExit("REFUSED (--mode reinit): 流入元分率 ξ が導けない (source に roXi も exhaust_fraction も純流入ラベル種も無い)")
        _check_finite(fails, "ξ", xi)
        if fails:
            print("[convert] FAILED (書き込みなし):"); [print("   - " + m) for m in fails]; sys.exit(1)
        print(f"[convert]   ξ = {xi_how}: min {xi.min():.6g} max {xi.max():.6g} mean {xi.mean():.6g}")
        print(f"[convert]   composition re-initialized from destination streams (source composition dropped): "
              f"inflow Y_t {np.round(dst['streams']['inflow'], 6).tolist()}, external Y_t "
              f"{np.round(dst['streams']['external'], 6).tolist() if 'external' in dst['streams'] else 'n/a'}")
        Ydst = transfer_reinit(dst, xi)
    nd = len(dst["names"])

    # ---- トレーサ (destination config が tracer: exhaust なら roXi を必ず書く: 持ち越し (conserve) か再生成 (reinit)) ----
    roXi_out = None
    if dst["tracer"]:
        if a.mode == "reinit":
            if xi is None:
                raise SystemExit("REFUSED: destination は tracer (roXi) を要るが ξ が導けない")
            roXi_out = ro * xi; how = f"regenerated roXi = ρ·ξ with ξ = {xi_how}"
        elif roXi is not None:
            roXi_out = roXi.copy(); how = "carried from source roXi"
        elif xi is not None:
            roXi_out = ro * xi; how = f"generated roXi = ρ·ξ with ξ = {xi_how}"
        else:
            raise SystemExit("REFUSED: destination は tracer (roXi) を要るが source に roXi も exhaust_fraction も純流入ラベル種も無い")
        print(f"[convert]   tracer: {how}")
    elif roXi is not None:
        print("[convert]   note: source roXi is dropped (destination config has no physProp.tracer)")

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

    # ---- 検査 1: 有限性, ΣY, (conserve) 実種保存・総水量, roXi ----
    _check_finite(fails, "Y_dst", Ydst)
    dsum = Ydst.sum(axis=0)
    err_sum = _amax(dsum - 1.0)
    print(f"[convert] check ΣY_dst−1: max |{err_sum:.3e}| (tol 1e-6)")
    if not (err_sum <= 1e-6):
        fails.append(f"ΣρY=ρ が破れる (max |ΣY−1| {err_sum:.3e})")
    if (Ydst < 0.0).any():
        fails.append(f"destination の組成に負値 (min {Ydst.min():.3e})")
    ms, md = real_masses(src, Ysrc), real_masses(dst, Ydst)
    worst = 0.0
    for r in set(ms) | set(md):
        worst = max(worst, _amax((ms.get(r, 0.0) - md.get(r, 0.0)) * ro))
    lossy = (a.mode == "reinit")
    if lossy:
        print(f"[convert] info 実種ごとの ρY の変化 (reinit は組成を destination の流れ組成に置き換えるので保存しない): max |Δ| {worst:.3e} kg/m³")
    else:
        tol_m = 1e-6 * float(np.max(ro))
        print(f"[convert] check 実種ごとの ρY 保存: max |Δ| {worst:.3e} kg/m³ (tol {tol_m:.1e})")
        if not (worst <= tol_m):
            fails.append(f"実種の質量が保存されない (max |Δ ρY| {worst:.3e})")
    cond_name = dst["condensing"] or src["condensing"]
    if cond_name and cond_name in ms and cond_name in md:
        wt = vol if vol is not None else np.ones(n)
        tot_s = float(np.sum(ms[cond_name] * ro * wt)); tot_d = float(np.sum(md[cond_name] * ro * wt))
        rel = abs(tot_d - tot_s) / max(abs(tot_s), 1e-300) if np.isfinite(tot_s) and np.isfinite(tot_d) else float("inf")
        if lossy:
            print(f"[convert] info 総水量 ρ·Y_{cond_name}: source {tot_s:.9e} destination {tot_d:.9e} rel diff {rel:.3e} (reinit では保存しない)")
        else:
            print(f"[convert] check 総水量 ρ·Y_{cond_name} (液相込みの総水分率; {'体積重み' if vol is not None else 'CV 単純和'}): "
                  f"source {tot_s:.9e} destination {tot_d:.9e} rel diff {rel:.3e} (tol 1e-09)")
            if not (rel <= 1e-9):
                fails.append(f"総水量が保存されない (rel {rel:.3e})")
    if roXi_out is not None:
        _check_finite(fails, "roXi", roXi_out)
        xr = roXi_out / ro
        print(f"[convert] check roXi/ρ ∈ [0,1]: min {np.nanmin(xr):.6g} max {np.nanmax(xr):.6g}")
        if not ((xr >= -1e-6).all() and (xr <= 1.0 + 1e-6).all()):
            fails.append("roXi/ρ が [0,1] を外れる")
        roXi_out = np.clip(roXi_out, 0.0, ro)

    # ---- roe ----
    differs, why = _db_differs(src, dst, eos_s, eos_d)
    do_rec = (a.reconstruct_roe == "always") or (a.reconstruct_roe == "auto" and (differs or a.mode == "reinit" or (moments and not moments_out)))
    if do_rec:
        if gd is None:
            raise SystemExit("roe 再構成に destination の熱物性 (記録 / speciesDBFile) が要る (--dst-run)")
        Yl = list(Ydst); Rd = gd.Rmix(Yl)
        e_dst = gd.h(Yl, Tsrc) - Rd * Tsrc
        if gs is not None and roe is not None and (moments_out or not moments):
            roe_new = roe + ro * roe_delta_per_mass(gs, gd, list(Ysrc), Yl, Tsrc, g_src, g_dst, eos_s, eos_d)
            how = "差分形 roe += ρ{e_gas,dst(T) − e_gas,src(T) + g_dst(R T − L_dst) − g_src(R T − L_src)}"
        else:
            liq = eos_d.e_liquid_term(Tsrc, g_dst, Rd) if eos_d is not None else 0.0
            roe_new = ro * (e_dst + liq + ke); how = "完全再構成 ρ(e_gas,dst(T) + 液相項 + u²/2)"
        de = (roe_new - roe) / ro if roe is not None else np.zeros(n)
        print(f"[convert] roe 再構成 ({why}; {how}): Δe max {_amax(de):.3e} J/kg, mean {np.mean(np.abs(de)) if np.isfinite(de).all() else float('nan'):.3e} J/kg")
        roe_out = roe_new
    else:
        if roe is None:
            raise SystemExit("source に roe が無く再構成も指定されていない (--reconstruct-roe always)")
        roe_out = roe
        print(f"[convert] roe はそのまま ({why})")

    # ---- 検査 2: 有限性と T 保存 (destination DB + 二相 EOS で反転, 全セル; NaN は必ず失敗) ----
    _check_finite(fails, "roe_out", roe_out)
    _check_finite(fails, "T (source)", Tsrc)
    if gd is None:
        fails.append("destination の熱物性が解決できず T 保存を検査できない")
    elif not fails:
        Tchk = _invert_checked(fails, "destination roe の反転", gd, list(Ydst), roe_out / ro - ke, Tsrc, g_dst, eos_d)
        _check_finite(fails, "T (destination)", Tchk)
        dT = np.abs(Tchk - Tsrc); wet = g_dst > 0
        print(f"[convert] check T 保存 (destination DB{' + 二相 EOS' if eos_d is not None else ''} で roe を反転): "
              f"max |ΔT| 全セル {_amax(dT):.3e} K, 乾き {_amax(dT[~wet]) if (~wet).any() else 0:.3e} K, "
              f"湿潤 {_amax(dT[wet]) if wet.any() else 0:.3e} K ({int(wet.sum())} セル) (tol {a.T_tol} K)")
        if not (_amax(dT) <= a.T_tol):
            i = int(np.nanargmax(dT)) if np.isfinite(dT).any() else 0
            fails.append(f"T が保存されない (max |ΔT| {_amax(dT):.3e} K at cell {i}: T_src {Tsrc[i]:.3f}, T_chk {Tchk[i]:.3f}, g {g_dst[i]:.3e})")

    # ---- 書き込む配列を 1 つの dict にまとめ、宛先 dtype に変換してから全配列の有限性・正値性を検査 (codex result-3 M3 / result-4 M1) ----
    # 宛先 dtype (既存データセットはその dtype、新規は VALUE/ro の dtype = run の flow_float, 通常 float32) を先に読む。
    # float64 で有限でも float32 へ落とすと Inf になり得る (例 roOmega=1e39) ので、変換後の配列で検査する。
    # ここは読み取りだけ (既存データセットの削除・再作成は全検査を通った後の書き込み段でしか行わない)。--dry-run も同じ検査を通る。
    out64 = {"ro": ro, "roUx": ro * Ux, "roUy": ro * Uy, "roUz": ro * Uz, "roe": roe_out}
    for j in range(nd):
        out64[f"roY{j}"] = ro * Ydst[j]
    out64.update(moments_out)
    if roK is not None: out64["roK"] = roK
    if roOm is not None: out64["roOmega"] = roOm
    if roXi_out is not None: out64["roXi"] = roXi_out
    with h5py.File(a.dst, "r") as d:
        nd_ = d["VALUE/ro"].shape[0]
        if nd_ != n:
            raise SystemExit(f"REFUSED: CV 数が違う (source {n}, destination {nd_}); 同一メッシュの input h5 を指定する")
        dt_default = d["VALUE/ro"].dtype
        dtypes = {k: (d["VALUE/" + k].dtype if ("VALUE/" + k) in d else dt_default) for k in out64}
    out = {}
    for k, v in out64.items():
        v = np.asarray(v, float)
        if v.shape != (n,):
            fails.append(f"{k}: 長さ {v.shape} が CV 数 {n} と違う")
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            vc = v.astype(dtypes[k])                      # 最終 dtype (通常 float32) へ変換してから検査
        out[k] = vc
        vchk = vc.astype(np.float64)
        if not _check_finite(fails, f"write:{k} ({np.dtype(dtypes[k]).name})", vchk):
            if np.isfinite(v).all():
                fails.append(f"write:{k}: float64 では有限だが {np.dtype(dtypes[k]).name} への変換で overflow (max |v| {np.max(np.abs(v)):.3e})")
            continue
        if k == "ro" and not (vchk > 0.0).all():
            fails.append("write:ro に ρ<=0 がある")
        if (k in ("roK", "roOmega") or k.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_"))) \
                and (vchk < -1e-9 * max(float(np.max(np.abs(vchk))), 1e-300)).any():
            fails.append(f"write:{k} に負値がある (min {vchk.min():.3e})")
    missing_req = [k for k in dst["required"] if k not in out64]
    if missing_req:
        fails.append(f"destination config が要求する保存量が揃っていない: {missing_req}")
    if fails:
        print("[convert] FAILED (書き込みなし):"); [print("   - " + m) for m in fails]
        sys.exit(1)
    print(f"[convert] arrays to write ({len(out)}, dtype {np.dtype(dt_default).name}): {list(out)} (checked after cast: all finite; ρ>0; roK/roOmega/moments >= 0)")
    if a.dry_run:
        print("[convert] --dry-run: all checks passed (書き込みなし)" + ("; reinit: composition re-initialized" if lossy else "")); return

    # ---- 書き込み (同一メッシュ index コピー; 全検査通過後にだけ既存データセットを削除/再作成) ----
    with h5py.File(a.dst, "r+") as d:
        fsp.write_species_attrs(d, None)          # 書き込み途中で失敗しても古い属性が残らないように先に消す
        for k in list(d["VALUE"].keys()):
            if (k.startswith("roY") and k[3:].isdigit()) or k.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or k == "roXi":
                del d["VALUE/" + k]
        for k, v in out.items():
            ds = "VALUE/" + k
            if ds in d:
                d[ds][...] = v
            else:
                d.create_dataset(ds, data=v)
        moved = list(out)
        fsp.write_species_attrs(d, species_plan["attrs"])
    print(f"[convert] wrote {a.dst}: {moved}")
    print("[convert] species attributes: " + (f"destination species_hash {species_plan['attrs']['species_hash'][:16]} "
          "(species_input_unverified=0)" if species_plan["attrs"] else "none (input unverified -> output unverified)"))
    print("[convert] SUMMARY: all checks passed (finite, ρ>0, ΣY, " + ("T, roXi range; reinit: composition re-initialized)" if lossy else "real-species mass, total water, T, roXi range)"))


def _db_differs(src, dst, eos_s=None, eos_d=None):
    # 液相 (潜熱モデル) だけの違いも再構成を発火させる (plan #10; codex diagnose 2026-09-27)
    ks = eos_s.latent_key if eos_s is not None else None
    kd = eos_d.latent_key if eos_d is not None else None
    if (eos_s is not None or eos_d is not None) and ks != kd:
        return True, "liquid phase / latent heat model differs (" + ("unknown" if ks is None else str(ks[0])[:60]) + " -> " \
            + ("unknown" if kd is None else str(kd[0])[:60]) + ")"
    if src["Tref"] != dst["Tref"]:
        return True, f"thermoHrefTemp {src['Tref']} -> {dst['Tref']}"
    if src["db"] is None or dst["db"] is None:
        return True, "熱物性が片方で解決できない"
    if src["names"] != dst["names"]:
        return True, "species set/order changed"
    for n in dst["names"]:
        a, b = src["db"].get(n), dst["db"].get(n)
        if a is None or b is None:
            return True, f"{n} not in one DB"
        if not np.array_equal(np.asarray(a["MW"], float), np.asarray(b["MW"], float)):
            return True, f"{n}.MW differs"
        (Ta, ca), (Tb_, cb) = fsp.nasa9_intervals(a), fsp.nasa9_intervals(b)
        if len(ca) == 2 and len(cb) == 2:
            for k in ("nasa9_low", "nasa9_high", "Tmid"):
                if not np.array_equal(np.asarray(a[k], float), np.asarray(b[k], float)):
                    return True, f"{n}.{k} differs"
        elif (Ta, ca) != (Tb_, cb):
            return True, f"{n}: temperature intervals or coefficients differ ({len(ca)} vs {len(cb)} intervals)"
    return False, "same DB and datum"


if __name__ == "__main__":
    main()
