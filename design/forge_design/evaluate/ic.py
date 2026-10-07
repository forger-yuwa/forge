"""準 1 次元等エントロピー初期値の貼り付け (uniform IC は圧力比で発散するため必須)。

case/29.bell_vs_conical/mesh/set_isentropic_ic.py と同方式だが、contour CSV
でなく NozzleWall から直接 r_w(x) を取る。x<0 亜音速根 / x>=0 超音速根。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np


def _forge_species():
    """solver_density_cuda/tools/forge_species.py (化学種の解決済み記録・属性の共通 API)。"""
    tools = str(Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3])) / "solver_density_cuda" / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import forge_species
    return forge_species


def area_ratio(M, g):
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
        (g + 1.0) / (2.0 * (g - 1.0))
    )


def invert_area_ratio(AR, supersonic, g):
    AR = np.maximum(AR, 1.0 + 1e-12)
    lo = np.where(supersonic, 1.0, 1e-4) * np.ones_like(AR)
    hi = np.where(supersonic, 50.0, 1.0) * np.ones_like(AR)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f = area_ratio(mid, g) - AR
        go_hi = np.where(supersonic, f < 0.0, f > 0.0)
        lo = np.where(go_hi, mid, lo)
        hi = np.where(go_hi, hi, mid)
    return 0.5 * (lo + hi)


def _throat_of(wall) -> tuple:
    """IC の 1D 等エントロピーが使うスロート (x, r) [r_t]。物理壁 (`PhysicalNozzleWall`・`SingleBSplinePhysicalWall`) は
    属性 x_throat・r_throat が必須で、欠けていたり非有限なら例外 (黙って設計スロート (0, 1) に戻さない;
    plan tooling-nozzle-wall-single-bspline §4.2)。設計壁など物理スロートを持たない壁は従来どおり (0, 1)。"""
    from ..geometry.wall_axismach import PhysicalNozzleWall, SingleBSplinePhysicalWall
    if isinstance(wall, (PhysicalNozzleWall, SingleBSplinePhysicalWall)):
        try:
            x_thr, r_thr = float(wall.x_throat), float(wall.r_throat)
        except AttributeError as e:
            raise ValueError(f"paste_isentropic_ic: 物理壁 {type(wall).__name__} にスロート属性が無い ({e})") from e
        if not (np.isfinite(x_thr) and np.isfinite(r_thr) and r_thr > 0.0):
            raise ValueError(f"paste_isentropic_ic: 物理壁のスロートが不正 (x_throat={x_thr!r}, r_throat={r_thr!r})")
        return x_thr, r_thr
    return float(getattr(wall, "x_throat", 0.0)), float(getattr(wall, "r_throat", 1.0))


def paste_isentropic_ic(h5path, wall, scale, Pt, Tt, gamma, cp,
                        k_init=1.0, omega_init=18000.0, gas=None,
                        h_ref_T: float | None = None,
                        species_Y: list | None = None) -> dict:
    """wall: NozzleWall (無次元), scale: r* [m]。出口諸元 dict を返す。

    SST の roK/roOmega も貼る (既定 0 のままだと omega=0 で序盤 NaN →
    EOS 床洗浄後プラトーの指紋 — run_0001 で確認)。

    `gas` (semi-perfect, 2026-08-17): 面積比→M・T・P・ρ・roe を NASA-9 テーブルで
    作る。**roe は TP の内部エネルギー e(T)=h(T)−RT** で構成しないと forge (TP) が
    step 0 で温度を誤再構成する ([[wys-tp-divergence-is-cold-not-multispecies]] の
    IC 不整合と同じ罠)。エンタルピー基準は forge 内蔵 DB (絶対基準, thermoHrefTemp=0)
    と同一の NASA-9 なので整合する。

    `species_Y` (2026-08-17, split_h2o): 多成分 TP の組成 [Y_0, Y_1, ...] を一様に貼る
    (`/VALUE/roY{s}` = ρ Y_s)。forge は roY が無いと Y0=1 と解釈するので、2 種以上では必須。
    """
    g = gamma
    R_gas = cp * (g - 1.0) / g
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn = cc[:, 0] / scale  # 無次元軸位置
        # 物理壁 (A13) はスロートが (x_throat, r_throat) ≠ (0, 1) に動く
        x_thr, r_thr = _throat_of(wall)
        AR = np.maximum(wall.r(xn) / r_thr, 1.0) ** 2
        if gas is not None and getattr(gas, "kind", "cpg") == "semiperfect":
            M = _invert_area_ratio_gas(AR, xn >= x_thr, gas)
            T = np.where(M >= 1.0, gas.T_of_M(np.maximum(M, 1.0)),
                         np.interp(np.minimum(M, 1.0), gas._Mu, gas._Tu))
            # 等エントロピー P/Pt = exp(∫ cp/(R T) dT) (Tt→T)
            P = Pt * _isentropic_pressure_ratio_gas(T, gas)
            R_gas = gas.R
            ro = P / (R_gas * T)
            gam_loc = gas.gamma(T)
            u = M * np.sqrt(gam_loc * R_gas * T)
            e_int = _ic_e_int(gas, T, h_ref_T)                # TP 内部エネルギー (stamp_isentropic_ic_species が同じ式を照合する)
            roe = ro * e_int + 0.5 * ro * u * u
        else:
            M = invert_area_ratio(AR, xn >= x_thr, g)
            fac = 1.0 + 0.5 * (g - 1.0) * M * M
            T = Tt / fac
            P = Pt / fac ** (g / (g - 1.0))
            ro = P / (R_gas * T)
            u = M * np.sqrt(g * R_gas * T)
            roe = P / (g - 1.0) + 0.5 * ro * u * u
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roUz"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        extra = [("roK", ro * k_init), ("roOmega", ro * omega_init)]
        if species_Y is not None:
            Ys = np.asarray(species_Y, dtype=float)
            if abs(float(Ys.sum()) - 1.0) > 1e-6:
                raise ValueError(f"paste_isentropic_ic: species_Y の和 {Ys.sum():.6f} != 1")
            extra += [(f"roY{i}", ro * float(y)) for i, y in enumerate(Ys)]
        for name, val in extra:
            v32 = val.astype(np.float32)
            if f"/VALUE/{name}" in f:
                f[f"/VALUE/{name}"][:] = v32
            else:
                f.create_dataset(f"/VALUE/{name}", data=v32)
    ie = int(np.argmax(xn))
    return {"M_exit_1d": float(M[ie]), "P_exit_1d": float(P[ie]), "T_exit_1d": float(T[ie])}

def _ic_e_int(gas, T, h_ref_T):
    """semi-perfect IC の内部エネルギー e(T)=h(T)−RT [J/kg]。h_ref_T>0 なら forge の thermoHrefTemp と同じ
    sensible datum (h(T_ref)=0)。基準が食い違うと step 0 で T が跳ぶ。"""
    e_int = gas.h_mass(T) - gas.R * np.asarray(T, dtype=float)
    if h_ref_T is not None and h_ref_T > 0.0:
        e_int = e_int - float(np.ravel(gas.h_mass(np.array([float(h_ref_T)])))[0])
    return e_int


def stamp_isentropic_ic_species(h5path, run_dir, gas, h_ref_T, species: list, species_MW: list,
                                species_Y: list | None = None, forge=None) -> str:
    """`paste_isentropic_ic` (semi-perfect) で作った新規初期場に化学種の属性を付ける (plans/active/
    thermophysics-solver-owned-species-db.md §4.3「新規初期場 = IC 生成処理」, #3b)。

    宛先 run_dir を `forge --resolve-species` で解決し、IC に使った datum (`h_ref_T`)・輸送種の順序 (`species`)・
    組成を作るのに使った MW (`species_MW`)・内部エネルギー式 (`_ic_e_int` と gas.R) が宛先の記録と一致したときだけ
    属性 (`species_input_unverified=0`) を付ける。一致しなければ属性を付けずに `forge_species.SpeciesCheckError`。
    宛先を解決できない (--resolve-species を持つ forge が無い) ときも既定で SpeciesCheckError (#3c; その実行だけ
    FORGE_ALLOW_UNVERIFIED_SPECIES=1 なら属性なしで 'unverified')。
    IC の書き込み (と壁速度 0 化など同じ物性での後処理) が終わってから呼ぶ。戻り値: 'verified' | 'cpg' | 'unverified'。"""
    fs = _forge_species()
    Y = list(species_Y) if species_Y is not None else [1.0]
    mixes = [("inflow", Y, gas.R, lambda T: _ic_e_int(gas, T, h_ref_T))]
    return fs.stamp_new_field(h5path, run_dir, list(species), list(species_MW), h_ref_T, mixes, forge=forge,
                              tool="paste_isentropic_ic")


def _invert_area_ratio_gas(AR, supersonic, gas):
    """A/A* → M (ガスモデルのテーブル、超音速/亜音速枝)。"""
    AR = np.asarray(AR, dtype=float)
    M = np.empty_like(AR)
    sup = np.asarray(supersonic, bool)
    # 超音速枝: gas._AR は M 増加で単調増
    M[sup] = np.interp(AR[sup], gas._AR, gas._M)
    # 亜音速枝: gas._ARu は M 増加 (0→1) で単調減 → 反転して補間
    o = np.argsort(gas._ARu)
    M[~sup] = np.interp(AR[~sup], gas._ARu[o], gas._Mu[o])
    return M


def _isentropic_pressure_ratio_gas(T, gas):
    """P/Pt = exp( ∫_{Tt}^{T} cp/(R T') dT' ) を NASA-9 で数値積分 (T ごと)。"""
    T = np.asarray(T, dtype=float)
    Tg = np.linspace(gas.Tt, float(np.min(T)) * 0.98, 4000)
    cp = gas.cp_mass(Tg)
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] / Tg[1:] + cp[:-1] / Tg[:-1])
                                             * np.diff(Tg) / gas.R)])
    return np.exp(np.interp(T, Tg[::-1], integ[::-1]))


def cpg_field_to_tp(src_h5, dst_h5, species_names: list, Y: list, cp_cpg: float, gamma_cpg: float,
                    h_ref_T: float = 298.15, db: dict | None = None) -> dict:
    r"""CPG (thermalMethod 0) の発達場を **多成分 TP (thermalMethod 2) 整合の初期値**に変換する
    (計画 tooling-nozzle-tp-split-h2o-condensation)。

    (ρ, ρu, k, ω, 既存 roY) はそのまま、T は CPG の $e=c_vT$ ($c_v=c_p/\gamma$) から復元し、
    roe を $\rho\,[\sum_s Y_s e_s(T) + e_k]$ (NASA-9, sensible datum $h(T_{ref})=0$ =
    forge `thermoHrefTemp`) で作り直す。roY{s} が無ければ `Y` を一様に貼る。
    `db`: 擬似種を含む speciesDBFile 相当の dict (`mixture_pseudo_species_split` の戻り)。None なら
    内蔵 `SPECIES_NASA9` の名前で引く。case/16 の `gen_tp_ic_wys.py` (N2+H2O 専用) の一般化。
    低温は forge と同じ Tlo クランプ (cp 凍結・h 線形接続) を適用する。
    戻り: 診断 (T 範囲、Y 範囲)。"""
    import shutil
    from ..gas.semiperfect import DESIGN_T_BOUNDS, RU, SPECIES_NASA9, _cp_R_raw, _h_RT_raw, check_design_T
    shutil.copy(src_h5, dst_h5)
    names = [n.upper() for n in species_names]

    def _entry(name):
        if db is not None and name in db:
            e = db[name]
            return (float(e["MW"]), np.asarray(e["nasa9_low"], float), np.asarray(e["nasa9_high"], float),
                    float(e.get("Tlo", 200.0)), float(e.get("Tmid", 1000.0)), float(e.get("Thi", 6000.0)))
        sp = SPECIES_NASA9[name]      # 内蔵種: 共通データの先頭 2 区間 (6000 K 超は評価しない; plan #13-3)
        return (float(sp["MW"]), np.asarray(sp["low"], float), np.asarray(sp["high"], float), *DESIGN_T_BOUNDS)

    def h_mass(name, T):
        MW, lo, hi, Tlo, Tmid, Thi = _entry(name)
        if not (db is not None and name in db):
            check_design_T(T, f"IC {name}")
        R = RU / MW
        T = np.asarray(T, dtype=float)
        Tc = np.clip(T, Tlo, Thi)
        a = np.where(Tc[..., None] < Tmid, lo, hi)
        h = np.array([_h_RT_raw(a[i], Tc[i]) * R * Tc[i] for i in range(len(Tc))]) if T.ndim else 0.0
        cpc = np.array([_cp_R_raw(a[i], Tc[i]) * R for i in range(len(Tc))]) if T.ndim else 0.0
        return h + cpc * (T - Tc), R                    # 域外は境界 cp で線形外挿 (forge と同じ)

    with h5py.File(dst_h5, "r+") as f:
        ro = f["/VALUE/ro"][:].astype(float)
        ux, uy, uz = (f[f"/VALUE/{k}"][:].astype(float) / ro for k in ("roUx", "roUy", "roUz"))
        ek = 0.5 * (ux * ux + uy * uy + uz * uz)
        e_cpg = f["/VALUE/roe"][:].astype(float) / ro - ek
        cv = cp_cpg / gamma_cpg
        T = e_cpg / cv
        Ys = []
        for i, name in enumerate(names):
            if f"/VALUE/roY{i}" in f:
                Ys.append(f[f"/VALUE/roY{i}"][:].astype(float) / ro)
            else:
                Ys.append(np.full_like(ro, float(Y[i])))
        Ys = np.asarray(Ys)
        Ys = Ys / np.maximum(Ys.sum(axis=0), 1e-30)
        e_int = np.zeros_like(ro)
        for i, name in enumerate(names):
            h, R = h_mass(name, T)
            href, _ = h_mass(name, np.full_like(T, h_ref_T))
            e_int += Ys[i] * ((h - href) - R * T)
        roe = ro * (e_int + ek)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for i in range(len(names)):
            v32 = (ro * Ys[i]).astype(np.float32)
            if f"/VALUE/roY{i}" in f:
                f[f"/VALUE/roY{i}"][:] = v32
            else:
                f.create_dataset(f"/VALUE/roY{i}", data=v32)
    return {"T_min": float(T.min()), "T_max": float(T.max()),
            "Y": {n: (float(Ys[i].min()), float(Ys[i].max())) for i, n in enumerate(names)}}


def zero_wall_velocity_ic(h5path, keep_energy: str = "internal") -> int:
    r"""node (median-dual) の cross-mesh IC 後処理: **壁ノード (wall_dist=0) の速度を 0** にし、roe から
    その運動エネルギー分を落とす (内部エネルギー・T を保つ)。

    背景 (2026-08-17, Wyslouzil node): cell 場を最近傍で node へ移すと壁ノードに隣接セルの
    (u≠0, roe=ρ(e+½u²)) がそのまま乗る。`nodeWallDirichlet` が u=0 に固定すると ½u² が内部エネルギー
    に化けて壁 T が +100 K 以上跳ね (u=400 m/s で 108 K)、壁から P/T が Pt/Tt を超えて発散した
    (run_0172–0174: 壁ノードで T 434 K > Tt 286 K, P 116 kPa > Pt 59 kPa)。
    [[node-isothermal-wall-thin-cell-mass-source]] の「原始変数で補間せよ」と同じ趣旨の最小処置。
    戻り: 処理した壁ノード数。"""
    with h5py.File(h5path, "r+") as f:
        wd = f["/VALUE/wall_dist"][:]
        wall = wd <= 0.0
        if not wall.any():
            return 0
        ro = f["/VALUE/ro"][:].astype(float)
        u = [f[f"/VALUE/{k}"][:].astype(float) / ro for k in ("roUx", "roUy", "roUz")]
        ek = 0.5 * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
        roe = f["/VALUE/roe"][:].astype(float)
        roe[wall] = roe[wall] - ro[wall] * ek[wall]
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for k in ("roUx", "roUy", "roUz"):
            v = f[f"/VALUE/{k}"][:]; v[wall] = 0.0; f[f"/VALUE/{k}"][:] = v
        if "/VALUE/roK" in f:
            v = f["/VALUE/roK"][:]; v[wall] = 0.0; f["/VALUE/roK"][:] = v
    return int(wall.sum())
