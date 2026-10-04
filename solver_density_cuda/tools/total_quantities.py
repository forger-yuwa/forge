#!/usr/bin/env python3
"""全温 T0・全圧 P0 を res_*.h5 の全エンタルピー `VALUE/h0` から作る (AGENTS.md「出力と後処理の原則」)。

ソルバが書く h0 = e + p/ρ + u²/2 (+ k は sstEnergyIncludesK のときだけ; 属性 h0_includes_k) を逆算するだけで、
スクリプト側で T + u²/2c_p を組まない (k を含めるかの判断をソルバ側に閉じる)。

  usage: total_quantities.py RUN_DIR [--res res_N.h5] [--write] [--Tt 286.65] [--species-source auto|record|speciesDBFile]
    --write : VALUE/T0, VALUE/P0 を res に追記 (属性 includes_k / method)。既にあれば上書き。
    --Tt    : 比較用の入口全温 (省略時は bcondConfig の inlet の Tt を探す)
    --species-source : TP の熱物性の読み元 (既定 auto; speciesDBFile は旧経路との照合用)

gas model は solverConfig.yaml の physProp から判定:
  thermalMethod 0 (CPG)   : T0 = h0/c_p, P0 = P (T0/T)^{γ/(γ-1)}
  thermalMethod 2 (TP)    : h_mix(T0) = h0 を Newton で逆算 (NASA-9 と thermoHrefTemp の datum は forge_species.run_thermo:
                            res の属性が指すソルバの解決済み記録 > run の記録 > 従来の speciesDBFile > forge --resolve-species),
                            P0 = P exp((s°(T0) − s°(T))/R_mix) (凍結組成)
  凝縮 (g>0)               : 凍結組成の気相逆算のみ (潜熱項は未対応 → 警告)。
Python API: total_state(run_dir, res_path, species_source="auto") -> dict(T0, P0, includes_k, method, species_source)
"""
import argparse, glob, os, sys
import numpy as np, h5py, yaml

RU = 8.314462618   # ソルバの THERMO_RU (cuda_forge/thermo_d.cuh) と同じ値にする (CODATA 8.31446261815324 との差 1.8e-11 が h の 1e-6 J/kg 差になる)


def _nasa9(a, T):
    """NASA-9: (cp/R, h/(RT), s/R)"""
    T = np.asarray(T, dtype=np.float64)
    h_RT = (-a[0] / T**2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3
            + a[5] * T**3 / 4 + a[6] * T**4 / 5 + a[7] / T)
    cp_R = a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4
    s_R = (-a[0] / (2 * T**2) - a[1] / T + a[2] * np.log(T) + a[3] * T + a[4] * T**2 / 2
           + a[5] * T**3 / 3 + a[6] * T**4 / 4 + a[8])
    return cp_R, h_RT, s_R


class _TPGas:
    """種エントリ {MW, Tlo, Tmid, Thi, nasa9_low, nasa9_high} または区間可変 {MW, Tbounds, nasa9_intervals} (記録 / speciesDBFile;
    forge_species.run_thermo) の凍結組成混合。質量基準の h, cp, s° (datum: thermoHrefTemp)。
    範囲外の扱いはソルバ (cuda_forge/thermo_d.cuh thermo_cp_molar / thermo_h_molar / thermo_s0_mass) と同じ:
    種ごとの Tlo/Thi の外では cp を端の値で固定し、h は線形外挿 h(T)=h(Tb)+cp(Tb)(T−Tb)、s° は s°(Tb)+cp(Tb) ln(T/Tb)。
    係数は区間 k = Tb[k] <= T < Tb[k+1] (区切りちょうどは上の区間; 2 区間では T<Tmid で low、それ以外 high = thermo_pick_coeffs)。
    codex 2026-09-16 result-3 M2、区間可変は plan thermophysics-solver-owned-species-db #13-1。"""
    def __init__(self, db, names, Tref):
        self.sp = [db[n] for n in names]; self.R = [RU / s["MW"] for s in self.sp]; self.Tref = Tref
        self.href = [self._h1(s, np.array([Tref]))[0] if Tref > 0 else 0.0 for s in self.sp]

    @staticmethod
    def _intervals(s):
        """(境界 [Tlo, 区切り..., Thi], 係数 ndarray (nInt, 9))。2 区間の書式は区切りの既定値 200/1000/6000 K (C++ 外部 DB 読込と同じ)。"""
        if s.get("nasa9_intervals") is not None:
            return [float(x) for x in s["Tbounds"]], np.asarray(s["nasa9_intervals"], float)
        Tb = [float(s.get("Tlo", 200.0)), float(s.get("Tmid", 1000.0)), float(s.get("Thi", 6000.0))]
        return Tb, np.asarray([s["nasa9_low"], s["nasa9_high"]], float)

    @classmethod
    def _bounds(cls, s):
        Tb = cls._intervals(s)[0]
        return Tb[0], Tb[-1]

    def _coef(self, s, T):
        Tb, co = self._intervals(s)
        if len(co) == 2:
            return np.where((T < Tb[1])[:, None], co[0], co[1])
        # 区切り Tb[1..n-1] のうち !(T < 区切り) の数 = 区間番号 (NaN は最後の区間; C++ thermo_interval と同じ)
        k = np.zeros(np.shape(T), dtype=int)
        for b in Tb[1:-1]:
            k += ~(T < b)
        return co[k]

    def _raw(self, s, Tc):
        """クランプ済み温度での (cp, h, s°) [質量基準] (thermo_*_clamped)。"""
        a = self._coef(s, Tc); R = RU / s["MW"]
        cp_R, h_RT, s_R = _nasa9(a.T, Tc)
        return cp_R * R, h_RT * R * Tc, s_R * R

    def _props(self, s, T):
        """範囲クランプ + 外挿込みの (cp, h, s°) [質量基準] (thermo_cp_molar / thermo_h_molar / thermo_s0_mass と同式)。"""
        T = np.asarray(T, dtype=np.float64)
        Tlo, Thi = self._bounds(s)
        Tc = np.clip(T, Tlo, Thi)
        cp, h, s0 = self._raw(s, Tc)
        out = (T < Tlo) | (T > Thi)
        if np.any(out):
            # 端の値 cp(Tb), h(Tb), s°(Tb) から線形 (h) / 対数 (s°) 外挿。cp は端の値で一定。
            h = np.where(out, h + cp * (T - Tc), h)
            s0 = np.where(out, s0 + cp * np.log(np.maximum(T, 1e-300) / Tc), s0)
        return cp, h, s0

    def _h1(self, s, T):
        return self._props(s, T)[1]

    def h(self, Y, T):
        return sum(Y[i] * (self._props(s, T)[1] - self.href[i]) for i, s in enumerate(self.sp))

    def cp(self, Y, T):
        return sum(Y[i] * self._props(s, T)[0] for i, s in enumerate(self.sp))

    def s0(self, Y, T):
        return sum(Y[i] * self._props(s, T)[2] for i, s in enumerate(self.sp))

    def Rmix(self, Y):
        return sum(Y[i] * self.R[i] for i in range(len(self.sp)))


def total_state(run_dir, res_path=None, species_source="auto"):
    run_dir = os.path.abspath(run_dir)
    if res_path is None:
        res_path = sorted(glob.glob(os.path.join(run_dir, "res_[0-9]*.h5")), key=lambda s: int(s.split("_")[-1][:-3]))[-1]
    cfg = yaml.safe_load(open(os.path.join(run_dir, "solverConfig.yaml")))
    pp = cfg["physProp"]; tm = int(pp.get("thermalMethod", 0))
    with h5py.File(res_path, "r") as f:
        if "VALUE/h0" not in f:
            raise SystemExit(f"{res_path}: VALUE/h0 が無い (output.level>=1 の新しいバイナリで出力した res が必要)")
        h0 = f["VALUE/h0"][:].astype(np.float64); inc_k = int(f["VALUE/h0"].attrs.get("h0_includes_k", 0))
        T = f["VALUE/T"][:].astype(np.float64); P = f["VALUE/P"][:].astype(np.float64)
        n = len(T)
        Y = None
        g = f["VALUE/g_0"][:].astype(np.float64) if "VALUE/g_0" in f else None
    warn = []
    if g is not None and np.nanmax(g) > 1e-9:
        warn.append(f"凝縮 g_max={np.nanmax(g):.3e}: 潜熱項は未対応、気相 (凍結組成) の逆算として扱う")
    if tm == 0:
        cp = float(pp["cp"]); ga = float(pp["gamma"])
        T0 = h0 / cp
        P0 = P * (T0 / T) ** (ga / (ga - 1.0)); method = "CPG"
    elif tm == 2:
        # 熱物性は forge_species.run_thermo (plan thermophysics-solver-owned-species-db #8): res の属性が指す記録を優先し、
        # 記録の無い旧 run は従来どおり speciesDBFile から読む (species_db.yaml を前提にしない)
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import forge_species as fsp
        th = fsp.run_thermo(run_dir, res_path=res_path, source=species_source)
        names = th["names"]
        with h5py.File(res_path, "r") as f:
            if len(names) == 1 and "VALUE/Y0" not in f:      # 単一種は Y を書かない → Y=1
                Y = [np.ones(n)]
            else:
                Y = [f[f"VALUE/Y{i}"][:].astype(np.float64) for i in range(len(names))]
        gas = fsp.thermo_gas(th)
        # Newton: h(T0) = h0 (初期値 T)
        T0 = T.copy()
        for _ in range(30):
            r = gas.h(Y, T0) - h0; d = gas.cp(Y, T0)
            dT = -r / np.maximum(d, 1.0); T0 = np.maximum(T0 + dT, 1.0)
            if np.max(np.abs(dT)) < 1e-6: break
        Rm = gas.Rmix(Y)
        P0 = P * np.exp((gas.s0(Y, T0) - gas.s0(Y, T)) / Rm); method = "TP-NASA9(frozen)"
    else:
        raise SystemExit(f"thermalMethod {tm} は未対応")
    src = None if tm != 2 else f"{th['source']} ({th['how']}{': ' + os.path.basename(th['path']) if th['path'] else ''})"
    return {"T0": T0, "P0": P0, "includes_k": inc_k, "method": method, "res": res_path, "warn": warn, "T": T, "P": P,
            "species_source": src}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--res", default=None)
    ap.add_argument("--write", action="store_true"); ap.add_argument("--Tt", type=float, default=None)
    ap.add_argument("--species-source", default="auto", choices=("auto", "record", "speciesDBFile"),
                    help="TP の熱物性の読み元 (既定 auto = 記録優先; speciesDBFile は旧経路との照合用)")
    a = ap.parse_args()
    st = total_state(a.run_dir, a.res, species_source=a.species_source)
    Tt = a.Tt
    if Tt is None:
        try:
            bc = yaml.safe_load(open(os.path.join(a.run_dir, "bcondConfig.yaml")))
            for v in bc.values():
                if isinstance(v, dict) and str(v.get("kind", "")).startswith("inlet") and "Tt" in (v.get("floats") or {}):
                    Tt = float(v["floats"]["Tt"]); break
        except Exception:
            pass
    T0, P0 = st["T0"], st["P0"]
    print(f"{st['res']}: method={st['method']} h0_includes_k={st['includes_k']}"
          + (f" species={st['species_source']}" if st["species_source"] else ""))
    for w in st["warn"]: print("  WARN:", w)
    print(f"  T0: min {T0.min():.3f} / median {np.median(T0):.3f} / max {T0.max():.3f} K" + (f"  (Tt={Tt}: max−Tt {T0.max()-Tt:+.3f} K, n(T0>Tt+1) {(T0>Tt+1).sum()})" if Tt else ""))
    print(f"  P0: min {P0.min():.1f} / max {P0.max():.1f} Pa")
    if a.write:
        with h5py.File(st["res"], "r+") as f:
            for nm, arr in (("T0", T0), ("P0", P0)):
                if f"VALUE/{nm}" in f: del f[f"VALUE/{nm}"]
                ds = f.create_dataset(f"VALUE/{nm}", data=arr.astype(np.float32))
                ds.attrs["includes_k"] = st["includes_k"]; ds.attrs["method"] = st["method"]
        print("  written VALUE/T0, VALUE/P0")


if __name__ == "__main__":
    main()
