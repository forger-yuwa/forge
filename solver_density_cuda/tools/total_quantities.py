#!/usr/bin/env python3
"""全温 T0・全圧 P0 を res_*.h5 の全エンタルピー `VALUE/h0` から作る (AGENTS.md「出力と後処理の原則」)。

ソルバが書く h0 = e + p/ρ + u²/2 (+ k は sstEnergyIncludesK のときだけ; 属性 h0_includes_k) を逆算するだけで、
スクリプト側で T + u²/2c_p を組まない (k を含めるかの判断をソルバ側に閉じる)。

  usage: total_quantities.py RUN_DIR [--res res_N.h5] [--write] [--Tt 286.65]
    --write : VALUE/T0, VALUE/P0 を res に追記 (属性 includes_k / method)。既にあれば上書き。
    --Tt    : 比較用の入口全温 (省略時は bcondConfig の inlet の Tt を探す)

gas model は solverConfig.yaml の physProp から判定:
  thermalMethod 0 (CPG)   : T0 = h0/c_p, P0 = P (T0/T)^{γ/(γ-1)}
  thermalMethod 2 (TP)    : h_mix(T0) = h0 を Newton で逆算 (species_db.yaml の NASA-9, thermoHrefTemp の datum),
                            P0 = P exp((s°(T0) − s°(T))/R_mix) (凍結組成)
  凝縮 (g>0)               : 凍結組成の気相逆算のみ (潜熱項は未対応 → 警告)。
Python API: total_state(run_dir, res_path) -> dict(T0, P0, includes_k, method)
"""
import argparse, glob, os, sys
import numpy as np, h5py, yaml

RU = 8.31446261815324


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
    """species_db.yaml (NASA-9) の凍結組成混合。質量基準の h, cp, s° (datum: thermoHrefTemp)。"""
    def __init__(self, db, names, Tref):
        self.sp = [db[n] for n in names]; self.R = [RU / s["MW"] for s in self.sp]; self.Tref = Tref
        self.href = [self._h1(s, np.array([Tref]))[0] if Tref > 0 else 0.0 for s in self.sp]

    def _coef(self, s, T):
        lo, hi = np.asarray(s["nasa9_low"], float), np.asarray(s["nasa9_high"], float)
        return np.where((T < s["Tmid"])[:, None], lo, hi)

    def _h1(self, s, T):
        a = self._coef(s, T); R = RU / s["MW"]
        _, h_RT, _ = _nasa9(a.T, T); return h_RT * R * T

    def _props(self, s, T):
        a = self._coef(s, T); R = RU / s["MW"]
        cp_R, h_RT, s_R = _nasa9(a.T, T); return cp_R * R, h_RT * R * T, s_R * R

    def h(self, Y, T):
        return sum(Y[i] * (self._props(s, T)[1] - self.href[i]) for i, s in enumerate(self.sp))

    def cp(self, Y, T):
        return sum(Y[i] * self._props(s, T)[0] for i, s in enumerate(self.sp))

    def s0(self, Y, T):
        return sum(Y[i] * self._props(s, T)[2] for i, s in enumerate(self.sp))

    def Rmix(self, Y):
        return sum(Y[i] * self.R[i] for i in range(len(self.sp)))


def total_state(run_dir, res_path=None):
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
        if tm == 2:
            names = pp.get("species", None)
            if names:
                if len(names) == 1 and "VALUE/Y0" not in f:      # 単一擬似種は Y を書かない → Y=1
                    Y = [np.ones(n)]
                else:
                    Y = [f[f"VALUE/Y{i}"][:].astype(np.float64) for i in range(len(names))]
            else:
                names = None
        g = f["VALUE/g_0"][:].astype(np.float64) if "VALUE/g_0" in f else None
    warn = []
    if g is not None and np.nanmax(g) > 1e-9:
        warn.append(f"凝縮 g_max={np.nanmax(g):.3e}: 潜熱項は未対応、気相 (凍結組成) の逆算として扱う")
    if tm == 0:
        cp = float(pp["cp"]); ga = float(pp["gamma"])
        T0 = h0 / cp
        P0 = P * (T0 / T) ** (ga / (ga - 1.0)); method = "CPG"
    elif tm == 2:
        db = yaml.safe_load(open(os.path.join(run_dir, pp.get("speciesDBFile", "species_db.yaml"))))
        if names is None:
            names = list(db.keys())[:1]; Y = [np.ones(n)]
        gas = _TPGas(db, names, float(pp.get("thermoHrefTemp", 0.0)))
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
    return {"T0": T0, "P0": P0, "includes_k": inc_k, "method": method, "res": res_path, "warn": warn, "T": T, "P": P}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--res", default=None)
    ap.add_argument("--write", action="store_true"); ap.add_argument("--Tt", type=float, default=None)
    a = ap.parse_args()
    st = total_state(a.run_dir, a.res)
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
    print(f"{st['res']}: method={st['method']} h0_includes_k={st['includes_k']}")
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
