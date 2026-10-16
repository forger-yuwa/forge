"""`contur_calibrate`・`wt_qualify` の試験用の合成の run (plan tooling-nozzle-card-local-mesh §5.1 #6a・#6b)。

`make_run(dir, ...)` は NS の run の較正に要るファイルだけを書く: prepare_info.json (scale_m・initializer・x_E・physical_wall.domain・Md)、
delta_r_initial.json・delta_r_initial.csv (積分法の関数の値)、δ_E の表 (delta_r_next.csv の形式)。δ は実寸 [m] で
  δ_C(x) = contur_m(x)、入れた δ = (c₀ + c₁·x)·δ_C + δ₀ (今の係数)、δ_E = (t₀ + t₁·x)·δ_C + t_d + bump(x) (目標の係数 + 形の残差)
とし、表は r_t 単位 (÷ S) で書く。`v2=True` なら停止規則 v2 の状態 (STOP_RULE_STATE.json の delta_E と stop_rule_v2/delta_E_<k>/) に置く。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

S_DEF = 0.1
X_IN_RT, X_OUT_RT = -10.0, 80.0


def x_columns(n_up: int = 300, n_down: int = 1500, clustered: bool = True, with_zero: bool = False) -> np.ndarray:
    """NS の列の x [r_t] (スロートの前後に寄せる)。既定は 0 ちょうどの列を置かない (実際の表と同じ; スロートの値は両隣の列の直線補間)。
    with_zero=True なら 0 の列を置く (スロートの補間が厳密になり、形の残差の無い δ_E から係数を丸めの差で取り戻せる)。"""
    up = np.linspace(X_IN_RT, -1e-3, n_up)
    s = np.linspace(0.0, 1.0, n_down)
    x0 = 0.0 if with_zero else 1.5e-3
    down = x0 + (X_OUT_RT - x0) * (s ** 2 if clustered else s)
    return np.r_[up, down]


def contur_m(x_m) -> np.ndarray:
    """係数前の CONTUR の δ の代わり [m] (上流は小さく、下流で単調に増える滑らかな形)。"""
    x = np.asarray(x_m, dtype=float)
    h = 0.5 * (1.0 + np.tanh(x / 0.3))
    return 2e-4 + 3e-3 * h + 4.5e-3 * x * h


def bump_m(x_m, amp: float = 3e-4, kind: str = "hump") -> np.ndarray:
    """δ_E に足す形の残差 [m] (両端で 0)。hump = x ≥ 0 の中ほどの山 (3 係数で一部を吸収できる)、wave = 3 周期の波
    (3 係数ではほとんど減らせない)。"""
    x = np.asarray(x_m, dtype=float)
    L = X_OUT_RT * S_DEF
    s = np.clip(x / L, 0.0, 1.0)
    if kind == "wave":
        return np.where(x > 0.0, amp * np.sin(6.0 * np.pi * s), 0.0)
    return np.where(x > 0.0, amp * np.sin(np.pi * s) ** 2, 0.0)


def _sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_delta_next(path, x_rt, delta_in_rt, delta_e_rt, held=None):
    held = np.zeros(len(x_rt)) if held is None else np.asarray(held, dtype=float)
    np.savetxt(path, np.c_[x_rt, delta_e_rt, delta_in_rt, delta_e_rt, held], delimiter=",", comments="",
               header="x_rt,delta_target,delta_in_prev,delta_E,held (synthetic)")


def make_run(d, *, S: float = S_DEF, coefs=(0.85, 0.0, -3.2e-4), target=(0.68, 0.021, -3.06e-4), bump: float = 0.0, bump_kind: str = "hump",
             x_rt=None, v2: bool = False, step: int = 46202, held_at=None, delta_r_csv: bool = False, initializer_in_prepare: bool = True,
             mismatch_json: bool = False) -> Path:
    """合成の run を d に作って返す。coefs = 今の係数 (c₀, c₁ [1/m], δ₀ [m]) か None (pass 0 = 係数なし)。"""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    x_rt = x_columns() if x_rt is None else np.asarray(x_rt, dtype=float)
    x_m = x_rt * S
    C = contur_m(x_m)
    c0, c1, d0 = coefs if coefs is not None else (1.0, 0.0, 0.0)
    din = (c0 + c1 * x_m) * C + d0
    t0, t1, td = target
    E = (t0 + t1 * x_m) * C + td + bump_m(x_m, bump, bump_kind)
    held = None
    if held_at is not None:
        held = np.zeros(len(x_rt))
        held[int(np.argmin(np.abs(x_rt - held_at)))] = 1
    ini = {"model": "contur_momentum_integral", "closure": "contur", "cf_scale": 1.0, "n_scale": 1.0, "rt_m": S,
           "delta_r_throat": float(np.interp(0.0, x_m, din) / S), "delta_r_exit": float(din[-1] / S)}
    if coefs is not None:
        ini.update(delta_scale=c0, delta_offset_m=d0)
        if c1 != 0.0:
            ini["delta_scale_slope_per_m"] = c1
        ini["two_coef"] = {"delta_scale": c0, "delta_offset_m": d0, "scale_m": S, **({"delta_scale_slope_per_m": c1} if c1 != 0.0 else {}),
                           "delta_r_throat_contur": float(np.interp(0.0, x_m, C) / S), "delta_r_exit_contur": float(C[-1] / S)}
    info = {"scale_m": S, "Md": 5.5, "x_E": 40.0, "physical_wall": {"domain": [float(x_rt[0]), float(x_rt[-1])]},
            "delta_r_csv": "some.csv" if delta_r_csv else None,
            "dstar_source": "integral_bl:contur_momentum_integral (radial)", "initializer": ini if initializer_in_prepare else None}
    (d / "prepare_info.json").write_text(json.dumps(info, indent=1))
    ini_js = dict(ini)
    if mismatch_json:
        ini_js["delta_scale"] = float(ini_js.get("delta_scale", 1.0)) * 1.001
    (d / "delta_r_initial.json").write_text(json.dumps(ini_js, indent=1))
    xf = np.linspace(x_rt[0], x_rt[-1], 1500)
    xfm = xf * S
    dfun = (c0 + c1 * xfm) * contur_m(xfm) + d0
    np.savetxt(d / "delta_r_initial.csv", np.c_[xf, dfun / S, np.zeros((len(xf), 8))], delimiter=",", comments="",
               header="x_rt,delta_r_rt,dstar_n_rt,theta_m,H,N,Cf,M_e,Tw,Taw")
    if v2:
        od = d / "stop_rule_v2" / f"delta_E_{step}"
        od.mkdir(parents=True, exist_ok=True)
        write_delta_next(od / "delta_r_next.csv", x_rt, din / S, E / S, held)
        st = {"version": 2, "method": "v2", "status": "REACHED", "stop_step": step, "reeval": {"verdict": "SETTLED"},
              "delta_E": {"status": "ok", "step": step, "dir": f"/elsewhere/{d.name}/stop_rule_v2/delta_E_{step}",
                          "res_sha256": "x" * 64, "delta_r_next_sha256": _sha(od / "delta_r_next.csv")}}
        (d / "STOP_RULE_STATE.json").write_text(json.dumps(st, indent=1))
    else:
        write_delta_next(d / "delta_r_next.csv", x_rt, din / S, E / S, held)
    return d
