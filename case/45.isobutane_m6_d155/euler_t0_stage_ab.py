"""plan verification-case45-euler-total-enthalpy §6 E1 (2026-10-07 登録): ISEN の起動前の prep と soft 段の後の状態を、
同じ保存量の復元で比べ、全温の超過が IC の生成で生じるか、起動後に生じるかを判別する (CFD 0 step、ファイルには書かない)。

A = 起動前の prep `_prep_moc_v5_mocG1_isen/nozzle.h5`、B = soft 段の後 `run_0153_euler_icdep_mocG1_isen/nozzle.h5`、
C = 本段の初期場 `run_0153_euler_icdep_mocG1_isen/res_0.h5` (B と保存量がビット一致かを確かめる)。
復元は V5c の経路 B (`moc_v5c_thermo_ab.path_b`) と同じ式: Y は負を 0 にして和で割り、e = roe/ro − |u|²/2 から
h_mix(T) − R_mix T = e を Newton で解き、h0 = (roe + P)/ro (Euler なので k は足さない)、h(T0) = h0 で全温。
判定 (§6 E1): A が全節点で ±1 K 以内かつ B が +315 K を再現 → 「IC の生成だけで説明する」説を退ける / A にも +100 K 以上 →
「起動後に初めて生じた」説を退ける / その中間 → 寄与の大きさを記録し、単独の原因とは認定しない。
usage: python3 euler_t0_stage_ab.py [case_dir] → _band_ab/euler_t0_stage_ab.json"""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import moc_v5c_thermo_ab as V  # noqa: E402

PLAN_REG_COMMIT = "7fa6295e"
REPORTED_B_MAX_T0 = 1915.38   # V5c: run_0153 res_0.h5 の最大全温 (= +315 K)
STATES = {"A_prep": ("_prep_moc_v5_mocG1_isen", "nozzle.h5"),
          "B_after_soft": ("run_0153_euler_icdep_mocG1_isen", "nozzle.h5"),
          "C_main_res0": ("run_0153_euler_icdep_mocG1_isen", "res_0.h5")}
CONS = ("ro", "roUx", "roUy", "roUz", "roe", "roY0", "roY1")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def recon(run_dir: Path, h5: Path) -> dict:
    th = V.FSP.run_thermo(str(run_dir), res_path=str(h5), source="auto")
    gas = V.FSP.thermo_gas(th)
    names = list(th["names"])
    with h5py.File(h5, "r") as f:
        g = {k: f[f"VALUE/{k}"][:].astype(np.float64) for k in ("ro", "roUx", "roUy", "roUz", "roe")}
        roY = [f[f"VALUE/roY{i}"][:].astype(np.float64) for i in range(len(names))]
    ro = g["ro"]
    Yc = [np.where(r / ro < 0.0, 0.0, r / ro) for r in roY]
    ys = np.sum(Yc, axis=0)
    Y = [y / np.maximum(ys, 1e-30) for y in Yc]
    ek = 0.5 * ((g["roUx"] / ro) ** 2 + (g["roUy"] / ro) ** 2 + (g["roUz"] / ro) ** 2)
    e = g["roe"] / ro - ek
    Rm = gas.Rmix(Y)
    T, sT, _ = V.newton_T(lambda t: gas.h(Y, t) - Rm * t, lambda t: np.maximum(gas.cp(Y, t) - Rm, 1e-2 * Rm), e,
                          np.full(ro.shape, V.NEWTON_T_INIT))
    P = ro * Rm * T
    h0 = (g["roe"] + P) / ro
    T0, sT0, _ = V.newton_T(lambda t: gas.h(Y, t), lambda t: np.maximum(gas.cp(Y, t), 1.0), h0, T)
    return {"T0": T0, "h0": h0, "Y": Y, "gas": gas, "src": f"{th['source']} ({th['how']})",
            "newton_ok": bool(np.all(sT < V.NEWTON_TOL) and np.all(sT0 < V.NEWTON_TOL))}


def main(argv=None):
    a = sys.argv[1:] if argv is None else argv
    case = Path(a[0]).resolve() if a else C
    Tt = V.read_Tt(case / STATES["B_after_soft"][0])
    out = {"plan": "plans/active/verification-case45-euler-total-enthalpy.md §6 E1", "plan_reg_commit": PLAN_REG_COMMIT,
           "evaluator_sha256": _sha(Path(__file__)), "Tt": Tt, "states": {}}
    R = {}
    for tag, (rd, fn) in STATES.items():
        run_dir, h5 = case / rd, case / rd / fn
        r = recon(run_dir, h5)
        X, Rr, ni, nj = V.load_geometry(run_dir)
        d = r["T0"] - Tt
        i_max = int(np.argmax(np.abs(d)))
        hmix_Tt = r["gas"].h(r["Y"], np.full_like(r["T0"], Tt))
        out["states"][tag] = {"file": f"{rd}/{fn}", "sha256": _sha(h5), "thermo_source": r["src"], "newton_ok": r["newton_ok"],
                              "T0_max": float(r["T0"].max()), "T0_min": float(r["T0"].min()),
                              "max_abs_T0_minus_Tt": float(np.abs(d).max()), "loc_max_abs": V._loc(i_max, X, Rr),
                              "n_above_Tt_plus_1K": int(np.sum(d > 1.0)), "n_below_Tt_minus_1K": int(np.sum(d < -1.0)),
                              "n_above_Tt_plus_100K": int(np.sum(d > 100.0)), "n_below_Tt_minus_100K": int(np.sum(d < -100.0)),
                              "h0_minus_hmixTt": V._stats(r["h0"] - hmix_Tt), "layers": V.layer_stats(r["T0"], X, Tt),
                              "n_nodes": int(r["T0"].size)}
        R[tag] = r
    # B と C の保存量のビット一致
    with h5py.File(case / STATES["B_after_soft"][0] / "nozzle.h5", "r") as fb, \
            h5py.File(case / STATES["C_main_res0"][0] / "res_0.h5", "r") as fc:
        out["B_equals_C_conserved"] = {k: bool(np.array_equal(fb[f"VALUE/{k}"][:], fc[f"VALUE/{k}"][:])) for k in CONS}
    # IC の記録との照合 (A)
    icm = case / STATES["A_prep"][0] / "IC_MAP.json"
    if icm.is_file():
        rec = json.loads(icm.read_text())
        out["A_ic_record"] = {k: rec.get(k) for k in ("VERDICT", "mode") if k in rec}
        shas = [v for k, v in rec.items() if "sha" in k.lower() and isinstance(v, str)]
        out["A_ic_record"]["sha256_fields"] = {k: v for k, v in rec.items() if "sha" in k.lower() and isinstance(v, str)}
        out["A_sha256_in_record"] = out["states"]["A_prep"]["sha256"] in shas
    A, B = out["states"]["A_prep"], out["states"]["B_after_soft"]
    a_ok = A["max_abs_T0_minus_Tt"] <= 1.0
    b_rep = abs(B["T0_max"] - REPORTED_B_MAX_T0) <= 1.0
    if not (A["newton_ok"] and B["newton_ok"]):
        v = "判定不能: Newton が収束しない節点がある"
    elif a_ok and b_rep:
        v = "退ける: 「数百 K の超過を IC の生成だけで説明する」説 → 起動後の処理 (soft 段の計算・段の引き継ぎ・本段の起動) を対象にする"
    elif A["n_above_Tt_plus_100K"] > 0:
        v = "退ける: 「起動後に初めて生じた」説 → IC の整合を先に扱う"
    else:
        v = "中間: 寄与の大きさを記録し、単独の原因とは認定しない"
    out["VERDICT"] = v
    out["VERDICT_detail"] = {"A_max_abs_T0_minus_Tt": A["max_abs_T0_minus_Tt"], "B_T0_max": B["T0_max"],
                             "B_reproduces_reported": b_rep, "reported_B_T0_max": REPORTED_B_MAX_T0}
    (case / "_band_ab").mkdir(exist_ok=True)
    (case / "_band_ab/euler_t0_stage_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for tag, s in out["states"].items():
        print(f"{tag:14s} T0 {s['T0_min']:.2f}〜{s['T0_max']:.2f} K  max|T0−Tt| {s['max_abs_T0_minus_Tt']:.3f} K at "
              f"x {s['loc_max_abs']['x_rt']:.3f} r_t j_w {s['loc_max_abs']['j_from_wall']}  "
              f">+1K {s['n_above_Tt_plus_1K']}  <−1K {s['n_below_Tt_minus_1K']}  >+100K {s['n_above_Tt_plus_100K']}  newton {s['newton_ok']}")
    print("B == C (保存量):", out["B_equals_C_conserved"])
    print("A の IC 記録:", out.get("A_ic_record"), "sha 一致:", out.get("A_sha256_in_record"))
    print("VERDICT:", v)


if __name__ == "__main__":
    main()
