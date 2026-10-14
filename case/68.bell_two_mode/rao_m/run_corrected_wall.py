#!/usr/bin/env python3
r"""③ベル (M) L100: δ* 補正壁 (plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-4) の事前登録どおり)。

段 (サブコマンド):

- `design` (約 25 分、手元は 1 本・nice): 補正後の非粘性壁 — CFD ピン版 (M) L100 と同じ初期線 (凍結源 run_0015 の res_30000、
  `run_pin_line.py` と同じ `moc_initial_line` n 321)・同じ L (`cone15_length(50)`)・円弧 kernel で、ε だけ
  ε_inv = (√50 − δ_r,e)² に変えて `design_rao` を回す。δ_r,e = kδ·C(x_E) (kδ は `kdelta_l100_calibration.json`、
  C = contur_v2 の `contur_l100_euler_run0018.npz` の `delta_r_ann`、凍結した一発補正)。出力 `rao_corr_parts/<tag>.json`・
  `<tag>_wall.npz` (`run_rao_t6.py` と同じ形式) とキャンペーンの要約 `rao_corr.json` (`setup_by_tag` に ε_inv)。
- `offset` (数秒): 物理壁 = 非粘性壁 + δ(x) (固定 x の半径方向) の δ(x) を 2 通り作る (`geometry.wall.DisplacementOffset`):
  V-ext (較正式を延長、x ∈ [0, X_EXT] の 5 次の滑らかな立ち上げでスロートの r・r′・r″ を保つ) と V-zero (x < X0 は 0、
  [X0, X1] で同じ立ち上げ)。kδ·C は x ∈ [0, x_E] で 5 次の平滑化スプライン。出力 `rao_corr_parts/offset_<variant>.npz`・`.json`。
- `check` (数秒): 補正壁の問題 YAML 2 本と狙いの Euler の YAML を runner の壁の口 (`_build_wall_rao`) で組み、
  壁の検査 (`validate`)・スロートの接合・出口半径 √50・δ の当てはめの残差・壁の r′・曲率の範囲を記録 (`rao_corr_check.json`)。

実行 (手元 PC はユーザの作業機): `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_corrected_wall.py design|offset|check`
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(HERE))
import run_pin_line as pl  # noqa: E402  (凍結源・ガス・初期線・記録の形式を共用)
import run_rao_t6 as t6  # noqa: E402
from forge_design.feedback.cfd_initial_line import moc_initial_line  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402

PARTS = HERE / "rao_corr_parts"
OUT = HERE / "rao_corr.json"
CAL = HERE / "kdelta_l100_calibration.json"
CONTUR = HERE / "contur_l100_euler_run0018.npz"
EPS_PHYS = pl.EPS
L100 = pl.L100
X_EXT = 0.41      # V-ext の立ち上げの幅: x < 0.41 は M3 で δ_r が時間的に定常でない区間 (§5.1 #7)
X0, X1 = 1.17, 2.0  # V-zero: 較正の区間の始まり x 1.17 から立ち上げ、2.0 で較正式に一致 (§5.1 #7 (7b-4))
VARIANTS = {"vext": (0.0, X_EXT), "vzero": (X0, X1)}
PROBLEMS = {"vext": CASE / "problem_bell_tp_m_ns_l100_corr_vext.yaml",
            "vzero": CASE / "problem_bell_tp_m_ns_l100_corr_vzero.yaml",
            "target": CASE / "problem_bell_tp_m_euler_l100_inv.yaml",
            "target_fine": CASE / "problem_bell_tp_m_euler_l100_inv_fine.yaml"}


def frozen_delta() -> dict:
    cal = json.loads(CAL.read_text())
    co = np.load(CONTUR, allow_pickle=True)
    x, C = np.asarray(co["x"], float), np.asarray(co["delta_r_ann"], float)
    kd = float(cal["k_delta"])
    de = kd * float(np.interp(L100, x, C))
    return {"k_delta": kd, "x": x, "C": C, "x_E": float(x[-1]), "C_E": float(C[-1]), "delta_E": de,
            "eps_inv": (np.sqrt(EPS_PHYS) - de) ** 2,
            "sha256_16": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in (CAL, CONTUR)}}


def tag_inv(eps_inv: float) -> str:
    return f"L100_inv_eps{eps_inv:.4f}".replace(".", "p") + "_pin_n321"


def cmd_design() -> None:
    t0 = time.time()
    fd = frozen_delta()
    if abs(fd["x_E"] - L100) > 1e-6:
        raise SystemExit(f"contur の終点 {fd['x_E']} が L100 {L100} と違う")
    ig = pl.gas()
    th, match = pl.pinned_throat(pl.FREEZE, ig)
    line = moc_initial_line(th, ig.moc, pl.N_LINE)
    tag = tag_inv(fd["eps_inv"])
    spec = {"tag": tag, "L_frac": 1.0, "variant": "pin_inv_corr"}
    res = ra.design_rao(ig, R_t=pl.R_T, L=L100, eps=fd["eps_inv"], kernel="arc", line=line, p_a=pl.P_A)
    rec = t6._record(res, spec, L100, time.time() - t0)
    lr = {"source": "cfd_pinned", "run": pl.FREEZE, "res": pl.RES, "n": pl.N_LINE, "match": match,
          **pl._setup(line)["initial_line_detail"]}
    rec.update({"line_source": res.get("line_source"), "initial_line": lr, "eps_inv": fd["eps_inv"],
                "delta_E": fd["delta_E"], "k_delta": fd["k_delta"],
                "note": "δ* 補正の非粘性壁 (§5.1 #7 (7b-4)): ピン版 L100 と同じ初期線・L で ε = ε_inv"})
    pl._save(PARTS, rec, res)
    print(t6._line(rec), flush=True)
    setup = {**pl._setup(line), "eps": fd["eps_inv"], "eps_phys": EPS_PHYS,
             "note": "ε は非粘性壁の ε_inv = (√50 − kδ·C(x_E))²。物理壁は offset_<variant>.npz を足す (§5.1 #7 (7b-4))"}
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #7 (7b-4)",
           "decision_record": "notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md",
           "setup": setup, "setup_by_tag": {tag: setup}, **pl.git_info(),
           "frozen": {k: fd[k] for k in ("k_delta", "x_E", "C_E", "delta_E", "eps_inv", "sha256_16")},
           "rows": [rec], "design_total_s": time.time() - t0}
    pl._dump(OUT, out)
    print(f"design 完了 ({time.time() - t0:.0f} s) → {PARTS / tag}_wall.npz", flush=True)


def cmd_offset() -> None:
    from forge_design.geometry.wall import DisplacementOffset
    fd = frozen_delta()
    m = (fd["x"] > 0.0)
    xs = np.r_[0.0, fd["x"][m]]
    xs[-1] = L100   # 終点を壁の x_E (生成器の L) に揃える (contur の終点との差 1.5e-9)
    ds = fd["k_delta"] * np.r_[np.interp(0.0, fd["x"], fd["C"]), fd["C"][m]]
    recs = {}
    for v, (a, b) in VARIANTS.items():
        off = DisplacementOffset.fit(xs, ds, ramp=(a, b), meta={
            "variant": v, "k_delta": fd["k_delta"], "contur": CONTUR.name, "calibration": CAL.name,
            "sha256_16": fd["sha256_16"], "eps_phys": EPS_PHYS, "eps_inv": fd["eps_inv"],
            "plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #7 (7b-4)"})
        PARTS.mkdir(exist_ok=True)
        off.save(PARTS / f"offset_{v}.npz")
        r = off.record()
        r["delta_at_xE_minus_frozen"] = float(off(np.array([fd["x_E"]]))[0] - fd["delta_E"])
        recs[v] = r
        (PARTS / f"offset_{v}.json").write_text(json.dumps(r, ensure_ascii=False, indent=1))
        print(v, json.dumps({k: r[k] for k in ("fit_rms", "fit_max", "fit_end_resid", "d0", "d1_0", "d2_0", "d1_max", "x_d1_max", "d2_absmax", "x_d2_absmax", "delta_at_xE_minus_frozen")}), flush=True)


def cmd_check() -> None:
    from forge_design.evaluate.runner import _build_wall_rao
    from forge_design.probdef import load_problem
    out = {}
    for k, p in PROBLEMS.items():
        if not p.exists():
            out[k] = {"skipped": f"{p.name} が無い"}
            continue
        w, rec = _build_wall_rao(load_problem(p))
        xs = np.linspace(0.0, w.x_e, 20001)
        r1, r2 = w.r(xs, 1), w.r(xs, 2)
        kap = r2 / (1 + r1 * r1) ** 1.5
        out[k] = {"problem": p.name, "validate": w.validate(), "throat_jump": w.throat_jump(),
                  "r_exit": float(w.r(np.array([w.x_e]))[0]), "eps_exit": float(w.r(np.array([w.x_e]))[0]) ** 2,
                  "rp_min": float(r1.min()), "kappa_min": float(kap.min()), "kappa_max": float(kap.max()),
                  "checks": rec.get("checks"), "offset": rec.get("offset")}
        print(k, json.dumps({kk: out[k][kk] for kk in ("validate", "r_exit", "eps_exit", "rp_min", "kappa_min")}), flush=True)
    (HERE / "rao_corr_check.json").write_text(json.dumps(t6._jsonable(out), ensure_ascii=False, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("design", "offset", "check"))
    a = ap.parse_args()
    {"design": cmd_design, "offset": cmd_offset, "check": cmd_check}[a.cmd]()


if __name__ == "__main__":
    main()
