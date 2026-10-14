#!/usr/bin/env python3
r"""③ベル (M): δ* 補正壁を生産の手順 (`forge_design.feedback.bell_deltastar_correction.correct_rao`、縁 = MOC の壁 + 凍結源 Euler、
V-ext) で L70〜L100 に作る — plan tooling-nozzle-bell-two-mode §5.1 #7c と V3 の補正列 (codex 2026-10-13 option-design 案 2)。

  design [--tags L100 L70 ...] [--workers 2] : 補正しない (M) の NS の問題 (`problem_bell_tp_m_ns_l{k}.yaml`) から補正壁を作り、
        `rao_corrprod_parts/` に壁・記録・offset、`rao_corrprod.json` に setup_by_tag、`../problem_bell_tp_m_ns_l{k}_corr.yaml` を書く (1 本 約 25 分)。
  check : (7c) の合格 — L100 の offset が旧手順 (Euler の縁、`rao_corr_parts/offset_vext.npz`) と δ の差 ≤ 0.5 % (x ≥ 0.41)、
        ε_inv の差 ≤ 1e-3。各 tag で補正壁の問題を runner で組み、出口の面積比が spec (相対 1e-9)、スロートの r・r′・r″ が補正しない壁と同じ、
        validate が空。出力 rao_corrprod_check.json。

手元 (並列 2 本まで、nice): nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_corrected_walls_prod.py design
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))

FREEZE = CASE / "run_0015_tp_m_euler_l100"
PARTS = HERE / "rao_corrprod_parts"
TAGS = ("L100", "L70", "L80", "L90")
OLD_OFFSET = HERE / "rao_corr_parts" / "offset_vext.npz"
OLD_RECORD = HERE / "rao_corr.json"


def src_problem(k: str) -> Path:
    return CASE / f"problem_bell_tp_m_ns_{k.lower()}.yaml"


def out_problem(k: str) -> Path:
    return CASE / f"problem_bell_tp_m_ns_{k.lower()}_corr.yaml"


def tag_of(k: str) -> str:
    return f"{k}_corr_vext_n321"


def job(k: str) -> dict:
    from forge_design.feedback.bell_deltastar_correction import correct_rao
    kd = float(json.loads((HERE / "kdelta_l100_calibration.json").read_text())["k_delta"])
    return correct_rao(src_problem(k), FREEZE, k_delta=kd, parts_dir=PARTS, tag=tag_of(k))


def write_problem(k: str, r: dict) -> None:
    src = src_problem(k).read_text()
    body = src[src.index("name:"):]
    wall_line = [ln for ln in body.splitlines() if ln.strip().startswith("rao_wall_file:")][0]
    rel = f"rao_m/{PARTS.name}/{tag_of(k)}"
    body = body.replace(wall_line, f"  rao_wall_file: {rel}_wall.npz   # 補正後の非粘性壁 (ε_inv {r['eps_inv']:.6f})\n"
                                   f"  rao_offset_file: {rel}_offset.npz   # 物理壁 = 非粘性壁 + δ(x) (V-ext)")
    name_line = [ln for ln in body.splitlines() if ln.startswith("name:")][0]
    body = body.replace(name_line, name_line + "_corr")
    hdr = (f"# forge_design 問題定義 — ③ベル (M) {k} の **δ* 補正壁** (生産の手順、plan tooling-nozzle-bell-two-mode §4.4・§5.1 #7c、V3 の補正列)。\n"
           f"# 元の問題 {src_problem(k).name} (補正しない) から forge_design.feedback.bell_deltastar_correction.correct_rao で作った:\n"
           f"# kδ {r['k_delta']}・縁 = MOC の壁 + 凍結源 Euler {FREEZE.name}・V-ext (x 0〜0.41 の立ち上げ)・δ_e {r['delta_E']:.6f} r_t・"
           f"ε_inv {r['eps_inv']:.6f}。台本 rao_m/run_corrected_walls_prod.py design。\n"
           f"# NS の設定・格子族・段は元の問題と同じ。\n")
    out_problem(k).write_text(hdr + body)


def cmd_design(tags, workers: int) -> None:
    from forge_design.feedback.bell_deltastar_correction import write_campaign
    if workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    done = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(job, k): k for k in tags}
        for fu in as_completed(futs):
            k = futs[fu]
            try:
                r = fu.result()
            except Exception as e:  # noqa: BLE001
                print(f"{k}: ERROR {type(e).__name__}: {e}\n{traceback.format_exc()}", flush=True)
                continue
            done.append(r)
            write_campaign(PARTS, [r], extra={"plan": "plans/active/tooling-nozzle-bell-two-mode.md §4.4・§5.1 #7c", "commit": head})
            write_problem(k, r)
            w = r["wall_record"]
            print(f"{k}: ok ε_inv {r['eps_inv']:.6f} δ_e {r['delta_E']:.6f} θ_B {w['theta_B_deg']:.4f}° θ_E {w['theta_E_deg']:.4f}° "
                  f"det_D {w['det_D']:+.4f} 推力差 {w['thrust_rel_diff']:+.2e} ({r['time_s']:.0f} s)", flush=True)
    print(f"design 完了 {len(done)}/{len(tags)} ({time.time() - t0:.0f} s)", flush=True)


def cmd_check() -> None:
    from forge_design.evaluate.runner import _build_wall_rao
    from forge_design.geometry.wall import DisplacementOffset
    from forge_design.probdef import load_problem
    out = {}
    if (PARTS / f"{tag_of('L100')}_offset.npz").exists():
        new, old = DisplacementOffset.load(PARTS / f"{tag_of('L100')}_offset.npz"), DisplacementOffset.load(OLD_OFFSET)
        xs = np.linspace(0.41, new.x_end, 4001)
        rel = np.abs(new(xs) / old(xs) - 1.0)
        eo = float(json.loads(OLD_RECORD.read_text())["frozen"]["eps_inv"])
        out["L100_vs_old"] = {"delta_rel_max_x_ge_0p41": float(rel.max()), "x_at_max": float(xs[int(np.argmax(rel))]),
                              "eps_inv_new": new.meta["eps_inv"], "eps_inv_old": eo, "eps_inv_diff": new.meta["eps_inv"] - eo,
                              "pass": bool(rel.max() <= 0.005 and abs(new.meta["eps_inv"] - eo) <= 1e-3)}
    for k in TAGS:
        if not out_problem(k).exists():
            continue
        p = load_problem(out_problem(k))
        w, rec = _build_wall_rao(p)
        w0, _ = _build_wall_rao(load_problem(src_problem(k)))
        z = np.array([0.0])
        tj = {d: float(w.r(z, d)[0] - w0.r(z, d)[0]) for d in (0, 1, 2)}
        eps_e = float(w.r(np.array([w.x_e]))[0]) ** 2
        out[k] = {"eps_exit": eps_e, "eps_rel": eps_e / float(p.spec["eps"]) - 1.0, "throat_vs_uncorrected": tj, "validate": w.validate(),
                  "offset": {kk: rec["offset"][kk] for kk in ("delta_end", "d1_max", "d2_absmax", "fit_rms", "ramp")},
                  "pass": bool(abs(eps_e / float(p.spec["eps"]) - 1.0) <= 1e-9 and max(abs(v) for v in tj.values()) <= 1e-12
                               and not w.validate())}
        print(k, json.dumps(out[k], ensure_ascii=False), flush=True)
    if "L100_vs_old" in out:
        print("L100 vs 旧手順", json.dumps(out["L100_vs_old"], ensure_ascii=False))
    (HERE / "rao_corrprod_check.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    d = sp.add_parser("design"); d.add_argument("--tags", nargs="*", default=list(TAGS)); d.add_argument("--workers", type=int, default=2)
    sp.add_parser("check")
    a = ap.parse_args()
    if a.cmd == "design":
        cmd_design(a.tags, a.workers)
    else:
        cmd_check()


if __name__ == "__main__":
    main()
