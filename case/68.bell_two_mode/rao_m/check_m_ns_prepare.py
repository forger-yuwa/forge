#!/usr/bin/env python3
r"""③ベル (M) の NS 評価用の問題 YAML (`../problem_bell_tp_m_ns_l{70,80,90,100}.yaml`) が runner の prepare を**変換器の手前まで**
通ることを確かめる (forge・変換器は回さない)。plan tooling-nozzle-bell-two-mode §5.1 #8・§6 V3。

方法は `design/tests/run_bell_cv_gate_tests.py` の (M) の prepare と同じ: `runner._run_converter` を偽物に替え、品質検査用の cell 変換
(`nozzle_qc.h5`) だけは生成した格子の四角形と座標を h5 に書いて返し、本変換 (`nozzle.h5`) の呼び出しで止める。したがって
本物の `check_mesh_quality.py --ar-max <mesh.ar_max>` (MESH_QUALITY.txt の VERDICT) と AR 緩和の条件の検査
(`runner.mesh_ar_relax_check`: AR > 1000 の全セルでスキュー ≤ 0.1、入口・出口面が平面 → MESH_AR_RELAX.json) が通る。
加えて WALL_DESIGN.json (照合・接合・moc_reference)、本変換に渡す solverConfig (node・SST 低 Re・viscMethod 2 + transport・soft 段)・
bcondConfig (断熱 no-slip 壁・outflow 出口) を読み、格子の第一セル比と壁法線の層の端の断面を記録する。

出力: `rao_pin_ns_prepare.json` (1 YAML 1 行)。run dir は一時ディレクトリ (`--tmp`、既定は mkdtemp) に作って残さない。
実行: `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/check_m_ns_prepare.py [L70 L80 ...]`
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner  # noqa: E402
from forge_design.meshing.mesh2d import generate_axisym_mesh  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

LABS = ("L70", "L80", "L90", "L100")
OUT = HERE / "rao_pin_ns_prepare.json"


class _StopAtConverter(Exception):
    pass


def write_qc_h5(path, coords, quads):
    """primal の四角形 (VIZMESH/CONNE、XDMF の四角形 = 5) と座標だけの h5 (check_mesh_quality.py と mesh_ar_relax_check が読む分;
    run_bell_cv_gate_tests.py の同名関数と同じ)。"""
    with h5py.File(path, "w") as f:
        f["MESH/COORD"] = coords.astype(np.float64).ravel()
        f["VIZMESH/CONNE"] = np.column_stack([np.full(len(quads), 5), quads]).astype(np.int32).ravel()
        f["VIZMESH"].attrs["nVizCells"] = len(quads)


def check_one(lab: str, tmp: Path) -> dict:
    t0 = time.time()
    prob = CASE / f"problem_bell_tp_m_ns_{lab.lower()}.yaml"
    p = load_problem(prob)
    wall, rec = runner.build_wall(p)
    mp = runner.mesh_params(p, float(p.spec["r_throat"]))
    cm, qm, bm = generate_axisym_mesh(wall, mp)
    seen = []

    def _fake(conv, msh, out, run_dir):
        seen.append((out, (Path(run_dir) / "solverConfig.yaml").read_text()))
        if out == "nozzle_qc.h5":
            write_qc_h5(Path(run_dir) / out, cm, qm)
            return
        raise _StopAtConverter(out)
    orig = runner._run_converter
    runner._run_converter = _fake
    rd = tmp / f"prep_{lab}"
    stopped = None
    try:
        runner.prepare(prob, rd)
    except _StopAtConverter as e:
        stopped = str(e)
    finally:
        runner._run_converter = orig
    mq = (rd / "MESH_QUALITY.txt").read_text()
    vq = [ln for ln in mq.splitlines() if ln.startswith("VERDICT")]
    rx = json.loads((rd / "MESH_AR_RELAX.json").read_text())
    wd = json.loads((rd / "WALL_DESIGN.json").read_text())
    sc_main = yaml.safe_load(seen[-1][1]) if seen else {}
    sc_qc = yaml.safe_load(seen[0][1]) if seen else {}
    bc = yaml.safe_load((rd / "bcondConfig.yaml").read_text())
    Rg = cm[:, 1].reshape(mp.ni, mp.nj) / 0.14
    Xg = cm[:, 0].reshape(mp.ni, mp.nj) / 0.14
    first = np.hypot(Xg[:, -1] - Xg[:, -2], Rg[:, -1] - Rg[:, -2]) / Rg[:, -1]
    out = {"problem": str(prob.relative_to(ROOT)), "stopped_at": stopped,
           "converter_calls": [s[0] for s in seen],
           "reached_converter": stopped == "nozzle.h5" and [s[0] for s in seen] == ["nozzle_qc.h5", "nozzle.h5"],
           "mesh_quality": {"verdict": vq[-1] if vq else None, "lines": [ln for ln in mq.splitlines() if ln.strip()][-6:]},
           "ar_relax": rx,
           "mesh": {"ni": mp.ni, "nj": mp.nj, "n_cells": int(len(qm)), "x_e": float(wall.x_e),
                    "first_cell_over_r_min": float(first.min()), "first_cell_over_r_max": float(first.max()),
                    "first_cell_over_r_throat": float(first[int(np.argmin(np.abs(Xg[:, -1])))]),
                    "first_cell_over_r_exit": float(first[-1])},
           "wall_design": {k: wd.get(k) for k in ("bell_mode", "design_status", "L", "theta_B_deg", "theta_e_deg", "theta_max_deg",
                                                  "wall_file", "wall_file_sha256", "checks", "joints")},
           "campaign_setup": (wd.get("campaign") or {}).get("setup"),
           "campaign_setup_source": (wd.get("campaign") or {}).get("setup_source"),
           "moc_reference": {k: wd["moc_reference"].get(k) for k in ("F_cs_N", "eta_F_cs", "F_wall_N", "eta_F_wall", "F_1D_N",
                                                                     "mdot_ratio_1d", "accounts_rel_diff")},
           "solverConfig_main": {"discretization": sc_main.get("mesh", {}).get("discretization"),
                                 "viscMethod": sc_main.get("physProp", {}).get("viscMethod"),
                                 "transport": "transport" in sc_main.get("physProp", {}),
                                 "turbulence": sc_main.get("turbulence"), "space": sc_main.get("space"),
                                 "nStepOuter": sc_main.get("time", {}).get("last", {}).get("nStepOuter")},
           "solverConfig_qc_discretization": sc_qc.get("mesh", {}).get("discretization"),
           "bcond": {k: (v or {}).get("kind") for k, v in bc.items()},
           "time_s": time.time() - t0}
    sl = out["solverConfig_main"]
    turb = sl.get("turbulence") or {}
    out["ok"] = bool(out["reached_converter"] and out["mesh_quality"]["verdict"] and out["mesh_quality"]["verdict"].startswith("VERDICT: PASS")
                     and rx["ok"] and rx["required"] and wd.get("design_status") == "OK"
                     and abs(wd["checks"]["L_rel_diff"]) <= 1e-9 and abs(wd["checks"]["eps_rel_diff"]) <= 1e-9
                     and wd["checks"]["gas"]["verified"] and sl["discretization"] == "node" and sl["viscMethod"] == 2
                     and sl["transport"] and str(turb.get("model")).lower() == "sst" and int(turb.get("wallTreatmentSST", -1)) == 0
                     and out["bcond"].get("wall") == "wall" and out["bcond"].get("outlet") == "outflow"
                     and out["solverConfig_qc_discretization"] == "cell")
    shutil.rmtree(rd, ignore_errors=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("labs", nargs="*", default=list(LABS))
    ap.add_argument("--tmp", default=None)
    a = ap.parse_args()
    tmp = Path(a.tmp) if a.tmp else Path(tempfile.mkdtemp(prefix="m_ns_prep_"))
    tmp.mkdir(parents=True, exist_ok=True)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    rows = dict(prev.get("rows") or {})
    for lab in a.labs:
        try:
            r = check_one(lab, tmp)
        except Exception as e:  # noqa: BLE001
            import traceback
            r = {"ok": False, "error": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
        rows[lab] = r
        if r.get("error"):
            print(f"{lab}: ERROR {r['error']}", flush=True)
            continue
        rx = r["ar_relax"]
        print(f"{lab}: {'OK' if r['ok'] else 'NG'} 変換器の手前まで {r['reached_converter']}・{r['mesh_quality']['verdict']}・"
              f"AR 緩和 {'PASS' if rx['ok'] else 'FAIL'} (AR > 1000 の {rx['n_high_ar']} セルでスキュー最大 {rx['max_skew_high_ar']}・"
              f"違反 {rx['n_violations']}・入口 {rx['inlet_planar']}・出口 {rx['outlet_planar']}・AR 最大 {rx['max_ar']:.1f}・"
              f"skew 最大 {rx['max_skew']:.3f})・η_cs {r['moc_reference']['eta_F_cs']:.6f} ({r['time_s']:.0f} s)", flush=True)
    import subprocess
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    OUT.write_text(json.dumps({"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #8・§6 V3 ((M) の NS 評価の準備)",
                               "note": "壁は δ* 補正なしの非粘性壁のまま (#7 未決)。prepare は変換器の手前まで (変換器は偽物で止めた)",
                               "commit": head, "rows": rows}, ensure_ascii=False, indent=1, default=float))
    print(OUT.name, flush=True)


if __name__ == "__main__":
    main()
