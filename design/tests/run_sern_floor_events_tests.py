#!/usr/bin/env python3
"""毎更新の EOS 床事象ゲートの単体テスト (design/.venv-opt の Python で実行)。

plans/active/tooling-sern-te-wake-grid.md §4「床の判定」(2026-10-08 のやり直し)・§5.1 #2 (codex plan レビュー M1・M2):
`sern_gates.floor_event_gate` / `evaluate_gates(require_floor_events=...)` が、生産系列で
「途中で 1 件・最終場は正常」「記録の欠落」「restart 後の区間の欠落」を拒否し、全 0 件だけを通すこと。
既存の呼び出し (必須にしない既定) の合否は変えないこと。runner の config・metrics の要約・driver の台帳までの接続。
判定の正本は solver_density_cuda/tools/check_floor_events.py (不正入力の試験は tools/test_gate_bad_input.py)。
"""
import csv
import json
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
import test_gate_bad_input as TGB  # noqa: E402  (合成記録の生成器; sern_gates が tools を sys.path に足す)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


tmp = Path(tempfile.mkdtemp(prefix="sern_floor_events_"))


def make_run(d: Path, steps=6000, rc=0):
    """健全な最終場 (floor_gate・field_health は合格) とプラトーの残差を持つ合成 run。"""
    d.mkdir(parents=True, exist_ok=True)
    n = 50
    with h5py.File(d / f"res_{steps}.h5", "w") as f:
        for k, val in (("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3.0e5), ("P", 1.0e5), ("T", 300.0)):
            f.create_dataset(f"VALUE/{k}", data=np.full(n, val))
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roe"])
        for s in range(0, steps + 1, 100):
            v = 10 ** (-1 - 1.5 * min(s / steps * 3, 1))
            w.writerow([s, -1, "outer_end", v, v * 10, v * 100])
    (d / "run_case_stdout.log").write_text(f"[run_case] forge exit={rc}\n")
    return d


steady_hist = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0} for i in range(8)]

# --- (1) ゲート ---------------------------------------------------------------------------------------
run_ok = make_run(tmp / "fe_ok"); TGB.write_floor_csv(str(run_ok), 6000)
g = G.evaluate_gates(run_ok, steady_hist, 0, require_floor_events=True)
check("必須: 全更新 0 件 → PASS (区間 (3000, 6000]、3000 更新)", g["verdict"] == "PASS" and g["floor_events"]["verdict"] == "PASS"
      and g["floor_events"]["window"] == [3000, 6000] and g["floor_events"]["n_updates_checked"] == 3000)
run_mid = make_run(tmp / "fe_mid"); TGB.write_floor_csv(str(run_mid), 6000, events={4321: (1, 0, 0)})
g = G.evaluate_gates(run_mid, steady_hist, 0, require_floor_events=True)
check("必須: 途中で 1 件・最終場は正常 (floor_gate・field_health は ok) → FAIL/FLOOR_EVENT",
      g["verdict"] == "FAIL" and g["fail_class"] == "FLOOR_EVENT" and g["floors"]["ok"] and g["field"]["ok"]
      and g["floor_events"]["events"]["T"]["first_q_index"] == 4321)
run_none = make_run(tmp / "fe_none")
g = G.evaluate_gates(run_none, steady_hist, 0, require_floor_events=True)
check("必須: 記録の欠落 (floor_events.csv 無し) → FAIL/FLOOR_UNVERIFIABLE (0 件として扱わない)",
      g["verdict"] == "FAIL" and g["fail_class"] == "FLOOR_UNVERIFIABLE")
run_drop = make_run(tmp / "fe_drop"); TGB.write_floor_csv(str(run_drop), 6000, drop=(5000,))
check("必須: 区間の step の欠け → FAIL/FLOOR_UNVERIFIABLE",
      G.evaluate_gates(run_drop, steady_hist, 0, require_floor_events=True)["fail_class"] == "FLOOR_UNVERIFIABLE")
run_rs = make_run(tmp / "fe_restart"); TGB.write_floor_csv(str(run_rs), 6000, sessions=[20000])
g = G.evaluate_gates(run_rs, steady_hist, 0, require_floor_events=True, floor_window_steps=10000)
check("必須: restart 後の区間の欠落 (区間 10000 step > restart 後の 6000 step) → FAIL/FLOOR_UNVERIFIABLE",
      g["fail_class"] == "FLOOR_UNVERIFIABLE" and "restart" in " ".join(g["reasons"]))
g = G.evaluate_gates(run_rs, steady_hist, 0, require_floor_events=True, floor_window_steps=5000)
check("必須: restart 後の session に区間が収まれば判定 → PASS", g["verdict"] == "PASS")
run_end = make_run(tmp / "fe_noend"); TGB.write_floor_csv(str(run_end), 6000, end=False)
check("必須: session_end 無し (異常終了) → FAIL/FLOOR_UNVERIFIABLE",
      G.evaluate_gates(run_end, steady_hist, 0, require_floor_events=True)["fail_class"] == "FLOOR_UNVERIFIABLE")
run_old = make_run(tmp / "fe_stale", steps=6000); TGB.write_floor_csv(str(run_old), 4000)
check("必須: 記録の最終 step (4000) と res_*.h5 (6000) が違う (古い記録) → FAIL/FLOOR_UNVERIFIABLE",
      G.evaluate_gates(run_old, steady_hist, 0, require_floor_events=True)["fail_class"] == "FLOOR_UNVERIFIABLE")
g = G.evaluate_gates(run_mid, steady_hist, 0)
check("既定 (必須でない): 既存の呼び出しの合否は変えない (PASS)、結果は notes に載る",
      g["verdict"] == "PASS" and not g["floor_events"]["required"] and g["floor_events"]["verdict"] == "FAIL"
      and g["floor_events"]["notes"] and not g["floor_events"]["reasons"] and g["floor_events"]["fail_class"] is None)
check("既定: 記録が無くても既存の合否は PASS のまま", G.evaluate_gates(run_none, steady_hist, 0)["verdict"] == "PASS")
g = G.evaluate_gates(run_mid, steady_hist, 1, require_floor_events=True)
check("必須: rc≠0 なら fail_class は DIVERGED が先 (理由には床事象も積む)",
      g["fail_class"] == "DIVERGED" and any(r.startswith("床事象") for r in g["reasons"]))

# --- (2) 問題 YAML からの経路 ---------------------------------------------------------------------------
check("floor_events_options: キー無し → 必須でない", G.floor_events_options({}) == {"require_floor_events": False, "floor_window_steps": None})
check("floor_events_options: floor_events 1 + 区間 10000",
      G.floor_events_options({"floor_events": 1, "floor_events_window_steps": 10000}) == {"require_floor_events": True, "floor_window_steps": 10000})
for bad in ({"floor_events": 2}, {"floor_events": 1, "floor_events_window_steps": 0}, {"floor_events": 1, "floor_events_window_steps": "1e4"}):
    try:
        G.floor_events_options(bad); okb = False
    except ValueError:
        okb = True
    check(f"floor_events_options: 不正値 {bad} は拒否", okb)
check("floor_events_config_line: 必須のときだけ output 節", G.floor_events_config_line({"floor_events": 1}) == "output: {floorEvents: 1}\n"
      and G.floor_events_config_line({}) == "")
sm = G.floor_events_summary(G.evaluate_gates(run_mid, steady_hist, 0, require_floor_events=True)["floor_events"])
check("floor_events_summary: 必須・判定・区間・種類別の件数・fail_class", sm["required"] and sm["verdict"] == "FAIL" and sm["n_events"]["T"] == 1
      and sm["window"] == [3000, 6000] and sm["fail_class"] == "FLOOR_EVENT")
check("floor_events_summary(None) は None (床事象を持たない旧 gates)", G.floor_events_summary(None) is None)
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402
prob = load_problem(ROOT.parent / "case" / "46.sern_design" / "problem_moo_frozen_tp_cycle3op.yaml")
R2.select_operating_point(prob, "m6_on")     # prepare と同じく作動点の組成を入れてから config を作る
c0 = R2._solver_config(prob, 100, 50, 0.5, 2851.0)
prob.evaluate["floor_events"] = 1
c1 = R2._solver_config(prob, 100, 50, 0.5, 2851.0)
check("runner: evaluate.floor_events 無し → config に floorEvents が現れない (従来の config)", "floorEvents" not in c0 and "output:" not in c0)
check("runner: evaluate.floor_events 1 → config に output.floorEvents 1 だけが増える (YAML として読める)",
      (yaml.safe_load(c1) or {}).get("output") == {"floorEvents": 1} and c1.replace("output: {floorEvents: 1}\n", "") == c0)
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
check("runner_sern3d: 2D と同じ config 生成を通る", "output: {floorEvents: 1}" in R3._solver_config(prob, 100, 50, 0.5, 2851.0))


# --- (3) driver の台帳 -----------------------------------------------------------------------------------
def fake_runner(fe):
    """R.prepare / run_staged / collect の差し替え。fe[op] = (gate_verdict, fail_class, floor_events_verdict)。"""
    m = types.SimpleNamespace()

    def prepare(prob, rd, op=None, **kw):
        rd = Path(rd); rd.mkdir(parents=True, exist_ok=False)
        info = {"design": {"L_ramp": 10.0, "warnings": []}}
        (rd / "prepare_info.json").write_text(json.dumps(info)); return info

    def run_staged(rd, stages, **kw):
        Path(rd, "res_6000.h5").write_bytes(b""); return 0

    def collect(prob, rd, out_dir=None, rc=None, require_residual_pass=False):
        op = Path(rd).name.split("_", 2)[-1].replace("_retry", "")
        verdict, fc, fev = fe[op]
        st = {k: {"verdict": "STEADY", "detail": ""} for k in ("C_T_with_shear", "C_T", "C_L", "C_M")}
        fe_gate = {"required": True, "verdict": fev, "fail_class": fc if fc and fc.startswith("FLOOR") else None,
                   "window": [10000, 20000], "n_updates_checked": 10000,
                   "events": {"T": {"n_events": 1 if fc == "FLOOR_EVENT" else 0}}, "first_events": [], "reasons": [], "notes": []}
        return {"forge_rc": 0, "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0, "step": 6000,
                "gates": {"verdict": verdict, "fail_class": fc, "reasons": [fc or ""], "objective": "C_T_with_shear",
                          "residual": {"verdict": "NOT CONVERGED (stalled/plateau)"}, "steadiness": {"series": st},
                          "floor_events": fe_gate}}

    m.prepare, m.run_staged, m.collect = prepare, run_staged, collect
    return m


BASE = {"type": "sern_2d", "dv": {k: {"min": 0.0, "max": 1.0} for k in D.DV_ORDER},
        "spec": {"inflow": {"M_in": 2.5, "p_in": 1e5}, "external": {"p_inf": 5e3}, "H_m": 0.1,
                 "operating_points": [{"name": "m6_on", "weight": 0.6, "external": {}}, {"name": "m4_off", "weight": 0.4, "external": {}}]},
        "gas": {"gamma": 1.4, "cp": 1004.5}, "geometry": {"L_ramp_max": 12.0}, "opt": {"cm_min": -7.0},
        "evaluate": {"floor_events": 1}}


def campaign(name):
    d = tmp / name; d.mkdir()
    (d / "problem.yaml").write_text(yaml.safe_dump(json.loads(json.dumps(BASE))))
    return D.SernCampaign(d / "problem.yaml", d / "camp")


X = [0.5] * 5
c = campaign("d_fail"); D.R = fake_runner({"m6_on": ("FAIL", "FLOOR_EVENT", "FAIL"), "m4_off": ("PASS", None, "PASS")})
row = c.evaluate(X, "doe_000")
check("driver: 床事象 (必須) で不合格 → 台帳 FAIL/FLOOR_EVENT、作動点の要約に floor_events、学習に入らない",
      row["status"] == "FAIL" and row["fail_class"] == "FLOOR_EVENT" and row["ops"]["m6_on"]["gate_fail_class"] == "FLOOR_EVENT"
      and row["ops"]["m6_on"]["floor_events"]["verdict"] == "FAIL" and len(c._XF()[0]) == 0)
led = [json.loads(line) for line in (c.dir / "ledger.jsonl").read_text().splitlines() if line.strip()]
check("driver: ledger.jsonl の行にも floor_events の要約と FAIL/FLOOR_EVENT が残る",
      led[-1]["fail_class"] == "FLOOR_EVENT" and led[-1]["ops"]["m6_on"]["floor_events"]["fail_class"] == "FLOOR_EVENT")
c = campaign("d_unv"); D.R = fake_runner({"m6_on": ("FAIL", "FLOOR_UNVERIFIABLE", "INDETERMINATE"), "m4_off": ("PASS", None, "PASS")})
row = c.evaluate(X, "doe_000")
check("driver: 記録が判定不能 → 台帳 FAIL/FLOOR_UNVERIFIABLE", row["status"] == "FAIL" and row["fail_class"] == "FLOOR_UNVERIFIABLE")
c = campaign("d_pass"); D.R = fake_runner({"m6_on": ("PASS", None, "PASS"), "m4_off": ("PASS", None, "PASS")})
row = c.evaluate(X, "doe_000")
check("driver: 全作動点 PASS → 台帳 PASS、作動点の要約に floor_events (必須・PASS・区間)",
      row["status"] == "PASS" and row["ops"]["m6_on"]["floor_events"] == {
          "required": True, "verdict": "PASS", "fail_class": None, "window": [10000, 20000], "n_updates_checked": 10000,
          "n_events": {"T": 0}, "first_events": [], "reasons": []})

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
