#!/usr/bin/env python3
r"""case/71 の δ* の補正の 1 回目の NS `run_0005_ns` を作って回す (plan tooling-nozzle-card-local-mesh §5.1 #3e (5)、
ユーザ決定 2026-10-15「aかな」)。

`deltastar_loop.run_from_table(..., prepare_only=True)` で補正の表の物理壁の NS を準備し (初期場は run_0003_ns_cont の最後の場から
`prepare_ns` の異なる格子への補間)、生成器の solverConfig に次だけを足してから `run_staged_ns(stages="full")` で回す:
  - `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1` (#3d と同じ全部 double; 条件 timeIntegration 11・
    blockDPLUR 1・lowMachPrecond < 2 を検査)
  - `time.outStepInterval: 2000` (停止規則の窓 2 万に 11 点、上限 20 万で全場 100 枚 ≈ 10 GB; 1,000 だと A の容量を超える)
  - 起動の環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1`・`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_BIN` (FP64 のビルド、sha256)
停止規則は明示 (窓 2 万・保持 4 万・上限 20 万)、出力の間隔 2,000 で設定を決め直して prepare_info.json に書く。
本段の後は `deltastar_loop._after_ns` (受け渡しの関門は上書きなし・抽出・設計の標準出力) まで。

使い方 (AWS A):
  cd ~/forge-cardmesh/design
  python3 ../case/71.m55_combustion_wt/_runs/run_0005_ns/p1_ns.py prepare --problem ../case/71.m55_combustion_wt/problem_m55_p1.yaml \
      --table ../case/71.m55_combustion_wt/problem_m55_p1_card_delta_r.csv --ic-from ../case/71.m55_combustion_wt/run_0003_ns_cont \
      --run-dir ../case/71.m55_combustion_wt/run_0005_ns \
      --euler-run ../case/71.m55_combustion_wt/run_0004_euler --euler-res res_12000.h5 --euler-sha256 <64 桁> \
      --forge-bin ~/forge-integ-dflow2-fp64/solver_density_cuda/build/forge --forge-sha256 <64 桁>
  FORGE_BIN=... FORGE_CONVERTER=... FORGE_CUDA_BLOCKSIZE=128 FORGE_TRANSPORT_TABLE=0 FORGE_DIAG_FACE_H_DOUBLE=1 \
      nohup python3 ../case/71.m55_combustion_wt/_runs/run_0005_ns/p1_ns.py run --run-dir ../case/71.m55_combustion_wt/run_0005_ns > <log> 2>&1 < /dev/null &
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]          # リポジトリ直下 (case/71.../_runs/run_0005_ns/ の 4 つ上)
sys.path.insert(0, str(ROOT / "design"))

PLAN = "plans/active/tooling-nozzle-card-local-mesh.md §5.1 #3e"
RECORD = "P1_NS.json"
ENV_REQUIRED = {"FORGE_TRANSPORT_TABLE": "0", "FORGE_DIAG_FACE_H_DOUBLE": "1", "FORGE_CUDA_BLOCKSIZE": "128"}
ADD = {("physProp", "thermoFloat"): 0, ("time", "deltaT", "implicitSolvePrecision"): 1}
OUT_INT = 2000


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def flatten(d, pre=()):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flatten(v, pre + (k,)))
    else:
        out[pre] = d
    return out


def edit_config(text: str) -> tuple[str, dict]:
    """生成器の solverConfig (flow 形式) に ADD を足し、outStepInterval を OUT_INT にする。読み直して差が要求ちょうどであることを検査。"""
    import yaml
    from forge_design.evaluate.runner_axismach import _cfg_set
    fb = flatten(yaml.safe_load(text))
    if any(k in fb for k in ADD):
        raise SystemExit(f"生成器の solverConfig に既に {[k for k in ADD if k in fb]} がある — 想定と違う")
    for old, new in (("thermoHrefTemp: 298.15}", "thermoHrefTemp: 298.15, thermoFloat: 0}"),
                     ("detectNaN: 1}", "detectNaN: 1, implicitSolvePrecision: 1}")):
        if text.count(old) != 1:
            raise SystemExit(f"solverConfig に '{old}' がちょうど 1 回ない ({text.count(old)} 回) — 書き換えない")
        text = text.replace(old, new)
    want = dict(ADD)
    if fb.get(("time", "outStepInterval")) != OUT_INT:
        text = _cfg_set(text, {"outStepInterval": OUT_INT})
        want[("time", "outStepInterval")] = OUT_INT
    fa = flatten(yaml.safe_load(text))
    diff = {k for k in set(fb) | set(fa) if fb.get(k, "<無>") != fa.get(k, "<無>")}
    if diff != set(want) or any(fa[k] != v for k, v in want.items()):
        raise SystemExit(f"solverConfig の差が要求と違う: {sorted(diff)} (要求 {sorted(want)})")
    return text, {".".join(k): {"before": fb.get(k), "after": v} for k, v in want.items()}


def check_implicit_conditions(cfg: dict) -> dict:
    """implicitSolvePrecision 1 の条件 (`solverConfig.cpp:1404-1419`): timeIntegration 11・blockDPLUR 1・lowMachPrecond < 2。"""
    t = cfg["time"]
    got = {"timeIntegration": t.get("timeIntegration"), "blockDPLUR": t["deltaT"].get("blockDPLUR"),
           "lowMachPrecond": t["deltaT"].get("lowMachPrecond", 0), "implicitSolvePrecision": t["deltaT"].get("implicitSolvePrecision")}
    ok = got["timeIntegration"] == 11 and got["blockDPLUR"] == 1 and int(got["lowMachPrecond"]) < 2 and got["implicitSolvePrecision"] == 1
    if not ok:
        raise SystemExit(f"implicitSolvePrecision 1 の条件を満たさない: {got}")
    return {**got, "ok": True}


def prepare(a) -> None:
    import yaml
    from forge_design.evaluate import stop_rule as SR
    from forge_design.evaluate.runner_axismach import PRESTAGE_CONFIG
    from forge_design.feedback import deltastar_loop as DL
    run = Path(a.run_dir).resolve()
    if run.exists():
        raise SystemExit(f"{run} が既にある (既存 run を上書きしない)")
    fb = Path(a.forge_bin).expanduser().resolve()
    if sha256(fb) != a.forge_sha256:
        raise SystemExit(f"{fb} の sha256 が指定と違う")
    eref = {"run": str(Path(a.euler_run).resolve()), "res": a.euler_res, "sha256": a.euler_sha256}
    SR.check_euler_ref(eref)                                   # 無い・sha が違えば ValueError
    ic = Path(a.ic_from).resolve()
    info = DL.run_from_table(a.problem, eref, a.table, run, ic_from=str(ic), prev_run=str(ic), prepare_only=True, stop_rule=True)
    if (run / PRESTAGE_CONFIG).exists():
        raise SystemExit(f"{run / PRESTAGE_CONFIG} がある (main_line の run は想定していない)")
    cfg0 = (run / "solverConfig.yaml").read_text()
    (run / "solverConfig.generated.yaml").write_text(cfg0)
    cfg, changes = edit_config(cfg0)
    (run / "solverConfig.yaml").write_text(cfg)
    cond = check_implicit_conditions(yaml.safe_load(cfg))
    # 停止規則を出力の間隔 2,000 で決め直す (明示・上限 20 万・保持 4 万; Euler の参照は固定のスナップショット)
    pinfo = json.loads((run / "prepare_info.json").read_text())
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_cap": int(a.cap), "stop_rule_hold_steps": int(a.hold)},
                             nsteps_arg=None, n_main=int(a.cap), n_main_source=f"停止規則の上限 ({PLAN})", out_int=OUT_INT,
                             euler_ref=eref)
    if not (sr["use"] and sr["explicit"]):
        raise SystemExit(f"停止規則を明示で使えない: {sr.get('reason')}")
    pinfo["stop_rule_generated"] = pinfo.get("stop_rule")
    pinfo["stop_rule"] = sr
    pinfo["p1_ns"] = {"plan": PLAN, "config_changes": changes, "env_required": ENV_REQUIRED,
                      "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
                      "note": "solverConfig.generated.yaml が生成器の出力、solverConfig.yaml はそれに config_changes を足したもの"}
    (run / "prepare_info.json").write_text(json.dumps(pinfo, indent=1, ensure_ascii=False, default=str))
    rec = {"plan": PLAN, "prepared": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "git_head": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
           "problem": str(Path(a.problem).resolve()), "problem_sha256": sha256(a.problem),
           "table": str(Path(a.table).resolve()), "table_sha256": sha256(a.table), "ic_from": str(ic),
           "euler_ref": eref, "nozzle_h5_sha256": sha256(run / "nozzle.h5"),
           "solverConfig_sha256": sha256(run / "solverConfig.yaml"),
           "solverConfig_generated_sha256": sha256(run / "solverConfig.generated.yaml"),
           "bcondConfig_sha256": sha256(run / "bcondConfig.yaml"), "config_changes": changes, "implicit_conditions": cond,
           "stop_rule": {k: sr[k] for k in ("explicit", "cap", "hold_steps", "window", "out_step_interval", "euler_ref", "use")},
           "env_required": ENV_REQUIRED, "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
           "prepare": {k: info.get(k) for k in ("throat_physical", "dstar_source", "delta_r_table", "nStepOuter", "cfl_main")}}
    (run / RECORD).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(rec, indent=1, ensure_ascii=False, default=str))


def run_main(a) -> int:
    import time
    import yaml
    run = Path(a.run_dir).resolve()
    rec = json.loads((run / RECORD).read_text())
    bad = {k: os.environ.get(k) for k, v in ENV_REQUIRED.items() if os.environ.get(k) != v}
    fb = os.environ.get("FORGE_BIN")
    if bad or not fb:
        raise SystemExit(f"環境変数が要求と違う: {bad} / FORGE_BIN={fb!r} (要求 {ENV_REQUIRED} と FORGE_BIN) — 起動しない")
    if sha256(fb) != rec["forge_sha256"]:
        raise SystemExit(f"FORGE_BIN {fb} の sha256 が準備の記録と違う — 起動しない")
    for name, key in (("solverConfig.yaml", "solverConfig_sha256"), ("nozzle.h5", "nozzle_h5_sha256")):
        if sha256(run / name) != rec[key]:
            raise SystemExit(f"{name} が準備の後に変わっている — 起動しない")
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    check_implicit_conditions(cfg)
    if cfg["physProp"].get("thermoFloat") != 0 or cfg["time"].get("outStepInterval") != OUT_INT:
        raise SystemExit("solverConfig の thermoFloat・outStepInterval が要求と違う — 起動しない")
    from forge_design.evaluate.runner_axismach import run_staged_ns, collect
    from forge_design.feedback import deltastar_loop as DL
    print(f"[p1_ns] {datetime.now(timezone.utc).isoformat(timespec='seconds')} 起動 (stages=full, stop_rule=True) pid {os.getpid()}",
          flush=True)
    t0 = time.time()
    rc = run_staged_ns(run, stages="full", stop_rule=True)
    print(f"[p1_ns] forge rc {int(rc)} 停止の状態 {getattr(rc, 'stop_rule', None)} ({time.time() - t0:.0f} s)", flush=True)
    m = collect(rec["problem"], run)
    m["ns_wall_time_s"] = time.time() - t0
    (run / "metrics.json").write_text(json.dumps(m, indent=1))
    # 受け渡しの関門 (上書きなし)・抽出・設計の標準出力 (run_from_table の後段と同じ)
    DL._after_ns(run, rec["euler_ref"], 0.5, m, "extract(new)", None, None, problem=rec["problem"])
    print(f"[p1_ns] {datetime.now(timezone.utc).isoformat(timespec='seconds')} 終了", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--problem", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--ic-from", required=True)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--euler-run", required=True)
    p.add_argument("--euler-res", required=True)
    p.add_argument("--euler-sha256", required=True)
    p.add_argument("--forge-bin", required=True)
    p.add_argument("--forge-sha256", required=True)
    p.add_argument("--cap", type=int, default=200000)
    p.add_argument("--hold", type=int, default=40000)
    r = sub.add_parser("run")
    r.add_argument("--run-dir", required=True)
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a)
        return 0
    return run_main(a)


if __name__ == "__main__":
    sys.exit(main())
