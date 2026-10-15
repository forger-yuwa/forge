#!/usr/bin/env python3
r"""case/71 の 3 係数の壁 (§4.3 の末尾、§5.1 #3t) の NS `run_0010_ns_b0` (本段 B0) を作って回す。`run_0008_ns_b0/p3_ns.py` の写しで、PLAN・記録の名前だけ違う。
(plan tooling-nozzle-card-local-mesh §5.1 #3p・§6 V6、ユーザ 2026-10-15「line implicitを使う方法でどうなるか、も同時に確かめてほしい」)。
`run_0007_ns/p2_ns.py` の写しで、違いは: 問題 YAML (`problem_m55_p2_b0.yaml`) の `evaluate.main_line` が本段の設定に B0 (lineImplicit 1・
lineDtDirectional 1・lineDtDirectionalCap 50・implicitThermalJacobian 5・cfl 4・implicitRelax 0.7・リミッタの基準値の固定) を書き、
前の段 (S1・S2) の元の設定は solverConfig_prestage.yaml に置かれる (`runner_axismach.PRESTAGE_CONFIG`)。outStepInterval 2,000 の
書き換えは本段と前の段の両方に当てる (run_0007_ns と同じ出力)。本段の設定に B0 のキーが入っていることを検査する。
使い方は p2_ns.py と同じ (`p3_ns.py prepare --problem .../problem_m55_p2_b0.yaml ... --run-dir .../run_0008_ns_b0` → `p3_ns.py run`)。
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

PLAN = "plans/active/tooling-nozzle-card-local-mesh.md §5.1 #3t (§4.3 の末尾、3 係数)"
RECORD = "P4_NS.json"
ENV = {"double": {"FORGE_TRANSPORT_TABLE": "0", "FORGE_DIAG_FACE_H_DOUBLE": "1", "FORGE_CUDA_BLOCKSIZE": "128"},
       "float": {"FORGE_TRANSPORT_TABLE": None, "FORGE_DIAG_FACE_H_DOUBLE": None, "FORGE_CUDA_BLOCKSIZE": "128"}}   # None = 無いこと
ADD = {"double": {("physProp", "thermoFloat"): 0, ("time", "deltaT", "implicitSolvePrecision"): 1}, "float": {}}
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


def edit_config(text: str, precision: str) -> tuple[str, dict]:
    """生成器の solverConfig (flow 形式) に ADD[precision] を足し、outStepInterval を OUT_INT にする。読み直して差が要求ちょうどであることを検査。"""
    import yaml
    from forge_design.evaluate.runner_axismach import _cfg_set
    fb = flatten(yaml.safe_load(text))
    add = ADD[precision]
    if any(k in fb for k in ADD["double"]):
        raise SystemExit(f"生成器の solverConfig に既に {[k for k in ADD['double'] if k in fb]} がある — 想定と違う (精度の切り替えは既定のはず)")
    for old, new in () if not add else (("thermoHrefTemp: 298.15}", "thermoHrefTemp: 298.15, thermoFloat: 0}"),
                     ("detectNaN: 1}", "detectNaN: 1, implicitSolvePrecision: 1}")):
        if text.count(old) != 1:
            raise SystemExit(f"solverConfig に '{old}' がちょうど 1 回ない ({text.count(old)} 回) — 書き換えない")
        text = text.replace(old, new)
    want = dict(add)
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


def implicit_settings(cfg: dict) -> dict:
    """float のときの記録だけ (検査はしない): 陰解法の設定と implicitSolvePrecision (既定 = キー無し)。"""
    t = cfg["time"]
    return {"timeIntegration": t.get("timeIntegration"), "blockDPLUR": t["deltaT"].get("blockDPLUR"),
            "lowMachPrecond": t["deltaT"].get("lowMachPrecond", 0), "implicitSolvePrecision": t["deltaT"].get("implicitSolvePrecision"),
            "thermoFloat": cfg["physProp"].get("thermoFloat")}


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
    info = DL.run_pass0_integral(a.problem, eref, run, ic_from=str(ic), initializer=None, prepare_only=True, stop_rule=True)
    ini = info.get("initializer") or {}
    if not (isinstance(ini, dict) and "delta_scale" in ini and "delta_offset_m" in ini):
        raise SystemExit(f"準備の initializer に 2 係数が無い ({ini!r}) — 問題 YAML の deltastar_initializer を確かめる")
    if not (run / PRESTAGE_CONFIG).exists():
        raise SystemExit(f"{run / PRESTAGE_CONFIG} が無い (main_line の B0 の run のはず)")
    pre0 = (run / PRESTAGE_CONFIG).read_text()
    (run / "solverConfig_prestage.generated.yaml").write_text(pre0)
    pre, pre_changes = edit_config(pre0, a.precision)
    (run / PRESTAGE_CONFIG).write_text(pre)
    cfg0 = (run / "solverConfig.yaml").read_text()
    (run / "solverConfig.generated.yaml").write_text(cfg0)
    cfg, changes = edit_config(cfg0, a.precision)
    (run / "solverConfig.yaml").write_text(cfg)
    cond = check_implicit_conditions(yaml.safe_load(cfg)) if a.precision == "double" else implicit_settings(yaml.safe_load(cfg))
    dtm = yaml.safe_load(cfg)["time"]["deltaT"]
    b0 = {k: dtm.get(k) for k in ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap", "implicitThermalJacobian", "implicitRelax", "cfl", "cfl_pseudo")}
    if not (b0["lineImplicit"] == 1 and b0["lineDtDirectional"] == 1 and b0["lineDtDirectionalCap"] == 50.0 and b0["implicitThermalJacobian"] == 5):
        raise SystemExit(f"本段の設定に B0 のキーが無い: {b0}")
    if any(k in yaml.safe_load(pre)["time"]["deltaT"] for k in ("lineImplicit", "lineDtDirectionalCap", "implicitThermalJacobian")):
        raise SystemExit("前の段の設定に B0 のキーがある")
    # 停止規則を出力の間隔 2,000 で決め直す (明示・上限 20 万・保持 4 万; Euler の参照は固定のスナップショット)
    pinfo = json.loads((run / "prepare_info.json").read_text())
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_cap": int(a.cap), "stop_rule_hold_steps": int(a.hold)},
                             nsteps_arg=None, n_main=int(a.cap), n_main_source=f"停止規則の上限 ({PLAN})", out_int=OUT_INT,
                             euler_ref=eref)
    if not (sr["use"] and sr["explicit"]):
        raise SystemExit(f"停止規則を明示で使えない: {sr.get('reason')}")
    pinfo["stop_rule_generated"] = pinfo.get("stop_rule")
    pinfo["stop_rule"] = sr
    pinfo["p4_ns"] = {"b0": b0, "prestage_changes": pre_changes, "plan": PLAN, "precision": a.precision, "config_changes": changes, "env_required": ENV[a.precision],
                      "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
                      "note": "solverConfig.generated.yaml が生成器の出力、solverConfig.yaml はそれに config_changes を足したもの"}
    (run / "prepare_info.json").write_text(json.dumps(pinfo, indent=1, ensure_ascii=False, default=str))
    rec = {"plan": PLAN, "prepared": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "git_head": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
           "problem": str(Path(a.problem).resolve()), "problem_sha256": sha256(a.problem),
           "initializer": ini, "ic_from": str(ic),
           "euler_ref": eref, "nozzle_h5_sha256": sha256(run / "nozzle.h5"),
           "solverConfig_sha256": sha256(run / "solverConfig.yaml"),
           "solverConfig_generated_sha256": sha256(run / "solverConfig.generated.yaml"),
           "bcondConfig_sha256": sha256(run / "bcondConfig.yaml"), "config_changes": changes, "implicit_conditions": cond,
           "b0": b0, "prestage_changes": pre_changes, "prestage_sha256": sha256(run / PRESTAGE_CONFIG),
           "stop_rule": {k: sr[k] for k in ("explicit", "cap", "hold_steps", "window", "out_step_interval", "euler_ref", "use")},
           "precision": a.precision, "env_required": ENV[a.precision], "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
           "prepare": {k: info.get(k) for k in ("throat_physical", "dstar_source", "initializer", "nStepOuter", "cfl_main")}}
    (run / RECORD).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(rec, indent=1, ensure_ascii=False, default=str))


def run_main(a) -> int:
    import time
    import yaml
    run = Path(a.run_dir).resolve()
    rec = json.loads((run / RECORD).read_text())
    prec = rec["precision"]
    bad = {k: os.environ.get(k) for k, v in ENV[prec].items() if os.environ.get(k) != v}
    fb = os.environ.get("FORGE_BIN")
    if bad or not fb:
        raise SystemExit(f"環境変数が要求と違う: {bad} / FORGE_BIN={fb!r} (要求 {ENV[prec]} (None = 無いこと) と FORGE_BIN) — 起動しない")
    if sha256(fb) != rec["forge_sha256"]:
        raise SystemExit(f"FORGE_BIN {fb} の sha256 が準備の記録と違う — 起動しない")
    from forge_design.evaluate.runner_axismach import PRESTAGE_CONFIG
    for name, key in (("solverConfig.yaml", "solverConfig_sha256"), ("nozzle.h5", "nozzle_h5_sha256"), (PRESTAGE_CONFIG, "prestage_sha256")):
        if sha256(run / name) != rec[key]:
            raise SystemExit(f"{name} が準備の後に変わっている — 起動しない")
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    want_tf = 0 if prec == "double" else None                     # float は既定 (キー無し = thermoFloat 1)
    if prec == "double":
        check_implicit_conditions(cfg)
    elif cfg["time"]["deltaT"].get("implicitSolvePrecision") is not None:
        raise SystemExit("float なのに implicitSolvePrecision がある — 起動しない")
    if cfg["physProp"].get("thermoFloat") != want_tf or cfg["time"].get("outStepInterval") != OUT_INT:
        raise SystemExit("solverConfig の thermoFloat・outStepInterval が要求と違う — 起動しない")
    from forge_design.evaluate.runner_axismach import run_staged_ns, collect
    from forge_design.feedback import deltastar_loop as DL
    print(f"[p4_ns] {datetime.now(timezone.utc).isoformat(timespec='seconds')} 起動 ({prec}, stages=full, stop_rule=True) pid {os.getpid()}",
          flush=True)
    t0 = time.time()
    rc = run_staged_ns(run, stages="full", stop_rule=True)
    print(f"[p4_ns] forge rc {int(rc)} 停止の状態 {getattr(rc, 'stop_rule', None)} ({time.time() - t0:.0f} s)", flush=True)
    m = collect(rec["problem"], run)
    m["ns_wall_time_s"] = time.time() - t0
    (run / "metrics.json").write_text(json.dumps(m, indent=1))
    # 受け渡しの関門 (上書きなし)・抽出・設計の標準出力 (run_from_table の後段と同じ)
    DL._after_ns(run, rec["euler_ref"], 0.5, m, "extract(new)", None, None, problem=rec["problem"])
    print(f"[p4_ns] {datetime.now(timezone.utc).isoformat(timespec='seconds')} 終了", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--problem", required=True)
    p.add_argument("--ic-from", required=True)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--euler-run", required=True)
    p.add_argument("--euler-res", required=True)
    p.add_argument("--euler-sha256", required=True)
    p.add_argument("--forge-bin", required=True)
    p.add_argument("--forge-sha256", required=True)
    p.add_argument("--cap", type=int, default=200000)
    p.add_argument("--hold", type=int, default=40000)
    p.add_argument("--precision", choices=("float", "double"), default="float",
                   help="float (既定、ユーザ決定 2026-10-15) = 生成器の精度の既定のまま; double = #3d の 4 つの切り替え")
    r = sub.add_parser("run")
    r.add_argument("--run-dir", required=True)
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a)
        return 0
    return run_main(a)


if __name__ == "__main__":
    sys.exit(main())
