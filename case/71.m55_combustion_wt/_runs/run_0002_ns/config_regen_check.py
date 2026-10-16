#!/usr/bin/env python3
r"""case/71 run_0002_ns の実効の設定が問題 YAML と生成器どおりかを照合する (forge は起動しない; plan tooling-nozzle-card-local-mesh §5.1 #3b)。

問題 YAML から生成器 (`runner_axismach` の `_ns_configs`・`_bcond`・`ns_stage_configs`、prepare_ns・run_staged_ns と同じ引数:
`deltastar_loop --init-integral --stages full` は nsteps・cfl_main・implicit_relax を渡さない) で S1・S2・本段の solverConfig と
bcondConfig を作り直し、run に残った記録と比べる:
  - 本段の solverConfig.yaml (停止規則の上限で nStepOuter を書き換えた後) と、作り直した本段の設定に同じ書き換えをした文字列の一致
  - MAIN_STAGE_STARTED の本段・前の段 (S1・S2) の設定の sha256、prepare_info.json の ns_stage_configs の sha256
  - forge_launches.jsonl の起動ごとの cfg_fnv・bcond_fnv (forge が実際に読んだ設定の FNV-1a 64)
  - bcondConfig.yaml・species_meta.yaml の一致、問題 YAML の sha256

使い方 (AWS A): python3 config_regen_check.py --design ~/forge-cardmesh/design --tools ~/forge-cardmesh/solver_density_cuda/tools \
    --problem <problem_m55_final.yaml> --run <run_0002_ns> --out-json <json>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--tools", required=True)
    ap.add_argument("--problem", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(Path(a.design).resolve()))
    sys.path.insert(0, str(Path(a.tools).resolve()))
    from forge_design.evaluate import runner_axismach as RA
    from stage_manifest import fnv1a64

    run = Path(a.run).resolve()
    info = json.loads((run / "prepare_info.json").read_text())
    started = json.loads((run / "MAIN_STAGE_STARTED").read_text())
    launches = [json.loads(l) for l in (run / "forge_launches.jsonl").read_text().splitlines() if l.strip()]
    p = RA.load_problem(a.problem)
    n, out_int, n_src = RA._ns_steps(p, None)
    ml = RA._main_line_spec(p, None, None)
    with tempfile.TemporaryDirectory() as td:
        cfg_ns, cfg_main = RA._ns_configs(p, Path(td), n, out_int, None, None, ml)
        meta = (Path(td) / "species_meta.yaml").read_text()
    bc = RA._bcond_with_species(RA._bcond(p, euler=False, profile_ok=True), RA._tp_species_Y(p))
    sha_prep = RA.ns_stage_sha256(cfg_main, cfg_ns if ml is not None else None)
    pre, _, _ = RA.ns_stage_configs(cfg_main, "full", None, 1000, cfg_ns if ml is not None else None)
    cap = int(info["stop_rule"]["cap"])
    main_run = RA._cfg_set(cfg_main, {"nStepOuter": cap}) if info["stop_rule"].get("use") else cfg_main
    run_cfg = (run / "solverConfig.yaml").read_text()
    run_bc = (run / "bcondConfig.yaml").read_text()
    run_meta = (run / "species_meta.yaml").read_text()
    cfgs = {tag: c for tag, c, _, _ in pre}
    cfgs["main"] = main_run
    exp_fnv = [fnv1a64(cfgs["S1_soft"]), fnv1a64(cfgs["S2_mid"]), fnv1a64(main_run)]
    checks = {
        "problem_sha256": sha(Path(a.problem).read_bytes()),
        "nStepOuter_fallback": {"n": n, "source": n_src, "out_int": out_int, "prepare_info": info.get("nStepOuter")},
        "main_line": ml,
        "stop_rule_cap": cap,
        "main_solverConfig_equal": run_cfg == main_run,
        "main_solverConfig_sha256": {"run": sha(run_cfg.encode()), "regen": sha(main_run.encode()),
                                     "MAIN_STAGE_STARTED": started["main_config"]["sha256"]},
        "prepare_ns_stage_configs": {"prepare_info": info.get("ns_stage_configs"), "regen": sha_prep,
                                     "equal": all(info["ns_stage_configs"].get(k) == sha_prep.get(k) for k in ("S1_soft", "S2_mid", "main"))},
        "previous_stages_sha256": {st["tag"]: {"MAIN_STAGE_STARTED": st["config_sha256"], "regen": sha(cfgs[st["tag"]].encode()),
                                               "equal": st["config_sha256"] == sha(cfgs[st["tag"]].encode())}
                                   for st in started["previous_stages"]},
        "forge_launches_cfg_fnv": {"recorded": [l["cfg_fnv"] for l in launches], "regen_S1_S2_main": exp_fnv,
                                   "equal": [l["cfg_fnv"] for l in launches] == exp_fnv},
        "bcond_equal": run_bc == bc,
        "bcond_fnv": {"recorded": sorted({l["bcond_fnv"] for l in launches}), "regen": fnv1a64(bc)},
        "species_meta_equal": run_meta == meta,
        "problem_values": {"Pt": p.spec.get("Pt"), "Tt": p.spec.get("Tt"), "p_ambient": p.spec.get("p_ambient"),
                           "wall_thermal": p.spec.get("wall_thermal"), "cfl_main": p.evaluate.get("cfl_main"),
                           "stop_rule_keys": {k: p.evaluate.get(k) for k in ("stop_rule", "stop_rule_cap", "stop_rule_hold_steps",
                                                                            "stop_rule_thresholds", "outStepInterval", "ns_steps")}},
        "run_bcond": run_bc,
    }
    if run_cfg != main_run:
        import difflib
        checks["main_solverConfig_diff"] = list(difflib.unified_diff(main_run.splitlines(), run_cfg.splitlines(), "regen", "run", lineterm=""))
    if run_bc != bc:
        import difflib
        checks["bcond_diff"] = list(difflib.unified_diff(bc.splitlines(), run_bc.splitlines(), "regen", "run", lineterm=""))
    ok = (checks["main_solverConfig_equal"] and checks["prepare_ns_stage_configs"]["equal"]
          and all(v["equal"] for v in checks["previous_stages_sha256"].values()) and checks["forge_launches_cfg_fnv"]["equal"]
          and checks["bcond_equal"] and checks["bcond_fnv"]["recorded"] == [checks["bcond_fnv"]["regen"]] and checks["species_meta_equal"]
          and checks["main_solverConfig_sha256"]["run"] == checks["main_solverConfig_sha256"]["MAIN_STAGE_STARTED"])
    checks["verdict"] = "一致 (生成器どおり)" if ok else "不一致あり (欄を見る)"
    Path(a.out_json).write_text(json.dumps(checks, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({k: v for k, v in checks.items() if k != "run_bcond"}, indent=1, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
