#!/usr/bin/env python3
"""§6 V1a (plan tooling-design-orchestrator §5.1 #6d、事前登録 2026-10-17): 実行役が case/71 の 1 pass のために作る起動の入力を、V6 の run
(run_0014_chain_v2) のものと比べる。実行役の側は recipe の NS の段と同じ CLI を --prepare-only で回した一時の run (forge は起動しない)。

使い方: python compare.py <実行役の run (prepare-only)> <run_0014_chain_v2> <run_0008_ns_b0 の写し> --out result.json
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import yaml

DESIGN = Path(__file__).resolve().parents[4] / "design"
sys.path.insert(0, str(DESIGN))
from forge_design.campaign import ic as IC  # noqa: E402
from forge_design.campaign import plan as EP  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402


def sha(p):
    p = Path(p)
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def flat(d, pre=()):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, pre + (str(k),)))
    elif isinstance(d, list) and d and all(isinstance(x, dict) for x in d):
        for i, v in enumerate(d):
            out.update(flat(v, pre + (str(i),)))
    else:
        out[pre] = d
    return out


def diff(a, b):
    fa, fb = flat(a), flat(b)
    return sorted(".".join(k) for k in set(fa) | set(fb) if fa.get(k, "<無>") != fb.get(k, "<無>"))


ap = argparse.ArgumentParser()
ap.add_argument("exec_run")
ap.add_argument("v6_run")
ap.add_argument("ic_src")
ap.add_argument("--problem", required=True)
ap.add_argument("--card", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
E, V = Path(a.exec_run), Path(a.v6_run)
res = {"exec_run": str(E), "v6_run": str(V), "items": {}}

# 計画の確定 (準備の後・起動の前)
ic = IC.resolve_ic(Path(a.ic_src) / "res_90000.h5", a.problem, for_ns=True)
appr = E.parent / "dummy_approval.json"
appr.write_text("{}")
pin = EP.fixed_inputs(card=a.card, problem=a.problem, approval=appr, mesh_msh=E / "nozzle.msh", ic=ic,
                      euler_ref=Path(json.loads((E / "prepare_info.json").read_text())["euler_ref"]["run"]) / "res_12000.h5")
if (E / EP.PLAN_FILE).exists():
    (E / EP.PLAN_FILE).unlink()
doc = EP.confirm(E, pin, stages="full")

# (1) 全段の設定
v6m = json.loads((V / "stage_manifest.json").read_text())
got = [(s["tag"], s["fnv"]) for s in doc["confirmed"]["stages"]]
want = [(s["tag"], s["cfg_fnv"]) for s in v6m["stages"]]
res["items"]["1_stage_configs"] = {"ok": got == want, "exec": got, "v6": want}
# (2)(3) バイト一致
for k, f in (("2_bcond", "bcondConfig.yaml"), ("3_prestage", "solverConfig_prestage.yaml")):
    res["items"][k] = {"ok": sha(E / f) == sha(V / f), "exec": sha(E / f), "v6": sha(V / f)}
# (4) 停止の基準
ce, cv = yaml.safe_load((E / "settle_criteria.yaml").read_text()), yaml.safe_load((V / "settle_criteria.yaml").read_text())
d4 = diff(ce, cv)
allowed4 = {"source.generator_sha256"}
res["items"]["4_settle_criteria"] = {"ok": not [x for x in d4 if x not in allowed4], "diffs": d4, "allowed": sorted(allowed4)}
# (5) 節点の表
nodes = lambda d: [ln for ln in (d / "settle_nodes.txt").read_text().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]  # noqa: E731
res["items"]["5_settle_nodes"] = {"ok": nodes(E) == nodes(V), "n": len(nodes(E))}
# (6) 格子
he, hv = RA.h5_mesh_hashes(E / "nozzle.h5"), RA.h5_mesh_hashes(V / "nozzle.h5")
res["items"]["6_mesh"] = {"ok": he == hv and sha(E / "nozzle.msh") == sha(V / "nozzle.msh"), "h5_exec": he, "h5_v6": hv,
                          "msh_exec": sha(E / "nozzle.msh"), "msh_v6": sha(V / "nozzle.msh")}
# (7) 準備の記録の主な欄
pe, pv = json.loads((E / "prepare_info.json").read_text()), json.loads((V / "prepare_info.json").read_text())
keys = ("initializer", "Md", "x_E", "scale_m", "physical_wall", "main_line", "stop_rule")
d7 = {}
for k in keys:
    dd = diff({k: pe.get(k)}, {k: pv.get(k)})
    dd = [x for x in dd if not any(t in x for t in ("path", "time", "created", "file", "run", "_dir", "sha256"))]
    if dd:
        d7[k] = dd
lim_e, lim_v = (yaml.safe_load((E / "solverConfig.yaml").read_text()).get("space") or {}), (yaml.safe_load((V / "solverConfig.yaml").read_text()).get("space") or {})
lim = {k: (lim_e.get(k), lim_v.get(k)) for k in ("limiterRoRef", "limiterPRef", "limiterARef", "limiterRefLength") if lim_e.get(k) != lim_v.get(k)}
res["items"]["7_prepare_info"] = {"ok": not d7 and not lim, "diffs": d7, "limiter_ref_diffs": lim,
                                  "note": "パス・時刻・ファイル名・run 名・sha256 の欄は正規化で除く"}
# (8) 初期場
res["items"]["8_ic"] = {"ok": (pe.get("ic_snapshot") or {}).get("sha256") == ic["snapshot"]["sha256"] and pe["ic_snapshot"]["step"] == 90000,
                        "exec_ic_snapshot": pe.get("ic_snapshot"), "v6_ic_from": pv.get("ic_from")}
# (9) バイナリ・変換器
res["items"]["9_binary"] = {"ok": bool((doc["inputs"].get("forge_bin") or {}).get("sha256")),
                            "exec_forge_bin": doc["inputs"].get("forge_bin"), "exec_converter": doc["inputs"].get("converter"),
                            "v6_provenance": [ln for ln in (V / "RUN_PROVENANCE.txt").read_text().splitlines() if ln.startswith(("forge_bin", "forge_sha256"))],
                            "note": "値の一致は求めない (手元と AWS B の FP64 のビルドは別; 事前登録の (9))"}
res["verdict"] = "PASS" if all(v["ok"] for v in res["items"].values()) else "FAIL"
Path(a.out).write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))
for k, v in res["items"].items():
    print(f"{'ok  ' if v['ok'] else 'FAIL'} {k}  {json.dumps({kk: vv for kk, vv in v.items() if kk in ('diffs', 'limiter_ref_diffs')}, ensure_ascii=False)[:300]}")
print(f"VERDICT: {res['verdict']}")
