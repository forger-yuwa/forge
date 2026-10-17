#!/usr/bin/env python3
"""§6 V1a′ (plan tooling-design-orchestrator、事前登録 2026-10-17): 2 つの準備した run の起動の入力を、完全なキーパスで比べる。

    python compare_ab.py noise <B1> <B2> --out noise.json           揺らぎの欄 (同じ経路で 2 回準備して違うキー) を記録する
    python compare_ab.py ab <A> <B1> --noise noise.json --out ab.json  A と B1 を比べる (揺らぎの欄と「別に記すもの」だけを除く)

比べるもの: 全段の起動直前の設定 (`staged_ns_launch_configs` の全文)・bcondConfig.yaml・solverConfig_prestage.yaml・solverConfig.yaml・settle_criteria.yaml・
settle_nodes.txt (コメントの行を除く)・nozzle.h5 の座標・トポロジ・VALUE の全データセットの値 (補間した後の初期場)・nozzle.msh・prepare_info.json の全キー。
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402

SEPARATE = ("prepare_info:stop_rule.explicit", "prepare_info:stop_rule.requested_source", "prepare_info:ic_snapshot.how", "prepare_info:ic_from")


def flat(d, pre):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{pre}.{k}" if not pre.endswith(":") else f"{pre}{k}"))
    elif isinstance(d, list) and any(isinstance(x, (dict, list)) for x in d):
        for i, v in enumerate(d):
            out.update(flat(v, f"{pre}[{i}]"))
    else:
        out[pre] = d
    return out


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def keys_of(run: Path) -> dict:
    out = {}
    lc = RA.staged_ns_launch_configs(run, "full")
    for s in lc["stages"]:
        out[f"stage:{s['tag']}"] = sha_bytes(s["config"].encode())
        out[f"stage:{s['tag']}.nsteps"] = s["nsteps"]
    for f in ("bcondConfig.yaml", "solverConfig_prestage.yaml", "nozzle.msh"):
        out[f"file:{f}"] = sha_bytes((run / f).read_bytes())
    out.update(flat(yaml.safe_load((run / "solverConfig.yaml").read_text()), "solverConfig:"))
    out.update(flat(yaml.safe_load((run / "settle_criteria.yaml").read_text()), "settle_criteria:"))
    out["settle_nodes"] = sha_bytes("\n".join(ln for ln in (run / "settle_nodes.txt").read_text().splitlines()
                                              if ln.strip() and not ln.lstrip().startswith("#")).encode())
    hh = RA.h5_mesh_hashes(run / "nozzle.h5")
    out["nozzle_h5:coord"], out["nozzle_h5:topology"] = hh["coord_sha256"], hh["topology_sha256"]
    with h5py.File(run / "nozzle.h5", "r") as h:
        for k in sorted(h["VALUE"]):
            a = np.ascontiguousarray(np.asarray(h["VALUE/" + k][:]))
            out[f"nozzle_h5:VALUE/{k}"] = sha_bytes(a.tobytes()) + f"|{a.dtype}|{a.shape}"
    out.update(flat(json.loads((run / "prepare_info.json").read_text()), "prepare_info:"))
    return out


def diff(a: dict, b: dict) -> list:
    return sorted(k for k in set(a) | set(b) if json.dumps(a.get(k, "<無>"), sort_keys=True, default=str) != json.dumps(b.get(k, "<無>"), sort_keys=True, default=str))


ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=("noise", "ab"))
ap.add_argument("run1")
ap.add_argument("run2")
ap.add_argument("--noise", default=None)
ap.add_argument("--out", required=True)
a = ap.parse_args()
k1, k2 = keys_of(Path(a.run1)), keys_of(Path(a.run2))
d = diff(k1, k2)
if a.mode == "noise":
    res = {"mode": "noise", "run1": a.run1, "run2": a.run2, "n_keys": len(k1), "volatile": d,
           "values": {k: [k1.get(k), k2.get(k)] for k in d}}
    print(f"揺らぎの欄 {len(d)} (全 {len(k1)} キー): {d}")
else:
    noise = set(json.loads(Path(a.noise).read_text())["volatile"])
    sep = [k for k in d if k in SEPARATE]
    rest = [k for k in d if k not in noise and k not in SEPARATE]
    res = {"mode": "ab", "A": a.run1, "B": a.run2, "n_keys": len(k1), "noise_excluded": sorted(set(d) & noise),
           "separately_reported": {k: [k1.get(k), k2.get(k)] for k in sep}, "remaining_diffs": {k: [k1.get(k), k2.get(k)] for k in rest},
           "verdict": "PASS" if not rest else "FAIL"}
    print(f"全 {len(k1)} キー、差 {len(d)} (揺らぎ {len(res['noise_excluded'])}・別に記す {len(sep)}・残り {len(rest)})")
    for k in rest[:20]:
        print("  残り:", k, str(res["remaining_diffs"][k])[:160])
    print(f"VERDICT: {res['verdict']}")
Path(a.out).write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))
