"""実行の計画 `execution_plan.json` (plan tooling-design-orchestrator §4.1.1、§5.1 #4b; 親 plan tooling-design-problem-campaign-recipe §4.3 の
今回の限定の実装の例外)。

1 つの run に 1 つ。3 段で作り、確定した後は上書きしない。

1. **準備の前 (入力の固定)** `fixed_inputs`: カード・問題 YAML・承認の記録・格子の msh・初期場 (`campaign.ic.resolve_ic` の戻り)・Euler の参照・
   δ の表・forge の実行ファイル・変換器の sha256 と、数値に効く環境変数の値。
2. **準備の後・起動の前 (実効の計画の確定)** `confirm`: `runner_axismach.staged_ns_launch_configs` で全段の起動直前の設定の全文を求め
   (本段の nStepOuter は停止規則の上限)、版付きの recipe (`RECIPE`) の段の列と、段の間で許す設定の差だけであることを確かめて書く。
   境界条件・停止規則の基準・節点の表・準備の後の格子と初期場 (`nozzle.h5`、元の snapshot とは別の生成物)・準備の記録の sha256 も書く。
3. **実行と照合**: `runner_axismach.run_staged_ns(plan=…)` が、最初の段の前に全段 (`assert_launch_set`)、各段の起動の直前に
   その段 (`assert_launch`) を計画と照合し、違えば forge を起動しない。終わった後 `verify_after` が段の記録・起動の記録・準備の記録を照合する。
   起動の後に走った設定を計画に書き足して一致とはしない。

試験: design/tests/run_execution_plan_tests.py (GPU・forge なし; 起動を記録するダミーの forge)。
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
from pathlib import Path

FORMAT = "forge.execution_plan/v1"
PLAN_FILE = "execution_plan.json"
PLAN_DOC = "plans/active/tooling-design-orchestrator.md §4.1.1"
# 数値に効く環境変数 (計画に値を書き、起動の記録と照合する)
NUMERIC_ENV = ("FORGE_CUDA_BLOCKSIZE",)


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_text(t: str) -> str:
    return _sha256_bytes(t.encode())


def sha256_file(p) -> str | None:
    p = Path(p)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _recipe() -> dict:
    """版付きの recipe `axismach.wt_card_chain/v1` の NS の段の部分: 段の列と、段の間で許す設定の差 (前の段と本段の差のパス)。
    物理 (physProp・turbulence・mesh) と境界条件の差は許さない (今の NS の段の生成は物理を保って次数・CFL を変える)。"""
    from ..evaluate import runner_axismach as RA
    allowed = ({("time", "deltaT", "cfl"), ("time", "deltaT", "cfl_pseudo"), ("space", "convMethod"), ("time", "nStepInner"),
                ("time", "last", "nStepOuter"), ("time", "outStepInterval"), ("time", "fieldOutStepInterval")}
               | set(RA._MAIN_LINE_PATHS) | set(RA._NS_OUTPUT_PATHS))
    return {"name": "axismach.wt_card_chain", "version": "v1",
            "ns_stages": {"full": ["S1_soft", "S2_mid", "main"], "none": ["main"]},
            "allowed_stage_diff_paths": sorted(".".join(p) for p in allowed),
            "forbidden_stage_diff_prefixes": ["physProp", "turbulence", "mesh"]}


RECIPE = None


def recipe() -> dict:
    global RECIPE
    if RECIPE is None:
        RECIPE = _recipe()
    return RECIPE


def fnv1a64(text: str) -> str:
    from ..evaluate import runner_axismach as RA
    RA._stage_manifest_cls()
    import stage_manifest as SM
    return SM.fnv1a64(text)


# ---------------------------------------------------------------------------------------------------- 1. 準備の前

def fixed_inputs(*, card, problem, approval, mesh_msh=None, ic=None, euler_ref=None, delta_table=None, forge_bin=None, converter=None,
                 expected_initializer=None) -> dict:
    """準備の前に固定する入力 (パスと sha256)。ic は `campaign.ic.resolve_ic` の戻り (初期場のファイル・継承する値・元の run の記録)。"""
    def ent(p):
        return None if p is None else {"path": str(Path(p).resolve()), "sha256": sha256_file(p)}
    fb = forge_bin or os.environ.get("FORGE_BIN")
    conv = converter or os.environ.get("FORGE_CONVERTER")
    out = {"card": ent(card), "problem": ent(problem), "approval": ent(approval), "mesh_msh": ent(mesh_msh), "ic": ic,
           "euler_ref": ent(euler_ref), "delta_table": ent(delta_table), "forge_bin": ent(fb), "converter": ent(conv),
           "env": {k: os.environ.get(k) for k in NUMERIC_ENV}, "fixed_at": _dt.datetime.now().isoformat(timespec="seconds"),
           "expected_initializer": dict(expected_initializer or {})}
    miss = [k for k in ("card", "problem", "approval", "forge_bin") if out[k] is None or out[k]["sha256"] is None]
    if miss:
        raise ValueError(f"計画の入力が無い・読めない {miss}")
    return out


# ---------------------------------------------------------------------------------------------------- 2. 準備の後・起動の前

def initializer_problems(expected: dict, ini: dict) -> list:
    """準備の記録の initializer に、計画で期待した初期壁の係数 (カードの `deltastar_initializer.*`) が入っているか。型 (数値)・有限・値
    (相対 1e-12) を確かめる。違い・欠け・非有限の列 (旧の台本 p5_ns.py の準備の時の検査; codex 2026-10-17 で NaN を通していた比較式を直した)。"""
    import math
    probs = []
    for k, v in (expected or {}).items():
        got = (ini or {}).get(k)
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            if got != v:
                probs.append(f"{k}: 準備の記録 {got!r} ≠ 期待 {v!r}")
            continue
        if isinstance(got, bool) or not isinstance(got, (int, float)):
            probs.append(f"{k}: 準備の記録に数値が無い ({got!r})")
        elif not (math.isfinite(float(got)) and math.isfinite(float(v))):
            probs.append(f"{k}: 有限でない (準備の記録 {got!r}、期待 {v!r})")
        elif abs(float(got) - float(v)) > 1e-12 * max(1.0, abs(float(v))):
            probs.append(f"{k}: 準備の記録 {got!r} ≠ 期待 {v!r}")
    return probs


def _diff_paths(a_text: str, b_text: str) -> list:
    from ..evaluate import runner_axismach as RA
    d = RA._cfg_diff_paths(RA._cfg_load(a_text), RA._cfg_load(b_text))
    return sorted(".".join(p) for p in d)


def stage_diff_problems(stages: list, rec: dict) -> list:
    """前の段と本段の差が recipe の許す差だけか (物理と境界条件の差は許さない)。"""
    probs = []
    main = next((s for s in stages if s["tag"] == "main"), None)
    if main is None:
        return ["本段 main が無い"]
    allowed = set(rec["allowed_stage_diff_paths"])
    for s in stages:
        if s["tag"] == "main":
            continue
        diffs = _diff_paths(s["config"], main["config"])
        bad = [d for d in diffs if not any(d == a or d.startswith(a + ".") for a in allowed)]
        phys = [d for d in diffs if d.split(".")[0] in rec["forbidden_stage_diff_prefixes"]]
        if bad or phys:
            probs.append(f"段 {s['tag']} と本段の差に recipe が許さないパス {sorted(set(bad) | set(phys))}")
    return probs


def confirm(run_dir, inputs: dict, *, stages: str = "full", attempt_id: str | None = None) -> dict:
    """準備の後・起動の前に計画を確定して run_dir/execution_plan.json に書く (既にあれば止める)。"""
    from ..evaluate import runner_axismach as RA
    run_dir = Path(run_dir)
    out = run_dir / PLAN_FILE
    if out.exists():
        raise ValueError(f"{out} が既にある — 確定した計画は上書きしない (新しい試行は新しい run で)")
    rec = recipe()
    want = rec["ns_stages"].get(stages)
    if want is None:
        raise ValueError(f"recipe {rec['name']}/{rec['version']} は段の列 {stages!r} を持たない ({list(rec['ns_stages'])})")
    lc = RA.staged_ns_launch_configs(run_dir, stages)
    tags = [s["tag"] for s in lc["stages"]]
    probs = []
    if tags != want:
        probs.append(f"段の列 {tags} が recipe の {want} と違う")
    probs += stage_diff_problems(lc["stages"], rec)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ic = inputs.get("ic") or {}
    snap = (ic.get("snapshot") or {})
    if snap and (info.get("ic_snapshot") or {}).get("sha256") != snap.get("sha256"):
        probs.append(f"準備の記録の初期場 {info.get('ic_snapshot')} が計画の初期場 {snap} と違う")
    probs += [f"初期壁の係数: {x}" for x in initializer_problems(inputs.get("expected_initializer"), info.get("initializer"))]
    msh = (inputs.get("mesh_msh") or {}).get("sha256")
    rec_msh = ((info.get("mesh") or {}).get("hashes") or {}).get("msh_sha256")
    if msh and rec_msh != msh:
        probs.append(f"準備の記録の格子 (msh_sha256 {str(rec_msh)[:16]}…) が計画の格子 {msh[:16]}… と違う")
    if probs:
        raise ValueError("計画を確定しない: " + "; ".join(probs))
    stages_rec = [{"order": i, "tag": s["tag"], "nsteps": s["nsteps"], "sha256": sha256_text(s["config"]), "fnv": fnv1a64(s["config"]),
                   "config": s["config"]} for i, s in enumerate(lc["stages"])]
    doc = {"format": FORMAT, "plan": PLAN_DOC, "recipe": {k: rec[k] for k in ("name", "version")}, "attempt_id": attempt_id,
           "run_dir": str(run_dir.resolve()), "inputs": inputs,
           "confirmed": {"at": _dt.datetime.now().isoformat(timespec="seconds"), "stages_mode": stages, "stages": stages_rec,
                         "bcond_sha256": sha256_text(lc["bcond"]), "stop_rule_cap": int(lc["stop_rule"]["cap"]),
                         "files_sha256": {f: sha256_file(run_dir / f) for f in ("settle_criteria.yaml", "settle_nodes.txt", "nozzle.h5",
                                                                                "prepare_info.json", "species_meta.yaml")}},
           "recipe_rules": {k: rec[k] for k in ("allowed_stage_diff_paths", "forbidden_stage_diff_prefixes")}}
    doc["sha256"] = sha256_text(json.dumps({k: v for k, v in doc.items() if k != "sha256"}, sort_keys=True, ensure_ascii=False, default=str))
    out.write_text(json.dumps(doc, indent=1, ensure_ascii=False, default=str))
    return doc


# ---------------------------------------------------------------------------------------------------- 3. 実行と照合

def load_confirmed(path) -> dict:
    p = Path(path)
    if p.is_dir():
        p = p / PLAN_FILE
    doc = json.loads(p.read_text())
    if doc.get("format") != FORMAT or "confirmed" not in doc:
        raise ValueError(f"{p} は確定した実行の計画でない")
    body = {k: v for k, v in doc.items() if k != "sha256"}
    if sha256_text(json.dumps(body, sort_keys=True, ensure_ascii=False, default=str)) != doc.get("sha256"):
        raise ValueError(f"{p} の中身と sha256 が合わない (確定した計画を書き換えた?)")
    return doc


# 最初の段の前に、計画の確定の後に変わっていないことを確かめるファイル (旧の台本 p5_ns.py の「準備の後に変わっていない」の検査; §5.1 #6d)。
# nozzle.h5 は段の間の restart で書き換わるので最初の段の前だけ
PRELAUNCH_FILES = ("settle_criteria.yaml", "settle_nodes.txt", "nozzle.h5", "prepare_info.json", "species_meta.yaml")


def assert_launch_set(doc: dict, planned: list, bc_text: str, run_dir=None) -> None:
    """起動しようとする全段 [(tag, 設定)] と境界条件が計画と同じか、run_dir を渡せば計画の確定の後に変わってはいけないファイル
    (停止規則の基準・節点の表・格子と初期場・準備の記録) が同じか。違えば ValueError (forge を 1 回も起動しない)。"""
    st = doc["confirmed"]["stages"]
    probs = []
    if run_dir is not None:
        for f in PRELAUNCH_FILES:
            want = (doc["confirmed"].get("files_sha256") or {}).get(f)
            got = sha256_file(Path(run_dir) / f)
            if want != got:
                probs.append(f"{f} が計画の確定の後に変わった (sha256 {str(got)[:12]} ≠ 計画 {str(want)[:12]})")
    if [t for t, _ in planned] != [s["tag"] for s in st]:
        probs.append(f"段の列 {[t for t, _ in planned]} が計画の {[s['tag'] for s in st]} と違う")
    for (t, c), s in zip(planned, st):
        if sha256_text(c) != s["sha256"]:
            probs.append(f"段 {t} の設定が計画と違う (差 {_diff_paths(s['config'], c)[:8]})")
    if sha256_text(bc_text) != doc["confirmed"]["bcond_sha256"]:
        probs.append("境界条件が計画と違う")
    if probs:
        raise ValueError("計画と違う — forge を起動しない: " + "; ".join(probs))


def assert_launch(doc: dict, tag: str, cfg: str, bc_text: str) -> None:
    """その段の起動の直前の照合。"""
    s = next((x for x in doc["confirmed"]["stages"] if x["tag"] == tag), None)
    if s is None:
        raise ValueError(f"段 {tag} は計画に無い — forge を起動しない")
    if sha256_text(cfg) != s["sha256"]:
        raise ValueError(f"段 {tag} の設定が計画と違う (差 {_diff_paths(s['config'], cfg)[:8]}) — forge を起動しない")
    if sha256_text(bc_text) != doc["confirmed"]["bcond_sha256"]:
        raise ValueError(f"段 {tag}: 境界条件が計画と違う — forge を起動しない")


def verify_after(run_dir) -> dict:
    """走った後の照合: 段の記録 (段の列・設定の FNV)、起動の記録 (forge の sha256・数値に効く環境変数)、準備の記録 (初期場・格子)。"""
    run_dir = Path(run_dir)
    doc = load_confirmed(run_dir)
    probs = []
    try:
        sm = json.loads((run_dir / "stage_manifest.json").read_text())
    except (OSError, ValueError) as e:
        sm = None
        probs.append(f"段の記録を読めない ({e})")
    if sm is not None:
        got = [(s.get("tag"), s.get("cfg_fnv")) for s in sm.get("stages") or []]
        want = [(s["tag"], s["fnv"]) for s in doc["confirmed"]["stages"]]
        if got != want:
            probs.append(f"段の記録 {got} が計画 {want} と違う")
    prov = {}
    try:
        for ln in (run_dir / "RUN_PROVENANCE.txt").read_text().splitlines():
            k, _, v = ln.partition(":")
            prov.setdefault(k.strip(), []).append(v.strip())
    except OSError as e:
        probs.append(f"起動の記録を読めない ({e})")
    fb = (doc["inputs"].get("forge_bin") or {}).get("sha256")
    if prov and fb and (prov.get("forge_sha256") or [None])[0] != fb:
        probs.append(f"起動の記録の forge の sha256 {(prov.get('forge_sha256') or [None])[0]} が計画 {fb} と違う")
    for k, v in (doc["inputs"].get("env") or {}).items():
        if v is not None and prov and f"{k}={v}" not in (prov.get("env") or []):
            probs.append(f"起動の記録に環境変数 {k}={v} が無い")
    info = json.loads((run_dir / "prepare_info.json").read_text())
    snap = ((doc["inputs"].get("ic") or {}).get("snapshot") or {}).get("sha256")
    if snap and (info.get("ic_snapshot") or {}).get("sha256") != snap:
        probs.append("準備の記録の初期場が計画と違う")
    return {"status": "OK" if not probs else "MISMATCH", "problems": probs}
