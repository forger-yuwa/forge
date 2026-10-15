#!/usr/bin/env python3
r"""case/71 の NS の延長 `run_0003_ns_cont` を作って回す (plan tooling-nozzle-card-local-mesh §5.1 #3d、
ユーザ決定 2026-10-15「どうせやるなら全部doubleかな」「thermoだけとは限らないよね」、条件は codex 諮問
`notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md` の「延長へ進む場合」)。

run_0002_ns の res_183000.h5 の保存量を**同じ格子** (run_0002_ns の nozzle.h5 の複製) に index restart し
(`restart_field.py --keep-src-dtype`、前段は回さない)、本段だけを `runner_axismach.run_staged_ns(stages="none", stop_rule=True)`
で停止規則の見張りごと回す。

run_0002_ns の本段の solverConfig・bcond から変えるのは次だけ (`prepare` が YAML を読み直して差を検査する):
  - `physProp.thermoFloat: 0`                 (TP の状態の逆算を double Newton に; 既定 1 は float Newton + double 研磨)
  - `time.deltaT.implicitSolvePrecision: 1`   (block-DPLUR の左辺の solve を double に; 条件 timeIntegration 11・blockDPLUR 1・
                                               lowMachPrecond < 2 を検査する。`solverConfig.cpp:1396-1420`)
  - `time.outStepInterval: 2000`              (全場・壁・出口の出力の間隔)
  - 起動の環境変数 `FORGE_TRANSPORT_TABLE=0` (輸送物性を表でなく double の `transport_mix_Y` で; `thermo_d.cu:165`)・
    `FORGE_DIAG_FACE_H_DOUBLE=1` (SLAU の TP 多成分の面の全エンタルピーを double の `thermo_h_mix` で; `convectiveFlux_slau_d.inc.cuh:397`)。
    `run` はこの 2 つと `FORGE_CUDA_BLOCKSIZE=128`・`FORGE_BIN` (FP64 のビルド、sha256 を照合) が無ければ起動しない。
  - 停止規則: 明示 (`evaluate.stop_rule: true` 相当、使えなければ起動しない)・窓 2 万・保持 4 万・上限 20 万 step (この run の
    step 数 = 追加分; 到達の予測ではなく予算の上限)・出力 2,000。Euler の参照は `run_0001_euler/res_12000.h5` を sha256 で固定。
    帯はこの run の最初の出力から取る (`stop_rule.fix_band`) ので θ_r は run_0002_ns とは別の系列。
残る float (この run でも float のまま): 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)、
リミッタの長さの `sqrtf` (`limiter_d.cu:398`)、XMF の型の宣言 (表示だけ)。

使い方 (AWS A):
  cd ~/forge-cardmesh/design
  python3 ../case/71.m55_combustion_wt/_runs/run_0003_ns_cont/cont_0003.py prepare \
      --src-run ../case/71.m55_combustion_wt/run_0002_ns --step 183000 \
      --run-dir ../case/71.m55_combustion_wt/run_0003_ns_cont \
      --euler-run ../case/71.m55_combustion_wt/run_0001_euler --euler-res res_12000.h5 --euler-sha256 <64 桁> \
      --forge-bin ~/forge-integ-dflow2-fp64/solver_density_cuda/build/forge --forge-sha256 <64 桁>
  FORGE_BIN=... FORGE_CUDA_BLOCKSIZE=128 FORGE_TRANSPORT_TABLE=0 FORGE_DIAG_FACE_H_DOUBLE=1 \
      nohup python3 ../case/71.m55_combustion_wt/_runs/run_0003_ns_cont/cont_0003.py run \
      --run-dir ../case/71.m55_combustion_wt/run_0003_ns_cont > <log> 2>&1 < /dev/null &
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]          # リポジトリ直下 (case/71.../_runs/run_0003_ns_cont/ の 4 つ上)
sys.path.insert(0, str(ROOT / "design"))
TOOLS = ROOT / "solver_density_cuda" / "tools"

PLAN = "plans/active/tooling-nozzle-card-local-mesh.md §5.1 #3d"
RECORD = "CONTINUATION.json"
# run_0002_ns から写すもの (格子・境界・物性の記録・形状の記録)。nozzle.h5 は写した後に保存量を restart で上書きする
COPY = ["nozzle.h5", "nozzle.xmf", "bcondConfig.yaml", "probe.yaml", "species_meta.yaml",
        "resolved_species_fc8ba1a8828147e3.yaml", "MESH_QUALITY.txt",
        "wall_design.csv", "wall_physical.csv", "wall_repr.json", "target_axis_M.csv",
        "delta_r_initial.csv", "delta_r_initial.json"]
ENV_REQUIRED = {"FORGE_TRANSPORT_TABLE": "0", "FORGE_DIAG_FACE_H_DOUBLE": "1", "FORGE_CUDA_BLOCKSIZE": "128"}
# 本段の solverConfig の変更 (YAML のパス → 値)。これ以外の値は run_0002_ns の本段と同じでなければならない
CHANGES = {("physProp", "thermoFloat"): 0, ("time", "deltaT", "implicitSolvePrecision"): 1, ("time", "outStepInterval"): 2000}


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


def edit_config(text: str) -> str:
    """run_0002_ns の本段の solverConfig (flow 形式) に CHANGES を書く。新しいキーは flow の mapping の末尾に足し、
    outStepInterval は値のトークンだけを置き換える。書いた後に YAML を読み直して、差が CHANGES ちょうどであることを確かめる。"""
    import yaml
    from forge_design.evaluate.runner_axismach import _cfg_set
    before = yaml.safe_load(text)
    fb = flatten(before)
    # 足すキーは元に無いこと (重複キーを作らない; 既定 thermoFloat 1・implicitSolvePrecision 0 のまま回っていた)・出力は 1,000
    if ("physProp", "thermoFloat") in fb or ("time", "deltaT", "implicitSolvePrecision") in fb or fb.get(("time", "outStepInterval")) != 1000:
        raise SystemExit("run_0002_ns の本段の solverConfig が想定 (thermoFloat・implicitSolvePrecision 無し、outStepInterval 1000) と違う")
    for old, new in (("thermoHrefTemp: 298.15}", "thermoHrefTemp: 298.15, thermoFloat: 0}"),
                     ("detectNaN: 1}", "detectNaN: 1, implicitSolvePrecision: 1}")):
        if text.count(old) != 1:
            raise SystemExit(f"solverConfig に '{old}' がちょうど 1 回ない ({text.count(old)} 回) — 書き換えない")
        text = text.replace(old, new)
    text = _cfg_set(text, {"outStepInterval": 2000})
    after = yaml.safe_load(text)
    fa = flatten(after)
    diff = {k for k in set(fb) | set(fa) if fb.get(k, "<無>") != fa.get(k, "<無>")}
    if diff != set(CHANGES) or any(fa[k] != v for k, v in CHANGES.items()):
        raise SystemExit(f"solverConfig の差が要求と違う: {sorted(diff)} (要求 {sorted(CHANGES)})")
    return text


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
    src = Path(a.src_run).resolve()
    run = Path(a.run_dir).resolve()
    res = src / f"res_{a.step}.h5"
    if not res.is_file():
        raise SystemExit(f"{res} が無い")
    fb = Path(a.forge_bin).expanduser().resolve()
    if sha256(fb) != a.forge_sha256:
        raise SystemExit(f"{fb} の sha256 が指定と違う")
    eu_run = Path(a.euler_run).resolve()
    eref = {"run": str(eu_run), "res": a.euler_res, "sha256": a.euler_sha256}
    SR.check_euler_ref(eref)                                   # 無い・sha が違えば ValueError
    run.mkdir(parents=False, exist_ok=False)                   # 既にあれば失敗 (既存 run を上書きしない)
    for name in COPY:
        shutil.copy2(src / name, run / name)
    # 本段の設定
    cfg0 = (src / "solverConfig.yaml").read_text()
    cfg = edit_config(cfg0)
    (run / "solverConfig.yaml").write_text(cfg)
    cond = check_implicit_conditions(yaml.safe_load(cfg))
    # 停止規則 (明示)・Euler の参照の固定・この run の段の記録
    info = json.loads((src / "prepare_info.json").read_text())
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_cap": int(a.cap), "stop_rule_hold_steps": int(a.hold)},
                             nsteps_arg=None, n_main=int(a.cap), n_main_source=f"延長の予算の上限 ({PLAN})", out_int=2000,
                             euler_ref=eref)
    if not (sr["use"] and sr["explicit"]):
        raise SystemExit(f"停止規則を明示で使えない: {sr.get('reason')}")
    info["stop_rule"] = sr
    info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}
    info["continuation"] = {
        "plan": PLAN, "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "from_run": str(src), "from_res": res.name, "from_res_sha256": sha256(res), "from_step": int(a.step),
        "restart": "restart_field.py --keep-src-dtype (同じ格子の保存量の index コピー、前段は回さない)",
        "config_changes": {".".join(k): v for k, v in CHANGES.items()},
        "env_required": ENV_REQUIRED, "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
        "note": ("prepare_info の他の項目 (nStepOuter・ns_stage_configs ほか) は run_0002_ns の準備の時のまま。"
                 "この run の本段の設定は MAIN_STAGE_STARTED と stage_manifest.json を見る"),
    }
    (run / "prepare_info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    # 保存量の restart (同じ格子: MESH/COORD の一致・保存量のビット一致を restart_field が検査する)
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
    r = subprocess.run([sys.executable, str(TOOLS / "restart_field.py"), str(res), str(run / "nozzle.h5"), "--keep-src-dtype",
                        "--dst-run", str(run), "--forge", str(fb)], env=env, capture_output=True, text=True)
    (run / "restart_field.log").write_text(r.stdout + r.stderr)
    print(r.stdout + r.stderr)
    want = "VERDICT: OK (9 量を移した、SRC とビット一致)"
    if r.returncode != 0 or want not in r.stdout:
        raise SystemExit(f"restart_field が要求どおりでない (rc {r.returncode}; 期待 '{want}')")
    rec = {"plan": PLAN, "prepared": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "git_head": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
           "from": {"run": str(src), "res": res.name, "sha256": info["continuation"]["from_res_sha256"]},
           "nozzle_h5_sha256_after_restart": sha256(run / "nozzle.h5"),
           "solverConfig_sha256": sha256(run / "solverConfig.yaml"), "solverConfig_src_sha256": sha256(src / "solverConfig.yaml"),
           "bcondConfig_sha256": sha256(run / "bcondConfig.yaml"), "config_changes": info["continuation"]["config_changes"],
           "implicit_conditions": cond, "stop_rule": {k: sr[k] for k in ("explicit", "cap", "hold_steps", "window", "out_step_interval",
                                                                          "euler_ref", "use")},
           "env_required": ENV_REQUIRED, "forge_bin": str(fb), "forge_sha256": a.forge_sha256,
           "restart_verdict": want, "copied": COPY}
    (run / RECORD).write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print(json.dumps(rec, indent=1, ensure_ascii=False))


def run_main(a) -> int:
    import yaml
    run = Path(a.run_dir).resolve()
    rec = json.loads((run / RECORD).read_text())
    bad = {k: os.environ.get(k) for k, v in ENV_REQUIRED.items() if os.environ.get(k) != v}
    fb = os.environ.get("FORGE_BIN")
    if bad or not fb:
        raise SystemExit(f"環境変数が要求と違う: {bad} / FORGE_BIN={fb!r} (要求 {ENV_REQUIRED} と FORGE_BIN) — 起動しない")
    if sha256(fb) != rec["forge_sha256"]:
        raise SystemExit(f"FORGE_BIN {fb} の sha256 が準備の記録と違う — 起動しない")
    if sha256(run / "solverConfig.yaml") != rec["solverConfig_sha256"]:
        raise SystemExit("solverConfig.yaml が準備の後に変わっている — 起動しない")
    if sha256(run / "nozzle.h5") != rec["nozzle_h5_sha256_after_restart"]:
        raise SystemExit("nozzle.h5 が restart の後に変わっている — 起動しない")
    cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
    check_implicit_conditions(cfg)
    if cfg["physProp"].get("thermoFloat") != 0 or cfg["time"].get("outStepInterval") != 2000:
        raise SystemExit("solverConfig の thermoFloat・outStepInterval が要求と違う — 起動しない")
    from forge_design.evaluate.runner_axismach import run_staged_ns
    print(f"[cont_0003] {datetime.now(timezone.utc).isoformat(timespec='seconds')} 本段を起動 (stages=none, stop_rule=True) "
          f"pid {os.getpid()}", flush=True)
    rc = run_staged_ns(run, stages="none", stop_rule=True)
    print(f"[cont_0003] 終了 rc {int(rc)} 停止の状態 {getattr(rc, 'stop_rule', None)}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--src-run", required=True)
    p.add_argument("--step", type=int, required=True)
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
