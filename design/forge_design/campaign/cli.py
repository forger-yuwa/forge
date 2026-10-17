"""設計の実行役の入口 `forge-design` (plan tooling-design-orchestrator §4.8・§4.5.1、§5.1 #6c・#9)。`design/bin/forge-design` が
`python -m forge_design.campaign.cli` を呼ぶ。

    forge-design run <campaign.yaml> --role fixed|fixer [--resume] [--foreground]
        driver を user service で起動する (動いていれば状態を返す)。台帳が無ければ新しく始め (--resume を付けない)、
        台帳があれば --resume が要る (付け忘れで別の campaign を上書きしない)。--role は campaign の role と同じでなければ止める
        (役割は起動のときに人が決める; モデルの自己申告で変えない、§4.0)。
    forge-design status <campaign.yaml>
        台帳を読んで状態と集計 (試行の数・人に返した回数・直した回数・消費した予算) を出す (書かない)。
    forge-design note <campaign.yaml> --path P --reason R
        直した入力の変更の要求 (fixer; driver が台帳に取り込む)。

- 終了コード: run は 0 起動した・動いている・(--foreground で) 終わった、1 起動できない・止まった、2 起動の前の検査・引数の照合で止めた。
  status・note は 0 (台帳が無くても 0)。所要: run・status・note は数秒 (run の driver は service で回り、段ごとに数分〜数時間)。
- 役割は campaign の `role` (台帳を作った後は変えられない)。`fixed` は段の処理を承認の commit の checkout から回し (`gate_env.ensure_checkout`)、
  `fixer` は作業ツリーから回す。AWS の実行担当 (`--executor aws`) は plan §5.1 #7 の後。
- 起動の前の検査: systemd の user service と `Linger=yes` (ssh を切っても止まらない)。`--foreground` は検査を飛ばして driver をこのプロセスで回す (試験用)。
- 走らせる環境 (FORGE_BIN・FORGE_CONVERTER・FORGE_CUDA_BLOCKSIZE・FREECADCMD・FORGE_MESH_PY・FORGE_PPTX_PY) は今のシェルから driver に渡し、driver が初回に
  台帳へ写す (再開で違えば止まる)。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from . import campaign as CP
from .ledger import put_request, read_ledger

DESIGN_DIR = Path(__file__).resolve().parents[2]
REPO_DIR = DESIGN_DIR.parent
ENV_KEYS = ("FORGE_BIN", "FORGE_CONVERTER", "FORGE_CUDA_BLOCKSIZE", "FREECADCMD", "FORGE_MESH_PY", "FORGE_PPTX_PY")


def precheck() -> list:
    probs = []
    r = subprocess.run(["systemctl", "--user", "is-system-running"], capture_output=True, text=True)
    if (r.stdout or "").strip() not in ("running", "degraded"):
        probs.append(f"systemd の user service が使えない ({(r.stdout or r.stderr).strip()})")
    lg = subprocess.run(["loginctl", "show-user", os.environ.get("USER", ""), "-p", "Linger"], capture_output=True, text=True)
    if (lg.stdout or "").strip() != "Linger=yes":
        probs.append(f"{(lg.stdout or lg.stderr).strip()}: ssh を切ると user service が止まる — `sudo loginctl enable-linger $USER` を一度 (人に確かめてから)")
    return probs


def repo_for(campaign: dict) -> Path:
    if campaign["role"] == "fixed":
        from ..feedback import gate_env as GE
        apr = json.loads(Path(campaign["approval"]["file"]).read_text())
        return GE.ensure_checkout(apr["gate"]["commit"])
    return REPO_DIR


def driver_unit(campaign: dict) -> str:
    return "forge-campaign-" + hashlib.sha256(str(Path(campaign["ledger"]).resolve()).encode()).hexdigest()[:12]


def run_args_problems(c: dict, role: str, resume: bool) -> list:
    """run の引数と campaign・台帳の照合 (plan §4.0・§4.8)。問題の列 (空なら起動してよい)。"""
    probs = []
    if role != c["role"]:
        probs.append(f"--role {role} が campaign の role {c['role']} と違う (役割は campaign に固定; 変えるなら新しい campaign)")
    exists = read_ledger(c["ledger"]) is not None
    if exists and not resume:
        probs.append(f"台帳 {c['ledger']} が既にある — 再開なら --resume を付ける (新しく始めるなら別の ledger)")
    if resume and not exists:
        probs.append(f"--resume だが台帳 {c['ledger']} が無い — 新しく始めるなら --resume を外す")
    return probs


def summary(d: dict) -> dict:
    """台帳からの集計 (plan §4.8): 試行の数 (状態ごと)・人に返した回数 (driver が止まった記録 stops)・直した回数 (変更の要求 note)・
    消費した予算。"""
    by_state = {}
    for x in d.get("attempts") or []:
        by_state[x.get("state")] = by_state.get(x.get("state"), 0) + 1
    reqs = d.get("requests") or []
    notes = [r for r in reqs if r.get("kind") == "note"]
    return {"attempts": len(d.get("attempts") or []), "attempts_by_state": by_state,
            "returned_to_human": len(d.get("stops") or []), "fixes": len(notes),
            "fix_reasons": [str(r.get("reason"))[:80] for r in notes],
            "budget_consumed": dict((d.get("budget") or {}).get("consumed") or {}),
            "budget_approved": dict((d.get("budget") or {}).get("approved") or {})}


def cmd_run(a) -> int:
    c = CP.load(a.campaign)
    probs = run_args_problems(c, a.role, a.resume)
    if probs:
        print("止めた — 引数の照合:\n- " + "\n- ".join(probs), file=sys.stderr)
        return 2
    repo = repo_for(c)
    design = repo / "design"
    if a.foreground:
        sys.path.insert(0, str(design))
        from forge_design.campaign.driver import Driver, recipe_for
        r = Driver(a.campaign, recipe_for(c, repo), design_dir=design).run()
        print(json.dumps({k: v for k, v in r.items() if k != "ledger"}, ensure_ascii=False))
        return {"done": 0, "busy": 0}.get(r["status"], 1)
    probs = precheck()
    if probs:
        print("止めた — 起動の前の検査:\n- " + "\n- ".join(probs), file=sys.stderr)
        return 2
    unit = driver_unit(c)
    st = subprocess.run(["systemctl", "--user", "is-active", unit], capture_output=True, text=True).stdout.strip()
    if st in ("active", "activating"):
        print(f"driver は動いている ({unit}) — 状態は status で")
        return 0
    envs = [f"--setenv={k}={os.environ[k]}" for k in ENV_KEYS if os.environ.get(k)]
    r = subprocess.run(["systemd-run", "--user", "--quiet", f"--unit={unit}", "--collect", "-p", "ExitType=main", "-p", "Restart=no",
                        f"--working-directory={design}", f"--setenv=PYTHONPATH={design}", "--setenv=PYTHONDONTWRITEBYTECODE=1", *envs,
                        sys.executable, "-m", "forge_design.campaign.driver", str(Path(a.campaign).resolve())], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"driver を起動できない: {r.stderr.strip()[-300:]}", file=sys.stderr)
        return 1
    print(f"driver を起動した ({unit}、段の処理は {repo} から)。状態: forge-design status {a.campaign}")
    return 0


def cmd_status(a) -> int:
    c = CP.load(a.campaign)
    d = read_ledger(c["ledger"])
    if d is None:
        print("台帳はまだ無い (run で始める)")
        return 0
    b = d["budget"]
    print(f"campaign {d['name']} ({d['campaign_id']})、役割 {d['role']}、recipe {d['recipe']['name']}/{d['recipe']['version']}")
    print(f"予算: GPU 時間 消費 {b['consumed']['gpu_hours']:.3f} / 承認 {b['approved'].get('gpu_hours')}、ディスク {b['consumed']['disk_gb']:.2f} / {b['approved'].get('disk_gb')} GB")
    for x in d["attempts"]:
        it = x.get("interpretation") or {}
        print(f"  {x['id']:<24} {x['state']:<9} {'OK ' if it.get('ok') else '   '} {str(it.get('detail') or '')[:120]}")
    sm = summary(d)
    print(f"集計: 試行 {sm['attempts']} ({', '.join(f'{k} {v}' for k, v in sm['attempts_by_state'].items()) or '無し'})・"
          f"人に返した {sm['returned_to_human']} 回・直した {sm['fixes']} 回")
    for st_ in (d.get("stops") or [])[-3:]:
        print(f"  止まった {st_.get('time')} {st_.get('status')}: {str(st_.get('why'))[:120]}")
    if a.json:
        print(json.dumps(sm, ensure_ascii=False))
    return 0


def cmd_note(a) -> int:
    c = CP.load(a.campaign)
    p = Path(a.path).resolve()
    sha = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    f = put_request(c["ledger"], "note", path=str(p), input_sha256=sha, reason=a.reason)
    print(f"変更の要求を置いた: {f} (driver が次に台帳へ取り込む)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="forge-design", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="driver を起動する (台帳があれば --resume)")
    r.add_argument("campaign", help="campaign の YAML (forge.campaign/v1: 承認の記録・recipe・役割・台帳の場所)")
    r.add_argument("--role", required=True, choices=("fixed", "fixer"), help="役割 (campaign の role と同じであること)")
    r.add_argument("--resume", action="store_true", help="既にある台帳から再開する")
    r.add_argument("--foreground", action="store_true", help="driver をこのプロセスで回す (試験用; systemd の検査を飛ばす)")
    s = sub.add_parser("status", help="状態と集計を出す (書かない)")
    s.add_argument("campaign")
    s.add_argument("--json", action="store_true", help="集計を JSON でも出す")
    n = sub.add_parser("note", help="直した入力の変更の要求 (fixer)")
    n.add_argument("campaign")
    n.add_argument("--path", required=True)
    n.add_argument("--reason", required=True)
    a = ap.parse_args(argv)
    return {"run": cmd_run, "status": cmd_status, "note": cmd_note}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
