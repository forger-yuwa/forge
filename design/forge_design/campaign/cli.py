"""設計の実行役の入口 (plan tooling-design-orchestrator §4.8・§4.5.1、§5.1 #6c)。

    python -m forge_design.campaign.cli run <campaign.yaml> [--foreground]   driver を user service で起動する (動いていれば状態を返す)
    python -m forge_design.campaign.cli status <campaign.yaml>               台帳を読んで状態を出す (書かない)
    python -m forge_design.campaign.cli note <campaign.yaml> --path P --reason R   直した入力の変更の要求 (fixer; driver が台帳に取り込む)

- 役割は campaign の `role` (台帳を作った後は変えられない)。`fixed` は段の処理を承認の commit の checkout から回し (`gate_env.ensure_checkout`)、
  `fixer` は作業ツリーから回す。
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


def cmd_run(a) -> int:
    c = CP.load(a.campaign)
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
    print(f"driver を起動した ({unit}、段の処理は {repo} から)。状態: python -m forge_design.campaign.cli status {a.campaign}")
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
    return 0


def cmd_note(a) -> int:
    c = CP.load(a.campaign)
    p = Path(a.path).resolve()
    sha = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    f = put_request(c["ledger"], "note", path=str(p), input_sha256=sha, reason=a.reason)
    print(f"変更の要求を置いた: {f} (driver が次に台帳へ取り込む)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m forge_design.campaign.cli", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="driver を起動する (再開も同じ)")
    r.add_argument("campaign")
    r.add_argument("--foreground", action="store_true", help="driver をこのプロセスで回す (試験用; 検査を飛ばす)")
    s = sub.add_parser("status", help="状態を出す")
    s.add_argument("campaign")
    n = sub.add_parser("note", help="直した入力の変更の要求 (fixer)")
    n.add_argument("campaign")
    n.add_argument("--path", required=True)
    n.add_argument("--reason", required=True)
    a = ap.parse_args(argv)
    return {"run": cmd_run, "status": cmd_status, "note": cmd_note}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
