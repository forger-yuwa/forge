"""判定専用の環境 (plan tooling-design-orchestrator §4.7.1 (5)、§5.1 #2b)。

判定 (`wt_qualify --approval`) を、承認の記録に書いた版 (commit) の別の checkout と、承認のときに記録した Python・外部の道具で回す。
作業側の作業ツリー・環境変数 (`PYTHONPATH`・`FORGE_*` の道具の選択) は判定に効かない。

- **checkout**: `~/.cache/forge-gates/<commit>` に `git worktree add --detach` で作る (`design`・`solver_density_cuda/tools` だけの sparse checkout)。判定の前に、HEAD が承認の commit であること、
  追跡しているファイルに変更が無いこと、コードの探索先 (`design/forge_design`・`solver_density_cuda/tools`) に未追跡・ignored のファイルが
  無いことを確かめる (Python は -B で回し、バイトコードを書かない)。
- **Python の指紋** (`python_fingerprint`): 実行ファイルの実体と sha256、版、`sys.prefix`、配布物ごとの名前・版・RECORD の sha256・editable の有無。
  editable の配布物があれば止める (作業ツリーを参照しうる)。`pip freeze` は使わない (別の Python・外部の実行ファイルを捉えない)。
- **外部の道具** (`TOOLS`): gmsh 用の Python (`FORGE_MESH_PY`、指紋)、FreeCAD (`FREECADCMD`、実体と sha256)、変換器 (`FORGE_CONVERTER`、実体と sha256)。
  判定のときは記録したパスをその環境変数に入れて渡す。
- **起動**: 隔離のモードの Python (`-I -B`) に、checkout の `design` を sys.path の先頭に入れて `wt_qualify` を回す。環境変数は許すものだけを渡し、
  `FORGE_GATE_ENV=1`・`FORGE_GATE_CHECKOUT` を足す。`wt_qualify` は承認のモードでこの印と自分の import の元・checkout の HEAD を確かめる
  (`assert_in_gate_env`)。

使い方:
  python -m forge_design.feedback.gate_env snapshot [--python PY]            承認のときに記録する環境 (JSON)
  python -m forge_design.feedback.gate_env judge --approval card_approval.json -- <wt_qualify の引数 (--approval は自動で足す)>
  python -m forge_design.feedback.gate_env check --approval card_approval.json   環境の照合だけ (判定しない)
終了コード: judge は wt_qualify の終了コード、環境の照合に失敗したら 1。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

GATE_ROOT_ENV = "FORGE_GATE_ROOT"
DEFAULT_GATE_ROOT = Path.home() / ".cache" / "forge-gates"
CODE_DIRS = ("design/forge_design", "solver_density_cuda/tools")
CHECKOUT_PATHS = ("design", "solver_density_cuda/tools")      # sparse checkout に入れる所
MARK_ENV, CHECKOUT_ENV = "FORGE_GATE_ENV", "FORGE_GATE_CHECKOUT"
PASS_ENV = ("HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ")   # 判定に渡してよい環境変数 (道具の選択は下の TOOLS で決める)
SAFE_PATH = "/usr/local/bin:/usr/bin:/bin"

FINGERPRINT_PY = r'''
import hashlib, importlib.metadata as m, json, os, sys
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
ds = []
for d in m.distributions():
    rec = d.read_text("RECORD")
    du = d.read_text("direct_url.json")
    ed = False
    if du:
        try:
            ed = bool(json.loads(du).get("dir_info", {}).get("editable"))
        except ValueError:
            ed = None
    ds.append({"name": (d.metadata["Name"] or "").lower(), "version": d.version,
               "record_sha256": None if rec is None else hashlib.sha256(rec.encode()).hexdigest(), "editable": ed})
ds.sort(key=lambda x: (x["name"], x["version"]))
exe = os.path.realpath(sys.executable)
print(json.dumps({"executable": sys.executable, "executable_real": exe, "executable_sha256": sha(exe), "version": sys.version,
                  "prefix": sys.prefix, "base_prefix": sys.base_prefix, "dists": ds}))
'''


class GateEnvError(ValueError):
    """判定専用の環境の誤り (wt_qualify の main は ValueError として「止めた」を返す)。"""


def _sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def python_fingerprint(python) -> dict:
    r = subprocess.run([str(python), "-I", "-B", "-c", FINGERPRINT_PY], capture_output=True, text=True,
                       env={"PATH": SAFE_PATH, "HOME": str(Path.home())}, timeout=120)
    if r.returncode != 0:
        raise GateEnvError(f"Python {python} の指紋を取れない (rc {r.returncode}): {r.stderr[-400:]}")
    fp = json.loads(r.stdout)
    fp["editable_dists"] = [d["name"] for d in fp["dists"] if d["editable"]]
    return fp


def file_identity(path) -> dict | None:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(path), "exists": False}
    real = Path(os.path.realpath(p))
    return {"path": str(path), "real": str(real), "sha256": _sha256_file(real), "size": real.stat().st_size, "exists": True}


def current_tools() -> dict:
    """今の規則で選ばれる外部の道具 (承認のときに記録する)。"""
    from ..evaluate import runner_axismach as RA
    from ..export import deliverables as DV
    from ..export import wall_step as WS
    conv = None
    try:
        conv = str(RA.converter_path())
    except Exception:  # noqa: BLE001 — 見つからなければ記録しない (判定のときに成果物の再検査が止まる)
        conv = None
    return {"FORGE_MESH_PY": DV.GMSH_PYTHON, "FREECADCMD": WS.FREECADCMD, "FORGE_CONVERTER": conv}


def tool_identity(env: str, path) -> dict | None:
    if path is None:
        return None
    if env == "FORGE_MESH_PY":
        ident = file_identity(path)
        if ident and ident.get("exists"):
            ident["fingerprint"] = python_fingerprint(path)
        return ident
    return file_identity(path)


def snapshot(python=None, tools: dict | None = None) -> dict:
    """判定専用の環境の記録 (承認の記録の gate.env に入れる)。"""
    py = python or sys.executable
    tl = current_tools() if tools is None else tools
    return {"python": python_fingerprint(py), "python_path": str(py),
            "tools": {k: tool_identity(k, v) for k, v in tl.items()}}


def compare_env(recorded: dict, now: dict) -> list:
    probs = []
    a, b = recorded.get("python") or {}, now.get("python") or {}
    for k in ("executable_real", "executable_sha256", "version", "prefix", "base_prefix"):
        if a.get(k) != b.get(k):
            probs.append(f"判定の Python の {k} が違う (記録 {str(a.get(k))[:60]!r} / 今 {str(b.get(k))[:60]!r})")
    if a.get("dists") != b.get("dists"):
        an = {(d["name"], d["version"], d["record_sha256"]) for d in a.get("dists") or []}
        bn = {(d["name"], d["version"], d["record_sha256"]) for d in b.get("dists") or []}
        probs.append(f"判定の Python の配布物が違う (記録だけ {sorted(an - bn)[:5]}、今だけ {sorted(bn - an)[:5]})")
    if b.get("editable_dists"):
        probs.append(f"判定の Python に editable の配布物がある {b['editable_dists']} (作業ツリーを参照しうる)")
    ta, tb = recorded.get("tools") or {}, now.get("tools") or {}
    for k in sorted(set(ta) | set(tb)):
        x, y = ta.get(k), tb.get(k)
        if x is None and y is None:
            continue
        if (x or {}).get("sha256") != (y or {}).get("sha256") or (x or {}).get("real") != (y or {}).get("real"):
            probs.append(f"外部の道具 {k} が違う (記録 {(x or {}).get('real')} {str((x or {}).get('sha256'))[:12]} / "
                         f"今 {(y or {}).get('real')} {str((y or {}).get('sha256'))[:12]})")
        fx, fy = (x or {}).get("fingerprint"), (y or {}).get("fingerprint")
        if fx is not None or fy is not None:
            sub = compare_env({"python": fx or {}}, {"python": fy or {}})
            probs += [f"{k}: {s}" for s in sub]
    return probs


# ---------------------------------------------------------------------------------------------------- checkout

def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def gate_root() -> Path:
    return Path(os.environ.get(GATE_ROOT_ENV) or DEFAULT_GATE_ROOT)


def _git(path, *args) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(path), *args], capture_output=True, text=True)


def checkout_problems(path, commit: str) -> list:
    probs = []
    path = Path(path)
    if not path.is_dir():
        return [f"checkout {path} が無い"]
    h = _git(path, "rev-parse", "HEAD")
    if h.returncode != 0 or h.stdout.strip() != commit:
        probs.append(f"checkout の HEAD {h.stdout.strip()!r} が承認の commit {commit!r} と違う")
    t = _git(path, "status", "--porcelain", "--untracked-files=no")
    if t.returncode != 0 or t.stdout.strip():
        probs.append(f"checkout の追跡ファイルに変更がある: {t.stdout.strip()[:300]}")
    u = _git(path, "status", "--porcelain", "--ignored", "--untracked-files=all", "--", *CODE_DIRS)
    if u.returncode != 0 or u.stdout.strip():
        probs.append(f"checkout のコードの探索先に未追跡・ignored のファイルがある: {u.stdout.strip()[:300]}")
    return probs


def ensure_checkout(commit: str, repo=None) -> Path:
    if not (isinstance(commit, str) and len(commit) == 40):
        raise GateEnvError(f"承認の commit が完全な hash でない ({commit!r})")
    repo = Path(repo) if repo else repo_root()
    dst = gate_root() / commit
    if not dst.exists():
        # 判定に要る所だけの sparse checkout (全体は約 170 MB、design と solver_density_cuda/tools だけなら約 8 MB)
        dst.parent.mkdir(parents=True, exist_ok=True)
        for args, where in ((("worktree", "add", "--no-checkout", "--detach", str(dst), commit), repo),
                            (("sparse-checkout", "set", "--cone", *CHECKOUT_PATHS), dst),
                            (("checkout", "-q", "--detach", commit), dst)):
            r = _git(where, *args)
            if r.returncode != 0:
                raise GateEnvError(f"checkout を作れない ({commit}; git {args[0]}): {r.stderr.strip()[:400]}")
    probs = checkout_problems(dst, commit)
    if probs:
        raise GateEnvError("判定の checkout が使えない: " + "; ".join(probs))
    return dst


# ---------------------------------------------------------------------------------------------------- 判定の起動・印の検査

def sanitized_env(checkout: Path, tools: dict) -> dict:
    env = {k: os.environ[k] for k in PASS_ENV if k in os.environ}
    env.update(PATH=SAFE_PATH, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", MPLBACKEND="Agg", QT_QPA_PLATFORM="offscreen")
    env[MARK_ENV] = "1"
    env[CHECKOUT_ENV] = str(checkout)
    for k, ident in (tools or {}).items():
        if ident and ident.get("exists"):
            env[k] = ident["path"]
    return env


LAUNCH_PY = r'''
import runpy, sys
co = sys.argv[1]
sys.path.insert(0, co + "/design")
sys.argv = [co + "/design/forge_design/feedback/wt_qualify.py"] + sys.argv[2:]
runpy.run_module("forge_design.feedback.wt_qualify", run_name="__main__", alter_sys=True)
'''


def assert_in_gate_env(apr: dict) -> dict:
    """wt_qualify の承認のモードが判定専用の環境で回っているか (印・import の元・checkout の HEAD)。違えば GateEnvError。"""
    gate = apr.get("gate") or {}
    commit = gate.get("commit")
    if os.environ.get(MARK_ENV) != "1":
        raise GateEnvError("判定専用の環境で回していない (python -m forge_design.feedback.gate_env judge --approval … で回す)")
    co = Path(os.environ.get(CHECKOUT_ENV, "")).resolve()
    here = Path(__file__).resolve()
    if co not in here.parents:
        raise GateEnvError(f"import の元 {here} が判定の checkout {co} の中でない")
    probs = checkout_problems(co, commit)
    if probs:
        raise GateEnvError("判定の checkout が承認と違う: " + "; ".join(probs))
    return {"checkout": str(co), "commit": commit, "module": str(here)}


def check(approval_path, python=None) -> dict:
    """承認の記録の gate と、今の checkout・Python・外部の道具を照合する (判定はしない)。"""
    from ..intake import approval as AP
    apr = AP.load_approval(approval_path)
    gate = apr.get("gate") or {}
    if not gate.get("commit"):
        raise GateEnvError("承認の記録に gate.commit が無い")
    if gate.get("tracked_clean") is not True:
        raise GateEnvError("承認のときの作業ツリーに追跡ファイルの変更があった (gate.tracked_clean が true でない) — commit してから承認し直す")
    rec_env = gate.get("env")
    if not rec_env:
        raise GateEnvError("承認の記録に gate.env (判定の Python・外部の道具) が無い")
    co = ensure_checkout(gate["commit"])
    py = python or rec_env.get("python_path")
    now = snapshot(py, {k: (v or {}).get("path") for k, v in (rec_env.get("tools") or {}).items()})
    probs = compare_env(rec_env, now)
    return {"status": "OK" if not probs else "MISMATCH", "checkout": str(co), "python": py, "problems": probs}


PROBE_PY = r'''
import json, sys
co, ap = sys.argv[1], sys.argv[2]
sys.path.insert(0, co + "/design")
from forge_design.feedback import gate_env as GE
from forge_design.intake import approval as AP
import forge_design
print(json.dumps({"assert": GE.assert_in_gate_env(AP.load_approval(ap)), "forge_design": forge_design.__file__,
                  "sys_path0": sys.path[0], "executable": sys.executable}))
'''


def probe(approval_path, python=None) -> dict:
    """判定専用の環境で import の元と印の検査だけを回す (判定はしない; 試験と事前の確認用)。"""
    r = check(approval_path, python)
    if r["status"] != "OK":
        return {"status": "MISMATCH", "problems": r["problems"]}
    from ..intake import approval as AP
    apr = AP.load_approval(approval_path)
    co = Path(r["checkout"])
    env = sanitized_env(co, apr["gate"]["env"].get("tools") or {})
    q = subprocess.run([r["python"], "-I", "-B", "-c", PROBE_PY, str(co), str(Path(approval_path).resolve())], env=env, cwd=str(co),
                       capture_output=True, text=True)
    if q.returncode != 0:
        return {"status": "FAIL", "problems": [q.stderr.strip()[-600:]]}
    return {"status": "OK", **json.loads(q.stdout)}


def judge(approval_path, wt_args: list, python=None) -> int:
    r = check(approval_path, python)
    if r["status"] != "OK":
        print("gate_env: 判定専用の環境が承認と違う — 判定しない:\n- " + "\n- ".join(r["problems"]), file=sys.stderr)
        return 1
    from ..intake import approval as AP
    apr = AP.load_approval(approval_path)
    co = Path(r["checkout"])
    env = sanitized_env(co, apr["gate"]["env"].get("tools") or {})
    args = [str(x) for x in wt_args]
    if "--approval" in args:
        raise GateEnvError("wt_qualify の引数に --approval を書かない (gate_env が足す)")
    cmd = [r["python"], "-I", "-B", "-c", LAUNCH_PY, str(co), *args, "--approval", str(Path(approval_path).resolve())]
    return subprocess.run(cmd, env=env, cwd=str(co)).returncode


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m forge_design.feedback.gate_env", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot", help="判定専用の環境の記録 (JSON)")
    s.add_argument("--python", default=None)
    c = sub.add_parser("check", help="承認の記録の gate と今の環境の照合")
    c.add_argument("--approval", required=True)
    pr = sub.add_parser("probe", help="判定専用の環境で import の元と印の検査だけを回す")
    pr.add_argument("--approval", required=True)
    j = sub.add_parser("judge", help="判定専用の環境で wt_qualify --approval を回す")
    j.add_argument("--approval", required=True)
    j.add_argument("wt_args", nargs=argparse.REMAINDER, help="-- の後に wt_qualify の引数")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "snapshot":
            print(json.dumps(snapshot(a.python), indent=1, ensure_ascii=False))
            return 0
        if a.cmd == "check":
            r = check(a.approval)
            print(json.dumps(r, indent=1, ensure_ascii=False))
            print(f"VERDICT: {r['status']}")
            return 0 if r["status"] == "OK" else 1
        if a.cmd == "probe":
            r = probe(a.approval)
            print(json.dumps(r, indent=1, ensure_ascii=False))
            return 0 if r["status"] == "OK" else 1
        wa = a.wt_args[1:] if a.wt_args[:1] == ["--"] else a.wt_args
        return judge(a.approval, wa)
    except GateEnvError as e:
        print(f"gate_env: 止めた — {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
