#!/usr/bin/env python3
"""凍結入力に対する EOS 1 回の全出力のビット比較を、準備済みの入力ディレクトリ 1 つについて回す (AWS の g3 用、ローカルでも可)。

  prepare: 既存 run から入力を作る (res_*・ログ・判定・時系列を除いて複製 → restart_field で同じ res から保存量を移す →
           nStepOuter 1・outStepInterval 1・output.floorEvents 1 にする。両版で同じファイルを使う: 旧版は floorEvents を読まない)
      python3 eos_ab_dir.py prepare --src-run RUN --res res_20000.h5 --input IN --tools TOOLS --forge NEWBIN
  run: 旧版 2 回・新版 2 回 (入力ディレクトリを毎回複製し、入力ファイルの sha256 が複製元と同じことを確かめてから回す)、
       5 組 (旧同士・新同士・旧対新 3 組) を compare_eos_dump.py で比べ、ダンプは比較が済んだ順に消す (--keep-dumps で残す)。同時に残るダンプは最大 3 本。
      python3 eos_ab_dir.py run --input IN --work W --old OLDBIN --new NEWBIN --tools TOOLS [--env FORGE_CUDA_BLOCKSIZE=128]
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

SKIP_PREFIX = ("res_",)
SKIP_NAMES = {"forge_run.log", "residual_history.csv", "residual_history.png", "floor_events.csv", "RUN_PROVENANCE.txt",
              "CONVERGENCE_VERDICT.txt", "forge_launches.jsonl", "eos_dump.h5"}
HDF5_LIB = "/usr/lib/x86_64-linux-gnu/hdf5/serial"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def value_file(d):
    c = open(os.path.join(d, "solverConfig.yaml")).read()
    return re.search(r'valueFileName:\s*"?([^",}\s]+)', c).group(1)


def set_output_on(cfg):
    m = re.search(r"(?m)^output:\s*\{([^}]*)\}", cfg)
    if m:
        inner = m.group(1).strip()
        if "floorEvents" in inner:
            return re.sub(r"floorEvents:\s*\d+", "floorEvents: 1", cfg)
        return cfg[:m.start()] + "output: {" + (inner + ", " if inner else "") + "floorEvents: 1}" + cfg[m.end():]
    if re.search(r"(?m)^output:", cfg):
        raise SystemExit("block 形式の output 節は未対応 (手で output.floorEvents: 1 を入れて --no-edit で回す)")
    return cfg.rstrip("\n") + "\noutput: {floorEvents: 1}\n"


def prepare(a):
    if os.path.exists(a.input):
        raise SystemExit(f"{a.input} は既にある (上書きしない)")
    os.makedirs(a.input)
    copied = []
    for f in sorted(os.listdir(a.src_run)):
        p = os.path.join(a.src_run, f)
        if os.path.isdir(p) or f.startswith(SKIP_PREFIX) or f.endswith(".xmf") or f in SKIP_NAMES or f.endswith(".msh"):
            continue
        shutil.copy2(p, a.input)
        copied.append(f)
    val = value_file(a.input)
    if val not in copied:
        raise SystemExit(f"入力 h5 {val} が {a.src_run} に無い")
    src = os.path.join(a.src_run, a.res)
    env = dict(os.environ, FORGE_BIN=a.forge, LD_LIBRARY_PATH=HDF5_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
    r = subprocess.run([sys.executable, os.path.join(a.tools, "restart_field.py"), src, os.path.join(a.input, val)],
                       env=env, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    print(out)
    if r.returncode != 0 or "VERDICT: OK" not in out:
        raise SystemExit("restart_field が OK でない")
    c = open(os.path.join(a.input, "solverConfig.yaml")).read()
    c2 = re.sub(r"nStepOuter:\s*\d+", "nStepOuter: 1", c)
    c2 = re.sub(r"(?m)^(\s*)outStepInterval:\s*\d+", r"\g<1>outStepInterval: 1", c2)
    c2 = set_output_on(c2)
    open(os.path.join(a.input, "solverConfig.yaml"), "w").write(c2)
    try:
        import yaml
        y = yaml.safe_load(c2)
        print("config: nStepOuter", y["time"]["last"]["nStepOuter"], "output", y.get("output"))
    except Exception as e:  # noqa: BLE001
        print("config の YAML 検査をしなかった:", e)
    files = {f: sha(os.path.join(a.input, f)) for f in sorted(os.listdir(a.input)) if os.path.isfile(os.path.join(a.input, f))}
    json.dump({"src_run": a.src_run, "res": a.res, "restart_field": out.splitlines()[-1] if out else "", "files": files},
              open(os.path.join(a.input, "INPUT.json"), "w"), indent=1)
    print(f"入力 {a.input}: {len(files)} ファイル (INPUT.json に sha256)")


def run_one(a, inp, files, name, binp, extra_env):
    rd = os.path.join(a.work, name)
    if os.path.exists(rd):
        shutil.rmtree(rd)
    shutil.copytree(inp, rd)
    bad = [f for f, h in files.items() if f != "INPUT.json" and sha(os.path.join(rd, f)) != h]
    if bad:
        raise SystemExit(f"{name}: 複製した入力が複製元と違う {bad}")
    dump = os.path.join(rd, "eos_dump.h5")
    env = dict(os.environ, FORGE_DUMP_EOS_STEP=str(a.step), FORGE_DUMP_EOS_FILE=dump,
               LD_LIBRARY_PATH=HDF5_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
    env.update(extra_env)
    with open(os.path.join(rd, "forge_run.log"), "w") as lf:
        r = subprocess.run([binp], cwd=rd, env=env, stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    log = open(os.path.join(rd, "forge_run.log")).read()
    ok = r.returncode == 0 and "[eos-dump] wrote" in log and os.path.exists(dump)
    fe = os.path.join(rd, "floor_events.csv")
    rows = [l for l in open(fe).read().splitlines() if l.startswith(("init,", "eos,"))] if os.path.exists(fe) else []
    # 入力の大きいファイルと初期出力は run から消す (正本は入力ディレクトリ)
    for f in list(files) + [x for x in os.listdir(rd) if x.startswith("res_")]:
        p = os.path.join(rd, f)
        if f not in ("solverConfig.yaml", "INPUT.json") and os.path.isfile(p) and os.path.getsize(p) > (1 << 20):
            os.remove(p)
    info = {"rc": r.returncode, "ok": ok, "dump_sha256": sha(dump) if os.path.exists(dump) else None,
            "dump_bytes": os.path.getsize(dump) if os.path.exists(dump) else 0, "floor_rows": rows}
    print(f"  {name}: rc={r.returncode} dump={'あり' if os.path.exists(dump) else '無し'} ({info['dump_bytes'] / 2**30:.2f} GiB) "
          f"floor_events {len(rows)} 行", flush=True)
    for l in rows:
        print("     ", l[:150])
    if not ok:
        raise SystemExit(f"{name}: forge が EOS ダンプを書いて終わらなかった (forge_run.log を見る)")
    return dump, info


def compare(a, x, y, dx, dy, res):
    out = os.path.join(a.work, f"cmp_{x}_vs_{y}")
    r = subprocess.run([sys.executable, os.path.join(a.tools, "compare_eos_dump.py"), dx, dy, "--json", out + ".json"],
                       capture_output=True, text=True)
    open(out + ".txt", "w").write(r.stdout + r.stderr)
    j = json.load(open(out + ".json"))
    t = j.get("totals", {})
    res[f"{x}_vs_{y}"] = {"verdict": j["verdict"], "eos_args_post_diff_bytes": t.get("eos_args_post_diff_bytes"),
                          "eos_args_pre_diff_bytes": t.get("eos_args_pre_diff_bytes"), "eos_args_bytes": t.get("eos_args_bytes"),
                          "other_pre_diff_bytes": t.get("other_pre_diff_bytes"), "other_post_diff_bytes": t.get("other_post_diff_bytes"),
                          "invalid": j.get("invalid"), "different": [d["array"] for d in j.get("different", [])]}
    print(f"  {x} vs {y}: VERDICT {j['verdict']}  EOS 引数 post 差分バイト {t.get('eos_args_post_diff_bytes')} / {t.get('eos_args_bytes')}"
          f"  pre 差分 {t.get('eos_args_pre_diff_bytes')}  その他 pre/post 差分 {t.get('other_pre_diff_bytes')}/{t.get('other_post_diff_bytes')}", flush=True)
    for s in j.get("invalid", [])[:4]:
        print("      不成立:", s[:160])


def run(a):
    os.makedirs(a.work, exist_ok=True)
    inp = os.path.abspath(a.input)
    files = json.load(open(os.path.join(inp, "INPUT.json")))["files"] if os.path.exists(os.path.join(inp, "INPUT.json")) else \
        {f: sha(os.path.join(inp, f)) for f in sorted(os.listdir(inp)) if os.path.isfile(os.path.join(inp, f))}
    extra = dict(kv.split("=", 1) for kv in (a.env or []))
    S = {"input": inp, "old": {"bin": a.old, "sha256": sha(a.old)}, "new": {"bin": a.new, "sha256": sha(a.new)},
         "env": extra, "step": a.step, "runs": {}, "compare": {}}
    print(f"old {S['old']['sha256'][:16]}… new {S['new']['sha256'][:16]}… env {extra}")
    d = {}
    d["old_r1"], S["runs"]["old_r1"] = run_one(a, inp, files, "old_r1", a.old, extra)
    d["new_r1"], S["runs"]["new_r1"] = run_one(a, inp, files, "new_r1", a.new, extra)
    compare(a, "old_r1", "new_r1", d["old_r1"], d["new_r1"], S["compare"])
    d["old_r2"], S["runs"]["old_r2"] = run_one(a, inp, files, "old_r2", a.old, extra)
    compare(a, "old_r1", "old_r2", d["old_r1"], d["old_r2"], S["compare"])
    compare(a, "old_r2", "new_r1", d["old_r2"], d["new_r1"], S["compare"])
    if not a.keep_dumps:
        os.remove(d["old_r2"])
    d["new_r2"], S["runs"]["new_r2"] = run_one(a, inp, files, "new_r2", a.new, extra)
    compare(a, "new_r1", "new_r2", d["new_r1"], d["new_r2"], S["compare"])
    compare(a, "old_r1", "new_r2", d["old_r1"], d["new_r2"], S["compare"])
    if not a.keep_dumps:
        for k in ("old_r1", "new_r1", "new_r2"):
            os.remove(d[k])
    json.dump(S, open(os.path.join(a.work, "SUMMARY.json"), "w"), indent=1)
    v = {k: c["verdict"] for k, c in S["compare"].items()}
    print("SUMMARY:", v)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--src-run", required=True)
    p.add_argument("--res", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--tools", required=True, help="restart_field.py のある tools ディレクトリ")
    p.add_argument("--forge", required=True, help="restart_field の化学種解決に使う forge (新版)")
    r = sub.add_parser("run")
    r.add_argument("--input", required=True)
    r.add_argument("--work", required=True)
    r.add_argument("--old", required=True)
    r.add_argument("--new", required=True)
    r.add_argument("--tools", required=True, help="compare_eos_dump.py のある tools ディレクトリ")
    r.add_argument("--step", type=int, default=1)
    r.add_argument("--env", action="append", help="KEY=VALUE (両版に同じ値を渡す。例 FORGE_CUDA_BLOCKSIZE=128)")
    r.add_argument("--keep-dumps", action="store_true")
    a = ap.parse_args()
    return prepare(a) if a.cmd == "prepare" else run(a)


if __name__ == "__main__":
    sys.exit(main() or 0)
