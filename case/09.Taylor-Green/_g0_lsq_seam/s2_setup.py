#!/usr/bin/env python3
"""S2 物理 A/B の双子 run を作る (plan gradient-scalar-lsq-unification §6 S2・§5.1 #5)。

起点 run の入力ファイルを新しい run ディレクトリへ複製し、solverConfig を編集して (`mesh.scalarGradient` を明示、
step 数・出力間隔、廃止キーの削除、任意の上書き)、起点の res を `restart_field.py` で保存量 index コピーする。
最後に `check_solver_config.py` を通す。gg と lsq の双子は同じ引数で `--grad` だけ変えて作る。

  python3 s2_setup.py --start DIR --res DIR/res_N.h5 --dst RUN_DIR --grad gg|lsq --nstep N --out-int K \\
      [--extra species_db.yaml inlet_profile_1.csv ...] [--drop mesh.nodeAxisDirichlet ...] [--set time.deltaT.cfl=2.0 ...]

- `--dst` が既にあれば失敗する (上書きしない)。
- YAML に anchor / alias / 重複キーがあれば失敗する (一括書き換えで節を取り違えないため)。
- `IC_FROM.txt` に起点 res のパスと sha256、変えた設定の差分を書く。
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
TOOLS = os.path.join(REPO, "solver_density_cuda", "tools")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class UniqueLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader, node, deep=False):
    keys = set()
    for k, _ in node.value:
        kk = loader.construct_object(k, deep=deep)
        if kk in keys:
            raise SystemExit(f"YAML に重複キー {kk!r} (一括書き換えを拒否)")
        keys.add(kk)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def parse_val(s):
    return yaml.safe_load(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--res", default=None, help="起点の res (保存量を index コピー)。--no-restart なら不要")
    ap.add_argument("--no-restart", action="store_true", help="起点の入力 h5 の VALUE をそのまま IC に使う (起点 run と同じ IC から回す双子)")
    ap.add_argument("--dst", required=True)
    ap.add_argument("--grad", choices=("gg", "lsq"), required=True)
    ap.add_argument("--nstep", type=int, required=True)
    ap.add_argument("--out-int", type=int, required=True)
    ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--drop", nargs="*", default=[])
    ap.add_argument("--set", nargs="*", default=[])
    a = ap.parse_args()

    if os.path.exists(a.dst):
        raise SystemExit(f"{a.dst} は既にある (上書きしない)")
    cfg_txt = open(os.path.join(a.start, "solverConfig.yaml")).read()
    if re.search(r"(^|[\s:\[{,])[&*][A-Za-z_]", cfg_txt) or "<<:" in cfg_txt:
        raise SystemExit("solverConfig.yaml に anchor/alias/merge key がある (一括書き換えを拒否)")
    cfg = yaml.load(cfg_txt, Loader=UniqueLoader)
    mesh_name = cfg["mesh"]["meshFileName"]
    value_name = cfg["mesh"].get("valueFileName", mesh_name)   # IC の書き込み先 (restart_field の DST)

    changes = []
    for key in a.drop:
        sec, k = key.split(".", 1) if "." in key else (None, key)
        node = cfg if sec is None else cfg.get(sec, {})
        # 2 段 (time.deltaT.x) まで
        if "." in k:
            s2, k = k.split(".", 1); node = node.get(s2, {})
        if k in node:
            changes.append(f"drop {key} (起点値 {node[k]!r})"); del node[k]
    for kv in a.set:
        key, v = kv.split("=", 1)
        parts = key.split(".")
        node = cfg
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        old = node.get(parts[-1], "<なし>")
        node[parts[-1]] = parse_val(v)
        changes.append(f"set {key}: {old!r} -> {node[parts[-1]]!r}")
    cfg["mesh"]["scalarGradient"] = a.grad
    changes.append(f"set mesh.scalarGradient: {a.grad} (明示)")
    old_n = cfg["time"]["last"].get("nStepOuter")
    cfg["time"]["last"]["nStepOuter"] = a.nstep
    changes.append(f"set time.last.nStepOuter: {old_n} -> {a.nstep}")
    old_o = cfg["time"].get("outStepInterval")
    cfg["time"]["outStepInterval"] = a.out_int
    cfg["time"]["outStepStart"] = 0
    changes.append(f"set time.outStepInterval: {old_o} -> {a.out_int}")

    os.makedirs(a.dst)
    for f in ["bcondConfig.yaml", "probe.yaml", mesh_name] + ([value_name] if value_name != mesh_name else []) + a.extra:
        src = os.path.join(a.start, f)
        if not os.path.exists(src):
            if f == "probe.yaml":
                continue
            raise SystemExit(f"起点に {f} が無い")
        shutil.copy2(src, a.dst)
    with open(os.path.join(a.dst, "solverConfig.yaml"), "w") as fp:
        fp.write(f"# S2 双子 (plan gradient-scalar-lsq-unification §6)。起点 {a.start}\n")
        yaml.safe_dump(cfg, fp, sort_keys=False, default_flow_style=None, allow_unicode=True)
    if a.no_restart:
        changes.append(f"IC: 起点の {value_name} の VALUE をそのまま使う (restart なし)")
    else:
        if not a.res:
            raise SystemExit("--res か --no-restart が要る")
        r = subprocess.run([sys.executable, os.path.join(TOOLS, "restart_field.py"), a.res, os.path.join(a.dst, value_name)],
                           capture_output=True, text=True)
        print(r.stdout[-1500:], r.stderr[-1500:])
        if r.returncode != 0:
            raise SystemExit("restart_field.py 失敗")
    r2 = subprocess.run([sys.executable, os.path.join(TOOLS, "check_solver_config.py"), os.path.join(a.dst, "solverConfig.yaml")],
                        capture_output=True, text=True)
    print(r2.stdout[-2000:], r2.stderr[-800:])
    with open(os.path.join(a.dst, "IC_FROM.txt"), "w") as fp:
        fp.write(f"start_run: {a.start}\nres: {a.res}\nres_sha256: {sha(a.res) if a.res else '-'}\n"
                 f"mesh_sha256 (起点の {mesh_name}): {sha(os.path.join(a.start, mesh_name))}\n"
                 f"value file (restart 先): {value_name}\n")
        fp.write("config changes vs start:\n" + "".join(f"  - {c}\n" for c in changes))
        fp.write(f"check_solver_config rc={r2.returncode}\n")
    print("\n".join(changes))
    if r2.returncode != 0:
        raise SystemExit("check_solver_config が失敗 (上の出力)")


if __name__ == "__main__":
    main()
