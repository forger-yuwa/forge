#!/usr/bin/env python3
"""case/66 の入力テンプレート inputs/<名前>/ を元 run から複製して作る (ローカルで実行)。

    python3 prepare_inputs.py [名前 ...]        # 省略時は matrix_spec.INPUTS の全部
    python3 prepare_inputs.py --list

- 元 run (matrix_spec.SRC_ROOT 配下、別セッションのワークツリー) は**読むだけ**。書き込みはこのディレクトリの inputs/ だけ。
- 設定の修正は matrix_spec の edits (正規表現・期待一致数) だけを行い、置換後に PyYAML で読んで expect の値を検査する。
  コメントと書式は元のまま残る (yaml を読み直して書き出さない: 重複キー・anchor の罠を避けるため)。
- `mesh.bndFirstOrder` と `wallTreatmentSST: 1` が修正後に残っていたら失敗させる (AGENTS.md・recommended-settings §2)。
- 同一メッシュ restart の種は seed_src.h5 として置くだけで、restart_field.py は AWS で掛ける (run_matrix.py seed)。
  化学種の記録の照合に forge (--resolve-species) が要るため。
- 各 inputs/<名前>/SOURCE.txt に元・複製したファイル・修正・sha256 を残す。FILES に run へ写すファイルの一覧を残す。
"""
import argparse
import hashlib
import os
import re
import shutil
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matrix_spec as ms  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def get_key(d, dotted):
    cur = d
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return "<無し>"
        cur = cur[k]
    return cur


def check_no_dup_keys(text, fname):
    """重複キー・anchor/alias を含む yaml は自動で扱わない (solver の yaml-cpp は先勝ち、PyYAML は後勝ち)。"""
    if re.search(r"(^|\s)[&*][A-Za-z0-9_]+", text) or "<<:" in text:
        raise SystemExit(f"{fname}: anchor/alias/merge key を含むので自動で扱わない")

    class Loader(yaml.SafeLoader):
        pass

    def construct_mapping(loader, node, deep=False):
        keys = set()
        for k, _ in node.value:
            kk = loader.construct_object(k, deep=deep)
            if kk in keys:
                raise SystemExit(f"{fname}: 重複キー {kk!r}")
            keys.add(kk)
        return yaml.SafeLoader.construct_mapping(loader, node, deep)

    Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)
    return yaml.load(text, Loader=Loader)


def prepare(name, spec):
    dst = os.path.join(HERE, "inputs", name)
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    src = os.path.join(ms.SRC_ROOT, spec["src"])
    log = [f"input: {name}", f"src  : {src}", ""]
    files = []

    # 派生 (変換器): 既に作った入力から config を持ってくる
    if "derive" in spec:
        dsrc = os.path.join(HERE, "inputs", spec["derive"]["input"])
        for fn in spec["derive"]["files"]:
            p = os.path.join(dsrc, fn)
            if not os.path.exists(p):
                raise SystemExit(f"{name}: 派生元 {p} が無い (先に {spec['derive']['input']} を作る)")
            shutil.copy2(p, os.path.join(dst, fn))
            files.append(fn)
            log.append(f"derive {fn} <- inputs/{spec['derive']['input']}/{fn}")

    for item in spec.get("copy", []):
        frm, to = (item["from"], item["to"]) if isinstance(item, dict) else (item, item)
        p = os.path.normpath(os.path.join(src, frm))
        if not os.path.exists(p):
            raise SystemExit(f"{name}: 元ファイル {p} が無い")
        shutil.copy2(p, os.path.join(dst, to))
        files.append(to)
        log.append(f"copy  {to} <- {p} (sha256 {sha256(p)[:16]})")

    if "seed" in spec:
        p = os.path.normpath(os.path.join(src, spec["seed"]["src"]))
        shutil.copy2(p, os.path.join(dst, "seed_src.h5"))
        log.append(f"seed  seed_src.h5 <- {p} (sha256 {sha256(p)[:16]}) → AWS で restart_field.py を {spec['seed']['dst']} へ")
        # 種の隣の解決済み記録 (化学種) があれば一緒に持っていく (restart_field の照合が読む)
        for fn in sorted(os.listdir(os.path.dirname(p))):
            if fn.startswith("resolved_species_") and fn.endswith(".yaml"):
                shutil.copy2(os.path.join(os.path.dirname(p), fn), os.path.join(dst, fn))
                log.append(f"copy  {fn} (種の隣の解決済み記録)")

    for fn, content in spec.get("files", {}).items():
        with open(os.path.join(dst, fn), "w") as f:
            f.write(content)
        if fn not in files:
            files.append(fn)
        log.append(f"write {fn} (matrix_spec で生成)")

    for fn, edits in spec.get("edits", {}).items():
        p = os.path.join(dst, fn)
        text = open(p).read()
        for pat, rep, cnt in edits:
            new, n = re.subn(pat, rep, text)
            if n != cnt:
                raise SystemExit(f"{name}/{fn}: 置換 {pat!r} の一致が {n} 件 (期待 {cnt})")
            text = new
            log.append(f"edit  {fn}: {pat!r} -> {rep!r} ({n} 件)")
        with open(p, "w") as f:
            f.write(text)

    if os.path.exists(os.path.join(dst, "solverConfig.yaml")):
        text = open(os.path.join(dst, "solverConfig.yaml")).read()
        cfg = check_no_dup_keys(text, f"{name}/solverConfig.yaml")
        if re.search(r"bndFirstOrder", text):
            raise SystemExit(f"{name}: mesh.bndFirstOrder が残っている (使用禁止)")
        if get_key(cfg, "turbulence.wallTreatmentSST") == 1:
            raise SystemExit(f"{name}: wallTreatmentSST 1 が残っている (使用禁止)")
        for k, v in spec.get("expect", {}).items():
            got = get_key(cfg, k)
            if got != v:
                raise SystemExit(f"{name}: {k} = {got!r} (期待 {v!r})")
        log.append("")
        log.append(f"nStepOuter = {get_key(cfg, 'time.last.nStepOuter')}, "
                   f"outStepInterval = {get_key(cfg, 'time.outStepInterval')}")
    for fn in ("bcondConfig.yaml",):
        if os.path.exists(os.path.join(dst, fn)):
            check_no_dup_keys(open(os.path.join(dst, fn)).read(), f"{name}/{fn}")

    if "ckpt" in spec:
        c = spec["ckpt"]
        files.append(c["dst"])
        log.append(f"ckpt  {c['dst']} は AWS で置く: base r1 の {c['cfg']} の {c['file']} (run_matrix.py set-ckpt)")
    if "seed" in spec and spec["seed"]["dst"] not in files:
        files.append(spec["seed"]["dst"])
    log.append(f"note : {spec.get('note', '')}")
    with open(os.path.join(dst, "SOURCE.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    with open(os.path.join(dst, "FILES"), "w") as f:
        f.write("\n".join(files) + "\n")
    print(f"[prepare] {name}: {len(files)} files")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for k, v in ms.INPUTS.items():
            print(f"{k:22s} {v['src']}")
        return
    names = a.names or list(ms.INPUTS)
    # 派生する入力は派生元より後に作る
    names = [n for n in names if "derive" not in ms.INPUTS[n]] + [n for n in names if "derive" in ms.INPUTS[n]]
    for n in names:
        prepare(n, ms.INPUTS[n])


if __name__ == "__main__":
    main()
