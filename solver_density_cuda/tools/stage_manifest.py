#!/usr/bin/env python3
r"""段階起動の **段ごとの実効設定**を記録・照合する (収束判定の区間を名前でなくデータで決めるため)。

なぜ必要か (2026-09-19 codex plan レビュー Major 1):

段階起動 (slip → 層流 → 等温 → SST → 2 次ランプ → 本段) の run では、
**本段だけを判定すると前段の収束場から始まるので低下桁数が小さく出る**一方、
**全段を連結すると別の方程式・BC の過渡を本段の基準にしてしまう**
(`check_convergence.py` は系列全体の最大値を低下桁数の基準に取る)。
段名のプレフィックスで選ぶ方式は**方程式・BC・離散化の同一性を保証しない**し、
`residual_history_all.csv` のような集約済みファイルを二重に拾う危険もある。

そこで、run 生成側が段ごとに

    {"order": n, "tag": ..., "history": "residual_history_<tag>.csv",
     "restart_from": <前段の tag or null>,
     "key": {"equations": ..., "bc": ..., "space": ...},   # ← **これが同じ区間だけ連結してよい**
     "soft": {"cfl": ..., "ninner": ...}}                   # ← 違っても連結してよい

を `stage_manifest.json` に書き、判定側はこれを読む。

使い方 (run 生成側):

    from stage_manifest import StageManifest
    sm = StageManifest(run_dir)
    sm.add("S0_slip", cfg_text, bcond_text, history="residual_history_S0_slip.csv")
    ...
    sm.write()

使い方 (判定側):

    python3 solver_density_cuda/tools/stage_manifest.py <run_dir> [--segments]

`--segments` は**連結してよい最大の区間**を順に出す。最後の区間が本段を含む判定区間になる。
"""
import argparse
import hashlib
import json
import os
import re
import sys

# 連結の可否を決める **hard キー** (これが変わったら別区間)
HARD_PATTERNS = [
    # (表示名, solverConfig の正規表現)
    ("convMethod", r"convMethod\s*:\s*(\S+)"),
    ("limiter", r"\blimiter\s*:\s*(\S+)"),
    # リミッタの**式そのもの**を変えるキー。これが変わった段を同一区間として連結してはいけない
    # (codex 2026-09-20 limiter-config-simplify plan レビュー M5)。
    ("limiterScaled", r"limiterScaled\s*:\s*(\S+)"),
    ("venkatK", r"venkatK\s*:\s*(\S+)"),
    ("turbulenceModel", r"turbulenceModel\s*:\s*(\S+)"),
    ("wallTreatmentSST", r"wallTreatmentSST\s*:\s*(\S+)"),
    ("viscMethod", r"viscMethod\s*:\s*(\S+)"),
    ("thermalMethod", r"thermalMethod\s*:\s*(\S+)"),
    ("unsteady", r"unsteady\s*:\s*(\S+)"),
    ("dualTime", r"dualTime\s*:\s*(\S+)"),
    ("axisymmetric", r"axisymmetric\s*:\s*(\S+)"),
]
# 変わっても連結してよい **soft キー**
SOFT_PATTERNS = [
    ("cfl", r"\bcfl\s*:\s*(\S+)"),
    ("cfl_pseudo", r"cfl_pseudo\s*:\s*(\S+)"),
    ("nStepInner", r"nStepInner\s*:\s*(\S+)"),
    ("nStepOuter", r"nStepOuter\s*:\s*(\S+)"),
    ("implicitRelax", r"implicitRelax\s*:\s*(\S+)"),
]


def _grab(text, pats):
    out = {}
    for name, pat in pats:
        m = re.search(pat, text or "")
        if m:
            out[name] = m.group(1).strip().strip('"\'')
    return out


def stage_key(cfg_text, bcond_text):
    """段の **hard キー** (方程式・BC・空間離散化)。BC は全文のハッシュで見る。"""
    k = _grab(cfg_text, HARD_PATTERNS)
    k["bcond_sha1"] = hashlib.sha1((bcond_text or "").encode()).hexdigest()[:12]
    return k


class StageManifest:
    def __init__(self, run_dir):
        self.run = str(run_dir)
        self.stages = []

    def add(self, tag, cfg_text, bcond_text, history=None, restart_from=None):
        self.stages.append({
            "order": len(self.stages),
            "tag": tag,
            "history": history or ("residual_history_%s.csv" % tag),
            "restart_from": restart_from if restart_from is not None else (
                self.stages[-1]["tag"] if self.stages else None),
            "key": stage_key(cfg_text, bcond_text),
            "soft": _grab(cfg_text, SOFT_PATTERNS),
        })

    def write(self):
        p = os.path.join(self.run, "stage_manifest.json")
        with open(p, "w") as f:
            json.dump({"stages": self.stages}, f, indent=2, ensure_ascii=False)
        return p


def segments(man):
    """hard キーが同じ連続区間に切る。戻り値 [[stage, ...], ...]"""
    segs, cur = [], []
    for st in man["stages"]:
        if cur and st["key"] != cur[-1]["key"]:
            segs.append(cur); cur = []
        cur.append(st)
    if cur:
        segs.append(cur)
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--segments", action="store_true", help="連結してよい区間を出す")
    a = ap.parse_args()
    p = os.path.join(a.run, "stage_manifest.json")
    if not os.path.exists(p):
        print("stage_manifest.json が無い: %s" % p)
        print("  -> 段階起動の run なら生成側で StageManifest を使うこと。"
              "無い run では**判定区間を人が明示する**しかない (AGENTS.md 収束確認)。")
        return 2
    man = json.load(open(p))
    segs = segments(man)
    print("=== %s : %d 段 / %d 区間 ===" % (a.run, len(man["stages"]), len(segs)))
    for i, sg in enumerate(segs):
        tags = " -> ".join(s["tag"] for s in sg)
        print("  区間 %d: %s" % (i, tags))
        diff = {k: v for k, v in sg[0]["key"].items()}
        print("        hard: %s" % json.dumps(diff, ensure_ascii=False))
        softs = {s["tag"]: s["soft"] for s in sg}
        if len(set(json.dumps(v, sort_keys=True) for v in softs.values())) > 1:
            print("        soft (区間内で変化・連結可): %s" % json.dumps(softs, ensure_ascii=False))
    print("\n**判定区間は最後の区間** (%s)。"
          % " -> ".join(s["tag"] for s in segs[-1]) if segs else "(段が無い)")
    print("  `check_convergence.py` にはこの区間の履歴だけを渡すこと。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
