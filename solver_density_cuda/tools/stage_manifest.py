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
    # **`turbulence: {model: ...}` の書式も拾う** (2026-09-20 codex result M3)。
    # `turbulenceModel:` だけを見ていたため、層流 (`model: "none"`) と SST (`model: "sst"`)
    # が**同一キー**になり、段階起動の層流段と SST 段が 1 区間に連結されていた。
    # 区間分離という本ツールの目的そのものが効いていなかった。
    # **正規表現は flow 形式・二重引用符しか拾えない**ので、下の `YAML_PATHS` による
    # 構造解析が正本で、この行は YAML が壊れている場合の保険 (2026-09-21 codex result M4)。
    ("turbulenceModel", r"turbulenceModel\s*:\s*(\S+)"),
    ("turbulence.model", r"turbulence\s*:\s*\{[^}]*?\bmodel\s*:\s*[\"']?([A-Za-z0-9_]+)"),
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


# **YAML として構造解析して拾う hard キー** (2026-09-21 codex result M4)。
# 正規表現は書式に依存する: `turbulence:\n  model: "sst"` (block 形式) や
# `turbulence: {model: 'sst'}` (単引用符) では層流と SST が同一キーになっていた。
# **構造解析を正本とし、正規表現は YAML が読めないときの保険**にする。
YAML_HARD_PATHS = [
    ("turbulence.model", ("turbulence", "model")),
    ("turbulence.wallTreatmentSST", ("turbulence", "wallTreatmentSST")),
    ("turbulence.transition", ("turbulence", "transition")),   # 遷移モデルの有無は別の方程式系
    ("space.convMethod", ("space", "convMethod")),
    ("space.limiter", ("space", "limiter")),
    ("space.limiterScaled", ("space", "limiterScaled")),
    ("space.venkatK", ("space", "venkatK")),
    ("physics.viscMethod", ("physics", "viscMethod")),
    ("physics.thermalMethod", ("physics", "thermalMethod")),
    ("time.unsteady", ("time", "unsteady")),
    ("time.dualTime", ("time", "dualTime")),
    # **ソルバ内 CHT** (2026-09-23)。連成の有無・固体モデル・界面熱量の定義・更新間隔・平均窓・D_f は
    # いずれも**解いている方程式を変える**ので、これらが違う段を 1 区間として連結してはいけない
    # (plan boundary-conjugate-heat-transfer §4.6a / §5.1 #67 ④)。
    # 固体そのもの (メッシュ・孔 Robin・k_s(T)) は下の `conjugate.solid_sha1` で見る。
    ("conjugate.mode", ("conjugate", "mode")),
    ("conjugate.flux", ("conjugate", "flux")),
    ("conjugate.interval", ("conjugate", "interval")),
    ("conjugate.flux_avg", ("conjugate", "flux_avg")),
    ("conjugate.Df_scale", ("conjugate", "Df_scale")),
    ("conjugate.thickness", ("conjugate", "thickness")),
    ("conjugate.k_solid", ("conjugate", "k_solid")),
    ("conjugate.T_b", ("conjugate", "T_b")),
    # 背面条件は固体抵抗 $R_{\rm tot}=t/k_s+R_{\rm back}$ を変える = 別の方程式
    # (codex result 2026-09-23 M7: `h_c` 100 と 10000 が同一キーになっていた)。
    ("conjugate.back", ("conjugate", "back")),
    ("conjugate.h_c", ("conjugate", "h_c")),
    ("conjugate.relax", ("conjugate", "relax")),
    # 壁隣接面の質量流束 χ を面法線マッハから取る opt-in (plan convection-slau-wall-normal-chi, 既定 0)。
    # 空間離散化を変えるので 0 と 1 を同一区間にしない (codex result-3 M5, 2026-09-25)。
    ("space.slauWallNormalChi", ("space", "slauWallNormalChi")),
    # 流束方式そのもの。SLAU→ROE の切替を 1 区間に連結しない (plan convection-slau-wall-normal-chi-usage-rule
    # codex plan m8, 2026-09-25)。solver は必須キーで既定値が無いので、大文字小文字だけ揃える。
    ("solver", ("solver",)),
]

# 既定値と同じなら**キーごと落とす** hard キー。省略と明示の既定値を同一区間にし、
# このキーを足す前に書かれた stage_manifest.json (キー無し) とも一致させるため。
# (`space.limiterScaled` / `space.venkatK` も省略と明示で割れるが、既定が変わった経緯があり別判断。)
YAML_HARD_DEFAULTS = {"space.slauWallNormalChi": "0"}


def _yaml_grab(text, paths):
    """solverConfig を YAML として読み、指定パスの値を拾う。読めなければ空 dict。"""
    try:
        import yaml
    except ImportError:
        return {}
    try:
        doc = yaml.safe_load(text or "")
    except Exception:
        return {}
    if not isinstance(doc, dict):
        return {}
    out = {}
    for name, path in paths:
        cur = doc
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                cur = None
                break
            cur = cur[k]
        if cur is not None and not isinstance(cur, (dict, list)):
            out[name] = str(cur).strip().strip('"\'')
    return out


def _grab(text, pats):
    out = {}
    for name, pat in pats:
        m = re.search(pat, text or "")
        if m:
            # flow 形式 `{convMethod: 1, limiter: 2}` では `\S+` が区切りの `,` `}` まで拾い、
            # 同じ値が末尾のキーか否かで `2}` / `2,` に割れていた (block 形式とも不一致)。
            # 区切り記号を落として書式に依存しない値にする (2026-09-25, plan convection-slau-wall-normal-chi #12)。
            out[name] = m.group(1).rstrip(",}]").strip().strip('"\'')
    return out


def _conjugate_solid_sha1(cfg_text, run_dir):
    """`conjugate.solid` が指す固体 HDF5 の sha1 (先頭 12 桁)。

    固体メッシュ・孔 Robin・$k_s(T)$ は**解いている方程式そのもの**なので、
    これが変わった段を同一区間として連結してはいけない。ファイルが見つからないときは
    `missing` を返す (黙って「同じ」にしない)。
    """
    try:
        import yaml
        cj = (yaml.safe_load(cfg_text) or {}).get("conjugate")
    except Exception:
        return None
    if not isinstance(cj, dict) or "solid" not in cj:
        return None
    if not run_dir:
        return "unknown"
    path = os.path.join(str(run_dir), str(cj["solid"]))
    if not os.path.exists(path):
        return "missing"
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def stage_key(cfg_text, bcond_text, run_dir=None):
    """段の **hard キー** (方程式・BC・空間離散化)。BC は全文のハッシュで見る。"""
    k = _grab(cfg_text, HARD_PATTERNS)
    # **構造解析の値で上書きする** (正規表現より優先。書式差で取りこぼさないため)。
    k.update(_yaml_grab(cfg_text, YAML_HARD_PATHS))
    if "solver" in k:
        k["solver"] = k["solver"].upper()
    for name, dflt in YAML_HARD_DEFAULTS.items():
        if name in k and k[name] == dflt:
            del k[name]
    k["bcond_sha1"] = hashlib.sha1((bcond_text or "").encode()).hexdigest()[:12]
    # ソルバ内 CHT: 連成の有無と固体の中身 (plan §5.1 #67 ④)。
    try:
        import yaml
        k["conjugate.enabled"] = "1" if isinstance((yaml.safe_load(cfg_text) or {}).get("conjugate"), dict) else "0"
    except Exception:
        k["conjugate.enabled"] = "1" if re.search(r"^conjugate\s*:", cfg_text or "", re.M) else "0"
    sha = _conjugate_solid_sha1(cfg_text, run_dir)
    if sha is not None:
        k["conjugate.solid_sha1"] = sha
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
            "key": stage_key(cfg_text, bcond_text, self.run),
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
