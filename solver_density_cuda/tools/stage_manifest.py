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
    # node のスカラー勾配の作用素 (gg / lsq)。明示 gg↔lsq を 1 区間に連結しない (plan gradient-scalar-lsq-unification
    # §5.1 #5j、codex result-1 M5)。**暫定**: 省略 (= 既定) と明示 gg は別キーになる (分けすぎる側)。起動順・バイナリ id・
    # legacy の区別は別 plan `tooling-stage-manifest-launch-binding` (#2g)。
    ("mesh.scalarGradient", ("mesh", "scalarGradient")),
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
    # space.slauWallNormalChi は既定が構成依存 (auto) になったので YAML の値でなく**実効値**で扱う
    # (下の infer_wall_normal_chi と segments()。plan convection-slau-wall-normal-chi-default §4.3)。
    # 流束方式そのもの。SLAU→ROE の切替を 1 区間に連結しない (plan convection-slau-wall-normal-chi-usage-rule
    # codex plan m8, 2026-09-25)。solver は必須キーで既定値が無いので、大文字小文字だけ揃える。
    ("solver", ("solver",)),
]

CHI_KEY = "space.slauWallNormalChi.effective"
MANIFEST_VERSION = 2

# 厚さ 0 の板の自由端の近傍で速度の再構成を節点値にする処置 (plan convection-zero-thickness-edge-reconstruction §4
# 「設計チェーンへの配線」、codex plan レビュー 2026-10-08 M3、codex diagnose 2026-10-08 edge-weight-preprocessing)。
# ソルバの設定は有効/無効とフィールド名だけ (`space.zeroThicknessEdgeVelocity: {field: w_recon_vel}`)。端の判定
# (タグ・rings) は前処理の道具 (`tools/mark_zero_thickness_edges.py`) が **meshFileName の `/AUX/<field>`** に書いた
# 節点の重み w とその属性にある。区間識別には有効状態・フィールド名・属性のタグ (並べ替え)・rings に加えて、
# **その段の forge が読む meshFileName から再計算した格子署名と w のハッシュ**を入れる (属性に書かれたハッシュは
# 転記しない。計算は道具の `mesh_signature(h5path)` / `field_hash(w)` を import して使い、同じ実装を二重に持たない)。
# **キーが無い・明示の無効の段には何も足さない** (既存 run の stage_key は不変)。
ZTE_CFG_KEY = "zeroThicknessEdgeVelocity"
ZTE_PREFIX = "space." + ZTE_CFG_KEY
ZTE_GROUP = "AUX"
ZTE_TOOL_MODULE = "mark_zero_thickness_edges"


ZTE_DEFAULT_FIELD = "w_recon_vel"


def zte_config(cfg_text):
    """solverConfig の `space.zeroThicknessEdgeVelocity` を読む (ソルバ solverConfig.cpp と同じ解釈)。
    戻り: None (キーが無い・`{enabled: 0}` = 無効) か {"field": 名前 or None, "invalid": 理由 or None}。
    有効は `{enabled: 1}` か `{field: ...}` (field の既定は w_recon_vel)。map でない値 (null・false・数・文字列)・未知キー・
    enabled が 0/1/真偽値でないものはソルバが起動時に拒否するので、invalid の有効扱い (無効の段と黙って連結しない)。
    YAML が読めないときも、キーの文字列が本文にあれば invalid の有効扱い。"""
    try:
        import yaml
    except ImportError:
        raise SystemExit("stage_manifest: PyYAML が無いので段の区間を判定できない (pip install pyyaml)")
    try:
        doc = yaml.safe_load(cfg_text or "")
    except Exception:
        return {"field": None, "invalid": "yaml_unparsable"} if ZTE_CFG_KEY in (cfg_text or "") else None
    sp = doc.get("space") if isinstance(doc, dict) else None
    if not isinstance(sp, dict) or ZTE_CFG_KEY not in sp:
        return None
    v = sp[ZTE_CFG_KEY]
    if not isinstance(v, dict):
        return {"field": None, "invalid": "not_a_map:%r" % (v,)}
    unknown = sorted(str(k) for k in v if k not in ("enabled", "field"))
    if unknown:
        return {"field": None, "invalid": "unknown_key:%s" % ",".join(unknown)}
    en = v.get("enabled", 1)
    if isinstance(en, bool):
        en = int(en)
    if en not in (0, 1) or isinstance(en, float):
        return {"field": None, "invalid": "enabled:%r" % (en,)}
    f = v.get("field", ZTE_DEFAULT_FIELD)
    if not (isinstance(f, str) and re.match(r"^[A-Za-z0-9_.-]+$", f.strip())):
        return {"field": None, "invalid": "field_not_a_name:%r" % (f,)}
    if en == 0:
        return None
    return {"field": f.strip(), "invalid": None}


def zte_signature(tags, rings):
    """処置の識別子 (設計 DB・学習の採否で使う)。タグは並べ替えて重複を落とす。"""
    return "tags=%s;rings=%d" % (",".join(sorted(set(str(t) for t in tags))), int(rings))


def _zte_json(v):
    """h5 の属性を JSON にできる値へ (来歴の記録用)。"""
    try:
        import numpy as np
    except ImportError:
        np = None
    if isinstance(v, bytes):
        return v.decode(errors="replace")
    if np is not None and isinstance(v, np.ndarray):
        return [_zte_json(x) for x in v.tolist()] if v.size <= 64 else "array(%s, %d)" % (v.dtype, v.size)
    if np is not None and isinstance(v, np.generic):
        return v.item()
    if isinstance(v, (list, tuple)):
        return [_zte_json(x) for x in v]
    return v


def _zte_tags(v):
    """属性のタグ (文字列 "a,b" / "a b"、または文字列の配列) を並べ替えたリストに。読めなければ None。"""
    v = _zte_json(v)
    if isinstance(v, str):
        items = [s for s in re.split(r"[,\s]+", v) if s]
    elif isinstance(v, list) and all(isinstance(s, str) for s in v):
        items = [s.strip() for s in v if s.strip()]
    else:
        return None
    return sorted(set(items)) or None


def zte_tool():
    """前処理の道具 `mark_zero_thickness_edges.py` (格子署名 `mesh_signature(h5path)` と w のハッシュ `field_hash(w)` の正本)。
    import できない・関数が無いときは止める (有効の段を識別できないまま他の段と連結しない)。"""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        import importlib
        t = importlib.import_module(ZTE_TOOL_MODULE)
    except ImportError as e:
        raise SystemExit("stage_manifest: 処置の前処理の道具 %s.py を import できない (%s)。格子署名と w のハッシュを"
                         "再計算できないので、処置を有効にした段を識別できない" % (ZTE_TOOL_MODULE, e))
    for fn in ("mesh_signature", "field_hash"):
        if not callable(getattr(t, fn, None)):
            raise SystemExit("stage_manifest: %s.py に %s() が無い" % (ZTE_TOOL_MODULE, fn))
    return t


def zte_mesh_file(cfg_text, run_dir):
    """処置の w を読む入力 h5 = **mesh.meshFileName** (valueFileName ではない)。決められなければ None。"""
    if not run_dir:
        return None
    try:
        import yaml
        ms = (yaml.safe_load(cfg_text or "") or {}).get("mesh") or {}
    except Exception:
        return None
    v = ms.get("meshFileName") if isinstance(ms, dict) else None
    return os.path.join(str(run_dir), str(v)) if v else None


def zte_identity(h5path, field):
    """meshFileName の `/AUX/<field>` を読み、識別に要る量を返す。ファイル・フィールドが無ければ None。
    戻り: {"path", "n" (長さ), "tags" (属性、並べ替え; 読めなければ None), "rings" (属性; int or None),
           "mesh_signature" (道具で再計算), "field_hash" (道具で再計算), "attrs" (全属性、来歴用)}。
    **ハッシュは属性から転記せず毎回再計算する** (同じ設定でも格子や w が変われば別物として扱うため)。"""
    if not h5path or not os.path.exists(str(h5path)):
        return None
    try:
        import h5py
        import numpy as np
    except ImportError:
        raise SystemExit("stage_manifest: h5py / numpy が無いので処置のフィールドを照合できない")
    path = "%s/%s" % (ZTE_GROUP, field)
    with h5py.File(str(h5path), "r") as f:
        if path not in f or not isinstance(f[path], h5py.Dataset):
            return None
        ds = f[path]
        attrs = {k: _zte_json(v) for k, v in ds.parent.attrs.items()}
        attrs.update({k: _zte_json(v) for k, v in ds.attrs.items()})
        w = np.asarray(ds[...])
    rings = attrs.get("rings")
    try:
        rings = int(rings) if rings is not None and float(rings) == int(rings) else None
    except (TypeError, ValueError):
        rings = None
    t = zte_tool()
    return {"path": path, "n": int(w.shape[0]) if w.shape else 0,
            "tags": _zte_tags(attrs.get("tags")) if "tags" in attrs else None, "rings": rings,
            "mesh_signature": str(t.mesh_signature(str(h5path))), "field_hash": str(t.field_hash(w)), "attrs": attrs}


def _zte_stage_key(cfg_text, zc, run_dir):
    """有効な段の区間識別 (zte_config が None でないときだけ呼ぶ)。h5 が読めないときは unknown / missing を値にして、
    他の段と黙って連結しない。"""
    out = {ZTE_PREFIX + ".enabled": "1",
           ZTE_PREFIX + ".field": zc["field"] if zc["field"] else "invalid:" + str(zc["invalid"])}
    names = ("tags", "rings", "mesh_signature", "field_hash")
    if not zc["field"]:
        val = dict.fromkeys(names, "invalid")
    elif not run_dir:
        val = dict.fromkeys(names, "unknown")
    else:
        a = zte_identity(zte_mesh_file(cfg_text, run_dir), zc["field"])
        if a is None:
            val = dict.fromkeys(names, "missing")
        else:
            val = {"tags": ",".join(a["tags"]) if a["tags"] else "missing",
                   "rings": str(a["rings"]) if a["rings"] is not None else "missing",
                   "mesh_signature": a["mesh_signature"], "field_hash": a["field_hash"]}
    for k, v in val.items():
        out[ZTE_PREFIX + "." + k] = v
    return out


def fnv1a64(data):
    """forge の appendLaunchRecord (main.cpp) と同じ FNV-1a 64。段と forge_launches.jsonl の起動を結び付ける。"""
    h = 14695981039346656037
    for b in (data.encode() if isinstance(data, str) else data):
        h ^= b
        h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return format(h, "x")


def infer_wall_normal_chi(cfg_text):
    """solverConfig から slauWallNormalChi の実効値を**推定**する (C++ solverConfig.cpp の auto 解決と同じ規則)。
    戻り (値 "0"/"1", 由来 "explicit"/"inferred")。起動の記録 (forge_launches.jsonl) があればそちらが正本。"""
    try:
        import yaml
        doc = yaml.safe_load(cfg_text or "") or {}
    except Exception:
        return "0", "inferred"
    sp = doc.get("space") or {}
    if isinstance(sp, dict) and "slauWallNormalChi" in sp:
        return str(int(sp["slauWallNormalChi"])), "explicit"
    ms = doc.get("mesh") or {}
    node = str(ms.get("discretization", "")).strip('"\'') == "node"
    nwd = int(ms.get("nodeWallDirichlet", 1)) == 1
    slau = str(doc.get("solver", "")).strip('"\'').upper() in ("SLAU", "SLAU2")
    return ("1" if (node and nwd and slau) else "0"), "inferred"


AMBIGUOUS = "ambiguous"


def load_launches(run_dir):
    """forge_launches.jsonl を読む (cfg_fnv → 実効値)。無ければ空。
    **同じ cfg_fnv の起動で実効値が食い違う**ときは値を `ambiguous` にする (旧バイナリ → 新バイナリで同じ YAML を
    起動し直した等。どの段がどの起動か分からないので、最後の起動の値で全段を上書きしない。codex result 2026-09-27
    chi-default M3。起動との結び付けは plan tooling-stage-manifest-launch-binding で後回し)。"""
    out = {}
    p = os.path.join(str(run_dir), "forge_launches.jsonl") if run_dir else None
    if p and os.path.exists(p):
        for line in open(p):
            try:
                r = json.loads(line)
                v = str(int(r["slauWallNormalChi"]))
            except Exception:
                continue
            k = r.get("cfg_fnv")
            out[k] = v if out.get(k, v) == v else AMBIGUOUS
    return out


def _yaml_grab(text, paths):
    """solverConfig を YAML として読み、指定パスの値を拾う。読めなければ空 dict。"""
    try:
        import yaml
    except ImportError:
        # PyYAML が無いと YAML 側だけにある hard キー (mesh.scalarGradient 等) が消えて段が黙って連結される
        # (codex 2026-09-26 gradient-scalar-lsq result-2 M1)。区間判定の必須依存として止める。
        raise SystemExit("stage_manifest: PyYAML が無いので段の区間を判定できない (pip install pyyaml)")
    try:
        doc = yaml.safe_load(text or "")
    except Exception:
        # 読めない YAML を空 dict にすると、hard キーが全部欠けた同士として**黙って連結**される
        # (codex 2026-09-26 gradient-scalar-lsq result-1 M5)。本文のハッシュを入れて別区間にする。
        return {"yaml_unparsable": hashlib.sha1((text or "").encode()).hexdigest()[:12]}
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


def _twophase_effective(cfg_text, bcond_text):
    """二相拡散の実効値 (0/1; tools/twophase_state.py)。凝縮 OFF かつ未指定、または YAML が読めなければ None (キーを足さない)。

    `condTwoPhaseDiffusion` は液・モーメントの輸送作用素を変えるので方程式の一部 (plan condensation-two-phase-default §4-4)。
    OFF の場から ON で継続した段を同一区間に連結しないため、実効値を hard キーに入れる。指定値でなく実効値を見るのは、
    不活性 (inactive-a/b) な指定の違いでは作用素が変わらないから。"""
    try:
        import yaml
        import importlib.util
        y = yaml.safe_load(cfg_text or "") or {}
        b = yaml.safe_load(bcond_text or "") if bcond_text else None
        if not isinstance(y, dict):
            return None
        cd = y.get("condensation") or {}
        if not isinstance(cd, dict) or (int(cd.get("condensation", 0) or 0) != 1 and "condTwoPhaseDiffusion" not in cd):
            return None
        here = os.path.dirname(os.path.abspath(__file__))
        spec = importlib.util.spec_from_file_location("twophase_state", os.path.join(here, "twophase_state.py"))
        tps = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(tps)
        return tps.classify(y, b if isinstance(b, dict) else None)["effective"]
    except Exception:
        return None


def stage_key(cfg_text, bcond_text, run_dir=None):
    """段の **hard キー** (方程式・BC・空間離散化)。BC は全文のハッシュで見る。"""
    k = _grab(cfg_text, HARD_PATTERNS)
    # **構造解析の値で上書きする** (正規表現より優先。書式差で取りこぼさないため)。
    k.update(_yaml_grab(cfg_text, YAML_HARD_PATHS))
    if "solver" in k:
        k["solver"] = k["solver"].upper()
    k[CHI_KEY] = infer_wall_normal_chi(cfg_text)[0]
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
    tp = _twophase_effective(cfg_text, bcond_text)
    if tp is not None:
        k["twophase_diffusion_effective"] = str(tp)   # 凝縮 run だけ (他の run の署名は従来どおり)
    # 厚さ 0 の板の自由端の速度再構成の処置: 有効な段だけ (無効・キー無しの段の署名は従来どおり)
    zc = zte_config(cfg_text)
    if zc is not None:
        k.update(_zte_stage_key(cfg_text, zc, run_dir))
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
            "cfg_fnv": fnv1a64(cfg_text or ""),
            "chi_source": infer_wall_normal_chi(cfg_text)[1],
            "soft": _grab(cfg_text, SOFT_PATTERNS),
        })

    def write(self):
        p = os.path.join(self.run, "stage_manifest.json")
        with open(p, "w") as f:
            json.dump({"manifest_version": MANIFEST_VERSION, "stages": self.stages}, f, indent=2, ensure_ascii=False)
        return p


def _normalize(man, launches):
    """段の key を実効値で揃える (旧形式の移行と起動記録での確定)。戻り: 段のコピーのリスト。"""
    out = []
    for st in man["stages"]:
        st = dict(st); k = dict(st["key"])
        if CHI_KEY not in k:
            # 旧形式 (manifest_version < 2): 明示 1 は "space.slauWallNormalChi": "1"、キー無しは当時の既定 0
            k[CHI_KEY] = str(k.pop("space.slauWallNormalChi", "0"))
            st["chi_source"] = "legacy"
        fnv = st.get("cfg_fnv")
        if launches and fnv in launches:
            k[CHI_KEY] = launches[fnv]
            st["chi_source"] = AMBIGUOUS if launches[fnv] == AMBIGUOUS else "launch"
        st["key"] = k
        out.append(st)
    return out


def segments(man, launches=None):
    """hard キーが同じ連続区間に切る。戻り [[stage, ...], ...]。
    実効値が**推定** (inferred) の段と、確定 (launch/explicit/legacy) の段は、値が同じでも連結しない
    (由来不明の推定値を既知と自動連結しない。plan convection-slau-wall-normal-chi-default §4.3、codex plan M1)。"""
    segs, cur = [], []
    for st in _normalize(man, launches):
        inf = st.get("chi_source") == "inferred"
        amb = st.get("chi_source") == AMBIGUOUS or (cur and cur[-1].get("chi_source") == AMBIGUOUS)
        if cur and (amb or st["key"] != cur[-1]["key"] or inf != (cur[-1].get("chi_source") == "inferred")):
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
    segs = segments(man, load_launches(a.run))
    print("=== %s : %d 段 / %d 区間 ===" % (a.run, len(man["stages"]), len(segs)))
    for i, sg in enumerate(segs):
        tags = " -> ".join(s["tag"] for s in sg)
        print("  区間 %d: %s" % (i, tags))
        diff = {k: v for k, v in sg[0]["key"].items()}
        print("        hard: %s" % json.dumps(diff, ensure_ascii=False))
        softs = {s["tag"]: s["soft"] for s in sg}
        if len(set(json.dumps(v, sort_keys=True) for v in softs.values())) > 1:
            print("        soft (区間内で変化・連結可): %s" % json.dumps(softs, ensure_ascii=False))
    if segs and any(st.get("chi_source") == AMBIGUOUS for st in segs[-1]):
        print("\n**判定不能**: 最後の区間の段は同じ設定で chi の実効値が違う起動が混ざっている (%s)。"
              "どの起動の履歴か分からないので、判定区間を人が明示すること。" % " -> ".join(s["tag"] for s in segs[-1]))
        return 3
    print("\n**判定区間は最後の区間** (%s)。"
          % " -> ".join(s["tag"] for s in segs[-1]) if segs else "(段が無い)")
    print("  `check_convergence.py` にはこの区間の履歴だけを渡すこと。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
