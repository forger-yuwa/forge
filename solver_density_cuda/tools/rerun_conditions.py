#!/usr/bin/env python3
r"""既存 run の形状 (格子) を固定したまま、入口条件・背圧・入口分率だけを変えた run を作る (run の準備まで; forge は起動しない)。

    python3 solver_density_cuda/tools/rerun_conditions.py REF_RUN NEW_RUN [--res res_N.h5]
        [--Pt P] [--Tt T] [--Y NAME=v ... | --Y1 v] [--balance NAME] [--k K] [--omega W]
        [--Ps P | --keep-Ps] [--Tw T | --keep-Tw] [--steps N] [--out-interval M] [--cfl C]
        [--scale-ic none|pt] [--override-recommended] [--forge BIN] [--dry-run]

plan: plans/active/tooling-rerun-conditions.md §4 (設計方針)。手順の正本は procedures/nozzle-design-workflow.md §3a、
仕様の解説は methods/design/overview.md「既存 run の条件変更 (rerun_conditions)」。

やること (順に; どこで止まっても NEW_RUN は残さない):
1. **入力契約の検査** (NEW_RUN を作る前): 単一の `inlet_Pressure` (floats に `Y{s}`、`inletProfile` なし)・
   単一の `outlet_statPress`・`wall`/`wall_isothermal`/`slip` (Euler の滑り壁)・`axis` だけ、bcond は 1 行 1 境界の flow 形式、
   `meshFileName == valueFileName == "nozzle.h5"`、solverConfig のファイル参照は許可リストの run 内相対だけ。
2. **書き換えの計画**: 指定値と参照の実効値の有限性・物理範囲 (Pt・Tt・Ps・Tw > 0、k ≥ 0、omega > 0、cfl > 0、
   steps・out-interval > 0) を検査。bcond は対象行の floats の当該トークンだけを置換し (他の行はバイト一致)、YAML で読み直して検証。
   Y は全種を書いて Σ=1 を 1e−12 で検査 (forge の起動検査は 1e−3 で、入口カーネルが黙って正規化する)。変更しないときも
   参照の入口 Y を各成分 [0, 1]・Σ=1 (1e−9) で検査する。
   Pt を変えたら `--Ps` か `--keep-Ps` が必須 (node の出口は壁列が亜音速で Ps を見る)。
   `recommended_stages` (plan §4.7) と生成 config の本段 cfl・step 数が食い違えば、必要な引数を示して作成前に停止
   (`--override-recommended` で明示的に通す)。
3. **初期場の検査**: config から必要保存量集合を決め、SRC (参照 res) と DST (参照 nozzle.h5) の存在・形・有限・ρ>0・
   0 ≤ ρY/ρ ≤ 1+1e−6・|Σ ρY − ρ| ≤ 1e−6 ρ を検査。DST の /VALUE に集合と `wall_dist` 以外があれば拒否。
4. NEW_RUN を作り、許可リストのファイルだけ複製 → 書き換え → `restart_field.py` (VERDICT OK 必須) →
   `--scale-ic pt` のときだけ必要保存量を全部 f = Pt_new/Pt_ref 倍 (T・U・Y・k・ω を保つ初期場変換; 変換後も場を検査) → 記録。
5. `RERUN_CONDITIONS.json` (変更前後・必要保存量・スケール検査・P_exit_ref・restart の VERDICT・forge の sha256・
   ツールの commit・recommended_stages) と `prepare_info.json` (`ic_from`・`rerun_of`) を書く。

v1 で止めるもの: 乾き成分 lump の組成変更 (`--lump`; restart_field・convert_species_field・forge の全経路が拒否する)、
`X{s}` 形式・`inletProfile: 1`・複数 inlet・外部参照、凝縮 block・Tt 変更・組成変更と `--scale-ic pt` の併用。
"""
import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import forge_species as fsp  # noqa: E402

REPO_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MESH_FILE = "nozzle.h5"
# 複製の許可リスト (系譜・ログ・VERDICT・series・report は持ち込まない; plan §4.2)
ALLOW_COPY = ("nozzle.h5", "nozzle.xmf", "bcondConfig.yaml", "solverConfig.yaml", "species_meta.yaml",
              "resolved_species_*.yaml", "species_db_external.yaml", "probe.yaml", "prepare_info.json",
              "wall_*.csv", "target_axis_M.csv", "delta_r_initial.*", "MESH_QUALITY.txt")
REQUIRED_FILES = ("nozzle.h5", "bcondConfig.yaml", "solverConfig.yaml")
INERT_WHEN_NO_TURB = ("roK", "roOmega")   # 乱流モデルなしの run に残る変換器既定の入れ物
ALLOWED_KINDS = ("inlet_Pressure", "outlet_statPress", "wall", "wall_isothermal", "slip", "axis")   # slip = Euler の滑り壁 (2026-10-06 追加)
Y_SUM_TOL = 1e-12          # 書いた Y の Σ=1 の許容 (forge の起動検査 1e-3 より厳しく; plan §4.3)
REF_Y_SUM_TOL = 1e-9       # 変更しない参照の入口 Y の Σ=1 の許容 (参照 BC の ΣY = 1.0001 などを拒否; codex result 段 #1)
FIELD_Y_TOL = 1e-6         # 場の 0 ≤ ρY/ρ ≤ 1+tol, |Σ ρY − ρ| ≤ tol·ρ (plan §4.4)
SCALE_RTOL = 1e-6          # スケール検査 allclose(d_new, f·d_ref, rtol, atol=0)
KEEP_FROM_DST = ("wall_dist",)
CONDITION_KEYS = ("Pt", "Tt", "Y", "k", "omega", "Ps", "Tw")


class RerunError(Exception):
    """入力契約・条件の不整合・初期場の検査で止めるとき (NEW_RUN は作らない / 消す)。"""


# ---------------------------------------------------------------------------
# 必要保存量集合と場の検査 (runner_axismach.run_staged_ns の段終了ゲートも使う)
# ---------------------------------------------------------------------------
def required_conserved_from_cfg(cfg):
    """solverConfig (dict) から restart に必須の保存量データセット名を返す (plan §4.4(1))。
    = forge_species.required_conserved (ro/roU/roe/roY*/roXi) + SST → roK・roOmega + 遷移 → roGamma・roReth
      + 凝縮 → rog_*・roQ2_*・roQ1_*・roQ0_* (variables::registerCondensation の命名)。"""
    pp = cfg.get("physProp") or {}
    tm = int(pp.get("thermalMethod", 0))
    names = [str(s["name"]) if isinstance(s, dict) else str(s) for s in (pp.get("species") or [])]
    if tm == 2 and not names:
        names = ["N2"]                      # solverConfig.cpp の既定 (単成分 N2)
    tracer = str(pp.get("tracer") or "").strip() or None
    if tracer in ("none", "0"):
        tracer = None
    req = list(fsp.required_conserved({"thermalMethod": tm, "names": names, "tracer": tracer}))
    tb = cfg.get("turbulence") or {}
    model = str(tb.get("model") or "none").strip().lower()
    if model.startswith("sst"):
        req += ["roK", "roOmega"]
    if str(tb.get("transition") or "none").strip().lower() not in ("none", "0"):
        req += ["roGamma", "roReth"]
    cd = cfg.get("condensation") or {}
    if isinstance(cd, dict) and int(cd.get("condensation", 0) or 0) == 1:
        for s in range(int(cd.get("nCondSpecies", 0) or 0)):
            req += [f"{b}_{s}" for b in ("rog", "roQ2", "roQ1", "roQ0")]
    return req


def field_problems(h5path, required, species_bounds=True):
    """h5 の /VALUE について必要保存量の存在・1 次元で同じ長さ・有限・ρ>0 (と species_bounds なら ρY の範囲・総和) を検査する。
    戻り (問題のリスト, 節点数)。問題が空なら合格。"""
    import h5py
    probs, n = [], None
    with h5py.File(h5path, "r") as f:
        if "VALUE" not in f:
            return [f"{h5path}: /VALUE が無い"], None
        V = f["VALUE"]
        data = {}
        for name in required:
            if name not in V:
                probs.append(f"{name}: 無い")
                continue
            a = np.asarray(V[name])
            if a.ndim != 1:
                probs.append(f"{name}: 1 次元でない (shape {a.shape})")
                continue
            if n is None:
                n = a.shape[0]
            elif a.shape[0] != n:
                probs.append(f"{name}: 長さ {a.shape[0]} != {n}")
                continue
            bad = int(np.count_nonzero(~np.isfinite(a)))
            if bad:
                probs.append(f"{name}: 非有限 {bad}/{a.size} 点")
                continue
            data[name] = a.astype(np.float64)
    if "ro" in data:
        nb = int(np.count_nonzero(data["ro"] <= 0.0))
        if nb:
            probs.append(f"ro: ρ ≤ 0 が {nb} 点")
    if species_bounds and "ro" in data and not probs:
        ro = data["ro"]
        ys = [k for k in required if re.fullmatch(r"roY\d+", k)]
        if ys:
            tot = np.zeros_like(ro)
            for k in ys:
                y = data[k] / ro
                if y.min() < 0.0 or y.max() > 1.0 + FIELD_Y_TOL:
                    probs.append(f"{k}/ro: 範囲外 [{y.min():.3e}, {y.max():.9f}] (0 ≤ Y ≤ 1+{FIELD_Y_TOL:g})")
                tot += data[k]
            err = float(np.max(np.abs(tot - ro) / ro))
            if err > FIELD_Y_TOL:
                probs.append(f"|Σ ρY − ρ|/ρ 最大 {err:.3e} > {FIELD_Y_TOL:g}")
    return probs, n


# ---------------------------------------------------------------------------
# bcond (1 行 1 境界の flow 形式) の解析とトークン置換
# ---------------------------------------------------------------------------
_LINE_RX = re.compile(r"^([A-Za-z_][\w.-]*)\s*:\s*(\{.*\})\s*(#.*)?$")
_FLOATS_RX = re.compile(r"\bfloats\s*:\s*(\{[^{}]*\})")


def parse_bcond(text):
    """bcond を 1 行 1 境界の flow 形式として読む。戻り {名前: {line, kind, floats, ints, entry}} (出現順)。
    flow 形式でない行・YAML として読めない・名前が重複するものは RerunError。"""
    try:
        doc = yaml.safe_load(text)
    except yaml.YAMLError as e:
        raise RerunError(f"bcondConfig.yaml が YAML として読めない: {e}")
    if not isinstance(doc, dict):
        raise RerunError("bcondConfig.yaml が境界名 → 設定の辞書でない")
    out = {}
    for i, line in enumerate(text.splitlines()):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        m = _LINE_RX.match(line)
        if not m:
            raise RerunError(f"bcondConfig.yaml {i + 1} 行目が 1 行 1 境界の flow 形式でない (v1 は flow 形式だけ): {line!r}")
        name = m.group(1)
        if name in out:
            raise RerunError(f"bcondConfig.yaml: 境界名 {name} が重複")
        ent = doc.get(name)
        if not isinstance(ent, dict):
            raise RerunError(f"bcondConfig.yaml: {name} の設定が辞書でない")
        out[name] = {"line": i, "kind": str(ent.get("kind")), "floats": dict(ent.get("floats") or {}),
                     "ints": dict(ent.get("ints") or {}), "entry": ent}
    if set(out) != set(doc):
        raise RerunError(f"bcondConfig.yaml の行と YAML のキーが一致しない ({sorted(out)} vs {sorted(doc)})")
    return out


def fmt(v):
    """書き込む数値の表記 (往復で同じ float に戻る最短表記)。"""
    return repr(float(v))


def replace_float_tokens(line, updates):
    """1 行の floats {…} の中で updates のキーの値だけを置換する。各キーはちょうど 1 回当たること (足さない)。"""
    m = _FLOATS_RX.search(line)
    if not m:
        raise RerunError(f"floats {{…}} が無い行: {line!r}")
    body = m.group(1)
    for key, val in updates.items():
        rx = re.compile(r"(?<![\w.])(" + re.escape(key) + r"\s*:\s*)([^,}\s]+)")
        hits = rx.findall(body)
        if len(hits) != 1:
            raise RerunError(f"floats のキー {key} が {len(hits)} 回当たる (ちょうど 1 回のはず): {line!r}")
        body = rx.sub(lambda mm: mm.group(1) + val, body, count=1)
    return line[:m.start(1)] + body + line[m.end(1):]


def _sub_once(rx, repl, text, what):
    hits = re.findall(rx, text)
    if len(hits) != 1:
        raise RerunError(f"solverConfig.yaml の {what} がちょうど 1 回当たらない ({len(hits)} 回)")
    return re.sub(rx, repl, text, count=1)


# ---------------------------------------------------------------------------
# 参照 run の読み取り
# ---------------------------------------------------------------------------
def _res_step(path):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(path))
    return int(m.group(1)) if m else None


def find_ref_res(ref, res_arg):
    if res_arg:
        p = res_arg if os.path.isabs(res_arg) or os.path.exists(res_arg) else os.path.join(ref, res_arg)
        if not os.path.exists(p):
            raise RerunError(f"--res {res_arg} が無い")
        return os.path.abspath(p)
    cands = [p for p in glob.glob(os.path.join(ref, "res_[0-9]*.h5")) if _res_step(p) is not None]
    if not cands:
        raise RerunError(f"{ref} に res_[0-9]*.h5 が無い")
    return os.path.abspath(max(cands, key=_res_step))


def exit_pressure_ref(ref, res_path, outlet_id, excluded_ids=()):
    """参照の出口圧 P_exit_ref = 参照 res の**出口断面 (最終 x の節点列) の内部節点の静圧 P の中央値**。
    壁節点 (wall / wall_isothermal の BC 節点) と軸節点 (axis の BC 節点) は除く (`excluded_ids` = その physID)。
    `res_outlet_*` の Ps は課した値そのものなので使わない (2026-10-06 主セッション決定; plan §4.1)。
    出口 BC の節点が単一の x に並ばない (節点列が取れない) ときは、出口 BC の節点と同じ x 座標を持つ節点で代替する。
    節点の座標と BC 節点は参照 run の nozzle.h5 (MESH/COORD・BCONDS/<physID>/iCells; node 離散化では節点番号) から読む。
    戻り (値 or None, 出所の説明)。取れなければ None (呼び手は停止せず警告として記録する)。"""
    import h5py
    mesh = os.path.join(ref, MESH_FILE)
    with h5py.File(res_path, "r") as f:
        if "VALUE/P" not in f:
            return None, f"取れない ({os.path.basename(res_path)} に VALUE/P が無い)"
        P = np.asarray(f["VALUE/P"], dtype=np.float64)
    with h5py.File(mesh, "r") as f:
        if "MESH/COORD" not in f:
            return None, f"取れない ({MESH_FILE} に MESH/COORD が無い)"
        C = np.asarray(f["MESH/COORD"], dtype=np.float64)

        def _nodes(pid):
            k = f"BCONDS/{pid}/iCells"
            return np.asarray(f[k], dtype=np.int64) if k in f else None
        out_nodes = _nodes(outlet_id)
        excl = [n for n in (_nodes(pid) for pid in excluded_ids) if n is not None]
    if C.size != 3 * P.size:
        return None, "取れない (COORD と P の点数が合わない; node 離散化でない)"
    x = C.reshape(-1, 3)[:, 0]
    tol = 1e-6 * max(float(x.max() - x.min()), 1e-300)
    if out_nodes is None or out_nodes.size == 0 or out_nodes.max() >= x.size:
        return None, f"取れない ({MESH_FILE} に出口 BCONDS/{outlet_id}/iCells が無い)"
    xo = x[out_nodes]
    if float(xo.max() - xo.min()) <= tol:
        sel = np.abs(x - float(np.mean(xo))) <= tol          # 最終 x の節点列
        how = f"出口断面 x={float(np.mean(xo)):.6g} の節点列"
    else:
        xs = np.unique(xo)                                    # 代替: 出口 BC の節点と同じ x を持つ節点
        k = np.clip(np.searchsorted(xs, x), 1, xs.size - 1)
        sel = np.minimum(np.abs(x - xs[k - 1]), np.abs(x - xs[k])) <= tol
        how = "出口 BC の節点と同じ x の節点 (出口が単一 x に並ばない)"
    n_all = int(sel.sum())
    for n in excl:
        sel[n[n < x.size]] = False
    if not sel.any():
        return None, f"取れない ({how} に壁・軸以外の節点が無い)"
    return (float(np.median(P[sel])),
            f"{os.path.basename(res_path)}: {how} の内部節点 {int(sel.sum())}/{n_all} 点 (壁・軸の BC 節点を除く) の P median")


def _file_refs(obj, path=""):
    """solverConfig 内のファイル参照らしき文字列値 (拡張子で判定) を列挙する。"""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            out += _file_refs(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += _file_refs(v, f"{path}[{i}]")
    elif isinstance(obj, str) and re.search(r"\.(ya?ml|h5|hdf5|csv|dat|txt|json|xmf|msh|inp)$", obj.strip(), re.I):
        out.append((path, obj.strip()))
    return out


def _allowed_name(name):
    import fnmatch
    return any(fnmatch.fnmatch(name, pat) for pat in ALLOW_COPY)


def _repo_rel(p):
    p = os.path.abspath(p)
    try:
        r = os.path.relpath(p, REPO_ROOT)
    except ValueError:
        return p
    return p if r.startswith("..") else r


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _num(v):
    """YAML の値を数値として読む (`5e+0` は PyYAML では文字列になるが forge は数値として読む)。数値でなければ None。"""
    if isinstance(v, bool) or v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _check_range(what, v, *, zero_ok=False):
    """v が有限で > 0 (zero_ok なら ≥ 0) であることを検査する。違反は RerunError。"""
    x = _num(v)
    if x is None or not np.isfinite(x) or (x < 0.0 if zero_ok else x <= 0.0):
        raise RerunError(f"{what} = {v!r} は有限で {'≥ 0' if zero_ok else '> 0'} でなければならない")
    return x


def _deltaT(cfg):
    return ((cfg.get("time") or {}).get("deltaT") or {})


# ---------------------------------------------------------------------------
# 計画 (検査と書き換え内容の決定; 何も書かない)
# ---------------------------------------------------------------------------
def build_plan(a):
    """引数から書き換え計画を作る。入力契約・条件の不整合・初期場の検査はすべてここで行い、違反は RerunError。"""
    ref = os.path.abspath(a.ref_run)
    new = os.path.abspath(a.new_run)
    if a.lump is not None:
        raise RerunError("--lump (乾き成分 lump の組成変更) は v1 で停止する: restart_field (互換ハッシュ不一致)・"
                         "convert_species_field conserve (実種ごとの ρY 保存検査)・reinit (ξ の出所なし)・"
                         "forge (旧ハッシュ属性) の全経路が拒否する (plan §4.4(4)、§5.1 #9)")
    # --- 指定値の有限性と物理範囲 (codex result 段 #1: 非物理な値で計画が通り、--Pt 0 で ZeroDivisionError だった) ---
    for what, v, zok in (("--Pt", a.Pt, False), ("--Tt", a.Tt, False), ("--Ps", a.Ps, False), ("--Tw", a.Tw, False),
                         ("--k", a.k, True), ("--omega", a.omega, False), ("--cfl", a.cfl, False),
                         ("--steps", a.steps, False), ("--out-interval", a.out_interval, False)):
        if v is not None:
            _check_range(what, v, zero_ok=zok)
    if os.path.exists(new):
        raise RerunError(f"NEW_RUN {new} が既にある (上書きしない)")
    if not os.path.isdir(ref):
        raise RerunError(f"REF_RUN {ref} が無い")
    for fn in REQUIRED_FILES:
        if not os.path.exists(os.path.join(ref, fn)):
            raise RerunError(f"REF_RUN に {fn} が無い")

    cfg_text = open(os.path.join(ref, "solverConfig.yaml"), encoding="utf-8").read()
    bc_text = open(os.path.join(ref, "bcondConfig.yaml"), encoding="utf-8").read()
    try:
        cfg = yaml.safe_load(cfg_text) or {}
    except yaml.YAMLError as e:
        raise RerunError(f"solverConfig.yaml が YAML として読めない: {e}")

    # --- 入力契約: メッシュと値のファイル ---
    mesh = cfg.get("mesh") or {}
    if str(mesh.get("meshFileName")) != MESH_FILE or str(mesh.get("valueFileName")) != MESH_FILE:
        raise RerunError(f"mesh.meshFileName / valueFileName が両方 {MESH_FILE!r} でない "
                         f"({mesh.get('meshFileName')!r} / {mesh.get('valueFileName')!r}); v1 の対応外")
    # --- 入力契約: solverConfig のファイル参照は許可リストの run 内相対だけ ---
    pp = cfg.get("physProp") or {}
    for path, val in _file_refs(cfg):
        if path in ("mesh.meshFileName", "mesh.valueFileName"):
            continue
        if path == "physProp.speciesDBFile" and not os.path.isabs(val) and os.sep not in os.path.normpath(val) \
                and _allowed_name(val) and os.path.exists(os.path.join(ref, val)):
            continue
        raise RerunError(f"solverConfig.yaml の {path} = {val!r} は外部参照・未知の入力依存 (v1 の対応外)")

    # --- 入力契約: bcond ---
    bc = parse_bcond(bc_text)
    kinds = {n: e["kind"] for n, e in bc.items()}
    bad = {n: k for n, k in kinds.items() if k not in ALLOWED_KINDS}
    if bad:
        raise RerunError(f"対応外の境界種別 {bad} (v1 は {ALLOWED_KINDS} だけ)")
    inlets = [n for n, k in kinds.items() if k == "inlet_Pressure"]
    outlets = [n for n, k in kinds.items() if k == "outlet_statPress"]
    if len(inlets) != 1:
        raise RerunError(f"inlet_Pressure が {len(inlets)} 本 (v1 は単一の inlet_Pressure だけ)")
    if len(outlets) != 1:
        raise RerunError(f"outlet_statPress が {len(outlets)} 本 (v1 は単一の outlet_statPress だけ)")
    walls_iso = [n for n, k in kinds.items() if k == "wall_isothermal"]
    for n, e in bc.items():
        if "inletProfile" in e["ints"] or "inletProfile" in e["floats"] or "inletProfile" in e["entry"]:
            raise RerunError(f"{n}: inletProfile は v1 の対応外 (入口分布の CSV を引き継げない)")
    inl, outl = bc[inlets[0]], bc[outlets[0]]
    fin = inl["floats"]
    if any(re.fullmatch(r"X\d+", str(k)) for k in fin):
        raise RerunError(f"{inlets[0]}: 入口組成が X{{s}} (モル分率) 形式 — v1 は Y{{s}} 形式だけ")
    for k in ("Pt", "Tt"):
        if k not in fin:
            raise RerunError(f"{inlets[0]}: floats に {k} が無い")
    if "Ps" not in outl["floats"] or "Pt" not in outl["floats"]:
        raise RerunError(f"{outlets[0]}: floats に Ps / Pt が無い")

    # --- 化学種 ---
    names = [str(s["name"]) if isinstance(s, dict) else str(s) for s in (pp.get("species") or [])]
    tm = int(pp.get("thermalMethod", 0))
    nsp = len(names) if (tm == 2 and len(names) >= 2) else 0
    ykeys = sorted((k for k in fin if re.fullmatch(r"Y\d+", str(k))), key=lambda k: int(k[1:]))
    if nsp:
        if ykeys != [f"Y{i}" for i in range(nsp)]:
            raise RerunError(f"{inlets[0]}: 入口 floats の Y{{s}} が {ykeys} (physProp.species {names} の全種 Y0..Y{nsp - 1} が要る)")
    elif ykeys:
        raise RerunError(f"{inlets[0]}: 単成分なのに Y{{s}} がある ({ykeys})")

    # --- 参照の実効値の有限性と物理範囲 (f = Pt_new/Pt_ref などの割り算の前に) ---
    Pt_ref = _check_range(f"参照 {inlets[0]}.Pt", fin["Pt"])
    Tt_ref = _check_range(f"参照 {inlets[0]}.Tt", fin["Tt"])
    for key, zok in (("k", True), ("omega", False)):
        if key in fin:
            _check_range(f"参照 {inlets[0]}.{key}", fin[key], zero_ok=zok)
    Ps_ref = _check_range(f"参照 {outlets[0]}.Ps", outl["floats"]["Ps"])
    for key in ("Pt", "Tt"):
        if key in outl["floats"]:
            _check_range(f"参照 {outlets[0]}.{key}", outl["floats"][key])
    for n in walls_iso:
        if "Ts" in bc[n]["floats"]:
            _check_range(f"参照 {n}.Ts", bc[n]["floats"]["Ts"])
    for key in ("cfl", "cfl_pseudo"):
        if key in _deltaT(cfg):
            _check_range(f"参照 solverConfig time.deltaT.{key}", _deltaT(cfg)[key])

    changes = {}               # 条件名 → (旧, 新)
    inlet_upd, outlet_upd, wall_upd = {}, {}, {}

    def _set(key, new, old, upd, tok):
        if new is None:
            return
        if float(new) != float(old):
            changes[key] = (float(old), float(new))
            upd[tok] = fmt(new)

    _set("Pt", a.Pt, Pt_ref, inlet_upd, "Pt")
    _set("Tt", a.Tt, Tt_ref, inlet_upd, "Tt")
    if a.k is not None:
        if "k" not in fin:
            raise RerunError(f"{inlets[0]}: floats に k が無い (--k を書けない)")
        _set("k", a.k, float(fin["k"]), inlet_upd, "k")
    if a.omega is not None:
        if "omega" not in fin:
            raise RerunError(f"{inlets[0]}: floats に omega が無い (--omega を書けない)")
        _set("omega", a.omega, float(fin["omega"]), inlet_upd, "omega")

    # --- Y (全種を書く; 吸収種は 2 種なら自動・3 種以上は --balance) ---
    Y_ref = []
    for i in range(nsp):
        y = _num(fin[f"Y{i}"])
        if y is None or not np.isfinite(y):
            raise RerunError(f"参照 {inlets[0]}.Y{i} ({names[i]}) = {fin[f'Y{i}']!r} が有限の数値でない")
        Y_ref.append(y)
    Y_new = list(Y_ref)
    yreq = {}
    for item in (a.Y or []):
        if "=" not in item:
            raise RerunError(f"--Y は NAME=v の形 (got {item!r})")
        nm, v = item.split("=", 1)
        hit = [i for i, s in enumerate(names) if s == nm] or [i for i, s in enumerate(names) if s.upper() == nm.upper()]
        if not hit:
            raise RerunError(f"--Y {nm}: physProp.species {names} に無い")
        yreq[hit[0]] = float(v)
    if a.Y1 is not None:
        yreq[1] = float(a.Y1)
    if yreq and not nsp:
        raise RerunError("単成分 (または thermalMethod != 2) の run に --Y / --Y1 は使えない")
    if a.balance is not None and not yreq:
        raise RerunError("--balance は --Y / --Y1 と一緒に使う")
    bal = None
    if yreq:
        if any(i >= nsp for i in yreq):
            raise RerunError(f"--Y1: 種数 {nsp} の範囲外")
        if a.balance is not None:
            hit = [i for i, s in enumerate(names) if s == a.balance] or [i for i, s in enumerate(names) if s.upper() == a.balance.upper()]
            if not hit:
                raise RerunError(f"--balance {a.balance}: physProp.species {names} に無い")
            bal = hit[0]
            if bal in yreq:
                raise RerunError(f"--balance {a.balance} に値を指定している (吸収種は残りから決まる)")
        else:
            rest = [i for i in range(nsp) if i not in yreq]
            if nsp == 2 and len(rest) == 1:
                bal = rest[0]
            elif rest:
                raise RerunError(f"{nsp} 種なので吸収種を --balance NAME で指定すること (候補 {[names[i] for i in rest]})")
        for i, v in yreq.items():
            Y_new[i] = v
        if bal is not None:
            Y_new[bal] = 1.0 - sum(Y_new[i] for i in range(nsp) if i != bal)
        for i, v in enumerate(Y_new):
            if not (np.isfinite(v) and 0.0 <= v <= 1.0):
                raise RerunError(f"Y{i} ({names[i]}) = {v} が [0, 1] の外")
        Y_back = [float(fmt(v)) for v in Y_new]
        if abs(sum(Y_back) - 1.0) > Y_SUM_TOL:
            raise RerunError(f"ΣY = {sum(Y_back)!r} が 1 から {Y_SUM_TOL:g} より外れる (--balance で吸収種を指定するか値を見直す)")
        if Y_back != Y_ref:
            changes["Y"] = (Y_ref, Y_back)
            for i, v in enumerate(Y_back):
                inlet_upd[f"Y{i}"] = fmt(v)
        Y_new = Y_back

    if nsp and "Y" not in changes:
        # 入口 Y を変えないときも、そのまま使う参照の Y を検査する (forge の起動検査は 1e-3 で、入口カーネルが黙って正規化する;
        # codex result 段 #1: 参照 BC の ΣY = 1.0001 を受理していた)。変えるときは上で 1e-12 で検査済み
        bad = [f"Y{i} ({names[i]}) = {y!r}" for i, y in enumerate(Y_ref) if not (0.0 <= y <= 1.0)]
        if bad:
            raise RerunError(f"参照 {inlets[0]} の入口 Y が [0, 1] の外: {bad} (--Y で全種を書き直すこと)")
        if abs(sum(Y_ref) - 1.0) > REF_Y_SUM_TOL:
            raise RerunError(f"参照 {inlets[0]} の ΣY = {sum(Y_ref)!r} が 1 から {REF_Y_SUM_TOL:g} より外れる "
                             "(参照 BC の入口組成が正規化されていない; --Y で全種を書き直すこと)")

    # --- 出口 Ps (Pt を変えたら --Ps / --keep-Ps 必須) ---
    res_path = find_ref_res(ref, a.res)
    f_pt = (changes["Pt"][1] / changes["Pt"][0]) if "Pt" in changes else 1.0
    outlet_id = outl["entry"].get("physID")
    P_exit_ref, P_exit_src = exit_pressure_ref(
        ref, res_path, outlet_id,
        excluded_ids=[bc[n]["entry"].get("physID") for n, k in kinds.items() if k in ("wall", "wall_isothermal", "slip", "axis")])
    if P_exit_ref is not None and not (np.isfinite(P_exit_ref) and P_exit_ref > 0.0):
        P_exit_src = f"取れない (出口断面の P の中央値 {P_exit_ref!r} が有限の正値でない; {P_exit_src})"
        P_exit_ref = None
    p_exit_warn = None if P_exit_ref is not None else f"P_exit_ref: null — {P_exit_src}"
    if "Pt" in changes and a.Ps is None and not a.keep_Ps:
        ratio = (Ps_ref / (f_pt * P_exit_ref)) if P_exit_ref else None
        raise RerunError(
            "Pt を変えるときは --Ps P か --keep-Ps が必須 (node の出口は壁列が常に亜音速で Ps を見る; Pt だけ下げて Ps 据え置きは"
            " 出口列の不安定要因。plan §4.1)\n"
            f"    f = Pt_new/Pt_ref = {f_pt:.6g}, P_exit_ref = {P_exit_ref} [{P_exit_src}], Ps_ref = {Ps_ref}\n"
            f"    参照の Ps_ref/P_exit_ref = {(Ps_ref / P_exit_ref) if P_exit_ref else None}\n"
            f"    --keep-Ps のとき Ps/(f·P_exit_ref) = {ratio if ratio is None else f'{ratio:.6g}'}\n"
            f"    参照と同じ比にするなら --Ps {f_pt * Ps_ref:.10g} (= f·Ps_ref)")
    if a.Ps is not None and float(a.Ps) != Ps_ref:
        changes["Ps"] = (Ps_ref, float(a.Ps))
        outlet_upd["Ps"] = fmt(a.Ps)
        outlet_upd["Pt"] = fmt(a.Ps)        # 出口の逆流用 Pt も Ps と同時に書く (Tt は据え置き; plan §4.3)
    Ps_new = float(a.Ps) if a.Ps is not None else Ps_ref
    ps_ratio = (Ps_new / (f_pt * P_exit_ref)) if P_exit_ref else None

    # --- 等温壁 ---
    if a.Tw is not None and not walls_iso:
        raise RerunError("--Tw を指定したが wall_isothermal の境界が無い")
    if walls_iso and "Tt" in changes and a.Tw is None and not a.keep_Tw:
        raise RerunError(f"等温壁 {walls_iso} があり Tt を変える: --Tw T か --keep-Tw を指定すること")
    wall_old = {}
    for n in walls_iso:
        if "Ts" not in bc[n]["floats"]:
            raise RerunError(f"{n}: wall_isothermal の floats に Ts が無い")
        wall_old[n] = float(bc[n]["floats"]["Ts"])
        if a.Tw is not None and float(a.Tw) != wall_old[n]:
            changes.setdefault("Tw", {})[n] = (wall_old[n], float(a.Tw))
            wall_upd[n] = {"Ts": fmt(a.Tw)}

    # --- 初期場スケール (opt-in) ---
    scale = a.scale_ic
    if scale == "pt":
        why = []
        if "Pt" not in changes:
            why.append("Pt を変えていない (f = 1)")
        if "condensation" in cfg:
            why.append("solverConfig に condensation block がある (凝縮モーメントは状態変換で保てない)")
        if "Tt" in changes:
            why.append("Tt を変える")
        if "Y" in changes:
            why.append("入口組成を変える")
        if why:
            raise RerunError("--scale-ic pt は使えない: " + "; ".join(why) + " (plan §4.4(3))")

    # --- 必要保存量と場の検査 (SRC = 参照 res, DST = 参照 nozzle.h5) ---
    required = required_conserved_from_cfg(cfg)
    probs, n_src = field_problems(res_path, required)
    if probs:
        raise RerunError(f"SRC {res_path} の必要保存量が不正:\n    " + "\n    ".join(probs))
    dst_ref = os.path.join(ref, MESH_FILE)
    probs, n_dst = field_problems(dst_ref, required)
    if probs:
        raise RerunError(f"DST {dst_ref} の必要保存量が不正:\n    " + "\n    ".join(probs))
    if n_src != n_dst:
        raise RerunError(f"SRC と DST の点数が違う ({n_src} vs {n_dst}) — 同一格子でない")
    import h5py
    with h5py.File(dst_ref, "r") as fh:
        extra = sorted(set(fh["VALUE"]) - set(required) - set(KEEP_FROM_DST))
    with h5py.File(res_path, "r") as fh:
        src_dtypes = {k: str(fh["VALUE"][k].dtype) for k in required}
        src_species_hash = fh.attrs.get("species_hash")
    # 乱流モデルなし (Euler・層流) の run では、変換器が既定で作る roK/roOmega の入れ物が残っている (forge は読まない)。
    # それだけは「使われない余りの量」として許し、記録する (スケールもしない)。他の未知の量は従来どおり拒否 (2026-10-06)。
    tb_model = str((cfg.get("turbulence") or {}).get("model") or "none").strip().lower()
    inert = [k for k in extra if k in INERT_WHEN_NO_TURB and not tb_model.startswith("sst")]
    extra = [k for k in extra if k not in inert]
    inert_note = (f"乱流モデルなしの run の未使用量 {inert} は restart_field がそのまま写す (forge は読まない・スケールしない)" if inert else None)
    if extra:
        raise RerunError(f"DST {dst_ref} の /VALUE に必要保存量と wall_dist 以外がある {extra} "
                         "(restart_field は SRC に同名があれば写すので、設定と食い違う量が持ち込まれうる; codex M1)")

    # --- solverConfig の書き換え (step 数・出力間隔・CFL) ---
    new_cfg_text = cfg_text
    cfg_changes = {}
    if a.steps is not None:
        new_cfg_text = _sub_once(r"(nStepOuter:\s*)\d+", lambda m: m.group(1) + str(int(a.steps)), new_cfg_text, "nStepOuter")
        cfg_changes["nStepOuter"] = int(a.steps)
    if a.out_interval is not None:
        new_cfg_text = _sub_once(r"(outStepInterval:\s*)\d+", lambda m: m.group(1) + str(int(a.out_interval)),
                                 new_cfg_text, "outStepInterval")
        cfg_changes["outStepInterval"] = int(a.out_interval)
    if a.cfl is not None:
        new_cfg_text = _sub_once(r"\bcfl:\s*[-+.\deE]+,\s*cfl_pseudo:\s*[-+.\deE]+",
                                 lambda m: f"cfl: {fmt(a.cfl)}, cfl_pseudo: {fmt(a.cfl)}", new_cfg_text, "`cfl: X, cfl_pseudo: X`")
        cfg_changes["cfl"] = float(a.cfl)
    new_cfg = yaml.safe_load(new_cfg_text) or {}
    eff_cfl = {k: _num(_deltaT(new_cfg).get(k)) for k in ("cfl", "cfl_pseudo")}
    if a.cfl is not None and any(v is not None and v != float(a.cfl) for v in eff_cfl.values()):
        raise RerunError(f"書き換え後の solverConfig を読み直したら time.deltaT の cfl/cfl_pseudo = {eff_cfl} (期待 {float(a.cfl)})")
    for k, v in eff_cfl.items():
        if k in _deltaT(new_cfg):
            _check_range(f"生成 solverConfig time.deltaT.{k}", v)
    tt = new_cfg.get("time") or {}
    n_outer = int(((tt.get("last") or {}).get("nStepOuter")) or 0)
    n_out = int(tt.get("outStepInterval") or 0)
    if n_out <= 0 or n_outer <= 0 or n_outer % n_out != 0:
        raise RerunError(f"nStepOuter {n_outer} が outStepInterval {n_out} の倍数でない (最終 res が書かれない)")
    ref_outer = int((((cfg.get("time") or {}).get("last") or {}).get("nStepOuter")) or 0)
    warnings = [p_exit_warn] if p_exit_warn else []
    if inert_note:
        warnings.append(inert_note)
    if _res_step(res_path) is not None and _res_step(res_path) != ref_outer:
        warnings.append(f"参照 res の step {_res_step(res_path)} が参照 config の nStepOuter {ref_outer} と違う "
                        "(途中の res を種にしている / config が延長後のものでない)")

    # --- bcond の書き換え (対象行の当該トークンだけ) ---
    lines = bc_text.splitlines(keepends=True)
    for name, upd in ((inlets[0], inlet_upd), (outlets[0], outlet_upd), *((n, u) for n, u in wall_upd.items())):
        if upd:
            i = bc[name]["line"]
            eol = lines[i][len(lines[i].rstrip("\r\n")):]
            lines[i] = replace_float_tokens(lines[i].rstrip("\r\n"), upd) + eol
    new_bc_text = "".join(lines)
    # 読み直して要求値と一致することを検査
    nb = parse_bcond(new_bc_text)
    want = {inlets[0]: inlet_upd, outlets[0]: outlet_upd, **wall_upd}
    for name, e in bc.items():
        for k, v in e["floats"].items():
            exp = float(want.get(name, {}).get(k, v)) if v is not None else None
            got = nb[name]["floats"].get(k)
            if (exp is None) != (got is None) or (exp is not None and float(got) != exp):
                raise RerunError(f"書き換え後の bcond を読み直したら {name}.{k} = {got} (期待 {exp})")
        if set(nb[name]["floats"]) != set(e["floats"]) or nb[name]["kind"] != e["kind"]:
            raise RerunError(f"書き換え後の bcond の {name} のキー・種別が変わった")
    if nsp:
        yb = [float(nb[inlets[0]]["floats"][f"Y{i}"]) for i in range(nsp)]
        tol = Y_SUM_TOL if "Y" in changes else REF_Y_SUM_TOL
        if any(not (0.0 <= y <= 1.0) for y in yb) or abs(sum(yb) - 1.0) > tol:
            raise RerunError(f"書き換え後の入口 Y = {yb} (各成分 [0, 1]・|ΣY − 1| ≤ {tol:g} でない)")

    # --- species_meta の同期 ---
    meta_text_new, meta_note = None, None
    meta_path = os.path.join(ref, "species_meta.yaml")
    if "Y" in changes and os.path.exists(meta_path):
        meta_text_new, meta_note = _sync_species_meta(open(meta_path, encoding="utf-8").read(), names, Y_new)

    cond_changed = bool(changes)
    info_path = os.path.join(ref, "prepare_info.json")
    pinfo = json.load(open(info_path)) if os.path.exists(info_path) else None
    # NS / Euler の分類は実効 config から (codex result-2 2026-10-06 Major 2): 壁 BC が粘着 (wall / wall_isothermal) なら NS、
    # 全部 slip なら Euler。乱流モデル (sst*) があれば NS。prepare_info.viscous は照合だけに使い、食い違えば作成前に止める。
    wall_kinds = [k for k in kinds.values() if k in ("wall", "wall_isothermal", "slip")]
    if not wall_kinds:
        raise RerunError("壁 BC (wall / wall_isothermal / slip) が無いので NS / Euler を判別できない")
    noslip = any(k in ("wall", "wall_isothermal") for k in wall_kinds)
    if noslip and any(k == "slip" for k in wall_kinds):
        raise RerunError(f"粘着壁と slip 壁が混在 ({wall_kinds}) — NS / Euler を判別できない (v1 対象外)")
    is_ns = noslip or tb_model.startswith("sst")
    if tb_model.startswith("sst") and not noslip:
        raise RerunError("乱流モデル (sst) なのに壁が全部 slip — config が不整合")
    if pinfo is not None and pinfo.get("viscous") is not None and bool(pinfo["viscous"]) != is_ns:
        raise RerunError(f"prepare_info.viscous={pinfo['viscous']} と実効 config (壁 {wall_kinds}・乱流 {tb_model}) の分類が食い違う")
    scale_allowed = not (("condensation" in cfg) or ("Tt" in changes) or ("Y" in changes) or ("Tw" in changes))
    if cond_changed and is_ns and "Pt" in changes and not scale_allowed:
        # Pt と Tt・組成・壁温・凝縮の複合変更: scale-ic pt は使えない (§4.4 の禁止条件) — 推奨は本段 cfl 1 だけ、整定の実績なし
        rec = {"stages": "full", "cfl": 1.0, "nStepOuter_min": 60000,
               "note": "Pt と Tt/組成/壁温/凝縮の複合変更 → run_staged_ns(stages='full')・本段 cfl 1・60000 step。"
                       "scale-ic pt は禁止条件に当たるので推奨しない。複合変更の整定の実績はない (plan §4.7)"}
        warnings.append("Pt と Tt/組成/壁温/凝縮を同時に変えた: この組合せの起動・整定は未検証 (plan §4.7)。量が STEADY になるまで延長して確認すること")
    elif cond_changed and is_ns and "Pt" in changes:
        # plan tooling-rerun-conditions §6 (ii′)・A3 (2026-10-06): Pt 0.8 倍は stages full でも本段 cfl 5 で出口壁際の角から発散
        # (scale あり step 468、なし step 2 で入口)、scale あり + full + 本段 cfl 1 は STEADY → Pt 変更の本段は cfl 1、scale-ic pt を推奨
        rec = {"stages": "full", "cfl": 1.0, "nStepOuter_min": 60000,
               "note": "Pt を変えた → run_staged_ns(stages='full') (soft→mid→本段) で本段 cfl 1・60000 step、"
                       "--scale-ic pt を推奨 (plan §4.7、§6 (ii′)・A3)"}
        if a.scale_ic != "pt":
            warnings.append("Pt 変更で --scale-ic none: 検証では stages full・本段 cfl 5 で入口から step 2 で発散した (plan §6 (ii′) 腕 B2)。"
                            "--scale-ic pt を推奨 (none は本段 cfl 1 でも入口配管の壁際に逆流域を残して別の状態に向かった: §6 (ii″) B3)")
    elif cond_changed:
        rec = {"stages": "full",
               "note": "条件を変えた → run_staged_ns(stages='full') (soft→mid→本段)。細分格子は本段 cfl 1 (plan §4.7)"}
    else:
        rec = {"stages": "none", "cfl": cfg_changes.get("cfl", _ref_cfl(cfg)),
               "note": "条件の変更なし → 参照場からの継続 (stages='none', 参照 cfl)"}
    # --- 推奨と生成 config の整合 (codex result 段 #3: 表示する run_staged_ns の行は生成 config の cfl で本段を回す) ---
    # §4.7 の推奨は NS (run_staged_ns) の実測に基づく。Euler 参照 (δ_E 用の対; plan §4.8) は run_staged で回し、
    # §6 (ii) 腕 E は stages none・cfl 2・6000 で合格している — §4.7 は Euler の推奨を定めていないので整合検査をしない。
    if cond_changed and not is_ns:
        # Euler 参照 (δ_E 用の対; plan §4.8)。実績があるのは Pt だけの変更: §6 (ii) 腕 E (Pt 0.8、scale-ic pt) は
        # run_staged(stages="none")・cfl 2・6000 step で STEADY (run_0120)。Tt・組成の変更 (run_0134: Tt 1500・H2O 0.10) は
        # 同じ設定で 6000 step では DRIFTING (出口 M 6.022、直前窓差 0.021) — 推奨は未確立 (2026-10-06 監査で判明)
        ps_scaled = ("Pt" in changes and abs(Ps_new / Ps_ref - f_pt) <= 1e-3 * abs(f_pt))
        if set(changes) <= {"Pt", "Ps"} and a.scale_ic == "pt" and ps_scaled:
            # 検証済みの条件に限る (codex result-2 Major 3): 腕 E は保存量と背圧をともに f 倍 (scale-ic pt・Ps = f·Ps_ref)
            rec = {"stages": "none", "cfl": 2.0,
                   "note": "Euler 参照・Pt のみ変更 (scale-ic pt・Ps も同じ比) → run_staged(stages='none')・cfl 2・6000 step (plan §4.7: 腕 E run_0120 の実績)"}
        elif set(changes) <= {"Pt", "Ps"}:
            rec = {"stages": "full", "cfl": 2.0,
                   "note": "Euler 参照・Pt 変更だが scale-ic pt でない、または Ps が Pt と同じ比でない → 未検証。run_staged(stages='full') で STEADY まで (plan §4.7)"}
            warnings.append("Euler 参照の Pt 変更で、検証済み条件 (scale-ic pt かつ Ps = f·Ps_ref) と違う: 起動の実績なし。STEADY を確認すること")
        else:
            rec = {"stages": "full", "cfl": 2.0,
                   "note": "Euler 参照・Tt/組成の変更 → 推奨は未確立 (run_0134 は stages none・cfl 2・6000 で DRIFTING)。"
                           "run_staged(stages='full') で回し、量が STEADY になるまで延長すること (plan §4.7)"}
            warnings.append("Euler 参照で Tt・組成を変えた: 検証では stages none・cfl 2・6000 step で準定常に達しなかった (run_0134)。"
                            "δ_E の参照に使う前に STEADY を確認すること")
    rec["runner"] = "run_staged_ns" if is_ns else "run_staged"
    mism, fix = [], []
    if rec.get("cfl") is not None and any(v is not None and v != _num(rec["cfl"]) for v in eff_cfl.values()):
        mism.append(f"本段 cfl/cfl_pseudo = {eff_cfl['cfl']}/{eff_cfl['cfl_pseudo']} (推奨 {rec['cfl']})")
        fix.append(f"--cfl {float(rec['cfl'])}")
    if rec.get("nStepOuter_min") is not None and n_outer < int(rec["nStepOuter_min"]):
        mism.append(f"nStepOuter = {n_outer} (推奨 ≥ {rec['nStepOuter_min']})")
        fix.append(f"--steps {int(rec['nStepOuter_min'])}")
        if int(rec["nStepOuter_min"]) % n_out:
            fix.append("--out-interval <nStepOuter の約数>")
    if mism and not is_ns:
        warnings.append("Euler 参照 (run_staged で回す): 推奨 (stages none・cfl 2、plan §4.7 の Euler 行) と生成 config が食い違う — "
                        "Euler は停止せず警告だけ (--cfl 2.0 で推奨どおり)。食い違い: " + "; ".join(mism))
        mism = []
    if mism:
        msg = ("推奨 (recommended_stages: " + json.dumps(rec, ensure_ascii=False) + ") と生成 config が食い違う: "
               + "; ".join(mism))
        if not a.override_recommended:
            raise RerunError(msg + f"\n    推奨どおりにするなら {' '.join(fix)} を足す。"
                             "意図して違う config で作るなら --override-recommended (記録に残る)")
        warnings.append(msg + " — --override-recommended で推奨と違う config のまま作った")
    rec["config_effective"] = {"cfl": eff_cfl["cfl"], "cfl_pseudo": eff_cfl["cfl_pseudo"], "nStepOuter": n_outer,
                               "outStepInterval": n_out, "override": bool(mism)}
    euler_note = None
    if cond_changed and pinfo and pinfo.get("viscous"):
        euler_note = ("δ_E の評価には同条件の Euler rerun を対で作ること (Euler 参照 run にも同じ条件引数で rerun_conditions)。"
                      "旧条件の Euler 参照に対する mdot_ratio_vs_euler は診断量として記録するだけ (plan §4.8)")

    return {
        "ref": ref, "new": new, "res": res_path, "cfg_text": cfg_text, "bc_text": bc_text,
        "new_cfg_text": new_cfg_text, "new_bc_text": new_bc_text, "meta_text_new": meta_text_new, "meta_note": meta_note,
        "changes": changes, "cfg_changes": cfg_changes, "species": names, "Y_ref": Y_ref, "Y_new": Y_new,
        "balance": names[bal] if bal is not None else None,
        "inlet": inlets[0], "outlet": outlets[0], "walls_isothermal": walls_iso,
        "bc_before": {n: e["floats"] for n, e in bc.items()}, "bc_after": {n: e["floats"] for n, e in nb.items()},
        "required": required, "src_dtypes": src_dtypes, "src_species_hash": src_species_hash,
        "scale_ic": scale, "f": f_pt, "P_exit_ref": P_exit_ref, "P_exit_ref_source": P_exit_src,
        "Ps": Ps_new, "Ps_over_fPexit": ps_ratio, "keep_Ps": bool(a.keep_Ps), "keep_Tw": bool(a.keep_Tw),
        "recommended_stages": rec, "euler_note": euler_note, "warnings": warnings, "prepare_info": pinfo,
        "forge": a.forge,
    }


def _ref_cfl(cfg):
    dT = ((cfg.get("time") or {}).get("deltaT") or {})
    return dT.get("cfl")


def _sync_species_meta(text, names, Y):
    """species_meta.yaml の streams.inflow.Y_transport・Y (実種; expansion で展開) を同期し、
    実種の MW が内蔵表で揃うときだけ X も作り直す (揃わなければ X を消す)。戻り (新しい本文, 注記)。"""
    meta = yaml.safe_load(text) or {}
    if [str(s) for s in (meta.get("species") or names)] != names:
        raise RerunError(f"species_meta.yaml の species {meta.get('species')} が physProp.species {names} と違う")
    inflow = ((meta.get("streams") or {}).get("inflow"))
    if not isinstance(inflow, dict):
        raise RerunError("species_meta.yaml に streams.inflow が無い (Y_transport を同期できない)")
    inflow["Y_transport"] = [float(v) for v in Y]
    exp = meta.get("expansion") or {}
    note = []
    if all(n in exp for n in names):
        real = {}
        for n, y in zip(names, Y):
            for k, w in (exp[n] or {}).items():
                real[str(k)] = real.get(str(k), 0.0) + float(y) * float(w)
        order = list((inflow.get("Y") or {}).keys()) or list(real)
        inflow["Y"] = {k: real.get(k, 0.0) for k in order + [k for k in real if k not in order]}
        mws = {k: fsp._find_ci(fsp.BUILTIN_MW, k) for k in inflow["Y"]}
        if all(v is not None for v in mws.values()):
            nmol = {k: inflow["Y"][k] / mws[k] for k in inflow["Y"]}
            tot = sum(nmol.values())
            inflow["X"] = {k: nmol[k] / tot for k in nmol}
            note.append("X は forge_species.BUILTIN_MW で作り直した")
        else:
            inflow.pop("X", None)
            note.append(f"実種の MW が揃わない ({[k for k, v in mws.items() if v is None]}) ので X を消した")
    else:
        inflow.pop("Y", None)
        inflow.pop("X", None)
        note.append("expansion が無い種があるので実種の Y・X を消した (Y_transport だけ同期)")
    return yaml.safe_dump(meta, sort_keys=False, allow_unicode=True), "; ".join(note)


# ---------------------------------------------------------------------------
# 実行 (NEW_RUN を作る; 失敗したら消す)
# ---------------------------------------------------------------------------
def scale_fields(dst_h5, src_h5, required, f):
    """dst の必要保存量を f·src にする (T・U・Y・k・ω を保つ初期場変換)。検査 allclose(d_new, f·d_ref, rtol, atol=0)
    (比で検査しない — ゼロ成分で NaN)。戻り {量: 最大相対誤差}。不合格は RerunError。"""
    import h5py
    out = {}
    with h5py.File(src_h5, "r") as s, h5py.File(dst_h5, "r+") as d:
        for name in required:
            ref = np.asarray(s["VALUE"][name], dtype=np.float64)
            ds = d["VALUE"][name]
            ds[...] = (f * ref).astype(ds.dtype)
        d.flush()
        for name in required:
            ref = np.asarray(s["VALUE"][name], dtype=np.float64)
            got = np.asarray(d["VALUE"][name], dtype=np.float64)
            want = f * ref
            if not np.all(np.isfinite(got)) or not np.allclose(got, want, rtol=SCALE_RTOL, atol=0.0):
                raise RerunError(f"スケール検査 NG: {name} が f·d_ref に rtol {SCALE_RTOL:g} で一致しない")
            nz = want != 0.0
            out[name] = float(np.max(np.abs(got[nz] - want[nz]) / np.abs(want[nz]))) if nz.any() else 0.0
    # allclose は掛け算の正確さしか見ないので、変換後の場の物理的妥当性 (有限・ρ>0・Y の範囲) を改めて検査する
    # (codex result 段 #1: f < 0 で全点負密度でも合格していた)
    probs, _ = field_problems(dst_h5, required)
    if probs:
        raise RerunError(f"スケール後 (f = {f!r}) の初期場が不正:\n    " + "\n    ".join(probs))
    return out


def _git_info():
    def run(*args):
        try:
            return subprocess.run(["git", "-C", REPO_ROOT, *args], capture_output=True, text=True).stdout.strip()
        except OSError:
            return ""
    return {"head": run("rev-parse", "HEAD") or "unknown",
            "tool_dirty": bool(run("status", "--porcelain", "--", os.path.relpath(os.path.abspath(__file__), REPO_ROOT)))}


def execute(plan):
    """計画どおりに NEW_RUN を作る。途中で失敗したら NEW_RUN ごと消して RerunError を投げ直す。"""
    ref, new = plan["ref"], plan["new"]
    os.makedirs(new)              # 既にあれば失敗 (build_plan でも検査済み)
    try:
        copied = []
        for fn in sorted(os.listdir(ref)):
            if _allowed_name(fn) and os.path.isfile(os.path.join(ref, fn)):
                shutil.copy2(os.path.join(ref, fn), os.path.join(new, fn))
                copied.append(fn)
        with open(os.path.join(new, "bcondConfig.yaml"), "w", encoding="utf-8") as f:
            f.write(plan["new_bc_text"])
        with open(os.path.join(new, "solverConfig.yaml"), "w", encoding="utf-8") as f:
            f.write(plan["new_cfg_text"])
        if plan["meta_text_new"] is not None:
            with open(os.path.join(new, "species_meta.yaml"), "w", encoding="utf-8") as f:
                f.write(plan["meta_text_new"])

        # --- 初期場: restart_field (保存量の index コピー + 化学種の照合) ---
        dst = os.path.join(new, MESH_FILE)
        cmd = [sys.executable, os.path.join(HERE, "restart_field.py"), plan["res"], dst, "--dst-run", new]
        if any(v == "float64" for v in plan["src_dtypes"].values()):
            cmd.append("--keep-src-dtype")
        if plan["forge"]:
            cmd += ["--forge", plan["forge"]]
        r = subprocess.run(cmd, capture_output=True, text=True)
        with open(os.path.join(new, "restart_field.log"), "w", encoding="utf-8") as f:
            f.write("$ " + " ".join(cmd) + "\n" + r.stdout + r.stderr)
        verdict = [ln for ln in r.stdout.splitlines() if ln.startswith("VERDICT: OK")]
        if r.returncode != 0 or not verdict:
            raise RerunError(f"restart_field.py が失敗 (rc={r.returncode}):\n{(r.stdout + r.stderr)[-1500:]}")
        import h5py
        with h5py.File(dst, "r") as fh:
            got_hash = fh.attrs.get("species_hash")
        if plan["src_species_hash"] is not None and got_hash != plan["src_species_hash"]:
            raise RerunError(f"restart 後の species_hash {got_hash!r} が SRC {plan['src_species_hash']!r} と違う "
                             "(未検証のまま写された; FORGE_ALLOW_UNVERIFIED_SPECIES 等で通していないか)")
        probs, _ = field_problems(dst, plan["required"])
        if probs:
            raise RerunError("restart 後の初期場が不正:\n    " + "\n    ".join(probs))

        scale_check = None
        if plan["scale_ic"] == "pt":
            scale_check = scale_fields(dst, plan["res"], plan["required"], plan["f"])

        # --- forge (化学種の解決に使ったもの) ---
        try:
            exe = fsp.find_forge(plan["forge"])
        except Exception as e:  # noqa: BLE001  記録だけ (restart_field は通っている)
            exe, _ = None, e
        forge_rec = {"path": exe, "sha256": _sha256(exe) if exe and os.path.exists(exe) else None}

        # --- prepare_info.json: 幾何は据え置き、ic_from と rerun_of ---
        if plan["prepare_info"] is not None:
            pinfo = dict(plan["prepare_info"])
            pinfo["ic_from"] = _repo_rel(plan["res"])
            pinfo["rerun_of"] = _repo_rel(ref)
            with open(os.path.join(new, "prepare_info.json"), "w") as f:
                json.dump(pinfo, f, indent=1)

        record = {
            "tool": "solver_density_cuda/tools/rerun_conditions.py", "plan": "plans/active/tooling-rerun-conditions.md",
            "argv": sys.argv[1:], "ref_run": _repo_rel(ref), "ref_res": _repo_rel(plan["res"]),
            "copied": copied,
            "changes": {k: v for k, v in plan["changes"].items()}, "solverConfig_changes": plan["cfg_changes"],
            "bcond_before": plan["bc_before"], "bcond_after": plan["bc_after"],
            "species": plan["species"], "Y_ref": plan["Y_ref"], "Y_new": plan["Y_new"], "balance": plan["balance"],
            "species_meta_sync": plan["meta_note"],
            "required_conserved": plan["required"],
            "scale_ic": plan["scale_ic"], "f": plan["f"], "scale_check_max_rel": scale_check,
            "P_exit_ref": plan["P_exit_ref"], "P_exit_ref_source": plan["P_exit_ref_source"],
            "Ps": plan["Ps"], "Ps_over_f_P_exit_ref": plan["Ps_over_fPexit"], "keep_Ps": plan["keep_Ps"],
            "keep_Tw": plan["keep_Tw"],
            "restart_field_verdict": verdict[-1], "species_hash": got_hash if got_hash is None else str(got_hash),
            "forge": forge_rec, "tool_commit": _git_info(),
            "recommended_stages": plan["recommended_stages"], "euler_reference": plan["euler_note"],
            "warnings": plan["warnings"],
        }
        with open(os.path.join(new, "RERUN_CONDITIONS.json"), "w", encoding="utf-8") as f:
            json.dump(record, f, indent=1, ensure_ascii=False, default=str)
        return record
    except BaseException:
        shutil.rmtree(new, ignore_errors=True)
        raise


def make_parser():
    ap = argparse.ArgumentParser(description="既存 run の形状を固定して入口条件・背圧・入口分率だけ変えた run を作る (forge は起動しない)")
    ap.add_argument("ref_run")
    ap.add_argument("new_run")
    ap.add_argument("--res", help="初期場の参照 res (既定: res_[0-9]*.h5 の最大番号)")
    ap.add_argument("--Pt", type=float, help="入口全圧 [Pa] (変えたら --Ps か --keep-Ps が必須)")
    ap.add_argument("--Tt", type=float, help="入口全温 [K]")
    yg = ap.add_mutually_exclusive_group()
    yg.add_argument("--Y", action="append", metavar="NAME=v", help="入口の輸送種の質量分率 (physProp.species の名前)")
    yg.add_argument("--Y1", type=float, help="入口 Y1 (輸送種 index 1) の質量分率")
    ap.add_argument("--balance", help="吸収種 (3 種以上で必須; 2 種なら自動)")
    ap.add_argument("--k", type=float, help="入口 k (既定は据え置き)")
    ap.add_argument("--omega", type=float, help="入口 omega (既定は据え置き)")
    pg = ap.add_mutually_exclusive_group()
    pg.add_argument("--Ps", type=float, help="出口背圧 [Pa] (出口の Ps と Pt を同時に書く)")
    pg.add_argument("--keep-Ps", action="store_true", help="Pt を変えても出口 Ps を据え置く")
    tg = ap.add_mutually_exclusive_group()
    tg.add_argument("--Tw", type=float, help="等温壁の壁温 Ts [K]")
    tg.add_argument("--keep-Tw", action="store_true", help="Tt を変えても等温壁の壁温を据え置く")
    ap.add_argument("--steps", type=int, help="nStepOuter")
    ap.add_argument("--out-interval", type=int, help="outStepInterval")
    ap.add_argument("--cfl", type=float, help="`cfl: X, cfl_pseudo: X` を書き換える")
    ap.add_argument("--scale-ic", choices=("none", "pt"), default="none",
                    help="pt: 必要保存量を全部 Pt_new/Pt_ref 倍 (T・U・Y・k・ω を保つ初期場変換; 明示 opt-in)")
    ap.add_argument("--override-recommended", action="store_true",
                    help="recommended_stages (plan §4.7) の本段 cfl・step 数と違う config のまま作る (既定は停止; 記録に残る)")
    ap.add_argument("--lump", default=None, help="(v1 では停止) 乾き成分 lump の組成変更")
    ap.add_argument("--forge", help="--resolve-species を持つ forge (restart_field に渡す; 既定 FORGE_BIN)")
    ap.add_argument("--dry-run", action="store_true", help="検査と書き換え内容の表示だけ (何も作らない)")
    return ap


def _summary(plan):
    out = {"ref_res": plan["res"], "changes": plan["changes"], "solverConfig_changes": plan["cfg_changes"],
           "required_conserved": plan["required"], "scale_ic": plan["scale_ic"], "f": plan["f"],
           "P_exit_ref": plan["P_exit_ref"], "P_exit_ref_source": plan["P_exit_ref_source"],
           "Ps_over_f_P_exit_ref": plan["Ps_over_fPexit"], "recommended_stages": plan["recommended_stages"],
           "euler_reference": plan["euler_note"], "warnings": plan["warnings"]}
    return json.dumps(out, indent=1, ensure_ascii=False, default=str)


def main(argv=None):
    a = make_parser().parse_args(argv)
    try:
        plan = build_plan(a)
    except RerunError as e:
        print(f"[rerun_conditions] STOP (NEW_RUN は作っていない): {e}", file=sys.stderr)
        return 2
    for w in plan["warnings"]:
        print(f"[rerun_conditions] warning: {w}", file=sys.stderr)
    if a.dry_run:
        print(_summary(plan))
        print("--dry-run: 何も作っていない")
        return 0
    try:
        rec = execute(plan)
    except RerunError as e:
        print(f"[rerun_conditions] FAILED (NEW_RUN を消した): {e}", file=sys.stderr)
        return 1
    print(_summary(plan))
    print(f"restart_field: {rec['restart_field_verdict']}")
    print(f"作成: {plan['new']} (記録 RERUN_CONDITIONS.json)")
    st = plan["recommended_stages"]["stages"]
    eff = plan["recommended_stages"]["config_effective"]
    print(f"回し方: {plan['recommended_stages']['runner']}(Path({_repo_rel(plan['new'])!r}), stages={st!r}) "
          "(design/forge_design/evaluate/runner_axismach.py) — 本段は生成 config の "
          f"cfl {eff['cfl']}・nStepOuter {eff['nStepOuter']} で回る")
    if plan["euler_note"]:
        print(f"注意: {plan['euler_note']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
