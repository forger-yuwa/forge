"""⑤ SERN 評価 runner (S3): 逆設計 → 2 バンド構造メッシュ → forge 平面 2D (Euler/SST) → 力係数。

plan: plans/active/tooling-nozzle-sern-chain.md §4.2, §4.7。問題型 `sern_2d`。
1 run = 1 作動点。作動点は spec.external で与える (多作動点束ねは S6 で driver 側)。
段階起動: soft (1 次 + cfl 0.5, 3000 step) → 本段 (2 次 + cfl_main)。IC は領域別一様
(中間線より上 = 燃焼器出口状態、下 = 外部流)。
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

from ..gas.frozen import FrozenGas, ideal_gross_thrust_frozen
from ..geometry.moc_sern import PlanarMOC, SernKernelSpec, wall_forces
from ..geometry.rao_planar import ideal_gross_thrust
from ..meshing.mesh_sern import PHYS_SERN, SernMeshParams, generate_sern_mesh, write_msh41_named
from ..metrics.sern_forces import force_history, write_force_history_csv
from ..metrics.sern_gates import evaluate_gates, forge_rc_from_log
from ..probdef import Problem, dv_value, load_problem
from .ic import _forge_species

# リポジトリ位置から導く (AWS など別マシンでも動くように。FORGE_ROOT で上書き可)
FORGE_ROOT = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
FORGE_TOOLS = FORGE_ROOT / "solver_density_cuda" / "tools"
FORGE_BUILD = FORGE_ROOT / "solver_density_cuda" / "build"
_ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
MESH = "sern.h5"



# 既定が変わった日。設計 DB で旧既定の評価と混ぜないための識別子。
#   2026-09-26: slauWallNormalChi の既定 0 → auto (node+SLAU で 1)
#   2026-09-27: mesh.scalarGradient の node 既定 gg → lsq (plan gradient-scalar-lsq-unification #6)。
#               日付の一致だけでは旧評価の混入を防げないので、学習側は実効 scalarGradient も見る (codex diagnose 2026-09-27)
FLAG_POLICY = "2026-09-27"


def _last_launch_value(run_dir, key, allowed=None):
    """forge_launches.jsonl の**最後の非空行** (= 最後の起動) の実効値 (文字列)。無い・壊れている・キーが無い・
    allowed に無い値は None (= 不明)。**前の起動の値に遡らない** (codex result 2026-09-27 M1: 最後の起動が
    読めないときに過去の lsq を引き継いで学習に採用していた)。"""
    p = Path(run_dir) / "forge_launches.jsonl"
    if not p.exists():
        return None
    lines = [l for l in p.read_text().splitlines() if l.strip()]
    if not lines:
        return None
    try:
        v = str(json.loads(lines[-1])[key])
    except Exception:
        return None
    return v if (allowed is None or v in allowed) else None


def _last_launch_chi(run_dir):
    """最後の起動の slauWallNormalChi 実効値 (0/1)。無い・壊れている・キーが無い・0/1 以外は None (= 不明)。
    前の起動の値に遡らない (codex result 2026-09-27 chi-default M4、scalarGradient と同じ厳格さ)。"""
    v = _last_launch_value(run_dir, "slauWallNormalChi", allowed=("0", "1"))
    return None if v is None else int(v)


# --- 厚さ 0 の板の自由端の近傍で速度の再構成を節点値にする処置 -------------------------------------------------------
# plan convection-zero-thickness-edge-reconstruction §4「設計チェーンへの配線」(codex plan レビュー 2026-10-08 M3)、
# 2026-10-08 のユーザ判断 (端の判定はソルバでなく前処理) と codex diagnose 2026-10-08 (edge-weight-preprocessing) の Major:
# 問題 YAML の `evaluate.zero_thickness_edge_velocity: {tags: [cowl_in, cowl_out], rings: 2}` は**前処理の入力**で、
# prepare が**最終の node 変換の後**に道具で meshFileName の `/AUX/w_recon_vel` に節点の重み w を書く (格子を作るたびに
# 作り直す)。各段の solverConfig の space には有効/フィールド名だけを出す。未指定 = 無効 (config は今とバイト一致)。
# 区間識別と来歴には、各段の起動前に meshFileName から**再計算した**格子署名と w のハッシュを入れる (属性のハッシュは
# 転記しない)。読み出しと再計算の正本は tools/stage_manifest.py の zte_* (計算そのものは道具の関数)。
ZTE_FIELD = "w_recon_vel"
ZTE_SPACE_FRAGMENT = f", zeroThicknessEdgeVelocity: {{field: {ZTE_FIELD}}}"
ZTE_TOOL = FORGE_TOOLS / "mark_zero_thickness_edges.py"
ZTE_LAUNCHES = "zte_launches.jsonl"     # 処置を有効にした段の起動前の記録 (runner が書く。forge_launches.jsonl とは cfg_fnv で対応)
ZTE_EFFECTIVE = "zero_thickness_edge_velocity_effective"     # 設計 DB の作動点要約で処置の識別に使うキー
_ZTE_TAG_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _stage_manifest():
    """solver_density_cuda/tools/stage_manifest.py (区間識別と処置のフィールドの読み出しの正本)。"""
    tools = str(FORGE_TOOLS)
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import stage_manifest
    return stage_manifest


def zte_spec_from_evaluate(ev: dict, known_tags=None) -> dict | None:
    """`evaluate.zero_thickness_edge_velocity` を正規化する。未指定・null・false → None (無効)。
    指定 → {"tags": 並べ替えたタグ, "rings": int}。known_tags を渡すとタグがその境界名に含まれるかも検査する。"""
    v = (ev or {}).get("zero_thickness_edge_velocity")
    if v is None or v is False:
        return None
    if not isinstance(v, dict):
        raise ValueError(f"evaluate.zero_thickness_edge_velocity は {{tags: [...], rings: N}} か null: {v!r}")
    unknown = set(v) - {"tags", "rings"}
    if unknown:
        raise ValueError(f"evaluate.zero_thickness_edge_velocity の未知キー: {sorted(unknown)} (tags | rings)")
    tags, rings = v.get("tags"), v.get("rings")
    # 道具 (mark_zero_thickness_edges.py) の契約に合わせる: タグは板の上下の壁面のちょうど 2 つ、rings は 1 以上
    if not isinstance(tags, (list, tuple)) or len(tags) != 2 or not all(isinstance(t, str) and _ZTE_TAG_RE.match(t) for t in tags):
        raise ValueError(f"evaluate.zero_thickness_edge_velocity.tags は板の上下の壁面の境界名 (英数・_.-) ちょうど 2 つ: {tags!r}")
    if len(set(tags)) != len(tags):
        raise ValueError(f"evaluate.zero_thickness_edge_velocity.tags に重複: {tags!r}")
    if isinstance(rings, bool) or not isinstance(rings, int) or rings < 1:
        raise ValueError(f"evaluate.zero_thickness_edge_velocity.rings は 1 以上の整数: {rings!r}")
    if known_tags is not None:
        bad = sorted(set(tags) - set(known_tags))
        if bad:
            raise ValueError(f"evaluate.zero_thickness_edge_velocity.tags に格子に無い境界名: {bad} (候補 {sorted(known_tags)})")
    return {"tags": sorted(tags), "rings": int(rings)}


def zte_spec(p: Problem, known_tags=None) -> dict | None:
    return zte_spec_from_evaluate(p.evaluate, known_tags)


def zte_signature_of(spec: dict | None) -> str:
    """処置の識別子。無効は "off"、有効は stage_manifest.zte_signature (タグ並べ替え・rings)。"""
    return "off" if spec is None else _stage_manifest().zte_signature(spec["tags"], spec["rings"])


def qc_cell_config(cfg: str, disc: str) -> str:
    """品質ゲート用の primal (cell) 変換の config。node 専用のキーを落とす (変換器も solverConfig を読んで検査する)。
    処置のキーも落とす: cell では起動時エラーで、w は最終の node 格子に書くのでこの段階には無い。"""
    return (cfg.replace(f'discretization: "{disc}"', 'discretization: "cell"')
            .replace(", nodeWallDirichlet: 1", "").replace(", nodeInletCornerWall: 1", "")
            .replace(ZTE_SPACE_FRAGMENT, ""))


def _mesh_nodes(h5path) -> int | None:
    with h5py.File(h5path, "r") as f:
        return int(f["VALUE/ro"].shape[0]) if "VALUE/ro" in f else None


def verify_zte_field(run_dir, spec: dict, mesh: str = MESH) -> dict:
    """meshFileName の `/AUX/<ZTE_FIELD>` を要求と照合する: 属性のタグ・rings が要求と一致、長さ = 節点数。
    戻り = stage_manifest.zte_identity (格子署名と w のハッシュは再計算した値)。"""
    h5 = Path(run_dir) / mesh
    a = _stage_manifest().zte_identity(h5, ZTE_FIELD)
    if a is None:
        raise ValueError(f"{h5}: 処置の重み /AUX/{ZTE_FIELD} が無い (前処理の道具が書いていない)")
    if a["tags"] != spec["tags"] or a["rings"] != spec["rings"]:
        raise ValueError(f"{h5}: /AUX/{ZTE_FIELD} の属性 tags={a['tags']} rings={a['rings']} が要求 {spec} と違う")
    n = _mesh_nodes(h5)
    if n is not None and a["n"] != n:
        raise ValueError(f"{h5}: w の長さ {a['n']} が節点数 {n} と違う (別の格子の w)")
    return a


def mark_zte_field(run_dir, spec: dict, mesh: str = MESH) -> dict:
    """前処理の道具で meshFileName に処置の重み w を書き、照合して識別量を返す。**道具の呼び出しはここ 1 か所**
    (道具の CLI が決まったらここだけ直す。試験はこの関数を差し替える)。格子を作るたびに prepare から呼ぶ。"""
    run_dir = Path(run_dir)
    if not ZTE_TOOL.exists():
        raise RuntimeError(f"処置の前処理の道具 {ZTE_TOOL} が無い (evaluate.zero_thickness_edge_velocity を使うには道具が要る)")
    cmd = [sys.executable, str(ZTE_TOOL), mesh, "--tags", *spec["tags"], "--bcond-config", "bcondConfig.yaml",
           "--rings", str(spec["rings"]), "--field", ZTE_FIELD]
    r = subprocess.run(cmd, cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "ZTE_MARK.txt").write_text(" ".join(cmd) + "\n" + (r.stdout or "") + (r.stderr or ""))
    if r.returncode != 0:
        raise RuntimeError(f"処置の前処理の道具が失敗 (rc={r.returncode})\n" + (r.stdout or "")[-1500:] + (r.stderr or "")[-1500:])
    return verify_zte_field(run_dir, spec, mesh)


def check_zte_stage(cfg_text: str, bcond_text: str, run_dir=None) -> dict | None:
    """処置を有効にした段の config が、ソルバの起動時エラーの条件 (node・SLAU/SLAU2・gpu 1・非周期・非軸対称、
    input/zeroThicknessEdge.cpp) に当たらないこと、
    run_dir を渡せば meshFileName に w があることを確かめ、**その段が読む w の識別量** (再計算) を返す。
    無効の段は None (何もしない)。違反は ValueError (driver では発散 [RuntimeError] でなく ERROR に分類される)。"""
    sm = _stage_manifest()
    zc = sm.zte_config(cfg_text)
    if zc is None:
        return None
    import yaml
    doc = yaml.safe_load(cfg_text) or {}
    ms = doc.get("mesh") or {}
    errs = []
    if zc["invalid"]:
        errs.append(f"space.zeroThicknessEdgeVelocity が不正 ({zc['invalid']})")
    if str(ms.get("discretization", "")).strip('"\'') != "node":
        errs.append(f"discretization {ms.get('discretization')!r} (node のみ)")
    if str(doc.get("solver", "")).strip('"\'') not in ("SLAU", "SLAU2"):
        errs.append(f"solver {doc.get('solver')!r} (SLAU / SLAU2 のみ)")
    if str(doc.get("gpu", "")).strip() != "1":
        errs.append(f"gpu {doc.get('gpu')!r} (gpu: 1 のみ)")
    if int(ms.get("isAxisymmetric", 0) or 0) or int((doc.get("physProp") or {}).get("isAxisymmetric", 0) or 0):
        errs.append("isAxisymmetric (軸対称は未検証で起動時エラー)")
    if re.search(r"\bkind\s*:\s*[\"']?periodic", bcond_text or ""):
        errs.append("周期境界 (初版は起動時エラー)")
    a = None
    if run_dir is not None and not errs:
        h5 = sm.zte_mesh_file(cfg_text, run_dir)
        try:
            a = sm.zte_identity(h5, zc["field"])
        except SystemExit as e:      # 道具・h5py が無い: 識別できないまま起動しない (driver では ERROR)
            raise ValueError(str(e)) from None
        if a is None:
            errs.append(f"meshFileName {h5} に /AUX/{zc['field']} が無い")
    if errs:
        raise ValueError("処置 (space.zeroThicknessEdgeVelocity) を有効にした段がソルバの起動条件に合わない: " + "; ".join(errs))
    return a


def _record_zte_launch(run_dir) -> None:
    """起動直前に、その段の config が処置を有効にしていれば条件を検査し、meshFileName から再計算した格子署名と
    w のハッシュを `zte_launches.jsonl` に追記する (無効の段は何も書かない = 旧 run とファイル構成も同じ)。"""
    rd = Path(run_dir)
    cfgp = rd / "solverConfig.yaml"
    cfg = cfgp.read_text() if cfgp.exists() else ""
    if "zeroThicknessEdgeVelocity" not in cfg:
        return
    bcp = rd / "bcondConfig.yaml"
    a = check_zte_stage(cfg, bcp.read_text() if bcp.exists() else "", rd)
    if a is None:
        return
    sm = _stage_manifest()
    rec = {"time": int(time.time()), "cfg_fnv": sm.fnv1a64(cfg), "field": sm.zte_config(cfg)["field"],
           "tags": a["tags"], "rings": a["rings"], "n": a["n"], "mesh_signature": a["mesh_signature"], "field_hash": a["field_hash"]}
    with open(rd / ZTE_LAUNCHES, "a") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")


def _last_jsonl(path):
    """jsonl の最後の非空行 (dict)。無い・壊れていれば None。前の行に遡らない。"""
    p = Path(path)
    if not p.exists():
        return None
    lines = [l for l in p.read_text().splitlines() if l.strip()]
    try:
        rec = json.loads(lines[-1]) if lines else None
    except Exception:
        return None
    return rec if isinstance(rec, dict) else None


_ZTE_LOG_EFF = re.compile(r"^'zeroThicknessEdgeVelocity' effective: (\d)(?: \(field (\S+), w<1 nodes (\d+), "
                          r"field_sha256 ([0-9a-f]+), mesh_signature ([0-9a-f]+)\))?")
_ZTE_LOG_FULL = re.compile(r"field_sha256[^(\n]*\(([0-9a-f]{64})\)[^(\n]*mesh_signature[^(\n]*\(([0-9a-f]{64})")
_ZTE_LOG_COUNTS = re.compile(r"(n_edge_nodes|n_marked_nodes|n_nodes)(?:\([A-Z0-9_]+\))?=(\d+)")


def _zte_from_log(run_dir) -> dict | None:
    """起動ログ forge_run.log (run_case.sh が起動ごとに上書き = 最後の起動のログ) から処置の記録を拾う。
    `'zeroThicknessEdgeVelocity' effective: 0|1 (field …, w<1 nodes N, field_sha256 <先頭16>, mesh_signature <先頭16>)` と、
    あれば全桁のハッシュと属性の数 (n_edge_nodes(E)= / n_marked_nodes(S)=)。行が無ければ None。"""
    p = Path(run_dir) / "forge_run.log"
    if not p.exists():
        return None
    info, after = None, 0
    with open(p, errors="replace") as fh:          # 長い run のログは大きいので、起動時の行の少し後で読むのをやめる
        for line in fh:
            m = _ZTE_LOG_EFF.match(line)
            if m:
                info = {"enabled": m.group(1) == "1", "source": "forge_run.log"}
                if m.group(2):
                    info.update({"field": m.group(2).rsplit("/", 1)[-1], "w_lt1_nodes": int(m.group(3)),
                                 "field_hash": m.group(4), "mesh_signature": m.group(5)})
                continue
            if info is not None:
                mf = _ZTE_LOG_FULL.search(line)
                if mf:
                    info["field_hash"], info["mesh_signature"] = mf.group(1), mf.group(2)
                for k, v in _ZTE_LOG_COUNTS.findall(line):
                    info[k] = int(v)
                after += 1
                if after >= 20:
                    break
    return info


def _last_launch_zte(run_dir):
    """最後の起動の処置の記録。forge_launches.jsonl の最後の非空行に処置のキー (入れ子 `zeroThicknessEdgeVelocity: {...}` /
    平坦 `zeroThicknessEdgeVelocity_<名>`) があればそれ、無ければ起動ログ (forge_run.log) の effective 行。
    戻り: (状態, 記録, cfg_fnv)。状態 = "absent" (起動記録が無い・壊れている) / "no_info" (処置の記録が無い) / "ok"。
    前の起動の値に遡らない (_last_launch_value と同じ厳格さ)。"""
    rec = _last_jsonl(Path(run_dir) / "forge_launches.jsonl")
    if rec is None:
        return "absent", None, None
    key = "zeroThicknessEdgeVelocity"
    v = rec.get(key)
    info = dict(v) if isinstance(v, dict) else ({} if v is None else {"enabled": v})
    for k, x in rec.items():
        if k.startswith(key) and k != key:
            info[k[len(key):].lstrip("_.")] = x
    if info:
        info["source"] = "forge_launches.jsonl"
        en = info.get("enabled")
        try:
            info["enabled"] = bool(int(en)) if en is not None else (True if info.get("field") else None)
        except (TypeError, ValueError):
            info["enabled"] = None
    else:
        info = _zte_from_log(run_dir)
    if not info:
        return "no_info", None, rec.get("cfg_fnv")
    return "ok", info, rec.get("cfg_fnv")


def _hash_agrees(recorded, full) -> bool:
    """ソルバの記録したハッシュ (全桁、または起動ログの先頭 16 桁) が再計算の全桁と一致するか。記録が無ければ True。"""
    if recorded in (None, ""):
        return True
    r = str(recorded).lower()
    return len(r) >= 16 and str(full).lower().startswith(r)


def zte_provenance(run_dir) -> dict:
    """評価の来歴: 処置の有効状態・フィールド・属性のタグ/rings・格子署名・w のハッシュ・節点数・起動の記録と、
    設計 DB の識別に使う `effective` ("off" / "tags=...;rings=N" / None = 不明)。
    有効の評価が確定する条件: 最後の起動 (forge_launches.jsonl) の config が solverConfig_main と同じ・runner の起動前の
    記録 (zte_launches.jsonl) がその起動に対応する・いま meshFileName から再計算した格子署名と w のハッシュが起動前と
    同じ・**ソルバが有効を記録** (forge_launches.jsonl に処置のキーがあればそれ、無ければ起動ログ forge_run.log の
    effective 行。格子署名・ハッシュを記録していれば一致)。旧バイナリは未知の space キーを
    黙って無視しうるので、設定だけでは処置が掛かったと言えない。無効は設定にキーが無ければ確定
    (最後の起動が有効を記録していたら矛盾で不明)。"""
    sm = _stage_manifest()
    rd = Path(run_dir)
    cfgp = (rd / "solverConfig_main.yaml") if (rd / "solverConfig_main.yaml").exists() else (rd / "solverConfig.yaml")
    st, launch, launch_fnv = _last_launch_zte(rd)
    pre = _last_jsonl(rd / ZTE_LAUNCHES)
    out = {"enabled": False, "field": None, "tags": None, "rings": None, "mesh_signature": None, "field_hash": None,
           "n_nodes": None, "attrs": None, "launch_state": st, "launch": launch, "pre_launch": pre, "effective": None, "reason": ""}
    if not cfgp.exists():
        out["reason"] = "solverConfig が無い"
        return out
    cfg = cfgp.read_text()
    zc = sm.zte_config(cfg)
    if zc is None:
        if launch is not None and launch.get("enabled"):
            out["reason"] = "設定は無効なのに最後の起動が有効を記録"
        else:
            out["effective"] = "off"
        return out
    out.update({"enabled": True, "field": zc["field"]})
    if zc["invalid"]:
        out["reason"] = f"設定が不正 ({zc['invalid']})"
        return out
    h5 = sm.zte_mesh_file(cfg, rd)
    try:
        a = sm.zte_identity(h5, zc["field"])
    except SystemExit as e:          # 道具・h5py が無い: 識別できないので不明 (collect・再判定を止めない)
        out["reason"] = str(e)
        return out
    if a is None:
        out["reason"] = f"meshFileName {h5} に /AUX/{zc['field']} が無い"
        return out
    out.update({"tags": a["tags"], "rings": a["rings"], "mesh_signature": a["mesh_signature"], "field_hash": a["field_hash"],
                "n_nodes": a["n"], "attrs": a["attrs"], "h5": Path(h5).name})
    n = _mesh_nodes(h5)
    fnv = sm.fnv1a64(cfg)
    if not a["tags"] or a["rings"] is None:
        out["reason"] = "w の属性 (tags / rings) が欠けている"
    elif n is not None and a["n"] != n:
        out["reason"] = f"w の長さ {a['n']} が節点数 {n} と違う"
    elif launch_fnv is None:
        out["reason"] = "ソルバの起動記録 (forge_launches.jsonl) が無い (どの config で起動したか結び付けられない)"
    elif launch_fnv != fnv:
        out["reason"] = "最後の起動の config が solverConfig_main と違う (本段で終わっていない)"
    elif pre is None or pre.get("cfg_fnv") != fnv:
        out["reason"] = f"runner の起動前の記録 ({ZTE_LAUNCHES}) が最後の段に対応していない"
    elif (pre.get("mesh_signature"), pre.get("field_hash")) != (a["mesh_signature"], a["field_hash"]):
        out["reason"] = "起動前と今で格子署名または w のハッシュが違う (起動後に h5 が変わった)"
    elif launch is None or launch.get("enabled") is not True:
        out["reason"] = f"最後の起動が処置の有効を記録していない (起動記録: {st})"
    elif launch.get("field") not in (None, zc["field"]):
        out["reason"] = f"最後の起動のフィールド {launch.get('field')} が設定 {zc['field']} と違う"
    elif not all(_hash_agrees(launch.get(k), a[k]) for k in ("mesh_signature", "field_hash")):
        out["reason"] = "ソルバが記録した格子署名または w のハッシュが再計算と違う"
    else:
        out["effective"] = sm.zte_signature(a["tags"], a["rings"])
    return out


def _dv(p: Problem, name, default=None) -> float:
    v = dv_value(p, name, default)
    return float(v["value"] if isinstance(v, dict) else v)


def design_snapshot(p: Problem) -> dict:
    """作動点で上書きされる**前**の設計点 (入口・外部流・ガス) を控える。逆設計はこれで固定する
    (plan §4.10: 形状は設計点で 1 つに決まる。作動点は CFD の BC/IC だけを変える)。"""
    return {"inflow": dict(p.spec["inflow"]), "external": dict(p.spec["external"]),
            "gamma": float(p.gamma), "cp": float(p.cp), "composition": p.raw.get("gas", {}).get("exhaust_composition")}


def select_operating_point(p: Problem, op: str | None) -> dict:
    """spec.operating_points[] から名前で 1 点を選び spec.external / inflow / ガスに反映する。
    op=None で operating_points が無ければ spec.external をそのまま使う。戻り値 = 選んだ点 (重み込み)。

    作動点は飛行条件 (external) だけでなく**燃焼器出口 (inflow) とガス (gamma/cp) も動く** —
    飛行 M が変わればインレット圧縮と燃焼加熱が変わり、powered と power-off では γ が 1.18 と 1.39 で違う
    (plan §4.10、NASA TM X-71972 TABLE 1 + CEA2)。"""
    ops = p.spec.get("operating_points")
    if not ops:
        if op not in (None, "", "default"):
            raise ValueError("spec.operating_points が無いのに --op が指定された")
        return {"name": "default", "weight": 1.0, "external": dict(p.spec["external"]),
                "inflow": dict(p.spec["inflow"]), "gas": {"gamma": p.gamma, "cp": p.cp}}
    names = [o["name"] for o in ops]
    if op in (None, "", "default"):
        op = names[0]
    if op not in names:
        raise ValueError(f"作動点 '{op}' が無い (候補: {names})")
    o = ops[names.index(op)]
    p.spec["external"] = dict(o["external"])
    if "inflow" in o:
        p.spec["inflow"] = {**p.spec["inflow"], **o["inflow"]}
    if "gas" in o:                      # 作動点ごとの γ / cp (powered 1.18 vs power-off 1.39) と排気組成 (R3, frozen_tp)
        unknown = set(o["gas"]) - {"gamma", "cp", "composition"}
        if unknown:
            raise ValueError(f"operating_points[{op}].gas の未知キー: {sorted(unknown)} (gamma | cp | composition)")
        gas = {k: float(v) for k, v in o["gas"].items() if k in ("gamma", "cp")}
        p.gamma = gas.get("gamma", p.gamma)
        p.cp = gas.get("cp", p.cp)
        p.raw.setdefault("gas", {}).update(gas)
        if "composition" in o["gas"]:   # モル分率 (CEA 凍結組成)。frozen_tp のときだけ効く
            if any(isinstance(k, bool) for k in o["gas"]["composition"]):
                raise ValueError(f"operating_points[{op}].gas.composition のキーに真偽値: 'NO' をクォートすること")
            p.raw["gas"]["exhaust_composition"] = {str(k): float(v) for k, v in o["gas"]["composition"].items()}
    return {"name": op, "weight": float(o.get("weight", 1.0)), "external": dict(p.spec["external"]),
            "inflow": dict(p.spec["inflow"]), "gas": {"gamma": p.gamma, "cp": p.cp, "composition": p.raw.get("gas", {}).get("exhaust_composition")}}


# --- R3: 凍結組成 TP (排気 = CEA 凍結組成の lump EXH, 外気 = 空気の lump AMB) ------------------------------------
# 外気 lump は 2026-09-30 に AIR → AMB へ改名 (ソルバ内蔵の擬似種 AIR と衝突するため。plan tooling-nozzle-sern-chain R8)。順序は不変
SPECIES_ORDER = ("EXH", "AMB")     # 既定 (evaluate.tp_species 省略時の別名 = lumps {EXH: stream inflow, AMB: stream external})


def frozen_gases(p: Problem) -> dict | None:
    """gas.model: frozen_tp のとき {"exhaust": FrozenGas, "ext": FrozenGas, "href_T", "layout": SpeciesLayout, "db",
    "transported": [FrozenGas (輸送種ごと)]}。cpg なら None。
    排気組成は作動点 `gas.composition` (select_operating_point が `gas.exhaust_composition` に写す) のモル分率、
    外気は `spec.external.composition` (モル分率) があればそれ、無ければ乾燥空気。
    輸送種の配置は統一スキーマ `evaluate.tp_species` (省略時 = `[EXH, AMB]` = 流れごとの lump; `{mode: full}` で実種の和集合、
    このとき排気率は受動スカラ `roXi` で輸送) を流れ {inflow, external} で解決する (plan cea-mole-fraction §4.5)。"""
    if not p.is_frozen_tp:
        return None
    from ..gas.composition import parse_tp_species, resolve_species_layout
    href = float(p.raw["gas"].get("thermo_href_temp", 298.15))
    comp = p.raw["gas"].get("exhaust_composition")
    if not comp:
        raise ValueError("gas.model: frozen_tp には gas.exhaust_composition か operating_points[].gas.composition (モル分率) が要る")
    db = p.species_db
    ext_comp = p.spec["external"].get("composition")
    from ..gas.frozen import AIR_MOLE
    exh = FrozenGas.from_mole(comp, "EXH", href, db)
    ext = FrozenGas.from_mole(ext_comp if ext_comp else AIR_MOLE, "AMB", href, db)
    ev = dict(p.evaluate)
    ev.setdefault("tp_species", ["EXH", "AMB"])
    layout = resolve_species_layout(parse_tp_species(ev), {"inflow": exh.Y, "external": ext.Y}, db,
                                    condensing_species=None, condensation=False)
    transported = []
    for sname in layout.species:
        e = layout.entries[sname]
        transported.append(FrozenGas(e.lump_mass if e.lump_mass else {sname: 1.0}, sname, href, db))
    return {"exhaust": exh, "ext": ext, "href_T": href, "layout": layout, "db": db, "transported": transported}


def frozen_transport(p: Problem, layout) -> dict:
    """frozen_tp の `gas.transport` をこの作動点の実種に絞って照合する ({実種: モデル}、実種の順)。
    SERN は作動点で排気の構成種が変わる (m4_off は燃料なしで N2/O2/AR/CO2 だけ) ので、YAML には全作動点の実種の和集合を書き、
    ここで作動点に無い種を落とす。作動点の実種の書き漏れはエラー (resolve_transport, required=True)。2026-10-01 (R8)。"""
    from ..gas.composition import resolve_transport, transport_real_species
    reals = set(transport_real_species(layout))
    sub = {k: v for k, v in (p.gas_transport or {}).items() if k in reals}
    return resolve_transport(layout, sub, required=True)


def write_species_db(p: Problem, run_dir, gases: dict | None) -> None:
    """`species_meta.yaml` を run dir に書く (cpg なら何も書かない)。lump の熱物性はソルバが起動時に合成するので
    合成済み擬似種の `species_db.yaml` は書かない (2026-09-30、plan thermophysics-solver-owned-species-db §5.2 / SERN R8)。
    ソルバ内蔵で解決できない実種 (外部 DB 由来など) があるときだけ、その生エントリを `species_db_external.yaml` に書く。"""
    if gases is None:
        return
    from ..gas.composition import solver_species_config, species_db_raw_yaml, write_species_meta
    _, external = solver_species_config(gases["layout"])
    tr = frozen_transport(p, gases["layout"]) if p.raw.get("gas", {}).get("transport") is not None else None
    write_species_meta(gases["layout"], run_dir, tr)
    if external:
        (Path(run_dir) / "species_db_external.yaml").write_text(species_db_raw_yaml(external))


def gas_states(p: Problem) -> dict:
    """入口 (燃焼器出口) と外気の一様状態。cpg: 単一 (γ, R) を両方に使う (旧; codex C2 が指摘した外部動圧 −15 % の原因)。
    frozen_tp (R3): 排気は CEA 凍結組成の擬似種、外気は空気で、ρ = P/(R T)・u = M a(T) をそれぞれの NASA-9 物性で計算する。"""
    g, cp = p.gamma, p.cp
    R = cp * (g - 1.0) / g
    fi, ex = p.spec["inflow"], p.spec["external"]
    if fi.get("mode", "supersonic") != "supersonic":
        raise ValueError("inflow.mode は現状 supersonic のみ (sonic_throat は S1 接続が未実装)")
    gases = frozen_gases(p)

    def turb(ro, u):
        k = 1.5 * (0.01 * u) ** 2
        return {"k": float(k), "omega": float(ro * k / (1.8e-5 * 10.0))}

    def st(M, P, T, gas: FrozenGas | None):
        if gas is None:
            ro = P / (R * T); u = M * np.sqrt(g * R * T)
            out = {"M": float(M), "P": float(P), "T": float(T), "ro": float(ro), "u": float(u), "R": R, "gamma_T": g, "gas": "cpg"}
        else:
            out = gas.state(M, P, T)
        out.update(turb(out["ro"], out["u"]))
        return out
    st_ex = st(fi["M_in"], fi["p_in"], fi["T_in"], gases["exhaust"] if gases else None)
    st_en = st(ex["M_inf"], ex["p_inf"], ex["T_inf"], gases["ext"] if gases else None)
    out = {"exhaust": st_ex, "ext": st_en, "R": R, "gas_model": "frozen_tp" if gases else "cpg",
           "q_inf": 0.5 * st_en["ro"] * st_en["u"] ** 2}
    if gases:
        L = gases["layout"]
        out["species"] = list(L.species); out["href_T"] = gases["href_T"]; out["tp_mode"] = L.mode; out["tracer"] = bool(L.tracer)
        out["exhaust"]["Y"] = L.Y_transport("inflow"); out["ext"]["Y"] = L.Y_transport("external")
        if L.tracer:                     # full: 排気率 ξ (排気入口 1 / 外気入口 0) を受動スカラ roXi で輸送
            out["exhaust"]["Xi"] = 1.0; out["ext"]["Xi"] = 0.0
        out["gas_summary"] = {"exhaust": gases["exhaust"].summary(), "ext": gases["ext"].summary()}
    return out


def ideal_thrust(p: Problem, st: dict) -> tuple:
    """理想総推力 F/(p_in H) と出口 M。frozen_tp は同じ NASA-9 物性の等エントロピー膨張 (R3)、cpg は解析式。"""
    ex, en = st["exhaust"], st["ext"]
    gases = frozen_gases(p)
    if gases:
        F_nd, M_e, _ = ideal_gross_thrust_frozen(gases["exhaust"], ex["M"], ex["T"], ex["P"], en["P"])
        return F_nd, M_e
    return ideal_gross_thrust(ex["M"], en["P"] / ex["P"], p.gamma)


def design_from_problem(p: Problem, design: dict | None = None):
    """逆設計は**設計点**で行う (`design` = 作動点適用前の `design_snapshot(p)`)。作動点 (operating_points) は
    CFD の境界条件・IC だけを変え、形状は変えない (2026-09-05 修正: それまで作動点ごとに p_ext が変わり kernel/形状が
    作動点依存になっていた — run_0010/0017 は作動点間で形状が一致していない)。
    2026-09-05 追補 (§4.10): 作動点が inflow / gas も動かすようになったので、入口状態と γ も設計点で固定する。"""
    geo = p.geometry
    d0 = design or design_snapshot(p)
    fi_d, ext_d, g_d = d0["inflow"], d0["external"], float(d0["gamma"])
    p_ext_ratio = float(ext_d["p_inf"]) / float(fi_d["p_in"])
    spec = SernKernelSpec(M_in=float(fi_d["M_in"]),
                          theta_r0=np.deg2rad(_dv(p, "theta_r0_deg")), theta_c0=np.deg2rad(_dv(p, "theta_c0_deg")),
                          L_cowl=_dv(p, "L_cowl"), gamma=g_d, p_ext_over_p_in=p_ext_ratio,
                          x_max=float(geo.get("x_max_kernel", 10.0)), nj=int(geo.get("nj_moc", 301)),
                          dx=float(geo.get("dx_moc", 2e-3)))
    if str(geo.get("mode", "keypoint")) == "straight":
        k = PlanarMOC(spec).march()
        d = k.straight_design(_dv(p, "L_ramp"))
    else:
        k = PlanarMOC(spec).march(stop_at=(_dv(p, "f"), _dv(p, "M_c")))   # c の少し先で打ち切る
        d = k.design_ramp(M_c=_dv(p, "M_c"), f=_dv(p, "f"))
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    fr = wall_forces(d, spec.M_in, g_d, pa_over_pin=p_ext_ratio, x_ref=float(xr), y_ref=float(yr))
    theta_b = float(k.TH[-1, 0])   # 自由境界の終端角 (せん断層に格子線を沿わせる)
    return k, d, fr, theta_b


def _solver_config(p: Problem, nsteps: int, out_int: int, cfl: float, p_ref: float) -> str:
    """`evaluate.implicit_relax` / `evaluate.p_min` を deltaT / space に挿入できる (2026-09-06)。

    SERN は元から `blockDPLUR: 1` (陰解法) だが `implicitRelax` を設定しておらず既定 1.0 だった。
    風洞チェーンは同じ陰解法に対し「cfl 6 + implicitRelax 0.7 が生産推奨」(case/45 run_0018) と配管済みで、
    メモリ [[implicit-cfl-ceiling-eos-floor]] も「NS 陰解法の上限は P 床洗浄律速、効くのは implicitRelax のみ」
    と記録している。run_0075 の発散 (M∞10 の boat-tail 膨張で 18 % のノードが pMin=1 Pa に着地 → 負密度) は
    まさにこの指紋なので、relax を効かせる。"""
    _lim = int(p.evaluate.get("limiter", 2))   # 2=Venkatakrishnan (既定), 1=Barth
    # リミッタの試行値を流束が適用する増分と同じ形・同じ点で評価する (既定 0 = 従来)。
    # plan convection-node-wall-reconstruction §4.8。Barth と組むと厳密有界になる
    # `limiter_match_recon` は廃止 (plan limiter-config-simplify §4.2)。`limiter_scaled` に内包した。
    if "limiter_match_recon" in p.evaluate:
        raise ValueError("evaluate.limiter_match_recon は廃止。limiter_scaled: 0 (旧経路) / 1 (評価点一致 + 無次元化) を使うこと "
                         "(中間の match_recon=1, scaled=0 は機能打ち切り)")
    _lsc = int(p.evaluate.get("limiter_scaled", 1))   # 既定 1 = 修正版 (2026-09-20)
    if _lsc not in (0, 1):
        raise ValueError(f"evaluate.limiter_scaled は 0 か 1 (比の形 2 は棄却済み): {_lsc}")
    # 既定はソルバと揃える: 修正版 (1) は 0.05、旧経路 (0) は 1.0 (旧経路の K は device 側で 1.f 固定)
    _vk = float(p.evaluate.get("venkat_k", 0.05 if _lsc == 1 else 1.0))
    # space.slauWallNormalChi (2026-09-26 既定化、plan convection-slau-wall-normal-chi-default §4.4): runner は既定でキーを書かない
    # (= auto。node+Dirichlet 壁+SLAU なら 1、品質検査用の cell 変換では 0 に静かに解決する)。明示 1 を書くと cell 変換が起動エラーになる。
    # 問題 YAML `mesh.slau_wall_normal_chi: 0` のときだけ明示 0 (旧挙動) を書く。
    _wnc = p.mesh.get("slau_wall_normal_chi", None)
    if _wnc is not None and int(_wnc) not in (0, 1):
        raise ValueError(f"mesh.slau_wall_normal_chi must be 0 or 1 (or omitted for auto): {_wnc}")
    _wnc_key = ", slauWallNormalChi: 0" if (_wnc is not None and int(_wnc) == 0) else ""
    # `evaluate.limiter_ref: {length, ro, p, a}` (2026-09-29): リミッタ基準値の明示固定。既定 (キーなし) は forge が
    # 領域の対角長と初期場の平均から自動で決めるので、**領域の大きさを変える比較 (幅系列) では離散化そのものが変わる**
    # (plan boundary-node-farfield-characteristic §5.1 #3 の V2a で確認)。幅系列では全幅で同じ値を書く。
    _lr = p.evaluate.get("limiter_ref")
    if _lr is not None:
        _wnc_key += (f", limiterRefLength: {float(_lr['length'])!r}, limiterRoRef: {float(_lr['ro'])!r}, "
                     f"limiterPRef: {float(_lr['p'])!r}, limiterARef: {float(_lr['a'])!r}")
    ir = p.evaluate.get("implicit_relax")
    _relax = f", implicitRelax: {float(ir)}" if ir is not None else ""
    pm = p.evaluate.get("p_min")
    _pmin = f", pMin: {float(pm)}" if pm is not None else ""   # physProp 配下 (space ではない)
    disc = p.mesh.get("discretization", "cell")
    model = p.evaluate.get("model", "euler")
    # 厚さ 0 の板の自由端の速度再構成 (plan convection-zero-thickness-edge-reconstruction §4): 有効/フィールド名だけを出す。
    # 段階起動の全段に同じ値が出る (run_staged の各段は cfg_main の置換なので)。1 次・層流暖機の段では効かないが、値は
    # 同じに出して区間識別と来歴をそろえる。node 以外はソルバの起動時エラーなので、ここで先に止める
    _zte_key = ""
    if zte_spec(p) is not None:
        if disc != "node":
            raise ValueError(f"evaluate.zero_thickness_edge_velocity は node のみ (mesh.discretization: {disc})")
        _zte_key = ZTE_SPACE_FRAGMENT
    node_keys = ", nodeWallDirichlet: 1" if (disc == "node" and model != "euler") else ""
    # R4b(i) (2026-09-13): 入口∩壁の角ノードの半割面所有を壁側に (converter が変換時に読む)。既定 0 = 旧 run とビット一致。
    # 生産 YAML は mesh.node_inlet_corner_wall: 1 (角ノードの壁圧 1.75 p_in 対策; plans/active/boundary-node-inlet-corner-wall.md)
    if disc == "node" and int(p.mesh.get("node_inlet_corner_wall", 0)):
        node_keys += ", nodeInletCornerWall: 1"
    # R3 (frozen_tp): 排気 EXH / 外気 AMB の 2 lump TP。thermoHrefTemp (sensible datum) は陰解法の χ_eos 桁違い対策で必須
    # ([[isobutane-wt-semiperfect]] / runner_axismach と同じ)。IC の roe も同じ基準で組む (paste_region_ic)
    # physProp.species は lump 記法 ({name, lump: {構成種: モル分率}, basis: mole}) で、NASA-9 はソルバが起動時に合成する (R8)
    if p.is_frozen_tp:
        from ..gas.composition import physprop_species_flow, solver_species_config
        gases = frozen_gases(p); L = gases["layout"]
        items, external = solver_species_config(L)
        _db = ', speciesDBFile: "species_db_external.yaml"' if external else ""
        _tp = f", species: {physprop_species_flow(items)}{_db}, thermoHrefTemp: {gases['href_T']}"
        if L.tracer:
            _tp += ", tracer: exhaust"
        _tm = 2
    else:
        _tp = ""; _tm = 0
    if model == "euler":
        phys = f"physProp: {{thermalMethod: {_tm}, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}{_pmin}{_tp}}}"
        turb = 'turbulence: {model: "none"}'
    else:
        # 輸送物性 (R8 段 (ii)、2026-09-30): frozen_tp で problem に `gas.transport` ({実種: モデル}) があるときだけ、
        # 実種ごとの輸送物性 (viscMethod 2 + physProp.transport、混合則 CEA frozen) にする。無ければ従来どおり
        # viscMethod 1 (空気の Sutherland) — 既存の run と同じ config。thermCondMethod は viscMethod 2 では読まれないので落とす。
        # visc (dt と陰解法対角の剛性見積り)・thermCond (必須キー)・prandtlLam (SST 壁関数の回復係数) は残す (runner_axismach と同じ)
        _transport = None
        if p.is_frozen_tp and p.raw.get("gas", {}).get("transport") is not None:
            from ..gas.composition import physprop_transport_flow
            _transport = frozen_transport(p, frozen_gases(p)["layout"])
        if _transport is not None:
            phys = (f"physProp: {{thermalMethod: {_tm}, viscMethod: 2, visc: 1.8e-5, thermCond: 0.0257, "
                    f"prandtlLam: 0.72, cp: {p.cp}, gamma: {p.gamma}{_pmin}{_tp}, transport: {physprop_transport_flow(_transport)}}}")
        else:
            phys = (f"physProp: {{thermalMethod: {_tm}, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, "
                    f"thermCondMethod: 1, prandtlLam: 0.72, cp: {p.cp}, gamma: {p.gamma}{_pmin}{_tp}}}")
        # 壁処理は**既定 0 (低 Re 壁解像)**。node の SST 壁関数は使わない方針 (2026-09-20)。
        # 壁関数を使うには問題 YAML に `evaluate.wall_treatment_sst: 1` を明示し、理由を run の README に書くこと。
        _wts = int(p.evaluate.get("wall_treatment_sst", 0))
        turb = ('turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, '
                f'katoLaunder: 1, wallTreatmentSST: {_wts}}}')
    return f"""mesh: {{discretization: "{disc}", isAxisymmetric: 0{node_keys}, meshFileName: "{MESH}", valueFileName: "{MESH}"}}
gpu: 1
solver: "SLAU"
{phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1{_relax}, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: {_lim}, pRef: {p_ref}, limiterScaled: {_lsc}, venkatK: {_vk}{_wnc_key}{_zte_key}}}
{turb}
initial: "uniform_p101325_u10"
"""


def _bcond_config(p: Problem, st: dict) -> str:
    model = p.evaluate.get("model", "euler")
    ex, en = st["exhaust"], st["ext"]
    # `evaluate.outlet_kind`: outflow (既定・全量外挿) / statPress。**既定は 2026-09-23 に statPress から変更**。
    # SERN の出口と bottom は設計上つねに超音速なので、静圧指定は node の壁列・後流の**亜音速ノード**に
    # Ps ≪ 実出口圧を課し、そこから圧力が育つ (procedures/recommended-settings.md「出口」/ [[node-supersonic-exit-outflow]])。
    # 実績: run_0121 で出口 P 7.5 → 128 kPa で発散、run_0430–0436 で出口の亜音速率 3.5 → 99.5 %・far_bottom 22 MPa。
    okind = str(p.evaluate.get("outlet_kind", "outflow"))

    def inlet(name, pid, s):
        return (f"{name}: {{physID: {pid}, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , "
                f"floats: {{ro: {s['ro']:.6g}, Ux: {s['u']:.6g}, Uy: 0.0, Uz: 0.0, Ps: {s['P']:.6g}, k: {s['k']:.6g}, omega: {s['omega']:.6g}{inlet_species_floats(s)}}}}}\n")

    def outlet(name, pid):
        if okind == "outflow":
            return f"{name}: {{physID: {pid}, kind: outflow, outputHDFflg: 0, ints: , floats: }}\n"
        return (f"{name}: {{physID: {pid}, kind: outlet_statPress, outputHDFflg: 0, ints: , "
                f"floats: {{Ps: {en['P']:.6g}, Pt: {en['P']:.6g}, Tt: {en['T']:.6g}}}}}\n")

    def wall(name, pid):
        # 壁行は spec.wall_thermal が単一ソース (断熱 wall / 等温 wall_isothermal+Ts)。Euler は slip
        return f"{name}: {p.wall_bcond_line(model == 'euler', phys_id=pid, output=1)}\n"
    P = PHYS_SERN
    return (inlet("inlet_nozzle", P["inlet_nozzle"], ex) + inlet("inlet_ext", P["inlet_ext"], en)
            + outlet("outlet", P["outlet"]) + wall("ramp", P["ramp"]) + wall("cowl_in", P["cowl_in"])
            + wall("cowl_out", P["cowl_out"]) + outlet("bottom", P["bottom"])
            + (outlet("top_out", P["top_out"]) if p.evaluate.get("top_out_kind", "outlet") == "outlet"
               else f"top_out: {{physID: {P['top_out']}, kind: slip, outputHDFflg: 0, ints: , floats: }}\n")
            # 機体上面 + base (§4.11)。機体の力なので帳簿外だが base 圧の診断のため壁出力する。
            # 既定は slip (2D 中立モデル)。**`evaluate.vehicle_kind: wall` で等温粘性壁**にできる
            # (3D の生産仕様 R4f/R4e と揃えるため。有限ベース `mesh.t_base > 0` の診断で使う —
            #  slip の base は鋭い 90° 角で wall_dist が 0 に落ち、SST の ω が発散する)
            + ((f"vehicle: {p.wall_bcond_line(model == 'euler', phys_id=P['vehicle'], output=1)}\n"
                if str(p.evaluate.get("vehicle_kind", "slip")) == "wall"
                else f"vehicle: {{physID: {P['vehicle']}, kind: slip, outputHDFflg: 1, ints: , floats: }}\n")
               if int(p.mesh.get("ext_top", 0)) else ""))


def inlet_species_floats(s: dict) -> str:
    """frozen_tp の入口組成 (Y0 = 排気, Y1 = 空気)。cpg (Y 無し) は空文字。"""
    Y = s.get("Y")
    if Y is None:
        return ""
    txt = "".join(f", Y{i}: {float(y):.6g}" for i, y in enumerate(Y))
    if "Xi" in s:
        txt += f", Xi: {float(s['Xi']):.6g}"
    return txt


def region_ic_arrays(upper, st: dict, gamma: float) -> dict:
    """領域マスク upper (排気側) から保存量 IC を作る。cpg: roe = P/(γ−1) + ½ρu²。
    frozen_tp: roe = ρ (h_sens(T) − R T) + ½ρu² を各領域のガスで (forge の thermoHrefTemp 基準と同一)、roY0/roY1 = 領域組成。"""
    ex, en = st["exhaust"], st["ext"]
    ro = np.where(upper, ex["ro"], en["ro"]); u = np.where(upper, ex["u"], en["u"]); P = np.where(upper, ex["P"], en["P"])
    out = {"ro": ro, "roUx": ro * u, "roUy": np.zeros_like(ro), "roUz": np.zeros_like(ro),
           "roK": ro * np.where(upper, ex["k"], en["k"]), "roOmega": ro * np.where(upper, ex["omega"], en["omega"])}
    if st.get("gas_model") == "frozen_tp":
        e_ex = ex["h_sens"] - ex["R"] * ex["T"]; e_en = en["h_sens"] - en["R"] * en["T"]
        out["roe"] = ro * (np.where(upper, e_ex, e_en) + 0.5 * u * u)
        for i, (ye, ya) in enumerate(zip(ex["Y"], en["Y"])):      # 輸送種ごとの領域組成 (lumped: [1,0]/[0,1], full: 実種ベクトル)
            out[f"roY{i}"] = ro * np.where(upper, float(ye), float(ya))
        if st.get("tracer"):
            out["roXi"] = np.where(upper, ro, 0.0)
    else:
        out["roe"] = P / (gamma - 1.0) + 0.5 * ro * u * u
    return out


def write_ic_arrays(v, arrays: dict) -> None:
    """VALUE グループへ IC を書く。roY{s} は無ければ作る (converter は化学種を知らない; forge は VALUE/roY を優先して読む)。"""
    for k, a in arrays.items():
        if k in ("roK", "roOmega") and k not in v:
            continue
        if k in v:
            v[k][:] = a
        else:
            v.create_dataset(k, data=np.asarray(a, dtype=np.float32))


def apply_wall_offset(design, wall_offset: dict, H: float):
    """壁を法線方向 (流体と反対側) に dn(x) [m] だけ動かした SernDesign を返す (無次元化して適用)。
    ramp: 上壁なので +n = 上、cowl: 下壁 (内面が上向き) なので +n = 下。表は (x_m, dn_m) の 2 列。"""
    import copy
    d = copy.deepcopy(design)
    for name, arr in (("ramp", d.ramp_xy), ("cowl", d.cowl_xy)):
        tbl = wall_offset.get(name)
        if tbl is None:
            continue
        tbl = np.asarray(tbl, dtype=float)
        dn = np.interp(arr[:, 0], tbl[:, 0] / H, tbl[:, 1] / H, left=tbl[0, 1] / H, right=tbl[-1, 1] / H)
        t = np.gradient(arr, axis=0); t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-30)[:, None]
        nrm = np.column_stack([-t[:, 1], t[:, 0]]) if name == "ramp" else np.column_stack([t[:, 1], -t[:, 0]])
        arr += dn[:, None] * nrm
    return d



def moc_ic_arrays(kern, xn, yn, upper, st: dict, gamma: float, gas=None) -> tuple:
    """**MOC 場を初期値にする** (2026-09-19, ユーザ提案)。

    現行の領域別一様 IC は、ノズル内を燃焼器出口状態 (例 101 kPa) で埋める。実際の解は出口で
    ~6.5 kPa まで膨張するので、**初期値が 17 倍ずれた状態**から始めることになり、
    梯子 12000 step を cfl 0.1 で這わせる主因になっている。MOC は同じ形状の非粘性解を
    station ごとに (Y, TH, M) で持っているので、それを内挿すれば初期値が解のすぐ近くから始まる。

    返り値: (arrays, n_moc) — n_moc は MOC を当てられたノード数 (残りは一様 IC のまま)。
    MOC の被覆外 (kernel の x 範囲外・上下境界の外) は `upper` による一様値に落とす。
    等エントロピー: よどみ量は排気の入口状態から作り、M(x,y)・θ(x,y) で静圧・静温・速度に展開する。
    """
    ex = st["exhaust"]
    g = float(gamma); gm = g - 1.0
    R = float(ex.get("R", ex["P"] / (ex["ro"] * ex["T"])))
    M_in = float(ex.get("M", 0.0))
    X = np.asarray(kern.X)
    M = np.full(len(xn), np.nan); TH = np.full(len(xn), np.nan)
    # **被覆外は外挿する** (2026-09-19)。x を kernel 範囲に、y を各 station の範囲にクランプして端の値を伸ばす。
    # 落とすと排気域の中に一様値 (入口状態) の塊が残り、その境界が 17 倍の圧力段差になって
    # 一様 IC より悪い初期値になる (実測: 排気域 32268 ノードのうち 13111 が一様のまま → mid 段 step 6 で発散)
    sel = np.flatnonzero(upper)
    xc = np.clip(xn[sel], X[0], X[-1])
    idx = np.clip(np.searchsorted(X, xc), 1, len(X) - 1)
    for j, (i1, x, y) in enumerate(zip(idx, xc, yn[sel])):
        i0 = i1 - 1
        t = (x - X[i0]) / max(X[i1] - X[i0], 1e-30)
        a = np.interp(y, kern.Y[i0], kern.M[i0]); b = np.interp(y, kern.Y[i1], kern.M[i1])   # 端はクランプ = 外挿
        ta = np.interp(y, kern.Y[i0], kern.TH[i0]); tb = np.interp(y, kern.Y[i1], kern.TH[i1])
        M[sel[j]] = a + t * (b - a)
        TH[sel[j]] = ta + t * (tb - ta)
    ok = np.isfinite(M) & (M > 0.0)
    base = region_ic_arrays(upper, st, gamma)
    if not ok.any():
        return base, 0
    # **NASA-9 に整合な等エントロピー展開** (2026-09-19, codex plan レビュー M4)。
    # 旧実装は一定 γ の式で T/P/q を作り `roe` だけ NASA-9 に置換していたため、入口に対して
    # 全エンタルピーが +1.20 %・エントロピーが +30.7 J/(kg·K) ずれていた。
    # ここでは同じ NASA-9 物性で次の 2 式を解く:
    #   h_sens(T) + ½ M² γ(T) R T = h0_in     (全エンタルピー保存)
    #   p = p_in · exp[(s°(T) − s°(T_in)) / R]  (等エントロピー)
    # 速度は q = M · a(T)。cpg のときは従来どおり一定 γ の式。
    Mo = M[ok]
    if st.get("gas_model") == "frozen_tp":
        if gas is None:
            raise ValueError("frozen_tp の MOC IC には FrozenGas が要る")
        T_in = float(ex["T"]); P_in = float(ex["P"])
        h0_in = float(np.ravel(gas.h_sens(T_in))[0]) + 0.5 * float(ex["u"]) ** 2
        T = np.full_like(Mo, T_in)
        for _ in range(40):                      # h0 一定から T を Newton で解く (γ(T) も更新)
            gT = np.asarray(gas.gamma(T)); hT = np.asarray(gas.h_sens(T))
            F = hT + 0.5 * Mo ** 2 * gT * R * T - h0_in
            cpT = np.asarray(gas.cp_mass(T))
            dF = cpT + 0.5 * Mo ** 2 * gT * R      # γ の T 依存は 2 次なので無視 (収束には十分)
            step = F / np.maximum(dF, 1e-30)
            step = np.clip(step, -0.3 * T, 0.3 * T)
            T = np.maximum(T - step, 1.0)
            if np.max(np.abs(step)) < 1e-8 * np.max(T):
                break
        s0 = np.asarray(gas.s0_mass(T)); s0_in = float(np.ravel(gas.s0_mass(T_in))[0])
        P = P_in * np.exp((s0 - s0_in) / R)
        ro = P / (R * T)
        q = Mo * np.asarray(gas.a(T))
    else:
        T0 = float(ex["T"]) * (1.0 + 0.5 * gm * M_in * M_in)
        P0 = float(ex["P"]) * (1.0 + 0.5 * gm * M_in * M_in) ** (g / gm)
        f = 1.0 + 0.5 * gm * Mo ** 2
        T = T0 / f; P = P0 / f ** (g / gm); ro = P / (R * T)
        q = Mo * np.sqrt(g * R * T)
    base["ro"][ok] = ro
    base["roUx"][ok] = ro * q * np.cos(TH[ok])
    base["roUy"][ok] = ro * q * np.sin(TH[ok])
    base["roK"][ok] = ro * float(ex["k"]); base["roOmega"][ok] = ro * float(ex["omega"])
    if st.get("gas_model") == "frozen_tp":
        # 内部エネルギーは **NASA-9 をそのまま使う** (定 cv の近似はしない)。forge の thermoHrefTemp 基準と同一
        if gas is None:
            raise ValueError("frozen_tp の MOC IC には FrozenGas が要る")
        base["roe"][ok] = ro * (gas.h_sens(T) - ex["R"] * T + 0.5 * q * q)
        if st.get("tracer"):
            base["roXi"][ok] = ro
        for i, ye in enumerate(ex["Y"]):
            base[f"roY{i}"][ok] = ro * float(ye)
    else:
        base["roe"][ok] = P / gm + 0.5 * ro * q * q
    return base, int(ok.sum())


def paste_region_ic(h5path, y_mid, y_top, scale: float, st: dict, gamma: float, kern=None, gas=None) -> int:
    """領域別一様 IC: 中間線とランプ/プルーム上線の間 = 燃焼器出口状態、それ以外 (カウル下・ランプ側外部流) = 外部流。
    `kern` を渡すと排気側を **MOC 場**で埋める (`moc_ic_arrays`)。戻り値 = MOC を当てたノード数。"""
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn, yn = cc[:, 0] / scale, cc[:, 1] / scale
        upper = (yn > y_mid(xn)) & (yn < y_top(xn))
        if kern is None:
            write_ic_arrays(f["/VALUE"], region_ic_arrays(upper, st, gamma)); return 0
        arrays, n = moc_ic_arrays(kern, xn, yn, upper, st, gamma, gas)
        write_ic_arrays(f["/VALUE"], arrays); return n


def stamp_region_ic_species(h5path, run_dir, st: dict, gases: dict | None) -> str | None:
    """新規初期場 (paste_region_ic) に化学種の属性を付ける (frozen_tp のみ; plan thermophysics-solver-owned-species-db §4.3 #3b)。
    IC が roe を作った排気・外気のガス (datum・輸送種組成・R・e_sens) を宛先の `forge --resolve-species` の記録と照合し、
    一致したときだけ付ける (違えば forge_species.SpeciesCheckError)。cpg は何もしない。"""
    if gases is None or st.get("gas_model") != "frozen_tp":
        return None
    L = gases["layout"]
    species = list(L.species)
    mixes = [("exhaust", st["exhaust"]["Y"], gases["exhaust"].R, gases["exhaust"].e_sens),
             ("external", st["ext"]["Y"], gases["ext"].R, gases["ext"].e_sens)]
    return _forge_species().stamp_new_field(h5path, run_dir, species, [float(L.entries[k].MW) for k in species],
                                            gases["href_T"], mixes, tool="paste_region_ic")


def convert_mesh(run_dir, msh: str, out: str) -> None:
    """gmsh msh → forge h5。**exit code で判定しない**: 一部の環境 (AWS g5 / CUDA 13) で converter は
    h5 を書き切ってから終了時に `GPUassert: invalid argument` を出して非零で抜ける (既知・無害)。
    成否は出力ファイルの存在とサイズで見る。"""
    r = subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), msh, out], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    f = Path(run_dir) / out
    if not f.exists() or f.stat().st_size < 1024:
        raise RuntimeError(f"convertGmshToForge が {out} を作れなかった (rc={r.returncode})\n"
                           + (r.stdout or "")[-1500:] + (r.stderr or "")[-1500:])


def prepare(problem_path, run_dir, nsteps=None, op: str | None = None, wall_offset=None) -> dict:
    """op: 作動点名 (spec.operating_points)。wall_offset: {"ramp": (x_m, dn_m), "cowl": (x_m, dn_m)} の
    法線オフセット表 [m] (S5 δ* 一発補正。壁を流体と反対側へ dn だけ動かす)。"""
    p = load_problem(problem_path)
    if p.type != "sern_2d":
        raise ValueError("runner_sern は sern_2d 専用")
    # 厚さ 0 の板の自由端の処置: 指定の検査 (タグは格子の境界名) を格子を作る前に済ませる
    zte = zte_spec(p, known_tags=PHYS_SERN.keys())
    if zte is not None and p.mesh.get("discretization", "cell") != "node":
        raise ValueError("evaluate.zero_thickness_edge_velocity は node のみ (mesh.discretization: node にすること)")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d0 = design_snapshot(p)                        # 設計点 (作動点で上書きされる前に保存)
    opinfo = select_operating_point(p, op)
    st = gas_states(p)
    kern, design, fr_moc, theta_b = design_from_problem(p, design=d0)
    H = float(p.spec["H_m"])
    m = p.mesh
    mp = SernMeshParams(ni_up=int(m.get("ni_up", 16)), ni_noz=int(m.get("ni_noz", 120)), ni_plume=int(m.get("ni_plume", 220)),
                        nj_top=int(m.get("nj_top", 101)), nj_bot=int(m.get("nj_bot", 61)), L_up=float(m.get("L_up", 0.5)),
                        x_out_extra=float(m.get("x_out_extra", 2.0)), bot_depth=float(m.get("bot_depth", 3.0)),
                        first_wall_frac=float(m.get("first_wall_frac", 2e-3)),
                        cowl_thickness=float(m.get("cowl_thickness", 2e-3 if m.get("discretization", "cell") == "node" else 0.0)),
                        interface_angle=float(m.get("interface_angle_rad", theta_b)),
                        top_ext_angle=float(np.deg2rad(m.get("top_ext_angle_deg", np.rad2deg(design.info["theta_e"])))),
                        ext_top=bool(int(m.get("ext_top", 0))), top_depth=float(m.get("top_depth", 2.0)),
                        nj_ext_top=int(m.get("nj_ext_top", 41)), nj_wake=int(m.get("nj_wake", 9)),
                        vehicle_clearance=float(m.get("vehicle_clearance", 0.06)), first_top_frac=float(m.get("first_top_frac", 0.02)),
                        vehicle_taper=float(m.get("vehicle_taper", 0.0)), t_base=float(m.get("t_base", 0.0)), first_wake_frac=float(m.get("first_wake_frac", 0.0)), split_plume_at_te=bool(m.get("split_plume_at_te", False)),
                        vehicle_wedge_deg=float(m.get("vehicle_wedge_deg", 3.0)), ramp_fillet=float(m.get("ramp_fillet", 0.0)),
                        scale=H)
    if wall_offset:
        design = apply_wall_offset(design, wall_offset, H)
    coords, quads, bedges, minfo, y_mid, y_top = generate_sern_mesh(design, mp)
    write_msh41_named(run_dir / "sern.msh", coords, quads, bedges, PHYS_SERN)
    np.savetxt(run_dir / "ramp_contour.csv", design.ramp_xy * H, delimiter=",", header="x_m,y_m", comments="")
    np.savetxt(run_dir / "cowl_contour.csv", design.cowl_xy * H, delimiter=",", header="x_m,y_m", comments="")
    n = int(nsteps or p.evaluate.get("nStepOuter", 6000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl = float(p.evaluate.get("cfl_main", 4.0))
    cfg = _solver_config(p, n, out_int, cfl, st["ext"]["P"])
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p, st))
    (run_dir / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    write_species_db(p, run_dir, frozen_gases(p))     # R3: species_meta.yaml (lump の NASA-9 はソルバが合成、cpg なら何も書かない)
    disc = p.mesh.get("discretization", "cell")
    # 品質ゲートは primal (cell) 変換で
    (run_dir / "solverConfig.yaml").write_text(qc_cell_config(cfg, disc))
    # 品質ゲートの primal (cell) 変換では farfield を slip に読み替える: farfield は node 専用で、変換器も境界の対応範囲を
    # 検査して止まる (2026-10-05、R7a で初めて farfield 入りの生産 YAML から格子を作って発覚)。どちらも壁ではないので
    # 壁距離・品質判定は同じ。node の本変換の前に元の bcond に戻す
    _bc = (run_dir / "bcondConfig.yaml").read_text()
    (run_dir / "bcondConfig.yaml").write_text(_bc.replace("kind: farfield", "kind: slip"))
    convert_mesh(run_dir, "sern.msh", "sern_qc.h5")
    (run_dir / "bcondConfig.yaml").write_text(_bc)
    # AR 上限は問題 YAML の `mesh.ar_max` で緩められる (既定 1000)。**壁法線に沿った構造格子の
    # 境界層セルに限り 5000 まで** (AGENTS.md「メッシュ品質チェック」2026-09-12 ユーザ決定)。
    # 他の設計チェーン (`runner_axismach` / `runner_wt`) は既にこの knob を持っている。
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "sern_qc.h5", "--mode", "2d",
                        "--ar-max", str(int(p.mesh.get("ar_max", 1000)))], cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    if disc == "cell":
        (run_dir / "sern_qc.h5").rename(run_dir / MESH)
    else:
        (run_dir / "sern_qc.h5").unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg)
        convert_mesh(run_dir, "sern.msh", MESH)
    for f in run_dir.glob("sern_qc.xmf"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg)
    if zte is not None:      # 処置の節点フィールドを変換した格子に書く (格子を作るたびに作り直す)
        mark_zte_field(run_dir, zte)
    # `mesh.ic: moc` で排気側を MOC 場から与える (既定 uniform)。一様 IC は入口状態を全域に置くので
    # 出口で 17 倍ずれており、梯子 12000 step の主因になっている (2026-09-19)
    _ic = str(p.mesh.get("ic", "uniform")).lower()
    _gs = frozen_gases(p)
    n_moc = paste_region_ic(run_dir / MESH, y_mid, y_top, H, st, p.gamma,
                           kern=(kern if _ic == "moc" else None),
                           gas=((_gs or {}).get("exhaust")))
    stamp_region_ic_species(run_dir / MESH, run_dir, st, _gs)     # 新規初期場の化学種属性 (TP のみ)
    ex = st["exhaust"]
    F_ideal_nd, M_e_id = ideal_thrust(p, st)
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n, "H_m": H, "states": st, "gas_model": st["gas_model"],
            "ic": {"mode": _ic, "n_moc_nodes": int(n_moc)},
            "operating_point": opinfo, "wall_offset": bool(wall_offset), "design_point": d0,
            "design": {"key_point": list(design.key_point), "foot_a": list(design.foot_a), "lip_e": list(design.lip_e),
                       "L_ramp": design.L_ramp, "mass_fraction_check": design.mass_fraction_check,
                       "theta_e_deg": float(np.rad2deg(design.info["theta_e"])), "theta_b_deg": float(np.rad2deg(theta_b)),
                       "p_te_over_p_in": kern.p_te_over_p_in, "warnings": design.info["warnings"]},
            "moc_forces": fr_moc, "F_ideal_N_per_m": F_ideal_nd * ex["P"] * H, "M_e_ideal": M_e_id,
            "mesh": minfo, "discretization": disc, "model": p.evaluate.get("model", "euler")}
    if zte is not None:      # IC・化学種属性を書いた後でも w が残っていることを確かめ、再計算した識別量を来歴に残す
        info["zero_thickness_edge_velocity"] = {"requested": zte, "signature": zte_signature_of(zte), "field": ZTE_FIELD,
                                                "identity": verify_zte_field(run_dir, zte)}
    check_zte_stage(cfg, (run_dir / "bcondConfig.yaml").read_text(), run_dir)     # 本段の起動条件 (無効なら何もしない)
    (run_dir / "solverConfig_main.yaml").write_text(cfg)
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def _species_signature(run_dir) -> dict | None:
    """run dir の輸送種の署名を**実 config + 解決済み熱物性** (forge_species.run_thermo) から作る (codex result-2 M2): 種順序・MW・全区間の NASA-9 係数・
    温度区切り・thermoHrefTemp・tracer 設定。`species_meta.yaml` があれば順序の矛盾を拒否。CPG (thermalMethod≠2) は None。
    TP なのに config/DB が読めなければ ValueError (照合不能)。"""
    from ..gas.composition import load_species_meta, load_yaml_str
    rd = Path(run_dir)
    cfgp = (rd / "solverConfig_main.yaml") if (rd / "solverConfig_main.yaml").exists() else (rd / "solverConfig.yaml")
    if not cfgp.exists():
        raise ValueError(f"{rd}: solverConfig.yaml が無く種配置を照合できない")
    cfg = load_yaml_str(cfgp.read_text()); pp = cfg.get("physProp", {})
    if int(pp.get("thermalMethod", 0)) != 2:
        return None
    names = [(str(k["name"]) if isinstance(k, dict) else str(k)).upper() for k in (pp.get("species") or ["N2"])]
    # 熱物性は共通の読み出し forge_species.run_thermo (plan thermophysics-solver-owned-species-db #8): ソルバの解決済み記録 >
    # 従来の speciesDBFile > forge --resolve-species。lump 記法 ({name, lump, basis}) の config も読める
    try:
        th = _forge_species().run_thermo(rd)
    except ValueError as e:
        raise ValueError(f"{rd}: 熱物性を解決できず種配置を照合できない ({e})") from None
    if [str(k).upper() for k in th["names"]] != names:
        raise ValueError(f"{rd}: 解決済み熱物性の種 {th['names']} が {cfgp.name} の {names} と違う")
    ents = {}
    for k, n in zip(names, th["names"]):
        e = th["species"][n]
        # 2 区間 (Tlo/Tmid/Thi, nasa9_low/high) と区間可変 (Tbounds, nasa9_intervals; 種 DB 段 3 の解決済み記録) の両方を読む
        Tb, co = _forge_species().nasa9_intervals(e)
        ents[k] = {"MW": float(e["MW"]), "coefs": co, "ranges": Tb}
    meta = load_species_meta(rd)
    if meta is not None and [str(k).upper() for k in meta["species"]] != names:
        raise ValueError(f"{rd}: species_meta.yaml の種順序 {meta['species']} が solverConfig の {names} と矛盾")
    tracer = str(pp.get("tracer", "none")).lower() not in ("none", "", "0")
    if meta is not None and bool(meta.get("tracer", {}).get("enabled", False)) != tracer:
        raise ValueError(f"{rd}: species_meta.yaml のトレーサ設定が solverConfig (tracer: {pp.get('tracer', 'none')}) と矛盾")
    return {"species": names, "entries": ents, "href": float(pp.get("thermoHrefTemp", 0.0)), "tracer": tracer}


def check_species_compatible(src_run_dir, dst_run_dir, what: str = "restart", allow_db_change: bool = False) -> None:
    """restart 経路の共通照合 (codex result M3 / result-2 M2): 実 config + DB の署名で種順序・MW・NASA-9 係数・温度区切り・
    datum・トレーサを照合し、一致しないと拒否。allow_db_change (作動点変更の warm start) は lump の MW/係数の変化だけ許す。
    両方 CPG なら通す。片方だけ TP は拒否。"""
    a, b = _species_signature(src_run_dir), _species_signature(dst_run_dir)
    if a is None and b is None:
        return
    if a is None or b is None:
        raise ValueError(f"{what}: 片方だけ TP (thermalMethod 2) で種配置を照合できない ({src_run_dir} → {dst_run_dir})")
    if a["species"] != b["species"]:
        raise ValueError(f"{what}: 輸送種の順序が違う (元 {a['species']} / 先 {b['species']}); tools/convert_species_field.py を使う")
    if not allow_db_change:
        for k in a["species"]:
            ea, eb = a["entries"][k], b["entries"][k]
            if abs(ea["MW"] / eb["MW"] - 1.0) > 1e-9:
                raise ValueError(f"{what}: 種 {k} の MW が違う ({ea['MW']} / {eb['MW']}) — DB が異なる")
            if ea["ranges"] != eb["ranges"]:
                raise ValueError(f"{what}: 種 {k} の温度区切りが違う ({ea['ranges']} / {eb['ranges']})")
            for i, (ca, cb) in enumerate(zip(ea["coefs"], eb["coefs"])):
                if any(abs(x - y) > 1e-12 * max(abs(x), abs(y), 1.0) for x, y in zip(ca, cb)):
                    raise ValueError(f"{what}: 種 {k} の NASA-9 係数 (区間 {i}) が違う — DB が異なる")
    if abs(a["href"] - b["href"]) > 1e-9:
        raise ValueError(f"{what}: thermoHrefTemp が違う ({a['href']} / {b['href']})")
    if a["tracer"] != b["tracer"]:
        raise ValueError(f"{what}: トレーサの有無が違う (元 {a['tracer']} / 先 {b['tracer']})")


def _require_datasets(h5, names, what: str) -> None:
    missing = [k for k in names if k not in h5["VALUE"]]
    if missing:
        raise ValueError(f"{what}: 元の VALUE に {missing} が無い")


def restart_by_index(res_h5, mesh_h5) -> None:
    """同一メッシュの stage 間移植: VALUE を index でコピーする (座標最近傍の `interp_field.py` は使わない)。
    理由 (2026-09-04, case/46 run_0009): スリットカウルの上下壁ノードは座標が一致し、最近傍補間が双子を同じ元
    ノードに写す → 排気側の壁ノードが外部流の圧力を持ち 2 次で発散した (interp_field の全 134 station で誤写像を確認)。"""
    check_species_compatible(Path(res_h5).parent, Path(mesh_h5).parent, "restart_by_index")
    sig = _species_signature(Path(mesh_h5).parent)
    # 化学種の属性 (§4.3): SRC の記録を検証し、宛先を --resolve-species で解決して一致なら継承 (不一致は書き込み前に停止)
    fsp = _forge_species()
    try:
        species_plan = fsp.plan_inherit(res_h5, Path(mesh_h5).parent, tool="restart_by_index")
    except fsp.SpeciesCheckError as e:
        raise ValueError(f"restart_by_index: {e}") from None
    fsp.write_species_attrs(mesh_h5, None)
    with h5py.File(res_h5, "r") as src, h5py.File(mesh_h5, "r+") as dst:
        n = len(dst["VALUE/ro"])
        if sig is not None:
            _require_datasets(src, ["ro", "roUx", "roUy", "roUz", "roe"] + [f"roY{i}" for i in range(len(sig["species"]))] + (["roXi"] if sig["tracer"] else []), "restart_by_index")
        keys = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]        # 状態量のみ (wall_dist は触らない)
        keys += sorted(k for k in src["VALUE"] if re.fullmatch(r"roY\d+", k)) + ["roXi"]   # 化学種・トレーサも引き継ぐ (旧: 7 変数のみで ΣY が壊れた; codex M4)
        for k in keys:
            if k in src["VALUE"] and len(src["VALUE"][k]) == n:
                if k not in dst["VALUE"]:
                    if k.startswith("roY") or k == "roXi":
                        dst["VALUE"].create_dataset(k, data=np.asarray(src["VALUE"][k][:], dtype=np.float32))
                    continue
                dst["VALUE"][k][:] = src["VALUE"][k][:]
        fsp.commit_inherit(dst, species_plan)


def warm_from_run(dst_run_dir, src_run_dir) -> dict:
    """**別作動点**の収束場を同一メッシュの run へ移植する (plan §4.7 / §5.1-1c、codex レビューの A′)。

    保存量の素コピー (`restart_by_index`) は**作動点間では使えない**: 同じ `roe` を別の γ で読むと
    P = (γ−1)(roe − ½ρ|u|²) が (γ_dst−1)/(γ_src−1) 倍ずれる (m6_on→m10_on で 1.24 倍、→m4_off で 2.18 倍)。
    そこで入口状態の比で相似スケールし、**目標作動点の γ で roe を組み直す**:

        ρ' = ρ·sρ,  (ρu)' = (ρu)·sρ·su,  P' = P·sP,  roe' = P'/(γ_dst−1) + ½|ρu'|²/ρ'
        (ρk)' = (ρk)·sρ·su²   (k ~ u²),   (ρω)' = (ρω)·sρ·su   (ω ~ u/L, L は同一メッシュで不変)

    sρ, su, sP は入口 (燃焼器出口) 状態の比。NPR が違えば内部の波の当たり方は変わるので、
    移植後に**短い適応段**を必ず入れる (run_staged の warm_adapt_steps)。戻り値 = 使ったスケール。
    """
    dst_run_dir, src_run_dir = Path(dst_run_dir), Path(src_run_dir)
    si = json.loads((src_run_dir / "prepare_info.json").read_text())
    di = json.loads((dst_run_dir / "prepare_info.json").read_text())
    if si["mesh"]["cells"] != di["mesh"]["cells"] or si["design"]["L_ramp"] != di["design"]["L_ramp"]:
        raise ValueError("warm_from_run: メッシュが同一でない (同一設計・同一 dv の run 間でのみ使う)")
    se, de = si["states"]["exhaust"], di["states"]["exhaust"]
    g_s = float(si["operating_point"]["gas"]["gamma"]); g_d = float(di["operating_point"]["gas"]["gamma"])
    s_ro, s_u, s_P = de["ro"] / se["ro"], de["u"] / se["u"], de["P"] / se["P"]
    res = sorted(src_run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not res:
        raise RuntimeError(f"warm_from_run: {src_run_dir} に res_*.h5 が無い")
    gases_d = None
    if di.get("gas_model") == "frozen_tp":
        # 順序・datum・トレーサは一致必須。MW (lump の組成 = 作動点の排気組成) は違ってよい: roe は目標作動点の DB で T から組み直す
        check_species_compatible(src_run_dir, dst_run_dir, "warm_from_run", allow_db_change=True)   # codex result M3
        # 作動点適用後の組成で擬似種を作る (prepare_info の problem は作動点未適用の YAML なので op を再選択)
        pd_ = load_problem(di["problem"]); select_operating_point(pd_, di["operating_point"]["name"]); gases_d = frozen_gases(pd_)
    # 化学種の属性 (§4.3「種変換」と同じ扱い): roe は目標作動点の物性で作り直すので、入口 (元 res) が検証済みなら
    # 書き込み後に宛先の記録と照合して付ける。記録が壊れていれば止める。元が未検証で宛先が TP なら既定で止める
    # (ソルバと同じ規約, #3c); その実行だけ FORGE_ALLOW_UNVERIFIED_SPECIES=1 で許可したときは宛先も未検証 (属性なし)
    fsp = _forge_species()
    src_state = fsp.source_species_state(res[-1])
    if src_state["state"] == "broken":
        raise ValueError(f"warm_from_run: 元 res の化学種記録を検証できない (照合不能): {src_state['why']}")
    if gases_d is not None and src_state["state"] in ("none", "unverified"):
        try:
            fsp.refuse_unverified("warm_from_run", f"SRC {res[-1]} is unverified ({src_state['why']}) and destination "
                                                   f"{dst_run_dir} is thermally perfect (frozen_tp)")
        except fsp.SpeciesCheckError as e:
            raise ValueError(f"warm_from_run: {e}") from None
    fsp.write_species_attrs(dst_run_dir / MESH, None)
    with h5py.File(res[-1], "r") as src, h5py.File(dst_run_dir / MESH, "r+") as dst:
        n = len(dst["VALUE/ro"])
        if len(src["VALUE/ro"]) != n:
            raise ValueError("warm_from_run: VALUE 長が一致しない")
        ro = src["VALUE/ro"][:].astype(np.float64)
        mom = [src[f"VALUE/{k}"][:].astype(np.float64) for k in ("roUx", "roUy", "roUz")]
        roe = src["VALUE/roe"][:].astype(np.float64)
        ro_n = ro * s_ro
        mom_n = [m * (s_ro * s_u) for m in mom]
        if gases_d is None:
            P = (g_s - 1.0) * (roe - 0.5 * sum(m * m for m in mom) / np.maximum(ro, 1e-30))
            roe_n = P * s_P / (g_d - 1.0) + 0.5 * sum(m * m for m in mom_n) / np.maximum(ro_n, 1e-30)
        else:
            # frozen_tp (R3): 圧力は出力の P を相似スケール、組成 (Y_EXH, Y_AMB) は場のまま持ち越し、
            # T' = P'/(ρ' R_mix(Y)) と目標作動点の擬似種 (排気組成が違う) で roe' = ρ'(Σ Y_s e_sens,s(T') + ½|u'|²) を組み直す
            # 組成の再初期化 (codex result-2 M1): 元の組成は排気率 ξ 以外捨て、**目標作動点**の入口ベクトルから
            # Y_t = ξ Y_in^dst + (1−ξ) Y_ext^dst を組む (lumped [EXH, AMB] では Y_EXH の持ち越しと同値、full / lumped+keep では
            # 実種分率が新作動点の排気組成に変わる)。ξ は元 run の exhaust_fraction (tracer なら roXi/ρ、無ければ流入元ラベル種)
            from ..gas.composition import exhaust_fraction, reinit_transport_vector
            tg = gases_d["transported"]; Ld = gases_d["layout"]; names = list(Ld.species)
            spec = exhaust_fraction(src_run_dir)
            _require_datasets(src, ["P", spec["conserved"]], "warm_from_run (frozen_tp)")
            xi = np.clip(src["VALUE/" + spec["conserved"]][:].astype(np.float64) / np.maximum(ro, 1e-30), 0.0, 1.0)
            Ys = reinit_transport_vector(xi, Ld)
            P = src["VALUE/P"][:].astype(np.float64) * s_P
            R_mix = sum(y * g.R for y, g in zip(Ys, tg))
            T_n = P / np.maximum(ro_n * R_mix, 1e-30)
            e_n = sum(y * g.e_sens(T_n) for y, g in zip(Ys, tg))
            roe_n = ro_n * e_n + 0.5 * sum(m * m for m in mom_n) / np.maximum(ro_n, 1e-30)
            if not (np.all(np.isfinite(roe_n)) and np.all(np.isfinite(T_n)) and np.all(ro_n > 0)):
                raise ValueError("warm_from_run: 非有限または非正の状態 (書込み前に中止)")
            for i_, y in enumerate(Ys):
                if f"roY{i_}" in dst["VALUE"]:
                    dst[f"VALUE/roY{i_}"][:] = ro_n * y
                else:
                    dst["VALUE"].create_dataset(f"roY{i_}", data=(ro_n * y).astype(np.float32))
            if Ld.tracer:
                if "roXi" in dst["VALUE"]:
                    dst["VALUE/roXi"][:] = ro_n * xi
                else:
                    dst["VALUE"].create_dataset("roXi", data=(ro_n * xi).astype(np.float32))
        dst["VALUE/ro"][:] = ro_n
        for k, m in zip(("roUx", "roUy", "roUz"), mom_n):
            dst["VALUE"][k][:] = m
        dst["VALUE/roe"][:] = roe_n
        if "roK" in src["VALUE"] and "roK" in dst["VALUE"]:
            dst["VALUE/roK"][:] = src["VALUE/roK"][:] * (s_ro * s_u * s_u)
            dst["VALUE/roOmega"][:] = src["VALUE/roOmega"][:] * (s_ro * s_u)
    species = "unverified"
    if gases_d is not None and src_state["state"] == "verified":
        tg = gases_d["transported"]; Ld = gases_d["layout"]; names = list(Ld.species)
        mixes = [(f"species {k}", [1.0 if j == i else 0.0 for j in range(len(names))], g.R, g.e_sens)
                 for i, (k, g) in enumerate(zip(names, tg))]
        species = fsp.stamp_new_field(dst_run_dir / MESH, dst_run_dir, names, [float(Ld.entries[k].MW) for k in names],
                                      gases_d["href_T"], mixes, tool="warm_from_run")
    return {"src": str(src_run_dir), "s_ro": s_ro, "s_u": s_u, "s_P": s_P, "gamma": [g_s, g_d], "species_attrs": species}


def run_forge(run_dir) -> int:
    _record_zte_launch(run_dir)     # 処置を有効にした段だけ: 起動条件の検査と、読む w の識別量の記録 (無効の段は何もしない)
    r = subprocess.run([str(FORGE_TOOLS / "run_case.sh"), str(Path(run_dir).resolve())], env=_ENV, capture_output=True, text=True)
    (Path(run_dir) / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode



def _archive_stage(run_dir, tag: str) -> None:
    """段の `residual_history.csv` を `residual_<tag>.csv` に退避する (2026-09-19, codex plan レビュー M5)。
    forge は段ごとに上書きするので、退避しないと**どの段で残差が上がったかを後から追えない**。
    設計 B の `rms_roY1` 上昇を「soft/mid で起きた」と断じた根拠が無かったのはこれが理由。"""
    src = Path(run_dir) / "residual_history.csv"
    if src.exists():
        try:
            (Path(run_dir) / f"residual_{tag}.csv").write_bytes(src.read_bytes())
        except OSError:
            pass


def run_staged(run_dir, stages: str = "full", soft_steps: int = 3000, soft_cfl: float = 0.5, soft_conv: int = 0,
               warm_lam_steps: int = 0, warm_lam_cfl: float = 0.2, mid_steps: int = 0,
               warm_lam_ramp=None, soft_ramp=None,
               warm_src=None, warm_adapt_steps: int = 500) -> int:
    """soft_cfl / soft_conv: soft 段の CFL と convMethod (既定 0.5 / 1 次)。3D SST の後縁 3 重点など、1 次でも
    立ち上がりが厳しいケースで下げる。

    warm_lam_steps > 0 (SST のみ): soft 段の**前に層流暖機段** (`turbulence: none` + 粘性あり、1 次、`warm_lam_cfl`) を
    入れる (plan §4.7)。カウル後縁のせん断層は排気と外部流の密度・温度比が大きく、平均場が立つ前に SST を回すと
    `roOmega` が発散する — 板厚では作動点ごとに要求が逆転して解けなかった (m4_off は 2e-3、m10_on は 5e-3 が必要)。
    暖機で平均場を作ってから乱流方程式を入れると 3 作動点とも通り、力係数は暖機なしの成功例と一致する (case/46 run_0048)。

    mid_steps > 0: soft と本段の間に **2 次 + soft CFL** の段を入れる。次数と CFL を同時に上げるとランプ膨張角部で
    `roOmega` が発散する (case/46 run_0054 doe_000, θ_r0 18.3°)。`procedures/solver-settings.md` の段階戦略に沿う。

    warm_src: **別作動点の収束 run ディレクトリ**。与えると層流暖機と soft を飛ばし、`warm_from_run` で
    熱力学整合リマップ → **適応段** (1 次, soft CFL, warm_adapt_steps) → mid → 本段 とする。
    NPR が大きく違う作動点間 (m6_on 35.4 ↔ m4_off 2.8) では使わないこと (§5.1-1c)。"""
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig_main.yaml").read_text() if (run_dir / "solverConfig_main.yaml").exists() \
        else (run_dir / "solverConfig.yaml").read_text()
    if stages == "main":  # soft 段の res_*.h5 が残っている状態から本段だけ回す
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg_main)
        return run_forge(run_dir)
    if stages == "none":
        return run_forge(run_dir)
    if warm_src is not None:      # 別作動点からの warm start: 暖機と soft を飛ばし、適応段だけ入れる
        info_warm = warm_from_run(run_dir, warm_src)
        (run_dir / "WARM_START.json").write_text(json.dumps(info_warm, indent=1))
        ad = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {soft_cfl}, cfl_pseudo: {soft_cfl}", cfg_main)
        ad = ad.replace("convMethod: 1", f"convMethod: {soft_conv}")
        ad = re.sub(r"nStepOuter: \d+", f"nStepOuter: {warm_adapt_steps}", ad)
        ad = re.sub(r"outStepInterval: \d+", f"outStepInterval: {warm_adapt_steps}", ad)
        (run_dir / "solverConfig.yaml").write_text(ad)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError("warm start の適応段が失敗 (res_nan_*.h5 / forge_run.log を見る)")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
        warm_lam_steps = 0      # 以降は mid → 本段
    if warm_lam_steps > 0 and 'model: "sst"' in cfg_main:      # 層流暖機段 (SST を後から入れる)
        lam0 = re.sub(r'turbulence: \{model: "sst"[^}]*\}', 'turbulence: {model: "none"}', cfg_main)
        if 'model: "none"' not in lam0:
            raise RuntimeError("層流暖機: turbulence 行の置換に失敗 (solverConfig の書式が変わった)")
        lam0 = lam0.replace("convMethod: 1", "convMethod: 0")
        # **CFL ramp** (2026-09-19 ユーザ提案): 暖機は固定 CFL だと 0.1 が上限だが、場が育つにつれ上げられる。
        # ソルバ側に ramp が無いので runner が forge を複数回起動して段階昇圧する
        # (風洞チェーン `runner_axismach.run_staged_ns(stages="ramp")` と同じ方式)。
        # `warm_lam_ramp` が空なら従来どおり `warm_lam_cfl` 固定の 1 段。
        legs = [(float(c), int(warm_lam_steps / max(len(warm_lam_ramp), 1))) for c in warm_lam_ramp] \
            if warm_lam_ramp else [(float(warm_lam_cfl), int(warm_lam_steps))]
        for c, n_leg in legs:
            lam = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {c}, cfl_pseudo: {c}", lam0)
            lam = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n_leg}", lam)
            lam = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n_leg}", lam)
            (run_dir / "solverConfig.yaml").write_text(lam)
            rc = run_forge(run_dir)
            res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
            if rc != 0 or not res:
                raise RuntimeError(f"層流暖機段が失敗 (cfl {c}; res_nan_*.h5 / forge_run.log を見る)")
            _archive_stage(run_dir, f"warm_cfl{c:g}")
            restart_by_index(res[-1], run_dir / MESH)
            for f in run_dir.glob("res_*"):
                f.unlink()
    if warm_src is not None:      # soft は適応段で代替済み → mid へ
        return _run_mid_and_main(run_dir, cfg_main, soft_cfl, mid_steps)
    # soft 段も **CFL ramp** できる (`opt.soft_ramp`, 2026-09-19)。暖機で ramp が効いた (固定 0.5 は step 30 で
    # 発散するのに ramp なら 1.0 まで到達) のと同じ理屈。固定 soft_cfl 1.0 は設計によって
    # `rms_roY1` (排気∩外気のせん断層) が上昇するので、段階昇圧で通す。
    _soft_legs = [(float(c), max(int(soft_steps / len(soft_ramp)), 1)) for c in soft_ramp] \
        if soft_ramp else [(float(soft_cfl), int(soft_steps))]
    for c, n_leg in _soft_legs:
        soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {c}, cfl_pseudo: {c}", cfg_main)
        soft = soft.replace("convMethod: 1", f"convMethod: {soft_conv}")
        soft = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n_leg}", soft)
        soft = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n_leg}", soft)
        (run_dir / "solverConfig.yaml").write_text(soft)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError(f"soft 段が失敗 (cfl {c}; res_nan_*.h5 / forge_run.log を見る)")
        _archive_stage(run_dir, f"soft_cfl{c:g}")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
    return _run_mid_and_main(run_dir, cfg_main, soft_cfl, mid_steps, soft_ramp)


def _run_mid_and_main(run_dir, cfg_main: str, soft_cfl: float, mid_steps: int, soft_ramp=None) -> int:
    """mid 段 (2 次 + soft CFL) → 本段。soft/暖機/warm start の後段として共用する。"""
    if mid_steps > 0:      # mid 段: 2 次に上げるが CFL は soft のまま (次数と CFL を同時に上げない)
        legs = [(float(c), max(int(mid_steps / len(soft_ramp)), 1)) for c in soft_ramp] \
            if soft_ramp else [(float(soft_cfl), int(mid_steps))]
        for _ci, (c, n_leg) in enumerate(legs):
            mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {c}, cfl_pseudo: {c}", cfg_main)
            mid = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n_leg}", mid)
            mid = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n_leg}", mid)
            (run_dir / "solverConfig.yaml").write_text(mid)
            rc = run_forge(run_dir)
            _archive_stage(run_dir, f"mid_cfl{c:g}")
            if _ci < len(legs) - 1:
                _r = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(ch for ch in f.stem if ch.isdigit())))
                if rc != 0 or not _r:
                    raise RuntimeError(f"mid 段が失敗 (cfl {c})")
                restart_by_index(_r[-1], run_dir / MESH)
                for f in run_dir.glob("res_*"):
                    f.unlink()
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res:
            raise RuntimeError("mid 段 (2 次 + soft CFL) が失敗 (res_nan_*.h5 / forge_run.log を見る)")
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)


def collect(problem_path, run_dir, out_dir=None, rc=None, require_residual_pass: bool = False) -> dict:
    """力係数履歴 + 受理ゲート (plan §5.1 R1)。metrics.json / force_history.csv は out_dir (既定 run_dir) に書く
    (再判定で元 run を汚さない)。rc は forge の終了コード (None なら run_case_stdout.log から読む)。"""
    p = load_problem(problem_path)
    run_dir = Path(run_dir); out_dir = Path(out_dir) if out_dir else run_dir; out_dir.mkdir(parents=True, exist_ok=True)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    st = info["states"]; ex, en = st["exhaust"], st["ext"]; H = info["H_m"]
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    hist = force_history(run_dir, p_a=en["P"], F_ideal=info["F_ideal_N_per_m"], H=H, x_ref=float(xr) * H, y_ref=float(yr) * H,
                         mdot_u_in=ex["ro"] * ex["u"] ** 2 * H, p_in=ex["P"],
                         twall_on_fluid=(info.get("discretization", "cell") == "cell"))
    if rc is None:
        rc = forge_rc_from_log(run_dir)
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-2:] if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else []
    gates = evaluate_gates(run_dir, hist, rc, require_residual_pass=require_residual_pass,
                           p_min=float(p.evaluate.get("p_min", 1.0)))
    out = {"convergence_verdict": verdict, "n_snapshots": len(hist), "history": hist, "forge_rc": rc,
           "operating_point": info.get("operating_point"), "L_ramp": info["design"]["L_ramp"],
           "gates": gates, "steadiness": gates["steadiness"]["series"], "objective": gates["objective"]}
    if hist:
        last = hist[-1]
        out.update({k: last[k] for k in ("step", "C_T", "C_T_wall", "C_L", "C_M", "T_wall", "L", "M_noseup")})
        for k in ("C_T_with_shear", "C_T_friction", "sep_frac_ramp", "sep_x_min_ramp"):
            if k in last:
                out[k] = last[k]
        out["moc_forces"] = info["moc_forces"]
        # MOC は**設計点**の値なので、作動点が設計点と一致するときだけ差を出す (2026-09-05, plan §4.10:
        # 作動点は inflow/gas も動かすので、オフデザイン run で差を取ると意味の無い数になる)。
        d0, o = info.get("design_point"), info.get("operating_point") or {}
        on_design = bool(d0) and all(o.get(k) == d0.get(k) for k in ("inflow", "external")) \
            and (o.get("gas", {}).get("gamma") == d0.get("gamma"))
        out["on_design_point"] = on_design
        if on_design:
            out["cfd_vs_moc"] = {k: (last[k] - info["moc_forces"][k]) for k in ("C_T", "C_L", "C_M")}
        write_force_history_csv(out_dir / "force_history.csv", hist)
    # 実効 slauWallNormalChi と設定方針 (plan convection-slau-wall-normal-chi-default §4.4、codex plan M3)。
    # 起動記録 forge_launches.jsonl の**最後の起動** (本段) の値。記録が無い run (旧バイナリ) は None = 不明。
    out["slau_wall_normal_chi_effective"] = _last_launch_chi(run_dir)
    # 実効 mesh.scalarGradient (plan gradient-scalar-lsq-unification #6、codex diagnose 2026-09-27)。記録が無ければ None = 不明。
    out["scalar_gradient_effective"] = _last_launch_value(run_dir, "scalarGradient", allowed=("gg", "lsq"))
    # 厚さ 0 の板の自由端の処置 (plan convection-zero-thickness-edge-reconstruction §4): 有効状態・タグ・rings・格子署名・
    # w のハッシュ・起動の記録と、設計 DB の識別子 (処置の有無・タグ・rings が違う評価を同じ設計点として混ぜない)
    _zp = zte_provenance(run_dir)
    out["zero_thickness_edge_velocity"] = _zp
    out[ZTE_EFFECTIVE] = _zp["effective"]
    out["flag_policy"] = FLAG_POLICY
    (out_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="forge_design ⑤ SERN 評価 (1 作動点)")
    ap.add_argument("problem"); ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--stages", default="full", choices=["full", "none", "main"])
    ap.add_argument("--op", default=None, help="作動点名 (spec.operating_points)")
    ap.add_argument("--wall-offset", default=None, help="δ* オフセット JSON ({ramp: [[x_m, dn_m],...], cowl: [...]})")
    a = ap.parse_args(argv)
    wo = json.loads(Path(a.wall_offset).read_text()) if a.wall_offset else None
    info = prepare(a.problem, a.run_dir, a.steps, op=a.op, wall_offset=wo)
    print(json.dumps({k: info[k] for k in ("design", "moc_forces", "mesh")}, indent=1))
    if a.prepare_only:
        return 0
    # 段階起動のパラメータは problem の `opt:` に置いてある (driver_sern と同じ正本)。
    # ここで読まないと既定 (層流暖機なし・soft CFL 0.5) で回り、生産と別物になる (run_0078 で 1 本無駄にした)。
    o = (load_problem(a.problem).raw.get("opt") or {})
    rc = run_staged(a.run_dir, a.stages,
                    soft_steps=int(o.get("soft_steps", 3000)), soft_cfl=float(o.get("soft_cfl", 0.5)),
                    warm_lam_steps=int(o.get("warm_lam_steps", 0)), warm_lam_cfl=float(o.get("warm_lam_cfl", 0.2)),
                    warm_lam_ramp=o.get("warm_lam_ramp"), soft_ramp=o.get("soft_ramp"),
                    mid_steps=int(o.get("mid_steps", 0)),
                    warm_adapt_steps=int(o.get("warm_adapt_steps", 500)))
    out = collect(a.problem, a.run_dir, rc=rc, require_residual_pass=bool(o.get("require_residual_pass", False)))
    print(json.dumps({k: v for k, v in out.items() if k not in ("history", "gates")}, indent=1))
    g = out["gates"]; print(f"GATES: {g['verdict']} fail_class={g['fail_class']} objective={g['objective']} reasons={g['reasons']}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
