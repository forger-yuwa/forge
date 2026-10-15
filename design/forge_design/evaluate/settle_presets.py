r"""収束の判定の基準 `settle_criteria.yaml` のプリセット (plan `plans/active/tooling-settle-convergence-and-monitor.md` §4.2・§5.1 #2〜#3・#5d)。

問題 YAML の `evaluate.settle` (`preset` + 足し引き・閾値の上書き) を、準備のときに run の設定 (prepare_info・solverConfig・
bcondConfig) で解決して `settle_criteria.yaml` の中身にする。判定は `solver_density_cuda/tools/settle.py` (`check_settled.py`)。

`evaluate.settle` のキー (知らないキーは ValueError):
  preset: wt_nozzle_v1                 必須 (wt_nozzle_v1 / wt_nozzle_v2)
  thresholds: {<量の名前>: {<閾値のキー>: 値}}   閾値の上書き (知らない量・キーは例外)
  quantities_remove: [<量の名前>]       外す (criteria の excluded に「問題 YAML で外した」と残す)
  quantities_add: [<量の指定>]          足す (settle_criteria.yaml の quantities の 1 要素の形; 検査は settle.validate_criteria)
  x_range_m: [x0, x1]                  試験部の線の x 範囲 [m] (既定 (x_E + 2, x_F − 1) r_t)
  line_eta: 0.03                       試験部の線の r/r_w (既定 0.03; ユーザ決定 2026-10-16)

## プリセット `wt_nozzle_v1` (要求カードの風洞ノズル; 閾値はすべて暫定 `provisional: true`、plan §6 V3 で決める)

量 6 つ: residual・mass_balance・exit_core_mach・line_mach (r/r_w = 0.03)・exit_flow_angle・wall_heat_nodes。
- 構造格子 ni・nj・scale_m は prepare_info、`structured.mesh_file` は solverConfig の `mesh.meshFileName` (節点の座標を
  `nozzle_report.load_field` と同じ格子の h5 の COORD から作る; float のビルドでも一致させるため)。
- 判定の格子 = solverConfig の `time.outStepInterval`、窓 W 20,000 (min_points 11)、保持 40,000 (停止規則用、記録だけ)。
  `grid.start_step` は 0 (forge の起動ごとに step は 0 から)、`grid.segment` は判定区間の識別 (設計チェーンは `main_stage_segment()` =
  本段の開始の通知 MAIN_STAGE_STARTED を指す; §5.1 #6)。
- 残差の必要な列: rms_ro・rms_roUx・rms_roUy・rms_roe、SST なら rms_roK・rms_roOmega、遷移 lm2009 なら rms_roGamma・rms_roReth、
  化学種 2 つ以上なら rms_roY0〜、トレーサ exhaust なら rms_roXi (forge `main.cpp` の `residualEquationNames` と同じ規則)。
  非活動を許す列: rms_roUz。凝縮 (condensation 1) はモーメントの列名を決められないので例外 (問題 YAML で residual を外して足す)。
  閾値: slope_dec log10(1.05)・amp_dec 0.1・rise_dec log10(1.10)・excursion_dec 0.5。
- 流量の収支: inlet = bcond の kind が `inlet_*`、outlet = `outlet_*`・`outflow`。`farfield` や知らない kind・入口か出口が無い → ValueError
  (判別できない)。壁・軸・slip は記録だけ、periodic は forge が CSV に書かない。閾値: mean_abs 5e−4・swing 5e−3・mdot_change_rel 1e−4。
- 出口コア M: η ∈ [0.05, 0.7]、tol = 0.1 × 0.05 % × M_d (絶対値; M_d は prepare_info の `Md`)。
- 試験部の線: η 0.03、x = (x_E + 2, x_F − 1) r_t を 2001 点 (`nozzle_report.metrics` の w_test と同じ範囲; x_F は
  prepare_info の `physical_wall.domain[1]` = 格子の最後の列の軸の x)、Md = prepare_info の `Md`。tol_maxabs 0.01・tol_halfamp 0.005・
  tol_pos 0.01 [%pt]。
- 出口の流れ角: η ∈ [0.05, 0.7]、tol 0.01°。
- 壁の熱流: bcond の `wall_isothermal` か `ints.wallProfile` の壁。無ければ外して excluded に理由を書く (黙って外さない)。
  q_tol 1.0e3 W/m²・eps_rel 1e−3。軸対称は solverConfig の mesh.isAxisymmetric。
- 量の前提 (forge の出力の設定; output.boundaryMassflux 1・output.nodeProbeFile・output.interfaceDiag 1・壁の outputHDFflg 1・定常) は
  `preconditions` に満たされているかを記録する (満たされていなければ判定は INDETERMINATE (データ不足) になる; 有効にするのは雛形 #6)。

## プリセット `wt_nozzle_v2` (区間 1 万 step・保持なし; plan §4.4「判定の区間を 1 万 step にする (2026-10-16 ユーザ決定)」、§5.1 #5d)

v1 との違いだけ (量・閾値の値・線・出口コア・入口出口・残差の列・前提の記録は v1 と同じ; 閾値はすべて暫定):
- 窓 W 10,000 (min_points 6 = 間隔 2,000 の端点込みの点)、保持 0 (停止規則 v2 は最初に SETTLED になった格子点で要求)。
  最初に判定できる格子点は開始 + 2W + 間隔 = 22,000。許容は区間の長さで縮めない (ユーザ「区間の中の幅 ≤ 許容の 1/8 これやめろ」)。
- 残差は `series: raw` (毎 step の outer_end の全行で傾き・振幅・上昇; 行の欠け・重複は判定不能)。
- 壁の熱流は `measure: range` (節点ごとの区間の幅 max − min; 相対の参照は区間の終端の値; codex diagnose 2026-10-16 settle-v2-span10k)。
- v1 は run_0013 の検証のために変えない (v1 の解決の結果は今と同じ)。

節点の表 (forge の `output.nodeProbeFile`): `probe_node_ids` が節点の量に要る節点の和集合を返し、`write_node_probe_file` が forge の
書式 (1 行 1 節点番号、`#` 以降と空行は無視、重複・範囲外は forge が起動時に拒否; `output/monitorOutputs.cpp` の readNodeList) で書く。
線は列ごとに 軸 j = 0・η を挟む 2 節点とその両側 (j−1〜j+2; η が節点の η と丸めの差で一致するとき列によって挟む対が
(j0, j1) か (j1, j2) になっても足りるように)・壁 j = nj − 1、線が使う列とその両隣の列。出口の列 i = ni − 1 は全節点。
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml


def _settle():
    """solver_density_cuda/tools/settle.py (判定の共通の部品)。"""
    tools = str(Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3])) / "solver_density_cuda" / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import settle
    return settle


SETTLE_KEYS = ("preset", "thresholds", "quantities_remove", "quantities_add", "x_range_m", "line_eta")

PRESETS = {
    "wt_nozzle_v1": {
        "version": "wt_nozzle_v1 (2026-10-16; 線 η 0.03)",
        "window": {"W": 20000, "min_points": 11},
        "hold": 40000,
        "thresholds": {
            "residual": {"slope_dec": math.log10(1.05), "amp_dec": 0.1, "rise_dec": math.log10(1.10), "excursion_dec": 0.5,
                         "provisional": True},
            "mass_balance": {"mean_abs": 5.0e-4, "swing": 5.0e-3, "mdot_change_rel": 1.0e-4, "provisional": True},
            "exit_core_mach": {"tol_over_Md": 0.1 * 0.0005, "provisional": True},   # tol = 0.1 × 0.05 % × M_d
            "line_mach": {"tol_maxabs": 0.01, "tol_halfamp": 0.005, "tol_pos": 0.01, "provisional": True},
            "exit_flow_angle": {"tol": 0.01, "provisional": True},
            "wall_heat_nodes": {"q_tol": 1.0e3, "eps_rel": 1.0e-3, "provisional": True},
        },
        "exit_eta_range": [0.05, 0.7],
        "line": {"eta": 0.03, "n": 2001, "x_offsets_rt": [2.0, -1.0]},
    },
}
# 区間 1 万 step・保持なし (plan §4.4「判定の区間を 1 万 step にする (2026-10-16 ユーザ決定)」、§5.1 #5d)。閾値は v1 と同じ値
# (区間の長さで許容を縮めない; ユーザ「区間の中の幅 ≤ 許容の 1/8 これやめろ 10000ステップはそのまま」)。残差は毎 step の系列
# (settle の residual.series raw)、壁の熱流は節点ごとの区間の幅 (wall_heat_nodes.measure range; codex diagnose 2026-10-16
# settle-v2-span10k)。v1 は run_0013 の検証のために変えない。
PRESETS["wt_nozzle_v2"] = {
    "version": "wt_nozzle_v2 (2026-10-16; 区間 1 万・保持なし (2026-10-16 ユーザ決定); 残差は毎 step の系列・壁の熱流は区間の幅; 線 η 0.03)",
    "window": {"W": 10000, "min_points": 6},
    "hold": 0,
    "residual_series": "raw",
    "wall_measure": "range",
    "thresholds": copy.deepcopy(PRESETS["wt_nozzle_v1"]["thresholds"]),
    "exit_eta_range": [0.05, 0.7],
    "line": {"eta": 0.03, "n": 2001, "x_offsets_rt": [2.0, -1.0]},
    "notes": ["区間 1 万 step・保持なし (2026-10-16 ユーザ決定; plan §4.4)。格子点 s ごとに区間 (s − 10,000, s] を 1 回判定し、"
              "最初に合格した格子点で止める。目的の量・流量は区間の中の幅を v1 と同じ許容で、残差は毎 step の系列 (outer_end の全行) で"
              "区間あたりの v1 と同じ閾値、上昇は直前の 1 万 step の区間と比べる。壁の熱流は節点ごとの区間の幅 (相対の参照は区間の終端の値)",
              "最初に判定できる格子点は開始 + 2 × 10,000 + 間隔 (間隔 2,000 なら 22,000)"],
}

INLET_PREFIX = "inlet_"
OUTLET_KINDS = ("outlet_statPress", "outflow")
OTHER_KINDS = ("wall", "wall_isothermal", "slip", "axis")


def _sha256_json(o) -> str:
    return hashlib.sha256(json.dumps(o, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def _get(d, *path, default=None):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def classify_bconds(bcond_cfg: dict) -> dict:
    """bcondConfig の境界を inlet / outlet / other に分ける。判別できない kind (farfield・知らない kind) は ValueError。"""
    out = {"inlet": [], "outlet": [], "other": [], "periodic": []}
    if not isinstance(bcond_cfg, dict) or not bcond_cfg:
        raise ValueError("bcondConfig が空・対応表でない")
    for name, b in bcond_cfg.items():
        if not isinstance(b, dict) or "physID" not in b or "kind" not in b:
            raise ValueError(f"bcondConfig の {name}: physID・kind が無い")
        kind = str(b["kind"])
        pid = int(b["physID"])
        if kind.startswith(INLET_PREFIX):
            out["inlet"].append(pid)
        elif kind in OUTLET_KINDS or kind.startswith("outlet_"):
            out["outlet"].append(pid)
        elif kind in OTHER_KINDS:
            out["other"].append(pid)
        elif kind == "periodic":
            out["periodic"].append(pid)
        else:
            raise ValueError(f"bcondConfig の {name} (physID {pid}) の kind '{kind}' は入口か出口か判別できない "
                             "(farfield は流入と流出の両方がありうる; 問題 YAML で mass_balance を外して inlet・outlet を明示して足す)")
    if not out["inlet"] or not out["outlet"]:
        raise ValueError(f"入口 ({out['inlet']}) か出口 ({out['outlet']}) が無い — 流量の収支を作れない")
    return out


def residual_columns(solver_cfg: dict) -> tuple[list, list]:
    """forge `main.cpp` の residualEquationNames と同じ規則で、必要な残差の列と非活動を許す列を返す。"""
    req = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roe"]
    inactive = ["rms_roUz"]
    turb = _get(solver_cfg, "turbulence") or {}
    model = str(turb.get("model", "none")).lower() if isinstance(turb, dict) else "none"
    if model == "sst":
        req += ["rms_roK", "rms_roOmega"]
    elif model not in ("none", ""):
        raise ValueError(f"turbulence.model '{model}' の残差の列を決められない (sst・none だけ)")
    tr = str(turb.get("transition", "none")).lower() if isinstance(turb, dict) else "none"
    if tr == "lm2009":
        req += ["rms_roGamma", "rms_roReth"]
    elif tr not in ("none", ""):
        raise ValueError(f"turbulence.transition '{tr}' の残差の列を決められない")
    species = _get(solver_cfg, "physProp", "species") or []
    if not isinstance(species, list):
        raise ValueError("physProp.species がリストでない")
    if len(species) >= 2:
        req += [f"rms_roY{s}" for s in range(len(species))]
    if str(_get(solver_cfg, "physProp", "tracer", default="")) == "exhaust":
        req += ["rms_roXi"]
    cond = _get(solver_cfg, "condensation")
    if isinstance(cond, dict) and int(cond.get("condensation", 0)) == 1:
        raise ValueError("凝縮 (condensation 1) のモーメントの残差の列はプリセットで決められない — 問題 YAML の evaluate.settle で "
                         "quantities_remove: [residual] と quantities_add に required の列を明示した residual を書く")
    return req, inactive


def _isothermal_walls(bcond_cfg: dict) -> list:
    walls = []
    for name, b in bcond_cfg.items():
        ints = b.get("ints") or {}
        if b.get("kind") == "wall_isothermal" or (isinstance(ints, dict) and int(ints.get("wallProfile", 0) or 0) == 1):
            walls.append({"name": str(name), "physID": int(b["physID"]), "outputHDFflg": int(b.get("outputHDFflg", 0) or 0)})
    return walls


def criteria_file_name() -> str:
    """run の判定の基準のファイルの名前 (settle.CRITERIA_FILE)。"""
    return _settle().CRITERIA_FILE


def main_stage_segment() -> dict:
    """設計チェーンの NS の run (runner_axismach.run_staged_ns) の判定区間の識別 (`grid.segment`; plan §4.2・§5.1 #6)。
    判定区間 = 本段の forge の起動 (forge の step は起動ごとに 0 から = grid.start_step)。どの起動が本段かは、本段の起動の直前に
    run_staged_ns が書く通知 MAIN_STAGE_STARTED (`main_config.fnv` = 本段の solverConfig.yaml の FNV-1a 64、`settle_criteria.sha256` =
    この基準) で決まり、stage_manifest.json の本段の段 (tag main、forge の後に足す) と forge_launches.jsonl の起動とは cfg_fnv で突き合わせる。
    リミッタの基準値はその段の limiter_ref (forge の limiter_ref.json の写し) にある。基準のファイルは準備の後に書き換えない。"""
    return {"stage": "main", "start_record": "MAIN_STAGE_STARTED", "config_fnv_key": "main_config.fnv",
            "criteria_sha256_key": "settle_criteria.sha256", "stage_manifest": "stage_manifest.json", "stage_tag": "main",
            "match": "cfg_fnv"}


def resolve(problem_cfg_evaluate_settle: dict, prepare_info: dict, solver_cfg: dict, bcond_cfg: dict, segment=None) -> dict:
    """問題 YAML の evaluate.settle を run の設定で解決して settle_criteria.yaml の中身 (dict) を返す。不正・判別できない入力は ValueError
    (settle.CriteriaError は ValueError の子)。segment は `grid.segment` (判定区間の識別; 設計チェーンは `main_stage_segment()`、
    None = 未記録)。"""
    st = _settle()
    es = problem_cfg_evaluate_settle
    if not isinstance(es, dict):
        raise ValueError(f"evaluate.settle が対応表でない ({es!r})")
    unknown = [k for k in es if k not in SETTLE_KEYS]
    if unknown:
        raise ValueError(f"evaluate.settle の知らないキー {unknown} (許されるキー: {list(SETTLE_KEYS)})")
    name = es.get("preset")
    if name not in PRESETS:
        raise ValueError(f"evaluate.settle.preset '{name}' を知らない (プリセット: {list(PRESETS)})")
    P = PRESETS[name]
    th = copy.deepcopy(P["thresholds"])

    # ---- run の設定から
    for k in ("Md", "scale_m"):
        if k not in prepare_info:
            raise ValueError(f"prepare_info に {k} が無い")
    Md = float(prepare_info["Md"])
    S = float(prepare_info["scale_m"])
    ni = _get(prepare_info, "mesh", "ni")
    nj = _get(prepare_info, "mesh", "nj")
    if ni is None or nj is None:
        raise ValueError("prepare_info.mesh に ni・nj が無い")
    dt = _get(solver_cfg, "time", "outStepInterval")
    if not isinstance(dt, int) or isinstance(dt, bool) or dt <= 0:
        raise ValueError(f"solverConfig の time.outStepInterval が正の整数でない ({dt!r})")
    if int(_get(solver_cfg, "time", "unsteady", default=0) or 0) != 0:
        raise ValueError(f"{name} は定常 (time.unsteady 0) の run だけ (非定常の判定は対象外)")
    axisym = int(_get(solver_cfg, "mesh", "isAxisymmetric", default=0) or 0) == 1
    out = _get(solver_cfg, "output") or {}
    req, inactive = residual_columns(solver_cfg)
    bc = classify_bconds(bcond_cfg)

    # ---- 試験部の線
    eta = float(es.get("line_eta", P["line"]["eta"]))
    if "x_range_m" in es:
        xr = es["x_range_m"]
        if not (isinstance(xr, list) and len(xr) == 2):
            raise ValueError(f"evaluate.settle.x_range_m が [x0, x1] でない ({xr!r})")
        x_range = [float(xr[0]), float(xr[1])]
        x_src = "evaluate.settle.x_range_m"
    else:
        dom = _get(prepare_info, "physical_wall", "domain")
        if "x_E" not in prepare_info or not (isinstance(dom, (list, tuple)) and len(dom) == 2):
            raise ValueError("prepare_info に x_E か physical_wall.domain が無い — 試験部の線の x 範囲を evaluate.settle.x_range_m で与える")
        xE = float(prepare_info["x_E"])
        xF = float(dom[1])
        o0, o1 = P["line"]["x_offsets_rt"]
        x_range = [(xE + o0) * S, (xF + o1) * S]
        x_src = (f"既定 (x_E + {o0:g}, x_F {o1:+g}) r_t × scale_m (x_E {xE!r}、x_F = prepare_info.physical_wall.domain[1] {xF!r}; "
                 "nozzle_report.metrics の w_test と同じ範囲)")

    q_res = {"kind": "residual", "name": "residual", "file": "residual_history.csv", "required": req, "inactive_ok": inactive,
             "thresholds": th["residual"]}
    if "residual_series" in P:          # v1 は書かない (settle の既定 grid; 基準の中身は今と同じ)
        q_res["series"] = P["residual_series"]
    q_mass = {"kind": "mass_balance", "name": "mass_balance", "file": "boundary_massflux.csv", "inlet": sorted(bc["inlet"]),
              "outlet": sorted(bc["outlet"]), "thresholds": th["mass_balance"],
              "notes": [f"入口・出口は bcondConfig の kind から (inlet_* / outlet_*・outflow); 記録だけの境界 physID {sorted(bc['other'])}"]}
    tec = dict(th["exit_core_mach"])
    tol_over = tec.pop("tol_over_Md")
    q_core = {"kind": "exit_core_mach", "name": "exit_core_mach", "nodes_file": st.DEFAULT_NODES_FILE, "eta_range": list(P["exit_eta_range"]),
              "thresholds": {"tol": tol_over * Md, **tec}, "notes": [f"tol = {tol_over:g} × M_d (0.1 × 0.05 % × M_d、M_d {Md!r})"]}
    q_line = {"kind": "line_mach", "name": "line_mach", "nodes_file": st.DEFAULT_NODES_FILE, "eta": eta, "x_range_m": x_range,
              "n": int(P["line"]["n"]), "Md": Md, "thresholds": th["line_mach"], "notes": [f"x 範囲: {x_src}"]}
    q_ang = {"kind": "exit_flow_angle", "name": "exit_flow_angle", "nodes_file": st.DEFAULT_NODES_FILE,
             "eta_range": list(P["exit_eta_range"]), "thresholds": th["exit_flow_angle"]}
    quantities = [q_res, q_mass, q_core, q_line, q_ang]
    excluded = []
    walls = _isothermal_walls(bcond_cfg)
    if walls:
        quantities.append({"kind": "wall_heat_nodes", "name": "wall_heat_nodes",
                           "walls": [{"name": w["name"], "physID": w["physID"]} for w in walls], "axisymmetric": axisym,
                           "file_pattern": st.DEFAULT_WALL_PATTERN, "dataset": st.DEFAULT_WALL_DATASET, "thresholds": th["wall_heat_nodes"]})
        if "wall_measure" in P:         # v1 は書かない (settle の既定 terminal_diff)
            quantities[-1]["measure"] = P["wall_measure"]
    else:
        excluded.append({"kind": "wall_heat_nodes", "name": "wall_heat_nodes",
                         "reason": "等温壁が無い (bcondConfig に wall_isothermal も ints.wallProfile も無い) — iface_Qf_eff は等温壁・定常の"
                                   "ときだけ値を持つ (conjugateWall.cpp fillInterfaceDiagnostics)"})

    # ---- 問題 YAML の足し引き・閾値の上書き
    names = {q["name"]: q for q in quantities}
    for nm in es.get("quantities_remove", []) or []:
        if nm not in names:
            raise ValueError(f"evaluate.settle.quantities_remove: 量 '{nm}' がプリセットに無い (量: {list(names)})")
        q = names.pop(nm)
        quantities.remove(q)
        excluded.append({"kind": q["kind"], "name": nm, "reason": "問題 YAML の evaluate.settle.quantities_remove で外した"})
    ov = es.get("thresholds", {}) or {}
    if not isinstance(ov, dict):
        raise ValueError("evaluate.settle.thresholds が対応表でない")
    for nm, kv in ov.items():
        if nm not in names:
            raise ValueError(f"evaluate.settle.thresholds: 量 '{nm}' が無い (量: {list(names)})")
        if not isinstance(kv, dict):
            raise ValueError(f"evaluate.settle.thresholds.{nm} が対応表でない")
        allowed = st.THRESHOLD_KEYS[names[nm]["kind"]] + ("provisional",)
        bad = [k for k in kv if k not in allowed]
        if bad:
            raise ValueError(f"evaluate.settle.thresholds.{nm}: 知らないキー {bad} (許されるキー: {list(allowed)})")
        names[nm]["thresholds"] = {**names[nm]["thresholds"], **kv}
        names[nm].setdefault("notes", []).append(f"閾値を問題 YAML で上書き: {kv}")
    for q in es.get("quantities_add", []) or []:
        quantities.append(copy.deepcopy(q))

    # 節点の座標は load_field と同じ格子の h5 から読む (float のビルドは res_nodes の座標が float32; settle.py の docstring)
    structured = {"ni": int(ni), "nj": int(nj), "scale_m": S}
    mesh_file = _get(solver_cfg, "mesh", "meshFileName")
    if not (isinstance(mesh_file, str) and mesh_file):
        raise ValueError("solverConfig の mesh.meshFileName が無い (節点の座標を読む格子の h5 を決められない)")
    structured["mesh_file"] = mesh_file
    pre = [
        {"quantity": "mass_balance", "requires": "solverConfig output.boundaryMassflux: 1",
         "satisfied": int(out.get("boundaryMassflux", 0) or 0) == 1},
        {"quantity": "exit_core_mach・line_mach・exit_flow_angle", "requires": "solverConfig output.nodeProbeFile (節点の表; probe_node_ids)",
         "satisfied": bool(out.get("nodeProbeFile"))},
    ]
    if walls:
        pre.append({"quantity": "wall_heat_nodes", "requires": "solverConfig output.interfaceDiag: 1",
                    "satisfied": int(out.get("interfaceDiag", 0) or 0) == 1})
        pre.append({"quantity": "wall_heat_nodes", "requires": "壁の bcond の outputHDFflg: 1 (" + ", ".join(w["name"] for w in walls) + ")",
                    "satisfied": all(w["outputHDFflg"] == 1 for w in walls)})
    crit = {
        "schema": st.SCHEMA,
        "source": {"preset": name, "preset_version": P["version"], "preset_sha256": _sha256_json(P),
                   "evaluate_settle": es, "evaluate_settle_sha256": _sha256_json(es),
                   "generator": "design/forge_design/evaluate/settle_presets.py",
                   "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "structured": structured,
        "grid": {"out_step_interval": int(dt), "start_step": 0, "segment": segment},
        "window": dict(P["window"]),
        "hold": int(P["hold"]),
        "quantities": quantities,
        "excluded": excluded,
        "preconditions": pre,
        "notes": ["grid.start_step は本段の forge の起動 (step 0)。segment は判定区間の識別 (設計チェーンでは本段の開始の通知 "
                  "MAIN_STAGE_STARTED を指す; main_stage_segment、plan §5.1 #6)。None = 未記録",
                  "閾値は暫定 (provisional: true; plan §6 V3 で決める)"] + list(P.get("notes", [])),
    }
    st.validate_criteria(crit)       # 誤用はここで CriteriaError (ValueError)
    return crit


def write_settle_criteria(run_dir, problem_cfg_evaluate_settle: dict, prepare_info: dict, solver_cfg: dict, bcond_cfg: dict,
                          segment=None) -> dict:
    """resolve して run の settle_criteria.yaml を書く (一時ファイル → rename)。出典に run の入力のファイルの sha256 を足す。
    戻り {'path', 'sha256', 'criteria'}。"""
    st = _settle()
    run = Path(run_dir)
    crit = resolve(problem_cfg_evaluate_settle, prepare_info, solver_cfg, bcond_cfg, segment=segment)
    inputs = {}
    for fn in ("prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml"):
        p = run / fn
        if p.is_file():
            inputs[fn] = hashlib.sha256(p.read_bytes()).hexdigest()
    crit["source"]["inputs_sha256"] = inputs
    crit["source"]["written"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    head = ("# settle_criteria.yaml — design/forge_design/evaluate/settle_presets.py が問題 YAML の evaluate.settle から作った\n"
            f"# (プリセット {crit['source']['preset_version']})。判定: python3 solver_density_cuda/tools/check_settled.py <run_dir>\n"
            "# 手で直さない (直すなら問題 YAML の evaluate.settle を直して作り直す)。\n")
    text = head + yaml.safe_dump(crit, allow_unicode=True, sort_keys=False, width=200)
    st.validate_criteria(st._yaml_load_strict(text))      # 往復して読めること
    p = run / st.CRITERIA_FILE
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, p)
    return {"path": str(p), "sha256": hashlib.sha256(text.encode()).hexdigest(), "criteria": crit}


# ------------------------------------------------------------------ 節点の表
def mesh_xr(nozzle_h5, ni: int, nj: int, scale_m: float):
    """格子の h5 から X, R [r_t 単位] (ni, nj) を作る (`nozzle_report.load_field` と同じ: COORD / scale_m を COORD の型のまま)。"""
    import h5py
    with h5py.File(nozzle_h5, "r") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    if nc.shape[0] != ni * nj:
        raise ValueError(f"{nozzle_h5}: 節点数 {nc.shape[0]} が ni × nj = {ni * nj} でない")
    S = float(scale_m)
    return (nc[:, 0] / S).reshape(ni, nj), (nc[:, 1] / S).reshape(ni, nj)


def probe_node_ids(ni: int, nj: int, X, R, criteria: dict) -> np.ndarray:
    """節点の量 (exit_core_mach・line_mach・exit_flow_angle) に要る節点の和集合 (昇順の構造格子の番号 i·nj + j)。
    X, R は `mesh_xr` の (ni, nj) [r_t]。series_csv の量は節点を使わない。"""
    X = np.asarray(X)
    R = np.asarray(R)
    if X.shape != (ni, nj) or R.shape != (ni, nj):
        raise ValueError(f"X, R の形 {X.shape}・{R.shape} が (ni, nj) = ({ni}, {nj}) でない")
    S = float(criteria["structured"]["scale_m"])
    ids = set()
    exit_needed = False
    for q in criteria["quantities"]:
        if "series_csv" in q:
            continue
        if q["kind"] in ("exit_core_mach", "exit_flow_angle"):
            exit_needed = True
        elif q["kind"] == "line_mach":
            exit_needed = True          # 出口の列は全部 (線が使わなくても入れる)
            xq = np.linspace(q["x_range_m"][0] / S, q["x_range_m"][1] / S, q["n"])
            xa = X[:, 0]
            i = np.clip(np.searchsorted(xa, xq) - 1, 0, ni - 2)
            cols = set(np.unique(i).tolist()) | set((np.unique(i) + 1).tolist())
            cols |= {c - 1 for c in cols if c > 0} | {c + 1 for c in cols if c < ni - 1}   # 両隣の列も (丸めの差の余裕)
            eta = float(q["eta"])
            for c in sorted(cols):
                xp = R[c] / R[c, -1]
                k = int(np.searchsorted(xp, eta, side="right")) - 1
                js = {0, nj - 1}
                for j in range(k - 1, k + 3):        # 挟む 2 節点 (k, k+1) とその両側
                    if 0 <= j < nj:
                        js.add(j)
                ids.update(c * nj + j for j in js)
    if exit_needed:
        ids.update((ni - 1) * nj + j for j in range(nj))
    return np.asarray(sorted(ids), dtype=np.int64)


def write_node_probe_file(path, ids, header: str | None = None) -> Path:
    """forge の output.nodeProbeFile の書式で書く: 1 行 1 節点番号 (10 進の整数)、`#` 以降はコメント。重複は ValueError。"""
    ids = [int(x) for x in ids]
    if len(set(ids)) != len(ids):
        raise ValueError("節点の番号に重複")
    if any(x < 0 for x in ids):
        raise ValueError("負の節点の番号")
    p = Path(path)
    lines = []
    if header:
        lines += ["# " + ln for ln in header.splitlines()]
    lines += [str(x) for x in ids]
    p.write_text("\n".join(lines) + "\n")
    return p


def read_node_probe_file(path) -> list:
    """forge の readNodeList と同じ規則で読む (試験用): `#` 以降と空行は無視、1 行 1 整数、重複は ValueError。"""
    out = []
    seen = set()
    for ln in Path(path).read_text().splitlines():
        t = ln.split("#", 1)[0].strip()
        if not t:
            continue
        v = int(t, 10)
        if v in seen:
            raise ValueError(f"重複 {v}")
        seen.add(v)
        out.append(v)
    return out
