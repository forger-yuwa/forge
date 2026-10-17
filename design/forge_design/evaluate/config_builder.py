"""段の solverConfig の仕上げと段の宣言 (plan tooling-generator-defaults-and-launch-gate §4.2・§4.3 の 4 (a)・4b・§5.1 #11)。

生成器 (runner_wt・runner_walldriven・runner_axismach の Euler/NS・runner の TP の段階起動・opt/driver の CPG の評価・
runner_sern の 2D/3D (driver_sern を含む)・feedback/euler_loop) は、forge に渡す段の config について次の 2 つを行う。

1. `finish_config` (段の config を書く前): 生成器の機種 (`FAMILIES` の family) の候補プロファイルから、config・bcond の適用条件
   (`recipe_settings.derive_applies`) に合うプロファイルを選ぶ。合えば、明示すべきキー (`recommended_settings.yaml` の
   `explicit_keys` の基準、`recipe_settings.explicit_for`) のうち config が書いていないものを**旧い実効値** (書き写しに使った
   ソルバの既定 = 今の実効値) で足す。書き換えは `yaml_strict.insert_scalars` (既にある写像へのキーの挿入だけ。既存の値・書式は
   1 文字も変えない)。合うプロファイルの無い config (風洞 NS v1・CPG の axis-Mach・断熱壁の axis-Mach NS など、書き写していない
   機種) は変えない。**数値の既定は変えない** (足すのは今の実効値だけ; 生成器の段の作り方・値は従来のまま)。
2. `declare_stage` (段の config を書いた後・forge の起動の前): 段の宣言 `RECIPE.json` を書く (`recipe_settings.declare`; 宣言の
   config の sha256 = 今の solverConfig.yaml)。段の tag ごとの写し `RECIPE_<tag>.json` も残す。
   **問題 YAML が明示した数値の上書き** (CFL の梯子・`cfl_main` など、プロファイルの `fixed` と違う値; plan §4.3 の 4b) は
   理由つきの逸脱に変換する: プロファイルとの照合で FAIL したキーのうち、生成器がそのキーを決めるのに読む問題 YAML のキー
   (`PROBLEM_SOURCES`) が問題 YAML に書いてあるもの・生成器が受け取った引数で決まったもの (`args`) を、出所 (問題 YAML のパスと
   キーと値・引数) つきの逸脱にする。問題 YAML に書いていないキーは**生成器の既定** (`GENERATOR_DEFAULTS`・`GENERATOR_WRITES`、
   出所 `generator_default`: 関数名と値。driver_sern の `opt.soft_steps` 1500、gas.transport の無い SERN の viscMethod 1 など;
   plan §5.1 #27) も出所にする。理由の 1 行は共用 (`DEVIATION_REASON`)。それでも通らない (出所の分からない違い・禁止キー・
   未知のキー・書き写していない段) ときは宣言を書かず、`RECIPE_<tag>.FAILED.json` に FAIL を残す (関門は宣言の無い run と
   して扱う: 移行期は WARN、強制の後は FAIL)。合うプロファイルの無い config は宣言しない (古い宣言は config の sha256 が
   違えば消す)。

宣言の照合は投入の関門 (`solver_density_cuda/tools/launch_gate.py`) と同じ `recipe_settings.check_config`。

段の値の書き換え (soft の CFL・次数・段の長さなど) は `set_scalars` (`yaml_strict.replace_scalars`: YAML 上の位置で値トークンだけを
置換; 正規表現・文字列置換の置き換え、書く文字列は従来と同じ)。plan §4.2 の `build_solver_config(profile, stage, problem)`
(プロファイルから config を組む) の形にはしていない: 段取り 1 のプロファイルは生成器の出力の書き写し (`recipes/regen_profiles.py`)
なので、プロファイルから組むと書き写しが循環する。生成器の段の作り方は従来のまま、仕上げと宣言だけをここに集めた (§5.1 #11)。
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

FORGE_ROOT = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
FORGE_TOOLS = FORGE_ROOT / "solver_density_cuda" / "tools"
PLAN = "plans/active/tooling-generator-defaults-and-launch-gate.md"

RECIPE = "RECIPE.json"
CONTEXT = "RECIPE_CONTEXT.json"
DEVIATION_REASON = ("問題 YAML・呼び出しの引数が明示した上書き、または問題 YAML に書いていないときの生成器の既定 "
                    "(生成器が逸脱に自動で変換、出所は source; plan tooling-generator-defaults-and-launch-gate §4.3 の 4b・§5.1 #27)")

# 生成器の機種 → 候補プロファイル (適用条件の合う最初のものを選ぶ)。書き写した機種だけ (plan §5.1 #6・棚卸し §1.1)
FAMILIES = {
    "wt_modef": ("wt_modef_euler_cell",),                       # runner_wt.run_staged (モード F、Euler cell)
    "wt_feedback": ("wt_feedback_euler_node",),                 # feedback/euler_loop (node)
    "walldriven": ("walldriven_euler_cpg",),                    # runner_walldriven.run_staged
    "axismach_euler": ("axismach_euler_tp",),                   # runner_axismach.run_staged
    "axismach_ns": ("axismach_ns_sst_tp",),                     # runner_axismach.prepare_ns → run_staged_ns
    "axismach_ns_b0": ("axismach_ns_sst_tp_b0",),               # 同 (evaluate.main_line の B0)
    "bell_tp": ("bell_tp_ns_sst", "bell_tp_euler"),             # runner.run_staged_tp・continue_tp
    "bell_cpg": ("bell_cpg_ns_sst",),                           # opt/driver の CPG の評価 (soft → 本段)
    "sern2d": ("sern2d_finite_base_frozen_tp_sst",),            # runner_sern.run_staged (2D、driver_sern)
    "sern3d": ("sern3d_frozen_tp_sst",),                        # runner_sern3d → runner_sern.run_staged
}

_CFL = ("time.deltaT.cfl", "time.deltaT.cfl_pseudo")


def _each(keys, srcs) -> dict:
    return {k: tuple(srcs) for k in keys}


# 出所の候補の書き方: "問題 YAML のキー" (値 = config の値のとき出所) か (キー, 検査)。検査:
#   "in"                       問題 YAML の列 (CFL の梯子) に config の値がある
#   "bool01"                   問題 YAML の真偽 → config の 0/1
#   ("legs", 梯子のキー, 床1)  段の長さ: 梯子 (列) があれば int(値 / 脚の数) (床1 なら最小 1)、無ければ int(値) — runner_sern の脚の分け方
#   ("venkat",)                limiter_scaled から決まる venkatK (1 → 0.05、0 → 1.0; runner_sern._solver_config)
# 値の合わない候補は出所にしない (手で変えた値を、たまたま書いてある問題 YAML のキーのせいにしない)。
_AXI_ALL = {"mesh.axisymMethod": ("evaluate.axisym_method",), "mesh.axisRFloor": ("evaluate.axis_r_floor",)}
_SERN_ALL = {"space.limiter": ("evaluate.limiter",), "space.limiterScaled": ("evaluate.limiter_scaled",),
             "space.venkatK": ("evaluate.venkat_k", ("evaluate.limiter_scaled", ("venkat",))),
             "space.slauWallNormalChi": ("mesh.slau_wall_normal_chi",),
             "time.deltaT.implicitRelax": ("evaluate.implicit_relax",),
             "turbulence.wallTreatmentSST": ("evaluate.wall_treatment_sst",),
             "mesh.nodeInletCornerWall": (("mesh.node_inlet_corner_wall", "bool01"),)}
_SERN_LADDER = {"main": _each(_CFL, ["evaluate.cfl_main"]),
                "warm_lam": {**_each(_CFL, ["opt.warm_lam_cfl", ("opt.warm_lam_ramp", "in")]),
                             "time.last.nStepOuter": (("opt.warm_lam_steps", ("legs", "opt.warm_lam_ramp", False)),)},
                "soft": {**_each(_CFL, ["opt.soft_cfl", ("opt.soft_ramp", "in")]),
                         "time.last.nStepOuter": (("opt.soft_steps", ("legs", "opt.soft_ramp", True)),)},
                "mid": {**_each(_CFL, ["opt.soft_cfl", ("opt.soft_ramp", "in")]),
                        "time.last.nStepOuter": (("opt.mid_steps", ("legs", "opt.soft_ramp", True)),)},
                "*": _SERN_ALL}

# family → 段 (`*` は全段) → config のキー → 生成器がそのキーを決めるのに読む問題 YAML のキー (点区切り) と検査。
# 各生成器の該当行: runner._solver_config_tp・tp_stage_plan (cfl_main・stage_steps)、runner_axismach.main・prepare・
# _apply_gas_to_config・_ns_configs・_main_line_spec (cfl_main・implicit_relax・axisym_method・axis_r_floor・main_line)、
# runner_sern._solver_config・main (cfl_main・limiter・limiter_scaled・venkat_k・slau_wall_normal_chi・implicit_relax・
# wall_treatment_sst・node_inlet_corner_wall、段の梯子は opt:)、runner_sern3d.main (opt.soft_conv も渡す)。
# runner_wt・walldriven・euler_loop・opt/driver の CPG は、問題 YAML が段の数値を変えるキー (問題のキー以外) を読まない。
PROBLEM_SOURCES = {
    "bell_tp": {"main": _each(_CFL, ["evaluate.cfl_main"]),
                "soft": {"time.last.nStepOuter": ("evaluate.stage_steps.soft",)},
                "mid": {"time.last.nStepOuter": ("evaluate.stage_steps.mid",)}},
    "axismach_euler": {"main": _each(_CFL, ["evaluate.cfl_main"]), "*": _AXI_ALL},
    "axismach_ns": {"main": _each(_CFL, ["evaluate.cfl_main"]),
                    "*": {**_AXI_ALL, "time.deltaT.implicitRelax": ("evaluate.implicit_relax",)}},
    "axismach_ns_b0": {"main": {**_each(_CFL, ["evaluate.main_line.cfl"]),
                                "time.deltaT.implicitRelax": ("evaluate.main_line.relax", "evaluate.implicit_relax"),
                                "time.deltaT.lineDtDirectionalCap": ("evaluate.main_line.cap",),
                                "time.deltaT.implicitThermalJacobian": ("evaluate.main_line.itj",)},
                       "*": {**_AXI_ALL, "time.deltaT.implicitRelax": ("evaluate.implicit_relax",)}},
    "sern2d": _SERN_LADDER,
    "sern3d": {**_SERN_LADDER, "soft": {**_SERN_LADDER["soft"], "space.convMethod": ("opt.soft_conv",)}},
}

# 生成器の既定 (generator_default; plan §5.1 #27・#11a (3))。問題 YAML にキーが書いていないとき生成器が使う値も逸脱の出所にする
# (driver_sern は `opt.soft_steps` が無ければ 1500 で回すが、書き写しの代表の問題は 4000 を書いていたので、出所が無く宣言が FAILED に
# なった)。出所の記録は {"generator_default": {"function": 既定を持つ関数, "key": 問題 YAML のキー, "value": 既定, ...}}。
# (1) GENERATOR_DEFAULTS: family → 問題 YAML のキー → [(関数, 生成器の名前の集合 (宣言の generator の先頭の語) か None = どれでも,
#     既定)]。既定は「問題 YAML のキーの値」として PROBLEM_SOURCES の同じキーの検査に当てる (梯子の脚の分け方なども同じ)。
#     PROBLEM_SOURCES の候補のキーが問題 YAML に無いときだけ使う (書いてあればそちらが出所)。
# (2) GENERATOR_WRITES: family → 段 (`*` は全段) → config のキー → [(関数, 問題 YAML のキー, 書く値)]: 問題 YAML のキーが無いとき
#     生成器が config に書く値 (None = 書かない; `finish_config` が旧い実効値で明示したキーはその値とも照合する)。
# 各値の該当行: runner_sern.main (opt: の既定)・opt/driver_sern.SernCampaign._eval_op (同)・runner_sern3d.run_staged と main
# (同、--soft-steps の既定 2000)・runner_sern.prepare (cfl_main 4.0)・runner_sern3d.prepare (1.0)・runner_sern._solver_config
# (limiter 2・limiter_scaled 1・venkat_k・wall_treatment_sst 0、gas.transport が無ければ viscMethod 1 + thermCondMethod 1、
# implicit_relax・slau_wall_normal_chi・node_inlet_corner_wall は無ければ書かない)・runner_axismach._euler_config
# (_config_euler_node の cfl 4.0)・runner_axismach._ns_configs (cfl_main 1.0、implicit_relax は無ければ書かない)・
# runner._solver_config_tp (cfl_main 4.0)・runner.TP_STAGE_STEPS (soft・mid 2000)。B0 の main_line の各キーは必須 (既定なし)。
_SERN_SC = "runner_sern._solver_config"


def _sern_defaults(entry: str, ladder: dict, prepare: str, cfl_main: float) -> dict:
    out = {f"opt.{k}": [(entry, {entry}, v)] for k, v in ladder.items()}
    out["evaluate.cfl_main"] = [(prepare, None, cfl_main)]
    out.update({"evaluate.limiter": [(_SERN_SC, None, 2)], "evaluate.limiter_scaled": [(_SERN_SC, None, 1)],
                "evaluate.wall_treatment_sst": [(_SERN_SC, None, 0)]})
    return out


_SERN_MAIN_OPT = {"soft_steps": 3000, "soft_cfl": 0.5, "warm_lam_steps": 0, "warm_lam_cfl": 0.2, "mid_steps": 0}
_SERN_DRV_OPT = {"soft_steps": 1500, "soft_cfl": 0.5, "warm_lam_steps": 0, "warm_lam_cfl": 0.2, "mid_steps": 0}
_SERN3D_OPT = {"soft_steps": 2000, "soft_cfl": 0.5, "soft_conv": 0, "warm_lam_steps": 0, "warm_lam_cfl": 0.2, "mid_steps": 0}
_SERN2D_DEFAULTS = _sern_defaults("runner_sern.main", _SERN_MAIN_OPT, "runner_sern.prepare", 4.0)
for _k, _v in _SERN_DRV_OPT.items():
    _SERN2D_DEFAULTS[f"opt.{_k}"].append(("opt.driver_sern._eval_op", {"opt.driver_sern._eval_op"}, _v))
GENERATOR_DEFAULTS = {
    "sern2d": _SERN2D_DEFAULTS,
    "sern3d": _sern_defaults("runner_sern3d.run_staged", _SERN3D_OPT, "runner_sern3d.prepare", 1.0),
    "axismach_euler": {"evaluate.cfl_main": [("runner_axismach._euler_config", None, 4.0)]},
    "axismach_ns": {"evaluate.cfl_main": [("runner_axismach._ns_configs", None, 1.0)]},
    "bell_tp": {"evaluate.cfl_main": [("runner._solver_config_tp", None, 4.0)],
                "evaluate.stage_steps.soft": [("runner.TP_STAGE_STEPS", None, 2000)],
                "evaluate.stage_steps.mid": [("runner.TP_STAGE_STEPS", None, 2000)]},
}
_SERN_WRITES = {"*": {"physProp.viscMethod": [(_SERN_SC, "gas.transport", 1)],
                      "physProp.thermCondMethod": [(_SERN_SC, "gas.transport", 1)],
                      "time.deltaT.implicitRelax": [(_SERN_SC, "evaluate.implicit_relax", None)],
                      "space.slauWallNormalChi": [(_SERN_SC, "mesh.slau_wall_normal_chi", None)],
                      "mesh.nodeInletCornerWall": [(_SERN_SC, "mesh.node_inlet_corner_wall", None)]}}
GENERATOR_WRITES = {
    "sern2d": _SERN_WRITES,
    "sern3d": _SERN_WRITES,
    "axismach_ns": {"*": {"time.deltaT.implicitRelax": [("runner_axismach._ns_configs", "evaluate.implicit_relax", None)]}},
    "axismach_ns_b0": {"*": {"time.deltaT.implicitRelax": [("runner_axismach._main_line_spec", "evaluate.implicit_relax", None)]}},
}


# ------------------------------------------------------------------------------------------------ 道具

def _tools():
    """recipe_settings・yaml_strict (solver_density_cuda/tools)。"""
    if str(FORGE_TOOLS) not in sys.path:
        sys.path.insert(0, str(FORGE_TOOLS))
    import recipe_settings as RS
    import yaml_strict as YS
    return RS, YS


def _bcond_doc(bcond):
    """bcond (文字列・dict・None) → dict (読めなければ None)。"""
    if bcond is None or isinstance(bcond, dict):
        return bcond
    _, YS = _tools()
    try:
        doc = YS.load(str(bcond))
    except Exception:
        return None
    return doc if isinstance(doc, dict) else None


def _token(v) -> str:
    """明示すべきキーの値 → YAML のトークン (数値はそのまま、文字列は二重引用符)。"""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    raise ValueError(f"明示すべきキーの値 {v!r} をトークンにできない")


def set_scalars(cfg_text: str, updates: dict) -> str:
    """段の config の値を YAML 上の位置で書き換える (`yaml_strict.replace_scalars`: その位置の値トークンだけを置換し、読み直して
    要求値・他の値の不変を検査する)。updates {完全修飾パス: トークン (文字列)}。段の作り方の正規表現・文字列置換の置き換え
    (plan §4.2・§5.1 #11): トークンは従来の書き方と同じ文字列を渡す (書かれる文字列は従来と同じ)。パスが無ければ ValueError
    (従来の置換は黙って当たらなかった)。"""
    _, YS = _tools()
    return YS.replace_scalars(cfg_text, {tuple(k.split(".")): str(v) for k, v in updates.items()})


def stage_of(tag: str) -> str | None:
    """段の tag (soft・S1_soft・S3_main_p・warm_cfl0.1_2・soft_cfl0.5・R1_cfl1 など) → プロファイルの段の名前。分からなければ None。
    (recipes/regen_profiles.py の `_stage_of` と同じ対応。)"""
    if re.fullmatch(r"R\d+_cfl.*", tag):
        return "ramp"
    t = re.sub(r"^S\d+_", "", tag)
    if t.startswith("warm_cfl"):
        return "warm_lam"
    if t.startswith("warm_adapt"):
        return "warm_adapt"
    for s in ("soft", "mid", "main"):
        if t.startswith(s):
            return s
    return None


def _nest(flat: dict) -> dict:
    """{点区切りのキー: 値} → 入れ子の dict。"""
    out = {}
    for k, v in flat.items():
        cur = out
        segs = k.split(".")
        for s_ in segs[:-1]:
            cur = cur.setdefault(s_, {})
        cur[segs[-1]] = v
    return out


def select_profile(rc, family: str, doc: dict, bcond, stage: str | None = None) -> str | None:
    """family の候補のうち、config・bcond の適用条件 (導けるもの: 離散化・軸対称・物理・気体・壁・時間) が合う最初のプロファイル。
    段 stage がプロファイルにあり、その段の `fixed` が適用条件を決めるキーを全部持つなら、その条件は適用条件の値でなく段の `fixed` の
    値が決める区分 (`derive_applies`) で比べる (check_config と同じ扱い: 例 SERN の層流暖機は SST の機種の段だが turbulence.model none を
    fixed に持つ)。fixed の値そのものとの違い (viscMethod 1 と 2 など) は宣言の照合 (check_config) と逸脱の出所が扱う。"""
    RS, _ = _tools()
    if family not in FAMILIES:
        raise ValueError(f"未知の生成器の機種 {family!r} ({sorted(FAMILIES)})")
    got = RS.derive_applies(doc, _bcond_doc(bcond))
    flat = RS.flatten(doc)
    for name in FAMILIES[family]:
        if name not in rc["profiles"]:
            continue
        res = RS.resolve(rc, name)
        fixed = (res["stages"].get(stage) or {}).get("fixed") or {}
        ok = True
        for k, v in res["applies"].items():
            cols = RS.APPLIES_KEYS.get(k)
            if cols and all(c in fixed for c in cols):
                # 段の fixed の値が決める区分と config の区分を比べる (2026-10-15、plan §5.1 #27: 値そのものを比べていたので、
                # gas.transport の無い SERN の viscMethod 1 (プロファイルは 2、どちらも ns_sst) がどのプロファイルにも合わず、
                # 宣言も FAILED の記録も書かなかった。値の違いは check_config の fixed の照合と逸脱の出所で扱う)
                ok &= RS.derive_applies(_nest({c: fixed[c] for c in cols})).get(k) == got.get(k, RS.derive_applies(doc).get(k))
            elif k in got:
                ok &= RS._applies_match(v, got[k])
        if ok:
            return name
    return None


# ------------------------------------------------------------------------------------------------ 1. 段の config の仕上げ

def finish_config(cfg_text: str, bcond, family: str, rc=None, tag: str | None = None) -> tuple[str, dict]:
    """段の config に明示すべきキーを旧い実効値で足す (合うプロファイルがあるときだけ)。戻り (config, 仕上げの記録
    {family, profile, explicit_added: {キー: 値}})。足すものが無ければ config はそのまま。tag は段の tag (適用条件の比べ方に使う;
    `select_profile`)。足すキーは段の config から決まる (`explicit_for` の when) ので段ごとに違いうる。"""
    RS, YS = _tools()
    rc = rc or RS.load_recipes()
    doc, fails = RS.parse_config(cfg_text)
    if fails:
        raise ValueError(f"生成器の config を読めない (明示すべきキーを足せない): {fails}")
    bc = _bcond_doc(bcond)
    prof = select_profile(rc, family, doc, bc, stage_of(tag) if tag else None)
    info = {"family": family, "profile": prof, "explicit_added": {}}
    if prof is None:
        return cfg_text, info
    # 出所の検査 (verify) は書き写しの台本の役目。生成器は今の実効値 (旧い実効値) を書くだけなので、ソルバの既定が後で変わっても
    # 生成器の出力は変わらない (それが明示する理由; plan §4.3 の 4 (a))
    add = RS.explicit_for(rc, cfg_text, bc, RS.explicit_table(rc, verify=False))
    if add:
        cfg_text, _ = YS.insert_scalars(cfg_text, {tuple(k.split(".")): _token(v) for k, v in add.items()})
    info["explicit_added"] = dict(add)
    return cfg_text, info


# ------------------------------------------------------------------------------------------------ 2. 段の宣言

def write_context(run_dir, generator: str, problem=None, args: dict | None = None) -> Path:
    """段階起動を別の呼び出しで回す生成器 (prepare_ns → run_staged_ns) のための記録 `RECIPE_CONTEXT.json`: 生成器、問題 YAML の
    所在、prepare が受け取った引数 (None は書かない)。`declare_stage` が逸脱の出所を決めるのに読む。"""
    p = Path(run_dir) / CONTEXT
    rec = {"generator": generator, "problem": (str(Path(problem).resolve()) if problem else None),
           "arguments": {k: v for k, v in (args or {}).items() if v is not None}, "plan": PLAN}
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    return p


def load_context(run_dir) -> dict:
    p = Path(run_dir) / CONTEXT
    try:
        return json.loads(p.read_text()) if p.is_file() else {}
    except (OSError, ValueError):
        return {}


def _dig(raw, dotted):
    cur = raw
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return False, None
        cur = cur[k]
    return True, cur


def _rel(path) -> str:
    p = Path(path).resolve()
    try:
        return str(p.relative_to(FORGE_ROOT.resolve()))
    except ValueError:
        return str(p)


def _problem_raw(problem):
    if not problem:
        return None
    _, YS = _tools()
    try:
        raw = YS.load(Path(problem).read_text(encoding="utf-8"))
    except Exception:
        import yaml
        try:
            raw = yaml.safe_load(Path(problem).read_text(encoding="utf-8"))
        except Exception:
            return None
    return raw if isinstance(raw, dict) else None


def _source_ok(check, pv, cv, raw) -> tuple[bool, dict]:
    """出所の候補の検査 (PROBLEM_SOURCES の書き方)。戻り (合うか, 記録に足す別のキー {キー: 値})。"""
    RS, _ = _tools()
    try:
        if check == "eq":
            return RS.same(pv, cv), {}
        if check == "in":
            return isinstance(pv, (list, tuple)) and any(RS.same(x, cv) for x in pv), {}
        if check == "bool01":
            return RS.same(1 if int(pv) else 0, cv), {}
        if isinstance(check, tuple) and check[0] == "venkat":
            return RS.same(0.05 if int(pv) == 1 else 1.0, cv), {}
        if isinstance(check, tuple) and check[0] == "legs":
            _, ramp_key, floor1 = check
            ok_r, ramp = _dig(raw, ramp_key)
            if ok_r and isinstance(ramp, (list, tuple)) and ramp:
                n = int(float(pv) / len(ramp))
                return RS.same(max(n, 1) if floor1 else n, cv), {ramp_key: ramp}
            return RS.same(int(float(pv)), cv), {}
    except (TypeError, ValueError):
        return False, {}
    raise ValueError(f"PROBLEM_SOURCES の検査 {check!r} が未知")


def _generator_id(generator) -> str | None:
    """宣言の生成器の名前 ("opt.driver_sern._eval_op (m6_on 標準レシピ)" など) → 関数名 (先頭の語)。"""
    return str(generator).split()[0] if generator and str(generator).split() else None


def _default_source(family, stage, key, cv, raw, cands, generator, explicit) -> dict | None:
    """問題 YAML に書いていないキーについて、生成器の既定 (GENERATOR_DEFAULTS・GENERATOR_WRITES) が config の値 cv を説明するか。
    戻り {"generator_default": {...}} か None。cands は PROBLEM_SOURCES の候補 [(問題 YAML のキー, 検査)]、explicit は finish_config が
    旧い実効値で明示したキー {キー: 値}。"""
    RS, _ = _tools()
    gid = _generator_id(generator)
    for k, chk in cands:
        if _dig(raw or {}, k)[0]:
            continue
        for fn, gens, dflt in GENERATOR_DEFAULTS.get(family, {}).get(k, ()):
            if gens is not None and gid not in gens:
                continue
            good, more = _source_ok(chk, dflt, cv, raw or {})
            if good:
                rec = {"function": fn, "key": k, "value": dflt, "note": "問題 YAML に書いていないので生成器の既定"}
                if more:
                    rec["with"] = more
                return {"generator_default": rec}
    table = GENERATOR_WRITES.get(family, {})
    for fn, k, val in list(table.get(stage, {}).get(key, ())) + list(table.get("*", {}).get(key, ())):
        if _dig(raw or {}, k)[0]:
            continue
        if RS.same(val, cv) if val is not None else (cv is None or (key in (explicit or {}) and RS.same(explicit[key], cv))):
            rec = {"function": fn, "key": k, "value": None, "config_value": val,
                   "note": (f"問題 YAML に {k} が無いので生成器は {key} を {RS.show(val)} と書く" if val is not None else
                            f"問題 YAML に {k} が無いので生成器は {key} を書かない"
                            + (f" (finish_config が旧い実効値 {RS.show(explicit[key])} で明示)" if key in (explicit or {}) else ""))}
            return {"generator_default": rec}
    return None


def _source_for(family, stage, key, cv, raw, problem, args, generator=None, explicit=None):
    """FAIL したキー (config の値 cv、書いていなければ None) の出所。引数 (args[段 or *][キー]、値が config と同じとき) が
    あればそれ (引数は問題 YAML より優先)、無ければ PROBLEM_SOURCES の候補のうち問題 YAML に書いてあって検査に合うもの
    (キーと値)、それも無ければ問題 YAML に書いていない候補についての生成器の既定 (`_default_source`; plan §5.1 #27)。
    どれも無ければ None。"""
    RS, _ = _tools()
    arg = ((args or {}).get(stage) or {}).get(key) or ((args or {}).get("*") or {}).get(key)
    if arg and ("value" not in arg or RS.same(arg["value"], cv)):
        return {"argument": arg}
    table = PROBLEM_SOURCES.get(family, {})
    cands, seen = [], set()
    for c in list(table.get(stage, {}).get(key, ())) + list(table.get("*", {}).get(key, ())):
        k, chk = (c, "eq") if isinstance(c, str) else c
        if k not in seen:
            seen.add(k)
            cands.append((k, chk))
    present = {}
    for k, chk in cands:
        ok, v = _dig(raw or {}, k)
        if not ok:
            continue
        good, more = _source_ok(chk, v, cv, raw or {})
        if good:
            present[k] = v
            present.update(more)
    if not present:
        dsrc = _default_source(family, stage, key, cv, raw, cands, generator, explicit)
        if dsrc is not None and problem:
            dsrc["problem"] = _rel(problem)
        return dsrc
    return {"problem": _rel(problem), "keys": present}


def _dev_token(v) -> str:
    return "~" if v is None else json.dumps(v, ensure_ascii=False)


def _unlink(*paths):
    for p in paths:
        try:
            Path(p).unlink()
        except FileNotFoundError:
            pass


def declare_stage(run_dir, tag: str, family: str, *, finished: dict | None = None, problem=None, args: dict | None = None,
                  generator: str | None = None, note: str | None = None, rc=None) -> dict:
    """段の config (run_dir/solverConfig.yaml) を書いた後に呼ぶ。合うプロファイルと段があれば宣言 RECIPE.json (と写し
    RECIPE_<tag>.json) を書く。戻り {declared, profile, stage, deviations, fails}。
    problem: 問題 YAML (逸脱の出所)。無ければ run の RECIPE_CONTEXT.json・prepare_info.json の problem。
    args: 引数による上書き {段 or '*': {config のキー: 出所の説明 (dict)}}。無ければ RECIPE_CONTEXT.json の arguments から
    (axis-Mach NS の cfl_main・implicit_relax)。finished: finish_config の記録 (宣言に `explicit_added` として写す)。"""
    RS, _ = _tools()
    rc = rc or RS.load_recipes()
    run_dir = Path(run_dir)
    cfg_p, bc_p = run_dir / "solverConfig.yaml", run_dir / "bcondConfig.yaml"
    cfg_raw = cfg_p.read_bytes()
    cfg = cfg_raw.decode("utf-8", errors="replace")
    bc = _bcond_doc(bc_p.read_text()) if bc_p.exists() else None
    stage = stage_of(tag)
    _unlink(run_dir / f"RECIPE_{tag}.json", run_dir / f"RECIPE_{tag}.FAILED.json")
    decl_p = run_dir / RECIPE
    out = {"declared": False, "profile": None, "stage": stage, "deviations": [], "fails": []}
    doc, pf = RS.parse_config(cfg)
    prof = select_profile(rc, family, doc, bc, stage) if doc is not None else None
    if prof is None:
        # 書き写していない機種: 宣言しない。前の段の宣言 (config の sha256 が違う) は関門で FAIL になるので消す
        if decl_p.exists():
            try:
                import hashlib
                old = json.loads(decl_p.read_text())
                keep = (old.get("config") or {}).get("sha256") == hashlib.sha256(cfg_raw).hexdigest()
            except (OSError, ValueError):
                keep = False
            if not keep:
                _unlink(decl_p)
        out["fails"] = pf or [("(profile)", f"生成器の機種 {family} の候補 {list(FAMILIES[family])} に適用条件の合うプロファイルが無い")]
        return out
    _unlink(decl_p)
    out["profile"] = prof
    ctx = load_context(run_dir)
    if problem is None:
        problem = ctx.get("problem")
    if problem is None and (run_dir / "prepare_info.json").is_file():
        try:
            problem = json.loads((run_dir / "prepare_info.json").read_text()).get("problem")
        except (OSError, ValueError):
            problem = None
    if args is None:
        args = _context_args(family, ctx)
    extra = {"generator": {"name": generator, "family": family, "tag": tag, "note": note,
                           "problem": (_rel(problem) if problem else None)},
             "explicit_added": dict((finished or {}).get("explicit_added") or {})}
    res = RS.resolve(rc, prof)
    if stage not in res["stages"]:
        out["fails"] = [("(stage)", f"段 {tag} ({stage}) はプロファイル {prof} に無い (書き写していない段: {list(res['stages'])})")]
        _write_failed(run_dir, tag, prof, stage, out["fails"], extra)
        return out
    r0 = RS.check_config(rc, cfg, prof, stage, {}, bc)
    devs, srcs = {}, {}
    if r0["fails"]:
        raw = _problem_raw(problem)
        flat = RS.flatten(doc)
        for k, _m in r0["fails"]:
            if k in devs or k.startswith("(") or k.startswith("applies."):
                continue
            src = _source_for(family, stage, k, flat.get(k), raw, problem, args, generator=generator,
                              explicit=extra["explicit_added"])
            if src is not None:
                devs[k], srcs[k] = flat.get(k), src
    code, r, path = RS.declare(run_dir, prof, stage, [f"{k}={_dev_token(v)}" for k, v in devs.items()],
                               DEVIATION_REASON if devs else None, recipes=rc, deviation_sources=srcs, extra=extra)
    if code != 0:
        out["fails"] = list(r["fails"])
        _write_failed(run_dir, tag, prof, stage, out["fails"], extra)
        return out
    (run_dir / f"RECIPE_{tag}.json").write_bytes(Path(path).read_bytes())
    out.update(declared=True, deviations=r["record"]["deviations"])
    return out


def _context_args(family: str, ctx: dict) -> dict:
    """RECIPE_CONTEXT.json の arguments → declare_stage の args (prepare_ns の cfl_main・implicit_relax)。"""
    a = ctx.get("arguments") or {}
    gen = ctx.get("generator") or "prepare"
    out = {}
    if a.get("cfl_main") is not None:
        out.setdefault("main", {}).update({k: {"name": f"{gen} の引数 cfl_main", "value": a["cfl_main"]} for k in _CFL})
    if a.get("implicit_relax") is not None:
        out.setdefault("*", {})["time.deltaT.implicitRelax"] = {"name": f"{gen} の引数 implicit_relax", "value": a["implicit_relax"]}
    return out


def _write_failed(run_dir: Path, tag, prof, stage, fails, extra):
    rec = {"declared": False, "profile": prof, "stage": stage, "tag": tag,
           "fails": [{"key": k, "message": m} for k, m in fails], **extra,
           "note": "宣言を書かなかった (関門は宣言の無い run として扱う: 移行期は WARN、強制の後は FAIL)。"
                   "問題 YAML の上書きなら PROBLEM_SOURCES に出所を足す、手で変えたなら recipe_settings.py declare で宣言する",
           "plan": PLAN}
    (run_dir / f"RECIPE_{tag}.FAILED.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(f"[recipe] {run_dir.name}/{tag}: 宣言を書かなかった (プロファイル {prof} 段 {stage}): "
          + "; ".join(f"{k}: {m}" for k, m in list(fails)[:3]), flush=True)


def launch_stage(run_dir, cfg_text: str, bcond, tag: str, family: str | None, **kw) -> str:
    """段の config を仕上げて書き、宣言する (forge の起動の直前に呼ぶ)。family が None なら config をそのまま書くだけ。
    戻り: 書いた config (stage_manifest.json に記録する文字列)。kw は declare_stage へ (problem・args・generator・note)。"""
    run_dir = Path(run_dir)
    if family is None:
        (run_dir / "solverConfig.yaml").write_text(cfg_text)
        return cfg_text
    RS, _ = _tools()
    rc = RS.load_recipes()
    cfg_text, fin = finish_config(cfg_text, bcond, family, rc=rc, tag=tag)
    (run_dir / "solverConfig.yaml").write_text(cfg_text)
    declare_stage(run_dir, tag, family, finished=fin, rc=rc, **kw)
    return cfg_text
