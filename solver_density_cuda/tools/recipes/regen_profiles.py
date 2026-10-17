#!/usr/bin/env python3
"""機種ごとのプロファイルの書き写し (plan tooling-generator-defaults-and-launch-gate §4.0 段取り 1・§5.1 #6)。

代表の問題で各生成器を実際に呼び、**起動の全段の solverConfig** を作って、そこから各段の `fixed` を書き写す (手で値を写さない)。
段の config は、各 runner の段階起動の関数 (`run_staged`・`run_staged_ns`・`run_staged_tp`・opt/driver の CPG 評価) を偽の forge で
回して、forge に渡る直前の文字列を取る。本段の config は各 runner の `prepare` が書く文字列を同じ関数の組み合わせで作る
(`prepare` はメッシュ・変換器・設計の MOC を回すので呼ばない。組み合わせは `prepare` の該当行を写した最小限: 下の `_fake_prepare_*`)。

書き写す値: 各段の config を完全修飾パスに平らにし、`problem_keys` (問題 YAML の物理条件・入出力) と `steady_base` の `ignore`
(定常では効かない dt) を除いたもの。Euler (viscMethod 0) の visc・thermCond は 0 でなければならないので残す。
明示すべきキー (plan §5.1 #16): YAML の `explicit_keys` の基準に当たるのに段の config が書いていないキーを、その段の `explicit` に
旧い実効値 (書き写しに使ったソルバの既定) で載せる (`recipe_settings.explicit_for`)。旧い実効値の出所は今のソースで確かめ、
合わなければ止まる (`recipe_settings.explicit_table`)。照合では `explicit` も `fixed` と同じく省略が FAIL。
#11 (2026-10-15) から生成器が段の config に明示すべきキーを旧い実効値で書く (`forge_design.evaluate.config_builder`) ので、
書き写すと `explicit` は空になり、キーは `fixed` に入る (値は #16 の `explicit` と同じ; design/tests/run_generator_recipe_tests.py)。
`explicit` は生成器がまた省略を始めたとき (または新しい明示すべきキーを足したとき) に表に出る。
§1 との差 (`differs_from_recommended`) は `recipe_settings.compute_differs` で共通プロファイルと比べて作る。

    design/.venv-opt/bin/python solver_density_cuda/tools/recipes/regen_profiles.py            # 生成区間を出す
    design/.venv-opt/bin/python solver_density_cuda/tools/recipes/regen_profiles.py --write    # recommended_settings.yaml を書き換える
    design/.venv-opt/bin/python solver_density_cuda/tools/recipes/regen_profiles.py --check    # YAML と一致するか (違えば終了 1)

forge・変換器・AWS は使わない。

代表の問題 YAML は**固定の写し** `recipes/fixtures/<機種>.yaml` を読む (plan §5.1 #17)。case の問題 YAML を直接読むと、case 側の
編集 (別の作業の比較・較正のやり直し) で書き写しの結果が黙って変わり、推奨値の変更が承認を経ずに入る。写しの元のパスと sha256 は
`fixtures/MANIFEST.yaml` にあり、`SPECS` の `problem` (= 元のパス、プロファイルの `transcribed_from.problem`) と食い違うか、
写しの sha256 が記録と違えば (写しを手で書き換えた) 止まる。
写しの中の相対パス (`gas.exhaust`・`geometry.rao_wall_file` など、問題 YAML の所在を基準に解決するもの) は**元の問題 YAML の所在**で
解決する (`materialize`: 作業ディレクトリに置く写しで、そのキーを元の所在の絶対パスに直す。既知のキーの外で元の所在に実在する
ファイルの相対パスがあれば止まる)。

    design/.venv-opt/bin/python solver_density_cuda/tools/recipes/regen_profiles.py --check-fixtures   # 写しと元の違い (WARN、終了 0)

**写しの更新の手順** (写しの更新は推奨値の変更の承認ではない):
1. `--check-fixtures` で元と違う写しを見る (WARN)。違いは case 側の編集であって、それだけでは書き写しの結果は何も変わらない。
2. 代表を新しい元に替えると決めたら `--refresh-fixture <機種> ...` (または `all`) で写しと MANIFEST を書き直す (元のバイト列の写し)。
3. `--check` を回す。一致すれば推奨値は変わらないので、写しの更新だけを commit してよい。
4. 不一致なら推奨値 (生成区間) が変わる。これは推奨値の変更なので、差を plan に書いて上位の判断を得てから `--write` する
   (判断を得るまで写しの更新も commit しない。戻すときは `git checkout` で写しと MANIFEST を戻す)。
"""
import argparse
import contextlib
import datetime
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent
ROOT = TOOLS.parents[1]
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(ROOT / "design"))

import yaml  # noqa: E402

import recipe_settings as RS  # noqa: E402

YAML_PATH = HERE / "recommended_settings.yaml"
FIXTURE_DIR = HERE / "fixtures"
MANIFEST_NAME = "MANIFEST.yaml"
# 問題 YAML の所在を基準に解決する相対パスのキーのうち、opt/driver の PROBLEM_PATH_KEYS (absolutize_problem_paths) に無いもの
EXTRA_PATH_KEYS = (("geometry", "initial_line_run"),)   # runner_axismach.py:447
BEGIN ="  # >>> BEGIN GENERATED: recipes/regen_profiles.py (手で編集しない) >>>"
END = "  # <<< END GENERATED <<<"
PLAN = "plans/active/tooling-generator-defaults-and-launch-gate.md"
DATE = "2026-10-14"
KEEP_FOR_EULER = ("physProp.visc", "physProp.thermCond")   # viscMethod 0 は定数粘性: Euler は 0 でなければならない


# ------------------------------------------------------------------------------------------------ 偽の forge と差し替え

class Capture:
    """forge の代わり: 起動ごとに solverConfig.yaml を控え、nStepOuter の res (中身は空でない印) と残差履歴を書く。"""

    def __init__(self):
        self.cfgs = []

    def __call__(self, run_dir):
        run_dir = Path(run_dir)
        cfg = (run_dir / "solverConfig.yaml").read_text()
        self.cfgs.append(cfg)
        n = int(re.search(r"nStepOuter:\s*(\d+)", cfg).group(1))
        (run_dir / f"res_{n}.h5").write_bytes(b"regen_profiles placeholder")
        (run_dir / "residual_history.csv").write_text(f"step,rms_ro\n0,1\n{n},0.001\n")
        (run_dir / "CONVERGENCE_VERDICT.txt").write_text("VERDICT: PASS (regen_profiles placeholder)\n")
        return 0


@contextlib.contextmanager
def patched(pairs):
    """[(module, 名前, 置き換え)] を差し替えて、抜けるときに戻す。"""
    saved = [(m, n, getattr(m, n)) for m, n, _ in pairs]
    try:
        for m, n, v in pairs:
            setattr(m, n, v)
        yield
    finally:
        for m, n, v in saved:
            setattr(m, n, v)


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _stage_of(tag):
    """manifest の段の tag (soft・S1_soft・warm_cfl0.1_2・R1_cfl1 など) → 段の名前。"""
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
    raise ValueError(f"段の tag {tag!r} を段の名前に写せない")


def _manifest_stages(rd, cfgs):
    """run の stage_manifest.json の段と、偽の forge が控えた config を順に対応させる (cfg_fnv で確かめる)。"""
    from stage_manifest import fnv1a64
    man = json.loads((Path(rd) / "stage_manifest.json").read_text())["stages"]
    if len(man) != len(cfgs):
        raise RuntimeError(f"{rd}: manifest の段 {len(man)} と起動 {len(cfgs)} の数が違う")
    out = []
    for st, cfg in zip(man, cfgs):
        if st["cfg_fnv"] != fnv1a64(cfg):
            raise RuntimeError(f"{rd}: 段 {st['tag']} の config が manifest の記録と違う")
        out.append((_stage_of(st["tag"]), cfg))
    return out


# ------------------------------------------------------------------------------------------------ 生成器ごとの起動

def _common_patches():
    from forge_design.evaluate import runner_axismach as A
    return [(A, "stage_gate", lambda res, cfg: []), (A, "_restart_same_mesh", lambda res, mesh, **k: None),
            (A, "next_stage_init_problems", lambda *a, **k: ([], {}))]


def run_bell_cpg(problem, work):
    """③ベル CPG: opt/driver の CPG 評価 (soft を別ディレクトリ `<tag>_stgA` で → 本段)。"""
    from forge_design.evaluate import runner as R
    from forge_design.opt import driver as DRV
    from forge_design.probdef import load_problem
    cap = Capture()

    def fake_prepare(prob, run_dir, nsteps=None, warm_from=None, stage=None):
        # runner.prepare の CPG の行 (n・出力の間隔・`_solver_config`・`_bcond_config`) と同じ
        p = load_problem(prob)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
        (run_dir / "solverConfig.yaml").write_text(R._solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
        (run_dir / "bcondConfig.yaml").write_text(R._bcond_config(p))
        return {}

    def fake_continue(prob, src_run, new_run, nsteps=None):
        # runner.continue_cpg の config の行 (n・出力の間隔・`_solver_config`) と同じ。bcond は継続元から写す (格子・restart は無し)
        p = load_problem(prob)
        new_run = Path(new_run)
        new_run.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
        (new_run / "solverConfig.yaml").write_text(R._solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
        (new_run / "bcondConfig.yaml").write_text((Path(src_run) / "bcondConfig.yaml").read_text())
        return {}

    raw = yaml.safe_load(Path(problem).read_text())
    camp = DRV.Campaign(problem, work / "campaign", seed_res=str(work / "seed_res.h5"))
    x = [float(raw["dv"][k]["value"]) for k in camp.dv_order]
    fake_collect = lambda prob, rd: {"eta_cf": 0.97, "mdot_kg_s": camp.mdot_1d, "thrust_N": 1.0,
                                     "res_file": sorted(Path(rd).glob("res_*.h5"))[-1].name}
    with patched([(R, "prepare", fake_prepare), (R, "continue_cpg", fake_continue), (R, "run_forge", cap), (R, "collect", fake_collect),
                  (DRV.Campaign, "_diverged", lambda self, rd, out: False),
                  (DRV.Campaign, "_quasisteady", lambda self, rd: "ALL STEADY"),
                  (DRV.Campaign, "_warm_for", lambda self, x: work / "seed_res.h5")]):
        _quiet(camp.evaluate, x, "regen")
    if len(cap.cfgs) != 2:
        raise RuntimeError(f"bell CPG: 起動が {len(cap.cfgs)} 回 (soft と本段の 2 回のはず)")
    bc = yaml.safe_load((work / "campaign" / "regen" / "bcondConfig.yaml").read_text())
    return [("soft", cap.cfgs[0]), ("main", cap.cfgs[1])], bc


def run_bell_tp(problem, work):
    """③ベル TP: runner.run_staged_tp (`--staged`)。"""
    from forge_design.evaluate import runner as R
    from forge_design.probdef import load_problem
    cap = Capture()
    rd = work / "run"

    def fake_prepare(prob, run_dir, nsteps=None, warm_from=None, stage=None):
        p = load_problem(prob)
        Path(run_dir).mkdir(parents=True, exist_ok=False)
        (Path(run_dir) / "bcondConfig.yaml").write_text(R._bcond_config_tp(p))
        (Path(run_dir) / "nozzle.h5").write_bytes(b"x")
        return {}

    with patched([(R, "prepare", fake_prepare), (R, "run_forge", cap)] + _common_patches()):
        _quiet(R.run_staged_tp, problem, rd)
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def run_wt_modef(problem, work):
    """① 風洞 モード F: runner_wt.main (prepare(euler=True) → run_staged)。"""
    from forge_design.evaluate import runner_wt as W
    from forge_design.probdef import load_problem
    cap = Capture()
    rd = work / "run"

    def fake_prepare(problem_path, run_dir, euler=True, nsteps=None, ic_from=None):
        # runner_wt.prepare の n・out_int・`_config_euler`/`_config_sst_node`・`_bcond` と同じ
        p = load_problem(problem_path)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
        out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
        cfg = W._config_euler(p, n, out_int, 4.0, 1) if euler else W._config_sst_node(p, n, out_int, 4.0)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        (run_dir / "bcondConfig.yaml").write_text(W._bcond(p, euler))
        return {}

    with patched([(W, "prepare", fake_prepare), (W, "run_forge", cap), (W, "collect", lambda *a: {})] + _common_patches()):
        _quiet(W.main, [str(problem), str(rd)])
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def run_wt_feedback(problem, work):
    """① 風洞 B10 の帰還 (feedback/euler_loop): パスごとに本段 1 段。cold のパスの config (`euler_loop.pass_config`、n = max(24000, 12000))。"""
    from forge_design.evaluate import runner_wt as W
    from forge_design.feedback import euler_loop as EL
    from forge_design.probdef import load_problem
    p = load_problem(problem)
    cfg, _ = EL.pass_config(p, max(24000, 12000))
    return [("main", cfg)], yaml.safe_load(W._bcond(p, euler=True))


def run_axismach_euler(problem, work):
    """① 風洞 axis-Mach Euler: runner_axismach.main (prepare → run_staged)。"""
    from forge_design.evaluate import runner_axismach as A
    from forge_design.evaluate import runner_wt as W
    from forge_design.probdef import load_problem
    cap = Capture()
    rd = work / "run"

    def fake_prepare(problem_path, run_dir, nsteps=None, ic_from=None, cfl_main=None, implicit_relax=None):
        # runner_axismach.prepare の n・out_int・config (`_euler_config` = `_config_euler_node` → `_apply_gas_to_config` → cfl・relax)・bcond と同じ
        p = load_problem(problem_path)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
        out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
        (run_dir / "bcondConfig.yaml").write_text(A._bcond_with_species(W._bcond(p, euler=True), A._tp_species_Y(p)))
        (run_dir / "solverConfig.yaml").write_text(A._euler_config(p, run_dir, n, out_int, cfl_main, implicit_relax))
        return {}

    with patched([(A, "prepare", fake_prepare), (A, "run_forge", cap), (A, "collect", lambda *a: {})] + _common_patches()):
        _quiet(A.main, [str(problem), str(rd)])
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def run_axismach_ns(problem, work):
    """① 風洞 axis-Mach NS: prepare_ns の config (`_ns_configs`、main_line なら PRESTAGE_CONFIG) → run_staged_ns(stages="full")。"""
    from forge_design.evaluate import runner_axismach as A
    from forge_design.evaluate import runner_wt as W
    from forge_design.probdef import load_problem
    cap = Capture()
    rd = work / "run"
    rd.mkdir(parents=True)
    p = load_problem(problem)
    # prepare_ns の該当行と同じ (n・out_int は `_ns_steps`、main_line は `_main_line_spec`、設定は `_ns_configs`、bcond、
    # 段の宣言の文脈 RECIPE_CONTEXT.json = `write_recipe_context`)
    n, out_int, _ = A._ns_steps(p, None)
    ml = A._main_line_spec(p, None, None)
    (rd / "bcondConfig.yaml").write_text(A._bcond_with_species(W._bcond(p, euler=False, profile_ok=True), A._tp_species_Y(p)))
    rec = {}
    cfg_ns, cfg_main = A._ns_configs(p, rd, n, out_int, None, None, ml, record=rec)
    A.write_recipe_context(rd, problem, {"cfl_main": None, "implicit_relax": None}, rec)
    (rd / "solverConfig.yaml").write_text(cfg_main if ml is not None else cfg_ns)
    if ml is not None:
        (rd / A.PRESTAGE_CONFIG).write_text(cfg_ns)
    (rd / "nozzle.h5").write_bytes(b"x")
    with patched([(A, "run_forge", cap)] + _common_patches()):
        _quiet(A.run_staged_ns, rd, stages="full")
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def run_walldriven(problem, work):
    """① 風洞 wall-driven: runner_walldriven.main (prepare → run_staged)。"""
    from forge_design.evaluate import runner_walldriven as D
    from forge_design.evaluate import runner_wt as W
    from forge_design.probdef import load_problem
    cap = Capture()
    rd = work / "run"

    def fake_prepare(problem_path, run_dir, nsteps=None):
        # runner_walldriven.prepare の n・out_int・`_config_euler_node`・`_bcond` と同じ
        p = load_problem(problem_path)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 16000))
        out_int = int(p.evaluate.get("outStepInterval", max(n // 4, 1)))
        (run_dir / "bcondConfig.yaml").write_text(W._bcond(p, euler=True))
        (run_dir / "solverConfig.yaml").write_text(W._config_euler_node(p, n, out_int, 4.0, 1))
        return {}

    with patched([(D, "prepare", fake_prepare), (D, "run_forge", cap)] + _common_patches()):
        _quiet(D.main, [str(problem), str(rd)])
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def _sern_fake_prepare(S, three_d):
    from forge_design.probdef import load_problem

    def fake_prepare(problem_path, run_dir, nsteps=None, op=None, wall_offset=None):
        # runner_sern(3d).prepare の作動点の選択・気体の状態・n・out_int・cfl・`_solver_config`・`_bcond_config` と同じ
        from forge_design.evaluate import runner_sern as R2
        p = load_problem(problem_path)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        R2.select_operating_point(p, op)
        st = R2.gas_states(p)
        n = int(nsteps or p.evaluate.get("nStepOuter", 4000 if three_d else 6000))
        out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
        cfl = float(p.evaluate.get("cfl_main", 1.0 if three_d else 4.0))
        cfg = S._solver_config(p, n, out_int, cfl, st["ext"]["P"])
        (run_dir / "bcondConfig.yaml").write_text(S._bcond_config(p, st))
        (run_dir / "solverConfig.yaml").write_text(cfg)
        (run_dir / "solverConfig_main.yaml").write_text(cfg)
        return {"design": {}, "moc_forces": {}, "mesh": {}}
    return fake_prepare


def _sern_collect(*a, **k):
    return {"gates": {"verdict": "PASS", "fail_class": None, "objective": None, "reasons": []}, "history": []}


def run_sern2d(problem, work):
    """⑤ SERN 2D: runner_sern.main (prepare → run_staged、段の梯子は問題の `opt:`)。"""
    from forge_design.evaluate import runner_sern as S
    cap = Capture()
    rd = work / "run"
    with patched([(S, "prepare", _sern_fake_prepare(S, False)), (S, "run_forge", cap), (S, "collect", _sern_collect),
                  (S, "restart_by_index", lambda res, mesh, **k: None)]):
        _quiet(S.main, [str(problem), str(rd)])
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


def run_sern3d(problem, work):
    """⑤ SERN 3D: runner_sern3d.main (prepare → runner_sern.run_staged)。"""
    from forge_design.evaluate import runner_sern as S
    from forge_design.evaluate import runner_sern3d as S3
    cap = Capture()
    rd = work / "run"
    with patched([(S3, "prepare", _sern_fake_prepare(S3, True)), (S, "run_forge", cap), (S3, "collect", _sern_collect),
                  (S, "restart_by_index", lambda res, mesh, **k: None)]):
        _quiet(S3.main, [str(problem), str(rd)])
    return _manifest_stages(rd, cap.cfgs), yaml.safe_load((rd / "bcondConfig.yaml").read_text())


# ------------------------------------------------------------------------------------------------ 機種の一覧

SPECS = [
    dict(name="bell_cpg_ns_sst", run=run_bell_cpg, problem="case/40.nozzle_design_tool/problem_bell_top.yaml",
         machine="③ベル CPG (node SST・断熱壁、opt/driver の CPG 評価)", dimension="2d_axisym", ref="common_node_ns_sst",
         generator="runner._solver_config + opt/driver の soft (別ディレクトリ)",
         call="opt.driver.Campaign.evaluate (prepare・forge・collect を差し替え)"),
    dict(name="bell_tp_ns_sst", run=run_bell_tp, problem="case/68.bell_two_mode/problem_bell_tp_m_ns_l100_corr.yaml",
         machine="③ベル TP (凍結組成の lump、node SST・低 Re・断熱壁)", dimension="2d_axisym", ref="common_node_ns_sst",
         generator="runner._solver_config_tp (TP_STAGES)", call="runner.run_staged_tp (prepare・forge を差し替え)"),
    dict(name="bell_tp_euler", run=run_bell_tp, problem="case/68.bell_two_mode/problem_bell_tp_m_euler_l100.yaml",
         machine="③ベル TP Euler (slip 壁、(M) 特性曲線モードの評価)", dimension="2d_axisym", ref="common_node_euler",
         generator="runner._solver_config_tp (TP_STAGES、evaluate.model: euler)", call="runner.run_staged_tp (prepare・forge を差し替え)"),
    dict(name="wt_modef_euler_cell", run=run_wt_modef, problem="case/41.wind_tunnel_design/problem_m4.yaml",
         machine="① 風洞 モード F (Euler、cell)", dimension="2d_axisym", ref="common_node_euler",
         generator="runner_wt._config_euler (cell 固定) + run_staged の文字列置換", call="runner_wt.main (prepare・forge を差し替え)"),
    dict(name="wt_feedback_euler_node", run=run_wt_feedback, problem="case/41.wind_tunnel_design/problem_m4_b10.yaml",
         machine="① 風洞 B10 の帰還 (feedback/euler_loop、Euler node、段階起動なし)", dimension="2d_axisym", ref="common_node_euler",
         generator="runner_wt._config_euler_node (euler_loop の cold のパス)", call="euler_loop の `_cfg(n, n // 4)` (n = 24000)"),
    dict(name="axismach_euler_tp", run=run_axismach_euler, problem="case/45.isobutane_m6_d155/problem_d155_euler_pin_G1_recal_mono.yaml",
         machine="① 風洞 axis-Mach Euler (TP、slip 壁)", dimension="2d_axisym", ref="common_node_euler",
         generator="runner_wt._config_euler_node → runner_axismach._apply_gas_to_config + run_staged (_cfg_set)",
         call="runner_axismach.main (prepare・forge を差し替え)"),
    dict(name="axismach_ns_sst_tp", run=run_axismach_ns, problem="case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300.yaml",
         machine="① 風洞 axis-Mach NS (TP、SST 低 Re、等温壁、main_line なし)", dimension="2d_axisym", ref="common_node_ns_sst",
         generator="runner_wt._config_sst_node → _apply_gas_to_config(viscous) + ns_stage_configs",
         call="runner_axismach._ns_configs → run_staged_ns(stages=full) (forge を差し替え)"),
    dict(name="axismach_ns_sst_tp_b0", run=run_axismach_ns,
         problem="case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml",
         machine="① 風洞 axis-Mach NS (TP、SST 低 Re、等温壁、本段に main_line の B0)", dimension="2d_axisym", ref="common_node_ns_sst",
         generator="runner_wt._config_sst_node → _apply_gas_to_config(viscous) → _apply_main_line + ns_stage_configs (PRESTAGE_CONFIG)",
         call="runner_axismach._ns_configs → run_staged_ns(stages=full) (forge を差し替え)"),
    dict(name="walldriven_euler_cpg", run=run_walldriven, problem="case/41.wind_tunnel_design/problem_m4_walldriven.yaml",
         machine="① 風洞 wall-driven (Euler、node)", dimension="2d_axisym", ref="common_node_euler",
         generator="runner_wt._config_euler_node + run_staged の文字列置換", call="runner_walldriven.main (prepare・forge を差し替え)"),
    dict(name="sern2d_finite_base_frozen_tp_sst", run=run_sern2d, problem="case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml",
         machine="⑤ SERN 2D 有限ベース (frozen TP、node SST、等温壁、既定の作動点)", dimension="2d_planar", ref="common_node_ns_sst",
         generator="runner_sern._solver_config + run_staged の正規表現・文字列置換 (問題の opt: の梯子)",
         call="runner_sern.main (prepare・forge・restart を差し替え)"),
    dict(name="sern3d_frozen_tp_sst", run=run_sern3d, problem="case/46.sern_design/problem_3d_prod_m6on_g2.yaml",
         machine="⑤ SERN 3D (frozen TP、node SST、等温壁、既定の作動点)", dimension="3d", ref="common_node_ns_sst",
         generator="runner_sern._solver_config + run_staged の正規表現・文字列置換 (問題の opt: の梯子)",
         call="runner_sern3d.main (prepare・forge・restart を差し替え)"),
]


# ------------------------------------------------------------------------------------------------ 代表の問題 YAML の固定の写し

def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixture_path(sp, fixture_dir=None):
    return Path(fixture_dir or FIXTURE_DIR) / f"{sp['name']}.yaml"


def load_manifest(fixture_dir=None):
    """MANIFEST.yaml → {機種: {source, sha256, copied}}。無い・読めない・形が違えば RuntimeError。"""
    path = Path(fixture_dir or FIXTURE_DIR) / MANIFEST_NAME
    if not path.is_file():
        raise RuntimeError(f"{path} が無い (代表の問題 YAML の写しの記録)")
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    fx = doc.get("fixtures")
    if not isinstance(fx, dict) or any(not isinstance(v, dict) or not {"source", "sha256"} <= set(v) for v in fx.values()):
        raise RuntimeError(f"{path}: `fixtures` の各項目に source・sha256 が要る")
    return fx


def verify_fixtures(specs=None, fixture_dir=None):
    """写しが MANIFEST の記録どおりか (元のパス = SPECS の problem、写しの sha256 = 記録)。違えば RuntimeError。
    元 (case の問題 YAML) との違いはここでは見ない (`fixture_status`、WARN)。"""
    man = load_manifest(fixture_dir)
    errs = []
    for sp in (specs or SPECS):
        ent, fp = man.get(sp["name"]), fixture_path(sp, fixture_dir)
        if ent is None:
            errs.append(f"{sp['name']}: MANIFEST に記録が無い")
        elif ent["source"] != sp["problem"]:
            errs.append(f"{sp['name']}: MANIFEST の元 {ent['source']} と SPECS の problem {sp['problem']} が違う")
        elif not fp.is_file():
            errs.append(f"{sp['name']}: 写し {fp} が無い")
        elif _sha256(fp) != ent["sha256"]:
            errs.append(f"{sp['name']}: 写し {fp.name} の sha256 が MANIFEST の記録と違う (写しを手で書き換えた; "
                        "元から取り直すなら --refresh-fixture)")
    if errs:
        raise RuntimeError("代表の問題 YAML の写しが記録と合わない:\n  " + "\n  ".join(errs))
    return man


def fixture_status(specs=None, fixture_dir=None, root=None):
    """写しと元 (case の問題 YAML) の比較: [(機種, 'same' | 'differs' | 'missing', 説明)]。differs・missing は WARN の扱い
    (元の編集だけでは書き写しの結果は変わらない。写しを更新するかは docstring の手順)。"""
    man = load_manifest(fixture_dir)
    out = []
    for sp in (specs or SPECS):
        ent = man.get(sp["name"]) or {}
        src = Path(root or ROOT) / sp["problem"]
        if not src.is_file():
            out.append((sp["name"], "missing", f"元 {sp['problem']} が無い (写しは {ent.get('copied', '?')} に取った)"))
        elif _sha256(src) != ent.get("sha256"):
            out.append((sp["name"], "differs", f"元 {sp['problem']} が写し (sha256 {str(ent.get('sha256'))[:12]}…、"
                                               f"{ent.get('copied', '?')}) と違う (今 {_sha256(src)[:12]}…)"))
        else:
            out.append((sp["name"], "same", sp["problem"]))
    return out


def refresh_fixtures(names, specs=None, fixture_dir=None, root=None, today=None):
    """元の問題 YAML のバイト列を写しにし、MANIFEST を書き直す。[(機種, 旧 sha256, 新 sha256)] を返す。"""
    specs = specs or SPECS
    fdir = Path(fixture_dir or FIXTURE_DIR)
    by_name = {sp["name"]: sp for sp in specs}
    if names == ["all"]:
        names = list(by_name)
    unknown = [n for n in names if n not in by_name]
    if unknown:
        raise RuntimeError(f"機種 {unknown} は SPECS に無い ({sorted(by_name)})")
    path = fdir / MANIFEST_NAME
    man = load_manifest(fdir) if path.is_file() else {}
    fdir.mkdir(parents=True, exist_ok=True)
    done = []
    for n in names:
        sp = by_name[n]
        data = (Path(root or ROOT) / sp["problem"]).read_bytes()
        fixture_path(sp, fdir).write_bytes(data)
        old = (man.get(n) or {}).get("sha256")
        man[n] = {"source": sp["problem"], "sha256": hashlib.sha256(data).hexdigest(),
                  "copied": str(today or datetime.date.today())}
        done.append((n, old, man[n]["sha256"]))
    order = [sp["name"] for sp in specs if sp["name"] in man] + sorted(set(man) - set(by_name))
    head = ("# 代表の問題 YAML の固定の写しの記録 (regen_profiles.py が読む; plan tooling-generator-defaults-and-launch-gate §5.1 #17)。\n"
            "# 手で編集しない: 写しの更新は `regen_profiles.py --refresh-fixture` (手順は regen_profiles.py の docstring)。\n"
            "# source = 元の問題 YAML (リポジトリ基準)、sha256 = 写し (= 写した時の元) のバイト列、copied = 写した日。\n")
    path.write_text(head + yaml.safe_dump({"version": 1, "fixtures": {k: man[k] for k in order}}, sort_keys=False,
                                          allow_unicode=True, width=10000), encoding="utf-8")
    return done


def _path_leaves(d, pre=()):
    if isinstance(d, dict):
        for k, v in d.items():
            yield from _path_leaves(v, pre + (str(k),))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            yield from _path_leaves(v, pre + (f"[{i}]",))
    else:
        yield pre, d


def materialize(sp, work, fixture_dir=None, root=None):
    """写しを作業ディレクトリに置き、問題 YAML の所在を基準に解決する相対パスを**元の問題 YAML の所在**の絶対パスに直す
    (opt/driver の `absolutize_problem_paths` の 4 キーと `EXTRA_PATH_KEYS`)。既知のキーの外に、元の所在に実在するファイルを指す
    相対パス (`/` か拡張子を含む文字列) があれば止まる (黙って作業ディレクトリ基準で解決しない)。置いたパスを返す。"""
    from forge_design.opt.driver import PROBLEM_PATH_KEYS, absolutize_problem_paths
    base = (Path(root or ROOT) / sp["problem"]).resolve().parent
    raw = yaml.safe_load(fixture_path(sp, fixture_dir).read_text(encoding="utf-8"))
    out = absolutize_problem_paths(raw, base)
    for path in EXTRA_PATH_KEYS:
        cur = out
        for k in path[:-1]:
            cur = cur.get(k) if isinstance(cur, dict) else None
        if isinstance(cur, dict) and isinstance(cur.get(path[-1]), str) and cur[path[-1]] \
                and not Path(cur[path[-1]]).is_absolute():
            cur[path[-1]] = str((base / cur[path[-1]]).resolve())
    known = {tuple(k) for k in PROBLEM_PATH_KEYS} | set(EXTRA_PATH_KEYS)
    stray = [(".".join(k), v) for k, v in _path_leaves(out)
             if k not in known and isinstance(v, str) and v and "\n" not in v and not Path(v).is_absolute()
             and ("/" in v or Path(v).suffix) and (base / v).exists()]
    if stray:
        raise RuntimeError(f"{sp['name']}: 写しに元の所在 ({base}) 基準と思われる相対パスがある {stray}。"
                           "解決の基準を確かめて EXTRA_PATH_KEYS に足す")
    dst = Path(work) / "problem" / Path(sp["problem"]).name
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(yaml.safe_dump(out, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return dst


# ------------------------------------------------------------------------------------------------ 書き写し

def stage_fixed(recipes, cfg_text, stage="main"):
    """段の config → fixed (problem_keys と steady_base の ignore を除く。Euler の visc・thermCond と、本段以外の段の長さ
    `time.last.nStepOuter` は残す)。"""
    doc, fails = RS.parse_config(cfg_text)
    if fails:
        raise RuntimeError(f"生成器の config を読めない: {fails}")
    flat = RS.flatten(doc)
    ignore = set(RS.resolve(recipes, "steady_base")["ignore"])
    euler = str(flat.get("physProp.viscMethod", "")) == "0"
    out = {}
    for k in sorted(flat):
        if k in ignore:
            continue
        if RS.is_problem_key(recipes, k) and not (euler and k in KEEP_FOR_EULER) \
                and not (stage != "main" and k == "time.last.nStepOuter"):
            continue
        out[k] = flat[k]
    return out


def build_profiles(recipes, specs=None, keep=False):
    """全機種を回して生成区間のプロファイル (dict) と、段ごとの config を返す。"""
    profiles, configs = {}, {}
    table = RS.explicit_table(recipes)                 # 旧い実効値の出所を今のソースで確かめる (合わなければ RecipeError)
    verify_fixtures(specs)                             # 写しが MANIFEST の記録どおり (元との違いは見ない: fixture_status)
    for sp in (specs or SPECS):
        work = Path(tempfile.mkdtemp(prefix=f"regen_{sp['name']}_"))
        problem = materialize(sp, work)                # 固定の写し (相対パスは元の問題 YAML の所在で解決)
        stages, bcond = sp["run"](problem, work)
        fixed, explicit, legs = {}, {}, {}
        for st, cfg in stages:
            fx = stage_fixed(recipes, cfg, st)
            ex = RS.explicit_for(recipes, cfg, bcond, table)
            legs.setdefault(st, []).append(cfg)
            if st in fixed and (fixed[st] != fx or explicit[st] != ex):
                raise RuntimeError(f"{sp['name']}: 段 {st} の脚の値が違う (梯子の代表の問題は脚を 1 つにする)")
            fixed[st], explicit[st] = fx, ex
        doc, _ = RS.parse_config(stages[-1][1])
        app = RS.derive_applies(doc, bcond)
        applies = {"dimension": sp["dimension"], "discretization": app["discretization"], "axisymmetric": app["axisymmetric"],
                   "physics": app["physics"], "gas": app["gas"], "wall": app.get("wall"), "mesh": "structured", "time": app["time"]}
        prof = {"inherits": "steady_base", "machine": sp["machine"],
                "applies": applies,
                "source": {"plan": PLAN, "section": "§4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6", "date": DATE},
                "transcribed_from": {"generator": sp["generator"], "problem": sp["problem"], "call": sp["call"]},
                "recommended_ref": sp["ref"],
                "stages": {st: ({"fixed": fx, "explicit": explicit[st]} if explicit[st] else {"fixed": fx})
                           for st, fx in fixed.items()}}
        prof["differs_from_recommended"] = RS.compute_differs(recipes, fixed, sp["ref"])
        profiles[sp["name"]] = prof
        configs[sp["name"]] = stages
        if not keep:
            import shutil
            shutil.rmtree(work, ignore_errors=True)
    return profiles, configs


def render_block(profiles):
    """生成区間の YAML (profiles の写像の中、2 字下げ)。"""
    body = yaml.safe_dump(profiles, sort_keys=False, allow_unicode=True, width=10000, default_flow_style=False)
    lines = [BEGIN] + ["  " + ln if ln else ln for ln in body.rstrip("\n").split("\n")] + [END]
    return "\n".join(lines) + "\n"


def current_block(text):
    i, j = text.find(BEGIN), text.find(END)
    if i < 0 or j < 0 or j < i:
        raise RuntimeError(f"{YAML_PATH} に生成区間の印が無い")
    return text[i:j + len(END)] + "\n", i, j + len(END) + 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--write", action="store_true", help="recommended_settings.yaml の生成区間を書き換える")
    g.add_argument("--check", action="store_true", help="生成区間が今の生成器の出力と一致するか (違えば終了 1)")
    g.add_argument("--check-fixtures", action="store_true",
                   help="代表の問題 YAML の写しと元 (case) の違いを WARN で出す (終了 0。写しが記録と合わなければ終了 1)")
    g.add_argument("--refresh-fixture", nargs="+", metavar="機種",
                   help="元の問題 YAML を写し直して MANIFEST を書く (`all` で全機種)。推奨値の変更の承認ではない (docstring の手順)")
    ap.add_argument("--keep", action="store_true", help="作業ディレクトリを消さない")
    a = ap.parse_args(argv)
    if a.check_fixtures:
        try:
            verify_fixtures()
        except RuntimeError as e:
            print(f"FAIL {e}")
            return 1
        st = fixture_status()
        for n, s, msg in st:
            print(f"{'ok  ' if s == 'same' else 'WARN'} {n}: {msg}")
        nw = sum(s != "same" for _, s, _ in st)
        print(f"写し {len(st)} 本、元と違う・元が無い {nw} 本" + (" (写しの更新は docstring の手順。推奨値の変更の承認とは別)" if nw else ""))
        return 0
    if a.refresh_fixture:
        for n, old, new in refresh_fixtures(a.refresh_fixture):
            print(f"写した: {n} ({(old or '記録なし')[:12]} → {new[:12]})")
        print("次に --check: 不一致なら推奨値の変更なので、差を plan に書いて判断を得てから --write (docstring の手順)")
        return 0
    text = YAML_PATH.read_text(encoding="utf-8")
    recipes = RS.load_recipes(text=text)
    profiles, _ = build_profiles(recipes, keep=a.keep)
    block = render_block(profiles)
    cur, i, j = current_block(text)
    if a.check:
        ok = cur == block
        print("一致" if ok else "不一致: 生成器の出力か代表の問題 YAML の写しが変わった。推奨値の変更なので、差を plan に書いて"
                                 "判断を得てから `regen_profiles.py --write` で書き直す")
        return 0 if ok else 1
    if a.write:
        new = text[:i] + block + text[j:]
        RS.load_recipes(text=new)                      # 書く前に読めることを確かめる
        YAML_PATH.write_text(new, encoding="utf-8")
        print(f"書いた: {YAML_PATH} ({len(profiles)} 機種)")
        return 0
    sys.stdout.write(block)
    return 0


if __name__ == "__main__":
    sys.exit(main())
