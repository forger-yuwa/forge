#!/usr/bin/env python3
"""起動の経路の段の記録と同一メッシュの段間の引き継ぎの試験 (plan tooling-generator-defaults-and-launch-gate §4.5・§5.1 #9)。

forge の代わりに res_<n>.h5 (実際の HDF5) と residual_history.csv を書く偽の forge を差し込み、各 runner の段階起動を回す:

(a) runner_wt.run_staged (Euler cell・NS node SST の 2 つの本段の config): 段の config の文字列が変更前 (HEAD 570790e8) と
    バイト同一 (sha256 を固定)、`stage_manifest.json` の段が soft → main、段の履歴 residual_history_<tag>.csv が残る、
    段間は `restart_field.py` (本物を子プロセスで呼ぶ): 初期場の保存量が soft の最終場とビット一致、`wall_dist` は宛先のまま
    (res 側の wall_dist は写さない)、`restart_field.log` が残る。
(b) restart の負例 (runner_wt・walldriven): soft の最終場に (1) 保存量 roe が無い、(2) roUx に Inf、(3) ρ < 0 が 1 点、
    (4) NS の config で roOmega が無い → RuntimeError (段終了ゲート) で本段を起動しない (forge の起動は 1 回だけ)、
    初期場 (nozzle.h5) は書き換わらない。(5) 節点数の違う res → restart_field.py が拒否して本段を起動しない。
(c) runner_walldriven.run_staged: (a) と同じ (段の config は変更前とバイト同一、manifest soft → main、restart_field)。
(d) runner_sern.run_staged (層流暖機・soft・mid の梯子、warm start の適応段、stages="main"・"none"): 段の config の文字列が
    変更前とバイト同一、`stage_manifest.json` の段の並び (warm_cfl* → soft_cfl* → mid_cfl* → main)、段の履歴 residual_<tag>.csv、
    同じ CFL を繰り返す梯子でも tag が重ならない、stages="main" は既存の記録に本段を足す。段間の引き継ぎ (restart_by_index) は
    **差し替えない** (2026-10-15、plan §5.1 #19: 段終了ゲート + 本物の restart_field.py。初期場が前の段の最終場とビット一致、
    restart_field.log の VERDICT の数 = 段間の数)。warm start の作動点間のリマップ (warm_from_run) だけは差し替える。
(g) SERN の restart の負例 (§5.1 #19): soft の最終場に #9 の 5 種 → 本段を起動しない (節点数違いは restart_field.py の拒否 = ValueError)、
    初期場は書き換わらない。層流暖機 → SST: 層流の res に roK・roOmega が無くても前の段のゲートは通り、足す量は宛先の初期場のまま
    (restart_field.log に記録)、足す量の ω ≤ 0・k の非有限・宛先に roOmega が無い → soft を起動しない。
(h) euler_loop の継続チャンク (`start_cont_chunk`): 同じ格子の写しへ段終了ゲート + restart_field.py (ビット一致、wall_dist は宛先、
    2 つ目のチャンクは 1 つ目の res から)。負例 (Euler なので roOmega 以外の 4 種 + 節点数違い) でチャンクを起動しない。
(i) CPG のベルの同一形状の継続 (`runner.continue_cpg`): 継続元の変換済みの格子を写して index コピー (作り直さない)、config は
    `_solver_config` と同じ。負例 5 種で書き込まず RuntimeError。opt/driver の CPG の評価 (soft → 本段) と opt/polish の CPG の
    チャンク継続を偽の forge で通し、ゲート落ちは本段を起動しない (driver は RETRYABLE)。

`--golden` は段の config の sha256 を出すだけ (変更前のコードで回して固定値を作るため; forge・restart は差し替え)。

#11 (2026-10-15、plan §5.1 #11) から、機種のプロファイルに合う段の config (ここでは runner_wt の Euler cell と walldriven) には生成器が
明示すべきキーを旧い実効値で足し、段の宣言 RECIPE.json・RECIPE_<tag>.json を書く。バイト同一の比較は、宣言の `explicit_added` の
断片を除いた文字列で行う (足したキーの値と例外表の一致は design/tests/run_generator_recipe_tests.py)。合うプロファイルの無い config
(NS v1 の node SST・CPG の SERN) は今までどおりバイト同一で、宣言も書かない。
"""
import argparse
import hashlib
import json
import os
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--golden", action="store_true", help="段の config の sha256 を出すだけ (変更前のコードで固定値を作る)")
ap.add_argument("--design", default=str(Path(__file__).resolve().parents[1]), help="forge_design のある design/ (既定: この試験の隣)")
args = ap.parse_args()
sys.path.insert(0, args.design)
sys.path.insert(0, str(Path(args.design).resolve().parent / "solver_density_cuda" / "tools" / "recipes"))
sys.path.insert(0, str(Path(args.design).resolve().parent / "solver_density_cuda" / "tools"))

from forge_design.evaluate import runner as RR  # noqa: E402
from forge_design.evaluate import runner_sern as RS  # noqa: E402
from forge_design.evaluate import runner_walldriven as RD  # noqa: E402
from forge_design.evaluate import runner_wt as RW  # noqa: E402

FAIL = 0
TMP = Path(tempfile.mkdtemp(prefix="launch_path_"))
N = 40


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def h5sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


P = types.SimpleNamespace(cp=1004.5, gamma=1.4)
CFG_EULER_CELL = RW._config_euler(P, 12000, 4000, 4.0, 1)        # runner_wt.prepare(euler=True) の本段
CFG_SST_NODE = RW._config_sst_node(P, 12000, 4000, 4.0)          # runner_wt.prepare(euler=False) の本段
CFG_EULER_NODE = RW._config_euler_node(P, 16000, 4000, 4.0, 1)   # runner_walldriven.prepare の本段
BCOND = "wall: {physID: 3, kind: slip, outputHDFflg: 1, ints: , floats: }\n"

# 変更前のコード (HEAD 570790e8) で回した段の config の sha256 (`--golden` で作った値)
GOLDEN = {
    "wt_euler_cell": "4c6acb6a5bd3da4f4560fc339ca9742e87067b75fe963cc3d5c7be0930010525",
    "wt_sst_node": "844c91c6ca260a1631edec6587be68f9807a0691a5e9f524d09199727ca6f0ab",
    "walldriven": "1217829500be9793ee20cca42b905d0e7dcd4bd7ae37d44b16d75b219531fa52",
    "sern_ladder": "67ecb0e57a105f90f21813ccd240859c68a495fddf51b156947a2084fb5cbde2",
    "sern_warm": "b075fb3549d679e7027372a9df17588305ca17944bd48da6b9ded6730c5c0cd3",
}


def fields(n=N, step=0, sst=False, wd_offset=0.0):
    """保存量の場 (一様に近い有限・正の値)。wall_dist は宛先と res で違う値にする (写さないことを見る)。"""
    x = np.linspace(0.0, 1.0, n)
    v = {"ro": 1.0 + 0.1 * x + 1e-4 * step, "roUx": 100.0 + x, "roUy": 0.5 * x, "roUz": np.zeros(n),
         "roe": 2.5e5 + 10.0 * x + step, "wall_dist": 1e-3 + x + wd_offset}
    if sst:
        v["roK"] = 1.0 + x
        v["roOmega"] = 1.8e4 + x
    return {k: np.asarray(a, dtype=np.float64) for k, a in v.items()}


def write_h5(path, v, n=N, mesh=False):
    """res (forge の出力と同じく /MESH/COORD・/MESH/CONNE・/VALUE) を書く。mesh=True は変換済みの格子の代わり (境界の所属も持つ;
    restart_field.py は res をその run の格子に座標・接続で結び付ける: plan §5.1 #20)。"""
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.column_stack([np.linspace(0, 1, n), np.zeros(n), np.zeros(n)]).ravel())
        f.create_dataset("MESH/CONNE", data=np.array([x for i in range(n) for x in (1, i)], dtype=np.int32))
        if mesh:
            f.create_dataset("BCONDS/3/iCells", data=np.arange(n, dtype=np.int32))
            f["BCONDS/3"].attrs["bcondKind"] = "slip"
        for k, a in v.items():
            f.create_dataset(f"VALUE/{k}", data=a)


def read_value(path):
    with h5py.File(path, "r") as f:
        return {k: np.asarray(f["VALUE"][k]) for k in f["VALUE"]}


class FakeForge:
    """forge の代わり: solverConfig.yaml を控え、nStepOuter の res を書く。`soft_edit` は最初の起動の res の場を壊す関数。
    `cfg_edit(cfg, v)` は config を見て res の場を変える関数 (層流暖機の res から roK・roOmega を落とすなど)。起動ごとに
    格子 (`mesh.meshFileName`) の初期場を `inits` に控え (段の初期場を見るため)、収束判定の印 CONVERGENCE_VERDICT.txt も書く。"""

    def __init__(self, sst=False, soft_edit=None, n_res=N, cfg_edit=None):
        self.cfgs, self.sst, self.soft_edit, self.n_res, self.cfg_edit = [], sst, soft_edit, n_res, cfg_edit
        self.inits = []

    def __call__(self, run_dir):
        run_dir = Path(run_dir)
        cfg = (run_dir / "solverConfig.yaml").read_text()
        self.cfgs.append(cfg)
        import re
        m = re.search(r'meshFileName: "([^"]+)"', cfg)
        mesh = run_dir / (m.group(1) if m else "nozzle.h5")
        self.inits.append(read_value(mesh) if mesh.exists() else None)
        n = int(re.search(r"nStepOuter: (\d+)", cfg).group(1))
        v = fields(self.n_res, step=n, sst=self.sst, wd_offset=5.0)
        if self.cfg_edit is not None:
            self.cfg_edit(cfg, v)
        if self.soft_edit is not None and len(self.cfgs) == 1:
            self.soft_edit(v)
        write_h5(run_dir / f"res_{n}.h5", v, self.n_res)
        (run_dir / "residual_history.csv").write_text(f"step,rms_ro\n0,1.0\n{n},{1e-3 * len(self.cfgs)}\n")
        (run_dir / "CONVERGENCE_VERDICT.txt").write_text("VERDICT: PASS (fake forge)\n")
        return 0


class FakeSub:
    """変更前のコードの interp_field.py の呼び出しを差し替える (--golden 用)。"""
    @staticmethod
    def run(*a, **k):
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")


def wt_dir(tag, cfg, sst=False):
    rd = TMP / tag
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(BCOND)
    write_h5(rd / "nozzle.h5", fields(sst=sst), mesh=True)
    return rd


def run_wt(tag, cfg, sst=False, soft_edit=None, n_res=N, mod=RW):
    rd = wt_dir(tag, cfg, sst)
    fake = FakeForge(sst=sst, soft_edit=soft_edit, n_res=n_res)
    orig = mod.run_forge
    mod.run_forge = fake
    orig_sub = getattr(mod, "subprocess", None)
    if args.golden:
        mod.subprocess = FakeSub
    err = None
    try:
        if mod is RW:
            RW.run_staged(None, rd)
        else:
            RD.run_staged(rd)
    except RuntimeError as e:
        err = e
    finally:
        mod.run_forge = orig
        if args.golden:
            mod.subprocess = orig_sub
    return rd, fake, err


def manifest_tags(rd):
    p = Path(rd) / "stage_manifest.json"
    return [s["tag"] for s in json.loads(p.read_text())["stages"]] if p.exists() else None


def golden_of(cfgs):
    return sha("\n=====\n".join(cfgs))


def strip_explicit(cfg, added):
    """#11 で生成器が足した明示すべきキー (宣言の explicit_added) の断片を除く (flow 形式の `, key: 値`)。断片が 1 か所に決まらなければ
    AssertionError。"""
    from forge_design.evaluate.config_builder import _token
    for k, v in (added or {}).items():
        frag = f", {k.split('.')[-1]}: {_token(v)}"
        assert cfg.count(frag) == 1, (k, frag)
        cfg = cfg.replace(frag, "", 1)
    return cfg


def stripped(rd, cfgs, tags):
    """段の config から、その段の宣言 (RECIPE_<tag>.json) の explicit_added を除いた文字列の列 (宣言の無い段はそのまま)。"""
    out = []
    for cfg, tag in zip(cfgs, tags):
        p = Path(rd) / f"RECIPE_{tag}.json"
        out.append(strip_explicit(cfg, json.loads(p.read_text()).get("explicit_added")) if p.exists() else cfg)
    return out


# ------------------------------------------------------------------------------------------------ SERN (段の梯子)
SERN_CFG = '''mesh: {discretization: "node", isAxisymmetric: 0, nodeWallDirichlet: 1, meshFileName: "sern.h5", valueFileName: "sern.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 12000}
  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5,
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 2000
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 1, limiter: 2, pRef: 2851.0, limiterScaled: 1, venkatK: 0.05}
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0}
initial: "uniform_p101325_u10"
'''


def n_restarts(rd):
    """restart_field.py の成功の数 (restart_field.log の VERDICT: OK の行)。"""
    log = Path(rd) / "restart_field.log"
    return log.read_text().count("VERDICT: OK") if log.exists() else 0


def run_sern(tag, call, pre=None, fake=None, sst_mesh=True):
    """SERN の段階起動を偽の forge で回す。段間の restart_by_index は**本物** (段終了ゲート + restart_field.py; §5.1 #19)。
    格子と初期場は sern.h5 (prepare が書くものの代わり; sst_mesh=False は roK・roOmega の無い初期場)。戻り
    (run, 偽の forge, restart の成功の数, rc, 例外)。"""
    rd = TMP / tag
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(SERN_CFG)
    (rd / "solverConfig_main.yaml").write_text(SERN_CFG)
    (rd / "bcondConfig.yaml").write_text(BCOND)
    write_h5(rd / "sern.h5", fields(sst=sst_mesh), mesh=True)
    if pre is not None:
        pre(rd)
    fake = fake or FakeForge(sst=True)
    orig = (RS.run_forge, RS.warm_from_run)
    RS.run_forge = fake
    RS.warm_from_run = lambda dst, src: {"src": str(src)}
    rc, err = None, None
    try:
        rc = call(rd)
    except (RuntimeError, ValueError) as e:
        err = e
    finally:
        RS.run_forge, RS.warm_from_run = orig
    return rd, fake, n_restarts(rd), rc, err


SERN_LADDER = dict(soft_steps=4000, soft_cfl=0.1, warm_lam_steps=4000, warm_lam_cfl=0.1, mid_steps=4000)
SERN_RAMP = dict(soft_steps=3000, soft_cfl=0.5, warm_lam_steps=2000, warm_lam_ramp=[0.05, 0.05], mid_steps=2000,
                 soft_ramp=[0.2, 0.5])

if args.golden:
    out = {}
    _, f, _ = run_wt("g_wt_e", CFG_EULER_CELL)
    out["wt_euler_cell"] = golden_of(f.cfgs)
    _, f, _ = run_wt("g_wt_s", CFG_SST_NODE, sst=True)
    out["wt_sst_node"] = golden_of(f.cfgs)
    _, f, _ = run_wt("g_wd", CFG_EULER_NODE, mod=RD)
    out["walldriven"] = golden_of(f.cfgs)
    _, f, _, _, _ = run_sern("g_sern", lambda rd: RS.run_staged(rd, "full", **SERN_LADDER))
    _, f2, _, _, _ = run_sern("g_sern_r", lambda rd: RS.run_staged(rd, "full", **SERN_RAMP))
    out["sern_ladder"] = golden_of(f.cfgs + ["--"] + f2.cfgs)
    _, f, _, _, _ = run_sern("g_sern_w", lambda rd: RS.run_staged(rd, "full", warm_src=TMP, soft_cfl=0.1, mid_steps=4000,
                                                                   warm_adapt_steps=500))
    out["sern_warm"] = golden_of(f.cfgs)
    print(json.dumps(out, indent=1))
    sys.exit(0)

# --- (a) runner_wt -------------------------------------------------------------------------------------------------
for tag, cfg, sst in (("wt_euler_cell", CFG_EULER_CELL, False), ("wt_sst_node", CFG_SST_NODE, True)):
    rd, f, err = run_wt(tag, cfg, sst=sst)
    check(f"(a) {tag}: 例外なしで soft → 本段 (forge 2 回)", err is None and len(f.cfgs) == 2, str(err)[:200])
    base = stripped(rd, f.cfgs, ["soft", "main"])
    check(f"(a) {tag}: 段の config は変更前とバイト同一 (#11 の明示すべきキーの断片を除いて)", golden_of(base) == GOLDEN[tag],
          golden_of(base)[:16])
    check(f"(a) {tag}: 本段の config = prepare の config (#11 の明示すべきキーを除いて書き換えない)", base[-1] == cfg)
    decl = (rd / "RECIPE_main.json").exists() and (rd / "RECIPE_soft.json").exists()
    check(f"(a) {tag}: 段の宣言 (#11): Euler cell (wt_modef_euler_cell に合う) は soft・main を宣言、NS v1 は宣言しない (config もそのまま)",
          decl if not sst else (not list(rd.glob("RECIPE*")) and f.cfgs[-1] == cfg))
    check(f"(a) {tag}: stage_manifest の段は soft → main", manifest_tags(rd) == ["soft", "main"], manifest_tags(rd))
    check(f"(a) {tag}: 段の履歴 residual_history_soft/main.csv",
          (rd / "residual_history_soft.csv").exists() and (rd / "residual_history_main.csv").exists())
    v = read_value(rd / "nozzle.h5")
    want = fields(step=3000, sst=sst, wd_offset=5.0)
    moved = [k for k in want if k != "wall_dist"]
    check(f"(a) {tag}: 初期場の保存量 = soft の最終場 (ビット一致、restart_field.py)",
          all(np.array_equal(v[k], want[k]) for k in moved), str(moved))
    check(f"(a) {tag}: wall_dist は宛先のまま (res の wall_dist は写さない)",
          np.array_equal(v["wall_dist"], fields(sst=sst)["wall_dist"]))
    log = (rd / "restart_field.log")
    check(f"(a) {tag}: restart_field.log が残る (interp_field.py でない)",
          log.exists() and "VERDICT: OK" in log.read_text())

# --- (b) restart の負例 -----------------------------------------------------------------------------------------------
NEG = [("roe が無い", lambda v: v.pop("roe"), False, "roe: 無い"),
       ("roUx に Inf", lambda v: v["roUx"].__setitem__(3, np.inf), False, "roUx: 非有限"),
       ("ρ < 0 が 1 点", lambda v: v["ro"].__setitem__(7, -1.0), False, "ρ ≤ 0"),
       ("NS で roOmega が無い", lambda v: v.pop("roOmega"), True, "roOmega: 無い")]
for mod, name in ((RW, "runner_wt"), (RD, "walldriven")):
    for i, (what, edit, sst, must) in enumerate(NEG):
        if mod is RD and sst:
            continue                    # walldriven は Euler だけ
        cfg = (CFG_SST_NODE if sst else (CFG_EULER_CELL if mod is RW else CFG_EULER_NODE))
        rd, f, err = run_wt(f"neg_{name}_{i}", cfg, sst=sst, soft_edit=edit, mod=mod)
        init = h5sha(TMP / f"neg_{name}_{i}" / "nozzle.h5")
        ok = (err is not None and "段終了ゲート" in str(err) and must in str(err) and len(f.cfgs) == 1)
        check(f"(b) {name}: {what} → 本段を起動しない", ok, str(err)[:160])
        rd0 = wt_dir(f"neg_ref_{name}_{i}", cfg, sst)
        check(f"(b) {name}: {what} → 初期場は書き換わらない", init == h5sha(rd0 / "nozzle.h5"))
        check(f"(b) {name}: {what} → manifest は soft だけ", manifest_tags(rd) == ["soft"], manifest_tags(rd))
    rd, f, err = run_wt(f"neg_{name}_n", CFG_EULER_CELL if mod is RW else CFG_EULER_NODE, n_res=N + 1, mod=mod)
    check(f"(b) {name}: 節点数の違う res → restart_field.py が拒否して本段を起動しない",
          err is not None and "restart_field.py が失敗" in str(err) and len(f.cfgs) == 1, str(err)[:160])

# --- (c) walldriven ----------------------------------------------------------------------------------------------------
rd, f, err = run_wt("walldriven", CFG_EULER_NODE, mod=RD)
check("(c) walldriven: 例外なしで soft → 本段", err is None and len(f.cfgs) == 2, str(err)[:200])
check("(c) walldriven: 段の config は変更前とバイト同一 (#11 の明示すべきキーの断片を除いて)",
      golden_of(stripped(rd, f.cfgs, ["soft", "main"])) == GOLDEN["walldriven"], golden_of(stripped(rd, f.cfgs, ["soft", "main"]))[:16])
check("(c) walldriven: 段の宣言 RECIPE_soft・RECIPE_main (#11)", (rd / "RECIPE_soft.json").exists() and (rd / "RECIPE_main.json").exists())
check("(c) walldriven: stage_manifest の段は soft → main", manifest_tags(rd) == ["soft", "main"], manifest_tags(rd))
v = read_value(rd / "nozzle.h5")
check("(c) walldriven: 初期場 = soft の最終場 (restart_field.py)、wall_dist は宛先のまま",
      np.array_equal(v["ro"], fields(step=3000)["ro"]) and np.array_equal(v["wall_dist"], fields()["wall_dist"])
      and "VERDICT: OK" in (rd / "restart_field.log").read_text())

# --- (d) SERN ----------------------------------------------------------------------------------------------------------
rd, f, rs, rc, err = run_sern("sern", lambda rd: RS.run_staged(rd, "full", **SERN_LADDER))
rd2, f2, rs2, rc2, err2 = run_sern("sern_r", lambda rd: RS.run_staged(rd, "full", **SERN_RAMP))
check("(d) SERN: 例外なしで全段を起動 (段間は本物の restart_by_index)", err is None and err2 is None and len(f.cfgs) == 4,
      f"{err} / {err2}"[:200])
check("(d) SERN: 段の config は変更前とバイト同一 (固定の梯子と CFL の梯子)",
      golden_of(f.cfgs + ["--"] + f2.cfgs) == GOLDEN["sern_ladder"], golden_of(f.cfgs + ["--"] + f2.cfgs)[:16])
check("(d) SERN: manifest の段 = 層流暖機 → soft → mid → 本段",
      manifest_tags(rd) == ["warm_cfl0.1", "soft_cfl0.1", "mid_cfl0.1", "main"], manifest_tags(rd))
check("(d) SERN: 段の履歴 residual_<tag>.csv が段ごとに残る",
      all((rd / f"residual_{t}.csv").exists() for t in manifest_tags(rd)))
check("(d) SERN: 段間の restart_by_index (restart_field.py の VERDICT: OK) は前の段の数だけ", rs == 3 and rs2 == 6, (rs, rs2))
check("(d) SERN: 各段の初期場の保存量 = 前の段の最終場 (ビット一致、restart_field.py)、wall_dist は宛先のまま",
      all(all(np.array_equal(f.inits[i + 1][k], fields(step=4000, sst=True, wd_offset=5.0)[k])
              for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"))
          and np.array_equal(f.inits[i + 1]["wall_dist"], fields()["wall_dist"]) for i in range(3)))
check("(d) SERN: 層流暖機 → SST で足す roK・roOmega は層流の res の値を写した (restart_field.log に記録)",
      "roK (前の段の res の値を写した), roOmega (前の段の res の値を写した)" in (rd / "restart_field.log").read_text())
tags2 = manifest_tags(rd2)
check("(d) SERN: 同じ CFL を繰り返す梯子でも tag が重ならない",
      tags2 == ["warm_cfl0.05", "warm_cfl0.05_2", "soft_cfl0.2", "soft_cfl0.5", "mid_cfl0.2", "mid_cfl0.5", "main"]
      and len(set(json.loads((rd2 / "stage_manifest.json").read_text())["stages"][i]["history"] for i in range(len(tags2))))
      == len(tags2), tags2)
man = json.loads((rd / "stage_manifest.json").read_text())
check("(d) SERN: 本段の段の鍵は層流暖機と違う (check_convergence --segment が暖機を本段に連結しない)",
      man["stages"][0]["key"] != man["stages"][-1]["key"])
rd3, f3, rs3, _, err3 = run_sern("sern_w", lambda rd: RS.run_staged(rd, "full", warm_src=TMP, soft_cfl=0.1, mid_steps=4000,
                                                                    warm_adapt_steps=500))
check("(d) SERN warm start: 段の config は変更前とバイト同一", golden_of(f3.cfgs) == GOLDEN["sern_warm"], golden_of(f3.cfgs)[:16])
check("(d) SERN warm start: manifest = 適応段 → mid → 本段", manifest_tags(rd3) == ["warm_adapt", "mid_cfl0.1", "main"],
      manifest_tags(rd3))
rd4, f4, rs4, _, err4 = run_sern("sern_main", lambda rd: RS.run_staged(rd, "main"),
                           pre=lambda rd: (write_h5(rd / "res_4000.h5", fields(sst=True)),
                                           (rd / "stage_manifest.json").write_text(json.dumps(
                                               {"manifest_version": 2, "stages": [{"order": 0, "tag": "soft_cfl0.1",
                                                                                    "history": "residual_soft_cfl0.1.csv",
                                                                                    "restart_from": None, "key": {}, "soft": {}}]}))))
check("(d) SERN stages=main: 既存の記録に本段を足す、本段の config = solverConfig_main.yaml",
      manifest_tags(rd4) == ["soft_cfl0.1", "main"] and f4.cfgs == [SERN_CFG] and rs4 == 1 and err4 is None, manifest_tags(rd4))
rd5, f5, _, _, _ = run_sern("sern_none", lambda rd: RS.run_staged(rd, "none"))
check("(d) SERN stages=none: solverConfig.yaml をそのまま回し、本段を記録", manifest_tags(rd5) == ["main"] and f5.cfgs == [SERN_CFG])

# --- (g) SERN の restart の負例 (§5.1 #19) -------------------------------------------------------------------------------
SERN_SOFT_MAIN = dict(soft_steps=4000, soft_cfl=0.1)
for i, (what, edit, _sst, must) in enumerate(NEG):
    rd, f, rs, _, err = run_sern(f"neg_sern_{i}", lambda rd: RS.run_staged(rd, "full", **SERN_SOFT_MAIN),
                                 fake=FakeForge(sst=True, soft_edit=edit))
    ref = TMP / f"neg_sern_ref_{i}"
    ref.mkdir()
    write_h5(ref / "sern.h5", fields(sst=True), mesh=True)
    check(f"(g) SERN: {what} → 本段を起動しない (段終了ゲート)",
          isinstance(err, RuntimeError) and "段終了ゲート" in str(err) and must in str(err) and len(f.cfgs) == 1 and rs == 0,
          str(err)[:160])
    check(f"(g) SERN: {what} → 初期場 sern.h5 は書き換わらない", h5sha(rd / "sern.h5") == h5sha(ref / "sern.h5"))
rd, f, rs, _, err = run_sern("neg_sern_n", lambda rd: RS.run_staged(rd, "full", **SERN_SOFT_MAIN),
                             fake=FakeForge(sst=True, n_res=N + 1))
check("(g) SERN: 節点数の違う res → restart_field.py が拒否 (ValueError) して本段を起動しない",
      isinstance(err, ValueError) and "restart_field.py が失敗" in str(err) and len(f.cfgs) == 1, str(err)[:160])


def lam_drop(cfg, v):
    """層流暖機の段 (turbulence: none) の res から roK・roOmega を落とす (前の段のゲートが足す量を求めないことを見る)。"""
    if 'model: "none"' in cfg:
        v.pop("roK", None)
        v.pop("roOmega", None)


def lam_edit(fn):
    return lambda cfg, v: fn(v) if 'model: "none"' in cfg else None


SERN_LAM = dict(soft_steps=4000, soft_cfl=0.1, warm_lam_steps=4000, warm_lam_cfl=0.1)
rd, f, rs, _, err = run_sern("lam_drop", lambda rd: RS.run_staged(rd, "full", **SERN_LAM), fake=FakeForge(sst=True, cfg_edit=lam_drop))
check("(g) SERN 層流 → SST: 層流の res に roK・roOmega が無くても前の段 (層流) のゲートは通り、soft・本段を起動",
      err is None and len(f.cfgs) == 3 and manifest_tags(rd) == ["warm_cfl0.1", "soft_cfl0.1", "main"], str(err)[:200])
check("(g) SERN 層流 → SST: soft の初期場の roK・roOmega = 宛先の初期場のまま、他の保存量 = 層流の最終場",
      err is None and np.array_equal(f.inits[1]["roOmega"], fields(sst=True)["roOmega"])
      and np.array_equal(f.inits[1]["roK"], fields(sst=True)["roK"])
      and np.array_equal(f.inits[1]["ro"], fields(step=4000, wd_offset=5.0)["ro"]))
check("(g) SERN 層流 → SST: 足す量の出所を restart_field.log に記録",
      "roK (宛先の初期場のまま), roOmega (宛先の初期場のまま)" in (rd / "restart_field.log").read_text())
LAM_NEG = [("足す roOmega に 0 が 1 点", lam_edit(lambda v: v["roOmega"].__setitem__(4, 0.0)), True,
            "roOmega (次の段で足す量、SRC の値): ≤ 0 が 1 点"),
           ("足す roK に Inf", lam_edit(lambda v: v["roK"].__setitem__(2, np.inf)), True, "roK (次の段で足す量、SRC の値): 非有限"),
           ("層流の res にも宛先にも roOmega が無い", lam_drop, False, "roOmega: 宛先の初期場 sern.h5 に無い"),
           ("層流の res に roe が無い (前の段のゲート)", lam_edit(lambda v: v.pop("roe")), True, "roe: 無い")]
for i, (what, edit, sst_mesh, must) in enumerate(LAM_NEG):
    rd, f, rs, _, err = run_sern(f"lam_neg_{i}", lambda rd: RS.run_staged(rd, "full", **SERN_LAM),
                                 fake=FakeForge(sst=True, cfg_edit=edit), sst_mesh=sst_mesh)
    check(f"(g) SERN 層流 → SST: {what} → soft を起動しない",
          isinstance(err, RuntimeError) and must in str(err) and len(f.cfgs) == 1 and rs == 0, str(err)[:200])

# --- (h) euler_loop の継続チャンク (§5.1 #19) ----------------------------------------------------------------------------
from forge_design.feedback import euler_loop as EL  # noqa: E402


def el_pass(tag, edit=None, n_res=N):
    """euler_loop のパスの dir (格子 nozzle.h5・config・bcond・probe) と、そのチャンクの最終場 res_6000.h5。"""
    rd = TMP / tag
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(CFG_EULER_NODE)
    (rd / "bcondConfig.yaml").write_text(BCOND)
    (rd / "probe.yaml").write_text("probes: []\n")
    write_h5(rd / "nozzle.h5", fields(), mesh=True)
    v = fields(n_res, step=6000, wd_offset=5.0)
    if edit is not None:
        edit(v)
    write_h5(rd / "res_6000.h5", v, n_res)
    return rd


rd = el_pass("el")
cd1 = EL.start_cont_chunk(rd, 1, rd / "res_6000.h5", CFG_EULER_NODE)
v = read_value(cd1 / "nozzle.h5")
want = fields(step=6000, wd_offset=5.0)
check("(h) euler_loop: 継続チャンクの初期場 = 直前のチャンクの最終場 (ビット一致、restart_field.py)、wall_dist は宛先のまま",
      all(np.array_equal(v[k], want[k]) for k in want if k != "wall_dist") and np.array_equal(v["wall_dist"], fields()["wall_dist"])
      and "VERDICT: OK" in (cd1 / "restart_field.log").read_text())
check("(h) euler_loop: チャンクの config・bcond・probe を置く", (cd1 / "solverConfig.yaml").read_text() == CFG_EULER_NODE
      and (cd1 / "bcondConfig.yaml").read_text() == BCOND and (cd1 / "probe.yaml").exists())
write_h5(cd1 / "res_6000.h5", fields(step=12000, wd_offset=5.0))
cd2 = EL.start_cont_chunk(rd, 2, cd1 / "res_6000.h5", CFG_EULER_NODE)
check("(h) euler_loop: 2 つ目のチャンクは 1 つ目の res から (格子は同じパスの写し)",
      np.array_equal(read_value(cd2 / "nozzle.h5")["roe"], fields(step=12000)["roe"]))
for i, (what, edit, sst, must) in enumerate(NEG):
    if sst:
        continue                    # euler_loop は Euler だけ
    rd = el_pass(f"el_neg_{i}", edit=edit)
    launched = []
    try:
        cd = EL.start_cont_chunk(rd, 1, rd / "res_6000.h5", CFG_EULER_NODE)
        launched.append(cd)          # run_loop はここで宣言を書いて起動する
        err = None
    except RuntimeError as e:
        err = e
    check(f"(h) euler_loop: {what} → チャンクを起動しない (段終了ゲート)",
          err is not None and "段終了ゲート" in str(err) and must in str(err) and not launched, str(err)[:160])
    check(f"(h) euler_loop: {what} → チャンクの初期場は格子の写しのまま",
          h5sha(rd / "cont_1" / "nozzle.h5") == h5sha(rd / "nozzle.h5"))
rd = el_pass("el_neg_n", n_res=N + 1)
try:
    EL.start_cont_chunk(rd, 1, rd / "res_6000.h5", CFG_EULER_NODE)
    err = None
except RuntimeError as e:
    err = e
check("(h) euler_loop: 節点数の違う res → restart_field.py が拒否してチャンクを起動しない",
      err is not None and "restart_field.py が失敗" in str(err), str(err)[:160])

# --- (i) CPG のベルの同一形状の継続 (§5.1 #19) ----------------------------------------------------------------------------
import regen_profiles as RG  # noqa: E402
import yaml  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

BELL = RG.materialize(next(sp for sp in RG.SPECS if sp["name"] == "bell_cpg_ns_sst"), TMP / "bell_problem")
PB = load_problem(BELL)
CFG_BELL = RR._solver_config(PB, 12000, 4000)
BC_BELL = RR._bcond_config(PB)


def bell_src(tag, edit=None, n_res=N, cfg=CFG_BELL, step=12000):
    """CPG の継続元の run (変換済みの格子 nozzle.h5・config・bcond・probe・prepare_info・品質の記録) と最終場 res_<step>.h5。"""
    rd = TMP / tag
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(BC_BELL)
    (rd / "probe.yaml").write_text("probes: []\n")
    (rd / "MESH_QUALITY.txt").write_text("VERDICT: PASS (fake)\n")
    (rd / "prepare_info.json").write_text(json.dumps({"problem": str(BELL), "run_dir": str(rd), "nsteps": step}))
    write_h5(rd / "nozzle.h5", fields(sst=True), mesh=True)
    v = fields(n_res, step=step, sst=True, wd_offset=5.0)
    if edit is not None:
        edit(v)
    write_h5(rd / f"res_{step}.h5", v, n_res)
    return rd


src = bell_src("bell_src")
info = RR.continue_cpg(BELL, src, TMP / "bell_p1")
v = read_value(TMP / "bell_p1" / "nozzle.h5")
want = fields(step=12000, sst=True, wd_offset=5.0)
check("(i) continue_cpg: 継続元の格子を写し、初期場 = 継続元の最終場 (ビット一致)、wall_dist は宛先のまま",
      all(np.array_equal(v[k], want[k]) for k in want if k != "wall_dist") and np.array_equal(v["wall_dist"], fields()["wall_dist"])
      and "VERDICT: OK" in (TMP / "bell_p1" / "restart_field.log").read_text())
check("(i) continue_cpg: config = runner._solver_config (prepare の CPG と同じ)、bcond・品質の記録を写す、res は写さない",
      (TMP / "bell_p1" / "solverConfig.yaml").read_text() == CFG_BELL and (TMP / "bell_p1" / "bcondConfig.yaml").read_text() == BC_BELL
      and (TMP / "bell_p1" / "MESH_QUALITY.txt").exists() and not list((TMP / "bell_p1").glob("res_*")))
check("(i) continue_cpg: prepare_info に継続元の res", info.get("continued_from", "").endswith("bell_src/res_12000.h5")
      and info["nsteps"] == 12000)
for i, (what, edit, _sst, must) in enumerate(NEG):
    s_ = bell_src(f"bell_neg_{i}", edit=edit)
    try:
        RR.continue_cpg(BELL, s_, TMP / f"bell_neg_{i}_p1")
        err = None
    except RuntimeError as e:
        err = e
    check(f"(i) continue_cpg: {what} → RuntimeError (段終了ゲート)", err is not None and "段終了ゲート" in str(err) and must in str(err),
          str(err)[:160])
    check(f"(i) continue_cpg: {what} → 写した格子の初期場は書き換わらない",
          h5sha(TMP / f"bell_neg_{i}_p1" / "nozzle.h5") == h5sha(s_ / "nozzle.h5"))
s_ = bell_src("bell_neg_n", n_res=N + 1)
try:
    RR.continue_cpg(BELL, s_, TMP / "bell_neg_n_p1")
    err = None
except RuntimeError as e:
    err = e
check("(i) continue_cpg: 節点数の違う res → restart_field.py が拒否", err is not None and "restart_field.py が失敗" in str(err),
      str(err)[:160])

# opt/driver の CPG の評価 (soft → 本段) を偽の forge で
from forge_design.opt import driver as DRV  # noqa: E402
from forge_design.opt import polish as POL  # noqa: E402


def drv_eval(tag, soft_edit=None):
    work = TMP / tag
    work.mkdir()
    fake = FakeForge(sst=True, soft_edit=soft_edit)

    def fake_prepare(prob, run_dir, nsteps=None, warm_from=None, stage=None):
        # runner.prepare の CPG (格子の生成・変換・IC の代わりに合成の格子と初期場)
        p = load_problem(prob)
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=False)
        n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
        (run_dir / "solverConfig.yaml").write_text(RR._solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
        (run_dir / "bcondConfig.yaml").write_text(RR._bcond_config(p))
        (run_dir / "probe.yaml").write_text("probes: []\n")
        write_h5(run_dir / "nozzle.h5", fields(sst=True), mesh=True)
        (run_dir / "prepare_info.json").write_text(json.dumps({"problem": str(prob), "run_dir": str(run_dir), "nsteps": n}))
        return {}

    camp = DRV.Campaign(BELL, work / "campaign", seed_res=str(work / "seed_res.h5"))
    raw = yaml.safe_load(BELL.read_text())
    x = [float(raw["dv"][k]["value"]) for k in camp.dv_order]
    fake_collect = lambda prob, rd: {"eta_cf": 0.97, "mdot_kg_s": camp.mdot_1d, "thrust_N": 1.0,   # noqa: E731
                                     "res_file": sorted(Path(rd).glob("res_[0-9]*.h5"))[-1].name}
    saved = (RR.prepare, RR.run_forge, RR.collect, DRV.Campaign._diverged, DRV.Campaign._quasisteady, DRV.Campaign._warm_for)
    RR.prepare, RR.run_forge, RR.collect = fake_prepare, fake, fake_collect
    DRV.Campaign._diverged = lambda self, rd, out: False
    DRV.Campaign._quasisteady = lambda self, rd: "ALL STEADY"
    DRV.Campaign._warm_for = lambda self, x: work / "seed_res.h5"
    try:
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            row = camp.evaluate(x, "e1")
    finally:
        (RR.prepare, RR.run_forge, RR.collect, DRV.Campaign._diverged, DRV.Campaign._quasisteady, DRV.Campaign._warm_for) = saved
    return work / "campaign", fake, row


cdir, f, row = drv_eval("drv")
v = f.inits[1] if len(f.inits) > 1 else {}
check("(i) driver CPG: soft → 本段 (forge 2 回)、PASS", len(f.cfgs) == 2 and row["status"] == "PASS", json.dumps(row, ensure_ascii=False)[:200])
check("(i) driver CPG: 本段の初期場 = soft の最終場 (ビット一致、soft の格子の写し)、wall_dist は宛先のまま",
      bool(v) and all(np.array_equal(v[k], fields(step=3000, sst=True, wd_offset=5.0)[k]) for k in ("ro", "roUx", "roe", "roK", "roOmega"))
      and np.array_equal(v["wall_dist"], fields()["wall_dist"]) and "VERDICT: OK" in (cdir / "e1" / "restart_field.log").read_text())
check("(i) driver CPG: 本段の config (明示すべきキーを除く) = runner._solver_config",
      strip_explicit(f.cfgs[1], json.loads((cdir / "e1" / "RECIPE_main.json").read_text()).get("explicit_added")) == CFG_BELL)
cdir, f, row = drv_eval("drv_neg", soft_edit=NEG[1][1])
check("(i) driver CPG: soft の最終場に Inf → 本段を起動しない、RETRYABLE (soft 段の失敗)",
      len(f.cfgs) == 1 and row["status"] == "FAIL" and row["fail_class"] == "RETRYABLE" and "段終了ゲート" in row["note"],
      json.dumps(row, ensure_ascii=False)[:240])

# opt/polish の CPG のチャンク継続
src = bell_src("pol_src")
fake = FakeForge(sst=True)
etas = iter([0.95, 0.9502])
saved = (RR.run_forge, RR.collect)
RR.run_forge, RR.collect = fake, (lambda prob, rd: {"eta_cf": next(etas)})
try:
    out = POL.polish_run(BELL, src, 0.90, tol=3e-4, max_chunks=2)
    err = None
except RuntimeError as e:
    out, err = None, e
finally:
    RR.run_forge, RR.collect = saved
check("(i) polish CPG: _p1 は継続元から、_p2 は _p1 から (forge 2 回、頭打ちで settled)",
      err is None and len(fake.cfgs) == 2 and out["settled"] and out["chunks"] == 2, str(err or out)[:200])
_ci = [json.loads((TMP / f"pol_src_p{k}" / "prepare_info.json").read_text()).get("continued_from", "") for k in (1, 2)] \
    if err is None else ["", ""]
check("(i) polish CPG: 各チャンクの初期場 = 前のチャンクの最終場 (ビット一致; _p1 ← pol_src、_p2 ← _p1 の res)",
      err is None and all(np.array_equal(fake.inits[i]["roe"], fields(step=12000, sst=True)["roe"]) for i in range(2))
      and _ci[0].endswith("pol_src/res_12000.h5") and _ci[1].endswith("pol_src_p1/res_12000.h5")
      and "VERDICT: OK" in (TMP / "pol_src_p2" / "restart_field.log").read_text(), _ci)
src = bell_src("pol_neg", edit=NEG[2][1])
fake = FakeForge(sst=True)
RR.run_forge, RR.collect = fake, (lambda prob, rd: {"eta_cf": 0.95})
try:
    POL.polish_run(BELL, src, 0.90, max_chunks=2)
    err = None
except RuntimeError as e:
    err = e
finally:
    RR.run_forge, RR.collect = saved
check("(i) polish CPG: 継続元の最終場に ρ < 0 → チャンクを起動しない", err is not None and "段終了ゲート" in str(err) and not fake.cfgs,
      str(err)[:160])

if FAIL == 0:
    import shutil
    shutil.rmtree(TMP, ignore_errors=True)
print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}  (一時ディレクトリ {TMP} を残した)'}")
sys.exit(1 if FAIL else 0)
