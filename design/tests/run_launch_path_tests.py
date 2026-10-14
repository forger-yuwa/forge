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
    差し替え (呼ばれた回数だけ見る)。

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
    """forge の代わり: solverConfig.yaml を控え、nStepOuter の res を書く。`soft_edit` は最初の起動の res の場を壊す関数。"""

    def __init__(self, sst=False, soft_edit=None, n_res=N):
        self.cfgs, self.sst, self.soft_edit, self.n_res = [], sst, soft_edit, n_res

    def __call__(self, run_dir):
        run_dir = Path(run_dir)
        cfg = (run_dir / "solverConfig.yaml").read_text()
        self.cfgs.append(cfg)
        import re
        n = int(re.search(r"nStepOuter: (\d+)", cfg).group(1))
        v = fields(self.n_res, step=n, sst=self.sst, wd_offset=5.0)
        if self.soft_edit is not None and len(self.cfgs) == 1:
            self.soft_edit(v)
        write_h5(run_dir / f"res_{n}.h5", v, self.n_res)
        (run_dir / "residual_history.csv").write_text(f"step,rms_ro\n0,1.0\n{n},{1e-3 * len(self.cfgs)}\n")
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


def run_sern(tag, call, pre=None):
    rd = TMP / tag
    rd.mkdir()
    (rd / "solverConfig.yaml").write_text(SERN_CFG)
    (rd / "solverConfig_main.yaml").write_text(SERN_CFG)
    (rd / "bcondConfig.yaml").write_text(BCOND)
    if pre is not None:
        pre(rd)
    fake = FakeForge(sst=True)
    restarts = []
    orig = (RS.run_forge, RS.restart_by_index, RS.warm_from_run)
    RS.run_forge = fake
    RS.restart_by_index = lambda res, mesh: restarts.append(Path(res).name)
    RS.warm_from_run = lambda dst, src: {"src": str(src)}
    try:
        rc = call(rd)
    finally:
        RS.run_forge, RS.restart_by_index, RS.warm_from_run = orig
    return rd, fake, restarts, rc


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
    _, f, _, _ = run_sern("g_sern", lambda rd: RS.run_staged(rd, "full", **SERN_LADDER))
    _, f2, _, _ = run_sern("g_sern_r", lambda rd: RS.run_staged(rd, "full", **SERN_RAMP))
    out["sern_ladder"] = golden_of(f.cfgs + ["--"] + f2.cfgs)
    _, f, _, _ = run_sern("g_sern_w", lambda rd: RS.run_staged(rd, "full", warm_src=TMP, soft_cfl=0.1, mid_steps=4000,
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
rd, f, rs, rc = run_sern("sern", lambda rd: RS.run_staged(rd, "full", **SERN_LADDER))
rd2, f2, rs2, rc2 = run_sern("sern_r", lambda rd: RS.run_staged(rd, "full", **SERN_RAMP))
check("(d) SERN: 段の config は変更前とバイト同一 (固定の梯子と CFL の梯子)",
      golden_of(f.cfgs + ["--"] + f2.cfgs) == GOLDEN["sern_ladder"], golden_of(f.cfgs + ["--"] + f2.cfgs)[:16])
check("(d) SERN: manifest の段 = 層流暖機 → soft → mid → 本段",
      manifest_tags(rd) == ["warm_cfl0.1", "soft_cfl0.1", "mid_cfl0.1", "main"], manifest_tags(rd))
check("(d) SERN: 段の履歴 residual_<tag>.csv が段ごとに残る",
      all((rd / f"residual_{t}.csv").exists() for t in manifest_tags(rd)))
check("(d) SERN: 段間の restart_by_index は前の段の数だけ", len(rs) == 3, rs)
tags2 = manifest_tags(rd2)
check("(d) SERN: 同じ CFL を繰り返す梯子でも tag が重ならない",
      tags2 == ["warm_cfl0.05", "warm_cfl0.05_2", "soft_cfl0.2", "soft_cfl0.5", "mid_cfl0.2", "mid_cfl0.5", "main"]
      and len(set(json.loads((rd2 / "stage_manifest.json").read_text())["stages"][i]["history"] for i in range(len(tags2))))
      == len(tags2), tags2)
man = json.loads((rd / "stage_manifest.json").read_text())
check("(d) SERN: 本段の段の鍵は層流暖機と違う (check_convergence --segment が暖機を本段に連結しない)",
      man["stages"][0]["key"] != man["stages"][-1]["key"])
rd3, f3, rs3, _ = run_sern("sern_w", lambda rd: RS.run_staged(rd, "full", warm_src=TMP, soft_cfl=0.1, mid_steps=4000,
                                                              warm_adapt_steps=500))
check("(d) SERN warm start: 段の config は変更前とバイト同一", golden_of(f3.cfgs) == GOLDEN["sern_warm"], golden_of(f3.cfgs)[:16])
check("(d) SERN warm start: manifest = 適応段 → mid → 本段", manifest_tags(rd3) == ["warm_adapt", "mid_cfl0.1", "main"],
      manifest_tags(rd3))
rd4, f4, rs4, _ = run_sern("sern_main", lambda rd: RS.run_staged(rd, "main"),
                           pre=lambda rd: (write_h5(rd / "res_4000.h5", fields(sst=True)),
                                           (rd / "stage_manifest.json").write_text(json.dumps(
                                               {"manifest_version": 2, "stages": [{"order": 0, "tag": "soft_cfl0.1",
                                                                                    "history": "residual_soft_cfl0.1.csv",
                                                                                    "restart_from": None, "key": {}, "soft": {}}]}))))
check("(d) SERN stages=main: 既存の記録に本段を足す、本段の config = solverConfig_main.yaml",
      manifest_tags(rd4) == ["soft_cfl0.1", "main"] and f4.cfgs == [SERN_CFG] and rs4 == ["res_4000.h5"], manifest_tags(rd4))
rd5, f5, _, _ = run_sern("sern_none", lambda rd: RS.run_staged(rd, "none"))
check("(d) SERN stages=none: solverConfig.yaml をそのまま回し、本段を記録", manifest_tags(rd5) == ["main"] and f5.cfgs == [SERN_CFG])

if FAIL == 0:
    import shutil
    shutil.rmtree(TMP, ignore_errors=True)
print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}  (一時ディレクトリ {TMP} を残した)'}")
sys.exit(1 if FAIL else 0)
