#!/usr/bin/env python3
"""run_0509–0511 (2026-09-27): va3 M4.19 L_c8 を入口 lump モル分率 (bcond X0/X1) で段階起動する。

    python3 run_lumpX_staged.py run_0509_... [run_0510_... ...]

run dir は `runner_axismach --prepare-only --cfl 6 --implicit-relax 0.7` で作ったもの。
- 入口 bcond の Y0/Y1 (質量分率) を lump (MIXDRY/H2O) のモル分率 X0/X1 に置き換える
  (値は prepare_info.json の正規化済みモル分率。forge が species_db の MW で Y に換算する)。
- 段階起動 (recommended-settings §1.2): S0_soft (1 次, cfl 0.5, nStepInner 10, 3000) →
  S1_mid (1 次, cfl 1.0, nStepInner 10, 3000) → S2_main (2 次, cfl 6 + implicitRelax 0.7, nStepInner 4)。
- 段間は同一メッシュなので restart_field.py (保存量の index コピー)。runner の run_staged は
  interp_field.py (cross-mesh 用) を使うため使わない (AGENTS.md「メッシュ変更後の restart」)。
- 段ごとの実効設定を stage_manifest.json に残す (check_convergence.py --segment 用)。
"""
import json, re, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))
from stage_manifest import StageManifest  # noqa: E402


def lump_X(run):
    sp = json.loads((run / "prepare_info.json").read_text())["species"]
    X, keep, tr = sp["X"], sp["keep"], sp["transported"]
    out = []
    for s in tr:
        out.append(X[s] if s in keep else sum(v for k, v in X.items() if k not in keep))
    tot = sum(out)
    return [x / tot for x in out]


def set_inlet_X(run):
    b = (run / "bcondConfig.yaml").read_text()
    Xs = lump_X(run)
    ins = ", ".join(f"X{i}: {x:.12g}" for i, x in enumerate(Xs))
    b2 = re.sub(r"Y0: [\d.eE+-]+, Y1: [\d.eE+-]+", ins, b, count=1)
    if b2 == b or "Y0" in b2:
        raise RuntimeError("inlet Y0/Y1 の置換に失敗")
    (run / "bcondConfig.yaml").write_text(b2)
    return Xs


def main_cfg(run):
    c = (run / "solverConfig.yaml").read_text()
    c = c.replace("nStepInner: 5", "nStepInner: 4")          # recommended-settings §1 (2026-09-12)
    if "condensation:" in c and "output:" not in c:
        c = c.rstrip("\n") + "\noutput: {level: 2}\n"          # §3: condLim_/condClampCorr_ の確認用
    return c


def stage_cfg(c, conv, cfl, ninner, nsteps):
    c = re.sub(r"convMethod: \d", f"convMethod: {conv}", c)
    c = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {cfl}, cfl_pseudo: {cfl}", c)
    c = re.sub(r"nStepInner: \d+", f"nStepInner: {ninner}", c)
    c = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", c)
    return c


def last_res(run):
    res = sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    return res[-1] if res else None


def run_one(run):
    run = Path(run).resolve()
    Xs = set_inlet_X(run)
    print(f"[{run.name}] inlet X = {Xs}", flush=True)
    cmain = main_cfg(run)
    n_main = int(re.search(r"nStepOuter: (\d+)", cmain).group(1))
    bc = (run / "bcondConfig.yaml").read_text()
    stages = [("S0_soft", stage_cfg(cmain, 0, 0.5, 10, 3000).replace("outStepInterval: 4000", "outStepInterval: 3000"), 3000),
              ("S1_mid", stage_cfg(cmain, 0, 1.0, 10, 3000).replace("outStepInterval: 4000", "outStepInterval: 3000"), 3000),
              ("S2_main", cmain, n_main)]
    sm = StageManifest(run)
    for i, (tag, cfg, n) in enumerate(stages):
        (run / "solverConfig.yaml").write_text(cfg)
        chk = subprocess.run([sys.executable, str(TOOLS / "check_solver_config.py"), str(run)],
                             capture_output=True, text=True)
        (run / f"CONFIG_CHECK_{tag}.txt").write_text(chk.stdout + chk.stderr)
        if chk.returncode != 0:
            raise RuntimeError(f"{tag}: check_solver_config FAIL\n{chk.stdout}")
        rc = subprocess.run([str(TOOLS / "run_case.sh"), str(run)], capture_output=True, text=True)
        (run / f"run_case_stdout_{tag}.log").write_text(rc.stdout + rc.stderr)
        shutil.copy(run / "residual_history.csv", run / f"residual_history_{tag}.csv")
        shutil.copy(run / "forge_run.log", run / f"forge_run_{tag}.log")
        sm.add(tag, cfg, bc, history=f"residual_history_{tag}.csv")
        sm.write()
        r = last_res(run)
        step = int(r.stem.split("_")[1]) if r else -1
        print(f"[{run.name}] {tag}: rc={rc.returncode} last res step={step}", flush=True)
        if rc.returncode != 0 or step < n:
            raise RuntimeError(f"{tag} 段が失敗 (rc {rc.returncode}, last step {step})")
        if i < len(stages) - 1:
            subprocess.run([sys.executable, str(TOOLS / "restart_field.py"), str(r), str(run / "nozzle.h5")],
                           check=True, capture_output=True, text=True)
            for f in list(run.glob("res_*")):
                f.unlink()


if __name__ == "__main__":
    for r in sys.argv[1:]:
        run_one(r)
