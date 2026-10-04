#!/usr/bin/env python3
"""case/60 (TN D-7275 試験 26、8-ft HTST 大型校正パネル) の forge run を生成して段階起動する。

    python3 gen_runs.py --run run_0001_t26 [--mesh fp_d7275_y3] [--forge-tools DIR] [--dry]

- 自由流: Table II の M∞ 6.64・p∞ 2117 Pa・Tt 1867 K。T∞ は燃焼ガスモデル (case/50 の CombustionProducts、
  φ は Tt から) で全エンタルピーを合わせて解く (Table III の T_l 228.3 K に対し 227.7 K、Re∞ −1 %、(ρVcp)∞ −1.4 %)。
- 壁: 等温 300 K (原報『ほぼ室温で投入』。Tw は表に無い。比較は St で行うので一次の影響は打ち消される)。
- 物性・熱力学・乱流の設定は case/56 の 2D 平板 (`tools/make_case.py` の CFG) と同じ。`qAccumulatorFP64: 1` を足す。
- 段: lam → soft → mid → 2 次ランプ 0.3/0.6/1.0 → main。段の引き継ぎは **restart_field.py** (同一メッシュ)、
  stage_manifest.json を書く。全域 FP64 で回すときは --forge-tools に倍精度ビルドの tools を渡す。
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import h5py, numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
C56 = ROOT / "case" / "56.gap_tp1187" / "tools"
sys.path.insert(0, str(TOOLS)); sys.path.insert(0, str(C56))
from stage_manifest import StageManifest                        # noqa: E402
import make_case as mc56                                        # noqa: E402  (CFG / BC / m50 / seed_turb)
from gas_htst import htst                                       # noqa: E402

ENV = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128",
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
M_INF, P_INF, TT = 6.64, 2117.0, 1867.0
TW = 300.0


def freestream():
    g = htst(TT)
    h0 = g.h(TT); lo, hi = 100.0, 600.0
    for _ in range(80):
        T = 0.5 * (lo + hi); U = M_INF * (g.gamma(T) * g.R * T) ** 0.5
        if g.h(T) + 0.5 * U * U > h0:
            hi = T
        else:
            lo = T
    ro = P_INF / (g.R * T)
    k = 1.5 * (0.005 * U) ** 2
    om = ro * k / (10.0 * g.mu(T))
    return g, dict(T=T, U=U, rho=ro, p=P_INF, k=k, om=om)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", default="fp_d7275_y3")
    ap.add_argument("--main-steps", type=int, default=40000)
    ap.add_argument("--out-int", type=int, default=5000)
    ap.add_argument("--cfl", type=float, default=1.5)
    ap.add_argument("--forge-tools", default=None)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    gas, fs = freestream()
    rd = HERE / a.run
    if rd.exists():
        raise SystemExit(f"{rd} は既にある。消さない")
    rd.mkdir()
    shutil.copy(HERE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    mc56.m50.write_species_db(rd)
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    ymf = ", ".join(f"Y{i}: {gas.Y[sp]:.8f}" for i, sp in enumerate(mc56.SPECIES))
    bc = mc56.BC.format(ro=fs["rho"], U=fs["U"], p=fs["p"], t=fs["T"], tw=TW, kinf=fs["k"], ominf=fs["om"],
                        ymf=ymf, ints="")
    (rd / "bcondConfig.yaml").write_text(bc)
    a_inf = float(np.sqrt(gas.gamma(fs["T"]) * gas.R * fs["T"]))
    common = dict(mu_inf=gas.mu(900.0), lam_inf=gas.lam(900.0), cp_ref=gas.cp(fs["T"]), gam_ref=gas.gamma(fs["T"]),
                  relax=1.0, kinf=fs["k"], ominf=fs["om"], turb="sst", ro_ref=fs["rho"], p_ref=fs["p"], a_ref=a_inf)
    stages = [("lam", dict(conv=0, lim=0, cfl=0.2, inner=10, nsteps=3000, out_int=3000, **{**common, "turb": "none"})),
              ("soft", dict(conv=0, lim=0, cfl=0.2, inner=10, nsteps=3000, out_int=3000, **common)),
              ("mid", dict(conv=0, lim=0, cfl=0.5, inner=10, nsteps=3000, out_int=3000, **common))]
    for i, c in enumerate((0.3, 0.6, 1.0)):
        stages.append((f"ramp{i}_cfl{c:g}", dict(conv=1, lim=2, cfl=c, inner=6, nsteps=2500, out_int=2500, **common)))
    stages.append(("main", dict(conv=1, lim=2, cfl=a.cfl, inner=4, nsteps=a.main_steps, out_int=a.out_int, **common)))

    def cfg(kw):
        return mc56.CFG.format(species=", ".join(mc56.SPECIES), **kw).replace("detectNaN: 1", "detectNaN: 1, qAccumulatorFP64: 1")

    sm = StageManifest(rd)
    for tag, st in stages:
        sm.add(tag, cfg(st), bc, history="residual_history.csv" if tag == "main" else f"residual_history_{tag}.csv")
    sm.write()
    (rd / "solverConfig.yaml").write_text(cfg(stages[0][1]))
    (rd / "case_setup.json").write_text(json.dumps(dict(
        d7275_test=26, M=M_INF, p_inf=fs["p"], T_inf=fs["T"], U_inf=fs["U"], ro_inf=fs["rho"], Tt=TT, Tw=TW,
        phi=gas.phi, Y=gas.Y, R=gas.R, Taw_table=1728.0, rhoVcp_star_l=0.413 * 68.89e3, mesh=a.mesh,
        cuda_blocksize=128, gen_args=sys.argv[1:]), indent=2, ensure_ascii=False))
    mc56.m50.IC_DELTA0 = mc56.IC_DELTA0
    mc56.m50.patch_ic(rd / "mesh.h5", dict(rho=fs["rho"], U=fs["U"], T=fs["T"], p=fs["p"]), gas, TW)
    mc56.seed_turb(rd / "mesh.h5", fs["k"], fs["om"])
    print(f"{a.run}: T∞ {fs['T']:.2f} K  U∞ {fs['U']:.1f}  ρ∞ {fs['rho']:.5f}  p∞ {fs['p']:.0f}  Tw {TW}")
    if a.dry:
        return
    run_tools = Path(a.forge_tools) if a.forge_tools else TOOLS
    for tag, st in stages:
        (rd / "solverConfig.yaml").write_text(cfg(st))
        print(f"--- {tag}: {st['nsteps']} step, cfl {st['cfl']}", flush=True)
        rc = subprocess.run([str(run_tools / "run_case.sh"), str(rd)], env=ENV).returncode
        cur = rd / f"res_{st['nsteps']}.h5"
        if rc != 0 or not cur.exists() or list(rd.glob("res_nan_*.h5")):
            raise SystemExit(f"stage {tag}: rc={rc} res={cur.exists()} → 中断")
        if tag == "main":
            break
        for suf in ("forge_run.log", "residual_history.csv", "CONVERGENCE_VERDICT.txt", "residual_history.png"):
            if (rd / suf).exists():
                (rd / suf).rename(rd / f"{Path(suf).stem}_{tag}{Path(suf).suffix}")
        subprocess.run([sys.executable, str(TOOLS / "restart_field.py"), str(cur), str(rd / "mesh.h5")], env=ENV, check=True)
        mc56.seed_turb(rd / "mesh.h5", fs["k"], fs["om"])
        for p in list(rd.glob("res_*")):
            p.unlink()                                  # 段の出力は残さない (AWS のディスク逼迫)
    print("main 終了")


if __name__ == "__main__":
    main()
