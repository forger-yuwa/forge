#!/usr/bin/env python3
"""case/60 (TN D-7275 試験 26・28、8-ft HTST 大型校正パネル) の forge run を生成して段階起動する。

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
TW = 300.0
# 試験ごとの Table II / III (digitize_*.json と同じ値)。alpha > 0 は斜め衝撃波の背後の局所状態を平板に平行な入口にする
TESTS = {
    26: dict(M=6.64, p=2117.0, Tt=1867.0, alpha=0.0, Taw=1728.0, rvcs=0.413 * 68.89e3),
    28: dict(M=6.60, p=2144.0, Tt=1817.0, alpha=9.8, Taw=1689.0, rvcs=1.46 * 69.91e3),
}


def oblique(g, T1, U1, p1, theta_deg):
    """熱的完全気体 (h(T)) の斜め衝撃波。弱い解の背後の (T, U, p, ρ, β)。"""
    from scipy.optimize import brentq
    r1 = p1 / (g.R * T1)

    def post(beta):
        un1, ut = U1 * np.sin(beta), U1 * np.cos(beta)

        def f(eps):
            un2 = un1 * eps; p2 = p1 + r1 * un1 ** 2 * (1 - eps)
            h2 = g.h(T1) + 0.5 * (un1 ** 2 - un2 ** 2)
            T2 = brentq(lambda T: g.h(T) - h2, 100, 3000)
            return p2 - r1 / eps * g.R * T2
        eps = brentq(f, 0.05, 0.99)
        un2 = un1 * eps; p2 = p1 + r1 * un1 ** 2 * (1 - eps)
        T2 = brentq(lambda T: g.h(T) - (g.h(T1) + 0.5 * (un1 ** 2 - un2 ** 2)), 100, 3000)
        return beta - np.arctan2(un2, ut), T2, float(np.hypot(un2, ut)), p2, r1 / eps
    M1 = U1 / (g.gamma(T1) * g.R * T1) ** 0.5
    b = brentq(lambda b: np.degrees(post(b)[0]) - theta_deg, np.arcsin(1 / M1) + 1e-3, np.radians(40))
    _, T2, U2, p2, r2 = post(b)
    return dict(T=T2, U=U2, p=p2, rho=r2, beta_deg=float(np.degrees(b)))


def freestream(test=26):
    c = TESTS[test]
    g = htst(c["Tt"])
    h0 = g.h(c["Tt"]); lo, hi = 100.0, 600.0
    for _ in range(80):
        T = 0.5 * (lo + hi); U = c["M"] * (g.gamma(T) * g.R * T) ** 0.5
        if g.h(T) + 0.5 * U * U > h0:
            hi = T
        else:
            lo = T
    free = dict(T=T, U=U, rho=c["p"] / (g.R * T), p=c["p"])
    loc = oblique(g, T, U, c["p"], c["alpha"]) if c["alpha"] > 0 else dict(free)
    k = 1.5 * (0.005 * loc["U"]) ** 2          # Tu 0.5 %、μt/μ 10 (入口 = 平板の縁の状態)
    om = loc["rho"] * k / (10.0 * g.mu(loc["T"]))
    return g, dict(loc, k=k, om=om, free=free)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--test", type=int, default=26, choices=sorted(TESTS))
    ap.add_argument("--mesh", default="fp_d7275_y3")
    ap.add_argument("--main-steps", type=int, default=40000)
    ap.add_argument("--out-int", type=int, default=5000)
    ap.add_argument("--cfl", type=float, default=1.5)
    ap.add_argument("--forge-tools", default=None)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    gas, fs = freestream(a.test)
    tc = TESTS[a.test]
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
        d7275_test=a.test, M=tc["M"], alpha_deg=tc["alpha"], p_inf=fs["free"]["p"], T_inf=fs["free"]["T"], U_inf=fs["free"]["U"],
        ro_inf=fs["free"]["rho"], inlet_local=dict(T=fs["T"], U=fs["U"], p=fs["p"], rho=fs["rho"], beta_deg=fs.get("beta_deg")),
        Tt=tc["Tt"], Tw=TW, phi=gas.phi, Y=gas.Y, R=gas.R, Taw_table=tc["Taw"], rhoVcp_star_l=tc["rvcs"], mesh=a.mesh,
        cuda_blocksize=128, gen_args=sys.argv[1:]), indent=2, ensure_ascii=False))
    mc56.m50.IC_DELTA0 = mc56.IC_DELTA0
    mc56.m50.patch_ic(rd / "mesh.h5", dict(rho=fs["rho"], U=fs["U"], T=fs["T"], p=fs["p"]), gas, TW)
    mc56.seed_turb(rd / "mesh.h5", fs["k"], fs["om"])
    print(f"{a.run}: 試験 {a.test}  入口 (平板の縁) T {fs['T']:.2f} K  U {fs['U']:.1f}  ρ {fs['rho']:.5f}  p {fs['p']:.0f}  Tw {TW}"
          f"  (自由流 T∞ {fs['free']['T']:.2f} K、β {fs.get('beta_deg')})")
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
