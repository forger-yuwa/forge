#!/usr/bin/env python3
"""case/56 の **3D** すきま run を組む (縦すきま × 横すきま交差)。

2D 平板の断面を `inletProfile` で入口に与え、初期場も同じ断面で埋める
(一様 IC だと入口分布と食い違って序盤で発散する)。段構成は 2D の `make_case.py` と
同じ梯子 (層流暖機 → SST soft → mid → 2 次ランプ → 本段)。

    python3 tools/make_case3d.py --run run_0049_ab3d_wall --mesh gap3d_c50 \\
        --inlet-run run_0002_fp_t8_long --inlet-x 1.7076 [--dry]

`--dry` は投入せず run ディレクトリだけ作る。段の実行は `run_case.sh` 経由。
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(HERE))
import importlib.util as _ilu                                      # noqa: E402
_sp = _ilu.spec_from_file_location("c56_make_case", HERE / "make_case.py")
mc = _ilu.module_from_spec(_sp); _sp.loader.exec_module(mc)
from gas_htst import htst                                          # noqa: E402
import extract_inlet_table as eit                                  # noqa: E402

SPECIES = mc.SPECIES
BC_SYM = """sym:    {physID: 7, kind: slip, outputHDFflg: 0, ints: , floats: }
"""


def patch_ic_profile(h5, tab, gas, Tw):
    """入口断面 (y の関数) を全域に敷く。すきま内 (y<0) は静止・壁温。

    `case/50` の `patch_ic` は一様自由流 + tanh で、$\\delta_{99}$ 2.6 cm の
    厚い乱流 BL には合わない (入口分布と食い違って序盤で発散する)。
    """
    y = tab["y"]
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        yy = np.clip(c[:, 1], y[0], y[-1])
        ro = np.interp(yy, y, tab["ro"])
        u = np.interp(yy, y, tab["Ux"])
        ps = np.interp(yy, y, tab["Ps"])
        kk = np.interp(yy, y, tab["k"])
        om = np.interp(yy, y, tab["omega"])
        ingap = c[:, 1] < 0.0                    # すきま内は静止・壁温から始める
        u[ingap] = 0.0
        T = ps / (gas.R * ro)
        T[ingap] = Tw
        ro[ingap] = ps[ingap] / (gas.R * Tw)
        wall = wd <= 0.0
        u[wall] = 0.0
        T[wall] = Tw
        ro[wall] = ps[wall] / (gas.R * Tw)
        uT = np.unique(T)
        emap = {t: gas.h(t) - gas.R * t for t in uT}
        e = np.array([emap[t] for t in T])
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros(n, np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = (ro * (e + 0.5 * u ** 2)).astype(np.float32)
        for name, v in (("roK", ro * kk), ("roOmega", ro * om)):
            key = "/VALUE/" + name
            if key in f:
                f[key][:] = v.astype(np.float32)
            else:
                f.create_dataset(key, data=v.astype(np.float32))
        for i, s in enumerate(SPECIES):
            key = f"/VALUE/roY{i}"
            v = (ro * gas.Y[s]).astype(np.float32)
            if key in f:
                f[key][:] = v
            else:
                f.create_dataset(key, data=v)
        print(f"  IC: 節点 {n}, すきま内 {int(ingap.sum())}, 壁 {int(wall.sum())}, "
              f"U 端 {u.max():.1f} m/s, T 範囲 {T.min():.1f}-{T.max():.1f} K")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", required=True, help="mesh/<name>.h5")
    ap.add_argument("--series", type=int, default=8)
    ap.add_argument("--tw", type=float, default=300.0)
    ap.add_argument("--inlet-run", default="run_0002_fp_t8_long")
    ap.add_argument("--inlet-x", type=float, default=1.7076)
    ap.add_argument("--soft-steps", type=int, default=2000)
    ap.add_argument("--mid-steps", type=int, default=3000)
    ap.add_argument("--ramp", default="0.3,0.6,1.0")
    ap.add_argument("--ramp-steps", type=int, default=2500)
    ap.add_argument("--soft-cfl", type=float, default=0.2)
    ap.add_argument("--mid-cfl", type=float, default=0.5)
    ap.add_argument("--cfl", type=float, default=1.5)
    ap.add_argument("--main-steps", type=int, default=60000)
    ap.add_argument("--out-int", type=int, default=10000)
    ap.add_argument("--blocksize", type=int, default=128,
                    help="FORGE_CUDA_BLOCKSIZE。既定 512 は node SLAU のレジスタ上限 "
                         "(REG 136 -> 481 スレッド) を超えて起動できない。**値は結果に効く**ので "
                         "A/B の両腕で必ず揃える (2D 生産 run も 128)")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    os.environ["FORGE_CUDA_BLOCKSIZE"] = str(a.blocksize)

    der = {o["run"]: o for o in json.loads(
        (CASE / "derived.json").read_text(encoding="utf-8"))["series"]}
    o = der[a.series]
    gas = htst(o["Tt_c"])
    T, U, rho, p = o["T_inf"], o["U_inf"], o["ro_inf"], o["p_inf"]
    mu = gas.mu(T)
    k_inf = 1.5 * (0.005 * U) ** 2
    om_inf = rho * k_inf / (10.0 * mu)

    rd = CASE / a.run
    rd.mkdir(exist_ok=True)
    for f in list(rd.glob("res_*")) + list(rd.glob("*.log")) + list(rd.glob("residual_history*")):
        f.unlink()
    shutil.copy(CASE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    mc.m50.write_species_db(rd)
    (rd / "probe.yaml").write_text(
        "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n", encoding="utf-8")

    # --- 入口分布 ---
    tabf = rd / "_inlet_table.txt"
    subprocess.run([sys.executable, str(HERE / "extract_inlet_table.py"),
                    "--run", a.inlet_run, "--x", str(a.inlet_x),
                    "--out", str(tabf), "--steady"], check=True)
    raw = np.loadtxt(tabf, skiprows=1)
    cols = tabf.read_text().splitlines()[0].split()
    tab = {c: raw[:, i] for i, c in enumerate(cols)}

    ymf = ", ".join(f"Y{i}: {gas.Y[sp]:.8f}" for i, sp in enumerate(SPECIES))
    bc = mc.BC.format(ro=rho, U=U, p=p, t=T, tw=a.tw, kinf=k_inf, ominf=om_inf,
                      ymf=ymf, ints="{inletProfile: 1}")
    bc += mc.BC_GAP.format(tw=a.tw) + BC_SYM
    (rd / "bcondConfig.yaml").write_text(bc, encoding="utf-8")
    patch_ic_profile(rd / "mesh.h5", tab, gas, a.tw)

    a_inf = float(np.sqrt(gas.gamma(T) * gas.R * T))
    common = dict(mu_inf=mu, lam_inf=gas.lam(T), cp_ref=gas.cp(T), gam_ref=gas.gamma(T),
                  relax=1.0, kinf=k_inf, ominf=om_inf, turb="sst",
                  ro_ref=rho, p_ref=p, a_ref=a_inf)
    stages = [("lam",  dict(conv=0, lim=0, cfl=0.2, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **{**common, "turb": "none"})),
              ("soft", dict(conv=0, lim=0, cfl=a.soft_cfl, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **common)),
              ("mid",  dict(conv=0, lim=0, cfl=a.mid_cfl, inner=10, nsteps=a.mid_steps,
                            out_int=a.mid_steps, **common))]
    for i, c in enumerate([float(x) for x in a.ramp.split(",") if x.strip()]):
        stages.append((f"ramp{i}_cfl{c:g}", dict(conv=1, lim=2, cfl=c, inner=6,
                                                 nsteps=a.ramp_steps,
                                                 out_int=a.ramp_steps, **common)))
    stages.append(("main", dict(conv=1, lim=2, cfl=a.cfl, inner=4, nsteps=a.main_steps,
                                out_int=a.out_int, **common)))
    (rd / "stages.json").write_text(json.dumps(
        [{**st, "tag": tag} for tag, st in stages], indent=2), encoding="utf-8")
    (rd / "case_setup.json").write_text(json.dumps(
        dict(tp1187_run=a.series, M=o["M"], T_inf=T, U_inf=U, rho_inf=rho, p_inf=p,
             mu_inf=mu, T_aw=o["T_aw"], Tt_c=o["Tt_c"], Tw=a.tw, mesh=a.mesh,
             inlet_run=a.inlet_run, inlet_x_m=a.inlet_x, dim="3D",
             cuda_blocksize=a.blocksize,
             phi=gas.phi, Y=gas.Y, R=gas.R),
        indent=2, ensure_ascii=False), encoding="utf-8")

    # **FP64 影アキュムレータを全段で入れる** (深部の float32 アーチファクトを消す)。
    cfg0 = mc.CFG.format(species=", ".join(SPECIES), **stages[0][1])
    (rd / "solverConfig.yaml").write_text(
        cfg0.replace("detectNaN: 1", "detectNaN: 1, qAccumulatorFP64: 1"), encoding="utf-8")
    # 入口分布の CSV は solverConfig.yaml が要る (gen が読む) のでここで作る
    subprocess.run([sys.executable, str(TOOLS / "gen_inlet_profile.py"), "gen",
                    "--run", str(rd), "--physID", "1", "--table", str(tabf),
                    "--axis", "y"], check=True)
    print(f"  段: {', '.join(t for t, _ in stages)}  -> {rd}")
    if a.dry:
        print("dry: 投入しない"); return
    for tag, kw in stages:
        cfg = mc.CFG.format(species=", ".join(SPECIES), **kw)
        (rd / "solverConfig.yaml").write_text(
            cfg.replace("detectNaN: 1", "detectNaN: 1, qAccumulatorFP64: 1"), encoding="utf-8")
        print(f"--- stage {tag}: cfl={kw['cfl']} conv={kw['conv']} turb={kw['turb']} "
              f"-> {kw['nsteps']} step")
        rc = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)]).returncode
        cur = rd / f"res_{kw['nsteps']}.h5"
        if rc != 0 or not cur.exists() or list(rd.glob("res_nan_*.h5")):
            raise SystemExit(f"stage {tag}: rc={rc} → 段の失敗。古い場を引き継がずに中断")
        for suf in ("forge_run.log", "residual_history.csv", "CONVERGENCE_VERDICT.txt"):
            src = rd / suf
            if src.exists():
                shutil.copy(src, rd / f"{Path(suf).stem}_{tag}{Path(suf).suffix}")
        mc.seed_turb(rd / "mesh.h5", k_inf, om_inf)


if __name__ == "__main__":
    main()
