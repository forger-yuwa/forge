#!/usr/bin/env python3
"""case/49 のメッシュ取り込みと段階起動 (plan §4.2 / §4.5 / §4.6)。

usage:
  python3 gen_runs.py convert --msh cad/forge.msh --out mesh/stageA.h5
  python3 gen_runs.py run --run run_0001_stageA_cpg --mesh mesh/stageA.h5 \
         --inlet-csv precursor/run_0002_.../inlet_profile_1.csv [--gas CPG] [--main-steps 20000]

設定はすべて manifest.json (setup.py --resolve) から。physID は manifest の `phys_id`。

重要 (plan §4.2 / §4.6, codex plan-1 M5):
  - **変換は最終形の壁タグで行う** (SST の wall_dist は no-slip 壁から作られる)。
    段階起動 S0 の「全面 slip」は実行時に bcondConfig.yaml を差し替えて実現する。
  - `space.pRef` = P_inf (非直交 float32 の自由流保存)。
  - IC はキャビティ内 (z<0) を **静止・壁温** から始める (一様 M5 のままだと大きな速度過渡)。
  - 同一メッシュの段間引き継ぎは index コピー (**3D で interp_field.py は使わない**)。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
BUILD = ROOT / "solver_density_cuda" / "build"
# FORGE_BIN は build-native に固定する (別セッションが build/ を使うため。
# 解像壁の qwall/utau 診断出力を入れたバイナリ)。
ENV = dict(os.environ,
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""),
           FORGE_BIN=os.environ.get("FORGE_BIN", str(ROOT / "solver_density_cuda" / "build-native" / "forge")))
sys.path.insert(0, str(HERE))
from setup import load as load_conditions  # noqa: E402

# manifest は環境変数で差し替えられる (塞ぎ形状 manifest_plug.json などの併存用)
sys.path.insert(0, str(HERE / "tools"))
import geom_common as gc  # noqa: E402

MAN = json.loads((HERE / os.environ.get("CASE49_MANIFEST", "manifest.json")).read_text())
PID, G = MAN["phys_id"], MAN["geometry"]
D = load_conditions()

# 群は manifest に**実在するものだけ**にする (塞ぎ形状 plug_cavity ではキャビティ壁が無い)
_iso = ["cav_outer", "cav_floor", "cyl_side"] + (
    ["cyl_top"] if D["cyl_top_thermal"] == "isothermal" else [])
_ad = ["plate", "plate_in"] + ([] if D["cyl_top_thermal"] == "isothermal" else ["cyl_top"])
WALLS_ISO = [g for g in _iso if g in PID]
WALLS_AD = [g for g in _ad if g in PID]
SLIPS = [g for g in (["top", "side", "sym"] + (["runup"] if G["has_runup"] else [])) if g in PID]


# ---------------------------------------------------------------- config
WALL_T_OVERRIDE = {}          # {壁グループ名: 壁温 [K]}。--wall-temps で埋める


def wall_T_of(group):
    """その壁グループの等温壁温度。**非一様壁温 (plan §5.1 #21) のための上書きに対応**。

    等温壁 3 本 (全壁同温) では熱回路の Q_i = G_i0(T_ext - T_i) + Σ_j G_ij(T_j - T_i) の
    第 2 項 (他の壁からガス越しに受ける熱) が恒等的に 0 になり同定できない
    (2026-09-19 codex)。壁ごとに違う温度を与えてこの項を測る。
    """
    return float(WALL_T_OVERRIDE.get(group, D["wall_T"]))


def bcond(stage, inlet_profile=False):
    """stage: 'convert'|'slip'|'adiabatic'|'isothermal'
    convert = 最終形 (wall_dist 用)。slip = S0 用に全壁 slip。"""
    L = []
    prof = "{inletProfile: 1}" if inlet_profile else ""
    L.append("inlet:     {physID: %d, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: %s, "
             "floats: {ro: %.8g, Ux: %.8g, Uy: 0.0, Uz: 0.0, Ps: %.8g, k: %.8g, omega: %.8g}}"
             % (PID["inlet"], prof, D["ro_inf"], D["U_inf"], D["P_inf"], D["k_inf"], D["omega_inf"]))
    L.append("outlet:    {physID: %d, kind: outlet_statPress, outputHDFflg: 0, ints: , "
             "floats: {Ps: %.8g, Pt: %.8g, Tt: %.8g}}"
             % (PID["outlet"], D["P_inf"], D["P_inf"], D["T_inf"]))
    for s in SLIPS:
        L.append("%-10s {physID: %d, kind: slip, outputHDFflg: 0, ints: , floats: }" % (s + ":", PID[s]))
    for w in WALLS_AD + WALLS_ISO:
        if stage == "slip":
            L.append("%-10s {physID: %d, kind: slip, outputHDFflg: 0, ints: , floats: }" % (w + ":", PID[w]))
        elif stage == "adiabatic" or (stage == "convert" and w in WALLS_AD) or \
                (stage in ("convert", "isothermal") and w in WALLS_AD):
            L.append("%-10s {physID: %d, kind: wall, outputHDFflg: 1, ints: , "
                     "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}" % (w + ":", PID[w]))
        else:                                  # convert / isothermal の等温壁
            L.append("%-10s {physID: %d, kind: wall_isothermal, outputHDFflg: 1, ints: , "
                     "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: %.8g}}" % (w + ":", PID[w], wall_T_of(w)))
    return "\n".join(L) + "\n"


def urans_cfg(nsteps, dt, *, cfl_pseudo=12.0, nsub=10, outint=200, gas=None, mesh_file="mesh.h5"):
    """dual-time の非定常設定 (plan §4.9(b))。
    **`deltaT.control: 0` (固定物理 dt) が必須** — 定常の control:1 を継承すると
    main.cpp が例外終了する。`dualTime: 1` と `blockDPLUR: 1` も必須。
    **`implicitRelax` は書かない** (物理時間項が安定性を担うので緩和は遅いモードを残すだけ)。"""
    base = solver_cfg(nsteps, cfl_pseudo, outint=outint, gas=gas, mesh_file=mesh_file)
    base = base.replace("  unsteady: 0\n  dualTime: 0", "  unsteady: 1\n  dualTime: 1")
    base = base.replace("control: 1, dt: 1e-9", "control: 0, dt: %.6g" % dt)
    base = base.replace(", implicitRelax: 0.7", "")
    base = base.replace("  nStepInner: 4", "  nStepInner: 4\n  nSubIterDualTime: %d" % nsub)
    return base


def perturb_asym(h5, amp=0.05):
    """キャビティ気体に**一様な周方向旋回**を与える (plan §4.9(a))。

    半割 (y=0 対称面) が禁じるのは「鏡像で符号が変わるモード」で、環状キャビティでは
    その代表が**周方向循環 (swirl)** である。鏡像は旋回の向きを反転させるので、
    対称解では循環が**厳密に 0**。したがって
      「旋回を与えて、それが減衰して 0 に戻るか」
    が半割の可否をそのまま判定する。**この指標はメッシュの左右非対称に汚染されない**
    (以前は開口の k を ±1 % だけ触っていたが、T に効かず指標が動かなかった)。

    amp はキャビティ気体の rms 速さに対する相対振幅。運動量を足し、対応する運動
    エネルギーを roe にも足して**温度を変えない**。
    """
    with h5py.File(h5, "r+") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        ro = np.array(f["/VALUE/ro"])
        rux, ruy = np.array(f["/VALUE/roUx"]), np.array(f["/VALUE/roUy"])
        ruz = np.array(f["/VALUE/roUz"])
        roe = np.array(f["/VALUE/roe"])
        m = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], MAN, shrink=0.0)
        if m.sum() < 100:
            print("  キャビティ節点が取れないので擾乱なし"); return
        u = np.stack([rux[m] / ro[m], ruy[m] / ro[m], ruz[m] / ro[m]], axis=1)
        uref = float(np.sqrt(np.mean(np.sum(u ** 2, axis=1))))
        r = np.hypot(c[m, 0], c[m, 1])
        # e_theta = (-y, x)/r  (右ねじ方向の一様旋回)
        et = np.stack([-c[m, 1] / np.maximum(r, 1e-12), c[m, 0] / np.maximum(r, 1e-12)], axis=1)
        du = amp * uref
        k0 = 0.5 * (rux[m] ** 2 + ruy[m] ** 2 + ruz[m] ** 2) / ro[m]
        rux[m] += ro[m] * du * et[:, 0]
        ruy[m] += ro[m] * du * et[:, 1]
        k1 = 0.5 * (rux[m] ** 2 + ruy[m] ** 2 + ruz[m] ** 2) / ro[m]
        roe[m] += (k1 - k0)                      # 内部エネルギー (= 温度) を保つ
        for nm, arr in (("roUx", rux), ("roUy", ruy), ("roe", roe)):
            f["/VALUE/" + nm][:] = arr.astype(f["/VALUE/" + nm].dtype)
        print("  一様旋回擾乱: %d ノードに u_theta = %.3g m/s (rms 速さ %.3g の %.0f %%)"
              % (int(m.sum()), du, uref, 100 * amp))


def solver_cfg(nsteps, cfl, *, conv=1, lim=2, ninner=4, relax=0.7, outint=2000,
               model="sst", mesh_file="mesh.h5", gas=None, extra_fields=None):
    # **診断勾配の追加出力** (残作業 #28)。`dTd*` は `output_cellValNames` に無く出せないので、
    # TP (単一擬似種 = R 一定) の EOS で ∇T = ∇P/(ρR) − P∇ρ/(ρ²R) と再構成するために
    # `dPd*` / `drod*` を、粘性仕事のために速度勾配を足す。level 2 より軽い。
    extra = ("" if not extra_fields else ", extraFields: [%s]" % ", ".join(extra_fields))
    gas = (gas or D["gas"]).upper()
    kc = D["mu_inf"] * D["cp"] / D["prandtl_lam"]
    if gas == "TP":
        # semi-perfect 乾燥空気 1 擬似種 (plan §4.10)。`thermoHrefTemp` は必須
        # (絶対基準 h だと chi_eos が桁違いになる)。cp/gamma も parser が要求するので残す。
        phys = ('physProp: {thermalMethod: 2, species: ["MIXDRY"], speciesDBFile: "species_db.yaml", '
                'thermoHrefTemp: 298.15, viscMethod: 1, visc: %.8g, thermCond: %.8g, thermCondMethod: 1, '
                'prandtlLam: %s, cp: %s, gamma: %s}'
                % (D["mu_inf"], kc, D["prandtl_lam"], D["cp"], D["gamma"]))
    else:
        phys = ('physProp: {thermalMethod: 0, viscMethod: 1, visc: %.8g, thermCond: %.8g, '
                'thermCondMethod: 1, prandtlLam: %s, cp: %s, gamma: %s}'
                % (D["mu_inf"], kc, D["prandtl_lam"], D["cp"], D["gamma"]))
    turb = ('turbulence: {model: "none"}' if model == "none" else
            'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 0, '
            'wallTreatmentSST: 0, turbulentPrandtl: %.4g, kInit: %.8g, omegaInit: %.8g}'
            % (D["prandtl_turb"], D["k_inf"], D["omega_inf"]))
    return f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "{mesh_file}", valueFileName: "{mesh_file}"}}
gpu: 1
solver: "SLAU"
{phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-9, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-11, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}, pRef: {D['P_inf']:.8g}}}
{turb}
output: {{level: 1{extra}}}
initial: "uniform_p101325_u10"
"""


# ---------------------------------------------------------------- 取り込み
def cmd_convert(a):
    conv = HERE / "mesh" / "_conv"
    conv.mkdir(parents=True, exist_ok=True)
    # 変換は幾何と wall_dist だけなので EOS に依らない -> CPG で通す (species DB 不要)
    (conv / "solverConfig.yaml").write_text(solver_cfg(10, 0.5, conv=0, lim=0, ninner=5,
                                                       outint=10, mesh_file="m.h5", gas="CPG"))
    (conv / "bcondConfig.yaml").write_text(bcond("convert"))
    r = subprocess.run([str(BUILD / "convertGmshToForge"), str(Path(a.msh).resolve()), "m.h5"],
                       cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / "convert.log").write_text(r.stdout + r.stderr)
    if not (conv / "m.h5").exists():
        print(r.stdout[-3000:], r.stderr[-2000:])
        raise SystemExit("convert failed")
    out = Path(a.out if os.path.isabs(a.out) else HERE / a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(conv / "m.h5"), str(out))
    for x in conv.glob("m.xmf"):
        x.unlink()
    for line in (conv / "convert.log").read_text().splitlines():
        if any(k in line for k in ("closure", "dual", "nBconds", "Number of", "bname")):
            print("   ", line)
    q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(out)],
                       capture_output=True, text=True)
    (out.parent / ("quality_%s.txt" % out.stem)).write_text(q.stdout + q.stderr)
    print("\n".join(l for l in q.stdout.splitlines() if "VERDICT" in l or "max" in l.lower())[:2000])
    print("wrote", out)


# ---------------------------------------------------------------- IC
def _tp_energy(T, ro, u2):
    """TP (semi-perfect 1 擬似種) の roe。datum は thermoHrefTemp=298.15
    (h(298.15)=0 の基準)。e = h(T) - h(Tref) - R T。"""
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.gas.semiperfect import GasSemiPerfect
    g = GasSemiPerfect(D["dry_air_Y"], Tt=2000.0)
    href = float(g.h_mass(298.15))
    h = np.array([float(g.h_mass(t)) for t in np.atleast_1d(T)]) - href
    e = h - D["R_tp"] * np.atleast_1d(T)
    return ro * (e + 0.5 * u2)


def patch_ic(h5, inlet_csv, gas="CPG"):
    """外部流は入口 BL 分布を z で写す。キャビティ内 (z<0) は静止・壁温・P_inf。"""
    prof = np.genfromtxt(inlet_csv, names=True)
    zc = np.asarray(prof["z"], float)
    o = np.argsort(zc)
    zc = zc[o]
    get = lambda n: np.asarray(prof[n], float)[o]      # noqa: E731
    with h5py.File(h5, "r+") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        n = f["/VALUE/ro"].shape[0]
        z = c[:, 2]
        ro = np.interp(z, zc, get("ro"))
        ux = np.interp(z, zc, get("Ux"))
        uz = np.interp(z, zc, get("Uz"))
        ps = np.interp(z, zc, get("Ps"))
        kk = np.interp(z, zc, get("k"))
        om = np.interp(z, zc, get("omega"))
        uy = np.zeros(n)
        incav = z < 0.0                                # キャビティ内は静止・壁温
        ro[incav] = D["P_inf"] / (D["R"] * D["wall_T"])
        ux[incav] = 0.0; uz[incav] = 0.0
        ps[incav] = D["P_inf"]
        kk[incav] = D["k_inf"] * 1e-3; om[incav] = D["omega_inf"]
        wd = np.array(f["/VALUE/wall_dist"]) if "/VALUE/wall_dist" in f else np.ones(n)
        wall = wd <= 0.0
        ux[wall] = 0.0; uy[wall] = 0.0; uz[wall] = 0.0
        u2 = ux ** 2 + uy ** 2 + uz ** 2
        if gas.upper() == "TP":
            # T を (ps, ro) から作り直して TP の e で roe を組む (CPG の ps/(γ-1) は datum が違う)
            Tfield = ps / np.maximum(ro * D["R_tp"], 1e-30)
            # 一意な温度だけ評価して展開 (NASA-9 の評価は node 数ぶん回すと遅い)
            Tq = np.round(Tfield, 2)
            uniq, inv = np.unique(Tq, return_inverse=True)
            sys.path.insert(0, str(ROOT / "design"))
            from forge_design.gas.semiperfect import GasSemiPerfect
            g = GasSemiPerfect(D["dry_air_Y"], Tt=2000.0)
            # **numpy 2 系では 0 次元でない配列の float() がエラーになる**ので、
            # setup.py と同じ `_f` で明示的にスカラー化する (AWS の numpy で実際に落ちる)。
            from setup import _f
            href = _f(g.h_mass(298.15))
            e_u = np.array([_f(g.h_mass(t)) - href - D["R_tp"] * t for t in uniq])
            roe = ro * (e_u[inv] + 0.5 * u2)
            print("  IC(TP): T %.1f..%.1f K, 一意温度 %d 点" % (Tfield.min(), Tfield.max(), len(uniq)))
        else:
            roe = ps / (D["gamma"] - 1.0) + 0.5 * ro * u2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * ux).astype(np.float32)
        f["/VALUE/roUy"][:] = (ro * uy).astype(np.float32)
        f["/VALUE/roUz"][:] = (ro * uz).astype(np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, v in (("roK", ro * kk), ("roOmega", ro * om)):
            ds = "/VALUE/" + name
            if ds in f:
                f[ds][:] = v.astype(np.float32)
            else:
                f.create_dataset(ds, data=v.astype(np.float32))
        print("IC: n=%d  cavity nodes=%d  wall nodes=%d" % (n, int(incav.sum()), int(wall.sum())))


def index_copy(res, mesh):
    """同一メッシュの段間引き継ぎ: res の VALUE を mesh.h5 に index コピー (wall_dist は残す)。"""
    with h5py.File(res, "r") as s, h5py.File(mesh, "r+") as d:
        for k in s["VALUE"]:
            if k == "wall_dist":
                continue
            arr = np.array(s["VALUE/" + k])
            ds = "VALUE/" + k
            if ds in d:
                if d[ds].shape != arr.shape:
                    continue
                d[ds][...] = arr.astype(d[ds].dtype)
            else:
                d.create_dataset(ds, data=arr)


def run_forge(rd):
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


_SM = {}      # run ごとの StageManifest


def stage(rd, tag, cfgtext, bctext, nsteps, keep=False):
    # **段ごとの実効設定を記録する** (収束判定の区間を段名でなくデータで決めるため。
    #  2026-09-19 codex: 段名のプレフィックスでは方程式・BC・離散化の同一性を保証できない)
    sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
    from stage_manifest import StageManifest
    sm = _SM.setdefault(str(rd), StageManifest(rd))
    sm.add(tag, cfgtext, bctext)
    sm.write()
    (rd / "solverConfig.yaml").write_text(cfgtext)
    (rd / "bcondConfig.yaml").write_text(bctext)
    rc = run_forge(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    # 完走判定は residual_history の最終 step で行う (出力間隔と nStepOuter が割り切れないと
    # 最後の res が出ないため)。段間引き継ぎに使う最終場が要るので outint=nsteps にしてある。
    last_step = -1
    rh = rd / "residual_history.csv"
    if rh.exists():
        tail = rh.read_text().strip().splitlines()
        if len(tail) > 1:
            last_step = int(float(tail[-1].split(",")[0]))
    ok = res and int(res[-1].stem.split("_")[1]) >= nsteps and last_step >= nsteps - 2
    print("  [%s] rc=%d last_res=%s last_step=%d" %
          (tag, rc, res[-1].name if res else "-", last_step), flush=True)
    if not ok:
        raise SystemExit("stage %s failed (rc=%d)" % (tag, rc))
    index_copy(res[-1], rd / "mesh.h5")
    if not keep:
        for f in rd.glob("res_*"):
            f.unlink()
    for nm in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
        if (rd / nm).exists():
            shutil.move(str(rd / nm), str(rd / ("%s_%s%s" % (Path(nm).stem, tag, Path(nm).suffix))))


def cmd_run(a):
    rd = HERE / a.run
    if rd.exists():
        raise SystemExit("%s exists" % rd)
    rd.mkdir()
    shutil.copy(HERE / a.mesh if not os.path.isabs(a.mesh) else a.mesh, rd / "mesh.h5")
    shutil.copy(a.inlet_csv, rd / ("inlet_profile_%d.csv" % PID["inlet"]))
    gas = (a.gas or D["gas"]).upper()
    if gas == "TP":
        sdb = HERE / "species_db.yaml"
        if not sdb.exists():
            raise SystemExit("species_db.yaml が無い (setup で生成: mixture_pseudo_species)")
        shutil.copy(sdb, rd / "species_db.yaml")
    print("  gas =", gas)
    (rd / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    # **使った manifest を run に残す** (run を自己記述にする)。偏心スイープのように幾何が
    # run ごとに違うとき、後から共有 manifest.json で評価すると別の CV マスクになる。
    shutil.copy(HERE / os.environ.get("CASE49_MANIFEST", "manifest.json"), rd / "manifest.json")
    # **解決済みの作動条件も run に固定する** (2026-09-19 codex Major 10: 評価側が共有 case.json を
    # 読み直していたため、M や壁温を変えると**過去 run の評価が変わって**しまっていた)。
    # 使った EOS も書く (回復温度 Taw を EOS に合わせて選ぶため)。
    if getattr(a, "wall_temps", ""):
        for kv in a.wall_temps.split(","):
            g, _, t = kv.partition("=")
            g = g.strip()
            if g not in PID:
                raise SystemExit("--wall-temps: 未知の壁グループ %r (候補 %s)"
                                 % (g, ",".join(sorted(PID))))
            WALL_T_OVERRIDE[g] = float(t)
        print("  壁温上書き: " + ", ".join("%s=%.2f K" % (k, v)
                                           for k, v in sorted(WALL_T_OVERRIDE.items())))
    cond = dict(D)
    if WALL_T_OVERRIDE:
        cond["wall_T_by_group"] = dict(WALL_T_OVERRIDE)
    cond["gas_used"] = gas
    (rd / "conditions.json").write_text(json.dumps(cond, indent=2, ensure_ascii=False,
                                                   default=float))
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    if a.ic_from:
        # **引き継ぎ元と条件が一致しているか必ず検査する** (2026-09-20)。
        # `CASE49_CASE` / `CASE49_MANIFEST` を export し忘れると既定の `case.json`
        # (別のマッハ数・別の壁温) で回ってしまい、IC だけ前 run という無意味な計算が
        # 静かに完走する。実際に M9/壁温 20・1000 degC のつもりで M5/500 degC を 2 本回した。
        pj = HERE / a.ic_from / "conditions.json"
        if pj.exists():
            import json as _json
            par = _json.loads(pj.read_text())
            keys = ("mach", "wall_T", "gas", "T_inf", "P_inf", "cyl_top_thermal",
                    "wall_T_by_group")
            bad = [(k, par.get(k), D.get(k)) for k in keys
                   if k in par and k in D
                   and (abs(par[k] - D[k]) > 1e-9 * max(1.0, abs(par[k]))
                        if isinstance(par[k], (int, float)) and isinstance(D[k], (int, float))
                        else par[k] != D[k])]
            if bad and getattr(a, "ic_force", ""):
                print("  [条件不一致だが --ic-force で続行: %s]" % a.ic_force)
                for k, x, y in bad:
                    print("     %-16s 親 %s  ->  今回 %s" % (k, x, y))
                bad = []
            if bad:
                raise SystemExit(
                    "引き継ぎ元 %s と条件が違う (CASE49_CASE / CASE49_MANIFEST の export 漏れ?):\n"
                    % a.ic_from
                    + "\n".join("  %-16s 親 %s  !=  今回 %s" % (k, x, y) for k, x, y in bad))
        else:
            print("  [警告] %s/conditions.json が無く、条件の一致を検査できない" % a.ic_from)
        src = sorted((HERE / a.ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int(f.stem.split("_")[1]))[-1]
        index_copy(src, rd / "mesh.h5")
        (rd / "CONTINUED_FROM").write_text(str(src) + "\n")
        print("  IC: %s から index コピー" % src.name)
    else:
        patch_ic(rd / "mesh.h5", a.inlet_csv, gas=gas)
    if a.urans:
        if a.perturb > 0:
            perturb_asym(rd / "mesh.h5", a.perturb)
        cfgu = urans_cfg(a.main_steps, a.dt, cfl_pseudo=a.cfl_pseudo, nsub=a.nsub,
                         outint=a.out_int, gas=gas)
        (rd / "solverConfig.yaml").write_text(cfgu)
        (rd / "bcondConfig.yaml").write_text(bcond("isothermal", inlet_profile=True))
        if a.dry:
            return
        rc = run_forge(rd)
        print("urans rc", rc)
        print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-600:])
        return
    if a.dry:
        (rd / "solverConfig.yaml").write_text(solver_cfg(a.main_steps, a.cfl, gas=gas))
        (rd / "bcondConfig.yaml").write_text(bcond("isothermal", inlet_profile=True))
        return
    if a.main_only:
        # **本段だけを回す** (収束済みの場を引き継いで伸ばす用)。`--ic-from` だけでは
        # 段階起動をやり直してしまい、**S0 の全 slip が解を壊す**ので専用の入口を用意する。
        if not a.ic_from:
            raise SystemExit("--main-only は --ic-from と併用する (引き継ぐ場が要る)")
        cfgtext = solver_cfg(a.main_steps, a.cfl, outint=a.out_int, gas=gas, relax=a.relax,
                             extra_fields=(["dPdx", "dPdy", "dPdz", "drodx", "drody", "drodz",
                                            "dUxdz", "dUydz", "dUzdz", "dUzdx", "dUzdy"]
                                           if a.grad_out else None))
        bctext = bcond("isothermal", inlet_profile=True)
        # **継続 run にも stage_manifest.json を書く**。以前は書いていなかったため
        # `check_convergence.py --segment` が「区間を決められない」で落ち、ゲートは
        # 弱い代替判定 (上昇の有無) に落ちていた (2026-09-19)。段は本段 1 つだけ。
        sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
        from stage_manifest import StageManifest
        sm = StageManifest(rd)
        sm.add("S6_main_ext", cfgtext, bctext)
        sm.write()
        (rd / "solverConfig.yaml").write_text(cfgtext)
        (rd / "bcondConfig.yaml").write_text(bctext)
        rc = run_forge(rd)
        print("main-only rc", rc)
        print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-600:])
        return
    ip = dict(inlet_profile=True)
    # S0: 全壁 slip・層流・1 次 (キャビティ内圧の平衡化)
    sc = float(getattr(a, "startup_cfl_scale", 1.0))
    if sc != 1.0:
        print("  段階起動の CFL を x%.3g" % sc)
    stage(rd, "S0_slip", solver_cfg(2000, 0.5 * sc, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("slip", **ip), 2000)
    # S1: no-slip 断熱 (層流)
    stage(rd, "S1_lam", solver_cfg(2000, 0.3 * sc, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("adiabatic", **ip), 2000)
    # S2: キャビティ等温壁 (層流)
    stage(rd, "S2_iso", solver_cfg(2000, 0.5 * sc, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("isothermal", **ip), 2000)
    # S3/S4: SST soft -> mid (1 次)
    stage(rd, "S3_sst_soft", solver_cfg(6000, 0.3 * sc, conv=0, lim=0, ninner=10, outint=3000, gas=gas), bcond("isothermal", **ip), 6000)
    stage(rd, "S4_sst_mid", solver_cfg(3000, 1.0 * sc, conv=0, lim=0, ninner=10, outint=3000, gas=gas), bcond("isothermal", **ip), 3000)
    # S5: 2 次ランプ
    for i, cv in enumerate([float(v) for v in a.ramp.split(",") if v]):
        stage(rd, "S5_ramp%d_cfl%g" % (i, cv), solver_cfg(2000, cv * sc, outint=2000, gas=gas), bcond("isothermal", **ip), 2000)
    # S6: 本段
    main_cfg = solver_cfg(a.main_steps, a.cfl, outint=a.out_int, gas=gas)
    main_bc = bcond("isothermal", **ip)
    sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
    from stage_manifest import StageManifest
    sm = _SM.setdefault(str(rd), StageManifest(rd))
    sm.add("S6_main", main_cfg, main_bc, history="residual_history.csv")
    sm.write()
    (rd / "solverConfig.yaml").write_text(main_cfg)
    (rd / "bcondConfig.yaml").write_text(main_bc)
    rc = run_forge(rd)
    print("main rc", rc)
    print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-900:])


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert")
    c.add_argument("--msh", required=True)
    c.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--run", required=True)
    r.add_argument("--mesh", required=True)
    r.add_argument("--inlet-csv", required=True)
    r.add_argument("--main-steps", type=int, default=20000)
    r.add_argument("--cfl", type=float, default=2.0)
    r.add_argument("--out-int", type=int, default=2000)
    r.add_argument("--ramp", default="0.5,1,2")
    r.add_argument("--gas", default=None, choices=["CPG", "TP"], help="既定は case.json の gas")
    r.add_argument("--urans", action="store_true", help="定常場から dual-time の非定常へ (plan §4.9)")
    r.add_argument("--dt", type=float, default=3.0e-7, help="URANS の物理 dt [s]")
    r.add_argument("--nsub", type=int, default=10)
    r.add_argument("--cfl-pseudo", type=float, default=12.0)
    r.add_argument("--ic-from", default=None, help="定常場の run (mesh.h5 を index コピーで引き継ぐ)")
    r.add_argument("--perturb", type=float, default=0.01, help="URANS の左右非対称擾乱 (相対)")
    r.add_argument("--startup-cfl-scale", type=float, default=1.0,
                   help="段階起動 (S0-S5) の CFL に掛ける係数。**第一層を薄くすると SST 投入段 "
                        "(S3) で omega 壁値が y1^-2 で増えて落ちる** ため 0.3-0.5 で緩める")
    r.add_argument("--grad-out", action="store_true",
                   help="開口収支をソルバ勾配で組むための診断場を追加出力する (#28)")
    r.add_argument("--relax", type=float, default=0.7,
                   help="implicitRelax。**1.0 = 緩和なし** (解が緩和で動いていないかの確認用)。"
                        "AGENTS の SERN 系は relax を使わない方針なので感度を取る")
    r.add_argument("--ic-force", default="",
                   help="引き継ぎ元と条件が違っても進める理由 (意図的な壁温変更など)。"
                        "空なら不一致で停止する")
    r.add_argument("--wall-temps", default="",
                   help="壁ごとの等温壁温度を上書き (例 cav_outer=1273.15,cyl_side=1273.15,cav_floor=293.15)")
    r.add_argument("--main-only", action="store_true",
                   help="段階起動を飛ばして本段だけ回す (--ic-from 必須。収束済みの場を伸ばす用)")
    r.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    (cmd_convert if a.cmd == "convert" else cmd_run)(a)


if __name__ == "__main__":
    main()
