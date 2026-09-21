#!/usr/bin/env python3
"""case/51 の run 生成 — T2' (NASA TN D-8233 平板 + 単独横すきま) の 2D 段階起動 run。

`case/56/tools/make_case.py` を雛型にし、**乾燥空気 CPG** 向けに単純化した
(多成分 TP の `species_db.yaml`・組成 `Y{i}` 列・NASA-9 経路が要らない)。
case/56 から引き継いだ仕組み:

  - 段階起動 (層流暖機 → soft → mid → 2 次ランプ → 本段)。既定値は case/55 の T3 実績
    (`case/55/_t3_chain.sh`) に合わせる。
  - **段間の引き継ぎはその段が書いた `res_<nsteps>.h5` を名前で取る** (2026-09-20 codex
    result M9: 全段の最大番号を取ると mid 3000 step が ramp 2500 step より大きく、
    ramp/本段が古い mid の場から再開して ramp 段が空振りしていた)。
  - **リミッタ基準値を自由流で固定する** (同 M10: auto だと開始場依存の作用素になり、
    分割実行が連続実行と一致しない)。
  - 段の失敗 (rc≠0 / 最終 res が無い / `res_nan_*`) では古い場を引き継がずに中断する。

判定区間は `stage_manifest.json` (`solver_density_cuda/tools/stage_manifest.py`) に書く。
段名ではなく**方程式・BC・空間離散化の同一性**で区間が決まるので、収束判定は
`check_convergence.py <run> --segment` を使うこと。

入力は case/51 に既にある台帳から読む (数値をハードコードしない):

  `derived.json`    Gate A の自由流 (既定 Re' = 1.47e6/m の系列)
  `conditions.json` M・模型寸法 (記録用)
  `geometry.json`   領域高さ H (入口 CSV の上端)
  `mesh/<tag>.h5`   `gen_mesh.py` が作った node 用メッシュ

入口 BL は `tools/gen_inlet_csv.py` をサブプロセスで呼んで `inlet_profile_1.csv` にする。
**`inlet_bl.json` は訂正前 (θ 誤読値で拘束した) の台帳なので使わない** — 理由は
`gen_inlet_csv.py` の docstring と README「訂正後の台帳」を参照。

usage:
  python3 tools/make_case.py --run run_0001_t2p_re147_tw300 --mesh t2p_y30um_g101_r1025 --dry
  python3 tools/make_case.py --run run_0001_t2p_re147_tw300 --mesh t2p_y30um_g101_r1025
  python3 tools/make_case.py --run run_0002_t2p_re147_tr161 --plate-Tw 483.0 # 壁温比 T_surf/T_gap=1.61
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path
import numpy as np, h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))
from stage_manifest import StageManifest                          # noqa: E402

# 乾燥空気 CPG。`gen_mesh.py` の変換 config と同じ値に揃える (physProp を 2 か所で持たない)。
CP, GAMMA, PRANDTL_LAM = 1004.5, 1.4, 0.72
R_GAS = CP * (GAMMA - 1.0) / GAMMA              # = 287.0 J/(kg·K)
# EOS 正値化フロア。既定 (pMin 1.0 Pa, roMin 1e-4, tMin 1e-4) は大気スケール想定であり、
# 本試験 (p∞ = 57 Pa, ρ∞ = 3.8e-3 kg/m³) では pMin が運転圧の 1/57、roMin が ρ∞ の 1/38 と
# 近すぎて物理的に正しい状態をクランプしうる。自由流比 1e-4 級になるよう下げる。
PMIN, ROMIN, TMIN = 1.0e-2, 1.0e-6, 1.0e-4

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp:
  thermalMethod: 0
  viscMethod: 1
  thermCondMethod: 1
  prandtlLam: {pr:.3f}
  visc: {mu_inf:.6e}
  thermCond: {kcond:.6e}
  cp: {cp:.2f}
  gamma: {gam:.3f}
  pMin: {pmin:.6g}
  roMin: {romin:.6g}
  tMin: {tmin:.6g}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1.0e-9, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1,
    dt_min: 1.0e-12, dt_max: 1.0, detectNaN: 1, monitorInterval: 100}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: {inner}
space: {{convMethod: {conv}, limiter: {lim}, limiterRoRef: {ro_ref:.10g}, limiterPRef: {p_ref:.10g}, limiterARef: {a_ref:.10g}}}
turbulence: {{model: "{turb}", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInit: {kinf:.6g}, omegaInit: {ominf:.6g}}}
initial: "uniform_p101325_u10"
output: {{level: 1, extraFields: [vis_lam, vis_turb]}}
"""

# physID は `gen_mesh.py` の Physical Curve と一致させる (1 inlet / 2 outlet / 3 top /
# 4 plate / 5 buffer / 6 gap)。入口は `inletProfile: 1` で CSV の分布が floats を上書きする。
BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {{inletProfile: 1}}, floats: {{ro: {ro:.8e}, Ux: {u:.4f}, Uy: 0.0, Uz: 0.0, Ps: {p:.6f}, k: {kinf:.6g}, omega: {ominf:.6g}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6f}, Pt: {p:.6f}, Tt: {tt:.4f}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
plate:  {{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw_plate:.2f}}}}}
buffer: {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
gap:    {{physID: 6, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw_gap:.2f}}}}}
"""


def patch_ic(h5, st, tw_plate, tw_gap, k_inf, om_inf, ic_delta):
    """一様自由流 + 壁近傍の速度ランプを input h5 に書く (case/56 `patch_ic` と同じ考え方)。

    CPG なので $\\rho e = p/(\\gamma-1) + \\rho u^2/2$ で直接組む。圧力は全域 p∞ 一様、
    壁ノード (wall_dist=0) だけ u=0・T=T_w にして密度を EOS で合わせる (等温 BC が step 0
    から効く)。すきま内部は wall_dist が小さいのでランプが自動的に u≈0 にする。
    """
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:].astype(float)
        u = st["U"] * np.tanh(np.maximum(wd, 0.0) / ic_delta)
        ro = np.full(n, st["ro"])
        wall = wd <= 0.0
        ingap = c[:, 1] < -1e-9                      # y<0 の壁ノード = すきま壁 (前壁/床/後壁)
        u[wall] = 0.0
        ro[wall & ~ingap] = st["p"] / (R_GAS * tw_plate)
        ro[wall & ingap] = st["p"] / (R_GAS * tw_gap)
        roe = st["p"] / (GAMMA - 1.0) + 0.5 * ro * u ** 2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros(n, np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, val in (("roK", k_inf), ("roOmega", om_inf)):
            key = f"/VALUE/{name}"
            v = (ro * val).astype(np.float32)
            if key in f:
                f[key][:] = v
            else:
                f.create_dataset(key, data=v)
        print(f"IC: n={n}, 壁ノード {int(wall.sum())} (うちすきま {int((wall & ingap).sum())}), "
              f"T∞={st['T']:.2f} K, U∞={st['U']:.1f} m/s, ramp δ_ic={ic_delta*1e3:.2f} mm, "
              f"Tw_plate={tw_plate:.1f} K, Tw_gap={tw_gap:.1f} K")


def seed_turb(h5, k_inf, om_inf):
    """roK / roOmega が無い or 全ゼロなら自由流値で埋める。

    層流暖機段の出力には **roK/roOmega がそもそも存在しない**ので、そのまま SST 段に渡すと
    k=0 から始まり初手で発散する (case/56 で 2026-09-20 に踏んだ)。段の引き継ぎ後に必ず呼ぶ。
    """
    with h5py.File(h5, "r+") as f:
        ro = f["/VALUE/ro"][:].astype(float)
        for name, val in (("roK", k_inf), ("roOmega", om_inf)):
            key = f"/VALUE/{name}"
            v = (ro * val).astype(np.float32)
            if key not in f:
                f.create_dataset(key, data=v); print(f"  seed {name} (新規)")
            elif float(np.max(np.abs(f[key][:]))) == 0.0:
                f[key][:] = v; print(f"  seed {name} (全ゼロだった)")


def build_stages(a, common):
    """(tag, cfg kwargs) の列。case/55 T3 と同じ並び: lam → soft → mid → ramp* → main。"""
    stages = [("lam",  dict(conv=0, lim=0, cfl=a.soft_cfl, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **{**common, "turb": "none"})),
              ("soft", dict(conv=0, lim=0, cfl=a.soft_cfl, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **common)),
              ("mid",  dict(conv=0, lim=0, cfl=a.mid_cfl, inner=10, nsteps=a.mid_steps,
                            out_int=a.mid_steps, **common))]
    for i, c in enumerate([float(x) for x in a.ramp.split(",") if x.strip()]):
        stages.append((f"ramp{i}_cfl{c:g}", dict(conv=1, lim=2, cfl=c, inner=6,
                                                 nsteps=a.ramp_steps, out_int=a.ramp_steps, **common)))
    stages.append(("main", dict(conv=1, lim=2, cfl=a.cfl, inner=4, nsteps=a.main_steps,
                                out_int=a.out_int, **{**common, "relax": a.relax})))
    return stages


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", default="t2p_y30um_g101_r1025",
                    help="mesh/<tag>.h5 (本命 t2p_y30um_g101_r1025 / 格子感度 t2p_y30um_g51_coarse)")
    ap.add_argument("--re", type=float, default=1.47e6, help="derived.json の Re'/m 系列")
    ap.add_argument("--plate-Tw", type=float, default=300.0,
                    help="平板 (physID 4) の壁温 T_surf [K]。入口 BL の再構成にも使う")
    ap.add_argument("--gap-Tw", type=float, default=300.0,
                    help="すきま壁 (physID 6) の壁温 T_gap [K]。壁温比は T_surf/T_gap")
    ap.add_argument("--ic-delta", type=float, default=2.0e-3, help="IC の壁近傍ランプ幅 [m]")
    ap.add_argument("--tu", type=float, default=0.005, help="自由流乱れ度 (k∞ = 1.5 (Tu U)²)")
    ap.add_argument("--mut-ratio", type=float, default=10.0, help="自由流 μ_t/μ (ω∞ = ρk/(比·μ))")
    ap.add_argument("--soft-steps", type=int, default=3000)
    ap.add_argument("--mid-steps", type=int, default=3000)
    ap.add_argument("--ramp", default="0.3,0.6,1.0")
    ap.add_argument("--ramp-steps", type=int, default=2500)
    ap.add_argument("--soft-cfl", type=float, default=0.2)
    ap.add_argument("--mid-cfl", type=float, default=0.5)
    ap.add_argument("--cfl", type=float, default=1.2)
    ap.add_argument("--relax", type=float, default=0.7, help="本段の implicitRelax")
    ap.add_argument("--main-steps", type=int, default=400000,
                    help="本段の step 数 (case/55 のキャビティが 400k step で静定した実績)")
    ap.add_argument("--out-int", type=int, default=5000)
    ap.add_argument("--prepare-only", action="store_true",
                    help="run ディレクトリ一式を作るだけで forge を起動しない")
    ap.add_argument("--dry", action="store_true",
                    help="--prepare-only に加えて生成した config を表示する")
    a = ap.parse_args()

    # --- 台帳 (case/51 の既存ファイル) から自由流・幾何を読む ---
    der = json.loads((CASE / "derived.json").read_text(encoding="utf-8"))
    try:
        o = next(x for x in der["series"] if abs(x["Re_m"] - a.re) / a.re < 1e-6)
    except StopIteration:
        raise SystemExit(f"derived.json に Re_m={a.re:g} が無い "
                         f"(あるのは {[x['Re_m'] for x in der['series']]})")
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    geo = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    assert geo["units"] == "mm", f"geometry.json の単位が mm でない: {geo['units']}"
    H = geo["tunnel"]["half_height"] * 1e-3
    T, U, ro, p, mu = o["T_inf"], o["U_inf"], o["ro_inf"], o["p_inf"], o["mu_inf"]
    st = dict(T=T, U=U, ro=ro, p=p, mu=mu)

    # CPG EOS (R = cp(γ-1)/γ) と Gate A の (p, T, ρ) の整合を確認する。
    ro_eos = p / (R_GAS * T)
    d_ro = abs(ro_eos - ro) / ro
    print(f"=== T2' 自由流 (Re' = {o['Re_m']:.3e}/m, M = {cond['M']}) ===")
    print(f"  T∞={T:.3f} K  U∞={U:.2f} m/s  p∞={p:.3f} Pa  ρ∞={ro:.6e} kg/m³  μ∞={mu:.4e} Pa·s")
    print(f"  CPG EOS 整合: ρ=p/(RT)={ro_eos:.6e} (差 {d_ro*100:.3f} %), R={R_GAS:.3f} J/(kg·K)")
    if d_ro > 0.01:
        raise SystemExit("CPG EOS と derived.json の (p,T,ρ) が 1 % 以上ずれている")
    a_inf = float(np.sqrt(GAMMA * R_GAS * T))
    print(f"  a∞={a_inf:.3f} m/s → M = {U/a_inf:.3f} (台帳 {cond['M']})")

    # 自由流の乱流量 (Tu と μ_t/μ から)。0 は不可 (node SST は k=0 で初手発散する)。
    k_inf = 1.5 * (a.tu * U) ** 2
    om_inf = ro * k_inf / (a.mut_ratio * mu)
    print(f"  乱流入口: Tu={a.tu*100:.2f} %, μ_t/μ={a.mut_ratio:g} → "
          f"k∞={k_inf:.6g} m²/s², ω∞={om_inf:.6g} 1/s")

    mesh = CASE / "mesh" / f"{a.mesh}.h5"
    if not mesh.exists():
        raise SystemExit(f"メッシュが無い: {mesh} (あるのは "
                         f"{[q.stem for q in (CASE / 'mesh').glob('*.h5')]})")

    # --- run ディレクトリ (旧 res / ログは必ず消す: 混在事故の防止) ---
    rd = CASE / a.run
    rd.mkdir(exist_ok=True)
    for f in (list(rd.glob("res_*")) + list(rd.glob("*.log"))
              + list(rd.glob("residual_history*")) + list(rd.glob("CONVERGENCE_VERDICT*"))
              + list(rd.glob("RUN_PROVENANCE*"))):
        f.unlink()
    shutil.copy(mesh, rd / "mesh.h5")
    (rd / "probe.yaml").write_text(
        "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n", encoding="utf-8")

    bc = BC.format(ro=ro, u=U, p=p, tt=o["Tt"], kinf=k_inf, ominf=om_inf,
                   tw_plate=a.plate_Tw, tw_gap=a.gap_Tw)
    (rd / "bcondConfig.yaml").write_text(bc, encoding="utf-8")

    # --- 入口 BL 分布 (自前で再構成せず gen_inlet_csv.py を呼ぶ) ---
    subprocess.run([sys.executable, str(HERE / "gen_inlet_csv.py"), "--re", repr(a.re),
                    "--Tw", repr(a.plate_Tw), "--H", repr(H),
                    "--out", str(rd / "inlet_profile_1.csv")], check=True)

    (rd / "case_setup.json").write_text(json.dumps(
        dict(re_m=o["Re_m"], M=cond["M"], T_inf=T, U_inf=U, ro_inf=ro, p_inf=p, mu_inf=mu,
             Tt=o["Tt"], T_aw=o["T_aw"], a_inf=a_inf, k_inf=k_inf, omega_inf=om_inf,
             tu=a.tu, mut_ratio=a.mut_ratio, Tw_plate=a.plate_Tw, Tw_gap=a.gap_Tw,
             Tw_ratio=a.plate_Tw / a.gap_Tw,        # conditions.json の T_surf/T_gap
             mesh=a.mesh, domain_H=H,
             gap_W_mm=geo["gap"]["width"], gap_D_mm=geo["gap"]["depth"],
             gas=dict(model="CPG", cp=CP, gamma=GAMMA, R=R_GAS, prandtlLam=PRANDTL_LAM,
                      viscMethod="1 (Sutherland)"),
             eos_floor=dict(pMin=PMIN, roMin=ROMIN, tMin=TMIN),
             ic_delta=a.ic_delta, inlet_profile="inlet_profile_1.csv (tools/gen_inlet_csv.py)"),
        indent=2, ensure_ascii=False), encoding="utf-8")

    patch_ic(rd / "mesh.h5", st, a.plate_Tw, a.gap_Tw, k_inf, om_inf, a.ic_delta)

    # --- 段構成 ---
    # リミッタ基準値を**自由流で固定**する (auto だと開始場依存の作用素になり、分割実行が
    # 連続実行と一致しない。2026-09-20 codex result M10)。
    # thermCond は thermCondMethod: 1 (Prandtl 則) でも必須キー。基準値として mu*cp/Pr を入れる
    common = dict(mu_inf=mu, cp=CP, gam=GAMMA, pr=PRANDTL_LAM, kcond=mu * CP / PRANDTL_LAM,
                  pmin=PMIN, romin=ROMIN, tmin=TMIN,
                  relax=1.0, kinf=k_inf, ominf=om_inf, turb="sst",
                  ro_ref=ro, p_ref=p, a_ref=a_inf)
    stages = build_stages(a, common)
    # 収束判定の区間は段名でなく**実効設定**で決まる。本段だけ標準名 (residual_history.csv,
    # CONVERGENCE_VERDICT.txt) を残し、前段は _<tag> を付けて退避する。
    sm = StageManifest(rd)
    for tag, kw in stages:
        hist = "residual_history.csv" if tag == "main" else f"residual_history_{tag}.csv"
        sm.add(tag, CFG.format(**kw), bc, history=hist)
    sm.write()
    (rd / "stages.json").write_text(json.dumps(
        [{**kw, "tag": tag} for tag, kw in stages], indent=2), encoding="utf-8")
    print(f"段構成 ({len(stages)} 段): " + " -> ".join(
        f"{tag}(cfl {kw['cfl']:g}, {kw['nsteps']} step)" for tag, kw in stages))

    (rd / "solverConfig.yaml").write_text(CFG.format(**stages[0][1]), encoding="utf-8")
    if a.dry or a.prepare_only:
        if a.dry:
            for tag, kw in stages:
                print(f"\n===== stage {tag} : solverConfig.yaml =====\n" + CFG.format(**kw))
            print("===== bcondConfig.yaml =====\n" + bc)
        print(f"{'dry' if a.dry else 'prepare-only'}: {rd} を用意した (forge は起動しない)")
        return

    for tag, kw in stages:
        (rd / "solverConfig.yaml").write_text(CFG.format(**kw), encoding="utf-8")
        print(f"--- stage {tag}: cfl={kw['cfl']} conv={kw['conv']} turb={kw['turb']} "
              f"→ {kw['nsteps']} step")
        # **失敗段から古い場を引き継がない**。run_case.sh の戻り値と、当該段の最終出力・
        # NaN の有無を確認してから次段へ進む (case/56 の 2026-09-20 codex result-2 M4)。
        rc = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)]).returncode
        cur = rd / f"res_{kw['nsteps']}.h5"
        if rc != 0 or not cur.exists() or list(rd.glob("res_nan_*.h5")):
            raise SystemExit(f"stage {tag}: rc={rc}, res_{kw['nsteps']}.h5={cur.exists()}, "
                             f"nan={bool(list(rd.glob('res_nan_*.h5')))} "
                             "→ 段の失敗。古い場を引き継がずに中断する")
        if tag == stages[-1][0]:
            print(f"    本段 {tag}: res_*.h5 {len(list(rd.glob('res_[0-9]*.h5')))} 個と"
                  " 標準名の残差・VERDICT を保持 (準定常判定に使う)")
            break
        # **その段が書いた res を名前で取る** (2026-09-20 codex result M9)。全段の最大番号を
        # 取ると mid 3000 step が ramp 2500 step より大きく、ramp/本段が古い mid の場から
        # 再開してしまう。
        subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(cur),
                        str(rd / "mesh.h5")], check=True)
        seed_turb(rd / "mesh.h5", k_inf, om_inf)
        for suf in ("forge_run.log", "residual_history.csv", "residual_history.png",
                    "CONVERGENCE_VERDICT.txt", "RUN_PROVENANCE.txt"):
            src = rd / suf
            if src.exists():
                src.rename(rd / f"{Path(suf).stem}_{tag}{Path(suf).suffix}")
        for f in rd.glob("res_[0-9]*.h5"):      # 前段の res は消す (本段の時系列に混ぜない)
            f.unlink()
    print("done:", rd)


if __name__ == "__main__":
    main()
