#!/usr/bin/env python3
r"""C3X 翼列 run の準備と段階起動 (V5 段 (a): **実測壁温を課して $h$ を当てる**)。

段階起動は procedures/divergence-and-startup.md の作法:
  warm (層流・1 次・低 CFL) → sst_soft (1 次) → sst_mid (1 次・CFL 上げ) → ramp (2 次) → main (2 次・本段)
各段は同じ run ディレクトリで forge を順に呼び、**前段の最終場を interp_field で引き継ぐ** (同一メッシュ)。
段ごとに `forge_run_<stage>.log` と `CONVERGENCE_VERDICT_<stage>.txt` を残す。

壁は `wall_isothermal` + `ints: {wallProfile: 1}` で、**実測壁温を弧長で写した CSV** を与える
(`make_wall_profile` で作る。マッピングは infer_internal_bc.py と同じ規約)。

usage: python3 case/53.c3x_vane_cht/tools/setup_run.py <run_dir> [--run run108] [--stages ...] [--go]
"""
import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
from run_data import TABLES                                    # noqa: E402
from gen_solid_mesh import parse_msh41                         # noqa: E402

H0, TREF = 1135.0, 811.0
# M1 は**報告値 (0.17) をそのまま入口に課すと閉塞流量を 4 % 超える** (実測 2026-09-20:
# 入口 Pt が 470 kPa まで上がり、出口が M=1.87 になって亜音速出口 BC と矛盾し発散した)。
# 入口は質量流束指定なので、**報告の出口条件 (M2, Pt, Tt, スロート面積) と整合する流量**から決める:
#   mdot = rho2 U2 A_throat,  rho1 U1 = mdot / pitch  ->  M1 = 0.1628 (報告の 0.17 に対し -4 %)
# 得られた**入口全圧が報告の 319.5 kPa になるか**を後段で照合する (これが整合性の検査)。
COND = {
    "run108": dict(vane="c3x", Pt=319500.0, Tt=786.0, M2=0.90, Tu=0.065, M1=0.1628,
                   M1_reported=0.17, throat_m=0.03292, pitch_m=0.11773),
    # Mark II (超音速出口)。報告 表 VIII: Pt1 48.89 psia, Tt1 788 K, M1 0.19, M2 1.04
    "run42": dict(vane="markii", Pt=337100.0, Tt=788.0, M2=1.04, Tu=0.065, M1=0.19,
                  M1_reported=0.19, throat_m=0.03983, pitch_m=0.12974),
}
CASE = {"c3x": "case/53.c3x_vane_cht", "markii": "case/54.markii_vane_cht"}

FORGE = ROOT / "solver_density_cuda/.build-native/relwithdebinfo/forge"
TOOLS = ROOT / "solver_density_cuda/tools"
ENV_LD = "/usr/lib/x86_64-linux-gnu/hdf5/serial"


def solver_cfg(stage, nstep, out_int, mesh="mesh.h5", precond=0):
    """段ごとの solverConfig。**1 次/2 次・乱流・CFL だけを変える** (区間判定のため)。"""
    model = '"none"' if stage == "warm" else '"sst"'
    # **convMethod 1 は (リミッタ無しでも) 2 次**。起動段は 0 = 1 次にする
    # (codex 診断 2026-09-20: 以前の「1 次で起動」は実際には無制限 2 次だった)。
    conv, lim = (1, 2) if stage.startswith(("ramp", "main")) else (0, 0)
    # main2 / main4 は本段 (2 次) で cfl だけ上げた区間 (区間判定は stage_manifest でなく
    # 「数値設定が同一の区間」なので、cfl だけの違いは連結してよい)。
    cfl = {"warm": 0.05, "sst_soft": 0.2, "sst_mid": 0.5, "ramp": 0.5,
           "main": 0.5, "main2": 2.0, "main4": 4.0}[stage]
    # **前縁よどみ点の k 過大生成**で SST が step 78 で落ちたので、SST 段は Kato-Launder を入れる
    # (実測 2026-09-20: KL 無しは cfl 0.02 まで、KL 有りで 0.05。KL 有りでも 0.2 は落ちる)。
    kl = 0 if stage == "warm" else 1
    return f"""# C3X cascade — stage {stage}
mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "{mesh}", valueFileName: "{mesh}"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1, visc: 3.3e-5, thermCond: 0.05, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nstep}}}
  deltaT: {{control: 1, dt: 1e-9, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: 1.0, blockDPLUR: 1, lowMachPrecond: {precond}, dt_min: 1e-12, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 4
space: {{convMethod: {conv}, limiter: {lim}, pRef: 250000.0}}
turbulence: {{model: {model}, scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: {kl}, wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInf: 57.4, omegaInf: 240000.0}}
output: {{level: 1, interfaceDiag: 1}}
initial: "uniform_p101325_u10"
"""


def bcond_cfg(pitch_m, Ps_exit, Pt, Tt, M1=0.17, wall_profile=True, inlet="pt"):
    wp = "{wallProfile: 1}" if wall_profile else " "
    # **入口は全圧入口 `inlet_Pressure_dir` を既定にした (2026-09-20)**。
    # 一様速度入口 (`inlet_uniformVelocity`) は「速度 3 成分 + config エントロピー」を課す
    # 正しく posed な BC だが、**全温を固定しない**: 実測 (run_0001_measTw) では内部静圧が
    # config アンカー (313.6 kPa) でなく 252 kPa に落ち着き、入口 Tt が 786 → 738 K と
    # 48 K 低くなった。h の比較には Tt が効くので使えない。
    # `inlet_Pressure_dir` が step 214 で NaN 化していた真因は BC カーネルの 2 か所の欠陥で、
    # 2026-09-20 に修正済 (boundaryCond_d.cu): (1) P_c>Pt で sqrt の引数が負 → NaN、
    # (2) 方向ベクトルを保持する bvar Ux/Uy/Uz をカーネルが次元付き速度で上書きするため
    #     M→0 で長さ 0 になり次 step が 0/0。
    ro1 = Pt / (1.0 + 0.2 * M1 ** 2) ** 3.5 / (287.0 * Tt / (1.0 + 0.2 * M1 ** 2))
    T1 = Tt / (1.0 + 0.2 * M1 ** 2)
    U1 = M1 * math.sqrt(1.4 * 287.0 * T1)
    Ps1 = Pt / (1.0 + 0.2 * M1 ** 2) ** 3.5
    if inlet == "pt":
        inlet_line = (f"inlet:   {{physID: 1, kind: inlet_Pressure_dir, outputHDFflg: 0, ints: , "
                      f"floats: {{Ux: 1.0, Uy: 0.0, Uz: 0.0, Pt: {Pt}, Tt: {Tt}, k: 57.4, omega: 240000.0}}}}")
    else:
        inlet_line = (f"inlet:   {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , "
                      f"floats: {{ro: {ro1:.6f}, Ux: {U1:.3f}, Uy: 0.0, Uz: 0.0, Ps: {Ps1:.1f}, "
                      f"k: 57.4, omega: 240000.0}}}}")
    return f"""{inlet_line}
outlet:  {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {Ps_exit:.1f}, Pt: {Pt}, Tt: {Tt}}}}}
per_low: {{physID: 3, kind: periodic, outputHDFflg: 0, ints: {{type: 0, partnerBCID: 4}}, floats: {{dx: 0.0, dy: {pitch_m:.9f}, dz: 0.0}}}}
per_up:  {{physID: 4, kind: periodic, outputHDFflg: 0, ints: {{type: 0, partnerBCID: 3}}, floats: {{dx: 0.0, dy: {-pitch_m:.9f}, dz: 0.0}}}}
wall:    {{physID: 5, kind: wall_isothermal, outputHDFflg: 1, ints: {wp}, floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 560.0}}}}
"""


def wall_profile_csv(vane, run, msh_path):
    """壁節点に実測壁温を写した CSV (x y z Ts) を作る。マッピングは infer_internal_bc.py と同規約。"""
    xy, tris, phys = parse_msh41(msh_path)
    wall_edges = phys[5]
    adj = {}
    for a_, b_ in wall_edges:
        adj.setdefault(int(a_), []).append(int(b_)); adj.setdefault(int(b_), []).append(int(a_))
    start = min(adj); loop = [start]; prev, cur = None, start
    while True:
        nxt = [v for v in adj[cur] if v != prev]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        if cur == start:
            break
        loop.append(cur)
    loop = np.array(loop, int)
    P = xy[loop]
    seg = np.hypot(*np.diff(np.vstack([P, P[:1]]), axis=0).T)
    n = len(loop)
    i_le = int(np.argmin(P[:, 0])); i_te = int(np.argmax(P[:, 0]))
    s = np.zeros(n); acc = 0.0
    for k in range(n):
        idx = (i_le + k) % n
        s[idx] = acc; acc += seg[idx]
    total = acc; s_te = s[i_te]
    side_fwd = s <= s_te
    arc_f, arc_b = s_te, total - s_te
    ss_is_fwd = arc_f > arc_b
    s_norm = np.where(side_fwd, s / arc_f, (total - s) / arc_b)
    is_ss = side_fwd if ss_is_fwd else ~side_fwd

    # **判読不能セル (None) は落とす**。壁温は弧長で内挿するので、欠測点を飛ばせばよい。
    rows = [r for r in TABLES[run]["rows"] if r[2] is not None]
    sd = np.array([r[0] for r in rows]); td = np.array([r[2] for r in rows]) * TREF
    i_stag = int(np.argmin(sd))
    dat = {"PS": (sd[:i_stag + 1][::-1], td[:i_stag + 1][::-1]), "SS": (sd[i_stag:], td[i_stag:])}
    Tw = np.array([np.interp(s_norm[k], *dat["SS" if is_ss[k] else "PS"]) for k in range(n)])
    return P, Tw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--vane", default="c3x")
    ap.add_argument("--stages", nargs="+",
                    default=["warm:2000", "sst_soft:2000", "sst_mid:3000", "ramp:3000", "main:20000"])
    ap.add_argument("--ic-p", type=float, default=250000.0)
    ap.add_argument("--ic-T", type=float, default=760.0)
    ap.add_argument("--ic-u", type=float, default=120.0)
    ap.add_argument("--resume", action="store_true", help="既存 run の mesh.h5 (前段の場) から続ける")
    ap.add_argument("--inlet", default="pt", choices=["pt", "u"], help="pt=全圧入口 / u=一様速度入口")
    ap.add_argument("--ps-exit", type=float, default=None,
                    help="出口静圧 [Pa] を直接指定 (既定は Pt/(1+0.2 M2^2)^3.5 = 全圧損失ゼロ仮定)")
    ap.add_argument("--precond", type=int, default=0, help="lowMachPrecond (低マッハ市松対策は 2)")
    ap.add_argument("--go", action="store_true", help="準備だけでなく実行もする")
    a = ap.parse_args()

    case = ROOT / CASE[a.vane]
    mesh_dir = case / "mesh"
    rd = Path(a.run_dir if Path(a.run_dir).is_absolute() else case / a.run_dir)
    rd.mkdir(parents=True, exist_ok=True)
    if not a.resume:
        shutil.copy(mesh_dir / f"fluid_{a.vane}.h5", rd / "mesh.h5")
    if a.resume:
        print("[setup] resume: 既存 mesh.h5 の場をそのまま使う")
    # **一様 IC を BC と整合させる** (既定の uniform_p101325_u10 は入口全圧 320 kPa と桁が合わず step 10 で発散した)。
    # 通路内の代表状態 (入口と出口の中間) を直接書き込む。
    import h5py
    p0, T0, u0 = a.ic_p, a.ic_T, a.ic_u
    if a.resume:
        p0 = None
    R = 1004.5 * (1.4 - 1.0) / 1.4
    ro = p0 / (R * T0) if p0 is not None else None
    with h5py.File(rd / "mesh.h5", "a") as f:
        n = f["CELLS/centCoords"].shape[0] // 3
        g = f.require_group("VALUE")
        if p0 is None:
            pass
        # **roK / roOmega も必ず入れる**: 変換直後の値は 0 で、SST を掛けると step 3 で
        # omega=0 → NaN になる (実測 2026-09-20)。自由流値をそのまま入れる。
        if p0 is not None:
            for name, val in (("ro", ro), ("roUx", ro * u0), ("roUy", 0.0), ("roUz", 0.0),
                              ("roe", ro * (1004.5 / 1.4 * T0 + 0.5 * u0 ** 2)),
                              ("roK", ro * 57.4), ("roOmega", ro * 240000.0)):
                if name in g:
                    del g[name]
                g.create_dataset(name, data=np.full(n, val, dtype=np.float32))
    if p0 is not None:
        print(f"[setup] IC: p={p0/1000:.1f} kPa, T={T0:.0f} K, u={u0:.0f} m/s (ro={ro:.4f})")
    (rd / "probe.yaml").write_text("outStepInterval: 100000\noutStepStart: 0\npoints:\nsurfaces:\n")

    c = COND[a.run]
    # 既定は損失ゼロ仮定。実際は翼列損失 (~2 %) の分だけ M2 が下がるので、
    # 1 回目の結果の出口 Pt を使って --ps-exit で追い込む。
    Ps_exit = a.ps_exit if a.ps_exit else c["Pt"] / (1 + 0.2 * c["M2"] ** 2) ** 3.5
    pitch_m = c["pitch_m"]          # COND は m で持つ (cm ではない)
    (rd / "bcondConfig.yaml").write_text(bcond_cfg(pitch_m, Ps_exit, c["Pt"], c["Tt"], c["M1"], inlet=a.inlet))

    P, Tw = wall_profile_csv(a.vane, a.run, mesh_dir / f"fluid_{a.vane}.msh")
    with open(rd / "wall_profile_5.csv", "w") as f:
        f.write("x y z Ts\n")
        for (x, y), t in zip(P, Tw):
            f.write(f"{x:.9e} {y:.9e} 0.0 {t:.6f}\n")
    print(f"[setup] wall profile: {len(P)} nodes, Tw {Tw.min():.1f} .. {Tw.max():.1f} K "
          f"(exit Ps {Ps_exit/1000:.1f} kPa)")

    (rd / "STAGES").write_text("\n".join(a.stages) + "\n")
    if not a.go:
        print(f"[setup] prepared {rd.relative_to(ROOT)} (実行は --go)")
        return

    import os
    env = dict(os.environ, LD_LIBRARY_PATH=ENV_LD + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
    for si, spec in enumerate(a.stages):
        # 段の指定は `name:steps` または `name:steps:ps_exit`。
        # **閉塞する翼列 (Mark II は M2=1.04) は背圧を段階的に下げる** — 一様 IC からいきなり
        # 超音速出口の背圧を課すと、出口ブロックが過膨張して圧力床に張り付き、
        # 「完走したが M=628・P=1 Pa」という壊れた場になる (実測 2026-09-20)。
        parts = spec.split(":")
        stage, nstep = parts[0], int(parts[1])
        if len(parts) > 2:
            ps_stage = float(parts[2])
            (rd / "bcondConfig.yaml").write_text(
                bcond_cfg(pitch_m, ps_stage, c["Pt"], c["Tt"], c["M1"], inlet=a.inlet))
            print(f"[setup]   back pressure -> {ps_stage/1000:.1f} kPa")
        (rd / "solverConfig.yaml").write_text(solver_cfg(stage, nstep, max(nstep // 4, 500), precond=a.precond))
        print(f"[setup] stage {stage} ({nstep} steps) ...", flush=True)
        with open(rd / f"forge_run_{stage}.log", "w") as log:
            r = subprocess.run([str(FORGE)], cwd=rd, stdout=log, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            sys.exit(f"[setup] stage {stage} failed (exit {r.returncode}); see forge_run_{stage}.log")
        subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(rd)],
                       stdout=open(rd / f"CONVERGENCE_VERDICT_{stage}.txt", "w"),
                       stderr=subprocess.STDOUT, env=env)
        res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
        if si + 1 < len(a.stages) and res:
            subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]),
                            str(rd / "mesh.h5")], check=True, env=env, stdout=subprocess.DEVNULL)
            for p in rd.glob("res_*"):
                p.rename(rd / f"_{stage}_{p.name}")
    print(f"[setup] done -> {rd.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
