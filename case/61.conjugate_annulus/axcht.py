#!/usr/bin/env python3
r"""軸対称 CHT 検証 (case/61 同心円環・case/62 同軸円板) の共通部品。

plan [`boundary-cht-axisymmetric-fem2d.md`](../../plans/accepted/boundary-cht-axisymmetric-fem2d.md) §6 V-ax2 / V-ax2b。
case/62 もこのファイルを import する (幾何だけが違い、手順と判定の作法は同じ)。

- 流体メッシュ: gmsh (.geo → msh4.1) → `convertGmshToForge` (node) → **静止・一様 IC をパッチ** (case/52 と同じ作法。
  変換器の `initial:` は `uniform_p101325_u10` のまま通し、`VALUE` を上書きする。`setInitial.hpp` は触らない)
- 固体帯: `case/58.conjugate_slot/gen_solid_strip.py` と**同じ帯トポロジ**
  `nid(j,i) = j*n_i + i` (j = 背面 0 → 界面 n_s、i = 界面に沿った並び)、各セルの対角は (j,i)–(j+1,i+1)。
  界面の並び i と背面の向きで対角の見かけの向きが決まる (case/61 は左上→右下、case/62 は左下→右上。
  各 case の gen_solid.py が登録文の向きを検査する)
- 評価: 流体の r 重み面積 $A^r_{\rm fluid}=|S_f|\max(r_f,r_{\rm floor})$ を `mesh.h5` の PLANES から組む
  (`conjugateWall.cpp` の `fillInterfaceDiagnostics` と同じ定義)、固体の $Q_{\rm sol}=(Ku-b)_{\rm iface}$ は
  `Fem2DOperator(axisym=True)` (case/58 `eval_v6p.py` と同じ作法)
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))

# 既定のビルド。変換器は幾何を書くだけなので float ビルドのもので良い (座標は h5 に倍精度で残らないが、
# 流体・固体・評価器が同じ値を読むので整合は保たれる)。FORGE_BUILD で差し替える。
BUILD = Path(os.environ.get("FORGE_BUILD", ROOT / "solver_density_cuda" / "build-ypls"))
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:"
           + os.environ.get("LD_LIBRARY_PATH", ""))

# ---- 共通の物性 (case/52 と同じ。登録値) ----
CP, GAMMA = 1004.5, 1.4
P0 = 1013.25            # 低圧で熱拡散率を上げる (case/52 README)。解析解は圧力に依らない
T_IC = 325.0            # 一様 IC (流体も固体も)。plan §6 V-ax2「流体残差」


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ============================================================ 流体メッシュ
def gmsh_and_convert(mdir: Path, name: str, geo: str, conv_cfg: str, conv_bc: str) -> Path:
    """.geo を書いて gmsh → msh4.1 → node 変換。戻り値は変換後の h5。"""
    mdir.mkdir(parents=True, exist_ok=True)
    work = mdir / f"_conv_{name}"
    work.mkdir(exist_ok=True)
    (mdir / f"{name}.geo").write_text(geo)
    subprocess.run(["gmsh", "-2", str(mdir / f"{name}.geo"), "-o", str(mdir / f"{name}.msh"),
                    "-format", "msh41", "-v", "1"], check=True, stdout=subprocess.DEVNULL)
    (work / "solverConfig.yaml").write_text(conv_cfg)
    (work / "bcondConfig.yaml").write_text(conv_bc)
    conv = BUILD / "convertGmshToForge"
    with open(work / "convert.log", "w") as log:
        subprocess.run([str(conv), str(mdir / f"{name}.msh"), "out.h5"], cwd=work, check=True,
                       env=ENV, stdout=log, stderr=subprocess.STDOUT)
    out = mdir / f"{name}.h5"
    (work / "out.h5").replace(out)
    return out


def patch_static_ic(h5path: Path, p0: float = P0, T0: float = T_IC) -> float:
    """静止・一様 (p0, T0) の IC を VALUE に書く (roe = ro cv T, ek = 0)。dtype は変換器の出力に揃える。"""
    R = CP * (GAMMA - 1.0) / GAMMA
    cv = CP / GAMMA
    ro = p0 / (R * T0)
    with h5py.File(h5path, "a") as f:
        g = f["VALUE"]
        n = g["ro"].shape[0]
        for name, val in (("ro", ro), ("roUx", 0.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", ro * cv * T0)):
            dt = g[name].dtype
            del g[name]
            g.create_dataset(name, data=np.full(n, val, dtype=dt))
        f["VALUE"].attrs["static_ic"] = f"p0={p0} T0={T0} u=0 (patched by axcht.patch_static_ic)"
    return ro


def mesh_quality(h5: Path, ar_max: float | None = None) -> str:
    cmd = [sys.executable, str(TOOLS / "check_mesh_quality.py"), str(h5)]
    if ar_max:
        cmd += ["--ar-max", str(ar_max)]
    p = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    out = p.stdout + p.stderr
    v = [l for l in out.splitlines() if "VERDICT" in l]
    return (v[-1].strip() if v else "(no VERDICT)") + "\n" + out


# ============================================================ 固体帯
def solid_strip(iface_xy: np.ndarray, normal_axis: int, back_coord: float, ns: int):
    r"""界面節点列 `iface_xy` (並び = 帯の i) から、法線座標 (`normal_axis` 0=x / 1=y) が `back_coord` の
    背面までを `ns` 層に一様分割した帯を作る (界面が法線座標一定の直線であること)。

    `gen_solid_strip.py` (case/58) と同じ帯トポロジ:
        節点 nid(j,i) = j*n_i + i、j = 0 (背面) … ns (界面)
        三角形 [(j,i),(j+1,i),(j+1,i+1)], [(j,i),(j+1,i+1),(j,i+1)]  (対角 = (j,i)–(j+1,i+1))
        界面 outer_edges = 行 j=ns、背面 hole1 = 行 j=0
    界面の行は壁ダンプの座標をそのまま使い (1 ulp も動かさない)、背面の法線座標は `back_coord` ちょうどにする。
    """
    iface_xy = np.asarray(iface_xy, float)
    ni = len(iface_xy)
    c0 = iface_xy[:, normal_axis]
    if np.ptp(c0) > 0.0:
        raise SystemExit(f"界面の法線座標が一定でない (ptp {np.ptp(c0):.3e} m)")
    nodes = np.repeat(iface_xy[None, :, :], ns + 1, axis=0).copy()
    for j in range(ns):                                   # j=ns (界面) は触らない
        nodes[j, :, normal_axis] = back_coord + (c0 - back_coord) * (j / ns)
    nodes = nodes.reshape(-1, 2)

    def nid(j, i):
        return j * ni + i

    tris = []
    for j in range(ns):
        for i in range(ni - 1):
            tris.append([nid(j, i), nid(j + 1, i), nid(j + 1, i + 1)])
            tris.append([nid(j, i), nid(j + 1, i + 1), nid(j, i + 1)])
    outer = np.array([[nid(ns, i), nid(ns, i + 1)] for i in range(ni - 1)], int)
    back = np.array([[nid(0, i), nid(0, i + 1)] for i in range(ni - 1)], int)
    return nodes, np.array(tris, int), outer, back


def diag_direction(nodes, ni, ns) -> str:
    """セル (j=0,i=0) の対角 (j,i)–(j+1,i+1) の見かけの向きを「左上→右下」「左下→右上」で返す。
    全セルで同じ対角なので 1 セルで決まる (solid_strip が保証する)。"""
    a, b = nodes[0 * ni + 0], nodes[1 * ni + 1]
    dx, dy = b[0] - a[0], b[1] - a[1]
    return "左上→右下" if dx * dy < 0 else "左下→右上"


def write_solid(out: Path, nodes, tris, outer, back, k_s, h, Tc, note, grid_meta: dict):
    """npz + json を書き、`solid_mesh_to_h5.py` で h5 にし、格子の記述 (`*.grid.json`) を書く。"""
    out.parent.mkdir(parents=True, exist_ok=True)
    npz = out.with_suffix(".npz")
    np.savez(npz, nodes=nodes, tris=tris, outer_edges=outer, hole1=back)
    spec = {"_note": note, "mesh_npz": npz.name, "k_solid": k_s,
            "holes": [{"h": h, "T_c": Tc}], "T_init": T_IC}
    js = out.with_suffix(".json")
    js.write_text(json.dumps(spec, ensure_ascii=False, indent=2))
    h5 = out.with_suffix(".h5")
    p = subprocess.run([sys.executable, str(TOOLS / "solid_mesh_to_h5.py"), "--solid", str(js),
                        "--npz", str(npz), "--out", str(h5)], capture_output=True, text=True, check=True)
    with h5py.File(h5, "r") as f:
        perm = np.asarray(f["MESH/PERM"][:], int)
        content = str(f.attrs["content_sha1"])
        iface_sha = str(f.attrs["iface_sha1"])
    inv = np.empty(len(perm), int)
    inv[perm] = np.arange(len(perm))
    meta = dict(grid_meta)
    meta.update({
        "solid_npz": npz.name, "solid_h5": h5.name,
        "npz_sha256": sha256(npz), "h5_content_sha1": content, "h5_iface_sha1": iface_sha,
        "n_nodes": int(len(nodes)), "n_tris": int(len(tris)),
        "node_order": "nid(j,i) = j*n_i + i (npz 順)。j = 0 背面 … n_s 界面、i = 界面に沿った並び。"
                      "h5 は RCM 並べ替え済み (h5 index = npz_to_h5[npz index])",
        "connectivity": "各セル [(j,i),(j+1,i),(j+1,i+1)] と [(j,i),(j+1,i+1),(j,i+1)] (gen_solid_strip.py と同じ)",
        "diagonal": diag_direction(nodes, grid_meta["n_iface"], grid_meta["n_layers"]),
        "npz_to_h5": inv.tolist(),
    })
    out.with_suffix(".grid.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1))
    return h5, meta, p.stdout


def read_wall_coords(dump: Path) -> np.ndarray:
    with h5py.File(dump, "r") as f:
        return np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)


# ============================================================ run の読み出し
def last_step(run: Path) -> int:
    st = [int(re.search(r"res_(\d+)\.h5$", p).group(1)) for p in glob.glob(str(run / "res_[0-9]*.h5"))]
    if not st:
        raise SystemExit(f"{run} に res_<step>.h5 が無い")
    return max(st)


def read_yaml(path: Path) -> dict:
    import yaml
    return yaml.safe_load(Path(path).read_text())


def axis_r_floor(run: Path) -> float:
    """`conjugateWall.cpp` の axisRFloorOf と同じ: mesh.axisRFloor > 0 ならその値、でなければ 1e-20。"""
    cfg = read_yaml(run / "solverConfig.yaml")
    m = cfg.get("mesh", {}) or {}
    if int(m.get("isAxisymmetric", 0)) != 1 or int(m.get("axisymMethod", 0)) != 0:
        raise SystemExit(f"{run}/solverConfig.yaml が mesh.isAxisymmetric: 1 / axisymMethod: 0 でない")
    v = float(m.get("axisRFloor", 0.0) or 0.0)
    return v if v > 0.0 else 1.0e-20


def fluid_wall_areas(run: Path, pid: int, rfloor: float):
    """壁 physID の各壁節点 (bc 順 = 壁ダンプの順) の平面面積 |S_f| と r 重み面積 |S_f| max(r_f, r_floor)。

    `fillInterfaceDiagnostics` と同じ定義 (面重心 r_f = PLANES/centCoords の y)。
    戻り値 (node_xyz[n,3], A_planar[n], A_r[n], r_face[n])。"""
    with h5py.File(run / "mesh.h5", "r") as f:
        ip = np.asarray(f[f"BCONDS/{pid}/iPlanes"][:], int)
        ic = np.asarray(f[f"BCONDS/{pid}/iCells"][:], int)
        sa = np.asarray(f["PLANES/surfArea"][:], float)
        pc = np.asarray(f["PLANES/centCoords"][:], float).reshape(-1, 3)
        # 節点座標は MESH/COORD (node 変換では CELLS/centCoords は双対重心で、実行時にノード座標へ置換される)
        cc = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)
    A_pl = sa[ip]
    rf = pc[ip, 1]
    A_r = A_pl * np.where(rf > rfloor, rf, rfloor)
    return cc[ic], A_pl, A_r, rf


def wall_dump(run: Path, name: str, pid: int, step: int) -> dict:
    fn = run / f"res_{name}_{pid}_{step}.h5"
    if not fn.exists():
        raise SystemExit(f"{fn} が無い (outputHDFflg: 1 と interfaceDiag: 1 の run か)")
    with h5py.File(fn, "r") as f:
        V = f["VALUE"]
        d = {"xyz": np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)}
        for k in ("Ts", "iface_q_eff", "iface_q_eff_raw", "iface_Qf_eff"):
            if k not in V:
                raise SystemExit(f"{fn} に VALUE/{k} が無い")
            d[k] = np.asarray(V[k][:], float)
        for k in ("ifaceRraw", "ifaceFw"):            # 恒等式 q_raw A = R_raw − F_w の照合用 (無ければ省く)
            if k in V:
                d[k] = np.asarray(V[k][:], float)
        d["dtype"] = str(V["iface_Qf_eff"].dtype)
    return d


def check_wall_order(dump_xyz, mesh_xyz, what: str, tol: float = 1e-9):
    if dump_xyz.shape != mesh_xyz.shape or np.abs(dump_xyz - mesh_xyz).max() > tol:
        raise SystemExit(f"{what}: 壁ダンプの節点順が mesh.h5 の BCONDS iCells 順と一致しない")


def solid_state(run: Path, pid: int, step: int):
    """固体ダンプと run の solid.h5 から Q_sol = (K u − b)_iface [W/rad] (軸対称 Fem2DOperator)。

    戻り値 dict: T (全節点), q_hole, q_iface, iface (h5 index), A_r (固体側集中量 m²/rad),
    Qsol_if, iface_xy。"""
    from solid_fem2d import Fem2DOperator
    fn = run / f"res_solid_{pid}_{step}.h5"
    if not fn.exists():
        raise SystemExit(f"{fn} が無い")
    with h5py.File(fn, "r") as g:
        sc = np.asarray(g["MESH/COORD"][:], float).reshape(-1, 3)
        T = np.asarray(g["VALUE/T"][:], float)
        q_hole = np.asarray(g["VALUE/q_hole"][:], float)
        q_iface = np.asarray(g["VALUE/q_iface"][:], float)
    with h5py.File(run / "solid.h5", "r") as g:
        sxy = np.asarray(g["MESH/COORD"][:], float)
        ifn = np.asarray(g["IFACE/NODES"][:], int)
        robin = [(int(e[0]), int(e[1]), float(h), float(t)) for e, h, t in
                 zip(g["ROBIN/EDGES"][:], g["ROBIN/H"][:], g["ROBIN/TC"][:])]
        kT, kV = np.asarray(g["SOLID/K_T"][:], float), np.asarray(g["SOLID/K_V"][:], float)
        op = Fem2DOperator(sxy, np.asarray(g["MESH/TRIS"][:], int), ifn,
                           np.asarray(g["IFACE/EDGES"][:], int), robin,
                           float(kV[0]) if len(kV) == 1 else (kT, kV), axisym=True)
    if np.abs(sxy - sc[:, :2]).max() > 1e-12:
        raise SystemExit("固体ダンプの節点順が solid.h5 と一致しない (Q_sol を組めない)")
    K, bs = op.assemble_full(T)
    Qsol = (K @ T - bs)
    return {"T": T, "q_hole": q_hole, "q_iface": q_iface, "iface": ifn, "A_r": op.area,
            "Qsol_if": Qsol[ifn], "iface_xy": sxy[ifn], "xy": sxy, "op": op}


def match_iface_to_wall(iface_xy, wall_xyz, tol=1e-7):
    """固体界面節点 → 壁節点 (座標一致 1 対 1)。ソルバの initSolidFem2d と同じ 1e-7 m。"""
    j = np.array([int(np.argmin(np.hypot(wall_xyz[:, 0] - x, wall_xyz[:, 1] - y))) for x, y in iface_xy])
    d = np.hypot(wall_xyz[j, 0] - iface_xy[:, 0], wall_xyz[j, 1] - iface_xy[:, 1])
    if d.max() > tol or len(set(j.tolist())) != len(j):
        raise SystemExit(f"界面節点と壁節点が 1 対 1 に対応しない (最大ずれ {d.max():.3e} m)")
    return j


# ============================================================ 判定の部品
def gate(rows, name, value, tol, unit="", fmt="{:.6g}"):
    ok = bool(np.isfinite(value) and value <= tol)
    rows.append((name, value, tol, unit, ok, fmt))
    return ok


def print_rows(rows):
    print(f"\n{'項目':<56}{'値':>14}{'許容':>12}  判定")
    for name, v, tol, unit, ok, fmt in rows:
        print(f"{name:<56}{fmt.format(v):>14}{fmt.format(tol):>12} {unit:<6}{'PASS' if ok else '**FAIL**'}")


def q_eff_check(q, q_star, frac=5e-3):
    """全節点 max |q − q*| と許容 (frac × q*)。戻り値 (max_err, tol, ok)。"""
    q = np.asarray(q, float)
    err = float(np.max(np.abs(q - q_star))) if np.all(np.isfinite(q)) else float("inf")
    tol = frac * q_star
    return err, tol, err <= tol


def write_series_csv(run: Path, pid: int, Tc: float, out: Path):
    """節点ログ (`conjugate.node_log: 1`) から全界面節点の毎更新系列を書く (準定常の判定用)。

    列: step, Tw_mTc_<i> (= T_w − T_c)、q_<i> (= Q_f,i / A_i^r、固体側集中量で割った熱流束 [W/m²])。
    **絶対温度は書かない** (plan §6 V-ax2 準定常)。戻り値は列名 (step を除く)。"""
    from check_cht_interface import load_node_log
    nd, lg, _, ni = load_node_log(run, pid)
    A = nd[:, 4]
    ups = np.unique(lg[:, 0].astype(int))
    steps = np.empty(len(ups), int)
    Tw = np.empty((len(ups), ni))
    Q = np.empty_like(Tw)
    for k, u in enumerate(ups):
        blk = lg[lg[:, 0].astype(int) == u]
        idx = blk[:, 2].astype(int)
        Tw[k, idx] = blk[:, 6]
        Q[k, idx] = blk[:, 4] / A[idx]
        steps[k] = int(blk[0, 1])
    cols = [f"Tw_mTc_{i}" for i in range(ni)] + [f"q_{i}" for i in range(ni)]
    data = np.column_stack([steps, Tw - Tc, Q])
    np.savetxt(out, data, delimiter=",", header=",".join(["step"] + cols), comments="",
               fmt=["%d"] + ["%.12e"] * (data.shape[1] - 1))
    return cols, len(ups)


# ============================================================ run ディレクトリの組立て
def make_run(template: Path, run: Path, mesh_h5: Path, solid_h5: Path | None = None,
             dry: bool = False, nstep: int | None = None, note: str = ""):
    """template の config と入力 h5 を新しい run ディレクトリに置く (既存なら止める)。

    dry=True は `solverConfig_dry.yaml` / `bcondConfig_dry.yaml` を使う (共役なし 1 step)。
    nstep はローカルの起動確認用 (nStepOuter だけ差し替え、run 名で区別すること)。"""
    import shutil
    if run.exists():
        raise SystemExit(f"{run} は既にある (run を使い回さない: AGENTS.md)")
    run.mkdir(parents=True)
    sfx = "_dry" if dry else ""
    cfg = (template / f"solverConfig{sfx}.yaml").read_text()
    if nstep is not None:
        cfg, n = re.subn(r"last: \{nStepOuter: \d+\}", f"last: {{nStepOuter: {nstep}}}", cfg)
        if n != 1:
            raise SystemExit("nStepOuter を差し替えられない (template の書式が変わった)")
    (run / "solverConfig.yaml").write_text(cfg)
    shutil.copy(template / f"bcondConfig{sfx}.yaml", run / "bcondConfig.yaml")
    shutil.copy(template / "probe.yaml", run / "probe.yaml")          # 無いとソルバが止まる (プローブ点は無し)
    shutil.copy(mesh_h5, run / "mesh.h5")
    prov = [f"mesh.h5  <- {mesh_h5} (sha256 {sha256(mesh_h5)})"]
    if not dry:
        if solid_h5 is None:
            raise SystemExit("共役 run には --solid が要る")
        shutil.copy(solid_h5, run / "solid.h5")
        g = Path(str(solid_h5)[:-3] + ".grid.json")
        if g.exists():
            shutil.copy(g, run / "solid.grid.json")
        prov.append(f"solid.h5 <- {solid_h5} (sha256 {sha256(solid_h5)})")
    prov.append(f"template: {template}{' (dry)' if dry else ''}" + (f", nStepOuter={nstep}" if nstep else ""))
    if note:
        prov.append(f"note: {note}")
    (run / "RUN_INPUTS.txt").write_text("\n".join(prov) + "\n")
    print(f"[make_run] {run}\n  " + "\n  ".join(prov))


# ============================================================ 評価器の自己試験 (流体計算なし)
def synth_run(dst: Path, mesh_h5: Path, solid_h5: Path, solver_cfg: Path, bcond_cfg: Path,
              walls: list, pid_cj: int, name_cj: str, q_cj: float, step: int = 1):
    r"""流体計算をせずに、評価器が読むファイル一式を**合成**した run を作る (評価器の正例・負例用)。

    - 共役壁: $Q_{f,i}=q_{cj}A^r_{{\rm fluid},i}$ (一様熱流束の厳密な面積分)、`iface_q_eff`=`raw`=$q_{cj}$、
      固体は `Fem2DOperator(axisym=True)` で $Ku=b+E^{\rm T}Q_f$ を直接解き、`Ts` = 固体の界面温度
    - `walls` の各 (pid, name, q) は非連成の壁: $Q_{f,i}=qA^r_{{\rm fluid},i}$、`Ts` はダミー
    - `q_hole` は組立てと同じ consistent 行列 $M_e(T-T_c)$ (plan §4.4b)
    出力は float64 (FP64 run と同じ)。戻り値は固体の界面温度 (壁節点順)。"""
    import shutil

    import scipy.sparse.linalg as spla
    from solid_fem2d import Fem2DOperator
    dst.mkdir(parents=True)
    shutil.copy(mesh_h5, dst / "mesh.h5")
    shutil.copy(solid_h5, dst / "solid.h5")
    shutil.copy(solver_cfg, dst / "solverConfig.yaml")
    shutil.copy(bcond_cfg, dst / "bcondConfig.yaml")
    rfloor = axis_r_floor(dst)
    (dst / f"res_{step}.h5").write_bytes(b"")          # last_step 用の目印 (中身は読まない)

    def dump(name, pid, Ts, q, Q):
        xyz, _, _, _ = fluid_wall_areas(dst, pid, rfloor)
        with h5py.File(dst / f"res_{name}_{pid}_{step}.h5", "w") as f:
            f.create_dataset("MESH/COORD", data=xyz.reshape(-1))
            # 強形式の定義 q_raw A = R_raw − F_w を満たすよう R_raw = Q・F_w = 0 と置く (蓄積項 0 = 定常)
            for k, v in (("Ts", Ts), ("iface_q_eff", q), ("iface_q_eff_raw", q), ("iface_Qf_eff", Q),
                         ("ifaceRraw", Q), ("ifaceFw", np.zeros_like(np.asarray(Q, float)))):
                f.create_dataset(f"VALUE/{k}", data=np.asarray(v, float))

    for pid, name, q in walls:
        xyz, _, Ar, _ = fluid_wall_areas(dst, pid, rfloor)
        dump(name, pid, np.full(len(Ar), np.nan), np.full(len(Ar), q), q * Ar)

    xyz, _, Ar, _ = fluid_wall_areas(dst, pid_cj, rfloor)
    Qf = q_cj * Ar
    with h5py.File(solid_h5, "r") as g:
        sxy = np.asarray(g["MESH/COORD"][:], float)
        ifn = np.asarray(g["IFACE/NODES"][:], int)
        robin = [(int(e[0]), int(e[1]), float(h), float(t)) for e, h, t in
                 zip(g["ROBIN/EDGES"][:], g["ROBIN/H"][:], g["ROBIN/TC"][:])]
        kV = float(np.asarray(g["SOLID/K_V"][:], float)[0])
        op = Fem2DOperator(sxy, np.asarray(g["MESH/TRIS"][:], int), ifn,
                           np.asarray(g["IFACE/EDGES"][:], int), robin, kV, axisym=True)
    jw = match_iface_to_wall(sxy[ifn], xyz)
    K, b = op.assemble_full(np.full(len(sxy), 300.0))
    rhs = b.copy()
    rhs[ifn] += Qf[jw]
    T = spla.spsolve(K.tocsc(), rhs)
    qh = np.zeros(len(sxy))
    for (n0, n1, h, Tc) in robin:
        L = float(np.linalg.norm(sxy[n1] - sxy[n0]))
        Me, _ = op._robin_edge(n0, n1, L)
        dT = np.array([T[n0] - Tc, T[n1] - Tc])
        qh[[n0, n1]] += h * (Me @ dT)
    qi = np.zeros(len(sxy))
    qi[ifn] = Qf[jw]
    with h5py.File(dst / f"res_solid_{pid_cj}_{step}.h5", "w") as f:
        f.create_dataset("MESH/COORD", data=np.column_stack([sxy, np.zeros(len(sxy))]))
        for k, v in (("T", T), ("q_hole", qh), ("q_iface", qi), ("k_s", np.full(len(sxy), kV))):
            f.create_dataset(f"VALUE/{k}", data=v)
    Tw = np.empty(len(Ar))
    Tw[jw] = T[ifn]
    dump(name_cj, pid_cj, Tw, np.full(len(Ar), q_cj), Qf)
    return Tw


def synthetic_ab(mesh_h5: Path, pid: int, solid_h5: Path, q_star: float, rfloor: float = 1.0e-20):
    r"""判別 A/B (plan §6 V-ax2b、流体計算不要): 合成荷重 $Q_i=q_*A_i^r$ を
    A = $A_i^r$ で割る / B = $A_{{\rm planar},i}$ で割る。A PASS・B FAIL でなければ評価器か面積定義に不備。

    $A_i^r$ は固体側の集中量 (Fem2DOperator axisym)、$A_{{\rm planar},i}$ は同じ界面の平面の集中辺長。
    流体側 (mesh.h5 PLANES の $|S_f|\max(r_f,r_{\rm floor})$ と $|S_f|$) でも同じ A/B を取る。"""
    from solid_fem2d import Fem2DOperator
    with h5py.File(solid_h5, "r") as g:
        sxy = np.asarray(g["MESH/COORD"][:], float)
        tris = np.asarray(g["MESH/TRIS"][:], int)
        ifn = np.asarray(g["IFACE/NODES"][:], int)
        ife = np.asarray(g["IFACE/EDGES"][:], int)
    Ar_s = Fem2DOperator(sxy, tris, ifn, ife, [], 1.0, axisym=True).area
    Apl_s = Fem2DOperator(sxy, tris, ifn, ife, [], 1.0, axisym=False).area
    with h5py.File(mesh_h5, "r") as f:
        ip = np.asarray(f[f"BCONDS/{pid}/iPlanes"][:], int)
        Apl_f = np.asarray(f["PLANES/surfArea"][:], float)[ip]
        rf = np.asarray(f["PLANES/centCoords"][:], float).reshape(-1, 3)[ip, 1]
    Ar_f = Apl_f * np.where(rf > rfloor, rf, rfloor)
    out = {}
    for side, Ar, Apl in (("solid", Ar_s, Apl_s), ("fluid", Ar_f, Apl_f)):
        Q = q_star * Ar
        eA, tol, okA = q_eff_check(Q / Ar, q_star)
        eB, _, okB = q_eff_check(Q / Apl, q_star)
        out[side] = {"A_err": eA, "B_err": eB, "tol": tol, "A_pass": okA, "B_pass": okB}
    return out


def identity_check(run: Path, pid: int, name: str, step: int, rfloor: float, use_floor: bool = True):
    r"""恒等式 (plan §6 V-ax2b): `iface_q_eff`·$A^r_{\rm fluid}$ = `iface_Qf_eff`、
    `iface_q_eff_raw`·$A^r_{\rm fluid}$ = $R^{\rm raw}-F_w$ (`fillInterfaceDiagnostics` の定義。強形式の等温壁) を全節点で照合する。

    $A^r_{\rm fluid}=|S_f|\max(r_f, r_{\rm floor})$。`use_floor=False` は床を外した面積 (|S_f| r_f) で照合する
    (軸に触れる負例で「床を使っていること」を示すための対照)。戻り値 dict。"""
    cfg = read_yaml(run / "solverConfig.yaml")
    weak = int((cfg.get("mesh", {}) or {}).get("nodeIsothermalEnergyBC", 0) or 0) != 0
    xyz, Apl, Ar, rf = fluid_wall_areas(run, pid, rfloor)
    d = wall_dump(run, name, pid, step)
    check_wall_order(d["xyz"], xyz, name)
    A = Ar if use_floor else Apl * rf
    Q = d["iface_Qf_eff"]

    def rel(lhs, rhs):
        den = np.abs(rhs)
        return np.where(den > 0, np.abs(lhs - rhs) / np.where(den > 0, den, 1.0),
                        np.where(lhs == rhs, 0.0, np.inf))

    r_eff = rel(d["iface_q_eff"] * A, Q)
    out = {"dtype": d["dtype"], "n": len(Q), "rel_q_eff": float(np.max(r_eff)),
           "n_floor": int(np.sum(~(rf > rfloor))), "r_face_min": float(rf.min()), "rfloor": rfloor}
    if "ifaceFw" in d and (weak or "ifaceRraw" in d):
        Rraw = np.zeros_like(Q) if weak else d["ifaceRraw"]
        out["rel_q_raw"] = float(np.max(rel(d["iface_q_eff_raw"] * A, Rraw - d["ifaceFw"])))
    else:
        out["rel_q_raw"] = None
    return out
