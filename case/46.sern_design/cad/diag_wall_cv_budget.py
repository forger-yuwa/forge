#!/usr/bin/env python3
"""接続模型の起動診断 (plan tooling-sern-mesh-blocking §5.1 B1d ②)。

指定した CV (節点) について、`res_*.h5` の場から **1 次 SLAU の面別質量流束**を再計算し、
移流項・圧力差項・$\\chi$・両側音速・面積を面ごとに出す。あわせて CV の密度履歴と
「面流束の和 × dt / V」と実際の Δρ を突き合わせる (床補正が入っていれば差が出る。codex M4)。

実装は `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh` の
mdot に合わせる (`convMethod: 0` = 再構成なしの 1 次、L = owner の節点値, R = neighbor):

    Vn_p = U_L·n,  Vn_m = U_R·n,  c_hat = (c_L + c_R)/2
    g    = -max(min(M_p,0),-1) * min(max(M_m,0),1)
    Vn_hat_abs = (ro_L|Vn_p| + ro_R|Vn_m|)/(ro_L+ro_R)
    M_hat = min(1, sqrt((|U_L|^2+|U_R|^2)/2)/c_hat),   chi = (1-M_hat)^2
    mdot  = A * 0.5 * ( ro_L(Vn_p + Vn_hat_p_abs) + ro_R(Vn_m - Vn_hat_m_abs) - chi/c_hat * (P_R-P_L) )

符号: res_ro[owner] -= mdot, res_ro[neighbor] += mdot なので、owner から見て mdot > 0 は流出。

使い方:
  python3 diag_wall_cv_budget.py sern.h5 --auto 5                  # 最終 dump で ρ 最小の 5 節点を追う
  python3 diag_wall_cv_budget.py sern.h5 --nodes 1.2006,0.009,1.0025 --h 0.1
  python3 diag_wall_cv_budget.py sern.h5 --ids 123456 --faces      # 面別の内訳を出す
"""
import argparse, glob, re, sys
from pathlib import Path
import numpy as np
import h5py


def parse_struct(st, nplane):
    """PLANES/STRUCT = [nNodes, nodes..., nCells, cells...] の繰り返し → (owner, neighbor)。neighbor < 0 は境界面"""
    own = np.full(nplane, -1, np.int64); nei = np.full(nplane, -1, np.int64)
    i = ip = 0; n = len(st)
    while i < n and ip < nplane:
        i += 1 + int(st[i]); nc = int(st[i]); i += 1
        own[ip] = st[i]
        if nc >= 2: nei[ip] = st[i + 1]
        i += nc; ip += 1
    return own, nei, ip


def res_files(run_dir):
    fs = sorted(glob.glob(str(run_dir / "res_[0-9]*.h5")), key=lambda f: int(re.findall(r"\d+", Path(f).name)[0]))
    return fs + sorted(glob.glob(str(run_dir / "res_nan_*.h5")))


def load_state(fn, keys=("ro", "P", "Ux", "Uy", "Uz", "sonic", "T")):
    with h5py.File(fn, "r") as f:
        V = f["VALUE"]
        miss = [k for k in keys if k not in V]
        if miss: raise SystemExit(f"{fn}: VALUE/{miss} が無い (output.level を 1 以上に)")
        return {k: np.asarray(V[k][...], dtype=np.float64) for k in keys}


def slau_mdot(A, nx, ny, nz, L, R, wall_face=False):
    """1 次 SLAU の質量流束 (kg/s)。L/R は dict of scalars。

    `wall_face=True` のとき、**質量流束の chi だけ**を面法線成分で組み直す
    (`space.slauWallNormalChi: 1` と同型。`convectiveFlux_slau_d.inc.cuh:544-556`)。
    圧力束は本ツールでは扱わないので、chi_pressure は返さない。

    戻り: (mdot, 移流項, 圧力差項, chi_mass, c_hat, Vn_L, Vn_R)"""
    VnL = L["Ux"] * nx + L["Uy"] * ny + L["Uz"] * nz
    VnR = R["Ux"] * nx + R["Uy"] * ny + R["Uz"] * nz
    c_hat = 0.5 * (L["sonic"] + R["sonic"])
    Mp, Mm = VnL / c_hat, VnR / c_hat
    g = -max(min(Mp, 0.0), -1.0) * min(max(Mm, 0.0), 1.0)
    vh = (L["ro"] * abs(VnL) + R["ro"] * abs(VnR)) / (L["ro"] + R["ro"])
    vhp = (1.0 - g) * vh + g * abs(VnL)
    vhm = (1.0 - g) * vh + g * abs(VnR)
    u2L = L["Ux"] ** 2 + L["Uy"] ** 2 + L["Uz"] ** 2
    u2R = R["Ux"] ** 2 + R["Uy"] ** 2 + R["Uz"] ** 2
    M_hat = min(1.0, np.sqrt(0.5 * (u2L + u2R)) / c_hat)
    chi = (1.0 - M_hat) ** 2
    if wall_face:                                   # space.slauWallNormalChi: 1
        M_hat_n = min(1.0, np.sqrt(0.5 * (VnL ** 2 + VnR ** 2)) / c_hat)
        chi = (1.0 - M_hat_n) ** 2
    adv = 0.5 * (L["ro"] * (VnL + vhp) + R["ro"] * (VnR - vhm))
    pre = -0.5 * chi / c_hat * (R["P"] - L["P"])
    return A * (adv + pre), A * adv, A * pre, chi, c_hat, VnL, VnR


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", help="sern.h5 (幾何と PLANES)")
    ap.add_argument("--run-dir", default=".", help="res_*.h5 のあるディレクトリ (既定: mesh と同じ)")
    ap.add_argument("--ids", default="", help="CV (節点) 番号をカンマ区切りで")
    ap.add_argument("--nodes", default="", help="節点座標 x,y,z を ';' 区切りで (--h で割った無次元)。最近傍の節点を取る")
    ap.add_argument("--auto", type=int, default=0, help="最終 dump で ρ 最小の N 節点を選ぶ")
    ap.add_argument("--h", type=float, default=0.1, help="--nodes の無次元化長さ (既定 0.1 m = H)")
    ap.add_argument("--faces", action="store_true", help="面別の内訳を出す (既定は最終 dump のみ)")
    ap.add_argument("--faces-all", action="store_true", help="全 dump で面別の内訳を出す")
    ap.add_argument("--summary", action="store_true", help="dump ごとの ro/P 最小値と位置・床到達数だけを出す")
    ap.add_argument("--equilibrium", action="store_true",
                    help="隣接状態を凍結し、全接続面の Σmdot = 0 となる壁 CV の P_w を 1 次元求根で出し、観測値と比べる "
                         "(plan convection-slau-wall-normal-chi §6 V1-c)。壁ノードは T=Tw ピンなので ro_w = P_w/(R Tw)")
    ap.add_argument("--tw", type=float, default=1000.0, help="--equilibrium の壁温 [K]")
    ap.add_argument("--wall-normal-chi", action="store_true",
                    help="`space.slauWallNormalChi: 1` の run を解析するときに指定する。壁隣接面 (いずれかの端点が壁ノード) "
                         "の質量流束の chi を面法線 Mach で組む。**指定を忘れると別スキームの流束を計算する** "
                         "(codex result M1: 指定なしで壁 0/接線 1000/法線 0/Δp 900 の面は mdot=0、指定ありは −0.643 kg/s)")
    ap.add_argument("--wall-phys-ids", default="",
                    help="壁の physID をカンマ区切りで (例 1,2,3,4,10,11,12,13,15)。--wall-normal-chi に必須")
    ap.add_argument("--ro-min", type=float, default=1e-4, help="roMin (床到達数の判定用)")
    ap.add_argument("--p-min", type=float, default=20.0, help="pMin (床到達数の判定用)")
    a = ap.parse_args()

    mesh = Path(a.mesh); run_dir = Path(a.run_dir) if a.run_dir != "." else mesh.parent
    with h5py.File(mesh, "r") as f:
        S = np.asarray(f["PLANES/surfVect"], dtype=np.float64).reshape(-1, 3)
        st = np.asarray(f["PLANES/STRUCT"])
        vol = np.asarray(f["CELLS/volume"], dtype=np.float64)
        cc = np.asarray(f["CELLS/centCoords"], dtype=np.float64).reshape(-1, 3)
        xyz = np.asarray(f["MESH/COORD"], dtype=np.float64).reshape(-1, 3)
    own, nei, got = parse_struct(st, len(S))
    if got != len(S): raise SystemExit(f"PLANES/STRUCT を {got}/{len(S)} 面しか読めない")
    ncv = len(vol)
    A = np.linalg.norm(S, axis=1)
    print(f"# mesh {mesh}  CV {ncv}  faces {len(S)} (内部 {(nei >= 0).sum()}, 境界半割 {(nei < 0).sum()})")
    if len(xyz) != ncv:
        print(f"# WARNING: MESH/COORD {len(xyz)} != CV {ncv} — 座標指定は centCoords (双対重心) で行う")

    wallf = np.zeros(ncv, bool)
    if a.wall_normal_chi:
        if not a.wall_phys_ids:
            raise SystemExit("--wall-normal-chi には --wall-phys-ids が要る (カーネルと同じ壁ノード集合を作るため)")
        with h5py.File(mesh, "r") as f:
            for k in a.wall_phys_ids.split(","):
                k = k.strip()
                if k not in f["BCONDS"]: raise SystemExit(f"BCONDS/{k} が無い")
                ic = np.asarray(f["BCONDS"][k]["iCells"]).ravel()
                wallf[ic[(ic >= 0) & (ic < ncv)]] = True
        print(f"# slauWallNormalChi: 1 として解析 (壁ノード {int(wallf.sum())} / {ncv})")

    def is_wall_face(ip):
        o, n = int(own[ip]), int(nei[ip])
        return bool(a.wall_normal_chi and ((0 <= o < ncv and wallf[o]) or (0 <= n < ncv and wallf[n])))

    fs = res_files(run_dir)
    if not fs: raise SystemExit(f"{run_dir} に res_*.h5 が無い")
    if a.summary:
        print(f"{'dump':>16} {'ro_min':>11} {'at (x,y,z)/H':>26} {'P_min':>11} {'at (x,y,z)/H':>26} {'ro<=roMin':>10} {'P<=pMin':>9} {'nonfinite':>10}")
        ref = xyz if len(xyz) == ncv else cc
        for fn in fs:
            V = load_state(fn, keys=("ro", "P"))
            ok = np.isfinite(V["ro"]) & np.isfinite(V["P"])
            r = np.where(ok, V["ro"], np.nan); q = np.where(ok, V["P"], np.nan)
            i = int(np.nanargmin(r)); j = int(np.nanargmin(q))
            print(f"{Path(fn).name:>16} {r[i]:11.4e} {str(np.round(ref[i] / a.h, 4)):>26} {q[j]:11.4e} {str(np.round(ref[j] / a.h, 4)):>26}"
                  f" {int(np.nansum(r <= a.ro_min)):10d} {int(np.nansum(q <= a.p_min)):9d} {int((~ok).sum()):10d}")
        if not (a.ids or a.nodes or a.auto): return 0
    last = load_state(fs[-1])

    ids = [int(s) for s in a.ids.split(",") if s.strip()]
    if a.nodes:
        ref = xyz if len(xyz) == ncv else cc
        for tok in a.nodes.split(";"):
            p = np.array([float(v) for v in tok.split(",")]) * a.h
            ids.append(int(np.argmin(((ref - p) ** 2).sum(1))))
    if a.auto:
        ro = np.where(np.isfinite(last["ro"]), last["ro"], np.inf)
        ids += list(np.argsort(ro)[:a.auto])
    ids = list(dict.fromkeys(int(i) for i in ids))
    if not ids: raise SystemExit("--ids / --nodes / --auto のどれかを指定する")

    # 節点ごとの面リスト
    face_of = {i: [] for i in ids}
    sel = np.zeros(ncv, bool); sel[ids] = True
    for arr, role in ((own, "own"), (nei, "nei")):
        idx = np.nonzero((arr >= 0) & sel[np.clip(arr, 0, ncv - 1)])[0]
        for ip in idx:
            c = int(arr[ip])
            if c in face_of: face_of[c].append((int(ip), role))

    states = [(Path(fn).name, load_state(fn)) for fn in fs]
    for i in ids:
        pos = xyz[i] / a.h if len(xyz) == ncv else cc[i] / a.h
        print(f"\n=== CV {i}  node (x,y,z)/H = {np.round(pos, 5)}  dual centroid/H = {np.round(cc[i] / a.h, 5)}  V = {vol[i]:.4e} m^3  faces {len(face_of[i])} ===")
        print(f"{'dump':>16} {'ro':>11} {'P':>11} {'|U|':>9} {'c':>8} {'sum mdot(out)':>14} {'dro/dt':>12} {'(dro/dt)/ro':>12}")
        for name, V in states:
            if not np.isfinite(V["ro"][i]):
                print(f"{name:>16}   NaN"); continue
            tot = 0.0
            for ip, role in face_of[i]:
                j = int(nei[ip]) if role == "own" else int(own[ip])
                if j < 0: continue
                L = {k: V[k][i if role == "own" else j] for k in V}
                R = {k: V[k][j if role == "own" else i] for k in V}
                if not (np.isfinite(L["ro"]) and np.isfinite(R["ro"])): continue
                n = S[ip] / A[ip]
                m, _, _, _, _, _, _ = slau_mdot(A[ip], *n, L, R, is_wall_face(ip))
                tot += m if role == "own" else -m       # CV i から見た流出
            u = np.hypot(np.hypot(V["Ux"][i], V["Uy"][i]), V["Uz"][i])
            drodt = -tot / vol[i]
            print(f"{name:>16} {V['ro'][i]:11.4e} {V['P'][i]:11.4e} {u:9.1f} {V['sonic'][i]:8.1f} {tot:14.5e} {drodt:12.4e} {drodt / max(V['ro'][i], 1e-30):12.4e}")
        if a.equilibrium:
            V = states[-1][1]
            if np.isfinite(V["ro"][i]):
                R_gas = float(V["P"][i] / max(V["ro"][i] * V["T"][i], 1e-30))
                fl = [(ip, role, int(nei[ip]) if role == "own" else int(own[ip])) for ip, role in face_of[i]]
                fl = [(ip, role, j) for ip, role, j in fl if j >= 0 and np.isfinite(V["ro"][j])]

                def net(Pw):
                    """壁 CV の圧力を Pw としたときの全接続面の正味流出 [kg/s] (隣接状態は凍結)"""
                    row = 1.0 if R_gas <= 0 else Pw / (R_gas * a.tw)
                    tot = 0.0
                    for ip, role, j in fl:
                        self_st = {"ro": row, "P": Pw, "Ux": 0.0, "Uy": 0.0, "Uz": 0.0,
                                   "sonic": V["sonic"][i]}
                        oth = {k: V[k][j] for k in V}
                        L, Rr = (self_st, oth) if role == "own" else (oth, self_st)
                        n = S[ip] / A[ip]
                        m, _, _, _, _, _, _ = slau_mdot(A[ip], *n, L, Rr, is_wall_face(ip))
                        tot += m if role == "own" else -m
                    return tot

                lo, hi = 1.0, max(V["P"][j] for _, _, j in fl) * 1.5 if fl else 1.0
                flo, fhi = net(lo), net(hi)
                obs = V["P"][i]
                print(f"  -- 平衡壁圧 (隣接凍結, 全 {len(fl)} 面, R={R_gas:.1f}, Tw={a.tw:g}) --")
                if flo * fhi > 0:
                    print(f"     求根できず (net({lo:.3g})={flo:+.3e}, net({hi:.3g})={fhi:+.3e} が同符号)")
                else:
                    for _ in range(80):
                        mid = 0.5 * (lo + hi)
                        if net(lo) * net(mid) <= 0: hi = mid
                        else: lo = mid
                    pred = 0.5 * (lo + hi)
                    r = obs / pred if pred > 0 else float("inf")
                    ok = (1 / 1.5) <= r <= 1.5
                    print(f"     予測 P_w = {pred:10.4e} Pa   観測 P_w = {obs:10.4e} Pa   観測/予測 = {r:.3f}"
                          f"   -> {'PASS (×/÷1.5 以内)' if ok else 'FAIL'}")
        if a.faces or a.faces_all:
            for name, V in (states if a.faces_all else states[-1:]):
                if not np.isfinite(V["ro"][i]): continue
                print(f"  -- 面別 ({name}) — mdot > 0 は CV {i} から流出 --")
                print(f"  {'face':>9} {'other':>9} {'A[m2]':>10} {'n·(out)':>24} {'mdot':>12} {'adv':>12} {'press':>12} {'chi':>7} {'c_hat':>8} {'Vn_self':>9} {'Vn_oth':>9} {'P_self':>10} {'P_oth':>10} {'ro_oth':>10}")
                rows = []
                for ip, role in face_of[i]:
                    j = int(nei[ip]) if role == "own" else int(own[ip])
                    if j < 0: continue
                    L = {k: V[k][i if role == "own" else j] for k in V}
                    R = {k: V[k][j if role == "own" else i] for k in V}
                    n = S[ip] / A[ip]
                    m, adv, pre, chi, ch, VnL, VnR = slau_mdot(A[ip], *n, L, R, is_wall_face(ip))
                    s = 1.0 if role == "own" else -1.0      # CV i から見た流出符号
                    rows.append((s * m, ip, j, A[ip], s * n, s * m, s * adv, s * pre, chi, ch,
                                 s * (VnL if role == "own" else VnR), s * (VnR if role == "own" else VnL),
                                 V["P"][i], V["P"][j], V["ro"][j]))
                for r in sorted(rows, key=lambda r: -r[0]):
                    _, ip, j, Af, nn, m, adv, pre, chi, ch, vs, vo, ps, po, roo = r
                    print(f"  {ip:9d} {j:9d} {Af:10.3e} ({nn[0]:6.3f},{nn[1]:6.3f},{nn[2]:6.3f}) {m:12.4e} {adv:12.4e} {pre:12.4e} {chi:7.4f} {ch:8.1f} {vs:9.1f} {vo:9.1f} {ps:10.3e} {po:10.3e} {roo:10.3e}")
                nb = sum(1 for ip, _ in face_of[i] if nei[ip] < 0)
                if nb: print(f"  (境界半割面 {nb} 面は convectiveFlux_boundary が別に扱う。no-slip 壁なら質量流束 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
