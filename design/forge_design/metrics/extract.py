"""res_*.h5 からの目的関数抽出 (③ベル最小セット)。

規約 (サーベイ B4.5 / case/29 postprocess_thrust.py):
- 出口面の状態は境界ダンプでなく res h5 を BCONDS/<outlet>/iCells で引く。
- 軸対称の PLANES/surfArea は 2D 辺長 (dr) — dA = 2πr·dr を手で構成する。
- P は roe から再構成 (CPG)。
- 全指標は最小化形でなく生値も併記 (η は最大化量なので optimizer へは -η)。

TP (③ベルの燃焼ガス) の推力は `thrust_cv` (閉じた CV の壁会計・出口面会計と CV 残差、plan
tooling-nozzle-bell-two-mode §4.1・§5.1 #2b。Euler の slip 壁は `wall_model="slip"`、同 §5.1 #6 (W1))。
`thrust_metrics` (CPG の出口面法) は従来どおり。
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from ..evaluate.ic import invert_area_ratio


def thrust_metrics(mesh_h5, res_h5, gamma, Pt, p_ambient, r_throat,
                   outlet_physid=2) -> dict:
    """推力・C_F・η (= C_F / C_F,ideal 同面積比)。

    出口面の運動量法 F = ∫(ρu² + (P - Pa)) dA。壁摩擦・壁圧の影響は出口状態に
    既に織り込まれているため**別加算しない** (加算すると二重計上 — サーベイ
    B4.5 の「+壁摩擦」はこの抽出法では不採用。子 plan §9 参照)。
    """
    g = gamma
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf = pc[ip, 1]
        dr = nz["/PLANES/surfArea"][:][ip]  # 軸対称: 2D 辺長 (dr)
        dA = 2.0 * np.pi * rf * dr
    with h5py.File(res_h5, "r") as f:
        ro = f["/VALUE/ro"][:][ic]
        u = f["/VALUE/roUx"][:][ic] / ro
        v = f["/VALUE/roUy"][:][ic] / ro
        roe = f["/VALUE/roe"][:][ic]
        P = (g - 1.0) * (roe - 0.5 * ro * (u * u + v * v))
    F = float(np.sum((ro * u * u + (P - p_ambient)) * dA))
    At = np.pi * r_throat ** 2
    Ae = float(np.sum(dA))
    CF = F / (Pt * At)
    CF_id = cf_ideal(Ae / At, g, p_ambient / Pt)
    return {
        "thrust_N": F,
        "CF": CF,
        "CF_ideal": CF_id,
        "eta_cf": CF / CF_id,
        "eps_measured": Ae / At,
        "mdot_kg_s": float(np.sum(ro * u * dA)),
    }


# ---------------------------------------------------------------- 推力の CV 会計 (plan tooling-nozzle-bell-two-mode §4.1, §5.1 #2b)
#
# ノズル内壁・入口面・出口面で閉じた検査体積 CV の全境界を、外向き法線 n で
#     Π_x = ∫[ρ u_x (u·n) + (p − p_a) n_x − (τ_eff n)_x] dA,   dA = 2π r ds (軸対称)
# に統一する。定常で Σ Π_x = 0 (軸対称の x 運動量には源項が無い。軸は r = 0 なので面積 0) なので
#     出口面会計  F_exit = Π_x(出口)                     (出口面の運動量流束 + 圧力 − 粘性応力)
#     壁会計      F_wall = −Π_x(入口) − Π_x(内壁)         (入口面の運動量流束・圧力・粘性応力 + 内壁の圧力と摩擦)
# で、CV 残差 = (F_wall − F_exit)/F_1D = −Σ Π_x / F_1D。推力は「流れが +x へ出るときに正」。
# 内壁の粘性 traction は forge の twall_* から作る。forge の規約は離散化で逆 (metrics/sern_forces.py:73 と
# viscousFlux_d.cu): cell は `tau = 面の粘性流束` を残差に足す = **流体に働く** traction、node (既定の
# `wallStressForOutput_node_d`) は内部双対面から壁 CV に入る粘性運動量 = **壁に働く**力 (流体に働く traction の逆符号)。
# Euler (slip 壁、2026-10-12 plan bell-two-mode §5.1 #6 (W1)): 式は同じで粘性応力が 0。forge の slip 境界は twall_* を持たない
# (boundaryCond.hpp の slip の valueTypes) ので、呼び出し側が `wall_model="slip"` を明示したときだけ内壁の traction を 0 とし
# 壁ダンプを読まない (NS で twall が欠けていたら止める — 黙って 0 にしない)。開境界の粘性項は式のまま (viscMethod 0 は
# vis_lam = visc = 0 なので 0 になる)。slip は場の μ_eff (vis_lam + vis_turb) が境界の DOF で 0 であることを確かめてから使う。
CV_IDS = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4}   # meshing.mesh2d.PHYS と同じ
WALL_MODELS = ("viscous", "slip")                           # thrust_cv の内壁の扱い (既定 viscous = 従来の NS)


def _step_of(path) -> int:
    """res_<step>.h5 / res_<群>_<id>_<step>.h5 の step (末尾の数字)。"""
    import re
    m = re.search(r"_(\d+)\.h5$", str(path))
    if not m:
        raise ValueError(f"step を読めないファイル名: {path}")
    return int(m.group(1))


def tp_pressure_from_conserved(ro, roe, ke, gas) -> np.ndarray:
    """TP (凍結組成の FrozenGas) の保存量から P を組む: e = (ρE − ρu²/2)/ρ を forge の thermoHrefTemp 基準の
    顕内部エネルギー e_sens(T) で引いて T、P = ρRT。CPG の (γ−1)(ρE − ρu²/2) は使わない。"""
    Tg = np.linspace(50.0, 5999.0, 60000)
    eg = gas.e_sens(Tg)
    e = (np.asarray(roe, float) - np.asarray(ke, float)) / np.maximum(np.asarray(ro, float), 1e-30)
    T = np.interp(e, eg, Tg)
    return np.asarray(ro, float) * float(gas.R) * T


def _pressure(V, idx, gas=None, gamma=None) -> tuple:
    """DOF idx の静圧と出所。ソルバの原始量 `VALUE/P` (output level ≥ 1) を優先し、無ければ保存量から組む
    (TP: `gas` の e_sens(T) の逆算 / CPG: `gamma`)。どちらも無ければ ValueError。"""
    if "P" in V:
        return np.asarray(V["P"][:][idx], float), "VALUE/P"
    ro = V["ro"][:][idx].astype(float)
    ke = 0.5 * (V["roUx"][:][idx].astype(float) ** 2 + V["roUy"][:][idx].astype(float) ** 2
                + V["roUz"][:][idx].astype(float) ** 2) / np.maximum(ro, 1e-30)
    roe = V["roe"][:][idx].astype(float)
    if gas is not None:
        return tp_pressure_from_conserved(ro, roe, ke, gas), "TP e_sens(T) の逆算 (VALUE/P なし)"
    if gamma is not None:
        return (float(gamma) - 1.0) * (roe - ke), "CPG (γ−1)(ρE − ρu²/2) (VALUE/P なし)"
    raise ValueError("VALUE/P が無く、保存量から P を組む気体 (gas / gamma) も与えられていない")


def _halfedge_geometry(m, phys_id: int) -> dict:
    """node (median-dual): 境界 physID の主格子の辺を両端の節点側の半辺に分け、半辺ごとに
    (所有節点, 外向き単位法線 n, r 重みの面積 ∫ r ds) を返す。直線の半辺では ∫ r ds = L·r_mid が厳密。
    外向きは 2 通りで決めて照合する: (a) 辺を含む主格子の四角形 (VIZMESH/CONNE) の重心から離れる向き、
    (b) 変換器が書いた境界双対面の外向きベクトル (`PLANES/surfVect` を bcond の `iPlanes`/`iCells` で端点に寄せたもの)。
    両方あって 1 辺でも食い違えば ValueError、片方だけならそれを使い、どちらも無ければ ValueError
    (双対 CV 重心からの推定は高 AR の曲がった壁で向きを誤るので使わない)。
    角の節点 (入口∩壁など) も半辺の所有節点の状態で評価するので、境界条件の帰属 (nodeInletCornerWall) に依らず
    閉曲面をちょうど 1 回覆う。"""
    xy = m["MESH/COORD"][:].reshape(-1, 3)[:, :2].astype(float)
    e = m[f"BCONDS/{phys_id}/vizBfaceNodes"][:].reshape(-1, 2).astype(np.int64)
    sizes = m[f"BCONDS/{phys_id}/vizBfaceSizes"][:]
    if not np.all(sizes == 2):
        raise ValueError(f"BCONDS/{phys_id}: 2 節点の辺以外の境界面がある (2D 軸対称の node 専用)")
    a, b = e[:, 0], e[:, 1]
    N = len(xy)
    mid = 0.5 * (xy[a] + xy[b])
    d = xy[b] - xy[a]
    L = np.hypot(d[:, 0], d[:, 1])
    if np.any(L <= 0):
        raise ValueError(f"BCONDS/{phys_id}: 長さ 0 の境界辺がある")
    n = np.column_stack([d[:, 1], -d[:, 0]]) / L[:, None]
    s_quad = None                                   # (a) 四角形の重心から離れる向き
    if "VIZMESH/CONNE" in m:
        cn = m["VIZMESH/CONNE"][:]
        if len(cn) % 5 == 0 and np.all(cn[::5] == 5):
            q = cn.reshape(-1, 5)[:, 1:].astype(np.int64)
            qc = xy[q].mean(axis=1)
            ek = np.concatenate([np.minimum(q[:, k], q[:, (k + 1) % 4]) * N + np.maximum(q[:, k], q[:, (k + 1) % 4])
                                 for k in range(4)])
            eq = np.concatenate([np.arange(len(q))] * 4)
            bk = np.minimum(a, b) * N + np.maximum(a, b)
            o = np.argsort(ek, kind="stable")
            eks = ek[o]
            pos = np.minimum(np.searchsorted(eks, bk), len(eks) - 1)
            if np.all(eks[pos] == bk):
                s_quad = np.sign(np.einsum("ij,ij->i", n, mid - qc[eq[o][pos]]))
    s_sv = None                                     # (b) 境界双対面の外向きベクトル (端点の和)
    g = f"BCONDS/{phys_id}"
    if "PLANES/surfVect" in m and f"{g}/iPlanes" in m and f"{g}/iCells" in m:
        ip, ic = m[f"{g}/iPlanes"][:], m[f"{g}/iCells"][:].astype(np.int64)
        svn = np.zeros((N, 2))
        np.add.at(svn, ic, m["PLANES/surfVect"][:].reshape(-1, 3)[ip, :2].astype(float))
        ref = svn[a] + svn[b]
        dot = np.einsum("ij,ij->i", n, ref)
        ok = np.abs(dot) > 0.1 * np.hypot(ref[:, 0], ref[:, 1])
        s_sv = np.where(ok, np.sign(dot), 0.0)
    if s_quad is None and (s_sv is None or np.any(s_sv == 0)):
        raise ValueError(f"BCONDS/{phys_id}: 境界辺の外向きを決められない (VIZMESH/CONNE も境界双対面のベクトルも使えない)")
    if s_quad is not None and np.any(s_quad == 0):
        raise ValueError(f"BCONDS/{phys_id}: 四角形の重心が境界辺の上にある (向きを決められない)")
    if s_quad is not None and s_sv is not None:
        bad = (s_sv != 0) & (s_sv != s_quad)
        if np.any(bad):
            raise ValueError(f"BCONDS/{phys_id}: 境界辺の外向きが四角形と境界双対面で {int(bad.sum())} 辺食い違う")
    n = n * (s_quad if s_quad is not None else s_sv)[:, None]
    owner = np.concatenate([a, b])
    nn = np.concatenate([n, n])
    # 半辺: [a, mid] と [mid, b]。r 重み = (L/2)·(r_端 + r_mid)/2
    rA = 0.5 * L * 0.5 * (xy[a, 1] + mid[:, 1])
    rB = 0.5 * L * 0.5 * (xy[b, 1] + mid[:, 1])
    return {"owner": owner, "n": nn, "rds": np.concatenate([rA, rB]), "n_edges": int(len(e))}


def _face_geometry_cell(m, phys_id: int) -> dict:
    """cell: 境界の主格子面 (iPlanes) ごとに (owner セル, 外向き単位法線, r_c·|S|)。面は直線なので r_c·|S| が厳密。"""
    ip = m[f"BCONDS/{phys_id}/iPlanes"][:]
    ic = m[f"BCONDS/{phys_id}/iCells"][:]
    sv = m["PLANES/surfVect"][:].reshape(-1, 3)[ip, :2].astype(float)
    pc = m["PLANES/centCoords"][:].reshape(-1, 3)[ip]
    S = np.hypot(sv[:, 0], sv[:, 1])
    return {"owner": ic.astype(np.int64), "n": sv / S[:, None], "rds": pc[:, 1].astype(float) * S, "n_edges": int(len(ip))}


def _wall_traction_on_fluid(wall_dump, dof_xyz, wall_dofs, node_mode: bool) -> tuple:
    """壁面ダンプ res_wall_<id>_<step>.h5 の twall_x/y を DOF に座標の厳密一致で対応づけ、**流体に働く** traction を
    返す (node: −twall、cell: +twall; 本節冒頭の規約)。戻り値 (traction (n,2) — wall_dofs の順, 規約の文字列)。"""
    from scipy.spatial import cKDTree
    with h5py.File(wall_dump, "r") as f:
        if "VALUE/twall_x" not in f or "VALUE/twall_y" not in f:
            raise ValueError(f"{wall_dump}: 壁ダンプに VALUE/twall_x・twall_y が無い (slip 壁のダンプ?) — 摩擦を 0 とみなさない。"
                             "Euler (slip 壁) の run なら thrust_cv(..., wall_model='slip') を明示する")
        wxyz = f["MESH/COORD"][:].reshape(-1, 3)[:, :2].astype(float)
        tw = np.column_stack([f["VALUE/twall_x"][:], f["VALUE/twall_y"][:]]).astype(float)
        if not node_mode:
            cn = f["MESH/CONNE"][:].reshape(-1, 4)[:, 2:4]
            wxyz = 0.5 * (wxyz[cn[:, 0]] + wxyz[cn[:, 1]])          # cell: 面値 → 面中心で対応づける
    if len(wxyz) != len(tw):
        raise ValueError(f"{wall_dump}: 座標 {len(wxyz)} と twall {len(tw)} の数が合わない")
    d, j = cKDTree(wxyz).query(dof_xyz[wall_dofs])
    span = float(np.ptp(dof_xyz[:, 0])) or 1.0
    if float(d.max()) > 1e-6 * span:
        raise ValueError(f"{wall_dump}: 壁ダンプの点を壁 DOF に対応づけられない (最大 {float(d.max()):.3e} m)")
    if len(np.unique(j)) != len(j):
        raise ValueError(f"{wall_dump}: 壁ダンプの同じ点に複数の壁 DOF が対応した (座標の重複)")
    sgn = -1.0 if node_mode else 1.0
    return sgn * tw[j], ("node: 流体に働く traction = −twall (壁 CV に入る粘性運動量)" if node_mode
                         else "cell: 流体に働く traction = +twall (残差に足す粘性流束)")


def thrust_cv(mesh_h5, res_h5, p_ambient: float, F_ideal: float, *, wall_dump=None, gas=None, gamma=None,
              iso_k: bool = False, ids: dict | None = None, mdot_ideal: float | None = None,
              wall_model: str = "viscous") -> dict:
    r"""閉じた CV の 2 会計 (壁会計・出口面会計) と CV 残差 (本節冒頭の定義)。

    - 状態は境界の DOF 値 (node: 半辺の所有節点、cell: 面の owner セル)。P は `VALUE/P` を優先 (`_pressure`)。
    - 粘性応力: 内壁は twall (`wall_dump`、無ければ res と同じ step の `res_wall_<id>_<step>.h5`)。入口・出口面は
      $\mu_{eff}=\mu+\mu_t$ (`vis_lam`+`vis_turb`) とソルバの勾配 (`dUxdx, dUxdy, dUydx, dUydy, axisym_divU`;
      output.extraFields) から $\tau_{xx}=\mu_{eff}(2\partial_x u-\tfrac23\nabla\!\cdot\!u)$、
      $\tau_{xr}=\mu_{eff}(\partial_r u+\partial_x v)$。`iso_k` (`sstIsotropicStress: 1`) なら $-\tfrac23\rho k$ を足す。
      勾配が res に無いときは開境界の粘性項を評価せず `open_viscous: "missing"`・`cv_closable: False` にする。
    - `F_ideal`: η と CV 残差の分母 (TP: 同じ凍結 TP の $F_{1D}(\varepsilon;p_a)$)。`mdot_ideal` は ṁ/ṁ_1D 用 (任意)。
    - `wall_model`: `"viscous"` (既定、従来どおり。壁ダンプに twall_* が無ければ ValueError) / `"slip"` (Euler の slip 壁:
      内壁の traction を 0 とし壁ダンプを読まない。場の μ_eff = vis_lam + vis_turb が境界の全 DOF で 0 でなければ ValueError —
      粘性の場で摩擦を黙って 0 にしない)。使った経路は戻り値の `wall_model` と `wall_traction_convention` に残す。
    戻り値の推力は [N] (2π を含む全周)。η は出口面会計 `eta_exit` と壁会計 `eta_wall` の両方を返す。"""
    if wall_model not in WALL_MODELS:
        raise ValueError(f"thrust_cv: wall_model '{wall_model}' は未知 ({list(WALL_MODELS)})")
    ids = {**CV_IDS, **(ids or {})}
    with h5py.File(mesh_h5, "r") as m, h5py.File(res_h5, "r") as f:
        V = f["VALUE"]
        n_val = len(V["ro"])
        node_mode = n_val == len(m["MESH/COORD"][:]) // 3
        geo = {k: (_halfedge_geometry(m, ids[k]) if node_mode else _face_geometry_cell(m, ids[k]))
               for k in ("inlet", "outlet", "wall")}
        dof_xyz = (m["MESH/COORD"][:].reshape(-1, 3) if node_mode else m["CELLS/centCoords"][:].reshape(-1, 3))[:, :2]
        wall_dofs = np.unique(m[f"BCONDS/{ids['wall']}/iCells"][:])
        allidx = np.unique(np.concatenate([geo[k]["owner"] for k in geo]))
        P_all, p_src = _pressure(V, allidx, gas=gas, gamma=gamma)
        P = np.full(n_val, np.nan); P[allidx] = P_all
        ro = V["ro"][:].astype(float)
        if "Ux" in V:
            ux, ur = V["Ux"][:].astype(float), V["Uy"][:].astype(float)
        else:
            ux, ur = V["roUx"][:].astype(float) / ro, V["roUy"][:].astype(float) / ro
        grads = ("dUxdx", "dUxdy", "dUydx", "dUydy", "axisym_divU")
        have_grad = all(g in V for g in grads) and "vis_lam" in V and "vis_turb" in V
        G = {g: V[g][:].astype(float) for g in grads} if have_grad else None
        mu = (V["vis_lam"][:].astype(float) + V["vis_turb"][:].astype(float)) if have_grad else None
        rk = (ro * V["k"][:].astype(float)) if (iso_k and "k" in V) else None
        if wall_model == "slip":
            # slip (Euler) の前提: 境界の全 DOF で μ_eff = 0 (viscMethod 0・visc 0・乱流なし)。粘性の場なら止める
            if "vis_lam" not in V or "vis_turb" not in V:
                raise ValueError("thrust_cv(wall_model='slip'): res に vis_lam・vis_turb が無く、μ_eff = 0 (非粘性) を確かめられない")
            mu_b = V["vis_lam"][:][allidx].astype(float) + V["vis_turb"][:][allidx].astype(float)
            if not np.all(mu_b == 0.0):
                raise ValueError(f"thrust_cv(wall_model='slip'): 境界の μ_eff が 0 でない (最大 {float(np.max(np.abs(mu_b))):.3e} Pa·s) — "
                                 "粘性の場で壁の摩擦を 0 にしない (slip は Euler の run だけ)")
    if iso_k and rk is None:
        raise ValueError("iso_k (sstIsotropicStress: 1) だが res に VALUE/k が無い")
    two_pi = 2.0 * np.pi
    out = {"res": str(res_h5), "node_mode": bool(node_mode), "p_ambient": float(p_ambient), "F_ideal": float(F_ideal),
           "pressure_source": p_src, "open_viscous": "evaluated" if have_grad else "missing", "iso_k": bool(iso_k),
           "wall_model": wall_model}

    def open_face(key):
        g = geo[key]; o = g["owner"]; n = g["n"]; dA = two_pi * g["rds"]
        un = ux[o] * n[:, 0] + ur[o] * n[:, 1]
        mom = float(np.sum(ro[o] * ux[o] * un * dA))
        prs = float(np.sum((P[o] - p_ambient) * n[:, 0] * dA))
        vis = None
        if have_grad:
            txx = mu[o] * (2.0 * G["dUxdx"][o] - (2.0 / 3.0) * G["axisym_divU"][o])
            txr = mu[o] * (G["dUxdy"][o] + G["dUydx"][o])
            if rk is not None:
                txx = txx - (2.0 / 3.0) * rk[o]
            vis = float(np.sum((txx * n[:, 0] + txr * n[:, 1]) * dA))   # ∫ (τ n)_x dA
        mdot = float(np.sum(ro[o] * un * dA))
        Pi = mom + prs - (vis if vis is not None else 0.0)
        return {"Pi_x": Pi, "momentum": mom, "pressure": prs, "viscous_tau_n_x": vis, "mdot_out": mdot,
                "area": float(np.sum(dA)), "n_faces": g["n_edges"]}

    inl = open_face("inlet")
    ext = open_face("outlet")
    # 内壁: u = 0 (no-slip、node は Dirichlet) なので運動量流束は 0 とし、圧力と粘性 traction だけを積む
    gw = geo["wall"]; o = gw["owner"]; n = gw["n"]; dA = two_pi * gw["rds"]
    if not np.all(np.isin(o, wall_dofs)):
        # 壁の辺の両端は壁 DOF のはず (入口∩壁・壁∩出口の角も壁 bcond の DOF)
        raise ValueError("内壁の半辺の所有 DOF が壁 bcond の DOF に無い (メッシュの境界構成が想定外)")
    w_prs = float(np.sum((P[o] - p_ambient) * n[:, 0] * dA))
    w_un = float(np.max(np.abs(ux[o] * n[:, 0] + ur[o] * n[:, 1]))) if len(o) else 0.0
    if wall_model == "slip":
        # Euler の slip 壁: 内壁の粘性 traction は 0 (壁ダンプは読まない — slip の壁ダンプに twall は無い)。運動量流束は no-slip と
        # 同じく積まない (slip の数値流束は圧力だけ。壁 DOF の |u·n| の最大は記録)
        w_vis = 0.0
        out["wall_traction_convention"] = "slip (Euler): 内壁の粘性 traction = 0 (壁ダンプを読まない、μ_eff = 0 を確認)"
    else:
        if wall_dump is None:
            cand = Path(res_h5).parent / f"res_wall_{ids['wall']}_{_step_of(res_h5)}.h5"
            wall_dump = cand if cand.exists() else None
        if wall_dump is not None:
            trac, conv = _wall_traction_on_fluid(wall_dump, dof_xyz, wall_dofs, node_mode)
            tmap = np.full((n_val, 2), np.nan); tmap[wall_dofs] = trac
            w_vis = float(np.sum(tmap[o, 0] * dA))      # ∫ (τ n)_x dA (twall は traction ベクトルそのもの: その x 成分)
            out["wall_traction_convention"] = conv
            out["wall_dump"] = str(wall_dump)
        else:
            w_vis = None
            out["wall_traction_convention"] = "壁面ダンプ無し (摩擦を評価できない)"
    wal = {"Pi_x": w_prs - (w_vis if w_vis is not None else 0.0), "pressure": w_prs, "viscous_tau_n_x": w_vis,
           "area": float(np.sum(dA)), "n_faces": gw["n_edges"], "max_abs_u_n": w_un}
    F_exit = ext["Pi_x"]
    F_wall = -inl["Pi_x"] - wal["Pi_x"]
    if wall_model == "slip":
        # μ_eff = 0 を確かめたので、開境界の粘性項も勾配の有無に依らず 0 (勾配があれば式で評価した値 = 0 を記録)
        closable = True
        if not have_grad:
            out["open_viscous"] = "zero (μ_eff = 0、勾配は res に無い)"
    else:
        closable = have_grad and (w_vis is not None)
    mdot_in = -inl["mdot_out"]
    out.update({
        "groups": {"inlet": inl, "outlet": ext, "wall": wal},
        "thrust_exit_N": F_exit, "thrust_wall_N": F_wall,
        "friction_N": w_vis,                                   # 内壁の摩擦の推力寄与 (= +∫(τn)_x dA、付着流では負)
        "eta_exit": F_exit / F_ideal, "eta_wall": F_wall / F_ideal,
        "cv_residual": (F_wall - F_exit) / F_ideal, "cv_closable": bool(closable),
        "mdot_in_kg_s": mdot_in, "mdot_exit_kg_s": ext["mdot_out"],
        "mass_imbalance": (ext["mdot_out"] - mdot_in) / max(abs(mdot_in), 1e-30),
        "exit_area_m2": ext["area"],
    })
    if mdot_ideal is not None:
        out["mdot_ideal_kg_s"] = float(mdot_ideal)
        out["mdot_ratio_1d"] = ext["mdot_out"] / float(mdot_ideal)
    return out


def cf_ideal(eps, gamma, pa_over_pt) -> float:
    """1D 等エントロピーの理想推力係数 (同面積比・同背圧)。"""
    g = gamma
    Me = float(invert_area_ratio(np.array([eps]), np.array([True]), g)[0])
    pe_pt = (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-g / (g - 1.0))
    gterm = np.sqrt(
        (2.0 * g * g / (g - 1.0))
        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
        * (1.0 - pe_pt ** ((g - 1.0) / g))
    )
    return float(gterm + (pe_pt - pa_over_pt) * eps)


def core_radius_traced(mesh_h5, res_h5, x_d, x_plane, r_start_frac=0.02,
                       ds_frac=2e-4, gamma=1.4) -> float:
    r"""**有効菱形の上流境界**を実場でトレースして求める (テストコア半径)。

    軸が $M_d$ に達する点 $(x_d,0)$ から出る C⁺ ($dr/dx=\tan(\theta+\mu)$) を
    **実際の CFD 場**の中で前進積分し、`x_plane` での半径を返す。この線より下が
    一様域 (有効菱形の内側)。

    直線近似 $(x-x_d)\tan\mu_d$ は特性線が**一様域でしか直線でない**ことを無視する
    ため近似にすぎない (本ケースでは 0.824 vs 実トレース 0.830 で 0.7% 差)。

    **cell/node の違い** (2026-08-17): cell は軸上に DOF が無く補間の凸包が
    $r_{min}$ で切れるため `r_start_frac`×r_max から開始せざるを得ない。
    node は軸上に DOF を持つので、補間点を **`MESH/COORD` (真の節点座標)** に
    すれば $r=0$ から開始できる。`/CELLS/centCoords` は node では双対 CV 重心に
    書き換わり最小 $r$ が 0 でない (実測 0.0089 $r_t$) ので、node では使わない。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        node_c = (nz["/MESH/COORD"][:].reshape(-1, 3) if "/MESH/COORD" in nz else None)
    with h5py.File(res_h5, "r") as f:
        Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
    from scipy.interpolate import LinearNDInterpolator
    if node_c is not None and len(node_c) == len(Ux):
        cc = node_c                      # node: 真の節点座標 (軸 r=0 を含む)
        r_start_frac = 0.0
    x, r = cc[:, 0], cc[:, 1]
    itp_M = LinearNDInterpolator(np.c_[x, r], np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    itp_th = LinearNDInterpolator(np.c_[x, r], np.arctan2(Uy, Ux))
    r_max = float(r.max())
    xx, rr = float(x_d), max(r_start_frac * r_max, float(r.min()) * 1.05)
    ds = ds_frac * float(x.max() - x.min())
    for _ in range(200000):
        if xx >= x_plane or rr >= r_max:
            break
        Mv, tv = float(itp_M(xx, rr)), float(itp_th(xx, rr))
        if not (np.isfinite(Mv) and np.isfinite(tv)) or Mv <= 1.0:
            break
        mu = np.arcsin(1.0 / max(Mv, 1.0 + 1e-12))
        rr += ds * np.tan(tv + mu)
        xx += ds
    return float(min(rr, r_max))


def core_radius_straight(x_plane: float, x_d: float, M_design: float,
                         r_wall: float) -> float:
    r"""有効菱形上流境界の**直線近似** $(x-x_d)\tan\mu_d$ (壁でクリップ)。

    一様域外では特性線が曲がるので近似。トレース版との差は本ケースで 0.7%。"""
    mu = np.arcsin(1.0 / max(float(M_design), 1.0 + 1e-12))
    return float(min(max((float(x_plane) - float(x_d)) * np.tan(mu), 0.0), float(r_wall)))


def exit_uniformity(mesh_h5, res_h5, M_design, x_d=None, outlet_physid=2,
                    core_mode: str = "traced", core_radius=None,
                    gamma=1.4, n_grid: int | None = None) -> dict:
    r"""出口一様性 $\varepsilon_M$ / $\varepsilon_\theta$ (サーベイ B4.5 定義)。

    **幾何は出口境界面の実 connectivity から取る** (`BCONDS/<id>/iPlanes` +
    `PLANES/surfArea`)。軸対称では面中心が半径中点なので
    $2\pi r_c\,dr = \pi(r_{out}^2-r_{in}^2)$ が**厳密**に成り立ち、環状面積は
    近似なしで得られる (セル中心から仮想境界を作る旧実装は近似だった)。
    壁半径も最外面の $r_{out}$ = **真の壁半径** (最外セル中心ではない)。

    - $\varepsilon_M$ = コア上の**質量流束重み RMS 偏差**
      $\sqrt{\langle (M-M_d)^2\rangle_w}/M_d$、$w_i=\rho_i u_{x,i} A_i$。
    - $\varepsilon_\theta$ = コア上の $\max|\theta|$ [deg] (独立の目的)。

    `core_mode`: `"traced"` (既定, 実場で C⁺ をトレース) / `"straight"` (直線近似)
    / `"full"` (コア制限なし) / `core_radius` 明示。
    `n_grid` を指定すると $r/r_{wall}$ 固定格子へ補間して積分し直す
    (サーベイの別法。厳密面積との一致を確認する交差検証用)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf, xf = pc[ip, 1], pc[ip, 0]
        dr = nz["/PLANES/surfArea"][:][ip]
    with h5py.File(res_h5, "r") as f:
        Ux, Uy = f["/VALUE/Ux"][:][ic], f["/VALUE/Uy"][:][ic]
        son, ro = f["/VALUE/sonic"][:][ic], f["/VALUE/ro"][:][ic]
    o = np.argsort(rf)
    rf, dr, Ux, Uy, son, ro = rf[o], dr[o], Ux[o], Uy[o], son[o], ro[o]
    A = 2.0 * np.pi * rf * dr              # = π(r_out²−r_in²) (厳密)
    r_wall = float(rf[-1] + 0.5 * dr[-1])  # 真の壁半径 (最外面の外縁)
    x_plane = float(np.mean(xf))
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th_deg = np.degrees(np.arctan2(Uy, Ux))
    Md = float(M_design)

    if core_radius is not None:
        rc, how = float(core_radius), "explicit"
    elif core_mode == "full" or x_d is None:
        rc, how = r_wall, "full(no_core_limit)"
    elif core_mode == "straight":
        rc, how = core_radius_straight(x_plane, x_d, Md, r_wall), "rhombus_straight"
    else:
        rc, how = core_radius_traced(mesh_h5, res_h5, x_d, x_plane, gamma=gamma), "rhombus_traced"
    core = rf <= rc
    if core.sum() < 4:
        raise ValueError(f"コア内の出口面が不足 ({int(core.sum())}, r_core={rc:.4g})")
    w = ro[core] * Ux[core] * A[core]
    eM = float(np.sqrt(np.average((M[core] - Md) ** 2, weights=w)) / Md)
    out = {"x_plane": x_plane, "r_core": rc, "core_def": how, "r_wall": r_wall,
           "eps_M_rms": eM,
           "eps_M_max": float(np.max(np.abs(M[core] - Md)) / Md),
           "eps_theta_max_deg": float(np.max(np.abs(th_deg[core]))),
           "M_core_massflux_avg": float(np.average(M[core], weights=w)),
           "n_core_faces": int(core.sum()), "n_exit_faces": int(len(rf))}
    if n_grid:   # 交差検証: 固定 r 格子へ補間して同じ量を積分し直す
        rg = np.linspace(0.0, rc, int(n_grid))
        Mg = np.interp(rg, rf, M)
        wg = np.interp(rg, rf, ro * Ux) * 2.0 * np.pi * rg
        out["eps_M_rms_grid"] = float(np.sqrt(np.average((Mg - Md) ** 2,
                                                         weights=np.maximum(wg, 0))) / Md)
    return out


def axis_mach(mesh_h5, res_h5, axis_band) -> tuple:
    r"""軸上/軸近傍の (x, M) を x 昇順で返す。**node では真の軸値**を返す。

    **cell/node で座標の意味が変わる罠** (2026-08-17 実測):
    node モードでは `/CELLS/centCoords` が**双対 CV 重心**に書き換わり VALUE と
    同じ長さになる。したがって cell 用の実装がエラーも出さずに「動いてしまう」が、
    軸上ノード (`MESH/COORD` で $r=0$) の双対 CV 重心は $r>0$ (実測 中央値
    0.022 $r_t$、最大 0.029) なので、`cy < axis_band` の帯選択には**非軸ノードが
    混入する** (実測: band=0.09 で 567 DOF 中 246 個が非軸、$M$ が最大 +0.001
    ずれる)。node は軸上に DOF を持つのが利点なので、それを使う。

    - node (`/MESH/COORD` の長さ = VALUE 長): $r=0$ の**軸ノードのみ**を返す。
      `axis_band` は無視される (真の軸値なので帯平均が不要)。
    - cell: 従来どおり `/CELLS/centCoords` の $r<$`axis_band` 帯。cell には軸上に
      DOF が無いため帯平均が避けられない (真の軸値が要る用途では
      `feedback.cfd_anchor.axis_curve_true` の偶関数フィットを使う)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        node_c = (nz["/MESH/COORD"][:].reshape(-1, 3) if "/MESH/COORD" in nz else None)
    with h5py.File(res_h5, "r") as f:
        Ux = f["/VALUE/Ux"][:]
        Uy = f["/VALUE/Uy"][:]
        son = f["/VALUE/sonic"][:]
    is_node = node_c is not None and len(node_c) == len(Ux)
    if is_node:
        m = node_c[:, 1] < 1e-12 * max(float(np.max(node_c[:, 1])), 1.0)
        x = node_c[m, 0]
    else:
        m = cc[:, 1] < axis_band
        x = cc[m, 0]
    M = np.hypot(Ux[m], Uy[m]) / np.maximum(son[m], 1e-9)
    o = np.argsort(x)
    return x[o], M[o]
