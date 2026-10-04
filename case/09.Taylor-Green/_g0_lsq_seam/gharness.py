#!/usr/bin/env python3
"""G0/G1/G2/G2' ハーネスの共通部 (plan boundary-node-periodic-gradient-fix §5.1 #5a, §6)。

ソルバ本体は触らず、変換済み h5 と res_*.h5 だけから次を CPU で再現する。

- 周期 group (root = group 内の最小 index。`mesh.cpp` buildPeriodicNodeGroups と同じ union-find の結果)。
  partner は周期半割面の重心 + 並進量で相手 bcond の面を探して決める (`setPeriodicPartner` と同じ幾何照合)。
- 合併 stencil の LSQ (G2 の参照): 座標は **格納 float32 のノード座標** (node は常に `nodeValueAtNode`=1 で
  ccx ← MESH/COORD、`main.cpp:1139`)、差分は double。同値類 (隣接の root 一致 ∧ |Δx−Δx'| ≤ 1e-4 h_min) の
  重複数 α=1/count、M_r = Σ α w d dᵀ、b = Σ α w d Δφ、スペクトル打ち切り (thresh 1e-2) は M_r に 1 回
  (`calcGradient_d.cu` lsqPre_mergePeriodic / lsqPre_pinv と同式)。
- 係数の float32 焼き込み版 (G0 の float32 反例用): c_mj = f32(M⁺ α w d) を焼き、f32 で Σ c Δφ。
- node GG (G1): `ransTransport_d.cu:41-67` calc_scalar_gradient_face_d と同じ面値規則 (内部面 fx=0.5、
  非周期境界面は φ[ic0]、周期半割面は除外)、合併体積で除算 (k/ω は除算、受動種は f32(1/V) を乗算)、
  gather は和 → broadcast。(a) float32 演算の再現、(b) double。
"""
import os
import shutil
import subprocess

import h5py
import numpy as np
import yaml

EPS32 = float(np.finfo(np.float32).eps)       # 2^-23 = 1.19e-7 (plan §6 の ε_f32)
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RUN_CASE = os.path.join(REPO, "solver_density_cuda", "tools", "run_case.sh")
LSQ_THRESH = 1.0e-2                            # gradLSQDegenThresh (固定内部定数。起動ログ "thresh=1.00e-02")


# ----------------------------------------------------------------------------------------------------
# メッシュ
# ----------------------------------------------------------------------------------------------------
def parse_plane_cells(strct, n_planes):
    """PLANES/STRUCT ([nn, nodes..., nc, cells...]) から各面の (ic0, ic1)。ic1 が無い面は -1。"""
    s = np.asarray(strct)
    ic = np.full((n_planes, 2), -1, dtype=np.int64)
    p = 0
    for ip in range(n_planes):
        nn = s[p]; p += 1 + nn
        nc = s[p]; p += 1
        ic[ip, 0] = s[p]
        if nc >= 2:
            ic[ip, 1] = s[p + 1]
        p += nc
    return ic


class Mesh:
    def __init__(self, h5path, bcond_yaml):
        with h5py.File(h5path, "r") as f:
            a = f["MESH"].attrs
            self.nCells = int(a["nCells"]); self.nNormal = int(a["nNormalPlanes"]); self.nPlanes = int(a["nPlanes"])
            self.xyz = np.asarray(f["MESH/COORD"], dtype=np.float32).reshape(-1, 3)[:self.nCells]
            self.vol = np.asarray(f["CELLS/volume"], dtype=np.float32)[:self.nCells]
            self.sv = np.asarray(f["PLANES/surfVect"], dtype=np.float32).reshape(-1, 3)
            self.pcent = np.asarray(f["PLANES/centCoords"], dtype=np.float32).reshape(-1, 3)
            self.pc = parse_plane_cells(f["PLANES/STRUCT"][()], self.nPlanes)
            self.bconds = {}
            for bid in f["BCONDS"]:
                g = f["BCONDS/" + bid]
                self.bconds[int(bid)] = dict(iPlanes=np.asarray(g["iPlanes"], dtype=np.int64),
                                             kind=g.attrs["bcondKind"] if "bcondKind" in g.attrs else "")
        with open(bcond_yaml) as fp:
            self.bcfg = {v["physID"]: v for v in yaml.safe_load(fp).values() if isinstance(v, dict) and "physID" in v}
        self.periodic_ids = {p for p, v in self.bcfg.items() if v.get("kind") == "periodic"}
        self.root = self._periodic_roots()
        self.groups = {}
        for c in np.nonzero(self.root != np.arange(self.nCells))[0]:
            self.groups.setdefault(int(self.root[c]), [int(self.root[c])]).append(int(c))
        self.nmember = np.ones(self.nCells, dtype=np.int64)
        for r, mem in self.groups.items():
            self.nmember[mem] = len(mem)
        # 内部面 incidence (両向き)
        i0, i1 = self.pc[:self.nNormal, 0], self.pc[:self.nNormal, 1]
        self.inc_m = np.concatenate([i0, i1]); self.inc_j = np.concatenate([i1, i0])
        x = self.xyz.astype(np.float64)
        self.inc_d = x[self.inc_j] - x[self.inc_m]        # GPU: double(ccx[jc]) - double(cx) (float 同士の差は double で厳密)
        self.inc_plane = np.concatenate([np.arange(self.nNormal), np.arange(self.nNormal)])
        # 周期半割面 / 非周期境界面
        per_planes = np.concatenate([self.bconds[p]["iPlanes"] for p in self.bconds if p in self.periodic_ids]) \
            if self.periodic_ids else np.empty(0, np.int64)
        self.is_periodic_plane = np.zeros(self.nPlanes, dtype=bool); self.is_periodic_plane[per_planes] = True
        self.bnd_planes = np.array([ip for ip in range(self.nNormal, self.nPlanes) if not self.is_periodic_plane[ip]],
                                   dtype=np.int64)
        # 合併体積 (float32、index 順の和: buildPeriodicNodeGroups と同じ)
        gv = np.zeros(self.nCells, dtype=np.float32)
        for c in range(self.nCells):
            gv[self.root[c]] = np.float32(gv[self.root[c]] + self.vol[c])
        self.vmerged32 = gv[self.root]
        gv64 = np.zeros(self.nCells); np.add.at(gv64, self.root, self.vol.astype(np.float64))
        self.vmerged64 = gv64[self.root]
        self._alpha = None

    # --- 周期 group ---
    def _periodic_roots(self):
        from scipy.spatial import cKDTree
        n = self.nCells
        par = np.arange(n)

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a

        self.partner_resid = 0.0
        for pid in sorted(self.periodic_ids):
            v = self.bcfg[pid]
            if int(v.get("ints", {}).get("type", 0)) != 0:
                raise SystemExit("回転周期は対象外")
            q = int(v["ints"]["partnerBCID"])
            fl = v.get("floats", {})
            off = np.array([fl.get("dx", 0.0), fl.get("dy", 0.0), fl.get("dz", 0.0)], dtype=np.float64)
            pa, pb = self.bconds[pid]["iPlanes"], self.bconds[q]["iPlanes"]
            tree = cKDTree(self.pcent[pb].astype(np.float64))
            dist, k = tree.query(self.pcent[pa].astype(np.float64) + off)
            self.partner_resid = max(self.partner_resid, float(dist.max()))
            for a_node, b_node in zip(self.pc[pa, 0], self.pc[pb[k], 0]):
                ra, rb = find(a_node), find(b_node)
                if ra != rb:
                    if rb < ra:
                        ra, rb = rb, ra
                    par[rb] = ra
        return np.array([find(c) for c in range(n)], dtype=np.int64)

    # --- LSQ 同値類の重複数 α (lsqPre_mergePeriodic と同じ規則) ---
    def alpha(self):
        if self._alpha is not None:
            return self._alpha
        al = np.ones(len(self.inc_m))
        by_root = {}
        gm = np.nonzero(self.nmember[self.inc_m] > 1)[0]
        for e in gm:
            by_root.setdefault(int(self.root[self.inc_m[e]]), []).append(e)
        self.n_cls_merged = 0
        for r, es in by_root.items():
            # root を先頭 (lsqPre_mergePeriodic: mem = [root, members...] の順に incidence を走査)
            es.sort(key=lambda e: (self.inc_m[e] != r, self.inc_m[e], e))
            d = self.inc_d[es]
            L = np.linalg.norm(d, axis=1)
            tol = 1e-4 * L[L > 0].min()
            jr = self.root[self.inc_j[es]]
            cls_rep, cls_cnt, lab = [], [], []
            for k in range(len(es)):
                found = -1
                for c, (cj, cd) in enumerate(cls_rep):
                    if cj == jr[k] and np.linalg.norm(d[k] - cd) <= tol:
                        found = c; break
                if found < 0:
                    cls_rep.append((jr[k], d[k])); cls_cnt.append(0); found = len(cls_rep) - 1
                cls_cnt[found] += 1; lab.append(found)
            for k, e in enumerate(es):
                al[e] = 1.0 / cls_cnt[lab[k]]
            self.n_cls_merged += sum(1 for c in cls_cnt if c > 1)
        self._alpha = al
        return al

    def face_count(self):
        """GG で node (group) に積算される面の数 (内部面 incidence + 非周期境界面、全 partial 合計)。"""
        cnt = np.bincount(self.root[self.inc_m], minlength=self.nCells).astype(np.int64)
        cnt += np.bincount(self.root[self.pc[self.bnd_planes, 0]], minlength=self.nCells)
        return cnt[self.root]


# ----------------------------------------------------------------------------------------------------
# LSQ
# ----------------------------------------------------------------------------------------------------
def pinv_sym3(M, thresh=LSQ_THRESH):
    """lsqPre_pinv と同じ打ち切り規則 (固有値 l0<=l1<=l2、cut=thresh·l2)。戻り値 (Minv, degen)。"""
    lam, V = np.linalg.eigh(M)
    out = np.zeros_like(M); degen = np.zeros(len(M), dtype=bool)
    for i in range(len(M)):
        l0, l1, l2 = lam[i]
        if l2 <= 1e-300:
            degen[i] = True; continue
        keep = [2] if l1 < thresh * l2 else ([1, 2] if l0 < thresh * l2 else [0, 1, 2])
        degen[i] = len(keep) < 3
        if not degen[i]:
            out[i] = np.linalg.inv(M[i])
        else:
            out[i] = sum(np.outer(V[i][:, k], V[i][:, k]) / lam[i][k] for k in keep)
    return out, degen


def _lsq_system(msh, merged=True):
    """LSQ の (M⁺ を持つ節点 index の root 配列, M⁺, incidence の重み w)。merged=False は継ぎ目を合併しない
    (root = 自分、α = 1: 軸対称 × 周期・回転周期で GPU が使う片側 LSQ と同じ。plan gradient-scalar-lsq-unification S0-a)。"""
    d = msh.inc_d
    al = msh.alpha() if merged else np.ones(len(msh.inc_m))
    w = al / np.maximum(np.sum(d * d, axis=1), 1e-300)
    root = msh.root if merged else np.arange(msh.nCells)
    r = root[msh.inc_m]
    n = msh.nCells
    M = np.zeros((n, 3, 3))
    for a in range(3):
        for c in range(a, 3):
            M[:, a, c] = np.bincount(r, weights=w * d[:, a] * d[:, c], minlength=n)
            M[:, c, a] = M[:, a, c]
    roots = np.unique(root)
    Minv, degen = pinv_sym3(M[roots])
    Mi = np.zeros((n, 3, 3)); Mi[roots] = Minv
    msh.lsq_ndegen = int(degen.sum())
    return root, Mi, w


def lsq_merged_ref(msh, phi, merged=True):
    """合併 stencil LSQ の double 参照 (G2)。phi は格納 float32 の場 (nCells)。戻り値 (n,3)、group 全員同値。
    merged=False: 継ぎ目を合併しない片側 LSQ (節点ごと、α=1)。"""
    root, Mi, w = _lsq_system(msh, merged)
    d = msh.inc_d
    dphi = phi.astype(np.float64)[msh.inc_j] - phi.astype(np.float64)[msh.inc_m]
    r = root[msh.inc_m]
    n = msh.nCells
    b = np.zeros((n, 3))
    for a in range(3):
        b[:, a] = np.bincount(r, weights=w * d[:, a] * dphi, minlength=n)
    g = np.einsum("nij,nj->ni", Mi, b)
    return g[root]


def lsq_partials(msh, phi):
    """合併係数での member ごとの部分和 p_m = M⁺_root · Σ_{m の incidence} α w d Δφ (double)。GPU の周期 gather 前の
    局所配列に当たる。和 Σ_m p_m が lsq_merged_ref。戻り値 (n,3) (節点ごと = member ごと)。S0-c の Σ|部分和| 用。"""
    root, Mi, w = _lsq_system(msh, True)
    d = msh.inc_d
    dphi = phi.astype(np.float64)[msh.inc_j] - phi.astype(np.float64)[msh.inc_m]
    n = msh.nCells
    b = np.zeros((n, 3))
    for a in range(3):
        b[:, a] = np.bincount(msh.inc_m, weights=w * d[:, a] * dphi, minlength=n)
    return np.einsum("nij,nj->ni", Mi[root], b)


def bitcmp(a, b):
    """(不一致数 = ビット比較, 最大絶対差 (NaN を除く), 不一致位置の bool 配列)。形状が違えば (None, None, None)。"""
    if a.shape != b.shape:
        return None, None, None
    ai = a.view(np.uint32) if a.dtype == np.float32 else a.view(np.uint64)
    bi = b.view(np.uint32) if b.dtype == np.float32 else b.view(np.uint64)
    m = ai != bi
    dd = np.abs(a.astype(np.float64) - b.astype(np.float64))
    dd = dd[np.isfinite(dd)]
    return int(np.count_nonzero(m)), (float(dd.max()) if dd.size else 0.0), m


def noise_rule(A, B):
    """plan gradient-scalar-lsq-unification §6 S1 の判定規則 (2026-09-26 訂正版)。A, B: 同一設定反復の配列のリスト (各 3 本)。
    - A 同士の全対がビット一致 → A–B の全対もビット一致を要求。
    - それ以外: ノイズ対 = A 同士 + B 同士 (各 3 対)。A–B (9 対) の最大差 ≤ ノイズ対の最大差の 2 倍 かつ
      A–B の不一致数の最大 ≤ ノイズ対の不一致数の最大の 2 倍。
    戻り値 dict(verdict, aa, bb, ab, mask) (aa/bb/ab は (不一致数, 最大差) のリスト、mask は A–B の不一致位置の和集合)。"""
    import itertools
    aa = [bitcmp(x, y) for x, y in itertools.combinations(A, 2)]
    bb = [bitcmp(x, y) for x, y in itertools.combinations(B, 2)]
    ab = [bitcmp(x, y) for x in A for y in B]
    if any(t[0] is None for t in aa + bb + ab):
        return dict(verdict="判定不能 (形状不一致)", aa=[], bb=[], ab=[], mask=None)
    mask = np.zeros(A[0].shape, bool)
    for t in ab:
        mask |= t[2]
    n_ab, m_ab = max(t[0] for t in ab), max(t[1] for t in ab)
    if all(t[0] == 0 for t in aa):
        v = "PASS" if n_ab == 0 else "FAIL (旧同士ビット一致・旧新不一致)"
    else:
        noise = aa + bb
        n_no, m_no = max(t[0] for t in noise), max(t[1] for t in noise)
        ok_a = m_ab <= 2.0 * m_no
        ok_b = n_ab <= 2.0 * n_no
        v = "PASS" if (ok_a and ok_b) else f"FAIL ({'' if ok_a else '最大差 > 2×ノイズ '}{'' if ok_b else '不一致数 > 2×ノイズ'})".replace("( ", "(")
    return dict(verdict=v, aa=[t[:2] for t in aa], bb=[t[:2] for t in bb], ab=[t[:2] for t in ab], mask=mask)


def periodic_sin(w, ext, c0, amp, ph):
    """継ぎ目中心座標 w (n,3) の滑らかな周期関数 c0 + amp·[sin(2πw_x/L_x+φ0) + 0.7 sin(2πw_y/L_y+φ1)·cos(2πw_z/L_z+φ2)]。
    L = 各軸の領域幅 (周期軸では周期長なので折返しで不連続が出ない)。幅 0 の軸 (2D) は定数扱い。"""
    k = np.array([2 * np.pi / e if e > 0 else 0.0 for e in ext])
    return c0 + amp * (np.sin(k[0] * w[:, 0] + ph[0]) + 0.7 * np.sin(k[1] * w[:, 1] + ph[1]) * np.cos(k[2] * w[:, 2] + ph[2]))


# ----------------------------------------------------------------------------------------------------
# Green–Gauss (node、合併)
# ----------------------------------------------------------------------------------------------------
def gg_merged(msh, phi, precision="f32", div="divide", include_periodic=False):
    """node GG の合併勾配。precision: 'f32' (GPU の float32 演算を再現) / 'f64'。
    div: 'divide' (k/ω: dK /= V) / 'mulinv' (受動種・化学種: dY *= f32(1.0/V))。
    include_periodic: 周期半割面も境界面 (面値 φ[ic0]) として積算する (診断用。plan §4.2 は除外。
    実行時の plane_cells では周期半割面の ic1 が ghost index になる (`mesh.cpp:443-456`) ので、
    `excludePeriodic` の条件 `ic1 < nCells` が成り立たず GPU は積算している)。"""
    n = msh.nCells
    ii = np.arange(msh.nNormal)
    i0, i1 = msh.pc[ii, 0], msh.pc[ii, 1]
    b = np.arange(msh.nNormal, msh.nPlanes) if include_periodic else msh.bnd_planes
    if precision == "f64":
        p = phi.astype(np.float64)
        pf = np.concatenate([0.5 * p[i0] + 0.5 * p[i1], p[msh.pc[b, 0]]])
        S = msh.sv[np.concatenate([ii, b])].astype(np.float64)
        acc = np.zeros((n, 3))
        for a in range(3):
            acc[:, a] = np.bincount(np.concatenate([i0, msh.pc[b, 0]]), weights=S[:, a] * pf, minlength=n) \
                - np.bincount(i1, weights=S[:len(ii), a] * pf[:len(ii)], minlength=n)
        g = np.zeros((n, 3)); np.add.at(g, msh.root, acc)
        return g[msh.root] / msh.vmerged64[:, None]
    # float32 再現: 面値 f32(0.5a+0.5b)、項 f32(s·φf)、節点ごとに f32 で逐次加算 (順序は GPU と異なる → G1-a の 4ε)
    p = phi.astype(np.float32)
    pf_in = (0.5 * p[i0].astype(np.float64) + 0.5 * p[i1].astype(np.float64)).astype(np.float32)
    pf_b = p[msh.pc[b, 0]]
    acc = np.zeros((n, 3), dtype=np.float32)
    for a in range(3):
        s_in = msh.sv[ii, a]; s_b = msh.sv[b, a]
        t_in = (s_in.astype(np.float64) * pf_in.astype(np.float64)).astype(np.float32)
        t_b = (s_b.astype(np.float64) * pf_b.astype(np.float64)).astype(np.float32)
        col = np.zeros(n, dtype=np.float32)
        np.add.at(col, i0, t_in); np.add.at(col, i1, -t_in); np.add.at(col, msh.pc[b, 0], t_b)
        acc[:, a] = col
    V = msh.vmerged32
    if div == "divide":
        part = (acc.astype(np.float64) / V.astype(np.float64)[:, None]).astype(np.float32)
    else:
        invv = (1.0 / np.maximum(V.astype(np.float64), 1e-30)).astype(np.float32)
        part = (acc.astype(np.float64) * invv.astype(np.float64)[:, None]).astype(np.float32)
    g = part.copy()
    for c in np.nonzero(msh.root != np.arange(n))[0]:
        g[msh.root[c]] = (g[msh.root[c]].astype(np.float64) + part[c].astype(np.float64)).astype(np.float32)
    return g[msh.root]


# ----------------------------------------------------------------------------------------------------
# 場の焼き込み・run
# ----------------------------------------------------------------------------------------------------
def wrapped_coords(msh, cut_frac=0.5):
    """各軸で継ぎ目中心の局所座標 w (group 内で同値)。周期軸だけ折り返す: x > x0 + cut·(xL−x0) を x−xL に。
    x0/xL は格納 float32 の最小/最大 (差は double で厳密)。戻り値 (n,3) double と、折返し不連続に触れない節点のマスク。"""
    x = msh.xyz.astype(np.float64)
    w = x.copy()
    per_axes = set()
    for pid in msh.periodic_ids:
        fl = msh.bcfg[pid].get("floats", {})
        for a, k in enumerate(("dx", "dy", "dz")):
            if abs(fl.get(k, 0.0)) > 0:
                per_axes.add(a)
    for a in per_axes:
        x0, xL = x[:, a].min(), x[:, a].max()
        up = x[:, a] > x0 + cut_frac * (xL - x0)
        w[up, a] = x[up, a] - xL
        w[~up, a] = x[~up, a] - x0
    w = w[msh.root]                 # group 内で root の値に揃える (root が上側でも下側でも同じ 0)
    # 折返し不連続: incidence の Δw が Δx と一致しない節点を group ごと除外
    bad = np.any(np.abs((w[msh.inc_j] - w[msh.inc_m]) - msh.inc_d) > 1e-9, axis=1)
    badn = np.zeros(msh.nCells, dtype=bool); badn[msh.root[msh.inc_m[bad]]] = True
    return w, ~badn[msh.root], sorted(per_axes)


def bake(h5path, fields):
    """VALUE の保存量を上書き。fields: dict 名前 -> 配列 (float64 可)。"""
    with h5py.File(h5path, "r+") as f:
        v = f["VALUE"]
        for name, a in fields.items():
            if name in v:
                v[name][...] = np.asarray(a, dtype=np.float32)
            else:
                v.create_dataset(name, data=np.asarray(a, dtype=np.float32))


def prim_state_fields(msh, ro, ux, uy, uz, p, k=None, om=None, xi=None, gamma=1.4):
    ro = np.broadcast_to(np.asarray(ro, dtype=np.float64), (msh.nCells,))
    e = p / (gamma - 1.0) + 0.5 * ro * (ux ** 2 + uy ** 2 + uz ** 2)
    out = {"ro": ro, "roUx": ro * ux, "roUy": ro * uy, "roUz": ro * uz, "roe": e}
    if k is not None:
        out["roK"] = ro * k
    if om is not None:
        out["roOmega"] = ro * om
    if xi is not None:
        out["roXi"] = ro * xi
    return out


def make_run(run_dir, src_h5, solver_cfg, bcond_cfg, h5name="mesh.h5"):
    """スクラッチに run を作る (既存なら消して作り直す: スクラッチ専用)。solver_cfg は dict。"""
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)
    shutil.copy(src_h5, os.path.join(run_dir, h5name))
    cfg = dict(solver_cfg)
    cfg["mesh"] = dict(cfg.get("mesh", {}), meshFileName=h5name, valueFileName=h5name, discretization="node")
    with open(os.path.join(run_dir, "solverConfig.yaml"), "w") as fp:
        yaml.safe_dump(cfg, fp, sort_keys=False)
    if isinstance(bcond_cfg, str):
        shutil.copy(bcond_cfg, os.path.join(run_dir, "bcondConfig.yaml"))
    else:
        with open(os.path.join(run_dir, "bcondConfig.yaml"), "w") as fp:
            yaml.safe_dump(bcond_cfg, fp, sort_keys=False)
    with open(os.path.join(run_dir, "probe.yaml"), "w") as fp:
        fp.write("outStepInterval: 1\noutStepStart: 0\npoints:\nsurfaces:\n")
    return os.path.join(run_dir, h5name)


def run_forge(run_dir, expect_merge=True, forge_bin=None):
    """expect_merge=False: 周期はあるが継ぎ目の合併が無効な構成 (軸対称 × 周期) で、合併ログの有無を検査しない。
    forge_bin: FORGE_BIN (既定は run_case.sh の build/forge)。"""
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_CUDA_BLOCKSIZE_SMALL="128")
    env.pop("FORGE_DUMP_SCALARGRAD", None)
    if forge_bin:
        env["FORGE_BIN"] = forge_bin
    r = subprocess.run(["bash", RUN_CASE, run_dir], env=env, capture_output=True, text=True)
    log = open(os.path.join(run_dir, "forge_run.log")).read() if os.path.exists(os.path.join(run_dir, "forge_run.log")) else ""
    if not os.path.exists(os.path.join(run_dir, "res_1.h5")):
        print(log[-3000:])
        raise SystemExit(f"forge run failed (res_1.h5 無し): {run_dir} rc={r.returncode}")
    if expect_merge and "gradLSQ=2 periodic seam" not in log and "periodic" in open(os.path.join(run_dir, "bcondConfig.yaml")).read():
        print(log[-3000:])
        raise SystemExit(f"forge run failed or seam merge inactive: {run_dir} rc={r.returncode}")
    return log


def read_res(path, names):
    with h5py.File(path, "r") as f:
        return {k: np.asarray(f["VALUE/" + k], dtype=np.float32) for k in names if "VALUE/" + k in f}


def base_solver_cfg(sst=True, tracer=False, visc=0.0, nstep=1):
    cfg = {
        "gpu": 1, "solver": "SLAU",
        "physProp": {"thermalMethod": 0, "viscMethod": 0, "visc": visc, "thermCond": 0.0, "cp": 1038.8, "gamma": 1.4},
        "time": {"unsteady": 0, "dualTime": 0, "last": {"nStepOuter": nstep},
                 "deltaT": {"control": 1, "dt": 1e-9, "cfl": 0.01, "cfl_pseudo": 0.01, "blockDPLUR": 1,
                            "dt_min": 1e-12, "dt_max": 1e-7, "detectNaN": 1},
                 "outStepStart": 0, "outStepInterval": 1, "timeIntegration": 11, "nStepInner": 5},
        "space": {"convMethod": 0, "limiter": 0},
        "turbulence": {"model": "sst" if sst else "none"},
        "initial": "uniform_p101325_u10",
        "output": {"level": 2, "extraFields": ["dt_local", "volume"]},
    }
    if tracer:
        cfg["physProp"]["tracer"] = "exhaust"
        cfg["time"]["deltaT"]["speciesFaceReconstruction"] = 1
    return cfg


def git_rev():
    try:
        h = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", REPO, "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return h + (" (dirty)" if dirty else "")
    except Exception:
        return "unknown"


def provenance(run_dir):
    p = os.path.join(run_dir, "RUN_PROVENANCE.txt")
    if not os.path.exists(p):
        return ""
    return "".join(l for l in open(p) if l.startswith(("forge_sha256", "git_head", "forge_mtime")))


# ----------------------------------------------------------------------------------------------------
# FORGE_DUMP_PREGATHER (plan gradient-scalar-lsq-unification §5.1 #4a) の読み出しと、周期 gather の順列和
# ----------------------------------------------------------------------------------------------------
def read_pregather(path_tag):
    """`<path>.<tag>` と `.names` を読む。戻り値 dict 名前 -> (nCells, 3) float32。"""
    with open(path_tag + ".names") as fp:
        nv, nc = map(int, fp.readline().split())
        names = [l.strip() for l in fp if l.strip()]
    a = np.fromfile(path_tag, dtype=np.float32).reshape(nv, nc, 3)
    return {nm: a[i] for i, nm in enumerate(names)}


def gather_perm_sums(msh, pre, max_members=8):
    """周期 gather (`periodicGather1ToRoot_d`: root の値に member を atomicAdd) の結果になりうる float32 値の集合。
    pre: (nCells, 3) float32 の gather 前配列。戻り値 dict root -> 成分ごとの set (frozenset of float32 bit pattern)。
    root の値から始め、member を任意の順で float32 で足す (member 同士の順列 (n−1)! 通り)。"""
    import itertools
    out = {}
    for r, mem in msh.groups.items():
        others = [m for m in mem if m != r]
        if len(mem) > max_members:
            continue
        sets = []
        for c in range(3):
            vals = set()
            for perm in itertools.permutations(others):
                acc = np.float32(pre[r, c])
                for m in perm:
                    acc = np.float32(acc + np.float32(pre[m, c]))
                vals.add(acc.view(np.uint32).item())
            sets.append(vals)
        out[r] = sets
    return out
