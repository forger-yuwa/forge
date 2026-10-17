#!/usr/bin/env python3
"""check_mesh_quality.py の壁直交性・AR の緩和の条件・結果の再利用の試験
(plan tooling-generator-defaults-and-launch-gate §5.1 #8・#8c、§6.2「壁直交性」・§6.3「再利用の鍵の試験」)。

#8c (codex diagnose 2026-10-14 `notes/reviews/2026-10-14-wall-orthogonality-judgments-diagnose.md`) の反例と事前登録の受入:
判別 A/B (直線壁の 5 セルを単一の physID と 2 分割で判定し、両方とも緩和の確認 2・警告 1・WARN)、AR 100・φ 40° の格子を
1 ID と 20 ID に分けても警告 20、真の角部を WARN から外さない、AR 6000.3 の 1 セルが --bcond で FAIL。

メッシュは小さい合成の h5 を直接書く (変換器は使わない)。node 変換と同じ配置: primal は `/VIZMESH/CONNE`
(XDMF の 5 = 四角形・4 = 三角形)、壁の面は `/BCONDS/<physID>/vizBfaceNodes`・`vizBfaceSizes`、座標は float32 の
`/MESH/COORD` (標準の変換器と同じ精度)。cell 変換の経路 (`iPlanes` → `/PLANES/STRUCT`) も 1 本試す。

合成の格子 (軸対称ノズルの縮流部の形、mesh2d と同じ作り方):
  壁 r_w(x) = 2 − x·tan40° (r_t 単位、x ∈ [0, 0.5/tan40°])、軸 r = 0、j 方向の線は半径方向 (x 一定)、壁側に等比で寄せる。
  `normal_layer=True` は mesh2d の wall_normal_layer と同じ式 (d ≤ d_n = 0.02 r_w を壁法線上に置き、d_b = 0.8 r_w まで
  smoothstep で半径方向へ戻す; forge-integ-1005 の mesh2d.py:310)。第一セルの AR は x の中央で ar_first。
"""
import json
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..', '..', 'tools', 'check_mesh_quality.py')
SCALE = 0.0767        # r_t [m] (case/45 と同じ桁; 座標の float32 の丸めを生産と同程度にする)

fails = 0


def check(cond, msg):
    global fails
    if not cond:
        fails += 1
        print('  FAIL:', msg)


# ---------------------------------------------------------------------------------------------------------------------
# 合成の格子
# ---------------------------------------------------------------------------------------------------------------------
def _radial_fracs(nj, f1):
    """s_0 = 0 (軸) … s_{nj−1} = 1 (壁)。壁側の第一の幅 f1、等比で軸へ広げる。"""
    lo, hi = 1.0 + 1e-12, 3.0
    for _ in range(200):                      # Σ_{k<nj−1} f1 q^k = 1 を q について二分法で解く
        q = 0.5 * (lo + hi)
        tot = f1 * (q ** (nj - 1) - 1.0) / (q - 1.0)
        lo, hi = (q, hi) if tot < 1.0 else (lo, q)
    gaps = f1 * q ** np.arange(nj - 1)
    gaps /= gaps.sum()
    s = 1.0 - np.r_[0.0, np.cumsum(gaps)]      # 壁 (1) から軸 (0) へ
    s[-1] = 0.0
    return s[::-1].copy()


def nozzle_mesh(theta_deg=40.0, ar_first=4000.0, normal_layer=False, ni=41, nj=61, dn=0.02, db=0.8, x_shift=0.0):
    """戻り (coords (N,3) [m], cells [list of node lists], 境界 {名前: (physID, [辺])})。"""
    th = np.radians(theta_deg)
    L = 0.5 / np.tan(th)
    xs = np.linspace(0.0, L, ni)
    rw = 2.0 - xs * np.tan(th)
    dx = xs[1] - xs[0]
    f1 = dx / (np.cos(th) * ar_first * 1.75)                  # 中央 (r_w = 1.75) で AR ≈ ar_first
    s = _radial_fracs(nj, f1)
    X = np.repeat(xs[:, None], nj, axis=1)
    R = rw[:, None] * s[None, :]
    if normal_layer:
        drw = np.full_like(xs, -np.tan(th))
        nrm = np.sqrt(1.0 + drw ** 2)
        d = rw[:, None] - R
        t = np.clip((d / rw[:, None] - dn) / (db - dn), 0.0, 1.0)
        beta = 1.0 - t * t * (3.0 - 2.0 * t)
        Xn = xs[:, None] + d * (drw / nrm)[:, None]
        Rn = rw[:, None] - d / nrm[:, None]
        X = X + beta * (Xn - X)
        R = R + beta * (Rn - R)
        assert np.all(np.diff(X, axis=0) > 0) and np.all(np.diff(R, axis=1) > 0), '格子の線が交差する'
    coords = np.zeros((ni * nj, 3))
    coords[:, 0] = X.ravel() * SCALE + x_shift          # x_shift [m]: 座標を大きくして float32 の丸めを効かせる
    coords[:, 1] = R.ravel() * SCALE
    nid = lambda i, j: i * nj + j                                  # noqa: E731
    cells = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    bnd = {'inlet': (1, [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)]),
           'outlet': (2, [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)]),
           'wall': (3, [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)]),
           'axis': (4, [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)])}
    return coords, cells, bnd


def split_wall_quads(cells, nj, ni_split):
    """壁に接する四角形のうち、x の先頭 ni_split 個を三角形 2 つに割る (四角形と三角形の混在)。"""
    out = []
    for k, c in enumerate(cells):
        i, j = divmod(k, nj - 1)
        if j == nj - 2 and i < ni_split:
            out += [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]
        else:
            out.append(c)
    return out


def strip_mesh(theta_deg=40.0, ar=100.0, ni=21, nj=6):
    """§6.3 の固定メッシュ: 下辺・上辺とも 40° 傾斜、線は鉛直、層厚一定 → 平行四辺形 (AR ar・skew 0.444・φ 40°)。"""
    th = np.radians(theta_deg)
    dx = 1.0e-3
    h = dx / (np.cos(th) * ar)
    nid = lambda i, j: i * nj + j                                  # noqa: E731
    coords = np.zeros((ni * nj, 3))
    for i in range(ni):
        for j in range(nj):
            coords[nid(i, j)] = (i * dx, i * dx * np.tan(th) + j * h, 0.0)
    cells = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    bnd = {'left': (1, [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)]),
           'right': (2, [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)]),
           'bottom': (3, [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)]),
           'top': (4, [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)])}
    return coords, cells, bnd


def extrude(coords, cells, bnd, nk=3, dz=None):
    """2D の四角形の格子を z 方向に nk 面で押し出して六面体の 3D 格子にする (z 平面が 3 枚以上 → 3D モード)。
    境界の辺は四角形の面に、前後 (z の両端) は physID 5・6 の面にする。"""
    n2 = len(coords)
    if dz is None:
        e = np.asarray(coords)[[c[1] for c in cells]] - np.asarray(coords)[[c[0] for c in cells]]
        dz = float(np.median(np.linalg.norm(e, axis=1)))
    co = np.concatenate([np.c_[np.asarray(coords)[:, :2], np.full(n2, k * dz)] for k in range(nk)])
    hexes = [[c[0] + k * n2, c[1] + k * n2, c[2] + k * n2, c[3] + k * n2,
              c[0] + (k + 1) * n2, c[1] + (k + 1) * n2, c[2] + (k + 1) * n2, c[3] + (k + 1) * n2]
             for k in range(nk - 1) for c in cells]
    b3 = {nm: (pid, [(a + k * n2, b + k * n2, b + (k + 1) * n2, a + (k + 1) * n2) for k in range(nk - 1) for a, b in edges])
          for nm, (pid, edges) in bnd.items()}
    b3['front'] = (5, [tuple(c) for c in cells])
    b3['back'] = (6, [tuple(x + (nk - 1) * n2 for x in c) for c in cells])
    return co, hexes, b3


def lcorner_mesh(ni=11, nj=11, h=1.0e-3):
    """直交の正方形格子。下辺 (physID 3) と左辺 (physID 1) を壁にすると原点が角部 (2 つの壁 physID が会う)。"""
    nid = lambda i, j: i * nj + j                                  # noqa: E731
    coords = np.array([(i * h, j * h, 0.0) for i in range(ni) for j in range(nj)])
    cells = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    bnd = {'left': (1, [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)]),
           'right': (2, [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)]),
           'bottom': (3, [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)]),
           'top': (4, [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)])}
    return coords, cells, bnd


def band_mesh():
    """codex の判別 A/B (diagnose 2026-10-14) の 5 セル: 直線壁 (y = 0) の上に 1 層、x の幅 1・1・0.1・0.1・0.1、層厚 0.00025
    → AR 4000・4000・400・400・400.8。最後のセルだけ外側の右の節点を h·tan40° ずらし、その横断辺の φ を 40° にする。
    戻り (coords, cells, 壁の辺 [5 本], 壁以外の境界 {名前: (physID, [辺])})。"""
    xs = np.array([0.0, 1.0, 2.0, 2.1, 2.2, 2.3])
    h = 0.00025
    shift = np.zeros(len(xs))
    shift[-1] = h * np.tan(np.radians(40.0))
    coords = np.array([(x + (shift[i] if j else 0.0), j * h, 0.0) for i, x in enumerate(xs) for j in range(2)])
    n = len(xs) - 1
    cells = [[2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1] for i in range(n)]
    wall = [(2 * i, 2 * i + 2) for i in range(n)]
    bnd = {'left': (1, [(0, 1)]), 'right': (2, [(2 * n, 2 * n + 1)]), 'top': (4, [(2 * i + 1, 2 * i + 3) for i in range(n)])}
    return coords, cells, wall, bnd


def band_kinds(wall_parts):
    """band_mesh の bcond: 壁以外は固定、壁は wall_parts {名前: (physID, 辺)} の分け方で kind wall。"""
    out = {'left': (1, 'inlet_Pressure'), 'right': (2, 'outlet_statPress'), 'top': (4, 'slip')}
    out.update({nm: (pid, 'wall') for nm, (pid, _) in wall_parts.items()})
    return out


STRIP_WALL_BOTTOM = {'left': (1, 'inlet_Pressure'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'wall'), 'top': (4, 'slip')}
STRIP_CORNER_KINDS = {'left': (1, 'wall'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'wall'), 'top': (4, 'slip')}


def strip_split(n_ids):
    """§6.3 の固定メッシュ (傾斜 40°・AR 100・φ 40°) の下辺の壁 20 面を n_ids 個の physID (1 または 20) に分ける。
    戻り (coords, cells, 境界, bcond の {名前: (physID, kind)})。"""
    co, ce, bn = strip_mesh()
    edges = bn.pop('bottom')[1]
    kinds = {k: v for k, v in STRIP_WALL_BOTTOM.items() if k != 'bottom'}
    if n_ids == 1:
        bn['bottom'] = (3, edges)
        kinds['bottom'] = (3, 'wall')
    else:
        assert len(edges) == n_ids
        for i, e in enumerate(edges):
            bn['bottom_%02d' % i] = (10 + i, [e])
            kinds['bottom_%02d' % i] = (10 + i, 'wall')
    return co, ce, bn, kinds


def ar6000_mesh():
    """codex の追加指摘の入力: 平らな壁の上の 1 層 2000 セル (AR 100)、最後のセルだけ外側の右の節点を 1/60 に下げて AR 6000.3。
    座標は float64 で書く (float32 だと 1.7e-7 の辺が丸めで消える)。"""
    co, ce, bn = strip_mesh(theta_deg=0.0, ar=100.0, ni=2001, nj=2)
    co = np.array(co)
    co[-1, 1] /= 60.0
    return co, ce, bn


def inproc_wall(coords, cells, walls, ar_max=5000.0):
    """check_mesh_quality.wall_orthogonality をプロセス内で回す (四角形だけ・2D)。walls: {physID: (名前, kind, [辺])}。
    戻り (結果の dict, ars)。セルごとの判定 (masks) を physID の分け方の間で比べるのに使う。"""
    sys.path.insert(0, os.path.dirname(TOOL))
    import argparse
    import check_mesh_quality as cmq                               # noqa: E402
    co = np.asarray(coords, dtype=np.float64)
    conn = np.asarray(cells, dtype=np.int64).ravel()
    offs = np.arange(1, len(cells) + 1, dtype=np.int64) * 4
    ty = np.full(len(cells), 9)
    ars, _ = cmq.metrics_vectorized(co[:, :2], conn, offs, ty, False)
    a = argparse.Namespace(ar_max=ar_max, phi_construct=cmq.PHI_CONSTRUCT_DEFAULT, theta_tol=cmq.THETA_TOL_DEFAULT,
                           corner_angle=cmq.CORNER_ANGLE_DEFAULT)
    w = {pid: (nm, kind, [np.asarray(e, dtype=np.int64) for e in faces]) for pid, (nm, kind, faces) in walls.items()}
    return cmq.wall_orthogonality(co, conn, offs, ty, w, ars, a, rnode=cmq.node_rounding_radius(co)), ars


def write_h5(path, coords, cells, bnd, layout='node', dtype=np.float32, kinds=None):
    """forge の入力 h5 のうち check_mesh_quality が読む部分だけを書く (BCONDS の bcondKind は kinds があれば書く)。"""
    kinds = kinds or {}
    conne = []
    for c in cells:
        conne += [{3: 4, 4: 5, 8: 9}[len(c)]] + list(c)          # XDMF: 三角形 4・四角形 5・六面体 9
    with h5py.File(path, 'w') as f:
        f.create_dataset('MESH/COORD', data=np.asarray(coords, dtype=dtype).ravel())
        if layout == 'node':
            g = f.create_group('VIZMESH')
            g.attrs['nVizCells'] = len(cells)
            g.attrs['vizCONNE_dim'] = len(conne)
            f.create_dataset('VIZMESH/CONNE', data=np.asarray(conne, dtype=np.int32))
            f.create_dataset('MESH/CONNE', data=np.zeros(1, dtype=np.int32))      # 双対 (この試験では読まない)
        else:
            f.create_dataset('MESH/CONNE', data=np.asarray(conne, dtype=np.int32))
            f.create_dataset('VALUE/ro', data=np.ones(len(cells), dtype=np.float32))
        st, ip = [], 0
        for _name, (pid, edges) in bnd.items():
            g = f.create_group('BCONDS/%d' % pid)
            g.attrs['bcondKind'] = kinds.get(pid, 'slip')
            n = len(edges)
            g.create_dataset('iCells', data=np.arange(n, dtype=np.int32))
            g.create_dataset('iBPlanes', data=np.arange(ip, ip + n, dtype=np.int32))
            g.create_dataset('iPlanes', data=np.arange(ip, ip + n, dtype=np.int32))
            if layout == 'node':
                g.create_dataset('vizBfaceNodes', data=np.asarray(edges, dtype=np.int32).ravel())
                g.create_dataset('vizBfaceSizes', data=np.asarray([len(e) for e in edges], dtype=np.int32))
            for e in edges:
                st += [len(e)] + [int(x) for x in e] + [2, 0, 0]
            ip += n
        if layout != 'node':
            f.create_dataset('PLANES/STRUCT', data=np.asarray(st, dtype=np.int32))


def write_bcond(path, kinds_by_name):
    """kinds_by_name: {名前: (physID, kind)}。"""
    with open(path, 'w') as fh:
        for nm, (pid, kind) in kinds_by_name.items():
            fh.write('%s: {physID: %d, kind: %s, outputHDFflg: 0, ints: , floats: }\n' % (nm, pid, kind))


def run(h5, *args):
    """check_mesh_quality を回す。戻り (終了コード, 出力, JSON の要約 or None)。"""
    js = h5 + '.summary.json'
    if os.path.exists(js):
        os.remove(js)
    extra = list(args)
    if '--bcond' in extra:
        extra += ['--json', js]
    p = subprocess.run([sys.executable, TOOL, h5] + extra, capture_output=True, text=True)
    summ = json.load(open(js)) if os.path.exists(js) else None
    return p.returncode, p.stdout + p.stderr, summ


def verdict_of(out):
    for ln in out.splitlines():
        if ln.startswith('VERDICT:'):
            return ln.split()[1]
    return None


NOZZLE_KINDS = {'inlet': (1, 'inlet_Pressure'), 'outlet': (2, 'outlet_statPress'), 'wall': (3, 'wall_isothermal'), 'axis': (4, 'axis')}
H5_KINDS = {pid: k for pid, k in NOZZLE_KINDS.values()}         # 変換時の bcondKind (bcond と同じ)


def main():
    tmp = tempfile.mkdtemp(prefix='cmq_wall_')
    bc = os.path.join(tmp, 'bcondConfig.yaml')
    write_bcond(bc, NOZZLE_KINDS)
    rows = []

    # ---- §6.2: 40° 傾斜壁の半径方向の格子 (AR 4000) は緩和を受けられず FAIL ----
    co, ce, bn = nozzle_mesh(ar_first=4000.0, normal_layer=False)
    h_rad = os.path.join(tmp, 'radial_ar4000.h5'); write_h5(h_rad, co, ce, bn, kinds=H5_KINDS)
    rc, out, sm = run(h_rad, '--ar-max', '5000', '--bcond', bc, '--no-cache')
    v = verdict_of(out)
    check(v == 'FAIL' and rc == 1, f'radial AR 4000 must FAIL with --bcond (got {v}, rc {rc})\n{out}')
    rel = (sm or {}).get('summary', {}).get('relax', {})
    check(rel.get('n_target', 0) > 0 and rel.get('n_verified', -1) == 0,
          f'radial AR 4000: every relaxation target must be unverified ({rel})')
    phimax = (sm or {}).get('summary', {}).get('phi_max')
    check(phimax is not None and abs(phimax - 40.0) < 0.5, f'radial: near-wall φ max must be ≈40° ({phimax})')
    rows.append(('radial 40° AR 4000 (--ar-max 5000)', v, rc, phimax, rel))

    # ---- §6.2: wall_normal_layer の格子 (AR 4000) は PASS ----
    co, ce, bn = nozzle_mesh(ar_first=4000.0, normal_layer=True)
    h_nrm = os.path.join(tmp, 'normal_ar4000.h5'); write_h5(h_nrm, co, ce, bn, kinds=H5_KINDS)
    rc, out, sm = run(h_nrm, '--ar-max', '5000', '--bcond', bc, '--no-cache')
    v = verdict_of(out)
    rel = (sm or {}).get('summary', {}).get('relax', {})
    phimax = (sm or {}).get('summary', {}).get('phi_max')
    check(v == 'PASS' and rc == 0, f'wall-normal layer AR 4000 must PASS (got {v}, rc {rc})\n{out}')
    check(phimax is not None and phimax < 5.0, f'wall-normal: near-wall φ (float32 storage) must stay below 5° ({phimax})')
    check(rel.get('n_target', 0) > 0 and rel.get('n_unverified', -1) == 0, f'wall-normal: all relaxation targets verified ({rel})')
    rows.append(('wall_normal_layer AR 4000 (--ar-max 5000)', v, rc, phimax, rel))

    # 同じ格子でも --ar-max を上げなければ従来どおり AR > 1000 は外れ値 (緩和の対象ではない)
    rc, out, sm = run(h_nrm, '--bcond', bc, '--no-cache')
    v = verdict_of(out)
    check(v == 'FAIL' and rc == 1, f'wall-normal AR 4000 with the default --ar-max 1000 must stay FAIL (got {v})')

    # ---- §6.2: AR 800 の半径方向の格子は WARN (φ の最大と位置) ----
    co, ce, bn = nozzle_mesh(ar_first=800.0, normal_layer=False)
    h_800 = os.path.join(tmp, 'radial_ar800.h5'); write_h5(h_800, co, ce, bn, kinds=H5_KINDS)
    rc, out, sm = run(h_800, '--bcond', bc, '--no-cache')
    v = verdict_of(out)
    wn = (sm or {}).get('summary', {}).get('warn_near_wall', {})
    check(v == 'WARN' and rc == 0, f'radial AR 800 must WARN with exit 0 (got {v}, rc {rc})\n{out}')
    check(wn.get('n', 0) > 0 and wn.get('phi_max') is not None and wn.get('at') is not None and 'WARN 近壁' in out,
          f'the WARN must give count, max φ and location ({wn})')
    rows.append(('radial 40° AR 800', v, rc, wn.get('phi_max'), {'warn_cells': wn.get('n'), 'at': wn.get('at')}))

    # ---- §6.2: 四角形と三角形の混在 → 三角形を対象外と出す ----
    # (AR の小さい格子で割る。AR 800 の四角形を割ると三角形の skew が 0.9 を超え、別の理由で FAIL になる)
    co, ce, bn = nozzle_mesh(ar_first=3.0, normal_layer=True)
    ce2 = split_wall_quads(ce, 61, 10)
    h_mix = os.path.join(tmp, 'mixed.h5'); write_h5(h_mix, co, ce2, bn, kinds=H5_KINDS)
    rc, out, sm = run(h_mix, '--bcond', bc, '--no-cache')
    w3 = [w for w in (sm or {}).get('summary', {}).get('walls', []) if w['physID'] == 3]
    check(bool(w3) and w3[0]['excluded_faces'] == 10 and w3[0]['excluded_types'] == {'三角形': 10} and '対象外 10 (三角形 10)' in out,
          f'mixed quad/tri: the 10 triangle wall faces must be reported as excluded ({w3})\n{out}')
    check(verdict_of(out) == 'PASS' and rc == 0, f'mixed quad/tri (wall-normal, AR <= 6) must PASS ({verdict_of(out)})')
    rows.append(('wall-normal AR 800 + 10 wall triangles', verdict_of(out), rc, (sm or {}).get('summary', {}).get('phi_max'),
                 {'excluded': w3[0]['excluded_types'] if w3 else None}))

    # ---- 3D (六面体、z 平面 3 枚 → 3D モード): 面は四角形、横断辺は 4 本 ----
    bc3 = os.path.join(tmp, 'bcond3d.yaml')
    write_bcond(bc3, dict(NOZZLE_KINDS, front=(5, 'slip'), back=(6, 'slip')))
    for nl, want, wrc in ((False, 'FAIL', 1), (True, 'PASS', 0)):
        co, ce, bn = nozzle_mesh(ar_first=4000.0, normal_layer=nl, ni=21, nj=41)
        co3, he, b3 = extrude(co, ce, bn, nk=3, dz=2.0e-3)
        h3 = os.path.join(tmp, 'hex_%s.h5' % ('normal' if nl else 'radial')); write_h5(h3, co3, he, b3, kinds=H5_KINDS)
        rc, out, sm = run(h3, '--ar-max', '5000', '--bcond', bc3, '--no-cache')
        s3 = (sm or {}).get('summary', {})
        check('mode: 3D' in out and verdict_of(out) == want and rc == wrc and s3.get('relax', {}).get('n_target', 0) > 0,
              f'3D hex ({"wall-normal" if nl else "radial"}): expected {want}/{wrc}, got {verdict_of(out)}/{rc}\n{out}')
        w3 = [w for w in s3.get('walls', []) if w['physID'] == 3]
        check(bool(w3) and w3[0]['edges_regular']['n'] == 4 * w3[0]['target_faces'], f'3D hex: 4 transverse edges per wall face ({w3})')
        rows.append(('3D hex %s AR 4000 (--ar-max 5000)' % ('wall-normal' if nl else 'radial'), verdict_of(out), rc, s3.get('phi_max'), s3.get('relax')))

    # ---- 角部: 2 つの壁 physID が会う節点に接するセルは分けて数える (直交の L 字なら φ 0・PASS) ----
    co, ce, bn = lcorner_mesh()
    hL = os.path.join(tmp, 'lcorner.h5'); write_h5(hL, co, ce, bn)
    bL = os.path.join(tmp, 'bcL.yaml')
    write_bcond(bL, {'left': (1, 'wall'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'wall_isothermal'), 'top': (4, 'slip')})
    rc, out, sm = run(hL, '--bcond', bL, '--no-cache')
    wl = {w['physID']: w for w in (sm or {}).get('summary', {}).get('walls', [])}
    check(verdict_of(out) == 'PASS' and set(wl) == {1, 3} and all('edges_corner' in w for w in wl.values())
          and all(w['edges_corner']['phi_max'] < 1e-6 for w in wl.values()) and '角部: 対象辺' in out,
          f'L corner: corner cells must be counted separately with φ 0 ({wl})\n{out}')
    cor = (sm or {}).get('summary', {}).get('corners', {})
    check(cor.get('n_nodes') == 1 and abs((cor.get('fold_max') or 0) - 90.0) < 1e-6 and cor.get('at') == [0.0, 0.0, 0.0]
          and cor.get('n_seam_nodes') == 0, f'L corner: one geometric corner (fold 90° at the origin), no label seam ({cor})')

    # ---- 格納精度: 座標が大きい (x + 1 m) 薄い層は float32 の丸めで φ が揺れる → 確かめられない (丸めの範囲) で FAIL、
    #      同じ格子を float64 で書けば PASS (判定は格納された形で行い、丸めの上限は理由の区別にだけ使う) ----
    co, ce, bn = nozzle_mesh(ar_first=4000.0, normal_layer=True, x_shift=1.0)
    for dt, want, wrc in ((np.float32, 'FAIL', 1), (np.float64, 'PASS', 0)):
        hs = os.path.join(tmp, 'normal_shift_%s.h5' % np.dtype(dt).name); write_h5(hs, co, ce, bn, dtype=dt, kinds=H5_KINDS)
        rc, out, sm = run(hs, '--ar-max', '5000', '--bcond', bc, '--no-cache')
        rel = (sm or {}).get('summary', {}).get('relax', {})
        ok = verdict_of(out) == want and rc == wrc
        if dt is np.float32:
            ok = ok and rel.get('reasons', {}).get('層の条件 (丸めの範囲)', 0) > 0 and '格納精度で確かめられない' in out
        check(ok, f'shifted wall-normal mesh in {np.dtype(dt).name}: expected {want}/{wrc} (float32: reason = rounding), '
                  f'got {verdict_of(out)}/{rc} {rel}\n{out}')
        rows.append(('wall-normal AR 4000, x+1 m, %s' % np.dtype(dt).name, verdict_of(out), rc,
                     (sm or {}).get('summary', {}).get('phi_max'), rel))

    # ---- cell 変換の経路 (iPlanes → PLANES/STRUCT) でも同じ判定 ----
    co, ce, bn = nozzle_mesh(ar_first=4000.0, normal_layer=False)
    h_cell = os.path.join(tmp, 'radial_ar4000_cell.h5'); write_h5(h_cell, co, ce, bn, layout='cell', kinds=H5_KINDS)
    rc, out, sm = run(h_cell, '--ar-max', '5000', '--bcond', bc, '--no-cache')
    check(verdict_of(out) == 'FAIL' and abs(((sm or {}).get('summary', {}).get('phi_max') or 0) - rows[0][3]) < 1e-6,
          f'cell layout must give the same φ and verdict as node layout\n{out}')

    # ---- 既存の挙動の回帰: --bcond なしは従来の規則 (AR・skew だけ、--ar-max 5000 は無条件) ----
    for h, args, want, wrc in ((h_rad, ['--ar-max', '5000'], 'PASS', 0), (h_nrm, ['--ar-max', '5000'], 'PASS', 0),
                               (h_rad, [], 'FAIL', 1), (h_800, [], 'PASS', 0)):
        rc, out, _ = run(h, *args)
        check(verdict_of(out) == want and rc == wrc and 'wall orthogonality' not in out and not os.path.exists(h + '.quality.json'),
              f'without --bcond {os.path.basename(h)} {args}: legacy verdict {want}/{wrc} expected, got {verdict_of(out)}/{rc}')
        # #8c (d): --bcond なしの判定は関門の合格の根拠にならない旨を必ず 1 行出す (VERDICT の行より前、語 VERDICT を含まない)
        note = [ln for ln in out.splitlines() if '投入の関門の合格の根拠にならない' in ln]
        check(len(note) == 1 and 'VERDICT' not in note[0] and out.splitlines().index(note[0]) < len(out.splitlines()) - 1,
              f'without --bcond {os.path.basename(h)} {args}: the "not evidence for the launch gate" note must appear once\n{out}')
    rc = subprocess.run([sys.executable, TOOL, h_rad, '--ar-max', '6000'], capture_output=True, text=True).returncode
    check(rc == 2, f'--ar-max above 5000 must be rejected (rc {rc})')

    # ---- §6.3: 再利用の鍵 (傾斜 40°・AR 100、試験用の閾値 φ 5°) ----
    co, ce, bn = strip_mesh()
    h63 = os.path.join(tmp, 'strip40_ar100.h5'); write_h5(h63, co, ce, bn)
    # 既存の cell_metrics() で AR 100・skew 0.444 (plan §6.3 の前提) を確かめる
    sys.path.insert(0, os.path.dirname(TOOL))
    import check_mesh_quality as cmq                               # noqa: E402
    ar0, sk0 = cmq.cell_metrics(np.asarray(co)[ce[0]][:, :2])
    check(abs(ar0 - 100.0) < 1e-6 and abs(sk0 - 0.444444) < 1e-5, f'strip cell must have AR 100 / skew 0.444 ({ar0}, {sk0})')
    cases = {'A': {'left': (1, 'inlet_Pressure'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'slip'), 'top': (4, 'slip')},
             'B': {'left': (1, 'inlet_Pressure'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'wall'), 'top': (4, 'slip')},
             'C': {'left': (1, 'inlet_Pressure'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'slip'), 'top': (4, 'wall')}}
    rec63 = {}
    for tag in ('A', 'B', 'C', 'B2'):
        b = os.path.join(tmp, 'bc63_%s.yaml' % tag)
        write_bcond(b, cases[tag[0]])
        rc, out, sm = run(h63, '--bcond', b, '--phi-construct', '5')
        s = (sm or {}).get('summary') or {}
        rec63[tag] = {'walls': s.get('n_walls'), 'phi_max': s.get('phi_max'), 'reused': (sm or {}).get('reused'),
                      'verdict': verdict_of(out), 'rc': rc}
    check(rec63['A'] == {'walls': 0, 'phi_max': None, 'reused': False, 'verdict': 'PASS', 'rc': 0}, f'§6.3 A: {rec63["A"]}')
    for tag in ('B', 'C'):
        r = rec63[tag]
        check(r['walls'] == 1 and r['phi_max'] is not None and abs(r['phi_max'] - 40.0) < 0.5 and r['reused'] is False
              and r['verdict'] == 'WARN' and r['rc'] == 0, f'§6.3 {tag}: must WARN with φ 40° and must not reuse the previous result ({r})')
    # 対照: 鍵が同じなら再利用する (B → C → B2 の B2 は C の結果を再利用しない、B2 の直後にもう一度回すと再利用)
    rc, out, sm = run(h63, '--bcond', os.path.join(tmp, 'bc63_B2.yaml'), '--phi-construct', '5')
    check((sm or {}).get('reused') is True and verdict_of(out) == 'WARN' and '結果の再利用: ' in out and '鍵が一致' in out,
          f'§6.3 control: an identical key must reuse the stored result ({sm})')
    rec63['B2-again'] = {'reused': (sm or {}).get('reused'), 'verdict': verdict_of(out)}
    # 閾値を変えたら再利用しない
    rc, out, sm = run(h63, '--bcond', os.path.join(tmp, 'bc63_B2.yaml'), '--phi-construct', '45')
    check((sm or {}).get('reused') is False and verdict_of(out) == 'PASS', f'§6.3: another threshold must not reuse ({sm})')
    # 保存量 (/VALUE) だけを書き換えても再利用する (メッシュの不変部分のハッシュに /VALUE を含めない)
    with h5py.File(h63, 'a') as f:
        f.create_dataset('VALUE/roe', data=np.ones(4, dtype=np.float32))
    rc, out, sm = run(h63, '--bcond', os.path.join(tmp, 'bc63_B2.yaml'), '--phi-construct', '45')
    check((sm or {}).get('reused') is True, f'§6.3: updating /VALUE only must keep the key ({sm})')
    # 座標を 1 節点だけ動かしたら再利用しない
    with h5py.File(h63, 'a') as f:
        c = np.array(f['MESH/COORD']); c[-1] += np.float32(1e-6); f['MESH/COORD'][...] = c
    rc, out, sm = run(h63, '--bcond', os.path.join(tmp, 'bc63_B2.yaml'), '--phi-construct', '45')
    check((sm or {}).get('reused') is False and 'mesh_invariant_sha256' in out, f'§6.3: moving a node must change the key ({sm})')

    # ==== #8c (codex diagnose 2026-10-14 の反例と事前登録の受入) ====
    rec8c = []

    def summ_of(sm):
        s = (sm or {}).get('summary') or {}
        rel, wn = s.get('relax', {}), s.get('warn_near_wall', {})
        return {'verified': rel.get('n_verified'), 'unverified': rel.get('n_unverified'), 'warn': wn.get('n'),
                'warn_corner': wn.get('n_corner'), 'corners': s.get('corners', {}).get('n_nodes'),
                'seams': s.get('corners', {}).get('n_seam_nodes')}

    # ---- 判別 A/B (事前登録): 直線壁の 5 セル (AR 4000・4000・400・400・400.8、下流の低 AR のセルに φ 40°)。
    #      A = 単一の physID、B = 上流 2 面と下流 3 面に分割、キャッシュなし。受入: 両方とも緩和の確認 2・警告 1・WARN
    #      (修正前は A = 確認 0・FAIL、B = 確認 2・WARN) ----
    co, ce, wall, bnd = band_mesh()
    ab = {}
    for tag, parts in (('A', {'bottom': (3, wall)}), ('B', {'up': (3, wall[:2]), 'down': (5, wall[2:])})):
        kinds = band_kinds(parts)
        hb = os.path.join(tmp, 'band_%s.h5' % tag)
        write_h5(hb, co, ce, dict(bnd, **parts), dtype=np.float64, kinds={pid: k for pid, k in kinds.values()})
        bb = os.path.join(tmp, 'band_%s.yaml' % tag); write_bcond(bb, kinds)
        rc, out, sm = run(hb, '--ar-max', '5000', '--bcond', bb, '--no-cache')
        r = dict(summ_of(sm), verdict=verdict_of(out), rc=rc, target=(sm or {}).get('summary', {}).get('relax', {}).get('n_target'))
        ab[tag] = r
        check(r['target'] == 2 and r['verified'] == 2 and r['unverified'] == 0 and r['warn'] == 1 and r['verdict'] == 'WARN'
              and r['rc'] == 0, f'#8c band {tag}: pre-registered acceptance is verified 2 / warn 1 / WARN (got {r})\n{out}')
        rec8c.append(('band 5 cells %s (%s)' % (tag, 'single physID' if tag == 'A' else 'split 2+3'), r))
    # セルごとにも同じ (緩和の可否・警告の対象が physID の分け方で変わらない)
    rA, _ = inproc_wall(co, ce, {3: ('bottom', 'wall', wall)})
    rB, _ = inproc_wall(co, ce, {3: ('up', 'wall', wall[:2]), 5: ('down', 'wall', wall[2:])})
    check(all(np.array_equal(rA['masks'][k], rB['masks'][k]) for k in ('verified', 'warn', 'layer'))
          and list(np.flatnonzero(rA['masks']['verified'])) == [0, 1] and list(np.flatnonzero(rA['masks']['warn'])) == [4],
          f'#8c band: per-cell verified / warn must not depend on the physID split '
          f'(A verified {np.flatnonzero(rA["masks"]["verified"])} warn {np.flatnonzero(rA["masks"]["warn"])}, '
          f'B verified {np.flatnonzero(rB["masks"]["verified"])} warn {np.flatnonzero(rB["masks"]["warn"])})')

    # ---- 境界のラベルの継ぎ目は角部ではない: AR 100・φ 40° の格子の下辺の壁を 1 ID と 20 ID に分けても警告は 20 ----
    warn_masks = {}
    for nids in (1, 20):
        co, ce, bn, kinds = strip_split(nids)
        hs = os.path.join(tmp, 'strip_%dids.h5' % nids)
        write_h5(hs, co, ce, bn, kinds={pid: k for pid, k in kinds.values()})
        bs = os.path.join(tmp, 'strip_%dids.yaml' % nids); write_bcond(bs, kinds)
        rc, out, sm = run(hs, '--bcond', bs, '--no-cache')
        r = dict(summ_of(sm), verdict=verdict_of(out), rc=rc)
        check(r['warn'] == 20 and r['warn_corner'] == 0 and r['corners'] == 0 and r['seams'] == (19 if nids == 20 else 0)
              and r['verdict'] == 'WARN' and r['rc'] == 0, f'#8c strip AR 100 φ 40° with {nids} wall IDs: warn 20 / WARN (got {r})\n{out}')
        rec8c.append(('strip AR 100 φ 40°, %d wall ID%s' % (nids, 's' if nids > 1 else ''), r))
        walls_in = {pid: (nm, k, [e for e in bn[nm][1]]) for nm, (pid, k) in kinds.items() if k == 'wall'}
        warn_masks[nids] = inproc_wall(co, ce, walls_in, ar_max=1000.0)[0]['masks']['warn']
    check(np.array_equal(warn_masks[1], warn_masks[20]) and int(warn_masks[1].sum()) == 20,
          '#8c strip: the warned cells must be the same set for 1 and 20 wall IDs')

    # ---- 真の角部も WARN から外さない: 40° 傾斜の下辺 + 左辺の壁 (内角 50° → 法線の折れ 130°) ----
    co, ce, bn = strip_mesh()
    hc = os.path.join(tmp, 'strip_corner.h5')
    write_h5(hc, co, ce, bn, kinds={pid: k for pid, k in STRIP_CORNER_KINDS.values()})
    bcc = os.path.join(tmp, 'strip_corner.yaml'); write_bcond(bcc, STRIP_CORNER_KINDS)
    rc, out, sm = run(hc, '--bcond', bcc, '--no-cache')
    r = dict(summ_of(sm), verdict=verdict_of(out), rc=rc)
    cor = (sm or {}).get('summary', {}).get('corners', {})
    wn = (sm or {}).get('summary', {}).get('warn_near_wall', {})
    check(r['warn'] == 24 and r['warn_corner'] == 1 and r['corners'] == 1 and abs((cor.get('fold_max') or 0) - 130.0) < 1e-3
          and cor.get('at') == [0.0, 0.0, 0.0] and wn.get('corner_phi_max') is not None and abs(wn['corner_phi_max'] - 40.0) < 0.5
          and r['verdict'] == 'WARN' and '折れ角 最大 130.0°' in out and 'うち角部 1 セル' in out,
          f'#8c true corner: 20 + 5 − 1 = 24 warned cells incl. the corner cell, fold 130° at the origin (got {r}, {cor})\n{out}')
    rec8c.append(('strip + left wall (fold 130°)', dict(r, fold=cor.get('fold_max'))))

    # ---- 3D の角部 (六面体の L 字、直交): 角の線の 3 節点が折れ 90°、PASS ----
    co, ce, bn = lcorner_mesh()
    co3, he, b3 = extrude(co, ce, bn, nk=3, dz=1.0e-3)
    kL3 = {'left': (1, 'wall'), 'right': (2, 'outlet_statPress'), 'bottom': (3, 'wall_isothermal'), 'top': (4, 'slip'),
           'front': (5, 'slip'), 'back': (6, 'slip')}
    hL3 = os.path.join(tmp, 'lcorner3d.h5'); write_h5(hL3, co3, he, b3, kinds={pid: k for pid, k in kL3.values()})
    bL3 = os.path.join(tmp, 'lcorner3d.yaml'); write_bcond(bL3, kL3)
    rc, out, sm = run(hL3, '--bcond', bL3, '--no-cache')
    cor = (sm or {}).get('summary', {}).get('corners', {})
    check('mode: 3D' in out and verdict_of(out) == 'PASS' and rc == 0 and cor.get('n_nodes') == 3
          and abs((cor.get('fold_max') or 0) - 90.0) < 1e-6 and cor.get('n_seam_nodes') == 0,
          f'#8c 3D L corner: 3 corner nodes with fold 90°, PASS (got {verdict_of(out)}/{rc} {cor})\n{out}')
    rec8c.append(('3D hex L corner', dict(summ_of(sm), verdict=verdict_of(out), rc=rc, fold=cor.get('fold_max'))))

    # ---- AR > 5000 が 1 セル (2000 セル中) は --bcond で件数によらず FAIL (修正前は SOFT-PASS・終了コード 0) ----
    co, ce, bn = ar6000_mesh()
    h6 = os.path.join(tmp, 'ar6000.h5')
    write_h5(h6, co, ce, bn, dtype=np.float64, kinds={pid: k for pid, k in STRIP_WALL_BOTTOM.values()})
    b6 = os.path.join(tmp, 'ar6000.yaml'); write_bcond(b6, STRIP_WALL_BOTTOM)
    for args in (['--ar-max', '5000'], []):
        rc, out, sm = run(h6, '--bcond', b6, '--no-cache', *args)
        s = (sm or {}).get('summary') or {}
        check(verdict_of(out) == 'FAIL' and rc != 0 and s.get('ar_bad') == 1 and s.get('ar_over_cap') == 1
              and 'SOFT-PASS を当てない' in out and abs(float(out.split('AR 最大 ')[1].split()[0]) - 6000.3) < 0.1,
              f'#8c one cell of AR 6000.3 among 2000 with --bcond {args}: FAIL / rc != 0 (got {verdict_of(out)}/{rc})\n{out}')
        rec8c.append(('AR 6000.3 x1 / 2000, --bcond %s' % ' '.join(args), {'verdict': verdict_of(out), 'rc': rc,
                                                                          'ar_bad': s.get('ar_bad')}))
    # --bcond なしは旧判定のまま (互換): SOFT-PASS・終了コード 0
    rc, out, _ = run(h6, '--ar-max', '5000')
    check(verdict_of(out) == 'SOFT-PASS' and rc == 0 and '投入の関門の合格の根拠にならない' in out,
          f'#8c AR 6000.3 without --bcond keeps the legacy SOFT-PASS / rc 0 and prints the note (got {verdict_of(out)}/{rc})')
    rec8c.append(('AR 6000.3 x1 / 2000, no --bcond (legacy)', {'verdict': verdict_of(out), 'rc': rc}))

    print('--- §6.2 壁直交性 (φ は格納された座標で測った近壁の最大)')
    for r in rows:
        print('  %-44s VERDICT %-5s rc %d  φ max %s  %s' % (r[0], r[1], r[2], 'None' if r[3] is None else '%.2f°' % r[3], r[4]))
    print('--- §6.3 再利用の鍵 (傾斜 40°・AR 100、試験用の閾値 5°)')
    for k, r in rec63.items():
        print('  %-9s %s' % (k, r))
    print('--- #8c 列ごとの全層・角部・AR の外れ値 (codex diagnose 2026-10-14 の反例)')
    for k, r in rec8c:
        print('  %-44s %s' % (k, r))
    print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
