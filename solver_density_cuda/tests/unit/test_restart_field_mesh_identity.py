#!/usr/bin/env python3
"""restart_field.py の同一メッシュの検査の試験 (plan tooling-generator-defaults-and-launch-gate §4.5・§5.1 #20)。

    python3 solver_density_cuda/tests/unit/test_restart_field_mesh_identity.py

合格条件 (§5.1 #20): 座標を格納の精度のまま比べる、接続・CV の順・境界の所属を照合する、座標・接続が無ければ FAIL。
codex の再現 (座標を約 1e−8 変えた別の入力; float32 に丸めると元と同じ) と接続だけ違う入力が FAIL、正常な同一メッシュの
restart (同じ run の res → その run の格子、親の res → 親の格子の複製、変換済みの入力 → 複製) が PASS。

変換器・forge は使わない。格子は変換済みの h5 と同じ置き方の小さな合成 (node は /MESH/CONNE が CV ごと・/VIZMESH/CONNE が
primal の四角形、cell は /MESH/CONNE が四角形; /PLANES・/CELLS・/BCONDS を持つ)、res は forge の出力と同じく
/MESH/COORD (float のビルドの型)・/MESH/CONNE (node は /VIZMESH/CONNE と同じ)・/VALUE だけを持つ (output.cpp:149・:191)。
宛先は CPG (thermalMethod 0) なので化学種の属性の照合は通る (ここでは見ない)。FAIL のときは DST のバイトが変わらないことも見る。
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.abspath(os.path.join(HERE, '..', '..', 'tools', 'restart_field.py'))
fails = 0
CONS = ('ro', 'roUx', 'roUy', 'roUz', 'roe')


def check(cond, msg, detail=''):
    global fails
    print(('  ok   ' if cond else '  FAIL ') + msg + (f'\n         {detail}' if detail and not cond else ''))
    if not cond:
        fails += 1


def grid(ni=7, nj=5):
    # float32 に丸めると値が変わる座標にする (丸めの照合が素通りしないように)
    xs = 0.1234567891 + np.linspace(0.0, 1.3, ni) * 1.0000001234
    ys = np.linspace(0.0, 0.7, nj) * 1.0000000123
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    coord = np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)]
    nid = lambda i, j: i * nj + j                                   # noqa: E731
    quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    wall = [nid(i, 0) for i in range(ni)]                           # 下辺 = 壁
    inlet = [nid(0, j) for j in range(nj)]                          # 左辺 = 入口
    return coord, quads, wall, inlet


def typed(quads):
    return np.array([x for q in quads for x in [5] + list(q)], dtype=np.int32)


def write_mesh(path, node=True, coord=None, quads=None, coord_dtype=np.float64, wall=None, inlet=None,
               volume_scale=1.0, wall_dist=None, drop=()):
    c0, q0, w0, i0 = grid()
    coord = c0 if coord is None else coord
    quads = q0 if quads is None else quads
    wall = w0 if wall is None else wall
    inlet = i0 if inlet is None else inlet
    nn = len(coord)
    ncv = nn if node else len(quads)
    with h5py.File(path, 'w') as f:
        g = f.create_group('MESH')
        g.attrs['nNodes'] = np.int32(nn)
        g.attrs['nCells'] = np.int32(ncv)
        f['MESH/COORD'] = np.asarray(coord, coord_dtype).ravel()
        if node:
            f['MESH/CONNE'] = np.array([x for n in range(nn) for x in (1, n)], dtype=np.int32)
            f['VIZMESH/CONNE'] = typed(quads)
            f['VIZMESH'].attrs['nVizCells'] = np.int32(len(quads))
        else:
            f['MESH/CONNE'] = typed(quads)
        f['PLANES/STRUCT'] = np.array([x for q in quads for x in (2, q[0], q[1])], dtype=np.int32)
        f['PLANES/surfArea'] = np.ones(len(quads), np.float32)
        f['CELLS/STRUCT'] = np.array([x for q in quads for x in (4, *q)], dtype=np.int32)
        f['CELLS/volume'] = np.full(ncv, 0.01 * volume_scale, np.float32)
        for pid, (kind, mem) in {1: ('wall', wall), 2: ('inlet_Pressure', inlet)}.items():
            b = f.create_group(f'BCONDS/{pid}')
            b['iCells'] = np.array(mem, dtype=np.int32)
            b.attrs['bcondKind'] = kind
            b['VALUE/Ps'] = np.full(len(mem), 1.0e5)
        f['VALUE/ro'] = np.ones(ncv, np.float32)
        for k in ('roUx', 'roUy', 'roUz'):
            f['VALUE/' + k] = np.zeros(ncv, np.float32)
        f['VALUE/roe'] = np.full(ncv, 2.5e5, np.float32)
        f['VALUE/wall_dist'] = (np.asarray(coord)[:, 1] if node else np.linspace(0, 1, ncv)) if wall_dist is None else wall_dist
        for k in drop:
            del f[k]


def write_res(path, mesh, dtype=np.float32, seed=1, coord_override=None, conne_override=None, drop=()):
    """forge の出力と同じ置き方の res (格子の構造は持たない)。保存量は宛先と違う値にする。"""
    rng = np.random.default_rng(seed)
    with h5py.File(mesh, 'r') as m, h5py.File(path, 'w') as f:
        n = m['VALUE/ro'].shape[0]
        f['MESH/COORD'] = (np.asarray(m['MESH/COORD']).astype(dtype) if coord_override is None else coord_override)
        ck = 'VIZMESH/CONNE' if 'VIZMESH/CONNE' in m else 'MESH/CONNE'
        f['MESH/CONNE'] = np.asarray(m[ck]) if conne_override is None else conne_override
        for k in CONS:
            f['VALUE/' + k] = (1.0 + rng.random(n)).astype(dtype)
        f['VALUE/P'] = np.full(n, 1.0e5, dtype)
        f['VALUE/wall_dist'] = np.asarray(m['VALUE/wall_dist']).astype(dtype) + 5.0
        for k in drop:
            del f[k]


def write_cfg(d, mesh='mesh.h5', value=None, node=True, text=None):
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as fh:
        fh.write(text if text is not None else
                 f'mesh: {{meshFileName: "{mesh}", valueFileName: "{value or mesh}", '
                 f'discretization: "{"node" if node else "cell"}"}}\n'
                 'gpu: 1\nsolver: "SLAU"\nphysProp: {thermalMethod: 0, cp: 1004.5, gamma: 1.4}\n')


class Env:
    def __init__(self):
        self.root = tempfile.mkdtemp(prefix='rf_ident_')
        self.n = 0

    def run_dir(self, node=True, cfg=True, **mesh_kw):
        self.n += 1
        d = os.path.join(self.root, f'run_{self.n:02d}')
        os.makedirs(d)
        if cfg:
            write_cfg(d, node=node)
        write_mesh(os.path.join(d, 'mesh.h5'), node=node, **mesh_kw)
        return d


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def tool(*args):
    env = dict(os.environ)
    env.pop('FORGE_ALLOW_UNVERIFIED_SPECIES', None)
    env['FORGE_BIN'] = '/nonexistent/forge'            # forge は起動しない (CPG の宛先では呼ばれない)
    r = subprocess.run([sys.executable, TOOL, *map(str, args)], capture_output=True, text=True, env=env)
    return r.returncode, r.stdout + r.stderr


def refused(dst, *args, want=()):
    """restart_field が同一メッシュと確かめられずに停止し、DST を書き換えないこと。"""
    before = sha(dst)
    rc, out = tool(*args)
    ok = (rc != 0 and 'REFUSED (nothing written)' in out and '同一メッシュと確かめられない' in out
          and all(w in out for w in want) and sha(dst) == before)
    return ok, f'rc={rc} {out.strip()[-500:]}'


def accepted(src, dst, *extra):
    rc, out = tool(src, dst, *extra)
    with h5py.File(src, 'r') as s, h5py.File(dst, 'r') as d:
        moved = all(np.array_equal(np.asarray(s['VALUE/' + k]).astype(d['VALUE/' + k].dtype), np.asarray(d['VALUE/' + k]))
                    for k in CONS)
    return rc == 0 and 'VERDICT: OK' in out and moved, out


def main():
    E = Env()
    try:
        print('=== 正例: 正常な同一メッシュの restart ===')
        # (P1) 同じ run の res (float32) → その run の格子 (倍精度の座標)。設計の runner の段間 (runner_axismach._restart_same_mesh)
        d = E.run_dir()
        mesh = os.path.join(d, 'mesh.h5')
        write_res(os.path.join(d, 'res_100.h5'), mesh)
        wd0 = np.asarray(h5py.File(mesh, 'r')['VALUE/wall_dist'])
        ok, out = accepted(os.path.join(d, 'res_100.h5'), mesh)
        wd1 = np.asarray(h5py.File(mesh, 'r')['VALUE/wall_dist'])
        check(ok and 'SRC と DST の格子は同じファイル' in out and 'float64 を float32 に丸めて' in out and np.array_equal(wd0, wd1),
              '(P1) 同じ run の float32 の res → その run の倍精度の格子: 通る、res は float32 で格子に結び付ける、wall_dist は宛先のまま',
              out[-600:])
        # (P2) 親の res → 親の格子の複製 (polish・rerun_conditions・case の継続)
        d2 = os.path.join(E.root, 'child')
        os.makedirs(d2)
        write_cfg(d2)
        shutil.copy(mesh, os.path.join(d2, 'mesh.h5'))
        ok, out = accepted(os.path.join(d, 'res_100.h5'), os.path.join(d2, 'mesh.h5'))
        check(ok and 'SRC と DST の格子が一致' in out, '(P2) 親の res → 親の格子の複製: 部分のハッシュが一致して通る', out[-600:])
        # (P3) cell の格子 (res の /MESH/CONNE = 格子の /MESH/CONNE)
        dc = E.run_dir(node=False)
        write_res(os.path.join(dc, 'res_100.h5'), os.path.join(dc, 'mesh.h5'))
        ok, out = accepted(os.path.join(dc, 'res_100.h5'), os.path.join(dc, 'mesh.h5'))
        check(ok and '接続 (/MESH/CONNE)' in out, '(P3) cell の格子: 通る (接続は /MESH/CONNE と照合)', out[-600:])
        # (P4) 変換済みの入力 h5 (格子を持つ) → その複製
        d3 = os.path.join(E.root, 'chain')
        os.makedirs(d3)
        write_cfg(d3)
        shutil.copy(mesh, os.path.join(d3, 'mesh.h5'))
        ok, out = accepted(os.path.join(d2, 'mesh.h5'), os.path.join(d3, 'mesh.h5'))
        check(ok and 'SRC の格子: ' in out and '(自身)' in out, '(P4) 変換済みの入力 h5 → 複製: 自身を格子として通る', out[-600:])
        # (P5) res だけを別の場所へ写した場合 → --src-mesh で格子を与えれば通る
        lone = os.path.join(E.root, 'lone')
        os.makedirs(lone)
        shutil.copy(os.path.join(d, 'res_100.h5'), lone)
        ok, out = accepted(os.path.join(lone, 'res_100.h5'), os.path.join(d2, 'mesh.h5'), '--src-mesh', mesh)
        check(ok and '(--src-mesh)' in out, '(P5) 隣に config の無い res + --src-mesh → 通る', out[-600:])
        # (P6) 倍精度のビルドの res (float64) と float32 の格子: 広げるのは厳密なのでそのまま一致。--keep-src-dtype の振る舞いは従来どおり
        d6 = E.run_dir(coord_dtype=np.float32)
        m6 = os.path.join(d6, 'mesh.h5')
        write_res(os.path.join(d6, 'res_100.h5'), m6, dtype=np.float64)
        shutil.copy(m6, os.path.join(d6, 'keep.h5'))
        rc, out = tool(os.path.join(d6, 'res_100.h5'), m6)
        with h5py.File(m6, 'r') as f:
            dt = f['VALUE/ro'].dtype
        check(rc == 0 and dt == np.float32 and '型の縮小' in out and '丸めあり' in out,
              '(P6) float64 の res → float32 の格子 (--keep-src-dtype なし): 通る、DST は float32 のまま・丸めを報告 (従来どおり)', out[-600:])
        shutil.copy(os.path.join(d6, 'keep.h5'), m6)
        ok, out = accepted(os.path.join(d6, 'res_100.h5'), m6, '--keep-src-dtype')
        with h5py.File(m6, 'r') as f:
            dt = f['VALUE/ro'].dtype
        check(ok and dt == np.float64 and 'SRC とビット一致' in out,
              '(P6) 同じく --keep-src-dtype: DST を float64 で作り直しビット一致 (従来どおり)', out[-600:])
        # (P7) wall_dist だけ違う格子 (壁距離の条件の違う変換) → 宛先の wall_dist を残すので通り、注記を出す
        d7 = E.run_dir(wall_dist=np.linspace(1.0, 2.0, 35))
        ok, out = accepted(os.path.join(d, 'res_100.h5'), os.path.join(d7, 'mesh.h5'))
        check(ok and 'wall_dist は SRC の格子と違う' in out, '(P7) wall_dist だけ違う → 通る (注記あり)', out[-600:])
        # (P8) DST が格子を持たない値のファイル (valueFileName が別): 宛先 run の meshFileName の格子に結び付けて通る
        d8 = os.path.join(E.root, 'valfile')
        os.makedirs(d8)
        write_cfg(d8, mesh='mesh.h5', value='val.h5')
        shutil.copy(mesh, os.path.join(d8, 'mesh.h5'))
        write_res(os.path.join(d8, 'val.h5'), mesh, seed=9)
        ok, out = accepted(os.path.join(d, 'res_100.h5'), os.path.join(d8, 'val.h5'))
        check(ok and 'DST の格子: ' in out and 'val.h5 の 座標' in out, '(P8) DST が値のファイル → 宛先 run の格子に結び付けて通る', out[-600:])

        print('=== 負例: 同一メッシュでない・確かめられない ===')
        c0, q0, w0, i0 = grid()
        # (N1) codex の再現: 座標を約 1e−8 (相対) 変えた別の入力。float32 に丸めると元と同じ成分だけ変える (旧実装は通した)
        pert = c0 * (1.0 + 1e-8)
        same32 = pert.astype(np.float32) == c0.astype(np.float32)
        pert = np.where(same32, pert, c0)
        nchg = int(np.count_nonzero(pert != c0))
        dn1 = E.run_dir(coord=pert)
        m1 = os.path.join(dn1, 'mesh.h5')
        check(nchg > 0 and np.array_equal(pert.astype(np.float32), c0.astype(np.float32)),
              f'(N1) 前提: 座標の {nchg} 成分を相対 1e-8 変え、float32 では元と同じ')
        ok, info = refused(m1, os.path.join(d, 'res_100.h5'), m1, want=('部分 coord', '最大 |Δ|'))
        check(ok, '(N1) codex の再現 (res の run の格子と座標が ~1e-8 違う宛先) → FAIL (coord)、DST 不変', info)
        ok, info = refused(m1, os.path.join(d, 'res_100.h5'), m1, '--dry-run', want=('部分 coord',))
        check(ok, '(N1) 同じく --dry-run でも停止', info)
        ok, info = refused(m1, mesh, m1, want=('部分 coord',))
        check(ok, '(N1) 変換済みの入力どうし (SRC が格子) でも FAIL', info)
        # (N2) 接続だけ違う: node は primal の四角形の節点の並び、cell はセルの順 (= CV の順) と接続
        q2 = [list(q) for q in q0]
        q2[0] = [q2[0][1], q2[0][2], q2[0][3], q2[0][0]]
        dn2 = E.run_dir(quads=q2)
        m2 = os.path.join(dn2, 'mesh.h5')
        ok, info = refused(m2, os.path.join(d, 'res_100.h5'), m2, want=('部分 topology', 'VIZMESH/CONNE'))
        check(ok, '(N2) node: 座標は同じで接続だけ違う格子 → FAIL (topology)', info)
        q3 = [q0[1], q0[0]] + q0[2:]
        dn3 = E.run_dir(node=False, quads=q3)
        m3 = os.path.join(dn3, 'mesh.h5')
        ok, info = refused(m3, os.path.join(dc, 'res_100.h5'), m3, want=('部分 topology', 'MESH/CONNE'))
        check(ok, '(N2) cell: 座標は同じでセルの順 (CV の順) と接続だけ違う格子 → FAIL (topology)', info)
        # (N3) CV の順: 節点の番号を入れ替えた同じ形の格子 (node は CV = 節点)
        perm = np.arange(len(c0))[::-1]
        inv = np.empty_like(perm)
        inv[perm] = np.arange(len(perm))
        dn4 = E.run_dir(coord=c0[perm], quads=[[int(inv[n]) for n in q] for q in q0], wall=[int(inv[n]) for n in w0],
                        inlet=[int(inv[n]) for n in i0])
        m4 = os.path.join(dn4, 'mesh.h5')
        ok, info = refused(m4, os.path.join(d, 'res_100.h5'), m4, want=('部分 coord',))
        check(ok, '(N3) 節点の番号を入れ替えた同じ形の格子 (CV の順が違う) → FAIL', info)
        # (N4) 境界の所属だけ違う
        dn5 = E.run_dir(wall=w0[:-1])
        m5 = os.path.join(dn5, 'mesh.h5')
        ok, info = refused(m5, os.path.join(d, 'res_100.h5'), m5, want=('部分 boundary', 'BCONDS/1/iCells'))
        check(ok, '(N4) 境界の所属だけ違う → FAIL (boundary)', info)
        # (N5) 座標の格納の型だけ違う (同じ値を float32 で格納) → 格納の精度のまま比べるので別の入力
        dn6 = E.run_dir(coord_dtype=np.float32)
        m6b = os.path.join(dn6, 'mesh.h5')
        ok, info = refused(m6b, os.path.join(d, 'res_100.h5'), m6b, want=('部分 coord', 'float64', 'float32'))
        check(ok, '(N5) 同じ座標を float32 で格納した格子 (res は float32 で両方に一致) → FAIL (coord の型)', info)
        # (N6) 幾何だけ違う (CV の体積)
        dn7 = E.run_dir(volume_scale=1.5)
        m7 = os.path.join(dn7, 'mesh.h5')
        ok, info = refused(m7, os.path.join(d, 'res_100.h5'), m7, want=('部分 geometry', 'CELLS/volume'))
        check(ok, '(N6) 幾何だけ違う → FAIL (geometry)', info)
        # (N7) 座標・接続が無い
        dn8 = E.run_dir(drop=('MESH/CONNE',))
        m8 = os.path.join(dn8, 'mesh.h5')
        ok, info = refused(m8, os.path.join(d, 'res_100.h5'), m8, want=('/MESH/CONNE が無い',))
        check(ok, '(N7) DST の格子に /MESH/CONNE が無い → FAIL', info)
        dn9 = E.run_dir(drop=('MESH/COORD',))
        m9 = os.path.join(dn9, 'mesh.h5')
        ok, info = refused(mesh, m9, mesh, want=('/MESH/COORD が無い',))
        check(ok, '(N7) SRC の格子に /MESH/COORD が無い → FAIL', info)
        for k in ('MESH/CONNE', 'MESH/COORD'):
            r = os.path.join(d, f'res_no_{k.split("/")[1]}.h5')
            write_res(r, mesh, drop=(k,))
            ok, info = refused(os.path.join(d2, 'mesh.h5'), r, os.path.join(d2, 'mesh.h5'), want=(f'/{k} が無い',))
            check(ok, f'(N7) res に /{k} が無い → FAIL', info)
        # (N8) res が別の格子の場 (結び付けの失敗): 座標が違う / 接続が違う
        r = os.path.join(d, 'res_wrong_coord.h5')
        write_res(r, mesh, coord_override=(c0 + [1e-3, 0, 0]).ravel().astype(np.float32))
        ok, info = refused(os.path.join(d2, 'mesh.h5'), r, os.path.join(d2, 'mesh.h5'), want=('座標が格子',))
        check(ok, '(N8) res の座標が自分の run の格子と違う → FAIL', info)
        r = os.path.join(d, 'res_wrong_conne.h5')
        write_res(r, mesh, conne_override=typed(q2))
        ok, info = refused(os.path.join(d2, 'mesh.h5'), r, os.path.join(d2, 'mesh.h5'), want=('接続 /MESH/CONNE',))
        check(ok, '(N8) res の接続が自分の run の格子と違う → FAIL', info)
        # (N9) 格子を決められない
        ok, info = refused(os.path.join(d2, 'mesh.h5'), os.path.join(lone, 'res_100.h5'), os.path.join(d2, 'mesh.h5'),
                           want=('格子を決められない', 'solverConfig.yaml が無い'))
        check(ok, '(N9) 隣に config の無い res (--src-mesh なし) → FAIL', info)
        for k, (tag, txt) in enumerate((('meshFileName なし', 'mesh: {discretization: "node"}\nphysProp: {thermalMethod: 0}\n'),
                                        ('重複キー', 'mesh: {meshFileName: "mesh.h5", meshFileName: "other.h5"}\nphysProp: {thermalMethod: 0}\n'),
                                        ('無いファイル', 'mesh: {meshFileName: "gone.h5", valueFileName: "gone.h5"}\nphysProp: {thermalMethod: 0}\n'))):
            dd = os.path.join(E.root, f'cfg_{k}')
            os.makedirs(dd)
            write_cfg(dd, text=txt)
            shutil.copy(os.path.join(d, 'res_100.h5'), dd)
            ok, info = refused(os.path.join(d2, 'mesh.h5'), os.path.join(dd, 'res_100.h5'), os.path.join(d2, 'mesh.h5'))
            check(ok, f'(N9) SRC の run の config が {tag} → FAIL', info)
        dd = os.path.join(E.root, 'cfg_res_as_mesh')
        os.makedirs(dd)
        write_cfg(dd, mesh='res_100.h5')
        shutil.copy(os.path.join(d, 'res_100.h5'), dd)
        ok, info = refused(os.path.join(d2, 'mesh.h5'), os.path.join(dd, 'res_100.h5'), os.path.join(d2, 'mesh.h5'),
                           want=('変換済みの格子でない',))
        check(ok, '(N9) meshFileName が格子の構造を持たない h5 → FAIL', info)
        ok, info = refused(os.path.join(d3, 'mesh.h5'), mesh, os.path.join(d3, 'mesh.h5'), '--src-mesh', mesh,
                           want=('--src-mesh は使わない',))
        check(ok, '(N9) 格子を持つ SRC に --src-mesh → FAIL', info)
        # (N10) --src-mesh で別の格子 (座標が float32 でも違う) を与えた → 結び付けで FAIL
        ok, info = refused(os.path.join(d2, 'mesh.h5'), os.path.join(lone, 'res_100.h5'), os.path.join(d2, 'mesh.h5'),
                           '--src-mesh', m4, want=('座標が格子',))
        check(ok, '(N10) --src-mesh に別の格子 → FAIL', info)
        # (N11) DST が値のファイルで、宛先 run の格子と座標が違う
        d11 = os.path.join(E.root, 'valfile_bad')
        os.makedirs(d11)
        write_cfg(d11, mesh='mesh.h5', value='val.h5')
        shutil.copy(m4, os.path.join(d11, 'mesh.h5'))
        write_res(os.path.join(d11, 'val.h5'), mesh, seed=9)
        ok, info = refused(os.path.join(d11, 'val.h5'), os.path.join(d, 'res_100.h5'), os.path.join(d11, 'val.h5'),
                           want=('DST', '座標が格子'))
        check(ok, '(N11) DST の値のファイルが宛先 run の格子の場でない → FAIL', info)
    finally:
        shutil.rmtree(E.root, ignore_errors=True)
    print('ALL PASS' if fails == 0 else f'{fails} FAIL')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
