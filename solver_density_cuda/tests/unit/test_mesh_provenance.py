#!/usr/bin/env python3
"""mesh_provenance.py (変換済みメッシュの来歴の記録と照合) の試験
(plan tooling-generator-defaults-and-launch-gate §4.3 の 2・§5.1 #7・§6.2「メッシュの実効条件の負例」)。

合格条件 (§5.1 #7): 記録なし・記録の取り違え (別メッシュの記録)・入口の分類だけの変更・`wallDistExtraPhysIDs` だけの変更・
変換時の `nodeInletCornerWall` 0 と run の 1 が FAIL、保存量だけを更新した正常な restart (`restart_field.py`) は PASS。

変換器のバイナリは使わない (CFD も回さない)。偽の変換器 (下の FAKE_CONVERTER) が作業ディレクトリの solverConfig.yaml・
bcondConfig.yaml を読み、本物と同じ置き方の h5 を書く。変換時の条件への依存は本物の該当箇所をまねる:
- node: primal のセルを /VIZMESH/CONNE に退避 (gmshReader.hpp:2659)、cell は /VIZMESH なし
- nodeInletCornerWall 1 (node): 入口の境界の壁節点を入口の所属から壁の所属へ移す (gmshReader.hpp:1674 の半割面の帰属)
- wall_dist: 壁 (wall・wall_isothermal) と wallDistExtraPhysIDs の境界節点からの最短距離 (calcWallDistance_kdtree.cpp:142)
- axisCentroidShift: CV の重心を動かす、renumber rcm: 節点の順を入れ替え /MESH/RENUMBER_PERM、node の 2D 軸対称: /PLANES/rSurfVect
- /BCONDS/<physID> の属性 bcondKind は bcondConfig の種類 (gmshReader.hpp:2748)

入力の .msh の代わりに JSON (ni・nj・座標のずらし) を渡す。偽の変換器の挙動は環境変数で壊せる (config を読まない・失敗する・
変換中に config を書き換える・書かずに 0 で終わる)。
"""
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, '..', '..', 'tools'))
sys.path.insert(0, TOOLS)
import mesh_provenance as mp  # noqa: E402

TOOL = os.path.join(TOOLS, 'mesh_provenance.py')
fails = 0


def check(cond, msg):
    global fails
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        fails += 1


FAKE_CONVERTER = r'''
import json, os, sys
import numpy as np, h5py, yaml

msh, out = sys.argv[1], sys.argv[2]
if os.environ.get('FAKE_CONV_RC'):
    print('*** Read Mesh ***'); sys.exit(int(os.environ['FAKE_CONV_RC']))
if os.environ.get('FAKE_CONV_NOWRITE'):
    print('*** Read Mesh ***'); sys.exit(0)        # 書かずに 0 で終わる (旧い h5 が残る)
cfg = yaml.safe_load(open('solverConfig.yaml')) or {}
bc = yaml.safe_load(open('bcondConfig.yaml')) or {}
mesh = cfg.get('mesh') or {}
if os.environ.get('FAKE_CONV_IGNORE_CFG'):
    mesh = {}                                   # config を読まない変換器 (cell・既定のまま)
if os.environ.get('FAKE_CONV_TOUCH'):
    with open('solverConfig.yaml', 'a') as fh:
        fh.write('# touched during conversion\n')
spec = json.load(open(msh))
ni, nj, L, H = spec['ni'], spec['nj'], spec.get('L', 1.0), spec.get('H', 0.2)
dx = spec.get('shift', 0.0)
node = mesh.get('discretization', 'cell') == 'node'
corner = int(mesh.get('nodeInletCornerWall', 0)) != 0
shiftc = int(mesh.get('axisCentroidShift', 1)) != 0
rcm = mesh.get('renumber', 'none') == 'rcm'
axi = int(mesh.get('isAxisymmetric', (cfg.get('physProp') or {}).get('isAxisymmetric', 0))) == 1
extra = set(int(x) for x in mesh.get('wallDistExtraPhysIDs', []) or [])
kinds = {int(v['physID']): v['kind'] for v in bc.values()}

xs = np.linspace(0.0, L, ni) + dx
ys = np.linspace(0.0, H, nj)
X, Y = np.meshgrid(xs, ys, indexing='ij')
coord = np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)]
nid = lambda i, j: i * nj + j
quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
# 境界の辺: 1 入口 (x=0)、2 出口 (x=L)、3 壁 (上, x < 0.6L)、5 slip の延長 (上, x >= 0.6L)、4 軸 (下)
edges = {1: [], 2: [], 3: [], 4: [], 5: []}
for j in range(nj - 1):
    edges[1].append((nid(0, j), nid(0, j + 1)))
    edges[2].append((nid(ni - 1, j), nid(ni - 1, j + 1)))
for i in range(ni - 1):
    edges[4].append((nid(i, 0), nid(i + 1, 0)))
    edges[3 if xs[i] - dx < 0.6 * L else 5].append((nid(i, nj - 1), nid(i + 1, nj - 1)))
for pid in edges:
    if pid not in kinds:
        print(f"Error: couldn't find physID {pid} in bcondConfig.yaml"); sys.exit(1)
perm = np.arange(ni * nj)
if rcm:
    perm = perm[::-1].copy()                    # new -> old
inv = np.empty_like(perm); inv[perm] = np.arange(len(perm))
coord = coord[perm]
quads = [[int(inv[n]) for n in q] for q in quads]
edges = {p: [(int(inv[a]), int(inv[b])) for a, b in e] for p, e in edges.items()}
walls = [p for p, k in kinds.items() if k in ('wall', 'wall_isothermal') and p in edges]
bnodes = {p: sorted({n for e in edges[p] for n in e}) for p in edges}
if node and corner:                             # 入口の壁節点の半割面を壁へ (gmshReader.hpp:1674)
    wn = {n for p in walls for n in bnodes[p]}
    for p in edges:
        if kinds[p].startswith('inlet_'):
            moved = [n for n in bnodes[p] if n in wn]
            bnodes[p] = [n for n in bnodes[p] if n not in wn]
            if walls:
                bnodes[walls[0]] = sorted(set(bnodes[walls[0]]) | set(moved))
wpts = [coord[n] for p in edges if (p in walls or p in extra) for n in bnodes[p]]
if wpts:
    W = np.array(wpts)
    wd = np.sqrt(((coord[:, None, :2] - W[None, :, :2]) ** 2).sum(-1)).min(1)
else:
    wd = np.zeros(len(coord))
cent = coord.copy()
if shiftc:
    cent[:, 0] += 1e-3 * (np.arange(len(coord)) % 3)
with h5py.File(out, 'w') as f:
    g = f.create_group('MESH')
    for k, v in (('nCells', len(coord) if node else len(quads)), ('nNodes', len(coord)), ('nBconds', len(edges))):
        g.attrs[k] = np.int64(v)
    f['MESH/COORD'] = coord.ravel()
    if rcm:
        f['MESH/RENUMBER_PERM'] = perm.astype(np.int64)
    conne = np.array([x for q in quads for x in [5] + q], dtype=np.int64)
    if node:
        f['MESH/CONNE'] = np.array([x for n in range(len(coord)) for x in (1, n)], dtype=np.int64)
        f['VIZMESH/CONNE'] = conne
        f['VIZMESH'].attrs['nVizCells'] = np.int64(len(quads))
    else:
        f['MESH/CONNE'] = conne
    planes = [x for p in sorted(bnodes) for n in bnodes[p] for x in (1, n, 1, n)]
    f['PLANES/STRUCT'] = np.array(planes, dtype=np.int64)
    f['PLANES/surfArea'] = np.ones(len(planes) // 4)
    if node and axi:
        f['PLANES/rSurfVect'] = np.ones(3 * (len(planes) // 4))
    f['CELLS/centCoords'] = cent.ravel()
    f['CELLS/volume'] = np.ones(len(coord))
    f['VALUE/ro'] = np.ones(len(coord), np.float32)
    for k in ('roUx', 'roUy', 'roUz'):
        f['VALUE/' + k] = np.zeros(len(coord), np.float32)
    f['VALUE/roe'] = np.full(len(coord), 2.5e5, np.float32)
    f['VALUE/wall_dist'] = wd
    for p in sorted(edges):
        b = f.create_group(f'BCONDS/{p}')
        b['iCells'] = np.array(bnodes[p], dtype=np.int64)
        b['iPlanes'] = np.arange(len(bnodes[p]), dtype=np.int64)
        b['iBPlanes'] = np.arange(len(bnodes[p]), dtype=np.int64)
        if node:
            b['vizBfaceNodes'] = np.array([n for e in edges[p] for n in e], dtype=np.int64)
            b['vizBfaceSizes'] = np.full(len(edges[p]), 2, dtype=np.int64)
        b.attrs['bcondKind'] = kinds[p]
        b['VALUE/Ps'] = np.full(len(bnodes[p]), 1.0e5)
print('*** Write Input HDF5 ***')
'''

BASE_BC = {'inlet': (1, 'inlet_Pressure'), 'outlet': (2, 'outlet_statPress'), 'wall': (3, 'wall'), 'axis': (4, 'axis'),
           'buffer': (5, 'slip')}


def write_cfg(d, mesh_extra='', name='mesh.h5', value=None, physprop=''):
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as fh:
        fh.write(f'mesh: {{meshFileName: "{name}", valueFileName: "{value or name}"{mesh_extra}}}\n'
                 f'gpu: 1\nsolver: "SLAU"\n{physprop}')


def write_bc(d, overrides=None):
    bcs = dict(BASE_BC)
    bcs.update(overrides or {})
    with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as fh:
        for nm, (pid, kind) in bcs.items():
            fh.write(f'{nm}: {{physID: {pid}, kind: {kind}, outputHDFflg: 0, ints: , floats: }}\n')


class Env:
    def __init__(self):
        self.root = tempfile.mkdtemp(prefix='mp_test_')
        self.conv_py = os.path.join(self.root, 'fake_convert.py')
        with open(self.conv_py, 'w') as fh:
            fh.write(FAKE_CONVERTER)
        self.conv = os.path.join(self.root, 'convertGmshToForge')
        with open(self.conv, 'w') as fh:
            fh.write(f'#!/bin/sh\nexec "{sys.executable}" "{self.conv_py}" "$@"\n')
        os.chmod(self.conv, os.stat(self.conv).st_mode | stat.S_IXUSR)
        self.n = 0

    def run_dir(self, mesh_extra=', discretization: "node"', bc=None, spec=None, physprop=''):
        self.n += 1
        d = os.path.join(self.root, f'run_{self.n:02d}')
        os.makedirs(d)
        write_cfg(d, mesh_extra, physprop=physprop)
        write_bc(d, bc)
        with open(os.path.join(d, 'mesh.json'), 'w') as fh:
            json.dump(spec or {'ni': 21, 'nj': 6}, fh)
        return d

    def convert(self, d, env=None):
        old = {k: os.environ.get(k) for k in (env or {})}
        os.environ.update(env or {})
        try:
            return mp.convert_and_record('mesh.json', 'mesh.h5', cwd=d, converter=self.conv)
        finally:
            for k, v in old.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

    def raw_convert(self, d):
        """記録を書かずに変換だけ (記録なしの負例・旧い記録の負例を作る)。"""
        r = subprocess.run([self.conv, 'mesh.json', 'mesh.h5'], cwd=d, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


def fail_keys(res):
    return [k for k, _ in res['fails']]


def cli(*args):
    r = subprocess.run([sys.executable, TOOL] + list(args), capture_output=True, text=True)
    out = r.stdout + r.stderr
    v = [ln for ln in out.splitlines() if ln.startswith('VERDICT:')]
    return r.returncode, (v[-1] if v else '(VERDICT なし)'), out


def main():
    E = Env()
    try:
        print('=== 正例: 変換して記録、同じ config で照合 ===')
        d0 = E.run_dir(', discretization: "node", nodeInletCornerWall: 1')
        path, rec = E.convert(d0)
        check(os.path.isfile(path) and path.endswith('mesh.h5.provenance.json'), f'記録 {os.path.basename(path)} を書いた')
        check(rec['conversion']['mesh.nodeInletCornerWall'] == 1 and rec['conversion']['mesh.wallDistExtraPhysIDs'] is None,
              f'変換時の要求値 {rec["conversion"]}')
        check(rec['converter']['sha256'] == mp._sha_file(E.conv) and rec['source']['sha256'] is not None,
              '変換器・入力の sha256 を記録')
        check([b['kind'] for b in rec['bconds']] == ['inlet_Pressure', 'outlet_statPress', 'wall', 'axis', 'slip'],
              '全境界の種類と physID を記録')
        rc, v, out = cli('check', d0)
        check(rc == 0 and v == 'VERDICT: PASS', f'同じ config で照合 → {v} (rc={rc})')

        print('=== 負例 1: 記録の無いメッシュ ===')
        d1 = E.run_dir()
        E.raw_convert(d1)
        res = mp.check_run(d1)
        check(fail_keys(res) == ['(provenance)'] and '記録の無いメッシュ' in res['fails'][0][1], f'FAIL {res["fails"]}')
        rc, v, _ = cli('check', d1)
        check(rc == 1 and v == 'VERDICT: FAIL', f'CLI → {v} (rc={rc})')

        print('=== 負例 2: 記録の取り違え (別メッシュの記録) ===')
        dA = E.run_dir(spec={'ni': 21, 'nj': 6})
        E.convert(dA)
        dB = E.run_dir(spec={'ni': 21, 'nj': 6, 'shift': 0.01})        # 同じ大きさ・同じ境界、座標だけ違う
        E.raw_convert(dB)
        shutil.copy(os.path.join(dA, 'mesh.h5.provenance.json'), os.path.join(dB, 'mesh.h5.provenance.json'))
        res = mp.check_run(dB)
        msg = dict(res['fails']).get('(provenance)', '')
        check(fail_keys(res) == ['(provenance)'] and 'coord' in msg and '別のメッシュ' in msg, f'FAIL {msg}')
        dC = E.run_dir(spec={'ni': 31, 'nj': 6})                         # 大きさも違う
        E.raw_convert(dC)
        shutil.copy(os.path.join(dA, 'mesh.h5.provenance.json'), os.path.join(dC, 'mesh.h5.provenance.json'))
        res = mp.check_run(dC)
        check(fail_keys(res) == ['(provenance)'], f'節点数の違う別メッシュ → FAIL {fail_keys(res)}')

        print('=== 負例 3: 入口の分類だけの変更 ===')
        # (a) 変換時 nodeInletCornerWall 1 の node、run の bcond で入口 (physID 1) を slip に
        write_bc(d0, {'inlet': (1, 'slip')})
        res = mp.check_run(d0)
        check(fail_keys(res) == ['bcond.1'] and '入口の分類が変わった' in res['fails'][0][1], f'FAIL {res["fails"]}')
        rc, v, _ = cli('check', d0)
        check(rc == 1 and v == 'VERDICT: FAIL', f'CLI → {v} (rc={rc})')
        # (b) 入口の分類を変えて変換し直したのに旧い記録が残っている (corner 1: 境界の所属が変わる)
        d3 = E.run_dir(', discretization: "node", nodeInletCornerWall: 1')
        E.convert(d3)
        old = open(os.path.join(d3, 'mesh.h5.provenance.json')).read()
        write_bc(d3, {'inlet': (1, 'slip')})
        E.raw_convert(d3)
        with open(os.path.join(d3, 'mesh.h5.provenance.json'), 'w') as fh:
            fh.write(old)
        write_bc(d3)                                                   # run の bcond は記録と同じに戻す (記録との照合だけを試す)
        res = mp.check_run(d3)
        msg = dict(res['fails']).get('(provenance)', '')
        check(fail_keys(res) == ['(provenance)'] and 'boundary' in msg and 'coord' not in msg.split('(違う部分:')[1].split(')')[0],
              f'変換し直した h5 に旧い記録 → FAIL {msg}')
        # (c) corner 0 を明記した node では入口の分類は変換に効かない → WARN (終了 0)
        d4 = E.run_dir(', discretization: "node", nodeInletCornerWall: 0')
        E.convert(d4)
        write_bc(d4, {'inlet': (1, 'slip')})
        res = mp.check_run(d4)
        check(not res['fails'] and [k for k, _ in res['warns']] == ['bcond.1'], f'corner 0 → WARN {res["warns"]}')
        # (d) corner を書いていない node (変換器の既定。値を確かめられない) → FAIL
        d5 = E.run_dir(', discretization: "node"')
        E.convert(d5)
        write_bc(d5, {'inlet': (1, 'outlet_statPress')})
        res = mp.check_run(d5)
        check(fail_keys(res) == ['bcond.1'], f'corner 省略 → FAIL {fail_keys(res)}')
        # (e) 同じ分類の中の種類の違いは注記だけ
        write_bc(d5, {'inlet': (1, 'inlet_uniformVelocity')})
        res = mp.check_run(d5)
        check(not res['fails'] and not res['warns'] and any(k == 'bcond.1' for k, _ in res['notes']),
              f'inlet_Pressure → inlet_uniformVelocity は注記 {res["notes"]}')
        # (f) 壁の分類の変更 (slip → wall) は corner によらず FAIL
        write_bc(d4, {'buffer': (5, 'wall')})
        res = mp.check_run(d4)
        check(fail_keys(res) == ['bcond.5'] and '壁の分類' in res['fails'][0][1], f'slip → wall → FAIL {res["fails"]}')
        write_bc(d4, {'wall': (3, 'wall_isothermal')})
        res = mp.check_run(d4)
        check(not res['fails'] and not res['warns'], 'wall → wall_isothermal (同じ壁の分類) は通す')
        bcs = dict(BASE_BC)
        del bcs['buffer']
        with open(os.path.join(d4, 'bcondConfig.yaml'), 'w') as fh:
            for nm, (pid, kind) in bcs.items():
                fh.write(f'{nm}: {{physID: {pid}, kind: {kind}, outputHDFflg: 0, ints: , floats: }}\n')
        res = mp.check_run(d4)
        check(fail_keys(res) == ['bcond.5'], f'変換時の physID 5 が run の bcond に無い → FAIL {fail_keys(res)}')

        print('=== 負例 4: wallDistExtraPhysIDs だけの変更 ===')
        d6 = E.run_dir(', discretization: "node", nodeInletCornerWall: 1')
        E.convert(d6)
        write_cfg(d6, ', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [5]')
        res = mp.check_run(d6)
        check(fail_keys(res) == ['mesh.wallDistExtraPhysIDs'], f'run だけ [5] → FAIL {res["fails"]}')
        # 変換し直した h5 (wall_dist だけ違う) に旧い記録
        old = open(os.path.join(d6, 'mesh.h5.provenance.json')).read()
        E.raw_convert(d6)                                              # 今の config ([5]) で変換
        with open(os.path.join(d6, 'mesh.h5.provenance.json'), 'w') as fh:
            fh.write(old)
        write_cfg(d6, ', discretization: "node", nodeInletCornerWall: 1')
        res = mp.check_run(d6)
        msg = dict(res['fails']).get('(provenance)', '')
        check(fail_keys(res) == ['(provenance)'] and '違う部分: wall_dist)' in msg, f'wall_dist だけ違う h5 → FAIL {msg}')
        # 順序・重複は効かない
        d7 = E.run_dir(', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [5, 2]')
        E.convert(d7)
        write_cfg(d7, ', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [2, 5, 5]')
        res = mp.check_run(d7)
        check(not res['fails'], f'[5, 2] と [2, 5, 5] は同じ → {res["fails"]}')
        write_cfg(d7, ', discretization: "node", nodeInletCornerWall: 1')
        res = mp.check_run(d7)
        check(fail_keys(res) == ['mesh.wallDistExtraPhysIDs'], f'変換時 [2, 5]・run は書いていない → FAIL {fail_keys(res)}')

        print('=== 負例 5: nodeInletCornerWall の不一致 ===')
        d8 = E.run_dir(', discretization: "node", nodeInletCornerWall: 0')
        E.convert(d8)
        write_cfg(d8, ', discretization: "node", nodeInletCornerWall: 1')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.nodeInletCornerWall'], f'変換時 0・run 1 → FAIL {res["fails"]}')
        rc, v, _ = cli('check', d8)
        check(rc == 1 and v == 'VERDICT: FAIL', f'CLI → {v} (rc={rc})')
        d9 = E.run_dir(', discretization: "node"')
        E.convert(d9)
        write_cfg(d9, ', discretization: "node", nodeInletCornerWall: 1')
        res = mp.check_run(d9)
        check(fail_keys(res) == ['mesh.nodeInletCornerWall'], f'変換時は書いていない・run 1 → FAIL {fail_keys(res)}')
        write_cfg(d9, ', discretization: "node", nodeInletCornerWall: 0')
        res = mp.check_run(d9)
        check(fail_keys(res) == ['mesh.nodeInletCornerWall'],
              f'変換時は書いていない・run 0 → FAIL (今のソースの既定 0 で代用しない) {fail_keys(res)}')

        print('=== その他の変換時の条件 ===')
        write_cfg(d8, ', discretization: "cell", nodeInletCornerWall: 0')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.discretization'], f'node で変換・run は cell → FAIL {fail_keys(res)}')
        write_cfg(d8, ', nodeInletCornerWall: 0')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.discretization'], f'node で変換・run は離散化を書いていない → FAIL {fail_keys(res)}')
        write_cfg(d8, ', discretization: "node", nodeInletCornerWall: 0, axisCentroidShift: 0')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.axisCentroidShift'], f'axisCentroidShift の違い → FAIL {fail_keys(res)}')
        write_cfg(d8, ', discretization: "node", nodeInletCornerWall: 0, renumber: "rcm"')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.renumber'], f'renumber の違い → FAIL {fail_keys(res)}')
        write_cfg(d8, ', discretization: "node", nodeInletCornerWall: 0, isAxisymmetric: 1')
        res = mp.check_run(d8)
        check(fail_keys(res) == ['mesh.isAxisymmetric'], f'平面で変換・run は軸対称 → FAIL {fail_keys(res)}')
        dA2 = E.run_dir(', discretization: "node", nodeInletCornerWall: 0, isAxisymmetric: 1, renumber: "rcm"')
        _, recA = E.convert(dA2)
        check(recA['artifact']['rSurfVect'] and recA['artifact']['renumber_perm'], f'成果物の事実 {recA["artifact"]}')
        write_cfg(dA2, ', discretization: "node", nodeInletCornerWall: 0, renumber: "rcm"')
        res = mp.check_run(dA2)
        check(not res['fails'] and [k for k, _ in res['warns']] == ['mesh.isAxisymmetric'], f'軸対称で変換・run は平面 → WARN {res["warns"]}')
        dA3 = E.run_dir(', discretization: "node", nodeInletCornerWall: 0', physprop='physProp: {isAxisymmetric: 1}\n')
        _, recA3 = E.convert(dA3)
        check(recA3['conversion']['mesh.isAxisymmetric'] == 1, '旧 physProp.isAxisymmetric も solver と同じ優先順で読む')
        dc = E.run_dir(', discretization: "cell"')
        _, recc = E.convert(dc)
        write_bc(dc, {'inlet': (1, 'slip')})
        res = mp.check_run(dc)
        check(not recc['artifact']['node'] and not res['fails'] and [k for k, _ in res['warns']] == ['bcond.1'],
              f'cell 変換: 入口の分類の変更は WARN {res["warns"]}')

        print('=== 正例: 保存量だけを更新した restart (restart_field.py) ===')
        dr = E.run_dir(', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [5]')
        E.convert(dr)
        mesh_h5 = os.path.join(dr, 'mesh.h5')
        res_h5 = os.path.join(E.root, 'res_00001000.h5')
        with h5py.File(mesh_h5, 'r') as f, h5py.File(res_h5, 'w') as g:
            n = f['VALUE/ro'].shape[0]
            rng = np.random.default_rng(7)
            g['MESH/COORD'] = np.asarray(f['MESH/COORD'])
            g['MESH/CONNE'] = np.asarray(f['VIZMESH/CONNE'])  # forge の node の出力は primal の接続 (output.cpp:154)
            for k in ('ro', 'roUx', 'roUy', 'roUz', 'roe'):
                g['VALUE/' + k] = (np.asarray(f['VALUE/' + k]) + rng.random(n)).astype(np.float32)
            g['VALUE/wall_dist'] = np.asarray(f['VALUE/wall_dist']).astype(np.float32)
            g['VALUE/P'] = np.ones(n, np.float32)
        # res は run の外に置いたので、その場の格子を --src-mesh で与える (restart_field.py の同一メッシュの検査、plan §5.1 #20)
        r = subprocess.run([sys.executable, os.path.join(TOOLS, 'restart_field.py'), res_h5, mesh_h5, '--src-mesh', mesh_h5],
                           capture_output=True, text=True)
        check(r.returncode == 0 and 'VERDICT: OK' in r.stdout, f'restart_field.py rc={r.returncode} {r.stdout.strip().splitlines()[-1:]}')
        with h5py.File(mesh_h5, 'r') as f, h5py.File(res_h5, 'r') as g:
            moved = all(np.array_equal(f['VALUE/' + k], g['VALUE/' + k]) for k in ('ro', 'roUx', 'roe'))
        check(moved, '保存量が写った (初期場は変わった)')
        res = mp.check_run(dr)
        check(not res['fails'] and not res['warns'] and any(k == '(initial_field)' for k, _ in res['notes']),
              f'restart の後のメッシュ → PASS、初期場の更新は注記 {res["notes"]}')
        rc, v, _ = cli('check', dr)
        check(rc == 0 and v == 'VERDICT: PASS', f'CLI → {v} (rc={rc})')
        # 化学種の属性 (ルートの属性) を足しても通る
        with h5py.File(mesh_h5, 'r+') as f:
            f.attrs['species_hash'] = 'abc'
            f['VALUE'].create_dataset('roY0', data=np.ones(n, np.float32))
        res = mp.check_run(dr)
        check(not res['fails'], 'ルートの属性・/VALUE の新しい量は照合に効かない')
        # wall_dist を書き換えたら FAIL (不変部分)
        with h5py.File(mesh_h5, 'r+') as f:
            w = np.asarray(f['VALUE/wall_dist'])
            w[0] += 1e-9
            f['VALUE/wall_dist'][...] = w
        res = mp.check_run(dr)
        check(fail_keys(res) == ['(provenance)'] and 'wall_dist' in res['fails'][0][1], f'wall_dist の書き換え → FAIL {fail_keys(res)}')
        w[0] -= 1e-9
        with h5py.File(mesh_h5, 'r+') as f:
            f['VALUE/wall_dist'][...] = w
        check(not mp.check_run(dr)['fails'], '元に戻せば通る')

        print('=== valueFileName が別ファイル (res から直接の restart) ===')
        shutil.copy(res_h5, os.path.join(dr, 'res_00001000.h5'))
        write_cfg(dr, ', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [5]', value='res_00001000.h5')
        res = mp.check_run(dr)
        check(not res['fails'], f'float32 の wall_dist (メッシュの double を丸めた値) → 通る {res["fails"]}')
        with h5py.File(os.path.join(dr, 'res_00001000.h5'), 'r+') as g:
            g['VALUE/wall_dist'][3] = g['VALUE/wall_dist'][3] * 1.5 + 0.01
        res = mp.check_run(dr)
        check(fail_keys(res) == ['mesh.valueFileName'] and 'wall_dist' in res['fails'][0][1], f'違う wall_dist の場 → FAIL {res["fails"]}')
        with h5py.File(os.path.join(dr, 'res_00001000.h5'), 'r+') as g:
            del g['VALUE/wall_dist']
        res = mp.check_run(dr)
        check(fail_keys(res) == ['mesh.valueFileName'], 'wall_dist の無い場 → FAIL')
        write_cfg(dr, ', discretization: "node", nodeInletCornerWall: 1, wallDistExtraPhysIDs: [5]', value='missing.h5')
        check(fail_keys(mp.check_run(dr)) == ['mesh.valueFileName'], '無い valueFileName → FAIL')

        print('=== 検査できない入力は合格にしない ===')
        dz = E.run_dir()
        E.convert(dz)
        pz = os.path.join(dz, 'mesh.h5.provenance.json')
        good = open(pz).read()
        with open(pz, 'w') as fh:
            fh.write(good[:len(good) // 2])
        check(fail_keys(mp.check_run(dz)) == ['(provenance)'], '壊れた記録 → FAIL')
        bad = json.loads(good)
        bad['version'] = 99
        with open(pz, 'w') as fh:
            json.dump(bad, fh)
        check(fail_keys(mp.check_run(dz)) == ['(provenance)'], '版の違う記録 → FAIL')
        bad = json.loads(good)
        bad['parts']['scheme'] = 'mesh-parts-v0'
        with open(pz, 'w') as fh:
            json.dump(bad, fh)
        check(fail_keys(mp.check_run(dz)) == ['(provenance)'], 'ハッシュの方式の違う記録 → FAIL')
        with open(pz, 'w') as fh:
            fh.write(good)
        check(not mp.check_run(dz)['fails'], '元の記録 → 通る')
        with open(os.path.join(dz, 'solverConfig.yaml'), 'w') as fh:
            fh.write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node", discretization: "cell"}\n')
        check(fail_keys(mp.check_run(dz)) == ['(config)'], '重複キーの config → FAIL')
        with open(os.path.join(dz, 'solverConfig.yaml'), 'w') as fh:
            fh.write('gpu: 1\n')
        check(fail_keys(mp.check_run(dz)) == ['mesh.meshFileName'], 'meshFileName の無い config → FAIL')
        write_cfg(dz, ', discretization: "node"', name='nothere.h5')
        check(fail_keys(mp.check_run(dz)) == ['mesh.meshFileName'], '無いメッシュ → FAIL')
        write_cfg(dz, ', discretization: "node"')
        os.remove(os.path.join(dz, 'bcondConfig.yaml'))
        check(fail_keys(mp.check_run(dz)) == ['(bcond)'], 'bcondConfig の無い run → FAIL')
        with open(os.path.join(dz, 'bcondConfig.yaml'), 'w') as fh:
            fh.write('inlet: {physID: 1, kind: inlet_Pressure}\noutlet: {physID: 1, kind: outlet_statPress}\n')
        check(fail_keys(mp.check_run(dz)) == ['(bcond)'], 'physID が重複した bcond → FAIL')
        with open(os.path.join(dz, 'mesh.h5'), 'wb') as fh:
            fh.write(b'not an hdf5 file')
        write_bc(dz)
        check(fail_keys(mp.check_run(dz)) == ['(mesh)'], '読めないメッシュ → FAIL')

        print('=== 変換の負例: 記録を書かない ===')
        dn = E.run_dir(', discretization: "node"')
        E.convert(dn)
        pn = os.path.join(dn, 'mesh.h5.provenance.json')
        try:
            E.convert(dn, {'FAKE_CONV_RC': '1'})
            check(False, '変換器の失敗で例外にならない')
        except mp.ProvenanceError as e:
            check(not os.path.exists(pn), f'変換器の失敗 → 記録しない・旧い記録を消す ({str(e)[:40]}…)')
        try:
            E.convert(dn, {'FAKE_CONV_IGNORE_CFG': '1'})
            check(False, 'config を読まない変換器で例外にならない')
        except mp.ProvenanceError as e:
            check('/VIZMESH' in str(e) and not os.path.exists(pn), 'config と成果物の食い違い (node の要求に cell の成果物) → 記録しない')
        try:
            E.convert(dn, {'FAKE_CONV_TOUCH': '1'})
            check(False, '変換中の書き換えで例外にならない')
        except mp.ProvenanceError as e:
            check('書き換わった' in str(e) and not os.path.exists(pn), '変換中に config が書き換わった → 記録しない')
        write_cfg(dn, ', discretization: "node"')
        try:
            E.convert(dn, {'FAKE_CONV_NOWRITE': '1'})
            check(False, '書かずに 0 で終わる変換器で例外にならない')
        except mp.ProvenanceError as e:
            check('書き換わっていない' in str(e) and not os.path.exists(pn), '変換器が h5 を書かずに 0 で終わった → 旧い h5 を記録しない')
        db = E.run_dir(', discretization: "node"')
        E.raw_convert(db)
        try:                                                            # 成果物の bcondKind と違う bcond で記録しようとする
            mp.build_record(os.path.join(db, 'mesh.h5'), open(os.path.join(db, 'solverConfig.yaml')).read(),
                            open(os.path.join(db, 'bcondConfig.yaml')).read().replace('kind: wall,', 'kind: wall_isothermal,'))
            check(False, 'bcondKind の食い違いで例外にならない')
        except mp.ProvenanceError as e:
            check('bcondKind' in str(e), 'bcondKind が変換時の bcond と違う → 記録しない')
        rc, v, out = cli('convert', 'mesh.json', 'mesh.h5', '--cwd', dn, '--converter', E.conv)
        check(rc == 0 and v == 'VERDICT: PASS' and os.path.exists(pn), f'CLI convert → {v} (rc={rc})')
        rc, v, out = cli('convert', 'mesh.json', 'mesh.h5', '--cwd', dn, '--converter', os.path.join(E.root, 'nope'))
        check(rc == 1 and v == 'VERDICT: FAIL' and not os.path.exists(pn), f'無い変換器 → {v} (rc={rc})、旧い記録は消す')

        print('=== 部分ごとのハッシュ (#20 が使う) ===')
        pa = mp.mesh_parts(os.path.join(dA, 'mesh.h5'))
        pb = mp.mesh_parts(os.path.join(dB, 'mesh.h5'))
        check(set(pa['parts']) == set(mp.INVARIANT_PARTS) | {'initial_field'}, f'部分 {sorted(pa["parts"])}')
        check(pa['parts']['coord'] != pb['parts']['coord'] and pa['parts']['topology'] == pb['parts']['topology'],
              '座標だけ違うメッシュ: coord が違い topology は同じ')
        with h5py.File(os.path.join(dA, 'mesh.h5'), 'r') as f:
            check(mp.mesh_parts(f)['invariant'] == pa['invariant'], '開いた h5py.File でも同じハッシュ')
    finally:
        shutil.rmtree(E.root, ignore_errors=True)

    print('ALL PASS' if fails == 0 else f'{fails} FAIL')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
