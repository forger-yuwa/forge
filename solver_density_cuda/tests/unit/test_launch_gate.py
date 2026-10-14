#!/usr/bin/env python3
"""投入の関門 (`launch_gate.py`) と `run_case.sh` の呼び出しの試験
(plan tooling-generator-defaults-and-launch-gate §5.1 #10、§6.2 の設定・メッシュの来歴・品質・宣言の負例)。

forge の代わりに偽のバイナリ (起動されたら印 `FORGE_STARTED` を書くだけの台本) を `FORGE_BIN` で渡し、本物の `run_case.sh` を
通しで回す。本物の forge は起動しない。合格条件 (plan §5.1 #10):

- 負例 (設定・メッシュの来歴・品質・宣言) で forge が起動しない (`run_case.sh` の終了コードが 0 でなく、印が無い)。
- 正例で起動する。
- `FORGE_LAUNCH_GATE=warn` は `fixed`・`explicit` の不一致・省略だけを WARN に落とし、`forbid`・未知のキー・宣言の欠落 (強制時)・
  メッシュの FAIL・検査不能は落とさない。
- バイナリの sha256 とプロファイルの内容のハッシュが `LAUNCH_GATE.txt` に記録され、`RUN_PROVENANCE.txt` に 1 行写る。
  初期場の h5 に `restart_field.py --allow-coord-upcast` の記録があれば両方に写る (判定は変えない; plan §5.1 #20b)。
- 移行期 (既定) は宣言の無い run・記録の無いメッシュを WARN で起動し、強制のモード (`FORGE_LAUNCH_GATE=enforce`・定数
  `ENFORCE`) では止める。
- 適用条件 (plan §5.1 #28): 宣言のプロファイルの次元・格子の種類を変換の成果物 (来歴の記録の `artifact`、無ければ h5) から照合する。
  2D/3D の取り違え (SERN 2D の宣言 + 六面体の格子、SERN 3D の宣言 + 四角形の格子) と構造/非構造の取り違え (SERN 2D の宣言 +
  三角形の格子) は名指しで FAIL、正しい組 (SERN 2D + 四角形、SERN 3D + 六面体) は通る。型を読めない格子は判定不能として
  宣言した run を FAIL (warn でも落とさない)。codex result 段 M2 の再現 (プロファイルの適用条件だけを `dimension: 3d` 等に) も FAIL。

メッシュは node 変換と同じ置き方の小さい合成の h5 (長方形の四角形の格子、上辺が壁。適用条件の試験は三角形・六面体・型を読めない
格子も) を直接書き、来歴の記録は `mesh_provenance.build_record` で作る (変換器は使わない)。宣言は `recipe_settings.declare` で
書く (プロファイル `common_node_ns_sst` の段 main。適用条件の試験は SERN 2D・3D の本段)。

実行: python3 solver_density_cuda/tests/unit/test_launch_gate.py   (約 1〜2 分、出力は一時ディレクトリ。全部合格なら消す)
"""
import contextlib
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, '..', '..', 'tools'))
sys.path.insert(0, TOOLS)
import launch_gate as LG          # noqa: E402
import mesh_provenance as MP      # noqa: E402
import recipe_settings as RS      # noqa: E402

RUN_CASE = os.path.join(TOOLS, 'run_case.sh')
TMP = tempfile.mkdtemp(prefix='launch_gate_test_')
RC = RS.load_recipes()
PROFILE, STAGE = 'common_node_ns_sst', 'main'
fails = 0


def check(cond, msg, detail=''):
    global fails
    print(('ok   ' if cond else 'FAIL ') + msg + (f'  [{detail}]' if detail and not cond else ''), flush=True)
    if not cond:
        fails += 1


# ------------------------------------------------------------------------------------------------ 合成の run

PROBLEM = {'mesh.meshFileName': 'mesh.h5', 'mesh.valueFileName': 'mesh.h5', 'gpu': 1, 'initial': 'uniform_p101325_u10',
           'physProp.cp': 1004.5, 'physProp.gamma': 1.4, 'physProp.visc': 1.8e-5, 'physProp.thermCond': 0.0257,
           'time.last.nStepOuter': 12000, 'time.outStepInterval': 2000, 'time.outStepStart': 0,
           'time.deltaT.dt': 1e-8, 'time.deltaT.dt_min': 1e-9, 'time.deltaT.dt_max': 1e-3, 'time.deltaT.cfl_pseudo': 4.0}
BCOND = {'inlet': (1, 'inlet_Pressure'), 'outlet': (2, 'outlet_statPress'), 'wall': (3, 'wall'), 'bottom': (4, 'slip')}


def unflatten(flat):
    out = {}
    for k, v in flat.items():
        cur = out
        parts = k.split('.')
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    return out


def cfg_flat(profile=PROFILE, stage=STAGE):
    res = RS.resolve(RC, profile)
    flat = dict(PROBLEM)
    flat.update(res['stages'][stage]['fixed'])
    flat.update(res['stages'][stage]['explicit'])
    return flat


def cfg_text(flat):
    return yaml.safe_dump(unflatten(flat), sort_keys=False)


def bcond_text(kinds=None):
    kinds = kinds or BCOND
    return ''.join(f'{nm}: {{physID: {pid}, kind: {k}, outputHDFflg: 0, ints: , floats: }}\n' for nm, (pid, k) in kinds.items())


def write_mesh(path, ni=21, nj=6, L=1.0, H=0.2, first=None, shift=0.0, kinds=None, cells='quad'):
    """node 変換と同じ置き方の合成の h5 (primal は /VIZMESH/CONNE、双対は節点ごとの 1 点の CV)。上辺 (physID 3) が壁。
    first を与えると壁 (上辺) の第一層の厚さをそれにする (AR = (L/(ni−1)) / first)。
    cells (plan §5.1 #28 の適用条件の試験): 'quad' (既定、2D の構造格子)・'tri' (四角形を対角で 2 つの三角形に割った 2D の
    非構造格子)・'hex' (z 方向に 2 層押し出した 3D の六面体; 境界は押し出した辺の四角形の面)・'polygon' (知らない型の code 3
    = 型を読めない格子)。"""
    if cells == 'hex':
        return _write_mesh_hex(path, ni, nj, L, H, kinds)
    kinds = kinds or BCOND
    xs = np.linspace(0.0, L, ni) + shift
    ys = np.linspace(0.0, H, nj)
    if first is not None:
        ys = np.r_[np.linspace(0.0, H - first, nj - 1), H]
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    coord = np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)]
    nid = lambda i, j: i * nj + j                                                       # noqa: E731
    quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    edges = {1: [(nid(0, j), nid(0, j + 1)) for j in range(nj - 1)],
             2: [(nid(ni - 1, j), nid(ni - 1, j + 1)) for j in range(nj - 1)],
             3: [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)],
             4: [(nid(i, 0), nid(i + 1, 0)) for i in range(ni - 1)]}
    pid_kind = {pid: k for _, (pid, k) in kinds.items()}
    with h5py.File(path, 'w') as f:
        g = f.create_group('MESH')
        g.attrs['nCells'] = np.int64(len(coord))
        f['MESH/COORD'] = coord.ravel()
        f['MESH/CONNE'] = np.array([x for n in range(len(coord)) for x in (1, n)], dtype=np.int64)
        if cells == 'tri':
            tris = [t for q in quads for t in ([q[0], q[1], q[2]], [q[0], q[2], q[3]])]
            f['VIZMESH/CONNE'] = np.array([x for t in tris for x in [4] + t], dtype=np.int64)
            f['VIZMESH'].attrs['nVizCells'] = np.int64(len(tris))
        else:
            code = [3, 4] if cells == 'polygon' else [5]
            f['VIZMESH/CONNE'] = np.array([x for q in quads for x in code + q], dtype=np.int64)
            f['VIZMESH'].attrs['nVizCells'] = np.int64(len(quads))
        f['CELLS/volume'] = np.ones(len(coord))
        for k in ('ro', 'roUx', 'roUy', 'roUz', 'roe'):
            f['VALUE/' + k] = np.ones(len(coord), np.float32)
        f['VALUE/wall_dist'] = (H - coord[:, 1])
        for pid, es in edges.items():
            b = f.create_group(f'BCONDS/{pid}')
            nodes = sorted({n for e in es for n in e})
            b['iCells'] = np.array(nodes, dtype=np.int64)
            b['vizBfaceNodes'] = np.array([n for e in es for n in e], dtype=np.int64)
            b['vizBfaceSizes'] = np.full(len(es), 2, dtype=np.int64)
            b.attrs['bcondKind'] = pid_kind[pid]


def _write_mesh_hex(path, ni, nj, L, H, kinds, nk=3, W=0.1):
    """write_mesh の 3D 版 (cells='hex'): 長方形の格子を z 方向に nk−1 層押し出した六面体。境界は押し出した 4 辺の四角形の面
    (z の両端の面は境界にしない; 品質の検査と来歴の記録に要るのは壁の面だけ)。"""
    kinds = kinds or BCOND
    xs, ys, zs = np.linspace(0.0, L, ni), np.linspace(0.0, H, nj), np.linspace(0.0, W, nk)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    coord = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    nid = lambda i, j, k: (i * nj + j) * nk + k                                        # noqa: E731
    hexes = [[nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
              nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
             for i in range(ni - 1) for j in range(nj - 1) for k in range(nk - 1)]
    faces = {1: [(nid(0, j, k), nid(0, j + 1, k), nid(0, j + 1, k + 1), nid(0, j, k + 1)) for j in range(nj - 1) for k in range(nk - 1)],
             2: [(nid(ni - 1, j, k), nid(ni - 1, j + 1, k), nid(ni - 1, j + 1, k + 1), nid(ni - 1, j, k + 1))
                 for j in range(nj - 1) for k in range(nk - 1)],
             3: [(nid(i, nj - 1, k), nid(i + 1, nj - 1, k), nid(i + 1, nj - 1, k + 1), nid(i, nj - 1, k + 1))
                 for i in range(ni - 1) for k in range(nk - 1)],
             4: [(nid(i, 0, k), nid(i + 1, 0, k), nid(i + 1, 0, k + 1), nid(i, 0, k + 1)) for i in range(ni - 1) for k in range(nk - 1)]}
    pid_kind = {pid: k for _, (pid, k) in kinds.items()}
    with h5py.File(path, 'w') as f:
        g = f.create_group('MESH')
        g.attrs['nCells'] = np.int64(len(coord))
        f['MESH/COORD'] = coord.ravel()
        f['MESH/CONNE'] = np.array([x for n in range(len(coord)) for x in (1, n)], dtype=np.int64)
        f['VIZMESH/CONNE'] = np.array([x for h in hexes for x in [9] + h], dtype=np.int64)
        f['VIZMESH'].attrs['nVizCells'] = np.int64(len(hexes))
        f['CELLS/volume'] = np.ones(len(coord))
        for k in ('ro', 'roUx', 'roUy', 'roUz', 'roe'):
            f['VALUE/' + k] = np.ones(len(coord), np.float32)
        f['VALUE/wall_dist'] = (H - coord[:, 1])
        for pid, fs in faces.items():
            b = f.create_group(f'BCONDS/{pid}')
            b['iCells'] = np.array(sorted({n for fc in fs for n in fc}), dtype=np.int64)
            b['vizBfaceNodes'] = np.array([n for fc in fs for n in fc], dtype=np.int64)
            b['vizBfaceSizes'] = np.full(len(fs), 4, dtype=np.int64)
            b.attrs['bcondKind'] = pid_kind[pid]


FAKE_FORGE = os.path.join(TMP, 'fake_forge.sh')
with open(FAKE_FORGE, 'w') as fh:
    fh.write('#!/bin/sh\n# 偽の forge (試験用): 起動されたら印を書くだけ\necho started > FORGE_STARTED\nexit 0\n')
os.chmod(FAKE_FORGE, os.stat(FAKE_FORGE).st_mode | stat.S_IXUSR)
FAKE_SHA = hashlib.sha256(open(FAKE_FORGE, 'rb').read()).hexdigest()
PROFILE_SHA = RS.resolve(RC, PROFILE)['sha256']

_n = [0]


def make_run(tag, flat=None, kinds=None, record=True, declare=True, conv_flat=None, mesh_kw=None, conv_kinds=None,
             deviate=None, reason=None, experiment=False, stage=STAGE, profile=PROFILE):
    """run を作る。record: 来歴の記録を書く (conv_flat・conv_kinds = 変換時の config・bcond、既定は run と同じ)。
    declare: recipe_settings.declare で宣言する (通らなければ例外)。"""
    _n[0] += 1
    d = os.path.join(TMP, f'{_n[0]:02d}_{tag}')
    os.makedirs(d)
    flat = cfg_flat() if flat is None else flat
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as fh:
        fh.write(cfg_text(flat))
    with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as fh:
        fh.write(bcond_text(kinds))
    write_mesh(os.path.join(d, 'mesh.h5'), kinds=conv_kinds or kinds, **(mesh_kw or {}))
    if record:
        h5 = os.path.join(d, 'mesh.h5')
        MP.write_record(h5, MP.build_record(h5, cfg_text(conv_flat or flat), bcond_text(conv_kinds or kinds)))
    if declare:
        rc, r, path = RS.declare(d, None if experiment else profile, None if experiment else stage, deviate, reason,
                                 experiment=experiment, recipes=RC)
        assert rc == 0, (tag, r['fails'])
    return d


def gate_inproc(d, **kw):
    """launch_gate.run_gate をプロセス内で回す (LAUNCH_GATE.txt の本文の表示は捨てる)。"""
    with contextlib.redirect_stdout(io.StringIO()):
        return LG.run_gate(d, **kw)


def edit_config(d, flat=None, text=None, rebind=True):
    """config を書き換える。rebind なら宣言の config の sha256 を新しい config に付け替える (生成器が段ごとに宣言を書き直した形;
    宣言はそのままに、照合の中身だけを試すため)。"""
    p = os.path.join(d, 'solverConfig.yaml')
    text = cfg_text(flat) if text is None else text
    with open(p, 'w') as fh:
        fh.write(text)
    if rebind:
        rp = os.path.join(d, 'RECIPE.json')
        rec = json.load(open(rp))
        rec['config']['sha256'] = hashlib.sha256(text.encode()).hexdigest()
        json.dump(rec, open(rp, 'w'), ensure_ascii=False)


def run_case(d, **env_extra):
    """本物の run_case.sh を偽の forge で回す。戻り (終了コード, stdout, 起動したか, LAUNCH_GATE.txt)。"""
    env = dict(os.environ, FORGE_BIN=FAKE_FORGE)
    env.pop('FORGE_LAUNCH_GATE', None)
    env.pop('FORGE_STRICT_CONVERGENCE', None)
    env.update(env_extra)
    for f in ('FORGE_STARTED', 'LAUNCH_GATE.txt'):
        if os.path.exists(os.path.join(d, f)):
            os.unlink(os.path.join(d, f))
    p = subprocess.run(['bash', RUN_CASE, d], capture_output=True, text=True, env=env)
    gate = open(os.path.join(d, 'LAUNCH_GATE.txt'), encoding='utf-8').read() if os.path.exists(os.path.join(d, 'LAUNCH_GATE.txt')) else ''
    return p.returncode, p.stdout + p.stderr, os.path.exists(os.path.join(d, 'FORGE_STARTED')), gate


def verdict(gate):
    lines = [ln for ln in gate.splitlines() if ln.startswith('VERDICT:')]
    return lines[-1].split()[1] if lines else None


def items(gate, level):
    """LAUNCH_GATE.txt の FAIL / WARN の行 → [(区分, キー)]。"""
    out = []
    for ln in gate.splitlines():
        if ln.startswith(level + ' ') and ': [' in ln:
            body = ln.split(': [', 1)[1]
            area, rest = body.split('] ', 1)
            out.append((area, rest.split(': ', 1)[0]))
    return out


def keys(gate, level):
    return {k for _, k in items(gate, level)}


def field(gate, name):
    for ln in gate.splitlines():
        if ln.startswith(name.ljust(12) + ':') or ln.startswith(name + ':') or ln.startswith(name.ljust(12)[:12] + ':'):
            return ln.split(':', 1)[1].strip()
    return None


def stopped(rc, out, started, gate, why):
    """負例: run_case.sh が 0 でなく、forge (偽物) を起動せず、関門が FAIL を書いた。"""
    check(rc != 0 and not started and verdict(gate) == 'FAIL' and '[run_case] launch gate exit=' in out
          and '[run_case] forge exit=' not in out, f'{why} -> forge を起動しない (rc {rc})',
          f'started {started} verdict {verdict(gate)} out {out[-300:]!r}')


def launched(rc, out, started, gate, why, want=('PASS', 'WARN')):
    check(rc == 0 and started and verdict(gate) in want and '[run_case] forge exit=0' in out,
          f'{why} -> 起動する (VERDICT {verdict(gate)})', f'rc {rc} started {started} gate:\n{gate[-1500:]}')


# ------------------------------------------------------------------------------------------------ 正例

print('== 正例 (宣言・来歴の記録・品質が揃った run)')
good = make_run('good')
rc, out, started, gate = run_case(good)
launched(rc, out, started, gate, '宣言 + 記録 + 品質 PASS')
check(field(gate, 'forge_sha256') == FAKE_SHA, 'LAUNCH_GATE.txt にバイナリの sha256 (偽の forge の実測と同じ)',
      field(gate, 'forge_sha256'))
check(f'(sha256 {PROFILE_SHA})' in (field(gate, 'profile') or '') and field(gate, 'profile').startswith(PROFILE),
      'LAUNCH_GATE.txt にプロファイルの内容のハッシュ (recipe_settings.resolve の sha256 と同じ)', field(gate, 'profile'))
check(field(gate, 'stage') == STAGE, 'LAUNCH_GATE.txt に段')
rec = json.load(open(os.path.join(good, 'mesh.h5.provenance.json')))
check((field(gate, 'mesh_hash') or '').startswith(rec['parts']['invariant']) and '来歴の記録と一致' in field(gate, 'mesh_hash'),
      'LAUNCH_GATE.txt にメッシュの不変部分のハッシュ (来歴の記録と同じ)', field(gate, 'mesh_hash'))
check((field(gate, 'quality') or '').startswith('PASS (再利用 なし'), '品質 PASS (初回は再利用なし)', field(gate, 'quality'))
check(keys(gate, 'FAIL') == set() and items(gate, 'WARN') == [('バイナリ', FAKE_FORGE)],
      'WARN は「対応の分からないバイナリ」だけ (拒否はしない)', str(items(gate, 'WARN')))
prov = open(os.path.join(good, 'RUN_PROVENANCE.txt'), encoding='utf-8').read().splitlines()
gl = [ln for ln in prov if ln.startswith('launch_gate :')]
check(len(gl) == 1 and 'VERDICT WARN' in gl[0] and f'profile {PROFILE}' in gl[0] and FAKE_SHA[:16] in gl[0]
      and any(ln.startswith('forge_sha256: ' + FAKE_SHA) for ln in prov),
      'RUN_PROVENANCE.txt に関門の 1 行 (run_case.sh の forge_sha256 と同じバイナリ)', str(gl))
rc, out, started, gate = run_case(good)
launched(rc, out, started, gate, '同じ run をもう一度')
check((field(gate, 'quality') or '').startswith('PASS (再利用 あり'), '2 回目は品質の結果を再利用 (鍵が一致)', field(gate, 'quality'))
qj = os.path.join(good, 'mesh.h5.quality.json')
q = json.load(open(qj))
q['key']['args']['ar_max'] = 1000.0
json.dump(q, open(qj, 'w'))
rc, out, started, gate = run_case(good)
check(started and (field(gate, 'quality') or '').startswith('PASS (再利用 なし'), '鍵の違う品質の結果は再利用しない',
      field(gate, 'quality'))

# 対応の分かるバイナリ (このリポジトリの build/forge) なら WARN も無く PASS (プロセス内で既定のパスを偽物に向ける)
_def = LG.DEFAULT_FORGE
LG.DEFAULT_FORGE = FAKE_FORGE
try:
    rc_, res_, info_ = gate_inproc(good, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
finally:
    LG.DEFAULT_FORGE = _def
check(rc_ == 0 and res_.verdict() == 'PASS', '既定の build/forge (と同じ実体) なら PASS', str(res_.items))
rc_, res_, info_ = gate_inproc(good, forge_bin=os.path.join(TMP, 'no_such_forge'), environ={}, dry_run=True)
check(rc_ == 0 and info_['forge_sha256'].startswith('不明') and ('バイナリ', os.path.join(TMP, 'no_such_forge')) in
      [(a, k) for lv, a, k, _ in res_.items if lv == 'WARN'], '無いバイナリは記録 (不明) と WARN (関門は止めない)')

# 初期場の来歴 (restart_field.py --allow-coord-upcast の記録、plan §5.1 #20b): 判定は変えず、LAUNCH_GATE.txt と RUN_PROVENANCE.txt に写す
up = make_run('coord_upcast')
rec_up = {'kind': 'coord_upcast', 'tool': 'restart_field.py --allow-coord-upcast', 'src_coord_dtype': 'float32',
          'src_invariant': 'a' * 64, 'dst_coord_dtype': 'float64', 'dst_invariant': 'b' * 64, 'parts_scheme': MP.PARTS_SCHEME,
          'date': '2026-10-15T00:00:00+09:00'}
with h5py.File(os.path.join(up, 'mesh.h5'), 'r+') as f:
    MP.write_coord_upcast(f, rec_up)
rc, out, started, gate = run_case(up)
launched(rc, out, started, gate, '初期場に座標の精度の引き上げの記録がある run')
check('座標の精度の引き上げ' in (field(gate, 'initial_field') or '') and 'a' * 16 + ' → float64 の格子 ' + 'b' * 16 in field(gate, 'initial_field')
      and '来歴の記録と一致' in (field(gate, 'mesh_hash') or ''),
      'LAUNCH_GATE.txt に初期場の来歴 (旧 → 新の不変部分の sha256)、ルートの記録は格子の来歴の照合を変えない', field(gate, 'initial_field'))
gl = [ln for ln in open(os.path.join(up, 'RUN_PROVENANCE.txt'), encoding='utf-8').read().splitlines() if ln.startswith('launch_gate :')]
check(len(gl) == 1 and '初期場 mesh.h5: 座標の精度の引き上げ' in gl[0], 'RUN_PROVENANCE.txt の関門の行に初期場の来歴', str(gl))
check(field(open(os.path.join(good, 'LAUNCH_GATE.txt'), encoding='utf-8').read(), 'initial_field') is None,
      '記録の無い run は initial_field の行を書かない')

# ------------------------------------------------------------------------------------------------ 設定の負例 (宣言あり)

print('== 設定の負例 (宣言あり、宣言の config の sha256 は書き換えた config に付け替える)')
base = cfg_flat()


def neg_setting(tag, flat=None, text=None, want_key=None, warn_drops=False):
    d = make_run(tag)
    edit_config(d, flat=flat, text=text)
    rc, out, started, gate = run_case(d)
    stopped(rc, out, started, gate, f'{tag}')
    if want_key is not None:
        check(want_key in keys(gate, 'FAIL'), f'{tag}: {want_key} を名指しで FAIL', str(items(gate, 'FAIL')))
    rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
    if warn_drops:
        launched(rc, out, started, gate, f'{tag} + FORGE_LAUNCH_GATE=warn (fixed の不一致は WARN に落とす)', want=('WARN',))
        check(any('FORGE_LAUNCH_GATE=warn で FAIL から WARN に落とした' in ln for ln in gate.splitlines() if ln.startswith('WARN')),
              f'{tag} + warn: 落とした旨を LAUNCH_GATE.txt に書く')
    else:
        stopped(rc, out, started, gate, f'{tag} + FORGE_LAUNCH_GATE=warn (落とさない)')
        if want_key is not None:
            check(want_key in keys(gate, 'FAIL'), f'{tag} + warn: {want_key} は FAIL のまま', str(items(gate, 'FAIL')))
    return gate


neg_setting('nStepInner5', flat=dict(base, **{'time.nStepInner': 5}), want_key='time.nStepInner', warn_drops=True)
fl = dict(base)
del fl['time.nStepInner']
neg_setting('omit_fixed', flat=fl, want_key='time.nStepInner', warn_drops=True)
neg_setting('bndFirstOrder', flat=dict(base, **{'mesh.bndFirstOrder': 1}), want_key='mesh.bndFirstOrder')
neg_setting('unknown_key', flat=dict(base, **{'time.deltaT.notAKey': 1}), want_key='time.deltaT.notAKey')
neg_setting('misplaced_key', flat=dict(base, **{'time.deltaT.nStepInner': 4}), want_key='time.deltaT.nStepInner')
g = neg_setting('dualtime_relax', flat=dict(base, **{'time.unsteady': 1, 'time.dualTime': 1, 'time.nSubIterDualTime': 20,
                                                     'time.deltaT.cfl_pseudo': 12.0}), want_key='time.deltaT.implicitRelax')
check({'time.unsteady', 'time.dualTime'} <= keys(g, 'WARN') and 'time.unsteady' not in keys(g, 'FAIL'),
      'dualtime_relax + warn: fixed の不一致 (unsteady・dualTime) だけ WARN に落ち、dual-time の緩和 (check_solver_config) は FAIL のまま',
      str(items(g, 'FAIL')))
neg_setting('off_profile_key', flat=dict(base, **{'space.slauWallNormalChi': 0}), want_key='space.slauWallNormalChi')
txt = cfg_text(base)
neg_setting('duplicate_key', text=txt + 'gpu: 0\n', want_key='gpu')
neg_setting('merge_key', text=txt.replace('turbulence:\n', 'turbulence:\n  <<: {sstEnergyIncludesK: 1}\n', 1), want_key='turbulence.<<')
d = make_run('stage_soft_main')
rp = os.path.join(d, 'RECIPE.json')
r_ = json.load(open(rp))
r_['stage'] = 'soft'
json.dump(r_, open(rp, 'w'), ensure_ascii=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '宣言した段 (soft) と config の段 (main) の不一致')
check({'space.convMethod', 'time.nStepInner'} <= keys(gate, 'FAIL'), '段の不一致は soft の fixed (次数・内部反復) を名指しで FAIL',
      str(items(gate, 'FAIL')))

# ------------------------------------------------------------------------------------------------ 宣言の負例・正例

print('== 宣言')
d = make_run('undeclared', declare=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '宣言の無い run (移行期)', want=('WARN',))
check(('宣言', 'RECIPE.json') in items(gate, 'WARN'), '宣言の欠落は WARN として記録')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='enforce')
stopped(rc, out, started, gate, '宣言の無い run (FORGE_LAUNCH_GATE=enforce)')
check(('宣言', 'RECIPE.json') in items(gate, 'FAIL'), '強制では宣言の欠落を FAIL として名指し')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='enforce,warn')
stopped(rc, out, started, gate, '宣言の無い run (enforce + warn: warn は宣言の欠落を落とさない)')
# 定数 ENFORCE の切り替え (強制の開始) — プロセス内で定数だけを変える
LG.ENFORCE = True
try:
    rc_, res_, _ = gate_inproc(d, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
    rc2, res2, _ = gate_inproc(good, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
finally:
    LG.ENFORCE = False
check(rc_ == 1 and ('FAIL', '宣言') in [(lv, a) for lv, a, _, _ in res_.items], '定数 ENFORCE = True: 宣言の無い run は FAIL')
check(rc2 == 0, '定数 ENFORCE = True: 宣言・記録の揃った run は通る', str(res2.items))
rc_, res_, _ = gate_inproc(d, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
check(rc_ == 0, '定数を戻せば移行期 (宣言の無い run は WARN)')
# 宣言が無くても check_solver_config の FAIL と未知のキーは止める
d = make_run('undeclared_bnd', flat=dict(base, **{'mesh.bndFirstOrder': 1}), declare=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '宣言の無い run の bndFirstOrder')
# 未知・誤配置のキーは、宣言が無いあいだは宣言の欠落と同じ水準 (移行期 WARN・強制 FAIL)。宣言した run では FAIL (上の unknown_key)
d = make_run('undeclared_unknown', flat=dict(base, **{'space.notAKey': 1}), declare=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '宣言の無い run の未知のキー (移行期)', want=('WARN',))
check(('設定', 'space.notAKey') in items(gate, 'WARN'), '宣言の無い run の未知のキーは名指しで WARN (移行期)')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='enforce')
stopped(rc, out, started, gate, '宣言の無い run の未知のキー (強制)')
check('space.notAKey' in keys(gate, 'FAIL'), '強制では宣言の無い run の未知のキーも名指しで FAIL')
# migrate_solver_config.py が直せるキー (time.last.time・旧 physProp.isAxisymmetric) は文言で migrate を案内する (plan §5.1 #25)
d = make_run('undeclared_lasttime', flat=dict(base, **{'time.last.time': 1.0, 'physProp.isAxisymmetric': 0}), declare=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '宣言の無い run の time.last.time・旧 physProp.isAxisymmetric (移行期)', want=('WARN',))
_mig = [ln for ln in gate.splitlines() if ln.startswith('WARN') and ('time.last.time' in ln or 'physProp.isAxisymmetric' in ln)]
check(len(_mig) == 2 and all('migrate_solver_config.py' in ln for ln in _mig),
      'LAUNCH_GATE.txt の time.last.time・physProp.isAxisymmetric の WARN が migrate_solver_config.py を案内する', str(_mig))
d = make_run('undeclared_rejected', flat=dict(base, **{'space.limiterMatchRecon': 1}), declare=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '宣言の無い run の solver が起動時に拒否するキー (移行期でも)')
d = make_run('decl_broken')
open(os.path.join(d, 'RECIPE.json'), 'w').write('{')
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '読めない宣言 (検査不能)')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
stopped(rc, out, started, gate, '読めない宣言 + warn (検査不能は落とさない)')
d = make_run('decl_stale')
edit_config(d, flat=dict(base, **{'time.deltaT.cfl_pseudo': 5.0}), rebind=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '宣言の後に config が変わった (sha256 の不一致)')
check(any('宣言の後に config が変わった' in ln for ln in gate.splitlines()), '宣言の config の sha256 の不一致を名指し')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
stopped(rc, out, started, gate, '宣言の後に config が変わった + warn (落とさない)')
d = make_run('decl_noreason', flat=dict(base, **{'time.nStepInner': 5}), deviate=['time.nStepInner=5'], reason='試し')
rp = os.path.join(d, 'RECIPE.json')
r_ = json.load(open(rp))
r_['deviations'][0]['reason'] = ''
json.dump(r_, open(rp, 'w'), ensure_ascii=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '理由の無い逸脱の宣言')
d = make_run('decl_unknown_profile')
r_ = json.load(open(os.path.join(d, 'RECIPE.json')))
r_['profile'] = 'no_such_profile'
json.dump(r_, open(os.path.join(d, 'RECIPE.json'), 'w'), ensure_ascii=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '無いプロファイルの宣言')
d = make_run('deviation_ok', flat=dict(base, **{'time.nStepInner': 5}), deviate=['time.nStepInner=5'], reason='試し: 内部反復 5')
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '理由つきの逸脱の宣言')
check(any(ln.startswith('deviation   : time.nStepInner = 5') and '理由 試し: 内部反復 5' in ln for ln in gate.splitlines()),
      'LAUNCH_GATE.txt に逸脱 (キー・値・理由)')
fl = dict(base)
del fl['time.nStepInner']
d = make_run('deviation_omit', flat=fl, deviate=['time.nStepInner=~'], reason='省略の試し')
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '「書かない」の逸脱 (key=~)')
check(any(ln.startswith('deviation   : time.nStepInner = ~') and '未確定' in ln for ln in gate.splitlines()),
      'LAUNCH_GATE.txt に「書かない」の逸脱と実効値 未確定')
d = make_run('experiment', flat=dict(base, **{'time.nStepInner': 7}), experiment=True, reason='実験: 内部反復の感度')
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '「実験」の宣言')
r_ = json.load(open(os.path.join(d, 'RECIPE.json')))
r_['reason'] = ' '
json.dump(r_, open(os.path.join(d, 'RECIPE.json'), 'w'), ensure_ascii=False)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '理由の無い「実験」の宣言')

# ------------------------------------------------------------------------------------------------ メッシュの来歴の負例

print('== メッシュの来歴')
d = make_run('no_record', record=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, '記録の無いメッシュ (移行期)', want=('WARN',))
check(('来歴', 'mesh.h5.provenance.json') in items(gate, 'WARN') and '来歴の記録なし' in (field(gate, 'mesh_hash') or ''),
      '記録の無いメッシュは WARN、メッシュのハッシュは今の h5 から記録', str(items(gate, 'WARN')))
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='enforce')
stopped(rc, out, started, gate, '記録の無いメッシュ (強制)')
d = make_run('other_record')
other = os.path.join(TMP, 'other_mesh.h5')
write_mesh(other, shift=0.5)
MP.write_record(other, MP.build_record(other, cfg_text(base), bcond_text()))
shutil.copy(other + '.provenance.json', os.path.join(d, 'mesh.h5.provenance.json'))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '別メッシュの記録')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
stopped(rc, out, started, gate, '別メッシュの記録 + warn (メッシュの FAIL は落とさない)')
d = make_run('corner_mismatch', conv_flat=dict(base, **{'mesh.nodeInletCornerWall': 0}))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '変換時の nodeInletCornerWall 0・run の config 1')
check('mesh.nodeInletCornerWall' in keys(gate, 'FAIL'), 'nodeInletCornerWall を名指しで FAIL')
d = make_run('inlet_class', kinds=dict(BCOND, inlet=(1, 'outlet_statPress')), conv_kinds=BCOND)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, '入口の分類だけを変えた bcond')

# symlink したメッシュ (plan §5.1 #23): 記録・品質の結果はメッシュの実体 (realpath) の隣を探す。以前は symlink の名前の隣を
# 探して「記録の無いメッシュ」(移行期 WARN、強制で FAIL) にしていた
src = make_run('symlink_src')
rc, out, started, gate = run_case(src)
launched(rc, out, started, gate, 'symlink の元の run (実体の隣に記録・品質の結果を書く)')
d = make_run('symlink_run', record=False)
os.unlink(os.path.join(d, 'mesh.h5'))
os.symlink(os.path.join(src, 'mesh.h5'), os.path.join(d, 'mesh.h5'))
before = sorted((f, os.path.getmtime(os.path.join(src, f))) for f in os.listdir(src))
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, 'メッシュを symlink で置いた run')
check(not [k for a_, k in items(gate, 'WARN') + items(gate, 'FAIL') if a_ == '来歴'] and '来歴の記録と一致' in (field(gate, 'mesh_hash') or ''),
      'symlink の run: 実体の隣の記録が見つかり来歴は PASS', str(items(gate, 'WARN')))
check((field(gate, 'quality') or '').startswith('PASS (再利用 あり'), 'symlink の run: 実体の隣の品質の結果を再利用', field(gate, 'quality'))
check(sorted((f, os.path.getmtime(os.path.join(src, f))) for f in os.listdir(src)) == before
      and not os.path.exists(os.path.join(d, 'mesh.h5.provenance.json')),
      'symlink の run の関門は実体のディレクトリに何も書かない (記録も写さない)')
other_real = os.path.join(TMP, 'symlink_other.h5')
write_mesh(other_real, shift=0.5)
os.unlink(os.path.join(d, 'mesh.h5'))
os.symlink(other_real, os.path.join(d, 'mesh.h5'))
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='enforce')
stopped(rc, out, started, gate, 'symlink の先を記録の無い別のメッシュに替えた (強制)')
check(('来歴', 'symlink_other.h5.provenance.json') in items(gate, 'FAIL'), '記録の無い実体を名指し', str(items(gate, 'FAIL')))
os.unlink(os.path.join(d, 'mesh.h5'))
os.symlink(os.path.join(src, 'mesh.h5'), os.path.join(d, 'mesh.h5'))
write_mesh(os.path.join(src, 'mesh.h5'), shift=0.25)         # 実体を別のメッシュに差し替える (記録は元のまま)
q_src, q_run = os.path.join(src, 'mesh.h5.quality.json'), os.path.join(d, 'mesh.h5.quality.json')
q_src_before, q_run_before = open(q_src).read(), open(q_run).read()     # q_run は上の symlink_other.h5 を測った結果
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'symlink の実体を別のメッシュに差し替えた (記録は元のまま)')
check(('来歴', '(provenance)') in items(gate, 'FAIL'), '記録と実体の不一致を名指しで FAIL', str(items(gate, 'FAIL')))
check((field(gate, 'quality') or '').startswith('PASS (再利用 なし') and open(q_src).read() == q_src_before
      and open(q_run).read() != q_run_before,
      '品質の結果は鍵が違えば実体の隣・run の中のどちらも使わず、測り直した結果は run の中 (symlink の名前の隣) にだけ書く',
      field(gate, 'quality'))

# ------------------------------------------------------------------------------------------------ メッシュの品質の負例

print('== メッシュの品質')
d = make_run('ar6000', mesh_kw={'first': (1.0 / 20) / 6000.0})
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'AR 6000 の壁の層 (緩和の上限 5000 を超える)')
check(('品質', 'mesh.h5') in items(gate, 'FAIL') and (field(gate, 'quality') or '').startswith('FAIL'), '品質の FAIL を名指し')
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
stopped(rc, out, started, gate, 'AR 6000 + warn (落とさない)')
d = make_run('ar2000_normal', mesh_kw={'first': (1.0 / 20) / 2000.0})
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, 'AR 2000 の壁法線の構造層 (緩和の条件を確かめて 5000 まで)')
d = make_run('no_bcond')
os.unlink(os.path.join(d, 'bcondConfig.yaml'))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'bcondConfig.yaml が無い (--bcond を渡せない)')
check(('品質', 'bcondConfig.yaml') in items(gate, 'FAIL'), '--bcond の欠落を名指しで FAIL')
d = make_run('bcond_type', declare=False)
open(os.path.join(d, 'bcondConfig.yaml'), 'w').write('wall: {physID: 3, type: wall_isothermal, Tw: 370.0}\n')
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'kind の無い bcondConfig (検査不能)')

# ------------------------------------------------------------------------------------------------ 適用条件 (次元・格子の種類)

print('== 適用条件 (次元・格子の種類・格子の node/cell; plan §5.1 #28)')
SERN2D, SERN3D = 'sern2d_finite_base_frozen_tp_sst', 'sern3d_frozen_tp_sst'
K_ISO = dict(BCOND, wall=(3, 'wall_isothermal'))
s2, s3 = cfg_flat(SERN2D, STAGE), cfg_flat(SERN3D, STAGE)


ok2 = make_run('sern2d_quad', flat=s2, kinds=K_ISO, profile=SERN2D)
rc, out, started, gate = run_case(ok2)
launched(rc, out, started, gate, 'SERN 2D の宣言 + 2D の四角形の格子 (正しい組)')
check(keys(gate, 'FAIL') == set() and items(gate, 'WARN') == [('バイナリ', FAKE_FORGE)],
      'SERN 2D の正しい組: FAIL 0、WARN は「対応の分からないバイナリ」だけ', str(items(gate, 'FAIL') + items(gate, 'WARN')))
check(field(gate, 'mesh_kind') == 'node; quad 100 (来歴の記録の成果物の事実)', 'LAUNCH_GATE.txt に格子の事実 (来歴の記録から)',
      field(gate, 'mesh_kind'))
dec = json.load(open(os.path.join(ok2, 'RECIPE.json')))
check({u['key'] for u in dec.get('unchecked_applies') or []} == {'applies.dimension', 'applies.mesh'},
      '宣言は格子を読まず、次元・格子の種類を照合していないと unchecked_applies に名指し (関門が照合する)',
      str(dec.get('unchecked_applies')))
ok3 = make_run('sern3d_hex', flat=s3, kinds=K_ISO, profile=SERN3D, mesh_kw={'cells': 'hex'})
rc, out, started, gate = run_case(ok3)
launched(rc, out, started, gate, 'SERN 3D の宣言 + 3D の六面体の格子 (正しい組)')
check(keys(gate, 'FAIL') == set() and field(gate, 'mesh_kind') == 'node; hex 200 (来歴の記録の成果物の事実)',
      'SERN 3D の正しい組: FAIL 0、格子の事実は六面体', field(gate, 'mesh_kind'))

# 2D/3D の取り違え (宣言は格子を読まないので書ける。関門が成果物から照合して止める)
d = make_run('sern2d_on_hex', flat=s2, kinds=K_ISO, mesh_kw={'cells': 'hex'}, profile=SERN2D)
rc_, r_, _ = RS.declare(d, SERN2D, STAGE, recipes=RC, write=False, mesh=MP.mesh_kind(os.path.join(d, 'mesh.h5')))
check(rc_ == 1 and 'applies.dimension' in {k for k, _ in r_['fails']},
      'declare に格子の事実 (mesh) を渡せば宣言の時点で止める (applies.dimension を名指し)', str(r_['fails']))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'SERN 2D の宣言 + 3D の六面体の格子')
fl_ = [ln for ln in gate.splitlines() if ln.startswith('FAIL') and 'applies.dimension' in ln]
check(len(fl_) == 1 and "2d_planar" in fl_[0] and "3d" in fl_[0], '2D/3D の取り違えを applies.dimension で名指し (プロファイル 2d_planar・格子 3d)',
      str(items(gate, 'FAIL')))
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
stopped(rc, out, started, gate, 'SERN 2D の宣言 + 3D の格子 + FORGE_LAUNCH_GATE=warn (適用条件の外は落とさない)')
d = make_run('sern3d_on_quad', flat=s3, kinds=K_ISO, profile=SERN3D)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'SERN 3D の宣言 + 2D の四角形の格子')
check('applies.dimension' in keys(gate, 'FAIL'), '3D/2D の取り違えを applies.dimension で名指し', str(items(gate, 'FAIL')))

# 構造/非構造の取り違え
d = make_run('sern2d_on_tri', flat=s2, kinds=K_ISO, mesh_kw={'cells': 'tri'}, profile=SERN2D)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'SERN 2D (structured) の宣言 + 三角形の格子')
fl_ = [ln for ln in gate.splitlines() if ln.startswith('FAIL') and 'applies.mesh' in ln]
check(len(fl_) == 1 and 'structured' in fl_[0] and 'unstructured' in fl_[0] and 'applies.dimension' not in keys(gate, 'FAIL'),
      '構造/非構造の取り違えを applies.mesh で名指し (次元は合う)', str(items(gate, 'FAIL')))

# 判定不能 (型を読めない格子) は一致扱いしない: 宣言した run は移行期でも FAIL、warn でも落とさない
d = make_run('sern2d_on_polygon', flat=s2, kinds=K_ISO, mesh_kw={'cells': 'polygon'}, record=False, profile=SERN2D)
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'SERN 2D の宣言 + 型を読めない格子 (判定不能)')
unc = [ln for ln in gate.splitlines() if ln.startswith('FAIL') and ('applies.dimension' in ln or 'applies.mesh' in ln)]
check(len(unc) == 2 and all('照合できない' in ln for ln in unc), '判定不能の次元・格子の種類を名指しで FAIL (一致扱いしない)', str(unc))
rc, out, started, gate = run_case(d, FORGE_LAUNCH_GATE='warn')
check({'applies.dimension', 'applies.mesh'} <= keys(gate, 'FAIL') and not started, '判定不能 + warn: 落とさない', str(items(gate, 'FAIL')))
d = make_run('exp_on_polygon', flat=s2, kinds=K_ISO, mesh_kw={'cells': 'polygon'}, record=False, declare=False)
assert RS.declare(d, SERN2D, STAGE, experiment=True, reason='実験: 型の試験', recipes=RC)[0] == 0
rc_, res_, _ = gate_inproc(d, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
check(('NOTE', '設定', 'applies.dimension') in [(lv, a, k) for lv, a, k, _ in res_.items]
      and ('FAIL', '設定', 'applies.dimension') not in [(lv, a, k) for lv, a, k, _ in res_.items],
      '実験の宣言 (プロファイルつき) では照合できない適用条件は NOTE', str(res_.items))

# 来歴の記録が無い・#28 より前の記録 (cell_types なし) は h5 から読む
d = make_run('sern2d_norecord', flat=s2, kinds=K_ISO, profile=SERN2D, record=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, 'SERN 2D の正しい組・記録の無いメッシュ (移行期)', want=('WARN',))
check(not {k for k in keys(gate, 'FAIL') if k.startswith('applies.')} and
      field(gate, 'mesh_kind') == 'node; quad 100 (h5 から読んだ; 来歴の記録が無い)', '記録が無ければ h5 から格子の事実を読む',
      field(gate, 'mesh_kind'))
d = make_run('sern2d_oldrecord', flat=s2, kinds=K_ISO, profile=SERN2D)
pj = os.path.join(d, 'mesh.h5.provenance.json')
rj = json.load(open(pj))
del rj['artifact']['cell_types']
json.dump(rj, open(pj, 'w'), ensure_ascii=False)
rc, out, started, gate = run_case(d)
launched(rc, out, started, gate, 'SERN 2D の正しい組・#28 より前の記録 (cell_types なし)')
check(field(gate, 'mesh_kind') == 'node; quad 100 (h5 から読んだ; 来歴の記録に cell_types が無い)',
      '記録に cell_types が無ければ h5 から読む', field(gate, 'mesh_kind'))

# codex result 段 M2 の再現: SERN 2D の正常な config に、プロファイルの適用条件だけを変える (以前は FAIL 0・WARN 0)
import copy                                                                                     # noqa: E402
_load = RS.load_recipes
for what, val, key in (('dimension', '3d', 'applies.dimension'), ('mesh', 'unstructured', 'applies.mesh'),
                       ('dimension', '4d', 'applies.dimension'), ('dimensions', '2d_planar', 'applies.dimensions')):
    rc2 = copy.deepcopy(RC)
    rc2['profiles'][SERN2D].setdefault('applies', {})[what] = val
    RS.load_recipes = lambda *a, _r=rc2, **k: _r                                                 # noqa: E731
    try:
        rc_, res_, _ = gate_inproc(ok2, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
    finally:
        RS.load_recipes = _load
    fk = [(a, k) for lv, a, k, _ in res_.items if lv == 'FAIL']
    check(rc_ == 1 and ('設定', key) in fk, f'プロファイルの適用条件だけ {what}: {val} → {key} を名指しで FAIL', str(fk))
rc_, res_, _ = gate_inproc(ok2, forge_bin=FAKE_FORGE, environ={}, dry_run=True)
check(rc_ == 0, '戻したプロファイルでは同じ run が通る', str(res_.items))
r_ = RS.check_config(RC, cfg_text(s2), SERN2D, STAGE, bcond=yaml.safe_load(bcond_text(K_ISO)))
check(not r_['fails'] and {k for k, _ in r_['unchecked']} == {'applies.dimension', 'applies.mesh'},
      'check_config に格子の事実を渡さなければ、次元・格子の種類は照合していないと名指し (黙って飛ばさない)', str(r_['unchecked']))
r_ = RS.check_config(RC, cfg_text(s2), SERN2D, STAGE, bcond=yaml.safe_load(bcond_text(K_ISO)),
                     mesh={'node': False, 'cell_types': {'quad': 100}})
check({k for k, _ in r_['fails']} == {'applies.discretization'} and not r_['unchecked'],
      'check_config: 格子そのものが cell (config は node) なら applies.discretization を名指しで FAIL', str(r_['fails']))

# ------------------------------------------------------------------------------------------------ モード・入力

print('== モード・入力')
rc, out, started, gate = run_case(good, FORGE_LAUNCH_GATE='wran')
stopped(rc, out, started, gate, 'FORGE_LAUNCH_GATE の未知の値 (綴り違い)')
d = make_run('cfg_list', declare=False)
open(os.path.join(d, 'solverConfig.yaml'), 'w').write('- 1\n- 2\n')
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'solverConfig.yaml が写像でない (検査不能)')
d = make_run('cfg_missing', declare=False)
os.unlink(os.path.join(d, 'solverConfig.yaml'))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'solverConfig.yaml が無い (検査不能)')
d = make_run('mesh_missing')
os.unlink(os.path.join(d, 'mesh.h5'))
rc, out, started, gate = run_case(d)
stopped(rc, out, started, gate, 'メッシュの h5 が無い')

print()
if fails:
    print(f'{fails} FAIL (出力を残した: {TMP})')
    sys.exit(1)
shutil.rmtree(TMP, ignore_errors=True)
print('ALL PASS')
