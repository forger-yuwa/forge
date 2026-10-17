#!/usr/bin/env python3
"""計算の後の判定 `judge_run.py` の性能の試験 (plan tooling-convergence-and-wall-resolution-gates §4.6-6・§5.1 #12 後半)。

**目標**: 45 万点・出力 80 枚の run で `judge_run.py <run>` (既定の監査の範囲: 残差の全履歴 + 最終の場; 準定常の snapshot の量は全出力を
読む) が **2 分以内**。条件を固定する (この台本の既定の引数が条件):

- 合成の run: node の構造格子 ni × nj = 1500 × 300 = 450,000 点 (MESH/COORD・CELLS/centCoords・PLANES/STRUCT ≈ 90 万面・壁 1 枚)、
  体積の出力 `res_<N>.h5` 80 枚 (N = 1000, 2000, …, 80000)。各出力は実際の run (case/45 run_0138 の res_6000: float32 の 23 量 +
  MESH/COORD float32・MESH/CONNE int32) と同じ形 (1 枚 約 56 MB、計 約 4.5 GB)、壁の出力 `res_wall_3_<N>.h5` 80 枚、
  `residual_history.csv` 16 万行 (8 万 step × outer_begin/outer_end、7 列)、`forge_run.log` 30 MB。
- 判定の仕様: 準定常は snapshot の machmax・pmax (全 80 枚を読む) と CSV の列 1 つ、壁解像は定数粘性の壁 1 枚 (適用)、収束は単段の履歴全体。
- 媒体: 作業ツリーのあるローカルのディスク (`--dir`、既定は一時ディレクトリ)。キャッシュ: **冷えた状態** は各ファイルに
  `posix_fadvise(POSIX_FADV_DONTNEED)` を掛けてページキャッシュから落としてから測る (root 不要。落ちたかは `--check-cold` の
  1 ファイルの読み直しの時間で確かめる)。落とせない環境では温かい状態として明記する。続けて温かい状態 (2 回目) も測る。

実行: python3 solver_density_cuda/tests/unit/perf_judge_run.py [--dir 置き場所] [--keep]   (生成 約 1〜2 分、判定 2 回)
終了コード: 冷えた状態 (落とせなければ温かい状態) の判定が 2 分以内なら 0、超えたら 1。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, '..', '..', 'tools'))
JUDGE = os.path.join(TOOLS, 'judge_run.py')
LIMIT_S = 120.0
CP, GAMMA, VISC = 1004.5, 1.4, 1.8e-5
FIELDS = ('P', 'T', 'Ux', 'Uy', 'Uz', 'Y0', 'Y1', 'h0', 'k', 'omega', 'ro', 'roK', 'roOmega', 'roUx', 'roUy', 'roUz', 'roY0', 'roY1',
          'roe', 'sonic', 'vis_lam', 'vis_turb', 'wall_dist')


def build(d, ni, nj, nout, step_int, log_mb):
    os.makedirs(d, exist_ok=True)
    n = ni * nj
    x = np.linspace(0.0, 1.0, ni)
    R = 0.1
    dr = 1.0e-5                                        # 第一層 (壁から 1 層目まで)
    yj = np.concatenate([[0.0], R - np.geomspace(0.999 * R, dr, nj - 2), [R]])   # 0 … R − dr, R (単調、壁で細かい)
    X, Y = np.meshgrid(x, yj, indexing='ij')
    coords = np.c_[X.ravel(), Y.ravel(), np.zeros(n)]
    nid = np.arange(n).reshape(ni, nj)
    a1, b1 = nid[:-1, :].ravel(), nid[1:, :].ravel()            # i 方向の面
    a2, b2 = nid[:, :-1].ravel(), nid[:, 1:].ravel()            # j 方向の面
    a = np.concatenate([a1, a2])
    b = np.concatenate([b1, b2])
    sv = np.concatenate([np.tile([1.0, 0.0, 0.0], (len(a1), 1)), np.tile([0.0, 1.0, 0.0], (len(a2), 1))])
    wall = nid[:, -1]
    dx = np.gradient(x)
    st_int = np.c_[np.full(len(a), 2), a, b, np.full(len(a), 2), a, b].ravel()
    st_bnd = np.c_[np.ones(ni, int), wall, np.ones(ni, int), wall].ravel()
    n_int = len(a)
    with h5py.File(os.path.join(d, 'mesh.h5'), 'w') as f:
        f.create_dataset('MESH/COORD', data=coords.ravel())
        f.create_dataset('CELLS/centCoords', data=coords.ravel())
        f.create_dataset('PLANES/STRUCT', data=np.concatenate([st_int, st_bnd]).astype(np.int32))
        f.create_dataset('PLANES/surfVect', data=np.concatenate([sv, np.c_[np.zeros(ni), dx, np.zeros(ni)]]).ravel())
        f.create_dataset('PLANES/surfArea', data=np.concatenate([np.ones(n_int), dx]))
        f.create_dataset('BCONDS/3/iCells', data=wall.astype(np.int32))
        f.create_dataset('BCONDS/3/iPlanes', data=(n_int + np.arange(ni)).astype(np.int32))
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as f:
        f.write('mesh: {discretization: "node", isAxisymmetric: 0, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}\n'
                'physProp: {thermalMethod: 0, cp: %r, gamma: %r, viscMethod: 0, visc: %r}\n'
                'turbulence: {model: "sst", wallTreatmentSST: 0}\n' % (CP, GAMMA, VISC))
    with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as f:
        f.write('inlet: {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 1.0e5, Tt: 300.0}}\n'
                'wall: {physID: 3, kind: wall, outputHDFflg: 1, ints: , floats: }\n')
    steps = [step_int * (k + 1) for k in range(nout)]
    with open(os.path.join(d, 'residual_history.csv'), 'w') as f:
        cols = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe', 'rms_roK', 'rms_roOmega']
        f.write('step,inner,phase,' + ','.join(cols) + '\n')
        ns = steps[-1]
        for s in range(ns):
            v = '%.6e' % (10.0 ** (-5.0 * s / (ns - 1)))
            row = ','.join([v] * len(cols))
            f.write('%d,-1,outer_begin,%s\n%d,-1,outer_end,%s\n' % (s, row, s, row))
    with open(os.path.join(d, 'forge_run.log'), 'w') as f:
        line = 'step %8d  dt 1.0e-07  cfl 4.0  rms_ro 1.0e-05  [synthetic log line padding ................................]\n'
        k = 0
        while f.tell() < log_mb * (1 << 20):
            f.write(line % k)
            k += 1
        f.write('Time = 1234.5 s (wall, %d steps, 15.4 ms/step)\n' % steps[-1])
    rng = np.random.default_rng(0)
    base = {k: np.full(n, v, np.float32) for k, v in (('P', 1.0e5), ('T', 300.0), ('Ux', 500.0), ('Uy', 0.0), ('Uz', 0.0),
                                                       ('Y0', 1.0), ('Y1', 0.0), ('k', 1.0), ('omega', 1e4), ('ro', 1.16),
                                                       ('sonic', 347.0), ('vis_lam', VISC), ('vis_turb', 0.0))}
    base['h0'] = (CP * base['T'] + 0.5 * base['Ux'] ** 2).astype(np.float32)
    base['wall_dist'] = (R - Y.ravel()).astype(np.float32)
    for k_, src in (('roK', 'k'), ('roOmega', 'omega'), ('roUx', 'Ux'), ('roUy', 'Uy'), ('roUz', 'Uz'), ('roY0', 'Y0'), ('roY1', 'Y1')):
        base[k_] = (base['ro'] * base[src]).astype(np.float32)
    base['roe'] = (base['P'] / (GAMMA - 1.0) + 0.5 * base['ro'] * base['Ux'] ** 2).astype(np.float32)
    quads = np.c_[nid[:-1, :-1].ravel(), nid[1:, :-1].ravel(), nid[1:, 1:].ravel(), nid[:-1, 1:].ravel()]
    conne = np.c_[np.full(len(quads), 5), quads].ravel().astype(np.int32)
    cf = coords.astype(np.float32).ravel()
    for s in steps:
        with h5py.File(os.path.join(d, 'res_%d.h5' % s), 'w') as f:
            f.create_dataset('MESH/COORD', data=cf)
            f.create_dataset('MESH/CONNE', data=conne)
            for k_ in FIELDS:
                v = base[k_]
                if k_ in ('P', 'Ux'):
                    v = v * np.float32(1.0 + 1e-6 * rng.standard_normal())
                f.create_dataset('VALUE/' + k_, data=v)
    # 壁の出力 (壁の節点、逆順; 接線の壁応力は y₁⁺ = 0.5 になる値)
    perm = wall[::-1]
    y1 = float(R - yj[-2])
    tau = (0.5 * VISC / y1) ** 2 / 1.16
    wc = np.c_[2 * np.ones(ni - 1, int), 2 * np.ones(ni - 1, int), np.arange(ni - 1), np.arange(1, ni)].ravel()
    for s in steps:
        with h5py.File(os.path.join(d, 'res_wall_3_%d.h5' % s), 'w') as f:
            f.create_dataset('MESH/COORD', data=coords[perm].ravel())
            f.create_dataset('MESH/CONNE', data=wc.astype(np.int64))
            for k_, v in (('ro', 1.16), ('Ts', 300.0), ('Ps', 1.0e5), ('twall_x', tau), ('twall_y', 0.0), ('twall_z', 0.0)):
                f.create_dataset('VALUE/' + k_, data=np.full(ni, v))
    with open(os.path.join(d, 'force_history.csv'), 'w') as f:
        f.write('step,C_T\n' + ''.join('%d,%r\n' % (s, 0.95) for s in steps))
    spec = {'judge_spec_version': 1, 'source': 'perf_judge_run.py', 'quasisteady': {'quantities': [
        {'name': n_, 'definition': '性能の試験', 'unit': '-', 'region': '全域', 'extractor': {'kind': 'snapshot', 'quantity': n_}}
        for n_ in ('machmax', 'pmax')] + [
        {'name': 'C_T', 'definition': '性能の試験', 'unit': '-', 'region': '壁',
         'extractor': {'kind': 'csv_column', 'file': 'force_history.csv', 'column': 'C_T'}}]}}
    with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
        json.dump(spec, f)


def drop_cache(d):
    """run の全ファイルをページキャッシュから落とす (posix_fadvise DONTNEED)。戻り 落としたファイルの数 (使えなければ None)。"""
    if not hasattr(os, 'posix_fadvise'):
        return None
    k = 0
    for fn in os.listdir(d):
        p = os.path.join(d, fn)
        if not os.path.isfile(p):
            continue
        fd = os.open(p, os.O_RDONLY)
        try:
            os.fsync(fd)
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            k += 1
        finally:
            os.close(fd)
    return k


def read_time(p):
    t = time.time()
    with open(p, 'rb') as f:
        while f.read(1 << 24):
            pass
    return time.time() - t


def run_judge(d):
    t = time.time()
    r = subprocess.run([sys.executable, JUDGE, d, '--solver-rc', '0'], capture_output=True, text=True)
    dt = time.time() - t
    v = json.load(open(os.path.join(d, 'RUN_VERDICT.json')))
    return dt, r.returncode, v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=None, help='合成の run の置き場所 (既定: 一時ディレクトリ)')
    ap.add_argument('--ni', type=int, default=1500)
    ap.add_argument('--nj', type=int, default=300)
    ap.add_argument('--nout', type=int, default=80)
    ap.add_argument('--step-int', type=int, default=1000)
    ap.add_argument('--log-mb', type=int, default=30)
    ap.add_argument('--keep', action='store_true', help='合成の run を消さない')
    a = ap.parse_args()
    base = a.dir or tempfile.mkdtemp(prefix='perf_judge_run_')
    d = os.path.join(base, 'run_perf')
    t0 = time.time()
    build(d, a.ni, a.nj, a.nout, a.step_int, a.log_mb)
    size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
    print('合成の run: %s  (%d 点・出力 %d 枚・計 %.2f GB、生成 %.0f s)' % (d, a.ni * a.nj, a.nout, size / 1e9, time.time() - t0), flush=True)
    probe = os.path.join(d, 'res_%d.h5' % (a.step_int * (a.nout // 2)))
    t_warm_read = read_time(probe)
    n_drop = drop_cache(d)
    t_cold_read = read_time(probe) if n_drop is not None else None
    if n_drop is not None:
        drop_cache(d)
    # 落とせたと確かめられるのは読み直しが温かい読みの 3 倍以上かかったときだけ (WSL2 では仮想ディスクがホスト側でキャッシュされ、
    # ゲストのページキャッシュを落としても冷えない; そのときは温かい状態として書く)
    cold = n_drop is not None and t_cold_read is not None and t_cold_read >= 3.0 * t_warm_read
    dt1, rc1, v1 = run_judge(d)
    dt2, rc2, v2 = run_judge(d)
    mb = os.path.getsize(probe) / 1e6
    state = ('冷えた状態 (posix_fadvise DONTNEED で落とした; 1 枚 %.0f MB の読み直し %.3f s vs 温かい %.3f s)' % (mb, t_cold_read, t_warm_read)
             if cold else '温かい状態 (ページキャッシュを落とせたと確かめられない: 1 枚 %.0f MB の読み直し %s s vs 温かい %.3f s)'
             % (mb, ('%.3f' % t_cold_read) if t_cold_read is not None else '—', t_warm_read))
    print('1 回目 (%s): %.1f s  (終了 %d、総合 %s)' % (state, dt1, rc1, v1['overall']))
    print('2 回目 (温かい状態): %.1f s  (終了 %d、総合 %s)' % (dt2, rc2, v2['overall']))
    for k, c in v1['checks'].items():
        print('  %-16s %-15s %s' % (k, c['status'], c['reason'][:150]))
    for nm, q in (v1['checks'].get('quasisteady', {}).get('quantities') or {}).items():
        print('    %-10s 点数 %s (%s)' % (nm, q.get('n'), q['extractor']['kind']))
    ok = dt1 <= LIMIT_S
    print('VERDICT: %s (1 回目 %.1f s %s %d s; 条件: %d 点・出力 %d 枚・%.2f GB・%s)'
          % ('PASS' if ok else 'FAIL', dt1, '≤' if ok else '>', LIMIT_S, a.ni * a.nj, a.nout, size / 1e9,
             '冷えた状態' if cold else '温かい状態'))
    if not a.keep:
        shutil.rmtree(base if a.dir is None else d, ignore_errors=True)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
