#!/usr/bin/env python3
"""計算の後の判定 `judge_run.py` と `run_case.sh` の結線の試験 (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12 前半)。

合成の run (forge なし) を一時ディレクトリに作って `judge_run.py` を回し、`RUN_VERDICT.json` の各項目と総合を見る。
本物の forge は起動しない (`run_case.sh` の通しは偽の forge を `FORGE_BIN` で渡す)。見るもの:

- 真理値表 (`_combine` の全組 + 実際の run で FAIL と判定不能の併存が FAIL になること)。
- 欠損: 段 (manifest の段の履歴が無い)・列 (系列の CSV の列が無い)・量 (抽出が適用不能)・壁 (最後の step の壁の出力が無い)・
  時刻 (壁の出力が古い step だけ、系列の step の重複・逆行)。空集合 (区間に点が無い)。窓の不足 (点数・末尾窓の step 幅)。
- 古い判定: 入力 (残差の履歴・仕様) が変わったら `--check-fresh` が 1。
- 定数粘性 (`viscMethod 0`) の壁で壁解像が適用になる (`wall_resolution_applies`) こと、仕様の applies と実効の食い違い。
- solver の終了コード (異常終了・監視の停止・不明)、非有限 (残差・最終の場・`res_nan_*`・`--full-audit`)、物理性、停止規則の状態の対応、
  収束の mode (非定常・単段と確かめられない・not_applicable・from_floor)、仕様の誤り、判定器の例外 (終了 3)。
- `run_case.sh` の通し: `RUN_VERDICT.*` ができる・関門で未起動なら判定しない・`LAUNCH_ECHO.txt` と `residual_history.png` が残る・
  終了コードは solver のまま・段階起動 (manifest / `FORGE_JUDGE_DEFER=1`) は `state: deferred`・`FORGE_STRICT_CONVERGENCE` は従来どおり。

実行: python3 solver_density_cuda/tests/unit/test_judge_run.py   (約 1 分、全部合格なら一時ディレクトリを消す)
"""
import contextlib
import io
import json
import os
import shutil
import socket
import stat
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.abspath(os.path.join(HERE, '..', '..', 'tools'))
sys.path.insert(0, TOOLS)
import judge_run as JR                    # noqa: E402
import stage_manifest as SM               # noqa: E402
import test_gate_bad_input as TG          # noqa: E402  (M6 と同じ壁の合成の run: make_wall_run・wall_cases)

JUDGE = os.path.join(TOOLS, 'judge_run.py')
RUN_CASE = os.path.join(TOOLS, 'run_case.sh')
TMP = tempfile.mkdtemp(prefix='judge_run_test_')
fails = 0


def check(cond, msg, detail=''):
    global fails
    print(('ok   ' if cond else 'FAIL ') + msg + (f'  [{detail}]' if detail and not cond else ''), flush=True)
    if not cond:
        fails += 1


# ------------------------------------------------------------------------------------------------ 合成の run

COLS = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe']
S_LAST = 2000


def write_hist(path, n=200, d0=0.0, d1=-5.0, nan_at=None, steps=None):
    """残差の履歴 (outer_end の行だけ)。d0 → d1 桁で落ちる (stalled にするなら d1 = d0)。"""
    with open(path, 'w') as f:
        f.write('step,inner,phase,' + ','.join(COLS) + '\n')
        for i in range(n):
            v = 10.0 ** (d0 + (d1 - d0) * i / (n - 1))
            vals = ['nan' if (nan_at is not None and i == nan_at) else '%.6e' % v for _ in COLS]
            f.write('%d,-1,outer_end,%s\n' % (i if steps is None else steps[i], ','.join(vals)))


CP, GAMMA = 1004.5, 1.4
H0_OK = CP * 300.0 + 0.5 * 10.0 ** 2          # h0 = c_p T + u²/2 (T 300 K・U 10 m/s): 全温 > 静温 → 全圧 > 静圧


def write_res(path, n, ro=1.2, nan_field=None, h0=H0_OK):
    """体積の出力。h0 (全エンタルピー、静圧 ≤ 全圧の判定が total_quantities で読む) は None なら書かない。"""
    with h5py.File(path, 'w') as f:
        for k, val in (('ro', ro), ('roUx', 10.0), ('roUy', 0.0), ('roUz', 0.0), ('roe', 2.5e5), ('P', 1.0e5), ('T', 300.0),
                       ('Ux', 10.0), ('Uy', 0.0), ('Uz', 0.0), ('sonic', 340.0)) + ((('h0', h0),) if h0 is not None else ()):
            a = np.full(n, val, float)
            if nan_field == k:
                a[0] = np.nan
            f.create_dataset('VALUE/' + k, data=a)


_n = [0]


def wall_run(tag, case='control', steps=None, hist=True, res=(1000, S_LAST), spec=None, solver_log=True, outlet=False):
    """2 枚の no-slip 壁 (定数粘性 viscMethod 0) の node の run (M6 と同じ make_wall_run) に、残差の履歴・体積の出力・
    forge_run.log・仕様を足す。"""
    _n[0] += 1
    d = os.path.join(TMP, '%02d_%s' % (_n[0], tag))
    walls = TG.wall_cases()[case]
    TG.make_wall_run(d, walls, steps if steps is not None else {'wall_a': [1000, S_LAST], 'wall_b': [1000, S_LAST]},
                     outlet=outlet)
    cp_ = os.path.join(d, 'solverConfig.yaml')          # 静圧 ≤ 全圧 (total_quantities の CPG) に要る熱物性を足す
    _t = open(cp_).read()
    assert 'physProp: {' in _t, _t
    open(cp_, 'w').write(_t.replace('physProp: {', 'physProp: {thermalMethod: 0, cp: %r, gamma: %r, ' % (CP, GAMMA), 1))
    with h5py.File(os.path.join(d, 'mesh.h5'), 'r') as f:
        n = len(f['MESH/COORD']) // 3
    if hist:
        write_hist(os.path.join(d, 'residual_history.csv'))
    for s in res:
        write_res(os.path.join(d, 'res_%d.h5' % s), n)
    if solver_log:
        with open(os.path.join(d, 'forge_run.log'), 'w') as f:
            f.write('[config] node scheme\nstep 2000 ...\nTime = 1.234 s (wall, 2000 steps, 0.62 ms/step)\n')
    if spec is not None:
        with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
            json.dump(spec, f)
    return d


QS_NA = {'judge_spec_version': 1, 'quasisteady': {'not_applicable': '試験: 報告する派生量が無い'}}


def judge(d, *args, env=None):
    e = dict(os.environ)
    e.pop('FORGE_JUDGE_DEFER', None)
    e.update(env or {})
    r = subprocess.run([sys.executable, JUDGE, d] + list(args), capture_output=True, text=True, env=e)
    v = None
    p = os.path.join(d, 'RUN_VERDICT.json')
    if os.path.exists(p):
        v = json.load(open(p))
    return r.returncode, v, r.stdout + r.stderr


def st(v, k):
    return (v or {}).get('checks', {}).get(k, {}).get('status')


def why(v, k):
    return (v or {}).get('checks', {}).get(k, {}).get('reason', '')


# ------------------------------------------------------------------------------------------------ 真理値表

print('== 真理値表 (_combine の全組)')
S = JR.STATUSES
bad = []
for a in S:
    for b in S:
        for c in S:
            got = JR._combine([a, b, c])
            want = ('FAIL' if 'FAIL' in (a, b, c) else 'INDETERMINATE' if 'INDETERMINATE' in (a, b, c)
                    else 'PASS' if 'PASS' in (a, b, c) else 'NOT_APPLICABLE')
            if got != want:
                bad.append((a, b, c, got))
check(not bad, '3 項目の 64 組: FAIL > INDETERMINATE > PASS > NOT_APPLICABLE (既知の FAIL を判定不能で隠さない)', bad)
check(JR._combine([]) == 'NOT_APPLICABLE' and JR.EXIT_OF['NOT_APPLICABLE'] == 2, '空・全部適用外は NOT_APPLICABLE で終了 2 (合格でない)')

# ------------------------------------------------------------------------------------------------ 正例と壁解像

print('== 正例 (2 枚の no-slip 壁・定数粘性・収束・仕様で準定常は適用外)')
d = wall_run('pass', spec=QS_NA)
rc, v, out = judge(d, '--solver-rc', '0')
check(rc == 0 and v['overall'] == 'PASS' and v['state'] == 'final', '総合 PASS・終了 0・state final', f'rc {rc} {out[-800:]}')
check(st(v, 'wall_resolution') == 'PASS' and v['checks']['wall_resolution']['step'] == S_LAST
      and sorted(v['checks']['wall_resolution']['rows']) == ['wall_a', 'wall_b'],
      '壁解像: 定数粘性 (viscMethod 0) の 2 壁を最後の step・面積の割合で PASS', why(v, 'wall_resolution'))
check('--weight' in v['checks']['wall_resolution']['cmd'] and 'area' in v['checks']['wall_resolution']['cmd'], '壁解像は --weight area で回す')
check(st(v, 'convergence') == 'PASS' and 'whole' in why(v, 'convergence') and 'single_stage' in v['checks']['convergence'],
      '収束: manifest 無し・単段と確かめて履歴全体で PASS', why(v, 'convergence'))
check(st(v, 'stop_rule') == 'NOT_APPLICABLE' and v['checks']['stop_rule']['raw'] == 'NOT_USED', '停止規則: NOT_USED → 適用外')
check(st(v, 'quasisteady') == 'NOT_APPLICABLE', '準定常: 仕様の not_applicable')
txt = open(os.path.join(d, 'RUN_VERDICT.txt'), encoding='utf-8').read()
check(txt.splitlines()[0].startswith('OVERALL: PASS'), 'RUN_VERDICT.txt の先頭の 1 行が総合', txt.splitlines()[0])
for k in ('tools', 'inputs', 'execution'):
    check(bool(v.get(k)), f'RUN_VERDICT.json に {k} (道具の sha・入力の識別・実行の識別子)')
check(v['tools'].get('check_convergence.py', '').__len__() == 64 and v['inputs'].get('residual_history.csv', '').startswith('sha256:')
      and v['execution'].get('run_id'), '道具の sha256・入力の sha256・run_id が入る')
check(not [f for f in os.listdir(d) if f.endswith('.tmp')], '一時ファイルが残らない (rename で書いた)')

print('== 壁解像の適用 (仕様側の関数)')
cfg0 = {'physProp': {'viscMethod': 0, 'visc': 1.8e-5}}
bc = {'wall': {'physID': 3, 'kind': 'wall'}, 'iso': {'physID': 5, 'kind': 'wall_isothermal'}, 'out': {'physID': 2, 'kind': 'outlet_statPress'},
      'sym': {'physID': 4, 'kind': 'slip'}}
a, w_, ws = JR.wall_resolution_applies(cfg0, bc)
check(a is True and ws == {'wall': 3, 'iso': 5} and '定数粘性' in w_, '定数粘性 (viscMethod 0・visc 1.8e-5) + no-slip 壁 → 適用 (slip・出口は壁に入れない)', w_)
check(JR.wall_resolution_applies({'physProp': {'viscMethod': 0, 'visc': 0.0}}, bc)[0] is False, 'viscMethod 0・visc 0 (Euler) → 適用外')
check(JR.wall_resolution_applies({'physProp': {'viscMethod': 1, 'visc': 0.0}}, bc)[0] is True, 'Sutherland (viscMethod 1) → 適用')
check(JR.wall_resolution_applies(cfg0, {'sym': {'physID': 4, 'kind': 'slip'}})[0] is False, 'no-slip 壁が無い → 適用外')
check(JR.wall_resolution_applies({'physProp': {}}, bc)[0] is None, 'viscMethod を読めない → 決められない (None)')
check(JR.wall_resolution_applies(cfg0, None)[0] is None, 'bcondConfig を読めない → 決められない (None)')

print('== 壁の欠損・時刻の欠損')
d = wall_run('wall_b_old', spec=QS_NA, steps={'wall_a': [1000, S_LAST], 'wall_b': [1000]})
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'wall_resolution') == 'INDETERMINATE' and 'wall_b' in why(v, 'wall_resolution') and v['overall'] == 'INDETERMINATE' and rc == 2,
      '片方の壁の最後の step の出力が無い (古い step だけ) → 判定不能 (最新に戻らない)・総合 INDETERMINATE・終了 2', why(v, 'wall_resolution'))
d = wall_run('wall_b_absent', spec=QS_NA, steps={'wall_a': [1000, S_LAST]})
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'wall_resolution') == 'INDETERMINATE' and 'wall_b' in why(v, 'wall_resolution'),
      '壁 wall_b の出力が 1 つも無い (outputHDFflg 0 相当) → 判定不能', why(v, 'wall_resolution'))
d = wall_run('area_fail', case='inv_area_fail', spec=QS_NA)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'wall_resolution') == 'FAIL' and v['overall'] == 'FAIL' and rc == 1, '面積の割合で y₁⁺ > 1 が許容を超える → FAIL・終了 1', why(v, 'wall_resolution'))
d = wall_run('spec_mismatch', spec=dict(QS_NA, wall_resolution={'applies': False, 'reason': '試験'}))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'wall_resolution') == 'INDETERMINATE' and '実効' in why(v, 'wall_resolution'),
      '仕様の applies false が実効 (定数粘性の壁) と違う → 判定不能 (仕様で壁解像を外せない)', why(v, 'wall_resolution'))
d = wall_run('inviscid', spec=QS_NA)
cfgp = os.path.join(d, 'solverConfig.yaml')
_t = open(cfgp).read()
assert 'visc: 1.8e-05' in _t, _t
open(cfgp, 'w').write(_t.replace('visc: 1.8e-05', 'visc: 0.0'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'wall_resolution') == 'NOT_APPLICABLE' and v['overall'] == 'PASS', 'Euler (visc 0) → 壁解像は適用外で総合 PASS', why(v, 'wall_resolution'))

# ------------------------------------------------------------------------------------------------ 準定常

print('== 準定常 (仕様の量だけ・要求と評価の完全一致)')


def qspec(*items):
    return {'judge_spec_version': 1, 'quasisteady': {'quantities': list(items)}}


def qitem(name, extractor, **kw):
    it = {'name': name, 'definition': '試験の量', 'unit': '-', 'region': '全体', 'extractor': extractor}
    it.update(kw)
    return it


def write_series(d, rows, cols=('step', 'C_T')):
    with open(os.path.join(d, 'force_history.csv'), 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join(str(x) for x in r) + '\n')


CSV_CT = {'kind': 'csv_column', 'file': 'force_history.csv', 'column': 'C_T'}
d = wall_run('qs_none')
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '対象未設定' in why(v, 'quasisteady') and v['overall'] == 'INDETERMINATE',
      '仕様も judge.yaml も無い → 準定常は「対象未設定」= 判定不能 (黙って省かない)', why(v, 'quasisteady'))

d = wall_run('qs_steady', spec=qspec(qitem('C_T', CSV_CT)))
write_series(d, [(s, 0.95) for s in range(0, 2001, 100)])
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'PASS' and v['checks']['quasisteady']['quantities']['C_T']['verdict'] == 'STEADY' and v['overall'] == 'PASS',
      'CSV の列 (一定) → STEADY・総合 PASS', why(v, 'quasisteady'))
write_series(d, [(s, 0.95 + 1e-4 * s) for s in range(0, 2001, 100)])
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'FAIL' and v['checks']['quasisteady']['quantities']['C_T']['verdict'] == 'DRIFTING' and v['overall'] == 'FAIL',
      'CSV の列 (ドリフト) → DRIFTING = FAIL', why(v, 'quasisteady'))
write_series(d, [(s, 0.95) for s in range(0, 2001, 100)], cols=('step', 'C_L'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '列' in why(v, 'quasisteady'), '系列の CSV に要求の列が無い → 判定不能 (列の欠損)', why(v, 'quasisteady'))
write_series(d, [(1, 0.95)] * 21)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '単調' in why(v, 'quasisteady'), 'step が全行同じ (時間が進まない)・値一定 → 判定不能', why(v, 'quasisteady'))
write_series(d, [(s, 0.95) for s in list(range(0, 1001, 100)) + list(range(500, 1501, 100))])
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '逆行' in why(v, 'quasisteady'), 'step が逆行する (継ぎ足しの系列) → 判定不能', why(v, 'quasisteady'))
write_series(d, [(s, 0.95 if s != 1500 else float('nan')) for s in range(0, 2001, 100)])
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'FAIL' and 'NONFINITE' in why(v, 'quasisteady'), '系列に NaN → NONFINITE = FAIL', why(v, 'quasisteady'))
write_series(d, [(s, 0.95) for s in range(0, 2001, 100)])
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('C_T', CSV_CT, segment={'step_min': 5000})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '空集合' in why(v, 'quasisteady'), '区間 [5000, ] に点が無い → 判定不能 (空集合)', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('C_T', CSV_CT, window={'min_points': 30})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '窓の不足' in why(v, 'quasisteady'), '点数 21 < min_points 30 → 判定不能 (窓の不足)', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('C_T', CSV_CT, window={'min_tail_span': 20000})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and 'min_tail_span' in why(v, 'quasisteady'),
      '末尾窓の step 幅 < min_tail_span (停止規則の窓 2 万 step 相当) → 判定不能', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('C_T', CSV_CT), qitem('C_L', dict(CSV_CT, column='C_L'), required=False)), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'PASS' and v['checks']['quasisteady']['quantities']['C_L']['status'] == 'INDETERMINATE'
      and '参考' in why(v, 'quasisteady'), '必須でない量の判定不能は参考 (総合に入れないが出す)', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('asym', {'kind': 'snapshot', 'quantity': 'asym'})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'INDETERMINATE' and '適用不能' in why(v, 'quasisteady'),
      'snapshot の asym (上下対称でないメッシュ) → 判定不能 (量の欠損; 評価量ゼロで STEADY にしない)', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('machmax', {'kind': 'snapshot', 'quantity': 'machmax'}, window={'min_points': 2})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') in ('PASS', 'FAIL') and v['checks']['quasisteady']['quantities']['machmax'].get('n') == 2,
      'snapshot の machmax は全 res_<N>.h5 (2 枚) から抽出して判定する', why(v, 'quasisteady'))

print('== 仕様の誤り・置き場所')
d = wall_run('spec_bad', spec={'judge_spec_version': 1, 'quasisteady': {'not_applicable': 'x'}, 'convergnce': {}})
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['error'] and st(v, 'quasisteady') == 'INDETERMINATE' and st(v, 'convergence') == 'INDETERMINATE'
      and st(v, 'wall_resolution') == 'INDETERMINATE' and v['overall'] == 'INDETERMINATE',
      '仕様の未知のキー (綴り違い) → 仕様に依る項目はすべて判定不能', v['spec'].get('error'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump({'quasisteady': {'not_applicable': 'x'}}, f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['error'] and 'judge_spec_version' in v['spec']['error'], 'JUDGE_SPEC.json に版が無い → 仕様の誤り', v['spec'].get('error'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['error'] and st(v, 'quasisteady') == 'INDETERMINATE', 'quantities が空 (空集合) → 仕様の誤り・判定不能', v['spec'].get('error'))
os.unlink(os.path.join(d, 'JUDGE_SPEC.json'))
case_yaml = os.path.join(os.path.dirname(d), 'judge.yaml')
with open(case_yaml, 'w') as f:
    f.write('quasisteady: {not_applicable: case の仕様}\n')
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['source'] == 'judge.yaml (case)' and st(v, 'quasisteady') == 'NOT_APPLICABLE', 'case の judge.yaml を読む (run に仕様が無いとき)')
with open(os.path.join(d, 'judge.yaml'), 'w') as f:
    f.write('quasisteady: {not_applicable: run の仕様}\nquasisteady: {not_applicable: 重複}\n')
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['source'] == 'judge.yaml (run)' and v['spec']['error'], 'run の judge.yaml が優先・重複キーは仕様の誤り', v['spec'].get('error'))
os.unlink(case_yaml)
with open(os.path.join(d, 'RUN_PROVENANCE.txt'), 'w') as f:
    f.write('date        : %s\n' % JR._now())
time.sleep(1.1)
with open(os.path.join(d, 'judge.yaml'), 'w') as f:
    f.write('quasisteady: {not_applicable: 後から書いた}\n')
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['written_after_start'] is True and '開始より後' in open(os.path.join(d, 'RUN_VERDICT.txt'), encoding='utf-8').read(),
      'run の開始より後に書いた judge.yaml は written_after_start を記録する')

# ------------------------------------------------------------------------------------------------ 収束

print('== 収束 (区間・単段の確かめ・mode)')
d = wall_run('stalled', spec=QS_NA, steps={'wall_a': [1000, S_LAST]})
write_hist(os.path.join(d, 'residual_history.csv'), d0=-2.0, d1=-2.0)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'FAIL' and st(v, 'wall_resolution') == 'INDETERMINATE' and v['overall'] == 'FAIL' and rc == 1,
      '収束 FAIL (停滞) と壁の判定不能の併存 → 総合 FAIL (既知の FAIL を隠さない)', v['overall_reason'])
d = wall_run('stage_files', spec=QS_NA)
shutil.copy(os.path.join(d, 'residual_history.csv'), os.path.join(d, 'residual_history_S1.csv'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'INDETERMINATE' and '単段と確かめられない' in why(v, 'convergence'),
      '段の履歴があるのに manifest が無い → 単段と確かめられず判定不能', why(v, 'convergence'))
d = wall_run('concat', spec=QS_NA)
write_hist(os.path.join(d, 'residual_history.csv'), steps=list(range(100)) + list(range(100)))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'INDETERMINATE' and '単調' in why(v, 'convergence'), '残差の履歴の step が戻る (継ぎ足し) → 判定不能', why(v, 'convergence'))
d = wall_run('unsteady', spec=QS_NA)
with open(cfgp := os.path.join(d, 'solverConfig.yaml'), 'a') as f:
    f.write('time: {unsteady: 1, dualTime: 1}\n')
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'INDETERMINATE' and '非定常' in why(v, 'convergence'), '非定常の run で仕様が無い → 収束は判定不能', why(v, 'convergence'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(dict(QS_NA, convergence={'mode': 'not_applicable', 'reason': '非定常 LES'}), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'NOT_APPLICABLE' and '非定常 LES' in why(v, 'convergence'), '仕様の not_applicable (理由つき) → 適用外')

print('== 段の欠損 (manifest の最後の区間)')
cfg_txt = 'turbulence: {model: "sst", wallTreatmentSST: 0}\nspace: {convMethod: 1, limiter: 2}\n'


def staged(tag, main_hist=True):
    d = wall_run(tag, spec=QS_NA)
    m = SM.StageManifest(d)
    m.add('S1', cfg_txt, 'wall: {}\n', history='residual_history_S1.csv')
    m.add('main', cfg_txt, 'wall: {}\n', history='residual_history_main.csv')
    m.write()
    write_hist(os.path.join(d, 'residual_history_S1.csv'), d0=0.0, d1=-2.0)
    write_hist(os.path.join(d, 'residual_history.csv'), d0=-2.0, d1=-5.0)
    if main_hist:
        shutil.copy(os.path.join(d, 'residual_history.csv'), os.path.join(d, 'residual_history_main.csv'))
    return d


d = staged('seg_ok')
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'PASS' and 'segment' in why(v, 'convergence')
      and any('S1 -> main' in s for s in v['checks']['convergence']['segment_lines']), 'manifest あり → --segment (S1 -> main) で PASS',
      why(v, 'convergence'))
d = staged('seg_missing', main_hist=False)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'INDETERMINATE' and st(v, 'nonfinite') == 'INDETERMINATE',
      '本段の段の履歴が無い → 収束も非有限の監査も判定不能 (前段だけで合格しない)', why(v, 'convergence'))

print('== from_floor')
ref = wall_run('ref', spec=QS_NA)
d = wall_run('cont', spec=dict(QS_NA, convergence={'mode': 'from_floor', 'reference': '../' + os.path.basename(ref)}))
write_hist(os.path.join(d, 'residual_history.csv'), d0=-5.0, d1=-5.0)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'convergence') == 'PASS' and '--from-floor' in v['checks']['convergence']['cmd'], '継続の run: 仕様の from_floor で参照の床に留まれば PASS',
      why(v, 'convergence'))

# ------------------------------------------------------------------------------------------------ solver・非有限・物理性・停止規則

print('== solver の終了コード')
d = wall_run('rc', spec=QS_NA)
rc, v, _ = judge(d, '--solver-rc', '1')
check(st(v, 'solver') == 'FAIL' and v['overall'] == 'FAIL', 'forge の異常終了 (終了 1) → FAIL')
rc, v, _ = judge(d)
check(st(v, 'solver') == 'FAIL' and '前の判定' in v['checks']['solver']['rc_source'],
      '--solver-rc なしの再判定は前の判定の終了コードを引き継ぐ (同じ forge_run.log; 既知の FAIL を消さない)', v['checks']['solver'].get('rc_source'))
os.unlink(os.path.join(d, 'RUN_VERDICT.json'))
rc, v, _ = judge(d)
check(st(v, 'solver') == 'PASS' and 'Time =' in v['checks']['solver']['rc_source'], '記録が無ければ forge_run.log の正常終了の行で PASS')
with open(os.path.join(d, 'forge_run.log'), 'w') as f:
    f.write('[config] node scheme\nstep 100 ...\n')
os.unlink(os.path.join(d, 'RUN_VERDICT.json'))
rc, v, _ = judge(d)
check(st(v, 'solver') == 'INDETERMINATE' and v['overall'] == 'INDETERMINATE', '終了コードが分からない (正常終了の行も無い) → 判定不能')
with open(os.path.join(d, 'STOP_RULE_STATE.json'), 'w') as f:
    json.dump({'version': 1, 'status': 'REACHED', 'stop_step': S_LAST, 'forge_rc': 143}, f)
rc, v, _ = judge(d, '--solver-rc', '143')
check(st(v, 'solver') == 'PASS' and v['checks']['solver'].get('monitor_stop') and st(v, 'stop_rule') == 'PASS',
      '停止規則の見張りが止めた (143・REACHED) → 正常な監視の停止 (solver PASS)', why(v, 'solver'))
os.unlink(os.path.join(d, 'STOP_RULE_STATE.json'))
rc, v, _ = judge(d, '--solver-rc', '143')
check(st(v, 'solver') == 'FAIL', '停止の記録が無い 143 → 異常終了 (FAIL)')

print('== 停止規則の状態 (そのまま記録)')
want = {'REACHED': 'PASS', 'CENSORED': 'FAIL', 'DIVERGED': 'FAIL', 'EXEC_ERROR': 'FAIL', 'DATA_ERROR': 'INDETERMINATE',
        'UNAVAILABLE': 'NOT_APPLICABLE', 'NOT_USED': 'NOT_APPLICABLE'}
got = {}
for s_, w in want.items():
    with open(os.path.join(d, 'STOP_RULE_STATE.json'), 'w') as f:
        json.dump({'version': 1, 'status': s_}, f)
    got[s_] = (JR.check_stop_rule(JR.read_stop_state(d))['status'], JR.read_stop_state(d)['raw'])
check(all(got[s_] == (w, s_) for s_, w in want.items()), '7 状態 → 項目の状態の対応、raw は状態をそのまま', got)
os.unlink(os.path.join(d, 'STOP_RULE_STATE.json'))
with open(os.path.join(d, 'stop_rule.json'), 'w') as f:
    json.dump({'status': 'RUNNING'}, f)
check(JR.check_stop_rule(JR.read_stop_state(d))['status'] == 'INDETERMINATE', '見張りの途中の記録 (stop_rule.json の RUNNING) → 判定不能')
os.unlink(os.path.join(d, 'stop_rule.json'))
with open(os.path.join(d, 'prepare_info.json'), 'w') as f:
    json.dump({'stop_rule': {'requested': True}}, f)
s_ = JR.read_stop_state(d)
check(s_['raw'] == 'MISSING' and JR.check_stop_rule(s_)['status'] == 'INDETERMINATE', 'prepare_info が要求したのに状態が無い → MISSING (判定不能)')

print('== 非有限・物理性')
d = wall_run('nan_hist', spec=QS_NA)
write_hist(os.path.join(d, 'residual_history.csv'), nan_at=150)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'FAIL' and 'step 150' in why(v, 'nonfinite') and v['overall'] == 'FAIL', '残差の履歴に NaN → FAIL (最初の step)', why(v, 'nonfinite'))
d = wall_run('nan_field', spec=QS_NA)
with h5py.File(os.path.join(d, 'mesh.h5'), 'r') as f:
    nn = len(f['MESH/COORD']) // 3
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, nan_field='P')
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'FAIL' and "'P'" in why(v, 'nonfinite'), '最終の場の P に NaN → FAIL', why(v, 'nonfinite'))
d = wall_run('nan_mid', spec=QS_NA)
write_res(os.path.join(d, 'res_1000.h5'), nn, nan_field='T')
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'PASS', '途中の場の NaN は通常の監査 (残差の全履歴 + 最終の場) では読まない (AGENTS.md の範囲)')
rc, v, _ = judge(d, '--solver-rc', '0', '--full-audit')
check(st(v, 'nonfinite') == 'FAIL' and 'res_1000.h5' in why(v, 'nonfinite'), '--full-audit で途中の場の NaN も FAIL', why(v, 'nonfinite'))
d = wall_run('resnan', spec=QS_NA)
write_res(os.path.join(d, 'res_nan_1500.h5'), nn)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'FAIL' and 'res_nan' in why(v, 'nonfinite'), 'forge の NaN のダンプ res_nan_* → FAIL', why(v, 'nonfinite'))
d = wall_run('neg_ro', spec=QS_NA)
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, ro=-1.0)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'physical') == 'FAIL', '最終の場の ro ≤ 0 → 物理性 FAIL', why(v, 'physical'))
print('== 静圧 ≤ 全圧 (total_quantities の h0 から; plan §5.1 #12a (6))')
d = wall_run('pt_ok', spec=QS_NA)
rc, v, _ = judge(d, '--solver-rc', '0')
ph = v['checks']['physical']
check(st(v, 'physical') == 'PASS' and ph['static_le_total']['n_violations'] == 0 and ph['static_le_total']['method'] == 'CPG'
      and 'P ≤ P0' in why(v, 'physical'), 'h0 = c_p T + u²/2 → P ≤ P0 で PASS (CPG)', why(v, 'physical'))
d = wall_run('pt_bad', spec=QS_NA)
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, h0=CP * 290.0)          # 全温 290 K < 静温 300 K → 全圧 < 静圧
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'physical') == 'FAIL' and '静圧 > 全圧' in why(v, 'physical') and v['overall'] == 'FAIL' and rc == 1,
      'h0 < c_p T (全温 < 静温) → 静圧 > 全圧で物理性 FAIL・総合 FAIL', why(v, 'physical'))
tiny = CP * 300.0 * (1.0 - 1.6e-7)                                            # float32 の丸め相当 (P/P0 − 1 ≈ 5.6e-7; 実測の最大 6.4e-7・TP 1.0e-6)
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, h0=tiny)
rc, v, _ = judge(d, '--solver-rc', '0')
mx = v['checks']['physical']['static_le_total']['max_P_over_P0_minus_1']
check(st(v, 'physical') == 'PASS' and 0.0 < mx < JR.PT_REL_TOL, '丸め相当の P/P0 − 1 = %.1e は許容 %.0e の内側 → PASS' % (mx, JR.PT_REL_TOL),
      why(v, 'physical'))
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, h0=None)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'physical') == 'INDETERMINATE' and '静圧 ≤ 全圧を確かめられない' in why(v, 'physical') and 'h0' in why(v, 'physical'),
      'h0 が無い → 静圧 ≤ 全圧は判定不能 (ρ・P・T > 0 は PASS と書く)', why(v, 'physical'))
write_res(os.path.join(d, 'res_%d.h5' % S_LAST), nn, h0=None, ro=-1.0)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'physical') == 'FAIL', 'ρ ≤ 0 と静圧 ≤ 全圧の判定不能の併存 → FAIL (既知の FAIL を隠さない)', why(v, 'physical'))
d = wall_run('pt_nocp', spec=QS_NA)
_c = os.path.join(d, 'solverConfig.yaml')
_t = open(_c).read()
open(_c, 'w').write(_t.replace('cp: %r, ' % CP, ''))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'physical') == 'INDETERMINATE' and 'cp' in why(v, 'physical'), '物性 (cp) が読めない → 静圧 ≤ 全圧は判定不能', why(v, 'physical'))

print('== 絶対の判定 (0 を横切る量: 質量の収支)')
CSV_R = {'kind': 'csv_column', 'file': 'force_history.csv', 'column': 'r'}
ABS = {'criterion': 'absolute', 'thresholds': {'abs_mean': 5e-4, 'abs_amp': 5e-3}}
d = wall_run('abs', spec=qspec(qitem('r', CSV_R, **ABS)))
write_series(d, [(s_, 1e-4 * (-1) ** (s_ // 100)) for s_ in range(100, 2001, 100)], cols=('step', 'r'))
rc, v, _ = judge(d, '--solver-rc', '0')
q = v['checks']['quasisteady']['quantities']['r']
check(st(v, 'quasisteady') == 'PASS' and q['verdict'] == 'STEADY' and q['criterion'] == 'absolute',
      '|平均| ≤ 5e-4・振幅 ≤ 5e-3 → STEADY (0 を横切っても相対の drift を使わない)', why(v, 'quasisteady'))
write_series(d, [(s_, 2e-3 + 1e-5 * (-1) ** (s_ // 100)) for s_ in range(100, 2001, 100)], cols=('step', 'r'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'FAIL' and v['checks']['quasisteady']['quantities']['r']['verdict'] == 'OFFSET',
      '|平均| 2e-3 > 5e-4 → OFFSET (受け入れられない) = FAIL', why(v, 'quasisteady'))
write_series(d, [(s_, 1e-2 * (-1) ** (s_ // 100)) for s_ in range(100, 2001, 100)], cols=('step', 'r'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'quasisteady') == 'FAIL' and v['checks']['quasisteady']['quantities']['r']['verdict'] == 'OSCILLATING',
      '振幅 1e-2 > 5e-3 → OSCILLATING (accept に無ければ FAIL)', why(v, 'quasisteady'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(qspec(qitem('r', CSV_R, criterion='absolute', thresholds={'drift': 0.1})), f)
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['spec']['error'] and 'thresholds' in v['spec']['error'], 'absolute に drift・osc を書く (abs_mean・abs_amp が無い) → 仕様の誤り',
      v['spec'].get('error'))
with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as f:
    json.dump(dict(qspec(qitem('r', CSV_R, **ABS)), notes=['試験の注記: 生成器が書いた']), f)
write_series(d, [(s_, 1e-4) for s_ in range(100, 2001, 100)], cols=('step', 'r'))
rc, v, _ = judge(d, '--solver-rc', '0')
check('仕様の注記: 試験の注記' in open(os.path.join(d, 'RUN_VERDICT.txt'), encoding='utf-8').read(), '仕様の notes を RUN_VERDICT.txt に写す')

d = wall_run('nores', spec=QS_NA, res=())
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'INDETERMINATE' and st(v, 'physical') == 'INDETERMINATE' and st(v, 'wall_resolution') == 'INDETERMINATE',
      '最終の場が無い → 非有限・物理性・壁解像は判定不能')
d = wall_run('nohist', spec=QS_NA, hist=False)
rc, v, _ = judge(d, '--solver-rc', '0')
check(st(v, 'nonfinite') == 'INDETERMINATE' and st(v, 'convergence') == 'INDETERMINATE', '残差の履歴が無い → 判定不能')

# ------------------------------------------------------------------------------------------------ 鮮度・例外・段の記録

print('== 古い判定 (--check-fresh)')
d = wall_run('fresh', spec=QS_NA)
judge(d, '--solver-rc', '0')
rc, _, out = judge(d, '--check-fresh')
check(rc == 0 and out.startswith('FRESH'), '判定の直後は FRESH (終了 0)', out)
write_hist(os.path.join(d, 'residual_history.csv'), d1=-6.0)
rc, _, out = judge(d, '--check-fresh')
check(rc == 1 and 'residual_history.csv' in out, '残差の履歴が変わった → STALE (終了 1)', out)
judge(d, '--solver-rc', '0')
with open(os.path.join(d, 'judge.yaml'), 'w') as f:
    f.write('quasisteady: {not_applicable: x}\n')
rc, _, out = judge(d, '--check-fresh')
check(rc == 1 and 'judge.yaml' in out, '仕様が増えた → STALE', out)
judge(d, '--solver-rc', '0')
write_res(os.path.join(d, 'res_3000.h5'), nn)
rc, _, out = judge(d, '--check-fresh')
check(rc == 1 and 'res_*' in out, '出力が増えた (最後の step が変わる) → STALE', out)
rc, _, out = judge(os.path.join(TMP, '01_pass'), '--check-fresh')
check(rc == 0, '他の run の判定は FRESH のまま', out)
dd = tempfile.mkdtemp(dir=TMP)
rc, _, out = judge(dd, '--check-fresh')
check(rc == 2, '判定が無い → 終了 2', out)

print('== 判定器の例外 (終了 3) と引数の誤り')
d = wall_run('boom', spec=QS_NA)
orig = JR.check_physical
JR.check_physical = lambda run: (_ for _ in ()).throw(RuntimeError('試験の例外'))
buf = io.StringIO()
with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
    rc = JR.main([d, '--solver-rc', '0'])
JR.check_physical = orig
v = json.load(open(os.path.join(d, 'RUN_VERDICT.json')))
check(rc == 3 and v['state'] == 'error' and v['overall'] == 'INDETERMINATE' and '試験の例外' in v['overall_reason'],
      '判定器の例外 → 終了 3・state error・INDETERMINATE を書く', f'rc {rc}')
r = subprocess.run([sys.executable, JUDGE, d, '--no-such-flag'], capture_output=True, text=True)
check(r.returncode == 3, '引数の誤り → 終了 3 (2 = INDETERMINATE と取り違えない)', r.returncode)

print('== 段の記録 (--defer・FORGE_JUDGE_DEFER)')
d = staged('defer')
rc, v, _ = judge(d, '--from-run-case', '--solver-rc', '0')
check(rc == 4 and v['state'] == 'deferred' and v['overall'] == 'INDETERMINATE' and len(v['stages']) == 1 and not v['checks'],
      'manifest のある run を run_case.sh から → 段の記録だけ (state deferred・終了 4・総合は確定しない)', f'rc {rc}')
rc, v, _ = judge(d, '--from-run-case', '--solver-rc', '1')
check(rc == 4 and len(v['stages']) == 2 and v['overall'] == 'FAIL', '2 段目の記録を足す・段の既知の FAIL は総合に出す (確定はしない)',
      v.get('overall_reason'))
rc, v, _ = judge(d, '--solver-rc', '0')
check(v['state'] == 'final' and len(v.get('stages') or []) == 2 and st(v, 'convergence') == 'PASS',
      '引数なし (親の確定) → state final・段の記録は参考として引き継ぐ', v['state'])
d = wall_run('defer_env', spec=QS_NA)
rc, v, _ = judge(d, '--from-run-case', '--solver-rc', '0', env={'FORGE_JUDGE_DEFER': '1'})
check(rc == 4 and v['state'] == 'deferred' and 'FORGE_JUDGE_DEFER' in v['state_note'], 'FORGE_JUDGE_DEFER=1 → manifest が無くても段の記録だけ')
rc, v, _ = judge(d, '--solver-rc', '0', env={'FORGE_JUDGE_DEFER': '1'})
check(v['state'] == 'final', 'FORGE_JUDGE_DEFER は run_case.sh からの呼び出しだけに効く (手で回せば確定する)')

# ------------------------------------------------------------------------------------------------ run_case.sh の通し

print('== run_case.sh の通し (偽の forge)')
FAKE = os.path.join(TMP, 'fake_forge.py')
with open(FAKE, 'w') as f:
    f.write('''#!%s
# 偽の forge (試験用): 残差の履歴・体積の出力を書き、正常終了の行を出す。FAKE_FORGE_RC で終了コード、FAKE_FORGE_STALL で停滞
import os, sys
import h5py, numpy as np
open('FORGE_STARTED', 'w').write('started')
stall = os.environ.get('FAKE_FORGE_STALL') == '1'
with open('residual_history.csv', 'w') as f:
    f.write('step,inner,phase,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\\n')
    for i in range(200):
        v = 1e-2 if stall else 10.0 ** (-5.0 * i / 199.0)
        for ph in ('outer_begin', 'outer_end'):
            f.write('%%d,-1,%%s,%%s\\n' %% (i, ph, ','.join(['%%.6e' %% v] * 5)))
n = len(h5py.File('mesh.h5', 'r')['MESH/COORD']) // 3
for s in (100, 200):
    with h5py.File('res_%%d.h5' %% s, 'w') as h:
        for k, val in (('ro', 1.2), ('P', 1e5), ('T', 300.0), ('roUx', 1.0), ('h0', 1004.5 * 300.0)):
            h.create_dataset('VALUE/' + k, data=np.full(n, val))
print('Time = 0.010 s (wall, 200 steps, 0.05 ms/step)')
sys.exit(int(os.environ.get('FAKE_FORGE_RC', '0')))
''' % sys.executable)
os.chmod(FAKE, os.stat(FAKE).st_mode | stat.S_IXUSR)


def rc_run(tag, bcond=True, viscous=False, wall_out=0):
    """関門を通る node の小さな run (test_gate_bad_input の lg_run と同じ置き方、上辺が壁)。既定は非粘性 (壁解像は適用外)。
    wall_out: 壁の outputHDFflg (粘性の run で 0 なら関門が止める)。"""
    _n[0] += 1
    d = os.path.join(TMP, '%02d_rc_%s' % (_n[0], tag))
    os.makedirs(d)
    ni, nj = 11, 4
    X, Y = np.meshgrid(np.linspace(0.0, 1.0, ni), np.linspace(0.0, 0.2, nj), indexing='ij')
    nid = lambda i, j: i * nj + j                                                     # noqa: E731
    quads = [[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for i in range(ni - 1) for j in range(nj - 1)]
    wall = [(nid(i, nj - 1), nid(i + 1, nj - 1)) for i in range(ni - 1)]
    with h5py.File(os.path.join(d, 'mesh.h5'), 'w') as f:
        f.create_dataset('MESH/COORD', data=np.c_[X.ravel(), Y.ravel(), np.zeros(ni * nj)].ravel())
        f.create_dataset('MESH/CONNE', data=np.array([x for n in range(ni * nj) for x in (1, n)], np.int64))
        f.create_dataset('VIZMESH/CONNE', data=np.array([x for q in quads for x in [5] + q], np.int64))
        f['VIZMESH'].attrs['nVizCells'] = len(quads)
        f.create_dataset('VALUE/wall_dist', data=0.2 - Y.ravel())
        f.create_dataset('BCONDS/3/vizBfaceNodes', data=np.array(wall, np.int64).ravel())
        f.create_dataset('BCONDS/3/vizBfaceSizes', data=np.full(len(wall), 2, np.int64))
        f['BCONDS/3'].attrs['bcondKind'] = 'wall'
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as fh:
        fh.write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5", discretization: "node"}\n'
                 'physProp: {thermalMethod: 0, cp: 1004.5, gamma: 1.4, viscMethod: 0, visc: %s, thermCond: 0.0}\n'
                 % ('1.8e-5' if viscous else '0.0'))
    if bcond:
        with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as fh:
            fh.write('wall: {physID: 3, kind: wall, outputHDFflg: %d, ints: , floats: }\n' % wall_out)
    with open(os.path.join(d, 'JUDGE_SPEC.json'), 'w') as fh:
        json.dump(QS_NA, fh)
    return d


def run_case(d, **env_extra):
    env = dict(os.environ, FORGE_BIN=FAKE)
    for k in ('FORGE_LAUNCH_GATE', 'FORGE_STRICT_CONVERGENCE', 'FORGE_JUDGE_DEFER', 'FAKE_FORGE_RC', 'FAKE_FORGE_STALL'):
        env.pop(k, None)
    env.update(env_extra)
    p = subprocess.run(['bash', RUN_CASE, d], capture_output=True, text=True, env=env)
    v = json.load(open(os.path.join(d, 'RUN_VERDICT.json'))) if os.path.exists(os.path.join(d, 'RUN_VERDICT.json')) else None
    return p.returncode, p.stdout + p.stderr, v


d = rc_run('ok')
rc, out, v = run_case(d)
check(rc == 0 and os.path.exists(os.path.join(d, 'FORGE_STARTED')) and v and v['state'] == 'final' and v['overall'] == 'PASS',
      '関門を通る run → forge (偽物) を起動し、RUN_VERDICT.* が総合 PASS で残る・終了 0', f'rc {rc} {out[-1500:]}')
check(os.path.exists(os.path.join(d, 'RUN_VERDICT.txt')) and '[run_case] judge exit=0' in out, 'RUN_VERDICT.txt・judge の終了コードの行')
check(v and v['execution']['solver_rc'] == 0 and v['checks']['solver']['rc_source'].startswith('run_case.sh'), 'solver の終了コードは run_case.sh が渡す')
check(os.path.exists(os.path.join(d, 'LAUNCH_ECHO.txt')) and os.path.exists(os.path.join(d, 'residual_history.png'))
      and os.path.exists(os.path.join(d, 'CONVERGENCE_VERDICT.txt')) and 'CONVERGENCE_EXIT: 0' in open(os.path.join(d, 'CONVERGENCE_VERDICT.txt')).read(),
      '起動の後のエコー LAUNCH_ECHO.txt・残差のプロット・CONVERGENCE_VERDICT.txt (副判定) は今のまま残る')
check(v and st(v, 'wall_resolution') == 'NOT_APPLICABLE', '非粘性の run は壁解像が適用外')

d = rc_run('visc_nowall_out', viscous=True)
rc, out, v = run_case(d)
gate_txt = open(os.path.join(d, 'LAUNCH_GATE.txt'), encoding='utf-8').read() if os.path.exists(os.path.join(d, 'LAUNCH_GATE.txt')) else ''
check(rc != 0 and not os.path.exists(os.path.join(d, 'FORGE_STARTED')) and v is None and '[壁の出力]' in gate_txt
      and 'judge exit' not in out,
      '定数粘性の壁で壁の出力が無い (outputHDFflg 0) → 投入の関門が止める (計算の後では壁解像を判定できない)・判定しない', f'rc {rc}')
d = rc_run('visc_wall_out', viscous=True, wall_out=1)
rc, out, v = run_case(d)
check(rc == 0 and v and st(v, 'wall_resolution') == 'INDETERMINATE' and 'res_wall_3_200.h5' in why(v, 'wall_resolution'),
      '同じ run で壁の outputHDFflg 1 → 起動する (偽の forge は壁の出力を書かないので壁解像は判定不能)', f'rc {rc} {why(v, "wall_resolution")}')

d = rc_run('nobc', bcond=False)
rc, out, v = run_case(d)
check(rc != 0 and not os.path.exists(os.path.join(d, 'FORGE_STARTED')) and v is None and not os.path.exists(os.path.join(d, 'RUN_VERDICT.txt'))
      and 'judge exit' not in out, '関門で未起動 (bcondConfig.yaml が無い) → 判定しない (RUN_VERDICT を書かない)', f'rc {rc}')

d = rc_run('forge_fail')
rc, out, v = run_case(d, FAKE_FORGE_RC='3')
check(rc == 3 and v and st(v, 'solver') == 'FAIL' and v['overall'] == 'FAIL', 'forge が 3 で終わる → run_case は 3 のまま・判定は solver FAIL', f'rc {rc}')

d = rc_run('stall')
rc, out, v = run_case(d, FAKE_FORGE_STALL='1')
check(rc == 0 and v and st(v, 'convergence') == 'FAIL' and v['overall'] == 'FAIL',
      '停滞の run → 判定は FAIL だが run_case の終了コードは 0 (FORGE_STRICT_CONVERGENCE なし)', f'rc {rc}')
d = rc_run('stall_strict')
rc, out, v = run_case(d, FAKE_FORGE_STALL='1', FORGE_STRICT_CONVERGENCE='1')
check(rc == 1 and v and v['overall'] == 'FAIL', 'FORGE_STRICT_CONVERGENCE=1 は従来どおり収束の副判定の終了コード (1) を返す', f'rc {rc}')

d = rc_run('staged')
m = SM.StageManifest(d)
m.add('S1', cfg_txt, 'wall: {}\n', history='residual_history_S1.csv')
m.write()
rc, out, v = run_case(d)
check(rc == 0 and v and v['state'] == 'deferred' and '[run_case] judge exit=4' in out,
      '段階起動の run (stage_manifest.json がある) → run_case は段の記録だけ (state deferred)', f'rc {rc}')
d = rc_run('defer_env')
rc, out, v = run_case(d, FORGE_JUDGE_DEFER='1')
check(rc == 0 and v and v['state'] == 'deferred', '親が FORGE_JUDGE_DEFER=1 を渡す → 段の記録だけ', f'rc {rc}')

# ------------------------------------------------------------------------------------------------ Stop フック

print('== Stop フック (hook_convergence_gate.py: 判定済みか・今の入力に対応するか、同じセッションの範囲)')
HOOK = os.path.join(TOOLS, 'hook_convergence_gate.py')
HROOT = os.path.join(TMP, 'hookroot')
HHOME = os.path.join(TMP, 'hookhome')
os.makedirs(os.path.join(HHOME, '.cache', 'forge-gate-claims'))
with open(os.path.join(HHOME, '.cache', 'forge-gate-claims', 'sess1.txt'), 'w') as f:
    f.write('01.mine\n')
SRC = os.path.join(TMP, '01_pass')                       # 上の正例 (PASS の判定がある run) を写して使う


def hook_run(case, name, verdict=True, age=200.0):
    dst = os.path.join(HROOT, 'case', case, name)
    shutil.copytree(SRC, dst)
    for fn in ('RUN_VERDICT.json', 'RUN_VERDICT.txt'):
        if os.path.exists(os.path.join(dst, fn)):
            os.unlink(os.path.join(dst, fn))
    if verdict:
        judge(dst, '--solver-rc', '0')
    t = time.time() - age
    os.utime(os.path.join(dst, 'residual_history.csv'), (t, t))
    return dst


def hook(session='sess1'):
    env = dict(os.environ, FORGE_GATE_ROOT=HROOT, HOME=HHOME)
    r = subprocess.run([sys.executable, HOOK], input=json.dumps({'session_id': session} if session else {}), capture_output=True,
                       text=True, env=env)
    try:
        return json.loads(r.stdout or '{}'), r
    except ValueError:
        return {'decision': 'parse-error'}, r


ok_run = hook_run('01.mine', 'run_0001_ok')
out, r = hook()
check(out == {} and r.returncode == 0, '判定済み・今の入力に対応 (FRESH) → block しない', r.stdout + r.stderr)
hook_run('01.mine', 'run_0002_fail')
_v = json.load(open(os.path.join(HROOT, 'case', '01.mine', 'run_0002_fail', 'RUN_VERDICT.json')))
write_hist(os.path.join(HROOT, 'case', '01.mine', 'run_0002_fail', 'residual_history.csv'), d0=-2.0, d1=-2.0)
judge(os.path.join(HROOT, 'case', '01.mine', 'run_0002_fail'), '--solver-rc', '0')
t_ = time.time() - 200
os.utime(os.path.join(HROOT, 'case', '01.mine', 'run_0002_fail', 'residual_history.csv'), (t_, t_))
_v = json.load(open(os.path.join(HROOT, 'case', '01.mine', 'run_0002_fail', 'RUN_VERDICT.json')))
out, r = hook()
check(out == {} and _v['overall'] == 'FAIL', 'FAIL を正しく判定した run → block しない (PASS を要求しない)', r.stdout)
nov = hook_run('01.mine', 'run_0003_nov', verdict=False)
out, r = hook()
check(out.get('decision') == 'block' and 'run_0003_nov' in out.get('reason', '') and 'RUN_VERDICT.json が無い' in out['reason']
      and 'run_0001_ok' not in out['reason'] and 'run_0002_fail' not in out['reason'], '判定の無い run → block (その run だけを名指し)',
      out.get('reason'))
judge(nov, '--solver-rc', '0')
os.utime(os.path.join(nov, 'residual_history.csv'), (time.time() - 200,) * 2)
out, r = hook()
check(out == {}, 'judge_run.py を回した後は block しない', out.get('reason'))
write_hist(os.path.join(nov, 'residual_history.csv'), d1=-6.0)
os.utime(os.path.join(nov, 'residual_history.csv'), (time.time() - 200,) * 2)
out, r = hook()
check(out.get('decision') == 'block' and 'STALE' in out.get('reason', '') and 'residual_history.csv' in out['reason'],
      '判定の後に入力 (残差の履歴) が変わった (古い判定) → block', out.get('reason'))
judge(nov, '--solver-rc', '0')
os.utime(os.path.join(nov, 'residual_history.csv'), (time.time() - 200,) * 2)
dfr = hook_run('01.mine', 'run_0004_deferred', verdict=False)
judge(dfr, '--defer', '--solver-rc', '0')
out, r = hook()
check(out.get('decision') == 'block' and 'run_0004_deferred' in out.get('reason', '') and 'deferred' in out['reason'],
      '段の記録 (state deferred) のまま親が確定せずに終わった → block', out.get('reason'))
with open(os.path.join(dfr, 'JUDGE_PENDING.json'), 'w') as f:
    json.dump({'pid': os.getpid(), 'host': socket.gethostname()}, f)
out, r = hook()
check(out == {}, '親 (runner) が生きている (JUDGE_PENDING.json の pid) → 確定待ちとして block しない', out.get('reason'))
with open(os.path.join(dfr, 'JUDGE_PENDING.json'), 'w') as f:
    json.dump({'pid': 2 ** 22 + 12345, 'host': socket.gethostname()}, f)
out, r = hook()
check(out.get('decision') == 'block' and 'run_0004_deferred' in out.get('reason', ''), '親が死んでいる (pid が無い) → block', out.get('reason'))
judge(dfr, '--solver-rc', '0')
os.utime(os.path.join(dfr, 'residual_history.csv'), (time.time() - 200,) * 2)
err = hook_run('01.mine', 'run_0005_error', verdict=False)
orig = JR.check_physical
JR.check_physical = lambda run: (_ for _ in ()).throw(RuntimeError('試験の例外'))
with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    JR.main([err, '--solver-rc', '0'])
JR.check_physical = orig
os.utime(os.path.join(err, 'residual_history.csv'), (time.time() - 200,) * 2)
out, r = hook()
check(out == {}, '判定器の例外 (state error = 判定不能) を報告して終えるのは block しない', out.get('reason'))
hook_run('01.mine', 'run_0006_active', verdict=False, age=10.0)
out, r = hook()
check(out == {}, '残差の履歴を書き込み中 (更新 120 s 未満) は見ない', out.get('reason'))
hook_run('02.other', 'run_0001_theirs', verdict=False)
out, r = hook()
check(out == {}, '他のセッションの case (claim していない) の判定の無い run → このセッションは block しない', out.get('reason'))
out, r = hook(session=None)
check(out.get('decision') == 'block' and '02.other' in out.get('reason', ''),
      'session_id が無ければ絞り込まず全件を見る (他の case の判定の無い run も block)', out.get('reason'))

print()
if fails:
    print('%d FAIL (出力を残した: %s)' % (fails, TMP))
    sys.exit(1)
shutil.rmtree(TMP, ignore_errors=True)
print('ALL PASS')
