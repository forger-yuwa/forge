#!/usr/bin/env python3
"""recipe_settings.py の試験 (plan tooling-generator-defaults-and-launch-gate §5.1 #6・§6.2 の設定の負例・正例)。

負例はそれぞれ対象のキーを名指しで FAIL: 宣言したプロファイルと違う nStepInner (逸脱の宣言なし)・bndFirstOrder・dual-time の
implicitRelax・未知のキー・誤った節のキー・YAML の重複キー・merge key・理由の無い逸脱・段の不一致 (soft の宣言で 2 次の config)・
プロファイルに無いキー。advisory の範囲外は WARN で終了 0。正例: プロファイルの値で組んだ config・#4 の実在の正当な config
(「実験」の宣言)。継承の衝突 (違う値・違う種類・2 つの親・循環) は読み込み時に例外。「実験」の宣言は理由が無いと書けない。
文書の同期: procedures/recommended-settings.md の生成区間 = `recipe_settings.py --doc` (plan §5.1 #5)。
生成器の出力そのものの正例と書き写しの台本の一致は design/tests/run_recipe_profiles_tests.py (design の venv が要る)。
"""
import json
import os
import shutil
import sys
import tempfile

here = os.path.dirname(os.path.abspath(__file__))
tools = os.path.abspath(os.path.join(here, '..', '..', 'tools'))
root = os.path.abspath(os.path.join(here, '..', '..', '..'))
sys.path.insert(0, tools)
import yaml  # noqa: E402

import recipe_settings as RS  # noqa: E402

fails = 0
TMP = tempfile.mkdtemp(prefix='recipe_settings_')


def check(cond, msg):
    global fails
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond:
        fails += 1


def fail_keys(r):
    return {k for k, _ in r['fails']}


def warn_keys(r):
    return {k for k, _ in r['warns']}


def unflatten(flat):
    out = {}
    for k, v in flat.items():
        cur = out
        parts = k.split('.')
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    return out


RC = RS.load_recipes()
PROBLEM = {'mesh.meshFileName': 'nozzle.h5', 'mesh.valueFileName': 'nozzle.h5', 'gpu': 1, 'initial': 'uniform_p101325_u10',
           'physProp.cp': 1004.5, 'physProp.gamma': 1.4, 'physProp.visc': 1.8e-5, 'physProp.thermCond': 0.0257,
           'time.last.nStepOuter': 12000, 'time.outStepInterval': 2000, 'time.outStepStart': 0,
           'time.deltaT.dt': 1e-8, 'time.deltaT.dt_min': 1e-9, 'time.deltaT.dt_max': 1e-3}


def cfg_of(profile, stage, edit=None, extra=None):
    """プロファイルの段の fixed (と問題のキー) で組んだ config の文字列。"""
    res = RS.resolve(RC, profile)
    flat = dict(PROBLEM)
    flat.update(res['stages'][stage]['fixed'])
    if 'physProp.viscMethod' in flat and int(flat['physProp.viscMethod']) == 0:
        flat['physProp.visc'] = flat['physProp.thermCond'] = 0.0
    if extra:
        flat.update(extra)
    doc = unflatten(flat)
    if edit:
        edit(doc)
    return yaml.safe_dump(doc, sort_keys=False)


def run_dir(tag, text, bcond=None):
    d = os.path.join(TMP, tag)
    os.makedirs(d)
    with open(os.path.join(d, 'solverConfig.yaml'), 'w') as f:
        f.write(text)
    if bcond is not None:
        with open(os.path.join(d, 'bcondConfig.yaml'), 'w') as f:
            f.write(bcond)
    return d


BC_NS = ('inlet: {physID: 1, kind: inlet_Pressure, outputHDFflg: 0, ints: , floats: {Pt: 1.0e6, Tt: 3000.0, k: 1.0, omega: 18000.0}}\n'
         'wall: {physID: 3, kind: wall, outputHDFflg: 1, ints: , floats: }\n'
         'axis: {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n')

# --- 読み込み・一覧 --------------------------------------------------------------------------------------------------
names = RS.declarable(RC)
check('common_node_ns_sst' in names and 'common_node_euler' in names, f'共通プロファイルがある ({names[:3]}…)')
check(all(n in names for n in ('bell_tp_ns_sst', 'axismach_euler_tp', 'axismach_ns_sst_tp', 'axismach_ns_sst_tp_b0',
                                'sern2d_finite_base_frozen_tp_sst', 'sern3d_frozen_tp_sst', 'walldriven_euler_cpg',
                                'wt_modef_euler_cell')), '機種ごとのプロファイル (書き写し) がある')
for n in names:
    r = RS.resolve(RC, n)
    check('main' in r['stages'] and r['source'], f'{n}: 本段と根拠 (source) がある')
check(all(d['acceptance'] == '未' for n in names for d in RS.resolve(RC, n)['differs_from_recommended']),
      '§1 との差の受入の状態は全部「未」')
check('time.nStepInner' in {d['key'] for d in RS.resolve(RC, 'axismach_ns_sst_tp')['differs_from_recommended']},
      'axismach_ns_sst_tp の §1 との差に nStepInner (5 vs 4) がある')

# --- 正例: プロファイルの値で組んだ config --------------------------------------------------------------------------------
for prof in ('common_node_ns_sst', 'common_node_euler', 'bell_tp_ns_sst', 'axismach_ns_sst_tp', 'sern2d_finite_base_frozen_tp_sst'):
    for st in RS.resolve(RC, prof)['stages']:
        r = RS.check_config(RC, cfg_of(prof, st), prof, st)
        check(not r['fails'] and not r['warns'], f'正例 {prof}/{st}: PASS (fails {sorted(fail_keys(r))} warns {sorted(warn_keys(r))})')

# --- 負例: nStepInner (逸脱の宣言なし / あり / 値の食い違い / 理由なし) -------------------------------------------------
txt = cfg_of('bell_tp_ns_sst', 'main', extra={'time.nStepInner': 5})
r = RS.check_config(RC, txt, 'bell_tp_ns_sst', 'main')
check(fail_keys(r) == {'time.nStepInner'}, f'nStepInner 5 (プロファイル 4) は名指しで FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, txt, 'bell_tp_ns_sst', 'main', deviations={'time.nStepInner': 5})
check(not r['fails'] and r['applied'] == [{'key': 'time.nStepInner', 'value': 5, 'profile_value': 4}],
      f'逸脱を宣言すれば通り、逸脱が記録される ({r["applied"]})')
r = RS.check_config(RC, txt, 'bell_tp_ns_sst', 'main', deviations={'time.nStepInner': 6})
check('time.nStepInner' in fail_keys(r), '逸脱の宣言 (6) と config (5) が違えば FAIL')
r = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main'), 'bell_tp_ns_sst', 'main', deviations={'time.nStepInner': 4})
check(not r['fails'] and 'time.nStepInner' in warn_keys(r), '使われていない逸脱の宣言は WARN (config はプロファイルどおり)')
d = run_dir('no_reason', txt)
rc, r, path = RS.declare(d, 'bell_tp_ns_sst', 'main', deviate=['time.nStepInner=5'], reason=None, recipes=RC)
check(rc == 2 and path is None and not os.path.exists(os.path.join(d, 'RECIPE.json')) and '(deviate)' in fail_keys(r),
      f'理由の無い逸脱は宣言できない (rc {rc})')
rc, r, path = RS.declare(d, 'bell_tp_ns_sst', 'main', deviate=['time.nStepInner=5'], reason='試し: 内部反復 5 の比較', recipes=RC)
rec = json.load(open(path)) if path else {}
import hashlib
check(rc == 0 and rec.get('kind') == 'profile' and rec.get('profile') == 'bell_tp_ns_sst' and rec.get('stage') == 'main'
      and rec['deviations'][0]['key'] == 'time.nStepInner' and rec['deviations'][0]['reason'].startswith('試し')
      and rec['deviations'][0]['basis'] is None
      and rec['config']['sha256'] == hashlib.sha256(txt.encode()).hexdigest()
      and rec['profile_sha256'] == RS.resolve(RC, 'bell_tp_ns_sst')['sha256'] and rec['profile_version'] == RC['version'],
      'RECIPE.json: プロファイル・版・段・逸脱 (理由、根拠は任意)・config の sha256')

# --- 負例: 禁止・dual-time の緩和・未知のキー・誤った節のキー ----------------------------------------------------------------
r = RS.check_config(RC, cfg_of('common_node_ns_sst', 'main', extra={'mesh.bndFirstOrder': 1}), 'common_node_ns_sst', 'main',
                    deviations={'mesh.bndFirstOrder': 1})
check('mesh.bndFirstOrder' in fail_keys(r) and any('逸脱として宣言できない' in m for k, m in r['fails']),
      'bndFirstOrder は FAIL、逸脱にもできない')
dual = cfg_of('common_node_ns_sst', 'main', extra={'time.unsteady': 1, 'time.dualTime': 1, 'time.nSubIterDualTime': 20,
                                                   'time.deltaT.cfl_pseudo': 12.0})
r = RS.check_config(RC, dual, 'common_node_ns_sst', 'main')
check('time.deltaT.implicitRelax' in fail_keys(r), f'dual-time の implicitRelax 0.7 は名指しで FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, cfg_of('common_node_ns_sst', 'main', extra={'time.deltaT.notAKey': 1}), 'common_node_ns_sst', 'main')
check('time.deltaT.notAKey' in fail_keys(r), '未知のキーは名指しで FAIL')
r = RS.check_config(RC, cfg_of('common_node_ns_sst', 'main', extra={'time.deltaT.nStepInner': 4}), 'common_node_ns_sst', 'main')
check('time.deltaT.nStepInner' in fail_keys(r) and any('time.nStepInner' in m for k, m in r['fails'] if k == 'time.deltaT.nStepInner'),
      '誤った節のキー (time.deltaT.nStepInner) は名指しで FAIL し、正しい節を示す')
r = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', extra={'space.slauWallNormalChi': 0}), 'bell_tp_ns_sst', 'main')
check(fail_keys(r) == {'space.slauWallNormalChi'}, f'プロファイルに無いキーは名指しで FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', extra={'space.slauWallNormalChi': 0}), 'bell_tp_ns_sst', 'main',
                    deviations={'space.slauWallNormalChi': 0})
check(not r['fails'], 'プロファイルに無いキーも逸脱を宣言すれば通る')
r2 = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', edit=lambda d: d['time']['deltaT'].pop('implicitRelax')),
                     'bell_tp_ns_sst', 'main')
check('time.deltaT.implicitRelax' in fail_keys(r2) and any('書いていない' in m for k, m in r2['fails']),
      'プロファイルの fixed を書いていない (solver 既定 1.0 ≠ 0.7) は FAIL')
r3 = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', edit=lambda d: d['time']['deltaT'].pop('implicitRelax')),
                     'bell_tp_ns_sst', 'main', deviations={'time.deltaT.implicitRelax': None})
check(not r3['fails'], '「書かない」(~) の逸脱を宣言すれば通る')

# --- 負例: YAML の重複キー・merge key ---------------------------------------------------------------------------------
base = cfg_of('bell_tp_ns_sst', 'main')
dup = base.replace('  nStepInner: 4\n', '  nStepInner: 4\n  nStepInner: 5\n', 1)
check(dup != base, '(重複キーの入力を作れた)')
r = RS.check_config(RC, dup, 'bell_tp_ns_sst', 'main')
check('time.nStepInner' in fail_keys(r), f'YAML の重複キーは完全修飾パスで FAIL ({sorted(fail_keys(r))})')
merge = 'common: &c {convMethod: 1, limiter: 2}\n' + base.replace('space:\n', 'space:\n  <<: *c\n', 1)
r = RS.check_config(RC, merge, 'bell_tp_ns_sst', 'main')
check(any('merge' in m for k, m in r['fails']), f'merge key は FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, 'time: [1, 2\n', 'bell_tp_ns_sst', 'main')
check('(config)' in fail_keys(r), '読めない YAML は検査できない入力として FAIL')

# --- 負例: 段の不一致・無い段・抽象・宣言の欠落 ----------------------------------------------------------------------------
r = RS.check_config(RC, cfg_of('axismach_ns_sst_tp', 'main'), 'axismach_ns_sst_tp', 'soft')
check('space.convMethod' in fail_keys(r), f'soft の宣言で 2 次の config は space.convMethod を名指しで FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, base, 'bell_tp_ns_sst', 'warm_lam')
check('(stage)' in fail_keys(r), 'プロファイルに無い段の宣言は FAIL')
r = RS.check_config(RC, base, 'common_node_steady', 'main')
check('(profile)' in fail_keys(r), '抽象のプロファイルは宣言できない')
r = RS.check_config(RC, base, None, None)
check('(profile)' in fail_keys(r), 'プロファイルの宣言が無ければ FAIL')

# --- advisory の範囲外は WARN で終了 0 --------------------------------------------------------------------------------------
txt = cfg_of('common_node_ns_sst', 'main', extra={'time.deltaT.cfl_pseudo': 8.0, 'time.deltaT.cfl': 8.0})
d = run_dir('advisory', txt, BC_NS)
rc, r, path = RS.declare(d, 'common_node_ns_sst', 'main', recipes=RC)
check(rc == 0 and not r['fails'] and warn_keys(r) == {'time.deltaT.cfl_pseudo'} and path is not None,
      f'cfl_pseudo 8 (目安 4〜6) は WARN で終了 0、宣言を書く ({sorted(warn_keys(r))})')
check(RS.main(['check', d, '--profile', 'common_node_ns_sst', '--stage', 'main']) == 0, 'CLI check: WARN は終了 0')
check(RS.main(['check', d, '--profile', 'common_node_ns_sst', '--stage', 'soft']) == 1, 'CLI check: 段の不一致は終了 1')
d2 = run_dir('cli_fail', cfg_of('bell_tp_ns_sst', 'main', extra={'time.nStepInner': 5}))
check(RS.main(['declare', d2, '--profile', 'bell_tp_ns_sst', '--stage', 'main']) == 1
      and not os.path.exists(os.path.join(d2, 'RECIPE.json')), 'CLI declare: FAIL なら書かずに終了 1')

# --- 「実験」の宣言 ----------------------------------------------------------------------------------------------------
txt = cfg_of('bell_tp_ns_sst', 'main', extra={'time.nStepInner': 5, 'space.slauWallNormalChi': 0})
d = run_dir('exp', txt)
rc, r, path = RS.declare(d, 'bell_tp_ns_sst', 'main', experiment=True, reason='', recipes=RC)
check(rc == 2 and path is None and not os.path.exists(os.path.join(d, 'RECIPE.json')), '「実験」は理由が無いと書けない')
rc, r, path = RS.declare(d, 'bell_tp_ns_sst', 'main', experiment=True, reason='内部反復 5 と χ 0 の組の試し', recipes=RC)
rec = json.load(open(path)) if path else {}
check(rc == 0 and rec.get('kind') == 'experiment' and rec.get('reason', '').startswith('内部反復')
      and {o['key'] for o in rec.get('off_profile', [])} == {'time.nStepInner', 'space.slauWallNormalChi'},
      f'「実験」は推奨への適合を求めず、外れたキーと理由を書く ({[o["key"] for o in rec.get("off_profile", [])]})')
d = run_dir('exp_bfo', cfg_of('bell_tp_ns_sst', 'main', extra={'mesh.bndFirstOrder': 1}))
rc, r, path = RS.declare(d, 'bell_tp_ns_sst', 'main', experiment=True, reason='x', recipes=RC)
check(rc == 1 and 'mesh.bndFirstOrder' in fail_keys(r), '「実験」でも bndFirstOrder は FAIL')
check(RS.main(['declare', d, '--experiment']) == 2, 'CLI: 「実験」に --reason が無ければ終了 2')

# --- 正例: #4 の実在の正当な config (「実験」の宣言で、未知・誤配置・起動時拒否が無いこと) -------------------------------------------
POSITIVE = ['case/46.sern_design/chi_evidence/46/0968_chidef_2d_m10_on_omit/solverConfig.yaml',
            'case/44.vitiated_air_wt/run_0164_va3_M4.19_Lc8_noneq_inletTt_cell_lim1d_cfl2/solverConfig.yaml',
            'case/46.sern_design/chi_evidence/39/0953_chidef_b1p_omit/solverConfig.yaml',
            'case/09.Taylor-Green/_g0_lsq_seam/solverConfig.yaml',
            'case/09.Taylor-Green/_g0_lsq_seam/s2_evidence/16/0961_sglsq_s2_fct_gg_nokey/solverConfig.yaml',
            'case/09.Taylor-Green/_g0_lsq_seam/s2_evidence/16/0955_sglsq_s2_fct_lsq/solverConfig.yaml']
for rel in POSITIVE:
    src = os.path.join(root, rel)
    if not os.path.exists(src):
        check(False, f'正例が無い: {rel}')
        continue
    d = os.path.join(TMP, 'pos_' + rel.replace('/', '_'))
    os.makedirs(d)
    shutil.copy(src, os.path.join(d, 'solverConfig.yaml'))
    bp = os.path.join(os.path.dirname(src), 'bcondConfig.yaml')
    if os.path.exists(bp):
        shutil.copy(bp, os.path.join(d, 'bcondConfig.yaml'))
    rc, r, path = RS.declare(d, experiment=True, reason='#4 の正例', recipes=RC)
    check(rc == 0 and path is not None, f'正例 (実験の宣言) {rel}: FAIL なし ({r["fails"][:2]})')

# --- 継承の衝突・鎖でない継承・循環 ----------------------------------------------------------------------------------------
BASE = '''version: 1
problem_keys: []
profiles:
  a: {abstract: true, all_stages: {fixed: {time.nStepInner: 4}}, forbid: {mesh.bndFirstOrder: x}}
  b: {inherits: a, stages: {main: {fixed: {space.convMethod: 1}}}}
'''


def raises(text, must):
    try:
        rc = RS.load_recipes(text=text)
        for n in rc['profiles']:
            RS.resolve(rc, n)
    except RS.RecipeError as e:
        return must in str(e), str(e)[:120]
    return False, '例外にならない'


ok, msg = raises(BASE + '  c: {inherits: b, stages: {main: {fixed: {time.nStepInner: 5}}}}\n', '継承の衝突')
check(ok, f'鎖の中で同じキーに違う値は例外 ({msg})')
ok, msg = raises(BASE + '  c: {inherits: b, stages: {main: {advisory: {time.nStepInner: [4, 6]}}}}\n', '継承の衝突')
check(ok, f'鎖の中で同じキーに違う種類 (fixed / advisory) は例外 ({msg})')
ok, msg = raises(BASE + '  c: {inherits: b, stages: {main: {fixed: {mesh.bndFirstOrder: 0}}}}\n', '継承の衝突')
check(ok, f'禁止キーを fixed にするのは例外 ({msg})')
ok, msg = raises(BASE + '  c: {inherits: [a, b], stages: {main: {fixed: {space.limiter: 2}}}}\n', '1 本の鎖')
check(ok, f'2 つの親は例外 ({msg})')
ok, msg = raises(BASE + '  c: {inherits: d, stages: {main: {}}}\n  d: {inherits: c}\n', '循環')
check(ok, f'継承の循環は例外 ({msg})')
ok, msg = raises(BASE + '  c: {inherits: b, stages: {main: {fixed: {time.nStepInner: 4}}}}\n', '')
check(not ok and msg == '例外にならない', '鎖の中で同じキーに同じ値は通る')
ok, msg = raises(BASE.replace('version: 1', 'version: 1\nversion: 2'), '重複キー')
check(ok, f'推奨設定のファイルの重複キーも拒否 ({msg})')

# --- 文書の同期: recommended-settings.md の生成区間 = YAML から作った表 (plan §5.1 #5・§6.2「文書の同期」) -----------------
doc = os.path.join(root, 'procedures', 'recommended-settings.md')
check(RS.doc_block(open(doc, encoding='utf-8').read()) is not None, 'recommended-settings.md に生成区間の印がある')
check(RS.update_doc(doc, RC, check_only=True),
      'recommended-settings.md の生成区間 = recipe_settings.py --doc (違えば `--doc --write procedures/recommended-settings.md`)')

if fails == 0:
    shutil.rmtree(TMP, ignore_errors=True)
print(f"\n{'ALL PASS' if fails == 0 else f'FAIL {fails}'}")
sys.exit(1 if fails else 0)
