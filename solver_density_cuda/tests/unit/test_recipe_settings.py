#!/usr/bin/env python3
"""recipe_settings.py の試験 (plan tooling-generator-defaults-and-launch-gate §5.1 #6・#16・§6.2 の設定の負例・正例)。

負例はそれぞれ対象のキーを名指しで FAIL: 宣言したプロファイルと違う nStepInner (逸脱の宣言なし)・bndFirstOrder・dual-time の
implicitRelax・未知のキー・誤った節のキー・YAML の重複キー・merge key・理由の無い逸脱・段の不一致 (soft の宣言で 2 次の config)・
プロファイルに無いキー。advisory の範囲外は WARN で終了 0。正例: プロファイルの値で組んだ config・#4 の実在の正当な config
(「実験」の宣言)。継承の衝突 (違う値・違う種類・2 つの親・循環) は読み込み時に例外。「実験」の宣言は理由が無いと書けない。
文書の同期: procedures/recommended-settings.md の生成区間 = `recipe_settings.py --doc` (plan §5.1 #5)。
#16 (省略の FAIL と明示すべきキー): `fixed`・`explicit` の省略は solver の既定が同じでも名指しで FAIL (判別 A/B = SERN 2D 本段の
`space.limiterScaled`)、省略を通すのは `key=~` + 理由の逸脱だけ (実効値は未確定として記録)、明示すべきキーの旧い実効値の出所は
今のソースで確かめる (合わなければ例外)。プロファイルの値で組む正例の config は `fixed` と `explicit` の両方を書く
(#16 で `explicit` を足したので、`fixed` だけで組んだ config は機種のプロファイルに FAIL する)。
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


def cfg_of(profile, stage, edit=None, extra=None, explicit=True):
    """プロファイルの段の fixed・explicit (明示すべきキー、#16) と問題のキーで組んだ config の文字列。
    explicit=False は fixed だけ (= 今の生成器と同じく明示すべきキーを省いた config)。"""
    res = RS.resolve(RC, profile)
    flat = dict(PROBLEM)
    flat.update(res['stages'][stage]['fixed'])
    if explicit:
        flat.update(res['stages'][stage]['explicit'])
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

# passiveFct の警告は check_solver_config 本体が「明示した有効化が作動しない」ときだけ出す (plan §5.1 #10)。
# 照合の側で「省略なら落とす」は撤去した: 定常の正例 (省略) は上で WARN 0、明示した 1 は WARN が残る
r = RS.check_config(RC, cfg_of('common_node_ns_sst', 'main', extra={'time.deltaT.passiveFct': 1}), 'common_node_ns_sst', 'main')
check('time.deltaT.passiveFct' in warn_keys(r), f'定常で passiveFct: 1 を明示すると「作動しない」の WARN が残る ({sorted(warn_keys(r))})')
check(r['fixed_fails'] == [], f'プロファイルに無いキーは fixed の不一致に数えない ({r["fixed_fails"]})')
r = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', extra={'time.nStepInner': 5}), 'bell_tp_ns_sst', 'main')
check([k for k, _ in r['fixed_fails']] == ['time.nStepInner'] and set(r['fixed_fails']) <= set(r['fails']),
      f'fixed_fails = fixed の値の不一致 (関門の FORGE_LAUNCH_GATE=warn が WARN に落とす分) ({r["fixed_fails"]})')

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
check(fail_keys(r2) == {'time.deltaT.implicitRelax'} and any('書いていない' in m for k, m in r2['fails']),
      f'プロファイルの fixed を書いていない (プロファイル 0.7) は名指しで FAIL ({sorted(fail_keys(r2))})')
r3 = RS.check_config(RC, cfg_of('bell_tp_ns_sst', 'main', edit=lambda d: d['time']['deltaT'].pop('implicitRelax')),
                     'bell_tp_ns_sst', 'main', deviations={'time.deltaT.implicitRelax': None})
check(not r3['fails'] and r3['applied'] == [{'key': 'time.deltaT.implicitRelax', 'value': None, 'profile_value': 0.7,
                                             'effective': RS.UNDETERMINED}],
      f'「書かない」(~) の逸脱を宣言すれば通り、実効値は未確定として記録 ({r3["applied"]})')

# --- #16: 省略の FAIL (solver の既定が同じでも) と明示すべきキー ------------------------------------------------------------
# 判別 A/B (2026-10-15 codex diagnose、forge 0 step): SERN 2D 本段のプロファイル・config で A = space.limiterScaled: 1 を明記、
# B = そのキーだけ削除。HEAD 79feaefb (#16 の前) では A・B とも FAIL/WARN 0 件だった (solver の既定 1 = プロファイル 1 で省略を通した)
SERN = 'sern2d_finite_base_frozen_tp_sst'
check(RS.resolve(RC, SERN)['stages']['main']['fixed'].get('space.limiterScaled') == 1
      and RS.solver_defaults().get('space.limiterScaled') == 1, '(A/B の前提: プロファイル 1、solver の既定も 1)')
rA = RS.check_config(RC, cfg_of(SERN, 'main'), SERN, 'main')
check(not rA['fails'] and not rA['warns'], f'A/B の A (limiterScaled: 1 を明記): PASS (fails {sorted(fail_keys(rA))})')
rB = RS.check_config(RC, cfg_of(SERN, 'main', edit=lambda d: d['space'].pop('limiterScaled')), SERN, 'main')
check(fail_keys(rB) == {'space.limiterScaled'} and not rB['warns'] and any('書いていない' in m for k, m in rB['fails']),
      f'A/B の B (limiterScaled だけ削除): 同じキーを名指しで FAIL ({sorted(fail_keys(rB))})')
r = RS.check_config(RC, cfg_of(SERN, 'main', edit=lambda d: d['space'].pop('limiterScaled')), SERN, 'main',
                    deviations={'space.limiterScaled': None})
check(not r['fails'] and r['applied'] and r['applied'][0].get('effective') == RS.UNDETERMINED,
      'B も space.limiterScaled=~ の逸脱を宣言すれば通り、実効値は未確定')
r = RS.check_config(RC, cfg_of('common_node_ns_sst', 'main', edit=lambda d: d['turbulence'].pop('scalarDiffusion')),
                    'common_node_ns_sst', 'main')
check(fail_keys(r) == {'turbulence.scalarDiffusion'},
      f'共通プロファイルでも fixed の省略は FAIL (scalarDiffusion: solver の既定 1 = プロファイル 1) ({sorted(fail_keys(r))})')
# 明示すべきキー (explicit) の省略・違う値・逸脱
ex_main = RS.resolve(RC, SERN)['stages']['main']['explicit']
check({'mesh.scalarGradient', 'time.deltaT.implicitRelax', 'turbulence.turbulentPrandtl'} <= set(ex_main),
      f'SERN 2D 本段の明示すべきキーに scalarGradient・implicitRelax・turbulentPrandtl がある ({sorted(ex_main)})')
r = RS.check_config(RC, cfg_of(SERN, 'main', edit=lambda d: d['mesh'].pop('scalarGradient')), SERN, 'main')
check(fail_keys(r) == {'mesh.scalarGradient'} and any('明示すべきキー' in m for k, m in r['fails']),
      f'明示すべきキーの省略は名指しで FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, cfg_of(SERN, 'main', extra={'mesh.scalarGradient': 'gg'}), SERN, 'main')
check(fail_keys(r) == {'mesh.scalarGradient'}, f'明示すべきキーの値が旧い実効値と違えば FAIL ({sorted(fail_keys(r))})')
r = RS.check_config(RC, cfg_of(SERN, 'main', extra={'mesh.scalarGradient': 'gg'}), SERN, 'main',
                    deviations={'mesh.scalarGradient': 'gg'})
check(not r['fails'] and r['applied'] == [{'key': 'mesh.scalarGradient', 'value': 'gg', 'profile_value': 'lsq'}],
      '明示すべきキーも理由つきの逸脱で外せる')
# 今の生成器と同じく明示すべきキーを省いた config (TP の温度反転の条件 thermoHrefTemp は問題のキーなので足す)
gen_like = cfg_of(SERN, 'main', explicit=False, extra={'physProp.thermoHrefTemp': 298.15})
r = RS.check_config(RC, gen_like, SERN, 'main')
check(fail_keys(r) == set(ex_main) and not r['warns'],
      f'fixed だけの config (今の生成器と同じ) は明示すべきキーの省略だけで FAIL ({len(fail_keys(r))} キー = 例外表 {len(ex_main)})')
BC_ISO = {'inlet': {'kind': 'inlet_Pressure'}, 'wall': {'kind': 'wall_isothermal'}, 'outlet': {'kind': 'outflow'}}
check(RS.explicit_for(RC, gen_like, BC_ISO) == ex_main, 'explicit_for (書き写しの台本の計算) = プロファイルの explicit')
d = run_dir('omit_decl', cfg_of(SERN, 'main', edit=lambda d: d['space'].pop('limiterScaled')))
rc, r, path = RS.declare(d, SERN, 'main', deviate=['space.limiterScaled=~'], reason='省略の試し', recipes=RC)
rec = json.load(open(path)) if path else {}
check(rc == 0 and rec.get('deviations') and rec['deviations'][0]['key'] == 'space.limiterScaled'
      and rec['deviations'][0]['value'] is None and rec['deviations'][0]['effective'] == RS.UNDETERMINED,
      f'RECIPE.json: 「書かない」の逸脱は実効値「未確定」で記録 ({rec.get("deviations")})')
# when の条件
ctx = {'a': 1, 'b': 'SLAU'}
check(RS.when_match({'a': 1, 'b': ['slau2', 'SLAU']}, ctx) and RS.when_match({'c': [None, 0]}, ctx)
      and not RS.when_match({'c': 0}, ctx) and RS.when_match({'c': {'not': [1]}}, ctx) and not RS.when_match({'a': {'not': [1]}}, ctx)
      and RS.when_match({'a': {'gt': 0}}, ctx) and not RS.when_match({'c': {'gt': 0}}, ctx) and RS.when_match({'a': {'min': 1}}, ctx),
      'when: 値・列・null (書いていない)・not・gt・min')
# 旧い実効値の出所の検査 (ソルバの既定が変わったら止まる)
check(len(RS.explicit_table(RC)) >= 20, f'explicit_keys の出所が今のソースと合う ({len(RS.explicit_table(RC))} キー)')
for tag, old, new in (('inventory', "source: {inventory: '1.0', at: solverConfig.cpp:370}", "source: {inventory: '0.9', at: solverConfig.cpp:370}"),
                      ('file', "text: 'int nodeWallDirichlet = 1;'", "text: 'int nodeWallDirichlet = 0;'")):
    txt0 = open(RS.RECIPES, encoding='utf-8').read()
    check(txt0.count(old) == 1, f'(出所の負例の入力を作れた: {tag})')
    try:
        RS.explicit_table(RS.load_recipes(text=txt0.replace(old, new)))
        check(False, f'出所 ({tag}) が今のソースと合わなければ例外')
    except RS.RecipeError as e:
        check('変わった' in str(e), f'出所 ({tag}) が今のソースと合わなければ例外 ({str(e)[:80]}…)')

# --- #21: 旧キー physProp.isAxisymmetric (solver と同じ優先順: mesh が優先、無ければ旧 physProp) ----------------------------
# 2026-10-15 までは derive_applies が mesh.isAxisymmetric だけを見て、旧キーだけの軸対称の config を平面と読んでいた
for tag, doc, want in (('mesh だけ 1', {'mesh': {'isAxisymmetric': 1}}, True),
                       ('旧 physProp だけ 1 (負例: 旧実装は平面)', {'physProp': {'isAxisymmetric': 1}}, True),
                       ('mesh 0・physProp 1 (mesh が優先)', {'mesh': {'isAxisymmetric': 0}, 'physProp': {'isAxisymmetric': 1}}, False),
                       ('mesh 1・physProp 0 (mesh が優先)', {'mesh': {'isAxisymmetric': 1}, 'physProp': {'isAxisymmetric': 0}}, True),
                       ('どちらも無し', {'mesh': {}}, False)):
    got = RS.derive_applies(doc)['axisymmetric']
    check(got is want, f'derive_applies の軸対称: {tag} → {want} (得た {got})')
AXI = 'axismach_ns_sst_tp'
BC_AXI = {'inlet': {'kind': 'inlet_Pressure'}, 'wall': {'kind': 'wall_isothermal'}, 'axis': {'kind': 'axis'}}


def to_legacy_axi(d):
    d['mesh'].pop('isAxisymmetric')
    d['mesh'].pop('axisymMethod', None)
    d['physProp']['isAxisymmetric'] = 1


leg = cfg_of(AXI, 'main', explicit=False, edit=to_legacy_axi)
ex_leg = RS.explicit_for(RC, leg, BC_AXI)
check('mesh.axisymMethod' in ex_leg and 'mesh.axisSegmentRWeight' in ex_leg,
      f'旧キーだけの軸対称の config でも明示すべきキーの軸対称の条件 (applies.axisymmetric) が当たる ({sorted(k for k in ex_leg if "axis" in k)})')
r = RS.check_config(RC, cfg_of(AXI, 'main', edit=to_legacy_axi), AXI, 'main')
check({'mesh.isAxisymmetric', 'physProp.isAxisymmetric'} <= fail_keys(r) and 'physProp.isAxisymmetric' in warn_keys(r),
      f'軸対称の機種で mesh.isAxisymmetric を旧キーに替えると、mesh 側の省略とプロファイルに無い旧キーで FAIL、旧キーは WARN '
      f'({sorted(fail_keys(r))}・warns {sorted(warn_keys(r))})')
r = RS.check_config(RC, cfg_of(SERN, 'main', edit=lambda d: d['physProp'].__setitem__('isAxisymmetric', 1)), SERN, 'main')
check(fail_keys(r) == {'physProp.isAxisymmetric'} and any('効かない' in m for k, m in r['warns'] if k == 'physProp.isAxisymmetric'),
      f'平面の機種に旧キー 1 を足しても mesh 0 が優先 (適用条件は平面のまま)、旧キーはプロファイルに無いキーで FAIL・効かない旨の WARN '
      f'({sorted(fail_keys(r))})')

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
ok, msg = raises(BASE + '  c: {inherits: b, stages: {main: {explicit: {time.nStepInner: 4}}}}\n', '継承の衝突')
check(ok, f'鎖の中で同じキーに違う種類 (fixed / explicit) は例外 ({msg})')
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
