#!/usr/bin/env python3
"""check_solver_config.py の試験: 残差では気づけない設定ミス (dual-time の緩和・不活性な S3・作動しない FCT・存在しないキー) を FAIL/WARN で止めること。"""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('csc', os.path.join(here, '..', '..', 'tools', 'check_solver_config.py'))
csc = importlib.util.module_from_spec(spec); spec.loader.exec_module(csc)

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)

def cfg(**kw):
    dT = {'dt': 8e-6, 'cfl_pseudo': kw.pop('cflp', 12.0), 'passiveScalarScheme': kw.pop('scheme', 1),
          'passiveFct': kw.pop('fct', 1), 'speciesFaceReconstruction': kw.pop('sfr', 2)}
    for k in ('implicitRelax', 'speciesImplicitCoupling'):
        if k in kw: dT[k] = kw.pop(k)
    y = {'solver': kw.pop('solver', 'SLAU'), 'mesh': {}, 'space': {'convMethod': kw.pop('conv', 1), 'limiter': 1},
         'physProp': {'species': ['A', 'B']}, 'condensation': {'condensation': kw.pop('cond', 1), 'nCondSpecies': 1},
         'time': {'unsteady': kw.pop('unsteady', 1), 'dualTime': kw.pop('dualTime', 1), 'timeIntegration': kw.pop('ti', 11),
                  'nSubIterDualTime': kw.pop('nsub', 20), 'last': {'nStepOuter': 100}, 'deltaT': dT}}
    for k, v in kw.items(): y[k] = v
    return y

def keys(res): return {k for k, _ in res}

f, w = csc.check(cfg())
check(not f and not w, f'recommended dual-time settings must be clean (fails {keys(f)}, warns {keys(w)})')

f, w = csc.check(cfg(implicitRelax=0.7))
check('time.deltaT.implicitRelax' in keys(f), 'implicitRelax < 1 in dual-time must FAIL')

f, w = csc.check(cfg(implicitRelax=0.7, unsteady=0, dualTime=0, ti=4))
check('time.deltaT.implicitRelax' not in keys(f), 'implicitRelax 0.7 in a steady run must be accepted')

f, w = csc.check(cfg(cflp=2.0))
check('time.deltaT.cfl_pseudo' in keys(w), 'small cfl_pseudo in dual-time must WARN')

f, w = csc.check(cfg(nsub=5))
check('time.nSubIterDualTime' in keys(w), 'too few sub-iterations must WARN')

f, w = csc.check(cfg(conv=0))
check('space.convMethod' in keys(f), 'speciesFaceReconstruction 2 with convMethod 0 must FAIL (S3 inactive)')

f, w = csc.check(cfg(scheme=0))
check('time.deltaT.passiveScalarScheme' in keys(f), 'condensation dual-time with scheme 0 must FAIL (no BDF term on the moments)')

f, w = csc.check(cfg(solver='KEEP'))
check('time.deltaT.passiveFct' in keys(w), 'passiveFct 1 with a non-SLAU solver must WARN (never activates)')

f, w = csc.check(cfg(sfr=0, fct=1))
check('time.deltaT.passiveFct' in keys(w), 'passiveFct 1 with speciesFaceReconstruction 0 must WARN')

# passiveFct の省略 (既定 1) は「作動する構成なら使う」既定: 定常・RK の通常の config で作動しないのは正しい動作なので警告しない。
# 明示した有効化 (passiveFct: 1) が作動しないときだけ警告する (2026-10-15 codex diagnose 案 5、plan §5.1 #10)。
def drop_fct(y):
    del y['time']['deltaT']['passiveFct']
    return y

for name, y in (('steady (unsteady 0)', cfg(unsteady=0, dualTime=0, ti=4)), ('RK (dualTime 0)', cfg(dualTime=0, ti=3)),
                ('dual-time with KEEP', cfg(solver='KEEP')), ('dual-time with SFR 0', cfg(sfr=0))):
    f, w = csc.check(drop_fct(y))
    check('time.deltaT.passiveFct' not in keys(w), f'omitted passiveFct ({name}) must not WARN (default, not a request)')
f, w = csc.check(cfg(unsteady=0, dualTime=0, ti=4, fct=1))
check('time.deltaT.passiveFct' in keys(w), 'explicit passiveFct 1 in a steady run must still WARN (the request never activates)')
f, w = csc.check(cfg(unsteady=0, dualTime=0, ti=4, fct=0))
check('time.deltaT.passiveFct' not in keys(w), 'explicit passiveFct 0 must not WARN')

f, w = csc.check(cfg(speciesImplicitCoupling=0))
check('time.deltaT.speciesImplicitCoupling' in keys(w), 'S3 with coupling 0 must WARN')

y = cfg(); y['turbulence'] = {'model': 'SST', 'kInf': 1.0, 'omegaInf': 1000.0}
f, w = csc.check(y)
check('turbulence.kInf' in keys(w) and 'turbulence.omegaInf' in keys(w), 'keys the solver never reads must WARN (kInf/omegaInf are not real keys)')

y = cfg(); y['turbulence'] = {'model': 'SST', 'kInit': 1.0, 'omegaInit': 1000.0}
f, w = csc.check(y)
check('turbulence.kInit' not in keys(w), 'the real initial-value keys must not be flagged')

# 誤配置の検出 (codex plan-2 M5 の実例)。末端名だけの照合では全部素通りしていた。
for path, val in [('lowMachPrecond', 2), ('physProp.speciesFaceReconstruction', 2), ('space.keepDissType', 1)]:
    y = cfg(); node = y
    parts = path.split('.')
    for seg in parts[:-1]: node = node.setdefault(seg, {})
    node[parts[-1]] = val
    f, w = csc.check(y)
    check(path in keys(w), f'misplaced key {path} must WARN (silently ignored where it is written)')

# 起動時に拒否されるキーは WARN でなく FAIL
y = cfg(); y['mesh']['gradLSQDegenThresh'] = 1e-6
f, w = csc.check(y)
check('mesh.gradLSQDegenThresh' in keys(f), 'a key the solver rejects at startup must FAIL')

# 読まれないが必須として全 run が書いているキーは、誤検出しない
y = cfg(); y['time']['last']['time'] = 1.0
f, w = csc.check(y)
check('time.last.time' in keys(w), 'time.last.time is not read by the solver and must WARN')
# migrate_solver_config.py が直せるキーの WARN・FAIL はそのツールを案内する (plan tooling-generator-defaults-and-launch-gate §5.1 #25)
wm = dict(w)
check('migrate_solver_config.py' in wm.get('time.last.time', '') and '(消す)' in wm.get('time.last.time', ''),
      f'the time.last.time WARN points to migrate_solver_config.py ({wm.get("time.last.time")})')
y = cfg(); y['time']['last']['control'] = 0
f, w = csc.check(y)
check('migrate_solver_config.py' in dict(f).get('time.last.control', ''), f'the time.last.control FAIL points to migrate ({dict(f).get("time.last.control")})')
y = cfg(); y['physProp']['isAxisymmetric'] = 1
f, w = csc.check(y)
check('migrate_solver_config.py' in dict(w).get('physProp.isAxisymmetric', '') and 'mesh.isAxisymmetric へ移す' in dict(w)['physProp.isAxisymmetric'],
      f'the legacy isAxisymmetric WARN points to migrate ({dict(w).get("physProp.isAxisymmetric")})')
y = cfg(); y['time']['last']['tiem'] = 1.0
f, w = csc.check(y)
check('migrate_solver_config.py' not in dict(w).get('time.last.tiem', 'x'), 'keys migrate does not handle get no migrate hint')

# 正しい位置のキーは素通りする (偽陽性の回帰)
y = cfg(); y['mesh']['primPack'] = 1; y['time']['deltaT']['blockDPLURDqPack'] = 1
f, w = csc.check(y)
check(not keys(f) and 'mesh.primPack' not in keys(w), 'the restored opt-in performance switches must be accepted')

y = cfg(); y['mesh']['bndFirstOrder'] = 1
f, w = csc.check(y)
check('mesh.bndFirstOrder' in keys(f), 'bndFirstOrder must FAIL (banned)')

# 起動時に拒否されるキー (if の条件で参照して "no longer supported" で落とす書き方) は FAIL
y = cfg(); y['space']['limiterMatchRecon'] = 1
f, w = csc.check(y)
check('space.limiterMatchRecon' in keys(f), 'limiterMatchRecon is rejected by the solver at startup and must FAIL')

# 旧キー physProp.isAxisymmetric (plan tooling-generator-defaults-and-launch-gate §5.1 #21): solver は deprecated として読む
# (mesh 側が優先: solverConfig.cpp:945・:953) ので FAIL にしない。WARN で移行を促し、mesh と両方あれば physProp 側が効かないと出す
for tag, mesh, pp, want in (('mesh only', {'isAxisymmetric': 1}, {}, 1), ('legacy only', {}, {'isAxisymmetric': 1}, 1),
                            ('mesh 0 wins over legacy 1', {'isAxisymmetric': 0}, {'isAxisymmetric': 1}, 0),
                            ('mesh 1 wins over legacy 0', {'isAxisymmetric': 1}, {'isAxisymmetric': 0}, 1),
                            ('neither', {}, {}, 0), ('legacy as a string', {}, {'isAxisymmetric': '1'}, 1),
                            ('empty mesh value falls back', {'isAxisymmetric': None}, {'isAxisymmetric': 1}, 1)):
    got = csc.effective_axisymmetric({'mesh': mesh, 'physProp': pp})
    check(got == want, f'effective_axisymmetric ({tag}): expected {want}, got {got}')
y = cfg(); y['physProp']['isAxisymmetric'] = 1
f, w = csc.check(y); wk = dict(w)
check('physProp.isAxisymmetric' not in keys(f) and 'physProp.isAxisymmetric' in wk and 'mesh.isAxisymmetric' in wk['physProp.isAxisymmetric'],
      f'legacy physProp.isAxisymmetric alone is accepted with a WARN that names mesh.isAxisymmetric ({wk.get("physProp.isAxisymmetric")})')
check('physProp.isAxisymmetric' not in {k for k, _ in csc.unknown_keys(y)[1]}, 'the legacy key is read by the solver, so it is not an unknown key')
y = cfg(); y['physProp']['isAxisymmetric'] = 1; y['mesh']['isAxisymmetric'] = 0
f, w = csc.check(y); wk = dict(w)
check('physProp.isAxisymmetric' in wk and '効かない' in wk['physProp.isAxisymmetric'] and '食い違う' in wk['physProp.isAxisymmetric'],
      f'legacy key overridden by a different mesh value must say it has no effect ({wk.get("physProp.isAxisymmetric")})')
y = cfg(); y['mesh']['isAxisymmetric'] = 1
f, w = csc.check(y)
check('physProp.isAxisymmetric' not in keys(w) and 'mesh.isAxisymmetric' not in keys(w), 'mesh.isAxisymmetric alone is clean')

# ---- 受理キーの一覧の完全性 (plan tooling-generator-defaults-and-launch-gate §5.1 #4、codex plan 段 M1) ----
# 抽出は手書きの対応表でなく、ソースの束縛 (`auto cj = config["conjugate"];`) から節を導く。
import re, yaml
root = os.path.abspath(os.path.join(here, '..', '..', '..'))
spec2 = importlib.util.spec_from_file_location('cki', os.path.join(here, '..', '..', 'tools', 'config_key_inventory.py'))
cki = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(cki)
cpp_raw = open(cki.SRC_CPP, errors='replace').read()
cpp_src = cki.strip_comments(cpp_raw)

def flagged(y):
    f, w = csc.unknown_keys(y)
    return {k for k, _ in f} | {k for k, _ in w}

# (a) CHT: ソルバが受理する conjugate.* (solverConfig.cpp の allowedCj・allowedGate) は一つも誤検出しない
m_cj = re.search(r'allowedCj\s*=\s*\{([^}]*)\}', cpp_src)
m_gate = re.search(r'allowedGate\s*=\s*\{([^}]*)\}', cpp_src)
check(m_cj and m_gate, 'allowedCj / allowedGate whitelists not found in solverConfig.cpp')
allowed_cj = re.findall(r'"([A-Za-z0-9_]+)"', m_cj.group(1)) if m_cj else []
allowed_gate = re.findall(r'"([A-Za-z0-9_]+)"', m_gate.group(1)) if m_gate else []
values, sections, _ = csc._spec()
for k in allowed_cj:
    p = 'conjugate.' + k
    check(p in values or p in sections, f'{p} is accepted by the solver (allowedCj) but missing from the extracted list')
for k in allowed_gate:
    p = 'conjugate.gate.' + k
    check(p in values, f'{p} is accepted by the solver (allowedGate) but missing from the extracted list')
for wrong in ('cj.mode', 'cj.flux', 'cj.back'):
    check(wrong not in values, f'{wrong} is an artefact of the old variable table and must not be in the list')
y = cfg(); y['conjugate'] = {'mode': 'fem2d', 'solid': 'solid.h5', 'flux': 'q_eff', 'gate': {'eps_rel': 0.01, 'eps_abs_Wm2': 150.0, 'dT_K': 0.01, 'n_consec': 80, 'tol_solid': 1e-6}}
check(not flagged(y), f'valid fem2d conjugate block must not be flagged ({flagged(y)})')
# local1d の実例 (/home/sano/work/forge の case/52.conjugate_slab/run_0900_absorb、追跡外なので値を写す)
y = cfg(); y['conjugate'] = yaml.safe_load('{mode: local1d, thickness: 2.0e-3, k_solid: 0.01, back: isothermal, T_b: 300.0, interval: 50, warmup: 0}')
check(not flagged(y), f'valid local1d conjugate block must not be flagged ({flagged(y)})')

# (b) 実在の正当な config (追跡済み・各種別の最近のもの)。旗が立つのは、ソルバが読まないキー (下の明示の集合) だけ。
#     読まないことは、末端名の文字列リテラルがソルバのソースに無いことで独立に確かめる。
POSITIVE = [
    # (種別, パス, 期待する機能, ソルバが読まない既知のキー)
    ('CHT fem2d+gate', 'case/53.c3x_vane_cht/_speed_on/solverConfig.yaml', {'conjugate', 'periodic'}, {'turbulence.kInf', 'turbulence.omegaInf'}),
    ('TP', 'case/46.sern_design/chi_evidence/46/0968_chidef_2d_m10_on_omit/solverConfig.yaml', {'tp'}, set()),
    ('TP+condensation', 'case/46.sern_design/chi_evidence/16/0964_chidef_b1c_omit/solverConfig.yaml', {'tp', 'cond'}, {'turbulence.kInf', 'turbulence.omegaInf'}),
    ('TP+condensation', 'case/44.vitiated_air_wt/run_0164_va3_M4.19_Lc8_noneq_inletTt_cell_lim1d_cfl2/solverConfig.yaml', {'tp', 'cond'}, set()),
    ('periodic', 'case/46.sern_design/chi_evidence/39/0953_chidef_b1p_omit/solverConfig.yaml', {'periodic'}, set()),
    ('periodic', 'case/09.Taylor-Green/_g0_lsq_seam/solverConfig.yaml', {'periodic'}, set()),
    ('dual-time TP', 'case/09.Taylor-Green/_g0_lsq_seam/s2_evidence/16/0961_sglsq_s2_fct_gg_nokey/solverConfig.yaml', {'tp', 'dual'}, set()),
    ('dual-time TP', 'case/09.Taylor-Green/_g0_lsq_seam/s2_evidence/16/0955_sglsq_s2_fct_lsq/solverConfig.yaml', {'tp', 'dual'}, set()),
]
for kind, rel, feats, ignored in POSITIVE:
    path = os.path.join(root, rel)
    if not os.path.exists(path):
        check(False, f'positive example missing: {rel}'); continue
    y = yaml.safe_load(open(path)) or {}
    pp = y.get('physProp') or {}; tt = y.get('time') or {}
    bp = os.path.join(os.path.dirname(path), 'bcondConfig.yaml')
    have = set()
    if 'conjugate' in y: have.add('conjugate')
    if pp.get('thermalMethod') == 2 and pp.get('species'): have.add('tp')
    if (y.get('condensation') or {}).get('condensation') == 1: have.add('cond')
    if tt.get('unsteady') == 1 and tt.get('dualTime') == 1: have.add('dual')
    if os.path.exists(bp) and 'periodic' in open(bp).read(): have.add('periodic')
    check(feats <= have, f'{rel}: expected features {sorted(feats)}, found {sorted(have)} (wrong example)')
    got = flagged(y)
    check(got == ignored, f'{kind} {rel}: flagged {sorted(got)}, expected only the solver-ignored {sorted(ignored)}')
    for p in ignored:
        check(f'"{p.split(".")[-1]}"' not in cpp_src, f'{p}: the leaf literal appears in solverConfig.cpp, so "not read" is not established')
    print(f'  positive {kind:16s} {rel}: flagged {sorted(got) or "-"}')

# (c) 負例: 存在しないキー・誤った節のキーは名指しで出る (今回は WARN のまま。FAIL 化は #10)
y = cfg(); y['conjugate'] = {'mode': 'fem2d', 'solid': 's.h5', 'modee': 'x', 'gate': {'eps_rell': 1e-3, 'eps_abs_Wm2': 1.0, 'dT_K': 1.0, 'n_consec': 1}}
f, w = csc.check(y); wk = dict(w)
check('conjugate.modee' in wk and 'conjugate.gate.eps_rell' in wk, f'non-existent conjugate keys must be named ({sorted(wk)})')
y = cfg(); y['conjugate'] = {'mode': 'fem2d', 'solid': 's.h5', 'eps_rel': 1e-3}
f, w = csc.check(y); wk = dict(w)
check('conjugate.eps_rel' in wk and 'conjugate.gate.eps_rel' in wk.get('conjugate.eps_rel', ''), f'misplaced gate key must name the right section ({wk.get("conjugate.eps_rel")})')
y = cfg(); y['time']['deltaT']['nStepInner'] = 4; y['physProp']['turbulentPrandtl'] = 0.9; y['notAKey'] = 1
f, w = csc.check(y); wk = dict(w)
check('time.deltaT.nStepInner' in wk and 'time.nStepInner' in wk['time.deltaT.nStepInner'], 'misplaced nStepInner must name time.nStepInner')
check('physProp.turbulentPrandtl' in wk and 'turbulence.turbulentPrandtl' in wk['physProp.turbulentPrandtl'], 'misplaced turbulentPrandtl must name turbulence.turbulentPrandtl')
check('notAKey' in wk, 'a top-level key the solver never reads must be named')

# (d) 抽出器そのもの: 束縛から節を導く (新しい別名・入れ子の別名・括弧を含む既定値・表示用の節名の誤り)
fake = cki.strip_comments(r"""
void f() {
    YAML::Node root = YAML::LoadFile(fn);
    auto sec = root["newSection"];
    const YAML::Node sub = sec["inner"];
    x = getOptionalValidatedValue<std::string>(sec, "mode", std::string("a,b"), "newSection");
    y = getValidatedValue<double>(sub, "eps", "WRONG.label");
    z = getOptionalValidatedValue<int>(root["time"], "n", (q == 1) ? 2 : 3, "time");
    if (sec["oldKey"]) { throw std::runtime_error("Key 'oldKey' is no longer supported."); }
    for (const char* k : {"gone1", "gone2"}) { if (sub[k]) { std::exit(1); } }
    for (const auto& e : sec["list"]) { w = e["name"]; }
}
""")
fk, frep = cki.extract_keys_with_report(fake)
for p in ('newSection.mode', 'newSection.inner.eps', 'time.n', 'newSection.list'):
    check(p in fk and fk[p]['kind'] == 'value', f'synthetic source: {p} must be extracted as a value ({fk.get(p)})')
check(fk.get('newSection.mode', {}).get('default') == '"a,b"', 'a default containing a comma inside std::string(...) must be kept whole')
for p in ('newSection.oldKey', 'newSection.inner.gone1', 'newSection.inner.gone2'):
    check(fk.get(p, {}).get('kind') == 'rejected', f'synthetic source: {p} must be rejected-only')
whys = {(r['why'], r['text'].split(' ')[0]) for r in frep}
check(any(w == 'label-mismatch' for w, _ in whys), 'a wrong display label must be reported (path comes from the binding)')
check(any(w == 'list-element' for w, _ in whys) and 'newSection.list.name' not in fk, 'subscripts on list elements are value contents, not config paths')

# (e) 実ソースで静的に解決できない参照は、値キーの中身 (list-element) と条件付きの拒否 (conditional-reject) だけ。
#     新しい書き方で解決できない参照が増えたらここで落ちる (一覧の取りこぼしを黙って生まない)。
_, rep = cki.extract_keys_with_report(cki.strip_comments(cpp_raw + '\n' + open(cki.SRC_HPP, errors='replace').read()))
bad = [r for r in rep if r['why'] not in ('list-element', 'conditional-reject')]
check(not bad, f'unresolvable references in solverConfig.cpp: {bad}')
import collections as _c
print('  static report:', dict(_c.Counter(r['why'] for r in rep)),
      ' conditional-reject:', sorted(r['text'] for r in rep if r['why'] == 'conditional-reject'))

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
