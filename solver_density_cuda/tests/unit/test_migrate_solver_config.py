#!/usr/bin/env python3
"""migrate_solver_config.py の試験 (plan config-key-pruning §4.1 S2)。

重点は **節を取り違えないこと**: `time.last.control` を消すつもりで `time.deltaT.control` を消すと、
必須キーが消えて起動しなくなる。フロー形・ブロック形・コメント付きの各形を見る。
`time.last.time` の削除と `physProp.isAxisymmetric` → `mesh.isAxisymmetric` の移送 (plan tooling-generator-defaults-and-launch-gate
§5.1 #25) も同じ形で見る。
"""
import importlib.util, os, sys, yaml

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('msc', os.path.join(here, '..', '..', 'tools', 'migrate_solver_config.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print('  FAIL:', msg)

def roundtrip(text):
    out, status = m.migrate_one(text)
    assert out is not None, f'unexpectedly refused: {status}'
    y = yaml.safe_load(text) or {}
    drops, move = m.plan_for(y)
    return yaml.safe_load(out) or {}, out, drops, move

BLOCK = """mesh:
  meshFormat: "hdf5"
  discretization: "node"
  meshFileName: "a.h5"
physProp:
  isCompressible: 1   # 圧縮性
  thermalMethod: 0
  ro : 1.2
  visc : 1.8e-5
time:
  last:
    control: 0
    nStepOuter: 60000
  deltaT:
    control: 1
    cfl: 20.0
"""
got, out, drops, move = roundtrip(BLOCK)
check(got['time']['deltaT']['control'] == 1, 'time.deltaT.control must survive (block form)')
check('control' not in got['time']['last'], 'time.last.control must be removed')
check(got['time']['last']['nStepOuter'] == 60000, 'nStepOuter must survive')
check('meshFormat' not in got['mesh'] and got['mesh']['discretization'] == 'node', 'mesh.meshFormat removed, siblings kept')
check('isCompressible' not in got['physProp'] and 'ro' not in got['physProp'], 'physProp.isCompressible and ro removed')
check(got['physProp']['visc'] == 1.8e-5 and got['physProp']['thermalMethod'] == 0, 'physProp siblings kept')

FLOW = """mesh: {meshFormat: "hdf5", discretization: "node", meshFileName: "a.h5", valueFileName: "a.h5"}
gpu: 1
physProp: {isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, cp: 1004.5}
time:
  last: {control: 0, nStepOuter: 200}
  deltaT: {control: 1, cfl: 2.0, detectNaNInterval: 7}
"""
got, out, drops, move = roundtrip(FLOW)
check(got['mesh'] == {'discretization': 'node', 'meshFileName': 'a.h5', 'valueFileName': 'a.h5'}, f"flow mesh: {got['mesh']}")
check(got['physProp'] == {'thermalMethod': 0, 'viscMethod': 1, 'visc': 1.8e-5, 'cp': 1004.5}, f"flow physProp: {got['physProp']}")
check(got['time']['last'] == {'nStepOuter': 200}, f"flow time.last: {got['time']['last']}")
check(got['time']['deltaT']['control'] == 1 and got['time']['deltaT']['cfl'] == 2.0, 'time.deltaT survives the flow edit')
check('detectNaNInterval' not in got['time']['deltaT'], 'the alias is removed')
check(got.get('detectNaNInterval') == 7, f'the alias value moves to the top level, got {got.get("detectNaNInterval")}')
check(got['gpu'] == 1, 'unrelated top-level keys survive')

# ブロック形で**末尾コメントに , が入る**ケース (実配置で 599 本が壊れた形)
COMMENTED = """mesh:
  meshFormat: "hdf5"
  meshFileName: "cavity.msh"
physProp:
  isCompressible: 0 # 0:incomp, 1:compressible
  ro : 1.2 # [kg/m3]
  visc : 0.000018 #[Pa*s]
time:
  last:
    control: 0 # 0:nStep , 1:time base
    nStepOuter: 1000
"""
got, out, drops, move = roundtrip(COMMENTED)
check(got == {'mesh': {'meshFileName': 'cavity.msh'},
              'physProp': {'visc': 1.8e-5},
              'time': {'last': {'nStepOuter': 1000}}}, f'block form with trailing commas in comments: {got}')
check('1:compressible' not in out and '1:time base' not in out, 'no comment fragment is left behind')

# トップレベルに既にあるなら、別名は捨ててよい (トップレベルが優先されるため)
BOTH = FLOW + "detectNaNInterval: 100\n"
got, out, drops, move = roundtrip(BOTH)
check(got.get('detectNaNInterval') == 100, f'the top-level value wins, got {got.get("detectNaNInterval")}')
check(out.count('detectNaNInterval') == 1, 'the alias line is gone and the top-level one is not duplicated')

# 何も無い config は触らない
NONE = "gpu: 1\ntime:\n  deltaT: {control: 1, cfl: 2.0}\n"
y = yaml.safe_load(NONE); drops, move = m.plan_for(y)
check(drops == [] and move is None, 'a config without staged keys needs no edit')

# --- 自動移送を拒否すべき入力 (codex plan-5 M1/M2/m6) ---
ANCHOR = """time:
  last: &shared
    control: 1
    nStepOuter: 200
  deltaT: *shared
"""
out, status = m.migrate_one(ANCHOR)
check(out is None and status.startswith('skip:') and 'anchor' in status,
      f'shared anchors must be refused, not silently dropped from both sections: {status}')

DUP = """time:
  deltaT:
    cfl: 1.0
    cfl: 0.5
  last:
    control: 0
    nStepOuter: 10
"""
out, status = m.migrate_one(DUP)
check(out is None and '重複' in status, f'duplicate keys must be refused (yaml-cpp takes the first): {status}')

GMSH = """mesh:
  meshFormat: "gmsh"
  meshFileName: "a.msh"
physProp:
  ro: 1.2
"""
out, status = m.migrate_one(GMSH)
check(out is None and 'meshFormat' in status, f'a non-hdf5 mesh format must be reported, not dropped: {status}')

# 移送は冪等: 2 回目は「既に清潔」
out, status = m.migrate_one(BLOCK)
out2, status2 = m.migrate_one(out)
check(status2 == 'clean' and out2 is None, f'a second pass must find nothing to do, got {status2}')

# トップレベルに null で存在する場合も「ある」扱い (値検証と存在判定を分ける)
NULLTOP = """detectNaNInterval:
time:
  deltaT: {control: 1, detectNaNInterval: 7}
"""
out, status = m.migrate_one(NULLTOP)
got = yaml.safe_load(out)
check(got.get('detectNaNInterval') is None and 'detectNaNInterval' not in got['time']['deltaT'],
      f'a null top-level key still wins; no duplicate is created: {got}')
check(out.count('detectNaNInterval') == 1, 'no duplicate top-level key is appended')

# --- time.last.time と physProp.isAxisymmetric (plan tooling-generator-defaults-and-launch-gate §5.1 #25) ---
LASTTIME = """time:
  last:
    control: 0 # 0:nStep , 1:time base
    nStepOuter: 10000
    time: 1.00 # [second]
  deltaT:
    control: 1
    cfl: 2.0
"""
got, out, drops, move = roundtrip(LASTTIME)
check(got['time'] == {'last': {'nStepOuter': 10000}, 'deltaT': {'control': 1, 'cfl': 2.0}},
      f'block: time.last.time and time.last.control removed, time.deltaT.control kept: {got}')
check('[second]' not in out and 'deltaT:' in out, 'the comment of time.last.time goes with its line')
got, out, drops, move = roundtrip('time:\n  last: {nStepOuter: 10000, time: 1.00}\n  deltaT: {control: 1, cfl: 2.0}\n')
check(got['time'] == {'last': {'nStepOuter': 10000}, 'deltaT': {'control': 1, 'cfl': 2.0}}, f'flow: time.last.time removed: {got}')
out, status = m.migrate_one('time:\n  last: {nStepOuter: 10, time: 1.0}\n')
out2, status2 = m.migrate_one(out)
check(status == 'ok' and status2 == 'clean', f'time.last.time removal is a plain drop and idempotent: {status} / {status2}')

AXB = """mesh:
  meshFileName: "a.h5"
  discretization: "node"   # node
physProp:
  isAxisymmetric: 1   # 旧キー
  thermalMethod: 0
time:
  last: {nStepOuter: 10}
"""
out, status = m.migrate_one(AXB)
got = yaml.safe_load(out)
check(status == 'moved' and got['mesh'] == {'meshFileName': 'a.h5', 'discretization': 'node', 'isAxisymmetric': 1}
      and got['physProp'] == {'thermalMethod': 0}, f'block: physProp.isAxisymmetric moves to mesh with its value: {status} {got}')
check(out.count('isAxisymmetric') == 1 and '  isAxisymmetric: 1\nphysProp:' in out, f'inserted once, at the end of the mesh block:\n{out}')
AXF = 'mesh: {meshFileName: "a.h5", discretization: "node"}\nphysProp: {thermalMethod: 0, isAxisymmetric: 1, cp: 1004.5}\n'
out, status = m.migrate_one(AXF)
got = yaml.safe_load(out)
check(status == 'moved' and got == {'mesh': {'meshFileName': 'a.h5', 'discretization': 'node', 'isAxisymmetric': 1},
                                    'physProp': {'thermalMethod': 0, 'cp': 1004.5}}, f'flow: moved: {got}')
out, status = m.migrate_one('mesh: {meshFileName: "a.h5", isAxisymmetric: 0}\nphysProp: {isAxisymmetric: 1, thermalMethod: 0}\n')
got = yaml.safe_load(out)
check(status == 'ok' and got == {'mesh': {'meshFileName': 'a.h5', 'isAxisymmetric': 0}, 'physProp': {'thermalMethod': 0}},
      f'mesh side present: only the old key goes (the effective value stays mesh 0): {status} {got}')
out, status = m.migrate_one('mesh: {meshFileName: "a.h5"}\nphysProp: {isAxisymmetric: , thermalMethod: 0}\n')
check(out is None and status.startswith('skip:') and '空' in status, f'an empty old key is reported, not moved: {status}')
out, status = m.migrate_one('mesh: {meshFileName: "a.h5", isAxisymmetric: }\nphysProp: {isAxisymmetric: 1}\n')
check(out is None and status.startswith('skip:'), f'an empty mesh key is reported: {status}')
out, status = m.migrate_one('physProp: {isAxisymmetric: 1, thermalMethod: 0}\n')
check(out is None and status.startswith('skip:') and 'mesh 節' in status, f'no mesh section: reported: {status}')
out, status = m.migrate_one(AXB)
out2, status2 = m.migrate_one(out)
check(status2 == 'clean' and out2 is None, f'the axisymmetric move is idempotent: {status2}')
hp = m.handled_paths()
check(hp.get('time.last.time') == '消す' and 'mesh.isAxisymmetric' in hp.get('physProp.isAxisymmetric', ''),
      f'handled_paths lists the new keys: {hp}')

print('ALL PASS' if fails == 0 else f'FAILED ({fails})')
sys.exit(1 if fails else 0)
