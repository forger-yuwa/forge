#!/usr/bin/env python3
"""solverConfig のキー棚卸し (plan config-key-pruning §5.1 #1): 定義・既定値と、case 配下の run config・procedures/methods での使われ方を突き合わせ、分類表の素を出す。\n使い方: config_key_inventory.py [--json OUT.json]"""
import re, glob, os, collections, json
ROOT='/home/sano/work/forge'
src=open(ROOT+'/solver_density_cuda/input/solverConfig.cpp').read()
# getOptionalValidatedValue<T>(sec, "key", default, ...) / getValidatedValue<T>(sec, "key", ...)
keys={}

def add(path, typ, dflt, required):
    """path は完全修飾 ('time.deltaT.dt' など)。同名キーを節ごとに別物として扱う (codex plan-1 M3)。"""
    e = keys.setdefault(path, {'path': path, 'sec': '.'.join(path.split('.')[:-1]) or '(top)', 'key': path.split('.')[-1],
                               'type': typ, 'default': (dflt or '').strip(), 'required': required, 'how': set()})
    if dflt and not e['default']: e['default'] = dflt.strip()
    return e

# 1) ヘルパ経由: get[Optional]ValidatedValue<T>(sec_var, "key", default, "section.path")
#    最後の引数が YAML 上の節パス (例 "time.deltaT")。無いものは変数名から引く。
VAR2PATH={'deltaT':'time.deltaT','last':'time.last','turb':'turbulence','cond':'condensation','ch':'physProp.chemistry',
          'physProp':'physProp','space':'space','out':'output','config':'(top)','time':'time','mesh':'mesh'}
for m in re.finditer(r'get(Optional)?ValidatedValue<([^>]+)>\(\s*(\w+)\s*,\s*"([A-Za-z0-9_]+)"\s*(?:,\s*([^,\)]+))?\s*(?:,\s*"([^"]*)")?', src):
    opt, typ, var, key, dflt, secstr = m.groups()
    sec = (secstr or VAR2PATH.get(var, var)).strip()
    path = key if sec in ('(top)', '') else f'{sec}.{key}'
    add(path, typ.strip(), dflt, not opt)['how'].add('helper')
# 2) 直接参照: config["a"]["b"]... / sec["key"]
for m in re.finditer(r'(config|[a-zA-Z_]\w*)((?:\[\s*"[A-Za-z0-9_]+"\s*\])+)', src):
    var, chain = m.groups()
    parts = re.findall(r'"([A-Za-z0-9_]+)"', chain)
    if var == 'config':
        path = '.'.join(parts)
    else:
        base = VAR2PATH.get(var)
        if base is None: continue
        path = ('' if base == '(top)' else base + '.') + '.'.join(parts)
    if not path or path.split('.')[0] not in ('mesh','time','space','turbulence','condensation','physProp','output','solver','gpu','initial','bodyForce','bodyForceCtrl','keepDissType','keepDissCoeff','keepDissJump','keepDissCprime','keepDissPrecond','keepDissFdBlend','keepDissCoeffMax','keepDissCbCoeff','keepDissCbEps','keepDissOpBlendRaw','bodyForceCtrlTarget','bodyForceCtrlRelax'):
        # top-level スカラは別に拾う
        if len(path.split('.')) == 1: add(path, '?', '', False)['how'].add('direct')
        continue
    add(path, '?', '', False)['how'].add('direct')

# run config での使用
used=collections.Counter(); vals=collections.defaultdict(set)
# 全 worktree を見る (指摘 2026-09-17: 本ツリーだけだと forge-cond / forge-chem の実使用を見落として誤判定する)
TREES=[ROOT]+[d for d in glob.glob('/home/sano/work/forge-*') if os.path.isdir(d+'/case')]
cfgs=[f for t in TREES for f in glob.glob(t+'/case/*/run_*/solverConfig.yaml')]
for p in cfgs:
    try: txt=open(p, errors='replace').read()
    except Exception: continue
    for k in {v['key'] for v in keys.values()}:
        for m in re.finditer(r'(?<![A-Za-z0-9_])'+re.escape(k)+r'\s*:\s*([^,}\n#]+)', txt):
            used[k]+=1
            vals[k].add(m.group(1).strip().strip('"\'').rstrip('}').strip())
# 文書での言及
docs=' '.join(open(f, errors='replace').read() for f in glob.glob(ROOT+'/procedures/*.md')+glob.glob(ROOT+'/methods/**/*.md', recursive=True))
rec=open(ROOT+'/procedures/recommended-settings.md', errors='replace').read()
rows=[]
for path, d in sorted(keys.items()):
    k = d['key']
    def same_as_default(v):
        d0=(d['default'] or '').strip().rstrip('f')
        if d0 == '': return False
        try: return abs(float(v)-float(d0)) <= 1e-12*max(1.0, abs(float(d0)))
        except ValueError: return v == d0
    nonz={v for v in vals[k] if not same_as_default(v)}
    rows.append({'path': path, 'key': k, 'sec': d['sec'], 'how': sorted(d['how']), 'default': d['default'], 'runs': used[k],
                 'values': sorted(vals[k])[:6], 'nondefault': len(nonz) > 0, 'nondefault_values': sorted(nonz)[:6],
                 'in_docs': k in docs, 'in_recommended': k in rec})
json.dump({'nconfig': len(cfgs), 'rows': rows}, open('/tmp/claude-1000/-home-sano-work-forge/96eabcfb-a7fa-40da-aae8-ef553208f9cd/scratchpad/key_inventory.json','w'), ensure_ascii=False, indent=1)
never=[r for r in rows if r['runs']==0]
docless=[r for r in rows if not r['in_docs']]
print(f'keys (完全修飾) {len(rows)} / run configs {len(cfgs)} (worktrees: ' + ', '.join(os.path.basename(t) for t in TREES) + ')')
print(f'どの run でも設定されていない: {len(never)}')
print(f'procedures・methods に一度も出てこない: {len(docless)}')
print(f'recommended-settings に無い: {len([r for r in rows if not r["in_recommended"]])}')
print('\n-- 削除候補 (run で未使用 かつ 文書に無い) --')
for r in never:
    if not r['in_docs']: print(f"  {r['sec']:10s} {r['key']:34s} default {r['default']}")
