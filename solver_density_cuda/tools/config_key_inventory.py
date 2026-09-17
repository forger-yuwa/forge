#!/usr/bin/env python3
"""solverConfig のキー棚卸し (plan config-key-pruning §5.1 #1): 定義・既定値と、case 配下の run config・procedures/methods での使われ方を突き合わせ、分類表の素を出す。\n使い方: config_key_inventory.py [--json OUT.json]"""
import re, glob, os, collections, json
ROOT='/home/sano/work/forge'
src=open(ROOT+'/solver_density_cuda/input/solverConfig.cpp').read()
# getOptionalValidatedValue<T>(sec, "key", default, ...) / getValidatedValue<T>(sec, "key", ...)
keys={}
for m in re.finditer(r'get(Optional)?ValidatedValue<([^>]+)>\(\s*(\w+)\s*,\s*"([A-Za-z0-9_]+)"\s*(?:,\s*([^,\)]+))?', src):
    opt, typ, sec, key, dflt = m.groups()
    keys.setdefault(key, {'sec': sec, 'type': typ.strip(), 'default': (dflt or '').strip(), 'required': not opt})
for m in re.finditer(r'(\w+)\["([A-Za-z0-9_]+)"\]', src):
    sec, key = m.groups()
    if key not in keys and sec in ('deltaT','cond','turb','physProp','space','last','out','ch','config','mesh','time'):
        keys.setdefault(key, {'sec': sec, 'type': '?', 'default': '', 'required': False})
# run config での使用
used=collections.Counter(); vals=collections.defaultdict(set)
cfgs=glob.glob(ROOT+'/case/*/run_*/solverConfig.yaml')
for p in cfgs:
    try: txt=open(p, errors='replace').read()
    except Exception: continue
    for k in keys:
        for m in re.finditer(r'(?<![A-Za-z0-9_])'+re.escape(k)+r'\s*:\s*([^,}\n]+)', txt):
            used[k]+=1; vals[k].add(m.group(1).strip().strip('"\''))
# 文書での言及
docs=' '.join(open(f, errors='replace').read() for f in glob.glob(ROOT+'/procedures/*.md')+glob.glob(ROOT+'/methods/**/*.md', recursive=True))
rec=open(ROOT+'/procedures/recommended-settings.md', errors='replace').read()
rows=[]
for k, d in sorted(keys.items()):
    nonz={v for v in vals[k] if v not in ('0','0.0','false','none','')}
    rows.append({'key': k, 'sec': d['sec'], 'default': d['default'], 'runs': used[k],
                 'values': sorted(vals[k])[:6], 'nondefault': len(nonz) > 0,
                 'in_docs': k in docs, 'in_recommended': k in rec})
json.dump({'nconfig': len(cfgs), 'rows': rows}, open('/tmp/claude-1000/-home-sano-work-forge/96eabcfb-a7fa-40da-aae8-ef553208f9cd/scratchpad/key_inventory.json','w'), ensure_ascii=False, indent=1)
never=[r for r in rows if r['runs']==0]
docless=[r for r in rows if not r['in_docs']]
print(f'keys {len(rows)} / run configs {len(cfgs)}')
print(f'どの run でも設定されていない: {len(never)}')
print(f'procedures・methods に一度も出てこない: {len(docless)}')
print(f'recommended-settings に無い: {len([r for r in rows if not r["in_recommended"]])}')
print('\n-- 削除候補 (run で未使用 かつ 文書に無い) --')
for r in never:
    if not r['in_docs']: print(f"  {r['sec']:10s} {r['key']:34s} default {r['default']}")
