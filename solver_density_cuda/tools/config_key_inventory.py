#!/usr/bin/env python3
"""solverConfig のキー棚卸し (plan config-key-pruning §5.1 #1)。

`solverConfig.{cpp,hpp}` が読むキーを**完全修飾パス**で抽出し、全 worktree の run config での使用実績と
`procedures/` `methods/` での言及を突き合わせて、分類表の素を出す。

数え方 (2026-09-18 改訂; 旧実装はキー名ベースで節をまたいで合算・重複計上し「使用 0」を誤って出していた):
  - run config は **PyYAML で読み、完全修飾パスに平坦化**して引く。
  - run は **`<case>/<run ディレクトリ>` で重複排除**する (同じ run が複数 worktree に複製されているため)。
  - 「非既定」は**そのパスの既定値と数値比較**して判定する (既定値を明記しただけの run は非既定に数えない)。
  - 起動時拒否テーブル (`removed[]`) にあるキーは live ではないので `removed` 印をつける。

使い方: config_key_inventory.py [--json OUT.json] [--csv OUT.csv]
"""
import collections
import glob
import json
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
SRC_CPP = os.path.join(ROOT, 'solver_density_cuda', 'input', 'solverConfig.cpp')
JSON_DEFAULT = '/tmp/claude-1000/-home-sano-work-forge/96eabcfb-a7fa-40da-aae8-ef553208f9cd/scratchpad/key_inventory.json'

# ヘルパ呼び出しの第 1 引数 (変数名) → YAML 上の節パス
VAR2PATH = {'deltaT': 'time.deltaT', 'last': 'time.last', 'turb': 'turbulence', 'cond': 'condensation',
            'ch': 'physProp.chemistry', 'physProp': 'physProp', 'space': 'space', 'out': 'output',
            'config': '', 'time': 'time', 'mesh': 'mesh'}
TOP_SECTIONS = ('mesh', 'time', 'space', 'turbulence', 'condensation', 'physProp', 'output')


def extract_keys(src):
    """戻り {path: {'default','type','required','how'}}。path は 'time.deltaT.dt' のような完全修飾。"""
    keys = {}

    def add(path, typ, dflt, required, how):
        e = keys.setdefault(path, {'type': typ, 'default': (dflt or '').strip(), 'required': required, 'how': set()})
        if dflt and not e['default']:
            e['default'] = dflt.strip()
        e['how'].add(how)
        return e

    # 1) get[Optional]ValidatedValue<T>(sec, "key", default, "section.path")
    #    第 1 引数は変数名でも config["time"] 形でもありうる。
    pat = (r'get(Optional)?ValidatedValue<([^>]+)>\(\s*(\w+(?:\[\s*"[A-Za-z0-9_]+"\s*\])*)\s*,\s*'
           r'"([A-Za-z0-9_]+)"\s*(?:,\s*([^,\)]+))?\s*(?:,\s*"([^"]*)")?')
    for m in re.finditer(pat, src):
        opt, typ, var, key, dflt, secstr = m.groups()
        if '[' in var:
            inner = re.findall(r'"([A-Za-z0-9_]+)"', var)
            base = '' if var.startswith('config') else VAR2PATH.get(var.split('[')[0], var.split('[')[0])
            sec = '.'.join([p for p in [base] + inner if p])
        else:
            sec = (secstr or VAR2PATH.get(var, var)).strip()
        add('.'.join([p for p in [sec, key] if p]), typ.strip(), dflt, not opt, 'helper')

    # 2) 直接参照 config["a"]["b"] / sec["key"]
    for m in re.finditer(r'(config|[a-zA-Z_]\w*)((?:\[\s*"[A-Za-z0-9_]+"\s*\])+)', src):
        var, chain = m.groups()
        parts = re.findall(r'"([A-Za-z0-9_]+)"', chain)
        if var == 'config':
            path = '.'.join(parts)
        else:
            base = VAR2PATH.get(var)
            if base is None:
                continue
            path = '.'.join([p for p in [base] + parts if p])
        if not path:
            continue
        if len(path.split('.')) > 1 and path.split('.')[0] not in TOP_SECTIONS:
            continue
        add(path, '?', '', False, 'direct')
    return keys


def removed_keys(src):
    """起動時拒否テーブル removed[] に載っているキー (= live でない) を {path} で返す。"""
    out = set()
    m = re.search(r'static const Removed removed\[\]\s*=\s*\{(.*?)\}\s*;', src, re.S)
    if not m:
        return out
    for e in re.finditer(r'\{\s*"([A-Za-z0-9_]+)"\s*,\s*(?:"([A-Za-z0-9_]+)"|nullptr)\s*,\s*"([A-Za-z0-9_]+)"', m.group(1)):
        s1, s2, key = e.groups()
        out.add('.'.join([p for p in [s1, s2, key] if p]))
    return out


def flatten(node, prefix=()):
    """YAML を {完全修飾パス: 値} に平坦化する (リストは葉として扱う)。"""
    out = {}
    if isinstance(node, dict):
        for k, v in node.items():
            out.update(flatten(v, prefix + (str(k),)))
    elif prefix:
        out['.'.join(prefix)] = node
    return out


def same_as_default(val, dflt):
    d = (dflt or '').strip().rstrip('f').strip('"')
    if d == '':
        return False
    try:
        return abs(float(val) - float(d)) <= 1e-12*max(1.0, abs(float(d)))
    except (TypeError, ValueError):
        return str(val).strip('"') == d


def main():
    import yaml
    src = open(SRC_CPP, errors='replace').read() + open(SRC_CPP.replace('.cpp', '.hpp'), errors='replace').read()
    keys = extract_keys(src)
    gone = removed_keys(src)
    # 直読み (config["mesh"][...] 形) は既定値が式に出てこないので、hpp のメンバ初期化子から補う
    hpp = open(SRC_CPP.replace('.cpp', '.hpp'), errors='replace').read()
    member_default = {m.group(1): m.group(2).strip()
                      for m in re.finditer(r'^\s*(?:int|flow_float|double|std::string)\s+(\w+)\s*=\s*([^;]+);', hpp, re.M)}
    for path, d in keys.items():
        if not d['default']:
            d['default'] = member_default.get(path.split('.')[-1], '')
    # 拒否済みキーは読みが消えているので抽出に出てこない。表には「removed」として並べる
    for path in gone:
        keys.setdefault(path, {'type': '-', 'default': '(removed)', 'required': False, 'how': {'removed'}})

    trees = [ROOT] + [d for d in sorted(glob.glob(os.path.dirname(ROOT) + '/forge-*')) if os.path.isdir(d + '/case')]
    seen, used, vals, nrun = set(), collections.Counter(), collections.defaultdict(set), 0
    for t in trees:
        for f in glob.glob(t + '/case/*/run_*/solverConfig.yaml'):
            rid = '/'.join(f.split('/')[-3:-1])          # <case>/<run> で重複排除
            if rid in seen:
                continue
            seen.add(rid)
            nrun += 1
            try:
                y = yaml.safe_load(open(f, errors='replace')) or {}
            except Exception:
                continue
            for path, v in flatten(y).items():
                if path in keys and not same_as_default(v, keys[path]['default']):
                    used[path] += 1
                    vals[path].add(str(v))

    docs = ' '.join(open(x, errors='replace').read()
                    for x in glob.glob(ROOT + '/procedures/*.md') + glob.glob(ROOT + '/methods/**/*.md', recursive=True))
    rec = open(ROOT + '/procedures/recommended-settings.md', errors='replace').read()

    rows = []
    for path, d in sorted(keys.items()):
        key = path.split('.')[-1]
        rows.append({'path': path, 'key': key, 'sec': '.'.join(path.split('.')[:-1]) or '(top)',
                     'default': d['default'], 'required': d['required'], 'how': sorted(d['how']),
                     'removed': path in gone, 'nondefault_runs': used[path],
                     'values': sorted(vals[path])[:8], 'in_docs': key in docs, 'in_recommended': key in rec})

    out = {'nrun': nrun, 'ntrees': len(trees), 'rows': rows}
    jpath = JSON_DEFAULT
    if '--json' in sys.argv:
        jpath = sys.argv[sys.argv.index('--json') + 1]
    try:
        json.dump(out, open(jpath, 'w'), ensure_ascii=False, indent=1)
    except OSError:
        jpath = None
    live = [r for r in rows if not r['removed']]
    print(f'live キー (完全修飾パス) {len(live)} / 拒否済み {len(rows)-len(live)} / run config {nrun} 本 ({len(trees)} worktree)')
    print(f"どの run でも既定以外の値で設定されていない: {len([r for r in live if r['nondefault_runs'] == 0])}")
    print(f"procedures・methods に一度も出てこない: {len([r for r in live if not r['in_docs']])}")
    print(f"recommended-settings に無い: {len([r for r in live if not r['in_recommended']])}")
    if jpath:
        print(f'詳細: {jpath}')
    print('\n-- 非既定の使用が無く、文書にも出てこない live キー (削除候補の一次ふるい) --')
    for r in live:
        if r['nondefault_runs'] == 0 and not r['in_docs']:
            print(f"  {r['path']:46s} default {r['default']}")


if __name__ == '__main__':
    main()
