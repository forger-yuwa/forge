#!/usr/bin/env python3
"""solverConfig のキー棚卸し (plan config-key-pruning §5.1 #1)。

`solverConfig.{cpp,hpp}` が読むキーを**完全修飾パス**で抽出し、全 worktree の run config での使用実績と
`procedures/` `methods/` での言及を突き合わせて、分類表の素を出す。

数え方 (2026-09-18 第 2 次改訂; codex plan レビュー 2 の指摘 2/3/4):
  - 抽出は**コメントを除いたソース**に対して行う (コメント中のキー名を live と数えない)。
  - **節そのもの** (`time.deltaT` など) と**値を読むキー**を分ける。節は `--all` でしか出さない。
  - **拒否専用の参照** (起動時に落とすためだけに見ているキー) は `rejected` 印をつけ live から外す。
  - **必須キーの既定値は「なし」**。ヘルパ第 4 引数はエラー表示用の節名であって既定値ではない。
  - run config は `case/` 配下を**入れ子まで**探索し、PyYAML で読んで完全修飾パスに平坦化して引く。
    識別子は `<case>/<case 以下の相対パス>` で、**同じ識別子で内容が違うもの**は別設定として数える。
    読み取りに失敗したファイルは黙って捨てず、集計不完全として報告する。
  - 「記載 run 数」と「非既定 run 数」を**別々に**出す。既定値が解決できないキーは非既定を判定しない (不明)。
  - 既定との比較は**数値表記を正規化した厳密比較**で行う (一律の絶対許容差を置かない。1e-20 と 1e-30 を
    同じと見なしてしまうため)。

使い方: config_key_inventory.py [--json OUT.json] [--all]
"""
import collections
import glob
import hashlib
import json
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
SRC_CPP = os.path.join(ROOT, 'solver_density_cuda', 'input', 'solverConfig.cpp')
SRC_HPP = SRC_CPP.replace('.cpp', '.hpp')
JSON_DEFAULT = os.path.join(os.environ.get('TMPDIR', '/tmp'), 'key_inventory.json')

# ヘルパ呼び出しの第 1 引数 (変数名) → YAML 上の節パス
VAR2PATH = {'deltaT': 'time.deltaT', 'last': 'time.last', 'turb': 'turbulence', 'cond': 'condensation',
            'ch': 'physProp.chemistry', 'physProp': 'physProp', 'space': 'space', 'out': 'output',
            'config': '', 'time': 'time', 'mesh': 'mesh'}
TOP_SECTIONS = ('mesh', 'time', 'space', 'turbulence', 'condensation', 'physProp', 'output')


def strip_comments(src):
    """C/C++ のコメントを空白に潰す (文字列リテラルの中は触らない)。行数は保つ。"""
    out, i, n = [], 0, len(src)
    while i < n:
        c = src[i]
        if c == '"':
            j = i + 1
            while j < n and src[j] != '"':
                j += 2 if src[j] == '\\' else 1
            out.append(src[i:j + 1]); i = j + 1
        elif src.startswith('//', i):
            j = src.find('\n', i)
            j = n if j < 0 else j
            out.append(' ' * (j - i)); i = j
        elif src.startswith('/*', i):
            j = src.find('*/', i + 2)
            j = n - 2 if j < 0 else j
            out.append(''.join(ch if ch == '\n' else ' ' for ch in src[i:j + 2])); i = j + 2
        else:
            out.append(c); i += 1
    return ''.join(out)


def extract_keys(src):
    """戻り {path: {...}}。path は 'time.deltaT.dt' のような完全修飾。

    kind: 'value' (値を読む) / 'section' (節そのもの) / 'rejected' (拒否専用の参照)。
    default: 既定値の文字列。required なら None、解決できなければ '?'。
    """
    keys = {}

    def add(path, typ, dflt, required, how, member=None):
        e = keys.setdefault(path, {'type': typ, 'default': dflt, 'required': required,
                                   'how': set(), 'member': member, 'kind': 'value'})
        if dflt is not None and e['default'] in (None, '?', ''):
            e['default'] = dflt
        if member and not e['member']:
            e['member'] = member
        e['required'] = e['required'] or required
        e['how'].add(how)
        return e

    def sec_of(var, secstr):
        if '[' in var:
            inner = re.findall(r'"([A-Za-z0-9_]+)"', var)
            base = '' if var.startswith('config') else VAR2PATH.get(var.split('[')[0], var.split('[')[0])
            return '.'.join([p for p in [base] + inner if p])
        return (secstr or VAR2PATH.get(var, var)).strip()

    # 1) 必須ヘルパ: getValidatedValue<T>(sec, "key", "section.path")  — 第 3 引数はエラー表示用の節名
    pat_req = (r'getValidatedValue<([^>]+)>\(\s*(\w+(?:\[\s*"[A-Za-z0-9_]+"\s*\])*)\s*,\s*'
               r'"([A-Za-z0-9_]+)"\s*(?:,\s*"([^"]*)")?\s*\)')
    for m in re.finditer(pat_req, src):
        typ, var, key, secstr = m.groups()
        sec = sec_of(var, secstr)
        add('.'.join([p for p in [sec, key] if p]), typ.strip(), None, True, 'helper-required')

    # 2) 任意ヘルパ: getOptionalValidatedValue<T>(sec, "key", default, "section.path")
    pat_opt = (r'getOptionalValidatedValue<([^>]+)>\(\s*(\w+(?:\[\s*"[A-Za-z0-9_]+"\s*\])*)\s*,\s*'
               r'"([A-Za-z0-9_]+)"\s*,\s*([^,\)]+?)\s*(?:,\s*"([^"]*)")?\s*\)')
    for m in re.finditer(pat_opt, src):
        typ, var, key, dflt, secstr = m.groups()
        sec = sec_of(var, secstr)
        add('.'.join([p for p in [sec, key] if p]), typ.strip(), dflt.strip(), False, 'helper-optional')

    # 3) 直接参照 config["a"]["b"] / sec["key"]
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
        if not path or (len(path.split('.')) > 1 and path.split('.')[0] not in TOP_SECTIONS):
            continue
        # 同じ文の代入先メンバ (this->meshRenumber = config["mesh"]["renumber"]...) を既定値の引き当てに使う
        line_start = src.rfind('\n', 0, m.start()) + 1
        stmt_end = src.find(';', m.end())
        stmt = src[line_start:stmt_end if stmt_end > 0 else m.end()]
        mm = re.search(r'this->(\w+)\s*=', stmt)
        e = add(path, '?', None, False, 'direct', mm.group(1) if mm else None)
        # 拒否専用の参照: 同じ文の直後に "no longer supported" の throw がある
        tail = src[m.end():m.end() + 400]
        if '.IsDefined()' in stmt and 'no longer supported' in tail:
            e['kind'] = 'rejected'

    # 4) ループ形の拒否表: for (const char* k : {"A","B"}) { if (sec[k]) { ...exit/throw... } }
    for m in re.finditer(r'for\s*\(\s*const char\*\s*(\w+)\s*:\s*\{([^}]*)\}\s*\)\s*\{\s*'
                         r'if\s*\(\s*(\w+(?:\[\s*"[A-Za-z0-9_]+"\s*\])*)\s*\[\s*\1\s*\]', src):
        _, lst, var = m.groups()
        sec = sec_of(var, None)
        for key in re.findall(r'"([A-Za-z0-9_]+)"', lst):
            add('.'.join([p for p in [sec, key] if p]), '-', None, False, 'rejected-loop')['kind'] = 'rejected'
    return keys


def removed_keys(src):
    """起動時拒否テーブル removed[] に載っているキー (= live でない) を {path: note} で返す。"""
    out = {}
    m = re.search(r'static const Removed removed\[\]\s*=\s*\{(.*?)\}\s*;', src, re.S)
    if not m:
        return out
    for e in re.finditer(r'\{\s*"([A-Za-z0-9_]+)"\s*,\s*(?:"([A-Za-z0-9_]+)"|nullptr)\s*,\s*"([A-Za-z0-9_]+)"\s*,\s*"([^"]*)"', m.group(1)):
        s1, s2, key, note = e.groups()
        out['.'.join([p for p in [s1, s2, key] if p])] = note
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


def norm_num(x):
    """数値表記を正規化して float にする。数値でなければ None。"""
    if isinstance(x, bool):
        return float(int(x))
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().strip('"').rstrip('fF')
    s = re.sub(r'^\(.*?\)\s*', '', s)          # (flow_float)0.0 のようなキャスト
    try:
        return float(s)
    except ValueError:
        return None


def same_as_default(val, dflt):
    """既定と同じなら True、違えば False、判定できなければ None。"""
    if dflt is None or dflt in ('', '?'):
        return None
    a, b = norm_num(val), norm_num(dflt)
    if a is not None and b is not None:
        return a == b                           # 正規化した厳密比較 (絶対許容差を置かない)
    return str(val).strip().strip('"') == str(dflt).strip().strip('"')


def main():
    import yaml
    raw = open(SRC_CPP, errors='replace').read() + '\n' + open(SRC_HPP, errors='replace').read()
    src = strip_comments(raw)
    keys = extract_keys(src)
    gone = removed_keys(src)

    # 節そのもの (他のパスの真の接頭辞) は値キーでない
    allpaths = set(keys)
    for path in allpaths:
        if any(p != path and p.startswith(path + '.') for p in allpaths):
            keys[path]['kind'] = 'section'

    # 直読みキーの既定値は hpp のメンバ初期化子から引く (末端名でなく**代入先メンバ名**で引き当てる)
    hpp = strip_comments(open(SRC_HPP, errors='replace').read())
    member_default = {m.group(1): m.group(2).strip()
                      for m in re.finditer(r'^\s*(?:int|flow_float|double|bool|std::string)\s+(\w+)\s*=\s*([^;]+);', hpp, re.M)}
    for path, d in keys.items():
        if d['default'] is None and not d['required'] and d['kind'] == 'value':
            d['default'] = member_default.get(d['member'], '?') if d['member'] else '?'
    for path, note in gone.items():
        keys.setdefault(path, {'type': '-', 'default': None, 'required': False,
                               'how': {'removed'}, 'member': None, 'kind': 'value'})['note'] = note

    trees = [ROOT] + [d for d in sorted(glob.glob(os.path.dirname(ROOT) + '/forge-*')) if os.path.isdir(d + '/case')]
    seen, used, written, vals, nconf, nfail = {}, collections.Counter(), collections.Counter(), collections.defaultdict(set), 0, []
    dup_conflict = 0
    for t in trees:
        for f in glob.glob(t + '/case/**/solverConfig.yaml', recursive=True):
            rel = os.path.relpath(f, t + '/case')          # <case>/<case 以下の相対パス>
            try:
                text = open(f, errors='replace').read()
                y = yaml.safe_load(text) or {}
            except Exception as e:
                nfail.append((os.path.relpath(f, os.path.dirname(ROOT)), str(e)[:80]))
                continue
            h = hashlib.md5(text.encode('utf-8', 'replace')).hexdigest()
            if rel in seen:
                if seen[rel] == h:
                    continue                                # 同一内容の worktree 複製は 1 本と数える
                dup_conflict += 1
                rel = rel + '#' + h[:8]                     # 内容が違うものは別設定として数える
                if rel in seen:
                    continue
            seen[rel] = h
            nconf += 1
            for path, v in flatten(y).items():
                if path not in keys:
                    continue
                written[path] += 1
                if same_as_default(v, keys[path]['default']) is False:
                    used[path] += 1
                    vals[path].add(str(v))

    docs = ' '.join(open(x, errors='replace').read()
                    for x in glob.glob(ROOT + '/procedures/*.md') + glob.glob(ROOT + '/methods/**/*.md', recursive=True))
    rec = open(ROOT + '/procedures/recommended-settings.md', errors='replace').read()

    rows = []
    for path, d in sorted(keys.items()):
        key = path.split('.')[-1]
        rows.append({'path': path, 'key': key, 'sec': '.'.join(path.split('.')[:-1]) or '(top)',
                     'kind': 'removed' if path in gone else d['kind'],
                     'default': ('(required)' if d['required'] else d['default']),
                     'required': d['required'], 'how': sorted(d['how']), 'member': d['member'],
                     'removed': path in gone, 'written_runs': written[path],
                     'nondefault_runs': (used[path] if ((not d['required']) and d['default'] not in (None, '?', '')) else None),
                     'default_known': (not d['required']) and d['default'] not in (None, '?', ''),
                     'values': sorted(vals[path])[:8], 'in_docs': key in docs, 'in_recommended': key in rec})

    out = {'nconf': nconf, 'ntrees': len(trees), 'nread_fail': len(nfail), 'ndup_conflict': dup_conflict, 'rows': rows}
    jpath = JSON_DEFAULT
    if '--json' in sys.argv:
        jpath = sys.argv[sys.argv.index('--json') + 1]
    try:
        json.dump(out, open(jpath, 'w'), ensure_ascii=False, indent=1)
    except OSError:
        jpath = None
    live = [r for r in rows if r['kind'] == 'value']
    print(f'live な値キー (完全修飾パス) {len(live)} / 節 {len([r for r in rows if r["kind"]=="section"])} / '
          f'拒否専用 {len([r for r in rows if r["kind"]=="rejected"])} / 起動時拒否 {len([r for r in rows if r["kind"]=="removed"])}')
    print(f'run config {nconf} 本 ({len(trees)} worktree, 内容違いの同名 {dup_conflict} 本を別設定として計上)')
    if nfail:
        print(f'**読み取り失敗 {len(nfail)} 本 (集計は不完全)**: ' + ', '.join(p for p, _ in nfail[:5]))
    known = [r for r in live if r['default_known']]
    print(f"必須キー (既定値なし): {len([r for r in live if r['required']])}")
    print(f"既定値が解決できず非既定を判定できない: {len([r for r in live if not r['default_known'] and not r['required']])}")
    print(f"どの run にも書かれていない (既定値が分かるもの): {len([r for r in known if r['written_runs'] == 0])}")
    print(f"書かれているが既定値のみ (既定値が分かるもの): {len([r for r in known if r['written_runs'] and r['nondefault_runs'] == 0])}")
    print(f"procedures・methods に一度も出てこない: {len([r for r in live if not r['in_docs']])}")
    if jpath:
        print(f'詳細: {jpath}')
    print('\n-- 非既定の使用が無く、文書にも出てこない live キー (削除候補の一次ふるい) --')
    for r in live:
        if r['default_known'] and r['nondefault_runs'] == 0 and not r['in_docs']:
            print(f"  {r['path']:46s} default {r['default']}  (記載 {r['written_runs']} run)")
    if '--all' in sys.argv:
        print('\n-- 全パス --')
        for r in rows:
            nd = '  不明' if r['nondefault_runs'] is None else f"{r['nondefault_runs']:6d}"
            print(f"  {r['kind']:8s} {r['path']:46s} default {str(r['default']):20s} 記載 {r['written_runs']:5d} 非既定 {nd}")


if __name__ == '__main__':
    main()
