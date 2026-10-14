#!/usr/bin/env python3
"""solverConfig のキー棚卸し (plan config-key-pruning §5.1 #1)。

`solverConfig.{cpp,hpp}` が読むキーを**完全修飾パス**で抽出し、全 worktree の run config での使用実績と
`procedures/` `methods/` での言及を突き合わせて、分類表の素を出す。

数え方 (2026-09-18 第 2 次改訂; codex plan レビュー 2 の指摘 2/3/4):
  - 抽出は**コメントを除いたソース**に対して行う (コメント中のキー名を live と数えない)。
  - 変数がどの節を指すかは**ソースの束縛から導く** (`auto cj = config["conjugate"];` → `conjugate`)。
    ヘルパの引数は括弧の対応で分け、第 4 引数 (表示用の節名) は束縛が解けないときだけ使う (2026-10-14)。
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

# YAML ノードの添字の連鎖 `var["a"]["b"]` (添字は文字列リテラルだけ)。変数がどの節を指すかはソースの束縛から導く
# (2026-10-14 まで手書きの対応表 VAR2PATH で引いていたため、表に無い `cj` (conjugate)・`g` (conjugate.gate) を
#  取りこぼし、ソルバが受理する `conjugate.mode`・`conjugate.gate` を誤検出した。plan tooling-generator-defaults-and-launch-gate #4)。
_CHAIN = r'(\w+)((?:\s*\[\s*"[A-Za-z0-9_]+"\s*\])*)'
_LIT = re.compile(r'"([A-Za-z0-9_]+)"')
_NOT_MEMBER = r'(?<![\w.])(?<!->)'                     # `a.b["x"]`・`a->b["x"]` のメンバ参照は対象外
_DECL = r'(?:\bconst\s+)?(?:\bauto\b|\bYAML::Node\b)\s*&?\s*'
# 拒否の文言 (起動時に落とすための参照を、値を読む参照と分ける)
_REMOVAL_MSG = re.compile(r'no longer supported|撤去|廃止')


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


def _match_close(src, i):
    """src[i] が開き括弧 ( [ { のとき、対応する閉じ括弧の位置を返す (文字列・文字リテラルは飛ばす)。無ければ -1。"""
    pairs = {'(': ')', '[': ']', '{': '}'}
    stack, j, n = [], i, len(src)
    while j < n:
        c = src[j]
        if c in '"\'':
            k = j + 1
            while k < n and src[k] != c:
                k += 2 if src[k] == '\\' else 1
            j = k + 1
            continue
        if c in pairs:
            stack.append(pairs[c])
        elif c in ')]}':
            if not stack or stack.pop() != c:
                return -1
            if not stack:
                return j
        j += 1
    return -1


def _split_args(src, i):
    """src[i] == '(' の呼び出しの引数を最上位のカンマで分ける。戻り (引数のリスト, 閉じ括弧の位置)。

    既定値に括弧を含む呼び出し (`std::string("local1d")`、`(a == 1) ? 0.05 : 1.0`) を正規表現で読むと
    第 4 引数を落とす (2026-10-14 まで `conjugate.mode` が節名なしで抽出されていた原因の一つ)。"""
    end = _match_close(src, i)
    if end < 0:
        return None, -1
    args, depth, start, j = [], 0, i + 1, i + 1
    while j < end:
        c = src[j]
        if c in '"\'':
            k = j + 1
            while k < end and src[k] != c:
                k += 2 if src[k] == '\\' else 1
            j = k + 1
            continue
        if c in '([{':
            depth += 1
        elif c in ')]}':
            depth -= 1
        elif c == ',' and depth == 0:
            args.append(src[start:j].strip())
            start = j + 1
        j += 1
    args.append(src[start:end].strip())
    return args, end


def _guarded_body(src, pos):
    """pos が `if (...)` の条件の中にあれば、その if が守る文 (ブロックか 1 文) を返す。無ければ None。"""
    lo = max(0, pos - 400)
    for m in reversed(list(re.finditer(r'\bif\s*\(', src[lo:pos]))):
        op = lo + m.end() - 1
        cl = _match_close(src, op)
        if cl < 0 or cl < pos:
            continue
        k = cl + 1
        while k < len(src) and src[k].isspace():
            k += 1
        if k < len(src) and src[k] == '{':
            e = _match_close(src, k)
            return src[k:e + 1] if e > 0 else None
        e = src.find(';', k)
        return src[k:e + 1] if e > 0 else None
    return None


class _NodeScope:
    """YAML ノードを指す変数 (`auto cj = config["conjugate"];` など) の束縛をソースから集めて、位置ごとに解決する。

    - 根: `YAML::Node X = YAML::LoadFile(...)`。
    - 別名: `auto X = Y["a"]["b"];` / `const YAML::Node X = Y["a"];` (Y は根か別名)。
    - 範囲 for の要素: `for (const auto& X : Y["a"])` は列の要素 (値の中身) で、設定のパスではない。
    同名の変数は**その位置より前の最後の束縛**で解決する (C++ の有効範囲の近似。違えば解決不能として報告に出る)。"""

    def __init__(self, src):
        self.defs = collections.defaultdict(list)          # 変数 → [(位置, ('path'|'elem', パスの部品) | None)]
        events = []
        for m in re.finditer(_DECL + r'(\w+)\s*=\s*([^;]*);', src):
            events.append((m.start(), m.group(1), 'path', m.group(2).strip()))
        for m in re.finditer(r'\bfor\s*\(\s*' + _DECL + r'(\w+)\s*:\s*([^)]*)\)', src):
            events.append((m.start(), m.group(1), 'elem', m.group(2).strip()))
        for pos, var, kind, rhs in sorted(events):
            val = None
            if kind == 'path' and re.match(r'YAML::Load(?:File)?\s*\(', rhs):
                val = ('path', ())
            else:
                mm = re.fullmatch(_CHAIN, rhs)
                if mm:
                    base = self.resolve(mm.group(1), pos)
                    if base is not None and base[0] == 'path':
                        val = (kind, base[1] + tuple(_LIT.findall(mm.group(2))))
                    elif base is not None:
                        val = base                          # 列の要素の中身の別名も、列の要素のまま
            self.defs[var].append((pos, val))

    def resolve(self, var, pos):
        """戻り ('path', 部品) / ('elem', 列のパスの部品) / None (YAML ノードと分からない・解決できない)。"""
        cand = [v for p, v in self.defs.get(var, []) if p < pos]
        return cand[-1] if cand else None

    def is_node_var(self, var, pos):
        """var が pos の時点で YAML ノードの束縛を持つ変数として宣言されたことがあるか (解決できたかは問わない)。"""
        return any(p < pos for p, _ in self.defs.get(var, []))


def extract_keys(src):
    """戻り {path: {...}}。path は 'time.deltaT.dt' のような完全修飾。

    kind: 'value' (値を読む) / 'section' (節そのもの) / 'rejected' (拒否専用の参照)。
    default: 既定値の文字列。required なら None、解決できなければ '?'。
    """
    return extract_keys_with_report(src)[0]


def extract_keys_with_report(src):
    """extract_keys と同じ表に加え、静的に解決できなかった参照の報告を返す。戻り (keys, report)。

    report は [{'line', 'why', 'text'}]。why の種類:
      - 'unresolved-node': YAML ノードの変数だが、どの節を指すかをソースから導けない (ヘルパの表示用の節名で代用した)。
      - 'list-element': 列の要素への添字 (値キーの中身。設定のパスではないので表に入れない)。
      - 'dynamic-subscript': 文字列リテラルでない添字 (ループ形の拒否表・削除表の書き方以外)。
      - 'conditional-reject': 値を読むキーを、別の設定値の条件付きで拒否している (静的な一覧では表せない)。
      - 'label-mismatch': ヘルパの表示用の節名と、束縛から導いた節が違う (表示の誤り。抽出は束縛を使う)。
    """
    keys, report = {}, []
    refs = collections.defaultdict(set)                     # path → {'read', 'reject', 'reject-loop'}
    reject_line = {}                                        # path → 拒否している行 (報告用)
    scope = _NodeScope(src)

    def line(pos):
        return src.count('\n', 0, pos) + 1

    def add(path, typ, dflt, required, how, member=None):
        e = keys.setdefault(path, {'type': typ, 'default': dflt, 'required': required,
                                   'how': set(), 'member': member, 'kind': 'value'})
        if dflt is not None and e['default'] in (None, '?', ''):
            e['default'] = dflt
        if member and not e['member']:
            e['member'] = member
        e['required'] = e['required'] or required
        e['how'].add(how)
        # 読むパスの親は必ず写像 (節)。親だけを読む参照が無くても節として載せる
        parts = path.split('.')
        for i in range(1, len(parts)):
            keys.setdefault('.'.join(parts[:i]), {'type': '-', 'default': None, 'required': False,
                                                  'how': set(), 'member': None, 'kind': 'value'})['how'].add('prefix')
        return e

    # 1) ヘルパ: getValidatedValue<T>(node, "key"[, "label"]) / getOptionalValidatedValue<T>(node, "key", default[, "label"])
    #    引数は括弧の対応で分ける。節は第 1 引数の束縛から導き、表示用の節名 (label) は束縛が解けないときだけ使う。
    for m in re.finditer(r'\b(getValidatedValue|getOptionalValidatedValue)\s*<', src):
        i, depth = m.end(), 1
        while i < len(src) and depth:
            depth += {'<': 1, '>': -1}.get(src[i], 0)
            i += 1
        typ = src[m.end():i - 1].strip()
        k = i
        while k < len(src) and src[k].isspace():
            k += 1
        if k >= len(src) or src[k] != '(':
            continue                                        # 定義側 (template の宣言) など
        args, _ = _split_args(src, k)
        required = m.group(1) == 'getValidatedValue'
        if not args or len(args) < (2 if required else 3):
            report.append({'line': line(m.start()), 'why': 'unresolved-node', 'text': src[m.start():k + 80].split('\n')[0]})
            continue
        km = re.fullmatch(r'"([A-Za-z0-9_]+)"', args[1])
        label_arg = args[2] if required and len(args) > 2 else (args[3] if not required and len(args) > 3 else None)
        lm = re.fullmatch(r'"([^"]*)"', label_arg) if label_arg is not None else None
        label = lm.group(1).strip() if lm else None
        nm = re.fullmatch(_CHAIN, args[0])
        val = scope.resolve(nm.group(1), m.start()) if nm else None
        if not km:
            report.append({'line': line(m.start()), 'why': 'dynamic-subscript', 'text': f'{m.group(1)}({args[0]}, {args[1]}, ...)'})
            continue
        if val is not None and val[0] == 'path':
            sec = '.'.join(val[1] + tuple(_LIT.findall(nm.group(2))))
            if label is not None and label != sec:
                report.append({'line': line(m.start()), 'why': 'label-mismatch',
                               'text': f'{args[0]}, "{km.group(1)}": 表示 "{label}" / 束縛 "{sec}"'})
        else:
            report.append({'line': line(m.start()), 'why': 'unresolved-node',
                           'text': f'{m.group(1)}({args[0]}, "{km.group(1)}") — 表示用の節名 "{label}" で代用'})
            if label is None:
                continue
            sec = label
        path = '.'.join([p for p in [sec, km.group(1)] if p])
        dflt = None
        if not required:
            dflt = args[2].strip()
            sm = re.fullmatch(r'std::string\s*\(\s*("[^"]*")\s*\)', dflt)
            if sm:
                dflt = sm.group(1)
        add(path, typ, dflt, required, 'helper-required' if required else 'helper-optional')
        refs[path].add('read')

    # 2) 直接参照 X["a"]["b"] (X は根・別名)。if の条件の中で、守る文が拒否の文言で落とすなら拒否専用。
    for m in re.finditer(_NOT_MEMBER + r'(\w+)((?:\s*\[\s*"[A-Za-z0-9_]+"\s*\])+)', src):
        var, chain = m.groups()
        val = scope.resolve(var, m.start())
        if val is None:
            if scope.is_node_var(var, m.start()):
                report.append({'line': line(m.start()), 'why': 'unresolved-node', 'text': m.group(0)})
            continue
        if val[0] == 'elem':
            report.append({'line': line(m.start()), 'why': 'list-element',
                           'text': f'{m.group(0)} (列 {".".join(val[1])} の要素)'})
            continue
        path = '.'.join(val[1] + tuple(_LIT.findall(chain)))
        if re.match(r'\s*\[', src[m.end():]):
            add(path, '?', None, False, 'prefix')          # `config["mesh"][k]` の親: 節であって、この参照は値も拒否も表さない
            continue
        # 同じ文の代入先メンバ (this->meshRenumber = config["mesh"]["renumber"]...) を既定値の引き当てに使う
        line_start = src.rfind('\n', 0, m.start()) + 1
        stmt_end = src.find(';', m.end())
        stmt = src[line_start:stmt_end if stmt_end > 0 else m.end()]
        mm = re.search(r'this->(\w+)\s*=', stmt)
        add(path, '?', None, False, 'direct', mm.group(1) if mm else None)
        body = _guarded_body(src, m.start())
        if body is not None and re.search(r'\bthrow\b|\bexit\s*\(', body) and _REMOVAL_MSG.search(body):
            refs[path].add('reject')
            reject_line.setdefault(path, line(m.start()))
        else:
            refs[path].add('read')

    # 3) ループ形の拒否表: for (const char* k : {"A","B"}) [{] if (X[k]) ...throw/exit...
    loop_spans = []
    for m in re.finditer(r'\bfor\s*\(\s*const\s+char\s*\*\s*(\w+)\s*:\s*\{([^}]*)\}\s*\)\s*\{?\s*'
                         r'if\s*\(\s*' + _CHAIN + r'\s*\[\s*\1\s*\]', src):
        _, lst, var, chain = m.groups()
        val = scope.resolve(var, m.start())
        if val is None or val[0] != 'path':
            report.append({'line': line(m.start()), 'why': 'unresolved-node', 'text': m.group(0)[:80]})
            continue
        body = _guarded_body(src, m.end() - 1)
        if body is None or not re.search(r'\bthrow\b|\bexit\s*\(', body):
            continue
        loop_spans.append((m.start(), m.end()))
        sec = val[1] + tuple(_LIT.findall(chain))
        for key in _LIT.findall(lst):
            path = '.'.join(sec + (key,))
            add(path, '-', None, False, 'rejected-loop')
            refs[path].add('reject-loop')
            reject_line.setdefault(path, line(m.start()))

    # 4) 文字列リテラルでない添字 (ループ形の拒否表・削除表 removed[] の走査は上と removed_keys() が読む)
    rm = re.search(r'\bfor\s*\(\s*const\s+auto\s*&\s*\w+\s*:\s*removed\s*\)', src)
    if rm:
        e = src.find('{', rm.end())
        loop_spans.append((rm.start(), _match_close(src, e) if e >= 0 else rm.end()))
    for m in re.finditer(_NOT_MEMBER + r'(\w+)((?:\s*\[\s*"[A-Za-z0-9_]+"\s*\])*)\s*\[\s*(?!["\s])([^\]]+)\]', src):
        val = scope.resolve(m.group(1), m.start())
        if val is None or val[0] != 'path' or any(a <= m.start() <= b for a, b in loop_spans):
            continue
        report.append({'line': line(m.start()), 'why': 'dynamic-subscript', 'text': m.group(0)})

    # 拒否専用の判定: どこかで値として読むパスは拒否にしない (条件付きの拒否は報告に回す)
    for path, r in refs.items():
        if r and r <= {'reject', 'reject-loop'}:
            keys[path]['kind'] = 'rejected'
        elif 'reject-loop' in r or 'reject' in r:
            report.append({'line': reject_line.get(path, 0), 'why': 'conditional-reject', 'text': path})
    return keys, report


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
    _, unresolved = extract_keys_with_report(src)

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
    if unresolved:
        cnt = collections.Counter(r['why'] for r in unresolved)
        print('静的に一覧へ入れなかった参照: ' + ', '.join(f'{k} {v}' for k, v in sorted(cnt.items()))
              + ' (list-element は値キーの中身、conditional-reject は別の設定値に依る拒否)')
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
