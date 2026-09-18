#!/usr/bin/env python3
"""段階移行キーを既存 run config から外す (plan config-key-pruning §4.2 S2 の前提)。

S1 で任意化した 5 キーは「書いてあると警告、S2 で起動時エラー」である。**S2 に進む前に、再実行したい既存 run の
config からこれらを外しておく**必要がある。外さないと過去の run を回し直せなくなる。

取り除くもの:
  `mesh.meshFormat` / `physProp.isCompressible` / `physProp.ro` / `time.last.control` — solver が読まないか合法値が 1 つ。
  `time.deltaT.detectNaNInterval` — **値は捨てない**。トップレベルに無ければそちらへ**移送**する
  (この別名は無効キーではなく、書いた値が実際に効く)。トップレベルに既にあればそちらが優先なので別名だけ捨てる。

**壊さないための取り決め** (codex plan-5):
  - **節を特定してから編集する**。`time.last.control` を消すつもりで `time.deltaT.control` を消さない。
  - **anchor / alias / merge key (`&a` `*a` `<<`) を含む入力は触らない**。節どうしが同じ実体を共有していると、
    片方から消したつもりが両方から消える (照合も共有のせいで気づけない)。手動移送の対象として報告する。
  - **重複キーを含む入力は触らない**。PyYAML は後勝ち、solver の yaml-cpp は**先勝ち**なので、
    PyYAML での照合が solver の実効設定を表さない。
  - **`meshFormat` が `hdf5` 以外の入力は触らない**。現行 solver では動かない旧入力 (`gmsh` = 旧 smac 系) であり、
    キーを消すだけでは移送にならない。
  - 書き戻す前に、**編集後の YAML が「元 − 消したキー (+ 移送)」と厳密一致**することを確かめる。1 か所でも違えば書かない。
  - **`--dry-run` も編集と照合まで同じことをする**。書き込みだけをしない (移送できない入力を dry-run で見つけられるように)。
  - 既定の対象は**このツリーだけ**。他 worktree は別ブランチの作業物なので、明示したパスでしか触らない。

使い方:
  migrate_solver_config.py [--dry-run] [PATH ...]      # PATH 省略時はこのツリーの case/ 配下
"""
import argparse
import copy
import glob
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
DROP = [('mesh', 'meshFormat'), ('physProp', 'isCompressible'), ('physProp', 'ro'), ('time.last', 'control')]
MOVE = ('time.deltaT', 'detectNaNInterval')


def get(node, path):
    for seg in path.split('.'):
        if not isinstance(node, dict) or seg not in node:
            return None
        node = node[seg]
    return node


def unsafe_reason(text):
    """自動移送してはいけない入力なら理由を返す。安全なら None。"""
    import yaml
    if re.search(r'(^|[\s\[{,])[&*][A-Za-z0-9_-]+', text) or re.search(r'^\s*<<\s*:', text, re.M):
        return 'anchor/alias/merge key を含む (節どうしが同じ実体を共有しうる)'

    dup = []

    class DupCheck(yaml.SafeLoader):
        pass

    def mapping(loader, node, deep=False):
        seen = set()
        for k, _ in node.value:
            key = loader.construct_object(k, deep=True)
            if key in seen:
                dup.append(str(key))
            seen.add(key)
        return yaml.SafeLoader.construct_mapping(loader, node, deep)

    DupCheck.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)
    try:
        y = yaml.load(text, Loader=DupCheck) or {}
    except Exception as e:
        return f'YAML として読めない: {str(e)[:60]}'
    if dup:
        return f'重複キー {sorted(set(dup))[:3]} を含む (solver の yaml-cpp は先勝ち、PyYAML は後勝ちで照合できない)'
    mf = get(y, 'mesh.meshFormat')
    if mf is not None and str(mf).strip('"') != 'hdf5':
        return f'mesh.meshFormat が {mf!r} (現行 solver では動かない旧入力。キーを消すだけでは移送にならない)'
    return None


def plan_for(y):
    """戻り: (消すもの [(節, キー)], 移送する値 or None)。"""
    drops = []
    for sec, key in DROP:
        n = get(y, sec)
        if isinstance(n, dict) and key in n:
            drops.append((sec, key))
    move = None
    n = get(y, MOVE[0])
    if isinstance(n, dict) and MOVE[1] in n:
        move = None if MOVE[1] in y else n[MOVE[1]]   # トップレベルに**キーがあれば**そちらが優先 (値が null でも)
        drops.append(MOVE)
    return drops, move


def section_span(text, sec):
    """節 `sec` の本文の範囲と形を `(開始, 終了, フロー形か)` で返す。無ければ None。"""
    lines = text.split('\n')
    off = [0]
    for l in lines:
        off.append(off[-1] + len(l) + 1)
    segs = sec.split('.')
    lo, hi, depth = 0, len(lines), -1
    for si, seg in enumerate(segs):
        pat = re.compile(r'^([ \t]*)' + re.escape(seg) + r'[ \t]*:(.*)$')
        found = None
        for i in range(lo, hi):
            m = pat.match(lines[i])
            if not m:
                continue
            ind = len(m.group(1).expandtabs(4))
            if ind <= depth:
                continue
            found = (i, ind, m.group(2))
            break
        if found is None:
            return None
        i, ind, rest = found
        depth = ind
        if '{' in rest and si == len(segs) - 1:
            a = off[i] + lines[i].index('{') + 1
            d, k = 1, a
            while k < len(text):
                if text[k] == '{':
                    d += 1
                elif text[k] == '}':
                    d -= 1
                    if d == 0:
                        return (a, k, True)
                k += 1
            return None
        lo = i + 1
        end = hi
        for j in range(i + 1, hi):
            if lines[j].strip() == '' or lines[j].lstrip().startswith('#'):
                continue
            if len(lines[j][:len(lines[j]) - len(lines[j].lstrip())].expandtabs(4)) <= ind:
                end = j
                break
        hi = end
    return (off[lo], off[hi] if hi < len(lines) else len(text), False)


def edit_text(text, drops, move):
    """節を特定したうえで、その範囲内のキーだけを最小編集で消す。"""
    for sec, key in drops:
        span = section_span(text, sec)
        if span is None:
            continue
        a, b, is_flow = span
        body = text[a:b]
        if is_flow:     # フロー形: , 区切りの 1 項目 (コメントは入らない)
            pats = [re.compile(r'[ \t]*' + re.escape(key) + r'[ \t]*:[^,}]*,[ \t\n]*'),
                    re.compile(r',[ \t\n]*' + re.escape(key) + r'[ \t]*:[^,}]*')]
        else:           # ブロック形: 行末まで (**末尾コメントに , があっても行ごと**)
            pats = [re.compile(r'^[ \t]*' + re.escape(key) + r'[ \t]*:[^\n]*(?:\n|$)', re.M)]
        for pat in pats:
            new = pat.sub('', body, count=1)
            if new != body:
                body = new
                break
        text = text[:a] + body + text[b:]
    if move is not None:
        text = text.rstrip('\n') + f'\ndetectNaNInterval: {move}\n'
    return text


def migrate_one(text):
    """戻り: (新しい本文 or None, 状態文字列)。`None` なら書き換えない。"""
    import yaml
    why = unsafe_reason(text)
    if why:
        return None, 'skip: ' + why
    y = yaml.safe_load(text) or {}
    drops, move = plan_for(y)
    if not drops:
        return None, 'clean'
    out = edit_text(text, drops, move)
    want = copy.deepcopy(y)
    for sec, key in drops:
        n = get(want, sec)
        if isinstance(n, dict):
            n.pop(key, None)
    if move is not None:
        want['detectNaNInterval'] = move
    try:
        got = yaml.safe_load(out) or {}
    except Exception as e:
        return None, f'fail: 編集後に YAML として壊れた ({str(e)[:50]})'
    if got != want:
        d = sorted(k for k in set(list(got) + list(want)) if got.get(k) != want.get(k))
        return None, f'fail: 期待と一致しない (差分節 {d})'
    return out, ('moved' if move is not None else 'ok')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='*', help='対象 (既定: このツリーの case/ 配下。他 worktree は明示が要る)')
    ap.add_argument('--dry-run', action='store_true', help='編集と照合まで行い、書き込みだけしない')
    a = ap.parse_args()                       # 未知オプションはここで拒否される

    files = []
    for p in (a.paths or [ROOT + '/case']):
        files += (glob.glob(p + '/**/solverConfig.yaml', recursive=True) if os.path.isdir(p) else [p])
    files = sorted(set(files))

    n = {'ok': 0, 'moved': 0, 'clean': 0}
    skipped, failed = [], []
    for f in files:
        try:
            text = open(f, errors='replace').read()
        except Exception as e:
            failed.append((f, f'読めない: {str(e)[:50]}'))
            continue
        out, status = migrate_one(text)
        if status.startswith('skip:'):
            skipped.append((f, status[6:]))
            continue
        if status.startswith('fail:'):
            failed.append((f, status[6:]))
            continue
        n[status] += 1
        if out is not None and not a.dry_run:
            with open(f, 'w') as fh:
                fh.write(out)

    rel = lambda p: os.path.relpath(p, os.path.dirname(ROOT))
    print(f'対象 {len(files)} 本: 移送 {n["ok"] + n["moved"]} 本 '
          f'(うち detectNaNInterval をトップレベルへ {n["moved"]} 本) / 既に清潔 {n["clean"]} 本 / '
          f'手動対応 {len(skipped)} 本 / 失敗 {len(failed)} 本'
          + ('  [--dry-run: 書き込みなし]' if a.dry_run else ''))
    for f, why in skipped:
        print(f'  手動: {rel(f)} — {why}')
    for f, why in failed:
        print(f'  失敗: {rel(f)} — {why}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
