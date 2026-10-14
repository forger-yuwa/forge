#!/usr/bin/env python3
"""機種と段ごとの推奨設定 (プロファイル) の読み込み・照合・宣言・文書の生成
(plan tooling-generator-defaults-and-launch-gate §4.1・§4.3・§5.1 #6・#16)。

値の正本は `recipes/recommended_settings.yaml`。プロファイル = 機種と段ごとの推奨設定で、各キーは次のどれか:

- `fixed`    宣言したプロファイルの値。違えば FAIL (理由つきの逸脱 `--deviate key=value --reason ...` を宣言すれば通す)。
- `explicit` 明示すべきキー (plan §5.1 #16): 値に依存するのに書き写しの出力 (生成器の config) が省略していたキーと、その旧い実効値
             (書き写しに使ったソルバの既定)。照合は `fixed` と同じ。#11 で生成器に書かせるときの「足してよいキーと旧い実効値」の
             例外表。基準は YAML の `explicit_keys` (`recipes/regen_profiles.py` が段ごとに作る)。
- `advisory` 推奨の範囲 (`[lo, hi]`・`{in: [...]}`・`any`)。外れたら WARN (終了 0)。
- `forbid`   書けば FAIL (逸脱にもできない)。

**`fixed`・`explicit` の省略は FAIL** (出所によらない。plan §4.3 の 4 (a)、2026-10-15 codex diagnose の Major)。省略したキーの値は
実行するバイナリの既定で決まり (`run_case.sh` は `FORGE_BIN` で別ビルドも回す)、今のソースの既定との一致は照合の代わりにならない。
省略そのものを試すときは `--deviate key=~ --reason ...` (「書かない」の逸脱) で宣言する。逸脱の記録の実効値は未確定とする。

どれにも無いキーは、問題 YAML の物理条件・入出力のキー (YAML の `problem_keys`) か段の値に効かないキー (`ignore`) でなければ
**プロファイルに無いキー**として FAIL (逸脱を宣言すれば通す)。加えて `check_solver_config.py` の FAIL、
未知のキー・誤った節のキー (`check_solver_config.unknown_keys`; そちらでは WARN)、YAML の重複キー・merge key (`yaml_strict.py`)、
検査できない入力は FAIL。

継承は 1 本の鎖 (`inherits: <親の名前>`)。鎖の中で同じ段・同じキーに違う値・違う種類を与えたら読み込み時に `RecipeError`
(黙って上書きしない)。`all_stages` はそのプロファイルの全段に当てる。`abstract: true` のプロファイルは宣言できない。

使い方:
    recipe_settings.py list
    recipe_settings.py show PROFILE
    recipe_settings.py check   RUN_DIR|solverConfig.yaml --profile P --stage S [--deviate key=value ... --reason R [--basis B]]
    recipe_settings.py declare RUN_DIR --profile P --stage S [--deviate key=value ... --reason R [--basis B]] [--out RECIPE.json]
    recipe_settings.py declare RUN_DIR --experiment --reason R [--profile P --stage S]
    recipe_settings.py --doc [--write MD | --check MD]

`declare` は照合が通ったときだけ `RECIPE.json` (プロファイル・版・段・逸脱・config の sha256) を書く (終了 0)。FAIL は書かずに終了 1。
「実験」の宣言 (`--experiment`) は推奨への適合を求めず、プロファイルから外れるキーと理由 (1 行の自由記述、必須) だけを書く
(禁止キー・check_solver_config の FAIL・未知のキー・重複キーは実験でも FAIL)。逸脱・実験の根拠 (`--basis`) は任意
(ユーザ決定 2026-10-14「試し計算は手軽にやりたい」)。逸脱の値に `~` を書くと「書かない」(キーを省いた) の宣言
(`RECIPE.json` の逸脱には `effective: 未確定` を付ける)。
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
RECIPES = os.path.join(HERE, 'recipes', 'recommended_settings.yaml')
DOC = os.path.join(ROOT, 'procedures', 'recommended-settings.md')
DOC_BEGIN = '<!-- BEGIN GENERATED: recipe_settings.py --doc (手で編集しない。正本は solver_density_cuda/tools/recipes/recommended_settings.yaml) -->'
DOC_END = '<!-- END GENERATED: recipe_settings.py --doc -->'
DECL_FORMAT = 'forge-recipe-declaration'
DECL_VERSION = 1
KINDS = ('fixed', 'advisory', 'explicit')
UNDETERMINED = '未確定 (省略: 実行するバイナリの既定で決まる)'

sys.path.insert(0, HERE)


class RecipeError(Exception):
    """推奨設定のファイルの不正 (継承の衝突・鎖でない継承・循環・未知のプロファイル・形式の誤り)。"""


# ------------------------------------------------------------------------------------------------ 値の扱い

def norm(v):
    """比較用の正規化: 数値は float、文字列は前後の空白と大文字小文字を無視、列・写像は再帰。"""
    if isinstance(v, bool):
        return ('bool', v)
    if isinstance(v, (int, float)):
        return ('num', float(v))
    if isinstance(v, str):
        return ('str', v.strip().strip('"\'').lower())
    if isinstance(v, (list, tuple)):
        return ('list', tuple(norm(x) for x in v))
    if isinstance(v, dict):
        return ('map', tuple(sorted((str(k), norm(x)) for k, x in v.items())))
    return ('other', repr(v))


def same(a, b):
    return norm(a) == norm(b)


def show(v):
    """表示用の値 (YAML の flow 形式に近い 1 行)。"""
    if v is None:
        return '書かない'
    if isinstance(v, str):
        return v
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return repr(v)
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(show(x) for x in v) + ']'
    if isinstance(v, dict):
        return '{' + ', '.join(f'{k}: {show(x)}' for k, x in v.items()) + '}'
    return repr(v)


def flatten(doc, prefix=''):
    """設定の dict を完全修飾パス → 値 に平らにする (列は値のまま)。"""
    out = {}
    if isinstance(doc, dict):
        for k, v in doc.items():
            p = f'{prefix}.{k}' if prefix else str(k)
            if isinstance(v, dict) and v:
                out.update(flatten(v, p))
            else:
                out[p] = v
    return out


def advisory_ok(spec, v):
    """advisory の範囲: [lo, hi] / {in: [...]} / {min, max} / 'any'。戻り (合否, 範囲の表示)。"""
    if spec == 'any' or spec is None:
        return True, '任意'
    if isinstance(spec, (list, tuple)) and len(spec) == 2:
        lo, hi = spec
        txt = f'{show(lo)} (目安)' if same(lo, hi) else f'{show(lo)}〜{show(hi)} (目安)'
        try:
            x = float(v)
        except (TypeError, ValueError):
            return False, txt
        return (float(lo) <= x <= float(hi)), txt
    if isinstance(spec, dict) and 'in' in spec:
        vals = list(spec['in'])
        txt = (f'{show(vals[0])} (書くなら)' if len(vals) == 1 else ' / '.join(show(c) for c in vals) + ' のどれか')
        return any(same(v, c) for c in vals), txt
    if isinstance(spec, dict) and ('min' in spec or 'max' in spec):
        try:
            x = float(v)
        except (TypeError, ValueError):
            return False, f'{show(spec.get("min"))}〜{show(spec.get("max"))}'
        ok = ('min' not in spec or x >= float(spec['min'])) and ('max' not in spec or x <= float(spec['max']))
        return ok, f'{show(spec.get("min", ""))}〜{show(spec.get("max", ""))}'
    raise RecipeError(f'advisory の範囲の形が不正: {spec!r} ([lo, hi] / {{in: [...]}} / {{min, max}} / any)')


def advisory_text(spec):
    return advisory_ok(spec, 0)[1] if spec != 'any' else '任意'


# ------------------------------------------------------------------------------------------------ solver の既定値

_DEFAULTS = None
_RAW_DEFAULTS = None


def solver_default_strings():
    """config_key_inventory が solverConfig.{cpp,hpp} から抜いた既定値の文字列 {キー: 文字列 or None} (失敗は例外)。"""
    global _RAW_DEFAULTS
    if _RAW_DEFAULTS is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location('cki_rs', os.path.join(HERE, 'config_key_inventory.py'))
        cki = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cki)
        raw = open(cki.SRC_CPP, errors='replace').read() + '\n' + open(cki.SRC_HPP, errors='replace').read()
        _RAW_DEFAULTS = {p: d.get('default') for p, d in cki.extract_keys(cki.strip_comments(raw)).items()}
    return _RAW_DEFAULTS


def solver_defaults():
    """solverConfig が省略時に使う値 (config_key_inventory の抽出から、数値・文字列に読めるものだけ)。
    表示 (§1 との差の「書かない (solver 既定 X)」) にだけ使う。照合には使わない (省略は FAIL、plan §4.3 の 4 (a))。"""
    global _DEFAULTS
    if _DEFAULTS is not None:
        return _DEFAULTS
    out = {}
    try:
        for p, d in solver_default_strings().items():
            v = _parse_default(d)
            if v is not None:
                out[p] = v
    except Exception:
        out = {}
    _DEFAULTS = out
    return out


def _parse_default(s):
    if s is None:
        return None
    s = str(s).strip()
    m = re.fullmatch(r'"([^"]*)"', s) or re.fullmatch(r'std::string\("([^"]*)"\)', s)
    if m:
        return m.group(1)
    s2 = re.sub(r'(?<=[0-9.])[fF]$', '', s)
    try:
        return int(s2)
    except ValueError:
        pass
    try:
        return float(s2)
    except ValueError:
        return None


# ------------------------------------------------------------------------------------------------ 明示すべきキー (#16)

WHEN_OPS = ('not', 'gt', 'min')


def _when_value_ok(spec, v):
    """`when` の 1 条件。spec: 値 (同じ値) / 列 (どれか; `null` = 書いていない) / {not: [...]} / {gt: x} / {min: x}。
    v は文脈の値 (None = 書いていない)。"""
    if isinstance(spec, dict):
        if len(spec) != 1 or set(spec) - set(WHEN_OPS):
            raise RecipeError(f'explicit_keys の when の条件の形が不正: {spec!r} (値 / 列 / {{not: [...]}} / {{gt: x}} / {{min: x}})')
        op, arg = next(iter(spec.items()))
        if op == 'not':
            return not _when_value_ok(list(arg) if isinstance(arg, (list, tuple)) else [arg], v)
        if v is None or isinstance(v, bool):
            return False
        try:
            x = float(v)
        except (TypeError, ValueError):
            return False
        return x > float(arg) if op == 'gt' else x >= float(arg)
    if isinstance(spec, (list, tuple)):
        return any((c is None and v is None) or (c is not None and v is not None and same(c, v)) for c in spec)
    return v is not None and same(spec, v)


def when_match(when, ctx):
    """`when` (写像。全部の条件が成り立てば真) を文脈 ctx (config の完全修飾パス・`applies.X`・`species_count`) に当てる。"""
    if when is not None and not isinstance(when, dict):
        raise RecipeError(f'explicit_keys の when は写像: {when!r}')
    return all(_when_value_ok(spec, ctx.get(k)) for k, spec in (when or {}).items())


def explicit_context(doc, bcond=None):
    """`when` を当てる文脈: config の完全修飾パス、`applies.X` (derive_applies; 壁の種類は bcond から)、`species_count` (化学種の数)。"""
    flat = flatten(doc)
    ctx = dict(flat)
    ctx.update({f'applies.{k}': v for k, v in derive_applies(doc, bcond).items()})
    sp = flat.get('physProp.species')
    ctx['species_count'] = len(sp) if isinstance(sp, (list, tuple)) else 0
    return flat, ctx


def _verify_source(key, case, raw, where):
    """旧い実効値の出所を今のソースで確かめる。合わなければ RecipeError (ソルバの既定が変わった: 旧い実効値を黙って更新しない)。"""
    src = case['source']
    if 'inventory' in src:
        got = raw.get(key)
        if key not in raw or str(got) != str(src['inventory']):
            raise RecipeError(f'{where}: config_key_inventory の抽出した既定 {got!r} が記録 {src["inventory"]!r} と違う '
                              '(ソルバの既定が書き写しの後に変わった。旧い実効値を見直し、plan に書いてから直す)')
        v = _parse_default(got)
        if v is not None and not same(v, case['value']):
            raise RecipeError(f'{where}: 旧い実効値 {show(case["value"])} が config_key_inventory の既定 {show(v)} と違う')
    for s in ([src] if 'file' in src else []) + list(src.get('also') or []):
        try:
            with open(os.path.join(ROOT, str(s['file'])), encoding='utf-8', errors='replace') as f:
                text = f.read()
        except (OSError, KeyError, TypeError) as e:
            raise RecipeError(f'{where}: 出所のファイルを読めない: {e}') from None
        for t in (s.get('text') if isinstance(s.get('text'), list) else [s.get('text')]):
            if not t or str(t) not in text:
                raise RecipeError(f'{where}: 出所 {s["file"]} に {t!r} が無い '
                                  '(ソルバの既定の決め方が書き写しの後に変わった。旧い実効値を見直す)')


def explicit_table(recipes, verify=True):
    """YAML の `explicit_keys` (明示すべきキーの基準、plan §5.1 #16) を検査して返す {キー: {why, basis, cases}}。
    各 case = {when (適用条件), value (旧い実効値), source (出所)}。source は `inventory` (config_key_inventory が抜いた既定値の
    文字列) か `file` + `text` (既定を決めるソースの文字列; inventory が抜けない既定・条件つきの既定。`also: [{file, text}]` で
    別のファイルも確かめる)、`at` (ファイル:行、表示用)。
    verify なら出所を今のソースで確かめる (`_verify_source`)。"""
    table = recipes.get('explicit_keys') or {}
    if not isinstance(table, dict):
        raise RecipeError('explicit_keys は写像 (キー → {why, basis, cases})')
    raw = solver_default_strings() if verify else {}
    out = {}
    for key, ent in table.items():
        key = str(key)
        if not isinstance(ent, dict) or not isinstance(ent.get('cases'), list) or not ent['cases']:
            raise RecipeError(f'explicit_keys.{key}: cases (列) が要る')
        if set(ent) - {'why', 'basis', 'cases'}:
            raise RecipeError(f'explicit_keys.{key}: 未知の項目 {sorted(set(ent) - {"why", "basis", "cases"})}')
        if not str(ent.get('why') or '').strip():
            raise RecipeError(f'explicit_keys.{key}: why (値に依存する理由) が要る')
        if is_problem_key(recipes, key):
            raise RecipeError(f'explicit_keys.{key}: 問題 YAML のキー (problem_keys) は明示すべきキーにしない')
        for i, c in enumerate(ent['cases']):
            where = f'explicit_keys.{key}.cases[{i}]'
            if not isinstance(c, dict) or 'value' not in c or not isinstance(c.get('source'), dict):
                raise RecipeError(f'{where}: value (旧い実効値) と source (出所) が要る')
            if set(c) - {'when', 'value', 'source'}:
                raise RecipeError(f'{where}: 未知の項目 {sorted(set(c) - {"when", "value", "source"})}')
            when_match(c.get('when'), {})                 # 形の検査
            src = c['source']
            if not ('inventory' in src or ('file' in src and 'text' in src)):
                raise RecipeError(f'{where}: source は inventory (抽出した既定値の文字列) か file + text (ソースの文字列)')
            if verify:
                _verify_source(key, c, raw, where)
        out[key] = ent
    return out


def explicit_for(recipes, cfg_text, bcond=None, table=None):
    """段の config が省略している明示すべきキー → 旧い実効値 {キー: 値}。`explicit_keys` の各キーのうち config に書いていないもので、
    最初に when が成り立つ case の値。書き写しの台本 (`recipes/regen_profiles.py`) が段ごとの `explicit` を作るのに使う。"""
    doc, fails = parse_config(cfg_text)
    if fails:
        raise RecipeError(f'config を読めない: {fails}')
    flat, ctx = explicit_context(doc, bcond)
    out = {}
    for key, ent in (explicit_table(recipes) if table is None else table).items():
        if key in flat or any(p.startswith(key + '.') for p in flat):
            continue
        for c in ent['cases']:
            if when_match(c.get('when'), ctx):
                out[key] = c['value']
                break
    return dict(sorted(out.items()))


def explicit_source(recipes, key, value):
    """明示すべきキーの旧い実効値の出所 (表示用の 1 行)。"""
    for c in ((recipes.get('explicit_keys') or {}).get(key) or {}).get('cases') or []:
        if same(c.get('value'), value):
            src = c.get('source') or {}
            if src.get('at'):
                return str(src['at'])
            return 'config_key_inventory' if 'inventory' in src else str(src.get('file', ''))
    return None


# ------------------------------------------------------------------------------------------------ 読み込みと継承

def load_recipes(path=None, text=None):
    """推奨設定のファイルを読む (重複キー・merge key は yaml_strict で拒否)。"""
    import yaml_strict
    if text is None:
        with open(path or RECIPES, encoding='utf-8') as f:
            text = f.read()
    try:
        doc = yaml_strict.load(text)
    except Exception as e:
        raise RecipeError(f'推奨設定のファイルを読めない: {e}') from None
    if not isinstance(doc, dict) or not isinstance(doc.get('profiles'), dict):
        raise RecipeError('推奨設定のファイルに profiles の写像が無い')
    doc['_sha256'] = hashlib.sha256(text.encode()).hexdigest()
    return doc


def problem_keys(recipes):
    return [str(e['key'] if isinstance(e, dict) else e) for e in (recipes.get('problem_keys') or [])]


def is_problem_key(recipes, key):
    """問題 YAML の物理条件・入出力のキー (その写像の下も) か。"""
    return any(key == p or key.startswith(p + '.') for p in problem_keys(recipes))


def _chain(recipes, name):
    profs = recipes['profiles']
    chain, seen, cur = [], [], name
    while cur is not None:
        if cur in seen:
            raise RecipeError(f'継承が循環している: {" → ".join(seen + [cur])}')
        if cur not in profs:
            raise RecipeError(f'プロファイル {cur!r} が無い' + (f' ({seen[-1]} の inherits)' if seen else ''))
        seen.append(cur)
        inh = profs[cur].get('inherits')
        if isinstance(inh, (list, tuple)):
            if len(inh) != 1:
                raise RecipeError(f'{cur}: 継承は 1 本の鎖に限る (inherits に {len(inh)} 個の親: {list(inh)})')
            inh = inh[0]
        cur = inh
    return list(reversed(seen))


def resolve(recipes, name):
    """プロファイルを継承の鎖で解決する。戻り {name, chain, abstract, stages: {段: {fixed, advisory, explicit}}, forbid, ignore, ...}。
    鎖の中で同じ段・同じキーに違う値・違う種類 (fixed/advisory/explicit/forbid/ignore) を与えたら RecipeError。"""
    chain = _chain(recipes, name)
    profs = recipes['profiles']
    per_stage = {}            # 段 → キー → (種類, 値, 出所)
    all_entries = []          # (種類, キー, 値, 出所) — all_stages は最後に全段へ
    forbid, ignore, applies = {}, {}, {}
    stages = []

    def put(table, key, kind, val, origin):
        if key in table:
            k0, v0, o0 = table[key]
            if k0 != kind or (kind != 'ignore' and not same(v0, val)):
                raise RecipeError(f'{name}: 継承の衝突 — {key} を {o0} は {k0} {show(v0)}、{origin} は {kind} {show(val)} '
                                  '(黙って上書きしない。片方を直すこと)')
            return
        table[key] = (kind, val, origin)

    for lvl in chain:
        pr = profs[lvl] or {}
        unknown = set(pr) - {'inherits', 'abstract', 'description', 'applies', 'source', 'forbid', 'ignore', 'all_stages',
                             'stages', 'differs_from_recommended', 'transcribed_from', 'recommended_ref', 'machine'}
        if unknown:
            raise RecipeError(f'{lvl}: 未知の項目 {sorted(unknown)}')
        for k, v in (pr.get('applies') or {}).items():
            if k in applies and not same(applies[k][0], v):
                raise RecipeError(f'{name}: 継承の衝突 — 適用条件 {k} を {applies[k][1]} は {v!r}、{lvl} は {applies[k][0]!r}')
            applies[k] = (v, lvl)
        for k, why in (pr.get('forbid') or {}).items():
            forbid.setdefault(str(k), (why, lvl))
        for k, why in (pr.get('ignore') or {}).items():
            ignore.setdefault(str(k), (why, lvl))
        for kind in KINDS:
            for k, v in ((pr.get('all_stages') or {}).get(kind) or {}).items():
                all_entries.append((kind, str(k), v, f'{lvl}.all_stages'))
        for st, body in (pr.get('stages') or {}).items():
            if st not in stages:
                stages.append(st)
            tbl = per_stage.setdefault(st, {})
            for kind in KINDS:
                for k, v in ((body or {}).get(kind) or {}).items():
                    put(tbl, str(k), kind, v, f'{lvl}.stages.{st}')
            extra = set(body or {}) - set(KINDS)
            if extra:
                raise RecipeError(f'{lvl}.stages.{st}: 未知の項目 {sorted(extra)} (fixed / advisory / explicit だけ)')
    for st in stages:
        tbl = per_stage.setdefault(st, {})
        for kind, k, v, origin in all_entries:
            put(tbl, k, kind, v, origin)
    out_stages = {}
    for st in stages:
        tbl = per_stage[st]
        for k in tbl:
            if k in forbid:
                raise RecipeError(f'{name}: 継承の衝突 — {k} は forbid ({forbid[k][1]}) なのに {tbl[k][2]} が {tbl[k][0]} にしている')
            if k in ignore:
                raise RecipeError(f'{name}: 継承の衝突 — {k} は ignore ({ignore[k][1]}) なのに {tbl[k][2]} が {tbl[k][0]} にしている')
        out_stages[st] = {kind: {k: v for k, (kd, v, _) in sorted(tbl.items()) if kd == kind} for kind in KINDS}
        for k, v in out_stages[st]['advisory'].items():
            advisory_text(v)                          # 形の検査
    leaf = profs[name] or {}
    res = {'name': name, 'chain': chain, 'abstract': bool(leaf.get('abstract', False)),
           'description': leaf.get('description'), 'machine': leaf.get('machine'),
           'applies': {k: v for k, (v, _) in applies.items()}, 'source': leaf.get('source'),
           'transcribed_from': leaf.get('transcribed_from'), 'recommended_ref': leaf.get('recommended_ref'),
           'differs_from_recommended': leaf.get('differs_from_recommended') or [],
           'forbid': {k: w for k, (w, _) in forbid.items()}, 'ignore': {k: w for k, (w, _) in ignore.items()},
           'stages': out_stages, 'version': recipes.get('version')}
    res['sha256'] = profile_hash(res)
    return res


def profile_hash(res):
    """解決したプロファイルの内容のハッシュ (宣言・関門の記録用)。"""
    body = {k: res[k] for k in ('name', 'chain', 'applies', 'forbid', 'ignore', 'stages', 'version')}
    return hashlib.sha256(json.dumps(body, sort_keys=True, ensure_ascii=False, default=str).encode()).hexdigest()


def declarable(recipes):
    return [n for n, p in recipes['profiles'].items() if not (p or {}).get('abstract')]


# ------------------------------------------------------------------------------------------------ 照合

def _dup_path(text):
    """重複キー・merge key の位置を完全修飾パスで探す (yaml_strict が見つけたものを名指しするため)。"""
    import yaml
    try:
        root = yaml.compose(text, Loader=yaml.SafeLoader)
    except Exception:
        return None

    def walk(node, path):
        if isinstance(node, yaml.MappingNode):
            seen = set()
            for k, v in node.value:
                kv = k.value if isinstance(k, yaml.ScalarNode) else str(k)
                if k.tag == 'tag:yaml.org,2002:merge':
                    return '.'.join(path + ['<<'])
                if kv in seen:
                    return '.'.join(path + [kv])
                seen.add(kv)
                r = walk(v, path + [kv])
                if r:
                    return r
        return None
    return walk(root, [])


def parse_config(text):
    """config の文字列を読む。戻り (dict, fails)。重複キー・merge key・読めない入力は fails (名指し)。"""
    import yaml_strict
    try:
        doc = yaml_strict.load(text)
    except yaml_strict.DuplicateKeyError as e:
        return None, [(_dup_path(text) or '(重複キー)', f'YAML の重複キー: {e}')]
    except yaml_strict.MergeKeyError as e:
        return None, [(_dup_path(text) or '<<', f'YAML の merge key: {e}')]
    except Exception as e:
        return None, [('(config)', f'YAML として読めない (検査できない入力): {e}')]
    if not isinstance(doc, dict):
        return None, [('(config)', 'YAML の写像でない (検査できない入力)')]
    return doc, []


def parse_deviations(items, reason=None, basis=None):
    """`key=value` の列 → {key: value}。value は YAML の値として読む (`~` は「書かない」)。理由の無い逸脱は ValueError。"""
    import yaml
    out = {}
    for it in items or []:
        if '=' not in it:
            raise ValueError(f'--deviate は key=value の形 ({it!r})')
        k, v = it.split('=', 1)
        k = k.strip()
        if not k:
            raise ValueError(f'--deviate のキーが空 ({it!r})')
        out[k] = yaml.safe_load(v) if v.strip() else None
    if out and not (reason or '').strip():
        raise ValueError('理由の無い逸脱は宣言できない (--reason に 1 行で理由を書く)')
    return out


def derive_applies(doc, bcond=None):
    """config (と bcond) から適用条件の値を導く (照合・書き写しに使う)。軸対称は solver と同じ優先順
    (`mesh.isAxisymmetric`、無ければ旧 `physProp.isAxisymmetric`: `check_solver_config.effective_axisymmetric`)。"""
    import check_solver_config as csc
    flat = flatten(doc)
    tb = str(flat.get('turbulence.model', 'none')).strip().lower()
    vm = flat.get('physProp.viscMethod')
    physics = ('ns_sst' if tb.startswith('sst') else
               ('euler' if (tb in ('none', '') and vm is not None and int(vm) == 0) else 'ns_laminar'))
    out = {'discretization': str(flat.get('mesh.discretization', 'cell')).strip().strip('"').lower(),
           'axisymmetric': csc.effective_axisymmetric(doc) == 1,
           'gas': 'tp' if int(flat.get('physProp.thermalMethod', 0) or 0) == 2 else 'cpg',
           'physics': physics,
           'time': 'dual_time' if int(flat.get('time.dualTime', 0) or 0) == 1 else
                   ('unsteady' if int(flat.get('time.unsteady', 0) or 0) == 1 else 'steady')}
    if isinstance(bcond, dict):
        kinds = {str((v or {}).get('kind', '')) for v in bcond.values() if isinstance(v, dict)}
        out['wall'] = ('isothermal' if 'wall_isothermal' in kinds else
                       ('noslip_adiabatic' if 'wall' in kinds else 'slip'))
    return out


# 適用条件 → それを決める config のキー (derive_applies と同じ対応)。軸対称は旧 physProp.isAxisymmetric でも決まるが、
# mesh.isAxisymmetric を書けばそちらが優先されるので、段の fixed が mesh.isAxisymmetric を持てば fixed の照合で足りる
# (書いていなければ fixed の省略で FAIL)
APPLIES_KEYS = {'discretization': ('mesh.discretization',), 'axisymmetric': ('mesh.isAxisymmetric',),
                'gas': ('physProp.thermalMethod',), 'physics': ('turbulence.model', 'physProp.viscMethod'),
                'time': ('time.unsteady', 'time.dualTime')}


def _applies_match(want, got):
    if isinstance(want, (list, tuple)):
        return any(same(w, got) for w in want)
    return same(want, got)


def check_config(recipes, cfg_text, profile=None, stage=None, deviations=None, bcond=None, experiment=False):
    """config をプロファイルの段と照合する。戻り dict:
    fails/warns [(キー, 説明)]、applied (宣言どおりの逸脱 [{key, value, profile_value}])、off_profile (実験で外れたキー)、
    profile (解決したプロファイル、無ければ None)、fixed_fails (fails のうち `fixed`・`explicit` の値の不一致・省略によるもの。
    投入の関門 `launch_gate.py` が `FORGE_LAUNCH_GATE=warn` のときだけ WARN に落とす; plan §4.3・§5.1 #10)。"""
    import check_solver_config as csc
    fails, warns, applied, off, fixed_fails = [], [], [], [], []
    deviations = dict(deviations or {})
    res = None
    doc, pf = parse_config(cfg_text)
    fails += pf
    if profile is not None:
        if profile not in recipes['profiles']:
            fails.append(('(profile)', f'プロファイル {profile!r} が無い (有るもの: {", ".join(declarable(recipes))})'))
        else:
            res = resolve(recipes, profile)
            if res['abstract']:
                fails.append(('(profile)', f'{profile} は抽象 (継承の親) で宣言できない'))
            elif stage not in res['stages']:
                fails.append(('(stage)', f'プロファイル {profile} に段 {stage!r} が無い (有る段: {", ".join(res["stages"])})'))
    elif not experiment:
        fails.append(('(profile)', 'プロファイルの宣言が無い (--profile / 実験なら --experiment)'))
    if doc is None:
        return {'fails': fails, 'warns': warns, 'applied': applied, 'off_profile': off, 'profile': res, 'fixed_fails': []}

    # --- check_solver_config (起動時拒否・禁止・dual-time の緩和など) と未知・誤配置のキー (ここでは FAIL) ---
    _, uk_w = csc.unknown_keys(doc)              # 起動時拒否 (unknown_keys の fails) は csc.check の fails に含まれる
    cf, cw = csc.check(doc, bcond)
    ukw_keys = {k for k, _ in uk_w}
    fails += [(k, m) for k, m in cf]
    fails += [(k, m + ' (宣言の照合では FAIL)') for k, m in uk_w]
    flat = flatten(doc)
    warns += [(k, m) for k, m in cw if k not in ukw_keys]
    if res is None or res['abstract'] or stage not in res['stages']:
        if experiment:
            for k in deviations:
                warns.append((k, '実験の宣言では逸脱を個別に宣言しない (外れたキーは自動で記録する)'))
        return {'fails': _dedupe(fails), 'warns': _dedupe(warns), 'applied': applied, 'off_profile': off, 'profile': res,
                'fixed_fails': []}

    st = res['stages'][stage]
    fixed, adv, expl = st['fixed'], st['advisory'], st['explicit']

    def off_or_fail(key, msg, actual, want, kind=None):
        if experiment:
            off.append({'key': key, 'value': actual, 'profile_value': want, 'note': msg})
            return
        if key in deviations:
            dv = deviations[key]
            if dv is None and actual is None:
                # 「書かない」(~) の逸脱: 実効値は実行するバイナリの既定で決まるので未確定として記録する (plan §4.3 の 4 (a))
                applied.append({'key': key, 'value': None, 'profile_value': want, 'effective': UNDETERMINED})
                return
            if dv is not None and actual is not None and same(dv, actual):
                applied.append({'key': key, 'value': actual, 'profile_value': want})
                return
            fails.append((key, f'逸脱の宣言 {show(dv)} と config の値 {show(actual)} が違う'))
            if kind in ('fixed', 'explicit'):
                fixed_fails.append(fails[-1])
            return
        fails.append((key, msg))
        if kind in ('fixed', 'explicit'):
            fixed_fails.append(fails[-1])

    # 禁止 (逸脱にもできない)
    for k, why in res['forbid'].items():
        hit = [p for p in flat if p == k or p.startswith(k + '.')]
        if hit:
            fails.append((k, f'禁止キー ({why})'))
        if k in deviations:
            fails.append((k, '禁止キーは逸脱として宣言できない'))
    # fixed と explicit (明示すべきキー、#16): 値を照合し、省略は FAIL。省略を solver の既定で代用しない
    # (既定は実行するバイナリで決まる。2026-10-15 codex diagnose の Major、plan §4.3 の 4 (a))
    for kind, table in (('fixed', fixed), ('explicit', expl)):
        for k, want in table.items():
            what = (f'プロファイル {profile} の段 {stage} は {show(want)}' if kind == 'fixed' else
                    f'明示すべきキー (plan §5.1 #16): プロファイル {profile} の段 {stage} は {show(want)} '
                    f'(生成器が省略していた値 = 書き写しに使ったソルバの既定、出所 {explicit_source(recipes, k, want) or "explicit_keys"}。'
                    '#11 で生成器に書かせる例外表のキー)')
            if k in flat:
                got = flat[k]
                if not same(got, want):
                    off_or_fail(k, f'{what} (config は {show(got)})', got, want, kind)
            else:
                off_or_fail(k, f'書いていない ({what}。省略の実効値は実行するバイナリで決まるので照合に使わない。'
                               f'省略を試すなら --deviate {k}=~ --reason "...")', None, want, kind)
    # advisory
    for k, spec in adv.items():
        if k in flat:
            ok, rng = advisory_ok(spec, flat[k])
            if not ok:
                if experiment:
                    off.append({'key': k, 'value': flat[k], 'profile_value': rng, 'note': 'advisory の範囲外'})
                elif k in deviations and same(deviations[k], flat[k]):
                    applied.append({'key': k, 'value': flat[k], 'profile_value': rng})
                else:
                    warns.append((k, f'推奨の範囲 {rng} の外 ({show(flat[k])}; プロファイル {profile} の段 {stage})'))
    # プロファイルに無いキー
    known = set(fixed) | set(expl) | set(adv) | set(res['forbid']) | set(res['ignore'])
    for k, v in flat.items():
        if k in known or any(k.startswith(f + '.') for f in res['forbid']) or is_problem_key(recipes, k):
            continue
        if k in ukw_keys:
            continue                                       # 未知・誤配置として既に FAIL
        off_or_fail(k, f'プロファイル {profile} の段 {stage} に無いキー (逸脱を宣言するか、プロファイルに足す)', v, None)
    # 宣言した逸脱の点検
    for k, dv in deviations.items():
        if is_problem_key(recipes, k) and k not in known:
            warns.append((k, '問題 YAML の物理条件・入出力のキーは逸脱の宣言が要らない (宣言を無視した)'))
        elif k not in flat and dv is not None:
            fails.append((k, f'逸脱の宣言 {show(dv)} だが config に書いていない'))
        elif k not in res['forbid'] and k not in {a['key'] for a in applied} and not any(fk == k for fk, _ in fails):
            warns.append((k, f'逸脱の宣言 {show(dv)} が使われていない (config はプロファイル {profile} の段 {stage} どおり)'))
    # 適用条件 (段の fixed がその条件を決めるキーを全部持つ段では fixed の照合に任せる: 例 SERN の層流暖機は SST の機種の段だが
    # turbulence.model none を fixed に持つ)
    got_app = derive_applies(doc, bcond)
    for k, want in res['applies'].items():
        if APPLIES_KEYS.get(k) and all(c in fixed for c in APPLIES_KEYS[k]):
            continue
        if k in got_app and not _applies_match(want, got_app[k]):
            msg = f'適用条件の外: {k} はプロファイル {show(want)}、config は {show(got_app[k])}'
            if experiment:
                off.append({'key': f'applies.{k}', 'value': got_app[k], 'profile_value': want, 'note': msg})
            else:
                fails.append((f'applies.{k}', msg))
    return {'fails': _dedupe(fails), 'warns': _dedupe(warns), 'applied': applied, 'off_profile': off, 'profile': res,
            'fixed_fails': _dedupe(fixed_fails)}


def _dedupe(items):
    seen, out = set(), []
    for k, m in items:
        if (k, m) not in seen:
            seen.add((k, m))
            out.append((k, m))
    return out


# ------------------------------------------------------------------------------------------------ 宣言

def _read_run(run):
    if os.path.isdir(run):
        cfg = os.path.join(run, 'solverConfig.yaml')
        bc = os.path.join(run, 'bcondConfig.yaml')
        return run, cfg, (bc if os.path.exists(bc) else None)
    d = os.path.dirname(os.path.abspath(run))
    bc = os.path.join(d, 'bcondConfig.yaml')
    return d, run, (bc if os.path.exists(bc) else None)


def _load_bcond(path):
    if path is None:
        return None
    try:
        import yaml_strict
        with open(path, encoding='utf-8') as f:
            return yaml_strict.load(f.read())
    except Exception:
        return None


def declare(run, profile=None, stage=None, deviate=None, reason=None, basis=None, experiment=False, out=None,
            recipes=None, write=True):
    """config を照合し、通れば宣言 (RECIPE.json) を書く。戻り (終了コード, 結果 dict, 書いた path or None)。"""
    recipes = recipes or load_recipes()
    run_dir, cfg_path, bc_path = _read_run(run)
    if experiment and not (reason or '').strip():
        return 2, {'fails': [('(reason)', '「実験」の宣言には理由 (--reason、1 行の自由記述) が要る')], 'warns': []}, None
    try:
        devs = parse_deviations(deviate, reason, basis) if not experiment else {}
    except ValueError as e:
        return 2, {'fails': [('(deviate)', str(e))], 'warns': []}, None
    if experiment and deviate:
        return 2, {'fails': [('(deviate)', '「実験」の宣言では --deviate を使わない (外れたキーは自動で記録する)')], 'warns': []}, None
    try:
        with open(cfg_path, 'rb') as f:
            raw = f.read()
    except OSError as e:
        return 1, {'fails': [('(config)', f'solverConfig.yaml を読めない: {e}')], 'warns': []}, None
    text = raw.decode('utf-8', errors='replace')
    r = check_config(recipes, text, profile, stage, devs, _load_bcond(bc_path), experiment=experiment)
    if r['fails']:
        return 1, r, None
    res = r['profile']
    rec = {'format': DECL_FORMAT, 'version': DECL_VERSION, 'kind': 'experiment' if experiment else 'profile',
           'profile': profile, 'stage': stage,
           'profile_version': recipes.get('version'), 'profile_sha256': (res or {}).get('sha256'),
           'recipes_file': os.path.relpath(RECIPES, ROOT), 'recipes_sha256': recipes.get('_sha256'),
           'config': {'file': os.path.basename(cfg_path), 'sha256': hashlib.sha256(raw).hexdigest()},
           'deviations': [{**a, 'reason': reason, 'basis': basis} for a in r['applied']],
           'warnings': [{'key': k, 'message': m} for k, m in r['warns']],
           'declared_at': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
           'tool': 'solver_density_cuda/tools/recipe_settings.py'}
    if experiment:
        rec.update({'reason': reason, 'basis': basis, 'off_profile': r['off_profile']})
    path = None
    if write:
        path = out or os.path.join(run_dir, 'RECIPE.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(rec, f, indent=1, ensure_ascii=False, default=str)
    r['record'] = rec
    return 0, r, path


# ------------------------------------------------------------------------------------------------ 推奨との差 (書き写しの台本が使う)

STAGE_TO_COMMON = {'soft': 'soft', 'mid': 'mid', 'main': 'main'}


def compute_differs(recipes, gen_stages, ref_name, ignore=()):
    """生成器の段の fixed ({段: {キー: 値}}) を共通プロファイル ref_name (§1) と比べ、違うキーの一覧を返す。
    同じ差が複数の段にあれば 1 行にまとめる。共通プロファイルに対応の無い段 (層流暖機・適応段・梯子) は比べない。"""
    ref = resolve(recipes, ref_name)
    defaults = solver_defaults()
    rows = {}
    for st, fixed in gen_stages.items():
        cst = STAGE_TO_COMMON.get(st)
        if cst is None or cst not in ref['stages']:
            continue
        rf, ra = ref['stages'][cst]['fixed'], ref['stages'][cst]['advisory']
        for k in sorted(set(rf) | set(ra)):
            if k in ignore or k in ref['ignore'] or is_problem_key(recipes, k) and k not in rf:
                continue
            got = fixed.get(k)
            if k in rf:
                want = rf[k]
                if got is None:
                    d = defaults.get(k)
                    if d is not None and same(d, want):
                        continue
                    this = '書かない' + (f' (solver 既定 {show(d)})' if d is not None else ' (solver 既定 不明)')
                elif same(got, want):
                    continue
                else:
                    this = show(got)
                rec_txt = show(want)
            else:
                if got is None:
                    continue
                ok, rng = advisory_ok(ra[k], got)
                if ok:
                    continue
                this, rec_txt = show(got), rng
            key = (k, rec_txt, this)
            rows.setdefault(key, []).append(st)
    order = {s: i for i, s in enumerate(['warm_lam', 'warm_adapt', 'soft', 'mid', 'ramp', 'main'])}
    out = []
    for (k, rec_txt, this), sts in sorted(rows.items(), key=lambda kv: (kv[0][0], order.get(kv[1][0], 9))):
        out.append({'key': k, 'stages': sorted(sts, key=lambda s: order.get(s, 9)), 'recommended': rec_txt, 'this': this,
                    'acceptance': '未'})
    return out


# ------------------------------------------------------------------------------------------------ 文書の生成

def _cell(v):
    return show(v).replace('|', '\\|')


def render_doc(recipes):
    """`procedures/recommended-settings.md` の生成区間 (機種と段ごとの値の表)。"""
    profs = recipes['profiles']
    L = [DOC_BEGIN, '']
    L.append(f'推奨設定のファイルの版 {recipes.get("version")}。プロファイルの一覧 (`python3 solver_density_cuda/tools/recipe_settings.py list`):')
    L.append('')
    L.append('| プロファイル | 機種 | 適用条件 | 段 | 根拠 |')
    L.append('| --- | --- | --- | --- | --- |')
    names = [n for n in profs if not (profs[n] or {}).get('abstract')]
    for n in names:
        r = resolve(recipes, n)
        src = r.get('source') or {}
        src_txt = ' '.join(str(src[k]) for k in ('plan', 'section', 'date') if src.get(k))
        app = ', '.join(f'{k} {show(v)}' for k, v in r['applies'].items())
        L.append(f'| `{n}` | {_cell(r.get("machine") or "")} | {_cell(app)} | {" → ".join(r["stages"])} | {_cell(src_txt)} |')
    L.append('')
    L.append('**読み方**: 各表の値は `fixed` (宣言したプロファイルと違えば FAIL、理由つきの逸脱で通す)、`〜` の範囲と「任意」は `advisory`'
             ' (外れたら WARN)。「(明示すべき)」の値は `explicit` (plan §5.1 #16: 生成器が省略していた、値に依存するキーと旧い実効値。'
             '照合は `fixed` と同じ)。**`fixed`・`explicit` の省略は FAIL** (省略の実効値は実行するバイナリで決まるので、solver の既定との'
             '一致で代用しない。省略を試すなら `--deviate key=~ --reason "..."`、実効値は未確定として記録)。'
             '本段 (main) 以外の列の「=」は本段と同じ値、空欄はその段では照合しない (本段の長さ'
             ' `time.last.nStepOuter` は問題 YAML のキー)。表に無いキーは、問題 YAML の物理条件・入出力のキー'
             ' (`problem_keys`) か無視するキー以外は「プロファイルに無いキー」として FAIL (逸脱を宣言すれば通す)。'
             '禁止キーは書けば FAIL。')
    L.append('')
    ek = recipes.get('explicit_keys') or {}
    if ek:
        L.append('**明示すべきキーの基準** (YAML の `explicit_keys`、plan §5.1 #16): 各機種・段で、下の適用条件が成り立ち、'
                 '生成器の config が書いていないキーを、その段の `explicit` に旧い実効値で載せる (`recipes/regen_profiles.py`)。'
                 '旧い実効値の出所は書き写しの時のソースで確かめてある (変わると書き写しの台本が止まる)。')
        L.append('')
        L.append('| キー | 適用条件 | 旧い実効値 | 出所 | 値に依存する理由 |')
        L.append('| --- | --- | --- | --- | --- |')
        for k, ent in ek.items():
            for c in ent.get('cases') or []:
                when = ', '.join(f'{a} {show(b)}' for a, b in (c.get('when') or {}).items()) or '常に'
                src = c.get('source') or {}
                at = src.get('at') or ('config_key_inventory' if 'inventory' in src else src.get('file', ''))
                L.append(f'| `{k}` | {_cell(when)} | {_cell(c.get("value"))} | {_cell(at)} | {_cell(ent.get("why", ""))} |')
        L.append('')
    for n in names:
        r = resolve(recipes, n)
        L.append(f'#### `{n}`')
        L.append('')
        if r.get('description'):
            L.append(str(r['description']).strip())
            L.append('')
        tf = r.get('transcribed_from')
        if tf:
            L.append('書き写し: ' + '、'.join(f'{k} `{v}`' for k, v in tf.items()))
            L.append('')
        sts = list(r['stages'])
        main = 'main' if 'main' in sts else sts[-1]
        keys = sorted({k for st in sts for kind in KINDS for k in r['stages'][st][kind]})
        L.append('| キー | ' + ' | '.join(sts) + ' |')
        L.append('| --- |' + ' --- |' * len(sts))
        for k in keys:
            cells = []
            mv = _val_txt(r['stages'][main], k)
            for st in sts:
                v = _val_txt(r['stages'][st], k)
                cells.append('=' if (st != main and v == mv and v != '') else v)
            L.append(f'| `{k}` | ' + ' | '.join(_cell(c) if c not in ('=', '') else c for c in cells) + ' |')
        L.append('')
        if r['forbid']:
            L.append('禁止: ' + '、'.join(f'`{k}` ({w})' for k, w in r['forbid'].items()))
            L.append('')
        if r['ignore']:
            L.append('照合しない: ' + '、'.join(f'`{k}`' for k in r['ignore']) + ' (' + '; '.join(sorted(set(map(str, r['ignore'].values())))) + ')')
            L.append('')
        dif = r['differs_from_recommended']
        if dif:
            L.append(f'§1 (`{r.get("recommended_ref")}`) との差:')
            L.append('')
            L.append('| キー | 段 | §1 | この機種 | 受入 |')
            L.append('| --- | --- | --- | --- | --- |')
            for d in dif:
                L.append(f'| `{d["key"]}` | {", ".join(d["stages"])} | {_cell(d["recommended"])} | {_cell(d["this"])} | {d["acceptance"]} |')
            L.append('')
    L.append(DOC_END)
    return '\n'.join(L) + '\n'


def _val_txt(stage_tbl, k):
    if k in stage_tbl['fixed']:
        return show(stage_tbl['fixed'][k])
    if k in stage_tbl['advisory']:
        return advisory_text(stage_tbl['advisory'][k])
    if k in stage_tbl.get('explicit', {}):
        return show(stage_tbl['explicit'][k]) + ' (明示すべき)'
    return ''


def doc_block(md_text):
    """文書の生成区間 (BEGIN〜END を含む) を返す。無ければ None。"""
    i, j = md_text.find(DOC_BEGIN), md_text.find(DOC_END)
    if i < 0 or j < 0 or j < i:
        return None
    return md_text[i:j + len(DOC_END)] + '\n'


def update_doc(md_path, recipes, check_only=False):
    """生成区間を YAML から作り直す。check_only なら一致を確かめるだけ。戻り True = 一致 / 書いた。"""
    text = open(md_path, encoding='utf-8').read()
    cur = doc_block(text)
    if cur is None:
        raise RecipeError(f'{md_path} に生成区間の印 ({DOC_BEGIN[:30]}…) が無い')
    new = render_doc(recipes)
    if check_only:
        return cur == new
    i = text.find(DOC_BEGIN)
    text = text[:i] + new + text[i + len(cur):]
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write(text)
    return True


# ------------------------------------------------------------------------------------------------ CLI

def _print_result(r, head):
    print(head)
    for k, m in r['fails']:
        print(f'  FAIL {k}: {m}')
    for k, m in r['warns']:
        print(f'  WARN {k}: {m}')
    for a in r.get('applied') or []:
        eff = f'、実効値 {a["effective"]}' if a.get('effective') else ''
        print(f'  DEVIATION {a["key"]}: {show(a["value"])} (プロファイル {show(a["profile_value"])}{eff})')
    for o in r.get('off_profile') or []:
        print(f'  OFF-PROFILE {o["key"]}: {show(o["value"])} (プロファイル {show(o["profile_value"])})')
    print(f'  VERDICT: {"FAIL" if r["fails"] else ("WARN" if r["warns"] else "PASS")}')


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == '--doc':
        ap = argparse.ArgumentParser(prog='recipe_settings.py --doc')
        ap.add_argument('--doc', action='store_true')
        g = ap.add_mutually_exclusive_group()
        g.add_argument('--write', metavar='MD', help='生成区間を書き換える')
        g.add_argument('--check', metavar='MD', help='生成区間と YAML の一致を確かめる (違えば終了 1)')
        a = ap.parse_args(argv)
        rc = load_recipes()
        if a.write:
            update_doc(a.write, rc)
            print(f'書いた: {a.write}')
            return 0
        if a.check:
            ok = update_doc(a.check, rc, check_only=True)
            print('一致' if ok else f'不一致: {a.check} の生成区間を `recipe_settings.py --doc --write {a.check}` で作り直す')
            return 0 if ok else 1
        sys.stdout.write(render_doc(rc))
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('list')
    sp = sub.add_parser('show')
    sp.add_argument('profile')
    for name in ('check', 'declare'):
        p = sub.add_parser(name)
        p.add_argument('run', help='run ディレクトリ (solverConfig.yaml) か solverConfig.yaml')
        p.add_argument('--profile')
        p.add_argument('--stage')
        p.add_argument('--deviate', action='append', default=[], metavar='KEY=VALUE',
                       help='プロファイルからの逸脱 (`~` は書かない)。--reason 必須')
        p.add_argument('--reason', help='逸脱・実験の理由 (1 行の自由記述)')
        p.add_argument('--basis', help='根拠 (plan の節・run パス; 任意)')
        p.add_argument('--experiment', action='store_true', help='「実験」の宣言 (推奨への適合を求めない; --reason 必須)')
        if name == 'declare':
            p.add_argument('--out', help='宣言の出力 (既定 <run>/RECIPE.json)')
    a = ap.parse_args(argv)
    try:
        recipes = load_recipes()
    except RecipeError as e:
        print(f'推奨設定のファイルが不正: {e}')
        return 2
    if a.cmd == 'list':
        for n in declarable(recipes):
            r = resolve(recipes, n)
            print(f'{n:<34} {" → ".join(r["stages"]):<40} {r.get("machine") or ""}')
        return 0
    if a.cmd == 'show':
        try:
            print(json.dumps(resolve(recipes, a.profile), indent=1, ensure_ascii=False, default=str))
        except RecipeError as e:
            print(f'ERROR: {e}')
            return 2
        return 0
    if not a.experiment and (not a.profile or not a.stage):
        ap.error('--profile と --stage が要る (推奨に従わない run は --experiment --reason "...")')
    if a.profile and not a.stage:
        ap.error('--profile には --stage が要る')
    rc, r, path = declare(a.run, a.profile, a.stage, a.deviate, a.reason, a.basis, a.experiment,
                          getattr(a, 'out', None), recipes=recipes, write=(a.cmd == 'declare'))
    head = f'--- {a.run}  profile {a.profile or "(なし)"}  stage {a.stage or "(なし)"}' + ('  [実験]' if a.experiment else '')
    _print_result(r, head)
    if path:
        print(f'  宣言を書いた: {path}')
    elif a.cmd == 'declare' and rc != 0:
        print('  宣言は書いていない')
    return rc


if __name__ == '__main__':
    sys.exit(main())
