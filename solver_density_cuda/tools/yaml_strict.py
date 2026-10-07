#!/usr/bin/env python3
"""重複キーを拒否する YAML ローダー (全階層・flow 形式 `{a: 1, a: 2}` と block 形式の両方)。

PyYAML (`yaml.safe_load`) は同じ mapping の重複キーを**後勝ち**で黙って読むが、solver の yaml-cpp は**先勝ち**である
(`migrate_solver_config.py` の `unsafe_reason` と同じ認識)。重複キーを含む config を PyYAML で読んで判断すると、
solver の実効設定と違うものを見ることになる (例: `turbulence: {model: "sst", model: "none"}` を PyYAML は none、
solver は sst と読む; plan tooling-rerun-conditions codex result 段 3 回目 #1)。そこで入力検査に使う読み込みは
本モジュールの `load` に揃え、重複キーがあれば `DuplicateKeyError` で止める。

    from yaml_strict import load, DuplicateKeyError
    doc = load(text)              # 重複キーがあれば DuplicateKeyError (yaml.YAMLError の派生)

merge key (`<<`) は**全階層で拒否**する (`MergeKeyError`)。solver の yaml-cpp は merge を展開せず、`node[key]` で
直接キーを探して未定義なら既定値を使う (`solverConfig.cpp` の `getOptionalValidatedValue`) ので、PyYAML が展開した
値を「検査した実効設定」とすると solver の方程式と食い違う (例: `turbulence: {<<: {sstEnergyIncludesK: 1}}` を PyYAML は
1、solver は既定 0 と読む)。また merge の中身の重複キー (`{<<: {model: sst, model: none}}`) も検査を迂回できた
(plan tooling-rerun-conditions codex result 段 4 回目 #1)。anchor / alias (`&a` / `*a`) 単体は yaml-cpp も同じ値に解決するので許す。
"""
import yaml


class DuplicateKeyError(yaml.YAMLError):
    """同じ mapping に同じキーが 2 回以上ある (PyYAML は後勝ち・yaml-cpp は先勝ちで解釈が割れる)。"""


class MergeKeyError(yaml.YAMLError):
    """merge key (`<<`) がある (solver の yaml-cpp は展開しないので PyYAML と実効設定が割れる)。"""


class _StrictLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader, node, deep=False):
    seen = {}
    for k, _ in node.value:
        if k.tag == "tag:yaml.org,2002:merge":
            raise MergeKeyError(
                f"merge key `<<` ({k.start_mark.line + 1} 行目; solver の yaml-cpp は merge を展開しないので "
                "PyYAML の解釈と実効設定が割れる。展開した明示キーで書き直すこと)")
        key = loader.construct_object(k, deep=True)
        try:
            hash(key)
        except TypeError:
            continue                                   # 非 hashable のキーは下の construct_mapping が拒否する
        if key in seen:
            raise DuplicateKeyError(
                f"重複キー {key!r} ({seen[key] + 1} 行目と {k.start_mark.line + 1} 行目; "
                "PyYAML は後勝ち・solver の yaml-cpp は先勝ちで解釈が割れる)")
        seen[key] = k.start_mark.line
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


_StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def load(text):
    """`yaml.safe_load` と同じ結果を返す。ただし全階層の mapping で重複キーがあれば DuplicateKeyError、
    merge key `<<` があれば MergeKeyError。"""
    return yaml.load(text, Loader=_StrictLoader)


def _find_value_node(root, path):
    """compose した木で path (キーの tuple) の値ノードを返す。無い・途中が mapping でなければ None。"""
    node = root
    for key in path:
        if not isinstance(node, yaml.MappingNode):
            return None
        hit = [v for k, v in node.value if isinstance(k, yaml.ScalarNode) and k.value == key]
        if len(hit) != 1:
            return None                                # 0 回 (無い) / 2 回以上 (load が DuplicateKeyError で先に止める)
        node = hit[0]
    return node


def _get(doc, path):
    v = doc
    for k in path:
        if not isinstance(v, dict) or k not in v:
            return None
        v = v[k]
    return v


def replace_scalars(text, updates):
    """updates {キーの tuple パス: 書き込むトークン (文字列)} を、YAML 上のその位置の**値スカラーのトークンだけ**置換する。

    正規表現で `key: 値` を探すとコメント (`# previous nStepOuter: 6000`) に当たり、`{nStepOuter : 6000}` (コロン前の空白)・
    引用符付きキー (`"nStepOuter": 6000`) に当たらない (plan tooling-rerun-conditions codex result 段 4 回目 #3)。
    ここでは compose した木で値ノードを特定し、その文字範囲 (start_mark.index〜end_mark.index) を書き換える
    (コメント・書式・flow/block の形式はそのまま)。書き換え後に読み直し、(1) 各パスの値が「トークンを単独で読んだ値」と一致、
    (2) それ以外の値が不変 であることを検査する。違反・パスが無い・値がスカラーでない・anchor/alias/tag 付きの値は ValueError
    (重複キー・merge key は load と同じく DuplicateKeyError / MergeKeyError)。"""
    import copy
    before = load(text)
    root = yaml.compose(text, Loader=_StrictLoader)
    spans = []
    for path, token in updates.items():
        path = tuple(path)
        node = _find_value_node(root, path)
        if node is None:
            raise ValueError(f"{'.'.join(path)} が無い (書き換えられない)")
        if not isinstance(node, yaml.ScalarNode):
            raise ValueError(f"{'.'.join(path)} の値がスカラーでない")
        s, e = node.start_mark.index, node.end_mark.index
        span = text[s:e]
        if span[:1] in ("&", "*", "!") or node.tag not in (
                "tag:yaml.org,2002:int", "tag:yaml.org,2002:float", "tag:yaml.org,2002:str"):
            raise ValueError(f"{'.'.join(path)} の値 {span!r} は anchor / alias / tag 付き・数値/文字列以外 (共有実体を書き換えうるので扱わない)")
        if "\n" in span:
            raise ValueError(f"{'.'.join(path)} の値が複数行にまたがる ({span!r})")
        spans.append((s, e, path, str(token)))
    spans.sort()
    for (s0, e0, p0, _), (s1, _e1, p1, _t) in zip(spans, spans[1:]):
        if s1 < e0:
            raise ValueError(f"{'.'.join(p0)} と {'.'.join(p1)} の値が同じ位置 (alias で共有している)")
    out = text
    for s, e, _p, token in reversed(spans):
        out = out[:s] + token + out[e:]
    after = load(out)
    want = copy.deepcopy(before)
    for path, token in updates.items():
        path = tuple(path)
        expect = load(str(token))
        got = _get(after, path)
        if got != expect or type(got) is not type(expect):
            raise ValueError(f"書き換え後に読み直したら {'.'.join(path)} = {got!r} (要求 {expect!r})")
        node = want
        for k in path[:-1]:
            node = node[k]
        node[path[-1]] = expect
    if want != after:
        raise ValueError("書き換えで要求したパス以外の値が変わった")
    return out
