#!/usr/bin/env python3
"""重複キーを拒否する YAML ローダー (全階層・flow 形式 `{a: 1, a: 2}` と block 形式の両方)。

PyYAML (`yaml.safe_load`) は同じ mapping の重複キーを**後勝ち**で黙って読むが、solver の yaml-cpp は**先勝ち**である
(`migrate_solver_config.py` の `unsafe_reason` と同じ認識)。重複キーを含む config を PyYAML で読んで判断すると、
solver の実効設定と違うものを見ることになる (例: `turbulence: {model: "sst", model: "none"}` を PyYAML は none、
solver は sst と読む; plan tooling-rerun-conditions codex result 段 3 回目 #1)。そこで入力検査に使う読み込みは
本モジュールの `load` に揃え、重複キーがあれば `DuplicateKeyError` で止める。

    from yaml_strict import load, DuplicateKeyError
    doc = load(text)              # 重複キーがあれば DuplicateKeyError (yaml.YAMLError の派生)

merge key (`<<: *a`) で取り込んだキーを明示キーで上書きするのは YAML の仕様どおりの上書きなので重複とみなさない
(検査は merge 展開の前に、明示キーどうしだけで行う)。
"""
import yaml


class DuplicateKeyError(yaml.YAMLError):
    """同じ mapping に同じキーが 2 回以上ある (PyYAML は後勝ち・yaml-cpp は先勝ちで解釈が割れる)。"""


class _StrictLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader, node, deep=False):
    seen = {}
    for k, _ in node.value:
        if k.tag == "tag:yaml.org,2002:merge":
            continue                                   # `<<` は merge 展開 (重複検査の対象外)
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
    """`yaml.safe_load` と同じ結果を返す。ただし全階層の mapping で重複キーがあれば DuplicateKeyError。"""
    return yaml.load(text, Loader=_StrictLoader)
