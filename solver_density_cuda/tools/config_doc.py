#!/usr/bin/env python3
"""solverConfig.yaml を「読める」ようにする: キーの意味・既定値・注意をソースから機械生成して付ける。

なぜツールなのか: run ディレクトリは複製で増える (case/16 だけで 100 以上)。config に手でコメントを書くと
複製で運ばれ、既定値や意味が変わっても更新されない (廃止キー `LESorRANS` が複製で生き残った前例)。
正本はコード側に置いたまま、必要なときに生成する。

抽出元 (どちらもコードと同じ場所にあるので構造は腐らない):
  input/solverConfig.cpp  … 読み出し式から (節, キー, 型, 既定値, 必須か) を取る。対応する書き方は 5 つ:
      1) getOptionalValidatedValue<T>(node, "key", default, "section")
      2) getValidatedValue<T>(node, "key"[, "section"])
      3) node["key"].as<T>()            (既定は if-else で書かれるのでコードにしか無い)
      4) for (... : node["key"])        (配列キー)
      5) if (node["key"])               (存在判定のみ)
     node は `auto out = config["output"];` のような束縛を追跡して節名に直す。
  input/solverConfig.hpp  … メンバ宣言のコメントを説明として取る (行末 + 桁を揃えた継続行、または同じ字下げの前置き)。
条件付きの必須 (`physProp.chemistry.mechanismFile` は `enabled: 1` のときだけ要る 等) はここでは検査しない。
ソルバ自身が起動時に理由付きで弾くので (solverConfig.cpp の `requires` / `is required` の throw)、
同じ条件をツールに写すと二重の正本になって腐るため。
廃止キーは solverConfig.cpp 自身の拒否メッセージ ("Key 'x' ... is no longer supported"、
`for (const char* old_key : {...})` の列挙を含む) と、ルール由来の内蔵リストから取る。

**説明の中身はコードから来ない** (人が書いたコメントなので陳腐化しうる)。`gpu` を「GPU 番号」と書いて
いた誤りが実際にあった (正しくは 0/1 フラグ)。説明を足すときは値の消費箇所まで当たること。

使い方:
  config_doc.py check    <run_dir|yaml> [...]  未知キー・廃止キー・節違い・必須欠落・既定からの逸脱 (問題があれば exit 1)
  config_doc.py annotate <run_dir|yaml>        注釈付きの写しを solverConfig.annotated.yaml に出す (元は触らない)
  config_doc.py template [--section s]         全キー入りの注釈付き雛形を stdout に出す (廃止キーは出さない)
  config_doc.py list     [--section s]         キー・節・既定値・説明の一覧 (docs 生成用)
  config_doc.py coverage                       抽出できなかった読み出し箇所を行番号付きで出す (0 件であるべき)
  config_doc.py selftest                       期待値つきの自己試験 (過去に壊れた入力をそのまま残してある)
"""
import argparse, os, re, sys, unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "input")
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))

# ルールで禁止されているキー (コードは受け付けるが使ってはいけないもの)。AGENTS.md / recommended-settings.md §9 由来。
DEPRECATED_FALLBACK = {
    "bndFirstOrder": "使用禁止 (粘性応力を破壊・疑似 2D で全域に効く)。段階起動で代替",
    "nodeAxisDirichlet": "廃止 (2026-08-16, 書くと起動エラー)",
    "nodeMidpointFx": "廃止 (2026-08-16)",
    "nodeValueAtNode": "廃止 (2026-08-16)",
    "nodeReconEdgeMidpoint": "廃止 (2026-08-16)",
    "nodeAxisUrDirichlet": "廃止 (2026-08-16)",
}

# メンバ宣言に説明を置けないキー (配列読み・ローカル変数受け) の説明はここに書く
KEY_DESC = {
    ("physProp", "species"): "混合を構成する化学種名のリスト。並び順が index s (roY{s}) を決める。物性は speciesDBFile から引く",
    ("space", "uRef"): "移流ゲージの参照速度 [m/s] ([ux, uy, uz])。space.roRef と対で KEEP の一様流保存に使う",
    ("bodyForce",): "一様体積力 [N/m³] ([fx, fy, fz])。周期チャネル等の駆動に使う",
    ("turbulence", "model"): "乱流モデル: none | wale | sigma | sst | sst-ddes | sst-iddes (旧 LESorRANS/RANSmodel は起動エラー)",
    ("mesh", "wallDistExtraPhysIDs"): "壁距離の計算に含める非 wall 境界の physical ID リスト (変換時, procedures/solver-settings.md)",
    ("output", "extraFields"): "output.level で落とされる場のうち個別に追加したい名前のリスト (output_cellValNames にあるもの)",
}

TYPE = r"<((?:[^<>]|<(?:[^<>]|<[^<>]*>)*>)*)>"        # std::vector<std::string> のような入れ子型
NODE = r"([A-Za-z_]\w*(?:\[\"\w+\"\])*)"               # config["a"]["b"] / out / ch


def _strip_line_comments(cpp):
    """行頭コメント (コメントアウトされた旧コード) を消す。行番号は保つ。"""
    return "\n".join("" if re.match(r"^\s*//", ln) else ln for ln in cpp.splitlines())


def _desc_map(hpp):
    """solverConfig.hpp のメンバ宣言に付いたコメントを説明として拾う。

    優先順は (1) 行末コメント + その桁に揃えて続く `//` 行、(2) 宣言と同じ字下げで直前に置かれた `//` ブロック。
    桁で見分けるのは、右端に流した継続行と「次のメンバの前置き」が同じ見た目になるため。
    """
    decl = r"^(\s+)(?:int|double|bool|std::string|float|flow_float)\s+(\w+)\s*(?:=[^;]*)?;\s*(?://\s*(.+?)\s*)?$"
    desc, lines = {}, hpp.splitlines()

    def cmt_col(t):
        if re.match(r"^\s*//\s*-{3,}", t):
            return None
        m = re.match(r"^(\s*)//\s*\S", t)
        return len(m.group(1)) if m else None

    for i, ln in enumerate(lines):
        m = re.match(decl, ln)
        if not m:
            continue
        indent, name, tail = len(m.group(1)), m.group(2), m.group(3)
        parts = []
        if tail:
            col = ln.index("//")
            parts.append(tail)
            j = i + 1
            while j < len(lines):
                c = cmt_col(lines[j])
                if c is None or abs(c - col) > 2:
                    break
                parts.append(re.sub(r"^\s*//\s?", "", lines[j]).strip())
                j += 1
        else:
            j, pre = i - 1, []
            while j >= 0:
                c = cmt_col(lines[j])
                if c is None or abs(c - indent) > 2:
                    break
                pre.append(re.sub(r"^\s*//\s?", "", lines[j]).strip())
                j -= 1
            parts = list(reversed(pre))
        if parts:
            desc[name] = " ".join(parts).strip()
    return desc


def _bindings(cpp):
    """`auto out = config["output"];` / `const YAML::Node ch = physProp["chemistry"];` を節パスに解く。"""
    b = {"config": ""}
    for _ in range(4):     # 束縛が束縛を参照するので数回まわす
        for m in re.finditer(r"(?:auto|(?:const\s+)?YAML::Node)\s+(\w+)\s*=\s*" + NODE + r"\s*;", cpp):
            name, expr = m.group(1), m.group(2)
            base = re.match(r"^[A-Za-z_]\w*", expr).group(0)
            if base not in b:
                continue
            parts = [p for p in [b[base]] if p] + re.findall(r"\[\"(\w+)\"\]", expr)
            b[name] = ".".join(parts)
    return b


def _resolve(expr, binds):
    base = re.match(r"^[A-Za-z_]\w*", expr).group(0)
    if base not in binds:
        return None
    parts = [p for p in [binds[base]] if p] + re.findall(r"\[\"(\w+)\"\]", expr)
    return ".".join(parts)


def _line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def _member_before(cpp, pos):
    m = re.search(r"this->(\w+)\s*=\s*$", cpp[max(0, pos - 140):pos])
    return m.group(1) if m else ""


def load_keys():
    """(節, キー) -> {type, default, member, desc, required, cond, lines} と、拾えなかった読み出しの一覧を返す。

    節はドット付きのフルパス ("time.deltaT")。トップレベルは ""。同じキーが複数箇所で読まれる場合
    (`detectNaN` は time.deltaT とトップレベルの両方) は、実在する読み出し位置をそれぞれ登録する。
    """
    raw_cpp = open(os.path.join(SRC, "solverConfig.cpp"), encoding="utf-8").read()
    hpp = open(os.path.join(SRC, "solverConfig.hpp"), encoding="utf-8").read()
    cpp = _strip_line_comments(raw_cpp)
    desc = _desc_map(hpp)
    binds = _bindings(cpp)

    keys, seen_pos, misses = {}, set(), []

    def put(section, key, typ, default, pos, required, cond):
        if section is None:
            return False
        path = tuple(section.split(".")) + (key,) if section else (key,)
        member = _member_before(cpp, pos)
        cur = keys.get(path)
        if cur is None:
            keys[path] = dict(
                type=typ, default=default, member=member,
                desc=KEY_DESC.get(path) or desc.get(member) or desc.get(key, ""),
                required=required, cond=cond, lines=[_line_of(cpp, pos)])
        else:
            cur["lines"].append(_line_of(cpp, pos))
            if cur["default"].startswith("(") and not default.startswith("("):
                cur.update(type=typ, default=default)
            if not cur["desc"]:
                cur["desc"] = KEY_DESC.get(path) or desc.get(member) or desc.get(key, "")
            if required and not cur["required"]:
                cur.update(required=True, cond=cond)
        return True

    def indent_of(pos):
        bol = cpp.rfind("\n", 0, pos) + 1
        return len(cpp[bol:]) - len(cpp[bol:].lstrip())

    # 1) 省略可能キー
    for m in re.finditer(r"getOptionalValidatedValue" + TYPE + r"\(\s*" + NODE +
                         r"\s*,\s*\"([^\"]+)\"\s*,\s*(.+?)\s*,\s*\"([^\"]*)\"\s*\)", cpp, re.S):
        typ, node, key, default, sec = m.groups()
        sec = sec or _resolve(node, binds)
        if not put(sec, key, typ, re.sub(r"\s+", " ", default).strip(), m.start(), False, False):
            misses.append((_line_of(cpp, m.start()), "節が解決できない: %s" % node))
        seen_pos.add(m.start())
    # 2) 必須キー (字下げ 8 = read() 直下 = 無条件、それより深ければ条件付き)
    for m in re.finditer(r"(?<!Optional)getValidatedValue" + TYPE + r"\(\s*" + NODE +
                         r"\s*,\s*\"([^\"]+)\"(?:\s*,\s*\"([^\"]*)\")?\s*\)", cpp, re.S):
        typ, node, key, sec = m.groups()
        sec = sec if sec is not None else _resolve(node, binds)
        cond = indent_of(m.start()) > 8
        if not put(sec, key, typ, "(必須)" if not cond else "(条件付き必須)", m.start(), True, cond):
            misses.append((_line_of(cpp, m.start()), "節が解決できない: %s" % node))
        seen_pos.add(m.start())
    # 3) 直接読み  4) 配列の for  5) 存在判定だけ
    for pat, kind in ((NODE + r"\[\"(\w+)\"\]\s*\.as" + TYPE + r"\(\)", "as"),
                      (r"for\s*\([^)]*:\s*" + NODE + r"\[\"(\w+)\"\]\s*\)", "for"),
                      (NODE + r"\[\"(\w+)\"\]", "exists")):
        for m in re.finditer(pat, cpp):
            g = m.groups()
            node, key = (g[0], g[1])
            typ = g[2] if kind == "as" else ("list" if kind == "for" else "?")
            sec = _resolve(node, binds)
            if sec is None:
                misses.append((_line_of(cpp, m.start()), "節が解決できない: %s[\"%s\"]" % (node, key)))
                continue
            put(sec, key, typ, "(コード側で分岐)", m.start(), False, False)
    # 存在判定だけで拾ったもののうち、(a) 節そのもの `if (config["mesh"])` と (b) 廃止キーの検出
    # `if (last["nStep"])` は「設定キー」ではないので落とす
    sections = {path[:-1] for path in keys if len(path) > 1}
    dep_keys = set(load_deprecated())
    for path in [pp for pp, i in keys.items()
                 if i["type"] == "?" and (pp in sections or pp[-1] in dep_keys)]:
        del keys[path]

    # 拾えなかった getValidated* 呼び出し (正規表現の取りこぼし検出)
    for m in re.finditer(r"get(?:Optional)?ValidatedValue\s*<", cpp):
        if m.start() not in seen_pos:
            misses.append((_line_of(cpp, m.start()), "読み出し式を解析できない"))
    return keys, misses


def load_deprecated():
    """廃止キー = (1) コード自身が拒否するもの + (2) ルールで禁止された内蔵リスト。

    (1) は `"Key 'x' ... is no longer supported"` と、`for (const char* old_key : {"A","B"})` のように
    実行時に連結して投げる形の両方を拾う。recommended-settings.md §9 の表は自然文なので機械抽出しない
    (「`thermoHrefTemp` 未指定が発散要因」のような行を『キーが廃止』と誤読するため)。
    """
    out = dict(DEPRECATED_FALLBACK)
    try:
        cpp = _strip_line_comments(open(os.path.join(SRC, "solverConfig.cpp"), encoding="utf-8").read())
    except OSError:
        return out
    for m in re.finditer(r"\"Key '(\w+)'(?: in '([^']*)')? is no longer supported\.?\s*([^\"]*)\"", cpp):
        key, _sec, tail = m.groups()
        out[key] = ("起動エラー: " + tail.strip()) if tail.strip() else "起動エラー (コードが拒否)"
    for m in re.finditer(r"for\s*\([^)]*:\s*\{([^}]*)\}\s*\)\s*\{(.{0,600}?)\n\s*\}", cpp, re.S):
        listed, body = m.groups()
        if "no longer supported" not in body:
            continue
        tail = re.search(r"is no longer supported\.?\s*([^\"]*)\"", body)
        note = ("起動エラー: " + tail.group(1).strip()) if tail and tail.group(1).strip() else "起動エラー (コードが拒否)"
        for k in re.findall(r"\"(\w+)\"", listed):
            out.setdefault(k, note)
    return out


def read_yaml(path):
    import yaml
    return yaml.safe_load(open(path, encoding="utf-8"))


def config_path(arg):
    return os.path.join(arg, "solverConfig.yaml") if os.path.isdir(arg) else arg


def disp(path):
    return ".".join(path)


def _width(t):
    """端末での表示幅 (日本語は 2 桁)。列を揃えるため。"""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in t)


def pad(t, n):
    return t + " " * max(1, n - _width(t))


def walk(node, path=()):
    """config を (パスのタプル, 値, 辞書か) に平坦化。辞書そのものも 1 件として出す
    (`mesh.bndFirstOrder: {}` のような空の節や、スカラーのキーに辞書を書いた事故を見逃さないため)。"""
    if not isinstance(node, dict):
        return
    for k, v in node.items():
        p = path + (str(k),)
        if isinstance(v, dict):
            yield (p, v, True)
            yield from walk(v, p)
        else:
            yield (p, v, False)


def fmt(v):
    """表示用。YAML として読み戻したとき値が変わらない表記にする。"""
    import yaml
    t = yaml.safe_dump(v, default_flow_style=True, allow_unicode=True, width=10**6).strip()
    return t.removesuffix("...").strip()


def same_as_default(val, default):
    d = default.strip().rstrip("f")
    try:
        return abs(float(val) - float(d)) <= 1e-12 * max(1.0, abs(float(d)))
    except (TypeError, ValueError):
        return str(val).strip('"') == d.strip('"')


def check_config(cfg, keys, dep, show_missing=False):
    """点検結果を (種別, 表示名, 値, 説明) の列で返す。種別が '非既定'/'未指定' 以外なら要修正。"""
    out = []
    sections = {path[:-1] for path in keys if len(path) > 1}
    have = {p for (p, _v, _d) in walk(cfg)}
    flagged = set()          # 廃止・型違いの辞書の下は掘らない (同じ誤りを子の数だけ繰り返さない)
    unknown_sec = set()      # 未知の節の下は「正しい節はここ」と言える子だけ出す
    for path, val, is_dict in walk(cfg):
        if any(path[:i] in flagged for i in range(1, len(path))):
            continue
        under_unknown = any(path[:i] in unknown_sec for i in range(1, len(path)))
        name, label = path[-1], disp(path)
        same_name = sorted({pp[:-1] for pp in keys if pp[-1] == name})
        if name in dep:
            out.append(("廃止", label, val, dep[name]))
            if is_dict:
                flagged.add(path)
        elif is_dict:
            if path in keys:
                out.append(("型違い", label, val, "スカラーのキーに辞書が書かれている (ソルバは読めない)"))
                flagged.add(path)
            elif path not in sections:
                if not under_unknown:
                    out.append(("未知", label, val, "ソルバはこの節を読まない (綴り間違い/旧キー/別ツール用)"))
                unknown_sec.add(path)
        elif path in keys:
            info = keys[path]
            if not info["required"] and not info["default"].startswith("(") \
                    and not same_as_default(val, info["default"]):
                out.append(("非既定", label, val, "(既定 %s) %s" % (info["default"], info["desc"][:60])))
        elif same_name:
            out.append(("節違い", label, val, "ここでは読まれない。正しい節は %s"
                        % " / ".join((disp(w + (name,)) if w else "トップレベルの " + name) for w in same_name)))
        elif not under_unknown:
            out.append(("未知", label, val, "ソルバは読まない (綴り間違い/旧キー/別ツール用)"))
    for path, info in sorted(keys.items()):
        if info["required"] and not info["cond"] and path not in have:
            out.append(("必須欠落", disp(path), None, "現ソルバでは必須 (無いと起動時エラー)  " + info["desc"][:40]))
    if show_missing:
        for path, info in sorted(keys.items()):
            if path not in have and info["desc"]:
                out.append(("未指定", disp(path), None, "既定 %-10s %s" % (info["default"], info["desc"][:60])))
    return out


OK_KINDS = ("非既定", "未指定")


def cmd_check(a):
    keys, misses = load_keys()
    dep = load_deprecated()
    bad = 0
    for arg in a.paths:
        p = config_path(arg)
        print("== %s" % os.path.relpath(p, REPO))
        for kind, label, val, note in check_config(read_yaml(p), keys, dep, a.missing):
            shown = "" if val is None else fmt(val)
            if isinstance(val, dict):
                shown = "{...}" if val else "{}"
            print("  %s%s%s%s" % (pad("[%s]" % kind, 12), pad(label, 30),
                                  pad(("= " + shown) if val is not None else "", 16), note))
            if kind not in OK_KINDS:
                bad += 1
    if misses:
        print("\n警告: 解析できなかった読み出しが %d 件ある (config_doc.py coverage で確認)。点検は不完全。" % len(misses))
        bad += 1
    return 1 if bad else 0


def wrap_comment(text, ind, width=96):
    out, cur = [], ""
    for w in text.split():
        if cur and len(cur) + 1 + len(w) > width:
            out.append("%s# %s" % (ind, cur)); cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        out.append("%s# %s" % (ind, cur))
    return out


def same_yaml(a, b):
    """2 つの YAML テキストが同じ値を表すか。ブロックスカラー内への挿入など「原文復元では見えない」破壊を捕まえる。"""
    import yaml
    try:
        return yaml.safe_load(a) == yaml.safe_load(b)
    except yaml.YAMLError:
        return False


def _notes_for(path, ind, keys, dep, sections):
    """1 キーの上に差し込むコメント行。節そのものには何も付けない。"""
    name = path[-1]
    if name in dep:
        return ["%s# [廃止] %s" % (ind, dep[name])]
    if path in keys:
        info = keys[path]
        lines = wrap_comment(info["desc"], ind) if info["desc"] else []
        return lines + ["%s# 既定 %s  (型 %s)" % (ind, info["default"], info["type"])]
    if path in sections:
        return []
    same_name = sorted({pp[:-1] for pp in keys if pp[-1] == name})
    if same_name:
        return ["%s# [節違い] 正しい節は %s" % (ind, " / ".join(disp(w) or "トップレベル" for w in same_name))]
    return ["%s# [未知] ソルバは読まない" % ind]


def _key_positions(src):
    """YAML を構文解析して (パス, 行, 桁) を集める。

    行を正規表現で見るだけだと、ブロックスカラー (`|-`) の中の `drive: mesh.h5` をキーと誤認して
    文字列の内側に `#` を挿し込む・複数行のフロー形式 (`deltaT: {..., \n  dt_min: ...}`) の入れ子を
    取り違える、という事故が起きる (codex result-2 M1/M2)。構文木の位置を使えばどちらも起きない。
    """
    import yaml
    out = []

    def rec(node, path):
        if isinstance(node, yaml.MappingNode):
            for k, v in node.value:
                if not isinstance(k, yaml.ScalarNode):
                    continue
                p = path + (str(k.value),)
                out.append((p, k.start_mark.line, k.start_mark.column))
                rec(v, p)

    rec(yaml.compose(src), ())
    return out


def annotate_text(src, keys, dep):
    """**元のテキストを保ったまま**、キーの行の前にコメント行だけを挿入する。

    - 行頭に単独で現れるキーには、その字下げで説明ブロックを差し込む。
    - 1 行に複数キーが並ぶフロー形式 (`space: {convMethod: 1, limiter: 2}`) では、その行の前に
      フルパス付きの 1 行要約をまとめて置く (行の中には手を入れない)。
    返り値は (出力行, 挿入した行かどうかのフラグ)。
    """
    sections = {path[:-1] for path in keys if len(path) > 1}
    lines = src.splitlines()
    per_line = {}
    for path, line, col in _key_positions(src):
        per_line.setdefault(line, []).append((col, path))

    out, added = [], []
    for i, raw in enumerate(lines):
        entries = sorted(per_line.get(i, []))
        if entries:
            indent = raw[:len(raw) - len(raw.lstrip())]
            own_line = len(entries) == 1 and entries[0][0] == len(indent)
            if own_line:
                notes = _notes_for(entries[0][1], indent, keys, dep, sections)
            else:
                notes = []
                for _col, path in entries:
                    for ln in _notes_for(path, indent, keys, dep, sections):
                        notes.append(ln.replace("# ", "# %s: " % disp(path), 1) if ln.strip().startswith("#") else ln)
            for ln in notes:
                out.append(ln); added.append(True)
        out.append(raw); added.append(False)
    return out, added


def cmd_annotate(a):
    keys, _ = load_keys()
    dep = load_deprecated()
    p = os.path.abspath(config_path(a.path))
    src = open(p, encoding="utf-8").read()
    dst = os.path.abspath(a.out or os.path.join(os.path.dirname(p), "solverConfig.annotated.yaml"))
    if dst == p or (os.path.exists(dst) and os.path.samefile(dst, p)):
        print("拒否: 出力先が入力と同じファイル (%s)。annotate は写しを作るもので、元 config は書き換えない。" % dst,
              file=sys.stderr)
        return 2
    lines, added = annotate_text(src, keys, dep)
    # 二重に確かめる: (1) 挿入行を外すと元テキストに戻る (2) YAML として読んだ値が同じ
    # (1) だけでは、ブロックスカラーの内側に `#` を挿した場合を見逃す
    kept = "\n".join(ln for ln, ad in zip(lines, added) if not ad)
    body = "\n".join(lines) + "\n"
    if kept.rstrip("\n") != src.rstrip("\n") or not same_yaml(body, src):
        print("拒否: 生成した写しが元 config と一致しない (原文復元または YAML 値)。出力しない。", file=sys.stderr)
        return 2
    head = ["# solverConfig.yaml の注釈つきの写し (tools/config_doc.py annotate が生成; 実行には使わない)",
            "# 説明と既定値は solver_density_cuda/input/solverConfig.{cpp,hpp} から抽出したもの。",
            "# 元: %s" % os.path.relpath(p, REPO), ""]
    open(dst, "w", encoding="utf-8").write("\n".join(head + lines) + "\n")
    print("wrote %s" % os.path.relpath(dst, REPO))
    return 0


def cmd_template(a):
    keys, _ = load_keys()
    dep = load_deprecated()
    print("# solverConfig.yaml の注釈つき雛形 (tools/config_doc.py template が生成)")
    print("# 全キーを列挙してある。**使うキーだけ残す** こと (既定のままのキーは書かない方が読みやすい)。")
    print("# 廃止キー・使用禁止キーは除いてある。推奨値は procedures/recommended-settings.md が正本 (ここはコード上の既定値)。")
    live = {path: i for path, i in keys.items() if path[-1] not in dep}
    sections = sorted({path[:-1] for path in live})
    if a.section:
        want = tuple(a.section.split("."))
        sections = [s2 for s2 in sections if s2[:len(want)] == want]
    printed = set()
    for sec in sections:
        for d in range(len(sec)):
            head = sec[:d + 1]
            if head not in printed:
                print("\n%s%s:" % ("  " * d, sec[d]))
                printed.add(head)
        if not sec:
            print("\n# --- トップレベル ---")
        ind = "  " * len(sec)
        for path, info in sorted(live.items()):
            if path[:-1] != sec:
                continue
            if info["desc"]:
                for ln in wrap_comment(info["desc"], ind):
                    print(ln)
            if info["default"].startswith("("):
                val = "REQUIRED" if info["required"] else "SEE_CODE"
                print("%s%s: %s   # %s (%s) — 値を入れてから使う" % (ind, path[-1], val, info["default"], info["type"]))
            else:
                print("%s%s: %s" % (ind, path[-1], info["default"]))
    return 0


def cmd_list(a):
    keys, misses = load_keys()
    print("| セクション | キー | 型 | 既定値 | 説明 |")
    print("| --- | --- | --- | --- | --- |")
    for path, info in sorted(keys.items()):
        sec = disp(path[:-1])
        if a.section and sec != a.section:
            continue
        print("| %s | `%s` | %s | `%s` | %s |" % (sec or "(top)", path[-1], info["type"], info["default"], info["desc"]))
    nd = sum(1 for i in keys.values() if not i["desc"])
    print("\n抽出 %d キー。説明なし %d キー。解析できなかった読み出し %d 件。" % (len(keys), nd, len(misses)))
    if a.missing_desc:
        print("\n説明が無いキー (solverConfig.hpp に行末コメントを足すと出る):")
        for path, info in sorted(keys.items()):
            if not info["desc"]:
                print("  %s" % disp(path))
    return 1 if misses else 0


def cmd_coverage(a):
    keys, misses = load_keys()
    print("抽出 %d キー / 解析できなかった読み出し %d 件" % (len(keys), len(misses)))
    for line, why in sorted(misses):
        print("  input/solverConfig.cpp:%d  %s" % (line, why))
    if a.verbose:
        print("\n読み出し位置 (input/solverConfig.cpp の行):")
        for path, info in sorted(keys.items()):
            print("  %-40s %s" % (disp(path), ",".join(str(x) for x in sorted(set(info["lines"])))))
    return 1 if misses else 0


# 期待値つきの自己試験。過去に実際に壊れた入力をそのまま残してある (codex レビュー 2 回分の指摘)。
SELFTEST = [
    ("節違い (同名キーがトップレベルにあっても見逃さない)",
     "space:\n  convMethod: 1\n  keepDissCoeff: 0.05\n", [("節違い", "space.keepDissCoeff")]),
    ("未知キー (綴り間違い)",
     "turbulence:\n  model: \"sst\"\n  kInf: 1.0\n", [("未知", "turbulence.kInf")]),
    ("配列キーは誤検知しない",
     "mesh:\n  wallDistExtraPhysIDs: [6]\n", []),
    ("使用禁止キー (値が空の辞書でも見逃さない)",
     "mesh:\n  bndFirstOrder: {}\n", [("廃止", "mesh.bndFirstOrder")]),
    ("ドット入りのキー名を入れ子と同一視しない",
     "time.deltaT:\n  cfl: 1\n", [("未知", "time.deltaT"), ("節違い", "time.deltaT.cfl")]),
    ("未知の節の下でも、正しい節が分かる子は出す",
     "time:\n  implicit:\n    nLoop: 10\n    cfl_pseudo: 0.5\n",
     [("未知", "time.implicit"), ("節違い", "time.implicit.cfl_pseudo")]),
    ("スカラーのキーに辞書",
     "space:\n  limiter: {a: 1}\n", [("型違い", "space.limiter")]),
]

SELFTEST_TEXT = [
    ("ブロックスカラーの中身をキーと誤認しない",
     "mesh:\n  meshFileName: |-\n    drive: mesh.h5\n"),
    ("複数行のフロー形式でも入れ子を取り違えない",
     "time:\n  deltaT: {control: 1, dt: 1e-5,\n           dt_min: 1e-8, dt_max: 1.0}\n"),
    ("キーの後ろのコメントを入れ子開始と誤認しない",
     "time: # 時間設定\n  unsteady: 0\n"),
    ("YAML 真偽値に見える化学種名 (NO) を壊さない",
     "physProp:\n  species: [\"MIXDRY\", NO]\n  cp: 1039.0\n"),
    ("空の節を潰さない",
     "output: {}\nturbulence:\n  model: \"sst\"\n"),
    ("引用符・バックスラッシュを含む文字列",
     "mesh:\n  meshFileName: \"a'b\\\\c.h5\"\n"),
]


# フロー形式の実 config でも、正しいキーに正しい注釈が付くこと (誤った [節違い] を付けない)
SELFTEST_NOTES = [
    ("フロー形式のキーにフルパス付きの注釈が付く",
     "time:\n  deltaT: {control: 1, dt: 1e-5,\n           dt_min: 1e-8, dt_max: 1.0}\n",
     ["time.deltaT.dt_min", "time.deltaT.dt_max"]),
]


def cmd_selftest(a):
    import io, yaml
    keys, misses = load_keys()
    dep = load_deprecated()
    ng = 0
    for title, text, expect in SELFTEST:
        got = [(k, lbl) for (k, lbl, _v, _n) in check_config(yaml.safe_load(text), keys, dep)
               if k not in OK_KINDS and k != "必須欠落"]
        ok = sorted(got) == sorted(expect)
        print("%s %s" % ("PASS" if ok else "FAIL", title))
        if not ok:
            print("     期待 %s / 実際 %s" % (sorted(expect), sorted(got)))
            ng += 1
    for title, text in SELFTEST_TEXT:
        lines, added = annotate_text(text, keys, dep)
        kept = "\n".join(ln for ln, ad in zip(lines, added) if not ad)
        body = "\n".join(lines) + "\n"
        ok = kept.rstrip("\n") == text.rstrip("\n") and same_yaml(body, text)
        print("%s %s (原文復元 + YAML 値一致)" % ("PASS" if ok else "FAIL", title))
        if not ok:
            ng += 1
    for title, text, expect in SELFTEST_NOTES:
        lines, added = annotate_text(text, keys, dep)
        got = [ln.strip() for ln, ad in zip(lines, added) if ad]
        ok = all(any(e in g for g in got) for e in expect) and not any("[節違い]" in g for g in got if "dt_min" in g)
        print("%s %s" % ("PASS" if ok else "FAIL", title))
        if not ok:
            print("     期待に含む %s / 実際 %s" % (expect, got))
            ng += 1
    print("%s coverage: 解析できなかった読み出し %d 件" % ("PASS" if not misses else "FAIL", len(misses)))
    ng += 1 if misses else 0
    nd = sum(1 for i in keys.values() if not i["desc"])
    print("%s 説明なし %d キー / 抽出 %d キー" % ("PASS" if nd == 0 else "FAIL", nd, len(keys)))
    ng += 1 if nd else 0
    print("\nVERDICT: %s" % ("PASS" if ng == 0 else "FAIL (%d 件)" % ng))
    return 1 if ng else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("check"); p.add_argument("paths", nargs="+"); p.add_argument("--missing", action="store_true"); p.set_defaults(fn=cmd_check)
    p = sp.add_parser("annotate"); p.add_argument("path"); p.add_argument("-o", "--out"); p.set_defaults(fn=cmd_annotate)
    p = sp.add_parser("template"); p.add_argument("--section"); p.set_defaults(fn=cmd_template)
    p = sp.add_parser("list"); p.add_argument("--section"); p.add_argument("--missing-desc", action="store_true"); p.set_defaults(fn=cmd_list)
    p = sp.add_parser("coverage"); p.add_argument("-v", "--verbose", action="store_true"); p.set_defaults(fn=cmd_coverage)
    p = sp.add_parser("selftest"); p.set_defaults(fn=cmd_selftest)
    a = ap.parse_args()
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
