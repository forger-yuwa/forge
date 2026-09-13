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
"""
import argparse, os, re, sys

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
    ("", "bodyForce"): "一様体積力 [N/m³] ([fx, fy, fz])。周期チャネル等の駆動に使う",
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
        member = _member_before(cpp, pos)
        cur = keys.get((section, key))
        if cur is None:
            keys[(section, key)] = dict(
                type=typ, default=default, member=member,
                desc=KEY_DESC.get((section, key)) or desc.get(member) or desc.get(key, ""),
                required=required, cond=cond, lines=[_line_of(cpp, pos)])
        else:
            cur["lines"].append(_line_of(cpp, pos))
            if cur["default"].startswith("(") and not default.startswith("("):
                cur.update(type=typ, default=default)
            if not cur["desc"]:
                cur["desc"] = KEY_DESC.get((section, key)) or desc.get(member) or desc.get(key, "")
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
    sections = {s for (s, _k) in keys}
    dep_keys = set(load_deprecated())
    for sk in [(s, k) for (s, k), i in keys.items()
               if i["type"] == "?" and ((("%s.%s" % (s, k)) if s else k) in sections or k in dep_keys)]:
        del keys[sk]

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


def walk(node, section=""):
    """config を (ドット付き節名, キー, 値) に平坦化。トップレベルのスカラーは節名 ''。"""
    if not isinstance(node, dict):
        return
    for k, v in node.items():
        if isinstance(v, dict):
            yield from walk(v, "%s.%s" % (section, k) if section else k)
        else:
            yield (section, k, v)


def fmt(v):
    """YAML として読み戻したとき値が変わらない表記にする (引用符・バックスラッシュを含む文字列対策)。"""
    import yaml
    return yaml.safe_dump(v, default_flow_style=True, allow_unicode=True, width=10**6).strip().rstrip("\n").removesuffix("...").strip()


def same_as_default(val, default):
    d = default.strip().rstrip("f")
    try:
        return abs(float(val) - float(d)) <= 1e-12 * max(1.0, abs(float(d)))
    except (TypeError, ValueError):
        return str(val).strip('"') == d.strip('"')


def cmd_check(a):
    keys, misses = load_keys()
    dep = load_deprecated()
    bad = 0
    for arg in a.paths:
        p = config_path(arg)
        print("== %s" % os.path.relpath(p, REPO))
        seen = list(walk(read_yaml(p)))
        have = {(s, k) for (s, k, _v) in seen}
        for section, key, val in seen:
            info = keys.get((section, key))
            label = "%s.%s" % (section, key) if section else key
            where = sorted({s2 for (s2, k2) in keys if k2 == key})
            if key in dep:
                print("  [廃止]   %-28s = %-12s %s" % (label, fmt(val), dep[key]))
                bad += 1
            elif info is not None:
                if not info["required"] and not info["default"].startswith("(") \
                        and not same_as_default(val, info["default"]):
                    print("  [非既定] %-28s = %-12s (既定 %s) %s" % (label, fmt(val), info["default"], info["desc"][:60]))
            elif where:
                print("  [節違い] %-28s = %-12s ここでは読まれない。正しい節は %s"
                      % (label, fmt(val), " / ".join(("%s.%s" % (w, key)) if w else ("トップレベルの " + key) for w in where)))
                bad += 1
            else:
                print("  [未知]   %-28s = %-12s ソルバは読まない (綴り間違い/旧キー/別ツール用)" % (label, fmt(val)))
                bad += 1
        for (s2, k2), info in sorted(keys.items()):
            if info["required"] and not info["cond"] and (s2, k2) not in have:
                print("  [必須欠落] %-26s 現ソルバでは必須 (無いと起動時エラー)  %s"
                      % ("%s.%s" % (s2, k2) if s2 else k2, info["desc"][:40]))
                bad += 1
        if a.missing:
            for (s, k), info in sorted(keys.items()):
                if (s, k) not in have and (s or k) and info["desc"]:
                    print("  [未指定] %-28s   既定 %-10s %s"
                          % ("%s.%s" % (s, k) if s else k, info["default"], info["desc"][:60]))
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


def annotated_lines(cfg, keys, dep, section="", depth=0):
    """元の構造と値を保ったまま、各キーの上に説明と既定値のコメントを差し込む。"""
    out, ind = [], "  " * depth
    for k, v in cfg.items():
        path = "%s.%s" % (section, k) if section else k
        if isinstance(v, dict):
            out.append("%s%s:" % (ind, k))
            out.extend(annotated_lines(v, keys, dep, path, depth + 1))
            if depth == 0:
                out.append("")
            continue
        info = keys.get((section, k))
        if k in dep:
            out.append("%s# [廃止] %s" % (ind, dep[k]))
        elif info is None:
            where = sorted({s2 for (s2, k2) in keys if k2 == k})
            out.append("%s# [%s] %s" % (ind, "節違い" if where else "未知",
                                        ("正しい節は " + " / ".join(w or "トップレベル" for w in where))
                                        if where else "ソルバは読まない"))
        else:
            out.extend(wrap_comment(info["desc"], ind) if info["desc"] else [])
            out.append("%s# 既定 %s  (型 %s)" % (ind, info["default"], info["type"]))
        out.append("%s%s: %s" % (ind, k, fmt(v)))
    return out


def cmd_annotate(a):
    import yaml
    keys, _ = load_keys()
    dep = load_deprecated()
    p = os.path.abspath(config_path(a.path))
    cfg = read_yaml(p)
    dst = os.path.abspath(a.out or os.path.join(os.path.dirname(p), "solverConfig.annotated.yaml"))
    if dst == p or (os.path.exists(dst) and os.path.samefile(dst, p)):
        print("拒否: 出力先が入力と同じファイル (%s)。annotate は写しを作るもので、元 config は書き換えない。" % dst,
              file=sys.stderr)
        return 2
    head = ["# solverConfig.yaml の注釈つきの写し (tools/config_doc.py annotate が生成; 実行には使わない)",
            "# 説明と既定値は solver_density_cuda/input/solverConfig.{cpp,hpp} から抽出したもの。",
            "# 元: %s" % os.path.relpath(p, REPO), ""]
    body = "\n".join(head + annotated_lines(cfg, keys, dep)) + "\n"
    # 生成物を読み戻して値が変わっていないことを確かめる (引用符・バックスラッシュ・数値表記の事故防止)
    back = yaml.safe_load(body)
    if back != cfg:
        print("拒否: 生成した写しを読み戻すと元の config と一致しない。出力しない。", file=sys.stderr)
        return 2
    open(dst, "w", encoding="utf-8").write(body)
    print("wrote %s" % os.path.relpath(dst, REPO))
    return 0


def cmd_template(a):
    keys, _ = load_keys()
    dep = load_deprecated()
    print("# solverConfig.yaml の注釈つき雛形 (tools/config_doc.py template が生成)")
    print("# 全キーを列挙してある。**使うキーだけ残す** こと (既定のままのキーは書かない方が読みやすい)。")
    print("# 廃止キー・使用禁止キーは除いてある。推奨値は procedures/recommended-settings.md が正本 (ここはコード上の既定値)。")
    live = {(s, k): i for (s, k), i in keys.items() if k not in dep}
    sections = sorted({s for (s, _k) in live})
    if a.section:
        sections = [s for s in sections if s == a.section or s.startswith(a.section + ".")]
    printed = set()
    for sec in sections:
        parts = sec.split(".") if sec else []
        for d in range(len(parts)):
            head = ".".join(parts[:d + 1])
            if head not in printed:
                print("\n%s%s:" % ("  " * d, parts[d]))
                printed.add(head)
        if not sec:
            print("\n# --- トップレベル ---")
        ind = "  " * len(parts)
        for (s2, k), info in sorted(live.items()):
            if s2 != sec:
                continue
            if info["desc"]:
                for ln in wrap_comment(info["desc"], ind):
                    print(ln)
            if info["default"].startswith("("):
                val = "REQUIRED" if info["required"] else "SEE_CODE"
                print("%s%s: %s   # %s (%s) — 値を入れてから使う" % (ind, k, val, info["default"], info["type"]))
            else:
                print("%s%s: %s" % (ind, k, info["default"]))
    return 0


def cmd_list(a):
    keys, misses = load_keys()
    print("| セクション | キー | 型 | 既定値 | 説明 |")
    print("| --- | --- | --- | --- | --- |")
    for (s, k), info in sorted(keys.items()):
        if a.section and s != a.section:
            continue
        print("| %s | `%s` | %s | `%s` | %s |" % (s or "(top)", k, info["type"], info["default"], info["desc"]))
    nd = sum(1 for i in keys.values() if not i["desc"])
    print("\n抽出 %d キー。説明なし %d キー。解析できなかった読み出し %d 件。" % (len(keys), nd, len(misses)))
    if a.missing_desc:
        print("\n説明が無いキー (solverConfig.hpp に行末コメントを足すと出る):")
        for (s, k), info in sorted(keys.items()):
            if not info["desc"]:
                print("  %s.%s" % (s, k) if s else "  %s" % k)
    return 1 if misses else 0


def cmd_coverage(a):
    keys, misses = load_keys()
    print("抽出 %d キー / 解析できなかった読み出し %d 件" % (len(keys), len(misses)))
    for line, why in sorted(misses):
        print("  input/solverConfig.cpp:%d  %s" % (line, why))
    if a.verbose:
        print("\n読み出し位置 (input/solverConfig.cpp の行):")
        for (s, k), info in sorted(keys.items()):
            print("  %-40s %s" % ("%s.%s" % (s, k) if s else k, ",".join(str(x) for x in sorted(set(info["lines"])))))
    return 1 if misses else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("check"); p.add_argument("paths", nargs="+"); p.add_argument("--missing", action="store_true"); p.set_defaults(fn=cmd_check)
    p = sp.add_parser("annotate"); p.add_argument("path"); p.add_argument("-o", "--out"); p.set_defaults(fn=cmd_annotate)
    p = sp.add_parser("template"); p.add_argument("--section"); p.set_defaults(fn=cmd_template)
    p = sp.add_parser("list"); p.add_argument("--section"); p.add_argument("--missing-desc", action="store_true"); p.set_defaults(fn=cmd_list)
    p = sp.add_parser("coverage"); p.add_argument("-v", "--verbose", action="store_true"); p.set_defaults(fn=cmd_coverage)
    a = ap.parse_args()
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
