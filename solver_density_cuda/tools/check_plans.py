#!/usr/bin/env python3
"""plans/active/*.md の構造 lint (AGENTS.md 「方針の反映トリガは『実装』ではなく『決定』」の実体化)。

チェック内容 (いずれも誤検知ゼロを狙った機械的な存在確認のみ。内容の良し悪しは見ない):

  1. 残作業表      — 見出しに「残作業」を含む節があり、実体 (表行 `|` か箇条書き) が 1 行以上ある
  2. 変更ログ      — `## ... 変更ログ` の節がある
  3. README 記載   — plans/README.md の一覧にファイル名が出てくる

使い方:
  python3 solver_density_cuda/tools/check_plans.py                 # plans/active/*.md を全部見る
  python3 solver_density_cuda/tools/check_plans.py PATH [PATH ...] # 指定ファイルだけ見る
  python3 solver_density_cuda/tools/check_plans.py --summary       # 件数だけ (移行状況の把握用)

終了コード: 0 = 全 PASS / 1 = FAIL あり。
"""
import sys, os, re, glob

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
ACTIVE = os.path.join(ROOT, "plans", "active")
README = os.path.join(ROOT, "plans", "README.md")

HEAD_RE = re.compile(r"^(#{2,4})\s+(.*)$")


def sections(text):
    """[(見出し文字列, 本文行リスト)] を返す。"""
    out, cur, body = [], None, []
    for line in text.splitlines():
        m = HEAD_RE.match(line)
        if m:
            if cur is not None:
                out.append((cur, body))
            cur, body = m.group(2), []
        elif cur is not None:
            body.append(line)
    if cur is not None:
        out.append((cur, body))
    return out


def check(path, readme_text):
    """[(ok, ラベル, 補足)] を返す。"""
    text = open(path, encoding="utf-8").read()
    secs = sections(text)
    res = []

    todo = [(h, b) for h, b in secs if "残作業" in h]
    if not todo:
        res.append((False, "残作業表", "見出しに「残作業」を含む節が無い (plans/_template.md §5.1 を参照)"))
    else:
        has_body = any(ln.lstrip().startswith(("|", "-", "*")) for _, b in todo for ln in b)
        res.append((has_body, "残作業表", "" if has_body else "節はあるが実体 (表行/箇条書き) が無い"))

    log = any("変更ログ" in h for h, _ in secs)
    res.append((log, "変更ログ", "" if log else "「変更ログ」の節が無い"))

    listed = os.path.basename(path) in readme_text
    res.append((listed, "README 記載", "" if listed else "plans/README.md の一覧に無い"))
    return res


def main(argv):
    summary = "--summary" in argv
    args = [a for a in argv if not a.startswith("--")]
    paths = args or sorted(glob.glob(os.path.join(ACTIVE, "*.md")))
    readme_text = open(README, encoding="utf-8").read() if os.path.exists(README) else ""

    n_fail = 0
    for p in paths:
        p = os.path.abspath(p)
        if not (os.path.isfile(p) and p.endswith(".md")):
            continue
        res = check(p, readme_text)
        bad = [r for r in res if not r[0]]
        if bad:
            n_fail += 1
        if summary:
            continue
        rel = os.path.relpath(p, ROOT)
        if bad:
            print(f"FAIL {rel}")
            for _, label, note in bad:
                print(f"       - {label}: {note}")
        else:
            print(f"PASS {rel}")

    total = len([p for p in paths if p.endswith(".md")])
    print(f"\nVERDICT: {'FAIL' if n_fail else 'PASS'}  ({total - n_fail}/{total} plans OK)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
