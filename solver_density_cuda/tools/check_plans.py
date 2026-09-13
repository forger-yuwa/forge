#!/usr/bin/env python3
"""plans/active/*.md の構造 lint (AGENTS.md 「方針の反映トリガは『実装』ではなく『決定』」の実体化)。

チェック内容 (いずれも誤検知ゼロを狙った機械的な存在確認のみ。内容の良し悪しは見ない):

  1. 残作業表      — 見出しに「残作業」を含む節があり、実体 (表行 `|` か箇条書き) が 1 行以上ある
  2. 変更ログ      — `## ... 変更ログ` の節がある
  3. README 記載   — plans/README.md の一覧にファイル名が出てくる
  4. レビュー記録  — 見出しに「レビュー記録」を含む節があり (雛型 §6.1)、status に応じた codex レビュー行がある:
                     in_progress 以上 → `plan` 行、done → `result` 行 (免除は `plan 免除` / `result 免除` + 理由)。
                     行の記録欄に書かれた notes/reviews/*.md は実在すること (AGENTS.md 「codex レビュー」)

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
STATUS_RE = re.compile(r"\*\*status\*\*\s*[:：]\s*`?([A-Za-z_\-]+)`?")
REVIEW_ROW_RE = re.compile(r"^\s*\|\s*(plan|result)\s*(免除)?\s*\|", re.IGNORECASE)
REVIEW_FILE_RE = re.compile(r"notes/reviews/[\w\-.]+\.md")


def plan_status(text):
    m = STATUS_RE.search(text)
    if not m:
        return ""
    return m.group(1).lower().replace("-", "_")


def review_rows(secs):
    """レビュー記録節の表行 → [(stage, exempt, cells)]"""
    rows = []
    for h, b in secs:
        if "レビュー記録" not in h:
            continue
        for ln in b:
            m = REVIEW_ROW_RE.match(ln)
            if not m:
                continue
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            rows.append((m.group(1).lower(), bool(m.group(2)), cells))
    return rows


def check_reviews(path, text, secs):
    """[(ok, ラベル, 補足)] — レビュー記録の存在と status 整合。"""
    has_sec = any("レビュー記録" in h for h, _ in secs)
    if not has_sec:
        return [(False, "レビュー記録", "「レビュー記録」の節が無い (plans/_template.md §6.1。codex レビューの記録表)")]
    rows = review_rows(secs)
    st = plan_status(text)
    out = []
    need = []
    if st in ("in_progress", "done"):
        need.append("plan")
    if st == "done":
        need.append("result")
    for stage in need:
        have = [r for r in rows if r[0] == stage]
        if not have:
            out.append((False, "レビュー記録",
                        f"status `{st}` なのに `{stage}` 段の codex レビュー行が無い "
                        f"(`{stage}` 行を足すか、`{stage} 免除` 行 + 理由。実行: codex_review.py <plan> --stage {stage})"))
    for stage, exempt, cells in rows:
        if exempt:
            reason = cells[-1] if len(cells) >= 4 else ""
            if not reason or reason.startswith("<"):
                out.append((False, "レビュー記録", f"`{stage} 免除` 行に理由が無い"))
            continue
        joined = " ".join(cells)
        m = REVIEW_FILE_RE.search(joined)
        if not m:
            out.append((False, "レビュー記録", f"`{stage}` 行の記録欄に notes/reviews/*.md のパスが無い"))
        elif not os.path.isfile(os.path.join(ROOT, m.group(0))):
            out.append((False, "レビュー記録", f"`{stage}` 行の記録 {m.group(0)} が実在しない"))
    if not out:
        out.append((True, "レビュー記録", ""))
    return out


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

    res.extend(check_reviews(path, text, secs))
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
