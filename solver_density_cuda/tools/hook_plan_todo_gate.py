#!/usr/bin/env python3
"""PostToolUse フック: plans/active/*.md を編集したとき、その plan に残作業表が無ければ知らせる。

AGENTS.md 「方針の反映トリガは『実装』ではなく『決定』」の実体化。判定は
`solver_density_cuda/tools/check_plans.py` に委譲し (残作業表 / 変更ログ / README 記載の存在確認のみ)、
**編集した plan だけ**を見る。既存 plan 全部を一度に FAIL にすると無視されるようになるので、
「触ったら直す」段階移行にしてある。

PostToolUse の decision=block は「編集は済んでいて、理由が Claude に返る」だけで、編集は取り消されない。
"""
import sys, os, json, subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
LINT = os.path.join(ROOT, "solver_density_cuda", "tools", "check_plans.py")


def main():
    try:
        ev = json.load(sys.stdin)
    except Exception:
        print("{}"); return

    ti = ev.get("tool_input") or {}
    tr = ev.get("tool_response") or {}
    path = tr.get("filePath") or ti.get("file_path") or ""
    if not path:
        print("{}"); return
    path = os.path.abspath(path)

    active = os.path.join(ROOT, "plans", "active") + os.sep
    if not (path.startswith(active) and path.endswith(".md")):
        print("{}"); return
    if not os.path.isfile(path):
        print("{}"); return

    try:
        p = subprocess.run([sys.executable, LINT, path], capture_output=True, text=True, timeout=20)
    except Exception:
        print("{}"); return
    if p.returncode == 0:
        print("{}"); return

    rel = os.path.relpath(path, ROOT)
    reason = (
        f"{rel} の構造チェックが FAIL です:\n\n{p.stdout.strip()}\n\n"
        "AGENTS.md 「方針の反映トリガは『実装』ではなく『決定』」より: 各 plans/active/*.md は\n"
        "**優先順つきの残作業表を持つ** (雛型は plans/_template.md §5.1)。残作業の実体はここに置き、\n"
        "notes/sessions/ の引き継ぎ文書には写しとポインタだけを置く。\n"
        "この plan を触ったついでに不足分を足すこと (編集自体は成功しています)。")
    print(json.dumps({"decision": "block", "reason": reason}, ensure_ascii=False))


if __name__ == "__main__":
    main()
