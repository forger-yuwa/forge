"""case の README の「計算 run 一覧」に campaign の行を同期する (plan tooling-design-orchestrator §4.5.1・§4.8、AGENTS.md「case の run 索引」)。

campaign ごとに印 `<!-- forge-design:<campaign の ID> begin -->` と `... end -->` の間の行を台帳から作り直す。ほかの行 (人が書いた行・別の campaign の行) は残す。
case ごとのロック (`<case>/.forge-readme.lock` の flock) の下で書く。README が無ければ、ケースの見出しと run 一覧の表の見出しを作る。
"""
from __future__ import annotations

import fcntl
import os
from pathlib import Path

HEADER = "| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |\n| --- | --- | --- | --- |\n"


def _rows(doc: dict) -> list:
    rows = []
    for a in doc.get("attempts") or []:
        if not a.get("run_dir"):
            continue
        it = a.get("interpretation") or {}
        st = "active" if a["id"] not in (doc.get("ingested") or []) else ("ref" if it.get("ok") else "停止")
        rows.append(f"| `{Path(a['run_dir']).name}` | campaign {doc.get('name')} の段 {a['step']} (試行 {a['id']}、recipe "
                    f"{doc['recipe']['name']}/{doc['recipe']['version']}) | {str(it.get('detail') or a.get('state'))[:200]} | {st} |")
    return rows


def sync_runs(case_dir: Path, doc: dict) -> Path:
    case_dir = Path(case_dir)
    rd = case_dir / "README.md"
    lk = os.open(case_dir / ".forge-readme.lock", os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(lk, fcntl.LOCK_EX)
        text = rd.read_text() if rd.is_file() else f"# {case_dir.name}\n\n## 計算 run 一覧\n\n{HEADER}"
        if "## 計算 run 一覧" not in text:
            text = text.rstrip("\n") + f"\n\n## 計算 run 一覧\n\n{HEADER}"
        b, e = f"<!-- forge-design:{doc['campaign_id']} begin -->", f"<!-- forge-design:{doc['campaign_id']} end -->"
        block = "\n".join([b, *_rows(doc), e])
        if b in text and e in text:
            i, j = text.index(b), text.index(e) + len(e)
            text = text[:i] + block + text[j:]
        else:
            lines = text.split("\n")
            h = next(n for n, ln in enumerate(lines) if ln.startswith("## 計算 run 一覧"))
            # 見出しの後の表 (| で始まる行の連続) の最後の行の後ろに足す。表が無ければ表の見出しを作る
            n = h + 1
            while n < len(lines) and not lines[n].startswith("|") and not lines[n].startswith("#"):
                n += 1
            if n < len(lines) and lines[n].startswith("|"):
                while n < len(lines) and (lines[n].startswith("|") or lines[n].startswith("<!-- forge-design:")):
                    n += 1
                lines[n:n] = block.split("\n")
            else:
                lines[h + 1:h + 1] = ["", *HEADER.rstrip("\n").split("\n"), *block.split("\n")]
            text = "\n".join(lines)
        tmp = rd.with_name(f".{rd.name}.tmp")
        tmp.write_text(text)
        os.replace(tmp, rd)
        return rd
    finally:
        fcntl.flock(lk, fcntl.LOCK_UN)
        os.close(lk)
