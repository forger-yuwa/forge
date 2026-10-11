#!/usr/bin/env python3
"""区間付きの収束 VERDICT (`CONVERGENCE_SEGMENT.txt`) が、その run の残差の履歴から作られたことの照合
(plan time_integration-line-viscous-jacobian-dt の result 段 codex M3、諮問 2026-10-11-faceh-dt-result-disposition-diagnose)。

保存のテキストだけを読む判定は、残差 CSV が無くても・別の CSV から作られたテキストでも通っていた。ここでは
  (1) stage_manifest.json と、そこに書かれた段の履歴 (段が 1 つで段の履歴が無ければ residual_history.csv) の実在、
  (2) それらと forge_launches.jsonl を **作業用のディレクトリへ複製して** `check_convergence.py --segment` を回し直す
      (run ディレクトリで回すと residual_history_segment.csv を書き、証拠を書き換える)、
  (3) 保存のテキストと回し直しの出力で、VERDICT のトークン・判定区間名・[last step N] の 3 つが一致し、
      N が履歴の最後の outer_end の step と一致すること
を要求する。食い違いは理由付きで返し、合わせに行かない。全文は比べない (ツールの表示の細部に依らないため)。
"""
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TOOL = Path(__file__).resolve().parent / "../../solver_density_cuda/tools/check_convergence.py"
RE_SEG = re.compile(r"判定区間 = (\S+)")
RE_LAST = re.compile(r"\[last step (\d+)\]\s*->\s*([A-Z][A-Z ]+[A-Z])")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def parse(text):
    """(判定区間名, last step, VERDICT) を返す。どれかが無ければ None を入れる。区間名・VERDICT が複数あれば None"""
    segs, lasts = RE_SEG.findall(text), RE_LAST.findall(text)
    seg = segs[0] if len(segs) == 1 else None
    last, verdict = (int(lasts[0][0]), lasts[0][1]) if len(lasts) == 1 else (None, None)
    return seg, last, verdict


def last_outer_end(p):
    n = None
    with open(p, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("phase") == "outer_end":
                n = int(r["step"])
    return n


def recheck(run_dir):
    """戻り値: (ok, 記録の dict)。ok が False なら記録の "why" に理由"""
    d = Path(run_dir)
    rec = {"tool_sha256": sha(TOOL)}
    stored_p = d / "CONVERGENCE_SEGMENT.txt"
    if not stored_p.is_file():
        return False, {**rec, "why": "CONVERGENCE_SEGMENT.txt が無い"}
    stored = parse(stored_p.read_text(errors="replace"))
    rec["stored"] = {"segment": stored[0], "last_step": stored[1], "verdict": stored[2]}
    if None in stored:
        return False, {**rec, "why": "保存のテキストから区間名・last step・VERDICT を 1 つずつ読めない"}
    mp = d / "stage_manifest.json"
    if not mp.is_file():
        return False, {**rec, "why": "stage_manifest.json が無い"}
    stages = json.loads(mp.read_text()).get("stages") or []
    hist = [s.get("history") for s in stages if s.get("history")]
    have = [h for h in hist if (d / h).is_file()]
    if not have:
        if len(stages) == 1 and (d / "residual_history.csv").is_file():
            have = ["residual_history.csv"]          # 段が 1 つの run の扱い (check_convergence.build_segment_csv と同じ)
        else:
            return False, {**rec, "why": f"段の履歴が無い ({hist})"}
    rec["csv_sha256"] = {h: sha(d / h) for h in have}
    with tempfile.TemporaryDirectory(prefix="lvc_segcheck_") as t:
        for f in ["stage_manifest.json", "forge_launches.jsonl"] + have:
            if (d / f).is_file():
                shutil.copy2(d / f, Path(t) / f)
        o = subprocess.run([sys.executable, str(TOOL), t, "--segment"], capture_output=True, text=True)
        fresh = parse(o.stdout + o.stderr)
    rec["fresh"] = {"segment": fresh[0], "last_step": fresh[1], "verdict": fresh[2], "rc": o.returncode}
    if fresh != stored:
        return False, {**rec, "why": f"回し直しが保存のテキストと違う (保存 {stored}、回し直し {fresh})"}
    # 区間が 1 段 (main) の run だけを扱う。段が複数なら step はオフセット付きで連結されるので最後の段の CSV と比べられない
    seg_csv = have[-1] if len(have) == 1 else None
    if seg_csv is None:
        return False, {**rec, "why": f"段の履歴が複数 ({have}) — この照合は 1 段の run 用"}
    n = last_outer_end(d / seg_csv)
    rec["csv_last_outer_end"] = n
    if n != stored[1]:
        return False, {**rec, "why": f"last step {stored[1]} が履歴の最後の outer_end の step {n} と違う"}
    return True, rec


if __name__ == "__main__":
    ok, r = recheck(sys.argv[1])
    print(json.dumps({"ok": ok, **r}, ensure_ascii=False, indent=1))
    sys.exit(0 if ok else 1)
