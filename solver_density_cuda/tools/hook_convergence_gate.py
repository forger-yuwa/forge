#!/usr/bin/env python3
"""Stop フック: forge を回したのに収束チェックしていない run があればターン終了を block する。

case/*/run_*/residual_history.csv のうち「最近 (既定 180 分) 更新された」ものを対象に、
同じ run の CONVERGENCE_VERDICT.txt が無い / residual より古い ものを検出したら block する。
これは forge の実行方法 (直接/ラッパー/バックグラウンド) に依らずファイル状態だけで判定するので、
「位置安定だけ見て収束と判断する」近道を構造的に防ぐ。

block 時は decision=block + reason を返し、解消法 (run_case.sh 経由 or check_convergence.py で
VERDICT を残す) と「収束/一致を主張するなら VERDICT 行を引用」を伝える。

**セッション別スコープ (2026-09-20)**: 共有ワークツリーで複数セッションが並行作業すると、
このフックは case/*/run_* をファイル状態だけで見るため **他セッションの run** で
終了を block してしまう (実際に起きた)。そこで PreToolUse(Bash) フックが claim した
「このセッションが forge を回した case」だけを検査する (forge_gate_claims)。
session_id が取れない場合は絞り込まず**従来どおり全件検査**する (ゲートを弱めない)。
"""
import sys, os, glob, json, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forge_gate_claims as claims

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RECENT_SEC = 180 * 60

try:
    _ev = json.load(sys.stdin)
except Exception:
    _ev = {}
_claimed = claims.load(_ev.get("session_id"))   # None = 絞り込み不能 (全件検査)

bad = []
now = time.time()
for rh in glob.glob(os.path.join(ROOT, "case", "*", "run_*", "residual_history.csv")):
    try:
        m = os.path.getmtime(rh)
    except OSError:
        continue
    if now - m > RECENT_SEC:
        continue  # 旧 run は対象外 (このセッションの活動のみ強制)
    if now - m < 120:
        continue  # forge 実行中 (residual を書き込み中) はスキップ。完了時に wrapper が VERDICT を書く
    rel = os.path.relpath(os.path.dirname(rh), ROOT)
    if _claimed is not None:
        case = rel.split(os.sep)[1] if rel.startswith("case" + os.sep) else ""
        if case not in _claimed:
            continue  # 他セッションが回した run — このセッションは検査しない
    v = os.path.join(os.path.dirname(rh), "CONVERGENCE_VERDICT.txt")
    if (not os.path.exists(v)) or (os.path.getmtime(v) < m):
        bad.append(rel)

if bad:
    reason = (
        "forge を実行したのに収束チェック未実施/古い run があります:\n  - "
        + "\n  - ".join(sorted(bad))
        + "\n\n各 run を `solver_density_cuda/tools/run_case.sh <run_dir>` 経由で実行する"
        "(自動で VERDICT を残す)か、`python3 solver_density_cuda/tools/check_convergence.py "
        "<run_dir> > <run_dir>/CONVERGENCE_VERDICT.txt` を実行して VERDICT を残すこと。\n"
        "**収束/一致/ショックレス等を主張する応答では check_convergence.py の VERDICT 行を必ず引用する。"
        "衝撃波位置の安定や rms_ro 単独で『収束』と判断しない。**")
    print(json.dumps({"decision": "block", "reason": reason}))
else:
    print("{}")
