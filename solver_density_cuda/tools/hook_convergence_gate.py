#!/usr/bin/env python3
"""Stop フック: forge を回したのに**計算の後の判定が無い・古い・確定していない** run があればターン終了を block する。

case/*/run_*/residual_history.csv のうち「最近 (既定 180 分) 更新された」ものを対象に、同じ run の判定 `RUN_VERDICT.json`
(`judge_run.py`、plan tooling-convergence-and-wall-resolution-gates §4.6-5・§5.1 #12) を読み、次のどれかなら block する:

- 判定が無い (`RUN_VERDICT.json` が無い・読めない)。
- 判定が今の入力に対応しない (`judge_run.check_fresh` が STALE: config・残差の履歴・manifest・停止の記録・出力の一覧・仕様などが
  判定の後に変わった)。
- 段階起動の段の記録 (`state: deferred`) のまま総合が確定していない。ただし親 (設計チェーンの runner) が実行中
  (`JUDGE_PENDING.json` の pid がこのホストで生きている) なら確定待ちとして見逃す。

**FAIL や判定不能 (INDETERMINATE・判定器の例外 `state: error`) を正しく報告して終えるのは許す** (PASS になるまで止める設計にしない;
§4.6-5)。block の理由には、確定のさせ方 (`judge_run.py <run>`) と報告の規則 (`RUN_VERDICT.txt` の総合と各判定の行を貼る) を書く。
これは forge の実行方法 (直接/ラッパー/バックグラウンド) に依らずファイル状態だけで判定する。

**セッション別スコープ (2026-09-20)**: 共有ワークツリーで複数セッションが並行作業すると、このフックは case/*/run_* を
ファイル状態だけで見るため **他セッションの run** で終了を block してしまう (実際に起きた)。そこで PreToolUse(Bash) フックが
claim した「このセッションが forge を回した case」だけを検査する (forge_gate_claims)。session_id が取れない場合は絞り込まず
**従来どおり全件検査**する (ゲートを弱めない)。

試験用: 環境変数 `FORGE_GATE_ROOT` で検査するリポジトリの根を差し替えられる (既定はこのファイルから 2 つ上)。
"""
import glob
import json
import os
import socket
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import forge_gate_claims as claims  # noqa: E402

ROOT = os.path.abspath(os.environ.get("FORGE_GATE_ROOT") or os.path.join(HERE, "..", ".."))
RECENT_SEC = 180 * 60
ACTIVE_SEC = 120
PENDING = "JUDGE_PENDING.json"           # 設計チェーンの runner が段階起動の間だけ置く (design/forge_design/evaluate/judge.py)


def _pending_alive(run):
    """JUDGE_PENDING.json の親がこのホストで生きているか。"""
    try:
        with open(os.path.join(run, PENDING), encoding="utf-8") as f:
            rec = json.load(f)
        if rec.get("host") != socket.gethostname():
            return False
        pid = int(rec["pid"])
        os.kill(pid, 0)
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run_problem(run):
    """run の判定の問題 → 理由の文字列 (問題が無ければ None)。"""
    import judge_run as JR
    p = os.path.join(run, JR.VERDICT_JSON)
    if not os.path.exists(p):
        return "計算の後の判定が無い (RUN_VERDICT.json が無い)"
    try:
        with open(p, encoding="utf-8") as f:
            v = json.load(f)
    except (OSError, ValueError) as e:
        return "RUN_VERDICT.json を読めない (%s)" % e
    state = v.get("state")
    if state == "error":
        return None                       # 判定器の例外 = 判定不能 (INDETERMINATE) として報告して終えてよい
    if state == "deferred":
        if _pending_alive(run):
            return None                   # 親 (runner) が実行中: 確定待ち
        return ("段階起動の段の記録 (state deferred、総合 %s) のまま総合の判定が確定していない "
                "(親が確定せずに終わった)" % v.get("overall"))
    if state != "final":
        return "RUN_VERDICT.json の state %r が既知でない" % (state,)
    rc, msg = JR.check_fresh(run)
    if rc == 0:
        return None
    return msg


def main():
    try:
        ev = json.load(sys.stdin)
    except Exception:                    # noqa: BLE001
        ev = {}
    claimed = claims.load(ev.get("session_id"))   # None = 絞り込み不能 (全件検査)
    bad = []
    now = time.time()
    for rh in glob.glob(os.path.join(ROOT, "case", "*", "run_*", "residual_history.csv")):
        try:
            m = os.path.getmtime(rh)
        except OSError:
            continue
        if now - m > RECENT_SEC:
            continue  # 旧 run は対象外 (このセッションの活動のみ強制)
        if now - m < ACTIVE_SEC:
            continue  # forge 実行中 (residual を書き込み中) はスキップ。完了時に run_case.sh・親が判定を書く
        run = os.path.dirname(rh)
        rel = os.path.relpath(run, ROOT)
        if claimed is not None:
            case = rel.split(os.sep)[1] if rel.startswith("case" + os.sep) else ""
            if case not in claimed:
                continue  # 他セッションが回した run — このセッションは検査しない
        try:
            why = run_problem(run)
        except Exception as e:           # noqa: BLE001 — フック自体は落とさない (判定できないことを理由にして止める)
            why = "判定の状態を確かめられない (%s: %s)" % (type(e).__name__, e)
        if why:
            bad.append("%s: %s" % (rel, why))
    if bad:
        reason = (
            "forge を実行した run の計算の後の判定が無い・古い・確定していない:\n  - "
            + "\n  - ".join(sorted(bad))
            + "\n\n`python3 solver_density_cuda/tools/judge_run.py <run_dir>` を回して RUN_VERDICT.{json,txt} を確定する "
            "(run_case.sh は forge の後に自動で回す。段階起動は親が最後に 1 度確定する)。\n"
            "**FAIL・判定不能のまま正しく報告して終えてよい** (PASS は要らない)。「収束した」「定常」と書く応答には "
            "RUN_VERDICT.txt の総合 (OVERALL) と各判定の行を貼ること。")
        print(json.dumps({"decision": "block", "reason": reason}, ensure_ascii=False))
    else:
        print("{}")


if __name__ == "__main__":
    main()
