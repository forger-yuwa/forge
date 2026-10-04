#!/usr/bin/env python3
"""フック間で共有する「このセッションが forge を回した case」の台帳。

共有ワークツリーで複数セッションが並行作業すると、Stop フック
(hook_convergence_gate.py) が **他セッションの run** でターン終了を block してしまう。
フックは case/*/run_* をファイル状態だけで見ており、誰が回したかを区別できないため。

そこで PreToolUse(Bash) フック (hook_forge_guard.py) が「forge を実行する
Bash コマンド」を見たときに、その case ディレクトリを session_id 単位で claim し、
Stop フックは **claim した case の run だけ**を検査する。

- 台帳: ~/.cache/forge-gate-claims/<session_id>.txt (1 行 1 case ディレクトリ名)
- session_id が取れない場合は claim 不能 → Stop 側は **従来どおり全件検査** (ゲートを弱めない)
- 台帳が空/不在 = このセッションは forge を回していない → 検査対象なし
"""
import os, re, time

CLAIM_DIR = os.path.join(os.path.expanduser("~"), ".cache", "forge-gate-claims")
MAX_AGE_SEC = 7 * 24 * 3600
CASE_RE = r"(\d+\.[A-Za-z0-9_.\-]+)"


def _path(session_id):
    if not session_id or not re.fullmatch(r"[A-Za-z0-9_.\-]{1,128}", session_id):
        return None
    return os.path.join(CLAIM_DIR, session_id + ".txt")


def _prune():
    now = time.time()
    try:
        for n in os.listdir(CLAIM_DIR):
            p = os.path.join(CLAIM_DIR, n)
            if now - os.path.getmtime(p) > MAX_AGE_SEC:
                os.remove(p)
    except OSError:
        pass


def claim(session_id, cases):
    """cases (case ディレクトリ名の集合) をこのセッションの台帳に追記する。"""
    p = _path(session_id)
    if not p or not cases:
        return
    try:
        os.makedirs(CLAIM_DIR, exist_ok=True)
        _prune()
        have = load(session_id) or set()
        new = set(cases) - have
        if not new:
            return
        with open(p, "a", encoding="utf-8") as f:
            for c in sorted(new):
                f.write(c + "\n")
    except OSError:
        pass


def load(session_id):
    """claim 集合を返す。session_id 不明なら None (= 絞り込み不能)。"""
    p = _path(session_id)
    if not p:
        return None
    try:
        with open(p, encoding="utf-8") as f:
            return {ln.strip() for ln in f if ln.strip()}
    except FileNotFoundError:
        return set()
    except OSError:
        return None


def cases_in(cmd, cwd):
    """Bash コマンドと cwd から case ディレクトリ名を拾う。"""
    out = set(re.findall(r"case/" + CASE_RE, cmd or ""))
    m = re.search(r"/case/" + CASE_RE, cwd or "")
    if m:
        out.add(m.group(1))
    return out


def runs_forge(cmd):
    """この Bash コマンドは forge を回すか。

    run_case.sh / build/forge の直接指定のほか、case ディレクトリ内のシェルスクリプト
    (_t3_chain.sh 等が内部で run_case.sh を呼ぶ) も対象に含める。
    """
    cmd = cmd or ""
    if re.search(r"run_case\.sh|build/forge\b", cmd):
        return True
    return re.search(r"(^|[;&|]\s*|\bnohup\s+)(\./|bash\s+|sh\s+)\S*\.sh\b", cmd) is not None
