#!/usr/bin/env python3
"""codex_review.py の判定役の切り替え・結果の分類の試験 (plan tooling-review-fallback-backend §6 V2 (a)〜(j)
と、既存の記録の後方互換)。**実モデルは呼ばない**。
result 段レビュー (2026-10-14) の指摘 1〜3 の負例・回帰も持つ: (o) 判定は明示の判定の行から・指摘数は最後の行、
(p) 証拠の欄 (「証拠: 足りる」/「証拠不足:」) が必須、(q) 証拠一式の必須の git の取得の失敗 → 子を起動せず failed。

仕組み:
  * 偽の `codex`・`claude` を一時ディレクトリに作って PATH の先頭に置く (本物の ~/.local/bin は PATH から外す)。
    偽の CLI は呼ばれた引数・cwd・stdin・プロンプトを calls.jsonl に書き、環境変数 FAKE_CODEX_MODE / FAKE_CLAUDE_MODE
    で決めた振る舞い (成功・NO-GO・rate limit・usage 上限・overloaded・認証・timeout・空・誤りの文の本文・証拠不足) をする。
  * スイッチ ~/.config/forge/diagnose-backend は一時の HOME で制御する。
  * codex_review.py は小さな一時の git リポジトリ (fixture) に写して走らせる。ROOT はスクリプトの場所で決まるので、
    本物の notes/reviews/ を汚さない。出力先は一時ディレクトリ (--out-dir)。
  * 後方互換は本物のリポジトリの notes/reviews/ の記録 (codex_review.py が書いた形式のもの) から判定・指摘数・結論が
    取れることを確かめる (証拠の欄の契約より前の記録なので、分類は incomplete になってよい)。
  * 証拠一式の git の timeout は、偽の git を PATH の先頭に置き、モジュールを直に読んで打ち切りを 1 s にして試す。

使い方:
  python3 solver_density_cuda/tools/test_codex_review_backend.py        # 全部 ok なら 0、失敗があれば 1
  python3 solver_density_cuda/tools/test_codex_review_backend.py -v     # 失敗した試験の stdout も出す
"""
import argparse, contextlib, glob, importlib.util, io, json, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
TOOL = os.path.join(HERE, "codex_review.py")
CHECK_PLANS = os.path.join(HERE, "check_plans.py")
VERBOSE = "-v" in sys.argv

FAKE_CODEX = r'''#!@PY@
import json, os, stat, sys, time
argv = sys.argv[1:]
if argv[:1] == ["--version"]:
    print("codex-cli 0.0.0-fake"); sys.exit(0)
def devnull_stdin():
    st, dn = os.fstat(0), os.stat("/dev/null")
    return stat.S_ISCHR(st.st_mode) and st.st_rdev == dn.st_rdev
mode = os.environ.get("FAKE_CODEX_MODE", "ok_plan")
with open(os.environ["FAKE_CALLS"], "a", encoding="utf-8") as f:
    f.write(json.dumps({"cli": "codex", "argv": argv, "cwd": os.getcwd(), "stdin_devnull": devnull_stdin(),
                        "mode": mode}, ensure_ascii=False) + "\n")
out = argv[argv.index("-o") + 1]
print("Reading additional input from stdin...")
print("OpenAI Codex v0.0.0-fake\n--------\nworkdir: %s\nmodel: fake-codex-model\nprovider: openai\n--------" % os.getcwd())
print("user\n<prompt>\n")
sys.stdout.flush()
BODY = {
    "ok_plan": "**判定: GO-with-changes**\n要約。\n\n1. Major — 何か (根拠 a.py:1)\n\n証拠: 足りる\n指摘数: Critical 0 / Major 1 / Minor 2",
    "ok_diag": "結論: 次は A/B を 1 つ回す\n第 1 仮説: x  確度: 中\n不足情報: 無し\n証拠: 足りる",
    "nogo": "**判定: NO-GO**\n直すまで進めない。\n\n証拠: 足りる\n指摘数: Critical 1 / Major 0 / Minor 0",
    # result 段レビュー 2026-10-14 指摘 1 の入力 (+ 証拠の欄): 本文の前の方の「GO」でなく判定の行の NO-GO を取る
    "nogo_prose": "GO は不適切です。\n判定: NO-GO\n証拠: 足りる\n指摘数: Critical 1 / Major 0 / Minor 0",
    "verdict_conflict": "**判定: GO**\n要約。\n\n判定: NO-GO\n証拠: 足りる\n指摘数: Critical 0 / Major 0 / Minor 0",
    "no_evidence": "**判定: GO**\n要約。\n\n指摘数: Critical 0 / Major 0 / Minor 0",
    # 同 指摘 2 の入力: diagnose の書式どおりの「結論: 判断できません」+「不足情報:」、証拠の欄なし
    "diag_lack": "結論: 判断できません\n不足情報: 最終 HDF5 の NaN 検査結果",
    "errbody": "Error: rate limit reached",
    "readfail": "レビューに必要なファイルを読み取れません。",
    "insufficient": "**判定: GO**\n証拠不足: run_0001 の residual_history.csv\n指摘数: Critical 0 / Major 0 / Minor 0",
    # 本文に「rate limit」「ERROR:」を含むだけの正常なレビュー (unavailable にしてはいけない)
    "ok_mentions_ratelimit": "**判定: GO**\nERROR: rate limit の扱いを確認した。429 も見た。\n\n証拠: 足りる\n指摘数: Critical 0 / Major 0 / Minor 0",
}
if mode in BODY:
    if mode == "ok_mentions_ratelimit":
        print("ERROR: rate limit reached (これは codex が読んだファイルの中身の行)")
    open(out, "w", encoding="utf-8").write(BODY[mode]); print("tokens used\n123"); sys.exit(0)
if mode == "empty":
    open(out, "w", encoding="utf-8").write(""); sys.exit(0)
if mode == "rate_limit":
    print("stream error: unexpected status 429 Too Many Requests; retrying 1/5")
    print("ERROR: exceeded retry limit, last status: 429 Too Many Requests"); sys.exit(1)
if mode == "usage_limit":
    print("ERROR: You've hit your usage limit. Visit https://chatgpt.com/codex/settings/usage to purchase more credits "
          "or try again at 8:40 PM."); sys.exit(1)
if mode == "overloaded":
    print("ERROR: unexpected status 529: Overloaded"); sys.exit(1)
if mode == "auth":
    print("ERROR: unexpected status 401 Unauthorized: Not logged in. Run `codex login`."); sys.exit(1)
if mode == "timeout":
    time.sleep(60); sys.exit(0)
if mode == "failed":   # 本文はあるが rc 2・unavailable の型に当たらない
    open(out, "w", encoding="utf-8").write("**判定: GO**\n指摘数: Critical 0 / Major 0 / Minor 0")
    print("ERROR: sandbox setup failed"); sys.exit(2)
print("unknown mode " + mode); sys.exit(9)
'''

FAKE_CLAUDE = r'''#!@PY@
import json, os, re, stat, sys
argv = sys.argv[1:]
if argv[:1] in (["--version"], ["-v"]):
    print("0.0.0-fake (Claude Code)"); sys.exit(0)
def devnull_stdin():
    st, dn = os.fstat(0), os.stat("/dev/null")
    return stat.S_ISCHR(st.st_mode) and st.st_rdev == dn.st_rdev
mode = os.environ.get("FAKE_CLAUDE_MODE", "ok_plan")
prompt = argv[-1]
m = re.search(r"`([^`]+\.claude\.prompt\.md)`", prompt)
pfile = None
if m:
    p = m.group(1) if os.path.isabs(m.group(1)) else os.path.join(os.getcwd(), m.group(1))
    pfile = {"path": p, "exists": os.path.isfile(p), "text": open(p, encoding="utf-8").read() if os.path.isfile(p) else ""}
with open(os.environ["FAKE_CALLS"], "a", encoding="utf-8") as f:
    f.write(json.dumps({"cli": "claude", "argv": argv, "cwd": os.getcwd(), "stdin_devnull": devnull_stdin(),
                        "mode": mode, "prompt_file": pfile}, ensure_ascii=False) + "\n")
model = argv[argv.index("--model") + 1] if "--model" in argv else "?"
res = {"type": "result", "subtype": "success", "is_error": False, "num_turns": 3, "total_cost_usd": 0.0123,
       "modelUsage": {"claude-fake-%s-9" % model: {"outputTokens": 100}, "claude-fake-haiku-9": {"outputTokens": 5}},
       "permission_denials": [], "api_error_status": None}
BODY = {
    "ok_plan": "**判定: GO**\n要約。\n\n証拠: 足りる\n指摘数: Critical 0 / Major 0 / Minor 1",
    "ok_diag": "結論: 根拠の行を確かめた\n第 1 仮説: y  確度: 高\n証拠: 足りる",
    "errbody": "Error: rate limit reached",
    "readfail": "必要なファイルを読み取れません",
    "insufficient": "結論: 判断できない\n証拠不足: git diff の全文 (証拠一式に無い)",
    "insufficient_none": "結論: 判断した\n証拠: 足りる\n証拠不足: 無し",
    "insufficient_none2": "結論: 判断した\n証拠不足: 該当なし",
    "diag_lack": "結論: 判断できません\n不足情報: 最終 HDF5 の NaN 検査結果",
    "pretty": "結論: 整形された JSON でも読める\n証拠: 足りる",
    "denial": "結論: 書き込みは拒否された\n証拠: 足りる",
    "empty": "",
}
rc = 0
if mode in BODY:
    res["result"] = BODY[mode]
    if mode == "denial":
        res["permission_denials"] = [{"tool_name": "Bash", "tool_use_id": "x", "tool_input": {"command": "touch a"}}]
elif mode == "unavail_429":
    res.update(is_error=True, subtype="error_during_execution", api_error_status=429,
               result="API Error: 429 {\"type\":\"error\",\"error\":{\"type\":\"rate_limit_error\"}}")
    rc = 1
elif mode == "unavail_text":   # status は無いが本文が usage 上限
    res.update(is_error=True, result="Claude AI usage limit reached|1760000000")
    rc = 1
elif mode == "is_error_other":
    res.update(is_error=True, subtype="error_max_turns", result="max turns")
    rc = 1
elif mode == "notjson":
    print("panic: something broke", file=sys.stderr); sys.exit(1)
else:
    print("unknown mode " + mode, file=sys.stderr); sys.exit(9)
print(json.dumps(res, ensure_ascii=False, indent=2 if mode == "pretty" else None))
sys.exit(rc)
'''

# 偽の git (証拠一式の取得の失敗の試験用)。FAKE_GIT_MODE で result 段の必須の取得 (`<base>...HEAD` / `<base>..HEAD`) だけを
# 止める (sleep: 打ち切りまで眠る、fail_log / fail_diff: その取得だけ rc 1)。それ以外は本物の git に渡す。
FAKE_GIT = r'''#!@PY@
import os, sys, time
args = sys.argv[1:]
mode = os.environ.get("FAKE_GIT_MODE", "")
ranged = any(".." in x and x.endswith("HEAD") for x in args)
if ranged:
    with open(os.environ["FAKE_GIT_LOG"], "a", encoding="utf-8") as f:
        f.write(" ".join(args) + "\n")
    if mode == "sleep":
        time.sleep(60)
    if (mode == "fail_log" and args[0] == "log") or (mode == "fail_diff" and args[0] == "diff" and "--stat=200" not in args):
        print("fatal: fake failure (" + mode + ")", file=sys.stderr); sys.exit(1)
os.execv("@GIT@", ["git"] + args)
'''

PLAN_TEXT = """# demo plan

## メタ

- **status**: `draft`

## 4. 設計方針

なにか。

## 5. 実装ステップ

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 実装 | x | O |

## 6. 検証

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 9. 変更ログ

- 起票。
"""


class Fixture:
    def __init__(self):
        self.tmp = tempfile.mkdtemp(prefix="codex_review_test_")
        self.repo = os.path.join(self.tmp, "repo")
        tools = os.path.join(self.repo, "solver_density_cuda", "tools")
        os.makedirs(tools)
        shutil.copy(TOOL, tools)
        shutil.copy(CHECK_PLANS, tools)
        self.tool = os.path.join(tools, "codex_review.py")
        self.check_plans = os.path.join(tools, "check_plans.py")
        os.makedirs(os.path.join(self.repo, "plans", "active"))
        os.makedirs(os.path.join(self.repo, "notes", "reviews", "briefs"))
        self.plan = os.path.join(self.repo, "plans", "active", "tooling-demo.md")
        open(self.plan, "w", encoding="utf-8").write(PLAN_TEXT)
        open(os.path.join(self.repo, "plans", "README.md"), "w", encoding="utf-8").write("- tooling-demo.md\n")
        self.brief = os.path.join(self.repo, "notes", "reviews", "briefs", "2026-01-01-demo.md")
        open(self.brief, "w", encoding="utf-8").write("## 観測事実\n\nx\n")
        self.big = os.path.join(self.repo, "notes", "big.txt")
        open(self.big, "w", encoding="utf-8").write(("0123456789abcdef" * 63 + "\n") * 110)   # 約 110 KiB
        self.g("init", "-q", "-b", "main")
        self.g("add", "-A")
        self.g("commit", "-q", "-m", "base")
        self.g("checkout", "-q", "-b", "feature")
        os.makedirs(os.path.join(self.repo, "design"))
        open(os.path.join(self.repo, "design", "foo.py"), "w").write("x = 1\n" * 50)
        open(os.path.join(self.repo, "notes", "other.txt"), "w").write("y\n" * 400)
        self.g("add", "-A")
        self.g("commit", "-q", "-m", "feature change")
        # 偽の CLI
        py = sys.executable
        self.bin_all = os.path.join(self.tmp, "bin_all")
        self.bin_noc = os.path.join(self.tmp, "bin_noc")   # codex が無い
        for d, names in ((self.bin_all, ("codex", "claude")), (self.bin_noc, ("claude",))):
            os.makedirs(d)
            for n in names:
                p = os.path.join(d, n)
                open(p, "w").write((FAKE_CODEX if n == "codex" else FAKE_CLAUDE).replace("@PY@", py))
                os.chmod(p, 0o755)
        # 偽の git (本物の git の絶対パスへ渡す)
        self.bin_git = os.path.join(self.tmp, "bin_git")
        os.makedirs(self.bin_git)
        p = os.path.join(self.bin_git, "git")
        open(p, "w").write(FAKE_GIT.replace("@PY@", py).replace("@GIT@", shutil.which("git")))
        os.chmod(p, 0o755)
        self.home = os.path.join(self.tmp, "home")
        os.makedirs(os.path.join(self.home, ".config", "forge"))
        self.calls = os.path.join(self.tmp, "calls.jsonl")
        self.n = 0
        sys_path = ":".join(d for d in ("/usr/local/bin", "/usr/bin", "/bin") if os.path.isdir(d))
        self.sys_path = sys_path
        for n in ("codex", "claude"):
            if shutil.which(n, path=sys_path):
                raise SystemExit(f"試験の前提が崩れている: {n} が {sys_path} にある (偽の CLI だけを見せられない)")

    def g(self, *args):
        subprocess.run(["git", *args], cwd=self.repo, check=True, capture_output=True,
                       env={**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
                            "GIT_COMMITTER_EMAIL": "t@t"})

    def run(self, args, codex="ok_plan", claude="ok_plan", switch=None, no_codex=False, out_dir=None, timeout=60):
        self.n += 1
        if out_dir is None:
            out_dir = os.path.join(self.tmp, f"out_{self.n:02d}")
        sw = os.path.join(self.home, ".config", "forge", "diagnose-backend")
        if switch is None:
            if os.path.exists(sw):
                os.remove(sw)
        else:
            open(sw, "w").write(switch + "\n")
        if os.path.exists(self.calls):
            os.remove(self.calls)
        env = {k: v for k, v in os.environ.items() if not k.startswith("FAKE_")}
        env.update(PATH=(self.bin_noc if no_codex else self.bin_all) + ":" + self.sys_path, HOME=self.home,
                   FAKE_CALLS=self.calls, FAKE_CODEX_MODE=codex, FAKE_CLAUDE_MODE=claude)
        cmd = [sys.executable, self.tool, *args]
        if out_dir != "default":   # "default" は --out-dir を付けない (既定の notes/reviews/)
            cmd += ["--out-dir", out_dir]
        p = subprocess.run(cmd, cwd=self.repo, env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                           timeout=timeout)
        calls = [json.loads(x) for x in open(self.calls, encoding="utf-8")] if os.path.exists(self.calls) else []
        real_out = os.path.join(self.repo, "notes", "reviews") if out_dir == "default" else out_dir
        mds = sorted(glob.glob(os.path.join(real_out, "*.md")))
        mds = [m for m in mds if not m.endswith(".prompt.md")]
        md = mds[-1] if mds else ""
        return Result(p.returncode, p.stdout + p.stderr, calls, real_out, md)


class Result:
    def __init__(self, rc, out, calls, out_dir, md):
        self.rc, self.out, self.calls, self.out_dir, self.md = rc, out, calls, out_dir, md
        self.text = open(md, encoding="utf-8").read() if md else ""

    def n(self, cli):
        return sum(1 for c in self.calls if c["cli"] == cli)

    def call(self, cli):
        return next(c for c in self.calls if c["cli"] == cli)

    def has_row(self, stage="plan"):
        return bool(re.search(rf"^\| {stage} \| \d{{4}}-\d\d-\d\d \| \[", self.out, re.M))

    def files(self, pat):
        return sorted(glob.glob(os.path.join(self.out_dir, pat)))


RESULTS = []


def check(name, cond, r=None, note=""):
    RESULTS.append((bool(cond), name))
    print(("ok  " if cond else "NG  ") + name + (f"  ({note})" if note and not cond else ""))
    if not cond and VERBOSE and r is not None:
        print("    --- stdout ---\n    " + r.out.replace("\n", "\n    "))
        print("    --- md ---\n    " + r.text[:3000].replace("\n", "\n    "))


def load_module():
    spec = importlib.util.spec_from_file_location("codex_review_mod", TOOL)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    fx = Fixture()
    P = [fx.plan, "--stage", "plan", "--timeout", "20"]
    D = ["--stage", "diagnose", "--brief", fx.brief, "--timeout", "20"]

    # (a) codex が成功 → claude を呼ばない。記録の本文は codex の最終メッセージのまま
    r = fx.run(P, codex="ok_plan")
    check("(a) codex 成功: rc 0・claude 0 回・codex 1 回", r.rc == 0 and r.n("claude") == 0 and r.n("codex") == 1, r)
    check("(a) 記録: 題「codex レビュー」・分類 ok・backend codex・実効のモデル (ログの見出し)・判定",
          r.text.startswith("# codex レビュー: tooling-demo (plan)") and "- **分類**: `ok`" in r.text
          and "- **backend**: `codex`" in r.text and "`fake-codex-model`" in r.text
          and "- **判定**: **GO-with-changes**, 指摘 C0/M1/m2" in r.text and "- **CLI の版**: `codex-cli 0.0.0-fake`" in r.text
          and "依頼側との関係**: 不明" in r.text, r)
    check("(a) 本文は codex の最終メッセージそのまま",
          r.text.split("\n---\n", 1)[1].strip() == "**判定: GO-with-changes**\n要約。\n\n1. Major — 何か (根拠 a.py:1)\n\n"
          "証拠: 足りる\n指摘数: Critical 0 / Major 1 / Minor 2", r)
    cp = r.call("codex")["argv"][-1]
    check("(a) codex のプロンプトの出力形式: 判定の行・証拠の欄 (「証拠: 足りる」か「証拠不足: <要るもの>」)・指摘数",
          "`判定: <GO / GO-with-changes / NO-GO のどれか 1 語>`" in cp and "`証拠: 足りる`" in cp
          and "`証拠不足: <要るもの>`" in cp and "`指摘数: Critical N / Major N / Minor N`" in cp, r)
    check("(a) §6.1 の行を出す ([Claude 代替] は付かない)", r.has_row() and "[Claude 代替]" not in r.out, r)
    check("(a) codex のときは証拠一式を作らない (今の振る舞い)", not r.files("*.evidence"), r)
    check("(a) codex の起動作法は今のまま (read-only・stdin /dev/null・cwd=ROOT・effort)",
          r.call("codex")["argv"][:4] == ["exec", "--sandbox", "read-only", "-C"] and r.call("codex")["stdin_devnull"]
          and r.call("codex")["cwd"] == fx.repo and 'model_reasoning_effort="high"' in r.call("codex")["argv"], r)

    # (a') 本文に「rate limit」「ERROR:」があるだけの正常なレビュー → ok (unavailable にしない)
    r = fx.run(P, codex="ok_mentions_ratelimit")
    check("(a') 本文に rate limit と書いてある正常なレビューは ok・切り替えない",
          r.rc == 0 and r.n("claude") == 0 and "- **分類**: `ok`" in r.text, r)

    # (b) codex の unavailable の各型 → claude 1 回・理由が記録に入る
    for mode, kind, frag in (("rate_limit", "rate_limit", "429"), ("usage_limit", "usage_limit", "usage limit"),
                             ("overloaded", "overloaded", "529"), ("auth", "auth", "401"),
                             ("timeout", "timeout", "timeout"), ("empty", "empty", "最終メッセージが空"),
                             ("cli_missing", "cli_missing", "PATH に無い")):
        args = P[:-1] + ["2"] if mode == "timeout" else P
        r = fx.run(args, codex=mode if mode != "cli_missing" else "ok_plan", no_codex=(mode == "cli_missing"))
        reason = re.search(r"^- \*\*backend\*\*: `claude` — 選んだ理由: (.+)$", r.text, re.M)
        check(f"(b) codex {mode} → claude 1 回・rc 0・ok",
              r.rc == 0 and r.n("claude") == 1 and r.n("codex") == (0 if mode == "cli_missing" else 1)
              and "- **分類**: `ok`" in r.text, r)
        check(f"(b) codex {mode} → 理由「自動の切り替え (codex が unavailable/{kind}: …{frag}…)」が記録に入る",
              bool(reason) and f"自動の切り替え (codex が unavailable/{kind}" in reason.group(1) and frag in reason.group(1), r,
              note=reason.group(1) if reason else "理由の行が無い")
        # (j) 試行ごとのログが残る
        check(f"(j) codex {mode} → .codex.log と .claude.log が両方残り、記録に両方が載る",
              os.path.isfile(r.md[:-3] + ".codex.log") and os.path.isfile(r.md[:-3] + ".claude.log")
              and f"(codex: unavailable/{kind})" in r.text and "(claude: ok)" in r.text, r)

    # (c) 有効な NO-GO → 切り替えない・ok
    r = fx.run(P, codex="nogo")
    check("(c) codex NO-GO → claude 0 回・rc 0・ok・行に NO-GO",
          r.rc == 0 and r.n("claude") == 0 and "- **分類**: `ok`" in r.text and r.has_row()
          and re.search(r"^\| plan .*\| NO-GO, C1/M0/m0 \|", r.out, re.M), r)

    # (d) 誤りの文の本文 → ok にならない・非 0・§6.1 の行を出さない (codex・claude とも)
    for mode in ("errbody", "readfail"):
        r = fx.run(P, codex=mode)
        check(f"(d) codex 本文「{mode}」(rc 0) → incomplete・非 0・行なし・切り替えない",
              r.rc == 2 and not r.has_row() and r.n("claude") == 0 and r.text.startswith("**未成立 (incomplete)**"), r)
        r = fx.run(P + ["--backend", "claude"], claude=mode)
        check(f"(d) claude 本文「{mode}」 → incomplete・非 0・行なし",
              r.rc == 2 and not r.has_row() and r.text.startswith("**未成立 (incomplete)**"), r)
        r = fx.run(D + ["--backend", "claude"], claude=mode)
        check(f"(d) claude diagnose 本文「{mode}」 → incomplete・非 0・行なし",
              r.rc == 2 and not r.has_row("diagnose") and r.text.startswith("**未成立 (incomplete)**"), r)
    r = fx.run(D, codex="errbody")
    check("(d) codex diagnose 本文「Error: rate limit reached」 → incomplete・非 0・行なし",
          r.rc == 2 and not r.has_row("diagnose") and r.n("claude") == 0 and r.text.startswith("**未成立 (incomplete)**"), r)

    # (e) 証拠不足 → incomplete
    r = fx.run(P, codex="insufficient")
    check("(e) codex「証拠不足:」 → incomplete・切り替えない", r.rc == 2 and r.n("claude") == 0
          and r.text.startswith("**未成立 (incomplete)** — incomplete: 証拠不足:"), r)
    r = fx.run(D + ["--backend", "claude"], claude="insufficient")
    check("(e) claude「証拠不足:」 → incomplete・非 0・行なし",
          r.rc == 2 and not r.has_row("diagnose") and "証拠不足: git diff の全文" in r.text.splitlines()[0], r)
    r = fx.run(D + ["--backend", "claude"], claude="insufficient_none")
    check("(e) 「証拠: 足りる」があれば「証拠不足: 無し」は不足として数えない → ok", r.rc == 0 and "- **分類**: `ok`" in r.text, r)
    r = fx.run(D + ["--backend", "claude"], claude="insufficient_none2")
    check("(e) 「証拠不足: 該当なし」だけでは証拠の欄の代わりにならない → incomplete・非 0・行なし",
          r.rc == 2 and not r.has_row("diagnose") and "証拠の欄" in r.text.splitlines()[0], r)

    # (f) 96 KiB 超のプロンプト → ファイル渡し (codex・claude とも)
    r = fx.run(P + ["--extra", fx.big], codex="ok_plan")
    pf = r.files("*.prompt.md")
    arg = r.call("codex")["argv"][-1]
    check("(f) codex: 96 KiB 超 → .prompt.md に書いて引数はその案内だけ",
          r.rc == 0 and len(pf) == 1 and os.path.getsize(pf[0]) > 96 * 1024 and len(arg.encode()) < 2048
          and os.path.basename(pf[0]) in arg, r)
    r = fx.run(P + ["--extra", fx.big, "--backend", "claude"])
    c = r.call("claude")
    check("(f) claude: 96 KiB 超 → .claude.prompt.md (前書きを含む) に書いて引数はその案内だけ",
          r.rc == 0 and c["prompt_file"] and c["prompt_file"]["exists"] and len(c["argv"][-1].encode()) < 2048
          and c["prompt_file"]["text"].startswith("## この依頼の前書き") and len(c["prompt_file"]["text"].encode()) > 96 * 1024, r)

    # (g) スイッチ claude → codex を呼ばない
    r = fx.run(P, switch="claude")
    check("(g) スイッチ claude → codex 0 回・claude 1 回・理由「スイッチ」",
          r.rc == 0 and r.n("codex") == 0 and r.n("claude") == 1 and "選んだ理由: スイッチ" in r.text, r)
    r = fx.run(P, switch="codex")
    check("(g') スイッチ codex (従来の値) → auto は codex を試す", r.rc == 0 and r.n("codex") == 1 and r.n("claude") == 0, r)

    # (h) --backend codex は失敗しても切り替えない
    for mode in ("rate_limit", "cli_missing"):
        r = fx.run(P + ["--backend", "codex"], codex=mode if mode != "cli_missing" else "ok_plan",
                   no_codex=(mode == "cli_missing"))
        check(f"(h) --backend codex・codex {mode} → claude 0 回・unavailable・rc 3・行なし・記録の先頭に未成立",
              r.rc == 3 and r.n("claude") == 0 and not r.has_row() and r.text.startswith("**未成立 (unavailable)**")
              and "選んだ理由: 明示 (`--backend codex`)" in r.text, r)
    r = fx.run(P, codex="failed")
    check("(h') codex が本文ありの rc 2 (unavailable の型に当たらない) → failed・切り替えない・rc 1",
          r.rc == 1 and r.n("claude") == 0 and r.text.startswith("**未成立 (failed)**"), r)

    # (i) result 段の証拠一式 (diff の stat) がプロンプトに入る
    r = fx.run([fx.plan, "--stage", "result", "--base", "main", "--backend", "claude"])
    stat = subprocess.run(["git", "diff", "--stat=200", "main...HEAD"], cwd=fx.repo, capture_output=True, text=True).stdout
    prompt = r.call("claude")["argv"][-1]
    ev = r.md[:-3] + ".evidence"
    check("(i) result 段: diff の stat がプロンプトに入り、証拠一式のファイルがある",
          r.rc == 0 and stat.strip() and stat.rstrip("\n") in prompt
          and all(os.path.isfile(os.path.join(ev, f)) for f in ("README.md", "head.txt", "diff_stat.txt", "diff.patch", "log.txt"))
          and "feature change" in open(os.path.join(ev, "log.txt")).read()
          and "design/foo.py" in open(os.path.join(ev, "diff.patch")).read(), r)
    check("(i) 前書き: 道具は Read・Grep・Glob だけ・証拠一式の場所・sed の読み替え・同調しない・証拠不足の行",
          prompt.startswith("## この依頼の前書き") and "道具は Read・Grep・Glob だけ" in prompt and ev in prompt
          and "sed -n 'A,Bp'" in prompt and "同調せず" in prompt and "`証拠不足: <要るもの>`" in prompt, r)
    r = fx.run([fx.plan, "--stage", "plan", "--backend", "claude"])
    ev = r.md[:-3] + ".evidence"
    check("(i') plan 段: 証拠一式は HEAD・ブランチ", os.path.isfile(os.path.join(ev, "head.txt"))
          and "branch: feature" in open(os.path.join(ev, "head.txt")).read() and not os.path.exists(os.path.join(ev, "diff.patch")), r)

    # 証拠一式の diff を切る (上限を小さくして関数を直に呼ぶ)
    mod = load_module()
    mod.ROOT = fx.repo
    evd = os.path.join(fx.tmp, "ev_trunc")
    e = mod.build_evidence(evd, "result", "main", diff_max_bytes=600)
    readme = open(os.path.join(evd, "README.md"), encoding="utf-8").read()
    patch = open(os.path.join(evd, "diff.patch"), encoding="utf-8").read()
    check("(i'') diff が上限を超えたら優先パス (design/) を入れ、入らなかったファイルを README に列挙",
          "design/foo.py" in patch and "notes/other.txt" not in patch and "省略: `notes/other.txt`" in readme
          and any("上限" in n for n in e["notes"]), note=readme[:400])
    check("(i'') 意図した差分の省略は取得の失敗 (errors) にしない", e["errors"] == [], note=str(e["errors"]))

    # claude の起動の契約 (plan §4.2)
    r = fx.run(P + ["--backend", "claude", "--effort", "xhigh", "--requester-model", "claude-fake-opus-9"])
    c = r.call("claude")
    a = c["argv"]
    want = ["-p", "--safe-mode", "--model", "opus", "--effort", "xhigh", "--tools", "Read,Grep,Glob",
            "--permission-mode", "default", "--no-session-persistence", "--strict-mcp-config"]
    check("(k) claude の引数: -p --safe-mode --model opus --effort (codex と同じ値) --tools Read,Grep,Glob --permission-mode default "
          "--no-session-persistence --strict-mcp-config --output-format json",
          a[:len(want)] == want and a[a.index("--output-format") + 1] == "json", r, note=str(a[:-1]))
    check("(k) Bash・Write・Edit を渡さない (プロンプト以外の引数に現れない)",
          not any(t in x for x in a[:-1] for t in ("Bash", "Write", "Edit")), r)
    check("(k) stdin は /dev/null・cwd はリポジトリのルート", c["stdin_devnull"] and c["cwd"] == fx.repo, r)
    check("(k) 記録の外 (--out-dir がリポジトリ外) なら --add-dir で証拠一式を読めるようにする",
          "--add-dir" in a and a[a.index("--add-dir") + 1] == r.out_dir, r)
    check("(k) 記録: 題「Claude 代替レビュー」・代替の明記・要求/実効のモデル・CLI の版・依頼側との関係 同一・permission_denials・行に [Claude 代替]",
          r.text.startswith("# Claude 代替レビュー: tooling-demo (plan)") and "**Claude の別インスタンスによる代替**" in r.text
          and "- **要求したモデル**: `opus`" in r.text and "- **実効のモデル**: `claude-fake-opus-9`" in r.text
          and "- **CLI の版**: `0.0.0-fake (Claude Code)`" in r.text and "- **依頼側との関係**: 同一" in r.text
          and "- **permission_denials**: 0 件" in r.text and "total_cost_usd 0.0123" in r.text
          and re.search(r"^\| plan \|.*\| \[Claude 代替\] GO, C0/M0/m1 \|", r.out, re.M), r)
    r = fx.run(P + ["--backend", "claude", "--requester-model", "gpt-x"])
    check("(k) 依頼側との関係: 別", "- **依頼側との関係**: 別 (依頼側 `gpt-x`" in r.text, r)
    r = fx.run(P + ["--backend", "claude"], out_dir="default")
    c = r.call("claude")
    check("(k) 既定の置き場 (リポジトリ内) なら --add-dir なし・証拠一式は相対パス notes/reviews/<名>.evidence",
          "--add-dir" not in c["argv"] and "`notes/reviews/" in c["argv"][-1] and ".evidence/`" in c["argv"][-1]
          and r.md.startswith(os.path.join(fx.repo, "notes", "reviews")), r)
    check("(k) 名前は codex と同じ規則 (YYYY-MM-DD-<stem>-<stage>.md)",
          re.fullmatch(r"\d{4}-\d\d-\d\d-tooling-demo-plan\.md", os.path.basename(r.md)) is not None, r)
    r = fx.run(D + ["--backend", "claude"], claude="denial")
    check("(k) permission_denials を記録に残す (道具名と入力)", "- **permission_denials**: 1 件 — Bash" in r.text
          and "touch a" in r.text, r)
    check("(k) diagnose (claude): 題「Claude 代替諮問」・結論・§6.1 の行に [Claude 代替]・痕跡の文",
          r.text.startswith("# Claude 代替諮問 (diagnose): demo") and "- **結論**: 書き込みは拒否された" in r.text
          and re.search(r"^\| diagnose \|.*\| \[Claude 代替\] 書き込みは拒否された \|", r.out, re.M)
          and "Claude の別インスタンス (同一系統の代替) に諮った" in r.out, r)

    # claude 自体が unavailable / failed
    r = fx.run(P + ["--backend", "claude"], claude="unavail_429")
    check("(l) claude is_error・api_error_status 429 → unavailable/rate_limit・rc 3・行なし",
          r.rc == 3 and not r.has_row() and r.text.startswith("**未成立 (unavailable)** — unavailable/rate_limit"), r)
    r = fx.run(P + ["--backend", "claude"], claude="unavail_text")
    check("(l) claude is_error・本文が usage 上限 → unavailable/usage_limit", r.rc == 3
          and "unavailable/usage_limit" in r.text.splitlines()[0], r)
    r = fx.run(P + ["--backend", "claude"], claude="is_error_other")
    check("(l) claude is_error (型に当たらない) → failed・rc 1", r.rc == 1 and r.text.startswith("**未成立 (failed)**"), r)
    r = fx.run(P + ["--backend", "claude"], claude="notjson")
    check("(l) claude の JSON が読めない → failed", r.rc == 1 and r.text.startswith("**未成立 (failed)**"), r)
    r = fx.run(D + ["--backend", "claude"], claude="pretty")
    check("(l) claude の JSON が複数行に整形されていても読める → ok", r.rc == 0 and "- **結論**: 整形された JSON でも読める" in r.text, r)
    r = fx.run(P + ["--backend", "claude"], claude="empty")
    check("(l) claude の result が空 → unavailable/empty", r.rc == 3 and "unavailable/empty" in r.text.splitlines()[0], r)
    r = fx.run(P, codex="rate_limit", claude="unavail_429")
    check("(l) auto: codex も claude も unavailable → rc 3・両方のログ・行なし",
          r.rc == 3 and r.n("codex") == 1 and r.n("claude") == 1 and not r.has_row()
          and "(codex: unavailable/rate_limit)" in r.text and "(claude: unavailable/rate_limit)" in r.text, r)
    r = fx.run(P, no_codex=True, claude="ok_plan")
    check("(l) codex が無く claude 成功の記録の本文は claude の result", r.text.split("\n---\n", 1)[1].strip().startswith("**判定: GO**"), r)

    # 表示された §6.1 の行で check_plans.py が通る (fixture の plan に貼って in_progress にする)
    r = fx.run(P + ["--backend", "claude"], out_dir="default")
    row = next((ln for ln in r.out.splitlines() if ln.startswith("| plan |")), "")
    text = open(fx.plan, encoding="utf-8").read().replace("- **status**: `draft`", "- **status**: `in_progress`")
    text = text.replace("| --- | --- | --- | --- | --- |\n", "| --- | --- | --- | --- | --- |\n" + row + "\n", 1)
    p2 = os.path.join(fx.tmp, "plancopy", "tooling-demo.md")
    os.makedirs(os.path.dirname(p2))
    open(p2, "w", encoding="utf-8").write(text)
    cp = subprocess.run([sys.executable, fx.check_plans, p2], capture_output=True, text=True)
    check("(m) 表示された §6.1 の行 ([Claude 代替]) を貼った plan で check_plans.py が PASS",
          bool(row) and cp.returncode == 0 and "VERDICT: PASS" in cp.stdout, r, note=cp.stdout[-300:])

    # (o) 判定は明示の判定の行から・指摘数は最後の行 (result 段レビュー 2026-10-14 指摘 1)
    mod = load_module()
    EV = "\n証拠: 足りる\n"
    CNT = "指摘数: Critical 0 / Major 0 / Minor 0"
    for name, body, want_v, want_c, want_cls in (
            ("指摘 1 の入力 (「GO は不適切です。」の後に「判定: NO-GO」、証拠の欄なし) → 判定は NO-GO (分類は証拠の欄なしで incomplete)",
             "GO は不適切です。\n判定: NO-GO\n指摘数: Critical 1 / Major 0 / Minor 0", "NO-GO", "C1/M0/m0", "incomplete"),
            ("指摘 1 の入力 + 証拠の欄 → ok・NO-GO",
             "GO は不適切です。\n判定: NO-GO" + EV + "指摘数: Critical 1 / Major 0 / Minor 0", "NO-GO", "C1/M0/m0", "ok"),
            ("判定の語だけの行 `**NO-GO**  ` (既存の記録 32 件の書式) → NO-GO", "**NO-GO**  \n要約。" + EV + CNT, "NO-GO", "C0/M0/m0", "ok"),
            ("`## 判定: **GO-with-changes**` → GO-with-changes", "## 判定: **GO-with-changes**\n要約。" + EV + CNT,
             "GO-with-changes", "C0/M0/m0", "ok"),
            ("全角コロン `**判定：GO**` → GO", "**判定：GO**\n要約。" + EV + CNT, "GO", "C0/M0/m0", "ok"),
            ("コードブロックの中の例文の判定の行・指摘数の行は拾わない (指摘数は最後の行)",
             "**判定: GO-with-changes**\n\n1. Major — 例:\n\n   ```text\n   GO は不適切です。\n   判定: NO-GO\n"
             "   指摘数: Critical 1 / Major 0 / Minor 0\n   ```\n" + EV + "指摘数: Critical 0 / Major 3 / Minor 2",
             "GO-with-changes", "C0/M3/m2", "ok"),
            ("判定の行が相反する (GO と NO-GO) → incomplete・判定 ?", "**判定: GO**\n\n判定: NO-GO" + EV + CNT, "?", "C0/M0/m0",
             "incomplete"),
            ("同じ判定の行が 2 回 → その判定", "**判定: NO-GO**\n\n判定: NO-GO" + EV + CNT, "NO-GO", "C0/M0/m0", "ok"),
            ("判定の行が無い (文中の GO・NO-GO だけ) → incomplete", "GO は不適切です。NO-GO とすべき。" + EV + CNT, "?", "C0/M0/m0",
             "incomplete"),
            ("箇条書きの「- 判定: GO」・後ろに文が続く「判定: GO は早い」は判定の行でない",
             "- 判定: GO の場合は進める\n判定: GO は早い\n**判定: NO-GO**" + EV + CNT, "NO-GO", "C0/M0/m0", "ok"),
            ("指摘数は最後の行 (本文中の例文 `指摘数: Critical 1 / …` を拾わない)",
             "**判定: GO**\n例文 `指摘数: Critical 1 / Major 0 / Minor 0` を拾っていた。" + EV + "指摘数: Critical 0 / Major 3 / Minor 2",
             "GO", "C0/M3/m2", "ok")):
        cls, why, info = mod.classify_body("result", body)
        check(f"(o) {name}", (cls, info.get("verdict"), info.get("counts")) == (want_cls, want_v, want_c),
              note=f"{cls} {why} {info}")
    # 回帰: 指摘 1〜3 を出したレビュー自身の記録。旧来の抽出は見出しを「C1/M0/m0」と誤った (本文の例文を拾った)
    rec = os.path.join(REAL_ROOT, "notes", "reviews", "2026-10-14-tooling-review-fallback-backend-result.md")
    if os.path.isfile(rec):
        cls, why, info = mod.classify_body("result", open(rec, encoding="utf-8").read().split("\n---\n", 1)[1])
        check("(o) 回帰: 2026-10-14-tooling-review-fallback-backend-result.md の本文 → GO-with-changes・C0/M3/m2 "
              "(コードブロックの「判定: NO-GO」・例文の指摘数を拾わない。証拠の欄の契約より前なので分類は incomplete)",
              (info.get("verdict"), info.get("counts"), cls) == ("GO-with-changes", "C0/M3/m2", "incomplete")
              and "証拠の欄" in why, note=f"{cls} {why} {info}")
    else:
        check("(o) 回帰の入力 2026-10-14-tooling-review-fallback-backend-result.md がある", False, note=rec)
    r = fx.run(P, codex="nogo_prose")
    check("(o) codex 本文「GO は不適切です。/判定: NO-GO」 → ok・rc 0・見出しと §6.1 の行が NO-GO, C1/M0/m0",
          r.rc == 0 and "- **判定**: **NO-GO**, 指摘 C1/M0/m0" in r.text
          and re.search(r"^\| plan .*\| NO-GO, C1/M0/m0 \|", r.out, re.M), r)
    r = fx.run(P, codex="verdict_conflict")
    check("(o) codex 本文の判定の行が相反 → incomplete・rc 2・行なし・切り替えない",
          r.rc == 2 and not r.has_row() and r.n("claude") == 0
          and r.text.startswith("**未成立 (incomplete)** — incomplete: 判定の行が相反する: GO / NO-GO"), r)

    # (p) 証拠の欄が必須 (同 指摘 2)。diagnose の「不足情報:」はそれだけでは incomplete にしない
    for name, body, want in (
            ("「不足情報:」+「証拠: 足りる」 → ok (不足情報は追加の希望の欄)", "結論: x\n不足情報: y\n証拠: 足りる", "ok"),
            ("「証拠: 足りない (…)」は証拠の欄として数えない → incomplete", "結論: x\n証拠: 足りない (run の VERDICT)", "incomplete"),
            ("書式の写し「証拠: 足りる | 証拠不足: <要るもの>」は数えない → incomplete",
             "結論: x\n証拠: 足りる | 証拠不足: <要るもの>", "incomplete"),
            ("「`証拠: 足りる`」(バッククォート)・「**証拠: 足りる (diff を確認)**」は数える → ok",
             "結論: x\n`証拠: 足りる`\n**証拠: 足りる (diff を確認)**", "ok"),
            ("「証拠: 足りる」と「証拠不足: <要るもの>」が両方 → 不足を優先して incomplete",
             "結論: x\n証拠: 足りる\n証拠不足: run_0001 の CONVERGENCE_VERDICT.txt", "incomplete")):
        cls, why, info = mod.classify_body("diagnose", body)
        check(f"(p) {name}", cls == want, note=f"{cls} {why}")
    r = fx.run(D, codex="diag_lack")
    check("(p) codex diagnose 指摘 2 の入力 (「結論: 判断できません」+「不足情報: …」、証拠の欄なし) → incomplete・rc 2・行なし・切り替えない",
          r.rc == 2 and not r.has_row("diagnose") and r.n("claude") == 0
          and r.text.startswith("**未成立 (incomplete)**") and "証拠の欄" in r.text.splitlines()[0], r)
    r = fx.run(D + ["--backend", "claude"], claude="diag_lack")
    check("(p) claude diagnose 指摘 2 の入力 → incomplete・rc 2・行なし",
          r.rc == 2 and not r.has_row("diagnose") and "証拠の欄" in r.text.splitlines()[0], r)
    r = fx.run(P, codex="no_evidence")
    check("(p) codex plan 判定・指摘数はあるが証拠の欄なし → incomplete・rc 2・行なし",
          r.rc == 2 and not r.has_row() and "証拠の欄" in r.text.splitlines()[0], r)
    r = fx.run(D, codex="ok_diag")
    dp = r.call("codex")["argv"][-1]
    check("(p) diagnose のプロンプトの出力形式に証拠の欄と「不足情報:」との区別がある",
          r.rc == 0 and "証拠: 足りる | 証拠不足: <要るもの>" in dp and "`不足情報:` は結論を妨げない" in dp, r)
    r = fx.run([fx.plan, "--stage", "result", "--base", "main"], codex="ok_plan")
    rp = r.call("codex")["argv"][-1]
    check("(p) result 段 (codex) のプロンプトにも証拠の欄", r.rc == 0 and "`証拠: 足りる`" in rp and "`証拠不足: <要るもの>`" in rp
          and r.has_row("result"), r)
    r = fx.run([fx.plan, "--stage", "result", "--base", "main", "--backend", "claude"])
    check("(p) claude の前書き: 足りているときは「証拠: 足りる」を書く (どちらか一方は必須)",
          r.rc == 0 and "足りているときは `証拠: 足りる` の行を書く" in r.call("claude")["argv"][-1], r)

    # (q) 証拠一式の必須の git の取得の失敗 → 子を起動せず failed: 証拠一式 (同 指摘 3)
    r = fx.run([fx.plan, "--stage", "result", "--base", "no-such-ref", "--backend", "claude", "--timeout", "20"])
    log = open(r.md[:-3] + ".claude.log", encoding="utf-8").read() if r.md and os.path.isfile(r.md[:-3] + ".claude.log") else ""
    evr = os.path.join(r.md[:-3] + ".evidence", "README.md") if r.md else ""
    check("(q) 不正な --base (claude) → claude を起動しない・rc 1・記録の先頭「未成立 (failed: 証拠一式)」・行なし",
          r.rc == 1 and r.n("claude") == 0 and not r.has_row("result")
          and r.text.startswith("**未成立 (failed: 証拠一式)** — failed: 証拠一式: 必須の git の取得に失敗"), r)
    check("(q) 不正な --base → 失敗した git のコマンドと stderr が claude のログ・証拠一式の README・記録に残る",
          "`git diff --stat=200 no-such-ref...HEAD` が失敗 (rc=128)" in log and os.path.isfile(evr)
          and "取得の失敗" in open(evr, encoding="utf-8").read() and "no-such-ref" in r.text
          and "(claude: failed: 証拠一式)" in r.text, r, note=log[:300])
    r = fx.run([fx.plan, "--stage", "result", "--base", "no-such-ref", "--timeout", "20"], codex="rate_limit")
    check("(q) auto: codex unavailable → claude への切り替えでも不正な --base なら起動しない・rc 1・両方の試行が記録に残る",
          r.rc == 1 and r.n("codex") == 1 and r.n("claude") == 0 and r.text.startswith("**未成立 (failed: 証拠一式)**")
          and "(codex: unavailable/rate_limit)" in r.text and "(claude: failed: 証拠一式)" in r.text, r)

    # 偽の git を PATH の先頭に置き、モジュールを直に呼ぶ (打ち切りを 1 s にする)
    gitlog = os.path.join(fx.tmp, "fake_git.log")

    def in_process(git_mode, fn):
        if os.path.exists(gitlog):
            os.remove(gitlog)
        if os.path.exists(fx.calls):
            os.remove(fx.calls)
        old = dict(os.environ)
        os.environ.update(PATH=fx.bin_git + ":" + fx.bin_all + ":" + fx.sys_path, HOME=fx.home, FAKE_CALLS=fx.calls,
                          FAKE_GIT_MODE=git_mode, FAKE_GIT_LOG=gitlog, FAKE_CLAUDE_MODE="ok_plan")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                return fn()
        finally:
            os.environ.clear()
            os.environ.update(old)

    def git_calls():
        return open(gitlog, encoding="utf-8").read().splitlines() if os.path.exists(gitlog) else []

    mq = load_module()
    mq.ROOT = fx.repo
    mq.EVIDENCE_GIT_TIMEOUT = 1
    outq = os.path.join(fx.tmp, "q_out")
    os.makedirs(outq)
    ns = argparse.Namespace(claude_model="opus", stage="result", base="main", effort="high", timeout=20)
    at = in_process("sleep", lambda: mq.run_claude(ns, "prompt", "q-timeout", outq))
    gc = git_calls()
    claude_calls = open(fx.calls, encoding="utf-8").read().splitlines() if os.path.exists(fx.calls) else []
    check("(q) git の timeout (偽の git が眠る) → failed: 証拠一式・claude を起動しない・最初の取得で止める (timeout を重ねない)",
          at.cls == "failed" and at.phase == "証拠一式" and not at.ran and "timeout 1s" in at.detail and not claude_calls
          and len(gc) == 1 and gc[0].startswith("diff --stat=200 main...HEAD"), note=f"{at.cls} {at.detail} git={gc}")
    for mode, cmd in (("fail_log", "git log --oneline main..HEAD"), ("fail_diff", "git diff main...HEAD")):
        e = in_process(mode, lambda: mq.build_evidence(os.path.join(fx.tmp, "q_" + mode), "result", "main"))
        check(f"(q) 2 つ目以降の必須の取得 ({cmd}) の失敗も errors に入る",
              len(e["errors"]) == 1 and f"`{cmd}` が失敗 (rc=1)" in e["errors"][0] and "fake failure" in e["errors"][0],
              note=str(e["errors"]))
    e = in_process("", lambda: mq.build_evidence(os.path.join(fx.tmp, "q_ok"), "result", "main"))
    check("(q) 偽の git でも取得が成功すれば errors なし (stat・log・diff の 3 回)", e["errors"] == [] and len(git_calls()) == 3,
          note=f"{e['errors']} {git_calls()}")

    # (n) 後方互換: 既存の記録 (codex_review.py が書いた形式、証拠の欄の契約より前) から判定・指摘数・結論が取れる。
    #     分類は新しい契約では incomplete になってよい (証拠の欄が無い)。判定・指摘数は記録の見出しと一致することも見る。
    mod = load_module()
    head_rev = re.compile(r"^# (?:codex レビュー|Claude 代替レビュー): \S+ \((plan|result)\)$")
    head_diag = re.compile(r"^# (?:codex 諮問|Claude 代替諮問) \(diagnose\): \S+$")
    hdr_v = re.compile(r"^- \*\*判定\*\*: \*\*(\S+?)\*\*, 指摘 (C\S+)$", re.M)
    # 旧来の抽出 (本文で最初の一致) が見出しを誤った記録。見出しとの一致は見ない (正しい値は (o) の回帰で確かめる)
    known_bad_header = {"2026-10-14-tooling-review-fallback-backend-result.md"}
    n_all, n_recent, bad, mismatch, other_why = 0, 0, [], [], []
    for f in sorted(glob.glob(os.path.join(REAL_ROOT, "notes", "reviews", "*.md"))):
        t = open(f, encoding="utf-8").read()
        first = t.splitlines()[0] if t else ""
        stage = "plan" if head_rev.match(first) else ("diagnose" if head_diag.match(first) else "")
        if not stage or "\n---\n" not in t:
            continue
        hdr, body = t.split("\n---\n", 1)
        cls, why, info = mod.classify_body(stage, body)
        b = os.path.basename(f)
        n_all += 1
        if re.match(r"2026-10-1\d-", b):
            n_recent += 1
        got = info["conclusion"] != "?" if stage == "diagnose" else (info["verdict"] != "?" and "?" not in info["counts"])
        if not got:
            bad.append(f"{b}: {why} {info}")
        if cls != "ok" and not why.startswith("必須の項目が無い: 証拠の欄"):
            other_why.append(f"{b}: {why}")
        h = hdr_v.search(hdr)
        if stage != "diagnose" and h and b not in known_bad_header and (h.group(1), h.group(2)) != (info["verdict"], info["counts"]):
            mismatch.append(f"{b}: 見出し {h.group(1)} {h.group(2)} / 抽出 {info['verdict']} {info['counts']}")
    check(f"(n) 後方互換 (証拠の欄の契約より前の記録、分類は incomplete でよい): 既存の記録 {n_all} 件 "
          f"(うち 2026-10-1x {n_recent} 件) から判定・指摘数・結論が取れる", not bad and n_recent > 0, note="; ".join(bad[:5]))
    check("(n) 後方互換: 判定・指摘数が記録の見出し (旧来の抽出) と一致 (旧来の抽出が誤った 1 件を除く)", not mismatch,
          note="; ".join(mismatch[:5]))
    check("(n) 後方互換: 既存の記録が ok でない理由は証拠の欄が無いことだけ (判定の相反・証拠不足は無い)", not other_why,
          note="; ".join(other_why[:5]))

    n_ng = sum(1 for ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_ng}/{len(RESULTS)} ok" + (f", {n_ng} NG" if n_ng else ""))
    if n_ng == 0:
        shutil.rmtree(fx.tmp, ignore_errors=True)
    else:
        print(f"一時ディレクトリを残した: {fx.tmp}")
    return 1 if n_ng else 0


if __name__ == "__main__":
    sys.exit(main())
