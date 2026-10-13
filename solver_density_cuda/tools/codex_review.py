#!/usr/bin/env python3
"""codex に plan の外部レビューを依頼し、結果を notes/reviews/ に残す (AGENTS.md 「codex レビュー」の実体化)。
codex が使えないときは、文脈を分けた Claude の別インスタンスで回し直す (下の「判定役」)。

段階は 3 つ:
  --stage diagnose 診断・設計判断の諮問 (AGENTS.md「モデル分担とエスカレーション」の 7 条件)。旧 `diagnostician`
                  サブエージェント (Fable) の置き換え (2026-09-26 ユーザ決定、usage 削減)。ブリーフを --brief で渡す。
  --stage plan    計画立案時 (設計方針 §4 が書けた時点、実装着手前)。方針・スコープ・検証計画の妥当性を問う。
  --stage result  検証結果が出た時 (status を done にして accepted/ へ移す前)。実装 diff と run の VERDICT を
                  突き合わせ、変更ログの主張が裏付けられているか・accepted にしてよいかを問う。

起動作法 (過去のハング事故 2 件の再発防止):
  * `codex exec --sandbox read-only` を stdin=/dev/null で起動する (stdin を開いたままにすると EOF 待ちで永久ブロック)。
  * stdout はパイプに通さずログファイルへ直書きする (tail 等を挟むと終了まで何も見えない)。
  * 所要 5〜15 分。Claude Code からは run_in_background + timeout 1200 s 以上で呼ぶ。

判定役 (--backend、plan tooling-review-fallback-backend §4.1・§4.2):
  auto   (既定) スイッチ ~/.config/forge/diagnose-backend の中身が `claude` なら最初から claude。そうでなければ codex を試し、
         codex が unavailable (CLI 無し・認証・usage 上限・rate limit・overloaded・timeout・最終メッセージが空) のときだけ、
         理由を記録して claude で回し直す。有効な NO-GO・incomplete・failed では切り替えない (弱い判定役で通し直さない)。
  codex  codex だけ (失敗しても切り替えない。今までの振る舞い)。
  claude 最初から claude。**文脈を共有しない Claude の別インスタンス** (同じ系統のモデルになりうる。独立性は文脈の分離だけ)。
         `claude -p --safe-mode` (CLAUDE.md・auto memory・hook・skill・plugin・MCP を読まない) で、道具は Read・Grep・Glob だけ
         (Bash を渡さない: 許可した git のコマンドにも `git diff --output=<file>` の書き込みの経路がある)。
         git の情報は親 (このスクリプト) が証拠一式 notes/reviews/<記録名>.evidence/ (git 追跡外) にして渡す。

結果の分類 (codex・claude 共通):
  ok           段ごとの必須の項目がある (diagnose: 「結論:」の行、plan・result: 判定 GO / GO-with-changes / NO-GO と「指摘数:」の行)
  unavailable  上の定義
  failed       それ以外の非 0 終了
  incomplete   本文はあるが必須の項目が無い、または「証拠不足:」の行がある
  **ok のときだけ** plan の §6.1 に貼る行を表示して 0 で終わる。それ以外は記録の先頭に「**未成立 (分類)**」と書き、
  非 0 で終わる (failed 1 / incomplete 2 / unavailable 3)。

使い方:
  python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage plan
  python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage result --base main
  python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --focus "§4.2 の作動点定義に集中" --extra case/46.sern_design/README.md
  python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --dry-run   # プロンプトだけ表示
  python3 solver_density_cuda/tools/codex_review.py --stage diagnose --brief notes/reviews/briefs/<日付>-<slug>.md \
      [plans/active/<plan>.md] [--extra case/NN/README.md ...]           # 諮問 (plan は任意)
  python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --backend claude [--claude-model opus] \
      [--requester-model claude-opus-5-5]                                # Claude の別インスタンスで

出力:
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.md          記録 (見出し + 判定役の最終メッセージ)
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.codex.log   codex の試行の生ログ (*.log は git 追跡外)
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.claude.log  claude の試行の生ログ (JSON の出力 + stderr)
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.evidence/   claude に渡した証拠一式 (git 追跡外)
  diagnose は stem にブリーフのファイル名を使う (notes/reviews/YYYY-MM-DD-<brief stem>-diagnose.md)
  試行のログは切り替えたときも全部残す。
最後に (ok のときだけ) plan の「レビュー記録」表へ貼る行を表示する。
"""
import argparse, datetime, json, os, re, shutil, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT_DIR = os.path.join(ROOT, "notes", "reviews")

# 判定役の切り替えスイッチ (AGENTS.md「モデル分担とエスカレーション」)。中身が `claude` なら最初から Claude の別インスタンス。
SWITCH_FILE = "~/.config/forge/diagnose-backend"
# これを超えるプロンプトは引数で渡さずファイルに書いて読ませる (execve の引数長上限、2026-09-20)
PROMPT_ARG_MAX_BYTES = 96 * 1024
# 子の claude に渡す道具 (読むだけ)。Bash は渡さない (plan §4.2、諮問 review-fallback-adoption)
CLAUDE_TOOLS = "Read,Grep,Glob"
# codex の --effort → claude の --effort (claude は low/medium/high/xhigh/max を受ける。同じ名前はそのまま写す)
CLAUDE_EFFORT = {"low": "low", "medium": "medium", "high": "high", "xhigh": "xhigh"}
# 管理ポリシー (managed settings) は --safe-mode でも外せない。あれば記録に残す
MANAGED_SETTINGS_PATHS = ("/etc/claude-code/managed-settings.json", "/etc/claude-code/managed-settings.d",
                          "/etc/claude-code/managed-mcp.json")
# 証拠一式: result 段の git diff がこれを超えたら、stat と上位のファイルの差分に切る
EVIDENCE_DIFF_MAX_BYTES = 512 * 1024
# 証拠一式: プロンプトに直に入れる diff の stat の上限 (超えたら頭だけ入れてファイルを読ませる)
EVIDENCE_STAT_INLINE_MAX_BYTES = 16 * 1024
# 差分を切るときに優先するパス (result 段のプロンプトが見る範囲 + plans)
EVIDENCE_DIFF_PRIORITY = ("solver_density_cuda/", "design/", "methods/", "procedures/", "plans/")
EXIT_CODE = {"ok": 0, "failed": 1, "incomplete": 2, "unavailable": 3}

# unavailable の型 (分類, 正規表現)。codex はログの末尾のエラー行、claude は is_error のときの本文・stderr に当てる。
# **成功した試行 (codex: rc 0 かつ最終メッセージあり、claude: is_error でない) の本文には当てない**
# (本文に「rate limit」と書いてあるだけの正常なレビューを unavailable にしないため)。
UNAVAILABLE_PATTERNS = (
    ("usage_limit", re.compile(r"usage limit|hit your usage|purchase more credits|insufficient_quota|quota exceeded"
                               r"|credit balance is too low|out of credits", re.I)),
    ("rate_limit", re.compile(r"rate[ _-]?limit|\b429\b|too many requests", re.I)),
    ("overloaded", re.compile(r"overloaded|\b529\b|\b503\b|service unavailable|server is busy", re.I)),
    ("auth", re.compile(r"\b401\b|\b403\b|unauthori[sz]ed|not logged in|/login\b|log ?in again|please (?:run )?login"
                        r"|authenticat|invalid api key|token (?:has )?expired|forbidden", re.I)),
)
# codex のログでエラーとして読む行 (codex 自身が出す `ERROR: ...` と stream の失敗)。ログの末尾だけを見る。
CODEX_ERROR_LINE_RE = re.compile(r"^\s*(?:ERROR|Error|error)\b.*$|^.*(?:stream error|unexpected status|stream disconnected).*$",
                                 re.M)
CODEX_ERROR_TAIL_LINES = 200
# codex のログの見出しの `model: <名前>` (実効のモデル)
CODEX_MODEL_RE = re.compile(r"^model:\s*(\S+)\s*$", re.M)
# claude の api_error_status → unavailable の型
CLAUDE_UNAVAILABLE_STATUS = {429: "rate_limit", 529: "overloaded", 503: "overloaded", 401: "auth", 403: "auth"}

# 段ごとの必須の項目。「結論:」は行頭の `**`・`-`・`#`・`>` を許す (`**結論: …**` の記録が既存にある)
CONCLUSION_RE = re.compile(r"^[ \t>*#\-]*\**[ \t]*結論[ \t]*\**[ \t]*[:：][ \t]*(.+)$", re.M)
COUNTS_RE = re.compile(r"指摘数\s*[:：]\s*Critical\s*(\d+)\s*/\s*Major\s*(\d+)\s*/\s*Minor\s*(\d+)")
VERDICT_RE = re.compile(r"\b(GO-with-changes|NO-GO|GO)\b")
# 「証拠不足: <要るもの>」の行 → incomplete。「証拠不足: 無し」は数えない
INSUFFICIENT_RE = re.compile(r"^[ \t>*\-]*\**[ \t]*証拠不足[ \t]*\**[ \t]*[:：][ \t]*(.*)$", re.M)
INSUFFICIENT_NONE_RE = re.compile(r"^\**\s*(?:無し|なし|特になし|特に無し|該当なし|該当無し|none|n/?a|-)\s*\**\s*[。.]?$", re.I)

COMMON_HEADER = """forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。
"""

STAGE_PLAN = """## 依頼: 計画立案時レビュー (stage = plan)

対象の plan は下に全文を貼る (`{plan_rel}`)。これから実装に入る前の段階なので、次を順に評価せよ。

1. **目的とスコープ** (§1, §2): 解こうとしている課題は正しく同定されているか。既に解決済み/別 plan と重複していないか
   (`plans/README.md` と `plans/accepted/` を確認)。
2. **設計方針** (§4): 数学的・数値的に健全か。forge の既存構造 (`methods/architecture/overview.md`、該当 `methods/<area>/`)
   と整合するか。node/cell 両離散化、float32、陰解法 (block-DPLUR)、周期・軸対称などの既知の落とし穴に抵触しないか。
   代替案と比べて費用対効果は妥当か。
3. **実装ステップと残作業表** (§5, §5.1): 順序・粒度は妥当か。抜けている前提 (メッシュ品質、IC、段階起動) はないか。
4. **検証計画** (§6): 判定基準は定量的か。検証ケースの選択は `procedures/verification/README.md` と整合するか。
   「収束」「一致」を何で判定するかが書かれているか。
5. **見落としているリスク**: 我々が気づいていない構造的問題があれば挙げよ。

最後に「この計画で実装に進んでよいか」を **GO / GO-with-changes / NO-GO** の 1 語で判定し、
GO-with-changes なら実装前に直すべき点を優先順で列挙すること。
"""

STAGE_RESULT = """## 依頼: 検証結果レビュー (stage = result)

対象の plan は下に全文を貼る (`{plan_rel}`)。実装と検証が終わり、`status: done` にして `plans/accepted/` へ移す直前の段階である。
次を順に評価せよ。

1. **実装 diff の検証**: `git diff {base}...HEAD -- solver_density_cuda design methods procedures` (必要なら `git log {base}..HEAD --oneline`)
   を自分で取り、plan §4 の設計方針どおりに実装されているか、符号・単位・境界 (node の境界半割面、周期 seam、軸)・
   float32 桁落ち・ゼロ割ガードの絶対閾値などの誤りがないかを見る。
2. **検証結果の裏付け**: plan の §6 / 変更ログに書かれた数値・主張 (収束、一致、改善率) を、挙げられている `run_*`
   ディレクトリの `CONVERGENCE_VERDICT.txt` / `residual_history.csv` / `README.md` の run 一覧で確認する。
   主張と実測が食い違う箇所、VERDICT が NOT CONVERGED / DRIFTING のまま「一致」と書いている箇所を挙げよ。
3. **回帰**: 既存機能 (既定値、他ケース) を壊していないか。既定挙動が変わった場合にそれが plan に明記されているか。
4. **文書整合**: `methods/<area>/` の現在仕様、`procedures/` の手順、`methods/index.md`、`plans/README.md` が
   実装と一致しているか。
5. **残作業表**: 未解決事項が §5.1 の残作業表に残っているか (対話で決めて書いていない、が無いか)。

最後に「accepted に移してよいか」を **GO / GO-with-changes / NO-GO** の 1 語で判定し、
GO-with-changes なら移す前に直すべき点を優先順で列挙すること。
"""

STAGE_DIAGNOSE = """## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。
"""

DIAG_OUTPUT = """## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
"""


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def git(*args):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


def build_diagnose_prompt(brief_path, plan_path, extras, focus):
    parts = [COMMON_HEADER, STAGE_DIAGNOSE]
    if focus:
        parts.append("## 重点\n\n" + focus.strip() + "\n")
    brief_rel = os.path.relpath(brief_path, ROOT)
    parts.append(f"## ブリーフ (`{brief_rel}`)\n\n{read(brief_path).strip()}\n")
    if plan_path:
        plan_rel = os.path.relpath(plan_path, ROOT)
        parts.append(f"## 関連 plan 全文 (`{plan_rel}`)\n\n```markdown\n{read(plan_path).strip()}\n```\n")
    for ex in extras:
        rel = os.path.relpath(os.path.abspath(ex), ROOT)
        if os.path.isfile(ex):
            parts.append(f"## 参考: `{rel}`\n\n```\n{read(ex).strip()}\n```\n")
        else:
            parts.append(f"## 参考: `{rel}` (ディレクトリ。中を自分で読むこと)\n")
    parts.append(DIAG_OUTPUT)
    return "\n".join(parts)


def build_prompt(plan_path, stage, base, extras, focus):
    plan_rel = os.path.relpath(plan_path, ROOT)
    body = STAGE_PLAN if stage == "plan" else STAGE_RESULT
    parts = [COMMON_HEADER, body.format(plan_rel=plan_rel, base=base)]
    if focus:
        parts.append("## 重点\n\n" + focus.strip() + "\n")
    parts.append(f"## plan 全文 (`{plan_rel}`)\n\n```markdown\n{read(plan_path).strip()}\n```\n")
    for ex in extras:
        rel = os.path.relpath(os.path.abspath(ex), ROOT)
        if os.path.isfile(ex):
            parts.append(f"## 参考: `{rel}`\n\n```\n{read(ex).strip()}\n```\n")
        else:
            parts.append(f"## 参考: `{rel}` (ディレクトリ。中を自分で読むこと)\n")
    parts.append(
        "## 出力形式\n\n"
        "1. 冒頭に **判定 (GO / GO-with-changes / NO-GO)** と 3 行以内の要約。\n"
        "2. 指摘一覧 (Critical → Major → Minor の順、番号付き。各項目に根拠と対案)。\n"
        "3. 推奨 (1 つに絞る)。\n"
        "4. 末尾に `指摘数: Critical N / Major N / Minor N` の 1 行。\n")
    return "\n".join(parts)


def git_raw(*args, timeout=120):
    """git を ROOT で走らせ (rc, stdout, stderr) を返す (strip しない。証拠一式用)。"""
    try:
        p = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=timeout,
                           errors="replace")
        return p.returncode, p.stdout, p.stderr
    except Exception as e:
        return 1, "", str(e)


def rel_or_abs(path):
    """ROOT の中なら ROOT からの相対パス、外なら絶対パス。"""
    path = os.path.abspath(path)
    rel = os.path.relpath(path, ROOT)
    return path if rel.startswith("..") else rel


def cli_version(cmd):
    try:
        p = subprocess.run([cmd, "--version"], stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=30)
        out = (p.stdout or p.stderr).strip().splitlines()
        return out[0].strip() if out else "不明"
    except Exception:
        return "不明"


def read_switch():
    """スイッチ ~/.config/forge/diagnose-backend の中身 (無ければ "")。"""
    path = os.path.expanduser(SWITCH_FILE)
    try:
        return read(path).strip().lower()
    except OSError:
        return ""


def match_unavailable(text):
    """エラーの文 → unavailable の型 (当たらなければ "")。"""
    for kind, rx in UNAVAILABLE_PATTERNS:
        if rx.search(text or ""):
            return kind
    return ""


def excerpt(text, n=300):
    text = " ".join((text or "").split())
    return text if len(text) <= n else text[:n] + " …"


def classify_body(stage, text):
    """判定役の最終メッセージ → (分類 ok/incomplete, 理由, 抽出した項目)。codex・claude 共通。"""
    info = {}
    missing = []
    if stage == "diagnose":
        m = CONCLUSION_RE.search(text)
        info["conclusion"] = m.group(1).strip().strip("*").strip() if m else "?"
        if not m:
            missing.append("「結論:」の行")
    else:
        m = COUNTS_RE.search(text)
        v = VERDICT_RE.search(text)
        info["counts"] = f"C{m.group(1)}/M{m.group(2)}/m{m.group(3)}" if m else "C?/M?/m?"
        info["verdict"] = v.group(1) if v else "?"
        if not v:
            missing.append("判定 (GO / GO-with-changes / NO-GO)")
        if not m:
            missing.append("「指摘数: Critical n / Major n / Minor n」の行")
    lacks = [x.strip() for x in INSUFFICIENT_RE.findall(text) if not INSUFFICIENT_NONE_RE.match(x.strip())]
    if lacks:
        return "incomplete", "証拠不足: " + " / ".join(excerpt(x, 120) or "(内容なし)" for x in lacks), info
    if missing:
        return "incomplete", "必須の項目が無い: " + "・".join(missing), info
    return "ok", "", info


class Attempt:
    """判定役 1 回の試行の結果。"""

    def __init__(self, backend, log):
        self.backend = backend
        self.log = log
        self.cls = ""            # ok / unavailable / failed / incomplete
        self.kind = ""           # unavailable の型 (cli_missing / timeout / empty / usage_limit / rate_limit / overloaded / auth)
        self.detail = ""         # 分類の理由 (1 行)
        self.result = ""         # 最終メッセージ
        self.info = {}           # 結論 / 判定・指摘数
        self.rc = None
        self.dt = 0.0
        self.requested_model = ""
        self.effective_model = "不明"
        self.model_usage = {}
        self.cli_version = "不明"
        self.permission_denials = None
        self.cost_usd = None
        self.evidence_dir = ""
        self.ran = False         # CLI を起動したか

    def label(self):
        return self.cls + (f"/{self.kind}" if self.kind else "")


def run_codex(a, prompt, name, out_dir):
    """codex exec を 1 回 (今までの起動作法のまま)。"""
    at = Attempt("codex", os.path.join(out_dir, name + ".codex.log"))
    at.requested_model = a.model or "既定 (~/.codex/config.toml)"
    if not shutil.which("codex"):
        with open(at.log, "w", encoding="utf-8") as lf:
            lf.write("[codex_review] codex CLI が PATH に無い (~/.local/bin/codex)\n")
        at.cls, at.kind, at.detail = "unavailable", "cli_missing", "codex CLI が PATH に無い"
        return at
    at.cli_version = cli_version("codex")
    last_msg = os.path.join(out_dir, name + ".last.txt")
    cmd = ["codex", "exec", "--sandbox", "read-only", "-C", ROOT, "--color", "never",
           "-c", f'model_reasoning_effort="{a.effort}"', "-o", last_msg]
    if a.model:
        cmd += ["-m", a.model]
    # **長いプロンプトは引数で渡せない** (2026-09-20)。plan が育つと `execve` の
    # 引数長上限に当たり `OSError: [Errno 7] Argument list too long` で落ちる
    # (case/49 の plan で 142 kB)。上限に近ければ**ファイルに落として読ませる**。
    # codex は read-only サンドボックスでもリポジトリ内のファイルを読める。
    if len(prompt.encode()) > PROMPT_ARG_MAX_BYTES:
        prompt_file = os.path.join(out_dir, name + ".prompt.md")
        with open(prompt_file, "w", encoding="utf-8") as pf:
            pf.write(prompt)
        rel = os.path.relpath(prompt_file, ROOT)
        cmd.append("レビュー依頼の全文は `%s` にある。**まずこのファイルを読み**、"
                   "その指示どおりにレビューして結果を出力せよ。"
                   "(プロンプトが %d kB と長く引数で渡せないためファイル渡しにしている)"
                   % (rel, len(prompt.encode()) // 1024))
        print("  prompt: %s (%d kB, 引数長上限のためファイル渡し)"
              % (rel, len(prompt.encode()) // 1024))
    else:
        cmd.append(prompt)
    print(f"  codex : log {rel_or_abs(at.log)}  (進捗はこのファイルのサイズ、または ~/.codex/sessions/ の rollout jsonl で確認)")
    sys.stdout.flush()

    t0 = datetime.datetime.now()
    timed_out = False
    with open(at.log, "w", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd[:-1]) + " '<prompt>'\n\n")
        lf.flush()
        try:
            p = subprocess.run(cmd, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=lf, stderr=subprocess.STDOUT,
                               timeout=a.timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            lf.write(f"\n[codex_review] timeout {a.timeout}s\n")
            rc = 124
            timed_out = True
    at.ran = True
    at.rc = rc
    at.dt = (datetime.datetime.now() - t0).total_seconds()
    at.result = read(last_msg).strip() if os.path.exists(last_msg) else ""
    if os.path.exists(last_msg):
        os.remove(last_msg)
    log_text = read(at.log) if os.path.exists(at.log) else ""
    m = CODEX_MODEL_RE.search("\n".join(log_text.splitlines()[:40]))
    if m:
        at.effective_model = m.group(1)

    if timed_out:
        at.cls, at.kind, at.detail = "unavailable", "timeout", f"timeout {a.timeout}s"
    elif rc != 0 or not at.result:
        # 失敗した試行だけ、ログの末尾のエラー行で型を見分ける
        tail = "\n".join(log_text.splitlines()[-CODEX_ERROR_TAIL_LINES:])
        err = "\n".join(x.group(0).strip() for x in CODEX_ERROR_LINE_RE.finditer(tail))
        kind = match_unavailable(err)
        if kind:
            at.cls, at.kind, at.detail = "unavailable", kind, excerpt(err.splitlines()[-1] if err else "")
        elif not at.result:
            at.cls, at.kind = "unavailable", "empty"
            at.detail = f"最終メッセージが空 (rc={rc})" + (f": {excerpt(err.splitlines()[-1])}" if err else "")
        else:
            at.cls, at.detail = "failed", f"codex が rc={rc} で終わった" + (f": {excerpt(err.splitlines()[-1])}" if err else "")
    else:
        at.cls, at.detail, at.info = classify_body(a.stage, at.result)
    return at


def build_evidence(ev_dir, stage, base, diff_max_bytes=EVIDENCE_DIFF_MAX_BYTES):
    """claude に渡す証拠一式を ev_dir に書く (plan §4.2 入力の契約)。

    全段: head.txt (HEAD・ブランチ)。result 段: diff_stat.txt・diff.patch・log.txt。
    diff が diff_max_bytes を超えたら、優先パス → 変更行数の多い順にファイル単位で入れ、入らなかったファイルは README.md に列挙する。
    返り値: {"files": [(名前, 説明)], "stat_inline": プロンプトに直に入れる stat, "notes": [注記]}
    """
    os.makedirs(ev_dir, exist_ok=True)
    files, notes = [], []
    head = git("rev-parse", "HEAD")
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    subject = git("log", "-1", "--format=%s")
    with open(os.path.join(ev_dir, "head.txt"), "w", encoding="utf-8") as f:
        f.write(f"HEAD: {head or '(取れない)'}\nbranch: {branch or '(取れない)'}\nsubject: {subject}\n")
    files.append(("head.txt", "HEAD のコミット・ブランチ (`git rev-parse HEAD` / `--abbrev-ref HEAD`)"))
    stat_inline = ""
    if stage == "result":
        rc, stat, err = git_raw("diff", "--stat=200", f"{base}...HEAD")
        if rc != 0:
            notes.append(f"`git diff --stat {base}...HEAD` が失敗 (rc={rc}): {excerpt(err, 200)}")
        with open(os.path.join(ev_dir, "diff_stat.txt"), "w", encoding="utf-8") as f:
            f.write(stat)
        files.append(("diff_stat.txt", f"`git diff --stat {base}...HEAD`"))
        rc, log, err = git_raw("log", "--oneline", f"{base}..HEAD")
        if rc != 0:
            notes.append(f"`git log --oneline {base}..HEAD` が失敗 (rc={rc}): {excerpt(err, 200)}")
        with open(os.path.join(ev_dir, "log.txt"), "w", encoding="utf-8") as f:
            f.write(log)
        files.append(("log.txt", f"`git log --oneline {base}..HEAD`"))
        rc, diff, err = git_raw("diff", f"{base}...HEAD")
        if rc != 0:
            notes.append(f"`git diff {base}...HEAD` が失敗 (rc={rc}): {excerpt(err, 200)}")
        size = len(diff.encode())
        if size <= diff_max_bytes:
            with open(os.path.join(ev_dir, "diff.patch"), "w", encoding="utf-8") as f:
                f.write(diff)
            files.append(("diff.patch", f"`git diff {base}...HEAD` (全文, {size // 1024} KiB)"))
        else:
            kept, dropped = truncate_diff(diff, diff_max_bytes)
            with open(os.path.join(ev_dir, "diff.patch"), "w", encoding="utf-8") as f:
                f.write("".join(c for _, c, _ in kept))
            files.append(("diff.patch", f"`git diff {base}...HEAD` を切ったもの ({len(kept)} ファイル。全文 {size // 1024} KiB は"
                                        f"上限 {diff_max_bytes // 1024} KiB を超える。入らなかったファイルは README.md)"))
            notes.append(f"diff は上限 {diff_max_bytes // 1024} KiB を超えた ({size // 1024} KiB) ので、"
                         f"{', '.join(EVIDENCE_DIFF_PRIORITY)} → その他 の順・変更行数の多い順にファイル単位で入れた。"
                         f"入らなかった {len(dropped)} ファイルは下に列挙 (今の作業ツリーの中身は Read で読める)")
            notes += [f"  - 省略: `{p}` ({len(c.encode()) // 1024} KiB, 変更 {n} 行)" for p, c, n in dropped]
        b = stat.encode()
        if len(b) <= EVIDENCE_STAT_INLINE_MAX_BYTES:
            stat_inline = stat
        else:
            lines = stat.splitlines()
            keep, total = [], 0
            for ln in lines[:-1]:
                total += len(ln.encode()) + 1
                if total > EVIDENCE_STAT_INLINE_MAX_BYTES:
                    break
                keep.append(ln)
            stat_inline = "\n".join(keep + [f" ... (残り {len(lines) - 1 - len(keep)} ファイルは diff_stat.txt)",
                                             lines[-1] if lines else ""]) + "\n"
    with open(os.path.join(ev_dir, "README.md"), "w", encoding="utf-8") as f:
        f.write("# 証拠一式 (codex_review.py が作成。子の claude は Read・Grep・Glob しか持たないので、git の情報はここで渡す)\n\n")
        f.write(f"- stage: `{stage}`" + (f" (diff base `{base}`)" if stage == "result" else "") + "\n")
        for fn, desc in files:
            f.write(f"- `{fn}`: {desc}\n")
        if notes:
            f.write("\n## 注記\n\n" + "\n".join(n if n.startswith("  ") else f"- {n}" for n in notes) + "\n")
    files.insert(0, ("README.md", "この証拠一式の目次と注記 (切った差分の一覧)"))
    return {"files": files, "stat_inline": stat_inline, "notes": notes}


def truncate_diff(diff, max_bytes):
    """git diff の全文をファイル単位に分け、優先パス → 変更行数の多い順に max_bytes まで入れる。
    返り値: (入れた [(パス, 差分, 変更行数)] (元の順), 入らなかった [(パス, 差分, 変更行数)])"""
    chunks = re.split(r"(?m)^(?=diff --git )", diff)
    items = []
    for i, c in enumerate(chunks):
        if not c.startswith("diff --git "):
            continue
        first = c.split("\n", 1)[0]
        path = first.rsplit(" b/", 1)[-1] if " b/" in first else first[len("diff --git "):]
        n = sum(1 for ln in c.splitlines() if (ln.startswith("+") and not ln.startswith("+++"))
                or (ln.startswith("-") and not ln.startswith("---")))
        grp = next((k for k, pre in enumerate(EVIDENCE_DIFF_PRIORITY) if path.startswith(pre)), len(EVIDENCE_DIFF_PRIORITY))
        items.append((grp, -n, i, path, c, n))
    kept, dropped, total = [], [], 0
    for grp, _, i, path, c, n in sorted(items):
        sz = len(c.encode())
        if total + sz <= max_bytes:
            kept.append((i, path, c, n))
            total += sz
        else:
            dropped.append((path, c, n))
    kept.sort()
    return [(p, c, n) for _, p, c, n in kept], dropped


def claude_preface(ev_rel, ev, stage, base):
    """claude に渡すプロンプトの前書き (codex と同じプロンプトの前に置く。plan §4.2)。"""
    lines = [
        "## この依頼の前書き (Claude の別インスタンス向け。下の本文の指示より優先する)",
        "",
        "- あなたは codex の代わりにこの依頼を受けている、**依頼側と文脈を共有しない Claude の別インスタンス**である。"
        "cwd はリポジトリのルート。",
        "- **道具は Read・Grep・Glob だけ**。シェル・git・python は使えず、ファイルを変更することもできない "
        "(本文の「読む・実行して確認するのは可」「`git diff` を自分で取り」「`python3` で統計量を出す」は、次の読み替えに従う)。",
        f"- **git の情報 (HEAD・ブランチ・差分・履歴) は証拠一式 `{ev_rel}/` のファイルにある**:",
    ]
    lines += [f"  - `{ev_rel}/{fn}`: {desc}" for fn, desc in ev["files"]]
    lines += [
        "- `sed -n 'A,Bp' f` のような行範囲の指示は、Read の offset/limit (offset=A, limit=B−A+1) で読む。長いファイルは分けて読む。",
        "- 実行しないと得られない数値 (HDF5 の統計・判定ツールの出力など) は、ブリーフ・証拠一式・リポジトリの既存のファイル "
        "(`CONVERGENCE_VERDICT.txt` など) にあるものだけを根拠にする。",
        "- 依頼側 (この依頼を書いたモデル) は**あなたと同じ系統のモデルかもしれない**。依頼側の結論・仮説・「除外済み」に同調せず、"
        "根拠の行 (`ファイル:行`) を自分で読んで確かめてから採る。確かめられないものは採らない。",
        "- 判定に要る証拠が足りなければ、推測で埋めずに `証拠不足: <要るもの>` の行を (足りないものごとに) 出す。"
        "この行があると、この依頼は不成立 (incomplete) として扱われる。足りているときはこの行を書かない (「証拠不足: 無し」とも書かない)。",
        "- 出力の形式は本文の「出力形式」に従う。",
    ]
    if stage == "result" and ev.get("stat_inline"):
        lines += ["", f"### `git diff --stat {base}...HEAD` (証拠一式 `diff_stat.txt` と同じ)", "", "```",
                  ev["stat_inline"].rstrip("\n"), "```"]
    return "\n".join(lines) + "\n\n---\n\n"


def managed_policy():
    found = [p for p in MANAGED_SETTINGS_PATHS if os.path.exists(p)]
    return ("有り: " + ", ".join(f"`{p}`" for p in found) + " (--safe-mode でも効く)") if found else \
        "無し (" + ", ".join(f"`{p}`" for p in MANAGED_SETTINGS_PATHS) + " が無い)"


def run_claude(a, prompt, name, out_dir):
    """claude -p を 1 回 (文脈を分けた別インスタンス、道具は Read・Grep・Glob だけ)。"""
    at = Attempt("claude", os.path.join(out_dir, name + ".claude.log"))
    at.requested_model = a.claude_model
    if not shutil.which("claude"):
        with open(at.log, "w", encoding="utf-8") as lf:
            lf.write("[codex_review] claude CLI が PATH に無い (~/.local/bin/claude)\n")
        at.cls, at.kind, at.detail = "unavailable", "cli_missing", "claude CLI が PATH に無い"
        return at
    at.cli_version = cli_version("claude")
    at.evidence_dir = os.path.join(out_dir, name + ".evidence")
    ev = build_evidence(at.evidence_dir, a.stage, a.base)
    ev_rel = rel_or_abs(at.evidence_dir)
    full = claude_preface(ev_rel, ev, a.stage, a.base) + prompt

    cmd = ["claude", "-p", "--safe-mode", "--model", a.claude_model, "--effort", CLAUDE_EFFORT[a.effort],
           "--tools", CLAUDE_TOOLS, "--permission-mode", "default", "--no-session-persistence", "--strict-mcp-config"]
    if rel_or_abs(out_dir) == os.path.abspath(out_dir):
        # 記録の置き場がリポジトリの外 (試験など) なら、証拠一式・プロンプトを読めるようにする
        cmd += ["--add-dir", os.path.abspath(out_dir)]
    cmd += ["--output-format", "json"]
    if len(full.encode()) > PROMPT_ARG_MAX_BYTES:
        prompt_file = os.path.join(out_dir, name + ".claude.prompt.md")
        with open(prompt_file, "w", encoding="utf-8") as pf:
            pf.write(full)
        rel = rel_or_abs(prompt_file)
        cmd.append("依頼の全文は `%s` にある。**まずこのファイルを Read で最後まで読み** (長いので offset/limit で分けて)、"
                   "その指示どおりに結果を出力せよ。道具は Read・Grep・Glob だけ。"
                   "(プロンプトが %d kB と長く引数で渡せないためファイル渡しにしている)"
                   % (rel, len(full.encode()) // 1024))
        print("  prompt: %s (%d kB, 引数長上限のためファイル渡し)" % (rel, len(full.encode()) // 1024))
    else:
        cmd.append(full)
    print(f"  claude: log {rel_or_abs(at.log)}  evidence {ev_rel}/  model {a.claude_model}")
    sys.stdout.flush()

    t0 = datetime.datetime.now()
    timed_out = False
    stderr = ""
    with open(at.log, "w", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd[:-1]) + " '<prompt>'\n\n")
        lf.flush()
        out_start = os.lseek(lf.fileno(), 0, os.SEEK_CUR)   # ここから子の stdout (JSON)
        try:
            p = subprocess.run(cmd, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=lf, stderr=subprocess.PIPE,
                               text=True, errors="replace", timeout=a.timeout)
            rc, stderr = p.returncode, p.stderr or ""
        except subprocess.TimeoutExpired as e:
            rc, timed_out = 124, True
            stderr = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else (e.stderr or "")
        out_end = os.lseek(lf.fileno(), 0, os.SEEK_CUR)
        lf.write("\n[stderr]\n" + stderr + ("\n" if stderr and not stderr.endswith("\n") else ""))
        if timed_out:
            lf.write(f"[codex_review] timeout {a.timeout}s\n")
    at.ran = True
    at.rc = rc
    at.dt = (datetime.datetime.now() - t0).total_seconds()

    with open(at.log, "rb") as f:
        f.seek(out_start)
        out = f.read(max(out_end - out_start, 0)).decode("utf-8", errors="replace")
    data = None
    try:   # stdout 全体が 1 つの JSON (整形されていてもよい)
        data = json.loads(out)
    except ValueError:   # 前後に余計な行があれば、`{` で始まる行を後ろから探す
        for ln in reversed(out.splitlines()):
            if ln.startswith("{"):
                try:
                    data = json.loads(ln)
                    break
                except ValueError:
                    continue
    if isinstance(data, dict):
        at.result = (data.get("result") or "").strip() if isinstance(data.get("result"), str) else ""
        at.model_usage = data.get("modelUsage") or {}
        if at.model_usage:
            at.effective_model = max(at.model_usage, key=lambda k: (at.model_usage[k] or {}).get("outputTokens", 0))
        at.permission_denials = data.get("permission_denials")
        at.cost_usd = data.get("total_cost_usd")

    if timed_out:
        at.cls, at.kind, at.detail = "unavailable", "timeout", f"timeout {a.timeout}s"
    elif not isinstance(data, dict):
        kind = match_unavailable(stderr + "\n" + out)
        if kind:
            at.cls, at.kind, at.detail = "unavailable", kind, excerpt((stderr or out).strip().splitlines()[-1])
        else:
            at.cls, at.detail = "failed", f"claude の JSON が読めない (rc={rc})" + (f": {excerpt(stderr)}" if stderr.strip() else "")
    elif data.get("is_error"):
        status = data.get("api_error_status")
        kind = CLAUDE_UNAVAILABLE_STATUS.get(status if isinstance(status, int) else -1, "") \
            or match_unavailable(f"{data.get('result') or ''}\n{stderr}")
        msg = excerpt(f"{data.get('subtype', '')}: {data.get('result') or ''}" + (f" (api_error_status {status})" if status else ""))
        if kind:
            at.cls, at.kind, at.detail = "unavailable", kind, msg
        else:
            at.cls, at.detail = "failed", f"claude が is_error で終わった (rc={rc}): {msg}"
    elif rc != 0:
        at.cls, at.detail = "failed", f"claude が rc={rc} で終わった" + (f": {excerpt(stderr)}" if stderr.strip() else "")
    elif not at.result:
        at.cls, at.kind, at.detail = "unavailable", "empty", "最終メッセージ (result) が空"
    else:
        at.cls, at.detail, at.info = classify_body(a.stage, at.result)
    return at


def requester_relation(requester, at):
    """依頼側との関係 (同一 / 別 / 不明)。--requester-model を渡されたときだけ判定し、推測で埋めない。"""
    if not requester:
        return "不明 (`--requester-model` 未指定)"
    eff = at.effective_model
    if not eff or eff == "不明":
        return f"不明 (依頼側 `{requester}`、実効のモデルが取れない)"
    r, e = requester.strip().lower(), eff.lower()
    if r == e:
        return f"同一 (依頼側 `{requester}` = 実効 `{eff}`)"
    if r in ("opus", "sonnet", "haiku", "fable") and r in re.split(r"[-_.]", e):
        return f"同一 (依頼側の別名 `{requester}` が実効 `{eff}` に一致)"
    claude_like = lambda s: "claude" in s or any(t in re.split(r"[-_.]", s) for t in ("opus", "sonnet", "haiku", "fable"))
    same_family = " (同じ Claude 系統)" if claude_like(r) and claude_like(e) else ""
    return f"別 (依頼側 `{requester}`、実効 `{eff}`){same_family}"


def write_record(a, ctx, attempts, reason):
    """記録の md を書く。返り値: 最終の試行。"""
    final = attempts[-1]
    claude = final.backend == "claude"
    stage, stem = a.stage, ctx["stem"]
    if stage == "diagnose":
        title = f"# {'Claude 代替諮問' if claude else 'codex 諮問'} (diagnose): {stem}"
    else:
        title = f"# {'Claude 代替レビュー' if claude else 'codex レビュー'}: {stem} ({stage})"
    L = []
    if final.cls != "ok":
        L += [f"**未成立 ({final.cls})** — {final.label()}: {final.detail}", ""]
    L += [title, ""]
    if claude:
        L += ["> **Claude の別インスタンスによる代替** (独立性は文脈の分離だけ。依頼側と同じ系統のモデルでありうる; "
              "codex・Fable と別系統の判定ではない)", ""]
    plan_rel, brief_rel = ctx["plan_rel"], ctx["brief_rel"]
    if stage == "diagnose":
        L.append(f"- **brief**: [`{brief_rel}`](../../{brief_rel})")
        if plan_rel:
            L.append(f"- **plan**: [`{plan_rel}`](../../{plan_rel})")
    else:
        L.append(f"- **plan**: [`{plan_rel}`](../../{plan_rel})")
        L.append(f"- **stage**: `{stage}`" + (f" (diff base `{a.base}`)" if stage == "result" else ""))
    L.append(f"- **date**: {ctx['today']}")
    L.append(f"- **commit**: `{ctx['head']}` ({ctx['branch']})")
    L.append(f"- **分類**: `{final.cls}`" + (f" ({final.kind})" if final.kind else ""))
    L.append(f"- **backend**: `{final.backend}` — 選んだ理由: {reason}")
    L.append(f"- **要求したモデル**: `{final.requested_model}`")
    eff = final.effective_model
    if claude and final.model_usage:
        eff_s = f"`{eff}` (modelUsage: " + ", ".join(
            f"`{k}` out {(v or {}).get('outputTokens', '?')}" for k, v in final.model_usage.items()) + ")"
    else:
        eff_s = f"`{eff}`" if eff != "不明" else "不明"
        if not claude and eff != "不明":
            eff_s += " (codex のログの見出し)"
    L.append(f"- **実効のモデル**: {eff_s}")
    L.append(f"- **CLI の版**: `{final.cli_version}`")
    L.append(f"- **依頼側との関係**: {requester_relation(a.requester_model, final)}")
    if claude:
        pd = final.permission_denials
        if pd is None:
            L.append("- **permission_denials**: 不明 (JSON が取れない)")
        elif not pd:
            L.append("- **permission_denials**: 0 件")
        else:
            L.append(f"- **permission_denials**: {len(pd)} 件 — " + "; ".join(
                excerpt(f"{d.get('tool_name', '?')} {json.dumps(d.get('tool_input', {}), ensure_ascii=False)}", 160)
                if isinstance(d, dict) else excerpt(str(d), 160) for d in pd))
        L.append(f"- **道具**: `{CLAUDE_TOOLS}` (`--safe-mode`、Bash なし) / 管理ポリシー: {managed_policy()}")
    else:
        L.append("- **permission_denials**: 対象外 (codex は read-only サンドボックス)")
    run = f"- **実行**: effort `{a.effort}`" + (f", model `{a.model}`" if (a.model and not claude) else "")
    run += f", {final.dt / 60:.1f} min, rc={final.rc if final.rc is not None else '-'}"
    if claude and final.cost_usd is not None:
        run += f", total_cost_usd {final.cost_usd:.4f}"
    L.append(run)
    L.append("- **試行のログ**: " + ", ".join(f"`{rel_or_abs(t.log)}` ({t.backend}: {t.label()})" for t in attempts))
    if claude and final.evidence_dir:
        L.append(f"- **証拠一式**: `{rel_or_abs(final.evidence_dir)}/` (git 追跡外)")
    if stage == "diagnose":
        L.append(f"- **結論**: {final.info.get('conclusion', '?')}")
    else:
        L.append(f"- **判定**: **{final.info.get('verdict', '?')}**, 指摘 {final.info.get('counts', 'C?/M?/m?')}")
    if a.focus:
        L.append(f"- **focus**: {a.focus}")
    if a.extra:
        L.append(f"- **extra**: {', '.join('`' + os.path.relpath(os.path.abspath(e), ROOT) + '`' for e in a.extra)}")
    who = "Claude の別インスタンスの最終メッセージ (JSON の `result`) " if claude else "codex の最終メッセージ"
    L.append("")
    if stage == "diagnose":
        L.append(f"本文は {who}をそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。"
                 "採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。")
    else:
        L.append(f"本文は {who}をそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。")
    L += ["", "---", ""]
    body = final.result
    if not body:
        tails = []
        for t in attempts:
            txt = read(t.log) if os.path.exists(t.log) else ""
            tails.append(f"### `{rel_or_abs(t.log)}` ({t.backend}: {t.label()}) の末尾\n\n```\n"
                         + "\n".join(txt.splitlines()[-20:]) + "\n```")
        body = "(本文なし)\n\n" + "\n\n".join(tails)
    with open(ctx["out_md"], "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n" + body + "\n")
    return final


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plan", nargs="?", default="", help="plans/active/<plan>.md (diagnose では任意)")
    ap.add_argument("--stage", choices=["plan", "result", "diagnose"], required=True)
    ap.add_argument("--brief", default="", help="diagnose: ブリーフ (観測事実/期待値と出典/再現条件/実施済み/仮説) の md")
    ap.add_argument("--base", default="main", help="result 段階の diff 基準 ref (既定 main)")
    ap.add_argument("--extra", nargs="*", default=[], help="プロンプトに全文を貼る補助ファイル (case README 等)")
    ap.add_argument("--focus", default="", help="重点的に見てほしい点 (自由文)")
    ap.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh"],
                    help="model_reasoning_effort の上書き (~/.codex/config.toml は low なので既定で high に上げる。"
                         "難所の diagnose は xhigh も可)。claude にも同じ値を --effort で渡す")
    ap.add_argument("--model", default="", help="codex の model 上書き (既定は config.toml のもの)")
    ap.add_argument("--backend", default="auto", choices=["auto", "codex", "claude"],
                    help="判定役 (既定 auto: スイッチが claude なら claude、そうでなければ codex を試し unavailable なら claude)")
    ap.add_argument("--claude-model", default="opus", help="claude の --model (既定 opus)")
    ap.add_argument("--requester-model", default="",
                    help="依頼側 (このスクリプトを呼んだモデル) の名前。渡したときだけ記録に 同一 / 別 を書く (無ければ「不明」)")
    ap.add_argument("--name-suffix", default="",
                    help="記録名の末尾に足す識別子 (例: claude → YYYY-MM-DD-<stem>-<stage>-claude.md。試験・同日の区別用)")
    ap.add_argument("--timeout", type=int, default=1800, help="秒 (既定 1800、試行ごと)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--dry-run", action="store_true", help="プロンプトを表示して終了")
    a = ap.parse_args()

    plan_path = os.path.abspath(a.plan) if a.plan else ""
    if plan_path and not os.path.isfile(plan_path):
        sys.exit(f"plan が無い: {a.plan}")
    brief_path = ""
    if a.stage == "diagnose":
        if not a.brief or not os.path.isfile(a.brief):
            sys.exit("diagnose には --brief <md> が要る (観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説)")
        brief_path = os.path.abspath(a.brief)
        prompt = build_diagnose_prompt(brief_path, plan_path, a.extra, a.focus)
        stem = os.path.splitext(os.path.basename(brief_path))[0]
        stem = re.sub(r"^\d{4}-\d{2}-\d{2}-", "", stem)
    else:
        if not plan_path:
            sys.exit("plan / result 段には plan のパスが要る")
        prompt = build_prompt(plan_path, a.stage, a.base, a.extra, a.focus)
        stem = os.path.splitext(os.path.basename(plan_path))[0]

    # 判定役の選択 (plan §4.1)
    switch = read_switch()
    if a.backend == "auto" and switch == "claude":
        first, reason = "claude", f"スイッチ (`{SWITCH_FILE}` = `claude`)"
    elif a.backend == "auto":
        first, reason = "codex", "既定 (auto: codex を試し、unavailable のときだけ claude)"
    else:
        first, reason = a.backend, f"明示 (`--backend {a.backend}`)"

    if a.dry_run:
        if first == "claude":
            ev = {"files": [("README.md", "(dry-run のため未作成)")], "stat_inline": ""}
            print(claude_preface("notes/reviews/<記録名>.evidence", ev, a.stage, a.base) + prompt)
        else:
            print(prompt)
        return 0

    out_dir = os.path.abspath(a.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    today = datetime.date.today().isoformat()
    base_name = f"{today}-{stem}-{a.stage}" + (f"-{a.name_suffix}" if a.name_suffix else "")
    name = base_name
    k = 1
    while os.path.exists(os.path.join(out_dir, name + ".md")):
        k += 1
        name = f"{base_name}-{k}"
    out_md = os.path.join(out_dir, name + ".md")

    head = git("rev-parse", "--short", "HEAD")
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    target = os.path.relpath(plan_path, ROOT) if plan_path else "-"
    print(f"review: stage={a.stage} plan={target} effort={a.effort} backend={a.backend} ({reason})"
          + (f" brief={os.path.relpath(a.brief, ROOT)}" if a.brief else ""))
    print(f"  out : {rel_or_abs(out_md)}")
    sys.stdout.flush()

    attempts = []
    if first == "codex":
        at = run_codex(a, prompt, name, out_dir)
        attempts.append(at)
        if at.cls == "unavailable" and a.backend == "auto":
            reason = f"自動の切り替え (codex が unavailable/{at.kind}: 「{excerpt(at.detail, 200)}」)"
            print(f"codex が unavailable ({at.kind}: {excerpt(at.detail, 200)}) → Claude の別インスタンスで回し直す")
            sys.stdout.flush()
            attempts.append(run_claude(a, prompt, name, out_dir))
    else:
        attempts.append(run_claude(a, prompt, name, out_dir))

    ctx = {"stem": stem, "today": today, "head": head, "branch": branch, "out_md": out_md,
           "plan_rel": os.path.relpath(plan_path, ROOT) if plan_path else "",
           "brief_rel": os.path.relpath(brief_path, ROOT) if brief_path else ""}
    final = write_record(a, ctx, attempts, reason)
    claude = final.backend == "claude"
    rel_md = rel_or_abs(out_md)

    if final.cls != "ok":
        print(f"\n未成立 ({final.label()}): {final.detail}")
        print(f"記録 (先頭に未成立と書いた): {rel_md}")
        print("試行のログ: " + ", ".join(f"{rel_or_abs(t.log)} ({t.backend}: {t.label()})" for t in attempts))
        print("§6.1 の行は出さない (レビュー・諮問は成立していない)。")
        return EXIT_CODE.get(final.cls, 1)

    tag = "[Claude 代替] " if claude else ""
    if a.stage == "diagnose":
        concl = final.info.get("conclusion", "?")
        print(f"\n完了 ({final.dt / 60:.1f} min, rc={final.rc}, backend={final.backend})")
        print(f"諮問記録: {rel_md}")
        print(f"結論: {concl}")
        who = "Claude の別インスタンス (同一系統の代替)" if claude else "codex (diagnose)"
        print(f"応答・plan に書く痕跡: 「{who} に諮った: {rel_md} — 結論 1 行」")
        link = os.path.relpath(out_md, os.path.dirname(plan_path)) if plan_path else rel_md
        print("\nplan の「レビュー記録」表に貼る行 (対応欄は採否を書いてから):")
        print(f"| diagnose | {today} | [{name}.md]({link}) | {tag}{excerpt(concl, 200)} | <採用/却下と理由> |")
        return 0

    verdict, counts = final.info.get("verdict", "?"), final.info.get("counts", "C?/M?/m?")
    print(f"\n完了 ({final.dt / 60:.1f} min, rc={final.rc}, backend={final.backend}) 判定={verdict} 指摘={counts}")
    print(f"レビュー: {rel_md}")
    print("\nplan の「レビュー記録」表に貼る行 (対応欄は採否を書いてから):")
    print(f"| {a.stage} | {today} | [{name}.md]({os.path.relpath(out_md, os.path.dirname(plan_path))}) | {tag}{verdict}, {counts} | <採用/却下と理由。Critical/Major は §5.1 残作業表へ> |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
