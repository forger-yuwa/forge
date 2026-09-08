#!/usr/bin/env python3
"""codex に plan の外部レビューを依頼し、結果を notes/reviews/ に残す (AGENTS.md 「codex レビュー」の実体化)。

段階は 2 つ:
  --stage plan    計画立案時 (設計方針 §4 が書けた時点、実装着手前)。方針・スコープ・検証計画の妥当性を問う。
  --stage result  検証結果が出た時 (status を done にして accepted/ へ移す前)。実装 diff と run の VERDICT を
                  突き合わせ、変更ログの主張が裏付けられているか・accepted にしてよいかを問う。

起動作法 (過去のハング事故 2 件の再発防止):
  * `codex exec --sandbox read-only` を stdin=/dev/null で起動する (stdin を開いたままにすると EOF 待ちで永久ブロック)。
  * stdout はパイプに通さずログファイルへ直書きする (tail 等を挟むと終了まで何も見えない)。
  * 所要 5〜15 分。Claude Code からは run_in_background + timeout 1200 s 以上で呼ぶ。

使い方:
  python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage plan
  python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage result --base main
  python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --focus "§4.2 の作動点定義に集中" --extra case/46.sern_design/README.md
  python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --dry-run   # プロンプトだけ表示

出力:
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.md   レビュー本文 (メタ情報ヘッダ + codex の最終メッセージ)
  notes/reviews/YYYY-MM-DD-<plan stem>-<stage>.log  生ログ (*.log は git 追跡外)
最後に plan の「レビュー記録」表へ貼る行を表示する。
"""
import argparse, datetime, os, re, shutil, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT_DIR = os.path.join(ROOT, "notes", "reviews")

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


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def git(*args):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plan", help="plans/active/<plan>.md")
    ap.add_argument("--stage", choices=["plan", "result"], required=True)
    ap.add_argument("--base", default="main", help="result 段階の diff 基準 ref (既定 main)")
    ap.add_argument("--extra", nargs="*", default=[], help="プロンプトに全文を貼る補助ファイル (case README 等)")
    ap.add_argument("--focus", default="", help="重点的に見てほしい点 (自由文)")
    ap.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh"],
                    help="model_reasoning_effort の上書き (~/.codex/config.toml は low なので既定で high に上げる)")
    ap.add_argument("--model", default="", help="codex の model 上書き (既定は config.toml のもの)")
    ap.add_argument("--timeout", type=int, default=1800, help="秒 (既定 1800)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--dry-run", action="store_true", help="プロンプトを表示して終了")
    a = ap.parse_args()

    plan_path = os.path.abspath(a.plan)
    if not os.path.isfile(plan_path):
        sys.exit(f"plan が無い: {a.plan}")
    prompt = build_prompt(plan_path, a.stage, a.base, a.extra, a.focus)
    if a.dry_run:
        print(prompt)
        return 0

    if not shutil.which("codex"):
        sys.exit("codex CLI が PATH に無い (~/.local/bin/codex)")

    os.makedirs(a.out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(plan_path))[0]
    today = datetime.date.today().isoformat()
    name = f"{today}-{stem}-{a.stage}"
    k = 1
    while os.path.exists(os.path.join(a.out_dir, name + ".md")):
        k += 1
        name = f"{today}-{stem}-{a.stage}-{k}"
    out_md = os.path.join(a.out_dir, name + ".md")
    out_log = os.path.join(a.out_dir, name + ".log")
    last_msg = os.path.join(a.out_dir, name + ".last.txt")

    cmd = ["codex", "exec", "--sandbox", "read-only", "-C", ROOT, "--color", "never",
           "-c", f'model_reasoning_effort="{a.effort}"', "-o", last_msg]
    if a.model:
        cmd += ["-m", a.model]
    cmd.append(prompt)

    head = git("rev-parse", "--short", "HEAD")
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    print(f"codex review: stage={a.stage} plan={os.path.relpath(plan_path, ROOT)} effort={a.effort}")
    print(f"  log : {os.path.relpath(out_log, ROOT)}  (進捗はこのファイルのサイズ、または ~/.codex/sessions/ の rollout jsonl で確認)")
    print(f"  out : {os.path.relpath(out_md, ROOT)}")
    sys.stdout.flush()

    t0 = datetime.datetime.now()
    with open(out_log, "w", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd[:-1]) + " '<prompt>'\n\n")
        lf.flush()
        try:
            p = subprocess.run(cmd, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=lf, stderr=subprocess.STDOUT,
                               timeout=a.timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            lf.write(f"\n[codex_review] timeout {a.timeout}s\n")
            rc = 124
    dt = (datetime.datetime.now() - t0).total_seconds()

    result = read(last_msg).strip() if os.path.exists(last_msg) else ""
    if os.path.exists(last_msg):
        os.remove(last_msg)
    if not result:
        print(f"codex の最終メッセージが取れなかった (rc={rc}, {dt:.0f}s)。ログ: {os.path.relpath(out_log, ROOT)}")
        return 1

    m = re.search(r"指摘数\s*[:：]\s*Critical\s*(\d+)\s*/\s*Major\s*(\d+)\s*/\s*Minor\s*(\d+)", result)
    counts = f"C{m.group(1)}/M{m.group(2)}/m{m.group(3)}" if m else "C?/M?/m?"
    verdict_m = re.search(r"\b(GO-with-changes|NO-GO|GO)\b", result)
    verdict = verdict_m.group(1) if verdict_m else "?"

    header = (
        f"# codex レビュー: {stem} ({a.stage})\n\n"
        f"- **plan**: [`{os.path.relpath(plan_path, ROOT)}`](../../{os.path.relpath(plan_path, ROOT)})\n"
        f"- **stage**: `{a.stage}`" + (f" (diff base `{a.base}`)" if a.stage == "result" else "") + "\n"
        f"- **date**: {today}\n"
        f"- **commit**: `{head}` ({branch})\n"
        f"- **codex**: effort `{a.effort}`" + (f", model `{a.model}`" if a.model else "") + f", {dt/60:.1f} min, rc={rc}\n"
        f"- **判定**: **{verdict}**, 指摘 {counts}\n"
        + (f"- **focus**: {a.focus}\n" if a.focus else "")
        + (f"- **extra**: {', '.join('`'+os.path.relpath(os.path.abspath(e), ROOT)+'`' for e in a.extra)}\n" if a.extra else "")
        + "\n本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。\n\n---\n\n"
    )
    with open(out_md, "w", encoding="utf-8") as f:
        f.write(header + result + "\n")

    rel_md = os.path.relpath(out_md, ROOT)
    print(f"\n完了 ({dt/60:.1f} min, rc={rc}) 判定={verdict} 指摘={counts}")
    print(f"レビュー: {rel_md}")
    print("\nplan の「レビュー記録」表に貼る行 (対応欄は採否を書いてから):")
    print(f"| {a.stage} | {today} | [{name}.md]({os.path.relpath(out_md, os.path.dirname(plan_path))}) | {verdict}, {counts} | <採用/却下と理由。Critical/Major は §5.1 残作業表へ> |")
    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
