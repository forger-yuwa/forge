# codex 諮問 (diagnose): review-fallback-adoption

- **brief**: [`notes/reviews/briefs/2026-10-14-review-fallback-adoption.md`](../../notes/reviews/briefs/2026-10-14-review-fallback-adoption.md)
- **plan**: [`plans/active/tooling-review-fallback-backend.md`](../../plans/active/tooling-review-fallback-backend.md)
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 1.4 min, rc=0
- **結論**: §4.2 の子プロセスから `Bash` を外し、必要な git の差分・履歴は親が取得して証拠一式に渡す設計への変更を推奨する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: §4.2 の子プロセスから `Bash` を外し、必要な git の差分・履歴は親が取得して証拠一式に渡す設計への変更を推奨する。

第 1 仮説: **Major — 許可した git コマンド自体の書き込み経路が残り、現行 V3 が合格しても「読み取り専用」を保証できない。** 確度: 高（検証の穴として。実際の CLI での通過は未確認）
  根拠: [plan:59](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:59) は `Bash(git diff:*)` を許可する。しかし `git diff` には `--output=<file>` があり、シェルのリダイレクトや `Write` を使わずファイルを書ける。[plan:110](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:110) の V3 は `touch`・リダイレクト・`Write`・`Edit` を試すだけで、この経路を必須にしていない。したがって M4 の採用方針は妥当だが、反映完了という評価は**要再検証**。
  反証条件: 宣言した CLI・起動引数で、`git diff --output=<file>` の実際のツール呼び出しが権限機構に拒否され、ファイル集合・ハッシュが不変であること。モデルが実行を控えただけでは反証にならない。

第 2・第 3 仮説: 無し。

判別 A/B: 呼び出し側で、一時コピーに既知の tracked 差分と検査用ファイルを用意し、現行の起動引数で各 1 回試す。A は `git diff HEAD`、B は **`--output=<検査用ファイル>` だけ追加**。見る量はツールの実行可否・`permission_denials`・前後のファイルハッシュ。
→ A が実行でき、B でファイルが変われば仮説を支持する／A が実行でき、B が権限拒否されて不変なら、この CLI での当該経路の懸念は棄却する。ツール呼び出し未実施なら判定不能。

やらない方がよいこと: `touch` が拒否された実測を、そのまま「許可された git の全引数が読み取り専用」という保証に拡張すること。

呼び出し側の前提への異議: [plan:62](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:62) の「書き込みの道具は無い」は受け入れない。`Write`・`Edit` が無いことと、`Bash` 経由の書き込みが不可能なことは別である。最大の穴は許可外 Bash の扱いより、**許可内コマンドの副作用**にある。

不足情報: 上記 B の権限判定とハッシュ比較の実測。本レビューではファイル変更・子 CLI の起動は行っていない。
