# codex レビュー: tooling-review-fallback-backend (result)

- **plan**: [`plans/active/tooling-review-fallback-backend.md`](../../plans/active/tooling-review-fallback-backend.md)
- **stage**: `result` (diff base `88675245`)
- **date**: 2026-10-14
- **commit**: `fc72a1df` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.7 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-tooling-review-fallback-backend-result.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C1/M0/m0

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

代替経路の基本設計は妥当で、V1 再試験の「6/6 正答・B 3/3 前提棄却」は保存データと一致しました。
ただし、判定の逆転と証拠不足の見逃しを再現しました。以下を直すまで `accepted` への移動は保留すべきです。

1. **Major — 実際の `NO-GO` を `GO` として記録できる。**

   根拠: [`codex_review.py:101`](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:101)、[抽出処理:327](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:327)。本文全体から最初の判定語を採っています。読み取りだけの関数実行で、次の入力が `ok / verdict=GO` になることを確認しました。

   ```text
   GO は不適切です。
   判定: NO-GO
   指摘数: Critical 1 / Major 0 / Minor 0
   ```

   この値が記録の見出しと §6.1 の行に使われます。既知事項でも、レビューの結論を逆転させるため持ち越せません。

   **対案:** 明示された判定行だけを抽出し、欠落・相反する複数判定は `incomplete` にする。説明中の判定語と最終判定が異なる回帰試験を追加し、§5.1 #5 を完了させる。

2. **Major — `diagnose` の指定書式に従った証拠不足が `ok` になる。**

   根拠: [出力テンプレート:209](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:209) は `不足情報:` を要求しますが、[分類処理:334](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:334) が検出するのは `証拠不足:` だけです。実行結果は次のとおりでした。

   ```text
   結論: 判断できません
   不足情報: 最終 HDF5 の NaN 検査結果
   → ok、conclusion=判断できません
   ```

   Claude 向け前書きには追加指示がありますが、codex 共通のテンプレートにはありません。「証拠不足は未成立」という §4.2 の契約を満たしていません。

   **対案:** 両 backend の出力契約を統一し、判定に必要な証拠の欠落を明示する欄を必須にする。既存の `不足情報:` も扱い、この入力が非 0・§6.1 行なしになることを試験する。新規残作業として §5.1 に追加する。

3. **Major — 必須の diff 取得に失敗しても result レビューを成立させられる。**

   根拠: [`build_evidence():463`](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:463) と [:475](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:475) は git の失敗を注記するだけです。[呼び出し側:588](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:588) は継続し、最終分類には証拠取得の成否を渡しません。

   ファイル操作をメモリに置換した確認で、git が `rc=128` を返しても関数は正常復帰し、空の `diff.patch` を作る経路を確認しました。その後の `GO` 本文は `ok` になります。親が把握済みの証拠欠落を、子の自己申告だけに委ねています。

   **対案:** 必須の git 取得失敗は親で未成立として終了する。意図した差分省略とは区別し、不正な `--base`・git timeout を負例に追加する。§5.1 に記載する。

4. **Minor — 独立した再試験の対策として `--out-dir` 変更だけでは足りない。**

   根拠: [plan:100](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:100) は「写しか別の `--out-dir`」を挙げています。しかし [起動処理:621](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:621) の `cwd` は常に元リポジトリで、先行記録への読み取りアクセスは残ります。対応する反復試験手順も手順書にありません。

   **対案:** 独立性が必要な反復では、先行回答・ログを除いた写しを使う、と手順を具体化する。`--out-dir` は出力先の分離に限定して説明する。今回の V1 再試験は `--tools ""` なので、この指摘で合格を取り消す必要はありません。

5. **Minor — 手順書と影響範囲の説明が実装と食い違う。**

   根拠: [`procedures/codex-review.md:23`](/home/sano/work/forge-deliv/procedures/codex-review.md:23) は diagnose を §6.1 に「書かない」としますが、[実装:903](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:903) は貼付行を表示し、本 plan にも記録されています。また [plan §7:151](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:151) の「codex が失敗したとき……回し直す」は、`unavailable` 限定の実装より広い説明です。

   **対案:** diagnose 行は plan/result の代用にはならないことを明記して記録方針を統一する。§7 は切り替え条件、終了コード、ログ名の変更まで現行実装に合わせる。

検証の裏付けとして、V1 の JSON 全6件、プロンプトの SHA-256、事前登録 commit と開始時刻を照合しました。V3 も一時領域に残るハッシュ表で、元の **9006 ファイルに変更・削除なし、追加4件は親の記録類**と確認できました。V4 の記録形式も一致しています。既存記録の分類は現在 **322/322 `ok`**、`check_plans.py` は `VERDICT: PASS` でした。ただし、この PASS は構造確認です。V2 全71件は一時ファイルを作成するため、今回の変更禁止条件では再実行していません。

指定 diff には別 plan のメッシュ規模制限も含まれますが、CUDA ソルバ本体の変更はありません。本 fallback plan に CFD の収束・一致の主張はなく、数値計算用 VERDICT を受入れ根拠に追加する必要はありません。

**推奨は、代替経路を維持したまま、上記1→2→3の成立判定を修正し、4→5の文書を同期してから移すことです。** 修正後の負例と V2 を通し、その結果・本レビューの採否を §5.1／§6.1 に反映してください。

指摘数: Critical 0 / Major 3 / Minor 2
