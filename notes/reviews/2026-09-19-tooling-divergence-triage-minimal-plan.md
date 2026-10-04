# codex レビュー: tooling-divergence-triage-minimal (v1 草案, plan 段)

- **plan**: [`plans/active/tooling-divergence-triage-minimal.md`](../../plans/active/tooling-divergence-triage-minimal.md)
  (本レビューの対象は現行 plan ではなく**その前身の草案 v1**。全文は末尾の付録に収めた)
- **stage**: `plan`
- **date**: 2026-09-19
- **codex**: `codex exec --sandbox read-only`, effort `high` (手動起動。`codex_review.py` と同じ作法:
  stdin `/dev/null`・stdout をパイプに通さない。plan がまだリポジトリ外だったため草案をプロンプトに直接貼った)
- **判定**: **NO-GO**, 指摘 C0/M11/m1
- **対象**: Jev (TypeSafe System One) を forge の開発・運用に組み込む草案 v1 (Jev クライアント層 + 発散トリアージを shadow mode で実装する案)
- **対応**: 全件採用し、順序を「証拠 → 規則 → 残余測定 → モデル」に反転した草案 v2 に改稿した。

Jev は 2026-09-15 早期アクセス開始の新しいモデルで codex の知識に無いため、公式ドキュメント
(`docs.typesafe.ai`) から起こした仕様・価格・弱点・信頼度運用を「所与の前提」としてプロンプトに添えた。

本文は codex の最終メッセージをそのまま転記 (指摘の採否は plan 側の §6.4 レビュー記録と §5.1 残作業表に書く)。

---

**判定: NO-GO**
現案のまま Jev クライアント・フック連携の実装には進むべきではありません。
診断支援の不足は実在しますが、Jev 固有の価値は未証明で、評価設計も採否を判断できる状態ではありません。まず決定論の診断要約と B1 を整備することを推奨します。

Jev の仕様は提示された事実を所与としました。以下の「草案 §…」は貼付された本文を指します。ファイル変更・外部 API 呼び出し・CFD 実行はしていません。

1. **Major — 埋めるべき隙間の同定より、外部モデル導入が先行しています。B1 が勝つ可能性は高いです。**

   **根拠:** [`check_solver_config.py:60`](/home/sano/work/forge/solver_density_cuda/tools/check_solver_config.py:60) は完全修飾パスで誤配置・未知キーを検出済みです。[`divergence-and-startup.md:10`](/home/sano/work/forge/procedures/divergence-and-startup.md:10) の主要チェックも、IC と入口の整合、出口逆流条件、初期 ω、`pRef` など、明示的な条件です。草案の「初手発散か」はモデルに聞く必要がありません。

   したがって、誤配置キー→正しい節、逆流条件不足→出口 BC 節、SST と非正の初期 ω→初期化節、という B1 は、該当する既知パターンでは確率モデルを介さず根拠付きで案内できます。ただし「uniform IC＋step 0/1 NaN→IC/BC 不整合」だけでは粗すぎます。`wall_dist=0` でも初手発散するためです。根拠は [`recommended-settings.md:71`](/home/sano/work/forge/procedures/recommended-settings.md:71)。

   **対案:** B1 を弱い文字列検索にせず、同じ構造化事実・既存 checker の結果を使う条件規則にしてください。既知規則で解決できない事例の頻度と調査時間を測ってから、モデルが必要か判断します。**B1 優位は現時点では構造からの予測であり、精度を実測した結論ではありません。**

2. **Major — A の価値はありますが、「誤答コスト低」は不適切です。最初の Jev 用途としても支持できません。**

   **根拠:** [`AGENTS.md:67`](/home/sano/work/forge/AGENTS.md:67) は原因と発生位置を切り分けてから対処することを要求しています。ところが草案 §4.4 は高信頼なら候補を一つだけ見せ、別の帯では追加 run を促します。誤誘導は GPU 時間、誤った設定変更、ソルバ不具合の見逃しにつながります。「読む順番だけ」ではありません。

   B〜E の評価も次のように修正すべきです。

   | 候補 | レビュー判断 |
   |---|---|
   | A | 診断支援として有用。ただし原因推定の誤答コストは中以上。まず決定論で実装する対象 |
   | B | 有用だが、既存 `check_solver_config.py` の拡張が自然。Jev の追加価値は小さい |
   | C | 意味解釈の隙間は実在。ただし直前の応答だけでは、引用された VERDICT や plan 反映の真偽を確認できない |
   | D | A より検証しやすいが、検索・見出し規則・既存エージェントに対する改善幅が小さい |
   | E | 不採用は妥当 |

   C の隙間は [`hook_convergence_gate.py:19`](/home/sano/work/forge/solver_density_cuda/tools/hook_convergence_gate.py:19) がファイル時刻だけを調べ、[`check_plans.py:4`](/home/sano/work/forge/solver_density_cuda/tools/check_plans.py:4) が内容の良否を見ないことから確認できます。

   **対案:** 最初の一本は **A の決定論版**にしてください。Jev の最初の用途を無理に選ぶ必要はありません。C も将来の自動判定には応答と証拠ファイルを対応付ける仕組みが先です。

3. **Major — 原因8択が、原因・症状・場所を混在させています。正解ラベル自体が一意になりません。**

   **根拠:** [`time_integration-update-positivity-guard.md:68`](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md:68) は、EOS 床 NaN を反復不安定の「終端症状」と明記しています。同じ事例が `cfl_too_high` と `eos_floor_washout` の両方に入り得ます。`axis_or_periodic_seam` は場所・境界構造であり、原因ではありません。

   また、[`boundary-node-inlet-corner-wall.md:16`](/home/sano/work/forge/plans/active/boundary-node-inlet-corner-wall.md:16) の入口・壁共有ノードの問題や、発散手順 §1(E) の free-stream 桁落ちには適切な専用分類がありません。未知原因も強制的に8択へ押し込まれます。`severity` の「設定変更／plan／バグ」も重大さの順序尺度ではありません。

   **対案:** 「観測事実」「確認すべき項目」「原因仮説」を分離します。仮説は複数選択・未知・証拠不足を許し、根本原因の確定と区別してください。学習・評価対象も、原因名より**次に行う確認の妥当性**に合わせる方が運用目的に近いです。

4. **Major — 数値を文章化するだけでは、弱点対策になりません。必要な情報を落とし、推論を前処理へ移しています。**

   **根拠:** [`divergence-and-startup.md:72`](/home/sano/work/forge/procedures/divergence-and-startup.md:72) は、最初の非有限値の位置を入口・出口・壁・内部に分けることを診断の中心にしています。しかし草案 §4.3 の5ブロックには、その情報を生成する具体的な契約がありません。メッシュ最悪位置と NaN 発生位置も別物です。

   [`check_convergence.py:39`](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:39) の読み込みは `outer_end` があれば他の phase を除外します。そのまま再利用して「最初の NaN」と呼ぶと、内反復中の先行異常を落とします。また [`check_mesh_quality.py:150`](/home/sano/work/forge/solver_density_cuda/tools/check_mesh_quality.py:150) は node メッシュでは primal セルを測り、同ファイル181行以降の出力には壁・軸・seam の位置分類がありません。

   **対案:** 自由文より、根拠付きの型付き事実を正本にします。最低限、観測 phase・検査間隔・欠測、NaN の領域、実効設定と既定値、メッシュ変換時設定、restart の出典、checker の版・引数を保持してください。未観測は `unknown` にします。数値比較は Python で済ませても、監査用の数値・単位・根拠を捨てる必要はありません。ケース名は分類入力から外し、匿名 ID と分離する方がよいです。

5. **Major — plan・notes・memory は、そのまま正解ラベル源にはできません。履歴の訂正と入力の事後変更が実在します。**

   **根拠:** [`architecture-axisym-axis-singularity.md:20`](/home/sano/work/forge/plans/accepted/architecture-axisym-axis-singularity.md:20) は、以前「真因」とした float 精度説を幾何不整合へ訂正しています。旧結論も本文に残っています。

   さらに [`config-key-pruning.md:477`](/home/sano/work/forge/plans/accepted/config-key-pruning.md:477) には **3,491 本の config を移送した**記録があります。現在置かれた config が、過去のログ生成時の入力と同一とは限りません。現在の solver ソースからキー仕様を抽出する checker も、旧バイナリの受理仕様とは一致しない可能性があります。

   **対案:** 各事例に実行時の config・バイナリ／commit・メッシュ・ログの対応を確立し、ラベルには訂正後の根拠箇所、確認介入、確定度を付けます。「変更後に回った」だけで因果確定としないこと。当時の入力を復元できない事例は、不確実事例として区別してください。memory は根拠探索の入口に留めます。

6. **Major — 母集団と分割方法に、選択バイアス・リーク・分布ずれがあります。**

   **根拠:** このツリーの保存済み判定1,679件を見出しで集計すると、`DIVERGED` 167件、`NOT CONVERGED` 1,448件、`PASS` 60件、分類不能4件でした。これは**保存された判定の集計で、再判定ではありません**。[`check_convergence.py:233`](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:233) は `NOT CONVERGED` に「まだ収束途中」も含めます。

   不明ラベルを全除外すれば、実運用で重要な未知原因・情報不足への対応を評価できません。同じ根本不具合の CFL 掃引や restart 派生を別事例として分割すると、独立性も失われます。さらに現行方針は [`recommended-settings.md:41`](/home/sano/work/forge/procedures/recommended-settings.md:41) の node 主体であり、歴史的な cell 事例を層化した比率が将来の利用比率とは限りません。

   **対案:** 調査エピソード・共通不具合・restart 系列を単位にまとめ、開発用・較正用・最終評価用を分離します。最終評価には将来側の期間を使い、未知・複合原因・収束途中・ログ欠損も残します。ラベルは事後証拠で確定してよいですが、モデル入力は**診断時点で利用可能だった情報だけ**に固定してください。`criteria` や B1 の調整にも最終評価事例を使わないことが必要です。

7. **Major — 60件では、提示された合格ラインを信頼できる形で確認できません。打ち切り条件も一致していません。**

   **根拠:** 草案 §6 の最小高信頼群は24件です。21件正解なら precision は87.5%で合格しますが、二項独立を仮定した95% Wilson 区間は約 **69–96%**です。精度85%以上を裏付ける証拠には弱すぎます。相関した派生 run があれば、さらに情報量が減ります。

   `+10 pt` は60件中6件差にすぎず、対応のある比較の不確実性を示していません。信頼度5ビンの単調性・ECE も少数標本では不安定です。同じ60件でしきい値調整と合否判定をすると楽観評価になります。また、表では `+10 pt` 必須なのに、打ち切りは「B1 を上回らない」場合だけで、`+1～9 pt` が宙に浮いています。

   **対案:** 60件は探索用に限定します。最終評価は独立事例で行い、対応のある精度差の区間、クラス別成績、棄権率と誤答率を示してください。しきい値は別集合で固定し、必要標本数は許容誤答率と検出したい改善幅から決めます。**高 confidence でも原因を確定扱いせず、必須チェックリストは省略しない**運用にします。

8. **Major — キャッシュ保存を、再現性の保証と混同しています。**

   **根拠:** 草案 §4.1 は `state+questions`、§4.5 は `state+questions+model` とキー定義が不一致です。後者でも `jev-latest` が同じ文字列のまま更新されれば、新旧モデルの結果が混在します。ログ保存は既存応答の再表示には使えますが、新規推論の再現を保証しません。

   実際、既存 checker にも判定モードや閾値があります。例は [`check_convergence.py:178`](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:178) の `--drop`、`--tail`、`--from-floor`。これらを落とした state は、同じ VERDICT 文字列でも意味が異なります。

   **対案:** 固定モデル版が利用可能かを導入条件として確認し、利用できなければ「保存応答の replay のみ再現可能」と明記します。入力・質問・要約器・文書ルール・閾値・実モデル版を版管理し、応答は上書きせず保存します。API 障害はモデルの「判定不能」と別の状態にしてください。

9. **Major — 想定したフック接続点が実装と合わず、設計ループにも依存が入り込みます。**

   **根拠:** [`run_case.sh:22`](/home/sano/work/forge/solver_density_cuda/tools/run_case.sh:22) に現在 `DIVERGED` 分岐はありません。checker の失敗を許容して出力し、最後は solver の終了コードを返します。

   Stop フックの `reason` は、VERDICT ファイルが未作成／古いときだけ生成されます。[`hook_convergence_gate.py:28`](/home/sano/work/forge/solver_density_cuda/tools/hook_convergence_gate.py:28)。**新しい DIVERGED 判定ファイルがある run は、その追加先を通りません。**

   加えて設計評価器は共通ラッパーを同期呼び出ししています。[`runner.py:242`](/home/sano/work/forge/design/forge_design/evaluate/runner.py:242)、[`runner_sern.py:605`](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:605)。自動組み込みすれば、目的関数に使わなくても設計ループに API 待ち・ログ・外部送信が入ります。5秒タイムアウトと複数回試行は、p95 1秒の保証にもなりません。

   **対案:** 本計画では同期フック連携を撤回し、診断は明示実行の CLI に限定します。将来連携する場合も、フックは既存のローカル結果だけを読む構成にし、終了コード・stdout の JSON・処理時間をモデルから隔離してください。回帰には欠損依存、壊れた応答／キャッシュ、書込不可、429/529、接続停滞、並行実行も必要です。

10. **Major — API キー管理とデータ送信管理が混同されています。shadow mode でも送信は発生します。**

   **根拠:** 草案 §7 自身が case 名・設定・残差を外部送信するとしています。一方、有効化条件は API キーの存在だけです。環境変数の継承は、そのプロジェクト・run の送信承認を意味しません。AWS では成果物を S3 へ回収する運用もあり、ログの保存範囲が広がります。[`tooling-cloud-gpu-env.md:27`](/home/sano/work/forge/plans/active/tooling-cloud-gpu-env.md:27)。

   なお、`.gitignore` 不足という指摘はしません。[`.gitignore:64`](/home/sano/work/forge/.gitignore:64) に `case/**/run_*/` があり、想定配置は対象です。ただし Git 除外は外部送信・キャッシュ保護ではありません。

   **対案:** キーと独立した明示的 opt-in、送信項目の許可リスト、ケース名・自由文・パスの除去を設けます。認証ヘッダはログ対象外とし、キャッシュの権限・保持期間・削除方法、提供者側の保存／学習利用条件を導入前に確定してください。オフライン試験は「無通信」を保証する必要があります。

11. **Major — 評価しているのは分類性能であり、目的である調査時間の削減ではありません。人手1–2時間も過小見積もりです。**

   **根拠:** 草案 §1・§4.2 の価値は人間の調査時間ですが、§6 は accuracy・API latency・料金だけです。上記の原因訂正や入力移送の確認を含め、60件を1–2時間で処理するには1件1–2分しかありません。

   B2 の「同じ state を一回」も分類器比較としては有効ですが、通常の既存エージェントは根拠ファイルを追加で読めます。実運用の代替案との比較が欠けています。また、価格が安くても、ラベル維持・要約器・質問文・API 障害対応の保守費は残ります。

   **対案:** 最初の10事例でラベル確定時間を実測してください。仮に1件5–15分なら60件で5–15時間、再確認や実装は別です。主要評価は、**B1＋根拠リンク＋必要時の既存エージェント**に対する、調査時間・不要な追加 run・誤変更の削減にします。API 単価は副次指標です。

12. **Minor — 導入・棄却の完了条件と、既存文書フローが揃っていません。**

   **根拠:** 草案 §8 は棄却でも `jev/` 実装と AGENTS 追記を要求する形です。§5 と §5.1 では API 実装／ベースライン作成の順序も一致していません。[`AGENTS.md:117`](/home/sano/work/forge/AGENTS.md:117) が求める現在仕様文書の整備も抜けています。

   §6.1 の「本セッションで実施」は草案の仮置きとしては理解できますが、正式記録としては [`check_plans.py:80`](/home/sano/work/forge/solver_density_cuda/tools/check_plans.py:80) の実在するレビュー記録要件を満たしません。

   **対案:** 探索終了・採用・棄却の完了条件を分け、棄却時に本番実装や利用規則を残さない形にします。正式化時にはレビュー記録と採否を保存し、採用時だけ仕様・運用文書を更新してください。

**推奨は一つです。Jev 導入を見送り、決定論の診断パッケージと B1 を先に作ってください。**

既存 checker の結果、実行時設定、最初に観測した異常、その証拠、対応する手順書をまとめ、規則で判断できない事例を既存エージェントへ渡します。この経路で残る調査負担を測ることが、追加モデルの必要性を判断する最短の方法です。現案の API クライアント・キャッシュ・フックを先に作る投資は支持しません。

変更禁止の依頼に従い、**plan 未反映・レビュー記録ファイル未作成**です。本回答がレビュー成果物です。

指摘数: Critical 0 / Major 11 / Minor 1


---

## 付録: レビュー対象の草案 v1 全文

現行 plan (v3) はこの草案を書き換えたものなので、上の指摘を読むための原文をここに残す。
節番号は本付録の見出しに対応する。

```markdown
# Jev (TypeSafe System One) を forge の開発・運用に組み込む

## メタ

- **area**: `tooling`
- **status**: `draft`
- **related_docs**:
  - `AGENTS.md` (VERDICT 必須ルール・plan 反映トリガ)
  - `procedures/divergence-and-startup.md` (発散原因の taxonomy)
  - `procedures/codex-review.md` (外部レビューの記録形式)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

forge の開発は「数値ゲート (決定論)」と「人間/LLM の判断」の 2 層でできている。前者は
`check_convergence.py` などで完全に自動化済みだが、後者 —— 発散したときにどの原因帯を先に疑うか、
この config は現行レシピと整合しているか、この応答は AGENTS.md のルールを踏んでいるか —— は毎回
エージェントか人間の注意に依存しており、取りこぼしが実際に起きている (rms_ro だけで「収束」と誤報告、
plan 未反映の方針決定、旧設定の混入)。

Jev は「文章を生成せず、型付きの質問に確率付きで答える」モデルで、70–500 ms・$0.042/MTok と
**ゲートに埋め込める速度と価格**を持つ。本計画は Jev を forge の判断層に**限定的に**導入し、
「決定論ゲートは正本のまま、その手前に高再現率のスクリーニングを 1 枚足す」状態を目指す。

完了時に得られる状態: 発散 run を投げると原因候補が確率付きで返り、どの手順書を読むべきかが即座に分かる。
かつ、その判定が過去 run の golden set で定量評価され、信頼度しきい値が根拠付きで決まっている。

## 2. スコープ

- **やる**
  - Jev クライアント層 (`solver_density_cuda/tools/jev/`): 認証・リトライ・キャッシュ・オフライン fallback・
    全リクエスト/レスポンスのログ。
  - 既存履歴 (2703 run / 1679 `CONVERGENCE_VERDICT.txt` / 2274 `forge_run.log` / 68 codex レビュー) から
    golden set を組む仕組み。
  - **ユースケース A「発散トリアージ」1 本**を shadow mode で実装・評価し、合否を判定する。
  - 合格した場合のみ advisory (助言) として `run_case.sh` / Stop フックに出す。
- **やらない**
  - ソルバ本体 (C++/CUDA) からの呼び出し。物理・数値の判断に Jev を入れることは**しない**。
  - `check_convergence.py` / `check_quasisteady.py` / `check_mesh_quality.py` の VERDICT を Jev で置換・
    補正すること (§4.1)。
  - 設計ループ (`design/forge_design/opt`) の目的関数・制約への組み込み (再現性が壊れる)。
  - 本計画では B〜E (§4.2) の横展開は行わない。A の評価結果を見てから別 plan にする。

## 3. 前提 (Jev の事実と、forge 側の既存ゲート)

Jev は 2026-09-15 早期アクセス開始、1.13 が 2026-09-18 リリースの新しいモデルである。公式ドキュメント
(`docs.typesafe.ai`) から、設計に効く事実だけを挙げる。

| 項目 | 事実 |
| --- | --- |
| エンドポイント | `POST https://api.typesafe.ai/v1/systemone`, `Authorization: Bearer`, `model: "jev-latest"` |
| リクエスト | `state` (string / object / array) + `questions` (id → 型付き質問のマップ)。複数質問を 1 リクエストで並列評価 |
| 質問型 | `noul` (真偽、`criteria: {true, false}` 任意) / `choice` (`criteria` 必須、最大 255 択) / `score` (順序付き水準の配列、2 以上) |
| レスポンス | `answers.<id>` に `noul|choice|score` の値 + `probabilities` + `confidence` (0–1)、`usage.input_tokens` |
| 制約 | テキストのみ (画像不可)、context 32k tokens、文字列生成をしない |
| 価格 | 入力 $0.042/MTok、出力無料 (補助金付きの可能性ありと明記) |
| レイテンシ | 70–500 ms |
| 較正 | RLCD で確率を較正。**集団としての較正であって個々の答えの正しさは保証しない** |
| 信頼度運用 | 公式推奨は 3 段: ≥0.7 自動実行 / 0.5–0.7 要確認 / <0.5 エスカレーション。結果の重大さでしきい値を変える |
| 弱点 (1.13 の jaggedness 文書) | **数え上げ・数値精度・日付比較が不得手**、間接推論 (hops) に弱い、無関係な文脈に引きずられる、指示を字義通りに読む、`P(true)+P(false)=1` を仮定してはいけない、説明 (理由) を返さない |

forge 側の既存ゲート (すべて決定論。本計画で触らない):
`check_convergence.py` (PASS / NOT CONVERGED / DIVERGED)、`check_quasisteady.py` (STEADY / DRIFTING /
OSCILLATING / TRANSIENT-UNSETTLED)、`check_mesh_quality.py` (AR・skew)、`check_solver_config.py`、
`check_plans.py`、`.claude/settings.json` の 3 フック (`hook_forge_guard.py` / `hook_convergence_gate.py` /
`hook_plan_todo_gate.py`)。

## 4. 設計方針

### 4.1 役割分担の原則 (これを破ると導入は失敗する)

1. **数値判定は Jev に渡さない**。Jev 自身が「数え上げ・数値精度・日付比較が不得手」と明記している。
   残差の低下桁数・NaN の有無・AR/skew・定常性はすべて既存の決定論ツールが出す。
2. **Jev の入力 state は「決定論コードが作った記述文」にする**。生の `residual_history.csv` を貼らず、
   前処理で「`rms_roOmega` が step 12 で NaN。その 3 step 前から `rms_roOmega` のみ 2 桁上昇、他列は横ばい」
   のような文に落としてから渡す (Jev の弱点を前処理で消す = §4.3)。
3. **Jev は決定を持たない**。返すのは「どの原因帯を先に見るか」「どの手順書を読むか」という**注意の配分**。
   AGENTS.md が要求する VERDICT の代替には**決してしない**。文書にも「Jev の出力は VERDICT ではない」と書く。
4. **再現性を壊さない**。Jev 応答はすべて `run_*/jev_triage.json` に生入力・生出力・モデル名・タイムスタンプ付きで
   保存し、キャッシュキー = state と questions のハッシュとする。ネットワーク不通・API エラー時は
   **黙って無効化して処理を続ける** (フックを落とさない)。

### 4.2 適用候補と優先順

| # | 候補 | 入力 state | 質問型 | 価値 | 誤答コスト | 判定 |
| --- | --- | --- | --- | --- | --- | --- |
| A | **発散トリアージ** (原因帯の当たり付け) | config 要約 + 残差の記述 + detectNaN 出力 + メッシュ品質サマリ | choice (原因 8 クラス) + noul 数個 | 高 (人間の時間が最も溶ける所) | 低 (見る順番が変わるだけ) | **Phase 1 で実装** |
| B | config スメル検査 | solverConfig/bcondConfig の要約 + 解析種別 | noul ×N (レシピ整合、旧設定混入、段階起動の要否) | 中 | 低〜中 | Phase 3 候補 |
| C | 応答ルール自己点検 (Stop フック補強) | 直前の応答テキスト | noul ×4 (VERDICT 引用なしで収束/一致を主張していないか、run パス明示があるか、plan 反映があるか 等) | 中〜高 | **中 (誤 block は実害)** → 高信頼のみ警告、block はしない | Phase 3 候補 |
| D | 文書ルーティング | 質問文 + 文書見出し一覧 | choice (procedures/plans/methods, 255 択以内) | 中 | 低 | Phase 3 候補 |
| E | 設計ループの run 健全性 | CFD 結果サマリ | score | 低 (決定論で足りる) | 高 (最適化結果が汚れる) | **やらない** |

A を選ぶ理由: (i) 入力が既に構造化されている、(ii) 誤答しても「読む順番」が変わるだけで結果は汚れない、
(iii) 正解ラベルが過去の plan 変更ログ・notes・memory に実在する、(iv) 効いたかどうかが体感で分かる。

### 4.3 state 設計 (Jev の弱点を前処理で消す)

`jev/state_builders.py` に決定論の要約器を置く。A の state は次の 5 ブロック (合計 2–4k tokens 目標、32k 上限に
対し十分な余裕を取る。無関係な文脈は Jev を引きずるので**削る方を優先**する):

1. `case` と離散化 (`node`/`cell`)、次元 (2D 平面 / 疑似 2D / 軸対称 / 3D)。
2. 設定要約: 対流スキーム・次数・リミッタ・陰解法の有無と `cfl`/`cfl_pseudo`/`implicitRelax`・
   乱流モデル・ガスモデル (CPG/TP/凝縮)・境界条件の種類一覧 (値ではなく種類と整合性の記述)。
3. IC の由来: uniform / 同一メッシュ restart / cross-mesh `interp_field` / `inletProfile` の有無。
4. 残差の記述文 (決定論生成): 最初に NaN になった step と列、その直前に上昇した列、列ごとのトレンド語
   (falling/flat/rising)。**数値の比較・桁数計算は Python 側で済ませ、結果の語だけ渡す**。
5. メッシュ品質 VERDICT と最悪 AR/skew の位置 (壁/軸/seam のどこか)。

質問 (1 リクエストにまとめる):

- `cause` (choice, 8 択): `ic_bc_mismatch` / `cfl_too_high` / `eos_floor_washout` / `mesh_quality` /
  `wall_function_or_walldist` / `axis_or_periodic_seam` / `species_or_condensation_stiffness` /
  `config_key_misplacement`。各 `criteria` は `procedures/divergence-and-startup.md` と memory の指紋を 1–2 文で書く。
- `is_step0_divergence` (noul): 初手発散か、しばらく回ってからか。
- `needs_staged_start` (noul): 段階起動 (1 次 → 2 次、層流 → 乱流) をやり直すべきか。
- `severity` (score, 3 水準): 設定変更で直る / plan を立てて調べる / ソルバのバグ疑い。

**`P(true)+P(false)=1` を仮定した派生量を作らない** (公式に禁止されている)。しきい値は生の
`noul` 値と `confidence` に対して引く。

### 4.4 信頼度ルーティング

公式推奨の 3 段を、forge の「誤答コスト」に合わせて次のように設定する (初期値。§6 で較正して更新)。

| confidence | 挙動 |
| --- | --- |
| ≥ 0.7 | 原因候補を 1 位だけ提示し、対応する手順書の該当節を指す |
| 0.5–0.7 | 上位 2 候補を確率付きで併記し「切り分け run を 1 本」促す |
| < 0.5 | 「Jev 判定不能」と出す (推測を出さない)。通常どおり `divergence-and-startup.md` のチェックリスト全体に回す |

### 4.5 実装アーキテクチャ

```
solver_density_cuda/tools/jev/
  client.py         # HTTP 1 本。429/529 は指数バックオフ、総リトライ上限 2、タイムアウト 5 s
  cache.py          # sha256(state+questions+model) → 応答 JSON。~/.cache/forge-jev/
  state_builders.py # 決定論の要約器 (A 用: build_divergence_state(run_dir))
  questions.py      # 質問定義 (criteria の文面。ここがプロンプト資産)
  triage.py         # CLI: python3 -m ... triage <run_dir> [--shadow] [--json]
  eval.py           # golden set に対する評価 (accuracy / high-conf precision / ECE / latency / cost)
```

- API キーは `FORGE_JEV_API_KEY` 環境変数のみ。**リポジトリに入れない**。未設定なら全機能を無効化して
  「Jev 無効 (キー未設定)」と 1 行出すだけにする (AWS / オフライン環境でも既存フローが壊れないこと)。
- `--shadow` では判定を**表示せず**ログにだけ残す (Phase 1 はこれで走らせて評価用データを貯める)。
- 出力は必ず `run_dir/jev_triage.json` に残す。`case/` 配下の run 成果物なので commit しない
  (AGENTS.md「コミット対象」)。

## 5. 実装ステップ

1. `jev/client.py` + `cache.py` + キー未設定時の no-op。疎通確認 1 件 (`usage.input_tokens` と課金額を実測)。
2. `jev/state_builders.py` の `build_divergence_state(run_dir)`。既存の `check_convergence.py` の系列読み込みを
   再利用し、**数値→記述文**の変換をここに閉じ込める。単体テストを 5 run 分の固定入力で置く。
3. `jev/questions.py` に A の質問セット。`criteria` の文面は `procedures/divergence-and-startup.md` と
   memory の「指紋」記述から起こす。
4. golden set 構築 (§6)。`eval.py` でベースラインと同時に評価。
5. 合格したら `triage.py` を `run_case.sh` の DIVERGED 時分岐に advisory として差す (shadow 解除)。
6. `AGENTS.md` / `procedures/` に「Jev の出力は VERDICT ではない」節を追記。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | API キー入手と疎通 | 早期アクセス登録。`usage` から実コストを 1 件で実測し、本 plan §3 の価格を確認 |
| 2 | golden set のラベル付け | DIVERGED run から層化サンプリング 60 件。ラベル源は plan 変更ログ / notes / memory。**人手 1–2 h が最初の関門** |
| 3 | state_builders の実装 | §4.3 の 5 ブロック。数値処理は全部ここ |
| 4 | ベースライン実装 | 最頻クラス + キーワード規則 (§6)。Jev より先に作る |
| 5 | 評価と採否判定 | §6 の合否ライン。落ちたら**打ち切って plan を archived へ** |
| 6 | advisory 組み込み | 合格時のみ。`run_case.sh` と `hook_convergence_gate.py` の reason 文に 1 行足す |
| 7 | 文書追記 | AGENTS.md の「Jev は VERDICT ではない」 |

## 6. 検証

**golden set**: `case/*/run_*` のうち `CONVERGENCE_VERDICT.txt` が DIVERGED または NOT CONVERGED の run から、
case・離散化・ガスモデルで層化して 60 件を抽出。正解ラベル (原因 8 クラス) は該当 plan の変更ログ・
`notes/investigations/` ・memory の記述から付け、判断がつかない run は golden set から落とす (無理にラベルしない)。
**ラベル付けは Jev の state を見る前に行う** (順序を逆にすると評価が汚れる)。

**ベースライン (Jev より先に作る)**:

- B0: 最頻クラス。
- B1: キーワード規則 (例: 「`step 0/1` で NaN かつ IC が uniform」→ `ic_bc_mismatch`)。
  **forge の発散は既に指紋が文書化されているので、B1 が強い可能性が高い。ここが本計画の本当の争点**。
- B2: 既存 LLM (Haiku 4.5) に同じ state と選択肢を 1 回投げる。

**判定基準** (Phase 1 の合否):

| 指標 | 合格ライン |
| --- | --- |
| top-1 accuracy | B1 を **+10 pt 以上**上回る |
| confidence ≥ 0.7 群の precision | ≥ 0.85、かつその群が全体の ≥ 40% を占める (カバレッジ) |
| 較正 | 信頼度ビン (0.5/0.6/0.7/0.8/0.9) の実測精度が単調増加。ECE ≤ 0.15 |
| レイテンシ | p95 ≤ 1 s (フックに入れられること) |
| コスト | 1 件 ≤ $0.005 (state 4k tokens 想定で $0.00017。桁が合うことの確認) |
| 対 B2 | 精度が同等以上、またはレイテンシ/コストが 1 桁以上有利で精度差 ≤ 5 pt |

**打ち切り条件 (kill criteria)**: 上記のうち top-1 が B1 を上回らない、または高信頼群の precision が 0.85 未満の
場合、**Phase 2 に進まず plan を `archived/` に送る**。「実装したから使う」はしない。

**回帰**: キー未設定・ネットワーク遮断の 2 条件で `run_case.sh` と 3 フックが従来どおり動くことを確認する
(Jev 導入でフックが落ちるのが最大の実害)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-19` | (本セッションで実施) | — | — |

## 7. 影響範囲

- 新規: `solver_density_cuda/tools/jev/` (6 ファイル)。既存ソルバ・ツールの挙動は不変。
- 変更 (Phase 2 で合格した場合のみ): `run_case.sh` に advisory 1 行、`hook_convergence_gate.py` の reason に 1 行、
  `AGENTS.md` に「Jev は VERDICT ではない」節、`procedures/divergence-and-startup.md` に triage の使い方。
- 秘密情報: API キーは環境変数のみ。`.gitignore` に `jev_triage.json` は不要 (`case/` 成果物は元々 commit しない)。
- 外部送信: state には case 名・設定・残差の記述が含まれる。**外部 API に送る**ことをユーザが承認していること。

## 8. 完了条件

- [ ] `jev/` 実装と単体テスト
- [ ] golden set 60 件 + ベースライン B0/B1/B2 の評価表
- [ ] §6 の合否ライン到達、または打ち切り判断の記録
- [ ] キー未設定・オフラインでの回帰確認
- [ ] codex レビュー (`plan` / `result`) を §6.1 に記録
- [ ] `AGENTS.md` への「Jev は VERDICT ではない」追記

## 9. 変更ログ

- `2026-09-19` — 初稿。
```
