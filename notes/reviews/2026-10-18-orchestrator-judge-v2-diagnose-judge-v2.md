# codex 諮問 (diagnose): orchestrator-judge-v2

- **brief**: [`notes/reviews/briefs/2026-10-18-orchestrator-judge-v2.md`](../../notes/reviews/briefs/2026-10-18-orchestrator-judge-v2.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-18
- **commit**: `b678772d` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-18-orchestrator-judge-v2-diagnose-judge-v2.codex.log` (codex: ok)
- **結論**: 停止時の入力まで照合する共通処理と、準定常の仕様検査を必須にする形へ plan を修正し、GPU 不要の入力変更 A/B を先に通す。
- **extra**: `plans/active/tooling-convergence-and-wall-resolution-gates.md`, `plans/accepted/tooling-settle-convergence-and-monitor.md`, `solver_density_cuda/tools/judge_run.py`, `solver_density_cuda/tools/check_settled.py`, `design/forge_design/evaluate/judge.py`, `design/forge_design/feedback/wt_qualify.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。**案 A は修正して採用**を推奨します。特に (2)・(5)・(6) は現案のままでは不十分です。ファイル変更・forge 起動は行っていません。

| 案 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| (1) 判定器の選択 | **採用** | 現在の `judge_run.py:856` は旧判定器を呼びます。`settle.py:1972` の選択規則へ統一し、終了コードだけでなく原判定・道具・区間・停止 step を保持してください。`convergence: PASS` は「設定された量・区間について SETTLED」と表示すればよく、残差の低下を意味するとは説明しないこと。 |
| (2) 停止時との照合 | **却下〔Major〕—合否の比較だけでは不足** | `stop_rule.py:1527` の `stop_verdict_record` は記録そのものの同一性を確かめます。現在の入力・起動・終端・量の数値まで照合する既存処理は **`match_reeval` (`:1626`)** です。これと停止状態・`STOPPED_AT` の照合を再利用してください。`reeval` が欠けた v2 の終了済み run を照合対象から外してはいけません。また、再判定の `DIVERGED` 等を照合不一致の `INDETERMINATE` で消さず、原判定の FAIL と照合不能を併記してください。既知の FAIL を保持する契約は `judge_run.py:1308` にあります。 |
| (3) 準定常判定の共用 | **採用。ただし (6) の無条件省略は不可〔Major〕** | 1 回の `check_settled` の結果を両項目で共有して構いません。ただし、**数値判定の共用と、仕様の完全性検査の省略は別**です。名前の欠け・余りで `quasisteady` が判定不能になっても、`convergence` 単独は SETTLED になり得ます。`wt_qualify` はこの不一致を必ず拒否してください。基準の sha256・量の種類・区間も参照で結び付け、別の閾値を持たせない設計を推奨します。仕様未設定を拒否する現契約は `judge_run.py:1099`。 |
| (4) `FAILED` | **採用〔Major：現実装の欠落〕** | `judge_run.py:151` に無く、関数を実行しても `FAILED → INDETERMINATE` でした。v2 の確定状態として FAIL に対応させ、途中の記録を確定と扱わない検査 (`:924`) にも反映してください。 |
| (5) 鮮度 | **却下〔Major〕—入力一覧の取得方法を修正** | 実行確認した **`check_settled --no-write --json` には `files` がありません**。一覧は `settle.py:1901` の書き込み処理内で別の dict に追加され、標準出力の JSON には入りません (`check_settled.py:262`, `:273`)。判定結果と入力一覧を作る共通処理を、書き込みから分離してください。鮮度確認時は保存済み hash を再利用せず現在のファイルを再計算し、欠損から出現した入力も検出すること。旧形式では実際に選ばれた `SETTLE_VERDICT.json` も依存入力です。 |
| (6) `wt_qualify` への統合 | **却下〔Major〕—固定版で再評価し、項目ごとの扱いを明示** | 現在の `approval_rechecks` は整定とメッシュ品質だけです (`wt_qualify.py:851`)。追加する非有限・物理性・solver 終了等の保証には不足します。**承認モードでは固定版の `judge_run` 全体を別出力先で再評価する方式**を選んでください。`check_fresh` は道具の hash が変わっても警告付きで FRESH を返します (`judge_run.py:1519`)。FRESH は承認した判定器で評価した証拠にはなりません。 |

(6) の壁解像については、**総合 FAIL を保ったまま `back_to_correction` に送ること自体は妥当**です。これは合格への迂回ではなく、失敗理由から次の処置を決める操作です。ただし、次を plan §4.2 に先に明記してください。

- 必須の健全性・整定・停止時照合・準定常の仕様検査が通り、壁解像だけが FAIL：`back_to_correction / mesh_first_layer`。解析全体は不合格のまま。
- 必須判定の欠落・照合不能・物理性等の FAIL：`incomplete`。
- `accepted`：壁解像を含め必要な検査がすべて合格。

現 plan の「解析の成立が通ったときだけ要求の達成を評価する」(`plans/active/tooling-design-orchestrator.md:209`) には、**失敗原因の診断と補正先の選択は許す**という整理が必要です。壁解像を別途評価する場合も、step・全壁・閾値を揃え、二つの判定が食い違ったら無視しないでください。また、対象の NS では `solver / nonfinite / physical / convergence / stop_rule` の一律の `NOT_APPLICABLE` 許容は不要です。

確認した証拠は以下です。これは整定の再評価結果であり、`RUN_VERDICT` 総合合格の報告ではありません。

| run | 書き込みなしの再評価 | 判定区間・停止観測 | 停止時との照合 |
|---|---|---|---|
| `case/71.m55_combustion_wt/run_0011_ns_v2/` | `VERDICT: SETTLED (閾値は暫定 (plan §6 V3 で決める))`、exit 0 | C: 36,000–46,000、P: 26,000–36,000、停止 46,202 | `match_reeval.matches = true` |
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 同上、exit 0 | C: 36,000–46,000、P: 26,000–36,000、停止 46,171 | `match_reeval.matches = true` |

両 run とも `JUDGE_SPEC.json`・`RUN_VERDICT.json` はありません。既存 run の検証では、後付けの仕様を使った移行試験と、起動前に仕様を確定する生成器の試験を区別してください。

結論: **停止時の入力まで照合する共通処理と、準定常の仕様検査を必須にする形へ plan を修正し、GPU 不要の入力変更 A/B を先に通す。**

第 1 仮説: 案 A の「記録の同一性＋合否の一致」だけでは、停止後に入力が変わっても同じ合否なら検出できない。確度: 高  
　根拠: `stop_rule.py:1639` と `:1654` は別の検査である。`run_0014_chain_v2` を使い、ファイルを変更せず hash の戻り値だけを変えたメモリ内負例では、`stop_verdict_record.exact = true`・合否一致のまま、`match_reeval.matches = false`、理由は `boundary_massflux.csv` の hash 不一致だった。  
　反証条件: 案 A の実装が合否比較とは独立に、停止時の全入力と現在の入力を照合し、この負例を拒否すること。

第 2 仮説: `quasisteady` を無条件に読まない設計では、仕様と監視量の不一致を見逃す。確度: 高。提案設計上の経路であり、未実装のため通し再現は未実施。

判別 A/B: 実装後、同じ合格 fixture の独立した写しを使い、A は変更なし、B は **`boundary_massflux.csv` のコメントを 1 文字だけ変更**する。forge は回さず判定を各 1 回行う。数値と SETTLED は両者で同じになる。→ A が PASS・B が入力照合不一致で拒否なら欠陥解消／B まで PASS なら第 1 仮説を支持。A が拒否されたら陽性対照の不成立。

やらない方がよいこと: 旧残差基準へ戻す、閾値を緩める、停止時の原本を再生成する、`quasisteady` の仕様不一致を「二重判定だから」で捨てること。

呼び出し側の前提への異議: 「1 回の数値計算を共用するから `quasisteady` を確認しなくてよい」は成立しない。また、停止規則の method v2 とプリセット `wt_nozzle_v2` は別であり、現行 `wt_nozzle_v3` からユーザ決定で外した残差・壁熱流を再導入しないこと (`plans/accepted/tooling-settle-convergence-and-monitor.md:80`)。

不足情報: 設計判断を妨げる不足は無し。修正後の実装については、上記 A/B と仕様欠落・壁解像 FAIL の分岐試験が必要。

証拠: 足りる
