# モデル分担: 次のセッションでのチェック手順 (2026-09-22 作成)

正本は [`AGENTS.md`](../../AGENTS.md)「モデル分担とエスカレーション」。ここは**初回運用の確認手順だけ**の使い捨てメモ。
確認が済んだ項目は結果を 1 行書き足し、全部済んだらこのファイルは消してよい。
経緯と codex レビューの採否は [`notes/reviews/2026-09-22-codex-model-tiering.md`](../reviews/2026-09-22-codex-model-tiering.md)。

## 0. 起動直後 (1 分)

| # | 確認 | やり方 | 期待 | 違ったら |
| --- | --- | --- | --- | --- |
| 0-1 | 主セッションが Opus | `/model` または `/status` | Opus 5 | `.claude/settings.local.json` の `"model": "opus"` を確認。設計・診断セッションなら意図的に `claude --model claude-fable-5-1` で起動 |
| 0-2 | 3 本が読み込まれている | `/agents` | `run-watcher` / `implementer` / `diagnostician` が出る | `.claude/agents/*.md` の frontmatter を確認。定義を直したら**セッションを立ち上げ直す** (同じセッションからは呼べない) |
| 0-3 | 収束ゲートのフックが生きている | `/hooks` | PreToolUse(Bash)・Stop・PostToolUse(Write\|Edit) の 3 本 | `.claude/settings.json` を確認 |

## 1. 最初の run で見る (下向きの委譲)

最初に forge を回す作業で、次が**指示しなくても**起きるかを見る。起きなければ「run-watcher に prepare させて」と明示し、
どこで起きなかったかを下の記録欄に書く (AGENTS.md の文言を強める材料になる)。

1. **prepare**: 主セッションが run 名 (`run_NNNN_<slug>`) を決めてから `run-watcher` を呼んだか。
   返却が `PREPARED` で、差分・MESH VERDICT (メッシュを変えた場合)・起動コマンドが入っているか。
   - 複製元に `wallTreatmentSST: 1` があった場合に、落とさず `FAILED` で返してきたか。
2. **起動**: 長時間 run を**主セッション自身が** `run_case.sh` を `run_in_background` で起動したか
   (`run-watcher` に完走まで待たせていないか)。応答に run パスが明示されているか。
3. **early-check**: 序盤の NaN 確認が返ってきたか。
4. **post-check**: 完了通知のあと `run-watcher` が呼ばれ、返却が `CHECKED` で、
   CONVERGENCE VERDICT (判定区間つき)・指定した量の QUASISTEADY VERDICT・`residual_history.png`・README の 1 行が揃ったか。
   - `CHECKED` は「実施した」の意味。VERDICT の中身は自分で読む。
   - **主セッションが、VERDICT が揃う前に結果を報告していないか**。
5. **Stop ゲートが子の run を見ているか** (1 回だけ確認すればよい):
   ```bash
   cat ~/.cache/forge-gate-claims/<このセッションの session_id>.txt   # 回した case 名が入っているか
   ```
   入っていなければ、サブエージェントの `session_id` の挙動が変わっている。AGENTS.md に書いた測り直し
   (`echo "run_case.sh case/99.claimprobe_child"` を子に打たせて台帳を見る) を行う。

## 2. エスカレーションが発火するか (上向き)

自動では発火しにくい側。次の 7 場面のどれかに当たったとき、応答に「`diagnostician` に諮った (結論 1 行)」が
書かれているかを見る。書かれずに先へ進んでいたら、その場で「diagnostician に聞いて」と止める。

1. plan §4・§6 を書く / 方針を変える
2. 発散対処を 2 回 (case README の run 一覧で 2 行) やって解けない
3. plan §6 に事前に書いた参照値・許容差との比較が FAIL・判定不能 (事前に書いていない比較も)
4. plan・手順に無い修正や再実行へ進む前、「原因は○○」と書く前
5. codex の Critical / Major の採否
6. `solver_density_cuda/cuda_forge/` の数値の振る舞いを変える編集の前
7. result 段の解釈を確定する前 (予想どおりに見えるときも)

諮ったときのブリーフも見る: **観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説** が分かれているか、
各試行に run パス・commit・設定差分・判定区間・VERDICT が付いているか。`diagnostician` の返答に
「呼び出し側の前提への異議」の欄があるか。

## 3. 費用 (1〜2 週間後)

- セッションの種類ごとに `/cost` を控える: (a) run 中心、(b) 実装中心、(c) 設計・診断 (Fable)。
- 比べるもの: 切り替え前の同種セッションの費用、`diagnostician` の呼び出し回数
  (1 回あたり固定で約 4.7 万トークン。細かく何度も呼んでいたら、判断単位にまとめるよう文言を強める)。
- 減っていなければ、どの種類のセッションが食っているかを見て分担を見直す。

## 記録欄

| 日付 | 項目 | 結果 |
| --- | --- | --- |
| 2026-09-22 | 0-2 (別プロセスで事前確認) | 3 本とも読み込み、Sonnet 5 / Opus 5 (1M) / Fable 5.1 で動作 |
| 2026-09-22 | 1-5 (事前確認) | 親・子とも同じ `session_id` の台帳に claim された |
| 2026-09-22 | 1-3 (read-only 試験) | 形式・境界は守った。終了済み run に `RUNNING` を付けた → 定義を修正済み |
| 2026-09-22 | 0-1 / 0-2 / 0-3 (新セッション実地) | 主 = Opus 5 (`settings.local.json` `"model": "opus"`)、3 本とも呼び出し可 (run-watcher=sonnet / implementer=opus / diagnostician=claude-fable-5-1)、フック 3 本 (PreToolUse Bash / Stop / PostToolUse Write\|Edit) 登録済み |
