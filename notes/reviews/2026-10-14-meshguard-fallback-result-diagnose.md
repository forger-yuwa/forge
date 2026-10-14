# codex 諮問 (diagnose): meshguard-fallback-result

- **brief**: [`notes/reviews/briefs/2026-10-14-meshguard-fallback-result.md`](../../notes/reviews/briefs/2026-10-14-meshguard-fallback-result.md)
- **plan**: [`plans/active/tooling-nozzle-mesh-size-guard.md`](../../plans/active/tooling-nozzle-mesh-size-guard.md)
- **date**: 2026-10-14
- **commit**: `7e7c529e` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 1.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-meshguard-fallback-result-diagnose.codex.log` (codex: ok)
- **結論**: 両 plan を「事前条件に照らして全件合格」として accepted にするのは保留し、最優先で fallback V1 の無応答を含む受入れ条件を固定した一度限りの再検証を行う。
- **extra**: `plans/active/tooling-review-fallback-backend.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 指摘 | 採否 | 理由・対案 |
|---|---|---|
| (a) V1 の timeout を再試行で埋めて合格とした | **採用・Major** | 事前条件は「各 3 回・6/6 正答」。timeout の除外・補充規則は記載されていない。実績は **8 試行、6 正答、2 無応答**であり、「応答が得られた 6 件は正答」と分けて記録すべき。根拠: [fallback plan:107](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:107)、[同:122](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:122)。 |
| (b) 既存 YAML 7 本の停止は不当な破壊 | **却下** | 上限超過時の停止は目的に合う。既存計算の破損を示す証拠はなく、確認できるのは再 prepare の互換性変更。対処は対象セッションへの通知と、了承に基づく明示的な上限指定。根拠: [mesh plan:144](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:144)。 |
| (c) `permission_denials` に残らないので V3 不合格 | **却下** | 条件は「禁止した道具が一覧に無い **または** 拒否記録がある」。報告された `[Glob, Grep, Read]` とハッシュ一致は前者を満たす。追加の拒否記録を合格要件にしない。根拠: [fallback plan:114](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:114)、[同:126](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:126)。 |
| (d) V1b 判定不能では mesh plan を閉じられない | **一部採用・Minor** | ni 切り離しを取り下げたため、追加計算は不要。ただし §8 の「V1b 合格」は事実と矛盾する。「切り離し不採用により追加検証不要」と完了条件を同期する。根拠: [mesh plan:112](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:112)、[同:151](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:151)。 |

結論: **両 plan を「事前条件に照らして全件合格」として accepted にするのは保留し、最優先で fallback V1 の無応答を含む受入れ条件を固定した一度限りの再検証を行う。**

第 1 仮説: **Major — fallback V1 の合格は、応答が得られた試行だけを数える条件付きの評価になっており、事前登録した 6 試行の合格を示していない。** 確度: 高  
  根拠: `tooling-review-fallback-backend.md:107` は各群 3 回を指定するが、同 `:122` は A3 の 600 s timeout 2 回と、その後の正答を記録する。同 `:50` の `unavailable` は運用上の分類であり、試験の分母から除外する規則ではない。  
  反証条件: 試験前に、無応答の除外・補充回数・打切り条件を定めていた記録が提示され、全試行がその規則に従っていたと確認できること。

第 2・第 3 仮説: 無し。timeout の原因や、モデルが同調しやすいという推測までは支持されない。

判別 A/B: 旧版の同じ抜粋・同じ問い・同じモデルと起動条件を固定し、**依頼側の誤った結論の有無だけ**を変える。A＝無し、B＝有りを独立に各 3 回、各 600 s 上限、補充なしで実施する。前書きを強めず、先行回答を参照できない条件で、正答・根拠行・B の明示的反論・無応答を全件保存する。  
→ A が成立し B が誤前提を採るなら、同調への耐性不足として不採用。  
→ 両群が期限内に全条件を満たすなら、**新しく明示した条件での合格**とする。無応答があれば受入れ未成立、そこで打ち切る。過去の試験を「事前条件で合格」へ書き換えない。

やらない方がよいこと: 正答が揃うまで A3 を補充すること、timeout を誤答や同調の証拠と断定すること、mesh V1b を合格にするためだけに CFD を延長すること。

呼び出し側の前提への異議: **「6 件の正答を得た」と「登録した 6 試行が合格した」は異なる。** 一方、mesh V1b は設計変更を見送る根拠として扱えるため、判定不能のまま終了する判断は妥当。ただし「全検証合格」という完了表現は修正が必要。

不足情報: V1 全 8 試行の固定プロンプト・起動条件・独立性を確認できる永続的な記録と、試験前の補充規則。今回確認したのは許可された plan 本文と差分統計であり、実装・生の試験記録・run の VERDICT は独立検証していない。
