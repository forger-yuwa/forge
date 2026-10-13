# codex 諮問 (diagnose): mknot-gate-output

- **brief**: [`notes/reviews/briefs/2026-10-14-mknot-gate-output.md`](../../notes/reviews/briefs/2026-10-14-mknot-gate-output.md)
- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **date**: 2026-10-14
- **commit**: `fa8450f9` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-mknot-gate-output-diagnose.codex.log` (codex: ok)
- **結論**: カード組込み前に、境界を全て残す距離判定を通った `chosen_N` と、その値で最終検証した組だけを出力する契約を §4.3・§4.4 に確定する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 穴 | 採否・重大度 | 理由と対案 |
|---|---|---|
| 1 | **採用・Major（補足必須）** | 出力は移動量によらず `chosen_N`。再固定点の判定は `moved_N` 単独でなく、**固定点計算に実際に使った M_K との差**で行う。差 ≤ 0.01 でも、出力する組の寸法・制約を再確認する（[plan:89](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:89)）。 |
| 2 | **採用・Minor** | 入力点と再選択した2点を両解像度で検査する実装は、V2(g) の棄却条件に沿う。入力点の制約違反を再選択で隠さず、距離は出力候補について判定する（[mknot_select.py:625](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:625)）。 |
| 3 | **却下・Major（現方式）** | 同方向・距離 ≤ 0.1 による境界の統合は、安全距離を過大評価しうる。**当面は N・2N の境界を全て残して距離を測る**。片側だけの境界を残す処置は採用（[mknot_select.py:505](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:505)）。 |
| 4 | **採用・Minor** | 全境界について成立側へ ΔM_safe の候補を追加する。上側境界では減算する現実装が正しい。評価後も全境界に対する距離で絞る（[mknot_select.py:345](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:345)）。 |

結論: カード組込み前に、境界を全て残す距離判定を通った `chosen_N` と、その値で最終検証した組だけを出力する契約を §4.3・§4.4 に確定する。

第 1 仮説: **Major — ゲート合格後に元の `chosen` を出力すると、安全距離を満たさない値を出力しうる。** 確度: 高  
  根拠: [mknot_select.py:647](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:647) の距離判定は再選択した2点だけに適用され、元の入力点の距離は記録のみ。提示された合成例では、元の 2.050625 は距離 0.03、`chosen_N` の 2.070625 は距離 0.05。したがって元の値を出力する案は棄却する。  
  反証条件: 元の値を出力する経路でも、その値について最終設定・共通境界で再検証し、距離不足なら必ず出力を拒否すること。ただし、その追加経路を設ける必要はない。

第 2 仮説: **Major — 境界統合が距離を過大評価する。** 合成入力で再現済み、実データでの発生は未確認。N の不成立区間を (2.00, 2.02)、2N を (2.08, 2.10) とすると、2.05 は両解像度で成立するが各境界まで 0.03。現実装は `upper=2.00` と `lower=2.10` だけを残し、絶対距離を 0.05 とする（[mknot_select.py:329](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:329)、[同:524](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:524)）。

判別 A/B: **上記合成入力で `pair_tol` だけを A=0.1／B=0 に変更する。forge は0 step、関数評価を各1回。** 事前の判別基準は「Aだけ距離合格なら境界統合による過大評価、両方不合格ならこの反例は不成立」。許可範囲の実関数を抽出して実行した結果、**A: 距離 0.05・合格、B: 距離 0.03・不合格**だった。これは距離条件の再現であり、ゲート全体の実行結果ではない。

やらない方がよいこと: `geometry.M_knot` だけを `chosen_N` に差し替え、古い寸法検証・選択行・記録を流用すること。特に `row_chosen_N` は名称に反して**入力の `chosen` の評価行**を保持する（[mknot_select.py:674](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:674)）。

呼び出し側の前提への異議: 「保守側の端を選べば、絶対距離も保守的になる」は成立しない。また、`moved_N` はゲート入力との差であり、固定点計算に使った値との差とは限らない。例えば固定点後の再選択で +0.008、ゲートでさらに +0.008 なら、総変化 +0.016 を見て再固定点が必要になる。差 ≤ 0.01 は寸法検証の免除条件ではない。

不足情報: カード統合実装は閲覧対象外のため、最終検証・記録・YAML の一致は未確認。M5 の提示数値も独立検証していない。今回の出力契約の判断と境界統合の反例には影響しない。
