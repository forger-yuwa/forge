# codex 諮問 (diagnose): deltastar-band-root-zone-result

- **brief**: [`notes/reviews/briefs/2026-10-13-deltastar-band-root-zone-result.md`](../../notes/reviews/briefs/2026-10-13-deltastar-band-root-zone-result.md)
- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **date**: 2026-10-13
- **commit**: `cd6a923b` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.3 min, rc=0
- **結論**: PR 前に、AWS の M6 生産で実際に使う NS/Euler の固定入力対で、`nonpos_ref` 検査の有無だけを変え、抽出から NS 補正入力の受入れまで通す A/B を 1 回行う。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: PR 前に、AWS の M6 生産で実際に使う NS/Euler の固定入力対で、`nonpos_ref` 検査の有無だけを変え、抽出から NS 補正入力の受入れまで通す A/B を 1 回行う。

第 1 仮説: **[Major・要再検証] `nonpos_ref` による hard 化が、生産入力では梯子内に留まらず、最終の補正表や寸法計算の停止へ伝播する可能性がある。**「欠陥の無い入力で同一」という結果だけでは、この変更の生産への影響を判定できない。確度: 中（生産での発現は未確認）
  根拠:
  - `design/forge_design/metrics/deltastar.py:277–278` は、D の符号分岐より前に非正・非有限の参照を NaN にする。`:589` で `nonpos_ref` は hard に分類される。
  - `design/forge_design/feedback/deltastar_loop.py:281–286` は、使用域の `held` が **1 点でも**あれば例外で止める。したがって、最終出力への伝播には運用上の影響がある。
  - `plans/active/tooling-nozzle-deltastar-band-root-zone.md:131` の記録では、手元の 12 対で各約 6,000 / 124,000 呼び出しが非正参照に該当する。ただし最終出力は不変。この事実は尊重するが、同ファイル `:140` の AWS 生産入力での確認は未実施である。
  反証条件: 対象の生産入力で、検査の有無によらず最終 CSV 2 本がバイト同一、使用域の hard/held がともに 0、NS 補正入力の受入れも両腕で成功する。この場合、この入力での伝播仮説を棄却する。

第 2・第 3 仮説: 追加しない。

判別 A/B: **A＝`nonpos_ref` 検査だけ無効、B＝現行の検査あり**。探索域修正・帯選択・平滑化・前回補正・保存済み設計壁を固定し、各腕で全断面を 1 回抽出する。forge は **0 step**。梯子／最終抽出別の該当件数、帯位置、最終 CSV 2 本、使用域の hard/held 数、`load_ns_correction` の成否を見る。
  → **両腕で最終出力・受入れが同一なら**、対象入力での伝播仮説を棄却する。
  → **B で最終表が変わる、または B だけ停止するなら**、「梯子内だけの変更で生産経路に影響しない」という仮説を棄却する。停止自体は修正の誤りとは限らず、旧抽出が不正な参照を通していた可能性を示す。

やらない方がよいこと: A の検査無効化を本番対策にすること、または V5 の未決着を今回のコード検証に持ち込んで長い NS 計算を追加すること。

呼び出し側の前提への異議: **V2 の合格から、生産入力への影響確認まで済んだとは言えない。**一方、M5 デモを合否対象外とする判断は、plan `:212` のユーザ決定に沿っており採用する。

不足情報: AWS 生産の正確な NS/Euler の run パス・抽出条件・上記 A/B の結果、および実行中とされた回帰 43 本の最終結果。今回の提案は **plan 未反映**（読み取り専用の諮問）。
