# 諮問: ⑤ SERN 設計チェーンの次の一手 — farfield・種 DB・main 統合が片付いた後の優先順 (2026-10-05)

関連 plan: `plans/active/tooling-nozzle-sern-chain.md` (§5.1 残作業表、§8 未確定事項)、`plans/active/tooling-nozzle-sern-3d.md` (§5.1、§4.46 = 現行の 3D 基準列)、
`plans/active/tooling-sern-mesh-blocking.md` (R5q の移管先、全ヘキサ接続模型)、`plans/accepted/boundary-node-farfield-characteristic.md` (2026-10-04 限定受理)。
メモリ相当の前提: ブランチ `feature/sern-design` = main 77318d0e + d22f4f04 (2026-10-05)。AWS g5 で run。

## 1. ここまでの到達点 (観測事実、各 plan に記録済み)
- 2D 評価チェーン: 評価ゲート修復 (R1)、力の集計分離 (R2)、凍結 TP (R3)、外部領域 (R4)、起動レシピ検証 (R-b)。生産 YAML = 種ごとの輸送 (R8)、種 DB 段 3 + #14 (R9/R9b)、species-transport 追加分 (R10)、main 統合 (R11) の回帰はすべて 2D m6_on でノイズ床内。
- 3D 基準列 (§4.46、2026-09-27): g1 `run_0973`・g3 `run_0971`・g4 `run_0972` (壁解像、chi 既定、lsq、出口 outflow)。g3→g4 の Δ は §8 許容内 (C_M −0.040 / 0.05)。残差はプラトー (ユーザ決定 §4.39 で受理の阻害要因にしない)。
- 側方境界: 2026-10-01 に生産 3D の side_far を farfield に (ユーザ決定)。新輸送・#14 後で必要幅 2.50 H、g3/g4 の G + D は §8 内 (限定付き、farfield plan #4f/#4h)。
- カウル後縁の低温スポット (R5h): g3 105 K・g4 95 K の 1 節点 (厚さ 0 の後縁ノード)。conv 1→0 の A/B と局所帳簿で「再構成依存を支持、EOS 単独原因は除外」(2026-09-27)、その後は未着手。
- R5q (接合部トポロジ) は全ヘキサ接続模型の plan に移管 (模型は PASS、生産規模 1250 万節点は載らない)。現行メッシャには `sz` テーパが残るが格子列内で一定。
- R6: (a)(b)(c)(e) 完了、(d) の許容値明記が残る。
- R7 (小規模探索で判別能力を確認してから MOO 再取得) は未着手。R1–R6 の後という順序 (codex 2026-09-09)。
- 3D の他の残: R4e-次 (W2 面単位フォールバックの A/B、別 plan `convection-node-wall-reconstruction.md`)、R5m (AR を下げる費用の測定)、R4c (旧トポロジ作り直し、R5q と重なる?)。

## 2. 問い
1. R7 に進むための残りの前提は何か (R6 (d)、R5h、R4e-次、R5m、R4c のどれが R7 の前に要るか、どれは R7 と並行・後回しでよいか)。
2. R7 の「小規模探索」の最小設計 (2D か 3D か、評価本数、判別能力の判定基準を事前にどう決めるか)。MOO は 2D 評価で回し 3D は代表点の確認に使う、という現行の役割分担でよいか。
3. 1・2 を plan §5.1 にどう書き直すべきか (完了済みだが取り消し線の無い行の整理を含む)。
