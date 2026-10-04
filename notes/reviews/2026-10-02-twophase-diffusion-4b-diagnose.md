# codex 諮問 (diagnose): twophase-diffusion-4b

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diffusion-4b.md`](../../notes/reviews/briefs/2026-10-02-twophase-diffusion-4b.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `22258d86` (feature/species-transport)
- **codex**: effort `high`, 7.3 min, rc=0
- **結論**: **初版を定常擬似時間に限定し、人工ソース S9 のパラメータ探索を打ち切って、実際の `condFloat` ソースと新更新を接続する最小検証へ進む。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 論点 | 重大度・採否 | 根拠と対案 |
|---|---|---|
| Q1：輸送の更新写像を除外 | **採用。ただし S8 の範囲限定** | 1000 更新を再実行した。判定から漏れていた総液量も含め、最大保存誤差は f32 停止で **A 4.38e−7／B 4.73e−7**、f64 停止で **A 6.96e−7／B 5.06e−7**。上限到達・状態補正はゼロ。輸送単独の当該懸念は閉じてよい。相変化・CUDA 実装まで除外したとは扱わない。 |
| Q1：蒸気残差の必須化 | **採用** | ブリーフの蒸気判定除外試験では両経路が不合格。`r_v = r_w − r_l` を独立した尺度で停止条件に含める。**Minor：判定コードも修正が必要**。[ハーネス:1081](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1081)の独立残差監査には `v` がなく、[1472](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1472)では保存集計から `rg` を除外している。輸送単独試験では液を必ず判定対象へ戻す。 |
| Q2：次の優先順位 | **Major／(a) を採用。ただし床の増大は未承認** | 最初の上限到達は液枯渇ではなく、**緩和された増分が格納丸めで消える状態**だった。下記の再測定では、同じ float32・同じ床で解ける。先に「許容床」と「反復写像の停滞」を区別する。`ρg_ref` の変更や保存許容の緩和を先行させない。 |
| Q3：2×2 前処理の必須化 | **Major／要再検証。現時点で要件化しない** | [ハーネス:1570](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1570)の対角比較は**緩和1・3更新だけ**。成功側は[1220](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1220)の空間結合込み **2N 系の double 求解**であり、局所2×2だけの検証でもない。実ソース・採用予定の緩和と前処理で必要性を判定する。 |
| 止めどころ | **定常擬似時間への限定を採用。「既存 limiter と同じだから受容」は却下** | [plan:88](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:88)は擬似反復ごとの保存を要求していない。S8 の輸送保存試験を維持し、相変化込み200物理更新の累積保存を初版の必須ゲートに追加しない。**Major：定常でも残差未達の停止は合格にしない**。既存経路も[推奨設定:172](/home/sano/work/forge-species/procedures/recommended-settings.md:172)で収束時の `condLim≈1`・補正ゼロを要求している。実ソース接続後の原方程式残差・補正監視を完了条件とし、dual-time 併用拒否は維持する。 |

結論: **初版を定常擬似時間に限定し、人工ソース S9 のパラメータ探索を打ち切って、実際の `condFloat` ソースと新更新を接続する最小検証へ進む。**

第 1 仮説: **S9 の少なくとも最初の上限到達は、緩和0.5と float32 格納丸めによる反復停滞である。** 確度: **高**

  根拠: [ハーネス:1227](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1227)で増分を緩和し、[978](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:978)で float32 状態へ加算する。変更なしのコードを追跡した結果、**物理更新 n=1（2回目）**で次を確認した。

  - θ は全セルで1、状態補正ゼロ。末尾の格納状態は変化しない。
  - 最大液残差のセル4は **ρg=1.32784e−4、g≈4.47e−4**。液枯渇近傍ではない。
  - 液増分 **−5.939e−12** に対し、格納液量の1 ULPは **1.455e−11**。増分は **0.408 ULP** で消える。
  - 液残差 **4.729e−11** は停止閾値 **3.577e−11** を超えたまま残る。

  反証条件: 同じ停止点で増分が0.5 ULPを十分超え、commit後も格納状態が継続して変化するなら、この停滞機構では説明できない。**残る全上限到達の原因までは特定していない。**

第 2 仮説: 実ソースでも連成前処理が有用である可能性はあるが、必須とは未確認。確度: **低**。S9 の簡易成長則と実際の核生成・モーメント依存ソースは異なる。

判別 A/B: **変更は対象更新の緩和係数だけ**。S9 の n=0 を緩和0.5で共通に計算し、その同一終状態から n=1 を A=0.5／B=1で再計算する。床6ε、ソース、前処理、limiterは固定し、上限3000回。格納値の変化、増分/ULP、原方程式残差を確認する。

  → **Aだけ停滞・Bが基準達成なら、この失敗を不可避な精度限界や液枯渇に帰す説明を棄却する。** Bも同じ残差床で停滞するなら、緩和0.5だけの説明は棄却する。

  **実施済み：Aは3000回で未達、Bは17回で基準達成。** 格納状態からEOSを再評価した独立残差比は **A 1.350／B 0.800**。これは局所的な判別であり、全200更新の合格や緩和1の一律採用を意味しない。

やらない方がよいこと: **`ρg_ref` を調整して合格させること、前処理対角で床を広げただけで原方程式を解いたと扱うこと、緩和1や2×2連成を一律の解決策にすること。**

呼び出し側の前提への異議:

- **Major：「float64で解けるので写像の問題ではない」は却下。** 有限精度の commit を含めて更新写像であり、今回まさにそこで止まった。
- **Major：S9 は実 `condFloat` の試験ではない。** [ハーネス:1067](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:1067)はソースをfloat64で評価する。実ソースは[condensationSourceF_d.cuh:121](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:121)で核生成・液滴モーメントに依存する。S9は人工問題の混合精度試験と記録し直すべき。
- 相変化込みS9のFAILは保持する。ただし、それを定常専用初版の無期限停止条件にする根拠はない。

不足情報: 実ソースと採用予定更新を接続した際の、残差・格納精度・モーメント補正・前処理の検証結果。CFDの収束・準定常性は今回評価していない。**ファイル変更なし、plan未反映**。呼び出し側で §4.2・§5.1 #4b・§6 に採否と定常限定の完了条件を反映すること。
