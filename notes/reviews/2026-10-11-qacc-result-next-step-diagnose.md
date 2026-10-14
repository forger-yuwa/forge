# codex 諮問 (diagnose): qacc-result-next-step

- **brief**: [`notes/reviews/briefs/2026-10-11-qacc-result-next-step.md`](../../notes/reviews/briefs/2026-10-11-qacc-result-next-step.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `8f7a56aa` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.5 min, rc=0
- **結論**: **同じ FP64 バイナリで、step 末尾の 9 保存量の float32 量子化だけを ON/OFF する、反復付き A/B を一つ行う。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/qa_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | §6.35 から平均流の丸めを原因候補から外す | **却下**。`qAccumulatorFP64` は累積の正本を double にするが、演算へ渡す `Q` は引き続き float。したがって、平均流の状態表現・残差入力の丸めまで除去した試験ではない。[implementation.md:248](/home/sano/work/forge-integ-1005/methods/time_integration/implementation.md:248)。記載は「今回の介入では大幅な抑制を観測しなかった。正式判定は判別不能」に限定する。 |
| **Major** | 「3 つの介入がどれも効かなかった」 | **却下**。§6.29・§6.33・§6.35 はいずれも総合判定が判別不能。特に SST 凍結では効果量条件も一部未達だった。[plan:760](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:760)。対案は、観測した差と、未確定の必要性・寄与を分けて記録すること。 |
| **Major** | step 末尾の丸めを「各 commit の保存精度だけを float にする」と呼ぶ | **要再検証**。平均流 commit と末尾の間には SST・化学種などの処理があり、末尾直前には等温壁ピンもある。[update_d.cu:45](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:45)、[main.cpp:2533](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2533)。途中で丸める演算と、すべて処理してから丸める演算は同じではない。**「step 末尾の 9 保存量の量子化の十分性試験」への限定を条件に採用**する。 |
| **Major** | G が離れれば丸めが「必要」、離れなければ原因は「演算側」 | **却下**。前者で示せるのは、この量子化介入が差を生む**十分性**。後者で棄却できるのは、この位置・頻度の量子化だけで差を生む仮説であり、途中の丸めや演算との相互作用は残る。§6.34 自体も二択を採用していない。[plan:859](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:859)。 |
| **Major** | B・B′・C・C′・E・E′・F をまとめた再実行差 | **却下**。異なる介入・実行環境の差を同条件の再実行差として扱えない。保存 JSON から再計算すると、N_E/N_B は **2.27〜4.43**。等分散も確認されていない。対案は、新しい対照・介入をそれぞれ反復し、腕内の分散だけを使うこと。 |
| **Minor** | 「FP64 ビルドなら線形の解きも double」 | **訂正して採用**。`implicitSolvePrecision: 0` では FP64 ビルドにも `ST=float` の経路がある。[plan:72](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:72)。今回の比較では設定を変えず、「既存の FP64 基準の演算経路を維持」と記す。 |

結論: **同じ FP64 バイナリで、step 末尾の 9 保存量の float32 量子化だけを ON/OFF する、反復付き A/B を一つ行う。**

第 1 仮説: **step 末尾の量子化だけでは観測規模の偏りを再現せず、step 内部の低精度演算・途中の丸めとの組合せが必要である。** 確度: **低**。  
　根拠: `case/45.isobutane_m6_d155/run_0592_qa_f32/`・`run_0593_qa_f32b/` の保存集計では、平均流の累積を double に変えても、A からの差の平均は基準 float の **97.68〜103.56 %**。ただし、これは上記仮説の傍証にとどまり、9 保存量の末尾量子化を直接検証していない。  
　反証条件: 下記 A/B で、量子化による差の区間下限が全 6 組で D の向きに **0.5|D| 以上**となること。

第 2 仮説: **9 保存量の末尾量子化だけでも、同方向・同規模の偏りを生む。** 未確認。平均流 5 本だけの累積変更では、ρk・ρω・ρY の影響を除けない。

判別 A/B:

- **対照 A\*** = 新しい FP64 バイナリ、量子化 OFF。**介入 G** = 同じバイナリ、量子化 ON。各 **3 回、30,000 step 固定**。H・I は今回は作らない。
- 共通入力・格子・実効設定・実行環境を揃え、各反復は同じ起点から独立に開始する。旧 A は参考資料とし、新しい A\* の反復を省略しない。
- ON は `applyNodeIsothermalWallPin` の後、`updateVariablesOuter` の前に、9 保存量を明示的に float32 へ丸めて double に戻す。`updateVariablesOuter_d` 本体のコピー対象は 7 本なので、**ρY0・ρY1 とそれらの履歴コピーも含むことを別途確認**する。[update_d.cu:63](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:63)。
- 投入前に、同じ入力配列に対する ON/OFF の単体照合で、丸め後の 9 本、履歴への反映、次の利用前の原始量・組成の整合を検証する。丸め以外の補正を新たに混ぜない。

**判定式は、再実行差に対する倍率試験から、事前指定した効果量に対する区間判定へ変える。** 過去の判定は変更しない。

各断面・各窓について、従来どおり起点からの変化 Δ を求め、

`e = sign(D)·(mean(ΔG) − mean(ΔA*))`

とする。独立な標本は **run ごとの窓平均**であり、窓内 11 出力や隣接する W1・W2 を反復数に数えない。

新試験では腕別の標本分散から Welch 型の差の区間を作り、6 組について Bonferroni 補正した全体 95 % の区間 `[L, U]` を事前登録する。これは独立性・run 平均の近似正規性を仮定するため、3 反復で精密な分布推定ができるとは扱わない。

- **結果 A：全 6 組で L ≥ 0.5|D|**  
  → 第 1 仮説を退け、**この量子化介入だけで、観測差の少なくとも半分の規模を生める**とする。必要性や実際の float の真因は確定しない。
- **結果 B：全 6 組で [L,U] が [−0.1|D|,+0.1|D|] に収まる**  
  → 第 2 仮説を退ける。次に調べる対象は途中の丸め・演算との組合せとなる。
- それ以外は判別不能。**有意差が出ないだけでは結果 B にしない。** 区間が広ければ情報不足であり、倍率や窓を変更して救済しない。

旧 A がほぼ動かないという前提も、新しい A\* で確認する。有限性、初期化後の保存量、介入の実効性を前提ゲートにし、`check_convergence`・`check_quasisteady` は対象区間と閾値を添えて記録する。Q_w は独立した副判定のままとする。

やらない方がよいこと: **§6.35 を結果 R に読み替えること、10 を結果に合わせて下げること、異なる介入の run をプールすること、最初から H・I に分岐すること。** `qAccumulatorFP64` の削除・既定変更を支持する証拠もない。

呼び出し側の前提への異議: §6.35 に書いてよいのは、次の範囲である。

> 登録された W1・W2 では、E・E′ の A からの差の平均は基準 float の 97.68〜103.56 % だった。大幅な抑制は観測されなかったが、事前登録判定は S 0/6・R 0/6 で判別不能。平均流の丸めの因果的関与は除外できない。

これは**有限期間の応答の記述**であり、収束解の一致や固定点の偏りの確定ではない。保存された収束判定は E・E′・F とも **`OVERALL: CHECK FAILURES ABOVE`**。また、F の合格範囲は B・B′ の生の最小〜最大ではなく、左右に N_B を加えた許容範囲である。

不足情報: このワークツリーには対象 run 本体、`tr_judge.json`、`ds_judge.json` がない。確認できた `qa_judge.json` はプロンプト内のものと一致し、保存系列・判定式の数値は照合できたが、元の HDF5・残差・ビルド成果物の独立検証はできていない。E・E′・F の `check_quasisteady` の VERDICT も提示資料にはない。

**ファイル変更・forge 起動は行っていない。提案は plan 未反映であり、反映先は `plans/active/architecture-float-state-double-geometry.md` §6.35 の解釈補足と次試験の事前登録節。**
