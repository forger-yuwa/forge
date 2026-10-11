# codex 諮問 (diagnose): float-after-freeze

- **brief**: [`notes/reviews/briefs/2026-10-11-float-after-freeze.md`](../../notes/reviews/briefs/2026-10-11-float-after-freeze.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `8fcef93a` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: ?

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major 1** | §6.29 の「判別不能」は採用。失敗理由の説明は訂正 | 提示表を判定式で再計算すると「残る」1/6・「抑えられる」0/6。ただし **W1 は再実行差の条件だけでなく、B′ が3断面すべてで 0.5｜D｜未満**（0.927<0.946、0.737<0.824、0.619<0.716）。「10倍条件だけで外れた」は誤り。また「抑制を支持できない」から SST 更新の必要性を棄却することもできない。[plan:766](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:766)、[判定式:77](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fz_an.py:77)。**両仮説とも採否未確定**と記す。 |
| **Major 2** | 「判定の前提を満たした」は要再検証 | `fz_an.py` は `check_convergence` だけを呼び、**`check_quasisteady` を呼ばない**。抽出後の θ_r・Q_w の有限性検査もない。[fz_an.py:37](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fz_an.py:37)、[同:66](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fz_an.py:66)。判定式へ NaN を入れると、入力異常と区別せず両条件が偽になることを確認した。**既存出力で抽出系列の有限性・点数・共通座標を確認し、対象量・区間・閾値付きの準定常 VERDICT を追加**する。元の判定閾値は変更しない。 |
| **Major 3** | 次の一手は候補(a)を条件付き採用、(b)は保留 | 既存診断は `fabs(d−act)` を集計し、領域も内部／境界だけなので、収縮部の符号付き偏りを示さない。[commitLossDiag.cu:59](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/commitLossDiag.cu:59)。さらに平均流 commit の基準は `Q_N` であり、単純な前後スナップショット差ではない。[update_d.cu:308](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:308)。**実際の基準・要求更新・commit直後を同時に記録する診断**を先に行う。 |
| **Major 4** | B–B′ の増大を原因特定や閾値緩和に使う案は却下 | 凍結時の B−B′ は、W1 の約 +0.18〜+0.24 ポイントから、W2 の約 −0.06〜−0.11 ポイントへ**全断面で符号反転**する。[plan:766](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:766)。凍結による減衰・位相・感度の変化とは整合するが、乱数的な分散の増加や SST の安定化作用は未証明。**既存判定を維持し、次の固定状態診断では両精度を反復評価**する。 |
| **Minor 1** | 記録文は修正して採用 | 「float が 0.6〜1.2 % 動いた」は Δ と δ の混同。表の **自分の起点からの ΔB・ΔB′ は約 +0.49〜+1.06 %**、**FP64 の移動を差し引いた δ は約 +0.62〜+1.17 ポイント**。[plan:764](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:764)。後者であることを明記する。「凍結なしの約5〜8割」は記述的な比に限り、原因の寄与率とは呼ばない。 |

**結論:** 候補(a)を採用し、共通 Q32・SST 輸送更新停止のまま、float／FP64 の最初の1更新を段階別に計測して、差が commit 前と commit 丸めのどちらに主に現れるかを調べる。

**第 1 仮説:** 共通 Q32 での最初の平均流更新差は、commit 前の演算差が主に説明する。**確度: 中。ただし今回の凍結条件では未確認。**

根拠: V3 には commit 前の残差差が残り、壁際の `res_roUy` は参照ノルム比 1.42、`res_roe` は 0.42 と記録されている。[plan §6.10](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:451)。これは候補を残す根拠であり、今回の起点・修正後の格子での実測ではない。凍結は SST 更新ブロックを飛ばすが、平均流更新と化学種更新は残る。[main.cpp:2366](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2366)、[同:2391](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2391)。

反証条件: 同じ起点の診断で、再実行差を十分上回る実更新差があり、その差の90%以上を commit の丸め差で説明でき、commit 前の寄与が10%以下なら、この**最初の更新についての仮説**を棄却する。

**第 2 仮説:** 平均流の commit 丸めが局所的に同じ向きへ偏る。確度: 低〜中。既存の `S_lostfrac` は符号を持たず、支持にも除外にも足りない。

**第 3 仮説:** 凍結で反復の減衰・位相が変わり、小さな演算順序差が増幅される。確度: 低。B−B′ の窓間の符号反転とは整合するが、発生源は未確認。

**判別 A/B:**

- **変える点:** ビルド精度だけ。A＝FP64、B＝float。共通 `_tr/q32_init.h5`、格子、実効設定、環境変数を固定し、両方で `FORGE_FREEZE_TURB=1`。
- **長さ:** 各1 step、各3回。同じ Q32 から毎回開始し、連続した3 stepにはしない。
- **記録:** 平均流5保存量について、組立後の残差 R、実際の commit 基準 b、緩和・ガード適用後の要求更新 d、commit直後 q。境界条件・EOSなど、その後の上書きは別に記録する。
- **比較:** 共通起点 Q₀ に対し、a＝q−Q₀、e＝q−b−d と置く。精度間差には  
  **Δa＝Δb＋Δd＋Δe**  
  が成立する。これを double で照合し、収縮部の内部・軸近傍・壁近傍を事前固定して、保存量ごとの符号付き集計と L2 を出す。異なる保存量を一つのノルムに混ぜない。

事前の判別は次のとおり。

- **結果A:** ｜｜Δe｜｜ ≤ 0.1｜｜Δa｜｜なら、最初の更新差は主に commit 前で生じる。第1仮説を支持し、**この更新差の主因としての**第2仮説を退ける。
- **結果B:** ｜｜Δb＋Δd｜｜ ≤ 0.1｜｜Δa｜｜なら、第1仮説を棄却し、commit 丸めを優先する。
- 判定対象の差は、両精度内の再実行差の最大の10倍超を要求する。中間・領域間の不一致・差が小さい場合は判別不能。10%・10倍は診断用の暫定基準で、信頼区間ではない。

**この1 stepで分かるのは差の発生段階までであり、30,000 step後の θ_r の偏りの主因までは確定しない。**

**やらない方がよいこと:** 今の結果から平均流だけの double commit を実装すること、`qAccumulatorFP64` の軸対称拒否だけを外すこと、30,000 stepを延長して合格窓を探すこと。軸対称の拒否には、commit基準の射影を FP64 正本へ反映する必要があるという理由がある。[main.cpp:3396](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3396)。

**呼び出し側の前提への異議:** §6.29 は「登録窓で残差・派生量を確認した上での記述」と「因果判定」を分ける必要がある。記載案は、次までに限定する。

> 提示された窓平均では、凍結した FP64 の移動を差し引いた float の θ_r の差は +0.62〜+1.17 ポイントだった。事前登録の判定は判別不能。SST 更新の必要性、平均流 commit の支配性は確定していない。

**不足情報:** 対象の `case/45.isobutane_m6_d155/run_0491_fz_fp64/`、`run_0492_fz_f32/`、`run_0493_fz_f32b/` と比較元の `run_0487〜0489`、`fz_judge.json`・`tr_judge.json` は手元に存在しない。このため本レビューで独立確認できたのは**コードと提示表の算術**まで。全期間の有限性・凍結・ビルド一致は未確認で、収束は plan 記載の **`CHECK FAILURES`** の引用に限る。凍結後の準定常 VERDICT は提示されていない。

ファイル変更・forge 起動は行っていない。**plan 未反映**。反映先は同 plan §6.29 と §5.1 #15。
