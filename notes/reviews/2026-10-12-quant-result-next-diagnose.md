# codex 諮問 (diagnose): quant-result-next

- **brief**: [`notes/reviews/briefs/2026-10-12-quant-result-next.md`](../../notes/reviews/briefs/2026-10-12-quant-result-next.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `62f02796` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **P1 を修正して採用し、G0c の継続場から `space.pRef: 0 / 5499000` を両精度・各2本・3,000 step で比較する。FP64＋pRef の対照も加える。**
- **extra**: `case/45.isobutane_m6_d155/rqs_an.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rqs_judge.json`, `case/45.isobutane_m6_d155/rpref.sh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **P1 を修正して採用し、G0c の継続場から `space.pRef: 0 / 5499000` を両精度・各2本・3,000 step で比較する。FP64＋pRef の対照も加える。**

第 1 仮説: 圧力ゲージ変更で抑えられる演算誤差が今回の差の発達にも寄与し、pRef だけで登録窓の温度差 RMS を少なくとも半減できる。確度: **中**。

  根拠: G0 の既存試験では、`case/45.isobutane_m6_d155/run_0703_rp_fp1`〜`run_0705_rp_fp3` の温度差 RMS は1.5〜1.6 Kで、無介入の3.7〜4.5 Kより約6割小さい。ただし正式判定は判別不能で、G0c・継続場への外挿は未確認（[plan §6.45](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1049)）。

  実装では pRef を SLAU の面圧力から差し引き、運動量流束に使い、hoop にも適用する。したがって試せるのは**この圧力ゲージ変更の十分性**である（[convectiveFlux_slau_d.inc.cuh:618](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:618)、[axisymmetricSource_d.cu:52](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:52)）。

  反証条件: 下記の対照条件を満たしても、pRef 後の差が両窓・両反復で無介入の半分を超えて残ること。

第 2 仮説: pRef が変更しない原始量・物性・残差・陰解法内部の精度差と、B0 の反復との相互作用が大きく、pRef 単独では半減できない。確度: 中、寄与箇所は未確認。

第 3 仮説: 精度依存の初期化差が主要な種になっている。確度: 低。共通入力でも C1/F1 の `res_0` の ρE 最大相対差は4.898e−5であり、除外されていない（[rqs_judge.json:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rqs_judge.json:75)）。

判別 A/B: **変更は `space.pRef` の1点だけ。** G0c・B0・上限50・`_ic/cont_q32.h5`・同一ソース版を固定し、C＝FP64、Cp＝FP64＋pRef、F＝float、Fp＝float＋pRef を各2本、3,000 step・250 step ごとに保存する。

- 評価域と W1/W2 は §6.54 を維持する。E₀ は F と C 平均の温度差 RMS、Eₚ は Fp と Cp 平均の温度差 RMS。N は C 同士・Cp 同士の再実行差 RMS の大きい方。G は Cp と C の全組合せの差 RMS の最大値。
- 前提は、両窓で `min E₀ ≥ 10N`、`G ≤ 0.05 × min E₀`。さらに G0c の実際の面ベクトルと hoop 面積で、§6.44 の閉性条件 `g ≤ 1e−12` を確認する。旧 G0 の検査結果を流用しない。
- 入力・バイナリ・設定差を照合し、同じ精度の pRef 有無で `res_0` の9保存量と座標が一致することを確認する。必要列の欠損・全残差と保存場の非有限値は前提不成立とする。
- **A：** 両窓で `max Eₚ + N + G ≤ 0.5 × min E₀`。第1仮説の有限期間の予測を支持し、「pRef 単独では半減できない」を退ける。
- **B：** 両窓で `min Eₚ − N − G > 0.5 × max E₀`。第1仮説を退ける。ただし圧力関連の誤差全般を除外する根拠にはしない。
- 中間・反復間の不一致・前提不成立は判別不能。以上は**新試験用の提案閾値**で、過去の判定は変更しない。全残差の `check_convergence` と対象量の準定常判定も添え、有限期間の効果と収束を分ける。

やらない方がよいこと:

- **P2・P4を先に回して、3,000 step 離れなければ生産に使えると判断すること。** 上限変更は平均流・SST の擬似 dt を縮めるため、短期の小差だけでは抑制と発達の遅延を分けられない（[solver-settings.md:76](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:76)）。生産設定の探索は、その後に行う。
- **P3 の部分 double 化を今すぐ実装すること。** §6.54 は右辺への局在化をしておらず、現時点では対象を選ぶ根拠が弱い。
- pRef 試験と `lowMachPrecond`・Thornber・左辺修正を同時に行うこと。pRef の成否から「低マッハ補正不足が真因」とは判定できない。

呼び出し側の前提への異議:

- **採用：§6.54 の結果 B と、その限定された解釈。** 保存済み JSON から再計算した `(max E_Q + N) / min E_F` は W1＝**0.02431**、W2＝**0.05762**で、両窓とも0.1以下だった。根拠 run は `case/45.isobutane_m6_d155/` 配下の C：`run_0737_qs_c1`・`run_0740_qs_c2`、Q：`run_0738_qs_q1`・`run_0741_qs_q2`、F：`run_0739_qs_f1`・`run_0742_qs_f2`。記録された収束判定は6本とも **`OVERALL: CHECK FAILURES ABOVE`**。したがって「この起点・評価域・3,000 step では末尾量子化だけで差を再現できない」までを採用する（[判定記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rqs_judge.json:82)、[run 索引](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:222)）。

- **Major／却下：「残りの大部分は残差演算の誤差」という寄与率の解釈。** Q の差も N の15〜20倍で成長している。量子化試験の差と float の差は、供給される摂動も反復経路も異なるため、RMS の比を原因の内訳にできない。対案は、右辺・途中の格納・陰解法・初期化を候補として残すこと（[plan §6.54 結果](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1356)）。

- **Major／却下：「左辺は両方単精度だから候補から外れる」。** `ST=float` でも、取り込む `gamma`・`dt_local`・`ro`・速度・音速・`Ht`・乱流粘性は各ビルドの計算結果であり、同一値とは限らない。該当箇所は試験 commit `b0303547` でも確認した。対案は、**組立精度が同じことと、組み立てた左辺が同じことを区別する**こと（[timeIntegration_d.cu:829](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:829)）。

- **Major／要修正：P1 に Cp がない。** pRef は内部面・境界面・hoop を変更する。G0c で FP64 の基準自体が動かないという確認が必要である。対案は上記の4腕と閉性検査。既存の §6.44 も Cp の関門を設けていた（[plan:1034](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1034)）。

- **Minor／要補強：「前提はすべて満たした」の検査範囲。** `rqs_an.py` は HDF5 の有限性を選択時点の `T` だけで検査し、残差の必須列の存在は検査していない。これだけでは全保存量・全保存時点の有限性を保証しない。対案は次の判定器で必要キーと全出力を検査し、既存結果にも確認範囲を明記すること（[rqs_an.py:42](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rqs_an.py:42)）。

不足情報: 対象6本の run ディレクトリはこの作業ツリーに存在しない。直接確認できたのは判定コード・保存済み JSON・ソースと計画の記録までで、CSV/HDF5 の再集計、実バイナリの照合、準定常判定の再実行はできていない。

**plan 未反映。** 読み取り専用の依頼に従い、ファイル変更・forge 起動は行っていない。
