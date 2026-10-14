# codex 諮問 (diagnose): core-grid-dtedge

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-dtedge.md`](../../notes/reviews/briefs/2026-10-11-core-grid-dtedge.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `27d0dee4` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.3 min, rc=0
- **結論**: **M-d の上限 50／1 を同じ条件で各 1 回だけ再実行し、同腕内の RHS 不一致と u5 の再現性を確認する。既存の「判定不能」は書き換えない。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **M-d の上限 50／1 を同じ条件で各 1 回だけ再実行し、同腕内の RHS 不一致と u5 の再現性を確認する。既存の「判定不能」は書き換えない。**

採否表:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 資格 (1) を直ちに「相対差 ≤ 1e−8」と読み替える | **却下**。判定器はハッシュ一致を要求し、実出力も `same_input: false`。7.0e−9 という小ささだけでは、加算順の非決定性も補正への無影響も確認できない。再実行による確認と、新しい資格の登録を先に行う。[判定器:16](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_sweep_judge.py:16)、[結果:6](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/_band_ab/core_grid/sweep_judge.json:6) |
| Major | 「集中域自身の dt 変化では説明できない」 | **撤回を採用**。根拠だった 0.977 は全域の中央値。ブリーフの残差最大帯では dt が約 3.3 倍違う。局所 dt・時間項の作用を候補に戻す。ただし、それだけで原因確定とはしない。[plan §4.36:746](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:746) |
| Major | H-edge と「上限 50→5 で境目だけを動かす」試験 | **要再検証／試験案は却下**。上限は dt の大きさと時間対角も変える。さらに拡大率 ≤ 5 の節点は初期状態で変わらず、拡大率が 1 になる境目も動かない。境目の位置だけの介入になっていない。[setDT_d.cu:194](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:194) |
| Major | M-d から「内部反復に問題なし」「H-LHS が真因」と結論する | **却下**。小さい u5 が否定的なのは、測定した状態・領域・右辺での「大きな 5→40 sweep の解き残し」。十分解いた近似左辺が外反復を減衰させるかは別問題。[plan §4.42:807](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:807) |
| Minor | 「両腕とも 5 回と 40 回が 0.1 % 以内で一致」 | **訂正を採用**。上限 50 の最大 u5 は 0.0917 %、上限 1 は 0.1643 %。後者は 0.1 % を超える。ただし登録閾値 u5 ≤ 0.1、すなわち 10 % は十分下回る。[集計結果:10](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/_band_ab/core_grid/sweep_judge.json:10) |

第 1 仮説: **十分解いた近似左辺と、実際の残差の応答との不整合が、境界層外縁のモードを減衰しにくくしている。** 確度: 中。ただし dx 依存の機構は未確定。

  根拠: `case/45.isobutane_m6_d155/run_0628_cg_sweep_cap50/` の島の u5 は、ρ が 3.56e−4、ρE が 9.06e−4。大きな内部解き残しを主因とする見立てには否定的である。一方、時間項 V/Δτ は対角へ入り、熱伝導の近似も自節点の対角へ加わる。[時間項:863](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:863)、[熱伝導項:1026](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1026)

  また、各 sweep 内の緩和により、解いている系は plan の表記で `(L/ω + N)δ = b₀`。ω = 0.7 なら、無緩和の `L + N` に対して約 `0.429 L` が加わる。**線形系がよく解けることと、外反復が速く進むことは同義ではない。** この追加項の当該モードへの作用は未測定であり、原因認定はできない。[緩和処理:2146](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2146)

  反証条件: 同一状態で十分解いた補正を使った残差の応答を調べ、上限 50 の当該モードが上限 1 と同等以上に減衰し、停滞を説明する悪い応答が認められなければ、この局所的な説明は支持されない。

第 2 仮説: **dt の空間的不均一性がモードの減衰を悪化させる H-edge。** 確度: 低。温度・μt・格子幅・残差の分布と交絡しており、現状は位置の相関にとどまる。

第 3 仮説: **μt・熱伝導係数の平均流への依存が関与する。** 確度: 低。F1 が止めたのは ρk・ρω の更新で、μt は 1.8〜4.2 % 変化した。F1 の結果を覆す必要はないが、「μt の経路まで除外済み」とは扱えない。[plan §4.39:780](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:780)

判別 A/B: **新規の機構試験より先に、既存 M-d の再現性を確認する。**

- 腕 A は上限 50、腕 B は上限 1。変更点は `lineDtDirectionalCap` だけ。既存と同じ起点・バイナリ・環境で、各 1 外反復・41 sweep、約 10 秒ずつ。新しい run に出力する。
- 各腕を既存の `run_0628_cg_sweep_cap50/`・`run_0629_cg_sweep_cap1/` と比較する。状態・節点の一致、成分別 RHS 差、後退誤差、g₄₀、u5 を見る。RHS 差は全ラインと島を分け、同じ体積重み付きノルムで測る。
- **判定 A**: 同腕でも微小な RHS 不一致が再現され、同腕・腕間とも登録した許容差内、資格 (2)(3) と u5 ≤ 0.1 が再現する → 「設定差がなければ RHS 不一致は生じない」を退ける。非決定的な数値差という説明を支持し、**新資格による再試験の結果 B**として記録する。
- **判定 B**: 許容差を超える、または u5・資格 (2)(3) が再現しない → 「微小な入力差だから判定に影響しない」という扱いを退け、入力・書き出し経路の診断を優先する。同腕がビット一致しただけの場合は非決定性を確認できておらず、原因を二択で断定しない。

記録は「旧登録: 判定不能」「再試験前に登録した数値許容差と変更理由」「再試験の結果」を分ける。**1e−8 は今回の差を見た後の基準であることを明記する。** 結果 B の意味も「この右辺に対する大きな内部解き残しを退ける」に限定する。

やらない方がよいこと: 上限 5 の 2 万 step、μt 固定の実装、Rhie–Chow 型補正の追加を先行させること。Grefc の発達した別の場との分布比較も、格子粗化の因果の判別にはならない。

呼び出し側の前提への異議:

- **方向別 dt は粘性 CFL だけを外す処理ではない。** 面の CFL は対流・音響・粘性の和で、ライン面ではその全体を除外する。「μt が大きいから拡大率が大きい」を確定するには、実際に最大値を取った面と各項の分解が必要。[setDT_d.cu:82](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:82)、[同:167](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:167)
- 0.977 に依存した「島の dt は不変」「壁側からの間接作用が必要」という推論は撤回する。D1 の交互更新・増分消失に対する否定的結果、F1・M-a の限定された結論、Grefc の残差改善の観測は、その訂正だけでは取り消されない。
- 上限が同じ状態の RHS を変えないことは、上限操作の作用が反復側にある根拠になる。しかし、**空間離散化の Jacobian が反復の難しさに関与しないという根拠にはならない。**

不足情報: 同腕の再現性、CFL の律速面と項別寄与、十分解いた補正に対する実残差の応答。巨大ファイルの禁止に従い、今回の数値確認は許可された集計 JSON と提示値に限った。ファイル変更・forge 起動は行っていない。**plan 未反映**であり、本提案の登録・反映は依頼どおり呼び出し側が行う。
