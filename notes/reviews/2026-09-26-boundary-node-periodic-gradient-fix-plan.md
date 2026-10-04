# codex レビュー: boundary-node-periodic-gradient-fix (plan)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `97d9c27a` (feature/sern-design)
- **codex**: effort `high`, 4.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **focus**: 実装前の方針レビュー。(i) §4.1 の合併 stencil LSQ (重複配分 α、M_r、打ち切り 1 回、係数焼き込み、毎 step の和) が全周 (周期なし) の一意 stencil LSQ と一致するか、実装箇所 calcGradient_d.cu:577-920 と mesh.cpp の group 情報で実現できるか、2/4/8 member の角・壁∩継ぎ目で破綻しないか。(ii) §4.2 の k/ω 勾配 gather の位置 (ransGradient 直後) と周期半割面の除外で正しくなるか、他に同じ型 (gather 後の作り直し) が無いか (§2 の点検表)。(iii) §6 の試験で修正の正しさと既存 run への影響を判定できるか。
- **extra**: `plans/active/boundary-node-rotational-periodic.md`, `notes/reviews/2026-09-26-boundary-node-rotational-periodic-plan-2.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
合併 stencil の LSQ と、`ransGradient` 直後の専用 gather を採用する方針は支持します。  
ただし、重複同定の精度仕様、軸対称の除外、作用素試験と回帰試験の合格条件を実装前に補う必要があります。

確認対象は `97d9c27a` と現在の作業ツリーです。`mesh.cpp/.hpp` には回転周期に関する未コミット変更があります。ファイルは変更していません。

目的の同定は正しいです。`plans/README.md` と `accepted/` を確認し、既存の cell 周期保存修正・node 境界疑似点撤去とは別の未解決問題と判断しました。

§4.1 の式も、**同一隣接の幾何と状態差が一致する**条件では正しいです。隣接の同値類を \(E\) とすると、

\[
\sum_{m,j}c_{mj}(\phi_j-\phi_m)
=M_{r,\tau}^{+}\sum_E w_Ed_E(\phi_E-\phi_r)
\]

となります。読み取り専用の CPU 数値実験でも、重複を持つ非対称 stencil を 2・4・8 member に配分した結果は、一意 stencil の解との差が最大 \(3.4\times10^{-16}\) でした。これは式の確認であり、forge の実装試験ではありません。

1. **Major — 重複同定に使う座標精度と一致判定が未定義**

   **根拠:** [plan:53](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:53) は `periodicRoot` と \(\Delta x\) の一致を要求します。一方、現行 LSQ は `flow_float` の座標を double に昇格して差を取っています。[`calcGradient_d.cu:653`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:653)

   float32 座標から現行と同じ順序で計算すると、同じ並進像の辺でも、

   ```text
   float32(0.2) − float32(0.1) = 0.10000000149011612
   float32(1.2) − float32(1.1) = 0.10000002384185791
   ```

   となり、完全一致しません。重複を取り逃がすと、二次場では丸め誤差より大きい差になります。今回の非対称 stencil の数値実験では、一意化ありの勾配 `(0.512338, 0.353896)` に対し、重複込みでは `(0.492580, 0.192580)` でした。

   **対案:** §4.1 に次を固定してください。

   - 周期像の識別と幾何の照合を分け、照合許容を局所辺長・座標精度から定義する。
   - 同定した隣接ごとに共通の \(d_E,w_E\) を決め、**行列組立と係数生成の両方で同じ値を使う**。
   - 配分は原則 \(\alpha=1/\text{重複数}\)。異なる周期像を統合しない。
   - 座標原点の移動、斜め並進、root 交換を G2 に追加する。

   `periodicRoot` と内部面 incidence から実現可能です。ただし、現在の group 情報を渡すだけでは重複配分まで完成しません。

2. **Major — 軸対称では「係数を合併して毎 step に和」の前提が成立しない**

   **根拠:** [`periodicNode_d.cu:174`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:174) は `isAxisymmetric == 1` で gather を終了します。LSQ の事前計算側には同じ除外がありません。提案どおり係数だけ合併すると、軸対称では**合併勾配の部分寄与だけが残る**危険があります。

   また、GG は軸対称で `A_planar` を使いますが、[`mesh.cpp:773`](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/mesh.cpp:773) が合併するのは `volume` です。既存 gather の軸対称除外を単に外す修正も不適切です。

   **対案:** 今回は**非軸対称の node 並進周期に限定**してください。LSQ の合併係数生成と SST gather に同じ適用条件を設け、軸対称の既存経路を変えないことを試験します。「すべての node 並進周期」「他スカラーは正しい」という記述にも、この限定を付けてください。

3. **Major — G1 の真値比較だけでは GG の合併を検証できず、G0 の角試験にも入力契約が必要**

   **根拠:** G1 は線形場の解析勾配に相対誤差 \(10^{-5}\) を要求します。しかし SST は内部面で節点値を補間し、壁境界では owner 値を使う GG です。[`ransTransport_d.cu:39`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransTransport_d.cu:39)  
   一般メッシュで線形厳密ではないことは、[現在仕様:78](/home/sano/work/forge-sern-design/methods/gradient.md:78) にも明記されています。

   具体的に、壁上の原点を共有する二つの三角形
   \[
   [(0,0),(2,0),(0,1)],\quad[(0,0),(0,1),(-1,0)]
   \]
   の median-dual を合併し、現行の面値規則で \(\phi=x\) を計算すると、GG 勾配は **\((1,0.5)\)**、真値は \((1,0)\) です。gather が正しくても G1 は落ちます。

   また、G0 の「局所作用素試験」という区別は適切ですが、§5.1 の線形場 HDF5 を通常の初期化へ渡すだけでは不十分です。[`main.cpp:1257`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1257) は周期状態を root からミラーします。三重周期の角では、非定数の大域線形場は周期条件を満たしません。

   **対案:** 試験を次の契約にしてください。

   - **G0/G2:** 各対象 group の展開した局所座標で場を作り、BC・時間更新による変更前の作用素を比較する。
   - **G1:** 一様直交格子では解析解、非対称・壁∩継ぎ目・2/4/8 member では **CPU double の合併 GG** を参照にする。
   - 非零定数場でゼロ勾配、ゼロ成分の絶対誤差も確認する。
   - 部分 stencil は退化するが、合併後には rank が回復する例を明示的に含める。

4. **Major — G2・R1・R2 に、実装の合否を確定できない条件が残っている**

   **根拠:** [plan:99–101](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:99) の「float32 丸め内」「同水準」「差を記録」は、閾値・尺度・比較区間を確定していません。R1 は VERDICT の併記を要求していますが、どの判定を受け入れるかも未定義です。

   `case/39.periodic_hills/run_0007_coarse_rans/` に判定ツールを実行した結果は、

   ```text
   NO residual_history.csv
   OVERALL: CHECK FAILURES ABOVE
   ```

   でした。この checkout の基準 run を収束済み対照として利用できません。G0 も現在確認できるのは **seam 平均 2.000000／内部平均 1.000000 の記録**までで、入力・出力 HDF5 はありません。

   **対案:** §5・§6 に以下を実装前の条件として追記してください。

   - G2 の比較対象を「格納済み float32 入力を double で評価した参照」と明確化し、係数丸め・積和・gather の誤差尺度を固定する。
   - R1 は同一の整備済み設定・メッシュ・IC で旧新を再生成する。`wallTreatmentSST: 0`、実際の `kInit/omegaInit`、段階起動、メッシュ品質・壁解像の確認を明記する。
   - R1 の定常評価は、同一設定区間の `check_convergence.py` **PASS** と対象量の準定常判定を受入条件にする。継ぎ目指標のノルム・比較列・許容倍率を固定する。
   - R2 はスキーム、粘性、共通 dt、終了時刻、保存誤差の正規化を指定する。非粘性 KEEP の基準は[既存手順:65](/home/sano/work/forge-sern-design/procedures/verification/09-taylor-green.md:65)から採用できるが、勾配修正の影響を見る粘性・再構成経路も明示する。
   - 固定状態で、非周期 node の LSQ と化学種・受動種の GG が変わらない回帰を追加する。

   `case/09` の選択と、TGV に定常 PASS を要求しない判断は妥当です。陰解法の行縮約を別計画に残す判断も支持します。

5. **Minor — 勾配以外にも、初回 F1 が生成後に上書きされる経路がある**

   **根拠:** 通常の残差組立では `ransBlendF1` の後に `ransTransport` を呼びます。しかし [`ransTransport_d.cu:104`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransTransport_d.cu:104) の `buildScalarDescs` は、初回に `sstF1` を全点 **1** で埋めます。したがって `sstSigmaBlend` 有効時の初回輸送は、直前に計算した F1 を使用しません。

   調べた四つの勾配生成経路では、SST 以外に同じ「gather 後の勾配再生成」は見つかりませんでした。ただし、F1 まで含めた生成・使用順序にはこの例外があります。

   **対案:** F1 の初期充填を変数初期化時へ移し、`ransBlendF1` 後には上書きしない構成にしてください。G1 では初回と二回目の残差組立を確認し、輸送が実際に読む F1 も検査します。

**推奨は、非軸対称の並進周期に限定して、共通幾何を使う一意 stencil LSQ と SST 専用 gather を実装することです。** 実装前の優先順は、①指摘1・2の適用契約、②指摘3の参照試験、③指摘4の受入条件、④指摘5の初期化整理です。`ransGradient → 専用 gather → ransBlendF1 → ransTransport` の順序は正しく、周期半割面の除外は node 周期だけに限定してください。

読み取り専用レビューのため、**plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
