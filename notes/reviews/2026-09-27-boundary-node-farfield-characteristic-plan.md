# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `f255f889` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M6/m1
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
目的は妥当ですが、現計画のままでは node 境界から組成・乱流量を正しく流入させられません。  
境界状態だけでなく、実際の流束・評価順序・検証基準を設計し直してから実装すべきです。

対象 checkout `../forge-sern-design` のコードを確認しました。`plans/README.md` と `accepted/` に同機能の完成済み計画はなく、既存の `outlet_statPress` とも役割は異なります。cell を起動時に拒否する範囲設定は、現行の検証方針と整合します。

参照 run の残差履歴・VERDICT 原本はローカルに存在せず、AWS 接続も失敗しました。`check_convergence.py` の実行結果は３件とも `NO residual_history.csv`、保存場への `check_quasisteady.py` は `NO mesh` で、**収束・準定常性は再判定不能**です。以下では、コードから確認できた欠陥、保存場の瞬時値、計画に記録された結果を区別します。ファイルは変更していません。

1. **Critical — §4.3 のゴースト更新では、node の化学種・SST 流入条件が効かない。さらに流向判定が１評価遅れる。**

   **根拠:** [scalarTransport_d.cu:166](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:166) は node 境界で外側のスカラー値を内部節点値に置き換えます。S3 経路も [passiveKernels_d.cuh:25](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:25) で `Yface`・ghost を使わず、内部の `roY/ro` を使います。従って、計画どおり ghost のみ更新しても `Y∞, k∞, omega∞` は輸送されません。

   また、[main.cpp:1458](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1458) ではスカラー境界処理が対流流束計算より先です。その位置で `massflux[ip]` を読むと、今回の境界状態から算出した流束ではありません。流向反転時に熱力学用組成と輸送用組成が食い違います。

   **対案:** ピンなしは維持し、`farfield` 専用の**面スカラー流束**を実装してください。同じ残差評価の境界状態から流向・面組成・面乱流量を決定し、一次輸送と S3 の両経路で `FρY = ṁY_upwind` を使います。`scalarTransport_d.cu`、`passiveKernels_d.cuh`、呼出順序を §5・§7 に追加する必要があります。

2. **Major — 「整合した `bvar` を埋めればよい」という §3 の前提は、実際の境界流束を説明していない。**

   **根拠:** [convectiveFlux_boundary_d.inc.cuh:192](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:192) では質量流束は境界状態ですが、流出時の運動量移流・エネルギー移流は **内部速度・内部全エンタルピー**です。境界状態全体の物理流束 `F(U_b)` ではありません。

   計画式の数値代入でも差が出ます。γ=1.4、両側 ρ=1・c=1、内部法線速度 0.3、自由流法線速度 0 とすると、構成後は `U_nb=0.15, c_b=1.03`。単位面積のエネルギー流束は、現カーネルで **0.442553**、`F(U_b)` で **0.463159**、差は **−4.45%**です。これは流束定義の差の実証であり、物理解に対する誤差率ではありません。

   参照する SU2 は構成状態と内部状態を数値流束へ渡します。forge の現境界処理と同じではありません。[SU2 `BC_Far_Field`](https://github.com/su2code/SU2/blob/master/SU2_CFD/src/solvers/CEulerSolver.cpp)

   **対案:** `farfield` の離散流束を §4 に明記してください。今回の推奨は、`farfield` 専用分岐で構成状態の `F(U_b)` を一貫して使い、既存境界は維持することです。状態構成だけを検証せず、保存量５成分の流束まで試験してください。

3. **Major — 超音速分岐の仕様が矛盾し、非物理状態への対処も未定義。**

   **根拠:** [対象 plan:57](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:57) は、本文で `U_nb ≥ c_b` を判定基準とし、括弧内では内部・自由流の判定が一致しなければ亜音速式としています。

   γ=1.4、`c_i=c∞=1, U_ni=2, U_n∞=0.5` なら、混合不変量から **`U_nb=1.25, c_b=1.15`**。本文では全量内部、括弧内では亜音速式となり、出力が一意に決まりません。また、大きな膨張では `R⁺−R⁻≤0` があり得ますが、そのまま `c_b²` を使うと不正な音速を隠します。

   **対案:** 特性速度の評価状態、分岐の優先順位、音速ゼロ近傍、両側の分類不一致を決定表にしてください。非有限・非正の状態は診断付きで拒否し、上記の反例、音速通過、流向反転を単体試験に追加します。

4. **Major — TP 近似の適用根拠を V1 では検証できない。実際の対象境界は一様外気ではない。**

   **根拠:** [対象 plan:49](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:49) は `γ_i≈γ∞` の妥当性を自由流保持で確かめるとしています。しかし一様流では最初から両者が等しく、異組成・異温度に対する近似誤差は検出できません。`c∞` を `γ∞` と `γ*` のどちらから作るかも明文化されていません。

   保存場 [r4d_view/z2p50H/res_20000.h5](/home/sano/work/forge-sern-design/case/46.sern_design/r4d_view/z2p50H/res_20000.h5) の `z=max(z)` にある 51,143 節点を読み取ると、**瞬時値**は `Y0=0〜0.128410`、`T=176.35〜596.97 K`、`P=2353.6〜41388.9 Pa`。出力から計算した `γ=sonic²ρ/P` は **1.37225〜1.40468**でした。これは定常性の主張ではなく、対象状態の範囲確認です。

   **対案:** frozen-γ 近似としての音速・エントロピー・EOS の契約を統一してください。その上で、内部と外気の温度・組成が異なる流入／流出試験を必須にし、圧力擾乱・反射率・種別質量収支で近似の許容範囲を定めます。V1 は自由流保持だけの試験と位置付けます。

5. **Major — `sstEnergyIncludesK` の内部 k 固定は、新しい流入 k と整合しない。**

   **根拠:** [対象 plan:72](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:72) は境界エネルギーに内部 k を使うと明記しています。[convectiveFlux_d.cu:413](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:413) でも `inlet` 接頭辞以外には境界 k を渡しません。一方、実流束は k から `H*` に `(5/3)k`、圧力に `(2/3)ρk` を加えます。

   指摘1を直して k∞ が流入するようになると、k 輸送は k∞、エネルギー・運動量は k_i という不整合が残ります。

   **対案:** `farfield` の流向で選んだ同じ k を、k 輸送と `H*`・`p*` に渡してください。`k_i≠k∞` の流入試験で全エネルギー収支を確認します。接頭辞判定の踏襲は設計根拠になりません。

6. **Major — V1・V2 が、今回追加する機能の主要な故障を検出しない。**

   **根拠:** [対象 plan:107](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:107) の V1 は一様組成なので、指摘1の「内部組成を使い続ける」実装でも組成条件を通過できます。V2 は斜め衝撃波１条件で、スカラー流入、流向反転、法線方向の音速通過を検証しません。参照領域の反射が評価域へ戻らない条件も、領域長・評価線・判定時刻が未指定です。

   **対案:** SERN の前に、以下を定量ゲートとして追加してください。

   - 内外で異なる Y・k・ω の流入／流出：面流束と体積積分収支を照合。
   - 亜音速の微小音響パルス：入射・反射振幅を分離し、反射率の上限を事前設定。
   - SLAU／ROE、陽解法／block-DPLUR：同じ境界試験で精度・安定性を確認。
   - V2：共通領域の格子、評価線、反射到達位置または物理時間窓を固定。

   ghostless `A⁺` の維持自体は合理的な初版方針です。[timeIntegration_d.cu:857](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:857) と整合しますが、新 BC での安定性まで保証するものではありません。周期・軸対称は検証するか、初版の非対応範囲として明示してください。

7. **Major — V3 は「遠方境界の正しさ」と「広い slip 領域への近さ」を混同している。**

   **根拠:** [対象 plan:114](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:114) は 2.50 H `farfield` と 6.19 H `slip` の一致を合格条件にしています。しかし親計画の R4d は、幅系列について **`check_convergence: NOT CONVERGED`** と記録し、slip 系列の頭打ちを無限遠との一致とは認定していません。[tooling-nozzle-sern-3d.md:1937](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1937)

   両者の差には BC と領域幅の２要因が入ります。また、広い slip 解に偶然近づいても、`farfield` 自身の領域独立性は示せません。

   **対案:** V3 はまず**同一格子・同一新バイナリで side_far だけを替える対照試験**にしてください。次に、合否にかかわらず `farfield` の幅系列を取り、同じ BC で領域感度を判定します。`top_out` などの実効 BC は生成 YAML で固定・照合します。

   `check_convergence.py` と `check_quasisteady.py` の判定区間・閾値・原本保存を §6 に明記し、未収束なら用途を「感度診断」に限定してください。20,000 step 完走や D≤ε だけを機能完成条件にしてはいけません。

8. **Minor — V0 の「更新後状態の md5 一致」は、既存の非決定性と衝突する。**

   **根拠:** [convectiveFlux_d.cu:424](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:424) は、残差の `atomicAdd` により状態は１ step でもビット再現しない一方、初回の面 `massflux` はビット比較できる、と明記しています。

   **対案:** 同一初期状態からの初回面流束・境界状態はビット比較し、更新後保存量は旧バイナリの反復実行で測った再現性幅を基準に比較してください。「既存コード経路を変えない」と「実行結果が常に md5 一致する」を分けます。

**推奨は１つです。** この plan を、**同一評価時点の面状態と保存量・化学種・SST の流束を一貫して構成する node 専用 `farfield`** に改訂し、再レビュー後に実装してください。優先順は **①スカラー経路と流束契約、②特性分岐と熱力学・k 整合、③小規模検証、④SERN 同一格子 A/B と farfield 幅系列**です。既存境界への修正や厳密境界 Jacobian の導入まで広げる必要はありません。

本レビューの提案は **plan 未反映**です。read-only の依頼に従い、記録・修正は行っていません。

指摘数: Critical 1 / Major 6 / Minor 1
