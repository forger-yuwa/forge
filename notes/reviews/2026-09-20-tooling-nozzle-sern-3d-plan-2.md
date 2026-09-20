# codex レビュー: tooling-nozzle-sern-3d (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-3d.md`](../../plans/active/tooling-nozzle-sern-3d.md)
- **stage**: `plan`
- **date**: 2026-09-20
- **commit**: `88eca004` (feature/sern-design)
- **codex**: effort `high`, 6.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: §4.35 R5h: カウル後縁が厚さ0で閉じている件の診断と案(a)有限厚+cowl_base vs 案(b)後縁だけ壁層を粗く、の選択。§4.34 R5g: プラトーをリミットサイクルと断じ受理条件を残差から力の平均±振幅へ書き直す判断。§4.32 側壁のy+1はARと両立しないという結論。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
有限厚カウル案 (a) は設計候補として妥当ですが、現状の診断と受理条件変更を根拠に実装へ進めません。  
生産条件との不一致、温度床の誤判定、リミットサイクルの未立証を先に直す必要があります。

指定プロンプト全文、関連計画・仕様・実装を確認しました。入力生成を独立実行し、既存メッシュ試験は `ALL PASS`。ただし、重点対象の `run_0415`〜`run_0417` の実データはこの checkout に無く、以下では計画掲載値と独立確認できた事実を区別しています。ファイルは変更していません。

1. **Major — R5 系の参照入力は、生産仕様の `frozen_tp`・等温壁になっていません。**

   **根拠:** `problem_3d_sst_cycle_m6on.yaml`、`..._wallres.yaml`、`..._wallres_y.yaml` を読み、現行生成関数で確認すると、すべて **`thermalMethod: 0`、断熱 `wall`、`outlet_statPress`** になります。実保存済みの [run_0400/solverConfig.yaml:4](/home/sano/work/forge/case/46.sern_design/run_0400_3d_base/solverConfig.yaml:4) と [bcondConfig.yaml:3](/home/sano/work/forge/case/46.sern_design/run_0400_3d_base/bcondConfig.yaml:3) も同じです。生成された外部動圧は **60,698.9 Pa** で、入力コメントの71,850 Paとも違います。

   原因は、ガスモデル省略時がCPG、`spec.wall_thermal` 省略時が断熱という既定値です（[probdef.py:62](/home/sano/work/forge/design/forge_design/probdef.py:62)、[同:149](/home/sano/work/forge/design/forge_design/probdef.py:149)）。「等温壁に統一済み」「生産TP条件を再取得済み」とは扱えません。

   **対案:** AWS側の各 run の実効config・BC・組成・バイナリ識別子を回収して条件表を確定してください。CPG・断熱の結果は診断用に限定し、修正した生産入力からTP・等温壁の基準を取り直す必要があります。

2. **Major — `FLOOR_STUCK` は実効温度床を誤認しており、37.6 K／49.8 Kは床到達の証拠になりません。**

   **根拠:** [sern_gates.py:151](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:151) は `thermalMethod` を見ず、全runの温度下限を最低50 Kにしています。一方、単相CPGの実装は `max(intE/(cp/gamma), tMin)` で、`tMin` の既定値は **`1e-4 K`** です（[dependentVariables_d.cu:288](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:288)、[solverConfig.hpp:483](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:483)）。

   したがって、前項の入力なら「50 K以下なので床に張り付いた」という診断は誤りです。低温・高密度の局所異常を調べる必要は残りますが、`FLOOR_STUCK` という分類から原因を逆算できません。

   **対案:** ゲートを実際のEOS経路に合わせ、**床の作動検知**と**物理的な温度の妥当性判定**を分離してください。対象点の保存量・温度の時系列と、EOS前後の保存量変更を確認してからR5hの診断を確定すべきです。ゲート修正だけで既存解を受理することも不適切です。

3. **Major — §4.34の測定では、リミットサイクルも物理的非定常性も立証できません。**

   **根拠:** [plan:1329](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:1329) の指標は、4000 step離れた場の差のL2ノルムです。これは符号と位相を失っており、**単調ドリフトでも一定値になる**ため、周期軌道への到達や振幅を示しません。掲載値でも `roK` の変化は減衰中です。

   また、5反復の終点間差は「4000 step後の再現性幅」であり、残差の丸め誤差床そのものではありません。「変化が反復間差の3〜62倍」から「残り約1桁が到達限界」とは導けません。定常局所擬似時間には物理時間の意味もありません（[solver-settings.md:37](/home/sano/work/forge/procedures/solver-settings.md:37)）。

   **対案:** R5gを未完了に戻し、符号付き局所量・力の密な時系列、複数の独立した末尾窓、CFL変更による平均・変動幅の感度で、ドリフトと反復振動を分離してください。物理的変動と主張するなら物理時間積分による検証が必要です。**現段階で残差条件を置き換えることには反対**します。限定的な目的量評価を認める場合も、`NOT CONVERGED` を保持し、反復誤差・収支誤差の数値上限を別途定めてください。

4. **Major — 案 (a) は、有限厚後縁から下流への接続設計が不足しています。**

   **根拠:** [mesh_sern3d.py:194](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:194) は上下の別節点を `i < i_te` にしか作らず、後縁では同じIDへ合流します。独立したメモリ上の生成でも、**x=0.12 mの25節点が `cowl_in`／`cowl_out` 共通ID**でした。これは後縁より上流にある「別IDの厚さゼロ板」と異なり、両者を同じ欠陥と断定できません。

   この構造で `tk` の終値だけを正にしてタグを追加しても、後縁上下の節点分離と、その間を下流で埋める流体領域が成立しません。さらに側端では `sz=0` に閉じるため、単純なベースquadは退化します（[同:200](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:200)）。

   **対案:** 実装前に、物理入力としての後縁厚さ、上下輪郭、後流ブロック、側端の閉じ方、共通節点IDを図と式で固定してください。`cowl_base` はノズル力に含め、BC・壁距離・初期場・圧力／摩擦／モーメント帳簿を定義する必要があります。検証にはベース面積、面所有数、符号付きJacobian、双対閉性、float32衝突、**生成後のカウル直後間隔**を含めてください。現在の `first_wake_frac` は機体ベース側の処理です。

5. **Major — カウル側端の形状がz方向格子に依存し、固定形状の格子感度試験を妨げます。**

   **根拠:** [mesh_sern3d.py:205](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:205) は、板厚を落とす範囲を「最後の2セル」で決めています。既存の分布関数で独立計算すると、`H=0.1 m` に対し、この物理的なテーパ幅は次のように変わります。

   | `nz_in` / `first_z_frac` | 側端テーパ幅 |
   |---|---:|
   | 25 / 0.004 | 0.868 mm |
   | 57 / 0.000044 | 0.00949 mm |

   約91倍の変化です。側壁を解像する変更がカウル形状も変えるため、その差を格子誤差だけに帰属できません。有限厚後縁を追加すると、この問題はベース面積にも及びます。

   **対案:** 側端テーパ幅・輪郭を格子非依存の物理入力へ移し、粗細格子で輪郭と面積が一致する試験を追加してください。後縁厚さの形状感度と、固定形状の格子感度は分けて評価すべきです。

6. **Major — §4.32は「試した分布では失敗」を「現行トポロジでは不可能」へ一般化し、未解像壁を受理上の残件から外しています。**

   **根拠:** [plan:1137](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:1137) の比較は数種類の点数・分布に限られます。ARは最長辺／最短辺なので、遠方の長辺をさらに分割する可能性と、そのメモリ費用を評価せずに数学的な非両立とは言えません。

   また掲載値では、`sidewall_in` の平均y₁⁺は50.9、`vehicle_side` は16.1、`cowl_in` は1.608、`vehicle_top` は2.354で、壁解像の **`VERDICT: FAIL`** が残っています（[plan:1243](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:1243)）。「場の変化が側壁に集中しない」は、側壁の離散化誤差が小さい証拠ではありません。粗細の推力差0.0033が許容0.002を超える結果も、細側の合格を保証しません。

   **対案:** 結論を「現行分布・試験した規模では不成立」に限定してください。全壁の局所y₁⁺と超過面積、固定形状での追加細分、係数ごとの格子・領域感度を受理の必須条件へ戻すべきです。局所クラスタ化の優先度は残差の発生位置だけで決めず、側壁の精度要件と計算費用で判断してください。

7. **Minor — §5.1が現状と同期しておらず、実装順序の正本として使えません。**

   **根拠:** [plan:1454](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:1454) 以降には、旧診断に基づく優先作業、「根治」という撤回済みの表現、残差PASS必須とR5gによる条件置換が併存しています。§4.23の「残作業表へ」とされた診断項目も、対応する行が明確ではありません。

   **対案:** 履歴は本文に残し、§5.1は現在の未完了事項・依存関係・完了条件に整理してください。共通ゲートの変更は既存の `tooling-convergence-and-wall-resolution-gates.md` と責務を揃え、`methods/design/overview.md` にも反映する必要があります。

**推奨は案 (a) です。ただし、生産条件と診断を訂正し、格子非依存の有限厚カウル＋適合する後流・側端接続を設計してから実装してください。** 優先順は、①実効入力と床ゲートの訂正、②有限厚形状・接続・帳簿の確定、③同一条件での原因検証、④全壁解像と格子・領域独立性の確認です。案 (b) の局所粗化を生産対策にすること、R5gを根拠に受理条件を緩めることは推奨しません。本応答はレビュー提案であり、plan未反映です。

指摘数: Critical 0 / Major 6 / Minor 1
