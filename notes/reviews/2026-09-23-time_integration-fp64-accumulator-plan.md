# codex レビュー: time_integration-fp64-accumulator (plan)

- **plan**: [`plans/active/time_integration-fp64-accumulator.md`](../../plans/active/time_integration-fp64-accumulator.md)
- **stage**: `plan`
- **date**: 2026-09-23
- **commit**: `75f01350` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

commit の丸め落ちは実測で確認でき、改善対象として妥当です。  
ただし、FP64 の値を次ステップまで保持する設計が未完で、倍精度ビルドの改善を「commit だけで再現できる」とする根拠も不足しています。  
さらに、現行 G2 は参照 run 自身が通りません。計画を改訂してから実装すべきです。

`plans/README.md` と関連 accepted plan を確認した範囲では、この保持精度の改善は解決済みではありません。既存の `implicitSolvePrecision` は線形解法内部の精度切替で、今回とは別です（[implementation.md:219](/home/sano/work/forge/methods/time_integration/implementation.md:219)）。以下はファイルを変更せず、コード確認と既存データの再計算で得た指摘です。

1. **Major — FP64 の増分が次ステップまで残る設計になっていない**

   **根拠:** 現行 commit は `Q += dq` ではなく `Q = Q_N + dq`（[update_d.cu:249](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:249)）。ステップ末尾で `Q` を `flow_float` の `Q_N` にコピーします（同ファイル [39行](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:39)）。`Q` 5 本だけの倍精度化では、このコピーで下位ビットを失います。

   計画の代表値を使った演算確認でも、`ρ=0.017755` に `dq=2.97e-10` を100回加える際、毎回 baseline を float32 に戻すと累積増分は **2.97e-10**、FP64 で保持すると **2.97e-8** でした。

   また、in-place 加算に変更するだけでも不十分です。`dependentVariables_d.cu` は float の `ro_temp` を密度へ書き戻し、エネルギーも再構成します（[76行](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)、[207行](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:207)、[296行](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:296)）。これは読出し専用処理ではありません。

   **対案:** 実装前に、FP64 状態の正本と各更新処理の関係を定義してください。`Q_N`、EOS、床処理、壁ピン、周期コピー、SST エネルギー補正を含む**書込み側**の棚卸しが必要です。微小増分を複数ステップ保持する試験を最初のゲートに置くべきです。

2. **Major — 「原因は commit に限局」「FP64 accumulator なら所定の床に届く」は未実証**

   **根拠:** `run_0020_double` は `flow_float` と `geom_float` をともに変更しています（[case README:387](/home/sano/work/forge/case/56.gap_tp1187/README.md:387)）。したがって、状態保持だけでなく流束、EOS、残差、線形解法、SST 更新なども変わっています。

   閉性誤差ゼロは、非一様流の流束評価・総和・勾配計算まで丸め誤差ゼロであることを意味しません。ブロックサイズ対照と sweep 数比較でも、この交絡は除けません。

   今回 `run_0017_resro_dump/res_1000.h5` の `z/W>10`、16,607 CV を再計算すると、平均 `|dq|/ULP=0.15977`、平均の実効更新／意図した更新は **0.259%** でした。これは丸め落ちの存在を支持しますが、唯一の律速要因である証明ではありません。

   **対案:** §3 を「commit 丸めが実在する。改善の十分性は未検証」に修正し、同一 IC・同一ブロックサイズによる **FP32／accumulator のみ FP64／全体 FP64** の比較を登録してください。§4.1 の `5.30e-15` は収束床の保証から外すべきです。更新を保持できても、FP32 の残差評価には別の限界が残ります。

3. **Major — G2 の参照値が終端挙動を表さず、合格条件も参照 run と矛盾する**

   **根拠:** `case/56.gap_tp1187/run_0020_double/` を再計算した結果です。

   | 確認 | 結果 |
   |---|---|
   | `mdot_decay.py`、既定 `z/W=15`、25k–400k | **べき乗則**。fit 残差 RMS: べき乗 0.8367、指数 1.0427 |
   | 300k の `|mdot|` | `4.5634e-12` |
   | 400k の `|mdot|` | **`1.4798e-11`** |
   | `check_convergence.py`、全400k区間 | **`NOT CONVERGED (stalled/plateau)`** |
   | `check_quasisteady.py`、保存済み `gap_series.csv` | `zW422`・`U_rms_deep`: **STEADY**、`zW844`・`zW1411`: **DRIFTING** |

   全残差判定の末尾値は `rms_ro=2.34e-6`、`rms_roUx=5.01e-3`、`rms_roUy=8.51e-4`、`rms_roe=5.07`。保存済み [CONVERGENCE_VERDICT.txt](/home/sano/work/forge/case/56.gap_tp1187/run_0020_double/CONVERGENCE_VERDICT.txt) と再実行結果は同じです。

   `mdot` 全時系列を `check_quasisteady.py` の `classify_series` に既定閾値で渡した結果も **DRIFTING**、drift **137.5%/tail** でした。300k の値を定常的な床として使えません。なお、この判定は「初期の指数的減衰が存在しない」という意味ではなく、床付近まで一括 fit する現在のゲートが不適切ということです。

   **対案:** 減衰区間の fit と末尾区間の床・変動を分離し、区間・指標・許容値を事前登録してください。G2 に全保存量の収束 VERDICT と対象量の準定常 VERDICT を追加し、局所的な改善と全体収束を区別してください。SU2 比較にも両者の判定が必要です。

4. **Major — runtime 型切替と restart の設計が欠落している**

   **根拠:** 状態管理は `std::map<std::string, flow_float*>`、host 側も `vector<flow_float>` です（[variables.hpp:22](/home/sano/work/forge/solver_density_cuda/variables.hpp:22)）。確保・転送もこの型を前提とします（[variables.cpp:312](/home/sano/work/forge/solver_density_cuda/variables.cpp:312)、[347行](/home/sano/work/forge/solver_density_cuda/variables.cpp:347)）。5配列の確保サイズだけでは runtime の型切替を実現できません。

   出力は `vector<flow_float>`（[output.cpp:156](/home/sano/work/forge/solver_density_cuda/output/output.cpp:156)）、保存量の読込みは `vector<geom_float>`（[variables.cpp:683](/home/sano/work/forge/solver_density_cuda/variables.cpp:683)）です。device だけ倍精度にすると、restart で蓄積した下位ビットを失います。

   **対案:** 型付きの専用格納領域、初期化・転送・解放、FP64 checkpoint、旧 FP32 入力からの初期化を設計に含めてください。連続実行と中断再開の比較をゲートに追加してください。

5. **Major — 既定 ON が未対応経路と他物理モデルへ及ぶ条件が未定義**

   **根拠:** §2 は `unsteady`／dual-time を対象外とする一方、§4.4 は無条件の既定1です。共有配列は陽解法・dual-time でも使われ、commit も別実装です（[update_d.cu:464](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:464)）。

   また、対象外とした SST 自体が平均流の `roe` を変更します（[ransTransport_d.cu:175](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransTransport_d.cu:175)）。周期状態コピー（[periodicNode_d.cu:119](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:119)）や等温壁ピン（[nodeWallDirichlet_d.cu:85](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85)）も整合対象です。「スカラー配列を倍精度化しない」と「連成経路を扱わない」は別です。

   **対案:** 初期版は既定0とし、対応する GPU・時間積分・物理モデルの組合せを明示してください。非対応で明示 ON なら起動時に拒否し、要求値と実効値を記録してください。既定変更は対応経路の検証後に判断すべきです。

6. **Major — G3 は現行検証方針と不整合で、不合格の解釈も誤っている**

   **根拠:** 主候補の `case/08.bump` は、現行手順で「まだ回帰の基準には使えない」と明記されています（[verification/README.md:23](/home/sano/work/forge/procedures/verification/README.md:23)）。また、「run-to-run ばらつき以内」は比較 run・報告量・数値閾値が未指定です。

   §4.3 の「ビット不変」が成立し得るのは、**同じ FP32 入力に対する同じ演算**についてです。蓄積によって入力状態が変われば、流束・リミッタの結果も変わります。さらに `updateGuardScale` の倍精度化は判定自体を変え得ます。したがって、§6 の「報告量が動いたら実装の切り方が誤り」という診断は成立しません。

   **対案:** G1 は OFF 経路の数値配列比較、G3 は ON 経路の精度・保存性・収束評価として定義し直してください。node の SST は `case/36`・`case/48`、対応を主張する軸対称 TP・凝縮は `case/44`、共有周期経路は `case/09` を選び、量ごとの許容値を登録してください。scalar/block 両 commit も試験対象です。

   cell は現行規則が実行回帰を禁止しています（同文書 [63行](/home/sano/work/forge/procedures/verification/README.md:63)）。コード上の整合を確認し、実測保証外と明記するのが適切です。メッシュ品質・同一 IC・段階起動区間の記録も検証前提へ追加してください。

7. **Major — +0.1% と +20 B/CV の見積りは、提示された設計から成立しない**

   **根拠:** §4.1 の演算回数には、FP64 化する `updateGuardScale`、変換、増えた読出し・書込み、baseline コピーが含まれていません。`Q` を倍精度配列へ置換するなら、その読出し側全体にも影響します。演算ピークで割った値は step 時間増加の予測ではありません。

   §4.5 のメモリ増分も、`Q_N` を倍精度化するか、別 accumulator と FP32 状態を併存させるかで変わります。前述の保持問題を解消する前に G5 を固定できません。

   **対案:** 格納方式を確定してから配列本数と転送量を再積算してください。G4 は同一 native 環境でウォームアップ、反復測定、ばらつきを含めて判定し、小規模 case/56 と代表的な大規模ケースで確認してください。

8. **Minor — 「全 CV が ½ ULP 未満」は実測と一致しない**

   **根拠:** `run_0017_resro_dump/res_1000.h5`、`z/W>10` の16,607 CVでは、最大 `|dq|/ULP` は **0.63836**、float32 加算で値が動く CV は **7個、0.0422%** でした。

   **対案:** 「ほぼ全 CV」に修正し、対象 snapshot・抽出条件・最大値・該当数を併記してください。丸め落ちが大半を占めるという主張は維持できます。

**推奨は、既存 FP32 配列を演算用に維持し、専用 FP64 accumulator を正本とする opt-in 案へ計画を改訂することです。** 全カーネルのポインタ型変更を避けつつ、更新・拘束・checkpoint の同期箇所を明示できます。ただし追加メモリは FP64 5本なら **40 B/CV** で、現行見積りの修正が必要です。

優先順は、①状態の保持・書込み・restart 契約、②対応範囲と OFF 保証、③微小増分保持試験、④三者比較と G2/G3 の再定義、⑤性能評価です。依頼どおりファイルは変更しておらず、この提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
