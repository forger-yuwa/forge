# codex レビュー: condensation-kantrowitz-gamma-twophase-sonic (plan)

- **plan**: [`plans/active/condensation-kantrowitz-gamma-twophase-sonic.md`](../../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)
- **stage**: `plan`
- **date**: 2026-09-10
- **commit**: `5db69d93` (feature/sern-design)
- **codex**: effort `high`, 4.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: 1 回目レビュー (NO-GO, C1/M5/m1) の指摘を全件採用してスコープを TP carrier H2O 非平衡/緩和形に縮小した改訂版。各指摘への対応 (§2 対象外の明示, §4.2 固定 g,Y frozen 近似としての LHS, §3/§7 境界 ghost, §6 凝縮残差列と時系列定常判定, node/cell・凝縮OFF・反復ノイズの回帰, 合否条件) が十分か。残る Critical/Major があれば根拠つきで
- **extra**: `methods/condensation.md`, `notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-plan.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
固定 `g,Y` の音速式と N2/CPG の除外は妥当で、前回の Critical は解消されています。  
ただし、既定適用範囲と検証に Major が残り、「前回指摘を全件解消」とは認められません。

1. **Major — 未対応境界への自動適用が残っている。緩和形も検証されていない。**

   **根拠:** [plan:104](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:104) の自動選択条件は境界を確認しません。したがって、対象外と書いた `outlet_statPress`・等温壁でも新音速が有効になります。

   [boundaryCond_d.cu:660](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:660) は亜音速出口 ghost を全蒸気 EOS で再構築し、[nodeWallDirichlet_d.cu:78](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:78) は壁の `roe/P/sonic` を上書きします。未対応と明記するだけでは、既定変更から除外できません。

   また、自動選択は `condEquilibrium:1` を含みますが、[検証表:166](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:166) は全て非平衡です。

   **対案:** 今回の自動適用を **TP carrier H2O・`condEquilibrium:0`・検証済み境界構成**に限定する。境界読込後に適用可否を確定し、解決値と理由をログへ出す。緩和形と未対応境界は今回の既定変更から外し、選択条件の回帰試験を追加する。

2. **Major — 「参照と同等以上」では、未収束のまま合格できる。残差列追加だけでも不十分。**

   **根拠:** [plan:188](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:188) は VERDICT の掲載を要求しますが、`PASS` を合格条件として明記していません。

   `case/16.nozzle_wys/run_0331_condfix_head_ref/` にツールを再実行した結果は、保存済み [CONVERGENCE_VERDICT.txt](/home/sano/work/forge/case/16.nozzle_wys/run_0331_condfix_head_ref/CONVERGENCE_VERDICT.txt) と同じです。

   ```text
   NOT CONVERGED (still converging — run more steps)
   rms_roe: final=1.02e-3, drop=2.6dec, falling
   ```

   凝縮列をメモリ上で検査対象に追加しても総合判定は `NOT CONVERGED`。凝縮４列はそれぞれ約 `4.9 / 5.6 / 6.4 / 7.1` 桁低下しており、この参照の未収束原因はエネルギー残差です。

   さらに、[check_convergence.py:90](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:90) の初期残差ゼロの特例は低下桁数を検査しません。実際の `analyze()` に合成系列 `[0,1,1,…,1]` を与えると、**低下ゼロなのに `ok=True`** になりました。これは run 実測とは別の、判定関数の再現試験です。

   **対案:** 効果比較に使う参照・新モデル双方で、全対象残差の `PASS` を必須にする。初期値ゼロでも途中から活動した列はピーク基準で判定し、全期間ゼロの列だけを除外する。R0 は延長してから比較基準にする。

3. **Major — 末尾２枚の判定では、報告量の定常性を保証できない。**

   **根拠:** [plan:163](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:163) と、確認時点の [compare_condfix.py:69](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:69) は末尾２点だけを比較します。振動の同位相や極値付近を拾えば、動いている量にも `STEADY` を返せます。

   実装上の判定対象は壁偏差の**平均**・onset・出口 `g` だけです。局所壁圧の変化は平均で相殺され得ます。表示する `M_exit/c_exit/h0err` と壁圧３地点は合否に入りません。また、[同:41](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:41) の onset は最初の閾値超過点なので、セル内の移動を検出できません。

   **対案:** `check_quasisteady.py` に対象量の抽出を追加し、複数スナップショットの末尾窓でトレンド・振幅を判定する。全報告量に許容値を設定し、onset は閾値交差を補間する。`DRIFTING`・`OSCILLATING`・データ不足を区別し、不合格時は非ゼロ終了にする。

4. **Major — cell run を追加しても、現行の後処理では評価できない。**

   **根拠:** [compare_condfix.py:22](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:22) は `MESH/COORD` を `VALUE` と同じ位置・順序の配列として扱い、[同:32](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:32) は `wall_dist<=0` で壁を抽出します。

   実データ `case/16.nozzle_wys/run_0197_user_cell_sst_cond_cont/res_24000.h5` では、

   - `MESH/COORD`: **83,280 点**
   - `VALUE/P`: **41,174 セル**
   - `wall_dist<=0`: **0 セル**

   同ファイルに実際の `lines()` を適用すると、次で失敗しました。

   ```text
   ValueError: operands could not be broadcast together
   with shapes (41174,) (83280,)
   ```

   **対案:** cell は入力メッシュのセル中心と境界面情報を使い、node と同じ物理位置へ補間して比較する。座標・場の対応と壁抽出の検証を、R4/E4 実行前の独立ステップにする。

5. **Major — frozen 近似として使う実 Jacobian の検証まで後送りしている。**

   **根拠:** [plan:149](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:149) は、固定 `g` の流束 FD と float32 行列作用を、厳密な `ρg` 列の実装と一緒に後続へ送っています。

   しかし今回変更する `gamma/sonic` は、[block_dplur_jacobian_d.cuh:38](/home/sano/work/forge/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:38) の `κ` と、[同:54](/home/sano/work/forge/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:54) の `χ` に直接入ります。音速スカラーの FD 一致と `sonic` 出力照合だけでは、`gam_array` の配線や実ブロックの作用を確認できません。残差が下がることも、その代用にはなりません。

   **対案:** 厳密な `ρg` 列の実装は後続のままでよいので、**今回使う固定 `g,Y` ブロックの試験は戻す**。実際の `accumulate_split_jacobian_cf` に wet 状態を渡し、固定 `g,Y` の流束 FD と比較する。float32 の方向微分と、実カーネルの `gamma`・`sonic` 双方を照合する。

6. **Minor — 単体試験の合格条件が `g=0` で矛盾する。**

   **根拠:** [plan:157](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:157) は掃引に `g=0` を含めながら、全点で `c₂φ<c全蒸気` を要求します。一方、直後の試験は `g=0` で厳密一致を要求しています。

   **対案:** `g=0` は等号、`g>0` の対象域では厳密不等号に分ける。

7. **Minor — carrier の「真値」に対する上下限という表現は、まだ強すぎる。**

   **根拠:** [plan:84](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:84) と §10 は `J_pure≤J_carrier≤J_iso` を真値の範囲として説明しています。しかし引用論文は、熱除去の改善と逆向きに作用する carrier の `pV` 仕事も区別しています。[Wedekind et al., 式9–12](https://arxiv.org/pdf/0804.1516)

   **対案:** この不等式を「同一の核生成障壁・前因子を固定し、衝突による熱除去だけを追加したモデル間の関係」に限定する。実際の carrier 中の核生成率の保証範囲とは呼ばない。

**推奨は、検証済み構成への限定導入を維持し、上記を修正してから実装することです。** 実装前の優先順は、①自動適用条件の制限、②残差判定と `PASS` 条件の修正、③全報告量の定常判定、④cell 後処理の整備、⑤実 Jacobian 試験の復帰、です。

目的は未解決の不整合を正しく捉えています。[plans/README.md](/home/sano/work/forge/plans/README.md) と関連 accepted plan を確認しましたが、完成済み機能との重複ではありません。固定 `g,Y` の導出は既存 EOS の `de/dT` と整合し、厳密な連成ブロック全体の実装を今回に要求する必要もありません。

ファイルは変更していません。指摘・推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
