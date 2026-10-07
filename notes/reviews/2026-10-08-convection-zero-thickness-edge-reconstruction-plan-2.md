# codex レビュー: convection-zero-thickness-edge-reconstruction (plan)

- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **stage**: `plan`
- **date**: 2026-10-08
- **commit**: `7fa7745c` (feature/sern-design)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
幾何で限定した速度再構成の抑制は、判別実験として実装する価値があります。  
ただし、§6.2 は冷却が継続・移動していても合格でき、恒久採用と旧評価の持ち越しの条件も混在しています。以下を実装前に修正してください。

目的の同定と基本構造は妥当です。`plans/README.md` と関連 `accepted/` を確認しましたが、既存の `slauWallNormalChi`・多成分面組成整合化とは変更対象が異なり、解決済み案件の重複ではありません。速度二乗・面エンタルピー計算前の介入位置と、内部面流束の逆符号加算は適切です。block-DPLUR はセル状態による近似ヤコビアンなので、今回 LHS まで変更する必要はありません。[SLAU:279](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:279)、[同:661](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:661)、[時間積分仕様:276](/home/sano/work/forge-sern-design/methods/time_integration/implementation.md:276)。

実測資料では、`case/46.sern_design/run_1066_r5h_m10_B_ledger100` の call 10 に、内部エネルギー更新 **−7.755×10³ J/kg**、EOS 差 **＋4.188×10⁻⁴ J/kg** が記録されています。局所処置を検証する根拠として十分です。ただし、対象 run の生 HDF5・残差 CSV・個別 VERDICT はこの checkout にないため、収束・定常性は再認定していません。[持ち帰り帳簿:16](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R5H_M10_LEDGER100.txt:16)。

1. **Major — 「孤立した冷点」の合格条件が、冷点の不存在も冷点移動の停止も保証しない。**

   **根拠:** §6.2 #2 は `min(T in R_TE) ≥ min(T outside R_TE)`、温度場の合否は #1・#2 だけです。領域外にもっと冷たい点が生じるほど合格しやすくなります。また、処置対象 E は自由側端まで含むのに、検査領域は後縁付近の x 範囲だけです。[plan:39](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:39)、[同:87](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:87)、[同:101](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:101)。

   読み取り実行した**合成系列による反例**では、局所最低温度が 20000 step で **130→90 K**、領域外最低温度が 80 K の場合、#2 は全保存時点で合格します。一定の力係数は `STEADY` ですが、同じ `check_quasisteady.classify()` で温度系列を判定すると **`DRIFTING`** です。力の定常性だけではこの穴を閉じられません。

   **対案:** 領域外の全域最低温度との比較を廃止してください。後縁・自由側端・マスク外周を覆う固定物理領域で、局所温度欠損、低温領域の体積と位置を追跡し、その閾値と参照場を事前登録すること。局所温度指標にも `check_quasisteady.py` を適用し、`DRIFTING` は恒久採用不可としてください。これは「正しい温度」の認定とは分けて扱えます。

2. **Major — 「リング数を選んで同じ物理幅にする」格子感度試験は、現格子では定義不足。**

   **根拠:** §6.2 #6 は、g4 のリング数だけで g3 と同じ物理幅を実現するとしています。しかし既存 g3→g4 は `first_wall_frac` が **1/4**、`nz_out` が **15→32**、x ブレンドも変更される異方的な変更です。単一の整数リング数では、壁法線・流れ方向・側端方向の広がりを同時に合わせられません。[plan:96](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:96)、[3D plan:1741](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1741)。

   **対案:** 生産候補の「2 リング固定」は維持し、幅感度の補助試験には、g3 の処置領域から**流れを参照せず事前に固定した物理的な領域**を両格子へ適用してください。方向別の広がり・対象体積・領域の対応誤差を記録すること。リング数だけを変えるなら「同じ物理幅」とは呼ばず、幅の不一致を含む感度試験としてください。

3. **Major — 床カウンタの「床分岐が作動した」という定義を、TP の実装に合わせて具体化する必要がある。**

   **根拠:** 今回の TP 温度床は外側の `tMin` 分岐ではなく、`DEPVAR_TMIN=50` を渡した温度反転関数の内部にあります。float Newton と double 研磨の**各反復**でクランプされ、最終的な `roe` は別の場所で再構成されます。分岐の呼出回数を数えるだけでは、物理状態の床違反件数・最終保存量補正と対応しません。[plan:57](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:57)、[dependentVariables_d.cu:172](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:172)、[thermo_d.cuh:660](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:660)、[同:748](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:748)。

   **対案:** 「Newton 中間反復の制限」と「更新された保存量が許容温度範囲外にあり、最終状態を床へ補正した事象」を分離してください。後者を**実節点・更新単位**で数え、ghost は別集計にすること。温度下限に対応する内部エネルギーと入力状態の比較、最終補正量、通常の丸め誤差の扱いを明記し、正常状態・床未満・下限近傍の試験を §5.3 に追加してください。

4. **Major — #4 の感度基準が、旧評価の持ち越し条件と修正自体の採用条件を兼ねている。**

   **根拠:** #4 は許容超過時に「修正の誤りと断定せず旧評価の持ち越しを止める」としています。一方、R7b 再開条件は **#4 合格が必須**、完了条件も全項目合格です。したがって、修正が必要なだけ評価値を変えた場合、旧評価を取り直しても完了できません。[plan:94](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:94)、[同:104](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:104)、[同:119](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:119)。

   **対案:** #4 の **5e−5／5e−4 は旧評価の持ち越しを判定する基準**として維持してください。超過時は旧評価・学習データを失効させ、修正後の基準列と設計差を取り直す分岐を明記すること。修正の採用は床・局所場・定常性・格子／領域感度で判断し、「旧結果に近いこと」を必要条件にしないでください。

5. **Major — 表にはなったが、#7・#9 と設計利得の最終判定はまだ事前登録として閉じていない。**

   **根拠:** #7 の「加重目的とその不確かさで判断」、#9 の「個別に評価」には、合否式・対象構成・比較量がありません。#5 の **5e−4／0.0125** の引用元は、有限厚模型の**設計差の格子感度**に対する暫定提案です。今回の処置有無による差の差に同じ数字を使っても、格子・領域・時間変動を含む総合的な判別能力の保証にはなりません。[plan:95](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:95)、[同:97](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:97)、[同:99](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:99)、[引用元:74](/home/sano/work/forge-sern-design/notes/reviews/2026-10-07-sern-coarse-design-mesh-diagnose.md:74)。

   **対案:** #5 の値は独立した暫定予算として扱い、#7 に重み・改善方向・不確かさの合成方法を記載してください。設計順位なら **改善量−総合感度幅 > 0**、0.002 を超える改善を主張するなら **>0.002** を要求する形です。各作動点の `C_M` 制約も別に判定してください。#9 は実際にマスクが発火する固定ケースと定量基準を指定し、該当ケースがなければ小規模な端付き検証ケースを登録してください。無効経路の `case/66` 回帰だけでは、有効経路の他格子への適用を検証できません。

6. **Minor — 現在仕様と参照先が plan に同期していない。**

   **根拠:** `methods/convection/implementation.md` は「対象が空なら旧経路」と記述していますが、plan は opt-in 時の空集合を起動時エラーとしています。また #6 が指す「親 plan §8」は未確定事項で、`C_M≤0.05` の正本は `tooling-nozzle-sern-3d.md` §8 です。[methods:72](/home/sano/work/forge-sern-design/methods/convection/implementation.md:72)、[plan:41](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:41)、[chain:547](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:547)、[3D plan:1974](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1974)。

   **対案:** 空集合の契約を同期し、許容値の参照先を正してください。#6 には推力2量・`C_L`・`C_M` の数値を省略せず記載してください。

**推奨は、幾何限定・速度のみ・opt-in の案を維持し、上記 1→2→3→4→5→6 の順で計画を修正してから、計測実装と §6.1 の判別 A/B に進むことです。** 周期・cell・非 SLAU・軸対称の拒否、計測を A/B より先に置く順序、同一格子／別格子で restart 手段を分ける方針は維持してください。短期 A/B の成功は恒久採用の承認ではありません。

ファイルは変更していません。修正提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
