# codex レビュー: tooling-nozzle-core-grid (plan)

- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `3de780a4` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m3
- **focus**: §4.4・§4.5 の事前登録 (格子 4 本・共通起動・到達後 4 万 step・τ と u・J0〜J4) が G1 の採否を決めるのに足りるか。特に §4.5 の『NOT CONVERGED の扱い』で諮問の推奨を一部却下した理由が妥当か (time_integration-line-implicit-speed.md §6.15 のユーザ決定)。§4.2 の分布の実装上の落とし穴 (壁法線の層との組み合わせ、c_i の x 方向の滑らかさ)。読むのは plans/active/tooling-nozzle-core-grid.md・design/forge_design/meshing/mesh2d.py の 20〜268 行・plans/active/time_integration-line-implicit-speed.md の 417〜436 行と 547〜600 行・notes/reviews/2026-10-10-core-grid-plan-diagnose.md・case/45.isobutane_m6_d155/cold_series.py の 1〜60 行に限る。res_*.h5・*.log・plans/README.md は読まない
- **extra**: `notes/reviews/2026-10-10-core-grid-plan-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
G1 の試験実装は妥当ですが、現行の J0・J1 では生産候補の採否を十分に判断できません。  
水準による停止は維持し、未収束時の採否条件、格子感度の相殺、壁法線変換後の検査を実装前に修正してください。

指定された閲覧範囲だけで評価しました。run 数値・VERDICT は文書記載値であり、独立検証していません。`plans/README.md`・`accepted/`・仕様／運用文書は閲覧対象外のため、計画の重複やリポジトリ全体との整合性は未確認です。以下の反例計算は実行しましたが、forge は起動していません。

1. **Major — 「NOT CONVERGED の扱い」の却下理由は、採否条件を緩める根拠として不十分です。**

   **根拠:** [速度計画:419–435](/home/sano/work/forge-coregrid/plans/active/time_integration-line-implicit-speed.md:419) の決定は、point の停止値を正解扱いする E2 を外し、水準までの時間を測るものです。格子変更後の反復バイアスを許容差内と認定した決定ではありません。実際、同計画の [557–565行](/home/sano/work/forge-coregrid/plans/active/time_integration-line-implicit-speed.md:557) では、M64 は水準に到達しても ω 残差中央値が B0 の約7倍で、採用を保留しています。

   [本計画:151–159](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:151) の u は時間変動しか測りません。例えば θ(t)=1＋10⁻⁸t という減衰しないドリフトでも、2万stepの変化は約0.020 %、4万stepの u は約0.040 %となり、θ の水準・J0 を通ります。4万stepという長さだけでは定常性を保証できません。`check_quasisteady.py` による判定も登録されていません。

   **対案:** 水準到達時間の比較は維持します。一方、生産候補に進めるには、全採否量の `check_quasisteady` VERDICT と、**同じ格子・同じ離散方程式で反復条件を変えた継続計算による感度確認**を追加してください。測定した反復感度にも許容差の予算を割り当て、未評価なら J1 を保留します。これは真の反復誤差上限ではなく、運用上の頑健性確認です。  
   **「NOT CONVERGED を一律失格にしない」は認めますが、「現行も未収束だから注記だけで通す」は認めません。**

2. **Major — Gc を追加しても、J1 は格子変更の効果の相殺を見逃します。**

   **根拠:** [本計画:153–158](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:153) は G1−G2 だけで合否を決め、Gc は不合格後の説明にしか使いません。τ＝1 %、u＝0 として、G1＝G2、Gc だけが両者から1.5 %離れる反例でも J1 は通ります。

   また、[133行](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:133) の効果分離は厳密ではありません。q と nj を変えると、解かれる c_i も変わります。例示入力 f＝1.4×10⁻⁵で計画式を解くと、Gc の c＝0.011197 に対して G2 は0.008998、**主流間隔も約19.6 %縮みます**。Gc→G2 は壁近く「だけ」の変更ではありません。

   **対案:** J1 に、各量について G1→Gc→G2 の変化の絶対値の和と時間変動の余裕を含める条件を追加してください。例えば  
   `|G1−Gc|＋|Gc−G2|＋u_G1＋2u_Gc＋u_G2 ≤ τ`。  
   J2 では符号付きの差と実際の c_i 分布を示し、「パラメータ変更への感度」として扱ってください。4本はこの感度試験には有用ですが、誤差上限や完全な原因分離を与える系列ではありません。

3. **Major — §4.2 の間隔上限は、壁法線変換後の物理的な間隔上限ではありません。**

   **根拠:** [mesh2d.py:238–244](/home/sano/work/forge-coregrid/design/forge_design/meshing/mesh2d.py:238) は、分布生成後に d と β(d) を使って X・R の両方を変更します。したがって g_k は変換前の間隔であり、遷移層内の半径方向間隔・辺長とは異なります。

   既存の変換式を、局所半径1・傾斜40°の直線壁、f＝3.5×10⁻⁷、G1 の条件に適用した数値例では、c＝0.029629 に対して最大半径方向間隔は **0.034588**、最大辺長は **0.037174**でした。これは case/45 の実測ではありませんが、`c_max=0.03` が物理間隔の保証にならない反例です。

   **対案:** c_max を「変換前の計算座標での上限」と明記し、変換後の全格子について物理的な半径方向間隔、辺長比、符号付き面積、AR・skew を検査してください。軸近傍では β＝0 なので、軸の上限とは分けて報告できます。  
   また、[本計画:149](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:149) の y₁⁺を記録だけにするのは不十分です。同じ第一層座標でも壁応力が変われば y₁⁺は変わるため、局所超過率・位置を採否時の壁解像確認に含めてください。

4. **Minor — c_i の隣接比1.006だけでは、x方向の滑らかさを確認できません。**

   **根拠:** [本計画:118](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:118) の評価に対し、第一セル表は [mesh2d.py:90–91](/home/sano/work/forge-coregrid/design/forge_design/meshing/mesh2d.py:90) で log 線形補間です。さらに、非飽和区間数を m とすると、固定した m の範囲では  
   `c＝[1−f(1＋q＋…＋qᵐ⁻¹)]/(n−m)`。  
   m が切り替わる点で c は連続でも導関数は一般に連続ではありません。計画式の数値例では、ほぼ同じ c のまま dc/df が約−17387から−13990へ切り替わりました。

   **対案:** 表の接続点と飽和区間数の切替点で、非一様な Δx を考慮した格子線の傾き変化・局所品質を検査してください。c_i を後から平滑化すると列の和や第一セルを壊すため、まず実座標の検査を追加するのが適切です。

5. **Minor — 後処理の新しい nj への対応と、判定量の定義が未完成です。**

   **根拠:** [cold_series.py:25–33](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cold_series.py:25) は `XC.mesh_info(len(ro))` と `reduce_fields` に依存します。閲覧範囲では nj＝92・122・170への対応を確認できず、[残作業#2](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:209) にもその検証がありません。また、[145–153行](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:145) は「波」を量に含めながら τ を定義せず、「全量」の範囲が曖昧です。

   **対案:** 長時間投入前に4格子すべてで後処理を通し、既知分布による積分・抽出の確認を行ってください。波を判定対象にするか記録だけにするか、相対差の分母、欠測・非有限値の扱いも確定させます。これは速度試験より先に必要です。

6. **Minor — 総時間の上下限の作り方が未登録です。**

   **根拠:** [本計画:139–141](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:139) は2500step間隔の到達判定と1000step×3の単価測定を定めていますが、[J1・J3](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:154) の時間上下限への変換規則がありません。32 ms/stepなら出力間隔だけで約80秒あります。[速度計画:433](/home/sano/work/forge-coregrid/plans/active/time_integration-line-implicit-speed.md:433) でも、この分解能が判別を左右しています。

   **対案:** 起動各段・本段の単価幅、出力時間、到達検出の分解能を含む上下限の式を事前登録してください。追加4万stepを除外するなら、評価名を「水準到達までの時間」とし、検証を含む実消費時間は別記します。

**推奨は、G1と4格子の構成を維持し、上記の順に採否条件・物理格子検査・後処理・時間区間を確定してから実装へ進むことです。** 固定 q・可変 c の分布自体は、節点を減らして主流へ配分する目的に合っています。FP64で得た採否を float32・cell中心・周期ケースへ一般化する根拠は、この試験にはありません。

ファイル変更はなく、**plan未反映**です。

指摘数: Critical 0 / Major 3 / Minor 3
