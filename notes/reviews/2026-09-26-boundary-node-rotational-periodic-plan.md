# codex レビュー: boundary-node-rotational-periodic (plan)

- **plan**: [`plans/active/boundary-node-rotational-periodic.md`](../../plans/active/boundary-node-rotational-periodic.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `df706869` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m2
- **focus**: (i) 変換の置き場所と向き (gather R^T / broadcast R、勾配テンソルの相似、DPLUR 対角の相似) に符号・転置の誤りが無いか。(ii) リミタで速度 2 成分の seam gather を止める方針の単調性への影響と、代替 (円筒成分 bound) を後回しにしてよいか。(iii) 角の bookkeeping (重み付き union-find、閉ループ、float cos/sin の精度)。(iv) 並進 (type 0) のビット不変が『同じカーネル』で担保されているか。(v) スカラー経路 (species/k-ω/凝縮/受動) を無変更にしてよい根拠。(vi) 同値性試験 (回転コピー全周と節点一致) と成立条件が測る前に固定されているか。(vii) slauWallNormalChi mask の seam 整合が回転でも成り立つか。(viii) 軸を含むセクタ・軸対称併用を対象外にした線引き。実装箇所 solver_density_cuda/cuda_forge/periodicNode_d.cu と mesh/mesh.cpp:512-775 を確認してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
seam 合算経路を拡張する目的と、回転の基本式は妥当です。  
ただし、陰解法の fold・リミタ・精度条件・同値性試験には実装前の修正が必要です。

レビュー基準は開始時の commit `df706869` と提示された plan です。調査中に別作業の変更を検出したため、それ以降は基準 commit のコードを参照しました。ファイル変更はしていません。

目的は未解決事項を正しく捉えています。既存の cell 周期修正とは重複せず、親 plan §4.5.8 の実装タスクとして妥当です。`R_m` を root→member と定義すれば、残差・スカラー勾配の gather `R_mᵀ`、速度勾配の `R_mᵀ G_m R_m`、対角ブロックの相似変換に転置の誤りはありません。問題は、それらを適用する**対象と順序が不足していること**です。

1. **Major — 「既存の対角 fold に回転を挿入する」という実装前提が誤っています。**

   **根拠:** 親 plan は、対角と LU 和の厳密 fold を「未実装」と明記しています。[discretization-median-dual-3d.md:152](/home/sano/work/forge-sern-design/plans/active/discretization-median-dual-3d.md:152)  
   実コードも CV ごとに対角・`neighbor_accum` を構築して solve し、その後 root の `dq` をミラーしています。[timeIntegration_d.cu:824](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:824)、[main.cpp:1591](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1591)

   対角だけを合算すると、member 側の非対角寄与が欠けます。また現在の対角には既に**合併体積**による `V/dt_local` と物理時間項が入るため、そのまま足すと時間項を重複計上します。

   **対案:** §4・§5を「回転 group の行縮約を新設する」に修正してください。部分 CV の空間対角を相似変換して合算し、各 sweep の近傍寄与も `Σ R_mᵀ neighbor_accum_m` として集約する。時間項は group に一度だけ加え、拘束行処理・対角キャッシュ・solve・`dq` broadcast の順序を明記する必要があります。小規模な明示行列との **1 sweep の比較**を追加し、残差履歴だけで正否を判定しないでください。

2. **Major — 速度 extrema の gather を止めるだけでは、リミタの回転整合が成立しません。円筒成分 bound の後回しには反対です。**

   **根拠:** 現実装は extrema の gather に加え、最後に `limiter_Q` 自体を成分別に min-gather しています。[limiter_d.cu:76](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:76)、[limiter_d.cu:90](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:90)

   90° 回転でも `ψ_y` と `ψ_z` は交換が必要です。例えば root の係数が `(0.2, 0.8)` なら、回転先は `(0.8, 0.2)` です。同名成分を min-gather すると両側 `(0.2, 0.2)` となり、全周側とは別の再構成になります。一般角では、成分別の制限操作そのものが回転と可換ではありません。

   片側 bound は通常より狭くなるため、直ちに「新極値を生む」とは断定できません。しかし、過剰制限や seam 依存の誤差が残ります。V6 の**節点値が全周の範囲内**という判定では、再構成面値の逸脱も検出できません。

   **対案:** 回転と整合する局所基底で、全近傍の bound と最終再構成を一体設計してください。推奨は円筒成分による方式で、全周対照にも同じ方式を適用することです。V6 は既存の面再構成検査を拡張し、局所 bound、回転後の面状態、滑らかな旋回場の精度を確認してください。[limiterPeriodic_d.cuh:168](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiterPeriodic_d.cuh:168)

3. **Major — 角度・座標の精度経路と、U0/U1 の閾値が両立していません。**

   **根拠:** `flow_float` と `geom_float` はともに float32 で、`dtheta` は YAML 読込時点で float32 に丸められます。[flowFormat.hpp:6](/home/sano/work/forge-sern-design/solver_density_cuda/flowFormat.hpp:6)、[boundaryCond.cpp:60](/home/sano/work/forge-sern-design/solver_density_cuda/boundaryCond.cpp:60)

   読み取り専用の数値確認では、次の値になりました。

   | 確認 | 誤差 |
   |---|---:|
   | float32 の `π/2` による単位半径の座標対応 | `4.371139e-8` |
   | 上記を4回合成した閉ループ角 | `1.748456e-7 rad` |
   | double から作った float32 の45°回転係数の直交性 | `3.422854e-8` |

   したがって、**後から角を double に格納しても精度は戻りません**。また U0 の `1e-12` を device の float32 行列にも要求するなら不適切です。

   **対案:** node 回転用の角度を入力から double で保持し、元の double 幾何での検査と、float32 格納後の検査を分けてください。`H` と誤差ノルムも定義する必要があります。閉ループは `abs(remainder(Δθ, 2π))` で検査し、逆向き BC、root 交換、連鎖、非90°角を含めてください。U0 は double の bookkeeping と実 CUDA 変換の許容を分離してください。

4. **Major — 勾配 gather の成立条件が Green–Gauss に限定されることを、受付条件に反映していません。**

   **根拠:** 現 gather の和が正しい理由は、各部分勾配が合併体積で割られていることです。[periodicNode_d.cu:166](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:166)  
   一方、`gradLSQ=1/2` は片側 stencil から勾配を直接解き、その後も呼出側が同じ gather を実行します。[calcGradient_d.cu:946](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:946)、[main.cpp:1459](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1459)

   両側が線形場の勾配をそれぞれ正しく再現した場合、回転して足すだけでは **2倍**になります。

   **対案:** 今回は回転周期を `gradLSQ=0` に限定し、非対応設定を起動エラーにすることを推奨します。軸を含むセクタ・軸対称併用の除外も妥当ですが、文章上の除外だけでなく負例試験を追加してください。`lineImplicit`・前処理付き陰解法なども、対応範囲を受付表で固定する必要があります。

5. **Major — V2 の「同 step・同 CFL・ノイズ床×2」は、全周との同値性を保証する条件になっていません。**

   **根拠:** `dt_local` は各 CV が持つ面集合の最大 CFL から決まります。部分 CV と全周の統合 CV では、同じ設定値でも同じ刻みになる保証がありません。[setDT_d.cu:157](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/setDT_d.cu:157)、[setDT_d.cu:216](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/setDT_d.cu:216)

   また `check_field_regress.py` は配列形状が同じことを要求し、座標対応・全周からの抽出・ベクトル回転を行いません。[check_field_regress.py:104](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_field_regress.py:104)  
   同一 mesh 内の反復ノイズは、seam 合算と全周内部合算の**異なる演算順による系統差**の上限でもありません。

   **対案:** V2をまず固定状態での空間作用素比較、その次に共通時間刻みでの更新比較へ分けてください。測定前に、以下を固定する必要があります。

   - 回転コピー後の節点統合、対応 ID、双対面・合併体積・壁距離の対応。
   - 全周の IC・BC・ソースがセクタ角の回転対称性を満たす条件。
   - 全周抽出と変換後の比較量、ゼロ近傍の絶対許容、丸め誤差の許容。
   - 非90°回転、非対角の速度勾配、非零の横運動量を持つ実カーネル試験。

   V1 の純軸流だけでは、横運動量の回転符号やテンソル変換を十分に励起できません。

6. **Major — V3 の保存収支式が誤っており、「各種スカラー対応」の検証も不足しています。**

   **根拠:** SST の残差には生成・消滅・交差拡散ソースが入ります。[ransSource_d.cu:248](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:248)  
   そのため、ソース込み残差について「総和＝境界流束」は成立しません。また gather 後は全 member に同じ残差が broadcast されるため、全節点を単純合計すると重複計上になります。[periodicNode_d.cu:62](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:62)

   **スカラー値の転送を恒等のままにすること自体は正しい**です。既存 species-DPLUR には近傍寄与の合算もあります。[speciesTransport_d.cu:529](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:529)  
   ただし、その根拠はスカラー勾配・面流束・ソース評価まで無変更でよいことを意味しません。

   **対案:** 輸送のみの収支と、体積ソース・拘束・更新補正込みの収支を分け、合算後は root のみを数えてください。float32 残差に相対 `1e-12` を要求せず、流束絶対値総和を尺度とする精度依存の許容を事前設定してください。さらに凝縮モーメント・受動種・遷移量を含む登録配列の転送試験と、非零ソースの小規模試験を追加してください。

7. **Major — V4/V7 は収束・準定常の受入ゲートとして不足し、V4 不合格時の処置も原因を決め打ちしています。**

   **根拠:** V7 の「全周の2倍以内・RISINGなし」では、両者が未収束のプラトーでも合格できます。`check_convergence.py` は残差低下も要求します。[check_convergence.py:337](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:337)  
   V4 には壁圧・`C_f` の準定常判定、局所壁解像、段階起動と判定区間がありません。[対象plan:103](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:103)

   `slauWallNormalChi` 自体は `u·n` から計算されるため、状態と法線が正しく回転し、mask が対応すれば回転不変です。[convectiveFlux_slau_d.inc.cuh:549](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:549)  
   従って、V4 不合格を直ちに「回転周期では auto→0」と処理する根拠はありません。

   **対案:** 両 run の全残差 `PASS`、各 run の報告量の準定常 VERDICT、局所 `y₁⁺`、同一設定区間を受入条件にしてください。flag の評価はセクタ／全周 × flag0／flag1 の比較と、対応面の流束照合で切り分けてください。

   なお V5 の参照先について実行した結果は、次のとおりです。再現用入力の確保も前提タスクに必要です。

   ```text
   check_convergence.py case/39.periodic_hills/run_0007_coarse_rans
   [case/39.periodic_hills/run_0007_coarse_rans] NO residual_history.csv
   OVERALL: CHECK FAILURES ABOVE
   ```

8. **Minor — 「同じカーネル＝並進のビット不変」という保証は強すぎます。**

   **根拠:** 既存 gather は `atomicAdd` を使用し、反復間のビット差が既存ツールでも明示されています。[periodicNode_d.cu:27](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:27)、[check_field_regress.py:4](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_field_regress.py:4)

   **対案:** 保証を「type0 の演算経路を維持する」と定義し、旧カーネルだけでなく root 配列・合併体積・面対応も旧版との一致対象にしてください。面流束のビット比較と、集約後の場のノイズ床比較を区別する V5 の構成は妥当です。

9. **Minor — 新設予定の case 番号が既存ケースと衝突しています。**

   **根拠:** `case/57.transition_flat_plate` が既に存在します。[README.md:1](/home/sano/work/forge-sern-design/case/57.transition_flat_plate/README.md:1)

   **対案:** 未使用番号を確保してから §5・§6・run 索引の参照先を統一してください。

**推奨は、既存の seam 合算構造を維持し、上記を反映した計画へ改訂してから実装することです。** 実装前の優先順は **①行縮約と時間項 → ②回転整合リミタ → ③入力精度と対応検査 → ④受付条件 → ⑤同値性・保存・収束ゲート**です。リミタを V6 の結果待ちにする進め方は推奨しません。

本レビューは **plan 未反映**です。依頼どおりファイルを変更していません。

指摘数: Critical 0 / Major 7 / Minor 2
