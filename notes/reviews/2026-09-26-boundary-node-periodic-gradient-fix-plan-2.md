# codex レビュー: boundary-node-periodic-gradient-fix (plan)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `5c891e18` (feature/sern-design)
- **codex**: effort `high`, 4.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2
- **focus**: 実装の区切り (§5.1 #3・#4 実装済、正式試験はこれから)。実装 diff (commit HEAD の calcGradient_d.cu の lsqPre_pinv 切り出しと lsqPre_mergePeriodic、periodicNode_d.cu の periodicSeamMergeActive、ransTransport_d.cu の excludePeriodic と ransGradient 直後の gather、ransSource_d.cu / ransTransport_d.cu の F1 初回充填) が §4.1・§4.2 の設計どおりか、誤りが無いかを見てほしい。特に (i) 同値類の同定 (jroot と Δx 1e-4 h_min)、代表値 (root 優先の最初の incidence)、α=1/count、M_r の組立と係数の焼き込みが codex plan M1 の式と一致するか、incidence の ilp と CSR の対応、float 化の順序。(ii) group の定義 (root[c]!=c の member を root ごとに集める) で root 自身の incidence が漏れないか。(iii) k/ω の gather の位置と、periodicGatherArray_d_wrapper が和と broadcast を両方行うこと。(iv) g_sstF1Computed の扱い。(v) §6 の正式試験に入る前に直すべき点。予備確認: 線形場で継ぎ目の速度・k・ω 勾配が真値と一致 (plan §6.2)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
HEAD `5c891e18` の合併 LSQ・`k/ω` gather は、主要部分で §4 の設計どおりです。  
ただし、代表変位への置換は G0 の精度条件を満たさない反例があり、正式試験前に設計修正が必要です。

1. **Major — 同値類の代表変位への置換と、G0 の線形精度保証が両立していない。**

   根拠: [`calcGradient_d.cu:764`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:764) で各 incidence の変位を取得しますが、770–796 行では許容差内の変位を最初の代表値へ置換します。一方、実行時の場の差分は各 incidence の実値です（[`同:878`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:878)）。

   この処理を CPU 上で再現した、非退化の局所 stencil に反例があります。float32 座標から得る接線変位を
   
   `d₀ = float32(0.03) − float32(0)`  
   `d₁ = float32(100.03) − float32(100)`
   
   とすると、差は `1.22003e-6`、照合許容差は `2.99988e-6` なので統合されます。格納座標に対する線形場を与え、係数・積算を float32 にすると、真値 `1` に対して **`0.99997973`、root 交換後は `1.0000203`**。いずれも G0 の `1e-5` を超えます。これは CFD run の実測ではなく、当該式の局所数値試験です。

   **対案:** 同値類は重複数と `α` の決定に使い、行列と係数には各 incidence の実変位を一貫して使ってください。
   \[
   M_r=\sum_{mj}\alpha_{mj}w_{mj}d_{mj}d_{mj}^{T},\qquad
   c_{mj}=M_{r,\tau}^{+}\alpha_{mj}w_{mj}d_{mj}.
   \]
   完全な重複では現在の一意 stencil 式と同じです。この式では上記反例は両順序とも `1.000000` でした。§4.1・G2 の参照定義を更新し、この反例を G0 に追加すべきです。

2. **Major — §6 の正式試験には、まだ合否を一意に決められない箇所がある。**

   根拠: [`plan:115`](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:115) の G2′ は、対称 stencil で二次場の勾配が厳密になる場合、誤差比・収束次数が丸め床で未定義になります。R2 の粘性 SLAU は今回の勾配変更を直接読む重要な回帰ですが、保存量を「記録」するだけで失敗条件がありません。また「初期総量で正規化」は、初期全運動量がゼロの TGV にそのまま適用できません。

   R1 も `C_f` の測定位置・`x_r` の抽出規則が未指定です。既存ツールは周期丘のこれらを直接抽出せず、外部生成した時系列を渡す必要があります（[`check_quasisteady.py:354`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:354)、同377–415行）。

   **対案:** 試験開始前に、G2′ の非対称 stencil と誤差床の扱い、SLAU の保存誤差上限、運動量の非零基準量 `ρ₀U₀V` を固定してください。R1 は定常 SST-RANS を明示し、`C_f` の3位置・接線応力の定義・再付着点の選択規則を固定した CSV を `--series-csv` で判定する形にしてください。既存の段階起動・メッシュ品質・同一区間の収束判定という方針は妥当です。

3. **Minor — `g_sstF1Computed` は通常経路では有効だが、計画した初期化責務と異なる。**

   根拠: [`ransSource_d.cu:363`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:363) のフラグはプロセス全体で共有され、配列の寿命に紐付きません。[`ransTransport_d.cu:113`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransTransport_d.cu:113) の `static inited` も残っています。§4.2 が指定した「変数初期化時へ移す」は未実施です。

   現行 `main` は毎回 `ransBlendF1 → ransTransport` の順なので、**通常実行で初回 F1 が上書きされる問題は防げています**。ここを現行計算の破綻とは判定しません。ただし、同一プロセスで別の `variables` を扱う試験では、前の配列の計算済み状態を引き継ぎます。

   **対案:** 計画どおり配列初期化時に初期値を設定し、`buildScalarDescs` を初期化副作用のない関数にしてください。初回・2回目とも `F1≠1` となる入力で、輸送が読む値を検査してください。

4. **Minor — 予備確認の再現物と、正式試験の作業表が不足している。**

   根拠: この作業ツリーの `case/09.Taylor-Green/_g0_lsq_seam/` には修正前の集計テキストと YAML 2個だけがあり、参照する `Taylor-Green.h5`、生成・集計スクリプト、修正後出力がありません。保存設定も [`solverConfig.yaml:18`](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/solverConfig.yaml:18) では `model: none` で、SST 予備確認を再現できません。さらに [`plan:105`](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:105) の検証作業には G0/G1・G2′・R3 が明示されていません。

   **対案:** 入力生成・集計手順、旧新 revision、修正後出力を保存し、全試験を作業表に登録してください。今回、§6.2 の修正後数値は独立に再集計できていません。

重点項目のうち、次は問題ありませんでした。

| 項目 | コード照合結果 |
|---|---|
| group の root incidence | [`calcGradient_d.cu:754`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:754) で root を先頭に追加。漏れていません。 |
| `α`・行列・CSR・float 化 | §4.1 の指定式どおり。`ilp` は実行時と同じ CSR 位置で、double 演算後に係数を float 化しています。指摘1は設計側の問題です。 |
| `lsqPre_pinv` の切り出し | 打ち切り分岐と計算式を維持しています。非周期係数を意図的に変える差分はありません。 |
| `k/ω` の合併 | 周期半割面を除外し、合併体積で除算後に gather。[`periodicNode_d.cu:226`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:226) は和と broadcast の両方を実行します。 |
| 適用条件・呼び出し順 | LSQ と SST は共通条件関数を使用。`ransGradient → gather → ransBlendF1 → ransTransport` が成立しています。 |

目的の同定と切り出しも妥当です。既存の cell 周期保存修正、node 境界 DOF-only 化とは対象が異なり、今回の欠陥が解決済みという根拠はありません。陰解法の行縮約とスカラー全面 LSQ 化を別計画に残す判断も支持します。

**推奨は、合併 LSQ＋SST gather の方針を維持し、① incidence 幾何を一貫使用する設計へ修正、②試験の合否を固定、③F1 初期化と再現物を整備、④G0/G1/G2・R3 を通してから R1/R2 に進む、の順です。** ファイルは変更しておらず、提案は **plan 未反映**です。CFD は実行していないため、収束・準定常の合格判定は行っていません。

指摘数: Critical 0 / Major 2 / Minor 2
