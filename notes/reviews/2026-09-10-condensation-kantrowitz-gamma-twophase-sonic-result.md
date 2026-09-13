# codex レビュー: condensation-kantrowitz-gamma-twophase-sonic (result)

- **plan**: [`plans/active/condensation-kantrowitz-gamma-twophase-sonic.md`](../../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-10
- **commit**: `5db69d93` (feature/sern-design)
- **codex**: effort `high`, 6.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 実装 diff (condensationEOS_d.cuh cond_twophase_sonic, dependentVariables_d.cu TP 分岐, condensationSource_d.cu の gamma_kw, condSonicResolve.hpp, check_convergence.py) が plan §4 と一致するか。§9 の検証主張 (回帰がノイズ床以内, E1 onset −0.45 mm, E2 が場に効かない理由=超音速で SLAU 全風上, PASS/STEADY) が run の実測 (case/16 run_0331–0342, compare_condfix*.txt, CONVERGENCE_VERDICT.txt) で裏付けられているか。accepted にしてよいか
- **extra**: `case/16.nozzle_wys/compare_condfix.txt`, `case/16.nozzle_wys/compare_condfix_cell.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

主要な音速・Kantrowitz 実装は §4 に沿っており、node 主要5 run の `PASS` と onset 約 −0.45 mm は再現できました。
しかし、cell の必須収束条件が未達で、回帰の合格根拠・SLAU の説明・判定ツールにも問題があります。現状では `accepted` に移せません。

指定の `git diff main...HEAD` を取得しましたが、対象実装は未コミットでした。そのため `git diff HEAD` の作業ツリー差分もレビューしました。ビルド・単体試験は再実行していません。保存場の全 `VALUE` の NaN/Inf 検査、収束判定、実カーネル照合は再実行しました。`--series` は描画初期化が read-only 環境で失敗するため、描画部分を除いて同じ計算・判定コードをメモリ上で実行しています。

1. **Major — cell の必須条件未達を「検証完了」にしている**

   [plan:204](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:204) は E4 を含めて `PASS` と `STEADY` の両方を要求しています。しかし再実行結果は次のとおりです。

   | run（`case/16.nozzle_wys/` 配下） | `check_convergence.py` | `compare_condfix.py --series` |
   |---|---|---|
   | `run_0331_condfix_head_ref`～`run_0335_condfix_new` | `PASS (converged)` | `STEADY` |
   | `run_0336_condfix_cell_head` | `NOT CONVERGED (stalled/plateau)` | `STEADY` |
   | `run_0337_condfix_cell_legacy` | `NOT CONVERGED (stalled/plateau)` | `STEADY` |
   | `run_0341_condfix_cell_new` | `NOT CONVERGED (stalled/plateau)` | `STEADY` |

   E4 の低下量は `rms_ro` 0.2 桁、`rms_roUy` 0.2 桁、`rms_roe` 0.2 桁、`rms_roOmega` 1.0 桁です。[保存 VERDICT](/home/sano/work/forge/case/16.nozzle_wys/run_0341_condfix_cell_new/CONVERGENCE_VERDICT.txt:2) も同じです。「既知の床」は規定の合格条件の代わりになりません。

   **対案:** §5.1 #5 と §8 の完了チェックを戻し、cell は「未収束の準定常比較」と記録する。restart 初期状態と残差床を切り分け、規定の収束条件を満たす検証を追加してください。単なる step 延長で解決するとは断定できません。

2. **Major — 「全 VALUE が反復ノイズ床以内」という回帰判定を裏付けられない**

   [plan:200](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:200) の基準と、[run 一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:336) の報告対象が一致していません。HDF5 から同じ `max|Δ|/max|ref|` を再計算すると、

   - node R1、`run_0331` 対 `run_0338`、24,000 step：`roQ0_0=1.905e-5`。既に「全場 ≤1.6e-5」を超えます。診断量まで含めれば `limiter_P=0.844`、`res_roK=1.164` です。
   - cell、`run_0336` を参照とした24,000 stepの `rog_0`：同一バイナリ反復 `run_0342` は **1.561e-5**、旧キー `run_0337` は **5.436e-5**。反復差の **3.48 倍**で、「≤1.5 倍」に収まりません。
   - node の反復 `run_0338` は24,000 stepで `NOT CONVERGED`。48,000 stepの回帰に対するノイズ評価としても条件が揃っていません。

   これは直ちにソルバの回帰を証明するものではありませんが、現在の合格宣言は成立しません。

   **対案:** 同じ終了条件で旧バイナリを複数回反復し、保存量・原始量・診断量それぞれに変数別の許容差を定める。凝縮モーメントを含む全対象列の結果を残し、§5.1 #4 を再判定してください。

3. **Major — E2 の「超音速だから SLAU が全風上で音速に依存しない」は誤り**

   SLAU は速度の大きさではなく、**面法線速度**から `M_p/M_m` を作ります。[法線速度:377](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:377)、[Mach 分岐:415](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:415)。法線方向が亜音速なら、`beta_p/beta_m` を通じて圧力流束に音速依存が残ります。[圧力流束:455](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:455)

   `case/16.nozzle_wys/run_0334_condfix_sonic/res_48000.h5` と入力メッシュから二次再構成すると、両側 `g>10⁻³` の内部面 **16,726 面中 8,304 面**は両側とも法線方向で亜音速でした。同じ場で音速だけ旧式に戻す後処理でも、これらの面の `p_tilde` は最大約 **0.024 Pa**変わります。

   また、参照との差は `g_0=1.434e-5`、`Uy=1.100e-5` で、README の「ρ,P,T,U,g が ≤4e-6」も誤りです。

   **対案:** 結論を「今回の検証条件では場の差が小さい」に修正する。原因説明には面法線 Mach と流束各項の感度を使い、「全風上」「影響0」という断定を削除してください。

4. **Major — 検証ツールが異常・評価不能を合格にする**

   [compare_condfix.py:99](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:99) は NaN を含む量を判定から除外します。実際の判定コードに `[1,2,3,NaN]`、全 NaN、全 Inf の系列を渡すと、いずれも **`VERDICT(series): STEADY`、`allok=True`** になりました。

   同様に [diff_res.py:7](/home/sano/work/forge/case/16.nozzle_wys/diff_res.py:7) は共通キーだけを比較し、形状不一致を失敗扱いにしません。[verify_sonic.py:52](/home/sano/work/forge/case/16.nozzle_wys/verify_sonic.py:52) は `gamma` の誤差を終了コードに反映しません。

   今回の保存場には NaN/Inf を検出していませんが、§6 の合否を保証するツールとしては不十分です。

   **対案:** 必須量の欠落・非有限値・形状不一致は非ゼロ終了にする。onset 未発生は明示的な評価不能状態として扱い、`sonic` と `gamma` の双方を合否条件に加えてください。

5. **Major — 現在仕様の既定適用範囲が実装と矛盾する**

   [methods/condensation.md:592](/home/sano/work/forge/methods/condensation.md:592) は `condEquilibrium 0/1` で自動 ON と記載し、境界条件の制限も欠いています。[同:169](/home/sano/work/forge/methods/condensation.md:169) は平衡拘束形の NS ブロックも二相 frozen に変更済みと説明しています。

   実装は [condSonicResolve.hpp:25](/home/sano/work/forge/solver_density_cuda/input/condSonicResolve.hpp:25) で `condEquilibrium!=0` を自動 OFF にし、未検証境界も除外します。`procedures/recommended-settings.md` は実装と整合していますが、plan §7 にも「非平衡/緩和形」が残っています。

   **対案:** `methods` と plan §7 を resolver の条件に統一する。CPG は明示 `1` でも変更されないことも、関連説明全体で揃えてください。

6. **Minor — 新既定の `h0` 誤差を参照値と取り違えている**

   [plan:255](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:255) と [README:341](/home/sano/work/forge/case/16.nozzle_wys/README.md:341) は新既定の中心線最大偏差を `0.267 kJ/kg、参照と同一` としています。`VALUE/h0` と正式比較コードでの再計算は、参照 `run_0331` が約 **0.2669**、新既定 `run_0335` が約 **0.2720 kJ/kg**でした。両者とも `VERDICT(series): STEADY` です。

   **対案:** 数値を訂正し、「同水準」の保存性許容差を明記する。約0.005 kJ/kgの増加だけから重大な非保存とは判断しませんが、「同一・悪化なし」の根拠にはできません。

7. **Minor — 前回採用した carrier 核生成率の限定が文書に反映されていない**

   plan §4.1 は上下限関係を「障壁・前因子を固定し、熱除去だけを追加するモデル間」に限定しています。しかし [methods/condensation.md:94](/home/sano/work/forge/methods/condensation.md:94) と plan §10 は依然として **carrier の真値を挟む上下限**として記述しています。§6.1 の「全件採用」と整合しません。

   **対案:** 両箇所を §4.1 の限定付き説明に統一してください。

推奨は、**`in_progress` を維持して検証をやり直すこと**です。まず判定ツールを修正し、cell の収束条件と回帰許容差の根拠を整え、その後に数値・因果説明・仕様文書を同期してください。`cond_twophase_sonic` の式、TP 分岐への接続、Kantrowitz の選択、CPG 分岐の維持については、今回確認した差分に重大な実装誤りは見つかっていません。

read-only 指示に従い、ファイル変更なし・plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 2
