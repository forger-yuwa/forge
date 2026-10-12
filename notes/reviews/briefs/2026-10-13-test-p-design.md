# 諮問ブリーフ: 試験 P (point の陰解法 + `lowMachPrecond 3`) の設計と事前登録 (2026-10-13)

- 依頼者: 主セッション (エスカレーション条件 1: §6 の検証計画を新規に書く、条件 4: 新しい run の前)
- plan: `plans/active/architecture-float-state-double-geometry.md` §5.1 #25 (ユーザの指示)・#26 (前処理の移植)、§6.57・§6.59・§6.60 と結果、§6.61 の結果
- ユーザ: 2026-10-12「低マッハすぎて収束性悪いと思ってる」「2 (低マッハの前処理の試験) はじゃあまずやってほしい」。2026-10-13「3 で試験を、あなたはするのかと思ってたけど。分担と状況把握はちゃんとしてくださいね」→ 試験 P はこのセッションの担当。`lowMachPrecond 2` (右辺 + 左辺) の試験は forge-25 が引き受けた (ユーザ決定、core-grid plan)。こちらは 2 を回さない。
- 前の諮問: `notes/reviews/2026-10-12-speedup-tests-diagnose.md` (P をそのまま 20,000 step は要再検証 → まず β = 1 の 1 step で境界の処理を確かめる (§6.57 → 結果 A、欠陥あり → §5.1 #26 で移植)。「FP64 同士なので前処理だけの比較」は却下 (P0 の `implicitSolvePrecision 0` は ST = float、P3 の行列・解法は flow_float = FP64 のビルドで double)。「停止値への残りが半分なら効く」は却下 (過渡の通過でも近づく)。加速は同じ収束・準定常の要件までの時間で判定、20,000 step は有限期間の応答試験に限る。速度は専有で)。

## 観測事実

- 移植 (`e224ae11`): 前処理のカーネルに node の行の処理 (境界の半割面の粘性の対角を除く、軸の ρv・壁 3 行・強いピンの等温壁 ρE を単位行) を入れた。R2 (§6.59): 結果 A (拘束の成分の max|d| が 0、入口・出口の更新の抑制が直った)。R1・R3 (§6.60): FP64 の 3 組は結果 A (前処理なしの経路の値は同等)、float の 3 組は判別不能 (main 統合は保留、#27)。前処理の経路 (`lowMachPrecond ≥ 2`) は β < 1 の多成分の `speciesPrecondDt` が未検証 (R2 は β = 1)。
- 移植後の FP64 のバイナリ: B の `~/forge-port/solver_density_cuda/build-post-fp64/forge` (sha256 75ee4675b0c90f41)。§6.58・§6.61 の run は移植前の 7e95e5bb (前処理なしの経路は SASS で同一)。
- §6.61 の結果: 登録した G0c の Q32 の起点・FP64・B0 (上限 50) は約 40.8 分 (占有) で停止し、T2 の停止の出力との比較の要件を満たした (`run_0804_sp_g0c`)。G0・B0 は止まるが θ_r が T2 から約 +1 %。G0c・FP64・B0 を登録の条件に限って採用した。
- forge-25 の結果 (core-grid plan §4.83、他セッション、こちらでは判定の出力を未読): Gref (forge-25 の格子) の `run_0608` res_60000 から point (cfl_pseudo 1) で 2 万 step、対照 `run_0688_pp_c` と `lowMachPrecond 3` の `run_0689_pp_p3`。前処理の資格 (島の刻み 6.04 倍・境界の凍結なし) は満たしたが、縮流部の ρ の残差は対照 217 → 49、=3 は 686 → 1296 (12〜21 倍悪い)。forge-25 は e224ae11 を cherry-pick した FP64 (A の `~/forge-fgeom7-fp64-lmp2`、sha d58ed787) を使った。B0 + `=1` (右辺だけ) は step 11 で発散 (§4.81)。
- 生産の point の設定 (`run_0713_pt_d1`): B0 からライン陰解法・方向別 dt・上限・ITJ のキーを外したもの (cfl 4・relax 0.7)。

## 設計の案

- 目的: point の陰解法に `lowMachPrecond 3` を足すと、G0c の Q32 の起点から停止規則の要件までが速くなるか。比べる相手は、同じ起点の point の対照 P0 と、いまの生産の候補 G0c・B0 (`run_0804`、40.8 分)。
- 腕 (B、GPU を占有して逐次、移植後の FP64 75ee4675 を両腕に、G0c の `_g0c/q32_init.h5`、設定は `run_0713_pt_d1` の point の設定を G0c に):
  - P0: point、`lowMachPrecond 0`
  - P3: point、`lowMachPrecond 3`、`precondEps` は既定
- 段 1 (有限期間の応答、各 20,000 step): 全残差列・θ_r・Q_w の推移と 1 step の時間を記録。P3 の縮流部の残差・全体の残差が P0 より明らかに悪い (例: 末尾 2,500 step の中央値で 2 倍以上) なら、速さの改善は無いとして段 2 に進まない。
- 段 2 (段 1 で P3 が劣らなければ): 両腕を §6.11 の停止規則 (上限 20 万 step) まで。停止の出力を T2 と M1〜M7 で比べ、T_reach を `run_0804` と比べる (§6.61 と同じ測り方)。
- 原因の帰属は分ける: P0 と P3 の差は「設定一式の効果」(前処理 + 解く精度 ST float 対 flow_float) と書く。前処理だけの効果を言うなら、P0 に `implicitSolvePrecision 1` の腕を足す (今回は入れない案)。

## 問い

1. この設計 (G0c・point・P0 対 P3・段 1 → 段 2) でユーザの問い (低マッハの前処理で速くなるか) に答えられるか。比べる相手 (P0 か `run_0804` の B0 か) の置き方。
2. 段 1 の打ち切りの基準を事前にどう書くか (残差の比・場所・どの列・何 step の窓)。forge-25 の結果 (Gref・cfl_pseudo 1) をこちらの判定にどう使うか、使わないか。
3. cfl (生産の point は 4、forge-25 は 1) と `precondEps` の選び方。cfl を振るなら何水準か。
4. 移植後のバイナリ 75ee4675 を両腕に使うこと (前処理なしの経路は FP64 の R1 で同等)、β < 1・多成分 (`speciesPrecondDt`) が未検証であることの扱い (case/45 は 2 成分)。
5. `implicitSolvePrecision` の違いを分ける腕を入れるべきか。

## 読んでよいファイル

- float plan §5.1 #21・#25・#26、§6.57・§6.59・§6.60・§6.61 と結果
- `notes/reviews/2026-10-12-speedup-tests-diagnose.md`
- `procedures/solver-settings.md` の lowMachPrecond・precondEps・speciesPrecondDt の節
- `solver_density_cuda/cuda_forge/timeIntegration_d.cu`・`lowMachPrecond_d.cuh`・`setDT_d.cu`
- `case/45.isobutane_m6_d155/sp_arm.sh`・`sp_judge.py`・`_band_ab/cold_pair/sp_judge.json`
