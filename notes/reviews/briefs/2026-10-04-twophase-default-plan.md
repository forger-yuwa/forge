# 諮問ブリーフ: 二相拡散の既定化 plan の §4 (設計方針)・§6 (検証計画)

エスカレーション条件 1 (後継 plan の §4・§6 を新規に書く)。ユーザ指定で diagnostician (Fable) に諮る。

## ユーザ決定と注意 (2026-10-03〜04)
- **決定**: 「物理的には二相拡散は ON にしないといけない。既定にしていってほしい」。理由の整理: 液が拡散するなら、その液が持つエンタルピー (h_l = h_v − L) を運ぶ −L·J_l が要る。既定経路 (OFF) は総水分 Y_w を乱流・分子拡散し (エネルギーは h_w J_w = 全部蒸気扱い)、液 ρg・モーメントは拡散しない → 液が混ざらず、潜熱の乱流輸送が欠ける (case/16 で欠けている潜熱流束/熱流束 p50 0.77・p95 4.8; `case/16.nozzle_wys/README.md` の `analyze_liquid_diffusion_error.py` 行)。
- **注意 (ユーザ)**: 「昔、relaxation 起因で dual-time の解が変わってしまったことがあったから気を付けて」。実例: `plans/accepted/species-passive-scalar-unification.md` §5.1 #26 (2026-09-17, case/44 `run_0357`–`0377`): dual-time の sub-iter 数依存の原因は `implicitRelax 0.7` 単独。緩和は不動点を変えないが遅いモードの収束を遅らせ、nSub 内で sub-iter が収束しきらず物理 step の解が nSub に依存する。残差ノルムには現れない。→ dual-time の時間次数・sub-iter 試験では緩和を使わない、と決めた。

## 現状 (plan `plans/active/condensation-two-phase-transport.md`、仕様 `methods/condensation.md` §7c)
- 実装: opt-in `condensation.condTwoPhaseDiffusion` (既定 0)。**定常擬似時間専用の初版**。起動時拒否: dual-time・陽解法・`speciesImplicitCoupling 2`・`passiveScalarScheme 0`・`condEquilibrium ≠ 0`・`condLimiterMode 0`・凝縮種 2 以上。CPG carrier (空気凝縮) と pure 凝縮は対象外 (液は拡散しない)。cell は未検証。
- 中身: 気相内の分子拡散 (z 基準、風上補正)、全相共通の乱流混合 μt/Sc_t、エネルギー Σh_kJ_k + h_vJ_w − L·J_l、非分割の全残差更新 R_v = R_w − R_g、θ 制限 (dg_max/dT_max、非負 θ は `condTwoPhaseNonnegLimit 0` で外せる)、再正規化係数の液・Q への適用、`condTwoPhaseRelax` (既定 1)、`condTwoPhaseSolver` (0 点対角 / 1 緩和整合 scalar-DPLUR)。
- 実装ファイル: `solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`、`condensationTransport_d.cu`、`speciesTransport_d.cu` (twophase_diffusion_d・監査・再正規化ゲート)、`renormGate_d.cuh`、`condensationRealizability_d.cuh`、`passiveKernels_d.cuh`、`main.cpp`、`input/solverConfig.{hpp,cpp}`。単体試験 `tests/unit/test_twophase_*`。
- 検証済み: case/16 (平面 2D node SST 定常、TP carrier) で ON (DPLUR・非負制限なし) が (A) 7 報告量 STEADY・(B) 負値/非有限 0・(C) RISING 0・(D) 補正 ≤κ を満たす (#4j、`run_0521` vs OFF `run_0520`)。ON−OFF: 出口 g 質量重み +3.2 %、壁圧 x42 +0.8 %、dev_pct 2.47→2.93 %、壁温 −2.4 K、onset 不変。**機構は未確定** (codex: OFF は ∇Y_w、ON の分子拡散は ∇z_v なので壁側の総水分低下は作用素差と整合、機構判別 #4k は未実施)。
- 回帰 (#5a): 既定経路 (OFF) は二相拡散の commit 前後でノイズ床内 (湿り凝縮 NS・乾き・Euler の 3 ケース PASS, `run_0540`–`0554`)。ON の軸対称・周期は入力が無く未検証。
- 関連の既知事項: dual-time の受動種は post-step FCT (`passiveFct`, 流束形 BDF2 履歴) で有界化 (species-passive-scalar-unification、accepted)。凝縮 θ リミッタの定常解の dt 依存は `condLimiterMode 1` で解決済み (condensation-source-limiter-steady、RK/dual-time は自動で 0 に降格・未検証)。残差の床はリミッタ既定由来で ON/OFF 共通 (limiter-inlet-column-oscillation、accepted)。

## 問い
1. 既定化の範囲と順序: どの制限から外すべきか (dual-time、陽解法、軸対称・周期、`speciesImplicitCoupling 2`、`passiveScalarScheme 0`、平衡凝縮、`condLimiterMode 0`、凝縮種 2 以上、CPG carrier、pure 凝縮)。全部を外してから既定にするのか、「使える構成では既定 ON・使えない構成では警告して OFF」の段階既定化か。後者の危険 (同じ物理が構成によって違う式になる) をどう扱うか。
2. dual-time への拡張の設計: 非分割更新・θ 制限・非負/実現可能性クランプ・再正規化・緩和 (`condTwoPhaseRelax`) を sub-iter の中でどう扱えば、**sub-iter が収束したとき BDF の解が nSub・緩和・θ に依存しない**ようにできるか (#26 の教訓)。既存の受動種 post-step FCT (物理 step 末尾の保存的 FCT) に二相の水 (蒸気・液・Q) をどう乗せるか。検証で何を固定・何を測るか (nSub 倍増差、3 水準の時間次数、緩和なしの対照、sub-iter 残差の低下桁)。
3. 機構判別 (#4k) は既定化の前提にすべきか。物理として正しい方向であることの根拠を何で示すか (解析解・文献ケース・SU2 等は無い。0 step 作用素 A/B と閉じた収支の提案あり)。
4. 既定変更の手順: 既存 run (config に `condTwoPhaseDiffusion` が無い凝縮 run) の結果が変わる。互換の扱い (旧挙動の明示キー、restart 時の扱い、記録)、検証ケースの範囲 (node 標準ケースのうち凝縮を含むもの: case/16・case/44・case/34 Arthur 等)、変化量の記録方法。
5. この plan の §6 として事前登録すべき合格条件の骨子。

読んでよいファイル: 上記すべて、`plans/active/condensation-two-phase-transport.md`、`plans/accepted/species-passive-scalar-unification.md` (§4.4・§5.1 #12・#26・#28)、`plans/accepted/condensation-source-limiter-steady.md`、`methods/condensation.md`、`methods/thermophysics.md`、`procedures/recommended-settings.md`、`solver_density_cuda/main.cpp` (`advanceImplicitDualTime`)、`solver_density_cuda/cuda_forge/` 配下。
