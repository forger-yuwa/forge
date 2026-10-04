# 諮問ブリーフ: G3 (閉じた収支) の設計と G2 の解釈

エスカレーション条件 1・6・7。plan: `plans/active/condensation-two-phase-default.md` §5.1 #4g3 (codex の要件)、#5 (G2 結果)、§6 G3。
前回諮問: `notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md` (D2 却下、作用素収支と更新写像収支に分ける)。
ユーザ判断 (2026-10-04): G3 を既定化前に実施する (「1」)。

## G2 の結果 (解釈の確認をお願いしたい)
`notes/investigations/2026-10-04-twophase-g2/evidence/G2_VERDICT.txt`。9 run 全て (A)〜(D) ok、7 量すべて ON−OFF が有意 (3s 超)、レシピ感度 abs(生産 − cfl2) は ON−OFF 差の 1e-4〜3e-2 倍、抽出照合一致。
私の解釈: 「生産レシピ (cfl_pseudo 4・implicitRelax 0.7・nStepInner 4) でも ON の準定常報告量は cfl2 と同じ値に落ち着き、既定化の効果 (ON−OFF 差) はレシピに依存しない」。異議があれば。

## G3 設計案
入力場: G2 の最終場 — ON = `case/16.nozzle_wys/run_0561_g2_onprod_r1/res_48000.h5`、OFF = `run_0567_g2_off_r1/res_48000.h5` (それぞれ自分の設定で評価)。

### G3-a 作用素の収支 (0 step、時間当たり)
- 成分 c ∈ {ρY_w, ρg, ρv = ρY_w − ρg, ρQ2, ρQ1, ρQ0}。単位は成分ごと (ρ 系は kg/s/m 奥行、Q 系はそれぞれの単位/s)。
- 制御体積 Ω (節点 CV の集合、事前固定): 全域、x 帯 [入口, 17 mm)・[17, 35)・[35, 70)・[70, 出口]、壁距離帯 wd < 0.4 mm と ≥ 0.4 mm (x 35〜70 mm 内)。
- 項 (すべて double で積算): 外向き正の境界流束 F_adv,∂Ω (移流、本番の移流カーネルの atomic 前の面流束を診断用の別配列にも書く)、F_diff,∂Ω (拡散: OFF は species_diffusion_d、ON は twophase_diffusion_d の atomic 前の面流束)、相変化ソース S⁺Ω・S⁻Ω (本番が実際に残差へ加える S·V を正負別に保存)、入口ピンの寄与 B_pin (ピン前後の差)、その他 R_other、最終残差 Σ_Ω R_final。
  恒等式: Σ_Ω R_final = −F_adv,∂Ω − F_diff,∂Ω + S⁺ + S⁻ + B_pin + R_other + E_assembly。蒸気は double で R_v = R_w − R_l 等。
- 記録点: wrapper 内部 (speciesTransport の移流後・拡散後、化学ソース、ピン前後、condensationTransport の移流後・二相拡散後、condensationSource 前後、passivePinResidual・2 度目の speciesPinResidual 前後、最終)。本番カーネルは書き換えず、診断有効時だけ atomic 前の値を別配列へも書く分岐を足す (既定経路は分岐のみ)、または D1 と同様に診断専用の写しで評価して本番との一致を検査する。
- **合否 (案)**: (1) 組立の整合: abs(E_assembly) ≤ 64·ε₃₂·Σ abs(各項) (float の atomic 加算の丸めの床)、全 Ω・全成分。(2) 物理の収支: 準定常の場で残る残差の寄与 abs(Σ_Ω R_final) が、液 (ρg) と総水分 (ρY_w) について、観測された液流束差 ΔF_l = abs(F_l,exit^ON − F_l,exit^OFF) (出口断面の移流 + 拡散) の 10 % 未満。Q 成分は対応する出口の Q 流束差の 10 % 未満。ρv は記録。

### G3-b 更新写像の収支 (1 更新、更新当たり)
- 同じ入力場から 1 外反復だけ更新し、成分ごとに Σ_Ω V(q_after − q_before) = Σ_Ω V·δq_update + Σ_a Σ_Ω V·C_q,a + E_update。
  δq_update = DPLUR・緩和・増分制限を通った実際の更新量。C_q,a = 各補正操作の符号付き格納差 (再正規化、commit の床、passive floor、上限クランプ、射影、液滴消滅、境界上書き; `twoPhaseHoldWater` の一時的な水更新と巻き戻しは二重計上しない)。符号付き積分と絶対量積分の両方。
- **合否 (案)**: (1) abs(E_update) ≤ 16·ε₃₂·Σ abs(V·Δq) (全 Ω・全成分)。(2) 補正の絶対量積分 Σ_a Σ abs(V·C_q,a) ≤ κ·Σ_Ω V·q (κ = 4.77e-7、既存 corr-gate と同じ基準) — ただし液滴消滅は物理的な除去として記録のみ。(3) ON と OFF の補正の内訳を比較して記録 (ON だけが液・Q に再正規化を掛ける)。

## 問い
1. G2 の解釈に異議は。
2. G3-a・G3-b の分け方・項・記録点・恒等式は前回の要件を満たすか。足りない項は (例: 境界半割面の移流は残る・拡散は skip、ghost、node の境界節点、周期は包絡外)。
3. 合否の尺度 (64ε₃₂ / 16ε₃₂、10 % ΔF_l、κ) は妥当か。単位の違う Q 成分の扱い、ΔF_l の定義 (出口断面の移流 + 拡散、ON−OFF) を具体化してほしい。
4. 実装方針: 本番カーネルに「診断時だけ別配列へも書く」分岐を足すのと、診断専用の写しで評価するのと、どちらがよいか (D1 は写し + 一致検査にした)。
5. 入力に G2 の最終場 (run_0561・run_0567) を使うのでよいか。

読んでよいファイル: 上記 plan と記録、`solver_density_cuda/main.cpp` (組立・更新・`twoPhaseHoldWater`)、`solver_density_cuda/cuda_forge/speciesTransport_d.cu`・`condensationTransport_d.cu`・`condensationSourceKernels_d.cuh`・`condensationRealizability_d.cuh`・`passiveKernels_d.cuh`、`notes/investigations/2026-10-04-twophase-g2/`、`notes/investigations/2026-10-03-twophase-onoff/wflux.py`。
