# 諮問ブリーフ: 凝縮域の二相拡散の再設計 — codex plan 段 NO-GO の採否と改訂 §4.2–§6 (2026-09-27)

作業ツリー `/home/sano/work/forge-species` (ブランチ `feature/species-transport`, HEAD 96f75d5a)。plan: `plans/active/condensation-two-phase-transport.md`。
AGENTS.md エスカレーション 1 (§4・§6 を新規に書く)・5 (codex Critical/Major の採否)。

## 問い (1 つ)

plan 段レビュー (`notes/reviews/2026-09-27-condensation-two-phase-transport-plan.md`, NO-GO C1/M6/m1) を全件採用して §4.2–§6 を改訂した。
**改訂版で #4 (カーネル実装) に進んでよいか。穴があれば最も重い 1 つに絞って直し方を示し、§6 の数値 (事前固定の合否) で直すべきものを挙げてほしい。**
採否 (全件採用) 自体に異論があれば指摘してほしい。

## 観測事実

- レビュー指摘 M3 はコードで確認: 既定の表引き経路 `gasProperties_d.cu:25-30` `gas_transport_cell_rY` は `roY` をそのまま渡すので、旧 §4.1 が対象にした `gas_transport_cell_Y` だけの修正では効かない。
- C1 は式で確認 (混合物基準 `−ρD_k∇Y_k` に気相組成の補正をかけると、気相一様・g 非一様で `ρ z_k (D_k − Σ z_j D_j)∇g` が残る)。
- run_0482/0483 は本体ワークツリー `/home/sano/work/forge/case/16.nozzle_wys/` にあり、species ワークツリーには無い (レビュー時に再検証できなかった理由)。
- ユーザ決定 (2026-09-27): 懸濁効果は無視、μ・λ は気相組成で、拡散は「蒸気の勾配で分子拡散・液は分子拡散なし・乱流は同じ Sc_t」、液 Sc は安定化用の選択肢として残したい。PC 負荷のため重い計算はローカルで回さない (AWS)。

## 改訂の要点 (本文は plan §4.2–§6)

- 分子拡散は気相基準: `j_k⁰ = −ρ_g D_k ∇z_k`、気相内補正 `j_k = j_k⁰ − z_k Σ j_j⁰`。液の分子拡散 0。有限 Sc_l は #8 に分離 (選択肢は残す)。
- 乱流は全輸送量 (Y_k, Y_w, g, Q0–Q2 の質量当たり) に共通の `μt/Sc_t`。
- エネルギー `Σ_{k≠水} h_k J_k + h_v J_w − L J_l`。
- 面流束を新カーネルで 1 回作り、化学種・液/Q・エネルギー残差と FCT 履歴で共有。呼び出しは液残差のゼロ初期化の後、境界ピン・周期集約の前。凝縮 TP では既存の水・気相種の分子拡散を置き換える。
- 非負: 蒸気 ρv と液 ρg を変数にしてそれぞれの対角で増分、`Δ(ρY_w) = Δρv + Δρg`。二点の正係数形。再正規化の係数を ρg と Q にも掛ける。
- 監視: 理由別 (上限違反・floor・増分制限・射影・消滅・再正規化) の区間値と累積。
- 対象外: CPG carrier (空気凝縮)、cell は未検証と記録。

## 仮説 / 当方の懸念 (棄却してよい)

- 二点の正係数形 (面法線差分) にすると、非直交補正を捨てることになり、現行の化学種拡散 (勾配使用?) と離散が変わる。g=0 の run でビット一致を保つには、凝縮 TP 以外で旧経路を残す必要がある。
- 蒸気/液変数での陰的対角は、化学種の他成分 (N2 等) の対角との整合 (ΣρY = ρ) を壊さないか。再正規化で吸収する設計でよいか。
- Q0–Q2 の乱流輸送の変数 (質量当たりか体積当たりか) と実現可能性 `Q2² ≤ Q1Q3` の保存: 同じ係数の線形拡散なら凸結合で保たれるはずだが、陰的更新で保てるか。
- §6 の保存許容 1e-6 (float32、1000 更新) は厳しすぎ/緩すぎないか。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `trans.inp`・`thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `plans/active/condensation-two-phase-transport.md` 全体
- `notes/reviews/2026-09-27-condensation-two-phase-transport-plan.md`
- `notes/reviews/2026-09-27-condensation-diffusion-error-diagnose.md`
- `sed -n '200,420p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`、`sed -n '1240,1300p;1690,1760p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`
- `sed -n '180,260p' solver_density_cuda/cuda_forge/passiveFct_d.cuh`
- `sed -n '80,200p' solver_density_cuda/cuda_forge/condensationRealizability_d.cuh`
- `sed -n '360,420p' solver_density_cuda/cuda_forge/condensationTransport_d.cu`
- `sed -n '340,400p' solver_density_cuda/cuda_forge/condensationEOS_d.cuh`
- `sed -n '1800,1850p;2070,2110p' solver_density_cuda/main.cpp`
