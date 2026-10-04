# 諮問ブリーフ: 二相拡散の蒸気・液・Q 更新を点対角から緩和整合 scalar-DPLUR へ (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD は git log -1)。plan `plans/active/condensation-two-phase-transport.md` §5.1 #4e (点対角前処理を選択)・#1b〜#1b-r2 (case/16 で受入未達、区切り済み)。設計メモ §14–§19。AGENTS.md エスカレーション 6 (cuda_forge の数値の振る舞いを変える編集の前)。**ユーザ指示 (2026-10-02)**: 「もうそれじゃん、すすめて」— 下の仮説で実装に進む。

## 観測事実と仮説

- forge の定常陰解法の 1 擬似時間ステップ (`main.cpp:2081` `implicitNonlinearUpdate`): (1) `assembleResidual` で流れ・乱流・化学種・液・Q の残差を同じ状態 Uⁿ から一度に計算 (凍結残差)、(2) 流れ 5 変数を block-DPLUR で `nStepInner` (case/16 は 5) 回 Jacobi sweep して commit、(3) 化学種を凍結残差のまま scalar-DPLUR (`speciesImplicitCoupling 1`; 同一 Δτ・`implicitRelax`・sweep 回数) で commit し ΣρY=ρ へ再正規化、(4) 液・Q を受動種 scalar-DPLUR (`passiveImplicitCoupling 1`)。
- `methods/thermophysics.md:102-108`: 化学種を流れと異なる緩和で進めると (e,Y) が各ステップで不整合になり T=T(e,Y) が跳ね、H2O (生成エンタルピー大) で発散を誘発する → 緩和整合 scalar-DPLUR を採用した経緯。scalar-DPLUR の 1 次風上移流の凍結 Jacobian は M 行列 (対角 = 流出質量流束、非対角 = 流入質量流束)、`nStepInner=1, ω=1, 非対角無視` で点陰的に一致。
- 二相拡散 (#4e) は蒸気・液・Q を**点対角**で更新 (`twoPhaseDiffusion_d.cuh` `tp_vl_update`、`condensationTransport_d.cu` `twophase_vl_update_d`; メモ §14.1) — 水だけ流れとの緩和整合を外している。
- case/16 (SST・凝縮、run_0482 res_48000 から): 緩和 1 で最初の 200 step に残差 2 桁上昇→平坦、θ<1 は全て非負制限 (蒸気 28874・液 43190、dg/dT 制限 0)、緩和 0.5 で 1/10〜1/60 に改善 (振動の型を支持) するが未収束。
- **仮説**: 点対角による水の緩和不整合 → (e,Y) の不整合・行き過ぎ → 非負制限 θ の多発 → 停滞/振動。

## 素案 (実装)

- 蒸気・液・Q の更新を、化学種と同じ緩和整合 scalar-DPLUR に載せる: 凍結全残差 (R_v = R_w − R_g、R_g、R_Q) を右辺、対角 = V/Δτ + 輸送対角 (移流流出 + 拡散) (+ ソース Jacobian)、非対角 = 隣接の流入質量流束 + 拡散係数、流れと同じ `nStepInner` 回の Jacobi sweep・同じ `implicitRelax`。sweep の後に θ (`vl_limit_commit`、安全側丸め) と commit、再正規化は現行どおり。`condTwoPhaseRelax` は sweep 後の増分に掛ける (既定 1)。
- 点対角は opt-in の比較用に残す (新キー例 `condTwoPhaseSolver: 0=点対角, 1=DPLUR`)。初版の既定をどちらにするかは結果を見てから。

## 問い

1. 素案に穴があれば最も重い 1 つ。特に (a) 二相拡散の拡散非対角 (気相内補正 −z_k Σj⁰ の非線形・種間結合、風上 z) を DPLUR の非対角にどう入れるか (入れずに移流 + 乱流拡散の非対角だけでよいか)、(b) DPLUR の sweep 中の中間状態での非負性 (M 行列なら Jacobi sweep は単調か)、(c) 凍結全残差 R_v = R_w − R_g と DPLUR の線形系の整合 (固定点が変わらないこと)。
2. 事前に固定する合格条件: 単体 (T1–T3 の再実行、1D 実ソースで DPLUR 版の収束と固定点一致)、case/16 の判別 A/B (同じ初期場 run_0482 res_48000、緩和 1 で 点対角 vs DPLUR、2000 step、#1b-r1 と同じ基準: Q 残差・θ=0 頻度・Qcut・再正規化が 1/10 以下 → 仮説支持)。数値を確定してほしい。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- plan、設計メモ §14–§19、`methods/thermophysics.md:95-115`、`solver_density_cuda/main.cpp:2081-2230`、`cuda_forge/twoPhaseDiffusion_d.cuh`、`condensationTransport_d.cu` (grep `twophase_vl_update`・`passive`)、`speciesTransport_d.cu` (grep `speciesImplicitDPLUR`・`passiveImplicit`)
