# 諮問ブリーフ: 二相更新の非負制限 θ を外し、負値は既存の実現可能性クランプで受ける (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD は git log -1)。plan `plans/active/condensation-two-phase-transport.md` §5.1 #4g (DPLUR 化の結果) と #4 (codex plan 段 M2 の判断: 非負は更新で保つ、bounds 吸収を非負保証と呼ばない)。AGENTS.md エスカレーション 6 (カーネル変更) と、過去の判断を覆すので 5 相当。**ユーザ提案・了承 (2026-10-02)**: 「時に負になっちゃうことは許容したら、その制限をしなくて済むんじゃないか。トラブルの元な気がする」→ 主セッションが下の案を提示し「ほい」。

## 観測事実

- θ = min(θ_しきい値 (dg_max・dT_max), θ_非負 (ρv/(−δρv), ρg/(−δρg))) を蒸気・液・Q に共通で掛ける (`twoPhaseDiffusion_d.cuh` `tp_vl_update`)。ρg = 0 で δρg < 0 なら θ = 0 でセル全体の更新が止まる。
- case/16 #4g A/B (run_0515 点対角 / run_0516 DPLUR、2000 step、末尾 200 更新): DPLUR は Q 残差を A の 1/2800 級、Qcut・再正規化を大幅改善したが、θ=0 頻度が 105 セル/更新 (A 359) で事前基準 (0.1 倍) 未達。run_0516 の θ=0 は 179 固有セル (61 セルが 90 % 以上の更新)、**理由は全て液の非負 (20977 件)、dg_max・dT_max は 0 件**、位置は x 17〜28 mm・y 2.7〜2.9 mm (凝縮開始手前の壁寄り)。
- #1b-r2 (run_0511、緩和 0.5、点対角): 乾燥停止セルの液残差は double でも負 (最小 −1.1e-34)、最も負の項は全行で移流; 停止セルは Q 残差を直接支配しない。
- 既存の実現可能性クランプ (`condensationRealizability_d.cuh`: 0 ≤ ρg ≤ ρY_w、モーメント射影、液滴消滅) と `[cond-corr]` の理由別計上、`[renorm-gate]` は実装済み。

## 案

- θ から**非負の項を外す** (θ = θ_しきい値 のみ)。dg_max・dT_max (硬いソース対策) は残す。
- 更新後に負になった蒸気・液・Q は既存の実現可能性クランプで受け、補正量を `[cond-corr]` で監視し受入条件に入れる (例: 末尾窓のクランプ補正量 / 総量 ≤ κ = 2n_sε₃₂)。
- EOS・ソースが負値を直接読まないこと (値を使う時点で max(·,0) かクランプ後の値) を確認・保証する。
- 新キー (例 `condTwoPhaseNonnegLimit`: 1 = 現行 (既定) / 0 = 外す) で A/B: case/16、DPLUR、同じ初期場 run_0482 res_48000、2000 step。

## 問い

1. この方針変更 (M2 の「非負は更新で保つ」を、定常擬似時間の途中では負を許しクランプで受ける、に変える) は妥当か。穴があれば最も重い 1 つ。特にクランプは保存を崩す (液を削る = 強制蒸発) ので、補正量の許容の定義と、固定点での補正 0 (またはκ以下) の要求の仕方。
2. 負値を読む箇所 (二相 EOS、核生成・成長ソース、モーメント実現可能性、物性の気相組成) の確認範囲。
3. A/B の事前基準 (θ=0 頻度は外す代わりに何を見るか: クランプ補正量・独立残差・非負の最終状態・再正規化)。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- plan、`notes/reviews/2026-09-27-condensation-two-phase-transport-plan.md` (M2)、設計メモ §14–§20、`cuda_forge/twoPhaseDiffusion_d.cuh`、`condensationRealizability_d.cuh`、`condensationTransport_d.cu` (grep `realiz`・`clamp`)、`condensationEOS_d.cuh` (grep `rog`・`fmax`)、`condensationSourceF_d.cuh` (grep 冒頭 100 行)
