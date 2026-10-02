# codex 諮問 (diagnose): twophase-dplur

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-dplur.md`](../../notes/reviews/briefs/2026-10-02-twophase-dplur.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `1f5a6953` (feature/species-transport)
- **codex**: effort `high`, 6.2 min, rc=0
- **結論**: **移流非対角だけを使う既存式の scalar-DPLUR に蒸気・液・Q の増分生成を接続し、単体条件を満たしてから、二相拡散 ON 同士の点対角対 DPLUR を各 2000 step 比較する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（レビュー対象 HEAD: `1f5a6953`）

| 論点 | 判定 | 根拠と対案 |
|---|---|---|
| 蒸気・液・Q の scalar-DPLUR 化 | **採用。ただし診断用の変更として** | 点対角経路は `implicitRelax` を使わず、`condTwoPhaseRelax` のみを掛けている。[twoPhaseDiffusion_d.cuh:159](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:159)。既存 scalar-DPLUR と同じ sweep 式を使う。 |
| 拡散非対角も同時に追加する | **初回は却下** | 既存の化学種 DPLUR の非対角は移流だけ。[speciesTransport_d.cu:408](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:408)。最初は移流非対角だけに限定し、分子・補正流束・乱流拡散は現在の残差と対角に残す。 |
| M 行列なら sweep 中も非負・単調 | **却下 — Major、最も重い穴** | 解く対象は符号が混在する全残差を右辺とした「増分」。M 行列だけでは更新後状態の非負性を保証しない。sweep 中のクリップは加えず、最後の共通 θ・安全側丸め・Q 非負化を維持する。 |
| 全残差 `R_v = R_w − R_g`、制限後に水を戻す | **採用** | 現行の残差変換と commit を維持する。[twoPhaseDiffusion_d.cuh:159](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:159)、[同:194](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:194)。輸送・ソースを別々に前処理して加え直さない。 |
| 点対角の比較経路、既定値 | **比較経路を残し、既定は当面点対角** | 二相拡散自体も既定 OFF のまま。2000 step の診断改善を正式受入・既定変更の根拠にはしない。 |

結論: **移流非対角だけを使う既存式の scalar-DPLUR に蒸気・液・Q の増分生成を接続し、単体条件を満たしてから、二相拡散 ON 同士の点対角対 DPLUR を各 2000 step 比較する。**

第 1 仮説: 点対角への置換による増分の緩和・空間伝播の違いが、過大更新と制限・再正規化の反復を増幅している。確度: **中**。

  根拠: 現行二相更新は各残差を対角で一度割る。一方、化学種は増分をゼロ初期化し、`implicitRelax` と `nStepInner` に従って近傍補正を伝播する。[speciesTransport_d.cu:974](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:974)。ブリーフの緩和感度はこの仮説を支持するが、**「エネルギーと組成の不整合が真因」という因果までは測っていない**。

  反証条件: 下記の実装同値性を確認した DPLUR 側でも事前基準を満たさなければ、**「この最小変更だけで 2000 step 内に所定の改善が得られる」仮説は棄却**する。DPLUR 全般の無効や、緩和差が無関係だとは結論しない。

第 2・第 3 仮説: 第 2 仮説は、輸送・相変化・温度応答の結合が残ること。確度は中以下、未確認。乾燥停止セルが Q 残差を**直接**支配する説は再提出しない。plan §5.1 #1b-r2 の記録では、200/200 更新で棄却されている。

判別 A/B:

**実装を次の形に固定する。**

- 右辺は `R_v = R_w − R_g`、`R_g`、各 `R_Q`。現在の拡散・ソースを全て含める。
- 対角は現在の `V/Δτ + transport_diag + V·src_jac` を維持する。蒸気の対角として既に計算された `transport_diag_Yw` を使う。
- 非対角は、既存の流入質量流束／組立時の隣接密度だけ。
- 増分は毎回ゼロから開始し、既存どおり  
  `δ⁽ᵏ⁺¹⁾ = ω D⁻¹(R + N_adv δ⁽ᵏ⁾)`  
  を `nStepInner` 回適用する。
- 最終増分に `condTwoPhaseRelax` を一度だけ掛け、その後に現行の θ・commit・再正規化を適用する。入口ピン行の増分ゼロ、周期の近傍集約と増分ミラーも既存経路を継承する。

ここでの ω は `implicitRelax`。**既存式には `(1−ω)δ⁽ᵏ⁾` がない**。[speciesTransport_d.cu:425](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:425)。通常の weighted Jacobi に置き換えると別の更新になる。非対角なし、`R/D=1`、ω=0.7 の代数確認では、既存式の 5 sweep 後は 0.7、weighted Jacobi は 0.99757 だった。

**単体の合格条件を先に固定する。**

1. `R=0`・初期増分ゼロなら、全 sweep・commit の増分と補正量が厳密にゼロ。`nStepInner=1`・`implicitRelax=1` では、同じ入力の点対角経路と一致する。これを非零残差でも確認する。
2. T1・T2 は §14.2 と §15 の修正後条件を維持する。特に独立残差による判定、3 セルの 1/1000 更新で全対象総量の相対変化 ≤1e−6、非負、上限到達ゼロを免除しない。
3. T3 の実ソース 1D は CFL 5、`DT_MAX` 1 K／0.01 K、上限 20000 外側反復。全対象成分について  
   `max|r_q| ≤ max(1e−7 r0_q, 6ε₃₂ max A_q)`、  
   非有限・負値ゼロ、実反復数の末尾 10% で Qcut・vround ゼロ、最終 θ・θ_src が全セル 1。0.01 K 条件では途中の θ<1 発生も必須。
4. 点対角の受入済み 1D 解との比較は、出口 g の相対差 ≤1e−4、出口 T の差 ≤0.01 K。**前処理を変えるので、旧 T3 の反復数 ±10% 条件は適用しない**。反復数は記録し、収束判定は独立残差で行う。

**case/16 の A/B は一組だけ。**

- 共通初期場: `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5`。
- 両側とも二相拡散 ON、`condTwoPhaseRelax=1`。
- A: 点対角。B: 上記 DPLUR。変更する設定は新設する solver 選択キーだけ。
- `implicitRelax` は流れ側の現行値を両側で維持する。「緩和 1」は `condTwoPhaseRelax` を指し、流れの緩和を同時変更しない。
- 各 2000 外側更新。判定窓は最後の 200 更新。履歴上の別ビルドの run を A の代用にしない。

事前基準は次のとおり。

| 指標 | B の改善条件 |
|---|---|
| Q2・Q1・Q0 の独立残差 | 各成分で、自身の初期値と A の末尾値の**両方の 0.1 倍以下**。単一最終値の偶然を避け、末尾 200 更新の最大値で比較する |
| θ=0 頻度 | 末尾窓のセル・更新数／200 が A の 0.1 倍以下 |
| Qcut | 同じ末尾窓・同じ集計定義で A の 0.1 倍以下。A がゼロなら B もゼロ |
| 再正規化 | 末尾窓の `max|f−1|` と各成分の `max_n C_q,n` を別々に比較し、各々 A の 0.1 倍以下。既に κ 以下の項目は κ 以下を維持 |
| 有限性・非負性 | 非有限ゼロ。commit 後、再正規化・bounds 前にも蒸気・液・各 Q が非負 |

κ は既存の **4.7683716e−7**。ブリーフの初期 Q 残差が再現すれば、初期値から決まる上限は Q2 **5.63e−3**、Q1 **2.76e5**、Q0 **2.25e13**。実際の判定では同じ初期場を現在の監査で評価した値を使う。

→ **B が全基準を満たすなら、最小 DPLUR 化による診断改善を支持する。満たさなければ「この変更だけで十分」を棄却する。** いずれの場合も、2000 step の結果から定常場の一致や正式収束は主張しない。

やらない方がよいこと:

- **M 行列を理由に θ を撤去すること。** 例えば M=`[[2,−1],[−1,2]]` は逆行列が非負だが、右辺 `[-1,1]` の補正は符号混在になる。既存式・ω=0.7 の 5 sweep を計算すると、初期状態 `[0,1]` は無制限なら `[-0.26062094,1.26062094]` になる。今回の全残差では非負右辺という前提がない。
- sweep 中に状態を更新して EOS・ソースを再評価したり、増分を途中でクリップしたりすること。凍結残差に対する線形反復という判別条件が変わる。
- 拡散非対角、θ 分割、ソース Jacobian の変更を一緒に入れること。
- DPLUR 化後の停止を、残差ゼロと同一視すること。`R=0 ⇒ δ=0` は保てるが、θ=0 や格納丸めによる停止では逆が成り立たない。

呼び出し側の前提への異議:

- **同じ Δτ・ω・sweep 数は、熱力学的整合の保証ではない。** 流れのブロックとスカラーでは前処理行列が異なり、凝縮ではさらに種間・温度との結合がある。今回確認できるのは更新方式を揃える効果まで。
- `scalarCflMax` が有効なら、スカラーの実効 Δτ は流れと異なる。[speciesTransport_d.cu:1228](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1228)。この値を確認せず「同一 Δτ」と記録しない。
- 再正規化は区間累積／現在総量を使わず、既存の更新ごとの `C_q,n` を使う。正式受入には従来どおり `check_convergence`、独立監査、再正規化ゲート、制限解除が必要。派生量の定常性は別途 `check_quasisteady` で判定する。

不足情報: 禁止指定に従い実測ファイルは読んでいないため、run 数値はブリーフ・plan の報告値として扱った。実効 `implicitRelax`、`scalarDtScale`、現在ビルドでの共通初期残差は未確認。ファイル変更・forge 起動は行っていない。**plan 未反映**。呼び出し側が `plans/active/condensation-two-phase-transport.md` の §4・§5.1・§6 に上記の実装範囲と事前基準を反映する。
