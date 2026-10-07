# 凝縮 dual-time の液滴モーメントが同じバイナリでも run ごとに 2 通りの状態に分かれる

## メタ

- **area**: `condensation / time_integration`
- **status**: `draft`
- **related_docs**:
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) (凝縮・受動種・dual-time の設定、`condLimiterMode`・`passiveScalarScheme`・`implicitRelax`)
- **related_plans**:
  - 起点: [`architecture-solver-host-memory.md`](architecture-solver-host-memory.md) §6.2–§6.4 (ホストメモリ削減の回帰で、この構成だけ変更前後を判別できなかった)
- **created**: `2026-10-07`
- **owner**: `CFD Dev (Claude セッション: SERN 3D)`

## 1. 目的

凝縮つき dual-time (case/44 の汚染空気風洞ノズル、軸対称 node、`condensation 1`・`passiveScalarScheme 1`・BDF2) では、**同じ入力・同じバイナリでも**、
少数の節点で液滴モーメント (Q0・Q1・Q2) が run ごとに 2 通りの状態に分かれる。回帰試験の物差しにならず、結果の再現性そのものにも関わる。
何が分岐を作っているかを特定し、再現する (または分岐が物理・数値的に許容できる範囲だと示す) 状態にする。

## 2. スコープ

- **やる**: 分岐の位置・時点・引き金の特定 (どの補正・制限が入るか)、設定 (内反復・擬似 CFL・更新クランプ・緩和) への依存の測定、対策の要否の判断。
- **やらない**: ホストメモリ削減の受入れ (起点 plan で範囲限定済み)。凝縮モデル自体の物理の変更 (要るなら別 plan)。

## 3. 関連 docs と前提

- 観測の原本: `case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/notes.txt`・`RESULT.txt`・`result/describe.txt`、
  入力 `case/66.hostmem_regression/inputs/c44dual_ckpt100/` (元 `case/44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float`、100 step、`outStepInterval 100`)。
- 入力の設定: dual-time BDF2・`nSubIterDualTime 20`・`cfl_pseudo 12`・`nStepInner 5`・block-DPLUR・`implicitRelax` 未指定 (既定 1.0 = 緩和なし)・
  `condLimiterMode 1` (更新クランプ θ_u)・`condFloat 1`・`passiveFct 1`・`speciesImplicitCoupling 1`・Euler (`visc 0`・`thermCond 0`)・`output.level 2`。
- 関連する過去の知見 (メモリ記録、要再確認): case/44 のモーメント時間次数は float32 のノイズ床で決まっていた (核生成の温度感度が強い)。
  凝縮の θ 制限は `dt_local` 比例で解が擬似 CFL に依存していた (`condLimiterMode` の導入理由)。

## 4. 設計方針

**未定** (調査の後に書く。§4・§6 を書く時点で上位に諮る)。

### 4.1 観測事実 (2026-10-07)

- `condClampCorrQ_0` = $\max_k |\Delta Q_k| / \max(|Q_{k,\mathrm{before}}|, 10^{-30})$ (`condensationRealizability_d.cuh:313`、補正の相対量)。
  一部の run で節点 19954 (x 3.804、y 0.792) が約 1.46e23 になる。補正後の `roQ2_0` = 1.46e-7 が 1.46e23 × 1e-30 と 7 桁一致 = 補正前の Q2 がほぼ 0 (≤ 1e-30) で、
  そこへ 1e-7 程度の値を入れた状態と整合 (補正前の値は出力に無い)。
- その状態の run では、同じ節点で `roQ1_0` が他の run の約 2.1 倍、`roQ2_0` が約 4.8 倍、`roQ0_0` が −10 %。一方 ρ・T・P は 5 桁以上一致、g は 0.5 % 以内。
  → 少数節点の**液滴の大きさの分布だけ**が 2 通りの状態に分かれる。別の節点 (18784、23680 など) でも同型の外れ (相対補正 1e10〜1e20) が出る。
- 出現頻度 (`res_100.h5` の時点): 同じ入力の新しい標本で変更前バイナリ 2/6・変更後 2/6、既存の run を含めると変更前 3/12・変更後 6/12 (構成に `FORGE_PIN_DIAG=1` の変種を含む)。
  6 本ずつでは頻度差を絞れない。

### 4.2 仮説 (未検証)

- **H1 (ユーザの見立て、2026-10-07)**: 緩和・切り詰めの類 (更新クランプ θ_u、実現性の補正、FCT の制限、あるいは緩和係数) が、内反復の途中の値に対して入るかどうかで節点ごとに分岐し、
  内反復が収束しきらないまま物理 step が進むので、その経路が解に残る。GPU の加算順 (atomicAdd) の揺れが引き金。
- H2: 核生成の温度感度による増幅 (float32 のノイズ床)。H1 と両立しうる。
- 判別の候補 (未登録): 内反復の回数・擬似 CFL・`condLimiterMode` を 1 つずつ変えたときの分岐の出現頻度と、内反復の残差の落ち方。分岐の時点の特定 (毎 step 出力)。

## 5. 実装ステップ

調査の後に書く。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 分岐の時点と引き金の特定 | 同じ入力で 1 step ごとに出力する run を数本 (変更後バイナリ)、節点 19954 等で Q0〜Q2・g・`condClampCorrQ_0`・θ_u・内反復の残差を時系列で。どの step のどの補正で状態が分かれるか | O |
| 2 | 設定依存の判別 A/B の設計 | §4.2 の H1/H2 を判別する A/B を設計 (変える因子・本数・合格条件を事前に)。上位に諮る | F |
| 3 | case/44 の担当との共有 | この現象は case/44 の dual-time 全般に関わる。case/44 の README・関連 plan にポインタを置く | O |

## 6. 検証

調査の後に書く (§5.1 #2)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

case/44 系の凝縮つき dual-time の結果の再現性。回帰試験の設計 (この構成は再実行で値が分かれる)。

## 8. 完了条件

分岐の引き金を特定し、(a) 対策で再現するようになる、または (b) 分岐が用途に照らして許容できる範囲だと根拠つきで示す。

## 未確定事項

- この分岐が定常解・時間平均に効くか (100 step の値しか見ていない)。

## 9. 変更ログ

- `2026-10-07` — 起票。ホストメモリ削減の回帰 (起点 plan §6.3) で観測。ユーザの見立て (緩和・切り詰めの類が効いていそう) を H1 として記録。
