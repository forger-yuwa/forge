# node のスカラー勾配を LSQ に統一する (k/ω・化学種・受動種・凝縮モーメント)

## メタ

- **area**: `gradient`
- **status**: `draft`
- **related_docs**:
  - [`methods/gradient.md`](../../methods/gradient.md) (「境界寄与」node × 周期の継ぎ目、スカラー勾配)
  - [`methods/discretization.md`](../../methods/discretization.md) §7.3 / §7.3.1 (「LSQ は NS だけ、スカラーは GG」の記述を更新する)
- **related_plans**:
  - [`boundary-node-periodic-gradient-fix.md`](../accepted/boundary-node-periodic-gradient-fix.md) (前提。§4.1 の合併 LSQ 係数は「変数に依らない形」で作ってある。§5.1 #7 の継続課題 (1)(3)(4) を本 plan に移す)
  - [`boundary-node-rotational-periodic.md`](boundary-node-rotational-periodic.md) (回転周期。本 plan では扱わない)
- **created**: `2026-09-26`
- **owner**: `sano`

## 1. 目的

node の勾配は、NS の原始量 ($\rho, u, P, T$) だけが LSQ (`gradLSQ: 2`、事前計算係数) で、スカラー ($k,\omega$、化学種 $Y_s$、受動種 $\xi$、凝縮モーメント) は Green–Gauss (GG) である。
2026-09-26 のユーザ決定 (「そろえるべき」) により、node のスカラー勾配も同じ LSQ に揃える。

根拠 (前提 plan の測定):

- GG は float32 で部分和 $\Sigma\phi_f S_f$ が $\sim\phi/h$ の打ち消しになり、**壁半割面込みの CV では定数場が閉じない** (channel 実測、壁節点で約 $5.3\,\varepsilon|\phi|/h$、継ぎ目に依らない。前提 plan §5.1 #7(3))。LSQ は定数場で厳密に 0。
- NS の LSQ は継ぎ目の合併係数 (`lsqPre_mergePeriodic`) を持ち、係数は変数に依らないので、スカラーも**同じ係数を流用**できる (並進周期の継ぎ目の扱いが NS と自動的に揃う)。
- 同じ場に NS は LSQ、スカラーは GG と作用素が混在していると、境界の扱い (LSQ は内部隣接のみ、GG は owner 値の半割面) が変数ごとに違う。

## 2. スコープ

- **やる**: node のスカラー勾配 4 系統 (`ransGradient` の $k,\omega$、`speciesGradient` の $Y_s$、`passiveGradient` の $\xi$・凝縮モーメント) を、NS と同じ事前計算 LSQ 係数 (`cInt`) による gather に置き換える。周期 gather の登録、出力配列名は現行どおり。検証ハーネス (`case/09.Taylor-Green/_g0_lsq_seam/`) の GG 参照を LSQ 参照 (`lsq_merged_ref`) に切り替える。前提 plan §5.1 #7 (1)(3)(4) の引き取り。
- **やらない**: cell モード (GG のまま。ユーザ方針で cell は使わない)。回転周期 (別 plan)。軸対称×周期の継ぎ目合併 (前提 plan #7(5)。LSQ 化で半割面閉包の問題は消えるが、合併は NS と同じく片側のまま)。スカラーのリミタの変更。遷移モデル (スカラー勾配を持たない)。

## 3. 関連 docs と前提

- 現状の事実 (2026-09-26 調査、`solver_density_cuda/`):
  - LSQ 係数 `cInt[3*ilp]` は `cell_planes` CSR の incidence ごとに float32 で 3 個、`calcGradient_d_wrapper` 内の static ローカル (`calcGradient_d.cu:959`)。境界 incidence (`ip >= nNormalPlanes`) は 0 (`:822`)、stencil は実在の内部隣接ノードだけ (ghost・bvar・境界点なし)。M⁺ は倍精度、打ち切り 1e-2。
  - 適用は `lsqPreGrad_internal_d` (`:850-894`) のノード並列 gather (atomic なし)。
  - スカラー GG: $k,\omega$ は `calc_scalar_gradient_face_d` (`ransTransport_d.cu:16-75`、境界面は owner 値 k[ic0]・ω[ic0])、$Y_s$・$\xi$・モーメントは `species_gradient_d` (`speciesTransport_d.cu:641-677`、境界面は実質 owner 値)。軸対称は `A_planar`・`s*_planar`。
  - 使用先: $k,\omega$ 勾配は SST の $CD_{k\omega}$・$F_1$ だけ (移流・拡散・粘性流束は勾配を使わない)。$Y_s$ は SLAU の面再構成 (`Yd_recon`、Venkat ψ_Y は nSpecies ≥ 2)。受動種は SLAU S3 の再構成 (`dPdx_recon`、無次元化 Venkat) と、それを通じて FCT の Pface。
  - 周期 gather: `periodicGradientGather_d_wrapper` に NS 18 本・divU・dY・受動種を登録 (k/ω は `ransGradient` 直後の専用 gather)。合併係数の部分和は、group 内で状態が同値 (root→member ミラー、スカラーもミラー済み) なら和で合併 LSQ になる。
  - 軸対称の NS LSQ は係数が純幾何 (planar LSQ) で、planar GG と同じ量を推定する。

## 4. 設計方針 (2026-09-26 `diagnostician` 確定)

**前提の事実**: 既定の設定で生きているスカラー勾配は $k,\omega$ だけ。`speciesFaceReconstruction` の既定は 0 (`solverConfig.hpp:110`) で、$Y_s$・$\xi$・モーメントの勾配は `main.cpp` の `if (speciesFaceReconstruction >= 1)` の中でしか計算されない (case/*/run_* に ≥ 1 の config は無い)。したがって物理 A/B の主役は $k,\omega$ で、化学種・受動種は作用素試験と opt-in の 1 ケースで確かめる。

### 4.1 係数の共有

- `cInt` を static ローカルから外に出し、アクセサで公開する。NS と同じ配列を使う (合併済み係数を含む)。係数の再計算・変数別の係数は作らない。

### 4.2 適用カーネル

- 汎用の多変数 LSQ gather カーネルを 1 本足す: 入力 `flow_float**` (N 変数)、出力 3N 配列。ノード並列 (atomic なし)、アキュムレータは最大 4 変数ずつのチャンク (CSR と `cInt` を読み直す)。**REG は `cuobjdump -res-usage` で記録**する。
- **差分形** $\nabla\phi_i=\sum_j c_{ij}(\phi_j-\phi_i)$ (NS と同じ) にし、定数場の勾配を厳密に 0 にする (前提 plan #7(3) の GG 壁閉包 5.3ε はこれで消える)。
- 境界の扱いは NS と完全に同一 (内部隣接のみ、境界 incidence の係数 0、疑似点なし)。根拠: accepted の node 境界勾配の原則 (DOF-only) と同じで、変数ごとに境界の扱いが違う現状の方が不整合。GG の内部面値が使う幾何 fx (高 AR 曲面壁で 0.07–0.96 に振れる既知欠陥) も使わない。
- 軸対称: 係数は planar LSQ のまま (NS と同じ)。GG の `A_planar` 除算は不要。
- 壁節点の ω: 低 Re + `sstNodeWallKPin` (既定 ON) では壁ノードの $k,\omega$ の残差・対角が 0 (`ransSource_d.cu:303-307`) なので壁ノードの $CD_{k\omega}$ は解に入らず、$F_1$ は wall_dist = 0 で 1 に張り付く。影響は第一内層ノードの $CD_{k\omega}$・$F_1$ と W–I 面の拡散係数 (面の F1 平均) に限られる → S2 の case/48 で記録。

### 4.3 周期

- 並進 (`periodicSeamMergeActive`): 合併係数 + 既存の gather (和)。$k,\omega$ の専用 gather の位置は現行どおり。
- それ以外 (軸対称×周期、回転周期): gather しない (各部分 CV の片側 LSQ。NS の軸対称と同じ)。軸対称×周期は既知の制約として methods に書く。
- **node + 回転周期 (type 1) は起動時に警告** (`[config] node 回転周期は未対応 (plan boundary-node-rotational-periodic)`)。現状は速度も勾配も回さずに走ってしまう。エラー化は回転周期 plan の #0。

### 4.4 切り替え

- **最初は opt-in**: `mesh.scalarGradient: gg|lsq`、既定 `gg`。**S0–S2 が通った commit で本 plan の中で既定を `lsq` に切り替える** (後継 plan は作らない)。cell は GG 固定。
  opt-in から入る理由: 共有ワークツリーで並行セッションの node SST run と設計 DB が検証途中で黙って変わるのを避ける。
- 付随 (`slauWallNormalChi` の既定化と同じ): 起動エコー `'scalarGradient' effective: <値> (default|explicit)`、`stage_manifest` の hard キー `mesh.scalarGradient` (実効値。キー無しの旧 manifest は `gg`)、`RUN_PROVENANCE` への記録、`procedures/recommended-settings.md` §9 に「〜既定化の日付までは gg」。

## 5. 実装ステップ

1. `methods/gradient.md`・`methods/discretization.md` §7.3 を更新 (opt-in の LSQ 経路、既知の制約) — **済 2026-09-26**。
2. `cInt` の公開、汎用 LSQ gather カーネル、3 wrapper の分岐、キー・エコー・manifest・provenance、回転周期の起動警告。
3. ハーネスに LSQ 参照を追加し S0/S1、続いて S2 (AWS)、S3。
4. 既定を `lsq` に切り替え、docs 同期、codex result。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 (**完了 2026-09-26**) | §4・§6 の確定 | 判断: 2026-09-26 `diagnostician`・§4 案 1–3 採用、周期は並進のみ合併、既定は opt-in から開始、§6 の合否を固定 | F |
| 2 | codex plan 段 | Critical/Major は `diagnostician` | F |
| 3 | 実装 | §4、§5 の 1–2 | O |
| 4 | S0/S1 | ハーネス (AWS) | O |
| 5 | S2/S3 | §6 の表の 4 ケース (**すべて AWS**、ユーザ指示で 2D も AWS) | O (結論 F) |
| 6 | 既定の切り替え + docs + codex result | | F |
| 7 (**完了 2026-09-26**: `boundary-node-rotational-periodic` §5.1 #0a) | 回転周期 plan への登録 | `boundary-node-rotational-periodic` §5.1 #0 に「node + type 1 周期は回転実装まで起動エラー (現状は速度も勾配も回さずに走る)」と、NS の `periodicGradientGather` が並進を判定しない件 | O |

## 6. 検証 (2026-09-26 `diagnostician` 確定、測る前に固定)

- **S0 作用素** (ハーネス 7 変種 + 軸対称×周期 1 本): $k,\omega,\xi,Y$ の LSQ 勾配が `lsq_merged_ref` (double) と ≤ 1e-5·S。**定数場は全節点 (壁・継ぎ目含む) で厳密に 0** (差分形なので == 0 を要求。前提 plan #7(3) を閉じる)。同じ場を NS の変数に入れたときの NS 勾配と**ビット同一** (同じ係数・同じ和の順序。順序が違うなら ≤ 2 ulp)。
- **S1 非干渉**: `scalarGradient: gg` で実装前バイナリと全配列が R3 規則で一致 (旧経路の保存)。`lsq` で NS の勾配・リミタ配列が `gg` と同じ step でビット一致 (1 step、初回ダンプ)。
- **S2 物理 A/B** — 規則: **各ケースの既存の検証ゲートを `lsq` でも通すことが合格**。GG との差は記録し、差の上限を超えたら自動で不合格にせず、ケース既存の粗/細メッシュの格子対で「差が h とともに縮むか」を見る (縮めば離散化差として合格、縮まなければ欠陥として停止)。分岐は同じ収束場から `restart_field.py` で `gg`/`lsq` の 2 本、同じ step 数。両側 `check_convergence --from-floor` PASS、比較量は STEADY が前提。

  | ケース | 起点 | 量 | 合格 (既存ゲート) | 差の上限 (超えたら格子対) |
  | --- | --- | --- | --- | --- |
  | case/48 冷却平板 (低 Re SST、主対象) | `run_0005_B_tw300_y3` 最終場、+24000 | Cf・St @ x 0.3/0.6/0.9、第一内層 F1 L∞ (記録)、壁ノード F1 が両経路とも 1 | Cf/VD-II 0.97–0.99 の帯を ±3 % 内で維持 | Cf・St 各 ±1 % |
  | case/40 軸対称ノズル (node y⁺1、SST) | `run_0045_node_yp1_outletfix_cont` res_24000、+12000 | η_CF、ṁ、壁温 L∞ | case の生産ゲート (η_CF・出口列) | η_CF ±0.1 %、壁温 ±15 K |
  | case/39 周期丘 (並進周期 SST) | `run_0039` 延長の最終場、+200k | 継ぎ目比 r_gradk・r_gradw、dF1_inf、Cf 3 点、x_r | 継ぎ目比 [0.9, 1.1]、両側 segment PASS | Cf・x_r は記録のみ |
  | case/16 S3 (opt-in の Y・ξ 経路) | `compare_s3_wall_pp0.py` の既存 S3 構成 | 壁 p/p0、onset | 実験との差が GG S3 より 0.5 % 以上悪化しない | 壁 p/p0 L∞ 0.5 %、onset 1 間隔 |

  case/39 の dF1_inf はゲートにしない。`lsq` で STEADY になれば前提 plan #7(1) を「GG 経路の性質」として閉じ、DRIFTING なら継続課題 (どちらでも記録)。
- **S3 性能**: case/48 と case/39 の step 時間で `lsq` の増分が全体の 5 % 以下。超えたら既定の切り替えの前に最適化の項目を立てる (不合格ではない)。
- **前提 plan #7(4) (ジッタ格子の LSQ 次数)**: 本 plan で「欠陥ではなく作用素の性質」として閉じる。根拠は前提 plan の固定形状試験 (合併作用素の 1 次整合、継ぎ目の誤差定数 ≤ 内部) と、S0 のビット同一 (スカラーは NS と同じ係数・同じ次数の性質)。`methods/gradient.md` に記録。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`ransTransport_d.cu`、`speciesTransport_d.cu`、`periodicNode_d.cu`、`input/solverConfig.{hpp,cpp}`。
- **既存の node SST・化学種再構成・受動種 run の結果が変わる** (勾配作用素の変更)。

## 8. 完了条件

- [ ] `methods/` 更新
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — §4・§6 を `diagnostician` 判断で確定 (opt-in から開始、S0–S2 通過後に本 plan 内で既定化、既存ゲート + 格子対の合否規則)。
- `2026-09-26` — 起票 (ユーザ決定「スカラー勾配も LSQ に揃える」、前提 plan の accepted を受けて)。
