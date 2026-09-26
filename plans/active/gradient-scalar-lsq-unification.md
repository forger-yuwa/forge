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

## 4. 設計方針 (案 — `diagnostician` 判断待ち)

### 4.1 係数の共有

- `cInt` を static ローカルから外に出し、アクセサ (例 `lsqPreCoeffDevicePtr()`) で公開する。NS と同じ配列を使う (合併済み係数を含む)。
- 係数の再計算・変数別の係数は作らない (変数に依らない前提を守る)。

### 4.2 適用カーネル

- 汎用の多変数 LSQ gather カーネルを 1 本足す: 入力 `flow_float**` (N 変数)、出力 3N 配列。ノード並列、アキュムレータは最大 4 変数ずつのチャンク (レジスタを抑えるため CSR と `cInt` を読み直す)。
- 境界の扱いは NS と同一 (内部隣接のみ、境界 incidence の係数 0)。GG の「境界半割面を owner 値で積算」は無くなる。
- 軸対称: 係数は planar LSQ のまま (NS と同じ)。GG の `A_planar` 除算は不要。

### 4.3 周期

- `periodicSeamMergeActive` のとき: 合併係数 + 既存の gather (和)。k/ω の専用 gather の位置は現行どおり。
- それ以外の周期 (軸対称×周期、回転周期): 案 = gather しない (各部分 CV の片側 LSQ のまま。NS の軸対称と同じ扱い)。※調査で、NS の `periodicGradientGather` は並進を判定しないので回転周期で片側 LSQ の和 (2 倍の疑い) になっていることが分かった → 回転周期 plan の §5.1 に登録する。

### 4.4 切り替え

- 案: node の既定を LSQ にし、A/B 用に `mesh.scalarGradient: gg` で旧経路を残す (opt-in で旧経路)。cell は GG 固定。

## 5. 実装ステップ

1. `methods/gradient.md`・`methods/discretization.md` §7.3 を計画後の仕様に更新。
2. `cInt` の公開、汎用 LSQ gather カーネル、3 wrapper の置換、キー追加。
3. ハーネスの参照切替と検証 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の確定 | `diagnostician` に諮る (条件 1) | F |
| 2 | codex plan 段 | Critical/Major は `diagnostician` | F |
| 3 | 実装 | §4 どおり | O |
| 4 | 検証 | §6 | O (結論 F) |
| 5 | docs + codex result | | F |

## 6. 検証 (案 — `diagnostician` 判断待ち、測る前に固定する)

- **S0 作用素**: ハーネス 7 変種で $k,\omega,\xi$ (と $Y_s$) の LSQ 勾配が `lsq_merged_ref` (double) と ≤ 1e-5·S。定数場は全節点 (壁・継ぎ目含む) で ≤ 4ε|φ|/h (前提 plan #7(3) を閉じる)。
- **S1 非干渉**: NS の勾配・リミタ配列が本変更の前後で R3 規則 (ビット一致 / atomicAdd 床)。
- **S2 A/B (GG → LSQ)**: case/48 (SST 平板、Cf・St)、case/16 (化学種・凝縮ノズル)、case/39 (周期丘、継ぎ目比・dF1 の定常性 = 前提 plan #7(1))、軸対称 1 ケース。合否ラインは未定。
- **S3 性能**: step 時間 (AWS、同一ブロック)。

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

- `2026-09-26` — 起票 (ユーザ決定「スカラー勾配も LSQ に揃える」、前提 plan の accepted を受けて)。§4・§6 は案。
