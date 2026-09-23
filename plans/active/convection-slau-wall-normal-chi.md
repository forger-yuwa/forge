# SLAU の圧力散逸を壁隣接面で面法線 Mach により評価する (`slauWallNormalChi`)

## メタ

- **area**: `convection`
- **status**: `draft`
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (SLAU の $\chi$ と質量流束)
  - [`methods/convection/implementation.md`](../../methods/convection/implementation.md)
  - [`methods/boundary/`](../../methods/boundary/) (node の壁ノード Dirichlet)
- **related_plans**:
  - [`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §4.13.1 (本件を発見した診断。run_0430–0434)
- **created**: `2026-09-23`
- **owner**: `CFD Dev`

## 1. 目的

node 方式の Dirichlet 壁ノード ($u=0$) が、隣接内点の $|\mathbf u|$ が $\sqrt2\,\hat c$ を超える面に接すると、
SLAU の質量流束の圧力差項 $-\chi(P_R-P_L)/\hat c$ が $\chi=0$ で消え、**壁 CV に質量を戻す経路が無くなる**。
結果として壁 CV が指数的に排出され、EOS 床以降は一定シンクとなって有限 step で負密度に達する
([`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §4.13.1)。
本計画は**壁隣接面に限って** $\chi$ を面法線 Mach で評価する opt-in を入れ、剥離縁を持つ有限厚の板で起動できるようにする。

## 2. スコープ

- **やる**: `space.slauWallNormalChi` (opt-in, 既定 0) の追加。node 方式・`nodeWallDirichlet: 1` で、
  内部面のうち少なくとも一端が壁ノードの面について $\chi$ を $\widehat M_n$ から作る。SLAU / SLAU2 の双方。
- **やらない**: 既定の変更 (§6 V3 が全部通ってから別途判断)。cell 方式の実装。
  EOS 床の意味の修正 ([`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B1e。**別 commit・別判断**)。
  Roe 側の改修。接続模型の格子解像度 ($\Delta x_1$) の振り直し。

## 3. 関連 docs と前提

- SLAU の $\chi$ の定義と質量流束は [`methods/convection/theory.md`](../../methods/convection/theory.md) の「SLAU」節。
  本計画は**差分のみ**を書く。
- 壁ノードのピン留めは `cuda_forge/nodeWallDirichlet_d.cu` (運動量・$\omega$・$k$・条件付きで `res_roe` を 0 化。
  **`res_ro` は触らない** = 密度行は自由)。
- 境界半割面は `convectiveFlux_boundary_d.inc.cuh:192` の $\dot m = A\rho_R(\mathbf U_b\!\cdot\!\mathbf n)$ で
  no-slip なら厳密 0。本計画の対象外。
- **前提の確認 (起票時点で未了)**: Shima–Kitamura (2011) が $\widehat M$ に速度ベクトルの大きさ $|\mathbf u|$ を使う理由
  (多次元性 / carbuncle 対策と理解しているが一次資料で確認する)。AUSM⁺-up は面法線 Mach で低 Mach スケーリングするので先例はある。

## 4. 設計方針

**キー**: `space.slauWallNormalChi` (int, 既定 0)。`space` 直下に置く ([[keepdiss-keys-toplevel-trap]] の罠を踏まないよう、
起動ログで実効値を出す)。0 のとき現行とビット同一。

**対象面**: node 方式 × `nodeWallDirichlet: 1` で、**内部面のうち少なくとも一端が `wall_flag==1` の面**。

- W↔W' (両側が壁ノード) は両側 $V_n=0$ で既に $\widehat M=0,\ \chi=1$ なので実効的な変化は無い。
- したがって実効は **W↔I 面** (壁ノード ↔ 内点)。
- 境界半割面は $\dot m=0$ なので対象外。cell 方式は実装しない。

**定義**: 現行の

$$\widehat M = \min\Big(1, \frac{\sqrt{\tfrac12(|\mathbf u_L|^2+|\mathbf u_R|^2)}}{\hat c}\Big), \qquad \chi=(1-\widehat M)^2$$

を、対象面でのみ**面法線成分**で置き換える:

$$\widehat M_n = \min\Big(1, \frac{\sqrt{\tfrac12(V_{nL}^2+V_{nR}^2)}}{\hat c}\Big), \qquad \chi_n=(1-\widehat M_n)^2 .$$

case/46 の面 575140 ($V_{nL}=0,\ V_{nR}=276.2,\ \hat c=700.5$) では $\widehat M_n=0.28 \Rightarrow \chi_n=0.52$
(現行は $|\mathbf u_R|\ge991$ より $\widehat M=1,\ \chi=0$)。補充は $\approx 0.52\times\tfrac12 A\,\Delta P/\hat c \approx 1.2$e-7 kg/s で、
排出 4.3e-9 kg/s を 2 桁上回る → 壁 CV は $P_w \to P_i$ へ向かう。

**適用箇所**: $\chi$ は質量流束 $\dot m$ と圧力束の第 3 項 $(1-\chi)(\beta_++\beta_--1)\tfrac12(P_L+P_R)$ の**両方**に入るので、
対象面では**両方**を $\chi_n$ で組む (SLAU は 1 つの $\chi$ で構成されている)。SLAU2 の第 3 項は $\chi$ を使わないので、
SLAU2 では $\dot m$ だけが変わる。

**付随リスクが小さいと考える理由**: 付着境界層では壁法線面の $\Delta P \approx 0$ なので、$\chi_n \ne 0$ でも
圧力差項はほぼ 0 のまま。効くのは**壁法線方向に大きな圧力段差がある剥離縁**だけで、これは狙った場所そのもの。
とはいえ壁散逸を増やす変更なので、§6 V3 で摩擦・熱流束の回帰を必ず見る。

**代替案 (併記。推奨は上記)**:

| 案 | 内容 | 評価 |
| --- | --- | --- |
| (b) | Roe 型の音響散逸を壁隣接面の質量流束にだけ加える | 効果は確認済み (run_0433 で符号反転) だが 2 つのスキームの混成になる |
| (c) | 既存 opt-in `time.deltaT.updateGuardAlpha` (commit 時に $dq$ を正値性で縮める) | **排出を止めず遅らせるだけ**。単独では不可 |
| (d) | EOS 床の意味を直す | 独立の欠陥。[`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B1e で別途 |

## 5. 実装ステップ

1. `input/solverConfig.{hpp,cpp}` に `slauWallNormalChi` を追加し、起動ログに実効値を出す。
2. `cuda_forge/convection/convectiveFlux_slau_d.inc.cuh` の $\widehat M$ 算出箇所 (`:538` 付近) に分岐を入れる。
   壁ノード判定は既存の `wall_flag` を流束カーネルへ渡す (無ければ `mesh` から引き回す)。
3. `methods/convection/theory.md` の SLAU 節に $\chi_n$ の節を追加し、`methods/index.md` と整合させる。
4. §6 の V1–V3 を回す。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | **codex plan レビュー** | §4 と §6 が書けた本状態で `codex_review.py <この plan> --stage plan`。Critical / Major の採否を §6.1 と本表に反映してから実装に入る | F |
| 2 | 出典の確認 | `papers/` に SLAU の一次資料が無い。$\widehat M$ に $|\mathbf u|$ を使う理由を確認し §3 に追記 (多次元性 / carbuncle 対策の理解が正しいか)。**実装前** | F |
| 3 | 実装 (§5 の 1–2) | 触るファイル: `input/solverConfig.{hpp,cpp}`, `cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`。合格: ビル ド成功 + flag 0 で V2 がビット同一 | O |
| 4 | V1 (本件の検証) | case/46 接続模型。合格条件は §6 のとおり | O |
| 5 | V2 / V3 (無害性と回帰) | case/48・case/16・SERN 2D 生産。1 つでも外れたら**既定化しない** | O |
| 6 | docs 同期 | `methods/convection/theory.md` に $\chi_n$、`methods/index.md` の目次 | O |
| 7 | codex result レビュー | V1–V3 の VERDICT が出そろってから `--stage result` | F |

## 6. 検証

- **単体 / ビルド**: `slauWallNormalChi: 0` で既存バイナリとビット同一を確認 (V2)。
- **検証ケース**: [`../../procedures/verification/README.md`](../../procedures/verification/README.md) の選定に従い、
  下の V1–V4。

**判定基準 (結果を見る前に固定する)**:

| # | 内容 | 合格 | 不合格 |
| --- | --- | --- | --- |
| **V1** | case/46 接続模型。`run_0430` の `res_1800` 起点、SLAU + $\chi_n$、**6000 step 以上**、200 step 出力 + 系列 CSV | $\rho$(CV 153797), $\rho$(CV 189814) が**常に $10\rho_{Min}$ 以上**を保ち、床到達 0、`check_quasisteady.py --series-csv` が **STEADY** (OSCILLATING なら振幅 < 20 %)、終端で $P_w/P_i \in [0.3, 1.5]$ | 下向き DRIFTING / NaN / $P_w > 2P_i$ (過補正) |
| **V2** | 無害性。`slauWallNormalChi: 0` で case/48・case/16 | `check_field_regression` で現行と**ビット同一** | 差があれば実装の誤り |
| **V3** | flag 1 の回帰。case/48 冷却平板 (y⁺≈1)・case/16・SERN 2D 生産 | case/48 の $C_f$・$q_w$ が flag 0 と **±1 %**、case/16 の壁 $p/p_0$ **±0.5 %**、SERN 2D の $C_T$ **±0.1 %** かつ起動レシピ通過 | 1 つでも外れたら**既定化しない** (opt-in のまま残す) |
| **V4** (任意) | [[base-wake-resolution-rule]] の旧ベース形状 | 壁 CV の排出率が 0 になる | — |

いずれも `check_convergence.py` の同一設定区間 VERDICT を併記する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | (未実施) | | | §5.1 #1。**実装前に回す** |

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{hpp,cpp}`
- `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`
- `methods/convection/theory.md`, `methods/index.md`
- 既定 0 なので既存ケース・手順への影響は無い (V2 で担保)。

## 8. 完了条件

- [ ] 関連 `methods/convection/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 の V1–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status` を `done` に変更し §9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-09-23` — 初稿。[`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §4.13.1 の診断 (run_0430–0434) と
  `diagnostician` の設計方針を受けて起票。**実装前に codex plan 段が必要** (§5.1 #1)。
