# SLAU の圧力散逸を壁隣接面で面法線 Mach により評価する (`slauWallNormalChi`)

## メタ

- **area**: `convection`
- **status**: `in_progress`
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
  内部面のうち少なくとも一端が壁ノードの面について、**質量流束 $\dot m$ の $\chi$ のみ** $\widehat M_n$ から作る
  (`chi_mass = χ_n`, `chi_pressure = χ`)。SLAU / SLAU2 の双方。
- **やらない**: **圧力束の $\chi$ の変更** (§4.2 で却下)。既定の変更 (§6 V3 が全部通ってから別途判断)。cell 方式の実装
  (明示的に無効化する)。EOS 床の意味の修正 ([`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B1e。
  **別 commit・別判断**)。Roe 側の改修。接続模型の格子解像度 ($\Delta x_1$) の振り直し。
- **制限事項 (flag 1 の適用範囲)**: 本計画で検証するのは **非周期・非軸対称の node 構成** (§6 の V1–V3 のケース) に限る。
  $\chi_n$ は面局所量で体積ソースを持たず、mask は `wall_flag` (周期 seam の両メンバに存在) なので保存性は変わらない
  見込みだが~~**未検証**。**周期・軸対称で flag 1 を使う前に小規模 node 試験を要する** (§5.1)。~~
  **小規模 node 試験は PASS** (§6 V6、2026-09-24、§5.1 #7)。**生産規模の周期・軸対称では未使用** (制限の解除ではなく範囲の明示。2026-09-25)。

## 3. 関連 docs と前提

- SLAU の $\chi$ の定義と質量流束は [`methods/convection/theory.md`](../../methods/convection/theory.md) の「SLAU」節。
  本計画は**差分のみ**を書く。
- 壁ノードのピン留めは `cuda_forge/nodeWallDirichlet_d.cu` (運動量・$\omega$・$k$・条件付きで `res_roe` を 0 化。
  **`res_ro` は触らない** = 密度行は自由)。
- 境界半割面は `convectiveFlux_boundary_d.inc.cuh:192` の $\dot m = A\rho_R(\mathbf U_b\!\cdot\!\mathbf n)$ で
  no-slip なら厳密 0。本計画の対象外。

### 3.1 文献 — 面法線マッハ版は既に研究され、**全面適用では悪化が報告されている**

**これは本計画の最大のリスクなので、先例の話より先に書く** (2026-09-23 調査。codex plan 段 M1 の指摘を受けて実施)。

- **[Furusawa & Kitamura (2023), IJNMF 95(6) 992–1010](https://onlinelibrary.wiley.com/doi/10.1002/fld.5183)**
  "Stability effect of multidimensional velocity components in numerical flux SLAU":
  **`mSLAU` = SLAU の速度を面法線成分だけにした版**が既に存在し (回転翼計算で使用)、本論文がその影響を調べている。
  結論は **SLAU の多次元速度成分が、質の悪い格子に対する安定性に寄与している** (等方的に十分な数値散逸を作るため)、
  **とくに低亜音速と超音速で**。`mSLAU` は minmod と組めば中程度のマッハ ($0.1<M<1.0$) は実用になるが、
  **面法線成分のみの使用は収束を悪化させ、格子形状に敏感になりうる**。
- **本計画との関係 (重要)**: `mSLAU` は**全面適用**である。本計画は**壁隣接面に限った opt-in の局所適用**なので、
  報告された悪化がそのまま当たるとは限らない。しかし**悪化が報告された条件 (超音速・質の悪い格子) は本件の条件そのもの**
  (case/46 は M6 級、壁層の AR は最大 676)。→ **§6 V3 の「収束の悪化」と「格子感度」を明示的に見る**。
  局所適用で回避できるかどうかが、本計画の賭けの中身である。
- **一次資料 (Shima–Kitamura 2011) の §III.K / Fig.20 が面法線マッハ版を比較している**と codex plan 段で指摘されたが、
  **原典は有料で未読** (ResearchGate / AIAA とも 403)。上記 2023 年論文が同じ著者グループによる追跡研究で、
  問いも同一なので、**当面はこれを出典として扱う**。原典は `result` 段の前に入手を試みる (§5.1)。
- なお `mSLAU` について「経験的に安定と主張されているが、面法線成分のみの使用は数値不安定を引き起こすと指摘されている」
  という二次的な記述もある。

### 3.2 SLAU2 は本件に効かない (forge 内の既存 A/B が裏づけ)

**SLAU2 (Kitamura–Shima 2013) が変えるのは圧力束の第 3 項だけで、質量流束 $\dot m$ は SLAU と完全に同一**
(`convectiveFlux_slau_d.inc.cuh:546-551`、[`methods/convection/theory.md`](../../methods/convection/theory.md) SLAU2 節)。
SLAU2 の狙いは衝撃波ロバスト性 (カーバンクル抑制) と $(1-\chi)$ 依存の除去であり、超音速域では両式は実質同等。
本件の排出は**質量流束の $\chi$** にあるので、**`solver: SLAU2` に替えても直らない**。

これは推測ではなく、**forge で既に同型の実験が済んでいる**:
[`convection-slau2-lowmach.md`](../accepted/convection-slau2-lowmach.md) §9 は、3D ノズル (陰解法・同条件・2000 step) で
SLAU と SLAU2 を A/B し、`rms_roe` フロア 5.05 対 5.06、チャンバー圧 std/mean 0.654 % 対 0.656 % と**ほぼ同一**で、
**低マッハの圧力–速度カップリングの問題は SLAU2 では解消しない**と結論している。理由も同じで、
そのカップリングを担うのが **SLAU/SLAU2 共通の質量流束項 $-\chi(P_R-P_L)/\hat c$** だから。

→ **症状は違う (低マッハ市松 / 壁 CV の排出) が、律速している項は同一**。本計画が質量流束の $\chi$ を狙うのは
この既存の切り分けとも整合する。**本計画の flag は SLAU / SLAU2 の双方で同じ効果を持つ**
(どちらも $\dot m$ の $\chi$ だけが変わる。SLAU2 の第 3 項は元々 $\chi$ を使わない)。

### 3.2 forge 内の先例 (構造としての先例。根拠の主役ではない)

AUSM⁺-up (Liou 2006, JCP 214) は質量流束の圧力拡散項の切替を面法線速度で組む。forge 自身がそのカーネルを持つ
(`cuda_forge/convection/AUSM_d.cu:277-290`):

```
U_L = (u_L·n),  U_R = (u_R·n)                 // 面法線成分
M_bar = ... (U_L, U_R から)                    // ← 速度ベクトルの大きさではない
M_p = -Kp/fa * max(1 - sig*M_bar^2, 0) * (P_R-P_L)/(ro_half*c_half^2)
```

「質量流束の圧力拡散を面法線マッハで切る」構造自体は AUSM 系の既存の作法である。
**ただし §3.1 のとおり、SLAU に対して同じことを全面で行った版 (`mSLAU`) には悪化報告がある**。
したがって**「AUSM⁺-up と同型だから妥当」とは主張しない** (codex M1 の指摘を採用)。

- **注意 (休眠コードの疑い)**: `AUSM_d.cu:277` の `M_bar` は Liou の $\bar M^2=(u_L^2+u_R^2)/(2a_{1/2}^2)$ と
  式形が合わない (sqrt の位置)。AUSM は `convectiveFlux_d.cu` の dispatch に無い休眠コードなので本計画に影響は無いが、
  **引くのは「面法線で切る」という構造であって、この実装の係数ではない**。AUSM を復活させるときは要確認。

## 4. 設計方針

**キー**: `space.slauWallNormalChi` (int, 既定 0)。`space` 直下に置く ([[keepdiss-keys-toplevel-trap]] の罠を踏まないよう、
起動ログで実効値を出す)。0 のとき現行とビット同一 (演算式として。場のビット一致は §6 V2 参照)。

**対象面**: node 方式 × `nodeWallDirichlet: 1` で、**内部面のうち少なくとも一端が `wall_flag==1` の面**。
境界半割面は $\dot m=0$ なので対象外。cell 方式は**明示的に無効化**する。

**定義**: 現行の $\widehat M = \min(1, \sqrt{\tfrac12(|\mathbf u_L|^2+|\mathbf u_R|^2)}/\hat c)$、$\chi=(1-\widehat M)^2$ を、
対象面でのみ**面法線成分**で置き換える:

$$\widehat M_n = \min\Big(1, \frac{\sqrt{\tfrac12(V_{nL}^2+V_{nR}^2)}}{\hat c}\Big), \qquad \chi_n=(1-\widehat M_n)^2 .$$

case/46 の面 575140 ($V_{nL}=0,\ V_{nR}=276.2,\ \hat c=700.5$) では $\widehat M_n=0.28 \Rightarrow \chi_n=0.52$
(現行は $|\mathbf u_R|\ge991$ より $\widehat M=1,\ \chi=0$)。

**$\widehat M_n$ の入力は「流束が実際に消費する面速度」**: ピン留めした節点値そのものではなく、
`interp_dispatch` による**再構成・補正後の面状態** (`convectiveFlux_slau_d.inc.cuh:185-198`) を使う (codex M5)。

### 4.1 適用箇所は**質量流束のみ** (圧力束は現行のまま)

$$\chi_{\text{mass}} = \chi_n, \qquad \chi_{\text{pressure}} = \chi \ (\text{現行のまま})$$

$\chi$ は質量流束 $\dot m$ と圧力束の第 3 項 $(1-\chi)(\beta_++\beta_--1)\tfrac12(P_L+P_R)$ の両方に入るが、
$\chi_n \ge \chi$ なので**圧力束側は係数 $(1-\chi)$ が小さくなり散逸が減る** — 補充機序 (質量) とは**別作用**である。
面 575140 の記録状態で代数評価すると、$\dot m$ は $+4.35$e-9 → $-1.17$e-7 kg/s と反転する一方、
**面圧力 $\tilde p$ は 807 → 1070 Pa (+32.5 %)** になる (codex の代数と当方の手計算が一致)。
これは V3 の $C_T \pm0.1$ % を確実に壊すので、**初回は質量流束に限定する**。

- **分離してよい根拠**: SLAU 族は「1 つの $\chi$」で閉じた性質に依存していない。**SLAU2 が既に第 3 項から $\chi$ を外している**
  (`convectiveFlux_slau_d.inc.cuh:546-551`)。質量流束と圧力束の散逸を独立に設計するのは同族の先例である。
- **適用範囲の正確な記述 (codex result M5)**: mask は `wall_flag` だけで判定しており、**剥離を検出していない**。
**壁隣接内部面すべて**が対象で、付着境界層でも壁法線面に接線速度差があれば $\chi_n \ne \chi$ になる。
flag の有無による流束差を決めるのは**面を跨ぐ圧力差 $\Delta p$** であって剥離の有無ではない
(単体試験 `cad/test_diag_wall_cv_budget.py`: $\Delta p$ = 1 / 10 / 100 Pa で差が 8.11e-4 / 8.11e-3 / 8.11e-2 と**正比例**、
$\Delta p$ = 0 なら接線速度があっても**差は厳密に 0**)。
→ 付着境界層で影響が小さいのは $\Delta p \approx 0$ だからであり、「剥離縁だけに効くから」ではない。

**付随する流束**: $\dot m$ が壁向きに反転すると運動量・エネルギーの風上項 $\tfrac12(\dot m-|\dot m|)U_R,\ h_R$ が壁 CV に入るが、
  運動量は Dirichlet 射影 (`nodeWallDirichlet_d.cu:57-59`) で捨てられ、等温壁は `roe` もピンで捨てる。
  **断熱壁は `roe` が自由なので質量と一緒にエネルギーが入る** (物理的に整合)。
  化学種は `massflux[ip] = mdot` (`:565`) 経由で自動的に同じ $\dot m$ を使う。

### 4.2 平衡壁圧 — 機序を説明する**近似モデル** (合否判定には使わない)

**この式は機序の説明用であり、V1 の合否には使わない** (codex M2 を採用)。壁 CV の排出と補充が釣り合う点を、
面 575140 だけの収支 (W↔W' を無視、$\rho_w \ll \rho_i$、$\rho_w=P_w/(RT_w)$) で解くと:

$$\frac{P_w}{P_i} \approx \frac{1}{1 + \dfrac{2\,V_{n,i}\,\hat c}{\chi_n\,R\,T_w}} \qquad (\chi_n=0.52 \Rightarrow 0.313)$$

**この近似は使えない**: (i) 移流項の正確形は $A\frac{\rho_w\rho_i}{\rho_w+\rho_i}V_{n,i}$ で、調和因子込みだと平衡比は **0.406**。
(ii) より重要なのは、**無視した W↔W' の補充が排出の約 65 %** ある (run_0430 記録: 内点面の排出 4.345e-9 に対し
壁同士の補充 合計 2.841e-9)。CV 189814 では壁同士の補充が内点面の排出を**上回る**。
→ **V1 は全接続面の収支で判定する** (§6 V1-c)。

**それでも残る含意 (記述として重要)**:

- **壁 CV は「満たされる」のではなく、$V_{n,i}/\hat c$ で決まる高さで止まる**。$V_{n,i}\to0$ (後流に再循環が立つ) になれば
  $P_w\to P_i$ に戻る。「$P_w \to P_i$ へ向かう」という書き方はしない。
- **リスクが残る場所**: 壁隣接面に $V_n\ne0$ と大きな $\Delta p$ が同居する所 = **衝撃の壁面衝突点**
  (SERN 2D のカウル衝撃がランプに当たる所)。§6 V3 で**衝撃足の壁圧分布**を明示的に見る。

### 4.3 検討して採らなかった案

| 案 | 判断 | 理由 |
| --- | --- | --- |
| **圧力束にも $\chi_n$ を適用** | **却下** (当初案から変更) | 面圧力が +32.5 % 変わり V3 の $C_T$ ±0.1 % を壊す。補充機序と別作用。**V1 が質量流束限定で通れば圧力側の変更は永久に不要**なので「将来の拡張」としても残さない |
| **B: 壁隣接面に KEEP の ES 行列散逸を重ねる** | **却下** | (1) σ が掛かるのは**散逸項だけ**なので正味流入の閾値は σ > 4.3e-9/2.4e-7 ≈ **0.018**。(2) より本質的に、§4.2 の式で σ=0.02 なら $P_w/P_i\approx0.02$ = **今の排出しきった状態 (0.016) を固定するだけ**。σ→1 にすれば A/Roe と同等だが、それは Roe を壁隣接面に貼るのと同じで KEEP を借りる意味が無い。(3) 5 成分すべての流束が変わる (`convectiveFlux_keep_d.inc.cuh:392-400,557-570`) ので $\chi$ だけを触る本案より侵襲的。(4) SLAU の $\dot m$ が既に持つ $-(S/2)\lvert\widehat V_n\rvert\Delta\rho$ と ES のエントロピー波が二重計上になる。(5) 本案が V3 で落ちるなら B も同じ面に散逸を足すので落ちる = **fallback にならない** |
| **C: 既存の改良 AUSM 族を丸ごと採用** | **保留 (出典として吸収する方向)** | 質量流束の圧力拡散の切替が面法線か $\lvert\mathbf u\rvert$ かを確認して §3 に加える。丸ごと採用はしない: forge の AUSM⁺/AUSM⁺-up は dispatch に無く (`convectiveFlux_d.cu:239-286` は SLAU/SLAU2/HLLE/ROE/KEEP のみ)、TP・node 対応も未確認。**1 つの隅の問題のためにスキームを替えない**。「AUSM⁺ の密度正値保存性」は壁 CV の補充とは別の性質なので根拠に使わない |
| **D: `time.deltaT.updateGuardAlpha`** | 併用可・単独不可 | commit 時に $dq$ を正値性で縮める既存 opt-in。**排出を止めず遅らせるだけ** |
| **E: EOS 床の意味を直す** | 別 plan | [`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B1e。独立の欠陥なので**本計画と同じ commit にしない** |

## 5. 実装ステップ

**docs 更新は実装前** (codex m7)。

1. `methods/convection/{theory,implementation}.md` に $\chi_n$ の仕様を書き、`methods/index.md` と整合させる。
2. `input/solverConfig.{hpp,cpp}` に `slauWallNormalChi` を追加し、**起動ログに実効値を出す**。
   `discretization != node` または `nodeWallDirichlet != 1` のとき 1 を指定したら**エラーで止める** (黙って無効にしない)。
3. 配線: `SLAU_d` の現在の引数に壁フラグが無く、呼出し元は `cuda_forge/convection/convectiveFlux_d.cu:260`。
   **wrapper を変更対象に含め**、`wall_flag_d` をカーネルへ渡す。
4. `cuda_forge/convection/convectiveFlux_slau_d.inc.cuh` の $\widehat M$ 算出 (`:538` 付近) に分岐を入れる。
   **質量流束の $\chi$ だけ**を $\chi_n$ にする (圧力束の $(1-\chi)$ は現行のまま。§4.1)。
5. 単体確認 (§6 V0) → V2 → V1 → V3 の順に回す。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| ~~1~~ (済 2026-09-23) | **codex plan レビュー** | **GO-with-changes, C0/M5/m2。全件採用**して §2–§7 を改訂 (§6.1) | F |
| ~~2~~ (一部済 2026-09-23、**閉 2026-09-25**) | 出典の確認 | (i) **済**: [Furusawa & Kitamura 2023](https://onlinelibrary.wiley.com/doi/10.1002/fld.5183) で `mSLAU` (**全面適用版**) の 「多次元速度成分が質の悪い格子への安定性に寄与、面法線のみは収束悪化と格子感度」を確認し §3.1 に反映。**本計画は壁隣接面への局所適用**なので直接は当たらないが、悪化が報告された条件 (超音速・質の悪い格子) は本件の条件そのもの → §6 V3 に収束悪化と格子感度を見る条件を追加した。(ii) ~~**未**: 原典 Shima–Kitamura 2011 §III.K / Fig.20 は **AIAA・ResearchGate とも 403 で入手できず**。`result` 段の前に再度試みる。~~ **決着 (2026-09-25, 入手不能で閉じる)**: result 段前の再試行も DOI 10.2514/1.J050905 (arc.aiaa.org)・ResearchGate 258474939 とも HTTP 403、ローカル `papers/` にも無し。**本 plan の設計は原典 §III.K に依存しない** — 本件は壁隣接面への局所適用で、原典は §3.2 のとおり「構造としての先例」止まり。定量的な裏付けは (i) Furusawa & Kitamura 2023 と本 plan の V1/V5 の実測。(iii) ~~**未**: 改良 AUSM 族の切替が面法線か $\lvert\mathbf u\rvert$ か。~~ **決着 (2026-09-25, 同じ扱い)**: 文献の分類であって本 plan の主張 (V0–V7 の実測) に依存しない。**実装は (i) を踏まえて先行した** (既定 0・ビット同一の opt-in なので後戻り可能) | F |
| ~~3~~ (済 2026-09-23, `5a4886ea`) | 実装 (§5 の 1–4) | 触るファイル: `methods/convection/{theory,implementation}.md`, `methods/index.md`, `input/solverConfig.{hpp,cpp}`, `cuda_forge/convection/convectiveFlux_d.cu` (wrapper), `cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`。合格: ビルド成功 + V0 + V2 | O |
| ~~4~~ (済 2026-09-23) | V0 (単体) → V2 (無害性) | §6 のとおり | O |
| ~~5~~ (a/c/e 済・d は #6b) | V1 (本件の検証) | case/46 接続模型。合格条件は §6 V1-a〜d | O |
| ~~6~~ (case/48・SERN 2D 済、case/16 は #6c) | V3 (回帰・受入) | case/48・case/16・SERN 2D 生産。**不合格なら受入保留・設計へ戻る** | O |
| ~~6b~~ (済 2026-09-23) | V1-d の STEADY 確認 | `run_0439` (通算 66000 step) で **`OVERALL: ALL STEADY`**。壁節点の $\rho$ は 5 桁目まで静止、床到達 0、`rms_ro` は通算 4.8 桁低下して 2.40e-10 | O |
| ~~6c~~ (**決着 2026-09-25: #10 で済**) | ~~V3 の case/16~~ | ~~run データが AWS・手元とも残っておらずメッシュ生成から要る。V3 は case/48 と SERN 2D の 2 件で通っているので優先度は低い~~ → #10 (V2・V3 とも PASS, §6 V3 case/16 条件) | O |
| ~~6d~~ (済 2026-09-23、**2026-09-25 限定: 定常差として未確定、#14**) | $C_L$ / $C_M$ の差の切り分け | **ノイズではなく実効果**と判明 (各側 3 本の反復で床 0.005–0.07 % に対し差は床の 5–10 倍)。$C_T$ 0.031 % / $C_L$ 0.433 % / $C_M$ 0.316 %。**V3 の合否は変わらない**が、既定化の判断には「どちらが正しいか」の独立した根拠が要る (§6.2)。**本計画では既定を変えない** | O |
| ~~**7**~~ (**済 2026-09-24, 全項目 PASS**) | **周期・軸対称の小規模 node 試験** | flag 1 を周期/軸対称で使う前に必要 (§2 制限事項)。**設計は §6 V6** (2026-09-24 `diagnostician`: 「カーネルに分岐は無い。変わるのは面集合の定義とツールの面積」)。V5 の P0–P5 + P6-per / P6-ax / 周期整合の前提 | O |
| ~~9a~~ (**起動成否 3 水準 済 2026-09-24**; 2026-09-25 に #9 を分割、#16) | **格子感度 — 起動成否** | **2 水準 (110 万 / 245 万) で済**: 対策なしは両方で NaN、対策ありは両方で完走 = **定性的に格子によらない**。NaN step は 3031 / 815 で**細かいほど早い** (壁 CV 体積比 1.98、局所 dt 込みで 3.72)。**残**: 3 水準目 (粗)。`--scale` は粗い側で使えない (成長率が跳ね上がり x 比が受入検査に落ちる、`--set` は `--scale` に上書きされる)。`--set` 単独で `HX` を詰めて x 比 1.71 → 1.29 まで来たが許容 1.224 に未達 | O |
| **9b** (**未**、#11 の前提へ移した、2026-09-25 #16) | **格子感度 — 定常解・壁圧・収束率の格子比較** | 受入条件から外した (codex result-3 M4)。flag 0 は細格子で解を持たないので、flag 1 の解の格子収束は生産の格子収束列 (R5n) と同じ作業になる (R5n `run_0419/0420` は flag 0 なので転用不可)。**測るときの合格 (事前に固定)**: 同一生成コマンドを記録・同一 blocksize で 3 水準、壁圧 L2 差と ΣCFL あたり低下桁数を表にする | F |
| ~~**10**~~ (**済 2026-09-25, V2・V3 とも PASS**、§6 V3 case/16 条件の結果欄) | **case/16 の V2・V3 — 現在の受入ゲート** (**判断: 2026-09-25 `diagnostician` — 自前の段構成をやめ既製レシピ `run_user_profile.py` + 1D 等エントロピー IC を使う。真因は一様亜音速 IC + `outflow` に圧力アンカーが無いこと**。**判断: 2026-09-25 `diagnostician` (2 回目) — flag 0 `run_0505` が旧前提ゲートを FAIL したのは参照値が 2D 層流でない誤指定のため。層流のまま前提 (a) に差し替え、run_0505 から 2 本分岐 (§6 V3)**) | **§6 V3 の受入条件なので「既定化の必須条件」へ格下げしない** (codex result-2 M1: 受入ゲートは結果を見てから動かさない)。node run が全部 Euler (`slip` 壁) で $\chi_n$ の対象外 → **NS 設定を新規に組む**。メッシュも手元・AWS とも無い。出口も `outlet_statPress` で B1f の見直し対象 | O |
| ~~**10b**~~ (**済 2026-09-25, PASS**) | **V3 SERN 2D の受入判定** | 終端差は許容内 ($C_T$ 0.031 %、衝撃足 0.842 %) だが**差の時系列が `DRIFTING`** (変動 53.5 % / 25.8 % / 1.3 %)。本段を延長して差が STEADY になるか確認するまで**合格と書かない** (§6.2) | O |
| ~~**10c**~~ (**限定合格 2026-09-25**、#15: 「測定区間の変化と CFL 間差は ε 以下」まで。固定点 (ドリフト減衰) は #11 の前提) | **固定点の判定** | `run_0442` と `run_0439` の終端場差 0.011 % は固定点の証明ではない。両 CFL で全保存量の収束判定 + 独立に延長した末尾区間の比較が要る (§6.2) | O |
| ~~**10d-1**~~ (済 2026-09-24) | **診断ダンプの実装** | env `FORGE_DUMP_MASSFLUX=<path>` で `var.p_d["massflux"]` を**最初の呼び出しだけ**ホストへ写して書く (`convectiveFlux_d.cu` の wrapper、カーネル本体は不変)。§6 V5 の前提。**数値の振る舞いを変えないのでカーネル変更ではない** (2026-09-24 `diagnostician` 判断) | O |
| ~~**10d-2**~~ (済 2026-09-24) | **人工状態の生成** | §6 V5 の $\Delta P$ 付き状態を `fp_y1_12um.h5` に書くスクリプト。$\Delta P=0$ の帯を含めること | O |
| ~~**10d**~~ (**済 2026-09-24, P0–P5 全て PASS**) | **実カーネル照合 (V5 P0–P5)** | 単体試験は**Python による式の確認**で float32・実 mask・非対象面のビット不変は未試験。**場では測らない** (1 step でもビット再現しない実測 + 陰解法では非対象ノードが動くのが正しい)。面流束 `massflux` で P0–P5 を判定する。設計: §6 V5 (2026-09-24 `diagnostician`: 「場でなく面で測れ」) | O |
| **11** | **既定化の判断** (**F**) | **前提**: #7・#9a・#10・#10b–#10d の完了 + **#9b (格子収束)・V7 の固定点 (ドリフトの減衰・漸近先の確認)・$C_L$/$C_M$ の準定常確認 (両 run STEADY) → その後に独立精度評価** (2026-09-25 #14–#16 で受入から移した) + $C_L$ 0.433 % / $C_M$ 0.316 % / 衝撃足の壁圧 0.842 % の**正否を独立に判定** (格子収束・別スキーム・実験/文献のいずれか)。これが無いうちは **opt-in のまま**。既定化は本 plan のスコープ外 (§2) | F |
| 8 | codex result レビュー (**3 回目の発火条件**) | **受入セットが揃ったら**回す: case/16 の flag 0/1 完了 + #10b の時系列 VERDICT + #10c の固定点判定 (延長 run)。**揃う前に出すと同じ M1 で NO-GO になりレビューを 1 周捨てる** (2 回目の教訓) | F |
| ~~**12**~~ (**済 2026-09-25**: `test_stage_manifest_wall_normal_chi.py` 7/7 PASS、`test_gate_bad_input.py` PASS。**同時に直した既存の欠陥**: 正規表現の保険経路が flow 形式の区切り `,` `}` を値に含め、flow/block や末尾キーの有無で区間が割れていた) | **codex result-3 M5: `stage_manifest` の hard キー漏れ** (優先 1) | `solver_density_cuda/tools/stage_manifest.py` の `YAML_HARD_PATHS` に `space.slauWallNormalChi` を追加し、`stage_key` でこのキーだけ省略 ≡ `0` にする (`limiterScaled`/`venkatK` の同じ穴は本 plan のスコープ外、備考のみ)。**合格 (測る前に固定)**: 単体試験 3 ケース — 省略 vs `0` → 同一 key・1 区間 / `0` vs `1` → 2 区間 / flow 形式と block 形式で同一 key。`test_gate_bad_input.py` も PASS。ローカル | O |
| ~~**13**~~ (**済 2026-09-25, V1-b PASS**: (a) `run_0437/res_0` で 3 CV とも $\Sigma\dot m$(流出正) −1.20e-7 / −1.35e-7 / −2.97e-7 kg/s = 正味流入; (b) `run_0439` res_24000/30000/36000 で $\lvert\Sigma\dot m\rvert\Delta t_l\cdot 6000/(\rho V)$ **最大 0.060 %/dump** (< 1 %)。$\Delta t_l$ は 1 step 継続で出力 (1.18–1.84e-8 s)。証拠 `case/46.sern_design/_r3_m1m3/`。**限定**: $\Delta t_l$ を出したバイナリは run_0439 当時の再ビルド (sha256 7a149acc、当時は c8c0b91b); ソルバ前処理 (床) は診断ツールに入っておらず「床から 15 倍以上離れている」ことで代える (§6.2 既記載)) | **codex result-3 M1: V1-b の未測定** (優先 3、AWS) | (i) **済 2026-09-25**: 受入結論の「V1-b 観測/予測 1.000」は V1-c の量なので直す。対応表 R2-M4 を「未完 → 測定後に再判定」に戻す。(ii) 測定: (a) `cad/diag_wall_cv_budget.py --faces-all` を `run_0437` の第 1・第 2 dump に当て、3 CV の $\Sigma\dot m$ の符号 (負 = 正味流入) を記録。(b) `run_0439` の res_24000/30000/36000 から `restart_field.py` で 1 step 継続、`output: {level: 2, extraFields: [dt_local]}` で $\Delta t_l$ を出し、$\lvert\Sigma\dot m\rvert \Delta t_l \cdot 6000/(\rho V) < 1\,\%$/dump を判定。**合格式は :233 のまま** (新設しない)。(a) で符号が正なら V1-b 不合格とし、補充機序の記述を V1-c/V5 の範囲に縮める。AWS は #14 と同じ 1 回の起動で (起動前に `aws_instance.sh status`/`busy` で他セッションの使用を確認) | O (判定 F) |
| ~~**14**~~ (**済 2026-09-25, 各 run STEADY**: `_v3sern/v3s_flag{0,1}_ext` 72 dump を `v3sern_foot_perrun.py` で系列化、`check_quasisteady.py --series-csv --drift 0.002 --osc 0.005 --tail 0.4` で **両 run とも `ALL STEADY`** (列 = 絶対 $x_{foot}$・固定窓平均 $p_w$・窓内 3 点 $p_w$)。末尾 29 dump の $x_{foot}$ span 0 (格子間隔 2.19e-4 m の 0 倍)。証拠 `case/46.sern_design/_r3_m1m3/PERRUN_QS_flag{0,1}.txt`。**判定列の変更 (理由つき)**: `diagnostician` 案の「各 dump vs 自 run 末尾平均の相対 L2」は偏差の系列で `fluct=span/mean` を当てると #10b と同じ operand 取り違えになるので、run の量そのもの (位置・圧力) に当て、自 run L2 は参考列 `selfL2_pct` に残した。$C_L$/$C_M$ の準定常は #11 の前提のまま) | **codex result-3 M3: SERN 衝撃足の各 run 準定常 + $C_L$ の「実効果」** (優先 3、AWS) | (i) **済 2026-09-25**: :974 の「実効果と確定」を撤回し「現区間の末尾平均差 0.437 % / 0.320 %、flag0 は `DRIFTING` で定常差として未確定」へ。#11 の前提に「$C_L$/$C_M$ の準定常確認 (両 run STEADY) → 独立精度評価」を明記。対応表 R2-M1 を「$C_T$ 済 / 衝撃足・$C_L$ 未」に分割。(ii) 測定: `case/46.sern_design/v3sern_series.py` に**各 run の絶対 $x_{foot}$** と**衝撃足窓 $[x_{foot}-2t, x_{foot}+5t]$ の $p_w/p_{ref}$ の「各 dump vs 自 run 末尾平均」相対 L2** を列追加し、`_v3sern/*_ext` の 72 dump に `check_quasisteady.py --series-csv --drift 0.002 --osc 0.005` を**各 run に**当ててから差ノルムを評価。**合格**: 両 run とも STEADY ($x_{foot}$ は格子 1 つ以内で不変)。OSCILLATING なら平均±振幅 | O |
| ~~**15**~~ (**済 2026-09-25**: V7 結果・V3 箇条・#10c・#11・対応表を書き換え) | **codex result-3 M2: V7 固定点の断定** (AWS 不要) | (i)「CFL に依存する固定点」を撤回し「**測定区間 (ΣCFL 13200–17200) の場の変化と CFL 間差は全列 ε 以下。ドリフト率は両 CFL で一致 (壁 ro 2.2e-7 vs 2.25e-7 %/ΣCFL) だが減衰・漸近先は未確認**」に限定。(ii) V3 の「固定点確認として `cfl_pseudo` 0.4 を 1 本」を「**ΣCFL を揃えた CFL 間差 ≤ ε**」(測った量) に書き換え、**固定点 (ドリフト減衰) の確認は #11 の前提へ移す**。(iii) #10c を「限定合格」に戻し、対応表 R2-M2 を「限定」に。**後続窓 run は足さない** (最悪列は ε 到達に 10^5 ΣCFL で、測れる長さで閉じる保証が無い) | F |
| ~~**16**~~ (**済 2026-09-25**: #9a/#9b に分割、V3 表の行を分割、撤回の残存 2 箇所に取り消し線。9b は未のまま #11 の前提) | **codex result-3 M4: 格子感度** (AWS 不要、今は測らない) | (i) #9 を「**9a 起動成否 3 水準 済**」と「**9b 定常解・壁圧・収束率の格子比較 未**」に分割。(ii) V3 表「収束の悪化と格子感度」行を分ける: 収束悪化 = 済 (case/16 (c) 比 ≤ 2・case/48・SERN 2D の VERDICT)、格子感度 = 9b = 未 → **受入条件から外して #11 の前提へ** (生産の格子収束列 R5n `run_0419/0420` は flag 0 なので転用不可)。(iii) 撤回済みの「細かいほど早い」の残存 (V3 格子感度節の「定量的には…」と「含意…」) を取り消し線に。9b を測るときの合格は「同一生成コマンドの記録・同一 blocksize・3 水準で壁圧 L2 差と ΣCFL あたり低下桁数を表にする」 | F (書き直し) + O (登録) |
| ~~**17**~~ (**済 2026-09-25**: plans/README・case/46 README・plan の「5 桁一致」・case/48 の最大差 (丸め前の値を正本に)・methods の検証範囲) | **codex result-3 m6: 台帳同期** | `plans/README.md` (「本命は未確定」→ 案 A 実装済・opt-in・in_progress)、`case/46.sern_design/README.md` (「5 桁一致」→ 0.0264 % 一致、旧固定点 PASS・格子非依存に撤回注記。plan の同じ記述も)、case/48 の最大差 (codex 再計算 $C_f$ 0.001804 % / $q_w$ 0.002388 % vs plan 0.000 / 0.011 %) は **plan 側の評価スクリプトで再計算して正本を決め、どちらを採ったか書く**、`methods/convection/theory.md` 付近に「検証範囲: 非周期・非軸対称の生産構成、周期・軸対称は小規模試験のみ」を 1 文 | O |
| ~~**18**~~ (**済 2026-09-25**: 表現を限定。$\Delta\chi\cdot\Delta p$ は未測定と明記) | **codex result-3 m7: case/16 の言い過ぎ** | 受入結論の「反復ノイズ以下」→「この 500 step 試験ではノイズ床比 0.80–1.47 (許容 2 倍以内) で V3 許容内」。$\chi_n \approx \chi$ は「期待値と整合」に留める (確定するなら run_0506 最終 dump の壁隣接面の $\Delta\chi\cdot\Delta p$ を測る。ツールが case/46 専用なら「未測定」と書いて終わり、任意) | O |
| 8b | **codex result 4 回目** | #12–#18 の後。対応表を「閉 / 限定 / 未」の 3 値に直し、**出す前に受入結論の各文が「閉」行だけを根拠にしているか 1 行ずつ照合する** (3 回とも「測れた範囲を超えて書いた」型で落ちた) | F |

## 6. 検証

- **検証ケース**: [`../../procedures/verification/README.md`](../../procedures/verification/README.md) の選定に従う。
- **基準 run は flag 0 の同一バイナリ・同一 commit で取り直す** (旧 run を基準にしない。codex M4)。

**判定基準 (結果を見る前に固定する)**:

### V0 単体 (実装直後)

- 面向き反転で $\dot m$ が符号反転のみ (SLAU / SLAU2 双方)。
- 等状態 ($L=R$) で $\dot m$ が厳密に移流のみ。
- **非対象面 (両端とも非壁) の流束が flag 0 とビット同一**。
- cell 方式で flag 1 を指定したらエラー終了。

### V1 case/46 接続模型 (本件の検証)

`run_0430` の `res_1800` 起点、SLAU + $\chi_n$、**6000 step 以上**、200 step 出力 + 毎 step probe
(CV 153797 / **153880** / 189814 / 1135684) + 系列 CSV。**入力の密度や床は変えない** (交絡。codex M3)。

| 記号 | 内容 | 合格 |
| --- | --- | --- |
| **V1-a** | 監視 CV の $\rho_w$ の時系列 | 単調増加 → プラトー |
| **V1-b** | 収支の符号。診断ツール `diag_wall_cv_budget.py` に**ソルバ前処理 (床・`pMin`・$T_{eos}$)** を足した版で、各 dump の監視 CV の**全接続面**の $\Sigma\dot m$ | 起点直後の dump で $\Sigma\dot m<0$ (正味流入)。プラトー後は $\lvert\Sigma\dot m\rvert\,\Delta t_l/V$ が $\rho_w$ の **1 %/dump 未満**。**これを $dq$ と同一視しない** (RHS 収支と陰解法の更新は別) |
| **V1-c** | 機序の定量。同じツールで、各 dump の隣接状態を凍結して $\Sigma_{\text{全面}}\dot m(P_w)=0$ となる $P_w$ を 1 次元求根 (調和因子込み・W↔W' 込み) | STEADY 後の最後の 3 dump で**観測 $P_w$ が予測の ×/÷1.5 以内**。外れたら「機序の理解が誤り」ではなく「**凍結近傍近似が破れた**」と読み、内点 $V_{n,i}$ の時系列を添えて再検討する |
| **V1-d** | 床と定常性 | **全壁ノードで床到達 0**。`check_quasisteady.py --series-csv` ($\rho_{153797},\rho_{153880},\rho_{189814},n_{floor},\rho_{min}$) が **STEADY** (既定 `--drift 5 % / --osc 10 %` でよい。V3 とは別) |
| **V1-e** | 回復の期限 | 153797・153880・189814 が **step 200 (第 1 dump) までに $10\rho_{Min}$ を超え、以後 V1 終了まで維持**。導出: 補充率 $\mathrm{d}\rho/\mathrm{d}t \approx (A/V)\tfrac12\chi_n(P_i-P_w)/\hat c \approx 1.0$e4 kg/m³/s → 1.7e-4 → 1e-3 に要する局所擬似時間 8e-8 s ≈ 17 step (陽的)。Roe の実測から陰解法の減衰 ≈2 倍、$\chi_n$ は Roe の約半分 → **≈50 step**、安全率 4 で 200 step。**Roe の飽和見積り 47 step から取らない** |
| (打切り) | — | **1000 step を超えても届かなければ補充係数の見積りが 1 桁違う**ので設計へ戻る (合否とは別の打切り条件) |

「常に $10\rho_{Min}$ 以上」とは**書かない** (起点が 1.711e-4 で要求 1e-3 の 17 % しかなく初期時点で不合格になる。codex M3)。

### V2 無害性 (flag 0)

- **面流束演算のビット同一性** (V0 で確認) と、**場の回帰**を分ける (codex m6)。
- 場は `solver_density_cuda/tools/check_field_regress.py` で、**同一バイナリの反復からノイズ床を作り、その何倍かで判定**
  (残差は `atomicAdd` で集積するので**無変更でもビット一致しない**)。**両側 3 本以上** ([[noise-floor-both-sides]])。
  case/48・case/16。`--boundary` で壁出力 (`twall_*`, `qwall`) も比較する。

### V3 回帰・受入 (flag 1)

**不合格なら「受入保留・設計へ戻る」** (「既定化しない」は受入ゲートにならない。既定化はもともとスコープ外。codex M4)。

| ケース | 量 | 許容差 | 準定常の要求 (**許容差の 1/5 以下**) |
| --- | --- | --- | --- |
| case/48 冷却平板 (y⁺≈1) | $C_f$, $q_w$ | ±1 % | `--drift 0.002 --osc 0.005` |
| case/16 | 壁 $p/p_0$ | ±0.5 % | `--drift 0.001 --osc 0.0025` |

**case/16 の条件 (2026-09-25 `diagnostician`、測る前に固定)**: **既製レシピ `case/16.nozzle_wys/run_user_profile.py`**
(`--disc node --phys sst --laminar`、出口 `outflow` へ `--bc-sub`) を使う。**自前の段構成を組まない**。

- **IC は `paste_isentropic_ic`** (`design/forge_design/evaluate/ic.py:31`、1D 等エントロピー = **発散部が step 0 から超音速**)。
  **一様亜音速 IC + `outflow` は使えない** — 「$u=0$・$P=P_t$ 全域」が入口 `inlet_Pressure` と `outflow` の
  両方を満たす定常解なので、そこへ落ちる (2026-09-25 に実際に踏んだ。§6.2)。
- **`initial:` キーは run 時に効かない**。`main.cpp:1193` が h5 の VALUE を無条件に読み、
  `setInitial()` は**変換時** (`convertGmshToForge.cpp:59`) にしか呼ばれない。IC は h5 に焼き込まれたものが全て。
- **入口 Pt 59070 / Tt 286.65** (`p_0` の正規化と実験 Fig.3 がこの値。`extract_wall_pp0.py:8` の既定も 59070)。
- **出口は全段 `outflow`**。途中で `outlet_statPress` から切り替えない (新しい過渡が flag 0/1 の起動経路に差を入れる)。
  全段 `outlet_statPress` にもしない (node の亜音速壁列に背圧を課す B1f のアーチファクトが
  **$\chi_n$ = 壁列の質量流束を変える機能と交絡する**)。
- **層流でよい** ($\chi_n$ は面局所の非粘性流束の話で乱流モデルに依存しない。case/48 が SST 低 Re で既に V3 を通している)。
  **ただし壁クラスタ $y_1$ 0.5–5 µm では壁隣接内点の $M \ll \sqrt2$ なので $\chi_n\approx\chi$、差は ≈0 が期待値**
  (case/48 の $C_f$ 差 0.0018 % と同型; 旧表記 0.000 % は表示値の丸め、#17)。**V3 は「害が無い」の回帰であって効果の検証ではない** (効果側は V1/V5 が担う)。
- **flag 0 と flag 1 を別々に起動しない**。flag 0 の収束場から 2 本を分岐する (起動経路の差を消す)。
- ~~**前提ゲート (合否条件ではなく比較の前提)**: flag 0 の壁 $p/p_0$ が README の既知値
  (x=16.4 mm: 0.367 / 45.6: 0.255 / 85: 0.187) と **±2 % 以内**、出口中心線 $M \gtrsim 1.7$。~~
  **撤回 (2026-09-25, `diagnostician`: 自身の誤指定)**: この参照値は README の `run_0221` (**3D 半幅層流**、側壁閉塞あり)
  と `run_0303` (**2D node SST**) の値で、**2D 層流の既知値は README に存在しない**。「層流でよい」と同時には満たせない
  ゲートだった (run_0505 で −2.4/−3.6/−7.3 % の FAIL を出して判明)。**変えたのは前提ゲートの参照値だけで、
  V3 の合否条件 (下のノルム ±0.5 %・定常性 `--drift 0.001 --osc 0.0025`) は動かしていない**。
  **流れていない場どうしの flag 0/1 比較は自明に通る**ので、差し替え後の前提 (a) を通さずに V3 を判定しない。
- **前提 (a) (差し替え版、2026-09-25、分岐 run の前に固定)**: flag 0 分岐 run で
  出口中心線 $M \ge 1.7$; 壁 $p/p_0$ @x=16.4/45.6/85 mm が **[1D 等エントロピー, 2D SST `run_0303`] の帯内**
  = [0.3398, 0.3657] / [0.2272, 0.2574] / [0.1566, 0.1874] (下限 = 非粘性、上限 = 乱流変位。層流 $\delta^*$ < 乱流 $\delta^*$
  の序列); 3 点の系列が `check_quasisteady.py --series-csv … --drift 0.001 --osc 0.0025` で STEADY (両 run)。
- **分岐の手順**: `run_0505_node_lam_flag0/res_12000.h5` を `restart_field.py` (同一メッシュ index コピー) で 2 本に写し、
  本段設定のみ (2 次 limiter 2, cfl 2, nStepInner 5, 出口 `outflow`, 層流) を 12000 step。flag 1 は
  `space: {convMethod: 1, limiter: 2, slauWallNormalChi: 1}` と同じ flow マップに書き、起動ログ
  `'slauWallNormalChi' in 'space': 1` を確認。両 run とも `FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128`
  (REG 上限。log に出ないので `RUN_PROVENANCE.txt` に明記)。
- **(c) 収束悪化の本ケース版 (2026-09-25、分岐 run の前に固定)**: プラトーからの分岐なので低下桁数は測れない。
  代わりに **flag 1 の各 `rms_*` 末尾平均 (末尾 40 %) が flag 0 の 2 倍以内、かつ RISING 列なし**。
- V3 の差は**末尾 40 % の dump の平均**で評価する。
- ~~(b) flag 差 $L^\infty$ の系列が同閾値で STEADY~~ **訂正 (2026-09-25, `diagnostician`)**: #10b (§6.2) で潰した operand 取り違えの再発
  (`fluct = span/|mean|` を平均 2e-6 の差系列に掛けても情報を持たない)。**STEADY の要求は各 run の量 (3 点 $p/p_0$・出口 $M$) に当て、
  差は #10b 式 $\lvert\Delta(\text{末尾平均})\rvert + \text{振幅} \le 0.5\,\%$ で判定する**。V3 表の閾値と 1/5 則は不変。
- **結果 (2026-09-25) — case/16 V3 PASS** (`diagnostician` 確認済み)。`case/16.nozzle_wys/run_0506_node_lam_flag0_br` /
  `run_0507_node_lam_flag1_br` (`run_0505/res_12000` から分岐)、評価 `case/16.nozzle_wys/v3_case16/`。
  (a) 出口中心線 $M$ 1.850 (両 flag)、$p/p_0$ 末尾平均 0.358276 / 0.245869 / 0.173292 = 帯内、3 点・$M$ の系列は全 STEADY
  (`QUASISTEADY_VERDICT.txt`; f1_x85 は単調減少で漸近値 0.173287、最終比 −0.001 %)。(b) #10b 式 0.00014 / 0.00095 / 0.00258 %
  (許容 0.5 %; `B_10b_RULE.txt`)、末尾 dump の輪郭壁 $L^\infty$ 1.1–3.3e-6 (x85 の差は float32 出力の最終桁 ±0〜4 個)。
  (c) rms 比 1.0003–1.0024、RISING なし (`C_CONVERGENCE_RATIO.txt`)。両 run の `check_convergence` は plateau (分岐なので 0.1–0.2 桁、:278 で許容済み)。
  **書き方の制約**: これは「害が無い」の合格で、「一致」「効果なし」ではない。**差は V3 許容差の 1/200 以下で、同一設定反復の
  ノイズ床と同程度と推定 (床は未測定)**。flag 1 の経路が有効な根拠は起動ログ `'slauWallNormalChi' in 'space': 1` と V5 (P0–P5) で、
  case/16 からは出ない (期待値 :268 のとおり壁隣接内点 $M \ll \sqrt2$)。
- **V2 (case/16) の組み方 (2026-09-25 固定、測る前)**: case/48 と同形。**旧バイナリ = commit `39526328`** を同一マシン・同一コンパイラで
  ビルドし、旧 ×3 / 新 ×3 (flag 0) を `run_0505/res_12000` から `restart_field.py` で 500 step ずつ。両側
  `FORGE_CUDA_BLOCKSIZE=128/SMALL=128`、`output: {level: 2}` (`--boundary` が壁出力を読むため、両側同一)。
  `check_field_regress.py --repeat 旧×3 --candidate 新×3 --boundary` で**全量ノイズ床比 ≤ 2**。
  副次 (V3 の合否に影響しない、V2 と混ぜて書かない): 新 flag 1 ×3 を足し `--repeat 新flag0×3 --candidate 新flag1×3` の床比で
  「本ケースで効果がノイズ以下か」を数値で言う。
  **訂正 (実施時)**: `--boundary` が読むのは `output.level` ではなく**壁 bcond の `outputHDFflg: 1`** が書く `res_wall_3_*.h5`。両側同一に付けた
  (`level: 2` も両側同一のまま残した)。
- **結果 (2026-09-25) — case/16 V2 `VERDICT: PASS`**。`case/16.nozzle_wys/_v2/v2_{old,new}_{a,b,c}` (旧 = `39526328` を同一マシン・
  RelWithDebInfo・sm_86 でビルド、`_v2/PROVENANCE.txt`)。全量ノイズ床比 **0.89–1.58** (L2 0.96–1.45 / L∞ 0.89–1.58; 最大は
  `wall_3/twall_x` L∞ 1.58)、NaN 0 (`_v2/V2_VERDICT.txt`)。**副次 (V2 と別物)**: 新 flag 0 ×3 vs 新 flag 1 ×3 の床比 **0.80–1.47**
  (`_v2/SIDE_new0_vs_new1.txt`) = **この 500 step 試験では flag 1 と flag 0 の差が反復ノイズと同程度 (床比 0.80–1.47、許容 2 倍以内)** (m7 で表現を限定、2026-09-25)。V3 の「差はノイズ床と同程度と推定」はこれで測定済みになった
  (期待値 :268 どおり。「効果が無い」一般論ではなく、壁隣接内点 $M \ll \sqrt2$ の本ケースで χ_n ≈ χ という意味)。
- **ノルム (測る前に固定)**: $p_0$ = **入口全圧 59070 Pa**。輪郭壁ノード $x\in[10, 94]$ mm で
  $\max_x \lvert (p/p_0)_1 - (p/p_0)_0 \rvert / (p/p_0)_0 \le 0.5\,\%$ (末尾平均で評価)。
- **量の定常性**: `extract_wall_pp0.py` を全 `res_*.h5` に当て、**固定 3 点の $p/p_0$ と flag 差の $L^\infty$** の
  時系列を `check_quasisteady.py --series-csv … --drift 0.001 --osc 0.0025` で判定 (#10b と同じ手順)。
  `NOT CONVERGED (plateau)` でもこれが STEADY なら比較可、`OSCILLATING` なら平均±振幅。
| SERN 2D 生産 | $C_T$ | ±0.1 % | `--drift 0.0002 --osc 0.0005` |
| SERN 2D 生産 | **カウル衝撃がランプに当たる衝撃足の壁圧分布** (§4.2 の残リスク) | **ノルムの定義 (2026-09-23 に固定。flag 0/1 の壁圧を並べる前に決めた)**: flag 0 の run で $\max\lvert \mathrm{d}p_w/\mathrm{d}x \rvert$ を与える位置を $x_{foot}$ とし、**ランプ上 $x \in [x_{foot}-2t,\ x_{foot}+5t]$** ($t$ = カウル板厚 2 mm) の $p_w/p_{ref}$ の**相対 L2 差 ≤ 1 %** かつ **$x_{foot}$ の差 ≤ 格子間隔 1 つ分** | 同上 |
| 全ケース | **収束の悪化** (§3.1 の `mSLAU` 報告) | flag 0 と比べ $\Sigma$ CFL あたりの低下桁数が悪化しないこと。**2026-09-25 (測る前に固定, `diagnostician`)**: SERN 2D は両 run が段階起動を各自回しているので登録量で判定する — `_v3sern/v3s_flag{0,1}` に `check_convergence.py --segment` (最終区間)、**合格 = flag 1 の各列の低下桁数 ≥ flag 0 − 0.2 桁**、ΣCFL を併記。case/16・case/48 は収束場からの分岐で低下桁数が定義できないので**代替量** (case/16 は事前固定 (c) の水準比 ≤ 2、case/48 は事後の水準比)。**結果 (2026-09-25)**: SERN 2D 本段 (12000 step × `cfl_pseudo` 0.5 = ΣCFL 6000、両 run 同一) で全 8 列 OK (flag0/flag1 の低下桁数 ro 0.4/0.4、roUx 0.4/0.4、roUy 0.3/0.4、roe 0.2/0.2、roK 0.4/0.4、roOmega 0.5/0.4、roY0 0.5/0.4、roY1 0.6/0.5) → **PASS**。証拠 `case/46.sern_design/_r3_m1m3/sern2d_conv/`。**判別力の限定**: 両 run ともプラトーで低下桁数自体が 0.2–0.6 桁と小さく、許容 0.2 桁と同程度 | — |
| ~~全ケース~~ | ~~**格子感度**~~ | **受入条件から外し #11 の前提へ (2026-09-25, codex result-3 M4, §5.1 #9b/#16)**。測ったのは起動成否 3 水準 (#9a) だけ | — |

- 両側 (flag 0 / flag 1) の `check_convergence.py` VERDICT を**同一設定区間**で併記する。
- `OSCILLATING` を許す場合は**平均 ± 振幅**で比較し、**終端 dump を代表値にしない**。
- ~~block-DPLUR の固定点確認として `cfl_pseudo` を 0.2 → 0.4 にした run を 1 本足す (`nStepInner 5` 固定)。~~ **書き換え (2026-09-25, codex result-3 M2, #15)**: `cfl_pseudo` 0.2 → 0.4 の run で **ΣCFL を揃えた CFL 間差 ≤ ε** を見る (測った量)。固定点 (ドリフトの減衰) の確認は受入条件から外し #11 の前提へ移す。
- **V1 の診断成功と、生産利用可能な検証完了を区別する**。

### V5 実カーネル照合 (面流束, #10d) — 2026-09-24 に設計、結果を見る前に固定

**場 (`res_*.h5`) では測らない。面流束 `massflux[ip]` で測る。** 理由は 2 つある。

1. **場は 1 step でもビット再現しない** (実測: `case/48.flat_plate_cooled_m4/run_0030_bitrep_a` と
   `run_0031_bitrep_b`、同一メッシュ・同一 IC・同一設定・同一ブロックサイズ 128・1 step で
   `roUy` 3181/89440・`roUx` 3・`ro` 2 が不一致、他は同一)。残差は `atomicAdd` で集める
   (`convectiveFlux_slau_d.inc.cuh:638`)。
2. **陰解法では「非対象ノードの場が動かないこと」は正しい要求ではない**。flag 1 で壁ノードの RHS が
   変われば block-DPLUR の sweep が差を隣接へ運ぶ。動くのが正しい。
   一方 `massflux[ip] = mdot` は **1 面 = 1 スレッドが非 atomic に書く** (同 :580) ので、
   面レベルではビット同一が定義どおり成立する。

**試験**: `case/48.flat_plate_cooled_m4/mesh/fp_y1_12um.h5` (89440 節点、`check_mesh_quality` PASS) に
**$\Delta P$ を設計で与えた人工状態**を置き、`space.slauWallNormalChi` 0 / 1 だけを変えて **RHS の第 1 評価**
(`iStep=0, inner=0`) の面流束をダンプして比較する。**収束場は使わない**: 収束した ZPG 平板は壁隣接面の
$\Delta P$ が最小のケース ($\partial p/\partial n\approx0$, $\mathrm{d}p/\mathrm{d}x\approx0$) で、
V3 case/48 の $C_f$ 差 0.0018 % (旧表記 0.000 % は表示値の丸め、#17) がその証拠。人工状態のほうが信号が大きく、かつ**閉形式で予測できる**。

- 人工状態: 内点 $U_x = U_t$ ($M_t\approx0.6$ = $\chi$ が 0/1 にクリップされない)、$U_y = 0.2U_t r_1$、
  $P = P_\infty(1+0.05 r_2)$、$T = T_\infty$ 一様、$\rho = P/(RT)$、壁ノード $u=0$。
  $r_1,r_2\in[-1,1]$ は節点 index からの決定的疑似乱数。**$x$ のある帯では $r_2=0$** として
  「対象面だが $\Delta P=0$」の集合を意図的に作る。
- 比較対象: 面配列 `m0`,`m1` (flag 0/1)、診断ツール `diag_wall_cv_budget.py` の **float64** 値
  `m_t(flag)` と $\Delta_t = m_t(1)-m_t(0)$。
- **許容差** $\tau_b = 10^{-5} A\bar\rho\hat c$。根拠: `mdot` は約 25 flop、各 ≤0.5 ulp32 (6e-8)、
  最大中間量は $M\le1$ で $\lesssim 3A\bar\rho\hat c$ → 最悪 ≈ $5\times10^{-6}A\bar\rho\hat c$。
  $10^{-5}$ はその 2 倍で「丸めでは決して落ちない」閾値。人工状態では
  $\Delta_t \approx A\tfrac12\Delta\chi\cdot0.05P_\infty/\hat c$ なので $\Delta_t/\tau_b \approx 1800\Delta\chi \gg 100$。

| 記号 | 内容 | 合格 |
| --- | --- | --- |
| **P0** (前提) | flag 0 を 2 本回す | `m0` が**全面ビット同一**。不成立なら**試験不能**として別の非決定源を探す (**ここで合否を判定しない**) |
| **P1** | 両端 `wall_flag=0` の面 | `m0 == m1` ビット同一。**1 面でも違えば不合格** |
| **P2** | 対象面のうち **カーネルがダンプした `Ps`** が両端で float32 同値の面、および両端が壁の面 | `m0 == m1` ビット同一。**母集団を `res_0.h5` の `P` で代理してはいけない** (1 ulp ずれる。§6.2 の経緯) |
| **P3** | 対象面 $S$ (`wall_flag` ∧ $\Delta P\ne0$ ∧ $\chi_n\ne\chi$) | $\lvert\Delta_k-\Delta_t\rvert \le 2\tau_b$。かつ $\lvert\Delta_t\rvert \ge 100\tau_b$ の面で相対差 ≤ 1e-2。$S$ で `m0==m1` になってよいのは $\lvert\Delta_t\rvert < \mathrm{ulp}(m_0)$ の面だけ |
| **P4** | flag 0・1 各々、全内部面 | $\lvert m_k - m_t\rvert \le \tau_b$。**flag 0 の非対象面で落ちたらカーネルの不合格ではなく状態不一致** (ダンプ時刻と `res_0` の書かれる位置を疑う) |
| **P5** | ブロックサイズ | flag 1 を `FORGE_CUDA_BLOCKSIZE=256` でもう 1 本。`m1` が**全面ビット同一**。→ 面流束照合は 512 へ外挿可 (ブロックサイズが効くのは `atomicAdd` の集積順だけ) |

**前提条件**: 試験 config の実効値ログで `convMethod 0`, `limiter 0`, `reconT 0`, `lowMachThornber 0`,
`contactBlend 0`, `slauContactFloor 0`, `lowMachPrecond 0`, `sstEnergyIncludesK 0`, `model none`, 単一化学種
を確認する (1 つでも非 0 ならツールの 1 次式と like-with-like でない)。ツールの mask には
**`wall` と `wall_isothermal` の両方の physID** を渡す (`mesh.cpp:968` で両方が `wall_flag` に入る)。

**やらないこと**: 場の差からノイズ床を作って P1 を主張する / `nStepInner 5` の末尾状態の `massflux` を
比べる (第 2 評価以降は flag 0/1 で状態が違う) / 収束場を作るために数時間回す /
$S$ の**全面**が変わることを要求する ($\lvert\Delta_t\rvert<\mathrm{ulp}$ の面は変わらないのが正しい)。

### V6 周期・軸対称 (#7) — 2026-09-24 に設計、結果を見る前に固定

**V5 の P0–P5 をそのまま使い、周期・軸対称に固有の項目だけ足す。** カーネル側に新しい分岐は無い
(2026-09-24 に確認: `convectiveFlux_slau_d.inc.cuh` の `periodic`/`axisym` への言及は 97 行のコメント 1 件のみで、
`chi_mass` の分岐は `wall_flag` と $V_{n}$ だけ)。変わるのは**面集合の定義**と**ツールが使う面積 $A$** である。

- **周期**: node の主対流ループは純内部双対面のみで、**周期半割面は主ループからも境界カーネルからも除外される**
  (`convectiveFlux_d.cu:174-182` のコメント、`:340` の `if (bc.bcondKind == "periodic") continue;`)。
  DOF 同一視は残差の gather 側なので `massflux` には現れない。よって対象面/非対象面は planar と同じく
  **各面の両端 CV の `wall_flag`** で読む (partner 参照は不要)。
- **軸対称** (`axisymMethod 0`): `variables.cpp:530-` が `sx,sy,sz,ss` に $r_{face}=\max(p_{cy}, r_{floor})$ を**掛ける**。
  $V_n = U\cdot S/|S|$ は $r$ が分子分母で消えて不変、$\dot m \propto |S|$ なので
  **$\tau_b$ の $A$ に同じ $r$ 重み付き面積を使えば退化しない** ($\Delta$ と $\tau_b$ が同じ因子で縮む)。
  `axisymMethod 1` は幾何が planar のままなので V5 で既にカバーされている。

**試験**: 各 4 本 (flag0 ×2 / flag1 / flag1 bs256)、**1 step**、`output.level 1`。

| | メッシュ | 規模 |
| --- | --- | --- |
| 周期 | `case/38.channel_wmles/mesh/channel_ret550.geo` を低解像度で生成 (ylo/yhi `wall`、x/z `periodic`) | 約 1 万節点 |
| 軸対称 | `case/27.axi_nozzle_plume/mesh/axi_nozzle_2d_publicrao_plume.geo` (`wall` physID 3 / `axis` 4) | 約 2.8 万節点 |

**判定 (P0–P5 は V5 の閾値のまま。追加は 3 つ)**:

| 記号 | 内容 | 合格 |
| --- | --- | --- |
| **前提** (P0 と同格) | 周期 run の `.state` で partner 間の 6 配列 (`ro,Ux,Uy,Uz,P,sonic`) | **ビット一致**。不成立なら**試験不能** (人工状態の生成器が周期整合でない) |
| **P6-per** (幾何、run 不要) | `BCONDS/<periodic>/iCells` と並進量から partner を座標一致で組み、union-find の各 group で `wall_flag` | 混在 group **0**。あわせて周期半割面を `iPlanes` で**明示的に**全集合から除外し (`nei<0` に頼らない)、`nei<0` の面数 = Σ 全 bcond の `iPlanes` 数を assert |
| **P6-ax** | (a) $A_{solver}=\lvert\text{surfVect}\rvert\times\max(\text{centCoords}[1], r_{floor})$ で $\tau_b$ を組む。(b) 内部面のうち $\text{centCoords}[1]\le r_{floor}$ の面数。(c) `wall_flag ∩ axis_flag` | (b) 期待 **0** (0 でなければその面は絶対判定 P4 のみに掛ける)。(c) 期待 **∅** (ノズル壁は軸に触れない)。いずれも**報告必須** |

**ツール側の落とし穴 (実装前に潰す)**:

1. `_v5_compare.py` は `PLANES/surfVect` の planar 面積を使う。**軸対称でそのまま流用すると全面で係数 $r$ ずれて
   P4 が全滅する** — カーネルの不合格ではない。$A_{solver}$ に直す。
2. `_v5_make_state.py` は $r_1,r_2$ を**節点 index** から作るので、周期 partner (別 index) に別の $P$/$U_y$ が乗る。
   **折り返し座標からハッシュする**よう直す。`BAND_X` も対象メッシュの $x$ 範囲に合わせる (平板の値がハードコード)。

**やらないこと**: 収束場や複数 step (第 2 評価以降は flag で状態が違う) / seam を跨ぐ面や partner 参照を面集合に
持ち込む (存在しない) / `res_0.h5` の `P` を母集団に使う / **軸対称で planar 面積の $\tau_b$ を使って「軸近傍で
退化する」と結論する** (退化するのはツールの $A$ の取り方) / 生産規模のメッシュを持ち出す (資源以前に情報が増えない)。

**V6 結果 (#7) — 周期・軸対称とも全項目 PASS** (2026-09-24)。

メッシュ (生成コマンドごと記録。どちらも `check_mesh_quality.py` **VERDICT: PASS**):

- **周期**: `case/38.channel_wmles/mesh/channel_ret550.geo` の `nx/ny/nz` を `32/24/12` に落として gmsh →
  node 変換 (x/z `periodic`、y 両壁 `wall`)。**10,725 節点 / 33,754 面**、AR 3.1 / skew 0.000。
- **軸対称**: **最小直管を新規に書いた** ($L$ 0.2 m / $R$ 0.02 m、60×30、`axis`/`wall` 分離)。
  **1,891 節点 / 3,873 面**、AR 5.0 / skew 0.000。**case/27 は使えない**:
  `axi_nozzle_2d_publicrao_plume.geo` が `axisym_nozzle_points_public_rao.geo` を `Include` するが、
  その点データも生成器もリポジトリに無い (輪郭 CSV のみ)。

| 判定 | 周期 | 軸対称 |
| --- | --- | --- |
| 前提 (partner 間の `.state` ビット一致) | **PASS** (group 1075、最大 4 member、6 配列とも不一致 0) | — |
| P0 flag0 ×2 | PASS 0/33754 | PASS 0/3873 |
| P1 非対象面 | PASS 0 (28114 面) | PASS 0 (3569 面) |
| P2 $\Delta P=0$ 面 | PASS 0 (494 面) | PASS 0 (25 面) |
| P3 対象面 vs ツール | PASS $0.001\tau_b$ / 相対 4.496e-06 | PASS $0.000\tau_b$ / 相対 8.387e-07 |
| P4 全内部面 | PASS max $e/\tau_b$ **0.009** (0/30596) | PASS **0.013** (0/3690) |
| P5 blocksize 128/256 | PASS 0/33754 | PASS 0/3873 |
| P6-per | **混在 group 0**、`nei<0` 3158 = Σ iPlanes 3158 | — |
| P6-ax | — | (b) 軸床面 **0**、(c) 壁∩軸 **0 節点** (`bcond` 定義・$\lvert y\rvert<10^{-9}$ の両方で) |

**途中で踏んだツールのバグ (記録)**: 軸対称で最初 P4 が 3690 面中 **3630 面**・$e/\tau_b$ 最大 **3.75e8** で全滅した。
原因は $r$ 重みを掛けた `A` を**法線の正規化 `n = S/A` にも使っていた**こと (`S` は planar なので法線が $1/r$ ずれる)。
面積を **法線用 (planar) と流束・$\tau_b$ 用 ($r$ 重み) に分離**したら 0.013 に収まった。
**`diagnostician` が事前に「P4 が軸対称だけ全面で落ちたらまずツールの $A$ を疑え」と書いていた指紋そのもの**だった。

**これで §2 の「周期・軸対称では未検証」は外せる** (ただし**試した 2 つのメッシュで**、という限定つき)。

### V7 固定点の判定 (#10c) — 2026-09-24 に設計、**測る前に**固定

**主判定は絶対値、比はラベルにしか使わない。** 比だけでは 2 通りに破綻する:
比 $\approx1$ でも両方大きければ「同じ過渡を追っている」、比 $\gg1$ でも両方が許容差以下なら
「CFL 依存はあるが許容内」。どちらも比では言い分けられない。

- $D_{intra}(\text{cfl})$ = 各 run の**種 dump → 延長終端**の場の差 (= その run が静止しているか)
- $D_{inter}$ = **両延長 run の終端場どうし**の差 (= 固定点が CFL に依るか)
- **主判定**: $\max\{D_{intra}(0.2),\,D_{intra}(0.4),\,D_{inter}\} \le \varepsilon$ を**列ごとに**満たす
- **副 (ラベルのみ)**: $D_{inter} \ge 5\max D_{intra}$ なら「CFL 依存を**検出**」、未満なら
  「ドリフト床以下で**検出されず**、検出限界 $=\max(D_{inter},\,5\max D_{intra})$」。
  5 倍は §6.2 で既に使っている「床の 5〜10 倍で有意」と揃える

**許容差は新設しない**。§6 V3 の表 × 既存の「許容差の 1/5」則から変数ごとに引く:

| 列 | 対応する V3 量 | 許容差 | $\varepsilon$ |
| --- | --- | --- | --- |
| `P` | $C_T$ (壁圧の積分) ±0.1 % | | **0.02 %** |
| `ro`, `roUx/y/z`, `roe` | $C_f$, $q_w$ (壁 $\rho$ が $\rho_w u_\tau$ に直結) ±1 % | | **0.2 %** |

**5 列を 1 スカラーに束ねない** (収束判定と同じ規律)。

**ノルム**: **壁隣接集合を主、全域を副**の 2 本立て。壁集合は V1-d で床到達を数えた**同じ節点集合**を使う
(定義を借り、新設しない)。分母は**同じ集合上のその列自身の RMS**、ただし**運動量 3 成分は共通分母
$\lVert\rho\lvert\mathbf u\rvert\rVert_2$** (壁集合で `roUz`≈0 だと自分の RMS で割ると発散する)。節点数等重み。
壁集合の $L^\infty$ と該当節点も報告する (ゲートにはしない)。

**累積 CFL は揃えない** (固定点は経路に依らない)。ただし**この延長は経路独立性の試験としてはほぼ無力**
(系譜で $0.2\times66000+2000$ 対 $0.4\times33000+4000$ = 13 % 差)。**価値は $D_{intra}$ (静止の証拠) にある**。
$D_{intra}$ は窓の $\Sigma$CFL が 2000 対 4000 と違うので、同じ物理ドリフトなら 0.4 側が約 2 倍出る →
$\max$ を取り $D_{intra}/\Sigma\mathrm{CFL}$ も併記する。

**残差側の「末尾が平坦」は窓を自作せず `--from-floor` を使う**。`check_convergence.py` に末尾窓の独立判定は無く、
`--from-floor` が求める判定そのもの (全列が参照床の 1.5 倍以内)。**参照は系譜を連結した CSV** を渡す
(単体 run は `--drop 3` を満たさず REFUSED になる。AGENTS.md「CFL・反復数だけの違いは連結可」に当たる):
cfl 0.2 = `run_0437+0438+0439` (通算 66000)、cfl 0.4 = `run_0441+0442` (通算 33000)。
**連結の前提は確認済み** (2026-09-24: `run_0441`/`0442` の差は `outStepInterval` のみで数値設定は同一)。
連結しても `--drop 3` に届かない列があれば、その列は**床判定不能 = 不合格**とし、延長を足す。
**順序**: 先に `--from-floor` を両 run で通し、**通ってから**場の距離を測る。

**$D_{intra}$ は 2 点で取る** (種→+6000、種→+10000)。$D(10000)/D(6000)\gtrsim1.5$ なら線形ドリフト = 過渡で不合格、
$1\sim1.4$ で飽和なら床の揺らぎ = 静止。中間 dump が無い側は「跳びとドリフトを分離できない」と明記する。

**報告の文言を先に固定する (4 通り。これ以外の書き方をしない)**:

| 条件 | 書く文言 |
| --- | --- |
| 全部 $\le\varepsilon$、比 $<5$ | 「**同一固定点**。検出限界 $\delta$」 |
| 全部 $\le\varepsilon$、比 $\ge5$ | 「**固定点は CFL に依存する** ($D_{inter}=x$、床の $r$ 倍で有意) が V3 許容差の 1/5 以下」。**「同一」とは書かない** |
| $D_{intra}>\varepsilon$ | 「**比較精度に足る静止に達していない**」。固定点の主張は保留 |
| $D_{intra}\le\varepsilon$ かつ $D_{inter}>\varepsilon$ | 「**CFL 依存の固定点**」。**先に倍精度対照で float32 の不定性を除外してから**設計へ戻る |

**V7 結果 (#10c) — ~~「固定点は CFL に依存する。ただし全列が V3 許容差の 1/5 以下」~~** (2026-09-25)。
**限定 (2026-09-25, codex result-3 M2, #15)**: 固定点の断定を撤回する。言えるのは「**測定区間 (ΣCFL 13200–17200) の場の変化と CFL 間差は全列 ε 以下。ドリフト率は両 CFL で一致 (壁 `ro` 2.2e-7 対 2.25e-7 %/ΣCFL) だが、減衰・漸近先は未確認**」まで。線形ドリフトを過渡として扱う :490 の規則に照らすと、最悪列は ε 到達に 10^5 ΣCFL 級で測れる長さではない。固定点の確認は #11 の前提へ移した。
**「同一固定点」とは書かない** (同一の主張は全列の論理積で、検出列がある以上否定される)。

run: `run_0447_cfl02_ext` / `run_0448_cfl04_ext` / `run_0449_cfl02_ext2` (いずれも `restart_field.py` で
**ビット一致検査つき**に継いだ 1 万 step)。前提はすべて満たした:

- **系譜連結で両 CFL とも `PASS (converged)`**: cfl 0.2 (`run_0437+0438+0439`, 通算 65997) 4.1–4.6 桁 /
  cfl 0.4 (`run_0441+0442`, 通算 32998) 4.3–4.8 桁。**単体 run の「NOT CONVERGED」は既に下がった場から
  再開したことによる見かけ**で、固定点の可否とは無関係 (この誤読を一度報告して撤回した)。
  連結の前提 (`run_0441`/`0442` の差は `outStepInterval` のみ) も diff で確認。
- **延長 3 本とも `--from-floor <系譜>` で `PASS`** (参照床の ×0.79–1.02)。

**主判定 (全列 $\le\varepsilon$) — 成立**。$D_{inter}$ は $\Sigma$CFL を 17200 に揃えた対で取った:

| 集合 | | `ro` | `roUx` | `roUy` | `roUz` | `roe` | `P` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 壁 | $\max D_{intra}$ | 0.00090 | — | — | — | 0.00054 | 0.00037 |
| | $D_{inter}$ | **0.00318** | — | — | — | **0.00317** | **0.00317** |
| | 比 / ラベル | 3.53 **未検出** | — | — | — | 5.87 **検出** | 8.57 **検出** |
| 第一内部 | $\max D_{intra}$ | 0.00153 | 0.00765 | 0.00170 | 0.00049 | 0.00076 | 0.00037 |
| | $D_{inter}$ | 0.01106 | **0.03752** | 0.00533 | 0.00297 | 0.00690 | 0.00317 |
| | 比 / ラベル | 7.23 **検出** | 4.90 未検出 | 3.14 未検出 | 6.06 **検出** | 9.08 **検出** | 8.57 **検出** |
| | $\varepsilon$ | 0.2 | 0.2 | 0.2 | 0.2 | 0.2 | **0.02** |

**最大でも $\varepsilon$ の 1/5.3** (第一内部 `roUx` 0.03752 % / 0.2 %)、`P` は 1/6.3。
未検出列の**検出限界** ($5\max D_{intra}$): 壁 `ro` 0.0045 %、第一内部 `roUx` 0.0383 %、`roUy` 0.0085 %。
**比が列ごとに 5 をまたぐ**ので、「5 列を束ねない」規律をラベルにもそのまま適用した (全体で 1 ラベルにしない)。

**$\Sigma$CFL を揃えた対で確認した (A/B)**: $D_{inter}^{(13200)}$ (種同士) 0.00310 % に対し
$D_{inter}^{(17200)}$ 0.00318 %、**$\lvert\Delta D\rvert/D = 0.026 < 0.1$** →
**$\Sigma$CFL を揃えても持続する静的オフセット**。ここで $\Sigma$CFL を揃えるのは §6 V7 の「揃えない」と
矛盾しない — **経路独立性の試験ではなく、ドリフトの進み (lead 2000 $\Sigma$CFL) を除くため**である。
なお終端同士 (15200 対 17200、揃えず) は 0.00322 % で、3 者は 3 桁目までしか違わない。

**$D_{intra}$ は「床の揺らぎ」ではなく線形ドリフト**。$D_{intra}/\Sigma$CFL が両 CFL で一致する
(壁 `ro` 2.2e-7 対 2.25e-7 %/$\Sigma$CFL、第一内部 `roUx` 1.92e-6 対 1.91e-6)。飽和した床なら窓長
(2000 対 4000) に依らず等しいはずで、0.4 側が 2 倍出たのは §6 V7 が予告した挙動そのもの。
**$\varepsilon$ 基準では静止だが、ドリフトが無いわけではない** (最悪列が $\varepsilon$ に届くには $10^5\,\Sigma$CFL 要る)。
2 点チェックの比 1.42 は規則の未定義帯 (1.4 未満=飽和 / 1.5 以上=過渡) だが、率の一致からドリフト側に読む。

**壁集合で `ro`/`roe`/`P` が 3 列とも一致する (0.00318/0.00317/0.00317) のは偶然ではなく等温ピンの拘束**。
全壁が `wall_isothermal` $T_w$=1000 K で、`nodeWallDirichlet_d.cu:84-86` が $T=T_w$, $\rho e=\rho(e+ek)$,
$P=\rho R(Y)T_w$ とピンするので $\Delta P/P = \Delta\rho/\rho + \Delta R/R$。**閉包試験で検証した**:
$\lvert\Delta P/P\rvert$ 上位 10 % の 3796 節点で $\lvert (r-1)+(\Delta R/R)/(\Delta P/P)\rvert$ の中央値 1.56e-4 に対し
$\lvert r-1\rvert$ の中央値 3.81e-3 = **拘束式が 95.9 % を説明**。ただし残差は float32 丸め (1e-6 級) より
2 桁大きく、**4 % は未説明**。また点ごとには **28 % の節点で比例が 1 % 以上破れる**ので
**厳密な恒等式ではない** (L2 が 3 桁一致するのは $\lvert\Delta P/P\rvert$ 大の節点が支配するため)。

**壁集合の運動量は「評価不能」ではなく「定義上対象外」**。`nodeWallDirichlet_d.cu:56-58` の Dirichlet 射影で
残差が 0 化され $\rho u \equiv 0$ になる。**§6 V7 は plan `:140` でこの射影を認識していながら見落としていた**。
**第一内部ノード集合 (壁ノードに面で隣接する内部ノード、38266 節点) は測定後に追加した副集合**であり、
$\varepsilon$ は変えていない。**第一内部の値は全列で壁より大きい**ので、追加が結果を有利にはしていない。

**その他の観測 (判定には使わない)**:

- **一様オフセット説は否定**: $\Delta P/P$ の $\lvert\text{mean}\rvert/\text{rms}$ は壁 0.564 / 第一内部 0.605 / 全域 0.362
  (一様なら 1 に近いはず)。
- $\cos(\delta, D)$ (ドリフト方向と run 間オフセット) は `ro` +0.170 / `roe` +0.403 / **`P` +0.616**。
  いずれも未決を強制する 0.7 未満だが、**`P` は予測 (≈0.2) よりかなり高い**。
- $L^\infty$ は 3 列とも**同一 CV 20751** (`cowl_in` のみに属する壁ノード、双対体積 3.6e-12 m³、
  RMS の 8.4 倍)。**複数 bcond の交線ではない**。局在の機序は未確認。
- **壁の遅いモードは組成にある**: `roY1` の $D_{intra}(0.2)$ は 0.04539 % で `ro` の 100 倍速い
  (2.3e-5 %/$\Sigma$CFL)。**V7 の $\varepsilon$ 表は組成列を持たない**。率からは 8800 $\Sigma$CFL
  (cfl 0.2 で 44000 step) で 0.2 % に達しうる。**後付けでゲート化はしない**が、次の plan 改訂の候補とする。

### V4 (任意)

[[base-wake-resolution-rule]] の旧ベース形状で壁 CV の排出率が 0 になるか。

### 6.2 結果 (2026-09-23)

**V5 実カーネル照合 (#10d) — P0–P5 全て PASS** (2026-09-24)。run: `case/48.flat_plate_cooled_m4/_v5/`
(`v5_flag0_a`, `v5_flag0_b`, `v5_flag1`, `v5_flag1_bs256`。各 1 step、メッシュ `mesh/fp_y1_12um.h5`
89440 節点・`check_mesh_quality` **VERDICT: PASS (AR≤1000, skew≤0.90)**)。
実効設定は `forge_run.log` で確認: `convMethod 0` / `limiter 0` / `slauContactFloor 0` / `reconT 0` /
`lowMachPrecond 0` / `model none` / `energyIncludesK 0` / 単一化学種 N2 / tracer none。
ダンプは `[FORGE_DUMP_MASSFLUX] ... (call 1)` = **第 1 評価**。

```
faces 180007 (内部 177754, 境界半割 2253)  壁ノード 1001/89440
P0 flag0 x2 全面ビット同一 : 不一致 0/180007                        -> PASS
P1 非対象面 175752 面      : 不一致 0                               -> PASS
P2 対象面かつ ΔP=0  270 面 : 不一致 0                               -> PASS
P3 対象面 S 868 面         : max|Δk−Δt|/τ_b = 0.003 (許容 2)        -> PASS
   |Δt|>=100τ_b の 795 面  : max 相対差 3.237e-06 (許容 1e-2)       -> PASS
P4 flag0 全内部面          : max e/τ_b 0.014, 超過 0/177754         -> PASS
   flag1 全内部面          : max e/τ_b 0.014, 超過 0/177754         -> PASS
P5 blocksize 128 vs 256    : 不一致 0/180007                        -> PASS
```

**P2 は初版で FAIL した。結果を見てから閾値を緩めたのではなく、代理量の誤りを直した** (経緯を残す):

1. 初版は母集団を「**`res_0.h5` の `P` が両端で float32 同値**」で代理し、273 面中 **3 面で FAIL**
   ($\Delta_k = -6.11/-7.28/-6.98\times10^{-10}$ kg/s)。
2. $\chi_{mass}$ が `mdot` に入る経路は `−chi_mass/c_diss·P_del` の 1 項だけ (`:578`) なので、
   **$\Delta_k\ne0$ は「カーネル内で $P_{del}\ne0$」の十分条件**である。したがって
   「$\Delta P=0$ なのに $\chi_n$ が効いた」という読みは論理的に成り立たず、問いは「なぜ $P_{del}\ne0$ か」に限られる
   (2026-09-24 `diagnostician`)。実際 $\Delta_k$ は $-A\tfrac12\Delta\chi/\hat c\cdot\mathrm{ulp}(P)$ の予測に
   **0.994 / 1.016 / 1.030** で乗った。
3. **カーネルが読む `Ps` を同じ第 1 呼び出しでダンプして確かめた** (#10d-1)。`res_0.h5` の `P` は
   カーネルの `Ps` と **1716/89440 節点で不一致、ずれは常にちょうど 1 ulp** (`sonic` は 409 節点、同じく 1 ulp。
   `ro`,`Ux`,`Uy`,`Uz` は完全一致)。3 面は**ダンプでは `5037.4` vs `5037.4004` = +1.00 ulp** だった。
   **代理が不適だった**のであってカーネルの不合格ではない。母集団を `Ps` で定義し直すと 270 面・不一致 0。
   外れた 3 面は $S$ へ移り ($865\to868$)、P3 は同じ閾値で PASS (相対差は 1.149e-05 → **3.237e-06** に改善)。
4. **機序をコードで特定した** (行は 2026-09-24 に確認): `dependentVariables_d.cu:291` が
   $P=\max((\gamma-1)(\rho e-\rho\,ek),\,p_{Min})$ を作り、**直後の `:297` が
   $\rho e = P/(\gamma-1)+\rho\,ek$ と書き戻す** (float32 の往復)。これが流束の前に **2 回**走る:
   `main.cpp:1225` (初期化) → **`:1344` が `res_0` に $P_1$ を書く** → `:1391` (step 0 の再更新) →
   **`:1438` の流束が読むのは $P_2$**。**配列は同じ `var.c_d["P"]` でも時点が違う**。
   「同じ配列だから同じ値」は時点を無視した推論だった。
5. **機序を測定でも独立に確認した**。人工状態の $\rho e = P/(\gamma-1)+\tfrac12\rho\lvert u\rvert^2$ を EOS で
   $P$ に戻すとき、$\lvert u\rvert^2$ が節点ごとに違うので丸めが節点依存になる。壁ノードは $u=0$ で
   $\rho e$ が厳密なので丸めが出ない — 帯内の**壁 137 節点は 0 個**、**内点 11645 節点中 741 個 (6.4 %)** がずれた。
   帯内の **I–I 面 23068 中 2778 面 (12.0 %)** で `Ps` が不一致 (非対象面なので P1 の合否には無関係)。
   `roe` は入力 h5 と `res_0` で 28867 節点ちがうが**最大 2 ulp** (相対 1.48e-7) で、同じ丸めの往復である。

**波及**: 「`res_0.h5` を読めば like-with-like」は**未確認の前提だった**。以後、面流束をツールと突き合わせる
ときは**ダンプした状態**をツールの入力にする (P3/P4 も今回はそうした。影響は $\tau_b$ の 0.003 で無視できるが、形を揃える)。


**V0 単体 — 合格** (2026-09-23。codex result M3「面反転・等状態・非対象面の試験記録を追えない」を受けて
`cad/test_diag_wall_cv_budget.py` に固定した。`python3 case/46.sern_design/cad/test_diag_wall_cv_budget.py` で `ALL PASS`):

| 試験 | 内容 | 結果 |
| --- | --- | --- |
| **カーネル照合** | カーネル式の独立な書き下しと項ごとに比較 (flag on/off 双方) | 一致 |
| **codex の反例** | 壁 0 / 内点の接線 1000 m/s / $\hat c$ 700 / 面法線 0 / $\Delta p$ 900 Pa / $A$=1 m² | off: $\dot m$=**0** / on: **−0.642857 kg/s** |
| **面向き反転** | $\mathbf n \to -\mathbf n$ かつ L↔R で $\dot m$ が符号反転のみ (flag on/off) | 一致 |
| **等状態** | $L=R$ で圧力差項が厳密 0、$\dot m = A\rho V_n$ (flag on/off) | 一致 |
| **$\Delta p$ 比例** | $\Delta p$ = 1/10/100 Pa で flag の差が 8.11e-4 / 8.11e-3 / 8.11e-2 | **正比例** |
| **$\Delta p$=0** | 接線速度 400 m/s があっても flag の差 | **厳密に 0** |

最後の 2 件が **M5 の実証**: flag の差を決めるのは $\Delta p$ であって剥離の有無ではない。

**非対象面の不変性**は config 検証と mask の構成から担保する (`is_wall_face` が偽なら `slau_mdot` は従来式と同一行を通る)。
**cell 方式での拒否**は下の config 検証に含む。

config 検証は 4 経路すべてで起動時に停止する (黙って無効化しない):
`discretization: cell` / `nodeWallDirichlet: 0` / `solver: ROE` / 値 2 で、それぞれ固有のエラーを出して終了。
正しい組合せでは実効値をログに出し、キー省略時 (既定 0) は何も出さない。

~~**V0 マスクの実効範囲** — 設計の「効くのは剥離縁だけ」を面数で確認した~~
**撤回 (2026-09-24, codex result-2 M5)**: mask は `wall_flag` だけで**剥離を検出していない**ので、面数比から
「剥離縁だけに効く」「そこは元々 $\chi=0$」は導けない。**3D の面数比を 2D の力係数変化の説明に転用するのも不可**。
正しい記述は「**対象は壁隣接の全内部面。凍結した状態での流束差は $\Delta\chi \times \Delta p$ に比例し、
位置は面数比からは言えない**」(単体試験で $\Delta p$ 比例と $\Delta p$=0 で差ゼロを確認済)。
以下の面数は**事実として残す**が、影響位置の根拠には使わない:

| | 面数 | 全内部面比 |
| --- | --- | --- |
| 内部面 (総数) | 7,271,119 | 100 % |
| 対象 (一端が `wall_flag`) | 113,852 | 1.57 % |
| **実際に $\chi$ が変わる面** | 38,335 | **0.53 %** |
| うち $\chi$ の増分 > 0.1 | 21,759 | **0.30 %** |

対象面の 3 分の 2 では $\chi$ が変わらない。壁ノードは 37,956 / 2,446,572 = 1.55 %。
**この面数は run_0437 (3D 接続模型) のもので、他ケース・他格子には転用できない。**
影響の位置を主張するなら**そのケース自身の $\Delta\chi$・$\Delta p$・流束差の分布**を測る必要がある (未実施)。

**V2 無害性 (flag 0) — `VERDICT: PASS`**。case/48 の収束場 (`run_0025_B_tw300_y3_fx05/res_48000.h5`) から
**旧バイナリ 3 本・新バイナリ 3 本**を 200 step ずつ (`_v2/v2_{old,new}_{a,b,c}`)。
旧バイナリは本 plan の直前 commit `39526328` を同一マシン・同一コンパイラでビルドしたもの。
`check_field_regress.py --boundary` で**全量がノイズ床比 0.86〜1.36** (許容 2 倍):

| 量 | 比 (L2 / L∞) | | 量 | 比 (L2 / L∞) |
| --- | --- | --- | --- | --- |
| `ro` | 1.03 / 1.18 | | `wall_4/qwall` | 1.03 / 1.10 |
| `roUx` | 1.03 / 1.16 | | `wall_4/twall_x` | 0.99 / 1.00 |
| `roe` | 1.03 / 1.16 | | `wall_4/twall_y` | 1.08 / 1.26 |
| `P` | 1.03 / 1.28 | | `wall_4/utau` | 1.01 / 1.36 |
| `roK` / `roOmega` | 0.99–1.02 / 1.03–1.33 | | `T` | 1.02 / 1.05 |

比が 1 前後 = **同一バイナリを 2 回回した差と同程度**。`ypls`・`roUz` 等は数値的にゼロで判定対象外。
なお**既定のブロックサイズ (512) では SLAU の node カーネルが起動できず** (`convectiveFlux_d.cu:316`,
`too many resources requested for launch`)、`FORGE_CUDA_BLOCKSIZE=64` を両バイナリ共通で指定した。
本変更とは無関係 (比較対象の基準バイナリでも同じ失敗が出る) が、~~ローカル GPU (RTX 3060) の制約~~
という当初の帰属は誤り (2026-09-24 撤回)。レジスタ 65,536/ブロックは **sm_70 以降どの GPU でも同じ定数**で、GPU の大小とは無関係である。
**2026-09-24 に単一 TU を各 commit でコンパイルして `cuobjdump -res-usage` の `REG` を測り、起点を特定した**
(nvcc 12.0 / sm_86、`-O3 --generate-code=arch=compute_86,code=[compute_86,sm_86]`):

| commit | `SLAU_d` REG | 512 起動 (R×512 ≤ 65536) |
| --- | --- | --- |
| `bcc68797` (2026-09-20) | 114 | 可 |
| `e7e4a0c1` | 114 | 可 |
| `f8fca224` / `a0426086` | 113 | 可 |
| `6430909d` (reconT) | **128** | 可 (65,536 = ちょうど上限) |
| `571e81df` (limiter_T) | **138** | **不可** ← 最初に越えた commit |
| `adcd5579` / `39526328` / HEAD | 136 | 不可 (上限 481 threads) |

**`6430909d` が余裕を使い切り `571e81df` が越えた** (手順と対策の正本は
[`procedures/development-environment.md`](../../procedures/development-environment.md) の
「カーネル起動が `too many resources requested for launch` で落ちるとき」)。どちらも既定 0 でビット不変の opt-in だが、
**レジスタは実行時フラグに依らず確保される**ので既定構成の起動が壊れた。
**本計画の `slauWallNormalChi` (`5a4886ea`) はレジスタを 1 つも増やしていない** (直前の `39526328` が既に 136)。
恒久対策 (既定ブロック 256 化 / `__launch_bounds__`) は**本計画の範囲外**。node の残差 gather は `atomicAdd` なので
ブロックサイズ変更はビット不変でない (別途判断が要る)。

**V3 case/48 (flag 1 の実害) — 合格**。同じ収束場から flag 0 / flag 1 を 20000 step
(`_v3/v3_flag0`, `_v3/v3_flag1`、新バイナリ・`slauWallNormalChi` のみ差分):

| 量 | 最大差 | 許容 |
| --- | --- | --- |
| $C_f$ (5 station) | ~~**0.000 %**~~ **0.0018 %** | ±1 % |
| $q_w$ (5 station) | ~~**0.011 %**~~ **0.0021 %** | ±1 % |

**訂正 (2026-09-25, codex result-3 m6, #17)**: 旧値 0.000 % / 0.011 % は評価表の**表示値 ($C_f\times10^3$ 小数 3 桁、$q_w$ 小数 2 桁) どうしの差**で、丸めの見かけだった。**正本は丸める前の値**: `_v3/v3_flag{0,1}/res_wall_4_20000.h5` の `twall_x`・`qwall` を 5 station の最寄り節点で比べた最大差 **0.0018 % (x=0.75) / 0.0021 % (x=0.90)**。codex の再計算 0.001804 % / 0.002388 % と同じ桁 (抽出位置の違い)。許容内という結論は変わらない。

§4 の予測どおり: **付着境界層では壁法線方向の $\Delta p \approx 0$ なので $\chi_n$ を上げても圧力差項は ≈0 のまま**。
残差水準も同等で、**`mSLAU` の収束悪化 (§3.1) は局所適用では現れていない**:

```
check_convergence.py: 両者とも NOT CONVERGED (stalled/plateau)
  ※ 収束済みの場からの再スタートなので残差が床から始まり低下桁数が出ない (両者同条件)。
     見るべきは水準の一致:
  rms_ro    flag0 1.94e-07 / flag1 1.83e-07      rms_roUx  flag0 2.71e-04 / flag1 2.56e-04
  rms_roe   flag0 2.43e-01 / flag1 2.29e-01      rms_roOmega flag0 2.12 / flag1 1.98
```

注意: flag 1 は `rms_roOmega` の初期値が 5.46e+02 (flag 0 は 2.85) と跳ねてから 2.4 桁下げて同水準に落ち着く。
壁隣接面の流束変更が $\omega$ に初期過渡を与えている。発散はしないが、SST を含む段階起動では留意する。

**V1 (case/46) — 進行中**。`run_0435_3d_junction_chi_n/`。`run_0430/res_1800` 起点、run_0432/0433 と同一設定
(差分は `slauWallNormalChi: 1` と step 数)。**初回投入はディスク不足で step ~2200 で書き込み失敗** (`outStepInterval: 200` で
31 dump = 10 GB に対し空き 9.8 GB) したが、**その時点までのプローブで排出の停止は確認できた**:

| 節点 | 起点 $\rho/\rho_{Min}$ | $10\rho_{Min}$ 到達 | ピーク | step ~2000 |
| --- | --- | --- | --- | --- |
| 153797 (側壁後端面) | 2.07 | **step 13** | 31.5 | 22.1 |
| 153880 (隣の壁節点) | 2.71 | **step 11** | 32.4 | 23.8 |
| 189814 (カウル後端面) | 5.88 | **step 8** | 33.7 | 15.5 |

SLAU (run_0432) では同じ節点が 1.71e-4 から減衰して床を割り step 2615 で NaN。
**期限 (V1-e, step 200) に対し 15〜25 倍の余裕**で到達し、**step 100 付近で頭打ち** (= §4.2 の「ある高さで止まる」と整合)。
内点の $|U|$ は 1316 → 1345 → 1310 と**維持**されており、駆動条件が消えたのではなく補充経路が回復したという読みと整合。
過補正 ($P_w>2P_i$) も無し (ピーク比 ≈0.33、近似式の予測 0.31〜0.38 の範囲)。
`outStepInterval: 1000` に変更して 6000 step を回し直し中。**正式判定 (V1-b/c/d) は完走後**。

#### 効果の分離 — 2×2 が揃った (2026-09-23, codex result M2 の採用結果)

codex result M2「**`outflow` + flag 0 が欠けており独立性は証明できない**」を採用し、
`run_0440_3d_junction_outflow_flag0` (run_0437 と**設定完全一致で flag だけ 0**、diff で確認、同一起点 `run_0430/res_1800`) を回した。

| | `outlet_statPress` | **`outflow`** |
| --- | --- | --- |
| **flag 0** | `run_0432`: **step 815 で NaN** | `run_0440`: **step 815 で NaN** |
| **flag 1** | `run_0435`: 6000 完走だが出口が破綻 (残差 +0.8 桁) | `run_0437–0439`: **通算 66000 完走・`ALL STEADY`** |

**壁節点の排出の軌跡が 5 桁目まで一致する** (CV 153797 の $\rho = P/(RT_w)$):

| step | run_0432 (statPress+flag0) | run_0440 (outflow+flag0) |
| --- | --- | --- |
| 1 | 1.7051e-04 | 1.7051e-04 |
| 200 | 1.1123e-04 | 1.1121e-04 |
| 400 | 8.1761e-05 | 8.1760e-05 |
| 600 | 5.1428e-05 | 5.1435e-05 |
| 800 | 2.6535e-05 | 2.6542e-05 |
| NaN | step 815 | step 815 |

→ **言えるのはここまで (codex result-2 m6)**: 「**2×2 で、step ≤ 800 の $\rho_w$ 軌跡が ~~5 桁一致し~~ 0.0264 % 以内で一致し (step 800: 2.6535e-5 対 2.6542e-5。5 有効桁の一致ではない; codex result-3 m6 で訂正 2026-09-25)、
NaN の step も同じだった。この区間では出口 BC の影響が検出されなかった**」。
**「独立」「確定」とは書かない**。理由: (i) 独立性は応答面の性質で、結果の一致は
「その応答量・その水準・その時間窓で交互作用が検出されなかった」までしか言えない。
(ii) 応答量の一つ (NaN step) は**同一入力でもばらつく非再現指標**なので一致は偶然に弱い。
(iii) 出口劣化は数千 step かけて上流へ育ち、排出は step ≤ 800 の局所現象 =
**設計上、交互作用が小さいはずの窓だけを見ている** (run_0435 と run_0437 が 6000 step で違う場になったのは、
窓の外では BC が効いた証拠)。
運用上の結論は変わらない: **出口だけ直しても step 815 で落ち、$\chi_n$ だけでも出口が壊れる。両方入れる。**
以前に根拠としていた「step 800 対 1000 の亜音速率 7.1 % 対 9.0 %」(同一時点でない 1 点比較) は**この直接比較に置き換える**。

#### V1 の再測定 (2026-09-23, 出口修正後)

**V1 の検証ケースに独立の欠陥が見つかった**: 出口 `outlet_statPress` (Ps 2851) が node の壁列・後流の亜音速ノードに
実圧と桁違いの背圧を課し、**出口から上流へ圧力が育つ**
([`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 **B1f**)。**本計画の $\chi_n$ とは独立**で、
step を揃えた比較 (対策なし step 800 / 対策あり step 1000) では出口の亜音速率が 7.1 % 対 9.0 % とほぼ同じ =
**対策は出口の劣化を加速していない**。

出口を `outflow` にした **`run_0437_3d_junction_chi_n_outflow`** (同一起点・同一設定、出口と `far_bottom` のみ変更) で測り直した:

| 指標 | `outlet_statPress` (run_0435) | **`outflow` (run_0437)** |
| --- | --- | --- |
| 出口の亜音速率 | 91.2 % | **4.4 %** (起点 3.5 %) |
| 出口の逆流 | 19.8 % | **0.1 %** |
| 全域で P>2e5 の節点 | 81,920 | **0** |
| `rms_ro` の変化 | **+0.8 桁 (悪化)** | **−1.3 桁** |
| 全域 $\rho$ 最小 | 9.33e-4 | 9.27e-4 (床の 9.3 倍) |

**V1-a / V1-e 合格** (判定区間 run_0435 step 0–2000、出口の汚染が後端面域 x/H 1.2–1.6 に届く前):
監視 CV は **step 8–13 で $10\rho_{Min}$ を超え**、期限 200 step に対し 15〜25 倍の余裕。床到達 0、過補正なし。

~~**V1-c 合格** (run_0437 `res_6000`)。予測 419.3/459.1/322.5 Pa に対し観測 565.9/604.7/449.3 Pa、観測/予測 1.317–1.393 で PASS。~~
**撤回 (2026-09-23, codex result M1)**: この判定に使った `diag_wall_cv_budget.py` の `slau_mdot` は
**変更前の $\chi$ (速度の大きさ基準) のまま**で、実装の `chi_mass` (`convectiveFlux_slau_d.inc.cuh:552`) を反映していなかった。
**検証したい補充項そのものが診断から抜けていた**ため、上の一致は意味を持たない。
codex の反例: 壁側 0 / 内点の接線 1000 m/s / $\hat c$ 700 / 面法線速度 0 / $\Delta p$ 900 Pa / $A$=1 m² で
**修正前ツール $\dot m = 0$、実装 $\dot m = -0.642857$ kg/s**。

**V1-b / V1-c 合格 (再測定、2026-09-23)** — ツールに `--wall-normal-chi` (カーネルと同じ壁ゲート + 面法線 $\widehat M_n$) を実装し、
**単体試験 `cad/test_diag_wall_cv_budget.py` でカーネル式と項ごとに一致することを固定**してから、
**`run_0439` (`ALL STEADY`) の最後の 3 dump** で判定した (plan §6 の「STEADY 後の最後の 3 dump」を満たす):

| CV | 予測 $P_w$ | 観測 $P_w$ | 観測/予測 | 判定 |
| --- | --- | --- | --- | --- |
| 153797 | 512.84 Pa | 512.84 Pa | **1.000** | PASS |
| 153880 | 512.70 Pa | 512.70 Pa | **1.000** | PASS |
| 189814 | 553.61 Pa | 553.61 Pa | **1.000** | PASS |

**V1-b も同時に合格**。最後の 3 dump の全接続面の正味流束 $\Sigma\dot m$:

| CV | res_24000 | res_30000 | res_36000 |
| --- | --- | --- | --- |
| 153797 | — | 3.02e-14 | 5.21e-14 |
| 153880 | 2.07e-14 | −5.39e-16 | 3.32e-14 |
| 189814 | −3.94e-13 | −8.03e-14 | −2.08e-14 |

起点の +1.49e-9 kg/s から **5 桁小さく**、$\lvert\Sigma\dot m\rvert\Delta t/V$ は要求の「$\rho_w$ の 1 %/dump 未満」を桁違いに下回る。
ツールが実装と同じ流束を計算するようにしたところ、予測と観測が 6 桁目まで一致した (上表)。**平衡壁圧の式はこの 3 CV・3 dump では機序を再現している**。

**V1-d 合格 — `run_0439_3d_junction_outflow_steady`** (run_0438 から 36000 step 継続、**通算 66000 step**):

```
check_quasisteady.py --series-csv ... --series-cols ro_153797,ro_153880,ro_189814,ro_min,n_floor
  ro_min   tail mean=0.0005005  drift=0.0%/tail  fluct=0.0%  [漸近値 0.000500518, 最終比 +0.000 %]  STEADY
  n_floor  tail mean=0           drift=0.0%/tail  fluct=0.0%                                        STEADY
  shock    tail mean=1.076       drift=0.0%  fluct=0.0%                                             STEADY
  machmax  tail mean=6.819       drift=0.0%  fluct=0.0%                                             STEADY
  pmax     tail mean=1.078e+05   drift=0.0%  fluct=0.0%                                             STEADY
OVERALL: ALL STEADY
```

壁節点の密度は 5 桁目まで静止:

| step | $\rho$(153797) | $\rho$(189814) | $\rho_{min}$ (全域) | 床到達数 |
| --- | --- | --- | --- | --- |
| 12000 | 1.5069e-03 | 1.6217e-03 | 4.9933e-04 | **0** |
| 24000 | 1.5068e-03 | 1.6265e-03 | 5.0049e-04 | **0** |
| 36000 | 1.5068e-03 | 1.6268e-03 | 5.0052e-04 | **0** |

`check_convergence` は最後の 36000 step 区間で `NOT CONVERGED (stalled/plateau)` (低下 0.6–1.2 桁) だが、
**絶対水準は `rms_ro` 2.40e-10** = 通算 **4.40 桁**低下 (6.02e-06 → 2.40e-10)。
**ここから「全方程式が収束した」とは言わない** (codex result M4)。`rms_ro` の絶対値だけでノイズ床を断定する根拠は無い。
**V1-d は「要求した密度系列が STEADY かつ全壁の床到達 0」という限定で合格**とし、方程式全体の収束とは分ける。
`pmax` 1.078e5 は入口の排気圧 101 kPa をわずかに上回る程度で、**出口起因の圧力の壁は消えている**。

#### V1-d と V3 (SERN 2D) の結果 (2026-09-23)

**V1-d — `run_0438_3d_junction_outflow_long`** (`run_0437/res_6000` から継続、通算 30000 step、出口 `outflow`):

```
check_convergence.py -> NOT CONVERGED (**still converging — run more steps**)
  rms_ro      2.65e-07 -> 1.25e-09  (2.1 dec) falling
  rms_roUx    4.87e-04 -> 2.13e-06  (2.1 dec) falling
  rms_roUy    1.38e-04 -> 4.12e-07  (2.2 dec) falling
  rms_roUz    8.54e-05 -> 6.81e-07  (1.9 dec) falling
  rms_roe     5.29e-01 -> 2.66e-03  (1.9 dec) falling
  rms_roY0/1  1.7-2.0 dec falling
```

**全列 `falling`** で、判定理由が `stalled/plateau` (頭打ち = スキーム変更が要る) から
**`still converging` (step を増やせ)** に変わった。通算では起点から **4.1 桁**低下
(run_0437 の 1.3 桁 + run_0438 の 2.1 桁)。**STEADY を名乗れる状態ではない**が、
「プラトーで止まっている」でもない。V1-d は**「さらに伸ばせば下がる」ところまで**とし、
STEADY の確認は残作業 (§5.1) に置く。

**V3 SERN 2D 生産 — PASS**。`problem_moo_frozen_tp_cycle3op.yaml` (設計点 m6_on)、
**両方とも新既定の `outflow`** で段階起動 (暖機 → soft 4000 → mid 4000 → 本段 12000) を rc=0 で完走
(`_v3sern/v3s_flag0`, `_v3sern/v3s_flag1`、差分は `slauWallNormalChi` のみ):

| 量 | flag 0 | flag 1 | 差 | 許容 |
| --- | --- | --- | --- | --- |
| **$C_T$** | 0.9293864 | 0.9291175 | **−0.029 %** | ±0.1 % → **PASS** |
| **$C_T$ (摩擦込)** | 0.8983631 | 0.8981576 | **−0.023 %** | ±0.1 % → **PASS** |
| $C_T$ 摩擦成分 | −0.0310233 | −0.0309599 | −0.205 % | (判定対象外) |
| $C_L$ | 0.2771498 | 0.2784627 | **+0.474 %** | (判定対象外) |
| $C_M$ | −7.1572 | −7.1823 | **+0.350 %** | (判定対象外) |
| `sep_frac_ramp` | 0.0 | 0.0 | 0 | — |

両者の `check_convergence` は `NOT CONVERGED (stalled/plateau)` だが**これは本ケースの性格**
(case/46 README run_0091: 「verdict は NOT CONVERGED stalled/plateau = 本ケースの性格」)。残差の水準は一致
(`rms_ro` fin 8.82e-05 対 8.09e-05、`rms_roOmega` 1.24e+02 対 1.02e+02) で、**flag 1 の方がわずかに低い**。

**$C_L$ / $C_M$ の差は run 間ノイズではなく実効果** (**2026-09-25 限定: flag0 が `DRIFTING` で定常差として未確定、#14**) (2026-09-23 に各側 3 本の反復で切り分け、
`_v3sern/v3s_flag{0,1}[_b,_c]`)。同一側の最大ばらつき (床) と flag 間の平均差:

| 量 | flag 0 (3 本平均±σ) | flag 1 (3 本平均±σ) | ノイズ床 | 差 | 床との比 |
| --- | --- | --- | --- | --- | --- |
| $C_T$ | 0.929415±0.000021 | 0.929125±0.000017 | 0.0049 % | **0.031 %** | 6.4 倍 |
| $C_T$ (摩擦込) | 0.898393±0.000021 | 0.898165±0.000020 | 0.0052 % | 0.025 % | 4.9 倍 |
| **$C_L$** | 0.277277±0.000090 | 0.278478±0.000060 | 0.0709 % | **0.433 %** | **6.1 倍** |
| **$C_M$** | −7.159973±0.001952 | −7.182615±0.001359 | 0.0591 % | **0.316 %** | **5.4 倍** |
| $C_T$ 摩擦成分 | −0.031022±0.000001 | −0.030960±0.000003 | 0.0196 % | 0.203 % | 10.3 倍 |

**全量が床の 5〜10 倍で有意に動いている**。V3 の合否 ($C_T$ ±0.1 %) は変わらないが、
**「ノイズかもしれない」という逃げ道は無い**。$\chi_n$ は確かに解を動かしており、$C_L$ で 0.43 %、$C_M$ で 0.32 %。

**評価が割れる点 (既定化の判断に直結)**:
- 動かしている面は run_0437 では全内部面の 0.53 % だった (ただし**これは 3D 接続模型の値で 2D には転用できない**。
  §6.2 の撤回参照)。$\chi=0$ だった面では**圧力–質量結合が欠落していた**ので**動いた方が物理的に正しい可能性がある**が、
  **どの面がそうだったかは 2D では測っていない**。
- 一方 0.43 % は設計最適化 (MOO) では無視できない量。既定化するなら「どちらが正しいか」の独立した根拠が要る
  (格子収束・別スキーム・実験/文献との比較のいずれか)。
- **本計画では既定を変えない** (§2 スコープどおり opt-in のまま)。この差の正否判定は既定化の検討時に別途行う。

**副産物**: この 2 本は**出口の既定変更 (B1f) が生産レシピを壊していない**ことの実証でもある
(段階起動 4 段を含めて rc=0、力係数が過去の生産値と同オーダー)。

#### #10b の再判定 — **手順と合格式を計算前に登録** (2026-09-25)

**2026-09-24 の時系列判定は operand の単位を取り違えていた** (許容差と 1/5 則は正しい)。
`check_quasisteady.py:280-282` の `fluct = span/|mean|` は**系列自身の平均に対する相対値**であり、
§6 V3 表 (:251-257) の `--drift 0.0002 --osc 0.0005` は**行の「量」= run の量** ($C_T$・壁圧分布) に
掛ける前提の数である。それを**差の系列** (`dFx_pct` 等、平均 0.13 %) に渡したので、
実効的に $F_x$ の $6.5\times10^{-5}$ % を要求していた — **osc で約 770 倍、drift で約 300 倍厳しい**
(当初「2000 倍」と書いたのは誤り)。**許容差 (0.1 % / 1 %) と 1/5 則は動かさない。直すのは当てる相手だけ。**

**本段を 36000 step に延長しても差は収まらなかった** (`_v3sern/v3s_flag0_ext` / `v3s_flag1_ext`、
`restart_field.py` で `res_12000` から継続、`VERDICT: OK`、500 step ごと 72 dump):
`foot_L2_pct` 変動 1.3 → **2.2 %**、`dFx_pct` 53.5 → **72.5 %**、`dFy_pct` 25.8 → **32.7 %**
(drift は 0.2→0.1 / 9.0→4.1 / 4.4→1.8 % と減少、末尾平均は 0.8385 / 0.1306 / 0.1931 で安定)。
**drift は既にほぼ 0 なので、さらに延長しても振幅は縮まない。**

**決め手の量 $C_T$ の時系列は一度も判定していない** (6d は終端 1 点の反復 3 本)。**次の手順で判定する**:

1. 各 dump に `design/forge_design/metrics/sern_momentum.py` の `check_run(run_dir, step=...)` を掛け、
   **run ごとに** $C_T$ / $C_T$(摩擦込) / $C_L$ / $C_M$ の時系列を作る。
2. **各 run を §6 の閾値そのもの** (`--drift 0.0002 --osc 0.0005`, `--tail 0.4`) で判定する。
3. 差は**末尾平均の差**、不確かさは**両 run の span/2 の大きい方** (振幅) とする。
4. **合格式 (計算前に固定)**: $\lvert\Delta(\text{末尾平均})\rvert + \text{振幅} \le 0.1\,\%$ ($C_T$ の許容差の単位で)。

- **A**: 両 run が STEADY、または `OSCILLATING` でも上式を満たす → **#10b PASS**
  (`foot_L2_pct` は 0.84 ± 0.01 % を併記)。
- **B**: `DRIFTING`、または上式を満たさない → **`flag` の差ではなく本ケースの限界サイクル振幅が
  受入判定の分解能を超えている**。「**判定不能 (限界サイクル振幅 > 1/5 則)**」と書き、
  振幅を下げる手 (`cfl_pseudo`・平均化窓) は §4 の設計判断として別途諮る。**延長では解けない。**

**撤回**: ~~「各 run 単体は `ALL STEADY` なので #10b を満たす」~~ — これは**既定閾値 (drift 5 % / osc 10 %) で
判定したもの**で、§6 の 0.02 % / 0.05 % ではない。$F_x$ の変動 0.1 % は §6 の osc を超えるので、
§6 の閾値では `OSCILLATING` である (2026-09-25 撤回。compare-like-with-like の型)。

**6d の「床の 5〜10 倍で実効果」は再判定が要る** (:814-826)。6d の床は**同一初期場・同一 step の反復 3 本の
散らばり**で、限界サイクルの**位相が揃う**ため**擬似時間内の振幅を含まない**。一方 `dFx` の変動が各 run の
変動とほぼ等しいことは、**flag0 と flag1 の限界サイクルが位相同期していない**ことを意味する (同期していれば
差は定常になる)。よって flag 間の 1 点差は最大で両振幅の和を含みうる。表の「ノイズ床」列は
**「反復床 (同 step)」**と改名し、**時系列振幅の列を足す**。$C_L$ 0.433 % が「実効果」かは同じ $C_L$ 系列で決める。

**#10b 結果 — PASS** (2026-09-25)。登録した手順と合格式の機械的な適用。
run: `_v3sern/v3s_flag0_ext` / `v3s_flag1_ext` (`res_12000` から `restart_field.py` で 36000 step、72 dump)。
力係数は**本番と同じ関数** (`forge_design.metrics.sern_forces.force_history`、`runner_sern.py:885` と同一引数・
`twall_on_fluid` の cell/node 規約込み) で作った (`v3sern_ct.py`)。**自前で係数を組まない**。

**各 run を §6 の閾値** (`--drift 0.0002 --osc 0.0005 --tail 0.4`) **で判定**:
$C_T$ と $C_T$(摩擦込) は**両 run とも `STEADY`** (drift 0.0 % / fluct 0.0 %)。
$C_L$ は flag0 `DRIFTING` (fluct 0.2 %) / flag1 `OSCILLATING` (0.2785 ± 0.00016)、
$C_M$ も同様 (−7.183 ± 0.0036) だが、**§6 V3 表のゲート量は $C_T$ のみ**である。

| 量 | flag0 末尾平均 | flag1 末尾平均 | Δ [%] | 振幅 [%] | $\lvert\Delta\rvert+\text{amp}$ | 許容 | 判定 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **$C_T$** | 0.9294199 | 0.9291329 | −0.0309 | 0.0073 | **0.0382** | 0.1 | **PASS** |
| $C_T$ (摩擦込) | 0.8983967 | 0.8981738 | −0.0248 | 0.0076 | 0.0324 | 0.1 | PASS |
| $C_L$ | 0.2772934 | 0.2785043 | +0.4367 | 0.0853 | 0.5220 | (対象外) | — |
| $C_M$ | −7.1602977 | −7.1832303 | −0.3203 | 0.0750 | 0.3953 | (対象外) | — |

`foot_L2_pct` は **0.8385 ± 0.0091 %** (振幅 0.009 pp、1/5 則 0.2 pp 以内)、平均+振幅 0.848 % ≤ 1 % で許容内。
`dx_foot` は 0 で `STEADY`。

**6d の床を「反復床 (同 step)」から「限界サイクル振幅」へ置き換えた再判定** (結論は変わらない):

| 量 | 差 | 反復床 (同 step) | 旧比 | **限界サイクル振幅** | **新比** |
| --- | --- | --- | --- | --- | --- |
| $C_L$ | 0.437 % | 0.0709 % | 6.1 倍 | **0.0853 %** | **5.1 倍** |
| $C_M$ | 0.320 % | 0.0591 % | 5.4 倍 | **0.0750 %** | **4.3 倍** |

反復床は**同一初期場・同一 step の反復**なので限界サイクルの**位相が揃い、擬似時間内の振幅を含まない**。
`dFx` の変動が各 run の変動とほぼ等しかったことは **flag0 と flag1 の限界サイクルが位相同期していない**
ことを意味する。正しい床は振幅側で、**比は 6.1→5.1 / 5.4→4.3 に下がるが~~「実効果」の結論は変わらない~~**。
**撤回 (2026-09-25, codex result-3 M3, #14)**: flag0 の $C_L$ は `DRIFTING` なので定常差として確定できない。言えるのは「**現区間の末尾平均差 $C_L$ 0.437 % / $C_M$ 0.320 %、flag0 は `DRIFTING` で定常差として未確定**」まで。準定常確認 (両 run STEADY) は #11 の前提。

**この項目で自分が犯した誤り (4 件。いずれも今日ほかの場面で踏んだ型の再発)**:

1. 「延長すれば収まる」という見込み → **外れ**。drift は減ったが fluct は 53.5→72.5 % と増えた。
   **drift がほぼ 0 の振動は延長では縮まない**。
2. 「意図の 2000 倍厳しい」→ **osc 約 770 倍 / drift 約 300 倍**の誤り。0.02 (drift) と 0.05 (osc) も混同した。
3. 「各 run は `ALL STEADY` だから #10b を満たす」→ それは**既定閾値 (5 %/10 %) での判定**で、
   §6 閾値 (0.02 %/0.05 %) ではない。**比べる相手を揃えていない** ([[compare-like-with-like]] の型)。撤回済み。
4. 力係数を**自前で組もうとした** → `sern_forces` が本番にあり、`twall_on_fluid` の cell/node 規約を
   落とすところだった ([[compare-like-with-like]] の 1 番目「ツールが実装と別のものを計算していた」の型)。

#### M3 の不足ゲート — 進捗 (2026-09-23)

**衝撃足の壁圧 — PASS**。ノルムは**比較する前に**固定した (§6 V3 の行。結果を見てから条件を作らない):
flag 0 の run で $\max\lvert\mathrm{d}p_w/\mathrm{d}x\rvert$ の位置を $x_{foot}$ とし、
ランプ上 $x \in [x_{foot}-2t,\ x_{foot}+5t]$ ($t$ = カウル板厚 2 mm) で比較。
`_v3sern/v3s_flag{0,1}` の最終壁 dump (`res_ramp_4_12000.h5`, 328 節点、帯内 34 節点):

| 判定項目 | 実測 | 合格条件 | 余裕 |
| --- | --- | --- | --- |
| 壁圧 $p_w$ の相対 L2 差 | **0.842 %** | ≤ 1 % | 1.19 倍 |
| $x_{foot}$ の差 | **0.00 格子間隔** | ≤ 1 格子間隔 | — |

$x_{foot}$ = −0.878 mm (ランプ上の格子間隔 中央値 1.081 mm)。位置は完全一致だが**壁圧の分布は 0.84 % 動いており、
余裕は 1.19 倍しかない**。$C_L$ の 0.433 % と整合する。**「影響は測定限界以下」ではない**。

**`cfl_pseudo` 0.4 — 「同じ累積 CFL での終端場差」まで (固定点 PASS は保留、codex result-2 M2)**。`run_0439` (cfl 0.2 × 66000 step, `ALL STEADY`) と
`run_0442_3d_junction_cfl04_long` (**cfl 0.4 × 33000 step = 同じ累積 CFL**) の比較:

| 量 | cfl 0.2 × 66000 | cfl 0.4 × 33000 | 相対差 |
| --- | --- | --- | --- |
| $\rho$(153797) | 1.506828e-03 | 1.507120e-03 | **0.019 %** |
| $\rho$(153880) | 1.506408e-03 | 1.506669e-03 | **0.017 %** |
| $\rho$(189814) | 1.626776e-03 | 1.627565e-03 | **0.049 %** |
| $\rho$ 全域 (相対 L2) | — | — | **0.011 %** |
| $P$ 全域 (相対 L2) | — | — | **0.0033 %** |

CFL を 2 倍にしても全域の相対 L2 差は **0.011 %** (~~「0.01 % 以内」~~は有利な側への丸めだった)。
**ただしこれは固定点の一致の証明ではない** (codex result-2 M2): `run_0442` の収束・準定常 VERDICT を出しておらず、
`run_0439` も `check_convergence` は `NOT CONVERGED`、密度系列に許した閾値は既定 drift 5 % / osc 10 % で
**0.011 % 規模の比較を支える精度ではない**。同じ累積 CFL で近い場になることは**同じ過渡を追っている場合にも起こる**。
固定点の判定には、両 CFL で比較精度より十分小さい時間変化を確認し、**全保存量の収束判定 + 独立に延長した末尾区間の比較**が要る (未実施)。

**最初の比較は不成立だった (記録として残す)**: `run_0441` (cfl 0.4, 6000 step) と `run_0437` (cfl 0.2, 6000 step) を
並べて 8〜16 % の差を得たが、**両者とも `still converging`** だったので、差は「固定点が違う」のではなく
**収束の進み具合の違い**だった。CFL が 2 倍なら同じ step 数での到達点が違うのは当然で、比較になっていない。
**累積 CFL を揃えたら 0.01 % 以内に収まった** = 読みが正しかったことになる。
**「比較するものが揃っているか」の確認不足**は本 plan で 3 回出ている (codex M1 のツール、M2 の 1 点比較、ここ)。

~~**未実施**: 格子感度 (3 水準)、case/16 の V2・V3 (run データが AWS・手元とも無く**メッシュ生成から**要る)。~~ **決着 (2026-09-25)**: 3 水準は済 (2026-09-24, `run_0445`/`0446`、下の V3 格子感度節。限定 3 つ付き)、case/16 は #10 で V2・V3 PASS (§6 V3 case/16 条件)。

#### V3 格子感度 (2026-09-23, codex result M3 の残り)

接続模型の格子を `cad/hex_junction_model.py --scale 0.71` で生成し (**1,108,434 節点**、採用格子 2,446,572 の 45 %、
`check_mesh_quality` **PASS** AR ≤ 1000 / skew ≤ 0.90)、**同一設定・出口 `outflow`・差分は flag のみ**で対にした
(`run_0443_3d_gs071_flag0` / `run_0444_3d_gs071_flag1`、暖機 1 次 CFL 0.2、6000 step)。

| 格子 | 節点数 | flag 0 | flag 1 |
| --- | --- | --- | --- |
| **粗 (新, 2026-09-24)** | **494,037** | **`VERDICT: DIVERGED (NaN/Inf)`** (step 1908) | **6000 完走、NaN 0** |
| 中 (scale 0.71) | 1,108,434 | **step 3031 で NaN** | **6000 完走 (rc=0)** |
| 採用 | 2,446,572 | step 815 で NaN | 完走 (通算 66000 で `ALL STEADY`) |

**3 水準で同じ成否が出た** (`run_0445_3d_gs050_flag0` / `run_0446_3d_gs050_flag1`)。`run_0446` の全残差列は
`rms_ro` −2.1 / `rms_roUx` −2.0 / `rms_roUy` −2.5 / `rms_roUz` −2.3 / `rms_roe` −2.3 / `rms_roY0` −2.5 桁 (falling)、
`rms_roY1` のみ −1.6 桁 plateau。**収束はしていない** (暖機 6000 step の段なので中・細の対と同じ扱い)。
最終 `res_6000.h5` は非有限 0、$\rho$ 3.92e-4〜0.322、$P$ 131〜1.10e5 Pa、$T$ 207〜2519 K。

**この比較の限定 (3 つ。外さないこと)**:

1. **一様な細分列ではない**。粗格子は流れ方向の刻み (`HX`) を採用格子と同じに据え置き、**壁法線 (`H1` 3.2e-4,
   `G` 1.44) と spanwise (`NZ` 16) だけ**を粗くした。`HX` を 2 倍にすると長さ 0.0574 の区間が 1 セルに潰れ、
   x 隣接比が 1.886 (許容 $G\times1.02=1.469$) で**メッシャの受入検査に落ちる** — 形状が決める区間なので、
   これは「粗ければ生成を失敗させる」という設計どおりの拒否である。
2. **中格子 `run_0443`/`0444` を作ったコマンドが記録に残っていない** (AWS 上で手打ち)。粗格子は
   `run_junction_model.py` の既定 (`problem_3d_prod_m6on_g1.yaml`) で組んだので、**設定が完全一致する保証はない**。
3. **ブロックサイズが揃っていない**。粗格子は `FORGE_CUDA_BLOCKSIZE=128` (現 HEAD の `SLAU_d` は `REG:136` で
   既定 512 が起動できない)、中・細格子は当時の既定 512。**面流束は 128 と 256 でビット同一** (§6.2 V5 P5) なので
   定性的な成否には効かないが、**NaN step の数値は比較対象にしない** (もともと非再現指標)。

**NaN step は 1908 / 3031 / 815 で単調でない。** 上の 3 つの限定があるので、この数から機序を組み立てない
(「細かいほど早い」という以前の説明は codex result-2 M3 で撤回済み)。**見たのは成否だけ**である。

**言えるのはここまで (codex result-2 M3)**: 「**試した 2 格子では同じ成否。採用格子の方が少ない反復で発散した**」。
**格子非依存の証明ではない** — 粗い側の flag 1 は 6000 step までで、**定常解・壁圧・収束率の格子比較をしていない**。
`mSLAU` の格子感度の懸念が解消した証拠にもならない。

~~**定量的には「細かいほど早く顕在化する」**~~ (**撤回の残存を処置 2026-09-25**, codex result-3 M4: NaN step は非再現指標で 3 水準では単調でもない。下の体積比は記録として残す)。側壁後端面の壁ノードの双対体積 (中央値):

| | 壁ノード数 | 体積の中央値 | 比 |
| --- | --- | --- | --- |
| 中 | 2,640 | 1.1731e-11 m³ | — |
| 採用 | 5,208 | 5.9322e-12 m³ | **1.98** |

~~NaN までの step 比 3.72 は体積比 1.98 と局所時間刻みで説明できる~~
**撤回 (2026-09-24, codex result-2 M3)**: 局所 $\Delta t$ を**実測していない**。更新は $\Sigma\dot m\,\Delta t/V$ で決まり、
格子を変えると面流束・時間刻み・初期補間誤差も変わるので**体積だけでは説明できない**。
加えて **NaN step は同一入力でもばらつく非再現指標**なので ([[base-wake-resolution-rule]])、比を取る指標として不適切。
機序を述べるなら対応 CV の面収支・$\Delta t/V$・初期状態の実測が要る (未実施)。

~~**含意 (plan の記述として重要)**: これまでの §4 の機序説明は格子依存性に触れていなかったが、
**壁解像度を上げるほどこの問題は早く顕在化する**。y⁺≈1 を要求する生産計算 (冷却壁など) ほど避けられない。
逆に粗い格子しか使わないなら顕在化しない可能性がある。~~ (**撤回 2026-09-25**, codex result-3 M4: 上と同じ。成否以外は言えない)**「どの解像度から要るか」は本 plan では決めていない**。

~~**未実施**: 3 水準目 (粗).~~ **決着 (2026-09-24)**: 粗格子は `--scale` を使わず壁法線と spanwise だけを粗くして作った (上の表・限定 1)。以下は当時の記録: `--scale` を粗い側に使うと**成長率が $G^{1/s}$ で跳ね上がり** (1.2 → 1.44)、
x 分割が減って**区間の継ぎ目で間隔比が 1.7–1.9 に跳ね**、生成器の受入検査 (`x_adjacent_ratio_max <= G*1.02`) に落ちる。
`--set` で直接指定しても `--scale` が**後から上書きする**ので併用できない。`HX` を 0.10 → 0.07 で x 比 1.71 → 1.29 まで改善したが未達。

#### M3 の現状まとめ (2026-09-23 時点、**2026-09-25 に状態欄を更新**)

| 項目 | 状態 |
| --- | --- |
| V0 面反転・等状態・非対象面 | **済** (`cad/test_diag_wall_cv_budget.py` が `ALL PASS`) |
| V1-b (収支の符号) | **済・合格** (STEADY 後 3 dump で $\Sigma\dot m$ が 5 桁縮小) |
| 衝撃足の壁圧ノルム | **済・PASS** (事前固定、L2 差 0.842 %、位置ずれ 0) |
| `cfl_pseudo` 0.4 の固定点 | **済・PASS** (累積 CFL を揃えて 0.01 % 一致) |
| 格子感度 | ~~**2 水準 済** (定性的に格子によらない、細かいほど早い)。**3 水準目は未**~~ **3 水準 済 (2026-09-24)**: 成否が同じ。**格子非依存の証明ではない**・「細かいほど早い」は撤回 (V3 格子感度節の限定) |
| case/16 の V2・V3 | **済 2026-09-25 (#10, V2・V3 PASS)**。以下は当時の記録: ~~**未**。~~ node run が**全部 Euler (`slip` 壁)** で $\chi_n$ の対象外。**NS 設定を新規に組む**必要があり、メッシュも手元・AWS とも無い。出口も `outlet_statPress` なので B1f の見直し対象 |
| 周期・軸対称の小規模試験 (#7) | **済 2026-09-24 (V6 全項目 PASS)**。以下は当時の記録: ~~**未**。~~case/23・case/43 ともメッシュが残っておらず生成から要る。case/43 の壁は `slip` で対象外 |

~~**未実施の 3 項目はいずれもメッシュ生成からになる**。~~ (3 項目とも 2026-09-25 までに済。)本 plan は既定 0 の opt-in で
`status: in_progress` なので、この状態でブランチに置いても既存ケースは踏まない (V2 で担保)。
**既定化の判断 (§5.1 #11) にはこの 3 項目が必須**。

#### 時系列判定 (2026-09-24, codex result-2 M1 の採用結果) — ~~**受入判定は保留**~~ **決着 (2026-09-25, #10b PASS — §6.2「#10b の再判定」節)**

M1「衝撃足 0.842 % は**最終壁 dump 1 枚**の比較で準定常判定が無い」を採用し、**共通 24 step の時系列**を作って
§6 で固定した閾値で判定した。**結果は `DRIFTING`**:

```
check_quasisteady.py --series-csv series_foot.csv --series-cols foot_L2_pct,dx_foot --drift 0.002 --osc 0.005
=== series_foot.csv [24 rows, steps 500..12000] -> DRIFTING ===
  foot_L2_pct : tail mean=0.8373  drift=0.2%/tail  fluct=1.3%   DRIFTING
  dx_foot     : tail mean=0       drift=0.0%/tail  fluct=0.0%   STEADY

check_quasisteady.py --series-csv series_force.csv --series-cols dFx_pct,dFy_pct --drift 0.0002 --osc 0.0005
=== series_force.csv [24 rows, steps 500..12000] -> DRIFTING ===
  dFx_pct : tail mean=0.149   drift=9.0%/tail  fluct=53.5%   DRIFTING
  dFy_pct : tail mean=0.2083  drift=4.4%/tail  fluct=25.8%   DRIFTING
```

| 量 | 終端差 | 末尾平均 | 変動 | 判定 |
| --- | --- | --- | --- | --- |
| 衝撃足の壁圧 L2 差 | 0.842 % | 0.8373 % | **1.3 %** | **DRIFTING** |
| 衝撃足の位置 $x_{foot}$ の差 | 0 | 0 | 0 % | STEADY |
| ランプ面 $F_x$ の差 | +0.123 % | 0.149 % | **53.5 %** | **DRIFTING** |
| ランプ面 $F_y$ の差 | +0.226 % | 0.208 % | **25.8 %** | **DRIFTING** |

→ **~~V3 SERN 2D PASS~~ を撤回**。正しくは「**終端差は許容内だが、差そのものが step ごとに 25–54 % 揺れており、
時系列判定は `DRIFTING`。受入判定は保留**」。終端 1 点の 0.031 % / 0.842 % は**揺れている量の瞬時値**である。
$x_{foot}$ だけは全 24 step で完全一致 (STEADY)。

**受入に必要なこと**: 本段をさらに延長して差の時系列が STEADY になるか確認する。
~~それまで V3 SERN 2D は**合格と書かない**。~~ **決着 (2026-09-25)**: この DRIFTING は差の系列に閾値を当てた operand 取り違えだった。#10b で各 run の $C_T$ に §6 の閾値を当て直して PASS。

#### V1-b/c の完全な記録 (2026-09-24, codex result-2 M4 の採用結果)

実行コマンド:
```
python3 ../cad/diag_wall_cv_budget.py sern.h5 --ids 153797,153880,189814 \
        --equilibrium --wall-normal-chi --wall-phys-ids 1,2,3,4,10,11,12,13,15
```

**床が不活性であることの証拠** (`roMin`=1e-4, `pMin`=20。§6 が求めた前処理は**実装していない**ので、
代わりに「床から十分離れている」ことを示す):

| dump | CV | $\rho$ | $P$ | $\rho/\rho_{Min}$ | $P/p_{Min}$ |
| --- | --- | --- | --- | --- | --- |
| res_24000 | 153797 | 1.50683e-03 | 5.12844e+02 | 15.1 | 25.6 |
| res_24000 | 153880 | 1.50641e-03 | 5.12701e+02 | 15.1 | 25.6 |
| res_24000 | 189814 | 1.62654e-03 | 5.53529e+02 | 16.3 | 27.7 |
| res_30000 | 153797 | 1.50683e-03 | 5.12844e+02 | 15.1 | 25.6 |
| res_30000 | 153880 | 1.50641e-03 | 5.12701e+02 | 15.1 | 25.6 |
| res_30000 | 189814 | 1.62677e-03 | 5.53608e+02 | 16.3 | 27.7 |
| res_36000 | 153797 | 1.50683e-03 | 5.12843e+02 | 15.1 | 25.6 |
| res_36000 | 153880 | 1.50641e-03 | 5.12700e+02 | 15.1 | 25.6 |
| res_36000 | 189814 | 1.62678e-03 | 5.53608e+02 | 16.3 | 27.7 |

**全接続面の正味流束と平衡予測** ($\Delta t_l$ を出力していないので正規化は $\lvert\Sigma\dot m\rvert/(\rho V)$ [1/s] で示す):

| dump | CV | $\Sigma\dot m$ [kg/s] | $\lvert\Sigma\dot m\rvert/(\rho V)$ [1/s] | 予測 $P_w$ | 観測 $P_w$ |
| --- | --- | --- | --- | --- | --- |
| res_24000 | 153797 | +4.1380e-14 | 2.34 | 5.1284e+02 | 5.12844e+02 |
| res_24000 | 153880 | +2.0696e-14 | 1.02 | 5.1270e+02 | 5.12701e+02 |
| res_24000 | 189814 | −3.9406e-13 | 8.55 | 5.5353e+02 | 5.53529e+02 |
| res_30000 | 153797 | +3.0223e-14 | 1.71 | 5.1284e+02 | 5.12844e+02 |
| res_30000 | 153880 | −5.3865e-16 | 0.027 | 5.1270e+02 | 5.12701e+02 |
| res_30000 | 189814 | −8.0313e-14 | 1.74 | 5.5361e+02 | 5.53608e+02 |
| res_36000 | 153797 | +5.2086e-14 | 2.95 | 5.1284e+02 | 5.12843e+02 |
| res_36000 | 153880 | +3.3228e-14 | 1.64 | 5.1270e+02 | 5.12700e+02 |
| res_36000 | 189814 | −2.0785e-14 | 0.45 | 5.5361e+02 | 5.53608e+02 |

起点 (`run_0432/res_0`) の $\Sigma\dot m$ = +1.488e-09 kg/s に対し、**最大残差 3.94e-13 は 3.58 桁の縮小**
(~~5 桁~~ は誤り。codex result-2 M4)。符号は CV・dump で揺れており、**正味流束は数値ゼロに張り付いている**。

**単体試験の範囲 (M4 の採用)**: `cad/test_diag_wall_cv_budget.py` の `kernel_mdot` は
**CUDA を実行せず Python でカーネルの式を書き直したもの**であり、**float64・mask は手で与えている**。
したがってこれは「**式の確認**」であって実カーネルの照合ではない。
**float32 の丸め・実際の `wall_flag`・非対象面のビット不変性は未試験**。

#### codex result 段の指摘への対応表 (2026-09-25, 3 回目の前に作成、4 回目の前に 3 回目の 7 件を追加)

前 2 回 (§6.1、原文 `notes/reviews/2026-09-23-convection-slau-wall-normal-chi-result{,-2}.md`) の Major 10 件。
**前 2 回とも同じ型「測れた範囲を超えて書いた」**だったので、状態欄は「何を測って閉じたか」で書く。

| 指摘 | 内容 (要約) | 対応 | 証拠 | 状態 |
| --- | --- | --- | --- | --- |
| R1-M1 | V1-c が別の流束で判定されていた | V1-c を撤回し、ツール `--wall-normal-chi` で再判定 → 実カーネルの面流束照合 (V5) | §6.2「V1-b/c の完全な記録」、V5 結果 (#10d, P0–P5 PASS) | 閉 (V5 で実カーネル) |
| R1-M2 | 出口 BC と $\chi_n$ の独立性が示せていない | `outflow`+flag 0 (`run_0440`) を足して 2×2 を揃えた | §6.2「効果の分離 — 2×2 が揃った」 | 閉 (「この窓では交絡は検出されなかった」に限定、m6) |
| R1-M3 | V3 の一部だけで受入全体を PASS にした | 受入ゲートを戻し、case/48・case/16・SERN 2D の 3 件をすべて判定 | case/48 V3 (§6.2)、case/16 (§6 V3 case/16 条件の結果、2026-09-25)、SERN 2D (#10b) | **限定** (case/48・case/16 は閉。SERN 2D は $C_T$ と衝撃足の各 run 準定常 (#14) が閉。格子感度は受入から外した → #9b) |
| R1-M4 | V1-d の準定常判定を方程式全体の収束証明に置き換えた | 4.40 桁に訂正・「ノイズ床」断定を削除し V1-d を限定合格に | §5.1 #6b、§6.2 V1-d | 閉 |
| R1-M5 | 「剥離縁だけに効く」は実装から成立しない | 撤回し「$\Delta p$ に比例」に訂正 (単体試験) | §6.2 (撤回箇所の取り消し線) | 閉 (R2-M5 で残存分も処置) |
| R2-M1 | SERN の準定常ゲートが未完、受入条件が既定化条件へ移った | case/16 を受入ゲートに戻し (#10)、SERN は時系列を #10b の登録式で再判定 | §6.2「#10b の再判定」(PASS)、#10 (V2・V3 PASS) | **閉** ($C_T$・case/16・衝撃足の各 run 準定常 (#14, 両 run ALL STEADY) は閉。$C_L$/$C_M$ の定常差はゲート量でなく #11 の前提) |
| R2-M2 | 累積 CFL を揃えただけでは固定点の一致を示せない | 両 CFL で全保存量の判定 + 独立延長区間で比較 (V7) | V7 結果 (#10c、限定) | **限定** (ΣCFL を揃えた CFL 間差 ≤ ε まで。固定点は #11 の前提へ; codex result-3 M2, #15) |
| R2-M3 | 2 水準から格子非依存と発現速度の一般則を導いた | 3 水準目 (粗) を追加。成否だけを言い、NaN step の機序説明を撤回 | §6.2「V3 格子感度」(限定 3 つ) | **限定** (起動成否 3 水準 #9a は閉。定常解・壁圧・収束率の格子比較 #9b は **未** で #11 の前提へ; codex result-3 M4, #16) |
| R2-M4 | V0・V1-b/c の全面合格を支える記録が不足 | 単体試験を「Python による式の確認」と明記し、実カーネル照合を V5 で実施。3 CV × 3 dump を記録 | V5 (#10d)、§6.2「V1-b/c の完全な記録」 | **閉** (V0・V1-c・V5 に加え、V1-b を #13 で測定: 起点で正味流入、プラトー後 最大 0.060 %/dump < 1 %) |
| R2-M5 | M5 の撤回が本文に反映しきれていない | 残っていた「剥離縁だけ」「元々 $\chi=0$」を取り消し線で撤回 | §6.2 該当箇所 | 閉 |
| R3-M1 | V1-b の $\Delta t_l$ 込み収支と起点直後の符号が未測定、受入結論で V1-c の量を V1-b と書いた | 文言を直し、AWS で測定 (#13) | §5.1 #13、`case/46.sern_design/_r3_m1m3/{run0437_budget.txt,V1B_b_VERDICT.txt}` | **閉** (起点で正味流入、プラトー後 最大 0.060 %/dump < 1 %) |
| R3-M2 | V7 は「CFL に依存する固定点」まで示していない | 固定点の断定を撤回し、測った範囲 (ΣCFL を揃えた CFL 間差 ≤ ε) に限定。固定点は #11 の前提へ (#15) | V7 結果 (限定)、§6 V3 箇条の書き換え | **限定** (受入条件から外した) |
| R3-M3 | SERN 衝撃足は差系列だけで、各 run の準定常が未判定。$C_L$ の「実効果」を DRIFTING のまま確定 | 各 run の量 (絶対 $x_{foot}$・窓内壁圧) を判定 (#14)、$C_L$ の確定を撤回 | §5.1 #14、`case/46.sern_design/_r3_m1m3/PERRUN_QS_flag{0,1}.txt` | **閉** (両 run ALL STEADY。$C_L$/$C_M$ は #11 の前提) |
| R3-M4 | 起動成否 3 水準を格子感度ゲートの完了として扱った。撤回済み機序が残存 | #9 を 9a (済)/9b (未) に分割、格子感度を受入から外し #11 の前提へ、残存 2 箇所に取り消し線 (#16) | §5.1 #9a/#9b、V3 表 | **限定** (受入条件から外した) |
| R3-M5 | `stage_manifest` の hard キーに新フラグが無い | キー追加 + 省略 ≡ 0、正規表現保険経路の区切り記号も除去 (#12) | `solver_density_cuda/tools/test_stage_manifest_wall_normal_chi.py` 7/7 PASS | **閉** |
| R3-m6 | 台帳・README・methods に撤回済みの数値・主張が残存 | 同期 (#17)。case/48 の最大差は丸め前の値 0.0018 % / 0.0021 % を正本に | §5.1 #17 | **閉** |
| R3-m7 | case/16 の「反復ノイズ以下」と機序の確定が言い過ぎ | 「500 step 試験でノイズ床比 0.80–1.47 (許容 2 倍以内)」「期待値と整合 (機序は未確定)」に限定 (#18) | 受入の範囲 | **閉** |

#### 受入の範囲 (result 段の結論, 2026-09-25, `diagnostician` 確認)

**言ってよいこと**:
- `space.slauWallNormalChi` (opt-in、既定 0) は **node・`nodeWallDirichlet: 1`・非周期/非軸対称**の構成で受入 PASS。
  - V0 と V5 (実カーネルの面流束、P0–P5)。
  - V2 (flag 0 がコード変更前 `39526328` と反復ノイズ床比 ≤ 2): case/48 と case/16。
  - V1 (case/46 接続模型): 壁 CV の排出が止まり床到達 0、要求系列 ALL STEADY、**V1-c** (平衡壁圧の観測/予測、合格範囲 ×/÷1.5; **全接続面収支の凍結近傍求根で、§4.2 の単面近似 (合否に使わない) ではない**) 1.000。**V1-b** (#13): 起点で 3 CV とも正味流入、プラトー後の $\lvert\Sigma\dot m\rvert\Delta t_l/(\rho V)$ 最大 0.060 %/dump (< 1 %、裕度 17 倍。$\Delta t_l$ は run_0439 後の再ビルドで出力、床前処理はツール外で「床の 15 倍以上」で代替)。
  - V3: case/48 は $C_f$/$q_w$ ±1 % 内。case/16 は壁 $p/p_0$ の差が V3 許容内で、**この 500 step 試験ではノイズ床比 0.80–1.47 (許容 2 倍以内)**。SERN 2D は $C_T$ が両 run STEADY で、#10b 式 (末尾平均差 0.031 % + 振幅 0.007 % =) **0.038 %** (許容 0.1 %)。衝撃足は**各 run が ALL STEADY** (#14: 絶対 $x_{foot}$ 不変・窓内壁圧 4 量) を確かめたうえで差ノルム L2 **0.8385 ± 0.0091 %** (平均+振幅 0.848 % ≤ 1 %)、$x_{foot}$ 差 0。
- 収束の悪化 (§3.1 `mSLAU` の報告): **SERN 2D は登録量で PASS** (本段 ΣCFL 6000 で全 8 列、flag 1 の低下桁数 ≥ flag 0 − 0.2 桁。ただし両 run ともプラトーで低下は 0.2–0.6 桁しかなく判別力は限られる)。**case/16・case/48 は収束場からの分岐で低下桁数が定義できないので代替量**: case/16 は事前固定 (c) の残差水準比 1.0003–1.0024 (≤ 2)、case/48 は事後の水準比 (`rms_ro` 1.83e-7 対 1.94e-7)。
- 周期・軸対称は V6 の**小規模試験 PASS まで** (生産規模では未使用)。
- 格子: 3 水準で「対策なし NaN / 対策あり完走」の**成否が同じ**。**格子非依存の証明ではない** (粗い側の flag 1 は 6000 step まで、定常解・壁圧・収束率の格子比較はしていない)。
- 固定点: **固定点の確認はしていない** (#11 の前提)。言えるのは、測定区間 (ΣCFL 13200–17200) の場の変化と、ΣCFL を揃えた CFL 0.2/0.4 間差が全列 ε 以下まで (差が検出される列はある: roe, P; 床の 5.9–8.6 倍)。
- case/16・case/48 で差が小さいことは、**壁隣接内点 $M \ll \sqrt2$ で $\chi_n \approx \chi$ という期待値と整合する** (当該ケースの $\Delta\chi$・$\Delta p$ は未測定なので機序は確定していない)。効果の根拠は V1/V5 のみ。

**言ってはいけないこと**:
- 「効果は無視できる」の一般化 (SERN 2D の $C_L$ / $C_M$ は現区間の末尾平均差 0.437 % / 0.320 % が振幅の 4–5 倍あるが、flag0 が `DRIFTING` で定常差として未確定)。
- 「$C_L$ への実効果を確定」(撤回、#14)。「固定点は CFL に依存する」(撤回、#15)。
- 「既定化できる/すべき」(§5.1 #11、スコープ外。前提は独立判定)。
- 「格子によらない」「固定点は同一」「剥離縁だけに効く」「3D の面数比で 2D の影響位置を語る」(撤回済み)。
- 「実験と一致」(case/16 の前提ゲートは流れている場の確認で、参照は非粘性 1D と 2D SST の run であって実験ではない)。
- 受入 = **「害が無い」の合格であって効果の検証ではない**。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-23` | [2026-09-23-convection-slau-wall-normal-chi-plan.md](../../notes/reviews/2026-09-23-convection-slau-wall-normal-chi-plan.md) | **GO-with-changes**, C0/M5/m2 | **全件採用** (却下なし)。M1 → §4.1 (質量流束限定) + §3.1 (文献) + §5.1 #2 を実装前へ。M2 → §4.2 を近似モデルへ降格 + §6 V1-b/c。M3 → §6 V1-e (期限方式) + 153880 を監視に追加。M4 → §6 V3 (受入保留・許容差の 1/5)。M5 → §4 (再構成後の面速度・条件付きの主張) + §6 V0 + §2 制限事項 + §5.1 #7。m6 → §6 V2 (ツール名 `check_field_regress.py` に訂正・ノイズ床比)。m7 → §5 の 1–3 (wrapper 配線・docs 先行) |
| **result** | `2026-09-23` | [2026-09-23-convection-slau-wall-normal-chi-result.md](../../notes/reviews/2026-09-23-convection-slau-wall-normal-chi-result.md) | **NO-GO**, C0/M5/m2 | **全件採用** (却下なし。2026-09-23 `diagnostician` 判断)。**M1 → V1-c を撤回**しツールに `--wall-normal-chi` を実装 + 単体試験でカーネル照合、`run_0439` の STEADY 後 3 dump で V1-b/c を再判定 (観測/予測 1.000)。**M2 → `outflow`+flag 0 (`run_0440`) を実施**し 2×2 を揃えた (排出の軌跡が 5 桁一致 = 独立)。**M3 → §5.1 に不足ゲートを戻し V3 総合 PASS を保留**。**M4 → 4.40 桁に訂正・「ノイズ床」断定を削除**し V1-d を限定合格に。**M5 → 「剥離縁だけ」を撤回**し「$\Delta p$ に比例」に訂正 (単体試験で実証)。m2 件も採用 |
| **result (2)** | `2026-09-23` | [2026-09-23-convection-slau-wall-normal-chi-result-2.md](../../notes/reviews/2026-09-23-convection-slau-wall-normal-chi-result-2.md) | **NO-GO**, C0/M5/m2 (ただし「**既定 0・`in_progress` で留める判断は妥当。`accepted` へ移すことは支持しない**」と明記) | **全件採用** (却下なし、2026-09-24 `diagnostician` 判断)。**全 5 件が同じ型「測れた範囲を超えて書いた」**。M1 → 衝撃足と力係数の**時系列を判定したら `DRIFTING`**、V3 SERN 2D の PASS を撤回し受入保留 (#10b)。case/16 を受入ゲートに戻した。M2 → 固定点 PASS を「終端場差 0.011 %」へ格下げ (#10c)。**「0.01 % 以内」は有利な側への丸めだった**。M3 → 「試した 2 格子で同じ成否」に限定し、NaN step 比の機序説明を撤回 (非再現指標で比を取っていた)。M4 → 単体試験を「Python による式の確認」と明記し実カーネル照合を #10d に。**「5 桁縮小」は 3.58 桁の誤り**。3 CV × 3 dump の完全な記録を §6.2 に。M5 → 本文に残っていた「剥離縁だけ」「元々 $\chi=0$」を取り消し線で撤回。m6 → 「独立と確定」を「この窓では検出されなかった」に限定 |
| result (3) | `2026-09-25` | [2026-09-25-convection-slau-wall-normal-chi-result.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-result.md) | **NO-GO**, C0/M5/m2 (実装・V5・case/16・case/48 の V2/V3 は裏付けありと明記、case/16 ゲート差し替えは不合格理由にしない) | **全件採用** (却下なし、2026-09-25 `diagnostician` 判断)。M5 → §5.1 #12 (事実を手元で確認済み)。M1 → #13 (V1-b/V1-c の取り違えを直し、AWS で $\Delta t_l$ 込み収支を測る)。M3 → #14 (各 run の衝撃足準定常を AWS で、$C_L$「実効果」は撤回)。**M2 → #15、M4 → #16 は「測って閉じる」でなく受入条件を測れた範囲に書き直し、固定点と格子収束を #11 の前提へ移す**。m6 → #17、m7 → #18。4 回目は #8b |

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{hpp,cpp}`
- `solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu` (**wrapper。`wall_flag_d` の受渡し**)
- `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`
- `methods/convection/{theory,implementation}.md`, `methods/index.md` (**実装前に更新**)
- 既定 0 なので既存ケース・手順への影響は無い (V2 で担保)。
- **制限事項**: flag 1 は**非周期・非軸対称の node 構成でのみ検証**する (§2)。周期・軸対称で使う前に §5.1 #7 が要る。
- block-DPLUR の Jacobian は近似 FVS のままでよい (`block_dplur_jacobian_d.cuh:22`)。
  ただし CFL・内反復数による固定点と収束性は V3 で確認する。

## 8. 完了条件

- [ ] 関連 `methods/convection/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 の V1–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status` を `done` に変更し §9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-09-24` — **codex result 段 2 回目も NO-GO (C0/M5/m2)、全件採用**。**5 件すべてが同じ型「測れた範囲を超えて書いた」**。衝撃足と力係数の**時系列を判定したら `DRIFTING`** (変動 25–54 %) で **V3 SERN 2D の PASS を撤回**し受入保留に。固定点 PASS を「終端場差 0.011 %」へ格下げ (**「0.01 % 以内」は有利な丸め**)。格子感度の結論を「試した 2 格子で同じ成否」に限定し NaN step 比の機序説明を撤回。単体試験を「Python による式の確認」と明記 (**「5 桁縮小」は 3.58 桁の誤り**)。本文に残っていた「剥離縁だけ」を撤回。case/16 を**受入ゲートに戻した**。3 回目のレビューは受入セットが揃うまで回さない (§5.1 #8)。
- `2026-09-23` — **codex result 段 NO-GO (C0/M5/m2) を全件採用**して復旧。**M1**: 診断ツールが変更前の $\chi$ を計算していたため V1-c を撤回 → `--wall-normal-chi` を実装し**単体試験でカーネル照合**、`run_0439` (ALL STEADY) の最後 3 dump で再判定し**観測/予測 1.000** (V1-b も $\Sigma\dot m$ が 5 桁縮小)。**M2**: 欠けていた `outflow`+flag 0 (`run_0440`) を実施し **2×2 が揃った** — 排出の軌跡が 5 桁一致し **2 つの欠陥は独立**と確定。**M3**: V0 単体・衝撃足 (事前固定ノルムで L2 差 0.842 % PASS)・`cfl_pseudo` 0.4 の固定点 (累積 CFL を揃えて 0.01 % 一致) を実施、格子感度と case/16 は §5.1 #9/#10 に残す。**M4**: 4.40 桁に訂正しノイズ床の断定を削除。**M5**: 「剥離縁だけ」を撤回し「$\Delta p$ に比例」に訂正 (単体試験で実証)。**status を `in_progress` に**。
- `2026-09-23` — 初稿。
- `2026-09-23` — 文献調査と `diagnostician` の判断で §3/§4/§6 を改訂。**A を本命のまま維持**し、根拠を「自作」から「**AUSM⁺-up の面法線マッハ切替則の局所適用**」に置き換え (`AUSM_d.cu:277-290`)。**§4.1 平衡壁圧の関係式**を追加し、V1 の合否をこれで判定する形に変更 (「完走」でも「$P_w\to P_i$」でもない)。**B (KEEP の ES 散逸流用) を却下** (σ=0.02 では平衡壁圧 0.02 で現状固定、A の fallback にもならない)、C は出典として吸収する方向で保留。[`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §4.13.1 の診断 (run_0430–0434) と
  `diagnostician` の設計方針を受けて起票。**実装前に codex plan 段が必要** (§5.1 #1)。
- `2026-09-23` — **codex plan 段 GO-with-changes (C0/M5/m2) を全件採用**して §2–§7 を全面改訂。**適用を質量流束のみに限定** (圧力束は現行のまま。面圧力が +32.5 % 変わるため)。**§4.2 平衡式を近似モデルへ降格**し V1 は全接続面の収支で判定 (無視していた W↔W' の補充が排出の 65 %)。V1 の合格条件を「常に $10\rho_{Min}$ 以上」から**期限方式 (step 200 までに超えて維持)** へ。V3 を受入ゲート化。**文献調査で `mSLAU` (全面適用版) の悪化報告を発見**し §3.1 に記載。
