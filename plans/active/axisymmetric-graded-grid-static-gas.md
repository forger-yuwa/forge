# 軸対称 node の静止ガスが非一様格子で静止を保てない — 偽流れの切り分けと修正

## メタ

- **area**: `axisymmetric`
- **status**: `draft`
- **related_docs**:
  - [`methods/axisymmetric/implementation.md`](../../methods/axisymmetric/implementation.md) — r 重み幾何・hoop ソース
  - [`methods/boundary.md`](../../methods/boundary.md) — slip / 等温壁 (node)
  - [`notes/investigations/node-slip-tangential-density-spurious-flow.md`](../../notes/investigations/node-slip-tangential-density-spurious-flow.md) — **既知の未修正欠陥**: node slip 境界 + 接線密度勾配で市松状の偽接線流
- **related_plans**:
  - [`boundary-cht-axisymmetric-fem2d.md`](boundary-cht-axisymmetric-fem2d.md) — **発注元**。§6 V-ax2b の非一様格子がこの障害で未達 (§5.1 #6)
  - [`axisymmetric-freestream-hoop-gauge.md`](axisymmetric-freestream-hoop-gauge.md) — hoop ソースの自由流保持 (pRef ゲージ・閉包面積)
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) — slip の作り直し予定 (既知欠陥の恒久修正先)
- **created**: `2026-09-27`
- **owner**: `sano`

## 1. 目的

軸対称 (`axisymMethod: 0`)・node で、**半径方向に非一様 (等比 1.1) な格子**の静止ガス層が、加熱 (両壁 350/325 K) のもとで
**静止を保てず偽の流れが立つ** (一様格子では立たない)。原因を切り分け、修正するか、使ってはいけない構成として拒否・文書化する。
完了時: 発注元 plan の V-ax2b (非一様格子) を登録どおりに判定できる状態。

## 2. スコープ

- **やる**: 原因の切り分け (A/B)、原因に応じた修正または構成の拒否・注意書き、回帰試験
- **やらない**: CHT 側の変更 (発注元 plan)、slip の全面作り直し (ghostless plan 側。本件が slip 起因ならそちらへ渡す)

## 3. 前提 — 観測事実 (2026-09-27、発注元 plan §5.1 #6 と `case/62.conjugate_disk/README.md`)

- 構成: 静止ガス層 $x\in[0,H{=}5\,\mathrm{mm}]$ × $r\in[5,20]$ mm (軸を含まない)、$x$=0 等温 350 K、$x=H$ 等温 325 K (連成前)、
  **$r$ 両端 slip**、`isAxisymmetric: 1`・`axisymMethod: 0`・node・`nodeWallDirichlet: 1`、定数物性、$p_0$ 1013.25 Pa、一様 IC 325 K・静止、
  陰解法 (block-DPLUR、`cfl_pseudo` 5、`implicitRelax` 0.5)、`convMethod` 1・`limiter` 2。**FP64 ビルド** (`6f64b873`)、`FORGE_CUDA_BLOCKSIZE=128`。
- 格子: $x$ 16 一様 × $r$ 32。**一様** (`case/62.conjugate_disk/run_0013_disk_r32u_w20k`) と **等比 1.1、$r$ が大きいほど広い** (`run_0014_disk_r32g_w20k`)。
- step 20000 (連成前) で: 一様 max|U| **4.3e-4 m/s**・P 1027 Pa / 非一様 max|U| **0.231 m/s** ($x$ 4.69, $r$ 6.19 mm)・**P 728–742 Pa**。
  非一様では $x=H$ 壁の熱流束が $r$≈5.9–6.6 mm で隣接節点ごとに符号反転 (±3000–4000 W/m²、伝導基準 120.5)。
- 起動直後 (step 250–500、`run_0015_hoop_A_hoop0`): max|U_y| **21.9 m/s**、半径方向圧力差 19.7 Pa、市松振幅 4.5e5 W/m² — 起動の過渡で最大。
- **棄却した仮説**: 入力メッシュの float32 閉包欠損 — `mesh.hoopAreaFromClosure` 0/1 の A/B で B/A = 1.001 / 0.998 / 0.996 (変化なし)。
- **既知の未修正欠陥** (2026-07-20、case/24): node の slip 境界に**接線方向の密度勾配**があると、slip 境界上の節点列に**市松状の偽接線流** (~0.5 m/s) が定在。
  slip → periodic で 1300 分の 1 に消える。本構成も $r$ 端の slip に沿って $x$ 方向の温度勾配がある。
  **ただし最大流速の位置は slip 境界上ではなく内部** ($r$ 6.19 mm = $r$=5 の slip から数セル) で、一様格子では出ていない — 既知欠陥と同一かは未確定。
- 発注元の同心円環 (`case/61`、一様格子、$x$ 端 slip・$r$ 方向の温度勾配 = slip の接線方向の勾配あり) は正常に収束し解析解と 0.04 % で一致した。

## 4. 設計方針 (草稿 — 上位に諮る前)

**切り分けを先にやり、修正方針は切り分けの後に決める**。候補:

- H1 **既知の slip 欠陥** (接線密度勾配 × slip 閉包) が非一様格子で増幅される
- H2 **加熱起動 × 低 Mach** の再構成・陰的更新 (等温壁ピンは密度を保って T・P・roe を変えるので起動直後は非平衡)
- H3 **r 重み幾何と非一様格子**の組合せ (閉包欠損は棄却済みだが、面重心 $r$・双対体積の扱いの非対称など)

判別 A/B の案 (1 変数ずつ):

- (i) **温度勾配を消す**: 両壁 325 K (加熱なし)。静止が保てれば「勾配駆動」(H1/H2)、保てなければ H3
- (ii) **平面 (`isAxisymmetric: 0`) の同じ非一様格子**: 出れば軸対称固有でない (H1/H2)、出なければ軸対称固有 (H3 か軸対称 × slip)
- (iii) **$r$ 端 slip → 等温壁 (no-slip)** または periodic (平面なら $r$ 方向 periodic が組める): 消えれば slip 起因 (H1)

## 5. 実装ステップ

1. 切り分け A/B (§4 の案を上位の判断で 1 つ〜2 つに絞る)
2. 原因に応じた修正 or 拒否・文書化
3. 回帰試験と発注元 V-ax2b の再判定

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 を上位に諮る | AGENTS.md 条件 1。判別 A/B を絞る | F |
| 2 | 切り分け A/B | #1 で決めたもの。AWS FP64、数千 step | O |
| 3 | 修正 or 拒否 | #2 の結果で決める (諮る) | F |

## 6. 検証 (草稿 — 上位に諮る前)

- **静止保持**: 対象構成 (非一様 32、加熱) で max|U| が一様格子と同程度 (≤1e-3 m/s を目安、要登録) かつ P が加熱に整合 (上昇)
- **発注元 V-ax2b**: 非一様 N_r=32 で全界面節点の壁温誤差 ≤0.5 % of 降下 (発注元の登録どおり)
- **回帰**: 一様格子・case/61・既知の slip ケース (case/24 `run_isoT_condL_node` 構成) で悪化しないこと

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- 修正が入る場合: node の slip 閉包 (`cuda_forge/boundaryCond_d.cu` 系) か軸対称幾何 (`variables.cpp`、`axisymmetricSource_d.cu`)。未定
- 既存ケース: slip と温度成層が同居する node ケース全般 (既知欠陥の含意)

## 8. 未確定事項

1. 原因 (H1/H2/H3)
2. 修正するか、構成の拒否・注意書きで済ませるか

## 9. 変更ログ

- `2026-09-27` — 起票 (ユーザ指示「じゃあ 1」= 流体側の別 plan を起票)。発注元 CHT plan の V-ax2b 非一様格子の障害。
