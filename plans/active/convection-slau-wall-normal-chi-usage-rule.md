# `slauWallNormalChi` の適用規則と、flag 1 が変える量の正否 (2D)

## メタ

- **area**: `convection`
- **status**: `draft` (codex plan 段 NO-GO → 改訂中)
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (「既知の限界」「対策 (opt-in)」節)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) (適用規則を書く先)
- **related_plans**:
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (前 plan。§5.1 #9b・#11、§6 V7、§6.2「未解決・適用限界」を引き継ぐ)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (3D 生産の格子収束列 R5n。3D の格子収束はこちらへ委譲)
- **created**: `2026-09-25`
- **owner**: `CFD Dev`

## 1. 目的

`space.slauWallNormalChi` を**いつ 1 にするか**の規則を決め、flag 1 が 2D 生産で動かす量の**大きさが生産許容 (R5n) に対してどこにあるか**と、
**その差が格子で減るか**を定量化する。**どちらが正しいかは決めない** (同一混合物性・同一 BC の参照手法が無い。§4.1)。**既定値は変えない**。

## 2. スコープ

- **やる**: 適用規則 (§4.2)、Q1 差の格子依存 (2D、3 水準 × 2 flag)、Q2 CFL 独立性 (2D、両フラグを独立に)、全 run の準定常 (§4.5)。
  ツール 2 点 (`v7_dist.py` の壁 ID・列の引数化、`stage_manifest` の `solver` hard キー)。
- **やらない**:
  - **既定化** (前 plan §5.1 #11 は「既定 0 のまま・規則で使う」で閉じた。§4.1)。
  - **正否の判定 (旧 Q3)**: **保留**。現行 Roe は EXH/AIR 2 成分 TP を `sp[0]` の単成分物性で再構成し (`convectiveFlux_roe_d.inc.cuh:159-162, 242-244`)、
    同一物理にならない。Roe の多成分化は別 plan (本 plan はコード変更なし)。**同一物理でない比較を「参考」として回さない**。
  - **3D の格子収束 (前 plan #9b)**: sern-3d の R5n 型を flag 1 で再取得する作業として [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 へ委譲。
  - **3D の固定点 (前 plan V7 のドリフト減衰)**: **2D の Q2 で解決した扱いにしない**。委譲先の未解決事項として残す。
  - SU2 を第三者にすること (frozen_tp 2 成分・等温壁 TP を同一 BC で組めない)。

## 3. 関連 docs と前提

- 前 plan の受入 (opt-in として、試験したケース・設定・区間で回帰許容内) と「未解決・適用限界」(§6.2)。
- 2D SERN 生産の起動レシピと判定: [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) §4.12 / §4.16 / §4.16.1、skill `sern-eval`。
- **2D 生産の収束状態 (訂正、codex plan M3)**: 2D 生産 (`case/46.sern_design/_r3_m1m3/sern2d_conv/flag{0,1}`) は flag 0/1 とも NaN なく完走するが、
  残差は**全列 `NOT CONVERGED (stalled/plateau)`、本段で 0.2–0.6 桁**。比較できるのは目的量 $C_T$ が §6 閾値で STEADY だから (前 plan #10b)。
  **「収束」とは書かない**。
- 許容は新設しない: 力係数は sern-3d §8 の R5n ($|\Delta C_T|\le0.002$, $|\Delta C_L|\le0.002$, $|\Delta C_M|\le0.02$、絶対値)、衝撃足は前 plan V3
  (L2 ≤ 1 %、$x_{foot}$ ≤ 格子 1 つ)、CFL 独立性は前 plan V7 の $\varepsilon$。
- 前 plan の実測: 2D 生産 (m6_on) で flag 1 − flag 0 = $C_L$ +0.437 % / $C_M$ −0.320 % (絶対 0.0012 / 0.023)、flag0 の $C_L$ は `DRIFTING`。
  **$C_M$ の差 0.023 は R5n 許容 0.02 を超える**。

## 4. 設計方針

### 4.1 既定化はしない / 正否は決めない (2026-09-25 `diagnostician` 判断、codex plan M1・M3 で改訂)

- flag 1 が**必要**な構成は、現状 1 つだけ: node + `nodeWallDirichlet: 1` で**側壁∩後端面の壁 CV が排出される 3D 接続構成**
  (flag 0 は 3 格子とも NaN、`case/46.sern_design/README.md` の run_0443–0446)。
- 2D 生産は flag 0/1 とも完走するが (残差はプラトー、§3)、flag 1 は $C_L$/$C_M$ を**正否不明のまま**動かす。既定 ON は既存 node 生産の答えを
  黙って変える操作で、利得が無い。規則 (opt-in を明示) の方が監査できる (`stage_manifest` の hard キー、`RUN_PROVENANCE`)。
- 本 plan で決めるのは (i) 適用規則、(ii) flag 差の大きさが R5n 許容に対してどこにあるか、(iii) その差が格子で減るか。**正否は決めない**。

### 4.2 適用規則 (codex plan M6 で改訂)

- **兆候 (診断を始める合図であって適用条件ではない)**: `diag_wall_cv_budget.py --summary` の床到達 > 0 (全域カウントで壁 CV に限らない)、
  壁隣接内点 $|\mathbf u|>\sqrt2\,\hat c$ (速度だけでは流束は変わらない。flag が変える項は再構成後の面状態での $(\chi_n-\chi)\Delta P$、
  `convectiveFlux_slau_d.inc.cuh:541`)。
- **適用判断 (3 条件すべて)**: (1) 対象壁 CV (`bcondConfig` の wall physID から作る) の $\rho_w$ が 3 dump 以上単調減少し $\rho_w/\rho_i<0.1$;
  (2) `diag_wall_cv_budget.py --faces-all` の全接続面正味流出 $\Sigma\dot m>0$ が 3 dump 以上持続; (3) 同 `--wall-normal-chi` で、再構成後の面状態に
  $\chi=0$ かつ $\Delta P\ne0$ の壁隣接面がある (= 補充項が消えている)。
- **適用範囲は既知の 3D node 側壁接続構成のみ**。周期・軸対称の生産構成へ一般化しない (V6 は小規模試験のみ)。**2D 生産は 0 のまま**
  (Q1 の結果で「差の大きさ」を注記する)。

### 4.3 Q1 格子 3 水準 (codex plan M2 で改訂)

**メッシャの `scale` は使わない** (`runner_sern.py:495` が `scale=H` を渡し、`mesh_sern.py:325` の `coords *= prm.scale` は**物理寸法**を変える。
縮めると Reynolds 数が変わり、固定形状の格子比較にならない)。`H_m`・形状・領域・`cowl_thickness`・`x_out_extra`・`bot_depth`・`top_depth`・
`vehicle_*`・`ramp_fillet` を固定し、`problem_moo_frozen_tp_cycle3op.yaml` の**解像度キーだけ**を変える:

| 水準 | `ni_up/ni_noz/ni_plume` | `nj_top/nj_bot/nj_ext_top/nj_wake` | `first_wall_frac` | `first_wake_frac` |
| --- | --- | --- | --- | --- |
| 粗 = 生産 | 16/160/180 | 81/51/41/9 | 0.004 | 既定 ($t_{base}/5$) |
| 中 (×√2) | 23/226/255 | 115/72/58/13 | 0.00283 | 既定/√2 |
| 細 (×2) | 32/320/360 | 161/101/81/17 | 0.002 | 既定/2 |

- **run ごとの記録**: 生成コマンド全文・節点数・衝撃足窓の局所接線間隔・3 station の第一内部ノード距離・**壁ノードの生産壁折れ線からの距離 max ≤ 1e-6 H**
  (形状不変)・`check_mesh_quality` VERDICT。
- **IC**: 生産 flag 0 の場を `interp_field.py` で各格子へ (cross-mesh、原始量から)、`warm_adapt` 500 step (1 次・soft CFL) → 本段 12000 を flag 0 で回し、
  その終端から `restart_field.py` で flag 0/1 を分岐して各 12000 (**別起動しない**)。

### 4.4 比較式 (codex plan M7、測る前に固定)

- 差は末尾 40 % の**平均差** $\Delta=\bar C^{(1)}-\bar C^{(0)}$。**上限 $U=\lvert\Delta\rvert+a_0+a_1$** ($a_f$ = 各 run の末尾 span/2。位相独立を仮定した保守側)。
  **合格は $U\le$ 許容**。
- 併記: **平均の安定性** $s_f=\lvert\bar C_{\text{前半}}-\bar C_{\text{後半}}\rvert$ (末尾窓の前半・後半) が許容の 1/5 以下でないと、平均差そのものを報告しない。
- 前 plan #10b の「両 run の半幅の大きい方」は過小 (位相独立なら両側が寄与) — 前 plan に訂正行を足した (§9)。

### 4.5 準定常 (codex plan M3 で改訂)

- **全 run (Q1・Q2) × 全量** ($C_T$, $C_L$, $C_M$, 衝撃足 L2, $x_{foot}$) に適用。閾値: $C_T$ `--drift 0.0002 --osc 0.0005`、$C_L$ `0.0014/0.0014`、
  $C_M$ `0.00056/0.00056`、衝撃足 `0.002/0.005`、`--tail 0.4`。
- **`DRIFTING` / `TRANSIENT-UNSETTLED` → 本段を 1 回だけ倍にして再判定。それでも同じなら「判定保留」** (比較しない・**ラベルを変えない**・予算上限を書く)。
  `OSCILLATING` → 平均±振幅で比較 (§4.4)。**未定常を「限界サイクル」と言い換えない**。
- 「残差プラトー」と「目的量が比較可能」を run ごとに**別欄**で記録する。

### 4.6 Q2 CFL 独立性 (codex plan M4 で改訂)

- 生産格子、各フラグ・各 CFL (`cfl_main` 0.5 と 0.25、**ΣCFL を揃える**)。**まず各 run の $D_{intra}$ (種→終端) と目的量の準定常を通し、それから $D_{inter}$**。
- 主判定は前 plan V7 と同じ「列ごと $\max\{D_{intra}(a), D_{intra}(b), D_{inter}\}\le\varepsilon$」を**フラグごとに独立に**。
  旧案の比 $D_{inter}(1)\le2D_{inter}(0)$ は削除 (丸め床付近で不安定)。比は検出限界つきの補助情報に留める。
- 壁集合は 2D の `bcondConfig` の wall physID から作る (`v7_dist.py` は 3D の壁 physID を固定しているので `--wall-ids` 引数化。第一内部ノードも同様)。
  列: `ro, roUx, roUy, roe, P` はゲート ($\varepsilon$: P 0.02 %、他 0.2 %)、`roY0/roY1` は ro と同じ 0.2 %、`roK/roOmega` は**監視のみ** (許容を新設しない)。

## 5. 実装ステップ

1. plan 改訂 + 前 plan #10b の訂正行 + codex plan 再レビュー。
2. ツール: `v7_dist.py --wall-ids/--cols`、`stage_manifest` に `solver` hard キー + 単体試験。
3. Q1 メッシュ 3 水準 (生成・記録・QC) + interp 起動 (flag 0)。
4. Q1 分岐 6 run + 準定常 + 判定。Q2 4 run。
5. 規則の本文化 (`recommended-settings.md`) + codex result。

### 5.1 残作業 (優先順)

**計算資源は AWS** (ユーザ指示 2026-09-25: ローカルで大きな計算をかけない。AWS は他セッションと共有なので起動前に `aws_instance.sh status`/`busy`、手動 stop はしない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | **plan 改訂 + codex plan 再レビュー** | 本改訂 (§4/§6)、前 plan #10b の訂正行。`codex_review.py --stage plan` を再度。Critical/Major は `diagnostician` に諮る | F |
| 2 | ツール 2 点 | `case/46.sern_design/cad/v7_dist.py` に `--wall-ids` / `--cols` 引数 (既定は現行の 3D 値でビット同一の出力)、`solver_density_cuda/tools/stage_manifest.py` に `solver` hard キー + 単体試験 (SLAU↔ROE が別区間・省略 ≡ 既定 SLAU)。合格: 単体試験 PASS、`test_gate_bad_input.py` PASS | O |
| 3 | Q1 メッシュ 3 水準 + interp 起動 (flag 0) | §4.3 の表のキーで生成、記録項目をすべて run に残す、`check_mesh_quality` PASS。`interp_field.py` → warm 500 → 本段 12000 (flag 0)。AWS、格子 ×2 で本段 ≈ 4 倍 → 計 ≈ 2–3 h | O |
| 4 | Q1 分岐 6 run + 準定常 + 判定 | `restart_field.py` で分岐、§4.5 の準定常を全量に、§6 Q1 の 3 結論で判定。AWS ≈ 2 h | O (結論は F) |
| 5 | Q2 4 run | 生産格子、`cfl_main` 0.5 / 0.25、ΣCFL 揃え、§4.6。AWS ≈ 1 h | O (結論は F) |
| 6 | 規則の本文化 + codex result | 結果を §6.2 に、規則を `procedures/recommended-settings.md` に (§4.2 の 3 条件)。`methods/convection/theory.md` からリンク | F |
| 7 | 3D の格子収束・固定点を sern-3d へ委譲 | `tooling-nozzle-sern-3d.md` §5.1 に「R5n を `slauWallNormalChi: 1` で再取得 (メッシュ生成コマンド・blocksize 記録)。3D 固定点は未解決」を 1 行 | F |
| — | 旧 Q3 (Roe) | **保留** (別 plan: Roe 多成分化)。同一物理でない比較は回さない | — |

## 6. 検証

**判定基準はすべて測る前に固定する** (2026-09-25 `diagnostician`、codex plan M3–M7 で改訂)。**新設の許容差は無い**。

| 項目 | 量 | 判定 (測る前に固定) |
| --- | --- | --- |
| Q1 差の格子依存 | 量ごとに (i) 各フラグの格子間変化 $\delta^{(f)}_k=C^{(f)}_{k+1}-C^{(f)}_k$、(ii) 各格子のフラグ差上限 $U_k$ (§4.4)。量 = $C_T, C_L, C_M$、衝撃足 L2、$x_{foot}$ | 結論は 3 つだけ。**検証範囲で許容内**: 3 格子すべてで $U_k\le$ 許容、かつ両フラグとも $\lvert\delta_{\text{細}}\rvert\le$ 許容。**差が残る**: 両フラグとも $\lvert\delta_{\text{細}}\rvert\le$ 許容なのに最細で $U>$ 許容 → 「flag は許容を超えて解を変える。正否は未判定 (参照手法なし)」、規則は変えない。**判定不能**: どちらかのフラグで $\lvert\delta_{\text{細}}\rvert>$ 許容、または $\delta$ の符号が水準間で反転 → 水準を 1 つ足す (1 回のみ)、それでも同じなら「格子収束が取れない」と書く。**比率 $\Delta_{細}/\Delta_{粗}$ は分岐に使わない** (表で推移を示すだけ)。**複数量は最悪の量の結論を全体の結論にする** |
| Q2 CFL 独立性 | 各フラグ・各 CFL の $D_{intra}$、ΣCFL を揃えた $D_{inter}$ (§4.6) | フラグごとに独立に、列ごと $\max\{D_{intra}(a), D_{intra}(b), D_{inter}\}\le\varepsilon$。比は補助情報 |
| 準定常 | 全 run × 全量 (§4.5) | STEADY / OSCILLATING (平均±振幅) のみ比較。DRIFTING / TRANSIENT-UNSETTLED は 1 回延長、それでも同じなら判定保留 |
| 収束 | 全 run | `check_convergence.py --segment` の VERDICT を同一設定区間で併記。「残差プラトー」と「目的量が比較可能」を別欄に |

- **参照値の出所**: Q1 は同格子の flag 0 と各フラグの格子間変化 (外部参照なし)、Q2 は前 plan V7 の $\varepsilon$、許容は sern-3d §8 の R5n 値と前 plan V3。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-25` | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用** (却下なし、2026-09-25 `diagnostician` 判断。根拠は各行をコードで確認)。M1 → Q3 (Roe) を保留、目的を「正否は決めない」に (§1・§2・§4.1)。M2 → `scale` をやめ解像度キーで 3 水準・形状不変の記録・IC 手順 (§4.3)。M3 → 「2D 生産は収束」を訂正 (§3)、準定常を全 run・全量に、未定常を言い換えない (§4.5)。M4 → Q2 を $D_{intra}$ 先行・フラグごと独立・比を削除、`v7_dist.py` 引数化 (§4.6、§5.1 #2)。M5 → Q1 を「格子間変化とフラグ差を別評価・3 結論」に (§6)。M6 → 兆候と適用 3 条件を分離、適用範囲を既知 3D 構成に限定 (§4.2)。M7 → 上限 $U=\lvert\Delta\rvert+a_0+a_1$ と平均の安定性 (§4.4)、前 plan #10b に訂正行。m8 → `stage_manifest` に `solver` hard キー (§5.1 #2) |

### 6.2 結果

(未実施)

## 7. 影響範囲

- コード変更なし (既定値は変えない)。`procedures/recommended-settings.md` に適用規則を追加。
- `methods/convection/theory.md` の「対策 (opt-in)」節に規則へのリンクを足す。
- 既存ケース・手順への影響なし (既定 0 のまま)。

## 8. 完了条件

- [ ] `procedures/recommended-settings.md` に適用規則を本文化 (暫定を外す)
- [ ] Q1・Q2 を §6 の基準で判定し §6.2 に記録 (旧 Q3 は保留)
- [ ] sern-3d §5.1 に R5n の flag 1 再取得と 3D 固定点の未解決を追加 (委譲)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-25` — codex plan 段 **NO-GO (C0/M7/m1)** を全件採用して改訂 (§6.1)。**目的を「適用規則 + flag 差の大きさと格子依存」に狭め、正否判定 (旧 Q3 Roe) は保留** (現行 Roe は 2 成分 TP を単成分物性で組むので同一物理でない)。格子 3 水準は `scale` でなく解像度キーで、判定規則・比較式・準定常の扱いを測る前に固定し直した。
- `2026-09-25` — 初稿。前 plan (accepted) の未解決 (#9b・#11・V7・$C_L$/$C_M$) を引き継ぐ。`diagnostician` の判断で**既定化を目的から外し**、
  適用規則の確定と 2D での正否判定に置き換えた (§4.1)。3D の格子収束は sern-3d へ委譲、3D 固定点は判断基準から外した。
