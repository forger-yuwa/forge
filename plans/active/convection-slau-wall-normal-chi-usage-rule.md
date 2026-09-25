# `slauWallNormalChi` の適用規則と、flag 1 が変える量の正否 (2D)

## メタ

- **area**: `convection`
- **status**: `draft` (codex plan 段 2 回 NO-GO → 縮小して改訂中)
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (「既知の限界」「対策 (opt-in)」節)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) (適用規則を書く先)
- **related_plans**:
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (前 plan。§5.1 #9b・#11、§6 V7、§6.2「未解決・適用限界」を引き継ぐ)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (3D 生産の格子収束列 R5n。3D の格子収束はこちらへ委譲)
- **created**: `2026-09-25`
- **owner**: `CFD Dev`

## 1. 目的

`space.slauWallNormalChi` を**いつ 1 にするか**の規則を書き、前 plan で未確定だった 2 点 — (a) 2D 生産での flag 差が生産許容の帯内か、
(b) カウル衝撃足での flag 差 (前 plan は膨張角を測っていた) — を**既存 dump と最小限の run** で確定して閉じる。**既定値は変えない**。
**どちらの flag が正しいかは決めない** (同一混合物性・同一 BC の参照手法が無い)。

## 2. スコープ

- **やる**: 適用規則 (§4.1)、既存 `_v3sern/*_ext` の $C_T/C_L/C_M$ 系列の再判定 (§4.3)、衝撃足の再同定 (§4.4)、前 plan の訂正、
  `stage_manifest` の `solver` hard キー。
- **やらない (理由つき)**:
  - **既定化**: 前 plan §5.1 #11 で「既定 0 のまま・規則で使う」と判断済み。
  - **Q1 差の格子依存 (2D 3 水準)**: 撤回 (2026-09-25 `diagnostician`、codex plan-2)。規則が「2D は 0」である以上**どの決定も変えない**うえ、
    `mesh_sern.py:347,359` が格子から機体形状を決めるので**形状固定の格子列が現状作れない** (codex が座標で確認: 3.97e-4 / 6.56e-4 H のずれ)。
    再開するときの要件は §7。
  - **Q2 CFL 独立性 (2D)**: 撤回。2D 生産で flag 1 を使わないので決定を変えない。3D (flag 1 必須) の CFL・格子は sern-3d の生産格子収束列で扱う。
  - **Q3 正否の第三者**: 保留。現行 Roe は EXH/AIR 2 成分 TP を `sp[0]` の単成分物性で再構成し (`convectiveFlux_roe_d.inc.cuh:159-162, 242-244`)、
    同一物理にならない。Roe の多成分化は別 plan。同一物理でない比較を「参考」として回さない。
  - **3D の格子収束・固定点 (前 plan #9b・V7)**: [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 へ委譲。**2D で解決した扱いにしない**。
  - SU2 を第三者にすること (同一 BC が組めない)。

## 3. 関連 docs と前提

- 前 plan の受入と「未解決・適用限界」、および 2026-09-25 の訂正 2 件 (#10b の不確かさ式、「衝撃足」の誤同定)。
- **生産許容 (現行)**: sern-3d §8 (`tooling-nozzle-sern-3d.md:1916-1920`、**2026-09-21 ユーザ決定で $C_M$ 0.02 → 0.05**):
  $|\Delta C_T|\le0.002$, $|\Delta C_L|\le0.002$, $|\Delta C_M|\le0.05$。**1 回目の設計 (2026-09-25) と codex 1 回目は旧値 0.02 を踏襲していた**
  (`diagnostician` の設計ミス、codex plan-2 M2 で訂正)。したがって「$C_M$ 差 0.023 は生産許容を超える」は撤回する。
- **2D 生産の収束状態**: flag 0/1 とも NaN なく完走するが、残差は**全列 `NOT CONVERGED (stalled/plateau)`、本段 0.2–0.6 桁**。比較できるのは目的量が
  §4.3 の閾値で STEADY/OSCILLATING のときだけ。「収束」とは書かない。
- 前 plan の実測 (末尾平均): flag 差 $C_T$ 0.0003 / $C_L$ 0.0012 / $C_M$ 0.023 (いずれも現行許容の帯内)。ただし $C_L/C_M$ の準定常は
  $C_T$ 用の閾値でしか判定されておらず、量別の閾値では**未判定** (flag0 は `DRIFTING` と記録)。

## 4. 設計方針

### 4.1 適用規則 (`procedures/recommended-settings.md` に書く文、codex plan M6・plan-2 M6 で改訂)

「`space.slauWallNormalChi: 1` は **node + `nodeWallDirichlet: 1` の 3D 側壁∩後端面接続構成 (case/46 接続模型・SERN 3D 生産) でのみ**使う。
2D 生産・case/16・case/48・周期・軸対称は 0 (周期・軸対称は小規模試験 V6 まで)。新しい構成で 1 にするのは、**`convMethod: 0` の起動区間**で
`diag_wall_cv_budget.py` (1 次 SLAU の再計算) により (i) 対象壁 CV (`bcondConfig` の wall physID から作る) の $\rho_w$ が 3 dump 以上単調減少し
$\rho_w/\rho_i<0.1$、(ii) 全接続面の正味流出 $\Sigma\dot m>0$ が 3 dump 以上持続、(iii) 元の $\chi=0$ かつ $\Delta P\neq0$ の壁隣接面があり、
同じ面状態で $\chi\to\chi_n$ に置換した再計算で当該 CV への補充 $\dot m$ が増える (元 $\chi$・$\chi_n$・$\Delta P$・置換前後の $\dot m$ を併記)、
の 3 つを満たしたときだけ。2 次生産場での診断は再構成後の面状態が取れるまで保留。flag は `stage_manifest` の hard キーなので段の途中で切り替えない。」

- `--summary` の床到達数 (全域) と「壁隣接内点 $|\mathbf u|>\sqrt2\hat c$」は**診断を始める兆候**であって適用条件ではない。
- `--wall-normal-chi` は「$\chi$ を $\chi_n$ に置換して再計算する」オプションであって「元の $\chi$ を確認する」ものではない
  (`diag_wall_cv_budget.py:75-77`)。docstring に明記する (§5.1 #5)。

### 4.2 比較式 — 差区間 (codex plan-2 M4・M5、測る前に固定)

- 各 run の末尾窓 (末尾 40 %) の観測範囲 $I_f=[\min C_f, \max C_f]$ から、**差区間** $I_\Delta=[\min C_1-\max C_0,\ \max C_1-\min C_0]$。
- 判定 (許容帯 $B=[-\text{tol}, +\text{tol}]$): **$I_\Delta\subset B$ → 帯内**、**$I_\Delta\cap B=\emptyset$ → 差が残る**、**それ以外 → 判定不能**。
  $I_\Delta$ は観測窓内の範囲であり、平均値の統計的不確かさ (末尾窓の前半/後半の平均差) とは**別欄**に記録する。
- 旧式 $\lvert\Delta\rvert+a_0+a_1$ (1 回目) と前 plan #10b の式は、非対称な変動で上限にならない (codex の反例: 旧式 0.019、実最大差 0.023)。
  前 plan #10b の $C_T$ と膨張角窓もこの差区間で再計算し、訂正行を足す (§5.1 #2。判定は計算してから書く)。

### 4.3 既存系列の再判定 (新規 run なし)

- 対象: 前 plan の `_v3sern/v3s_flag{0,1}_ext` (72 dump、36000 step)。量 $C_T, C_L, C_M$ (`forge_design.metrics.sern_forces.force_history` で作った
  既存 `ct_flag{0,1}.csv`)。
- **量別の準定常閾値** (絶対許容の 1/5 を各量の末尾平均で相対化。出所を併記): $C_T$ `--drift 0.0002 --osc 0.0005` (前 plan V3 表)、
  $C_L$ `0.0014/0.0014` (= 0.002/0.277/5)、$C_M$ `0.0014/0.0014` (= 0.05/7.16/5)、`--tail 0.4`。
- **`DRIFTING` / `TRANSIENT-UNSETTLED` → 「判定保留」**(ラベルを変えない。延長 run は本 plan ではしない: 2D は flag 0 で生産するので決定を変えない)。
  `STEADY` / `OSCILLATING` の量だけ §4.2 の差区間で判定する。

### 4.4 衝撃足の再同定 (codex plan-2 M7、測る前に固定)

- **探索区間**: ランプ上の物理座標 $x\in[L_{cowl}, L_{ramp}]$ (カウル唇より下流のランプ)。**条件**: $dp_w/dx>0$ (圧縮) の最大位置 $x_f$、
  かつ窓 $[x_f-2t, x_f+5t]$ で $p_{後}/p_{前}\ge1.1$ ($t$ = カウル板厚 2 mm)。圧縮・膨張を区別しない旧抽出 (`v3sern_foot_perrun.py:42`,
  `v3sern_series.py:39-40`) は使わない。
- **比較**: 両 run の壁圧を**共通の物理 x 格子** (窓内の局所間隔) に補間し、相対 L2 ≤ 1 %、$x_f$ 差 ≤ 局所間隔 1 つ (前 plan V3 と同じ許容)。
- **各 run の準定常**: $x_f$ と窓内壁圧を `0.002/0.005` で判定 (§4.3 の規則)。差は §4.2 の差区間。
- **圧縮足がランプ上に無ければ**「設計点 m6_on ではカウル衝撃がランプに当たらない」と記録し、**m4_off の flag 0/1 対** (2D、生産レシピ、flag 0 から分岐) を
  1 組だけ回して同じ手順で判定する (AWS)。

## 5. 実装ステップ

1. plan 縮小 + 前 plan の訂正 + codex plan 3 回目。
2. §4.3 の再判定と §4.4 の再同定 (既存 dump。必要なら m4_off の 1 組)。
3. `stage_manifest` の `solver` キー、`diag_wall_cv_budget.py` の docstring。
4. 規則の本文化、codex result、accepted。

### 5.1 残作業 (優先順)

**計算資源は AWS** (ユーザ指示 2026-09-25: ローカルで大きな計算をかけない。AWS は他セッションと共有なので起動前に `aws_instance.sh status`/`busy`、手動 stop はしない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | **plan 縮小 + codex plan 3 回目** | 本改訂。Critical/Major は `diagnostician` に諮る | F |
| 2 | 既存系列の再判定 (§4.3) + 前 plan #10b の差区間での再計算 | `_v3sern/ct_flag{0,1}.csv` (AWS から csv だけ取得) に量別閾値で `check_quasisteady.py --series-csv`、`STEADY`/`OSCILLATING` の量だけ差区間 (§4.2)。前 plan #10b の $C_T$・膨張角窓も差区間で再計算し訂正行 | O |
| 3 | 衝撃足の再同定 (§4.4) | 新スクリプト (探索区間・圧縮条件・共通 x 補間) を**結果を見る前に**書いて commit、既存 72 dump に当てる。圧縮足が無ければ m4_off の 1 組 (AWS ≈ 30 min)。結論は `diagnostician` | O (結論 F) |
| 4 | 前 plan の訂正行 | **済 2026-09-25**: V3 表・§6.2・#14・受入の範囲・§9 を「膨張角窓の壁圧」に改名し、カウル衝撃足を「未測定」に | O |
| 5 | ツール | `stage_manifest.py` に `solver` hard キー + 単体試験 (SLAU↔ROE が別区間・省略 ≡ 既定 SLAU)、`diag_wall_cv_budget.py` の docstring に `--wall-normal-chi` の意味。合格: 単体試験 PASS、`test_gate_bad_input.py` PASS | O |
| 6 | 規則の本文化 + codex result → accepted | §4.1 の文を `procedures/recommended-settings.md` に、`methods/convection/theory.md` の対策節からリンク | F |
| 7 | sern-3d への委譲 | `tooling-nozzle-sern-3d.md` §5.1 に「R5n を `slauWallNormalChi: 1` で再取得 (メッシュ生成コマンド・blocksize 記録)。3D 固定点は未解決」を 1 行 | F |

## 6. 検証

**判定基準はすべて測る前に固定する**。**新設の許容差は無い** (生産許容は sern-3d §8 の現行値、衝撃足は前 plan V3)。

| 項目 | 量 | 判定 (測る前に固定) |
| --- | --- | --- |
| 既存系列 (§4.3) | 各 run の $C_T, C_L, C_M$ | 量別閾値で `STEADY`/`OSCILLATING` の量だけ差区間 $I_\Delta$ を判定: $I_\Delta\subset B$ → 帯内 / $I_\Delta\cap B=\emptyset$ → 差が残る / それ以外 → 判定不能。`DRIFTING`/`TRANSIENT-UNSETTLED` は「判定保留」 |
| 衝撃足 (§4.4) | 各 run の $x_f$・窓内壁圧、共通 x 格子での相対 L2、$x_f$ 差 | 各 run が `STEADY`/`OSCILLATING` のときだけ比較。相対 L2 は差区間で ≤ 1 % → 帯内 (他は §4.2 の 3 値)、$x_f$ 差 ≤ 局所間隔 1 つ。圧縮足が見つからなければ m4_off で同じ判定 |
| ツール | `stage_manifest` 単体試験 | PASS |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-25` | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用** (却下なし、2026-09-25 `diagnostician` 判断。根拠は各行をコードで確認)。M1 → Q3 (Roe) を保留、目的を「正否は決めない」に (§1・§2・§4.1)。M2 → `scale` をやめ解像度キーで 3 水準・形状不変の記録・IC 手順 (§4.3)。M3 → 「2D 生産は収束」を訂正 (§3)、準定常を全 run・全量に、未定常を言い換えない (§4.5)。M4 → Q2 を $D_{intra}$ 先行・フラグごと独立・比を削除、`v7_dist.py` 引数化 (§4.6、§5.1 #2)。M5 → Q1 を「格子間変化とフラグ差を別評価・3 結論」に (§6)。M6 → 兆候と適用 3 条件を分離、適用範囲を既知 3D 構成に限定 (§4.2)。M7 → 上限 $U=\lvert\Delta\rvert+a_0+a_1$ と平均の安定性 (§4.4)、前 plan #10b に訂正行。m8 → `stage_manifest` に `solver` hard キー (§5.1 #2) |
| plan | `2026-09-25` (2 回目) | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md) | **NO-GO**, C0/M7/m0 | **全件採用** (2026-09-25 `diagnostician` 判断)。**M2 で前提が変わった** ($C_M$ 許容は現行 0.05、flag 差は全量帯内) ため **Q1/Q2 を撤回して plan を縮小** (§2)。M1・M3 → Q1 撤回で解消し、再開時の必須要件として §7 に残す。M4・M5 → 差区間 $I_\Delta$ と 3 値判定 (§4.2)。M6 → 適用規則を `convMethod: 0` 起動区間の 3 条件に (§4.1)。M7 → 衝撃足の再同定 (§4.4) と前 plan の訂正 (「膨張角窓の壁圧」に改名、カウル衝撃足は未測定) |

### 6.2 結果

(未実施)

## 7. 影響範囲

- コード変更なし (既定値は変えない)。`procedures/recommended-settings.md` に適用規則、`methods/convection/theory.md` からリンク。
- ツール: `solver_density_cuda/tools/stage_manifest.py` (`solver` キー)、`case/46.sern_design/cad/diag_wall_cv_budget.py` (docstring)、衝撃足の新抽出スクリプト。
- **Q1 (2D 格子 3 水準) を再開するときの必須要件** (codex plan M2/plan-2 M1・M3、撤回で未実施):
  生産格子の壁折れ線を固定入力として節点を足す格子生成器 (`mesh_sern.py:347,359` が格子から機体形状を決めるため現状不可)、
  実在するブロックだけの解像度表 (`first_top_frac` を含む。`t_base` 未指定では wake ブロックが無い)、
  優先順つき if/elif/else の判定表 + codex plan-2 の反例 2 つを判定器の試験に入れること。

## 8. 完了条件

- [ ] 適用規則を `procedures/recommended-settings.md` に本文化
- [ ] §4.3 の再判定と §4.4 の再同定を §6 の基準で判定し §6.2 に記録、前 plan に結果の訂正行
- [ ] sern-3d §5.1 に委譲の行を追加
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-25` — codex plan 段 2 回目 **NO-GO (C0/M7/m0)** を全件採用。**$C_M$ の生産許容は現行 0.05** (旧値 0.02 を設計と 1 回目レビューが踏襲していた) で flag 差は全量帯内 → **Q1 (格子 3 水準)・Q2 (CFL) を撤回**し、既存系列の再判定・衝撃足の再同定・適用規則だけに縮小。比較は差区間 $I_\Delta$ の 3 値判定に。前 plan の「衝撃足」は x≈0 膨張角を測っていたので訂正した。
- `2026-09-25` — codex plan 段 **NO-GO (C0/M7/m1)** を全件採用して改訂 (§6.1)。**目的を「適用規則 + flag 差の大きさと格子依存」に狭め、正否判定 (旧 Q3 Roe) は保留** (現行 Roe は 2 成分 TP を単成分物性で組むので同一物理でない)。格子 3 水準は `scale` でなく解像度キーで、判定規則・比較式・準定常の扱いを測る前に固定し直した。
- `2026-09-25` — 初稿。前 plan (accepted) の未解決 (#9b・#11・V7・$C_L$/$C_M$) を引き継ぐ。`diagnostician` の判断で**既定化を目的から外し**、
  適用規則の確定と 2D での正否判定に置き換えた (§4.1)。3D の格子収束は sern-3d へ委譲、3D 固定点は判断基準から外した。
