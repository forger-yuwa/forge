# `slauWallNormalChi` の適用規則と、flag 1 が変える量の正否 (2D)

## メタ

- **area**: `convection`
- **status**: `in_progress` (codex plan 3 回目 GO-with-changes を反映)
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
2D 生産・case/16・case/48・周期・軸対称は 0 (周期・軸対称は小規模試験 V6 まで)。新しい構成で 1 にするのは、**診断ツール `diag_wall_cv_budget.py` の質量流束が実カーネルと一致する設定** (§4.1.1 の表を機械的に検査し、
かつ同 dump の 1 step `FORGE_DUMP_MASSFLUX` で対象 CV の全接続面 $\dot m$ がツール値と前 plan V5 の許容内で一致) の**`convMethod: 0` の起動区間**で、
`diag_wall_cv_budget.py` (1 次 SLAU の再計算) により (i) 対象壁 CV (`bcondConfig` の wall physID から作る) の $\rho_w$ が 3 dump 以上単調減少し
$\rho_w/\rho_i<0.1$、(ii) 全接続面の正味流出 $\Sigma\dot m>0$ が 3 dump 以上持続、(iii) 元の $\chi=0$ かつ $\Delta P\neq0$ の壁隣接面があり、
同じ面状態で $\chi\to\chi_n$ に置換した再計算で当該 CV への補充 $\dot m$ が増える (元 $\chi$・$\chi_n$・$\Delta P$・置換前後の $\dot m$ を併記)、
の 3 つを満たしたときだけ。**検査に落ちる設定は「診断不能」で、既知構成 (3D 側壁接続) 以外では 1 にしない** (診断不能 ≠ 不要)。2 次生産場での診断は再構成後の面状態が取れるまで保留。flag は `stage_manifest` の hard キーなので段の途中で切り替えない。」

- `--summary` の床到達数 (全域) と「壁隣接内点 $|\mathbf u|>\sqrt2\hat c$」は**診断を始める兆候**であって適用条件ではない。
- `--wall-normal-chi` は「$\chi$ を $\chi_n$ に置換して再計算する」オプションであって「元の $\chi$ を確認する」ものではない
  (`diag_wall_cv_budget.py:75-77`)。docstring に明記する (§5.1 #5)。

#### 4.1.1 診断可能性の検査 (codex plan-3 M3)

ツールには実カーネルの `slauContactFloor` の質量流束追加 (`convectiveFlux_slau_d.inc.cuh:620`) と `sstEnergyIncludesK` の $p^*$ 差 (:516) が無く、
`slauContactFloor` 0.01 の人工状態でツール +0.099 (流出) / 実カーネル −0.191 (流入) と**符号が逆**になる (codex の再計算)。そこで診断区間の実効値を
`forge_run.log` の起動エコーと `solverConfig.yaml` の両方から取り (食い違えば診断不能、省略 ≡ 既定値に正規化)、次を要求する:

| キー | 要求 | 理由 |
| --- | --- | --- |
| `mesh.discretization` / `mesh.nodeWallDirichlet` | `node` / `1` | 対象構成 |
| `solver` | `SLAU` / `SLAU2` | 質量流束式が同じ (SLAU2 は圧力束のみ差)。Roe・KEEP は不可 |
| `space.convMethod` | `0` (診断区間) | ツールは再構成なし |
| `space.slauContactFloor` | `0` | :620 の追加項がツールに無い |
| `turbulence.sstEnergyIncludesK` | `0` | :516 の $p^*$ 差 |
| `space.lowMachPrecond` | `0` | 散逸スケール $c'$ |
| `space.slauWallNormalChi` | `0` (診断する run 側) | `--wall-normal-chi` は置換予測用 |
| `space.badReconFallback` 等の面状態差し替え | `0` | 面状態が変わる |
| 凝縮・二相 | なし | 面エンタルピー・$P=\rho(1-g)RT$ |
| `physProp.thermalMethod` | ツールの EOS 経路と同じ (CPG / TP を明記) | $\hat c$・$T$ |
| 周期・軸対称 | なし | 適用範囲外 |

- **キー検査に加えて面流束の照合**: 同 dump から 1 step、`FORGE_DUMP_MASSFLUX` で対象 CV の全接続面の $\dot m$ を取り、ツール値と前 plan V5 の許容内で一致すること。
- $\rho_i$ = 対象壁 CV の**壁でない隣接ノード全部の平均** (面選択に依らない)。
- 対応手順: 同一 dump・同一 CV id・同一面 id (`PLANES` index) で、面ごとに `chi, chi_n, ΔP, mdot(chi), mdot(chi_n)` を表にし、
  条件 (iii) は $\Sigma_{\text{面}}[\dot m(\chi_n)-\dot m(\chi)]<0$ (当該 CV への正味補充が増える) を同じ表で判定する。

### 4.2 比較式 — 差区間 (codex plan-2 M4・M5、測る前に固定)

- 各 run の末尾窓 (末尾 40 %) の観測範囲 $I_f=[\min C_f, \max C_f]$ から、**差区間** $I_\Delta=[\min C_1-\max C_0,\ \max C_1-\min C_0]$。
- 判定 (許容帯 $B=[-\text{tol}, +\text{tol}]$): **$I_\Delta\subset B$ → 帯内**、**$I_\Delta\cap B=\emptyset$ → 差が残る**、**それ以外 → 判定不能**。
  $I_\Delta$ は観測窓内の範囲であり、平均値の統計的不確かさ (末尾窓の前半/後半の平均差) とは**別欄**に記録する。
- **壁圧分布どうしの比較 (codex plan-3 M1)**: スカラーの差区間は分布に当てられない (反例 $p_0=(100,102)$, $p_1=(102,100)$ で相対 L2 1.98 %)。
  共通固定窓 $W=[x_f-2t,\ x_f+5t]$ (物理座標、$x_f$ = flag0 末尾 dump 群の中央値、$t$ = 2 mm)、座標は **`MESH/COORD` のノード座標**
  (`CELLS/centCoords` は使わない。`procedures/solver-settings.md` の node の値位置)。両 run の壁圧を $W$ 内の共通格子 (同一メッシュならランプ壁ノード) に
  線形補間、重みは台形則の区間長。末尾 40 % の**全 dump 対** $(i,j)$ で $\ell_{ij}=\lVert p_{1,j}-p_{0,i}\rVert_w/\lVert\bar p_0\rVert_w$。
  **帯内**: $\max\ell_{ij}\le0.01$ / **差が残る**: $\min\ell_{ij}>0.01$ / 他は**判定不能**。位置差 $x_f$ はスカラー差区間 (許容 = 窓内ノード間隔の中央値 1 つ)。
- 旧式 $\lvert\Delta\rvert+a_0+a_1$ (1 回目) と前 plan #10b の式は、非対称な変動で上限にならない (codex の反例: 旧式 0.019、実最大差 0.023)。
  前 plan #10b の $C_T$ と膨張角窓もこの差区間で再計算し、訂正行を足す (§5.1 #2。判定は計算してから書く)。

### 4.3 既存系列の再判定 (新規 run なし。codex plan-3 M4・m5 で改訂)

- 対象: 前 plan の `_v3sern/v3s_flag{0,1}_ext` (72 dump、36000 step)、既存 `ct_flag{0,1}.csv` (`forge_design.metrics.sern_forces.force_history`、`sern_forces.py:85-100`)。

| 列 | 定義 | 判定 | 相対 drift / osc (`--osc` は $(\max-\min)/\lvert\text{mean}\rvert$ の全幅) | 出所 | 絶対換算 (末尾平均) |
| --- | --- | --- | --- | --- | --- |
| `C_T` | 圧力のみ $F_{gross}/F_{ideal}$ | 帯判定 (許容 0.002) | 0.0002 / 0.0005 | 前 plan V3 から継承 (1/5 則ではない) | 0.000186 / 0.000465 @0.929 |
| `C_T_with_shear` | $(F_{gross}+F_{x,\tau})/F_{ideal}$ | **帯判定 (許容 0.002)** — sern-3d に $C_T$ 内・摩擦込み外の実例 | 0.0002 / 0.0005 | 同上 | 同上 |
| `C_L` / `C_L_with_shear` | $F_y/F_{ideal}$ / 摩擦込み | 帯判定 (0.002) | 0.0014 / 0.0014 | 0.002/0.277/5 | 0.00039 |
| `C_M` | $M_{noseup}/(F_{ideal}H)$ | 帯判定 (0.05) | 0.0014 / 0.0014 | 0.05/7.16/5 | 0.010 |
| `C_T_wall`, `C_T_friction` | 内訳 | 報告のみ | — | — | — |

- **「生産許容内」= 帯判定 5 列すべてが帯内**。1 列でも「差が残る」→ 差が残る、他 → 判定不能 (最悪列の結論)。
- `--tail 0.4`。**`DRIFTING` / `TRANSIENT-UNSETTLED` → 「判定保留」** (ラベルを変えない。延長 run はしない)。`STEADY` / `OSCILLATING` の列だけ §4.2 の差区間で判定。
  `OSCILLATING` は比較禁止ではなく、平均 ± span/2 を報告して差区間へ。

### 4.4 衝撃足の再同定 (codex plan-2 M7・plan-3 M2、測る前に固定)

旧抽出 (`v3sern_foot_perrun.py:42`, `v3sern_series.py:39-40`) は $\max|dp/dx|$ をランプ全域から取り圧縮・膨張を区別しない。`v3sern_series.py:29` の
`CELLS/centCoords` も node の値位置と違う。新しい抽出は**結果を見る前に**書いて commit し、座標は `MESH/COORD` のノード座標を使う。

1. **探索区間**: ランプ壁ノードのうち $x\in[x_{lip},\ L_{ramp}]$ ($x_{lip}$ = カウル後縁の x)。
2. **候補**: $dp_w/dx>0$ の局所最大すべて。各候補 $x_c$ で上昇比 $r=p(x_c+2t)/p(x_c-2t)$ (壁ノード線形補間)。**登録条件 $r\ge1.05$**
   (1.1 は 8 % 上昇を落とす。codex の反例)。
3. **唇由来の識別**: 候補 $(x_c, y_w)$ から唇 $(x_{lip}, y_{lip})$ への線分上に 50 点を取り、場の $\nabla p$ を線分の法線方向 (壁側→外側) に射影した符号が
   **80 % 以上の点で同符号** (連続した圧縮帯) なら「唇から続く圧縮構造」と同定。満たさない候補は外す。
4. **分岐表 (上から優先)**:

| 条件 | 結論 |
| --- | --- |
| flag0・flag1 とも識別済み候補が 1 つで $\lvert x_c^{(1)}-x_c^{(0)}\rvert\le$ 局所間隔 | 比較へ (§4.2 の分布比較と位置差) |
| 複数候補 | flag0 で $r$ 最大の候補を採る。flag1 の $r$ 最大候補と位置が 1 間隔超で違う → **判定不能 (候補不一致)** |
| 窓 $W$ がランプ端を越える | 端で切り、残りが窓長の 50 % 未満なら **判定保留** |
| 片側のみ検出 | **判定不能** (両分布を図で残す。「差が残る」とは書かない) |
| m6_on で両側未検出 | 「**登録条件では同定できない (m6_on)**」→ m4_off の flag 0/1 対 (2D 生産レシピ、flag 0 から分岐、各 1 本、AWS) で 1–3 を再実行 |
| m4_off でも未検出 | **判定保留 (本手順では衝撃足を同定できず)**。「衝撃が当たらない」とは書かない |

- 各 run の準定常は固定座標上の系列 ($x_f$、窓平均壁圧、$x_f-t$ / $x_f+t$ / $x_f+3t$ の 3 点) に `--drift 0.002 --osc 0.005` で別途判定。

## 5. 実装ステップ

1. plan 縮小 + 前 plan の訂正 + codex plan 3 回目。
2. §4.3 の再判定と §4.4 の再同定 (既存 dump。必要なら m4_off の 1 組)。
3. `stage_manifest` の `solver` キー、`diag_wall_cv_budget.py` の docstring。
4. 規則の本文化、codex result、accepted。

### 5.1 残作業 (優先順)

**計算資源は AWS** (ユーザ指示 2026-09-25: ローカルで大きな計算をかけない。AWS は他セッションと共有なので起動前に `aws_instance.sh status`/`busy`、手動 stop はしない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| ~~1~~ (**済 2026-09-25**: plan 3 回目 GO-with-changes C0/M4/m1、全件採用して §4 に反映) | **plan 縮小 + codex plan 3 回目** | 本改訂。Critical/Major は `diagnostician` に諮る | F |
| 2 | 既存系列の再判定 (§4.3) + 前 plan #10b の差区間での再計算 | `_v3sern/ct_flag{0,1}.csv` (AWS から csv だけ取得) の **5 列 (`C_T`, `C_T_with_shear`, `C_L`, `C_L_with_shear`, `C_M`)** に §4.3 の閾値表で `check_quasisteady.py --series-csv`、`STEADY`/`OSCILLATING` の量だけ差区間 (§4.2)。前 plan #10b の $C_T$・膨張角窓も差区間で再計算し訂正行 | O |
| 3 | 衝撃足の再同定 (§4.4) | 新スクリプト (§4.4 の手順 1–3・分岐表、座標は `MESH/COORD`、分布比較は §4.2 の全 dump 対 $\ell_{ij}$) を**結果を見る前に**書いて commit、既存 72 dump に当てる (AWS)。m6_on で両側未検出なら m4_off の 1 組 (AWS ≈ 30 min)。結論は `diagnostician` | O (結論 F) |
| 3b | 診断可能性の検査スクリプト (§4.1.1) | `case/46.sern_design/cad/diag_applicability.py`: 表のキー検査 (起動エコーと yaml の両方、省略 ≡ 既定値) + `FORGE_DUMP_MASSFLUX` 1 step 照合。合格 = 全キー OK かつ対象 CV 全接続面の $\lvert\dot m_{tool}-\dot m_{kernel}\rvert$ が前 plan V5 の許容内。試験: 既知の 3D 側壁接続 (`run_0437/res_0` 相当、AWS) で PASS、`slauContactFloor` 0.01 の config で「診断不能」を返すこと | O |
| 4 | 前 plan の訂正行 | **済 2026-09-25**: V3 表・§6.2・#14・受入の範囲・§9 を「膨張角窓の壁圧」に改名し、カウル衝撃足を「未測定」に | O |
| 5 | ツール | `stage_manifest.py` に `solver` hard キー + 単体試験 (SLAU↔ROE が別区間・省略 ≡ 既定 SLAU)、`diag_wall_cv_budget.py` の docstring に `--wall-normal-chi` の意味。合格: 単体試験 PASS、`test_gate_bad_input.py` PASS | O |
| 6 | 規則の本文化 + codex result → accepted | §4.1 の文を `procedures/recommended-settings.md` に、`methods/convection/theory.md` の対策節からリンク | F |
| 7 | sern-3d への委譲 | `tooling-nozzle-sern-3d.md` §5.1 に「R5n を `slauWallNormalChi: 1` で再取得 (メッシュ生成コマンド・blocksize 記録)。3D 固定点は未解決」を 1 行 | F |

## 6. 検証

**判定基準はすべて測る前に固定する**。**新設の許容差は無い** (生産許容は sern-3d §8 の現行値、衝撃足は前 plan V3)。

| 項目 | 量 | 判定 (測る前に固定) |
| --- | --- | --- |
| 既存系列 (§4.3) | 各 run の `C_T`, `C_T_with_shear`, `C_L`, `C_L_with_shear`, `C_M` (帯判定 5 列、全列帯内で「生産許容内」) | 量別閾値で `STEADY`/`OSCILLATING` の量だけ差区間 $I_\Delta$ を判定: $I_\Delta\subset B$ → 帯内 / $I_\Delta\cap B=\emptyset$ → 差が残る / それ以外 → 判定不能。`DRIFTING`/`TRANSIENT-UNSETTLED` は「判定保留」 |
| 衝撃足 (§4.4) | 各 run の $x_f$・窓内壁圧、全 dump 対の $\ell_{ij}$、$x_f$ 差 | §4.4 の分岐表の結論に従う。比較に進んだ場合、各 run が `STEADY`/`OSCILLATING` のときだけ: $\max\ell_{ij}\le0.01$ → 帯内 / $\min\ell_{ij}>0.01$ → 差が残る / 他は判定不能。$x_f$ 差はスカラー差区間 (許容 = 局所間隔 1 つ) |
| 診断可能性 (§4.1.1) | キー検査 + 面流束照合 | 合格 = 全キー OK かつ照合が V5 許容内。既知構成で PASS・`slauContactFloor` 0.01 で「診断不能」 |
| ツール | `stage_manifest` 単体試験 | PASS |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-25` | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用** (却下なし、2026-09-25 `diagnostician` 判断。根拠は各行をコードで確認)。M1 → Q3 (Roe) を保留、目的を「正否は決めない」に (§1・§2・§4.1)。M2 → `scale` をやめ解像度キーで 3 水準・形状不変の記録・IC 手順 (§4.3)。M3 → 「2D 生産は収束」を訂正 (§3)、準定常を全 run・全量に、未定常を言い換えない (§4.5)。M4 → Q2 を $D_{intra}$ 先行・フラグごと独立・比を削除、`v7_dist.py` 引数化 (§4.6、§5.1 #2)。M5 → Q1 を「格子間変化とフラグ差を別評価・3 結論」に (§6)。M6 → 兆候と適用 3 条件を分離、適用範囲を既知 3D 構成に限定 (§4.2)。M7 → 上限 $U=\lvert\Delta\rvert+a_0+a_1$ と平均の安定性 (§4.4)、前 plan #10b に訂正行。m8 → `stage_manifest` に `solver` hard キー (§5.1 #2) |
| plan | `2026-09-25` (2 回目) | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-2.md) | **NO-GO**, C0/M7/m0 | **全件採用** (2026-09-25 `diagnostician` 判断)。**M2 で前提が変わった** ($C_M$ 許容は現行 0.05、flag 差は全量帯内) ため **Q1/Q2 を撤回して plan を縮小** (§2)。M1・M3 → Q1 撤回で解消し、再開時の必須要件として §7 に残す。M4・M5 → 差区間 $I_\Delta$ と 3 値判定 (§4.2)。M6 → 適用規則を `convMethod: 0` 起動区間の 3 条件に (§4.1)。M7 → 衝撃足の再同定 (§4.4) と前 plan の訂正 (「膨張角窓の壁圧」に改名、カウル衝撃足は未測定) |
| plan | `2026-09-25` (3 回目) | [2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-3.md](../../notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan-3.md) | **GO-with-changes**, C0/M4/m1 (縮小と「膨張角窓」訂正は支持) | **全件採用** (2026-09-25 `diagnostician` 判断)。M1 → 分布比較は全 dump 対 $\ell_{ij}$ の 3 値判定、座標は `MESH/COORD` (§4.2)。M2 → 衝撃足の同定手順 (候補 $r\ge1.05$・唇由来の識別) と分岐表、未検出は「同定できない」で「当たらない」とは書かない (§4.4)。M3 → 規則を「ツールと実カーネルが一致する設定のときだけ」に、診断可能性の検査表と 1 step 面流束照合 (§4.1.1、§5.1 #3b)。M4 → `C_T_with_shear`・`C_L_with_shear` を帯判定に追加し列定義表 (§4.3)。m5 → 閾値の出所と絶対換算を表に (§4.3) |

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

- `2026-09-25` — codex plan 段 3 回目 **GO-with-changes (C0/M4/m1)** を全件採用して §4–§6 に反映。分布比較の式、衝撃足の同定手順と分岐表、診断可能性の検査 (設定表 + 1 step 面流束照合)、摩擦込み係数の追加。`in_progress` へ。
- `2026-09-25` — codex plan 段 2 回目 **NO-GO (C0/M7/m0)** を全件採用。**$C_M$ の生産許容は現行 0.05** (旧値 0.02 を設計と 1 回目レビューが踏襲していた) で flag 差は全量帯内 → **Q1 (格子 3 水準)・Q2 (CFL) を撤回**し、既存系列の再判定・衝撃足の再同定・適用規則だけに縮小。比較は差区間 $I_\Delta$ の 3 値判定に。前 plan の「衝撃足」は x≈0 膨張角を測っていたので訂正した。
- `2026-09-25` — codex plan 段 **NO-GO (C0/M7/m1)** を全件採用して改訂 (§6.1)。**目的を「適用規則 + flag 差の大きさと格子依存」に狭め、正否判定 (旧 Q3 Roe) は保留** (現行 Roe は 2 成分 TP を単成分物性で組むので同一物理でない)。格子 3 水準は `scale` でなく解像度キーで、判定規則・比較式・準定常の扱いを測る前に固定し直した。
- `2026-09-25` — 初稿。前 plan (accepted) の未解決 (#9b・#11・V7・$C_L$/$C_M$) を引き継ぐ。`diagnostician` の判断で**既定化を目的から外し**、
  適用規則の確定と 2D での正否判定に置き換えた (§4.1)。3D の格子収束は sern-3d へ委譲、3D 固定点は判断基準から外した。
