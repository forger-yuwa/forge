# `slauWallNormalChi` の適用規則と、flag 1 が変える量の正否 (2D)

## メタ

- **area**: `convection`
- **status**: `draft`
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (「既知の限界」「対策 (opt-in)」節)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) (適用規則を書く先)
- **related_plans**:
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (前 plan。§5.1 #9b・#11、§6 V7、§6.2「未解決・適用限界」を引き継ぐ)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (3D 生産の格子収束列 R5n。3D の格子収束はこちらへ委譲)
- **created**: `2026-09-25`
- **owner**: `CFD Dev`

## 1. 目的

`space.slauWallNormalChi` を**いつ 1 にするか**の規則を決め、flag 1 が動かす量 ($C_L$ / $C_M$ / 衝撃足の壁圧) の
**正否を 2D で独立に判定する**。**既定値は変えない** (既定 0 のまま。§4.1)。

## 2. スコープ

- **やる**: 2D SERN 生産構成での (Q1) 差の格子依存、(Q2) CFL 独立性の両フラグ比較、(Q3) 必要なときだけ Roe を第三者にした正否、
  (Q4) 各 run の $C_L$/$C_M$ 準定常。結果から `procedures/recommended-settings.md` の適用規則を確定する。
- **やらない**:
  - **既定化** (前 plan §5.1 #11 の「既定化の判断」は、本 plan で「既定 0 のまま・規則で使う」に読み替えて閉じる。§4.1)。
  - **3D 接続模型の格子収束 (前 plan #9b)**: 3D 生産は flag 0 で解を持たず常に flag 1 なので、決めるべきは「flag 1 の生産解の格子収束」=
    sern-3d の R5n 型の再取得。[`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 へ委譲する (§5.1 #5)。
  - **3D の固定点 (前 plan V7 のドリフト減衰)**: flag 0 の解が無いのでフラグに帰属できない量。判断基準から外し、2D の Q2 で置き換える。
  - **SU2 を第三者にすること**: frozen_tp 2 成分・等温壁 TP を SU2 で同一条件にできず、[`su2-cross-check.md`](../../procedures/su2-cross-check.md) の
    前提「同一 BC」を満たせない。

## 3. 関連 docs と前提

- 前 plan の受入 (opt-in として、試験したケース・設定・区間で回帰許容内) と、その「未解決・適用限界」(§6.2)。
- 2D SERN 生産の起動レシピと収束判定: [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) §4.12 / §4.16 / §4.16.1、skill `sern-eval`。
- 許容差は新設しない: 力係数は sern-3d §8 の R5n 格子収束許容 ($|\Delta C_T|\le0.002$, $|\Delta C_L|\le0.002$, $|\Delta C_M|\le0.02$、絶対値)、
  CFL 独立性は前 plan V7 の $\varepsilon$ (P 0.02 %、ro/roU/roe 0.2 %、壁隣接集合) を借りる。
- 前 plan の実測: 2D 生産 (設計点 m6_on) で flag 1 − flag 0 = $C_L$ +0.437 % / $C_M$ −0.320 % (絶対 0.0012 / 0.023)、
  flag0 の $C_L$ は `DRIFTING`。**$C_M$ の差 0.023 は R5n 許容 0.02 を超える**。

## 4. 設計方針

### 4.1 既定化はしない (2026-09-25 `diagnostician` 判断)

- flag 1 が**必要**な構成は、現状 1 つだけ: node + `nodeWallDirichlet: 1` で**側壁∩後端面の壁 CV が排出される 3D 接続構成**
  (flag 0 は 3 格子とも NaN、`case/46.sern_design/README.md` の run_0443–0446)。
- 2D 生産は flag 0 で収束し、flag 1 は $C_L$/$C_M$ を**正否不明のまま**動かす。既定 ON は既存 node 生産 (2D SERN・case/16・case/48) の
  答えを黙って変える操作で、利得が無い。
- 規則 (opt-in を明示して使う) の方が監査できる: `stage_manifest` の hard キー (前 plan #12) と `RUN_PROVENANCE` にフラグが残る。

**適用規則の草案** (§5.1 #1 で「暫定」として `recommended-settings.md` に書き、Q1–Q3 の結果で確定する):
「node + `nodeWallDirichlet: 1` で**壁 CV の排出**が出る構成 (`diag_wall_cv_budget.py --summary` で床到達 > 0、または壁隣接内点
$|\mathbf u|>\sqrt2\,\hat c$ の面がある。現状は 3D 側壁接続) では `space.slauWallNormalChi: 1`。**2D 生産は 0 のまま** (Q1 の結果で更新)。」

### 4.2 Q1 差の格子依存 (2D SERN 生産)

- 構成: 設計点 m6_on、frozen_tp、`problem_moo_frozen_tp_cycle3op.yaml`。**3 水準** (生産 / 壁法線・流れ方向 ×1/√2 / ×1/2) を
  メッシャの scale で作り、**生成コマンドを run に記録**する。各水準 flag 0/1。同一起動レシピ (chain §4.16)、`FORGE_CUDA_BLOCKSIZE` 同一。
- 量: $C_T$, $C_L$, $C_M$ (`forge_design.metrics.sern_forces.force_history`、自前で組まない)、衝撃足 L2 (前 plan のノルム)。
- $\Delta_k = C^{(1)}_k - C^{(0)}_k$ を**末尾平均の差** (#10b 式: $|\Delta\text{末尾平均}|+\text{振幅}$) で取る。

### 4.3 Q2 CFL 独立性 (2D、両フラグ同型)

- 生産格子で `cfl_main` 5 と 2.5 を **ΣCFL を揃えて**各フラグ 1 対。前 plan V7 と同じ $D_{inter}$ を両フラグで取る。
- 3D V7 のドリフトは flag 0 の解が無いので flag に帰属できず、判断基準から外す (§2)。

### 4.4 Q3 正否の第三者 (Q1 が (b) のときだけ)

- 同じ 2D 最細格子で **Roe** (forge 内) を flag 0 の解を起点に同じ本段長で回す。Roe は壁 CV の補充で符号が反転する
  (前 plan §1: 「同じ場で流束を Roe に替えると 1 step 目から符号が反転する」) ので、$\chi$ の定義に依らない第三者になる。

### 4.5 Q4 $C_L$/$C_M$ の準定常

- Q1 の全 run に掛ける。閾値は R5n 許容の 1/5 を相対化したもの (§6 の表)。STEADY にならなければ本段を 1 回だけ倍にし、
  それでも駄目なら「限界サイクル」と書いて平均±振幅で進む。

## 5. 実装ステップ

1. `procedures/recommended-settings.md` に適用規則の草案を「暫定」で書く (§4.1)。
2. 2D SERN のメッシャで 3 水準を作り、生成コマンドと `check_mesh_quality` の VERDICT を記録する (`case/46.sern_design/`)。
3. Q1 の 6 run、Q2 の 4 run を AWS で回し、run-watcher で post-check。
4. Q1 が (b) なら Q3 (Roe) を 1–2 run。
5. 結果を §6.2 に書き、規則を確定する。sern-3d §5.1 に R5n の flag 1 再取得を 1 行入れる。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 0 | **codex plan 段** | `codex_review.py <この plan> --stage plan`。Critical/Major の採否は `diagnostician` に諮る | F |
| 1 | 規則の草案を「暫定」で書く | `procedures/recommended-settings.md` に §4.1 の草案 + 「Q1–Q3 の結果で更新」 | F |
| 2 | **Q1: 2D 3 水準 × 2 flag (6 run)** | 触る: `case/46.sern_design/` (run 生成は既存の SERN 2D runner)。メッシュ 3 水準の生成コマンドを記録・`check_mesh_quality` PASS。各 run `check_convergence --segment` + Q4 の準定常。判定は §6 の (a)/(b)/(c)。**計算は AWS** (ユーザ指示: ローカルで大きな計算をかけない)。見積もり: 1 run 3–10 min (ローカル実績 463 s/3 作動点から)、計 ≈ 1 h | O |
| 3 | **Q2: CFL 対 (4 run)** | 生産格子、`cfl_main` 5 / 2.5、ΣCFL を揃える。V7 の $\varepsilon$ で列ごと判定 (§6)。AWS ≈ 40 min | O |
| 4 | Q3: Roe (Q1 が (b) のときのみ、1–2 run) | 最細格子、flag 0 の解を起点。符号と距離で判定 (§6)。Roe が起動しなければ F に戻す。AWS ≈ 20 min | O (判定 F) |
| 5 | 3D の格子収束を sern-3d へ委譲 | `tooling-nozzle-sern-3d.md` §5.1 に「R5n を `slauWallNormalChi: 1` で再取得 (メッシュ生成コマンド・blocksize 記録)」を 1 行入れ、本 plan からリンク。参考費用: 粗 ≈ 12 min / 中 ≈ 25 min / 採用 245 万節点 ≈ 50 min (66k step、40–50 ms/step の外挿)、計 ≈ 1.5–2 h + メッシュ生成 | F |
| 6 | 規則の確定 + codex result 段 | 結果を §6.2 に、規則を `recommended-settings.md` に本文化。`stage_manifest`/`RUN_PROVENANCE` にフラグが残ることを明記 | F |
| 7 (任意) | 前 plan V7 の窓 B | AWS を他用途で起動したときだけ。`run_0447` の終端から +10000。$D_{intra}(B)/D_{intra}(A)$ ≤ 0.7 → 減衰する過渡 / 0.7–1.3 → 持続ドリフト (ケース固有、flag に帰属しない) / それ以外 → 判定不能。**どの結果でも規則は変えない** (記録のみ) | O |

## 6. 検証

**判定基準はすべて測る前に固定する** (2026-09-25 `diagnostician`)。**新設の許容差は無い**。

| 項目 | 量 | 合格・分岐 (測る前に固定) |
| --- | --- | --- |
| Q1 差の格子依存 | 3 水準それぞれの $\Delta C_T, \Delta C_L, \Delta C_M$ (末尾平均の差、#10b 式) と衝撃足 L2 | **(a)** $\lvert\Delta\rvert$ が細分で単調減少し、最細で R5n 許容 ($\lvert\Delta C_T\rvert\le0.002$, $\lvert\Delta C_L\rvert\le0.002$, $\lvert\Delta C_M\rvert\le0.02$) 以内 → 「差は離散化誤差の範囲。2D は 0 のまま (安定性で選ぶ)」。**(b)** $\lvert\Delta_{\text{fine}}\rvert/\lvert\Delta_{\text{coarse}}\rvert>0.7$ かつ許容超え → 「flag が収束解を変える」→ Q3 へ。**(c)** 単調でない → 判定不能、水準を 1 つ足す |
| Q2 CFL 独立性 | 各フラグの $D_{inter}$ (前 plan V7 と同定義、ΣCFL 揃え) | 列ごと $D_{inter}(\text{flag})\le\varepsilon$ かつ $D_{inter}(1)\le2\,D_{inter}(0)$ → 「$\chi_n$ は CFL 依存を足していない」 |
| Q3 正否の第三者 | Roe・flag 0・flag 1 の $C_L$, $C_M$ | $\mathrm{sign}(C_{Roe}-C^{(0)})=\mathrm{sign}(C^{(1)}-C^{(0)})$ かつ $\lvert C_{Roe}-C^{(1)}\rvert<\lvert C_{Roe}-C^{(0)}\rvert$ ($C_L$, $C_M$ とも) → flag 1 を支持。片方でも逆 → 支持しない (規則は「安定性が要るときだけ 1」に留める) |
| Q4 準定常 | Q1 各 run の $C_L$, $C_M$ 系列 | `check_quasisteady.py --series-csv`: $C_L$ `--drift 0.0014 --osc 0.0014`、$C_M$ `--drift 0.00056 --osc 0.00056`、`--tail 0.4`。OSCILLATING は平均±振幅。STEADY にならなければ本段を 1 回だけ倍、それでも駄目なら「限界サイクル」と明記 |
| 収束 | 全 run | `check_convergence.py --segment` の VERDICT を同一設定区間で併記 (本段) |

- **参照値の出所**: Q1 は同格子の flag 0 (外部参照なし)、Q2 は前 plan V7 の $\varepsilon$、Q3 は Roe、許容は sern-3d §8 の R5n 値。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

### 6.2 結果

(未実施)

## 7. 影響範囲

- コード変更なし (既定値は変えない)。`procedures/recommended-settings.md` に適用規則を追加。
- `methods/convection/theory.md` の「対策 (opt-in)」節に規則へのリンクを足す。
- 既存ケース・手順への影響なし (既定 0 のまま)。

## 8. 完了条件

- [ ] `procedures/recommended-settings.md` に適用規則を本文化 (暫定を外す)
- [ ] Q1・Q2 (と必要なら Q3) を §6 の基準で判定し §6.2 に記録
- [ ] sern-3d §5.1 に R5n の flag 1 再取得を追加 (委譲)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-25` — 初稿。前 plan (accepted) の未解決 (#9b・#11・V7・$C_L$/$C_M$) を引き継ぐ。`diagnostician` の判断で**既定化を目的から外し**、
  適用規則の確定と 2D での正否判定に置き換えた (§4.1)。3D の格子収束は sern-3d へ委譲、3D 固定点は判断基準から外した。
