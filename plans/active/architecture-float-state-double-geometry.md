# 状態を float にして、差を取る幾何の量を double の座標から作る

## メタ

- **area**: `architecture`
- **status**: `draft`
- **related_docs**:
  - `methods/architecture/overview.md` (§6.2a 幾何の量の精度)
  - `methods/time_integration/implementation.md` (commit の丸め、`qAccumulatorFP64`)
- **related_plans**: [`time_integration-line-implicit-speed`](time_integration-line-implicit-speed.md) §5.1 #23・§6.21・§6.22 (発端)、[`tooling-nozzle-isothermal-wall-chain`](tooling-nozzle-isothermal-wall-chain.md) (冷却壁の格子と FP64 のビルド)
- **created**: `2026-10-10`
- **owner**: Claude (Opus 5.5)

## 1. 目的

case/45 の冷却壁の M6 ノズル (第一層厚 / 半径 ≈ 1e-7〜1e-6) は、float の座標では第一層が数 ulp にしかならないので、FP64 のビルドで回している。
float のビルドは同じ構成で 1 step が FP64 の 0.64 倍 (20.5 対 32.0 ms、line-implicit-speed §6.21) だが、残差は約 10 倍で使えない。
座標の差を取る量を、読み込み時に double の座標から作って float で渡す。これで状態は float のまま、θ_r・Q_w の到達を FP64 と同等に保ち、1 step を 3 割ほど短くする。

## 2. スコープ

- **やる**:
  - メッシュの座標・面重心を double でも読み込み (ホスト)、面ごとの差の量 (節点間ベクトル、拡散の幾何係数、再構成のベクトル) を double で作って `flow_float` の配列に入れる。
  - 毎 step 絶対座標の差を取っているカーネルを、それらの配列を読むように書き換える (調査 [`notes/investigations/2026-10-10-geometry-precision-inventory.md`](../../notes/investigations/2026-10-10-geometry-precision-inventory.md) の §1)。
  - 読み込み時に作る量 (LSQ の係数、軸対称の closure、node の壁関数の代表距離、弱形式の等温壁の d1) を double の座標から作る。
  - case/45 で float の到達を FP64 と比べ、標準ケースの回帰を確かめる。
- **やらない**:
  - 変換器の float 経路の修正。冷却壁の格子は従来どおり FP64 の変換器 (`stod`・17 桁) で作る。
  - `geom_float` と `flow_float` を別の型にする (調査 §0-1 のとおり、配列の表とカーネルの引数の作り直しになる)。
  - 状態の commit の精度 (`qAccumulatorFP64` の軸対称対応)。§6 の測定で丸めによる停滞が見えたときに別に立てる (§5.1 #11)。
  - cell 離散化の正しさの確認 (使っていない)。ただし cell の経路もビルドと既存の回帰は壊さない。

## 3. 関連 docs と前提

- 現状の精度の扱い: `methods/architecture/overview.md` §6.2a。
- 調査: `notes/investigations/2026-10-10-geometry-precision-inventory.md`。
- 冷却壁の格子と FP64 の経緯: メモ [cooled-wall-mesh-precision]、`plans/active/tooling-nozzle-isothermal-wall-chain.md`。
- 発端の測定: `plans/active/time_integration-line-implicit-speed.md` §6.21 (float 20.5 対 FP64 32.0 ms/step、float の残差は約 10 倍) と §6.22 (codex の諮問。「混成 typedef では試験にならない、差を double で作る案を先に」)。
- B0 の水準での 1 step の変化と float の刻み (2026-10-10、`case/45.isobutane_m6_d155/run_0387_dqulp_B0`、FP64 で B0 の最終状態から 10 step、9 → 10 step の変化):
  - |Δq| が float32 の刻みの半分に満たない節点は、ρ 17.7 %、ρu 8.9 %、ρv 3.5 %、ρE 11.9 %、ρk 3.4 %、ρω 5.9 %。
  - 比の中央値は 2〜31 (刻みより大きく揺れている)。

## 4. 設計方針 (案、上位の判断の前)

### 4.1 どこを直すか

調査 §5 の順。**毎 step カーネルが絶対座標の差を作っている箇所**が本体で、読み込み時の量はその次。

### 4.2 差の量の作り方

1. 読み込み (`mesh.cpp`): `/MESH/COORD`・`/PLANES/centCoords`・`/CELLS/centCoords` を、`geom_float` の配列とは別に **double のホスト配列**にも読む。HDF5 が double ならそのまま、float なら float から広げる (精度は上がらない)。node 離散化の中心 (= 節点座標) も double で持つ。
2. 面ごとの量を double で作り、`flow_float` の面の配列 (`p`/`p_d`) に入れる:
   - `ex, ey, ez` = cc[ic1] − cc[ic0] (節点間ベクトル。node は節点の差)
   - `elen` = |e|
   - `ediff` = δ/|e| = |S|² / |e·S| (拡散の幾何係数。粘性の対角と transport_diag の `geo` もこれ)
   - `r0x..r0z`・`r1x..r1z` = pc − cc[ic0]、pc − cc[ic1] (辺中点を使わない再構成とリミタ用。node は辺中点なので 0.5·e で足りる経路が多い。どの経路が読むかは実装のときに洗い出す)
   - 境界の半割面 (ic1 がゴースト) の差: node はゴーストを使わない経路が多い。cell のゴーストの中心は double で作る。
3. 読み込み時に作る量:
   - LSQ の係数 (`gradLSQ` 1/2): d を double の座標から作る。
   - 軸対称の closure `A_closure = Σ±S_f r_f`: r は double の面重心から作る。
   - node の壁関数の代表距離 y = −(cc_I − cc_W)·n: 壁節点ごとに double で作って配列にする。
   - 弱形式の等温壁の d1/d2: double の座標から作る (今は丸めた座標の double の差)。
   - `delta_les` (DES): double の差から作る。
4. 書き換えるカーネル (調査 §1c〜§1j): 粘性 (`viscousFlux_d`)、スカラー・化学種・受動種の拡散、再構成 (Roe・HLLE・SLAU・KEEP と include される legacy)、リミタ (fused と periodic)、スカラー DPLUR・block DPLUR・前処理版の粘性の対角、壁関数 (`ransWallFunction_d`・`wmlesWallModel_d`)。読むのは e・|e|・δ/|e|・r で、座標の差の式は消す。

### 4.3 既定の扱い (上位に諮る点)

- 案 A (切り替えなし、常に新しい経路):
  - FP64 のビルドでは、差を double で作るのは今も同じなので、結果は丸めの順序の違いの範囲で一致するはず (§6 V2 で確かめる)。
  - float のビルドでは結果が変わる (精度が上がる方向)。標準ケースの回帰が変わる。
  - コードの経路は 1 本で済む。
- 案 B (`mesh.geometryDifferences: 1` の opt-in):
  - 既存の float の回帰は変わらない。
  - カーネル 20 余りに経路の分岐が要る。
- 呼び出し側の推奨は案 A。ユーザの方針は「標準にしたものは既定にする」(メモ [new-standard-becomes-default])。

### 4.4 メモリと速さ

- 面ごとに float 5〜11 個 (e 3・|e| 1・δ/|e| 1・r 6)。2D の node では面 ≈ 2 × 節点なので、+40〜90 B/節点。座標 6 個の読み込みが 3〜5 個に置き換わるので、帯域はほぼ変わらない見込み。
- 速さの目標: float の新しい経路で、FP64 の 1 step の 0.70 倍以下 (32.0 → 22.4 ms 以下)。

### 4.5 状態の commit の精度

- float の保存量では |dq| < ½ ULP の更新は消える。B0 の水準では、刻みの半分に満たない節点は量により 3〜18 %。比の中央値は 2〜31 (§3)。
- 場全体が止まる心配は小さいが、ゆっくりした漂い (θ_r のドリフト) が揺れに埋もれて進むかは分からない。§6 V4 で |dq|/ULP の統計と到達を見て、停滞が見えたら §5.1 #11 を立てる。

## 5. 実装ステップ

1. 読み込み: double のホスト配列と面の差の量の作成 (`mesh.cpp`・`mesh.hpp`・`variables.cpp`・`variables.hpp`)。
2. カーネルの書き換え: 粘性・拡散・陰解法の対角 (`viscousFlux_d.cu`・`scalarTransport_d.cu`・`speciesTransport_d.cu`・`passiveKernels_d.cuh`・`passiveFct_d.cuh`・`timeIntegration_d.cu`)。
3. カーネルの書き換え: 再構成・リミタ (`convection/*.inc.cuh`・`limiter_d.cu`・`limiterPeriodic_d.cuh`・`passiveLimiter_d.cuh`)。
4. 読み込み時の量: LSQ (`calcGradient_d.cu`)、closure (`variables.cpp`)、壁関数の代表距離 (`ransWallFunction_d.cu`・`wmlesWallModel_d.cu`)、d1 (`conjugateWall.cpp`)、`delta_les`。
5. 検証 (§6)、文書の更新 (`methods/architecture/overview.md` §6.2a、`procedures/` の FP64 の手順)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 設計の諮問 | §4 (特に 4.2 の量の選び方、4.3 の案 A/B) と §6 の基準を上位に諮る | F |
| 2 | codex plan 段 | §4・§6 が固まったら `codex_review.py --stage plan` | F |
| 3 | 読み込みと差の量 | §5 の 1。合格: FP64 のビルドで、作った e・δ/\|e\| が今のカーネルの式の値と一致 (相対 1e-14 以内) | O |
| 4 | 粘性・拡散・対角 | §5 の 2。合格: §6 V2 | O |
| 5 | 再構成・リミタ | §5 の 3。合格: §6 V2 | O |
| 6 | 読み込み時の量 | §5 の 4。合格: §6 V1・V2 | O |
| 7 | FP64 の同一性 | §6 V2 | O |
| 8 | float の幾何と到達 | §6 V1・V3・V4 | O |
| 9 | 速さ | §6 V5 | O |
| 10 | 標準ケースの回帰 | §6 V6 | O |
| 11 | 状態の commit の精度 (条件付き) | V4 で停滞が見えたら、`qAccumulatorFP64` の軸対称対応を別 plan に | F |

## 6. 検証 (案、上位の判断の前に事前登録の形にする)

- **V0 ビルド**: float と FP64 の両方のビルドが通る。
- **V1 幾何の係数** (float のビルド): case/45 の格子で、作った e・\|e\|・δ/\|e\| を FP64 のビルドの値と比べる。
  - 全面で相対差 ≤ 1e-6 (float の丸めの範囲)。特に壁の第一層の辺を見る。
  - 今の float の経路 (丸めた座標の差) の誤差も同じ物差しで記録する (2026-10-08 の実測では最大 28 %)。
- **V2 FP64 の同一性** (FP64 のビルド): 新しい経路と今のバイナリ (lineM_fp64) で、同じ状態 (B0 の最終状態) の 1 step の全残差を比べる。
  - 相対差 ≤ 1e-12 (丸めの順序の違いの範囲)。20 step の軌道も比べる。
- **V3 同じ状態の残差** (float 対 FP64): 同じ状態から 1 step の全残差の場。
  - 新しい float と今の float の両方を、FP64 との距離で比べる。新しい float のほうが近いこと (メモ [float-regress-double-truth] の物差し)。壁際の行を分けて見る。
- **V4 到達** (float の新しい経路、case/45 の B0 の構成):
  - `run_0183` の res_100000 から、水準 (2 出力連続、2500 step ごと、最大 20 万 step) まで。
  - 合格: 到達時の θ_r(40/70/94)・Q_w が FP64 の B0 の到達 (系列の 125000) から \|Δθ_r\| ≤ 0.05 %、\|ΔQ_w\| ≤ 0.1 % (line-implicit-speed §6.18 の R3 と同じ閾値)。
  - 記録: step 数、`check_convergence`、到達窓の `check_quasisteady`、末尾 2 万 step の全残差の中央値、|dq|/ULP の統計、壁の熱流束の分布の FP64 との差。
- **V5 速さ**: 専有・交互に 1000 step × 3 本。新しい float が FP64 の 0.70 倍以下。
- **V6 標準ケースの回帰**: float で動いている標準ケース (例: case/08 の bump の node + 陰解法)。
  - 新しい float と今の float の両方を FP64 (同じコミットの double のビルド) と比べ、新しい float が今の float より FP64 に近いか同等であること。
  - 物差しと対象の量は上位の判断で決める。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- すべてのケースのカーネル (粘性・拡散・再構成・リミタ・陰解法の対角・壁関数)。float のビルドの回帰の値が変わる (案 A の場合)。
- メモリ: 面ごとに float 5〜11 個。

## 変更ログ

- 2026-10-10: 起票 (draft)。ユーザ「float化よろ」。調査と §4・§6 の案。
