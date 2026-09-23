# 影アキュムレータの横展開と、残差側の丸めの帰属

## メタ

- **area**: `time_integration`
- **status**: `draft`
- **related_docs**:
  - [`methods/time_integration/implementation.md`](../../methods/time_integration/implementation.md) の「commit の丸め — 定常解の到達限界を決める」
- **related_plans**:
  - [`plans/active/time_integration-fp64-accumulator.md`](time_integration-fp64-accumulator.md) (v1a。本計画の発注元)
  - [`plans/active/case-hypersonic-gap-heating-validation.md`](case-hypersonic-gap-heating-validation.md) (発見の経緯)
- **created**: `2026-09-24`
- **owner**: gap-heating セッション

## 1. 目的

v1a ([`time_integration-fp64-accumulator.md`](time_integration-fp64-accumulator.md)) は
**定常陰解法 (`tI=11`, `unsteady=0`)・node・軸対称なし・周期なし・`sstEnergyIncludesK=0`** に
範囲を絞って閉じた。本計画はそこで**意図的に外した 3 つ**を引き取る:

1. **拒否したままの経路**を外す (軸対称・SST 交換エネルギー・周期・陽解法・dual-time)
2. **残差側の丸めの帰属**を決める (v1a で運動量だけ未解決のまま残った)
3. **深部の量の定常化** (v1a の合否から外したもの。ユーザ決定、v1a §6.0a)

## 2. スコープ

- **やる**: 上の 1〜3。および v1a で保留にした double-float の再評価条件の判定。
- **やらない**: 面流束**評価**そのものの精度向上 (相殺の構造を変える話で、SLAU の定式化に踏み込む)。
  本計画は「総和を FP64 にすれば済むのか」までを決める。

## 3. 関連 docs と前提

- v1a の §4.3 (影アキュムレータの設計)・§4.4 (拒否リストと理由)・§6.1 (正本ゲート表)。
- v1a の実測: §5.2 (4 残差の `E_order`/`E_eval`)・§6.1 の G2/G3 実測表。
- **前提**: v1a の実装 (`cuda_forge/qAccumulator_d.cuh`、両 commit の FP64 分岐、
  `updateVariablesOuter_d` への reconcile 融合、起動時の拒否リスト) が入っていること。

## 4. 設計方針

### 4.1 拒否したままの経路 (v1a の S1b-①〜④)

**いずれも「原理的に無理」ではなく「Qacc をまだ扱っていない」だけ**である (v1a §4.4)。

| 経路 | v1a で拒否した理由 | 要る作業 |
| --- | --- | --- |
| node 軸対称 | `axisymmetricSource_d.cu:312-323` が commit の**基準** `roeN`/`roUyN` を射影し、呼出が commit より前 (`main.cpp:1379`) なので `Q` 側の射影は上書きされて消える | `enforceAxisSymmetry_d` に `Qacc` 5 本を渡し、軸セルで `Qacc_roUy=0` と `Qacc_roe -= 0.5·Qacc_roUy²/Qacc_ro` を当てる。`applyBlockImplicitCorrection_d:253` の axis 分岐は現状 `nullptr` で inert なので、有効化するなら同時に対応 |
| `sstEnergyIncludesK=1` | `ransTransport_d.cu:175` の `roe -= (roK-roKref)` は**増分**。陽解法は毎段 `roK−roKN` (`main.cpp:1756`, `fromN=1`) なので毎段引くと重複 | 陰解法は増分を 1 回 `Qacc_roe -= (double)Δ`、陽解法は**最終段のみ**。ON 経路では補正後に必ず `roe=(flow_float)Qacc_roe` を生成 |
| node 周期 | `periodicNode_d.cu:119` は FP32 値だけを root→member に配るので、member の `Qacc` の下位ビットが root と食い違う | `periodicMirrorQacc` を FP32 ミラー (`main.cpp:1221`, `:1753`) と**一組**にする。reconcile に root/member 不一致カウンタ (期待 0) を入れる |
| 陽解法 `tI=1`/`tI=4` | 未配線なだけ | `tI=1` は `Q_N + res·dt/v` で commit と同型。`tI=4` は**最終段 (`loop==3`) のみ** `Q_N + res_m` が commit で、中間段は捨てる暫定値 |
| 陽解法 `tI=3` / dual-time | commit の基準が**別配列の FP32 値** (`Q_N`/`Q_M`、`Q_N`/`Q_NN`) | 履歴配列と並べて **FP64 影を 1 段ずつ shift** する (機構は同一・コピーのみ)。費用は履歴 1 段につき +40 B/CV |

**実使用の分布** (全 3172 run、v1a §4.4): `tI=11 unsteady=0` 2131 / `tI=11 dual-time` 456 /
**`tI=3` 516** / `tI=4` 49 / `tI=1` 14。**`tI=3` と dual-time で 970 run** なので、
「陽解法に対応する」と言うなら `tI=3` を含めないと誤解を生む。

### 4.2 残差側の丸めの帰属 (v1a で未解決)

v1a §5.2 の実測 (深部 z/W>10、16607 CV):

| 残差 | `E_order` (再実行の変動) | `E_eval` (FP32/FP64 の評価差) | 比 |
| --- | --- | --- | --- |
| `res_ro` | 1.111e-17 | 2.915e-11 | 3.8e-07 |
| **`res_roUx`** | 6.728e-09 | 2.133e-08 | **0.315** |
| **`res_roUy`** | 5.583e-10 | 1.743e-09 | **0.320** |
| `res_roe` | 3.640e-12 | 1.747e-05 | 2.1e-07 |

**質量・エネルギーは並べ替え上界で総和丸めを除外できる** (3608/3791 CV で並べ替えが起きているのに
max|Δ| が 2.1e-16 / 5.8e-11)。**運動量は除外できない**。

⚠ **`E_order` は再実行の変動しか測らない**ので、同じ順序で毎回出る総和の丸めは現れない
(反例: `[1e8, 1, -1e8]` と `[-1e8, 1, 1e8]` は**どちらも FP32 で 0**、厳密和 1)。
**だから比 0.32 は「誤差寄与率」ではない**。

**決める方法**: **同じ FP32 面流束を固定して、FP32 総和と FP64 総和を比べる**。
面流束を配列に保存し (質量束 `massflux[ip]` は既存、残り 4 本を追加)、
セル中心に固定順で gather する経路を診断用に作る。

### 4.3 gather の選択肢と費用

| 案 | 中身 | 得るもの | 代償 |
| --- | --- | --- | --- |
| FP64 `atomicAdd` | 残差配列を double に (sm_60 以降ハード対応) | 総和誤差が相対 1e-16 級へ | 残差配列 2 倍。**決定論にはならない** |
| 面流束保存 + 固定順 gather | `cell.iPlanes`/`iPlanesDir` (`mesh.hpp:41-42`) と device CSR (`mesh.cpp:872`) が既にある | **ビット再現性**。run 間ノイズが消える | 対流 5 変数で概ね 16 B/面 + 1 パス。面並列の連続読みが散在読みに変わる |
| 圧力の基準差分 | 面の圧力力から局所基準を引く。現状 `space.pRef` は**既定 0 で実効的に何も引いていない** | 総和以前の相殺が減る | 面の反対符号共有・閉性・軸対称ソースとの整合を壊しうる |

⚠ **対流だけ決定論化しても run 間ノイズは消えない**: `viscousFlux_d.cu:379` と
`periodicNode_d.cu:27` にも残差への `atomicAdd` がある。

### 4.4 double-float の再評価条件

v1a §4.3a で**見送った**。理由はメモリ半分が成立しないこと (検出のために上位を別に持つと
FP64 影と同じ 8 B/変数、`hi==Q` にすると `dependentVariables` の commit 前書き込みで黙って壊れる)。
**再評価する条件は「速度の実測で FP64 加算が効いていると分かったとき」だけ**。
v1a の G5 は `case/56` (65k CV) で差が区別できなかったので、**代表的 3D 規模で測り直したとき**が判断点。
試作 (`cuda_forge/qAccumulatorDF_d.cuh`) と試験 (`tests/unit/test_qacc_df.cu`) は残してある。

## 5. 実装ステップ

### 5.1 残作業 (優先順)

| # | 担当 | 項目 | 内容 |
| --- | --- | --- | --- |
| 1 | F | **残差側の帰属を決める A/B** | §4.2 のとおり面流束を固定して FP32/FP64 総和を比較。**対象は運動量**だが**他の保存量も測る** (質量・エネルギーを除外する理由は並べ替え上界であって比ではない)。<br>**事前登録**: FP64 総和で `res_roU*` の `E_order` 級の差が消えれば **A** (gather の精度が効く) / 消えなければ **B** (面流束値の評価が効く)。どちらでも**実用影響が無いこと**は v1a §6.1 G2 で確認済みなので、結果は記録して優先度を決めるだけ |
| 2 | O | **軸対称** (v1a S1b-①) | §4.1。合格条件: v1a の反例 (ρ=1, ρu_r=0.125, ρe=2, dq=0 → `Q` と `Qacc` の**両方**で ρu_r=0, ρe=1.9921875) を単体試験に追加して PASS。`case/44` で採用カウンタ 0、ON/OFF の場差が事前登録の許容内 |
| 3 | O | **SST 交換エネルギー** (S1b-②) | §4.1。合格条件: RK 反例 (roKN=8、各段 9/10/11/12 → 正しい補正 −4、全段累積の −10 でないこと) の単体試験 |
| 4 | O | **周期** (S1b-③) | §4.1。合格条件: 陰解法・陽解法とも root/member 不一致カウンタ 0。`case/09` で確認 |
| 5 | O | **陽解法 `tI=1`/`tI=4`** (S1b-④) | §4.1。`unsteady=1` も許可 (陽解法に BDF 履歴は無い)。合格条件: `case/09 run_0046` 設定で ON/OFF の場距離 ≤ OFF の run-to-run ノイズ床 (**両側測定**)、採用カウンタ 0 |
| 6 | F | **`tI=3` と dual-time** | §4.1 の「履歴配列と並べて FP64 影を shift」。**970 run が該当**するので優先度は低くない。費用は履歴 1 段につき +40 B/CV |
| 7 | F | **深部の量の定常化** | v1a が合否から外したもの (§6.0a)。全域 FP64 でも `DRIFTING` なので、**まず「到達可能な量なのか」を決める**。到達不能なら指標を変える |
| 8 | O | ~~`mdot_decay.py` の量の混在~~ **済 (2026-09-24)** | 保存量 `roUy` を直接積分する形に修正。`run_0026/res_400000` で 5.903983e-10 → 5.880652e-10 (**0.40 %**)、減衰率の結論は不変。`signed=True` も足した (床付近は符号が零をまたぐので絶対値だけ見ると「張り付いている」ように見える) |
| 9 | O | 3D 規模での G5/G6 | v1a は `case/56` (65k CV) のみ。**代表的 3D 規模**で速度とメモリを測る。これが §4.4 の double-float 再評価の判断点になる |
| 10 | F | codex レビュー (plan 段) | §4 の設計方針と §6 の検証計画が書けた時点、実装着手前 |

## 6. 検証

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 対応 |
| --- | --- | --- | --- | --- |

(まだ無い。`status: draft` のうちは不要。実装着手前に `plan` 段を回す → §5.1 #10)

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/` の `axisymmetricSource_d.cu` / `ransTransport_d.cu` /
  `periodicNode_d.cu` / `timeIntegration_d.cu` / `update_d.cu`
- `solver_density_cuda/main.cpp` の起動時拒否リスト (経路を外すたびに縮む)
- `case/56.gap_tp1187/tools/mdot_decay.py` (#8)

## 8. 完了条件

- [ ] §5.1 の #1〜#9 が完了し、それぞれの合格条件を満たしている
- [ ] 起動時の拒否リストが `tI=3` と dual-time 以外を通す (または #6 も完了して空になる)
- [ ] `methods/time_integration/implementation.md` の現在仕様が対応範囲を反映している
- [ ] codex レビューを plan 段と result 段で受け、指摘の採否を §6.1 に記録している

## 変更ログ

- `2026-09-24`: 起票。v1a ([`time_integration-fp64-accumulator.md`](time_integration-fp64-accumulator.md)) が
  範囲を絞って閉じるにあたり、外した経路・未解決の帰属・合否から外した量を引き取る。
