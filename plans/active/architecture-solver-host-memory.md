# forge 本体のホストメモリ削減 — GPU 経路で使わないホスト側の写しと行列を確保しない

## メタ

- **area**: `architecture`
- **status**: `draft`
- **related_docs**:
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) §6.4 (メモリの置き場所と計測)
- **related_plans**:
  - 起点: [`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B4b-1・B4b-3 (SERN 全体の全ヘキサ化で、g5 の律速が本体のホスト RAM だった)
  - 同型の先例: 同 plan §5.1 B4-5 (変換器の省メモリ化。3716 → 1579 B/節点、出力はデータセット単位で不変)
- **created**: `2026-10-07`
- **owner**: `CFD Dev (Claude セッション: SERN 3D)`

## 1. 目的

GPU 経路 (`gpu: 1`) の forge 本体は、計算に使わないホスト側の配列を初期化の後まで持ち続けており、
ホストのピークが約 2.7 kB/節点になっている。AWS g5 (ホスト 16 GB、他セッションと共有) では、これが格子の上限 (約 500 万節点) を決めている。
計算結果を変えずにホストのピークを下げ、同じ環境で回せる格子を大きくする。

## 2. スコープ

- **やる (段階 1)**: R1 使っていない行列 `mat_ns` を作らない、R2 GPU 経路でホストの面変数 `p` を確保しない、
  R3 GPU 経路でホストのセル変数 `c` をホストで読み書きする変数だけに絞る。いずれも `gpu: 1` のときだけ。
- **やらない (段階 2 以降、段階 1 の結果と SERN の規模判断の後に別途決める)**: R4 初期化の後の `msh.planes`/`cells` の解放
  (常駐は下がるがピークは `readMesh` に移るので、R5 と組みでないと効かない)、R5 `readMesh` の平坦な読込 (h5 から vector-of-structs を作らない)、
  R6 GPU 側で使っていない変数の確保をやめる (約 50 本、約 220 B/節点の候補。カーネルが未確保のポインタを読まないかの確認が要る)。
- **やらない**: CPU 経路 (`gpu: 0`) の挙動変更、数値スキーム・カーネルの変更、出力形式の変更。

## 3. 関連 docs と前提

- 計測の原本: [`notes/investigations/2026-10-07-forge-memlog/`](../../notes/investigations/2026-10-07-forge-memlog/) (`memlog_table.txt`、
  縮小した生産 SERN 3D 格子 25.4 / 67.3 / 118.9 万節点、生産と同じ本段設定、10 step、`FORGE_MEMLOG=1`)。
  AWS の g3 192 万節点の実測 (`case/46.sern_design/run_1071_b4b1_mem`、ホスト VmHWM 5.04 GiB、GPU 2974 MiB) と外挿が一致。
- 計測で分かったこと (2026-10-07):
  - ホストのピークは一時的な山でなく**常駐配列**で決まる (HWM と終了時 VmRSS の差は 3 % 以内)。最終 HWM は `setStructuralVariables` の時点 (常駐 + 約 150 MB)。
    一時的な山は `readMesh` の平坦な読込配列 (約 246 B/節点) と `setStructuralVariables_d` の作業配列 (+118 MB) だけ。
  - 常駐の内訳 (B/節点): ホストの `c` 964 (223 本)、`msh.planes` 705、`mat_ns` 312、`msh.cells` 254、ホストの `p` 195 (16 本)、`msh.nodes` 104、`vizCONNE` 35、出力用の局所格子 36。
    傾き 2626 B/節点 + 約 200 MB。GPU は 1395 B/節点 + 約 300 MB (context 別)。
  - 計算時間は律速でない (g3 で 38.8 ms/step、約 20 ns/節点・step)。
- 計測用のコード (`mesh/memlog.hpp`、`main.cpp` の `MEMLOG_SOLVER`、`mesh/mesh.cpp`) は既定で無出力・計算に触れない。

## 4. 設計方針

**原則**: 削るのは「GPU 経路でホストが読み書きしない配列」だけ。数値の経路 (デバイス側の配列・カーネル) には触れない。
**確保しなかった配列にホストから触ったら、黙って空の配列を使うのでなく、変数名つきのエラーで止める** (誤って 0 長の配列を H2D/D2H すると
GPU 側を壊すか、出力が黙って欠けるため)。

### 4.1 R1 `mat_ns` を作らない (見込み −312 B/節点)

`main.cpp` の `mat_ns.initMatrix(msh)` (起動ログ "Init Matrix (but not used now)") を `gpu: 1` では呼ばない。
2026-10-07 の確認: ビルド対象 (`CMakeLists.txt` の `forge_cppfiles`) で `mat_ns` のメンバ (`structure`・`lhs`・`rhs`・`localPlnOfCell`) を読む箇所は無く
(読むのはビルド対象外の `solvePoisson_amgx.cpp`・`solveNavierStokes.cpp` だけ)、各 wrapper は参照を受け渡すだけ。
`matrix` オブジェクト自体は残す (関数の引数として渡っているため)。

### 4.2 R2 ホストの面変数 `p` を確保しない (見込み −195 B/節点)

`variables::allocVariables` (`variables.cpp`) で `useGPU == 1` のときホスト側の `p[name]` を resize しない (デバイス側 `p_d` は今までどおり)。
確認: ホスト側の `p` を読むのは CPU 経路 (`variables::setStructuralVariables` の `gpu==0` 部、`setStructualVariables.cpp`、`gradient.cpp` の `gpu==0` 部)
とビルド対象外のファイルだけで、`copyVariables_plane_H2D/D2H` の呼び出し元は無い。ガード: `copyVariables_plane_*` は 0 長のホスト配列を受けたら停止。

### 4.3 R3 ホストのセル変数 `c` を必要な名前だけにする (見込み約 −850 B/節点)

`useGPU == 1` のとき、ホスト側の `c[name]` を `nCells_all` 長で確保するのは **ホスト確保集合 H** に入る名前だけにする。
H は起動時に設定から決める (固定リストにしない):
- 出力で D2H する名前 (`output.level`・`extraFields`・境界出力・`res_wall` 等で実際に書く名前。level 2 は名前が増える)
- 初期場の読込 (`readValueHDF5`) で読む名前、restart で読む名前
- 初期化でホストが書いて H2D する名前 (入口分布・壁距離・共役伝熱・凝縮/種の初期化など。**監査で列挙する**)
- probe (`pprobes`・`probe.yaml`) が読む名前、環境変数で有効になる診断 (`FORGE_IMPLICIT_DIAG_CSV` 等) が読む名前
- 実行中にホストで読む名前 (残差・力・積分量の集計で D2H するもの。監査で列挙する)

**H に無い名前のホスト配列は長さ 0 のまま**にし、`copyVariables_cell_H2D/D2H` と `readValueHDF5` は長さが `nCells_all` でない配列を受けたら
変数名つきで停止する (遅延確保はしない — どこで要るかを監査で確定させ、漏れは停止で見つける)。`gpu: 0` は今までどおり全確保。

監査の方法: ビルド対象の全ソースで、ホスト側 `c` への参照 (`c[`・`c.at(`・`.c.count(`・`copyVariables_cell_*` の名前リスト・`output_cellValNames` 系) を列挙し、
GPU 経路で到達するものを H に入れる。列挙は試験 (§6) の構成行列で裏付ける (停止しなければ漏れなし、ではなく、構成ごとに出力・残差が一致すること)。

### 4.4 期待値

R1–R3 で常駐は約 2626 → 約 1270 B/節点、ピークは `readMesh` の末尾 (約 1365 B/節点 + 定数) に移る (計測の外挿)。
g5 (ホスト 16 GB、共有) で回せる格子は約 500 万 → **約 1000 万節点**の見込み。GPU (約 1.4 kB/節点、24 GB で約 1600 万) はまだ律速でない。
それ以上 (SERN 全体の全ヘキサ見積もり 1400〜1700 万節点) には段階 2 (R4+R5、R6) が要る。

## 5. 実装ステップ

1. R1: `main.cpp` の `initMatrix` を `gpu: 0` のときだけに。
2. R2: `variables.cpp` の面変数ループで `useGPU == 1` のときホスト resize を省く。`copyVariables_plane_*` に長さガード。
3. R3: ホスト確保集合 H を作る関数 (`variables` か `main.cpp` の初期化) と、`allocVariables` での条件付き確保、`copyVariables_cell_*`・`readValueHDF5` の長さガード。監査結果を本 plan に表で残す。
4. 計測と回帰 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan レビュー | 本 plan の §4・§6 を `codex_review.py --stage plan` に。エスカレーション条件 1 の諮問を兼ねる (設計の選択肢は計測で 1 つに絞れているため) | F (採否) |
| 2 | R3 の監査 | ホスト側 `c` の参照を列挙し、GPU 経路で到達する名前と理由の表を §4.3 に追記。触るファイル: なし (読むだけ) | O |
| 3 | R1–R3 の実装 | §5 の 1〜3。触るファイル: `main.cpp`、`variables.cpp`/`.hpp`、必要なら `output/output.cpp`。合格: §6 の全項目 | O |
| 4 | 計測と回帰 | §6 | O |

## 6. 検証

事前に決める合格条件 (結果を見てから作らない):

- **メモリ**: 計測と同じ縮小格子 3 点 (25.4 / 67.3 / 118.9 万節点、`notes/investigations/2026-10-07-forge-memlog/` と同じ入力) を `FORGE_MEMLOG=1` で回し、
  ホスト VmHWM の傾きが **≤ 1500 B/節点** (計測前 2723)、GPU の傾きが計測前 (1395 B/節点) から **±2 % 以内** (デバイス側は変えていないことの確認)。
- **計算結果が変わらない**: 次の各構成で、変更前後のビルドを同じ入力で回し、
  (a) step 0 の残差行と初期出力 (`res_0` 相当) の全データセットが一致 (残差は既知の 1 ulp の非決定性 [atomicAdd、計測前ビルド同士でも出る] を除きビット一致)、
  (b) 数十〜数百 step 後の残差・出力が、変更前ビルド同士の再実行のばらつきの範囲内 (両側 3 回以上で幅を取る、[[noise-floor-both-sides]])。
  構成: 3D node SERN (縮小格子、生産設定) / 2D node (case/08.bump) / cell モード 2D / 軸対称 2D node (case/45) /
  共役伝熱 (CHT、case/63 か 64 の小さい入力) / 周期境界 (case/09 TGV 小格子) / 凝縮か多成分 (case/44 か 16 の小さい入力) /
  `output.level: 2` と `extraFields` / probe あり / restart (前 run の `res_*.h5` から)。
- **ガード**: H から 1 つ名前を外したビルドで、該当構成が変数名つきで停止すること (負例。確認後に戻す)。
- **ビルド**: `forge`・`convertGmshToForge` とも成功、構造体レイアウトを変えたらクリーンビルド ([[stale-build-struct-layout-trap]])。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

`gpu: 1` の全ケースのホスト側確保 (出力・restart・probe・診断・初期化のホスト計算)。数値カーネルとデバイス側の確保は変えない。
`gpu: 0` (CPU 経路) は変えない。

## 8. 完了条件

§6 の全項目が合格し、codex result レビューを経ること。g5 で回せる格子の上限の新しい見積もりを起点 plan (tooling-sern-mesh-blocking §5.1 B4b-2) に返す。

## 未確定事項

- 段階 2 (R4+R5、R6) に進むか — 段階 1 の結果と、SERN 全体の規模のユーザ判断を見て決める。

## 9. 変更ログ

- `2026-10-07` — 起票。ユーザ決定「forge 本体のメモリを小さくする余地を、計測から調べる」(tooling-sern-mesh-blocking §5.1 B4b-3) の計測結果を受け、段階 1 (R1–R3) を設計。
