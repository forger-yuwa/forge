# forge 本体のホストメモリ削減 — GPU 経路で使わないホスト側の写しと行列を確保しない

## メタ

- **area**: `architecture`
- **status**: `in_progress`
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

**初期化の順序 (2026-10-07 plan レビュー M2)**: 現状は `allocVariables` (`main.cpp:1715`) の後、初期化の終わりで環境変数による出力登録
(`FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP`、`main.cpp:3417`) が走り、出力側は `h0` を `Ht`・`k` から組み (`output.cpp:117`)、
`/CHECKPOINT` の履歴変数は通常の出力名リストと別 (`output.cpp:217` 以降) になっている。これを
**変数・診断の登録 → 出力と checkpoint が要る名前の確定 → H の構築 → 確保** の順にする。
出力が要る名前は**出力側と確保側で同じ関数**から作る (`output.cpp` を共通化の対象にする。別実装で名前を再現しない — 出力を変えたときに食い違うため)。

**ホスト参照の保護 (同 M1)**: 転送関数 (`copyVariables_cell_*`) や `readValueHDF5` の入口だけで長さを見ても、その前にホストへ書く経路
(例: dual-time の checkpoint 復元 `main.cpp:1148` は `var.c.at(nm)` に `v[i]` で書いてから H2D する) は守れない。
GPU 経路のホスト参照は、**変数名と期待長を検査する共通アクセサ** (例 `var.hostCell(name)`、長さが `nCells_all` でなければ名前つきで停止)
を通す。添字アクセス・イテレータ演算・生ポインタの受け渡しより前に検査する。`map::at()` はキーの存在しか見ないので保護にならない。

### 4.1 R1 `mat_ns` を作らない (見込み −312 B/節点)

`main.cpp` の `mat_ns.initMatrix(msh)` (起動ログ "Init Matrix (but not used now)") を `gpu: 1` では呼ばない。
2026-10-07 の確認: ビルド対象 (`CMakeLists.txt` の `forge_cppfiles`) で `mat_ns` のメンバ (`structure`・`lhs`・`rhs`・`localPlnOfCell`) を読む箇所は無く
(読むのはビルド対象外の `solvePoisson_amgx.cpp`・`solvePoisson_amgcl.cpp` だけ。`solveNavierStokes.cpp` は参照を受け取るだけ — 2026-10-07 実装時の再確認で訂正)、各 wrapper は参照を受け渡すだけ。
`matrix` オブジェクト自体は残す (関数の引数として渡っているため)。

### 4.2 R2 ホストの面変数 `p` を確保しない (見込み −195 B/節点)

`variables::allocVariables` (`variables.cpp`) で `useGPU == 1` のときホスト側の `p[name]` を resize しない (デバイス側 `p_d` は今までどおり)。
**例外 (監査 2026-10-07)**: 診断 `FORGE_DIAG_PSI_DUALEVAL` の退避 (`main.cpp:2293-2303` の `pdeSize(s.var.p, k, 0)`) はホスト `p` の長さをデバイス配列の長さとして使い、fallback が 0 — R2 で面配列が**黙って退避から外れる**ので、fallback を `msh.nPlanes` に直す (R2 に含める)。確認: ホスト側の `p` を読むのは CPU 経路 (`variables::setStructuralVariables` の `gpu==0` 部、`setStructualVariables.cpp`、`gradient.cpp` の `gpu==0` 部)
とビルド対象外のファイルだけで、`copyVariables_plane_H2D/D2H` の呼び出し元は無い。ガード: `copyVariables_plane_*` は 0 長のホスト配列を受けたら停止。

### 4.3 R3 ホストのセル変数 `c` を必要な名前だけにする (見込み約 −850 B/節点)

`useGPU == 1` のとき、ホスト側の `c[name]` を `nCells_all` 長で確保するのは **ホスト確保集合 H** に入る名前だけにする。
H は起動時に設定から決める (固定リストにしない。出力・checkpoint の分は出力側と共通の関数で作る):
- 出力で D2H する名前 (`output.level`・`extraFields`・境界出力・`res_wall` 等で実際に書く名前。level 2 は名前が増える)
- 初期場の読込 (`readValueHDF5`) で読む名前、restart で読む名前
- 初期化でホストが書いて H2D する名前 (入口分布・壁距離・共役伝熱・凝縮/種の初期化など。**監査で列挙する**)
- probe (`pprobes`・`probe.yaml`) が読む名前、環境変数で有効になる診断 (`FORGE_IMPLICIT_DIAG_CSV` 等) が読む名前
- 実行中にホストで読む名前 (残差・力・積分量の集計で D2H するもの。監査で列挙する)
- dual-time の checkpoint 履歴変数 (復元・書出し)、`h0` を組むための `Ht`・`k` のような出力の依存名

**H に無い名前のホスト配列は長さ 0 のまま**にし、`copyVariables_cell_H2D/D2H` と `readValueHDF5` は長さが `nCells_all` でない配列を受けたら
変数名つきで停止する (遅延確保はしない — どこで要るかを監査で確定させ、漏れは停止で見つける)。`gpu: 0` は今までどおり全確保。

**監査の結果 (2026-10-07、全文 [`notes/investigations/2026-10-07-host-memory-audit.md`](../../notes/investigations/2026-10-07-host-memory-audit.md)、HEAD 9c9f623c の行番号)**:
GPU 経路でホスト `c` を読み書きするのは 5 か所だけ — (1) 初期場の読込 `readValueHDF5` (保存量・`wall_dist`・k/ω・種・遷移・トレーサ・凝縮モーメントを書いて H2D、`variables.cpp:787-939`)、
(2) dual-time の checkpoint 復元 (`main.cpp:1148`、書いて H2D) と書出し (`output.cpp:227`)、(3) 出力 (`output.cpp`: level 0/1/2・`extraFields`・`h0` の依存 `Ht`/`k`・`/CHECKPOINT`、NaN ダンプも同じ関数)、
(4) 環境変数で有効な診断 2 つ (`FORGE_IMPLICIT_DIAG_CSV` の `main.cpp:847-856`、`FORGE_PIN_DIAG` の `main.cpp:2210-2238`)、(5) `lineImplicit` の `ccx/ccy/ccz` (`main.cpp:3503`)。
probe・壁出力・CHT・残差/NaN 集計・入口分布・壁距離・種/凝縮/受動種の初期化は、デバイスか自前の局所バッファで済み**ホスト `c` を経由しない**
(§4.3 当初の「probe が読む名前を H に入れる」は不要)。GPU 経路にホスト `c` の `operator[]` (暗黙生成) は無い。名前の集合を変える環境変数は 10 個
(登録・除去: `FORGE_WI_FORCE_DIAG`・`FORGE_WF_OMEGA_SOURCE`・`FORGE_WF_CLOSURE_DIAG`・`FORGE_WF_REP_DIAG`・`FORGE_OMEGA_BUDGET`・`FORGE_SPECIES_RAW_DIAG`、出力名追加: `FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP`、ホストで読む: `FORGE_IMPLICIT_DIAG_CSV`・`FORGE_PIN_DIAG`)。
**共通関数の設計 (監査 §3 の案を採用)**: `output/outputFieldNames.hpp` (新規) に `outputFieldPlan` (出力名の列・`h0` 依存・checkpoint 履歴)、`dualTimeHistoryNames`、`initialValueNames`、
`hostCellSet` (H = それらの和 ∪ 診断 ∪ lineImplicit)、`registerOutputDiagnostics` (`main.cpp:3417-3431` を確保の前へ移す)、`applyEnvGatedRemovals` (`variables.cpp:301-356` を切り出す)。
アクセサ `variables::hostCell(name)` (キー無し・長さ違いで名前と長さを出して停止) を、`copyVariables_cell_*`・`output.cpp:207/230/263`・`variables.cpp:787-931`・`main.cpp:850-856/1148/2211-2238/3503` に通す。
`output.cpp:196` の `for (auto& v : var.c)` は h5 の書込み順を変えないため残し、名前で絞った後に `hostCell` を通す。
**守る順序依存**: H の配列は `nCells_all` 長で 0 初期化を保つ (H2D はゴーストを含めて写すので、縮めるとビット一致しない)。`FORGE_RESID_SNAP` の登録条件 (`e != nullptr`) は今のまま写す。
`extraFields` の受付は「登録済み (cellValNames)」で判定 (`c_d` は実行中にキーが増えうる)。`variables.hpp` にメンバを足すとレイアウトが変わる — クリーンビルド。
`sstF1 = 1` の H2D がループ内で繰り返される件 (`variables.cpp:381-384`、最終状態は同じ) は**ついでに直さない** (別件)。
監査表は **名前・有効になる条件 (設定・環境変数)・最初の利用箇所 (ファイル:行)** を列にする (全文の §1)。監査の方法: ビルド対象の全ソースで、ホスト側 `c` への参照 (`c[`・`c.at(`・`.c.count(`・`copyVariables_cell_*` の名前リスト・`output_cellValNames` 系) を列挙し、
GPU 経路で到達するものを H に入れる。列挙は試験 (§6) の構成行列で裏付ける (停止しなければ漏れなし、ではなく、構成ごとに出力・残差が一致すること)。

### 4.4 期待値

**ピークは工程別の「傾き + 切片」の最大値で見積もる** (同 m1)。`setStructuralVariables` の作業配列は定数でなく
約 148 B/節点 + 3.6 MiB (`variables.cpp:567` が面数・セル数に比例して確保、`memlog_table.txt` の HWM − RSS が 40/98/172 MiB)。
その時点の HWM の傾き 2751 B/節点から R1–R3 (312 + 195 + 850) を引くと **約 1394 B/節点**で、`readMesh` 末尾 (1365) より大きい —
ピークが `readMesh` に移るとは言えず、`setStructuralVariables` 付近に残る見込み。g5 (ホスト 16 GB、他セッションと共有) で回せる格子は
約 500 万 → **約 1000 万節点 (暫定)**。正常終了した run の実測と、共有環境の実効予算から更新する。GPU (約 1.4 kB/節点、24 GB で約 1600 万) はまだ律速でない。
R1 → R2 → R3 の各段階で効果を測る。
それ以上 (SERN 全体の全ヘキサ見積もり 1400〜1700 万節点) には段階 2 (R4+R5、R6) が要る。

## 5. 実装ステップ

1. 監査 (§5.1 #2) と基準入力の確定 (§5.1 #3) — 実装の前に。
2. R1: `main.cpp` の `initMatrix` を `gpu: 0` のときだけに。測る。
3. R2: `variables.cpp` の面変数ループで `useGPU == 1` のときホスト resize を省く。`copyVariables_plane_*` に長さガード。測る。
4. R3: (a) 出力・checkpoint が要る名前を作る共通関数 (`output/output.cpp` と確保側で共用)、(b) 初期化の順序の入れ替え (環境変数の出力登録を確保の前へ)、
   (c) H の構築と `allocVariables` での条件付き確保、(d) GPU 経路のホスト参照を共通アクセサへ。触るファイル: `main.cpp`、`variables.cpp`/`.hpp`、`output/output.cpp` ほか監査で挙がった所。測る。
5. 回帰 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| ~~1~~ | ~~codex plan レビュー~~ (**済 2026-10-07**: GO-with-changes M5/m1、全件採用 — §6.1。エスカレーション条件 1 の諮問を兼ねた [設計の選択肢は計測で 1 つに絞れているため]。判断: 2026-10-07・R1–R3 の方針は維持、監査と受入条件を先に確定) | | F |
| ~~2~~ (**済 2026-10-07**: [監査全文](../../notes/investigations/2026-10-07-host-memory-audit.md)、要点は §4.3。R2 に `pdeSize` の fallback 修正を追加) | **監査** (M1・M2・M5) | ホスト側 `c`/`p` の全参照を、**名前・有効条件・最初の利用箇所 (ファイル:行)・書くか読むか・保護の要否**の表に (§4.3)。出力と checkpoint の依存名 (`h0` ← `Ht`/`k`、`/CHECKPOINT` 履歴、`FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP` 等の環境変数) を含める。表の各行に対応する**試験 (入力・確認する成果物・判定)** を割り当てる (§6 の構成表の元)。触るファイル: なし (読むだけ)。合格: 表が本 plan に入り、§6 の構成表と対応している | O |
| 3 | **基準入力の確定と回帰ハーネス** (M3・M5。**2026-10-07 途中**: `case/66.hostmem_regression/` に入力 23・構成 30 [forge 25・変換器 5] とハーネス `matrix_spec.py`・`prepare_inputs.py`・`run_matrix.py`・`compare_runs.py`・`memlog_summary.py`。base r1 は 28 構成が完走・NaN 0、SERN g3 base 3 回 [`case/46` run_1072–1074] 完走。**base r2/r3 の途中で AWS が idle 自動停止** — base バイナリを `forge_9c9f623c` の名前で起動していたため `idle_autostop.sh` の `pgrep -x forge` に掛からず、短い 2D run の合間を idle と数えられた。ハーネスは `.bin/<ビルド>/forge` の symlink 経由に修正済み (未実行)。**g4 は最終場が AWS に残っておらず不可** → メモリの傾きはローカル s070/s085 の 2 点 (変更前後とも同じ入力) で取り、g3 は外挿の確認点とする [§6 の「g3・g4」からの変更、理由: 入力の消失]。保留: NaN ダンプ (入力が別セッションの scratch)、case/13 の cell) | **専用 case `case/66.hostmem_regression/`** に、§6 の構成表の各構成の入力を (元 case から**複製して**。元の case・他セッションの作業ディレクトリには書かない) 置き、変更前バイナリ (AWS `~/bin-hostmem/forge_9c9f623c`) で各 3 回・N step を回すハーネス (投入と比較のスクリプト) を作る。各構成が現行バイナリで起動し NaN 無く N step 走ることを確認 (古い config の廃止キーは複製側だけ直し、直した内容を README に書く)。SERN g3/g4 は `case/46.sern_design` に run_1072 以降で。2026-10-07 の縮小格子 3 点は**発散前までのメモリ観測**として保持し、回帰の基準には使わない (`s050` は step 9 で `roe` 非有限、3 点とも `check_convergence` は NOT CONVERGED)。変更後も同じハーネスで回す | O |
| 4 | R1–R3 の実装 (**R1・R2 は済 2026-10-07**: gpu: 1 で `mat_ns` とホスト `p` を確保しない、`pdeSize` の面 fallback を `nPlanes` に、面配列の長さガード `requireHostPlaneLength`。ローカル縮小格子 s070/s085 でホスト HWM の傾き 2737 → 2227 B/節点 [−510、期待 −507]、GPU は不変、s070 の step 0 は既知の 1 ulp 列を除き一致・`res_0.h5` は全 25 データセット一致、`FORGE_DIAG_PSI_DUALEVAL` の退避は 543/543 で base と同じ。`variables.hpp` にメンバ `nPlanesAlloc` を追加 = **構造体レイアウトが変わったので他のビルドはクリーンビルド**。本格回帰は #5 の AWS ハーネスで。R3 も済 2026-10-07 (ローカル): 共通関数 `output/outputFieldNames.{hpp,cpp}` [出力名・h0 依存・checkpoint・初期場・診断・H]、順序「登録 → `registerOutputDiagnostics` → `applyEnvGatedRemovals` → `hostCellSet` → `allocVariables(…, H)`」、アクセサ `hostCell` [キー無し・長さ違いで名前つき停止] を転送・出力・初期場読込・checkpoint 復元・診断・lineImplicit に。s070/s085 でホスト HWM の傾き 2737 → **1371 B/節点** [R3 単独 −856、合計 −1366]、ピークは `setStructuralVariables` (1371) と `readMesh` 末尾 (1362) がほぼ並ぶ、GPU 不変。既定構成 (level 1) の H は 223 本中 23 本、**level 2 では大きい** (case/09 で 241 本中 98 本)。s070 の step 0・`res_0.h5` (25 本) 一致、case/09 dual-time の `res_0` (level 2、100 本) 一致・checkpoint 14 本と属性一致・再開直後の `/CHECKPOINT` が入力とビット一致、env 変種の `res_0` (106 本) 一致。負例 2 つ (`P`・`roN` を H から外す) は書込み/転送の前に名前つきで停止。変換器は旧 (7637528e) と新で s070 の変換結果 398 データセットが一致。**本格回帰 (#5、AWS ハーネス) が残り**) | §5 の 2〜4。各段階でメモリを測る (g3/g4、`FORGE_MEMLOG=1`) | O |
| 5 | 回帰 (**2026-10-07 実施: 登録判定は (a) と (b) 4 構成 7 量が FAIL、他は PASS — §6.2**。追加診断 B は全量 PASS だが 7 量を区別できない・new 側の異常も幅を広げる → §6.3 の固定幅の独立 A/B は**判定不能** (新規 base も超過)。判断: 2026-10-07 codex diagnose・限定事項つき `done` も却下。次案 (経路の一致) は方向は採用だが CUPTI と読戻し再供給の基盤が要り、それでも `done` 不可 — **検証の深さをユーザ判断へ (2026-10-07)**) | §6 の全項目。負例 2 つ (出力用変数を H から外す・checkpoint 履歴変数を外す → いずれも書込み/転送の前に名前つきで停止) | O |

## 6. 検証

事前に決める合格条件 (結果を見てから作らない。2026-10-07 plan レビュー M3–M5・m1 を反映):

- **メモリ**: 正常な基準入力 (§5.1 #3 の g3・g4) を `FORGE_MEMLOG=1` で変更前後に回し、工程別の傾き + 切片を出す。
  ホスト VmHWM の傾き **≤ 1500 B/節点** (変更前 2723)、GPU の傾きは変更前から **±2 % 以内** (デバイス側を変えていないことの確認)。
- **比較の時点を変数ごとに固定する (M4)**: 初期出力 (`writeInitialOutputs`、残差の組立より前) では**初期化済みの保存量・原始量・幾何量だけ**を比べる。
  残差 (`res_*`)・補正量 (`dq_*`) などは最初に定義される組立・更新の後の出力で比べる。level 2 の初期出力には確保時に 0 初期化されない平均流残差
  (`variables.cpp:365` の対象外) が入るので、そのビット列は合格条件にしない (未定義値を初期出力に書く既存の問題として §未確定事項に記録)。
  出力互換性は**データセットの集合・shape・dtype・属性**の一致も含む。
- **計算結果が変わらない**:
  (a) step 0 の残差行は全列ビット一致 (既知の 1 ulp の非決定性 [atomicAdd、変更前ビルド同士でも同じ 2 値が出る] がある列は、その 2 値のどちらかであること)。
  (b) N step 後の残差・出力の差は、**変更前ビルド 3 回・変更後ビルド 3 回の同ビルド内ペア差 (計 6 対) の最大の 2 倍以内** (データセットごとに
  max|差| / max|基準|)。[[repeat-range-as-limit-rule-trap]] (同じ分布でも「旧 3 回の範囲」では 8 割が不合格になる) を避けるため、両側をまとめ係数 2 を置く。
  (c) 全保存量・全出力に NaN/Inf が無い (各 run で確認)。
  (d) 定常ケースで「一致」と言うときは両側の `check_convergence.py` と、派生量 (推力・力係数など) の `check_quasisteady.py` の VERDICT を貼る。
  短い試験の合格は「回帰差が許容内」に限定し、収束の主張に使わない。非定常ケースは同じ物理時刻で比べる。
- **構成** (監査表の到達分岐ごとに具体化する。下は初期の割り当て):

  | 構成 | 入力 (元。複製して使う) | 確認する成果物 |
  | --- | --- | --- |
  | SERN 3D node 生産設定 + メモリ | AWS `case/46` g3 (`run_1068` 最終場から restart)、g4 (`run_1045`/`run_1046` 最終場) | 残差 CSV・`res_*.h5`・力の時系列・`FORGE_MEMLOG` |
  | 2D node 標準 (A・G・H・probe・壁出力) | `case/36` `run_sym_H_2up_node` の最終場から (probe を 2〜3 点足す) | `res_*.h5` (集合・shape・dtype・属性 `h0_includes_k` 含む)、残差、壁 h5、probe 出力 |
  | 同 + `FORGE_IMPLICIT_DIAG_CSV` / `FORGE_DIAG_PSI_DUALEVAL` | 同上 | `diag.csv`、`psi_dualeval.csv` と log の `snapshot: X of Y device arrays` (X が減ったら R2 の fallback 漏れ) |
  | dual-time の checkpoint (F・I・J・K・L・B・D・N) | `case/09` `run_0160`/`0162`/`0164` 系 (node 周期・種・トレーサ・FCT・BDF2・level 2・extraFields)。100 + 再開 100 と連続 200。env 変種: `FORGE_OUT_RESIDUALS=1`+`FORGE_RESID_SNAP=0`、`FORGE_SPECIES_RAW_DIAG=1`、`FORGE_PIN_DIAG=1` | `/CHECKPOINT` の全データセットと属性、log の履歴復元行、`res_*_m`・`dq_*_new`・`roYraw*` の有無、分割と連続の差 |
  | 軸対称・多成分・凝縮 | `case/44` の定常 active run と dual-time 版 (`run_0376` 型) | 種・モーメント、`/CHECKPOINT/<cons>P`、pin 診断 |
  | 共役伝熱 | `case/52` `run_0007_fxhalf` (node・`interfaceDiag 1`) | 壁 h5 (`iface*`)・CHT の CSV |
  | cell モード | `case/36` `run_sym_H_2up_cell`、`case/20` `001.test/run_slau`、cell dual-time `run_case04_les_unsteady_dualtime` | 出力・残差 (cell は CONNE を cells から組む) |
  | 遷移モデル (C・lm*) | `case/57` `run_0014_t3a_lm_unitcheck` (roGamma を読む経路) と、SST だけの場から LM を始める経路 | `lm*`、log の遷移初期化行 |
  | line-implicit と extraFields (O・J) | `case/56` `run_0019_lineimplicit`、`run_0027_s6_f32_b`、`case/48` `run_0903_absorb2` | 出力の集合、ライン構築の log |
  | 環境変数で登録が変わる診断 (I) | `case/26` `run_0086_optin_base` を env 無しと `FORGE_WI_FORCE_DIAG`・`WF_CLOSURE_DIAG`・`OMEGA_BUDGET`・`WF_REP_DIAG`=1 で | データセット集合、警告行 |
  | 変換器 | 上の各入力の msh (変換器も `variables.cpp` を共有する) | 変換結果のデータセット単位比較 |

- **ガード (負例 2 つ)**: (1) `case/36` で H から `P` を外したビルド → `writeInitialOutputs` の D2H より前で、(2) `case/09` の再開で H から `roN` を外したビルド → `main.cpp:1148` の書込みより前で、いずれも**変数名つきで停止**すること (確認後に戻す)。
- **ビルド**: `forge`・`convertGmshToForge` とも成功、構造体レイアウトを変えたらクリーンビルド ([[stale-build-struct-layout-trap]])。

### 6.2 回帰の結果と追加診断 (2026-10-07)

**登録判定の結果 (base 9c9f623c 対 new 93e55957、各 3 回、`case/66.hostmem_regression/`、原本 `results/2026-10-07_base9c9f623c_vs_new93e55957/`)**:
(c) NaN・初期出力 (保存量・原始量・幾何量のビット一致)・出力互換・ログ行・変換器 5 構成は全構成 PASS。**(a) step 0 は 1 構成を除き FAIL**、**(b) D ≤ 2S は 4 構成・7 量で FAIL**
(c44dual_ckpt100 `condClampCorrQ_0`・`condR30_0`、c44dual_restart100 `condClampCorrQ_0`、c20cell_rk3 `rms_roUz`、c20cell_dual `Uz`・`roUz`・`CHECKPOINT/roUzN`)。
SERN g3 のホスト VmHWM 5157 → 2630 MiB (2816 → 1436 B/節点)、GPU 不変 (メモリの合格条件は満たす)。**登録判定は FAIL のまま記録する** (書き換えない)。

**判断: 2026-10-07 codex (diagnose) [記録](../../notes/reviews/2026-10-07-hostmem-regression-fail-diagnose.md) — 全件採用**:
(1) 比較器に欠陥がある (呼び出し側で確認済み): `compare_runs.py` の `metric` は `max|A−B|/max|A|` で**分母がペアの片側**、同ビルド内は `i<j` の片方向・ビルド間は base→new の全組合せなので**S と D の尺度が揃わず、反復の並び順で判定が変わる**
(最小再現: base `[100,2,1]` new `[100,2,1]` が FAIL、new を `[1,2,100]` に並べ替えると PASS)。さらに `inf ≤ 2·inf` が真になり**比較不能を PASS にする** (base `[0,1,1]` new `[100,100,100]` が PASS)。
(2) (a) の前提 (1 ulp の 2 値) は誤りで、base 自身が 19–75 ulp (2D 定常)、52 ulp (SERN)、1e6 ulp 級 (dual-time は step 0 に内反復 22 行) に割れていた — **(a) は判別力を持たなかった**。(a) は「初回組立」と「最初の物理 step 全体」を混同している。
(3) 提案した「初期化直後の全デバイス配列のハッシュ」は却下 (一部配列は 0 初期化されず未初期化領域を含む、初期化中の境界勾配に `atomicAdd` がある [`calcGradient_d.cu:251`]、`qacc_d` 等 `c_d` 外の状態もある)。反復を増やした順位和検定も今は却下 (比較指標の修正が先、非有意は同等の証明にならない)。
(4) 「差はすべて非決定性」はまだ解釈。`roUz`・`CHECKPOINT/roUzN` は名目ゼロでも保存量で、診断量として除外できない。

**追加診断 (事前登録 2026-10-07、追加の forge 実行なし)**: 既存の base 3 本・new 3 本の同じ保存時点・同じ配列・同じ CSV 行で、尺度だけを替える。
- **A**: 登録済みの `m(A,B)` の結果 (上)。
- **B**: `d(A,B) = max|A − B|` (絶対 L∞)。同ビルド内 6 対の最大を `S_abs`、ビルド間 9 対の最大を `D_abs`。**FAIL の 7 量だけでなく従来 PASS の全量を再評価**。
- 条件: 全入力が有限で、shape・列・行キーが対応すること (比較不能は FAIL、PASS にしない)。`D_abs ≤ 2·S_abs`、`S_abs = 0` なら `D_abs = 0`。整数・構造情報は厳密一致。
  **反復の並べ替えと base/new の交換で結果が変わらない**こと (比較器の自己検査)。
- B で超過が消える量: その量の登録 FAIL は**尺度依存で説明できる**と記録 (「追加診断で反復内差の 2 倍以内」)。B でも超過する量: 比較器だけでは説明できない → 変更起因の差を候補に戻す (諮る)。
- **この追加診断の PASS は、登録判定の書換えや「非決定性だけだった」証明には使わない**。R1–R3 を合格にするかは、B の結果を見て result 段 (codex `--stage result`) と上位の判断で決める。
- やらない: run の順序を選んで PASS にする、FAIL 量を事後に除外する、S が増えるまで反復を足す、全配列を 0 初期化してハッシュを揃える。

**追加診断 B の結果 (2026-10-07、原本 `case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957_abs/`)**: 30 構成 3147 量 + SERN g3 255 量で超過 0・比較不能 0、整数・文字列 185 量は全 run で厳密一致、
並べ替え・交換の自己検査 (122472 回) で判定の変化 0、単体試験 8 件 PASS。ただし A の FAIL 7 量は B で全部 D/2S = 0.50 (S と D を同じ外れ 1 本が決める) = **B はこの 7 量で base と new を区別できていない**。
SERN g3 `res_vehicle_base_18_100:twall_z` は D/2S = 1.000 (S_abs 2⁻⁷、D_abs 2⁻⁶) の境界上 PASS。

**判断: 2026-10-07 codex (diagnose) [記録](../../notes/reviews/2026-10-07-hostmem-result-interpretation-diagnose.md) — 全件採用。R1–R3 の `done` は保留**:
(1) B の S は new 内の差を含むので、**変更後だけの異常も許容幅を広げる** (人工入力 base [0,0,0]・new [0,0,1e23] が PASS)。B の全量 PASS を受入れの根拠にしない。
(2) 7 量を「ゼロ近傍の診断量だから無害」として除外しない: `condClampCorrQ_0` は $\max_k |\Delta Q_k| / \max(|Q_{k,\mathrm{before}}|, 10^{-30})$ の最大 (`condensationRealizability_d.cuh:313`) で絶対補正量ではない (微小分母で増幅しうる)、`roUz`・履歴 `roUzN` は保存量。
(3) `twall_z` は B の境界上 PASS としては正しいが「float32 の刻みの問題」とは断定しない (最悪差の位置・6 本の値・その位置の ULP 距離を記録する)。
(4) メモリの実測は採用: g3 で 5157 → 2630 MiB、GPU 不変。1436 B/節点は**切片込みの単点値**、1371 B/節点はローカル 2 格子の**傾き**として区別して書く。
**記録の仕方**: 「A は FAIL を保存 / B は定義した条件で全量 PASS / 7 量の変更起因性は未解決 / `twall_z` は境界上 PASS・原因未確定」。
**段階 2 の注意**: 固定した許容幅・入力・バイナリ・対照 run を先に確定する。ホストピークが `readMesh` と `setStructuralVariables` で並ぶので、進めるなら R4+R5 を先に、R6 は別の変更・検証。B の方式を恒久的な受入れゲートにしない。

### 6.3 固定幅の独立 A/B (事前登録 2026-10-07、codex diagnose 同記録)

**変える因子は実行バイナリだけ**。構成 `c44dual_ckpt100` (A の FAIL 7 量のうち 2 量を持ち、外れ値の出る凝縮 dual-time) の同一入力から、base = `9c9f623c`・new = `93e55957` を**各 6 本・各 100 step**、
**順序は事前に固定して交互** (b, n, b, n, …)、新しい run ディレクトリに保存。
- **実行前に**、既存の同構成 base 3 本だけから各量の $S_0 = \max_{i<j} \|B_i - B_j\|_\infty$ を計算し、診断幅 $T = 2 S_0$ を**凍結して commit する** ($S_0 = 0$ の量は差 0 を要求)。new は幅の算定に使わない。
- 新しい各 run と既存 base 各 run の差 $\|X - B_i\|_\infty$ を同じ $T$ で評価。対象は既存の全 196 量。初期出力・構造は厳密一致、非有限は失敗。残差の行キーと物理時刻を固定。

| 結果 | 判定 |
| --- | --- |
| 新規 base・new とも全量で幅内 | 「追加標本でもこの固定幅を超える差は検出されなかった」— **限定した支持材料**として result レビューへ |
| 新規 base は幅内、new だけ超過 | 変更起因の差を優先して追う (諮る) |
| 新規 base も超過 | 基準 3 本では再現性を捉えられていない → **判定不能** (その場で $T$ や反復数を増やして合格にしない) |

これは受入れ幅の妥当性を確認する診断で、物理的に許せる誤差の上限は未定義 — 全量が幅内でも包括的な同等性は宣言しない。
やらない: S が増えるまで反復を足す、7 量を除外する、`condClampCorrQ_0` の分母床を変えて差を消す、未収束の 100–200 step 比較を定常解の一致と呼ぶ。
**凍結 (2026-10-07、段階 2 の投入前)**: `case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/T_frozen.tsv` (sha256 `45168ae5adf4f0f3896b1f85aa1aef110857bf30f4049b0818084be3b444451e`、既存 base run_0029・0038・0068 だけ)、
最終出力 196 量 (幅 121・差 0 を要求 72・厳密一致 3) + 行キー + 初期出力 44 量。投入計画 `PLAN.txt`/`plan.json` (run_0194–0205、b/n 交互、base `forge_9c9f623c` sha c65dbb39…・new `forge_93e55957` sha 7121a8e5…、入力 7 ファイルの sha を固定)。
`condClampCorrQ_0` の T = 7.25e19。**投入前の観察 (既存 run、判定に使わない)**: 約 1.46e23 の状態 (節点 19954、`roQ2_0` = `condClampCorrQ_0` × 1e-30 と 7 桁一致 = 分母が床 1e-30 だった状態と整合) は、
c44dual_ckpt100 で base 0/3・new 1/3、同じ入力に `FORGE_PIN_DIAG=1` を足した c44dual_pindiag で base 1/3 (run_0077)・new 3/3 — base にも出るので変更でしか起きない状態ではないが、new に偏る兆候。
段階 2 では**この状態に入った本数 (`condClampCorrQ_0` > 1e22) をビルドごとに記述として記録する** (判定は上の 3 行の表のまま、この本数で合否を決めない)。
あわせて (追加実行なし): `twall_z` の最悪差の位置・6 本の値・ULP 距離、`condClampCorrQ_0` が最大になった run・位置での補正前後のモーメント値 (出力にあれば) を記録する。


**段階 2 の結果 (2026-10-07、`fixedwidth_c44dual_ckpt100/RESULT.txt`、run_0194–0205)**: **表の 3 行目「新規 base も超過 → 判定不能」** (新規 base 4/6・new 4/6 が凍結幅を超過、全量幅内は base 2/6・new 2/6、構造・行キー・初期出力・整数は 12 本とも厳密一致)。
`condClampCorrQ_0` > 1e22 の状態 (`res_100.h5` の時点): 新しい標本で base 2/6・new 2/6。new だけ超えた量は 4 つ・各 1 本 (`condClampCorr_0`、`Q0_0`・`roQ0_0`・`CHECKPOINT/roQ1_0_fctH`)。
**判断: 2026-10-07 codex (diagnose) [記録](../../notes/reviews/2026-10-07-hostmem-fixedwidth-undecidable-diagnose.md) — 全件採用**: 「判定不能」で終了 (幅・本数を変えない)。**限定事項つきの `done` は却下** (物理的な許容差が未定義、§8 未達)。
「2/6 対 2/6 で頻度差は解消」も却下 (2/6 の 95 % 正確区間は約 4–78 %)。巨大値は分母だけの問題ではない (run_0106 の節点 19954 で `roQ1_0` が他の約 2.1 倍、`roQ2_0` が約 4.8 倍 = 保存量も動いている)。
ばらつきの小さい別構成の合格では代替できない (凝縮 dual-time 固有の経路)。「状態に入った本数」は `res_100.h5` 時点の値で、100 step 中に一度でも入った割合ではない (診断配列は毎 step 上書き)。
codex の次案: 同構成で base/new 各 128 本の終点事象頻度 A/B (Newcombe–Wilson 90 % 区間・±20 pt) — ただし A でも `done` にはならず、用途に基づく許容差が要る。
**現時点の記録** (codex の文案): 「ホストメモリ削減と、試験した初期化・出力・checkpoint の互換性を確認した。凝縮 dual-time の固定幅試験は変更前対照も超過して判定不能だった。step 100 の巨大診断値は両ビルドで観測されたが、出現確率・保存量への影響・物理的に許容できる差は未確定である。」
**経路の一致で決着させる案の判断: 2026-10-07 codex (diagnose) [記録](../../notes/reviews/2026-10-07-hostmem-h2d-channel-diagnose.md) — 全件採用**: 方向は採用 (`9c9f623c..93e55957` に CUDA ソースの差分は無く、`c_d`・`p_d` の確保長も維持)。
ただし提案の記録項目は不十分: 定数への転送 (`convectiveFlux_d.cu:186`)・`cudaMemset`・D2D・カーネル引数 (`dt`・係数) を含む**全操作列** (確保・解放と寿命、全方向の転送と範囲、symbol 更新、memset の値、全引数・grid/block・shared memory、stream/context と同期) が要る。
「base 同士で変わる H2D は中身を免除」は却下 (D2H → ホスト処理 → H2D の経路がある [`calcGradient_d.cu:737/799`]; base が `r`、new が `2r` を送っても免除される) — **共通の読戻し系列を両ビルドに再供給**してホスト処理を比べる。
デバイスポインタを含む転送 (`boundaryCond_d.cu:911`) は「確保 ID・世代・offset」で対応付け、読み取る領域の初期化は別に確認。記録は **CUPTI の Runtime/Driver callback** を基礎に (直接の `cudaMemcpy` の個別ラップより漏れを管理しやすい)。
これが通っても「当該入力系列でのホスト側操作生成の同等性」に限られ、通常実行の頻度差の解消でも §8 の完了でもない。**この提案だけでは `done` 不可**。一方、用途に基づく許容差は必須ではなく、デバイスコード・意味のある入力・ホスト処理・初期化/同期・出力処理の同等性を十分に立証できれば、数値演算を変えない変更として受け入れる道はある。
**以下は判断前の次案の記録**: 統計でなく**経路で決着させる** — R1–R3 が GPU の結果に影響しうる経路は (1) ホスト → デバイスの転送 (H2D) の中身と (2) デバイス側の確保 (名前・大きさ・順序) だけ (カーネル・起動設定は不変)。両ビルドに H2D の (名前・長さ・内容ハッシュ) と `cudaMalloc` の (名前・大きさ) の列を記録する診断を載せ、全構成で列が一致することを確かめる (D2H 由来で base 同士でも変わる H2D は base 同士の比較で識別し、名前・長さの一致だけを見る)。

**今後の回帰の作法 (同判断)**: 初期化済みの決定的状態のビット一致と、最初の組立後の残差比較を**別ゲート**にする。残差は `(step, inner, phase)` を固定し、
方程式別の絶対・相対許容差を、独立した base の校正と必要な検出幅から事前に決める (「観測した 2 値のどちらか」や全列共通の最大 ulp 幅を使わない)。
比較尺度はペアに対称で、反復の並べ替えに不変なもの (絶対差か、全ペア共通の分母) にし、比較不能は FAIL にする。完了時に `procedures/verification/` へ移す。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-07 | [2026-10-07-architecture-solver-host-memory-plan.md](../../notes/reviews/2026-10-07-architecture-solver-host-memory-plan.md) | GO-with-changes, C0/M5/m1 (R1・R2 は妥当、R3 は保護と順序の補強、§6 は受入判定できない) | **全件採用** (根拠は自分で確認: `main.cpp:1148` の checkpoint 復元・`main.cpp:3417` の出力登録が確保より後・`procedures/verification/README.md` の case/08 除外・case/63/64 はこのツリーに無い)。M1 → §4 原則 (共通アクセサ) + §5.1 #5 負例、M2 → §4 原則 (順序・共通関数) + §5.1 #2、M3 → §5.1 #3 (基準入力)、M4 → §6 (比較時点)、M5 → §6 (構成表・VERDICT) + §5.1 #2/#3、m1 → §4.4 (工程別モデル、1000 万は暫定)。条件 1 の諮問を兼ねる |

## 7. 影響範囲

`gpu: 1` の全ケースのホスト側確保 (出力・restart・probe・診断・初期化のホスト計算)。数値カーネルとデバイス側の確保は変えない。
`gpu: 0` (CPU 経路) は変えない。

## 8. 完了条件

§6 の全項目が合格し、codex result レビューを経ること。g5 で回せる格子の上限の新しい見積もりを起点 plan (tooling-sern-mesh-blocking §5.1 B4b-2) に返す。

## 未確定事項

- 段階 2 (R4+R5、R6) に進むか — 段階 1 の結果と、SERN 全体の規模のユーザ判断を見て決める。
- level 2 の初期出力に、確保時に 0 初期化されない平均流残差が書かれている (未定義値の出力。2026-10-07 plan レビュー M4 で判明、本 plan の範囲外の既存問題)。
- 監査で見つかった既存の問題 (本 plan の範囲外、挙動は今のまま再現する): (a) `FORGE_OUT_RESIDUALS`/`FORGE_RESID_SNAP` は level < 2 では出力に効かない (printf は「追加した」と言う)、level 2 では `res_ro` 等が `output_cellValNames` と重複し XDMF の Attribute と D2H が二重になる。
  (b) `lineImplicit` はホストの `ccx/ccy/ccz` を読むが GPU 経路では一度も書かれない (cell で `nodes.size() < nCells` のとき 0 座標を読む。node はノード座標を使うので影響なし)。
  (c) `readValueHDF5` はファイル側のデータセット長を検査しない。(d) `sstF1 = 1` の H2D がループ内で繰り返される (時間の無駄、最終状態は同じ)。

## 9. 変更ログ

- `2026-10-07` — 経路の一致で決着させる案: codex は方向を採用、ただし全操作列 (CUPTI)・読戻しの再供給・ポインタの対応付けが要り、通っても `done` 不可。検証の深さをユーザ判断へ。
- `2026-10-07` — 固定幅 A/B は判定不能 (新規 base も超過)。codex: 限定事項つき `done` も却下、物理的許容差が未定義。経路 (H2D・確保) の一致で決着させる案を諮問。
- `2026-10-07` — 追加診断 B: 3402 量で超過 0 だが、7 量で base/new を区別できず、S が new 内の差を含むので変更後だけの異常も吸収しうる (codex diagnose)。`done` を保留し、§6.3 の固定幅の独立 A/B (c44dual_ckpt100、各 6 本、幅は既存 base 3 本だけで凍結) を事前登録。
- `2026-10-07` — 回帰 (各 3 回): 登録判定は (a) step 0 と (b) 4 構成 7 量が FAIL、NaN・初期出力・出力互換・変換器は PASS、SERN g3 ホスト 5157 → 2630 MiB。codex diagnose: 比較器の尺度の欠陥 (非対称な分母・inf を PASS) を確認、登録 FAIL は維持し、絶対 L∞ の追加診断 B を事前登録 (§6.2)。
- `2026-10-07` — R3 を実装 (ローカル: ホスト HWM 2737 → 1371 B/節点、出力・checkpoint 一致、負例 2 つ停止、変換器の出力一致)。回帰ハーネス (`case/66.hostmem_regression/`) を作成、base r1 完走、AWS の idle 自動停止で r2/r3 が中断。
- `2026-10-07` — R1・R2 を実装 (ローカル計測でホスト HWM −510 B/節点、GPU 不変、step 0 と初期出力は一致)。本格回帰は AWS ハーネスで。
- `2026-10-07` — §5.1 #2 監査済み (ホスト `c` の GPU 経路の参照は 5 か所、probe 等は不要、R2 に `pdeSize` の fallback 修正を追加、試験構成を具体化)。回帰は専用 case `case/66.hostmem_regression/` で。
- `2026-10-07` — codex plan レビュー (GO-with-changes, M5/m1) を全件採用: 初期化の順序と出力/確保の共通関数、ホスト参照の共通アクセサ、正常な基準入力、比較時点の固定、構成表と VERDICT、工程別メモリモデル。
- `2026-10-07` — 起票。ユーザ決定「forge 本体のメモリを小さくする余地を、計測から調べる」(tooling-sern-mesh-blocking §5.1 B4b-3) の計測結果を受け、段階 1 (R1–R3) を設計。
