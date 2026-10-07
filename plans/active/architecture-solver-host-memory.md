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
| 5 | 回帰 | §6 の全項目。負例 2 つ (出力用変数を H から外す・checkpoint 履歴変数を外す → いずれも書込み/転送の前に名前つきで停止) | O |

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

- `2026-10-07` — R3 を実装 (ローカル: ホスト HWM 2737 → 1371 B/節点、出力・checkpoint 一致、負例 2 つ停止、変換器の出力一致)。回帰ハーネス (`case/66.hostmem_regression/`) を作成、base r1 完走、AWS の idle 自動停止で r2/r3 が中断。
- `2026-10-07` — R1・R2 を実装 (ローカル計測でホスト HWM −510 B/節点、GPU 不変、step 0 と初期出力は一致)。本格回帰は AWS ハーネスで。
- `2026-10-07` — §5.1 #2 監査済み (ホスト `c` の GPU 経路の参照は 5 か所、probe 等は不要、R2 に `pdeSize` の fallback 修正を追加、試験構成を具体化)。回帰は専用 case `case/66.hostmem_regression/` で。
- `2026-10-07` — codex plan レビュー (GO-with-changes, M5/m1) を全件採用: 初期化の順序と出力/確保の共通関数、ホスト参照の共通アクセサ、正常な基準入力、比較時点の固定、構成表と VERDICT、工程別メモリモデル。
- `2026-10-07` — 起票。ユーザ決定「forge 本体のメモリを小さくする余地を、計測から調べる」(tooling-sern-mesh-blocking §5.1 B4b-3) の計測結果を受け、段階 1 (R1–R3) を設計。
