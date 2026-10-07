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
(読むのはビルド対象外の `solvePoisson_amgx.cpp`・`solveNavierStokes.cpp` だけ)、各 wrapper は参照を受け渡すだけ。
`matrix` オブジェクト自体は残す (関数の引数として渡っているため)。

### 4.2 R2 ホストの面変数 `p` を確保しない (見込み −195 B/節点)

`variables::allocVariables` (`variables.cpp`) で `useGPU == 1` のときホスト側の `p[name]` を resize しない (デバイス側 `p_d` は今までどおり)。
確認: ホスト側の `p` を読むのは CPU 経路 (`variables::setStructuralVariables` の `gpu==0` 部、`setStructualVariables.cpp`、`gradient.cpp` の `gpu==0` 部)
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

監査表は **名前・有効になる条件 (設定・環境変数)・最初の利用箇所 (ファイル:行)** を列にする。監査の方法: ビルド対象の全ソースで、ホスト側 `c` への参照 (`c[`・`c.at(`・`.c.count(`・`copyVariables_cell_*` の名前リスト・`output_cellValNames` 系) を列挙し、
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
| 2 | **監査** (M1・M2・M5) | ホスト側 `c`/`p` の全参照を、**名前・有効条件・最初の利用箇所 (ファイル:行)・書くか読むか・保護の要否**の表に (§4.3)。出力と checkpoint の依存名 (`h0` ← `Ht`/`k`、`/CHECKPOINT` 履歴、`FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP` 等の環境変数) を含める。表の各行に対応する**試験 (入力・確認する成果物・判定)** を割り当てる (§6 の構成表の元)。触るファイル: なし (読むだけ)。合格: 表が本 plan に入り、§6 の構成表と対応している | O |
| 3 | **基準入力の確定** (M3・M5) | 発散しない基準入力を固定し、恒久的な `run_*` と case README に残す: SERN 3D 生産設定 = AWS の g3 (`run_1068` の最終場から restart) と g4 (`run_1045`/`run_1046` の最終場) — 2 サイズでメモリの傾きも取る。2D node = `case/36` (手順書の標準)。軸対称・多成分・凝縮 = `case/44`。CHT = `case/52.conjugate_slab`。cell モード = `case/20` か `case/13` (全 run が cell)。dual-time (checkpoint 復元) = 監査で選ぶ。遷移モデル・probe・`output.level: 2`・`extraFields`・環境変数診断は、監査表で到達する構成に割り当てる。**他セッションの case ディレクトリには書かない** (入力を自分の run へ複製)。2026-10-07 の縮小格子 3 点は**発散前までのメモリ観測**として保持し、回帰の基準には使わない (`s050` は step 9 で `roe` 非有限、3 点とも `check_convergence` は NOT CONVERGED) | O |
| 4 | R1–R3 の実装 | §5 の 2〜4。各段階でメモリを測る (g3/g4、`FORGE_MEMLOG=1`) | O |
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

  | 構成 | 入力 | 確認する成果物 |
  | --- | --- | --- |
  | SERN 3D node 生産設定 | AWS g3 (`run_1068` 最終場から restart)、g4 | 残差 CSV・`res_*.h5`・力の時系列・メモリ |
  | 2D node (標準) | `case/36` | 残差・出力 |
  | 軸対称・多成分・凝縮 | `case/44` | 残差・出力・種/モーメント |
  | 共役伝熱 | `case/52.conjugate_slab` | 残差・出力・CHT 出力 |
  | cell モード | `case/20` または `case/13` | 残差・出力 (cell は毎回 `cells[].iNodes` で CONNE を組む) |
  | dual-time の checkpoint | 監査で選ぶ非定常ケース | **連続実行と checkpoint 経由の分割実行**の一致、`/CHECKPOINT` の履歴属性 |
  | `output.level: 2`・`extraFields`・`FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP` | SERN 3D か `case/36` | 出力データセット集合・値 |
  | probe・遷移モデル・その他の環境変数診断 | 監査表で割り当て | 該当成果物 |

- **ガード (負例 2 つ)**: (1) 出力用の変数を H から外したビルド、(2) checkpoint 履歴変数を外したビルドで、いずれも**ホストへの書込み・転送の前に**変数名つきで停止すること (確認後に戻す)。
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

## 9. 変更ログ

- `2026-10-07` — codex plan レビュー (GO-with-changes, M5/m1) を全件採用: 初期化の順序と出力/確保の共通関数、ホスト参照の共通アクセサ、正常な基準入力、比較時点の固定、構成表と VERDICT、工程別メモリモデル。
- `2026-10-07` — 起票。ユーザ決定「forge 本体のメモリを小さくする余地を、計測から調べる」(tooling-sern-mesh-blocking §5.1 B4b-3) の計測結果を受け、段階 1 (R1–R3) を設計。
