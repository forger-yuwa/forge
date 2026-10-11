forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-10-07-hostmem-regression-fail.md`)

# 諮問: ホストメモリ削減 (R1–R3) の回帰で、事前登録した (a) step 0 一致と (b) D ≤ 2S が FAIL — 解釈と、判別のための追加試験 (2026-10-07)

関連 plan: `plans/active/architecture-solver-host-memory.md` (§4 設計、§6 合格条件 [事前登録]、§5.1 #3–#5)。エスカレーション条件 3 (事前登録の比較が FAIL)。
結果の原本: `case/66.hostmem_regression/README.md` (判定・run 一覧)、`case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/*.txt` (構成ごとの比較)、比較スクリプト `case/66.hostmem_regression/compare_runs.py`。
コード: base = 9c9f623c (計測ログのみ追加、挙動は変更前と同一)、new = 93e55957 (R1: `mat_ns` を作らない、R2: ホスト面変数 `p` を確保しない、R3: ホストのセル変数を H だけ確保 + 共通関数 + `hostCell` アクセサ + 初期化の順序変更)。設計上、**デバイス側の確保・転送・カーネルは変えていない**。

## 観測事実
- 30 構成 (forge 25・変換器 5) + SERN g3、base 3 回・new 3 回 (AWS g5、`FORGE_CUDA_BLOCKSIZE=128`、同一入力)。全 run rc 0・NaN/Inf なし。
- **PASS**: (c) NaN、初期出力 (保存量・原始量・幾何量のビット一致、level 2 の 100 データセットを含む)、出力互換 (ファイル集合・データセット集合・shape・dtype・属性)、ログ行 (行の集合で比較)、変換器 5 構成 (base 内・new 内・base 対 new の全データセットで差 0)、dual-time の分割 (100 + checkpoint 再開 100) 対連続 200 (`split_vs_cont_c09.txt`: 保存量はビット一致、トレーサ roXi の m は base 3.31e-6/2.81e-6/6.0e-7、new 2.71e-6/2.76e-6/6.0e-7 = 反復ばらつきと同程度)。再開直後の `/CHECKPOINT` は入力とビット一致 (ローカル確認)。ψ 退避 `snapshot: 269 of 269`。
- **(a) step 0 FAIL** (c20cell_impdiag 以外の全構成): §6 は「既知の 1 ulp の 2 値のどちらか」を前提にしたが、**base 自身が 3 回で 3 値に割れる列がある** (幅: 2D 定常 19–75 ulp、SERN 52 ulp、dual-time は step 0 に内反復 22 行が入り 1e6 ulp 級)。new と最も近い base 値との距離は大半の構成で base 内の幅以下。上回ったのは c44steady (21961 / 13561 ulp)、c52cht (29 / 19)、c26optin (239 / 199)、c26optin_env (146 / 94)、c44dual_restart100 (1.45M / 1.30M)。c52cht は追加で各 4 回回すと base 7 回中 4 回が new と同じ値を出した (`rms_roUx`)。
- **(b) D ≤ 2S FAIL (4 構成・7 量)** — 保存量・原始量 (ro・roUx・roe・P・T・ρY・凝縮モーメント・k/ω・γ) は全構成で D ≤ 2S (D/2S 最大 0.55–0.95、SERN g3 の最悪は `res_vehicle_base_18_100:twall_z` で 1.00 = 境界ちょうど)。FAIL はいずれも名目ゼロか疎な診断量で、base 自身の run 間で桁が動く:
  - c44dual_ckpt100 `condClampCorrQ_0` (S 44.5 / 1.0、D 1.79e5): run ごとの最大 base 8.2e17・8.1e17・3.6e19、new 1.5e23・8.1e17・1.7e17 (c44dual_pindiag の base r3 も 1.45e23)。`condR30_0` D 3.14・S 1.04。
  - c44dual_restart100 `condClampCorrQ_0`: base 1.4e10・8.0e6・1.8e3、new 7.9e6・7.9e6・3.1e2 (D 4.39e3、S 1.0)。
  - c20cell_rk3 `rms_roUz` (擬似 2D の spanwise 残差 1e-9〜3e-7、rms_roUx は 6e-3): D/2S 4.0。
  - c20cell_dual `Uz`・`roUz`・`CHECKPOINT/roUzN`: roUz の最大 base 1.7e-8・4.3e-7・5.7e-9、new 4.2e-7・5.8e-9・5.7e-9 (roUx は 156)。D/2S 1.68。
- メモリ: SERN g3 ホスト VmHWM 5157 → 2630 MiB (2816 → 1436 B/節点、−49 %)、GPU 不変。ローカル縮小格子 2 点の傾き 2737 → 1371 B/節点 (§6 の ≤ 1500 は満たす。g4 は入力消失で g3 単点 + ローカル傾きで判定)。

## 期待値と出典
- plan §6 (2026-10-07 事前登録): (a) step 0 の残差行は全列ビット一致 (1 ulp の 2 値の列は例外)、(b) 各データセットで「base 3 回・new 3 回の同ビルド内ペア差 (計 6 対) の最大 S の 2 倍以内」(D = base 対 new のペア差の最大)。[[repeat-range-as-limit-rule-trap]] を意識して係数 2・両側プール。

## 仮説 (呼び出し側)
- 差は全て run 間の非決定性 (atomicAdd 集計) で、変更起因の系統差は無い。根拠: デバイス側を変えていない設計、初期出力と checkpoint のビット一致、保存量・原始量は全構成で D ≤ 2S、FAIL 量は base 内で桁が動く重い裾の診断量。
- ただし (a) の前提 (1 ulp の 2 値) は事前に誤っていたので、(a) は判別力を持たなかった。

## 問い
1. この FAIL をどう扱うか (登録判定は FAIL のまま記録し、追加の判別試験で結論を出す、が呼び出し側の案)。
2. 判別試験の案と、その事前登録すべき合格条件:
   (i) **初期化直後の全デバイス配列のハッシュ**: env で有効な小さな診断 (例 `FORGE_DEVICE_HASH=1`: 初期化の終わりに全 `c_d`/`p_d`/`bvar_d`/格子マップを D2H して名前ごとのハッシュを出力) を base コミットと new コミットの両方に載せた 2 バイナリで、全構成のハッシュが一致すること (計算は決定的なので一致すべき)。一致すれば、以後の差は非決定性だけ、と言えるか。
   (ii) FAIL した構成・量だけ反復を増やす (例 各 10 回) 統計比較 (順位和検定や max の分布比較)。必要か、(i) で十分か。
   (iii) 1 step 目の組立直後のデバイス状態のハッシュ (atomicAdd の前) は取れるか / 要るか。
3. (a) の基準をどう作り直すべきだったか (今後の回帰の作法として plan と procedures に残すため)。

## 関連 plan 全文 (`plans/active/architecture-solver-host-memory.md`)

```markdown
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
```

## 参考: `case/66.hostmem_regression/README.md`

```
# case/66 — ホストメモリ削減の回帰ハーネス

## ケース概要

plan [architecture-solver-host-memory](../../plans/active/architecture-solver-host-memory.md) §5.1 #3 (基準入力の確定と回帰ハーネス) の実体。
§6 の構成表の各構成を、変更前ビルド (base) と変更後ビルド (new) で**同じ入力・同じ手順で 3 回ずつ**回し、
§6 の (a) step 0 のビット一致・(b) N step 後の差が同ビルド内ばらつきの 2 倍以内・(c) NaN/Inf 無し・出力互換性 (データセット集合・shape・dtype・属性) を
`compare_runs.py` で判定する。

- **base** = `9c9f623c` (AWS `~/bin-hostmem/forge_9c9f623c`、sha256 `c65dbb39f772…`、変換器 `~/bin-hostmem/convert_9c9f623c`)。
- **new** = `93e55957` (R1–R3。AWS `~/bin-hostmem/forge_93e55957` sha256 `7121a8e51e7b…`、変換器 `~/bin-hostmem/convert_93e55957` `144db0f8b20f…`。
  `~/forge-b4` を 93e55957 へ `git checkout -f` し、新しいビルドディレクトリ `solver_density_cuda/build_93e55957` でクリーンビルド。
  base と同じ `RelWithDebInfo`・`CMAKE_CUDA_ARCHITECTURES=86`・`CMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl`)。
  ハーネスはどちらも `.bin/<ビルド名>/forge` (シンボリックリンク) 経由で起動する。
- **run の置き場所は AWS** `~/forge-b4/case/66.hostmem_regression/run_NNNN_<構成>_<ビルド>_r<k>/`。比較も AWS で行い、手元にはテキストだけを持ち帰る
  (ローカルでは forge を回さない。重い h5 を手元へ引かない)。このディレクトリ (ローカル) には入力 config・スクリプト・README だけを置く
  (h5 は `.gitignore` で追跡外)。
- 元入力は別セッションのワークツリー `/home/sano/work/forge/case/...` にあり、**読むだけ**で複製して使う (`prepare_inputs.py`)。
- SERN 3D (g3) は `case/46.sern_design` の base `run_1072`〜`run_1074`・new `run_1075`〜`run_1077` (AWS `~/forge-r8/case/46.sern_design/`、[case/46 README](../46.sern_design/README.md) の run 一覧)。

## ファイル

| ファイル | 役割 |
| --- | --- |
| `matrix_spec.py` | **構成表の正本** (入力の元 run・複製するファイル・設定の修正・種・checkpoint・構成ごとの環境変数) |
| `prepare_inputs.py` | (ローカル) 元 run から `inputs/<入力名>/` を作る。修正は正規表現 + 期待一致数で行い、PyYAML で期待値を検査。`mesh.bndFirstOrder`・`wallTreatmentSST: 1` が残れば失敗 |
| `run_matrix.py` | (AWS) `seed` (restart_field.py を種に掛ける)・`set-ckpt`・`launch` (run を作って 1 本ずつ順に回すワーカーを裏で起動)・`resume`・`status`・`verify` (RUN_PROVENANCE の forge_bin/sha256 照合)・`note` (README の状態・比較から除外)・`table` (下の run 一覧の行) |
| `compare_runs.py` | (AWS) §6 の判定。`--cfg X` / `--all` (registry.tsv から run を選ぶ)、`--base … --new …` (直接指定)、`--diff2 A.h5 B.h5` (分割と連続など 2 ファイル) |
| `memlog_summary.py` | (AWS) `FORGE_MEMLOG=1` の工程別 RSS/HWM/GPU と `--memwatch` の 1 s 採取 (`mem_samples.csv`) を表にする |
| `split_vs_cont.py` | (AWS) dual-time の分割 (ckpt100 → 再開 100) と連続 200 の差を、ビルドごとに反復内の差と並べる (判定はしない) |
| `results/<日付>_<base>_vs_<new>/` | AWS から持ち帰った判定のテキスト (`summary.txt`・構成ごとの報告・`sern_g3.txt`・`sern_g3_memlog.txt`・`split_vs_cont_c09.txt`・`registry.tsv`) |
| `inputs/<入力名>/` | 入力テンプレート。`SOURCE.txt` (元・sha256・修正の全記録)、`FILES` (run に写すファイル)、AWS では `SEEDED.txt` (種の適用記録)・`CKPT_FROM.txt` |

## 構成表

N は入力の `nStepOuter`。出力は `outStepInterval = N` で初期出力 `res_0` と最終 `res_N` の 2 回 (c09cont200 だけ 100 ごと)。2D/小格子だけなので N は全部 100〜200。

| 構成 | 入力 (`inputs/`) ← 元 run (`/home/sano/work/forge/case/` 相対) | 直した点 | N | 環境変数 | 確認する成果物 |
| --- | --- | --- | --- | --- | --- |
| `c36node` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | - | `res_*.h5` (集合・shape・dtype・属性)・残差・壁 h5 `res_wall_4_*`・`point_probe_*.out` |
| `c36node_impdiag` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `diag.csv` |
| `c36node_psidual` | `c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | 旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと) | 200 | `FORGE_DIAG_PSI_DUALEVAL=2,4` | `psi_dualeval.csv`・log `snapshot: X of Y device arrays` (base は 269 of 269) |
| `c09ckpt100` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | - | `/CHECKPOINT` 14 データセットと属性 (layout・nHistoryValid・totalTime) |
| `c09restart100` | `c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | 無修正。`ic_ckpt.h5` = base r1 の c09ckpt100 の res_100.h5 (全ビルド・全反復で同じファイル) | 100 | - | log の履歴復元行・`/CHECKPOINT`・分割 (本構成の res_100) と連続 (c09cont200 の res_200) の差 |
| `c09cont200` | `c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | 無修正 (連続 200 step、res_100/res_200) | 200 | - | res_100・res_200 |
| `c09ckpt100_outres` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_OUT_RESIDUALS=1`, `FORGE_RESID_SNAP=0` | `res_*`・`dq_*`・`res_*_m`・`dq_*_new` の有無 |
| `c09ckpt100_rawdiag` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_SPECIES_RAW_DIAG=1` | `roYraw0/1` の有無 |
| `c09ckpt100_pindiag` | `c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | 無修正 (100 step、res_100 に /CHECKPOINT) | 100 | `FORGE_PIN_DIAG=1` | log のピン診断行 |
| `c44steady` | `c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | 無修正 (定常 active run、res_24000 から restart)。凝縮モーメント rog_0/roQ*_0 は変換直後の nozzle.h5 に無いので 0 で作ってから写す。種の記録が無い場なので restart_field は --force-species、solver は FORGE_ALLOW_UNVERIFIED_SPECIES=1。species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須、現行 solver が拒否する) | 200 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | 種・モーメント (`rog_0`・`roQ*_0`)、`res_outlet_2_*` |
| `c44dual_ckpt100` | `c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須)。run_0376 型 (plan §6) のうち推奨設定に合う run_0468 (run_0376 は dual-time で implicitRelax 0.7・SFR 2 なのに convMethod 0 で check_solver_config が FAIL)。200 step を 100 + 100 に分割 (前半) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `/CHECKPOINT/<cons>P` (26 データセット) |
| `c44dual_restart100` | `c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す。後半。`valueFileName` を ic.h5 → ic_ckpt.h5 (= base r1 の c44dual_ckpt100 の res_100.h5) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | log の履歴復元行 |
| `c44dual_pindiag` | `c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須)。run_0376 型 (plan §6) のうち推奨設定に合う run_0468 (run_0376 は dual-time で implicitRelax 0.7・SFR 2 なのに convMethod 0 で check_solver_config が FAIL)。200 step を 100 + 100 に分割 (前半) | 100 | `FORGE_ALLOW_UNVERIFIED_SPECIES=1`, `FORGE_PIN_DIAG=1` | log のピン診断行 (入口あり) |
| `c52cht` | `c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | 無修正 (res_8000 から restart。CHT の内部状態は res に無いので界面は再初期化) | 200 | - | 壁 h5 (`iface*`)・`conjugate_*.csv` |
| `c36cell` | `c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | 旧キー → `model: sst`、`wallTreatmentSST` 1→0 | 200 | - | 出力・残差 (cell は CONNE を cells から組む) |
| `c20cell_rk3` | `c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | 無修正 (cell RK3、変換時の初期場から) | 200 | - | 出力・残差・壁 h5 |
| `c20cell_dual` | `c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | 旧キー LESorRANS 1 / LESmodel 1 → `model: wale` | 100 | - | 出力・残差・壁 h5 |
| `c20cell_impdiag` | `c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | 旧キー LESorRANS 0 → `model: none` (cell 陰解法。FORGE_IMPLICIT_DIAG_CSV の cell 版) | 200 | `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `diag.csv` |
| `c57lm` | `c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | 無修正 (mesh.h5 は run_0011 の場を種にしたもの = roGamma を読む経路) | 200 | - | `lm*`、log `[variables] transition … read from input` |
| `c57lm_fromsst` | `c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | run_0013 の LM 設定 + run_0012 (SST だけ) の mesh.h5 に res_20000 を写した場 (roGamma が無い = 初期化する経路) | 200 | - | log `… not in input: will be initialised` |
| `c56lineimp` | `c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | 無修正 (res_100000 から restart) | 200 | - | 出力の集合 (extraFields)、log `[lineImplicit] lines=…` |
| `c56extra` | `c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | 無修正 (extraFields に res_ro 等。元は 1 step → 200 step) | 200 | - | extraFields の `res_*` |
| `c48absorb` | `c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | 無修正 (extraFields に roN.. と dq_block_old_*。元は 1 step → 200 step) | 200 | - | extraFields の `roN..`・`dq_block_old_*` |
| `c26optin` | `c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | `wallTreatmentSST` 1→0 (env による名前の登録・除去は壁処理に依らない: variables.cpp:301-356) | 200 | - | データセット集合、警告 `… 確保されていない変数なので無視する` |
| `c26optin_env` | `c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | `wallTreatmentSST` 1→0 (env による名前の登録・除去は壁処理に依らない: variables.cpp:301-356) | 200 | `FORGE_WI_FORCE_DIAG=1`, `FORGE_WF_CLOSURE_DIAG=1`, `FORGE_OMEGA_BUDGET=1`, `FORGE_WF_REP_DIAG=1` | データセット集合 (wi_*・wf_g・omg_*・rep_*) |
| `v36node` | `v36node` ← `36.passive_pseudoshock_control` | c36node の config (node) で mesh/passive_solid.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v36cell` | `v36cell` ← `36.passive_pseudoshock_control` | c36cell の config (cell) で mesh/passive_solid.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v09` | `v09` ← `09.Taylor-Green` | c09ckpt100 の config (node 周期・2 種・トレーサ) で mesh/Taylor-Green.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v52` | `v52` ← `52.conjugate_slab` | case/52 mesh/ の変換 config のまま | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |
| `v44` | `v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | c44steady の config (軸対称・2 種・凝縮) で nozzle.msh を変換 | 変換 | - | 変換結果のデータセット単位比較 (h5 はバイト単位で再現しない) |

**元から差し替えた入力**: §6 の「`run_0376` 型」は `run_0376` 自体が dual-time で `implicitRelax: 0.7`・`speciesFaceReconstruction 2` なのに `convMethod 0`
(`check_solver_config.py` が FAIL) なので、同じ型 (軸対称 node・[MIXDRY, H2O]・凝縮・scheme 1・level 2・inletProfile) で推奨設定に合う
`run_0468_sweep_cflp12_nsub20_float` (PASS) を使った。

**起動のために複製側だけ直したもの** (全部 `inputs/*/SOURCE.txt` に記録):
- 旧キー `LESorRANS`/`LESmodel`/`RANSmodel` (現行 solver が起動時に拒否) → `turbulence.model` (c36node・c36cell・c20cell_dual・c20cell_impl)。
- `wallTreatmentSST: 1` (使用禁止) → 0 (c36node・c36cell・c26optin)。c26 の env 診断は名前の登録・除去が壁処理に依らない (`variables.cpp:301-356`) ので試験の目的は保たれる。
- case/44: `species_db.yaml` から H2O を外した (凝縮 ON では凝縮する H2O は組込みでなければならず、DB の H2O は MW 0.0180153 vs 組込み 0.01801528 で拒否される)。
  元の場に化学種の記録 (`species_hash`) が無いので、restart は `--force-species`、solver は `FORGE_ALLOW_UNVERIFIED_SPECIES=1` (監査の env 一覧で host `c` の名前集合に効かない側)。
  c44steady は変換直後の `nozzle.h5` に凝縮モーメントが無いので `rog_0`・`roQ*_0` を 0 で作ってから写した (13 量をビット一致で移送)。
- 同一メッシュ restart の種: c36node (res_80000)・c36cell (res_60000)・c44steady (res_24000)・c52cht (res_8000)・c56lineimp (res_100000)・
  c57lm_fromsst (run_0012 res_20000 を run_0012 の mesh.h5 へ)。いずれも `restart_field.py` (9c9f623c) の VERDICT OK (保存量ビット一致)。
- dual-time の再開 (c09restart100・c44dual_restart100) は **base r1 の ckpt100 の `res_100.h5` を全ビルド・全反復で同じ checkpoint として使う**
  (`ic_ckpt.h5`、`CKPT_FROM.txt` に sha256)。c09 は `run_0001` (sha256 `04a0f362…`)、c44 は `run_0029` (`1473d282…`)。

## 判定 (`compare_runs.py`、§6 を実装したもの。結果を見てから変えない)

- **(a)** `residual_history.csv` の最初の step の全行を全列ビット一致 (文字列一致) で比べる。base 3 回で値が割れた列は「new の値が base の観測値のどれか」で可。
  **初期出力** (最初の `res_*.h5`) は保存量・原始量・幾何量 (`T0_VALUE` の正規表現 + `/VALUE` 以外の全データセット) だけをビット一致で比べ、
  残差・補正量・勾配・リミッタ・診断は比べない (情報として不一致の名前だけ出す)。
- **(b)** 最終出力 (最終 `res_N.h5`・境界出力・CSV/probe 出力・残差履歴の各列) で m(A,B) = max|A−B| / max|A|。
  S = 同ビルド内ペア (base 3 対 + new 3 対) の最大、D = base×new 9 対の最大。合格は **D ≤ 2·S** (S = 0 なら D = 0)。new が無いときは S_base を出す。
- **(c)** 各 run の `NANCHECK.txt` (ワーカーが残差 CSV の数値列と run が書いた全 h5 の浮動小数データセットを検査)。
- **出力互換**: 出力ファイル集合、全 h5 のデータセット集合・shape・dtype・属性 (値まで) を全 run で比べる。
- **ログ**: 名前登録・checkpoint 復元・遷移初期化・ψ 退避件数 (`snapshot: X of Y`)・警告の行が全 run で一致すること。

**比較スクリプトの対象の直し (2026-10-07、new の比較で誤 FAIL・比較漏れが出た所。§6 の判定条件そのものは変えていない)**:
- ログは**行の集合**で比べる (順序は見ない)。変更後は環境変数の出力登録を確保の前へ移すので、`[FORGE_OUT_RESIDUALS]` 等の行の位置が変わるのは設計どおり (監査 §2 (i))。
- `[memlog]` 行 (FORGE_MEMLOG の計測値。"lineImplicit" で拾われ、RSS の実測値が run ごとに違う) と、変更後が設計どおり新しく出す
  `[variables] host cell arrays (gpu: 1): N of M registered` 行はログ比較から外す (後者の N は下の「結果」に別表で載せる)。
- 空白区切りの CSV (`conjugate_Tw_*.csv`) を列に分けられず比較から漏れていた → `,` が無い表は空白で切る。
- `mem_samples.csv` (ハーネスの採取で forge の出力ではない) は比較しない。
- 判定には使わない情報を足した: (a) の base 内 ulp 幅と new の最近傍 base 値までの ulp 距離、(b) の D/2S 最大のデータセット、FAIL の量の max|A|・max|A−B|。
## 使い方 (変更後のビルドで回す手順)

前提: AWS が `running` (`bash solver_density_cuda/tools/aws_instance.sh status`、start/stop は自分でしない)、`pgrep -x forge` の cwd で他セッションを確認。

```bash
# 0) AWS 側の case/66 は ~/forge-b4 の未追跡ファイル。case/66 を含むコミットへ checkout するときは
#    同名の追跡ファイルと衝突するので `git checkout -f <commit>` にする (run_*・h5 は .gitignore 済みで残る)。
#    変更後のバイナリは ~/bin-hostmem/forge_<hash> に置く。
cd ~/forge-b4/case/66.hostmem_regression
python3 run_matrix.py status                     # running? (中断) / queued を確認
python3 run_matrix.py note <中断した run ...> --state 破棄予定 --text "AWS 停止で中断" --exclude
python3 run_matrix.py resume                     # queued のまま残った run を登録どおりに回す
# 1) base の不足分を足す (各構成 3 本そろうまで。番号は自動で連番)
python3 run_matrix.py launch <構成 ...> --build base --bin ~/bin-hostmem/forge_9c9f623c \
    --convert-bin ~/bin-hostmem/convert_9c9f623c --reps 2 3
# 2) 変更後 (ビルド名は段階ごとに変える: r1 / r2 / r3 / new 等)
python3 run_matrix.py launch $(python3 -c "import matrix_spec as m; print(' '.join(m.CONFIGS))") \
    --build new --bin ~/bin-hostmem/forge_<hash> --convert-bin ~/bin-hostmem/convert_<hash> --reps 1 2 3
python3 run_matrix.py status --build new; python3 run_matrix.py verify
# 3) 判定 (構成ごとの報告 + summary.txt)
python3 compare_runs.py --all --new-build new --out-dir compare_new
python3 compare_runs.py --diff2 run_<c09restart100>/res_100.h5 run_<c09cont200>/res_200.h5   # 分割と連続
# 4) テキストだけ持ち帰る (除外を先に書く)
rsync -az --exclude='*.h5' --exclude='residual_history.csv' --include='*/' --include='*.txt' --include='registry.tsv' \
    --include='summary.txt' --exclude='*' ubuntu@<IP>:forge-b4/case/66.hostmem_regression/ case/66.hostmem_regression/_aws/
python3 run_matrix.py table   # → 下の「計算 run 一覧」を置き換える
```

SERN g3 (メモリ): `case/46` に `run_1072` と同じ作り方で run を作り (sern.h5 はハードリンクで共有・読み取り専用)、
`python3 run_matrix.py launch-dir <run ...> --bin <bin> --label <ビルド名> --env FORGE_MEMLOG=1 --memwatch`、
`python3 memlog_summary.py <run ...> --items` で工程別の表。比較は `compare_runs.py --base <base 3 本> --new <new 3 本>`。

## 結果 (2026-10-07、base 9c9f623c 3 回 対 new 93e55957 3 回)

原本: [`results/2026-10-07_base9c9f623c_vs_new93e55957/`](results/2026-10-07_base9c9f623c_vs_new93e55957/summary.txt)。全 run `VERIFY: PASS` (RUN_PROVENANCE の forge_bin・sha256 が登録どおり)。
短い回帰 run なので収束の主張には使わない (全 forge run の check_convergence は NOT CONVERGED、回帰差の判定だけに使う)。

| 構成 | (a) step 0 | step 0 の ulp (base 内の幅 / new の最近傍距離) | 初期出力 | (b) D ≤ 2S | (b) D/2S 最大のデータセット | (c) NaN | 出力互換 | ログ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `c36node` | FAIL (9 値) | 44 / 13 | PASS | PASS | `res_wall_4_200.h5:VALUE/utau` 0.67 (D 4.35e-04, S 3.27e-04) | PASS | PASS | PASS |
| `c36node_impdiag` | FAIL (18 値) | 75 / 33 | PASS | PASS | `point_probe_2.out:T` 0.75 (D 3.22e-07, S 2.15e-07) | PASS | PASS | PASS |
| `c36node_psidual` | FAIL (15 値) | 42 / 5 | PASS | PASS | `psi_dualeval.csv:S_ApmA` 0.92 (D 1.10e+00, S 5.97e-01) | PASS | PASS | PASS |
| `c09ckpt100` | FAIL (300 値) | 4476783 / 3664948 | PASS | PASS | `res_100.h5:VALUE/dUydx` 0.60 (D 5.19e-01, S 4.32e-01) | PASS | PASS | PASS |
| `c09restart100` | FAIL (315 値) | 3575316 / 3435028 | PASS | PASS | `residual_history.csv:rms_roUx` 0.55 (D 5.97e-02, S 5.46e-02) | PASS | PASS | PASS |
| `c09cont200` | FAIL (316 値) | 3744452 / 3408329 | PASS | PASS | `residual_history.csv:rms_roXi` 0.75 (D 1.31e-07, S 8.75e-08) | PASS | PASS | PASS |
| `c09ckpt100_outres` | FAIL (312 値) | 3952886 / 2146423 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.75 (D 3.00e+00, S 2.00e+00) | PASS | PASS | PASS |
| `c09ckpt100_rawdiag` | FAIL (311 値) | 4208959 / 2727567 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.58 (D 1.75e+00, S 1.50e+00) | PASS | PASS | PASS |
| `c09ckpt100_pindiag` | FAIL (315 値) | 3959569 / 2485060 | PASS | PASS | `res_100.h5:VALUE/res_roUy` 0.67 (D 2.00e+00, S 1.50e+00) | PASS | PASS | PASS |
| `c44steady` | FAIL (54 値) | 13561 / 21961 | PASS | PASS | `res_200.h5:VALUE/res_ro` 0.95 (D 2.73e+00, S 1.43e+00) | PASS | PASS | PASS |
| `c44dual_ckpt100` | FAIL (627 値) | 1850076 / 1176704 | PASS | FAIL (2 量) | `res_100.h5:VALUE/condClampCorrQ_0` 2016.58 (D 1.79e+05, S 4.45e+01) | PASS | PASS | PASS |
| `c44dual_restart100` | FAIL (628 値) | 1301738 / 1454616 | PASS | FAIL (1 量) | `res_100.h5:VALUE/condClampCorrQ_0` 2192.98 (D 4.39e+03, S 1.00e+00) | PASS | PASS | PASS |
| `c44dual_pindiag` | FAIL (626 値) | 2482237 / 1193940 | PASS | PASS | `res_100.h5:CHECKPOINT/rog_0_fctH` 0.87 (D 5.23e-04, S 3.02e-04) | PASS | PASS | PASS |
| `c52cht` | FAIL (18 値) | 19 / 29 | PASS | PASS | `res_wall_bot_3_200.h5:VALUE/ifaceRro` 0.72 (D 4.82e-05, S 3.35e-05) | PASS | PASS | PASS |
| `c36cell` | FAIL (21 値) | 28 / 9 | PASS | PASS | `res_200.h5:VALUE/roUx` 0.60 (D 3.95e-04, S 3.28e-04) | PASS | PASS | PASS |
| `c20cell_rk3` | FAIL (14 値) | 18051 / 7617 | PASS | FAIL (1 量) | `residual_history.csv:rms_roUz` 4.00 (D 1.02e+03, S 1.28e+02) | PASS | PASS | PASS |
| `c20cell_dual` | FAIL (172 値) | 1025677 / 984352 | PASS | FAIL (3 量) | `res_100.h5:VALUE/Uz` 1.68 (D 8.87e+01, S 2.65e+01) | PASS | PASS | PASS |
| `c20cell_impdiag` | PASS | 0 / 0 | PASS | PASS | `residual_history.csv:rms_roUx` 0.61 (D 1.69e-06, S 1.39e-06) | PASS | PASS | PASS |
| `c57lm` | FAIL (72 値) | 38642 / 24679 | PASS | PASS | `res_200.h5:VALUE/lmPgamma` 0.75 (D 2.10e-05, S 1.40e-05) | PASS | PASS | PASS |
| `c57lm_fromsst` | FAIL (54 値) | 21209 / 15991 | PASS | PASS | `residual_history.csv:rms_roOmega` 0.65 (D 1.76e-01, S 1.35e-01) | PASS | PASS | PASS |
| `c56lineimp` | FAIL (27 値) | 48 / 40 | PASS | PASS | `res_200.h5:VALUE/roOmega` 0.85 (D 2.35e-05, S 1.38e-05) | PASS | PASS | PASS |
| `c56extra` | FAIL (27 値) | 66 / 33 | PASS | PASS | `res_200.h5:VALUE/roOmega` 0.76 (D 8.69e-06, S 5.72e-06) | PASS | PASS | PASS |
| `c48absorb` | FAIL (54 値) | 4409 / 3672 | PASS | PASS | `residual_history.csv:rms_roK` 0.68 (D 4.41e-03, S 3.24e-03) | PASS | PASS | PASS |
| `c26optin` | FAIL (21 値) | 199 / 239 | PASS | PASS | `res_200.h5:VALUE/k` 0.57 (D 1.83e-05, S 1.60e-05) | PASS | PASS | PASS |
| `c26optin_env` | FAIL (24 値) | 94 / 146 | PASS | PASS | `res_200.h5:VALUE/roUy` 0.74 (D 1.44e-03, S 9.75e-04) | PASS | PASS | PASS |
| `v36node` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v36cell` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v09` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v52` | - | - | - | PASS | - | PASS | PASS | PASS |
| `v44` | - | - | - | PASS | - | PASS | PASS | PASS |
| `SERN g3 (case/46 run_1072–1077)` | FAIL (45 値) | 52 / 15 | PASS | PASS | `res_vehicle_base_18_100.h5:VALUE/twall_z` 1.00 (D 4.11e-07, S 2.05e-07) | PASS | PASS | PASS |

**まとめ**: 初期出力 (保存量・原始量・幾何量のビット一致)・(c) NaN・出力互換 (ファイル集合・データセット集合・shape・dtype・属性値)・ログ行は**全構成 PASS**。
変換器 5 構成は出力 h5 が base 同士・new 同士・base 対 new とも全データセットで差 0。
ψ 二重評価の退避件数は new 3 回とも `snapshot: 269 of 269 device arrays` (base と同じ = R2 の pdeSize の fallback 修正で面配列が漏れていない)。
dual-time の再開 (c09restart100・c44dual_restart100) は new でも `history restored from ic_ckpt.h5 (/CHECKPOINT: 10 levels …)` と FCT 履歴の復元行が出る。

**(a) step 0 は §6 の文言どおりでは c20cell_impdiag 以外の全構成 FAIL**。§6 の前提「既知の 1 ulp の非決定性で、変更前ビルド同士でも同じ 2 値が出る」が
base 3 回の時点で成り立っていない: 最初の組立の残差でも 3 回で 3 値に割れる列があり、幅は 2D node 定常で 19〜75 ulp、SERN g3 で 52 ulp、
dual-time は step 0 に内反復 (22 行) が入るので 10⁶ ulp 級。new の値と最近傍 base 値の距離は大半の構成で base 内の幅以下、
上回るのは c44steady (21961 / 13561)・c52cht (29 / 19)・c26optin (239 / 199)・c26optin_env (146 / 94)・c44dual_restart100 (1454616 / 1301738)。
c52cht の `rms_roUx` は base 3 回・new 3 回がそれぞれ揃って 1 ulp 違ったので、追加で base・new 各 4 回 (run_0186–0193、比較から除外) を回した:
base 7 回中 4 回が new と同じ値を出す (= 系統差でなく 2 値の非決定性)。判定条件の扱い (ulp 幅で読み替えるか) は plan 側の決定事項 (ここでは変えない)。

**(b) の FAIL (4 構成・7 量)** — いずれも名目上ゼロか疎な診断量で、base 自身の反復の間で桁が動いている:

| 構成 | データセット | S (base / new) | D | run ごとの最大絶対値 (base r1–r3 / new r1–r3) |
| --- | --- | --- | --- | --- |
| c44dual_ckpt100 | `res_100:condClampCorrQ_0` (凝縮モーメントのクランプ補正、非零 ~1400 節点) | 44.5 / 1.0 | 1.79e5 | 8.2e17, 8.1e17, 3.6e19 / **1.5e23**, 8.1e17, 1.7e17 (c44dual_pindiag の base r3 も 1.45e23) |
| c44dual_ckpt100 | `res_100:condR30_0` | 1.04 / 1.0 | 3.14 | 1.6e-5, 1.3e-5, 1.4e-5 / 4.2e-5, 2.5e-5, 1.3e-5 |
| c44dual_restart100 | `res_100:condClampCorrQ_0` | 1.0 / 1.0 | 4.39e3 | 1.4e10, 8.0e6, 1.8e3 / 7.9e6, 7.9e6, 3.1e2 |
| c20cell_rk3 | `residual_history.csv:rms_roUz` (擬似 2D の spanwise 残差、`rms_roUx` 6e-3 に対し 1e-9〜3e-7) | 1.0 / 128 | 1.02e3 | (残差履歴の列) |
| c20cell_dual | `res_100:Uz`・`roUz`・`CHECKPOINT/roUzN` (擬似 2D の spanwise 速度) | 26.5 / 1.3 | 88.7 | roUz: 1.7e-8, **4.3e-7**, 5.7e-9 / **4.2e-7**, 5.8e-9, 5.7e-9 (roUx は 156) |

保存量・原始量の主要データセット (ro・roUx・roe・P・T・ρY・凝縮モーメント・k/ω・γ/Re_θt) は全構成で D ≤ 2S。

**dual-time の分割と連続** ([`split_vs_cont_c09.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957/split_vs_cont_c09.txt)): c09 は一様流 (ro・roUx・roe・P・T・ρY は
分割・連続・反復のどれでもビット一致) で、動くのはトレーサ roXi と名目ゼロの roUy・roUz。roXi の m (分割 vs 連続 / 連続内 / 分割内) は
base 3.31e-6 / 2.81e-6 / 6.0e-7、new 2.71e-6 / 2.76e-6 / 6.0e-7 (checkpoint の roXiP も同程度) — 分割と連続の差は反復のばらつきと同じ大きさで、ビルド間で変わらない。

**SERN g3 のメモリ** (192 万節点、本段設定 100 step、`FORGE_MEMLOG=1`、[`sern_g3_memlog.txt`](results/2026-10-07_base9c9f623c_vs_new93e55957/sern_g3_memlog.txt)。3 回とも工程別に同値):

| 工程 | base RSS / HWM [MiB] | new RSS / HWM [MiB] | GPU 使用 [MiB] (両方) |
| --- | --- | --- | --- |
| readMesh 末尾 (平坦な読込配列が生存) | 2614 / 2614 | 2614 / 2614 | — |
| initMatrix の後 | 2721 / 2721 | 2143 / 2614 | 260 |
| allocVariables の後 | 4864 / 4864 | 2344 / 2614 | 2440 |
| setStructuralVariables の後 | 4865 / **5157** | 2338 / **2630** | 2602 |
| 主ループの後 / 終了 | 5001 / 5157 | 2474 / 2630 | 2982 |

- ホスト VmHWM **5157 → 2630 MiB (2816 → 1436 B/節点、−49 %)**。終了時の内訳 (new): host `c` 23/223 本 183 MB (base 1774)、host `p` 0 (357)、`mat_ns` 0 (572)。
  ピークは readMesh 末尾 (2614) のすぐ上の setStructuralVariables (2630) に移った。
- GPU は工程別に base と同値 (プロセス 2974 MiB、初回 memlog から +2724 MB) → §6 の「GPU 傾き ±2 % 以内」を満たす。
- ローカルの縮小格子 s070/s085 からの外挿 (HWM 傾き 2737 → 1371 B/節点、切片 163 → 133 MiB) では g3 で base 5175・new 2644 MiB の見込みで、
  実測 5157・2630 MiB と −0.3 %・−0.5 % で合う。§6 の「ホスト VmHWM の傾き ≤ 1500 B/節点」は g3 単点の平均 1436 B/節点 (切片込み) で下回る
  (傾きの判定は 2 格子が要り、g4 が無いので縮小格子の傾き 1371 を併用)。

**new の host cell 確保数** (`[variables] host cell arrays (gpu: 1): N of M registered`、new r1 の起動ログ):
c36node 19/191・impdiag 20/191・c09ckpt100 98/241 (outres 114/241・rawdiag 100/243・pindiag 102/241)・c44steady 120/303・c44dual 133/303 (pindiag 140/303)・
c52cht・c36cell・c20cell_rk3 19/191・c20cell_dual 26/191・c57lm 92/211・c57lm_fromsst 24/202・c56lineimp 23/191・c56extra 25/191・c48absorb 29/191・
c26optin 20/191 (env 23/212)・SERN g3 23/223。

## 既知の注意

- **AWS の自動停止 (2026-10-07 01:05 UTC 頃 = JST 10:05、base 投入中に停止した)**: idle 自動停止 (`solver_density_cuda/tools/cloud/idle_autostop.sh`、5 分おき・6 回連続で停止) は
  「GPU 使用率 0・`pgrep -x forge` 0 件・ログイン 0・load < 1」を idle と数える。base は `forge_9c9f623c` という名前で起動していたので
  **`pgrep -x forge` に一致せず**、短い 2D run の合間 (GPU 使用率の瞬間値 0) と準備作業の時間で 30 分が数えられたと推定する (停止の記録は見ていない)。
  対策として `run_matrix.py` は `--bin` を `.bin/<ビルド名>/forge` というシンボリックリンク経由で起動するようにした (comm が `forge` になる)。
  再起動後の最初の投入で `pgrep -x forge` が自分の run に一致する (comm=forge、exe は実ファイル) ことを確認した (2026-10-07)。
- **`~/forge-b4` の `git checkout -f` で `solver_density_cuda/third_party/HighFive` が空になった** (この worktree の HighFive は
  `.git` ファイルがローカル PC のパスを指す複製だった)。同じ pin (`dfc06537`) の `~/forge-sern` から中身を複製して戻した
  (`third_party/HighFive.COPIED_FROM.txt`)。次に checkout するときも空になったら同じ手順で戻す。
- 変換器は終了時の `cudaFree` で GPUassert (invalid argument) を出して exit 1 になる既知の罠がある (出力 h5 は完全)。ワーカーは rc でなく
  出力 h5 が開けて `/VALUE`・`/MESH` を持つかで完了を判定する (`written=1`)。cell モードは `/VIZMESH` を書かない。
- case/36 の c36node は 80223 節点・200 step で約 0.3 s、全 2D 構成で 1 本数秒〜十数秒。ワーカーは 1 本ずつ順に回す。
- `run_case.sh` は `~/forge-b4` のものを使うので `RUN_PROVENANCE.txt` の `git_head` は `~/forge-b4` の HEAD (base では 9c9f623c)。

## 保留

- **SERN g4**: `run_1045`/`run_1046` (と延長 `run_1048`/`run_1049`) の最終場・g4 の `sern.h5` は AWS に残っていない (場の h5 は削除済み、
  g4 メッシュも無い)。g4 を使うには格子の再生成 (`problem_3d_prod_m6on_g4.yaml`) と g3 からの cross-mesh restart か段階起動が要る。
  2 格子での傾きは 2026-10-07 の縮小格子 3 点 (`notes/investigations/2026-10-07-forge-memlog/`) でしか出せない。
- **NaN ダンプ (表の P)**: 入力 (s050 縮小格子) が別セッションの scratch にあるので入れていない。
- case/13 `run_slau` (cell) は入れていない (cell は case/36・case/20 の 3 構成で代表)。
- (a) step 0 の判定条件の扱いと (b) の FAIL 7 量 (名目ゼロ・疎な診断量) の扱いは plan の決定事項 (上の「結果」)。

## 計算 run 一覧

AWS `~/forge-b4/case/66.hostmem_regression/` (2026-10-07)。base = `.bin/base/forge` → `~/bin-hostmem/forge_9c9f623c` (r1 の 34 本だけは実ファイル名のまま起動)、
new = `.bin/new/forge` → `~/bin-hostmem/forge_93e55957`、全部 `FORGE_CUDA_BLOCKSIZE=128`。N step の短い回帰 run なので収束判定は NOT CONVERGED
(回帰差の基準であって収束の主張には使わない)。判定済みの run からは入力 h5 の複製を削除した (`INPUTS_REMOVED.txt`、テンプレートは `inputs/` に残る)。
出力は `res_0`・最終 `res_N`・境界出力・CSV だけなので消していない。表は `python3 run_matrix.py table` の出力。

| run | 目的・設定差分 | 元の入力 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0001_c09ckpt100_base_r1` | dual-time 100 step (checkpoint 書出し) (base r1) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0002_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0003_c36node_base_r1` | 2D node 標準 (base r1) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0004_c36node_impdiag_base_r1` | 2D node + 陰解法診断 CSV (表 M) (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0005_c36node_psidual_base_r1` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r1) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0006_c09cont200_base_r1` | dual-time 連続 200 step (base r1) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0007_c09ckpt100_outres_base_r1` | dual-time + 残差出力 (表 K・L) (base r1) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0008_c09ckpt100_rawdiag_base_r1` | dual-time + roYraw (表 I) (base r1) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0009_c09ckpt100_pindiag_base_r1` | dual-time + ピン診断 (表 N) (base r1) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0010_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0011_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 1, NaN PASS, ?。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0012_c52cht_base_r1` | 共役伝熱 (base r1) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0013_c36cell_base_r1` | cell 定常 SST (base r1) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0014_c20cell_rk3_base_r1` | cell RK3 (base r1) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0015_c20cell_dual_base_r1` | cell LES dual-time (base r1) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0016_c20cell_impdiag_base_r1` | cell 陰解法 + 診断 CSV (base r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0017_c57lm_base_r1` | 遷移 (roGamma を読む) (base r1) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0018_c57lm_fromsst_base_r1` | 遷移 (SST の場から初期化) (base r1) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0019_c56lineimp_base_r1` | line-implicit + extraFields (base r1) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0020_c56extra_base_r1` | extraFields (res_* を含む) + FP64 アキュムレータ (base r1) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0021_c48absorb_base_r1` | extraFields (roN・dq_block_old_*) (base r1) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0022_c26optin_base_r1` | env 診断なし (名前が消えて警告) (base r1) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0023_c26optin_env_base_r1` | env 診断 4 種 ON (表 I) (base r1) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0024_v36node_base_r1` | 変換器 node (base r1) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0025_v36cell_base_r1` | 変換器 cell (base r1) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0026_v09_base_r1` | 変換器 node 周期・種 (base r1) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0027_v52_base_r1` | 変換器 node (CHT 用スラブ) (base r1) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0028_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS。起動拒否: species_db.yaml の H2O が組込みと違い凝縮 ON で拒否 (入力修正前) | 破棄予定 |
| `run_0029_c44dual_ckpt100_base_r1` | 軸対称・凝縮 dual-time 前半 (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0030_c44steady_base_r1` | 軸対称・多成分・凝縮 (定常) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0031_c44dual_pindiag_base_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0032_v44_base_r1` | 変換器 軸対称・種・凝縮 (base r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0033_c09restart100_base_r1` | dual-time 再開 100 step (checkpoint 復元) (base r1) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0034_c44dual_restart100_base_r1` | 軸対称・凝縮 dual-time 後半 (再開) (base r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0035_c09ckpt100_base_r2` | dual-time 100 step (checkpoint 書出し) (base r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | done rc=? (中断) nan=-, NOT CONVERGED。AWS 停止 (2026-10-07 01:05 UTC 頃) で中断 | 破棄予定 |
| `run_0036_c09restart100_base_r2` | dual-time 再開 100 step (checkpoint 復元) (base r2) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0037_c09cont200_base_r2` | dual-time 連続 200 step (base r2) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0038_c44dual_ckpt100_base_r2` | 軸対称・凝縮 dual-time 前半 (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0039_c44dual_restart100_base_r2` | 軸対称・凝縮 dual-time 後半 (再開) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0040_c36node_base_r2` | 2D node 標準 (base r2) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0041_c36node_impdiag_base_r2` | 2D node + 陰解法診断 CSV (表 M) (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0042_c36node_psidual_base_r2` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r2) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0043_c09ckpt100_outres_base_r2` | dual-time + 残差出力 (表 K・L) (base r2) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0044_c09ckpt100_rawdiag_base_r2` | dual-time + roYraw (表 I) (base r2) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0045_c09ckpt100_pindiag_base_r2` | dual-time + ピン診断 (表 N) (base r2) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0046_c44steady_base_r2` | 軸対称・多成分・凝縮 (定常) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0047_c44dual_pindiag_base_r2` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0048_c52cht_base_r2` | 共役伝熱 (base r2) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0049_c36cell_base_r2` | cell 定常 SST (base r2) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0050_c20cell_rk3_base_r2` | cell RK3 (base r2) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0051_c20cell_dual_base_r2` | cell LES dual-time (base r2) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0052_c20cell_impdiag_base_r2` | cell 陰解法 + 診断 CSV (base r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0053_c57lm_base_r2` | 遷移 (roGamma を読む) (base r2) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0054_c57lm_fromsst_base_r2` | 遷移 (SST の場から初期化) (base r2) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0055_c56lineimp_base_r2` | line-implicit + extraFields (base r2) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0056_c56extra_base_r2` | extraFields (res_* を含む) + FP64 アキュムレータ (base r2) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0057_c48absorb_base_r2` | extraFields (roN・dq_block_old_*) (base r2) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0058_c26optin_base_r2` | env 診断なし (名前が消えて警告) (base r2) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0059_c26optin_env_base_r2` | env 診断 4 種 ON (表 I) (base r2) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0060_v36node_base_r2` | 変換器 node (base r2) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0061_v36cell_base_r2` | 変換器 cell (base r2) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0062_v09_base_r2` | 変換器 node 周期・種 (base r2) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0063_v52_base_r2` | 変換器 node (CHT 用スラブ) (base r2) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0064_v44_base_r2` | 変換器 軸対称・種・凝縮 (base r2) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0065_c09ckpt100_base_r3` | dual-time 100 step (checkpoint 書出し) (base r3) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0066_c09restart100_base_r3` | dual-time 再開 100 step (checkpoint 復元) (base r3) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0067_c09cont200_base_r3` | dual-time 連続 200 step (base r3) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0068_c44dual_ckpt100_base_r3` | 軸対称・凝縮 dual-time 前半 (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0069_c44dual_restart100_base_r3` | 軸対称・凝縮 dual-time 後半 (再開) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0070_c36node_base_r3` | 2D node 標準 (base r3) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0071_c36node_impdiag_base_r3` | 2D node + 陰解法診断 CSV (表 M) (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0072_c36node_psidual_base_r3` | 2D node + ψ 二重評価 (R2 の pdeSize) (base r3) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0073_c09ckpt100_outres_base_r3` | dual-time + 残差出力 (表 K・L) (base r3) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0074_c09ckpt100_rawdiag_base_r3` | dual-time + roYraw (表 I) (base r3) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0075_c09ckpt100_pindiag_base_r3` | dual-time + ピン診断 (表 N) (base r3) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0076_c44steady_base_r3` | 軸対称・多成分・凝縮 (定常) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0077_c44dual_pindiag_base_r3` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (base r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0078_c52cht_base_r3` | 共役伝熱 (base r3) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0079_c36cell_base_r3` | cell 定常 SST (base r3) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0080_c20cell_rk3_base_r3` | cell RK3 (base r3) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0081_c20cell_dual_base_r3` | cell LES dual-time (base r3) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0082_c20cell_impdiag_base_r3` | cell 陰解法 + 診断 CSV (base r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0083_c57lm_base_r3` | 遷移 (roGamma を読む) (base r3) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0084_c57lm_fromsst_base_r3` | 遷移 (SST の場から初期化) (base r3) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0085_c56lineimp_base_r3` | line-implicit + extraFields (base r3) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0086_c56extra_base_r3` | extraFields (res_* を含む) + FP64 アキュムレータ (base r3) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0087_c48absorb_base_r3` | extraFields (roN・dq_block_old_*) (base r3) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0088_c26optin_base_r3` | env 診断なし (名前が消えて警告) (base r3) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0089_c26optin_env_base_r3` | env 診断 4 種 ON (表 I) (base r3) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0090_v36node_base_r3` | 変換器 node (base r3) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0091_v36cell_base_r3` | 変換器 cell (base r3) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0092_v09_base_r3` | 変換器 node 周期・種 (base r3) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0093_v52_base_r3` | 変換器 node (CHT 用スラブ) (base r3) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0094_v44_base_r3` | 変換器 軸対称・種・凝縮 (base r3) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0095_c09ckpt100_base_r2` | dual-time 100 step (checkpoint 書出し) (base r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0096_c36node_new_r1` | 2D node 標準 (new r1) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0097_c36node_impdiag_new_r1` | 2D node + 陰解法診断 CSV (表 M) (new r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0098_c36node_psidual_new_r1` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r1) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0099_c09ckpt100_new_r1` | dual-time 100 step (checkpoint 書出し) (new r1) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0100_c09restart100_new_r1` | dual-time 再開 100 step (checkpoint 復元) (new r1) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0101_c09cont200_new_r1` | dual-time 連続 200 step (new r1) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0102_c09ckpt100_outres_new_r1` | dual-time + 残差出力 (表 K・L) (new r1) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0103_c09ckpt100_rawdiag_new_r1` | dual-time + roYraw (表 I) (new r1) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0104_c09ckpt100_pindiag_new_r1` | dual-time + ピン診断 (表 N) (new r1) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0105_c44steady_new_r1` | 軸対称・多成分・凝縮 (定常) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0106_c44dual_ckpt100_new_r1` | 軸対称・凝縮 dual-time 前半 (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0107_c44dual_restart100_new_r1` | 軸対称・凝縮 dual-time 後半 (再開) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0108_c44dual_pindiag_new_r1` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r1) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0109_c52cht_new_r1` | 共役伝熱 (new r1) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0110_c36cell_new_r1` | cell 定常 SST (new r1) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0111_c20cell_rk3_new_r1` | cell RK3 (new r1) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0112_c20cell_dual_new_r1` | cell LES dual-time (new r1) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0113_c20cell_impdiag_new_r1` | cell 陰解法 + 診断 CSV (new r1) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0114_c57lm_new_r1` | 遷移 (roGamma を読む) (new r1) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0115_c57lm_fromsst_new_r1` | 遷移 (SST の場から初期化) (new r1) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0116_c56lineimp_new_r1` | line-implicit + extraFields (new r1) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0117_c56extra_new_r1` | extraFields (res_* を含む) + FP64 アキュムレータ (new r1) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0118_c48absorb_new_r1` | extraFields (roN・dq_block_old_*) (new r1) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0119_c26optin_new_r1` | env 診断なし (名前が消えて警告) (new r1) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0120_c26optin_env_new_r1` | env 診断 4 種 ON (表 I) (new r1) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0121_v36node_new_r1` | 変換器 node (new r1) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0122_v36cell_new_r1` | 変換器 cell (new r1) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0123_v09_new_r1` | 変換器 node 周期・種 (new r1) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0124_v52_new_r1` | 変換器 node (CHT 用スラブ) (new r1) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0125_v44_new_r1` | 変換器 軸対称・種・凝縮 (new r1) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0126_c36node_new_r2` | 2D node 標準 (new r2) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0127_c36node_impdiag_new_r2` | 2D node + 陰解法診断 CSV (表 M) (new r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0128_c36node_psidual_new_r2` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r2) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0129_c09ckpt100_new_r2` | dual-time 100 step (checkpoint 書出し) (new r2) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0130_c09restart100_new_r2` | dual-time 再開 100 step (checkpoint 復元) (new r2) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0131_c09cont200_new_r2` | dual-time 連続 200 step (new r2) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0132_c09ckpt100_outres_new_r2` | dual-time + 残差出力 (表 K・L) (new r2) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0133_c09ckpt100_rawdiag_new_r2` | dual-time + roYraw (表 I) (new r2) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0134_c09ckpt100_pindiag_new_r2` | dual-time + ピン診断 (表 N) (new r2) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0135_c44steady_new_r2` | 軸対称・多成分・凝縮 (定常) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0136_c44dual_ckpt100_new_r2` | 軸対称・凝縮 dual-time 前半 (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0137_c44dual_restart100_new_r2` | 軸対称・凝縮 dual-time 後半 (再開) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0138_c44dual_pindiag_new_r2` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r2) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0139_c52cht_new_r2` | 共役伝熱 (new r2) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0140_c36cell_new_r2` | cell 定常 SST (new r2) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0141_c20cell_rk3_new_r2` | cell RK3 (new r2) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0142_c20cell_dual_new_r2` | cell LES dual-time (new r2) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0143_c20cell_impdiag_new_r2` | cell 陰解法 + 診断 CSV (new r2) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0144_c57lm_new_r2` | 遷移 (roGamma を読む) (new r2) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0145_c57lm_fromsst_new_r2` | 遷移 (SST の場から初期化) (new r2) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0146_c56lineimp_new_r2` | line-implicit + extraFields (new r2) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0147_c56extra_new_r2` | extraFields (res_* を含む) + FP64 アキュムレータ (new r2) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0148_c48absorb_new_r2` | extraFields (roN・dq_block_old_*) (new r2) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0149_c26optin_new_r2` | env 診断なし (名前が消えて警告) (new r2) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0150_c26optin_env_new_r2` | env 診断 4 種 ON (表 I) (new r2) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0151_v36node_new_r2` | 変換器 node (new r2) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0152_v36cell_new_r2` | 変換器 cell (new r2) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0153_v09_new_r2` | 変換器 node 周期・種 (new r2) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0154_v52_new_r2` | 変換器 node (CHT 用スラブ) (new r2) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0155_v44_new_r2` | 変換器 軸対称・種・凝縮 (new r2) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0156_c36node_new_r3` | 2D node 標準 (new r3) | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0157_c36node_impdiag_new_r3` | 2D node + 陰解法診断 CSV (表 M) (new r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0158_c36node_psidual_new_r3` | 2D node + ψ 二重評価 (R2 の pdeSize) (new r3) + `FORGE_DIAG_PSI_DUALEVAL=2,4` | `inputs/c36node` ← `36.passive_pseudoshock_control/run_sym_H_2up_node` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0159_c09ckpt100_new_r3` | dual-time 100 step (checkpoint 書出し) (new r3) | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0160_c09restart100_new_r3` | dual-time 再開 100 step (checkpoint 復元) (new r3) | `inputs/c09restart100` ← `09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0161_c09cont200_new_r3` | dual-time 連続 200 step (new r3) | `inputs/c09cont200` ← `09.Taylor-Green/run_0162_passiveG_fct_cont200` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0162_c09ckpt100_outres_new_r3` | dual-time + 残差出力 (表 K・L) (new r3) + `FORGE_OUT_RESIDUALS=1 FORGE_RESID_SNAP=0` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0163_c09ckpt100_rawdiag_new_r3` | dual-time + roYraw (表 I) (new r3) + `FORGE_SPECIES_RAW_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0164_c09ckpt100_pindiag_new_r3` | dual-time + ピン診断 (表 N) (new r3) + `FORGE_PIN_DIAG=1` | `inputs/c09ckpt100` ← `09.Taylor-Green/run_0160_passiveG_fct_ckpt100` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0165_c44steady_new_r3` | 軸対称・多成分・凝縮 (定常) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44steady` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0166_c44dual_ckpt100_new_r3` | 軸対称・凝縮 dual-time 前半 (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0167_c44dual_restart100_new_r3` | 軸対称・凝縮 dual-time 後半 (再開) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1` | `inputs/c44dual_restart100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0168_c44dual_pindiag_new_r3` | 軸対称・凝縮 dual-time + ピン診断 (入口あり) (new r3) + `FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_PIN_DIAG=1` | `inputs/c44dual_ckpt100` ← `44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0169_c52cht_new_r3` | 共役伝熱 (new r3) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0170_c36cell_new_r3` | cell 定常 SST (new r3) | `inputs/c36cell` ← `36.passive_pseudoshock_control/run_sym_H_2up_cell` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0171_c20cell_rk3_new_r3` | cell RK3 (new r3) | `inputs/c20cell_rk3` ← `20.naca_ml/001.test/run_slau` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0172_c20cell_dual_new_r3` | cell LES dual-time (new r3) | `inputs/c20cell_dual` ← `20.naca_ml/001.test/run_case04_les_unsteady_dualtime` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0173_c20cell_impdiag_new_r3` | cell 陰解法 + 診断 CSV (new r3) + `FORGE_IMPLICIT_DIAG_CSV=diag.csv` | `inputs/c20cell_impl` ← `20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0174_c57lm_new_r3` | 遷移 (roGamma を読む) (new r3) | `inputs/c57lm` ← `57.transition_flat_plate/run_0014_t3a_lm_unitcheck` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0175_c57lm_fromsst_new_r3` | 遷移 (SST の場から初期化) (new r3) | `inputs/c57lm_fromsst` ← `57.transition_flat_plate/run_0013_t3b_lm` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0176_c56lineimp_new_r3` | line-implicit + extraFields (new r3) | `inputs/c56lineimp` ← `56.gap_tp1187/run_0019_lineimplicit` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0177_c56extra_new_r3` | extraFields (res_* を含む) + FP64 アキュムレータ (new r3) | `inputs/c56extra` ← `56.gap_tp1187/run_0027_s6_f32_b` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0178_c48absorb_new_r3` | extraFields (roN・dq_block_old_*) (new r3) | `inputs/c48absorb` ← `48.flat_plate_cooled_m4/run_0903_absorb2` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0179_c26optin_new_r3` | env 診断なし (名前が消えて警告) (new r3) | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0180_c26optin_env_new_r3` | env 診断 4 種 ON (表 I) (new r3) + `FORGE_WI_FORCE_DIAG=1 FORGE_WF_CLOSURE_DIAG=1 FORGE_OMEGA_BUDGET=1 FORGE_WF_REP_DIAG=1` | `inputs/c26optin` ← `26.flat_plate_sst/run_0086_optin_base` | rc 0, NaN PASS, NOT CONVERGED | active |
| `run_0181_v36node_new_r3` | 変換器 node (new r3) | `inputs/v36node` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0182_v36cell_new_r3` | 変換器 cell (new r3) | `inputs/v36cell` ← `36.passive_pseudoshock_control` | rc 1, NaN PASS | active |
| `run_0183_v09_new_r3` | 変換器 node 周期・種 (new r3) | `inputs/v09` ← `09.Taylor-Green` | rc 1, NaN PASS | active |
| `run_0184_v52_new_r3` | 変換器 node (CHT 用スラブ) (new r3) | `inputs/v52` ← `52.conjugate_slab` | rc 1, NaN PASS | active |
| `run_0185_v44_new_r3` | 変換器 軸対称・種・凝縮 (new r3) | `inputs/v44` ← `44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX` | rc 1, NaN PASS | active |
| `run_0186_c52cht_base_r4` | 共役伝熱 (base r4) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0187_c52cht_base_r5` | 共役伝熱 (base r5) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0188_c52cht_base_r6` | 共役伝熱 (base r6) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0189_c52cht_base_r7` | 共役伝熱 (base r7) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0190_c52cht_new_r4` | 共役伝熱 (new r4) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0191_c52cht_new_r5` | 共役伝熱 (new r5) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0192_c52cht_new_r6` | 共役伝熱 (new r6) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
| `run_0193_c52cht_new_r7` | 共役伝熱 (new r7) | `inputs/c52cht` ← `52.conjugate_slab/run_0007_fxhalf` | rc 0, NaN PASS, NOT CONVERGED。ノイズ確認の追加反復 (step 0 の rms_roUx が base でも 2 値に割れることの確認、判定に使わない) | ref |
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
