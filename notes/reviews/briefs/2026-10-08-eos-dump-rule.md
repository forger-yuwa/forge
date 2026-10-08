# 諮問: 凍結入力の EOS ビット比較 — 全配列で差分 0 バイトだが TP が規則上 INVALID / B 枝が step の入力にならない (2026-10-08)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §5.1 #2 の「追加の受入れ試験 (事前登録 2026-10-08)」。判断の記録 `notes/reviews/2026-10-08-floor-counter-regression-diagnose.md`。
エスカレーション条件 3 (事前登録の比較が判定不能) と 1 (規則の変更)。記録の写し `notes/investigations/2026-10-08-eos-dump/` (SUMMARY.json・cmp_*.txt・DUMP_SHA256.txt・フックの patch・harness)。

## 実施 (実装担当、commit b1448e54 のフック `FORGE_DUMP_EOS_STEP/FILE` と比較器 `tools/compare_eos_dump.py`)
- フックは `assembleResidualPre` の中で `dependentVariables` を挟む (境界ピン [no-slip・軸・等温壁] と化学種・受動種の原始量の後)。`/pre`・`/post` に var.c_d の全配列 (float32 のまま、nCells_all 長 = 壁・ghost を含む)、`/db` に device の種の熱物性表、属性に EOS が読む設定・経路・blocksize・solverConfig 本文。
- 旧版 = 9f35e3e7 (0ffc1a21 の親) + 同じフック (patch は index・hunk 行以外同一)、新版 = 6717f7f9 + フック、`output.floorEvents: 1`。バイナリ sha256 旧 3fb0e9e4…・新 7ea2a8f9…。フック無効時は HEAD と決定的な比較で一致・ログの追加行 0。
- 入力 7 種 (各入力の sha256 を run ごとに確認): TP 初期場・TP 発達場 (restart_field、保存量ビット一致)・TP の A 枝/B 枝 (小型試験の res_2 に注入後の roe を 1 節点入れて restart_field)・CPG c52cht・CPG の A 枝/B 枝。全入力 step 1 でダンプ。各入力 旧 2 回・新 2 回、6 組比較。
- EOS の書き込み先 (dependentVariables_d.cu から事前に列挙): 速度・TP 分岐 (T・P・ro・roe・Ht・sonic・gamma・cp・Rmix)・CPG 単相 (gamma・cp を除く)・k・ω。読む入力: ro・roUx/y/z・roe・roK・roOmega・T (温度反転の初期推定)・roY。gamma・cp は書くだけで読まない。全 32 run で「EOS の前後で変わった配列 ⊆ 列挙」= 取りこぼしなし。

## 結果
| 入力 | 判定 (6 組すべて) | post 差分バイト / 対象 |
| --- | --- | --- |
| TP 初期場・発達場・A 枝・B 枝 | INVALID | 0 / 5,388,516 |
| CPG c52cht・A 枝・B 枝 | IDENTICAL | 0 / 32,908 |
- 全入力・全組で pre の差分も 0、EOS が読まない配列 (TP 202 本・CPG 172 本) も pre・post とも差分 0。同じ秒に書かれた旧版と新版のダンプはファイル全体でもバイト一致 (tp_init の old_r2 と new_r1 の sha256 が同じ)。
- **TP の INVALID の理由は 1 つ**: pre の gamma・cp が ghost 2082 節点 (= nCells_all − nCells) で NaN。post は全配列が有限。推定 (未検証): 初期化の EOS (main.cpp:1846) が ghost の保存量が 0 の状態で走り (readValueHDF5 は nCells までしか埋めない、variables.cpp:794-805)、境界条件は ghost の ro・roe・T・P を埋めるが gamma・cp は埋めない (ghost の Rmix が厳密に 0、T・P は境界値)。
- **B 枝**: restart を経由すると初期化の EOS が床へ戻し `updateVariablesOuter` (:1876) が commit するので、床の下の状態は step k の EOS の入力にならない (新版の記録: init 行 nT 1 [TP 節点 18259、Δ(ρE) 449.74]、step 1 の eos 行は nT 0・床近傍 1)。得られたのは「床の下の状態に対する初期化の EOS の出力が step 1 の pre に残り、それが旧新でビット一致した」という間接の証拠。

## 問い
1. 有限性の要求を「EOS が読む入力の pre」と「EOS の書き込み先の post」に絞る規則の変更を認めるか (この場合、全 7 入力が IDENTICAL)。読まない ghost の gamma・cp の NaN は別の問題 (初期化の欠落) として記録・別件にしてよいか。それとも ghost の gamma・cp を何かが読んでいないか (境界流束など) を先に確かめるべきか。
2. B 枝 (床の下の入力) の扱い: (a) 初期化の EOS の位置でもダンプする (境界ピンの前の位置)、(b) step k の位置で凍結入力を読み込む (出力専用でなくなる)、(c) 間接の証拠で足りる、のどれか。
3. g3 の継続場 (AWS) は同じ規則で回してよいか。
