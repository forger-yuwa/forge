# 諮問ブリーフ: 移植の回帰 R1・R3 (前処理なしの値が変わらない) の回し直しの設計 (2026-10-13)

- 依頼者: 主セッション (エスカレーション条件 3: 事前登録の比較が FAIL、条件 4: 回し直しの前)
- plan: `plans/active/architecture-float-state-double-geometry.md` §5.1 #26・§6.59 と結果。統合の記録 `notes/investigations/2026-10-13-main-integration/README.md` (確認 6 が R1〜R3)
- ユーザ: main へ統合してよい (2026-10-13)。計算を速くしたい。

## 観測事実
- §6.59 の R1 (case/45 G0、point と B0、10 step、float と FP64) と R3 (case/48 平板 CPG、10 step、float と FP64): 各組で移植の前 2 本 (a・b、同じバイナリ) と後 1 本。登録の基準「全量で max|後 − 前a| ≤ max|前b − 前a|」で全組が不合格。例: R1 B0 float の ρv 5.14e-3 対 2.64e-3、ω 2.6e5 対 2.0e5、R1 B0 FP64 の h0 0.1944 対 0.1943・dt_local 1.282e-12 対 1.282e-12 (同値で >)、R1 point FP64 の P 0.5075 対 0.5033、R3 FP64 の P 6.9e-6 対 5.0e-6・T 3.0e-7 対 2.1e-7。いずれも前の再実行の差と同じ桁で 1〜2 倍。
- 前の 2 本はどの組でもビット一致しない (node の残差の atomicAdd の加算順)。ソルバに加算順を固定する設定は無い。
- 移植 (`e224ae11`) の実装の確認 (implementer): cuda_forge の全 .cu.o を関数ごとに SASS で比べ、float は 296 関数・FP64 は 292 関数が命令列まで同一。違うのは前処理のカーネル (`lowMachPrecond ≥ 2` だけが使う) と新しいカーネル `setDTlocal_species_beta_d` (同じく ≥ 2 だけ) と、FP64 の `dependentVariables_d` の PTX の 1 か所 (NaN 定数の mov の位置。ソースのパスに依存し、移植後のヘッダを移植前のパスでコンパイルすると移植前と一致)。今回の回帰の前後のバイナリは同じパス (`~/forge-port`) でビルドした。`solverConfig.hpp` の構造体に `int speciesPrecondDt` が足されて配置が変わった。`dt_local_sp` は多成分のときだけ確保 (case/45 は 2 成分なので確保される、R3 の case/48 は CPG で確保されない)。
- メモリの記録 (repeat-range-as-limit-rule-trap、2026-10-03): 「新旧の差の最大 ≤ 旧の再実行の差の最大」型は同じ分布でも多くが FAIL する。codex の 2026-10-03 の提案は「群の差 D + 群の中のばらつき E_A + E_B ≤ 許容」か、物理量の許容を先に決める形。

## 問い
1. R1・R3 の不合格をどう扱うか (登録の判定は読み替えない前提で)。
2. 回し直すならどの設計がよいか。案:
   - Q1 (群の方式): 各組で前 3 本・後 3 本。各量 (VALUE の全量・残差の rms 列) について、群の平均の差の最大 D と、群の中のばらつき (各群の中の全組の差の最大) E_A・E_B を取り、D ≤ E_A + E_B なら同等。
   - Q2 (物理の許容): 10 step 後の相対差の許容を量ごとに先に決める (例: 保存量・P・T の相対 1e-6 (FP64)・1e-4 (float))。
   - Q3 (コードの同等性を主にする): 前処理なしの経路で起動されるカーネルの SASS の同一性 (済み) と、構造体の配置の変化が読み違いを起こさないこと (設定の読み込みの単体試験 `test_check_solver_config.py` は PASS) で R1・R3 の目的 (値が変わらない) を満たすとし、run の比較は群の方式で補う。
3. 統合 (main) は、R1・R3 の回し直しの結果を待つべきか。

## 読んでよいファイル
- float plan §5.1 #26、§6.57、§6.59 と結果、`case/45.isobutane_m6_d155/port_reg.sh`・`port_reg_an.py`・`_band_ab/cold_pair/port_reg_judge.json`
- `notes/investigations/2026-10-13-main-integration/README.md`
- `solver_density_cuda/tests/unit/test_precond_node_rows.cu`
