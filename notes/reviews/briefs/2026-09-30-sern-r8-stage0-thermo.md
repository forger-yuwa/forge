# 諮問: SERN R8 回帰の照合 (0) が事前基準をわずかに超えた件 (2026-09-30、条件 3)

plan: `plans/active/tooling-nozzle-sern-chain.md` の R8 行 (事前登録: 回帰 (i) の (0)・(1))。
コード: `feature/sern-design` 772cecd8 (マージ 237357fd = feature/species-transport 36685217、ランナー 5eee3ff2)。
AWS: 旧 = `~/forge-pgrad-new` (ランナー c758e59d、バイナリ build-ff = a6ceee0b)、新 = `~/forge-r8` (8b87c0e7 のビルド)。

## 観測事実
- (0) 熱物性の照合 (事前基準: MW・NASA-9 係数が相対 1e-15 以内):
  旧ランナーが書いた `run_1002_r8A_1/species_db.yaml` の EXH/AIR と、新ランナーで準備した run (`/tmp/r8chk`) を `forge_species.run_thermo` (新バイナリの `--resolve-species`) で解決した EXH/AMB を比較。
  - AMB vs AIR: MW 相対 1.2e-16、係数は全て完全一致。
  - EXH: MW 相対 2.8e-16、係数は 18 個中 10 個が 1〜5 ulp 違い、最大 `nasa9_low[0]` 1560.6275407330086 (新) vs 1560.627540733007 (旧) = 相対 **1.02e-15** (基準超え)。次は `nasa9_low[8]` 6.1e-16、`nasa9_low[1]` 4.5e-16。
  - 旧は設計側 Python (`composition.lump_entry` 相当) が合成して書いた値、新はソルバ (C++ `speciesDB.cpp` の lump 合成) が起動時に合成した値。構成種は 11 種 (N2, H2O, H2, AR, OH, O2, NO, H, O, CO2, CO)、lump 内モル分率は同じ全桁の値 (config に repr で書く)。
- (1) 回帰 (事前基準: 全 run GATES PASS かつ 4 量とも |平均_B − 平均_A| ≤ N、N = max(3 × 群内最大差, 1e-6)): **PASS**。
  問題 `case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml`、m6_on、2D node SST、段階起動 + 本段 12000。A = run_1002–1004 (旧)、B = run_1005–1007 (新)、各 3 反復。本段末尾 6000 step 平均:
  C_T |Δ| 1.5e-6 (N 1.3e-5)、C_T_with_shear 2.0e-6 (1.2e-5)、C_L 8.6e-6 (4.1e-5)、C_M 1.8e-4 (9.4e-4)。6 本とも GATES PASS。

## 期待値と出典
- (0) の 1e-15 は私が事前に置いた閾値で、根拠は「同じ入力の合成なら丸め程度」。計算順序の違いは想定していなかった。

## 仮説
- H1: 差は合成の計算順序 (Python と C++ で和の順序・演算の組み方が違う) による丸めで、物理的な差ではない。(1) の回帰がノイズ床内であることとも整合。

## 問い
1. (0) を「FAIL (1.02e-15、丸め) として記録し、(1) の PASS をもって段 (i) を合格とする」でよいか。それとも (0) の閾値の置き方を改めた上で、別の確認 (例: 係数ごとの ulp 差の上限、あるいは Python 側と C++ 側の合成を同じ入力で直接比べる単体試験) を足すべきか。
2. 段 (ii) (輸送物性の切替: viscMethod 1 → 2 + gas.transport、変化量を記録・合否なし) に進んでよいか。記録すべき量や注意点は。
