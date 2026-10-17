# §6 V1a: 実行役の起動の入力と V6 の run の比較、旧の台本の検査と実行役の関門の対応表 (2026-10-17)

plan [`tooling-design-orchestrator`](../../../../plans/active/tooling-design-orchestrator.md) §5.1 #6d・§6 の「V1a の事前登録」。GPU・forge なし (準備だけ、手元)。

## 比較 (事前登録の 9 項目)

- 実行役の側: recipe の NS の段と同じ CLI (`deltastar_loop --init-integral --ic-file run_0008_ns_b0/res_90000.h5 --ic-sha256 c64a73f1… --euler-ref run_0009_euler
  --euler-ref-res res_12000.h5 --euler-ref-sha256 2b2d02ce…`、問題 `problem_m55_p5_settle_v2.yaml`) を `--prepare-only` で一時の run に準備し、`campaign.plan.confirm` で計画を確定した。
  `run_0008_ns_b0` の最後の場と設定は AWS B から写した (sha256 は B の上と同じ)。
- 比べる相手: V6 `run_0014_chain_v2` (台本 `../run_0014_chain_v2/p6_chain.sh`、AWS B、git_head f20cded9、FP64 のビルド)。
- 台本 [`compare.py`](compare.py)、結果 [`result.json`](result.json)。

| 項目 | 結果 |
| --- | --- |
| (1) 全段の起動直前の設定の FNV (S1_soft・S2_mid・main) | 一致 |
| (2) bcondConfig.yaml・(3) solverConfig_prestage.yaml | バイト一致 |
| (4) settle_criteria.yaml | **FAIL** (差は `source.generator_sha256` (許した) と `source.written` (書いた時刻; 事前登録の正規化に書き漏らした)) |
| (5) settle_nodes.txt・(6) 格子・(7) 準備の記録の初期壁の数値 | **FAIL** (x_E 40.633257 vs 40.630797、L_c の差 0.00246 r_t; 下の切り分けを参照) |
| (8) 初期場の固定 | 一致 (ic_snapshot の sha256 c64a73f1…・step 90000) |
| (9) バイナリ・変換器の sha256 の記録 | 記録あり (値の一致は求めない; 手元と B の FP64 のビルドは別) |

**V1a は事前登録の文言どおり FAIL**。訂正: (6) はトポロジが一致し、違うのは座標と `nozzle.msh` (codex の再計算; 最初の報告で「トポロジも違う」と書いたのは誤り)。節点の表は実行役 3,830 行・V6 3,831 行。`compare.py` は部分文字列で除外していて緩い (計測時間・配置のパスを除き、停止規則を事前登録より広く比べた; codex 諮問 `notes/reviews/2026-10-17-orchestrator-v1a-fail-diagnose.md`)。記録として残し、次の V1a′ では使わない。

## V1a′: 同じ今のコードで、旧の台本の呼び方 (A) と実行役の CLI (B) の比較 (事前登録 2026-10-17)

- 因子は呼び出し経路だけ (commit ea5b4f6a、手元の card4a のビルド、同じ問題 YAML・Euler の参照・初期場)。A = `run_pass0_integral(..., ic_from=<run_0008 の写し>, stop_rule=True, prepare_only=True)` (`p6_chain.sh` と同じ)、B = `deltastar_loop --init-integral --ic-file … --ic-sha256 … --prepare-only`。A が選んだ最大番号の場は res_90000.h5 (sha256 c64a73f1…)。
- 台本 [`compare_ab.py`](compare_ab.py) (完全なキーパスで 1,706 キー: 全段の設定の全文・境界条件・前の段の設定・solverConfig・停止の基準・節点の表・格子の座標とトポロジ・**補間した後の初期場 (VALUE の全データセット)**・msh・prepare_info の全キー)。
- 揺らぎの欄 (B を 2 回準備した差、A と比べる前に記録; [`v1a_prime_noise.json`](v1a_prime_noise.json)): `prepare_info:physical_wall.fit.build_seconds` (計測時間)・`prepare_info:settle.sha256` (書いた時刻を含む基準のファイルの sha256)・`settle_criteria:source.written` (書いた時刻) の 3 つ。
- 結果 ([`v1a_prime_ab.json`](v1a_prime_ab.json)): **PASS** — 差 6 = 揺らぎ 3 + 別に記すもの 3 (`ic_snapshot.how` 最大番号 / 固定したファイル、`stop_rule.explicit` true / false、`stop_rule.requested_source` 引数 / 既定)、残りの差 0。この入力について、呼び出し経路の差 (ディレクトリの最大番号の場 / 固定したファイル、停止規則の明示 / 既定) は起動の入力に効かない。
- 解釈の限定 (事前登録): A = B なので V1a の V6 との差は呼び出し経路によらない。ただし、その差の唯一の原因が #6e (eabbab97) とは確定しない (V6 とは環境・変換器も違う)。A を「V6 の再現」とは呼ばない。

## 旧の台本の検査と実行役の関門の対応表

旧の台本: `../run_0011_ns_v2/p5_ns.py`・`../run_0011_ns_v2/p5_all.sh`・`../run_0014_chain_v2/p6_chain.sh`。

| 旧の台本の検査 | 出典 | 実行役の関門 | 試験 |
| --- | --- | --- | --- |
| run が既にあれば止める | p5_all.sh・p6_chain.sh (`test -e`)、p5_ns.py:104 | driver: 出力先が既にあれば起動の前に止める | `run_driver_tests.py` (T12) |
| FORGE_BIN の sha256 の固定 (準備と起動で同じ) | p5_all.sh (FBSHA)、p5_ns.py:106・184 | 走らせる環境の固定 (driver が初回に台帳へ写し、再開で照合)、計画の `inputs.forge_bin`、走った後の照合 (`RUN_PROVENANCE.txt` の forge_sha256) | `run_driver_tests.py` (T13)、`run_execution_plan_tests.py` (N5) |
| Euler の参照の sha256 の固定 | p5_all.sh・p6_chain.sh・p5_ns.py:108〜109 | recipe の NS の段が Euler の段の結果の sha256 を `--euler-ref-sha256` で渡し、`deltastar_loop` が照合 | `run_recipe_tests.py` (R9)、本物は §6 V1b |
| 初期場の場がある | p5_all.sh (`test -f res_90000.h5`) | 初期場をファイルと sha256 で固定 (`select_ic`) | `run_campaign_ic_tests.py` (N1〜N3) |
| settle_criteria.yaml がある | p5_ns.py:114 | 計画は停止規則 v2 の run だけ確定する (`staged_ns_launch_configs`) | `run_execution_plan_tests.py` (P1) |
| 初期壁に較正の係数がある | p5_ns.py:117 | recipe の NS の段: カードの `deltastar_initializer.*` が準備の記録の initializer にある (`initializer_mismatch`、28945473) | `run_recipe_tests.py` (I1) |
| PRESTAGE_CONFIG がある・本段に B0 のキー・前の段に B0 のキーが無い | p5_ns.py:119・132〜135 | runner の `ns_stage_configs` の検査と、計画の全段の設定の全文の固定 | `run_execution_plan_tests.py` (P2: FNV が V6 と一致) |
| 停止規則を明示で使える | p5_ns.py:142 | `_stop_rule_v2_for_main` が使えなければ例外 (計画を確定しない・起動しない) | `run_execution_plan_tests.py` (P1 の前提) |
| 本段の nStepOuter が指定どおり | p5_ns.py:152 | 計画の本段の step 数 = 停止規則の上限、起動の直前の照合 | `run_execution_plan_tests.py` (P1・B) |
| 環境変数の要求 (FORGE_CUDA_BLOCKSIZE ほか) | p5_ns.py:183 | 計画の `inputs.env` と走った後の照合 (起動の記録の env)、走らせる環境の固定 | `run_execution_plan_tests.py` (N5)、`run_driver_tests.py` (T13) |
| 準備の後に設定・格子・前の段の設定・停止の基準が変わっていない | p5_ns.py:187〜189・201 | 起動の直前に全段の設定と境界条件を計画と照合、最初の段の前に停止の基準・節点の表・格子と初期場・準備の記録の sha256 を照合 (28945473) | `run_execution_plan_tests.py` (B・N4・N6) |
| 倍精度の implicitSolvePrecision・thermoFloat の条件 | p5_ns.py:61〜85・193〜197 | 対象外 (recipe v1 は float の生産; 精度の切り替えを持たない) | — |
| 別の run (run_0012) との設定の差が nStepOuter だけ | p5_all.sh | 対象外 (2 本の組の比較; 計画が全段の設定の全文を固定する) | — |
| V6 の確認: 停止規則 REACHED・関門 proceed・settled・非有限 0 | p6_chain.sh | recipe の NS の段の合格の印 (停止規則 REACHED・関門 proceed・SETTLED・非有限 0・計画の照合 OK・初期壁の係数) | `run_recipe_tests.py` (分岐)、本物は §6 V1b |
| V6 の確認: delta_r_next.csv が停止の抽出とバイト同一 | p6_chain.sh | 対象外 (recipe v1 の次の壁は較正の係数で作り、delta_r_next.csv を使わない) | — |
| 成果物 (`--skip cad`) | p6_chain.sh | recipe の成果物の段 (判定専用の環境、CAD も作る) | `run_recipe_gate_tests.py` (A) |

未対応 0 (対象外 3 は理由つき)。表を作る中で見つけた穴 2 つ (準備の後のファイルの変更・初期壁の係数) は 28945473 で関門にし、試験を足した。
