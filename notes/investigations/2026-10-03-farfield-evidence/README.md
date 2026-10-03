# farfield 受理範囲の一次証拠 (2026-10-03、codex result M3)

plan [`boundary-node-farfield-characteristic.md`](../../../plans/active/boundary-node-farfield-characteristic.md) §2 の限定受理の裏付け。run 本体は AWS
(`~/forge-pgrad-new/case/58.farfield_verification/`、`~/forge-pgrad-new/case/46.sern_design/`、`~/forge-r8/case/46.sern_design/`) にあり、ここには
小さなテキストだけを写した (回収スクリプトは AWS で実行、`collect_ff_evidence.py` の出力)。

## 索引 `INDEX.tsv` (147 run)

| 列 | 意味 |
| --- | --- |
| `bin` | `RUN_PROVENANCE.txt` の forge sha256 先頭 12 桁 (無い run は `forge_launches.jsonl` の実行ファイルサイズ) |
| `head` | 同 git_head |
| `conv` | `CONVERGENCE_VERDICT.txt` の判定行 (**残差の収束**。SERN と V1 は全列プラトーで NOT CONVERGED が既知の性質) |
| `gates` | `metrics.json` の SERN ゲート判定 (**残差収束とは別**: NaN・床・rising・量の定常性) |
| `steady` | 力係数の `check_quasisteady` 判定 |
| `ff_hll/vac` | `forge_run.log` の farfield 累積カウンタ (HLL 退避 / 真空・非物理の置換)。表示が無ければ増加なし |
| `nan` | ログ中の NaN 検出 |

受理範囲の run はすべて `ff_hll/vac` = 0/0・`nan` なし。例外は破棄済みの不具合版 `run_0002_v1b`・`run_0003_v1c` (置換 2.6e6 / 1.3e6、GPU の構造体代入の欠陥、修正後の `run_0005_v1b_fix`・`run_0006_v1c_fix` で置換 0)。

## 受理項目と判定原本

| 項目 | 判定原本 (このフォルダ) | run |
| --- | --- | --- |
| V0 初回評価のビット一致 | case/46 README の run_0991 行 (md5 同一) | `run_0991_ff_v0_old/new` |
| V0 後半 (多 step 再現性) | `c46_V0B_REPRO_20261003.txt` (**事前規則で FAIL、D/S 1.08–1.84**; 判断は plan §5.1 #5a) | `run_1032_v0b_{o,n}{1,2,3}` |
| V0u (v) 起動拒否 + 組成必須 | `c58_V0U_REJECT_VERDICT.txt`、`c58_V0U_REJECT_20261003.txt` | 一時ディレクトリ |
| V0u/V0k 単体 | `solver_density_cuda/tests/unit/test_farfield_flux.cu` (リポジトリ) | — |
| V1 自由流保持 | plan §5.1 #3 (VERDICT ファイルなし、eval_v1*.py の出力) | `run_0001`, `run_0004`–`0006` |
| V2a 音響反射 | `c58_V2A_VERDICT.txt` | `run_0040`–`0048` |
| V2b 保存収支 | `c58_V2B_VERDICT.txt` | `run_0030`–`0035` (+`_bal`) |
| V2d-1 接触波 | `c58_V2D_VERDICT.txt` | `run_0070`/`71`・`75`/`76`・`80`/`81` |
| V2e 陽解法・SST | `c58_V2E_VERDICT.txt` | `run_0050`–`0055` |
| V2f 局所逆流 | `c58_V2F_VERDICT.txt` | `run_0036` (+`_bal`) |
| V2c (未達・判定不能) | `c58_V2C_*`、`c58_V2C_ORD_20261003.log`、`c58_V2C_EXT_20261003.log` | `run_0063`–`0107`、`run_0150`–`0166` |
| V2d-2 (未達・測定不能) | `c58_V2D2_*`、`c58_V2_AB_20261003.log`、`c58_V2_REP_20261003.log` | `run_0110`–`0161` |
| 独立参照 | `c58_REF1D_COMPARE_20261003.txt` | AWS `~/forge-r8/case/58.farfield_verification/ref1d/` |
| V3 幅・初期場履歴 (旧輸送) | `c46_V3_EVAL.txt`、`c46_V3IC_EVAL.txt` | `run_0992`–`0998` |
| 窓条件 (3D、単体評価の再計算) | `V3_WINDOWS.txt` | case/46 の 3D run |
| #4c/#4e/#4f/#4h・R8–R10 | `INDEX.tsv` の該当行 + plan §5.1 の数値 | `run_0999`–`1031` |
