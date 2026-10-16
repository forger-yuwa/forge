# 収束の判定の新方式の result 段レビューの証拠の写し (2026-10-16、M6ノズルの設計)

plan `tooling-settle-convergence-and-monitor` §5.1 #13 (result 段レビュー 2 回目の「証拠不足」への対応)。原本は AWS B の `~/m55_p5/` と
`~/forge-cardmesh/case/71.m55_combustion_wt/run_*` (B のディスクにある)。

| ファイル | 中身 |
| --- | --- |
| `v3_replay/`・`v3b_replay/` | §5.1 #8 の再生 (プリセット `wt_nozzle_v3`; `v3b` は節点ごとの壁の熱流を外した後): 再生に使った基準 `criteria_v2_<run>.yaml`、作り方のログ `make_<run>.log`、再生のログ `replay_<run>.log`、run ごとの結果 (`<run>/`: `VERDICT.txt`・`grid_verdicts.csv`・`span_deviation.csv` など) |
| `replay_v3.sh`・`replay_v2.sh` | 再生を回した台本 (B の `~/m55_p5/`) |
| `listing_<run>.tsv` | run の全ファイルの一覧 (大きさ [bytes]・更新時刻・パス)。V4 の出力の数と大きさの根拠 |

run の写し (手元、git に入れない; 全場 `res_<n>.h5`・`nozzle.h5`・`*.msh` は写していない):
- `case/71.m55_combustion_wt/run_0015_ns_float_a/`・`run_0016_ns_float_b/` (float の V8 の run): 判定の記録 (`SETTLE_VERDICT.json`・`CONVERGENCE_VERDICT.txt`・`STOP_RULE_STATE.json`・`STOPPED_AT`・`MAIN_STAGE_STARTED`・`forge_launches.jsonl`・`settle_criteria.yaml`)、節点の出力 `res_nodes_*.h5`、境界の出力、`boundary_massflux.csv`、`residual_history*.csv`
- `run_0011_ns_v2/`・`run_0014_chain_v2/` に、`STOPPED_AT` の `files` にある停止の時の入口・出口・壁・軸の h5・xmf と `res_nodes_<停止の step>.h5` を足した

一覧から数え直した出力 (V4 (vi) の照合):

| run | 全場 `res_<n>.h5` | 監視の出力 (節点・境界の h5 と xmf) | run 全体 |
| --- | --- | --- | --- |
| run_0011_ns_v2 (V4) | 4 枚 (0・20,000・40,000・46,202)、403.3 MB | 216 ファイル、68.9 MB | 768.7 MB |
| run_0014_chain_v2 (V6) | 4 枚 (0・20,000・40,000・46,171)、403.3 MB | 216 ファイル、68.9 MB | 1,024.1 MB (成果物・報告を含む) |
| run_0015_ns_float_a (V8) | 5 枚 (0〜60,000・61,758)、274.0 MB | 279 ファイル、48.5 MB | 618.8 MB |
| run_0016_ns_float_b (V8) | 5 枚 (0〜60,000・60,720)、274.0 MB | 279 ファイル、48.5 MB | 617.5 MB |

V4 の結果に書いた「全場以外の出力 76.4 MB」は、数え方の記録が無く、この一覧からは再現できない (監視の出力 68.9 MB に `forge_run.log` 5.2 MB を足すと 74.1 MB)。
どの数え方でも事前登録の上限 (≤ 300 MB) の合格は変わらない。
