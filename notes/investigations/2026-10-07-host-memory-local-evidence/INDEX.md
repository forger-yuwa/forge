# ホストメモリ削減 R1–R3 のローカル証跡 (2026-10-07)

plan [`plans/accepted/architecture-solver-host-memory.md`](../../../plans/accepted/architecture-solver-host-memory.md) §5.1 #4 の「傾き 2737 → 2227 → 1371 B/節点」と負例 2 つの原本
(2026-10-07 result レビュー M2 で、plan から辿れるように scratchpad から移した)。実行はローカル WSL (RTX、sm_86)、`FORGE_MEMLOG=1`、`FORGE_ALLOW_UNVERIFIED_SPECIES=1`、10 step、`outStepInterval 5`。

## 入力
- 縮小した生産 SERN 3D 格子 (生産の問題 YAML `case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml` を `runner_sern3d --prepare-only --op m6_on` で ni/nj/nz を 0.7・0.85 倍、計測専用に第一層厚を 1/s 倍、品質 PASS)。
  s070 = 672,871 節点、s085 = 1,188,542 節点。設定は生産の本段と同じ (node・SST・frozen_tp 2 成分・SLAU・陰解法 block-DPLUR・2 次・level 1)。
  **一様 IC からの開始なので値の回帰には使わない** (「発散前までのメモリ観測」、s050 は step 9 で `roe` 非有限)。
- 負例と checkpoint の確認: `case/09.Taylor-Green/run_0160_passiveG_fct_ckpt100` の入力を複製 (42,471 節点、dual-time・level 2)。

## バイナリ (sha256 の先頭 8 桁、`binaries_sha256.txt`)
| 名前 | コミット | sha | 使った run |
| --- | --- | --- | --- |
| base | 6ca9db49 (R1–R3 の前、計測ログ入り) | 958a65e5 | `runs/mem_s070_base`・`mem_s085_base` |
| R1+R2 | 作業ツリー (c1ea79ca と同内容) | f0bc9585 | `runs/mem_s070_new`・`mem_s085_new` |
| R3 (計測時) | 作業ツリー (93e55957 の 1 つ前の状態) | 86ad3c7e | `runs/mem_s070_r3`・`mem_s085_r3` |
| R3 最終 | 93e55957 と同内容 | a61d1b09 | s070 で 1 回 (step 0・HWM が 86ad3c7e と同じことを確認、run は scratch のみ) |
| 負例 P / roN | R3 + H から 1 名前を外した一時ビルド | b90efd8b / 7f0c8bcf | `runs/ck_neg_P`・`mem_s070_neg_P`・`ck_neg_roN` |

**注意**: s085 の R3 計測は最終バイナリ (a61d1b09) でなく 86ad3c7e。最終版との差は空白・コメントの修正だけとの申告 (実装者)、s070 で step 0 と HWM が同じことは確認済み。

## 傾きの出し方
各工程の値 M (MB) から、傾き = (M_s085 − M_s070) / (N_s085 − N_s070) [B/節点、1 MB = 2^20 B]、切片 = M_s070 − 傾き × N_s070。表は `memlog_cmp_r1r2.txt` (base 対 R1+R2) と `memlog_cmp_r3.txt` (base 対 R3)、
生成は `cmp_memlog.py`、各 run の `[memlog]` 行は `runs/<run>/memlog.txt`、`/usr/bin/time -v` は `runs/<run>/time_v.txt`。
- ホスト HWM のピーク (`after setStructuralVariables`): base 2737 → R1+R2 2227 → R3 1371 B/節点。GPU の増分は base と同じ (s085 で 1872 MB)。

## 負例 (`runs/*/key_lines.txt`・`rc.txt`)
- H から `P` を外す: `ck_neg_P` (case/09、新規開始) と `mem_s070_neg_P` (s070) とも、初期出力の D2H で
  `[variables] ERROR: host cell variable 'P' has length 0 but nCells_all = …` を出して rc 1、`res_0.h5` は作られない。
- H から `roN` を外す: `ck_neg_roN` (case/09 の res_100 から再開) で、checkpoint 復元の書込みより前に
  `[variables] ERROR: host cell variable 'roN' has length 0 …` を出して rc 1 (「history restored」行は出ていない)。
