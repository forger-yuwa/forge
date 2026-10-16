# #6e の試作: 終端 C⁺ の補間 (2026-10-16、M6ノズルの設計)

plan `tooling-nozzle-card-local-mesh` §4.5・§6 V9。forge-7d の掃引 (`../lcscan_forge7d/`) と同じ問題・同じ形。

- コード: feature/card-local-mesh e580a1fb の `design/forge_design` と `solver_density_cuda/data` を `git archive` で展開し、`prototype.diff` を当てたもの
  (環境変数 `INTERP=1` で補間、`INTERP=0` で今の最寄り)。Python は手元の `design/.venv-opt`、BLAS 1 スレッド、8 並列。
- 問題: forge-7d の `problem_B.yaml` (sha256 675339395af9…、カード B、r_t 0.10490242069948211 m、M_K 2.7)。
- `INTERP=0` で L_c 37.2020873623039 → x_F 76.24832719272557 (forge-7d の基準とビット一致)。
- `interp2400.json`: `INTERP=1`・n_axis_inv 2400・L_c 37.17〜37.26 の 181 点 (`scan_interp.py`)。比べ方 `cmp.py ../lcscan_forge7d/fine2400.json interp2400.json`:

| | 局所の傾き | 直線の当てはめからの残差の最大 |
| --- | --- | --- |
| 今 (最寄り、forge-7d `fine2400.json`) | 中央値 0.2920 (段で 105.7) | 2.84e−2 r_t |
| 補間 | 0.9999〜1.0001 | 1.34e−6 r_t |
