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

## V9 (plan §6、事前登録 2026-10-16) の実施 — 実装 (`moc_inverse._terminal_bracket`・`_cplus_wall_cross`) で

入力 `problem_B.yaml` (sha256 675339395af9336536eb63d708ee45eaa90619381be25784840549cbf2df6d35)。

| 項目 | 結果 | 記録 |
| --- | --- | --- |
| (a)・(a′) 節点に乗る x_E・端・異常な入力 | PASS (節点に置くと最寄りの値とビット同一、中点は 2 本の平均、範囲の外・NaN・inf・昇順・重複は ValueError) | `design/tests/run_terminal_bracket_tests.py` |
| (b′) 解像度の A/B (7 点、2400 と 4800) | PASS: |x_F(2400) − x_F(4800)| 最大 1.99e−6 r_t = 0.00021 mm (許容 0.25 mm)、終端の Mach 線との整合 最大 1.1e−9 r_t、切り口の差 0、壁の検査の違反 0 | `ab_res.py`・`ab2400.json`・`ab4800.json` |
| (b) なめらかさ (181 点) | PASS: 2400 は傾き 0.99990〜1.00010・残差 1.34e−6 r_t (試作とビット同一)、4800 は 0.99995〜1.00005・3.67e−7 r_t | `scan_interp.py`・`v9b_2400.json`・`v9b_4800.json` |
| (c)・(d) 4 枚の `wt_card` | PASS: 4 枚とも自動調整なし (n_axis_inv 2400、L_total_tol はカードの公差から)、寸法 PASS、壁の検査で止まらない。拡大部の長さのずれ case/71 +1.2e−7 m・A −1.2e−7 m・B −1.5e−7 m・C −4.5e−8 m。M_knot 2.80・2.65・2.70・2.85 (前と同じ)。各 57〜59 s (手元 12 コア) | `v9c/*.log` |
| (e) 既存の試験 | 下の表 | |
| (f) ほかの分岐の不変 | PASS: streamline・flux の characteristic と cplus × lip は変更の前後で x_F・切る前の壁の sha256 が同一。cplus × characteristic だけ変わる (x_F 19.190001 → 19.192146; flux の経路は 19.192157) | `invariance.py`・`inv_old.json`・`inv_new.json` |
