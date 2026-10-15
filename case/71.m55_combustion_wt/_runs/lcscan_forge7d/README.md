# from_length の x_F(L_c) の掃引 (カード B、2026-10-16、forge-7d)

- コード: feature/card-local-mesh 2871a489 の `design/` と `solver_density_cuda/data` を `git archive` で `src/` に展開。Python: `/home/sano/work/forge/design/.venv-opt` (numba 0.68.0)、BLAS 1 スレッド、手元 WSL。
- 問題: M6 の試走のカード B の出力 `problem_p0_card.yaml` の写し `problem_B.yaml` (sha256 675339395af9…) = r_t 0.10490242069948211 m・M_K 2.7・x_end_margin 2.3・axis_dx0 0.03・analytic/converge。
  L_c を固定 (`Lc_mode: explicit`) して `design_chain` を回し、`f = x_F − L_total` (L_total 76.26134789508718 r_t) を記録 (`scan.py`)。
  再現の確認: L_c 37.2020873623039・N 2400 で x_F 76.24832719272557 (試走の記録とビット一致)。
- `fine2400.json` (37.17〜37.26、0.0005 刻み 181 点)、`wide2400.json` (37.0〜37.5、0.01 刻み 51 点)、`fine4800.json` (0.001 刻み 91 点)、
  `align2400.json` (診断の A/B: x_E に最も近い軸の節点を x_E に乗せた網、`ALIGN_XE=1`、91 点)。解析 `analyze.py`、節点との照合 `grid_check.py`。
