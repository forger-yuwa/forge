# V4 実経路の形式の試験: `radial_node_count` の axis_gap_frac のときの nj (plan tooling-review-fallback-backend §6 V4)

これは `codex_review.py --backend claude` の実経路の形式 (記録の見出し・分類・§6.1 の行) を確かめるための小さな諮問である。
答えは呼び出し側が既に知っている (形式の試験であって、判断を仰ぐものではない)。

## 観測事実

- `design/forge_design/meshing/mesh2d.py` の `radial_node_count(prm, max_nj)` は、`generate_axisym_mesh` が作る径方向の点数 (実効の nj) を返す関数である。
- この作業ツリーでは同じファイルを別の作業が編集中なので、行番号は自分で確かめること。

## 問い

`prm.axis_gap_frac` を指定し (`axis_cap_frac` は指定しない)、`max_nj` を渡さないとき、`radial_node_count` が返す nj は
入力の `prm.nj` そのものか、`_capped_layout` から導いた値か。根拠の行 (`ファイル:行`) を示せ。

## 期待値と出典

- 出典: `design/forge_design/meshing/mesh2d.py` の `radial_node_count` 本体。

## 再現条件・実施済みの操作

- 無し (コードを読むだけで答えが出る)。

## 仮説

- 無し。
