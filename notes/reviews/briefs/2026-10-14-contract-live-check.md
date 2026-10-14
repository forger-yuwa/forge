# 諮問ブリーフ: 出力の契約の実機確認 (plan tooling-review-fallback-backend §5.1 #5・#7)

ツールの出力の契約 (判定の行・`証拠:` の欄) を実物の判定役が守るかを確かめるための小さな問い。

## 観測事実
- `design/forge_design/meshing/mesh2d.py` の `radial_node_count` は、`axis_gap_frac` のとき指定の nj から実効の nj を数える。

## 問い
`radial_node_count` で、`max_nj` を渡し、かつ指定の nj − 1 が `max_nj` を超えるとき、`_capped_layout` (配列を作る関数) は呼ばれるか。根拠の行を示して、推奨を 1 つ (「このままでよい」も可)。

## 禁止事項 (厳守)
- ファイルを変更しない。
- 下の範囲以外のファイル読みをしない。grep は可。
- 推奨は 1 つに絞る。根拠は `ファイル:行` で示す。

## 読んでよいもの
- `sed -n '150,190p' design/forge_design/meshing/mesh2d.py`
