# codex 諮問 (diagnose): length-match-interp

- **brief**: [`notes/reviews/briefs/2026-10-16-length-match-interp.md`](../../notes/reviews/briefs/2026-10-16-length-match-interp.md)
- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-16
- **commit**: `45d161ea` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.0 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-length-match-interp-diagnose.codex.log` (codex: ok)
- **結論**: 案 b の既定化に向け、まず補間を固定した2400／4800の小規模A/Bで、出口位置の解像度依存を確認する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否は **案 b を採用、既定化は追加検証後、対象経路の #6d 自動倍化・許容差緩和は外す**です。ファイル変更・forge 起動は行っていません。

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 案 b：終端 C⁺ の交点補間 | **採用** | 現行の最寄り選択は [moc_inverse.py:773](/home/sano/work/forge-cardmesh/design/forge_design/geometry/moc_inverse.py:773)。保存された181点を再集計し、補間後の傾き **0.9998967〜1.0001037**、当てはめ残差最大 **1.34×10⁻⁶ r_t**、失敗0を確認した。節点を E に合わせた案 a の成功88点との出口差も最大 **1.37×10⁻⁶ r_t**。網を作り替える案 a より、今回の修正として案 b を推奨する。 |
| V9だけで既定化する | **要再検証・Major** | [plan:192](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-card-local-mesh.md:192) は各解像度での滑らかさを測るが、**2400と4800の出口位置の差を制限していない**。一定の偏りがあっても合格する。解像度間の差と終端 Mach 線との整合を追加する。既存の傾き・残差閾値はカードBの回帰条件として維持してよいが、他の軸則まで傾き1を要求しない。 |
| #6dを保険として残し、r_tごとにリセットする | **却下・Major** | [wt_card.py:1969](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:1969) が毎回のカード由来設定を `auto` で上書きし、[同:1997](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:1997) が許容差0.05を保持する。r_tごとのリセットは持ち越しを直すだけで、そのr_tでも公差より緩い結果を許す。補間の検証後、**`cplus × characteristic` の長さ未収束に対する自動倍化・緩和を外す**。未収束行の検出・記録と停止は残し、谷を黙って捨てない。場の延長 `x_end_margin` は別原因への対処なので維持する。 |
| 試作の節点一致・範囲端処理 | **要再検証・Minor** | [prototype.diff:41](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/lcscan_e6/prototype.diff:41) は無条件に `m_lo=m_hi+1` を参照する。式を再実行すると、降順軸 `[3,2,1]` の `x_E=1` は節点上でも `IndexError`。**節点一致ならその1本を直接返す**。非有限・範囲外・補間に必要な特性線の不足は明示的に拒否する。V9(a)に最上流節点と異常入力を加える。 |
| 「同じL_cなら壁は変わらない」「移動量最大0.019 r_t」 | **訂正・Minor** | 不変なのは切断前の網と `wall_full`。切断後の表から壁を再構築するため、スプラインの端条件も変わる（[wall_axismach.py:182](/home/sano/work/forge-cardmesh/design/forge_design/geometry/wall_axismach.py:182)）。保存181点の旧方式との差は最大 **0.0262356 r_t＝2.75218 mm**。旧値からの変化をカード公差で制限せず、長さを解き直した後の寸法と壁検査で判定する。 |

物理面では、上流側のC⁺が非一様域を通ることだけで補間を退ける理由にはなりません。使用中の `knot` 則は E で M′＝M″＝0 として一様域へ接続しています（[axis_law.py:180](/home/sano/work/forge-cardmesh/design/forge_design/geometry/axis_law.py:180)）。実際、`lcscan_e6/interp2400.json` の全181点で、一様流の終端 Mach 線との整合誤差  
`|x_F − x_E − r_F√(M_d²−1)|` は最大 **1.18×10⁻⁹ r_t** でした。ただし、これでMOC全体の二次精度やCFD出口の一様性まで証明したことにはしません。

`_cut_wall_at_exit` は新しいx_Fで **全列を同じ壁区間から補間**するため、方式として整合しています（[moc_inverse.py:737](/home/sano/work/forge-cardmesh/design/forge_design/geometry/moc_inverse.py:737)）。V9には既存の `run_moc_exit_cut_tests.py` を明示し、`0 ≤ x_F − wall_inv[-1,0] < 1e−9 r_t`、半径の一致、再構築後の `wall.validate()` を含めてください。

既定化の変更範囲は `cplus × characteristic` に限定し、キーや試作の `INTERP` 環境変数は製品コードに残さない方針でよいです。`streamline`／`flux` の `terminal_exit` と `cplus × lip` は別分岐です（[moc_inverse.py:1030](/home/sano/work/forge-cardmesh/design/forge_design/geometry/moc_inverse.py:1030)、[同:1072](/home/sano/work/forge-cardmesh/design/forge_design/geometry/moc_inverse.py:1072)）。これらの不変性も回帰対象にします。

結論: **案 b の既定化に向け、まず補間を固定した2400／4800の小規模A/Bで、出口位置の解像度依存を確認する。**

第 1 仮説: 補間後は2400で出口位置の誤差をカード公差内に抑えられ、自動倍化・許容差緩和は不要になる。確度: **中**。最寄り選択がのこぎり歯の主因である確度は高い。  
根拠: `case/71.m55_combustion_wt/_runs/lcscan_e6/interp2400.json` の181点、上記の傾き・終端線整合、および `lcscan_forge7d/align2400.json` との88点比較。  
反証条件: 下記A/Bで出口差が事前の許容を超える、または元の公差・2400のままではV9(c)の4枚を通せない。

第 2 仮説: 補間で段差を消しても、MOCの網の誤差や隣接C⁺の不足が実用上の制約として残る。確度: **低・未確認**。

判別 A/B: **変更は `n_axis_inv` のみ、A＝2400、B＝4800**。同一のカードB入力・r_t・M_knot・補間方式を固定し、L_c＝37.170、37.1885、37.2070、37.2075、37.2255、37.2440、37.260の7点を各1評価する。forgeは0 step、4800は逐次実行とする。見る量はx_F、終端線整合、切断端整合、壁検査。追加の判別基準として、**全点で出口差≤カード長さ公差の1/4＝0.25 mm**、終端線整合≤1e−4 r_t、切断端整合と壁検査の違反0を事前登録する。  
→ **A：すべて以内**なら第1仮説を支持し、V9全体を完了して既定化・対象経路の自動調整撤去へ進む。  
→ **B：いずれか超過**なら「2400で十分」を棄却し、第2仮説の切り分けへ戻す。許容差を緩めて通さない。

やらない方がよいこと: 傾きが1になっただけで物理的精度まで合格とすること、#6dの持ち越しだけ直して保険を恒久化すること、既存試験の期待値を理由なく更新すること。

呼び出し側の前提への異議: O(Δx²)は現時点では見込み。「網と壁はx_Eに依らない」も不正確で、`x_end` 自体がx_Eに依存する（[runner_axismach.py:586](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:586)）。正しくは「同じ設計入力で出口選択方式だけ変えた場合、切断前の網と壁を保存する」。

不足情報: 指定ディレクトリには `problem_B.yaml` 本体がなく、補間4800と4枚の通し結果も未保存。次の検証時に入力本体・完全なhash・結果を一緒に残すこと。これは上記の次手の判断を妨げないが、既定化の完了証拠には必要。

証拠: 足りる
