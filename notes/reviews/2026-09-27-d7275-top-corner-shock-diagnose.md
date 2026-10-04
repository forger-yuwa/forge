# codex 諮問 (diagnose): d7275-top-corner-shock

- **brief**: [`notes/reviews/briefs/2026-09-27-d7275-top-corner-shock.md`](../../notes/reviews/briefs/2026-09-27-d7275-top-corner-shock.md)
- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **date**: 2026-09-27
- **commit**: `3c4d0609` (feature/gap-heating-precision)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **上端 `slip` を維持し、既存格子の上に8層だけ追加して H≈0.794506 m にする高さ A/B を選ぶ。G15 は未達のまま維持する。**
- **extra**: `case/60.flatplate_d7275_m7/README.md`, `case/60.flatplate_d7275_m7/acceptance.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **上端 `slip` を維持し、既存格子の上に8層だけ追加して H≈0.794506 m にする高さ A/B を選ぶ。G15 は未達のまま維持する。**

第 1 仮説: **前縁からの圧力擾乱と、近すぎる上端 `slip` の離散境界処理との干渉が、上端下流の残差を維持している。** 確度: **中。ただし「共有角ノードに定常解がない」は未確認。**

  根拠:
  
- `case/60.flatplate_d7275_m7/run_0012_t26_wall` の台帳は、流れ残差二乗和の99.7–99.9%が上端下流に局在すると記録している。記載された判定は **`NOT CONVERGED (stalled/plateau)`**、St は **`ALL STEADY`**。これらは原本を再判定した値ではない。[acceptance.json:145](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:145)
- メッシュ定義の H=0.5、160節点、伸長率1.05959255から独立計算すると、**上端直下の節点は y=0.471876666 m、最上層厚さは0.028123334 m**。報告された「衝撃波位置 y≈0.472 m」は、この節点列と一致する。圧力擾乱が上端のごく近くまで来ていることは支持するが、その位置精度は粗い。[fp_d7275_y3_wall.geo:22](/home/sano/work/forge/case/60.flatplate_d7275_m7/mesh/fp_d7275_y3_wall.geo:22)
- `slip` は法線速度を反転する閉包なので、上端を遠ざける操作には実際の作用がある。[boundaryCond_d.cu:68](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68)

  反証条件: **上端と圧力擾乱の間に十分な未擾乱域ができ、移植過渡も落ち着いたのに、旧領域の絶対残差と局所振幅が対照と同水準なら、「上端が近いことが主因」を棄却する。**

第 2 仮説: **衝撃波付近の空間離散化と陰的更新による反復振動。** 確度: 低。緩和 A/B は記録上「判別保留」であり、除外されていない。[README.md:37](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:37)

第 3 仮説: **初期場移植の非忠実性による過渡の混入。** 確度: 低。`run_0012` は同じ幾何・分割なのに `interp_field.py` を使用しており、同ツールは運動量・乱流量を原始量から再構成する。ただし、それが40k後の局在残差を説明する証拠はない。[interp_field.py:157](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:157)

判別 A/B:

- **対照**: `run_0012` を対照の出発点にする。同じ最終保存量から、A=H=0.5 m の新規継続、B=上方8層追加の新規 run を各20,000 step。過去の `run_0012` 末尾だけとの比較では、継続時間と再起動の影響を分離できない。
- **唯一の設計変更**: 上端位置。B は元の x 座標・既存 y 座標・接続を保持し、最上層から同じ伸長率で8層追加する。`.geo` 定義からの高さは **0.794506364 m**。壁終端2.8 m、全BC、物性、SST、空間スキーム、CFL、緩和、内反復数、実バイナリを固定する。生成後に共通領域の座標・接続とメッシュ品質を検査する。
- **初期化**: A は `restart_field.py --keep-src-dtype`。B は cross-mesh 移植後、共通節点の全12保存量がAと一致することを確認する。追加領域の初期値は入口自由流に固定する。旧上端の擾乱を最近傍で新上端まで引き延ばす初期化は避ける。壁距離は各メッシュの値を保持する。
- **測定**: 全残差列、旧上端付近・新上端・出口の成分別残差、局所 `ρ,P,k,ω` の時系列、固定10点のSt、位置IIと固定区間のパネル平均熱流束。局所振幅は連続反復の時系列で測り、離れた2スナップショットの差で代用しない。
- **判定窓**: まず12–16kと16–20kを比較する。流れ・化学種の対象指標の中央値変化が5%を超えるなら過渡として、両腕とも40kまで延長して同じ末尾窓で判定する。それでも落ち着かなければ判別保留。
- **事前の読み**:
  - **Bで旧上端域の各未達流れ・化学種残差がAの1/10以下となり、局所振幅も減衰し、別の場所へ同水準の残差が移っていない** → 上端近接説を支持。「その場所の格子・更新だけで上端位置によらず維持される」説を退ける。
  - **Bで上端から圧力擾乱が離れても、同じ指標がAの0.8–1.2倍に残る** → 上端近接主因説を棄却する。
  - 中間的改善・残差の移動・過渡継続は保留。k/ωの未収束は別途残し、流れ残差の改善で免除しない。

  **1/10は診断上の支持条件に限る。G15の代替合格条件は置かない。** B自身の同一設定区間で `check_convergence` のPASS、対象量の `check_quasisteady`、壁解像のゲートを確認する。未収束の `run_0012` を合格済み参照床として使わない。

やらない方がよいこと:

**Major — 現行生成器に単に `--H 0.8`、または `--y-file --wall-dn` を渡して「高さだけの変更」と扱うこと。** 前者は伸長率を再計算して既存のy分布も変える。後者は `wall_dn` を受け取らず、下流底辺をslipへ戻し、x方向伸長率も1.06／1.03へ固定する。今回の元格子は1.0／1.01である。関数をメモリ上で実行して出力を確認した。**対案は既存格子保持の追加層方式と、生成物の差分検査。** [gen_mesh.py:102](/home/sano/work/forge/case/56.gap_tp1187/gen_mesh.py:102)、[gen_mesh.py:152](/home/sano/work/forge/case/56.gap_tp1187/gen_mesh.py:152)、[gen_mesh.py:163](/home/sano/work/forge/case/56.gap_tp1187/gen_mesh.py:163)

呼び出し側の前提への異議:

- **Major — 「衝撃波が共有角に当たり、定常解を持たない」は観測を越えている。** 報告された残差最大位置はx=2.758 mで、出口角x=2.8 mではない。圧力閾値の最上点も格子列に量子化されている。対案は「上端付近の圧力擾乱との干渉」という仮説に限定し、上記A/Bで判別する。[acceptance.json:146](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:146)
- **Major — 高さ変更後のRMS比だけでは改善を測れない。** 実装は√(Σres²/N)なので、節点追加だけでも値が下がる。今回の160→168列では約0.976倍。対案は全域の√Σres²と、同じ物理領域の体積当たり残差・最大値も併記する。[residualMonitor_d.cu:151](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:151)
- **Major — 「6回とも5e−6以下」は提示資料では立証されていない。** 後半のEXT／WALLの記録は「表示桁で一致」であり、抽出器のR表示は小数3桁。対案は丸め前の点別差を再集計する。確認できても、それは試した変更への感度であり、真値との誤差上限やG15免除の根拠にはならない。[acceptance.json:125](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:125)、[d7275_compare.py:52](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/d7275_compare.py:52)、[plan:1030](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1030)
- **Major — BC変更だけなら、別HDF5ファイルになったことはcross-mesh移植の理由にならない。** 今回確認した元／壁延長 `.geo` の差はPhysical指定の2行だけ。対案は実HDF5の節点対応を確認し、同一なら保存量を直接コピーする。[restart_field.py:63](/home/sano/work/forge/solver_density_cuda/tools/restart_field.py:63)

不足情報: **対象runの実効設定、実バイナリ識別情報、残差履歴、VERDICT原本、残差集計コード、上端付近の時系列がローカルにない。** 数値結果は台帳の引用であり、独立再判定ではない。独立確認したのはコード、生成器の出力、`.geo`の差分と幾何計算。ファイル変更・`forge`起動なし。**plan未反映。**
