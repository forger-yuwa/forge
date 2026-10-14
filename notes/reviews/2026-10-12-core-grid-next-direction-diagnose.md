# codex 諮問 (diagnose): core-grid-next-direction

- **brief**: [`notes/reviews/briefs/2026-10-12-core-grid-next-direction.md`](../../notes/reviews/briefs/2026-10-12-core-grid-next-direction.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-12
- **commit**: `6cd2c912` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.2 min, rc=0
- **結論**: **(a) を採用し、Gref の上限 50 の双子を同じ起点から各 3,000 step だけ再実行して、250 step ごとの列別振幅から波束自身の増幅と移動を測る。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **(a) を採用し、Gref の上限 50 の双子を同じ起点から各 3,000 step だけ再実行して、250 step ごとの列別振幅から波束自身の増幅と移動を測る。**

元の問いへの現在の答えは、「上限 50 に依存する微小差の増幅と、縮流部の粗化による停滞の緩和は観測されている」です。「細分すると増幅域に長く滞在するため停滞する」は有力な候補ですが、**その因果も、反復が差を増幅する数値的な原因も未確定**です。

| 重大度 | 採否 | 指摘・根拠・対案 |
|---|---|---|
| Major | 採用 | **機構の調査を続ける。** Gref の追跡データ不足が明記されており、元の依頼は未完了。まず同じ測定方法で比較可能にする。[plan:1040](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1040) |
| Major | 要再検証 | **「停滞の正体は非線形に飽和した揺れ」は断定が早い。** Grefc は一度頭打ちした後、再上昇している。双子の差は定常解からの誤差でもない。「有限時間の差の増幅」と「残差停滞を維持する原因」を分ける。[plan:1029](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1029)、[plan:1037](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1037) |
| Major | 却下 | **「増大の約 1/4 は流入」という寄与率の認定。** 追従窓と固定窓の対数傾きの差だけでは、形の変化や窓内の分布変化と流入を分離できない。「追従窓の傾きは固定窓の 0.73 倍」と記載する。正式判定も「判定不能」のまま維持する。[plan:1046](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1046) |
| Major | 却下 | **「収束した参照」という呼称。** Gref point の `check_convergence` は通算 0〜200k で **NOT CONVERGED**。160k〜200k の対象 8 量は **STEADY** なので、「対象量の参照資格を満たした場」とする。[plan:663](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:663) |

第 1 仮説: **Gref でも、固定領域への流入だけでなく、移動する波束自身の有限時間増幅が起きている。細分化による滞在反復数の増加は、その先に検証すべき格子差の説明候補である。** 確度: 中。

- 根拠: Gref の固定領域では λ = 3.83e−3/step、上限 1 では λ = −4e−5/step。一方、Grefc の追従窓では λ = 3.30e−3/step、R² = 0.995。ただし、Gref の追従窓の値は未測定。[plan:1008](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1008)、[plan:1046](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1046)
- 反証条件: 十分な信号と追跡資格を保った Gref の追従窓で、増幅がほぼ消えること。その場合、現在の λ を波束の滞在時間に掛ける説明は退ける。

第 2 仮説: Gref の固定領域の増大は、波束の流入・形の変化による寄与が大きく、Grefc と同程度の波束自身の増幅を表していない。

第 3 仮説: 格子間の差には、別々に発達した背景場や局所の反復作用の違いが強く寄与し、滞在反復数だけでは説明できない。未確認であり、今回の試験だけでは除外できない。

判別 A/B: **初期の ρE 摂動の有無だけを変える、Gref の基準・摂動の 2 本。**

- 起点は既存の `case/45.isobutane_m6_d155/run_0619_cg_gref_b0_ctl/` の 60,000 step。§4.62 と同じバイナリ・設定・摂動配列を使い、上限 50 のまま各 3,000 step。今回の追加計算はこの 2 本に限る。
- 250 step ごとに列別振幅を保存し、**step 1500〜2750 の 6 点**を主判定区間として事前固定する。固定窓・追従窓は同じ振幅定義とし、§4.67 の線形範囲・捕捉率・窓幅感度・相互相関の資格を使う。信号が測定限界を十分上回ることも確認する。
- **結果 A:** 追従窓で λ ≥ 1e−3/step、R² ≥ 0.95、前半・後半とも正 → Gref の波束自身の増幅を支持し、「流入だけで増大した」を退ける。
- **結果 B:** 固定窓の増大が再現する一方、追従窓の λ が固定窓の 20 % 以下 → 第 1 仮説を退け、滞在時間による説明を組み直す。
- 中間・資格不足は判定不能。これは**新しい事前判定**であり、§4.68 の旧判定は変更しない。

同じ出力で、追跡した中心の実際の列移動から速度を求め、既存 Grefc と比較する。物理座標速度を代表 Δx で割る方法は使わない。波束上の実効 `dt_local` と律速条件も記録し、「Δτ ∝ Δx」を確認事項にする。**結果 A でも、背景場の交絡、総増幅の列数依存、残差・θ_r への因果は未解決**と明記する。

やらない方がよいこと: この段階で上限 1 の生産既定化、中間上限の探索、左辺の別項の変更へ進むこと。Gref の上限 1 は平均流を改善した一方、20k の判定は **NOT CONVERGED（RISING）** で、連成系の解決は示していない。[plan:708](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:708) 熱伝導の近傍結合だけでは停滞が緩和しなかった試験も、今は繰り返さない。

呼び出し側の前提への異議: 「対策が見つかったので原因調査を閉じられる」という整理は受け入れない。また、約 0.04 列/step、頭打ちの約 30 倍差、固定領域の λ を組み合わせても、同じ波束の総増幅を測ったことにはならない。

不足情報: Gref の線形区間の列別振幅系列、波束上の実効 dt と律速条件、波束増幅が残差停滞を維持することの直接的な証拠。今回のレビューは指定範囲の文書記録に基づき、生データ・判定ツールの再実行は行っていない。**ファイル変更なし・plan 未反映。**
