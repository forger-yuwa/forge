# codex 諮問 (diagnose): sern3d-r5h-te-coldspot

- **brief**: [`../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-te-coldspot.md`](../../../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-te-coldspot.md)
- **plan**: [`../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md`](../../../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md)
- **date**: 2026-09-27
- **commit**: `4f6cf102` (feature/gap-heating-precision)
- **codex**: effort `high`, 2.9 min, rc=0
- **結論**: **R5h は現形状で診断を続け、同一保存場から対流再構成だけを変える A/B で局所の質量・エネルギー収支を確認する。有限厚後縁への移管と監視項目への格下げは、現時点では採用しない。**
- **extra**: `../forge-sern-design/plans/active/tooling-sern-mesh-blocking.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **R5h は現形状で診断を続け、同一保存場から対流再構成だけを変える A/B で局所の質量・エネルギー収支を確認する。有限厚後縁への移管と監視項目への格下げは、現時点では採用しない。**

採否表:

| 重大度 | 主張・判断 | 採否と根拠・対案 |
|---|---|---|
| **Major** | 厚さ 0 の後縁が原因なので、有限厚化に委ねる | **要再検証**。後縁で板厚が 0 になることと、上下が同一ノードを共有することはコードで確認できる。ただし、これは正しく閉じた鋭い後縁であり、離散化の破綻を証明しない。根拠: [mesh_sern3d.py:180](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:180)、[同:239](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:239)。**形状と数値処理の相互作用として診断する。** |
| **Major** | 等温壁化で改善したので、回復温度が寄与した | **却下〔因果の帰属〕**。比較は断熱 CPG → 等温 TP で、熱物性も同時に変わり、最小点の位置・壁距離も違う。根拠: [対象 plan:1533](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1533)、[同:1555](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1555)。**「条件一式の変更で最低温度が上昇した」に訂正する。** |
| **Major** | 力係数の格子差が許容内なので、低温点の影響も小さい | **却下〔上限評価として〕**。g3・g4 とも低温点を含み、共通の誤差は差分から見えない。記録上も `check_convergence: NOT CONVERGED`、力係数のみ `check_quasisteady: STEADY`。根拠: [対象 plan:1857](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1857)、[同:1878](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1878)。**格子差・準定常性・低温点による誤差を分けて扱う。** |

第 1 仮説: **後縁近傍の高次対流再構成が、質量輸送に対して内部エネルギーの供給を不足させ、低温極小を維持している。** 確度: **低〔機序は未確認〕**。  
　根拠: ブリーフ記載の `case/46.sern_design/run_0971_3d_g3_chidef_cont16k/` は 105.4 K、`run_0972_3d_g4_chidef_cont40k/` は 94.7 K で、どちらも後縁直下に局在する。ただし、局在は再構成原因の証明ではなく、優先して調べる根拠に留まる。  
　反証条件: 共通保存場で高次再構成を外しても、低温を維持する局所収支が変わらず、その収支を一次流束・粘性項・境界処理で説明できること。

第 2 仮説: **高次再構成に依存しない、後縁の壁処理・一次流束・粘性輸送の結合に問題がある。** 確度: 低、未確認。有限厚化で消えても、この可能性は否定されない。  
第 3 仮説: **保存量からの温度反転、組成処理、または出力の対応付けによる見かけの低温。** 確度: 低、未確認。保存量と EOS 前後の照合が不足している。

判別 A/B: **同一 g3 メッシュ・同一保存場で、`space.convMethod` だけを現行値から `0` に変える。**

- 現行の実効値が `1` であることを確認し、A＝現行、B＝`0`。値の意味は [solver-settings.md:9](/home/sano/work/forge/procedures/solver-settings.md:9) に従う。既に `0` なら第 1 仮説はこの時点で棄却する。
- `run_0971` の同一スナップショットから、`restart_field.py` で保存量をビット一致させた新規 run を二つ作る。BC・ガス・SST・χ・勾配法・リミッタ・CFL は固定する。
- **各 100 step を診断用に回し、最初の共通状態での残差評価を主判定にする。** 対象は低温点、接する後縁壁点、その隣接 CV。質量・運動量・エネルギー・化学種の面別輸送、粘性寄与、更新量、EOS 前後の保存量と T を記録する。内部エネルギーの評価では運動エネルギー、組成、設定に応じた k の寄与も整合させる。
- **分岐 A:** 再構成差による流束収支の変化が、低温点を加熱側へ動かす更新差を説明し、EOS 書き換えでは説明されない → 再構成依存を支持し、「EOS 処理だけが原因」を除外する。
- **分岐 B:** 再構成差が更新差を説明せず、共通項または EOS 処理が低温を説明する → 第 1 仮説を棄却する。収支が閉じなければ判定不能とする。

この短い A/B は**機序の切り分け専用**であり、温度の定常値や力係数への影響上限は判定しない。一次化で最低温度が上がったという結果だけでも、原因確定にはしない。

やらない方がよいこと: **温度床を上げる、後縁だけ粗くして受理する、旧メッシャと新ブロックメッシャの力係数差を低温点の影響上限と呼ぶこと。** 新方式は側壁厚・フィレット・後流接続も変えるため交絡する。有限厚 plan 自体も、根治確認を B6 に、力の帰属を後続作業に残している。根拠: [mesh-blocking plan:5](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:5)、[同:114](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:114)。

呼び出し側の前提への異議: **「形状起因」と「数値起因」は排他的ではない。** また、外部流温度・粗格子の最低温度・全場百分位は、対象節点の物理的な温度下限ではない。現行の等温壁では、上下の流体温度差だけから「壁温条件が競合している」とも言えない。

不足情報: 対象 `run_0971`〜`run_0973` は両 checkout に存在せず、上記数値・VERDICT はブリーフと plan の記録であり、今回の独立再測定ではない。実効 config、低温点と近傍の保存量・組成・時系列、局所収支が必要。**ファイル変更・forge 起動はしていない。plan 未反映であり、反映先は呼び出し側の `tooling-nozzle-sern-3d.md` §5.1 R5h。**
