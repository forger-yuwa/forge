# codex 諮問 (diagnose): g4-initial-field

- **brief**: [`notes/reviews/briefs/2026-10-11-g4-initial-field.md`](../../notes/reviews/briefs/2026-10-11-g4-initial-field.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-10
- **commit**: `925b441a` (feature/sern-design)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: **案 (a) を「現行 3D の領域別初期場から段階起動」へ訂正し、まず実格子上の内外ラベルを使う 0-step A/B で初期化の対応付けを検証する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **Major** | **(a) は初期場と起動手順を訂正して採用** | 現行 3D は MOC 初期場ではない。[`runner_sern3d.py:108`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:108) は節点 index で排気・外気を分ける**領域別一様初期場**で、[`:204`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:204) がこれを無条件に呼ぶ。[`:416`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:416) は暖機 ramp も渡さない。`run_0970` の台帳記載版 `2fa3826c` でも同じ経路だった。**既存 3D の初期化と段階起動を使い、その実効設定を記録する**。R-b は [2D・m6_on 限定](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-startup.md:14)で、今回の起動保証には使えない。 |
| **Major** | **(b) の重心への単純置換は却下** | [`interp_field.py:46`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/interp_field.py:46) は、高 AR の双対重心が節点から外れ、近壁で壁値・内部値を交互に拾った既往を理由に節点座標を優先している。**重心が違うことは、正しい側の最近傍になる保証ではない**。今回は補間器の変更を持ち込まない。将来対応するなら、接続・境界面から識別した側で探索候補を制限する。 |
| **Major** | **(c) は却下** | [`interp_field.py:212`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/interp_field.py:212) は座標だけで検索して全転送量へ同じ index を適用する。[旧事故の記録](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:134)は134 station の双子の合流と step 7 の発散。**検査は時間積分前に行う**。事後に回復しても初期場の正しさは証明できない。 |
| **Major** | **「両側 STEADY なら履歴の違いを除外できる」は却下** | [plan §6 #4](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:115) の定常性条件は時間変動を制限するだけで、同じ解への到達を検証していない。**§6 #4 に初期履歴の差を明記し、比較を「異なる起動履歴を含む生産候補の頑健性確認」と位置付ける**。差を純粋な格子効果へ帰属しない。 |
| **Minor** | **g4 の 0 格子は流れを回さない方針を採用** | [既決事項](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:144)どおり、同じ g4 仕様の投入条件の基準として使う。**1.0 格子の品質・双対閉性・`--admission` を確認する**。この比較から「g4 でも 0→1.0 が冷点を消した」とは主張しない。 |

結論: **案 (a) を「現行 3D の領域別初期場から段階起動」へ訂正し、まず実格子上の内外ラベルを使う 0-step A/B で初期化の対応付けを検証する。**

第 1 仮説: **座標だけの cross-mesh 最近傍転送は、カウル・側壁の双子で内外の状態を混同する。** 確度: **高**
  
根拠: [`mesh_sern3d.py:284`](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:284) はカウルの `dup1` に加えて側壁の `dup2` も作り、[`:326`](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:326) は側壁コピーへ同一座標を与える。今回、同じ `cKDTree.query` を3断面・双子3組の合成入力で実行し、**3/3組が同じ donor に合流、6節点中3節点が誤側ラベル**になった。これは合成試験であり、実 g4 の誤写像数は未測定。

反証条件: 実 g3→g4 の全対象双子について、接続から独立に定めた内外ラベルを現行転送が誤り0件で保持すること。実格子に対象の双子がない場合も、今回への適用仮説は成立しない。

第 2 仮説: **起動履歴が異なる準定常状態への到達に影響する。** 確度: **低・未確認**。`STEADY` だけでは未除外。

判別 A/B: **変えるのは初期値の割当方式だけ。時間積分は0 step。**

- g3 のカウル・側壁について、接続から内外を識別した合成ラベルを作る。自由端で共有する単一節点は双子検査から除く。
- **A**：現行 `interp_field` の座標最近傍写像で g4 へ転送。
- **B**：`paste_region_ic3d` の index 分類で g4 に同じ内外ラベルを直接割り当てる。
- g4 の境界面接続を正解として、側ごとの誤割当数・双子の同一 donor 化数を数える。

→ **Aに誤り、Bが0件なら**、転送の曖昧さを確認でき、訂正した (a) を進める。  
→ **Bにも誤りがあれば**、「既存 3D 初期化なら内外を保持する」を棄却し、時間積分へ進まない。  
→ **両方0件なら**、実格子での誤転送仮説を退ける。合成試験の結果を実格子へ読み替えない。

やらない方がよいこと: `mesh.ic: moc` や ramp のキーを書くだけで 3D に効いたと扱うこと、重心の最近傍なら安全と判断すること、壊れた初期場を回して自然回復に期待すること。

呼び出し側の前提への異議: **観測として確認できたのは双子を識別しない探索実装であり、今回の実 g4 の誤写像・発散はまだ観測されていない。** また `run_1079` は、保存判定が **`NOT CONVERGED (stalled/plateau)`**。力CSVを `check_quasisteady.py --tail 0.5` で再判定すると、step 500～20000 の末尾半分で **`OVERALL: ALL STEADY`** だった。「収束解」とは呼ばない。

不足情報: `run_0970/0972` の実初期場・段別設定・VERDICT、今回の g4 入力格子と投入条件判定、plan §5.1 #2 (4) の g3 継続場に対する EOS 比較の合格記録。**後者は既存の時間積分開始条件として維持する。**

**plan 未反映。** 呼び出し側で `tooling-sern-te-wake-grid.md` §5 手順3・§5.1 #4・§6 #4 に訂正と解釈上の制限を反映すること。ファイル変更・forge 起動は行っていない。
