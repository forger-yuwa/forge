# codex レビュー: convection-zero-thickness-edge-reconstruction (plan)

- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **stage**: `plan`
- **date**: 2026-10-08
- **commit**: `cd5fa5d4` (feature/sern-design)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
幾何で限定して速度再構成だけを抑える案は、判別実験として妥当です。  
ただし、対象集合・診断計測・設計チェーンへの配線・恒久採用の判定基準を補ってから実装してください。

`plans/README.md` と関連する `accepted/` を確認しました。既存の `slauWallNormalChi` は質量流束の補正であり、本件の速度再構成とは別です。解決済み案件の重複とは判断しません。指定の挿入位置は速度二乗・面エンタルピー計算前で、内部面の逆符号加算も維持できます。block-DPLUR はセル状態による近似ヤコビアンを使うため、今回の変更に合わせて LHS を変更する必然性もありません。[SLAU:279](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:279)、[同:661](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:661)、[時間積分仕様:276](/home/sano/work/forge-sern-design/methods/time_integration/implementation.md:276)。

なお、対象 `run_1055`・診断 run の生 HDF5、残差・準定常 VERDICT はローカルにありません。実測根拠として確認できたのは持ち帰り帳簿で、収束・定常性の再認定はしていません。

1. **Major — E の抽出仕様が、実際の node メッシュ表現と 2D に対して不足している。**

   **根拠:** §4 は「上下壁が共有する辺」、§5 は「変換済み h5」からの生成を指定しています。しかし node 変換後の境界面の `iNodes` は **1 節点だけ**です。元の境界面接続は別の `vizBfaceNodes` に保存されています。また、2D の板端は共有する**頂点**なので、共有辺という定義をそのまま使うと検出できません。[変換器:2339](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:2339)、[同:2686](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:2686)、[plan:39](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:39)。

   **対案:** §4 に、3D は primal 境界面の共有辺、2D は境界線の共有端点から E を作ると明記してください。入力は `vizBfaceNodes` 等に固定し、欠落・タグ誤記を「正常な空集合」と混同しないこと。実格子で **517160・517199 ∈ S₂** を必須条件にし、共有端・座標一致別 ID・有限厚・端なしの小規模検査を実装の最初に置くべきです。

2. **Major — 周期境界を含む場合、現在の S₂ 定義では同じ DOF に異なる処置を掛け得る。**

   **根拠:** 適用条件は node＋SLAU のみで、周期の除外がありません。しかし forge の周期節点は別 ID・別接続面を持ちながら同一 DOF として更新されます。通常の内部面グラフだけで距離を測ると、周期対応点のマスクや継ぎ目を越える距離が一致する保証がありません。[plan:40](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:40)、[periodicNode_d.cuh:14](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cuh:14)。

   **対案:** 初版は目的に必要な**非周期構成に限定し、有効化＋周期境界を起動時エラー**にしてください。cell・非 SLAU も明示的に拒否する契約を追加してください。軸対称で有効な処置まで保証するなら、2D の端定義に加えて軸ピンとの交差を検査する必要があります。未検証構成を空マスク回帰で受入れ済みにしないことです。

3. **Major — ソルバの設定追加だけでは、判定区間と設計評価の新旧分離が成立しない。**

   **根拠:** §5 の変更対象に `stage_manifest.py`、runner、設計 DB の扱いがありません。読み取り専用で現行 `stage_key()` を実行した結果は、次のとおりでした。

   ```text
   stage_key OFF == ON: True
   stage_key rings2 == rings3: True
   ```

   現行 runner は固定テンプレートから `space` を生成し、評価来歴には `slauWallNormalChi`・`scalarGradient` を記録しています。新設定の配線と来歴追加がなければ、手動 A/B が成功しても設計ループで処置が抜けたり、新旧評価を混在させたりできます。[stage_manifest.py:266](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:266)、[runner_sern.py:346](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:346)、[runner_sern3d.py:359](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:359)。

   **対案:** §5 に「問題 YAML → 各段の solver config → 起動時実効値 → manifest／metrics／設計 DB」の配線を追加してください。有効状態・正規化したタグ・リング数を区間識別に含め、マスクとメッシュの識別情報も保存すること。**OFF→ON とリング変更で区間が分離する検査、旧評価の持ち越し禁止が働く検査**を必須にしてください。

4. **Major — 要求する帳簿と毎更新の床監視を、既存計測だけでは満たせない。**

   **根拠:** §6.1 は EOS・境界補正を別々に要求しますが、既存帳簿は `entry` と **`after_eos_bc`** の間に壁・軸・等温ピン、EOS、BC をまとめています。面帳簿も上限 200000 件で、超過時は切り詰めます。既存 `floor_gate()` は最終保存場を調べるため、§6.2 の「末尾区間の毎更新で補填なし」の証明にはなりません。[main.cpp:1975](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1975)、[同:2018](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2018)、[帳簿:665](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:665)、[同:721](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:721)、[sern_gates.py:200](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:200)。

   **対案:** A/B より前に計測実装を独立ステップとして追加してください。状態更新・ピン・EOS・BC の差分を分離し、面帳簿の欠落は判定不能にすること。長期 run には全域の床違反件数・床由来の保存量補正・発生位置を毎更新で集計する軽量カウンタが必要です。TP の EOS は通常時にも `roe` を再構成するので、**単なる `roe` の丸め差を床補填と数えない定義**も必要です。[dependentVariables_d.cu:208](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:208)。

5. **Major — 恒久採用の一部が「測る項目」に留まり、合否を事前に決められない。**

   **根拠:** §6.2 には「冷点移動なし」「温度場・収支・格子・リング感度」「設計差を変動幅込みで判定」とありますが、冷点の定義、比較領域、細分系列、許容差、保留条件がありません。一方、完了条件は「全項目が合格」です。これでは結果を見て合格基準を選べます。[plan:76](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:76)、[同:80](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:80)、[同:98](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:98)。

   **対案:** §6.2 を「量・領域・窓・許容差・判定不能条件」の表にしてください。特に固定リング数の細分では処置の物理幅も縮むため、**格子感度と処置幅感度を分ける系列**を事前指定すべきです。力の `d` は前後窓平均差の絶対値と明記し、差の差 δ の採否基準も固定してください。温度精度の認定を保留する方針は妥当ですが、設計用途として受け入れる感度の上限は必要です。細分格子の品質 VERDICT と cross-mesh restart も、この検証ステップに組み込んでください。

6. **Minor — §1 の「機序は診断済み」は m10_on まで含めると過剰な断定。**

   **根拠:** 持ち帰り帳簿では、`run_1066_r5h_m10_B_ledger100` の 517199 は call 10 で実内部エネルギー更新 **−7.755e3 J/kg**、EOS 差 **＋4.188e−4 J/kg**。これは当該更新の冷却が EOS 上書きで説明されない根拠ですが、未介入面の速度再構成が主因であることまでは分離していません。[帳簿:16](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R5H_M10_LEDGER100.txt:16)、[plan:20](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:20)。

   **対案:** 「m6_on の局所残差維持機序は診断済み。m10_on の移動先に対する S₂ 処置の十分性は未確認」と修正してください。§6.1 の反証可能な構成と揃います。

**推奨は、幾何限定・速度のみ・opt-in の案を維持し、上記 1→2→3→4→5 の順で §4・§5・§6 を具体化してから実装することです。** 空マスク回帰と短期 A/B の成功だけでは恒久採用せず、計測可能になった受入れ条件を通して判断してください。

ファイルは変更していません。上記修正案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
