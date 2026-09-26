# codex 諮問 (diagnose): stage-manifest-launch-binding-design

- **brief**: [`notes/reviews/briefs/2026-09-27-stage-manifest-launch-binding-design.md`](../../notes/reviews/briefs/2026-09-27-stage-manifest-launch-binding-design.md)
- **plan**: [`plans/active/tooling-stage-manifest-launch-binding.md`](../../plans/active/tooling-stage-manifest-launch-binding.md)
- **date**: 2026-09-27
- **commit**: `9607db76` (feature/sern-design)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: **呼び出し側で §4.2 を「起動試行 ID と残差履歴の直接対応」に改め、その欠落耐性を T4 の必須条件としてから §4・§6 を確定する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**現案のままの確定は不可**。Critical 0／Major 4／Minor 1）

| 対象 | 採否・重大度 | 根拠と修正案 |
|---|---|---|
| §4.2・T4：設定ハッシュ＋単調対応 | **却下・Major** | [plan:39](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:39) のアルゴリズムでは、同一設定の起動 `L1,L2,L3` から `L2` が欠落すると、対応は **`s1→L1, s2→L3, s3→不明`** になる。メモリ上で再現した。**起動試行ごとに一意な `launch_id` を発行し、起動記録と、その試行が生成した残差履歴を manifest で直接結ぶ方式へ変更する。** ハッシュは照合用にする。現在の [StageManifest.add:259](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:259) は `bcond_fnv` も保存していないため、これも追加が必要。 |
| §4.4・T3：キー欠如を `gg (legacy)` とする | **却下・Major** | 親 [plan:114](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:114) 自体が、既存 manifest の保存済みキーは修復されないと認めている。**記録していなかったことは、gg で実行した証拠ではない。** 段に対応する旧バイナリ・設定原本で裏付けられる場合だけ `legacy gg` とし、裏付け不能は `unknown`。不明段をまたぐ自動連結は停止する。 |
| §4.5・T6：DB に最後の起動値を付ける | **却下・Major** | [runner_sern.py:41](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:41) は最後に解析できた値を返すだけで、評価した出力との対応を確認しない。最後が失敗した再起動なら、前の成果物に別の実効値が付く設計になる。**DB 行は評価対象成果物の `launch_id`・`exe_sha256`・実効値を参照する。** T6 に「成果物生成後の別設定の失敗起動で行の由来が変わらない」を追加する。 |
| T1・T3・T4・T5：検証の検出力 | **要再検証・Major** | [T1:63](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:63) と T5 は、**バイナリ SHA だけで分割しても合格**でき、実効値の誤対応を見逃す。既存 [chi テスト:40](/home/sano/work/forge-sern-design/solver_density_cuda/tools/test_stage_manifest_wall_normal_chi.py:40) も「同一 YAML」に異なる `cfg_fnv` を与えている。区間数に加えて、**各段の対応 ID・実効値・由来を完全一致で検査**する。T4 は欠落・余分な失敗起動・同一段の再試行を含め、誤対応 0 件を合格条件にする。 |
| §4.1：バイナリ識別 | **採用、記録形式は修正・Minor** | [plan:38](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:38) の正本を実行ファイル SHA にする判断は妥当。**常に SHA-256 全64桁を正本、git rev＋dirty は補助情報**とする。16桁への切り詰めは不要。起動ごとに1回読む設計を採用するが、コスト許容は未測定なので対象環境で読み取り時間を記録する。取得失敗を共通の空値として連結してはいけない。 |
| §4.3：異なるバイナリは別区間 | **採用** | [plan:40](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:40) の保守側の分割でよい。ただし「再ビルドしたか」ではなく、**実行ファイルの内容 SHA が違うか**で判断する。実効値が同じでも SHA が違えば分割する試験を入れる。 |
| T1：既定だけ変えた試験ビルド／T2 | **採用** | 親 [plan §4.4:67](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:67) の既定化前検証として適切。試験ビルドは生産用と分離し、差分・SHA・起動記録の実効値を保存する。T2 は維持し、T3 には **cell で明示 lsq でも実効 gg** を含める（[仕様:371](/home/sano/work/forge-sern-design/procedures/solver-settings.md:371)）。 |

結論: **呼び出し側で §4.2 を「起動試行 ID と残差履歴の直接対応」に改め、その欠落耐性を T4 の必須条件としてから §4・§6 を確定する。**

第 1 仮説: 最大の欠陥は、設定の同一性を起動の同一性として扱っていることにある。バイナリ ID を追加しても、別の起動を段へ割り当てれば防げない。 **確度: 高**  
  根拠: [stage_manifest.py:152](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:152) と [:289](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:289) を実行し、同じ `cfg_fnv` の起動実効値 `[0,1]` が段では `[1,1]`、**1区間**になることを再現した。提案の単調対応でも、中間記録欠落時に上述の誤対応が発生した。  
  反証条件: 記録欠落・余分な起動・再試行があっても、各残差履歴の生成起動を独立した証拠で一意に特定でき、誤対応が0件になること。現案にはその証拠を記録する仕組みがない。

第 2・第 3 仮説: 追加しない。今回の判断には第1仮説で十分。

判別 A/B: **同一設定の3段・3起動を持つメモリ上の fixture で、中間起動記録 `L2` の有無だけを変える。forge は0 step。** 見る量は各段の対応 ID と判定不能状態。欠落時に **`s1→L1, s2→不明, s3→L3`** を保持し、不明段をまたぐ自動連結を拒否できる結果 A なら、欠落による対応ずれの仮説をその実装について棄却する。`s2→L3` へずれる結果 B なら、対応方式を却下する。**現案を忠実に再現したアルゴリズムは B だった。**

やらない方がよいこと: バイナリ SHA の追加と「2区間になった」だけで M7 を完了扱いにすること。再試行を履歴と結び付けず読み飛ばすこと。1つの履歴が複数起動を含む場合も、起動別の範囲を記録せず単一の実効値を付けてはいけない。

呼び出し側の前提への異議: 「途中で記録が欠けた段だけ inferred にできる」は現案からは導けない。また、**記録の根拠**（`launch`／`inferred` 等）と**設定値の決まり方**（`default`／`explicit`）は別項目にするべき。[main.cpp:115](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:115) の実効値と由来を両方保持し、`launch` で後者を潰さない。

不足情報: 許可された範囲には、manifest を確定する呼び出し箇所、再試行時の履歴の保存方法、DB の評価成果物・選別処理が含まれない。このため具体的な接続箇所と T6 の現状は未確認。既存テストは scalarGradient **5件 PASS**、chi **17件 PASS**だが、上記欠陥を検出しない。流れの収束・定常性は今回判定していない。**ファイル変更なし、plan 未反映。**
