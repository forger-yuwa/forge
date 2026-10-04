# codex レビュー: boundary-node-farfield-characteristic (result)

- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `result` (diff base `bfdfd040^`)
- **date**: 2026-10-03
- **commit**: `c111796d` (feature/sern-design)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2
- **focus**: §2 の 2026-10-03 ユーザ決定で受理範囲を縮小した (汎用の遠方境界としては受理せず、V0–V2 の合格項目と SERN 3D 側方の用途に限る)。この範囲での accepted の可否を判定してほしい。diff は farfield 関係 (solver_density_cuda/cuda_forge/convection/farfieldFlux_d.inc.cuh, convectiveFlux_d.cu, boundaryCond.cpp, ransBoundary_d.cu ほか §7 の影響範囲, tests/unit/test_farfield_flux.cu, design/forge_design/evaluate/runner_sern3d.py, design/forge_design/metrics/sern_momentum.py, methods/boundary.md) に絞ってよい。species-transport ブランチのマージ (種 DB・輸送・凝縮・リミッタ調査) はこの plan の範囲外
- **extra**: `case/58.farfield_verification/README.md`, `case/46.sern_design/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

§2 の縮小された受理範囲は妥当です。V2c・V2d-2 の未達だけを理由に却下しません。
ただし、入力検査・使用手順・検証証拠の不足を解消するまでは、`accepted` への移動を認めません。

1. **Major — 多成分 `farfield` の必須組成が、省略しても受理される。**

   [boundaryCond.cpp:175](/home/sano/work/forge-sern-design/solver_density_cuda/boundaryCond.cpp:175) は `farfield` にも入口用の既定補完を適用し、組成未指定を `Y0=1`、残りを 0 にします。[plan §4.1:52](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:52) の「多成分は Y/X 必須」と不一致です。SERN の EXH/AIR 配列では、外気組成の書き忘れが純 EXH の自由流として通り、EOS と流入流束が変わります。

   **対案:** `farfield` では組成解決結果が空なら起動時に拒否する。既存入口の補完は維持し、組成省略の拒否試験と、明示した Y/X の正常系を追加する。

2. **Major — 文書どおりの設定では、側方境界が `slip` のままになる。**

   [methods/boundary.md:167](/home/sano/work/forge-sern-design/methods/boundary.md:167) は `mesh3d.side_far_kind: farfield` と案内していますが、[runner_sern3d.py:62](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:62) が読むのは `evaluate.side_far_kind` です。関数を抽出した実行確認でも、前者では `kind: slip`、後者では `kind: farfield` が生成されました。今回受理する用途そのものに影響します。

   **対案:** 文書を `evaluate.side_far_kind` に訂正し、生成 YAML の実効 BC を検査する試験を残す。併せて、`top_out_kind: outflow` は現在 `ValueError` になるため、`outlet` と `outlet_kind: outflow` を使う仕様を明記する。

3. **Major — result レビューに必要な一次証拠を、このレビュー環境で監査できない。**

   [case/58 README:5](/home/sano/work/forge-sern-design/case/58.farfield_verification/README.md:5) のとおり、対象 run は AWS にあり、手元には対象 run の残差・VERDICT 原本がありません。AWS 接続も失敗しました。

   ローカルの [V3_EVAL.txt](/home/sano/work/forge-sern-design/case/46.sern_design/V3_EVAL.txt) と [V3IC_EVAL.txt](/home/sano/work/forge-sern-design/case/46.sern_design/V3IC_EVAL.txt) は、plan の旧輸送幅比較・履歴比較の数値と整合します。しかし、これだけでは V0–V2 の合格、新輸送の #4f/#4h、全区間の NaN・置換ゼロを独立に確認できません。**結果が誤りだという指摘ではなく、受理の裏付けが未監査という指摘です。**

   **対案:** 受理対象の run ごとに、実効設定・バイナリ識別子・残差 CSV・収束／準定常 VERDICT・専用試験の判定原本・置換カウンタ記録を回収し、判定区間付きで索引化する。`GATES PASS` と残差収束は別欄にする。V0 は初回流束だけでなく、§6 が要求する更新後保存量と旧版反復幅の比較記録も提示する。

4. **Minor — 診断ダンプの実装が、完了扱いの §4.3 を満たしていない。**

   [farfieldFlux_d.inc.cuh:225](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/farfieldFlux_d.inc.cuh:225) の出力は流れの5流束と外側の `Y0/k/omega` で、化学種・乱流の実際の面流束はありません。[convectiveFlux_d.cu:697](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:697) の帳簿も種 0/1 の固定列で、転送は節点・変数ごとです。「実種数から列挙・一括転送」と一致しません。

   **対案:** 限定受理では「2種試験で必要量を再構成した診断」と正確に記述し、一般種数の流束出力・一括転送を未完了項目として残す。

5. **Minor — 縮小後の受理範囲と最新結果が、文書間で同期されていない。**

   [methods/boundary.md:176](/home/sano/work/forge-sern-design/methods/boundary.md:176) は独立参照の精度を未達としていますが、plan §5.1 #3e は完了です。[case/46 README:367](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:367) の「#4f の結論を持ち越し」は、plan の「広幅への持ち越しは未証明」という訂正を反映していません。[plans/README.md:28](/home/sano/work/forge-sern-design/plans/README.md:28) も `draft・plan-9 レビュー待ち` のままです。

   **対案:** `methods`・case 台帳・計画一覧・§8 完了条件を§2の限定受理に統一する。V2c・V2d-2 と広幅の LJ 応答は、後継 plan の具体的な残作業へリンクする。2水準の格子差を精度保証と読める表現も避ける。

**推奨は、限定受理を維持したまま、①一次証拠の回収・照合、②組成検査と設定手順の修正、③診断仕様・残作業・文書の同期、の順に完了してから移すことです。** 専用流束への分岐、一次／S3 の流入面値、混合後 `k_R`、運動量ダンプの圧力基準補正には、確認した範囲で設計との重大な相違を認めませんでした。

ファイルは変更していません。以上の提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
