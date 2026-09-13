# codex レビュー: tooling-config-self-documentation (result)

- **plan**: [`plans/active/tooling-config-self-documentation.md`](../../plans/active/tooling-config-self-documentation.md)
- **stage**: `result` (diff base `HEAD~1`)
- **date**: 2026-09-13
- **commit**: `c7e5e96b` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**  
数値計算経路・既定値の変更はありませんが、注釈生成による値の改変と、設定点検の見逃しを再現しました。  
以下を修正・再検証してから `accepted` に移してください。ファイルは変更していません。

1. **Major — `annotate` が化学種名を改変し、往復検証も通してしまう**

   根拠: [config_doc.py:253](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:253) の `yaml.safe_load` は、引用符なしの `NO` を真偽値として読みます。[同ファイル:383](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:383) の比較は、この変換後の値同士です。

   メモリ上で入出力を代替して `cmd_annotate` を実行した結果、`physProp: {species: [N2, NO]}` が **`species: [N2, false]`、exit 0** になりました。ソルバは [solverConfig.cpp:635](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:635) で各要素を `as<std::string>()` として読むため、元の `NO` と生成後の `false` は別の化学種名になります。

   **対案:** YAML の scalar 表記・タグ・構造を保持して注釈を挿入し、ソルバが読む値の保存を検証してください。`NO`／`on`／数値に見える文字列を回帰試験に追加します。

2. **Major — `check` が異なる入れ子構造を同一視し、空の禁止キーも見逃す**

   根拠: [config_doc.py:262](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:262) は節名をドットで連結し、辞書自体は照合しません。

   `case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/solverConfig.yaml` から既知の未知キー2個をメモリ上で除いた入力を基準に、次を再現しました。

   - `time.deltaT` の入れ子を、トップレベルの単一キー `"time.deltaT"` に移動：**exit 0、指摘0件**。
   - `mesh.bndFirstOrder: {}` と未知の空節を追加：**exit 0、指摘0件**。

   前者をソルバは [solverConfig.cpp:340](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:340) の `config["time"]["deltaT"]` として読めません。後者も使用禁止キーの検出要件を満たしません。

   **対案:** パスを文字列のタプルとして保持し、既知の節・設定キーを各階層で照合してください。空辞書も検査対象にします。

3. **Minor — 正当な空セクションを含む config を注釈生成できない**

   根拠: [config_doc.py:347](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:347) は `output: {}` を `output:` と出力し、空辞書を null に変えます。実行結果は **exit 2、出力なし**。破損防止は働きますが、plan §6 の構造保存要件は未達です。

   **対案:** 空辞書を `{}` として保持し、空セクションを含む往復試験を追加してください。

4. **Minor — `keepDissCoeffMax` の説明が適用領域を逆に示す**

   根拠: 新しい [solverConfig.hpp:155](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.hpp:155) は「DES シールド外で使う」と説明しています。しかし [convectiveFlux_keep_d.inc.cuh:169](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh:169) の実装は、`ransFrac` を使った `max(keepDissCoeff, ransFrac*keepDissCoeffMax)` です。上限側が効くのは RANS 側です。

   **対案:** 「RANS 領域で使う散逸係数の上限。`keepDissFdBlend=1` で有効」と訂正してください。

5. **Minor — 検証完了を追跡できる記録が不足している**

   根拠: [plan §6:95](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:95) は検証項目の列挙で、40本の対象一覧・判定内訳・ハーネス実行結果への参照がありません。[変更ログ:145](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:145) にも、それらの実測がありません。

   今回は `coverage` の **159キー／未解析0件／説明なし0件**、禁止キーなしの template、同一ファイルへの出力拒否、C++ 構文検査を再確認しました。ただしハーネスはリンク時の一時ファイル作成が read-only 制限で拒否され、起動ログの実行確認はできていません。

   **対案:** 試験入力・期待値・実行結果・40本の対象一覧を保存し、§6 から参照してください。本レビューの修正事項と採否も §5.1／§6.1 に反映します。

指定の `c7e5e96b` の diff では、境界面・周期 seam・軸・浮動小数点演算・ゼロ割ガードに変更はありません。対象 plan に流れの収束や改善率の主張はなく、CFD の VERDICT 比較は該当しません。既知の古い config の問題は §5.1 に残っており、`methods/index.md` と `plans/README.md` の現状も整合しています。

**推奨は、上記1→5の順に修正・記録し、result 再レビュー後に `accepted` へ移すことです。** 本応答の指摘は plan 未反映です（read-only）。

指摘数: Critical 0 / Major 2 / Minor 3
