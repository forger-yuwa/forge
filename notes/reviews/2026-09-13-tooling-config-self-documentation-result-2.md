# codex レビュー: tooling-config-self-documentation (result)

- **plan**: [`plans/active/tooling-config-self-documentation.md`](../../plans/active/tooling-config-self-documentation.md)
- **stage**: `result` (diff base `HEAD~2`)
- **date**: 2026-09-13
- **commit**: `483de03b` (feature/perf-3d-speedup)
- **codex**: effort `high`, 3.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

抽出 159 キー・未解析 0 件・自己試験 12 項目 PASS、既存 config 32 本の主要指摘件数、起動ログの 53／76／0 行は再現できました。
ただし、`annotate` に値の破壊と誤注釈が残っています。現状のまま `accepted` へ移すべきではありません。

1. **Major — `annotate` がブロックスカラーの値を変える。安全性検査も通過する。**

   根拠: [config_doc.py:413](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:413)、[config_doc.py:442](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:442)。

   次の入力をメモリ上で注釈生成すると、`drive:` の直前にインデント付きの `# [未知] …` が挿入されました。

   ```yaml
   mesh:
     meshFileName: |-
       drive: mesh.h5
   ```

   この位置の `#` は YAML コメントではなく**文字列の一部**です。生成前後の YAML 値は不一致なのに、挿入行を取り除く独自検査は PASS します。「原文を復元できる」ことは「生成物の値が同じ」ことを保証しません。

   **対案:** YAML の構文木・トークン位置から安全なコメント挿入位置を求め、スカラー内部への挿入を禁止する。原文保持検査に加え、スカラーの表記を保持した構文比較を行い、この再現例を試験に追加する。

2. **Major — 実際の検証 config に誤った `[節違い]` を付け、主要設定の注釈を落とす。**

   根拠: [config_doc.py:414](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:414)、[config_doc.py:426](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:426)。実入力は [run_0456 の solverConfig.yaml:28](/home/sano/work/forge-perf/case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/solverConfig.yaml:28)。

   `case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/` の注釈生成を再現すると、正しい `time.deltaT.dt_min` に `[節違い] 正しい節は time.deltaT` が付付きました。複数行のフロー形式を入れ子として認識できないためです。また、`last`・`space`・`turbulence` のフロー形式内のキーには説明も既定値も付きません。`time: # コメント` でも入れ子を取り違えます。

   **対案:** 1 と同じ YAML 構文解析からフルパスを取得する。フロー形式には、その行の前にフルパス付きの注釈をまとめて付ける。実 config について「値が保持される」だけでなく、注釈対象キーと警告の期待値を検証する。

3. **Minor — `detectNaN` の最終値をログで確認できるという説明が誤り。**

   根拠: [solver-settings.md:54](/home/sano/work/forge-perf/procedures/solver-settings.md:54)、[overview.md:217](/home/sano/work/forge-perf/methods/architecture/overview.md:217)。実装の [solverConfig.cpp:382](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:382) は上書き時に出力しません。

   保存済みハーネスへ `time.deltaT.detectNaN: 0` とトップレベル `detectNaN: 1` を渡すと、終了コードは 0、ログには前者の `0` だけが残りました。コード上の最終値は `1` です。

   **対案:** 文書を「この上書きの最終値は現ログでは確認できない」に訂正する。実効値ログの追加は、既存の §5.1 #4 に残す。

4. **Minor — plan の設計方針と検証条件が、対応済みとする実装に追随していない。**

   根拠: [plan §4:67](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:67) は現在も「YAML シリアライザで書く」、[§6:102](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:102) は「読み戻して構造・値が完全一致」としています。実装はコメント挿入と挿入行除去の比較です。

   **対案:** §4・§6 を修正後の実方式と試験に合わせ、今回の未解決事項を §5.1 に登録する。既存の残作業 5 件は記載されていますが、別途対応するものは移動前に後継の active plan への参照を残してください。

指定 diff の C++ 変更はログ追加とコメントに限定され、既定値や数値計算経路の変更はありません。`g++ -std=c++17 -fsyntax-only` は通過しました。`template` は 26 トップレベル項目・禁止キーを除く 158 キーを保持し、入力自身への出力は exit 2 で拒否しました。本 plan は流れ場の収束・定常性を主張していないため、それらの VERDICT を根拠とする数値評価は対象外です。

**推奨:** 原文を保持し、YAML の構文位置に基づいて注釈する方式へ修正してください。移動前の優先順は **1 → 2 → 3 → 4**。追加試験と文書同期を終えてから `accepted` へ移すことを推奨します。今回は read-only レビューのため、ファイル変更なし・**plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 2
