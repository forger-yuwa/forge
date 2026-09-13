# codex レビュー: tooling-config-self-documentation (plan)

- **plan**: [`plans/active/tooling-config-self-documentation.md`](../../plans/active/tooling-config-self-documentation.md)
- **stage**: `plan`
- **date**: 2026-09-13
- **commit**: `513f36d9` (feature/perf-3d-speedup)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
課題設定は妥当で、`plans/accepted/` に同機能の既存計画は見当たりません。ただし、現案では誤った設定説明・見逃し・誤検知を生成します。  
既存の実装途中ファイルを読み取り専用で検証しました。ファイル変更・CFD 実行は行っておらず、数値結果の不変性は未検証です。

1. **Major — 抽出方式が現行コードを網羅していない。**

   根拠: [config_doc.py:129](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:129) の束縛抽出は `auto ... = config[...]` に限定され、同ファイル153行の型抽出はネストした `std::vector<T>` を扱えません。実行結果は「157キー、読み出し163箇所、説明なし0」ですが、次が欠落・誤登録されています。

   - `mesh.wallDistExtraPhysIDs`: 抽出されない。実際には [solverConfig.cpp:216](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:216) で読む。
   - `output.extraFields`: 抽出されない。実際には [solverConfig.cpp:431](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:431) で読む。
   - `physProp.chemistry.mechanismFile`: `ch.mechanismFile` として登録される。[solverConfig.cpp:683](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:683) の `const YAML::Node ch = physProp["chemistry"]` が原因。

   `case/16.nozzle_wys/_aws_perf_evidence/ref_run_0234_evidence/solverConfig.yaml` に対する実行でも、正しい `wallDistExtraPhysIDs: [6]` を「ソルバは読まない」と誤診しました。廃止キーも文字列連結で拒否する `LESmodel`／`DESmode` が抽出されません。

   **対案:** 対応構文を明文化し、ネスト型・束縛の連鎖・動的な廃止キー列挙を扱うこと。未対応の読み出しは所在付きで抽出エラーにし、完全な一覧として公開しないこと。

2. **Major — フルパス照合という設計と実装が逆で、必須条件も見逃す。**

   根拠: [config_doc.py:207](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:207) は親節からトップレベルまで検索します。メモリ上の入力試験では、`space.keepDissCoeff: 0.05` を加えても `check` は **exit 0、節違い指摘0件**でした。ソルバが読むのは [solverConfig.cpp:282](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:282) のトップレベルだけです。

   また、`chemistry.enabled: 1` で `mechanismFile` が無くても **exit 0**。ソルバ側には [solverConfig.cpp:692](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:692) の必須検査があります。

   **対案:** パスは完全一致にする。`detectNaN` 等の互換パスは、実在する読み出し位置だけを登録する。必須性には条件を持たせ、無条件の `getValidatedValue` 検出だけで「必須欠落を点検できる」としないこと。

3. **Major — 「説明なし0」は説明の正しさを保証しない。実際に危険な誤説明が入った。**

   根拠: 追加された [solverConfig.hpp:24](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.hpp:24) は `gpu` を「使用するGPU番号、0起点」と説明しています。しかし [variables.cpp:232](/home/sano/work/forge-perf/solver_density_cuda/variables.cpp:232) は `useGPU == 1` の場合にGPUメモリを確保します。`gpu: 0` はGPU番号0の指定ではありません。生成一覧にもこの誤説明がそのまま出ました。

   **対案:** `gpu` をGPU使用フラグとして訂正し、追加する説明は値の消費箇所まで照合すること。コード内の自然言語コメントも陳腐化するため、「コードから生成するので腐らない」という前提は撤回し、意味のレビューを完了条件に追加すること。

4. **Major — 全キー入り雛形が禁止キーを再投入する。**

   根拠: [config_doc.py:335](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:335) の `template` は禁止・廃止情報を参照せず全キーを出します。実行結果には有効なYAML項目として `mesh.bndFirstOrder: (コード側で分岐)` が含まれました。これは関連する [architecture-bndfirstorder-removal.md:20](/home/sano/work/forge-perf/plans/active/architecture-bndfirstorder-removal.md:20) の使用禁止方針と逆行します。

   **対案:** 雛形から禁止・廃止キーを除外すること。未解決の既定値や必須値のプレースホルダーはコメントとして出し、値として解釈される形で出さないこと。

5. **Major — `annotate` に「写し」の保証と元ファイル保護がない。**

   根拠: [config_doc.py:218](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:218) は文字列をエスケープせず引用符で囲みます。メモリ上の往復試験で、引用符を含むパスはYAML構文エラー、バックスラッシュを含む文字列は値が変化しました。

   さらに [config_doc.py:330](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:330) は出力先を無条件に書き込みモードで開くため、`-o` に入力自身を指定すると元configを上書きします。

   **対案:** YAMLシリアライザを使い、生成後に再読込して値・型の同一性を検証すること。入力と出力が同一ファイルの場合は、シンボリックリンク・ハードリンクを含め拒否すること。

6. **Major — §6の検証では、上記の不具合を合格させてしまう。**

   根拠: [plan:74](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:74) は説明件数と既存40本の指摘確認が中心です。しかし説明なし0でも指摘1～5は再現しました。[config_doc.py:133](/home/sano/work/forge-perf/solver_density_cuda/tools/config_doc.py:133) の分母は関数定義や配列要素の変換も数えるため、キー網羅率の基準になりません。また「指摘が正しい」確認だけでは見逃しを測れません。

   **対案:** 読み出し箇所と抽出結果の対応表を作り、未対応0件を検査すること。正しい入力に加え、節移動・条件付き必須欠落・配列・廃止キー・文字列特殊文字の期待結果付き試験を追加すること。

   `g++ -std=c++17 -fsyntax-only -Isolver_density_cuda .../solverConfig.cpp` は成功しましたが、ログ動作や数値不変性の検証にはなりません。C++読込試験でログ抑止・明示値・省略値と解決後設定の同一性を確認し、§6(4)を残すなら変更前後の標準ケース比較方法も明記してください。

7. **Minor — `[default]` は実効設定の完全な記録にならない。**

   根拠: [solverConfig.cpp:201](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:201) の `discretization` 省略経路はヘルパーを通りません。また [solverConfig.cpp:381](/home/sano/work/forge-perf/solver_density_cuda/input/solverConfig.cpp:381) はネスト側の既定値をログした後、トップレベルの `detectNaN` で上書きします。node の `gradLSQ=2` 強制も別経路です。

   **対案:** 読込時の既定候補と最終解決値を区別すること。最終値の記録を追加するまでは、ログの説明を「ヘルパー経由の省略記録」に限定してください。

8. **Minor — 死にキーの除去対象から、再流入元の推奨レシピが漏れている。**

   根拠: [plan:68](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:68) は生産configの掃除を挙げますが、[recommended-settings.md:92](/home/sano/work/forge-perf/procedures/recommended-settings.md:92) 自体が `kInf`／`omegaInf` を推奨しています。実際に読むのは `kInit`／`omegaInit` で、[variables.cpp:642](/home/sano/work/forge-perf/solver_density_cuda/variables.cpp:642) に初期場欠落時の適用条件があります。

   **対案:** 推奨レシピと生成元を先に修正すること。既存runの削除だけで終わらせず、意図した初期乱流値が必要なら正しいキーと適用条件を記載してください。

9. **Minor — plan/resultレビュー統合はリポジトリの免除条件に該当しない。**

   根拠: [plan:83](/home/sano/work/forge-perf/plans/active/tooling-config-self-documentation.md:83) の統合理由は、[AGENTS.md:168](/home/sano/work/forge-perf/AGENTS.md:168) の実装前レビュー要件を満たしません。起票日は免除期限より後です。`check_plans.py` も **`VERDICT: FAIL (0/1 plans OK)`**、理由はレビュー記録パス欠落でした。§5には `methods/` の仕様更新もありません。

   **対案:** 今回をplanレビューとして記録し、resultレビューは別途実施すること。既存実装を試作として明示し、`methods/architecture/overview.md` §7.1、plan §4・§5.1・§6を整えてから実装を継続してください。

**推奨は、コード由来の生成方針を維持しつつ、対応構文を限定して未対応を検出し、C++読込との照合で保証する方式に改訂することです。** 汎用C++解析器の導入より、この規模では費用対効果が高いと判断します。ログ追加そのものに数値式やCUDAカーネルの変更はなく、主なリスクは設定を誤案内することです。

実装継続前の優先順は次のとおりです。

1. 指摘1～3を反映し、抽出・照合・説明の正確性を仕様化する。
2. 指摘4～5の生成物と元ファイル保護の契約を決める。
3. 指摘6の期待結果付き検証を§6へ追加する。
4. ログの保証範囲、推奨レシピ、仕様文書、レビュー記録を同期する。

読み取り専用の依頼に従い、以上の提案は **plan未反映** です。

指摘数: Critical 0 / Major 6 / Minor 3
