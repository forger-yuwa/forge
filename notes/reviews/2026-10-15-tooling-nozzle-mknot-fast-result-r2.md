# codex レビュー: tooling-nozzle-mknot-fast (result)

- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-15
- **commit**: `24f5a17c` (feature/mknot-fast)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.0 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-mknot-fast-result-r2.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

最終版のカード2枚は再照合で `VERDICT: PASS`。選択値・寸法・所要 81.7 s／72.1 s は保存記録と一致しました。
今回の変更範囲で Critical／Major の不具合は確認していません。ただし、保証の説明・試験・現在仕様の3点を直してから移してください。

1. **Minor — 「真の境目まで 0.05 以上」は、この探索では保証できません。**  
   根拠：[plan:76](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:76)、[現在仕様:407](/home/sano/work/forge-mknot/methods/design/overview.md:407)。実装は格子点と両隣の成立を確認しますが、格子点間の不成立域は検査しません。  
   合成反例を実装に通すと、不成立域を `2.015 < M_K < 2.035` に置いても、選択 `2.00`、最終3点 `1.95／2.00／2.05` は全部成立、推定距離は `0.50` でした。実際の最寄り境界は `0.015` です。これは実ノズルでの発生報告ではなく、保証文への反例です。  
   **対案：** 現行アルゴリズムを維持し、「標本上の成立区間から推定した距離」と記述してください。受け入れた限界に、2N 固有の不成立だけでなく、格子点間の狭い不成立域も含めます。

2. **Minor — 点堆積の「対の全集合」の比較になっていません。**  
   根拠：[試験:556](/home/sano/work/forge-mknot/design/tests/run_moc_jit_tests.py:556) は参照側の対 `(a,b)` を `(段, x_a, r_a)` に縮約し、[試験:586](/home/sano/work/forge-mknot/design/tests/run_moc_jit_tests.py:586) で比較しています。[実装:421](/home/sano/work/forge-mknot/design/forge_design/geometry/moc_fill_nb.py:421) も片端の位置だけを返します。同じ `a` に対して相手 `b` を取り違えても、件数が同じなら通ります。§4.6・§6.4 の「対の全集合が一致」という主張より弱い試験です。  
   **対案：** 試験用の取得経路で `(段, 元の添字a, 元の添字b)` を返し、切り詰め前の集合を比較してください。公開出力の先頭20件は、現在の順序契約のままで構いません。

3. **Minor — 現在仕様と受入条件に旧記述が残っています。**  
   根拠：
   - [現在仕様:413](/home/sano/work/forge-mknot/methods/design/overview.md:413) の case/71「約55 s」は、最終記録の **72.1 s** と不一致。
   - [現在仕様:687](/home/sano/work/forge-mknot/methods/design/overview.md:687) と [手順:31](/home/sano/work/forge-mknot/procedures/nozzle-design-workflow.md:31) の「出力が完全一致」は、`pileup_locs_interior` の選択・順序変更という例外を落としています。
   - [plan:167](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:167) は「3倍以上のまま」、[試験:27](/home/sano/work/forge-mknot/design/tests/run_moc_jit_tests.py:27) は「AWS A・3倍」と記載。§6 の改訂済み条件と矛盾します。  
   **対案：** 現在仕様を §6.6 の最終値・例外・改訂条件に統一し、旧条件には失効表示を付けてください。`plans/README.md` の「同点の規則を外す」も「旧版の同点規則を置換」とすると実装と対応します。

実装・証拠について確認した内容は次のとおりです。

- 指定の `main...HEAD` を取得しました。ただしローカル `main=77318d0e` は起点より古く、過去の統合変更を含みます。本計画の実装照合は **`4b1708ac..HEAD`**。この範囲に CUDA ソルバ変更はなく、node 境界半割面・周期 seam・float32 流束の変更はありません。MOC は float64 で、軸処理・幾何的棄却・MOC の seam 除外を確認しました。
- 整数格子、区間ごとの谷、確認後の比較、一度だけの同点集合、実効入力変更時の評価破棄、再固定点の計時は §4 と対応しています。既定 numba と対応ガスの制限も明記されています。
- [最終証拠](/home/sano/work/forge-mknot/case/70.m5_air_wt_demo/_mknot_fast_final/README.md) のカード2枚を既存ツールで再照合し、両方 **`VERDICT: PASS`**。保存評価表から選択を再計算して **2.05／2.80** を再現し、最終3点と YAML の実効入力、計時内訳の合計も一致しました。
- [基準計時](/home/sano/work/forge-mknot/case/70.m5_air_wt_demo/_mknot_fast_final/time_base.txt:1) と [最終計時](/home/sano/work/forge-mknot/case/70.m5_air_wt_demo/_mknot_fast_final/time_final.txt:3) は、約 **2.40／2.42倍**を裏付けます。証拠版 `cba2914f` と HEAD の生産コードに差分はありません。
- V1 合成試験は `ALL OK`。JIT 合成試験は、キャッシュ書き込みと参照ファイル生成を避ける実行方法で **253項目、`ALL OK`**。実データ V5 は保存結果を確認し、今回再実行はしていません。
- 本件は CFD 0 step の設計検証です。カードの `PASS` を CFD の収束・準定常性とは扱っていません。本計画の主張に、`NOT CONVERGED`／`DRIFTING` を「一致」の根拠にした箇所は確認していません。
- §5.1 に `L_c` 再利用の後回しは残っています。今回の3指摘を未完了項目として追加し、修正後に閉じてください。`methods/index.md` の既存リンクは有効です。

**推奨：現行の選択規則と JIT を採用し、上記1→2→3を修正した後、レビュー記録・完了条件・索引を同期して `accepted` へ移す。**

証拠: 足りる

指摘数: Critical 0 / Major 0 / Minor 3
