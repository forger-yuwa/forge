**未成立 (incomplete)** — incomplete: 証拠不足: `cba2914f` に対応する M5・M5.5 のカード全体の検証記録（最終三点、寸法、計時内訳、壁時計、実効入力）。 / AWS の新規則 V3 比較表、および §6.4 の倍率・81秒・メモリピークを直接照合できる元記録。確認できた実行出力はカード全体117秒／停止まで71秒等に限られます。

# codex レビュー: tooling-nozzle-mknot-fast (result)

- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-15
- **commit**: `cba2914f` (feature/mknot-fast)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.0 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-mknot-fast-result.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M1/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

選択規則と JIT の実装に、今回確認した範囲で新たな数値的欠陥は見つかりませんでした。
ただし、Newton 修正後の最終版に対するカード全体の受入確認が閉じておらず、現状での `accepted` 移行は認められません。
3 倍という速度条件の撤回は、明記されたユーザ決定として扱い、指摘対象にしません。

1. **Major — 最終版の受入証拠が揃っていない。**  
   根拠: [plan:218](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:218) は「最終の構成で V4 を回し直す」としています。しかし、掲載された速度測定は修正前の `f7e185c4`（[plan:285](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:285)）。M5.5 はカード生成が途中で停止しています（[plan:288](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:288)）。寸法 PASS の代替記録も、`7713f84c`・numpy 版です。  
   最終 `cba2914f` の V5・既存試験ログは確認できましたが、V5 は固定した M_K の比較であり、**掃引全体の選択・固定点・寸法・計時を通した確認にはなりません**。§5.1 にこの未完了項目も明確には残っていません。  
   **対案:** 最終 commit で M5・M5.5 のカード生成を通し、最終三点の成立、寸法、選択値、計時内訳と壁時計を保存する。M5.5 の格子設定変更は検証入力に明示し、§5.1 に完了まで残す。

2. **Minor — 条件が違う境界からの距離を報告している。**  
   根拠: [plan:210](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:210) は、旧境界との距離を示す条件を「r_t と実効 geometry が同じときだけ」と限定しています。一方、[plan:224](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:224) の新旧 r_t は `0.024974004` / `0.024973968 m` と異なり、そのまま距離 `0.078` を掲げています。旧カードの実値 `0.024973968248658388 m` も確認しました。  
   **対案:** この距離を削除する。同一条件で再評価しない限り、最終出力の安全距離の裏付けには使わない。

3. **Minor — 点の堆積位置に関する受入条件を、試験側だけで弱めている。**  
   根拠: [plan:154](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:154) は位置一覧の「集合と件数が同じ」を要求しています。しかし [run_moc_jit_tests.py:554](/home/sano/work/forge-mknot/design/tests/run_moc_jit_tests.py:554) は、堆積が20件を超えると集合の比較を免除します。実装では先頭20件の選び方も変わっています。  
   **対案:** 診断用の代表20件は変更可と契約を明記し、試験では切り詰め前の検出対集合を比較する。現在の「完全一致」という説明も対象を限定する。

4. **Minor — 現在仕様・残作業表に失効した説明が残っている。**  
   根拠: [methods/design/overview.md:413](/home/sano/work/forge-mknot/methods/design/overview.md:413) は、反例によって使い回しを見送った後も「M_K は長さを変えない」と一般化し、V3・V4 を未了としています。[plans/README.md:28](/home/sano/work/forge-mknot/plans/README.md:28) は寸法ゲート待ち、§5.1 の #3b・#4 も古い状態です。[plan:103](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:103) の「JIT は後回し」も残っています。  
   **対案:** 現行仕様、完了事項、今回残る検証、後回しの L_c 再利用を揃える。`methods/index.md` は既存ファイルへの入口が維持されており、今回の追加修正は不要です。

確認できた裏付けは次のとおりです。

- 指定の `main...HEAD` 差分を取得。過去の変更も含むため、計画の分岐点 `4b1708ac..cba2914f` で今回の範囲を特定しました。今回の差分に CUDA ソルバ本体の変更はありません。
- V1 合成試験は `ALL OK`。JIT 合成試験は、書込みを避けるためキャッシュ無効・参照コードをメモリ内で読み込んで実行し、**253項目成功、exit 0**。
- 修正後 V5 の保存ログは `ALL OK`。最大位置差約 `2.05e−12 r_t`、壁差約 `1.56e−12 r_t`。既存19本の exit 0 も確認。ただし一部に入力不足による省略があります。
- M5.5 寸法記録の SHA-256 は提示値と一致。拡大部長差 `4.4465861e−5 m`、物理スロート径差 `1.1534500e−8 m`、`dimensions.pass=true` を確認しました。
- 本計画は CFD 0 step の設計処理の検証です。`CONVERGENCE_VERDICT`／準定常 VERDICT の欠如自体は欠陥に数えていません。

**推奨:** 現行の単純化・JIT 方針を維持し、最終版のカード全体の検証と証拠保存を先に完了させ、上記の記録を直してから再レビューしてください。

証拠不足: `cba2914f` に対応する M5・M5.5 のカード全体の検証記録（最終三点、寸法、計時内訳、壁時計、実効入力）。
証拠不足: AWS の新規則 V3 比較表、および §6.4 の倍率・81秒・メモリピークを直接照合できる元記録。確認できた実行出力はカード全体117秒／停止まで71秒等に限られます。

指摘数: Critical 0 / Major 1 / Minor 3
