# codex レビュー: tooling-nozzle-mknot-fast (plan)

- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `21efe401` (feature/mknot-fast)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.5 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-mknot-fast-plan.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M3/m3
- **focus**: §4.4 (最初の r_t のゲートを N_p だけにし 2N_p は最終の組だけ) の可否、§4.2 の捕捉を最後の 2 つに限ることの正しさ (from_length の最良は最後か最後から 2 番目か)、§6 V2〜V4 の合格条件

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

§4.4 の初回 `2N_p` 省略と、§4.2 の捕捉を最後の 2 件に制限する案は採用可能です。
ただし、粗い探索の候補 1 点と局所ゲートだけでは、既存の選択規則を維持できません。V2 の統合試験と V4 の完了条件も実装前に修正が必要です。

**指摘一覧**

1. **Major — 生産網のゲートは、粗い探索による別の谷の選択を検出できない。**

   根拠: [plan:51](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:51) は粗い網の選択 1 点だけを渡します。既存の [`select_mknot`:519](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:519) は全探索範囲の局所極小を拾いますが、[`resolution_gate`:741](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:741) は入力点を中心とする局所窓しか評価しません。

   現行コードで次の**合成反例**を再現しました。全点で制約を満たし、生産網の `N_p` と `2N_p` は同じ評価関数です。

   | 経路 | 選択 M_K | 生産網の θ_max | ゲート |
   |---|---:|---:|---|
   | 粗い探索 → 生産網の局所ゲート | 2.00 | 15.0° | 全項目 PASS |
   | 生産網の既存探索 | 3.45 | 14.5° | — |

   評価関数は、粗い網で `min(15 + 2(M−2)², 15.02 + 400(M−3.45)²)`、生産網で第 2 項の定数を `14.5` に変更したものです。粗い探索でも局所極小は 2 個見つかっていますが、選択 1 点に縮約した時点でもう一方を捨てます。

   **対案:** 今回の既定経路では `N_p` の既存 `select_mknot` を維持してください。粗い候補 1 点への置換は外します。V2 にこの順位逆転例を追加します。V3 の 2 ケースで θ 差 ≤0.01°を満たしても、一般の入力についてこの欠陥は解消されません。

2. **Major — §4.4 は成立するが、初回の仮選択と最終保証を分ける実装・試験が不足している。**

   根拠: [`wt_card.py`:1952](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:1952) は条件次第で `gate_first` を最終ゲートとして再利用します。また、[`mknot_select.py`:735](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:735) は渡された表を解像度検証なしで seed します。初回を `N_p` だけに変更する場合、その結果を従来の完全な `gate_first` と同じ扱いにはできません。

   さらに、現行コードは移動量 ≤0.01 でも、物理スロート径を解く経路では出力組の固定点残差を再検査します（[`wt_card.py`:1987](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:1987)）。§4.4 の「>0.01 ならやり直す」だけではこの条件を表現しきれていません。

   **対案:** 初回は「仮選択」と明記し、最終ゲートと別の状態として保持してください。V2 は `resolution_gate` 単体ではなくカード生成の制御経路を通し、最低限、以下を確認します。

   - `r_t` 不変でも完全な `N_p/2N_p` ゲートを必ず実行する。
   - `r_t`・実効 geometry の変更後は旧評価を再利用しない。粗い網の行も生産網へ seed しない。
   - 移動量 >0.01 の再固定点、≤0.01 の残差再検査、再試行上限での例外を通す。
   - 最終 YAML・card_log・最終ゲートの評価組が対応する。

   **初回 `2N_p` は不要です。** ただし「最後に 1 回だけ」ではなく、固定点をやり直すたびに、その更新後の組で完全ゲートが必要です。

3. **Major — V4 が目的達成を判定せず、既存の計時も再固定点の費用を除外している。**

   根拠: [plan:120](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:120) は 10 分未達を不合格にしません。このままでは短縮できなくても完了条件を満たせます。

   また、[`wt_card.py`:1972](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:1972) は `mk["t"]` を加算した**後**に `solve_loop` を呼びます。§4.4 によって再固定点が増えると、その追加時間が `M_knot_rule.elapsed_s` から抜けます。既存 M5 の記録は全体 5404 s、M_knot 5339 s、やり直し 0 回なので、この欠落を露呈しません（[check_final.txt:26](/home/sano/work/forge-mknot/case/70.m5_air_wt_demo/_mknot_card_v4_3d/check_final.txt:26)）。

   **対案:** M_knot に起因する再固定点・再ゲートを含めて計時し、カード全体の時間も併記します。V4 は同一 AWS 環境・同一入力で測り、両ケースで 600 s 以下を完了条件にしてください。未達なら機能上の検証成功と性能目標未達を分け、plan は未完了にします。

4. **Minor — 並列数の式に、0 worker・メモリ取得不能時の扱いがない。**

   根拠: [plan:68](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:68) の式は、例えば使用可能 7 GiB・1 評価 5.6 GiB なら 0 を返します。一方、[`Evaluator`:299](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:299) は 0 を 1 に補正します。これでは確保したつもりの 30% の余裕を失います。`MemAvailable` は取得失敗時に `None` になります。

   **対案:** 算出値 0 は開始前に明示的に停止し、取得不能時の方針と手元 PC の既定上限を固定してください。捕捉を 2 件にしても、新しい網の構築中には旧 2 件が残るため、並列数は変更後のプロセスピーク実測から決めます。

5. **Minor — V2(iv) の「行がビット同一」「例外」は現行 API と一致しない。**

   根拠: [`evaluate_mknot`:215](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:215) は捕捉失敗を含む例外を行の `error` に変換し、221 行で実行ごとに変わる `time_s` を追加します。行全体のビット同一も、呼び出し元への直接の例外伝播も期待できません。

   **対案:** 時間などを除いた数値・判定フィールドの完全一致を条件にします。捕捉対象外の場合は `error`、`ok=False`、明示的な理由を確認し、選択経路がそれを採用しないことまで試験してください。

6. **Minor — 設計文書に相互矛盾と旧方針が残っている。**

   根拠: [plan:55](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:55) と実装項目 #3 は初回も完全ゲート、§4.4 は初回 `N_p` のみです。[plans/README.md:28](/home/sano/work/forge-mknot/plans/README.md:28) には、既に却下された「生産網の 11 点だけ確認・失敗時 fallback」が残っています。

   **対案:** レビュー採否の決定時に §4・§5・§6、`methods/design/overview.md`、索引を同じ方針へ揃えてください。

**推奨**

**`N_p` の既存探索を維持し、捕捉を最後の 2 件に制限してメモリに応じたプロセス並列を使い、初回 `2N_p` を省略する経路**を推奨します。実装前の修正順は、①候補探索の維持、②初回・最終状態の分離と V2、③再固定点を含む V4、④並列条件と文書同期です。

捕捉 2 件は妥当です。厳密には、現行 `from_length` の**正常返却は最後の評価**です。改善停止時は直前が `best` ですが、その直前が許容差以内なら次の評価へ進まないため、改善停止後は未収束例外になります（[`runner_axismach.py`:653](/home/sano/work/forge-mknot/design/forge_design/evaluate/runner_axismach.py:653)）。現行分岐を抽出した 256 通りの残差系列でも、正常返却 125 通りはすべて最後、残り 131 通りは例外でした。

V3 の最終完全ゲート・寸法 PASS・θ 差 ≤0.01°は、2 ケースの受入条件として妥当です。目的も未解決の短縮課題として正しく、既存の交差検査高速化 plan との重複ではありません。本変更は CPU 側の設計処理なので、CUDA/float32、node/cell、block-DPLUR、周期境界の回帰計算を追加する必要はありません。CFD の収束・定常性を主張する変更でもありません。

**証拠の確認範囲**

現行コード、既存 M5 の生成物・ゲート値・計時を照合し、YAML と card_log の整合検査は `ok=True` でした。既存試験は書き込みを要する 1 件を除いてメモリ内で実行し、96 チェック成功、失敗 0。試験 (n)(o) と上記の反例も確認しました。ファイル変更はありません。これは計画段階の判定であり、未実装の新経路や AWS の新旧性能比較を合格としたものではありません。

証拠: 足りる

指摘数: Critical 0 / Major 3 / Minor 3
