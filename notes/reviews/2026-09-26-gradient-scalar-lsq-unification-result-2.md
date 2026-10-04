# codex レビュー: gradient-scalar-lsq-unification (result)

- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **stage**: `result` (diff base `36d8ba03`)
- **date**: 2026-09-26
- **commit**: `926837c1` (feature/sern-design)
- **codex**: effort `high`, 4.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: Phase 1 (opt-in) の result 2 回目。1 回目 NO-GO の 7 件は §5.1 #5f–#5m で対応済み。S2 の収束ゲートの基準はユーザ決定で『起点 run の床』から『同じ現行バイナリの gg 双子』へ差し替えた (§6 の差し替え規則、判定前に固定)。この基準での判定と、Phase 1 を done にして Phase 2 (既定化、前提 #2g) に進めるかを判断してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
新基準による S2 の残差比 PASS は保存記録と整合し、LSQ の中核実装にも明確な欠陥は見つかりません。  
ただし、区間判定の不具合が残り、一次記録の検証と Phase 2 が未完了のため、現状の plan を `done`／`accepted` にはできません。

指定の `36d8ba03...926837c1` の diff を取得して確認しました。係数共有、差分形、境界 incidence の skip、最大4変数のチャンク、周期二重 gather の回避、軸対称での体積除算除去は §4 と整合しています。新たな絶対ゼロ割閾値の導入はなく、既定 `gg`・cell の GG 固定も維持されています。`f99f236d` 以降、対象の solver 本体に追加変更はありません。

保存された [S2_GGREF.txt](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/S2_GGREF.txt:1) では、指定された最終段の全活動列の最大比は次のとおりです。

| ケース | 最大の lsq／gg 末尾残差比 | 新基準 |
|---|---:|---|
| case/48 | 1.419 | PASS |
| case/40 | 1.024 | PASS |
| case/39 | 1.000 | PASS |
| case/16・0476 | 1.007 | PASS |
| case/16・0482 | 1.002 | PASS |

同ファイルの gg 通常判定は全件 `NOT CONVERGED (stalled/plateau…)` です。これは今回の plateau 許容規則とは矛盾しません。**旧床への復帰を再び要求する理由にはしませんが、通常判定で収束したという意味でもありません。**

1. **Major — 前回 M5 の修正後も、明示 gg→lsq が同一区間になる経路が残っています。**

   [stage_manifest.py:166](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:166) は、`PyYAML` の読み込み失敗時に空辞書を返します。`mesh.scalarGradient` は YAML 解析側にしか追加されていないため、この経路で識別が消えます。

   ファイルを書かないメモリ上の再現試験で確認しました。

   | 条件 | 明示 gg／lsq のキー | `segments()` |
   |---|---|---:|
   | 通常 | `gg`／`lsq` | 2区間 |
   | `import yaml` が `ImportError` | 両方欠落 | **1区間** |

   通常環境では追加4試験と既存 chi 試験は PASS しましたが、この条件は試験されていません。[solver-settings.md:374](/home/sano/work/forge-sern-design/procedures/solver-settings.md:374) の「別区間として扱う」という保証も成立しません。これは既定化後だけの問題ではなく、現在の opt-in 切替にも影響します。

   **対案:** `PyYAML` を区間判定の必須依存にし、読み込めない場合は明示的に判定を停止してください。同条件の負の試験を追加し、#5j を再度閉じる必要があります。

2. **Major — 新しい S2 合格条件の全項目を独立に検証できる証拠が不足しています。**

   対象 `run_095*` はローカルに存在せず、[case/48 の索引](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:55) と [case/16 の索引](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393) も AWS のみと明記しています。今回、AWS API への接続は失敗しました。したがって残差 CSV、`CONVERGENCE_VERDICT.txt`、準定常系列・VERDICT の原本を再検査できていません。

   特に FCT の新条件は「同じ判定、かつ lsq の remainder が gg 以下」です。しかし `8.45e-2` という丸め値の一致だけでは大小関係を証明できません。[M3_fct.txt:29](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/M3_fct.txt:29) が保存しているのは感度条件の `B` であり、収支 checker の完全な出力ではありません。

   **対案:** 既存 run から、判定区間付き残差 CSV、通常判定、新床比判定、準定常系列と実行引数・VERDICT、FCT の丸め前 remainder と項目別判定、入力・バイナリ provenance を回収してください。追加 CFD 計算より先に、既存証拠をレビュー可能にするべきです。現時点で数値の誤りを断定しているわけではありません。

3. **Major — Phase 1 の完了と、この plan 全体の完了が混同されています。**

   [plan §4.4](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:67) は、Phase 2 の既定化・S4・設計 DB・result レビューまでを同じ plan に含めています。[§5.1 #2g・#2h](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:92) と [#6](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:118) は未完了で、参照先 `plans/active/tooling-stage-manifest-launch-binding.md` 自体も存在しません。

   **対案:** Phase 1 の判定を §6.1 に記録し、plan 全体は `in_progress` に維持してください。#2g を起票・完了し、S4・設計 DB・既定化後の検証を終えてから `accepted` に移すのが、現在の「分割しない・2段ゲート」という設計に沿います。

4. **Minor — 最新判定と、現在仕様・run 索引・残作業表が同期していません。**

   [methods/gradient.md:80](/home/sano/work/forge-sern-design/methods/gradient.md:80) と [solver-settings.md:372](/home/sano/work/forge-sern-design/procedures/solver-settings.md:372) は S2 を「条件付き・未合格」と記載し、case 索引にも旧基準の「判断待ち」が残っています。§5.1 #2a–#2e は後続行に実施記録があるのに状態未更新、#2i は旧 S3 数値のままです。既知の chi 文書不整合 #5e も未解消です。`methods/index.md` の既存リンクには問題ありません。

   **対案:** 旧 FAIL 記録を保存したまま、「旧規則の結果」と「差し替え後の現行判定」を明示し、各索引と残作業行を同期してください。FCT の精度未検証は、保存収支の非悪化判定と別に残す必要があります。

**推奨は、既定 `gg`・plan 全体の `in_progress` を維持し、①区間判定の修正、②既存一次記録の回収・照合、③Phase 1 の正式な完了記録、の順に進めることです。** その後に #2g を前提として Phase 2 へ進めてください。今回の新基準に対して、旧床復帰や FCT 精度審査を追加の合格条件として復活させる必要はありません。

床比ツールの9例はファイル入出力をメモリに置き換えて再実行し、すべて PASS しました。ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
