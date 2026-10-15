# codex レビュー: tooling-nozzle-card-local-mesh (plan)

- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `0a92b976` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan-2.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M7/m0
- **focus**: §4.4 (NS の本段・判定・出力の新しい既定) とその「修正」「判定の基準をケースごとに・未設定は判定不能」、§5.1 #4a (forge の出力の追加: fieldOutStepInterval・協調停止・節点の出力・境界の質量流束の積算) と #4b (check_settled.py・settle_criteria.yaml・停止規則 v2・生成器の既定)。#4a は実装を始めている。検証計画 (§6) に #4a・#4b の合格条件が足りているかも見ること
- **extra**: `solver_density_cuda/tools/check_settled.py`, `notes/reviews/2026-10-15-m55-settle-and-monitor-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

§4.4 の方向性は妥当ですが、誤合格を防ぐ入力検証、壁熱量の取得、協調停止の完了保証が受入条件まで具体化されていません。残差判定の試作では、短すぎる履歴と末尾の急増を `SETTLED` にすることを再現しました。ファイル変更・forge 起動は行っていません。

1. **Major — 窓の時間幅・欠測・必須残差列を検証する契約が不足している。**

   根拠: [check_settled.py:45](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:45)、[同:68](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:68)、[同:84](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:84)。現行関数への読み取り専用の再現結果は次のとおりです。

   | 入力 | 実際の判定 |
   |---|---|
   | step 0〜10 の一定値、計11点だけ | `SETTLED` |
   | 約4万stepの欠測を挟み、末尾11点だけ一定 | `SETTLED` |
   | `rms_ro` しか存在しない一定系列 | `SETTLED` |
   | 2,000step間隔、0〜40,000の完全な一定系列 | `INDETERMINATE` |

   最後の例は窓の左端を除外するため、2万step窓が10点になり、仕様の「端点込み11点」と食い違います。直前窓が不足すると上昇検査を省略して合格でき、恒等ゼロの列も設定との照合なしに除外されます。

   **対案:** `settle_criteria.yaml` に必要列・非活動を許す列・観測間隔・判定区間を解決して記録する。W・2Wの被覆、欠測、stepの単調性、必要列を先に検証し、不足は `INDETERMINATE` にする。上表と「活動中の乱流列が全ゼロ」を必須の負例試験に追加してください。

2. **Major — 分位幅と回帰傾きだけでは、有限値の残差爆発を合格させる。**

   根拠: [check_settled.py:82](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:82)。保存量・乱流・成分の9列を持つ40,001点の系列で、最後の `rms_ro` だけを **10⁸倍**にしても、総合判定は **`SETTLED`** でした。末尾窓の傾きは約0.00240桁、p95−p5は0、窓平均の上昇は0.0004桁となり、すべて閾値内に入ります。

   「非有限を拒否する」「`run_0002` を落とす」だけでは、この誤合格を防げません。

   **対案:** 平滑な準定常判定とは別に、生系列の末尾急増・短時間のバーストを拒否する検査を設け、閾値を基準ファイルに明記する。単発急増、短い反復バースト、持続上昇を負例にし、監視中に発生した場合は保持時間を取り消してください。

3. **Major — 必須にした `wall_heat_nodes` のデータ取得条件が #4a・#4b に接続されていない。**

   根拠: [plan:93](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:93) は `iface_Qf_eff` を要求しますが、#4a の節点出力には熱量がなく、生成器の変更項目にも `output.interfaceDiag: 1` がありません。[conjugateWall.cpp:206](/home/sano/work/forge-cardmesh/solver_density_cuda/conjugateWall.cpp:206) は同設定がなければ生成を省略します。さらに [同:265](/home/sano/work/forge-cardmesh/solver_density_cuda/conjugateWall.cpp:265) では、実効熱量の適用条件が定常・等温壁などに限られます。

   実際に `case/71.m55_combustion_wt/run_0010_ns_b0/res_wall_3_86000.h5` を確認すると、`ifaceFw`・`ifaceRraw` 等はありますが、**`iface_Qf_eff` はありません**。既存runの出力をそのまま新規則の完全な再生試験には使えません。

   **対案:** 各量について「必要な出力設定・入力ファイル・対応ID・単位・適用条件」を定義し、準備時に自動設定・検査する。等温壁では `interfaceDiag` を有効化し、断熱壁などはプリセット解決時に適用可否を明示する。欠測を理由に必須量を黙って除外しないこと、総熱量一定でも局所分布が変わる負例を落とすことを受入条件にしてください。

4. **Major — 協調停止の成功条件が、出力失敗まで含めて定義されていない。**

   根拠: [plan:138](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:138) の試験は、主に「全場の存在と残差stepの対応」です。実装中の [monitorOutputs.cpp:325](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:325) は流量CSVの `flush()` 後に書込み失敗を確認しません。また [同:271](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:271) は完了通知の書込み・rename失敗をログだけにし、[main.cpp:3595](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:3595) はそのまま終了コード0へ進みます。

   容量不足が既に問題になった運用で、全場だけが存在して監視系列や完了通知が欠けるケースを取りこぼします。

   **対案:** 全場・必要な境界・節点・流量・残差の書込み成功を確認した後だけ完了通知を発行する。通常出力も読取り可能な完成状態を識別できるようにし、監視側は実際の停止stepまで再評価する。非出力stepでの停止、書込み失敗、通知失敗、要求後に条件が崩れる場合を試験し、失敗時は `REACHED` にしない契約を追加してください。

5. **Major — #4a の流量試験は異なる定義の量を正解としており、節点監視の同値性試験も不足している。**

   根拠: [plan:138](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:138) は数値流束和と全場の折れ線積分の相対差1e−4を合格条件にしています。しかし [convectiveFlux_boundary_d.inc.cuh:188](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:188) の流束は境界状態と双対面を使います。節点状態の折れ線積分とは、離散化・境界状態が異なり、定常になっても一致は保証されません。

   また `eta_line` は [nozzle_report.py:138](/home/sano/work/forge-cardmesh/design/forge_design/report/nozzle_report.py:138) の補間です。「節点と重みが同じ」だけでは、節点でMachを計算して補間するか、速度・音速を補間してからMachを計算するかの差が残ります。

   **対案:** 流量出力の正解は、同じ流束評価時点の面別 `massflux` の独立集計とする。符号・2π・軸端の有限双対面・境界角・周期面除外を検査し、折れ線積分との差は別の診断として扱う。節点値は同時出力の全場とビット比較し、Mach線・出口平均・流れ角は既存の抽出関数との同値性を、演算順まで固定して検証してください。

6. **Major — 新判定への移行範囲と、既存の区間管理・呼び出し側の移行が未設計。**

   根拠: [plan:103](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:103) はリポジトリ全体を「`SETTLED` だけが収束」に替える一方、プリセットと検証対象は主にノズルです。[run_case.sh:66](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/run_case.sh:66) は旧判定を実行し、[deltastar_loop.py:653](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:653) も旧判定を直接呼びます。AGENTSとStopフックだけ変更すると判定が分裂します。

   既存の [区間管理plan:76](/home/sano/work/forge-cardmesh/plans/active/tooling-convergence-and-wall-resolution-gates.md:76) が防いだ、異なる方程式・BCの履歴混在も、新ツールへ引き継ぐ仕様が必要です。現在の本段manifestは [runner_axismach.py:2405](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2405) で監視終了後に書かれます。

   **対案:** 共通の判定APIと結果形式を先に定義し、オンライン監視・終了後判定・ラッパー・受け渡しゲートを同じAPIへ移行する。結果には基準のhashだけでなく、起動・メッシュ・判定区間・終端step・使用データの識別を含める。本段開始時に区間を記録し、旧データ、未設定ケース、入口出口のない周期問題、非定常問題の扱いを明文化してください。ノズルの検証完了だけで全ケースの移行完了とはしないことが必要です。

7. **Major — 採用済み修正が現行仕様・検証・完了条件へ統合されておらず、計画どおり実装しても旧案が残る。**

   根拠: §6.1では前回指摘を採用していますが、本文には次が残っています。

   - [plan:52](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:52)：縮流部係数1.0、拡大部0.5。後続の決定は0.6・0.49。
   - [plan:53](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:53)：実寸の下包絡を半径で割る旧手順。採用済み修正は正規化した要求の下包絡と区間内部の検査。
   - [plan:150](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:150)：存在しない `check_mesh_quality.py --bcond`。
   - [plan:183](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:183)：完了条件が依然「V1〜V3、V4は記録」。#4a・#4bの合格が不要。

   `mesh_euler` を対象外とする記述も、niの結合による変更を説明できていません。§6には新しい判定・出力・B0既定化・リミッタ継承の受入項目がありません。

   **対案:** 修正追記を現行仕様へ統合し、旧案は履歴へ移す。§6に上記1〜6の試験、出力有効化前後の数値不変、B0がNS本段だけに効くこと、ρ・P・a・`L_ref`の往復精度と継承、オンラインと再生判定の一致を追加する。§8もその合格を必須にし、関連する生成器・収束ゲートplanとの担当範囲を明記してください。

推奨は、**§4・§5.1・§6・§8を単一の現行仕様へ書き直し、不正入力と出力完了の契約を確定してから #4a・#4b の実装を進めること**です。B0採用、保持4万step、最終全場からのδ抽出という決定は維持できます。

今回の実行では `case/71.m55_combustion_wt/_runs/run_0002_ns/residual_history.csv.gz` に対する試作判定 **`VERDICT: NOT_SETTLED`** を再現しました。ただし、これだけでは上記の合成負例や新しい監視量を検証できません。既存判定も `run_0008_ns_b0` は `PASS`、`run_0010_ns_b0` は `NOT CONVERGED` のままで、新規則の合格へ読み替えていません。保存記録の索引は [case README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21) です。

証拠: 足りる

指摘数: Critical 0 / Major 7 / Minor 0
