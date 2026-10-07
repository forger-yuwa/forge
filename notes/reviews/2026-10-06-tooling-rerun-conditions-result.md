# codex レビュー: tooling-rerun-conditions (result)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `7d7214f0` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

入力検査と段階起動に、再現できる不具合があります。正常な Pt スケール計算自体は確認できましたが、主要な検証 run の原データを参照できず、結果の受入れは認められません。ファイルは変更していません。

1. **Major — 非物理的な入力を受理し、負密度の初期場がスケール検査に合格する**

   根拠: [`rerun_conditions.py:393`](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:393)、[同:695](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:695)、[同:749](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:749)。

   実在する `run_0094_ns_c2pin_pass2_ext6k` を読み、保存先だけをメモリ内 HDF5 にして再現しました。

   - `--Pt -1 --keep-Ps --scale-ic pt`：受理され、**全 121,250 点が負密度**になっても `scale_fields` は合格。
   - `--Pt 0 --keep-Ps --scale-ic pt`：`ZeroDivisionError`。
   - `--Tt -1`、`--Ps -1`、`--cfl -1`、`--cfl inf`：生成計画が受理される。

   `allclose` は掛け算の正確さしか保証せず、物理的妥当性は保証しません。`field_problems` はスケール前にしか実行されません。また、入口 Y の検査も変更時に限定され、参照 BC の ΣY＝1.0001 を受理しました。

   **対案:** 作成前に指定値・実効値の有限性と物理範囲、常時の Y 範囲・総和を検査する。スケール後にも `field_problems` を実行する。

2. **Major — 受理した YAML の書式によって段階起動が無効になる**

   根拠: [`runner_axismach.py:1006`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1006)、[同:1076](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1076)。

   参照 config の表記だけをメモリ内で変更して確認しました。正当な YAML の `convMethod:  2` と `cfl: 5.0e+0, cfl_pseudo: 5.0e+0` は `build_plan` を通過しますが、soft 段の置換後も **convMethod＝2、両 CFL＝5** のままです。要求した一次化・CFL 低減が黙って失敗します。

   **対案:** YAML の構造に基づいて変更し、各段の実効値を起動前に検査する。空白数・指数表記・block 形式を含む回帰試験を追加する。

3. **Major — 表示する起動手順が推奨 CFL を適用しない**

   根拠: [`rerun_conditions.py:552`](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:552)、[同:604](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:604)、[同:859](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:859)。

   `run_0094` に `--Pt 4400000 --Ps 1789.6 --scale-ic pt` を指定すると、推奨は `full / cfl 1`、生成予定 config は **cfl＝cfl_pseudo＝5** でした。表示される `run_staged_ns(..., stages='full')` は推奨 JSON を読みません。表示どおり実行すると、計画が発散を報告した本段条件になります。

   **対案:** 推奨値と実効値が異なる場合は、必要な `--cfl 1 --steps 60000` を明示して作成前に停止する。表示する実行手順と生成 config の整合を試験する。

4. **Major — 検証結果を独立に監査できる証拠が不足している**

   根拠: [`plan:176`](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:176)、[同:179](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:179)、[case README:96](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:96)。

   対象ツリーと別作業ツリーを調べましたが、検証用 `run_0119`〜`run_0139` の原データがありません。したがって、A3 の STEADY、B3 の DRIFTING、Euler の流量差、AWS の単体試験 SKIP 0 は**未確認**です。数値が誤りと断定する根拠もありません。

   ローカルにある参照 `run_0094` の保存済み判定は **`NOT CONVERGED (stalled/plateau)`**。再実行した判定器は `NO residual_history.csv` となり、こちらも残差から再検証できませんでした。

   また、計画自身が B3 を DRIFTING とする以上、[plan:65](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:65) の「偽のはく離」という原因の断定は裏付け不足です。有限時間内の逆流残存だけでは、長い過渡と別の定常解を区別できません。

   **対案:** 各 run の config、残差全列、manifest、量の CSV、判定コマンドと出力、必要な HDF5 を監査可能にする。台帳は省略名を廃して run ごとに成果物へリンクする。B3 の結論は「規定時間内に準定常未達」に限定する。

5. **Minor — 「DRIFTING だが漸近値あり」という判定分岐に到達できない**

   根拠: [`plan:135`](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:135)、[同:144](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:144)、[`check_quasisteady.py:293`](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:293)。

   判定器は DRIFTING なら即 return し、漸近値計算は STEADY の経路にしかありません。合成系列 `100＋10×0.8ⁿ` でも DRIFTING の detail に漸近値は出ず、別途 `_monotone_limit` を呼ぶと 100 が得られました。**「漸近値が出ない」ことを線形ドリフトの証拠にはできません。**

   **対案:** 準定常判定と外挿診断を分け、計画の解釈を修正する。ユーザが承認した整定長の持ち越しは維持し、「残り約 0.1 %」は予測値との差として記載する。

6. **Minor — 現在仕様と完了管理に更新漏れがある**

   根拠と不整合:

   - [`methods/design/overview.md:927`](/home/sano/work/forge-integ-1005/methods/design/overview.md:927)：廃止した `res_outlet_*` の Ps を出口圧の定義として記載。
   - [同:903](/home/sano/work/forge-integ-1005/methods/design/overview.md:903)、[同:935](/home/sano/work/forge-integ-1005/methods/design/overview.md:935)：`slip` と非乱流時の `roK/roOmega` 許容が未反映。
   - [`plan:88`](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:88)：B3 完了後も「本段 cfl 1 は未検証」。
   - [同:157](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:157)：`run_staged_ns` を変更したのに「既存ツールは変えない」。
   - [`plans/README.md:66`](/home/sano/work/forge-integ-1005/plans/README.md:66)：状態が `draft`。

   **対案:** 現在仕様・影響範囲・残作業の状態を同期する。§5.1 #9/#11 の持ち越し自体は記録されていますが、移動時には継続先の active plan を明示する。

正常経路では、参照出口圧 **2241.9975586 Pa、内部 95/97 点**を再現しました。9 保存量の 0.8 倍変換も最大相対誤差 **4.77×10⁻⁸**、変換後の場検査は問題なしでした。`run_mesh_params_tests.py` は FAIL 0。作成を伴う `test_rerun_conditions.py` 全体と実 forge 経路は今回再実行していません。

**推奨は、active に留めて指摘 1→2→3 を修正し、4 の原データで結果を再監査してから移動することです。** スケール機能を削除する根拠はありません。今回はレビューのみで、plan 未反映です。

指摘数: Critical 0 / Major 4 / Minor 2
