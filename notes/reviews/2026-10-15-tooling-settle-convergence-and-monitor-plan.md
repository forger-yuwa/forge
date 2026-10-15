# codex レビュー: tooling-settle-convergence-and-monitor (plan)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `4a7374a3` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-settle-convergence-and-monitor-plan.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **extra**: `solver_density_cuda/tools/check_settled.py`, `design/forge_design/evaluate/stop_rule.py`, `notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan-2.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

目的と分割方針は妥当で、前回の主要指摘も計画に取り込まれています。ただし、停止時の判定窓、再生試験の成立性、既存ゲートの移行には未確定の契約が残ります。
以下を優先順に計画へ反映してから実装してください。ファイル変更・forge 起動は行っていません。

1. **Major — 非出力 step での協調停止と、判定窓の完全性が整合していない。**

   根拠: [plan:43](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:43)、[plan:58](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:58)、[plan:65](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:65)。

   2,000 step 間隔の系列を **80,001 step** で停止すると、直前窓 `[40,001, 60,001]` にある通常観測は **10点**です。停止時の追加出力を含む末尾窓は11点ですが、直前窓は増えません。端点込み11点の検査をそのまま適用すると、一定系列でも停止後の再判定が成立しません。

   また、[main.cpp:2571](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:2571) は場を `iStep+1` で書いた後、残差を `iStep` で記録します。残差はさらにバッファリングされるため、通常監視でもファイルの存在だけでは同じ終端まで揃ったことを保証できません。

   **対案:** 通常観測の格子上で W・2W と保持期間を判定し、停止時の追加観測には別途、有限性・急変・最終値の条件を適用する仕様を定めてください。残差番号の対応と更新前後の評価時点も明記し、必要な出力一式の完成を識別できるようにすること。V4 に「一定系列を80,001 stepで停止して正常に合格する正例」を追加してください。

2. **Major — V3 の既存データでは、新しいプリセット全体を再生・校正できない。**

   根拠: [plan:103](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:103)、[case README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21)。

   読み取り確認では、次の制約があります。

   - 対象の既存 run に `boundary_massflux.csv` はありません。旧 `stop_rule.json` の流量は全場の折れ線積分で、新規の数値流束和とは別定義です。
   - `run_0010_ns_b0/res_wall_3_86000.h5` に **`iface_Qf_eff` はありません**。
   - 手元の `run_0005_ns`・`run_0010_ns_b0` の全場は各1枚です。READMEにも途中全場の削除が記録されています。
   - 保存済み `v5c_series.csv` は最大偏差と片振幅のみで、線の位置ごとの変化を再生できません。

   この状態で全条件を AND すれば、既存 run はデータ不足になります。それを「未整定の負例を落とせた」と数えても判定性能の検証にはなりません。

   **対案:** V3 を「残存データで可能な成分別再生」と「新出力を備えた run での全条件検証」に分けてください。新規 run で数値流束・節点線・壁熱量を取得し、閾値を固定した後、別の系列で正例・負例を確認すること。各試験で、不合格理由が閾値超過かデータ不足かも検査してください。

3. **Major — δ 時系列廃止後の受け渡しゲートの置換条件が不足している。**

   根拠: [plan:65](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:65)、[plan:77](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:77)、[deltastar_loop.py:1006](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1006)。

   現行 `handoff_gate` は残差判定だけでなく、δ_E・δ_target の `two_window_check` を必須にしています。記録欠落は上書きでも通せません。旧収束判定の呼出しだけを差し替えると、新方式で `SETTLED` になっても必ずここで止まります。一方、ゲート全体を `SETTLED` に置き換えると、壁解像・Euler参照・抽出被覆・次の壁の健全性まで失うおそれがあります。

   **対案:** 新方式では旧残差判定・旧流量判定・δの2窓判定を何で置き換えるかを条件ごとに明記してください。壁解像、参照の同一性、最終δ抽出の被覆、次の壁の検査は維持し、旧方式との分岐と上書き規則も定義すること。V4/V5の後に、受け渡しゲートまで通す移行試験を追加してください。

4. **Major — リミッタ基準値の「記録」だけでは、異なる作用素の履歴連結を防げない。**

   根拠: [plan:73](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:73)、[stage_manifest.py:266](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/stage_manifest.py:266)、[main.cpp:1793](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:1793)。

   既存APIで、他の設定を同一にして `limiterRoRef` を **1→10** と変更した2段を作ると、**hard key は同一、区間数は1**でした。現行の区間判定には基準値4個が含まれていません。記録欄へ追加するだけでは、この挙動は変わりません。

   **対案:** forge が確定した ρ・P・a・L_ref を区間比較の hard key に含め、未確定値を同一扱いしない仕様を追加してください。自動決定後に本段の記録を確定する順序も必要です。V5には「同値の継承は同一区間」「各基準値の変更は別区間」「自動値が不明なら連結不可」を追加してください。起動識別は既存の [launch-binding plan](/home/sano/work/forge-cardmesh/plans/active/tooling-stage-manifest-launch-binding.md:30) との担当範囲も明示すべきです。

5. **Major — 必須の壁熱量判定に、ゼロ近傍で成立する尺度がない。**

   根拠: [plan:60](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:60)。

   同節は固定尺度を使うとした直後に、壁熱量だけ「節点ごとの相対変化」としています。分母・絶対許容差・符号反転の扱いが未定義です。局所熱量がゼロなら完全に一定でも相対変化を計算できず、ゼロ近傍では小さな揺らぎが判定を支配します。等温壁であることは、全節点の熱量がゼロから離れる保証にはなりません。

   **対案:** 節点ごとに `|ΔQ| ≤ Q_abs,i + ε_rel·Q_ref,i` のような絶対・相対許容差を定義してください。固定の熱流束公差から局所面積で `Q_abs,i` を作れば、局所性も保持できます。一定ゼロを正例、許容内のゼロ交差と許容を超える局所変動を別々の試験にしてください。

6. **Minor — V5 の1 GB条件には、対象格子と容量の集計範囲が必要。**

   根拠: [plan:105](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:105)。

   実測した `case/71.m55_combustion_wt/run_0010_ns_b0/res_86000.h5` は **100,814,972 bytes**。同サイズなら10枚だけで **1,008,149,720 bytes**となり、境界・監視・残差出力を含める前に1 GBを超えます。

   **対案:** 節点数・精度・`output.level`・初期出力の扱い・GB/GiB・残差CSVを含めるかを固定した試験条件にしてください。case/71の規模にも1 GBを要求するなら、出力間隔か保持枚数の調整が必要です。

推奨は、**基本方針を維持し、上記の契約と受入試験を§4・§5.1・§6へ統合してから実装すること**です。新しい判定自体は既存の準定常ツールと目的が異なり、追加する理由があります。数値カーネルの変更やδ時系列の復活は必要ありません。

証拠として、コード・作業差分・保存データを確認し、区間連結と停止窓の問題を再現しました。`run_0010_ns_b0` の保存系列を既存 `check_quasisteady.py` で再実行した結果は `OVERALL: ALL STEADY`。保存済み `check_convergence` 判定は `run_0008_ns_b0` が `PASS`、`run_0010_ns_b0` が `NOT CONVERGED` です。これらを新プリセット全体の合格とは読み替えていません。

証拠: 足りる

指摘数: Critical 0 / Major 5 / Minor 1
