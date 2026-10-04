# codex レビュー: case-plate-annular-cavity-m5 (result)

- **plan**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-20
- **commit**: `b64717a7` (feature/sern-design)
- **codex**: effort `high`, 6.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M9/m1
- **extra**: `case/49.plate_annular_cavity_m5/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`GATES: PASS` が plan の必須条件を保証しておらず、実 run でも誤合格を再現しました。  
入口 EOS、収支評価、格子不確かさ、FEM 出力に未解消事項があり、`accepted` への移動は認められません。

指定の diff を取得し、case 関連実装・共有判定ツール・文書・ローカル成果物を照合しました。`run_0401_..._ext3` など AWS 側の最終正本は本ツリーになく、その最終値の独立再検証はできていません。以下では、実行して確認した事項と文書上の未達を区別します。

1. **Major — TP 生産計算の入口に CPG 前駆分布を使っている**

   根拠: [precursor/gen_runs.py:39](/home/sano/work/forge/case/49.plate_annular_cavity_m5/precursor/gen_runs.py:39) は `thermalMethod: 0` 固定です。[gen_runs.py:359](/home/sano/work/forge/case/49.plate_annular_cavity_m5/gen_runs.py:359) は入口 CSV をそのままコピーします。plan §4.10 の「前駆・入口も TP」と一致しません。

   `case/49.plate_annular_cavity_m5/_local_s2/run_0401_s2_off000_tw20_ext1/inlet_profile_1.csv` を採用 TP 物性で再評価すると、壁際静温は **3232.48 K**、入口分布の最大総温は **3258.45 K**。条件から求めた主流総温 **3165.65 K** を約 **93 K** 上回ります。これは既定条件と整合する TP 境界層を与えた検証になっていません。

   **対案:** TP 前駆を作り直し、plan §4.3 の塞ぎ 3D 検証で入口移送・下流発達を確認してから、生産結果と熱回路係数を再評価する。

2. **Major — 継続 run なら残差の不合格を通してしまう**

   根拠: [check_case_gates.py:92](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_case_gates.py:92) は終了コードを無視し、`RISING` / `DIVERGED` という文字列がなければ合格にします。参照場が収束済みかも確認していません。メモリ内の反例でも、必須残差列が欠けた入力を正本ツールは拒否する一方、この判定式は通します。

   上記 `_local_s2/run_0401_s2_off000_tw20_ext1/` の再実行結果は次のとおりです。

   ```text
   check_convergence.py:
     NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)
   check_cavity_steady.py:
     VERDICT: STEADY (全量)
   check_wall_resolution.py:
     VERDICT: FAIL
   check_case_gates.py:
     GATES: PASS (報告してよい)
   ```

   **対案:** `CONTINUED_FROM` の存在で免除せず、収束確認済み参照に対する `--from-floor`、または同一設定区間の正規判定を採用する。入力不備・判定不能・非ゼロ終了コードを確実に不合格とする。

3. **Major — 判定区間の識別が層流→SST の変更を検出しない**

   根拠: [stage_manifest.py:43](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:43) は `turbulenceModel` を検索しますが、生成 config は `turbulence: {model: ...}` です。

   `model: "none"` と `model: "sst"` だけを変えた config を `stage_key()` に渡すと、**同一キーになることを再現**しました。異なる方程式の履歴を分離するという共有ツールの保証が成立していません。

   **対案:** YAML を構造として解析し、実際の `turbulence.model` と方程式・離散化を変える設定を比較する。層流→SST を別区間にする回帰試験を追加する。

4. **Major — エネルギー収支ゲートが実際には作動していない**

   根拠: [check_case_gates.py:149](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_case_gates.py:149) は JSON の `budget_residual` を読み、欠けていれば検査を省略します。しかし [cavity_eval.py:847](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:847) が保存するのは `field` と `wall` だけで、収支残差は標準出力にしか出ません。確認した既存 JSON にも当該キーはありません。

   また、質量不釣合いの許容は plan §6.4 の **0.1%** に対して実装既定が **2%**。上記 run は末尾平均 **1.334%** で合格しました。

   **対案:** 評価面・壁積分・対流・伝導・粘性仕事・正規化残差を JSON に保存し、欠損や非有限を拒否する。閾値を plan と共通化し、修正後に正本 run を再判定する。

5. **Major — TP の開口伝導流束を CPG の定数比熱で計算している**

   根拠: [cavity_eval.py:300](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:300) は `conditions.json` の `cp=1004.5` を使います。一方、ソルバは [gasProperties_d.cu:65](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:65) で局所 `cp(T)` を使います。

   同じ乾燥空気モデルで再計算すると、1273.15 K の `cp` は **1184.81 J/(kg·K)**。同じ粘性・温度勾配なら後処理の伝導係数はソルバより **15.2% 小さい**値です。「全項を入れて収支が閉じた」という §4.8.1 の根拠を、そのまま採用できません。

   **対案:** ソルバの `thermCond`・`cp` を出力して用いるか、採用 EOS と同じ局所物性で再構成する。近似的な開口収支と、未実装の離散流束収支も区別する。

6. **Major — M5 の格子評価を M9・別形状の不確かさに転用している**

   根拠: [plan §4.11.6:1400](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1400) は、旧形状・M5・CPG 系列の **GCI 0.31～1.6%** を、形状2・M9・TP の不確かさ表に入れ、「全部合わせても数%」と結論しています。さらに、その GCI にはリップ除外積分が含まれ、総入熱の精度とは同一ではありません。

   同じ plan §4.11.3 には `cyl_top` の細分で熱流束 **+14.6%**、なお壁解像未達とあります。キャビティ側の第一層を変えない A/B で、キャビティ側の壁解像まで保証することもできません。

   **対案:** 報告する M9 条件で、キャビティ壁と `cyl_top` を含む細分比較を行う。それまでは総入熱・熱回路係数を暫定値とし、「数%の不確かさ」「確定値」を撤回する。

7. **Major — 半割の妥当性と物理的非定常性の断定が検証範囲を超えている**

   根拠: [plan §4.9.1:1120](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1120) は、時間刻み半減・内部反復倍増・格子細分・観測窓倍増を未実施と明記しています。それでも §5.1 #9 は「半割は妥当」で完了扱いです。M5 の粗い全周 tet での旋回減衰は、M9・形状2・偏心条件の保証にもなりません。

   また、`case/49.plate_annular_cavity_m5/_local_off150/run_0205_off150_s14_ext/` の再判定は **`VERDICT: DRIFTING`**。config は定常計算なので、擬似時間の振動から物理的リミットサイクルと確定することはできません。

   **対案:** plan §6.4 の四条件を対象条件で検証し、それまでは「当該設定で旋回擾乱が減衰」「定常反復が未収束」と限定する。

8. **Major — FEM 出力はリップ帯を均していないのに、均し済みと説明している**

   根拠: [export_fem_bc.py:148](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/export_fem_bc.py:148) の `--eps-mm` はメタデータと説明文に使われるだけで、`qpp_grid()` の積分・分布加工に渡されません。[出力 README:16](/home/sano/work/forge/case/49.plate_annular_cavity_m5/_fem_bc/fem_bc_conc/README.md:16) は「開口端1 mmは均してある」と記しています。

   実際の `qpp_cav_outer_Tw293K.csv` には、その帯内に18層が残り、同じ下流側方位で **172～1267 kW/m²** の値が並んでいます。

   **対案:** 指定幅で面を切断し、方位別の帯積分熱量を保存する均しを実装する。幅を変えたときの分布変化と積分熱量保存を検証し、加工前後を明示する。

9. **Major — 残作業の正本と「完了」表示が実態を保持していない**

   根拠: [plan §5.1:1512](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1512) の表は **#20 で終わる**一方、本文は未実施作業を **#28・#29・#30** に送っています。離散保存検算が未実装でも #4 は完了、URANS 必須比較が未実施でも #9 は完了です。

   数値の同期も不十分です。§4.7.7 の底面 `G_i0=0.00091 W/K` に対して、現在の [network.json](/home/sano/work/forge/case/49.plate_annular_cavity_m5/_fem_bc/fem_bc_conc/network.json) は **0.00108691 W/K**、対応する `h` は **4.11→4.908 W/m²K**。どの再評価で更新した係数を正本とするか追えません。

   **対案:** 未完了ゲートを優先順付きの実在行へ戻し、正本 run・評価版・VERDICT・係数を結び付ける。plan、case README、FEM 成果物を同じ評価結果から同期する。

10. **Minor — 旋回の回帰試験が本番判定を呼んでいない**

    根拠: [test_symmetry_gate.py:17](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/test_symmetry_gate.py:17) は判定ロジックを別実装しています。付属試験は `VERDICT: PASS` でしたが、本番コードの変更や「3点未満を拒否」の退行を保証しません。

    **対案:** 本番の判定を関数化して試験から直接呼び、短系列・非有限・終点だけ小さい振動を検査する。

**推奨:** `active` に留め、まず入口 EOS と判定・収支ツールを修正し、その後に対象条件の検証と FEM 出力の再生成を行って、result レビューを再実施してください。診断用 `qwall` 追加の確認した分岐には直接の残差変更を認めませんでしたが、それだけで解析結果全体の妥当性は保証されません。ファイルは変更しておらず、本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 9 / Minor 1
