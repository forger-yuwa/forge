# codex レビュー: tooling-rerun-conditions (result)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `5baa5892` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

主要な検証数値と準定常判定は再現できました。ただし、入力検査と設定書き換えに実効条件を誤認する穴が残っています。
`accepted` への移行は、以下の実装修正と回帰確認が終わるまで保留を推奨します。

1. **Major — YAML merge によって設定解釈の不一致と重複キー検査の迂回が残る。**

   根拠: [yaml_strict.py:27](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/yaml_strict.py:27) は `<<` を検査対象から外し、後から PyYAML に展開させています。メモリ内試験で、次が拒否されず `model: none` になることを確認しました。

   ```yaml
   turbulence: {<<: {model: sst, model: none}}
   ```

   また、既存の SST 設定に `<<: {sstEnergyIncludesK: 1}` を追加しても `build_plan` は受理します。検査側は値を 1 と解釈しますが、ソルバは merge を展開せず直接キーを検索し、未定義なら 0 を使います（[solverConfig.cpp:26](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:26)、[同:840](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:840)）。「検査した実効設定」と実際の方程式が食い違います。

   **対案:** v1 では merge key を全階層で拒否する。merge 内の重複と、merge で与えた任意設定の試験を追加してください。

2. **Major — 使用禁止の `mesh.bndFirstOrder` を新しい run に持ち込める。**

   根拠: [rerun_conditions.py:376](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:376) の入力検査に禁止キーの確認がありません。`run_0094` の設定へ `bndFirstOrder: 1` を追加したメモリ内 fixture は、Pt 変更・スケール・推奨 CFL／step 数を指定しても受理され、生成設定に同キーが残りました。

   ソルバは現在もこの値を読みます（[solverConfig.cpp:196](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:196)）。有効時には速度勾配をゼロ化し、その勾配を粘性流束が使用します（[calcGradient_d.cu:471](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/calcGradient_d.cu:471)、[viscousFlux_d.cu:130](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:130)）。単なる記録上の問題ではありません。

   **対案:** キーが存在する参照 run は作成前に拒否する。参照設定を是正してから再実行させ、禁止設定を複製しないことを回帰試験で保証してください。

3. **Major — `--steps` が実効値を変えず、変更済みとして記録される。**

   根拠: [rerun_conditions.py:618](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:618) はコメントも対象に正規表現で置換します。その後は倍数条件だけを検査し、要求した step 数との一致を確認しません（[同:637](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:637)）。

   次の有効な YAML 表記で `--steps 12000` を与えると再現します。

   ```yaml
   last: {nStepOuter : 6000}
   # previous nStepOuter: 6000
   ```

   結果は、コメントだけが `12000` に変わり、`cfg_changes` は `12000`、実効 `nStepOuter` は **6000**、例外なしでした。`outStepInterval` にも同じ欠陥があります。

   **対案:** 対象 YAML パスの値トークンを変更し、再読込後に要求値との完全一致と他設定の不変を検査する。コメント・引用符付きキー・コロン前の空白を含む回帰試験が必要です。

4. **Minor — `species_meta.yaml` の重複キー拒否は、組成変更時にしか働かない。**

   根拠: [rerun_conditions.py:677](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:677) では `"Y" in changes` のときだけ同ファイルを読みます。重複した `species` キーを追加して無変更 rerun を計画すると、拒否されませんでした。「solverConfig・bcondConfig・species_meta の全階層で拒否」という仕様と不一致です。

   **対案:** ファイルが存在すれば変更内容に関係なく strict loader で検査し、書き換えだけを組成変更時に限定してください。

5. **Minor — plan の現行まとめに、訂正済みの「整定余量」が残っている。**

   根拠: [plan:147](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:147) の「§4.7 最終形」は、Tt／Y 変更について依然として「漸近中［残り約 0.1 %］」と記載しています。しかし `run_0138_rerun_fullpath_blk4` は再判定でも **流量・δ_E が DRIFTING**。0.1 % は簡易予想流量との差であり、真の整定余量ではありません。変更ログと README の訂正が現行まとめまで届いていません。

   **対案:** 「流量・δ_E は DRIFTING、簡易予想値との差は約 0.1 %、整定余量は未確定」に統一してください。

検証結果の照合では、次を確認しました。run はすべて `case/45.isobutane_m6_d155/` 配下、実データは `/home/sano/work/forge/` 側です。[run 一覧](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:96)とも照合しました。

| run | 再計算した末尾5枚平均など | `check_convergence` | `check_quasisteady` |
|---|---|---|---|
| `run_0119_rerun_ctrl` | δ_E = 0.725282185、M = 5.999261065 | NOT CONVERGED | 3量 STEADY |
| `run_0120_rerun_euler_pt08` | M = 5.999997716、流量 = 13344.965039 | NOT CONVERGED | M・流量 STEADY |
| `run_0132_rerun_pt08_scale_cfl1_ext` | δ_E = 0.749448414、M = 5.992124197 | NOT CONVERGED | 3量 STEADY |
| `run_0134_rerun_euler_tt1500` | M の直前窓差 = +0.020794759 | NOT CONVERGED | M・流量 DRIFTING |
| `run_0138_rerun_fullpath_blk4` | 流量の直前窓差 = +0.821875 | NOT CONVERGED | M STEADY、流量・δ_E DRIFTING |
| `run_0139_rerun_pt08_noscale_cfl1_ext` | M = 5.993976711、最終場の Ux < 0 は779節点 | NOT CONVERGED | M・流量 DRIFTING |

したがって、A3 の量の準定常到達と B3 の規定時間内未達という限定した結論は支持できます。残差収束を示した結果ではありません。Tt／組成変更の整定長と Euler 参照の再整定は、§5.1 #11 に残っています。

対象 plan によるソルバ数値カーネルの変更はありません。段設定の空白・指数表記・block 形式の処理、および負密度・NaN・Inf を拒否する段終了ゲートは、ファイルを作らない試験で確認できました。ファイル生成を伴う既存テストスイート全体は再実行していません。

**推奨:** `in_progress` を維持し、指摘1〜3を修正・回帰確認した後、4〜5の入力契約と文書を同期して再レビューしてください。今回の NO-GO は、持ち越し済みの物理受入れではなく、準備ツールの実装欠陥によるものです。ファイルは変更しておらず、本レビューの提案は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 2
